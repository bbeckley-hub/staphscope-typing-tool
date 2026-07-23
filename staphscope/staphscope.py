#!/usr/bin/env python3
"""
StaphScope Main Orchestrator - v1.3.2
All module writes happen in /tmp, final results are copied to user output. HPC/ Docker-friendly
Author: Brown Beckley <brownbeckley94@gmail.com>
Affiliation: University of Ghana Medical School
Version: 1.3.2
Date: 2026-07-18
MIT
"""

import os
import sys
import glob
import argparse
import subprocess
import shutil
import tempfile
import logging
import traceback
import signal
from pathlib import Path
from datetime import datetime
from typing import Dict, List

try:
    from .core.banner import StaphScopeBanner
except (ImportError, SystemError):
    sys.path.insert(0, str(Path(__file__).parent))
    from core.banner import StaphScopeBanner


class ColoredHelpFormatter(argparse.HelpFormatter):
    """Custom help formatter with colored output for terminal readability."""
    HEADER = '\033[95m'
    BLUE = '\033[94m'
    GREEN = '\033[92m'
    YELLOW = '\033[93m'
    RED = '\033[91m'
    CYAN = '\033[96m'
    BOLD = '\033[1m'
    RESET = '\033[0m'

    def _format_usage(self, usage, actions, groups, prefix):
        usage_str = super()._format_usage(usage, actions, groups, prefix)
        if usage_str:
            usage_str = self.GREEN + self.BOLD + "usage: " + self.RESET + usage_str
        return usage_str

    def _format_action(self, action):
        action_str = super()._format_action(action)
        if not action_str:
            return action_str
        lines = action_str.split('\n')
        colored_lines = []
        for line in lines:
            if line.strip():
                if line.lstrip().startswith('-'):
                    parts = line.split('  ', 1)
                    if len(parts) == 2:
                        options = parts[0].strip()
                        help_text = parts[1]
                        colored_line = f"  {self.CYAN}{options}{self.RESET}  {help_text}"
                    else:
                        colored_line = f"  {self.CYAN}{line.strip()}{self.RESET}"
                else:
                    colored_line = f"  {self.YELLOW}{line}{self.RESET}"
                colored_lines.append(colored_line)
            else:
                colored_lines.append(line)
        return '\n'.join(colored_lines)

    def _format_text(self, text):
        if not text:
            return text
        return f"{self.BLUE}{text}{self.RESET}"

    def start_section(self, heading):
        heading = f"{self.BOLD}{self.GREEN}{heading}{self.RESET}"
        super().start_section(heading)


class StaphScopeOrchestrator:
    """
    Main orchestrator for StaphScope pipeline.
    Manages temporary directories, runs modules in isolation, and collects results.
    """

    def __init__(self):
        self.banner = StaphScopeBanner()
        self.base_dir = Path(__file__).parent
        self.user_output_dir = None
        self.logger = None
        self.keep_temp = False
        self.temp_dirs = []

        signal.signal(signal.SIGINT, self._signal_handler)
        signal.signal(signal.SIGTERM, self._signal_handler)

    def _signal_handler(self, signum, frame):
        self.banner.display_warning("Interrupted by user. Cleaning up temporary directories...")
        for temp_dir in self.temp_dirs:
            if Path(temp_dir).exists():
                shutil.rmtree(temp_dir, ignore_errors=True)
                self.logger.info(f"Removed temporary directory: {temp_dir}")
        sys.exit(1)

    def _register_temp_dir(self, path: Path):
        self.temp_dirs.append(str(path))

    def setup_logging(self, output_dir: Path):
        log_file = output_dir / "staphscope_run.log"
        logging.basicConfig(
            level=logging.INFO,
            format='%(asctime)s - %(levelname)s - %(message)s',
            handlers=[logging.FileHandler(log_file, mode='w')]
        )
        self.logger = logging.getLogger("StaphScope")
        self.logger.info(f"Logging to {log_file}")
        self.user_output_dir = output_dir

    def get_module_path(self, module_name: str) -> Path:
        if hasattr(sys, 'prefix'):
            share_path = Path(sys.prefix) / "share" / "staphscope" / "modules" / module_name
            if share_path.exists():
                return share_path
        return self.base_dir / "modules" / module_name

    def run_module_in_temp(self, module_name: str, fasta_files: List[Path],
                           cmd_str: str, result_subdir: str = None) -> bool:
        module_orig = self.get_module_path(module_name)
        if not module_orig.exists():
            self.banner.display_error(f"Module directory not found: {module_orig}")
            return False

        temp_dir = Path(tempfile.mkdtemp(prefix=f"staphscope_{module_name}_"))
        self._register_temp_dir(temp_dir)
        self.logger.info(f"Temporary directory for {module_name}: {temp_dir}")

        try:
            shutil.copytree(module_orig, temp_dir / module_name, dirs_exist_ok=True)
            for f in fasta_files:
                shutil.copy2(f, temp_dir / f.name)

            self.banner.display_info(f"Copied {len(fasta_files)} files to {module_name} module")
            pattern = self.get_file_pattern(fasta_files)
            self.banner.display_info(f"Running {module_name} analysis with pattern: {pattern}")

            result = subprocess.run(cmd_str, shell=True, cwd=temp_dir, capture_output=True, text=True)
            if result.stdout:
                self.logger.info(result.stdout)
            if result.stderr:
                self.logger.warning(result.stderr)
            if result.returncode != 0:
                self.logger.error(f"{module_name} failed with return code {result.returncode}")
                return False

            if result_subdir:
                src = temp_dir / result_subdir
                if src.exists():
                    dst = self.user_output_dir / result_subdir
                    if dst.exists():
                        shutil.rmtree(dst)
                    shutil.copytree(src, dst)
                    self.logger.info(f"Results copied to {dst}")

            for extra in ["mutation_summary.html", "mutation_summary.tsv", "mutation_master_summary.json"]:
                extra_src = temp_dir / extra
                if extra_src.exists():
                    shutil.copy2(extra_src, self.user_output_dir / extra)
                    self.logger.info(f"Copied {extra} to output directory")

            return True

        except Exception as e:
            self.logger.error(f"Exception in {module_name}: {e}\n{traceback.format_exc()}")
            return False
        finally:
            if not self.keep_temp:
                shutil.rmtree(temp_dir, ignore_errors=True)
                self.logger.info(f"Removed temporary directory: {temp_dir}")

    # --- Module-specific run methods ---

    def run_fasta_qc_analysis(self, fasta_files: List[Path], output_dir: Path, threads: int) -> bool:
        pattern = self.get_file_pattern(fasta_files)
        cmd = f"{sys.executable} fasta_qc_module/staph_fasta_qc.py {pattern} -o fasta_qc_results -c {threads}"
        return self.run_module_in_temp("fasta_qc_module", fasta_files, cmd, "fasta_qc_results")

    def run_mlst_analysis(self, fasta_files: List[Path], output_dir: Path, threads: int) -> bool:
        pattern_with_quotes = self.get_file_pattern(fasta_files)
        pattern_unquoted = pattern_with_quotes.strip('"')
        temp_dir = Path(tempfile.mkdtemp(prefix="staphscope_mlst_"))
        self._register_temp_dir(temp_dir)
        try:
            module_orig = self.get_module_path("mlst_module")
            shutil.copytree(module_orig, temp_dir, dirs_exist_ok=True)
            for f in fasta_files:
                shutil.copy2(f, temp_dir / f.name)

            self.banner.display_info(f"Copied {len(fasta_files)} files to MLST module")
            self.banner.display_info(f"Running MLST analysis with pattern: {pattern_with_quotes}")

            script = temp_dir / "mlst_module.py"
            cmd = [sys.executable, str(script), "-i", pattern_unquoted, "-o", "mlst_results", "-db", "db", "-sc", "bin", "--batch"]
            result = subprocess.run(cmd, cwd=temp_dir, capture_output=True, text=True)

            if result.stdout:
                self.logger.info(result.stdout)
            if result.stderr:
                self.logger.warning(result.stderr)

            if result.returncode != 0:
                self.logger.error(f"MLST failed: {result.stderr}")
                return False

            src = temp_dir / "mlst_results"
            if src.exists():
                dst = self.user_output_dir / "mlst_results"
                if dst.exists():
                    shutil.rmtree(dst)
                shutil.copytree(src, dst)
                self.logger.info(f"MLST results copied to {dst}")
            return True
        except Exception as e:
            self.logger.error(f"MLST exception: {e}")
            return False
        finally:
            if not self.keep_temp:
                shutil.rmtree(temp_dir, ignore_errors=True)

    def run_spa_typing(self, fasta_files: List[Path], output_dir: Path, threads: int) -> bool:
        temp_dir = Path(tempfile.mkdtemp(prefix="staphscope_spa_"))
        self._register_temp_dir(temp_dir)
        try:
            module_orig = self.get_module_path("spa_module")
            shutil.copytree(module_orig, temp_dir / "spa_module", dirs_exist_ok=True)
            for f in fasta_files:
                shutil.copy2(f, temp_dir / f.name)

            self.banner.display_info(f"Copied {len(fasta_files)} files to spa module")
            self.banner.display_info("Running spa typing analysis...")

            spatyper_dir = temp_dir / "spa_module" / "spatyper"
            script = spatyper_dir / "spa_typing" / "spa_typing_module.py"
            env = os.environ.copy()
            env['PYTHONPATH'] = str(spatyper_dir) + ":" + env.get('PYTHONPATH', '')
            cmd = [sys.executable, str(script), "-i", ".", "-o", "spa_results"]
            result = subprocess.run(cmd, cwd=temp_dir, env=env, capture_output=True, text=True)

            if result.stdout:
                self.logger.info(result.stdout)
            if result.stderr:
                self.logger.warning(result.stderr)

            if result.returncode == 0:
                src = temp_dir / "spa_results"
                if src.exists():
                    dst = self.user_output_dir / "spa_results"
                    if dst.exists():
                        shutil.rmtree(dst)
                    shutil.copytree(src, dst)
                    self.logger.info(f"Spa results copied to {dst}")
                return True
            else:
                self.logger.error(f"spa typing failed: {result.stderr}")
                return False
        except Exception as e:
            self.logger.error(f"spa typing exception: {e}")
            return False
        finally:
            if not self.keep_temp:
                shutil.rmtree(temp_dir, ignore_errors=True)

    def run_sccmec_analysis(self, fasta_files: List[Path], output_dir: Path, threads: int) -> bool:
        pattern_with_quotes = self.get_file_pattern(fasta_files)
        pattern_unquoted = pattern_with_quotes.strip('"')
        temp_dir = Path(tempfile.mkdtemp(prefix="staphscope_sccmec_"))
        self._register_temp_dir(temp_dir)
        try:
            module_orig = self.get_module_path("sccmec_module")
            shutil.copytree(module_orig, temp_dir, dirs_exist_ok=True)
            for f in fasta_files:
                shutil.copy2(f, temp_dir / f.name)

            self.banner.display_info(f"Copied {len(fasta_files)} files to SCCmec module")
            self.banner.display_info(f"Running SCCmec analysis with pattern: {pattern_with_quotes}")

            batch_script = temp_dir / "run_sccmec_batch.sh"
            summary_script = temp_dir / "generate_staphscope_summary.sh"

            cmd_batch = f"bash {batch_script} {pattern_unquoted}"
            result_batch = subprocess.run(cmd_batch, shell=True, cwd=temp_dir, capture_output=True, text=True)
            if result_batch.stdout:
                self.logger.info(result_batch.stdout)
            if result_batch.stderr:
                self.logger.warning(result_batch.stderr)

            if result_batch.returncode != 0:
                self.logger.error(f"SCCmec batch failed: {result_batch.stderr}")
                return False

            if summary_script.exists():
                subprocess.run(f"bash {summary_script}", shell=True, cwd=temp_dir, capture_output=True, text=True)

            target_dir = self.user_output_dir / "sccmec_results"
            target_dir.mkdir(parents=True, exist_ok=True)

            for s_dir in temp_dir.glob("s_*"):
                if s_dir.is_dir():
                    shutil.copytree(s_dir, target_dir / s_dir.name, dirs_exist_ok=True)

            for fname in ["staphscope_summary.html", "staphscope_summary.tsv", "staphscope_detailed_results.csv"]:
                src = temp_dir / fname
                if src.exists():
                    shutil.copy2(src, target_dir / fname)

            self.logger.info(f"SCCmec results copied to {target_dir}")
            return True
        except Exception as e:
            self.logger.error(f"SCCmec exception: {e}")
            return False
        finally:
            if not self.keep_temp:
                shutil.rmtree(temp_dir, ignore_errors=True)

    def run_amrfinder_analysis(self, fasta_files: List[Path], output_dir: Path, threads: int,
                               min_identity: float = None, min_coverage: float = None,
                               skip_mutations: bool = False, force_update: bool = False) -> bool:
        if not self.ensure_amr_database():
            self.banner.display_error("AMR database is missing and could not be updated automatically.")
            return False

        pattern = self.get_file_pattern(fasta_files)
        cmd = f"{sys.executable} amr_module/amrfinder_standalone.py {pattern}"
        if min_identity is not None:
            cmd += f" --min-identity {min_identity}"
        if min_coverage is not None:
            cmd += f" --min-coverage {min_coverage}"
        if skip_mutations:
            cmd += " --skip-mutations"
        if force_update:
            self.banner.display_info("Forcing AMR database update before analysis...")
            self.update_amr_database(force=True)

        return self.run_module_in_temp("amr_module", fasta_files, cmd, "staph_amrfinder_results")

    def run_abricate_analysis(self, fasta_files: List[Path], output_dir: Path, threads: int,
                              min_identity: int = 80, min_coverage: int = 80) -> bool:
        pattern = self.get_file_pattern(fasta_files)
        cmd = (f"{sys.executable} abricate_module/abricate_standalone.py {pattern} "
               f"--minid {min_identity} --mincov {min_coverage}")
        return self.run_module_in_temp("abricate_module", fasta_files, cmd, "abricate_results")

    def run_agr_analysis(self, fasta_files: List[Path], output_dir: Path, threads: int) -> bool:
        module_name = "agr_module"
        module_orig = self.get_module_path(module_name)
        if not module_orig.exists():
            self.banner.display_error(f"Agr module not found: {module_orig}")
            return False

        temp_dir = Path(tempfile.mkdtemp(prefix="staphscope_agr_"))
        self._register_temp_dir(temp_dir)
        self.logger.info(f"Temporary directory for agr: {temp_dir}")

        try:
            shutil.copytree(module_orig, temp_dir / module_name, dirs_exist_ok=True)
            for f in fasta_files:
                shutil.copy2(f, temp_dir / f.name)

            self.banner.display_info(f"Copied {len(fasta_files)} files to agr module")
            self.banner.display_info("Running agr typing analysis...")

            script = temp_dir / module_name / "agr_typing_module.py"
            if not script.exists():
                script = temp_dir / "agr_typing_module.py"
                if not script.exists():
                    self.logger.error("agr_typing_module.py not found in agr_module")
                    return False

            cmd = [sys.executable, str(script), "-i", ".", "-o", "agr_results"]
            result = subprocess.run(cmd, cwd=temp_dir, capture_output=True, text=True)

            if result.stdout:
                self.logger.info(result.stdout)
            if result.stderr:
                self.logger.warning(result.stderr)

            if result.returncode != 0:
                self.logger.error(f"agr typing failed with return code {result.returncode}")
                return False

            src = temp_dir / "agr_results"
            if src.exists():
                dst = self.user_output_dir / "agr_results"
                if dst.exists():
                    shutil.rmtree(dst)
                shutil.copytree(src, dst)
                self.logger.info(f"agr results copied to {dst}")
                self.banner.display_success("✅ Agr typing completed successfully.")
                return True
            else:
                self.logger.error("agr_results directory not found after running agr module")
                return False

        except Exception as e:
            self.logger.error(f"Agr exception: {e}\n{traceback.format_exc()}")
            return False
        finally:
            if not self.keep_temp:
                shutil.rmtree(temp_dir, ignore_errors=True)

    def run_sample_centric_analysis(self, output_dir: Path) -> bool:
        """
        Run the sample‑centric reporter (staphscope_ultimate_samplecentric_reporter.py).
        This module depends on all TSV summary files, mutation_summary.tsv, and comprehensive report files.
        It runs in a temporary directory where all required files are copied.
        """
        module_name = "sample_centric_module"
        module_orig = self.get_module_path(module_name)
        if not module_orig.exists():
            self.banner.display_error(f"Sample‑centric module not found: {module_orig}")
            return False

        temp_dir = Path(tempfile.mkdtemp(prefix="staphscope_sample_centric_"))
        self._register_temp_dir(temp_dir)
        self.logger.info(f"Temporary directory for sample‑centric: {temp_dir}")

        try:
            # Copy the entire module directory
            shutil.copytree(module_orig, temp_dir, dirs_exist_ok=True)

            # Define required files with source paths (relative to output_dir)
            required_files = [
                ("mlst_results", "mlst_summary.tsv"),
                ("spa_results", "spa_summary.tsv"),
                ("sccmec_results", "staphscope_summary.tsv"),
                ("agr_results", "agr_summary.tsv"),
                ("staph_amrfinder_results", "staph_amrfinder_summary.tsv"),
                ("fasta_qc_results", "FASTA_QC_summary.html"),
                # mutation_summary.tsv is inside staph_amrfinder_results
                ("staph_amrfinder_results", "mutation_summary.tsv"),
            ]

            # All ABRicate summary TSVs
            abricate_sources = [
                ("abricate_results", "staph_argannot_abricate_summary.tsv"),
                ("abricate_results", "staph_bacmet2_abricate_summary.tsv"),
                ("abricate_results", "staph_card_abricate_summary.tsv"),
                ("abricate_results", "staph_megares_abricate_summary.tsv"),
                ("abricate_results", "staph_ncbi_abricate_summary.tsv"),
                ("abricate_results", "staph_plasmidfinder_abricate_summary.tsv"),
                ("abricate_results", "staph_resfinder_abricate_summary.tsv"),
                ("abricate_results", "staph_vfdb_abricate_summary.tsv"),
            ]
            required_files.extend(abricate_sources)

            # Copy all required files from their subdirectories
            for subdir, filename in required_files:
                if subdir:
                    src = output_dir / subdir / filename
                else:
                    src = output_dir / filename
                if src.exists():
                    shutil.copy2(src, temp_dir / filename)
                    self.logger.info(f"Copied {filename} to sample‑centric temp dir")
                else:
                    self.logger.warning(f"Required file not found: {src} (skipping)")

            # Copy comprehensive report files (HTML, JSON, TSV) from top-level output
            comprehensive_extensions = [".html", ".json", ".tsv"]
            for ext in comprehensive_extensions:
                src = output_dir / f"staphscope_comprehensive_report{ext}"
                if src.exists():
                    shutil.copy2(src, temp_dir / f"staphscope_comprehensive_report{ext}")
                    self.logger.info(f"Copied staphscope_comprehensive_report{ext} to sample‑centric temp dir")
                else:
                    self.logger.warning(f"staphscope_comprehensive_report{ext} not found; sample‑centric reporter may fail")

            self.banner.display_info("Running sample‑centric reporter...")
            script = temp_dir / "staphscope_ultimate_samplecentric_reporter.py"
            if not script.exists():
                self.logger.error("staphscope_ultimate_samplecentric_reporter.py not found in sample_centric_module")
                return False

            cmd = [sys.executable, str(script), "-i", "."]
            result = subprocess.run(cmd, cwd=temp_dir, capture_output=True, text=True)

            if result.stdout:
                self.logger.info(result.stdout)
            if result.stderr:
                self.logger.warning(result.stderr)

            if result.returncode != 0:
                self.logger.error(f"sample‑centric reporter failed with return code {result.returncode}")
                return False

            # Copy output directory back
            src_dir = temp_dir / "STAPHSCOPE_ULTIMATE_SAMPLE_CENTRIC_REPORTS"
            if src_dir.exists():
                dst_dir = self.user_output_dir / "STAPHSCOPE_ULTIMATE_SAMPLE_CENTRIC_REPORTS"
                if dst_dir.exists():
                    shutil.rmtree(dst_dir)
                shutil.copytree(src_dir, dst_dir)
                self.logger.info(f"Sample‑centric reports copied to {dst_dir}")
                self.banner.display_success("✅ Sample‑centric reporter completed successfully.")
                return True
            else:
                self.logger.error("STAPHSCOPE_ULTIMATE_SAMPLE_CENTRIC_REPORTS not found after running sample‑centric reporter")
                return False

        except Exception as e:
            self.logger.error(f"Sample‑centric exception: {e}\n{traceback.format_exc()}")
            return False
        finally:
            if not self.keep_temp:
                shutil.rmtree(temp_dir, ignore_errors=True)

    def run_lineage_analysis(self, output_dir: Path) -> bool:
        temp_dir = Path(tempfile.mkdtemp(prefix="staphscope_lineage_"))
        self._register_temp_dir(temp_dir)
        try:
            module_orig = self.get_module_path("lineage_module")
            shutil.copytree(module_orig, temp_dir / "lineage_module", dirs_exist_ok=True)
            self.banner.display_info("Generating lineage reference database...")
            cmd = f"{sys.executable} lineage_module/html_reference.py"
            result = subprocess.run(cmd, shell=True, cwd=temp_dir, capture_output=True, text=True)
            if result.stdout:
                self.logger.info(result.stdout)
            if result.stderr:
                self.logger.warning(result.stderr)
            if result.returncode == 0:
                src = temp_dir / "staphscope_lineage_reference.html"
                if src.exists():
                    dst_dir = self.user_output_dir / "lineage_results"
                    dst_dir.mkdir(exist_ok=True)
                    shutil.copy2(src, dst_dir / "staphscope_lineage_reference.html")
                    self.logger.info("Lineage reference copied")
                return True
            else:
                self.logger.error(f"Lineage failed: {result.stderr}")
                return False
        except Exception as e:
            self.logger.error(f"Lineage exception: {e}")
            return False
        finally:
            if not self.keep_temp:
                shutil.rmtree(temp_dir, ignore_errors=True)

    # --- Helper functions ---

    def find_fasta_files(self, input_path: str) -> List[Path]:
        self.banner.display_info(f"Searching for files with pattern: {input_path}")
        if '*' in input_path or '?' in input_path:
            matched_files = glob.glob(input_path)
            fasta_files = [Path(f) for f in matched_files if Path(f).is_file() and
                          f.lower().endswith(('.fna', '.fasta')) and
                          not Path(f).name.startswith('.')]
            self.banner.display_success(f"Found {len(fasta_files)} FASTA files")
            return sorted(fasta_files)
        input_path_obj = Path(input_path)
        if input_path_obj.is_file() and input_path_obj.suffix.lower() in ['.fna', '.fasta']:
            self.banner.display_success(f"Found single FASTA file: {input_path_obj.name}")
            return [input_path_obj]
        if input_path_obj.is_dir():
            patterns = [f"{input_path}/*.fna", f"{input_path}/*.fasta"]
            fasta_files = []
            for pattern in patterns:
                for file_path in glob.glob(pattern):
                    path = Path(file_path)
                    if path.is_file() and not path.name.startswith('.'):
                        fasta_files.append(path)
            fasta_files = sorted(list(set(fasta_files)))
            if fasta_files:
                self.banner.display_success(f"Found {len(fasta_files)} FASTA files in directory")
            else:
                self.banner.display_warning(f"No FASTA files found in directory: {input_path}")
            return fasta_files
        self.banner.display_error(f"Input path not found: {input_path}")
        return []

    def get_file_pattern(self, fasta_files: List[Path]) -> str:
        if not fasta_files:
            return '"*.fna"'
        extensions = set(f.suffix.lower() for f in fasta_files)
        if len(extensions) == 1:
            ext = list(extensions)[0]
            return f'"*{ext}"'
        return '"*"'

    def update_amr_database(self, force: bool = False) -> bool:
        amr_module_path = self.get_module_path("amr_module")
        amr_script = amr_module_path / "amrfinder_standalone.py"
        if not amr_script.exists():
            self.banner.display_error(f"AMR script not found at: {amr_script}")
            return False

        # Ensure logger exists
        if self.logger is None:
            import logging
            logging.basicConfig(level=logging.INFO, format='%(message)s')
            self.logger = logging.getLogger("StaphScope")

        self.banner.display_info("Updating AMRfinderPlus database...")
        flag = "--force-update" if force else "--update-db"
        cmd = [sys.executable, str(amr_script), flag]
        result = subprocess.run(cmd, capture_output=True, text=True, cwd=amr_module_path)
        if result.stdout:
            self.logger.info(result.stdout)
        if result.stderr:
            self.logger.warning(result.stderr)
        if result.returncode == 0:
            self.banner.display_success("AMR database updated successfully.")
            version_cmd = [sys.executable, str(amr_script), "--db-version"]
            version_result = subprocess.run(version_cmd, capture_output=True, text=True, cwd=amr_module_path)
            if version_result.returncode == 0:
                self.banner.display_info(f"New database version: {version_result.stdout.strip()}")
            return True
        else:
            self.banner.display_error("AMR database update failed.")
            if result.stderr:
                self.logger.warning(result.stderr)
            return False

    def ensure_amr_database(self) -> bool:
        amr_module_path = self.get_module_path("amr_module")
        amr_script = amr_module_path / "amrfinder_standalone.py"
        if not amr_script.exists():
            self.banner.display_error("AMR script not found, cannot check database.")
            return False

        if self.logger is None:
            import logging
            logging.basicConfig(level=logging.INFO, format='%(message)s')
            self.logger = logging.getLogger("StaphScope")

        cmd = [sys.executable, str(amr_script), "--db-version"]
        result = subprocess.run(cmd, capture_output=True, text=True, cwd=amr_module_path)
        if result.stdout:
            self.logger.info(result.stdout)
        if result.stderr:
            self.logger.warning(result.stderr)
        if result.returncode == 0 and "Unknown" not in result.stdout and "No database" not in result.stdout:
            self.banner.display_success(f"AMR database already present: {result.stdout.strip()}")
            return True
        else:
            self.banner.display_warning("AMR database not found or outdated. Attempting automatic update...")
            return self.update_amr_database(force=False)

    def run_comprehensive_and_ultimate_reports(self, output_dir: Path) -> bool:
        """
        Run comprehensive_report.py and staphscope_ultimate_reporter.py (gene‑centric).
        Copies all required TSV and HTML files into a temporary directory.
        Now includes agr_summary.tsv.
        """
        temp_dir = Path(tempfile.mkdtemp(prefix="staphscope_summary_"))
        self._register_temp_dir(temp_dir)
        self.logger.info(f"Temporary directory for summary reports: {temp_dir}")

        try:
            summary_module_path = self.get_module_path("summary_module")
            shutil.copytree(summary_module_path, temp_dir, dirs_exist_ok=True)

            # Required TSV files – now includes agr_summary.tsv
            required_tsvs = [
                ("mlst_results", "mlst_summary.tsv"),
                ("spa_results", "spa_summary.tsv"),
                ("sccmec_results", "staphscope_summary.tsv"),
                ("agr_results", "agr_summary.tsv"),
            ]
            for subdir, filename in required_tsvs:
                src = output_dir / subdir / filename
                if src.exists():
                    shutil.copy2(src, temp_dir / filename)
                    self.logger.info(f"Copied {filename} to temporary summary_module directory")
                else:
                    self.logger.warning(f"Required TSV not found: {src}")

            # Required HTML files (all abricate, amr, qc, mutation)
            required_htmls = [
                ("staph_amrfinder_results", "staph_amrfinder_summary_report.html"),
                ("staph_amrfinder_results", "mutation_summary.html"),
                ("abricate_results", "staph_card_summary_report.html"),
                ("abricate_results", "staph_plasmidfinder_summary_report.html"),
                ("abricate_results", "staph_ncbi_summary_report.html"),
                ("abricate_results", "staph_vfdb_summary_report.html"),
                ("abricate_results", "staph_megares_summary_report.html"),
                ("abricate_results", "staph_resfinder_summary_report.html"),
                ("abricate_results", "staph_argannot_summary_report.html"),
                ("abricate_results", "staph_bacmet2_summary_report.html"),
                ("fasta_qc_results", "FASTA_QC_summary.html"),
            ]
            for subdir, filename in required_htmls:
                src = output_dir / subdir / filename
                if src.exists():
                    shutil.copy2(src, temp_dir / filename)
                    self.logger.info(f"Copied {filename} to temporary summary_module directory")
                else:
                    self.logger.warning(f"Required HTML not found: {src}")

            self.banner.display_info("Running comprehensive report...")
            cmd_comp = [sys.executable, "comprehensive_report.py"]
            result_comp = subprocess.run(cmd_comp, cwd=temp_dir, capture_output=True, text=True)
            if result_comp.stdout:
                self.logger.info(result_comp.stdout)
            if result_comp.stderr:
                self.logger.warning(result_comp.stderr)

            if result_comp.returncode != 0:
                self.logger.error(f"Comprehensive report failed:\n{result_comp.stderr}")
                self.banner.display_warning("Comprehensive report failed – continuing with ultimate reporter")
            else:
                self.banner.display_success("Comprehensive report generated successfully!")
                # Copy comprehensive report outputs to user output (top level)
                for ext in [".html", ".json", ".tsv"]:
                    src = temp_dir / f"staphscope_comprehensive_report{ext}"
                    if src.exists():
                        dst = self.user_output_dir / src.name
                        shutil.copy2(src, dst)
                        self.logger.info(f"Copied {src.name} to output directory")

            self.banner.display_info("Running ultimate reporter (gene‑centric)...")
            cmd_ultimate = [sys.executable, "staphscope_ultimate_reporter.py", "-i", "."]
            result_ultimate = subprocess.run(cmd_ultimate, cwd=temp_dir, capture_output=True, text=True)
            if result_ultimate.stdout:
                self.logger.info(result_ultimate.stdout)
            if result_ultimate.stderr:
                self.logger.warning(result_ultimate.stderr)

            if result_ultimate.returncode != 0:
                self.logger.error(f"Ultimate reporter failed:\n{result_ultimate.stderr}")
                self.banner.display_error(f"Ultimate reporter failed with exit code {result_ultimate.returncode}")
                return False
            else:
                self.banner.display_success("Ultimate reporter completed successfully!")
                src_dir = temp_dir / "STAPHSCOPE_ULTIMATE_GENE_CENTRIC_REPORTS"
                if src_dir.exists():
                    dst_dir = self.user_output_dir / "STAPHSCOPE_ULTIMATE_GENE_CENTRIC_REPORTS"
                    if dst_dir.exists():
                        shutil.rmtree(dst_dir)
                    shutil.copytree(src_dir, dst_dir)
                    self.logger.info(f"Ultimate reports copied to {dst_dir}")
                return True

        except Exception as e:
            self.logger.error(f"Summary reports exception: {e}\n{traceback.format_exc()}")
            return False
        finally:
            if not self.keep_temp:
                shutil.rmtree(temp_dir, ignore_errors=True)
                self.logger.info(f"Removed temporary directory: {temp_dir}")

    def copy_summary_results_to_final_directory(self, output_dir: Path):
        """
        Copy only the ultimate reports and comprehensive files into Staphscope_final_report.
        No individual module directories are copied.
        """
        try:
            self.banner.display_info("Copying summary results to final directory...")
            final_report_dir = output_dir / "Staphscope_final_report"
            final_report_dir.mkdir(parents=True, exist_ok=True)

            # Copy comprehensive report files
            comprehensive_files = [
                "staphscope_comprehensive_report.html",
                "staphscope_comprehensive_report.json",
                "staphscope_comprehensive_report.tsv"
            ]
            for file_name in comprehensive_files:
                source_file = output_dir / file_name
                if source_file.exists():
                    shutil.copy2(source_file, final_report_dir / file_name)
                    self.banner.display_info(f"  ✓ Copied: {file_name}")

            # Copy gene‑centric ultimate reports
            ultimate_reports_dir = output_dir / "STAPHSCOPE_ULTIMATE_GENE_CENTRIC_REPORTS"
            if ultimate_reports_dir.exists() and ultimate_reports_dir.is_dir():
                target_ultimate_dir = final_report_dir / "STAPHSCOPE_ULTIMATE_GENE_CENTRIC_REPORTS"
                if target_ultimate_dir.exists():
                    shutil.rmtree(target_ultimate_dir)
                shutil.copytree(ultimate_reports_dir, target_ultimate_dir)
                self.banner.display_info(f"  ✓ Copied: STAPHSCOPE_ULTIMATE_GENE_CENTRIC_REPORTS directory")

            # Copy sample‑centric ultimate reports
            sample_centric_dir = output_dir / "STAPHSCOPE_ULTIMATE_SAMPLE_CENTRIC_REPORTS"
            if sample_centric_dir.exists() and sample_centric_dir.is_dir():
                target_sample_dir = final_report_dir / "STAPHSCOPE_ULTIMATE_SAMPLE_CENTRIC_REPORTS"
                if target_sample_dir.exists():
                    shutil.rmtree(target_sample_dir)
                shutil.copytree(sample_centric_dir, target_sample_dir)
                self.banner.display_info(f"  ✓ Copied: STAPHSCOPE_ULTIMATE_SAMPLE_CENTRIC_REPORTS directory")

            self.banner.display_success(f"✅ All results copied to: {final_report_dir}")
            self.banner.display_info("Summary Reports Generated:")
            for file_path in sorted(final_report_dir.glob("*")):
                if file_path.is_dir():
                    dir_files = list(file_path.glob("*"))
                    self.banner.display_info(f"  📁 {file_path.name} ({len(dir_files)} files)")
                else:
                    self.banner.display_info(f"  📄 {file_path.name}")
        except Exception as e:
            self.banner.display_error(f"Error copying summary results: {str(e)}")

    def run_visualization_analysis(self, output_dir: Path) -> bool:
        temp_dir = Path(tempfile.mkdtemp(prefix="staphscope_visualization_"))
        self._register_temp_dir(temp_dir)
        self.logger.info(f"Temporary directory for visualization: {temp_dir}")
        try:
            vis_module_orig = self.get_module_path("visualization_module")
            shutil.copytree(vis_module_orig, temp_dir, dirs_exist_ok=True)

            required_files = [
                (output_dir / "staph_amrfinder_results", "staph_amrfinder_summary_report.html"),
                (output_dir / "abricate_results", "staph_card_summary_report.html"),
                (output_dir / "abricate_results", "staph_plasmidfinder_summary_report.html"),
                (output_dir / "abricate_results", "staph_ncbi_summary_report.html"),
                (output_dir / "abricate_results", "staph_vfdb_summary_report.html"),
                (output_dir / "abricate_results", "staph_megares_summary_report.html"),
                (output_dir / "abricate_results", "staph_resfinder_summary_report.html"),
                (output_dir / "abricate_results", "staph_argannot_summary_report.html"),
                (output_dir / "agr_results", "agr_summary.tsv"),
                (output_dir / "Staphscope_final_report", "staphscope_comprehensive_report.html"),
            ]
            copied_count = 0
            for source_dir, filename in required_files:
                src = source_dir / filename
                if src.exists():
                    shutil.copy2(src, temp_dir / filename)
                    copied_count += 1
                    self.logger.info(f"Copied {filename} to temporary visualization directory")
                else:
                    self.logger.warning(f"Required file not found: {src}")
            self.banner.display_info(f"Copied {copied_count} files to temporary visualization module")

            vis_script = temp_dir / "staphscope_visualizer.py"
            if not vis_script.exists():
                self.banner.display_error(f"Visualization script not found at: {vis_script}")
                return False

            self.banner.display_info("Running visualization module...")
            cmd = [sys.executable, str(vis_script)]
            result = subprocess.run(cmd, cwd=temp_dir, capture_output=True, text=True)

            if result.stdout:
                self.logger.info(result.stdout)
            if result.stderr:
                self.logger.warning(result.stderr)

            if result.returncode != 0:
                self.logger.error(f"Visualization failed:\n{result.stderr}")
                self.banner.display_warning("Visualization had issues")
                return False

            self.banner.display_success("Visualization completed successfully!")

            vis_output_dir = temp_dir / "STAPHSCOPE_VISUALIZATIONS"
            if vis_output_dir.exists() and vis_output_dir.is_dir():
                target_vis_dir = self.user_output_dir / "STAPHSCOPE_VISUALIZATIONS"
                if target_vis_dir.exists():
                    shutil.rmtree(target_vis_dir)
                shutil.copytree(vis_output_dir, target_vis_dir)
                vis_files = list(vis_output_dir.rglob("*"))
                html_files = [f for f in vis_files if f.suffix == '.html']
                image_files = [f for f in vis_files if f.suffix in ['.png', '.jpg', '.jpeg', '.svg']]
                self.banner.display_success(f"✅ Visualizations copied to: {target_vis_dir}")
                self.banner.display_info(f"   📊 {len(html_files)} HTML reports")
                self.banner.display_info(f"   🖼️  {len(image_files)} visualization images")
            return True

        except Exception as e:
            self.logger.error(f"Visualization exception: {e}\n{traceback.format_exc()}")
            return False
        finally:
            if not self.keep_temp:
                shutil.rmtree(temp_dir, ignore_errors=True)
                self.logger.info(f"Removed temporary directory: {temp_dir}")

    def run_sequential_analyses(self, fasta_files: List[Path], output_dir: Path, threads: int,
                               skip_modules: Dict[str, bool],
                               amr_min_identity: float, amr_min_coverage: float,
                               amr_skip_mutations: bool, amr_force_update: bool,
                               abricate_min_identity: int = 80, abricate_min_coverage: int = 80) -> Dict[str, bool]:
        analysis_functions = [
            ("FASTA QC", self.run_fasta_qc_analysis, "FASTA QC Analysis", "Sequence Quality Control & Statistics", not skip_modules.get('fasta_qc', False)),
            ("MLST", self.run_mlst_analysis, "MLST Analysis", "Multi-Locus Sequence Typing", not skip_modules.get('mlst', False)),
            ("spa typing", self.run_spa_typing, "SPA TYPING ANALYSIS", "Staphylococcal Protein A Typing", not skip_modules.get('spa', False)),
            ("SCCmec", self.run_sccmec_analysis, "SCCMEC ANALYSIS", "Methicillin Resistance Cassette Typing", not skip_modules.get('sccmec', False)),
            ("Agr", self.run_agr_analysis, "AGR TYPING", "Agr Accessory Gene Typing", not skip_modules.get('agr', False)),
            ("AMRFinderPlus", lambda f, o, t: self.run_amrfinder_analysis(f, o, t, amr_min_identity, amr_min_coverage, amr_skip_mutations, amr_force_update),
             "AMR ANALYSIS", "Antimicrobial Resistance Gene Detection", not skip_modules.get('amr', False)),
            ("ABRicate", lambda f, o, t: self.run_abricate_analysis(f, o, t, abricate_min_identity, abricate_min_coverage),
             "ABRICATE ANALYSIS", "Comprehensive Resistance, Plasmid & Virulence Gene Screening", not skip_modules.get('abricate', False))
        ]
        active = [(name, func, header, desc) for name, func, header, desc, enabled in analysis_functions if enabled]
        if not active:
            self.banner.display_warning("All analyses were skipped! Nothing to run.")
            return {}
        self.banner.display_info(f"Running {len(active)} analyses")
        results = {}
        for name, func, header, desc in active:
            self.banner.display_module_header(header, desc)
            try:
                if name == "AMRFinderPlus":
                    success = func(fasta_files, output_dir, max(1, threads // len(active)))
                elif name == "ABRicate":
                    success = func(fasta_files, output_dir, max(1, threads // len(active)))
                else:
                    success = func(fasta_files, output_dir, max(1, threads // len(active)))
                results[name] = success
                if success:
                    self.banner.display_success(f"✅ {name} completed")
                else:
                    self.banner.display_error(f"❌ {name} failed")
            except Exception as e:
                self.banner.display_error(f"❌ {name} failed with exception: {str(e)}")
                results[name] = False
            print()
        return results

    def run_complete_analysis(self, input_path: str, output_dir: str, threads: int = 1,
                             skip_modules: Dict[str, bool] = None,
                             skip_comprehensive: bool = False,
                             skip_visualization: bool = False,
                             update_amr_db_only: bool = False,
                             amr_min_identity: float = None,
                             amr_min_coverage: float = None,
                             amr_skip_mutations: bool = False,
                             amr_force_update: bool = False,
                             clean_output: bool = False,
                             skip_sample_centric: bool = False,
                             abricate_min_identity: int = 80,
                             abricate_min_coverage: int = 80):
        if skip_modules is None:
            skip_modules = {}
        if update_amr_db_only:
            self.update_amr_database(force=False)
            return

        start_time = datetime.now()
        try:
            self.banner.display_startup_sequence()
            self.banner.display_banner(show_quote=True, show_author=True)
            output_path = Path(output_dir)
            if clean_output and output_path.exists():
                shutil.rmtree(output_path)
            output_path.mkdir(parents=True, exist_ok=True)
            self.setup_logging(output_path)

            fasta_files = self.find_fasta_files(input_path)
            if not fasta_files:
                self.banner.display_error("No FASTA files found! Analysis stopped.")
                return

            extensions = set(f.suffix.lower() for f in fasta_files)
            self.banner.display_success(f"Starting analysis of {len(fasta_files)} samples")
            self.banner.display_info(f"File formats detected: {', '.join(extensions)}")

            # Create output subdirectories (only for intermediate results)
            subdirs = ["fasta_qc_results", "mlst_results", "spa_results", "sccmec_results",
                       "abricate_results", "staph_amrfinder_results", "lineage_results",
                       "agr_results"]
            for subdir in subdirs:
                (output_path / subdir).mkdir(exist_ok=True)

            # Display enabled modules
            self.banner.display_module_header("Analysis Plan", "Modules to be executed")
            analyses_to_run = [
                ("FASTA QC", not skip_modules.get('fasta_qc', False)),
                ("MLST", not skip_modules.get('mlst', False)),
                ("spa typing", not skip_modules.get('spa', False)),
                ("SCCmec", not skip_modules.get('sccmec', False)),
                ("Agr typing", not skip_modules.get('agr', False)),
                ("AMRFinderPlus", not skip_modules.get('amr', False)),
                ("ABRicate", not skip_modules.get('abricate', False)),
                ("Lineage Reference", not skip_modules.get('lineage', False)),
                ("Comprehensive Report", not skip_comprehensive),
                ("Ultimate Reporter (gene‑centric)", not skip_comprehensive),
                ("Sample‑Centric Reporter", not skip_sample_centric),
                ("Visualization", not skip_visualization)
            ]
            for analysis, enabled in analyses_to_run:
                status = "✅ ENABLED" if enabled else "⏸️  SKIPPED"
                print(f"   {status} - {analysis}")
            sys.stdout.flush()

            # Run primary modules
            analysis_results = self.run_sequential_analyses(fasta_files, output_path, threads, skip_modules,
                                                           amr_min_identity, amr_min_coverage,
                                                           amr_skip_mutations, amr_force_update,
                                                           abricate_min_identity, abricate_min_coverage)

            # Lineage
            if not skip_modules.get('lineage', False):
                self.banner.display_module_header("Lineage Database", "S. aureus Lineage Reference Generation")
                lineage_success = self.run_lineage_analysis(output_path)
                analysis_results["Lineage Reference"] = lineage_success
                print()

            # Comprehensive + Ultimate (gene‑centric)
            if not skip_comprehensive:
                self.banner.display_module_header("Comprehensive & Ultimate Reports", "Unified MLST, spa, SCCmec, agr and gene‑centric integration")
                summary_success = self.run_comprehensive_and_ultimate_reports(output_path)
                analysis_results["Comprehensive & Ultimate Reports"] = summary_success
                if not summary_success:
                    self.banner.display_warning("Summary reports had issues")
                print()

            # Sample‑centric reporter – depends on comprehensive HTML and all TSVs
            if not skip_sample_centric:
                if skip_comprehensive:
                    self.banner.display_warning("Sample‑centric reporter requires staphscope_comprehensive_report.html.")
                    self.banner.display_warning("Since --skip-comprehensive was used, sample‑centric reporter will be skipped.")
                    analysis_results["Sample‑Centric Reporter"] = False
                else:
                    self.banner.display_module_header("Sample‑Centric Reporter", "Interactive isolate‑centric report")
                    sample_success = self.run_sample_centric_analysis(output_path)
                    analysis_results["Sample‑Centric Reporter"] = sample_success
                    print()

            # Now copy all results into Staphscope_final_report (only the required items)
            self.copy_summary_results_to_final_directory(output_path)

            # Clean up top-level duplicates and intermediate directories
            for dup in ["staphscope_comprehensive_report.html",
                        "staphscope_comprehensive_report.json",
                        "staphscope_comprehensive_report.tsv"]:
                (output_path / dup).unlink(missing_ok=True)
            shutil.rmtree(output_path / "STAPHSCOPE_ULTIMATE_GENE_CENTRIC_REPORTS", ignore_errors=True)
            shutil.rmtree(output_path / "STAPHSCOPE_ULTIMATE_SAMPLE_CENTRIC_REPORTS", ignore_errors=True)

            # Visualization (runs after final cleanup, using files from Staphscope_final_report)
            if not skip_visualization:
                self.banner.display_module_header("Visualization", "Interactive Visualizations & Dashboard")
                visualization_success = self.run_visualization_analysis(output_path)
                analysis_results["Visualization"] = visualization_success
                print()

            # Final summary
            analysis_time = datetime.now() - start_time
            analysis_time_str = str(analysis_time).split('.')[0]
            successful_count = sum(analysis_results.values())
            total_count = len(analysis_results)

            self.banner.display_footer(analysis_time=analysis_time_str, samples_processed=len(fasta_files))

            if successful_count == total_count:
                self.banner.display_success(f"🎉 All {total_count} analyses completed successfully!")
                self.banner.display_success("🧹 All module directories have been cleaned up")
                print("\n📁 FINAL OUTPUT STRUCTURE:")
                # Show only Staphscope_final_report and any extra directories
                for subdir in sorted(output_path.glob("*")):
                    if subdir.is_dir():
                        files_count = len(list(subdir.rglob("*")))
                        if subdir.name == "Staphscope_final_report":
                            print(f"   📂 {subdir.name}/")
                            for item in sorted(subdir.glob("*")):
                                if item.is_file():
                                    print(f"      📄 {item.name}")
                                elif item.is_dir():
                                    sub_items = len(list(item.glob("*")))
                                    print(f"      📂 {item.name}/ ({sub_items} items)")
                        elif subdir.name not in ["STAPHSCOPE_VISUALIZATIONS", "staphscope_run.log"]:
                            print(f"   📂 {subdir.name}/ ({files_count} items)")
            else:
                self.banner.display_warning(f"⚠️  {successful_count}/{total_count} analyses completed successfully.")

            print("\n📚 Please cite our StaphScope paper:")
            print("   Beckley, B., Amarh, V. StaphScope: a species-optimized computational pipeline for rapid and accessible Staphylococcus aureus genotyping and surveillance. BMC Genomics (2026).")
            print("   https://doi.org/10.1186/s12864-026-12609-x")

        except KeyboardInterrupt:
            self.banner.display_error("Analysis interrupted by user")
        except Exception as e:
            self.banner.display_error(f"Critical error in analysis pipeline: {str(e)}")
            import traceback
            traceback.print_exc()


def main():
    parser = argparse.ArgumentParser(
        description="StaphScope: Advanced Staphylococcus aureus Typing & Lineage Analysis Platform with FASTA QC, Agr typing, and Visualization",
        formatter_class=ColoredHelpFormatter,
        epilog=f"""
{ColoredHelpFormatter.BOLD}{ColoredHelpFormatter.GREEN}Examples:{ColoredHelpFormatter.RESET}
  # Basic analysis
  staphscope -i genome.fna -o results/

  # Batch with glob pattern
  staphscope -i "*.fna" -o batch_results --threads 8

  # Skip some modules
  staphscope -i "*.fasta" -o analysis --threads 16 --skip-lineage --skip-visualization

  # AMR with custom thresholds and no mutation reporting
  staphscope -i "*.fna" -o results --amr-min-identity 0.95 --amr-min-coverage 0.9 --skip-amr-mutations

  # Force update AMR database before analysis
  staphscope -i "*.fna" -o results --amr-force-update

  # Update AMR database only (standalone)
  staphscope --update-amr-db

  # Force update AMR database only
  staphscope --force-update-amr-db

  # Skip the new sample‑centric reporter
  staphscope -i "*.fna" -o results --skip-sample-centric

  # Skip agr typing
  staphscope -i "*.fna" -o results --skip-agr

  # ABRicate custom thresholds
  staphscope -i "*.fna" -o results --abricate-minid 85 --abricate-mincov 90

{ColoredHelpFormatter.BOLD}{ColoredHelpFormatter.GREEN}Supported FASTA formats:{ColoredHelpFormatter.RESET} .fna, .fasta

{ColoredHelpFormatter.BOLD}{ColoredHelpFormatter.GREEN}Analysis Modules:{ColoredHelpFormatter.RESET}
  • FASTA QC (Quality Control & Statistics)
  • MLST (Multi-Locus Sequence Typing)
  • spa typing (Staphylococcal Protein A)
  • SCCmec typing (Methicillin Resistance Cassette)
  • Agr typing (agrVATE) – NEW
  • AMR profiling (AMRFinderPlus) – point mutations reported by default (use --skip-amr-mutations to disable)
  • ABRicate (Comprehensive resistance/plasmid/virulence)
  • Lineage reference database
  • Comprehensive report (MLST + spa + SCCmec + agr summary)
  • Ultimate reporter (Gene-centric integrated analysis)
  • Sample‑centric reporter (Isolate‑centric interactive report) 
  • Visualization (Interactive dashboards & plots)

{ColoredHelpFormatter.BOLD}{ColoredHelpFormatter.GREEN}Output:{ColoredHelpFormatter.RESET}
  Comprehensive results for all analyses in organized directories.
  A detailed log file (staphscope_run.log) is written to the output directory.

{ColoredHelpFormatter.BOLD}{ColoredHelpFormatter.GREEN}Citation:{ColoredHelpFormatter.RESET}
  Beckley B, Amarh V. StaphScope: a species-optimized computational pipeline for rapid and accessible
  Staphylococcus aureus genotyping and surveillance. BMC Genomics. 2026;27:261.
  doi:10.1186/s12864-026-12609-x

{ColoredHelpFormatter.BOLD}{ColoredHelpFormatter.GREEN}Support & Contributions:{ColoredHelpFormatter.RESET}
  • Issues & feature requests: https://github.com/bbeckley-hub/staphscope/issues
  • Email: brownbeckley94@gmail.com

{ColoredHelpFormatter.YELLOW}⭐ Star us on GitHub if you find this tool useful! ⭐{ColoredHelpFormatter.RESET}
        """
    )
    parser.add_argument('-i', '--input', help='Input FASTA file(s) - can use glob patterns like "*.fna" or "*.fasta"')
    parser.add_argument('-o', '--output', help='Output directory for all results')
    parser.add_argument('-t', '--threads', type=int, default=2, help='Number of threads (default: 2)')
    parser.add_argument('--skip-fasta-qc', action='store_true', help='Skip FASTA QC analysis')
    parser.add_argument('--skip-mlst', action='store_true', help='Skip MLST analysis')
    parser.add_argument('--skip-spa', action='store_true', help='Skip spa typing analysis')
    parser.add_argument('--skip-sccmec', action='store_true', help='Skip SCCmec analysis')
    parser.add_argument('--skip-agr', action='store_true', help='Skip agr typing analysis')
    parser.add_argument('--skip-amr', action='store_true', help='Skip AMR analysis (AMRfinderPlus)')
    parser.add_argument('--skip-abricate', action='store_true', help='Skip ABRicate analysis')
    parser.add_argument('--skip-lineage', action='store_true', help='Skip lineage reference generation')
    parser.add_argument('--skip-comprehensive', action='store_true', help='Skip comprehensive report AND ultimate reporter (gene‑centric)')
    parser.add_argument('--skip-sample-centric', action='store_true', help='Skip sample‑centric reporter (isolate‑centric)')
    parser.add_argument('--skip-visualization', action='store_true', help='Skip visualization module (dashboards & plots)')
    parser.add_argument('--amr-min-identity', type=float, help='Minimum identity for AMR hits (0..1)')
    parser.add_argument('--amr-min-coverage', type=float, help='Minimum coverage for AMR hits (0..1)')
    parser.add_argument('--skip-amr-mutations', action='store_true', help='Disable point mutation reporting in AMR (enabled by default)')
    parser.add_argument('--amr-force-update', action='store_true', help='Force update AMR database before analysis')
    parser.add_argument('--update-amr-db', action='store_true', help='Update AMRfinderPlus database (incremental) and exit')
    parser.add_argument('--force-update-amr-db', action='store_true', help='Force complete AMR database update (overwrites old) and exit')
    parser.add_argument('--keep-temp', action='store_true', help='Do not delete temporary directories (for debugging)')
    parser.add_argument('--clean-output', action='store_true', help='Delete output directory before analysis (prevents mixing results from different runs)')
    # New ABRicate threshold flags
    parser.add_argument('--abricate-minid', type=int, default=80,
                        help='Minimum identity for ABRicate hits (default: 80)')
    parser.add_argument('--abricate-mincov', type=int, default=80,
                        help='Minimum coverage for ABRicate hits (default: 80)')

    args = parser.parse_args()

    if args.update_amr_db or args.force_update_amr_db:
        orch = StaphScopeOrchestrator()
        if args.force_update_amr_db:
            orch.update_amr_database(force=True)
        else:
            orch.update_amr_database(force=False)
        sys.exit(0)

    if not args.input or not args.output:
        parser.error("When not using --update-amr-db or --force-update-amr-db, both -i/--input and -o/--output are required.")

    skip_modules = {
        'fasta_qc': args.skip_fasta_qc,
        'mlst': args.skip_mlst,
        'spa': args.skip_spa,
        'sccmec': args.skip_sccmec,
        'agr': args.skip_agr,
        'amr': args.skip_amr,
        'abricate': args.skip_abricate,
        'lineage': args.skip_lineage,
    }

    orch = StaphScopeOrchestrator()
    orch.keep_temp = args.keep_temp
    orch.run_complete_analysis(
        input_path=args.input,
        output_dir=args.output,
        threads=args.threads,
        skip_modules=skip_modules,
        skip_comprehensive=args.skip_comprehensive,
        skip_visualization=args.skip_visualization,
        amr_min_identity=args.amr_min_identity,
        amr_min_coverage=args.amr_min_coverage,
        amr_skip_mutations=args.skip_amr_mutations,
        amr_force_update=args.amr_force_update,
        clean_output=args.clean_output,
        skip_sample_centric=args.skip_sample_centric,
        abricate_min_identity=args.abricate_minid,
        abricate_min_coverage=args.abricate_mincov
    )


if __name__ == "__main__":
    main()
#!/usr/bin/env python3
"""
StaphScope Main Orchestrator - v2.0.0
All module writes happen in /tmp, final results are copied to user output.
HPC / Docker-friendly.

Modules:
  analysis      : fasta_qc, mlst, spa, sccmec_cge, sccmec_rpet, capsule, agr, amr, abricate, mge
  reporting     : gene_centric_module (comprehensive_report.py + ultimate gene-centric reporter)
                  sample_centric_module (interactive isolate-centric reporter)
  support       : lineage_module, visualization_module
Find a bug? Please reach out!!!
Author: Brown Beckley <brownbeckley94@gmail.com>
Affiliation: University of Ghana Medical School
Version: 2.0.0
Date: 2026-09-19
MIT License
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
import urllib.request
import zipfile
from pathlib import Path
from datetime import datetime
from typing import Dict, List

try:
    from .core.banner import StaphScopeBanner
except (ImportError, SystemError):
    sys.path.insert(0, str(Path(__file__).parent))
    from core.banner import StaphScopeBanner

__version__ = "2.0.0"


class ColoredHelpFormatter(argparse.HelpFormatter):
    """Custom help formatter with ANSI color codes for better readability."""

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
    Main orchestrator for the StaphScope pipeline.
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
        env_override = os.environ.get(f'STAPHSCOPE_{module_name.upper()}_MODULE_PATH')
        if env_override:
            return Path(env_override)

        if hasattr(sys, 'prefix'):
            share_path = Path(sys.prefix) / "share" / "staphscope" / "modules" / module_name
            if share_path.exists():
                return share_path
        return self.base_dir / "modules" / module_name

    def get_mlst_db_dir(self, for_update: bool = False) -> Path:
        """Return the MLST database directory.

        Resolution order:
          1. $STAPHSCOPE_MLST_DB env override (wins in all cases)
          2. User-local module db/ if it contains a populated scheme
          3. Packaged module db/ as fallback

        When for_update=True, always return a writable path (never the
        packaged tree) so updates can't fail on read-only installs.
        """
        env_db = os.environ.get("STAPHSCOPE_MLST_DB")
        if env_db:
            env_path = Path(env_db)
            if for_update:
                env_path.mkdir(parents=True, exist_ok=True)
                return env_path
            if (env_path / "pubmlst" / "saureus").exists():
                return env_path

        user_module = Path(os.environ.get(
            "STAPHSCOPE_MLST_MODULE_PATH",
            Path.home() / ".local" / "share" / "staphscope" / "mlst_module",
        ))
        user_db = user_module / "db"

        if for_update:
            user_db.mkdir(parents=True, exist_ok=True)
            return user_db

        if (user_db / "pubmlst" / "saureus").exists():
            self.banner.display_info(f"Using user-local MLST database: {user_db}")
            return user_db

        module_db = self.get_module_path("mlst_module") / "db"
        if module_db.exists():
            return module_db

        self.banner.display_warning(
            "No populated MLST database found. Run `staphscope --pull-mlst-db`."
        )
        user_db.mkdir(parents=True, exist_ok=True)
        return user_db

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

    # ------------------------------------------------------------------
    # Analysis modules
    # ------------------------------------------------------------------

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

            db_dir = self.get_mlst_db_dir()
            env = os.environ.copy()
            env['STAPHSCOPE_MLST_DB'] = str(db_dir)

            self.banner.display_info(f"Copied {len(fasta_files)} files to MLST module")
            self.banner.display_info(f"Running MLST analysis with pattern: {pattern_with_quotes}")
            self.banner.display_info(f"Using MLST database: {db_dir}")

            script = temp_dir / "mlst_module.py"
            cmd = [sys.executable, str(script), "-i", pattern_unquoted,
                   "-o", "mlst_results", "-db", str(db_dir), "-sc", "bin", "--batch"]
            result = subprocess.run(cmd, cwd=temp_dir, env=env, capture_output=True, text=True)

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
        """Run the CGE SCCmec caller in an isolated temp directory."""
        pattern_with_quotes = self.get_file_pattern(fasta_files)
        pattern_unquoted = pattern_with_quotes.strip('"')
        temp_dir = Path(tempfile.mkdtemp(prefix="staphscope_sccmec_"))
        self._register_temp_dir(temp_dir)
        try:
            module_orig = self.get_module_path("sccmec_module_cge")
            if not module_orig.exists():
                self.banner.display_error(f"sccmec_module_cge not found: {module_orig}")
                return False
            shutil.copytree(module_orig, temp_dir, dirs_exist_ok=True)
            for f in fasta_files:
                shutil.copy2(f, temp_dir / f.name)

            self.banner.display_info(f"Copied {len(fasta_files)} files to SCCmec (CGE) module")
            self.banner.display_info(f"Running SCCmec (CGE) analysis with pattern: {pattern_with_quotes}")

            batch_script = temp_dir / "run_sccmec_batch.sh"
            summary_script = temp_dir / "generate_staphscope_summary.sh"

            cmd_batch = f"bash {batch_script} {pattern_unquoted}"
            result_batch = subprocess.run(cmd_batch, shell=True, cwd=temp_dir, capture_output=True, text=True)
            if result_batch.stdout:
                self.logger.info(result_batch.stdout)
            if result_batch.stderr:
                self.logger.warning(result_batch.stderr)

            if result_batch.returncode != 0:
                self.logger.error(f"SCCmec (CGE) batch failed: {result_batch.stderr}")
                return False

            if summary_script.exists():
                subprocess.run(f"bash {summary_script}", shell=True, cwd=temp_dir,
                               capture_output=True, text=True)

            target_dir = self.user_output_dir / "sccmec_cge_results"
            target_dir.mkdir(parents=True, exist_ok=True)

            for s_dir in temp_dir.glob("s_*"):
                if s_dir.is_dir():
                    shutil.copytree(s_dir, target_dir / s_dir.name, dirs_exist_ok=True)

            for fname in ["staphscope_sccmec_cge_summary.html",
                          "staphscope_sccmec_cge_summary.tsv",
                          "staphscope_sccmec_cge_detailed_results.csv"]:
                src = temp_dir / fname
                if src.exists():
                    shutil.copy2(src, target_dir / fname)

            self.logger.info(f"SCCmec (CGE) results copied to {target_dir}")
            return True
        except Exception as e:
            self.logger.error(f"SCCmec (CGE) exception: {e}")
            return False
        finally:
            if not self.keep_temp:
                shutil.rmtree(temp_dir, ignore_errors=True)

    def run_sccmec_rpet_analysis(self, fasta_files: List[Path], output_dir: Path, threads: int) -> bool:
        """Run the RPet SCCmec caller (Robert Petit III) in an isolated temp directory."""
        pattern_with_quotes = self.get_file_pattern(fasta_files)
        pattern_unquoted = pattern_with_quotes.strip('"')
        temp_dir = Path(tempfile.mkdtemp(prefix="staphscope_sccmec_rpet_"))
        self._register_temp_dir(temp_dir)
        self.logger.info(f"Temporary directory for sccmec_rpet: {temp_dir}")

        try:
            module_orig = self.get_module_path("sccmec_module_rpet")
            if not module_orig.exists():
                self.banner.display_error(f"sccmec_module_rpet not found: {module_orig}")
                return False

            shutil.copytree(module_orig, temp_dir, dirs_exist_ok=True,
                            ignore=shutil.ignore_patterns('results',
                                                          'sccmec_rpet_results',
                                                          '*.fna',
                                                          '__pycache__'))
            for f in fasta_files:
                shutil.copy2(f, temp_dir / f.name)

            self.banner.display_info(f"Copied {len(fasta_files)} files to SCCmec RPet module")
            self.banner.display_info(f"Running SCCmec RPet analysis with pattern: {pattern_with_quotes}")

            script = temp_dir / "SCCmecFinder_RPet.py"
            if not script.exists():
                self.logger.error("SCCmecFinder_RPet.py not found in sccmec_module_rpet")
                return False

            cmd = f"{sys.executable} {script} -i {pattern_unquoted} -d sccmec_rpet_results -db_dir database"
            result = subprocess.run(cmd, shell=True, cwd=temp_dir, capture_output=True, text=True)

            if result.stdout:
                self.logger.info(result.stdout)
            if result.stderr:
                self.logger.warning(result.stderr)

            if result.returncode != 0:
                self.logger.error(f"SCCmec RPet failed with return code {result.returncode}")
                return False

            src = temp_dir / "sccmec_rpet_results"
            if src.exists():
                dst = self.user_output_dir / "sccmec_rpet_results"
                if dst.exists():
                    shutil.rmtree(dst)
                shutil.copytree(src, dst)
                self.logger.info(f"SCCmec RPet results copied to {dst}")
                self.banner.display_success("✅ SCCmec RPet analysis completed successfully.")
                return True
            else:
                self.logger.error("sccmec_rpet_results directory not found after running the module")
                return False

        except Exception as e:
            self.logger.error(f"SCCmec RPet exception: {e}\n{traceback.format_exc()}")
            return False
        finally:
            if not self.keep_temp:
                shutil.rmtree(temp_dir, ignore_errors=True)

    def run_capsule_analysis(self, fasta_files: List[Path], output_dir: Path, threads: int) -> bool:
        """Run the capsule typing module in an isolated temp directory."""
        pattern_with_quotes = self.get_file_pattern(fasta_files)
        pattern_unquoted = pattern_with_quotes.strip('"')
        temp_dir = Path(tempfile.mkdtemp(prefix="staphscope_capsule_"))
        self._register_temp_dir(temp_dir)
        self.logger.info(f"Temporary directory for capsule: {temp_dir}")

        try:
            module_orig = self.get_module_path("capsule_module")
            if not module_orig.exists():
                self.banner.display_error(f"capsule_module not found: {module_orig}")
                return False

            shutil.copytree(module_orig, temp_dir, dirs_exist_ok=True,
                            ignore=shutil.ignore_patterns('capsule_results',
                                                          '*.fna',
                                                          '__pycache__'))
            for f in fasta_files:
                shutil.copy2(f, temp_dir / f.name)

            self.banner.display_info(f"Copied {len(fasta_files)} files to capsule module")
            self.banner.display_info(f"Running capsule typing analysis with pattern: {pattern_with_quotes}")

            script = temp_dir / "CapsuleFinder.py"
            if not script.exists():
                self.logger.error("CapsuleFinder.py not found in capsule_module")
                return False

            cmd = f"{sys.executable} {script} -i {pattern_unquoted} -d capsule_results -db_dir database"
            result = subprocess.run(cmd, shell=True, cwd=temp_dir, capture_output=True, text=True)

            if result.stdout:
                self.logger.info(result.stdout)
            if result.stderr:
                self.logger.warning(result.stderr)

            if result.returncode != 0:
                self.logger.error(f"Capsule typing failed with return code {result.returncode}")
                return False

            src = temp_dir / "capsule_results"
            if src.exists():
                dst = self.user_output_dir / "capsule_results"
                if dst.exists():
                    shutil.rmtree(dst)
                shutil.copytree(src, dst)
                self.logger.info(f"Capsule results copied to {dst}")
                self.banner.display_success("✅ Capsule typing completed successfully.")
                return True
            else:
                self.logger.error("capsule_results directory not found after running the module")
                return False

        except Exception as e:
            self.logger.error(f"Capsule exception: {e}\n{traceback.format_exc()}")
            return False
        finally:
            if not self.keep_temp:
                shutil.rmtree(temp_dir, ignore_errors=True)

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

    def run_mge_analysis(self, fasta_files: List[Path], output_dir: Path, threads: int) -> bool:
        """Run the mobileOG-based MGE profiler in an isolated temp directory."""
        pattern_with_quotes = self.get_file_pattern(fasta_files)
        pattern_unquoted = pattern_with_quotes.strip('"')
        temp_dir = Path(tempfile.mkdtemp(prefix="staphscope_mge_"))
        self._register_temp_dir(temp_dir)
        self.logger.info(f"Temporary directory for mge: {temp_dir}")

        try:
            module_orig = self.get_module_path("mge_module")
            if not module_orig.exists():
                self.banner.display_error(f"mge_module not found: {module_orig}")
                return False

            shutil.copytree(module_orig, temp_dir, dirs_exist_ok=True,
                            ignore=shutil.ignore_patterns('mge_results',
                                                          'test_hits.tsv',
                                                          '*.fna',
                                                          '__pycache__'))
            for f in fasta_files:
                shutil.copy2(f, temp_dir / f.name)

            self.banner.display_info(f"Copied {len(fasta_files)} files to MGE module")
            self.banner.display_info(f"Running MGE analysis with pattern: {pattern_with_quotes}")

            script = temp_dir / "MGEFinder.py"
            if not script.exists():
                self.logger.error("MGEFinder.py not found in mge_module")
                return False

            cmd = f"{sys.executable} {script} -i {pattern_unquoted} -d mge_results -db_dir mobileOG-db/beatrix-1-6_v1_all"
            result = subprocess.run(cmd, shell=True, cwd=temp_dir, capture_output=True, text=True)

            if result.stdout:
                self.logger.info(result.stdout)
            if result.stderr:
                self.logger.warning(result.stderr)

            if result.returncode != 0:
                self.logger.error(f"MGE analysis failed with return code {result.returncode}")
                return False

            src = temp_dir / "mge_results"
            if src.exists():
                dst = self.user_output_dir / "mge_results"
                if dst.exists():
                    shutil.rmtree(dst)
                shutil.copytree(
                    src, dst,
                    ignore=shutil.ignore_patterns('annotations'),
                )
                self.logger.info(f"MGE results copied to {dst} (annotations cache skipped)")
                self.banner.display_success("✅ MGE analysis completed successfully.")
                return True
            else:
                self.logger.error("mge_results directory not found after running the module")
                return False

        except Exception as e:
            self.logger.error(f"MGE exception: {e}\n{traceback.format_exc()}")
            return False
        finally:
            if not self.keep_temp:
                shutil.rmtree(temp_dir, ignore_errors=True)

    # ------------------------------------------------------------------
    # Reporting modules
    # ------------------------------------------------------------------

    def run_comprehensive_and_ultimate_reports(self, output_dir: Path) -> bool:
        """Run comprehensive_report.py and the gene-centric ultimate reporter from gene_centric_module."""
        temp_dir = Path(tempfile.mkdtemp(prefix="staphscope_gene_centric_"))
        self._register_temp_dir(temp_dir)
        self.logger.info(f"Temporary directory for gene-centric reports: {temp_dir}")

        try:
            gene_centric_path = self.get_module_path("gene_centric_module")
            if not gene_centric_path.exists():
                self.banner.display_error(f"gene_centric_module not found: {gene_centric_path}")
                return False

            shutil.copytree(
                gene_centric_path, temp_dir, dirs_exist_ok=True,
                ignore=shutil.ignore_patterns(
                    'STAPHSCOPE_ULTIMATE_GENE_CENTRIC_REPORTS',
                    '*.fna', '__pycache__', '*.pyc',
                ),
            )

            required_tsvs = [
                ("mlst_results", "mlst_summary.tsv"),
                ("spa_results", "spa_summary.tsv"),
                ("sccmec_cge_results", "staphscope_sccmec_cge_summary.tsv"),
                ("sccmec_rpet_results", "staphscope_sccmec_rpet_summary.tsv"),
                ("capsule_results", "staphscope_capsule_summary.tsv"),
                ("agr_results", "agr_summary.tsv"),
            ]
            for subdir, filename in required_tsvs:
                src = output_dir / subdir / filename
                if src.exists():
                    shutil.copy2(src, temp_dir / filename)
                    self.logger.info(f"Copied {filename} to gene-centric temp dir")
                else:
                    self.logger.warning(f"Required TSV not found: {src}")

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
                ("mge_results", "staphscope_mge_summary.html"),
            ]
            for subdir, filename in required_htmls:
                src = output_dir / subdir / filename
                if src.exists():
                    shutil.copy2(src, temp_dir / filename)
                    self.logger.info(f"Copied {filename} to gene-centric temp dir")
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
                for ext in [".html", ".json", ".tsv"]:
                    src = temp_dir / f"staphscope_comprehensive_report{ext}"
                    if src.exists():
                        dst = self.user_output_dir / src.name
                        shutil.copy2(src, dst)
                        self.logger.info(f"Copied {src.name} to output directory")

            self.banner.display_info("Running ultimate reporter (gene-centric)...")
            cmd_ultimate = [sys.executable, "staphscope_ultimate_gene_centric_reporter.py", "-i", "."]
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
            self.logger.error(f"Gene-centric reports exception: {e}\n{traceback.format_exc()}")
            return False
        finally:
            if not self.keep_temp:
                shutil.rmtree(temp_dir, ignore_errors=True)
                self.logger.info(f"Removed temporary directory: {temp_dir}")

    def run_sample_centric_analysis(self, output_dir: Path) -> bool:
        module_name = "sample_centric_module"
        module_orig = self.get_module_path(module_name)
        if not module_orig.exists():
            self.banner.display_error(f"Sample-centric module not found: {module_orig}")
            return False

        temp_dir = Path(tempfile.mkdtemp(prefix="staphscope_sample_centric_"))
        self._register_temp_dir(temp_dir)
        self.logger.info(f"Temporary directory for sample-centric: {temp_dir}")

        try:
            shutil.copytree(module_orig, temp_dir, dirs_exist_ok=True)

            required_files = [
                ("staph_amrfinder_results", "staph_amrfinder_summary.tsv"),
                ("staph_amrfinder_results", "mutation_summary.tsv"),
                ("fasta_qc_results", "FASTA_QC_summary.html"),
                ("mge_results", "staphscope_mge_summary.html"),
            ]
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

            for subdir, filename in required_files:
                if subdir:
                    src = output_dir / subdir / filename
                else:
                    src = output_dir / filename
                if src.exists():
                    shutil.copy2(src, temp_dir / filename)
                    self.logger.info(f"Copied {filename} to sample-centric temp dir")
                else:
                    self.logger.warning(f"Required file not found: {src} (skipping)")

            comprehensive_extensions = [".html", ".json", ".tsv"]
            for ext in comprehensive_extensions:
                src = output_dir / f"staphscope_comprehensive_report{ext}"
                if src.exists():
                    shutil.copy2(src, temp_dir / f"staphscope_comprehensive_report{ext}")
                    self.logger.info(f"Copied staphscope_comprehensive_report{ext} to sample-centric temp dir")
                else:
                    self.logger.warning(f"staphscope_comprehensive_report{ext} not found; sample-centric reporter may fail")

            self.banner.display_info("Running sample-centric reporter...")
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
                self.logger.error(f"Sample-centric reporter failed with return code {result.returncode}")
                return False

            src_dir = temp_dir / "STAPHSCOPE_ULTIMATE_SAMPLE_CENTRIC_REPORTS"
            if src_dir.exists():
                dst_dir = self.user_output_dir / "STAPHSCOPE_ULTIMATE_SAMPLE_CENTRIC_REPORTS"
                if dst_dir.exists():
                    shutil.rmtree(dst_dir)
                shutil.copytree(src_dir, dst_dir)
                self.logger.info(f"Sample-centric reports copied to {dst_dir}")
                self.banner.display_success("✅ Sample-centric reporter completed successfully.")
                return True
            else:
                self.logger.error("STAPHSCOPE_ULTIMATE_SAMPLE_CENTRIC_REPORTS not found after running sample-centric reporter")
                return False

        except Exception as e:
            self.logger.error(f"Sample-centric exception: {e}\n{traceback.format_exc()}")
            return False
        finally:
            if not self.keep_temp:
                shutil.rmtree(temp_dir, ignore_errors=True)

    # ------------------------------------------------------------------
    # Support modules
    # ------------------------------------------------------------------

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

    def run_visualization_analysis(self, output_dir: Path) -> bool:
        """Run the visualization module against the CSVs and master TSV it expects."""
        temp_dir = Path(tempfile.mkdtemp(prefix="staphscope_visualization_"))
        self._register_temp_dir(temp_dir)
        self.logger.info(f"Temporary directory for visualization: {temp_dir}")
        try:
            vis_module_orig = self.get_module_path("visualization_module")
            # Skip stale CSVs / TSVs / HTMLs that may live in the module folder
            shutil.copytree(
                vis_module_orig, temp_dir, dirs_exist_ok=True,
                ignore=shutil.ignore_patterns(
                    '*.csv', '*.tsv', '*.html', '*.json',
                    'STAPHSCOPE_VISUALIZATIONS', '__pycache__',
                ),
            )

            # The visualizer reads from two final-output subfolders
            final_report_dir = output_dir / "Staphscope_final_report"
            gene_centric_dir = final_report_dir / "STAPHSCOPE_ULTIMATE_GENE_CENTRIC_REPORTS"

            # (source_path, filename_to_copy_into_temp_dir)
            required_files = [
                # Master typing table
                (final_report_dir / "staphscope_comprehensive_report.tsv",
                 "staphscope_comprehensive_report.tsv"),
                # Gene tables exported by the gene-centric module
                (gene_centric_dir / "amr_genes.csv",         "amr_genes.csv"),
                (gene_centric_dir / "virulence_genes.csv",   "virulence_genes.csv"),
                (gene_centric_dir / "bacmet_genes.csv",      "bacmet_genes.csv"),
                (gene_centric_dir / "plasmid_replicons.csv", "plasmid_replicons.csv"),
                (gene_centric_dir / "mutations.csv",         "mutations.csv"),
                (gene_centric_dir / "mge_profile.csv",       "mge_profile.csv"),
                (gene_centric_dir / "fasta_qc.csv",          "fasta_qc.csv"),
            ]

            copied_count = 0
            missing = []
            for src, dest_name in required_files:
                if src.exists():
                    shutil.copy2(src, temp_dir / dest_name)
                    copied_count += 1
                    self.logger.info(f"Copied {dest_name} → visualization temp dir")
                else:
                    missing.append(src.name)
                    self.logger.warning(f"Required file not found: {src}")

            self.banner.display_info(
                f"Copied {copied_count} file(s) to temporary visualization module"
            )
            if missing:
                self.banner.display_warning(
                    f"Missing files (visualizer will run degraded): {', '.join(missing)}"
                )

            vis_script = temp_dir / "staphscope_visualizer.py"
            if not vis_script.exists():
                self.banner.display_error(f"Visualization script not found at: {vis_script}")
                return False

            self.banner.display_info("Running visualization module...")
            cmd = [sys.executable, str(vis_script), "-i", ".", "-o", "STAPHSCOPE_VISUALIZATIONS"]
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
                image_files = [f for f in vis_files
                               if f.suffix in ('.png', '.jpg', '.jpeg', '.svg')]
                zip_files = [f for f in vis_files if f.suffix == '.zip']
                self.banner.display_success(f"✅ Visualizations copied to: {target_vis_dir}")
                self.banner.display_info(f"   📊 {len(html_files)} HTML report(s)")
                self.banner.display_info(f"   🖼️  {len(image_files)} visualization image(s)")
                if zip_files:
                    self.banner.display_info(f"   📦 {len(zip_files)} export bundle(s)")
            return True

        except Exception as e:
            self.logger.error(f"Visualization exception: {e}\n{traceback.format_exc()}")
            return False
        finally:
            if not self.keep_temp:
                shutil.rmtree(temp_dir, ignore_errors=True)
                self.logger.info(f"Removed temporary directory: {temp_dir}")

    def copy_summary_results_to_final_directory(self, output_dir: Path):
        try:
            self.banner.display_info("Copying summary results to final directory...")
            final_report_dir = output_dir / "Staphscope_final_report"
            final_report_dir.mkdir(parents=True, exist_ok=True)

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

            ultimate_reports_dir = output_dir / "STAPHSCOPE_ULTIMATE_GENE_CENTRIC_REPORTS"
            if ultimate_reports_dir.exists() and ultimate_reports_dir.is_dir():
                target_ultimate_dir = final_report_dir / "STAPHSCOPE_ULTIMATE_GENE_CENTRIC_REPORTS"
                if target_ultimate_dir.exists():
                    shutil.rmtree(target_ultimate_dir)
                shutil.copytree(ultimate_reports_dir, target_ultimate_dir)
                self.banner.display_info("  ✓ Copied: STAPHSCOPE_ULTIMATE_GENE_CENTRIC_REPORTS directory")

            sample_centric_dir = output_dir / "STAPHSCOPE_ULTIMATE_SAMPLE_CENTRIC_REPORTS"
            if sample_centric_dir.exists() and sample_centric_dir.is_dir():
                target_sample_dir = final_report_dir / "STAPHSCOPE_ULTIMATE_SAMPLE_CENTRIC_REPORTS"
                if target_sample_dir.exists():
                    shutil.rmtree(target_sample_dir)
                shutil.copytree(sample_centric_dir, target_sample_dir)
                self.banner.display_info("  ✓ Copied: STAPHSCOPE_ULTIMATE_SAMPLE_CENTRIC_REPORTS directory")

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

    # ------------------------------------------------------------------
    # Helpers
    # ------------------------------------------------------------------

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

        if self.logger is None:
            import logging
            logging.basicConfig(level=logging.INFO, format='%(message)s')
            self.logger = logging.getLogger("StaphScope")

        self.banner.display_info("Updating AMRFinderPlus database...")
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

    def run_sequential_analyses(self, fasta_files: List[Path], output_dir: Path, threads: int,
                                skip_modules: Dict[str, bool],
                                amr_min_identity: float, amr_min_coverage: float,
                                amr_skip_mutations: bool, amr_force_update: bool,
                                abricate_min_identity: int = 80,
                                abricate_min_coverage: int = 80) -> Dict[str, bool]:
        analysis_functions = [
            ("FASTA QC", self.run_fasta_qc_analysis,
             "FASTA QC Analysis", "Sequence Quality Control & Statistics",
             not skip_modules.get('fasta_qc', False)),
            ("MLST", self.run_mlst_analysis,
             "MLST Analysis", "Multi-Locus Sequence Typing",
             not skip_modules.get('mlst', False)),
            ("spa typing", self.run_spa_typing,
             "SPA TYPING ANALYSIS", "Staphylococcal Protein A Typing",
             not skip_modules.get('spa', False)),
            ("SCCmec (CGE)", self.run_sccmec_analysis,
             "SCCMEC (CGE) ANALYSIS", "Methicillin Resistance Cassette Typing (CGE)",
             not skip_modules.get('sccmec', False)),
            ("SCCmec (RPet)", self.run_sccmec_rpet_analysis,
             "SCCMEC (RPET) ANALYSIS", "Methicillin Resistance Cassette Typing (RPet)",
             not skip_modules.get('sccmec_rpet', False)),
            ("Capsule", self.run_capsule_analysis,
             "CAPSULE TYPING", "Capsular Polysaccharide Typing",
             not skip_modules.get('capsule', False)),
            ("Agr", self.run_agr_analysis,
             "AGR TYPING", "Agr Accessory Gene Typing",
             not skip_modules.get('agr', False)),
            ("AMRFinderPlus",
             lambda f, o, t: self.run_amrfinder_analysis(
                 f, o, t, amr_min_identity, amr_min_coverage, amr_skip_mutations, amr_force_update),
             "AMR ANALYSIS", "Antimicrobial Resistance Gene Detection",
             not skip_modules.get('amr', False)),
            ("ABRicate",
             lambda f, o, t: self.run_abricate_analysis(
                 f, o, t, abricate_min_identity, abricate_min_coverage),
             "ABRICATE ANALYSIS", "Comprehensive Resistance, Plasmid & Virulence Gene Screening",
             not skip_modules.get('abricate', False)),
            ("MGE", self.run_mge_analysis,
             "MGE ANALYSIS", "Mobile Genetic Element Profiling",
             not skip_modules.get('mge', False)),
        ]
        active = [(name, func, header, desc)
                  for name, func, header, desc, enabled in analysis_functions if enabled]
        if not active:
            self.banner.display_warning("All analyses were skipped! Nothing to run.")
            return {}
        self.banner.display_info(f"Running {len(active)} analyses")
        results = {}
        for name, func, header, desc in active:
            self.banner.display_module_header(header, desc)
            try:
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

            subdirs = ["fasta_qc_results", "mlst_results", "spa_results", "sccmec_cge_results",
                       "sccmec_rpet_results", "capsule_results", "abricate_results",
                       "staph_amrfinder_results", "lineage_results", "agr_results",
                       "mge_results"]
            for subdir in subdirs:
                (output_path / subdir).mkdir(exist_ok=True)

            self.banner.display_module_header("Analysis Plan", "Modules to be executed")
            analyses_to_run = [
                ("FASTA QC", not skip_modules.get('fasta_qc', False)),
                ("MLST", not skip_modules.get('mlst', False)),
                ("spa typing", not skip_modules.get('spa', False)),
                ("SCCmec (CGE)", not skip_modules.get('sccmec', False)),
                ("SCCmec (RPet)", not skip_modules.get('sccmec_rpet', False)),
                ("Capsule typing", not skip_modules.get('capsule', False)),
                ("Agr typing", not skip_modules.get('agr', False)),
                ("AMRFinderPlus", not skip_modules.get('amr', False)),
                ("ABRICATE", not skip_modules.get('abricate', False)),
                ("MGE profiling", not skip_modules.get('mge', False)),
                ("Lineage Reference", not skip_modules.get('lineage', False)),
                ("Comprehensive Report", not skip_comprehensive),
                ("Ultimate Reporter (gene-centric)", not skip_comprehensive),
                ("Sample-Centric Reporter", not skip_sample_centric),
                ("Visualization", not skip_visualization),
            ]
            for analysis, enabled in analyses_to_run:
                status = "✅ ENABLED" if enabled else "⏸️  SKIPPED"
                print(f"   {status} - {analysis}")
            sys.stdout.flush()

            analysis_results = self.run_sequential_analyses(
                fasta_files, output_path, threads, skip_modules,
                amr_min_identity, amr_min_coverage,
                amr_skip_mutations, amr_force_update,
                abricate_min_identity, abricate_min_coverage)

            if not skip_modules.get('lineage', False):
                self.banner.display_module_header("Lineage Database",
                                                  "S. aureus Lineage Reference Generation")
                lineage_success = self.run_lineage_analysis(output_path)
                analysis_results["Lineage Reference"] = lineage_success
                print()

            if not skip_comprehensive:
                self.banner.display_module_header(
                    "Comprehensive & Ultimate Reports",
                    "Unified MLST, spa, SCCmec, agr, capsule and gene-centric integration")
                summary_success = self.run_comprehensive_and_ultimate_reports(output_path)
                analysis_results["Comprehensive & Ultimate Reports"] = summary_success
                if not summary_success:
                    self.banner.display_warning("Summary reports had issues")
                print()

            if not skip_sample_centric:
                if skip_comprehensive:
                    self.banner.display_warning(
                        "Sample-centric reporter requires staphscope_comprehensive_report.html.")
                    self.banner.display_warning(
                        "Since --skip-comprehensive was used, sample-centric reporter will be skipped.")
                    analysis_results["Sample-Centric Reporter"] = False
                else:
                    self.banner.display_module_header("Sample-Centric Reporter",
                                                      "Interactive isolate-centric report")
                    sample_success = self.run_sample_centric_analysis(output_path)
                    analysis_results["Sample-Centric Reporter"] = sample_success
                    print()

            self.copy_summary_results_to_final_directory(output_path)

            for dup in ["staphscope_comprehensive_report.html",
                        "staphscope_comprehensive_report.json",
                        "staphscope_comprehensive_report.tsv"]:
                (output_path / dup).unlink(missing_ok=True)
            shutil.rmtree(output_path / "STAPHSCOPE_ULTIMATE_GENE_CENTRIC_REPORTS", ignore_errors=True)
            shutil.rmtree(output_path / "STAPHSCOPE_ULTIMATE_SAMPLE_CENTRIC_REPORTS", ignore_errors=True)

            if not skip_visualization:
                self.banner.display_module_header("Visualization",
                                                  "Interactive Visualizations & Dashboard")
                visualization_success = self.run_visualization_analysis(output_path)
                analysis_results["Visualization"] = visualization_success
                print()

            analysis_time = datetime.now() - start_time
            analysis_time_str = str(analysis_time).split('.')[0]
            successful_count = sum(analysis_results.values())
            total_count = len(analysis_results)

            self.banner.display_footer(analysis_time=analysis_time_str,
                                       samples_processed=len(fasta_files))

            if successful_count == total_count:
                self.banner.display_success(f"🎉 All {total_count} analyses completed successfully!")
                self.banner.display_success("🧹 All module directories have been cleaned up")
                print("\n📁 FINAL OUTPUT STRUCTURE:")
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
                self.banner.display_warning(
                    f"⚠️  {successful_count}/{total_count} analyses completed successfully.")

            print("\n📚 Please cite our StaphScope paper:")
            print("   Beckley, B., Amarh, V. StaphScope: a species-optimized computational pipeline for rapid and accessible Staphylococcus aureus genotyping and surveillance. BMC Genomics (2026).")
            print("   https://doi.org/10.1186/s12864-026-12609-x")

        except KeyboardInterrupt:
            self.banner.display_error("Analysis interrupted by user")
        except Exception as e:
            self.banner.display_error(f"Critical error in analysis pipeline: {str(e)}")
            import traceback
            traceback.print_exc()

    # ------------------------------------------------------------------
    # MLST Database Management
    # ------------------------------------------------------------------

    def _get_system_module_path(self, module_name: str) -> Path:
        """Return the original system module path (without fallback)."""
        if hasattr(sys, 'prefix'):
            share_path = Path(sys.prefix) / "share" / "staphscope" / "modules" / module_name
            if share_path.exists():
                return share_path
        return self.base_dir / "modules" / module_name

    def pull_mlst_database(self) -> bool:
        """Download the prebuilt S. aureus MLST database from GitHub.

        Writes into the packaged module directory if writable (user-local conda),
        otherwise falls back to ~/.local/share/staphscope/mlst_module/ (Docker,
        Apptainer, shared conda). The user-local path is exposed via
        STAPHSCOPE_MLST_MODULE_PATH and STAPHSCOPE_MLST_DB so subsequent calls
        in this process use it, and future processes pick it up automatically
        because get_mlst_db_dir() prefers a populated user-local db/.
        """
        module_path = self.get_module_path("mlst_module")

        if os.access(module_path, os.W_OK):
            target_path = module_path
            self.banner.display_info(f"Using writable module directory: {target_path}")
        else:
            target_path = Path(os.environ.get(
                "STAPHSCOPE_MLST_MODULE_PATH",
                Path.home() / ".local" / "share" / "staphscope" / "mlst_module",
            ))
            target_path.mkdir(parents=True, exist_ok=True)
            self.banner.display_info(
                f"System MLST module is read-only; using user-local path: {target_path}"
            )

        system_path = self._get_system_module_path("mlst_module")

        repo_url = "https://github.com/bbeckley-hub/mlst/archive/refs/heads/master.zip"
        temp_zip = Path(tempfile.mktemp(suffix=".zip"))
        extract_dir = None

        try:
            self.banner.display_info("Downloading MLST data from GitHub...")
            urllib.request.urlretrieve(repo_url, temp_zip)

            extract_dir = Path(tempfile.mkdtemp(prefix="staphscope_mlst_extract_"))
            with zipfile.ZipFile(temp_zip, 'r') as zip_ref:
                zip_ref.extractall(extract_dir)

            src_root = extract_dir / "mlst-master"
            if not src_root.exists():
                self.banner.display_error("Extracted repository root not found.")
                return False

            for dirname in ['bin', 'db', 'perl5', 'scripts']:
                src = src_root / dirname
                if src.exists():
                    dst = target_path / dirname
                    if dst.exists():
                        shutil.rmtree(dst)
                    shutil.copytree(src, dst)
                    self.logger.info(f"Copied {dirname} to {dst}")
                else:
                    self.banner.display_warning(f"'{dirname}' not found in repository, skipping.")

            # Copy helper scripts when writing to a user-local target
            if target_path != system_path:
                for script in ['mlst_module.py', 'mlst_database.sh']:
                    src_script = system_path / script
                    if src_script.exists():
                        shutil.copy2(src_script, target_path / script)
                        self.logger.info(f"Copied {script} to {target_path / script}")

            # Verify the scheme landed
            db_check = target_path / "db" / "pubmlst" / "saureus"
            if not db_check.exists():
                self.banner.display_error(f"MLST scheme missing after pull: {db_check}")
                return False

            self.banner.display_success(f"MLST module pulled to {target_path}")
            os.environ['STAPHSCOPE_MLST_MODULE_PATH'] = str(target_path)
            os.environ['STAPHSCOPE_MLST_DB'] = str(target_path / "db")
            return True

        except Exception as e:
            self.banner.display_error(f"Failed to download/extract MLST data: {e}")
            return False
        finally:
            if temp_zip.exists():
                temp_zip.unlink()
            if extract_dir and extract_dir.exists():
                shutil.rmtree(extract_dir, ignore_errors=True)

    def update_mlst_database(self) -> bool:
        """Refresh the S. aureus MLST scheme from PubMLST (API key required).

        Targets the same module tree that --pull-mlst-db populates:
          - packaged module dir if writable, otherwise
          - ~/.local/share/staphscope/mlst_module/

        Does NOT pull scripts or seed metadata — it only updates
        <db_dir>/pubmlst/saureus/ and rebuilds <db_dir>/blast/.
        Run --pull-mlst-db first on a fresh install.
        """
        module_path = self.get_module_path("mlst_module")
        if os.access(module_path, os.W_OK):
            target_module = module_path
        else:
            target_module = Path(os.environ.get(
                "STAPHSCOPE_MLST_MODULE_PATH",
                Path.home() / ".local" / "share" / "staphscope" / "mlst_module",
            ))

        db_dir = target_module / "db"
        script_dir = target_module / "bin"
        scheme_map = db_dir / "scheme_species_map.tab"

        if not scheme_map.exists():
            self.banner.display_error(
                f"MLST module tree is not populated at {target_module}"
            )
            self.banner.display_info(
                "Run `staphscope --pull-mlst-db` first to install the module, "
                "then run `staphscope --update-mlst-db` to refresh the scheme."
            )
            return False

        if not script_dir.exists():
            self.banner.display_error(f"MLST bin/ not found at {script_dir}")
            return False

        os.environ['STAPHSCOPE_MLST_MODULE_PATH'] = str(target_module)
        os.environ['STAPHSCOPE_MLST_DB'] = str(db_dir)

        try:
            from staphscope.modules.mlst_module.mlst_module import update_mlst_database as _update
        except ImportError:
            sys.path.insert(0, str(self.base_dir / "modules" / "mlst_module"))
            from mlst_module import update_mlst_database as _update

        self.banner.display_info(f"Using MLST database directory: {db_dir}")
        return _update(db_dir, script_dir)


def main():
    parser = argparse.ArgumentParser(
        description="StaphScope: Advanced Staphylococcus aureus Typing & Lineage Analysis Platform",
        formatter_class=ColoredHelpFormatter,
        epilog=f"""
{ColoredHelpFormatter.BOLD}{ColoredHelpFormatter.GREEN}First-Time Setup (Non-Docker Users){ColoredHelpFormatter.RESET}
  Before your first analysis, run these recommended commands to prepare the databases:

  1. AMR database (AMRFinderPlus):
       staphscope --update-amr-db
     (or use --force-update-amr-db to overwrite an existing database)

  2. MLST database (choose one):
       # Quick start – download pre-built S. aureus database from GitHub (no API key)
       staphscope --pull-mlst-db

       # Advanced – fetch the latest S. aureus data from PubMLST (requires API key)
       # First, set up your PubMLST API key (see instructions below), then:
       staphscope --update-mlst-db

  3. ABRicate databases (if you plan to use the ABRicate module):
       abricate --setupdb

{ColoredHelpFormatter.BOLD}{ColoredHelpFormatter.GREEN}Examples:{ColoredHelpFormatter.RESET}
  # Basic single-sample analysis
  staphscope -i genome.fna -o results/

  # Batch analysis with glob pattern
  staphscope -i "*.fna" -o batch_results --threads 8

  # Skip selected modules
  staphscope -i "*.fasta" -o analysis --threads 16 --skip-lineage --skip-visualization

  # AMR with custom thresholds and no mutation reporting
  staphscope -i "*.fna" -o results --amr-min-identity 0.95 --amr-min-coverage 0.9 --skip-amr-mutations

  # Force update AMR database before analysis
  staphscope -i "*.fna" -o results --amr-force-update

  # AMR database maintenance
  staphscope --update-amr-db                 # incremental update
  staphscope --force-update-amr-db           # complete overwrite

  # MLST database maintenance
  staphscope --pull-mlst-db                  # clone from GitHub (no credentials)
  staphscope --update-mlst-db                # update from PubMLST (API key required)

  # ABRicate with custom thresholds
  staphscope -i "*.fna" -o results --abricate-minid 85 --abricate-mincov 90

{ColoredHelpFormatter.BOLD}{ColoredHelpFormatter.GREEN}MLST Database Management:{ColoredHelpFormatter.RESET}
  --pull-mlst-db     : Download a pre-built S. aureus MLST database from GitHub.
                       Recommended for first-time users (no API key needed).
  --update-mlst-db   : Update the S. aureus MLST database directly from PubMLST.
                       Requires one-time API key setup (see instructions below).
                       Download may take 10-20 minutes depending on internet speed.

{ColoredHelpFormatter.BOLD}{ColoredHelpFormatter.GREEN}PubMLST API Key Setup (for --update-mlst-db):{ColoredHelpFormatter.RESET}
  1. Register/login at https://pubmlst.org/site-accounts
  2. Generate an API key from your profile (under "API keys").
  3. Run: mlstdb connect --db pubmlst --api-key
  4. Paste your API key when prompted.
  5. After setup, run staphscope --update-mlst-db whenever you want to refresh the database.

{ColoredHelpFormatter.BOLD}{ColoredHelpFormatter.GREEN}Supported FASTA formats:{ColoredHelpFormatter.RESET} .fna, .fasta

{ColoredHelpFormatter.BOLD}{ColoredHelpFormatter.GREEN}Analysis Modules:{ColoredHelpFormatter.RESET}
  • FASTA QC (Quality Control & Statistics)
  • MLST (Multi-Locus Sequence Typing)
  • spa typing (Staphylococcal Protein A)
  • SCCmec typing — CGE caller (Methicillin Resistance Cassette)
  • SCCmec typing — RPet caller (Robert Petit III's sccmec)
  • Capsule typing (cap5/cap8 serotype)
  • Agr typing (agrVATE)
  • AMR profiling (AMRFinderPlus) – point mutations reported by default
  • ABRICATE (Comprehensive resistance/plasmid/virulence)
  • MGE profiling (mobileOG-db)
  • Lineage reference database
  • Comprehensive report (MLST + spa + SCCmec + agr + capsule)
  • Ultimate reporter (Gene-centric integrated analysis)
  • Sample-centric reporter (Isolate-centric interactive report)
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

    parser.add_argument('-i', '--input',
                        help='Input FASTA file(s) - can use glob patterns like "*.fna" or "*.fasta"')
    parser.add_argument('-o', '--output', help='Output directory for all results')
    parser.add_argument('-t', '--threads', type=int, default=2,
                        help='Number of threads (default: 2)')

    parser.add_argument('--skip-fasta-qc', action='store_true', help='Skip FASTA QC analysis')
    parser.add_argument('--skip-mlst', action='store_true', help='Skip MLST analysis')
    parser.add_argument('--skip-spa', action='store_true', help='Skip spa typing analysis')
    parser.add_argument('--skip-sccmec', action='store_true', help='Skip SCCmec CGE analysis')
    parser.add_argument('--skip-sccmec-rpet', action='store_true', help='Skip SCCmec RPet typing analysis')
    parser.add_argument('--skip-capsule', action='store_true', help='Skip capsule typing analysis')
    parser.add_argument('--skip-agr', action='store_true', help='Skip agr typing analysis')
    parser.add_argument('--skip-amr', action='store_true', help='Skip AMR analysis (AMRFinderPlus)')
    parser.add_argument('--skip-abricate', action='store_true', help='Skip ABRicate analysis')
    parser.add_argument('--skip-mge', action='store_true', help='Skip MGE (mobileOG) profiling')
    parser.add_argument('--skip-lineage', action='store_true', help='Skip lineage reference generation')
    parser.add_argument('--skip-comprehensive', action='store_true',
                        help='Skip comprehensive report AND ultimate reporter (gene-centric)')
    parser.add_argument('--skip-sample-centric', action='store_true',
                        help='Skip sample-centric reporter (isolate-centric)')
    parser.add_argument('--skip-visualization', action='store_true',
                        help='Skip visualization module (dashboards & plots)')

    parser.add_argument('--amr-min-identity', type=float,
                        help='Minimum identity for AMR hits (0..1)')
    parser.add_argument('--amr-min-coverage', type=float,
                        help='Minimum coverage for AMR hits (0..1)')
    parser.add_argument('--skip-amr-mutations', action='store_true',
                        help='Disable point mutation reporting in AMR (enabled by default)')
    parser.add_argument('--amr-force-update', action='store_true',
                        help='Force update AMR database before analysis')
    parser.add_argument('--update-amr-db', action='store_true',
                        help='Update AMRFinderPlus database (incremental) and exit')
    parser.add_argument('--force-update-amr-db', action='store_true',
                        help='Force complete AMR database update (overwrites old) and exit')

    parser.add_argument('--pull-mlst-db', action='store_true',
                        help='Pull/refresh MLST database from GitHub fork (Recommended for first-time setup)')
    parser.add_argument('--update-mlst-db', action='store_true',
                        help='Update S. aureus MLST database from PubMLST (API key required)')

    parser.add_argument('--abricate-minid', type=int, default=80,
                        help='Minimum identity for ABRicate hits (default: 80)')
    parser.add_argument('--abricate-mincov', type=int, default=80,
                        help='Minimum coverage for ABRicate hits (default: 80)')

    parser.add_argument('--keep-temp', action='store_true',
                        help='Do not delete temporary directories (for debugging)')
    parser.add_argument('--clean-output', action='store_true',
                        help='Delete output directory before analysis (prevents mixing results from different runs)')

    parser.add_argument('-v', '--version', action='version', version=f'StaphScope version {__version__}')

    args = parser.parse_args()

    if args.update_amr_db or args.force_update_amr_db:
        orch = StaphScopeOrchestrator()
        if args.force_update_amr_db:
            orch.update_amr_database(force=True)
        else:
            orch.update_amr_database(force=False)
        sys.exit(0)

    if args.pull_mlst_db or args.update_mlst_db:
        orch = StaphScopeOrchestrator()
        temp_log_dir = Path(tempfile.mkdtemp(prefix="staphscope_mlst_db_"))
        orch.setup_logging(temp_log_dir)
        if args.pull_mlst_db:
            success = orch.pull_mlst_database()
        else:
            success = orch.update_mlst_database()
        if not orch.keep_temp:
            shutil.rmtree(temp_log_dir, ignore_errors=True)
        sys.exit(0 if success else 1)

    if not args.input or not args.output:
        parser.error("When not using --update-amr-db, --force-update-amr-db, "
                     "--pull-mlst-db, or --update-mlst-db, both -i/--input and -o/--output are required.")

    skip_modules = {
        'fasta_qc': args.skip_fasta_qc,
        'mlst': args.skip_mlst,
        'spa': args.skip_spa,
        'sccmec': args.skip_sccmec,
        'sccmec_rpet': args.skip_sccmec_rpet,
        'capsule': args.skip_capsule,
        'agr': args.skip_agr,
        'amr': args.skip_amr,
        'abricate': args.skip_abricate,
        'mge': args.skip_mge,
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
        abricate_min_coverage=args.abricate_mincov,
    )


if __name__ == "__main__":
    main()
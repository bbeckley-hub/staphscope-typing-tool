#!/usr/bin/env python3
"""
StaphScope Main Orchestrator - UPDATED WITH FASTA QC AND VISUALIZATION MODULES
Complete S. aureus typing pipeline - FASTA QC, MLST, spa, SCCmec, AMR, Virulence, Lineage, Visualization
Author: Brown Beckley <brownbeckley94@gmail.com>
Date: 2025
Send a quick mail for any issues or further explanations.
Affiliation: University of Ghana Medical School-Department of Medical Biochemistry
version 1.1.0
"""

import os
import sys
import glob
import argparse
import subprocess
import shutil
from pathlib import Path
from datetime import datetime
from typing import Dict, List, Set

# Import banner
try:
    from .core.banner import StaphScopeBanner
except (ImportError, SystemError):
    sys.path.insert(0, str(Path(__file__).parent))
    from core.banner import StaphScopeBanner

class StaphScopeOrchestrator:
    """Final StaphScope orchestrator - with FASTA QC first and Visualization last"""
    
    def __init__(self):
        self.banner = StaphScopeBanner()
        self.base_dir = Path(__file__).parent
        
        # Dictionary for visualization module HTML files (from Brown's directory listing)
        self.visualization_html_files = {
            "amrfinder": "staph_amrfinder_summary_report.html",
            "argannot": "staph_argannot_summary_report.html",
            "card": "staph_card_summary_report.html",
            "megares": "staph_megares_summary_report.html",
            "ncbi": "staph_ncbi_summary_report.html",
            "plasmidfinder": "staph_plasmidfinder_summary_report.html",
            "resfinder": "staph_resfinder_summary_report.html",
            "comprehensive": "staphscope_comprehensive_report.html",
            "vfdb": "staph_vfdb_summary_report.html"
        }
    
    def find_fasta_files(self, input_path: str) -> List[Path]:
        """Find all FASTA files using glob patterns"""
        self.banner.display_info(f"Searching for files with pattern: {input_path}")
        
        # Handle quoted wildcards properly
        if '*' in input_path or '?' in input_path:
            matched_files = glob.glob(input_path)
            fasta_files = [Path(f) for f in matched_files if Path(f).is_file() and 
                          f.lower().endswith(('.fna', '.fasta', '.fa', '.fn', '.faa', '.gb', '.gbk', '.gbff')) and
                          not Path(f).name.startswith('.')]
            self.banner.display_success(f"Found {len(fasta_files)} FASTA files")
            return sorted(fasta_files)
        
        # Handle direct file path
        input_path_obj = Path(input_path)
        if input_path_obj.is_file() and input_path_obj.suffix.lower() in ['.fna', '.fasta', '.fa', '.fn', '.faa', '.gb', '.gbk', '.gbff']:
            self.banner.display_success(f"Found single FASTA file: {input_path_obj.name}")
            return [input_path_obj]
        
        # Handle directory
        if input_path_obj.is_dir():
            patterns = [
                f"{input_path}/*.fna", f"{input_path}/*.fasta",
                f"{input_path}/*.fa", f"{input_path}/*.fn",
                f"{input_path}/*.faa", f"{input_path}/*.gb",
                f"{input_path}/*.gbk", f"{input_path}/*.gbff"
            ]
            fasta_files = []
            for pattern in patterns:
                matched_files = glob.glob(pattern)
                for file_path in matched_files:
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
        """Get the correct file pattern based on actual file extensions"""
        if not fasta_files:
            return '"*.fna"'  # Default fallback
        
        # Get all unique extensions from the input files
        extensions = set(f.suffix.lower() for f in fasta_files)
        
        # If all files have the same extension, use that
        if len(extensions) == 1:
            ext = list(extensions)[0]
            return f'"*{ext}"'
        
        # If mixed extensions, use a pattern that matches all FASTA files
        return '"*"'

    def cleanup_module_directory(self, module_path: Path, fasta_files: List[Path]):
        """COMPREHENSIVE cleanup of module directory after analysis"""
        try:
            self.banner.display_info(f"Cleaning up {module_path.name}...")
            
            # 1. Remove copied input files
            for fasta_file in fasta_files:
                temp_file = module_path / fasta_file.name
                if temp_file.exists():
                    temp_file.unlink()
            
            # 2. Remove common output directories
            output_dirs = [
                "mlst_results", "spa_results", "abricate_results", 
                "staph_amrfinder_results", "abricate_results", "fasta_qc_results"
            ]
            for output_dir in output_dirs:
                dir_path = module_path / output_dir
                if dir_path.exists():
                    shutil.rmtree(dir_path)
            
            # 3. Remove SCCmec-specific directories and files
            if module_path.name == "sccmec_module":
                # Remove s_* directories
                for s_dir in module_path.glob("s_*"):
                    if s_dir.is_dir():
                        shutil.rmtree(s_dir)
                # Remove summary files
                for summary_file in ["staphscope_summary.html", "staphscope_summary.tsv"]:
                    summary_path = module_path / summary_file
                    if summary_path.exists():
                        summary_path.unlink()
            
            # 4. Remove any other common temporary files
            temp_patterns = ["*.txt", "*.log", "*.tmp", "temp_*"]
            for pattern in temp_patterns:
                for temp_file in module_path.glob(pattern):
                    if temp_file.is_file():
                        temp_file.unlink()
            
            self.banner.display_success(f"✅ {module_path.name} cleaned up successfully")
            
        except Exception as e:
            self.banner.display_warning(f"⚠️  Partial cleanup issue in {module_path.name}: {str(e)}")

    def run_fasta_qc_analysis(self, fasta_files: List[Path], output_dir: Path, threads: int) -> bool:
        """Run FASTA QC analysis - first step in pipeline"""
        qc_module_path = self.base_dir / "modules" / "fasta_qc_module"
        
        try:
            qc_script = qc_module_path / "staph_fasta_qc.py"
            
            if not qc_script.exists():
                self.banner.display_error(f"FASTA QC script not found at: {qc_script}")
                return False
            
            # Copy files to QC module directory
            for fasta_file in fasta_files:
                target_file = qc_module_path / fasta_file.name
                shutil.copy2(fasta_file, target_file)
            
            self.banner.display_info(f"Copied {len(fasta_files)} files to FASTA QC module")
            
            # Get correct file pattern based on actual files
            file_pattern = self.get_file_pattern(fasta_files)
            
            # Build command - using the correct pattern
            cmd = [
                sys.executable, str(qc_script),
                file_pattern,
                "-o", "fasta_qc_results",
                "-c", str(threads)
            ]
            
            self.banner.display_info(f"Running FASTA QC analysis with pattern: {file_pattern}")
            result = subprocess.run(" ".join(cmd), shell=True, capture_output=True, text=True, cwd=qc_module_path)
            
            if result.returncode == 0:
                self.banner.display_success("FASTA QC analysis completed!")
                
                # Copy results to output directory
                qc_source = qc_module_path / "fasta_qc_results"
                qc_target = output_dir / "fasta_qc_results"
                
                if qc_source.exists():
                    if qc_target.exists():
                        shutil.rmtree(qc_target)
                    shutil.copytree(qc_source, qc_target)
                    self.banner.display_success(f"FASTA QC results copied to: {qc_target}")
                
                return True
            else:
                self.banner.display_warning("FASTA QC analysis had warnings/errors")
                if result.stderr:
                    error_lines = result.stderr.strip().split('\n')
                    for line in error_lines[:5]:  # Show first 5 errors
                        if line.strip():
                            print(f"  ❌ {line.strip()}")
                # Still return True to continue pipeline unless critical error
                return True
                
        except Exception as e:
            self.banner.display_error(f"FASTA QC analysis failed: {str(e)}")
            # Continue pipeline even if QC fails
            return True
        finally:
            # ALWAYS cleanup, even if analysis fails
            self.banner.display_info("Cleaning up fasta_qc_module...")
            self.cleanup_module_directory(qc_module_path, fasta_files)
            self.banner.display_success("✅ fasta_qc_module cleaned up successfully")

    def run_mlst_analysis(self, fasta_files: List[Path], output_dir: Path, threads: int) -> bool:
        """Run MLST analysis"""
        mlst_module_path = self.base_dir / "modules" / "mlst_module"
        
        try:
            mlst_script = mlst_module_path / "mlst_module.py"
            
            if not mlst_script.exists():
                self.banner.display_error(f"MLST script not found at: {mlst_script}")
                return False
            
            # Copy files to MLST module directory
            for fasta_file in fasta_files:
                target_file = mlst_module_path / fasta_file.name
                shutil.copy2(fasta_file, target_file)
            
            self.banner.display_info(f"Copied {len(fasta_files)} files to MLST module")
            
            # Get correct file pattern based on actual files
            file_pattern = self.get_file_pattern(fasta_files)
            
            # Build command
            cmd = [
                sys.executable, str(mlst_script),
                "-i", file_pattern,
                "-o", "mlst_results",
                "-db", "db",
                "-sc", "bin", 
                "--batch"
            ]
            
            self.banner.display_info(f"Running MLST analysis with pattern: {file_pattern}")
            result = subprocess.run(" ".join(cmd), shell=True, capture_output=True, text=True, cwd=mlst_module_path)
            
            if result.returncode == 0:
                self.banner.display_success("MLST analysis completed!")
                
                # Copy results to output directory
                mlst_source = mlst_module_path / "mlst_results"
                mlst_target = output_dir / "mlst_results"
                
                if mlst_source.exists():
                    if mlst_target.exists():
                        shutil.rmtree(mlst_target)
                    shutil.copytree(mlst_source, mlst_target)
                    self.banner.display_success(f"MLST results copied to: {mlst_target}")
                
                return True
            else:
                self.banner.display_warning("MLST analysis had warnings")
                if result.stderr:
                    error_lines = result.stderr.strip().split('\n')
                    for line in error_lines[:3]:
                        if line.strip():
                            print(f"  ⚠️  {line.strip()}")
                return True
                
        except Exception as e:
            self.banner.display_error(f"MLST analysis failed: {str(e)}")
            return False
        finally:
            # ALWAYS cleanup, even if analysis fails
            self.banner.display_info("Cleaning up mlst_module...")
            self.cleanup_module_directory(mlst_module_path, fasta_files)
            self.banner.display_success("✅ mlst_module cleaned up successfully")

    def run_spa_typing(self, fasta_files: List[Path], output_dir: Path, threads: int) -> bool:
        """Run spa typing"""
        spa_module_path = self.base_dir / "modules" / "spa_module" / "spatyper" / "spa_typing"
        
        try:
            spa_script = spa_module_path / "spa_typing_module.py"
            
            if not spa_script.exists():
                self.banner.display_error(f"spa script not found at: {spa_script}")
                return False
            
            # Copy files to spa module directory
            for fasta_file in fasta_files:
                target_file = spa_module_path / fasta_file.name
                shutil.copy2(fasta_file, target_file)
            
            self.banner.display_info(f"Copied {len(fasta_files)} files to spa module")
            
            # Set Python path for spa typing dependencies
            spa_env = os.environ.copy()
            spa_pythonpath = str(self.base_dir / "modules" / "spa_module" / "spatyper")
            current_pythonpath = spa_env.get('PYTHONPATH', '')
            spa_env['PYTHONPATH'] = f"{spa_pythonpath}:{current_pythonpath}" if current_pythonpath else spa_pythonpath
            
            # Build command
            cmd = [
                sys.executable, str(spa_script),
                "-i", ".",
                "-o", "spa_results"
            ]
            
            self.banner.display_info("Running spa typing analysis...")
            result = subprocess.run(cmd, capture_output=True, text=True, cwd=spa_module_path, env=spa_env)
            
            if result.returncode == 0:
                self.banner.display_success("spa typing analysis completed!")
                
                # Copy results to output directory
                spa_source = spa_module_path / "spa_results"
                spa_target = output_dir / "spa_results"
                
                if spa_source.exists():
                    if spa_target.exists():
                        shutil.rmtree(spa_target)
                    shutil.copytree(spa_source, spa_target)
                    self.banner.display_success(f"spa results copied to: {spa_target}")
                
                return True
            else:
                self.banner.display_warning("spa typing analysis had warnings")
                if result.stderr:
                    error_lines = result.stderr.strip().split('\n')
                    for line in error_lines[:3]:
                        if line.strip():
                            print(f"  ⚠️  {line.strip()}")
                return True
                
        except Exception as e:
            self.banner.display_error(f"spa typing analysis failed: {str(e)}")
            return False
        finally:
            # ALWAYS cleanup, even if analysis fails
            self.banner.display_info("Cleaning up spa_typing...")
            self.cleanup_module_directory(spa_module_path, fasta_files)
            self.banner.display_success("✅ spa_typing cleaned up successfully")

    def run_sccmec_analysis(self, fasta_files: List[Path], output_dir: Path, threads: int) -> bool:
        """Run SCCmec analysis"""
        sccmec_module_path = self.base_dir / "modules" / "sccmec_module"
        
        try:
            sccmec_script = sccmec_module_path / "run_sccmec_batch.sh"
            summary_script = sccmec_module_path / "generate_staphscope_summary.sh"
            
            if not sccmec_script.exists():
                self.banner.display_error(f"SCCmec script not found at: {sccmec_script}")
                return False
            
            # Copy files to SCCmec module directory
            for fasta_file in fasta_files:
                target_file = sccmec_module_path / fasta_file.name
                shutil.copy2(fasta_file, target_file)
            
            self.banner.display_info(f"Copied {len(fasta_files)} files to SCCmec module")
            
            # Get correct file pattern based on actual files
            file_pattern = self.get_file_pattern(fasta_files)
            file_pattern_clean = file_pattern.strip('"')
            
            # Build command
            cmd = ["bash", str(sccmec_script), file_pattern_clean]
            
            self.banner.display_info(f"Running SCCmec analysis with pattern: {file_pattern}")
            result = subprocess.run(" ".join(cmd), shell=True, capture_output=True, text=True, cwd=sccmec_module_path)
            
            if result.returncode == 0:
                self.banner.display_success("SCCmec analysis completed!")
                
                # Run summary generation
                self.banner.display_info("Generating SCCmec summary...")
                summary_result = subprocess.run(
                    ["bash", str(summary_script)], 
                    capture_output=True, text=True, cwd=sccmec_module_path
                )
                
                if summary_result.returncode == 0:
                    self.banner.display_success("SCCmec summary generated!")
                else:
                    self.banner.display_warning("SCCmec summary generation had issues")
                
                # Create SCCmec output directory
                sccmec_output = output_dir / "sccmec_results"
                sccmec_output.mkdir(parents=True, exist_ok=True)
                
                # Copy all s_* directories
                s_dirs_copied = 0
                for s_dir in sccmec_module_path.glob("s_*"):
                    if s_dir.is_dir():
                        target_dir = sccmec_output / s_dir.name
                        if target_dir.exists():
                            shutil.rmtree(target_dir)
                        shutil.copytree(s_dir, target_dir)
                        s_dirs_copied += 1
                
                # Copy summary files
                summary_files = ["staphscope_summary.html", "staphscope_summary.tsv", "staphscope_detailed_results.csv"]
                for summary_file in summary_files:
                    source_file = sccmec_module_path / summary_file
                    if source_file.exists():
                        target_file = sccmec_output / summary_file
                        shutil.copy2(source_file, target_file)
                
                self.banner.display_success(f"SCCmec results copied to: {sccmec_output}")
                self.banner.display_info(f"Copied {s_dirs_copied} sample directories and summary files")
                
                return True
            else:
                self.banner.display_warning("SCCmec analysis had warnings")
                if result.stderr:
                    error_lines = result.stderr.strip().split('\n')
                    for line in error_lines[:3]:
                        if line.strip():
                            print(f"  ⚠️  {line.strip()}")
                return True
                
        except Exception as e:
            self.banner.display_error(f"SCCmec analysis failed: {str(e)}")
            return False
        finally:
            # ALWAYS cleanup, even if analysis fails
            self.banner.display_info("Cleaning up sccmec_module...")
            self.cleanup_module_directory(sccmec_module_path, fasta_files)
            self.banner.display_success("✅ sccmec_module cleaned up successfully")

    def run_amrfinder_analysis(self, fasta_files: List[Path], output_dir: Path, threads: int) -> bool:
        """Run AMRFinderPlus analysis"""
        amr_module_path = self.base_dir / "modules" / "amr_module"
        
        try:
            amr_script = amr_module_path / "amrfinder_standalone.py"
            
            if not amr_script.exists():
                self.banner.display_error(f"AMRFinderPlus script not found at: {amr_script}")
                return False
            
            # Copy files to AMR module directory
            for fasta_file in fasta_files:
                target_file = amr_module_path / fasta_file.name
                shutil.copy2(fasta_file, target_file)
            
            self.banner.display_info(f"Copied {len(fasta_files)} files to AMR module")
            
            # Get correct file pattern based on actual files
            file_pattern = self.get_file_pattern(fasta_files)
            
            # Build command
            cmd = [
                sys.executable, str(amr_script),
                file_pattern
            ]
            
            self.banner.display_info(f"Running AMRFinderPlus analysis with pattern: {file_pattern}")
            result = subprocess.run(" ".join(cmd), shell=True, capture_output=True, text=True, cwd=amr_module_path)
            
            if result.returncode == 0:
                self.banner.display_success("AMRFinderPlus analysis completed!")
                
                # Copy results to output directory
                amr_source = amr_module_path / "staph_amrfinder_results"
                amr_target = output_dir / "amr_results"
                
                if amr_source.exists():
                    if amr_target.exists():
                        shutil.rmtree(amr_target)
                    shutil.copytree(amr_source, amr_target)
                    self.banner.display_success(f"AMR results copied to: {amr_target}")
                
                return True
            else:
                self.banner.display_warning("AMRFinderPlus analysis had warnings")
                if result.stderr:
                    error_lines = result.stderr.strip().split('\n')
                    for line in error_lines[:3]:
                        if line.strip():
                            print(f"  ⚠️  {line.strip()}")
                return True
                
        except Exception as e:
            self.banner.display_error(f"AMRFinderPlus analysis failed: {str(e)}")
            return False
        finally:
            # ALWAYS cleanup, even if analysis fails
            self.banner.display_info("Cleaning up amr_module...")
            self.cleanup_module_directory(amr_module_path, fasta_files)
            self.banner.display_success("✅ amr_module cleaned up successfully")

    def run_abricate_analysis(self, fasta_files: List[Path], output_dir: Path, threads: int) -> bool:
        """Run Abricate analysis"""
        abricate_module_path = self.base_dir / "modules" / "abricate_module"
        
        try:
            abricate_script = abricate_module_path / "abricate_standalone.py"
            
            if not abricate_script.exists():
                self.banner.display_error(f"ABRicate script not found at: {abricate_script}")
                return False
            
            # Copy files to Abricate module directory
            for fasta_file in fasta_files:
                target_file = abricate_module_path / fasta_file.name
                shutil.copy2(fasta_file, target_file)
            
            self.banner.display_info(f"Copied {len(fasta_files)} files to ABRicate module")
            
            # Get correct file pattern based on actual files
            file_pattern = self.get_file_pattern(fasta_files)
            
            # Build command
            cmd = [
                sys.executable, str(abricate_script),
                file_pattern
            ]
            
            self.banner.display_info(f"Running ABRicate analysis with pattern: {file_pattern}")
            result = subprocess.run(" ".join(cmd), shell=True, capture_output=True, text=True, cwd=abricate_module_path)
            
            if result.returncode == 0:
                self.banner.display_success("ABRicate analysis completed!")
                
                # Copy results to output directory
                abricate_source = abricate_module_path / "abricate_results"
                abricate_target = output_dir / "abricate_results"
                
                if abricate_source.exists():
                    if abricate_target.exists():
                        shutil.rmtree(abricate_target)
                    shutil.copytree(abricate_source, abricate_target)
                    self.banner.display_success(f"ABRicate results copied to: {abricate_target}")
                
                return True
            else:
                self.banner.display_warning("ABRicate analysis had warnings")
                if result.stderr:
                    error_lines = result.stderr.strip().split('\n')
                    for line in error_lines[:3]:
                        if line.strip():
                            print(f"  ⚠️  {line.strip()}")
                return True
                
        except Exception as e:
            self.banner.display_error(f"ABRicate analysis failed: {str(e)}")
            return False
        finally:
            # ALWAYS cleanup, even if analysis fails
            self.banner.display_info("Cleaning up abricate_module...")
            self.cleanup_module_directory(abricate_module_path, fasta_files)
            self.banner.display_success("✅ abricate_module cleaned up successfully")

    def run_lineage_analysis(self, output_dir: Path) -> bool:
        """Run lineage database generation"""
        try:
            lineage_module_path = self.base_dir / "modules" / "lineage_module"
            lineage_script = lineage_module_path / "html_reference.py"
            
            if not lineage_script.exists():
                self.banner.display_error(f"Lineage script not found at: {lineage_script}")
                return False
            
            # Build command
            cmd = [sys.executable, str(lineage_script)]
            
            self.banner.display_info("Generating lineage reference database...")
            result = subprocess.run(cmd, capture_output=True, text=True, cwd=lineage_module_path)
            
            if result.returncode == 0:
                self.banner.display_success("Lineage reference database generated!")
                
                # Create lineage output directory
                lineage_output = output_dir / "lineage_results"
                lineage_output.mkdir(parents=True, exist_ok=True)
                
                # Copy lineage HTML to output directory
                lineage_html = lineage_module_path / "staphscope_lineage_reference.html"
                if lineage_html.exists():
                    target_html = lineage_output / "staphscope_lineage_reference.html"
                    shutil.copy2(lineage_html, target_html)
                    self.banner.display_success(f"Lineage reference copied to: {target_html}")
                return True
            else:
                self.banner.display_warning("Lineage database generation had warnings")
                if result.stderr:
                    error_lines = result.stderr.strip().split('\n')
                    for line in error_lines[:3]:
                        if line.strip():
                            print(f"  ⚠️  {line.strip()}")
                return True
                
        except Exception as e:
            self.banner.display_error(f"Lineage database generation failed: {str(e)}")
            return False

    # =========================================================================================
    # EXTRACTED FROM WORKING ORCHESTRATOR - COMPREHENSIVE REPORT & ULTIMATE REPORTER
    # =========================================================================================

    def copy_required_files_for_comprehensive_report(self, output_dir: Path) -> bool:
        """Copy required files to summary_module for comprehensive report - EXACT FROM WORKING VERSION"""
        try:
            self.banner.display_info("Copying required files for comprehensive report...")
            
            summary_module_path = self.base_dir / "modules" / "summary_module"
            
            # Check if comprehensive report script exists
            comprehensive_script = summary_module_path / "comprehensive_report.py"
            if not comprehensive_script.exists():
                self.banner.display_error(f"Comprehensive report script not found at: {comprehensive_script}")
                return False
            
            # Files to copy from each module's results directory - EXACT FROM WORKING VERSION
            files_to_copy = {
                "mlst_results": ["mlst_summary.tsv"],
                "spa_results": ["spa_summary.tsv"],
                "sccmec_results": ["staphscope_summary.tsv"]
            }
            
            copied_files = 0
            
            for source_dir, files in files_to_copy.items():
                source_path = output_dir / source_dir
                if not source_path.exists():
                    self.banner.display_warning(f"Source directory not found: {source_path}")
                    continue
                    
                for file in files:
                    source_file = source_path / file
                    if source_file.exists():
                        target_file = summary_module_path / file
                        shutil.copy2(source_file, target_file)
                        copied_files += 1
                        self.banner.display_info(f"  ✓ Copied: {file}")
                    else:
                        self.banner.display_warning(f"  ✗ File not found: {source_file}")
            
            if copied_files == 3:
                self.banner.display_success(f"All required files copied to summary_module")
                return True
            else:
                self.banner.display_warning(f"Only {copied_files}/3 required files copied")
                # Still try to run it - EXACT FROM WORKING VERSION
                return True
                
        except Exception as e:
            self.banner.display_error(f"Error copying files for comprehensive report: {str(e)}")
            return False

    def copy_html_reports_for_ultimate_reporter(self, output_dir: Path) -> bool:
        """Copy HTML reports needed by ultimate reporter to summary_module - EXACT FROM WORKING VERSION"""
        try:
            self.banner.display_info("Copying HTML reports for ultimate reporter...")
            
            summary_module_path = self.base_dir / "modules" / "summary_module"
            
            # List of HTML reports needed by ultimate reporter - EXACT FROM WORKING VERSION
            html_reports_needed = [
                # AMR reports
                ("amr_results", "staph_amrfinder_summary_report.html"),
                # ABRicate reports
                ("abricate_results", "staph_card_summary_report.html"),
                ("abricate_results", "staph_plasmidfinder_summary_report.html"),
                ("abricate_results", "staph_ncbi_summary_report.html"),
                ("abricate_results", "staph_vfdb_summary_report.html"),
                ("abricate_results", "staph_megares_summary_report.html"),
                ("abricate_results", "staph_resfinder_summary_report.html"),
                ("abricate_results", "staph_argannot_summary_report.html"),
                ("abricate_results", "staph_bacmet2_summary_report.html"),
            ]
            
            copied_count = 0
            
            for source_dir, html_file in html_reports_needed:
                source_path = output_dir / source_dir / html_file
                if source_path.exists():
                    target_path = summary_module_path / html_file
                    shutil.copy2(source_path, target_path)
                    copied_count += 1
                    self.banner.display_info(f"  ✓ Copied: {html_file}")
                else:
                    self.banner.display_warning(f"  ✗ Not found: {html_file}")
            
            # Also ensure comprehensive report is already there (it should be from previous step)
            comprehensive_file = summary_module_path / "staphscope_comprehensive_report.html"
            if comprehensive_file.exists():
                self.banner.display_info(f"  ✓ Comprehensive report already exists")
            else:
                self.banner.display_error("Comprehensive report not found - required for ultimate reporter")
                return False
            
            self.banner.display_success(f"Copied {copied_count} HTML reports to summary_module")
            return True
            
        except Exception as e:
            self.banner.display_error(f"Error copying HTML reports: {str(e)}")
            return False

    def run_comprehensive_report(self, output_dir: Path) -> bool:
        """Run comprehensive report generation - EXACT FROM WORKING VERSION"""
        try:
            summary_module_path = self.base_dir / "modules" / "summary_module"
            comprehensive_script = summary_module_path / "comprehensive_report.py"
            
            if not comprehensive_script.exists():
                self.banner.display_error(f"Comprehensive report script not found: {comprehensive_script}")
                return False
            
            # Copy required files to summary_module - EXACT FROM WORKING VERSION
            if not self.copy_required_files_for_comprehensive_report(output_dir):
                self.banner.display_warning("Skipping comprehensive report due to missing files")
                return False
            
            # Run comprehensive report - EXACT FROM WORKING VERSION
            self.banner.display_info("Running comprehensive report...")
            
            cmd = [sys.executable, str(comprehensive_script)]
            result = subprocess.run(cmd, capture_output=True, text=True, cwd=summary_module_path)
            
            if result.returncode == 0:
                self.banner.display_success("Comprehensive report generated successfully!")
                return True
            else:
                self.banner.display_warning("Comprehensive report generation had issues")
                if result.stderr:
                    error_lines = result.stderr.strip().split('\n')
                    for line in error_lines:
                        if "error" in line.lower() or "failed" in line.lower():
                            print(f"  ⚠️  Comprehensive report error: {line}")
                # Still return True to continue pipeline - EXACT FROM WORKING VERSION
                return True
                
        except Exception as e:
            self.banner.display_error(f"Comprehensive report generation failed: {str(e)}")
            return False

    def run_ultimate_reporter(self, output_dir: Path) -> bool:
        """Run ultimate reporter - EXACT FROM WORKING VERSION"""
        try:
            summary_module_path = self.base_dir / "modules" / "summary_module"
            ultimate_script = summary_module_path / "staphscope_ultimate_reporter.py"
            
            if not ultimate_script.exists():
                self.banner.display_error(f"Ultimate reporter script not found: {ultimate_script}")
                return False
            
            # Copy HTML reports needed by ultimate reporter
            if not self.copy_html_reports_for_ultimate_reporter(output_dir):
                self.banner.display_warning("Skipping ultimate reporter due to missing HTML reports")
                return False
            
            # Run ultimate reporter silently - EXACT FROM WORKING VERSION
            self.banner.display_info("Running ultimate reporter (silent mode)...")
            
            cmd = [sys.executable, str(ultimate_script), "-i", "."]
            result = subprocess.run(cmd, capture_output=True, text=True, cwd=summary_module_path)
            
            if result.returncode == 0:
                self.banner.display_success("Ultimate reporter completed successfully!")
                return True
            else:
                self.banner.display_warning("Ultimate reporter had issues")
                if result.stderr:
                    error_lines = result.stderr.strip().split('\n')
                    for line in error_lines:
                        if "error" in line.lower() or "failed" in line.lower():
                            print(f"  ⚠️  Ultimate reporter error: {line}")
                return True
                
        except Exception as e:
            self.banner.display_error(f"Ultimate reporter failed: {str(e)}")
            return False

    def copy_summary_results_to_final_directory(self, output_dir: Path):
        """Copy only comprehensive report files and ultimate reporter directory to final location - EXACT FROM WORKING VERSION"""
        try:
            self.banner.display_info("Copying summary results to final directory...")
            
            summary_module_path = self.base_dir / "modules" / "summary_module"
            final_report_dir = output_dir / "Staphscope_final_report"
            final_report_dir.mkdir(parents=True, exist_ok=True)
            
            # Files to copy from comprehensive reporter - EXACT FROM WORKING VERSION
            comprehensive_files = [
                "staphscope_comprehensive_report.html",
                "staphscope_comprehensive_report.json",
                "staphscope_comprehensive_report.tsv"
            ]
            
            # Copy comprehensive report files
            for file_name in comprehensive_files:
                source_file = summary_module_path / file_name
                if source_file.exists():
                    target_file = final_report_dir / file_name
                    shutil.copy2(source_file, target_file)
                    self.banner.display_info(f"  ✓ Copied: {file_name}")
                else:
                    self.banner.display_warning(f"  ✗ Not found: {file_name}")
            
            # Copy ultimate reporter directory - EXACT FROM WORKING VERSION
            ultimate_reports_dir = summary_module_path / "STAPHSCOPE_ULTIMATE_REPORTS"
            if ultimate_reports_dir.exists() and ultimate_reports_dir.is_dir():
                target_ultimate_dir = final_report_dir / "STAPHSCOPE_ULTIMATE_REPORTS"
                if target_ultimate_dir.exists():
                    shutil.rmtree(target_ultimate_dir)
                shutil.copytree(ultimate_reports_dir, target_ultimate_dir)
                
                # Count files in ultimate directory
                ultimate_files = list(ultimate_reports_dir.glob("*"))
                self.banner.display_info(f"  ✓ Copied: STAPHSCOPE_ULTIMATE_REPORTS directory ({len(ultimate_files)} files)")
            
            self.banner.display_success(f"✅ All summary results copied to: {final_report_dir}")
            
            # Display what was copied - EXACT FROM WORKING VERSION
            self.banner.display_info("Summary Reports Generated:")
            all_files = list(final_report_dir.glob("*"))
            for file_path in sorted(all_files):
                if file_path.is_dir():
                    dir_files = list(file_path.glob("*"))
                    self.banner.display_info(f"  📁 {file_path.name} ({len(dir_files)} files)")
                else:
                    self.banner.display_info(f"  📄 {file_path.name}")
            
        except Exception as e:
            self.banner.display_error(f"Error copying summary results: {str(e)}")

    # =========================================================================================================================== #
    # VISUALIZATION METHODS (KEEPING MY EXISTING LOGIC)-BROWN BECKLEY                                                             #
    # =========================================================================================================================== #

    def copy_html_files_to_visualization_module(self, output_dir: Path) -> bool:
        """Copy HTML files from other modules to visualization module"""
        try:
            self.banner.display_info("Copying HTML files to visualization module...")
            
            vis_module_path = self.base_dir / "modules" / "visualization_module"
            
            # Check if visualization script exists
            vis_script = vis_module_path / "staphscope_visualizer.py"
            if not vis_script.exists():
                self.banner.display_error(f"Visualization script not found at: {vis_script}")
                return False
            
            # Files to copy from various output directories
            files_to_copy = [
                # From amr_results
                (output_dir / "amr_results", "staph_amrfinder_summary_report.html"),
                # From abricate_results
                (output_dir / "abricate_results", "staph_card_summary_report.html"),
                (output_dir / "abricate_results", "staph_plasmidfinder_summary_report.html"),
                (output_dir / "abricate_results", "staph_ncbi_summary_report.html"),
                (output_dir / "abricate_results", "staph_vfdb_summary_report.html"),
                (output_dir / "abricate_results", "staph_megares_summary_report.html"),
                (output_dir / "abricate_results", "staph_resfinder_summary_report.html"),
                (output_dir / "abricate_results", "staph_argannot_summary_report.html"),
                # From final report
                (output_dir / "Staphscope_final_report", "staphscope_comprehensive_report.html"),
            ]
            
            copied_count = 0
            missing_count = 0
            
            for source_dir, html_file in files_to_copy:
                source_path = source_dir / html_file
                target_path = vis_module_path / html_file
                
                if source_path.exists():
                    shutil.copy2(source_path, target_path)
                    copied_count += 1
                    self.banner.display_info(f"  ✓ Copied: {html_file}")
                else:
                    missing_count += 1
                    self.banner.display_warning(f"  ⚠️  Not found: {html_file}")
            
            self.banner.display_info(f"Copied {copied_count} files, {missing_count} files missing")
            
            # Check if we have at least the comprehensive report (most important)
            comprehensive_path = vis_module_path / "staphscope_comprehensive_report.html"
            if not comprehensive_path.exists():
                self.banner.display_error("Comprehensive report not found - required for visualization")
                return False
            
            return True
            
        except Exception as e:
            self.banner.display_error(f"Error copying HTML files to visualization module: {str(e)}")
            return False

    def run_visualization_analysis(self, output_dir: Path) -> bool:
        """Run visualization analysis - last step in pipeline"""
        vis_module_path = self.base_dir / "modules" / "visualization_module"
        
        try:
            vis_script = vis_module_path / "staphscope_visualizer.py"
            
            if not vis_script.exists():
                self.banner.display_error(f"Visualization script not found at: {vis_script}")
                return False
            
            # Copy HTML files from other modules to visualization module
            if not self.copy_html_files_to_visualization_module(output_dir):
                self.banner.display_warning("Some HTML files missing for visualization, but proceeding anyway...")
            
            # Build command - just python staphscope_visualizer.py
            cmd = [sys.executable, str(vis_script)]
            
            self.banner.display_info("Running visualization module...")
            result = subprocess.run(cmd, capture_output=True, text=True, cwd=vis_module_path)
            
            if result.returncode == 0:
                self.banner.display_success("Visualization completed successfully!")
                
                # Check if visualizations were generated
                vis_output_dir = vis_module_path / "STAPHSCOPE_VISUALIZATIONS"
                if vis_output_dir.exists() and vis_output_dir.is_dir():
                    # Copy STAPHSCOPE_VISUALIZATIONS to output directory
                    target_vis_dir = output_dir / "STAPHSCOPE_VISUALIZATIONS"
                    if target_vis_dir.exists():
                        shutil.rmtree(target_vis_dir)
                    shutil.copytree(vis_output_dir, target_vis_dir)
                    
                    # Count files in visualization directory
                    vis_files = list(vis_output_dir.rglob("*"))
                    html_files = [f for f in vis_files if f.suffix == '.html']
                    image_files = [f for f in vis_files if f.suffix in ['.png', '.jpg', '.jpeg', '.svg']]
                    
                    self.banner.display_success(f"✅ Visualizations copied to: {target_vis_dir}")
                    self.banner.display_info(f"   📊 {len(html_files)} HTML reports")
                    self.banner.display_info(f"   🖼️  {len(image_files)} visualization images")
                    
                    # List the generated files
                    print("\n📁 Generated visualization files:")
                    for file_path in sorted(vis_output_dir.glob("*")):
                        if file_path.is_file():
                            print(f"   📄 {file_path.name}")
                        elif file_path.is_dir():
                            sub_files = list(file_path.glob("*"))
                            print(f"   📂 {file_path.name}/ ({len(sub_files)} files)")
                else:
                    self.banner.display_warning("Visualization output directory not found")
                
                return True
            else:
                self.banner.display_warning("Visualization had issues")
                if result.stderr:
                    error_lines = result.stderr.strip().split('\n')
                    for line in error_lines[:5]:
                        if line.strip():
                            print(f"  ⚠️  {line.strip()}")
                return True
                
        except Exception as e:
            self.banner.display_error(f"Visualization failed: {str(e)}")
            return False

    # ========================================================================
    # SEQUENTIAL ANALYSES METHOD
    # ========================================================================

    def run_sequential_analyses(self, fasta_files: List[Path], output_dir: Path, threads: int, 
                               skip_modules: Dict[str, bool]) -> Dict[str, bool]:
        """Run analyses SEQUENTIALLY to keep module messages together"""
        # Define analysis functions in order with their display headers
        analysis_functions = [
            ("FASTA QC", self.run_fasta_qc_analysis, "FASTA QC Analysis", "Sequence Quality Control & Statistics", not skip_modules.get('fasta_qc', False)),
            ("MLST", self.run_mlst_analysis, "MLST Analysis", "Multi-Locus Sequence Typing", not skip_modules.get('mlst', False)),
            ("spa typing", self.run_spa_typing, "SPA TYPING ANALYSIS", "Staphylococcal Protein A Typing", not skip_modules.get('spa', False)),
            ("SCCmec", self.run_sccmec_analysis, "SCCMEC ANALYSIS", "Methicillin Resistance Cassette Typing", not skip_modules.get('sccmec', False)),
            ("AMRFinderPlus", self.run_amrfinder_analysis, "AMR ANALYSIS", "Antimicrobial Resistance Gene Detection", not skip_modules.get('amr', False)),
            ("ABRicate", self.run_abricate_analysis, "ABRICATE ANALYSIS", "Comprehensive Resistance, Plasmid & Virulence Gene Screening", not skip_modules.get('abricate', False))
        ]
        
        # Filter out skipped analyses
        active_analyses = [(name, func, header, description) for name, func, header, description, enabled in analysis_functions if enabled]
        
        if not active_analyses:
            self.banner.display_warning("All analyses were skipped! Nothing to run.")
            return {}
        
        self.banner.display_info(f"Running {len(active_analyses)} analyses")
        
        results = {}
        
        # Run analyses SEQUENTIALLY
        for analysis_name, analysis_func, header, description in active_analyses:
            # Display module header before each analysis
            self.banner.display_module_header(header, description)
            
            try:
                success = analysis_func(fasta_files, output_dir, max(1, threads // len(active_analyses)))
                results[analysis_name] = success
                
                if success:
                    self.banner.display_success(f"✅ {analysis_name} completed")
                else:
                    self.banner.display_error(f"❌ {analysis_name} failed")
                    
            except Exception as e:
                self.banner.display_error(f"❌ {analysis_name} failed with exception: {str(e)}")
                results[analysis_name] = False
            
            # Add spacing between modules
            print()
        
        return results

    # ========================================================================
    # MAIN COMPLETE ANALYSIS METHOD
    # ========================================================================

    def run_complete_analysis(self, input_path: str, output_dir: str, threads: int = 1, 
                             skip_modules: Dict[str, bool] = None, skip_comprehensive: bool = False,
                             skip_visualization: bool = False):
        """Run complete StaphScope analysis pipeline with FASTA QC first and Visualization last"""
        if skip_modules is None:
            skip_modules = {}
        
        start_time = datetime.now()
        
        try:
            # Display beautiful startup sequence
            self.banner.display_startup_sequence()
            self.banner.display_banner(show_quote=True, show_author=True)
            
            # Create output directory
            output_path = Path(output_dir)
            output_path.mkdir(parents=True, exist_ok=True)
            
            # Find input files
            fasta_files = self.find_fasta_files(input_path)
            
            if not fasta_files:
                self.banner.display_error("No FASTA files found! Analysis stopped.")
                return
            
            # Show file formats detected
            extensions = set(f.suffix.lower() for f in fasta_files)
            self.banner.display_success(f"Starting analysis of {len(fasta_files)} samples")
            self.banner.display_info(f"File formats detected: {', '.join(extensions)}")
            
            # Create output structure with new FASTA QC directory
            subdirs = [
                "fasta_qc_results", "mlst_results", "spa_results", "sccmec_results",
                "abricate_results", "amr_results", "lineage_results",
                "Staphscope_final_report"
            ]
            for subdir in subdirs:
                (output_path / subdir).mkdir(exist_ok=True)
            
            # Display analysis plan
            self.banner.display_module_header("Analysis Plan", "Modules to be executed")
            analyses_to_run = [
                ("FASTA QC", not skip_modules.get('fasta_qc', False)),
                ("MLST", not skip_modules.get('mlst', False)),
                ("spa typing", not skip_modules.get('spa', False)),
                ("SCCmec", not skip_modules.get('sccmec', False)),
                ("AMRFinderPlus", not skip_modules.get('amr', False)),
                ("ABRicate", not skip_modules.get('abricate', False)),
                ("Lineage Reference", not skip_modules.get('lineage', False)),
                ("Comprehensive Report", not skip_comprehensive),
                ("Ultimate Reporter", not skip_comprehensive),
                ("Visualization", not skip_visualization)
            ]
            
            for analysis, enabled in analyses_to_run:
                status = "✅ ENABLED" if enabled else "⏸️  SKIPPED"
                print(f"   {status} - {analysis}")
            
            # Force output flush before starting analyses
            sys.stdout.flush()
            
            # Run main analyses SEQUENTIALLY (FASTA QC first, then others)
            analysis_results = self.run_sequential_analyses(fasta_files, output_path, threads, skip_modules)
            
            # Run lineage analysis if not skipped
            if not skip_modules.get('lineage', False):
                self.banner.display_module_header("Lineage Database", "S. aureus Lineage Reference Generation")
                lineage_success = self.run_lineage_analysis(output_path)
                analysis_results["Lineage Reference"] = lineage_success
                print()
            
            # Run comprehensive report if not skipped
            if not skip_comprehensive:
                self.banner.display_module_header("Comprehensive Report", "Unified MLST, spa & SCCmec Analysis")
                comprehensive_success = self.run_comprehensive_report(output_path)
                analysis_results["Comprehensive Report"] = comprehensive_success
                
                # Run ultimate reporter immediately after comprehensive report - EXACT FROM WORKING VERSION
                if comprehensive_success:
                    self.banner.display_module_header("Ultimate Reporter", "Gene-centric Integrated Analysis")
                    ultimate_success = self.run_ultimate_reporter(output_path)
                    analysis_results["Ultimate Reporter"] = ultimate_success
                    
                    # Copy results to final directory
                    self.copy_summary_results_to_final_directory(output_path)
                else:
                    self.banner.display_warning("Skipping ultimate reporter due to comprehensive report failure")
                
                print()
            
            # Run visualization if not skipped (LAST step)
            if not skip_visualization:
                self.banner.display_module_header("Visualization", "Interactive Visualizations & Dashboard")
                visualization_success = self.run_visualization_analysis(output_path)
                analysis_results["Visualization"] = visualization_success
                print()
            
            # Calculate analysis time
            analysis_time = datetime.now() - start_time
            analysis_time_str = str(analysis_time).split('.')[0]
            
            # Display beautiful completion footer
            successful_count = sum(analysis_results.values())
            total_count = len(analysis_results)
            
            self.banner.display_footer(
                analysis_time=analysis_time_str,
                samples_processed=len(fasta_files)
            )
            
            # Final status
            if successful_count == total_count:
                self.banner.display_success(f"🎉 All {total_count} analyses completed successfully!")
                self.banner.display_success("🧹 All module directories have been cleaned up")
                
                # Show final output structure
                print("\n📁 FINAL OUTPUT STRUCTURE:")
                for subdir in sorted(output_path.glob("*")):
                    if subdir.is_dir():
                        files_count = len(list(subdir.rglob("*")))
                        if subdir.name == "Staphscope_final_report":
                            # Show detailed contents of final report
                            print(f"   📂 {subdir.name}/")
                            for item in sorted(subdir.glob("*")):
                                if item.is_file():
                                    print(f"      📄 {item.name}")
                                elif item.is_dir():
                                    sub_items = len(list(item.glob("*")))
                                    print(f"      📂 {item.name}/ ({sub_items} items)")
                        else:
                            print(f"   📂 {subdir.name}/ ({files_count} items)")
            
            else:
                self.banner.display_warning(
                    f"⚠️  {successful_count}/{total_count} analyses completed successfully."
                )

    # ===== NEW CITATION BLOCK =====
            print("\n📚 Please cite our StaphScope paper:")
            print("   Beckley, B., Amarh, V. StaphScope: a species-optimized computational pipeline for rapid and accessible Staphylococcus aureus genotyping and surveillance. BMC Genomics (2026).")
            print("   https://doi.org/10.1186/s12864-026-12609-x")
            
    # ===== END CITATION BLOCK =====
                
        except KeyboardInterrupt:
            self.banner.display_error("Analysis interrupted by user")
        except Exception as e:
            self.banner.display_error(f"Critical error in analysis pipeline: {str(e)}")
            import traceback
            traceback.print_exc()

def main():
    """Main entry point for StaphScope"""
    parser = argparse.ArgumentParser(
        description="StaphScope: Advanced Staphylococcus aureus Typing & Lineage Analysis Platform with FASTA QC & Visualization",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  staphscope -i genome.fna -o results/
  staphscope -i "*.fna" -o batch_results --threads 8
  staphscope -i "*.fasta" -o analysis --threads 16 --skip-lineage
  staphscope -i "genome*.fa" -o results/ --threads 4 --skip-comprehensive
  staphscope -i "samples/*.fasta" -o output/ --skip-visualization

Supported FASTA formats: .fna, .fasta, .fa, .fn, .faa

Analysis Modules:
  • FASTA QC (Quality Control & Statistics) 
  • MLST (Multi-Locus Sequence Typing)
  • spa typing (Staphylococcal Protein A)  
  • SCCmec typing (Methicillin Resistance Cassette)
  • AMR profiling (Antimicrobial Resistance)
  • ABRicate (Comprehensive resistance/Plasmid/virulence)
  • Lineage reference database
  • Comprehensive report (MLST + spa + SCCmec summary)
  • Ultimate reporter (Gene-centric integrated analysis) - runs with comprehensive report
  • Visualization (Interactive dashboards & plots) 

Output: Comprehensive results for all analyses in organized directories

Please run abricate --setupdb for recent gene annotations!!!

⭐ Star us on GitHub if you find this tool useful!

Transforming fragmented genomic data into coherent biological narratives 🧬✨
        """
    )
    
    parser.add_argument('-i', '--input', required=True,
                       help='Input FASTA file(s) - can use glob patterns like "*.fna" or "*.fasta"')
    parser.add_argument('-o', '--output', required=True,
                       help='Output directory for all results')
    parser.add_argument('-t', '--threads', type=int, default=2,
                       help='Number of threads (default: 2)')
    
    # Skip options
    parser.add_argument('--skip-fasta-qc', action='store_true', 
                       help='Skip FASTA QC analysis')
    parser.add_argument('--skip-amr', action='store_true', 
                       help='Skip AMR analysis (AMRfinderPlus)')
    parser.add_argument('--skip-abricate', action='store_true',
                       help='Skip ABRicate analysis')
    parser.add_argument('--skip-mlst', action='store_true',
                       help='Skip MLST analysis')
    parser.add_argument('--skip-spa', action='store_true',
                       help='Skip spa typing analysis')
    parser.add_argument('--skip-sccmec', action='store_true',
                       help='Skip SCCmec analysis')
    parser.add_argument('--skip-lineage', action='store_true',
                       help='Skip lineage reference generation')
    parser.add_argument('--skip-comprehensive', action='store_true',
                       help='Skip comprehensive report generation (MLST + spa + SCCmec) AND ultimate reporter')
    parser.add_argument('--skip-visualization', action='store_true',
                       help='Skip visualization module (dashboards & plots)')
    
    args = parser.parse_args()
    
    # Create skip modules dictionary
    skip_modules = {
        'fasta_qc': args.skip_fasta_qc,
        'amr': args.skip_amr,
        'abricate': args.skip_abricate,
        'mlst': args.skip_mlst,
        'spa': args.skip_spa,
        'sccmec': args.skip_sccmec,
        'lineage': args.skip_lineage
    }
    
    # Create and run StaphScope
    staphscope = StaphScopeOrchestrator()
    
    try:
        staphscope.run_complete_analysis(
            input_path=args.input,
            output_dir=args.output,
            threads=args.threads,
            skip_modules=skip_modules,
            skip_comprehensive=args.skip_comprehensive,
            skip_visualization=args.skip_visualization
        )
    except KeyboardInterrupt:
        print("\n❌ Analysis interrupted by user")
    except Exception as e:
        print(f"\n💥 Critical error: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)

if __name__ == "__main__":
    main()
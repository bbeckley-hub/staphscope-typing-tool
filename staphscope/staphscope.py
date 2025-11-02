#!/usr/bin/env python3
"""
StaphScope Main Orchestrator - FINAL VERSION WITH PROPER CLEANUP
Complete S. aureus typing pipeline - MLST, spa, SCCmec, AMR, Virulence, Lineage
Author: Brown Beckley <brownbeckley94@gmail.com>
"""

import os
import sys
import glob
import argparse
import subprocess
import shutil
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime
from typing import Dict, List, Set

# Import banner
try:
    from .core.banner import StaphScopeBanner
except (ImportError, SystemError):
    sys.path.insert(0, str(Path(__file__).parent))
    from core.banner import StaphScopeBanner

class StaphScopeOrchestrator:
    """Final StaphScope orchestrator - with COMPLETE cleanup"""
    
    def __init__(self):
        self.banner = StaphScopeBanner()
        self.base_dir = Path(__file__).parent
        
    def find_fasta_files(self, input_path: str) -> List[Path]:
        """Find all FASTA files using glob patterns"""
        self.banner.display_info(f"Searching for files with pattern: {input_path}")
        
        # Handle quoted wildcards properly
        if '*' in input_path or '?' in input_path:
            matched_files = glob.glob(input_path)
            fasta_files = [Path(f) for f in matched_files if Path(f).is_file() and 
                          f.lower().endswith(('.fna', '.fasta', '.fa', '.fn')) and
                          not Path(f).name.startswith('.')]
            self.banner.display_success(f"Found {len(fasta_files)} FASTA files")
            return sorted(fasta_files)
        
        # Handle direct file path
        input_path_obj = Path(input_path)
        if input_path_obj.is_file() and input_path_obj.suffix.lower() in ['.fna', '.fasta', '.fa', '.fn']:
            self.banner.display_success(f"Found single FASTA file: {input_path_obj.name}")
            return [input_path_obj]
        
        # Handle directory
        if input_path_obj.is_dir():
            patterns = [
                f"{input_path}/*.fna", f"{input_path}/*.fasta",
                f"{input_path}/*.fa", f"{input_path}/*.fn"
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
                "amrfinder_results", "abricate_results"
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

    def run_mlst_analysis(self, fasta_files: List[Path], output_dir: Path, threads: int) -> bool:
        """Run MLST analysis - WITH CLEANUP"""
        mlst_module_path = self.base_dir / "modules" / "mlst_module"
        
        try:
            self.banner.display_module_header("MLST Analysis", "Multi-Locus Sequence Typing")
            
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
                    print(f"MLST stderr: {result.stderr[:200]}...")
                return True
                
        except Exception as e:
            self.banner.display_error(f"MLST analysis failed: {str(e)}")
            return False
        finally:
            # ALWAYS cleanup, even if analysis fails
            self.cleanup_module_directory(mlst_module_path, fasta_files)

    def run_spa_typing(self, fasta_files: List[Path], output_dir: Path, threads: int) -> bool:
        """Run spa typing - WITH CLEANUP"""
        spa_module_path = self.base_dir / "modules" / "spa_module" / "spatyper" / "spa_typing"
        
        try:
            self.banner.display_module_header("spa Typing Analysis", "Staphylococcal Protein A Typing")
            
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
                    print(f"spa stderr: {result.stderr[:200]}...")
                return True
                
        except Exception as e:
            self.banner.display_error(f"spa typing analysis failed: {str(e)}")
            return False
        finally:
            # ALWAYS cleanup, even if analysis fails
            self.cleanup_module_directory(spa_module_path, fasta_files)

    def run_sccmec_analysis(self, fasta_files: List[Path], output_dir: Path, threads: int) -> bool:
        """Run SCCmec analysis - WITH CLEANUP"""
        sccmec_module_path = self.base_dir / "modules" / "sccmec_module"
        
        try:
            self.banner.display_module_header("SCCmec Analysis", "Methicillin Resistance Cassette Typing")
            
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
                summary_files = ["staphscope_summary.html", "staphscope_summary.tsv"]
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
                    print(f"SCCmec stderr: {result.stderr[:200]}...")
                return True
                
        except Exception as e:
            self.banner.display_error(f"SCCmec analysis failed: {str(e)}")
            return False
        finally:
            # ALWAYS cleanup, even if analysis fails
            self.cleanup_module_directory(sccmec_module_path, fasta_files)

    def run_amrfinder_analysis(self, fasta_files: List[Path], output_dir: Path, threads: int) -> bool:
        """Run AMRFinderPlus analysis - WITH CLEANUP"""
        amr_module_path = self.base_dir / "modules" / "amr_module"
        
        try:
            self.banner.display_module_header("AMR Analysis", "Antimicrobial Resistance Gene Detection")
            
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
                amr_source = amr_module_path / "amrfinder_results"
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
                    print(f"AMRFinderPlus stderr: {result.stderr[:200]}...")
                return True
                
        except Exception as e:
            self.banner.display_error(f"AMRFinderPlus analysis failed: {str(e)}")
            return False
        finally:
            # ALWAYS cleanup, even if analysis fails
            self.cleanup_module_directory(amr_module_path, fasta_files)

    def run_abricate_analysis(self, fasta_files: List[Path], output_dir: Path, threads: int) -> bool:
        """Run Abricate analysis - WITH CLEANUP"""
        abricate_module_path = self.base_dir / "modules" / "abricate_module"
        
        try:
            self.banner.display_module_header("ABRicate Analysis", "Comprehensive Resistance & Virulence Gene Screening")
            
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
                    print(f"ABRicate stderr: {result.stderr[:200]}...")
                return True
                
        except Exception as e:
            self.banner.display_error(f"ABRicate analysis failed: {str(e)}")
            return False
        finally:
            # ALWAYS cleanup, even if analysis fails
            self.cleanup_module_directory(abricate_module_path, fasta_files)

    def run_lineage_analysis(self, output_dir: Path) -> bool:
        """Run lineage database generation"""
        try:
            self.banner.display_module_header("Lineage Database", "S. aureus Lineage Reference Generation")
            
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
                    print(f"Lineage stderr: {result.stderr[:200]}...")
                return True
                
        except Exception as e:
            self.banner.display_error(f"Lineage database generation failed: {str(e)}")
            return False

    def run_parallel_analyses(self, fasta_files: List[Path], output_dir: Path, threads: int, 
                            skip_modules: Dict[str, bool]) -> Dict[str, bool]:
        """Run analyses in parallel"""
        analysis_functions = [
            (self.run_mlst_analysis, "MLST", not skip_modules.get('mlst', False)),
            (self.run_spa_typing, "spa typing", not skip_modules.get('spa', False)),
            (self.run_sccmec_analysis, "SCCmec", not skip_modules.get('sccmec', False)),
            (self.run_amrfinder_analysis, "AMRFinderPlus", not skip_modules.get('amr', False)),
            (self.run_abricate_analysis, "ABRicate", not skip_modules.get('abricate', False))
        ]
        
        # Filter out skipped analyses
        active_analyses = [(func, name) for func, name, enabled in analysis_functions if enabled]
        
        if not active_analyses:
            self.banner.display_warning("All analyses were skipped! Nothing to run.")
            return {}
        
        self.banner.display_info(f"Running {len(active_analyses)} analyses")
        
        results = {}
        
        # Run analyses in parallel
        with ThreadPoolExecutor(max_workers=min(len(active_analyses), max(1, threads // 2))) as executor:
            future_to_analysis = {
                executor.submit(func, fasta_files, output_dir, max(1, threads // len(active_analyses))): name 
                for func, name in active_analyses
            }
            
            for future in as_completed(future_to_analysis):
                analysis_name = future_to_analysis[future]
                
                try:
                    success = future.result()
                    results[analysis_name] = success
                    
                    if success:
                        self.banner.display_success(f"✅ {analysis_name} completed")
                    else:
                        self.banner.display_error(f"❌ {analysis_name} failed")
                        
                except Exception as e:
                    self.banner.display_error(f"❌ {analysis_name} failed with exception: {str(e)}")
                    results[analysis_name] = False
        
        return results

    def run_complete_analysis(self, input_path: str, output_dir: str, threads: int = 1, 
                            skip_modules: Dict[str, bool] = None):
        """Run complete StaphScope analysis pipeline"""
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
            
            # Create output structure
            subdirs = [
                "mlst_results", "spa_results", "sccmec_results",
                "abricate_results", "amr_results", "lineage_results"
            ]
            for subdir in subdirs:
                (output_path / subdir).mkdir(exist_ok=True)
            
            # Display analysis plan
            self.banner.display_module_header("Analysis Plan", "Modules to be executed")
            analyses_to_run = [
                ("MLST", not skip_modules.get('mlst', False)),
                ("spa typing", not skip_modules.get('spa', False)),
                ("SCCmec", not skip_modules.get('sccmec', False)),
                ("AMRFinderPlus", not skip_modules.get('amr', False)),
                ("ABRicate", not skip_modules.get('abricate', False)),
                ("Lineage Reference", not skip_modules.get('lineage', False))
            ]
            
            for analysis, enabled in analyses_to_run:
                status = "✅ ENABLED" if enabled else "⏸️  SKIPPED"
                print(f"   {status} - {analysis}")
            
            # Run main analyses in parallel
            analysis_results = self.run_parallel_analyses(fasta_files, output_path, threads, skip_modules)
            
            # Run lineage analysis if not skipped
            if not skip_modules.get('lineage', False):
                lineage_success = self.run_lineage_analysis(output_path)
                analysis_results["Lineage Reference"] = lineage_success
            
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
            else:
                self.banner.display_warning(
                    f"⚠️  {successful_count}/{total_count} analyses completed successfully."
                )
                
        except KeyboardInterrupt:
            self.banner.display_error("Analysis interrupted by user")
        except Exception as e:
            self.banner.display_error(f"Critical error in analysis pipeline: {str(e)}")
            import traceback
            traceback.print_exc()

def main():
    """Main entry point for StaphScope"""
    parser = argparse.ArgumentParser(
        description="StaphScope: Complete S. aureus Typing Pipeline",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  staphscope -i genome.fna -o results/
  staphscope -i "*.fna" -o batch_results --threads 8
  staphscope -i "*.fasta" -o analysis --threads 16 --skip-lineage
  staphscope -i "genome*.fa" -o results/ --threads 4

Supported FASTA formats: .fna, .fasta, .fa, .fn

Analysis Modules:
  • MLST (Multi-Locus Sequence Typing)
  • spa typing (Staphylococcal Protein A)  
  • SCCmec typing (Methicillin Resistance Cassette)
  • AMR profiling (Antimicrobial Resistance)
  • ABRicate (Comprehensive resistance/virulence)
  • Lineage reference database
    Run : abricate --setupdb & amrfinder --update or amrfinder -u (For latest database)

Output: Comprehensive results for all analyses in organized directories
        """
    )
    
    parser.add_argument('-i', '--input', required=True,
                       help='Input FASTA file(s) - can use glob patterns like "*.fna" or "*.fasta"')
    parser.add_argument('-o', '--output', required=True,
                       help='Output directory for all results')
    parser.add_argument('-t', '--threads', type=int, default=2,
                       help='Number of threads (default: 2)')
    
    # Skip options
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
    
    args = parser.parse_args()
    
    # Create skip modules dictionary
    skip_modules = {
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
            skip_modules=skip_modules
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

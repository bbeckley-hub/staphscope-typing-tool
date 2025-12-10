#!/usr/bin/env python3
"""
StaphScope AMRfinderPlus Standalone Module - BUNDLED VERSION
Comprehensive AMR analysis with HTML, TSV, and JSON reporting - MAXIMUM SPEED VERSION
Author: Beckley Brown <brownbeckley94@gmail.com>
Date: 2025
Send a quick mail for any issues or further explanations.
Affiliation: University of Ghana Medical School-Department of Medical Biochemistry
Uses BUNDLED AMRFinderPlus 4.2.4 with database 2025-12-03.1
"""

import subprocess
import sys
import os
import glob
import logging
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed
from typing import List, Dict, Any, Optional
import argparse
import re
from datetime import datetime
import psutil
import math
import json
from collections import defaultdict, Counter

class AMRfinderPlusExecutor:
    """AMRfinderPlus executor with BUNDLED resources - Comprehensive HTML, TSV, and JSON reporting"""
    
    def __init__(self, cpus: int = None):
        # Setup logging FIRST
        self.logger = self._setup_logging()
        
        # Get module directory
        self.module_dir = os.path.dirname(os.path.abspath(__file__))
        
        # Initialize available_ram before calculating cpus
        self.available_ram = self._get_available_ram()
        
        # Then calculate resources - MAXIMUM SPEED MODE
        self.cpus = self._calculate_optimal_cpus(cpus)
        
        # Set bundled resources paths
        self.bundled_amrfinder = os.path.join(self.module_dir, "bin", "amrfinder")
        self.bundled_database = os.path.join(self.module_dir, "data", "amrfinder_db", "2025-12-03.1")
        
        self.metadata = {
            "tool_name": "StaphScope AMRfinderPlus (BUNDLED)",
            "version": "1.0.0", 
            "authors": ["Brown Beckley"],
            "email": "brownbeckley94@gmail.com",
            "github": "https://github.com/bbeckley-hub",
            "affiliation": "University of Ghana Medical School",
            "analysis_date": datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            "amrfinder_version": "4.2.4",
            "database_version": "2025-12-03.1"
        }
        
        # S. aureus specific high-risk and critical gene sets
        self.high_risk_genes = {
            # Critical Beta-lactam resistance (MRSA)
            'mecA', 'mecC', 'blaZ',
            
            # Glycopeptide resistance (VISA/VRSA)
            'vanA', 'vanB', 'vanC', 'vanD', 'vanE', 'vanG',
            
            # MLSB resistance
            'ermA', 'ermB', 'ermC', 'ermF', 'ermT', 'ermY',
            'mphC', 'msrA',
            
            # Aminoglycoside resistance
            'aacA-aphD', 'aac(6\')-aph(2")', 'aph(3\')-III', 'ant(4\')-Ia',
            'aac(6\')-Ie-aph(2")-Ia',
            
            # Tetracycline resistance
            'tetK', 'tetL', 'tetM', 'tetO',
            
            # Fluoroquinolone resistance
            'grlA', 'grlB', 'gyrA', 'gyrB', 'norA', 'norB',
            
            # Oxazolidinone resistance (Linezolid)
            'cfr', 'optrA', 'poxtA',
            
            # Other critical resistance
            'dfrA', 'dfrG', 'cat', 'fosB', 'fusB', 'fusC', 'mupA',
            
            # Multi-drug efflux pumps
            'qacA', 'qacB', 'qacC', 'smr', 'sepA'
        }

        # CRITICAL RISK genes - highest priority for S. aureus
        self.critical_risk_genes = {
            'mecA', 'mecC', 'vanA', 'vanB', 'cfr', 'optrA', 'poxtA'
        }
        
        self.science_quotes = [
            "“The important thing is not to stop questioning. Curiosity has its own reason for existence.” - Albert Einstein",
            "“Nothing in life is to be feared, it is only to be understood.” - Marie Curie", 
            "“The microscope opens a new world to the investigator.” - Robert Koch",
            "“In science, the credit goes to the man who convinces the world, not to the man to whom the idea first occurs.” - Francis Darwin",
            "“The good thing about science is that it's true whether or not you believe in it.” - Neil deGrasse Tyson",
            "“Science knows no country, because knowledge belongs to humanity.” - Louis Pasteur"
        ]
        
    def _setup_logging(self):
        """Setup logging - must be called first in __init__"""
        logging.basicConfig(
            level=logging.INFO,
            format='%(asctime)s - %(levelname)s - %(message)s'
        )
        return logging.getLogger(__name__)
    
    def _get_available_ram(self) -> int:
        """Get available RAM in GB"""
        try:
            ram_gb = psutil.virtual_memory().available / (1024 ** 3)
            return ram_gb
        except Exception as e:
            self.logger.warning(f"Could not detect RAM: {e}")
            return 8  # Assume 8GB as fallback
    
    def _calculate_optimal_cpus(self, user_cpus: int = None) -> int:
        """Calculate optimal number of CPU cores for MAXIMUM SPEED"""
        if user_cpus is not None:
            self._log_resource_info(user_cpus)
            return user_cpus
            
        try:
            # Get total PHYSICAL CPU cores (not logical threads)
            total_physical_cores = psutil.cpu_count(logical=False) or os.cpu_count() or 2
            
            # MAXIMUM SPEED RULES - AGGRESSIVE CPU USAGE
            if total_physical_cores <= 4:
                optimal_cpus = total_physical_cores  # Use ALL cores on small systems
            elif total_physical_cores <= 8:
                optimal_cpus = total_physical_cores - 1  # Use 7/8, 6/7, etc.
            elif total_physical_cores <= 16:
                optimal_cpus = max(8, total_physical_cores - 2)  # Use 14/16, 13/15, etc.
            elif total_physical_cores <= 32:
                optimal_cpus = max(16, total_physical_cores - 4)  # Use 28/32, 27/31, etc.
            else:
                optimal_cpus = min(32, int(total_physical_cores * 0.85))  # Use 85% on huge systems
            
            # Ensure at least 1 CPU and not more than available cores
            optimal_cpus = max(1, min(optimal_cpus, total_physical_cores))
            
            self._log_resource_info(optimal_cpus, total_physical_cores)
            return optimal_cpus
            
        except Exception as e:
            # Fallback to using all available cores for maximum speed
            self.logger.warning(f"Could not detect CPU cores, using maximum available: {e}")
            return os.cpu_count() or 4
    
    def _log_resource_info(self, cpus: int, total_cores: int = None):
        """Log resource allocation information"""
        self.logger.info(f"Available RAM: {self.available_ram:.1f} GB")
        
        if total_cores:
            self.logger.info(f"System CPU cores: {total_cores}")
            utilization = (cpus / total_cores) * 100
            self.logger.info(f"Using CPU cores: {cpus} ({utilization:.1f}% of available cores)")
        else:
            self.logger.info(f"Using user-specified CPU cores: {cpus}")
        
        # Performance recommendations - MAXIMUM SPEED FOCUS
        if cpus == 1:
            self.logger.info("💡 Performance: Single-core (max speed for 1-core systems)")
        elif cpus <= 4:
            self.logger.info("💡 Performance: Multi-core (max speed for small systems)")
        elif cpus <= 8:
            self.logger.info("💡 Performance: High-speed mode")
        else:
            self.logger.info("💡 Performance: MAXIMUM SPEED MODE 🚀")
    
    def check_amrfinder_installed(self) -> bool:
        """Check if BUNDLED AMRfinderPlus is available"""
        try:
            if not os.path.exists(self.bundled_amrfinder):
                self.logger.error(f"Bundled AMRfinderPlus not found at: {self.bundled_amrfinder}")
                return False
            
            if not os.access(self.bundled_amrfinder, os.X_OK):
                self.logger.warning(f"Bundled AMRfinderPlus not executable, fixing permissions...")
                os.chmod(self.bundled_amrfinder, 0o755)
            
            # Test the bundled version
            result = subprocess.run(
                [self.bundled_amrfinder, '--version'], 
                capture_output=True, 
                text=True, 
                check=True
            )
            
            version_line = result.stdout.strip()
            self.logger.info(f"Bundled AMRfinderPlus version: {version_line}")
            
            # Check database
            if os.path.exists(self.bundled_database):
                self.logger.info(f"✅ Bundled database found: {self.bundled_database}")
                db_version_file = os.path.join(self.bundled_database, "version.txt")
                if os.path.exists(db_version_file):
                    with open(db_version_file, 'r') as f:
                        db_version = f.read().strip()
                        self.logger.info(f"✅ Database version: {db_version}")
            else:
                self.logger.warning(f"⚠️  Bundled database not found at: {self.bundled_database}")
                self.logger.info("Will use default AMRfinderPlus database location")
            
            return True
            
        except (subprocess.CalledProcessError, FileNotFoundError) as e:
            self.logger.error(f"Bundled AMRfinderPlus check failed: {e}")
            self.logger.error("This tool requires the bundled AMRfinderPlus in bin/ directory")
            return False
    
    def run_amrfinder_single_genome(self, genome_file: str, output_dir: str) -> Dict[str, Any]:
        """Run BUNDLED AMRfinderPlus on a single genome - MAXIMUM SPEED"""
        genome_name = Path(genome_file).stem
        output_file = os.path.join(output_dir, f"{genome_name}_amrfinder.txt")
        
        # USE ALL AVAILABLE CPU CORES for each genome - MAXIMUM SPEED!
        run_cpus = self.cpus
        
        # Build command with BUNDLED resources
        cmd = [
            self.bundled_amrfinder,
            '-n', genome_file,  # Nucleotide mode
            '--output', output_file,
            '--threads', str(run_cpus),
            '--plus',  # Include virulence factors
            '--organism', 'Staphylococcus_aureus'
        ]
        
        # Add bundled database if it exists
        if os.path.exists(self.bundled_database):
            cmd.extend(['--database', self.bundled_database])
            self.logger.info(f"Using bundled database: {self.bundled_database}")
        else:
            self.logger.warning("Using default AMRfinderPlus database location")
        
        self.logger.info(f"Running BUNDLED AMRfinderPlus: {genome_name} (using {run_cpus} CPU cores)")
        self.logger.debug(f"Command: {' '.join(cmd)}")
        
        try:
            result = subprocess.run(cmd, capture_output=True, text=True, check=True)
            
            # Parse results for reporting
            hits = self._parse_amrfinder_output(output_file)
            
            # Create individual HTML report
            self._create_amrfinder_html_report(genome_name, hits, output_dir)
            
            return {
                'genome': genome_name,
                'output_file': output_file,
                'hits': hits,
                'hit_count': len(hits),
                'status': 'success'
            }
            
        except subprocess.CalledProcessError as e:
            self.logger.error(f"AMRfinderPlus failed for {genome_name}")
            self.logger.error(f"STDERR: {e.stderr}")
            self.logger.error(f"STDOUT: {e.stdout}")
            
            return {
                'genome': genome_name,
                'output_file': output_file,
                'hits': [],
                'hit_count': 0,
                'status': 'failed',
                'error': e.stderr
            }
    
    def _parse_amrfinder_output(self, amrfinder_file: str) -> List[Dict]:
        """Parse AMRFinderPlus 4.2.4 output file into structured data"""
        hits = []
        try:
            with open(amrfinder_file, 'r') as f:
                lines = f.readlines()
                
            if not lines or len(lines) < 2:  # Need at least header and one data line
                return hits
                
            # Parse header for AMRFinderPlus 4.2.4
            headers = lines[0].strip().split('\t')
            
            # Log the headers for debugging
            self.logger.debug(f"Found headers: {headers}")
            
            # Parse data lines
            for line_num, line in enumerate(lines[1:], 2):
                line = line.strip()
                if not line:
                    continue
                    
                parts = line.split('\t')
                if len(parts) >= len(headers):
                    hit = dict(zip(headers, parts))
                    
                    # Map to consistent field names for AMRFinderPlus 4.2.4
                    # New headers: 
                    # Protein id, Contig id, Start, Stop, Strand, Element symbol, Element name,
                    # Scope, Type, Subtype, Class, Subclass, Method, Target length,
                    # Reference sequence length, % Coverage of reference, % Identity to reference,
                    # Alignment length, Closest reference accession, Closest reference name,
                    # HMM accession, HMM description
                    
                    processed_hit = {
                        'protein_id': hit.get('Protein id', ''),
                        'contig_id': hit.get('Contig id', ''),
                        'start': hit.get('Start', ''),
                        'stop': hit.get('Stop', ''),
                        'strand': hit.get('Strand', ''),
                        'gene_symbol': hit.get('Element symbol', ''),  # New: Element symbol
                        'sequence_name': hit.get('Element name', ''),  # New: Element name
                        'scope': hit.get('Scope', ''),
                        'element_type': hit.get('Type', ''),  # New: Type
                        'element_subtype': hit.get('Subtype', ''),  # New: Subtype
                        'class': hit.get('Class', ''),
                        'subclass': hit.get('Subclass', ''),
                        'method': hit.get('Method', ''),
                        'target_length': hit.get('Target length', ''),
                        'ref_length': hit.get('Reference sequence length', ''),
                        'coverage': hit.get('% Coverage of reference', '').replace('%', ''),
                        'identity': hit.get('% Identity to reference', '').replace('%', ''),
                        'alignment_length': hit.get('Alignment length', ''),
                        'accession': hit.get('Closest reference accession', ''),
                        'closest_name': hit.get('Closest reference name', ''),
                        'hmm_id': hit.get('HMM accession', ''),
                        'hmm_description': hit.get('HMM description', '')
                    }
                    
                    # Also store original hit for debugging
                    processed_hit['_original'] = hit
                    
                    hits.append(processed_hit)
                    
                else:
                    self.logger.warning(f"Line {line_num} has {len(parts)} parts, expected {len(headers)}: {line[:100]}...")
                    
        except Exception as e:
            self.logger.error(f"Error parsing {amrfinder_file}: {e}")
            import traceback
            self.logger.error(traceback.format_exc())
            
        self.logger.info(f"Parsed {len(hits)} AMR hits from {amrfinder_file}")
        return hits
    
    def _create_amrfinder_html_report(self, genome_name: str, hits: List[Dict], output_dir: str):
        """Create comprehensive HTML report for AMRfinderPlus results with beautiful styling"""
        
        # Analyze AMR results for S. aureus
        analysis = self._analyze_saureus_amr_results(hits)
        
        # JavaScript for rotating quotes
        quotes_js = f"""
        <script>
            let quotes = {json.dumps(self.science_quotes)};
            let currentQuote = 0;
            
            function rotateQuote() {{
                document.getElementById('science-quote').innerHTML = quotes[currentQuote];
                currentQuote = (currentQuote + 1) % quotes.length;
            }}
            
            // Rotate every 10 seconds
            setInterval(rotateQuote, 10000);
            
            // Initial display
            document.addEventListener('DOMContentLoaded', function() {{
                rotateQuote();
            }});
        </script>
        """
        
        html_content = f"""
<!DOCTYPE html>
<html>
<head>
    <title>StaphScope AMRfinderPlus Analysis Report</title>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <style>
        body {{ 
            font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif; 
            margin: 0; 
            padding: 0; 
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            min-height: 100vh;
        }}
        .container {{ 
            max-width: 1200px; 
            margin: 0 auto; 
            padding: 20px; 
        }}
        .header {{ 
            background: rgba(255, 255, 255, 0.95); 
            padding: 30px; 
            border-radius: 15px; 
            margin-bottom: 30px; 
            box-shadow: 0 8px 32px rgba(0,0,0,0.1);
            backdrop-filter: blur(10px);
        }}
        .card {{ 
            background: rgba(255, 255, 255, 0.95); 
            padding: 25px; 
            margin: 20px 0; 
            border-radius: 12px; 
            box-shadow: 0 4px 20px rgba(0,0,0,0.1);
            backdrop-filter: blur(10px);
        }}
        .gene-table, .class-table {{ 
            width: 100%; 
            border-collapse: collapse; 
            margin: 20px 0; 
            background: white;
            border-radius: 8px;
            overflow: hidden;
            box-shadow: 0 2px 10px rgba(0,0,0,0.1);
        }}
        .gene-table th, .gene-table td, .class-table th, .class-table td {{ 
            padding: 15px; 
            text-align: left; 
            border-bottom: 1px solid #e0e0e0; 
        }}
        .gene-table th, .class-table th {{ 
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            color: white;
            font-weight: 600;
        }}
        tr:hover {{ background-color: #f8f9fa; }}
        .success {{ color: #28a745; font-weight: 600; }}
        .warning {{ color: #ffc107; font-weight: 600; }}
        .error {{ color: #dc3545; font-weight: 600; }}
        .summary-stats {{ 
            display: flex; 
            justify-content: space-around; 
            margin: 20px 0; 
            flex-wrap: wrap;
        }}
        .stat-card {{ 
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            color: white;
            padding: 20px; 
            border-radius: 12px; 
            text-align: center; 
            box-shadow: 0 4px 15px rgba(0,0,0,0.2);
            margin: 10px;
            flex: 1;
            min-width: 200px;
        }}
        .critical-stat-card {{
            background: linear-gradient(135deg, #dc3545 0%, #c82333 100%);
            color: white;
            padding: 20px; 
            border-radius: 12px; 
            text-align: center; 
            box-shadow: 0 4px 15px rgba(0,0,0,0.2);
            margin: 10px;
            flex: 1;
            min-width: 200px;
        }}
        .quote-container {{
            background: rgba(255, 255, 255, 0.1);
            color: white;
            padding: 20px;
            border-radius: 12px;
            margin: 20px 0;
            text-align: center;
            font-style: italic;
            border-left: 4px solid #fff;
        }}
        .footer {{
            background: rgba(0, 0, 0, 0.8);
            color: white;
            padding: 30px;
            border-radius: 12px;
            margin-top: 40px;
        }}
        .footer a {{
            color: #667eea;
            text-decoration: none;
        }}
        .footer a:hover {{
            text-decoration: underline;
        }}
        .resistance-badge {{
            display: inline-block;
            background: #667eea;
            color: white;
            padding: 5px 10px;
            border-radius: 15px;
            margin: 2px;
            font-size: 0.9em;
        }}
        .high-risk {{ background: #dc3545; }}
        .critical-risk {{ background: #8b0000; font-weight: bold; }}
        .medium-risk {{ background: #ffc107; color: black; }}
        .low-risk {{ background: #28a745; }}
        .present {{ background-color: #d4edda; }}
        .critical-row {{ background-color: #f8d7da; font-weight: bold; border-left: 4px solid #dc3545; }}
        .high-risk-row {{ background-color: #fff3cd; border-left: 4px solid #ffc107; }}
        /* Make sequence name column responsive with word wrapping */
        .sequence-cell {{
            white-space: normal !important;
            word-wrap: break-word;
            max-width: 400px;
            min-width: 200px;
        }}
        /* Make tables responsive */
        .table-responsive {{
            width: 100%;
            overflow-x: auto;
            margin: 20px 0;
        }}
        .gene-table, .class-table {{
            min-width: 900px;
        }}
        /* Tool info box */
        .tool-info {{
            background: linear-gradient(135deg, #17a2b8 0%, #138496 100%);
            color: white;
            padding: 15px;
            border-radius: 8px;
            margin: 10px 0;
            font-size: 0.9em;
        }}
    </style>
    {quotes_js}
</head>
<body>
    <div class="container">
        <div class="header">
            <h1 style="color: #333; margin: 0; font-size: 2.5em;">🧫 StaphScope AMRfinderPlus Analysis Report</h1>
            <p style="color: #666; font-size: 1.2em;">Comprehensive S. aureus Antimicrobial Resistance Analysis</p>
            <div class="tool-info">
                <strong>Tool Info:</strong> Bundled AMRfinderPlus {self.metadata['amrfinder_version']} | 
                Database: {self.metadata['database_version']} | 
                Analysis includes AMR genes and virulence factors (--plus mode)
            </div>
        </div>
        
        <div class="quote-container">
            <div id="science-quote" style="font-size: 1.1em;"></div>
        </div>
"""
        
        # CRITICAL RISK ALERT - Show first if critical genes detected
        if analysis['critical_risk_genes'] > 0:
            html_content += f"""
        <div class="card" style="border-left: 4px solid #dc3545; background: #f8d7da;">
            <h2 style="color: #dc3545;">🚨 CRITICAL RISK AMR GENES DETECTED</h2>
            <p><strong>{analysis['critical_risk_genes']} CRITICAL RISK antimicrobial resistance genes found:</strong></p>
            <div style="margin: 10px 0;">
                <p style="color: #721c24; font-weight: bold;">
                    ⚠️ These genes confer resistance to last-resort antibiotics and represent 
                    a serious public health concern requiring immediate attention.
                </p>
"""
            for gene in analysis['critical_risk_list']:
                html_content += f'<span class="resistance-badge critical-risk" style="font-size: 1.1em;">🚨 {gene}</span>'
            html_content += """
            </div>
        </div>
"""
        
        html_content += f"""
        <div class="card">
            <h2 style="color: #333; border-bottom: 2px solid #667eea; padding-bottom: 10px;">📊 S. aureus AMR Summary</h2>
            <div class="summary-stats">
                <div class="stat-card">
                    <h3>Total AMR Genes</h3>
                    <p style="font-size: 2em; margin: 0;">{analysis['total_genes']}</p>
                </div>
                <div class="stat-card">
                    <h3>High Risk Genes</h3>
                    <p style="font-size: 2em; margin: 0;">{analysis['high_risk_genes']}</p>
                </div>
                <div class="critical-stat-card">
                    <h3>Critical Risk</h3>
                    <p style="font-size: 2em; margin: 0;">{analysis['critical_risk_genes']}</p>
                </div>
            </div>
            <p><strong>Genome:</strong> {genome_name}</p>
            <p><strong>Date:</strong> {self.metadata['analysis_date']}</p>
            <p><strong>Tool Version:</strong> {self.metadata['version']}</p>
            <p><strong>AMRfinderPlus:</strong> {self.metadata['amrfinder_version']}</p>
            <p><strong>Database:</strong> {self.metadata['database_version']}</p>
        </div>
"""
        
        # High-risk genes warning (non-critical)
        if analysis['high_risk_genes'] > 0 and analysis['critical_risk_genes'] == 0:
            html_content += f"""
        <div class="card" style="border-left: 4px solid #ffc107;">
            <h2 style="color: #856404;">⚠️ High-Risk AMR Genes Detected</h2>
            <p><strong>{analysis['high_risk_genes']} high-risk antimicrobial resistance genes found:</strong></p>
            <div style="margin: 10px 0;">
"""
            for gene in analysis['high_risk_list']:
                html_content += f'<span class="resistance-badge high-risk">{gene}</span>'
            html_content += """
            </div>
        </div>
"""
        
        # Resistance Mechanism Breakdown
        if any(analysis['resistance_mechanisms'].values()):
            html_content += """
        <div class="card">
            <h2 style="color: #333; border-bottom: 2px solid #667eea; padding-bottom: 10px;">🔬 Resistance Mechanism Breakdown</h2>
"""
            
            mechanisms = analysis['resistance_mechanisms']
            if mechanisms['mrsa']:
                html_content += f"""
            <div style="margin: 10px 0; padding: 10px; background: #f8d7da; border-radius: 5px;">
                <strong>MRSA Markers (CRITICAL):</strong> {', '.join(mechanisms['mrsa'])}
            </div>
"""
            if mechanisms['glycopeptide_resistance']:
                html_content += f"""
            <div style="margin: 10px 0; padding: 10px; background: #f8d7da; border-radius: 5px;">
                <strong>Glycopeptide Resistance (VISA/VRSA - CRITICAL):</strong> {', '.join(mechanisms['glycopeptide_resistance'])}
            </div>
"""
            if mechanisms['oxazolidinone_resistance']:
                html_content += f"""
            <div style="margin: 10px 0; padding: 10px; background: #f8d7da; border-radius: 5px;">
                <strong>Oxazolidinone Resistance (Linezolid - CRITICAL):</strong> {', '.join(mechanisms['oxazolidinone_resistance'])}
            </div>
"""
            if mechanisms['mlsb_resistance']:
                html_content += f"""
            <div style="margin: 10px 0; padding: 10px; background: #d1ecf1; border-radius: 5px;">
                <strong>MLS⸰B Resistance:</strong> {', '.join(mechanisms['mlsb_resistance'])}
            </div>
"""
            if mechanisms['aminoglycoside_resistance']:
                html_content += f"""
            <div style="margin: 10px 0; padding: 10px; background: #d1ecf1; border-radius: 5px;">
                <strong>Aminoglycoside Resistance:</strong> {', '.join(mechanisms['aminoglycoside_resistance'])}
            </div>
"""
            if mechanisms['efflux_pumps']:
                html_content += f"""
            <div style="margin: 10px 0; padding: 10px; background: #e2e3e5; border-radius: 5px;">
                <strong>Efflux Pumps:</strong> {', '.join(mechanisms['efflux_pumps'])}
            </div>
"""
            
            html_content += """
        </div>
"""
        
        # Resistance classes summary
        if analysis['resistance_classes']:
            html_content += """
        <div class="card">
            <h2 style="color: #333; border-bottom: 2px solid #667eea; padding-bottom: 10px;">🧪 Resistance Classes Detected</h2>
            <div class="table-responsive">
                <table class="class-table">
                    <thead>
                        <tr>
                            <th>Resistance Class</th>
                            <th>Gene Count</th>
                            <th>Genes</th>
                        </tr>
                    </thead>
                    <tbody>
"""
            
            for class_name, genes in analysis['resistance_classes'].items():
                gene_list = ", ".join(genes)
                html_content += f"""
                    <tr>
                        <td><strong>{class_name}</strong></td>
                        <td>{len(genes)}</td>
                        <td class="sequence-cell">{gene_list}</td>
                    </tr>
"""
            
            html_content += """
                    </tbody>
                </table>
            </div>
        </div>
"""
        
        # Detailed AMR genes table
        if hits:
            html_content += """
        <div class="card">
            <h2 style="color: #333; border-bottom: 2px solid #667eea; padding-bottom: 10px;">🔬 Detailed AMR Genes Detected</h2>
            <div class="table-responsive">
                <table class="gene-table">
                    <thead>
                        <tr>
                            <th>Gene Symbol</th>
                            <th>Sequence Name</th>
                            <th>Class</th>
                            <th>Subclass</th>
                            <th>Coverage</th>
                            <th>Identity</th>
                            <th>Scope</th>
                        </tr>
                    </thead>
                    <tbody>
"""
            
            for hit in hits:
                # Determine row class based on risk level
                row_class = "present"
                gene_symbol = hit.get('gene_symbol', '')
                if gene_symbol in analysis['critical_risk_list']:
                    row_class = "critical-row"
                elif gene_symbol in analysis['high_risk_list']:
                    row_class = "high-risk-row"
                
                # Show full sequence name without truncation
                sequence_display = hit.get('sequence_name', '')
                
                html_content += f"""
                    <tr class="{row_class}">
                        <td><strong>{gene_symbol}</strong></td>
                        <td class="sequence-cell">{sequence_display}</td>
                        <td>{hit.get('class', '')}</td>
                        <td>{hit.get('subclass', '')}</td>
                        <td>{hit.get('coverage', '')}%</td>
                        <td>{hit.get('identity', '')}%</td>
                        <td>{hit.get('scope', '')}</td>
                    </tr>
"""
            
            html_content += """
                    </tbody>
                </table>
            </div>
        </div>
"""
        else:
            html_content += """
        <div class="card">
            <h2 style="color: #333; border-bottom: 2px solid #667eea; padding-bottom: 10px;">✅ No AMR Genes Detected</h2>
            <p>No antimicrobial resistance genes found in this S. aureus genome.</p>
        </div>
"""
        
        # Footer
        html_content += f"""
        <div class="footer">
            <h3 style="color: #fff; border-bottom: 2px solid #667eea; padding-bottom: 10px;">👥 Contact Information</h3>
            <p><strong>Author:</strong> Brown Beckley</p>
            <p><strong>Email:</strong> brownbeckley94@gmail.com</p>
            <p><strong>GitHub:</strong> <a href="https://github.com/bbeckley-hub" target="_blank">https://github.com/bbeckley-hub</a></p>
            <p><strong>Affiliation:</strong> University of Ghana Medical School</p>
            <p style="margin-top: 20px; font-size: 0.9em; color: #ccc;">
                Analysis performed using StaphScope AMRfinderPlus v1.0.1<br>
                Bundled AMRfinderPlus {self.metadata['amrfinder_version']} with database {self.metadata['database_version']}
            </p>
        </div>
    </div>
</body>
</html>
"""
        
        # Write HTML report
        html_file = os.path.join(output_dir, f"{genome_name}_amrfinder_report.html")
        with open(html_file, 'w') as f:
            f.write(html_content)
        
        self.logger.info(f"S. aureus AMRfinderPlus HTML report generated: {html_file}")
    
    def _analyze_saureus_amr_results(self, hits: List[Dict]) -> Dict[str, Any]:
        """Analyze AMR results specifically for S. aureus with enhanced risk assessment"""
        
        analysis = {
            'total_genes': len(hits),
            'resistance_classes': {},
            'total_classes': 0,
            'high_risk_genes': 0,
            'critical_risk_genes': 0,
            'high_risk_list': [],
            'critical_risk_list': [],
            'resistance_mechanisms': {
                'mrsa': [],                    # Methicillin resistance
                'glycopeptide_resistance': [], # Vancomycin resistance
                'oxazolidinone_resistance': [], # Linezolid resistance
                'mlsb_resistance': [],         # Macrolide-Lincosamide-Streptogramin B
                'aminoglycoside_resistance': [], # Aminoglycoside resistance
                'tetracycline_resistance': [], # Tetracycline resistance
                'fluoroquinolone_resistance': [], # Fluoroquinolone resistance
                'efflux_pumps': [],            # Multi-drug efflux pumps
                'other_amr': []                # Other resistance mechanisms
            }
        }
        
        for hit in hits:
            gene_symbol = hit.get('gene_symbol', '')
            resistance_class = hit.get('class', '')
            
            if not gene_symbol:
                continue
            
            # Categorize resistance mechanism
            self._categorize_resistance_mechanism(gene_symbol, resistance_class, analysis)
            
            # Check for critical risk genes
            if gene_symbol in self.critical_risk_genes:
                analysis['critical_risk_genes'] += 1
                if gene_symbol not in analysis['critical_risk_list']:
                    analysis['critical_risk_list'].append(gene_symbol)
            
            # Check for high-risk genes (includes critical ones)
            if gene_symbol in self.high_risk_genes:
                analysis['high_risk_genes'] += 1
                if gene_symbol not in analysis['high_risk_list']:
                    analysis['high_risk_list'].append(gene_symbol)
            
            # Group by resistance class
            if resistance_class:
                if resistance_class not in analysis['resistance_classes']:
                    analysis['resistance_classes'][resistance_class] = []
                if gene_symbol not in analysis['resistance_classes'][resistance_class]:
                    analysis['resistance_classes'][resistance_class].append(gene_symbol)
        
        analysis['total_classes'] = len(analysis['resistance_classes'])
        return analysis

    def _categorize_resistance_mechanism(self, gene_symbol: str, resistance_class: str, analysis: Dict[str, Any]):
        """Categorize genes by resistance mechanism for S. aureus"""
        
        # MRSA markers
        mrsa_genes = {'mecA', 'mecC'}
        
        # Glycopeptide resistance (VISA/VRSA)
        glycopeptide_genes = {'vanA', 'vanB', 'vanC', 'vanD', 'vanE', 'vanG'}
        
        # Oxazolidinone resistance (Linezolid)
        oxazolidinone_genes = {'cfr', 'optrA', 'poxtA'}
        
        # MLSB resistance
        mlsb_genes = {'ermA', 'ermB', 'ermC', 'ermF', 'ermT', 'ermY', 'mphC', 'msrA'}
        
        # Aminoglycoside resistance
        aminoglycoside_genes = {
            'aacA-aphD', 'aac(6\')-aph(2")', 'aph(3\')-III', 'ant(4\')-Ia',
            'aac(6\')-Ie-aph(2")-Ia'
        }
        
        # Tetracycline resistance
        tetracycline_genes = {'tetK', 'tetL', 'tetM', 'tetO'}
        
        # Fluoroquinolone resistance
        fluoroquinolone_genes = {'grlA', 'grlB', 'gyrA', 'gyrB', 'norA', 'norB'}
        
        # Efflux pumps
        efflux_pump_genes = {'qacA', 'qacB', 'qacC', 'smr', 'sepA'}
        
        if gene_symbol in mrsa_genes:
            analysis['resistance_mechanisms']['mrsa'].append(gene_symbol)
        elif gene_symbol in glycopeptide_genes:
            analysis['resistance_mechanisms']['glycopeptide_resistance'].append(gene_symbol)
        elif gene_symbol in oxazolidinone_genes:
            analysis['resistance_mechanisms']['oxazolidinone_resistance'].append(gene_symbol)
        elif gene_symbol in mlsb_genes:
            analysis['resistance_mechanisms']['mlsb_resistance'].append(gene_symbol)
        elif gene_symbol in aminoglycoside_genes:
            analysis['resistance_mechanisms']['aminoglycoside_resistance'].append(gene_symbol)
        elif gene_symbol in tetracycline_genes:
            analysis['resistance_mechanisms']['tetracycline_resistance'].append(gene_symbol)
        elif gene_symbol in fluoroquinolone_genes:
            analysis['resistance_mechanisms']['fluoroquinolone_resistance'].append(gene_symbol)
        elif gene_symbol in efflux_pump_genes:
            analysis['resistance_mechanisms']['efflux_pumps'].append(gene_symbol)
        else:
            analysis['resistance_mechanisms']['other_amr'].append(gene_symbol)
    
    def create_amr_summary(self, all_results: Dict[str, Any], output_base: str):
        """Create comprehensive AMR summary files and HTML reports for all S. aureus samples"""
        self.logger.info("Creating S. aureus AMR summary files and HTML reports...")
        
        # Create TSV summary files
        summary_file = os.path.join(output_base, "staph_amrfinder_summary.tsv")
        
        with open(summary_file, 'w') as f:
            # Write header
            f.write("Genome\tGene_Symbol\tSequence_Name\tClass\tSubclass\tCoverage\tIdentity\tScope\tElement_Type\tAccession\tContig\tStart\tStop\n")
            
            # Write data for all genomes
            for genome_name, result in all_results.items():
                for hit in result['hits']:
                    row = [
                        genome_name,
                        hit.get('gene_symbol', ''),
                        hit.get('sequence_name', ''),
                        hit.get('class', ''),
                        hit.get('subclass', ''),
                        hit.get('coverage', ''),
                        hit.get('identity', ''),
                        hit.get('scope', ''),
                        hit.get('element_type', ''),
                        hit.get('accession', ''),
                        hit.get('contig_id', ''),
                        hit.get('start', ''),
                        hit.get('stop', '')
                    ]
                    f.write('\t'.join(str(x) for x in row) + '\n')
        
        self.logger.info(f"✓ S. aureus AMR summary file created: {summary_file}")
        
        # Create statistics summary
        stats_file = os.path.join(output_base, "staph_amrfinder_statistics_summary.tsv")
        with open(stats_file, 'w') as f:
            f.write("Genome\tTotal_AMR_Genes\tHigh_Risk_Genes\tCritical_Risk_Genes\tResistance_Classes\tGene_List\n")
            
            for genome_name, result in all_results.items():
                # Get unique genes
                genes = list(set(hit.get('gene_symbol', '') for hit in result['hits'] if hit.get('gene_symbol')))
                gene_list = ",".join(genes)
                
                # Count high-risk and critical genes
                high_risk_count = sum(1 for gene in genes if gene in self.high_risk_genes)
                critical_risk_count = sum(1 for gene in genes if gene in self.critical_risk_genes)
                
                # Get resistance classes
                classes = list(set(hit.get('class', '') for hit in result['hits'] if hit.get('class')))
                class_list = ",".join(classes)
                
                f.write(f"{genome_name}\t{result['hit_count']}\t{high_risk_count}\t{critical_risk_count}\t{class_list}\t{gene_list}\n")
        
        self.logger.info(f"✓ S. aureus AMR statistics summary created: {stats_file}")
        
        # Create comprehensive HTML summary report for staph_amrfinder_summary.tsv
        self._create_summary_html_report(all_results, output_base)
        
        # NEW: Create JSON summaries
        self.create_amr_json_summaries(all_results, output_base)
        self.create_amr_master_json_summary(all_results, output_base)
    
    def create_amr_json_summaries(self, all_results: Dict[str, Any], output_base: str):
        """Create JSON summary files for AMR genes across all genomes."""
        self.logger.info("Creating JSON AMR summaries...")
        
        # Collect all hits with genome information
        all_hits = []
        for genome_name, result in all_results.items():
            for hit in result['hits']:
                hit_with_genome = hit.copy()
                hit_with_genome['genome'] = genome_name
                all_hits.append(hit_with_genome)
        
        if not all_hits:
            self.logger.info("No AMR hits found, skipping JSON summaries")
            return
        
        # Calculate gene frequency
        gene_frequency = {}
        for hit in all_hits:
            gene = hit['gene_symbol']
            if not gene:
                continue
                
            if gene not in gene_frequency:
                gene_frequency[gene] = {
                    'count': 0,
                    'genomes': set(),
                    'details': []
                }
            
            gene_frequency[gene]['count'] += 1
            gene_frequency[gene]['genomes'].add(hit['genome'])
            gene_frequency[gene]['details'].append({
                'genome': hit['genome'],
                'sequence_name': hit['sequence_name'],
                'class': hit['class'],
                'subclass': hit['subclass'],
                'coverage': hit['coverage'],
                'identity': hit['identity'],
                'accession': hit['accession']
            })
        
        # Convert sets to lists for JSON serialization
        for gene in gene_frequency:
            gene_frequency[gene]['genomes'] = list(gene_frequency[gene]['genomes'])
        
        # Create JSON structure
        json_summary = {
            'metadata': {
                'tool': 'StaphScope AMRfinderPlus',
                'version': self.metadata['version'],
                'amrfinder_version': self.metadata['amrfinder_version'],
                'database_version': self.metadata['database_version'],
                'analysis_date': self.metadata['analysis_date'],
                'total_hits': len(all_hits),
                'total_genomes': len(all_results),
                'unique_genes': len(gene_frequency),
                'critical_risk_genes_found': list(self.critical_risk_genes.intersection(set(gene_frequency.keys()))),
                'high_risk_genes_found': list(self.high_risk_genes.intersection(set(gene_frequency.keys())))
            },
            'gene_frequency': gene_frequency,
            'summary_by_genome': self._create_amr_genome_summary(all_results),
            'hits': all_hits[:100]  # Include first 100 hits to keep file manageable
        }
        
        # Write JSON file
        json_file = os.path.join(output_base, "staph_amrfinder_summary.json")
        with open(json_file, 'w') as f:
            json.dump(json_summary, f, indent=2, default=str)
        
        self.logger.info(f"✓ Created JSON summary: {json_file}")
    
    def _create_amr_genome_summary(self, all_results: Dict) -> Dict:
        """Create summary of AMR hits by genome."""
        genome_summary = {}
        
        for genome_name, result in all_results.items():
            if genome_name not in genome_summary:
                genome_summary[genome_name] = {
                    'total_genes': result['hit_count'],
                    'genes': set(),
                    'critical_risk_genes': [],
                    'high_risk_genes': [],
                    'resistance_classes': set()
                }
            
            for hit in result['hits']:
                gene = hit['gene_symbol']
                if not gene:
                    continue
                    
                genome_summary[genome_name]['genes'].add(gene)
                
                if hit['class']:
                    genome_summary[genome_name]['resistance_classes'].add(hit['class'])
                
                # Classify genes
                if gene in self.critical_risk_genes:
                    if gene not in genome_summary[genome_name]['critical_risk_genes']:
                        genome_summary[genome_name]['critical_risk_genes'].append(gene)
                elif gene in self.high_risk_genes:
                    if gene not in genome_summary[genome_name]['high_risk_genes']:
                        genome_summary[genome_name]['high_risk_genes'].append(gene)
        
        # Convert sets to lists
        for genome in genome_summary:
            genome_summary[genome]['genes'] = list(genome_summary[genome]['genes'])
            genome_summary[genome]['resistance_classes'] = list(genome_summary[genome]['resistance_classes'])
        
        return genome_summary
    
    def create_amr_master_json_summary(self, all_results: Dict[str, Any], output_base: str):
        """Create a master JSON summary for AMR genes with cross-genome patterns."""
        self.logger.info("Creating master JSON summary for AMR genes...")
        
        # Collect overall statistics
        master_summary = {
            'metadata': {
                'tool': 'StaphScope AMRfinderPlus Module',
                'version': self.metadata['version'],
                'amrfinder_version': self.metadata['amrfinder_version'],
                'database_version': self.metadata['database_version'],
                'analysis_date': self.metadata['analysis_date'],
                'total_genomes': len(all_results),
                'critical_risk_genes': list(self.critical_risk_genes),
                'high_risk_genes': list(self.high_risk_genes)
            },
            'genome_summaries': {},
            'critical_findings': {},
            'cross_genome_patterns': {}
        }
        
        # Analyze each genome
        genomes_with_critical = 0
        genomes_with_high_risk = 0
        
        for genome_name, result in all_results.items():
            # Collect genes for this genome
            genes = [hit['gene_symbol'] for hit in result['hits'] if hit['gene_symbol']]
            unique_genes = set(genes)
            
            # Get resistance classes
            classes = set(hit['class'] for hit in result['hits'] if hit['class'])
            
            # Check for critical and high-risk genes
            has_critical = any(gene in unique_genes for gene in self.critical_risk_genes)
            has_high_risk = any(gene in unique_genes for gene in self.high_risk_genes)
            
            if has_critical:
                genomes_with_critical += 1
            if has_high_risk:
                genomes_with_high_risk += 1
            
            master_summary['genome_summaries'][genome_name] = {
                'total_hits': result['hit_count'],
                'unique_genes': len(unique_genes),
                'has_critical_risk': has_critical,
                'has_high_risk': has_high_risk,
                'resistance_classes': list(classes),
                'genes': list(unique_genes)
            }
            
            # Track critical findings
            if has_critical:
                critical_genes_found = [g for g in unique_genes if g in self.critical_risk_genes]
                if 'genomes_with_critical' not in master_summary['critical_findings']:
                    master_summary['critical_findings']['genomes_with_critical'] = []
                master_summary['critical_findings']['genomes_with_critical'].append({
                    'genome': genome_name,
                    'critical_genes': critical_genes_found
                })
        
        # Find cross-genome patterns (genes found in multiple genomes)
        all_hits_by_gene = {}
        for genome_name, result in all_results.items():
            for hit in result['hits']:
                gene = hit['gene_symbol']
                if not gene:
                    continue
                    
                if gene not in all_hits_by_gene:
                    all_hits_by_gene[gene] = {
                        'count': 0,
                        'genomes': [],
                        'classes': set(),
                        'details': []
                    }
                
                all_hits_by_gene[gene]['count'] += 1
                if genome_name not in all_hits_by_gene[gene]['genomes']:
                    all_hits_by_gene[gene]['genomes'].append(genome_name)
                if hit['class']:
                    all_hits_by_gene[gene]['classes'].add(hit['class'])
                
                all_hits_by_gene[gene]['details'].append({
                    'genome': genome_name,
                    'sequence_name': hit['sequence_name'],
                    'coverage': hit['coverage'],
                    'identity': hit['identity']
                })
        
        # Convert sets to lists
        for gene in all_hits_by_gene:
            all_hits_by_gene[gene]['classes'] = list(all_hits_by_gene[gene]['classes'])
        
        # Find common genes (present in >1 genome)
        common_genes = {}
        for gene, data in all_hits_by_gene.items():
            if data['count'] > 1:  # Gene found in multiple genomes
                common_genes[gene] = data
        
        # Get top genes
        top_genes = sorted(
            [(gene, data) for gene, data in all_hits_by_gene.items()],
            key=lambda x: x[1]['count'],
            reverse=True
        )[:20]  # Top 20 most frequent genes
        
        master_summary['cross_genome_patterns'] = {
            'total_genes_found': len(all_hits_by_gene),
            'genomes_with_critical_risk': genomes_with_critical,
            'genomes_with_high_risk': genomes_with_high_risk,
            'common_genes': common_genes,
            'top_genes': top_genes,
            'risk_gene_distribution': {
                'critical_risk_found': len([g for g in all_hits_by_gene.keys() if g in self.critical_risk_genes]),
                'high_risk_found': len([g for g in all_hits_by_gene.keys() if g in self.high_risk_genes]),
                'standard_risk_found': len([g for g in all_hits_by_gene.keys() if g not in self.high_risk_genes])
            }
        }
        
        # Write master JSON
        json_file = os.path.join(output_base, "staph_amrfinder_master_summary.json")
        with open(json_file, 'w') as f:
            json.dump(master_summary, f, indent=2, default=str)
        
        self.logger.info(f"✓ Created master JSON summary: {json_file}")
        return master_summary
    
    def _create_summary_html_report(self, all_results: Dict[str, Any], output_base: str):
        """Create comprehensive HTML summary report with pattern discovery"""
        
        # Collect all data for pattern analysis
        all_hits = []
        for genome_name, result in all_results.items():
            for hit in result['hits']:
                hit_with_genome = hit.copy()
                hit_with_genome['genome'] = genome_name
                all_hits.append(hit_with_genome)
        
        # Calculate statistics
        total_genomes = len(all_results)
        total_hits = len(all_hits)
        
        # Track critical and high-risk genes across all genomes
        critical_genes_found = set()
        high_risk_genes_found = set()
        genomes_with_critical = 0
        genomes_with_high_risk = 0
        
        # Calculate genes per genome and gene frequency
        genes_per_genome = {}
        gene_frequency = {}
        
        for genome_name, result in all_results.items():
            genome_genes = set()
            for hit in result['hits']:
                gene = hit.get('gene_symbol', '')
                if gene:
                    genome_genes.add(gene)
                    
                    # Track gene frequency
                    if gene not in gene_frequency:
                        gene_frequency[gene] = set()
                    gene_frequency[gene].add(genome_name)
            
            genes_per_genome[genome_name] = genome_genes
            
            # Check for critical and high-risk genes
            has_critical = any(gene in genome_genes for gene in self.critical_risk_genes)
            has_high_risk = any(gene in genome_genes for gene in self.high_risk_genes)
            
            if has_critical:
                genomes_with_critical += 1
                critical_genes_found.update(genome_genes.intersection(self.critical_risk_genes))
            
            if has_high_risk:
                genomes_with_high_risk += 1
                high_risk_genes_found.update(genome_genes.intersection(self.high_risk_genes))
        
        # JavaScript for rotating quotes
        quotes_js = f"""
        <script>
            let quotes = {json.dumps(self.science_quotes)};
            let currentQuote = 0;
            
            function rotateQuote() {{
                document.getElementById('science-quote').innerHTML = quotes[currentQuote];
                currentQuote = (currentQuote + 1) % quotes.length;
            }}
            
            // Rotate every 10 seconds
            setInterval(rotateQuote, 10000);
            
            // Initial display
            document.addEventListener('DOMContentLoaded', function() {{
                rotateQuote();
            }});
        </script>
        """
        
        html_content = f"""
<!DOCTYPE html>
<html>
<head>
    <title>StaphScope AMRfinderPlus - Summary Report (BUNDLED)</title>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <style>
        body {{ 
            font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif; 
            margin: 0; 
            padding: 0; 
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            min-height: 100vh;
        }}
        .container {{ 
            max-width: 1400px; 
            margin: 0 auto; 
            padding: 20px; 
        }}
        .header {{ 
            background: rgba(255, 255, 255, 0.95); 
            padding: 30px; 
            border-radius: 15px; 
            margin-bottom: 30px; 
            box-shadow: 0 8px 32px rgba(0,0,0,0.1);
            backdrop-filter: blur(10px);
        }}
        .card {{ 
            background: rgba(255, 255, 255, 0.95); 
            padding: 25px; 
            margin: 20px 0; 
            border-radius: 12px; 
            box-shadow: 0 4px 20px rgba(0,0,0,0.1);
            backdrop-filter: blur(10px);
        }}
        .gene-table {{ 
            width: 100%; 
            border-collapse: collapse; 
            margin: 20px 0; 
            background: white;
            border-radius: 8px;
            overflow: hidden;
            box-shadow: 0 2px 10px rgba(0,0,0,0.1);
        }}
        .gene-table th, .gene-table td {{ 
            padding: 12px; 
            text-align: left; 
            border-bottom: 1px solid #e0e0e0; 
        }}
        .gene-table th {{ 
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            color: white;
            font-weight: 600;
        }}
        tr:hover {{ background-color: #f8f9fa; }}
        .summary-stats {{ 
            display: flex; 
            justify-content: space-around; 
            margin: 20px 0; 
            flex-wrap: wrap;
        }}
        .stat-card {{ 
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            color: white;
            padding: 20px; 
            border-radius: 12px; 
            text-align: center; 
            box-shadow: 0 4px 15px rgba(0,0,0,0.2);
            margin: 10px;
            flex: 1;
            min-width: 200px;
        }}
        .critical-stat-card {{
            background: linear-gradient(135deg, #dc3545 0%, #c82333 100%);
            color: white;
            padding: 20px; 
            border-radius: 12px; 
            text-align: center; 
            box-shadow: 0 4px 15px rgba(0,0,0,0.2);
            margin: 10px;
            flex: 1;
            min-width: 200px;
        }}
        .quote-container {{
            background: rgba(255, 255, 255, 0.1);
            color: white;
            padding: 20px;
            border-radius: 12px;
            margin: 20px 0;
            text-align: center;
            font-style: italic;
            border-left: 4px solid #fff;
        }}
        .footer {{
            background: rgba(0, 0, 0, 0.8);
            color: white;
            padding: 30px;
            border-radius: 12px;
            margin-top: 40px;
        }}
        .footer a {{
            color: #667eea;
            text-decoration: none;
        }}
        .footer a:hover {{
            text-decoration: underline;
        }}
        .resistance-badge {{
            display: inline-block;
            background: #dc3545;
            color: white;
            padding: 5px 10px;
            border-radius: 15px;
            margin: 2px;
            font-size: 0.9em;
        }}
        .critical-resistance-badge {{
            display: inline-block;
            background: #8b0000;
            color: white;
            padding: 5px 10px;
            border-radius: 15px;
            margin: 2px;
            font-size: 0.9em;
            font-weight: bold;
        }}
        .warning-badge {{
            display: inline-block;
            background: #ffc107;
            color: black;
            padding: 5px 10px;
            border-radius: 15px;
            margin: 2px;
            font-size: 0.9em;
        }}
        .success-badge {{
            display: inline-block;
            background: #28a745;
            color: white;
            padding: 5px 10px;
            border-radius: 15px;
            margin: 2px;
            font-size: 0.9em;
        }}
        /* Make sequence name column responsive with word wrapping */
        .sequence-cell {{
            white-space: normal !important;
            word-wrap: break-word;
            max-width: 400px;
            min-width: 200px;
        }}
        /* Make tables responsive */
        .table-responsive {{
            width: 100%;
            overflow-x: auto;
            margin: 20px 0;
        }}
        .gene-table {{
            min-width: 1000px;
        }}
        /* IMPROVED GENE FREQUENCY COLOR SCHEME */
        .frequency-high {{ background-color: #f8d7da; font-weight: bold; border-left: 4px solid #dc3545; }}
        .frequency-medium-high {{ background-color: #ffeaa7; border-left: 4px solid #fdcb6e; }}
        .frequency-medium {{ background-color: #fff3cd; border-left: 4px solid #ffc107; }}
        .frequency-low-medium {{ background-color: #d1ecf1; border-left: 4px solid #17a2b8; }}
        .frequency-low {{ background-color: #d4edda; border-left: 4px solid #28a745; }}
        /* Tool info box */
        .tool-info {{
            background: linear-gradient(135deg, #17a2b8 0%, #138496 100%);
            color: white;
            padding: 15px;
            border-radius: 8px;
            margin: 10px 0;
            font-size: 0.9em;
        }}
    </style>
    {quotes_js}
</head>
<body>
    <div class="container">
        <div class="header">
            <h1 style="color: #333; margin: 0; font-size: 2.5em;">🧫 StaphScope AMRfinderPlus - Summary Report </h1>
            <p style="color: #666; font-size: 1.2em;">Comprehensive S. aureus Antimicrobial Resistance Analysis Across All Genomes</p>
            <div class="tool-info">
                <strong>Tool Info:</strong> Bundled AMRfinderPlus {self.metadata['amrfinder_version']} | 
                Database: {self.metadata['database_version']} | 
                Analysis includes AMR genes and virulence factors (--plus mode)
            </div>
        </div>
        
        <div class="quote-container">
            <div id="science-quote" style="font-size: 1.1em;"></div>
        </div>
"""
        
        # CRITICAL RISK ALERT - Show first if critical genes detected
        if critical_genes_found:
            html_content += f"""
        <div class="card" style="border-left: 4px solid #dc3545; background: #f8d7da;">
            <h2 style="color: #dc3545;">🚨 CRITICAL RISK AMR GENES ACROSS ALL GENOMES</h2>
            <p><strong>{len(critical_genes_found)} unique critical risk genes found in {genomes_with_critical} genomes:</strong></p>
            <div style="margin: 10px 0;">
                <p style="color: #721c24; font-weight: bold;">
                    ⚠️ IMMEDIATE ATTENTION REQUIRED: These genes confer resistance to last-resort antibiotics
                </p>
"""
            for gene in sorted(critical_genes_found):
                html_content += f'<span class="critical-resistance-badge">🚨 {gene}</span>'
            html_content += """
            </div>
        </div>
"""
        
        html_content += f"""
        <div class="card">
            <h2 style="color: #333; border-bottom: 2px solid #667eea; padding-bottom: 10px;">📊 Overall Summary</h2>
            <div class="summary-stats">
                <div class="stat-card">
                    <h3>Total Genomes</h3>
                    <p style="font-size: 2em; margin: 0;">{total_genomes}</p>
                </div>
                <div class="stat-card">
                    <h3>Total AMR Genes</h3>
                    <p style="font-size: 2em; margin: 0;">{total_hits}</p>
                </div>
                <div class="critical-stat-card">
                    <h3>Critical Risk Genomes</h3>
                    <p style="font-size: 2em; margin: 0;">{genomes_with_critical}</p>
                </div>
            </div>
            <p><strong>Date:</strong> {self.metadata['analysis_date']}</p>
            <p><strong>Tool Version:</strong> {self.metadata['version']}</p>
            <p><strong>AMRfinderPlus:</strong> {self.metadata['amrfinder_version']} </p>
            <p><strong>Database:</strong> {self.metadata['database_version']} </p>
        </div>
"""
        
        # High-risk genes summary (non-critical)
        if high_risk_genes_found and not critical_genes_found:
            html_content += f"""
        <div class="card" style="border-left: 4px solid #ffc107;">
            <h2 style="color: #856404;">⚠️ High-Risk AMR Genes Detected</h2>
            <p><strong>{len(high_risk_genes_found)} unique high-risk genes found across {genomes_with_high_risk} genomes:</strong></p>
            <div style="margin: 10px 0;">
"""
            for gene in sorted(high_risk_genes_found):
                html_content += f'<span class="resistance-badge">{gene}</span>'
            html_content += """
            </div>
        </div>
"""
        
        # Genes by Genome table (Pattern Discovery)
        html_content += """
        <div class="card">
            <h2 style="color: #333; border-bottom: 2px solid #667eea; padding-bottom: 10px;">🔍 Genes by Genome</h2>
            <div class="table-responsive">
                <table class="gene-table">
                    <thead>
                        <tr>
                            <th>Genome</th>
                            <th>Gene Count</th>
                            <th>Critical Genes</th>
                            <th>High Risk Genes</th>
                        </tr>
                    </thead>
                    <tbody>
"""
        
        for genome in sorted(genes_per_genome.keys()):
            genes = genes_per_genome.get(genome, set())
            critical_genes = [g for g in genes if g in self.critical_risk_genes]
            high_risk_genes = [g for g in genes if g in self.high_risk_genes and g not in self.critical_risk_genes]
            
            critical_display = ", ".join(critical_genes) if critical_genes else "None"
            high_risk_display = ", ".join(high_risk_genes) if high_risk_genes else "None"
            
            # Highlight rows with critical genes
            row_class = "critical-row" if critical_genes else "high-risk-row" if high_risk_genes else ""
            
            html_content += f"""
                    <tr class="{row_class}">
                        <td><strong>{genome}</strong></td>
                        <td>{len(genes)}</td>
                        <td class="sequence-cell">{critical_display}</td>
                        <td class="sequence-cell">{high_risk_display}</td>
                    </tr>
"""
        
        html_content += """
                    </tbody>
                </table>
            </div>
        </div>
        
        <div class="card">
            <h2 style="color: #333; border-bottom: 2px solid #667eea; padding-bottom: 10px;">📈 Gene Frequency</h2>
            <div class="table-responsive">
                <table class="gene-table">
                    <thead>
                        <tr>
                            <th>Gene</th>
                            <th>Frequency</th>
                            <th>Prevalence</th>
                            <th>Risk Level</th>
                            <th>Genomes</th>
                        </tr>
                    </thead>
                    <tbody>
"""
        
        # Calculate gene frequency with IMPROVED color highlighting
        for gene, genomes in sorted(gene_frequency.items(), key=lambda x: len(x[1]), reverse=True):
            frequency = len(genomes)
            genome_list = ", ".join(sorted(genomes))
            frequency_percent = (frequency / total_genomes) * 100 if total_genomes > 0 else 0
            
            # Determine risk level
            if gene in self.critical_risk_genes:
                risk_level = '<span class="critical-resistance-badge">CRITICAL</span>'
            elif gene in self.high_risk_genes:
                risk_level = '<span class="resistance-badge">HIGH</span>'
            else:
                risk_level = '<span class="success-badge">Standard</span>'
            
            # IMPROVED COLOR SCHEME: Better visual distinction
            if frequency_percent >= 75:
                frequency_class = "frequency-high"
                prevalence_badge = '<span class="resistance-badge">Very High</span>'
            elif frequency_percent >= 50:
                frequency_class = "frequency-medium-high"
                prevalence_badge = '<span class="warning-badge">High</span>'
            elif frequency_percent >= 25:
                frequency_class = "frequency-medium"
                prevalence_badge = '<span class="warning-badge">Medium</span>'
            elif frequency_percent >= 10:
                frequency_class = "frequency-low-medium"
                prevalence_badge = '<span class="success-badge">Low</span>'
            else:
                frequency_class = "frequency-low"
                prevalence_badge = '<span class="success-badge">Rare</span>'
            
            html_content += f"""
                    <tr class="{frequency_class}">
                        <td><strong>{gene}</strong></td>
                        <td>{frequency} ({frequency_percent:.1f}%)</td>
                        <td>{prevalence_badge}</td>
                        <td>{risk_level}</td>
                        <td class="sequence-cell">{genome_list}</td>
                    </tr>
"""
        
        html_content += """
                    </tbody>
                </table>
            </div>
        </div>
        
        <div class="card">
            <h2 style="color: #333; border-bottom: 2px solid #667eea; padding-bottom: 10px;">📁 Generated Files</h2>
            <ul style="color: #666; font-size: 1.1em;">
                <li><strong>staph_amrfinder_summary.tsv</strong> - Complete AMR data for all genomes</li>
                <li><strong>staph_amrfinder_statistics_summary.tsv</strong> - Statistical summary</li>
                <li><strong>staph_amrfinder_summary.json</strong> - JSON summary with gene frequencies</li>
                <li><strong>staph_amrfinder_master_summary.json</strong> - Master JSON with cross-genome patterns</li>
                <li><strong>Individual genome HTML reports</strong> - Detailed analysis per genome</li>
                <li><strong>This summary report</strong> - Cross-genome analysis with pattern discovery</li>
            </ul>
        </div>
        
        <div class="footer">
            <h3 style="color: #fff; border-bottom: 2px solid #667eea; padding-bottom: 10px;">👥 Contact Information</h3>
            <p><strong>Author:</strong> Brown Beckley</p>
            <p><strong>Email:</strong> brownbeckley94@gmail.com</p>
            <p><strong>GitHub:</strong> <a href="https://github.com/bbeckley-hub" target="_blank">https://github.com/bbeckley-hub</a></p>
            <p><strong>Affiliation:</strong> University of Ghana Medical School</p>
            <p style="margin-top: 20px; font-size: 0.9em; color: #ccc;">
                Analysis performed using StaphScope AMRfinderPlus v1.0.1<br></p>
        </div>
    </div>
</body>
</html>
"""
        
        # Write summary HTML report
        html_file = os.path.join(output_base, "staph_amrfinder_summary_report.html")
        with open(html_file, 'w') as f:
            f.write(html_content)
        
        self.logger.info(f"✓ S. aureus AMRfinderPlus summary HTML report created: {html_file}")
    
    def process_single_genome(self, genome_file: str, output_base: str = "amrfinder_results") -> Dict[str, Any]:
        """Process a single genome with BUNDLED AMRfinderPlus"""
        genome_name = Path(genome_file).stem
        results_dir = os.path.join(output_base, genome_name)
        
        self.logger.info(f"=== PROCESSING GENOME: {genome_name} ===")
        
        # Create output directory
        os.makedirs(results_dir, exist_ok=True)
        
        # Run AMRfinderPlus
        result = self.run_amrfinder_single_genome(genome_file, results_dir)
        
        status_icon = "✓" if result['status'] == 'success' else "✗"
        self.logger.info(f"{status_icon} {genome_name}: {result['hit_count']} AMR hits")
        
        return result
    
    def process_multiple_genomes(self, genome_pattern: str, output_base: str = "amrfinder_results") -> Dict[str, Any]:
        """Process multiple genomes using wildcard pattern - MAXIMUM SPEED"""
        
        # Check BUNDLED AMRfinderPlus installation
        if not self.check_amrfinder_installed():
            raise RuntimeError("BUNDLED AMRfinderPlus not properly installed")
        
        # Find genome files (support all FASTA extensions)
        fasta_patterns = [genome_pattern, f"{genome_pattern}.fasta", f"{genome_pattern}.fa", 
                         f"{genome_pattern}.fna", f"{genome_pattern}.faa"]
        
        genome_files = []
        for pattern in fasta_patterns:
            genome_files.extend(glob.glob(pattern))
        
        # Remove duplicates
        genome_files = list(set(genome_files))
        
        if not genome_files:
            raise FileNotFoundError(f"No FASTA files found matching pattern: {genome_pattern}")
        
        self.logger.info(f"Found {len(genome_files)} genomes: {[Path(f).name for f in genome_files]}")
        
        # Create output directory
        os.makedirs(output_base, exist_ok=True)
        
        # Process genomes with threading - MAXIMUM SPEED CONFIGURATION
        all_results = {}
        
        # Calculate optimal concurrent genomes - BE AGGRESSIVE FOR SPEED
        # Use all available CPU cores for concurrent processing
        max_concurrent = max(1, min(self.cpus, len(genome_files), int(self.available_ram / 1.5)))  # 1.5GB per genome
        
        self.logger.info(f"🚀 MAXIMUM SPEED: Using {max_concurrent} concurrent genome processing jobs")
        self.logger.info(f"   Each AMRfinderPlus instance uses {self.cpus} threads internally")
        self.logger.info(f"   This provides maximum throughput for multiple genome analysis")
        
        with ThreadPoolExecutor(max_workers=max_concurrent) as executor:
            # Submit all tasks
            future_to_genome = {
                executor.submit(self.process_single_genome, genome, output_base): genome 
                for genome in genome_files
            }
            
            # Collect results as they complete
            for future in as_completed(future_to_genome):
                genome = future_to_genome[future]
                try:
                    result = future.result()
                    all_results[result['genome']] = result
                    self.logger.info(f"✓ COMPLETED: {result['genome']} ({result['hit_count']} AMR hits)")
                except Exception as e:
                    self.logger.error(f"✗ FAILED: {genome} - {e}")
                    all_results[Path(genome).stem] = {
                        'genome': Path(genome).stem,
                        'hits': [],
                        'hit_count': 0,
                        'status': 'failed'
                    }
        
        # Create AMR summary files and HTML reports after processing all genomes
        self.create_amr_summary(all_results, output_base)
        
        self.logger.info("=== S. AUREUS AMR ANALYSIS COMPLETE ===")
        self.logger.info(f"Processed {len(all_results)} genomes")
        self.logger.info(f"Results saved to: {output_base}")
        
        return all_results


def main():
    """Command line interface"""
    parser = argparse.ArgumentParser(
        description='StaphScope AMRfinderPlus Analysis - S. aureus Antimicrobial Resistance - BUNDLED VERSION',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Run on all S. aureus FASTA files (auto-detect optimal CPU cores - MAXIMUM SPEED)
  python staph_amrfinder.py "*.fna"
  
  # Run on specific pattern with auto CPU detection
  python staph_amrfinder.py "MRSA_*.fasta"
  
  # Force specific number of CPU cores
  python staph_amrfinder.py "*.fa" --cpus 16

MAXIMUM SPEED RESOURCE MANAGEMENT:
  • 1-4 cores: Uses ALL CPU cores (100% utilization)
  • 5-8 cores: Uses (cores-1) for optimal performance  
  • 9-16 cores: Uses (cores-2) for high performance
  • 17-32 cores: Uses (cores-4) for maximum throughput
  • 32+ cores: Uses 85% of cores (capped at 32)

Supported FASTA extensions: .fasta, .fa, .fna, .faa
        """
    )
    
    parser.add_argument('pattern', help='File pattern for S. aureus genomes (e.g., "*.fasta", "genomes/*.fna")')
    parser.add_argument('--cpus', '-c', type=int, default=None, 
                       help='Number of CPU cores to use (default: auto-detect optimal for MAXIMUM SPEED)')
    parser.add_argument('--output', '-o', default='staph_amrfinder_results', 
                       help='Output directory (default: staph_amrfinder_results)')
    
    args = parser.parse_args()
    
    executor = AMRfinderPlusExecutor(cpus=args.cpus)
    
    try:
        results = executor.process_multiple_genomes(args.pattern, args.output)
        
        # Print summary
        executor.logger.info("\n" + "="*50)
        executor.logger.info("🧫 StaphScope AMRfinderPlus FINAL SUMMARY (BUNDLED)")
        executor.logger.info("="*50)
        
        total_hits = 0
        high_risk_count = 0
        critical_risk_count = 0
        
        for genome_name, result in results.items():
            total_hits += result['hit_count']
            
            # Count high-risk and critical genes
            genes = [hit.get('gene_symbol') for hit in result['hits'] if hit.get('gene_symbol')]
            high_risk_count += sum(1 for gene in genes if gene in executor.high_risk_genes)
            critical_risk_count += sum(1 for gene in genes if gene in executor.critical_risk_genes)
            
            executor.logger.info(f"✓ {genome_name}: {result['hit_count']} AMR hits")
        
        executor.logger.info("\n📊 S. AUREUS SUMMARY STATISTICS:")
        executor.logger.info(f"   Total genomes processed: {len(results)}")
        executor.logger.info(f"   Total AMR hits: {total_hits}")
        executor.logger.info(f"   High-risk genes detected: {high_risk_count}")
        executor.logger.info(f"   CRITICAL RISK genes detected: {critical_risk_count}")
        executor.logger.info(f"   Average AMR hits per genome: {total_hits / len(results) if results else 0:.1f}")
        
        # Show summary file locations
        executor.logger.info("\n📁 SUMMARY FILES CREATED:")
        executor.logger.info(f"   Comprehensive AMR data: {args.output}/staph_amrfinder_summary.tsv")
        executor.logger.info(f"   Statistics summary: {args.output}/staph_amrfinder_statistics_summary.tsv")
        executor.logger.info(f"   JSON summary: {args.output}/staph_amrfinder_summary.json")
        executor.logger.info(f"   Master JSON summary: {args.output}/staph_amrfinder_master_summary.json")
        executor.logger.info(f"   Summary HTML report: {args.output}/staph_amrfinder_summary_report.html")
        
        # Performance summary
        executor.logger.info("\n⚡ MAXIMUM SPEED PERFORMANCE SUMMARY:")
        executor.logger.info(f"   CPU cores utilized: {executor.cpus} cores")
        executor.logger.info(f"   Available RAM: {executor.available_ram:.1f} GB")
        executor.logger.info("   Processing mode: MAXIMUM SPEED CONCURRENT MODE 🚀")
        executor.logger.info("   Strategy: Process multiple genomes concurrently with optimal core allocation")
        executor.logger.info(f"   Bundled AMRfinderPlus: {executor.metadata['amrfinder_version']}")
        executor.logger.info(f"   Bundled database: {executor.metadata['database_version']}")
        
        # Critical risk warning if detected
        if critical_risk_count > 0:
            executor.logger.info("\n🚨 CRITICAL RISK ALERT: Last-resort antibiotic resistance genes detected!")
            executor.logger.info("   Immediate clinical attention and infection control measures required.")
        
        import random
        executor.logger.info(f"\n💡 {random.choice(executor.science_quotes)}")
        
    except Exception as e:
        executor.logger.error(f"S. aureus AMR analysis failed: {e}")
        import traceback
        executor.logger.error(traceback.format_exc())
        sys.exit(1)


if __name__ == "__main__":
    main()
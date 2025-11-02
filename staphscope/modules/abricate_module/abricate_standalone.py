#!/usr/bin/env python3
"""
StaphScope ABRicate Standalone Module
Comprehensive ABRicate analysis with HTML and summary.tsv reporting only - MAXIMUM SPEED VERSION
Author: Beckley Brown <brownbeckley94@gmail.com>
"""

import subprocess
import sys
import os
import glob
import logging
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed
from typing import List, Dict, Any
import argparse
import re
from datetime import datetime
import psutil
import math

class AbricateExecutor:
    """ABRicate executor with comprehensive HTML reporting - MAXIMUM SPEED"""
    
    def __init__(self, cpus: int = None):
        # Setup logging FIRST
        self.logger = self._setup_logging()
        
        # Initialize available_ram before calculating cpus
        self.available_ram = self._get_available_ram()
        
        # Then calculate resources - MAXIMUM SPEED MODE
        self.cpus = self._calculate_optimal_cpus(cpus)
        
        # All databases to be used - including the ones you specified
        self.required_databases = [
            'ncbi', 'card', 'resfinder', 'vfdb', 'argannot', 
            'plasmidfinder', 'megares', 'ecoh', 'ecoli_vf'
        ]
        
        # S. aureus specific genes for filtering
        self.sa_genes = ['mecA', 'mecC', 'lukS-PV', 'lukF-PV', 'tst', 'eta', 'etb', 'etd', 
                        'sea', 'seb', 'sec', 'sed', 'see', 'hlg', 'hla', 'hlb', 'hld',
                        'clf', 'fnb', 'spa', 'sdr', 'isd', 'cap', 'geh', 'aur', 'ssp', 
                        'map', 'vWbp', 'sak', 'scn', 'chp', 'adsA', 'esx', 'ess', 'esa']
    
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

    def check_abricate_installed(self) -> bool:
        """Check if ABRicate is installed and meets version requirements"""
        try:
            result = subprocess.run(['abricate', '--version'], 
                                  capture_output=True, text=True, check=True)
            version_line = result.stdout.strip()
            self.logger.info("ABRicate version: %s", version_line)
            
            version_match = re.search(r'(\d+\.\d+\.\d+)', version_line)
            if version_match:
                version_str = version_match.group(1)
                if version_str >= "1.0.1":
                    self.logger.info("✓ ABRicate version meets requirement (>=1.0.1)")
                    return True
                else:
                    self.logger.error("ABRicate version too old: %s. Required >=1.0.1", version_str)
                    return False
            self.logger.info("✓ ABRicate installed (version check skipped)")
            return True
            
        except (subprocess.CalledProcessError, FileNotFoundError):
            self.logger.error("ABRicate not found. Please install with: conda install -c bioconda abricate")
            return False
    
    def setup_abricate_databases(self):
        """Setup ABRicate databases if they don't exist"""
        self.logger.info("Setting up ABRicate databases...")
        
        available_dbs = []
        missing_dbs = []
        
        try:
            # Check which databases exist
            check_result = subprocess.run(['abricate', '--list'], 
                                        capture_output=True, text=True, check=True)
            
            for db in self.required_databases:
                if db in check_result.stdout:
                    self.logger.info("✓ Database available: %s", db)
                    available_dbs.append(db)
                else:
                    self.logger.warning("Database not available: %s", db)
                    missing_dbs.append(db)
            
            # Setup missing databases
            for db in missing_dbs:
                self.logger.info("Attempting to setup database: %s", db)
                try:
                    result = subprocess.run(
                        ['abricate', '--setupdb', '--db', db],
                        capture_output=True, text=True, check=True
                    )
                    self.logger.info("✓ Database setup completed: %s", db)
                    available_dbs.append(db)
                except subprocess.CalledProcessError as e:
                    self.logger.error("Failed to setup database %s: %s", db, e.stderr)
            
            # Update required databases to only include available ones
            self.required_databases = available_dbs
            self.logger.info("Using databases: %s", ", ".join(self.required_databases))
            
        except subprocess.CalledProcessError as e:
            self.logger.error("Error checking ABRicate databases: %s", e.stderr)
        except Exception as e:
            self.logger.error("Unexpected error setting up databases: %s", e)
    
    def run_abricate_single_db(self, genome_file: str, database: str, output_dir: str) -> Dict[str, Any]:
        """Run ABRicate on a single genome with specific database"""
        genome_name = Path(genome_file).stem
        output_file = os.path.join(output_dir, f"abricate_{database}.txt")
        
        cmd = [
            'abricate', 
            genome_file, 
            '--db', database,
            '--minid', '80',
            '--mincov', '80'
        ]
        
        self.logger.info("Running ABRicate: %s --db %s", genome_name, database)
        
        try:
            with open(output_file, 'w') as outfile:
                result = subprocess.run(cmd, stdout=outfile, stderr=subprocess.PIPE, text=True, check=True)
            
            # Parse results for reporting
            hits = self._parse_abricate_output(output_file)
            
            # Create individual database HTML report
            self._create_database_html_report(genome_name, database, hits, output_dir)
            
            return {
                'database': database,
                'genome': genome_name,
                'output_file': output_file,
                'hits': hits,
                'hit_count': len(hits),
                'status': 'success'
            }
            
        except subprocess.CalledProcessError as e:
            self.logger.error("ABRicate failed for %s on %s: %s", database, genome_name, e.stderr)
            return {
                'database': database,
                'genome': genome_name,
                'output_file': output_file,
                'hits': [],
                'hit_count': 0,
                'status': 'failed'
            }
    
    def _parse_abricate_output(self, abricate_file: str) -> List[Dict]:
        """Parse ABRicate output file - ROBUST VERSION that handles tabs in fields"""
        hits = []
        try:
            with open(abricate_file, 'r') as f:
                lines = f.readlines()
                
            if not lines:
                return hits
                
            # Find header line
            headers = []
            data_lines = []
            
            for line in lines:
                if line.startswith('#FILE') and not headers:
                    # This is the header line
                    headers = line.strip().replace('#', '').split('\t')
                elif line.strip() and not line.startswith('#'):
                    data_lines.append(line.strip())
            
            if not headers:
                self.logger.warning("No headers found in %s", abricate_file)
                return hits
            
            # Expected column count based on standard ABRicate output
            expected_columns = len(headers)
                
            # Parse data lines with robust tab handling
            for line_num, line in enumerate(data_lines, 1):
                # Split by tab but be careful about fields that contain tabs
                parts = line.split('\t')
                
                # Handle cases where there are more parts than headers due to tabs in fields
                if len(parts) > expected_columns:
                    # Combine extra fields into the last column (usually PRODUCT)
                    combined_parts = parts[:expected_columns-1]  # Take all but the last expected column
                    combined_parts.append('\t'.join(parts[expected_columns-1:]))  # Combine the rest into PRODUCT
                    parts = combined_parts
                elif len(parts) < expected_columns:
                    # Pad with empty strings if fewer columns
                    parts.extend([''] * (expected_columns - len(parts)))
                
                if len(parts) == expected_columns:
                    hit = {}
                    for i, header in enumerate(headers):
                        hit[header] = parts[i] if i < len(parts) else ''
                    
                    # Map to consistent field names
                    processed_hit = {
                        'file': hit.get('FILE', ''),
                        'sequence': hit.get('SEQUENCE', ''),
                        'start': hit.get('START', ''),
                        'end': hit.get('END', ''),
                        'strand': hit.get('STRAND', ''),
                        'gene': hit.get('GENE', ''),
                        'coverage': hit.get('COVERAGE', ''),
                        'coverage_map': hit.get('COVERAGE_MAP', ''),
                        'gaps': hit.get('GAPS', ''),
                        'coverage_percent': hit.get('%COVERAGE', ''),
                        'identity_percent': hit.get('%IDENTITY', ''),
                        'database': hit.get('DATABASE', ''),
                        'accession': hit.get('ACCESSION', ''),
                        'product': hit.get('PRODUCT', ''),
                        'resistance': hit.get('RESISTANCE', '')
                    }
                    hits.append(processed_hit)
                else:
                    self.logger.warning("Line %d has %d parts, expected %d: %s", 
                                      line_num, len(parts), expected_columns, line[:100] + "...")
                    
        except Exception as e:
            self.logger.error("Error parsing %s: %s", abricate_file, e)
            
        self.logger.info("Parsed %d hits from %s", len(hits), abricate_file)
        return hits
    
    def _create_database_html_report(self, genome_name: str, database: str, hits: List[Dict], output_dir: str):
        """Create individual HTML report for each database"""
        html_content = f"""
<!DOCTYPE html>
<html>
<head>
    <title>StaphScope ABRicate Report - {database.upper()}</title>
    <style>
        body {{ font-family: Arial, sans-serif; margin: 20px; }}
        .header {{ background: #2c3e50; color: white; padding: 20px; border-radius: 5px; }}
        .summary {{ background: #ecf0f1; padding: 15px; border-radius: 5px; margin: 10px 0; }}
        .gene-table {{ width: 100%; border-collapse: collapse; margin: 10px 0; }}
        .gene-table th, .gene-table td {{ border: 1px solid #ddd; padding: 8px; text-align: left; }}
        .gene-table th {{ background-color: #34495e; color: white; }}
        .present {{ background-color: #d4edda; }}
    </style>
</head>
<body>
    <div class="header">
        <h1>StaphScope ABRicate Analysis Report</h1>
        <p>Genome: {genome_name} | Database: {database.upper()} | Generated on: {datetime.now().strftime('%Y-%m-%d %H:%M')}</p>
    </div>
    
    <div class="summary">
        <h2>Database Summary</h2>
        <p><strong>Total genes detected:</strong> {len(hits)}</p>
        <p><strong>Database:</strong> {database.upper()}</p>
    </div>
"""
        
        if hits:
            html_content += """
    <h2>Genes Detected</h2>
    <table class="gene-table">
        <tr>
            <th>Gene</th>
            <th>Product</th>
            <th>Coverage</th>
            <th>Identity</th>
            <th>Accession</th>
        </tr>
"""
            
            for hit in hits:
                # Truncate very long product descriptions for display
                product_display = hit['product']
                if len(product_display) > 150:
                    product_display = product_display[:147] + "..."
                
                html_content += f"""
        <tr class="present">
            <td>{hit['gene']}</td>
            <td title="{hit['product']}">{product_display}</td>
            <td>{hit['coverage_percent']}%</td>
            <td>{hit['identity_percent']}%</td>
            <td>{hit['accession']}</td>
        </tr>
"""
            
            html_content += "</table>"
        else:
            html_content += """
    <div class="summary">
        <h2>No Genes Detected</h2>
        <p>No significant hits found in this database.</p>
    </div>
"""
        
        html_content += """
</body>
</html>
"""
        
        # Write individual database HTML report
        html_file = os.path.join(output_dir, f"abricate_{database}_report.html")
        with open(html_file, 'w') as f:
            f.write(html_content)
        
        self.logger.info("Individual database report: %s", html_file)
    
    def create_database_summaries(self, all_results: Dict[str, Any], output_base: str):
        """Create ABRicate summary files for each database across all genomes"""
        self.logger.info("Creating database summary files...")
        
        # Group results by database
        db_results = {}
        for genome_name, genome_result in all_results.items():
            for db, db_result in genome_result['results'].items():
                if db not in db_results:
                    db_results[db] = []
                
                # Add hits with genome name
                for hit in db_result['hits']:
                    hit['genome'] = genome_name
                    db_results[db].append(hit)
        
        # Create summary file for each database
        for db, hits in db_results.items():
            if hits:
                summary_file = os.path.join(output_base, f"{db}_abricate_summary.tsv")
                
                # Get headers from first hit
                headers = list(hits[0].keys())
                
                with open(summary_file, 'w') as f:
                    # Write header
                    f.write('\t'.join(headers) + '\n')
                    
                    # Write data
                    for hit in hits:
                        row = [str(hit.get(header, '')) for header in headers]
                        f.write('\t'.join(row) + '\n')
                
                self.logger.info("✓ Created %s summary: %s (%d hits)", db, summary_file, len(hits))
            else:
                self.logger.info("No hits for database %s, skipping summary", db)
    
    def analyze_saureus_genes(self, all_hits: List[Dict]) -> Dict[str, Any]:
        """Analyze S. aureus specific genes with improved detection and partial PVL handling"""
        analysis = {
            'mrsa_status': 'negative',
            'pvl_status': 'negative',  # negative, partial, or positive
            'pvl_genes_found': [],
            'sa_genes_found': [],
            'total_sa_genes': 0,
            'total_hits': len(all_hits)
        }
        
        # Expanded S. aureus gene patterns
        sa_gene_patterns = [
            'mecA', 'mecC', 'lukS-PV', 'lukF-PV', 'spa', 'sdr', 'clf', 'fnb',
            'hly', 'hla', 'hlb', 'hld', 'hlg', 'sel', 'se', 'et', 'tst',
            'isd', 'cap', 'geh', 'aur', 'ssp', 'map', 'vWbp', 'sak', 'scn', 'chp',
            'adsA', 'esx', 'ess', 'esa'
        ]
        
        pvl_genes_found = []
        
        for hit in all_hits:
            gene = hit['gene']
            
            # Check if this is an S. aureus gene
            is_sa_gene = any(pattern.lower() in gene.lower() for pattern in sa_gene_patterns)
            
            if is_sa_gene:
                # Determine function
                function = 'Resistance'
                if any(res in gene.lower() for res in ['mec']):
                    function = 'Resistance'
                    analysis['mrsa_status'] = 'positive'
                else:
                    function = 'Virulence'
                
                # Track PVL genes
                if gene in ['lukS-PV', 'lukF-PV']:
                    pvl_genes_found.append(gene)
                    analysis['pvl_genes_found'].append(gene)
                
                analysis['sa_genes_found'].append({
                    'gene': gene,
                    'product': hit['product'],  # Now correctly gets the full product description
                    'database': hit['database'],
                    'coverage': hit['coverage_percent'],  # Correct %coverage
                    'identity': hit['identity_percent'],  # Correct %identity
                    'function': function
                })
        
        # Set PVL status with partial detection
        unique_pvl_genes = set(pvl_genes_found)
        if len(unique_pvl_genes) == 2:  # Both lukS-PV and lukF-PV present
            analysis['pvl_status'] = 'positive'
        elif len(unique_pvl_genes) == 1:  # Only one PVL gene present
            analysis['pvl_status'] = 'partial'
        else:  # No PVL genes
            analysis['pvl_status'] = 'negative'
            
        analysis['total_sa_genes'] = len(analysis['sa_genes_found'])
        return analysis
    
    def create_comprehensive_html_report(self, genome_name: str, results: Dict, output_dir: str):
        """Create comprehensive HTML report with correct data mapping"""
        
        # Collect all hits
        all_hits = []
        for db_result in results.values():
            all_hits.extend(db_result['hits'])
        
        # Analyze S. aureus genes
        analysis = self.analyze_saureus_genes(all_hits)
        
        # Create HTML content
        html_content = f"""
<!DOCTYPE html>
<html>
<head>
    <title>StaphScope ABRicate Analysis Report</title>
    <style>
        body {{ font-family: Arial, sans-serif; margin: 20px; }}
        .header {{ background: #2c3e50; color: white; padding: 20px; border-radius: 5px; }}
        .summary {{ background: #ecf0f1; padding: 15px; border-radius: 5px; margin: 10px 0; }}
        .gene-table {{ width: 100%; border-collapse: collapse; margin: 10px 0; }}
        .gene-table th, .gene-table td {{ border: 1px solid #ddd; padding: 8px; text-align: left; }}
        .gene-table th {{ background-color: #34495e; color: white; }}
        .present {{ background-color: #d4edda; }}
        .absent {{ background-color: #f8d7da; }}
        .mrsa {{ background-color: #fff3cd; font-weight: bold; }}
        .pvl-partial {{ background-color: #ffeaa7; font-weight: bold; }}
    </style>
</head>
<body>
    <div class="header">
        <h1>StaphScope ABRicate Analysis Report</h1>
        <p>Genome: {genome_name} | Generated on: {datetime.now().strftime('%Y-%m-%d %H:%M')}</p>
    </div>
    
    <div class="summary">
        <h2>Summary</h2>
        <p><strong>Total genes detected:</strong> {analysis['total_hits']}</p>
        <p><strong>S. aureus virulence/resistance genes:</strong> {analysis['total_sa_genes']}</p>
        <p><strong>Databases analyzed:</strong> {len(results)}</p>
    </div>
"""
        
        # MRSA status
        if analysis['mrsa_status'] == 'positive':
            html_content += f"""
    <div class="summary mrsa">
        <h2>🔴 MRSA DETECTED</h2>
        <p>mecA/mecC genes found in database results</p>
    </div>
"""
        else:
            html_content += """
    <div class="summary">
        <h2>🟢 MSSA (No MRSA markers detected)</h2>
        <p>No mecA or mecC genes found</p>
    </div>
"""
        
        # PVL status - WITH PARTIAL DETECTION
        if analysis['pvl_status'] == 'positive':
            html_content += f"""
    <div class="summary">
        <h2>🦠 PVL POSITIVE</h2>
        <p>Both Panton-Valentine Leukocidin genes detected: {', '.join(analysis['pvl_genes_found'])}</p>
    </div>
"""
        elif analysis['pvl_status'] == 'partial':
            html_content += f"""
    <div class="summary pvl-partial">
        <h2>🟡 PARTIAL PVL DETECTED</h2>
        <p>Only one PVL gene detected: {', '.join(analysis['pvl_genes_found'])}</p>
        <p><em>Note: Complete PVL requires both lukS-PV and lukF-PV genes</em></p>
    </div>
"""
        else:
            html_content += """
    <div class="summary">
        <h2>🔵 PVL NEGATIVE</h2>
        <p>No PVL genes detected</p>
    </div>
"""
        
        # S. aureus genes table - WITH CORRECT DATA MAPPING
        if analysis['sa_genes_found']:
            html_content += """
    <h2>S. aureus Specific Genes Detected</h2>
    <table class="gene-table">
        <tr>
            <th>Gene</th>
            <th>Product</th>
            <th>Database</th>
            <th>Coverage</th>
            <th>Identity</th>
            <th>Function</th>
        </tr>
"""
            
            for gene_info in analysis['sa_genes_found']:
                # Truncate long product descriptions for display
                product_display = gene_info['product']
                if len(product_display) > 100:
                    product_display = product_display[:97] + "..."
                
                html_content += f"""
        <tr class="present">
            <td>{gene_info['gene']}</td>
            <td title="{gene_info['product']}">{product_display}</td>
            <td>{gene_info['database']}</td>
            <td>{gene_info['coverage']}%</td>
            <td>{gene_info['identity']}%</td>
            <td>{gene_info['function']}</td>
        </tr>
"""
            
            html_content += "</table>"
        
        # Database summary table
        html_content += """
    <h2>Database Results Summary</h2>
    <table class="gene-table">
        <tr>
            <th>Database</th>
            <th>Hits</th>
            <th>Status</th>
        </tr>
"""
        
        for db, result in results.items():
            status_icon = "✅" if result['status'] == 'success' else "❌"
            html_content += f"""
        <tr>
            <td>{db}</td>
            <td>{result['hit_count']}</td>
            <td>{status_icon} {result['status']}</td>
        </tr>
"""
        
        html_content += """
    </table>
</body>
</html>
"""
        
        # Write comprehensive HTML report
        html_file = os.path.join(output_dir, f"{genome_name}_comprehensive_abricate_report.html")
        with open(html_file, 'w') as f:
            f.write(html_content)
        
        self.logger.info("Comprehensive HTML report generated: %s", html_file)
    
    def process_single_genome(self, genome_file: str, output_base: str = "abricate_results") -> Dict[str, Any]:
        """Process a single genome with all databases and HTML reporting"""
        genome_name = Path(genome_file).stem
        results_dir = os.path.join(output_base, genome_name)
        
        self.logger.info("=== PROCESSING GENOME: %s ===", genome_name)
        
        # Create output directory
        os.makedirs(results_dir, exist_ok=True)
        
        databases = self.required_databases
        
        # Run ABRicate on all databases
        results = {}
        for db in databases:
            result = self.run_abricate_single_db(genome_file, db, results_dir)
            results[db] = result
            status_icon = "✓" if result['status'] == 'success' else "✗"
            self.logger.info("%s %s: %d hits", status_icon, db, result['hit_count'])
        
        # Create comprehensive HTML report
        self.create_comprehensive_html_report(genome_name, results, results_dir)
        
        return {
            'genome': genome_name,
            'results': results,
            'total_hits': sum(r['hit_count'] for r in results.values())
        }
    
    def process_multiple_genomes(self, genome_pattern: str, output_base: str = "abricate_results") -> Dict[str, Any]:
        """Process multiple genomes using wildcard pattern - MAXIMUM SPEED"""
        
        # Check ABRicate installation
        if not self.check_abricate_installed():
            raise RuntimeError("ABRicate not properly installed")
        
        # Setup databases
        self.setup_abricate_databases()
        
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
        
        self.logger.info("Found %d genomes: %s", len(genome_files), [Path(f).name for f in genome_files])
        
        # Create output directory
        os.makedirs(output_base, exist_ok=True)
        
        # Process genomes with parallel execution - MAXIMUM SPEED CONFIGURATION
        all_results = {}
        
        if len(genome_files) > 1 and self.cpus > 1:
            # Use ThreadPoolExecutor for parallel processing of multiple genomes
            self.logger.info("Using parallel processing with %d CPU cores (MAXIMUM SPEED)", self.cpus)
            
            with ThreadPoolExecutor(max_workers=self.cpus) as executor:
                # Submit all genomes for processing
                future_to_genome = {
                    executor.submit(self.process_single_genome, genome, output_base): genome 
                    for genome in genome_files
                }
                
                # Collect results as they complete
                for future in as_completed(future_to_genome):
                    genome = future_to_genome[future]
                    try:
                        result = future.result()
                        all_results[Path(genome).stem] = result
                        self.logger.info("✓ Completed: %s (%d total hits)", result['genome'], result['total_hits'])
                    except Exception as e:
                        self.logger.error("✗ Failed: %s - %s", genome, e)
        else:
            # Process genomes sequentially
            for genome in genome_files:
                try:
                    result = self.process_single_genome(genome, output_base)
                    all_results[Path(genome).stem] = result
                    self.logger.info("✓ Completed: %s (%d total hits)", result['genome'], result['total_hits'])
                except Exception as e:
                    self.logger.error("✗ Failed: %s - %s", genome, e)
        
        # Create database summary files after processing all genomes
        self.create_database_summaries(all_results, output_base)
        
        self.logger.info("=== ANALYSIS COMPLETE ===")
        self.logger.info("Processed %d genomes", len(all_results))
        self.logger.info("Results saved to: %s", output_base)
        
        return all_results


def main():
    """Command line interface"""
    parser = argparse.ArgumentParser(
        description='StaphScope ABRicate Analysis - MAXIMUM SPEED VERSION',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Run on all FASTA files (auto-detect optimal CPU cores - MAXIMUM SPEED)
  python abricate_standalone.py "*.fasta"
  
  # Run on specific pattern with auto CPU detection
  python abricate_standalone.py "MRSA_*.fna"
  
  # Run with custom output directory and auto CPUs
  python abricate_standalone.py "*.fa" --output my_results

  # Force specific number of CPU cores
  python abricate_standalone.py "*.fna" --cpus 4

MAXIMUM SPEED RESOURCE MANAGEMENT:
  • 1-4 cores: Uses ALL CPU cores (100% utilization)
  • 5-8 cores: Uses (cores-1) for optimal performance  
  • 9-16 cores: Uses (cores-2) for high performance
  • 17-32 cores: Uses (cores-4) for maximum throughput
  • 32+ cores: Uses 85% of cores (capped at 32)

Supported FASTA extensions: .fasta, .fa, .fna, .faa
        """
    )
    
    parser.add_argument('pattern', help='File pattern for genomes (e.g., "*.fasta", "genomes/*.fna")')
    parser.add_argument('--cpus', '-c', type=int, default=None, 
                       help='Number of CPU cores to use (default: auto-detect optimal for MAXIMUM SPEED)')
    parser.add_argument('--output', '-o', default='abricate_results', 
                       help='Output directory (default: abricate_results)')
    
    args = parser.parse_args()
    
    executor = AbricateExecutor(cpus=args.cpus)
    
    try:
        results = executor.process_multiple_genomes(args.pattern, args.output)
        
        # Print summary
        executor.logger.info("\n" + "="*50)
        executor.logger.info("📊 FINAL SUMMARY")
        executor.logger.info("="*50)
        for genome_name, result in results.items():
            executor.logger.info("✓ %s: %d total hits across %d databases", 
                               genome_name, result['total_hits'], len(result['results']))
        
        # Database usage summary
        executor.logger.info("\n" + "="*50)
        executor.logger.info("🗃️  DATABASE USAGE SUMMARY")
        executor.logger.info("="*50)
        executor.logger.info("Used databases: %s", ", ".join(executor.required_databases))
        
        # Performance summary
        executor.logger.info("\n" + "="*50)
        executor.logger.info("⚡ MAXIMUM SPEED PERFORMANCE SUMMARY")
        executor.logger.info("="*50)
        executor.logger.info("CPU cores utilized: %d cores", executor.cpus)
        executor.logger.info("Available RAM: %.1f GB", executor.available_ram)
        executor.logger.info("Total genomes processed: %d", len(results))
        executor.logger.info("Processing mode: MAXIMUM SPEED 🚀")
        
    except Exception as e:
        executor.logger.error("Analysis failed: %s", e)
        sys.exit(1)


if __name__ == "__main__":
    main()

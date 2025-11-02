#!/usr/bin/env python3
"""
StaphScope AMRfinderPlus Standalone Module
Comprehensive AMR analysis with HTML reporting
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

class AMRfinderPlusExecutor:
    """AMRfinderPlus executor with comprehensive HTML reporting"""
    
    def __init__(self, threads: int = 4):
        self.threads = threads
        self.logger = self._setup_logging()
        
    def _setup_logging(self):
        logging.basicConfig(
            level=logging.INFO,
            format='%(asctime)s - %(levelname)s - %(message)s'
        )
        return logging.getLogger(__name__)
    
    def check_amrfinder_installed(self) -> bool:
        """Check if AMRfinderPlus is installed and meets requirements"""
        try:
            result = subprocess.run(['amrfinder', '--version'], 
                                  capture_output=True, text=True, check=True)
            version_line = result.stdout.strip()
            self.logger.info("AMRfinderPlus version: %s", version_line)
            
            # Check for version number in output
            version_match = re.search(r'(\d+\.\d+\.\d+)', version_line)
            if version_match:
                version_str = version_match.group(1)
                self.logger.info("✓ AMRfinderPlus version: %s", version_str)
            else:
                self.logger.info("✓ AMRfinderPlus installed (version check skipped)")
            
            return True
            
        except (subprocess.CalledProcessError, FileNotFoundError):
            self.logger.error("AMRfinderPlus not found. Installing with conda...")
            return self._install_amrfinder()
    
    def _install_amrfinder(self) -> bool:
        """Install AMRfinderPlus using conda"""
        try:
            self.logger.info("Installing AMRfinderPlus with conda...")
            # Run conda install silently
            result = subprocess.run([
                'conda', 'install', '-c', 'bioconda', 'ncbi-amrfinderplus', '-y'
            ], capture_output=True, text=True, check=True)
            
            self.logger.info("✓ AMRfinderPlus installed successfully")
            return True
            
        except subprocess.CalledProcessError as e:
            self.logger.error("Failed to install AMRfinderPlus with conda. Trying pip...")
            try:
                # Try pip as fallback
                result = subprocess.run([
                    'pip', 'install', 'ncbi-amrfinderplus'
                ], capture_output=True, text=True, check=True)
                self.logger.info("✓ AMRfinderPlus installed successfully with pip")
                return True
            except subprocess.CalledProcessError as e2:
                self.logger.error("Failed to install AMRfinderPlus: %s", e2.stderr)
                return False
    
    def setup_amrfinder_database(self):
        """Setup and update AMRfinderPlus database"""
        try:
            self.logger.info("Updating AMRfinderPlus database...")
            # Run with --force to ensure database is setup
            result = subprocess.run([
                'amrfinder', '--update', '--force'
            ], capture_output=True, text=True, check=True)
            
            self.logger.info("✓ AMRfinderPlus database updated successfully")
            return True
            
        except subprocess.CalledProcessError as e:
            self.logger.error("Failed to update AMRfinderPlus database: %s", e.stderr)
            
            # Try alternative method
            try:
                self.logger.info("Trying alternative database setup method...")
                result = subprocess.run([
                    'amrfinder', '-u'
                ], capture_output=True, text=True, check=True)
                self.logger.info("✓ AMRfinderPlus database setup completed")
                return True
            except subprocess.CalledProcessError as e2:
                self.logger.error("Database setup failed: %s", e2.stderr)
                return False
    
    def run_amrfinder_single_genome(self, genome_file: str, output_dir: str) -> Dict[str, Any]:
        """Run AMRfinderPlus on a single genome"""
        genome_name = Path(genome_file).stem
        output_file = os.path.join(output_dir, f"{genome_name}_amrfinder.txt")
        
        cmd = [
            'amrfinder',
            '--nucleotide', genome_file,
            '--output', output_file,
            '--threads', str(self.threads),
            '--plus'
        ]
        
        self.logger.info("Running AMRfinderPlus: %s", genome_name)
        
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
            self.logger.error("AMRfinderPlus failed for %s: %s", genome_name, e.stderr)
            return {
                'genome': genome_name,
                'output_file': output_file,
                'hits': [],
                'hit_count': 0,
                'status': 'failed'
            }
    
    def _parse_amrfinder_output(self, amrfinder_file: str) -> List[Dict]:
        """Parse AMRfinderPlus output file into structured data"""
        hits = []
        try:
            with open(amrfinder_file, 'r') as f:
                lines = f.readlines()
                
            if not lines or len(lines) < 2:  # Need at least header and one data line
                return hits
                
            # Parse header
            headers = lines[0].strip().split('\t')
            
            # Parse data lines
            for line_num, line in enumerate(lines[1:], 2):
                line = line.strip()
                if not line:
                    continue
                    
                parts = line.split('\t')
                if len(parts) >= len(headers):
                    hit = {}
                    for i, header in enumerate(headers):
                        if i < len(parts):
                            hit[header] = parts[i]
                        else:
                            hit[header] = ''
                    
                    # Map to consistent field names
                    processed_hit = {
                        'protein_id': hit.get('Protein identifier', ''),
                        'contig_id': hit.get('Contig id', ''),
                        'start': hit.get('Start', ''),
                        'stop': hit.get('Stop', ''),
                        'strand': hit.get('Strand', ''),
                        'gene_symbol': hit.get('Gene symbol', ''),
                        'sequence_name': hit.get('Sequence name', ''),
                        'scope': hit.get('Scope', ''),
                        'element_type': hit.get('Element type', ''),
                        'element_subtype': hit.get('Element subtype', ''),
                        'class': hit.get('Class', ''),
                        'subclass': hit.get('Subclass', ''),
                        'method': hit.get('Method', ''),
                        'target_length': hit.get('Target length', ''),
                        'ref_length': hit.get('Reference sequence length', ''),
                        'coverage': hit.get('% Coverage of reference sequence', ''),
                        'identity': hit.get('% Identity to reference sequence', ''),
                        'alignment_length': hit.get('Alignment length', ''),
                        'accession': hit.get('Accession of closest sequence', ''),
                        'closest_name': hit.get('Name of closest sequence', ''),
                        'hmm_id': hit.get('HMM id', ''),
                        'hmm_description': hit.get('HMM description', '')
                    }
                    hits.append(processed_hit)
                else:
                    self.logger.warning("Line %d has %d parts, expected %d: %s", 
                                      line_num, len(parts), len(headers), line[:100] + "...")
                    
        except Exception as e:
            self.logger.error("Error parsing %s: %s", amrfinder_file, e)
            
        self.logger.info("Parsed %d AMR hits from %s", len(hits), amrfinder_file)
        return hits
    
    def _create_amrfinder_html_report(self, genome_name: str, hits: List[Dict], output_dir: str):
        """Create comprehensive HTML report for AMRfinderPlus results"""
        
        # Analyze AMR results
        analysis = self._analyze_amr_results(hits)
        
        html_content = f"""
<!DOCTYPE html>
<html>
<head>
    <title>StaphScope AMRfinderPlus Analysis Report</title>
    <style>
        body {{ font-family: Arial, sans-serif; margin: 20px; }}
        .header {{ background: #2c3e50; color: white; padding: 20px; border-radius: 5px; }}
        .summary {{ background: #ecf0f1; padding: 15px; border-radius: 5px; margin: 10px 0; }}
        .gene-table {{ width: 100%; border-collapse: collapse; margin: 10px 0; }}
        .gene-table th, .gene-table td {{ border: 1px solid #ddd; padding: 8px; text-align: left; }}
        .gene-table th {{ background-color: #34495e; color: white; }}
        .present {{ background-color: #d4edda; }}
        .mrsa {{ background-color: #fff3cd; font-weight: bold; }}
        .resistance-class {{ background-color: #e8f4fd; }}
        .stat-table {{ width: auto; border-collapse: collapse; margin: 10px 0; }}
        .stat-table th, .stat-table td {{ border: 1px solid #ddd; padding: 8px; text-align: left; }}
        .stat-table th {{ background-color: #34495e; color: white; }}
    </style>
</head>
<body>
    <div class="header">
        <h1>StaphScope AMRfinderPlus Analysis Report</h1>
        <p>Genome: {genome_name} | Generated on: {datetime.now().strftime('%Y-%m-%d %H:%M')}</p>
    </div>
    
    <div class="summary">
        <h2>Summary</h2>
        <p><strong>Total AMR genes detected:</strong> {analysis['total_genes']}</p>
        <p><strong>Resistance classes:</strong> {analysis['total_classes']}</p>
        <p><strong>MRSA status:</strong> {analysis['mrsa_status']}</p>
    </div>
"""
        
        # MRSA status
        if analysis['mrsa_status'] == 'positive':
            html_content += f"""
    <div class="summary mrsa">
        <h2>🔴 MRSA DETECTED</h2>
        <p>mecA gene found in AMR results</p>
    </div>
"""
        else:
            html_content += """
    <div class="summary">
        <h2>🟢 MSSA (No MRSA markers detected)</h2>
        <p>No mecA gene found</p>
    </div>
"""
        
        # Resistance classes summary
        if analysis['resistance_classes']:
            html_content += """
    <h2>Resistance Classes Detected</h2>
    <table class="stat-table">
        <tr>
            <th>Resistance Class</th>
            <th>Genes Count</th>
            <th>Genes</th>
        </tr>
"""
            
            for class_name, genes in analysis['resistance_classes'].items():
                gene_list = ", ".join(genes)
                html_content += f"""
        <tr class="resistance-class">
            <td>{class_name}</td>
            <td>{len(genes)}</td>
            <td>{gene_list}</td>
        </tr>
"""
            
            html_content += "</table>"
        
        # Detailed AMR genes table
        if hits:
            html_content += """
    <h2>Detailed AMR Genes Detected</h2>
    <table class="gene-table">
        <tr>
            <th>Gene Symbol</th>
            <th>Sequence Name</th>
            <th>Class</th>
            <th>Subclass</th>
            <th>Coverage</th>
            <th>Identity</th>
            <th>Scope</th>
        </tr>
"""
            
            for hit in hits:
                # Determine row class based on gene type
                row_class = "present"
                if hit['gene_symbol'] == 'mecA':
                    row_class = "mrsa"
                
                html_content += f"""
        <tr class="{row_class}">
            <td>{hit['gene_symbol']}</td>
            <td title="{hit['sequence_name']}">{hit['sequence_name'][:80]}{'...' if len(hit['sequence_name']) > 80 else ''}</td>
            <td>{hit['class']}</td>
            <td>{hit['subclass']}</td>
            <td>{hit['coverage']}%</td>
            <td>{hit['identity']}%</td>
            <td>{hit['scope']}</td>
        </tr>
"""
            
            html_content += "</table>"
            
            # Full data table (collapsible)
            html_content += """
    <h2>
        <button onclick="toggleFullTable()" style="background: #34495e; color: white; border: none; padding: 10px; border-radius: 5px; cursor: pointer;">
            📊 Show Full AMR Data Table
        </button>
    </h2>
    <div id="fullTable" style="display: none;">
        <table class="gene-table">
            <tr>
                <th>Protein ID</th>
                <th>Contig</th>
                <th>Start</th>
                <th>Stop</th>
                <th>Strand</th>
                <th>Gene Symbol</th>
                <th>Sequence Name</th>
                <th>Scope</th>
                <th>Element Type</th>
                <th>Class</th>
                <th>Subclass</th>
                <th>Coverage</th>
                <th>Identity</th>
                <th>Accession</th>
            </tr>
"""
            
            for hit in hits:
                html_content += f"""
            <tr class="present">
                <td>{hit['protein_id']}</td>
                <td>{hit['contig_id']}</td>
                <td>{hit['start']}</td>
                <td>{hit['stop']}</td>
                <td>{hit['strand']}</td>
                <td>{hit['gene_symbol']}</td>
                <td title="{hit['sequence_name']}">{hit['sequence_name'][:60]}{'...' if len(hit['sequence_name']) > 60 else ''}</td>
                <td>{hit['scope']}</td>
                <td>{hit['element_type']}</td>
                <td>{hit['class']}</td>
                <td>{hit['subclass']}</td>
                <td>{hit['coverage']}%</td>
                <td>{hit['identity']}%</td>
                <td>{hit['accession']}</td>
            </tr>
"""
            
            html_content += """
        </table>
    </div>
"""
        else:
            html_content += """
    <div class="summary">
        <h2>No AMR Genes Detected</h2>
        <p>No antimicrobial resistance genes found in this genome.</p>
    </div>
"""
        
        # JavaScript for toggling full table
        html_content += """
    <script>
        function toggleFullTable() {
            var fullTable = document.getElementById("fullTable");
            var button = document.querySelector("button");
            if (fullTable.style.display === "none") {
                fullTable.style.display = "block";
                button.innerHTML = "📊 Hide Full AMR Data Table";
            } else {
                fullTable.style.display = "none";
                button.innerHTML = "📊 Show Full AMR Data Table";
            }
        }
    </script>
</body>
</html>
"""
        
        # Write HTML report
        html_file = os.path.join(output_dir, f"{genome_name}_amrfinder_report.html")
        with open(html_file, 'w') as f:
            f.write(html_content)
        
        self.logger.info("AMRfinderPlus HTML report generated: %s", html_file)
    
    def _analyze_amr_results(self, hits: List[Dict]) -> Dict[str, Any]:
        """Analyze AMR results for summary reporting"""
        analysis = {
            'total_genes': len(hits),
            'mrsa_status': 'negative',
            'resistance_classes': {},
            'total_classes': 0
        }
        
        for hit in hits:
            gene_symbol = hit['gene_symbol']
            resistance_class = hit['class']
            
            # Check for MRSA
            if gene_symbol == 'mecA':
                analysis['mrsa_status'] = 'positive'
            
            # Group by resistance class
            if resistance_class:
                if resistance_class not in analysis['resistance_classes']:
                    analysis['resistance_classes'][resistance_class] = []
                if gene_symbol not in analysis['resistance_classes'][resistance_class]:
                    analysis['resistance_classes'][resistance_class].append(gene_symbol)
        
        analysis['total_classes'] = len(analysis['resistance_classes'])
        return analysis
    
    def create_amr_summary(self, all_results: Dict[str, Any], output_base: str):
        """Create comprehensive AMR summary file for all samples"""
        self.logger.info("Creating AMR summary file for all samples...")
        
        summary_file = os.path.join(output_base, "amrfinder_summary.tsv")
        
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
        
        self.logger.info("✓ AMR summary file created: %s", summary_file)
        
        # Also create a simplified summary with statistics
        stats_file = os.path.join(output_base, "amrfinder_statistics_summary.tsv")
        with open(stats_file, 'w') as f:
            f.write("Genome\tTotal_AMR_Genes\tMRSA_Status\tResistance_Classes\tGene_List\n")
            
            for genome_name, result in all_results.items():
                # Get unique genes
                genes = list(set(hit.get('gene_symbol', '') for hit in result['hits'] if hit.get('gene_symbol')))
                gene_list = ",".join(genes)
                
                # Check MRSA status
                mrsa_status = "Positive" if any(hit.get('gene_symbol') == 'mecA' for hit in result['hits']) else "Negative"
                
                # Get resistance classes
                classes = list(set(hit.get('class', '') for hit in result['hits'] if hit.get('class')))
                class_list = ",".join(classes)
                
                f.write(f"{genome_name}\t{result['hit_count']}\t{mrsa_status}\t{class_list}\t{gene_list}\n")
        
        self.logger.info("✓ AMR statistics summary created: %s", stats_file)
    
    def process_single_genome(self, genome_file: str, output_base: str = "amrfinder_results") -> Dict[str, Any]:
        """Process a single genome with AMRfinderPlus"""
        genome_name = Path(genome_file).stem
        results_dir = os.path.join(output_base, genome_name)
        
        self.logger.info("=== PROCESSING GENOME: %s ===", genome_name)
        
        # Create output directory
        os.makedirs(results_dir, exist_ok=True)
        
        # Run AMRfinderPlus
        result = self.run_amrfinder_single_genome(genome_file, results_dir)
        
        status_icon = "✓" if result['status'] == 'success' else "✗"
        self.logger.info("%s %s: %d AMR hits", status_icon, genome_name, result['hit_count'])
        
        return result
    
    def process_multiple_genomes(self, genome_pattern: str, output_base: str = "amrfinder_results") -> Dict[str, Any]:
        """Process multiple genomes using wildcard pattern"""
        
        # Check AMRfinderPlus installation
        if not self.check_amrfinder_installed():
            raise RuntimeError("AMRfinderPlus not properly installed")
        
        # Setup database
        if not self.setup_amrfinder_database():
            raise RuntimeError("AMRfinderPlus database setup failed")
        
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
        
        # Process genomes with threading
        all_results = {}
        
        with ThreadPoolExecutor(max_workers=self.threads) as executor:
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
                    self.logger.info("✓ Completed: %s (%d AMR hits)", result['genome'], result['hit_count'])
                except Exception as e:
                    self.logger.error("✗ Failed: %s - %s", genome, e)
                    all_results[Path(genome).stem] = {
                        'genome': Path(genome).stem,
                        'hits': [],
                        'hit_count': 0,
                        'status': 'failed'
                    }
        
        # Create AMR summary files after processing all genomes
        self.create_amr_summary(all_results, output_base)
        
        self.logger.info("=== AMR ANALYSIS COMPLETE ===")
        self.logger.info("Processed %d genomes", len(all_results))
        self.logger.info("Results saved to: %s", output_base)
        
        return all_results


def main():
    """Command line interface"""
    parser = argparse.ArgumentParser(
        description='StaphScope AMRfinderPlus Analysis - Comprehensive AMR Reporting',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Run on all FASTA files in current directory
  python amrfinder_standalone.py "*.fna"
  
  # Run on specific pattern with more threads
  python amrfinder_standalone.py "MRSA_*.fasta" --threads 8
  
  # Run with custom output directory
  python amrfinder_standalone.py "*.fa" --output my_amr_results

Supported FASTA extensions: .fasta, .fa, .fna, .faa
        """
    )
    
    parser.add_argument('pattern', help='File pattern for genomes (e.g., "*.fasta", "genomes/*.fna")')
    parser.add_argument('--threads', '-t', type=int, default=4, help='Number of threads (default: 4)')
    parser.add_argument('--output', '-o', default='amrfinder_results', help='Output directory (default: amrfinder_results)')
    
    args = parser.parse_args()
    
    executor = AMRfinderPlusExecutor(threads=args.threads)
    
    try:
        results = executor.process_multiple_genomes(args.pattern, args.output)
        
        # Print summary
        executor.logger.info("\n" + "="*50)
        executor.logger.info("🧬 AMRfinderPlus FINAL SUMMARY")
        executor.logger.info("="*50)
        
        total_hits = 0
        mrsa_count = 0
        
        for genome_name, result in results.items():
            total_hits += result['hit_count']
            # Check if MRSA was detected (mecA in hits)
            if any(hit.get('gene_symbol') == 'mecA' for hit in result['hits']):
                mrsa_count += 1
            
            executor.logger.info("✓ %s: %d AMR hits", genome_name, result['hit_count'])
        
        executor.logger.info("\n📊 SUMMARY STATISTICS:")
        executor.logger.info("   Total genomes processed: %d", len(results))
        executor.logger.info("   Total AMR hits: %d", total_hits)
        executor.logger.info("   MRSA detected in: %d genomes", mrsa_count)
        executor.logger.info("   Average AMR hits per genome: %.1f", total_hits / len(results) if results else 0)
        
        # Show summary file locations
        executor.logger.info("\n📁 SUMMARY FILES CREATED:")
        executor.logger.info("   Comprehensive AMR data: %s/amrfinder_summary.tsv", args.output)
        executor.logger.info("   Statistics summary: %s/amrfinder_statistics_summary.tsv", args.output)
        
    except Exception as e:
        executor.logger.error("AMR analysis failed: %s", e)
        sys.exit(1)


if __name__ == "__main__":
    main()

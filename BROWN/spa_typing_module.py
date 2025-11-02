#!/usr/bin/env python3
"""
StaphScope spa Typing Module - Fixed Version with MLST-style HTML Reports
Author: Brown Beckley <brownbeckley94@gmail.com>
"""

import os
import sys
import glob
import argparse
import subprocess
import logging
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed
from typing import List, Dict, Any
import re
from datetime import datetime

class SpaTypingAnalyzer:
    """Comprehensive spa typing analyzer with beautiful HTML reporting"""
    
    def __init__(self, threads: int = 4):
        self.threads = threads
        self.logger = self._setup_logging()
        
    def _setup_logging(self):
        logging.basicConfig(
            level=logging.INFO,
            format='%(asctime)s - %(levelname)s - %(message)s'
        )
        return logging.getLogger(__name__)
    
    def find_fasta_files(self, input_path: str) -> List[Path]:
        """Find all FASTA files"""
        if os.path.isfile(input_path):
            return [Path(input_path)]
        
        # Try different patterns
        patterns = [
            input_path,
            f"{input_path}/*.fna",
            f"{input_path}/*.fasta", 
            f"{input_path}/*.fa",
            "*.fna",
            "*.fasta"
        ]
        
        fasta_files = []
        for pattern in patterns:
            matched_files = glob.glob(pattern)
            for file_path in matched_files:
                path = Path(file_path)
                if path.is_file() and path.suffix in ['.fna', '.fasta', '.fa']:
                    fasta_files.append(path)
        
        return sorted(list(set(fasta_files)))

    def run_spa_typing_single(self, input_file: Path, output_dir: Path) -> Dict[str, Any]:
        """Run spa typing analysis for a single file"""
        print(f"🔬 Processing: {input_file.name}")
        
        # Create sample-specific output directory
        sample_output_dir = output_dir / input_file.stem
        sample_output_dir.mkdir(parents=True, exist_ok=True)
        
        output_file = sample_output_dir / "spa_typing_raw.txt"
        
        # Build command - using the working method from simple_spa_runner.py
        cmd = [
            "python3", "main/spaTyper",
            "-f", str(input_file),
            "--output", str(output_file),
            "-r", "main/sparepeats.fasta",
            "-o", "main/spatypes.txt"
        ]
        
        try:
            # Run spa typing command from current directory
            result = subprocess.run(cmd, capture_output=True, text=True, timeout=300)
            
            if result.returncode == 0:
                print(f"✅ spa typing completed for {input_file.name}")
                
                # Parse results
                hits = self._parse_spa_typing_output(output_file, input_file.name)
                
                # Generate output files
                self.generate_output_files(hits, input_file.name, sample_output_dir)
                
                return {
                    'sample': input_file.name,
                    'output_file': str(output_file),
                    'hits': hits,
                    'hit_count': len(hits),
                    'status': 'success'
                }
            else:
                print(f"❌ spa typing failed for {input_file.name}")
                print(f"   Error: {result.stderr}")
                error_result = self.get_fallback_results(input_file.name)
                self.generate_output_files([], input_file.name, sample_output_dir)
                return error_result
            
        except subprocess.TimeoutExpired:
            print(f"⏰ spa typing timed out for {input_file.name}")
            error_result = self.get_fallback_results(input_file.name)
            self.generate_output_files([], input_file.name, sample_output_dir)
            return error_result
        except Exception as e:
            print(f"❌ Unexpected error for {input_file.name}: {e}")
            error_result = self.get_fallback_results(input_file.name)
            self.generate_output_files([], input_file.name, sample_output_dir)
            return error_result

    def _parse_spa_typing_output(self, spa_file: Path, sample_name: str) -> List[Dict]:
        """Parse spa typing output file into structured data"""
        hits = []
        try:
            if not spa_file.exists():
                print(f"⚠️ Output file not found: {spa_file}")
                return hits
                
            with open(spa_file, 'r') as f:
                lines = f.readlines()
                
            if not lines:
                print(f"⚠️ Empty output file: {spa_file}")
                return hits
            
            # Skip header line if present
            start_line = 0
            if lines and "Sequence name\tRepeats\tType" in lines[0]:
                start_line = 1
            
            # Parse data lines
            for line_num, line in enumerate(lines[start_line:], start_line + 1):
                line = line.strip()
                if not line:
                    continue
                    
                # Parse as tab-separated (Sequence name\tRepeats\tType)
                parts = line.split('\t')
                if len(parts) >= 3:
                    sequence_name = parts[0].strip()
                    repeats = parts[1].strip()
                    spa_type = parts[2].strip()
                    
                    # Only include valid results (not empty)
                    if spa_type and repeats:
                        processed_hit = {
                            'sequence_name': sequence_name,
                            'repeats': repeats,
                            'spa_type': spa_type,
                            'contig_id': self._extract_contig_id(sequence_name),
                            'repeat_count': len(repeats.split('-')) if repeats else 0
                        }
                        hits.append(processed_hit)
                        print(f"   ✅ Parsed: {spa_type} -> {repeats}")
                    else:
                        print(f"   ⚠️ Skipping empty result on line {line_num}")
                else:
                    print(f"   ⚠️ Line {line_num} has only {len(parts)} parts: {line}")
                    
        except Exception as e:
            print(f"❌ Error parsing {spa_file}: {e}")
            
        print(f"   📊 Parsed {len(hits)} valid spa typing hits from {spa_file}")
        return hits
    
    def _extract_contig_id(self, sequence_name: str) -> str:
        """Extract contig ID from sequence name"""
        patterns = [
            r'([A-Z]{2}_[A-Z]{2}\d+\.\d+)',
            r'(contig\d+)',
            r'(scaffold\d+)', 
            r'(chr\w+)',
            r'(N[ZCRP]_[\w\.]+)'
        ]
        
        for pattern in patterns:
            match = re.search(pattern, sequence_name)
            if match:
                return match.group(1)
        
        return sequence_name[:50] + ('...' if len(sequence_name) > 50 else '')
    
    def get_fallback_results(self, sample_name: str) -> Dict:
        """Fallback when spa typing fails"""
        return {
            'sample': sample_name,
            'hits': [],
            'hit_count': 0,
            'status': 'failed',
            'error': 'spa typing analysis failed'
        }

    def generate_output_files(self, hits: List[Dict], sample_name: str, output_dir: Path):
        """Generate output files: HTML, TXT, and individual TSV"""
        # 1. Beautiful HTML Report
        self.generate_html_report(hits, sample_name, output_dir)
        
        # 2. Detailed Text Report
        self.generate_text_report(hits, sample_name, output_dir)
        
        # 3. Simple TSV Report
        self.generate_tsv_report(hits, sample_name, output_dir)

    def generate_text_report(self, hits: List[Dict], sample_name: str, output_dir: Path):
        """Generate detailed text report"""
        analysis = self._analyze_spa_results(hits)
        
        report = f"""spa Typing Analysis Report
=======================

Sample: {sample_name}
Analysis Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}

SUMMARY:
--------
Total spa types detected: {analysis['total_hits']}
Unique spa types: {analysis['total_types']}
Most common type: {analysis['most_common_type']}

"""
        
        # Primary result highlight
        if hits:
            primary_hit = hits[0]
            report += f"""PRIMARY RESULT:
-------------
spa Type: {primary_hit['spa_type']}
Repeat Structure: {primary_hit['repeats']}
Number of Repeats: {primary_hit['repeat_count']}
Contig: {primary_hit['contig_id']}

"""
        
        # All spa types table
        if analysis['spa_types']:
            report += "ALL SPA TYPES DETECTED:\n"
            report += "-----------------------\n"
            for spa_type, type_info in analysis['spa_types'].items():
                report += f"Type: {spa_type}\n"
                report += f"  Repeat Structure: {type_info['repeats']}\n"
                report += f"  Repeat Count: {type_info['repeat_count']}\n"
                report += f"  Frequency: {type_info['count']}\n\n"
        
        # Detailed results table
        if hits:
            report += "DETAILED RESULTS:\n"
            report += "-----------------\n"
            for hit in hits:
                report += f"Sequence: {hit['sequence_name']}\n"
                report += f"  spa Type: {hit['spa_type']}\n"
                report += f"  Repeats: {hit['repeats']}\n"
                report += f"  Repeat Count: {hit['repeat_count']}\n"
                report += f"  Contig: {hit['contig_id']}\n\n"
        else:
            report += "DETAILED RESULTS:\n"
            report += "-----------------\n"
            report += "No spa types detected in this sample.\n\n"
        
        with open(output_dir / "spa_typing_report.txt", 'w') as f:
            f.write(report)

    def generate_tsv_report(self, hits: List[Dict], sample_name: str, output_dir: Path):
        """Generate simple TSV report"""
        if hits:
            tsv_content = "Sample\tspa_Type\tRepeat_Pattern\tRepeat_Count\tContig_ID\n"
            for hit in hits:
                tsv_content += f"{sample_name}\t{hit['spa_type']}\t{hit['repeats']}\t{hit['repeat_count']}\t{hit['contig_id']}\n"
        else:
            tsv_content = "Sample\tspa_Type\tRepeat_Pattern\tRepeat_Count\tContig_ID\n"
            tsv_content += f"{sample_name}\tNo_spa_type_detected\t\t\t\n"
        
        with open(output_dir / "spa_typing_report.tsv", 'w') as f:
            f.write(tsv_content)

    def generate_html_report(self, hits: List[Dict], sample_name: str, output_dir: Path):
        """Generate beautiful HTML report in MLST style"""
        analysis = self._analyze_spa_results(hits)
        
        html_content = f'''<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>STAPHSCOPE - spa Typing Analysis Report</title>
    <style>
        * {{
            margin: 0;
            padding: 0;
            box-sizing: border-box;
        }}
        
        body {{
            background: linear-gradient(135deg, #1e3c72 0%, #2a5298 50%, #7e22ce 100%);
            font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif;
            color: #ffffff;
            padding: 20px;
            min-height: 100vh;
        }}
        
        .container {{
            max-width: 1400px;
            margin: 0 auto;
        }}
        
        .header {{
            text-align: center;
            margin-bottom: 30px;
        }}
        
        .ascii-container {{
            background: rgba(0, 0, 0, 0.7);
            padding: 20px;
            border-radius: 15px;
            margin-bottom: 20px;
            box-shadow: 0 8px 32px rgba(0, 0, 0, 0.4);
            border: 2px solid rgba(0, 255, 0, 0.3);
        }}
        
        .ascii-art {{
            font-family: 'Courier New', monospace;
            font-size: 10px;
            line-height: 1.1;
            white-space: pre;
            color: #00ff00;
            text-shadow: 0 0 10px rgba(0, 255, 0, 0.5);
            overflow-x: auto;
        }}
        
        .report-section {{
            background: rgba(255, 255, 255, 0.95);
            color: #1f2937;
            padding: 25px;
            border-radius: 10px;
            margin-bottom: 20px;
            box-shadow: 0 4px 15px rgba(0, 0, 0, 0.2);
        }}
        
        .report-section h2 {{
            color: #1e3a8a;
            border-bottom: 3px solid #3b82f6;
            padding-bottom: 10px;
            margin-bottom: 20px;
            font-size: 24px;
        }}
        
        .report-section h3 {{
            color: #1e40af;
            margin-top: 20px;
            margin-bottom: 10px;
            font-size: 18px;
        }}
        
        .metrics-grid {{
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(250px, 1fr));
            gap: 20px;
            margin-top: 15px;
        }}
        
        .metric-card {{
            background: linear-gradient(135deg, #8b5cf6 0%, #6d28d9 100%);
            color: white;
            padding: 20px;
            border-radius: 8px;
            box-shadow: 0 4px 12px rgba(0, 0, 0, 0.15);
        }}
        
        .metric-label {{
            font-size: 14px;
            opacity: 0.9;
            margin-bottom: 5px;
        }}
        
        .metric-value {{
            font-size: 24px;
            font-weight: bold;
        }}
        
        .primary-result {{
            background: linear-gradient(135deg, #10b981 0%, #059669 100%);
            color: white;
            padding: 25px;
            border-radius: 10px;
            margin: 20px 0;
            box-shadow: 0 4px 15px rgba(0, 0, 0, 0.2);
        }}
        
        .repeat-structure {{
            font-family: 'Courier New', monospace;
            background: #1f2937;
            color: #fbbf24;
            padding: 12px 16px;
            border-radius: 6px;
            margin: 10px 0;
            font-size: 16px;
            font-weight: bold;
            border-left: 4px solid #f59e0b;
        }}
        
        .summary-table {{
            width: 100%;
            border-collapse: collapse;
            margin-top: 20px;
            font-size: 14px;
        }}
        
        .summary-table th {{
            background: linear-gradient(135deg, #3b82f6 0%, #1e40af 100%);
            color: white;
            padding: 12px;
            text-align: left;
            font-weight: bold;
        }}
        
        .summary-table td {{
            padding: 12px;
            border-bottom: 1px solid #e5e7eb;
        }}
        
        .summary-table tr:nth-child(even) {{
            background-color: #f8fafc;
        }}
        
        .summary-table tr:hover {{
            background-color: #e0f2fe;
        }}
        
        .detail-table {{
            width: 100%;
            border-collapse: collapse;
            margin-top: 20px;
            font-size: 13px;
        }}
        
        .detail-table th {{
            background: linear-gradient(135deg, #6366f1 0%, #4f46e5 100%);
            color: white;
            padding: 10px;
            text-align: left;
            font-weight: bold;
        }}
        
        .detail-table td {{
            padding: 10px;
            border-bottom: 1px solid #e5e7eb;
        }}
        
        .detail-table tr:nth-child(even) {{
            background-color: #f8fafc;
        }}
        
        .spa-type-cell {{
            font-weight: bold;
            color: #1e40af;
        }}
        
        .repeat-cell {{
            font-family: 'Courier New', monospace;
            background-color: #f0f9ff;
            color: #0369a1;
            font-weight: bold;
        }}
        
        .footer {{
            text-align: center;
            margin-top: 30px;
            padding: 20px;
            background: rgba(0, 0, 0, 0.3);
            border-radius: 10px;
            font-size: 14px;
        }}
        
        .timestamp {{
            color: #fbbf24;
            font-weight: bold;
        }}
        
        .authorship {{
            margin-top: 15px;
            padding: 15px;
            background: rgba(255, 255, 255, 0.1);
            border-radius: 8px;
            font-size: 12px;
        }}
        
        @media (max-width: 768px) {{
            .ascii-art {{
                font-size: 6px;
            }}
            .metrics-grid {{
                grid-template-columns: 1fr;
            }}
        }}
    </style>
</head>
<body>
    <div class="container">
        <div class="header">
            <div class="ascii-container">
                <div class="ascii-art">███████╗████████╗ █████╗ ██████╗ ██╗  ██╗███████╗ ██████╗ ██████╗ ██████╗ ███████╗
██╔════╝╚══██╔══╝██╔══██╗██╔══██╗██║  ██║██╔════╝██╔════╝██╔═══██╗██╔══██╗██╔════╝
███████╗   ██║   ███████║██████╔╝███████║███████╗██║     ██║   ██║██████╔╝█████╗  
╚════██║   ██║   ██╔══██║██╔═══╝ ██╔══██║╚════██║██║     ██║   ██║██╔═══╝ ██╔══╝  
███████║   ██║   ██║  ██║██║     ██║  ██║███████║╚██████╗╚██████╔╝██║     ███████╗
╚══════╝   ╚═╝   ╚═╝  ╚═╝╚═╝     ╚═╝  ╚═╝╚══════╝ ╚═════╝ ╚═════╝ ╚═╝     ╚══════╝</div>
            </div>
        </div>
        
        <div class="report-section">
            <h2>📊 Sample Information</h2>
            <div class="metrics-grid">
                <div class="metric-card">
                    <div class="metric-label">Sample Name</div>
                    <div class="metric-value">{sample_name}</div>
                </div>
                <div class="metric-card">
                    <div class="metric-label">Analysis Date</div>
                    <div class="metric-value">{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}</div>
                </div>
                <div class="metric-card">
                    <div class="metric-label">Total Hits</div>
                    <div class="metric-value">{analysis['total_hits']}</div>
                </div>
                <div class="metric-card">
                    <div class="metric-label">Unique Types</div>
                    <div class="metric-value">{analysis['total_types']}</div>
                </div>
            </div>
        </div>
'''
        
        if hits:
            primary_hit = hits[0]
            html_content += f'''
        <div class="report-section">
            <h2>🎯 Primary spa Typing Result</h2>
            <div class="primary-result">
                <div style="font-size: 14px; opacity: 0.9; margin-bottom: 10px;">PRIMARY SPA TYPE</div>
                <div style="font-size: 36px; font-weight: bold; margin-bottom: 15px;">{primary_hit['spa_type']}</div>
                
                <div style="margin: 15px 0;">
                    <div style="font-size: 14px; opacity: 0.9;">Repeat Structure</div>
                    <div class="repeat-structure">{primary_hit['repeats']}</div>
                </div>
                
                <div class="metrics-grid">
                    <div class="metric-card">
                        <div class="metric-label">Repeat Count</div>
                        <div class="metric-value">{primary_hit['repeat_count']}</div>
                    </div>
                    <div class="metric-card">
                        <div class="metric-label">Contig ID</div>
                        <div class="metric-value">{primary_hit['contig_id']}</div>
                    </div>
                </div>
            </div>
        </div>
'''
        
        if analysis['spa_types']:
            html_content += '''
        <div class="report-section">
            <h2>📊 All spa Types Detected</h2>
            <table class="summary-table">
                <thead>
                    <tr>
                        <th>spa Type</th>
                        <th>Repeat Structure</th>
                        <th>Repeat Count</th>
                        <th>Frequency</th>
                    </tr>
                </thead>
                <tbody>
'''
            for spa_type, type_info in analysis['spa_types'].items():
                html_content += f'''
                    <tr>
                        <td class="spa-type-cell">{spa_type}</td>
                        <td class="repeat-cell">{type_info['repeats']}</td>
                        <td>{type_info['repeat_count']}</td>
                        <td>{type_info['count']}</td>
                    </tr>
'''
            html_content += '''
                </tbody>
            </table>
        </div>
'''
        
        if hits:
            html_content += '''
        <div class="report-section">
            <h2>🔍 Detailed Results</h2>
            <table class="detail-table">
                <thead>
                    <tr>
                        <th>Sequence Name</th>
                        <th>spa Type</th>
                        <th>Repeat Structure</th>
                        <th>Repeat Count</th>
                        <th>Contig ID</th>
                    </tr>
                </thead>
                <tbody>
'''
            for hit in hits:
                html_content += f'''
                    <tr>
                        <td>{hit['sequence_name'][:80]}{'...' if len(hit['sequence_name']) > 80 else ''}</td>
                        <td class="spa-type-cell">{hit['spa_type']}</td>
                        <td class="repeat-cell">{hit['repeats']}</td>
                        <td>{hit['repeat_count']}</td>
                        <td>{hit['contig_id']}</td>
                    </tr>
'''
            html_content += '''
                </tbody>
            </table>
        </div>
'''
        
        html_content += f'''
        <div class="footer">
            <p><strong>STAPHSCOPE</strong> - spa Typing Analysis Module</p>
            <p class="timestamp">Generated: {datetime.now().strftime("%Y-%m-%d %H:%M:%S")}</p>
            <div class="authorship">
                <p><strong>Technical Support & Inquiries:</strong></p>
                <p>Author: Brown Beckley | GitHub: bbeckley-hub</p>
                <p>Email: brownbeckley94@gmail.com</p>
                <p>Affiliation: University of Ghana Medical School - Department of Medical Biochemistry</p>
            </div>
        </div>
    </div>
</body>
</html>'''
        
        with open(output_dir / "spa_typing_report.html", 'w', encoding='utf-8') as f:
            f.write(html_content)

    def _analyze_spa_results(self, hits: List[Dict]) -> Dict[str, Any]:
        """Analyze spa typing results for summary reporting"""
        analysis = {
            'total_hits': len(hits),
            'total_types': 0,
            'most_common_type': 'None',
            'spa_types': {},
        }
        
        for hit in hits:
            spa_type = hit['spa_type']
            repeats = hit['repeats']
            
            if spa_type not in analysis['spa_types']:
                analysis['spa_types'][spa_type] = {
                    'count': 0,
                    'repeats': repeats,
                    'repeat_count': hit['repeat_count']
                }
            analysis['spa_types'][spa_type]['count'] += 1
        
        analysis['total_types'] = len(analysis['spa_types'])
        
        if analysis['spa_types']:
            analysis['most_common_type'] = max(analysis['spa_types'].items(), 
                                             key=lambda x: x[1]['count'])[0]
        
        return analysis

    def create_spa_summary(self, all_results: Dict[str, Dict], output_dir: Path):
        """Create comprehensive spa typing summary files for all samples"""
        print("📊 Creating spa typing summary files...")
        
        # Create TSV summary
        self.create_spa_tsv_summary(all_results, output_dir)
        
        # Create HTML summary
        self.create_spa_html_summary(all_results, output_dir)
        
        print("✅ spa typing summary files created successfully!")

    def create_spa_tsv_summary(self, all_results: Dict[str, Dict], output_dir: Path):
        """Create TSV summary file with all samples"""
        summary_file = output_dir / "spa_summary.tsv"
        
        with open(summary_file, 'w') as f:
            # Write header
            f.write("Sample\tspa_Type\tRepeat_Pattern\tRepeat_Count\tStatus\n")
            
            # Write data for each sample
            for sample_name, result in all_results.items():
                if result['hits']:
                    primary_hit = result['hits'][0]
                    f.write(f"{sample_name}\t{primary_hit['spa_type']}\t{primary_hit['repeats']}\t{primary_hit['repeat_count']}\t{result['status']}\n")
                else:
                    f.write(f"{sample_name}\tNo_spa_type_detected\t\t\t{result['status']}\n")
        
        print(f"📄 TSV summary created: {summary_file}")

    def create_spa_html_summary(self, all_results: Dict[str, Dict], output_dir: Path):
        """Create beautiful HTML summary in MLST style"""
        summary_file = output_dir / "spa_summary.html"
        
        html_content = f'''<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>STAPHSCOPE - Batch spa Typing Summary</title>
    <style>
        * {{
            margin: 0;
            padding: 0;
            box-sizing: border-box;
        }}
        
        body {{
            background: linear-gradient(135deg, #1e3c72 0%, #2a5298 50%, #7e22ce 100%);
            font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif;
            color: #ffffff;
            padding: 20px;
            min-height: 100vh;
        }}
        
        .container {{
            max-width: 1800px;
            margin: 0 auto;
        }}
        
        .header {{
            text-align: center;
            margin-bottom: 30px;
        }}
        
        .ascii-container {{
            background: rgba(0, 0, 0, 0.7);
            padding: 20px;
            border-radius: 15px;
            margin-bottom: 20px;
            box-shadow: 0 8px 32px rgba(0, 0, 0, 0.4);
            border: 2px solid rgba(0, 255, 0, 0.3);
        }}
        
        .ascii-art {{
            font-family: 'Courier New', monospace;
            font-size: 10px;
            line-height: 1.1;
            white-space: pre;
            color: #00ff00;
            text-shadow: 0 0 10px rgba(0, 255, 0, 0.5);
            overflow-x: auto;
        }}
        
        .report-section {{
            background: rgba(255, 255, 255, 0.95);
            color: #1f2937;
            padding: 25px;
            border-radius: 10px;
            margin-bottom: 20px;
            box-shadow: 0 4px 15px rgba(0, 0, 0, 0.2);
        }}
        
        .report-section h2 {{
            color: #1e3a8a;
            border-bottom: 3px solid #3b82f6;
            padding-bottom: 10px;
            margin-bottom: 20px;
            font-size: 24px;
        }}
        
        .summary-table {{
            width: 100%;
            border-collapse: collapse;
            margin-top: 20px;
            font-size: 14px;
        }}
        
        .summary-table th {{
            background: linear-gradient(135deg, #3b82f6 0%, #1e40af 100%);
            color: white;
            padding: 12px;
            text-align: left;
            font-weight: bold;
            position: sticky;
            top: 0;
        }}
        
        .summary-table td {{
            padding: 10px;
            border-bottom: 1px solid #e5e7eb;
        }}
        
        .summary-table tr:nth-child(even) {{
            background-color: #f8fafc;
        }}
        
        .summary-table tr:hover {{
            background-color: #e0f2fe;
        }}
        
        .spa-type-cell {{
            font-weight: bold;
            color: #1e40af;
        }}
        
        .repeat-cell {{
            font-family: 'Courier New', monospace;
            background-color: #f0f9ff;
            color: #0369a1;
            font-weight: bold;
        }}
        
        .success {{
            color: #16a34a;
            font-weight: bold;
        }}
        
        .failed {{
            color: #dc2626;
            font-weight: bold;
        }}
        
        .stats-grid {{
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(200px, 1fr));
            gap: 15px;
            margin-bottom: 20px;
        }}
        
        .stat-card {{
            background: linear-gradient(135deg, #8b5cf6 0%, #6d28d9 100%);
            color: white;
            padding: 15px;
            border-radius: 8px;
            text-align: center;
            box-shadow: 0 4px 12px rgba(0, 0, 0, 0.15);
        }}
        
        .stat-value {{
            font-size: 24px;
            font-weight: bold;
            margin-bottom: 5px;
        }}
        
        .stat-label {{
            font-size: 12px;
            opacity: 0.9;
        }}
        
        .footer {{
            text-align: center;
            margin-top: 30px;
            padding: 20px;
            background: rgba(0, 0, 0, 0.3);
            border-radius: 10px;
            font-size: 14px;
        }}
        
        .timestamp {{
            color: #fbbf24;
            font-weight: bold;
        }}
        
        @media (max-width: 768px) {{
            .ascii-art {{
                font-size: 6px;
            }}
            .summary-table {{
                font-size: 12px;
            }}
            .summary-table th,
            .summary-table td {{
                padding: 6px;
            }}
        }}
    </style>
</head>
<body>
    <div class="container">
        <div class="header">
            <div class="ascii-container">
                <div class="ascii-art">███████╗████████╗ █████╗ ██████╗ ██╗  ██╗███████╗ ██████╗ ██████╗ ██████╗ ███████╗
██╔════╝╚══██╔══╝██╔══██╗██╔══██╗██║  ██║██╔════╝██╔════╝██╔═══██╗██╔══██╗██╔════╝
███████╗   ██║   ███████║██████╔╝███████║███████╗██║     ██║   ██║██████╔╝█████╗  
╚════██║   ██║   ██╔══██║██╔═══╝ ██╔══██║╚════██║██║     ██║   ██║██╔═══╝ ██╔══╝  
███████║   ██║   ██║  ██║██║     ██║  ██║███████║╚██████╗╚██████╔╝██║     ███████╗
╚══════╝   ╚═╝   ╚═╝  ╚═╝╚═╝     ╚═╝  ╚═╝╚══════╝ ╚═════╝ ╚═════╝ ╚═╝     ╚══════╝</div>
            </div>
        </div>
        
        <div class="report-section">
            <h2>📊 Batch spa Typing Summary - All Samples</h2>
            
            <div class="stats-grid">
                <div class="stat-card">
                    <div class="stat-value">{len(all_results)}</div>
                    <div class="stat-label">SAMPLES PROCESSED</div>
                </div>
                <div class="stat-card">
                    <div class="stat-value">{sum(1 for r in all_results.values() if r['hits'])}</div>
                    <div class="stat-label">SUCCESSFUL</div>
                </div>
                <div class="stat-card">
                    <div class="stat-value">{len(set(hit['spa_type'] for r in all_results.values() for hit in r['hits']))}</div>
                    <div class="stat-label">UNIQUE TYPES</div>
                </div>
                <div class="stat-card">
                    <div class="stat-value">{datetime.now().strftime('%Y-%m-%d')}</div>
                    <div class="stat-label">ANALYSIS DATE</div>
                </div>
            </div>
            
            <table class="summary-table">
                <thead>
                    <tr>
                        <th>Sample</th>
                        <th>spa Type</th>
                        <th>Repeat Pattern</th>
                        <th>Repeat Count</th>
                        <th>Status</th>
                    </tr>
                </thead>
                <tbody>
'''
        
        for sample_name, result in sorted(all_results.items()):
            if result['hits']:
                primary_hit = result['hits'][0]
                html_content += f'''
                    <tr>
                        <td><strong>{sample_name}</strong></td>
                        <td class="spa-type-cell">{primary_hit['spa_type']}</td>
                        <td class="repeat-cell">{primary_hit['repeats']}</td>
                        <td>{primary_hit['repeat_count']}</td>
                        <td class="success">✓ Success</td>
                    </tr>
'''
            else:
                html_content += f'''
                    <tr>
                        <td><strong>{sample_name}</strong></td>
                        <td colspan="2">No spa type detected</td>
                        <td>-</td>
                        <td class="failed">✗ Failed</td>
                    </tr>
'''
        
        html_content += f'''
                </tbody>
            </table>
        </div>
        
        <div class="footer">
            <p><strong>STAPHSCOPE</strong> - Batch spa Typing Analysis Summary</p>
            <p class="timestamp">Generated: {datetime.now().strftime("%Y-%m-%d %H:%M:%S")}</p>
        </div>
    </div>
</body>
</html>'''
        
        with open(summary_file, 'w', encoding='utf-8') as f:
            f.write(html_content)

    def run_spa_typing_batch(self, input_path: str, output_dir: str):
        """Run spa typing analysis on all FASTA files"""
        print("🚀 Starting spa typing analysis...")
        
        # Create output directory
        output_path = Path(output_dir)
        output_path.mkdir(parents=True, exist_ok=True)
        
        # Find all FASTA files
        fasta_files = self.find_fasta_files(input_path)
        
        if not fasta_files:
            print("❌ No FASTA files found!")
            return
        
        print(f"📁 Found {len(fasta_files)} FASTA files to process")
        
        # Process files in parallel
        all_results = {}
        with ThreadPoolExecutor(max_workers=self.threads) as executor:
            future_to_file = {
                executor.submit(self.run_spa_typing_single, file, output_path): file 
                for file in fasta_files
            }
            
            for future in as_completed(future_to_file):
                file = future_to_file[future]
                try:
                    result = future.result()
                    all_results[result['sample']] = result
                except Exception as e:
                    print(f"❌ Error processing {file}: {e}")
                    all_results[file.name] = self.get_fallback_results(file.name)
        
        # Create summary files
        self.create_spa_summary(all_results, output_path)
        
        print(f"✅ spa typing analysis completed! Processed {len(all_results)} samples")
        print(f"📁 Results saved to: {output_path}")

def main():
    parser = argparse.ArgumentParser(
        description='STAPHSCOPE spa Typing Module - Comprehensive spa typing analysis with beautiful HTML reports',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python spa_typing_module.py -i sample.fasta -o results/
  python spa_typing_module.py -i /path/to/fasta/files -o results/ -t 8
  python spa_typing_module.py -i *.fasta -o spa_results/
        """
    )
    
    parser.add_argument('-i', '--input', required=True,
                       help='Input FASTA file or directory containing FASTA files')
    parser.add_argument('-o', '--output', required=True,
                       help='Output directory for results')
    parser.add_argument('-t', '--threads', type=int, default=4,
                       help='Number of threads to use (default: 4)')
    
    args = parser.parse_args()
    
    if not os.path.exists(args.input):
        print(f"❌ Error: Input path '{args.input}' does not exist!")
        sys.exit(1)
    
    analyzer = SpaTypingAnalyzer(threads=args.threads)
    analyzer.run_spa_typing_batch(args.input, args.output)

if __name__ == "__main__":
    main()

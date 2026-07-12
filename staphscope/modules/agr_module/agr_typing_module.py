#!/usr/bin/env python3
"""
StaphScope Agr Typing Module - Using agrVATE
Author: Brown Beckley <brownbeckley94@gmail.com>
Date: 2026-07-10
Affiliation: University of Ghana Medical School - Department of Medical Biochemistry
MIT
"""

import os
import sys
import glob
import argparse
import subprocess
import logging
import shutil
import json
import re
import random
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed
from typing import List, Dict, Any, Optional
from datetime import datetime


class AgrTypingAnalyzer:
    """Comprehensive agr typing analyzer using agrVATE with beautiful HTML reporting"""

    def __init__(self, threads: int = 4):
        self.threads = threads
        self.logger = self._setup_logging()

        # Science quotes for rotation
        self.science_quotes = [
            {"text": "The important thing is not to stop questioning. Curiosity has its own reason for existing.",
             "author": "Albert Einstein"},
            {"text": "Science is not only a disciple of reason but also one of romance and passion.",
             "author": "Stephen Hawking"},
            {"text": "Somewhere, something incredible is waiting to be known.", "author": "Carl Sagan"},
            {"text": "The good thing about science is that it's true whether or not you believe in it.",
             "author": "Neil deGrasse Tyson"},
            {"text": "In science, there are no shortcuts to truth.", "author": "Karl Popper"},
            {"text": "Science knows no country, because knowledge belongs to humanity.", "author": "Louis Pasteur"},
            {"text": "The science of today is the technology of tomorrow.", "author": "Edward Teller"},
            {"text": "Nothing in life is to be feared, it is only to be understood.", "author": "Marie Curie"},
            {"text": "Research is what I'm doing when I don't know what I'm doing.",
             "author": "Wernher von Braun"},
            {"text": "The universe is not required to be in perfect harmony with human ambition.",
             "author": "Carl Sagan"},
            {"text": "Staphscope represents the convergence of genomic surveillance and clinical diagnostics, "
                    "transforming raw sequences into actionable insights for infection control.",
             "author": "Brown Beckley"},
            {"text": "In the battle against antimicrobial resistance, tools like Staphscope are our eyes and ears, "
                    "revealing the genetic blueprints of resistant pathogens.", "author": "Brown Beckley"},
            {"text": "Staphscope isn't just a tool; it's a comprehensive system that bridges the gap between "
                    "sequencing data and public health action.", "author": "Brown Beckley"},
            {"text": "Through Staphscope, we turn the complexity of bacterial genomes into clear, interpretable "
                    "reports, empowering clinicians and researchers alike.", "author": "Brown Beckley"},
            {"text": "Staphscope is a testament to the power of bioinformatics in the modern era, making advanced "
                    "pathogen typing accessible to all.", "author": "Brown Beckley"}
        ]

    def get_random_quote(self):
        """Get a random science quote"""
        return random.choice(self.science_quotes)

    def _setup_logging(self):
        logging.basicConfig(
            level=logging.INFO,
            format='%(asctime)s - %(levelname)s - %(message)s'
        )
        return logging.getLogger(__name__)

    def find_agrvate(self) -> Optional[Path]:
        """Find agrVATE executable in PATH or common system locations"""
        print("🔍 Searching for agrVATE...")

        # First check if agrVATE is in PATH
        agrvate_path = shutil.which('agrvate')
        if agrvate_path:
            print(f"✅ Found agrVATE at: {agrvate_path}")
            return Path(agrvate_path)

        # Check common system locations (non-conda)
        system_locations = [
            '/usr/local/bin/agrvate',
            '/usr/bin/agrvate',
            Path.home() / '.local' / 'bin' / 'agrvate',
        ]
        for loc in system_locations:
            loc_path = Path(loc)
            if loc_path.exists():
                print(f"✅ Found agrVATE at: {loc_path}")
                return loc_path

        print("❌ agrVATE not found. Please ensure agrVATE is installed and in your PATH.")
        print("💡 Installation: conda install -c bioconda agrvate")
        return None

    def find_fasta_files(self, input_path: str) -> List[Path]:
        """Find all FASTA files using glob patterns"""
        if os.path.isfile(input_path):
            return [Path(input_path)]

        patterns = [
            input_path,
            f"{input_path}/*.fna",
            f"{input_path}/*.fasta",
            f"{input_path}/*.fa",
            "*.fna",
            "*.fasta",
            "*.fa"
        ]

        fasta_files = []
        for pattern in patterns:
            matched_files = glob.glob(pattern)
            for file_path in matched_files:
                path = Path(file_path)
                if path.is_file() and path.suffix in ['.fna', '.fasta', '.fa']:
                    fasta_files.append(path)

        return sorted(list(set(fasta_files)))

    def run_agr_typing_single(self, input_file: Path, output_dir: Path) -> Dict[str, Any]:
        """Run agr typing analysis for a single file using agrVATE"""
        sample_stem = input_file.stem  # without extension
        print(f"🔬 Processing: {input_file.name}")

        # Create sample-specific output directory
        sample_output_dir = output_dir / sample_stem
        sample_output_dir.mkdir(parents=True, exist_ok=True)

        agrvate = self.find_agrvate()
        if not agrvate:
            error_result = self.get_fallback_results(input_file.name, "agrVATE not found")
            self.generate_output_files([], input_file.name, sample_output_dir, error=True)
            return error_result

        # Copy input file to sample output directory
        input_copy = sample_output_dir / input_file.name
        shutil.copy2(input_file, input_copy)

        # Build command (typing-only, force overwrite)
        cmd = [
            str(agrvate),
            "-t",
            "-f",
            "-i", input_copy.name
        ]

        print(f"   Running: {' '.join(cmd)} in {sample_output_dir}")

        try:
            result = subprocess.run(cmd, capture_output=True, text=True, timeout=300, cwd=sample_output_dir)

            # Clean up copied input file
            if input_copy.exists():
                input_copy.unlink()

            if result.returncode == 0:
                print(f"✅ agr typing completed for {input_file.name}")
                hits = self._parse_agr_output(sample_output_dir, sample_stem)
                self.generate_output_files(hits, input_file.name, sample_output_dir)
                return {
                    'sample': input_file.name,
                    'hits': hits,
                    'hit_count': len(hits),
                    'status': 'success'
                }
            else:
                print(f"❌ agr typing failed for {input_file.name}")
                print(f"   Error: {result.stderr}")
                error_result = self.get_fallback_results(input_file.name, "agrVATE failed")
                self.generate_output_files([], input_file.name, sample_output_dir, error=True)
                return error_result

        except subprocess.TimeoutExpired:
            print(f"⏰ agr typing timed out for {input_file.name}")
            if input_copy.exists():
                input_copy.unlink()
            error_result = self.get_fallback_results(input_file.name, "Timeout")
            self.generate_output_files([], input_file.name, sample_output_dir, error=True)
            return error_result
        except Exception as e:
            print(f"❌ Unexpected error for {input_file.name}: {e}")
            if input_copy.exists():
                input_copy.unlink()
            error_result = self.get_fallback_results(input_file.name, str(e))
            self.generate_output_files([], input_file.name, sample_output_dir, error=True)
            return error_result

    def _parse_agr_output(self, sample_dir: Path, sample_stem: str) -> List[Dict]:
        """Parse agrVATE summary file into structured data"""
        hits = []
        results_dir = sample_dir / f"{sample_stem}-results"
        if not results_dir.exists():
            print(f"⚠️ Results directory not found: {results_dir}")
            return hits

        summary_file = results_dir / f"{sample_stem}-summary.tab"
        if not summary_file.exists():
            print(f"⚠️ Summary file not found: {summary_file}")
            return hits

        try:
            with open(summary_file, 'r') as f:
                lines = f.readlines()

            if not lines:
                print(f"⚠️ Empty summary file: {summary_file}")
                return hits

            # Parse header and data
            header = lines[0].strip().split('\t')
            if header[0].startswith('#'):
                header[0] = header[0][1:]  # remove leading '#'

            for line in lines[1:]:
                line = line.strip()
                if not line:
                    continue
                parts = line.split('\t')
                if len(parts) >= 6:
                    record = dict(zip(header, parts))
                    group = record.get('agr_group', '').strip()
                    type_map = {'gp1': 'I', 'gp2': 'II', 'gp3': 'III', 'gp4': 'IV'}
                    agr_type = type_map.get(group, 'NA')

                    hit = {
                        'sample': sample_stem,
                        'agr_group': group,
                        'agr_type': agr_type,
                        'match_score': record.get('match_score', ''),
                        'canonical_agrD': record.get('canonical_agrD', ''),
                        'multiple_agr': record.get('multiple_agr', ''),
                        'frameshifts': record.get('frameshifts', '')
                    }
                    hits.append(hit)
                    print(f"   ✅ Parsed: {agr_type} (group {group})")
                else:
                    print(f"   ⚠️ Skipping malformed line: {line}")

            # Clean up agrVATE leftovers (keep only our generated reports)
            shutil.rmtree(results_dir, ignore_errors=True)
            # Remove error report if it exists
            error_report = sample_dir / f"{sample_stem}.fna-error-report.tab"
            if error_report.exists():
                error_report.unlink()

        except Exception as e:
            print(f"❌ Error parsing {summary_file}: {e}")

        print(f"   📊 Parsed {len(hits)} agr typing hit(s)")
        return hits

    def get_fallback_results(self, sample_name: str, error_msg: str = "Unknown error") -> Dict:
        """Fallback when agr typing fails"""
        return {
            'sample': sample_name,
            'hits': [],
            'hit_count': 0,
            'status': 'failed',
            'error': error_msg
        }

    def generate_output_files(self, hits: List[Dict], sample_name: str, output_dir: Path, error: bool = False):
        """Generate output files: HTML, TXT, and TSV"""
        self.generate_html_report(hits, sample_name, output_dir, error)
        self.generate_text_report(hits, sample_name, output_dir, error)
        self.generate_tsv_report(hits, sample_name, output_dir)

    def generate_text_report(self, hits: List[Dict], sample_name: str, output_dir: Path, error: bool = False):
        """Generate detailed text report"""
        analysis = self._analyze_agr_results(hits)

        report = f"""agr Typing Analysis Report
==========================

Sample: {sample_name}
Analysis Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}

SUMMARY:
--------
Total agr hits: {analysis['total_hits']}
Unique agr types: {analysis['total_types']}
Most common type: {analysis['most_common_type']}

"""

        if not error and hits:
            primary = hits[0]
            report += f"""PRIMARY RESULT:
-------------
agr Type: {primary['agr_type']} (Group: {primary['agr_group']})
Match Score: {primary['match_score']}
Canonical AgrD: {primary['canonical_agrD']}
Multiple Agr: {primary['multiple_agr']}
Frameshifts: {primary['frameshifts']}

"""
        elif error:
            report += f"ERROR: {analysis.get('error', 'Agr typing failed')}\n\n"

        if analysis['agr_types']:
            report += "DETAILED RESULTS:\n"
            report += "-----------------\n"
            for ag_type, info in analysis['agr_types'].items():
                report += f"Type: {ag_type}\n"
                report += f"  Group: {info['agr_group']}\n"
                report += f"  Match Score: {info['match_score']}\n"
                report += f"  Frameshifts: {info['frameshifts']}\n"
                report += f"  Count: {info['count']}\n\n"

        with open(output_dir / "agr_typing_report.txt", 'w') as f:
            f.write(report)

    def generate_tsv_report(self, hits: List[Dict], sample_name: str, output_dir: Path):
        """Generate simple TSV report"""
        if hits:
            tsv_content = "Sample\tagr_Type\tagr_Group\tMatch_Score\tCanonical_AgrD\tMultiple_Agr\tFrameshifts\n"
            for hit in hits:
                tsv_content += f"{sample_name}\t{hit['agr_type']}\t{hit['agr_group']}\t{hit['match_score']}\t{hit['canonical_agrD']}\t{hit['multiple_agr']}\t{hit['frameshifts']}\n"
        else:
            tsv_content = "Sample\tagr_Type\tagr_Group\tMatch_Score\tCanonical_AgrD\tMultiple_Agr\tFrameshifts\n"
            tsv_content += f"{sample_name}\tNA\tNA\tNA\tNA\tNA\tNA\n"

        with open(output_dir / "agr_typing_report.tsv", 'w') as f:
            f.write(tsv_content)

    def generate_html_report(self, hits: List[Dict], sample_name: str, output_dir: Path, error: bool = False):
        """Generate beautiful HTML report"""
        analysis = self._analyze_agr_results(hits)
        random_quote = self.get_random_quote()

        html_content = f'''<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>STAPHSCOPE - agr Typing Analysis Report</title>
    <style>
        * {{ margin: 0; padding: 0; box-sizing: border-box; }}
        body {{
            background: linear-gradient(135deg, #1e3c72 0%, #2a5298 50%, #7e22ce 100%);
            font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif;
            color: #ffffff;
            padding: 20px;
            min-height: 100vh;
        }}
        .container {{ max-width: 1400px; margin: 0 auto; }}
        .header {{ text-align: center; margin-bottom: 30px; }}
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
        .quote-container {{
            background: rgba(255, 255, 255, 0.1);
            backdrop-filter: blur(10px);
            padding: 20px;
            border-radius: 10px;
            margin-bottom: 30px;
            text-align: center;
            min-height: 100px;
            display: flex;
            flex-direction: column;
            justify-content: center;
            box-shadow: 0 4px 20px rgba(0, 0, 0, 0.3);
            border: 1px solid rgba(255, 255, 255, 0.2);
            transition: opacity 0.5s ease-in-out;
        }}
        .quote-text {{ font-size: 18px; font-style: italic; margin-bottom: 10px; color: #ffffff; }}
        .quote-author {{ font-size: 14px; color: #fbbf24; font-weight: bold; }}
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
        .metric-label {{ font-size: 14px; opacity: 0.9; margin-bottom: 5px; }}
        .metric-value {{ font-size: 24px; font-weight: bold; }}
        .primary-result {{
            background: linear-gradient(135deg, #10b981 0%, #059669 100%);
            color: white;
            padding: 25px;
            border-radius: 10px;
            margin: 20px 0;
            box-shadow: 0 4px 15px rgba(0, 0, 0, 0.2);
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
        .summary-table tr:nth-child(even) {{ background-color: #f8fafc; }}
        .summary-table tr:hover {{ background-color: #e0f2fe; }}
        .type-badge {{
            display: inline-block;
            padding: 4px 12px;
            border-radius: 20px;
            font-weight: bold;
            font-size: 14px;
        }}
        .type-I {{ background-color: #16a34a; color: white; }}
        .type-II {{ background-color: #2563eb; color: white; }}
        .type-III {{ background-color: #f59e0b; color: white; }}
        .type-IV {{ background-color: #dc2626; color: white; }}
        .type-NA {{ background-color: #6b7280; color: white; }}
        .footer {{
            text-align: center;
            margin-top: 30px;
            padding: 20px;
            background: rgba(0, 0, 0, 0.3);
            border-radius: 10px;
            font-size: 14px;
        }}
        .timestamp {{ color: #fbbf24; font-weight: bold; }}
        .authorship {{ margin-top: 15px; padding: 15px; background: rgba(255, 255, 255, 0.1); border-radius: 8px; font-size: 12px; }}
        .error-box {{
            background: #fef2f2;
            border-left: 4px solid #dc2626;
            padding: 15px;
            border-radius: 8px;
            margin: 20px 0;
            color: #991b1b;
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
            <div class="quote-container" id="quoteContainer">
                <div class="quote-text" id="quoteText">"{random_quote['text']}"</div>
                <div class="quote-author" id="quoteAuthor">— {random_quote['author']}</div>
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

        if error:
            html_content += f'''
        <div class="report-section">
            <div class="error-box">
                <strong>⚠️ Error:</strong> {analysis.get('error', 'Agr typing failed')}
            </div>
        </div>
'''
        elif hits:
            primary = hits[0]
            type_class = f"type-{primary['agr_type']}" if primary['agr_type'] != 'NA' else 'type-NA'
            html_content += f'''
        <div class="report-section">
            <h2>🎯 Primary agr Typing Result</h2>
            <div class="primary-result">
                <div style="font-size: 14px; opacity: 0.9; margin-bottom: 10px;">PRIMARY AGR TYPE</div>
                <div style="font-size: 36px; font-weight: bold; margin-bottom: 15px;">
                    <span class="type-badge {type_class}">{primary['agr_type']}</span>
                </div>
                <div class="metrics-grid">
                    <div class="metric-card">
                        <div class="metric-label">agr Group</div>
                        <div class="metric-value">{primary['agr_group']}</div>
                    </div>
                    <div class="metric-card">
                        <div class="metric-label">Match Score</div>
                        <div class="metric-value">{primary['match_score']}</div>
                    </div>
                    <div class="metric-card">
                        <div class="metric-label">Canonical AgrD</div>
                        <div class="metric-value">{primary['canonical_agrD']}</div>
                    </div>
                    <div class="metric-card">
                        <div class="metric-label">Multiple Agr</div>
                        <div class="metric-value">{primary['multiple_agr']}</div>
                    </div>
                </div>
            </div>
        </div>
'''

        if analysis['agr_types']:
            html_content += '''
        <div class="report-section">
            <h2>📊 All agr Types Detected</h2>
            <table class="summary-table">
                <thead>
                    <tr>
                        <th>agr Type</th>
                        <th>agr Group</th>
                        <th>Match Score</th>
                        <th>Canonical AgrD</th>
                        <th>Multiple Agr</th>
                        <th>Frameshifts</th>
                        <th>Frequency</th>
                    </tr>
                </thead>
                <tbody>
'''
            for ag_type, info in analysis['agr_types'].items():
                type_class = f"type-{ag_type}" if ag_type != 'NA' else 'type-NA'
                html_content += f'''
                    <tr>
                        <td><span class="type-badge {type_class}">{ag_type}</span></td>
                        <td>{info['agr_group']}</td>
                        <td>{info['match_score']}</td>
                        <td>{info['canonical_agrD']}</td>
                        <td>{info['multiple_agr']}</td>
                        <td>{info['frameshifts']}</td>
                        <td>{info['count']}</td>
                    </tr>
'''
            html_content += '''
                </tbody>
            </table>
        </div>
'''

        html_content += f'''
        <div class="footer">
            <p><strong>STAPHSCOPE</strong> - agr Typing Analysis Module</p>
            <p class="timestamp">Generated: {datetime.now().strftime("%Y-%m-%d %H:%M:%S")}</p>
            <div class="authorship">
                <p><strong>Technical Support & Inquiries:</strong></p>
                <p>Author: Brown Beckley | GitHub: bbeckley-hub</p>
                <p>Email: brownbeckley94@gmail.com</p>
                <p>Affiliation: University of Ghana Medical School - Department of Medical Biochemistry</p>
            </div>
        </div>
    </div>

    <script>
        const quotes = [
            {{"text": "The important thing is not to stop questioning. Curiosity has its own reason for existing.", "author": "Albert Einstein"}},
            {{"text": "Science is not only a disciple of reason but also one of romance and passion.", "author": "Stephen Hawking"}},
            {{"text": "Somewhere, something incredible is waiting to be known.", "author": "Carl Sagan"}},
            {{"text": "The good thing about science is that it's true whether or not you believe in it.", "author": "Neil deGrasse Tyson"}},
            {{"text": "In science, there are no shortcuts to truth.", "author": "Karl Popper"}},
            {{"text": "Science knows no country, because knowledge belongs to humanity.", "author": "Louis Pasteur"}},
            {{"text": "The science of today is the technology of tomorrow.", "author": "Edward Teller"}},
            {{"text": "Nothing in life is to be feared, it is only to be understood.", "author": "Marie Curie"}},
            {{"text": "Research is what I'm doing when I don't know what I'm doing.", "author": "Wernher von Braun"}},
            {{"text": "The universe is not required to be in perfect harmony with human ambition.", "author": "Carl Sagan"}}
        ];

        const quoteContainer = document.getElementById('quoteContainer');
        const quoteText = document.getElementById('quoteText');
        const quoteAuthor = document.getElementById('quoteAuthor');

        function getRandomQuote() {{
            return quotes[Math.floor(Math.random() * quotes.length)];
        }}

        function displayQuote() {{
            quoteContainer.style.opacity = '0';
            setTimeout(() => {{
                const quote = getRandomQuote();
                quoteText.textContent = '"' + quote.text + '"';
                quoteAuthor.textContent = '— ' + quote.author;
                quoteContainer.style.opacity = '1';
            }}, 500);
        }}

        setInterval(displayQuote, 10000);
    </script>
</body>
</html>'''

        with open(output_dir / "agr_typing_report.html", 'w', encoding='utf-8') as f:
            f.write(html_content)

    def _analyze_agr_results(self, hits: List[Dict]) -> Dict[str, Any]:
        """Analyze agr typing results for summary reporting"""
        analysis = {
            'total_hits': len(hits),
            'total_types': 0,
            'most_common_type': 'None',
            'agr_types': {},
            'error': None
        }

        if not hits:
            return analysis

        for hit in hits:
            ag_type = hit['agr_type']
            group = hit['agr_group']
            match_score = hit['match_score']
            canonical_agrD = hit['canonical_agrD']
            multiple_agr = hit['multiple_agr']
            frameshifts = hit['frameshifts']

            if ag_type not in analysis['agr_types']:
                analysis['agr_types'][ag_type] = {
                    'count': 0,
                    'agr_group': group,
                    'match_score': match_score,
                    'canonical_agrD': canonical_agrD,
                    'multiple_agr': multiple_agr,
                    'frameshifts': frameshifts
                }
            analysis['agr_types'][ag_type]['count'] += 1

        analysis['total_types'] = len(analysis['agr_types'])
        if analysis['agr_types']:
            analysis['most_common_type'] = max(analysis['agr_types'].items(),
                                               key=lambda x: x[1]['count'])[0]

        return analysis

    def create_agr_summary(self, all_results: Dict[str, Dict], output_dir: Path):
        """Create comprehensive agr typing summary files for all samples"""
        print("📊 Creating agr typing summary files...")

        self.create_agr_tsv_summary(all_results, output_dir)
        self.create_agr_html_summary(all_results, output_dir)
        self.create_agr_json_summary(all_results, output_dir)

        print("✅ agr typing summary files created successfully!")

    def create_agr_tsv_summary(self, all_results: Dict[str, Dict], output_dir: Path):
        """Create TSV summary file with all samples (without frameshift column)"""
        summary_file = output_dir / "agr_summary.tsv"

        with open(summary_file, 'w') as f:
            f.write("Sample\tagr_Type\tagr_Group\tMatch_Score\tCanonical_AgrD\tMultiple_Agr\tStatus\n")

            for sample_name, result in sorted(all_results.items()):
                if result['hits']:
                    primary = result['hits'][0]
                    f.write(f"{sample_name}\t{primary['agr_type']}\t{primary['agr_group']}\t{primary['match_score']}\t{primary['canonical_agrD']}\t{primary['multiple_agr']}\t{result['status']}\n")
                else:
                    f.write(f"{sample_name}\tNA\tNA\tNA\tNA\tNA\t{result['status']}\n")

        print(f"📄 TSV summary created: {summary_file}")

    def create_agr_html_summary(self, all_results: Dict[str, Dict], output_dir: Path):
        """Create beautiful HTML summary with statistics"""
        summary_file = output_dir / "agr_summary.html"

        total_samples = len(all_results)
        successful = sum(1 for r in all_results.values() if r['hits'])
        failed = total_samples - successful

        type_counts = {'I': 0, 'II': 0, 'III': 0, 'IV': 0, 'NA': 0}
        for result in all_results.values():
            for hit in result['hits']:
                typ = hit['agr_type']
                if typ in type_counts:
                    type_counts[typ] += 1

        random_quote = self.get_random_quote()

        html_content = f'''<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>STAPHSCOPE - Batch agr Typing Summary</title>
    <style>
        * {{ margin: 0; padding: 0; box-sizing: border-box; }}
        body {{
            background: linear-gradient(135deg, #1e3c72 0%, #2a5298 50%, #7e22ce 100%);
            font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif;
            color: #ffffff;
            padding: 20px;
            min-height: 100vh;
        }}
        .container {{ max-width: 1800px; margin: 0 auto; }}
        .header {{ text-align: center; margin-bottom: 30px; }}
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
        .quote-container {{
            background: rgba(255, 255, 255, 0.1);
            backdrop-filter: blur(10px);
            padding: 20px;
            border-radius: 10px;
            margin-bottom: 30px;
            text-align: center;
            min-height: 100px;
            display: flex;
            flex-direction: column;
            justify-content: center;
            box-shadow: 0 4px 20px rgba(0, 0, 0, 0.3);
            border: 1px solid rgba(255, 255, 255, 0.2);
            transition: opacity 0.5s ease-in-out;
        }}
        .quote-text {{ font-size: 18px; font-style: italic; margin-bottom: 10px; color: #ffffff; }}
        .quote-author {{ font-size: 14px; color: #fbbf24; font-weight: bold; }}
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
        .stat-value {{ font-size: 24px; font-weight: bold; margin-bottom: 5px; }}
        .stat-label {{ font-size: 12px; opacity: 0.9; }}
        .type-distribution {{
            display: flex;
            flex-wrap: wrap;
            gap: 15px;
            margin: 15px 0;
        }}
        .type-pill {{
            padding: 8px 16px;
            border-radius: 20px;
            font-weight: bold;
            font-size: 14px;
            color: white;
        }}
        .type-I {{ background-color: #16a34a; }}
        .type-II {{ background-color: #2563eb; }}
        .type-III {{ background-color: #f59e0b; }}
        .type-IV {{ background-color: #dc2626; }}
        .type-NA {{ background-color: #6b7280; }}
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
        .summary-table tr:nth-child(even) {{ background-color: #f8fafc; }}
        .summary-table tr:hover {{ background-color: #e0f2fe; }}
        .type-badge {{
            display: inline-block;
            padding: 3px 10px;
            border-radius: 12px;
            font-weight: bold;
            font-size: 12px;
            color: white;
        }}
        .type-I {{ background-color: #16a34a; }}
        .type-II {{ background-color: #2563eb; }}
        .type-III {{ background-color: #f59e0b; }}
        .type-IV {{ background-color: #dc2626; }}
        .type-NA {{ background-color: #6b7280; }}
        .success {{ color: #16a34a; font-weight: bold; }}
        .failed {{ color: #dc2626; font-weight: bold; }}
        .footer {{
            text-align: center;
            margin-top: 30px;
            padding: 20px;
            background: rgba(0, 0, 0, 0.3);
            border-radius: 10px;
            font-size: 14px;
        }}
        .timestamp {{ color: #fbbf24; font-weight: bold; }}
        @media (max-width: 768px) {{
            .ascii-art {{ font-size: 6px; }}
            .summary-table {{ font-size: 12px; }}
            .summary-table th, .summary-table td {{ padding: 6px; }}
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
            <div class="quote-container" id="quoteContainer">
                <div class="quote-text" id="quoteText">"{random_quote['text']}"</div>
                <div class="quote-author" id="quoteAuthor">— {random_quote['author']}</div>
            </div>
        </div>

        <div class="report-section">
            <h2>📊 Batch agr Typing Summary - All Samples</h2>

            <div class="stats-grid">
                <div class="stat-card">
                    <div class="stat-value">{total_samples}</div>
                    <div class="stat-label">SAMPLES PROCESSED</div>
                </div>
                <div class="stat-card">
                    <div class="stat-value">{successful}</div>
                    <div class="stat-label">SUCCESSFUL</div>
                </div>
                <div class="stat-card">
                    <div class="stat-value">{failed}</div>
                    <div class="stat-label">FAILED</div>
                </div>
                <div class="stat-card">
                    <div class="stat-value">{datetime.now().strftime('%Y-%m-%d')}</div>
                    <div class="stat-label">ANALYSIS DATE</div>
                </div>
            </div>

            <h3>Agr Type Distribution</h3>
            <div class="type-distribution">
                <span class="type-pill type-I">Type I: {type_counts['I']}</span>
                <span class="type-pill type-II">Type II: {type_counts['II']}</span>
                <span class="type-pill type-III">Type III: {type_counts['III']}</span>
                <span class="type-pill type-IV">Type IV: {type_counts['IV']}</span>
                <span class="type-pill type-NA">NA: {type_counts['NA']}</span>
            </div>

            <table class="summary-table">
                <thead>
                    <tr>
                        <th>Sample</th>
                        <th>agr Type</th>
                        <th>agr Group</th>
                        <th>Match Score</th>
                        <th>Canonical AgrD</th>
                        <th>Multiple Agr</th>
                        <th>Status</th>
                    </tr>
                </thead>
                <tbody>
'''

        for sample_name, result in sorted(all_results.items()):
            if result['hits']:
                primary = result['hits'][0]
                type_class = f"type-{primary['agr_type']}" if primary['agr_type'] != 'NA' else 'type-NA'
                html_content += f'''
                    <tr>
                        <td><strong>{sample_name}</strong></td>
                        <td><span class="type-badge {type_class}">{primary['agr_type']}</span></td>
                        <td>{primary['agr_group']}</td>
                        <td>{primary['match_score']}</td>
                        <td>{primary['canonical_agrD']}</td>
                        <td>{primary['multiple_agr']}</td>
                        <td class="success">✓ Success</td>
                    </tr>
'''
            else:
                html_content += f'''
                    <tr>
                        <td><strong>{sample_name}</strong></td>
                        <td><span class="type-badge type-NA">NA</span></td>
                        <td>NA</td>
                        <td>NA</td>
                        <td>NA</td>
                        <td>NA</td>
                        <td class="failed">✗ Failed</td>
                    </tr>
'''

        html_content += f'''
                </tbody>
            </table>
        </div>

        <div class="footer">
            <p><strong>STAPHSCOPE</strong> - Batch agr Typing Analysis Summary</p>
            <p class="timestamp">Generated: {datetime.now().strftime("%Y-%m-%d %H:%M:%S")}</p>
            <p>Github: bbeckley-hub | Email: brownbeckley94@gmail.com</p>
        </div>
    </div>

    <script>
        const quotes = [
            {{"text": "The important thing is not to stop questioning. Curiosity has its own reason for existing.", "author": "Albert Einstein"}},
            {{"text": "Science is not only a disciple of reason but also one of romance and passion.", "author": "Stephen Hawking"}},
            {{"text": "Somewhere, something incredible is waiting to be known.", "author": "Carl Sagan"}},
            {{"text": "The good thing about science is that it's true whether or not you believe in it.", "author": "Neil deGrasse Tyson"}},
            {{"text": "In science, there are no shortcuts to truth.", "author": "Karl Popper"}},
            {{"text": "Science knows no country, because knowledge belongs to humanity.", "author": "Louis Pasteur"}},
            {{"text": "The science of today is the technology of tomorrow.", "author": "Edward Teller"}},
            {{"text": "Nothing in life is to be feared, it is only to be understood.", "author": "Marie Curie"}},
            {{"text": "Research is what I'm doing when I don't know what I'm doing.", "author": "Wernher von Braun"}},
            {{"text": "The universe is not required to be in perfect harmony with human ambition.", "author": "Carl Sagan"}}
        ];

        const quoteContainer = document.getElementById('quoteContainer');
        const quoteText = document.getElementById('quoteText');
        const quoteAuthor = document.getElementById('quoteAuthor');

        function getRandomQuote() {{
            return quotes[Math.floor(Math.random() * quotes.length)];
        }}

        function displayQuote() {{
            quoteContainer.style.opacity = '0';
            setTimeout(() => {{
                const quote = getRandomQuote();
                quoteText.textContent = '"' + quote.text + '"';
                quoteAuthor.textContent = '— ' + quote.author;
                quoteContainer.style.opacity = '1';
            }}, 500);
        }}

        setInterval(displayQuote, 10000);
    </script>
</body>
</html>'''

        with open(summary_file, 'w', encoding='utf-8') as f:
            f.write(html_content)

        print(f"🌐 HTML summary created: {summary_file}")

    def create_agr_json_summary(self, all_results: Dict[str, Dict], output_dir: Path):
        """Create JSON summary file with all agr typing results"""
        print("📊 Creating agr typing JSON summary...")

        summary_file = output_dir / "agr_summary.json"

        json_summary = {
            "metadata": {
                "analysis_date": datetime.now().isoformat(),
                "total_samples": len(all_results),
                "analysis_type": "agr typing",
                "version": "1.0",
                "analysis_method": "agrVATE with typing-only mode"
            },
            "statistics": self._calculate_agr_json_statistics(all_results),
            "samples": {},
            "summary_by_agr_type": self._create_agr_type_summary(all_results)
        }

        for sample_name, result in all_results.items():
            if result['hits']:
                primary = result['hits'][0]
                json_summary["samples"][sample_name] = {
                    "status": result['status'],
                    "total_hits": result['hit_count'],
                    "primary_result": {
                        "agr_type": primary['agr_type'],
                        "agr_group": primary['agr_group'],
                        "match_score": primary['match_score'],
                        "canonical_agrD": primary['canonical_agrD'],
                        "multiple_agr": primary['multiple_agr'],
                        "frameshifts": primary['frameshifts']
                    },
                    "all_hits": result['hits']
                }
            else:
                json_summary["samples"][sample_name] = {
                    "status": result['status'],
                    "total_hits": 0,
                    "error": result.get('error', 'No agr type detected')
                }

        with open(summary_file, 'w', encoding='utf-8') as f:
            json.dump(json_summary, f, indent=2, ensure_ascii=False)

        print(f"📄 JSON summary created: {summary_file}")

    def _calculate_agr_json_statistics(self, all_results: Dict[str, Dict]) -> Dict:
        """Calculate statistics for JSON summary"""
        total_samples = len(all_results)
        successful = sum(1 for r in all_results.values() if r['hits'])
        failed = total_samples - successful

        type_counts = {'I': 0, 'II': 0, 'III': 0, 'IV': 0, 'NA': 0}
        group_counts = {}
        match_scores = []

        for result in all_results.values():
            for hit in result['hits']:
                typ = hit['agr_type']
                if typ in type_counts:
                    type_counts[typ] += 1
                group = hit['agr_group']
                group_counts[group] = group_counts.get(group, 0) + 1
                try:
                    match_scores.append(int(hit['match_score']))
                except ValueError:
                    pass

        return {
            "total_samples": total_samples,
            "successful_samples": successful,
            "failed_samples": failed,
            "success_rate": (successful / total_samples * 100) if total_samples > 0 else 0,
            "type_distribution": type_counts,
            "group_distribution": group_counts,
            "match_score_stats": {
                "average": sum(match_scores) / len(match_scores) if match_scores else 0,
                "min": min(match_scores) if match_scores else 0,
                "max": max(match_scores) if match_scores else 0
            },
            "most_common_type": max(type_counts.items(), key=lambda x: x[1])[0] if any(type_counts.values()) else "NA"
        }

    def _create_agr_type_summary(self, all_results: Dict[str, Dict]) -> Dict:
        """Create summary grouped by agr type"""
        type_summary = {}
        for result in all_results.values():
            for hit in result['hits']:
                typ = hit['agr_type']
                if typ not in type_summary:
                    type_summary[typ] = {
                        "samples": [],
                        "count": 0,
                        "groups": set(),
                        "match_scores": []
                    }
                type_summary[typ]["samples"].append(result['sample'])
                type_summary[typ]["count"] += 1
                type_summary[typ]["groups"].add(hit['agr_group'])
                try:
                    type_summary[typ]["match_scores"].append(int(hit['match_score']))
                except ValueError:
                    pass

        for typ in type_summary:
            type_summary[typ]["groups"] = list(type_summary[typ]["groups"])
            if type_summary[typ]["match_scores"]:
                type_summary[typ]["average_match_score"] = sum(type_summary[typ]["match_scores"]) / len(type_summary[typ]["match_scores"])
            else:
                type_summary[typ]["average_match_score"] = 0

        return type_summary

    def run_agr_typing_batch(self, input_path: str, output_dir: str):
        """Run agr typing analysis on all FASTA files"""
        print("🚀 Starting agr typing analysis...")

        output_path = Path(output_dir)
        output_path.mkdir(parents=True, exist_ok=True)

        fasta_files = self.find_fasta_files(input_path)

        if not fasta_files:
            print("❌ No FASTA files found!")
            return

        print(f"📁 Found {len(fasta_files)} FASTA files to process")

        all_results = {}
        with ThreadPoolExecutor(max_workers=self.threads) as executor:
            future_to_file = {
                executor.submit(self.run_agr_typing_single, file, output_path): file
                for file in fasta_files
            }

            for future in as_completed(future_to_file):
                file = future_to_file[future]
                try:
                    result = future.result()
                    all_results[result['sample']] = result
                except Exception as e:
                    print(f"❌ Error processing {file}: {e}")
                    all_results[file.name] = self.get_fallback_results(file.name, str(e))

        self.create_agr_summary(all_results, output_path)

        print(f"✅ agr typing analysis completed! Processed {len(all_results)} samples")
        print(f"📁 Results saved to: {output_path}")


def main():
    parser = argparse.ArgumentParser(
        description='STAPHSCOPE Agr Typing Module - Comprehensive agr typing analysis with beautiful HTML reports',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python agr_typing_module.py -i sample.fna -o results/
  python agr_typing_module.py -i /path/to/fasta/files -o results/ -t 8
  python agr_typing_module.py -i . -o agr_results/
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

    analyzer = AgrTypingAnalyzer(threads=args.threads)
    analyzer.run_agr_typing_batch(args.input, args.output)


if __name__ == "__main__":
    main()
#!/usr/bin/env python3
"""
MLST Module for StaphScope - Modular with wildcard support
Author: Brown Beckley <brownbeckley94@gmail.com>
"""

import os
import sys
import json
import glob
import argparse
import subprocess
from pathlib import Path
from typing import Dict, List, Optional
import pandas as pd
from datetime import datetime

class ModularMLSTAnalyzer:
    def __init__(self, database_dir: Path, script_dir: Path):
        self.database_dir = database_dir
        self.script_dir = script_dir
        # Correct MLST binary path - it's in bin/mlst
        self.mlst_bin = script_dir / "mlst"
        
    def find_fasta_files(self, input_path: str) -> List[Path]:
        """Find all FASTA files using glob patterns"""
        # If it's a specific file, return it directly
        if os.path.isfile(input_path):
            return [Path(input_path)]
        
        fasta_patterns = [
            input_path,
            f"{input_path}/*.fna", f"{input_path}/*.fasta",
            f"{input_path}/*.fa", f"{input_path}/*.fn",
            f"{input_path}/*.fna.gz", f"{input_path}/*.fasta.gz",
            f"{input_path}/*.fa.gz", f"{input_path}/*.gb",
            f"{input_path}/*.gbk", f"{input_path}/*.gbff"
        ]
        
        fasta_files = []
        for pattern in fasta_patterns:
            matched_files = glob.glob(pattern)
            for file_path in matched_files:
                path = Path(file_path)
                if path.is_file():
                    fasta_files.append(path)
        
        # Remove duplicates and sort
        fasta_files = sorted(list(set(fasta_files)))
        return fasta_files

    def run_mlst_single(self, input_file: Path, output_dir: Path, scheme: str = "saureus") -> Dict:
        """Run MLST analysis for a single file"""
        print(f"🔬 Processing: {input_file.name}")
        
        # Create sample-specific output directory
        sample_output_dir = output_dir / input_file.stem
        sample_output_dir.mkdir(parents=True, exist_ok=True)
        
        # Run MLST command - use perl to run the mlst script
        mlst_cmd = [
            "perl", str(self.mlst_bin),
            str(input_file),
            "--scheme", scheme,
            "--csv",
            "--nopath",
            "--quiet"
        ]
        
        try:
            result = subprocess.run(mlst_cmd, capture_output=True, text=True, check=True)
            
            # Parse the CSV output from stdout
            mlst_results = self.parse_mlst_stdout(result.stdout, input_file.name)
            
            # Add lineage information
            mlst_results.update(self.get_lineage_info(mlst_results.get('st', 'ND')))
            
            # Generate all report formats
            self.generate_mlst_reports(mlst_results, sample_output_dir, input_file.name)
            
            print(f"✅ Completed: {input_file.name} -> ST{mlst_results.get('st', 'ND')}")
            return mlst_results
            
        except subprocess.CalledProcessError as e:
            print(f"❌ MLST failed for {input_file.name}: {e}")
            error_result = self.get_fallback_results(input_file.name)
            self.generate_mlst_reports(error_result, sample_output_dir, input_file.name)
            return error_result

    def run_mlst_batch(self, input_path: str, output_dir: Path, scheme: str = "saureus") -> Dict[str, Dict]:
        """Run MLST analysis for multiple files"""
        print("🔍 Searching for FASTA files...")
        fasta_files = self.find_fasta_files(input_path)
        
        if not fasta_files:
            print("❌ No FASTA files found!")
            return {}
        
        print(f"📁 Found {len(fasta_files)} FASTA files")
        
        results = {}
        for fasta_file in fasta_files:
            result = self.run_mlst_single(fasta_file, output_dir, scheme)
            results[fasta_file.name] = result
        
        # Generate batch summary
        self.generate_batch_summary(results, output_dir)
        
        return results

    def parse_mlst_stdout(self, stdout: str, sample_name: str) -> Dict:
        """Parse MLST stdout in CSV format"""
        lines = stdout.strip().split('\n')
        if len(lines) < 2:
            return self.get_empty_results(sample_name)
        
        # Parse CSV header and data
        header = lines[0].split('\t')
        data = lines[1].split('\t')
        
        if len(header) != len(data):
            return self.get_empty_results(sample_name)
        
        result_dict = dict(zip(header, data))
        
        # Extract alleles
        alleles = {}
        allele_profile = result_dict.get('alleles', '')
        if allele_profile:
            # Format: arcC(1)-aroE(4)-glpF(1)-gmk(4)-pta(12)-tpi(1)-yqiL(10)
            allele_parts = allele_profile.split('-')
            for part in allele_parts:
                if '(' in part and ')' in part:
                    gene = part.split('(')[0]
                    allele = part.split('(')[1].rstrip(')')
                    alleles[gene] = allele
        
        return {
            "sample": sample_name,
            "st": result_dict.get('ST', 'ND'),
            "scheme": result_dict.get('scheme', 'saureus'),
            "alleles": alleles,
            "allele_profile": allele_profile,
            "confidence": "HIGH" if result_dict.get('ST', 'ND') != 'ND' else "LOW"
        }

    def get_lineage_info(self, st: str) -> Dict:
        """Get lineage information based on ST"""
        lineage_db = {
            '5': {
                "clonal_complex": "CC5",
                "classification": "Healthcare-associated MRSA",
                "geographic_distribution": "Global",
                "clinical_significance": "Major healthcare-associated lineage, often multidrug-resistant",
                "common_virulence": ["Enterotoxins", "Immune evasion cluster"],
                "outbreak_potential": "HIGH"
            },
            '8': {
                "clonal_complex": "CC8", 
                "classification": "Community-associated MRSA",
                "geographic_distribution": "Global",
                "clinical_significance": "Includes USA300 clone, often PVL-positive",
                "common_virulence": ["Panton-Valentine Leukocidin", "Enterotoxins"],
                "outbreak_potential": "HIGH"
            },
            '22': {
                "clonal_complex": "CC22",
                "classification": "Healthcare-associated MRSA",
                "geographic_distribution": "Europe, Global", 
                "clinical_significance": "Epidemic MRSA-15, major healthcare clone",
                "common_virulence": ["Enterotoxins", "Immune evasion cluster"],
                "outbreak_potential": "HIGH"
            },
            '30': {
                "clonal_complex": "CC30",
                "classification": "Healthcare and Community-associated",
                "geographic_distribution": "Global",
                "clinical_significance": "Includes EMRSA-16, often PVL-positive",
                "common_virulence": ["Panton-Valentine Leukocidin", "Enterotoxins"],
                "outbreak_potential": "HIGH"
            }
        }
        
        return lineage_db.get(st, {
            "clonal_complex": f"CC{st}" if st != 'ND' else 'Unknown',
            "classification": "Further characterization needed",
            "geographic_distribution": "Variable",
            "clinical_significance": "Lineage requires additional analysis",
            "common_virulence": ["Unknown"],
            "outbreak_potential": "UNKNOWN"
        })

    def get_empty_results(self, sample_name: str) -> Dict:
        """Return empty results structure"""
        return {
            "sample": sample_name,
            "st": "ND",
            "scheme": "saureus",
            "alleles": {},
            "allele_profile": "",
            "confidence": "LOW"
        }

    def get_fallback_results(self, sample_name: str) -> Dict:
        """Fallback when MLST fails"""
        return {
            "sample": sample_name,
            "st": "UNKNOWN",
            "scheme": "saureus",
            "alleles": {},
            "allele_profile": "",
            "confidence": "LOW",
            "error": "MLST analysis failed"
        }

    def generate_mlst_reports(self, mlst_results: Dict, output_dir: Path, sample_name: str):
        """Generate all report formats for a single sample"""
        # JSON report
        with open(output_dir / "mlst_report.json", 'w') as f:
            json.dump(mlst_results, f, indent=2)
        
        # Text report
        self.generate_text_report(mlst_results, output_dir)
        
        # TSV report
        self.generate_tsv_report(mlst_results, output_dir)
        
        # Beautiful HTML report
        self.generate_html_report(mlst_results, output_dir)

    def generate_text_report(self, mlst_results: Dict, output_dir: Path):
        """Generate detailed text report"""
        report = f"""MLST Analysis Report
===================

Sample: {mlst_results['sample']}
Analysis Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}

MLST TYPING RESULTS:
-------------------
Sequence Type (ST): {mlst_results['st']}
Scheme: {mlst_results['scheme']}
Confidence: {mlst_results['confidence']}

Allele Profile:
{mlst_results['allele_profile']}

Detailed Alleles:
"""
        for gene, allele in mlst_results['alleles'].items():
            report += f"- {gene}: {allele}\n"
        
        # Add lineage information if available
        if 'clonal_complex' in mlst_results:
            report += f"""
LINEAGE ANALYSIS:
-----------------
Clonal Complex: {mlst_results['clonal_complex']}
Classification: {mlst_results['classification']}
Geographic Distribution: {mlst_results['geographic_distribution']}
Clinical Significance: {mlst_results['clinical_significance']}
Outbreak Potential: {mlst_results['outbreak_potential']}

Common Virulence Factors:
"""
            for virulence in mlst_results['common_virulence']:
                report += f"- {virulence}\n"
        
        with open(output_dir / "mlst_detailed_report.txt", 'w') as f:
            f.write(report)

    def generate_tsv_report(self, mlst_results: Dict, output_dir: Path):
        """Generate TSV report"""
        tsv_content = f"Sample\tST\tScheme\tConfidence\tAllele_Profile\tClonal_Complex\tClassification\tOutbreak_Potential\n"
        tsv_content += f"{mlst_results['sample']}\t{mlst_results['st']}\t{mlst_results['scheme']}\t{mlst_results['confidence']}\t{mlst_results['allele_profile']}\t{mlst_results.get('clonal_complex', 'Unknown')}\t{mlst_results.get('classification', 'Unknown')}\t{mlst_results.get('outbreak_potential', 'UNKNOWN')}\n"
        
        with open(output_dir / "mlst_report.tsv", 'w') as f:
            f.write(tsv_content)

    def generate_html_report(self, mlst_results: Dict, output_dir: Path):
        """Generate beautiful HTML report in StaphScope style"""
        html_content = f'''<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>STAPHSCOPE - MLST Analysis Report</title>
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
        
        .allele-grid {{
            display: grid;
            grid-template-columns: repeat(auto-fill, minmax(150px, 1fr));
            gap: 15px;
            margin-top: 15px;
        }}
        
        .allele-card {{
            background: linear-gradient(135deg, #3b82f6 0%, #1e40af 100%);
            color: white;
            padding: 15px;
            border-radius: 8px;
            box-shadow: 0 2px 8px rgba(0, 0, 0, 0.1);
            text-align: center;
            font-weight: bold;
        }}
        
        .confidence-high {{
            color: #16a34a;
            font-weight: bold;
        }}
        
        .confidence-medium {{
            color: #f59e0b;
            font-weight: bold;
        }}
        
        .confidence-low {{
            color: #dc2626;
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
                    <div class="metric-value">{mlst_results['sample']}</div>
                </div>
                <div class="metric-card">
                    <div class="metric-label">Analysis Date</div>
                    <div class="metric-value">{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}</div>
                </div>
                <div class="metric-card">
                    <div class="metric-label">MLST Scheme</div>
                    <div class="metric-value">{mlst_results['scheme']}</div>
                </div>
            </div>
        </div>
        
        <div class="report-section">
            <h2>🧬 MLST Typing Results</h2>
            <div class="metrics-grid">
                <div class="metric-card">
                    <div class="metric-label">Sequence Type</div>
                    <div class="metric-value">ST{mlst_results['st']}</div>
                </div>
                <div class="metric-card">
                    <div class="metric-label">Confidence</div>
                    <div class="metric-value confidence-{mlst_results['confidence'].lower()}">{mlst_results['confidence']}</div>
                </div>
                <div class="metric-card">
                    <div class="metric-label">Clonal Complex</div>
                    <div class="metric-value">{mlst_results.get('clonal_complex', 'Unknown')}</div>
                </div>
            </div>
            
            <h3>Allele Profile</h3>
            <div style="background: #f8fafc; padding: 15px; border-radius: 8px; margin: 15px 0;">
                <code style="font-size: 16px; color: #1e40af;">{mlst_results['allele_profile']}</code>
            </div>
            
            <h3>Individual Alleles</h3>
            <div class="allele-grid">
'''
        
        # Add allele cards
        for gene, allele in mlst_results['alleles'].items():
            html_content += f'''                <div class="allele-card">
                    <div style="font-size: 12px; opacity: 0.9;">{gene}</div>
                    <div style="font-size: 18px;">{allele}</div>
                </div>
'''
        
        html_content += f'''            </div>
        </div>
        
        <div class="report-section">
            <h2>🌍 Lineage Information</h2>
            <div class="metrics-grid">
                <div class="metric-card">
                    <div class="metric-label">Classification</div>
                    <div class="metric-value">{mlst_results.get('classification', 'Unknown')}</div>
                </div>
                <div class="metric-card">
                    <div class="metric-label">Geographic Distribution</div>
                    <div class="metric-value">{mlst_results.get('geographic_distribution', 'Unknown')}</div>
                </div>
                <div class="metric-card">
                    <div class="metric-label">Outbreak Potential</div>
                    <div class="metric-value">{mlst_results.get('outbreak_potential', 'UNKNOWN')}</div>
                </div>
            </div>
            
            <h3>Clinical Significance</h3>
            <p style="margin: 15px 0; line-height: 1.6;">{mlst_results.get('clinical_significance', 'Further analysis required.')}</p>
            
            <h3>Common Virulence Factors</h3>
            <div style="display: grid; grid-template-columns: repeat(auto-fill, minmax(200px, 1fr)); gap: 10px; margin-top: 15px;">
'''
        
        for virulence in mlst_results.get('common_virulence', []):
            html_content += f'''                <div style="background: #e0f2fe; color: #0369a1; padding: 10px; border-radius: 6px; text-align: center; font-weight: bold;">
                    {virulence}
                </div>
'''
        
        html_content += '''            </div>
        </div>
        
        <div class="footer">
            <p><strong>STAPHSCOPE</strong> - MLST Analysis Module</p>
            <p class="timestamp">Generated: ''' + datetime.now().strftime("%Y-%m-%d %H:%M:%S") + '''</p>
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
        
        with open(output_dir / "mlst_report.html", 'w', encoding='utf-8') as f:
            f.write(html_content)

    def generate_batch_summary(self, results: Dict[str, Dict], output_dir: Path):
        """Generate batch summary report"""
        summary_data = []
        for sample, result in results.items():
            summary_data.append({
                'Sample': sample,
                'ST': result['st'],
                'Scheme': result['scheme'],
                'Confidence': result['confidence'],
                'Allele_Profile': result['allele_profile'],
                'Clonal_Complex': result.get('clonal_complex', 'Unknown'),
                'Classification': result.get('classification', 'Unknown'),
                'Outbreak_Potential': result.get('outbreak_potential', 'UNKNOWN')
            })
        
        # JSON summary
        with open(output_dir / "batch_mlst_summary.json", 'w') as f:
            json.dump(summary_data, f, indent=2)
        
        # TSV summary
        df = pd.DataFrame(summary_data)
        df.to_csv(output_dir / "batch_mlst_summary.tsv", sep='\t', index=False)
        
        print(f"📊 Batch summary generated: {len(results)} samples processed")

def main():
    parser = argparse.ArgumentParser(description='StaphScope Modular MLST Analyzer')
    parser.add_argument('-i', '--input', required=True, 
                       help='Input FASTA file or directory (supports wildcards)')
    parser.add_argument('-o', '--output-dir', required=True, 
                       help='Output directory')
    parser.add_argument('-db', '--database-dir', required=True,
                       help='Database directory')
    parser.add_argument('-sc', '--script-dir', required=True,
                       help='Script directory (contains mlst binary)')
    parser.add_argument('-s', '--scheme', default='saureus',
                       help='MLST scheme (default: saureus)')
    parser.add_argument('--batch', action='store_true',
                       help='Process multiple files')
    
    args = parser.parse_args()
    
    # Initialize analyzer
    analyzer = ModularMLSTAnalyzer(
        database_dir=Path(args.database_dir),
        script_dir=Path(args.script_dir)
    )
    
    # Create output directory
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    if args.batch:
        # Batch processing
        results = analyzer.run_mlst_batch(args.input, output_dir, args.scheme)
        print(f"🎉 Batch MLST completed! Processed {len(results)} samples")
    else:
        # Single file processing
        input_file = Path(args.input)
        if input_file.exists():
            result = analyzer.run_mlst_single(input_file, output_dir, args.scheme)
            print(f"🎉 MLST completed for {input_file.name}: ST{result.get('st', 'ND')}")
        else:
            print(f"❌ Input file not found: {args.input}")

if __name__ == "__main__":
    main()

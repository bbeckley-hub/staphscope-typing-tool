#!/usr/bin/env python3
"""
MLST Module for StaphScope - Complete with Beautiful HTML Reports
Author: Brown Beckley <brownbeckley94@gmail.com>
GitHub: bbeckley-hub
Affiliation: University of Ghana Medical School - Department of Medical Biochemistry
Date: 2025
Send a quick mail for any issues or further explanations.
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
        self.mlst_bin = script_dir / "mlst"
        
    def find_fasta_files(self, input_path: str) -> List[Path]:
        """Find all FASTA files using glob patterns"""
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
        
        return sorted(list(set(fasta_files)))

    def run_mlst_single(self, input_file: Path, output_dir: Path, scheme: str = "saureus") -> Dict:
        """Run MLST analysis for a single file"""
        print(f"🔬 Processing: {input_file.name}")
        
        # Create sample-specific output directory
        sample_output_dir = output_dir / input_file.stem
        sample_output_dir.mkdir(parents=True, exist_ok=True)
        
        # Save raw MLST output first
        raw_output_file = sample_output_dir / "mlst_raw_output.txt"
        
        # Run MLST command
        mlst_cmd = [
            "perl", str(self.mlst_bin),
            str(input_file),
            "--scheme", scheme,
            "--csv",
            "--nopath"
        ]
        
        try:
            # Run and capture output
            result = subprocess.run(mlst_cmd, capture_output=True, text=True, check=True)
            
            # Save raw output
            with open(raw_output_file, 'w') as f:
                f.write("STDOUT:\n")
                f.write(result.stdout)
                f.write("\nSTDERR:\n")
                f.write(result.stderr)
            
            print(f"Raw MLST output: {result.stdout.strip()}")
            
            # Parse the CSV output (it's comma-separated!)
            mlst_results = self.parse_mlst_csv(result.stdout, input_file.name)
            
            # Add lineage information
            mlst_results.update(self.get_lineage_info(mlst_results.get('st', 'ND')))
            
            # Generate only 3 output files
            self.generate_output_files(mlst_results, sample_output_dir)
            
            print(f"✅ Completed: {input_file.name} -> ST{mlst_results.get('st', 'ND')}")
            return mlst_results
            
        except subprocess.CalledProcessError as e:
            print(f"❌ MLST failed for {input_file.name}")
            error_result = self.get_fallback_results(input_file.name)
            self.generate_output_files(error_result, sample_output_dir)
            return error_result

    def parse_mlst_csv(self, stdout: str, sample_name: str) -> Dict:
        """Parse MLST CSV output - it's comma-separated!"""
        print(f"Parsing CSV output for {sample_name}")
        
        lines = stdout.strip().split('\n')
        if not lines:
            return self.get_empty_results(sample_name)
        
        # Find the result line (usually the last line with data)
        result_line = None
        for line in reversed(lines):
            if line.strip() and ',' in line and not line.startswith('['):
                result_line = line.strip()
                break
        
        if not result_line:
            return self.get_empty_results(sample_name)
        
        print(f"CSV result line: {result_line}")
        
        # Split by COMMA, not tab!
        parts = result_line.split(',')
        print(f"CSV parts: {parts}")
        
        if len(parts) < 3:
            return self.get_empty_results(sample_name)
        
        # Extract components - format: filename,scheme,ST,allele1,allele2,...
        filename = parts[0]
        scheme = parts[1]
        st = parts[2]
        
        # Extract alleles from remaining parts
        alleles = {}
        allele_parts = []
        
        for i in range(3, len(parts)):
            allele_str = parts[i]
            if '(' in allele_str and ')' in allele_str:
                # Format: arcC(1)
                gene = allele_str.split('(')[0]
                allele = allele_str.split('(')[1].rstrip(')')
                alleles[gene] = allele
                allele_parts.append(f"{gene}({allele})")
        
        allele_profile = '-'.join(allele_parts) if allele_parts else ""
        
        return {
            "sample": sample_name,
            "st": st,
            "scheme": scheme,
            "alleles": alleles,
            "allele_profile": allele_profile,
            "confidence": "HIGH" if st and st != '-' and st != 'ND' else "LOW"
        }

    def get_lineage_info(self, st: str) -> Dict:
        """Get comprehensive lineage information based on ST"""
        lineage_db = {
            '1': {
                "clonal_complex": "CC1",
                "classification": "Community-associated MRSA",
                "geographic_distribution": "Global",
                "clinical_significance": "Includes USA400 clone (MW2), often PVL-positive, associated with community infections",
                "common_virulence": ["Panton-Valentine Leukocidin", "Enterotoxins SEH/SEK", "Hemolysins"],
                "outbreak_potential": "HIGH",
                "typical_sccmec": ["IVa", "IVc"],
                "typical_spa": ["t128", "t127", "t174"],
                "resistance_profile": ["Often community-associated resistance patterns"]
            },
            '5': {
                "clonal_complex": "CC5",
                "classification": "Healthcare-associated MRSA",
                "geographic_distribution": "Global",
                "clinical_significance": "Major healthcare-associated lineage, includes New York/Japan clone, often multidrug-resistant with extensive virulence arsenal",
                "common_virulence": ["Enterotoxins (SEC, SEL, SEI, SEM, SEO)", "Immune evasion cluster (scn, chp, sak)", "Hemolysins", "Proteases"],
                "outbreak_potential": "VERY HIGH",
                "typical_sccmec": ["II", "IV"],
                "typical_spa": ["t002", "t242", "t548", "t688"],
                "resistance_profile": ["Methicillin", "Multiple aminoglycosides", "Macrolides", "Tetracyclines"]
            },
            '8': {
                "clonal_complex": "CC8", 
                "classification": "Community and Healthcare-associated MRSA",
                "geographic_distribution": "Global",
                "clinical_significance": "Includes pandemic USA300 clone, highly virulent and transmissible, often PVL-positive with enhanced fitness",
                "common_virulence": ["Panton-Valentine Leukocidin", "Enterotoxins", "Phenol-soluble modulins", "Arginine catabolic mobile element"],
                "outbreak_potential": "VERY HIGH",
                "typical_sccmec": ["IVa", "IV"],
                "typical_spa": ["t008", "t064", "t121", "t024"],
                "resistance_profile": ["Methicillin", "Often fluoroquinolone resistance", "Community-associated resistance patterns"]
            },
            '22': {
                "clonal_complex": "CC22",
                "classification": "Healthcare-associated MRSA",
                "geographic_distribution": "Europe, Middle East, Global", 
                "clinical_significance": "Epidemic MRSA-15 (EMRSA-15), major healthcare clone with high transmission in hospital settings",
                "common_virulence": ["Enterotoxins", "Immune evasion cluster", "Hemolysins"],
                "outbreak_potential": "VERY HIGH",
                "typical_sccmec": ["IV", "IVh"],
                "typical_spa": ["t032", "t022", "t005", "t852"],
                "resistance_profile": ["Methicillin", "Multiple drug classes", "Often gentamicin resistant"]
            },
            '30': {
                "clonal_complex": "CC30",
                "classification": "Healthcare and Community-associated",
                "geographic_distribution": "Global",
                "clinical_significance": "Includes EMRSA-16 and Southwest Pacific clone, often PVL-positive, associated with both hospital and community settings",
                "common_virulence": ["Panton-Valentine Leukocidin", "Enterotoxins", "Hemolysins"],
                "outbreak_potential": "HIGH",
                "typical_sccmec": ["IV", "IVc"],
                "typical_spa": ["t019", "t021", "t318", "t018"],
                "resistance_profile": ["Methicillin", "Variable resistance patterns"]
            },
            '45': {
                "clonal_complex": "CC45",
                "classification": "Community and Healthcare-associated",
                "geographic_distribution": "Europe, North America",
                "clinical_significance": "Includes USA600 clone, often associated with both community and healthcare settings",
                "common_virulence": ["Enterotoxins", "Hemolysins"],
                "outbreak_potential": "MEDIUM",
                "typical_sccmec": ["IV", "II"],
                "typical_spa": ["t015", "t026", "t038"],
                "resistance_profile": ["Methicillin", "Variable resistance"]
            },
            '80': {
                "clonal_complex": "CC80",
                "classification": "Community-associated MRSA",
                "geographic_distribution": "Europe, Middle East, Asia",
                "clinical_significance": "European CA-MRSA clone, often PVL-positive, associated with skin and soft tissue infections",
                "common_virulence": ["Panton-Valentine Leukocidin", "Enterotoxins"],
                "outbreak_potential": "HIGH",
                "typical_sccmec": ["IV", "IVc"],
                "typical_spa": ["t044", "t131", "t186"],
                "resistance_profile": ["Methicillin", "Often fusidic acid resistant"]
            },
            '93': {
                "clonal_complex": "CC93",
                "classification": "Community-associated MRSA", 
                "geographic_distribution": "Australia, Asia",
                "clinical_significance": "Queensland clone, often associated with community infections in Australia and Asia",
                "common_virulence": ["Enterotoxins", "Hemolysins"],
                "outbreak_potential": "MEDIUM",
                "typical_sccmec": ["IV", "V"],
                "typical_spa": ["t202", "t1340"],
                "resistance_profile": ["Methicillin", "Variable resistance"]
            },
            '121': {
                "clonal_complex": "CC121",
                "classification": "Community-associated MRSA",
                "geographic_distribution": "Europe, Global",
                "clinical_significance": "Often associated with exotoxin production and skin infections",
                "common_virulence": ["Exfoliative toxins", "Enterotoxins"],
                "outbreak_potential": "MEDIUM",
                "typical_sccmec": ["IV", "V"],
                "typical_spa": ["t159", "t314", "t645"],
                "resistance_profile": ["Methicillin", "Often fusidic acid resistant"]
            },
            '152': {
                "clonal_complex": "CC152",
                "classification": "Community-associated MRSA",
                "geographic_distribution": "Europe, Middle East",
                "clinical_significance": "Often PVL-positive, associated with community-acquired infections",
                "common_virulence": ["Panton-Valentine Leukocidin", "Enterotoxins"],
                "outbreak_potential": "MEDIUM",
                "typical_sccmec": ["IV", "V"],
                "typical_spa": ["t355", "t657"],
                "resistance_profile": ["Methicillin", "Variable resistance"]
            },
            '398': {
                "clonal_complex": "CC398",
                "classification": "Livestock-associated MRSA", 
                "geographic_distribution": "Global",
                "clinical_significance": "Zoonotic transmission from livestock, emerging public health concern with human infections",
                "common_virulence": ["Limited virulence arsenal", "Adapted to animal hosts"],
                "outbreak_potential": "HIGH (in agricultural settings)",
                "typical_sccmec": ["IV", "V"],
                "typical_spa": ["t011", "t034", "t108"],
                "resistance_profile": ["Methicillin", "Tetracycline", "Often multidrug-resistant"]
            }
        }
        
        return lineage_db.get(st, {
            "clonal_complex": f"CC{st}" if st != 'ND' and st != '-' else 'Unknown',
            "classification": "Novel or uncommon lineage",
            "geographic_distribution": "Further characterization needed",
            "clinical_significance": "This sequence type requires additional epidemiological analysis to determine clinical significance",
            "common_virulence": ["Virulence profile not well characterized"],
            "outbreak_potential": "UNKNOWN",
            "typical_sccmec": ["Variable"],
            "typical_spa": ["Variable"],
            "resistance_profile": ["Further testing required"]
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

    def generate_output_files(self, mlst_results: Dict, output_dir: Path):
        """Generate only 3 output files: HTML, TXT, and TSV"""
        # 1. Beautiful HTML Report
        self.generate_html_report(mlst_results, output_dir)
        
        # 2. Detailed Text Report
        self.generate_text_report(mlst_results, output_dir)
        
        # 3. Simple TSV Report
        self.generate_tsv_report(mlst_results, output_dir)

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
        
        # Lineage information
        report += f"""
LINEAGE ANALYSIS:
-----------------
Clonal Complex: {mlst_results.get('clonal_complex', 'Unknown')}
Classification: {mlst_results.get('classification', 'Unknown')}
Geographic Distribution: {mlst_results.get('geographic_distribution', 'Unknown')}
Clinical Significance: {mlst_results.get('clinical_significance', 'Unknown')}
Outbreak Potential: {mlst_results.get('outbreak_potential', 'UNKNOWN')}

Typical SCCmec Types: {', '.join(mlst_results.get('typical_sccmec', ['Unknown']))}
Typical spa Types: {', '.join(mlst_results.get('typical_spa', ['Unknown']))}
Resistance Profile: {', '.join(mlst_results.get('resistance_profile', ['Unknown']))}

Common Virulence Factors:
"""
        for virulence in mlst_results.get('common_virulence', []):
            report += f"- {virulence}\n"
        
        with open(output_dir / "mlst_report.txt", 'w') as f:
            f.write(report)

    def generate_tsv_report(self, mlst_results: Dict, output_dir: Path):
        """Generate simple TSV report"""
        tsv_content = f"Sample\tST\tScheme\tClonal_Complex\tClassification\tOutbreak_Potential\tTypical_SCCmec\tTypical_spa\tResistance_Profile\n"
        tsv_content += f"{mlst_results['sample']}\t{mlst_results['st']}\t{mlst_results['scheme']}\t{mlst_results.get('clonal_complex', 'Unknown')}\t{mlst_results.get('classification', 'Unknown')}\t{mlst_results.get('outbreak_potential', 'UNKNOWN')}\t{','.join(mlst_results.get('typical_sccmec', ['Unknown']))}\t{','.join(mlst_results.get('typical_spa', ['Unknown']))}\t{','.join(mlst_results.get('resistance_profile', ['Unknown']))}\n"
        
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
        
        .virulence-grid {{
            display: grid;
            grid-template-columns: repeat(auto-fill, minmax(200px, 1fr));
            gap: 10px;
            margin-top: 15px;
        }}
        
        .virulence-card {{
            background: #e0f2fe;
            color: #0369a1;
            padding: 12px;
            border-radius: 6px;
            text-align: center;
            font-weight: bold;
            border-left: 4px solid #0ea5e9;
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
        
        .profile-box {{
            background: #f8fafc;
            padding: 15px;
            border-radius: 8px;
            margin: 15px 0;
            border-left: 4px solid #3b82f6;
        }}
        
        @media (max-width: 768px) {{
            .ascii-art {{
                font-size: 6px;
            }}
            .allele-grid {{
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
                    <div class="metric-value">{mlst_results['sample']}</div>
                </div>
                <div class="metric-card">
                    <div class="metric-label">Analysis Date</div>
                    <div class="metric-value">{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}</div>
                </div>
                <div class="metric-card">
                    <div class="metric-label">MLST Scheme</div>
                    <div class="metric-value">{mlst_results['scheme'].title()}</div>
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
            <div class="profile-box">
                <code style="font-size: 16px; color: #1e40af; font-weight: bold;">{mlst_results['allele_profile']}</code>
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
            <div class="profile-box">
                <p style="margin: 0; line-height: 1.6; color: #374151;">{mlst_results.get('clinical_significance', 'Further analysis required.')}</p>
            </div>
            
            <div style="display: grid; grid-template-columns: repeat(auto-fit, minmax(300px, 1fr)); gap: 20px; margin-top: 20px;">
                <div>
                    <h3>Typical Characteristics</h3>
                    <p><strong>SCCmec Types:</strong> {', '.join(mlst_results.get('typical_sccmec', ['Unknown']))}</p>
                    <p><strong>spa Types:</strong> {', '.join(mlst_results.get('typical_spa', ['Unknown']))}</p>
                    <p><strong>Resistance Profile:</strong> {', '.join(mlst_results.get('resistance_profile', ['Unknown']))}</p>
                </div>
                <div>
                    <h3>Common Virulence Factors</h3>
                    <div class="virulence-grid">
'''
        
        for virulence in mlst_results.get('common_virulence', []):
            html_content += f'''                        <div class="virulence-card">{virulence}</div>
'''
        
        html_content += '''                    </div>
                </div>
            </div>
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

    def create_mlst_summary(self, all_results: Dict[str, Dict], output_dir: Path):
        """Create comprehensive MLST summary files for all samples"""
        print("📊 Creating MLST summary files...")
        
        # Create TSV summary
        self.create_mlst_tsv_summary(all_results, output_dir)
        
        # Create HTML summary
        self.create_mlst_html_summary(all_results, output_dir)
        
        print("✅ MLST summary files created successfully!")

    def create_mlst_tsv_summary(self, all_results: Dict[str, Dict], output_dir: Path):
        """Create TSV summary file with all samples"""
        summary_file = output_dir / "mlst_summary.tsv"
        
        with open(summary_file, 'w') as f:
            # Write header
            f.write("Sample\tST\tAllele_Profile\t")
            
            # Get all unique gene names from all samples
            all_genes = set()
            for result in all_results.values():
                all_genes.update(result['alleles'].keys())
            
            # Write gene headers
            for gene in sorted(all_genes):
                f.write(f"{gene}\t")
            f.write("\n")
            
            # Write data for each sample
            for sample_name, result in all_results.items():
                f.write(f"{sample_name}\t{result['st']}\t{result['allele_profile']}\t")
                
                # Write allele values for each gene
                for gene in sorted(all_genes):
                    allele = result['alleles'].get(gene, '')
                    f.write(f"{allele}\t")
                f.write("\n")
        
        print(f"📄 TSV summary created: {summary_file}")

    def create_mlst_html_summary(self, all_results: Dict[str, Dict], output_dir: Path):
        """Create beautiful HTML summary with MLST styling"""
        summary_file = output_dir / "mlst_summary.html"
        
        # Get all unique gene names
        all_genes = set()
        for result in all_results.values():
            all_genes.update(result['alleles'].keys())
        sorted_genes = sorted(all_genes)
        
        html_content = f'''<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>STAPHSCOPE - MLST Summary Report</title>
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
        
        .st-cell {{
            font-weight: bold;
            color: #1e40af;
        }}
        
        .allele-cell {{
            font-family: 'Courier New', monospace;
            background-color: #f0f9ff;
            color: #0369a1;
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
            <h2>📊 MLST Summary - All Samples</h2>
            
            <div class="stats-grid">
                <div class="stat-card">
                    <div class="stat-value">{len(all_results)}</div>
                    <div class="stat-label">SAMPLES PROCESSED</div>
                </div>
                <div class="stat-card">
                    <div class="stat-value">{len(set(result['st'] for result in all_results.values()))}</div>
                    <div class="stat-label">UNIQUE STs</div>
                </div>
                <div class="stat-card">
                    <div class="stat-value">{len(sorted_genes)}</div>
                    <div class="stat-label">MLST GENES</div>
                </div>
            </div>
            
            <div style="overflow-x: auto;">
                <table class="summary-table">
                    <thead>
                        <tr>
                            <th>Sample</th>
                            <th>ST</th>
                            <th>Allele Profile</th>
'''
        
        # Add gene headers
        for gene in sorted_genes:
            html_content += f'                            <th>{gene}</th>\n'
        
        html_content += '''                        </tr>
                    </thead>
                    <tbody>
'''
        
        # Add data rows
        for sample_name, result in all_results.items():
            html_content += f'''                        <tr>
                            <td><strong>{sample_name}</strong></td>
                            <td class="st-cell">ST{result['st']}</td>
                            <td class="allele-cell">{result['allele_profile']}</td>
'''
            
            # Add allele values for each gene
            for gene in sorted_genes:
                allele = result['alleles'].get(gene, '')
                html_content += f'                            <td>{allele}</td>\n'
            
            html_content += '                        </tr>\n'
        
        html_content += '''                    </tbody>
                </table>
            </div>
        </div>
        
        <div class="footer">
            <p><strong>STAPHSCOPE</strong> - MLST Summary Report</p>
            <p class="timestamp">Generated: ''' + datetime.now().strftime("%Y-%m-%d %H:%M:%S") + '''</p>
            <p>Github: bbeckley-hub</p>
        </div>
    </div>
</body>
</html>'''
        
        with open(summary_file, 'w', encoding='utf-8') as f:
            f.write(html_content)
        
        print(f"🌐 HTML summary created: {summary_file}")

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
        
        # Create summary files after processing all samples
        self.create_mlst_summary(results, output_dir)
        
        return results

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
    
    analyzer = ModularMLSTAnalyzer(
        database_dir=Path(args.database_dir),
        script_dir=Path(args.script_dir)
    )
    
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    if args.batch:
        results = analyzer.run_mlst_batch(args.input, output_dir, args.scheme)
        print(f"🎉 Batch MLST completed! Processed {len(results)} samples")
    else:
        input_file = Path(args.input)
        if input_file.exists():
            result = analyzer.run_mlst_single(input_file, output_dir, args.scheme)
            print(f"🎉 MLST completed for {input_file.name}: ST{result.get('st', 'ND')}")
        else:
            print(f"❌ Input file not found: {args.input}")

if __name__ == "__main__":
    main()

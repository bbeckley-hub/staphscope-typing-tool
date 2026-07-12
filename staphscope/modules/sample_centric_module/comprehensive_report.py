#!/usr/bin/env python3
"""
StaphScope Comprehensive Report Generator
Creates unified report from MLST, spa, and SCCmec results
Author: Brown Beckley <brownbeckley94@gmail.com>
Date: 2025
Affiliation: University of Ghana Medical School-Department of Medical Biochemsitry
Please kindly send a quick email for any technical issues or assistant
"""

import os
import sys
import json
import pandas as pd
from pathlib import Path
from datetime import datetime
import random

# Science quotes for rotation
SCIENCE_QUOTES = [
    {"text": "The important thing is not to stop questioning. Curiosity has its own reason for existing.", "author": "Albert Einstein"},
    {"text": "Science is not only a disciple of reason but also one of romance and passion.", "author": "Stephen Hawking"},
    {"text": "Somewhere, something incredible is waiting to be known.", "author": "Carl Sagan"},
    {"text": "The good thing about science is that it's true whether or not you believe in it.", "author": "Neil deGrasse Tyson"},
    {"text": "In science, there are no shortcuts to truth.", "author": "Karl Popper"},
    {"text": "Science knows no country, because knowledge belongs to humanity.", "author": "Louis Pasteur"},
    {"text": "The science of today is the technology of tomorrow.", "author": "Edward Teller"},
    {"text": "Nothing in life is to be feared, it is only to be understood.", "author": "Marie Curie"},
    {"text": "Research is what I'm doing when I don't know what I'm doing.", "author": "Wernher von Braun"},
    {"text": "The universe is not required to be in perfect harmony with human ambition.", "author": "Carl Sagan"},
    {"text": "Thank you for using STAPHSCOPE and don't forget to share with us any feature suggestions in mind.", "author": "Brown Beckley"},
    {"text": "Staphscope represents the convergence of genomic surveillance and clinical diagnostics, transforming raw sequences into actionable insights for infection control.", "author": "Brown Beckley"},
    {"text": "In the battle against antimicrobial resistance, tools like Staphscope are our eyes and ears, revealing the genetic blueprints of resistant pathogens.", "author": "Brown Beckley"},
    {"text": "Staphscope isn't just a tool; it's a comprehensive system that bridges the gap between sequencing data and public health action.", "author": "Brown Beckley"},
    {"text": "Through Staphscope, we turn the complexity of bacterial genomes into clear, interpretable reports, empowering clinicians and researchers alike.", "author": "Brown Beckley"},
    {"text": "Staphscope is a testament to the power of bioinformatics in the modern era, making advanced pathogen typing accessible to all.", "author": "Brown Beckley"}
]

def normalize_sample_name(sample_name):
    """Remove .fna extension if present for consistent matching"""
    if not sample_name:
        return ""
    if sample_name.endswith('.fna'):
        return sample_name[:-4]
    if sample_name.endswith('.fasta'):
        return sample_name[:-6]
    if sample_name.endswith('.fa'):
        return sample_name[:-3]
    if sample_name.endswith('.fn'):
        return sample_name[:-3]
    return sample_name

def load_mlst_data(mlst_file):
    """Load MLST data from TSV file"""
    print(f"📊 Loading MLST data from: {mlst_file}")
    
    try:
        df = pd.read_csv(mlst_file, sep='\t')
        mlst_data = {}
        
        for _, row in df.iterrows():
            sample = str(row['Sample']).strip()
            st = str(row['ST']).strip()
            
            # Normalize sample name
            norm_sample = normalize_sample_name(sample)
            
            if st and st not in ['-', 'ND', '', 'UNKNOWN', 'Not Assigned']:
                mlst_data[norm_sample] = f"ST{st}"
            else:
                mlst_data[norm_sample] = "Not typed"
        
        print(f"   ✓ Loaded {len(mlst_data)} samples")
        return mlst_data
        
    except Exception as e:
        print(f"   ✗ Error loading MLST data: {e}")
        return {}

def load_spa_data(spa_file):
    """Load spa typing data from TSV file"""
    print(f"📊 Loading spa data from: {spa_file}")
    
    try:
        df = pd.read_csv(spa_file, sep='\t')
        spa_data = {}
        
        for _, row in df.iterrows():
            sample = str(row['Sample']).strip()
            spa_type = str(row['spa_Type']).strip()
            
            # Normalize sample name
            norm_sample = normalize_sample_name(sample)
            
            if spa_type and spa_type not in ['-', 'ND', '', 'UNKNOWN', 'Not typed']:
                spa_data[norm_sample] = spa_type
            else:
                spa_data[norm_sample] = "Not typed"
        
        print(f"   ✓ Loaded {len(spa_data)} samples")
        return spa_data
        
    except Exception as e:
        print(f"   ✗ Error loading spa data: {e}")
        return {}

def load_sccmec_data(sccmec_file):
    """Load SCCmec data from TSV file"""
    print(f"📊 Loading SCCmec data from: {sccmec_file}")
    
    try:
        df = pd.read_csv(sccmec_file, sep='\t')
        sccmec_data = {}
        
        for _, row in df.iterrows():
            sample = str(row['Sample_Name']).strip()
            sccmec_type = str(row['SCCmec_Type']).strip()
            mrsa_status = str(row['MRSA_Status']).strip()
            
            # Normalize sample name (SCCmec file doesn't have .fna extension)
            norm_sample = normalize_sample_name(sample)
            
            # Handle SCCmec type
            if sccmec_type and sccmec_type not in ['-', 'ND', '', 'UNKNOWN', 'Not typed']:
                final_sccmec = sccmec_type
            else:
                final_sccmec = "Not typed"
            
            # Handle MRSA status
            if mrsa_status and mrsa_status not in ['-', 'ND', '', 'UNKNOWN']:
                final_mrsa = mrsa_status
            else:
                final_mrsa = "Unknown"
            
            sccmec_data[norm_sample] = {
                'sccmec_type': final_sccmec,
                'mrsa_status': final_mrsa
            }
        
        print(f"   ✓ Loaded {len(sccmec_data)} samples")
        return sccmec_data
        
    except Exception as e:
        print(f"   ✗ Error loading SCCmec data: {e}")
        return {}

def combine_data(mlst_data, spa_data, sccmec_data):
    """Combine data from all three sources"""
    print("\n🔗 Combining data from all sources...")
    
    # Get all unique samples from all sources
    all_samples = set()
    all_samples.update(mlst_data.keys())
    all_samples.update(spa_data.keys())
    all_samples.update(sccmec_data.keys())
    
    combined = []
    
    for sample in sorted(all_samples):
        # Get data from each source
        mlst = mlst_data.get(sample, "Not typed")
        spa_type = spa_data.get(sample, "Not typed")
        
        sccmec_info = sccmec_data.get(sample, {})
        sccmec_type = sccmec_info.get('sccmec_type', "Not typed")
        mrsa_status = sccmec_info.get('mrsa_status', "Unknown")
        
        combined.append({
            'sample': sample,
            'mlst': mlst,
            'spa_type': spa_type,
            'sccmec_type': sccmec_type,
            'mrsa_status': mrsa_status
        })
    
    print(f"   ✓ Combined data for {len(combined)} samples")
    return combined

def generate_tsv_report(data, output_file):
    """Generate TSV report"""
    print(f"\n📄 Generating TSV report: {output_file}")
    
    try:
        df = pd.DataFrame(data)
        # Reorder columns
        df = df[['sample', 'mlst', 'spa_type', 'sccmec_type', 'mrsa_status']]
        df.to_csv(output_file, sep='\t', index=False)
        print(f"   ✓ TSV report generated successfully")
        return True
    except Exception as e:
        print(f"   ✗ Error generating TSV: {e}")
        return False

def generate_json_report(data, output_file):
    """Generate JSON report"""
    print(f"📄 Generating JSON report: {output_file}")
    
    try:
        # Calculate statistics
        total_samples = len(data)
        mrsa_count = sum(1 for d in data if d['mrsa_status'] == 'MRSA')
        mssa_count = sum(1 for d in data if d['mrsa_status'] == 'MSSA')
        
        report_data = {
            'metadata': {
                'generated_date': datetime.now().isoformat(),
                'total_samples': total_samples,
                'report_type': 'StaphScope Comprehensive Report',
                'version': '1.0',
                'author': 'Brown Beckley',
                'email': 'brownbeckley94@gmail.com',
                'affiliation': 'University of Ghana Medical School - Department of Medical Biochemistry'
            },
            'statistics': {
                'mrsa_count': mrsa_count,
                'mrsa_percentage': (mrsa_count / total_samples * 100) if total_samples > 0 else 0,
                'mssa_count': mssa_count,
                'mssa_percentage': (mssa_count / total_samples * 100) if total_samples > 0 else 0,
                'samples_with_mlst': sum(1 for d in data if d['mlst'] != 'Not typed'),
                'samples_with_spa': sum(1 for d in data if d['spa_type'] != 'Not typed'),
                'samples_with_sccmec': sum(1 for d in data if d['sccmec_type'] != 'Not typed')
            },
            'samples': data
        }
        
        with open(output_file, 'w', encoding='utf-8') as f:
            json.dump(report_data, f, indent=2, ensure_ascii=False)
        
        print(f"   ✓ JSON report generated successfully")
        return True
    except Exception as e:
        print(f"   ✗ Error generating JSON: {e}")
        return False

def generate_html_report(data, output_file):
    """Generate beautiful HTML report with ASCII art and science quotes"""
    print(f"📄 Generating HTML report: {output_file}")
    
    try:
        # Get random science quote
        random_quote = random.choice(SCIENCE_QUOTES)
        
        # Calculate statistics
        total_samples = len(data)
        mrsa_count = sum(1 for d in data if d['mrsa_status'] == 'MRSA')
        mssa_count = sum(1 for d in data if d['mrsa_status'] == 'MSSA')
        unknown_count = total_samples - mrsa_count - mssa_count
        
        # Generate table rows
        table_rows = ""
        for row in data:
            # Color code MRSA status
            if row['mrsa_status'] == 'MRSA':
                mrsa_class = 'mrsa-positive'
            elif row['mrsa_status'] == 'MSSA':
                mrsa_class = 'mrsa-negative'
            else:
                mrsa_class = 'unknown-status'
            
            table_rows += f"""
                <tr>
                    <td>{row['sample']}</td>
                    <td>{row['mlst']}</td>
                    <td>{row['spa_type']}</td>
                    <td>{row['sccmec_type']}</td>
                    <td class="{mrsa_class}">{row['mrsa_status']}</td>
                </tr>"""
        
        html_content = f'''<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>StaphScope Comprehensive Report</title>
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
        
        .quote-text {{
            font-size: 18px;
            font-style: italic;
            margin-bottom: 10px;
            color: #ffffff;
        }}
        
        .quote-author {{
            font-size: 14px;
            color: #fbbf24;
            font-weight: bold;
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
        
        table {{
            width: 100%;
            border-collapse: collapse;
            margin-top: 20px;
            font-size: 14px;
        }}
        
        th {{
            background: linear-gradient(135deg, #3b82f6 0%, #1e40af 100%);
            color: white;
            padding: 12px;
            text-align: left;
            position: sticky;
            top: 0;
        }}
        
        td {{
            padding: 10px;
            border-bottom: 1px solid #e5e7eb;
        }}
        
        tr:nth-child(even) {{
            background-color: #f8fafc;
        }}
        
        tr:hover {{
            background-color: #e0f2fe;
        }}
        
        .mrsa-positive {{
            color: #dc2626;
            font-weight: bold;
        }}
        
        .mrsa-negative {{
            color: #10b981;
            font-weight: bold;
        }}
        
        .unknown-status {{
            color: #6b7280;
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
        
        @media (max-width: 768px) {{
            .ascii-art {{
                font-size: 6px;
            }}
            table {{
                font-size: 12px;
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
            
            <div class="quote-container" id="quoteContainer">
                <div class="quote-text" id="quoteText">"{random_quote['text']}"</div>
                <div class="quote-author" id="quoteAuthor">— {random_quote['author']}</div>
            </div>
        </div>
        
        <div class="report-section">
            <h2>📊 Comprehensive Typing Report</h2>
            
            <div class="stats-grid">
                <div class="stat-card">
                    <div class="stat-value">{total_samples}</div>
                    <div class="stat-label">TOTAL SAMPLES</div>
                </div>
                <div class="stat-card">
                    <div class="stat-value">{mrsa_count}</div>
                    <div class="stat-label">MRSA</div>
                </div>
                <div class="stat-card">
                    <div class="stat-value">{mssa_count}</div>
                    <div class="stat-label">MSSA</div>
                </div>
                <div class="stat-card">
                    <div class="stat-value">{unknown_count}</div>
                    <div class="stat-label">UNKNOWN</div>
                </div>
            </div>
            
            <p><strong>Generated:</strong> {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}</p>
            <p><strong>Report includes:</strong> Sample, MLST, spa type, SCCmec type, MRSA/MSSA classification</p>
        </div>
        
        <div class="report-section">
            <h2>🧬 Sample Typing Data ({total_samples} samples)</h2>
            
            <div style="overflow-x: auto;">
                <table>
                    <thead>
                        <tr>
                            <th>Sample</th>
                            <th>MLST</th>
                            <th>spa Type</th>
                            <th>SCCmec Type</th>
                            <th>MRSA/MSSA Status</th>
                        </tr>
                    </thead>
                    <tbody>
                        {table_rows}
                    </tbody>
                </table>
            </div>
        </div>
        
        <div class="footer">
            <p><strong>STAPHSCOPE</strong> - Comprehensive Analysis Report</p>
            <p class="timestamp">Generated: {datetime.now().strftime("%Y-%m-%d %H:%M:%S")}</p>
            <p><strong>Author:</strong> Brown Beckley | <strong>GitHub:</strong> bbeckley-hub</p>
            <p><strong>Email:</strong> brownbeckley94@gmail.com</p>
            <p><strong>Affiliation:</strong> University of Ghana Medical School - Department of Medical Biochemistry</p>
            <p><em>Report generated from MLST, spa, and SCCmec typing results.</em></p>
            <p><em>⭐ Star us on GitHub if you find this tool useful!</em></p>
            <p><em>Transforming fragmented genomic data into coherent biological narratives 🧬✨</em></p>
        </div>
    </div>

    <script>
        // Science quotes for rotation
        const quotes = {json.dumps(SCIENCE_QUOTES)};
        
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
        
        // Rotate quotes every 10 seconds
        setInterval(displayQuote, 10000);
    </script>
</body>
</html>'''
        
        with open(output_file, 'w', encoding='utf-8') as f:
            f.write(html_content)
        
        print(f"   ✓ HTML report generated successfully")
        return True
    except Exception as e:
        print(f"   ✗ Error generating HTML: {e}")
        return False

def main():
    """Main function"""
    print("\n" + "="*60)
    print("STAPHSCOPE COMPREHENSIVE REPORT GENERATOR")
    print("="*60)
    
    # Check if we're in the right directory
    current_dir = Path.cwd()
    
    # Look for the required files
    mlst_file = current_dir / "mlst_summary.tsv"
    spa_file = current_dir / "spa_summary.tsv"
    sccmec_file = current_dir / "staphscope_summary.tsv"
    
    # Check if files exist
    if not mlst_file.exists():
        print(f"❌ Error: {mlst_file} not found!")
        print("   Please run this script from the directory containing the summary files.")
        sys.exit(1)
    
    if not spa_file.exists():
        print(f"❌ Error: {spa_file} not found!")
        sys.exit(1)
    
    if not sccmec_file.exists():
        print(f"❌ Error: {sccmec_file} not found!")
        sys.exit(1)
    
    print("✓ Found all required input files")
    
    # Load data
    mlst_data = load_mlst_data(mlst_file)
    spa_data = load_spa_data(spa_file)
    sccmec_data = load_sccmec_data(sccmec_file)
    
    if not mlst_data or not spa_data or not sccmec_data:
        print("❌ Failed to load one or more data files")
        sys.exit(1)
    
    # Combine data
    combined_data = combine_data(mlst_data, spa_data, sccmec_data)
    
    if not combined_data:
        print("❌ No data to generate report")
        sys.exit(1)
    
    # Generate reports
    base_name = "staphscope_comprehensive_report"
    
    tsv_success = generate_tsv_report(combined_data, current_dir / f"{base_name}.tsv")
    json_success = generate_json_report(combined_data, current_dir / f"{base_name}.json")
    html_success = generate_html_report(combined_data, current_dir / f"{base_name}.html")
    
    print("\n" + "="*60)
    if tsv_success and json_success and html_success:
        print("✅ COMPREHENSIVE REPORT GENERATED SUCCESSFULLY!")
        print("="*60)
        print(f"📁 Reports saved in: {current_dir}")
        print(f"   • {base_name}.html (HTML report)")
        print(f"   • {base_name}.json (JSON data)")
        print(f"   • {base_name}.tsv (TSV data)")
        print(f"📊 Summary: {len(combined_data)} samples")
        
        # Count MRSA/MSSA
        mrsa_count = sum(1 for d in combined_data if d['mrsa_status'] == 'MRSA')
        mssa_count = sum(1 for d in combined_data if d['mrsa_status'] == 'MSSA')
        print(f"🦠 MRSA: {mrsa_count}, MSSA: {mssa_count}")
    else:
        print("❌ Some reports failed to generate")
    
    print("="*60)

if __name__ == "__main__":
    main()

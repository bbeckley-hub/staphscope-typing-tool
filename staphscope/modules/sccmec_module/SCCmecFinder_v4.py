#!/usr/bin/env python3
"""
SCCmecFinder_v4 - Python 3 Implementation with Enhanced JSON Reporting
Original Author: Hulya Kaya
Rewritten by: Beckley Brown <brownbeckley94@gmail.com>
Date: 2025-08-18
Updated: 2025-08-20 (Added enhanced JSON reporting)
Updated: 2025-08-21 (Enhanced HTML reporting with detailed results)
Updated: 2025-08-22 (Improved k-mer results display with colored table)
Send a quick mail for any issues or further explanations.
Affiliation: University of Ghana Medical School-Department of Medical Bioichemistry
"""

import os
import sys
import time
import glob
import argparse
import subprocess
import json
from pathlib import Path
from typing import Dict, List, Set, Tuple, Optional, Any

# SCCmec type definitions
SCCMEC_DEFINITIONS = {
    'SCCmec_type_I(1B)': {'ccrA1', 'ccrB1', 'mecA', 'dmecR1', 'IS1272'},
    'SCCmec_type_II(2A)': {'ccrA2', 'ccrB2', 'mecA', 'mecR1', 'mecI'},
    'SCCmec_type_III(3A)': {'ccrA3', 'ccrB3', 'mecA', 'mecR1', 'mecI'},
    'SCCmec_type_IV(2B)': {'ccrA2', 'ccrB2', 'mecA', 'dmecR1', 'IS1272'},
    'SCCmec_type_IV(2B&5)': {'ccrA2', 'ccrB2', 'ccrC11', 'mecA', 'dmecR1', 'IS1272'},
    'SCCmec_type_V(5C2)': {'ccrC11', 'mec-class-C2', 'mecA'},
    'SCCmec_type_V(5C2&5)': {'ccrC11', 'ccrC12', 'mec-class-C2', 'mecA'},
    'SCCmec_type_VI(4B)': {'ccrA4', 'ccrB4', 'mecA', 'dmecR1', 'IS1272'},
    'SCCmec_type_VII(5C1)': {'ccrC11', 'mec-class-C1', 'mecA'},
    'SCCmec_type_VIII(4A)': {'ccrA4', 'ccrB4', 'mecA', 'mecR1', 'mecI'},
    'SCCmec_type_IX(1C2)': {'ccrA1', 'ccrB1', 'mec-class-C2', 'mecA'},
    'SCCmec_type_X(7C1)': {'ccrA1', 'ccrB6', 'mec-class-C1', 'mecA'},
    'SCCmec_type_XI(8E)': {'ccrA1', 'ccrB3', 'mecC', 'mecR1', 'mecI'},
    'SCCmec_type_XII(9C2)': {'ccrC2', 'mec-class-C2', 'mecA'},
    'SCCmec_type_XIII(9A)': {'ccrC2', 'mecA', 'dmecR1', 'IS1272'}
}

SCCMEC_CLASSES = {
    'SCCmec_type_I(1B)': {'ccr class 1', 'mec class B'},
    'SCCmec_type_II(2A)': {'ccr class 2', 'mec class A'},
    'SCCmec_type_III(3A)': {'ccr class 3', 'mec class A'},
    'SCCmec_type_IV(2B)': {'ccr class 2', 'mec class B'},
    'SCCmec_type_IV(2B&5)': {'ccr class 5', 'ccr class 2', 'mec class B'},
    'SCCmec_type_V(5C2)': {'ccr class 5', 'mec class C2'},
    'SCCmec_type_V(5C2&5)': {'ccr class 5&5', 'mec class C2'},
    'SCCmec_type_VI(4B)': {'ccr class 4', 'mec class B'},
    'SCCmec_type_VII(5C1)': {'ccr class 5', 'mec class C1'},
    'SCCmec_type_VIII(4A)': {'ccr class 4', 'mec class A'},
    'SCCmec_type_IX(1C2)': {'ccr class 1', 'mec class C2'},
    'SCCmec_type_X(7C1)': {'ccr class 7', 'mec class C1'},
    'SCCmec_type_XI(8E)': {'ccr class 8', 'mec class E'},
    'SCCmec_type_XII(9C2)': {'ccr class 9', 'mec class C2'},
    'SCCmec_type_XIII(9A)': {'ccr class 9', 'mec class A'}
}

def perform_ccr_gene_complex_typing(ccrAB_genes: Set[str], ccrC_genes: Set[str]) -> List[str]:
    """Determine CCR gene complex classes based on detected genes"""
    classes = []
    
    if {"ccrA1", "ccrB1"}.issubset(ccrAB_genes):
        classes.append("ccr class 1")
    if {"ccrA2", "ccrB2"}.issubset(ccrAB_genes):
        classes.append("ccr class 2")
    if {"ccrC11", "ccrC12"}.issubset(ccrC_genes):
        classes.append("ccr class 5&5")
    elif "ccrC11" in ccrC_genes:
        classes.append("ccr class 5")
    if {"ccrA3", "ccrB3"}.issubset(ccrAB_genes):
        classes.append("ccr class 3")
    if {"ccrA4", "ccrB4"}.issubset(ccrAB_genes):
        classes.append("ccr class 4")
    if {"ccrA5", "ccrB3"}.issubset(ccrAB_genes):
        classes.append("ccr class 6")
    if {"ccrA1", "ccrB6"}.issubset(ccrAB_genes):
        classes.append("ccr class 7")
    if {"ccrA1", "ccrB3"}.issubset(ccrAB_genes):
        classes.append("ccr class 8")
    if "ccrC2" in ccrC_genes:
        classes.append("ccr class 9")
    
    return classes

def perform_mec_gene_complex_typing(mec_genes: Set[str]) -> List[str]:
    """Determine MEC gene complex classes based on detected genes"""
    classes = []
    
    if {"mecA", "mecR1", "mecI"}.issubset(mec_genes):
        classes.append("mec class A")
    if {"mecA", "dmecR1", "IS1272"}.issubset(mec_genes):
        classes.append("mec class B")
    if {"mecA", "mec-class-C1"}.issubset(mec_genes):
        classes.append("mec class C1")
    if {"mecA", "mec-class-C2"}.issubset(mec_genes):
        classes.append("mec class C2")
    if {"mecC", "mecR1", "mecI"}.issubset(mec_genes):
        classes.append("mec class E")
    
    return classes

def perform_sccmec_typing(classes: Set[str]) -> List[str]:
    """Determine SCCmec types based on gene classes"""
    sccmec_types = []
    
    if {"ccr class 1", "mec class B"}.issubset(classes):
        sccmec_types.append("SCCmec_type_I(1B)")
    if {"ccr class 2", "mec class A"}.issubset(classes):
        sccmec_types.append("SCCmec_type_II(2A)")
    if {"ccr class 3", "mec class A"}.issubset(classes):
        sccmec_types.append("SCCmec_type_III(3A)")
    if {"ccr class 2", "ccr class 5", "mec class B"}.issubset(classes):
        sccmec_types.append("SCCmec_type_IV(2B&5)")
    elif {"ccr class 2", "mec class B"}.issubset(classes):
        sccmec_types.append("SCCmec_type_IV(2B)")
    if {"ccr class 5", "mec class C2"}.issubset(classes):
        sccmec_types.append("SCCmec_type_V(5C2)")
    if {"ccr class 5&5", "mec class C2"}.issubset(classes):
        sccmec_types.append("SCCmec_type_V(5C2&5)")
    if {"ccr class 4", "mec class B"}.issubset(classes):
        sccmec_types.append("SCCmec_type_VI(4B)")
    if {"ccr class 5", "mec class C1"}.issubset(classes):
        sccmec_types.append("SCCmec_type_VII(5C1)")
    if {"ccr class 4", "mec class A"}.issubset(classes):
        sccmec_types.append("SCCmec_type_VIII(4A)")
    if {"ccr class 1", "mec class C2"}.issubset(classes):
        sccmec_types.append("SCCmec_type_IX(1C2)")
    if {"ccr class 7", "mec class C1"}.issubset(classes):
        sccmec_types.append("SCCmec_type_X(7C1)")
    if {"ccr class 8", "mec class E"}.issubset(classes):
        sccmec_types.append("SCCmec_type_XI(8E)")
    if {"ccr class 9", "mec class C2"}.issubset(classes):
        sccmec_types.append("SCCmec_type_XII(9C2)")
    if {"ccr class 9", "mec class A"}.issubset(classes):
        sccmec_types.append("SCCmec_type_XIII(9A)")
    
    return sccmec_types

def run_command(cmd: List[str], cwd: Optional[Path] = None) -> None:
    """Run a system command with error handling"""
    try:
        subprocess.run(cmd, cwd=str(cwd) if cwd else None, check=True, 
                      stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    except subprocess.CalledProcessError as e:
        print(f"Error running command: {' '.join(cmd)}")
        print(f"Error message: {e.stderr.decode() if e.stderr else 'Unknown error'}")
        sys.exit(1)

def parse_detailed_gene_results(db_finder_file: Path) -> List[Dict]:
    """Parse detailed gene results from MyDbFinder output"""
    results = []
    if not db_finder_file.exists() or os.path.getsize(db_finder_file) == 0:
        return results
    
    with open(db_finder_file, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 4:
                gene_full = parts[0]
                identity = parts[1]
                length_info = parts[2]
                contig = parts[3]
                position = parts[4] if len(parts) > 4 else ""
                
                # Parse length info
                length_parts = length_info.split('/')
                query_length = length_parts[0] if len(length_parts) > 0 else ""
                hsp_length = length_parts[1] if len(length_parts) > 1 else ""
                
                results.append({
                    'gene': gene_full.split(':')[0] if ':' in gene_full else gene_full,
                    'gene_full': gene_full,
                    'identity': identity,
                    'query_length': query_length,
                    'hsp_length': hsp_length,
                    'contig': contig,
                    'position': position
                })
    
    return results

def parse_kmer_detailed_results(kmer_file: Path) -> List[Dict]:
    """Parse detailed k-mer results from MyKmerFinder output"""
    results = []
    if not kmer_file.exists() or os.path.getsize(kmer_file) == 0:
        return results
    
    with open(kmer_file, 'r') as f:
        next(f)  # Skip header
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 10:
                results.append({
                    'template': parts[0],
                    'score': parts[1],
                    'expected': parts[2],
                    'z': parts[3],
                    'p_value': parts[4],
                    'query_coverage': parts[5],
                    'template_coverage': parts[6],
                    'depth': parts[7],
                    'kmers_in_template': parts[8],
                    'description': parts[9] if len(parts) > 9 else ""
                })
    
    return results

def parse_raw_kmer_content_to_table(raw_content: str) -> Tuple[List[str], List[List[str]]]:
    """Parse raw k-mer content into header and rows for HTML table"""
    lines = raw_content.strip().split('\n')
    if not lines:
        return [], []
    
    # Extract header (remove # if present)
    header_line = lines[0]
    if header_line.startswith('#'):
        header_line = header_line[1:].strip()
    
    # Split header by tabs
    headers = [h.strip() for h in header_line.split('\t')]
    
    # Parse data rows
    rows = []
    for line in lines[1:]:
        if not line.strip():
            continue
        parts = line.strip().split('\t')
        rows.append([p.strip() for p in parts])
    
    return headers, rows

def get_top_kmer_predictions(kmer_results: List[Dict], has_mecA: bool) -> Dict[str, Any]:
    """
    Get top k-mer based predictions when gene-based fails.
    Returns top 2 hits with highest template coverage.
    """
    if not kmer_results or not has_mecA:
        return {
            'should_use_kmer': False,
            'top_hits': [],
            'contradiction': False,
            'message': ''
        }
    
    # Sort by template coverage (highest first)
    sorted_results = sorted(kmer_results, 
                          key=lambda x: float(x['template_coverage']), 
                          reverse=True)
    
    # Get top 2 hits
    top_hits = sorted_results[:2] if len(sorted_results) >= 2 else sorted_results[:1]
    
    # Check for contradiction if we have 2 hits
    contradiction = False
    message = ''
    
    if len(top_hits) >= 2:
        cov1 = float(top_hits[0]['template_coverage'])
        cov2 = float(top_hits[1]['template_coverage'])
        
        # If top 2 are close in coverage (within 10%), show contradiction
        if abs(cov1 - cov2) < 10.0 and cov1 >= 50.0:
            contradiction = True
            type1 = top_hits[0]['template'].split('|')[0]
            type2 = top_hits[1]['template'].split('|')[0]
            message = f"⚠️ K-mer based contradiction: Similar coverage between {type1} ({cov1:.1f}%) and {type2} ({cov2:.1f}%). Isolate could be either type."
        elif cov1 >= 50.0:
            # Good single prediction
            message = f"✅ K-mer based prediction: {top_hits[0]['template'].split('|')[0]} with {cov1:.1f}% coverage"
    
    return {
        'should_use_kmer': len(top_hits) > 0 and float(top_hits[0]['template_coverage']) >= 50.0,
        'top_hits': top_hits,
        'contradiction': contradiction,
        'message': message
    }

def generate_enhanced_json_report(sample_name: str, mrsa_gene: str, sccmec_types: List[str], 
                                 kmer_hits: List, kmer_prediction: Dict, ccrAB_genes: Set[str], 
                                 ccrC_genes: Set[str], mec_genes: Set[str], subtyping_genes: Set[str], 
                                 total_genes: Set[str], output_path: Path):
    """
    Generate reliable JSON report - only include what actually works
    """
    # Determine primary type - use the first gene-based prediction (most reliable)
    primary_type = sccmec_types[0] if sccmec_types else "NA"
    
    # Build MRSA evidence based on reliable gene detection
    mrsa_evidence = []
    if "mecA" in mec_genes or "mecA" in mrsa_gene:
        mrsa_evidence.append("mecA_detected")
    if "mecR1" in mec_genes:
        mrsa_evidence.append("mecR1_detected")
    
    # Simple confidence based on gene evidence
    if primary_type != "NA" and len(mrsa_evidence) >= 2:
        mrsa_confidence = "VERY_HIGH"
        sccmec_confidence = "HIGH"
    elif primary_type != "NA" and mrsa_evidence:
        mrsa_confidence = "HIGH" 
        sccmec_confidence = "HIGH"
    else:
        mrsa_confidence = "LOW"
        sccmec_confidence = "LOW"
    
    # Add SCCmec evidence if we have a type
    if primary_type != "NA":
        type_part = primary_type.split('_')[-1].split('(')[0]
        mrsa_evidence.append(f"sccmec_{type_part}_present")
    
    enhanced_report = {
        "sample": sample_name,
        "mrsa_status": {
            "classification": "CONFIRMED_MRSA" if mrsa_gene else "MSSA",
            "evidence": mrsa_evidence,
            "confidence": mrsa_confidence
        },
        "sccmec_typing": {
            "primary_type": primary_type,
            "all_predicted_types": sccmec_types,
            "confidence": sccmec_confidence,
            "typing_method": "gene_based"
        },
        "kmer_prediction": kmer_prediction,
        "gene_detection": {
            "mec_genes": [{"gene": gene, "detected": True} for gene in sorted(mec_genes)],
            "ccr_genes": [{"gene": gene, "detected": True} for gene in sorted(ccrAB_genes.union(ccrC_genes))],
            "regulator_genes": [{"gene": gene, "detected": True} for gene in sorted([g for g in total_genes if any(x in g for x in ['mecR', 'mecI', 'IS1272', 'dmec'])] if subtyping_genes else [])]
        },
        "quality_metrics": {
            "gene_completeness": assess_gene_completeness(ccrAB_genes, ccrC_genes, mec_genes),
            "contradictions": len(sccmec_types) > 1
        }
    }
    
    # Save JSON report
    json_file = output_path / "sccmec_enhanced_report.json"
    with open(json_file, 'w') as f:
        json.dump(enhanced_report, f, indent=2)
    
    return enhanced_report

def generate_staphscope_html_report(sample_name: str, mrsa_gene: str, sccmec_types: List[str],
                                   enhanced_report: Dict, ccrAB_genes: Set[str], ccrC_genes: Set[str],
                                   mec_genes: Set[str], kmer_hits: List, 
                                   detailed_gene_results: List[Dict], detailed_kmer_results: List[Dict],
                                   output_path: Path, args, kmer_prediction: Dict):
    """Generate beautiful StaphScope HTML report with ALL detailed results"""
    
    # Read the ACTUAL results_MyKmerFinder.txt file content
    full_kmer_content = ""
    kmer_file = output_path / "results_MyKmerFinder.txt"
    if kmer_file.exists():
        with open(kmer_file, 'r') as f:
            full_kmer_content = f.read()
    
    # Parse raw content for table display
    kmer_headers, kmer_rows = parse_raw_kmer_content_to_table(full_kmer_content)
    
    # Read additional files for comprehensive reporting
    html_kmer_content = ""
    html_kmer_file = output_path / "HTML_MyKmerFinder.txt"
    if html_kmer_file.exists():
        with open(html_kmer_file, 'r') as f:
            html_kmer_content = f.read()
    
    detailed_text_content = ""
    detailed_text_file = output_path / args.output_file
    if detailed_text_file.exists():
        with open(detailed_text_file, 'r') as f:
            detailed_text_content = f.read()
    
    html_content = '''<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>STAPHSCOPE - SCCmec Analysis Report</title>
    <style>
        * {
            margin: 0;
            padding: 0;
            box-sizing: border-box;
        }
        
        body {
            background: linear-gradient(135deg, #1e3c72 0%, #2a5298 50%, #7e22ce 100%);
            font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif;
            color: #ffffff;
            padding: 20px;
            min-height: 100vh;
        }
        
        .container {
            max-width: 1400px;
            margin: 0 auto;
        }
        
        .header {
            text-align: center;
            margin-bottom: 30px;
        }
        
        .ascii-container {
            background: rgba(0, 0, 0, 0.7);
            padding: 20px;
            border-radius: 15px;
            margin-bottom: 20px;
            box-shadow: 0 8px 32px rgba(0, 0, 0, 0.4);
            border: 2px solid rgba(0, 255, 0, 0.3);
        }
        
        .ascii-art {
            font-family: 'Courier New', monospace;
            font-size: 10px;
            line-height: 1.1;
            white-space: pre;
            color: #00ff00;
            text-shadow: 0 0 10px rgba(0, 255, 0, 0.5);
            overflow-x: auto;
        }
        
        .quote-container {
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
        }
        
        .quote-text {
            font-size: 18px;
            font-style: italic;
            margin-bottom: 10px;
            color: #ffffff;
        }
        
        .quote-author {
            font-size: 14px;
            color: #fbbf24;
            font-weight: bold;
        }
        
        .report-section {
            background: rgba(255, 255, 255, 0.95);
            color: #1f2937;
            padding: 25px;
            border-radius: 10px;
            margin-bottom: 20px;
            box-shadow: 0 4px 15px rgba(0, 0, 0, 0.2);
        }
        
        .report-section h2 {
            color: #1e3a8a;
            border-bottom: 3px solid #3b82f6;
            padding-bottom: 10px;
            margin-bottom: 20px;
            font-size: 24px;
        }
        
        .report-section h3 {
            color: #1e40af;
            margin-top: 20px;
            margin-bottom: 10px;
            font-size: 18px;
        }
        
        .status-badge {
            display: inline-block;
            padding: 8px 16px;
            border-radius: 20px;
            font-weight: bold;
            font-size: 16px;
            margin: 10px 0;
        }
        
        .status-mrsa {
            background: #dc2626;
            color: white;
        }
        
        .status-mssa {
            background: #16a34a;
            color: white;
        }
        
        .confidence-high {
            color: #16a34a;
            font-weight: bold;
        }
        
        .confidence-medium {
            color: #f59e0b;
            font-weight: bold;
        }
        
        .confidence-low {
            color: #dc2626;
            font-weight: bold;
        }
        
        .gene-grid {
            display: grid;
            grid-template-columns: repeat(auto-fill, minmax(200px, 1fr));
            gap: 15px;
            margin-top: 15px;
        }
        
        .gene-card {
            background: linear-gradient(135deg, #3b82f6 0%, #1e40af 100%);
            color: white;
            padding: 15px;
            border-radius: 8px;
            box-shadow: 0 2px 8px rgba(0, 0, 0, 0.1);
            text-align: center;
            font-weight: bold;
        }
        
        .detailed-table {
            width: 100%;
            border-collapse: collapse;
            margin-top: 15px;
            font-size: 14px;
        }
        
        .detailed-table th {
            background: #1e40af;
            color: white;
            padding: 12px;
            text-align: left;
            font-weight: bold;
            position: sticky;
            top: 0;
            border: 1px solid #e5e7eb;
        }
        
        .detailed-table td {
            padding: 10px;
            border-bottom: 1px solid #e5e7eb;
            border-right: 1px solid #e5e7eb;
        }
        
        .detailed-table tr:hover {
            background: #f3f4f6;
        }
        
        .highlight-red {
            background-color: #fee2e2 !important;
            color: #dc2626 !important;
            font-weight: bold;
        }
        
        .highlight-orange {
            background-color: #ffedd5 !important;
            color: #ea580c !important;
        }
        
        .highlight-yellow {
            background-color: #fef3c7 !important;
            color: #d97706 !important;
        }
        
        .kmer-table {
            width: 100%;
            border-collapse: collapse;
            margin-top: 15px;
            font-size: 14px;
        }
        
        .kmer-table th {
            background: #1e40af;
            color: white;
            padding: 12px;
            text-align: left;
            font-weight: bold;
        }
        
        .kmer-table td {
            padding: 10px;
            border-bottom: 1px solid #e5e7eb;
        }
        
        .kmer-table tr:hover {
            background: #f3f4f6;
        }
        
        .metrics-grid {
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(250px, 1fr));
            gap: 20px;
            margin-top: 15px;
        }
        
        .metric-card {
            background: linear-gradient(135deg, #8b5cf6 0%, #6d28d9 100%);
            color: white;
            padding: 20px;
            border-radius: 8px;
            box-shadow: 0 4px 12px rgba(0, 0, 0, 0.15);
        }
        
        .metric-label {
            font-size: 14px;
            opacity: 0.9;
            margin-bottom: 5px;
        }
        
        .metric-value {
            font-size: 24px;
            font-weight: bold;
        }
        
        .text-content {
            background: #f8fafc;
            border: 1px solid #e2e8f0;
            border-radius: 8px;
            padding: 20px;
            margin-top: 15px;
            font-family: 'Courier New', monospace;
            font-size: 14px;
            line-height: 1.5;
            white-space: pre-wrap;
            max-height: 500px;
            overflow-y: auto;
        }
        
        .file-source {
            font-size: 14px;
            color: #666;
            font-style: italic;
            margin-bottom: 10px;
            background: #f3f4f6;
            padding: 8px 12px;
            border-radius: 4px;
            display: inline-block;
        }
        
        .legend {
            display: flex;
            flex-wrap: wrap;
            gap: 15px;
            margin-top: 10px;
            margin-bottom: 15px;
            padding: 10px;
            background: #f8fafc;
            border-radius: 6px;
            border: 1px solid #e2e8f0;
        }
        
        .legend-item {
            display: flex;
            align-items: center;
            gap: 5px;
            font-size: 12px;
        }
        
        .legend-color {
            width: 15px;
            height: 15px;
            border-radius: 3px;
        }
        
        .footer {
            text-align: center;
            margin-top: 30px;
            padding: 20px;
            background: rgba(0, 0, 0, 0.3);
            border-radius: 10px;
            font-size: 14px;
        }
        
        .timestamp {
            color: #fbbf24;
            font-weight: bold;
        }
        
        .authorship {
            margin-top: 15px;
            padding: 15px;
            background: rgba(255, 255, 255, 0.1);
            border-radius: 8px;
            font-size: 12px;
        }
        
        .kmer-prediction-box {
            background: #f0f9ff;
            border-left: 4px solid #3b82f6;
            padding: 15px;
            margin: 15px 0;
            border-radius: 6px;
        }
        
        .kmer-warning-box {
            background: #fef3c7;
            border-left: 4px solid #f59e0b;
            padding: 15px;
            margin: 15px 0;
            border-radius: 6px;
        }
        
        @media (max-width: 768px) {
            .ascii-art {
                font-size: 6px;
            }
            .gene-grid {
                grid-template-columns: 1fr;
            }
            .detailed-table {
                font-size: 12px;
            }
            .detailed-table th, .detailed-table td {
                padding: 8px;
            }
        }
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
                <div class="quote-text" id="quoteText"></div>
                <div class="quote-author" id="quoteAuthor"></div>
            </div>
        </div>
'''
    
    # Add sample information
    html_content += f'''
        <div class="report-section">
            <h2>📊 Sample Information</h2>
            <div class="metrics-grid">
                <div class="metric-card">
                    <div class="metric-label">Sample Name</div>
                    <div class="metric-value">{sample_name}</div>
                </div>
                <div class="metric-card">
                    <div class="metric-label">Analysis Date</div>
                    <div class="metric-value">{time.strftime("%Y-%m-%d %H:%M:%S")}</div>
                </div>
                <div class="metric-card">
                    <div class="metric-label">Total Genes Detected</div>
                    <div class="metric-value">{len(ccrAB_genes) + len(ccrC_genes) + len(mec_genes)}</div>
                </div>
                <div class="metric-card">
                    <div class="metric-label">Database Used</div>
                    <div class="metric-value">{args.database_choice.title()}</div>
                </div>
            </div>
        </div>
'''
    
    # Add MRSA/MSSA status
    mrsa_status = enhanced_report['mrsa_status']
    status_class = "status-mrsa" if mrsa_status['classification'] == "CONFIRMED_MRSA" else "status-mssa"
    confidence_class = f"confidence-{mrsa_status['confidence'].lower()}"
    
    html_content += f'''
        <div class="report-section">
            <h2>🦠 MRSA/MSSA Classification</h2>
            <div class="status-badge {status_class}">{mrsa_status['classification']}</div>
            <p style="margin-top: 15px;"><strong>Confidence:</strong> <span class="{confidence_class}">{mrsa_status['confidence']}</span></p>
            <h3>Evidence:</h3>
            <ul>
'''
    for evidence in mrsa_status['evidence']:
        html_content += f'                <li>{evidence.replace("_", " ").title()}</li>\n'
    
    html_content += '''            </ul>
        </div>
'''
    
    # Add SCCmec typing results
    sccmec = enhanced_report['sccmec_typing']
    confidence_class = f"confidence-{sccmec['confidence'].lower()}"
    
    html_content += f'''
        <div class="report-section">
            <h2>🧬 SCCmec Typing Results</h2>
            <p><strong>Primary Type:</strong> <span style="font-size: 20px; color: #1e40af; font-weight: bold;">{sccmec['primary_type']}</span></p>
            <p><strong>Confidence:</strong> <span class="{confidence_class}">{sccmec['confidence']}</span></p>
            <p><strong>Typing Method:</strong> {sccmec['typing_method'].replace("_", " ").title()}</p>
'''
    
    # Add K-mer based prediction if gene-based failed
    if not sccmec_types and mrsa_gene and kmer_prediction.get('should_use_kmer', False):
        html_content += f'''
            <div class="kmer-prediction-box">
                <h3>🔍 K-mer Based Prediction (Gene-based typing inconclusive)</h3>
                <p><strong>{kmer_prediction.get('message', '')}</strong></p>
'''
        
        if kmer_prediction.get('top_hits'):
            html_content += '''
                <table class="detailed-table">
                    <thead>
                        <tr>
                            <th>Rank</th>
                            <th>Template</th>
                            <th>Template Coverage [%]</th>
                            <th>Query Coverage [%]</th>
                            <th>Score</th>
                        </tr>
                    </thead>
                    <tbody>
'''
            for i, hit in enumerate(kmer_prediction['top_hits'], 1):
                row_class = "highlight-red"  # ALWAYS HIGHLIGHT TOP 2 IN RED
                html_content += f'''                    <tr class="{row_class}">
                        <td>{i}</td>
                        <td><strong>{hit['template'].split('|')[0]}</strong></td>
                        <td>{hit['template_coverage']}</td>
                        <td>{hit['query_coverage']}</td>
                        <td>{hit['score']}</td>
                    </tr>
'''
            html_content += '''                    </tbody>
                </table>
'''
        
        if kmer_prediction.get('contradiction'):
            html_content += f'''
                <div class="kmer-warning-box">
                    <strong>⚠️ Interpretation Note:</strong> Top 2 k-mer hits have similar coverage. 
                    Consider both possibilities or perform additional testing.
                </div>
'''
        
        html_content += '''            </div>
'''
    
    if len(sccmec['all_predicted_types']) > 1:
        html_content += '''
            <h3>All Predicted Types:</h3>
            <ul>
'''
        for stype in sccmec['all_predicted_types']:
            html_content += f'                <li>{stype}</li>\n'
        html_content += '''            </ul>
'''
    
    html_content += '''        </div>
'''
    
    # Add Detailed Text Results (from sccmec_detailed_results.txt)
    if detailed_text_content:
        html_content += f'''
        <div class="report-section">
            <h2>📄 Detailed Text Results</h2>
            <p class="file-source">From: {args.output_file}</p>
            <div class="text-content">
{detailed_text_content}
            </div>
        </div>
'''
    
    # Add detailed gene detection results
    if detailed_gene_results:
        html_content += '''
        <div class="report-section">
            <h2>🧪 Detailed Gene Detection Results</h2>
            <p class="file-source">From: results_tab_MyDbFinder.txt</p>
            <div style="max-height: 600px; overflow-y: auto;">
                <table class="detailed-table">
                    <thead>
                        <tr>
                            <th>Gene</th>
                            <th>% Identity</th>
                            <th>Query/HSP Length</th>
                            <th>Contig</th>
                            <th>Position</th>
                        </tr>
                    </thead>
                    <tbody>
'''
        for result in detailed_gene_results:  # Show ALL results
            html_content += f'''
                        <tr>
                            <td><strong>{result['gene_full']}</strong></td>
                            <td>{result['identity']}</td>
                            <td>{result['query_length']}/{result['hsp_length']}</td>
                            <td>{result['contig']}</td>
                            <td>{result['position']}</td>
                        </tr>
'''
        html_content += f'''
                    </tbody>
                </table>
            </div>
            <p style="margin-top: 10px; font-size: 14px; color: #666;">
                Total genes detected: {len(detailed_gene_results)}
            </p>
        </div>
'''
    
    # Add COMPLETE results_MyKmerFinder.txt content in BEAUTIFUL COLORED TABLE
    if kmer_headers and kmer_rows:
        # Sort rows by template coverage (8th column, index 7)
        try:
            sorted_rows = sorted(kmer_rows, key=lambda x: float(x[7] if len(x) > 7 else 0), reverse=True)
        except (ValueError, IndexError):
            sorted_rows = kmer_rows
        
        html_content += '''
        <div class="report-section">
            <h2>📋 Complete K-mer Results (Formatted Table)</h2>
            <p class="file-source">From: results_MyKmerFinder.txt (All results in colored table)</p>
            
            <div class="legend">
                <div class="legend-item">
                    <div class="legend-color" style="background-color: #fee2e2;"></div>
                    <span>Top 2 hits with highest template coverage</span>
                </div>
                <div class="legend-item">
                    <div class="legend-color" style="background-color: #ffedd5;"></div>
                    <span>High coverage (≥ 70%)</span>
                </div>
                <div class="legend-item">
                    <div class="legend-color" style="background-color: #fef3c7;"></div>
                    <span>Medium coverage (≥ 50%)</span>
                </div>
            </div>
            
            <div style="max-height: 600px; overflow-y: auto;">
                <table class="detailed-table">
                    <thead>
                        <tr>
'''
        # Add headers
        for i, header in enumerate(kmer_headers):
            html_content += f'                            <th>{header}</th>\n'
        
        html_content += '''                        </tr>
                    </thead>
                    <tbody>
'''
        # Add rows with coloring
        for i, row in enumerate(sorted_rows):
            # Determine row class based on template coverage
            row_class = ""
            try:
                if len(row) > 7:
                    template_cov = float(row[7])
                    if i < 2:  # Top 2 rows
                        row_class = "highlight-red"
                    elif template_cov >= 70.0:
                        row_class = "highlight-orange"
                    elif template_cov >= 50.0:
                        row_class = "highlight-yellow"
            except (ValueError, IndexError):
                pass
            
            html_content += f'                        <tr class="{row_class}">\n'
            
            for j, cell in enumerate(row):
                if j == 0:  # Template column - format for better display
                    if '|' in cell:
                        parts = cell.split('|')
                        if len(parts) > 1 and not parts[1].startswith('gb'):
                            display_cell = f"{parts[0]} | {parts[1]}"
                        else:
                            display_cell = parts[0]
                        html_content += f'                            <td><strong>{display_cell}</strong></td>\n'
                    else:
                        html_content += f'                            <td><strong>{cell}</strong></td>\n'
                else:
                    html_content += f'                            <td>{cell}</td>\n'
            
            html_content += '                        </tr>\n'
        
        html_content += f'''                    </tbody>
                </table>
            </div>
            <p style="margin-top: 10px; font-size: 14px; color: #666;">
                Total templates analyzed: {len(sorted_rows)}
                <br><span style="color: #dc2626; font-weight: bold;">● Top 2 hits highlighted in RED</span>
                <br><span style="color: #ea580c;">● High coverage (≥70%) in ORANGE</span>
                <br><span style="color: #d97706;">● Medium coverage (≥50%) in YELLOW</span>
            </p>
        </div>
'''
    
    # Add parsed k-mer results table - HIGHLIGHTING TOP 2
    if detailed_kmer_results:
        # Sort by template coverage for highlighting
        sorted_kmer = sorted(detailed_kmer_results, 
                           key=lambda x: float(x.get('template_coverage', 0)), 
                           reverse=True)
        
        html_content += '''
        <div class="report-section">
            <h2>🔍 Detailed K-mer Homology Analysis (Parsed View)</h2>
            <p class="file-source">Parsed from: results_MyKmerFinder.txt</p>
            <div style="max-height: 600px; overflow-y: auto;">
                <table class="detailed-table">
                    <thead>
                        <tr>
                            <th>Rank</th>
                            <th>Template</th>
                            <th>Score</th>
                            <th>Expected</th>
                            <th>z</th>
                            <th>p-value</th>
                            <th>query coverage [%]</th>
                            <th>template coverage [%]</th>
                            <th>depth</th>
                            <th>Kmers in Template</th>
                        </tr>
                    </thead>
                    <tbody>
'''
        for i, result in enumerate(sorted_kmer):  # Show ALL results
            # HIGHLIGHT TOP 2 WITH HIGHEST COVERAGE IN RED
            template_cov = float(result.get('template_coverage', 0))
            row_class = "highlight-red" if i < 2 else ""  # TOP 2 HIGHLIGHTED
            
            # Format template name for cleaner display
            template = result['template']
            if '|' in template:
                parts = template.split('|')
                if len(parts) > 1 and not parts[1].startswith('gb'):
                    display_template = f"{parts[0]} | {parts[1]}"
                else:
                    display_template = parts[0]
            else:
                display_template = template
            
            html_content += f'''
                        <tr class="{row_class}">
                            <td>{i+1}</td>
                            <td><strong>{display_template}</strong></td>
                            <td>{result['score']}</td>
                            <td>{result['expected']}</td>
                            <td>{result['z']}</td>
                            <td>{result['p_value']}</td>
                            <td>{result['query_coverage']}</td>
                            <td>{result['template_coverage']}</td>
                            <td>{result['depth']}</td>
                            <td>{result['kmers_in_template']}</td>
                        </tr>
'''
        html_content += f'''
                    </tbody>
                </table>
            </div>
            <p style="margin-top: 10px; font-size: 14px; color: #666;">
                Total templates analyzed: {len(sorted_kmer)}
            </p>
        </div>
'''
    
    # Add gene detection summary
    html_content += '''
        <div class="report-section">
            <h2>🧬 Gene Detection Summary</h2>
'''
    
    # MEC genes
    if enhanced_report['gene_detection']['mec_genes']:
        html_content += '''
            <h3>MEC Genes</h3>
            <div class="gene-grid">
'''
        for gene_data in enhanced_report['gene_detection']['mec_genes']:
            html_content += f'''                <div class="gene-card">{gene_data['gene']}</div>
'''
        html_content += '''            </div>
'''
    
    # CCR genes
    if enhanced_report['gene_detection']['ccr_genes']:
        html_content += '''
            <h3>CCR Genes</h3>
            <div class="gene-grid">
'''
        for gene_data in enhanced_report['gene_detection']['ccr_genes']:
            html_content += f'''                <div class="gene-card">{gene_data['gene']}</div>
'''
        html_content += '''            </div>
'''
    
    # Regulator genes
    if enhanced_report['gene_detection']['regulator_genes']:
        html_content += '''
            <h3>Regulator Genes</h3>
            <div class="gene-grid">
'''
        for gene_data in enhanced_report['gene_detection']['regulator_genes']:
            html_content += f'''                <div class="gene-card">{gene_data['gene']}</div>
'''
        html_content += '''            </div>
'''
    
    html_content += '''        </div>
'''
    
    # Add HTML K-mer Report
    if html_kmer_content:
        html_content += f'''
        <div class="report-section">
            <h2>📋 HTML K-mer Report</h2>
            <p class="file-source">From: HTML_MyKmerFinder.txt</p>
            <div class="text-content">
{html_kmer_content}
            </div>
        </div> 
'''
    
    # Add k-mer summary results if available
    if kmer_hits:
        html_content += '''
        <div class="report-section">
            <h2>🎯 K-mer Homology Summary (Top 10 Hits)</h2>
            <table class="kmer-table">
                <thead>
                    <tr>
                        <th>Rank</th>
                        <th>Template</th>
                        <th>Score</th>
                        <th>Template Coverage (%)</th>
                        <th>Query Coverage (%)</th>
                        <th>Depth</th>
                    </tr>
                </thead>
                <tbody>
'''
        for i, (hit, data) in enumerate(kmer_hits[:10], 1):
            template_cov = float(data['template_coverage'])
            row_class = "highlight-red" if i <= 2 and template_cov >= 50.0 else ""
            
            html_content += f'''                    <tr class="{row_class}">
                        <td>{i}</td>
                        <td><strong>{hit}</strong></td>
                        <td>{data['Score']}</td>
                        <td>{data['template_coverage']}</td>
                        <td>{data['query_coverage']}</td>
                        <td>{data['depth']}</td>
                    </tr>
'''
        html_content += '''                </tbody>
            </table>
        </div>
'''
    
    # Add quality metrics
    quality = enhanced_report['quality_metrics']
    html_content += f'''
        <div class="report-section">
            <h2>📈 Quality Metrics</h2>
            <div class="metrics-grid">
                <div class="metric-card">
                    <div class="metric-label">Gene Completeness</div>
                    <div class="metric-value">{quality['gene_completeness'].upper()}</div>
                </div>
                <div class="metric-card">
                    <div class="metric-label">Type Contradictions</div>
                    <div class="metric-value">{"YES" if quality['contradictions'] else "NO"}</div>
                </div>
                <div class="metric-card">
                    <div class="metric-label">Analysis Parameters</div>
                    <div class="metric-value">ID: {args.id_threshold}% Len: {args.len_threshold}%</div>
                </div>
                <div class="metric-card">
                    <div class="metric-label">K-mer Prediction</div>
                    <div class="metric-value">{"AVAILABLE" if kmer_prediction.get('should_use_kmer') else "NOT USED"}</div>
                </div>
            </div>
        </div>
'''
    
    # Add footer with JavaScript for rotating quotes and authorship
    html_content += f'''
        <div class="footer">
            <p><strong>STAPHSCOPE</strong> - SCCmec Cassette Finder</p>
            <p class="timestamp">Generated: {time.strftime("%Y-%m-%d %H:%M:%S")}</p>
            <p style="margin-top: 10px; font-size: 12px;">Powered by MyDbFinder & MyKmerFinder</p>
            
            <div class="authorship">
                <p><strong>Technical Support & Inquiries:</strong></p>
                <p>Author: Brown Beckley</p>
                <p>GitHub: <a href="https://github.com/bbeckley-hub" style="color: #fbbf24;">bbeckley-hub</a></p>
                <p>Email: <a href="mailto:brownbeckley94@gmail.com" style="color: #fbbf24;">brownbeckley94@gmail.com</a></p>
                <p>Affiliation: University of Ghana Medical School - Department of Medical Biochemistry</p>
            </div>
        </div>
    </div>

    <script>
        const quotes = [
            {{ text: "The important thing is not to stop questioning. Curiosity has its own reason for existing.", author: "Albert Einstein" }},
            {{ text: "Science is not only a disciple of reason but also one of romance and passion.", author: "Stephen Hawking" }},
            {{ text: "Somewhere, something incredible is waiting to be known.", author: "Carl Sagan" }},
            {{ text: "The good thing about science is that it's true whether or not you believe in it.", author: "Neil deGrasse Tyson" }},
            {{ text: "In science, there are no shortcuts to truth.", author: "Karl Popper" }},
            {{ text: "Science knows no country, because knowledge belongs to humanity.", author: "Louis Pasteur" }},
            {{ text: "The science of today is the technology of tomorrow.", author: "Edward Teller" }},
            {{ text: "Nothing in life is to be feared, it is only to be understood.", author: "Marie Curie" }},
            {{ text: "Research is what I'm doing when I don't know what I'm doing.", author: "Wernher von Braun" }},
            {{ text: "The universe is not required to be in perfect harmony with human ambition.", author: "Carl Sagan" }},
            {{"text": "Staphscope represents the convergence of genomic surveillance and clinical diagnostics, transforming raw sequences into actionable insights for infection control.", "author": "Brown Beckley"}},
            {{"text": "In the battle against antimicrobial resistance, tools like Staphscope are our eyes and ears, revealing the genetic blueprints of resistant pathogens.", "author": "Brown Beckley"}},
            {{"text": "Staphscope isn't just a tool; it's a comprehensive system that bridges the gap between sequencing data and public health action.", "author": "Brown Beckley"}},    
            {{"text": "Through Staphscope, we turn the complexity of bacterial genomes into clear, interpretable reports, empowering clinicians and researchers alike.", "author": "Brown Beckley"}},
            {{"text": "Staphscope is a testament to the power of bioinformatics in the modern era, making advanced pathogen typing accessible to all.", "author": "Brown Beckley"}}
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
                quoteText.textContent = `"${{quote.text}}"`;
                quoteAuthor.textContent = `— ${{quote.author}}`;
                quoteContainer.style.opacity = '1';
            }}, 500);
        }}

        displayQuote();
        setInterval(displayQuote, 10000);
    </script>
</body>
</html>
'''
    
    # Save HTML report
    html_file = output_path / "staphscope_comprehensive_report.html"
    with open(html_file, 'w', encoding='utf-8') as f:
        f.write(html_content)
    
    return html_file
    
def determine_primary_type_with_confidence(sccmec_types: List[str], kmer_results: List[Dict],
                                         mec_genes: Set[str], ccrAB_genes: Set[str], 
                                         ccrC_genes: Set[str]) -> Tuple[str, str, List[str]]:
    """Determine primary type with confidence scoring and evidence tracking"""
    evidence = []
    
    if not sccmec_types:
        # Check if we have k-mer evidence even without gene-based types
        if kmer_results:
            best_kmer = kmer_results[0]
            coverage = float(best_kmer.get("template_coverage", 0))
            primary_type = best_kmer["template"].split("|")[0]
            evidence.append(f"kmer_only_evidence_{coverage:.1f}%")
            
            if coverage >= 85.0:
                return primary_type, "HIGH", evidence
            elif coverage >= 70.0:
                return primary_type, "MEDIUM", evidence
            else:
                return primary_type, "LOW", evidence
        return "NA", "LOW", ["no_sccmec_types_detected"]
    
    if len(sccmec_types) == 1:
        # Single prediction - check k-mer support
        primary_type = sccmec_types[0]
        evidence.append("single_gene_based_prediction")
        
        if kmer_results:
            best_kmer = kmer_results[0]
            coverage = float(best_kmer.get("template_coverage", 0))
            kmer_type = best_kmer["template"].split("|")[0]
            
            if kmer_type == primary_type:
                evidence.append(f"kmer_support_{coverage:.1f}%_coverage")
                if coverage >= 85.0:
                    return primary_type, "VERY_HIGH", evidence
                elif coverage >= 70.0:
                    return primary_type, "HIGH", evidence
                else:
                    return primary_type, "MEDIUM", evidence
            else:
                evidence.append(f"kmer_conflict_{kmer_type}_{coverage:.1f}%")
                # Still use gene-based but with lower confidence
                return primary_type, "MEDIUM", evidence
        
        # Gene-based only
        if len(mec_genes) >= 2 and len(ccrAB_genes) >= 2:
            evidence.append("complete_gene_set")
            return primary_type, "HIGH", evidence
        else:
            evidence.append("partial_gene_set")
            return primary_type, "MEDIUM", evidence
    
    else:
        # Multiple predictions - use k-mer to resolve
        evidence.append("multiple_gene_based_predictions")
        
        if kmer_results:
            best_kmer = kmer_results[0]
            coverage = float(best_kmer.get("template_coverage", 0))
            primary_type = best_kmer["template"].split("|")[0]
            
            if coverage >= 70.0 and primary_type in sccmec_types:
                evidence.append(f"kmer_resolution_{coverage:.1f}%_coverage")
                confidence = "VERY_HIGH" if coverage >= 90.0 else "HIGH" if coverage >= 85.0 else "MEDIUM"
                return primary_type, confidence, evidence
            elif primary_type in sccmec_types:
                evidence.append(f"weak_kmer_support_{coverage:.1f}%")
                return primary_type, "LOW", evidence
        
        # Fallback to first gene-based prediction
        evidence.append("first_gene_prediction_fallback")
        return sccmec_types[0], "LOW", evidence

def assess_gene_completeness(ccrAB_genes: Set[str], ccrC_genes: Set[str], mec_genes: Set[str]) -> str:
    """Assess completeness of gene detection"""
    mec_count = len(mec_genes)
    ccr_count = len(ccrAB_genes) + len(ccrC_genes)
    
    if mec_count >= 2 and ccr_count >= 2:
        return "complete"
    elif mec_count >= 1 or ccr_count >= 1:
        return "partial"
    else:
        return "none"

def main():
    start_time = time.time()
    
    # Parse command-line arguments
    parser = argparse.ArgumentParser(
        description='Prediction of SCCmec cassette in S. aureus',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("-iDb", "--input-db", dest="fasta_file_db", required=True,
                        help="Fasta input file for MyDbFinder")
    parser.add_argument("-iKm", "--input-km", dest="fasta_file_km", required=True,
                        help="Fasta input file for MyKmerFinder")
    parser.add_argument("-k", "--id-threshold", dest="id_threshold", default="90",
                        help="Minimum identity threshold")
    parser.add_argument("-l", "--len-threshold", dest="len_threshold", default="60",
                        help="Minimum length threshold")
    parser.add_argument("-o", "--output", dest="output_file", required=True,
                        help="Output file name")
    parser.add_argument("-d", "--output-dir", dest="output_dir", required=True,
                        help="Output directory")
    parser.add_argument("-db_dir", "--database-dir", dest="database_dir", required=True,
                        help="Database directory")
    parser.add_argument("-sc_dir", "--script-dir", dest="script_dir", required=True,
                        help="Script directory")
    parser.add_argument("-db_choice", "--db-choice", dest="database_choice", required=True,
                        choices=['reference', 'extended'],
                        help="Database choice: reference or extended")
    args = parser.parse_args()
    
    print(f"Finding the SCCmec cassette of: {args.fasta_file_db}")
    
    # Create output directory if needed
    output_path = Path(args.output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    
    # Initialize variables
    ccrAB_genes = set()
    ccrC_genes = set()
    mec_genes = set()
    subtyping_genes = set()
    mrsa_gene = ""
    total_genes = set()
    
    # Run MyDbFinder
    print("Running the BLAST-based approach")
    db_finder_cmd = [
        "python3",
        str(Path(args.script_dir) / "CGE_MyDbFinder-1.1.py"),
        "-o", str(output_path),
        "-r", str(Path(args.database_dir) / "single_genes_database_20171117.fasta"),
        "-k", args.id_threshold,
        "-l", args.len_threshold,
        "-i", args.fasta_file_db
    ]
    run_command(db_finder_cmd)
    
    # Parse MyDbFinder results
    db_finder_file = output_path / "results_tab_MyDbFinder.txt"
    detailed_gene_results = []
    if not db_finder_file.exists():
        print("Error: MyDbFinder did not produce a result file")
        sys.exit(1)
    
    # Check if file is empty before processing
    if os.path.getsize(db_finder_file) == 0:
        print("Warning: MyDbFinder results file is empty")
    else:
        detailed_gene_results = parse_detailed_gene_results(db_finder_file)
        with open(db_finder_file, 'r') as f:
            # Safely skip header only if it exists
            header = next(f, None)
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) < 4:
                    continue
                
                gene = parts[0].split(':')[0]
                if gene.startswith('ccr'):
                    if gene.startswith(('ccrA', 'ccrB')):
                        ccrAB_genes.add(gene)
                    elif gene.startswith('ccrC'):
                        ccrC_genes.add(gene)
                elif gene.startswith(('mec', 'IS', 'dme')):
                    mec_genes.add(gene)
                    if gene.startswith('mecALGA251'):
                        mrsa_gene = "mecALGA251"
                    elif gene.startswith('mecA'):
                        mrsa_gene = "mecA"
                    elif gene.startswith('mec-class'):
                        mrsa_gene = "mecA"
                elif gene.startswith('sub'):
                    subtyping_genes.add(gene.split("|")[0])

    # Run MyKmerFinder
    print("Running the kmer-based approach")
    template_db = "MyKmerFinder_reference_template" if args.database_choice == 'reference' else "MyKmerFinder_extended_template"
    kmer_finder_cmd = [
        "python3",
        str(Path(args.script_dir) / "findtemplate.py"),
        "-i", args.fasta_file_km,
        "-t", str(Path(args.database_dir) / "template_db" / template_db),
        "-o", str(output_path / "results_MyKmerFinder.txt")
    ]
    run_command(kmer_finder_cmd)
    
    # Parse MyKmerFinder results
    kmer_results = {}
    kmer_file = output_path / "results_MyKmerFinder.txt"
    detailed_kmer_results = []
    if kmer_file.exists() and os.path.getsize(kmer_file) > 0:
        detailed_kmer_results = parse_kmer_detailed_results(kmer_file)
        with open(kmer_file, 'r') as f:
            next(f)  # Skip header
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) < 10:
                    continue
                
                template_coverage = float(parts[6])
                if template_coverage >= 50.0:  # Cutoff value
                    kmer_results[parts[0]] = {
                        'Score': parts[1],
                        'Expected': parts[2],
                        'z': parts[3],
                        'p_value': parts[4],
                        'query_coverage': parts[5],
                        'template_coverage': parts[6],
                        'depth': parts[7],
                        'Kmers_in_Template': parts[8],
                        'Description': parts[9]
                    }
    
    # Perform gene typing
    total_genes = ccrAB_genes | ccrC_genes | mec_genes
    
    ccr_classes = perform_ccr_gene_complex_typing(ccrAB_genes, ccrC_genes)
    mec_classes = perform_mec_gene_complex_typing(mec_genes)
    all_classes = set(ccr_classes + mec_classes)
    
    sccmec_types = perform_sccmec_typing(all_classes)
    
    # Determine if subtyping should be performed
    perform_subtyping = 'No'
    if len(sccmec_types) == 1 and sccmec_types[0] in {
        'SCCmec_type_IV(2B)', 
        'SCCmec_type_V(5C2&5)', 
        'SCCmec_type_V(5C2)'
    }:
        perform_subtyping = 'Yes'
    
    # Process MyKmerFinder results
    kmer_hits = sorted(
        kmer_results.items(),
        key=lambda x: float(x[1]['Score']),
        reverse=True
    ) if kmer_results else []
    best_kmer_hit = kmer_hits[0][0].split("|")[0] if kmer_hits else ""
    kmer_subtype = kmer_hits[0][0].split("|")[1] if kmer_hits and "|" in kmer_hits[0][0] else ""
    
    # Get k-mer based prediction for isolates with mecA but no gene-based type
    has_mecA = "mecA" in mec_genes or "mecA" in mrsa_gene
    kmer_prediction = get_top_kmer_predictions(detailed_kmer_results, has_mecA)
    
    # Generate final text results
    result_file = output_path / args.output_file
    with open(result_file, 'w') as out:
        # MRSA/MSSA determination
        if mrsa_gene:
            out.write(f"The input organism was predicted as a MRSA isolate\nThe {mrsa_gene} gene was detected\n\n")
        else:
            out.write("The input organism was predicted as a MSSA isolate\nThe mecA/mecC gene was not detected\n\n")
        
        # SCCmec typing results
        if sccmec_types:
            if len(sccmec_types) > 1:
                out.write(f"Alert! Possible {len(sccmec_types)} SCCmec cassettes were predicted:\n")
                out.write("\n".join(sccmec_types) + "\n\n")
            else:
                db_hit = sccmec_types[0]
                
                # Handle subtyping if needed
                if perform_subtyping == 'Yes' and kmer_subtype:
                    if len(subtyping_genes) > 1:
                        out.write("Alert! Multiple subtype target genes detected.\n")
                    elif len(subtyping_genes) == 1:
                        subtype_gene = list(subtyping_genes)[0]
                        db_subtype = f"{db_hit.split('(')[0]}_{subtype_gene.split('-')[1]}"
                        
                        if db_subtype != kmer_subtype:
                            out.write("Alert! Contradicting subtype predictions.\n")
                        out.write(f"Prediction based on genes: {db_subtype}\n")
                
                out.write(f"Prediction based on genes: {db_hit}\n")
                
                if best_kmer_hit:
                    coverage = kmer_hits[0][1]['template_coverage']
                    out.write(f"Prediction based on homology: {best_kmer_hit} ({coverage}% coverage)\n")
                
                # Find additional complexes/genes
                additional_complexes = all_classes - SCCMEC_CLASSES.get(db_hit, set())
                additional_genes = total_genes - SCCMEC_DEFINITIONS.get(db_hit, set())
                
                if additional_complexes:
                    out.write("\nAdditional complexes found:\n")
                    out.write("\n".join(additional_complexes) + "\n")
                
                if additional_genes:
                    out.write("\nAdditional genes found:\n")
                    out.write("\n".join(additional_genes) + "\n")
        else:
            out.write("No SCCmec element was detected by gene-based prediction\n")
            
            # Add k-mer based prediction if available
            if kmer_prediction['should_use_kmer']:
                out.write("\n=== K-MER BASED PREDICTION ===\n")
                out.write(f"{kmer_prediction['message']}\n\n")
                out.write("Top k-mer hits:\n")
                for i, hit in enumerate(kmer_prediction['top_hits'], 1):
                    out.write(f"{i}. {hit['template']} - Template Coverage: {hit['template_coverage']}%, Query Coverage: {hit['query_coverage']}%\n")
                
                if kmer_prediction['contradiction']:
                    out.write("\n⚠️ NOTE: Top 2 hits have similar coverage. Consider both possibilities.\n")
            
            if all_classes:
                out.write("\nDetected gene complexes:\n")
                out.write("\n".join(all_classes) + "\n")
            
            if kmer_hits:
                best_hit = kmer_hits[0][0]
                coverage = kmer_hits[0][1]['template_coverage']
                out.write(f"\nBest homology match: {best_hit} ({coverage}% coverage)\n")
                
                # Find missing genes for best hit
                missing_genes = SCCMEC_DEFINITIONS.get(best_hit.split("|")[0], set()) - total_genes
                if missing_genes:
                    out.write("\nMissing genes for this cassette:\n")
                    out.write("\n".join(missing_genes) + "\n")
    
    # Generate enhanced JSON report
    try:
        enhanced_report = generate_enhanced_json_report(
            sample_name=Path(args.fasta_file_db).name,
            mrsa_gene=mrsa_gene,
            sccmec_types=sccmec_types,
            kmer_hits=kmer_hits,
            kmer_prediction=kmer_prediction,
            ccrAB_genes=ccrAB_genes,
            ccrC_genes=ccrC_genes,
            mec_genes=mec_genes,
            subtyping_genes=subtyping_genes,
            total_genes=total_genes,
            output_path=output_path
        )
        print(f"Enhanced JSON report saved to: {output_path / 'sccmec_enhanced_report.json'}")
    except Exception as e:
        print(f"Warning: Could not generate enhanced JSON report: {e}")
    
    # Generate beautiful StaphScope HTML report with detailed results
    try:
        html_file = generate_staphscope_html_report(
            sample_name=Path(args.fasta_file_db).name,
            mrsa_gene=mrsa_gene,
            sccmec_types=sccmec_types,
            enhanced_report=enhanced_report,
            ccrAB_genes=ccrAB_genes,
            ccrC_genes=ccrC_genes,
            mec_genes=mec_genes,
            kmer_hits=kmer_hits,
            detailed_gene_results=detailed_gene_results,
            detailed_kmer_results=detailed_kmer_results,
            output_path=output_path,
            args=args,
            kmer_prediction=kmer_prediction
        )
        print(f"Beautiful StaphScope HTML report saved to: {html_file}")
    except Exception as e:
        print(f"Warning: Could not generate StaphScope HTML report: {e}")
    
    # Generate original HTML report for MyKmerFinder
    html_file = output_path / "HTML_MyKmerFinder.txt"
    with open(html_file, 'w') as html:
        html.write("SCCmec elements\n")
        if kmer_hits:
            html.write("Template\tScore\tExpected\tz\tp_value\tQuery Coverage\tTemplate Coverage\tDepth\tKmers\tDescription\n")
            for hit, data in kmer_hits[:5]:  # Top 5 hits
                html.write(f"{hit}\t{data['Score']}\t{data['Expected']}\t{data['z']}\t{data['p_value']}\t")
                html.write(f"{data['query_coverage']}\t{data['template_coverage']}\t{data['depth']}\t")
                html.write(f"{data['Kmers_in_Template']}\t{data['Description']}\n")
        else:
            html.write("no whole SCCmec cassette was found\n")
    
    # Cleanup temporary files
    for pattern in ["MyKmerFinder_template.*", "Hit_in_genome_seq.fsa", 
                    "results.txt", "Database_gene_seq.fsa"]:
        for file in output_path.glob(pattern):
            try:
                file.unlink()
            except OSError:
                pass
    
    runtime = time.time() - start_time
    print(f"Prediction completed in {runtime:.2f} seconds")
    print(f"Results saved to: {result_file}")

if __name__ == "__main__":
    main()

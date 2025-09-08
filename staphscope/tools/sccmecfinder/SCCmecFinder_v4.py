#!/usr/bin/env python3
"""
SCCmecFinder_v4 - Python 3 Implementation
Original Author: Hulya Kaya
Rewritten by: Beckley Brown <brownbeckley94@gmail.com>
Date: 2025-08-18
Updated: 2025-08-19 (Fixed StopIteration error and improved error handling)
"""

import os
import sys
import time
import glob
import argparse
import subprocess
from pathlib import Path
from typing import Dict, List, Set, Tuple, Optional

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
    'SCCmec_type_XII(9C2)': {'ccrC21', 'mec-class-C2', 'mecA'},
    'SCCmec_type_XIII(9A)': {'ccrC21', 'mecA', 'dmecR1', 'IS1272'}
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
    if "ccrC21" in ccrC_genes:
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
        "perl",
        str(Path(args.script_dir) / "CGE_MyDbFinder-1.1.pl"),
        "-o", str(output_path),
        "-r", str(Path(args.database_dir) / "single_genes_database_20171117.fasta"),
        "-k", args.id_threshold,
        "-l", args.len_threshold,
        "-i", args.fasta_file_db
    ]
    run_command(db_finder_cmd)
    
    # Parse MyDbFinder results - FIXED SECTION
    db_finder_file = output_path / "results_tab_MyDbFinder.txt"
    if not db_finder_file.exists():
        print("Error: MyDbFinder did not produce a result file")
        sys.exit(1)
    
    # Check if file is empty before processing
    if os.path.getsize(db_finder_file) == 0:
        print("Warning: MyDbFinder results file is empty")
    else:
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
    if kmer_file.exists() and os.path.getsize(kmer_file) > 0:
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
    
    # Generate final results
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
            out.write("No SCCmec element was detected\n")
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
    
    # Generate HTML report
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

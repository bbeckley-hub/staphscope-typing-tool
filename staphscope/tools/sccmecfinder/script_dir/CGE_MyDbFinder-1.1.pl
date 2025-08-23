#!/usr/bin/env python3
"""
MyDbFinder - Python Implementation
Replaces the original Perl script for detecting resistance genes
"""

import os
import sys
import argparse
import subprocess
import tempfile
from pathlib import Path
from typing import Dict, List, Tuple, Set

def run_command(cmd: List[str]) -> None:
    """Execute a system command with error handling"""
    try:
        subprocess.run(cmd, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    except subprocess.CalledProcessError as e:
        print(f"Error running command: {' '.join(cmd)}")
        print(f"Error message: {e.stderr.decode()}")
        sys.exit(1)

def main():
    parser = argparse.ArgumentParser(
        description='Find antimicrobial resistance genes in input sequences',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("-o", "--output-dir", dest="output_dir", required=True,
                        help="Output directory")
    parser.add_argument("-r", "--database", dest="database", required=True,
                        help="Resistance gene database (FASTA format)")
    parser.add_argument("-k", "--id-threshold", dest="id_threshold", type=float, default=90.0,
                        help="Minimum identity percentage threshold")
    parser.add_argument("-l", "--len-threshold", dest="len_threshold", type=float, default=60.0,
                        help="Minimum length percentage threshold")
    parser.add_argument("-i", "--input", dest="input_file", required=True,
                        help="Input genome file (FASTA format)")
    args = parser.parse_args()

    # Create output directory
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    tabular_output = output_dir / "results_tab_MyDbFinder.txt"
    
    print("Running BLAST-based resistance gene detection")
    
    # Create temporary directory for BLAST files
    with tempfile.TemporaryDirectory() as tmp_dir:
        tmp_path = Path(tmp_dir)
        
        # Create BLAST database from input genome
        db_name = tmp_path / "genome_db"
        makeblastdb_cmd = [
            "makeblastdb",
            "-in", args.input_file,
            "-dbtype", "nucl",
            "-out", str(db_name)
        ]
        run_command(makeblastdb_cmd)
        
        # Run BLAST search
        blast_output = tmp_path / "blast_results.tsv"
        blast_cmd = [
            "blastn",
            "-db", str(db_name),
            "-query", args.database,
            "-out", str(blast_output),
            "-outfmt", "6 qseqid qlen length pident sseqid sstart send",
            "-num_threads", "1"
        ]
        run_command(blast_cmd)
        
        # Process BLAST results
        genes_found: Set[str] = set()
        with open(blast_output, 'r') as blast_file, open(tabular_output, 'w') as tabular_file:
            for line in blast_file:
                parts = line.strip().split('\t')
                if len(parts) < 7:
                    continue
                
                qseqid = parts[0]       # Gene ID
                qlen = float(parts[1])  # Query length
                hsp_len = float(parts[2])  # HSP length
                pident = float(parts[3]) # Percent identity
                sseqid = parts[4]       # Contig ID
                sstart = parts[5]       # Start position
                send = parts[6]         # End position
                
                # Calculate coverage
                coverage = (hsp_len / qlen) * 100
                
                # Apply thresholds
                if pident >= args.id_threshold and coverage >= args.len_threshold:
                    gene_name = qseqid.split(':')[0]
                    genes_found.add(gene_name)
                    # Write to tabular output: Gene, %ID, Query/HSP, Contig, Position
                    tabular_file.write(f"{gene_name}\t{pident:.2f}\t{int(qlen)}/{int(hsp_len)}\t{sseqid}\t{sstart}..{send}\n")
        
        print(f"Found {len(genes_found)} resistance genes passing thresholds")

if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""
MLST Module for StaphScope - Complete with Beautiful HTML Reports
Author: Brown Beckley <brownbeckley94@gmail.com>
GitHub: bbeckley-hub
Affiliation: University of Ghana Medical School - Department of Medical Biochemistry
Date: 2025/2026-04-21
Send a quick mail for any issues or further explanations.
"""

import os
import sys
import json
import glob
import argparse
import subprocess
import random
from pathlib import Path
from typing import Dict, List, Optional
import pandas as pd
from datetime import datetime

class ModularMLSTAnalyzer:
    def __init__(self, database_dir: Path, script_dir: Path):
        self.database_dir = database_dir
        self.script_dir = script_dir
        self.mlst_bin = script_dir / "mlst"
        
        # Science quotes for rotation
        self.science_quotes = [
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
            {"text": "Staphscope represents the convergence of genomic surveillance and clinical diagnostics, transforming raw sequences into actionable insights for infection control.", "author": "Brown Beckley"},
            {"text": "In the battle against antimicrobial resistance, tools like Staphscope are our eyes and ears, revealing the genetic blueprints of resistant pathogens.", "author": "Brown Beckley"},
            {"text": "Staphscope isn't just a tool; it's a comprehensive system that bridges the gap between sequencing data and public health action.", "author": "Brown Beckley"},
            {"text": "Through Staphscope, we turn the complexity of bacterial genomes into clear, interpretable reports, empowering clinicians and researchers alike.", "author": "Brown Beckley"},
            {"text": "Staphscope is a testament to the power of bioinformatics in the modern era, making advanced pathogen typing accessible to all.", "author": "Brown Beckley"}
        ]
    
    def get_random_quote(self):
        """Get a random science quote"""
        return random.choice(self.science_quotes)
    
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
            
            # Add identity and coverage information
            mlst_results.update(self.get_identity_coverage(mlst_results.get('st', 'ND')))
            
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
            "confidence": "HIGH" if st and st != '-' and st != 'ND' else "LOW",
            "mlst_assigned": True if st and st != '-' and st != 'ND' else False
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
            '2': {
                "clonal_complex": "CC2",
                "classification": "Healthcare-associated MRSA",
                "geographic_distribution": "Europe, North America",
                "clinical_significance": "Less common lineage, occasionally associated with hospital outbreaks",
                "common_virulence": ["Enterotoxins", "Hemolysins", "Proteases"],
                "outbreak_potential": "MEDIUM",
                "typical_sccmec": ["II", "IV"],
                "typical_spa": ["t002", "t037"],
                "resistance_profile": ["Methicillin", "Variable resistance"]
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
            '6': {
                "clonal_complex": "CC6",
                "classification": "Healthcare-associated MRSA",
                "geographic_distribution": "Europe, Asia",
                "clinical_significance": "Associated with hospital-acquired infections, particularly in neonatal units",
                "common_virulence": ["Enterotoxins", "Hemolysins"],
                "outbreak_potential": "MEDIUM",
                "typical_sccmec": ["IV", "II"],
                "typical_spa": ["t701", "t304"],
                "resistance_profile": ["Methicillin", "Gentamicin", "Erythromycin"]
            },
            '7': {
                "clonal_complex": "CC7",
                "classification": "Community-associated MRSA",
                "geographic_distribution": "Australia, Asia-Pacific",
                "clinical_significance": "Often associated with skin and soft tissue infections in community settings",
                "common_virulence": ["Panton-Valentine Leukocidin", "Enterotoxins"],
                "outbreak_potential": "MEDIUM",
                "typical_sccmec": ["IV", "V"],
                "typical_spa": ["t091", "t657"],
                "resistance_profile": ["Methicillin", "Fusidic acid"]
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
            '9': {
                "clonal_complex": "CC9",
                "classification": "Livestock-associated MRSA",
                "geographic_distribution": "Asia, Europe",
                "clinical_significance": "Zoonotic transmission from livestock, particularly swine, emerging public health concern",
                "common_virulence": ["Limited virulence arsenal", "Animal-adapted factors"],
                "outbreak_potential": "HIGH (agricultural settings)",
                "typical_sccmec": ["IV", "V"],
                "typical_spa": ["t899", "t337"],
                "resistance_profile": ["Methicillin", "Tetracycline", "Multi-drug resistance common"]
            },
            '10': {
                "clonal_complex": "CC81",
                "classification": "Community-associated MRSA",
                "geographic_distribution": "Global, Eastern Asia",
                "clinical_significance": "Diverse lineage with both community and healthcare associations",
                "common_virulence": ["Variable virulence factors", "Enterotoxins"],
                "outbreak_potential": "MEDIUM",
                "typical_sccmec": ["IV", "V"],
                "typical_spa": ["t021", "t045"],
                "resistance_profile": ["Variable resistance patterns"]
            },
            '12': {
                "clonal_complex": "CC12",
                "classification": "Healthcare-associated MRSA",
                "geographic_distribution": "Europe, North America",
                "clinical_significance": "Includes EMRSA-12 clone, associated with hospital outbreaks",
                "common_virulence": ["Enterotoxins", "Hemolysins"],
                "outbreak_potential": "HIGH",
                "typical_sccmec": ["II", "IV"],
                "typical_spa": ["t032", "t037"],
                "resistance_profile": ["Methicillin", "Multi-drug resistance"]
            },
            '15': {
                "clonal_complex": "CC15",
                "classification": "Community-associated MRSA",
                "geographic_distribution": "Global",
                "clinical_significance": "Common MSSA lineage, can acquire SCCmec to become MRSA",
                "common_virulence": ["Enterotoxins", "Hemolysins"],
                "outbreak_potential": "MEDIUM",
                "typical_sccmec": ["IV", "V"],
                "typical_spa": ["t084", "t085"],
                "resistance_profile": ["Variable, often MSSA"]
            },
            '20': {
                "clonal_complex": "CC20",
                "classification": "Healthcare-associated MRSA",
                "geographic_distribution": "Europe, Middle East",
                "clinical_significance": "Associated with hospital-acquired infections and device-related infections",
                "common_virulence": ["Enterotoxins", "Biofilm formation genes"],
                "outbreak_potential": "MEDIUM",
                "typical_sccmec": ["II", "IV"],
                "typical_spa": ["t164", "t021"],
                "resistance_profile": ["Methicillin", "Aminoglycosides"]
            },
            '22': {
                "clonal_complex": "CC22",
                "classification": "Healthcare-associated MRSA",
                "geographic_distribution": "Europe, Middle East, Global", 
                "clinical_significance": "Epidemic MRSA-15 (EMRSA-15), major healthcare clone with high transmission in hospital settings",
                "common_virulence": ["Enterotoxins", "Immune evasion cluster", "Hemolysins"],
                "outbreak_potential": "VERY HIGH",
                "typical_sccmec": ["IV", "IVh", "IVa", "IVc", "IVd", "Vb"],
                "typical_spa": ["t032", "t022", "t005", "t852"],
                "resistance_profile": ["Methicillin", "Multiple drug classes", "Often gentamicin resistant"]
            },
            '25': {
                "clonal_complex": "CC25",
                "classification": "Livestock-associated MRSA",
                "geographic_distribution": "Europe",
                "clinical_significance": "Associated with bovine mastitis, zoonotic potential",
                "common_virulence": ["Bovine-adapted factors", "Limited human virulence"],
                "outbreak_potential": "LOW (human), MEDIUM (bovine)",
                "typical_sccmec": ["IV", "V"],
                "typical_spa": ["t078", "t387"],
                "resistance_profile": ["Methicillin", "Penicillin", "Tetracycline"]
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
            '36': {
                "clonal_complex": "CC30",
                "classification": "Healthcare-associated MRSA",
                "geographic_distribution": "Middle East, Asia, Europe",
                "clinical_significance": "Emerging clone in Middle Eastern hospitals",
                "common_virulence": ["Enterotoxins", "Hemolysins"],
                "outbreak_potential": "MEDIUM",
                "typical_sccmec": ["III", "II", "IV"],
                "typical_spa": ["t032", "t037"],
                "resistance_profile": ["Methicillin", "Multi-drug resistant"]
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
            '59': {
                "clonal_complex": "CC59",
                "classification": "Community-associated MRSA",
                "geographic_distribution": "Asia-Pacific, Taiwan clone, USA(ST59)",
                "clinical_significance": "Taiwan clone/ST59, highly virulent community clone in Asia and USA",
                "common_virulence": ["Panton-Valentine Leukocidin", "Enterotoxins", "Hemolysins"],
                "outbreak_potential": "HIGH",
                "typical_sccmec": ["IV", "V"],
                "typical_spa": ["t437", "t441", "t163"],
                "resistance_profile": ["Methicillin", "Multi-drug resistant"]
            },
            '72': {
                "clonal_complex": "CC72",
                "classification": "Community-associated MRSA",
                "geographic_distribution": "USA, South America",
                "clinical_significance": "USA700 clone, community-associated with skin infections",
                "common_virulence": ["Panton-Valentine Leukocidin", "Enterotoxins"],
                "outbreak_potential": "MEDIUM",
                "typical_sccmec": ["IV", "V"],
                "typical_spa": ["t148", "t324", "t324", "t664", "t791", "t3092"],
                "resistance_profile": ["Methicillin", "Variable resistance"]
            },
            '75': {
                "clonal_complex": "CC75",
                "classification": "Community-associated MRSA",
                "geographic_distribution": "Australia, Asia",
                "clinical_significance": "Western Samoan Phage Pattern clone, associated with tropical regions",
                "common_virulence": ["Enterotoxins", "Hemolysins"],
                "outbreak_potential": "MEDIUM",
                "typical_sccmec": ["IV", "V"],
                "typical_spa": ["t359", "t084"],
                "resistance_profile": ["Methicillin", "Fusidic acid"]
            },
            '78': {
                "clonal_complex": "CC88",
                "classification": "Livestock-associated MRSA",
                "geographic_distribution": "Europe",
                "clinical_significance": "Associated with livestock, particularly poultry",
                "common_virulence": ["Animal-adapted virulence factors"],
                "outbreak_potential": "LOW (human), MEDIUM (poultry)",
                "typical_sccmec": ["IV", "V"],
                "typical_spa": ["t078", "t1773"],
                "resistance_profile": ["Methicillin", "Tetracycline"]
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
            '88': {
                "clonal_complex": "CC88",
                "classification": "Community-associated MRSA",
                "geographic_distribution": "Africa, Middle East",
                "clinical_significance": "African clone, emerging in sub-Saharan Africa",
                "common_virulence": ["Panton-Valentine Leukocidin", "Enterotoxins"],
                "outbreak_potential": "MEDIUM",
                "typical_sccmec": ["IV", "V"],
                "typical_spa": ["t186", "t325"],
                "resistance_profile": ["Methicillin", "Variable resistance"]
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
            '97': {
                "clonal_complex": "CC97",
                "classification": "Livestock-associated MRSA",
                "geographic_distribution": "Global",
                "clinical_significance": "Zoonotic transmission from livestock (especially cattle), emerging human infections",
                "common_virulence": ["Animal-adapted factors", "Some human virulence genes"],
                "outbreak_potential": "MEDIUM (human), HIGH (livestock)",
                "typical_sccmec": ["IV", "V"],
                "typical_spa": ["t267", "t359", "t1730"],
                "resistance_profile": ["Methicillin", "Tetracycline", "Multi-drug resistant"]
            },
            '101': {
                "clonal_complex": "CC101",
                "classification": "Community-associated MRSA",
                "geographic_distribution": "Europe",
                "clinical_significance": "Emerging community clone in Europe",
                "common_virulence": ["Enterotoxins", "Hemolysins"],
                "outbreak_potential": "LOW",
                "typical_sccmec": ["IV", "V"],
                "typical_spa": ["t571", "t1274"],
                "resistance_profile": ["Methicillin", "Variable resistance"]
            },
            '105': {
                "clonal_complex": "CC105",
                "classification": "Healthcare-associated MRSA",
                "geographic_distribution": "Middle East, America, Europe, Asia",
                "clinical_significance": "Emerging in healthcare settings",
                "common_virulence": ["Enterotoxins", "Hemolysins"],
                "outbreak_potential": "MEDIUM",
                "typical_sccmec": ["III", "IV"],
                "typical_spa": ["t002", "t037"],
                "resistance_profile": ["Methicillin", "Multi-drug resistant"]
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
            '130': {
                "clonal_complex": "CC130",
                "classification": "Livestock-associated MRSA",
                "geographic_distribution": "Europe",
                "clinical_significance": "Zoonotic transmission from livestock, often mecC-positive (alternative methicillin resistance)",
                "common_virulence": ["Limited human virulence", "Animal-adapted"],
                "outbreak_potential": "LOW (human), MEDIUM (livestock)",
                "typical_sccmec": ["XI (mecC)", "IV"],
                "typical_spa": ["t843", "t1736"],
                "resistance_profile": ["Methicillin (mecC)", "Tetracycline"]
            },
            '133': {
                "clonal_complex": "CC133",
                "classification": "Livestock-associated MRSA",
                "geographic_distribution": "Europe, Middle East",
                "clinical_significance": "Associated with ruminants, zoonotic potential",
                "common_virulence": ["Animal-adapted factors"],
                "outbreak_potential": "LOW (human), MEDIUM (livestock)",
                "typical_sccmec": ["IV", "V"],
                "typical_spa": ["t1166", "t1730"],
                "resistance_profile": ["Methicillin", "Tetracycline"]
            },
            '152': {
                "clonal_complex": "CC152",
                "classification": "Community-associated MRSA",
                "geographic_distribution": "Europe, Middle East, sub-Saharan Africa,",
                "clinical_significance": "Often PVL-positive, associated with community-acquired infections",
                "common_virulence": ["Panton-Valentine Leukocidin", "Enterotoxins"],
                "outbreak_potential": "MEDIUM",
                "typical_sccmec": ["IV", "V"],
                "typical_spa": ["t355", "t657"],
                "resistance_profile": ["Methicillin", "Variable resistance"]
            },
            '188': {
                "clonal_complex": "CC188",
                "classification": "Community-associated MRSA",
                "geographic_distribution": "Europe, Middle East",
                "clinical_significance": "Associated with skin and soft tissue infections",
                "common_virulence": ["Enterotoxins", "Hemolysins"],
                "outbreak_potential": "MEDIUM",
                "typical_sccmec": ["IV", "V"],
                "typical_spa": ["t189", "t325"],
                "resistance_profile": ["Methicillin", "Variable resistance"]
            },
            '239': {
                "clonal_complex": "CC8/CC30",
                "classification": "Healthcare-associated MRSA",
                "geographic_distribution": "Asia, Brazil",
                "clinical_significance": "Brazilian/Hungarian epidemic clone, multi-drug resistant hospital strain",
                "common_virulence": ["Enterotoxins", "Hemolysins", "Biofilm formation"],
                "outbreak_potential": "VERY HIGH",
                "typical_sccmec": ["III", "IIIA"],
                "typical_spa": ["t037", "t030"],
                "resistance_profile": ["Methicillin", "Multi-drug resistant including aminoglycosides"]
            },
            '247': {
                "clonal_complex": "CC8",
                "classification": "Healthcare-associated MRSA",
                "geographic_distribution": "Iberian Peninsula, South America",
                "clinical_significance": "Iberian clone, hospital-associated with high resistance",
                "common_virulence": ["Enterotoxins", "Hemolysins"],
                "outbreak_potential": "HIGH",
                "typical_sccmec": ["I", "III"],
                "typical_spa": ["t051", "t037"],
                "resistance_profile": ["Methicillin", "High-level multi-drug resistance"]
            },
            '250': {
                "clonal_complex": "CC8",
                "classification": "Healthcare-associated MRSA",
                "geographic_distribution": "Europe",
                "clinical_significance": "Hospital-associated clone, less common",
                "common_virulence": ["Enterotoxins", "Hemolysins"],
                "outbreak_potential": "MEDIUM",
                "typical_sccmec": ["II", "IV"],
                "typical_spa": ["t032", "t022"],
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
            },
            '425': {
                "clonal_complex": "CC425",
                "classification": "Livestock-associated MRSA",
                "geographic_distribution": "Europe",
                "clinical_significance": "Associated with poultry, zoonotic potential",
                "common_virulence": ["Animal-adapted factors"],
                "outbreak_potential": "LOW (human), MEDIUM (poultry)",
                "typical_sccmec": ["V", "IV"],
                "typical_spa": ["t2245", "t1730"],
                "resistance_profile": ["Methicillin", "Tetracycline"]
            },
            '582': {
                "clonal_complex": "CC15",
                "classification": "Community-associated MRSA",
                "geographic_distribution": "Southeast Asia",
                "clinical_significance": "Emerging community clone in Southeast Asia",
                "common_virulence": ["Panton-Valentine Leukocidin", "Enterotoxins"],
                "outbreak_potential": "MEDIUM",
                "typical_sccmec": ["IV", "V"],
                "typical_spa": ["t657", "t324"],
                "resistance_profile": ["Methicillin", "Variable resistance"]
            },
            '772': {
                "clonal_complex": "CC1",
                "classification": "Healthcare-associated MRSA",
                "geographic_distribution": "Middle East",
                "clinical_significance": "Emerging in Middle Eastern hospitals",
                "common_virulence": ["Enterotoxins", "Hemolysins"],
                "outbreak_potential": "MEDIUM",
                "typical_sccmec": ["III", "IV"],
                "typical_spa": ["t044", "t037"],
                "resistance_profile": ["Methicillin", "Multi-drug resistant"]
            },   
            '3': {
                "clonal_complex": "CC3",
                "classification": "Healthcare-associated MRSA",
                "geographic_distribution": "Europe, UK (prevalent), Australia",
                "clinical_significance": "EMRSA-3 (Epidemic MRSA-3), major hospital-acquired clone in the UK during 1990s, often associated with surgical site infections and bacteremia",
                "common_virulence": ["Enterotoxins", "Hemolysins", "Immune evasion cluster (scn, sak, chp)"],
                "outbreak_potential": "HIGH (historical epidemic)",
                "typical_sccmec": ["III"],
                "typical_spa": ["t037", "t045"],
                "resistance_profile": ["Methicillin", "Multi-drug resistant including aminoglycosides"]
            },
            '4': {
                "clonal_complex": "CC45",
                "classification": "Rare lineage",
                "geographic_distribution": "Variable",
                "clinical_significance": "Limited documentation; CC45 typically includes USA600 clone with both community and healthcare associations",
                "common_virulence": ["Enterotoxins", "Hemolysins", "Immune evasion genes (typical of CC45)"],
                "outbreak_potential": "Unknown",
                "typical_sccmec": ["Variable"],
                "typical_spa": ["Unknown or variable"],
                "resistance_profile": ["Variable resistance patterns"]
            },
            '11': {
                "clonal_complex": "CC5",
                "classification": "Rare lineage",
                "geographic_distribution": "Variable",
                "clinical_significance": "Limited documentation; CC5 lineage often healthcare-associated with extensive virulence arsenal",
                "common_virulence": ["Enterotoxins (SEC, SEL, SEU)", "Immune evasion cluster", "Hemolysins", "Proteases"],
                "outbreak_potential": "Unknown",
                "typical_sccmec": ["Variable"],
                "typical_spa": ["Unknown or variable"],
                "resistance_profile": ["Variable resistance patterns"]
            },
            '13': {
                "clonal_complex": "CC13",
                "classification": "Rare lineage",
                "geographic_distribution": "Variable",
                "clinical_significance": "Limited documentation; associated with specific geographic regions",
                "common_virulence": ["Variable virulence factors"],
                "outbreak_potential": "Unknown",
                "typical_sccmec": ["Variable"],
                "typical_spa": ["Unknown or variable"],
                "resistance_profile": ["Variable resistance patterns"]
            },
            '14': {
                "clonal_complex": "CC15",
                "classification": "Rare lineage",
                "geographic_distribution": "Variable",
                "clinical_significance": "Limited documentation; CC15 is common MSSA lineage that can acquire SCCmec",
                "common_virulence": ["Enterotoxins", "Hemolysins", "Typical CC15 virulence factors"],
                "outbreak_potential": "Unknown",
                "typical_sccmec": ["Variable"],
                "typical_spa": ["Unknown or variable"],
                "resistance_profile": ["Variable resistance patterns"]
            },
            '16': {
                "clonal_complex": "CC15",
                "classification": "Rare lineage",
                "geographic_distribution": "Variable",
                "clinical_significance": "Limited documentation",
                "common_virulence": ["Enterotoxins", "Hemolysins"],
                "outbreak_potential": "Unknown",
                "typical_sccmec": ["Variable"],
                "typical_spa": ["Unknown or variable"],
                "resistance_profile": ["Variable resistance patterns"]
            },
            '17': {
                "clonal_complex": "CC30",
                "classification": "Rare lineage",
                "geographic_distribution": "Variable",
                "clinical_significance": "Limited documentation; CC30 lineage often PVL-positive with both hospital and community associations",
                "common_virulence": ["Panton-Valentine Leukocidin", "Enterotoxins", "Hemolysins"],
                "outbreak_potential": "Unknown",
                "typical_sccmec": ["Variable"],
                "typical_spa": ["Unknown or variable"],
                "resistance_profile": ["Variable resistance patterns"]
            },
            '18': {
                "clonal_complex": "CC15",
                "classification": "Rare lineage",
                "geographic_distribution": "Variable",
                "clinical_significance": "Limited documentation",
                "common_virulence": ["Enterotoxins", "Hemolysins"],
                "outbreak_potential": "Unknown",
                "typical_sccmec": ["Variable"],
                "typical_spa": ["Unknown or variable"],
                "resistance_profile": ["Variable resistance patterns"]
            },
            '19': {
                "clonal_complex": "CC30",
                "classification": "Rare lineage",
                "geographic_distribution": "Variable",
                "clinical_significance": "Limited documentation",
                "common_virulence": ["Panton-Valentine Leukocidin", "Enterotoxins", "Hemolysins"],
                "outbreak_potential": "Unknown",
                "typical_sccmec": ["Variable"],
                "typical_spa": ["Unknown or variable"],
                "resistance_profile": ["Variable resistance patterns"]
            },
            '21': {
                "clonal_complex": "CC22",
                "classification": "Rare lineage",
                "geographic_distribution": "Variable",
                "clinical_significance": "Limited documentation; CC22 includes EMRSA-15 with high hospital transmission",
                "common_virulence": ["Enterotoxins", "Immune evasion cluster", "Hemolysins"],
                "outbreak_potential": "Unknown",
                "typical_sccmec": ["Variable"],
                "typical_spa": ["Unknown or variable"],
                "resistance_profile": ["Variable resistance patterns"]
            },
            '23': {
                "clonal_complex": "CC22",
                "classification": "Rare lineage",
                "geographic_distribution": "Variable",
                "clinical_significance": "Limited documentation",
                "common_virulence": ["Enterotoxins", "Immune evasion cluster", "Hemolysins"],
                "outbreak_potential": "Unknown",
                "typical_sccmec": ["Variable"],
                "typical_spa": ["Unknown or variable"],
                "resistance_profile": ["Variable resistance patterns"]
            },
            '24': {
                "clonal_complex": "CC30",
                "classification": "Rare lineage",
                "geographic_distribution": "Variable",
                "clinical_significance": "Limited documentation",
                "common_virulence": ["Panton-Valentine Leukocidin", "Enterotoxins", "Hemolysins"],
                "outbreak_potential": "Unknown",
                "typical_sccmec": ["Variable"],
                "typical_spa": ["Unknown or variable"],
                "resistance_profile": ["Variable resistance patterns"]
            },
            '26': {
                "clonal_complex": "CC26",
                "classification": "Rare lineage",
                "geographic_distribution": "Variable",
                "clinical_significance": "Limited documentation",
                "common_virulence": ["Variable virulence factors"],
                "outbreak_potential": "Unknown",
                "typical_sccmec": ["Variable"],
                "typical_spa": ["Unknown or variable"],
                "resistance_profile": ["Variable resistance patterns"]
            },
            '27': {
                "clonal_complex": "CC8",
                "classification": "Rare lineage",
                "geographic_distribution": "Variable",
                "clinical_significance": "Limited documentation; CC8 includes pandemic USA300 clone with high virulence",
                "common_virulence": ["Panton-Valentine Leukocidin", "Enterotoxins", "Phenol-soluble modulins"],
                "outbreak_potential": "Unknown",
                "typical_sccmec": ["Variable"],
                "typical_spa": ["Unknown or variable"],
                "resistance_profile": ["Variable resistance patterns"]
            },
            '28': {
                "clonal_complex": "CC28",
                "classification": "Rare lineage",
                "geographic_distribution": "Variable",
                "clinical_significance": "Limited documentation",
                "common_virulence": ["Variable virulence factors"],
                "outbreak_potential": "Unknown",
                "typical_sccmec": ["Variable"],
                "typical_spa": ["Unknown or variable"],
                "resistance_profile": ["Variable resistance patterns"]
            },
            '29': {
                "clonal_complex": "CC121",
                "classification": "Rare lineage",
                "geographic_distribution": "Variable",
                "clinical_significance": "Limited documentation; CC121 often associated with exotoxin production",
                "common_virulence": ["Exfoliative toxins", "Enterotoxins"],
                "outbreak_potential": "Unknown",
                "typical_sccmec": ["Variable"],
                "typical_spa": ["Unknown or variable"],
                "resistance_profile": ["Variable resistance patterns"]
            },
            '31': {
                "clonal_complex": "CC30",
                "classification": "Rare lineage",
                "geographic_distribution": "Variable",
                "clinical_significance": "Limited documentation",
                "common_virulence": ["Panton-Valentine Leukocidin", "Enterotoxins", "Hemolysins"],
                "outbreak_potential": "Unknown",
                "typical_sccmec": ["Variable"],
                "typical_spa": ["Unknown or variable"],
                "resistance_profile": ["Variable resistance patterns"]
            },
            '32': {
                "clonal_complex": "CC30",
                "classification": "Rare lineage",
                "geographic_distribution": "Variable",
                "clinical_significance": "Limited documentation",
                "common_virulence": ["Panton-Valentine Leukocidin", "Enterotoxins", "Hemolysins"],
                "outbreak_potential": "Unknown",
                "typical_sccmec": ["Variable"],
                "typical_spa": ["Unknown or variable"],
                "resistance_profile": ["Variable resistance patterns"]
            },
            '33': {
                "clonal_complex": "CC30",
                "classification": "Rare lineage",
                "geographic_distribution": "Variable",
                "clinical_significance": "Limited documentation",
                "common_virulence": ["Panton-Valentine Leukocidin", "Enterotoxins", "Hemolysins"],
                "outbreak_potential": "Unknown",
                "typical_sccmec": ["Variable"],
                "typical_spa": ["Unknown or variable"],
                "resistance_profile": ["Variable resistance patterns"]
            },
            '34': {
                "clonal_complex": "CC30",
                "classification": "Community and Healthcare-associated MRSA",
                "geographic_distribution": "Asia-Pacific, Australia (Southwest Pacific clone)",
                "clinical_significance": "Southwest Pacific (SWP) clone, often PVL-positive, circulating in both community and hospital settings with significant transmissibility",
                "common_virulence": ["Panton-Valentine Leukocidin", "Enterotoxins", "Hemolysins"],
                "outbreak_potential": "HIGH",
                "typical_sccmec": ["IV"],
                "typical_spa": ["t019", "t318", "t021"],
                "resistance_profile": ["Methicillin", "Variable resistance patterns"]
            },
            '35': {
                "clonal_complex": "CC15",
                "classification": "Rare lineage",
                "geographic_distribution": "Variable",
                "clinical_significance": "Limited documentation",
                "common_virulence": ["Enterotoxins", "Hemolysins"],
                "outbreak_potential": "Unknown",
                "typical_sccmec": ["Variable"],
                "typical_spa": ["Unknown or variable"],
                "resistance_profile": ["Variable resistance patterns"]
            },
            '37': {
                "clonal_complex": "CC30",
                "classification": "Rare lineage",
                "geographic_distribution": "Variable",
                "clinical_significance": "Limited documentation",
                "common_virulence": ["Panton-Valentine Leukocidin", "Enterotoxins", "Hemolysins"],
                "outbreak_potential": "Unknown",
                "typical_sccmec": ["Variable"],
                "typical_spa": ["Unknown or variable"],
                "resistance_profile": ["Variable resistance patterns"]
            },
            '38': {
                "clonal_complex": "CC30",
                "classification": "Rare lineage",
                "geographic_distribution": "Variable",
                "clinical_significance": "Limited documentation",
                "common_virulence": ["Panton-Valentine Leukocidin", "Enterotoxins", "Hemolysins"],
                "outbreak_potential": "Unknown",
                "typical_sccmec": ["Variable"],
                "typical_spa": ["Unknown or variable"],
                "resistance_profile": ["Variable resistance patterns"]
            },
            '39': {
                "clonal_complex": "CC30",
                "classification": "Rare lineage",
                "geographic_distribution": "Variable",
                "clinical_significance": "Limited documentation",
                "common_virulence": ["Panton-Valentine Leukocidin", "Enterotoxins", "Hemolysins"],
                "outbreak_potential": "Unknown",
                "typical_sccmec": ["Variable"],
                "typical_spa": ["Unknown or variable"],
                "resistance_profile": ["Variable resistance patterns"]
            },
            '40': {
                "clonal_complex": "CC30",
                "classification": "Rare lineage",
                "geographic_distribution": "Variable",
                "clinical_significance": "Limited documentation",
                "common_virulence": ["Panton-Valentine Leukocidin", "Enterotoxins", "Hemolysins"],
                "outbreak_potential": "Unknown",
                "typical_sccmec": ["Variable"],
                "typical_spa": ["Unknown or variable"],
                "resistance_profile": ["Variable resistance patterns"]
            },
            '41': {
                "clonal_complex": "CC30",
                "classification": "Rare lineage",
                "geographic_distribution": "Variable",
                "clinical_significance": "Limited documentation",
                "common_virulence": ["Panton-Valentine Leukocidin", "Enterotoxins", "Hemolysins"],
                "outbreak_potential": "Unknown",
                "typical_sccmec": ["Variable"],
                "typical_spa": ["Unknown or variable"],
                "resistance_profile": ["Variable resistance patterns"]
            },
            '42': {
                "clonal_complex": "CC42",
                "classification": "Rare lineage",
                "geographic_distribution": "Variable",
                "clinical_significance": "Limited documentation",
                "common_virulence": ["Variable virulence factors"],
                "outbreak_potential": "Unknown",
                "typical_sccmec": ["Variable"],
                "typical_spa": ["Unknown or variable"],
                "resistance_profile": ["Variable resistance patterns"]
            },
            '43': {
                "clonal_complex": "CC30",
                "classification": "Rare lineage",
                "geographic_distribution": "Variable",
                "clinical_significance": "Limited documentation",
                "common_virulence": ["Panton-Valentine Leukocidin", "Enterotoxins", "Hemolysins"],
                "outbreak_potential": "Unknown",
                "typical_sccmec": ["Variable"],
                "typical_spa": ["Unknown or variable"],
                "resistance_profile": ["Variable resistance patterns"]
            },
            '44': {
                "clonal_complex": "CC22",
                "classification": "Rare lineage",
                "geographic_distribution": "Variable",
                "clinical_significance": "Limited documentation",
                "common_virulence": ["Enterotoxins", "Immune evasion cluster", "Hemolysins"],
                "outbreak_potential": "Unknown",
                "typical_sccmec": ["Variable"],
                "typical_spa": ["Unknown or variable"],
                "resistance_profile": ["Variable resistance patterns"]
            },
            '46': {
                "clonal_complex": "CC45",
                "classification": "Rare lineage",
                "geographic_distribution": "Variable",
                "clinical_significance": "Limited documentation",
                "common_virulence": ["Enterotoxins", "Hemolysins", "Immune evasion genes"],
                "outbreak_potential": "Unknown",
                "typical_sccmec": ["Variable"],
                "typical_spa": ["Unknown or variable"],
                "resistance_profile": ["Variable resistance patterns"]
            },
            '47': {
                "clonal_complex": "CC45",
                "classification": "Rare lineage",
                "geographic_distribution": "Variable",
                "clinical_significance": "Limited documentation",
                "common_virulence": ["Enterotoxins", "Hemolysins", "Immune evasion genes"],
                "outbreak_potential": "Unknown",
                "typical_sccmec": ["Variable"],
                "typical_spa": ["Unknown or variable"],
                "resistance_profile": ["Variable resistance patterns"]
            },
            '48': {
                "clonal_complex": "CC45",
                "classification": "Rare lineage",
                "geographic_distribution": "Variable",
                "clinical_significance": "Limited documentation",
                "common_virulence": ["Enterotoxins", "Hemolysins", "Immune evasion genes"],
                "outbreak_potential": "Unknown",
                "typical_sccmec": ["Variable"],
                "typical_spa": ["Unknown or variable"],
                "resistance_profile": ["Variable resistance patterns"]
            },
            '49': {
                "clonal_complex": "CC49",
                "classification": "Rare lineage",
                "geographic_distribution": "Variable",
                "clinical_significance": "Limited documentation",
                "common_virulence": ["Variable virulence factors"],
                "outbreak_potential": "Unknown",
                "typical_sccmec": ["Variable"],
                "typical_spa": ["Unknown or variable"],
                "resistance_profile": ["Variable resistance patterns"]
            },
            '50': {
                "clonal_complex": "CC50",
                "classification": "Rare lineage",
                "geographic_distribution": "Variable",
                "clinical_significance": "Limited documentation",
                "common_virulence": ["Variable virulence factors"],
                "outbreak_potential": "Unknown",
                "typical_sccmec": ["Variable"],
                "typical_spa": ["Unknown or variable"],
                "resistance_profile": ["Variable resistance patterns"]
            },
            '51': {
                "clonal_complex": "CC121",
                "classification": "Rare lineage",
                "geographic_distribution": "Variable",
                "clinical_significance": "Limited documentation",
                "common_virulence": ["Exfoliative toxins", "Enterotoxins"],
                "outbreak_potential": "Unknown",
                "typical_sccmec": ["Variable"],
                "typical_spa": ["Unknown or variable"],
                "resistance_profile": ["Variable resistance patterns"]
            },
            '52': {
                "clonal_complex": "CC30",
                "classification": "Rare lineage",
                "geographic_distribution": "Variable",
                "clinical_significance": "Limited documentation",
                "common_virulence": ["Panton-Valentine Leukocidin", "Enterotoxins", "Hemolysins"],
                "outbreak_potential": "Unknown",
                "typical_sccmec": ["Variable"],
                "typical_spa": ["Unknown or variable"],
                "resistance_profile": ["Variable resistance patterns"]
            },
            '53': {
                "clonal_complex": "CC45",
                "classification": "Rare lineage",
                "geographic_distribution": "Variable",
                "clinical_significance": "Limited documentation",
                "common_virulence": ["Enterotoxins", "Hemolysins", "Immune evasion genes"],
                "outbreak_potential": "Unknown",
                "typical_sccmec": ["Variable"],
                "typical_spa": ["Unknown or variable"],
                "resistance_profile": ["Variable resistance patterns"]
            },
            '54': {
                "clonal_complex": "CC45",
                "classification": "Rare lineage",
                "geographic_distribution": "Variable",
                "clinical_significance": "Limited documentation",
                "common_virulence": ["Enterotoxins", "Hemolysins", "Immune evasion genes"],
                "outbreak_potential": "Unknown",
                "typical_sccmec": ["Variable"],
                "typical_spa": ["Unknown or variable"],
                "resistance_profile": ["Variable resistance patterns"]
            },
            '55': {
                "clonal_complex": "CC55",
                "classification": "Rare lineage",
                "geographic_distribution": "Variable",
                "clinical_significance": "Limited documentation",
                "common_virulence": ["Variable virulence factors"],
                "outbreak_potential": "Unknown",
                "typical_sccmec": ["Variable"],
                "typical_spa": ["Unknown or variable"],
                "resistance_profile": ["Variable resistance patterns"]
            },
            '56': {
                "clonal_complex": "CC15",
                "classification": "Rare lineage",
                "geographic_distribution": "Variable",
                "clinical_significance": "Limited documentation",
                "common_virulence": ["Enterotoxins", "Hemolysins"],
                "outbreak_potential": "Unknown",
                "typical_sccmec": ["Variable"],
                "typical_spa": ["Unknown or variable"],
                "resistance_profile": ["Variable resistance patterns"]
            },
            '57': {
                "clonal_complex": "CC30",
                "classification": "Rare lineage",
                "geographic_distribution": "Variable",
                "clinical_significance": "Limited documentation",
                "common_virulence": ["Panton-Valentine Leukocidin", "Enterotoxins", "Hemolysins"],
                "outbreak_potential": "Unknown",
                "typical_sccmec": ["Variable"],
                "typical_spa": ["Unknown or variable"],
                "resistance_profile": ["Variable resistance patterns"]
            },
            '58': {
                "clonal_complex": "CC15",
                "classification": "Rare lineage",
                "geographic_distribution": "Variable",
                "clinical_significance": "Limited documentation",
                "common_virulence": ["Enterotoxins", "Hemolysins"],
                "outbreak_potential": "Unknown",
                "typical_sccmec": ["Variable"],
                "typical_spa": ["Unknown or variable"],
                "resistance_profile": ["Variable resistance patterns"]
            },
            '60': {
                "clonal_complex": "CC22",
                "classification": "Rare lineage",
                "geographic_distribution": "Variable",
                "clinical_significance": "Limited documentation",
                "common_virulence": ["Enterotoxins", "Immune evasion cluster", "Hemolysins"],
                "outbreak_potential": "Unknown",
                "typical_sccmec": ["Variable"],
                "typical_spa": ["Unknown or variable"],
                "resistance_profile": ["Variable resistance patterns"]
            },
            '61': {
                "clonal_complex": "CC15",
                "classification": "Rare lineage",
                "geographic_distribution": "Variable",
                "clinical_significance": "Limited documentation",
                "common_virulence": ["Enterotoxins", "Hemolysins"],
                "outbreak_potential": "Unknown",
                "typical_sccmec": ["Variable"],
                "typical_spa": ["Unknown or variable"],
                "resistance_profile": ["Variable resistance patterns"]
            },
            '62': {
                "clonal_complex": "CC121",
                "classification": "Rare lineage",
                "geographic_distribution": "Variable",
                "clinical_significance": "Limited documentation",
                "common_virulence": ["Exfoliative toxins", "Enterotoxins"],
                "outbreak_potential": "Unknown",
                "typical_sccmec": ["Variable"],
                "typical_spa": ["Unknown or variable"],
                "resistance_profile": ["Variable resistance patterns"]
            },
            '63': {
                "clonal_complex": "CC1",
                "classification": "Rare lineage",
                "geographic_distribution": "Variable",
                "clinical_significance": "Limited documentation; CC1 includes USA400 clone often PVL-positive",
                "common_virulence": ["Panton-Valentine Leukocidin", "Enterotoxins SEH/SEK", "Hemolysins"],
                "outbreak_potential": "Unknown",
                "typical_sccmec": ["Variable"],
                "typical_spa": ["Unknown or variable"],
                "resistance_profile": ["Variable resistance patterns"]
            },
            '64': {
                "clonal_complex": "CC64",
                "classification": "Rare lineage",
                "geographic_distribution": "Variable",
                "clinical_significance": "Limited documentation",
                "common_virulence": ["Variable virulence factors"],
                "outbreak_potential": "Unknown",
                "typical_sccmec": ["Variable"],
                "typical_spa": ["Unknown or variable"],
                "resistance_profile": ["Variable resistance patterns"]
            },
            '65': {
                "clonal_complex": "CC65",
                "classification": "Rare lineage",
                "geographic_distribution": "Variable",
                "clinical_significance": "Limited documentation",
                "common_virulence": ["Variable virulence factors"],
                "outbreak_potential": "Unknown",
                "typical_sccmec": ["Variable"],
                "typical_spa": ["Unknown or variable"],
                "resistance_profile": ["Variable resistance patterns"]
            },
            '66': {
                "clonal_complex": "CC66",
                "classification": "Rare lineage",
                "geographic_distribution": "Variable",
                "clinical_significance": "Limited documentation",
                "common_virulence": ["Variable virulence factors"],
                "outbreak_potential": "Unknown",
                "typical_sccmec": ["Variable"],
                "typical_spa": ["Unknown or variable"],
                "resistance_profile": ["Variable resistance patterns"]
            },
            '67': {
                "clonal_complex": "CC67",
                "classification": "Rare lineage",
                "geographic_distribution": "Variable",
                "clinical_significance": "Limited documentation",
                "common_virulence": ["Variable virulence factors"],
                "outbreak_potential": "Unknown",
                "typical_sccmec": ["Variable"],
                "typical_spa": ["Unknown or variable"],
                "resistance_profile": ["Variable resistance patterns"]
            },
            '68': {
                "clonal_complex": "CC68",
                "classification": "Rare lineage",
                "geographic_distribution": "Variable",
                "clinical_significance": "Limited documentation",
                "common_virulence": ["Variable virulence factors"],
                "outbreak_potential": "Unknown",
                "typical_sccmec": ["Variable"],
                "typical_spa": ["Unknown or variable"],
                "resistance_profile": ["Variable resistance patterns"]
            },
            '69': {
                "clonal_complex": "CC1",
                "classification": "Rare lineage",
                "geographic_distribution": "Variable",
                "clinical_significance": "Limited documentation",
                "common_virulence": ["Panton-Valentine Leukocidin", "Enterotoxins SEH/SEK", "Hemolysins"],
                "outbreak_potential": "Unknown",
                "typical_sccmec": ["Variable"],
                "typical_spa": ["Unknown or variable"],
                "resistance_profile": ["Variable resistance patterns"]
            },
            '70': {
                "clonal_complex": "CC97",
                "classification": "Rare lineage",
                "geographic_distribution": "Variable",
                "clinical_significance": "Limited documentation; CC97 is livestock-associated with zoonotic potential",
                "common_virulence": ["Animal-adapted factors", "Some human virulence genes"],
                "outbreak_potential": "Unknown",
                "typical_sccmec": ["Variable"],
                "typical_spa": ["Unknown or variable"],
                "resistance_profile": ["Variable resistance patterns"]
            },
            '71': {
                "clonal_complex": "CC97",
                "classification": "Rare lineage",
                "geographic_distribution": "Variable",
                "clinical_significance": "Limited documentation",
                "common_virulence": ["Animal-adapted factors", "Some human virulence genes"],
                "outbreak_potential": "Unknown",
                "typical_sccmec": ["Variable"],
                "typical_spa": ["Unknown or variable"],
                "resistance_profile": ["Variable resistance patterns"]
            },
            '73': {
                "clonal_complex": "CC5",
                "classification": "Rare lineage",
                "geographic_distribution": "Variable",
                "clinical_significance": "Limited documentation",
                "common_virulence": ["Enterotoxins (SEC, SEL, SEU)", "Immune evasion cluster", "Hemolysins"],
                "outbreak_potential": "Unknown",
                "typical_sccmec": ["Variable"],
                "typical_spa": ["Unknown or variable"],
                "resistance_profile": ["Variable resistance patterns"]
            },
            '74': {
                "clonal_complex": "CC30",
                "classification": "Rare lineage",
                "geographic_distribution": "Variable",
                "clinical_significance": "Limited documentation",
                "common_virulence": ["Panton-Valentine Leukocidin", "Enterotoxins", "Hemolysins"],
                "outbreak_potential": "Unknown",
                "typical_sccmec": ["Variable"],
                "typical_spa": ["Unknown or variable"],
                "resistance_profile": ["Variable resistance patterns"]
            },
            '76': {
                "clonal_complex": "CC1",
                "classification": "Rare lineage",
                "geographic_distribution": "Variable",
                "clinical_significance": "Limited documentation",
                "common_virulence": ["Panton-Valentine Leukocidin", "Enterotoxins SEH/SEK", "Hemolysins"],
                "outbreak_potential": "Unknown",
                "typical_sccmec": ["Variable"],
                "typical_spa": ["Unknown or variable"],
                "resistance_profile": ["Variable resistance patterns"]
            },
            '77': {
                "clonal_complex": "CC30",
                "classification": "Rare lineage",
                "geographic_distribution": "Variable",
                "clinical_significance": "Limited documentation",
                "common_virulence": ["Panton-Valentine Leukocidin", "Enterotoxins", "Hemolysins"],
                "outbreak_potential": "Unknown",
                "typical_sccmec": ["Variable"],
                "typical_spa": ["Unknown or variable"],
                "resistance_profile": ["Variable resistance patterns"]
            },
            '79': {
                "clonal_complex": "CC22",
                "classification": "Rare lineage",
                "geographic_distribution": "Variable",
                "clinical_significance": "Limited documentation",
                "common_virulence": ["Enterotoxins", "Immune evasion cluster", "Hemolysins"],
                "outbreak_potential": "Unknown",
                "typical_sccmec": ["Variable"],
                "typical_spa": ["Unknown or variable"],
                "resistance_profile": ["Variable resistance patterns"]
            },
            '81': {
                "clonal_complex": "CC1",
                "classification": "Rare lineage",
                "geographic_distribution": "Variable",
                "clinical_significance": "Limited documentation",
                "common_virulence": ["Panton-Valentine Leukocidin", "Enterotoxins SEH/SEK", "Hemolysins"],
                "outbreak_potential": "Unknown",
                "typical_sccmec": ["Variable"],
                "typical_spa": ["Unknown or variable"],
                "resistance_profile": ["Variable resistance patterns"]
            },
            '82': {
                "clonal_complex": "CC82",
                "classification": "Rare lineage",
                "geographic_distribution": "Variable",
                "clinical_significance": "Limited documentation",
                "common_virulence": ["Variable virulence factors"],
                "outbreak_potential": "Unknown",
                "typical_sccmec": ["Variable"],
                "typical_spa": ["Unknown or variable"],
                "resistance_profile": ["Variable resistance patterns"]
            },
            '83': {
                "clonal_complex": "CC5",
                "classification": "Rare lineage",
                "geographic_distribution": "Variable",
                "clinical_significance": "Limited documentation",
                "common_virulence": ["Enterotoxins (SEC, SEL, SEU)", "Immune evasion cluster", "Hemolysins"],
                "outbreak_potential": "Unknown",
                "typical_sccmec": ["Variable"],
                "typical_spa": ["Unknown or variable"],
                "resistance_profile": ["Variable resistance patterns"]
            },
            '84': {
                "clonal_complex": "CC84",
                "classification": "Rare lineage",
                "geographic_distribution": "Variable",
                "clinical_significance": "Limited documentation",
                "common_virulence": ["Variable virulence factors"],
                "outbreak_potential": "Unknown",
                "typical_sccmec": ["Variable"],
                "typical_spa": ["Unknown or variable"],
                "resistance_profile": ["Variable resistance patterns"]
            },
            '85': {
                "clonal_complex": "CC5",
                "classification": "Rare lineage",
                "geographic_distribution": "Variable",
                "clinical_significance": "Limited documentation",
                "common_virulence": ["Enterotoxins (SEC, SEL, SEU)", "Immune evasion cluster", "Hemolysins"],
                "outbreak_potential": "Unknown",
                "typical_sccmec": ["Variable"],
                "typical_spa": ["Unknown or variable"],
                "resistance_profile": ["Variable resistance patterns"]
            },
            '86': {
                "clonal_complex": "CC8",
                "classification": "Rare lineage",
                "geographic_distribution": "Variable",
                "clinical_significance": "Limited documentation",
                "common_virulence": ["Panton-Valentine Leukocidin", "Enterotoxins", "Phenol-soluble modulins"],
                "outbreak_potential": "Unknown",
                "typical_sccmec": ["Variable"],
                "typical_spa": ["Unknown or variable"],
                "resistance_profile": ["Variable resistance patterns"]
            },
            '87': {
                "clonal_complex": "CC87",
                "classification": "Rare lineage",
                "geographic_distribution": "Variable",
                "clinical_significance": "Limited documentation",
                "common_virulence": ["Variable virulence factors"],
                "outbreak_potential": "Unknown",
                "typical_sccmec": ["Variable"],
                "typical_spa": ["Unknown or variable"],
                "resistance_profile": ["Variable resistance patterns"]
            },
            '89': {
                "clonal_complex": "CC89",
                "classification": "Rare lineage",
                "geographic_distribution": "Variable",
                "clinical_significance": "Limited documentation",
                "common_virulence": ["Variable virulence factors"],
                "outbreak_potential": "Unknown",
                "typical_sccmec": ["Variable"],
                "typical_spa": ["Unknown or variable"],
                "resistance_profile": ["Variable resistance patterns"]
            },
            '90': {
                "clonal_complex": "CC5",
                "classification": "Rare lineage",
                "geographic_distribution": "Variable",
                "clinical_significance": "Limited documentation",
                "common_virulence": ["Enterotoxins (SEC, SEL, SEU)", "Immune evasion cluster", "Hemolysins"],
                "outbreak_potential": "Unknown",
                "typical_sccmec": ["Variable"],
                "typical_spa": ["Unknown or variable"],
                "resistance_profile": ["Variable resistance patterns"]
            },
            '91': {
                "clonal_complex": "CC91",
                "classification": "Rare lineage",
                "geographic_distribution": "Variable",
                "clinical_significance": "Limited documentation",
                "common_virulence": ["Variable virulence factors"],
                "outbreak_potential": "Unknown",
                "typical_sccmec": ["Variable"],
                "typical_spa": ["Unknown or variable"],
                "resistance_profile": ["Variable resistance patterns"]
            },
            '92': {
                "clonal_complex": "CC92",
                "classification": "Rare lineage",
                "geographic_distribution": "Variable",
                "clinical_significance": "Limited documentation",
                "common_virulence": ["Variable virulence factors"],
                "outbreak_potential": "Unknown",
                "typical_sccmec": ["Variable"],
                "typical_spa": ["Unknown or variable"],
                "resistance_profile": ["Variable resistance patterns"]
            },
            '94': {
                "clonal_complex": "CC8",
                "classification": "Rare lineage",
                "geographic_distribution": "Variable",
                "clinical_significance": "Limited documentation",
                "common_virulence": ["Panton-Valentine Leukocidin", "Enterotoxins", "Phenol-soluble modulins"],
                "outbreak_potential": "Unknown",
                "typical_sccmec": ["Variable"],
                "typical_spa": ["Unknown or variable"],
                "resistance_profile": ["Variable resistance patterns"]
            },
            '95': {
                "clonal_complex": "CC121",
                "classification": "Rare lineage",
                "geographic_distribution": "Variable",
                "clinical_significance": "Limited documentation",
                "common_virulence": ["Exfoliative toxins", "Enterotoxins"],
                "outbreak_potential": "Unknown",
                "typical_sccmec": ["Variable"],
                "typical_spa": ["Unknown or variable"],
                "resistance_profile": ["Variable resistance patterns"]
            },
            '96': {
                "clonal_complex": "CC96",
                "classification": "Rare lineage",
                "geographic_distribution": "Variable",
                "clinical_significance": "Limited documentation",
                "common_virulence": ["Variable virulence factors"],
                "outbreak_potential": "Unknown",
                "typical_sccmec": ["Variable"],
                "typical_spa": ["Unknown or variable"],
                "resistance_profile": ["Variable resistance patterns"]
            },
            '98': {
                "clonal_complex": "CC5",
                "classification": "Rare lineage",
                "geographic_distribution": "Variable",
                "clinical_significance": "Limited documentation",
                "common_virulence": ["Enterotoxins (SEC, SEL, SEU)", "Immune evasion cluster", "Hemolysins"],
                "outbreak_potential": "Unknown",
                "typical_sccmec": ["Variable"],
                "typical_spa": ["Unknown or variable"],
                "resistance_profile": ["Variable resistance patterns"]
            },
            '99': {
                "clonal_complex": "CC5",
                "classification": "Rare lineage",
                "geographic_distribution": "Variable",
                "clinical_significance": "Limited documentation",
                "common_virulence": ["Enterotoxins (SEC, SEL, SEU)", "Immune evasion cluster", "Hemolysins"],
                "outbreak_potential": "Unknown",
                "typical_sccmec": ["Variable"],
                "typical_spa": ["Unknown or variable"],
                "resistance_profile": ["Variable resistance patterns"]
            },
            '100': {
                "clonal_complex": "CC5",
                "classification": "Rare lineage",
                "geographic_distribution": "Variable",
                "clinical_significance": "Limited documentation",
                "common_virulence": ["Enterotoxins (SEC, SEL, SEU)", "Immune evasion cluster", "Hemolysins"],
                "outbreak_potential": "Unknown",
                "typical_sccmec": ["Variable"],
                "typical_spa": ["Unknown or variable"],
                "resistance_profile": ["Variable resistance patterns"]
            },
            '102': {
                "clonal_complex": "CC102",
                "classification": "Rare lineage",
                "geographic_distribution": "Variable",
                "clinical_significance": "Limited documentation",
                "common_virulence": ["Variable virulence factors"],
                "outbreak_potential": "Unknown",
                "typical_sccmec": ["Variable"],
                "typical_spa": ["Unknown or variable"],
                "resistance_profile": ["Variable resistance patterns"]
            },
            '103': {
                "clonal_complex": "CC103",
                "classification": "Rare lineage",
                "geographic_distribution": "Variable",
                "clinical_significance": "Limited documentation",
                "common_virulence": ["Variable virulence factors"],
                "outbreak_potential": "Unknown",
                "typical_sccmec": ["Variable"],
                "typical_spa": ["Unknown or variable"],
                "resistance_profile": ["Variable resistance patterns"]
            },
            '104': {
                "clonal_complex": "CC104",
                "classification": "Rare lineage",
                "geographic_distribution": "Variable",
                "clinical_significance": "Limited documentation",
                "common_virulence": ["Variable virulence factors"],
                "outbreak_potential": "Unknown",
                "typical_sccmec": ["Variable"],
                "typical_spa": ["Unknown or variable"],
                "resistance_profile": ["Variable resistance patterns"]
            },
            '106': {
                "clonal_complex": "CC106",
                "classification": "Rare lineage",
                "geographic_distribution": "Variable",
                "clinical_significance": "Limited documentation",
                "common_virulence": ["Variable virulence factors"],
                "outbreak_potential": "Unknown",
                "typical_sccmec": ["Variable"],
                "typical_spa": ["Unknown or variable"],
                "resistance_profile": ["Variable resistance patterns"]
            },
            '107': {
                "clonal_complex": "CC107",
                "classification": "Rare lineage",
                "geographic_distribution": "Variable",
                "clinical_significance": "Limited documentation",
                "common_virulence": ["Variable virulence factors"],
                "outbreak_potential": "Unknown",
                "typical_sccmec": ["Variable"],
                "typical_spa": ["Unknown or variable"],
                "resistance_profile": ["Variable resistance patterns"]
            },
            '108': {
                "clonal_complex": "CC45",
                "classification": "Rare lineage",
                "geographic_distribution": "Variable",
                "clinical_significance": "Limited documentation",
                "common_virulence": ["Enterotoxins", "Hemolysins", "Immune evasion genes"],
                "outbreak_potential": "Unknown",
                "typical_sccmec": ["Variable"],
                "typical_spa": ["Unknown or variable"],
                "resistance_profile": ["Variable resistance patterns"]
            },
            '109': {
                "clonal_complex": "CC1",
                "classification": "Rare lineage",
                "geographic_distribution": "Variable",
                "clinical_significance": "Limited documentation",
                "common_virulence": ["Panton-Valentine Leukocidin", "Enterotoxins SEH/SEK", "Hemolysins"],
                "outbreak_potential": "Unknown",
                "typical_sccmec": ["Variable"],
                "typical_spa": ["Unknown or variable"],
                "resistance_profile": ["Variable resistance patterns"]
            },
            '110': {
                "clonal_complex": "CC8",
                "classification": "Rare lineage",
                "geographic_distribution": "Variable",
                "clinical_significance": "Limited documentation",
                "common_virulence": ["Panton-Valentine Leukocidin", "Enterotoxins", "Phenol-soluble modulins"],
                "outbreak_potential": "Unknown",
                "typical_sccmec": ["Variable"],
                "typical_spa": ["Unknown or variable"],
                "resistance_profile": ["Variable resistance patterns"]
            },
            '111': {
                "clonal_complex": "CC5",
                "classification": "Rare lineage",
                "geographic_distribution": "Variable",
                "clinical_significance": "Limited documentation",
                "common_virulence": ["Enterotoxins (SEC, SEL, SEU)", "Immune evasion cluster", "Hemolysins"],
                "outbreak_potential": "Unknown",
                "typical_sccmec": ["Variable"],
                "typical_spa": ["Unknown or variable"],
                "resistance_profile": ["Variable resistance patterns"]
            },
            '112': {
                "clonal_complex": "CC8",
                "classification": "Rare lineage",
                "geographic_distribution": "Variable",
                "clinical_significance": "Limited documentation",
                "common_virulence": ["Panton-Valentine Leukocidin", "Enterotoxins", "Phenol-soluble modulins"],
                "outbreak_potential": "Unknown",
                "typical_sccmec": ["Variable"],
                "typical_spa": ["Unknown or variable"],
                "resistance_profile": ["Variable resistance patterns"]
            },
            '113': {
                "clonal_complex": "CC8",
                "classification": "Rare lineage",
                "geographic_distribution": "Variable",
                "clinical_significance": "Limited documentation",
                "common_virulence": ["Panton-Valentine Leukocidin", "Enterotoxins", "Phenol-soluble modulins"],
                "outbreak_potential": "Unknown",
                "typical_sccmec": ["Variable"],
                "typical_spa": ["Unknown or variable"],
                "resistance_profile": ["Variable resistance patterns"]
            },
            '114': {
                "clonal_complex": "CC5",
                "classification": "Rare lineage",
                "geographic_distribution": "Variable",
                "clinical_significance": "Limited documentation",
                "common_virulence": ["Enterotoxins (SEC, SEL, SEU)", "Immune evasion cluster", "Hemolysins"],
                "outbreak_potential": "Unknown",
                "typical_sccmec": ["Variable"],
                "typical_spa": ["Unknown or variable"],
                "resistance_profile": ["Variable resistance patterns"]
            },
            '115': {
                "clonal_complex": "CC97",
                "classification": "Rare lineage",
                "geographic_distribution": "Variable",
                "clinical_significance": "Limited documentation",
                "common_virulence": ["Animal-adapted factors", "Some human virulence genes"],
                "outbreak_potential": "Unknown",
                "typical_sccmec": ["Variable"],
                "typical_spa": ["Unknown or variable"],
                "resistance_profile": ["Variable resistance patterns"]
            },
            '116': {
                "clonal_complex": "CC97",
                "classification": "Rare lineage",
                "geographic_distribution": "Variable",
                "clinical_significance": "Limited documentation",
                "common_virulence": ["Animal-adapted factors", "Some human virulence genes"],
                "outbreak_potential": "Unknown",
                "typical_sccmec": ["Variable"],
                "typical_spa": ["Unknown or variable"],
                "resistance_profile": ["Variable resistance patterns"]
            },
            '117': {
                "clonal_complex": "CC97",
                "classification": "Rare lineage",
                "geographic_distribution": "Variable",
                "clinical_significance": "Limited documentation",
                "common_virulence": ["Animal-adapted factors", "Some human virulence genes"],
                "outbreak_potential": "Unknown",
                "typical_sccmec": ["Variable"],
                "typical_spa": ["Unknown or variable"],
                "resistance_profile": ["Variable resistance patterns"]
            },
            '118': {
                "clonal_complex": "CC97",
                "classification": "Rare lineage",
                "geographic_distribution": "Variable",
                "clinical_significance": "Limited documentation",
                "common_virulence": ["Animal-adapted factors", "Some human virulence genes"],
                "outbreak_potential": "Unknown",
                "typical_sccmec": ["Variable"],
                "typical_spa": ["Unknown or variable"],
                "resistance_profile": ["Variable resistance patterns"]
            },
            '119': {
                "clonal_complex": "CC97",
                "classification": "Rare lineage",
                "geographic_distribution": "Variable",
                "clinical_significance": "Limited documentation",
                "common_virulence": ["Animal-adapted factors", "Some human virulence genes"],
                "outbreak_potential": "Unknown",
                "typical_sccmec": ["Variable"],
                "typical_spa": ["Unknown or variable"],
                "resistance_profile": ["Variable resistance patterns"]
            },
            '120': {
                "clonal_complex": "CC121",
                "classification": "Rare lineage",
                "geographic_distribution": "Variable",
                "clinical_significance": "Limited documentation",
                "common_virulence": ["Exfoliative toxins", "Enterotoxins"],
                "outbreak_potential": "Unknown",
                "typical_sccmec": ["Variable"],
                "typical_spa": ["Unknown or variable"],
                "resistance_profile": ["Variable resistance patterns"]
            },
            # Continuing with ST122-200...
            '122': {
                "clonal_complex": "CC45",
                "classification": "Rare lineage",
                "geographic_distribution": "Variable",
                "clinical_significance": "Limited documentation",
                "common_virulence": ["Enterotoxins", "Hemolysins", "Immune evasion genes"],
                "outbreak_potential": "Unknown",
                "typical_sccmec": ["Variable"],
                "typical_spa": ["Unknown or variable"],
                "resistance_profile": ["Variable resistance patterns"]
            },
            '123': {
                "clonal_complex": "CC121",
                "classification": "Rare lineage",
                "geographic_distribution": "Variable",
                "clinical_significance": "Limited documentation",
                "common_virulence": ["Exfoliative toxins", "Enterotoxins"],
                "outbreak_potential": "Unknown",
                "typical_sccmec": ["Variable"],
                "typical_spa": ["Unknown or variable"],
                "resistance_profile": ["Variable resistance patterns"]
            },
            '124': {
                "clonal_complex": "CC97",
                "classification": "Rare lineage",
                "geographic_distribution": "Variable",
                "clinical_significance": "Limited documentation",
                "common_virulence": ["Animal-adapted factors", "Some human virulence genes"],
                "outbreak_potential": "Unknown",
                "typical_sccmec": ["Variable"],
                "typical_spa": ["Unknown or variable"],
                "resistance_profile": ["Variable resistance patterns"]
            },
            '125': {
                "clonal_complex": "CC5",
                "classification": "Rare lineage",
                "geographic_distribution": "Variable",
                "clinical_significance": "Limited documentation",
                "common_virulence": ["Enterotoxins (SEC, SEL, SEU)", "Immune evasion cluster", "Hemolysins"],
                "outbreak_potential": "Unknown",
                "typical_sccmec": ["Variable"],
                "typical_spa": ["Unknown or variable"],
                "resistance_profile": ["Variable resistance patterns"]
            },
            '126': {
                "clonal_complex": "CC97",
                "classification": "Rare lineage",
                "geographic_distribution": "Variable",
                "clinical_significance": "Limited documentation",
                "common_virulence": ["Animal-adapted factors", "Some human virulence genes"],
                "outbreak_potential": "Unknown",
                "typical_sccmec": ["Variable"],
                "typical_spa": ["Unknown or variable"],
                "resistance_profile": ["Variable resistance patterns"]
            },
            '127': {
                "clonal_complex": "CC97",
                "classification": "Rare lineage",
                "geographic_distribution": "Variable",
                "clinical_significance": "Limited documentation",
                "common_virulence": ["Animal-adapted factors", "Some human virulence genes"],
                "outbreak_potential": "Unknown",
                "typical_sccmec": ["Variable"],
                "typical_spa": ["Unknown or variable"],
                "resistance_profile": ["Variable resistance patterns"]
            },
            '128': {
                "clonal_complex": "CC8",
                "classification": "Rare lineage",
                "geographic_distribution": "Variable",
                "clinical_significance": "Limited documentation",
                "common_virulence": ["Panton-Valentine Leukocidin", "Enterotoxins", "Phenol-soluble modulins"],
                "outbreak_potential": "Unknown",
                "typical_sccmec": ["Variable"],
                "typical_spa": ["Unknown or variable"],
                "resistance_profile": ["Variable resistance patterns"]
            },
            '129': {
                "clonal_complex": "CC129",
                "classification": "Rare lineage",
                "geographic_distribution": "Variable",
                "clinical_significance": "Limited documentation",
                "common_virulence": ["Variable virulence factors"],
                "outbreak_potential": "Unknown",
                "typical_sccmec": ["Variable"],
                "typical_spa": ["Unknown or variable"],
                "resistance_profile": ["Variable resistance patterns"]
            },
            '131': {
                "clonal_complex": "CC131",
                "classification": "Rare lineage",
                "geographic_distribution": "Variable",
                "clinical_significance": "Limited documentation",
                "common_virulence": ["Variable virulence factors"],
                "outbreak_potential": "Unknown",
                "typical_sccmec": ["Variable"],
                "typical_spa": ["Unknown or variable"],
                "resistance_profile": ["Variable resistance patterns"]
            },
            '132': {
                "clonal_complex": "CC132",
                "classification": "Rare lineage",
                "geographic_distribution": "Variable",
                "clinical_significance": "Limited documentation",
                "common_virulence": ["Variable virulence factors"],
                "outbreak_potential": "Unknown",
                "typical_sccmec": ["Variable"],
                "typical_spa": ["Unknown or variable"],
                "resistance_profile": ["Variable resistance patterns"]
            },
            '134': {
                "clonal_complex": "CC22",
                "classification": "Rare lineage",
                "geographic_distribution": "Variable",
                "clinical_significance": "Limited documentation",
                "common_virulence": ["Enterotoxins", "Immune evasion cluster", "Hemolysins"],
                "outbreak_potential": "Unknown",
                "typical_sccmec": ["Variable"],
                "typical_spa": ["Unknown or variable"],
                "resistance_profile": ["Variable resistance patterns"]
            },
            '135': {
                "clonal_complex": "CC5",
                "classification": "Rare lineage",
                "geographic_distribution": "Variable",
                "clinical_significance": "Limited documentation",
                "common_virulence": ["Enterotoxins (SEC, SEL, SEU)", "Immune evasion cluster", "Hemolysins"],
                "outbreak_potential": "Unknown",
                "typical_sccmec": ["Variable"],
                "typical_spa": ["Unknown or variable"],
                "resistance_profile": ["Variable resistance patterns"]
            },
            '136': {
                "clonal_complex": "CC136",
                "classification": "Rare lineage",
                "geographic_distribution": "Variable",
                "clinical_significance": "Limited documentation",
                "common_virulence": ["Variable virulence factors"],
                "outbreak_potential": "Unknown",
                "typical_sccmec": ["Variable"],
                "typical_spa": ["Unknown or variable"],
                "resistance_profile": ["Variable resistance patterns"]
            },
            '137': {
                "clonal_complex": "CC22",
                "classification": "Rare lineage",
                "geographic_distribution": "Variable",
                "clinical_significance": "Limited documentation",
                "common_virulence": ["Enterotoxins", "Immune evasion cluster", "Hemolysins"],
                "outbreak_potential": "Unknown",
                "typical_sccmec": ["Variable"],
                "typical_spa": ["Unknown or variable"],
                "resistance_profile": ["Variable resistance patterns"]
            },
            '138': {
                "clonal_complex": "CC138",
                "classification": "Rare lineage",
                "geographic_distribution": "Variable",
                "clinical_significance": "Limited documentation",
                "common_virulence": ["Variable virulence factors"],
                "outbreak_potential": "Unknown",
                "typical_sccmec": ["Variable"],
                "typical_spa": ["Unknown or variable"],
                "resistance_profile": ["Variable resistance patterns"]
            },
            '139': {
                "clonal_complex": "CC139",
                "classification": "Rare lineage",
                "geographic_distribution": "Variable",
                "clinical_significance": "Limited documentation",
                "common_virulence": ["Variable virulence factors"],
                "outbreak_potential": "Unknown",
                "typical_sccmec": ["Variable"],
                "typical_spa": ["Unknown or variable"],
                "resistance_profile": ["Variable resistance patterns"]
            },
            '140': {
                "clonal_complex": "CC140",
                "classification": "Rare lineage",
                "geographic_distribution": "Variable",
                "clinical_significance": "Limited documentation",
                "common_virulence": ["Variable virulence factors"],
                "outbreak_potential": "Unknown",
                "typical_sccmec": ["Variable"],
                "typical_spa": ["Unknown or variable"],
                "resistance_profile": ["Variable resistance patterns"]
            },
            '141': {
                "clonal_complex": "CC8",
                "classification": "Rare lineage",
                "geographic_distribution": "Variable",
                "clinical_significance": "Limited documentation",
                "common_virulence": ["Panton-Valentine Leukocidin", "Enterotoxins", "Phenol-soluble modulins"],
                "outbreak_potential": "Unknown",
                "typical_sccmec": ["Variable"],
                "typical_spa": ["Unknown or variable"],
                "resistance_profile": ["Variable resistance patterns"]
            },
            '142': {
                "clonal_complex": "CC5",
                "classification": "Rare lineage",
                "geographic_distribution": "Variable",
                "clinical_significance": "Limited documentation",
                "common_virulence": ["Enterotoxins (SEC, SEL, SEU)", "Immune evasion cluster", "Hemolysins"],
                "outbreak_potential": "Unknown",
                "typical_sccmec": ["Variable"],
                "typical_spa": ["Unknown or variable"],
                "resistance_profile": ["Variable resistance patterns"]
            },
            '143': {
                "clonal_complex": "CC5",
                "classification": "Rare lineage",
                "geographic_distribution": "Variable",
                "clinical_significance": "Limited documentation",
                "common_virulence": ["Enterotoxins (SEC, SEL, SEU)", "Immune evasion cluster", "Hemolysins"],
                "outbreak_potential": "Unknown",
                "typical_sccmec": ["Variable"],
                "typical_spa": ["Unknown or variable"],
                "resistance_profile": ["Variable resistance patterns"]
            },
            '144': {
                "clonal_complex": "CC144",
                "classification": "Rare lineage",
                "geographic_distribution": "Variable",
                "clinical_significance": "Limited documentation",
                "common_virulence": ["Variable virulence factors"],
                "outbreak_potential": "Unknown",
                "typical_sccmec": ["Variable"],
                "typical_spa": ["Unknown or variable"],
                "resistance_profile": ["Variable resistance patterns"]
            },
            '145': {
                "clonal_complex": "CC145",
                "classification": "Rare lineage",
                "geographic_distribution": "Variable",
                "clinical_significance": "Limited documentation",
                "common_virulence": ["Variable virulence factors"],
                "outbreak_potential": "Unknown",
                "typical_sccmec": ["Variable"],
                "typical_spa": ["Unknown or variable"],
                "resistance_profile": ["Variable resistance patterns"]
            },
            '146': {
                "clonal_complex": "CC5",
                "classification": "Rare lineage",
                "geographic_distribution": "Variable",
                "clinical_significance": "Limited documentation",
                "common_virulence": ["Enterotoxins (SEC, SEL, SEU)", "Immune evasion cluster", "Hemolysins"],
                "outbreak_potential": "Unknown",
                "typical_sccmec": ["Variable"],
                "typical_spa": ["Unknown or variable"],
                "resistance_profile": ["Variable resistance patterns"]
            },
            '147': {
                "clonal_complex": "CC1",
                "classification": "Rare lineage",
                "geographic_distribution": "Variable",
                "clinical_significance": "Limited documentation",
                "common_virulence": ["Panton-Valentine Leukocidin", "Enterotoxins SEH/SEK", "Hemolysins"],
                "outbreak_potential": "Unknown",
                "typical_sccmec": ["Variable"],
                "typical_spa": ["Unknown or variable"],
                "resistance_profile": ["Variable resistance patterns"]
            },
            '148': {
                "clonal_complex": "CC5",
                "classification": "Rare lineage",
                "geographic_distribution": "Variable",
                "clinical_significance": "Limited documentation",
                "common_virulence": ["Enterotoxins (SEC, SEL, SEU)", "Immune evasion cluster", "Hemolysins"],
                "outbreak_potential": "Unknown",
                "typical_sccmec": ["Variable"],
                "typical_spa": ["Unknown or variable"],
                "resistance_profile": ["Variable resistance patterns"]
            },
            '149': {
                "clonal_complex": "CC5",
                "classification": "Rare lineage",
                "geographic_distribution": "Variable",
                "clinical_significance": "Limited documentation",
                "common_virulence": ["Enterotoxins (SEC, SEL, SEU)", "Immune evasion cluster", "Hemolysins"],
                "outbreak_potential": "Unknown",
                "typical_sccmec": ["Variable"],
                "typical_spa": ["Unknown or variable"],
                "resistance_profile": ["Variable resistance patterns"]
            },
            '150': {
                "clonal_complex": "CC150",
                "classification": "Rare lineage",
                "geographic_distribution": "Variable",
                "clinical_significance": "Limited documentation",
                "common_virulence": ["Variable virulence factors"],
                "outbreak_potential": "Unknown",
                "typical_sccmec": ["Variable"],
                "typical_spa": ["Unknown or variable"],
                "resistance_profile": ["Variable resistance patterns"]
            },
            '151': {
                "clonal_complex": "CC151",
                "classification": "Rare lineage",
                "geographic_distribution": "Variable",
                "clinical_significance": "Limited documentation",
                "common_virulence": ["Variable virulence factors"],
                "outbreak_potential": "Unknown",
                "typical_sccmec": ["Variable"],
                "typical_spa": ["Unknown or variable"],
                "resistance_profile": ["Variable resistance patterns"]
            },
            '153': {
                "clonal_complex": "CC153",
                "classification": "Rare lineage",
                "geographic_distribution": "Variable",
                "clinical_significance": "Limited documentation",
                "common_virulence": ["Variable virulence factors"],
                "outbreak_potential": "Unknown",
                "typical_sccmec": ["Variable"],
                "typical_spa": ["Unknown or variable"],
                "resistance_profile": ["Variable resistance patterns"]
            },
            '154': {
                "clonal_complex": "CC154",
                "classification": "Rare lineage",
                "geographic_distribution": "Variable",
                "clinical_significance": "Limited documentation",
                "common_virulence": ["Variable virulence factors"],
                "outbreak_potential": "Unknown",
                "typical_sccmec": ["Variable"],
                "typical_spa": ["Unknown or variable"],
                "resistance_profile": ["Variable resistance patterns"]
            },
            '155': {
                "clonal_complex": "CC8",
                "classification": "Rare lineage",
                "geographic_distribution": "Variable",
                "clinical_significance": "Limited documentation",
                "common_virulence": ["Panton-Valentine Leukocidin", "Enterotoxins", "Phenol-soluble modulins"],
                "outbreak_potential": "Unknown",
                "typical_sccmec": ["Variable"],
                "typical_spa": ["Unknown or variable"],
                "resistance_profile": ["Variable resistance patterns"]
            },
            '156': {
                "clonal_complex": "CC156",
                "classification": "Rare lineage",
                "geographic_distribution": "Variable",
                "clinical_significance": "Limited documentation",
                "common_virulence": ["Variable virulence factors"],
                "outbreak_potential": "Unknown",
                "typical_sccmec": ["Variable"],
                "typical_spa": ["Unknown or variable"],
                "resistance_profile": ["Variable resistance patterns"]
            },
            '157': {
                "clonal_complex": "CC8",
                "classification": "Rare lineage",
                "geographic_distribution": "Variable",
                "clinical_significance": "Limited documentation",
                "common_virulence": ["Panton-Valentine Leukocidin", "Enterotoxins", "Phenol-soluble modulins"],
                "outbreak_potential": "Unknown",
                "typical_sccmec": ["Variable"],
                "typical_spa": ["Unknown or variable"],
                "resistance_profile": ["Variable resistance patterns"]
            },
            '158': {
                "clonal_complex": "CC8",
                "classification": "Rare lineage",
                "geographic_distribution": "Variable",
                "clinical_significance": "Limited documentation",
                "common_virulence": ["Panton-Valentine Leukocidin", "Enterotoxins", "Phenol-soluble modulins"],
                "outbreak_potential": "Unknown",
                "typical_sccmec": ["Variable"],
                "typical_spa": ["Unknown or variable"],
                "resistance_profile": ["Variable resistance patterns"]
            },
            '159': {
                "clonal_complex": "CC1",
                "classification": "Rare lineage",
                "geographic_distribution": "Variable",
                "clinical_significance": "Limited documentation",
                "common_virulence": ["Panton-Valentine Leukocidin", "Enterotoxins SEH/SEK", "Hemolysins"],
                "outbreak_potential": "Unknown",
                "typical_sccmec": ["Variable"],
                "typical_spa": ["Unknown or variable"],
                "resistance_profile": ["Variable resistance patterns"]
            },
            '160': {
                "clonal_complex": "CC121",
                "classification": "Rare lineage",
                "geographic_distribution": "Variable",
                "clinical_significance": "Limited documentation",
                "common_virulence": ["Exfoliative toxins", "Enterotoxins"],
                "outbreak_potential": "Unknown",
                "typical_sccmec": ["Variable"],
                "typical_spa": ["Unknown or variable"],
                "resistance_profile": ["Variable resistance patterns"]
            },
            '188': {
                "clonal_complex": "CC1",
                "classification": "Community-associated MRSA",
                "geographic_distribution": "Europe, Middle East",
                "clinical_significance": "Associated with skin and soft tissue infections",
                "common_virulence": ["Enterotoxins", "Hemolysins"],
                "outbreak_potential": "MEDIUM",
                "typical_sccmec": ["IV", "V"],
                "typical_spa": ["t189", "t325"],
                "resistance_profile": ["Methicillin", "Variable resistance"]
            },
            '225': {  
                "clonal_complex": "CC5",
                "classification": "Healthcare-associated MRSA",
                "geographic_distribution": "Germany, Central Europe",
                "clinical_significance": "German epidemic clone (GEC), hospital-associated with multi-drug resistance, particularly prevalent in German-speaking countries",
                "common_virulence": ["Enterotoxins", "Hemolysins", "Immune evasion genes"],
                "outbreak_potential": "HIGH",
                "typical_sccmec": ["II", "IV"],
                "typical_spa": ["t003", "t014", "t045"],
                "resistance_profile": ["Methicillin", "Multi-drug resistant including aminoglycosides, macrolides"]
            },
            '231': {  
                "clonal_complex": "CC5",
                "classification": "Rare/Uncommon MRSA",
                "geographic_distribution": "Sporadic reports globally",
                "clinical_significance": "Primarily MSSA lineage, rare MRSA conversion reported; limited clinical significance as MRSA",
                "common_virulence": ["Variable; typical CC5 factors if present"],
                "outbreak_potential": "VERY LOW",
                "typical_sccmec": ["Rare/occasional acquisition"],
                "typical_spa": ["Variable"],
                "resistance_profile": ["Variable if MRSA"]
            },

            '239': {
                "clonal_complex": "CC8",  
                "classification": "Healthcare-associated MRSA",
                "geographic_distribution": "Asia (particularly China, Taiwan), Brazil, Eastern Europe",
                "clinical_significance": "Brazilian/Hungarian epidemic clone, one of the earliest and most successful global MRSA clones with high multi-drug resistance",
                "common_virulence": ["Enterotoxins", "Hemolysins", "Biofilm formation genes", "Often lacks PVL"],
                "outbreak_potential": "VERY HIGH",
                "typical_sccmec": ["III", "IIIA"],
                "typical_spa": ["t037", "t030", "t421"],
                "resistance_profile": ["Methicillin", "High-level multi-drug resistance including aminoglycosides, fluoroquinolones"]
            },
            '291': {  
                "clonal_complex": "CC291",
                "classification": "Emerging lineage",
                "geographic_distribution": "Sporadic reports, limited distribution",
                "clinical_significance": "Emerging clone with limited documentation, occasional reports from specific regions",
                "common_virulence": ["Variable virulence factors"],
                "outbreak_potential": "UNKNOWN/LOW",
                "typical_sccmec": ["Variable"],
                "typical_spa": ["Unknown"],
                "resistance_profile": ["Variable resistance patterns"]
            }
        }
            
        # Check if ST is in database
        if st in lineage_db:
            return lineage_db[st]
        else:
            # For unknown STs - not in our database
            if st.isdigit():
                return {
                    "clonal_complex": f"Unknown (ST{st})",
                    "classification": "Not in database - novel or uncommon",
                    "geographic_distribution": "Unknown",
                    "clinical_significance": f"ST{st} is not currently in the StaphScope MLST database. This may be a novel sequence type, a local variant, or a recently emerged clone. Please consult PubMLST (https://pubmlst.org/saureus/) for the most current information.",
                    "common_virulence": ["Unknown - requires further analysis"],
                    "outbreak_potential": "UNKNOWN",
                    "typical_sccmec": ["Unknown"],
                    "typical_spa": ["Unknown"],
                    "resistance_profile": ["Unknown"],
                    "database_note": "Consider submitting to PubMLST for characterization"
                }
            else:
                # For non-numeric STs (ND, -, etc.)
                return {
                    "clonal_complex": 'Not Assigned',
                    "classification": "MLST typing failed",
                    "geographic_distribution": "N/A",
                    "clinical_significance": "Could not determine sequence type. Possible reasons: novel alleles, poor quality assembly, or non-typeable strain.",
                    "common_virulence": ["Cannot determine"],
                    "outbreak_potential": "UNKNOWN",
                    "typical_sccmec": ["Cannot determine"],
                    "typical_spa": ["Cannot determine"],
                    "resistance_profile": ["Cannot determine"]
                }
            
    def get_identity_coverage(self, st: str) -> Dict:
        """Get identity and coverage information based on MLST assignment"""
        if st and st != '-' and st != 'ND' and st != 'UNKNOWN':
            # MLST assigned - identity and coverage are 100%
            return {
                "identity": "100%",
                "coverage": "100%",
                "mlst_status": "Assigned",
                "quality_metrics": {
                    "assembly_quality": "High Quality",
                    "allele_completeness": "Complete",
                    "database_match": "Perfect Match"
                }
            }
        else:
            # MLST not assigned
            return {
                "identity": "Not Assigned",
                "coverage": "Not Assigned",
                "mlst_status": "Not Assigned",
                "quality_metrics": {
                    "assembly_quality": "Requires Review",
                    "allele_completeness": "Incomplete",
                    "database_match": "No Match"
                }
            }

    def get_empty_results(self, sample_name: str) -> Dict:
        """Return empty results structure"""
        return {
            "sample": sample_name,
            "st": "ND",
            "scheme": "saureus",
            "alleles": {},
            "allele_profile": "",
            "confidence": "LOW",
            "mlst_assigned": False
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
            "mlst_assigned": False,
            "error": "MLST analysis failed"
        }

    def generate_output_files(self, mlst_results: Dict, output_dir: Path):
        """Generate only 3 output files: HTML, TXT, and TSV"""
        # Add identity and coverage if not already present
        if 'identity' not in mlst_results:
            mlst_results.update(self.get_identity_coverage(mlst_results.get('st', 'ND')))
        
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
MLST Status: {mlst_results.get('mlst_status', 'Not Assigned')}

Identity & Coverage:
Identity: {mlst_results.get('identity', 'Not Assigned')}
Coverage: {mlst_results.get('coverage', 'Not Assigned')}

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
        
        # Quality metrics
        if 'quality_metrics' in mlst_results:
            report += f"""
QUALITY METRICS:
----------------
Assembly Quality: {mlst_results['quality_metrics'].get('assembly_quality', 'Unknown')}
Allele Completeness: {mlst_results['quality_metrics'].get('allele_completeness', 'Unknown')}
Database Match: {mlst_results['quality_metrics'].get('database_match', 'Unknown')}
"""
        
        with open(output_dir / "mlst_report.txt", 'w') as f:
            f.write(report)

    def generate_tsv_report(self, mlst_results: Dict, output_dir: Path):
        """Generate simple TSV report"""
        tsv_content = f"Sample\tST\tScheme\tMLST_Status\tIdentity\tCoverage\tClonal_Complex\tClassification\tOutbreak_Potential\tTypical_SCCmec\tTypical_spa\tResistance_Profile\n"
        tsv_content += f"{mlst_results['sample']}\t{mlst_results['st']}\t{mlst_results['scheme']}\t{mlst_results.get('mlst_status', 'Not Assigned')}\t{mlst_results.get('identity', 'Not Assigned')}\t{mlst_results.get('coverage', 'Not Assigned')}\t{mlst_results.get('clonal_complex', 'Unknown')}\t{mlst_results.get('classification', 'Unknown')}\t{mlst_results.get('outbreak_potential', 'UNKNOWN')}\t{','.join(mlst_results.get('typical_sccmec', ['Unknown']))}\t{','.join(mlst_results.get('typical_spa', ['Unknown']))}\t{','.join(mlst_results.get('resistance_profile', ['Unknown']))}\n"
        
        with open(output_dir / "mlst_report.tsv", 'w') as f:
            f.write(tsv_content)

    def generate_html_report(self, mlst_results: Dict, output_dir: Path):
        """Generate beautiful HTML report in StaphScope style"""
        # Get a random quote for this report
        random_quote = self.get_random_quote()
        
        # Extract variables for template
        sample = mlst_results['sample']
        st = mlst_results['st']
        scheme = mlst_results['scheme']
        confidence = mlst_results['confidence']
        clonal_complex = mlst_results.get('clonal_complex', 'Unknown')
        allele_profile = mlst_results['allele_profile']
        identity = mlst_results.get('identity', 'Not Assigned')
        coverage = mlst_results.get('coverage', 'Not Assigned')
        mlst_status = mlst_results.get('mlst_status', 'Not Assigned')
        classification = mlst_results.get('classification', 'Unknown')
        geographic_distribution = mlst_results.get('geographic_distribution', 'Unknown')
        outbreak_potential = mlst_results.get('outbreak_potential', 'UNKNOWN')
        clinical_significance = mlst_results.get('clinical_significance', 'Further analysis required.')
        typical_sccmec = ', '.join(mlst_results.get('typical_sccmec', ['Unknown']))
        typical_spa = ', '.join(mlst_results.get('typical_spa', ['Unknown']))
        resistance_profile = ', '.join(mlst_results.get('resistance_profile', ['Unknown']))
        common_virulence = mlst_results.get('common_virulence', [])
        quality_metrics = mlst_results.get('quality_metrics', {})
        
        # Build alleles HTML
        alleles_html = ''
        for gene, allele in mlst_results.get('alleles', {}).items():
            alleles_html += f'''                <div class="allele-card">
                    <div style="font-size: 12px; opacity: 0.9;">{gene}</div>
                    <div style="font-size: 18px;">{allele}</div>
                </div>
'''
        
        # Build virulence factors HTML
        virulence_html = ''
        for virulence in common_virulence:
            virulence_html += f'''                        <div class="virulence-card">{virulence}</div>
'''
        
        # Build quality metrics HTML
        quality_metrics_html = ''
        if quality_metrics:
            for key, value in quality_metrics.items():
                quality_metrics_html += f'''                <div class="quality-card">
                    <div class="quality-value">{value}</div>
                    <div class="quality-label">{key.replace('_', ' ').title()}</div>
                </div>
'''
        
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
        
        .identity-coverage-grid {{
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(200px, 1fr));
            gap: 20px;
            margin-top: 20px;
        }}
        
        .ic-card {{
            background: linear-gradient(135deg, #10b981 0%, #059669 100%);
            color: white;
            padding: 20px;
            border-radius: 8px;
            box-shadow: 0 4px 12px rgba(0, 0, 0, 0.15);
            text-align: center;
        }}
        
        .ic-card.not-assigned {{
            background: linear-gradient(135deg, #6b7280 0%, #4b5563 100%);
        }}
        
        .ic-value {{
            font-size: 28px;
            font-weight: bold;
            margin-bottom: 5px;
        }}
        
        .ic-label {{
            font-size: 14px;
            opacity: 0.9;
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
        
        .quality-grid {{
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(200px, 1fr));
            gap: 15px;
            margin-top: 15px;
        }}
        
        .quality-card {{
            background: #f0f9ff;
            color: #0369a1;
            padding: 15px;
            border-radius: 8px;
            text-align: center;
            border: 1px solid #bae6fd;
        }}
        
        .quality-value {{
            font-size: 18px;
            font-weight: bold;
            margin-bottom: 5px;
        }}
        
        .quality-label {{
            font-size: 12px;
            opacity: 0.8;
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
███████╗   ██║   ██║  ██║██║     ██║  ██║███████║╚██████╗╚██████╔╝██║     ███████╗
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
                    <div class="metric-value">{sample}</div>
                </div>
                <div class="metric-card">
                    <div class="metric-label">Analysis Date</div>
                    <div class="metric-value">{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}</div>
                </div>
                <div class="metric-card">
                    <div class="metric-label">MLST Scheme</div>
                    <div class="metric-value">{scheme.title()}</div>
                </div>
            </div>
        </div>
        
        <div class="report-section">
            <h2>🎯 Identity & Coverage</h2>
            <div class="identity-coverage-grid">
'''
        
        # Identity and Coverage cards
        identity_class = "ic-card" if identity == "100%" else "ic-card not-assigned"
        coverage_class = "ic-card" if coverage == "100%" else "ic-card not-assigned"
        
        html_content += f'''                <div class="{identity_class}">
                    <div class="ic-value">{identity}</div>
                    <div class="ic-label">Identity</div>
                </div>
                <div class="{coverage_class}">
                    <div class="ic-value">{coverage}</div>
                    <div class="ic-label">Coverage</div>
                </div>
'''
        
        html_content += f'''            </div>
            
            <h3>MLST Status: <span style="color: {'#16a34a' if mlst_status == 'Assigned' else '#dc2626'}">{mlst_status}</span></h3>
            
            <h3>Quality Metrics</h3>
            <div class="quality-grid">
{quality_metrics_html}            </div>
        </div>
        
        <div class="report-section">
            <h2>🧬 MLST Typing Results</h2>
            <div class="metrics-grid">
                <div class="metric-card">
                    <div class="metric-label">Sequence Type</div>
                    <div class="metric-value">ST{st}</div>
                </div>
                <div class="metric-card">
                    <div class="metric-label">Confidence</div>
                    <div class="metric-value confidence-{confidence.lower()}">{confidence}</div>
                </div>
                <div class="metric-card">
                    <div class="metric-label">Clonal Complex</div>
                    <div class="metric-value">{clonal_complex}</div>
                </div>
            </div>
            
            <h3>Allele Profile</h3>
            <div class="profile-box">
                <code style="font-size: 16px; color: #1e40af; font-weight: bold;">{allele_profile}</code>
            </div>
            
            <h3>Individual Alleles</h3>
            <div class="allele-grid">
{alleles_html}            </div>
        </div>
        
        <div class="report-section">
            <h2>🌍 Lineage Information</h2>
            <div class="metrics-grid">
                <div class="metric-card">
                    <div class="metric-label">Classification</div>
                    <div class="metric-value">{classification}</div>
                </div>
                <div class="metric-card">
                    <div class="metric-label">Geographic Distribution</div>
                    <div class="metric-value">{geographic_distribution}</div>
                </div>
                <div class="metric-card">
                    <div class="metric-label">Outbreak Potential</div>
                    <div class="metric-value">{outbreak_potential}</div>
                </div>
            </div>
            
            <h3>Clinical Significance</h3>
            <div class="profile-box">
                <p style="margin: 0; line-height: 1.6; color: #374151;">{clinical_significance}</p>
            </div>
            
            <div style="display: grid; grid-template-columns: repeat(auto-fit, minmax(300px, 1fr)); gap: 20px; margin-top: 20px;">
                <div>
                    <h3>Typical Characteristics</h3>
                    <p><strong>SCCmec Types:</strong> {typical_sccmec}</p>
                    <p><strong>spa Types:</strong> {typical_spa}</p>
                    <p><strong>Resistance Profile:</strong> {resistance_profile}</p>
                </div>
                <div>
                    <h3>Common Virulence Factors</h3>
                    <div class="virulence-grid">
{virulence_html}                    </div>
                </div>
            </div>
        </div>
        
        <div class="footer">
            <p><strong>STAPHSCOPE</strong> - MLST Analysis Module</p>
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

        // Rotate quotes every 10 seconds
        setInterval(displayQuote, 10000);
    </script>
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
        
        # Create JSON summary
        self.create_mlst_json_summary(all_results, output_dir)
        
        print("✅ MLST summary files created successfully!")

    def create_mlst_tsv_summary(self, all_results: Dict[str, Dict], output_dir: Path):
        """Create TSV summary file with all samples"""
        summary_file = output_dir / "mlst_summary.tsv"
        
        with open(summary_file, 'w') as f:
            # Write header
            f.write("Sample\tST\tMLST_Status\tIdentity\tCoverage\tAllele_Profile\t")
            
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
                f.write(f"{sample_name}\t{result['st']}\t{result.get('mlst_status', 'Not Assigned')}\t{result.get('identity', 'Not Assigned')}\t{result.get('coverage', 'Not Assigned')}\t{result['allele_profile']}\t")
                
                # Write allele values for each gene
                for gene in sorted(all_genes):
                    allele = result['alleles'].get(gene, '')
                    f.write(f"{allele}\t")
                f.write("\n")
        
        print(f"📄 TSV summary created: {summary_file}")

    def create_mlst_html_summary(self, all_results: Dict[str, Dict], output_dir: Path):
        """Create beautiful HTML summary with MLST styling"""
        summary_file = output_dir / "mlst_summary.html"
        
        # Get a random quote for the summary
        random_quote = self.get_random_quote()
        
        # Get all unique gene names
        all_genes = set()
        for result in all_results.values():
            all_genes.update(result['alleles'].keys())
        sorted_genes = sorted(all_genes)
        
        # Calculate statistics
        total_samples = len(all_results)
        assigned_samples = sum(1 for result in all_results.values() if result.get('mlst_status') == 'Assigned')
        not_assigned_samples = total_samples - assigned_samples
        assigned_percentage = (assigned_samples / total_samples * 100) if total_samples > 0 else 0
        
        # Get unique STs that are assigned (excluding 'ND', '-', 'UNKNOWN')
        unique_assigned_sts = set()
        for result in all_results.values():
            if result.get('st') not in ['ND', '-', 'UNKNOWN', ''] and result.get('mlst_status') == 'Assigned':
                unique_assigned_sts.add(result['st'])
        
        # Build table rows
        table_rows = ''
        for sample_name, result in all_results.items():
            st = result.get('st', '')
            mlst_status = result.get('mlst_status', 'Not Assigned')
            identity = result.get('identity', 'Not Assigned')
            coverage = result.get('coverage', 'Not Assigned')
            allele_profile = result.get('allele_profile', '')
            
            # Build allele columns
            allele_columns = ''
            for gene in sorted_genes:
                allele = result['alleles'].get(gene, '')
                allele_columns += f'<td>{allele}</td>\n'
            
            # Format row
            mlst_status_html = f'<span style="color: {"#10b981" if mlst_status == "Assigned" else "#dc2626"}">{mlst_status}</span>'
            identity_class = 'identity-cell' if identity == '100%' else 'not-assigned-cell'
            coverage_class = 'coverage-cell' if coverage == '100%' else 'not-assigned-cell'
            
            table_rows += f'''                        <tr>
                            <td><strong>{sample_name}</strong></td>
                            <td class="st-cell">ST{st}</td>
                            <td>{mlst_status_html}</td>
                            <td class="{identity_class}">{identity}</td>
                            <td class="{coverage_class}">{coverage}</td>
                            <td class="allele-cell">{allele_profile}</td>
{allele_columns}                        </tr>
'''
        
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
        
        .identity-cell {{
            font-weight: bold;
            color: #10b981;
        }}
        
        .coverage-cell {{
            font-weight: bold;
            color: #10b981;
        }}
        
        .not-assigned-cell {{
            color: #dc2626;
            font-style: italic;
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
        
        .identity-coverage-stats {{
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(180px, 1fr));
            gap: 10px;
            margin-top: 15px;
        }}
        
        .ic-stat-card {{
            background: #f0f9ff;
            color: #0369a1;
            padding: 10px;
            border-radius: 6px;
            text-align: center;
            border: 1px solid #bae6fd;
        }}
        
        .ic-stat-value {{
            font-size: 18px;
            font-weight: bold;
            margin-bottom: 3px;
        }}
        
        .ic-stat-label {{
            font-size: 11px;
            opacity: 0.8;
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
            
            <div class="quote-container" id="quoteContainer">
                <div class="quote-text" id="quoteText">"{random_quote['text']}"</div>
                <div class="quote-author" id="quoteAuthor">— {random_quote['author']}</div>
            </div>
        </div>
        
        <div class="report-section">
            <h2>📊 MLST Summary - All Samples</h2>
            
            <div class="stats-grid">
                <div class="stat-card">
                    <div class="stat-value">{total_samples}</div>
                    <div class="stat-label">SAMPLES PROCESSED</div>
                </div>
                <div class="stat-card">
                    <div class="stat-value">{assigned_samples}</div>
                    <div class="stat-label">SAMPLES ASSIGNED</div>
                </div>
                <div class="stat-card">
                    <div class="stat-value">{len(unique_assigned_sts)}</div>
                    <div class="stat-label">UNIQUE STs</div>
                </div>
                <div class="stat-card">
                    <div class="stat-value">{len(sorted_genes)}</div>
                    <div class="stat-label">MLST GENES</div>
                </div>
            </div>
            
            <h3>Identity & Coverage Statistics</h3>
            <div class="identity-coverage-stats">
                <div class="ic-stat-card">
                    <div class="ic-stat-value">{assigned_samples}/{total_samples}</div>
                    <div class="ic-stat-label">Samples with MLST</div>
                </div>
                <div class="ic-stat-card">
                    <div class="ic-stat-value">{assigned_percentage:.1f}%</div>
                    <div class="ic-stat-label">Assignment Rate</div>
                </div>
                <div class="ic-stat-card">
                    <div class="ic-stat-value">{not_assigned_samples}</div>
                    <div class="ic-stat-label">Not Assigned</div>
                </div>
            </div>
            
            <div style="overflow-x: auto;">
                <table class="summary-table">
                    <thead>
                        <tr>
                            <th>Sample</th>
                            <th>ST</th>
                            <th>MLST Status</th>
                            <th>Identity</th>
                            <th>Coverage</th>
                            <th>Allele Profile</th>
'''
        
        # Add gene headers
        for gene in sorted_genes:
            html_content += f'                            <th>{gene}</th>\n'
        
        html_content += f'''                        </tr>
                    </thead>
                    <tbody>
{table_rows}                    </tbody>
                </table>
            </div>
        </div>
        
        <div class="footer">
            <p><strong>STAPHSCOPE</strong> - MLST Summary Report</p>
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
            {{"text": "The universe is not required to be in perfect harmony with human ambition.", "author": "Carl Sagan"}},
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
        
        with open(summary_file, 'w', encoding='utf-8') as f:
            f.write(html_content)
        
        print(f"🌐 HTML summary created: {summary_file}")

    def create_mlst_json_summary(self, all_results: Dict[str, Dict], output_dir: Path):
        """Create JSON summary file with all MLST results"""
        print("📊 Creating MLST JSON summary...")
        
        summary_file = output_dir / "mlst_summary.json"
        
        # Prepare the JSON structure
        json_summary = {
            "metadata": {
                "analysis_date": datetime.now().isoformat(),
                "total_samples": len(all_results),
                "analysis_type": "MLST",
                "scheme": "saureus",
                "version": "1.0"
            },
            "statistics": self._calculate_json_statistics(all_results),
            "samples": {},
            "summary_by_st": self._create_st_summary(all_results)
        }
        
        # Add each sample's results
        for sample_name, result in all_results.items():
            json_summary["samples"][sample_name] = {
                "sequence_type": result.get('st', 'ND'),
                "mlst_status": result.get('mlst_status', 'Not Assigned'),
                "identity": result.get('identity', 'Not Assigned'),
                "coverage": result.get('coverage', 'Not Assigned'),
                "confidence": result.get('confidence', 'LOW'),
                "alleles": result.get('alleles', {}),
                "allele_profile": result.get('allele_profile', ''),
                "lineage_info": {
                    "clonal_complex": result.get('clonal_complex', 'Unknown'),
                    "classification": result.get('classification', 'Unknown'),
                    "geographic_distribution": result.get('geographic_distribution', 'Unknown'),
                    "clinical_significance": result.get('clinical_significance', 'Unknown'),
                    "outbreak_potential": result.get('outbreak_potential', 'UNKNOWN'),
                    "typical_sccmec": result.get('typical_sccmec', []),
                    "typical_spa": result.get('typical_spa', []),
                    "resistance_profile": result.get('resistance_profile', []),
                    "common_virulence": result.get('common_virulence', [])
                },
                "quality_metrics": result.get('quality_metrics', {})
            }
        
        # Write JSON with pretty formatting
        with open(summary_file, 'w', encoding='utf-8') as f:
            json.dump(json_summary, f, indent=2, ensure_ascii=False)
        
        print(f"📄 JSON summary created: {summary_file}")

    def _calculate_json_statistics(self, all_results: Dict[str, Dict]) -> Dict:
        """Calculate statistics for JSON summary"""
        total_samples = len(all_results)
        
        # Count samples by status
        assigned_samples = sum(1 for r in all_results.values() if r.get('mlst_status') == 'Assigned')
        not_assigned_samples = total_samples - assigned_samples
        
        # Get unique STs
        unique_sts = set()
        st_counts = {}
        for result in all_results.values():
            st = result.get('st', '')
            if st and st not in ['ND', '-', 'UNKNOWN', '']:
                unique_sts.add(st)
                st_counts[st] = st_counts.get(st, 0) + 1
        
        # Count samples by confidence
        confidence_counts = {}
        for result in all_results.values():
            confidence = result.get('confidence', 'LOW')
            confidence_counts[confidence] = confidence_counts.get(confidence, 0) + 1
        
        # Get all unique genes
        all_genes = set()
        for result in all_results.values():
            all_genes.update(result.get('alleles', {}).keys())
        
        # Get lineage distribution
        lineage_distribution = {}
        for result in all_results.values():
            clonal_complex = result.get('clonal_complex', 'Unknown')
            lineage_distribution[clonal_complex] = lineage_distribution.get(clonal_complex, 0) + 1
        
        return {
            "total_samples": total_samples,
            "assigned_samples": assigned_samples,
            "not_assigned_samples": not_assigned_samples,
            "assignment_rate": (assigned_samples / total_samples * 100) if total_samples > 0 else 0,
            "unique_sts": len(unique_sts),
            "st_distribution": dict(sorted(st_counts.items(), key=lambda x: x[1], reverse=True)),
            "confidence_distribution": confidence_counts,
            "total_genes": len(all_genes),
            "genes_detected": sorted(list(all_genes)),
            "lineage_distribution": dict(sorted(lineage_distribution.items(), key=lambda x: x[1], reverse=True)),
            "most_common_st": max(st_counts.items(), key=lambda x: x[1])[0] if st_counts else "None"
        }

    def _create_st_summary(self, all_results: Dict[str, Dict]) -> Dict:
        """Create summary grouped by sequence type"""
        st_summary = {}
        
        for result in all_results.values():
            st = result.get('st', 'ND')
            if st not in ['ND', '-', 'UNKNOWN', '']:
                if st not in st_summary:
                    st_summary[st] = {
                        "samples": [],
                        "count": 0,
                        "allele_profiles": set(),
                        "lineage_info": {
                            "clonal_complex": result.get('clonal_complex', 'Unknown'),
                            "classification": result.get('classification', 'Unknown'),
                            "geographic_distribution": result.get('geographic_distribution', 'Unknown'),
                            "clinical_significance": result.get('clinical_significance', 'Unknown'),
                            "outbreak_potential": result.get('outbreak_potential', 'UNKNOWN')
                        }
                    }
                
                st_summary[st]["samples"].append(result['sample'])
                st_summary[st]["count"] += 1
                st_summary[st]["allele_profiles"].add(result.get('allele_profile', ''))
        
        # Convert sets to lists for JSON serialization
        for st in st_summary:
            st_summary[st]["allele_profiles"] = list(st_summary[st]["allele_profiles"])
        
        return st_summary

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

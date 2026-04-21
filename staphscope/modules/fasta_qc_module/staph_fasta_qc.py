#!/usr/bin/env python3
"""
StaphScope FASTA QC - Comprehensive Quality Control with Beautiful HTML Reports
Author: Brown Beckley <brownbeckley94@gmail.com>
Affiliation: University of Ghana Medical School - Department of Medical Biochemistry
Date: 2026-01-28
Send a quick mail for any issues or further explanations.
"""

import os
import sys
import glob
import json
import math
import statistics
from pathlib import Path
from datetime import datetime
from collections import Counter, defaultdict
import argparse
import logging
import subprocess
from typing import List, Dict, Any, Optional
from concurrent.futures import ThreadPoolExecutor, as_completed

# BioPython imports
from Bio import SeqIO
from Bio.SeqUtils import gc_fraction

class StaphFASTAQC:
    """Comprehensive FASTA QC for Staphylococcus aureus with beautiful HTML"""
    
    def __init__(self, cpus: int = None):
        # Setup logging
        self.logger = self._setup_logging()
        
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
        
        # Metadata
        self.metadata = {
            "tool_name": "StaphScope FASTA QC Analysis",
            "version": "1.2.1", 
            "authors": ["Brown Beckley"],
            "email": "brownbeckley94@gmail.com",
            "github": "https://github.com/bbeckley-hub",
            "affiliation": "University of Ghana Medical School - Department of Medical Biochemistry",
            "analysis_date": datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            "biopython_version": self._get_biopython_version()
        }
        
        # S. aureus specific thresholds
        self.thresholds = {
            'gc_normal': (32, 38),  # S. aureus specific GC range
            'short_seq': 100,
            'long_seq': 1000000,
            'max_n_run': 100,
            'max_homopolymer': 20,
            'ambiguous_critical': 5.0,
            'ambiguous_warning': 1.0,
        }
        
        self.cpus = cpus or os.cpu_count() or 1
    
    def _setup_logging(self):
        """Setup logging"""
        logging.basicConfig(
            level=logging.INFO,
            format='%(asctime)s - %(levelname)s - %(message)s'
        )
        return logging.getLogger(__name__)
    
    def _get_biopython_version(self) -> str:
        """Get BioPython version"""
        try:
            import Bio
            return Bio.__version__
        except:
            return "Unknown"
    
    def get_random_quote(self):
        """Get a random science quote"""
        import random
        return random.choice(self.science_quotes)
    
    def analyze_file(self, fasta_file: str) -> Dict[str, Any]:
        """Comprehensive analysis of a FASTA file for S. aureus"""
        try:
            self.logger.info(f"🔬 Analyzing {os.path.basename(fasta_file)}...")
            
            # Read sequences
            sequences = list(SeqIO.parse(fasta_file, "fasta"))
            if not sequences:
                return {
                    'filename': os.path.basename(fasta_file),
                    'filepath': fasta_file,
                    'status': 'error',
                    'error': 'No sequences found'
                }
            
            # Basic stats
            seq_lengths = [len(seq) for seq in sequences]
            total_length = sum(seq_lengths)
            
            # Sort for N statistics
            sorted_lengths = sorted(seq_lengths, reverse=True)
            
            # N statistics
            n50 = self._calculate_nx(sorted_lengths, total_length, 50)
            n75 = self._calculate_nx(sorted_lengths, total_length, 75)
            n90 = self._calculate_nx(sorted_lengths, total_length, 90)
            
            # L statistics
            l50 = self._calculate_lx(sorted_lengths, total_length, 50)
            l75 = self._calculate_lx(sorted_lengths, total_length, 75)
            l90 = self._calculate_lx(sorted_lengths, total_length, 90)
            
            # Nucleotide composition
            base_counts = Counter()
            gc_contents = []
            
            for seq in sequences:
                seq_str = str(seq.seq).upper()
                base_counts.update(seq_str)
                gc_contents.append(gc_fraction(seq_str) * 100)
            
            total_bases = sum(base_counts.values())
            
            # Calculate percentages
            a_count = base_counts.get('A', 0)
            t_count = base_counts.get('T', 0)
            g_count = base_counts.get('G', 0)
            c_count = base_counts.get('C', 0)
            
            a_percent = (a_count / total_bases) * 100 if total_bases > 0 else 0
            t_percent = (t_count / total_bases) * 100 if total_bases > 0 else 0
            g_percent = (g_count / total_bases) * 100 if total_bases > 0 else 0
            c_percent = (c_count / total_bases) * 100 if total_bases > 0 else 0
            gc_percent = g_percent + c_percent
            at_percent = a_percent + t_percent
            
            # Ambiguous bases
            ambiguous_bases = sum(base_counts.get(b, 0) for b in ['N', 'Y', 'R', 'W', 'S', 'K', 'M', 'B', 'D', 'H', 'V'])
            ambiguous_percent = (ambiguous_bases / total_bases) * 100 if total_bases > 0 else 0
            
            # N statistics
            sequences_with_n = sum(1 for seq in sequences if 'N' in str(seq.seq).upper())
            total_n_bases = base_counts.get('N', 0)
            
            # Find longest N-run
            max_n_run = 0
            n_runs = []
            for seq in sequences:
                seq_str = str(seq.seq).upper()
                import re
                for match in re.finditer(r'N+', seq_str):
                    run_len = len(match.group())
                    n_runs.append(run_len)
                    if run_len > max_n_run:
                        max_n_run = run_len
            
            # Homopolymers
            homopolymers = []
            for seq in sequences:
                seq_str = str(seq.seq).upper()
                for base in ['A', 'T', 'G', 'C']:
                    for match in re.finditer(f'{base}+', seq_str):
                        if len(match.group()) > 4:
                            homopolymers.append({
                                'base': base,
                                'length': len(match.group()),
                                'position': match.start()
                            })
            
            max_homopolymer = max([h['length'] for h in homopolymers]) if homopolymers else 0
            
            # Duplicate sequences
            seq_hashes = set()
            duplicate_sequences = 0
            for seq in sequences:
                seq_hash = hash(str(seq.seq))
                if seq_hash in seq_hashes:
                    duplicate_sequences += 1
                else:
                    seq_hashes.add(seq_hash)
            
            # Short and long sequences
            short_sequences = sum(1 for length in seq_lengths if length < self.thresholds['short_seq'])
            long_sequences = sum(1 for length in seq_lengths if length > self.thresholds['long_seq'])
            
            # Length distribution
            length_distribution = self._create_length_bins(seq_lengths)
            
            # S. aureus specific checks
            s_aureus_status = self._check_s_aureus_specific(gc_percent, gc_contents, total_length, seq_lengths)
            
            # Compile results
            results = {
                'filename': os.path.basename(fasta_file),
                'filepath': fasta_file,
                'file_size_mb': os.path.getsize(fasta_file) / (1024 ** 2),
                'analysis_date': datetime.now().strftime('%Y-%m-%d %H:%M:%S'),
                'status': 'success',
                'total_sequences': len(sequences),
                'total_length': total_length,
                'total_bases': total_bases,
                'longest_sequence': max(seq_lengths),
                'shortest_sequence': min(seq_lengths),
                'mean_length': statistics.mean(seq_lengths) if seq_lengths else 0,
                'median_length': statistics.median(seq_lengths) if seq_lengths else 0,
                'n50': n50,
                'n75': n75,
                'n90': n90,
                'l50': l50,
                'l75': l75,
                'l90': l90,
                'gc_percent': gc_percent,
                'at_percent': at_percent,
                'a_percent': a_percent,
                't_percent': t_percent,
                'g_percent': g_percent,
                'c_percent': c_percent,
                'ambiguous_percent': ambiguous_percent,
                'sequences_with_n': sequences_with_n,
                'total_n_bases': total_n_bases,
                'max_n_run': max_n_run,
                'homopolymers_count': len(homopolymers),
                'max_homopolymer': max_homopolymer,
                'duplicate_sequences': duplicate_sequences,
                'short_sequences': short_sequences,
                'long_sequences': long_sequences,
                'gc_distribution': {
                    'mean': statistics.mean(gc_contents) if gc_contents else 0,
                    'median': statistics.median(gc_contents) if gc_contents else 0,
                    'min': min(gc_contents) if gc_contents else 0,
                    'max': max(gc_contents) if gc_contents else 0,
                },
                'length_distribution': length_distribution,
                'base_counts': dict(base_counts),
                'n_runs': n_runs,
                's_aureus_status': s_aureus_status,
                'warnings': self._generate_warnings(
                    gc_percent, ambiguous_percent, max_n_run, max_homopolymer,
                    short_sequences, long_sequences, duplicate_sequences, len(sequences),
                    s_aureus_status
                )
            }
            
            self.logger.info(f"✅ {os.path.basename(fasta_file)}: {len(sequences)} sequences, {total_length:,} bp, N50: {n50:,}, GC: {gc_percent:.1f}%")
            return results
            
        except Exception as e:
            self.logger.error(f"❌ Error analyzing {fasta_file}: {e}")
            return {
                'filename': os.path.basename(fasta_file),
                'status': 'error',
                'error': str(e)
            }
    
    def _check_s_aureus_specific(self, gc_percent: float, gc_contents: List[float], 
                                total_length: int, seq_lengths: List[int]) -> Dict:
        """Check S. aureus specific characteristics"""
        # S. aureus genome typically 2.7-3.0 Mbp
        genome_size_ok = total_length >= 2500000 and total_length <= 3200000
        
        # GC content should be 32-38%
        gc_ok = self.thresholds['gc_normal'][0] <= gc_percent <= self.thresholds['gc_normal'][1]
        
        # Check for typical contig count (assembly quality)
        if len(seq_lengths) <= 100:
            assembly_quality = "Good" if len(seq_lengths) <= 50 else "Moderate"
        else:
            assembly_quality = "Fragmented"
        
        return {
            'genome_size_mbp': total_length / 1000000,
            'genome_size_status': 'Normal' if genome_size_ok else 'Atypical',
            'gc_status': 'Normal' if gc_ok else 'Atypical',
            'assembly_quality': assembly_quality,
            'contig_count': len(seq_lengths),
            'expected_gc_range': f"{self.thresholds['gc_normal'][0]}-{self.thresholds['gc_normal'][1]}%",
            'expected_genome_size': "2.7-3.2 Mbp"
        }
    
    def _calculate_nx(self, sorted_lengths: List[int], total_length: int, x: int) -> int:
        """Calculate Nx (N50, N75, N90)"""
        if not sorted_lengths:
            return 0
        
        target = total_length * (x / 100)
        cumulative = 0
        
        for length in sorted_lengths:
            cumulative += length
            if cumulative >= target:
                return length
        
        return sorted_lengths[-1]
    
    def _calculate_lx(self, sorted_lengths: List[int], total_length: int, x: int) -> int:
        """Calculate Lx (L50, L75, L90)"""
        if not sorted_lengths:
            return 0
        
        target = total_length * (x / 100)
        cumulative = 0
        
        for i, length in enumerate(sorted_lengths, 1):
            cumulative += length
            if cumulative >= target:
                return i
        
        return len(sorted_lengths)
    
    def _create_length_bins(self, lengths: List[int]) -> Dict[str, int]:
        """Create length distribution bins"""
        bins = {
            '< 100 bp': 0,
            '100-500 bp': 0,
            '500-1k bp': 0,
            '1k-5k bp': 0,
            '5k-10k bp': 0,
            '10k-50k bp': 0,
            '50k-100k bp': 0,
            '100k-500k bp': 0,
            '500k-1M bp': 0,
            '> 1M bp': 0
        }
        
        for length in lengths:
            if length < 100:
                bins['< 100 bp'] += 1
            elif length < 500:
                bins['100-500 bp'] += 1
            elif length < 1000:
                bins['500-1k bp'] += 1
            elif length < 5000:
                bins['1k-5k bp'] += 1
            elif length < 10000:
                bins['5k-10k bp'] += 1
            elif length < 50000:
                bins['10k-50k bp'] += 1
            elif length < 100000:
                bins['50k-100k bp'] += 1
            elif length < 500000:
                bins['100k-500k bp'] += 1
            elif length < 1000000:
                bins['500k-1M bp'] += 1
            else:
                bins['> 1M bp'] += 1
        
        return bins
    
    def _generate_warnings(self, gc_percent: float, ambiguous_percent: float, 
                          max_n_run: int, max_homopolymer: int,
                          short_sequences: int, long_sequences: int, 
                          duplicate_sequences: int, total_sequences: int,
                          s_aureus_status: Dict) -> List[Dict]:
        """Generate warning messages for S. aureus"""
        warnings = []
        
        # GC content for S. aureus
        low, high = self.thresholds['gc_normal']
        if gc_percent < low:
            warnings.append({
                'level': 'warning',
                'message': f'Low GC content ({gc_percent:.1f}%) - below S. aureus range ({low}-{high}%)'
            })
        elif gc_percent > high:
            warnings.append({
                'level': 'warning',
                'message': f'High GC content ({gc_percent:.1f}%) - above S. aureus range ({low}-{high}%)'
            })
        
        # Ambiguous bases
        if ambiguous_percent > self.thresholds['ambiguous_critical']:
            warnings.append({
                'level': 'danger',
                'message': f'High ambiguous bases ({ambiguous_percent:.2f}%) - may indicate poor quality'
            })
        elif ambiguous_percent > self.thresholds['ambiguous_warning']:
            warnings.append({
                'level': 'warning',
                'message': f'Elevated ambiguous bases ({ambiguous_percent:.2f}%)'
            })
        
        # N-runs
        if max_n_run > 100:
            warnings.append({
                'level': 'danger',
                'message': f'Very long N-run detected ({max_n_run} bases) - may indicate assembly gaps'
            })
        elif max_n_run > 10:
            warnings.append({
                'level': 'warning',
                'message': f'Long N-run detected ({max_n_run} bases)'
            })
        
        # Homopolymers
        if max_homopolymer > 20:
            warnings.append({
                'level': 'danger',
                'message': f'Very long homopolymer ({max_homopolymer} bases) - may cause sequencing errors'
            })
        elif max_homopolymer > 10:
            warnings.append({
                'level': 'warning',
                'message': f'Long homopolymer ({max_homopolymer} bases)'
            })
        
        # Short sequences
        if short_sequences > total_sequences * 0.5:
            warnings.append({
                'level': 'danger',
                'message': f'Many short sequences ({short_sequences}) - may indicate poor assembly'
            })
        elif short_sequences > total_sequences * 0.1:
            warnings.append({
                'level': 'warning',
                'message': f'Some short sequences ({short_sequences})'
            })
        
        # Long sequences
        if long_sequences > 0:
            warnings.append({
                'level': 'warning',
                'message': f'Very long sequences detected ({long_sequences}) - may indicate contamination'
            })
        
        # Duplicate sequences
        duplicate_percent = (duplicate_sequences / total_sequences) * 100 if total_sequences > 0 else 0
        if duplicate_percent > 10:
            warnings.append({
                'level': 'warning',
                'message': f'Duplicate sequences detected ({duplicate_sequences}, {duplicate_percent:.1f}%)'
            })
        
        # S. aureus specific warnings
        if s_aureus_status.get('genome_size_status') == 'Atypical':
            warnings.append({
                'level': 'warning',
                'message': f'Atypical genome size ({s_aureus_status["genome_size_mbp"]:.2f} Mbp) for S. aureus (expected: {s_aureus_status["expected_genome_size"]})'
            })
        
        if s_aureus_status.get('assembly_quality') == 'Fragmented':
            warnings.append({
                'level': 'warning',
                'message': f'High contig count ({s_aureus_status["contig_count"]}) - consider additional assembly polishing'
            })
        
        return warnings
    
    def create_individual_html_report(self, results: Dict[str, Any], output_dir: str) -> str:
        """Create comprehensive HTML report for a single FASTA file in StaphScope style"""
        # Create individual folder for this FASTA file
        filename_no_ext = Path(results['filename']).stem
        sample_dir = os.path.join(output_dir, filename_no_ext)
        os.makedirs(sample_dir, exist_ok=True)
        
        html_file = os.path.join(sample_dir, f"{filename_no_ext}_fasta_qc_report.html")
        
        # Get a random quote
        random_quote = self.get_random_quote()
        
        # Extract variables for template
        sample = results['filename']
        total_sequences = results['total_sequences']
        total_length = results['total_length']
        n50 = results['n50']
        gc_percent = results['gc_percent']
        at_percent = results['at_percent']
        ambiguous_percent = results['ambiguous_percent']
        s_aureus_status = results.get('s_aureus_status', {})
        
        # Build warnings HTML
        warnings_html = ''
        if results.get('warnings'):
            for warning in results['warnings']:
                level_color = '#dc2626' if warning['level'] == 'danger' else '#f59e0b'
                warnings_html += f'''                <div style="background: rgba(220, 38, 38, 0.1); padding: 10px; border-radius: 6px; margin: 5px 0; border-left: 4px solid {level_color};">
                    <div style="font-weight: bold; color: {level_color};">{warning['level'].upper()}</div>
                    <div>{warning['message']}</div>
                </div>
'''
        
        # Build nucleotide composition HTML
        composition_html = f'''
                <div style="display: grid; grid-template-columns: repeat(auto-fit, minmax(150px, 1fr)); gap: 15px; margin-top: 15px;">
                    <div style="background: linear-gradient(135deg, #10b981 0%, #059669 100%); padding: 15px; border-radius: 8px; text-align: center; color: white;">
                        <div style="font-size: 20px; font-weight: bold;">{results['gc_percent']:.1f}%</div>
                        <div style="font-size: 12px; opacity: 0.9;">GC Content</div>
                    </div>
                    <div style="background: linear-gradient(135deg, #3b82f6 0%, #1e40af 100%); padding: 15px; border-radius: 8px; text-align: center; color: white;">
                        <div style="font-size: 20px; font-weight: bold;">{results['at_percent']:.1f}%</div>
                        <div style="font-size: 12px; opacity: 0.9;">AT Content</div>
                    </div>
                    <div style="background: #f0f9ff; padding: 15px; border-radius: 8px; text-align: center; color: #0369a1; border: 1px solid #bae6fd;">
                        <div style="font-size: 20px; font-weight: bold;">{results['a_percent']:.1f}%</div>
                        <div style="font-size: 12px; opacity: 0.9;">Adenine (A)</div>
                    </div>
                    <div style="background: #f0f9ff; padding: 15px; border-radius: 8px; text-align: center; color: #0369a1; border: 1px solid #bae6fd;">
                        <div style="font-size: 20px; font-weight: bold;">{results['t_percent']:.1f}%</div>
                        <div style="font-size: 12px; opacity: 0.9;">Thymine (T)</div>
                    </div>
                    <div style="background: #f0f9ff; padding: 15px; border-radius: 8px; text-align: center; color: #0369a1; border: 1px solid #bae6fd;">
                        <div style="font-size: 20px; font-weight: bold;">{results['g_percent']:.1f}%</div>
                        <div style="font-size: 12px; opacity: 0.9;">Guanine (G)</div>
                    </div>
                    <div style="background: #f0f9ff; padding: 15px; border-radius: 8px; text-align: center; color: #0369a1; border: 1px solid #bae6fd;">
                        <div style="font-size: 20px; font-weight: bold;">{results['c_percent']:.1f}%</div>
                        <div style="font-size: 12px; opacity: 0.9;">Cytosine (C)</div>
                    </div>
                </div>
'''
        
        # Build length distribution HTML
        length_dist_html = ''
        for length_range, count in results['length_distribution'].items():
            percentage = (count / total_sequences * 100) if total_sequences > 0 else 0
            length_dist_html += f'''
                        <tr>
                            <td>{length_range}</td>
                            <td>{count:,}</td>
                            <td>{percentage:.1f}%</td>
                        </tr>
'''
        
        # Build S. aureus status HTML
        s_aureus_html = ''
        if s_aureus_status:
            genome_color = '#10b981' if s_aureus_status['genome_size_status'] == 'Normal' else '#f59e0b'
            gc_color = '#10b981' if s_aureus_status['gc_status'] == 'Normal' else '#f59e0b'
            assembly_color = '#10b981' if s_aureus_status['assembly_quality'] == 'Good' else '#f59e0b' if s_aureus_status['assembly_quality'] == 'Moderate' else '#dc2626'
            
            s_aureus_html = f'''
            <div style="margin-top: 20px; padding: 20px; background: linear-gradient(135deg, #1e293b 0%, #0f172a 100%); border-radius: 10px; color: white;">
                <h3 style="color: white; margin-bottom: 15px;">🦠 Staphylococcus aureus Specific Analysis</h3>
                <div style="display: grid; grid-template-columns: repeat(auto-fit, minmax(200px, 1fr)); gap: 15px;">
                    <div>
                        <div style="font-size: 12px; opacity: 0.8;">Genome Size</div>
                        <div style="font-size: 18px; font-weight: bold; color: {genome_color};">{s_aureus_status['genome_size_mbp']:.2f} Mbp</div>
                        <div style="font-size: 11px; opacity: 0.7;">Status: {s_aureus_status['genome_size_status']}</div>
                    </div>
                    <div>
                        <div style="font-size: 12px; opacity: 0.8;">GC Content</div>
                        <div style="font-size: 18px; font-weight: bold; color: {gc_color};">{results['gc_percent']:.1f}%</div>
                        <div style="font-size: 11px; opacity: 0.7;">Expected: {s_aureus_status['expected_gc_range']}</div>
                    </div>
                    <div>
                        <div style="font-size: 12px; opacity: 0.8;">Assembly Quality</div>
                        <div style="font-size: 18px; font-weight: bold; color: {assembly_color};">{s_aureus_status['assembly_quality']}</div>
                        <div style="font-size: 11px; opacity: 0.7;">Contigs: {s_aureus_status['contig_count']}</div>
                    </div>
                    <div>
                        <div style="font-size: 12px; opacity: 0.8;">Expected Size</div>
                        <div style="font-size: 18px; font-weight: bold;">{s_aureus_status['expected_genome_size']}</div>
                        <div style="font-size: 11px; opacity: 0.7;">Typical S. aureus</div>
                    </div>
                </div>
            </div>
'''
        
        html_content = f'''<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>STAPHSCOPE - FASTA QC Report</title>
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
            max-width: 1600px;
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
        
        .stat-table {{
            width: 100%;
            border-collapse: collapse;
            margin: 20px 0;
            background: white;
            border-radius: 8px;
            overflow: hidden;
            box-shadow: 0 2px 10px rgba(0,0,0,0.1);
        }}
        
        .stat-table th {{
            background: linear-gradient(135deg, #3b82f6 0%, #1e40af 100%);
            color: white;
            padding: 12px;
            text-align: left;
            font-weight: bold;
        }}
        
        .stat-table td {{
            padding: 10px;
            border-bottom: 1px solid #e5e7eb;
        }}
        
        .stat-table tr:nth-child(even) {{
            background-color: #f8fafc;
        }}
        
        .stat-table tr:hover {{
            background-color: #e0f2fe;
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
        
        .controls {{
            background: #f8f9fa;
            padding: 15px;
            border-radius: 8px;
            margin: 15px 0;
            display: flex;
            flex-wrap: wrap;
            gap: 10px;
            align-items: center;
        }}
        
        .export-buttons {{
            display: flex;
            gap: 8px;
            flex-wrap: wrap;
        }}
        
        .export-buttons button {{
            padding: 8px 16px;
            background: #3b82f6;
            color: white;
            border: none;
            border-radius: 4px;
            cursor: pointer;
            font-size: 14px;
            transition: background 0.3s;
        }}
        
        .export-buttons button:hover {{
            background: #2563eb;
        }}
        
        .export-buttons button.print {{
            background: #10b981;
        }}
        
        .export-buttons button.print:hover {{
            background: #059669;
        }}
        
        .export-buttons button.json {{
            background: #8b5cf6;
        }}
        
        .export-buttons button.json:hover {{
            background: #7c3aed;
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
                    <div class="metric-value">{results['analysis_date']}</div>
                </div>
                <div class="metric-card">
                    <div class="metric-label">File Size</div>
                    <div class="metric-value">{results['file_size_mb']:.2f} MB</div>
                </div>
            </div>
        </div>
        
        <div class="report-section">
            <h2>🎯 FASTA QC Summary</h2>
            <div class="metrics-grid">
                <div class="metric-card">
                    <div class="metric-label">Total Sequences</div>
                    <div class="metric-value">{total_sequences:,}</div>
                </div>
                <div class="metric-card">
                    <div class="metric-label">Total Length</div>
                    <div class="metric-value">{total_length:,}</div>
                    <div class="metric-label">base pairs</div>
                </div>
                <div class="metric-card">
                    <div class="metric-label">N50</div>
                    <div class="metric-value">{n50:,}</div>
                    <div class="metric-label">base pairs</div>
                </div>
                <div class="metric-card">
                    <div class="metric-label">GC Content</div>
                    <div class="metric-value">{gc_percent:.1f}%</div>
                    <div class="metric-label">S. aureus: 32-38%</div>
                </div>
                <div class="metric-card">
                    <div class="metric-label">AT Content</div>
                    <div class="metric-value">{at_percent:.1f}%</div>
                </div>
                <div class="metric-card">
                    <div class="metric-label">Ambiguous Bases</div>
                    <div class="metric-value">{ambiguous_percent:.2f}%</div>
                </div>
            </div>
        </div>
        
        <div class="report-section">
            <h2>📈 Basic Statistics</h2>
            <table class="stat-table">
                <thead>
                    <tr>
                        <th>Metric</th>
                        <th>Value</th>
                        <th>Description</th>
                    </tr>
                </thead>
                <tbody>
                    <tr>
                        <td><strong>Total Sequences</strong></td>
                        <td>{results['total_sequences']:,}</td>
                        <td>Number of sequences in the file</td>
                    </tr>
                    <tr>
                        <td><strong>Total Length</strong></td>
                        <td>{results['total_length']:,} bp</td>
                        <td>Total number of bases</td>
                    </tr>
                    <tr>
                        <td><strong>Total Bases</strong></td>
                        <td>{results['total_bases']:,} bp</td>
                        <td>Total bases excluding ambiguous characters</td>
                    </tr>
                    <tr>
                        <td><strong>Longest Sequence</strong></td>
                        <td>{results['longest_sequence']:,} bp</td>
                        <td>Length of the longest sequence</td>
                    </tr>
                    <tr>
                        <td><strong>Shortest Sequence</strong></td>
                        <td>{results['shortest_sequence']:,} bp</td>
                        <td>Length of the shortest sequence</td>
                    </tr>
                    <tr>
                        <td><strong>Mean Length</strong></td>
                        <td>{results['mean_length']:,.0f} bp</td>
                        <td>Average sequence length</td>
                    </tr>
                    <tr>
                        <td><strong>Median Length</strong></td>
                        <td>{results['median_length']:,} bp</td>
                        <td>Median sequence length</td>
                    </tr>
                    <tr>
                        <td><strong>N50</strong></td>
                        <td>{results['n50']:,} bp</td>
                        <td>Length for which 50% of total bases are in longer sequences</td>
                    </tr>
                    <tr>
                        <td><strong>N75</strong></td>
                        <td>{results['n75']:,} bp</td>
                        <td>Length for which 75% of total bases are in longer sequences</td>
                    </tr>
                    <tr>
                        <td><strong>N90</strong></td>
                        <td>{results['n90']:,} bp</td>
                        <td>Length for which 90% of total bases are in longer sequences</td>
                    </tr>
                    <tr>
                        <td><strong>L50</strong></td>
                        <td>{results['l50']:,}</td>
                        <td>Number of sequences that make up 50% of total length</td>
                    </tr>
                    <tr>
                        <td><strong>L75</strong></td>
                        <td>{results['l75']:,}</td>
                        <td>Number of sequences that make up 75% of total length</td>
                    </tr>
                    <tr>
                        <td><strong>L90</strong></td>
                        <td>{results['l90']:,}</td>
                        <td>Number of sequences that make up 90% of total length</td>
                    </tr>
                </tbody>
            </table>
        </div>
        
        <div class="report-section">
            <h2>🧬 Nucleotide Composition</h2>
            {composition_html}
        </div>
        
        <div class="report-section">
            <h2>📊 Length Distribution</h2>
            <table class="stat-table">
                <thead>
                    <tr>
                        <th>Length Range</th>
                        <th>Count</th>
                        <th>Percentage</th>
                    </tr>
                </thead>
                <tbody>
                    {length_dist_html}
                </tbody>
            </table>
        </div>
        
        {s_aureus_html}
        
        <div class="report-section">
            <h2>⚠️ Quality Warnings</h2>
            {warnings_html if warnings_html else '<p style="color: #10b981; font-weight: bold;">✅ No quality warnings detected</p>'}
        </div>
        
        <div class="report-section">
            <h2>📋 Additional Statistics</h2>
            <div style="display: grid; grid-template-columns: repeat(auto-fit, minmax(200px, 1fr)); gap: 15px; margin-top: 15px;">
                <div style="background: #f0f9ff; padding: 15px; border-radius: 8px; color: #0369a1; border: 1px solid #bae6fd;">
                    <div style="font-size: 16px; font-weight: bold;">{results['sequences_with_n']:,}</div>
                    <div style="font-size: 12px;">Sequences with Ns</div>
                </div>
                <div style="background: #f0f9ff; padding: 15px; border-radius: 8px; color: #0369a1; border: 1px solid #bae6fd;">
                    <div style="font-size: 16px; font-weight: bold;">{results['total_n_bases']:,}</div>
                    <div style="font-size: 12px;">Total N Bases</div>
                </div>
                <div style="background: #f0f9ff; padding: 15px; border-radius: 8px; color: #0369a1; border: 1px solid #bae6fd;">
                    <div style="font-size: 16px; font-weight: bold;">{results['max_n_run']:,}</div>
                    <div style="font-size: 12px;">Max N-run</div>
                </div>
                <div style="background: #f0f9ff; padding: 15px; border-radius: 8px; color: #0369a1; border: 1px solid #bae6fd;">
                    <div style="font-size: 16px; font-weight: bold;">{results['homopolymers_count']:,}</div>
                    <div style="font-size: 12px;">Homopolymers</div>
                </div>
                <div style="background: #f0f9ff; padding: 15px; border-radius: 8px; color: #0369a1; border: 1px solid #bae6fd;">
                    <div style="font-size: 16px; font-weight: bold;">{results['max_homopolymer']:,}</div>
                    <div style="font-size: 12px;">Max Homopolymer</div>
                </div>
                <div style="background: #f0f9ff; padding: 15px; border-radius: 8px; color: #0369a1; border: 1px solid #bae6fd;">
                    <div style="font-size: 16px; font-weight: bold;">{results['duplicate_sequences']:,}</div>
                    <div style="font-size: 12px;">Duplicate Sequences</div>
                </div>
            </div>
        </div>
        
        <div class="footer">
            <p><strong>STAPHSCOPE</strong> - FASTA QC Analysis Module</p>
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
        const quotes = {json.dumps(self.science_quotes)};
        let currentQuote = 0;

        function getRandomQuote() {{
            return quotes[Math.floor(Math.random() * quotes.length)];
        }}

        function displayQuote() {{
            const quoteContainer = document.getElementById('quoteContainer');
            const quoteText = document.getElementById('quoteText');
            const quoteAuthor = document.getElementById('quoteAuthor');
            
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
        
        // Print report
        function printReport() {{
            window.print();
        }}
        
        // Export to JSON
        function exportToJSON() {{
            const reportData = {json.dumps(results, indent=2)};
            const dataStr = JSON.stringify(reportData, null, 2);
            const dataBlob = new Blob([dataStr], {{ type: 'application/json' }});
            const url = URL.createObjectURL(dataBlob);
            const a = document.createElement('a');
            a.href = url;
            a.download = '{filename_no_ext}_fasta_qc_report.json';
            document.body.appendChild(a);
            a.click();
            document.body.removeChild(a);
            URL.revokeObjectURL(url);
        }}
    </script>
</body>
</html>'''
        
        # Write HTML file
        with open(html_file, 'w', encoding='utf-8') as f:
            f.write(html_content)
        
        # Create JSON report
        json_file = os.path.join(sample_dir, f"{filename_no_ext}_fasta_qc_report.json")
        with open(json_file, 'w', encoding='utf-8') as f:
            json.dump(results, f, indent=2, default=str)
        
        self.logger.info(f"✅ HTML report generated: {html_file}")
        self.logger.info(f"✅ JSON report generated: {json_file}")
        
        return html_file
    
    def create_summary_report(self, all_results: List[Dict[str, Any]], output_dir: str):
        """Create summary report for multiple FASTA files in StaphScope style"""
        successful_results = [r for r in all_results if r.get('status') == 'success']
        
        if not successful_results:
            self.logger.warning("No successful analyses to create summary report")
            return
        
        # Create TSV summary
        tsv_file = os.path.join(output_dir, "FASTA_QC_summary.tsv")
        with open(tsv_file, 'w', encoding='utf-8') as f:
            # Header
            headers = [
                'Filename', 'Total Sequences', 'Total Length', 'Total Bases',
                'GC Content (%)', 'AT Content (%)', 'N50', 'N75', 'N90',
                'Median Length', 'Mean Length', 'Longest Sequence', 'Shortest Sequence',
                'Ambiguous Bases (%)', 'Sequences with Ns', 'Max N-run',
                'Homopolymers', 'Max Homopolymer', 'Duplicate Sequences',
                'Short Sequences (<100 bp)', 'Long Sequences (>1M bp)',
                'File Size (MB)', 'Warnings', 'S. aureus Status'
            ]
            f.write('\t'.join(headers) + '\n')
            
            for result in successful_results:
                s_aureus_status = result.get('s_aureus_status', {})
                status_str = f"GC:{s_aureus_status.get('gc_status', 'Unknown')},Size:{s_aureus_status.get('genome_size_status', 'Unknown')}"
                
                row = [
                    result['filename'],
                    str(result['total_sequences']),
                    str(result['total_length']),
                    str(result['total_bases']),
                    f"{result['gc_percent']:.2f}",
                    f"{result['at_percent']:.2f}",
                    str(result['n50']),
                    str(result['n75']),
                    str(result['n90']),
                    str(result['median_length']),
                    f"{result['mean_length']:.0f}",
                    str(result['longest_sequence']),
                    str(result['shortest_sequence']),
                    f"{result['ambiguous_percent']:.2f}",
                    str(result['sequences_with_n']),
                    str(result['max_n_run']),
                    str(result['homopolymers_count']),
                    str(result['max_homopolymer']),
                    str(result['duplicate_sequences']),
                    str(result['short_sequences']),
                    str(result['long_sequences']),
                    f"{result['file_size_mb']:.2f}",
                    str(len(result.get('warnings', []))),
                    status_str
                ]
                f.write('\t'.join(row) + '\n')
        
        # Create JSON summary
        json_summary = self._create_json_summary(successful_results)
        json_file = os.path.join(output_dir, "FASTA_QC_summary.json")
        with open(json_file, 'w', encoding='utf-8') as f:
            json.dump(json_summary, f, indent=2, default=str)
        
        # Create HTML summary report
        html_file = os.path.join(output_dir, "FASTA_QC_summary.html")
        self._create_summary_html_report(successful_results, html_file, json_summary)
        
        self.logger.info(f"✅ TSV summary created: {tsv_file}")
        self.logger.info(f"✅ JSON summary created: {json_file}")
        self.logger.info(f"✅ HTML summary created: {html_file}")
    
    def _create_json_summary(self, successful_results: List[Dict[str, Any]]) -> Dict[str, Any]:
        """Create JSON summary structure"""
        summary_data = {
            'metadata': {
                'tool': self.metadata['tool_name'],
                'version': self.metadata['version'],
                'biopython_version': self.metadata['biopython_version'],
                'analysis_date': self.metadata['analysis_date'],
                'total_files': len(successful_results),
                'total_sequences': sum(r['total_sequences'] for r in successful_results),
                'total_length': sum(r['total_length'] for r in successful_results)
            },
            'statistics': {
                'gc_content_range': {
                    'min': min(r['gc_percent'] for r in successful_results) if successful_results else 0,
                    'max': max(r['gc_percent'] for r in successful_results) if successful_results else 0,
                    'mean': statistics.mean(r['gc_percent'] for r in successful_results) if successful_results else 0,
                    'median': statistics.median(r['gc_percent'] for r in successful_results) if successful_results else 0,
                },
                'n50_range': {
                    'min': min(r['n50'] for r in successful_results) if successful_results else 0,
                    'max': max(r['n50'] for r in successful_results) if successful_results else 0,
                    'mean': statistics.mean(r['n50'] for r in successful_results) if successful_results else 0,
                    'median': statistics.median(r['n50'] for r in successful_results) if successful_results else 0,
                },
                'files_with_warnings': sum(1 for r in successful_results if r.get('warnings')),
                'total_warnings': sum(len(r.get('warnings', [])) for r in successful_results),
                'files_with_high_gc': sum(1 for r in successful_results if r['gc_percent'] > self.thresholds['gc_normal'][1]),
                'files_with_low_gc': sum(1 for r in successful_results if r['gc_percent'] < self.thresholds['gc_normal'][0]),
            },
            'files': successful_results
        }
        
        return summary_data
    
    def _create_summary_html_report(self, successful_results: List[Dict[str, Any]], 
                                  html_file: str, json_summary: Dict[str, Any]):
        """Create HTML summary report with StaphScope styling"""
        
        # Get a random quote
        random_quote = self.get_random_quote()
        
        # Build table rows
        table_rows = ''
        for result in successful_results:
            warning_count = len(result.get('warnings', []))
            warning_class = 'none' if warning_count == 0 else 'low' if warning_count < 3 else 'high'
            
            table_rows += f'''
                        <tr>
                            <td><strong>{result['filename']}</strong></td>
                            <td>{result['total_sequences']:,}</td>
                            <td>{result['total_length']:,}</td>
                            <td>{result['total_bases']:,}</td>
                            <td>{result['gc_percent']:.1f}</td>
                            <td>{result['at_percent']:.1f}</td>
                            <td>{result['n50']:,}</td>
                            <td>{result['n75']:,}</td>
                            <td>{result['n90']:,}</td>
                            <td>{result['median_length']:,}</td>
                            <td>{result['mean_length']:.0f}</td>
                            <td>{result['longest_sequence']:,}</td>
                            <td>{result['shortest_sequence']:,}</td>
                            <td>{result['ambiguous_percent']:.2f}</td>
                            <td>{result['sequences_with_n']:,}</td>
                            <td>{result['max_n_run']:,}</td>
                            <td>{result['homopolymers_count']:,}</td>
                            <td>{result['max_homopolymer']:,}</td>
                            <td>{result['duplicate_sequences']:,}</td>
                            <td>{result['short_sequences']:,}</td>
                            <td>{result['long_sequences']:,}</td>
                            <td>{result['file_size_mb']:.2f}</td>
                            <td><span style="display: inline-block; padding: 2px 8px; border-radius: 12px; font-size: 11px; font-weight: bold; background: {'#d4edda' if warning_class == 'none' else '#fff3cd' if warning_class == 'low' else '#f8d7da'}; color: {'#155724' if warning_class == 'none' else '#856404' if warning_class == 'low' else '#721c24'};">{warning_count}</span></td>
                        </tr>
'''
        
        # Overall statistics
        total_sequences = sum(r['total_sequences'] for r in successful_results)
        total_length = sum(r['total_length'] for r in successful_results)
        total_warnings = sum(len(r.get('warnings', [])) for r in successful_results)
        files_with_warnings = sum(1 for r in successful_results if r.get('warnings'))
        
        html_content = f'''<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>STAPHSCOPE - FASTA QC Summary Report</title>
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
            font-size: 12px;
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
        
        .table-container {{
            overflow-x: auto;
            margin: 20px 0;
            max-height: 600px;
            overflow-y: auto;
        }}
        
        @media (max-width: 768px) {{
            .ascii-art {{
                font-size: 6px;
            }}
            .summary-table {{
                font-size: 11px;
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
            <h2>📊 FASTA QC Summary - All Samples</h2>
            
            <div class="stats-grid">
                <div class="stat-card">
                    <div class="stat-value">{len(successful_results)}</div>
                    <div class="stat-label">SAMPLES PROCESSED</div>
                </div>
                <div class="stat-card">
                    <div class="stat-value">{total_sequences:,}</div>
                    <div class="stat-label">TOTAL SEQUENCES</div>
                </div>
                <div class="stat-card">
                    <div class="stat-value">{total_length:,}</div>
                    <div class="stat-label">TOTAL BASES</div>
                </div>
                <div class="stat-card">
                    <div class="stat-value">{total_warnings}</div>
                    <div class="stat-label">TOTAL WARNINGS</div>
                </div>
            </div>
            
            <h3>Staphylococcus aureus Analysis Statistics</h3>
            <p><strong>Analysis Date:</strong> {self.metadata['analysis_date']}</p>
            <p><strong>Tool Version:</strong> {self.metadata['version']}</p>
            <p><strong>BioPython Version:</strong> {self.metadata['biopython_version']}</p>
        </div>
        
        <div class="report-section">
            <h2>📈 Detailed FASTA QC Results</h2>
            <p style="color: #666; margin-bottom: 15px;">Scroll horizontally to view all columns.</p>
            
            <div class="table-container">
                <table class="summary-table" id="qc-summary-table">
                    <thead>
                        <tr>
                            <th>Filename</th>
                            <th>Total Sequences</th>
                            <th>Total Length</th>
                            <th>Total Bases</th>
                            <th>GC Content (%)</th>
                            <th>AT Content (%)</th>
                            <th>N50</th>
                            <th>N75</th>
                            <th>N90</th>
                            <th>Median Length</th>
                            <th>Mean Length</th>
                            <th>Longest Sequence</th>
                            <th>Shortest Sequence</th>
                            <th>Ambiguous Bases (%)</th>
                            <th>Sequences with Ns</th>
                            <th>Max N-run</th>
                            <th>Homopolymers</th>
                            <th>Max Homopolymer</th>
                            <th>Duplicate Sequences</th>
                            <th>Short Sequences (&lt;100 bp)</th>
                            <th>Long Sequences (&gt;1M bp)</th>
                            <th>File Size (MB)</th>
                            <th>Warnings</th>
                        </tr>
                    </thead>
                    <tbody>
                        {table_rows}
                    </tbody>
                </table>
            </div>
        </div>
        
        <div class="footer">
            <p><strong>STAPHSCOPE</strong> - FASTA QC Summary Report</p>
            <p class="timestamp">Generated: {datetime.now().strftime("%Y-%m-%d %H:%M:%S")}</p>
            <p>Github: bbeckley-hub | Email: brownbeckley94@gmail.com</p>
        </div>
    </div>

    <script>
        const quotes = {json.dumps(self.science_quotes)};
        let currentQuote = 0;

        function getRandomQuote() {{
            return quotes[Math.floor(Math.random() * quotes.length)];
        }}

        function displayQuote() {{
            const quoteContainer = document.getElementById('quoteContainer');
            const quoteText = document.getElementById('quoteText');
            const quoteAuthor = document.getElementById('quoteAuthor');
            
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
        
        # Write HTML file
        with open(html_file, 'w', encoding='utf-8') as f:
            f.write(html_content)
    
    def process_files(self, pattern: str, output_dir: str = "fasta_qc_results") -> List[Dict[str, Any]]:
        """Process multiple FASTA files with parallel execution and comprehensive reporting"""
        print("\n" + "="*80)
        print("🔬 STAPHSCOPE FASTA QC - Staphylococcus aureus Quality Control")
        print("="*80)
        print(f"Author: Brown Beckley | Email: brownbeckley94@gmail.com")
        print(f"Affiliation: University of Ghana Medical School - Department of Medical Biochemistry")
        print("="*80)
        print(f"Output directory: {output_dir}")
        print(f"Using {self.cpus} CPU cores")
        print("="*80)
        
        # Find FASTA files
        fasta_extensions = ['.fna', '.fasta', '.fa', '.faa', '.fn', '.fna.gz', '.fasta.gz', '.fa.gz']
        files = []
        
        # Try the pattern as-is first
        files.extend(glob.glob(pattern))
        
        # If no files found, try with extensions
        if not files:
            for ext in fasta_extensions:
                files.extend(glob.glob(f"{pattern}{ext}"))
        
        # Remove duplicates and non-existent files
        files = list(set([f for f in files if os.path.exists(f)]))
        
        if not files:
            self.logger.error(f"❌ No FASTA files found matching pattern: {pattern}")
            self.logger.info(f"🔍 Supported extensions: {', '.join(fasta_extensions)}")
            return []
        
        self.logger.info(f"📁 Found {len(files)} FASTA files: {[Path(f).name for f in files]}")
        
        # Create output directory
        os.makedirs(output_dir, exist_ok=True)
        
        # Process files in parallel
        all_results = []
        successful_files = 0
        
        with ThreadPoolExecutor(max_workers=self.cpus) as executor:
            future_to_file = {executor.submit(self.analyze_file, f): f for f in files}
            
            for future in as_completed(future_to_file):
                file = future_to_file[future]
                try:
                    result = future.result()
                    all_results.append(result)
                    
                    if result.get('status') == 'success':
                        successful_files += 1
                        # Create individual HTML report
                        self.create_individual_html_report(result, output_dir)
                    else:
                        self.logger.error(f"❌ {result['filename']}: {result.get('error', 'Unknown error')}")
                        
                except Exception as e:
                    self.logger.error(f"❌ {os.path.basename(file)}: ERROR - {str(e)}")
                    all_results.append({
                        'filename': os.path.basename(file),
                        'status': 'error',
                        'error': str(e)
                    })
        
        # Create summary reports if we have successful analyses
        if successful_files > 0:
            self.create_summary_report(all_results, output_dir)
        
        # Print final summary
        print("\n" + "="*80)
        print("🎉 ANALYSIS COMPLETE")
        print("="*80)
        print(f"Total files: {len(files)}")
        print(f"Successful analyses: {successful_files}")
        print(f"Failed analyses: {len(files) - successful_files}")
        print("\n📁 OUTPUT FILES:")
        print(f"   Individual HTML reports: {output_dir}/*/*_fasta_qc_report.html")
        print(f"   Individual JSON reports: {output_dir}/*/*_fasta_qc_report.json")
        if successful_files > 1:
            print(f"   Summary TSV: {output_dir}/FASTA_QC_summary.tsv")
            print(f"   Summary JSON: {output_dir}/FASTA_QC_summary.json")
            print(f"   Summary HTML: {output_dir}/FASTA_QC_summary.html")
        
        print("\n" + "="*80)
        
        return all_results

def main():
    """Command line interface for FASTA QC analysis"""
    parser = argparse.ArgumentParser(
        description='StaphScope FASTA QC - Staphylococcus aureus Quality Control with HTML Reports',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Run on all FASTA files
  python staph_fasta_qc.py "*.fna"
  
  # Run on specific pattern
  python staph_fasta_qc.py "genomes/*.fasta"
  
  # Run with custom output directory
  python staph_fasta_qc.py "*.fa" --output my_qc_results

Supported FASTA extensions: .fasta, .fa, .fna, .faa, .fn, .gb, .gbk, .gbff
        """
    )
    
    parser.add_argument('pattern', help='File pattern for FASTA files (e.g., "*.fasta", "genomes/*.fna")')
    parser.add_argument('--output', '-o', default='fasta_qc_results', 
                       help='Output directory (default: fasta_qc_results)')
    parser.add_argument('--cpus', '-c', type=int, default=None,
                       help='Number of CPU cores to use (default: all available)')
    
    args = parser.parse_args()
    
    try:
        qc = StaphFASTAQC(cpus=args.cpus)
        results = qc.process_files(args.pattern, args.output)
        
    except Exception as e:
        print(f"❌ FASTA QC analysis failed: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()

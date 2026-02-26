#!/usr/bin/env python3
"""
STAPHSCOPE VISUALIZATION MODULE - ULTIMATE S. AUREUS EDITION
Parses StaphScope HTML reports and generates publication-ready visualizations
Author: Brown Beckley <brownbeckley94@gmail.com>
Affiliation: University of Ghana Medical School
Version: 1.0.0 - Ultimate Staphylococcus Edition
Date: 26-01-2026
"""

import os
import sys
import re
import json
import argparse
import pandas as pd
import numpy as np
from pathlib import Path
from typing import Dict, List, Set, Tuple, Any, Optional, Union
from datetime import datetime
from collections import defaultdict, Counter
import warnings
warnings.filterwarnings('ignore')

# HTML parsing
from bs4 import BeautifulSoup

# Visualization
import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn as sns
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.patches as mpatches
from scipy import stats
import scipy.cluster.hierarchy as sch

# Set publication-quality style
plt.style.use('seaborn-v0_8-whitegrid')
mpl.rcParams['font.family'] = 'sans-serif'
mpl.rcParams['font.sans-serif'] = ['Arial', 'DejaVu Sans', 'Helvetica']
mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42
mpl.rcParams['savefig.dpi'] = 300
mpl.rcParams['savefig.bbox'] = 'tight'
mpl.rcParams['savefig.pad_inches'] = 0.1
mpl.rcParams['figure.max_open_warning'] = 50

# Staphylococcus-specific color palettes
COLOR_PALETTES = {
    'mrsa': ['#DC143C', '#4682B4', '#808080'],  # MRSA, MSSA, Unknown
    'sccmec': ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd',
               '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf',
               '#9edae5', '#dbdb8d'],
    'mlst': ['#4c78a8', '#f58518', '#e45756', '#72b7b2', '#54a24b',
             '#eeca3b', '#b279a2', '#ff9da6', '#9d755d', '#bab0ac',
             '#6a9fb5', '#d4b9da', '#00cc96', '#ab63fa', '#ffa15a'],
    'spa': ['#8dd3c7', '#ffffb3', '#bebada', '#fb8072', '#80b1d3',
            '#fdb462', '#b3de69', '#fccde5', '#d9d9d9', '#bc80bd',
            '#ccebc5', '#ffed6f', '#e78ac3', '#a6d854', '#ffd92f'],
    'database': ['#e41a1c', '#377eb8', '#4daf4a', '#984ea3', '#ff7f00',
                 '#ffff33', '#a65628', '#f781bf', '#999999', '#66c2a5'],
    'risk': ['#DC143C', '#FF4500', '#FF8C00', '#32CD32']  # Critical, High, Medium, Low
}


class StaphHTMLParser:
    """Parses ALL StaphScope HTML reports with BeautifulSoup"""
    
    def __init__(self):
        self.normalized_samples = {}
        
    def normalize_sample_name(self, sample_name: str) -> str:
        """Normalize sample name by removing extensions and paths"""
        sample = str(sample_name)
        
        # Remove common extensions
        extensions = ['.fna', '.fasta', '.fa', '.gb', '.gbk', '.gbff', 
                     '.fna.fna', '.fasta.fasta', '.fa.fa']
        for ext in extensions:
            if sample.endswith(ext):
                sample = sample[:-len(ext)]
        
        # Remove path if present
        if '/' in sample or '\\' in sample:
            sample = Path(sample).name
        
        # Remove any remaining extensions
        for ext in ['.fna', '.fasta', '.fa']:
            if sample.endswith(ext):
                sample = sample[:-len(ext)]
        
        return sample.strip()
    
    def parse_comprehensive_html(self, file_path: Path) -> pd.DataFrame:
        """Parse StaphScope comprehensive typing report"""
        print(f"  📄 Parsing Comprehensive HTML: {file_path.name}")
        
        try:
            with open(file_path, 'r', encoding='utf-8') as f:
                html_content = f.read()
            
            soup = BeautifulSoup(html_content, 'html.parser')
            
            # Find the main typing table
            tables = soup.find_all('table')
            typing_table = None
            
            for table in tables:
                if table.find('th') and any('MLST' in th.get_text() for th in table.find_all('th')):
                    typing_table = table
                    break
            
            if not typing_table:
                # Try to find by text content
                for table in tables:
                    if table.get_text().strip() and 'Sample' in table.get_text() and 'MLST' in table.get_text():
                        typing_table = table
                        break
            
            if not typing_table:
                raise ValueError("Could not find typing table in HTML")
            
            # Parse table headers
            headers = []
            header_row = typing_table.find('tr')
            if header_row:
                for th in header_row.find_all(['th', 'td']):
                    headers.append(th.get_text().strip())
            
            # Parse table rows
            rows = []
            for tr in typing_table.find_all('tr')[1:]:  # Skip header
                cols = tr.find_all(['td', 'th'])
                if cols:
                    row_data = [col.get_text().strip() for col in cols]
                    if len(row_data) >= 2:
                        rows.append(row_data)
            
            # Create DataFrame
            if len(headers) > 0 and len(rows) > 0:
                # Ensure headers match columns
                if len(headers) > len(rows[0]):
                    headers = headers[:len(rows[0])]
                elif len(headers) < len(rows[0]):
                    headers = headers + [f'Col{i}' for i in range(len(headers), len(rows[0]))]
                
                df = pd.DataFrame(rows, columns=headers)
                
                # Clean column names
                df.columns = [col.strip() for col in df.columns]
                
                # Standardize column names
                column_mapping = {
                    'Sample': 'Sample',
                    'sample': 'Sample',
                    'Genome': 'Sample',
                    'MLST': 'MLST',
                    'MLST Type': 'MLST',
                    'ST': 'MLST',
                    'spa Type': 'spa_Type',
                    'spa': 'spa_Type',
                    'SCCmec Type': 'SCCmec_Type',
                    'SCCmec': 'SCCmec_Type',
                    'MRSA/MSSA Status': 'MRSA_Status',
                    'MRSA Status': 'MRSA_Status',
                    'Status': 'MRSA_Status'
                }
                
                for old_col, new_col in column_mapping.items():
                    if old_col in df.columns:
                        df.rename(columns={old_col: new_col}, inplace=True)
                
                # Normalize sample names
                if 'Sample' in df.columns:
                    df['Sample'] = df['Sample'].apply(self.normalize_sample_name)
                
                # Clean MLST values
                if 'MLST' in df.columns:
                    df['MLST'] = df['MLST'].astype(str).str.replace('ST', '').str.strip()
                
                # Clean MRSA status
                if 'MRSA_Status' in df.columns:
                    df['MRSA_Status'] = df['MRSA_Status'].str.upper()
                
                # Keep only essential columns
                keep_cols = ['Sample', 'MLST', 'spa_Type', 'SCCmec_Type', 'MRSA_Status']
                df = df[[col for col in keep_cols if col in df.columns]]
                
                print(f"    ✓ Found {len(df)} samples, {df['MRSA_Status'].value_counts().to_dict()}")
                return df
            
            return pd.DataFrame()
            
        except Exception as e:
            print(f"    ✗ Error parsing comprehensive HTML: {e}")
            import traceback
            traceback.print_exc()
            return pd.DataFrame()
    
    def parse_amrfinder_html(self, file_path: Path) -> Dict[str, pd.DataFrame]:
        """Parse AMRfinder HTML report with dual tables"""
        print(f"  📄 Parsing AMRfinder HTML: {file_path.name}")
        
        try:
            with open(file_path, 'r', encoding='utf-8') as f:
                html_content = f.read()
            
            soup = BeautifulSoup(html_content, 'html.parser')
            
            # Find all tables
            tables = soup.find_all('table')
            if len(tables) < 2:
                print(f"    ⚠️ Only {len(tables)} table(s) found, expected at least 2")
                return {'genes_by_genome': pd.DataFrame(), 'gene_frequency': pd.DataFrame()}
            
            # Parse first table: Genes by Genome
            df1 = self._parse_html_table(tables[1])  # Usually second table (index 1)
            if not df1.empty:
                # Normalize column names
                df1.columns = [col.strip() for col in df1.columns]
                
                # Identify genome column
                genome_col = None
                for col in df1.columns:
                    if 'genome' in col.lower():
                        genome_col = col
                        break
                
                if genome_col:
                    df1['Sample'] = df1[genome_col].apply(self.normalize_sample_name)
                
                # Extract gene count
                count_col = None
                for col in df1.columns:
                    if 'count' in col.lower():
                        count_col = col
                        break
                
                if count_col:
                    df1['Gene_Count'] = pd.to_numeric(df1[count_col], errors='coerce').fillna(0)
                
                # Extract critical genes
                if 'Critical Genes' in df1.columns:
                    df1['Critical_Genes'] = df1['Critical Genes'].fillna('')
                
                # Extract high risk genes
                if 'High Risk Genes' in df1.columns:
                    df1['High_Risk_Genes'] = df1['High Risk Genes'].fillna('')
                
                # Keep essential columns
                keep_cols = ['Sample', 'Gene_Count', 'Critical_Genes', 'High_Risk_Genes']
                df1 = df1[[col for col in keep_cols if col in df1.columns]]
            
            # Parse second table: Gene Frequency
            df2 = self._parse_html_table(tables[2])  # Usually third table (index 2)
            if not df2.empty:
                # Normalize column names
                df2.columns = [col.strip() for col in df2.columns]
                
                # Identify gene column
                gene_col = None
                for col in df2.columns:
                    if 'gene' in col.lower():
                        gene_col = col
                        break
                
                if gene_col:
                    df2['Gene'] = df2[gene_col].astype(str).str.strip()
                
                # Extract frequency and count
                if 'Frequency' in df2.columns:
                    df2['Frequency_Text'] = df2['Frequency']
                    df2['Count'] = df2['Frequency'].apply(self._extract_count_from_frequency)
                
                # Extract prevalence
                if 'Prevalence' in df2.columns:
                    df2['Prevalence'] = df2['Prevalence'].astype(str).str.strip()
                
                # Extract risk level
                if 'Risk Level' in df2.columns:
                    df2['Risk_Level'] = df2['Risk Level'].astype(str).str.strip()
                
                # Extract genomes
                if 'Genomes' in df2.columns:
                    df2['Genomes_Text'] = df2['Genomes']
                
                # Keep essential columns
                keep_cols = ['Gene', 'Count', 'Frequency_Text', 'Prevalence', 'Risk_Level', 'Genomes_Text']
                df2 = df2[[col for col in keep_cols if col in df2.columns]]
            
            result = {
                'genes_by_genome': df1,
                'gene_frequency': df2
            }
            
            print(f"    ✓ AMRfinder: {len(df1)} genome entries, {len(df2)} gene entries")
            return result
            
        except Exception as e:
            print(f"    ✗ Error parsing AMRfinder HTML: {e}")
            import traceback
            traceback.print_exc()
            return {'genes_by_genome': pd.DataFrame(), 'gene_frequency': pd.DataFrame()}
    
    def parse_abricate_html(self, file_path: Path) -> Dict[str, pd.DataFrame]:
        """Parse ABRicate HTML report (for all databases: VFDB, CARD, ResFinder, etc.)"""
        print(f"  📄 Parsing ABRicate HTML: {file_path.name}")
        
        try:
            with open(file_path, 'r', encoding='utf-8') as f:
                html_content = f.read()
            
            soup = BeautifulSoup(html_content, 'html.parser')
            
            # Determine database name from file name
            db_name = file_path.stem.lower()
            
            # Clean up the database name
            if 'vfdb' in db_name:
                db_name = 'vfdb'
            elif 'card' in db_name:
                db_name = 'card'
            elif 'resfinder' in db_name:
                db_name = 'resfinder'
            elif 'plasmidfinder' in db_name:
                db_name = 'plasmidfinder'
            elif 'argannot' in db_name:
                db_name = 'argannot'
            elif 'megares' in db_name:
                db_name = 'megares'
            elif 'ncbi' in db_name:
                db_name = 'ncbi'
            
            # Find all tables
            tables = soup.find_all('table')
            if len(tables) < 2:
                print(f"    ⚠️ Only {len(tables)} table(s) found, expected at least 2")
                return {'database': db_name, 'genes_by_genome': pd.DataFrame(), 'gene_frequency': pd.DataFrame()}
            
            # Parse first table: Genes by Genome
            df1 = self._parse_html_table(tables[0])
            if not df1.empty:
                # Normalize column names
                df1.columns = [col.strip() for col in df1.columns]
                
                # Identify genome column
                genome_col = None
                for col in df1.columns:
                    if 'genome' in col.lower():
                        genome_col = col
                        break
                
                if genome_col:
                    df1['Sample'] = df1[genome_col].apply(self.normalize_sample_name)
                
                # Extract gene count
                count_col = None
                for col in df1.columns:
                    if 'count' in col.lower():
                        count_col = col
                        break
                
                if count_col:
                    df1['Gene_Count'] = pd.to_numeric(df1[count_col], errors='coerce').fillna(0)
                
                # Extract detected genes
                if 'Genes Detected' in df1.columns:
                    df1['Genes_Detected'] = df1['Genes Detected'].fillna('')
                
                # Keep essential columns
                keep_cols = ['Sample', 'Gene_Count', 'Genes_Detected']
                df1 = df1[[col for col in keep_cols if col in df1.columns]]
            
            # Parse second table: Gene Frequency
            df2 = self._parse_html_table(tables[1])
            if not df2.empty:
                # Normalize column names
                df2.columns = [col.strip() for col in df2.columns]
                
                # Identify gene column
                gene_col = None
                for col in df2.columns:
                    if 'gene' in col.lower():
                        gene_col = col
                        break
                
                if gene_col:
                    df2['Gene'] = df2[gene_col].astype(str).str.strip()
                
                # Extract frequency and count
                if 'Frequency' in df2.columns:
                    df2['Frequency_Text'] = df2['Frequency']
                    df2['Count'] = df2['Frequency'].apply(self._extract_count_from_frequency)
                
                # Extract genomes
                if 'Genomes' in df2.columns:
                    df2['Genomes_Text'] = df2['Genomes']
                
                # Keep essential columns
                keep_cols = ['Gene', 'Count', 'Frequency_Text', 'Genomes_Text']
                df2 = df2[[col for col in keep_cols if col in df2.columns]]
            
            result = {
                'database': db_name,
                'genes_by_genome': df1,
                'gene_frequency': df2
            }
            
            print(f"    ✓ {db_name.upper()}: {len(df1)} genome entries, {len(df2)} gene entries")
            return result
            
        except Exception as e:
            print(f"    ✗ Error parsing ABRicate HTML: {e}")
            import traceback
            traceback.print_exc()
            return {'database': 'unknown', 'genes_by_genome': pd.DataFrame(), 'gene_frequency': pd.DataFrame()}
    
    def _parse_html_table(self, table) -> pd.DataFrame:
        """Parse a single HTML table"""
        try:
            # Parse headers
            headers = []
            header_row = table.find('tr')
            if header_row:
                for th in header_row.find_all(['th', 'td']):
                    headers.append(th.get_text().strip())
            
            # Parse rows
            rows = []
            for tr in table.find_all('tr')[1:]:  # Skip header
                cols = tr.find_all(['td', 'th'])
                if cols:
                    row_data = [col.get_text().strip() for col in cols]
                    rows.append(row_data)
            
            # Create DataFrame
            if len(headers) > 0 and len(rows) > 0:
                # Ensure headers match columns
                if len(headers) > len(rows[0]):
                    headers = headers[:len(rows[0])]
                elif len(headers) < len(rows[0]):
                    headers = headers + [f'Col{i}' for i in range(len(headers), len(rows[0]))]
                
                return pd.DataFrame(rows, columns=headers)
            
            return pd.DataFrame()
            
        except Exception as e:
            print(f"    ⚠️ Error parsing table: {e}")
            return pd.DataFrame()
    
    def _extract_count_from_frequency(self, freq_str: str) -> int:
        """Extract count from frequency string like '24 (100.0%)'"""
        try:
            if pd.isna(freq_str):
                return 0
            
            freq_str = str(freq_str)
            # Look for numbers in the string (first number before space or parenthesis)
            match = re.search(r'(\d+)\s*\(', freq_str)
            if match:
                return int(match.group(1))
            
            # Try to find any number
            match = re.search(r'(\d+)', freq_str)
            if match:
                return int(match.group(1))
            
            return 0
        except:
            return 0


class StaphVisualizer:
    """Generates ALL requested visualizations for S. aureus"""
    
    def __init__(self, output_dir: Path):
        self.output_dir = Path(output_dir)
        
        # Create subdirectories
        self.subdirs = {
            'png': self.output_dir / 'PNG',
            'pdf': self.output_dir / 'PDF',
            'svg': self.output_dir / 'SVG',
            'data': self.output_dir / 'DATA'
        }
        
        for subdir in self.subdirs.values():
            subdir.mkdir(parents=True, exist_ok=True)
        
        # Staphylococcus critical genes
        self.critical_amr_genes = {'meca', 'mecA', 'mecc', 'mecC'}
        self.critical_virulence_genes = {
            'luks-pv', 'lukS-PV', 'lukf-pv', 'lukF-PV',  # PVL toxin
            'tst',  # TSST-1
            'sea', 'seb', 'sec', 'sed', 'see',  # Enterotoxins
            'eta', 'etb',  # Exfoliative toxins
        }
    
    def _save_figure(self, fig, name: str, formats: List[str] = ['png', 'pdf', 'svg']):
        """Save figure in multiple formats"""
        for fmt in formats:
            if fmt == 'png':
                fig.savefig(self.subdirs['png'] / f"{name}.png", dpi=300, bbox_inches='tight')
            elif fmt == 'pdf':
                fig.savefig(self.subdirs['pdf'] / f"{name}.pdf", bbox_inches='tight')
            elif fmt == 'svg':
                fig.savefig(self.subdirs['svg'] / f"{name}.svg", bbox_inches='tight')
        
        plt.close(fig)
    
    def _get_colors(self, n_colors: int, palette: str = 'categorical') -> List[str]:
        """Get colors from palette"""
        if palette in COLOR_PALETTES:
            return COLOR_PALETTES[palette][:n_colors]
        else:
            # Use tab20c colormap
            return plt.cm.tab20c(np.linspace(0, 1, n_colors))
    
    def create_mrsa_analysis(self, typing_data: pd.DataFrame):
        """Create MRSA-specific visualizations"""
        print("\n📊 Creating MRSA analysis plots...")
        
        if typing_data.empty or 'MRSA_Status' not in typing_data.columns:
            print("  ⚠️ No MRSA status data available")
            return
        
        # Calculate MRSA distribution
        mrsa_dist = typing_data['MRSA_Status'].value_counts().reset_index()
        mrsa_dist.columns = ['Status', 'Count']
        mrsa_dist['Percentage'] = (mrsa_dist['Count'] / mrsa_dist['Count'].sum() * 100).round(2)
        
        # Save data
        mrsa_dist.to_csv(self.subdirs['data'] / "mrsa_distribution.csv", index=False)
        
        # Get colors
        colors = self._get_colors(len(mrsa_dist), 'mrsa')
        
        # ========== PIE CHART ==========
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(18, 9))
        
        # Create pie chart
        wedges, texts, autotexts = ax1.pie(
            mrsa_dist['Count'],
            labels=mrsa_dist['Status'],
            colors=colors,
            autopct=lambda pct: f'{pct:.1f}%\n({int(pct/100.*sum(mrsa_dist["Count"]))})',
            startangle=90,
            textprops={'fontsize': 11}
        )
        
        ax1.set_title('MRSA vs MSSA Distribution - Pie Chart\nTotal Samples: {}'.format(
            mrsa_dist['Count'].sum()), fontsize=14, fontweight='bold', pad=20)
        
        # Improve label readability
        for text in texts:
            text.set_fontsize(10)
            text.set_fontweight('bold')
        for autotext in autotexts:
            autotext.set_fontsize(10)
            autotext.set_fontweight('bold')
            autotext.set_color('white')
        
        # ========== DONUT CHART ==========
        # Create donut chart
        wedges, texts, autotexts = ax2.pie(
            mrsa_dist['Count'],
            labels=mrsa_dist['Status'],
            colors=colors,
            autopct=lambda pct: f'{pct:.1f}%\n({int(pct/100.*sum(mrsa_dist["Count"]))})',
            startangle=90,
            textprops={'fontsize': 11},
            wedgeprops={'width': 0.6}  # This creates the donut hole
        )
        
        # Add circle in the middle to complete donut look
        centre_circle = plt.Circle((0, 0), 0.3, fc='white')
        ax2.add_artist(centre_circle)
        
        ax2.set_title('MRSA vs MSSA Distribution - Donut Chart\nTotal Samples: {}'.format(
            mrsa_dist['Count'].sum()), fontsize=14, fontweight='bold', pad=20)
        
        plt.suptitle('S. aureus MRSA Analysis', fontsize=16, fontweight='bold', y=0.98)
        plt.tight_layout()
        
        # Save figure
        self._save_figure(fig, "mrsa_distribution")
        print("    ✓ Created mrsa_distribution.[png/pdf/svg]")
        
        # Create additional MRSA analysis if SCCmec data is available
        if 'SCCmec_Type' in typing_data.columns:
            self._create_sccmec_analysis(typing_data)
        
        # Create MLST-MRSA correlation if MLST data is available
        if 'MLST' in typing_data.columns:
            self._create_mlst_mrsa_correlation(typing_data)
    
    def _create_sccmec_analysis(self, typing_data: pd.DataFrame):
        """Create SCCmec type analysis"""
        print("  Creating SCCmec analysis...")
        
        # Filter for MRSA only
        mrsa_data = typing_data[typing_data['MRSA_Status'] == 'MRSA']
        
        if mrsa_data.empty:
            print("    ⚠️ No MRSA data for SCCmec analysis")
            return
        
        # Calculate SCCmec distribution
        sccmec_dist = mrsa_data['SCCmec_Type'].value_counts().reset_index()
        sccmec_dist.columns = ['SCCmec_Type', 'Count']
        sccmec_dist = sccmec_dist[sccmec_dist['SCCmec_Type'] != 'Not Assigned']
        
        if sccmec_dist.empty:
            print("    ⚠️ No SCCmec types found in MRSA samples")
            return
        
        sccmec_dist['Percentage'] = (sccmec_dist['Count'] / sccmec_dist['Count'].sum() * 100).round(2)
        
        # Save data
        sccmec_dist.to_csv(self.subdirs['data'] / "sccmec_distribution.csv", index=False)
        
        # Get colors
        colors = self._get_colors(len(sccmec_dist), 'sccmec')
        
        # Plot
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(20, 10))
        
        # Bar chart
        y_pos = np.arange(len(sccmec_dist))
        bars = ax1.barh(y_pos, sccmec_dist['Count'], color=colors, edgecolor='black', linewidth=0.5)
        ax1.set_yticks(y_pos)
        
        # Create y-tick labels
        y_labels = [f"{row['SCCmec_Type']}\n(n={row['Count']}, {row['Percentage']}%)" 
                   for _, row in sccmec_dist.iterrows()]
        ax1.set_yticklabels(y_labels, fontsize=10)
        ax1.invert_yaxis()  # Highest count on top
        
        # Add count labels on bars
        for bar, count, pct in zip(bars, sccmec_dist['Count'], sccmec_dist['Percentage']):
            width = bar.get_width()
            ax1.text(width + max(sccmec_dist['Count']) * 0.01, bar.get_y() + bar.get_height()/2,
                    f'{count} ({pct}%)', va='center', fontsize=10, fontweight='bold')
        
        ax1.set_xlabel('Number of MRSA Samples', fontsize=12, fontweight='bold')
        ax1.set_title('SCCmec Type Distribution in MRSA\nTotal MRSA: {}'.format(
            sccmec_dist['Count'].sum()), fontsize=14, fontweight='bold', pad=20)
        ax1.grid(True, alpha=0.3, linestyle='--', axis='x')
        
        # Pie chart
        wedges, texts, autotexts = ax2.pie(
            sccmec_dist['Count'],
            labels=sccmec_dist['SCCmec_Type'],
            colors=colors,
            autopct=lambda pct: f'{pct:.1f}%\n({int(pct/100.*sum(sccmec_dist["Count"]))})',
            startangle=90,
            textprops={'fontsize': 9}
        )
        
        ax2.set_title('SCCmec Type Distribution - Pie Chart', 
                     fontsize=14, fontweight='bold', pad=20)
        
        plt.suptitle('SCCmec Typing Analysis for S. aureus MRSA', 
                    fontsize=16, fontweight='bold', y=0.98)
        plt.tight_layout()
        
        # Save figure
        self._save_figure(fig, "sccmec_analysis")
        print("    ✓ Created sccmec_analysis.[png/pdf/svg]")
    
    def _create_mlst_mrsa_correlation(self, typing_data: pd.DataFrame):
        """Create MLST-MRSA correlation plot"""
        print("  Creating MLST-MRSA correlation...")
        
        # Create cross-tabulation
        cross_tab = pd.crosstab(typing_data['MLST'], typing_data['MRSA_Status'])
        
        # Remove 'ND' if exists
        if 'ND' in cross_tab.index:
            cross_tab = cross_tab.drop('ND')
        if 'ND' in cross_tab.columns:
            cross_tab = cross_tab.drop(columns='ND')
        
        # Sort by total count
        cross_tab = cross_tab.loc[cross_tab.sum(axis=1).sort_values(ascending=False).index]
        
        # Plot stacked bar chart
        fig, ax = plt.subplots(figsize=(16, 10))
        
        # Get colors
        colors = self._get_colors(len(cross_tab.columns), 'mrsa')
        
        # Create stacked bar
        bottom = np.zeros(len(cross_tab))
        for i, col in enumerate(cross_tab.columns):
            ax.bar(cross_tab.index, cross_tab[col], bottom=bottom, 
                  label=col, color=colors[i], edgecolor='black', linewidth=0.5)
            bottom += cross_tab[col]
        
        ax.set_xlabel('MLST Sequence Type', fontsize=12, fontweight='bold')
        ax.set_ylabel('Number of Samples', fontsize=12, fontweight='bold')
        ax.set_title('MLST × MRSA Status Relationship\nTotal Samples: {}'.format(
            len(typing_data)), fontsize=14, fontweight='bold', pad=20)
        
        # Rotate x-tick labels for readability
        plt.xticks(rotation=45, ha='right', fontsize=10)
        plt.yticks(fontsize=10)
        
        # Add count labels on top of each bar
        for i, total in enumerate(cross_tab.sum(axis=1)):
            ax.text(i, total + max(cross_tab.sum(axis=1)) * 0.01, f'{int(total)}', 
                   ha='center', va='bottom', fontsize=9, fontweight='bold')
        
        # Add legend
        ax.legend(title='MRSA Status', bbox_to_anchor=(1.05, 1), loc='upper left', 
                 fontsize=10, title_fontsize=11)
        
        # Add grid
        ax.grid(True, alpha=0.3, linestyle='--', axis='y')
        
        plt.tight_layout()
        
        # Save figure
        self._save_figure(fig, "mlst_mrsa_correlation")
        print("    ✓ Created mlst_mrsa_correlation.[png/pdf/svg]")
    
    def create_typing_distributions(self, typing_data: pd.DataFrame):
        """Create distributions for MLST, spa, and SCCmec typing"""
        print("\n📊 Creating typing distribution plots...")
        
        if typing_data.empty:
            print("  ⚠️ No typing data available")
            return
        
        # MLST Distribution
        if 'MLST' in typing_data.columns:
            self._create_distribution_plot(
                typing_data, 'MLST',
                'MLST Sequence Type Distribution',
                'mlst_distribution',
                'mlst'
            )
        
        # spa Type Distribution
        if 'spa_Type' in typing_data.columns:
            self._create_distribution_plot(
                typing_data, 'spa_Type',
                'spa Typing Distribution',
                'spa_distribution',
                'spa'
            )
        
        # SCCmec Type Distribution (all samples)
        if 'SCCmec_Type' in typing_data.columns:
            self._create_distribution_plot(
                typing_data[typing_data['SCCmec_Type'] != 'Not Assigned'], 'SCCmec_Type',
                'SCCmec Type Distribution',
                'sccmec_distribution',
                'sccmec'
            )
    
    def _create_distribution_plot(self, data: pd.DataFrame, column: str, 
                                 title: str, filename: str, palette: str):
        """Create a distribution plot for a specific typing column"""
        if data.empty or column not in data.columns:
            print(f"    ⚠️ No data for {title}")
            return
        
        # Calculate distribution
        distribution = data[column].value_counts().reset_index()
        distribution.columns = [column, 'Count']
        distribution['Percentage'] = (distribution['Count'] / distribution['Count'].sum() * 100).round(2)
        
        # Sort by count
        distribution = distribution.sort_values('Count', ascending=False)
        
        # Save data
        distribution.to_csv(self.subdirs['data'] / f"{filename}.csv", index=False)
        
        # Get colors
        colors = self._get_colors(len(distribution), palette)
        
        # Create figure
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(20, 10))
        
        # Bar chart
        y_pos = np.arange(len(distribution))
        bars = ax1.barh(y_pos, distribution['Count'], color=colors, edgecolor='black', linewidth=0.5)
        ax1.set_yticks(y_pos)
        
        # Create y-tick labels
        y_labels = [f"{row[column]}\n(n={row['Count']}, {row['Percentage']}%)" 
                   for _, row in distribution.iterrows()]
        ax1.set_yticklabels(y_labels, fontsize=9)
        ax1.invert_yaxis()  # Highest count on top
        
        # Add count labels on bars
        for bar, count, pct in zip(bars, distribution['Count'], distribution['Percentage']):
            width = bar.get_width()
            ax1.text(width + max(distribution['Count']) * 0.01, bar.get_y() + bar.get_height()/2,
                    f'{count} ({pct}%)', va='center', fontsize=9, fontweight='bold')
        
        ax1.set_xlabel('Number of Samples', fontsize=12, fontweight='bold')
        ax1.set_title(f'{title} - Bar Chart\nTotal: {distribution["Count"].sum()} samples',
                     fontsize=14, fontweight='bold', pad=20)
        ax1.grid(True, alpha=0.3, linestyle='--', axis='x')
        
        # Pie chart
        wedges, texts, autotexts = ax2.pie(
            distribution['Count'],
            labels=distribution[column],
            colors=colors,
            autopct=lambda pct: f'{pct:.1f}%\n({int(pct/100.*sum(distribution["Count"]))})',
            startangle=90,
            textprops={'fontsize': 8}
        )
        
        ax2.set_title(f'{title} - Pie Chart', 
                     fontsize=14, fontweight='bold', pad=20)
        
        plt.suptitle(title, fontsize=16, fontweight='bold', y=0.98)
        plt.tight_layout()
        
        # Save figure
        self._save_figure(fig, filename)
        print(f"    ✓ Created {filename}.[png/pdf/svg]")
    
    def create_stacked_combinations(self, typing_data: pd.DataFrame):
        """
        Create stacked combination plots for S. aureus:
        1. MLST × SCCmec
        2. MLST × spa Type
        3. MLST × MRSA Status (already done)
        """
        print("\n📊 Creating stacked combination plots...")
        
        if typing_data.empty:
            print("  ⚠️ No typing data for stacked plots")
            return
        
        # Plot 1: MLST × SCCmec (MRSA only)
        if 'MLST' in typing_data.columns and 'SCCmec_Type' in typing_data.columns:
            mrsa_data = typing_data[typing_data['MRSA_Status'] == 'MRSA']
            mrsa_data = mrsa_data[mrsa_data['SCCmec_Type'] != 'Not Assigned']
            
            if not mrsa_data.empty:
                self._create_stacked_plot(
                    mrsa_data,
                    primary_col='MLST',
                    secondary_col='SCCmec_Type',
                    title='MLST × SCCmec Type Relationship (MRSA only)',
                    filename='stacked_mlst_sccmec'
                )
        
        # Plot 2: MLST × spa Type
        if 'MLST' in typing_data.columns and 'spa_Type' in typing_data.columns:
            self._create_stacked_plot(
                typing_data,
                primary_col='MLST',
                secondary_col='spa_Type',
                title='MLST × spa Type Relationship',
                filename='stacked_mlst_spa'
            )
        
        # Plot 3: SCCmec × spa Type (MRSA only)
        if 'SCCmec_Type' in typing_data.columns and 'spa_Type' in typing_data.columns:
            mrsa_data = typing_data[typing_data['MRSA_Status'] == 'MRSA']
            mrsa_data = mrsa_data[mrsa_data['SCCmec_Type'] != 'Not Assigned']
            
            if not mrsa_data.empty:
                self._create_stacked_plot(
                    mrsa_data,
                    primary_col='SCCmec_Type',
                    secondary_col='spa_Type',
                    title='SCCmec × spa Type Relationship (MRSA only)',
                    filename='stacked_sccmec_spa'
                )
    
    def _create_stacked_plot(self, data: pd.DataFrame, primary_col: str, 
                            secondary_col: str, title: str, filename: str):
        """Create a single stacked bar plot"""
        print(f"  Creating {filename}...")
        
        # Create cross-tabulation
        cross_tab = pd.crosstab(data[primary_col], data[secondary_col])
        
        # Remove 'ND' or empty values if they exist
        if 'ND' in cross_tab.columns:
            cross_tab = cross_tab.drop(columns='ND')
        if 'Not Assigned' in cross_tab.columns:
            cross_tab = cross_tab.drop(columns='Not Assigned')
        
        if cross_tab.empty:
            print(f"    ⚠️ No data for {filename}")
            return
        
        # Sort by total count
        cross_tab = cross_tab.loc[cross_tab.sum(axis=1).sort_values(ascending=False).index]
        
        # Sort columns by total
        cross_tab = cross_tab[cross_tab.sum().sort_values(ascending=False).index]
        
        # Plot
        fig, ax = plt.subplots(figsize=(16, 10))
        
        # Get colors
        n_colors = len(cross_tab.columns)
        colors = self._get_colors(n_colors, 'categorical')
        
        # Create stacked bar
        bottom = np.zeros(len(cross_tab))
        for i, col in enumerate(cross_tab.columns):
            ax.bar(cross_tab.index, cross_tab[col], bottom=bottom, 
                  label=str(col), color=colors[i], edgecolor='black', linewidth=0.5)
            bottom += cross_tab[col]
        
        ax.set_xlabel(primary_col, fontsize=12, fontweight='bold')
        ax.set_ylabel('Number of Samples', fontsize=12, fontweight='bold')
        ax.set_title(f'{title}\nTotal Samples: {len(data)}', 
                    fontsize=14, fontweight='bold', pad=20)
        
        # Rotate x-tick labels for readability
        plt.xticks(rotation=45, ha='right', fontsize=10)
        plt.yticks(fontsize=10)
        
        # Add count labels on top of each bar
        for i, total in enumerate(cross_tab.sum(axis=1)):
            ax.text(i, total + max(cross_tab.sum(axis=1)) * 0.01, f'{int(total)}', 
                   ha='center', va='bottom', fontsize=9, fontweight='bold')
        
        # Add legend
        ax.legend(title=secondary_col, bbox_to_anchor=(1.05, 1), loc='upper left', 
                 fontsize=9, title_fontsize=10)
        
        # Add grid
        ax.grid(True, alpha=0.3, linestyle='--', axis='y')
        
        plt.tight_layout()
        
        # Save figure
        self._save_figure(fig, filename)
        print(f"    ✓ Created {filename}.[png/pdf/svg]")
    
    def create_amr_analysis(self, amr_data: Dict[str, pd.DataFrame], 
                           databases_data: Dict[str, Dict[str, pd.DataFrame]]):
        """Create AMR gene analysis visualizations"""
        print("\n📊 Creating AMR gene analysis plots...")
        
        # 1. Critical gene analysis from AMRfinder
        if amr_data and 'gene_frequency' in amr_data and not amr_data['gene_frequency'].empty:
            self._create_critical_gene_analysis(amr_data['gene_frequency'])
        
        # 2. Database comparison
        if databases_data:
            self._create_database_comparison(databases_data)
        
        # 3. Gene frequency distributions
        if databases_data:
            self._create_gene_frequency_analysis(databases_data)
        
        # 4. Risk level analysis (if available)
        if amr_data and 'gene_frequency' in amr_data and 'Risk_Level' in amr_data['gene_frequency'].columns:
            self._create_risk_level_analysis(amr_data['gene_frequency'])
    
    def _create_critical_gene_analysis(self, gene_freq: pd.DataFrame):
        """Create analysis of critical AMR genes (mecA, mecC)"""
        print("  Creating critical gene analysis...")
        
        # Filter for critical genes
        critical_genes = gene_freq[gene_freq['Gene'].str.contains('mec', case=False, na=False)]
        
        if critical_genes.empty:
            print("    ⚠️ No critical AMR genes found")
            return
        
        # Plot
        fig, ax = plt.subplots(figsize=(12, 8))
        
        # Sort by count
        critical_genes = critical_genes.sort_values('Count', ascending=True)
        
        # Create horizontal bar chart
        y_pos = np.arange(len(critical_genes))
        bars = ax.barh(y_pos, critical_genes['Count'], 
                      color=['#DC143C' if 'cri' in str(risk).lower() else '#FF4500' 
                            for risk in critical_genes.get('Risk_Level', ['HIGH']*len(critical_genes))],
                      edgecolor='black', linewidth=1)
        
        ax.set_yticks(y_pos)
        ax.set_yticklabels(critical_genes['Gene'], fontsize=11, fontweight='bold')
        ax.set_xlabel('Number of Genomes', fontsize=12, fontweight='bold')
        ax.set_title('Critical AMR Genes in S. aureus\n(mecA/mecC Distribution)', 
                    fontsize=14, fontweight='bold', pad=20)
        
        # Add count labels
        for bar, count in zip(bars, critical_genes['Count']):
            width = bar.get_width()
            ax.text(width + max(critical_genes['Count']) * 0.01, bar.get_y() + bar.get_height()/2,
                   f'{count} genomes', va='center', fontsize=10, fontweight='bold')
        
        # Add risk level annotations
        for i, (_, row) in enumerate(critical_genes.iterrows()):
            risk_level = row.get('Risk_Level', 'UNKNOWN')
            ax.text(0, i, f'  {risk_level}', va='center', fontsize=9, 
                   fontweight='bold', color='white')
        
        ax.grid(True, alpha=0.3, linestyle='--', axis='x')
        plt.tight_layout()
        
        # Save figure
        self._save_figure(fig, "critical_amr_genes")
        print("    ✓ Created critical_amr_genes.[png/pdf/svg]")
    
    def _create_database_comparison(self, databases_data: Dict[str, Dict[str, pd.DataFrame]]):
        """Create database comparison plots"""
        print("  Creating database comparison...")
        
        # Prepare data for comparison
        comparison_data = []
        
        for db_name, db_data in databases_data.items():
            if 'genes_by_genome' in db_data and not db_data['genes_by_genome'].empty:
                df = db_data['genes_by_genome']
                if 'Gene_Count' in df.columns:
                    for count in df['Gene_Count']:
                        comparison_data.append({
                            'Database': db_name.upper(),
                            'Gene_Count': float(count)
                        })
        
        if not comparison_data:
            print("    ⚠️ No gene count data for database comparison")
            return
        
        # Create DataFrame
        comparison_df = pd.DataFrame(comparison_data)
        
        # Save data
        comparison_df.to_csv(self.subdirs['data'] / "database_comparison.csv", index=False)
        
        # Calculate statistics
        stats_df = comparison_df.groupby('Database')['Gene_Count'].agg([
            'count', 'mean', 'median', 'std', 'min', 'max'
        ]).round(2)
        stats_df.to_csv(self.subdirs['data'] / "database_statistics.csv")
        
        # ========== BOX PLOT ==========
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(20, 10))
        
        # Get colors
        n_databases = len(stats_df)
        colors = self._get_colors(n_databases, 'database')
        
        # Box plot
        box_data_by_db = [comparison_df[comparison_df['Database'] == db]['Gene_Count'].values 
                         for db in stats_df.index]
        
        bp = ax1.boxplot(box_data_by_db, patch_artist=True, showfliers=True, 
                        labels=stats_df.index, medianprops={'color': 'black', 'linewidth': 2})
        
        # Color boxes
        for patch, color in zip(bp['boxes'], colors):
            patch.set_facecolor(color)
            patch.set_alpha(0.7)
        
        ax1.set_xticklabels(stats_df.index, rotation=45, ha='right', fontsize=10)
        ax1.set_ylabel('Number of Genes Detected', fontsize=12, fontweight='bold')
        ax1.set_title('Gene Detection by Database - Box Plot', 
                     fontsize=14, fontweight='bold', pad=20)
        ax1.grid(True, alpha=0.3, linestyle='--')
        
        # ========== VIOLIN PLOT ==========
        sns.violinplot(data=comparison_df, x='Database', y='Gene_Count', 
                      ax=ax2, palette=colors, inner='quartile', cut=0)
        
        # Add swarm plot for individual points
        sns.swarmplot(data=comparison_df, x='Database', y='Gene_Count', 
                     ax=ax2, color='black', alpha=0.5, size=3)
        
        ax2.set_xticklabels(ax2.get_xticklabels(), rotation=45, ha='right', fontsize=10)
        ax2.set_ylabel('Number of Genes Detected', fontsize=12, fontweight='bold')
        ax2.set_title('Gene Detection by Database - Violin Plot', 
                     fontsize=14, fontweight='bold', pad=20)
        ax2.grid(True, alpha=0.3, linestyle='--')
        
        plt.suptitle('S. aureus - Database Comparison Analysis', 
                    fontsize=16, fontweight='bold', y=0.98)
        plt.tight_layout()
        
        # Save figure
        self._save_figure(fig, "database_comparison")
        print("    ✓ Created database_comparison.[png/pdf/svg]")
        
        # Create statistics table
        self._create_statistics_table(stats_df, "database_statistics")
    
    def _create_gene_frequency_analysis(self, databases_data: Dict[str, Dict[str, pd.DataFrame]]):
        """Create gene frequency distribution plots"""
        print("  Creating gene frequency distributions...")
        
        # Extract frequency data
        freq_data = {}
        for db_name, db_data in databases_data.items():
            if 'gene_frequency' in db_data and not db_data['gene_frequency'].empty:
                df = db_data['gene_frequency']
                if 'Count' in df.columns:
                    counts = pd.to_numeric(df['Count'], errors='coerce').dropna().values
                    if len(counts) > 0:
                        freq_data[db_name] = counts
        
        if not freq_data:
            print("    ⚠️ No frequency data found")
            return
        
        # Create subplots
        n_databases = len(freq_data)
        n_cols = min(3, n_databases)
        n_rows = (n_databases + n_cols - 1) // n_cols
        
        fig, axes = plt.subplots(n_rows, n_cols, figsize=(6 * n_cols, 5 * n_rows))
        if n_databases == 1:
            axes = np.array([axes])
        axes = axes.flatten()
        
        # Get colors
        colors = self._get_colors(n_databases, 'database')
        
        for i, (db_name, counts) in enumerate(freq_data.items()):
            if i >= len(axes):
                break
                
            ax = axes[i]
            
            # Plot histogram with KDE
            n, bins, patches = ax.hist(counts, bins=30, alpha=0.7, density=True,
                                      color=colors[i], edgecolor='black', linewidth=0.5)
            
            try:
                # Plot KDE
                if len(counts) > 1:
                    kde = stats.gaussian_kde(counts)
                    x_range = np.linspace(min(counts), max(counts), 1000)
                    ax.plot(x_range, kde(x_range), 'r-', linewidth=2, label='KDE')
            except:
                pass
            
            # Add statistics lines
            mean_val = np.mean(counts)
            median_val = np.median(counts)
            
            ax.axvline(mean_val, color='blue', linestyle='--', linewidth=2, 
                      label=f'Mean: {mean_val:.2f}')
            ax.axvline(median_val, color='green', linestyle='--', linewidth=2,
                      label=f'Median: {median_val:.2f}')
            
            ax.set_xlabel('Gene Frequency (Count)', fontsize=11)
            ax.set_ylabel('Density', fontsize=11)
            ax.set_title(f'{db_name.upper()}\nGene Frequency Distribution', 
                        fontsize=12, fontweight='bold', pad=15)
            ax.legend(fontsize=9)
            ax.grid(True, alpha=0.3, linestyle='--')
            
            # Add statistics text
            stats_text = f'n = {len(counts)}\nMean = {mean_val:.2f}\nMedian = {median_val:.2f}\nStd = {np.std(counts):.2f}'
            ax.text(0.95, 0.95, stats_text, transform=ax.transAxes, fontsize=9,
                   verticalalignment='top', horizontalalignment='right',
                   bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
        
        # Hide unused subplots
        for j in range(i + 1, len(axes)):
            axes[j].axis('off')
        
        plt.suptitle('Gene Frequency Distributions Across Databases', 
                    fontsize=16, fontweight='bold', y=0.98)
        plt.tight_layout()
        
        # Save figure
        self._save_figure(fig, "gene_frequency_distributions")
        print("    ✓ Created gene_frequency_distributions.[png/pdf/svg]")
    
    def _create_risk_level_analysis(self, gene_freq: pd.DataFrame):
        """Create risk level distribution analysis"""
        print("  Creating risk level analysis...")
        
        if 'Risk_Level' not in gene_freq.columns:
            return
        
        # Calculate risk level distribution
        risk_dist = gene_freq['Risk_Level'].value_counts().reset_index()
        risk_dist.columns = ['Risk_Level', 'Count']
        risk_dist['Percentage'] = (risk_dist['Count'] / risk_dist['Count'].sum() * 100).round(2)
        
        # Sort by predefined order
        risk_order = ['CRITICAL', 'HIGH', 'STANDARD', 'LOW', 'MEDIUM']
        risk_dist['Risk_Level'] = pd.Categorical(risk_dist['Risk_Level'], categories=risk_order, ordered=True)
        risk_dist = risk_dist.sort_values('Risk_Level')
        
        # Get colors
        colors = self._get_colors(len(risk_dist), 'risk')
        
        # Plot
        fig, ax = plt.subplots(figsize=(12, 8))
        
        # Create bar chart
        bars = ax.bar(range(len(risk_dist)), risk_dist['Count'], color=colors, 
                     edgecolor='black', linewidth=1)
        
        ax.set_xticks(range(len(risk_dist)))
        ax.set_xticklabels(risk_dist['Risk_Level'], fontsize=11, fontweight='bold')
        ax.set_ylabel('Number of Genes', fontsize=12, fontweight='bold')
        ax.set_title('AMR Gene Risk Level Distribution\n(AMRfinder Analysis)', 
                    fontsize=14, fontweight='bold', pad=20)
        
        # Add count labels
        for bar, count, pct in zip(bars, risk_dist['Count'], risk_dist['Percentage']):
            height = bar.get_height()
            ax.text(bar.get_x() + bar.get_width()/2., height + max(risk_dist['Count']) * 0.01,
                   f'{count} ({pct}%)', ha='center', va='bottom', fontsize=10, fontweight='bold')
        
        ax.grid(True, alpha=0.3, linestyle='--', axis='y')
        plt.tight_layout()
        
        # Save figure
        self._save_figure(fig, "risk_level_distribution")
        print("    ✓ Created risk_level_distribution.[png/pdf/svg]")
    
    def _create_statistics_table(self, stats_df: pd.DataFrame, filename: str):
        """Create a visual table of statistics"""
        fig, ax = plt.subplots(figsize=(12, 8))
        ax.axis('tight')
        ax.axis('off')
        
        # Create table data
        table_data = []
        for db in stats_df.index:
            row = [
                db,
                f"{stats_df.loc[db, 'count']}",
                f"{stats_df.loc[db, 'mean']:.2f}",
                f"{stats_df.loc[db, 'median']:.2f}",
                f"{stats_df.loc[db, 'std']:.2f}",
                f"{stats_df.loc[db, 'min']:.0f}",
                f"{stats_df.loc[db, 'max']:.0f}"
            ]
            table_data.append(row)
        
        # Create table
        table = ax.table(cellText=table_data,
                        colLabels=['Database', 'N', 'Mean', 'Median', 'Std', 'Min', 'Max'],
                        cellLoc='center',
                        loc='center',
                        colWidths=[0.15] * 7)
        
        # Style table
        table.auto_set_font_size(False)
        table.set_fontsize(10)
        table.scale(1, 2)
        
        # Color header
        for i in range(7):
            table[(0, i)].set_facecolor('#4C78A8')
            table[(0, i)].set_text_props(weight='bold', color='white')
        
        # Alternate row colors
        for i in range(1, len(table_data) + 1):
            if i % 2 == 0:
                for j in range(7):
                    table[(i, j)].set_facecolor('#f0f0f0')
        
        ax.set_title('Database Statistics Summary', fontsize=14, fontweight='bold', pad=20)
        plt.tight_layout()
        
        # Save figure
        self._save_figure(fig, filename)
        print(f"    ✓ Created {filename}.[png/pdf/svg]")
    
    def create_virulence_analysis(self, vfdb_data: Dict[str, pd.DataFrame]):
        """Create virulence gene analysis"""
        print("\n📊 Creating virulence gene analysis...")
        
        if not vfdb_data or 'gene_frequency' not in vfdb_data or vfdb_data['gene_frequency'].empty:
            print("  ⚠️ No VFDB data available")
            return
        
        gene_freq = vfdb_data['gene_frequency']
        
        # Identify critical virulence genes
        critical_virulence = []
        for gene in gene_freq['Gene']:
            if any(crit in str(gene).lower() for crit in self.critical_virulence_genes):
                critical_virulence.append(gene)
        
        # Filter for top genes (by count)
        top_genes = gene_freq.nlargest(20, 'Count') if len(gene_freq) > 20 else gene_freq
        
        # Plot top virulence genes
        fig, ax = plt.subplots(figsize=(16, 10))
        
        # Sort by count
        top_genes = top_genes.sort_values('Count', ascending=True)
        
        # Create horizontal bar chart
        y_pos = np.arange(len(top_genes))
        bars = ax.barh(y_pos, top_genes['Count'], 
                      color=['#DC143C' if any(crit in str(gene).lower() for crit in self.critical_virulence_genes) 
                            else '#4daf4a' for gene in top_genes['Gene']],
                      edgecolor='black', linewidth=0.5)
        
        ax.set_yticks(y_pos)
        ax.set_yticklabels(top_genes['Gene'], fontsize=10)
        ax.set_xlabel('Number of Genomes', fontsize=12, fontweight='bold')
        ax.set_title('Top Virulence Genes in S. aureus\n(VFDB Analysis)', 
                    fontsize=14, fontweight='bold', pad=20)
        
        # Add count labels
        for bar, count in zip(bars, top_genes['Count']):
            width = bar.get_width()
            ax.text(width + max(top_genes['Count']) * 0.01, bar.get_y() + bar.get_height()/2,
                   f'{count}', va='center', fontsize=9, fontweight='bold')
        
        # Add critical gene annotations
        for i, gene in enumerate(top_genes['Gene']):
            if any(crit in str(gene).lower() for crit in self.critical_virulence_genes):
                ax.text(0, i, ' ⚠️ CRITICAL', va='center', fontsize=9, 
                       fontweight='bold', color='white')
        
        ax.grid(True, alpha=0.3, linestyle='--', axis='x')
        plt.tight_layout()
        
        # Save figure
        self._save_figure(fig, "virulence_gene_analysis")
        print("    ✓ Created virulence_gene_analysis.[png/pdf/svg]")
        
        # Create critical virulence genes plot if found
        if critical_virulence:
            self._create_critical_virulence_plot(gene_freq, critical_virulence)
    
    def _create_critical_virulence_plot(self, gene_freq: pd.DataFrame, critical_genes: List[str]):
        """Create plot for critical virulence genes"""
        critical_data = gene_freq[gene_freq['Gene'].isin(critical_genes)]
        
        if not critical_data.empty:
            fig, ax = plt.subplots(figsize=(14, 8))
            
            # Sort by count
            critical_data = critical_data.sort_values('Count', ascending=True)
            
            # Create horizontal bar chart
            y_pos = np.arange(len(critical_data))
            bars = ax.barh(y_pos, critical_data['Count'], 
                          color='#DC143C', edgecolor='black', linewidth=1)
            
            ax.set_yticks(y_pos)
            ax.set_yticklabels(critical_data['Gene'], fontsize=11, fontweight='bold')
            ax.set_xlabel('Number of Genomes', fontsize=12, fontweight='bold')
            ax.set_title('Critical Virulence Genes in S. aureus\n(PVL, TSST-1, Enterotoxins)', 
                        fontsize=14, fontweight='bold', pad=20)
            
            # Add count labels
            for bar, count in zip(bars, critical_data['Count']):
                width = bar.get_width()
                ax.text(width + max(critical_data['Count']) * 0.01, bar.get_y() + bar.get_height()/2,
                       f'{count} genomes', va='center', fontsize=10, fontweight='bold')
            
            ax.grid(True, alpha=0.3, linestyle='--', axis='x')
            plt.tight_layout()
            
            # Save figure
            self._save_figure(fig, "critical_virulence_genes")
            print("    ✓ Created critical_virulence_genes.[png/pdf/svg]")


class StaphVisualizationReporter:
    """Main class to orchestrate the complete visualization pipeline for S. aureus"""
    
    def __init__(self, input_dir: Path = Path.cwd(), output_dir: Path = None):
        self.input_dir = Path(input_dir)
        
        if output_dir:
            self.output_dir = Path(output_dir)
        else:
            self.output_dir = self.input_dir / "STAPHSCOPE_VISUALIZATIONS"
        
        self.parser = StaphHTMLParser()
        self.visualizer = StaphVisualizer(self.output_dir)
        
        # Data storage
        self.typing_data = pd.DataFrame()
        self.amrfinder_data = {}
        self.databases_data = {}  # All ABRicate databases
        self.vfdb_data = {}  # Special handling for virulence
    
    def run_pipeline(self):
        """Run the complete visualization pipeline"""
        print("=" * 70)
        print("🧬 STAPHSCOPE VISUALIZATION MODULE - S. AUREUS EDITION")
        print("=" * 70)
        print(f"Input directory: {self.input_dir}")
        print(f"Output directory: {self.output_dir}")
        print("-" * 70)
        
        start_time = datetime.now()
        
        # Step 1: Parse all HTML files
        self._parse_all_files()
        
        # Step 2: Generate MRSA Analysis
        print("\n" + "=" * 70)
        print("🎯 CATEGORY 1: MRSA & TYPING ANALYSIS")
        print("=" * 70)
        self._generate_mrsa_typing_analysis()
        
        # Step 3: Generate AMR Analysis
        print("\n" + "=" * 70)
        print("🎯 CATEGORY 2: AMR GENE ANALYSIS")
        print("=" * 70)
        self._generate_amr_analysis()
        
        # Step 4: Generate Virulence Analysis
        print("\n" + "=" * 70)
        print("🎯 CATEGORY 3: VIRULENCE GENE ANALYSIS")
        print("=" * 70)
        self._generate_virulence_analysis()
        
        # Step 5: Generate summary report
        self._generate_summary_report(start_time)
        
        print("\n" + "=" * 70)
        print("✅ STAPHSCOPE VISUALIZATION PIPELINE COMPLETED!")
        print("=" * 70)
    
    def _parse_all_files(self):
        """Parse all StaphScope HTML files in the input directory"""
        print("\n📁 Parsing StaphScope HTML files...")
        
        # Parse comprehensive typing report
        comp_files = list(self.input_dir.glob("*comprehensive*.html"))
        if comp_files:
            self.typing_data = self.parser.parse_comprehensive_html(comp_files[0])
        
        # Parse AMRfinder report
        amrfinder_files = list(self.input_dir.glob("*amrfinder*.html"))
        if amrfinder_files:
            self.amrfinder_data = self.parser.parse_amrfinder_html(amrfinder_files[0])
        
        # Parse all ABRicate database files
        database_patterns = [
            "*vfdb*.html",
            "*card*.html",
            "*resfinder*.html",
            "*plasmidfinder*.html",
            "*argannot*.html",
            "*megares*.html",
            "*ncbi*.html"
        ]
        
        parsed_files = set()
        
        for pattern in database_patterns:
            for db_file in self.input_dir.glob(pattern):
                # Skip files already parsed
                if (db_file in comp_files or db_file in amrfinder_files or 
                    str(db_file).lower() in parsed_files):
                    continue
                    
                db_data = self.parser.parse_abricate_html(db_file)
                db_name = db_data.get('database', 'unknown')
                
                # Only add if we have valid data
                if (not db_data['genes_by_genome'].empty or 
                    not db_data['gene_frequency'].empty):
                    self.databases_data[db_name] = db_data
                    parsed_files.add(str(db_file).lower())
                    
                    # Store VFDB separately for virulence analysis
                    if db_name == 'vfdb':
                        self.vfdb_data = db_data
    
    def _generate_mrsa_typing_analysis(self):
        """Generate MRSA and typing analysis visualizations"""
        # MRSA Analysis
        if not self.typing_data.empty:
            self.visualizer.create_mrsa_analysis(self.typing_data)
            
            # Typing Distributions
            self.visualizer.create_typing_distributions(self.typing_data)
            
            # Stacked Combinations
            self.visualizer.create_stacked_combinations(self.typing_data)
        else:
            print("  ⚠️ No typing data available for MRSA analysis")
    
    def _generate_amr_analysis(self):
        """Generate AMR gene analysis visualizations"""
        if self.amrfinder_data or self.databases_data:
            self.visualizer.create_amr_analysis(self.amrfinder_data, self.databases_data)
        else:
            print("  ⚠️ No AMR data available for analysis")
    
    def _generate_virulence_analysis(self):
        """Generate virulence gene analysis visualizations"""
        if self.vfdb_data:
            self.visualizer.create_virulence_analysis(self.vfdb_data)
        else:
            print("  ⚠️ No VFDB data available for virulence analysis")
    
    def _generate_summary_report(self, start_time):
        """Generate a comprehensive summary report"""
        end_time = datetime.now()
        duration = end_time - start_time
        
        report_path = self.output_dir / "staphscope_visualization_report.txt"
        
        with open(report_path, 'w') as f:
            f.write("=" * 70 + "\n")
            f.write("🧬 STAPHSCOPE VISUALIZATION MODULE - SUMMARY REPORT\n")
            f.write("=" * 70 + "\n\n")
            
            f.write("📋 EXECUTION DETAILS\n")
            f.write("-" * 40 + "\n")
            f.write(f"Start time: {start_time.strftime('%Y-%m-%d %H:%M:%S')}\n")
            f.write(f"End time: {end_time.strftime('%Y-%m-%d %H:%M:%S')}\n")
            f.write(f"Duration: {duration}\n")
            f.write(f"Input directory: {self.input_dir}\n")
            f.write(f"Output directory: {self.output_dir}\n\n")
            
            f.write("📊 DATA SUMMARY\n")
            f.write("-" * 40 + "\n")
            
            # Typing data summary
            if not self.typing_data.empty:
                total_samples = len(self.typing_data)
                mrsa_count = len(self.typing_data[self.typing_data['MRSA_Status'] == 'MRSA']) if 'MRSA_Status' in self.typing_data.columns else 0
                unique_mlst = self.typing_data['MLST'].nunique() if 'MLST' in self.typing_data.columns else 0
                unique_spa = self.typing_data['spa_Type'].nunique() if 'spa_Type' in self.typing_data.columns else 0
                unique_sccmec = self.typing_data['SCCmec_Type'].nunique() if 'SCCmec_Type' in self.typing_data.columns else 0
                
                f.write(f"Typing Data: {total_samples} samples\n")
                f.write(f"  • MRSA: {mrsa_count} ({mrsa_count/total_samples*100:.1f}%)\n")
                f.write(f"  • Unique MLST: {unique_mlst}\n")
                f.write(f"  • Unique spa types: {unique_spa}\n")
                f.write(f"  • Unique SCCmec types: {unique_sccmec}\n")
            else:
                f.write("Typing Data: No data found\n")
            
            # AMRfinder data summary
            if self.amrfinder_data and 'gene_frequency' in self.amrfinder_data:
                amr_genes = len(self.amrfinder_data['gene_frequency'])
                f.write(f"AMRfinder Data: {amr_genes} AMR genes\n")
            
            # Database summary
            f.write(f"Databases Analyzed: {len(self.databases_data)}\n")
            for db_name in sorted(self.databases_data.keys()):
                f.write(f"  • {db_name.upper()}\n")
            
            f.write("\n📈 GENERATED VISUALIZATIONS\n")
            f.write("-" * 40 + "\n")
            
            # List generated files
            for fmt in ['PNG', 'PDF', 'SVG', 'DATA']:
                fmt_dir = self.output_dir / fmt
                if fmt_dir.exists():
                    files = list(fmt_dir.glob("*"))
                    if files:
                        f.write(f"\n{fmt} Files ({len(files)}):\n")
                        for file in sorted(files):
                            f.write(f"  • {file.name}\n")
            
            f.write("\n" + "=" * 70 + "\n")
            f.write("✅ REPORT COMPLETE\n")
            f.write("=" * 70 + "\n")
        
        print(f"\n📋 Summary report saved: {report_path}")
        
        # Print summary to console
        print("\n📋 FINAL SUMMARY:")
        print("-" * 40)
        print(f"Total time: {duration}")
        print(f"Total plots generated: {len(list((self.output_dir / 'PNG').glob('*.png')))}")
        print(f"Output directory: {self.output_dir}")
        print(f"\n🎯 Key visualizations created:")
        print("  • MRSA vs MSSA distribution charts")
        print("  • SCCmec type analysis")
        print("  • MLST and spa typing distributions")
        print("  • Critical AMR gene analysis (mecA/mecC)")
        print("  • Database comparison plots")
        print("  • Virulence gene analysis")
        print("  • Statistical summary tables")


def main():
    """Main function to run the StaphScope visualization pipeline"""
    parser = argparse.ArgumentParser(
        description="STAPHSCOPE Visualization Module - Parse HTML reports and generate publication-quality visualizations for S. aureus",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python staphscope_visualizer.py
  python staphscope_visualizer.py --input /path/to/html/files --output /path/to/results
  python staphscope_visualizer.py --force
        
This script will:
  1. Parse all StaphScope HTML files in the input directory
  2. Generate MRSA vs MSSA analysis with pie/donut charts
  3. Create SCCmec type distribution plots
  4. Generate MLST and spa typing visualizations
  5. Create critical AMR gene analysis (mecA/mecC)
  6. Generate database comparison plots
  7. Create virulence gene analysis
  8. Save results in PNG, PDF, SVG, and CSV formats
        """
    )
    
    parser.add_argument("--input", type=str, default=".",
                       help="Directory containing StaphScope HTML files (default: current directory)")
    parser.add_argument("--output", type=str, default=None,
                       help="Output directory (default: STAPHSCOPE_VISUALIZATIONS in input directory)")
    parser.add_argument("--force", action="store_true",
                       help="Overwrite existing output directory")
    
    args = parser.parse_args()
    
    # Create reporter
    reporter = StaphVisualizationReporter(
        input_dir=Path(args.input),
        output_dir=Path(args.output) if args.output else None
    )
    
    # Run pipeline
    try:
        reporter.run_pipeline()
    except Exception as e:
        print(f"\n❌ ERROR: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main()

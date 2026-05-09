#!/usr/bin/env python3
"""
STAPHSCOPE ULTIMATE REPORTER - GENE-CENTRIC S. AUREUS ANALYSIS
Enhanced with FASTA QC, AI Guide, Filter Buttons, Combination Tables, Sortable Headers
Author: Beckley Brown <brownbeckley94@gmail.com>
Affiliation: University of Ghana Medical School
Version: 2.1.0 - Ultimate Edition with Scrollable Genome Lists, Triple Typing, BACMET, Citations and Call to Action!
Date: 2026-05-08
"""

import os
import sys
import json
import re
import glob
import argparse
import pandas as pd
import numpy as np
from pathlib import Path
from typing import Dict, List, Set, Tuple, Any, Optional
from datetime import datetime
from collections import defaultdict, Counter
import warnings
warnings.filterwarnings('ignore')

# HTML parsing
from bs4 import BeautifulSoup

# Visualization
try:
    import plotly.graph_objects as go
    import plotly.express as px
    PLOTLY_AVAILABLE = True
except ImportError:
    PLOTLY_AVAILABLE = False

class StaphHTMLParser:
    """Ultimate HTML parser for all StaphScope reports"""
    
    def __init__(self):
        self.abricate_databases = [
            'card', 'resfinder', 'vfdb', 'argannot',
            'plasmidfinder', 'megares', 'ncbi', 'bacmet2'
        ]
    
    def normalize_sample_id(self, sample_id: str) -> str:
        """Normalize sample ID"""
        sample = str(sample_id)
        extensions = ['.fna', '.fasta', '.fa', '.gb', '.gbk', '.gbff', '.txt', '.tsv', '.csv']
        for ext in extensions:
            if sample.endswith(ext):
                sample = sample[:-len(ext)]
        
        if '/' in sample or '\\' in sample:
            sample = Path(sample).name
        
        return sample.strip()
    
    def parse_html_table(self, html_content: str, table_index: int = 0) -> pd.DataFrame:
        """Parse HTML table"""
        try:
            soup = BeautifulSoup(html_content, 'html.parser')
            tables = soup.find_all('table')
            
            if not tables or table_index >= len(tables):
                return pd.DataFrame()
            
            table = tables[table_index]
            rows = table.find_all('tr')
            
            headers = []
            for th in rows[0].find_all(['th', 'td']):
                headers.append(th.get_text().strip())
            
            data = []
            for row in rows[1:]:
                cols = row.find_all(['td', 'th'])
                if cols:
                    row_data = [col.get_text().strip() for col in cols]
                    if len(row_data) == len(headers):
                        data.append(row_data)
            
            if not data:
                return pd.DataFrame()
            
            return pd.DataFrame(data, columns=headers)
            
        except Exception as e:
            print(f"  ⚠️ Table parsing error: {e}")
            return pd.DataFrame()
    
    def parse_qc_report(self, file_path: Path) -> Dict[str, Dict]:
        """Parse FASTA QC summary HTML – extract all available metrics."""
        print(f"  🧬 Parsing FASTA QC: {file_path.name}")
        try:
            with open(file_path, 'r', encoding='utf-8') as f:
                html = f.read()
            df = self.parse_html_table(html, 0)
            if df.empty:
                return {}

            # Find sample column (usually 'Filename')
            sample_col = None
            for col in df.columns:
                if 'filename' in col.lower() or 'sample' in col.lower() or col == df.columns[0]:
                    sample_col = col
                    break
            if not sample_col:
                return {}

            results = {}
            for _, row in df.iterrows():
                sample_raw = row[sample_col]
                if not sample_raw:
                    continue
                sample = self.normalize_sample_id(sample_raw)
                # Extract all available columns
                qc_data = {}
                for col in df.columns:
                    if col == sample_col:
                        continue
                    val = row[col]
                    if pd.isna(val) or val == '' or val == 'ND':
                        qc_data[col] = 'ND'
                    else:
                        cleaned = str(val).replace('%', '').replace(',', '').strip()
                        try:
                            qc_data[col] = float(cleaned)
                        except:
                            qc_data[col] = str(val)
                results[sample] = qc_data
            print(f"    ✓ Parsed {len(results)} samples")
            return results
        except Exception as e:
            print(f"    ❌ Error parsing QC: {e}")
            return {}
    
    def parse_comprehensive_report(self, file_path: Path) -> Dict[str, Dict]:
        """Parse StaphScope comprehensive typing report"""
        print(f"  🧬 Parsing Comprehensive Report: {file_path.name}")
        
        try:
            with open(file_path, 'r', encoding='utf-8') as f:
                html_content = f.read()
            
            soup = BeautifulSoup(html_content, 'html.parser')
            tables = soup.find_all('table')
            
            typing_table = None
            for table in tables:
                if table.find(string=re.compile(r'Sample|MLST|spa|SCCmec|MRSA', re.I)):
                    typing_table = table
                    break
            
            if not typing_table:
                typing_table = soup.find('table')
            
            if not typing_table:
                return {}
            
            data = []
            rows = typing_table.find_all('tr')
            
            if len(rows) < 2:
                return {}
            
            headers = []
            header_cells = rows[0].find_all(['th', 'td'])
            for cell in header_cells:
                headers.append(cell.get_text().strip())
            
            for row in rows[1:]:
                cols = row.find_all(['td', 'th'])
                if cols:
                    row_data = [col.get_text().strip() for col in cols]
                    if len(row_data) >= 2:
                        data.append(row_data)
            
            if not data:
                return {}
            
            df = pd.DataFrame(data)
            if len(df.columns) > len(headers):
                df = df.iloc[:, :len(headers)]
            df.columns = headers[:len(df.columns)]
            
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
            
            df.rename(columns=column_mapping, inplace=True)
            
            # Ensure required columns exist
            if 'Sample' not in df.columns and len(df.columns) > 0:
                df.rename(columns={df.columns[0]: 'Sample'}, inplace=True)
            
            df['normalized_sample'] = df['Sample'].apply(self.normalize_sample_id)
            
            results = {}
            for _, row in df.iterrows():
                sample = row['normalized_sample']
                
                mlst = row.get('MLST', 'ND') if 'MLST' in df.columns else 'ND'
                spa_type = row.get('spa_Type', 'ND') if 'spa_Type' in df.columns else 'ND'
                sccmec_type = row.get('SCCmec_Type', 'ND') if 'SCCmec_Type' in df.columns else 'ND'
                mrsa_status = row.get('MRSA_Status', 'ND') if 'MRSA_Status' in df.columns else 'ND'
                
                # Clean up values
                mlst = str(mlst).strip() if pd.notna(mlst) else 'ND'
                spa_type = str(spa_type).strip() if pd.notna(spa_type) else 'ND'
                sccmec_type = str(sccmec_type).strip() if pd.notna(sccmec_type) else 'ND'
                mrsa_status = str(mrsa_status).strip() if pd.notna(mrsa_status) else 'ND'
                
                results[sample] = {
                    'MLST': mlst,
                    'spa_Type': spa_type,
                    'SCCmec_Type': sccmec_type,
                    'MRSA_Status': mrsa_status
                }
            
            print(f"    ✓ Found {len(results)} samples")
            return results
            
        except Exception as e:
            print(f"    ❌ Error parsing comprehensive report: {e}")
            return {}
    
    def parse_amrfinder_report(self, file_path: Path) -> Tuple[Dict[str, List], Dict[str, Dict]]:
        """Parse AMRfinder HTML report for S. aureus"""
        print(f"  🧬 Parsing AMRfinder: {file_path.name}")
        
        try:
            with open(file_path, 'r', encoding='utf-8') as f:
                html_content = f.read()
            
            soup = BeautifulSoup(html_content, 'html.parser')
            tables = soup.find_all('table')
            
            print(f"    Found {len(tables)} tables in the HTML")
            
            genes_by_genome = {}
            gene_frequencies = {}
            
            # Look for the "Genes by Genome" table
            genes_by_genome_table = None
            gene_frequency_table = None
            
            # Find tables by their content
            for i, table in enumerate(tables):
                table_text = table.get_text()
                # Check if this is the "Genes by Genome" table
                if 'Genome' in table_text and 'Critical Genes' in table_text:
                    genes_by_genome_table = table
                    print(f"    Found 'Genes by Genome' table at index {i}")
                # Check if this is the "Gene Frequency" table
                elif 'Gene' in table_text and 'Frequency' in table_text and 'Prevalence' in table_text:
                    gene_frequency_table = table
                    print(f"    Found 'Gene Frequency' table at index {i}")
            
            # Parse "Genes by Genome" table
            if genes_by_genome_table:
                print(f"    Parsing 'Genes by Genome' table...")
                try:
                    table_html = str(genes_by_genome_table)
                    df_genomes = pd.read_html(table_html)[0]
                    
                    print(f"    Table columns: {list(df_genomes.columns)}")
                    
                    genome_col = None
                    for col in df_genomes.columns:
                        if 'genome' in col.lower() or 'Genome' in col:
                            genome_col = col
                            break
                    
                    if not genome_col and len(df_genomes.columns) > 0:
                        genome_col = df_genomes.columns[0]
                    
                    if genome_col:
                        for _, row in df_genomes.iterrows():
                            sample = self.normalize_sample_id(row[genome_col])
                            
                            critical_genes = []
                            crit_col = None
                            for col in df_genomes.columns:
                                if 'critical' in col.lower():
                                    crit_col = col
                                    break
                            
                            if crit_col and pd.notna(row.get(crit_col)):
                                crit_val = str(row[crit_col]).strip()
                                if crit_val.lower() not in ['none', 'nan', '']:
                                    critical_genes = [g.strip() for g in crit_val.split(',') if g.strip()]
                            
                            high_risk_genes = []
                            high_risk_col = None
                            for col in df_genomes.columns:
                                if 'high risk' in col.lower() or 'high' in col.lower():
                                    high_risk_col = col
                                    break
                            
                            if high_risk_col and pd.notna(row.get(high_risk_col)):
                                high_val = str(row[high_risk_col]).strip()
                                if high_val.lower() not in ['none', 'nan', '']:
                                    high_risk_genes = [g.strip() for g in high_val.split(',') if g.strip()]
                            
                            genes_by_genome[sample] = {
                                'critical_genes': critical_genes,
                                'high_risk_genes': high_risk_genes,
                                'all_genes': []
                            }
                        
                        print(f"    Parsed {len(genes_by_genome)} samples from 'Genes by Genome' table")
                except Exception as e:
                    print(f"    Error parsing 'Genes by Genome' table with pandas: {e}")
            
            # Parse "Gene Frequency" table
            if gene_frequency_table:
                print(f"    Parsing 'Gene Frequency' table...")
                try:
                    table_html = str(gene_frequency_table)
                    df_genes = pd.read_html(table_html)[0]
                    
                    print(f"    Gene table columns: {list(df_genes.columns)}")
                    
                    gene_col = None
                    for col in df_genes.columns:
                        if 'gene' in col.lower() or 'Gene' in col:
                            gene_col = col
                            break
                    
                    if gene_col:
                        for _, row in df_genes.iterrows():
                            if pd.isna(row[gene_col]):
                                continue
                            
                            gene = str(row[gene_col]).strip()
                            
                            frequency = '0'
                            freq_col = None
                            for col in df_genes.columns:
                                if 'freq' in col.lower():
                                    freq_col = col
                                    break
                            
                            if freq_col and pd.notna(row.get(freq_col)):
                                frequency = str(row[freq_col]).strip()
                            
                            count = 0
                            match = re.search(r'(\d+)', frequency)
                            if match:
                                count = int(match.group(1))
                            
                            prevalence = 'ND'
                            prev_col = None
                            for col in df_genes.columns:
                                if 'prev' in col.lower():
                                    prev_col = col
                                    break
                            
                            if prev_col and pd.notna(row.get(prev_col)):
                                prevalence = str(row[prev_col]).strip()
                            
                            risk_level = 'STANDARD'
                            risk_col = None
                            for col in df_genes.columns:
                                if 'risk' in col.lower():
                                    risk_col = col
                                    break
                            
                            if risk_col and pd.notna(row.get(risk_col)):
                                risk_level = str(row[risk_col]).strip()
                            
                            genomes = []
                            genomes_col = None
                            for col in df_genes.columns:
                                if 'genome' in col.lower() or 'Genomes' in col:
                                    genomes_col = col
                                    break
                            
                            if genomes_col and pd.notna(row.get(genomes_col)):
                                genomes_str = str(row[genomes_col])
                                genomes_str = genomes_str.replace('tet(38)', 'tet38')
                                genomes_list = [g.strip() for g in genomes_str.split(',') if g.strip()]
                                genomes = [self.normalize_sample_id(g.replace('tet38', 'tet(38)')) for g in genomes_list]
                            
                            gene_frequencies[gene] = {
                                'frequency': frequency,
                                'count': count,
                                'prevalence': prevalence,
                                'risk_level': risk_level,
                                'genomes': genomes,
                                'database': 'amrfinder'
                            }
                            
                            for genome in genomes:
                                if genome in genes_by_genome:
                                    if gene not in genes_by_genome[genome]['all_genes']:
                                        genes_by_genome[genome]['all_genes'].append(gene)
                                else:
                                    genes_by_genome[genome] = {
                                        'critical_genes': [],
                                        'high_risk_genes': [],
                                        'all_genes': [gene]
                                    }
                        
                        print(f"    Parsed {len(gene_frequencies)} genes from 'Gene Frequency' table")
                except Exception as e:
                    print(f"    Error parsing 'Gene Frequency' table with pandas: {e}")
                    import traceback
                    traceback.print_exc()
            
            # If pandas failed, try alternative parsing
            if not gene_frequencies and gene_frequency_table:
                print(f"    Trying alternative parsing for gene frequency table...")
                rows = gene_frequency_table.find_all('tr')
                
                if len(rows) > 1:
                    headers = []
                    header_cells = rows[0].find_all(['th', 'td'])
                    for cell in header_cells:
                        headers.append(cell.get_text().strip())
                    
                    print(f"    Headers found: {headers}")
                    
                    gene_idx = -1
                    freq_idx = -1
                    prev_idx = -1
                    risk_idx = -1
                    genomes_idx = -1
                    
                    for i, header in enumerate(headers):
                        if 'gene' in header.lower():
                            gene_idx = i
                        elif 'freq' in header.lower():
                            freq_idx = i
                        elif 'prev' in header.lower():
                            prev_idx = i
                        elif 'risk' in header.lower():
                            risk_idx = i
                        elif 'genome' in header.lower():
                            genomes_idx = i
                    
                    for row in rows[1:]:
                        cols = row.find_all(['td', 'th'])
                        if cols and len(cols) >= len(headers):
                            row_data = [col.get_text().strip() for col in cols]
                            
                            if gene_idx >= 0 and gene_idx < len(row_data):
                                gene = row_data[gene_idx]
                                
                                frequency = '0'
                                if freq_idx >= 0 and freq_idx < len(row_data):
                                    frequency = row_data[freq_idx]
                                
                                count = 0
                                match = re.search(r'(\d+)', frequency)
                                if match:
                                    count = int(match.group(1))
                                
                                prevalence = 'ND'
                                if prev_idx >= 0 and prev_idx < len(row_data):
                                    prevalence = row_data[prev_idx]
                                
                                risk_level = 'STANDARD'
                                if risk_idx >= 0 and risk_idx < len(row_data):
                                    risk_level = row_data[risk_idx]
                                
                                genomes = []
                                if genomes_idx >= 0 and genomes_idx < len(row_data):
                                    genomes_str = row_data[genomes_idx]
                                    genomes_list = [g.strip() for g in genomes_str.split(',') if g.strip()]
                                    genomes = [self.normalize_sample_id(g) for g in genomes_list]
                                
                                gene_frequencies[gene] = {
                                    'frequency': frequency,
                                    'count': count,
                                    'prevalence': prevalence,
                                    'risk_level': risk_level,
                                    'genomes': genomes,
                                    'database': 'amrfinder'
                                }
                    
                    print(f"    Alternative parsing found {len(gene_frequencies)} genes")
            
            print(f"    ✓ Found {len(genes_by_genome)} samples, {len(gene_frequencies)} genes")
            return genes_by_genome, gene_frequencies
            
        except Exception as e:
            print(f"    ❌ Error parsing AMRfinder: {e}")
            import traceback
            traceback.print_exc()
            return {}, {}    
    
    def parse_abricate_report(self, file_path: Path) -> Tuple[str, Dict[str, List], Dict[str, Dict]]:
        """Parse ANY ABRicate database HTML report"""
        print(f"  🧬 Parsing ABRicate: {file_path.name}")
        
        try:
            with open(file_path, 'r', encoding='utf-8') as f:
                html_content = f.read()
            
            soup = BeautifulSoup(html_content, 'html.parser')
            tables = soup.find_all('table')
            
            if len(tables) < 2:
                return 'unknown', {}, {}
            
            # Determine database name from filename
            db_name = 'unknown'
            filename = file_path.name.lower()
            
            for db in self.abricate_databases:
                if db in filename:
                    db_name = db
                    break
            
            # Also try to get from title
            title_tag = soup.find('title')
            if title_tag:
                title_text = title_tag.get_text().lower()
                for db in self.abricate_databases:
                    if db in title_text:
                        db_name = db
                        break
            
            # Parse table 1: Genes by Genome
            genes_by_genome = {}
            df1 = self.parse_html_table(html_content, 0)
            if not df1.empty:
                df1.columns = [col.strip() for col in df1.columns]
                
                # Find genome column
                genome_col = None
                for col in df1.columns:
                    if 'genome' in col.lower() or 'sample' in col.lower():
                        genome_col = col
                        break
                
                if not genome_col and len(df1.columns) > 0:
                    genome_col = df1.columns[0]
                
                if genome_col:
                    for _, row in df1.iterrows():
                        sample = self.normalize_sample_id(row[genome_col])
                        genes = []
                        
                        # Find genes column
                        genes_col = None
                        for col in df1.columns:
                            if 'genes' in col.lower() or 'detected' in col.lower():
                                genes_col = col
                                break
                        
                        if genes_col and pd.notna(row.get(genes_col)):
                            gene_str = str(row[genes_col])
                            genes = [g.strip() for g in gene_str.split(',') if g.strip()]
                        
                        genes_by_genome[sample] = genes
            
            # Parse table 2: Gene Frequency
            gene_frequencies = {}
            df2 = self.parse_html_table(html_content, 1)
            if not df2.empty:
                df2.columns = [col.strip() for col in df2.columns]
                
                if 'Gene' in df2.columns:
                    for _, row in df2.iterrows():
                        gene = str(row['Gene']).strip()
                        frequency = str(row.get('Frequency', '0')).strip()
                        
                        genomes = []
                        if 'Genomes' in df2.columns and pd.notna(row.get('Genomes')):
                            genomes_str = str(row['Genomes'])
                            genomes = [self.normalize_sample_id(g.strip()) 
                                      for g in genomes_str.split(',') if g.strip()]
                        
                        count = 0
                        match = re.search(r'(\d+)', frequency)
                        if match:
                            count = int(match.group(1))
                        
                        gene_frequencies[gene] = {
                            'frequency': frequency,
                            'count': count,
                            'genomes': genomes,
                            'database': db_name
                        }
            
            print(f"    ✓ {db_name.upper()}: {len(genes_by_genome)} samples, {len(gene_frequencies)} genes")
            return db_name, genes_by_genome, gene_frequencies
            
        except Exception as e:
            print(f"    ❌ Error parsing ABRicate report: {e}")
            return 'unknown', {}, {}


class StaphDataAnalyzer:
    """Analyzes data for S. aureus gene-centric reporting"""
    
    def __init__(self):
        # S. aureus critical AMR genes
        self.critical_amr_genes = {
            'meca', 'mecA', 'mecc', 'mecC', 'vana', 'vanA', 'vanb', 'vanB',
            'vanc', 'vanC', 'erma', 'ermA', 'ermb', 'ermB', 'ermc', 'ermC',
            'msra', 'msrA', 'mphc', 'mphC', 'tetk', 'tetK', 'tetm', 'tetM', 'tetl', 'tetL'
        }
        
        # High-priority AMR genes for filter buttons 
        self.high_priority_amr = [
            'mecA', 'mecC', 'vanA', 'vanB', 'ermA', 'ermB', 'ermC',
            'msrA', 'mphC', 'tetK', 'tetM', 'aacA-aphD', 'aac(6\')-aph(2\'\')',
            'ant(4\')-Ia', 'ant(6)-Ia', 'aph(3\')-IIIa', 'satA', 'dfrA', 'dfrG', 'cat'
        ]
        
        # S. aureus critical virulence genes
        self.critical_virulence_genes = {
            'luks-pv', 'lukS-PV', 'lukf-pv', 'lukF-PV',  # PVL toxin
            'tsst',  # TSST-1
            'sea', 'seb', 'sec', 'sed', 'see',  # Enterotoxins A-E
            'seg', 'seh', 'sei', 'sej', 'sek', 'sel', 'sem', 'sen', 'seo', 'sep', 'seq', 'ser', 'seu',  # Other enterotoxins
            'eta', 'etb',  # Exfoliative toxins
            'hla', 'hlb', 'hlg', 'hld',  # Hemolysins
        }
        
        # High-priority virulence genes for filter buttons 
        self.high_priority_virulence = [
            'lukF-PV', 'lukS-PV', 'tsst', 'sea', 'seb', 'sec', 'sed', 'see',
            'seg', 'seh', 'sei', 'sej', 'sek', 'sel', 'sem', 'sen', 'seo', 'sep',
            'eta', 'etb', 'hla', 'hlb', 'hlg', 'hld'
        ]
        
        # SCCmec types for MRSA classification
        self.sccmec_types = {
            'I', 'II', 'III', 'IV', 'V', 'VI', 'VII', 'VIII', 'IX', 'X', 'XI', 'XII'
        }
    
    def create_gene_centric_tables(self, integrated_data: Dict[str, Any]) -> Dict[str, Any]:
        """Create gene-centric tables showing genes with their genomes"""
        gene_centric = {
            'amr_databases': {},
            'virulence_databases': {},
            'plasmid_databases': {},
            'bacmet_databases': {},  # New for biocide/heavy metal
            'combined_gene_frequencies': []
        }
        
        # Process AMRfinder
        if 'amrfinder' in integrated_data.get('gene_frequencies', {}):
            amr_data = integrated_data['gene_frequencies']['amrfinder']
            gene_list = []
            
            for gene, data in amr_data.items():
                gene_list.append({
                    'gene': gene,
                    'database': 'AMRfinder',
                    'frequency': data.get('frequency', '0'),
                    'count': data.get('count', 0),
                    'prevalence': data.get('prevalence', 'ND'),
                    'risk_level': data.get('risk_level', 'ND'),
                    'genomes': data.get('genomes', [])
                })
            
            gene_centric['amr_databases']['amrfinder'] = sorted(gene_list, key=lambda x: x['count'], reverse=True)
        
        # Process ABRicate databases
        if 'abricate' in integrated_data.get('gene_frequencies', {}):
            abricate_data = integrated_data['gene_frequencies']['abricate']
            
            for db_name, db_genes in abricate_data.items():
                gene_list = []
                
                for gene, data in db_genes.items():
                    gene_list.append({
                        'gene': gene,
                        'database': db_name.upper(),
                        'frequency': data.get('frequency', '0'),
                        'count': data.get('count', 0),
                        'genomes': data.get('genomes', [])
                    })
                
                # Sort by count and store in appropriate category
                if gene_list:
                    gene_list.sort(key=lambda x: x['count'], reverse=True)
                    
                    if db_name == 'vfdb':
                        gene_centric['virulence_databases'][db_name] = gene_list
                    elif db_name == 'plasmidfinder':
                        gene_centric['plasmid_databases'][db_name] = gene_list
                    elif db_name == 'bacmet2':
                        gene_centric['bacmet_databases'][db_name] = gene_list
                    else:
                        gene_centric['amr_databases'][db_name] = gene_list
        
        # Create combined gene frequencies for pattern discovery
        all_genes = []
        
        for db_type in ['amr_databases', 'virulence_databases', 'plasmid_databases', 'bacmet_databases']:
            for db_name, genes in gene_centric.get(db_type, {}).items():
                for gene_data in genes:
                    all_genes.append(gene_data)
        
        # Sort combined list by count
        all_genes.sort(key=lambda x: x['count'], reverse=True)
        gene_centric['combined_gene_frequencies'] = all_genes
        
        return gene_centric
    
    def create_cross_genome_patterns(self, integrated_data: Dict[str, Any]) -> Dict[str, Any]:
        """Create cross-genome patterns for S. aureus, including combination tables and triple typing"""
        patterns = {
            'mlst_distribution': Counter(),
            'spa_type_distribution': Counter(),
            'sccmec_distribution': Counter(),
            'mrsa_status_distribution': Counter(),
            'mlst_spa_combinations': defaultdict(list),
            'mlst_sccmec_combinations': defaultdict(list),
            'spa_sccmec_combinations': defaultdict(list),
            'triple_combinations': defaultdict(list),   # NEW: ST-spa-SCCmec
            'gene_cooccurrence': defaultdict(Counter),
            'high_risk_combinations': []
        }
        
        samples_data = integrated_data.get('samples', {})
        gene_centric = integrated_data.get('gene_centric', {})
        
        # Collect all genes per sample for co-occurrence
        sample_genes = defaultdict(list)
        for db_type in ['amr_databases', 'virulence_databases']:
            for db_name, genes in gene_centric.get(db_type, {}).items():
                for gene_data in genes:
                    for genome in gene_data['genomes']:
                        if gene_data['gene'] not in sample_genes[genome]:
                            sample_genes[genome].append(gene_data['gene'])
        
        # Analyze each sample
        for sample, data in samples_data.items():
            mlst = data.get('typing', {}).get('MLST', 'ND')
            spa_type = data.get('typing', {}).get('spa_Type', 'ND')
            sccmec_type = data.get('typing', {}).get('SCCmec_Type', 'ND')
            mrsa_status = data.get('typing', {}).get('MRSA_Status', 'ND')
            
            # Basic distributions
            if mlst != 'ND':
                patterns['mlst_distribution'][mlst] += 1
            if spa_type != 'ND':
                patterns['spa_type_distribution'][spa_type] += 1
            if sccmec_type != 'ND' and sccmec_type != 'Not Assigned':
                patterns['sccmec_distribution'][sccmec_type] += 1
            if mrsa_status != 'ND':
                patterns['mrsa_status_distribution'][mrsa_status] += 1
            
            # MLST-spa combinations
            if mlst != 'ND' and spa_type != 'ND':
                patterns['mlst_spa_combinations'][f"{mlst} - {spa_type}"].append(sample)
            
            # MLST-SCCmec combinations
            if mlst != 'ND' and sccmec_type != 'ND' and sccmec_type != 'Not Assigned':
                patterns['mlst_sccmec_combinations'][f"{mlst} - {sccmec_type}"].append(sample)
            
            # spa-SCCmec combinations
            if spa_type != 'ND' and sccmec_type != 'ND' and sccmec_type != 'Not Assigned':
                patterns['spa_sccmec_combinations'][f"{spa_type} - {sccmec_type}"].append(sample)
            
            # TRIPLE combination (MLST - spa - SCCmec)
            if mlst != 'ND' and spa_type != 'ND' and sccmec_type != 'ND' and sccmec_type != 'Not Assigned':
                patterns['triple_combinations'][f"{mlst} - {spa_type} - {sccmec_type}"].append(sample)
            
            # Gene co-occurrence
            genes = sample_genes.get(sample, [])
            for i, gene1 in enumerate(genes):
                for gene2 in genes[i+1:]:
                    patterns['gene_cooccurrence'][gene1][gene2] += 1
            
            # High-risk combinations (critical AMR + critical virulence)
            amr_genes = data.get('amrfinder', {}).get('all_genes', [])
            virulence_genes = data.get('abricate_databases', {}).get('vfdb', [])
            
            # Also include virulence factors from AMRfinder
            if 'amrfinder' in integrated_data.get('gene_frequencies', {}):
                amrfinder_genes = integrated_data['gene_frequencies']['amrfinder']
                for gene in amrfinder_genes:
                    if any(vir_gene in gene.lower() for vir_gene in self.critical_virulence_genes):
                        if gene not in virulence_genes:
                            virulence_genes.append(gene)
            
            critical_amr = [g for g in amr_genes if any(crit in str(g).lower() for crit in self.critical_amr_genes)]
            critical_vir = [g for g in virulence_genes if any(crit in str(g).lower() for crit in self.critical_virulence_genes)]
            
            if critical_amr and critical_vir:
                patterns['high_risk_combinations'].append({
                    'sample': sample,
                    'mlst': mlst,
                    'spa_type': spa_type,
                    'sccmec_type': sccmec_type,
                    'mrsa_status': mrsa_status,
                    'critical_amr_genes': critical_amr,
                    'critical_virulence_genes': critical_vir
                })
        
        return patterns


class StaphHTMLGenerator:
    """Generates ultimate HTML reports for S. aureus with gene-centric approach"""
    
    def __init__(self, data_analyzer: StaphDataAnalyzer):
        self.data_analyzer = data_analyzer
        self.tab_colors = {
            'summary': '#4CAF50',
            'sample_overview': '#2196F3',
            'qc': '#607D8B',
            'mlst': '#FF9800',
            'spa': '#9C27B0',
            'sccmec': '#009688',
            'mrsa': '#795548',
            'amr': '#F44336',
            'virulence': '#E91E63',
            'bacmet': '#FF5722',
            'plasmids': '#673AB7',
            'patterns': '#3F51B5',
            'aiguide': '#00BCD4',
            'citation': '#8BC34A',
            'calltoaction': '#FFC107',
            'export': '#9E9E9E'
        }
    
    def generate_main_report(self, integrated_data: Dict[str, Any], output_dir: Path) -> str:
        """Generate the ultimate HTML report for S. aureus"""
        print("\n🎨 Generating STAPHSCOPE ULTIMATE HTML report...")
        
        # Extract data
        samples_data = integrated_data.get('samples', {})
        patterns = integrated_data.get('patterns', {})
        gene_centric = integrated_data.get('gene_centric', {})
        metadata = integrated_data.get('metadata', {})
        
        # Create HTML
        html = self._create_ultimate_html(
            metadata=metadata,
            samples_data=samples_data,
            patterns=patterns,
            gene_centric=gene_centric,
            integrated_data=integrated_data
        )
        
        # Save HTML file
        output_file = output_dir / "staphscope_ultimate_report.html"
        with open(output_file, 'w', encoding='utf-8') as f:
            f.write(html)
        
        print(f"    ✅ HTML report saved: {output_file}")
        return str(output_file)
    
    def _create_ultimate_html(self, **kwargs) -> str:
        """Create ultimate HTML with all sections for S. aureus"""
        
        # CSS Styles - Updated with scrollable genome list and new tab colors
        css = """
        <style>
        :root {
            --summary-color: #4CAF50;
            --sample_overview-color: #2196F3;
            --qc-color: #607D8B;
            --mlst-color: #FF9800;
            --spa-color: #9C27B0;
            --sccmec-color: #009688;
            --mrsa-color: #795548;
            --amr-color: #F44336;
            --virulence-color: #E91E63;
            --bacmet-color: #FF5722;
            --plasmids-color: #673AB7;
            --patterns-color: #3F51B5;
            --aiguide-color: #00BCD4;
            --citation-color: #8BC34A;
            --calltoaction-color: #FFC107;
            --export-color: #9E9E9E;
        }
        
        * {
            margin: 0;
            padding: 0;
            box-sizing: border-box;
        }
        
        body {
            font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif;
            line-height: 1.6;
            color: #333;
            background: #f5f5f5;
            min-width: 1200px;
        }
        
        .container {
            max-width: none;
            margin: 0 auto;
            padding: 20px;
            width: 100%;
            overflow-x: auto;
        }
        
        .main-header {
            background: linear-gradient(135deg, #006400 0%, #228B22 100%);
            color: white;
            padding: 30px;
            border-radius: 15px;
            box-shadow: 0 10px 30px rgba(0,0,0,0.2);
            margin-bottom: 30px;
            text-align: center;
        }
        
        .main-header h1 {
            font-size: 2.8em;
            margin-bottom: 10px;
            color: white;
        }
        
        .metadata-bar {
            background: rgba(255,255,255,0.1);
            padding: 15px;
            border-radius: 10px;
            margin: 20px 0;
            display: flex;
            justify-content: space-around;
            flex-wrap: wrap;
            gap: 15px;
            backdrop-filter: blur(10px);
        }
        
        .metadata-item {
            display: flex;
            align-items: center;
            gap: 8px;
            font-size: 0.95em;
        }
        
        .dashboard-grid {
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(250px, 1fr));
            gap: 20px;
            margin-bottom: 30px;
        }
        
        .dashboard-card {
            background: white;
            padding: 25px;
            border-radius: 12px;
            box-shadow: 0 5px 20px rgba(0,0,0,0.1);
            text-align: center;
            transition: all 0.3s ease;
            cursor: pointer;
            border-left: 5px solid;
            position: relative;
            overflow: hidden;
        }
        
        .dashboard-card::before {
            content: '';
            position: absolute;
            top: 0;
            left: 0;
            right: 0;
            height: 3px;
            background: linear-gradient(90deg, transparent, rgba(255,255,255,0.8), transparent);
        }
        
        .dashboard-card:hover {
            transform: translateY(-10px);
            box-shadow: 0 15px 30px rgba(0,0,0,0.2);
        }
        
        .card-summary { border-left-color: var(--summary-color); }
        .card-sample_overview { border-left-color: var(--sample_overview-color); }
        .card-qc { border-left-color: var(--qc-color); }
        .card-mlst { border-left-color: var(--mlst-color); }
        .card-spa { border-left-color: var(--spa-color); }
        .card-sccmec { border-left-color: var(--sccmec-color); }
        .card-mrsa { border-left-color: var(--mrsa-color); }
        .card-amr { border-left-color: var(--amr-color); }
        .card-virulence { border-left-color: var(--virulence-color); }
        .card-bacmet { border-left-color: var(--bacmet-color); }
        .card-plasmids { border-left-color: var(--plasmids-color); }
        .card-patterns { border-left-color: var(--patterns-color); }
        .card-aiguide { border-left-color: var(--aiguide-color); }
        .card-citation { border-left-color: var(--citation-color); }
        .card-calltoaction { border-left-color: var(--calltoaction-color); }
        .card-export { border-left-color: var(--export-color); }
        
        .card-number {
            font-size: 3em;
            font-weight: bold;
            margin: 15px 0;
            background: linear-gradient(90deg, #006400, #228B22);
            -webkit-background-clip: text;
            -webkit-text-fill-color: transparent;
        }
        
        .tab-navigation {
            display: flex;
            gap: 5px;
            margin-bottom: 20px;
            flex-wrap: wrap;
            background: white;
            padding: 15px;
            border-radius: 12px;
            box-shadow: 0 5px 20px rgba(0,0,0,0.1);
            position: sticky;
            top: 10px;
            z-index: 100;
        }
        
        .tab-button {
            padding: 12px 20px;
            background: #f5f5f5;
            border: none;
            border-radius: 8px;
            cursor: pointer;
            font-weight: 600;
            color: #666;
            transition: all 0.3s ease;
            display: flex;
            align-items: center;
            gap: 8px;
            position: relative;
            overflow: hidden;
            font-size: 0.9em;
        }
        
        .tab-button::after {
            content: '';
            position: absolute;
            bottom: 0;
            left: 50%;
            right: 50%;
            height: 3px;
            background: currentColor;
            transition: all 0.3s ease;
        }
        
        .tab-button:hover::after {
            left: 10%;
            right: 10%;
        }
        
        .tab-button.active {
            color: white;
        }
        
        .tab-button.active::after {
            left: 10%;
            right: 10%;
        }
        
        .tab-button.summary.active { background: var(--summary-color); }
        .tab-button.sample_overview.active { background: var(--sample_overview-color); }
        .tab-button.qc.active { background: var(--qc-color); }
        .tab-button.mlst.active { background: var(--mlst-color); }
        .tab-button.spa.active { background: var(--spa-color); }
        .tab-button.sccmec.active { background: var(--sccmec-color); }
        .tab-button.mrsa.active { background: var(--mrsa-color); }
        .tab-button.amr.active { background: var(--amr-color); }
        .tab-button.virulence.active { background: var(--virulence-color); }
        .tab-button.bacmet.active { background: var(--bacmet-color); }
        .tab-button.plasmids.active { background: var(--plasmids-color); }
        .tab-button.patterns.active { background: var(--patterns-color); }
        .tab-button.aiguide.active { background: var(--aiguide-color); }
        .tab-button.citation.active { background: var(--citation-color); }
        .tab-button.calltoaction.active { background: var(--calltoaction-color); }
        .tab-button.export.active { background: var(--export-color); }
        
        .tab-content {
            display: none;
            background: white;
            padding: 30px;
            border-radius: 15px;
            box-shadow: 0 10px 30px rgba(0,0,0,0.1);
            margin-bottom: 30px;
            animation: fadeIn 0.5s ease;
            width: 100%;
            overflow-x: auto;
        }
        
        .tab-content.active {
            display: block;
        }
        
        @keyframes fadeIn {
            from { opacity: 0; transform: translateY(20px); }
            to { opacity: 1; transform: translateY(0); }
        }
        
        .section-header {
            color: #2c3e50;
            margin-bottom: 25px;
            padding-bottom: 15px;
            border-bottom: 3px solid;
            font-size: 1.8em;
            display: flex;
            align-items: center;
            justify-content: space-between;
        }
        
        .summary-header { border-color: var(--summary-color); }
        .sample_overview-header { border-color: var(--sample_overview-color); }
        .qc-header { border-color: var(--qc-color); }
        .mlst-header { border-color: var(--mlst-color); }
        .spa-header { border-color: var(--spa-color); }
        .sccmec-header { border-color: var(--sccmec-color); }
        .mrsa-header { border-color: var(--mrsa-color); }
        .amr-header { border-color: var(--amr-color); }
        .virulence-header { border-color: var(--virulence-color); }
        .bacmet-header { border-color: var(--bacmet-color); }
        .plasmids-header { border-color: var(--plasmids-color); }
        .patterns-header { border-color: var(--patterns-color); }
        .aiguide-header { border-color: var(--aiguide-color); }
        .citation-header { border-color: var(--citation-color); }
        .calltoaction-header { border-color: var(--calltoaction-color); }
        .export-header { border-color: var(--export-color); }
        
        .data-table {
            width: 100%;
            border-collapse: collapse;
            margin: 20px 0;
            font-size: 0.95em;
            box-shadow: 0 2px 10px rgba(0,0,0,0.1);
            border-radius: 8px;
            overflow: hidden;
            table-layout: auto;
        }
        
        .data-table th {
            background: #2c3e50;
            color: white;
            padding: 15px;
            text-align: left;
            font-weight: 600;
            position: sticky;
            top: 0;
            white-space: nowrap;
            cursor: pointer;
        }
        
        .data-table th:hover {
            background: #1a252f;
        }
        
        .data-table td {
            padding: 12px;
            border-bottom: 1px solid #e0e0e0;
            vertical-align: top;
            word-wrap: break-word;
            word-break: break-word;
            white-space: nowrap;
        }
        
        .data-table tr:hover {
            background: #f8f9fa;
        }
        
        .scrollable-table {
            max-height: none;
            overflow-y: auto;
            border: 1px solid #e0e0e0;
            border-radius: 8px;
            margin: 20px 0;
            width: 100%;
        }
        
        .master-scrollable-container {
            width: 100%;
            overflow-x: auto;
            border: 1px solid #e0e0e0;
            border-radius: 8px;
            margin: 20px 0;
        }
        
        /* Scrollable genome list - wrapping + vertical scroll */
        .genome-list {
            display: flex;
            flex-wrap: wrap;
            gap: 5px;
            max-height: 200px;
            overflow-y: auto;
            padding: 5px;
            background: #f8f9fa;
            border-radius: 5px;
        }

        .genome-tag {
            display: inline-block;
            background: #e6ffe6;
            color: #006400;
            padding: 3px 10px;
            border-radius: 12px;
            font-size: 0.85em;
            border: 1px solid #b3ffb3;
            white-space: nowrap;
            margin: 2px;
        }
        
        .genome-tag.highlight {
            background-color: #ffff99 !important;
            color: #000 !important;
            border: 1px solid #ffc107;
        }
        
        .search-box {
            width: 100%;
            padding: 12px;
            margin-bottom: 20px;
            border: 2px solid #e0e0e0;
            border-radius: 8px;
            font-size: 1em;
            transition: all 0.3s ease;
        }
        
        .search-box:focus {
            outline: none;
            border-color: #006400;
            box-shadow: 0 0 0 3px rgba(0, 100, 0, 0.1);
        }
        
        .badge {
            display: inline-block;
            padding: 5px 15px;
            border-radius: 20px;
            font-size: 0.85em;
            font-weight: 600;
            margin: 2px;
        }
        
        .badge-mrsa { background: #8B0000; color: white; }
        .badge-mssa { background: #4682B4; color: white; }
        .badge-critical { background: #DC143C; color: white; }
        .badge-high { background: #FF4500; color: white; }
        .badge-medium { background: #FF8C00; color: black; }
        .badge-low { background: #32CD32; color: white; }
        
        .alert-box {
            padding: 20px;
            border-radius: 10px;
            margin: 20px 0;
            display: flex;
            align-items: center;
            gap: 20px;
            border-left: 5px solid;
        }
        
        .alert-success { background: #d4edda; color: #155724; border-left-color: #28a745; }
        .alert-warning { background: #fff3cd; color: #856404; border-left-color: #ffc107; }
        .alert-danger { background: #f8d7da; color: #721c24; border-left-color: #dc3545; }
        .alert-info { background: #d1ecf1; color: #0c5460; border-left-color: #17a2b8; }
        
        .action-buttons {
            display: flex;
            gap: 10px;
            margin: 20px 0;
            flex-wrap: wrap;
        }
        
        .action-btn {
            padding: 10px 20px;
            border: none;
            border-radius: 8px;
            cursor: pointer;
            font-weight: 600;
            display: flex;
            align-items: center;
            gap: 8px;
            transition: all 0.3s ease;
        }
        
        .action-btn:hover {
            transform: translateY(-2px);
            box-shadow: 0 5px 15px rgba(0,0,0,0.2);
        }
        
        .btn-primary { background: #006400; color: white; }
        .btn-success { background: #28a745; color: white; }
        .btn-danger { background: #dc3545; color: white; }
        .btn-warning { background: #ffc107; color: black; }
        .btn-info { background: #17a2b8; color: white; }
        .btn-secondary { background: #6c757d; color: white; }
        .btn-light { background: #f8f9fa; color: #212529; border: 1px solid #dee2e6; }
        
        .database-section {
            margin: 30px 0;
            padding: 25px;
            border-radius: 12px;
            background: #f8f9fa;
            box-shadow: 0 3px 15px rgba(0,0,0,0.08);
        }
        
        .database-header {
            font-size: 1.4em;
            color: #2c3e50;
            margin-bottom: 20px;
            padding-bottom: 10px;
            border-bottom: 2px solid #006400;
            display: flex;
            align-items: center;
            justify-content: space-between;
        }
        
        .print-section-btn {
            background: #006400;
            color: white;
            border: none;
            border-radius: 5px;
            padding: 8px 15px;
            cursor: pointer;
            display: flex;
            align-items: center;
            gap: 5px;
            font-size: 0.9em;
        }
        
        .print-section-btn:hover {
            background: #228B22;
        }
        
        .footer {
            text-align: center;
            padding: 30px;
            color: white;
            margin-top: 40px;
            border-radius: 15px;
            background: linear-gradient(135deg, #2c3e50 0%, #34495e 100%);
        }
        
        .mrsa-highlight {
            background-color: #ffe6e6 !important;
            border-left: 3px solid #8B0000 !important;
        }
        
        .sort-icon {
            margin-left: 5px;
            font-size: 0.8em;
            opacity: 0.6;
        }
        
        @media print {
            body * {
                visibility: hidden;
            }
            .tab-content.active,
            .tab-content.active * {
                visibility: visible;
            }
            .tab-content.active {
                position: absolute;
                left: 0;
                top: 0;
                width: 100%;
                padding: 20px;
                box-shadow: none;
                border-radius: 0;
            }
            .print-section-btn,
            .tab-navigation,
            .dashboard-grid,
            .search-box,
            .action-buttons {
                display: none !important;
            }
            
            .data-table {
                page-break-inside: auto;
            }
            .data-table tr {
                page-break-inside: avoid;
                page-break-after: auto;
            }
            .data-table td, .data-table th {
                page-break-inside: avoid;
            }
        }
        
        @media (max-width: 768px) {
            .container {
                padding: 10px;
            }
            
            .main-header {
                padding: 20px;
            }
            
            .main-header h1 {
                font-size: 2em;
            }
            
            .tab-button {
                padding: 8px 12px;
                font-size: 0.8em;
            }
            
            .dashboard-grid {
                grid-template-columns: repeat(auto-fit, minmax(180px, 1fr));
            }
            
            .data-table {
                font-size: 0.8em;
            }
            
            body {
                min-width: auto;
                overflow-x: auto;
            }
        }
        </style>
        """
        
        # JavaScript - updated with genome highlight function
        js = """
        <script>
        // Tab switching
        function switchTab(tabName) {
            document.querySelectorAll('.tab-content').forEach(tab => {
                tab.classList.remove('active');
            });
            document.querySelectorAll('.tab-button').forEach(button => {
                button.classList.remove('active');
            });
            document.getElementById(tabName + '-tab').classList.add('active');
            event.currentTarget.classList.add('active');
            window.location.hash = tabName;
        }
        
        // Search functionality for tables
        function searchTable(tableId, searchId) {
            const input = document.getElementById(searchId);
            const filter = input.value.toUpperCase();
            const table = document.getElementById(tableId);
            const rows = table.getElementsByTagName('tr');
            
            for (let i = 1; i < rows.length; i++) {
                const cells = rows[i].getElementsByTagName('td');
                let found = false;
                for (let j = 0; j < cells.length; j++) {
                    const cell = cells[j];
                    if (cell) {
                        const txtValue = cell.textContent || cell.innerText;
                        if (txtValue.toUpperCase().indexOf(filter) > -1) {
                            found = true;
                            break;
                        }
                    }
                }
                rows[i].style.display = found ? '' : 'none';
            }
        }
        
        // Highlight genome tags
        function highlightGenome(tableId, searchId) {
            const filter = document.getElementById(searchId).value.toUpperCase().trim();
            const table = document.getElementById(tableId);
            const allTags = table.querySelectorAll('.genome-tag');
            allTags.forEach(tag => tag.classList.remove('highlight'));
            if (filter === '') return;
            allTags.forEach(tag => {
                if (tag.textContent.toUpperCase().indexOf(filter) > -1) {
                    tag.classList.add('highlight');
                }
            });
        }
        
        // Sort table function
        function sortTable(tableId, colIndex, type) {
            const table = document.getElementById(tableId);
            const tbody = table.tBodies[0];
            const rows = Array.from(tbody.rows);
            const isAscending = table.getAttribute('data-sort-dir') !== 'asc';
            
            rows.sort((a, b) => {
                let aVal = a.cells[colIndex].innerText.trim();
                let bVal = b.cells[colIndex].innerText.trim();
                
                if (type === 'number') {
                    aVal = parseFloat(aVal.replace(/,/g, '')) || 0;
                    bVal = parseFloat(bVal.replace(/,/g, '')) || 0;
                    return isAscending ? aVal - bVal : bVal - aVal;
                } else {
                    return isAscending ? aVal.localeCompare(bVal) : bVal.localeCompare(aVal);
                }
            });
            
            tbody.append(...rows);
            table.setAttribute('data-sort-dir', isAscending ? 'asc' : 'desc');
            
            // Update sort icons
            const headers = table.querySelectorAll('th');
            headers.forEach((th, idx) => {
                const icon = th.querySelector('.sort-icon');
                if (icon) icon.innerHTML = '⇅';
            });
            const currentHeader = headers[colIndex];
            const icon = currentHeader.querySelector('.sort-icon');
            if (icon) icon.innerHTML = isAscending ? '↑' : '↓';
        }
        
        // Print current section
        function printSection(sectionId) {
            const content = document.getElementById(sectionId);
            const printWindow = window.open('', '_blank');
            printWindow.document.write('<html><head><title>Print Section</title>');
            printWindow.document.write('<style>' + document.querySelector('style').textContent + '</style>');
            printWindow.document.write('</head><body>');
            printWindow.document.write(content.innerHTML);
            printWindow.document.write('</body></html>');
            printWindow.document.close();
            printWindow.print();
        }
        
        // Export table to CSV
        function exportTableToCSV(tableId, filename) {
            const table = document.getElementById(tableId);
            const rows = table.querySelectorAll('tr');
            const csv = [];
            for (let i = 0; i < rows.length; i++) {
                const row = [], cols = rows[i].querySelectorAll('td, th');
                for (let j = 0; j < cols.length; j++) {
                    row.push('"' + (cols[j].innerText || '').replace(/"/g, '""') + '"');
                }
                csv.push(row.join(','));
            }
            const csvFile = new Blob([csv.join('\\n')], {type: 'text/csv'});
            const downloadLink = document.createElement('a');
            downloadLink.download = filename;
            downloadLink.href = window.URL.createObjectURL(csvFile);
            downloadLink.style.display = 'none';
            document.body.appendChild(downloadLink);
            downloadLink.click();
            document.body.removeChild(downloadLink);
        }
        
        // Initialize from URL hash and add sort listeners
        document.addEventListener('DOMContentLoaded', function() {
            const hash = window.location.hash.substring(1);
            if (hash) {
                const tabButton = document.querySelector(`.tab-button.${hash}`);
                if (tabButton) {
                    tabButton.click();
                }
            } else {
                document.querySelector('.tab-button').click();
            }
            
            // Add sort listeners to all tables with sortable headers
            document.querySelectorAll('.data-table').forEach(table => {
                const headers = table.querySelectorAll('th');
                headers.forEach((header, idx) => {
                    const type = header.getAttribute('data-sort') || 'string';
                    header.style.cursor = 'pointer';
                    header.addEventListener('click', () => {
                        sortTable(table.id, idx, type);
                    });
                    // Add sort icon
                    const icon = document.createElement('span');
                    icon.className = 'sort-icon';
                    icon.innerHTML = '⇅';
                    header.appendChild(icon);
                });
            });
        });
        </script>
        """
        
        # Calculate total AMR genes count
        total_amr_genes = sum(len(genes) for genes in kwargs['gene_centric'].get('amr_databases', {}).values())
        total_virulence_genes = sum(len(genes) for genes in kwargs['gene_centric'].get('virulence_databases', {}).values())
        total_bacmet_genes = sum(len(genes) for genes in kwargs['gene_centric'].get('bacmet_databases', {}).values())
        
        # Build HTML
        html = f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>STAPHSCOPE Ultimate S. aureus Report</title>
    <link rel="stylesheet" href="https://cdnjs.cloudflare.com/ajax/libs/font-awesome/6.4.0/css/all.min.css">
    {css}
    {js}
</head>
<body>
    <div class="container">
        <!-- Main Header -->
        <div class="main-header">
            <h1><i class="fas fa-bacteria"></i> STAPHSCOPE Ultimate S. aureus Analysis Report</h1>
            <p>Gene-Centric Cross-Genome Analysis with Complete Genome Lists</p>
            
            <div class="metadata-bar">
                <div class="metadata-item">
                    <i class="fas fa-calendar"></i>
                    <span>Generated: {kwargs['metadata'].get('analysis_date', 'Unknown')}</span>
                </div>
                <div class="metadata-item">
                    <i class="fas fa-database"></i>
                    <span>Samples: {len(kwargs['samples_data'])}</span>
                </div>
                <div class="metadata-item">
                    <i class="fas fa-user-md"></i>
                    <span>Tool: STAPHSCOPE Ultimate v2.1.0</span>
                </div>
                <div class="metadata-item">
                    <i class="fas fa-university"></i>
                    <span>University of Ghana Medical School</span>
                </div>
            </div>
        </div>
        
        <!-- Dashboard -->
        <div class="dashboard-grid">
            <div class="dashboard-card card-summary" onclick="switchTab('summary')">
                <div class="card-number">{len(kwargs['samples_data'])}</div>
                <div class="card-label">Total Samples</div>
                <i class="fas fa-vial fa-2x" style="color: var(--summary-color); margin-top: 10px;"></i>
            </div>
            
            <div class="dashboard-card card-mlst" onclick="switchTab('mlst')">
                <div class="card-number">{len(kwargs['patterns'].get('mlst_distribution', {}))}</div>
                <div class="card-label">Unique STs</div>
                <i class="fas fa-code-branch fa-2x" style="color: var(--mlst-color); margin-top: 10px;"></i>
            </div>
            
            <div class="dashboard-card card-spa" onclick="switchTab('spa')">
                <div class="card-number">{len(kwargs['patterns'].get('spa_type_distribution', {}))}</div>
                <div class="card-label">spa Types</div>
                <i class="fas fa-dna fa-2x" style="color: var(--spa-color); margin-top: 10px;"></i>
            </div>
            
            <div class="dashboard-card card-amr" onclick="switchTab('amr')">
                <div class="card-number">{total_amr_genes}</div>
                <div class="card-label">AMR Genes</div>
                <i class="fas fa-biohazard fa-2x" style="color: var(--amr-color); margin-top: 10px;"></i>
            </div>
            
            <div class="dashboard-card card-virulence" onclick="switchTab('virulence')">
                <div class="card-number">{total_virulence_genes}</div>
                <div class="card-label">Virulence Genes</div>
                <i class="fas fa-virus fa-2x" style="color: var(--virulence-color); margin-top: 10px;"></i>
            </div>
            
            <div class="dashboard-card card-patterns" onclick="switchTab('patterns')">
                <div class="card-number">{len(kwargs['patterns'].get('high_risk_combinations', []))}</div>
                <div class="card-label">High-Risk Combos</div>
                <i class="fas fa-project-diagram fa-2x" style="color: var(--patterns-color); margin-top: 10px;"></i>
            </div>
        </div>
        
        <!-- Tab Navigation (reordered) -->
        <div class="tab-navigation">
            <button class="tab-button summary active" onclick="switchTab('summary')"><i class="fas fa-chart-pie"></i> Summary</button>
            <button class="tab-button sample_overview" onclick="switchTab('sample_overview')"><i class="fas fa-list-alt"></i> Sample Overview</button>
            <button class="tab-button qc" onclick="switchTab('qc')"><i class="fas fa-chart-line"></i> FASTA QC</button>
            <button class="tab-button mlst" onclick="switchTab('mlst')"><i class="fas fa-code-branch"></i> MLST</button>
            <button class="tab-button spa" onclick="switchTab('spa')"><i class="fas fa-dna"></i> spa Typing</button>
            <button class="tab-button sccmec" onclick="switchTab('sccmec')"><i class="fas fa-shield-alt"></i> SCCmec</button>
            <button class="tab-button mrsa" onclick="switchTab('mrsa')"><i class="fas fa-skull-crossbones"></i> MRSA Analysis</button>
            <button class="tab-button amr" onclick="switchTab('amr')"><i class="fas fa-biohazard"></i> AMR</button>
            <button class="tab-button virulence" onclick="switchTab('virulence')"><i class="fas fa-virus"></i> Virulence</button>
            <button class="tab-button bacmet" onclick="switchTab('bacmet')"><i class="fas fa-flask"></i> BACMET</button>
            <button class="tab-button plasmids" onclick="switchTab('plasmids')"><i class="fas fa-plug"></i> Plasmids</button>
            <button class="tab-button patterns" onclick="switchTab('patterns')"><i class="fas fa-project-diagram"></i> Patterns</button>
            <button class="tab-button aiguide" onclick="switchTab('aiguide')"><i class="fas fa-robot"></i> AI Guide</button>
            <button class="tab-button citation" onclick="switchTab('citation')"><i class="fas fa-book"></i> Citation</button>
            <button class="tab-button calltoaction" onclick="switchTab('calltoaction')"><i class="fas fa-bullhorn"></i> Call to Action</button>
            <button class="tab-button export" onclick="switchTab('export')"><i class="fas fa-download"></i> Export</button>
        </div>
        
        <!-- Summary Tab -->
        <div id="summary-tab" class="tab-content active">
            <h2 class="section-header summary-header">
                <i class="fas fa-chart-pie"></i> Executive Summary
                <button class="print-section-btn" onclick="printSection('summary-tab')">
                    <i class="fas fa-print"></i> Print
                </button>
            </h2>
            {self._generate_summary_section(kwargs)}
        </div>
        
        <!-- Sample Overview Tab -->
        <div id="sample_overview-tab" class="tab-content">
            <h2 class="section-header sample_overview-header">
                <i class="fas fa-list-alt"></i> Sample Overview
                <button class="print-section-btn" onclick="printSection('sample_overview-tab')">
                    <i class="fas fa-print"></i> Print
                </button>
            </h2>
            {self._generate_sample_overview_section(kwargs)}
        </div>
        
        <!-- FASTA QC Tab -->
        <div id="qc-tab" class="tab-content">
            <h2 class="section-header qc-header">
                <i class="fas fa-chart-line"></i> FASTA Quality Control Metrics
                <button class="print-section-btn" onclick="printSection('qc-tab')">
                    <i class="fas fa-print"></i> Print
                </button>
            </h2>
            {self._generate_qc_section(kwargs)}
        </div>
        
        <!-- MLST Tab -->
        <div id="mlst-tab" class="tab-content">
            <h2 class="section-header mlst-header">
                <i class="fas fa-code-branch"></i> MLST Analysis
                <button class="print-section-btn" onclick="printSection('mlst-tab')">
                    <i class="fas fa-print"></i> Print
                </button>
            </h2>
            {self._generate_mlst_section(kwargs)}
        </div>
        
        <!-- spa Typing Tab -->
        <div id="spa-tab" class="tab-content">
            <h2 class="section-header spa-header">
                <i class="fas fa-dna"></i> spa Typing Analysis
                <button class="print-section-btn" onclick="printSection('spa-tab')">
                    <i class="fas fa-print"></i> Print
                </button>
            </h2>
            {self._generate_spa_section(kwargs)}
        </div>
        
        <!-- SCCmec Tab -->
        <div id="sccmec-tab" class="tab-content">
            <h2 class="section-header sccmec-header">
                <i class="fas fa-shield-alt"></i> SCCmec Typing Analysis
                <button class="print-section-btn" onclick="printSection('sccmec-tab')">
                    <i class="fas fa-print"></i> Print
                </button>
            </h2>
            {self._generate_sccmec_section(kwargs)}
        </div>
        
        <!-- MRSA Analysis Tab -->
        <div id="mrsa-tab" class="tab-content">
            <h2 class="section-header mrsa-header">
                <i class="fas fa-skull-crossbones"></i> MRSA Analysis
                <button class="print-section-btn" onclick="printSection('mrsa-tab')">
                    <i class="fas fa-print"></i> Print
                </button>
            </h2>
            {self._generate_mrsa_section(kwargs)}
        </div>
        
        <!-- AMR Genes Tab -->
        <div id="amr-tab" class="tab-content">
            <h2 class="section-header amr-header">
                <i class="fas fa-biohazard"></i> Antimicrobial Resistance Genes
                <button class="print-section-btn" onclick="printSection('amr-tab')">
                    <i class="fas fa-print"></i> Print
                </button>
            </h2>
            {self._generate_amr_section(kwargs)}
        </div>
        
        <!-- Virulence Genes Tab -->
        <div id="virulence-tab" class="tab-content">
            <h2 class="section-header virulence-header">
                <i class="fas fa-virus"></i> Virulence Genes
                <button class="print-section-btn" onclick="printSection('virulence-tab')">
                    <i class="fas fa-print"></i> Print
                </button>
            </h2>
            {self._generate_virulence_section(kwargs)}
        </div>
        
        <!-- BACMET (Biocides) Tab -->
        <div id="bacmet-tab" class="tab-content">
            <h2 class="section-header bacmet-header">
                <i class="fas fa-flask"></i> BACMET: Biocide & Heavy Metal Resistance
                <button class="print-section-btn" onclick="printSection('bacmet-tab')">
                    <i class="fas fa-print"></i> Print
                </button>
            </h2>
            {self._generate_bacmet_section(kwargs)}
        </div>
        
        <!-- Plasmids Tab -->
        <div id="plasmids-tab" class="tab-content">
            <h2 class="section-header plasmids-header">
                <i class="fas fa-plug"></i> Plasmid Replicon Analysis
                <button class="print-section-btn" onclick="printSection('plasmids-tab')">
                    <i class="fas fa-print"></i> Print
                </button>
            </h2>
            {self._generate_plasmids_section(kwargs)}
        </div>
        
        <!-- Pattern Discovery Tab -->
        <div id="patterns-tab" class="tab-content">
            <h2 class="section-header patterns-header">
                <i class="fas fa-project-diagram"></i> Cross-Genome Pattern Discovery
                <button class="print-section-btn" onclick="printSection('patterns-tab')">
                    <i class="fas fa-print"></i> Print
                </button>
            </h2>
            {self._generate_pattern_discovery_section(kwargs)}
        </div>
        
        <!-- AI Guide Tab -->
        <div id="aiguide-tab" class="tab-content">
            <h2 class="section-header aiguide-header">
                <i class="fas fa-robot"></i> AI Assistant Guide
                <button class="print-section-btn" onclick="printSection('aiguide-tab')">
                    <i class="fas fa-print"></i> Print
                </button>
            </h2>
            {self._generate_aiguide_section(kwargs)}
        </div>
        
        <!-- Citation Tab -->
        <div id="citation-tab" class="tab-content">
            <h2 class="section-header citation-header">
                <i class="fas fa-book"></i> Citations & References
                <button class="print-section-btn" onclick="printSection('citation-tab')">
                    <i class="fas fa-print"></i> Print
                </button>
            </h2>
            {self._generate_citation_section(kwargs)}
        </div>
        
        <!-- Call to Action Tab -->
        <div id="calltoaction-tab" class="tab-content">
            <h2 class="section-header calltoaction-header">
                <i class="fas fa-bullhorn"></i> Call to Action
                <button class="print-section-btn" onclick="printSection('calltoaction-tab')">
                    <i class="fas fa-print"></i> Print
                </button>
            </h2>
            {self._generate_calltoaction_section(kwargs)}
        </div>
        
        <!-- Export Tab -->
        <div id="export-tab" class="tab-content">
            <h2 class="section-header export-header">
                <i class="fas fa-download"></i> Export Data
                <button class="print-section-btn" onclick="printSection('export-tab')">
                    <i class="fas fa-print"></i> Print
                </button>
            </h2>
            {self._generate_export_section(kwargs)}
        </div>
        
        <!-- Footer -->
        <div class="footer">
            <h3>STAPHSCOPE Ultimate S. aureus Reporter v2.1.0</h3>
            <p>University of Ghana Medical School | Brown Beckley <brownbeckley94@gmail.com></p>
            <p>Generated on {kwargs['metadata'].get('analysis_date', 'Unknown')}</p>
            <p>Please give a big STAR on GitHub if you found the integrated report useful.</p>
        </div>
    </div>
</body>
</html>
        """
        
        return html
    
    def _generate_summary_section(self, kwargs: Dict) -> str:
        """Generate summary section for S. aureus with biology info"""
        samples_data = kwargs['samples_data']
        patterns = kwargs['patterns']
        gene_centric = kwargs['gene_centric']
        
        total_samples = len(samples_data)
        total_amr_genes = sum(len(genes) for genes in gene_centric.get('amr_databases', {}).values())
        total_virulence_genes = sum(len(genes) for genes in gene_centric.get('virulence_databases', {}).values())
        total_plasmids = sum(len(genes) for genes in gene_centric.get('plasmid_databases', {}).values())
        total_bacmet = sum(len(genes) for genes in gene_centric.get('bacmet_databases', {}).values())
        critical_findings = len(patterns.get('high_risk_combinations', []))
        
        # Count MRSA vs MSSA
        mrsa_count = 0
        mssa_count = 0
        for sample_data in samples_data.values():
            status = sample_data.get('typing', {}).get('MRSA_Status', '')
            if 'MRSA' in status:
                mrsa_count += 1
            elif 'MSSA' in status:
                mssa_count += 1
        
        html = f"""
        <div class="alert-box alert-info">
            <i class="fas fa-info-circle fa-2x"></i>
            <div>
                <h3>📊 Analysis Overview – Gene‑Centric Approach</h3>
                <p>This report analyzes <strong>{total_samples}</strong> <em>Staphylococcus aureus</em> genomes using a <strong>gene‑centric</strong> approach: each resistance or virulence gene is displayed with <strong>all genomes</strong> that carry it. This makes outbreak tracking and co‑occurrence analysis immediate and transparent.</p>
                <p><strong>Why gene‑centric?</strong> Traditional sample‑centric tables force you to search sample by sample. Here, you see the full distribution of every gene in seconds.</p>
            </div>
        </div>
        
        <div class="action-buttons">
            <button class="action-btn btn-primary" onclick="switchTab('amr')">
                <i class="fas fa-biohazard"></i> View AMR Genes
            </button>
            <button class="action-btn btn-success" onclick="switchTab('virulence')">
                <i class="fas fa-virus"></i> View Virulence Genes
            </button>
            <button class="action-btn btn-danger" onclick="switchTab('patterns')">
                <i class="fas fa-exclamation-triangle"></i> Check High-Risk Combos
            </button>
        </div>
        
        <h3><i class="fas fa-chart-bar"></i> Key Statistics</h3>
        <div class="scrollable-table">
            <table class="data-table">
                <thead>
                    <tr>
                        <th>Metric</th>
                        <th>Count</th>
                        <th>Details</th>
                    </tr>
                </thead>
                <tbody>
                    <tr><td>Total Samples Analyzed</td><td><strong>{total_samples}</strong></td><td>Complete genomic analysis with all databases</td></tr>
                    <tr><td>MRSA Samples</td><td><span class="badge badge-mrsa">{mrsa_count}</span></td><td>Methicillin‑resistant S. aureus (carry mecA/mecC)</td></tr>
                    <tr><td>MSSA Samples</td><td><span class="badge badge-mssa">{mssa_count}</span></td><td>Methicillin‑sensitive S. aureus</td></tr>
                    <tr><td>Unique MLST Types</td><td><strong>{len(patterns.get('mlst_distribution', {}))}</strong></td><td>Sequence types (population structure)</td></tr>
                    <tr><td>Unique spa Types</td><td><strong>{len(patterns.get('spa_type_distribution', {}))}</strong></td><td>Protein A gene typing – outbreak tracking</td></tr>
                    <tr><td>Unique SCCmec Types</td><td><strong>{len(patterns.get('sccmec_distribution', {}))}</strong></td><td>SCCmec cassette – MRSA lineage marker</td></tr>
                    <tr><td>AMR Genes</td><td><strong>{total_amr_genes}</strong></td><td>Across AMRfinder, CARD, ResFinder, etc.</td></tr>
                    <tr><td>Virulence Genes</td><td><strong>{total_virulence_genes}</strong></td><td>From VFDB (toxins, adhesins, immune evasion)</td></tr>
                    <tr><td>BACMET (Biocide/Metal) Genes</td><td><strong>{total_bacmet}</strong></td><td>Biocide and heavy metal resistance (hospital environment)</td></tr>
                    <tr><td>Plasmid Replicons</td><td><strong>{total_plasmids}</strong></td><td>Plasmid types – horizontal gene transfer potential</td></tr>
                    <tr><td>High‑Risk AMR+Virulence Combos</td><td><span class="badge {'badge-critical' if critical_findings > 0 else 'badge-low'}">{critical_findings}</span></td><td>Samples with both critical AMR and virulence genes</td></tr>
                </tbody>
            </table>
        </div>
        
        <h3 style="margin-top: 30px;"><i class="fas fa-lightbulb"></i> How to Use This Report</h3>
        <div style="display: grid; grid-template-columns: repeat(auto-fit, minmax(300px, 1fr)); gap: 20px; margin: 20px 0;">
            <div class="database-section">
                <h4><i class="fas fa-gene"></i> Gene‑Centric Tables</h4>
                <p>In AMR, Virulence, and BACMET tabs, each gene is shown with <strong>all genomes</strong> that contain it. Use the <strong>genome search box</strong> to highlight specific isolates across multiple genes.</p>
            </div>
            <div class="database-section">
                <h4><i class="fas fa-print"></i> Section‑Specific Printing</h4>
                <p>Click the <strong>Print</strong> button on any section header to generate a clean printout of only that tab – perfect for reports.</p>
            </div>
            <div class="database-section">
                <h4><i class="fas fa-search"></i> Smart Filter Buttons</h4>
                <p>Use the coloured buttons to instantly filter for <strong>mecA, PVL, TSST‑1, enterotoxins</strong> and other key genes.</p>
            </div>
            <div class="database-section">
                <h4><i class="fas fa-chart-line"></i> FASTA QC</h4>
                <p>Assess assembly quality (N50, contig count, GC%) before drawing biological conclusions.</p>
            </div>
        </div>
        """
        return html
    
    def _generate_sample_overview_section(self, kwargs: Dict) -> str:
        """Generate sample overview with biology info"""
        samples_data = kwargs['samples_data']
        
        html = f"""
        <div class="alert-box alert-info">
            <i class="fas fa-info-circle fa-2x"></i>
            <div>
                <h3>🧬 Sample Overview – Population Snapshot</h3>
                <p>This table summarises the three main typing methods for each <em>S. aureus</em> isolate:</p>
                <ul>
                    <li><strong>MLST</strong> (Multi‑Locus Sequence Typing) – gold standard for global epidemiology.</li>
                    <li><strong>spa typing</strong> – high‑resolution typing based on the <em>spa</em> gene (protein A repeat region).</li>
                    <li><strong>SCCmec type</strong> – identifies the staphylococcal cassette chromosome mec, which carries <em>mecA</em> and defines MRSA.</li>
                </ul>
                <p><strong>MRSA</strong> samples are highlighted in red. Click on column headers to sort.</p>
            </div>
        </div>
        
        <input type="text" class="search-box" id="search-samples" 
               onkeyup="searchTable('samples-table', 'search-samples')" 
               placeholder="🔍 Search samples by ID, ST, spa type, or SCCmec...">
        
        <div class="action-buttons">
            <button class="action-btn btn-primary" onclick="exportTableToCSV('samples-table', 'sample_overview.csv')">
                <i class="fas fa-download"></i> Export to CSV
            </button>
            <button class="action-btn btn-success" onclick="document.getElementById('search-samples').value=''; searchTable('samples-table', 'search-samples')">
                <i class="fas fa-sync"></i> Clear Search
            </button>
        </div>
        
        <div class="scrollable-table">
            <table id="samples-table" class="data-table">
                <thead>
                    <tr>
                        <th data-sort="string">Sample</th>
                        <th data-sort="string">MLST</th>
                        <th data-sort="string">spa Type</th>
                        <th data-sort="string">SCCmec Type</th>
                        <th data-sort="string">MRSA Status</th>
                        <th data-sort="number">Virulence Gene Count</th>
                    </tr>
                </thead>
                <tbody>
        """
        
        for sample, data in samples_data.items():
            mlst = data.get('typing', {}).get('MLST', 'ND')
            spa_type = data.get('typing', {}).get('spa_Type', 'ND')
            sccmec_type = data.get('typing', {}).get('SCCmec_Type', 'ND')
            mrsa_status = data.get('typing', {}).get('MRSA_Status', 'ND')
            virulence_count = len(data.get('abricate_databases', {}).get('vfdb', []))
            
            row_class = 'class="mrsa-highlight"' if 'MRSA' in mrsa_status else ''
            status_badge = '<span class="badge badge-mrsa">MRSA</span>' if 'MRSA' in mrsa_status else ('<span class="badge badge-mssa">MSSA</span>' if 'MSSA' in mrsa_status else mrsa_status)
            
            html += f"""
                    <tr {row_class}>
                        <td><strong>{sample}</strong></td>
                        <td>{mlst}</td><td>{spa_type}</td><td>{sccmec_type}</td>
                        <td>{status_badge}</td><td>{virulence_count}</td>
                    </tr>
            """
        
        html += """
                </tbody>
            </table>
        </div>
        """
        return html
    
    def _generate_qc_section(self, kwargs: Dict) -> str:
        """Generate FASTA QC section with biology explanation"""
        qc_data = kwargs.get('integrated_data', {}).get('qc_data', {})
        
        if not qc_data:
            return """
            <div class="alert-box alert-warning">
                <i class="fas fa-exclamation-circle fa-2x"></i>
                <div>
                    <h3>No QC Data Available</h3>
                    <p>The FASTA_QC_summary.html file was not found or could not be parsed. Please ensure it is present in the input directory.</p>
                </div>
            </div>
            """
        
        # Collect metrics
        all_metrics = set()
        for metrics in qc_data.values():
            all_metrics.update(metrics.keys())
        metric_list = sorted(all_metrics)
        
        html = f"""
        <div class="alert-box alert-info">
            <i class="fas fa-chart-line fa-2x"></i>
            <div>
                <h3>📏 FASTA Quality Control – Why It Matters</h3>
                <p>Assembly quality directly affects the reliability of typing and gene detection.</p>
                <ul>
                    <li><strong>Number of contigs</strong> – lower is better; fragmented assemblies may miss genes or break up SCCmec cassettes.</li>
                    <li><strong>N50</strong> – the contig length such that 50% of the assembly is in contigs ≥ N50; higher N50 indicates better continuity.</li>
                    <li><strong>GC%</strong> – S. aureus is typically 32‑33%; large deviations may indicate contamination.</li>
                    <li><strong>Total length</strong> – should be ~2.8 Mbp for a complete genome.</li>
                </ul>
            </div>
        </div>
        
        <input type="text" class="search-box" id="search-qc" onkeyup="searchTable('qc-table', 'search-qc')" placeholder="🔍 Search sample...">
        
        <div class="action-buttons">
            <button class="action-btn btn-primary" onclick="exportTableToCSV('qc-table', 'fasta_qc.csv')">
                <i class="fas fa-download"></i> Export QC Data
            </button>
        </div>
        
        <div class="master-scrollable-container">
            <table id="qc-table" class="data-table">
                <thead><tr><th data-sort="string">Sample</th>
        """
        for metric in metric_list:
            html += f'<th data-sort="number">{metric}</th>'
        html += "</tr></thead><tbody>"
        
        for sample, metrics in sorted(qc_data.items()):
            html += f"<tr><td><strong>{sample}</strong></td>"
            for metric in metric_list:
                val = metrics.get(metric, 'ND')
                if isinstance(val, float):
                    if val > 1e6:
                        val = f"{val:,.0f}"
                    elif val > 1000:
                        val = f"{val:,.0f}"
                    else:
                        val = f"{val:.2f}"
                html += f"<td>{val}</td>"
            html += "</tr>"
        
        html += """
                </tbody>
            </table>
        </div>
        """
        return html
    
    def _generate_mlst_section(self, kwargs: Dict) -> str:
        """Generate MLST analysis section with ST-spa and ST-SCCmec combinations"""
        patterns = kwargs['patterns']
        mlst_dist = patterns.get('mlst_distribution', Counter())
        mlst_spa_combos = patterns.get('mlst_spa_combinations', {})
        mlst_sccmec_combos = patterns.get('mlst_sccmec_combinations', {})
        
        html = f"""
        <div class="alert-box alert-info">
            <i class="fas fa-code-branch fa-2x"></i>
            <div>
                <h3>🔬 MLST (Multi‑Locus Sequence Typing)</h3>
                <p>MLST indexes internal fragments of seven housekeeping genes (<em>arcC, aroE, glpF, gmk, pta, tpi, yqiL</em>). Each unique combination of alleles defines a Sequence Type (ST).</p>
                <p><strong>Why MLST?</strong> It is highly reproducible and portable, enabling global surveillance. Closely related STs belong to the same clonal complex (CC). Major epidemic clones include CC1, CC5, CC8, CC22, CC30 and CC45.</p>
                <p><strong>{len(mlst_dist)} unique STs</strong> identified.</p>
            </div>
        </div>
        
        <h3>📊 ST Distribution</h3>
        <div class="scrollable-table">
            <table id="mlst-table" class="data-table">
                <thead><tr><th data-sort="string">ST</th><th data-sort="number">Count</th><th data-sort="number">Percentage</th><th data-sort="string">Associated spa Types</th><th data-sort="string">Associated SCCmec Types</th></tr></thead>
                <tbody>
        """
        total = sum(mlst_dist.values())
        for mlst, count in mlst_dist.most_common():
            if mlst == 'ND': continue
            pct = (count/total)*100
            # Find associated spa
            associated_spa = []
            for combo in mlst_spa_combos:
                if f"{mlst} - " in combo:
                    spa = combo.split(" - ")[1]
                    if spa not in associated_spa:
                        associated_spa.append(spa)
            # Find associated SCCmec
            associated_sccmec = []
            for combo in mlst_sccmec_combos:
                if f"{mlst} - " in combo:
                    scc = combo.split(" - ")[1]
                    if scc not in associated_sccmec:
                        associated_sccmec.append(scc)
            spa_list = ', '.join(associated_spa) if associated_spa else 'ND'
            scc_list = ', '.join(associated_sccmec) if associated_sccmec else 'ND'
            html += f"<tr><td><strong>{mlst}</strong></td><td>{count}</td><td>{pct:.1f}%</td><td>{spa_list}</td><td>{scc_list}</td></tr>"
        html += "</tbody></table></div>"
        
        # ST-spa combinations
        html += f"""
        <h3>🔗 ST–spa Combinations</h3>
        <input type="text" class="search-box" id="search-mlst-spa" onkeyup="searchTable('mlst-spa-table','search-mlst-spa')" placeholder="🔍 Search ST-spa...">
        <div class="master-scrollable-container"><table id="mlst-spa-table" class="data-table"><thead><tr><th data-sort="string">ST-spa Combination</th><th data-sort="number">Count</th><th data-sort="string">Samples</th></tr></thead><tbody>
        """
        for combo, samples in sorted(mlst_spa_combos.items(), key=lambda x: len(x[1]), reverse=True):
            sample_tags = ''.join(f'<span class="genome-tag">{s}</span>' for s in samples)
            html += f"<tr><td><strong>{combo}</strong></td><td>{len(samples)}</td><td><div class='genome-list'>{sample_tags}</div></td></tr>"
        html += "</tbody></table></div>"
        
        # ST-SCCmec combinations
        html += f"""
        <h3>🔗 ST–SCCmec Combinations</h3>
        <input type="text" class="search-box" id="search-mlst-sccmec" onkeyup="searchTable('mlst-sccmec-table','search-mlst-sccmec')" placeholder="🔍 Search ST-SCCmec...">
        <div class="master-scrollable-container"><table id="mlst-sccmec-table" class="data-table"><thead><tr><th data-sort="string">ST-SCCmec Combination</th><th data-sort="number">Count</th><th data-sort="string">Samples</th></tr></thead><tbody>
        """
        for combo, samples in sorted(mlst_sccmec_combos.items(), key=lambda x: len(x[1]), reverse=True):
            sample_tags = ''.join(f'<span class="genome-tag">{s}</span>' for s in samples)
            html += f"<tr><td><strong>{combo}</strong></td><td>{len(samples)}</td><td><div class='genome-list'>{sample_tags}</div></td></tr>"
        html += "</tbody></table></div>"
        return html
    
    def _generate_spa_section(self, kwargs: Dict) -> str:
        """Generate spa typing section"""
        patterns = kwargs['patterns']
        spa_dist = patterns.get('spa_type_distribution', Counter())
        mlst_spa_combos = patterns.get('mlst_spa_combinations', {})
        spa_sccmec_combos = patterns.get('spa_sccmec_combinations', {})
        
        html = f"""
        <div class="alert-box alert-info">
            <i class="fas fa-dna fa-2x"></i>
            <div>
                <h3>🧬 spa Typing – High‑Resolution Outbreak Tracking</h3>
                <p>The <em>spa</em> gene encodes protein A, an important virulence factor. The repeat region is highly polymorphic; differences in repeat patterns define spa types.</p>
                <p><strong>Why spa typing?</strong> It is faster and cheaper than MLST, with higher discrimination. It is widely used for outbreak investigation and surveillance of MRSA.</p>
                <p><strong>{len(spa_dist)} unique spa types</strong> identified.</p>
            </div>
        </div>
        
        <h3>📊 spa Type Distribution</h3>
        <div class="scrollable-table"><table id="spa-table" class="data-table"><thead><tr><th data-sort="string">spa Type</th><th data-sort="number">Count</th><th data-sort="number">Percentage</th><th data-sort="string">Common STs</th></tr></thead><tbody>
        """
        total = sum(spa_dist.values())
        for spa, count in spa_dist.most_common():
            if spa == 'ND': continue
            pct = (count/total)*100
            sts = []
            for combo in mlst_spa_combos:
                if spa in combo:
                    st = combo.split(" - ")[0]
                    if st not in sts:
                        sts.append(st)
            st_list = ', '.join(sts) if sts else 'None'
            html += f"<tr><td><strong>{spa}</strong></td><td>{count}</td><td>{pct:.1f}%</td><td>{st_list}</td></tr>"
        html += "</tbody></table></div>"
        
        # spa-ST combinations (reverse order)
        html += f"""
        <h3>🔗 spa–ST Combinations</h3>
        <input type="text" class="search-box" id="search-spa-mlst" onkeyup="searchTable('spa-mlst-table','search-spa-mlst')" placeholder="🔍 Search spa-ST...">
        <div class="master-scrollable-container"><table id="spa-mlst-table" class="data-table"><thead><tr><th data-sort="string">spa-ST Combination</th><th data-sort="number">Count</th><th data-sort="string">Samples</th></tr></thead><tbody>
        """
        for combo, samples in sorted(mlst_spa_combos.items(), key=lambda x: len(x[1]), reverse=True):
            parts = combo.split(" - ")
            rev_combo = f"{parts[1]} - {parts[0]}" if len(parts)==2 else combo
            sample_tags = ''.join(f'<span class="genome-tag">{s}</span>' for s in samples)
            html += f"<tr><td><strong>{rev_combo}</strong></td><td>{len(samples)}</td><td><div class='genome-list'>{sample_tags}</div></td></tr>"
        html += "</tbody></table></div>"
        
        # spa-SCCmec combinations
        html += f"""
        <h3>🔗 spa–SCCmec Combinations</h3>
        <input type="text" class="search-box" id="search-spa-sccmec" onkeyup="searchTable('spa-sccmec-table','search-spa-sccmec')" placeholder="🔍 Search spa-SCCmec...">
        <div class="master-scrollable-container"><table id="spa-sccmec-table" class="data-table"><thead><tr><th data-sort="string">spa-SCCmec Combination</th><th data-sort="number">Count</th><th data-sort="string">Samples</th></tr></thead><tbody>
        """
        for combo, samples in sorted(spa_sccmec_combos.items(), key=lambda x: len(x[1]), reverse=True):
            sample_tags = ''.join(f'<span class="genome-tag">{s}</span>' for s in samples)
            html += f"<tr><td><strong>{combo}</strong></td><td>{len(samples)}</td><td><div class='genome-list'>{sample_tags}</div></td></tr>"
        html += "</tbody></table></div>"
        return html
    
    def _generate_sccmec_section(self, kwargs: Dict) -> str:
        """Generate SCCmec section"""
        patterns = kwargs['patterns']
        sccmec_dist = patterns.get('sccmec_distribution', Counter())
        mlst_sccmec_combos = patterns.get('mlst_sccmec_combinations', {})
        spa_sccmec_combos = patterns.get('spa_sccmec_combinations', {})
        
        html = f"""
        <div class="alert-box alert-info">
            <i class="fas fa-shield-alt fa-2x"></i>
            <div>
                <h3>🧬 SCCmec Typing – The MRSA Cassette</h3>
                <p>The staphylococcal cassette chromosome <em>mec</em> (SCCmec) carries the <em>mecA</em> or <em>mecC</em> gene, conferring methicillin resistance. SCCmec types are classified by their <em>ccr</em> (cassette chromosome recombinase) gene complex and <em>mec</em> gene complex.</p>
                <p>Major types: <strong>I, II, III</strong> (hospital‑acquired MRSA, often large), <strong>IV, V</strong> (community‑acquired MRSA, smaller, more mobile).</p>
                <p><strong>{len(sccmec_dist)} unique SCCmec types</strong> identified.</p>
            </div>
        </div>
        
        <h3>📊 SCCmec Type Distribution</h3>
        <div class="scrollable-table"><table id="sccmec-table" class="data-table"><thead><tr><th data-sort="string">SCCmec Type</th><th data-sort="number">Count</th><th data-sort="number">Percentage</th><th data-sort="string">Common STs</th><th data-sort="string">Common spa Types</th></tr></thead><tbody>
        """
        total = sum(sccmec_dist.values())
        for scc, count in sccmec_dist.most_common():
            if scc in ['ND', 'Not Assigned']: continue
            pct = (count/total)*100
            sts = [c.split(" - ")[0] for c in mlst_sccmec_combos if scc in c]
            spas = [c.split(" - ")[0] for c in spa_sccmec_combos if scc in c]
            st_list = ', '.join(set(sts)) if sts else 'None'
            spa_list = ', '.join(set(spas)) if spas else 'None'
            html += f"<tr><td><strong>{scc}</strong></td><td>{count}</td><td>{pct:.1f}%</td><td>{st_list}</td><td>{spa_list}</td></tr>"
        html += "</tbody></table></div>"
        
        # SCCmec-ST combinations
        html += f"""
        <h3>🔗 SCCmec–ST Combinations</h3>
        <input type="text" class="search-box" id="search-sccmec-mlst" onkeyup="searchTable('sccmec-mlst-table','search-sccmec-mlst')" placeholder="🔍 Search SCCmec-ST...">
        <div class="master-scrollable-container"><table id="sccmec-mlst-table" class="data-table"><thead><tr><th data-sort="string">SCCmec-ST Combination</th><th data-sort="number">Count</th><th data-sort="string">Samples</th></tr></thead><tbody>
        """
        for combo, samples in sorted(mlst_sccmec_combos.items(), key=lambda x: len(x[1]), reverse=True):
            parts = combo.split(" - ")
            rev_combo = f"{parts[1]} - {parts[0]}" if len(parts)==2 else combo
            sample_tags = ''.join(f'<span class="genome-tag">{s}</span>' for s in samples)
            html += f"<tr><td><strong>{rev_combo}</strong></td><td>{len(samples)}</td><td><div class='genome-list'>{sample_tags}</div></td></tr>"
        html += "</tbody></table></div>"
        
        # SCCmec-spa combinations
        html += f"""
        <h3>🔗 SCCmec–spa Combinations</h3>
        <input type="text" class="search-box" id="search-sccmec-spa" onkeyup="searchTable('sccmec-spa-table','search-sccmec-spa')" placeholder="🔍 Search SCCmec-spa...">
        <div class="master-scrollable-container"><table id="sccmec-spa-table" class="data-table"><thead><tr><th data-sort="string">SCCmec-spa Combination</th><th data-sort="number">Count</th><th data-sort="string">Samples</th></tr></thead><tbody>
        """
        for combo, samples in sorted(spa_sccmec_combos.items(), key=lambda x: len(x[1]), reverse=True):
            parts = combo.split(" - ")
            rev_combo = f"{parts[1]} - {parts[0]}" if len(parts)==2 else combo
            sample_tags = ''.join(f'<span class="genome-tag">{s}</span>' for s in samples)
            html += f"<tr><td><strong>{rev_combo}</strong></td><td>{len(samples)}</td><td><div class='genome-list'>{sample_tags}</div></td></tr>"
        html += "</tbody></table></div>"
        return html
    
    def _generate_mrsa_section(self, kwargs: Dict) -> str:
        """Generate MRSA-specific analysis with triple combinations"""
        patterns = kwargs['patterns']
        mrsa_status_dist = patterns.get('mrsa_status_distribution', Counter())
        mlst_spa_combos = patterns.get('mlst_spa_combinations', {})
        mlst_sccmec_combos = patterns.get('mlst_sccmec_combinations', {})
        spa_sccmec_combos = patterns.get('spa_sccmec_combinations', {})
        samples_data = kwargs['samples_data']
        
        # Filter MRSA samples
        mrsa_samples = [s for s, d in samples_data.items() if 'MRSA' in d.get('typing', {}).get('MRSA_Status', '')]
        
        # Build MRSA-specific combinations
        mrsa_mlst_spa = defaultdict(list)
        mrsa_mlst_sccmec = defaultdict(list)
        mrsa_spa_sccmec = defaultdict(list)
        for sample in mrsa_samples:
            data = samples_data[sample]
            mlst = data.get('typing', {}).get('MLST', 'ND')
            spa = data.get('typing', {}).get('spa_Type', 'ND')
            scc = data.get('typing', {}).get('SCCmec_Type', 'ND')
            if mlst!='ND' and spa!='ND':
                mrsa_mlst_spa[f"{mlst} - {spa}"].append(sample)
            if mlst!='ND' and scc!='ND' and scc!='Not Assigned':
                mrsa_mlst_sccmec[f"{mlst} - {scc}"].append(sample)
            if spa!='ND' and scc!='ND' and scc!='Not Assigned':
                mrsa_spa_sccmec[f"{spa} - {scc}"].append(sample)
        
        html = f"""
        <div class="alert-box alert-danger">
            <i class="fas fa-skull-crossbones fa-2x"></i>
            <div>
                <h3>⚠️ MRSA (Methicillin‑Resistant S. aureus) – A Clinical Priority</h3>
                <p>MRSA is resistant to all β‑lactam antibiotics, including penicillins and cephalosporins. It is a major cause of hospital‑ and community‑acquired infections. The presence of <em>mecA</em> or <em>mecC</em> (carried on SCCmec) defines MRSA.</p>
                <p><strong>{mrsa_status_dist.get('MRSA',0)} MRSA samples</strong> identified.</p>
            </div>
        </div>
        
        <h3>📊 MRSA vs MSSA Distribution</h3>
        <div class="scrollable-table"><table id="mrsa-status-table" class="data-table"><thead><tr><th data-sort="string">Status</th><th data-sort="number">Count</th><th data-sort="number">Percentage</th><th data-sort="string">Common STs</th><th data-sort="string">Common SCCmec Types</th></tr></thead><tbody>
        """
        total = sum(mrsa_status_dist.values())
        for status, count in mrsa_status_dist.most_common():
            if status == 'ND': continue
            pct = (count/total)*100
            sts = set()
            sccs = set()
            for s, d in samples_data.items():
                if d.get('typing', {}).get('MRSA_Status') == status:
                    mlst = d.get('typing', {}).get('MLST')
                    scc = d.get('typing', {}).get('SCCmec_Type')
                    if mlst and mlst!='ND':
                        sts.add(mlst)
                    if scc and scc!='ND' and scc!='Not Assigned':
                        sccs.add(scc)
            st_list = ', '.join(sorted(sts)) if sts else 'None'
            scc_list = ', '.join(sorted(sccs)) if sccs else 'None'
            badge = '<span class="badge badge-mrsa">MRSA</span>' if 'MRSA' in status else '<span class="badge badge-mssa">MSSA</span>'
            html += f"<tr><td>{badge}</td><td>{count}</td><td>{pct:.1f}%</td><td>{st_list}</td><td>{scc_list}</td></tr>"
        html += "</tbody></table></div>"
        
        # MRSA ST-spa combos
        html += f"""
        <h3>🔗 MRSA ST–spa Combinations ({len(mrsa_mlst_spa)} combinations)</h3>
        <input type="text" class="search-box" id="search-mrsa-mlst-spa" onkeyup="searchTable('mrsa-mlst-spa-table','search-mrsa-mlst-spa')" placeholder="🔍 Search...">
        <div class="master-scrollable-container"><table id="mrsa-mlst-spa-table" class="data-table"><thead><tr><th data-sort="string">ST-spa Combination</th><th data-sort="number">Count</th><th data-sort="string">Samples</th></tr></thead><tbody>
        """
        for combo, samples in sorted(mrsa_mlst_spa.items(), key=lambda x: len(x[1]), reverse=True):
            sample_tags = ''.join(f'<span class="genome-tag">{s}</span>' for s in samples)
            html += f"<tr><td><strong>{combo}</strong></td><td>{len(samples)}</td><td><div class='genome-list'>{sample_tags}</div></td></tr>"
        html += "</tbody></table></div>"
        
        # MRSA ST-SCCmec combos
        html += f"""
        <h3>🔗 MRSA ST–SCCmec Combinations ({len(mrsa_mlst_sccmec)} combinations)</h3>
        <input type="text" class="search-box" id="search-mrsa-mlst-sccmec" onkeyup="searchTable('mrsa-mlst-sccmec-table','search-mrsa-mlst-sccmec')" placeholder="🔍 Search...">
        <div class="master-scrollable-container"><table id="mrsa-mlst-sccmec-table" class="data-table"><thead><tr><th data-sort="string">ST-SCCmec Combination</th><th data-sort="number">Count</th><th data-sort="string">Samples</th></tr></thead><tbody>
        """
        for combo, samples in sorted(mrsa_mlst_sccmec.items(), key=lambda x: len(x[1]), reverse=True):
            sample_tags = ''.join(f'<span class="genome-tag">{s}</span>' for s in samples)
            html += f"<tr><td><strong>{combo}</strong></td><td>{len(samples)}</td><td><div class='genome-list'>{sample_tags}</div></td></tr>"
        html += "</tbody></table></div>"
        
        # MRSA spa-SCCmec combos
        html += f"""
        <h3>🔗 MRSA spa–SCCmec Combinations ({len(mrsa_spa_sccmec)} combinations)</h3>
        <input type="text" class="search-box" id="search-mrsa-spa-sccmec" onkeyup="searchTable('mrsa-spa-sccmec-table','search-mrsa-spa-sccmec')" placeholder="🔍 Search...">
        <div class="master-scrollable-container"><table id="mrsa-spa-sccmec-table" class="data-table"><thead><tr><th data-sort="string">spa-SCCmec Combination</th><th data-sort="number">Count</th><th data-sort="string">Samples</th></tr></thead><tbody>
        """
        for combo, samples in sorted(mrsa_spa_sccmec.items(), key=lambda x: len(x[1]), reverse=True):
            sample_tags = ''.join(f'<span class="genome-tag">{s}</span>' for s in samples)
            html += f"<tr><td><strong>{combo}</strong></td><td>{len(samples)}</td><td><div class='genome-list'>{sample_tags}</div></td></tr>"
        html += "</tbody></table></div>"
        return html
    
    def _generate_amr_section(self, kwargs: Dict) -> str:
        """Generate comprehensive AMR genes section with filter buttons, genome search, and database summary"""
        gene_centric = kwargs['gene_centric']
        amr_databases = gene_centric.get('amr_databases', {})
        total_samples = len(kwargs.get('samples_data', {}))
        
        html = """
        <div class="alert-box alert-info">
            <i class="fas fa-biohazard fa-2x"></i>
            <div>
                <h3>🧬 Antimicrobial Resistance (AMR) Genes – Gene‑Centric View</h3>
                <p>Each AMR gene is shown with <strong>all genomes</strong> that contain it. Critical S. aureus AMR genes (<em>mecA</em>, <em>mecC</em>, <em>vanA</em>, <em>vanB</em>) are highlighted. Use the filter buttons to focus on specific resistance classes.</p>
                <p>The genome list is <strong>scrollable</strong> scollable, you can see every sample.</p>
            </div>
        </div>
        
        <!-- Gene search -->
        <input type="text" class="search-box" id="search-amr" onkeyup="searchTable('amr-table','search-amr')" placeholder="🔍 Search AMR genes by name or database...">
        
        <!-- Genome highlight search -->
        <input type="text" class="search-box" id="search-amr-genome" onkeyup="highlightGenome('amr-table','search-amr-genome')" placeholder="🔍 Highlight genomes containing specific text (e.g., sample ID)">
        
        <div class="action-buttons">
            <button class="action-btn btn-primary" onclick="exportTableToCSV('amr-table', 'amr_genes.csv')"><i class="fas fa-download"></i> Export All AMR Genes</button>
            <button class="action-btn btn-danger" onclick="document.getElementById('search-amr').value='mecA'; searchTable('amr-table','search-amr')"><i class="fas fa-skull-crossbones"></i> mecA (MRSA)</button>
            <button class="action-btn btn-danger" onclick="document.getElementById('search-amr').value='mecC'; searchTable('amr-table','search-amr')"><i class="fas fa-skull-crossbones"></i> mecC</button>
            <button class="action-btn btn-warning" onclick="document.getElementById('search-amr').value='vanA'; searchTable('amr-table','search-amr')"><i class="fas fa-biohazard"></i> vanA</button>
            <button class="action-btn btn-warning" onclick="document.getElementById('search-amr').value='vanB'; searchTable('amr-table','search-amr')"><i class="fas fa-biohazard"></i> vanB</button>
            <button class="action-btn btn-info" onclick="document.getElementById('search-amr').value='erm'; searchTable('amr-table','search-amr')"><i class="fas fa-pills"></i> erm (Macrolides)</button>
            <button class="action-btn btn-info" onclick="document.getElementById('search-amr').value='msrA'; searchTable('amr-table','search-amr')"><i class="fas fa-pills"></i> msrA</button>
            <button class="action-btn btn-info" onclick="document.getElementById('search-amr').value='mphC'; searchTable('amr-table','search-amr')"><i class="fas fa-pills"></i> mphC</button>
            <button class="action-btn btn-secondary" onclick="document.getElementById('search-amr').value='tet'; searchTable('amr-table','search-amr')"><i class="fas fa-capsules"></i> tet (Tetracycline)</button>
            <button class="action-btn btn-secondary" onclick="document.getElementById('search-amr').value='aac'; searchTable('amr-table','search-amr')"><i class="fas fa-syringe"></i> aac (Aminoglycosides)</button>
            <button class="action-btn btn-secondary" onclick="document.getElementById('search-amr').value='aph'; searchTable('amr-table','search-amr')"><i class="fas fa-syringe"></i> aph</button>
            <button class="action-btn btn-secondary" onclick="document.getElementById('search-amr').value='ant'; searchTable('amr-table','search-amr')"><i class="fas fa-syringe"></i> ant</button>
            <button class="action-btn btn-light" onclick="document.getElementById('search-amr').value='dfr'; searchTable('amr-table','search-amr')"><i class="fas fa-tablets"></i> dfr (Trimethoprim)</button>
            <button class="action-btn btn-light" onclick="document.getElementById('search-amr').value='cat'; searchTable('amr-table','search-amr')"><i class="fas fa-tablets"></i> cat (Chloramphenicol)</button>
            <button class="action-btn btn-light" onclick="document.getElementById('search-amr').value='bla'; searchTable('amr-table','search-amr')"><i class="fas fa-tablets"></i> bla (Beta‑lactamase)</button>
            <button class="action-btn btn-light" onclick="document.getElementById('search-amr').value='pco'; searchTable('amr-table','search-amr')"><i class="fas fa-tablets"></i> pco (Copper)</button>
            <button class="action-btn btn-light" onclick="document.getElementById('search-amr').value=''; searchTable('amr-table','search-amr')"><i class="fas fa-sync"></i> Clear Search</button>
        </div>
        
        <div style="margin: 10px 0 20px 0; background: #f8f9fa; padding: 15px; border-radius: 8px; font-size: 0.9em; border-left: 4px solid #F44336;">
            <strong><i class="fas fa-info-circle"></i> Role of each AMR gene family in S. aureus:</strong><br>
            • <strong>mecA/mecC</strong> – Methicillin resistance (MRSA), confers resistance to all β‑lactam antibiotics.<br>
            • <strong>vanA/vanB</strong> – Vancomycin resistance, last‑line antibiotic for MRSA.<br>
            • <strong>erm</strong> – Macrolide, lincosamide, streptogramin B resistance (MLS<sub>B</sub>).<br>
            • <strong>msrA</strong> – Macrolide efflux pump.<br>
            • <strong>mphC</strong> – Macrolide phosphotransferase.<br>
            • <strong>tet</strong> – Tetracycline resistance (efflux or ribosomal protection).<br>
            • <strong>aac/aph/ant</strong> – Aminoglycoside modifying enzymes (gentamicin, kanamycin, tobramycin).<br>
            • <strong>dfr</strong> – Trimethoprim resistance.<br>
            • <strong>cat</strong> – Chloramphenicol resistance.<br>
        </div>
        
        <h3><i class="fas fa-shield-virus"></i> All AMR Genes Across Databases</h3>
        <div class="master-scrollable-container">
            <table id="amr-table" class="data-table">
                <thead>
                    <tr>
                        <th data-sort="string">Gene</th>
                        <th data-sort="string">Database</th>
                        <th data-sort="number">Count</th>
                        <th data-sort="number">Percentage</th>
                        <th data-sort="string">Genomes (scrollable)</th>
                    </tr>
                </thead>
                <tbody>
        """
        
        # Combine all AMR genes from all databases
        all_amr_genes = []
        for db_name, genes in amr_databases.items():
            for gene_data in genes:
                all_amr_genes.append(gene_data)
        
        # Sort by count descending
        all_amr_genes.sort(key=lambda x: x['count'], reverse=True)
        
        for gene_data in all_amr_genes:
            gene = gene_data['gene']
            database = gene_data['database']
            frequency = gene_data.get('frequency', str(gene_data['count']))
            count = gene_data['count']
            genomes = gene_data.get('genomes', [])
            
            # Calculate percentage string
            percentage = "0%"
            if '(' in frequency:
                # Already formatted like "5 (2.3%)"
                percentage = frequency.split('(')[-1].replace(')', '').strip()
            elif count > 0 and total_samples > 0:
                percentage = f"{(count/total_samples)*100:.1f}%"
            
            # Critical gene indicator (warning icon)
            is_critical = any(crit in gene.lower() for crit in self.data_analyzer.critical_amr_genes)
            gene_display = f"<strong>{gene}</strong>" + (" ⚠️" if is_critical else "")
            
            # Build scrollable genome tags
            genome_tags = ''.join(f'<span class="genome-tag">{g}</span>' for g in genomes)
            
            html += f"""
                <tr>
                    <td>{gene_display}</td>
                    <td>{database}</td>
                    <td><strong>{count}</strong></td>
                    <td>{percentage}</td>
                    <td><div class="genome-list">{genome_tags}</div></td>
                </tr>
            """
        
        html += """
                </tbody>
            </table>
        </div>
        
        <h3 style="margin-top: 30px;"><i class="fas fa-database"></i> AMR Databases Summary</h3>
        <div style="display: grid; grid-template-columns: repeat(auto-fit, minmax(300px, 1fr)); gap: 20px; margin: 20px 0;">
        """
        
        for db_name, genes in amr_databases.items():
            db_display = db_name.upper() if db_name != 'amrfinder' else 'AMRfinder'
            gene_count = len(genes)
            total_count = sum(g['count'] for g in genes) if genes else 0
            top_genes = ', '.join([f"{g['gene']} ({g['count']})" for g in genes[:3]])
            
            html += f"""
            <div class="database-section">
                <h4>{db_display}</h4>
                <p><strong>{gene_count} unique AMR genes</strong> (Total occurrences: {total_count})</p>
                <p>Top genes: {top_genes}...</p>
            </div>
            """
        
        html += """
        </div>
        """
        
        return html
    
    def _generate_virulence_section(self, kwargs: Dict) -> str:
        """Generate comprehensive virulence genes section with filter buttons, genome search, and database summary"""
        gene_centric = kwargs['gene_centric']
        virulence_databases = gene_centric.get('virulence_databases', {})
        total_samples = len(kwargs.get('samples_data', {}))
        
        html = """
        <div class="alert-box alert-info">
            <i class="fas fa-virus fa-2x"></i>
            <div>
                <h3>🧬 Virulence Factors – Gene‑Centric View</h3>
                <p>Each virulence gene is shown with <strong>all genomes</strong> that contain it. Critical S. aureus virulence genes (PVL, TSST‑1, enterotoxins) are highlighted. Use the filter buttons to focus on specific virulence mechanisms.</p>
                <p>The genome list is <strong>scrollable</strong> you can see every sample.</p>
            </div>
        </div>
        
        <!-- Gene search -->
        <input type="text" class="search-box" id="search-vir" onkeyup="searchTable('vir-table','search-vir')" placeholder="🔍 Search virulence genes by name or database...">
        
        <!-- Genome highlight search -->
        <input type="text" class="search-box" id="search-vir-genome" onkeyup="highlightGenome('vir-table','search-vir-genome')" placeholder="🔍 Highlight genomes containing specific text (e.g., sample ID)">
        
        <div class="action-buttons">
            <button class="action-btn btn-primary" onclick="exportTableToCSV('vir-table', 'virulence_genes.csv')"><i class="fas fa-download"></i> Export All Virulence Genes</button>
            <button class="action-btn btn-danger" onclick="document.getElementById('search-vir').value='luk'; searchTable('vir-table','search-vir')"><i class="fas fa-skull"></i> PVL (lukF/S-PV)</button>
            <button class="action-btn btn-danger" onclick="document.getElementById('search-vir').value='tsst'; searchTable('vir-table','search-vir')"><i class="fas fa-biohazard"></i> TSST‑1 (tsst)</button>
            <button class="action-btn btn-warning" onclick="document.getElementById('search-vir').value='se'; searchTable('vir-table','search-vir')"><i class="fas fa-virus"></i> Enterotoxin Genes</button>
            <button class="action-btn btn-warning" onclick="document.getElementById('search-vir').value='cap'; searchTable('vir-table','search-vir')"><i class="fas fa-virus"></i> Capsule</button>
            <button class="action-btn btn-warning" onclick="document.getElementById('search-vir').value='isd'; searchTable('vir-table','search-vir')"><i class="fas fa-virus"></i> Iron acquisition</button>
            <button class="action-btn btn-warning" onclick="document.getElementById('search-vir').value='esa'; searchTable('vir-table','search-vir')"><i class="fas fa-virus"></i> ESAT‑6 secretion (esa)</button>
            <button class="action-btn btn-warning" onclick="document.getElementById('search-vir').value='ess'; searchTable('vir-table','search-vir')"><i class="fas fa-virus"></i> ESAT‑6 secretion (ess)</button>
            <button class="action-btn btn-info" onclick="document.getElementById('search-vir').value='set'; searchTable('vir-table','search-vir')"><i class="fas fa-hand-peace"></i> SET exotoxin‑like proteins</button>
            <button class="action-btn btn-info" onclick="document.getElementById('search-vir').value='esx'; searchTable('vir-table','search-vir')"><i class="fas fa-hand-peace"></i> ESAT‑6 secretion (esx)</button>
            <button class="action-btn btn-secondary" onclick="document.getElementById('search-vir').value='hl'; searchTable('vir-table','search-vir')"><i class="fas fa-tint"></i> Hemolysins</button>
            <button class="action-btn btn-secondary" onclick="document.getElementById('search-vir').value='scn'; searchTable('vir-table','search-vir')"><i class="fas fa-tint"></i> Immune evasion (SCIN)</button>
            <button class="action-btn btn-secondary" onclick="document.getElementById('search-vir').value='ica'; searchTable('vir-table','search-vir')"><i class="fas fa-tint"></i> Biofilm (ica)</button>
            <button class="action-btn btn-secondary" onclick="document.getElementById('search-vir').value='ssp'; searchTable('vir-table','search-vir')"><i class="fas fa-tint"></i> Serine protease</button>
            <button class="action-btn btn-danger" onclick="document.getElementById('search-vir').value='ads'; searchTable('vir-table','search-vir')"><i class="fas fa-biohazard"></i> Adenosine synthase (immune modulation)</button>
            <button class="action-btn btn-danger" onclick="document.getElementById('search-vir').value='eap'; searchTable('vir-table','search-vir')"><i class="fas fa-biohazard"></i> Immune evasion (Eap)</button>
            <button class="action-btn btn-danger" onclick="document.getElementById('search-vir').value='hys'; searchTable('vir-table','search-vir')"><i class="fas fa-biohazard"></i> Hyaluronidase</button>
            <button class="action-btn btn-secondary" onclick="document.getElementById('search-vir').value='clf'; searchTable('vir-table','search-vir')"><i class="fas fa-tint"></i> Adhesins (ClfA/B)</button>
            <button class="action-btn btn-light" onclick="document.getElementById('search-vir').value=''; searchTable('vir-table','search-vir')"><i class="fas fa-sync"></i> Clear Search</button>
        </div>
        
        <div style="margin: 10px 0 20px 0; background: #f8f9fa; padding: 15px; border-radius: 8px; font-size: 0.9em; border-left: 4px solid #E91E63;">
            <strong><i class="fas fa-info-circle"></i> Role of each virulence factor in S. aureus:</strong><br>
            • <strong>PVL (lukF/S-PV)</strong> – Panton‑Valentine leukocidin, causes leukocyte destruction, linked to severe skin/soft tissue infections.<br>
            • <strong>TSST‑1 (tsst)</strong> – Toxic shock syndrome toxin.<br>
            • <strong>Enterotoxins (sea‑see, seg‑seu)</strong> – Superantigens causing food poisoning and toxic shock.<br>
            • <strong>Exfoliative toxins (eta, etb)</strong> – Staphylococcal scalded skin syndrome.<br>
            • <strong>Hemolysins (hla, hlb, hlg, hld)</strong> – Lyse red blood cells and contribute to tissue damage.<br>
            • <strong>Biofilm (icaADBC)</strong> – Polysaccharide intercellular adhesin, critical for device‑related infections.<br>
            • <strong>Immune evasion (scn, eap, ads)</strong> – Inhibit complement, neutrophil chemotaxis, and adenosine signalling.<br>
            • <strong>Adhesins (clfA/B, fnbA/B)</strong> – Bind to fibrinogen and fibronectin, promoting attachment to host tissues.<br>
            • <em>NB: If your favourite button is missing, please reach out to <strong>brownbeckley94@gmail.com</strong>. All suggestions are welcome!</em>
        </div>
        
        <h3><i class="fas fa-virus"></i> All Virulence Genes Across Databases</h3>
        <div class="master-scrollable-container">
            <table id="vir-table" class="data-table">
                <thead>
                    <tr>
                        <th data-sort="string">Gene</th>
                        <th data-sort="string">Database</th>
                        <th data-sort="number">Count</th>
                        <th data-sort="number">Percentage</th>
                        <th data-sort="string">Genomes (scrollable)</th>
                    </tr>
                </thead>
                <tbody>
        """
        
        # Combine all virulence genes from all databases
        all_vir_genes = []
        for db_name, genes in virulence_databases.items():
            for gene_data in genes:
                all_vir_genes.append(gene_data)
        
        # Sort by count descending
        all_vir_genes.sort(key=lambda x: x['count'], reverse=True)
        
        for gene_data in all_vir_genes:
            gene = gene_data['gene']
            database = gene_data['database']
            frequency = gene_data.get('frequency', str(gene_data['count']))
            count = gene_data['count']
            genomes = gene_data.get('genomes', [])
            
            # Calculate percentage string
            percentage = "0%"
            if '(' in frequency:
                percentage = frequency.split('(')[-1].replace(')', '').strip()
            elif count > 0 and total_samples > 0:
                percentage = f"{(count/total_samples)*100:.1f}%"
            
            # Critical gene indicator (warning icon)
            is_critical = any(crit in gene.lower() for crit in self.data_analyzer.critical_virulence_genes)
            gene_display = f"<strong>{gene}</strong>" + (" ⚠️" if is_critical else "")
            
            # Build scrollable genome tags
            genome_tags = ''.join(f'<span class="genome-tag">{g}</span>' for g in genomes)
            
            html += f"""
                <tr>
                    <td>{gene_display}</td>
                    <td>{database}</td>
                    <td><strong>{count}</strong></td>
                    <td>{percentage}</td>
                    <td><div class="genome-list">{genome_tags}</div></td>
                </tr>
            """
        
        html += """
                </tbody>
            </table>
        </div>
        
        <h3 style="margin-top: 30px;"><i class="fas fa-database"></i> Virulence Databases Summary</h3>
        <div style="display: grid; grid-template-columns: repeat(auto-fit, minmax(300px, 1fr)); gap: 20px; margin: 20px 0;">
        """
        
        for db_name, genes in virulence_databases.items():
            db_display = db_name.upper()
            gene_count = len(genes)
            total_count = sum(g['count'] for g in genes) if genes else 0
            top_genes = ', '.join([f"{g['gene']} ({g['count']})" for g in genes[:3]])
            
            html += f"""
            <div class="database-section">
                <h4>{db_display}</h4>
                <p><strong>{gene_count} unique virulence genes</strong> (Total occurrences: {total_count})</p>
                <p>Top genes: {top_genes}...</p>
            </div>
            """
        
        html += """
        </div>
        """
        
        return html
    
    def _generate_bacmet_section(self, kwargs: Dict) -> str:
        """Generate comprehensive BACMET section with many filter buttons and database summary"""
        gene_centric = kwargs['gene_centric']
        bacmet_databases = gene_centric.get('bacmet_databases', {})
        total_samples = len(kwargs.get('samples_data', {}))
        
        if not bacmet_databases or not any(bacmet_databases.values()):
            return """
            <div class="alert-box alert-warning">
                <i class="fas fa-flask fa-2x"></i>
                <div>
                    <h3>No BACMET Data Available</h3>
                    <p>Biocide and heavy metal resistance genes (BACMET2 database) were not detected in the input HTML files.</p>
                    <p><strong>Why include BACMET?</strong> These genes confer resistance to disinfectants (quaternary ammonium compounds, chlorhexidine) and heavy metals (mercury, copper, arsenic). They can co‑select for antibiotic resistance in hospital environments. If you have ABRicate results with the bacmet2 database, please ensure the HTML report is present.</p>
                </div>
            </div>
            """
        
        html = f"""
        <div class="alert-box alert-info">
            <i class="fas fa-flask fa-2x"></i>
            <div>
                <h3>🧪 BACMET: Biocide & Heavy Metal Resistance – Environmental Co‑Selection</h3>
                <p>BACMET2 genes provide resistance to hospital disinfectants (e.g., quaternary ammonium compounds – <em>qac</em> genes) and heavy metals (<em>mer</em> for mercury, <em>ars</em> for arsenic, <em>cop</em> for copper). These markers can co‑select for antibiotic resistance and promote persistence in healthcare settings.</p>
                <p>Each gene is shown with <strong>all genomes</strong> that carry it. Use the filter buttons to explore specific resistance categories.</p>
                <p>The genome list is <strong>scrollable</strong> – all samples are visible.</p>
            </div>
        </div>
        
        <!-- Gene search -->
        <input type="text" class="search-box" id="search-bac" onkeyup="searchTable('bac-table','search-bac')" placeholder="🔍 Search BACMET genes by name...">
        
        <!-- Genome highlight search -->
        <input type="text" class="search-box" id="search-bac-genome" onkeyup="highlightGenome('bac-table','search-bac-genome')" placeholder="🔍 Highlight genomes containing specific text">
        
        <div class="action-buttons">
            <button class="action-btn btn-primary" onclick="exportTableToCSV('bac-table', 'bacmet_genes.csv')"><i class="fas fa-download"></i> Export BACMET Genes</button>
            
            <!-- Biocide / Disinfectant resistance -->
            <button class="action-btn btn-info" onclick="document.getElementById('search-bac').value='qac'; searchTable('bac-table','search-bac')">qac (Quat. ammonium)</button>
            <button class="action-btn btn-info" onclick="document.getElementById('search-bac').value='qacEdelta1'; searchTable('bac-table','search-bac')">qacEdelta1 (Truncated)</button>
            <button class="action-btn btn-info" onclick="document.getElementById('search-bac').value='cep'; searchTable('bac-table','search-bac')">cep (Chlorhexidine)</button>
            <button class="action-btn btn-info" onclick="document.getElementById('search-bac').value='form'; searchTable('bac-table','search-bac')">form (Formaldehyde)</button>
            
            <!-- Heavy metals – Mercury -->
            <button class="action-btn btn-warning" onclick="document.getElementById('search-bac').value='mer'; searchTable('bac-table','search-bac')">mer (Mercury)</button>
            <!-- Heavy metals – Arsenic -->
            <button class="action-btn btn-warning" onclick="document.getElementById('search-bac').value='ars'; searchTable('bac-table','search-bac')">ars (Arsenic)</button>
            <button class="action-btn btn-warning" onclick="document.getElementById('search-bac').value='arsT'; searchTable('bac-table','search-bac')">arsT (Arsenic)</button>
            <!-- Heavy metals – Copper -->
            <button class="action-btn btn-warning" onclick="document.getElementById('search-bac').value='cop'; searchTable('bac-table','search-bac')">cop (Copper)</button>
            <!-- Heavy metals – Silver -->
            <button class="action-btn btn-warning" onclick="document.getElementById('search-bac').value='sil'; searchTable('bac-table','search-bac')">sil (Silver)</button>
            <!-- Heavy metals – Chromate -->
            <button class="action-btn btn-warning" onclick="document.getElementById('search-bac').value='chr'; searchTable('bac-table','search-bac')">chr (Chromate)</button>
            <!-- Heavy metals – Cadmium -->
            <button class="action-btn btn-warning" onclick="document.getElementById('search-bac').value='cad'; searchTable('bac-table','search-bac')">cad (Cadmium)</button>
            <!-- Heavy metals – Zinc -->
            <button class="action-btn btn-warning" onclick="document.getElementById('search-bac').value='znt'; searchTable('bac-table','search-bac')">znt (Zinc)</button>
            <!-- Heavy metals – Cobalt‑Zinc‑Cadmium efflux -->
            <button class="action-btn btn-warning" onclick="document.getElementById('search-bac').value='czc'; searchTable('bac-table','search-bac')">czc (Co‑Zn‑Cd efflux)</button>
            <!-- Heavy metals – Lead -->
            <button class="action-btn btn-warning" onclick="document.getElementById('search-bac').value='pbr'; searchTable('bac-table','search-bac')">pbr (Lead)</button>
            <!-- Heavy metals – Nickel transport -->
            <button class="action-btn btn-warning" onclick="document.getElementById('search-bac').value='nik'; searchTable('bac-table','search-bac')">nik (Nickel)</button>
            <!-- Heavy metals – Magnesium/Cobalt transport -->
            <button class="action-btn btn-warning" onclick="document.getElementById('search-bac').value='cor'; searchTable('bac-table','search-bac')">cor (Mg/Co transport)</button>
            
            <!-- Stress response & efflux -->
            <button class="action-btn btn-secondary" onclick="document.getElementById('search-bac').value='soxR'; searchTable('bac-table','search-bac')">soxR (Oxidative stress)</button>
            <button class="action-btn btn-secondary" onclick="document.getElementById('search-bac').value='cpxR'; searchTable('bac-table','search-bac')">cpxR (Envelope stress)</button>
            <button class="action-btn btn-secondary" onclick="document.getElementById('search-bac').value='baeR'; searchTable('bac-table','search-bac')">baeR (MDR efflux reg.)</button>
            <button class="action-btn btn-secondary" onclick="document.getElementById('search-bac').value='emr'; searchTable('bac-table','search-bac')">emr (MFS efflux)</button>
            <button class="action-btn btn-secondary" onclick="document.getElementById('search-bac').value='sme'; searchTable('bac-table','search-bac')">sme (MFS efflux)</button>
            <button class="action-btn btn-secondary" onclick="document.getElementById('search-bac').value='norA'; searchTable('bac-table','search-bac')">norA (MFS efflux)</button>
            
            <!-- Clear -->
            <button class="action-btn btn-light" onclick="document.getElementById('search-bac').value=''; searchTable('bac-table','search-bac')"><i class="fas fa-sync"></i> Clear Search</button>
        </div>
        
        <div style="margin: 10px 0 20px 0; background: #f8f9fa; padding: 15px; border-radius: 8px; font-size: 0.9em; border-left: 4px solid #FF5722;">
            <strong><i class="fas fa-info-circle"></i> Environmental co‑selection – why these genes matter:</strong><br>
            • <strong>qac family</strong> – Quaternary ammonium compounds (disinfectants used in hospitals).<br>
            • <strong>cep</strong> – Chlorhexidine resistance (antiseptic).<br>
            • <strong>mer/ars/cop/sil</strong> – Resistance to mercury, arsenic, copper, silver – often linked to metal‑based antimicrobials.<br>
            • <strong>czc/cad/znt</strong> – Efflux of zinc, cadmium, cobalt – co‑selection with antibiotic resistance.<br>
            • <strong>soxR/cpxR/baeR</strong> – Stress response regulators that also upregulate multidrug efflux pumps.<br>
            • <strong>emr/sme/norA</strong> – Multidrug efflux pumps that can export both biocides and antibiotics.
        </div>
        
        <h3>📋 BACMET Genes (Biocide / Heavy Metal Resistance)</h3>
        <div class="master-scrollable-container">
            <table id="bac-table" class="data-table">
                <thead>
                    <tr>
                        <th data-sort="string">Gene</th>
                        <th data-sort="string">Database</th>
                        <th data-sort="number">Count</th>
                        <th data-sort="number">Percentage</th>
                        <th data-sort="string">Genomes (scrollable)</th>
                    </tr>
                </thead>
                <tbody>
        """
        
        # Combine all BACMET genes
        all_bac = []
        for db_name, genes in bacmet_databases.items():
            for g in genes:
                all_bac.append(g)
        all_bac.sort(key=lambda x: x['count'], reverse=True)
        
        for g in all_bac:
            gene = g['gene']
            database = g['database']
            count = g['count']
            freq = g.get('frequency', f"{count} ({count/total_samples*100:.1f}%)" if total_samples>0 else "0")
            if '(' not in freq and total_samples>0:
                pct = (count/total_samples)*100
                freq = f"{count} ({pct:.1f}%)"
            genomes = g.get('genomes', [])
            genome_tags = ''.join(f'<span class="genome-tag">{gen}</span>' for gen in genomes)
            html += f"""
                <tr>
                    <td><strong>{gene}</strong></td>
                    <td>{database}</td>
                    <td>{count}</td>
                    <td>{freq}</td>
                    <td><div class='genome-list'>{genome_tags}</div></td>
                </tr>"""
        
        html += """
                </tbody>
            </table>
        </div>
        
        <h3 style="margin-top: 30px;"><i class="fas fa-database"></i> BACMET Databases Summary</h3>
        <div style="display: grid; grid-template-columns: repeat(auto-fit, minmax(300px, 1fr)); gap: 20px; margin: 20px 0;">
        """
        
        for db_name, genes in bacmet_databases.items():
            db_display = db_name.upper()
            gene_count = len(genes)
            total_count = sum(g['count'] for g in genes) if genes else 0
            top_genes = ', '.join([f"{g['gene']} ({g['count']})" for g in genes[:3]])
            html += f"""
            <div class="database-section">
                <h4>{db_display}</h4>
                <p><strong>{gene_count} unique BACMET genes</strong> (Total occurrences: {total_count})</p>
                <p>Top genes: {top_genes}...</p>
            </div>
            """
        
        html += """
        </div>
        """
        return html

    def _generate_plasmids_section(self, kwargs: Dict) -> str:
        """Generate plasmid replicon section with scrollable list"""
        gene_centric = kwargs['gene_centric']
        plasmid_databases = gene_centric.get('plasmid_databases', {})
        total_samples = len(kwargs.get('samples_data', {}))
        
        html = """
        <div class="alert-box alert-info">
            <i class="fas fa-plug fa-2x"></i>
            <div>
                <h3>🧬 Plasmid Replicons – Horizontal Gene Transfer</h3>
                <p>Plasmid replicons indicate the presence of specific plasmid types. Plasmids are major vehicles for spreading AMR and virulence genes (e.g., <em>mecA</em> on SCCmec is not a plasmid, but many resistance genes like <em>erm</em>, <em>tet</em>, and <em>vanA</em> are plasmid‑borne).</p>
                <p>Each replicon is shown with <strong>all genomes</strong> that carry it (scrollable list).</p>
            </div>
        </div>
        
        <input type="text" class="search-box" id="search-plasmid" onkeyup="searchTable('plasmid-table','search-plasmid')" placeholder="🔍 Search plasmids...">
        <input type="text" class="search-box" id="search-plasmid-genome" onkeyup="highlightGenome('plasmid-table','search-plasmid-genome')" placeholder="🔍 Highlight genomes...">
        
        <div class="action-buttons">
            <button class="action-btn btn-primary" onclick="exportTableToCSV('plasmid-table', 'plasmid_replicons.csv')"><i class="fas fa-download"></i> Export</button>
        </div>
        
        <h3>📋 Plasmid Replicons Detected</h3>
        <div class="master-scrollable-container"><table id="plasmid-table" class="data-table"><thead><tr><th data-sort="string">Replicon</th><th data-sort="string">Database</th><th data-sort="number">Count</th><th data-sort="number">Percentage</th><th data-sort="string">Genomes (scrollable)</th></tr></thead><tbody>
        """
        all_plas = []
        for db_name, genes in plasmid_databases.items():
            for g in genes:
                all_plas.append(g)
        all_plas.sort(key=lambda x: x['count'], reverse=True)
        
        for g in all_plas:
            gene = g['gene']
            database = g['database']
            count = g['count']
            freq = g.get('frequency', f"{count} ({count/total_samples*100:.1f}%)" if total_samples>0 else "0")
            if '(' not in freq and total_samples>0:
                pct = (count/total_samples)*100
                freq = f"{count} ({pct:.1f}%)"
            genomes = g.get('genomes', [])
            genome_tags = ''.join([f'<span class="genome-tag">{gen}</span>' for gen in genomes])
            html += f"<tr><td><strong>{gene}</strong></td><td>{database}</td><td>{count}</td><td>{freq}</td><td><div class='genome-list'>{genome_tags}</div></td></tr>"
        html += "</tbody></table></div>"
        return html
    
    def _generate_pattern_discovery_section(self, kwargs: Dict) -> str:
        """Generate pattern discovery including triple typing combinations"""
        patterns = kwargs['patterns']
        gene_centric = kwargs['gene_centric']
        total_samples = len(kwargs.get('samples_data', {}))
        
        # Triple typing combination (ST - spa - SCCmec)
        triple_combos = patterns.get('triple_combinations', {})
        
        html = f"""
        <div class="alert-box alert-info">
            <i class="fas fa-project-diagram fa-2x"></i>
            <div>
                <h3>🔍 Cross‑Genome Pattern Discovery</h3>
                <p>This tab reveals associations between typing results, gene co‑occurrence, and high‑risk combinations.</p>
                <p><strong>Triple Typing Combination (MLST – spa – SCCmec)</strong> – The most informative molecular fingerprint for S. aureus epidemiology. Use this table to track dominant clones.</p>
            </div>
        </div>
        
        <h3>🔗 Triple Typing Combination (ST – spa – SCCmec) – Outbreak Signature</h3>
        <input type="text" class="search-box" id="search-triple" onkeyup="searchTable('triple-table','search-triple')" placeholder="🔍 Search combination...">
        <div class="master-scrollable-container"><table id="triple-table" class="data-table"><thead><tr><th data-sort="string">ST - spa - SCCmec</th><th data-sort="number">Count</th><th data-sort="string">Samples</th></tr></thead><tbody>
        """
        for combo, samples in sorted(triple_combos.items(), key=lambda x: len(x[1]), reverse=True):
            sample_tags = ''.join(f'<span class="genome-tag">{s}</span>' for s in samples)
            html += f"<tr><td><strong>{combo}</strong></td><td>{len(samples)}</td><td><div class='genome-list'>{sample_tags}</div></td></tr>"
        html += "</tbody></table></div>"
        
        # High-risk combinations
        high_risk = patterns.get('high_risk_combinations', [])
        if high_risk:
            html += f"""
            <h3>⚠️ High‑Risk Combinations (Critical AMR + Critical Virulence)</h3>
            <div class="alert-box alert-danger"><i class="fas fa-radiation"></i><div><strong>{len(high_risk)} samples</strong> carry both critical AMR and virulence genes.</div></div>
            <div class="master-scrollable-container"><table id="highrisk-table" class="data-table"><thead><tr><th data-sort="string">Sample</th><th data-sort="string">MLST</th><th data-sort="string">spa</th><th data-sort="string">SCCmec</th><th data-sort="string">Critical AMR</th><th data-sort="string">Critical Virulence</th></tr></thead><tbody>
            """
            for c in high_risk:
                html += f"<tr><td><strong>{c['sample']}</strong></td><td>{c['mlst']}</td><td>{c['spa_type']}</td><td>{c['sccmec_type']}</td><td>{', '.join(c['critical_amr_genes'])}</td><td>{', '.join(c['critical_virulence_genes'])}</td></tr>"
            html += "</tbody></table></div>"
        
        # Gene co-occurrence (top 500)
        cooc = patterns.get('gene_cooccurrence', {})
        if cooc:
            cooc_list = []
            for g1, partners in cooc.items():
                for g2, cnt in partners.items():
                    cooc_list.append((g1, g2, cnt))
            cooc_list.sort(key=lambda x: x[2], reverse=True)
            html += "<h3>📈 Gene Co‑occurrence (Top 500)</h3><div class='master-scrollable-container'><table class='data-table'><thead><tr><th data-sort='string'>Gene 1</th><th data-sort='string'>Gene 2</th><th data-sort='number'>Co‑occurrence Count</th></tr></thead><tbody>"
            for g1, g2, cnt in cooc_list[:500]:
                html += f"<tr><td>{g1}</td><td>{g2}</td><td>{cnt}</td></tr>"
            html += "</tbody></table></div>"
        
        return html
    
    def _generate_aiguide_section(self, kwargs: Dict) -> str:
        """Generate AI assistant guide"""
        return """
        <div class="alert-box alert-info">
            <i class="fas fa-robot fa-2x"></i>
            <div>
                <h3>🤖 AI Assistant Guide – How to Use with ChatGPT, Claude, or Gemini</h3>
                <p>You can upload this HTML report (or its JSON data) to AI tools to ask detailed questions about your S. aureus dataset.</p>
            </div>
        </div>
        
        <div style="margin: 20px 0;">
            <div class="database-section">
                <h4><i class="fas fa-chart-line"></i> Example Questions</h4>
                <ul><li>What are the most common MLST sequence types in this dataset?</li>
                <li>Which spa types are dominant in MRSA vs MSSA?</li>
                <li>How many samples carry mecA? What are their STs and spa types?</li>
                <li>Are there any samples with vancomycin resistance genes (vanA/vanB)?</li>
                <li>Which samples carry the PVL toxin genes (lukF/S-PV)?</li>
                <li>List all samples with TSST-1 (tsst).</li>
                <li>Show me the triple typing combinations (ST-spa-SCCmec) with more than 2 isolates.</li>
                <li>Which resistance genes co-occur most frequently?</li>
                </ul>
            </div>
            <div class="database-section">
                <h4><i class="fas fa-upload"></i> How to Upload Data</h4>
                <p>Save the <strong>staphscope_ultimate_report.json</strong> file and upload it to the AI. Most AI tools accept file uploads, or you can copy‑paste relevant tables.</p>
            </div>
        </div>
        """
    
    def _generate_citation_section(self, kwargs: Dict) -> str:
        """Generate citation tab with colored references (larger font)"""
        return """
        <div class="alert-box alert-info">
            <i class="fas fa-quote-right fa-2x"></i>
            <div>
                <h3>📚 How to Cite StaphScope and Its Dependencies</h3>
                <p>If you use StaphScope in your research, please cite the main tool and the relevant third‑party tools and databases.</p>
            </div>
        </div>

        <style>
        .citation-card {
            margin: 20px 0;
            padding: 15px;
            border-radius: 12px;
            border-left: 5px solid;
            transition: transform 0.2s;
            font-size: 1rem;  /* Increased base font size */
        }
        .citation-card:hover {
            transform: translateX(5px);
        }
        .citation-card h4 {
            margin-bottom: 12px;
            display: flex;
            align-items: center;
            gap: 10px;
            font-size: 1.2rem;  /* Larger heading */
        }
        .citation-card p {
            font-size: 1rem;
            line-height: 1.4;
            margin: 8px 0;
        }
        .citation-card pre {
            background: rgba(0,0,0,0.05);
            padding: 12px;
            border-radius: 8px;
            overflow-x: auto;
            font-size: 0.95rem;  /* Increased from 0.85em */
            font-family: 'Courier New', Courier, monospace;
            margin-top: 12px;
            line-height: 1.4;
        }
        .citation-card a {
            font-size: 1rem;
        }
        /* Tool colors */
        .citation-main { background: #FFF3E0; border-left-color: #FF9800; }
        .citation-mlst { background: #E8F5E9; border-left-color: #4CAF50; }
        .citation-pubmlst { background: #E3F2FD; border-left-color: #2196F3; }
        .citation-abricate { background: #F3E5F5; border-left-color: #9C27B0; }
        .citation-amrfinder { background: #FFEBEE; border-left-color: #F44336; }
        .citation-sccmecfinder { background: #E0F2F1; border-left-color: #009688; }
        .citation-spa { background: #FFF8E1; border-left-color: #FFC107; }
        .citation-card { background: #F1F8E9; border-left-color: #8BC34A; }
        .citation-resfinder { background: #FBE9E7; border-left-color: #FF5722; }
        .citation-vfdb { background: #E8EAF6; border-left-color: #3F51B5; }
        .citation-plasmidfinder { background: #EFEBE9; border-left-color: #795548; }
        .citation-megares { background: #FCE4EC; border-left-color: #E91E63; }
        </style>

        <div class="citation-card citation-main">
            <h4><i class="fas fa-star"></i> Main StaphScope Citation</h4>
            <p>Beckley, B., Amarh, V. (2026). StaphScope: a species‑optimized computational pipeline for rapid and accessible <em>Staphylococcus aureus</em> genotyping and surveillance. <em>BMC Genomics</em>, 27:123.</p>
            <p><strong>DOI</strong>: <a href="https://doi.org/10.1186/s12864-026-12609-x" target="_blank">10.1186/s12864-026-12609-x</a></p>
            <pre>
    @article{beckley2026staphscope,
    title={StaphScope: a species‑optimized computational pipeline for rapid and accessible Staphylococcus aureus genotyping and surveillance},
    author={Beckley, Brown and Amarh, Vincent},
    journal={BMC Genomics},
    volume={27},
    pages={123},
    year={2026},
    doi={10.1186/s12864-026-12609-x}
    }
            </pre>
        </div>

        <div class="citation-card citation-mlst">
            <h4><i class="fas fa-tools"></i> MLST (Torsten Seemann)</h4>
            <p><a href="https://github.com/tseemann/mlst" target="_blank">https://github.com/tseemann/mlst</a></p>
            <pre>
    @software{seemann_mlst_2018,
    author = {Seemann, T.},
    title = {MLST: Scan contig files against traditional PubMLST typing schemes},
    year = {2018},
    publisher = {GitHub},
    url = {https://github.com/tseemann/mlst}
    }
            </pre>
        </div>

        <div class="citation-card citation-pubmlst">
            <h4><i class="fas fa-tools"></i> PubMLST (Jolley et al.)</h4>
            <p><a href="https://pubmlst.org" target="_blank">https://pubmlst.org</a></p>
            <pre>
    @article{jolley_pubmlst_2018,
    author = {Jolley, K. A. and Bray, J. E. and Maiden, M. C. J.},
    title = {Open-access bacterial population genomics: {BIGSdb} software, the {PubMLST.org} website and their applications},
    journal = {Wellcome Open Research},
    volume = {3},
    pages = {124},
    year = {2018},
    doi = {10.12688/wellcomeopenres.14826.1}
    }
            </pre>
        </div>

        <div class="citation-card citation-abricate">
            <h4><i class="fas fa-tools"></i> ABRicate (Torsten Seemann)</h4>
            <p><a href="https://github.com/tseemann/abricate" target="_blank">https://github.com/tseemann/abricate</a></p>
            <pre>
    @software{seemann_abricate_2018,
    author = {Seemann, T.},
    title = {ABRicate: Mass screening of contigs for antimicrobial resistance and virulence genes},
    year = {2018},
    publisher = {GitHub},
    url = {https://github.com/tseemann/abricate}
    }
            </pre>
        </div>

        <div class="citation-card citation-amrfinder">
            <h4><i class="fas fa-tools"></i> AMRFinderPlus (NCBI)</h4>
            <p><a href="https://www.ncbi.nlm.nih.gov/pathogens/amr/" target="_blank">NCBI AMRFinderPlus</a></p>
            <pre>
    @article{feldgarden_amrfinderplus_2019,
    author = {Feldgarden, M. et al.},
    title = {AMRFinderPlus and the Reference Gene Catalog facilitate examination of the genomic links among antimicrobial resistance, stress response, and virulence},
    journal = {Scientific Reports},
    volume = {11},
    pages = {12728},
    year = {2019},
    doi = {10.1038/s41598-021-91456-0}
    }
            </pre>
        </div>

    <div class="citation-card citation-spa">
        <h4><i class="fas fa-tools"></i> spa Typing (Ridom)</h4>
        <p><a href="https://spa.ridom.de" target="_blank">https://spa.ridom.de</a></p>
        <p><strong>Database statistics (May 2026):</strong> 22,638 spa types, 863 repeats, 472,120 total strains, from 181 countries.</p>
        <pre>
    @article{mellmann_spa_typing_2005,
    author = {Mellmann, A. et al.},
    title = {Evidenzbasierte Hygienemassnahmen mittels spa-Typisierung bei MRSA-Häufungen im Krankenhaus},
    journal = {Deutsche Medizinische Wochenschrift},
    volume = {130},
    number = {22},
    pages = {1364-1368},
    year = {2005},
    doi = {10.1055/s-2005-868351}
    }
        </pre>
    </div>
    
        <div class="citation-card citation-sccmecfinder">
            <h4><i class="fas fa-tools"></i> SCCmecFinder (CGE)</h4>
            <p><a href="https://cge.cbs.dtu.dk/services/SCCmecFinder/" target="_blank">CGE SCCmecFinder</a></p>
            <pre>
    @article{kaya_sccmecfinder_2018,
    author = {Kaya, H. et al.},
    title = {SCCmecFinder, a Web-Based Tool for Typing of Staphylococcal Cassette Chromosome mec in Staphylococcus aureus Using Whole-Genome Sequence Data},
    journal = {mSphere},
    volume = {3},
    number = {1},
    pages = {e00612-17},
    year = {2018},
    doi = {10.1128/mSphere.00612-17}
    }
            </pre>
        </div>

        <div class="citation-card citation-spa">
            <h4><i class="fas fa-tools"></i> spa Typing (Ridom)</h4>
            <p><a href="https://spa.ridom.de" target="_blank">https://spa.ridom.de</a></p>
            <pre>
    @article{mellmann_spa_typing_2005,
    author = {Mellmann, A. et al.},
    title = {Evidenzbasierte Hygienemassnahmen mittels spa-Typisierung bei MRSA-Häufungen im Krankenhaus},
    journal = {Deutsche Medizinische Wochenschrift},
    volume = {130},
    number = {22},
    pages = {1364-1368},
    year = {2005},
    doi = {10.1055/s-2005-868351}
    }
            </pre>
        </div>

        <div class="citation-card citation-card">
            <h4><i class="fas fa-database"></i> CARD</h4>
            <p>Comprehensive Antibiotic Resistance Database</p>
            <pre>
    @article{alcock_card_2023,
    author = {Alcock, B. P. et al.},
    title = {CARD 2023: expanded curation, support for machine learning, and resistome prediction at the Comprehensive Antibiotic Resistance Database},
    journal = {Nucleic Acids Research},
    volume = {51},
    number = {D1},
    pages = {D690-D699},
    year = {2023},
    doi = {10.1093/nar/gkac920}
    }
            </pre>
        </div>

        <div class="citation-card citation-resfinder">
            <h4><i class="fas fa-database"></i> ResFinder</h4>
            <pre>
    @article{bortolaia_resfinder_2020,
    author = {Bortolaia, V. et al.},
    title = {ResFinder 4.0 for predictions of phenotypes from genotypes},
    journal = {Journal of Antimicrobial Chemotherapy},
    volume = {75},
    number = {12},
    pages = {3491-3500},
    year = {2020},
    doi = {10.1093/jac/dkaa345}
    }
            </pre>
        </div>

        <div class="citation-card citation-vfdb">
            <h4><i class="fas fa-database"></i> VFDB – Virulence Factor Database</h4>
            <pre>
    @article{chen_vfdb_2016,
    author = {Chen, L. et al.},
    title = {VFDB 2016: hierarchical and refined dataset for big data analysis—10 years on},
    journal = {Nucleic Acids Research},
    volume = {44},
    number = {D1},
    pages = {D694-D697},
    year = {2016},
    doi = {10.1093/nar/gkv1239}
    }
            </pre>
        </div>

        <div class="citation-card citation-plasmidfinder">
            <h4><i class="fas fa-database"></i> PlasmidFinder</h4>
            <pre>
    @article{carattoli_plasmidfinder_2014,
    author = {Carattoli, A. et al.},
    title = {In silico detection and typing of plasmids using PlasmidFinder and plasmid multilocus sequence typing},
    journal = {Antimicrobial Agents and Chemotherapy},
    volume = {58},
    number = {7},
    pages = {3895-3903},
    year = {2014},
    doi = {10.1128/AAC.02412-14}
    }
            </pre>
        </div>

        <div class="citation-card citation-megares">
            <h4><i class="fas fa-database"></i> MEGARes 3.0</h4>
            <pre>
    @article{bonin_megares_2023,
    author = {Bonin, N. et al.},
    title = {MEGARes and AMR++, v3.0: an updated comprehensive database of antimicrobial resistance determinants and an improved software pipeline for classification using high-throughput sequencing},
    journal = {Nucleic Acids Research},
    volume = {51},
    number = {D1},
    pages = {D744-D752},
    year = {2023},
    doi = {10.1093/nar/gkac1047}
    }
            </pre>
        </div>

        <div class="alert-box alert-success">
            <i class="fas fa-hand-peace"></i>
            <div>
                <strong>Suggested acknowledgement:</strong><br>
                “Genomic analysis was performed using StaphScope [Beckley &amp; Amarh, 2026], which integrates MLST [Seemann, 2018] using the Pubmlst database [Jolley et al., 2018], ABRicate [Seemann, 2018], AMRFinderPlus [Feldgarden et al., 2019], and SCCmecFinder [Kaya et al., 2018] for comprehensive <em>S. aureus</em> characterization. Antimicrobial resistance genes were identified using the CARD [Alcock et al., 2023] and ResFinder [Bortolaia et al., 2020] databases.”
            </div>
        </div>
        """
        

    def _generate_calltoaction_section(self, kwargs: Dict) -> str:
        """Generate Call to Action tab with humorous ESKAPE tool lineup"""
        return """
        <div class="alert-box alert-info">
            <i class="fas fa-bullhorn fa-2x"></i>
            <div>
                <h3>🚀 Call to Action – Help Us Fight AMR!</h3>
                <p>If you find this report useful, please <strong>star our GitHub repository</strong> and share it with your colleagues. Your support drives further development (and keeps us caffeinated).</p>
            </div>
        </div>

        <div style="text-align: center; margin: 40px 0;">
            <i class="fas fa-star" style="font-size: 3em; color: #ffc107;"></i>
            <i class="fas fa-star" style="font-size: 2.5em; color: #ffc107; margin-left: 5px;"></i>
            <i class="fas fa-star" style="font-size: 2em; color: #ffc107; margin-left: 5px;"></i>
            <h3>⭐ Star us on GitHub</h3>
            <p><a href="https://github.com/bbeckley-hub/staphscope-typing-tool" target="_blank">https://github.com/bbeckley-hub/staphscope-typing-tool</a></p>
        </div>

        <div class="database-section">
            <h4><i class="fas fa-laugh-beam"></i> The ESKAPE Rogues' Gallery – Our Other Tools</h4>
            <p>We’ve built a whole zoo of pathogen‑specific pipelines – each one more obsessive about your favourite bug. Behold, the <strong>ESKAPE</strong> lineup (because one bug report is never enough):</p>
            <ul style="margin-left: 20px; list-style-type: none;">
                <li><strong><i class="fab fa-github"></i> EnteroMark</strong> – for <em>Enterococcus faecium</em> (VRE surveillance, the "E" in ESKAPE): <a href="https://github.com/bbeckley-hub/EnteroMark" target="_blank">github.com/bbeckley-hub/EnteroMark</a></li>
                <li><strong>🧫 <i class="fab fa-github"></i> StaphScope</strong> – for <em>Staphylococcus aureus</em> (MRSA, the "S" – you're using it right now!): <a href="https://github.com/bbeckley-hub/staphscope-typing-tool" target="_blank">github.com/bbeckley-hub/staphscope-typing-tool</a></li>
                <li><strong>🧬 <i class="fab fa-github"></i> Kleboscope</strong> – for <em>Klebsiella pneumoniae</em> (carbapenemase hunter, the "K"): <a href="https://github.com/bbeckley-hub/Kleboscope" target="_blank">github.com/bbeckley-hub/Kleboscope</a></li>
                <li><strong>🧪 <i class="fab fa-github"></i> Acinetoscope</strong> – for <em>Acinetobacter baumannii</em> (WHO Critical Priority, the "A"): <a href="https://github.com/bbeckley-hub/acinetoscope" target="_blank">github.com/bbeckley-hub/acinetoscope</a></li>
                <li><strong>🧴 <i class="fab fa-github"></i> Pseudoscope</strong> – for the mighty <em>Pseudomonas aeruginosa</em> (biofilm king, the "P"): <a href="https://github.com/bbeckley-hub/pseudoscope" target="_blank">github.com/bbeckley-hub/pseudoscope</a></li>
                <li><strong>🧫 <i class="fab fa-github"></i> Enteroscope</strong> – for <em>Enterobacter cloacae</em> (the forgotten ESKAPE member, the final "E"): <a href="https://github.com/bbeckley-hub/enteroscope" target="_blank">github.com/bbeckley-hub/enteroscope</a></li>
            </ul>
            <p>And because we like to go the extra mile (and bug you), we also have <strong><i class="fab fa-github"></i> EcoliTyper</strong> for <em>E. coli</em> serotyping – because even non‑ESKAPE bugs deserve love: <a href="https://github.com/bbeckley-hub/EcoliTyper" target="_blank">github.com/bbeckley-hub/EcoliTyper</a></p>
            <p>If you like what we’re doing, give them a ⭐ too – it helps us stay motivated (and caffeinated). We put the "fun" in "fungus" (even though most of these are bacteria).</p>
        </div>

        <div class="database-section">
            <h4><i class="fas fa-chalkboard-user"></i> Collaborate & Contribute – Join the Resistance!</h4>
            <p>We welcome bug reports, feature requests, and collaborations. Open an issue on GitHub or email <strong>brownbeckley94@gmail.com</strong>.</p>
            <p><i class="fas fa-globe-africa"></i> <strong>University of Ghana Medical School</strong> – Putting African genomics on the map, one genome at a time.</p>
            <p><i class="fas fa-microbe"></i> <strong>Fun fact:</strong> The ESKAPE pathogens are responsible for the majority of hospital‑acquired infections and are notoriously multi‑drug resistant. Our tools help you fight back – with code, not just antibiotics.</p>
        </div>

        <div class="alert-box alert-success">
            <i class="fas fa-hand-holding-heart"></i>
            <div>
                <strong>Together we can beat AMR – one genome at a time. Now go forth and type!</strong>
            </div>
        </div>
        """

    def _generate_export_section(self, kwargs: Dict) -> str:
        """Generate export section with links to CSV/JSON files"""
        return """
        <div class="alert-box alert-info">
            <i class="fas fa-download fa-2x"></i>
            <div>
                <h3>📥 Export Data – Download Your Results</h3>
                <p>Download all analysis tables in CSV format or the complete JSON data for downstream use.</p>
            </div>
        </div>

        <div style="display: grid; grid-template-columns: repeat(auto-fit, minmax(280px, 1fr)); gap: 20px; margin: 30px 0;">
            <div class="dashboard-card card-export" onclick="exportTableToCSV('samples-table', 'sample_overview.csv')">
                <div style="font-size: 2em; color: var(--export-color);"><i class="fas fa-table"></i></div>
                <div class="card-label">Sample Overview CSV</div>
                <p>All samples with MLST, spa, SCCmec, MRSA status</p>
            </div>
            <div class="dashboard-card card-export" onclick="exportTableToCSV('amr-table', 'amr_genes.csv')">
                <div style="font-size: 2em; color: var(--export-color);"><i class="fas fa-biohazard"></i></div>
                <div class="card-label">AMR Genes CSV</div>
                <p>Gene‑centric AMR table with genome lists</p>
            </div>
            <div class="dashboard-card card-export" onclick="exportTableToCSV('vir-table', 'virulence_genes.csv')">
                <div style="font-size: 2em; color: var(--export-color);"><i class="fas fa-virus"></i></div>
                <div class="card-label">Virulence Genes CSV</div>
                <p>Gene‑centric virulence table</p>
            </div>
            <div class="dashboard-card card-export" onclick="exportTableToCSV('bac-table', 'bacmet_genes.csv')">
                <div style="font-size: 2em; color: var(--export-color);"><i class="fas fa-flask"></i></div>
                <div class="card-label">BACMET Genes CSV</div>
                <p>Biocide & heavy metal resistance genes</p>
            </div>
            <div class="dashboard-card card-export" onclick="exportTableToCSV('plasmid-table', 'plasmid_replicons.csv')">
                <div style="font-size: 2em; color: var(--export-color);"><i class="fas fa-plug"></i></div>
                <div class="card-label">Plasmid Replicons CSV</div>
                <p>Plasmid replicon distribution</p>
            </div>
            <div class="dashboard-card card-export" onclick="exportTableToCSV('qc-table', 'fasta_qc.csv')">
                <div style="font-size: 2em; color: var(--export-color);"><i class="fas fa-chart-line"></i></div>
                <div class="card-label">FASTA QC CSV</div>
                <p>Assembly quality metrics</p>
            </div>
            <div class="dashboard-card card-export" onclick="location.href='staphscope_ultimate_report.json'">
                <div style="font-size: 2em; color: var(--export-color);"><i class="fas fa-file-code"></i></div>
                <div class="card-label">Complete JSON Data</div>
                <p>All analysis data in structured JSON format</p>
            </div>
        </div>

        <div class="alert-box alert-warning">
            <i class="fas fa-save"></i>
            <div>
                <strong>Note:</strong> The JSON file (<code>staphscope_ultimate_report.json</code>) is saved in the same directory as this HTML report. You can use it for custom scripts or upload to AI tools.
            </div>
        </div>
        """
# =============================================================================
# STAPH ULTIMATE REPORTER CLASS
# =============================================================================

class StaphUltimateReporter:
    """MASTER CLASS: Generates ultimate gene-centric reports for S. aureus"""
    
    def __init__(self, input_dir: Path):
        self.input_dir = Path(input_dir)
        self.output_dir = self.input_dir / "STAPHSCOPE_ULTIMATE_REPORTS"
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
        # Initialize components
        self.parser = StaphHTMLParser()
        self.analyzer = StaphDataAnalyzer()
        self.html_generator = StaphHTMLGenerator(self.analyzer)
        
        # Metadata
        self.metadata = {
            "tool_name": "STAPHSCOPE Ultimate S. aureus Reporter",
            "version": "2.1.0",
            "author": "Brown Beckley <brownbeckley94@gmail.com>",
            "affiliation": "University of Ghana Medical School",
            "analysis_date": datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            "input_directory": str(self.input_dir)
        }
    
    def find_html_files(self) -> Dict[str, List[Path]]:
        """Find all StaphScope HTML report files"""
        print("🔍 Searching for StaphScope HTML reports...")
        
        html_files = {
            'comprehensive': [],
            'amrfinder': [],
            'abricate': [],
            'qc': []
        }
        
        for html_file in self.input_dir.glob("**/*.html"):
            filename = html_file.name.lower()
            
            if 'comprehensive' in filename:
                html_files['comprehensive'].append(html_file)
            elif 'staph_amrfinder_summary_report.html' in filename:
                html_files['amrfinder'].append(html_file)
                print(f"    🎯 Found exact AMRfinder file: {html_file.name}")
            elif 'amrfinder' in filename and html_file not in html_files['amrfinder']:
                html_files['amrfinder'].append(html_file)
            elif 'fasta_qc_summary' in filename or 'fasta_qc' in filename or 'qc_summary' in filename:
                html_files['qc'].append(html_file)
                print(f"    🎯 Found FASTA QC file: {html_file.name}")
            elif any(db in filename for db in self.parser.abricate_databases + ['abricate']):
                html_files['abricate'].append(html_file)
        
        # Print findings
        for file_type, files in html_files.items():
            if files:
                if file_type == 'amrfinder' and len(files) > 1:
                    for f in files:
                        if 'staph_amrfinder_summary_report.html' in f.name.lower():
                            print(f"  📁 {file_type.upper()}: {len(files)} files found (Using: {f.name})")
                            break
                    else:
                        print(f"  📁 {file_type.upper()}: {len(files)} files found (Using: {files[0].name})")
                else:
                    print(f"  📁 {file_type.upper()}: {len(files)} files found")
        
        return html_files
    
    def integrate_all_data(self, html_files: Dict[str, List[Path]]) -> Dict[str, Any]:
        """Integrate data from all StaphScope reports"""
        print("\n🔗 Integrating data from all reports...")
        
        integrated_data = {
            'metadata': self.metadata,
            'samples': {},
            'patterns': {},
            'gene_centric': {},
            'qc_data': {}
        }
        
        # Parse FASTA QC report
        if html_files['qc']:
            integrated_data['qc_data'] = self.parser.parse_qc_report(html_files['qc'][0])
            print(f"  ✅ QC data parsed for {len(integrated_data['qc_data'])} samples")
        
        # Parse comprehensive typing report
        typing_data = {}
        if html_files['comprehensive']:
            typing_data = self.parser.parse_comprehensive_report(html_files['comprehensive'][0])
        
        # Parse AMRfinder report
        amr_by_sample, amr_gene_freq = {}, {}
        if html_files['amrfinder']:
            amr_by_sample, amr_gene_freq = self.parser.parse_amrfinder_report(html_files['amrfinder'][0])
        
        # Parse ABRicate databases (including bacmet2)
        abricate_by_sample = defaultdict(dict)
        abricate_gene_freq = {}
        
        for abricate_file in html_files['abricate']:
            db_name, genes_by_sample, gene_freq = self.parser.parse_abricate_report(abricate_file)
            
            if db_name != 'unknown':
                for sample, genes in genes_by_sample.items():
                    if sample not in abricate_by_sample:
                        abricate_by_sample[sample] = {}
                    abricate_by_sample[sample][db_name] = genes
                
                abricate_gene_freq[db_name] = gene_freq
        
        # Combine all samples
        all_samples = set()
        all_samples.update(typing_data.keys())
        all_samples.update(amr_by_sample.keys())
        all_samples.update(abricate_by_sample.keys())
        all_samples.update(integrated_data['qc_data'].keys())
        
        all_samples = sorted(list(all_samples))
        
        if not all_samples:
            print("❌ No samples found in any report!")
            return {}
        
        print(f"\n📊 Found {len(all_samples)} unique samples")
        
        # Integrate data for each sample
        for sample in all_samples:
            # Get typing data
            typing_info = typing_data.get(sample, {
                'MLST': 'ND',
                'spa_Type': 'ND',
                'SCCmec_Type': 'ND',
                'MRSA_Status': 'ND'
            })
            
            # Get AMRfinder data
            amr_info = amr_by_sample.get(sample, {
                'critical_genes': [],
                'high_risk_genes': [],
                'all_genes': []
            })
            
            # Get ABRicate database data
            abricate_info = abricate_by_sample.get(sample, {})
            
            sample_data = {
                'typing': typing_info,
                'amrfinder': amr_info,
                'abricate_databases': abricate_info
            }
            
            integrated_data['samples'][sample] = sample_data
        
        # Store gene frequencies
        integrated_data['gene_frequencies'] = {
            'amrfinder': amr_gene_freq,
            'abricate': abricate_gene_freq
        }
        
        # Process gene-centric data and patterns
        print("\n🧠 Processing gene-centric analysis...")
        integrated_data['gene_centric'] = self.analyzer.create_gene_centric_tables(integrated_data)
        integrated_data['patterns'] = self.analyzer.create_cross_genome_patterns(integrated_data)
        
        return integrated_data
    
    def generate_json_report(self, integrated_data: Dict[str, Any]) -> Path:
        """Generate comprehensive JSON report"""
        print("\n📝 Generating JSON report...")
        
        output_file = self.output_dir / "staphscope_ultimate_report.json"
        
        # Create serializable copy
        def make_serializable(obj):
            if obj is None:
                return None
            elif isinstance(obj, (str, int, float, bool)):
                return obj
            elif isinstance(obj, (list, tuple)):
                return [make_serializable(item) for item in obj]
            elif isinstance(obj, dict):
                return {str(k): make_serializable(v) for k, v in obj.items()}
            elif isinstance(obj, set):
                return [make_serializable(item) for item in obj]
            elif isinstance(obj, (Counter, defaultdict)):
                return {str(k): make_serializable(v) for k, v in obj.items()}
            elif isinstance(obj, Path):
                return str(obj)
            elif hasattr(obj, 'isoformat'):
                return obj.isoformat()
            else:
                try:
                    return str(obj)
                except:
                    return None
        
        serializable_data = make_serializable(integrated_data)
        
        with open(output_file, 'w', encoding='utf-8') as f:
            json.dump(serializable_data, f, indent=2, ensure_ascii=False)
        
        print(f"    ✅ JSON report saved: {output_file}")
        return output_file
    
    def generate_csv_reports(self, integrated_data: Dict[str, Any]):
        """Generate multiple CSV reports for S. aureus"""
        print("\n📊 Generating CSV reports...")
        
        # 1. Sample summary
        samples_data = []
        for sample, data in integrated_data['samples'].items():
            row = {
                'Sample': sample,
                'MLST': data['typing']['MLST'],
                'spa_Type': data['typing']['spa_Type'],
                'SCCmec_Type': data['typing']['SCCmec_Type'],
                'MRSA_Status': data['typing']['MRSA_Status'],
                'Critical_AMR_Genes': ';'.join(data['amrfinder']['critical_genes']),
                'High_Risk_AMR_Genes': ';'.join(data['amrfinder']['high_risk_genes']),
                'Total_AMR_Genes': len(data['amrfinder']['all_genes']),
                'VFDB_Genes': ';'.join(data.get('abricate_databases', {}).get('vfdb', []))
            }
            samples_data.append(row)
        
        df_samples = pd.DataFrame(samples_data)
        samples_file = self.output_dir / "sample_overview.csv"
        df_samples.to_csv(samples_file, index=False)
        
        # 2. AMR genes (gene-centric)
        amr_data = []
        gene_centric = integrated_data.get('gene_centric', {})
        
        for db_name, genes in gene_centric.get('amr_databases', {}).items():
            for gene_info in genes:
                amr_data.append({
                    'Gene': gene_info['gene'],
                    'Database': gene_info['database'],
                    'Count': gene_info['count'],
                    'Frequency': gene_info['frequency'],
                    'Percentage': f"{(gene_info['count']/len(integrated_data['samples']))*100:.1f}%" if len(integrated_data['samples']) > 0 else "0%",
                    'Genomes': ';'.join(gene_info.get('genomes', []))
                })
        
        if amr_data:
            df_amr = pd.DataFrame(amr_data)
            amr_file = self.output_dir / "amr_genes.csv"
            df_amr.to_csv(amr_file, index=False)
        
        # 3. Virulence genes (gene-centric)
        virulence_data = []
        for db_name, genes in gene_centric.get('virulence_databases', {}).items():
            for gene_info in genes:
                virulence_data.append({
                    'Gene': gene_info['gene'],
                    'Database': gene_info['database'],
                    'Count': gene_info['count'],
                    'Frequency': gene_info['frequency'],
                    'Percentage': f"{(gene_info['count']/len(integrated_data['samples']))*100:.1f}%" if len(integrated_data['samples']) > 0 else "0%",
                    'Genomes': ';'.join(gene_info.get('genomes', []))
                })
        
        if virulence_data:
            df_virulence = pd.DataFrame(virulence_data)
            virulence_file = self.output_dir / "virulence_genes.csv"
            df_virulence.to_csv(virulence_file, index=False)
        
        # 4. BACMET genes (gene-centric)
        bacmet_data = []
        for db_name, genes in gene_centric.get('bacmet_databases', {}).items():
            for gene_info in genes:
                bacmet_data.append({
                    'Gene': gene_info['gene'],
                    'Database': gene_info['database'],
                    'Count': gene_info['count'],
                    'Frequency': gene_info['frequency'],
                    'Percentage': f"{(gene_info['count']/len(integrated_data['samples']))*100:.1f}%" if len(integrated_data['samples']) > 0 else "0%",
                    'Genomes': ';'.join(gene_info.get('genomes', []))
                })
        
        if bacmet_data:
            df_bacmet = pd.DataFrame(bacmet_data)
            bacmet_file = self.output_dir / "bacmet_genes.csv"
            df_bacmet.to_csv(bacmet_file, index=False)
        
        # 5. Plasmid replicons (gene-centric)
        plasmid_data = []
        for db_name, genes in gene_centric.get('plasmid_databases', {}).items():
            for gene_info in genes:
                plasmid_data.append({
                    'Plasmid_Replicon': gene_info['gene'],
                    'Database': gene_info['database'],
                    'Count': gene_info['count'],
                    'Frequency': gene_info['frequency'],
                    'Percentage': f"{(gene_info['count']/len(integrated_data['samples']))*100:.1f}%" if len(integrated_data['samples']) > 0 else "0%",
                    'Genomes': ';'.join(gene_info.get('genomes', []))
                })
        
        if plasmid_data:
            df_plasmids = pd.DataFrame(plasmid_data)
            plasmids_file = self.output_dir / "plasmid_replicons.csv"
            df_plasmids.to_csv(plasmids_file, index=False)
        
        # 6. FASTA QC
        qc_data = integrated_data.get('qc_data', {})
        if qc_data:
            qc_rows = []
            for sample, metrics in qc_data.items():
                row = {'Sample': sample}
                row.update(metrics)
                qc_rows.append(row)
            df_qc = pd.DataFrame(qc_rows)
            qc_file = self.output_dir / "fasta_qc.csv"
            df_qc.to_csv(qc_file, index=False)
        
        # 7. Pattern discovery (including triple combinations)
        pattern_data = []
        patterns = integrated_data['patterns']
        
        for mlst, count in patterns.get('mlst_distribution', Counter()).items():
            pattern_data.append({
                'Pattern_Type': 'MLST_Distribution',
                'MLST': mlst,
                'Count': count
            })
        
        for combo, samples in patterns.get('mlst_spa_combinations', {}).items():
            pattern_data.append({
                'Pattern_Type': 'MLST_spa_Combination',
                'Combination': combo,
                'Samples': ';'.join(samples),
                'Count': len(samples)
            })
        
        for combo, samples in patterns.get('triple_combinations', {}).items():
            pattern_data.append({
                'Pattern_Type': 'Triple_Typing_(ST_spa_SCCmec)',
                'Combination': combo,
                'Samples': ';'.join(samples),
                'Count': len(samples)
            })
        
        for combo in patterns.get('high_risk_combinations', []):
            pattern_data.append({
                'Pattern_Type': 'High_Risk_Combination',
                'Sample': combo['sample'],
                'MLST': combo['mlst'],
                'spa_Type': combo['spa_type'],
                'SCCmec_Type': combo['sccmec_type'],
                'Critical_AMR_Genes': ';'.join(combo['critical_amr_genes']),
                'Critical_Virulence_Genes': ';'.join(combo['critical_virulence_genes'])
            })
        
        if pattern_data:
            df_patterns = pd.DataFrame(pattern_data)
            patterns_file = self.output_dir / "pattern_discovery.csv"
            df_patterns.to_csv(patterns_file, index=False)
        
        print(f"    ✅ CSV reports generated: sample_overview.csv, amr_genes.csv, virulence_genes.csv, bacmet_genes.csv, plasmid_replicons.csv, fasta_qc.csv, pattern_discovery.csv")
    
    def run(self):
        """Run the complete analysis for S. aureus"""
        print("=" * 80)
        print("🧬 STAPHSCOPE ULTIMATE S. AUREUS REPORTER v2.1.0")
        print("=" * 80)
        print(f"📁 Input directory: {self.input_dir}")
        
        # Find HTML files
        html_files = self.find_html_files()
        
        if not any(html_files.values()):
            print("❌ No HTML report files found!")
            return False
        
        # Integrate data
        integrated_data = self.integrate_all_data(html_files)
        if not integrated_data:
            return False
        
        # Generate reports
        print("\n" + "=" * 80)
        print("📊 GENERATING ULTIMATE STAPHSCOPE REPORTS")
        print("=" * 80)
        
        # Generate JSON
        json_file = self.generate_json_report(integrated_data)
        
        # Generate CSV
        self.generate_csv_reports(integrated_data)
        
        # Generate HTML
        html_file = self.html_generator.generate_main_report(integrated_data, self.output_dir)
        
        # Print summary
        total_samples = len(integrated_data['samples'])
        patterns = integrated_data['patterns']
        high_risk = len(patterns.get('high_risk_combinations', []))
        gene_centric = integrated_data['gene_centric']
        
        total_amr_genes = sum(len(genes) for genes in gene_centric.get('amr_databases', {}).values())
        total_virulence_genes = sum(len(genes) for genes in gene_centric.get('virulence_databases', {}).values())
        total_bacmet_genes = sum(len(genes) for genes in gene_centric.get('bacmet_databases', {}).values())
        total_plasmids = sum(len(genes) for genes in gene_centric.get('plasmid_databases', {}).values())
        
        # Count MRSA
        mrsa_count = 0
        for sample_data in integrated_data['samples'].values():
            if 'MRSA' in sample_data.get('typing', {}).get('MRSA_Status', ''):
                mrsa_count += 1
        
        print("\n" + "=" * 80)
        print("✅ ULTIMATE ANALYSIS COMPLETE!")
        print("=" * 80)
        print(f"📁 Output directory: {self.output_dir}")
        print(f"📄 Files generated:")
        print(f"   • staphscope_ultimate_report.html (Interactive report)")
        print(f"   • staphscope_ultimate_report.json (Complete data)")
        print(f"   • sample_overview.csv (Sample data)")
        print(f"   • amr_genes.csv (Gene-centric AMR data)")
        print(f"   • virulence_genes.csv (Gene-centric virulence data)")
        print(f"   • bacmet_genes.csv (Biocide/heavy metal genes)")
        print(f"   • plasmid_replicons.csv (Plasmid replicon data)")
        if integrated_data.get('qc_data'):
            print(f"   • fasta_qc.csv (FASTA QC metrics)")
        print(f"   • pattern_discovery.csv (Pattern analysis including triple typing)")
        
        print(f"\n🔬 KEY FEATURES FOR S. AUREUS:")
        print(f"   • Gene-centric tables with scrollable genome lists (no truncation)")
        print(f"   • MRSA focused analysis: {mrsa_count} MRSA samples identified")
        print(f"   • Complete spa typing: {len(patterns.get('spa_type_distribution', {}))} unique spa types")
        print(f"   • SCCmec analysis: {len(patterns.get('sccmec_distribution', {}))} SCCmec types")
        print(f"   • BACMET (biocide & heavy metal) resistance tracking")
        print(f"   • Triple typing combination (ST – spa – SCCmec)")
        print(f"   • FASTA QC metrics integrated")
        print(f"   • Detailed biology explanations for each tab")
        
        print(f"\n📈 ANALYSIS SUMMARY:")
        print(f"   • {total_samples} total S. aureus samples analyzed")
        print(f"   • {mrsa_count} MRSA samples ({mrsa_count/total_samples*100:.1f}%)")
        print(f"   • {total_amr_genes} AMR genes across all databases")
        print(f"   • {total_virulence_genes} virulence genes")
        print(f"   • {total_bacmet_genes} BACMET (biocide/metal) genes")
        print(f"   • {total_plasmids} plasmid replicons")
        print(f"   • {high_risk} high-risk AMR+virulence combinations")
        
        print("\n🎯 Next steps:")
        print("   1. Open staphscope_ultimate_report.html in your browser")
        print("   2. Use AMR, Virulence, and BACMET tabs to see genes with ALL their genomes (scrollable lists)")
        print("   3. Check MRSA analysis tab for MRSA-specific insights")
        print("   4. Use filter buttons to focus on key resistance/virulence genes")
        print("   5. Examine triple typing combinations under Pattern Discovery tab")
        print("   6. Use print buttons on each section header to print specific sections")
        print("   7. Export data using the Export tab or individual CSV buttons")
        
        print("\n" + "=" * 80)
        return True


# =============================================================================
# MAIN FUNCTION
# =============================================================================

def main():
    """Main function for StaphScope Ultimate Reporter"""
    parser = argparse.ArgumentParser(
        description='STAPHSCOPE Ultimate S. aureus Reporter - Gene‑Centric Analysis',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python staphscope_ultimate_reporter.py -i /path/to/staphscope/reports
  
Author: Brown Beckley <brownbeckley94@gmail.com>
Affiliation: University of Ghana Medical School
        """
    )
    
    parser.add_argument('-i', '--input-dir', required=True,
                       help='Directory containing StaphScope HTML report files')
    parser.add_argument('-o', '--output-dir',
                       help='Custom output directory (optional)')
    
    args = parser.parse_args()
    
    input_dir = Path(args.input_dir)
    
    if not input_dir.exists():
        print(f"❌ Input directory not found: {input_dir}")
        sys.exit(1)
    
    # Create and run reporter
    reporter = StaphUltimateReporter(input_dir)
    
    if args.output_dir:
        reporter.output_dir = Path(args.output_dir)
        reporter.output_dir.mkdir(parents=True, exist_ok=True)
    
    success = reporter.run()
    
    if not success:
        print("❌ Report generation failed!")
        sys.exit(1)


if __name__ == "__main__":
    main()
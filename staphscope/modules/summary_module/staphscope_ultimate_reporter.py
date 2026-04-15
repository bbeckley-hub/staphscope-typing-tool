#!/usr/bin/env python3
"""
STAPHSCOPE ULTIMATE REPORTER - GENE-CENTRIC S. AUREUS ANALYSIS
Enhanced with FASTA QC, AI Guide, Filter Buttons, Combination Tables, Sortable Headers
Author: Beckley Brown <brownbeckley94@gmail.com>
Affiliation: University of Ghana Medical School
Version: 2.0.0 - Kleboscope-inspired Ultimate Edition
Date: 2026-04-14
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
                        # try to convert to numeric if it looks like a number
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
        
        # High-priority AMR genes for filter buttons (20)
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
                    else:
                        gene_centric['amr_databases'][db_name] = gene_list
        
        # Create combined gene frequencies for pattern discovery
        all_genes = []
        
        for db_type in ['amr_databases', 'virulence_databases', 'plasmid_databases']:
            for db_name, genes in gene_centric.get(db_type, {}).items():
                for gene_data in genes:
                    all_genes.append(gene_data)
        
        # Sort combined list by count
        all_genes.sort(key=lambda x: x['count'], reverse=True)
        gene_centric['combined_gene_frequencies'] = all_genes
        
        return gene_centric
    
    def create_cross_genome_patterns(self, integrated_data: Dict[str, Any]) -> Dict[str, Any]:
        """Create cross-genome patterns for S. aureus, including combination tables"""
        patterns = {
            'mlst_distribution': Counter(),
            'spa_type_distribution': Counter(),
            'sccmec_distribution': Counter(),
            'mrsa_status_distribution': Counter(),
            'mlst_spa_combinations': defaultdict(list),
            'mlst_sccmec_combinations': defaultdict(list),
            'spa_sccmec_combinations': defaultdict(list),
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
            'samples': '#2196F3',
            'mlst': '#FF9800',
            'spa': '#9C27B0',
            'sccmec': '#009688',
            'mrsa': '#795548',
            'amr': '#F44336',
            'virulence': '#E91E63',
            'plasmids': '#673AB7',
            'patterns': '#FF5722',
            'qc': '#607D8B',
            'aiguide': '#3F51B5',
            'export': '#3F51B5'
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
        
        # CSS Styles - Updated for S. aureus with DEEP GREEN theme
        css = """
        <style>
        :root {
            --summary-color: #4CAF50;
            --samples-color: #2196F3;
            --mlst-color: #FF9800;
            --spa-color: #9C27B0;
            --sccmec-color: #009688;
            --mrsa-color: #795548;
            --amr-color: #F44336;
            --virulence-color: #E91E63;
            --plasmids-color: #673AB7;
            --patterns-color: #FF5722;
            --qc-color: #607D8B;
            --aiguide-color: #3F51B5;
            --export-color: #3F51B5;
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
        .card-samples { border-left-color: var(--samples-color); }
        .card-mlst { border-left-color: var(--mlst-color); }
        .card-spa { border-left-color: var(--spa-color); }
        .card-sccmec { border-left-color: var(--sccmec-color); }
        .card-mrsa { border-left-color: var(--mrsa-color); }
        .card-amr { border-left-color: var(--amr-color); }
        .card-virulence { border-left-color: var(--virulence-color); }
        .card-plasmids { border-left-color: var(--plasmids-color); }
        .card-patterns { border-left-color: var(--patterns-color); }
        .card-qc { border-left-color: var(--qc-color); }
        .card-aiguide { border-left-color: var(--aiguide-color); }
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
            padding: 12px 25px;
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
        .tab-button.samples.active { background: var(--samples-color); }
        .tab-button.mlst.active { background: var(--mlst-color); }
        .tab-button.spa.active { background: var(--spa-color); }
        .tab-button.sccmec.active { background: var(--sccmec-color); }
        .tab-button.mrsa.active { background: var(--mrsa-color); }
        .tab-button.amr.active { background: var(--amr-color); }
        .tab-button.virulence.active { background: var(--virulence-color); }
        .tab-button.plasmids.active { background: var(--plasmids-color); }
        .tab-button.patterns.active { background: var(--patterns-color); }
        .tab-button.qc.active { background: var(--qc-color); }
        .tab-button.aiguide.active { background: var(--aiguide-color); }
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
        .samples-header { border-color: var(--samples-color); }
        .mlst-header { border-color: var(--mlst-color); }
        .spa-header { border-color: var(--spa-color); }
        .sccmec-header { border-color: var(--sccmec-color); }
        .mrsa-header { border-color: var(--mrsa-color); }
        .amr-header { border-color: var(--amr-color); }
        .virulence-header { border-color: var(--virulence-color); }
        .plasmids-header { border-color: var(--plasmids-color); }
        .patterns-header { border-color: var(--patterns-color); }
        .qc-header { border-color: var(--qc-color); }
        .aiguide-header { border-color: var(--aiguide-color); }
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
        
        .genome-list {
            display: flex;
            flex-wrap: wrap;
            gap: 5px;
            margin-top: 5px;
            max-width: none;
        }
        
        .genome-tag {
            background: #e6ffe6;
            color: #006400;
            padding: 3px 10px;
            border-radius: 12px;
            font-size: 0.85em;
            border: 1px solid #b3ffb3;
            white-space: nowrap;
            margin: 2px;
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
                padding: 10px 15px;
                font-size: 0.9em;
            }
            
            .dashboard-grid {
                grid-template-columns: repeat(auto-fit, minmax(200px, 1fr));
            }
            
            .data-table {
                font-size: 0.85em;
            }
            
            body {
                min-width: auto;
                overflow-x: auto;
            }
        }
        </style>
        """
        
        # JavaScript
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
        
        // Search functionality
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
                    <span>Tool: STAPHSCOPE Ultimate v2.0.0</span>
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
                <div class="card-number">{sum(len(genes) for genes in kwargs['gene_centric'].get('virulence_databases', {}).values())}</div>
                <div class="card-label">Virulence Genes</div>
                <i class="fas fa-virus fa-2x" style="color: var(--virulence-color); margin-top: 10px;"></i>
            </div>
            
            <div class="dashboard-card card-patterns" onclick="switchTab('patterns')">
                <div class="card-number">{len(kwargs['patterns'].get('high_risk_combinations', []))}</div>
                <div class="card-label">High-Risk Combos</div>
                <i class="fas fa-project-diagram fa-2x" style="color: var(--patterns-color); margin-top: 10px;"></i>
            </div>
        </div>
        
        <!-- Tab Navigation -->
        <div class="tab-navigation">
            <button class="tab-button summary active" onclick="switchTab('summary')">
                <i class="fas fa-chart-pie"></i> Summary
            </button>
            <button class="tab-button samples" onclick="switchTab('samples')">
                <i class="fas fa-list-alt"></i> Sample Overview
            </button>
            <button class="tab-button mlst" onclick="switchTab('mlst')">
                <i class="fas fa-code-branch"></i> MLST Analysis
            </button>
            <button class="tab-button spa" onclick="switchTab('spa')">
                <i class="fas fa-dna"></i> spa Typing
            </button>
            <button class="tab-button sccmec" onclick="switchTab('sccmec')">
                <i class="fas fa-shield-alt"></i> SCCmec
            </button>
            <button class="tab-button mrsa" onclick="switchTab('mrsa')">
                <i class="fas fa-skull-crossbones"></i> MRSA Analysis
            </button>
            <button class="tab-button amr" onclick="switchTab('amr')">
                <i class="fas fa-biohazard"></i> AMR Genes
            </button>
            <button class="tab-button virulence" onclick="switchTab('virulence')">
                <i class="fas fa-virus"></i> Virulence Genes
            </button>
            <button class="tab-button plasmids" onclick="switchTab('plasmids')">
                <i class="fas fa-plug"></i> Plasmids
            </button>
            <button class="tab-button patterns" onclick="switchTab('patterns')">
                <i class="fas fa-project-diagram"></i> Pattern Discovery
            </button>
            <button class="tab-button qc" onclick="switchTab('qc')">
                <i class="fas fa-chart-line"></i> FASTA QC
            </button>
            <button class="tab-button aiguide" onclick="switchTab('aiguide')">
                <i class="fas fa-robot"></i> AI Guide
            </button>
            <button class="tab-button export" onclick="switchTab('export')">
                <i class="fas fa-download"></i> Export
            </button>
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
        <div id="samples-tab" class="tab-content">
            <h2 class="section-header samples-header">
                <i class="fas fa-list-alt"></i> Complete Sample Overview
                <button class="print-section-btn" onclick="printSection('samples-tab')">
                    <i class="fas fa-print"></i> Print
                </button>
            </h2>
            {self._generate_sample_overview_section(kwargs)}
        </div>
        
        <!-- MLST Analysis Tab -->
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
            <h3>STAPHSCOPE Ultimate S. aureus Reporter v2.0.0</h3>
            <p>University of Ghana Medical School | Brown Beckley <brownbeckley94@gmail.com></p>
            <p>Generated on {kwargs['metadata'].get('analysis_date', 'Unknown')}</p>
             <p>Please give a big STAR on github if you found the integrated report useful.</p>
        </div>
    </div>
</body>
</html>
        """
        
        return html
    
    def _generate_summary_section(self, kwargs: Dict) -> str:
        """Generate summary section for S. aureus"""
        samples_data = kwargs['samples_data']
        patterns = kwargs['patterns']
        gene_centric = kwargs['gene_centric']
        
        total_samples = len(samples_data)
        total_amr_genes = sum(len(genes) for genes in gene_centric.get('amr_databases', {}).values())
        total_virulence_genes = sum(len(genes) for genes in gene_centric.get('virulence_databases', {}).values())
        total_plasmids = sum(len(genes) for genes in gene_centric.get('plasmid_databases', {}).values())
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
                <h3>Analysis Overview</h3>
                <p>This ultimate gene-centric report analyzes <strong>{total_samples}</strong> S. aureus genomes. 
                Instead of listing genes per genome, we show each gene with all genomes that contain it - 
                making it easy to track gene distribution across samples.</p>
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
                    <tr>
                        <td>Total Samples Analyzed</th>
                        <td><strong>{total_samples}</strong></th>
                        <td>Complete genomic analysis with all databases</th>
                    </tr>
                    <tr>
                        <td>MRSA Samples</th>
                        <td><span class="badge badge-mrsa">{mrsa_count}</span></th>
                        <td>Methicillin-resistant S. aureus</th>
                    </tr>
                    <tr>
                        <td>MSSA Samples</th>
                        <td><span class="badge badge-mssa">{mssa_count}</span></th>
                        <td>Methicillin-sensitive S. aureus</th>
                    </tr>
                    <tr>
                        <td>Unique Sequence Types</th>
                        <td><strong>{len(patterns.get('mlst_distribution', {}))}</strong></th>
                        <td>MLST typing results</th>
                    </tr>
                    <tr>
                        <td>Unique spa Types</th>
                        <td><strong>{len(patterns.get('spa_type_distribution', {}))}</strong></th>
                        <td>Protein A gene typing</th>
                    </tr>
                    <tr>
                        <td>Unique SCCmec Types</th>
                        <td><strong>{len(patterns.get('sccmec_distribution', {}))}</strong></th>
                        <td>SCCmec cassette typing</th>
                    </tr>
                    <tr>
                        <td>AMR Genes Identified</th>
                        <td><strong>{total_amr_genes}</strong></th>
                        <td>Across all AMR databases (AMRfinder, CARD, ResFinder, etc.)</th>
                    </tr>
                    <tr>
                        <td>Virulence Genes Identified</th>
                        <td><strong>{total_virulence_genes}</strong></th>
                        <td>From virulence databases (VFDB)</th>
                    </tr>
                    <tr>
                        <td>Plasmid Replicons</th>
                        <td><strong>{total_plasmids}</strong></th>
                        <td>Plasmid replicon types detected</th>
                    </tr>
                    <tr>
                        <td>High-Risk Combinations</th>
                        <td><span class="badge {'badge-critical' if critical_findings > 0 else 'badge-low'}">{critical_findings}</span></th>
                        <td>Samples with both critical AMR and virulence genes</th>
                    </tr>
                </tbody>
            </table>
        </div>
        
        <h3 style="margin-top: 30px;"><i class="fas fa-lightbulb"></i> Report Features</h3>
        <div style="display: grid; grid-template-columns: repeat(auto-fit, minmax(300px, 1fr)); gap: 20px; margin: 20px 0;">
            <div class="database-section">
                <h4><i class="fas fa-gene"></i> Gene-Centric Approach</h4>
                <p>Each gene is shown with ALL genomes that contain it. No more searching through sample lists!</p>
            </div>
            
            <div class="database-section">
                <h4><i class="fas fa-print"></i> Section-Specific Printing</h4>
                <p>Print button on each section header prints only that section. Report any issues by quick mail!</p>
            </div>
            
            <div class="database-section">
                <h4><i class="fas fa-search"></i> Comprehensive Search</h4>
                <p>Search any table by gene name, sample name, ST, spa type, or any field.</p>
            </div>
            
            <div class="database-section">
                <h4><i class="fas fa-bacteria"></i> S. aureus Specific</h4>
                <p>Special focus on MRSA markers (mecA/mecC), SCCmec types, and S. aureus virulence factors.</p>
            </div>
        </div>
        """
        
        return html
    
    def _generate_sample_overview_section(self, kwargs: Dict) -> str:
        """Generate complete sample overview for S. aureus with sortable headers"""
        samples_data = kwargs['samples_data']
        
        html = f"""
        <input type="text" class="search-box" id="search-samples" 
               onkeyup="searchTable('samples-table', 'search-samples')" 
               placeholder="🔍 Search samples by any field...">
        
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
            
            # Highlight MRSA samples
            row_class = 'class="mrsa-highlight"' if 'MRSA' in mrsa_status else ''
            
            # Format status badge
            status_badge = ''
            if 'MRSA' in mrsa_status:
                status_badge = '<span class="badge badge-mrsa">MRSA</span>'
            elif 'MSSA' in mrsa_status:
                status_badge = '<span class="badge badge-mssa">MSSA</span>'
            else:
                status_badge = mrsa_status
            
            html += f"""
                    <tr {row_class}>
                        <td><strong>{sample}</strong></td>
                        <td>{mlst}</td>
                        <td>{spa_type}</td>
                        <td>{sccmec_type}</td>
                        <td>{status_badge}</td>
                        <td>{virulence_count}</td>
                    </tr>
            """
        
        html += """
                </tbody>
            </table>
        </div>
        
        <div class="alert-box alert-info" style="margin-top: 20px;">
            <i class="fas fa-info-circle"></i>
            <div>
                <p><strong>Note:</strong> MRSA samples are highlighted in red. SCCmec type "Not Assigned" indicates no SCCmec cassette detected.</p>
                <p><strong>Sortable columns:</strong> Click on any column header to sort the table.</p>
            </div>
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
            <i class="fas fa-info-circle fa-2x"></i>
            <div>
                <h3>MLST Analysis</h3>
                <p><strong>{len(mlst_dist)} unique sequence types</strong> identified. Each ST is shown with its associated spa types and SCCmec types.</p>
            </div>
        </div>
        
        <h3><i class="fas fa-chart-bar"></i> Sequence Type Distribution</h3>
        <div class="scrollable-table">
            <table id="mlst-table" class="data-table">
                <thead>
                    <tr>
                        <th data-sort="string">Sequence Type</th>
                        <th data-sort="number">Count</th>
                        <th data-sort="number">Percentage</th>
                        <th data-sort="string">Associated spa Types</th>
                        <th data-sort="string">Associated SCCmec Types</th>
                    </tr>
                </thead>
                <tbody>
        """
        
        total = sum(mlst_dist.values())
        for mlst, count in mlst_dist.most_common():
            if mlst == 'ND':
                continue
                
            percentage = (count / total) * 100
            
            # Get associated spa types for this ST
            associated_spa = []
            for combo, samples in mlst_spa_combos.items():
                if f"{mlst} - " in combo:
                    spa_type = combo.split(" - ")[1]
                    if spa_type not in associated_spa:
                        associated_spa.append(spa_type)
            
            # Get associated SCCmec types for this ST
            associated_sccmec = []
            for combo, samples in mlst_sccmec_combos.items():
                if f"{mlst} - " in combo:
                    sccmec_type = combo.split(" - ")[1]
                    if sccmec_type not in associated_sccmec:
                        associated_sccmec.append(sccmec_type)
            
            spa_list = ', '.join(associated_spa) if associated_spa else 'ND'
            sccmec_list = ', '.join(associated_sccmec) if associated_sccmec else 'ND'
            
            html += f"""
                    <tr>
                        <td><strong>{mlst}</strong></td>
                        <td>{count}</td>
                        <td>{percentage:.1f}%</td>
                        <td>{spa_list}</td>
                        <td>{sccmec_list}</td>
                    </tr>
            """
        
        html += """
                </tbody>
            </table>
        </div>
        
        <h3 style="margin-top: 30px;"><i class="fas fa-project-diagram"></i> ST-spa Combinations</h3>
        <input type="text" class="search-box" id="search-mlst-spa" onkeyup="searchTable('mlst-spa-table', 'search-mlst-spa')" placeholder="🔍 Search ST-spa combinations...">
        <div class="master-scrollable-container">
            <table id="mlst-spa-table" class="data-table">
                <thead>
                    <tr>
                        <th data-sort="string">ST-spa Combination</th>
                        <th data-sort="number">Count</th>
                        <th data-sort="string">Samples</th>
                    </tr>
                </thead>
                <tbody>
        """
        
        for combo, samples in sorted(mlst_spa_combos.items(), key=lambda x: len(x[1]), reverse=True):
            sample_list = ', '.join(samples)  # No truncation
            html += f"""
                    <tr>
                        <td><strong>{combo}</strong></td>
                        <td>{len(samples)}</td>
                        <td>{sample_list}</td>
                    </tr>
            """
        
        html += """
                </tbody>
            </table>
        </div>
        
        <h3 style="margin-top: 30px;"><i class="fas fa-project-diagram"></i> ST-SCCmec Combinations</h3>
        <input type="text" class="search-box" id="search-mlst-sccmec" onkeyup="searchTable('mlst-sccmec-table', 'search-mlst-sccmec')" placeholder="🔍 Search ST-SCCmec combinations...">
        <div class="master-scrollable-container">
            <table id="mlst-sccmec-table" class="data-table">
                <thead>
                    <tr>
                        <th data-sort="string">ST-SCCmec Combination</th>
                        <th data-sort="number">Count</th>
                        <th data-sort="string">Samples</th>
                    </tr>
                </thead>
                <tbody>
        """
        
        for combo, samples in sorted(mlst_sccmec_combos.items(), key=lambda x: len(x[1]), reverse=True):
            sample_list = ', '.join(samples)  # No truncation
            html += f"""
                    <tr>
                        <td><strong>{combo}</strong></td>
                        <td>{len(samples)}</td>
                        <td>{sample_list}</td>
                    </tr>
            """
        
        html += """
                </tbody>
            </table>
        </div>
        """
        
        return html
    
    def _generate_spa_section(self, kwargs: Dict) -> str:
        """Generate spa typing analysis section with spa-ST and spa-SCCmec combinations"""
        patterns = kwargs['patterns']
        spa_dist = patterns.get('spa_type_distribution', Counter())
        mlst_spa_combos = patterns.get('mlst_spa_combinations', {})
        spa_sccmec_combos = patterns.get('spa_sccmec_combinations', {})
        
        html = f"""
        <div class="alert-box alert-info">
            <i class="fas fa-info-circle fa-2x"></i>
            <div>
                <h3>spa Typing Analysis</h3>
                <p><strong>{len(spa_dist)} unique spa types</strong> identified across all samples.</p>
            </div>
        </div>
        
        <h3><i class="fas fa-chart-bar"></i> spa Type Distribution</h3>
        <div class="scrollable-table">
            <table id="spa-table" class="data-table">
                <thead>
                    <tr>
                        <th data-sort="string">spa Type</th>
                        <th data-sort="number">Count</th>
                        <th data-sort="number">Percentage</th>
                        <th data-sort="string">Common STs</th>
                    </tr>
                </thead>
                <tbody>
        """
        
        total = sum(spa_dist.values())
        for spa_type, count in spa_dist.most_common():
            if spa_type == 'ND':
                continue
                
            percentage = (count / total) * 100
            
            # Find STs with this spa type
            sts = []
            for combo, samples in mlst_spa_combos.items():
                if spa_type in combo:
                    st = combo.split(" - ")[0]
                    if st not in sts:
                        sts.append(st)
            
            st_list = ', '.join(sts) if sts else 'None'
            
            html += f"""
                    <tr>
                        <td><strong>{spa_type}</strong></td>
                        <td>{count}</td>
                        <td>{percentage:.1f}%</td>
                        <td>{st_list}</td>
                    </tr>
            """
        
        html += """
                </tbody>
            </table>
        </div>
        
        <h3 style="margin-top: 30px;"><i class="fas fa-project-diagram"></i> spa-ST Combinations</h3>
        <input type="text" class="search-box" id="search-spa-mlst" onkeyup="searchTable('spa-mlst-table', 'search-spa-mlst')" placeholder="🔍 Search spa-ST combinations...">
        <div class="master-scrollable-container">
            <table id="spa-mlst-table" class="data-table">
                <thead>
                    <tr>
                        <th data-sort="string">spa-ST Combination</th>
                        <th data-sort="number">Count</th>
                        <th data-sort="string">Samples</th>
                    </tr>
                </thead>
                <tbody>
        """
        
        # Reuse mlst_spa_combos but show as spa-ST
        for combo, samples in sorted(mlst_spa_combos.items(), key=lambda x: len(x[1]), reverse=True):
            # Convert "ST123 - spa123" to "spa123 - ST123" for this view
            parts = combo.split(" - ")
            if len(parts) == 2:
                rev_combo = f"{parts[1]} - {parts[0]}"
            else:
                rev_combo = combo
            sample_list = ', '.join(samples)
            html += f"""
                    <tr>
                        <td><strong>{rev_combo}</strong></td>
                        <td>{len(samples)}</td>
                        <td>{sample_list}</td>
                    </tr>
            """
        
        html += """
                </tbody>
            </table>
        </div>
        
        <h3 style="margin-top: 30px;"><i class="fas fa-project-diagram"></i> spa-SCCmec Combinations</h3>
        <input type="text" class="search-box" id="search-spa-sccmec" onkeyup="searchTable('spa-sccmec-table', 'search-spa-sccmec')" placeholder="🔍 Search spa-SCCmec combinations...">
        <div class="master-scrollable-container">
            <table id="spa-sccmec-table" class="data-table">
                <thead>
                    <tr>
                        <th data-sort="string">spa-SCCmec Combination</th>
                        <th data-sort="number">Count</th>
                        <th data-sort="string">Samples</th>
                    </tr>
                </thead>
                <tbody>
        """
        
        for combo, samples in sorted(spa_sccmec_combos.items(), key=lambda x: len(x[1]), reverse=True):
            sample_list = ', '.join(samples)
            html += f"""
                    <tr>
                        <td><strong>{combo}</strong></td>
                        <td>{len(samples)}</td>
                        <td>{sample_list}</td>
                    </tr>
            """
        
        html += """
                </tbody>
            </table>
        </div>
        """
        
        return html
    
    def _generate_sccmec_section(self, kwargs: Dict) -> str:
        """Generate SCCmec typing analysis section with SCCmec-ST and SCCmec-spa combinations"""
        patterns = kwargs['patterns']
        sccmec_dist = patterns.get('sccmec_distribution', Counter())
        mlst_sccmec_combos = patterns.get('mlst_sccmec_combinations', {})
        spa_sccmec_combos = patterns.get('spa_sccmec_combinations', {})
        
        html = f"""
        <div class="alert-box alert-info">
            <i class="fas fa-info-circle fa-2x"></i>
            <div>
                <h3>SCCmec Typing Analysis</h3>
                <p><strong>{len(sccmec_dist)} unique SCCmec types</strong> identified. SCCmec cassettes are critical for MRSA classification.</p>
            </div>
        </div>
        
        <h3><i class="fas fa-chart-bar"></i> SCCmec Type Distribution</h3>
        <div class="scrollable-table">
            <table id="sccmec-table" class="data-table">
                <thead>
                    <tr>
                        <th data-sort="string">SCCmec Type</th>
                        <th data-sort="number">Count</th>
                        <th data-sort="number">Percentage</th>
                        <th data-sort="string">Common STs</th>
                        <th data-sort="string">Common spa Types</th>
                    </tr>
                </thead>
                <tbody>
        """
        
        total = sum(sccmec_dist.values())
        for sccmec_type, count in sccmec_dist.most_common():
            if sccmec_type == 'ND' or 'Not Assigned' in sccmec_type:
                continue
                
            percentage = (count / total) * 100
            
            # Find STs with this SCCmec type
            sts = []
            for combo, samples in mlst_sccmec_combos.items():
                if sccmec_type in combo:
                    st = combo.split(" - ")[0]
                    if st not in sts:
                        sts.append(st)
            
            # Find spa types with this SCCmec type
            spas = []
            for combo, samples in spa_sccmec_combos.items():
                if sccmec_type in combo:
                    spa = combo.split(" - ")[0]
                    if spa not in spas:
                        spas.append(spa)
            
            st_list = ', '.join(sts) if sts else 'None'
            spa_list = ', '.join(spas) if spas else 'None'
            
            html += f"""
                    <tr>
                        <td><strong>{sccmec_type}</strong></td>
                        <td>{count}</td>
                        <td>{percentage:.1f}%</td>
                        <td>{st_list}</td>
                        <td>{spa_list}</td>
                    </tr>
            """
        
        html += """
                </tbody>
            </table>
        </div>
        
        <h3 style="margin-top: 30px;"><i class="fas fa-project-diagram"></i> SCCmec-ST Combinations</h3>
        <input type="text" class="search-box" id="search-sccmec-mlst" onkeyup="searchTable('sccmec-mlst-table', 'search-sccmec-mlst')" placeholder="🔍 Search SCCmec-ST combinations...">
        <div class="master-scrollable-container">
            <table id="sccmec-mlst-table" class="data-table">
                <thead>
                    <tr>
                        <th data-sort="string">SCCmec-ST Combination</th>
                        <th data-sort="number">Count</th>
                        <th data-sort="string">Samples</th>
                    </tr>
                </thead>
                <tbody>
        """
        
        # Reuse mlst_sccmec_combos but show as SCCmec-ST
        for combo, samples in sorted(mlst_sccmec_combos.items(), key=lambda x: len(x[1]), reverse=True):
            parts = combo.split(" - ")
            if len(parts) == 2:
                rev_combo = f"{parts[1]} - {parts[0]}"
            else:
                rev_combo = combo
            sample_list = ', '.join(samples)
            html += f"""
                    <tr>
                        <td><strong>{rev_combo}</strong></td>
                        <td>{len(samples)}</td>
                        <td>{sample_list}</td>
                    </tr>
            """
        
        html += """
                </tbody>
            </table>
        </div>
        
        <h3 style="margin-top: 30px;"><i class="fas fa-project-diagram"></i> SCCmec-spa Combinations</h3>
        <input type="text" class="search-box" id="search-sccmec-spa" onkeyup="searchTable('sccmec-spa-table', 'search-sccmec-spa')" placeholder="🔍 Search SCCmec-spa combinations...">
        <div class="master-scrollable-container">
            <table id="sccmec-spa-table" class="data-table">
                <thead>
                    <tr>
                        <th data-sort="string">SCCmec-spa Combination</th>
                        <th data-sort="number">Count</th>
                        <th data-sort="string">Samples</th>
                    </tr>
                </thead>
                <tbody>
        """
        
        # Reuse spa_sccmec_combos but show as SCCmec-spa
        for combo, samples in sorted(spa_sccmec_combos.items(), key=lambda x: len(x[1]), reverse=True):
            parts = combo.split(" - ")
            if len(parts) == 2:
                rev_combo = f"{parts[1]} - {parts[0]}"
            else:
                rev_combo = combo
            sample_list = ', '.join(samples)
            html += f"""
                    <tr>
                        <td><strong>{rev_combo}</strong></td>
                        <td>{len(samples)}</td>
                        <td>{sample_list}</td>
                    </tr>
            """
        
        html += """
                </tbody>
            </table>
        </div>
        """
        
        return html
    
    def _generate_mrsa_section(self, kwargs: Dict) -> str:
        """Generate MRSA analysis section with combination tables filtered by MRSA status"""
        patterns = kwargs['patterns']
        mrsa_status_dist = patterns.get('mrsa_status_distribution', Counter())
        mlst_spa_combos = patterns.get('mlst_spa_combinations', {})
        mlst_sccmec_combos = patterns.get('mlst_sccmec_combinations', {})
        spa_sccmec_combos = patterns.get('spa_sccmec_combinations', {})
        samples_data = kwargs['samples_data']
        
        # Filter only MRSA samples
        mrsa_samples = []
        for sample, data in samples_data.items():
            if 'MRSA' in data.get('typing', {}).get('MRSA_Status', ''):
                mrsa_samples.append(sample)
        
        # Build MRSA-specific combinations
        mrsa_mlst_spa = defaultdict(list)
        mrsa_mlst_sccmec = defaultdict(list)
        mrsa_spa_sccmec = defaultdict(list)
        
        for sample in mrsa_samples:
            data = samples_data[sample]
            mlst = data.get('typing', {}).get('MLST', 'ND')
            spa = data.get('typing', {}).get('spa_Type', 'ND')
            sccmec = data.get('typing', {}).get('SCCmec_Type', 'ND')
            
            if mlst != 'ND' and spa != 'ND':
                mrsa_mlst_spa[f"{mlst} - {spa}"].append(sample)
            if mlst != 'ND' and sccmec != 'ND' and sccmec != 'Not Assigned':
                mrsa_mlst_sccmec[f"{mlst} - {sccmec}"].append(sample)
            if spa != 'ND' and sccmec != 'ND' and sccmec != 'Not Assigned':
                mrsa_spa_sccmec[f"{spa} - {sccmec}"].append(sample)
        
        html = f"""
        <div class="alert-box alert-danger">
            <i class="fas fa-skull-crossbones fa-2x"></i>
            <div>
                <h3>MRSA Analysis</h3>
                <p><strong>{mrsa_status_dist.get('MRSA', 0)} MRSA samples</strong> identified. These carry the mecA or mecC gene conferring methicillin resistance.</p>
            </div>
        </div>
        
        <h3><i class="fas fa-chart-bar"></i> MRSA/MSSA Distribution</h3>
        <div class="scrollable-table">
            <table id="mrsa-status-table" class="data-table">
                <thead>
                    <tr>
                        <th data-sort="string">Status</th>
                        <th data-sort="number">Count</th>
                        <th data-sort="number">Percentage</th>
                        <th data-sort="string">Common STs</th>
                        <th data-sort="string">Common SCCmec Types</th>
                    </tr>
                </thead>
                <tbody>
        """
        
        total = sum(mrsa_status_dist.values())
        for status, count in mrsa_status_dist.most_common():
            if status == 'ND':
                continue
                
            percentage = (count / total) * 100
            
            # Find STs and SCCmec types for this status
            sts = []
            sccmec_types = []
            
            for sample, data in samples_data.items():
                if data.get('typing', {}).get('MRSA_Status') == status:
                    mlst = data.get('typing', {}).get('MLST', '')
                    sccmec = data.get('typing', {}).get('SCCmec_Type', '')
                    
                    if mlst and mlst != 'ND' and mlst not in sts:
                        sts.append(mlst)
                    
                    if sccmec and sccmec != 'ND' and sccmec != 'Not Assigned' and sccmec not in sccmec_types:
                        sccmec_types.append(sccmec)
            
            st_list = ', '.join(sts) if sts else 'None'
            sccmec_list = ', '.join(sccmec_types) if sccmec_types else 'None'
            
            status_badge = '<span class="badge badge-mrsa">MRSA</span>' if 'MRSA' in status else '<span class="badge badge-mssa">MSSA</span>'
            
            html += f"""
                    <tr>
                        <td>{status_badge}</td>
                        <td>{count}</td>
                        <td>{percentage:.1f}%</td>
                        <td>{st_list}</td>
                        <td>{sccmec_list}</td>
                    </tr>
            """
        
        html += f"""
                </tbody>
            </table>
        </div>
        
        <h3 style="margin-top: 30px;"><i class="fas fa-project-diagram"></i> MRSA ST-spa Combinations ({len(mrsa_mlst_spa)} combinations)</h3>
        <input type="text" class="search-box" id="search-mrsa-mlst-spa" onkeyup="searchTable('mrsa-mlst-spa-table', 'search-mrsa-mlst-spa')" placeholder="🔍 Search MRSA ST-spa combinations...">
        <div class="master-scrollable-container">
            <table id="mrsa-mlst-spa-table" class="data-table">
                <thead>
                    <tr>
                        <th data-sort="string">ST-spa Combination</th>
                        <th data-sort="number">Count</th>
                        <th data-sort="string">Samples</th>
                    </tr>
                </thead>
                <tbody>
        """
        
        for combo, samples in sorted(mrsa_mlst_spa.items(), key=lambda x: len(x[1]), reverse=True):
            sample_list = ', '.join(samples)
            html += f"""
                    <tr>
                        <td><strong>{combo}</strong></td>
                        <td>{len(samples)}</td>
                        <td>{sample_list}</td>
                    </tr>
            """
        
        html += f"""
                </tbody>
            </table>
        </div>
        
        <h3 style="margin-top: 30px;"><i class="fas fa-project-diagram"></i> MRSA ST-SCCmec Combinations ({len(mrsa_mlst_sccmec)} combinations)</h3>
        <input type="text" class="search-box" id="search-mrsa-mlst-sccmec" onkeyup="searchTable('mrsa-mlst-sccmec-table', 'search-mrsa-mlst-sccmec')" placeholder="🔍 Search MRSA ST-SCCmec combinations...">
        <div class="master-scrollable-container">
            <table id="mrsa-mlst-sccmec-table" class="data-table">
                <thead>
                    <tr>
                        <th data-sort="string">ST-SCCmec Combination</th>
                        <th data-sort="number">Count</th>
                        <th data-sort="string">Samples</th>
                    </tr>
                </thead>
                <tbody>
        """
        
        for combo, samples in sorted(mrsa_mlst_sccmec.items(), key=lambda x: len(x[1]), reverse=True):
            sample_list = ', '.join(samples)
            html += f"""
                    <tr>
                        <td><strong>{combo}</strong></td>
                        <td>{len(samples)}</td>
                        <td>{sample_list}</td>
                    </tr>
            """
        
        html += f"""
                </tbody>
            </table>
        </div>
        
        <h3 style="margin-top: 30px;"><i class="fas fa-project-diagram"></i> MRSA spa-SCCmec Combinations ({len(mrsa_spa_sccmec)} combinations)</h3>
        <input type="text" class="search-box" id="search-mrsa-spa-sccmec" onkeyup="searchTable('mrsa-spa-sccmec-table', 'search-mrsa-spa-sccmec')" placeholder="🔍 Search MRSA spa-SCCmec combinations...">
        <div class="master-scrollable-container">
            <table id="mrsa-spa-sccmec-table" class="data-table">
                <thead>
                    <tr>
                        <th data-sort="string">spa-SCCmec Combination</th>
                        <th data-sort="number">Count</th>
                        <th data-sort="string">Samples</th>
                    </tr>
                </thead>
                <tbody>
        """
        
        for combo, samples in sorted(mrsa_spa_sccmec.items(), key=lambda x: len(x[1]), reverse=True):
            sample_list = ', '.join(samples)
            html += f"""
                    <tr>
                        <td><strong>{combo}</strong></td>
                        <td>{len(samples)}</td>
                        <td>{sample_list}</td>
                    </tr>
            """
        
        html += """
                </tbody>
            </table>
        </div>
        """
        
        return html
    
    def _generate_amr_section(self, kwargs: Dict) -> str:
        """Generate AMR genes section with gene-centric approach and filter buttons"""
        gene_centric = kwargs['gene_centric']
        amr_databases = gene_centric.get('amr_databases', {})
        total_samples = len(kwargs.get('samples_data', {}))
        
        html = """
        <div class="alert-box alert-info">
            <i class="fas fa-info-circle fa-2x"></i>
            <div>
                <h3>AMR Gene Analysis</h3>
                <p>Each AMR gene is shown with ALL genomes that contain it. Critical S. aureus AMR genes (mecA, mecC) are highlighted. Click buttons in colors to filter genes.</p>
            </div>
        </div>
        
        <input type="text" class="search-box" id="search-amr" 
            onkeyup="searchTable('amr-table', 'search-amr')" 
            placeholder="🔍 Search AMR genes by name or database...">
        
        <div class="action-buttons">
            <button class="action-btn btn-primary" onclick="exportTableToCSV('amr-table', 'amr_genes.csv')">
                <i class="fas fa-download"></i> Export All AMR Genes
            </button>
            <button class="action-btn btn-danger" onclick="document.getElementById('search-amr').value='mecA'; searchTable('amr-table','search-amr')">
                <i class="fas fa-skull-crossbones"></i> mecA (MRSA)
            </button>
            <button class="action-btn btn-danger" onclick="document.getElementById('search-amr').value='mecC'; searchTable('amr-table','search-amr')">
                <i class="fas fa-skull-crossbones"></i> mecC
            </button>
            <button class="action-btn btn-warning" onclick="document.getElementById('search-amr').value='vanA'; searchTable('amr-table','search-amr')">
                <i class="fas fa-biohazard"></i> vanA
            </button>
            <button class="action-btn btn-warning" onclick="document.getElementById('search-amr').value='vanB'; searchTable('amr-table','search-amr')">
                <i class="fas fa-biohazard"></i> vanB
            </button>
            <button class="action-btn btn-info" onclick="document.getElementById('search-amr').value='erm'; searchTable('amr-table','search-amr')">
                <i class="fas fa-pills"></i> erm (Macrolides)
            </button>
            <button class="action-btn btn-info" onclick="document.getElementById('search-amr').value='msrA'; searchTable('amr-table','search-amr')">
                <i class="fas fa-pills"></i> msrA
            </button>
            <button class="action-btn btn-info" onclick="document.getElementById('search-amr').value='mphC'; searchTable('amr-table','search-amr')">
                <i class="fas fa-pills"></i> mphC
            </button>
            <button class="action-btn btn-secondary" onclick="document.getElementById('search-amr').value='tet'; searchTable('amr-table','search-amr')">
                <i class="fas fa-capsules"></i> tet (Tetracycline)
            </button>
            <button class="action-btn btn-secondary" onclick="document.getElementById('search-amr').value='aac'; searchTable('amr-table','search-amr')">
                <i class="fas fa-syringe"></i> aac (Aminoglycosides)
            </button>
            <button class="action-btn btn-secondary" onclick="document.getElementById('search-amr').value='aph'; searchTable('amr-table','search-amr')">
                <i class="fas fa-syringe"></i> aph
            </button>
            <button class="action-btn btn-secondary" onclick="document.getElementById('search-amr').value='ant'; searchTable('amr-table','search-amr')">
                <i class="fas fa-syringe"></i> ant
            </button>
            <button class="action-btn btn-light" onclick="document.getElementById('search-amr').value='dfr'; searchTable('amr-table','search-amr')">
                <i class="fas fa-tablets"></i> dfr (Trimethoprim)
            </button>
            <button class="action-btn btn-light" onclick="document.getElementById('search-amr').value='cat'; searchTable('amr-table','search-amr')">
                <i class="fas fa-tablets"></i> cat (Chloramphenicol)
            </button>
            <button class="action-btn btn-light" onclick="document.getElementById('search-amr').value='bla'; searchTable('amr-table','search-amr')">
                <i class="fas fa-tablets"></i> bla (Beta-lactamase)
            </button>    
            <button class="action-btn btn-light" onclick="document.getElementById('search-amr').value='pco'; searchTable('amr-table','search-amr')">
                <i class="fas fa-tablets"></i> pco (Copper)
            <button class="action-btn btn-light" onclick="document.getElementById('search-amr').value=''; searchTable('amr-table','search-amr')">
                <i class="fas fa-sync"></i> Clear Search
            </button>
        </div>
        
        <div style="margin: 10px 0 20px 0; background: #f8f9fa; padding: 15px; border-radius: 8px; font-size: 0.9em; border-left: 4px solid #F44336;">
            <strong><i class="fas fa-info-circle"></i> Role of each AMR gene family in S. aureus:</strong><br>
            • <strong>mecA/mecC</strong> – Methicillin resistance (MRSA), confers resistance to all β-lactam antibiotics.<br>
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
                        <th data-sort="string">Genomes</th>
                    </tr>
                </thead>
                <tbody>
        """
        
        # Combine all AMR genes from all databases
        all_amr_genes = []
        for db_name, genes in amr_databases.items():
            for gene_data in genes:
                all_amr_genes.append(gene_data)
        
        # Sort by count
        all_amr_genes.sort(key=lambda x: x['count'], reverse=True)
        
        for gene_data in all_amr_genes:
            gene = gene_data['gene']
            database = gene_data['database']
            frequency = gene_data['frequency']
            count = gene_data['count']
            genomes = gene_data.get('genomes', [])
            
            # Calculate percentage
            percentage = "0%"
            if '(' in frequency:
                percentage = frequency.split('(')[-1].replace(')', '').strip()
            elif count > 0 and total_samples > 0:
                percentage = f"{(count/total_samples)*100:.1f}%"
            
            # Check if critical AMR gene - show warning icon
            is_critical = any(crit_gene in gene.lower() for crit_gene in self.data_analyzer.critical_amr_genes)
            gene_display = f"<strong>{gene}</strong>" + (" ⚠️" if is_critical else "")
            
            # Create genome tags - NO TRUNCATION
            genome_tags = ''.join([f'<span class="genome-tag">{g}</span>' for g in genomes])
            
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
            
            html += f"""
            <div class="database-section">
                <h4>{db_display}</h4>
                <p><strong>{gene_count} unique AMR genes</strong> (Total occurrences: {total_count})</p>
                <p>Top genes: {', '.join([f"{g['gene']} ({g['count']})" for g in genes[:3]])}...</p>
            </div>
            """
        
        html += """
        </div>
        """
        
        return html
    
    def _generate_virulence_section(self, kwargs: Dict) -> str:
        """Generate virulence genes section with filter buttons and info box"""
        gene_centric = kwargs['gene_centric']
        virulence_databases = gene_centric.get('virulence_databases', {})
        total_samples = len(kwargs.get('samples_data', {}))
        
        html = """
        <div class="alert-box alert-info">
            <i class="fas fa-info-circle fa-2x"></i>
            <div>
                <h3>Virulence Gene Analysis</h3>
                <p>Each virulence gene is shown with ALL genomes that contain it. Critical S. aureus virulence genes (PVL, TSST-1, enterotoxins) are highlighted. Click buttons in colors to filter virulence factors.</p>
            </div>
        </div>
        
        <input type="text" class="search-box" id="search-virulence" 
               onkeyup="searchTable('virulence-table', 'search-virulence')" 
               placeholder="🔍 Search virulence genes by name or database...">
        
        <div class="action-buttons">
            <button class="action-btn btn-primary" onclick="exportTableToCSV('virulence-table', 'virulence_genes.csv')">
                <i class="fas fa-download"></i> Export All Virulence Genes
            </button>
            <button class="action-btn btn-danger" onclick="document.getElementById('search-virulence').value='luk'; searchTable('virulence-table','search-virulence')">
                <i class="fas fa-skull"></i> PVL (lukF/S-PV)
            </button>
            <button class="action-btn btn-danger" onclick="document.getElementById('search-virulence').value='tsst'; searchTable('virulence-table','search-virulence')">
                <i class="fas fa-biohazard"></i> TSST-1 (tsst)
            </button>
            <button class="action-btn btn-warning" onclick="document.getElementById('search-virulence').value='se'; searchTable('virulence-table','search-virulence')">
                <i class="fas fa-virus"></i> Enterotoxin Genes
            </button>
            <button class="action-btn btn-warning" onclick="document.getElementById('search-virulence').value='cap'; searchTable('virulence-table','search-virulence')">
                <i class="fas fa-virus"></i> Capsule 
            </button>
            <button class="action-btn btn-warning" onclick="document.getElementById('search-virulence').value='isd'; searchTable('virulence-table','search-virulence')">
                <i class="fas fa-virus"></i> Iron acquisition
            </button>
            <button class="action-btn btn-warning" onclick="document.getElementById('search-virulence').value='esa'; searchTable('virulence-table','search-virulence')">
                <i class="fas fa-virus"></i> ESAT‑6 secretion system
            </button>
            <button class="action-btn btn-warning" onclick="document.getElementById('search-virulence').value='ess'; searchTable('virulence-table','search-virulence')">
                <i class="fas fa-virus"></i> ESAT‑6 secretion system
            </button>
            <button class="action-btn btn-info" onclick="document.getElementById('search-virulence').value='set'; searchTable('virulence-table','search-virulence')">
                <i class="fas fa-hand-peace"></i> SET (exotoxin‑like) proteins
            </button>
            <button class="action-btn btn-info" onclick="document.getElementById('search-virulence').value='esx'; searchTable('virulence-table','search-virulence')">
                <i class="fas fa-hand-peace"></i> ESAT‑6 secretion system
            </button>
            </button>
            <button class="action-btn btn-secondary" onclick="document.getElementById('search-virulence').value='hl'; searchTable('virulence-table','search-virulence')">
                <i class="fas fa-tint"></i> Hemolysin
            </button>
            <button class="action-btn btn-secondary" onclick="document.getElementById('search-virulence').value='scn'; searchTable('virulence-table','search-virulence')">
                <i class="fas fa-tint"></i> Immune evasion (staphylococcal complement inhibitor)
            </button>
            <button class="action-btn btn-secondary" onclick="document.getElementById('search-virulence').value='ica'; searchTable('virulence-table','search-virulence')">
                <i class="fas fa-tint"></i> Biofilm
            </button>
            <button class="action-btn btn-secondary" onclick="document.getElementById('search-virulence').value='ssp'; searchTable('virulence-table','search-virulence')">
                <i class="fas fa-tint"></i> Serine protease
            </button>
            <button class="action-btn btn-danger" onclick="document.getElementById('search-virulence').value='ads'; searchTable('virulence-table','search-virulence')">
                <i class="fas fa-biohazard"></i> Adenosine synthase (immune modulation)
            </button>
            <button class="action-btn btn-danger" onclick="document.getElementById('search-virulence').value='eap'; searchTable('virulence-table','search-virulence')">
                <i class="fas fa-biohazard"></i> Immune evasion (extracellular adherence protein)
            </button>
            <button class="action-btn btn-danger" onclick="document.getElementById('search-virulence').value='hys'; searchTable('virulence-table','search-virulence')">
                <i class="fas fa-biohazard"></i> Hyaluronidase
            </button>
            <button class="action-btn btn-secondary" onclick="document.getElementById('search-virulence').value='clf'; searchTable('virulence-table','search-virulence')">
                <i class="fas fa-tint"></i> Adhesin
            </button>
            <button class="action-btn btn-light" onclick="document.getElementById('search-virulence').value=''; searchTable('virulence-table','search-virulence')">
                <i class="fas fa-sync"></i> Clear Search
            </button>
        </div>
        
        <div style="margin: 10px 0 20px 0; background: #f8f9fa; padding: 15px; border-radius: 8px; font-size: 0.9em; border-left: 4px solid #E91E63;">
            <strong><i class="fas fa-info-circle"></i> Role of each virulence factor in S. aureus:</strong><br>
            • <strong>PVL (lukF/S-PV)</strong> – Panton-Valentine leukocidin, causes leukocyte destruction and necrosis, associated with severe skin and soft tissue infections.<br>
            • <strong>TSST-1 (tsst)</strong> – Toxic shock syndrome toxin, causes toxic shock syndrome.<br>
            • <strong>Enterotoxins (sea-see, seg-seu)</strong> – Superantigens causing food poisoning and toxic shock.<br>
            • <strong>Exfoliative toxins (eta, etb)</strong> – Cause staphylococcal scalded skin syndrome.<br>
            • <strong>Hemolysins (hla, hlb, hlg, hld)</strong> – Lyse red blood cells and contribute to tissue damage.<br>
             • <strong>NB:</strong> Please you can reach out to brownbeckley94@gmail.com if your favourite button is missing. All suggestions are WeLcOmE.<br>
        </div>
        
        <h3><i class="fas fa-virus"></i> All Virulence Genes Across Databases</h3>
        <div class="master-scrollable-container">
            <table id="virulence-table" class="data-table">
                <thead>
                    <tr>
                        <th data-sort="string">Gene</th>
                        <th data-sort="string">Database</th>
                        <th data-sort="number">Count</th>
                        <th data-sort="number">Percentage</th>
                        <th data-sort="string">Genomes</th>
                    </tr>
                </thead>
                <tbody>
        """
        
        # Combine all virulence genes from all databases
        all_virulence_genes = []
        for db_name, genes in virulence_databases.items():
            for gene_data in genes:
                all_virulence_genes.append(gene_data)
        
        # Sort by count
        all_virulence_genes.sort(key=lambda x: x['count'], reverse=True)
        
        for gene_data in all_virulence_genes:
            gene = gene_data['gene']
            database = gene_data['database']
            frequency = gene_data['frequency']
            count = gene_data['count']
            genomes = gene_data.get('genomes', [])
            
            # Calculate percentage
            percentage = "0%"
            if '(' in frequency:
                percentage = frequency.split('(')[-1].replace(')', '').strip()
            elif count > 0 and total_samples > 0:
                percentage = f"{(count/total_samples)*100:.1f}%"
            
            # Check if critical virulence gene
            is_critical = any(crit_gene in gene.lower() for crit_gene in self.data_analyzer.critical_virulence_genes)
            gene_display = f"<strong>{gene}</strong>" + (" ⚠️" if is_critical else "")
            
            # Create genome tags - NO TRUNCATION
            genome_tags = ''.join([f'<span class="genome-tag">{g}</span>' for g in genomes])
            
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
            
            html += f"""
            <div class="database-section">
                <h4>{db_display}</h4>
                <p><strong>{gene_count} unique virulence genes</strong> (Total occurrences: {total_count})</p>
                <p>Top genes: {', '.join([f"{g['gene']} ({g['count']})" for g in genes[:3]])}...</p>
            </div>
            """
        
        html += """
        </div>
        """
        
        return html
    
    def _generate_plasmids_section(self, kwargs: Dict) -> str:
        """Generate plasmid replicon analysis section"""
        gene_centric = kwargs['gene_centric']
        plasmid_databases = gene_centric.get('plasmid_databases', {})
        total_samples = len(kwargs.get('samples_data', {}))
        
        html = """
        <div class="alert-box alert-info">
            <i class="fas fa-info-circle fa-2x"></i>
            <div>
                <h3>Plasmid Replicon Analysis</h3>
                <p>Each plasmid replicon is shown with ALL genomes that contain it. Plasmid replicons indicate horizontal gene transfer potential.</p>
            </div>
        </div>
        
        <input type="text" class="search-box" id="search-plasmids" 
               onkeyup="searchTable('plasmids-table', 'search-plasmids')" 
               placeholder="🔍 Search plasmid replicons...">
        
        <div class="action-buttons">
            <button class="action-btn btn-primary" onclick="exportTableToCSV('plasmids-table', 'plasmid_replicons.csv')">
                <i class="fas fa-download"></i> Export Plasmid Replicons
            </button>
        </div>
        
        <h3><i class="fas fa-plug"></i> Plasmid Replicons Detected</h3>
        <div class="master-scrollable-container">
            <table id="plasmids-table" class="data-table">
                <thead>
                    <tr>
                        <th data-sort="string">Plasmid Replicon</th>
                        <th data-sort="string">Database</th>
                        <th data-sort="number">Count</th>
                        <th data-sort="number">Percentage</th>
                        <th data-sort="string">Genomes</th>
                    </tr>
                </thead>
                <tbody>
        """
        
        # Combine all plasmid replicons
        all_plasmids = []
        for db_name, genes in plasmid_databases.items():
            for gene_data in genes:
                all_plasmids.append(gene_data)
        
        # Sort by count
        all_plasmids.sort(key=lambda x: x['count'], reverse=True)
        
        for gene_data in all_plasmids:
            gene = gene_data['gene']
            database = gene_data['database']
            frequency = gene_data['frequency']
            count = gene_data['count']
            genomes = gene_data.get('genomes', [])
            
            # Calculate percentage
            percentage = "0%"
            if '(' in frequency:
                percentage = frequency.split('(')[-1].replace(')', '').strip()
            elif count > 0 and total_samples > 0:
                percentage = f"{(count/total_samples)*100:.1f}%"
            
            # Create genome tags
            genome_tags = ''.join([f'<span class="genome-tag">{g}</span>' for g in genomes])
            
            html += f"""
                    <tr>
                        <td><strong>{gene}</strong></td>
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
        
        <h3 style="margin-top: 30px;"><i class="fas fa-database"></i> Plasmid Databases Summary</h3>
        <div style="display: grid; grid-template-columns: repeat(auto-fit, minmax(300px, 1fr)); gap: 20px; margin: 20px 0;">
        """
        
        for db_name, genes in plasmid_databases.items():
            db_display = db_name.upper()
            gene_count = len(genes)
            total_count = sum(g['count'] for g in genes) if genes else 0
            
            html += f"""
            <div class="database-section">
                <h4>{db_display}</h4>
                <p><strong>{gene_count} unique plasmid replicons</strong> (Total occurrences: {total_count})</p>
                <p>Top replicons: {', '.join([f"{g['gene']} ({g['count']})" for g in genes[:3]])}...</p>
            </div>
            """
        
        html += """
        </div>
        """
        
        return html
    
    def _generate_pattern_discovery_section(self, kwargs: Dict) -> str:
        """Generate pattern discovery section for S. aureus"""
        patterns = kwargs['patterns']
        gene_centric = kwargs['gene_centric']
        
        # Get critical genes from gene-centric data
        critical_amr_genes = []
        critical_virulence_genes = []
        
        for db_type in ['amr_databases', 'virulence_databases']:
            for db_name, genes in gene_centric.get(db_type, {}).items():
                for gene_data in genes:
                    gene = gene_data['gene'].lower()
                    if any(crit in gene for crit in self.data_analyzer.critical_amr_genes):
                        critical_amr_genes.append(gene_data)
                    if any(crit in gene for crit in self.data_analyzer.critical_virulence_genes):
                        critical_virulence_genes.append(gene_data)
        
        html = """
        <div class="alert-box alert-info">
            <i class="fas fa-info-circle fa-2x"></i>
            <div>
                <h3>Cross-Genome Pattern Discovery</h3>
                <p>Discover associations between genes and identify high-risk combinations across all S. aureus samples.</p>
            </div>
        </div>
        
        <div class="action-buttons">
            <button class="action-btn btn-primary" onclick="exportTableToCSV('high-risk-table', 'high_risk_combinations.csv')">
                <i class="fas fa-download"></i> Export High-Risk Combos
            </button>
        </div>
        """
        
        # High-risk combinations
        high_risk_combinations = patterns.get('high_risk_combinations', [])
        if high_risk_combinations:
            html += f"""
            <h3><i class="fas fa-exclamation-triangle"></i> High-Risk Combinations ({len(high_risk_combinations)})</h3>
            <div class="alert-box alert-danger">
                <i class="fas fa-radiation fa-2x"></i>
                <div>
                    <h3>⚠️ Critical Alert</h3>
                    <p><strong>{len(high_risk_combinations)} samples</strong> contain dangerous combinations of critical AMR and virulence genes.</p>
                </div>
            </div>
            
            <div class="master-scrollable-container">
                <table id="high-risk-table" class="data-table">
                    <thead>
                        <tr>
                            <th data-sort="string">Sample</th>
                            <th data-sort="string">MLST</th>
                            <th data-sort="string">spa Type</th>
                            <th data-sort="string">SCCmec Type</th>
                            <th data-sort="string">Critical AMR Genes</th>
                            <th data-sort="string">Critical Virulence Genes</th>
                        </tr>
                    </thead>
                    <tbody>
            """
            
            for combo in high_risk_combinations:
                amr_genes = ', '.join(combo['critical_amr_genes'])
                vir_genes = ', '.join(combo['critical_virulence_genes'])
                
                html += f"""
                        <tr>
                            <td><strong>{combo['sample']}</strong></td>
                            <td>{combo['mlst']}</td>
                            <td>{combo['spa_type']}</td>
                            <td>{combo['sccmec_type']}</td>
                            <td>{amr_genes}</td>
                            <td>{vir_genes}</td>
                        </tr>
                """
            
            html += """
                    </tbody>
                </table>
            </div>
            """
        
        # Critical genes summary - SHOW ALL, NO LIMIT
        html += f"""
        <h3 style="margin-top: 30px;"><i class="fas fa-skull-crossbones"></i> Critical Genes Summary</h3>
        <div style="display: grid; grid-template-columns: repeat(auto-fit, minmax(400px, 1fr)); gap: 20px; margin: 20px 0;">
            <div class="database-section">
                <h4>Critical AMR Genes</h4>
                <p><strong>{len(critical_amr_genes)} critical AMR genes</strong> identified</p>
                <div style="margin-top: 10px;">
        """
        
        # Show ALL critical AMR genes, no limit
        for gene_data in critical_amr_genes:
            html += f"""
                    <div style="padding: 5px; border-bottom: 1px solid #eee;">
                        <strong>{gene_data['gene']}</strong> - {gene_data['count']} genomes
                    </div>
            """
        
        html += """
                </div>
            </div>
            
            <div class="database-section">
                <h4>Critical Virulence Genes</h4>
                <p><strong>critical virulence genes</strong> identified</p>
                <div style="margin-top: 10px;">
        """
        
        # Show ALL critical virulence genes, no limit
        for gene_data in critical_virulence_genes:
            html += f"""
                    <div style="padding: 5px; border-bottom: 1px solid #eee;">
                        <strong>{gene_data['gene']}</strong> - {gene_data['count']} genomes
                    </div>
            """
        
        html += """
                </div>
            </div>
        </div>
        
        <h3 style="margin-top: 30px;"><i class="fas fa-project-diagram"></i> Gene Co-occurrence</h3>
        <p>Genes that frequently occur together across samples (top 100 associations):</p>
        <div class="master-scrollable-container">
            <table class="data-table">
                <thead>
                    <tr>
                        <th data-sort="string">Gene 1</th>
                        <th data-sort="string">Gene 2</th>
                        <th data-sort="number">Co-occurrence Count</th>
                    </tr>
                </thead>
                <tbody>
        """
        
        # Get top co-occurrences
        gene_cooccurrence = patterns.get('gene_cooccurrence', {})
        cooccurrence_list = []
        for gene1, partners in gene_cooccurrence.items():
            for gene2, count in partners.items():
                cooccurrence_list.append((gene1, gene2, count))
        
        cooccurrence_list.sort(key=lambda x: x[2], reverse=True)
        
        for gene1, gene2, count in cooccurrence_list[:100]:
            html += f"""
                    <tr>
                        <td><strong>{gene1}</strong></td>
                        <td><strong>{gene2}</strong></td>
                        <td>{count}</td>
                    </tr>
            """
        
        html += """
                </tbody>
            </table>
        </div>
        """
        
        return html
    
    def _generate_qc_section(self, kwargs: Dict) -> str:
        """Generate FASTA Quality Control section with sortable headers"""
        qc_data = kwargs.get('integrated_data', {}).get('qc_data', {})
        
        if not qc_data:
            return """
            <div class="alert-box alert-warning">
                <i class="fas fa-exclamation-circle fa-2x"></i>
                <div>
                    <h3>No QC Data Available</h3>
                    <p>The FASTA_QC_summary.html file was not found or could not be parsed.</p>
                </div>
            </div>
            """
        
        # Collect all metric columns from first sample
        all_metrics = set()
        for sample, metrics in qc_data.items():
            all_metrics.update(metrics.keys())
        
        metric_list = sorted(all_metrics)
        
        # Build table
        html = f"""
        <div class="alert-box alert-info">
            <i class="fas fa-chart-line fa-2x"></i>
            <div>
                <h3>FASTA Quality Control Metrics</h3>
                <p>Assembly quality metrics for each genome. Click on any column header to sort.</p>
            </div>
        </div>
        
        <input type="text" class="search-box" id="search-qc" 
               onkeyup="searchTable('qc-table', 'search-qc')" 
               placeholder="🔍 Search by sample name...">
        
        <div class="action-buttons">
            <button class="action-btn btn-primary" onclick="exportTableToCSV('qc-table', 'fasta_qc.csv')">
                <i class="fas fa-download"></i> Export QC Data
            </button>
        </div>
        
        <div class="master-scrollable-container">
            <table id="qc-table" class="data-table">
                <thead>
                    <tr>
                        <th data-sort="string">Sample</th>
        """
        
        for metric in metric_list:
            html += f'                        <th data-sort="number">{metric}</th>\n'
        
        html += """
                    </tr>
                </thead>
                <tbody>
        """
        
        for sample, metrics in sorted(qc_data.items()):
            html += f"""
                    <tr>
                        <td><strong>{sample}</strong></td>
            """
            for metric in metric_list:
                val = metrics.get(metric, 'ND')
                if isinstance(val, float):
                    if val > 1e6:
                        val = f"{val:,.0f}"
                    elif val > 1000:
                        val = f"{val:,.0f}"
                    else:
                        val = f"{val:.2f}"
                html += f"                        <td>{val}</td>\n"
            html += "                    </tr>\n"
        
        html += """
                </tbody>
            </table>
        </div>
        """
        
        return html
    
    def _generate_aiguide_section(self, kwargs: Dict) -> str:
        """Generate AI Assistant Guide tab"""
        return """
        <div class="alert-box alert-info">
            <i class="fas fa-lightbulb fa-2x"></i>
            <div>
                <h3>How to Use AI with This Report</h3>
                <p>You can upload this HTML report (or its JSON data) to AI tools like ChatGPT, Claude, or Gemini to ask detailed questions about your S. aureus dataset.</p>
            </div>
        </div>
        
        <div style="margin: 20px 0;">
            <div class="database-section">
                <h4><i class="fas fa-chart-line"></i> General Questions</h4>
                <ul style="margin-left: 20px;">
                    <li>What are the most common MLST sequence types in this dataset?</li>
                    <li>Which spa types are dominant in MRSA vs MSSA?</li>
                    <li>How many samples carry mecA? What are their STs and spa types?</li>
                    <li>Are there any samples with vancomycin resistance genes (vanA/vanB)?</li>
                    <li>What is the overall assembly quality (N50, GC%)?</li>
                </ul>
            </div>
            
            <div class="database-section">
                <h4><i class="fas fa-virus"></i> Virulence Analysis</h4>
                <ul style="margin-left: 20px;">
                    <li>Which samples carry the PVL toxin genes (lukF/S-PV)?</li>
                    <li>List all samples with TSST-1 (tsst).</li>
                    <li>Which enterotoxins are most prevalent (sea, seb, sec, sed, see)?</li>
                    <li>Are there samples with exfoliative toxins (eta, etb)?</li>
                    <li>Which STs are associated with PVL-positive MRSA?</li>
                </ul>
            </div>
            
            <div class="database-section">
                <h4><i class="fas fa-biohazard"></i> Resistance Patterns</h4>
                <ul style="margin-left: 20px;">
                    <li>What is the prevalence of macrolide resistance genes (erm, msrA, mphC)?</li>
                    <li>Which STs are most associated with SCCmec type II vs IV?</li>
                    <li>Are there samples with multidrug resistance profiles (e.g., mecA + erm + tet)?</li>
                    <li>Which resistance genes co-occur most frequently?</li>
                </ul>
            </div>
            
            <div class="database-section">
                <h4><i class="fas fa-project-diagram"></i> Pattern Discovery</h4>
                <ul style="margin-left: 20px;">
                    <li>Identify samples with high‑risk combinations (MRSA + PVL).</li>
                    <li>Which ST‑spa combinations are associated with CA-MRSA vs HA-MRSA?</li>
                    <li>Show me the distribution of SCCmec types across different STs.</li>
                </ul>
            </div>
            
            <div class="database-section">
                <h4><i class="fas fa-upload"></i> How to Upload Data</h4>
                <p>If you are using a text-based AI, you can copy the relevant tables or the JSON data and paste them into the conversation. For example, you can say: "I'm sending you the sample overview table from the StaphScope report. Can you analyze the MRSA distribution?"</p>
                <p>Alternatively, you can save the <strong>staphscope_ultimate_report.json</strong> file and ask the AI to parse it. Many AI tools accept file uploads.</p>
            </div>
        </div>
        """
    
    def _generate_export_section(self, kwargs: Dict) -> str:
        """Generate export section"""
        return """
        <div class="alert-box alert-info">
            <i class="fas fa-info-circle fa-2x"></i>
            <div>
                <h3>Export Data and Reports</h3>
                <p>Download comprehensive data in various formats for further analysis and reporting.</p>
            </div>
        </div>
        
        <div style="display: grid; grid-template-columns: repeat(auto-fit, minmax(300px, 1fr)); gap: 20px; margin: 30px 0;">
            <div class="dashboard-card card-export" onclick="exportTableToCSV('samples-table', 'sample_overview.csv')">
                <div style="font-size: 2.5em; color: var(--export-color);"><i class="fas fa-file-csv"></i></div>
                <div class="card-label">Sample Overview CSV</div>
                <p style="font-size: 0.9em; margin-top: 10px;">All samples with MLST, spa type, SCCmec, MRSA status</p>
            </div>
            
            <div class="dashboard-card card-export" onclick="exportTableToCSV('amr-table', 'amr_genes.csv')">
                <div style="font-size: 2.5em; color: var(--export-color);"><i class="fas fa-biohazard"></i></div>
                <div class="card-label">AMR Genes CSV</div>
                <p style="font-size: 0.9em; margin-top: 10px;">All AMR genes with genomes and frequencies</p>
            </div>
            
            <div class="dashboard-card card-export" onclick="exportTableToCSV('virulence-table', 'virulence_genes.csv')">
                <div style="font-size: 2.5em; color: var(--export-color);"><i class="fas fa-virus"></i></div>
                <div class="card-label">Virulence Genes CSV</div>
                <p style="font-size: 0.9em; margin-top: 10px;">All virulence genes with genomes and frequencies</p>
            </div>
            
            <div class="dashboard-card card-export" onclick="exportTableToCSV('plasmids-table', 'plasmid_replicons.csv')">
                <div style="font-size: 2.5em; color: var(--export-color);"><i class="fas fa-plug"></i></div>
                <div class="card-label">Plasmid Replicons CSV</div>
                <p style="font-size: 0.9em; margin-top: 10px;">All plasmid replicons with genomes</p>
            </div>
            
            <div class="dashboard-card card-export" onclick="exportTableToCSV('qc-table', 'fasta_qc.csv')">
                <div style="font-size: 2.5em; color: var(--export-color);"><i class="fas fa-chart-line"></i></div>
                <div class="card-label">FASTA QC CSV</div>
                <p style="font-size: 0.9em; margin-top: 10px;">Assembly quality metrics</p>
            </div>
            
            <div class="dashboard-card card-export" onclick="location.href='staphscope_ultimate_report.json'">
                <div style="font-size: 2.5em; color: var(--export-color);"><i class="fas fa-file-code"></i></div>
                <div class="card-label">Complete JSON Data</div>
                <p style="font-size: 0.9em; margin-top: 10px;">All analysis data in structured JSON format</p>
            </div>
        </div>
        
        <h3><i class="fas fa-download"></i> Available Export Files</h3>
        <div class="master-scrollable-container">
            <table class="data-table">
                <thead>
                    <tr>
                        <th>File</th>
                        <th>Description</th>
                        <th>Format</th>
                        <th>Contents</th>
                    </tr>
                </thead>
                <tbody>
                    <tr>
                        <td><strong>staphscope_ultimate_report.html</strong></td>
                        <td>This interactive HTML report</th>
                        <td>HTML</th>
                        <td>Complete analysis with all sections</th>
                    </tr>
                    <tr>
                        <td><strong>staphscope_ultimate_report.json</strong></th>
                        <td>Complete structured data</th>
                        <td>JSON</th>
                        <td>All analysis data for programmatic use</th>
                    </tr>
                    <tr>
                        <td><strong>sample_overview.csv</strong></th>
                        <td>Sample overview data</th>
                        <td>CSV</th>
                        <td>All samples with typing results</th>
                    </tr>
                    <tr>
                        <td><strong>amr_genes.csv</strong></th>
                        <td>AMR gene analysis</th>
                        <td>CSV</th>
                        <td>All AMR genes with genomes and frequencies</th>
                    </tr>
                    <tr>
                        <td><strong>virulence_genes.csv</strong></th>
                        <td>Virulence gene analysis</th>
                        <td>CSV</th>
                        <td>All virulence genes with genomes and frequencies</th>
                    </tr>
                    <tr>
                        <td><strong>plasmid_replicons.csv</strong></th>
                        <td>Plasmid replicon analysis</th>
                        <td>CSV</th>
                        <td>All plasmid replicons with genomes</th>
                    </tr>
                    <tr>
                        <td><strong>fasta_qc.csv</strong></th>
                        <td>FASTA QC metrics</th>
                        <td>CSV</th>
                        <td>Assembly quality metrics</th>
                    </tr>
                    <tr>
                        <td><strong>pattern_discovery.csv</strong></th>
                        <td>Pattern discovery results</th>
                        <td>CSV</th>
                        <td>Cross-genome patterns and associations</th>
                    </tr>
                </tbody>
            </table>
        </div>
        """


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
            "version": "2.0.0",
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
        
        # Parse ABRicate databases
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
        
        # 4. Plasmid replicons (gene-centric)
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
        
        # 5. FASTA QC
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
        
        # 6. Pattern discovery
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
        
        print(f"    ✅ CSV reports generated: sample_overview.csv, amr_genes.csv, virulence_genes.csv, plasmid_replicons.csv, fasta_qc.csv, pattern_discovery.csv")
    
    def run(self):
        """Run the complete analysis for S. aureus"""
        print("=" * 80)
        print("🧬 STAPHSCOPE ULTIMATE S. AUREUS REPORTER v2.0.0")
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
        print(f"   • plasmid_replicons.csv (Plasmid replicon data)")
        if integrated_data.get('qc_data'):
            print(f"   • fasta_qc.csv (FASTA QC metrics)")
        print(f"   • pattern_discovery.csv (Pattern analysis)")
        
        print(f"\n🔬 KEY FEATURES FOR S. AUREUS:")
        print(f"   • Gene-centric approach: Genes shown with all genomes")
        print(f"   • MRSA focused analysis: {mrsa_count} MRSA samples identified")
        print(f"   • Complete spa typing: {len(patterns.get('spa_type_distribution', {}))} unique spa types")
        print(f"   • SCCmec analysis: {len(patterns.get('sccmec_distribution', {}))} SCCmec types")
        print(f"   • NO TRUNCATION: Complete genome lists for each gene")
        print(f"   • FASTA QC metrics integrated")
        print(f"   • Combination tables: ST-spa, ST-SCCmec, spa-SCCmec")
        
        print(f"\n📈 ANALYSIS SUMMARY:")
        print(f"   • {total_samples} total S. aureus samples analyzed")
        print(f"   • {mrsa_count} MRSA samples ({mrsa_count/total_samples*100:.1f}%)")
        print(f"   • {total_amr_genes} AMR genes across all databases")
        print(f"   • {total_virulence_genes} virulence genes")
        print(f"   • {total_plasmids} plasmid replicons")
        print(f"   • {high_risk} high-risk AMR+virulence combinations")
        
        print("\n🎯 Next steps:")
        print("   1. Open staphscope_ultimate_report.html in your browser")
        print("   2. Use AMR and Virulence tabs to see genes with ALL their genomes")
        print("   3. Check MRSA analysis tab for MRSA-specific insights")
        print("   4. Use filter buttons to focus on key resistance/virulence genes")
        print("   5. Explore combination tables under MLST, spa, SCCmec, and MRSA tabs")
        print("   6. Use print buttons on each section header to print specific sections")
        print("   7. Export data using the Export tab or individual CSV buttons")
        
        print("\n" + "=" * 80)
        return True


def main():
    """Main function for StaphScope Ultimate Reporter"""
    parser = argparse.ArgumentParser(
        description='STAPHSCOPE Ultimate S. aureus Reporter - Gene-Centric Analysis',
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
                       help='Custom output directory')
    
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
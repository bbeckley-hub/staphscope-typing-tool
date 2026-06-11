#!/usr/bin/env python3
"""
STAPHSCOPE ULTIMATE REPORTER - GENE-CENTRIC S. AUREUS ANALYSIS
Enhanced with FASTA QC, AI Guide, Filter Buttons, Combination Tables, Sortable Headers
and DYNAMIC GROUPING BY MLST/spa/SCCmec FOR PATTERN DISCOVERY!
Now with MUTATION TAB (AMRfinderPlus point mutations)
Author: Beckley Brown <brownbeckley94@gmail.com>
Affiliation: University of Ghana Medical School
Version: 2.2.0 
Date: 2026-06-10
"""

import os
import sys
import json
import re
import glob
import argparse
import io  
import pandas as pd
import numpy as np
from pathlib import Path
from typing import Dict, List, Set, Tuple, Any, Optional
from datetime import datetime
from collections import defaultdict, Counter
import warnings
warnings.filterwarnings('ignore')

from bs4 import BeautifulSoup

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
        sample = str(sample_id)
        extensions = ['.fna', '.fasta', '.fa', '.gb', '.gbk', '.gbff', '.txt', '.tsv', '.csv']
        for ext in extensions:
            if sample.endswith(ext):
                sample = sample[:-len(ext)]
        if '/' in sample or '\\' in sample:
            sample = Path(sample).name
        return sample.strip()
    
    def parse_html_table(self, html_content: str, table_index: int = 0) -> pd.DataFrame:
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
        print(f"  🧬 Parsing FASTA QC: {file_path.name}")
        try:
            with open(file_path, 'r', encoding='utf-8') as f:
                html = f.read()
            df = self.parse_html_table(html, 0)
            if df.empty:
                return {}
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
            column_mapping = {
                'Sample': 'Sample', 'sample': 'Sample', 'Genome': 'Sample',
                'MLST': 'MLST', 'MLST Type': 'MLST', 'ST': 'MLST',
                'spa Type': 'spa_Type', 'spa': 'spa_Type',
                'SCCmec Type': 'SCCmec_Type', 'SCCmec': 'SCCmec_Type',
                'MRSA/MSSA Status': 'MRSA_Status', 'MRSA Status': 'MRSA_Status', 'Status': 'MRSA_Status'
            }
            df.rename(columns=column_mapping, inplace=True)
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
        print(f"  🧬 Parsing AMRfinder: {file_path.name}")
        try:
            with open(file_path, 'r', encoding='utf-8') as f:
                html_content = f.read()
            soup = BeautifulSoup(html_content, 'html.parser')
            tables = soup.find_all('table')
            print(f"    Found {len(tables)} tables in the HTML")
            genes_by_genome = {}
            gene_frequencies = {}
            genes_by_genome_table = None
            gene_frequency_table = None
            for i, table in enumerate(tables):
                table_text = table.get_text()
                if 'Genome' in table_text and 'Critical Genes' in table_text:
                    genes_by_genome_table = table
                    print(f"    Found 'Genes by Genome' table at index {i}")
                elif 'Gene' in table_text and 'Frequency' in table_text and 'Prevalence' in table_text:
                    gene_frequency_table = table
                    print(f"    Found 'Gene Frequency' table at index {i}")
            if genes_by_genome_table:
                print(f"    Parsing 'Genes by Genome' table...")
                try:
                    table_html = str(genes_by_genome_table)
                    df_genomes = pd.read_html(io.StringIO(table_html))[0]
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
            if gene_frequency_table:
                print(f"    Parsing 'Gene Frequency' table...")
                try:
                    table_html = str(gene_frequency_table)
                    df_genes = pd.read_html(io.StringIO(table_html))[0]
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
        print(f"  🧬 Parsing ABRicate: {file_path.name}")
        try:
            with open(file_path, 'r', encoding='utf-8') as f:
                html_content = f.read()
            soup = BeautifulSoup(html_content, 'html.parser')
            tables = soup.find_all('table')
            if len(tables) < 2:
                return 'unknown', {}, {}
            db_name = 'unknown'
            filename = file_path.name.lower()
            for db in self.abricate_databases:
                if db in filename:
                    db_name = db
                    break
            title_tag = soup.find('title')
            if title_tag:
                title_text = title_tag.get_text().lower()
                for db in self.abricate_databases:
                    if db in title_text:
                        db_name = db
                        break
            genes_by_genome = {}
            df1 = self.parse_html_table(html_content, 0)
            if not df1.empty:
                df1.columns = [col.strip() for col in df1.columns]
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
                        genes_col = None
                        for col in df1.columns:
                            if 'genes' in col.lower() or 'detected' in col.lower():
                                genes_col = col
                                break
                        if genes_col and pd.notna(row.get(genes_col)):
                            gene_str = str(row[genes_col])
                            genes = [g.strip() for g in gene_str.split(',') if g.strip()]
                        genes_by_genome[sample] = genes
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
    
    def parse_mutation_summary_html(self, file_path: Path) -> Dict[str, Any]:
        """Parse mutation_summary.html into gene-centric mutation data."""
        print(f"  🧬 Parsing mutation summary HTML: {file_path.name}")
        try:
            with open(file_path, 'r', encoding='utf-8') as f:
                html_content = f.read()
            soup = BeautifulSoup(html_content, 'html.parser')
            
            # Find the mutation table – look for a table containing the text "Gene"
            mutation_table = None
            for table in soup.find_all('table'):
                if table.find(string=re.compile(r'Gene', re.I)) and table.find(string=re.compile(r'Mutation', re.I)):
                    mutation_table = table
                    break
            
            if not mutation_table:
                print("    ⚠️ Could not find mutation table in HTML")
                return {}
            
            # Parse headers – they may be in <th> or in the first row's <td>
            header_row = None
            thead = mutation_table.find('thead')
            if thead:
                header_row = thead.find('tr')
            if not header_row:
                # Try the first row of the table (could be in <tbody> or direct)
                first_row = mutation_table.find('tr')
                if first_row:
                    header_row = first_row
            
            if not header_row:
                print("    ⚠️ Could not find header row in mutation table")
                return {}
            
            headers = []
            for cell in header_row.find_all(['th', 'td']):
                text = cell.get_text().strip()
                if text:
                    headers.append(text)
            
            # Find column indices based on header text
            col_idx = {}
            for idx, h in enumerate(headers):
                h_lower = h.lower()
                if 'gene' in h_lower:
                    col_idx['gene'] = idx
                elif 'mutation' in h_lower:
                    col_idx['mutation'] = idx
                elif 'count' in h_lower:
                    col_idx['count'] = idx
                elif 'genome' in h_lower:
                    col_idx['genomes'] = idx
                elif 'class' in h_lower:
                    col_idx['class'] = idx
                elif 'subclass' in h_lower:
                    col_idx['subclass'] = idx
            
            required = ['gene', 'mutation', 'count', 'genomes']
            for req in required:
                if req not in col_idx:
                    print(f"    ⚠️ Missing required column: {req}. Found headers: {headers}")
                    return {}
            
            # Parse data rows – skip the header row
            tbody = mutation_table.find('tbody')
            if tbody:
                rows = tbody.find_all('tr')
            else:
                rows = mutation_table.find_all('tr')[1:]  # skip first row (header)
            
            mutations_list = []
            genome_counts = defaultdict(int)
            
            for row in rows:
                cells = row.find_all('td')
                if len(cells) <= max(col_idx.values()):
                    continue
                
                gene = cells[col_idx['gene']].get_text().strip()
                mutation = cells[col_idx['mutation']].get_text().strip()
                count_str = cells[col_idx['count']].get_text().strip()
                count_match = re.search(r'(\d+)', count_str)
                count = int(count_match.group(1)) if count_match else 0
                
                genomes_str = cells[col_idx['genomes']].get_text().strip()
                genomes = [g.strip() for g in genomes_str.split(',') if g.strip()]
                if not genomes:
                    continue
                
                for g in genomes:
                    genome_counts[g] += 1
                
                class_name = ''
                if 'class' in col_idx:
                    class_name = cells[col_idx['class']].get_text().strip()
                subclass = ''
                if 'subclass' in col_idx:
                    subclass = cells[col_idx['subclass']].get_text().strip()
                
                mutations_list.append({
                    'gene': gene,
                    'mutation': mutation,
                    'class': class_name,
                    'subclass': subclass,
                    'count': count,
                    'genomes': genomes
                })
            
            mutations_list.sort(key=lambda x: x['count'], reverse=True)
            print(f"    ✓ Parsed {len(mutations_list)} unique mutations across {len(genome_counts)} genomes")
            return {
                'mutations': mutations_list,
                'genome_mutation_counts': dict(genome_counts)
            }
            
        except Exception as e:
            print(f"    ❌ Error parsing mutation summary HTML: {e}")
            import traceback
            traceback.print_exc()
            return {}


class StaphDataAnalyzer:
    """Analyzes data for S. aureus gene-centric reporting"""
    
    def __init__(self):
        self.critical_amr_genes = {
            'meca', 'mecA', 'mecc', 'mecC', 'vana', 'vanA', 'vanb', 'vanB',
            'vanc', 'vanC', 'erma', 'ermA', 'ermb', 'ermB', 'ermc', 'ermC',
            'msra', 'msrA', 'mphc', 'mphC', 'tetk', 'tetK', 'tetm', 'tetM', 'tetl', 'tetL'
        }
        self.high_priority_amr = [
            'mecA', 'mecC', 'vanA', 'vanB', 'ermA', 'ermB', 'ermC',
            'msrA', 'mphC', 'tetK', 'tetM', 'aacA-aphD', 'aac(6\')-aph(2\'\')',
            'ant(4\')-Ia', 'ant(6)-Ia', 'aph(3\')-IIIa', 'satA', 'dfrA', 'dfrG', 'cat'
        ]
        self.critical_virulence_genes = {
            'luks-pv', 'lukS-PV', 'lukf-pv', 'lukF-PV',
            'tsst',
            'sea', 'seb', 'sec', 'sed', 'see',
            'seg', 'seh', 'sei', 'sej', 'sek', 'sel', 'sem', 'sen', 'seo', 'sep', 'seq', 'ser', 'seu',
            'eta', 'etb',
            'hla', 'hlb', 'hlg', 'hld',
        }
        self.high_priority_virulence = [
            'lukF-PV', 'lukS-PV', 'tsst', 'sea', 'seb', 'sec', 'sed', 'see',
            'seg', 'seh', 'sei', 'sej', 'sek', 'sel', 'sem', 'sen', 'seo', 'sep',
            'eta', 'etb', 'hla', 'hlb', 'hlg', 'hld'
        ]
        self.sccmec_types = {
            'I', 'II', 'III', 'IV', 'V', 'VI', 'VII', 'VIII', 'IX', 'X', 'XI', 'XII'
        }
    
    def create_gene_centric_tables(self, integrated_data: Dict[str, Any]) -> Dict[str, Any]:
        gene_centric = {
            'amr_databases': {},
            'virulence_databases': {},
            'plasmid_databases': {},
            'bacmet_databases': {},
            'combined_gene_frequencies': []
        }
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
        all_genes = []
        for db_type in ['amr_databases', 'virulence_databases', 'plasmid_databases', 'bacmet_databases']:
            for db_name, genes in gene_centric.get(db_type, {}).items():
                for gene_data in genes:
                    all_genes.append(gene_data)
        all_genes.sort(key=lambda x: x['count'], reverse=True)
        gene_centric['combined_gene_frequencies'] = all_genes
        return gene_centric
    
    def create_cross_genome_patterns(self, integrated_data: Dict[str, Any]) -> Dict[str, Any]:
        patterns = {
            'mlst_distribution': Counter(),
            'spa_type_distribution': Counter(),
            'sccmec_distribution': Counter(),
            'mrsa_status_distribution': Counter(),
            'mlst_spa_combinations': defaultdict(list),
            'mlst_sccmec_combinations': defaultdict(list),
            'spa_sccmec_combinations': defaultdict(list),
            'triple_combinations': defaultdict(list),
            'gene_cooccurrence': defaultdict(Counter),
            'high_risk_combinations': []
        }
        samples_data = integrated_data.get('samples', {})
        gene_centric = integrated_data.get('gene_centric', {})
        sample_genes = defaultdict(list)
        for db_type in ['amr_databases', 'virulence_databases']:
            for db_name, genes in gene_centric.get(db_type, {}).items():
                for gene_data in genes:
                    for genome in gene_data['genomes']:
                        if gene_data['gene'] not in sample_genes[genome]:
                            sample_genes[genome].append(gene_data['gene'])
        for sample, data in samples_data.items():
            mlst = data.get('typing', {}).get('MLST', 'ND')
            spa_type = data.get('typing', {}).get('spa_Type', 'ND')
            sccmec_type = data.get('typing', {}).get('SCCmec_Type', 'ND')
            mrsa_status = data.get('typing', {}).get('MRSA_Status', 'ND')
            if mlst != 'ND':
                patterns['mlst_distribution'][mlst] += 1
            if spa_type != 'ND':
                patterns['spa_type_distribution'][spa_type] += 1
            if sccmec_type != 'ND' and sccmec_type != 'Not Assigned':
                patterns['sccmec_distribution'][sccmec_type] += 1
            if mrsa_status != 'ND':
                patterns['mrsa_status_distribution'][mrsa_status] += 1
            if mlst != 'ND' and spa_type != 'ND':
                patterns['mlst_spa_combinations'][f"{mlst} - {spa_type}"].append(sample)
            if mlst != 'ND' and sccmec_type != 'ND' and sccmec_type != 'Not Assigned':
                patterns['mlst_sccmec_combinations'][f"{mlst} - {sccmec_type}"].append(sample)
            if spa_type != 'ND' and sccmec_type != 'ND' and sccmec_type != 'Not Assigned':
                patterns['spa_sccmec_combinations'][f"{spa_type} - {sccmec_type}"].append(sample)
            if mlst != 'ND' and spa_type != 'ND' and sccmec_type != 'ND' and sccmec_type != 'Not Assigned':
                patterns['triple_combinations'][f"{mlst} - {spa_type} - {sccmec_type}"].append(sample)
            genes = sample_genes.get(sample, [])
            for i, gene1 in enumerate(genes):
                for gene2 in genes[i+1:]:
                    patterns['gene_cooccurrence'][gene1][gene2] += 1
            amr_genes = data.get('amrfinder', {}).get('all_genes', [])
            virulence_genes = data.get('abricate_databases', {}).get('vfdb', [])
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
            'mutation': '#00BCD4',
            'patterns': '#3F51B5',
            'aiguide': '#00BCD4',
            'citation': '#8BC34A',
            'funding': '#FFC107',
            'export': '#9E9E9E'
        }
    
    def generate_main_report(self, integrated_data: Dict[str, Any], output_dir: Path) -> str:
        print("\n🎨 Generating STAPHSCOPE ULTIMATE HTML report...")
        samples_data = integrated_data.get('samples', {})
        patterns = integrated_data.get('patterns', {})
        gene_centric = integrated_data.get('gene_centric', {})
        metadata = integrated_data.get('metadata', {})
        html = self._create_ultimate_html(
            metadata=metadata,
            samples_data=samples_data,
            patterns=patterns,
            gene_centric=gene_centric,
            integrated_data=integrated_data
        )
        output_file = output_dir / "staphscope_ultimate_report.html"
        with open(output_file, 'w', encoding='utf-8') as f:
            f.write(html)
        print(f"    ✅ HTML report saved: {output_file}")
        return str(output_file)
    
    def _create_ultimate_html(self, **kwargs) -> str:
        samples_data = kwargs.get('samples_data', {})
        sample_typing_js = {}
        for sample, data in samples_data.items():
            typing = data.get('typing', {})
            sample_typing_js[sample] = {
                "MLST": typing.get('MLST', 'ND'),
                "spa": typing.get('spa_Type', 'ND'),
                "SCCmec": typing.get('SCCmec_Type', 'ND'),
                "MRSA": typing.get('MRSA_Status', 'ND')
            }
        sample_typing_json = json.dumps(sample_typing_js)
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
            --mutation-color: #00BCD4;
            --patterns-color: #3F51B5;
            --aiguide-color: #00BCD4;
            --citation-color: #8BC34A;
            --funding-color: #FFC107;
            --export-color: #9E9E9E;
        }
        * { margin: 0; padding: 0; box-sizing: border-box; }
        body { font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif; line-height: 1.6; color: #333; background: #f5f5f5; min-width: 1200px; }
        .container { max-width: none; margin: 0 auto; padding: 20px; width: 100%; overflow-x: auto; }
        .main-header { background: linear-gradient(135deg, #006400 0%, #228B22 100%); color: white; padding: 30px; border-radius: 15px; box-shadow: 0 10px 30px rgba(0,0,0,0.2); margin-bottom: 30px; text-align: center; }
        .main-header h1 { font-size: 2.8em; margin-bottom: 10px; color: white; }
        .metadata-bar { background: rgba(255,255,255,0.1); padding: 15px; border-radius: 10px; margin: 20px 0; display: flex; justify-content: space-around; flex-wrap: wrap; gap: 15px; backdrop-filter: blur(10px); }
        .metadata-item { display: flex; align-items: center; gap: 8px; font-size: 0.95em; }
        .dashboard-grid { display: grid; grid-template-columns: repeat(auto-fit, minmax(250px, 1fr)); gap: 20px; margin-bottom: 30px; }
        .dashboard-card { background: white; padding: 25px; border-radius: 12px; box-shadow: 0 5px 20px rgba(0,0,0,0.1); text-align: center; transition: all 0.3s ease; cursor: pointer; border-left: 5px solid; position: relative; overflow: hidden; }
        .dashboard-card::before { content: ''; position: absolute; top: 0; left: 0; right: 0; height: 3px; background: linear-gradient(90deg, transparent, rgba(255,255,255,0.8), transparent); }
        .dashboard-card:hover { transform: translateY(-10px); box-shadow: 0 15px 30px rgba(0,0,0,0.2); }
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
        .card-funding { border-left-color: var(--funding-color); }
        .card-export { border-left-color: var(--export-color); }
        .card-number { font-size: 3em; font-weight: bold; margin: 15px 0; background: linear-gradient(90deg, #006400, #228B22); -webkit-background-clip: text; -webkit-text-fill-color: transparent; }
        .tab-navigation { display: flex; gap: 5px; margin-bottom: 20px; flex-wrap: wrap; background: white; padding: 15px; border-radius: 12px; box-shadow: 0 5px 20px rgba(0,0,0,0.1); position: sticky; top: 10px; z-index: 100; }
        .tab-button { padding: 12px 20px; background: #f5f5f5; border: none; border-radius: 8px; cursor: pointer; font-weight: 600; color: #666; transition: all 0.3s ease; display: flex; align-items: center; gap: 8px; position: relative; overflow: hidden; font-size: 0.9em; }
        .tab-button::after { content: ''; position: absolute; bottom: 0; left: 50%; right: 50%; height: 3px; background: currentColor; transition: all 0.3s ease; }
        .tab-button:hover::after { left: 10%; right: 10%; }
        .tab-button.active { color: white; }
        .tab-button.active::after { left: 10%; right: 10%; }
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
        .tab-button.mutation.active { background: var(--mutation-color); }
        .tab-button.patterns.active { background: var(--patterns-color); }
        .tab-button.aiguide.active { background: var(--aiguide-color); }
        .tab-button.citation.active { background: var(--citation-color); }
        .tab-button.funding.active { background: var(--funding-color); }
        .tab-button.export.active { background: var(--export-color); }
        .tab-content { display: none; background: white; padding: 30px; border-radius: 15px; box-shadow: 0 10px 30px rgba(0,0,0,0.1); margin-bottom: 30px; animation: fadeIn 0.5s ease; width: 100%; overflow-x: auto; }
        .tab-content.active { display: block; }
        @keyframes fadeIn { from { opacity: 0; transform: translateY(20px); } to { opacity: 1; transform: translateY(0); } }
        .section-header { color: #2c3e50; margin-bottom: 25px; padding-bottom: 15px; border-bottom: 3px solid; font-size: 1.8em; display: flex; align-items: center; justify-content: space-between; }
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
        .mutation-header { border-color: var(--mutation-color); }
        .patterns-header { border-color: var(--patterns-color); }
        .aiguide-header { border-color: var(--aiguide-color); }
        .citation-header { border-color: var(--citation-color); }
        .funding-header { border-color: var(--funding-color); }
        .export-header { border-color: var(--export-color); }
        .data-table { width: 100%; border-collapse: collapse; margin: 20px 0; font-size: 0.95em; box-shadow: 0 2px 10px rgba(0,0,0,0.1); border-radius: 8px; overflow: hidden; table-layout: auto; }
        .data-table th { background: #2c3e50; color: white; padding: 15px; text-align: left; font-weight: 600; position: sticky; top: 0; white-space: nowrap; cursor: pointer; }
        .data-table th:hover { background: #1a252f; }
        .data-table td { padding: 12px; border-bottom: 1px solid #e0e0e0; vertical-align: top; word-wrap: break-word; white-space: nowrap; }
        .data-table tr:hover { background: #f8f9fa; }
        .scrollable-table { max-height: none; overflow-y: auto; border: 1px solid #e0e0e0; border-radius: 8px; margin: 20px 0; width: 100%; }
        .master-scrollable-container { width: 100%; overflow-x: auto; border: 1px solid #e0e0e0; border-radius: 8px; margin: 20px 0; }
        .genome-list { display: flex; flex-wrap: wrap; gap: 5px; max-height: 200px; overflow-y: auto; padding: 5px; background: #f8f9fa; border-radius: 5px; }
        .genome-group { margin-bottom: 10px; width: 100%; }
        .genome-group-header { font-weight: bold; background: #e0e0e0; padding: 4px 8px; border-radius: 4px; margin: 5px 0; font-size: 0.85em; display: inline-block; }
        .genome-group-tags { display: flex; flex-wrap: wrap; gap: 5px; margin-left: 10px; }
        .genome-tag { display: inline-block; background: #e6ffe6; color: #006400; padding: 3px 10px; border-radius: 12px; font-size: 0.85em; border: 1px solid #b3ffb3; white-space: nowrap; margin: 2px; }
        .genome-tag.highlight { background-color: #ffff99 !important; color: #000 !important; border: 1px solid #ffc107; }
        .search-box { width: 100%; padding: 12px; margin-bottom: 20px; border: 2px solid #e0e0e0; border-radius: 8px; font-size: 1em; transition: all 0.3s ease; }
        .search-box:focus { outline: none; border-color: #006400; box-shadow: 0 0 0 3px rgba(0, 100, 0, 0.1); }
        .badge { display: inline-block; padding: 5px 15px; border-radius: 20px; font-size: 0.85em; font-weight: 600; margin: 2px; }
        .badge-mrsa { background: #8B0000; color: white; }
        .badge-mssa { background: #4682B4; color: white; }
        .badge-critical { background: #DC143C; color: white; }
        .badge-high { background: #FF4500; color: white; }
        .badge-medium { background: #FF8C00; color: black; }
        .badge-low { background: #32CD32; color: white; }
        .alert-box { padding: 20px; border-radius: 10px; margin: 20px 0; display: flex; align-items: center; gap: 20px; border-left: 5px solid; }
        .alert-success { background: #d4edda; color: #155724; border-left-color: #28a745; }
        .alert-warning { background: #fff3cd; color: #856404; border-left-color: #ffc107; }
        .alert-danger { background: #f8d7da; color: #721c24; border-left-color: #dc3545; }
        .alert-info { background: #d1ecf1; color: #0c5460; border-left-color: #17a2b8; }
        .action-buttons { display: flex; gap: 10px; margin: 20px 0; flex-wrap: wrap; }
        .action-btn { padding: 10px 20px; border: none; border-radius: 8px; cursor: pointer; font-weight: 600; display: flex; align-items: center; gap: 8px; transition: all 0.3s ease; }
        .action-btn:hover { transform: translateY(-2px); box-shadow: 0 5px 15px rgba(0,0,0,0.2); }
        .btn-primary { background: #006400; color: white; }
        .btn-success { background: #28a745; color: white; }
        .btn-danger { background: #dc3545; color: white; }
        .btn-warning { background: #ffc107; color: black; }
        .btn-info { background: #17a2b8; color: white; }
        .btn-secondary { background: #6c757d; color: white; }
        .btn-light { background: #f8f9fa; color: #212529; border: 1px solid #dee2e6; }
        .database-section { margin: 30px 0; padding: 25px; border-radius: 12px; background: #f8f9fa; box-shadow: 0 3px 15px rgba(0,0,0,0.08); }
        .database-header { font-size: 1.4em; color: #2c3e50; margin-bottom: 20px; padding-bottom: 10px; border-bottom: 2px solid #006400; display: flex; align-items: center; justify-content: space-between; }
        .print-section-btn { background: #006400; color: white; border: none; border-radius: 5px; padding: 8px 15px; cursor: pointer; display: flex; align-items: center; gap: 5px; font-size: 0.9em; }
        .print-section-btn:hover { background: #228B22; }
        .footer { text-align: center; padding: 30px; color: white; margin-top: 40px; border-radius: 15px; background: linear-gradient(135deg, #2c3e50 0%, #34495e 100%); }
        .mrsa-highlight { background-color: #ffe6e6 !important; border-left: 3px solid #8B0000 !important; }
        .sort-icon { margin-left: 5px; font-size: 0.8em; opacity: 0.6; }
        .grouping-controls { background: #f0f7f0; padding: 12px; border-radius: 8px; margin: 15px 0; display: flex; flex-wrap: wrap; gap: 10px; align-items: center; border-left: 4px solid #006400; }
        .grouping-controls label { font-weight: bold; margin-right: 5px; }
        .group-btn { background: white; border: 1px solid #006400; color: #006400; padding: 6px 12px; border-radius: 20px; cursor: pointer; font-size: 0.85em; transition: all 0.2s; }
        .group-btn:hover { background: #006400; color: white; }
        .group-btn.active { background: #006400; color: white; }
        .accordion { margin: 20px 0; }
        .accordion-item { background: #f8f9fa; border: 1px solid #dee2e6; margin-bottom: 10px; border-radius: 8px; overflow: hidden; }
        .accordion-header { background: #e9ecef; padding: 12px 20px; cursor: pointer; font-weight: bold; color: #1e3a8a; display: flex; justify-content: space-between; align-items: center; }
        .accordion-header:hover { background: #dee2e6; }
        .accordion-content { padding: 15px 20px; border-top: 1px solid #dee2e6; background: white; }
        .citation-list { list-style: none; padding-left: 0; }
        .citation-list li { margin-bottom: 15px; padding-bottom: 10px; border-bottom: 1px solid #e0e0e0; }
        .citation-list li:last-child { border-bottom: none; }
        .copy-btn { background: #006400; color: white; border: none; padding: 4px 12px; border-radius: 20px; cursor: pointer; font-size: 0.8em; margin-left: 10px; }
        .copy-btn:hover { background: #228B22; }
        @media print { body * { visibility: hidden; } .tab-content.active, .tab-content.active * { visibility: visible; } .tab-content.active { position: absolute; left: 0; top: 0; width: 100%; padding: 20px; box-shadow: none; border-radius: 0; } .print-section-btn, .tab-navigation, .dashboard-grid, .search-box, .action-buttons, .grouping-controls { display: none !important; } .data-table { page-break-inside: auto; } .data-table tr { page-break-inside: avoid; page-break-after: auto; } }
        @media (max-width: 768px) { .container { padding: 10px; } .main-header { padding: 20px; } .main-header h1 { font-size: 2em; } .tab-button { padding: 8px 12px; font-size: 0.8em; } .dashboard-grid { grid-template-columns: repeat(auto-fit, minmax(180px, 1fr)); } .data-table { font-size: 0.8em; } body { min-width: auto; overflow-x: auto; } }
        </style>
        """
        js = f"""
        <script>
        var sampleTyping = {sample_typing_json};
        var originalGenomeLists = {{}};
        function switchTab(tabName) {{
            document.querySelectorAll('.tab-content').forEach(tab => tab.classList.remove('active'));
            document.querySelectorAll('.tab-button').forEach(button => button.classList.remove('active'));
            document.getElementById(tabName + '-tab').classList.add('active');
            event.currentTarget.classList.add('active');
            window.location.hash = tabName;
        }}
        function searchTable(tableId, searchId) {{
            const input = document.getElementById(searchId);
            const filter = input.value.toUpperCase();
            const table = document.getElementById(tableId);
            const rows = table.getElementsByTagName('tr');
            for (let i = 1; i < rows.length; i++) {{
                const cells = rows[i].getElementsByTagName('td');
                let found = false;
                for (let j = 0; j < cells.length; j++) {{
                    const cell = cells[j];
                    if (cell) {{
                        const txtValue = cell.textContent || cell.innerText;
                        if (txtValue.toUpperCase().indexOf(filter) > -1) {{
                            found = true;
                            break;
                        }}
                    }}
                }}
                rows[i].style.display = found ? '' : 'none';
            }}
        }}
        function highlightGenome(tableId, searchId) {{
            const filter = document.getElementById(searchId).value.toUpperCase().trim();
            const table = document.getElementById(tableId);
            const allTags = table.querySelectorAll('.genome-tag');
            allTags.forEach(tag => tag.classList.remove('highlight'));
            if (filter === '') return;
            allTags.forEach(tag => {{
                if (tag.textContent.toUpperCase().indexOf(filter) > -1) {{
                    tag.classList.add('highlight');
                }}
            }});
        }}
        function getTypingValue(genome, groupBy) {{
            var info = sampleTyping[genome];
            if (!info) return "Unknown";
            if (groupBy === "MLST") return info.MLST;
            if (groupBy === "spa") return info.spa;
            if (groupBy === "SCCmec") return info.SCCmec;
            if (groupBy === "ST-spa") return info.MLST + " - " + info.spa;
            if (groupBy === "ST-SCCmec") return info.MLST + " - " + info.SCCmec;
            if (groupBy === "spa-SCCmec") return info.spa + " - " + info.SCCmec;
            if (groupBy === "triple") return info.MLST + " - " + info.spa + " - " + info.SCCmec;
            return "Unknown";
        }}
        function groupRowGenomes(row, groupBy, originalList) {{
            let genomesCell = null;
            for (let i = 0; i < row.cells.length; i++) {{
                if (row.cells[i].querySelector('.genome-list')) {{
                    genomesCell = row.cells[i];
                    break;
                }}
            }}
            if (!genomesCell) {{
                console.warn("Could not find genomes cell in row");
                return;
            }}
            var genomes = originalList.slice();
            if (genomes.length === 0) {{
                genomesCell.innerHTML = '<div class="genome-list">None</div>';
                return;
            }}
            var groups = {{}};
            genomes.forEach(function(genome) {{
                var key = getTypingValue(genome, groupBy);
                if (!groups[key]) groups[key] = [];
                groups[key].push(genome);
            }});
            var html = '<div class="genome-list">';
            for (var key in groups) {{
                var tags = groups[key].map(g => `<span class="genome-tag">${{g}}</span>`).join('');
                html += `<div class="genome-group"><div class="genome-group-header">${{key}}</div><div class="genome-group-tags">${{tags}}</div></div>`;
            }}
            html += '</div>';
            genomesCell.innerHTML = html;
        }}
        function groupGenomesByTyping(tableId, groupBy) {{
            var table = document.getElementById(tableId);
            if (!table) {{
                console.error("Table not found:", tableId);
                return;
            }}
            var tbody = table.tBodies[0];
            if (!tbody) {{
                console.error("No tbody found in table", tableId);
                return;
            }}
            var rows = tbody.rows;
            for (var i = 0; i < rows.length; i++) {{
                var row = rows[i];
                var geneNameCell = row.cells[0];
                if (!geneNameCell) continue;
                var geneName = geneNameCell.textContent.trim().replace(/⚠️/g, '').trim();
                if (!originalGenomeLists[geneName]) {{
                    var genomesCell = null;
                    for (var j = 0; j < row.cells.length; j++) {{
                        if (row.cells[j].querySelector('.genome-list')) {{
                            genomesCell = row.cells[j];
                            break;
                        }}
                    }}
                    if (genomesCell) {{
                        var tags = genomesCell.querySelectorAll('.genome-tag');
                        var genomes = Array.from(tags).map(tag => tag.textContent.trim());
                        originalGenomeLists[geneName] = genomes;
                    }} else {{
                        originalGenomeLists[geneName] = [];
                    }}
                }}
            }}
            for (var i = 0; i < rows.length; i++) {{
                var row = rows[i];
                var geneNameCell = row.cells[0];
                if (!geneNameCell) continue;
                var geneName = geneNameCell.textContent.trim().replace(/⚠️/g, '').trim();
                var original = originalGenomeLists[geneName] || [];
                groupRowGenomes(row, groupBy, original);
            }}
            var container = table.closest('.tab-content');
            if (container) {{
                var btns = container.querySelectorAll('.group-btn');
                btns.forEach(btn => btn.classList.remove('active'));
                var activeBtn = container.querySelector(`.group-btn[data-group="${{groupBy}}"]`);
                if (activeBtn) activeBtn.classList.add('active');
            }}
        }}
        function resetGenomeList(tableId) {{
            var table = document.getElementById(tableId);
            if (!table) return;
            var tbody = table.tBodies[0];
            if (!tbody) return;
            var rows = tbody.rows;
            for (var i = 0; i < rows.length; i++) {{
                var row = rows[i];
                var geneNameCell = row.cells[0];
                if (!geneNameCell) continue;
                var geneName = geneNameCell.textContent.trim().replace(/⚠️/g, '').trim();
                var original = originalGenomeLists[geneName] || [];
                var genomesCell = null;
                for (var j = 0; j < row.cells.length; j++) {{
                    if (row.cells[j].querySelector('.genome-list')) {{
                        genomesCell = row.cells[j];
                        break;
                    }}
                }}
                if (genomesCell) {{
                    var tags = original.map(g => `<span class="genome-tag">${{g}}</span>`).join('');
                    genomesCell.innerHTML = `<div class="genome-list">${{tags}}</div>`;
                }}
            }}
            var container = table.closest('.tab-content');
            if (container) {{
                var btns = container.querySelectorAll('.group-btn');
                btns.forEach(btn => btn.classList.remove('active'));
            }}
        }}
        function sortTable(tableId, colIndex, type) {{
            const table = document.getElementById(tableId);
            const tbody = table.tBodies[0];
            const rows = Array.from(tbody.rows);
            const isAscending = table.getAttribute('data-sort-dir') !== 'asc';
            rows.sort((a, b) => {{
                let aVal = a.cells[colIndex].innerText.trim();
                let bVal = b.cells[colIndex].innerText.trim();
                if (type === 'number') {{
                    aVal = parseFloat(aVal.replace(/,/g, '')) || 0;
                    bVal = parseFloat(bVal.replace(/,/g, '')) || 0;
                    return isAscending ? aVal - bVal : bVal - aVal;
                }} else {{
                    return isAscending ? aVal.localeCompare(bVal) : bVal.localeCompare(aVal);
                }}
            }});
            tbody.append(...rows);
            table.setAttribute('data-sort-dir', isAscending ? 'asc' : 'desc');
            const headers = table.querySelectorAll('th');
            headers.forEach((th, idx) => {{
                const icon = th.querySelector('.sort-icon');
                if (icon) icon.innerHTML = '⇅';
            }});
            const currentHeader = headers[colIndex];
            const icon = currentHeader.querySelector('.sort-icon');
            if (icon) icon.innerHTML = isAscending ? '↑' : '↓';
        }}
        function printSection(sectionId) {{
            const content = document.getElementById(sectionId);
            const printWindow = window.open('', '_blank');
            printWindow.document.write('<html><head><title>Print Section</title>');
            printWindow.document.write('<style>' + document.querySelector('style').textContent + '</style>');
            printWindow.document.write('</head><body>');
            printWindow.document.write(content.innerHTML);
            printWindow.document.write('</body></html>');
            printWindow.document.close();
            printWindow.print();
        }}
        function exportTableToCSV(tableId, filename) {{
            const table = document.getElementById(tableId);
            const rows = table.querySelectorAll('tr');
            const csv = [];
            for (let i = 0; i < rows.length; i++) {{
                const row = [], cols = rows[i].querySelectorAll('td, th');
                for (let j = 0; j < cols.length; j++) {{
                    row.push('"' + (cols[j].innerText || '').replace(/"/g, '""') + '"');
                }}
                csv.push(row.join(','));
            }}
            const csvFile = new Blob([csv.join('\\n')], {{type: 'text/csv'}});
            const downloadLink = document.createElement('a');
            downloadLink.download = filename;
            downloadLink.href = window.URL.createObjectURL(csvFile);
            downloadLink.style.display = 'none';
            document.body.appendChild(downloadLink);
            downloadLink.click();
            document.body.removeChild(downloadLink);
        }}
        document.addEventListener('DOMContentLoaded', function() {{
            const hash = window.location.hash.substring(1);
            if (hash) {{
                const tabButton = document.querySelector(`.tab-button.${{hash}}`);
                if (tabButton) tabButton.click();
            }} else {{
                document.querySelector('.tab-button').click();
            }}
            document.querySelectorAll('.data-table').forEach(table => {{
                const headers = table.querySelectorAll('th');
                headers.forEach((header, idx) => {{
                    const type = header.getAttribute('data-sort') || 'string';
                    header.style.cursor = 'pointer';
                    header.addEventListener('click', () => sortTable(table.id, idx, type));
                    const icon = document.createElement('span');
                    icon.className = 'sort-icon';
                    icon.innerHTML = '⇅';
                    header.appendChild(icon);
                }});
            }});
            document.querySelectorAll('.accordion-header').forEach(header => {{
                header.addEventListener('click', function() {{
                    const content = this.nextElementSibling;
                    content.style.display = content.style.display === 'block' ? 'none' : 'block';
                }});
            }});
            document.querySelectorAll('.copy-btn').forEach(btn => {{
                btn.addEventListener('click', function() {{
                    const citation = this.getAttribute('data-citation');
                    navigator.clipboard.writeText(citation).then(() => {{
                        const originalText = this.innerHTML;
                        this.innerHTML = '✓ Copied!';
                        setTimeout(() => {{ this.innerHTML = originalText; }}, 2000);
                    }});
                }});
            }});
        }});
        </script>
        """
        total_amr_genes = sum(len(genes) for genes in kwargs['gene_centric'].get('amr_databases', {}).values())
        total_virulence_genes = sum(len(genes) for genes in kwargs['gene_centric'].get('virulence_databases', {}).values())
        total_bacmet_genes = sum(len(genes) for genes in kwargs['gene_centric'].get('bacmet_databases', {}).values())
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
            <div class="main-header">
                <h1><i class="fas fa-bacteria"></i> STAPHSCOPE Ultimate S. aureus Analysis Report</h1>
                <p>Gene-Centric Cross-Genome Analysis with Genome Grouping by Typing</p>
                <div class="metadata-bar">
                    <div class="metadata-item"><i class="fas fa-calendar"></i><span>Generated: {kwargs['metadata'].get('analysis_date', 'Unknown')}</span></div>
                    <div class="metadata-item"><i class="fas fa-database"></i><span>Samples: {len(kwargs['samples_data'])}</span></div>
                    <div class="metadata-item"><i class="fas fa-user-md"></i><span>Tool: STAPHSCOPE Ultimate v2.2.0</span></div>
                    <div class="metadata-item"><i class="fas fa-university"></i><span>University of Ghana Medical School</span></div>
                </div>
            </div>
            <div class="dashboard-grid">
                <div class="dashboard-card card-summary" onclick="switchTab('summary')"><div class="card-number">{len(kwargs['samples_data'])}</div><div class="card-label">Total Samples</div><i class="fas fa-vial fa-2x" style="color: var(--summary-color); margin-top: 10px;"></i></div>
                <div class="dashboard-card card-mlst" onclick="switchTab('mlst')"><div class="card-number">{len(kwargs['patterns'].get('mlst_distribution', {}))}</div><div class="card-label">Unique STs</div><i class="fas fa-code-branch fa-2x" style="color: var(--mlst-color); margin-top: 10px;"></i></div>
                <div class="dashboard-card card-spa" onclick="switchTab('spa')"><div class="card-number">{len(kwargs['patterns'].get('spa_type_distribution', {}))}</div><div class="card-label">spa Types</div><i class="fas fa-dna fa-2x" style="color: var(--spa-color); margin-top: 10px;"></i></div>
                <div class="dashboard-card card-amr" onclick="switchTab('amr')"><div class="card-number">{total_amr_genes}</div><div class="card-label">AMR Genes</div><i class="fas fa-biohazard fa-2x" style="color: var(--amr-color); margin-top: 10px;"></i></div>
                <div class="dashboard-card card-virulence" onclick="switchTab('virulence')"><div class="card-number">{total_virulence_genes}</div><div class="card-label">Virulence Genes</div><i class="fas fa-virus fa-2x" style="color: var(--virulence-color); margin-top: 10px;"></i></div>
                <div class="dashboard-card card-patterns" onclick="switchTab('patterns')"><div class="card-number">{len(kwargs['patterns'].get('high_risk_combinations', []))}</div><div class="card-label">High-Risk Combos</div><i class="fas fa-project-diagram fa-2x" style="color: var(--patterns-color); margin-top: 10px;"></i></div>
            </div>
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
                <button class="tab-button mutation" onclick="switchTab('mutation')"><i class="fas fa-dna"></i> Mutations</button>
                <button class="tab-button patterns" onclick="switchTab('patterns')"><i class="fas fa-project-diagram"></i> Patterns</button>
                <button class="tab-button aiguide" onclick="switchTab('aiguide')"><i class="fas fa-robot"></i> AI Guide</button>
                <button class="tab-button citation" onclick="switchTab('citation')"><i class="fas fa-book"></i> Citation</button>
                <button class="tab-button funding" onclick="switchTab('funding')"><i class="fas fa-coffee"></i> Funding</button>
                <button class="tab-button export" onclick="switchTab('export')"><i class="fas fa-download"></i> Export</button>
            </div>
            <div id="summary-tab" class="tab-content active">
                <h2 class="section-header summary-header"><i class="fas fa-chart-pie"></i> Executive Summary<button class="print-section-btn" onclick="printSection('summary-tab')"><i class="fas fa-print"></i> Print</button></h2>
                {self._generate_summary_section(kwargs)}
            </div>
            <div id="sample_overview-tab" class="tab-content">
                <h2 class="section-header sample_overview-header"><i class="fas fa-list-alt"></i> Sample Overview<button class="print-section-btn" onclick="printSection('sample_overview-tab')"><i class="fas fa-print"></i> Print</button></h2>
                {self._generate_sample_overview_section(kwargs)}
            </div>
            <div id="qc-tab" class="tab-content">
                <h2 class="section-header qc-header"><i class="fas fa-chart-line"></i> FASTA Quality Control Metrics<button class="print-section-btn" onclick="printSection('qc-tab')"><i class="fas fa-print"></i> Print</button></h2>
                {self._generate_qc_section(kwargs)}
            </div>
            <div id="mlst-tab" class="tab-content">
                <h2 class="section-header mlst-header"><i class="fas fa-code-branch"></i> MLST Analysis<button class="print-section-btn" onclick="printSection('mlst-tab')"><i class="fas fa-print"></i> Print</button></h2>
                {self._generate_mlst_section(kwargs)}
            </div>
            <div id="spa-tab" class="tab-content">
                <h2 class="section-header spa-header"><i class="fas fa-dna"></i> spa Typing Analysis<button class="print-section-btn" onclick="printSection('spa-tab')"><i class="fas fa-print"></i> Print</button></h2>
                {self._generate_spa_section(kwargs)}
            </div>
            <div id="sccmec-tab" class="tab-content">
                <h2 class="section-header sccmec-header"><i class="fas fa-shield-alt"></i> SCCmec Typing Analysis<button class="print-section-btn" onclick="printSection('sccmec-tab')"><i class="fas fa-print"></i> Print</button></h2>
                {self._generate_sccmec_section(kwargs)}
            </div>
            <div id="mrsa-tab" class="tab-content">
                <h2 class="section-header mrsa-header"><i class="fas fa-skull-crossbones"></i> MRSA Analysis<button class="print-section-btn" onclick="printSection('mrsa-tab')"><i class="fas fa-print"></i> Print</button></h2>
                {self._generate_mrsa_section(kwargs)}
            </div>
            <div id="amr-tab" class="tab-content">
                <h2 class="section-header amr-header"><i class="fas fa-biohazard"></i> Antimicrobial Resistance Genes<button class="print-section-btn" onclick="printSection('amr-tab')"><i class="fas fa-print"></i> Print</button></h2>
                {self._generate_amr_section(kwargs)}
            </div>
            <div id="virulence-tab" class="tab-content">
                <h2 class="section-header virulence-header"><i class="fas fa-virus"></i> Virulence Genes<button class="print-section-btn" onclick="printSection('virulence-tab')"><i class="fas fa-print"></i> Print</button></h2>
                {self._generate_virulence_section(kwargs)}
            </div>
            <div id="bacmet-tab" class="tab-content">
                <h2 class="section-header bacmet-header"><i class="fas fa-flask"></i> BACMET: Biocide & Heavy Metal Resistance<button class="print-section-btn" onclick="printSection('bacmet-tab')"><i class="fas fa-print"></i> Print</button></h2>
                {self._generate_bacmet_section(kwargs)}
            </div>
            <div id="plasmids-tab" class="tab-content">
                <h2 class="section-header plasmids-header"><i class="fas fa-plug"></i> Plasmid Replicon Analysis<button class="print-section-btn" onclick="printSection('plasmids-tab')"><i class="fas fa-print"></i> Print</button></h2>
                {self._generate_plasmids_section(kwargs)}
            </div>
            <div id="mutation-tab" class="tab-content">
                <h2 class="section-header mutation-header"><i class="fas fa-dna"></i> Point Mutations (AMRfinderPlus)<button class="print-section-btn" onclick="printSection('mutation-tab')"><i class="fas fa-print"></i> Print</button></h2>
                {self._generate_mutation_section(kwargs)}
            </div>
            <div id="patterns-tab" class="tab-content">
                <h2 class="section-header patterns-header"><i class="fas fa-project-diagram"></i> Cross-Genome Pattern Discovery<button class="print-section-btn" onclick="printSection('patterns-tab')"><i class="fas fa-print"></i> Print</button></h2>
                {self._generate_pattern_discovery_section(kwargs)}
            </div>
            <div id="aiguide-tab" class="tab-content">
                <h2 class="section-header aiguide-header"><i class="fas fa-robot"></i> AI Assistant Guide<button class="print-section-btn" onclick="printSection('aiguide-tab')"><i class="fas fa-print"></i> Print</button></h2>
                {self._generate_aiguide_section(kwargs)}
            </div>
            <div id="citation-tab" class="tab-content">
                <h2 class="section-header citation-header"><i class="fas fa-book"></i> Citations & References<button class="print-section-btn" onclick="printSection('citation-tab')"><i class="fas fa-print"></i> Print</button></h2>
                {self._generate_citation_section(kwargs)}
            </div>
            <div id="funding-tab" class="tab-content">
                <h2 class="section-header funding-header"><i class="fas fa-coffee"></i> Funding & Support<button class="print-section-btn" onclick="printSection('funding-tab')"><i class="fas fa-print"></i> Print</button></h2>
                {self._generate_funding_section(kwargs)}
            </div>
            <div id="export-tab" class="tab-content">
                <h2 class="section-header export-header"><i class="fas fa-download"></i> Export Data<button class="print-section-btn" onclick="printSection('export-tab')"><i class="fas fa-print"></i> Print</button></h2>
                {self._generate_export_section(kwargs)}
            </div>
            <div class="footer">
                <h3>STAPHSCOPE Ultimate S. aureus Reporter v2.2.0</h3>
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
        samples_data = kwargs['samples_data']
        patterns = kwargs['patterns']
        gene_centric = kwargs['gene_centric']
        total_samples = len(samples_data)
        total_amr_genes = sum(len(genes) for genes in gene_centric.get('amr_databases', {}).values())
        total_virulence_genes = sum(len(genes) for genes in gene_centric.get('virulence_databases', {}).values())
        total_plasmids = sum(len(genes) for genes in gene_centric.get('plasmid_databases', {}).values())
        total_bacmet = sum(len(genes) for genes in gene_centric.get('bacmet_databases', {}).values())
        critical_findings = len(patterns.get('high_risk_combinations', []))
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
                <h3>📊 Analysis Overview – Gene‑Centric Approach + NEW Dynamic Grouping</h3>
                <p>This report analyzes <strong>{total_samples}</strong> <em>Staphylococcus aureus</em> genomes using a <strong>gene‑centric</strong> approach: each resistance or virulence gene is displayed with <strong>all genomes</strong> that carry it. This makes outbreak tracking and co‑occurrence analysis immediate and transparent.</p>
                <p><strong>NEW FEATURE: Group Genomes by Typing!</strong> In the AMR, Virulence, BACMET, Plasmids, and Mutations tabs, you will find a set of buttons that let you dynamically regroup the genome list by MLST, spa, SCCmec, or any combination (e.g., ST‑spa‑SCCmec). This allows you to instantly see, for example, which STs carry a specific gene/mutation, and exactly which genomes belong to each ST.</p>
                <p><strong>FULL DISCLOSURE:</strong> This idea was born while I was desperately trying to make <strong>Kleboscope</strong> do something useful (it's still in therapy). I liked it so much I stole… I mean, ported it to <strong>StaphScope</strong> first. Don't tell the <em>Klebsiella</em> people. They'll make me rename it "KleboGroup" and submit a competing paper.</p>
            </div>
        </div>
        <div class="alert-box alert-success">
            <i class="fas fa-magic fa-2x"></i>
            <div>
                <h3>📘 How to use the Grouping Feature</h3>
                <ol><li>Go to any gene‑centric table (AMR, Virulence, BACMET, Plasmids, or Mutations).</li><li>Above the table, you will see a bar labeled <strong>"Group genomes by:"</strong> with buttons like MLST, spa, SCCmec, ST‑spa, etc.</li><li>Click any button – the genome column will instantly reorganize, grouping genomes under sub‑headers (e.g., ST5, ST8).</li><li>Click <strong>"Reset"</strong> to return to the flat, scrollable list.</li><li>You can still use the search and highlight boxes to find specific genomes across all groups.</li></ol>
                <p><strong>Why is this powerful?</strong> Instead of manually scanning a long list of genome names, you can now see at a glance: "Is <em>mecA</em> confined to ST5 and ST8?" or "Does PVL appear only in spa type t008?" This turns your gene table into an instant epidemiological tool.</p>
                <p style="margin-top: 15px; font-style: italic; border-top: 1px dashed #28a745; padding-top: 10px;"><i class="fas fa-microbe"></i> <strong>P.S.</strong> If you don’t love this new feature, we won’t ask you to drop a star – we’ll just send <strong><em>mecA</em></strong> to your lab. And trust us, once it arrives, your β‑lactams will never be the same again. 😉 <strong>⭐ Star us on GitHub – keep MRSA where it belongs: in the report, not in your hospital.</strong></p>
            </div>
        </div>
        <div class="action-buttons">
            <button class="action-btn btn-primary" onclick="switchTab('amr')"><i class="fas fa-biohazard"></i> Try Grouping on AMR Genes</button>
            <button class="action-btn btn-success" onclick="switchTab('virulence')"><i class="fas fa-virus"></i> Try Grouping on Virulence Genes</button>
            <button class="action-btn btn-info" onclick="switchTab('mutation')"><i class="fas fa-dna"></i> Try Grouping on Mutations</button>
            <button class="action-btn btn-danger" onclick="switchTab('patterns')"><i class="fas fa-exclamation-triangle"></i> Check High-Risk Combos</button>
        </div>
        <h3><i class="fas fa-chart-bar"></i> Key Statistics</h3>
        <div class="scrollable-table"><table class="data-table"><thead><tr><th>Metric</th><th>Count</th><th>Details</th></tr></thead><tbody>
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
        </tbody></table></div>
        <h3 style="margin-top: 30px;"><i class="fas fa-lightbulb"></i> ✨ Key Features of This Report</h3>
        <div style="display: grid; grid-template-columns: repeat(auto-fit, minmax(300px, 1fr)); gap: 20px; margin: 20px 0;">
            <div class="database-section"><h4><i class="fas fa-gene"></i> Gene‑Centric Tables</h4><p>Each AMR, virulence, BACMET, or mutation is shown with <strong>all genomes</strong> that carry it – no more sample‑by‑sample searching.</p></div>
            <div class="database-section"><h4><i class="fas fa-layer-group"></i> Dynamic Grouping by Typing</h4><p>Reorganise genome lists instantly by MLST, spa, SCCmec, or any combination. Discover clones carrying specific genes/mutations in one click.</p></div>
            <div class="database-section"><h4><i class="fas fa-skull-crossbones"></i> MRSA Focus</h4><p>Dedicated MRSA analysis tab with ST‑spa, ST‑SCCmec, spa‑SCCmec, and triple combinations.</p></div>
            <div class="database-section"><h4><i class="fas fa-flask"></i> BACMET</h4><p>Biocide and heavy metal resistance genes – essential for hospital hygiene studies.</p></div>
            <div class="database-section"><h4><i class="fas fa-chart-line"></i> FASTA QC</h4><p>Assembly quality metrics (N50, contig count, GC%) with typical ranges for S. aureus.</p></div>
            <div class="database-section"><h4><i class="fas fa-dna"></i> Mutations</h4><p>AMRfinderPlus point mutations in gene‑centric view – detect resistance due to gyrA, parC, rpoB, 23S, etc.</p></div>
            <div class="database-section"><h4><i class="fas fa-print"></i> Section‑Specific Printing</h4><p>Print any tab individually with the “Print” button in each section header.</p></div>
            <div class="database-section"><h4><i class="fas fa-download"></i> Full Data Export</h4><p>All tables can be exported as CSV, and the complete dataset as JSON for downstream analysis or AI upload.</p></div>
            <div class="database-section"><h4><i class="fas fa-robot"></i> AI Assistant Guide</h4><p>Upload the JSON file to ChatGPT, Claude, or Gemini for interactive questions – with ethical guidelines.</p></div>
        </div>
        """
        return html
    
    def _generate_sample_overview_section(self, kwargs: Dict) -> str:
            """Generate sample overview with helpful info and sortable columns note"""
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
                    <p><strong>MRSA samples</strong> are highlighted in red.</p>
                    <p><strong>Note:</strong> "Not Assigned" in SCCmec type means no SCCmec cassette was detected (common in MSSA).</p>
                    <p><strong>Sortable columns:</strong> Click on any column header to sort the table.</p>
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
            """
            return html
    
    def _generate_qc_section(self, kwargs: Dict) -> str:
            """Generate FASTA QC section with typical ranges for S. aureus"""
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
                        <li><strong>Number of contigs</strong> – lower is better; fragmented assemblies may miss genes or break up SCCmec cassettes. <em>Typical S. aureus: &lt;200 contigs (good), &lt;100 (excellent).</em></li>
                        <li><strong>N50</strong> – the contig length such that 50% of the assembly is in contigs ≥ N50; higher N50 indicates better continuity. <em>Typical S. aureus: &gt;50 kb (good), &gt;100 kb (excellent).</em></li>
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
            html += "<tr></thead><tbody>"
            
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
                <div class="master-scrollable-container"><table id="mlst-sccmec-table" class="data-table"><thead><tr><th data-sort="string">ST-SCCmec Combination</th><th data-sort="number">Count</th><th data-sort="string">Samples</th><tr></thead><tbody>
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
            """Generate SCCmec section – no MRSA content here"""
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
        samples_data = kwargs['samples_data']
        
        # Check if we have any typing data at all
        has_typing = any(
            d.get('typing', {}).get('MLST', 'ND') != 'ND' or 
            d.get('typing', {}).get('spa_Type', 'ND') != 'ND' or
            d.get('typing', {}).get('SCCmec_Type', 'ND') != 'ND'
            for d in samples_data.values()
        )
        
        # If no typing data, show warning
        if not has_typing:
            return """
            <div class="alert-box alert-warning">
                <i class="fas fa-exclamation-triangle fa-2x"></i>
                <div>
                    <h3>⚠️ No Typing Data Available</h3>
                    <p>The comprehensive typing report (<code>*comprehensive*.html</code>) was not found or could not be parsed. MRSA analysis requires MLST, spa, SCCmec, and MRSA status.</p>
                    <p>Please ensure the StaphScope comprehensive report is present in the input directory and re-run the analysis.</p>
                </div>
            </div>
            """
        
        # Filter MRSA samples
        mrsa_samples = [s for s, d in samples_data.items() if 'MRSA' in d.get('typing', {}).get('MRSA_Status', '')]
        
        # Build MRSA-specific combinations (only if MRSA samples exist)
        mrsa_mlst_spa = defaultdict(list)
        mrsa_mlst_sccmec = defaultdict(list)
        mrsa_spa_sccmec = defaultdict(list)
        for sample in mrsa_samples:
            data = samples_data[sample]
            mlst = data.get('typing', {}).get('MLST', 'ND')
            spa = data.get('typing', {}).get('spa_Type', 'ND')
            scc = data.get('typing', {}).get('SCCmec_Type', 'ND')
            if mlst != 'ND' and spa != 'ND':
                mrsa_mlst_spa[f"{mlst} - {spa}"].append(sample)
            if mlst != 'ND' and scc != 'ND' and scc != 'Not Assigned':
                mrsa_mlst_sccmec[f"{mlst} - {scc}"].append(sample)
            if spa != 'ND' and scc != 'ND' and scc != 'Not Assigned':
                mrsa_spa_sccmec[f"{spa} - {scc}"].append(sample)
        
        # Start HTML
        html = f"""
        <div class="alert-box alert-danger">
            <i class="fas fa-skull-crossbones fa-2x"></i>
            <div>
                <h3>⚠️ MRSA (Methicillin‑Resistant S. aureus) – A Clinical Priority</h3>
                <p>MRSA is resistant to all β‑lactam antibiotics, including penicillins and cephalosporins. It is a major cause of hospital‑ and community‑acquired infections. The presence of <em>mecA</em> or <em>mecC</em> (carried on SCCmec) defines MRSA.</p>
                <p><strong>{len(mrsa_samples)} MRSA samples</strong> identified.</p>
            </div>
        </div>
        
        <h3>📊 MRSA vs MSSA Distribution</h3>
        <div class="scrollable-table">
            <table id="mrsa-status-table" class="data-table">
                <thead><tr><th data-sort="string">Status</th><th data-sort="number">Count</th><th data-sort="number">Percentage</th><th data-sort="string">Common STs</th><th data-sort="string">Common SCCmec Types</th></tr></thead>
                <tbody>
        """
        total = sum(mrsa_status_dist.values())
        if total == 0:
            html += "<tr><td colspan='5' style='text-align:center;'>No MRSA status data available</td></tr>"
        else:
            for status, count in mrsa_status_dist.most_common():
                if status == 'ND': continue
                pct = (count/total)*100
                sts = set()
                sccs = set()
                for s, d in samples_data.items():
                    if d.get('typing', {}).get('MRSA_Status') == status:
                        mlst = d.get('typing', {}).get('MLST')
                        scc = d.get('typing', {}).get('SCCmec_Type')
                        if mlst and mlst != 'ND':
                            sts.add(mlst)
                        if scc and scc != 'ND' and scc != 'Not Assigned':
                            sccs.add(scc)
                st_list = ', '.join(sorted(sts)) if sts else 'None'
                scc_list = ', '.join(sorted(sccs)) if sccs else 'None'
                badge = '<span class="badge badge-mrsa">MRSA</span>' if 'MRSA' in status else '<span class="badge badge-mssa">MSSA</span>'
                html += f"<tr><td>{badge}</td><td>{count}</td><td>{pct:.1f}%</td><td>{st_list}</td><td>{scc_list}</td></tr>"
        html += "</tbody></table></div>"
        
        # If no MRSA samples, show a message and skip the combination tables
        if not mrsa_samples:
            html += """
            <div class="alert-box alert-info">
                <i class="fas fa-info-circle fa-2x"></i>
                <div>
                    <strong>No MRSA samples detected</strong> – The combination tables below are empty because no methicillin‑resistant isolates were found.
                </div>
            </div>
            """
            return html
        
        # MRSA ST‑spa combinations
        html += f"""
        <h3>🔗 MRSA ST–spa Combinations ({len(mrsa_mlst_spa)} combinations)</h3>
        <input type="text" class="search-box" id="search-mrsa-mlst-spa" onkeyup="searchTable('mrsa-mlst-spa-table','search-mrsa-mlst-spa')" placeholder="🔍 Search ST‑spa...">
        <div class="master-scrollable-container">
            <table id="mrsa-mlst-spa-table" class="data-table">
                <thead><tr><th data-sort="string">ST‑spa Combination</th><th data-sort="number">Count</th><th data-sort="string">Samples</th></tr></thead>
                <tbody>
        """
        for combo, samples in sorted(mrsa_mlst_spa.items(), key=lambda x: len(x[1]), reverse=True):
            sample_tags = ''.join(f'<span class="genome-tag">{s}</span>' for s in samples)
            html += f"<tr><td><strong>{combo}</strong></td><td>{len(samples)}</td><td><div class='genome-list'>{sample_tags}</div></td></tr>"
        html += "</tbody></tr></div>"
        
        # MRSA ST‑SCCmec combinations
        html += f"""
        <h3>🔗 MRSA ST–SCCmec Combinations ({len(mrsa_mlst_sccmec)} combinations)</h3>
        <input type="text" class="search-box" id="search-mrsa-mlst-sccmec" onkeyup="searchTable('mrsa-mlst-sccmec-table','search-mrsa-mlst-sccmec')" placeholder="🔍 Search ST‑SCCmec...">
        <div class="master-scrollable-container">
            <table id="mrsa-mlst-sccmec-table" class="data-table">
                <thead><tr><th data-sort="string">ST‑SCCmec Combination</th><th data-sort="number">Count</th><th data-sort="string">Samples</th></tr></thead>
                <tbody>
        """
        for combo, samples in sorted(mrsa_mlst_sccmec.items(), key=lambda x: len(x[1]), reverse=True):
            sample_tags = ''.join(f'<span class="genome-tag">{s}</span>' for s in samples)
            html += f"<tr><td><strong>{combo}</strong></td><td>{len(samples)}</td><td><div class='genome-list'>{sample_tags}</div></td></tr>"
        html += "</tbody><tr></div>"
        
        # MRSA spa‑SCCmec combinations
        html += f"""
        <h3>🔗 MRSA spa–SCCmec Combinations ({len(mrsa_spa_sccmec)} combinations)</h3>
        <input type="text" class="search-box" id="search-mrsa-spa-sccmec" onkeyup="searchTable('mrsa-spa-sccmec-table','search-mrsa-spa-sccmec')" placeholder="🔍 Search spa‑SCCmec...">
        <div class="master-scrollable-container">
            <table id="mrsa-spa-sccmec-table" class="data-table">
                <thead><tr><th data-sort="string">spa‑SCCmec Combination</th><th data-sort="number">Count</th><th data-sort="string">Samples</th></tr></thead>
                <tbody>
        """
        for combo, samples in sorted(mrsa_spa_sccmec.items(), key=lambda x: len(x[1]), reverse=True):
            sample_tags = ''.join(f'<span class="genome-tag">{s}</span>' for s in samples)
            html += f"<tr><td><strong>{combo}</strong></td><td>{len(samples)}</td><td><div class='genome-list'>{sample_tags}</div></td></tr>"
        html += "</tbody></table></div>"
        
        return html

    def _generate_amr_section(self, kwargs: Dict) -> str:
        """Generate comprehensive AMR genes section with filter buttons and grouping controls"""
        gene_centric = kwargs['gene_centric']
        amr_databases = gene_centric.get('amr_databases', {})
        total_samples = len(kwargs.get('samples_data', {}))
        
        html = """
        <div class="alert-box alert-info">
            <i class="fas fa-biohazard fa-2x"></i>
            <div>
                <h3>🧬 Antimicrobial Resistance (AMR) Genes – Gene‑Centric View + Grouping by Typing</h3>
                <p>Each AMR gene is shown with <strong>all genomes</strong> that contain it. Use the grouping buttons below to reorganise the genome list by MLST, spa, SCCmec, or combinations – this reveals which clones carry specific resistance genes.</p>
                <p>The genome list is <strong>scrollable</strong> and searchable.</p>
            </div>
        </div>
        
        <!-- Grouping controls -->
        <div class="grouping-controls">
            <strong><i class="fas fa-layer-group"></i> Group genomes by:</strong>
            <button class="group-btn" data-group="MLST" onclick="groupGenomesByTyping('amr-table', 'MLST')">MLST (ST)</button>
            <button class="group-btn" data-group="spa" onclick="groupGenomesByTyping('amr-table', 'spa')">spa type</button>
            <button class="group-btn" data-group="SCCmec" onclick="groupGenomesByTyping('amr-table', 'SCCmec')">SCCmec</button>
            <button class="group-btn" data-group="ST-spa" onclick="groupGenomesByTyping('amr-table', 'ST-spa')">ST‑spa</button>
            <button class="group-btn" data-group="ST-SCCmec" onclick="groupGenomesByTyping('amr-table', 'ST-SCCmec')">ST‑SCCmec</button>
            <button class="group-btn" data-group="spa-SCCmec" onclick="groupGenomesByTyping('amr-table', 'spa-SCCmec')">spa‑SCCmec</button>
            <button class="group-btn" data-group="triple" onclick="groupGenomesByTyping('amr-table', 'triple')">Triple (ST‑spa‑SCCmec)</button>
            <button class="group-btn" onclick="resetGenomeList('amr-table')">Reset (flat list)</button>
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
                        <th data-sort="string">Genomes (scrollable, groupable)</th>
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
                percentage = frequency.split('(')[-1].replace(')', '').strip()
            elif count > 0 and total_samples > 0:
                percentage = f"{(count/total_samples)*100:.1f}%"
            
            # Critical gene indicator (warning icon)
            is_critical = any(crit in gene.lower() for crit in self.data_analyzer.critical_amr_genes)
            gene_display = f"<strong>{gene}</strong>" + (" ⚠️" if is_critical else "")
            
            # Build scrollable genome tags (flat list initially)
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
        """Generate comprehensive virulence genes section with grouping controls"""
        gene_centric = kwargs['gene_centric']
        virulence_databases = gene_centric.get('virulence_databases', {})
        total_samples = len(kwargs.get('samples_data', {}))
        
        html = """
        <div class="alert-box alert-info">
            <i class="fas fa-virus fa-2x"></i>
            <div>
                <h3>🧬 Virulence Factors – Gene‑Centric View + Grouping by Typing</h3>
                <p>Each virulence gene is shown with <strong>all genomes</strong> that contain it. Use the grouping buttons to reorganise the genome list by MLST, spa, SCCmec, or combinations – this reveals which clones carry specific virulence factors.</p>
                <p>The genome list is <strong>scrollable</strong> and searchable.</p>
            </div>
        </div>
        
        <!-- Grouping controls -->
        <div class="grouping-controls">
            <strong><i class="fas fa-layer-group"></i> Group genomes by:</strong>
            <button class="group-btn" data-group="MLST" onclick="groupGenomesByTyping('vir-table', 'MLST')">MLST (ST)</button>
            <button class="group-btn" data-group="spa" onclick="groupGenomesByTyping('vir-table', 'spa')">spa type</button>
            <button class="group-btn" data-group="SCCmec" onclick="groupGenomesByTyping('vir-table', 'SCCmec')">SCCmec</button>
            <button class="group-btn" data-group="ST-spa" onclick="groupGenomesByTyping('vir-table', 'ST-spa')">ST‑spa</button>
            <button class="group-btn" data-group="ST-SCCmec" onclick="groupGenomesByTyping('vir-table', 'ST-SCCmec')">ST‑SCCmec</button>
            <button class="group-btn" data-group="spa-SCCmec" onclick="groupGenomesByTyping('vir-table', 'spa-SCCmec')">spa‑SCCmec</button>
            <button class="group-btn" data-group="triple" onclick="groupGenomesByTyping('vir-table', 'triple')">Triple (ST‑spa‑SCCmec)</button>
            <button class="group-btn" onclick="resetGenomeList('vir-table')">Reset (flat list)</button>
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
                        <th data-sort="string">Genomes (scrollable, groupable)</th>
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
            
            # Build scrollable genome tags (flat list initially)
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
        """Generate comprehensive BACMET section with grouping controls"""
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
                <h3>🧪 BACMET: Biocide & Heavy Metal Resistance – Groupable by Typing</h3>
                <p>BACMET2 genes provide resistance to hospital disinfectants (e.g., quaternary ammonium compounds – <em>qac</em> genes) and heavy metals (<em>mer</em> for mercury, <em>ars</em> for arsenic, <em>cop</em> for copper). These markers can co‑select for antibiotic resistance and promote persistence in healthcare settings.</p>
                <p>Each gene is shown with <strong>all genomes</strong> that carry it. Use the grouping buttons to reorganise by typing.</p>
                <p>The genome list is <strong>scrollable</strong> – all samples are visible.</p>
            </div>
        </div>
        
        <!-- Grouping controls -->
        <div class="grouping-controls">
            <strong><i class="fas fa-layer-group"></i> Group genomes by:</strong>
            <button class="group-btn" data-group="MLST" onclick="groupGenomesByTyping('bac-table', 'MLST')">MLST (ST)</button>
            <button class="group-btn" data-group="spa" onclick="groupGenomesByTyping('bac-table', 'spa')">spa type</button>
            <button class="group-btn" data-group="SCCmec" onclick="groupGenomesByTyping('bac-table', 'SCCmec')">SCCmec</button>
            <button class="group-btn" data-group="ST-spa" onclick="groupGenomesByTyping('bac-table', 'ST-spa')">ST‑spa</button>
            <button class="group-btn" data-group="ST-SCCmec" onclick="groupGenomesByTyping('bac-table', 'ST-SCCmec')">ST‑SCCmec</button>
            <button class="group-btn" data-group="spa-SCCmec" onclick="groupGenomesByTyping('bac-table', 'spa-SCCmec')">spa‑SCCmec</button>
            <button class="group-btn" data-group="triple" onclick="groupGenomesByTyping('bac-table', 'triple')">Triple (ST‑spa‑SCCmec)</button>
            <button class="group-btn" onclick="resetGenomeList('bac-table')">Reset (flat list)</button>
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
                        <th data-sort="string">Genomes (scrollable, groupable)</th>
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
        """Generate plasmid replicon section with grouping controls"""
        gene_centric = kwargs['gene_centric']
        plasmid_databases = gene_centric.get('plasmid_databases', {})
        total_samples = len(kwargs.get('samples_data', {}))
        
        html = """
        <div class="alert-box alert-info">
            <i class="fas fa-plug fa-2x"></i>
            <div>
                <h3>🧬 Plasmid Replicons – Horizontal Gene Transfer + Grouping by Typing</h3>
                <p>Plasmid replicons indicate the presence of specific plasmid types. Use grouping buttons to see which clones carry which plasmids.</p>
                <p>Each replicon is shown with <strong>all genomes</strong> that carry it (scrollable list).</p>
            </div>
        </div>
        
        <!-- Grouping controls -->
        <div class="grouping-controls">
            <strong><i class="fas fa-layer-group"></i> Group genomes by:</strong>
            <button class="group-btn" data-group="MLST" onclick="groupGenomesByTyping('plasmid-table', 'MLST')">MLST (ST)</button>
            <button class="group-btn" data-group="spa" onclick="groupGenomesByTyping('plasmid-table', 'spa')">spa type</button>
            <button class="group-btn" data-group="SCCmec" onclick="groupGenomesByTyping('plasmid-table', 'SCCmec')">SCCmec</button>
            <button class="group-btn" data-group="ST-spa" onclick="groupGenomesByTyping('plasmid-table', 'ST-spa')">ST‑spa</button>
            <button class="group-btn" data-group="ST-SCCmec" onclick="groupGenomesByTyping('plasmid-table', 'ST-SCCmec')">ST‑SCCmec</button>
            <button class="group-btn" data-group="spa-SCCmec" onclick="groupGenomesByTyping('plasmid-table', 'spa-SCCmec')">spa‑SCCmec</button>
            <button class="group-btn" data-group="triple" onclick="groupGenomesByTyping('plasmid-table', 'triple')">Triple (ST‑spa‑SCCmec)</button>
            <button class="group-btn" onclick="resetGenomeList('plasmid-table')">Reset (flat list)</button>
        </div>
        
        <input type="text" class="search-box" id="search-plasmid" onkeyup="searchTable('plasmid-table','search-plasmid')" placeholder="🔍 Search plasmids...">
        <input type="text" class="search-box" id="search-plasmid-genome" onkeyup="highlightGenome('plasmid-table','search-plasmid-genome')" placeholder="🔍 Highlight genomes...">
        
        <div class="action-buttons">
            <button class="action-btn btn-primary" onclick="exportTableToCSV('plasmid-table', 'plasmid_replicons.csv')"><i class="fas fa-download"></i> Export</button>
        </div>
        
        <h3>📋 Plasmid Replicons Detected</h3>
        <div class="master-scrollable-container"><table id="plasmid-table" class="data-table"><thead><tr><th data-sort="string">Replicon</th><th data-sort="string">Database</th><th data-sort="number">Count</th><th data-sort="number">Percentage</th><th data-sort="string">Genomes (scrollable, groupable)</th></tr></thead><tbody>
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
    
    def _generate_mutation_section(self, kwargs: Dict) -> str:
        integrated_data = kwargs.get('integrated_data', {})
        mutation_data = integrated_data.get('mutation_data', {})
        mutations = mutation_data.get('mutations', [])
        total_samples = len(kwargs.get('samples_data', {}))
        if not mutations:
            return """<div class="alert-box alert-warning"><i class="fas fa-dna fa-2x"></i><div><h3>No Mutation Data Available</h3><p>The <code>mutation_summary.html</code> file was not found or could not be parsed. Make sure you have run the AMRfinderPlus module with mutation reporting enabled (default) and that <code>mutation_summary.html</code> is present in the input directory.</p></div></div>"""
        html = """
        <div class="alert-box alert-info"><i class="fas fa-dna fa-2x"></i><div><h3>🧬 Point Mutations – Gene‑Centric View + Grouping by Typing</h3><p>Each unique mutation (gene + element name) is shown with <strong>all genomes</strong> that carry it. Mutations in key genes (e.g., <em>gyrA, parC, rpoB, 23S rRNA</em>) can confer resistance even without acquired resistance genes. Use grouping buttons to reorganise by typing.</p></div></div>
        <div class="grouping-controls"><strong><i class="fas fa-layer-group"></i> Group genomes by:</strong><button class="group-btn" data-group="MLST" onclick="groupGenomesByTyping('mutation-table', 'MLST')">MLST (ST)</button><button class="group-btn" data-group="spa" onclick="groupGenomesByTyping('mutation-table', 'spa')">spa type</button><button class="group-btn" data-group="SCCmec" onclick="groupGenomesByTyping('mutation-table', 'SCCmec')">SCCmec</button><button class="group-btn" data-group="ST-spa" onclick="groupGenomesByTyping('mutation-table', 'ST-spa')">ST‑spa</button><button class="group-btn" data-group="ST-SCCmec" onclick="groupGenomesByTyping('mutation-table', 'ST-SCCmec')">ST‑SCCmec</button><button class="group-btn" data-group="spa-SCCmec" onclick="groupGenomesByTyping('mutation-table', 'spa-SCCmec')">spa‑SCCmec</button><button class="group-btn" data-group="triple" onclick="groupGenomesByTyping('mutation-table', 'triple')">Triple (ST‑spa‑SCCmec)</button><button class="group-btn" onclick="resetGenomeList('mutation-table')">Reset (flat list)</button></div>
        <input type="text" class="search-box" id="search-mutation" onkeyup="searchTable('mutation-table','search-mutation')" placeholder="🔍 Search mutation by gene or mutation name...">
        <input type="text" class="search-box" id="search-mutation-genome" onkeyup="highlightGenome('mutation-table','search-mutation-genome')" placeholder="🔍 Highlight genomes containing specific text (e.g., sample ID)">
        <div class="action-buttons"><button class="action-btn btn-primary" onclick="exportTableToCSV('mutation-table', 'mutations.csv')"><i class="fas fa-download"></i> Export All Mutations</button><button class="action-btn btn-danger" onclick="document.getElementById('search-mutation').value='LINEZOLID'; searchTable('mutation-table','search-mutation')"><i class="fas fa-skull-crossbones"></i> Linezolid‑related</button><button class="action-btn btn-warning" onclick="document.getElementById('search-mutation').value='QUINOLONE'; searchTable('mutation-table','search-mutation')"><i class="fas fa-biohazard"></i> Quinolone‑related</button><button class="action-btn btn-warning" onclick="document.getElementById('search-mutation').value='RIFAMPIN'; searchTable('mutation-table','search-mutation')"><i class="fas fa-biohazard"></i> Rifampin‑related</button><button class="action-btn btn-info" onclick="document.getElementById('search-mutation').value='gyrA'; searchTable('mutation-table','search-mutation')"><i class="fas fa-dna"></i> gyrA</button><button class="action-btn btn-info" onclick="document.getElementById('search-mutation').value='parC'; searchTable('mutation-table','search-mutation')"><i class="fas fa-dna"></i> parC</button><button class="action-btn btn-info" onclick="document.getElementById('search-mutation').value='rpoB'; searchTable('mutation-table','search-mutation')"><i class="fas fa-dna"></i> rpoB</button><button class="action-btn btn-info" onclick="document.getElementById('search-mutation').value='23S'; searchTable('mutation-table','search-mutation')"><i class="fas fa-dna"></i> 23S rRNA</button><button class="action-btn btn-light" onclick="document.getElementById('search-mutation').value=''; searchTable('mutation-table','search-mutation')"><i class="fas fa-sync"></i> Clear Search</button></div>
        <div style="margin: 10px 0 20px 0; background: #f8f9fa; padding: 15px; border-radius: 8px; font-size: 0.9em; border-left: 4px solid #00BCD4;"><strong><i class="fas fa-info-circle"></i> Clinical relevance of key mutations:</strong><br>• <strong>23S rRNA (linezolid)</strong> – Mutations (e.g., G2576T, T2500A) confer resistance to linezolid, a last‑line antibiotic for MRSA/VRE.<br>• <strong>gyrA/parC (quinolones)</strong> – Mutations in the quinolone resistance‑determining region (QRDR) reduce susceptibility to fluoroquinolones (ciprofloxacin, levofloxacin).<br>• <strong>rpoB (rifampin)</strong> – Mutations cause high‑level rifampin resistance, often used in combination therapy.<br>• <strong>mprF (daptomycin)</strong> – Mutations can lead to daptomycin non‑susceptibility.<br>• <strong>rplC/rplD (linezolid)</strong> – Ribosomal protein mutations also confer linezolid resistance.</div>
        <h3><i class="fas fa-dna"></i> All Detected Point Mutations</h3><div class="master-scrollable-container"><table id="mutation-table" class="data-table"><thead><tr>
            <th data-sort="string">Gene</th>
            <th data-sort="string">Mutation (Element name)</th>
            <th data-sort="string">Class</th>
            <th data-sort="string">Subclass</th>
            <th data-sort="number">Genome Count</th>
            <th data-sort="string">Genomes (scrollable, groupable)</th>
        </tr></thead><tbody>"""
        for mut in mutations:
            gene = mut['gene']
            mutation = mut['mutation']
            class_name = mut['class']
            subclass = mut['subclass']
            genomes = mut['genomes']
            # Use actual number of distinct genomes
            count = len(genomes)
            mutation_display = mutation[:80] + '…' if len(mutation) > 80 else mutation
            genome_tags = ''.join(f'<span class="genome-tag">{g}</span>' for g in genomes)
            html += f"""
                <tr>
                    <td><strong>{gene}</strong></td>
                    <td><span title='{mutation}'>{mutation_display}</span></td>
                    <td>{class_name}</td>
                    <td>{subclass}</td>
                    <td><strong>{count}</strong> ({count/total_samples*100:.1f}%)</td>
                    <td><div class="genome-list">{genome_tags}</div></td>
                </tr>"""
        html += """
                </tbody>
            </table>
        </div>
        """
        return html
    
    def _generate_pattern_discovery_section(self, kwargs: Dict) -> str:
        patterns = kwargs['patterns']
        total_samples = len(kwargs.get('samples_data', {}))
        triple_combos = patterns.get('triple_combinations', {})
        html = f"""<div class="alert-box alert-info"><i class="fas fa-project-diagram fa-2x"></i><div><h3>🔍 Cross‑Genome Pattern Discovery</h3><p>This tab reveals associations between typing results, gene co‑occurrence, and high‑risk combinations. <strong>Triple Typing Combination (MLST – spa – SCCmec)</strong> – the most informative molecular fingerprint for S. aureus epidemiology.</p></div></div>
        <h3>🔗 Triple Typing Combination (ST – spa – SCCmec) – Outbreak Signature</h3><input type="text" class="search-box" id="search-triple" onkeyup="searchTable('triple-table','search-triple')" placeholder="🔍 Search combination..."><div class="master-scrollable-container"><table id="triple-table" class="data-table"><thead><tr><th data-sort="string">ST - spa - SCCmec</th><th data-sort="number">Count</th><th data-sort="string">Samples</th></tr></thead><tbody>"""
        for combo, samples in sorted(triple_combos.items(), key=lambda x: len(x[1]), reverse=True):
            sample_tags = ''.join(f'<span class="genome-tag">{s}</span>' for s in samples)
            html += f"<tr><td><strong>{combo}</strong></td><td>{len(samples)}</td><td><div class='genome-list'>{sample_tags}</div></td></tr>"
        html += "</tbody></table></div>"
        high_risk = patterns.get('high_risk_combinations', [])
        if high_risk:
            html += f"""<h3>⚠️ High‑Risk Combinations (Critical AMR + Critical Virulence)</h3><div class="alert-box alert-danger"><i class="fas fa-radiation"></i><div><strong>{len(high_risk)} samples</strong> carry both critical AMR and virulence genes.</div></div><div class="master-scrollable-container"><table id="highrisk-table" class="data-table"><thead><tr><th data-sort="string">Sample</th><th data-sort="string">MLST</th><th data-sort="string">spa</th><th data-sort="string">SCCmec</th><th data-sort="string">Critical AMR</th><th data-sort="string">Critical Virulence</th></tr></thead><tbody>"""
            for c in high_risk:
                html += f"<tr><td><strong>{c['sample']}</strong></td><td>{c['mlst']}</td><td>{c['spa_type']}</td><td>{c['sccmec_type']}</td><td>{', '.join(c['critical_amr_genes'])}</td><td>{', '.join(c['critical_virulence_genes'])}</td></tr>"
            html += "</tbody></table></div>"
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
        return """
        <div class="alert-box alert-info"><i class="fas fa-robot fa-2x"></i><div><h3>🤖 AI Assistant Guide – How to Use with ChatGPT, Claude, or Gemini</h3><p>You can upload this HTML report to AI tools to ask detailed questions about your S. aureus dataset.</p></div></div>
        <div style="margin: 20px 0;"><div class="database-section"><h4><i class="fas fa-chart-line"></i> Example Questions</h4><ul><li>What are the most common MLST sequence types in this dataset?</li><li>Which spa types are dominant in MRSA vs MSSA?</li><li>How many samples carry mecA? What are their STs and spa types?</li><li>Are there any samples with vancomycin resistance genes (vanA/vanB)?</li><li>Which samples carry the PVL toxin genes (lukF/S-PV)?</li><li>List all samples with TSST-1 (tsst).</li><li>Show me the triple typing combinations (ST-spa-SCCmec) with more than 2 isolates.</li><li>Which resistance genes co-occur most frequently?</li><li>What are the most frequent point mutations in gyrA or parC?</li></ul></div><div class="database-section"><h4><i class="fas fa-upload"></i> How to Upload Data</h4><p>Save the <strong>staphscope_ultimate_report.json</strong> file and upload it to the AI. Most AI tools accept file uploads, or you can copy‑paste relevant tables.</p></div><div class="database-section"><h4><i class="fas fa-balance-scale"></i> Ethical Guidelines for AI Use</h4><ul><li><strong>AI assists, but experts decide:</strong> Always interpret AI-generated insights in the context of your local epidemiology and clinical knowledge.</li><li><strong>No patient data:</strong> Do not upload raw patient data or identifiable information. Only use aggregated, de‑identified genomic data.</li><li><strong>Verify critical findings:</strong> Important resistance/virulence calls should be confirmed by human review or secondary tools.</li><li><strong>Stay transparent:</strong> When publishing results, disclose that AI tools were used for exploratory analysis.</li></ul><p><i class="fas fa-smile-wink"></i> Remember: AI won’t take your job – but a microbiologist who knows how to use AI might! 😉</p></div></div>
        """
    
    def _generate_citation_section(self, kwargs: Dict) -> str:
        return """
        <div class="alert-box alert-info"><i class="fas fa-quote-right fa-2x"></i><div><h3>📚 How to Cite StaphScope and Its Dependencies</h3><p>If you use StaphScope in your research, please cite the main tool and the relevant third‑party tools and databases.</p></div></div>
        <div class="accordion"><div class="accordion-item"><div class="accordion-header"><span>📄 From the ESKAPE AMR Platform</span><i class="fas fa-chevron-down"></i></div><div class="accordion-content" style="display: none;"><ul class="citation-list"><li><strong>StaphScope</strong> – Beckley B, Amarh V. StaphScope: a species‑optimized computational pipeline for rapid and accessible <em>Staphylococcus aureus</em> genotyping and surveillance. <em>BMC Genomics</em>. 2026;27:261. doi:10.1186/s12864-026-12609-x<button class="copy-btn" data-citation='Beckley B, Amarh V. StaphScope: a species‑optimized computational pipeline for rapid and accessible Staphylococcus aureus genotyping and surveillance. BMC Genomics. 2026;27:261. doi:10.1186/s12864-026-12609-x'>📋 Copy citation</button></li></ul></div></div><div class="accordion-item"><div class="accordion-header"><span>🔧 Key Databases & Methods</span><i class="fas fa-chevron-down"></i></div><div class="accordion-content" style="display: none;"><ul class="citation-list"><li><strong>MLST</strong> – Seemann T. MLST: Scan contig files against PubMLST typing schemes. GitHub. 2018. https://github.com/tseemann/mlst<button class="copy-btn" data-citation='Seemann T. MLST: Scan contig files against PubMLST typing schemes. GitHub. 2018. https://github.com/tseemann/mlst'>📋 Copy citation</button></li><li><strong>PubMLST / BIGSdb</strong> – Jolley KA, Bray JE, Maiden MCJ. Open‑access bacterial population genomics: BIGSdb software, the PubMLST.org website and their applications. <em>Wellcome Open Res</em>. 2018;3:124. doi:10.12688/wellcomeopenres.14826.1<button class="copy-btn" data-citation='Jolley KA, Bray JE, Maiden MCJ. Open-access bacterial population genomics: BIGSdb software, the PubMLST.org website and their applications. Wellcome Open Res. 2018;3:124. doi:10.12688/wellcomeopenres.14826.1'>📋 Copy citation</button></li><li><strong>spa typing</strong> – Harmsen D, et al. Typing of methicillin‑resistant <em>Staphylococcus aureus</em> in a university hospital setting by using novel software for spa repeat determination and database management. <em>J Clin Microbiol</em>. 2003;41(12):5442-8. doi:10.1128/JCM.41.12.5442-5448.2003<button class="copy-btn" data-citation='Harmsen D, et al. Typing of methicillin-resistant Staphylococcus aureus in a university hospital setting by using novel software for spa repeat determination and database management. J Clin Microbiol. 2003;41(12):5442-8. doi:10.1128/JCM.41.12.5442-5448.2003'>📋 Copy citation</button></li><li><strong>SCCmecFinder</strong> – Kaya H, et al. SCCmecFinder, a Web‑Based Tool for Typing of Staphylococcal Cassette Chromosome mec in <em>Staphylococcus aureus</em> Using Whole‑Genome Sequence Data. <em>mSphere</em>. 2018;3(1):e00612-17. doi:10.1128/mSphere.00612-17<button class="copy-btn" data-citation='Kaya H, et al. SCCmecFinder, a Web‑Based Tool for Typing of Staphylococcal Cassette Chromosome mec in Staphylococcus aureus Using Whole‑Genome Sequence Data. mSphere. 2018;3(1):e00612-17. doi:10.1128/mSphere.00612-17'>📋 Copy citation</button></li><li><strong>AMRFinderPlus</strong> – Feldgarden M, et al. AMRFinderPlus and the Reference Gene Catalog facilitate examination of the genomic links among antimicrobial resistance, stress response, and virulence. <em>Sci Rep</em>. 2021;11(1):12728. doi:10.1038/s41598-021-91456-0<button class="copy-btn" data-citation='Feldgarden M, et al. AMRFinderPlus and the Reference Gene Catalog facilitate examination of the genomic links among antimicrobial resistance, stress response, and virulence. Sci Rep. 2021;11(1):12728. doi:10.1038/s41598-021-91456-0'>📋 Copy citation</button></li><li><strong>ABRicate</strong> – Seemann T. ABRicate: mass screening of contigs for antibiotic resistance genes. GitHub. 2024. https://github.com/tseemann/abricate<button class="copy-btn" data-citation='Seemann T. ABRicate: mass screening of contigs for antibiotic resistance genes. GitHub. 2024. https://github.com/tseemann/abricate'>📋 Copy citation</button></li><li><strong>CARD</strong> – McArthur AG, et al. The comprehensive antibiotic resistance database. <em>Antimicrob Agents Chemother</em>. 2013;57(7):3348-57. doi:10.1128/AAC.00419-13<button class="copy-btn" data-citation='McArthur AG, et al. The comprehensive antibiotic resistance database. Antimicrob Agents Chemother. 2013;57(7):3348-57. doi:10.1128/AAC.00419-13'>📋 Copy citation</button></li><li><strong>ResFinder</strong> – Florensa AF, et al. ResFinder – an open online resource for identification of antimicrobial resistance genes in next‑generation sequencing data and prediction of phenotypes from genotypes. <em>Microb Genom</em>. 2022;8(1):000748. doi:10.1099/mgen.0.000748<button class="copy-btn" data-citation='Florensa AF, et al. ResFinder – an open online resource for identification of antimicrobial resistance genes in next-generation sequencing data and prediction of phenotypes from genotypes. Microb Genom. 2022;8(1):000748. doi:10.1099/mgen.0.000748'>📋 Copy citation</button></li><li><strong>VFDB</strong> – Chen L, et al. VFDB 2012 update: toward the genetic diversity and molecular evolution of bacterial virulence factors. <em>Nucleic Acids Res</em>. 2012;40(Database issue):D641-5. doi:10.1093/nar/gkr989<button class="copy-btn" data-citation='Chen L, et al. VFDB 2012 update: toward the genetic diversity and molecular evolution of bacterial virulence factors. Nucleic Acids Res. 2012;40(Database issue):D641-5. doi:10.1093/nar/gkr989'>📋 Copy citation</button></li><li><strong>PlasmidFinder</strong> – Carattoli A, et al. <em>In silico</em> detection and typing of plasmids using PlasmidFinder and plasmid multilocus sequence typing. <em>Antimicrob Agents Chemother</em>. 2014;58(7):3895-903. doi:10.1128/AAC.02412-14<button class="copy-btn" data-citation='Carattoli A, et al. In silico detection and typing of plasmids using PlasmidFinder and plasmid multilocus sequence typing. Antimicrob Agents Chemother. 2014;58(7):3895-903. doi:10.1128/AAC.02412-14'>📋 Copy citation</button></li><li><strong>BacMet</strong> – Pal C, et al. BacMet: antibacterial biocide and metal resistance genes database. <em>Nucleic Acids Res</em>. 2014;42(Database issue):D737-43. doi:10.1093/nar/gkt1252<button class="copy-btn" data-citation='Pal C, et al. BacMet: antibacterial biocide and metal resistance genes database. Nucleic Acids Res. 2014;42(Database issue):D737-43. doi:10.1093/nar/gkt1252'>📋 Copy citation</button></li><li><strong>MEGARes</strong> – Doster E, et al. MEGARes 2.0: a database for classification of antimicrobial drug, biocide and metal resistance determinants in metagenomic sequence data. <em>Nucleic Acids Res</em>. 2020;48(D1):D561-D569. doi:10.1093/nar/gkz1010<button class="copy-btn" data-citation='Doster E, et al. MEGARes 2.0: a database for classification of antimicrobial drug, biocide and metal resistance determinants in metagenomic sequence data. Nucleic Acids Res. 2020;48(D1):D561-D569. doi:10.1093/nar/gkz1010'>📋 Copy citation</button></li><li><strong>ARG-ANNOT</strong> – Gupta SK, et al. ARG-ANNOT, a new bioinformatic tool to discover antibiotic resistance genes in bacterial genomes. <em>Antimicrob Agents Chemother</em>. 2014;58(1):212-20. doi:10.1128/AAC.01310-13<button class="copy-btn" data-citation='Gupta SK, et al. ARG-ANNOT, a new bioinformatic tool to discover antibiotic resistance genes in bacterial genomes. Antimicrob Agents Chemother. 2014;58(1):212-20. doi:10.1128/AAC.01310-13'>📋 Copy citation</button></li><li><strong>Biopython</strong> – Cock PJ, et al. Biopython: freely available Python tools for computational molecular biology and bioinformatics. <em>Bioinformatics</em>. 2009;25(11):1422-3. doi:10.1093/bioinformatics/btp163<button class="copy-btn" data-citation='Cock PJ, et al. Biopython: freely available Python tools for computational molecular biology and bioinformatics. Bioinformatics. 2009;25(11):1422-3. doi:10.1093/bioinformatics/btp163'>📋 Copy citation</button></li></ul></div></div></div>
        <div class="alert-box alert-success" style="margin-top: 20px;"><i class="fas fa-hand-peace"></i><div><strong>Suggested acknowledgement:</strong><br>"Genomic analysis was performed using StaphScope [Beckley &amp; Amarh, 2026], which integrates MLST [Seemann, 2018] using the PubMLST database [Jolley et al., 2018], ABRicate [Seemann, 2018], AMRFinderPlus [Feldgarden et al., 2021], and SCCmecFinder [Kaya et al., 2018] for comprehensive <em>S. aureus</em> characterization. Antimicrobial resistance genes were identified using the CARD [McArthur et al., 2013] and ResFinder [Florensa et al., 2022] databases. For biocide and heavy metal resistance genes, BacMet [Pal et al., 2014] was used. Virulence and plasmid screening were performed with ABRicate using the VFDB [Chen et al., 2012] and PlasmidFinder [Carattoli et al., 2014] databases."</div></div>
        """
    
    def _generate_funding_section(self, kwargs: Dict) -> str:
        return """
        <div class="alert-box alert-info"><i class="fas fa-coffee fa-2x"></i><div><h3>☕ Funding & Support – Keeping the Lights On (with code and caffeine)</h3><p>StaphScope is an <strong>independent, unfunded project</strong> born out of passion for genomic surveillance and AMR research at the University of Ghana Medical School.</p><p>No grants, no sponsors, no institutional backing – just a laptop, a lot of coffee, and a burning desire to help researchers fight antimicrobial resistance.</p></div></div>
        <div class="alert-box alert-warning"><i class="fas fa-heart fa-2x"></i><div><h3>💡 How You Can Help (Without Opening Your Wallet)</h3><ul><li><strong>⭐ Star us on GitHub</strong> – It takes two seconds and makes us feel like rockstars.</li><li><strong>🐛 Report bugs</strong> – If something breaks, let us know. We’ll fix it with joy.</li><li><strong>💡 Suggest features</strong> – Have an idea? We’re all ears (and we actually implement them, as you’ve seen!).</li><li><strong>🧬 Share your data</strong> – If you’ve used StaphScope and want to collaborate, we’d love to hear your story.</li><li><strong>📢 Spread the word</strong> – Tell your colleagues, tweet about it, or mention it in your next Zoom call.</li><li><strong>👋 Say hello!</strong> – Seriously, just drop an email to <strong>brownbeckley94@gmail.com</strong>. It makes our day (Literally My Day).</li></ul><p><i class="fas fa-microbe"></i> <strong>Fun fact:</strong> This project runs on 100% volunteer tears, 0% grant money. But we’re not bitter – we’re just caffeinated.</p></div></div>
        <div class="alert-box alert-success"><i class="fas fa-hand-holding-heart"></i><div><h3>🤝 Contribute to the ESKAPE AMR Platform</h3><p>We also maintain pipelines for other ESKAPE pathogens (AcinetoScope, Kleboscope, Pseudoscope, etc.). If you’re a developer, bioinformatician, or just someone who loves clean code and bacteria, we welcome: pull requests, issues, documentation improvements, ideas for new databases. Visit our GitHub: <a href="https://github.com/bbeckley-hub" target="_blank">https://github.com/bbeckley-hub</a> – star, fork, and let’s fight AMR together!</p><p><strong>Brown Beckley</strong> – <i class="fas fa-envelope"></i> brownbeckley94@gmail.com</p><p><i class="fas fa-laugh-beam"></i> <strong>P.S.</strong> If you ever meet Brown in person, buy him a coffee or I'll buy you a bug myself. He’ll probably talk your ear off about SCCmec types, but it’s worth it.</p></div></div>
        """
    
    def _generate_export_section(self, kwargs: Dict) -> str:
        return """
        <div class="alert-box alert-info"><i class="fas fa-download fa-2x"></i><div><h3>📥 Export Data – Download Your Results</h3><p>Download all analysis tables in CSV format or the complete JSON data for downstream use.</p></div></div>
        <div style="display: grid; grid-template-columns: repeat(auto-fit, minmax(280px, 1fr)); gap: 20px; margin: 30px 0;"><div class="dashboard-card card-export" onclick="exportTableToCSV('samples-table', 'sample_overview.csv')"><div style="font-size: 2em; color: var(--export-color);"><i class="fas fa-table"></i></div><div class="card-label">Sample Overview CSV</div><p>All samples with MLST, spa, SCCmec, MRSA status</p></div><div class="dashboard-card card-export" onclick="exportTableToCSV('amr-table', 'amr_genes.csv')"><div style="font-size: 2em; color: var(--export-color);"><i class="fas fa-biohazard"></i></div><div class="card-label">AMR Genes CSV</div><p>Gene‑centric AMR table with genome lists</p></div><div class="dashboard-card card-export" onclick="exportTableToCSV('vir-table', 'virulence_genes.csv')"><div style="font-size: 2em; color: var(--export-color);"><i class="fas fa-virus"></i></div><div class="card-label">Virulence Genes CSV</div><p>Gene‑centric virulence table</p></div><div class="dashboard-card card-export" onclick="exportTableToCSV('bac-table', 'bacmet_genes.csv')"><div style="font-size: 2em; color: var(--export-color);"><i class="fas fa-flask"></i></div><div class="card-label">BACMET Genes CSV</div><p>Biocide & heavy metal resistance genes</p></div><div class="dashboard-card card-export" onclick="exportTableToCSV('plasmid-table', 'plasmid_replicons.csv')"><div style="font-size: 2em; color: var(--export-color);"><i class="fas fa-plug"></i></div><div class="card-label">Plasmid Replicons CSV</div><p>Plasmid replicon distribution</p></div><div class="dashboard-card card-export" onclick="exportTableToCSV('mutation-table', 'mutations.csv')"><div style="font-size: 2em; color: var(--export-color);"><i class="fas fa-dna"></i></div><div class="card-label">Mutations CSV</div><p>Point mutations table</p></div><div class="dashboard-card card-export" onclick="exportTableToCSV('qc-table', 'fasta_qc.csv')"><div style="font-size: 2em; color: var(--export-color);"><i class="fas fa-chart-line"></i></div><div class="card-label">FASTA QC CSV</div><p>Assembly quality metrics</p></div><div class="dashboard-card card-export" onclick="location.href='staphscope_ultimate_report.json'"><div style="font-size: 2em; color: var(--export-color);"><i class="fas fa-file-code"></i></div><div class="card-label">Complete JSON Data</div><p>All analysis data in structured JSON format</p></div></div>
        <div class="alert-box alert-warning"><i class="fas fa-save"></i><div><strong>Note:</strong> The JSON file (<code>staphscope_ultimate_report.json</code>) is saved in the same directory as this HTML report. You can use it for custom scripts or upload to AI tools.</div></div>
        """


class StaphUltimateReporter:
    def __init__(self, input_dir: Path):
        self.input_dir = Path(input_dir)
        self.output_dir = self.input_dir / "STAPHSCOPE_ULTIMATE_REPORTS"
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.parser = StaphHTMLParser()
        self.analyzer = StaphDataAnalyzer()
        self.html_generator = StaphHTMLGenerator(self.analyzer)
        self.metadata = {
            "tool_name": "STAPHSCOPE Ultimate S. aureus Reporter",
            "version": "2.2.0",
            "author": "Brown Beckley <brownbeckley94@gmail.com>",
            "affiliation": "University of Ghana Medical School",
            "analysis_date": datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            "input_directory": str(self.input_dir)
        }
    
    def find_html_files(self) -> Dict[str, List[Path]]:
        print("🔍 Searching for StaphScope HTML reports...")
        html_files = {'comprehensive': [], 'amrfinder': [], 'abricate': [], 'qc': []}
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
        print("\n🔗 Integrating data from all reports...")
        integrated_data = {'metadata': self.metadata, 'samples': {}, 'patterns': {}, 'gene_centric': {}, 'qc_data': {}}
        if html_files['qc']:
            integrated_data['qc_data'] = self.parser.parse_qc_report(html_files['qc'][0])
            print(f"  ✅ QC data parsed for {len(integrated_data['qc_data'])} samples")
        # Parse mutation summary HTML
        mutation_data = {}
        mutation_html = self.input_dir / "mutation_summary.html"
        if not mutation_html.exists():
            mutation_html = self.input_dir / "staph_amrfinder_results" / "mutation_summary.html"
        if mutation_html.exists():
            mutation_data = self.parser.parse_mutation_summary_html(mutation_html)
        else:
            print("  ⚠️ mutation_summary.html not found; mutation tab will be empty")
        integrated_data['mutation_data'] = mutation_data
        typing_data = {}
        if html_files['comprehensive']:
            typing_data = self.parser.parse_comprehensive_report(html_files['comprehensive'][0])
        amr_by_sample, amr_gene_freq = {}, {}
        if html_files['amrfinder']:
            amr_by_sample, amr_gene_freq = self.parser.parse_amrfinder_report(html_files['amrfinder'][0])
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
        for sample in all_samples:
            typing_info = typing_data.get(sample, {'MLST': 'ND', 'spa_Type': 'ND', 'SCCmec_Type': 'ND', 'MRSA_Status': 'ND'})
            amr_info = amr_by_sample.get(sample, {'critical_genes': [], 'high_risk_genes': [], 'all_genes': []})
            abricate_info = abricate_by_sample.get(sample, {})
            integrated_data['samples'][sample] = {'typing': typing_info, 'amrfinder': amr_info, 'abricate_databases': abricate_info}
        integrated_data['gene_frequencies'] = {'amrfinder': amr_gene_freq, 'abricate': abricate_gene_freq}
        print("\n🧠 Processing gene-centric analysis...")
        integrated_data['gene_centric'] = self.analyzer.create_gene_centric_tables(integrated_data)
        integrated_data['patterns'] = self.analyzer.create_cross_genome_patterns(integrated_data)
        return integrated_data
    
    def generate_json_report(self, integrated_data: Dict[str, Any]) -> Path:
        print("\n📝 Generating JSON report...")
        output_file = self.output_dir / "staphscope_ultimate_report.json"
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
        print("\n📊 Generating CSV reports...")
        samples_data = []
        for sample, data in integrated_data['samples'].items():
            samples_data.append({
                'Sample': sample,
                'MLST': data['typing']['MLST'],
                'spa_Type': data['typing']['spa_Type'],
                'SCCmec_Type': data['typing']['SCCmec_Type'],
                'MRSA_Status': data['typing']['MRSA_Status'],
                'Critical_AMR_Genes': ';'.join(data['amrfinder']['critical_genes']),
                'High_Risk_AMR_Genes': ';'.join(data['amrfinder']['high_risk_genes']),
                'Total_AMR_Genes': len(data['amrfinder']['all_genes']),
                'VFDB_Genes': ';'.join(data.get('abricate_databases', {}).get('vfdb', []))
            })
        pd.DataFrame(samples_data).to_csv(self.output_dir / "sample_overview.csv", index=False)
        gene_centric = integrated_data.get('gene_centric', {})
        amr_data = []
        for db_name, genes in gene_centric.get('amr_databases', {}).items():
            for gene_info in genes:
                amr_data.append({'Gene': gene_info['gene'], 'Database': gene_info['database'], 'Count': gene_info['count'], 'Frequency': gene_info['frequency'], 'Percentage': f"{(gene_info['count']/len(integrated_data['samples']))*100:.1f}%" if len(integrated_data['samples']) > 0 else "0%", 'Genomes': ';'.join(gene_info.get('genomes', []))})
        if amr_data:
            pd.DataFrame(amr_data).to_csv(self.output_dir / "amr_genes.csv", index=False)
        virulence_data = []
        for db_name, genes in gene_centric.get('virulence_databases', {}).items():
            for gene_info in genes:
                virulence_data.append({'Gene': gene_info['gene'], 'Database': gene_info['database'], 'Count': gene_info['count'], 'Frequency': gene_info['frequency'], 'Percentage': f"{(gene_info['count']/len(integrated_data['samples']))*100:.1f}%" if len(integrated_data['samples']) > 0 else "0%", 'Genomes': ';'.join(gene_info.get('genomes', []))})
        if virulence_data:
            pd.DataFrame(virulence_data).to_csv(self.output_dir / "virulence_genes.csv", index=False)
        bacmet_data = []
        for db_name, genes in gene_centric.get('bacmet_databases', {}).items():
            for gene_info in genes:
                bacmet_data.append({'Gene': gene_info['gene'], 'Database': gene_info['database'], 'Count': gene_info['count'], 'Frequency': gene_info['frequency'], 'Percentage': f"{(gene_info['count']/len(integrated_data['samples']))*100:.1f}%" if len(integrated_data['samples']) > 0 else "0%", 'Genomes': ';'.join(gene_info.get('genomes', []))})
        if bacmet_data:
            pd.DataFrame(bacmet_data).to_csv(self.output_dir / "bacmet_genes.csv", index=False)
        plasmid_data = []
        for db_name, genes in gene_centric.get('plasmid_databases', {}).items():
            for gene_info in genes:
                plasmid_data.append({'Plasmid_Replicon': gene_info['gene'], 'Database': gene_info['database'], 'Count': gene_info['count'], 'Frequency': gene_info['frequency'], 'Percentage': f"{(gene_info['count']/len(integrated_data['samples']))*100:.1f}%" if len(integrated_data['samples']) > 0 else "0%", 'Genomes': ';'.join(gene_info.get('genomes', []))})
        if plasmid_data:
            pd.DataFrame(plasmid_data).to_csv(self.output_dir / "plasmid_replicons.csv", index=False)
        mutation_data = integrated_data.get('mutation_data', {})
        mutations = mutation_data.get('mutations', [])
        if mutations:
            mut_csv = []
            for m in mutations:
                mut_csv.append({'Gene': m['gene'], 'Mutation': m['mutation'], 'Class': m['class'], 'Subclass': m['subclass'], 'Count': m['count'], 'Genomes': ';'.join(m['genomes'])})
            pd.DataFrame(mut_csv).to_csv(self.output_dir / "mutations.csv", index=False)
        qc_data = integrated_data.get('qc_data', {})
        if qc_data:
            qc_rows = [{'Sample': sample, **metrics} for sample, metrics in qc_data.items()]
            pd.DataFrame(qc_rows).to_csv(self.output_dir / "fasta_qc.csv", index=False)
        pattern_data = []
        patterns = integrated_data['patterns']
        for mlst, count in patterns.get('mlst_distribution', Counter()).items():
            pattern_data.append({'Pattern_Type': 'MLST_Distribution', 'MLST': mlst, 'Count': count})
        for combo, samples in patterns.get('mlst_spa_combinations', {}).items():
            pattern_data.append({'Pattern_Type': 'MLST_spa_Combination', 'Combination': combo, 'Samples': ';'.join(samples), 'Count': len(samples)})
        for combo, samples in patterns.get('triple_combinations', {}).items():
            pattern_data.append({'Pattern_Type': 'Triple_Typing_(ST_spa_SCCmec)', 'Combination': combo, 'Samples': ';'.join(samples), 'Count': len(samples)})
        for combo in patterns.get('high_risk_combinations', []):
            pattern_data.append({'Pattern_Type': 'High_Risk_Combination', 'Sample': combo['sample'], 'MLST': combo['mlst'], 'spa_Type': combo['spa_type'], 'SCCmec_Type': combo['sccmec_type'], 'Critical_AMR_Genes': ';'.join(combo['critical_amr_genes']), 'Critical_Virulence_Genes': ';'.join(combo['critical_virulence_genes'])})
        if pattern_data:
            pd.DataFrame(pattern_data).to_csv(self.output_dir / "pattern_discovery.csv", index=False)
        print(f"    ✅ CSV reports generated: sample_overview.csv, amr_genes.csv, virulence_genes.csv, bacmet_genes.csv, plasmid_replicons.csv, fasta_qc.csv, pattern_discovery.csv" + (", mutations.csv" if mutations else ""))
    
    def run(self):
        print("=" * 80)
        print("🧬 STAPHSCOPE ULTIMATE S. AUREUS REPORTER v2.2.0")
        print("=" * 80)
        print(f"📁 Input directory: {self.input_dir}")
        html_files = self.find_html_files()
        if not any(html_files.values()):
            print("❌ No HTML report files found!")
            return False
        integrated_data = self.integrate_all_data(html_files)
        if not integrated_data:
            return False
        print("\n" + "=" * 80)
        print("📊 GENERATING ULTIMATE STAPHSCOPE REPORTS")
        print("=" * 80)
        self.generate_json_report(integrated_data)
        self.generate_csv_reports(integrated_data)
        self.html_generator.generate_main_report(integrated_data, self.output_dir)
        total_samples = len(integrated_data['samples'])
        patterns = integrated_data['patterns']
        high_risk = len(patterns.get('high_risk_combinations', []))
        gene_centric = integrated_data['gene_centric']
        total_amr_genes = sum(len(genes) for genes in gene_centric.get('amr_databases', {}).values())
        total_virulence_genes = sum(len(genes) for genes in gene_centric.get('virulence_databases', {}).values())
        total_bacmet_genes = sum(len(genes) for genes in gene_centric.get('bacmet_databases', {}).values())
        total_plasmids = sum(len(genes) for genes in gene_centric.get('plasmid_databases', {}).values())
        mrsa_count = 0
        for sample_data in integrated_data['samples'].values():
            if 'MRSA' in sample_data.get('typing', {}).get('MRSA_Status', ''):
                mrsa_count += 1
        mutation_count = len(integrated_data.get('mutation_data', {}).get('mutations', []))
        print("\n" + "=" * 80)
        print("✅ ULTIMATE ANALYSIS COMPLETE!")
        print("=" * 80)
        print(f"📁 Output directory: {self.output_dir}")
        print(f"📄 Files generated:")
        print(f"   • staphscope_ultimate_report.html (Interactive report with grouping)")
        print(f"   • staphscope_ultimate_report.json (Complete data)")
        print(f"   • sample_overview.csv (Sample data)")
        print(f"   • amr_genes.csv (Gene-centric AMR data)")
        print(f"   • virulence_genes.csv (Gene-centric virulence data)")
        print(f"   • bacmet_genes.csv (Biocide/heavy metal genes)")
        print(f"   • plasmid_replicons.csv (Plasmid replicon data)")
        if mutation_count > 0:
            print(f"   • mutations.csv (Point mutations)")
        if integrated_data.get('qc_data'):
            print(f"   • fasta_qc.csv (FASTA QC metrics)")
        print(f"   • pattern_discovery.csv (Pattern analysis including triple typing)")
        print(f"\n🔬 NEW FEATURE: Dynamic Grouping by Typing")
        print(f"   • In AMR, Virulence, BACMET, Plasmids, and Mutations tabs, use the grouping buttons to reorganise genome lists by MLST, spa, SCCmec, or combinations.")
        print(f"   • Instantly see which clones carry specific genes/mutations – powerful for outbreak and evolutionary studies.")
        print(f"\n🔬 KEY FEATURES FOR S. AUREUS:")
        print(f"   • Gene-centric tables with scrollable genome lists")
        print(f"   • MRSA focused analysis: {mrsa_count} MRSA samples identified")
        print(f"   • Complete spa typing: {len(patterns.get('spa_type_distribution', {}))} unique spa types")
        print(f"   • SCCmec analysis: {len(patterns.get('sccmec_distribution', {}))} SCCmec types")
        print(f"   • BACMET (biocide & heavy metal) resistance tracking")
        print(f"   • Triple typing combination (ST – spa – SCCmec)")
        print(f"   • FASTA QC metrics integrated")
        print(f"   • {mutation_count} unique point mutations across all genomes")
        print(f"   • Detailed biology explanations for each tab")
        print(f"\n📈 ANALYSIS SUMMARY:")
        print(f"   • {total_samples} total S. aureus samples analyzed")
        print(f"   • {mrsa_count} MRSA samples ({mrsa_count/total_samples*100:.1f}%)")
        print(f"   • {total_amr_genes} AMR genes across all databases")
        print(f"   • {total_virulence_genes} virulence genes")
        print(f"   • {total_bacmet_genes} BACMET (biocide/metal) genes")
        print(f"   • {total_plasmids} plasmid replicons")
        print(f"   • {mutation_count} unique point mutations")
        print(f"   • {high_risk} high-risk AMR+virulence combinations")
        print("\n🎯 Next steps:")
        print("   1. Open staphscope_ultimate_report.html in your browser")
        print("   2. Go to AMR, Virulence, BACMET, Plasmids, or Mutations tabs")
        print("   3. Click a grouping button (e.g., 'MLST' or 'Triple') to reorganise genomes by typing")
        print("   4. Use filter buttons to focus on key genes/mutations")
        print("   5. Examine triple typing combinations under Pattern Discovery tab")
        print("   6. Use print buttons on each section header to print specific sections")
        print("   7. Export data using the Export tab or individual CSV buttons")
        print("\n" + "=" * 80)
        return True


def main():
    parser = argparse.ArgumentParser(description='STAPHSCOPE Ultimate S. aureus Reporter - Gene‑Centric Analysis with Grouping by Typing', formatter_class=argparse.RawDescriptionHelpFormatter, epilog="""Examples:\n  python staphscope_ultimate_reporter.py -i /path/to/staphscope/reports\n\nAuthor: Brown Beckley <brownbeckley94@gmail.com>\nAffiliation: University of Ghana Medical School""")
    parser.add_argument('-i', '--input-dir', required=True, help='Directory containing StaphScope HTML report files')
    parser.add_argument('-o', '--output-dir', help='Custom output directory (optional)')
    args = parser.parse_args()
    input_dir = Path(args.input_dir)
    if not input_dir.exists():
        print(f"❌ Input directory not found: {input_dir}")
        sys.exit(1)
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
#!/usr/bin/env python3
"""
STAPHSCOPE ULTIMATE REPORTER - HYBRID GENE-CENTRIC & SAMPLE-CENTRIC
===================================================================
Version 1.0.0

Gene-centric for MLST/spa/SCCmec/Patterns.
Sample-centric interactive boxes for AMR, Virulence, BACMET, Plasmids, and Mutations.

Each sample-centric tab displays per-isolate boxes with:
- Sample name, total count, and typing badges (MLST, spa, SCCmec, MRSA/MSSA, agr).
- Horizontally scrollable tables with full details.
- Filtering by sample name and by database (where applicable).

Author: Brown Beckley <brownbeckley94@gmail.com>
Affiliation: University of Ghana Medical School
Date: 2026-07-12
"""

import os
import sys
import json
import re
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


# -----------------------------------------------------------------------------
# PARSER CLASS
# -----------------------------------------------------------------------------
class StaphHTMLParser:
    """Parser for StaphScope HTML and TSV reports."""

    def __init__(self):
        self.abricate_databases = [
            'card', 'resfinder', 'vfdb', 'argannot',
            'plasmidfinder', 'megares', 'ncbi', 'bacmet2'
        ]
        self.abricate_tsv_files = {
            'resfinder': 'staph_resfinder_abricate_summary.tsv',
            'vfdb': 'staph_vfdb_abricate_summary.tsv',
            'card': 'staph_card_abricate_summary.tsv',
            'argannot': 'staph_argannot_abricate_summary.tsv',
            'bacmet2': 'staph_bacmet2_abricate_summary.tsv',
            'plasmidfinder': 'staph_plasmidfinder_abricate_summary.tsv',
            'megares': 'staph_megares_abricate_summary.tsv',
            'ncbi': 'staph_ncbi_abricate_summary.tsv',
        }

    def normalize_sample_id(self, sample_id: str) -> str:
        sample = str(sample_id)
        extensions = ['.fna', '.fasta', '.fa', '.gb', '.gbk', '.gbff', '.txt', '.tsv', '.csv']
        for ext in extensions:
            if sample.endswith(ext):
                sample = sample[:-len(ext)]
        if '/' in sample or '\\' in sample:
            sample = Path(sample).name
        return sample.strip()

    # --------------------------------------------------------------------------
    # TSV LOADERS
    # --------------------------------------------------------------------------
    def load_typing_from_tsv(self, input_dir: Path) -> Dict[str, Dict]:
        tsv_path = input_dir / 'staphscope_comprehensive_report.tsv'
        if not tsv_path.exists():
            return {}
        df = pd.read_csv(tsv_path, sep='\t')
        typing = {}
        for _, row in df.iterrows():
            sample = row['sample']
            typing[sample] = {
                'MLST': row.get('mlst', 'ND'),
                'spa_Type': row.get('spa_type', 'ND'),
                'SCCmec_Type': row.get('sccmec_type', 'ND'),
                'MRSA_Status': row.get('mrsa_status', 'ND')
            }
        return typing

    def load_agr_from_tsv(self, input_dir: Path) -> Dict[str, Dict]:
        tsv_path = input_dir / 'agr_summary.tsv'
        if not tsv_path.exists():
            return {}
        df = pd.read_csv(tsv_path, sep='\t')
        agr_data = {}
        for _, row in df.iterrows():
            sample = row['Sample']
            sample = re.sub(r'\.(fna|fasta|fa)$', '', sample)
            agr_type = row.get('agr_Type', 'NA')
            agr_group = row.get('agr_Group', 'NA')
            match_score = row.get('Match_Score', '')
            canonical_agrD = row.get('Canonical_AgrD', '')
            multiple_agr = row.get('Multiple_Agr', '')
            status = row.get('Status', 'failed')
            agr_data[sample] = {
                'agr_Type': agr_type,
                'agr_Group': agr_group,
                'match_score': match_score,
                'canonical_agrD': canonical_agrD,
                'multiple_agr': multiple_agr,
                'status': status
            }
        return agr_data

    def load_amrfinder_from_tsv(self, input_dir: Path) -> Tuple[Dict[str, List[Dict]], Dict[str, int]]:
        tsv_path = input_dir / 'staph_amrfinder_summary.tsv'
        if not tsv_path.exists():
            return {}, {}
        df = pd.read_csv(tsv_path, sep='\t')
        amr_details = defaultdict(list)
        gene_counts = Counter()
        for _, row in df.iterrows():
            sample = row['Genome']
            gene_dict = row.to_dict()
            gene_dict = {k: (None if pd.isna(v) else v) for k, v in gene_dict.items()}
            if 'Gene_Symbol' in gene_dict:
                gene_dict['gene'] = gene_dict.pop('Gene_Symbol')
            amr_details[sample].append(gene_dict)
            gene_counts[gene_dict['gene']] += 1
        return dict(amr_details), dict(gene_counts)

    def load_abricate_from_tsv(self, input_dir: Path) -> Tuple[Dict[str, Dict[str, List[Dict]]], Dict[str, Dict[str, int]]]:
        abricate_details = defaultdict(lambda: defaultdict(list))
        abricate_gene_counts = defaultdict(lambda: defaultdict(int))

        for db, fname in self.abricate_tsv_files.items():
            path = input_dir / fname
            if not path.exists():
                continue
            df = pd.read_csv(path, sep='\t')
            sample_col = 'genome' if 'genome' in df.columns else 'file'
            for _, row in df.iterrows():
                sample = str(row[sample_col])
                sample = re.sub(r'\.(fasta|fna)$', '', sample)
                gene_dict = row.to_dict()
                gene_dict = {k: (None if pd.isna(v) else v) for k, v in gene_dict.items()}
                abricate_details[sample][db].append(gene_dict)
                gene = row['gene']
                abricate_gene_counts[db][gene] += 1

        return dict(abricate_details), dict(abricate_gene_counts)

    def load_mutations_from_tsv(self, input_dir: Path) -> Dict[str, List[Dict]]:
        tsv_path = input_dir / 'mutation_summary.tsv'
        if not tsv_path.exists():
            return {}
        df = pd.read_csv(tsv_path, sep='\t')
        mutations_by_sample = defaultdict(list)
        for _, row in df.iterrows():
            sample = row['genome']
            def clean_val(v):
                if pd.isna(v):
                    return ''
                return str(v)
            mut_dict = {
                'gene': clean_val(row.get('gene_symbol', '')),
                'mutation': clean_val(row.get('element_name', '')),
                'class': clean_val(row.get('class', '')),
                'subclass': clean_val(row.get('subclass', '')),
                'contig': clean_val(row.get('contig_id', '')),
                'start': clean_val(row.get('start', '')),
                'stop': clean_val(row.get('stop', '')),
                'strand': clean_val(row.get('strand', '')),
                'coverage': clean_val(row.get('coverage', '')),
                'identity': clean_val(row.get('identity', '')),
                'accession': clean_val(row.get('accession', ''))
            }
            mutations_by_sample[sample].append(mut_dict)
        return dict(mutations_by_sample)

    # --------------------------------------------------------------------------
    # HTML FALLBACK PARSERS
    # --------------------------------------------------------------------------
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

    def load_qc_from_html(self, file_path: Path) -> Dict[str, Dict]:
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
        print(f"  🧬 Parsing Comprehensive Report (HTML fallback): {file_path.name}")
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
            rows = typing_table.find_all('tr')
            if len(rows) < 2:
                return {}
            headers = []
            header_cells = rows[0].find_all(['th', 'td'])
            for cell in header_cells:
                headers.append(cell.get_text().strip())
            data = []
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
                results[sample] = {
                    'MLST': str(mlst).strip() if pd.notna(mlst) else 'ND',
                    'spa_Type': str(spa_type).strip() if pd.notna(spa_type) else 'ND',
                    'SCCmec_Type': str(sccmec_type).strip() if pd.notna(sccmec_type) else 'ND',
                    'MRSA_Status': str(mrsa_status).strip() if pd.notna(mrsa_status) else 'ND'
                }
            print(f"    ✓ Found {len(results)} samples")
            return results
        except Exception as e:
            print(f"    ❌ Error parsing comprehensive report: {e}")
            return {}

    def parse_amrfinder_report(self, file_path: Path) -> Tuple[Dict[str, List], Dict[str, Dict]]:
        print(f"  🧬 Parsing AMRfinder (HTML fallback): {file_path.name}")
        try:
            with open(file_path, 'r', encoding='utf-8') as f:
                html_content = f.read()
            soup = BeautifulSoup(html_content, 'html.parser')
            tables = soup.find_all('table')
            genes_by_genome = {}
            gene_frequencies = {}
            if tables:
                df = self.parse_html_table(html_content, 0)
                if not df.empty:
                    for _, row in df.iterrows():
                        sample = self.normalize_sample_id(row.get('Genome', row.get('genome', '')))
                        if sample:
                            genes_by_genome[sample] = {'all_genes': [], 'critical_genes': [], 'high_risk_genes': []}
            return genes_by_genome, gene_frequencies
        except Exception as e:
            print(f"    ❌ Error parsing AMRfinder: {e}")
            return {}, {}

    def parse_abricate_report(self, file_path: Path) -> Tuple[str, Dict[str, List], Dict[str, Dict]]:
        print(f"  🧬 Parsing ABRicate (HTML fallback): {file_path.name}")
        try:
            with open(file_path, 'r', encoding='utf-8') as f:
                html_content = f.read()
            db_name = 'unknown'
            filename = file_path.name.lower()
            for db in self.abricate_databases:
                if db in filename:
                    db_name = db
                    break
            df = self.parse_html_table(html_content, 0)
            genes_by_genome = {}
            if not df.empty:
                for _, row in df.iterrows():
                    sample = self.normalize_sample_id(row.get('genome', row.get('file', '')))
                    if sample:
                        genes_by_genome[sample] = []
            return db_name, genes_by_genome, {}
        except Exception as e:
            print(f"    ❌ Error parsing ABRicate: {e}")
            return 'unknown', {}, {}

    def parse_mutation_summary_html(self, file_path: Path) -> Dict[str, Any]:
        print(f"  🧬 Parsing mutation summary HTML: {file_path.name}")
        try:
            with open(file_path, 'r', encoding='utf-8') as f:
                html_content = f.read()
            soup = BeautifulSoup(html_content, 'html.parser')
            mutation_table = None
            for table in soup.find_all('table'):
                if table.find(string=re.compile(r'Gene', re.I)) and table.find(string=re.compile(r'Mutation', re.I)):
                    mutation_table = table
                    break
            if not mutation_table:
                print("    ⚠️ Could not find mutation table in HTML")
                return {}
            rows = mutation_table.find_all('tr')
            if len(rows) < 2:
                return {}
            headers = [th.get_text().strip() for th in rows[0].find_all(['th', 'td'])]
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
                    return {}
            mutations_list = []
            for row in rows[1:]:
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
                class_name = cells[col_idx['class']].get_text().strip() if 'class' in col_idx else ''
                subclass = cells[col_idx['subclass']].get_text().strip() if 'subclass' in col_idx else ''
                mutations_list.append({
                    'gene': gene,
                    'mutation': mutation,
                    'class': class_name,
                    'subclass': subclass,
                    'count': count,
                    'genomes': genomes
                })
            mutations_list.sort(key=lambda x: x['count'], reverse=True)
            return {'mutations': mutations_list}
        except Exception as e:
            print(f"    ❌ Error parsing mutation summary: {e}")
            return {}


# -----------------------------------------------------------------------------
# DATA ANALYZER
# -----------------------------------------------------------------------------
class StaphDataAnalyzer:
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
            'tsst', 'sea', 'seb', 'sec', 'sed', 'see',
            'seg', 'seh', 'sei', 'sej', 'sek', 'sel', 'sem', 'sen', 'seo', 'sep', 'seq', 'ser', 'seu',
            'eta', 'etb', 'hla', 'hlb', 'hlg', 'hld',
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
                if isinstance(data, dict):
                    count = data.get('count', 0)
                    frequency = data.get('frequency', str(count))
                    prevalence = data.get('prevalence', 'ND')
                    risk_level = data.get('risk_level', 'ND')
                    genomes = data.get('genomes', [])
                else:
                    count = int(data) if isinstance(data, (int, float)) else 0
                    frequency = str(count)
                    prevalence = 'ND'
                    risk_level = 'ND'
                    genomes = []
                gene_list.append({
                    'gene': gene,
                    'database': 'AMRfinder',
                    'frequency': frequency,
                    'count': count,
                    'prevalence': prevalence,
                    'risk_level': risk_level,
                    'genomes': genomes
                })
            gene_centric['amr_databases']['amrfinder'] = sorted(gene_list, key=lambda x: x['count'], reverse=True)

        if 'abricate' in integrated_data.get('gene_frequencies', {}):
            abricate_data = integrated_data['gene_frequencies']['abricate']
            for db_name, db_genes in abricate_data.items():
                gene_list = []
                for gene, data in db_genes.items():
                    if isinstance(data, dict):
                        count = data.get('count', 0)
                        frequency = data.get('frequency', str(count))
                        genomes = data.get('genomes', [])
                    else:
                        count = int(data) if isinstance(data, (int, float)) else 0
                        frequency = str(count)
                        genomes = []
                    gene_list.append({
                        'gene': gene,
                        'database': db_name.upper(),
                        'frequency': frequency,
                        'count': count,
                        'genomes': genomes
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


# -----------------------------------------------------------------------------
# HTML GENERATOR
# -----------------------------------------------------------------------------
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
            'export': '#9E9E9E',
            'agr': '#8B5CF6',
            'calltoaction': '#F472B6'
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
            agr = data.get('agr', {})
            sample_typing_js[sample] = {
                "MLST": typing.get('MLST', 'ND'),
                "spa": typing.get('spa_Type', 'ND'),
                "SCCmec": typing.get('SCCmec_Type', 'ND'),
                "MRSA": typing.get('MRSA_Status', 'ND'),
                "agr": agr.get('agr_Type', 'ND')
            }
        sample_typing_json = json.dumps(sample_typing_js)

        css = self._get_css()
        js = self._get_js(sample_typing_json)

        total_amr_genes = sum(len(genes) for genes in kwargs['gene_centric'].get('amr_databases', {}).values())
        total_virulence_genes = sum(len(genes) for genes in kwargs['gene_centric'].get('virulence_databases', {}).values())
        total_bacmet_genes = sum(len(genes) for genes in kwargs['gene_centric'].get('bacmet_databases', {}).values())
        total_plasmids = sum(len(genes) for genes in kwargs['gene_centric'].get('plasmid_databases', {}).values())

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
        <p>Gene‑Centric for Typing / Interactive Sample‑Centric for AMR, Virulence, BACMET, Plasmids & Mutations</p>
        <div class="metadata-bar">
            <div class="metadata-item"><i class="fas fa-calendar"></i><span>Generated: {kwargs['metadata'].get('analysis_date', 'Unknown')}</span></div>
            <div class="metadata-item"><i class="fas fa-database"></i><span>Samples: {len(kwargs['samples_data'])}</span></div>
            <div class="metadata-item"><i class="fas fa-user-md"></i><span>Tool: STAPHSCOPE Ultimate v1.0.0</span></div>
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
        <div class="dashboard-card card-agr" onclick="switchTab('agr')"><div class="card-number">{len(kwargs['integrated_data'].get('agr_data', {}))}</div><div class="card-label">Agr Typed</div><i class="fas fa-dna fa-2x" style="color: var(--agr-color); margin-top: 10px;"></i></div>
    </div>

    <div class="tab-navigation">
        <button class="tab-button summary active" onclick="switchTab('summary')"><i class="fas fa-chart-pie"></i> Summary</button>
        <button class="tab-button sample_overview" onclick="switchTab('sample_overview')"><i class="fas fa-list-alt"></i> Sample Overview</button>
        <button class="tab-button qc" onclick="switchTab('qc')"><i class="fas fa-chart-line"></i> FASTA QC</button>
        <button class="tab-button mlst" onclick="switchTab('mlst')"><i class="fas fa-code-branch"></i> MLST</button>
        <button class="tab-button spa" onclick="switchTab('spa')"><i class="fas fa-dna"></i> spa Typing</button>
        <button class="tab-button sccmec" onclick="switchTab('sccmec')"><i class="fas fa-shield-alt"></i> SCCmec</button>
        <button class="tab-button mrsa" onclick="switchTab('mrsa')"><i class="fas fa-skull-crossbones"></i> MRSA</button>
        <button class="tab-button agr" onclick="switchTab('agr')"><i class="fas fa-dna"></i> agr Typing</button>
        <button class="tab-button amr" onclick="switchTab('amr')"><i class="fas fa-biohazard"></i> AMR</button>
        <button class="tab-button virulence" onclick="switchTab('virulence')"><i class="fas fa-virus"></i> Virulence</button>
        <button class="tab-button bacmet" onclick="switchTab('bacmet')"><i class="fas fa-flask"></i> BACMET</button>
        <button class="tab-button plasmids" onclick="switchTab('plasmids')"><i class="fas fa-plug"></i> Plasmids</button>
        <button class="tab-button mutation" onclick="switchTab('mutation')"><i class="fas fa-dna"></i> Mutations</button>
        <button class="tab-button patterns" onclick="switchTab('patterns')"><i class="fas fa-project-diagram"></i> Patterns</button>
        <button class="tab-button aiguide" onclick="switchTab('aiguide')"><i class="fas fa-robot"></i> AI Guide</button>
        <button class="tab-button calltoaction" onclick="switchTab('calltoaction')"><i class="fas fa-globe"></i> Call to Action</button>
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
    <div id="agr-tab" class="tab-content">
        <h2 class="section-header agr-header"><i class="fas fa-dna"></i> Agr Typing Distribution<button class="print-section-btn" onclick="printSection('agr-tab')"><i class="fas fa-print"></i> Print</button></h2>
        {self._generate_agr_section(kwargs)}
    </div>
    <div id="amr-tab" class="tab-content">
        <h2 class="section-header amr-header"><i class="fas fa-biohazard"></i> Antimicrobial Resistance – Interactive Isolate Boxes<button class="print-section-btn" onclick="printSection('amr-tab')"><i class="fas fa-print"></i> Print</button></h2>
        {self._generate_sample_centric_boxes(kwargs, 'amr', 'AMR', ['amrfinder', 'resfinder', 'card', 'argannot', 'megares', 'ncbi'])}
    </div>
    <div id="virulence-tab" class="tab-content">
        <h2 class="section-header virulence-header"><i class="fas fa-virus"></i> Virulence Factors – Interactive Isolate Boxes<button class="print-section-btn" onclick="printSection('virulence-tab')"><i class="fas fa-print"></i> Print</button></h2>
        {self._generate_sample_centric_boxes(kwargs, 'virulence', 'Virulence', ['vfdb'])}
    </div>
    <div id="bacmet-tab" class="tab-content">
        <h2 class="section-header bacmet-header"><i class="fas fa-flask"></i> BACMET – Interactive Isolate Boxes<button class="print-section-btn" onclick="printSection('bacmet-tab')"><i class="fas fa-print"></i> Print</button></h2>
        {self._generate_sample_centric_boxes(kwargs, 'bacmet', 'BACMET', ['bacmet2'])}
    </div>
    <div id="plasmids-tab" class="tab-content">
        <h2 class="section-header plasmids-header"><i class="fas fa-plug"></i> Plasmid Replicons – Interactive Isolate Boxes<button class="print-section-btn" onclick="printSection('plasmids-tab')"><i class="fas fa-print"></i> Print</button></h2>
        {self._generate_sample_centric_boxes(kwargs, 'plasmids', 'Plasmids', ['plasmidfinder'])}
    </div>
    <div id="mutation-tab" class="tab-content">
        <h2 class="section-header mutation-header"><i class="fas fa-dna"></i> Point Mutations – Interactive Isolate Boxes<button class="print-section-btn" onclick="printSection('mutation-tab')"><i class="fas fa-print"></i> Print</button></h2>
        {self._generate_mutation_boxes(kwargs)}
    </div>
    <div id="patterns-tab" class="tab-content">
        <h2 class="section-header patterns-header"><i class="fas fa-project-diagram"></i> Cross-Genome Pattern Discovery<button class="print-section-btn" onclick="printSection('patterns-tab')"><i class="fas fa-print"></i> Print</button></h2>
        {self._generate_pattern_discovery_section(kwargs)}
    </div>
    <div id="aiguide-tab" class="tab-content">
        <h2 class="section-header aiguide-header"><i class="fas fa-robot"></i> AI Assistant Guide<button class="print-section-btn" onclick="printSection('aiguide-tab')"><i class="fas fa-print"></i> Print</button></h2>
        {self._generate_aiguide_section(kwargs)}
    </div>
    <div id="calltoaction-tab" class="tab-content">
        <h2 class="section-header calltoaction-header"><i class="fas fa-globe"></i> Call to Action – Fight AMR Together<button class="print-section-btn" onclick="printSection('calltoaction-tab')"><i class="fas fa-print"></i> Print</button></h2>
        {self._calltoaction_section()}
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
        <h3>STAPHSCOPE Ultimate S. aureus Reporter v1.0.0</h3>
        <p>University of Ghana Medical School | Brown Beckley &lt;brownbeckley94@gmail.com&gt;</p>
        <p>Generated on {kwargs['metadata'].get('analysis_date', 'Unknown')}</p>
        <p>⭐ Please give a big STAR on GitHub if you find this useful!</p>
    </div>
</div>
</body>
</html>
        """
        return html

    # --------------------------------------------------------------------------
    # CSS and JavaScript
    # --------------------------------------------------------------------------
    def _get_css(self) -> str:
        return """
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
            --agr-color: #8B5CF6;
            --calltoaction-color: #F472B6;
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
        .card-agr { border-left-color: var(--agr-color); }
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
        .tab-button.agr.active { background: var(--agr-color); }
        .tab-button.calltoaction.active { background: var(--calltoaction-color); }
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
        .agr-header { border-color: var(--agr-color); }
        .calltoaction-header { border-color: var(--calltoaction-color); }
        .data-table { width: 100%; border-collapse: collapse; margin: 20px 0; font-size: 0.95em; box-shadow: 0 2px 10px rgba(0,0,0,0.1); border-radius: 8px; overflow: hidden; table-layout: auto; }
        .data-table th { background: #2c3e50; color: white; padding: 15px; text-align: left; font-weight: 600; position: sticky; top: 0; white-space: nowrap; cursor: pointer; }
        .data-table th:hover { background: #1a252f; }
        .data-table td { padding: 12px; border-bottom: 1px solid #e0e0e0; vertical-align: top; white-space: nowrap; }
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

        /* Sample-centric box styles */
        .isolate-box {
            border: 1px solid #ddd;
            border-radius: 12px;
            margin-bottom: 30px;
            padding: 20px;
            background: #fafafa;
            box-shadow: 0 2px 8px rgba(0,0,0,0.06);
            border-bottom: 2px dashed #ccc;
        }
        .isolate-box .sample-header { display: flex; align-items: center; gap: 15px; flex-wrap: wrap; margin-bottom: 15px; border-bottom: 2px solid #e0e0e0; padding-bottom: 10px; }
        .isolate-box .sample-header h3 { font-size: 1.4em; margin: 0; }
        .isolate-box .sample-header .total-badge { background: #006400; color: white; padding: 4px 16px; border-radius: 20px; font-weight: bold; font-size: 0.9em; }
        .isolate-box .sample-header .typing-info { display: flex; flex-wrap: wrap; gap: 8px; margin-left: auto; }
        .isolate-box .sample-header .typing-badge { display: inline-block; padding: 3px 10px; border-radius: 12px; font-size: 0.8em; font-weight: 600; background: #e0e0e0; color: #333; border: 1px solid #ccc; }
        .isolate-box .sample-header .typing-badge.badge-mrsa { background: #8B0000; color: white; border-color: #8B0000; }
        .isolate-box .sample-header .typing-badge.badge-mssa { background: #4682B4; color: white; border-color: #4682B4; }
        .isolate-box .sample-header .typing-badge.agr-I { background: #16a34a; color: white; border-color: #16a34a; }
        .isolate-box .sample-header .typing-badge.agr-II { background: #2563eb; color: white; border-color: #2563eb; }
        .isolate-box .sample-header .typing-badge.agr-III { background: #f59e0b; color: white; border-color: #f59e0b; }
        .isolate-box .sample-header .typing-badge.agr-IV { background: #dc2626; color: white; border-color: #dc2626; }
        .isolate-box .sample-header .typing-badge.agr-NA { background: #6b7280; color: white; border-color: #6b7280; }
        .isolate-box .database-table-wrapper { margin: 15px 0; overflow-x: auto; }
        .isolate-box .database-table-wrapper table { width: 100%; border-collapse: collapse; font-size: 0.85em; min-width: 800px; }
        .isolate-box .database-table-wrapper table th { background: #2c3e50; color: white; padding: 8px 12px; text-align: left; white-space: nowrap; }
        .isolate-box .database-table-wrapper table td { padding: 8px 12px; border-bottom: 1px solid #e0e0e0; white-space: nowrap; }
        .isolate-box .database-table-wrapper table tr:hover { background: #f1f1f1; }
        .isolate-box .db-title { font-weight: bold; color: #006400; margin: 10px 0 5px 0; font-size: 1.1em; border-left: 4px solid #006400; padding-left: 10px; }
        .filter-controls { display: flex; flex-wrap: wrap; gap: 10px; align-items: center; background: #f8f9fa; padding: 15px; border-radius: 8px; margin-bottom: 20px; }
        .filter-controls select { padding: 10px; border-radius: 8px; border: 2px solid #ddd; background: white; min-width: 150px; }

        @media print { body * { visibility: hidden; } .tab-content.active, .tab-content.active * { visibility: visible; } .tab-content.active { position: absolute; left: 0; top: 0; width: 100%; padding: 20px; box-shadow: none; border-radius: 0; } .print-section-btn, .tab-navigation, .dashboard-grid, .search-box, .action-buttons, .grouping-controls, .filter-controls { display: none !important; } .data-table { page-break-inside: auto; } .data-table tr { page-break-inside: avoid; page-break-after: auto; } }
        @media (max-width: 768px) { .container { padding: 10px; } .main-header { padding: 20px; } .main-header h1 { font-size: 2em; } .tab-button { padding: 8px 12px; font-size: 0.8em; } .dashboard-grid { grid-template-columns: repeat(auto-fit, minmax(180px, 1fr)); } .data-table { font-size: 0.8em; } body { min-width: auto; overflow-x: auto; } }
        </style>
        """

    def _get_js(self, typing_json: str) -> str:
        return f"""
        <script>
        var sampleTyping = {typing_json};
        var originalGenomeLists = {{}};

        function switchTab(tabName) {{
            document.querySelectorAll('.tab-content').forEach(tab => tab.classList.remove('active'));
            document.querySelectorAll('.tab-button').forEach(button => button.classList.remove('active'));
            document.getElementById(tabName + '-tab').classList.add('active');
            event.currentTarget.classList.add('active');
            window.location.hash = tabName;
        }}

        // ---------- Gene‑centric grouping functions ----------
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
            if (groupBy === "agr") return info.agr;
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

        // ---------- Sample‑centric box filtering ----------
        function filterBoxes(tabId) {{
            var search = document.getElementById('search-' + tabId).value.toUpperCase();
            var dbFilter = document.getElementById('dbFilter-' + tabId).value;
            var boxes = document.querySelectorAll('#' + tabId + '-tab .isolate-box');
            boxes.forEach(function(box) {{
                var sample = box.getAttribute('data-sample') || '';
                var show = true;
                if (search && sample.toUpperCase().indexOf(search) === -1) show = false;
                var dbWrappers = box.querySelectorAll('.database-table-wrapper');
                if (dbFilter !== 'all') {{
                    dbWrappers.forEach(function(wrapper) {{
                        var dbName = wrapper.getAttribute('data-db') || '';
                        if (dbName === dbFilter) {{
                            wrapper.style.display = '';
                        }} else {{
                            wrapper.style.display = 'none';
                        }}
                    }});
                    var dbTitles = box.querySelectorAll('.db-title');
                    dbTitles.forEach(function(title) {{
                        var wrapper = title.nextElementSibling;
                        if (wrapper && wrapper.style.display === 'none') {{
                            title.style.display = 'none';
                        }} else {{
                            title.style.display = '';
                        }}
                    }});
                }} else {{
                    dbWrappers.forEach(function(wrapper) {{ wrapper.style.display = ''; }});
                    box.querySelectorAll('.db-title').forEach(function(title) {{ title.style.display = ''; }});
                }}
                box.style.display = show ? '' : 'none';
            }});
        }}

        function resetBoxFilters(tabId) {{
            document.getElementById('search-' + tabId).value = '';
            document.getElementById('dbFilter-' + tabId).value = 'all';
            filterBoxes(tabId);
        }}

        // ---------- General utilities ----------
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

    # --------------------------------------------------------------------------
    # SECTION GENERATORS
    # --------------------------------------------------------------------------
    def _generate_summary_section(self, kwargs: Dict) -> str:
        samples_data = kwargs['samples_data']
        patterns = kwargs['patterns']
        gene_centric = kwargs['gene_centric']
        integrated_data = kwargs['integrated_data']
        total_samples = len(samples_data)
        total_amr_genes = sum(len(genes) for genes in gene_centric.get('amr_databases', {}).values())
        total_virulence_genes = sum(len(genes) for genes in gene_centric.get('virulence_databases', {}).values())
        total_plasmids = sum(len(genes) for genes in gene_centric.get('plasmid_databases', {}).values())
        total_bacmet = sum(len(genes) for genes in gene_centric.get('bacmet_databases', {}).values())
        critical_findings = len(patterns.get('high_risk_combinations', []))
        mrsa_count = 0
        mssa_count = 0
        agr_count = len(integrated_data.get('agr_data', {}))
        agr_types = Counter()
        mutation_count = len(integrated_data.get('mutation_details', {}))
        for sample, data in samples_data.items():
            status = data.get('typing', {}).get('MRSA_Status', '')
            if 'MRSA' in status:
                mrsa_count += 1
            elif 'MSSA' in status:
                mssa_count += 1
            agr = data.get('agr', {})
            if agr.get('agr_Type') and agr['agr_Type'] != 'NA':
                agr_types[agr['agr_Type']] += 1

        agr_dist_str = ', '.join([f"{t}: {c}" for t, c in sorted(agr_types.items())]) if agr_types else 'None'

        html = f"""
        <div class="alert-box alert-info">
            <i class="fas fa-info-circle fa-2x"></i>
            <div>
                <h3>📊 Hybrid Analysis Overview – Gene‑Centric for Typing, Interactive Sample‑Centric for AMR/Virulence/BACMET/Plasmids/Mutations</h3>
                <p>This report analyzes <strong>{total_samples}</strong> <em>Staphylococcus aureus</em> genomes using a hybrid approach:</p>
                <ul>
                    <li><strong>Gene‑centric</strong> for typing tabs (MLST, spa, SCCmec, MRSA) – each gene/marker shown with all genomes that carry it.</li>
                    <li><strong>Interactive Sample‑Centric</strong> for AMR, Virulence, BACMET, Plasmids, and Mutations – each isolate gets its own box with detailed tables and typing badges.</li>
                </ul>
                <p><strong>NEW:</strong> Mutations tab now shows per‑isolate boxes, just like AMR/Virulence/BACMET/Plasmids, with all mutation details (gene, mutation, class, subclass, contig, start, stop, coverage, identity, accession).</p>
            </div>
        </div>
        <div class="alert-box alert-success">
            <i class="fas fa-magic fa-2x"></i>
            <div>
                <h3>📘 How to use the interactive isolate boxes</h3>
                <ol>
                    <li>Go to AMR, Virulence, BACMET, Plasmids, or Mutations.</li>
                    <li>Each isolate is displayed in its own box with a header showing sample name, total count, and typing badges.</li>
                    <li>Inside each box, you will find a table per database (or mutation table) with all details – scroll horizontally to see all columns.</li>
                    <li>Use the <strong>search bar</strong> to filter boxes by sample name.</li>
                    <li>For AMR/Virulence/BACMET/Plasmids, use the <strong>database dropdown</strong> to show only specific databases.</li>
                </ol>
            </div>
        </div>
        <div class="action-buttons">
            <button class="action-btn btn-primary" onclick="switchTab('amr')"><i class="fas fa-biohazard"></i> Explore AMR (Isolate Boxes)</button>
            <button class="action-btn btn-success" onclick="switchTab('virulence')"><i class="fas fa-virus"></i> Explore Virulence (Isolate Boxes)</button>
            <button class="action-btn btn-info" onclick="switchTab('bacmet')"><i class="fas fa-flask"></i> Explore BACMET (Isolate Boxes)</button>
            <button class="action-btn btn-warning" onclick="switchTab('mutation')"><i class="fas fa-dna"></i> Explore Mutations (Isolate Boxes)</button>
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
        <tr><td>Agr Typed Samples</td><td><strong>{agr_count}</strong></td><td>Distribution: {agr_dist_str}</td></tr>
        <tr><td>Samples with Mutations</td><td><strong>{mutation_count}</strong></td><td>Point mutations detected</td></tr>
        <tr><td>AMR Genes</td><td><strong>{total_amr_genes}</strong></td><td>Across AMRfinder, CARD, ResFinder, etc.</td></tr>
        <tr><td>Virulence Genes</td><td><strong>{total_virulence_genes}</strong></td><td>From VFDB (toxins, adhesins, immune evasion)</td></tr>
        <tr><td>BACMET (Biocide/Metal) Genes</td><td><strong>{total_bacmet}</strong></td><td>Biocide and heavy metal resistance (hospital environment)</td></tr>
        <tr><td>Plasmid Replicons</td><td><strong>{total_plasmids}</strong></td><td>Plasmid types – horizontal gene transfer potential</td></tr>
        <tr><td>High‑Risk AMR+Virulence Combos</td><td><span class="badge {'badge-critical' if critical_findings > 0 else 'badge-low'}">{critical_findings}</span></td><td>Samples with both critical AMR and virulence genes</td></tr>
        </tbody></table></div>
        <h3 style="margin-top: 30px;"><i class="fas fa-lightbulb"></i> ✨ Key Features of This Report</h3>
        <div style="display: grid; grid-template-columns: repeat(auto-fit, minmax(300px, 1fr)); gap: 20px; margin: 20px 0;">
            <div class="database-section"><h4><i class="fas fa-gene"></i> Gene‑Centric Typing</h4><p>MLST, spa, SCCmec, MRSA shown with all genomes – ideal for outbreak tracking.</p></div>
            <div class="database-section"><h4><i class="fas fa-box"></i> Interactive Isolate Boxes</h4><p>Each isolate has its own box with detailed gene/mutation tables per database, including all columns from the TSV files (Coverage, Identity, Contig, Start, Stop, etc.). Horizontally scrollable – no truncation. Displays typing badges including agr.</p></div>
            <div class="database-section"><h4><i class="fas fa-skull-crossbones"></i> MRSA Focus</h4><p>Dedicated MRSA analysis tab with ST‑spa, ST‑SCCmec, spa‑SCCmec, and triple combinations.</p></div>
            <div class="database-section"><h4><i class="fas fa-flask"></i> BACMET</h4><p>Biocide and heavy metal resistance genes – essential for hospital hygiene studies.</p></div>
            <div class="database-section"><h4><i class="fas fa-dna"></i> Mutations (Sample‑Centric)</h4><p>Per‑isolate point mutation tables with full details (gene, mutation, class, subclass, contig, start, stop, coverage, identity, accession).</p></div>
            <div class="database-section"><h4><i class="fas fa-print"></i> Section‑Specific Printing</h4><p>Print any tab individually with the “Print” button in each section header.</p></div>
            <div class="database-section"><h4><i class="fas fa-download"></i> Full Data Export</h4><p>All tables can be exported as CSV, and the complete dataset as JSON for downstream analysis or AI upload.</p></div>
            <div class="database-section"><h4><i class="fas fa-robot"></i> AI Assistant Guide</h4><p>Upload the JSON file to ChatGPT, Claude, or Gemini for interactive questions – with ethical guidelines.</p></div>
        </div>
        """
        return html

    def _generate_sample_overview_section(self, kwargs: Dict) -> str:
        samples_data = kwargs['samples_data']
        html = f"""
        <div class="alert-box alert-info">
            <i class="fas fa-info-circle fa-2x"></i>
            <div>
                <h3>🧬 Sample Overview – Population Snapshot</h3>
                <p>This table summarises the main typing methods for each <em>S. aureus</em> isolate:</p>
                <ul>
                    <li><strong>MLST</strong> – gold standard for global epidemiology.</li>
                    <li><strong>spa typing</strong> – high‑resolution outbreak tracking.</li>
                    <li><strong>SCCmec type</strong> – identifies MRSA lineage.</li>
                    <li><strong>agr type</strong> – accessory gene regulator (virulence regulation).</li>
                </ul>
                <p><strong>MRSA samples</strong> are highlighted in red.</p>
                <p><strong>Sortable columns:</strong> Click on any column header to sort.</p>
            </div>
        </div>
        <input type="text" class="search-box" id="search-samples" onkeyup="searchTable('samples-table', 'search-samples')" placeholder="🔍 Search samples by ID, ST, spa type, SCCmec, or agr...">
        <div class="action-buttons">
            <button class="action-btn btn-primary" onclick="exportTableToCSV('samples-table', 'sample_overview.csv')"><i class="fas fa-download"></i> Export to CSV</button>
            <button class="action-btn btn-success" onclick="document.getElementById('search-samples').value=''; searchTable('samples-table', 'search-samples')"><i class="fas fa-sync"></i> Clear Search</button>
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
                        <th data-sort="string">agr Type</th>
                        <th data-sort="number">Virulence Gene Count</th>
                    </tr>
                </thead>
                <tbody>
        """
        for sample, data in samples_data.items():
            typing = data.get('typing', {})
            agr = data.get('agr', {})
            mlst = typing.get('MLST', 'ND')
            spa_type = typing.get('spa_Type', 'ND')
            sccmec_type = typing.get('SCCmec_Type', 'ND')
            mrsa_status = typing.get('MRSA_Status', 'ND')
            agr_type = agr.get('agr_Type', 'ND')
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
                            <td>{agr_type}</td>
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
        qc_data = kwargs.get('integrated_data', {}).get('qc_data', {})
        if not qc_data:
            return """
            <div class="alert-box alert-warning">
                <i class="fas fa-exclamation-circle fa-2x"></i>
                <div><h3>No QC Data Available</h3><p>The FASTA_QC_summary.html file was not found or could not be parsed.</p></div>
            </div>
            """
        all_metrics = set()
        for metrics in qc_data.values():
            all_metrics.update(metrics.keys())
        metric_list = sorted(all_metrics)
        html = f"""
        <div class="alert-box alert-info">
            <i class="fas fa-chart-line fa-2x"></i>
            <div>
                <h3>📏 FASTA Quality Control</h3>
                <ul>
                    <li><strong>Number of contigs</strong> – lower is better; <em>Typical S. aureus: &lt;200 contigs (good), &lt;100 (excellent).</em></li>
                    <li><strong>N50</strong> – higher is better; <em>Typical S. aureus: &gt;50 kb (good), &gt;100 kb (excellent).</em></li>
                    <li><strong>GC%</strong> – S. aureus is typically 32‑33%.</li>
                    <li><strong>Total length</strong> – ~2.8 Mbp for a complete genome.</li>
                </ul>
            </div>
        </div>
        <input type="text" class="search-box" id="search-qc" onkeyup="searchTable('qc-table', 'search-qc')" placeholder="🔍 Search sample...">
        <div class="action-buttons">
            <button class="action-btn btn-primary" onclick="exportTableToCSV('qc-table', 'fasta_qc.csv')"><i class="fas fa-download"></i> Export QC Data</button>
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
        patterns = kwargs['patterns']
        mlst_dist = patterns.get('mlst_distribution', Counter())
        mlst_spa_combos = patterns.get('mlst_spa_combinations', {})
        mlst_sccmec_combos = patterns.get('mlst_sccmec_combinations', {})
        html = f"""
        <div class="alert-box alert-info">
            <i class="fas fa-code-branch fa-2x"></i>
            <div>
                <h3>🔬 MLST (Multi‑Locus Sequence Typing)</h3>
                <p>MLST indexes internal fragments of seven housekeeping genes. Each unique combination defines a Sequence Type (ST).</p>
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
            if mlst == 'ND':
                continue
            pct = (count / total) * 100 if total > 0 else 0
            associated_spa = []
            for combo in mlst_spa_combos:
                if f"{mlst} - " in combo:
                    spa = combo.split(" - ")[1]
                    if spa not in associated_spa:
                        associated_spa.append(spa)
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
        patterns = kwargs['patterns']
        spa_dist = patterns.get('spa_type_distribution', Counter())
        mlst_spa_combos = patterns.get('mlst_spa_combinations', {})
        spa_sccmec_combos = patterns.get('spa_sccmec_combinations', {})
        html = f"""
        <div class="alert-box alert-info">
            <i class="fas fa-dna fa-2x"></i>
            <div>
                <h3>🧬 spa Typing – High‑Resolution Outbreak Tracking</h3>
                <p>The <em>spa</em> gene encodes protein A. Repeat region polymorphisms define spa types.</p>
                <p><strong>{len(spa_dist)} unique spa types</strong> identified.</p>
            </div>
        </div>
        <h3>📊 spa Type Distribution</h3>
        <div class="scrollable-table"><table id="spa-table" class="data-table"><thead><tr><th data-sort="string">spa Type</th><th data-sort="number">Count</th><th data-sort="number">Percentage</th><th data-sort="string">Common STs</th></tr></thead><tbody>
        """
        total = sum(spa_dist.values())
        for spa, count in spa_dist.most_common():
            if spa == 'ND':
                continue
            pct = (count / total) * 100 if total > 0 else 0
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
            rev_combo = f"{parts[1]} - {parts[0]}" if len(parts) == 2 else combo
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
        patterns = kwargs['patterns']
        sccmec_dist = patterns.get('sccmec_distribution', Counter())
        mlst_sccmec_combos = patterns.get('mlst_sccmec_combinations', {})
        spa_sccmec_combos = patterns.get('spa_sccmec_combinations', {})
        html = f"""
        <div class="alert-box alert-info">
            <i class="fas fa-shield-alt fa-2x"></i>
            <div>
                <h3>🧬 SCCmec Typing – The MRSA Cassette</h3>
                <p>The staphylococcal cassette chromosome <em>mec</em> (SCCmec) carries <em>mecA</em> or <em>mecC</em>.</p>
                <p><strong>{len(sccmec_dist)} unique SCCmec types</strong> identified.</p>
            </div>
        </div>
        <h3>📊 SCCmec Type Distribution</h3>
        <div class="scrollable-table"><table id="sccmec-table" class="data-table"><thead><tr><th data-sort="string">SCCmec Type</th><th data-sort="number">Count</th><th data-sort="number">Percentage</th><th data-sort="string">Common STs</th><th data-sort="string">Common spa Types</th></tr></thead><tbody>
        """
        total = sum(sccmec_dist.values())
        for scc, count in sccmec_dist.most_common():
            if scc in ['ND', 'Not Assigned']:
                continue
            pct = (count / total) * 100 if total > 0 else 0
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
            rev_combo = f"{parts[1]} - {parts[0]}" if len(parts) == 2 else combo
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
            rev_combo = f"{parts[1]} - {parts[0]}" if len(parts) == 2 else combo
            sample_tags = ''.join(f'<span class="genome-tag">{s}</span>' for s in samples)
            html += f"<tr><td><strong>{rev_combo}</strong></td><td>{len(samples)}</td><td><div class='genome-list'>{sample_tags}</div></td></tr>"
        html += "</tbody></table></div>"
        return html

    def _generate_mrsa_section(self, kwargs: Dict) -> str:
        patterns = kwargs['patterns']
        mrsa_status_dist = patterns.get('mrsa_status_distribution', Counter())
        samples_data = kwargs['samples_data']
        has_typing = any(
            d.get('typing', {}).get('MLST', 'ND') != 'ND' or
            d.get('typing', {}).get('spa_Type', 'ND') != 'ND' or
            d.get('typing', {}).get('SCCmec_Type', 'ND') != 'ND'
            for d in samples_data.values()
        )
        if not has_typing:
            return """
            <div class="alert-box alert-warning">
                <i class="fas fa-exclamation-triangle fa-2x"></i>
                <div><h3>⚠️ No Typing Data Available</h3><p>The comprehensive typing report was not found or could not be parsed.</p></div>
            </div>
            """
        mrsa_samples = [s for s, d in samples_data.items() if 'MRSA' in d.get('typing', {}).get('MRSA_Status', '')]
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

        html = f"""
        <div class="alert-box alert-danger">
            <i class="fas fa-skull-crossbones fa-2x"></i>
            <div>
                <h3>⚠️ MRSA (Methicillin‑Resistant S. aureus) – A Clinical Priority</h3>
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
                if status == 'ND':
                    continue
                pct = (count / total) * 100
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

        if not mrsa_samples:
            html += """
            <div class="alert-box alert-info">
                <i class="fas fa-info-circle fa-2x"></i>
                <div><strong>No MRSA samples detected</strong> – combination tables are empty.</div>
            </div>
            """
            return html

        html += f"""
        <h3>🔗 MRSA ST–spa Combinations ({len(mrsa_mlst_spa)} combinations)</h3>
        <input type="text" class="search-box" id="search-mrsa-mlst-spa" onkeyup="searchTable('mrsa-mlst-spa-table','search-mrsa-mlst-spa')" placeholder="🔍 Search ST‑spa...">
        <div class="master-scrollable-container"><table id="mrsa-mlst-spa-table" class="data-table"><thead><tr><th data-sort="string">ST‑spa Combination</th><th data-sort="number">Count</th><th data-sort="string">Samples</th></tr></thead><tbody>
        """
        for combo, samples in sorted(mrsa_mlst_spa.items(), key=lambda x: len(x[1]), reverse=True):
            sample_tags = ''.join(f'<span class="genome-tag">{s}</span>' for s in samples)
            html += f"<tr><td><strong>{combo}</strong></td><td>{len(samples)}</td><td><div class='genome-list'>{sample_tags}</div></td></tr>"
        html += "</tbody></table></div>"

        html += f"""
        <h3>🔗 MRSA ST–SCCmec Combinations ({len(mrsa_mlst_sccmec)} combinations)</h3>
        <input type="text" class="search-box" id="search-mrsa-mlst-sccmec" onkeyup="searchTable('mrsa-mlst-sccmec-table','search-mrsa-mlst-sccmec')" placeholder="🔍 Search ST‑SCCmec...">
        <div class="master-scrollable-container"><table id="mrsa-mlst-sccmec-table" class="data-table"><thead><tr><th data-sort="string">ST‑SCCmec Combination</th><th data-sort="number">Count</th><th data-sort="string">Samples</th></tr></thead><tbody>
        """
        for combo, samples in sorted(mrsa_mlst_sccmec.items(), key=lambda x: len(x[1]), reverse=True):
            sample_tags = ''.join(f'<span class="genome-tag">{s}</span>' for s in samples)
            html += f"<tr><td><strong>{combo}</strong></td><td>{len(samples)}</td><td><div class='genome-list'>{sample_tags}</div></td></tr>"
        html += "</tbody></table></div>"

        html += f"""
        <h3>🔗 MRSA spa–SCCmec Combinations ({len(mrsa_spa_sccmec)} combinations)</h3>
        <input type="text" class="search-box" id="search-mrsa-spa-sccmec" onkeyup="searchTable('mrsa-spa-sccmec-table','search-mrsa-spa-sccmec')" placeholder="🔍 Search spa‑SCCmec...">
        <div class="master-scrollable-container"><table id="mrsa-spa-sccmec-table" class="data-table"><thead><tr><th data-sort="string">spa‑SCCmec Combination</th><th data-sort="number">Count</th><th data-sort="string">Samples</th></tr></thead><tbody>
        """
        for combo, samples in sorted(mrsa_spa_sccmec.items(), key=lambda x: len(x[1]), reverse=True):
            sample_tags = ''.join(f'<span class="genome-tag">{s}</span>' for s in samples)
            html += f"<tr><td><strong>{combo}</strong></td><td>{len(samples)}</td><td><div class='genome-list'>{sample_tags}</div></td></tr>"
        html += "</tbody></table></div>"
        return html

    # --------------------------------------------------------------------------
    # AGR SECTION
    # --------------------------------------------------------------------------
    def _generate_agr_section(self, kwargs: Dict) -> str:
        integrated_data = kwargs.get('integrated_data', {})
        agr_data = integrated_data.get('agr_data', {})
        samples_data = kwargs.get('samples_data', {})

        # Count agr types
        agr_counts = Counter()
        for sample, data in samples_data.items():
            agr = data.get('agr', {})
            if agr.get('agr_Type') and agr['agr_Type'] != 'NA':
                agr_counts[agr['agr_Type']] += 1
            else:
                agr_counts['NA'] += 1

        # Build agr type -> samples mapping
        agr_to_samples = defaultdict(list)
        for sample, data in samples_data.items():
            agr = data.get('agr', {}).get('agr_Type', 'NA')
            if agr != 'NA':
                agr_to_samples[agr].append(sample)
        for sample, data in samples_data.items():
            agr = data.get('agr', {}).get('agr_Type', 'NA')
            if agr == 'NA':
                agr_to_samples['NA'].append(sample)

        # Build combination maps
        combo_functions = {
            'agr-MLST': lambda s, d: f"{d.get('agr', {}).get('agr_Type', 'NA')} - {d.get('typing', {}).get('MLST', 'ND')}",
            'agr-spa': lambda s, d: f"{d.get('agr', {}).get('agr_Type', 'NA')} - {d.get('typing', {}).get('spa_Type', 'ND')}",
            'agr-SCCmec': lambda s, d: f"{d.get('agr', {}).get('agr_Type', 'NA')} - {d.get('typing', {}).get('SCCmec_Type', 'ND')}",
            'agr-MLST-spa': lambda s, d: f"{d.get('agr', {}).get('agr_Type', 'NA')} - {d.get('typing', {}).get('MLST', 'ND')} - {d.get('typing', {}).get('spa_Type', 'ND')}",
            'agr-MLST-SCCmec': lambda s, d: f"{d.get('agr', {}).get('agr_Type', 'NA')} - {d.get('typing', {}).get('MLST', 'ND')} - {d.get('typing', {}).get('SCCmec_Type', 'ND')}",
            'agr-spa-SCCmec': lambda s, d: f"{d.get('agr', {}).get('agr_Type', 'NA')} - {d.get('typing', {}).get('spa_Type', 'ND')} - {d.get('typing', {}).get('SCCmec_Type', 'ND')}",
            'agr-MLST-spa-SCCmec': lambda s, d: f"{d.get('agr', {}).get('agr_Type', 'NA')} - {d.get('typing', {}).get('MLST', 'ND')} - {d.get('typing', {}).get('spa_Type', 'ND')} - {d.get('typing', {}).get('SCCmec_Type', 'ND')}",
        }

        combo_data = {}
        for name, func in combo_functions.items():
            combos = defaultdict(list)
            for sample, data in samples_data.items():
                agr = data.get('agr', {}).get('agr_Type', 'NA')
                if agr == 'NA':
                    continue
                key = func(sample, data)
                if 'ND' not in key and 'Not Assigned' not in key:
                    combos[key].append(sample)
            combo_data[name] = dict(combos)

        html = f"""
        <div class="alert-box alert-info">
            <i class="fas fa-dna fa-2x"></i>
            <div>
                <h3>🧬 Agr Typing – Explore agr with All Typing Schemes</h3>
                <p>The <strong>accessory gene regulator (agr)</strong> system is a quorum‑sensing circuit that controls virulence gene expression in <em>S. aureus</em>. Four agr types (I‑IV) are recognised, with distinct epidemiological and clinical associations.</p>
                <p><strong>{len(agr_data)} samples</strong> were successfully typed. Use the tables below to explore how agr type correlates with MLST, spa, SCCmec, and combinations.</p>
            </div>
        </div>
        """

        # 1. Agr Type Distribution
        html += """
        <h3>📊 Agr Type Distribution</h3>
        <div class="scrollable-table">
            <table class="data-table">
                <thead><tr><th>agr Type</th><th>Count</th><th>Percentage</th></tr></thead>
                <tbody>
        """
        total = sum(agr_counts.values())
        for typ in ['I', 'II', 'III', 'IV', 'NA']:
            count = agr_counts.get(typ, 0)
            pct = (count / total * 100) if total > 0 else 0
            color_class = f"agr-{typ}" if typ != 'NA' else 'agr-NA'
            html += f"""
                <tr>
                    <td><span class="typing-badge {color_class}">{typ}</span></td>
                    <td>{count}</td>
                    <td>{pct:.1f}%</td>
                </tr>
            """
        html += """
                </tbody>
            </table>
        </div>
        """

        # 2. Samples by Agr Type 
        html += """
        <h3>📋 Samples by agr Type</h3>
        <input type="text" class="search-box" id="search-agr-samples" onkeyup="searchTable('agr-samples-table','search-agr-samples')" placeholder="🔍 Search agr type...">
        <input type="text" class="search-box" id="highlight-agr-samples" onkeyup="highlightGenome('agr-samples-table','highlight-agr-samples')" placeholder="🔍 Highlight genomes containing specific text...">
        <div class="master-scrollable-container">
            <table id="agr-samples-table" class="data-table">
                <thead>
                    <tr>
                        <th data-sort="string">agr Type</th>
                        <th data-sort="number">Count</th>
                        <th data-sort="string">Genomes</th>
                    </tr>
                </thead>
                <tbody>
        """
        for typ in ['I', 'II', 'III', 'IV', 'NA']:
            samples = agr_to_samples.get(typ, [])
            if samples:
                color_class = f"agr-{typ}" if typ != 'NA' else 'agr-NA'
                sample_tags = ''.join(f'<span class="genome-tag">{s}</span>' for s in sorted(samples))
                html += f"""
                    <tr>
                        <td><span class="typing-badge {color_class}">{typ}</span></td>
                        <td>{len(samples)}</td>
                        <td><div class="genome-list">{sample_tags}</div></td>
                    </tr>
                """
        html += """
                </tbody>
            </table>
        </div>
        """

        # 3. All combination tables
        table_configs = [
            ('agr-MLST', 'agr-mlst', 'agr-MLST', 'agr_MLST'),
            ('agr-spa', 'agr-spa', 'agr-spa', 'agr_spa'),
            ('agr-SCCmec', 'agr-sccmec', 'agr-SCCmec', 'agr_sccmec'),
            ('agr-MLST-spa', 'agr-mlst-spa', 'agr-MLST-spa', 'agr_mlst_spa'),
            ('agr-MLST-SCCmec', 'agr-mlst-sccmec', 'agr-MLST-SCCmec', 'agr_mlst_sccmec'),
            ('agr-spa-SCCmec', 'agr-spa-sccmec', 'agr-spa-SCCmec', 'agr_spa_sccmec'),
            ('agr-MLST-spa-SCCmec', 'agr-four-way', 'agr-MLST-spa-SCCmec', 'agr_four_way'),
        ]

        for key, table_id, display_name, search_prefix in table_configs:
            combos = combo_data.get(key, {})
            if not combos:
                continue

            html += f"""
        <h3>🔗 {display_name} Combinations</h3>
        <input type="text" class="search-box" id="search-{table_id}" onkeyup="searchTable('{table_id}','search-{table_id}')" placeholder="🔍 Search {display_name} combination...">
        <input type="text" class="search-box" id="highlight-{table_id}" onkeyup="highlightGenome('{table_id}','highlight-{table_id}')" placeholder="🔍 Highlight genomes containing specific text...">
        <div class="master-scrollable-container">
            <table id="{table_id}" class="data-table">
                <thead>
                    <tr>
                        <th data-sort="string">{display_name}</th>
                        <th data-sort="number">Count</th>
                        <th data-sort="string">Samples</th>
                    </tr>
                </thead>
                <tbody>
            """
            for combo, samples in sorted(combos.items(), key=lambda x: len(x[1]), reverse=True):
                sample_tags = ''.join(f'<span class="genome-tag">{s}</span>' for s in samples)
                html += f"""
                    <tr>
                        <td><strong>{combo}</strong></td>
                        <td>{len(samples)}</td>
                        <td><div class="genome-list">{sample_tags}</div></td>
                    </tr>
                """
            html += """
                </tbody>
            </table>
        </div>
        """

        return html

    # --------------------------------------------------------------------------
    # SAMPLE‑CENTRIC BOX GENERATOR 
    # --------------------------------------------------------------------------
    def _generate_sample_centric_boxes(self, kwargs: Dict, tab_id: str, title: str, db_list: List[str]) -> str:
        samples_data = kwargs.get('samples_data', {})
        amr_details = kwargs.get('integrated_data', {}).get('amrfinder_details', {})
        abricate_details = kwargs.get('integrated_data', {}).get('abricate_details', {})

        relevant_samples = []
        for sample in samples_data:
            has_data = False
            if 'amrfinder' in db_list and sample in amr_details and amr_details[sample]:
                has_data = True
            for db in db_list:
                if db != 'amrfinder' and sample in abricate_details and db in abricate_details[sample] and abricate_details[sample][db]:
                    has_data = True
            if has_data:
                relevant_samples.append(sample)
        relevant_samples.sort()

        if not relevant_samples:
            return f"""
            <div class="alert-box alert-warning">
                <i class="fas fa-exclamation-circle fa-2x"></i>
                <div><h3>No {title} Data Available</h3><p>No samples with {title} genes were found.</p></div>
            </div>
            """

        html = f"""
        <div class="alert-box alert-info">
            <i class="fas fa-info-circle fa-2x"></i>
            <div>
                <h3>🧬 {title} – Interactive Isolate Boxes</h3>
                <p>Each box represents one isolate. Inside, you will find separate tables for each database with full gene details – horizontally scrollable.</p>
                <p>Use the filters below to search by sample name or to show only a specific database.</p>
            </div>
        </div>
        <div class="filter-controls">
            <input type="text" class="search-box" id="search-{tab_id}" onkeyup="filterBoxes('{tab_id}')" placeholder="🔍 Search sample...">
            <select id="dbFilter-{tab_id}" onchange="filterBoxes('{tab_id}')">
                <option value="all">All Databases</option>
        """
        for db in db_list:
            display = db.upper()
            if db == 'amrfinder':
                display = 'AMRfinder'
            html += f'<option value="{db}">{display}</option>'
        html += f"""
            </select>
            <button class="action-btn btn-success" onclick="resetBoxFilters('{tab_id}')"><i class="fas fa-sync"></i> Clear Filters</button>
        </div>
        <div id="box-container-{tab_id}">
        """

        for sample in relevant_samples:
            typing = samples_data.get(sample, {}).get('typing', {})
            agr = samples_data.get(sample, {}).get('agr', {})
            mlst = typing.get('MLST', 'ND')
            spa = typing.get('spa_Type', 'ND')
            sccmec = typing.get('SCCmec_Type', 'ND')
            mrsa = typing.get('MRSA_Status', 'ND')
            mrsa_class = 'badge-mrsa' if 'MRSA' in mrsa else 'badge-mssa' if 'MSSA' in mrsa else ''
            agr_type = agr.get('agr_Type', 'ND')
            agr_class = f"agr-{agr_type}" if agr_type != 'ND' else 'agr-NA'

            html += f'<div class="isolate-box" data-sample="{sample}">'
            total_genes = 0
            if 'amrfinder' in db_list and sample in amr_details:
                total_genes += len(amr_details[sample])
            for db in db_list:
                if db != 'amrfinder' and sample in abricate_details and db in abricate_details[sample]:
                    total_genes += len(abricate_details[sample][db])
            html += f"""
                <div class="sample-header">
                    <h3><i class="fas fa-microbe"></i> {sample}</h3>
                    <span class="total-badge">Total Genes: {total_genes}</span>
                    <div class="typing-info">
                        <span class="typing-badge">ST: {mlst}</span>
                        <span class="typing-badge">spa: {spa}</span>
                        <span class="typing-badge">SCCmec: {sccmec}</span>
                        <span class="typing-badge {mrsa_class}">{mrsa}</span>
                        <span class="typing-badge {agr_class}">agr: {agr_type}</span>
                    </div>
                </div>
            """

            if 'amrfinder' in db_list and sample in amr_details and amr_details[sample]:
                html += self._make_database_table('AMRfinder', amr_details[sample], 'amrfinder')

            for db in db_list:
                if db == 'amrfinder':
                    continue
                if sample in abricate_details and db in abricate_details[sample] and abricate_details[sample][db]:
                    html += self._make_database_table(db.upper(), abricate_details[sample][db], db)

            html += '</div>'

        html += """
        </div>
        """
        return html

    def _make_database_table(self, db_name: str, genes: List[Dict], db_key: str) -> str:
        if not genes:
            return ''
        keys = list(genes[0].keys()) if genes else []
        priority = ['gene', 'product', 'coverage_percent', 'identity_percent', 'accession', 'contig', 'start', 'stop', 'class', 'subclass', 'scope', 'resistance']
        ordered = []
        for p in priority:
            if p in keys:
                ordered.append(p)
        for k in keys:
            if k not in ordered:
                ordered.append(k)

        html = f"""
                <div class="database-table-wrapper" data-db="{db_key}">
                    <div class="db-title">{db_name}</div>
                    <table>
                        <thead>
                            <tr>
        """
        for col in ordered:
            display = col.replace('_', ' ').title()
            html += f'<th>{display}</th>'
        html += '</tr></thead><tbody>'
        for gene_dict in genes:
            html += '<tr>'
            for col in ordered:
                val = gene_dict.get(col, '')
                if val is None:
                    val = ''
                html += f'<td>{val}</td>'
            html += '</tr>'
        html += """
                        </tbody>
                    </table>
                </div>
        """
        return html

    # --------------------------------------------------------------------------
    # MUTATION BOXES (SAMPLE‑CENTRIC)
    # --------------------------------------------------------------------------
    def _generate_mutation_boxes(self, kwargs: Dict) -> str:
        samples_data = kwargs.get('samples_data', {})
        integrated_data = kwargs.get('integrated_data', {})
        mutation_details = integrated_data.get('mutation_details', {})
        relevant_samples = [s for s in samples_data if s in mutation_details and mutation_details[s]]
        relevant_samples.sort()

        if not relevant_samples:
            return """
            <div class="alert-box alert-warning">
                <i class="fas fa-exclamation-circle fa-2x"></i>
                <div><h3>No Mutation Data Available</h3><p>No mutations found for any sample.</p></div>
            </div>
            """

        html = f"""
        <div class="alert-box alert-info">
            <i class="fas fa-dna fa-2x"></i>
            <div>
                <h3>🧬 Point Mutations – Sample‑Centric View</h3>
                <p>Each box represents one isolate. Inside, you will find a table of all detected point mutations with full details – horizontally scrollable.</p>
                <p>Use the search bar below to filter boxes by sample name.</p>
            </div>
        </div>
        <div class="filter-controls">
            <input type="text" class="search-box" id="search-mutation" onkeyup="filterBoxes('mutation')" placeholder="🔍 Search sample...">
            <button class="action-btn btn-success" onclick="resetBoxFilters('mutation')"><i class="fas fa-sync"></i> Clear Filters</button>
        </div>
        <div id="box-container-mutation">
        """

        for sample in relevant_samples:
            typing = samples_data.get(sample, {}).get('typing', {})
            agr = samples_data.get(sample, {}).get('agr', {})
            mlst = typing.get('MLST', 'ND')
            spa = typing.get('spa_Type', 'ND')
            sccmec = typing.get('SCCmec_Type', 'ND')
            mrsa = typing.get('MRSA_Status', 'ND')
            mrsa_class = 'badge-mrsa' if 'MRSA' in mrsa else 'badge-mssa' if 'MSSA' in mrsa else ''
            agr_type = agr.get('agr_Type', 'ND')
            agr_class = f"agr-{agr_type}" if agr_type != 'ND' else 'agr-NA'

            mutations = mutation_details.get(sample, [])
            total_mutations = len(mutations)

            html += f'<div class="isolate-box" data-sample="{sample}">'
            html += f"""
                <div class="sample-header">
                    <h3><i class="fas fa-microbe"></i> {sample}</h3>
                    <span class="total-badge">Total Mutations: {total_mutations}</span>
                    <div class="typing-info">
                        <span class="typing-badge">ST: {mlst}</span>
                        <span class="typing-badge">spa: {spa}</span>
                        <span class="typing-badge">SCCmec: {sccmec}</span>
                        <span class="typing-badge {mrsa_class}">{mrsa}</span>
                        <span class="typing-badge {agr_class}">agr: {agr_type}</span>
                    </div>
                </div>
            """

            if mutations:
                html += self._make_mutation_table(mutations)

            html += '</div>'

        html += """
        </div>
        """
        return html

    def _make_mutation_table(self, mutations: List[Dict]) -> str:
        if not mutations:
            return ''
        cols = ['gene', 'mutation', 'class', 'subclass', 'contig', 'start', 'stop', 'strand', 'coverage', 'identity', 'accession']
        html = """
        <div class="database-table-wrapper" data-db="mutations">
            <div class="db-title">Mutations</div>
            <table>
                <thead>
                    <tr>
        """
        for col in cols:
            display = col.replace('_', ' ').title()
            html += f'<th>{display}</th>'
        html += '</tr></thead><tbody>'
        for mut in mutations:
            html += '<tr>'
            for col in cols:
                val = mut.get(col, '')
                if val is None or (isinstance(val, float) and val != val):
                    val = ''
                html += f'<td>{val}</td>'
            html += '</tr>'
        html += """
                </tbody>
            </table>
        </div>
        """
        return html

    # --------------------------------------------------------------------------
    # OTHER SECTIONS 
    # --------------------------------------------------------------------------
    def _generate_pattern_discovery_section(self, kwargs: Dict) -> str:
        patterns = kwargs['patterns']
        triple_combos = patterns.get('triple_combinations', {})
        html = f"""<div class="alert-box alert-info"><i class="fas fa-project-diagram fa-2x"></i><div><h3>🔍 Cross‑Genome Pattern Discovery</h3><p>This tab reveals associations between typing results, gene co‑occurrence, and high‑risk combinations.</p></div></div>
        <h3>🔗 Triple Typing Combination (ST – spa – SCCmec)</h3>
        <input type="text" class="search-box" id="search-triple" onkeyup="searchTable('triple-table','search-triple')" placeholder="🔍 Search combination...">
        <div class="master-scrollable-container"><table id="triple-table" class="data-table"><thead><tr><th data-sort="string">Combination</th><th data-sort="number">Count</th><th data-sort="string">Samples</th></tr></thead><tbody>"""
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
            <div class="alert-box alert-info">
                <i class="fas fa-robot fa-2x"></i>
                <div>
                    <h3>🤖 AI Assistant Guide – Unleash the Power of AI for Genomic Epidemiology</h3>
                    <p>This guide shows you how to leverage large language models (LLMs) like <strong>ChatGPT, Claude, or Gemini</strong> to interact with your <em>S. aureus</em> dataset – turning static reports into dynamic conversations.</p>
                </div>
            </div>

            <div style="margin: 20px 0;">
                <!-- Why AI for genomics -->
                <div class="database-section">
                    <h4><i class="fas fa-brain"></i> Why Use AI for Genomic Data Analysis?</h4>
                    <p>Modern AI models excel at:</p>
                    <ul>
                        <li><strong>Pattern recognition</strong> – spotting epidemiological trends, clone associations, and co‑occurrence networks.</li>
                        <li><strong>Natural language queries</strong> – ask in plain English, get instant answers without writing code.</li>
                        <li><strong>Hypothesis generation</strong> – uncover unexpected correlations that merit experimental follow‑up.</li>
                        <li><strong>Literature synthesis</strong> – connect your findings with published resistance mechanisms and clinical guidelines.</li>
                    </ul>
                    <p><span style="background: #fff3cd; padding: 2px 8px; border-radius: 4px;"><i class="fas fa-lightbulb"></i> <strong>Scientific note:</strong> AI is <em>pattern‑finding</em>, not <em>causal‑inferring</em>. Use it to suggest, then verify with wet‑lab or clinical correlation.</span></p>
                </div>

                <!-- How to use AI with this report -->
                <div class="database-section">
                    <h4><i class="fas fa-upload"></i> How to Feed This Report to AI</h4>
                    <p>You have <strong>three powerful options</strong>:</p>
                    <ol>
                        <li><strong>Upload the JSON file</strong> – The file <code>staphscope_ultimate_report.json</code> contains all structured data. Upload it to ChatGPT (Advanced Data Analysis), Claude, or Gemini. <em>Best for precise quantitative queries.</em></li>
                        <li><strong>Upload the HTML report</strong> – Modern AI tools can read HTML and extract tables. You can upload the <code>.html</code> file directly – the AI will parse the tables and text. <em>Great for visual context.</em></li>
                        <li><strong>Copy‑paste specific tables</strong> – If you only need a quick insight, copy a table (e.g., AMR gene list) and paste it into the chat. <em>Instant, no file upload needed.</em></li>
                    </ol>
                    <p style="margin-top: 10px; background: #e8f5e9; padding: 10px; border-radius: 5px;">
                        <i class="fas fa-info-circle"></i> <strong>Pro tip:</strong> For best results, tell the AI: <em>"You are a bioinformatician analysing <em>S. aureus</em> genomes. The attached data contains typing, AMR, virulence, BACMET, and mutation information. Answer my questions with references to the data."</em>
                    </p>
                </div>

                <!-- Example questions (categorized) -->
                <div class="database-section">
                    <h4><i class="fas fa-chart-line"></i> Scientifically Relevant Questions to Ask</h4>
                    <p>Use these as starting points – adapt to your specific research question:</p>
                    <div style="display: grid; grid-template-columns: 1fr 1fr; gap: 10px;">
                        <div style="background: #f8f9fa; padding: 10px; border-radius: 8px;">
                            <strong>🧬 Epidemiology &amp; Clonality</strong>
                            <ul style="margin-top: 5px; font-size: 0.9em;">
                                <li>What are the most common MLST sequence types in this dataset?</li>
                                <li>Which spa types are dominant in MRSA vs MSSA?</li>
                                <li>Are there any ST‑spa‑SCCmec combinations with >2 isolates?</li>
                                <li>What is the agr type distribution, and does it correlate with MRSA status?</li>
                                <li>Which clones carry the most resistance genes?</li>
                            </ul>
                        </div>
                        <div style="background: #f8f9fa; padding: 10px; border-radius: 8px;">
                            <strong>💊 Antimicrobial Resistance</strong>
                            <ul style="margin-top: 5px; font-size: 0.9em;">
                                <li>How many samples carry mecA? What are their STs and spa types?</li>
                                <li>Are there any vanA/vanB positive samples? What is their SCCmec type?</li>
                                <li>Which AMR genes co‑occur most frequently?</li>
                                <li>What is the distribution of tetracycline (tet) resistance genes?</li>
                                <li>Do any samples have combined β‑lactam + macrolide resistance?</li>
                            </ul>
                        </div>
                        <div style="background: #f8f9fa; padding: 10px; border-radius: 8px;">
                            <strong>🦠 Virulence &amp; Toxins</strong>
                            <ul style="margin-top: 5px; font-size: 0.9em;">
                                <li>Which samples carry PVL (lukF/S-PV)? Are they associated with specific STs or agr types?</li>
                                <li>List all samples with TSST‑1 (tsst).</li>
                                <li>Which enterotoxin genes are most prevalent?</li>
                                <li>Is there a correlation between biofilm (ica) genes and MRSA?</li>
                                <li>Do any isolates carry both immune evasion and cytotoxin genes?</li>
                            </ul>
                        </div>
                        <div style="background: #f8f9fa; padding: 10px; border-radius: 8px;">
                            <strong>🧪 Mutations &amp; Biocides</strong>
                            <ul style="margin-top: 5px; font-size: 0.9em;">
                                <li>What are the most frequent point mutations in gyrA or parC?</li>
                                <li>Are there any linezolid‑related mutations (23S rRNA)?</li>
                                <li>Which samples carry qac genes (disinfectant resistance)?</li>
                                <li>Is mer (mercury) resistance linked to specific STs?</li>
                                <li>Which isolates have both efflux pumps and biocide resistance?</li>
                            </ul>
                        </div>
                    </div>
                    <p style="margin-top: 10px;"><i class="fas fa-arrow-right"></i> <strong>Beyond tables:</strong> Ask the AI to <em>“write a summary of the resistance profile for ST5”</em> or <em>“compare virulence gene carriage between MRSA and MSSA”</em> – the report contains all data.</p>
                </div>

                <!-- Scientific verification & ethics -->
                <div class="database-section">
                    <h4><i class="fas fa-balance-scale"></i> Scientific Rigour &amp; Ethical AI Use</h4>
                    <ul>
                        <li><strong>AI is your co‑pilot, not the pilot.</strong> Always interpret AI‑generated insights in the context of your local epidemiology, clinical guidelines, and laboratory validation.</li>
                        <li><strong>Verify, verify, verify.</strong> Cross‑check critical calls (e.g., MRSA status, vancomycin resistance) with primary literature, genome browsers, or secondary tools.</li>
                        <li><strong>No patient‑identifiable data.</strong> Only upload aggregated, de‑identified genomic data. This report contains no clinical metadata.</li>
                        <li><strong>Transparency in publications.</strong> If you use AI for exploratory analysis, mention it in the methods (e.g., “AI‑assisted pattern discovery was performed using a large language model, followed by manual curation”).</li>
                        <li><strong>AI hallucination is real.</strong> If the AI confidently tells you that <em>mecA</em> is found in <em>E. coli</em> or that TSST‑1 causes indigestion – <strong>don’t believe it</strong>. Treat every AI statement as a hypothesis, not a fact.</li>
                    </ul>
                </div>

                <!-- Humorous warnings -->
                <div class="database-section" style="background: #fff3cd; border-left: 6px solid #ffc107;">
                    <h4><i class="fas fa-smile-wink"></i> A (Mostly Serious) AI Survival Guide</h4>
                    <ul>
                        <li><strong>If the AI says “I don’t know”</strong> – trust it. It’s being honest.</li>
                        <li><strong>If the AI says “It is widely known”</strong> – ask for a reference. It may have made it up.</li>
                        <li><strong>If the AI offers a completely novel evolutionary theory</strong> – check if your coffee is spiked.</li>
                        <li><strong>Remember:</strong> AI won’t take your job – but a microbiologist who knows how to use AI might! So learn it, use it, and always keep a healthy dose of scepticism. 😉</li>
                    </ul>
                    <p style="margin-top: 10px;"><i class="fas fa-microbe"></i> <strong>Final thought:</strong> The best AI‑human partnership is one where the AI does the heavy pattern‑lifting, and you do the heavy thinking. Happy (and careful) exploring!</p>
                </div>
            </div>
            """

    def _generate_citation_section(self, kwargs: Dict) -> str:
            return """
            <div class="alert-box alert-info"><i class="fas fa-quote-right fa-2x"></i><div><h3>📚 How to Cite StaphScope and Its Dependencies</h3><p>If you use StaphScope in your research, please cite the main tool and the relevant third‑party tools and databases.</p></div></div>
            <div class="accordion"><div class="accordion-item"><div class="accordion-header"><span>📄 From the ESKAPE AMR Platform</span><i class="fas fa-chevron-down"></i></div><div class="accordion-content" style="display: none;"><ul class="citation-list"><li><strong>StaphScope</strong> – Beckley B, Amarh V. StaphScope: a species‑optimized computational pipeline for rapid and accessible <em>Staphylococcus aureus</em> genotyping and surveillance. <em>BMC Genomics</em>. 2026;27:261. doi:10.1186/s12864-026-12609-x<button class="copy-btn" data-citation='Beckley B, Amarh V. StaphScope: a species‑optimized computational pipeline for rapid and accessible Staphylococcus aureus genotyping and surveillance. BMC Genomics. 2026;27:261. doi:10.1186/s12864-026-12609-x'>📋 Copy citation</button></li></ul></div></div><div class="accordion-item"><div class="accordion-header"><span>🔧 Key Databases & Methods</span><i class="fas fa-chevron-down"></i></div><div class="accordion-content" style="display: none;"><ul class="citation-list">
            <li><strong>MLST</strong> – Seemann T. MLST: Scan contig files against PubMLST typing schemes. GitHub. 2018. https://github.com/tseemann/mlst<button class="copy-btn" data-citation='Seemann T. MLST: Scan contig files against PubMLST typing schemes. GitHub. 2018. https://github.com/tseemann/mlst'>📋 Copy citation</button></li>
            <li><strong>PubMLST / BIGSdb</strong> – Jolley KA, Bray JE, Maiden MCJ. Open‑access bacterial population genomics: BIGSdb software, the PubMLST.org website and their applications. <em>Wellcome Open Res</em>. 2018;3:124. doi:10.12688/wellcomeopenres.14826.1<button class="copy-btn" data-citation='Jolley KA, Bray JE, Maiden MCJ. Open-access bacterial population genomics: BIGSdb software, the PubMLST.org website and their applications. Wellcome Open Res. 2018;3:124. doi:10.12688/wellcomeopenres.14826.1'>📋 Copy citation</button></li>
            <li><strong>spa typing</strong> – Harmsen D, et al. Typing of methicillin‑resistant <em>Staphylococcus aureus</em> in a university hospital setting by using novel software for spa repeat determination and database management. <em>J Clin Microbiol</em>. 2003;41(12):5442-8. doi:10.1128/JCM.41.12.5442-5448.2003<button class="copy-btn" data-citation='Harmsen D, et al. Typing of methicillin-resistant Staphylococcus aureus in a university hospital setting by using novel software for spa repeat determination and database management. J Clin Microbiol. 2003;41(12):5442-8. doi:10.1128/JCM.41.12.5442-5448.2003'>📋 Copy citation</button></li>
            <li><strong>SCCmecFinder</strong> – Kaya H, et al. SCCmecFinder, a Web‑Based Tool for Typing of Staphylococcal Cassette Chromosome mec in <em>Staphylococcus aureus</em> Using Whole‑Genome Sequence Data. <em>mSphere</em>. 2018;3(1):e00612-17. doi:10.1128/mSphere.00612-17<button class="copy-btn" data-citation='Kaya H, et al. SCCmecFinder, a Web‑Based Tool for Typing of Staphylococcal Cassette Chromosome mec in Staphylococcus aureus Using Whole‑Genome Sequence Data. mSphere. 2018;3(1):e00612-17. doi:10.1128/mSphere.00612-17'>📋 Copy citation</button></li>
            <li><strong>agrVATE</strong> – Raghuram V, Alexander AM, Loo HQ, Petit RA 3rd, Goldberg JB, Read TD. Species-Wide Phylogenomics of the Staphylococcus aureus Agr Operon Revealed Convergent Evolution of Frameshift Mutations. <em>Microbiol Spectr</em>. 2022;10(1):e0133421. doi:10.1128/spectrum.01334-21<button class="copy-btn" data-citation='Raghuram V, Alexander AM, Loo HQ, Petit RA 3rd, Goldberg JB, Read TD. Species-Wide Phylogenomics of the Staphylococcus aureus Agr Operon Revealed Convergent Evolution of Frameshift Mutations. Microbiol Spectr. 2022;10(1):e0133421. doi:10.1128/spectrum.01334-21'>📋 Copy citation</button></li>
            <li><strong>AMRFinderPlus</strong> – Feldgarden M, et al. AMRFinderPlus and the Reference Gene Catalog facilitate examination of the genomic links among antimicrobial resistance, stress response, and virulence. <em>Sci Rep</em>. 2021;11(1):12728. doi:10.1038/s41598-021-91456-0<button class="copy-btn" data-citation='Feldgarden M, et al. AMRFinderPlus and the Reference Gene Catalog facilitate examination of the genomic links among antimicrobial resistance, stress response, and virulence. Sci Rep. 2021;11(1):12728. doi:10.1038/s41598-021-91456-0'>📋 Copy citation</button></li>
            <li><strong>ABRicate</strong> – Seemann T. ABRicate: mass screening of contigs for antibiotic resistance genes. GitHub. 2024. https://github.com/tseemann/abricate<button class="copy-btn" data-citation='Seemann T. ABRicate: mass screening of contigs for antibiotic resistance genes. GitHub. 2024. https://github.com/tseemann/abricate'>📋 Copy citation</button></li>
            <li><strong>CARD</strong> – McArthur AG, et al. The comprehensive antibiotic resistance database. <em>Antimicrob Agents Chemother</em>. 2013;57(7):3348-57. doi:10.1128/AAC.00419-13<button class="copy-btn" data-citation='McArthur AG, et al. The comprehensive antibiotic resistance database. Antimicrob Agents Chemother. 2013;57(7):3348-57. doi:10.1128/AAC.00419-13'>📋 Copy citation</button></li>
            <li><strong>ResFinder</strong> – Florensa AF, et al. ResFinder – an open online resource for identification of antimicrobial resistance genes in next‑generation sequencing data and prediction of phenotypes from genotypes. <em>Microb Genom</em>. 2022;8(1):000748. doi:10.1099/mgen.0.000748<button class="copy-btn" data-citation='Florensa AF, et al. ResFinder – an open online resource for identification of antimicrobial resistance genes in next-generation sequencing data and prediction of phenotypes from genotypes. Microb Genom. 2022;8(1):000748. doi:10.1099/mgen.0.000748'>📋 Copy citation</button></li>
            <li><strong>VFDB</strong> – Chen L, et al. VFDB 2012 update: toward the genetic diversity and molecular evolution of bacterial virulence factors. <em>Nucleic Acids Res</em>. 2012;40(Database issue):D641-5. doi:10.1093/nar/gkr989<button class="copy-btn" data-citation='Chen L, et al. VFDB 2012 update: toward the genetic diversity and molecular evolution of bacterial virulence factors. Nucleic Acids Res. 2012;40(Database issue):D641-5. doi:10.1093/nar/gkr989'>📋 Copy citation</button></li>
            <li><strong>PlasmidFinder</strong> – Carattoli A, et al. <em>In silico</em> detection and typing of plasmids using PlasmidFinder and plasmid multilocus sequence typing. <em>Antimicrob Agents Chemother</em>. 2014;58(7):3895-903. doi:10.1128/AAC.02412-14<button class="copy-btn" data-citation='Carattoli A, et al. In silico detection and typing of plasmids using PlasmidFinder and plasmid multilocus sequence typing. Antimicrob Agents Chemother. 2014;58(7):3895-903. doi:10.1128/AAC.02412-14'>📋 Copy citation</button></li>
            <li><strong>BacMet</strong> – Pal C, et al. BacMet: antibacterial biocide and metal resistance genes database. <em>Nucleic Acids Res</em>. 2014;42(Database issue):D737-43. doi:10.1093/nar/gkt1252<button class="copy-btn" data-citation='Pal C, et al. BacMet: antibacterial biocide and metal resistance genes database. Nucleic Acids Res. 2014;42(Database issue):D737-43. doi:10.1093/nar/gkt1252'>📋 Copy citation</button></li>
            <li><strong>MEGARes</strong> – Doster E, et al. MEGARes 2.0: a database for classification of antimicrobial drug, biocide and metal resistance determinants in metagenomic sequence data. <em>Nucleic Acids Res</em>. 2020;48(D1):D561-D569. doi:10.1093/nar/gkz1010<button class="copy-btn" data-citation='Doster E, et al. MEGARes 2.0: a database for classification of antimicrobial drug, biocide and metal resistance determinants in metagenomic sequence data. Nucleic Acids Res. 2020;48(D1):D561-D569. doi:10.1093/nar/gkz1010'>📋 Copy citation</button></li>
            <li><strong>ARG-ANNOT</strong> – Gupta SK, et al. ARG-ANNOT, a new bioinformatic tool to discover antibiotic resistance genes in bacterial genomes. <em>Antimicrob Agents Chemother</em>. 2014;58(1):212-20. doi:10.1128/AAC.01310-13<button class="copy-btn" data-citation='Gupta SK, et al. ARG-ANNOT, a new bioinformatic tool to discover antibiotic resistance genes in bacterial genomes. Antimicrob Agents Chemother. 2014;58(1):212-20. doi:10.1128/AAC.01310-13'>📋 Copy citation</button></li>
            <li><strong>Biopython</strong> – Cock PJ, et al. Biopython: freely available Python tools for computational molecular biology and bioinformatics. <em>Bioinformatics</em>. 2009;25(11):1422-3. doi:10.1093/bioinformatics/btp163<button class="copy-btn" data-citation='Cock PJ, et al. Biopython: freely available Python tools for computational molecular biology and bioinformatics. Bioinformatics. 2009;25(11):1422-3. doi:10.1093/bioinformatics/btp163'>📋 Copy citation</button></li>
        </ul></div></div></div>
        <div class="alert-box alert-success" style="margin-top: 20px;"><i class="fas fa-hand-peace"></i><div><strong>Suggested acknowledgement:</strong><br>
        "Genomic analysis was performed using StaphScope [Beckley &amp; Amarh, 2026], which integrates MLST [Seemann, 2018] using the PubMLST database [Jolley et al., 2018], ABRicate [Seemann, 2018], AMRFinderPlus [Feldgarden et al., 2021], SCCmecFinder [Kaya et al., 2018], and agrVATE [Raghuram et al., 2022] for comprehensive <em>S. aureus</em> characterization. Antimicrobial resistance genes were identified using the CARD [McArthur et al., 2013], ResFinder [Florensa et al., 2022], MEGARes [Doster et al., 2020], and ARG-ANNOT [Gupta et al., 2014] databases. For biocide and heavy metal resistance genes, BacMet [Pal et al., 2014] was used. Virulence and plasmid screening were performed with ABRicate using the VFDB [Chen et al., 2012] and PlasmidFinder [Carattoli et al., 2014] databases. Mutation detection was performed using AMRFinderPlus. FASTA QC was performed using Biopython [Cock et al., 2009]."
        </div></div>
            """

    def _generate_funding_section(self, kwargs: Dict) -> str:
        return """
        <div class="alert-box alert-info"><i class="fas fa-coffee fa-2x"></i><div><h3>☕ Funding & Support – Keeping the Lights On (with code and caffeine)</h3><p>StaphScope is an <strong>independent, unfunded project</strong> born out of passion for genomic surveillance and AMR research at the University of Ghana Medical School.</p><p>No grants, no sponsors, no institutional backing – just a laptop, a lot of coffee, and a burning desire to help researchers fight antimicrobial resistance.</p></div></div>
        <div class="alert-box alert-warning"><i class="fas fa-heart fa-2x"></i><div><h3>💡 How You Can Help (Without Opening Your Wallet)</h3><ul><li><strong>⭐ Star us on GitHub</strong> – It takes two seconds and makes us feel like rockstars.</li><li><strong>🐛 Report bugs</strong> – If something breaks, let us know. We’ll fix it with joy.</li><li><strong>💡 Suggest features</strong> – Have an idea? We’re all ears (and we actually implement them, as you’ve seen!).</li><li><strong>🧬 Share your data</strong> – If you’ve used StaphScope and want to collaborate, we’d love to hear your story.</li><li><strong>📢 Spread the word</strong> – Tell your colleagues, tweet about it, or mention it in your next Zoom call.</li><li><strong>👋 Say hello!</strong> – Seriously, just drop an email to <strong>brownbeckley94@gmail.com</strong>. It makes our day (Literally My Day).</li></ul><p><i class="fas fa-microbe"></i> <strong>Fun fact:</strong> This project runs on 100% volunteer tears, 0% grant money. But we’re not bitter – we’re just caffeinated.</p></div></div>
        <div class="alert-box alert-success"><i class="fas fa-hand-holding-heart"></i><div><h3>🤝 Contribute to the ESKAPE AMR Platform</h3><p>We also maintain pipelines for other ESKAPE pathogens (AcinetoScope, Kleboscope, Pseudoscope, etc.). If you’re a developer, bioinformatician, or just someone who loves clean code and bacteria, we welcome: pull requests, issues, documentation improvements, ideas for new databases. Visit our GitHub: <a href="https://github.com/bbeckley-hub" target="_blank">https://github.com/bbeckley-hub</a> – star, fork, and let’s fight AMR together!</p><p><strong>Brown Beckley</strong> – <i class="fas fa-envelope"></i> brownbeckley94@gmail.com</p><p><i class="fas fa-laugh-beam"></i> <strong>P.S.</strong> If you ever meet Brown in person, buy him a coffee or I'll buy you a bug myself. He’ll probably talk your ear off about SCCmec types, but it’s worth it.</p></div></div>
        """

    def _calltoaction_section(self):
        return """
        <div class="alert-box alert-info"><i class="fas fa-globe fa-2x"></i><div>
        <h3>The Global Burden of AMR and Our Call to Action</h3>
        <p>Antimicrobial resistance (AMR) is one of the top global public health threats, with an estimated <strong>1.27 million direct deaths annually</strong>. <em>Staphylococcus aureus</em> is a major contributor, causing skin and soft tissue infections, bloodstream infections, pneumonia, and endocarditis. Tracking AMR and virulence determinants is essential to inform treatment guidelines and infection control.</p>
        <p>We developed <strong>StaphScope</strong> to empower researchers and clinicians – especially in low‑resource settings – to analyse their own sequencing data without extensive bioinformatics expertise.</p>
        </div></div>
        
        <div style="background:#e8f5e9; padding:20px; border-radius:12px; margin:20px 0;">
        <h3><i class="fas fa-bacterium"></i> ESCAPE AMR – Our Ongoing Project (ESKAPE Pathogens)</h3>
        <p><strong>StaphScope</strong> is one of the first modules of a larger initiative called <strong>ESCAPE AMR</strong> (formerly ESKAPE). We target the notorious <strong>ESKAPE pathogens</strong>:</p>
        <ul>
            <li><strong>E</strong>nterococcus faecium</li>
            <li><strong>S</strong>taphylococcus aureus</li>
            <li><strong>K</strong>lebsiella pneumoniae</li>
            <li><strong>A</strong>cinetobacter baumannii</li>
            <li><strong>P</strong>seudomonas aeruginosa</li>
            <li><strong>E</strong>nterobacter species</li>
        </ul>
        <p>These bacteria “escape” the effects of antibiotics – hence the name. But we believe the name is also a global call to action:</p>
        <div style="background:#fff3e0; padding:15px; border-radius:8px; margin:15px 0;">
            <p><strong>🔹 E</strong>veryone must join forces – researchers, clinicians, policymakers, and citizens.<br>
            <strong>🔹 S</strong>mart surveillance is our first line of defence. No more guessing – we need genomic data.<br>
            <strong>🔹 K</strong>nowledge must be shared openly. No paywalls, no closed silos.<br>
            <strong>🔹 A</strong>frica bears a heavy AMR burden, but African solutions are already emerging.<br>
            <strong>🔹 P</strong>revention is cheaper than cure. Let's stop resistant infections before they spread.<br>
            <strong>🔹 E</strong>very day we delay, more lives are at stake. The time to act is now, not tomorrow.</p>
        </div>
        <p><i class="fas fa-laugh-squint"></i> <strong>“We didn’t choose the name ESKAPE because it sounds cool (though it does). We chose it because it reminds us every single day: we must ESCAPE the AMR crisis – together, urgently, and with the best science we have.”</strong><br>
        — Brown Beckley, lead developer (who secretly hopes this pun makes you smile, not roll your eyes 😉)</p>
        </div>
        
        <div style="text-align:center; margin:40px 0;">
            <i class="fas fa-star" style="font-size:3em; color:#ffc107;"></i>
            <h3>🤝 We Invite You to Contribute!</h3>
            <p><strong>If you find this tool useful, please:</strong></p>
            <div class="action-buttons" style="justify-content:center;">
                <a href="https://github.com/bbeckley-hub/staphscope-typing-tool" target="_blank" class="action-btn btn-primary" style="text-decoration:none;"><i class="fab fa-github"></i> ⭐ Star us on GitHub</a>
                <a href="mailto:brownbeckley94@gmail.com" class="action-btn btn-success" style="text-decoration:none;"><i class="fas fa-envelope"></i> Contact the team</a>
                <a href="https://github.com/bbeckley-hub/staphscope-typing-tool/issues" target="_blank" class="action-btn btn-warning" style="text-decoration:none;"><i class="fas fa-bug"></i> Report issues</a>
            </div>
            <p style="margin-top:20px;"><i class="fas fa-chalkboard-user"></i> <strong>We welcome collaborations</strong> to adapt this tool for other pathogens and to improve AMR surveillance in Africa and beyond.</p>
            <p><i class="fas fa-hand-holding-heart"></i> If you are a funder or organisation interested in supporting the <strong>ESCAPE AMR</strong> project, please reach out. Together we can build a free, open‑source ecosystem for genomic surveillance of all ESKAPE pathogens.</p>
        </div>
        
        <div class="info-text" style="background:#f8f9fa;">
        <i class="fas fa-quote-left"></i> “AMR is a silent pandemic, but we have the tools to fight it – if we share them, if we teach each other, and if we act with urgency. Let's escape the era of untreatable infections.”<br>
        <strong>— The ESCAPE AMR Team, University of Ghana Medical School</strong>
        </div>
        """

    def _generate_export_section(self, kwargs: Dict) -> str:
        return """
        <div class="alert-box alert-info"><i class="fas fa-download fa-2x"></i><div><h3>📥 Export Data</h3><p>Download tables as CSV using the buttons in each tab, or get the complete JSON data.</p></div></div>
        <div style="display: grid; grid-template-columns: repeat(auto-fit, minmax(280px, 1fr)); gap: 20px; margin: 30px 0;">
            <div class="dashboard-card card-export" onclick="exportTableToCSV('samples-table', 'sample_overview.csv')"><div style="font-size:2em;color:var(--export-color);"><i class="fas fa-table"></i></div><div class="card-label">Sample Overview CSV</div></div>
            <div class="dashboard-card card-export" onclick="exportTableToCSV('amr-table', 'amr_genes.csv')"><div style="font-size:2em;color:var(--export-color);"><i class="fas fa-biohazard"></i></div><div class="card-label">AMR Genes CSV</div></div>
            <div class="dashboard-card card-export" onclick="exportTableToCSV('vir-table', 'virulence_genes.csv')"><div style="font-size:2em;color:var(--export-color);"><i class="fas fa-virus"></i></div><div class="card-label">Virulence Genes CSV</div></div>
            <div class="dashboard-card card-export" onclick="exportTableToCSV('bac-table', 'bacmet_genes.csv')"><div style="font-size:2em;color:var(--export-color);"><i class="fas fa-flask"></i></div><div class="card-label">BACMET Genes CSV</div></div>
            <div class="dashboard-card card-export" onclick="exportTableToCSV('plasmid-table', 'plasmid_replicons.csv')"><div style="font-size:2em;color:var(--export-color);"><i class="fas fa-plug"></i></div><div class="card-label">Plasmid Replicons CSV</div></div>
            <div class="dashboard-card card-export" onclick="exportTableToCSV('mutation-table', 'mutations.csv')"><div style="font-size:2em;color:var(--export-color);"><i class="fas fa-dna"></i></div><div class="card-label">Mutations CSV</div></div>
            <div class="dashboard-card card-export" onclick="exportTableToCSV('qc-table', 'fasta_qc.csv')"><div style="font-size:2em;color:var(--export-color);"><i class="fas fa-chart-line"></i></div><div class="card-label">FASTA QC CSV</div></div>
            <div class="dashboard-card card-export" onclick="location.href='staphscope_ultimate_report.json'"><div style="font-size:2em;color:var(--export-color);"><i class="fas fa-file-code"></i></div><div class="card-label">Complete JSON Data</div></div>
        </div>
        """


# -----------------------------------------------------------------------------
# MAIN REPORTER CLASS
# -----------------------------------------------------------------------------
class StaphUltimateReporter:
    def __init__(self, input_dir: Path):
        self.input_dir = Path(input_dir)
        self.output_dir = self.input_dir / "STAPHSCOPE_ULTIMATE_SAMPLE_CENTRIC_REPORTS"
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.parser = StaphHTMLParser()
        self.analyzer = StaphDataAnalyzer()
        self.html_generator = StaphHTMLGenerator(self.analyzer)
        self.metadata = {
            "tool_name": "STAPHSCOPE Ultimate S. aureus Reporter",
            "version": "1.0.0",
            "author": "Brown Beckley <brownbeckley94@gmail.com>",
            "affiliation": "University of Ghana Medical School",
            "analysis_date": datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            "input_directory": str(self.input_dir)
        }

    def find_html_files(self) -> Dict[str, List[Path]]:
        print("🔍 Searching for StaphScope HTML reports (for QC, Mutations, fallback)...")
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
                print(f"  📁 {file_type.upper()}: {len(files)} files found")
        return html_files

    def integrate_all_data(self, html_files: Dict[str, List[Path]]) -> Dict[str, Any]:
        print("\n🔗 Integrating data – primary source: TSV files.")
        integrated_data = {
            'metadata': self.metadata,
            'samples': {},
            'patterns': {},
            'gene_centric': {},
            'qc_data': {},
            'amrfinder_details': {},
            'abricate_details': {},
            'mutation_details': {},   # per‑sample mutations
            'agr_data': {}
        }

        if html_files['qc']:
            integrated_data['qc_data'] = self.parser.load_qc_from_html(html_files['qc'][0])
        else:
            print("  ⚠️ QC HTML not found; QC tab will be empty.")

        # Load agr data
        agr_data = self.parser.load_agr_from_tsv(self.input_dir)
        if agr_data:
            integrated_data['agr_data'] = agr_data
            print(f"  ✅ Loaded agr typing for {len(agr_data)} samples")
        else:
            print("  ⚠️ agr_summary.tsv not found; agr tab will be empty.")

        # Load mutations – only per‑sample details
        mutation_tsv_data = self.parser.load_mutations_from_tsv(self.input_dir)
        if mutation_tsv_data:
            integrated_data['mutation_details'] = mutation_tsv_data
            print(f"  ✅ Loaded per‑sample mutations for {len(mutation_tsv_data)} samples")
        else:
            print("  ⚠️ No mutation data found; mutation tab will be empty.")

        typing_data = self.parser.load_typing_from_tsv(self.input_dir)
        amr_details, amr_gene_freq = self.parser.load_amrfinder_from_tsv(self.input_dir)
        abricate_details, abricate_gene_freq = self.parser.load_abricate_from_tsv(self.input_dir)

        if not typing_data and not amr_details and not abricate_details:
            print("  ⚠️ No TSV data found; falling back to HTML parsing.")
            if html_files['comprehensive']:
                typing_data = self.parser.parse_comprehensive_report(html_files['comprehensive'][0])
            if html_files['amrfinder']:
                amr_by_sample, amr_gene_freq = self.parser.parse_amrfinder_report(html_files['amrfinder'][0])
                amr_details = {}
                for sample, data in amr_by_sample.items():
                    amr_details[sample] = [{'gene': g} for g in data.get('all_genes', [])]
            if html_files['abricate']:
                abricate_details = defaultdict(lambda: defaultdict(list))
                abricate_gene_freq = defaultdict(dict)
                for abricate_file in html_files['abricate']:
                    db_name, genes_by_sample, gene_freq = self.parser.parse_abricate_report(abricate_file)
                    if db_name != 'unknown':
                        for sample, genes in genes_by_sample.items():
                            abricate_details[sample][db_name] = [{'gene': g} for g in genes]
                        abricate_gene_freq[db_name] = gene_freq
            if not typing_data and html_files['comprehensive']:
                typing_data = self.parser.parse_comprehensive_report(html_files['comprehensive'][0])

        all_samples = set(typing_data.keys()) | set(amr_details.keys()) | set(abricate_details.keys())
        all_samples.update(integrated_data['qc_data'].keys())
        all_samples.update(agr_data.keys())
        all_samples = sorted(list(all_samples))
        if not all_samples:
            print("❌ No samples found in any source!")
            return {}

        print(f"\n📊 Found {len(all_samples)} unique samples")

        for sample in all_samples:
            typing_info = typing_data.get(sample, {'MLST': 'ND', 'spa_Type': 'ND', 'SCCmec_Type': 'ND', 'MRSA_Status': 'ND'})
            agr_info = agr_data.get(sample, {
                'agr_Type': 'ND',
                'agr_Group': 'ND',
                'match_score': '',
                'canonical_agrD': '',
                'multiple_agr': '',
                'status': 'not_found'
            })
            amr_list = amr_details.get(sample, [])
            amr_gene_names = [d['gene'] for d in amr_list]
            amr_info = {
                'critical_genes': [],
                'high_risk_genes': [],
                'all_genes': amr_gene_names
            }
            for gene in amr_gene_names:
                if gene.lower() in self.analyzer.critical_amr_genes:
                    amr_info['critical_genes'].append(gene)
                if gene.lower() in self.analyzer.high_priority_amr:
                    amr_info['high_risk_genes'].append(gene)
            abricate_info = {}
            for db, genes in abricate_details.get(sample, {}).items():
                abricate_info[db] = [g['gene'] for g in genes]
            integrated_data['samples'][sample] = {
                'typing': typing_info,
                'agr': agr_info,
                'amrfinder': amr_info,
                'abricate_databases': abricate_info
            }

        integrated_data['amrfinder_details'] = amr_details
        integrated_data['abricate_details'] = abricate_details

        integrated_data['gene_frequencies'] = {
            'amrfinder': amr_gene_freq,
            'abricate': abricate_gene_freq
        }

        print("\n🧠 Processing gene‑centric and pattern analysis...")
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
            agr = data.get('agr', {})
            samples_data.append({
                'Sample': sample,
                'MLST': data['typing']['MLST'],
                'spa_Type': data['typing']['spa_Type'],
                'SCCmec_Type': data['typing']['SCCmec_Type'],
                'MRSA_Status': data['typing']['MRSA_Status'],
                'agr_Type': agr.get('agr_Type', 'ND'),
                'agr_Group': agr.get('agr_Group', 'ND'),
                'agr_Status': agr.get('status', 'not_found'),
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

        # Export per‑sample mutation details as CSV (if any)
        mutation_details = integrated_data.get('mutation_details', {})
        if mutation_details:
            mut_rows = []
            for sample, muts in mutation_details.items():
                for m in muts:
                    mut_rows.append({
                        'Sample': sample,
                        'Gene': m.get('gene', ''),
                        'Mutation': m.get('mutation', ''),
                        'Class': m.get('class', ''),
                        'Subclass': m.get('subclass', ''),
                        'Contig': m.get('contig', ''),
                        'Start': m.get('start', ''),
                        'Stop': m.get('stop', ''),
                        'Strand': m.get('strand', ''),
                        'Coverage': m.get('coverage', ''),
                        'Identity': m.get('identity', ''),
                        'Accession': m.get('accession', '')
                    })
            if mut_rows:
                pd.DataFrame(mut_rows).to_csv(self.output_dir / "mutations.csv", index=False)

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

        print(f"    ✅ CSV reports generated.")

    def run(self):
        print("=" * 80)
        print("🧬 STAPHSCOPE ULTIMATE S. AUREUS REPORTER v1.0.0")
        print("   (Hybrid: Gene‑centric for typing, Interactive Sample‑Centric for AMR/Virulence/BACMET/Plasmids/Mutations)")
        print("   (Primary data source: TSV summaries)")
        print("=" * 80)
        print(f"📁 Input directory: {self.input_dir}")

        html_files = self.find_html_files()
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
        agr_data = integrated_data.get('agr_data', {})
        mutation_samples = len(integrated_data.get('mutation_details', {}))
        print("\n" + "=" * 80)
        print("✅ ULTIMATE ANALYSIS COMPLETE!")
        print("=" * 80)
        print(f"📁 Output directory: {self.output_dir}")
        print(f"📄 Files generated:")
        print(f"   • staphscope_ultimate_report.html (Hybrid report with interactive isolate boxes and agr typing)")
        print(f"   • staphscope_ultimate_report.json (Complete data)")
        print(f"\n🔬 Interactive Isolate Boxes with Typing Badges including agr")
        print(f"   • Each isolate box shows MLST, spa, SCCmec, MRSA/MSSA, and agr type.")
        print(f"   • Horizontally scrollable tables – no truncation.")
        print(f"   • Filter by sample name or database.")
        print(f"   • Agr typing distribution tab available.")
        print(f"   • Mutations tab now shows per‑isolate boxes with full mutation details.")
        print(f"\n📈 ANALYSIS SUMMARY:")
        print(f"   • {total_samples} total S. aureus samples analyzed")
        print(f"   • {len(patterns.get('mlst_distribution', {}))} unique STs")
        print(f"   • {len(patterns.get('spa_type_distribution', {}))} unique spa types")
        print(f"   • {len(patterns.get('sccmec_distribution', {}))} SCCmec types")
        print(f"   • {len(agr_data)} samples with agr typing")
        print(f"   • {mutation_samples} samples with point mutations")
        print(f"   • {len(patterns.get('high_risk_combinations', []))} high‑risk AMR+virulence combos")
        print("\n🎯 Next steps:")
        print("   1. Open staphscope_ultimate_report.html in your browser")
        print("   2. Explore the interactive isolate boxes in AMR/Virulence/BACMET/Plasmids tabs")
        print("   3. Check the agr typing tab for distribution")
        print("   4. Use the new Mutations tab to see per‑isolate mutation tables")
        print("   5. Use the filters to focus on specific samples or databases")
        print("   6. Export any table as CSV for further analysis")
        print("\n" + "=" * 80)
        return True


def main():
    parser = argparse.ArgumentParser(
        description='STAPHSCOPE Ultimate S. aureus Reporter v1.0.0 - Hybrid Gene‑centric / Sample‑centric with agr',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""Examples:\n  python staphscope_ultimate_samplecentric_reporter.py -i /path/to/staphscope/reports\n\nAuthor: Brown Beckley <brownbeckley94@gmail.com>\nAffiliation: University of Ghana Medical School"""
    )
    parser.add_argument('-i', '--input-dir', required=True,
                        help='Directory containing StaphScope TSV summaries and/or HTML reports')
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
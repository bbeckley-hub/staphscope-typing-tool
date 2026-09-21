#!/usr/bin/env python3
"""
STAPHSCOPE ULTIMATE REPORTER - GENE-CENTRIC S. AUREUS ANALYSIS
===================================================================
Version 3.0.0

Gene-centric cross-genome analysis with dynamic grouping by typing.
Single master TSV for all typing (MLST, spa, agr, capsule, SCCmec CGE/RPet/Subtype, MRSA).
Gene-centric tabs (AMR, Virulence, BACMET, Plasmids, Mutations) with dropdown grouping,
filter buttons, gene-family info boxes, and rich acknowledgements.
New: Capsule tab, merged SCCmec tab, MGE tab, comparison tab.

Find a superbug? Say Hi or Open an ISSUE!!!!!

Author: Brown Beckley <brownbeckley94@gmail.com>
Affiliation: University of Ghana Medical School
MIT
"""

import os
import sys
import json
import re
import io
import glob
import argparse
import pandas as pd
from pathlib import Path
from typing import Dict, List, Tuple, Any, Optional
from datetime import datetime
from collections import defaultdict, Counter
import warnings
warnings.filterwarnings('ignore')

from bs4 import BeautifulSoup


def esc(v) -> str:
    """Escape a value for safe HTML embedding."""
    if v is None:
        return ""
    return (str(v).replace("&", "&amp;").replace("<", "&lt;")
            .replace(">", "&gt;").replace('"', "&quot;"))


# =============================================================================
# PARSER
# =============================================================================
class StaphHTMLParser:
    """Parses StaphScope HTML/TSV outputs into structured dictionaries."""

    def __init__(self):
        self.abricate_databases = [
            'card', 'resfinder', 'vfdb', 'argannot',
            'plasmidfinder', 'megares', 'ncbi', 'bacmet2'
        ]

    def normalize_sample_id(self, sample_id: str) -> str:
        """Strip file extension and directory path from a sample identifier."""
        sample = str(sample_id)
        for ext in ('.fna', '.fasta', '.fa', '.gb', '.gbk', '.gbff', '.txt', '.tsv', '.csv'):
            if sample.endswith(ext):
                sample = sample[:-len(ext)]
        if '/' in sample or '\\' in sample:
            sample = Path(sample).name
        return sample.strip()

    def parse_html_table(self, html_content: str, table_index: int = 0) -> pd.DataFrame:
        """Parse the Nth HTML table in a string into a DataFrame."""
        try:
            soup = BeautifulSoup(html_content, 'html.parser')
            tables = soup.find_all('table')
            if not tables or table_index >= len(tables):
                return pd.DataFrame()
            table = tables[table_index]
            rows = table.find_all('tr')
            if not rows:
                return pd.DataFrame()
            headers = [th.get_text().strip() for th in rows[0].find_all(['th', 'td'])]
            data = []
            for row in rows[1:]:
                cols = row.find_all(['td', 'th'])
                if cols and len(cols) == len(headers):
                    data.append([c.get_text().strip() for c in cols])
            return pd.DataFrame(data, columns=headers) if data else pd.DataFrame()
        except Exception as e:
            print(f"  ⚠️ Table parsing error: {e}")
            return pd.DataFrame()

    def load_master_tsv(self, path: Path) -> Dict[str, Dict]:
        """Load the master typing TSV that aggregates every typing field."""
        print(f"  🧬 Loading master typing TSV: {path.name}")
        if not path.exists():
            print(f"    ❌ File not found: {path}")
            return {}
        try:
            df = pd.read_csv(path, sep='\t', dtype=str).fillna('Not Assigned')
            results = {}
            for _, row in df.iterrows():
                sample_raw = row.get('Sample', '')
                if not sample_raw or sample_raw == 'Not Assigned':
                    continue
                sample = self.normalize_sample_id(sample_raw)
                results[sample] = {
                    'MLST':           str(row.get('MLST', 'Not Assigned')).strip(),
                    'spa_Type':       str(row.get('spa Type', 'Not Assigned')).strip(),
                    'agr_Type':       str(row.get('agr Type', 'Not Assigned')).strip(),
                    'capsule_type':   str(row.get('Capsule Type', 'Not Assigned')).strip(),
                    'SCCmec_CGE':     str(row.get('SCCmec Type (CGE)', 'Not Assigned')).strip(),
                    'SCCmec_RPet':    str(row.get('SCCmec Type (RPet)', 'Not Assigned')).strip(),
                    'SCCmec_Subtype': str(row.get('SCCmec Subtype', 'Not Assigned')).strip(),
                    'MRSA_Status':    str(row.get('MRSA/MSSA Status', 'Not Assigned')).strip(),
                }
            print(f"    ✓ Loaded {len(results)} samples")
            return results
        except Exception as e:
            print(f"    ❌ Error reading master TSV: {e}")
            return {}

    def parse_qc_report(self, file_path: Path) -> Dict[str, Dict]:
        """Parse the FASTA QC HTML summary into per-sample QC metrics."""
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
                        except Exception:
                            qc_data[col] = str(val)
                results[sample] = qc_data
            print(f"    ✓ Parsed {len(results)} samples")
            return results
        except Exception as e:
            print(f"    ❌ Error parsing QC: {e}")
            return {}

    def parse_amrfinder_report(self, file_path: Path) -> Tuple[Dict[str, List], Dict[str, Dict]]:
        """Parse the AMRFinderPlus HTML summary."""
        print(f"  🧬 Parsing AMRfinder: {file_path.name}")
        try:
            with open(file_path, 'r', encoding='utf-8') as f:
                html_content = f.read()
            soup = BeautifulSoup(html_content, 'html.parser')
            tables = soup.find_all('table')
            genes_by_genome, gene_frequencies = {}, {}
            genes_by_genome_table = gene_frequency_table = None
            for table in tables:
                t = table.get_text()
                if 'Genome' in t and 'Critical Genes' in t:
                    genes_by_genome_table = table
                elif 'Gene' in t and 'Frequency' in t and 'Prevalence' in t:
                    gene_frequency_table = table

            if genes_by_genome_table:
                try:
                    df_genomes = pd.read_html(io.StringIO(str(genes_by_genome_table)))[0]
                    genome_col = next((c for c in df_genomes.columns
                                       if 'genome' in c.lower()), df_genomes.columns[0])
                    crit_col = next((c for c in df_genomes.columns
                                     if 'critical' in c.lower()), None)
                    hr_col = next((c for c in df_genomes.columns
                                   if 'high risk' in c.lower() or 'high' in c.lower()), None)
                    for _, row in df_genomes.iterrows():
                        sample = self.normalize_sample_id(row[genome_col])
                        crit, hr = [], []
                        if crit_col and pd.notna(row.get(crit_col)):
                            v = str(row[crit_col]).strip()
                            if v.lower() not in ('none', 'nan', ''):
                                crit = [g.strip() for g in v.split(',') if g.strip()]
                        if hr_col and pd.notna(row.get(hr_col)):
                            v = str(row[hr_col]).strip()
                            if v.lower() not in ('none', 'nan', ''):
                                hr = [g.strip() for g in v.split(',') if g.strip()]
                        genes_by_genome[sample] = {
                            'critical_genes': crit,
                            'high_risk_genes': hr,
                            'all_genes': []
                        }
                except Exception as e:
                    print(f"    Error parsing genes-by-genome table: {e}")

            if gene_frequency_table:
                try:
                    df_genes = pd.read_html(io.StringIO(str(gene_frequency_table)))[0]
                    gene_col = next((c for c in df_genes.columns if 'gene' in c.lower()), None)
                    freq_col = next((c for c in df_genes.columns if 'freq' in c.lower()), None)
                    prev_col = next((c for c in df_genes.columns if 'prev' in c.lower()), None)
                    risk_col = next((c for c in df_genes.columns if 'risk' in c.lower()), None)
                    genomes_col = next((c for c in df_genes.columns if 'genome' in c.lower()), None)
                    if gene_col:
                        for _, row in df_genes.iterrows():
                            if pd.isna(row[gene_col]):
                                continue
                            gene = str(row[gene_col]).strip()
                            frequency = str(row[freq_col]).strip() if freq_col and pd.notna(row.get(freq_col)) else '0'
                            prevalence = str(row[prev_col]).strip() if prev_col and pd.notna(row.get(prev_col)) else 'ND'
                            risk = str(row[risk_col]).strip() if risk_col and pd.notna(row.get(risk_col)) else 'STANDARD'
                            genomes = []
                            if genomes_col and pd.notna(row.get(genomes_col)):
                                gstr = str(row[genomes_col]).replace('tet(38)', 'tet38')
                                genomes = [self.normalize_sample_id(g.replace('tet38', 'tet(38)'))
                                           for g in gstr.split(',') if g.strip()]
                            match = re.search(r'(\d+)', frequency)
                            count = int(match.group(1)) if match else 0
                            gene_frequencies[gene] = {
                                'frequency': frequency,
                                'count': count,
                                'prevalence': prevalence,
                                'risk_level': risk,
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
                except Exception as e:
                    print(f"    Error parsing gene frequency table: {e}")

            print(f"    ✓ Found {len(genes_by_genome)} samples, {len(gene_frequencies)} genes")
            return genes_by_genome, gene_frequencies
        except Exception as e:
            print(f"    ❌ Error parsing AMRfinder: {e}")
            return {}, {}

    def parse_abricate_report(self, file_path: Path) -> Tuple[str, Dict[str, List], Dict[str, Dict]]:
        """Parse an ABRicate HTML summary for one database."""
        print(f"  🧬 Parsing ABRicate: {file_path.name}")
        try:
            with open(file_path, 'r', encoding='utf-8') as f:
                html_content = f.read()
            soup = BeautifulSoup(html_content, 'html.parser')
            tables = soup.find_all('table')
            if len(tables) < 2:
                return 'unknown', {}, {}
            db_name = 'unknown'
            fname = file_path.name.lower()
            for db in self.abricate_databases:
                if db in fname:
                    db_name = db
                    break
            title = soup.find('title')
            if title:
                tl = title.get_text().lower()
                for db in self.abricate_databases:
                    if db in tl:
                        db_name = db
                        break

            genes_by_genome = {}
            df1 = self.parse_html_table(html_content, 0)
            if not df1.empty:
                df1.columns = [c.strip() for c in df1.columns]
                genome_col = next((c for c in df1.columns
                                   if 'genome' in c.lower() or 'sample' in c.lower()),
                                  df1.columns[0])
                genes_col = next((c for c in df1.columns
                                  if 'genes' in c.lower() or 'detected' in c.lower()), None)
                if genome_col and genes_col:
                    for _, row in df1.iterrows():
                        sample = self.normalize_sample_id(row[genome_col])
                        if pd.notna(row.get(genes_col)):
                            genes_by_genome[sample] = [
                                g.strip() for g in str(row[genes_col]).split(',') if g.strip()]

            gene_frequencies = {}
            df2 = self.parse_html_table(html_content, 1)
            if not df2.empty:
                df2.columns = [c.strip() for c in df2.columns]
                if 'Gene' in df2.columns:
                    for _, row in df2.iterrows():
                        gene = str(row['Gene']).strip()
                        frequency = str(row.get('Frequency', '0')).strip()
                        genomes = []
                        if 'Genomes' in df2.columns and pd.notna(row.get('Genomes')):
                            genomes = [self.normalize_sample_id(g.strip())
                                       for g in str(row['Genomes']).split(',') if g.strip()]
                        match = re.search(r'(\d+)', frequency)
                        count = int(match.group(1)) if match else 0
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
        """Parse the mutation summary HTML for point-mutation data."""
        print(f"  🧬 Parsing mutation summary: {file_path.name}")
        try:
            with open(file_path, 'r', encoding='utf-8') as f:
                html_content = f.read()
            soup = BeautifulSoup(html_content, 'html.parser')
            mutation_table = None
            for table in soup.find_all('table'):
                if (table.find(string=re.compile(r'Gene', re.I)) and
                        table.find(string=re.compile(r'Mutation', re.I))):
                    mutation_table = table
                    break
            if not mutation_table:
                return {}
            header_row = None
            thead = mutation_table.find('thead')
            if thead:
                header_row = thead.find('tr')
            if not header_row:
                header_row = mutation_table.find('tr')
            if not header_row:
                return {}
            headers = [c.get_text().strip()
                       for c in header_row.find_all(['th', 'td'])
                       if c.get_text().strip()]
            col_idx = {}
            for idx, h in enumerate(headers):
                hl = h.lower()
                if 'gene' in hl: col_idx['gene'] = idx
                elif 'mutation' in hl: col_idx['mutation'] = idx
                elif 'count' in hl: col_idx['count'] = idx
                elif 'genome' in hl: col_idx['genomes'] = idx
                elif 'class' in hl: col_idx['class'] = idx
                elif 'subclass' in hl: col_idx['subclass'] = idx
            for req in ['gene', 'mutation', 'count', 'genomes']:
                if req not in col_idx:
                    return {}
            tbody = mutation_table.find('tbody')
            rows = tbody.find_all('tr') if tbody else mutation_table.find_all('tr')[1:]
            mutations_list, genome_counts = [], defaultdict(int)
            for row in rows:
                cells = row.find_all('td')
                if len(cells) <= max(col_idx.values()):
                    continue
                gene = cells[col_idx['gene']].get_text().strip()
                mutation = cells[col_idx['mutation']].get_text().strip()
                count_str = cells[col_idx['count']].get_text().strip()
                m = re.search(r'(\d+)', count_str)
                count = int(m.group(1)) if m else 0
                genomes = [g.strip() for g in
                           cells[col_idx['genomes']].get_text().strip().split(',') if g.strip()]
                if not genomes:
                    continue
                for g in genomes:
                    genome_counts[g] += 1
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
            print(f"    ✓ {len(mutations_list)} mutations across {len(genome_counts)} genomes")
            return {
                'mutations': mutations_list,
                'genome_mutation_counts': dict(genome_counts)
            }
        except Exception as e:
            print(f"    ❌ Error parsing mutation summary: {e}")
            return {}

    def parse_mge_summary_html(self, file_path: Path) -> Dict[str, Any]:
        """Parse the mobileOG MGE summary HTML (stat cards + per-sample table)."""
        print(f"  🧬 Parsing MGE summary: {file_path.name}")
        try:
            with open(file_path, 'r', encoding='utf-8') as f:
                html_content = f.read()
            soup = BeautifulSoup(html_content, 'html.parser')

            stats = {}
            for label in ('Total Samples', 'With mobileOG hits', 'Without mobileOG hits',
                          'Total mobileOG hits', 'Key MGE-associated genes', 'Runtime'):
                for tag in soup.find_all(string=re.compile(re.escape(label), re.I)):
                    parent = tag.find_parent()
                    text_block = parent.get_text(" ", strip=True) if parent else ''
                    m = re.search(rf'{re.escape(label)}\s*[:]?\s*([\d,\.]+s?)\b',
                                  text_block, re.I)
                    if m:
                        stats[label] = m.group(1)
                        break
                    if parent:
                        for sib in parent.find_all_next(limit=3):
                            s = sib.get_text(strip=True)
                            if re.fullmatch(r'[\d,\.]+s?', s):
                                stats[label] = s
                                break
                    if label in stats:
                        break

            per_sample = {}
            mge_table = None
            for table in soup.find_all('table'):
                header_text = ' '.join(th.get_text().strip()
                                       for th in table.find_all('th'))
                if 'Sample' in header_text and 'mobileOG' in header_text:
                    mge_table = table
                    break
            if mge_table:
                try:
                    df = pd.read_html(io.StringIO(str(mge_table)))[0]
                    df.columns = [str(c).strip() for c in df.columns]
                    for _, row in df.iterrows():
                        sample = self.normalize_sample_id(row.get('Sample', ''))
                        if not sample:
                            continue
                        entry = {}
                        for col in df.columns:
                            if col == 'Sample':
                                continue
                            v = row.get(col, 0)
                            try:
                                entry[col] = int(str(v).replace(',', '').strip())
                            except Exception:
                                entry[col] = v
                        per_sample[sample] = entry
                except Exception as e:
                    print(f"    Error parsing MGE table: {e}")

            note = ""
            for tag in soup.find_all(string=re.compile(r'Interpretation note', re.I)):
                parent = tag.find_parent()
                if parent:
                    raw = parent.get_text(" ", strip=True)
                    raw = re.sub(r'^Interpretation note\s*:?\s*', '', raw, flags=re.I)
                    if len(raw) > 40:
                        note = raw
                        break

            print(f"    ✓ MGE stats: {len(stats)} cards, {len(per_sample)} samples")
            return {'stats': stats, 'per_sample': per_sample, 'note': note}
        except Exception as e:
            print(f"    ❌ Error parsing MGE summary: {e}")
            return {'stats': {}, 'per_sample': {}, 'note': ''}


# =============================================================================
# ANALYZER
# =============================================================================
class StaphDataAnalyzer:
    """Cross-genome patterns, gene-centric tables, MGE statistics."""

    def __init__(self):
        self.critical_amr_genes = {
            'meca', 'mecc', 'vana', 'vanb', 'vanc',
            'erma', 'ermb', 'ermc', 'msra', 'mphc',
            'tetk', 'tetm', 'tetl'
        }
        self.critical_virulence_genes = {
            'luks-pv', 'lukf-pv', 'tsst', 'sea', 'seb', 'sec', 'sed', 'see',
            'seg', 'seh', 'sei', 'sej', 'sek', 'sel', 'sem', 'sen', 'seo', 'sep',
            'seq', 'ser', 'seu', 'eta', 'etb', 'hla', 'hlb', 'hlg', 'hld',
        }

    def create_gene_centric_tables(self, integrated_data: Dict[str, Any]) -> Dict[str, Any]:
        """Group gene frequencies by database category for gene-centric tables."""
        gene_centric = {
            'amr_databases': {},
            'virulence_databases': {},
            'plasmid_databases': {},
            'bacmet_databases': {},
            'combined_gene_frequencies': []
        }

        if 'amrfinder' in integrated_data.get('gene_frequencies', {}):
            amr_data = integrated_data['gene_frequencies']['amrfinder']
            gene_centric['amr_databases']['amrfinder'] = sorted(
                [{'gene': g, 'database': 'AMRfinder',
                  'frequency': d.get('frequency', '0'),
                  'count': d.get('count', 0),
                  'prevalence': d.get('prevalence', 'ND'),
                  'risk_level': d.get('risk_level', 'ND'),
                  'genomes': d.get('genomes', [])}
                 for g, d in amr_data.items()],
                key=lambda x: x['count'], reverse=True)

        if 'abricate' in integrated_data.get('gene_frequencies', {}):
            for db_name, db_genes in integrated_data['gene_frequencies']['abricate'].items():
                gene_list = [{'gene': g, 'database': db_name.upper(),
                              'frequency': d.get('frequency', '0'),
                              'count': d.get('count', 0),
                              'genomes': d.get('genomes', [])}
                             for g, d in db_genes.items()]
                if not gene_list:
                    continue
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
        for db_type in ('amr_databases', 'virulence_databases',
                        'plasmid_databases', 'bacmet_databases'):
            for genes in gene_centric.get(db_type, {}).values():
                all_genes.extend(genes)
        all_genes.sort(key=lambda x: x['count'], reverse=True)
        gene_centric['combined_gene_frequencies'] = all_genes
        return gene_centric

    def create_cross_genome_patterns(self, integrated_data: Dict[str, Any]) -> Dict[str, Any]:
        """Build every combination table and distribution for typing fields."""
        samples_data = integrated_data.get('samples', {})

        P = {
            'mlst_distribution': Counter(),
            'spa_type_distribution': Counter(),
            'agr_type_distribution': Counter(),
            'capsule_distribution': Counter(),
            'sccmec_cge_distribution': Counter(),
            'sccmec_rpet_distribution': Counter(),
            'sccmec_subtype_distribution': Counter(),
            'mrsa_status_distribution': Counter(),
            'high_risk_combinations': [],
            'gene_cooccurrence': defaultdict(Counter),
        }

        COMBO_FIELDS = {
            'mlst_spa':              ['MLST', 'spa_Type'],
            'mlst_sccmec_cge':       ['MLST', 'SCCmec_CGE'],
            'mlst_sccmec_rpet':      ['MLST', 'SCCmec_RPet'],
            'mlst_subtype':          ['MLST', 'SCCmec_Subtype'],
            'mlst_agr':              ['MLST', 'agr_Type'],
            'mlst_capsule':          ['MLST', 'capsule_type'],
            'spa_sccmec_cge':        ['spa_Type', 'SCCmec_CGE'],
            'spa_sccmec_rpet':       ['spa_Type', 'SCCmec_RPet'],
            'spa_subtype':           ['spa_Type', 'SCCmec_Subtype'],
            'spa_agr':               ['spa_Type', 'agr_Type'],
            'spa_capsule':           ['spa_Type', 'capsule_type'],
            'agr_capsule':           ['agr_Type', 'capsule_type'],
            'agr_sccmec_cge':        ['agr_Type', 'SCCmec_CGE'],
            'agr_sccmec_rpet':       ['agr_Type', 'SCCmec_RPet'],
            'agr_subtype':           ['agr_Type', 'SCCmec_Subtype'],
            'capsule_sccmec_cge':    ['capsule_type', 'SCCmec_CGE'],
            'capsule_sccmec_rpet':   ['capsule_type', 'SCCmec_RPet'],
            'capsule_subtype':       ['capsule_type', 'SCCmec_Subtype'],
            'mlst_spa_agr':          ['MLST', 'spa_Type', 'agr_Type'],
            'mlst_spa_sccmec_cge':   ['MLST', 'spa_Type', 'SCCmec_CGE'],
            'mlst_spa_subtype':      ['MLST', 'spa_Type', 'SCCmec_Subtype'],
            'mlst_capsule_agr':      ['MLST', 'capsule_type', 'agr_Type'],
            'spa_capsule_agr':       ['spa_Type', 'capsule_type', 'agr_Type'],
            'mlst_spa_agr_capsule':  ['MLST', 'spa_Type', 'agr_Type', 'capsule_type'],
            'mlst_spa_sccmec_agr':   ['MLST', 'spa_Type', 'SCCmec_CGE', 'agr_Type'],
            'mlst_spa_subtype_agr':  ['MLST', 'spa_Type', 'SCCmec_Subtype', 'agr_Type'],
        }

        combo_dicts = {k: defaultdict(list) for k in COMBO_FIELDS}

        sample_genes = defaultdict(list)
        gene_centric = integrated_data.get('gene_centric', {})
        for db_type in ('amr_databases', 'virulence_databases'):
            for genes in gene_centric.get(db_type, {}).values():
                for g in genes:
                    for genome in g['genomes']:
                        if g['gene'] not in sample_genes[genome]:
                            sample_genes[genome].append(g['gene'])

        for sample, data in samples_data.items():
            typing = data.get('typing', {})
            values = {
                'MLST':            typing.get('MLST', 'Not Assigned'),
                'spa_Type':        typing.get('spa_Type', 'Not Assigned'),
                'agr_Type':        typing.get('agr_Type', 'Not Assigned'),
                'capsule_type':    typing.get('capsule_type', 'Not Assigned'),
                'SCCmec_CGE':      typing.get('SCCmec_CGE', 'Not Assigned'),
                'SCCmec_RPet':     typing.get('SCCmec_RPet', 'Not Assigned'),
                'SCCmec_Subtype':  typing.get('SCCmec_Subtype', 'Not Assigned'),
                'MRSA_Status':     typing.get('MRSA_Status', 'Not Assigned'),
            }

            def ok(v):
                return v and v not in ('Not Assigned', 'ND', '', 'nan')

            if ok(values['MLST']):         P['mlst_distribution'][values['MLST']] += 1
            if ok(values['spa_Type']):     P['spa_type_distribution'][values['spa_Type']] += 1
            if ok(values['agr_Type']):     P['agr_type_distribution'][values['agr_Type']] += 1
            if ok(values['capsule_type']): P['capsule_distribution'][values['capsule_type']] += 1
            if ok(values['SCCmec_CGE']):   P['sccmec_cge_distribution'][values['SCCmec_CGE']] += 1
            if ok(values['SCCmec_RPet']):  P['sccmec_rpet_distribution'][values['SCCmec_RPet']] += 1
            if ok(values['SCCmec_Subtype']): P['sccmec_subtype_distribution'][values['SCCmec_Subtype']] += 1
            if ok(values['MRSA_Status']):  P['mrsa_status_distribution'][values['MRSA_Status']] += 1

            for key, fields in COMBO_FIELDS.items():
                if all(ok(values[f]) for f in fields):
                    label = ' - '.join(values[f] for f in fields)
                    combo_dicts[key][label].append(sample)

            amr_genes = data.get('amrfinder', {}).get('all_genes', [])
            vir_genes = data.get('abricate_databases', {}).get('vfdb', [])
            critical_amr = [g for g in amr_genes
                            if any(c in str(g).lower() for c in self.critical_amr_genes)]
            critical_vir = [g for g in vir_genes
                            if any(c in str(g).lower() for c in self.critical_virulence_genes)]
            if critical_amr and critical_vir:
                P['high_risk_combinations'].append({
                    'sample': sample,
                    'mlst': values['MLST'],
                    'spa_type': values['spa_Type'],
                    'sccmec_type': values['SCCmec_CGE'],
                    'mrsa_status': values['MRSA_Status'],
                    'agr_type': values['agr_Type'],
                    'critical_amr_genes': critical_amr,
                    'critical_virulence_genes': critical_vir,
                })

            genes = sample_genes.get(sample, [])
            for i, g1 in enumerate(genes):
                for g2 in genes[i + 1:]:
                    P['gene_cooccurrence'][g1][g2] += 1

        for key in combo_dicts:
            P[key] = dict(combo_dicts[key])
        return P

    def compute_mge_stats(self, mge_data: Dict[str, Any]) -> Dict[str, Any]:
        """Aggregate MGE per-sample counts into category totals and means."""
        per_sample = mge_data.get('per_sample', {})
        if not per_sample:
            return {'category_totals': [], 'total_per_sample': 0}
        all_columns = set()
        for s in per_sample.values():
            all_columns.update(s.keys())
        numeric_cols = [c for c in all_columns
                        if all(isinstance(s.get(c, 0), (int, float))
                               for s in per_sample.values())]
        category_totals = []
        n = len(per_sample)
        grand_total = 0
        for col in numeric_cols:
            vals = {s: per_sample[s].get(col, 0) for s in per_sample}
            total = sum(vals.values())
            grand_total += total
            max_sample = max(vals, key=vals.get) if vals else 'ND'
            category_totals.append({
                'category': col,
                'total': total,
                'mean': round(total / n, 2) if n else 0,
                'max': vals[max_sample] if vals else 0,
                'max_sample': max_sample,
            })
        category_totals.sort(key=lambda x: x['total'], reverse=True)
        for c in category_totals:
            c['pct_of_all'] = round((c['total'] / grand_total) * 100, 2) if grand_total else 0
        return {
            'category_totals': category_totals,
            'grand_total': grand_total,
            'n_samples': n,
        }


# =============================================================================
# HTML GENERATOR
# =============================================================================
class StaphHTMLGenerator:
    """Builds the interactive multi-tab HTML report."""

    def __init__(self, analyzer: StaphDataAnalyzer):
        self.analyzer = analyzer
        self._current_samples_data = {}
        self.tab_colors = {
            'summary': '#4CAF50', 'sample_overview': '#2196F3', 'qc': '#607D8B',
            'mlst': '#FF9800', 'spa': '#9C27B0', 'sccmec': '#009688',
            'capsule': '#00ACC1', 'mrsa': '#795548', 'agr': '#8B5CF6',
            'amr': '#F44336', 'virulence': '#E91E63', 'bacmet': '#FF5722',
            'plasmids': '#673AB7', 'mutation': '#00BCD4', 'mge': '#16A085',
            'patterns': '#3F51B5', 'compare': '#0d9488', 'aiguide': '#00BCD4', 'calltoaction': '#F472B6',
            'citation': '#8BC34A', 'funding': '#FFC107', 'export': '#9E9E9E'
        }

    # -------------------------------------------------------------------------
    # Reusable HTML helpers
    # -------------------------------------------------------------------------
    def _credit_bar(self, color: str, icon: str, title: str, body: str) -> str:
        """Colored acknowledgment strip used at the top of every tool-driven tab."""
        return f'''
        <div class="scientific-note" style="background: linear-gradient(135deg, #f8f9fa 0%, #f0f4f8 100%); border-left: 6px solid {color}; margin-bottom: 20px; padding: 15px; border-radius: 8px;">
            <div style="display: flex; align-items: center; gap: 12px; flex-wrap: wrap;">
                <span style="font-size: 1.4em;">{icon}</span>
                <div>
                    <strong style="font-size: 1.1em; color: {color};">{title}</strong><br>
                    <span style="font-size: 0.95em; color: #333;">{body}</span>
                </div>
            </div>
        </div>'''

    def _alert(self, kind: str, icon: str, html_body: str) -> str:
        """Standard alert box. kind ∈ {info, success, warning, danger}."""
        return f'''
        <div class="alert-box alert-{kind}">
            <i class="fas {icon} fa-2x"></i>
            <div>{html_body}</div>
        </div>'''

    def _stat_card(self, value, label: str, color: str = '#4CAF50', icon: str = '') -> str:
        """Colored stat card used by the MGE tab and dashboard."""
        return f'''
        <div class="stat-card" style="background: linear-gradient(135deg, {color} 0%, {color}dd 100%);">
            {f'<i class="fas {icon} fa-2x" style="opacity:0.9; margin-bottom:8px;"></i>' if icon else ''}
            <div class="stat-value">{value}</div>
            <div class="stat-label">{label}</div>
        </div>'''

    def _filter_buttons(self, table_id: str, buttons: list) -> str:
        """Row of quick-filter buttons that populate the table's search box.

        buttons: list of (label, search_value, css_class, icon) tuples.
        """
        html = f'''<div class="action-buttons">
            <button class="action-btn btn-primary"
                onclick="exportTableToCSV('{table_id}', '{table_id}.csv')">
                <i class="fas fa-download"></i> Export</button>'''
        for btn in buttons:
            label, search_val = btn[0], btn[1]
            css_class = btn[2] if len(btn) > 2 else 'btn-info'
            icon = btn[3] if len(btn) > 3 else 'fa-filter'
            html += (f'''<button class="action-btn {css_class}"
                onclick="document.getElementById('search-{table_id}').value='{search_val}';'''
                     f'''searchTable('{table_id}','search-{table_id}')">
                <i class="fas {icon}"></i> {label}</button>''')
        html += (f'''<button class="action-btn btn-light"
            onclick="document.getElementById('search-{table_id}').value='';'''
                 f'''searchTable('{table_id}','search-{table_id}')">
            <i class="fas fa-sync"></i> Clear</button></div>''')
        return html

    def _gene_family_info(self, border_color: str, title: str, items: list) -> str:
        """Explanatory box listing the biological role of each gene family.

        items: list of (name, description) tuples.
        """
        html = (f'<div style="margin:10px 0 20px 0;background:#f8f9fa;padding:15px;'
                f'border-radius:8px;font-size:.9em;border-left:4px solid {border_color};">'
                f'<strong><i class="fas fa-info-circle"></i> {title}</strong><br>')
        for name, desc in items:
            html += f'• <strong>{name}</strong> – {desc}<br>'
        html += '</div>'
        return html

    def _database_cards(self, gene_dict: dict, db_labels: dict = None) -> str:
        """Per-database summary cards showing gene count and top hits."""
        db_labels = db_labels or {}
        html = ('<h3 style="margin-top:30px;"><i class="fas fa-database"></i> '
                'Database Summary</h3>'
                '<div style="display:grid;grid-template-columns:'
                'repeat(auto-fit,minmax(300px,1fr));gap:20px;margin:20px 0;">')
        for db, genes in gene_dict.items():
            label = db_labels.get(db, db.upper() if db != 'amrfinder' else 'AMRfinder')
            top = ', '.join(f"{g['gene']} ({g['count']})" for g in genes[:3])
            total = sum(g['count'] for g in genes)
            html += f'''<div class="database-section">
                <h4>{label}</h4>
                <p><strong>{len(genes)} unique genes</strong>
                (Total occurrences: {total})</p>
                <p>Top genes: {top}</p></div>'''
        html += '</div>'
        return html

    def _grouping_dropdown(self, table_id: str) -> str:
        """Single-select dropdown that triggers client-side genome grouping."""
        return f'''
        <div class="grouping-controls">
            <strong><i class="fas fa-layer-group"></i> Group genomes by:</strong>
            <select class="group-select" onchange="groupGenomesByTyping('{table_id}', this.value)">
                <option value="">— None (flat list) —</option>
                <optgroup label="Single field">
                    <option value="MLST">MLST</option>
                    <option value="spa">spa</option>
                    <option value="SCCmec_CGE">SCCmec (CGE)</option>
                    <option value="SCCmec_RPet">SCCmec (RPet)</option>
                    <option value="SCCmec_Subtype">SCCmec Subtype</option>
                    <option value="Capsule">Capsule</option>
                    <option value="agr">agr</option>
                </optgroup>
                <optgroup label="Two fields">
                    <option value="MLST-spa">MLST + spa</option>
                    <option value="MLST-Subtype">MLST + SCCmec Subtype</option>
                    <option value="MLST-Capsule">MLST + Capsule</option>
                    <option value="MLST-agr">MLST + agr</option>
                    <option value="spa-Capsule">spa + Capsule</option>
                    <option value="spa-agr">spa + agr</option>
                </optgroup>
                <optgroup label="Three fields">
                    <option value="MLST-spa-SCCmec_CGE">MLST + spa + SCCmec (CGE)</option>
                    <option value="MLST-spa-Subtype">MLST + spa + SCCmec Subtype</option>
                    <option value="MLST-Capsule-agr">MLST + Capsule + agr</option>
                    <option value="spa-Capsule-agr">spa + Capsule + agr</option>
                </optgroup>
                <optgroup label="Four fields">
                    <option value="MLST-spa-SCCmec_CGE-agr">MLST + spa + SCCmec (CGE) + agr</option>
                    <option value="MLST-spa-Subtype-agr">MLST + spa + SCCmec Subtype + agr</option>
                </optgroup>
                <optgroup label="Five fields">
                    <option value="MLST-spa-SCCmec_CGE-agr-Capsule">MLST + spa + SCCmec (CGE) + agr + Capsule</option>
                </optgroup>
            </select>
            <button class="group-btn" onclick="resetGenomeList('{table_id}')"><i class="fas fa-undo"></i> Reset</button>
        </div>'''

    def _search_pair(self, table_id: str, placeholder: str) -> str:
        """Two search boxes: one to filter rows, one to highlight genome tags."""
        return f'''
        <input type="text" class="search-box" id="search-{table_id}"
               onkeyup="searchTable('{table_id}', 'search-{table_id}')"
               placeholder="🔍 {placeholder}">
        <input type="text" class="search-box" id="highlight-{table_id}"
               onkeyup="highlightGenome('{table_id}', 'highlight-{table_id}')"
               placeholder="🔍 Highlight genomes containing specific text">'''

    def _combo_table(self, table_id: str, title: str, combo_dict: dict,
                     sample_key: str = 'Samples') -> str:
        """A single 'X – Y Combination | Count | Samples' table."""
        if not combo_dict:
            return ''
        rows = ''
        for combo, samples in sorted(combo_dict.items(),
                                     key=lambda x: len(x[1]), reverse=True):
            tags = ''.join(f'<span class="genome-tag">{esc(s)}</span>' for s in samples)
            rows += (f'<tr><td><strong>{esc(combo)}</strong></td>'
                     f'<td>{len(samples)}</td>'
                     f'<td><div class="genome-list">{tags}</div></td></tr>')
        return f'''
        <h3>🔗 {title}</h3>
        {self._search_pair(table_id, f'Search {title}...')}
        <div class="master-scrollable-container">
            <table id="{table_id}" class="data-table">
                <thead><tr>
                    <th data-sort="string">{title}</th>
                    <th data-sort="number">Count</th>
                    <th data-sort="string">{sample_key}</th>
                </tr></thead>
                <tbody>{rows}</tbody>
            </table>
        </div>'''

    def _gene_table(self, table_id: str, title: str, gene_list: list,
                    total_samples: int, critical_set: set = None) -> str:
        """Generic 'Gene | Database | Count | % | Genomes' table."""
        critical_set = critical_set or set()
        rows = ''
        for g in gene_list:
            is_crit = any(c in g['gene'].lower() for c in critical_set)
            gene_display = f"<strong>{esc(g['gene'])}</strong>" + (" ⚠️" if is_crit else "")
            freq = g.get('frequency', str(g['count']))
            pct = '0%'
            if '(' in freq:
                pct = freq.split('(')[-1].replace(')', '').strip()
            elif g['count'] > 0 and total_samples > 0:
                pct = f"{(g['count'] / total_samples) * 100:.1f}%"
            tags = ''.join(f'<span class="genome-tag">{esc(x)}</span>' for x in g.get('genomes', []))
            rows += (f'<tr><td>{gene_display}</td><td>{esc(g["database"])}</td>'
                     f'<td><strong>{g["count"]}</strong></td><td>{pct}</td>'
                     f'<td><div class="genome-list">{tags}</div></td></tr>')
        return f'''
        <h3><i class="fas fa-list"></i> {title}</h3>
        {self._grouping_dropdown(table_id)}
        {self._search_pair(table_id, f'Search {title}...')}
        <div class="master-scrollable-container">
            <table id="{table_id}" class="data-table">
                <thead><tr>
                    <th data-sort="string">Gene</th>
                    <th data-sort="string">Database</th>
                    <th data-sort="number">Count</th>
                    <th data-sort="number">Percentage</th>
                    <th data-sort="string">Genomes (scrollable, groupable)</th>
                </tr></thead>
                <tbody>{rows}</tbody>
            </table>
        </div>'''

    def _distribution_table(self, table_id: str, dist: Counter,
                            label: str, extra_headers: list = None) -> str:
        """Simple 'Value | Count | %' table for categorical distributions."""
        extra_headers = extra_headers or []
        total = sum(dist.values())
        rows = ''
        for val, cnt in dist.most_common():
            if val in ('Not Assigned', 'ND', ''):
                continue
            pct = (cnt / total * 100) if total else 0
            extras = ''.join('<td>—</td>' for _ in extra_headers)
            rows += f'<tr><td><strong>{esc(val)}</strong></td><td>{cnt}</td><td>{pct:.1f}%</td>{extras}</tr>'
        header_extra = ''.join(f'<th data-sort="string">{h}</th>' for h in extra_headers)
        return f'''
        <div class="scrollable-table">
            <table id="{table_id}" class="data-table">
                <thead><tr>
                    <th data-sort="string">{label}</th>
                    <th data-sort="number">Count</th>
                    <th data-sort="number">Percentage</th>
                    {header_extra}
                </tr></thead>
                <tbody>{rows}</tbody>
            </table>
        </div>'''

    def _colorize_capsule_cell(self, value: str) -> str:
        """Render a capsule type value with a colour badge (Type 5 = green, Type 8 = red)."""
        v = str(value).strip()
        if v == 'Type 5':
            return ('<span style="background:#d4edda;color:#155724;'
                    'font-weight:bold;padding:3px 10px;border-radius:10px;'
                    'display:inline-block;">Type 5</span>')
        if v == 'Type 8':
            return ('<span style="background:#f8d7da;color:#721c24;'
                    'font-weight:bold;padding:3px 10px;border-radius:10px;'
                    'display:inline-block;">Type 8</span>')
        return esc(v)

    # -------------------------------------------------------------------------
    # CSS + JS
    # -------------------------------------------------------------------------
    def _get_css(self) -> str:
        """Global stylesheet embedded in the report head."""
        return """
<style>
:root{
  --summary-color:#4CAF50;--sample_overview-color:#2196F3;--qc-color:#607D8B;
  --mlst-color:#FF9800;--spa-color:#9C27B0;--sccmec-color:#009688;--capsule-color:#00ACC1;
  --mrsa-color:#795548;--agr-color:#8B5CF6;--amr-color:#F44336;
  --virulence-color:#E91E63;--bacmet-color:#FF5722;--plasmids-color:#673AB7;
  --mutation-color:#00BCD4;--mge-color:#16A085;--patterns-color:#3F51B5;
  --aiguide-color:#00BCD4;--citation-color:#8BC34A;--funding-color:#FFC107;
  --export-color:#9E9E9E;--calltoaction-color:#F472B6;
}
*{margin:0;padding:0;box-sizing:border-box}
body{font-family:'Segoe UI',Tahoma,Geneva,Verdana,sans-serif;line-height:1.6;color:#333;background:#f5f5f5;min-width:1200px}
.container{max-width:none;margin:0 auto;padding:20px;width:100%;overflow-x:auto}
.main-header{background:linear-gradient(135deg,#006400 0%,#228B22 100%);color:white;padding:30px;border-radius:15px;box-shadow:0 10px 30px rgba(0,0,0,.2);margin-bottom:30px;text-align:center}
.main-header h1{font-size:2.8em;margin-bottom:10px;color:white}
.metadata-bar{background:rgba(255,255,255,.1);padding:15px;border-radius:10px;margin:20px 0;display:flex;justify-content:space-around;flex-wrap:wrap;gap:15px;backdrop-filter:blur(10px)}
.metadata-item{display:flex;align-items:center;gap:8px;font-size:.95em}
.dashboard-grid{display:grid;grid-template-columns:repeat(auto-fit,minmax(200px,1fr));gap:20px;margin-bottom:30px}
.dashboard-card{background:white;padding:20px;border-radius:12px;box-shadow:0 5px 20px rgba(0,0,0,.1);text-align:center;transition:all .3s ease;cursor:pointer;border-left:5px solid;position:relative;overflow:hidden}
.dashboard-card:hover{transform:translateY(-10px);box-shadow:0 15px 30px rgba(0,0,0,.2)}
.dashboard-card .card-number{font-size:2.5em;font-weight:bold;margin:10px 0;background:linear-gradient(90deg,#006400,#228B22);-webkit-background-clip:text;-webkit-text-fill-color:transparent}
.dashboard-card .card-label{font-size:.9em;color:#555;font-weight:600}
.card-summary{border-left-color:var(--summary-color)}
.card-mlst{border-left-color:var(--mlst-color)}
.card-spa{border-left-color:var(--spa-color)}
.card-sccmec{border-left-color:var(--sccmec-color)}
.card-capsule{border-left-color:var(--capsule-color)}
.card-mrsa{border-left-color:var(--mrsa-color)}
.card-agr{border-left-color:var(--agr-color)}
.card-amr{border-left-color:var(--amr-color)}
.card-virulence{border-left-color:var(--virulence-color)}
.card-mge{border-left-color:var(--mge-color)}
.card-patterns{border-left-color:var(--patterns-color)}
.tab-navigation{display:flex;gap:5px;margin-bottom:20px;flex-wrap:wrap;background:white;padding:15px;border-radius:12px;box-shadow:0 5px 20px rgba(0,0,0,.1);position:sticky;top:10px;z-index:100}
.tab-button{padding:10px 16px;background:#f5f5f5;border:none;border-radius:8px;cursor:pointer;font-weight:600;color:#666;transition:all .3s ease;display:flex;align-items:center;gap:6px;font-size:.88em}
.tab-button.active{color:white}
.tab-button.summary.active{background:var(--summary-color)}
.tab-button.sample_overview.active{background:var(--sample_overview-color)}
.tab-button.qc.active{background:var(--qc-color)}
.tab-button.mlst.active{background:var(--mlst-color)}
.tab-button.spa.active{background:var(--spa-color)}
.tab-button.sccmec.active{background:var(--sccmec-color)}
.tab-button.capsule.active{background:var(--capsule-color)}
.tab-button.mrsa.active{background:var(--mrsa-color)}
.tab-button.agr.active{background:var(--agr-color)}
.tab-button.amr.active{background:var(--amr-color)}
.tab-button.virulence.active{background:var(--virulence-color)}
.tab-button.bacmet.active{background:var(--bacmet-color)}
.tab-button.plasmids.active{background:var(--plasmids-color)}
.tab-button.mutation.active{background:var(--mutation-color)}
.tab-button.mge.active{background:var(--mge-color)}
.tab-button.patterns.active{background:var(--patterns-color)}
.tab-button.aiguide.active{background:var(--aiguide-color)}
.tab-button.calltoaction.active{background:var(--calltoaction-color)}
.tab-button.citation.active{background:var(--citation-color)}
.tab-button.funding.active{background:var(--funding-color)}
.tab-button.export.active{background:var(--export-color)}
.tab-content{display:none;background:white;padding:30px;border-radius:15px;box-shadow:0 10px 30px rgba(0,0,0,.1);margin-bottom:30px;animation:fadeIn .5s ease;width:100%;overflow-x:auto}
.tab-content.active{display:block}
@keyframes fadeIn{from{opacity:0;transform:translateY(20px)}to{opacity:1;transform:translateY(0)}}
.section-header{color:#2c3e50;margin-bottom:25px;padding-bottom:15px;border-bottom:3px solid;font-size:1.8em;display:flex;align-items:center;justify-content:space-between}
.summary-header{border-color:var(--summary-color)}
.sample_overview-header{border-color:var(--sample_overview-color)}
.qc-header{border-color:var(--qc-color)}
.mlst-header{border-color:var(--mlst-color)}
.spa-header{border-color:var(--spa-color)}
.sccmec-header{border-color:var(--sccmec-color)}
.capsule-header{border-color:var(--capsule-color)}
.mrsa-header{border-color:var(--mrsa-color)}
.agr-header{border-color:var(--agr-color)}
.amr-header{border-color:var(--amr-color)}
.virulence-header{border-color:var(--virulence-color)}
.bacmet-header{border-color:var(--bacmet-color)}
.plasmids-header{border-color:var(--plasmids-color)}
.mutation-header{border-color:var(--mutation-color)}
.mge-header{border-color:var(--mge-color)}
.patterns-header{border-color:var(--patterns-color)}
.aiguide-header{border-color:var(--aiguide-color)}
.calltoaction-header{border-color:var(--calltoaction-color)}
.citation-header{border-color:var(--citation-color)}
.funding-header{border-color:var(--funding-color)}
.export-header{border-color:var(--export-color)}
.data-table{width:100%;border-collapse:collapse;margin:20px 0;font-size:.95em;box-shadow:0 2px 10px rgba(0,0,0,.1);border-radius:8px;overflow:hidden}
.data-table th{background:#2c3e50;color:white;padding:14px;text-align:left;font-weight:600;position:sticky;top:0;white-space:nowrap;cursor:pointer}
.data-table th:hover{background:#1a252f}
.data-table td{padding:12px;border-bottom:1px solid #e0e0e0;vertical-align:top}
.data-table tr:hover td{background:#f8f9fa}
.scrollable-table{max-height:600px;overflow-y:auto;border:1px solid #e0e0e0;border-radius:8px;margin:20px 0;width:100%}
.master-scrollable-container{width:100%;overflow-x:auto;border:1px solid #e0e0e0;border-radius:8px;margin:20px 0}
.gene-centric-table{table-layout:fixed;width:auto}
.gene-centric-table th:first-child,
.gene-centric-table td:first-child{width:200px;box-sizing:border-box;word-break:break-word}
.gene-centric-table th:nth-child(2),
.gene-centric-table td:nth-child(2){width:120px;box-sizing:border-box;word-break:break-word}
.gene-centric-table th:nth-child(3),
.gene-centric-table td:nth-child(3){width:80px;box-sizing:border-box;text-align:center}
.gene-centric-table th:nth-child(4),
.gene-centric-table td:nth-child(4){width:110px;box-sizing:border-box;text-align:center}
.gene-centric-table th:last-child,
.gene-centric-table td:last-child{width:700px}
.data-table:not(.gene-centric-table) th{white-space:normal;word-break:break-word}
.genome-list{display:flex;flex-wrap:wrap;gap:6px;padding:8px;background:#f8f9fa;border-radius:5px;width:100%;max-height:280px;overflow:auto;box-sizing:border-box}
.genome-group{margin-bottom:10px;width:100%;box-sizing:border-box}
.genome-group-header{font-weight:bold;background:#e0e0e0;padding:5px 10px;border-radius:4px;margin:6px 0;font-size:.9em;display:inline-block;max-width:100%;box-sizing:border-box;word-break:break-word}
.genome-group-tags{display:flex;flex-wrap:wrap;gap:6px;margin-left:10px;max-width:calc(100% - 10px);box-sizing:border-box}
.genome-tag{display:inline-block;background:#e6ffe6;color:#006400;padding:4px 12px;border-radius:12px;font-size:.9em;border:1px solid #b3ffb3;white-space:nowrap;margin:2px;max-width:100%;box-sizing:border-box}
.genome-tag.highlight{background:#ffff99 !important;color:#000 !important;border:1px solid #ffc107}
.search-box{width:100%;padding:12px;margin-bottom:12px;border:2px solid #e0e0e0;border-radius:8px;font-size:1em;transition:all .3s ease}
.search-box:focus{outline:none;border-color:#006400;box-shadow:0 0 0 3px rgba(0,100,0,.1)}
.badge{display:inline-block;padding:5px 15px;border-radius:20px;font-size:.85em;font-weight:600;margin:2px}
.badge-mrsa{background:#8B0000;color:white}
.badge-mssa{background:#4682B4;color:white}
.badge-critical{background:#DC143C;color:white}
.alert-box{padding:20px;border-radius:10px;margin:20px 0;display:flex;align-items:flex-start;gap:20px;border-left:5px solid}
.alert-success{background:#d4edda;color:#155724;border-left-color:#28a745}
.alert-warning{background:#fff3cd;color:#856404;border-left-color:#ffc107}
.alert-danger{background:#f8d7da;color:#721c24;border-left-color:#dc3545}
.alert-info{background:#d1ecf1;color:#0c5460;border-left-color:#17a2b8}
.scientific-note{border-radius:8px}
.action-buttons{display:flex;gap:10px;margin:20px 0;flex-wrap:wrap}
.action-btn{padding:10px 20px;border:none;border-radius:8px;cursor:pointer;font-weight:600;display:flex;align-items:center;gap:8px;transition:all .3s ease;text-decoration:none}
.action-btn:hover{transform:translateY(-2px);box-shadow:0 5px 15px rgba(0,0,0,.2)}
.btn-primary{background:#006400;color:white}
.btn-success{background:#28a745;color:white}
.btn-danger{background:#dc3545;color:white}
.btn-warning{background:#ffc107;color:black}
.btn-info{background:#17a2b8;color:white}
.btn-secondary{background:#6c757d;color:white}
.btn-light{background:#f8f9fa;color:#212529;border:1px solid #dee2e6}
.database-section{margin:30px 0;padding:25px;border-radius:12px;background:#f8f9fa;box-shadow:0 3px 15px rgba(0,0,0,.08)}
.print-section-btn{background:#006400;color:white;border:none;border-radius:5px;padding:8px 15px;cursor:pointer;display:flex;align-items:center;gap:5px;font-size:.9em}
.print-section-btn:hover{background:#228B22}
.footer{text-align:center;padding:30px;color:white;margin-top:40px;border-radius:15px;background:linear-gradient(135deg,#2c3e50 0%,#34495e 100%)}
.mrsa-highlight{background:#ffe6e6 !important;border-left:3px solid #8B0000 !important}
.sort-icon{margin-left:5px;font-size:.8em;opacity:.6}
.grouping-controls{background:#f0f7f0;padding:12px;border-radius:8px;margin:15px 0;display:flex;flex-wrap:wrap;gap:10px;align-items:center;border-left:4px solid #006400}
.grouping-controls label,.grouping-controls strong{font-weight:bold;margin-right:5px;color:#2c3e50}
.group-select{padding:8px 14px;border:2px solid #006400;border-radius:6px;background:white;color:#333;font-size:.9em;font-weight:600;cursor:pointer;min-width:240px}
.group-select:focus{outline:none;box-shadow:0 0 0 3px rgba(0,100,0,.15)}
.group-btn{background:white;border:1px solid #006400;color:#006400;padding:6px 12px;border-radius:20px;cursor:pointer;font-size:.85em;transition:all .2s}
.group-btn:hover{background:#006400;color:white}
#qc-table{width:max-content;min-width:100%}
.stat-card{color:white;padding:18px;border-radius:10px;text-align:center;box-shadow:0 4px 15px rgba(0,0,0,.15);transition:transform .2s}
.stat-card:hover{transform:translateY(-3px)}
.stat-card .stat-value{font-size:1.9em;font-weight:bold;margin-bottom:4px}
.stat-card .stat-label{font-size:.85em;opacity:.95;text-transform:uppercase;letter-spacing:.5px}
.stats-grid{display:grid;grid-template-columns:repeat(auto-fit,minmax(180px,1fr));gap:15px;margin:20px 0}
.typing-badge{display:inline-block;padding:3px 10px;border-radius:12px;font-size:.8em;font-weight:600;background:#e0e0e0;color:#333;border:1px solid #ccc}
.typing-badge.agr-I{background:#16a34a;color:white;border-color:#16a34a}
.typing-badge.agr-II{background:#2563eb;color:white;border-color:#2563eb}
.typing-badge.agr-III{background:#f59e0b;color:white;border-color:#f59e0b}
.typing-badge.agr-IV{background:#dc2626;color:white;border-color:#dc2626}
.typing-badge.agr-NA{background:#6b7280;color:white;border-color:#6b7280}
.vir-details{margin-top:8px}
.vir-details summary{cursor:pointer;color:#006400;font-weight:600;font-size:.85em;list-style:none}
.vir-details summary::-webkit-details-marker{display:none}
.vir-details summary::before{content:'▸ ';transition:transform .2s}
.vir-details[open] summary::before{content:'▾ '}
.vir-details .vir-gene-tags{margin-top:6px;display:flex;flex-wrap:wrap;gap:4px}
.accordion{margin:20px 0}
.accordion-item{background:#f8f9fa;border:1px solid #dee2e6;margin-bottom:10px;border-radius:8px;overflow:hidden}
.accordion-header{background:#e9ecef;padding:12px 20px;cursor:pointer;font-weight:bold;color:#1e3a8a;display:flex;justify-content:space-between;align-items:center}
.accordion-header:hover{background:#dee2e6}
.accordion-content{padding:15px 20px;border-top:1px solid #dee2e6;background:white}
.citation-list{list-style:none;padding-left:0}
.citation-list li{margin-bottom:15px;padding-bottom:10px;border-bottom:1px solid #e0e0e0}
.citation-list li:last-child{border-bottom:none}
.copy-btn{background:#006400;color:white;border:none;padding:4px 12px;border-radius:20px;cursor:pointer;font-size:.8em;margin-left:10px}
.copy-btn:hover{background:#228B22}
@media print{body *{visibility:hidden}.tab-content.active,.tab-content.active *{visibility:visible}.tab-content.active{position:absolute;left:0;top:0;width:100%;padding:20px;box-shadow:none;border-radius:0}.print-section-btn,.tab-navigation,.dashboard-grid,.search-box,.action-buttons,.grouping-controls{display:none !important}}
/* Compare tab */
.compare-header { border-color: #0d9488; }
/* Mode toggle */
.mode-toggle { display: flex; gap: 8px; margin-bottom: 18px; }
.mode-btn { padding: 10px 22px; border: 2px solid #cbd5e1; background: white; color: #475569;
            border-radius: 8px; font-weight: 600; cursor: pointer; font-size: 0.95em;
            display: flex; align-items: center; gap: 8px; transition: all 0.2s; }
.mode-btn:hover { border-color: #006400; color: #006400; }
.mode-btn.active { background: #006400; border-color: #006400; color: white; }
.compare-mode { display: none; }
.compare-mode.active { display: block; animation: fadeIn 0.3s ease; }

/* Visual panel: gauge + metrics */
.compare-visuals-grid { display: grid; grid-template-columns: 180px 1fr; gap: 16px;
                        margin: 20px 0; align-items: center; }
.visual-card { background: white; border-radius: 12px; padding: 20px;
               box-shadow: 0 2px 8px rgba(0,0,0,0.06); border: 1px solid #e2e8f0; }
.gauge-card { display: flex; justify-content: center; align-items: center; padding: 12px; }
.similarity-gauge { width: 160px; height: 160px; display: block; }
.metrics-card { display: flex; flex-direction: column; gap: 18px; }
.metric-block { display: flex; flex-direction: column; gap: 6px; }
.metric-label { font-size: 0.85em; font-weight: 700; color: #475569;
                text-transform: uppercase; letter-spacing: 0.5px; }
.metric-bar-track { background: #f1f5f9; height: 22px; border-radius: 11px;
                    overflow: hidden; position: relative; }
.metric-bar-fill { height: 100%; border-radius: 11px; transition: width 0.8s ease; min-width: 4px; }
.metric-value { font-size: 0.88em; color: #334155; font-weight: 600; }

/* Typing pills */
.typing-pill { display: inline-block; padding: 5px 12px; border-radius: 12px;
               font-size: 0.85em; font-weight: 600; border: 1.5px solid;
               font-family: 'Consolas', monospace; }

/* Gene proportion bars */
.gene-proportion-bar { display: flex; height: 12px; border-radius: 6px;
                       overflow: hidden; margin-bottom: 14px; background: #f1f5f9; }
.gene-prop-segment { height: 100%; transition: width 0.6s ease; min-width: 2px; }
.gene-prop-segment.shared { background: #10b981; }
.gene-prop-segment.only-a { background: #3b82f6; }
.gene-prop-segment.only-b { background: #ec4899; }

/* Cluster cards */
.cluster-card { background: white; border-radius: 10px; padding: 16px 20px;
                margin-bottom: 12px; border-left: 5px solid;
                box-shadow: 0 2px 8px rgba(0,0,0,0.06); transition: transform 0.15s; }
.cluster-card:hover { transform: translateX(4px); }
.cluster-header { display: flex; align-items: center; gap: 14px; margin-bottom: 12px; flex-wrap: wrap; }
.cluster-badge { color: white; padding: 4px 14px; border-radius: 14px;
                 font-size: 0.85em; font-weight: 700; }
.cluster-stats { font-size: 0.9em; color: #64748b; font-weight: 600; }
.cluster-samples { display: flex; flex-wrap: wrap; gap: 6px; }
.cluster-sample-pill { background: #f1f5f9; color: #1e293b; padding: 4px 12px;
                       border-radius: 12px; font-size: 0.85em; font-weight: 600;
                       font-family: 'Consolas', monospace; border: 1px solid #cbd5e1; }

/* Cluster controls */
.cluster-controls { display: flex; gap: 14px; align-items: flex-end; flex-wrap: wrap; margin-bottom: 20px; }
.cluster-controls label { display: flex; flex-direction: column; gap: 4px; }
.cluster-controls label span { font-size: 0.85em; font-weight: 600; color: #475569; }
.cluster-controls select { padding: 10px 14px; border: 2px solid #cbd5e1;
                           border-radius: 8px; background: white; font-size: 0.95em; min-width: 220px; }

/* Heatmap */
.heatmap-wrapper { background: white; border-radius: 12px; padding: 20px;
                   margin-top: 24px; box-shadow: 0 2px 8px rgba(0,0,0,0.06); }
.heatmap-wrapper h3 { margin-bottom: 8px; color: #1e293b; }
.heatmap-grid { display: grid; gap: 2px; margin-top: 16px;
                overflow-x: auto; padding: 4px; }
.heatmap-cell { width: 32px; height: 32px; display: flex;
                align-items: center; justify-content: center;
                font-size: 0.7em; font-weight: 700; border-radius: 3px;
                transition: transform 0.15s; cursor: default; }
.heatmap-cell:hover { transform: scale(1.15); z-index: 2; box-shadow: 0 2px 8px rgba(0,0,0,0.2); }
.heatmap-col-header, .heatmap-row-header { background: transparent; color: #64748b;
                                           font-size: 0.75em; }
.heatmap-row-header { justify-content: flex-start; padding-left: 6px;
                      font-family: 'Consolas', monospace; }
.heatmap-col-header { writing-mode: vertical-rl; text-orientation: mixed;
                      height: 70px; width: 32px; padding: 4px 0;
                      font-family: 'Consolas', monospace; }
.heatmap-header-label { font-size: 0.85em; }
.heatmap-corner { background: transparent; }
.compare-panel { background: #f8f9fa; border-radius: 12px; padding: 20px; margin-bottom: 24px; }
.compare-select-row { display: flex; gap: 12px; align-items: flex-end; flex-wrap: wrap; margin-bottom: 16px; }
.compare-select-row label { display: flex; flex-direction: column; gap: 4px; }
.compare-select-row label span { font-size: 0.85em; font-weight: 600; color: #475569; }
.compare-select-row select { padding: 10px 14px; border: 2px solid #cbd5e1; border-radius: 8px; background: white; font-size: 0.95em; min-width: 220px; }
.compare-select-row select:focus { outline: none; border-color: #006400; box-shadow: 0 0 0 3px rgba(0,100,0,0.15); }

.compare-verdict-box { padding: 16px 20px; border-radius: 10px; margin-bottom: 16px; display: flex; gap: 12px; align-items: flex-start; font-size: 0.95em; }
.compare-verdict-box.danger  { background: #fef2f2; border-left: 5px solid #dc2626; color: #7f1d1d; }
.compare-verdict-box.warning { background: #fffbeb; border-left: 5px solid #f59e0b; color: #78350f; }
.compare-verdict-box.info    { background: #eff6ff; border-left: 5px solid #2563eb; color: #1e3a8a; }

.compare-warning { padding: 20px; background: #fef2f2; color: #7f1d1d; border-radius: 8px; border-left: 4px solid #dc2626; }

.compare-section { background: white; border-radius: 12px; padding: 20px; margin-bottom: 20px; box-shadow: 0 2px 8px rgba(0,0,0,0.06); }
.compare-section h3 { margin-bottom: 14px; color: #1e293b; font-size: 1.15em; padding-bottom: 8px; border-bottom: 2px solid #e2e8f0; }

.compare-table { display: flex; flex-direction: column; gap: 2px; }
.compare-row { display: grid; grid-template-columns: 180px 1fr 1fr; gap: 12px; padding: 10px 14px; border-radius: 6px; align-items: center; font-size: 0.92em; }
.compare-row.header { background: #f1f5f9; font-weight: 700; border-bottom: 2px solid #cbd5e1; }
.compare-row.diff { background: #fef2f2; border-left: 3px solid #dc2626; }
.compare-row.match { background: #f0fdf4; border-left: 3px solid #16a34a; }
.compare-row .row-label { font-weight: 600; color: #475569; font-size: 0.88em; text-transform: uppercase; letter-spacing: 0.3px; }
.compare-row .row-val { color: #1e293b; word-break: break-word; }

.compare-gene-block { margin-bottom: 24px; padding-bottom: 20px; border-bottom: 1px dashed #e2e8f0; }
.compare-gene-block:last-child { border-bottom: none; margin-bottom: 0; padding-bottom: 0; }
.compare-gene-block h4 { margin-bottom: 12px; color: #1e293b; font-size: 1.02em; display: flex; align-items: center; gap: 10px; flex-wrap: wrap; }
.gene-count { font-size: 0.82em; font-weight: 400; color: #64748b; background: #f1f5f9; padding: 3px 10px; border-radius: 12px; }

.gene-columns { display: grid; grid-template-columns: repeat(auto-fit, minmax(220px, 1fr)); gap: 14px; }
.gene-col { background: #f8fafc; border-radius: 8px; padding: 12px; border: 1px solid #e2e8f0; }
.gene-col-head { font-size: 0.85em; font-weight: 700; color: #475569; margin-bottom: 8px; text-transform: uppercase; letter-spacing: 0.3px; }

.compare-gene { display: inline-block; font-size: 0.82em; padding: 3px 10px; border-radius: 12px; margin: 2px; font-family: 'Consolas', monospace; }
.compare-gene.shared { background: #dcfce7; color: #14532d; border: 1px solid #86efac; }
.compare-gene.only-a { background: #dbeafe; color: #1e3a8a; border: 1px solid #93c5fd; }
.compare-gene.only-b { background: #fce7f3; color: #831843; border: 1px solid #f9a8d4; }

/* Compare — full profile expanders */
.compare-full-details { background: #f8fafc; border: 1px solid #e2e8f0; border-radius: 10px; padding: 12px 16px; margin-bottom: 10px; }
.compare-full-details[open] { background: #fff; border-color: #006400; }
.compare-full-details summary { cursor: pointer; font-weight: 700; color: #006400; list-style: none; display: flex; justify-content: space-between; align-items: center; gap: 10px; flex-wrap: wrap; padding: 4px 0; }
.compare-full-details summary::-webkit-details-marker { display: none; }
.compare-full-details summary::before { content: '▸'; display: inline-block; margin-right: 8px; color: #006400; }
.compare-full-details[open] summary::before { content: '▾'; }
.details-cat-name { font-size: 1em; }
.details-counts { font-size: 0.85em; font-weight: 500; color: #64748b; background: #f1f5f9; padding: 3px 10px; border-radius: 12px; margin-left: auto; }

.compare-full-profile { display: grid; grid-template-columns: 1fr 1fr; gap: 16px; margin-top: 14px; padding-top: 14px; border-top: 1px dashed #e2e8f0; }
.profile-col { background: #f8fafc; border-radius: 10px; padding: 14px; border: 1px solid #e2e8f0; }
.profile-head { display: flex; justify-content: space-between; align-items: center; gap: 8px; margin-bottom: 10px; padding-bottom: 8px; border-bottom: 2px solid #e2e8f0; flex-wrap: wrap; }
.profile-sample-label { font-weight: 700; color: #1e293b; font-size: 0.92em; word-break: break-word; }
.profile-sample-count { font-size: 0.78em; font-weight: 600; color: #006400; background: #dcfce7; padding: 2px 10px; border-radius: 12px; white-space: nowrap; }
.empty-note { font-size: 0.85em; color: #94a3b8; font-style: italic; padding: 6px 0; display: inline-block; }

@media (max-width: 900px) {
    .compare-full-profile { grid-template-columns: 1fr; }
    .compare-row { grid-template-columns: 1fr; }
}
    .compare-visuals-grid { grid-template-columns: 1fr; }
    .heatmap-grid { grid-template-columns: 120px repeat(var(--n, 10), 28px) !important; }
</style>"""

    def _get_js(self, typing_json: str) -> str:
        """Client-side scripting for tabs, search, sort, grouping, CSV export."""
        return f"""
<script>
var sampleTyping = {typing_json};
var originalGenomeLists = {{}};
function escapeHtml(str) {{
    return String(str)
        .replace(/&/g, '&amp;')
        .replace(/</g, '&lt;')
        .replace(/>/g, '&gt;')
        .replace(/"/g, '&quot;');
}}

function switchTab(tabName){{
    document.querySelectorAll('.tab-content').forEach(t=>t.classList.remove('active'));
    document.querySelectorAll('.tab-button').forEach(b=>b.classList.remove('active'));
    var content = document.getElementById(tabName+'-tab');
    var button = document.querySelector('.tab-button.'+tabName);
    if(content) content.classList.add('active');
    if(button) button.classList.add('active');
    if(event && event.currentTarget) event.currentTarget.classList.add('active');
    window.location.hash = tabName;
}}
function searchTable(tableId, searchId){{
    var input = document.getElementById(searchId);
    if(!input) return;
    var filter = input.value.toUpperCase();
    var table = document.getElementById(tableId);
    if(!table) return;
    var rows = table.tBodies[0].rows;
    for(var i=0;i<rows.length;i++){{
        var cells = rows[i].getElementsByTagName('td');
        var found = false;
        for(var j=0;j<cells.length;j++){{
            if(cells[j] && (cells[j].textContent||cells[j].innerText).toUpperCase().indexOf(filter)>-1){{
                found = true; break;
            }}
        }}
        rows[i].style.display = found ? '' : 'none';
    }}
}}
function highlightGenome(tableId, searchId){{
    var el = document.getElementById(searchId);
    if(!el) return;
    var filter = el.value.toUpperCase().trim();
    var table = document.getElementById(tableId);
    if(!table) return;
    table.querySelectorAll('.genome-tag').forEach(function(t){{
        t.classList.remove('highlight');
        if(filter && t.textContent.toUpperCase().indexOf(filter)>-1) t.classList.add('highlight');
    }});
}}
function getTypingValue(genome, groupBy){{
    var info = sampleTyping[genome];
    if(!info) return "Unknown";
    function pick(k){{ return info[k] || 'Not Assigned'; }}
    switch(groupBy){{
        case 'MLST': return pick('MLST');
        case 'spa': return pick('spa');
        case 'SCCmec_CGE': return pick('SCCmec_CGE');
        case 'SCCmec_RPet': return pick('SCCmec_RPet');
        case 'SCCmec_Subtype': return pick('SCCmec_Subtype');
        case 'Capsule': return pick('Capsule');
        case 'agr': return pick('agr');
        case 'MLST-spa': return pick('MLST')+' - '+pick('spa');
        case 'MLST-Subtype': return pick('MLST')+' - '+pick('SCCmec_Subtype');
        case 'MLST-Capsule': return pick('MLST')+' - '+pick('Capsule');
        case 'MLST-agr': return pick('MLST')+' - '+pick('agr');
        case 'spa-Capsule': return pick('spa')+' - '+pick('Capsule');
        case 'spa-agr': return pick('spa')+' - '+pick('agr');
        case 'MLST-spa-SCCmec_CGE': return pick('MLST')+' - '+pick('spa')+' - '+pick('SCCmec_CGE');
        case 'MLST-spa-Subtype': return pick('MLST')+' - '+pick('spa')+' - '+pick('SCCmec_Subtype');
        case 'MLST-Capsule-agr': return pick('MLST')+' - '+pick('Capsule')+' - '+pick('agr');
        case 'spa-Capsule-agr': return pick('spa')+' - '+pick('Capsule')+' - '+pick('agr');
        case 'MLST-spa-SCCmec_CGE-agr': return pick('MLST')+' - '+pick('spa')+' - '+pick('SCCmec_CGE')+' - '+pick('agr');
        case 'MLST-spa-Subtype-agr': return pick('MLST')+' - '+pick('spa')+' - '+pick('SCCmec_Subtype')+' - '+pick('agr');
        case 'MLST-spa-SCCmec_CGE-agr-Capsule': return pick('MLST')+' - '+pick('spa')+' - '+pick('SCCmec_CGE')+' - '+pick('agr')+' - '+pick('Capsule');
    }}
    return "Unknown";
}}
function groupRowGenomes(row, groupBy, originalList){{
    var genomesCell = null;
    for(var i=0;i<row.cells.length;i++){{
        if(row.cells[i].querySelector('.genome-list')){{ genomesCell = row.cells[i]; break; }}
    }}
    if(!genomesCell) return;
    var genomes = originalList.slice();
    if(genomes.length===0){{ genomesCell.innerHTML='<div class="genome-list">None</div>'; return; }}
    var groups = {{}};
    genomes.forEach(function(g){{
        var key = getTypingValue(g, groupBy);
        if(!groups[key]) groups[key]=[];
        groups[key].push(g);
    }});
    var html = '<div class="genome-list">';
    for(var key in groups){{
        var tags = groups[key].map(g=>'<span class="genome-tag">'+g+'</span>').join('');
        html += '<div class="genome-group"><div class="genome-group-header">'+key+'</div><div class="genome-group-tags">'+tags+'</div></div>';
    }}
    html += '</div>';
    genomesCell.innerHTML = html;
}}
function groupGenomesByTyping(tableId, groupBy){{
    var table = document.getElementById(tableId);
    if(!table) return;
    var tbody = table.tBodies[0];
    if(!tbody) return;
    var rows = tbody.rows;
    for(var i=0;i<rows.length;i++){{
        var row = rows[i];
        var nameCell = row.cells[0];
        if(!nameCell) continue;
        var name = nameCell.textContent.trim().replace(/⚠️/g,'').trim();
        if(!originalGenomeLists[name]){{
            var genomesCell = null;
            for(var j=0;j<row.cells.length;j++){{
                if(row.cells[j].querySelector('.genome-list')){{ genomesCell = row.cells[j]; break; }}
            }}
            if(genomesCell){{
                var tags = genomesCell.querySelectorAll('.genome-tag');
                originalGenomeLists[name] = Array.from(tags).map(t=>t.textContent.trim());
            }} else {{
                originalGenomeLists[name] = [];
            }}
        }}
    }}
    if(!groupBy){{ resetGenomeList(tableId); return; }}
    for(var i=0;i<rows.length;i++){{
        var row = rows[i];
        var nameCell = row.cells[0];
        if(!nameCell) continue;
        var name = nameCell.textContent.trim().replace(/⚠️/g,'').trim();
        groupRowGenomes(row, groupBy, originalGenomeLists[name]||[]);
    }}
}}
function resetGenomeList(tableId){{
    var table = document.getElementById(tableId);
    if(!table) return;
    var tbody = table.tBodies[0];
    if(!tbody) return;
    var rows = tbody.rows;
    for(var i=0;i<rows.length;i++){{
        var row = rows[i];
        var nameCell = row.cells[0];
        if(!nameCell) continue;
        var name = nameCell.textContent.trim().replace(/⚠️/g,'').trim();
        var original = originalGenomeLists[name] || [];
        for(var j=0;j<row.cells.length;j++){{
            if(row.cells[j].querySelector('.genome-list')){{
                var tags = original.map(g=>'<span class="genome-tag">'+g+'</span>').join('');
                row.cells[j].innerHTML = '<div class="genome-list">'+tags+'</div>';
                break;
            }}
        }}
    }}
    var sel = document.querySelector('#'+tableId)?.closest('.tab-content')?.querySelector('.group-select');
    if(sel) sel.value = '';
}}
function sortTable(tableId, colIndex, type){{
    var table = document.getElementById(tableId);
    if(!table) return;
    var tbody = table.tBodies[0];
    var rows = Array.from(tbody.rows);
    var asc = table.getAttribute('data-sort-dir') !== 'asc';
    rows.sort(function(a,b){{
        var av = a.cells[colIndex].innerText.trim();
        var bv = b.cells[colIndex].innerText.trim();
        if(type==='number'){{
            av = parseFloat(av.replace(/,/g,''))||0;
            bv = parseFloat(bv.replace(/,/g,''))||0;
            return asc ? av-bv : bv-av;
        }}
        return asc ? av.localeCompare(bv) : bv.localeCompare(av);
    }});
    tbody.append.apply(tbody, rows);
    table.setAttribute('data-sort-dir', asc ? 'asc' : 'desc');
}}
function printSection(id){{
    var content = document.getElementById(id);
    if(!content) return;
    var w = window.open('', '_blank');
    var style = document.querySelector('style');
    w.document.write('<html><head><title>Print</title>');
    if(style) w.document.write('<style>'+style.textContent+'</style>');
    w.document.write('</head><body>'+content.innerHTML+'</body></html>');
    w.document.close(); w.print();
}}
function exportTableToCSV(tableId, filename){{
    var table = document.getElementById(tableId);
    if(!table) return;
    var rows = table.querySelectorAll('tr');
    var csv = [];
    for(var i=0;i<rows.length;i++){{
        var row=[], cols=rows[i].querySelectorAll('td, th');
        for(var j=0;j<cols.length;j++){{
            row.push('"'+(cols[j].innerText||'').replace(/"/g,'""')+'"');
        }}
        csv.push(row.join(','));
    }}
    var blob = new Blob([csv.join('\\n')],{{type:'text/csv'}});
    var a = document.createElement('a');
    a.download = filename; a.href = URL.createObjectURL(blob);
    document.body.appendChild(a); a.click(); document.body.removeChild(a);
}}
// ---------- Compare tool ----------
function setCompareMode(mode) {{
    document.querySelectorAll('.mode-btn').forEach(b => {{
        b.classList.toggle('active', b.dataset.mode === mode);
    }});
    document.querySelectorAll('.compare-mode').forEach(m => m.classList.remove('active'));
    const el = document.getElementById('mode-' + mode);
    if (el) el.classList.add('active');
    if (mode === 'cluster') runCluster();
}}

function computeSimilarity(A, B) {{
    // Typing score
    var typingKeys = Object.keys(A.typing);
    var matches = 0;
    typingKeys.forEach(function(k) {{
        if ((A.typing[k] || '') === (B.typing[k] || '')) matches++;
    }});
    var typingPct = typingKeys.length ? (matches / typingKeys.length) * 100 : 0;

    // Gene score per category (Jaccard index)
    var cats = ['AMR', 'Virulence', 'BACMET', 'Plasmids', 'Mutations'];
    var catScores = [];
    cats.forEach(function(cat) {{
        var sa = new Set(A.genes[cat] || []);
        var sb = new Set(B.genes[cat] || []);
        var union = new Set([...sa, ...sb]);
        if (union.size === 0) return;
        var shared = 0;
        union.forEach(function(g) {{ if (sa.has(g) && sb.has(g)) shared++; }});
        catScores.push((shared / union.size) * 100);
    }});
    var genePct = catScores.length
        ? catScores.reduce((a, b) => a + b, 0) / catScores.length
        : 0;

    var overall = (typingPct * 0.5) + (genePct * 0.5);
    return {{
        typing: Math.round(typingPct),
        genes: Math.round(genePct),
        overall: Math.round(overall),
        typingMatches: matches,
        typingTotal: typingKeys.length
    }};
}}

function similarityColor(pct) {{
    // Red (0) → Amber (50) → Green (100)
    if (pct >= 95) return '#059669';
    if (pct >= 90) return '#10b981';
    if (pct >= 80) return '#84cc16';
    if (pct >= 70) return '#eab308';
    if (pct >= 50) return '#f59e0b';
    if (pct >= 30) return '#f97316';
    return '#dc2626';
}}

function drawGauge(pct) {{
    var color = similarityColor(pct);
    var radius = 60;
    var circumference = 2 * Math.PI * radius;
    var offset = circumference * (1 - pct / 100);
    var svg = '<svg viewBox="0 0 160 160" class="similarity-gauge">' +
        '<circle cx="80" cy="80" r="' + radius + '" fill="none" ' +
        'stroke="#e2e8f0" stroke-width="14"/>' +
        '<circle cx="80" cy="80" r="' + radius + '" fill="none" ' +
        'stroke="' + color + '" stroke-width="14" ' +
        'stroke-dasharray="' + circumference + '" ' +
        'stroke-dashoffset="' + offset + '" ' +
        'stroke-linecap="round" ' +
        'transform="rotate(-90 80 80)"/>' +
        '<text x="80" y="76" text-anchor="middle" ' +
        'font-size="28" font-weight="700" fill="' + color + '">' + pct + '%</text>' +
        '<text x="80" y="98" text-anchor="middle" ' +
        'font-size="11" fill="#64748b">overall</text>' +
        '</svg>';
    return svg;
}}

function typingValuePill(va, vb) {{
    var match = (va === vb);
    var color = match ? '#10b981' : '#dc2626';
    var bg = match ? '#dcfce7' : '#fee2e2';
    var icon = match ? '✓' : '✗';
    return '<span class="typing-pill" style="background:' + bg + ';color:' + color + ';border-color:' + color + ';">' +
        icon + ' ' + escapeHtml(va) + '</span>';
}}

function runCompare() {{
    var data = window.STAPHSCOPE_COMPARE_DATA || {{}};
    var a = document.getElementById('compare-a').value;
    var b = document.getElementById('compare-b').value;
    if (!a || !b) return;
    if (a === b) {{
        document.getElementById('compare-output').innerHTML =
            '<div class="compare-warning">Pick two different samples.</div>';
        document.getElementById('compare-visuals').innerHTML = '';
        document.getElementById('compare-verdict').style.display = 'none';
        return;
    }}
    if (!data[a] || !data[b]) {{
        document.getElementById('compare-output').innerHTML =
            '<div class="compare-warning">Sample data not found.</div>';
        return;
    }}

    var A = data[a], B = data[b];
    var sim = computeSimilarity(A, B);

    // ---------- VISUALS PANEL ----------
    var visualsHtml =
        '<div class="compare-visuals-grid">' +
            '<div class="visual-card gauge-card">' + drawGauge(sim.overall) + '</div>' +
            '<div class="visual-card metrics-card">' +
                '<div class="metric-block">' +
                    '<div class="metric-label">Typing match</div>' +
                    '<div class="metric-bar-track"><div class="metric-bar-fill" ' +
                    'style="width:' + sim.typing + '%;background:' + similarityColor(sim.typing) + ';"></div></div>' +
                    '<div class="metric-value">' + sim.typingMatches + '/' + sim.typingTotal + ' fields (' + sim.typing + '%)</div>' +
                '</div>' +
                '<div class="metric-block">' +
                    '<div class="metric-label">Gene content match</div>' +
                    '<div class="metric-bar-track"><div class="metric-bar-fill" ' +
                    'style="width:' + sim.genes + '%;background:' + similarityColor(sim.genes) + ';"></div></div>' +
                    '<div class="metric-value">' + sim.genes + '% Jaccard similarity</div>' +
                '</div>' +
            '</div>' +
        '</div>';
    document.getElementById('compare-visuals').innerHTML = visualsHtml;

    // ---------- TYPING TABLE (colored pills) ----------
    var typingKeys = Object.keys(A.typing);
    var typingRows = '';
    var typingDiff = 0;
    typingKeys.forEach(function(k) {{
        var va = A.typing[k] || 'Not Assigned';
        var vb = B.typing[k] || 'Not Assigned';
        var diff = va !== vb;
        if (diff) typingDiff++;
        typingRows += '<div class="compare-row ' + (diff ? 'diff' : 'match') + '">' +
            '<span class="row-label">' + k + '</span>' +
            '<span class="row-val">' + typingValuePill(va, vb) + '</span>' +
            '<span class="row-val">' + typingValuePill(vb, va) + '</span>' +
            '</div>';
    }});

    // ---------- GENE CONTENT ----------
    var geneCategories = ['AMR', 'Virulence', 'BACMET', 'Plasmids', 'Mutations'];
    var geneHtml = '';
    geneCategories.forEach(function(cat) {{
        var ga = new Set(A.genes[cat] || []);
        var gb = new Set(B.genes[cat] || []);
        var union = new Set([...ga, ...gb]);
        if (union.size === 0) return;
        var shared = [], onlyA = [], onlyB = [];
        union.forEach(function(g) {{
            if (ga.has(g) && gb.has(g)) shared.push(g);
            else if (ga.has(g)) onlyA.push(g);
            else onlyB.push(g);
        }});
        var total = union.size;
        var sharedPct = (shared.length / total * 100);
        var aPct = (onlyA.length / total * 100);
        var bPct = (onlyB.length / total * 100);

        geneHtml += '<div class="compare-gene-block">' +
            '<h4>' + cat + ' <span class="gene-count">' +
            shared.length + ' shared · ' + onlyA.length + ' only A · ' +
            onlyB.length + ' only B (' + sharedPct.toFixed(0) + '% match)</span></h4>' +
            // Visual proportion bar
            '<div class="gene-proportion-bar">' +
                '<div class="gene-prop-segment shared" style="width:' + sharedPct + '%;" ' +
                'title="Shared: ' + shared.length + '"></div>' +
                '<div class="gene-prop-segment only-a" style="width:' + aPct + '%;" ' +
                'title="Only A: ' + onlyA.length + '"></div>' +
                '<div class="gene-prop-segment only-b" style="width:' + bPct + '%;" ' +
                'title="Only B: ' + onlyB.length + '"></div>' +
            '</div>' +
            '<div class="gene-columns">';
        geneHtml += '<div class="gene-col">' +
            '<div class="gene-col-head">Shared (' + shared.length + ')</div>' +
            (shared.map(g => '<span class="compare-gene shared">' + escapeHtml(g) + '</span>').join('') || '<em>None</em>') +
            '</div>';
        geneHtml += '<div class="gene-col">' +
            '<div class="gene-col-head">Only in A (' + onlyA.length + ')</div>' +
            (onlyA.map(g => '<span class="compare-gene only-a">' + escapeHtml(g) + '</span>').join('') || '<em>None</em>') +
            '</div>';
        geneHtml += '<div class="gene-col">' +
            '<div class="gene-col-head">Only in B (' + onlyB.length + ')</div>' +
            (onlyB.map(g => '<span class="compare-gene only-b">' + escapeHtml(g) + '</span>').join('') || '<em>None</em>') +
            '</div>';
        geneHtml += '</div></div>';
    }});

    // ---------- FULL PROFILE EXPANDERS ----------
    var fullProfileHtml = '<div class="compare-section">' +
        '<h3>📋 Full Gene Profiles</h3>' +
        '<p style="font-size:0.9em;color:#64748b;margin-bottom:12px;">' +
        'Click any category to expand both samples\\' complete gene lists side by side.</p>';

    geneCategories.forEach(function(cat) {{
        var ga = A.genes[cat] || [];
        var gb = B.genes[cat] || [];
        if (ga.length === 0 && gb.length === 0) return;

        var tagsA = ga.length
            ? ga.map(g => '<span class="genome-tag">' + escapeHtml(g) + '</span>').join('')
            : '<em class="empty-note">None detected</em>';
        var tagsB = gb.length
            ? gb.map(g => '<span class="genome-tag">' + escapeHtml(g) + '</span>').join('')
            : '<em class="empty-note">None detected</em>';

        fullProfileHtml +=
            '<details class="vir-details compare-full-details">' +
                '<summary>' +
                    '<span class="details-cat-name">' + cat + '</span>' +
                    '<span class="details-counts">' +
                        'Sample A: ' + ga.length + ' · Sample B: ' + gb.length +
                    '</span>' +
                '</summary>' +
                '<div class="compare-full-profile">' +
                    '<div class="profile-col">' +
                        '<div class="profile-head">' +
                            '<span class="profile-sample-label">' + escapeHtml(a) + '</span>' +
                            '<span class="profile-sample-count">' + ga.length + ' gene(s)</span>' +
                        '</div>' +
                        '<div class="genome-list">' + tagsA + '</div>' +
                    '</div>' +
                    '<div class="profile-col">' +
                        '<div class="profile-head">' +
                            '<span class="profile-sample-label">' + escapeHtml(b) + '</span>' +
                            '<span class="profile-sample-count">' + gb.length + ' gene(s)</span>' +
                        '</div>' +
                        '<div class="genome-list">' + tagsB + '</div>' +
                    '</div>' +
                '</div>' +
            '</details>';
    }});
    fullProfileHtml += '</div>';

    // ---------- VERDICT BANNER ----------
    var verdict = document.getElementById('compare-verdict');
    if (sim.overall >= 95) {{
        verdict.innerHTML = '<div class="compare-verdict-box danger">' +
            '<i class="fas fa-exclamation-triangle"></i> ' +
            '<strong>Near-identical profile (' + sim.overall + '%)</strong> — ' +
            'these samples are likely the same strain. Investigate as a possible <em>transmission pair</em> ' +
            'or a duplicate sample.</div>';
    }} else if (sim.overall >= 85) {{
        verdict.innerHTML = '<div class="compare-verdict-box warning">' +
            '<i class="fas fa-info-circle"></i> ' +
            '<strong>High similarity (' + sim.overall + '%)</strong> — ' +
            'very likely the same clonal lineage with a small number of differing markers.</div>';
    }} else if (sim.overall >= 60) {{
        verdict.innerHTML = '<div class="compare-verdict-box warning">' +
            '<i class="fas fa-info-circle"></i> ' +
            '<strong>Moderate similarity (' + sim.overall + '%)</strong> — ' +
            'related lineage, but not a recent transmission pair.</div>';
    }} else {{
        verdict.innerHTML = '<div class="compare-verdict-box info">' +
            '<i class="fas fa-info-circle"></i> ' +
            '<strong>Low similarity (' + sim.overall + '%)</strong> — ' +
            'distinct lineages.</div>';
    }}
    verdict.style.display = 'block';

    var tableHtml = '<div class="compare-section">' +
        '<h3>🧬 Typing Profile</h3>' +
        '<div class="compare-table">' +
        '<div class="compare-row header">' +
        '<span class="row-label">Field</span>' +
        '<span class="row-val"><strong>' + escapeHtml(a) + '</strong></span>' +
        '<span class="row-val"><strong>' + escapeHtml(b) + '</strong></span>' +
        '</div>' + typingRows + '</div></div>';

    document.getElementById('compare-output').innerHTML =
        tableHtml +
        '<div class="compare-section"><h3>💊 Gene Content</h3>' +
        (geneHtml || '<p><em>No gene data available for these samples.</em></p>') +
        '</div>' + fullProfileHtml;
}}

function swapCompare() {{
    var a = document.getElementById('compare-a');
    var b = document.getElementById('compare-b');
    var tmp = a.value;
    a.value = b.value;
    b.value = tmp;
    runCompare();
}}

function resetCompare() {{
    var sel = document.getElementById('compare-a');
    if (!sel) return;
    var opts = sel.options;
    if (opts.length >= 2) {{
        document.getElementById('compare-a').value = opts[0].value;
        document.getElementById('compare-b').value = opts[1].value;
    }}
    document.getElementById('compare-output').innerHTML = '';
    document.getElementById('compare-visuals').innerHTML = '';
    document.getElementById('compare-verdict').style.display = 'none';
}}

// ---------- CLUSTER MODE ----------
function runCluster() {{
    var data = window.STAPHSCOPE_COMPARE_DATA || {{}};
    var samples = Object.keys(data).sort();
    if (samples.length < 2) {{
        document.getElementById('cluster-heatmap').innerHTML =
            '<div class="compare-warning">Need at least two samples.</div>';
        return;
    }}

    var thresholdEl = document.getElementById('cluster-threshold');
    var threshold = thresholdEl ? parseInt(thresholdEl.value) : 90;

    // Compute N×N similarity matrix
    var matrix = {{}};
    samples.forEach(function(s) {{ matrix[s] = {{}}; }});
    samples.forEach(function(s1) {{
        samples.forEach(function(s2) {{
            if (s1 === s2) {{
                matrix[s1][s2] = 100;
            }} else {{
                var sim = computeSimilarity(data[s1], data[s2]);
                matrix[s1][s2] = sim.overall;
            }}
        }});
    }});

    // Detect clusters using union-find
    var parent = {{}};
    samples.forEach(function(s) {{ parent[s] = s; }});
    function find(x) {{
        while (parent[x] !== x) {{
            parent[x] = parent[parent[x]];
            x = parent[x];
        }}
        return x;
    }}
    function union(a, b) {{
        var ra = find(a), rb = find(b);
        if (ra !== rb) parent[ra] = rb;
    }}
    for (var i = 0; i < samples.length; i++) {{
        for (var j = i + 1; j < samples.length; j++) {{
            if (matrix[samples[i]][samples[j]] >= threshold) {{
                union(samples[i], samples[j]);
            }}
        }}
    }}

    // Group samples by root
    var groups = {{}};
    samples.forEach(function(s) {{
        var r = find(s);
        if (!groups[r]) groups[r] = [];
        groups[r].push(s);
    }});

    // Only keep clusters with ≥2 samples
    var clusters = Object.values(groups).filter(function(g) {{ return g.length >= 2; }});
    clusters.sort(function(a, b) {{ return b.length - a.length; }});

    // ---- Cluster summary cards ----
    var clusterHtml = '';
    if (clusters.length === 0) {{
        clusterHtml = '<div class="no-alerts">' +
            'No clusters detected at ≥' + threshold + '% similarity. ' +
            'Every sample is sufficiently distinct.</div>';
    }} else {{
        clusters.forEach(function(cluster, idx) {{
            // Compute mean intra-cluster similarity
            var sum = 0, count = 0;
            for (var i = 0; i < cluster.length; i++) {{
                for (var j = i + 1; j < cluster.length; j++) {{
                    sum += matrix[cluster[i]][cluster[j]];
                    count++;
                }}
            }}
            var mean = count ? (sum / count).toFixed(1) : '100';
            var color = similarityColor(parseFloat(mean));
            clusterHtml += '<div class="cluster-card" style="border-left-color:' + color + ';">' +
                '<div class="cluster-header">' +
                    '<span class="cluster-badge" style="background:' + color + ';">' +
                    'Cluster ' + (idx + 1) + '</span>' +
                    '<span class="cluster-stats">' + cluster.length + ' samples · mean ' +
                    mean + '% similar</span>' +
                '</div>' +
                '<div class="cluster-samples">' +
                    cluster.map(s => '<span class="cluster-sample-pill">' + escapeHtml(s) + '</span>').join('') +
                '</div></div>';
        }});
    }}

    // ---- Heatmap ----
    var heatmapHtml = '<div class="heatmap-wrapper">' +
        '<h3>🔷 Pairwise Similarity Matrix</h3>' +
        '<p style="font-size:0.9em;color:#64748b;margin-bottom:12px;">' +
        'Colored by overall similarity: ' +
        '<span style="color:#059669;font-weight:700;">green ≥95%</span>, ' +
        '<span style="color:#eab308;font-weight:700;">yellow 70–80%</span>, ' +
        '<span style="color:#dc2626;font-weight:700;">red &lt;50%</span>. ' +
        'Hover any cell for details.</p>' +
        '<div class="heatmap-grid" style="grid-template-columns: 180px repeat(' +
        samples.length + ', 32px);">';

    // Header row
    heatmapHtml += '<div class="heatmap-cell heatmap-corner"></div>';
    samples.forEach(function(s) {{
        heatmapHtml += '<div class="heatmap-cell heatmap-col-header" title="' + escapeHtml(s) + '">' +
            '<span class="heatmap-header-label">' + escapeHtml(s.slice(-6)) + '</span></div>';
    }});

    // Rows
    samples.forEach(function(s1) {{
        heatmapHtml += '<div class="heatmap-cell heatmap-row-header" title="' + escapeHtml(s1) + '">' +
            '<span class="heatmap-header-label">' + escapeHtml(s1) + '</span></div>';
        samples.forEach(function(s2) {{
            var pct = matrix[s1][s2];
            var color = (s1 === s2) ? '#1e293b' : similarityColor(pct);
            var text = (s1 === s2) ? '—' : Math.round(pct);
            heatmapHtml += '<div class="heatmap-cell" style="background:' + color + ';color:white;" ' +
                'title="' + escapeHtml(s1) + ' vs ' + escapeHtml(s2) + ': ' + pct + '%">' +
                text + '</div>';
        }});
    }});
    heatmapHtml += '</div></div>';

    document.getElementById('cluster-clusters').innerHTML = clusterHtml;
    document.getElementById('cluster-heatmap').innerHTML = heatmapHtml;

    var summary = document.getElementById('cluster-summary');
    if (summary) {{
        summary.textContent = clusters.length + ' cluster(s) detected at ≥' + threshold + '%';
    }}
}}
document.addEventListener('DOMContentLoaded', function(){{
    var hash = window.location.hash.substring(1);
    var target = hash ? document.querySelector('.tab-button.'+hash) : document.querySelector('.tab-button');
    if(target) target.click();
    document.querySelectorAll('.data-table').forEach(function(table){{
        var headers = table.querySelectorAll('th');
        headers.forEach(function(h, idx){{
            var type = h.getAttribute('data-sort') || 'string';
            h.style.cursor = 'pointer';
            h.addEventListener('click', function(){{ sortTable(table.id, idx, type); }});
            var icon = document.createElement('span');
            icon.className='sort-icon'; icon.innerHTML='⇅';
            h.appendChild(icon);
        }});
    }});
    document.querySelectorAll('.accordion-header').forEach(function(header){{
        header.addEventListener('click', function(){{
            var content = this.nextElementSibling;
            content.style.display = content.style.display === 'block' ? 'none' : 'block';
        }});
    }});
    document.querySelectorAll('.copy-btn').forEach(function(b){{
        b.addEventListener('click', function(){{
            var c = this.getAttribute('data-citation') || '';
            navigator.clipboard.writeText(c).then(()=>{{
                var t = this.innerHTML;
                this.innerHTML='✓ Copied!';
                setTimeout(()=>{{this.innerHTML=t;}},2000);
            }});
        }});
    }});
}});
</script>"""

    # -------------------------------------------------------------------------
    # Main report assembly
    # -------------------------------------------------------------------------
    def generate_main_report(self, integrated_data: Dict[str, Any], output_dir: Path) -> str:
        """Assemble the full HTML report and write it to disk."""
        print("\n🎨 Generating STAPHSCOPE ULTIMATE HTML report...")
        samples_data = integrated_data.get('samples', {})
        patterns = integrated_data.get('patterns', {})
        gene_centric = integrated_data.get('gene_centric', {})
        metadata = integrated_data.get('metadata', {})
        self._current_samples_data = samples_data

        sample_typing_js = {}
        for sample, data in samples_data.items():
            t = data.get('typing', {})
            sample_typing_js[sample] = {
                "MLST": t.get('MLST', 'Not Assigned'),
                "spa": t.get('spa_Type', 'Not Assigned'),
                "SCCmec_CGE": t.get('SCCmec_CGE', 'Not Assigned'),
                "SCCmec_RPet": t.get('SCCmec_RPet', 'Not Assigned'),
                "SCCmec_Subtype": t.get('SCCmec_Subtype', 'Not Assigned'),
                "Capsule": t.get('capsule_type', 'Not Assigned'),
                "agr": t.get('agr_Type', 'Not Assigned'),
                "MRSA": t.get('MRSA_Status', 'Not Assigned'),
            }

        css = self._get_css()
        js = self._get_js(json.dumps(sample_typing_js))

        TABS = [
            ('summary', 'Summary', 'chart-pie'),
            ('sample_overview', 'Sample Overview', 'list-alt'),
            ('qc', 'FASTA QC', 'chart-line'),
            ('mlst', 'MLST', 'code-branch'),
            ('spa', 'spa Typing', 'dna'),
            ('sccmec', 'SCCmec', 'shield-alt'),
            ('capsule', 'Capsule', 'capsules'),
            ('mrsa', 'MRSA Analysis', 'skull-crossbones'),
            ('agr', 'agr Typing', 'dna'),
            ('amr', 'AMR', 'biohazard'),
            ('virulence', 'Virulence', 'virus'),
            ('bacmet', 'BACMET', 'flask'),
            ('plasmids', 'Plasmids', 'plug'),
            ('mutation', 'Mutations', 'dna'),
            ('mge', 'MGE', 'mobile-alt'),
            ('patterns', 'Patterns', 'project-diagram'),
            ('compare', 'Compare', 'balance-scale'),
            ('aiguide', 'AI Guide', 'robot'),
            ('calltoaction', 'Call to Action', 'globe'),
            ('citation', 'Citation', 'book'),
            ('funding', 'Funding', 'coffee'),
            ('export', 'Export', 'download'),
        ]

        nav_html = ''
        for i, (tid, title, icon) in enumerate(TABS):
            active = ' active' if i == 0 else ''
            nav_html += (f'<button class="tab-button {tid}{active}" onclick="switchTab(\'{tid}\')">'
                         f'<i class="fas fa-{icon}"></i> {title}</button>')

        section_methods = {
            'summary': self._sec_summary,
            'sample_overview': self._sec_sample_overview,
            'qc': self._sec_qc,
            'mlst': self._sec_mlst,
            'spa': self._sec_spa,
            'sccmec': self._sec_sccmec,
            'capsule': self._sec_capsule,
            'mrsa': self._sec_mrsa,
            'agr': self._sec_agr,
            'amr': self._sec_amr,
            'virulence': self._sec_virulence,
            'bacmet': self._sec_bacmet,
            'plasmids': self._sec_plasmids,
            'mutation': self._sec_mutation,
            'mge': self._sec_mge,
            'patterns': self._sec_patterns,
            'compare': self._sec_compare,
            'aiguide': self._sec_aiguide,
            'calltoaction': self._sec_calltoaction,
            'citation': self._sec_citation,
            'funding': self._sec_funding,
            'export': self._sec_export,
        }

        tabs_html = ''
        for i, (tid, title, icon) in enumerate(TABS):
            active = ' active' if i == 0 else ''
            color = self.tab_colors.get(tid, '#4CAF50')
            content = section_methods[tid](integrated_data)
            tabs_html += f'''
            <div id="{tid}-tab" class="tab-content{active}">
                <h2 class="section-header {tid}-header" style="border-color:{color};">
                    <span><i class="fas fa-{icon}"></i> {title}</span>
                    <button class="print-section-btn" onclick="printSection('{tid}-tab')">
                        <i class="fas fa-print"></i> Print</button>
                </h2>
                {content}
            </div>'''

        total_samples = len(samples_data)
        n_mlst = len(patterns.get('mlst_distribution', {}))
        n_spa = len(patterns.get('spa_type_distribution', {}))
        n_amr = sum(len(g) for g in gene_centric.get('amr_databases', {}).values())
        n_vir = sum(len(g) for g in gene_centric.get('virulence_databases', {}).values())
        n_agr = len(patterns.get('agr_type_distribution', {}))
        n_mge = len(integrated_data.get('mge_data', {}).get('per_sample', {}))

        dash_html = f'''
        <div class="dashboard-grid">
            <div class="dashboard-card card-summary" onclick="switchTab('summary')">
                <i class="fas fa-vial fa-2x" style="color:#4CAF50;"></i>
                <div class="card-number">{total_samples}</div>
                <div class="card-label">Total Samples</div></div>
            <div class="dashboard-card card-mlst" onclick="switchTab('mlst')">
                <i class="fas fa-code-branch fa-2x" style="color:#FF9800;"></i>
                <div class="card-number">{n_mlst}</div>
                <div class="card-label">Unique STs</div></div>
            <div class="dashboard-card card-spa" onclick="switchTab('spa')">
                <i class="fas fa-dna fa-2x" style="color:#9C27B0;"></i>
                <div class="card-number">{n_spa}</div>
                <div class="card-label">spa Types</div></div>
            <div class="dashboard-card card-amr" onclick="switchTab('amr')">
                <i class="fas fa-biohazard fa-2x" style="color:#F44336;"></i>
                <div class="card-number">{n_amr}</div>
                <div class="card-label">AMR Genes</div></div>
            <div class="dashboard-card card-virulence" onclick="switchTab('virulence')">
                <i class="fas fa-virus fa-2x" style="color:#E91E63;"></i>
                <div class="card-number">{n_vir}</div>
                <div class="card-label">Virulence Genes</div></div>
            <div class="dashboard-card card-agr" onclick="switchTab('agr')">
                <i class="fas fa-dna fa-2x" style="color:#8B5CF6;"></i>
                <div class="card-number">{n_agr}</div>
                <div class="card-label">agr Types</div></div>
            <div class="dashboard-card card-mge" onclick="switchTab('mge')">
                <i class="fas fa-mobile-alt fa-2x" style="color:#16A085;"></i>
                <div class="card-number">{n_mge}</div>
                <div class="card-label">MGE-profiled Samples</div></div>
            <div class="dashboard-card card-patterns" onclick="switchTab('patterns')">
                <i class="fas fa-project-diagram fa-2x" style="color:#3F51B5;"></i>
                <div class="card-number">{len(patterns.get('high_risk_combinations', []))}</div>
                <div class="card-label">High-Risk Combos</div></div>
        </div>'''

        html = f'''<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>STAPHSCOPE Ultimate S. aureus Report v3.0.0</title>
    <link rel="stylesheet" href="https://cdnjs.cloudflare.com/ajax/libs/font-awesome/6.4.0/css/all.min.css">
    {css}
    {js}
</head>
<body>
<div class="container">
    <div class="main-header">
        <h1><i class="fas fa-bacteria"></i> STAPHSCOPE Ultimate S. aureus Analysis Report</h1>
        <p>Gene-Centric Cross-Genome Analysis with Full Typing, MGE, and Interaction</p>
        <div class="metadata-bar">
            <div class="metadata-item"><i class="fas fa-calendar"></i><span>Generated: {metadata.get('analysis_date', 'Unknown')}</span></div>
            <div class="metadata-item"><i class="fas fa-database"></i><span>Samples: {total_samples}</span></div>
            <div class="metadata-item"><i class="fas fa-code-branch"></i><span>Tool: STAPHSCOPE Ultimate v3.0.0</span></div>
            <div class="metadata-item"><i class="fas fa-university"></i><span>University of Ghana Medical School</span></div>
        </div>
    </div>
    {dash_html}
    <div class="tab-navigation">{nav_html}</div>
    {tabs_html}
    <div class="footer">
        <h3>STAPHSCOPE Ultimate S. aureus Reporter v3.0.0</h3>
        <p>University of Ghana Medical School | Brown Beckley &lt;brownbeckley94@gmail.com&gt;</p>
        <p>Generated on {metadata.get('analysis_date', 'Unknown')}</p>
        <p>⭐ Please give a big STAR on GitHub if you found this useful!</p>
    </div>
</div>
</body>
</html>'''
        output_file = output_dir / "staphscope_ultimate_gene_centric_report.html"
        with open(output_file, 'w', encoding='utf-8') as f:
            f.write(html)
        print(f"    ✅ HTML report saved: {output_file}")
        return str(output_file)

    # =========================================================================
    # SUMMARY
    # =========================================================================
    def _sec_summary(self, d):
        """Executive summary with key statistics and feature highlights."""
        samples_data = d['samples']
        patterns = d['patterns']
        gene_centric = d['gene_centric']
        total = len(samples_data)
        total_amr = sum(len(g) for g in gene_centric.get('amr_databases', {}).values())
        total_vir = sum(len(g) for g in gene_centric.get('virulence_databases', {}).values())
        total_bac = sum(len(g) for g in gene_centric.get('bacmet_databases', {}).values())
        total_pla = sum(len(g) for g in gene_centric.get('plasmid_databases', {}).values())
        high_risk = len(patterns.get('high_risk_combinations', []))
        mrsa = sum(1 for s in samples_data.values()
                   if 'MRSA' in s.get('typing', {}).get('MRSA_Status', ''))
        mssa = sum(1 for s in samples_data.values()
                   if 'MSSA' in s.get('typing', {}).get('MRSA_Status', ''))
        agr_dist = patterns.get('agr_type_distribution', {})
        agr_str = ', '.join(f"{k}: {v}" for k, v in sorted(agr_dist.items())) or 'None'

        features = [
            ('Gene-Centric Tables', 'fa-gene',
             'Each AMR, virulence, BACMET, or mutation is shown with <strong>all genomes</strong> that carry it.'),
            ('Dynamic Grouping by Typing', 'fa-layer-group',
             'Group by MLST, spa, SCCmec, agr, capsule, subtype, or any combination.'),
            ('Capsule &amp; SCCmec Subtype', 'fa-capsules',
             'Two additional typing layers with full cross-tabulation.'),
            ('MGE Profiling', 'fa-mobile-alt',
             'mobileOG-db hit counts per sample, plus aggregate profiles by typing.'),
            ('MRSA Focus', 'fa-skull-crossbones',
             'Dedicated MRSA analysis with all typing combinations.'),
            ('BACMET', 'fa-flask',
             'Biocide and heavy-metal resistance genes — hospital-environment adaptation.'),
            ('Section-Specific Printing', 'fa-print',
             'Print any tab individually with one click.'),
            ('Full Data Export', 'fa-download',
             'All tables as CSV; complete dataset as JSON for AI-assisted analysis.'),
            ('AI Assistant Guide', 'fa-robot',
             'Upload JSON to ChatGPT, Claude, or Gemini for natural-language exploration.'),
        ]
        feature_tiles = ''
        for label, icon, desc in features:
            feature_tiles += f'''<div class="database-section">
                <h4><i class="fas {icon}"></i> {label}</h4><p>{desc}</p></div>'''

        return f'''
        <div class="alert-box alert-info">
            <i class="fas fa-info-circle fa-2x"></i>
            <div>
                <h3>📊 Analysis Overview – Gene‑Centric Approach with Full Typing Support</h3>
                <p>This report analyses <strong>{total}</strong> <em>Staphylococcus aureus</em> genomes using a <strong>gene‑centric</strong> approach: each resistance or virulence gene is displayed with <strong>all genomes</strong> that carry it — making outbreak tracking and co‑occurrence analysis immediate and transparent.</p>
                <p>Includes MLST, spa, SCCmec (CGE + RPet + Subtype), agr, capsule, MRSA status, AMR, virulence, BACMET, plasmids, point mutations, and mobile genetic element profiles.</p>
            </div>
        </div>
        <div class="alert-box alert-success">
            <i class="fas fa-magic fa-2x"></i>
            <div>
                <h3>📘 How to Use Grouping</h3>
                <ol>
                    <li>Go to any gene‑centric table (AMR, Virulence, BACMET, Plasmids, Mutations).</li>
                    <li>Use the <strong>"Group genomes by"</strong> dropdown above the table.</li>
                    <li>Pick any single field or combination — the genome column re‑organises instantly.</li>
                    <li>Click <strong>"Reset"</strong> to return to a flat list.</li>
                </ol>
                <p><strong>Why it matters:</strong> you can instantly ask “Does <em>mecA</em> only appear in ST5?” or “Which agr types carry PVL?”</p>
            </div>
        </div>
        <h3><i class="fas fa-chart-bar"></i> Key Statistics</h3>
        <div class="scrollable-table"><table class="data-table">
            <thead><tr><th>Metric</th><th>Count</th><th>Details</th></tr></thead>
            <tbody>
                <tr><td>Total Samples Analysed</td><td><strong>{total}</strong></td><td>Complete genomic analysis with all databases</td></tr>
                <tr><td>MRSA Samples</td><td><span class="badge badge-mrsa">{mrsa}</span></td><td>Methicillin‑resistant S. aureus</td></tr>
                <tr><td>MSSA Samples</td><td><span class="badge badge-mssa">{mssa}</span></td><td>Methicillin‑sensitive S. aureus</td></tr>
                <tr><td>Unique MLST Types</td><td><strong>{len(patterns.get('mlst_distribution', {}))}</strong></td><td>Sequence types (population structure)</td></tr>
                <tr><td>Unique spa Types</td><td><strong>{len(patterns.get('spa_type_distribution', {}))}</strong></td><td>Protein A gene typing</td></tr>
                <tr><td>Unique SCCmec (CGE)</td><td><strong>{len(patterns.get('sccmec_cge_distribution', {}))}</strong></td><td>SCCmec cassette (CGE caller)</td></tr>
                <tr><td>Unique SCCmec (RPet)</td><td><strong>{len(patterns.get('sccmec_rpet_distribution', {}))}</strong></td><td>SCCmec cassette (RPet caller)</td></tr>
                <tr><td>Unique SCCmec Subtypes</td><td><strong>{len(patterns.get('sccmec_subtype_distribution', {}))}</strong></td><td>Fine-grained cassette subtypes</td></tr>
                <tr><td>Unique Capsule Types</td><td><strong>{len(patterns.get('capsule_distribution', {}))}</strong></td><td>Serotype / vaccine markers</td></tr>
                <tr><td>agr Types Detected</td><td><strong>{len(agr_dist)}</strong></td><td>Distribution: {agr_str}</td></tr>
                <tr><td>AMR Genes</td><td><strong>{total_amr}</strong></td><td>Across all AMR databases</td></tr>
                <tr><td>Virulence Genes</td><td><strong>{total_vir}</strong></td><td>From VFDB</td></tr>
                <tr><td>BACMET Genes</td><td><strong>{total_bac}</strong></td><td>Biocide and heavy metal resistance</td></tr>
                <tr><td>Plasmid Replicons</td><td><strong>{total_pla}</strong></td><td>Plasmid families</td></tr>
                <tr><td>High‑Risk AMR+Virulence Combos</td><td><span class="badge badge-critical">{high_risk}</span></td><td>Samples with both critical AMR and virulence genes</td></tr>
            </tbody>
        </table></div>
        <h3 style="margin-top:30px;"><i class="fas fa-lightbulb"></i> ✨ Key Features</h3>
        <div style="display:grid;grid-template-columns:repeat(auto-fit,minmax(300px,1fr));gap:20px;margin:20px 0;">
            {feature_tiles}
        </div>'''

    # =========================================================================
    # SAMPLE OVERVIEW
    # =========================================================================
    def _sec_sample_overview(self, d):
        """Master per-sample table with all typing columns and gene expanders."""
        samples_data = d['samples']
        rows = ''
        for sample, data in sorted(samples_data.items()):
            t = data.get('typing', {})
            mlst = t.get('MLST', 'Not Assigned')
            spa = t.get('spa_Type', 'Not Assigned')
            agr = t.get('agr_Type', 'Not Assigned')
            cap = t.get('capsule_type', 'Not Assigned')
            cge = t.get('SCCmec_CGE', 'Not Assigned')
            rpet = t.get('SCCmec_RPet', 'Not Assigned')
            sub = t.get('SCCmec_Subtype', 'Not Assigned')
            mrsa = t.get('MRSA_Status', 'Not Assigned')

            vir_genes = data.get('abricate_databases', {}).get('vfdb', [])
            vir_count = len(vir_genes)
            mge_hits = data.get('mge_hits', 0)

            row_class = 'class="mrsa-highlight"' if 'MRSA' in mrsa else ''
            status_badge = ('<span class="badge badge-mrsa">MRSA</span>' if 'MRSA' in mrsa
                            else ('<span class="badge badge-mssa">MSSA</span>' if 'MSSA' in mrsa
                                  else mrsa))
            agr_class = f"agr-{agr}" if agr in ('I', 'II', 'III', 'IV') else 'agr-NA'

            vir_tags = ''.join(f'<span class="genome-tag">{esc(g)}</span>' for g in vir_genes)
            vir_details = (
                f'<details class="vir-details"><summary>{vir_count} gene(s)</summary>'
                f'<div class="vir-gene-tags">{vir_tags or "<em>None</em>"}</div></details>'
            )

            rows += f'''<tr {row_class}>
                <td><strong>{esc(sample)}</strong></td>
                <td>{esc(mlst)}</td>
                <td>{esc(spa)}</td>
                <td><span class="typing-badge {agr_class}">{esc(agr)}</span></td>
                <td>{self._colorize_capsule_cell(cap)}</td>
                <td>{esc(cge)}</td>
                <td>{esc(rpet)}</td>
                <td>{esc(sub)}</td>
                <td>{status_badge}</td>
                <td>{vir_details}</td>
                <td>{mge_hits}</td>
            </tr>'''

        return f'''
        <div class="alert-box alert-info">
            <i class="fas fa-info-circle fa-2x"></i>
            <div>
                <h3>🧬 Sample Overview – The Population Snapshot</h3>
                <p>This is your <strong>master reference table</strong>. Every downstream tab breaks down one slice of the dataset — this one shows all typing layers for every isolate on a single row, so you can spot clonal clusters, resistance-linked lineages, and outbreak candidates at a glance.</p>
            </div>
        </div>

        <div class="database-section">
            <h4><i class="fas fa-book-medical"></i> What Each Column Means</h4>
            <div style="display:grid;grid-template-columns:repeat(auto-fit,minmax(340px,1fr));gap:15px;margin-top:10px;">

                <div style="background:#fff;padding:14px;border-radius:8px;border-left:4px solid #FF9800;">
                    <strong style="color:#FF9800;">🔬 MLST (Multi-Locus Sequence Typing)</strong>
                    <p style="font-size:.92em;margin-top:6px;">Sequences seven conserved <em>housekeeping genes</em>. Each unique allele combination defines a <strong>Sequence Type (ST)</strong>. STs are the <em>lingua franca</em> of global <em>S. aureus</em> epidemiology — they reveal the deep population structure and let you compare your isolates to thousands worldwide. <strong>Example:</strong> ST5, ST8, ST239, ST398, ST22.</p>
                </div>

                <div style="background:#fff;padding:14px;border-radius:8px;border-left:4px solid #9C27B0;">
                    <strong style="color:#9C27B0;">🧬 spa Typing</strong>
                    <p style="font-size:.92em;margin-top:6px;">Sequences the variable repeat region of the <em>spa</em> gene (protein A). Gives <strong>outbreak-resolution discrimination</strong> within a single ST — two isolates from the same ST5 lineage might carry different spa types (e.g. t002 vs t688), telling you they came from different introductions. Used worldwide by reference labs and PulseNet-style surveillance.</p>
                </div>

                <div style="background:#fff;padding:14px;border-radius:8px;border-left:4px solid #009688;">
                    <strong style="color:#009688;">🛡️ SCCmec (CGE &amp; RPet callers)</strong>
                    <p style="font-size:.92em;margin-top:6px;">The Staphylococcal Cassette Chromosome <em>mec</em> carries <strong>mecA</strong> or <strong>mecC</strong> — the genetic determinant of methicillin resistance. Its <strong>type</strong> (I–XIII) and <strong>subtype</strong> (IIa, IVc, …) trace MRSA lineage origin. Two callers (CGE and RPet) are shown side-by-side so you can see where they agree or disagree. Subtype is the finest grain — useful for tracking regional clones.</p>
                </div>

                <div style="background:#fff;padding:14px;border-radius:8px;border-left:4px solid #8B5CF6;">
                    <strong style="color:#8B5CF6;">🧬 agr Typing</strong>
                    <p style="font-size:.92em;margin-top:6px;">The <em>accessory gene regulator</em> system — a quorum-sensing circuit that globally controls virulence gene expression. Four types (I–IV) exist. agr type influences <strong>toxin production, biofilm formation, and infection severity</strong>, and correlates with certain clonal complexes. It's the key axis linking <em>genotype</em> to <em>virulence phenotype</em>.</p>
                </div>

                <div style="background:#fff;padding:14px;border-radius:8px;border-left:4px solid #00ACC1;">
                    <strong style="color:#00ACC1;">💊 Capsule Type</strong>
                    <p style="font-size:.92em;margin-top:6px;">Capsular polysaccharide serotype — mainly <strong>Type 5</strong> and <strong>Type 8</strong>. The capsule helps the bacterium evade phagocytosis and is a major <strong>vaccine target</strong>. Type 5 is more common in MSSA, Type 8 in MRSA ST8 / ST239 lineages. Knowing the capsule type helps interpret vaccine coverage and immune-evasion potential.</p>
                </div>

                <div style="background:#fff;padding:14px;border-radius:8px;border-left:4px solid #795548;">
                    <strong style="color:#795548;">🦠 MRSA / MSSA Status</strong>
                    <p style="font-size:.92em;margin-top:6px;">Methicillin-<strong>R</strong>esistant vs Methicillin-<strong>S</strong>usceptible <em>S. aureus</em>. This is the <strong>clinical bottom line</strong> — MRSA excludes standard β-lactam therapy (methicillin, oxacillin, cefazolin) and forces clinicians to second-line agents (vancomycin, daptomycin, linezolid). MRSA rows are highlighted in red throughout the report.</p>
                </div>
            </div>
        </div>

        <div class="alert-box alert-success">
            <i class="fas fa-microscope fa-2x"></i>
            <div>
                <h3>💡 Why This Matters for AMR Research</h3>
                <ul style="margin-top:8px;">
                    <li><strong>Lineage–resistance linkage.</strong> Antibiotic resistance is not randomly distributed — it clusters with specific clones. Combined typing reveals <em>"ST239-II is our dominant MRSA clone"</em> or <em>"PVL is confined to ST121-ST152"</em>. That is epidemiological intelligence you cannot get from AMR gene presence alone.</li>
                    <li><strong>Outbreak detection.</strong> When two isolates share the same <strong>MLST + spa + SCCmec subtype + agr</strong>, they are almost certainly recent transmissions. Mismatches flag independent introductions.</li>
                    <li><strong>MRSA surveillance.</strong> Which STs/SCCmec types dominate in your setting? Do they match regional or global trends? Is a new clone emerging?</li>
                    <li><strong>Vaccine and therapy design.</strong> Capsule type + agr type + MRSA status → understand which prevention strategies and treatment options are plausible for each lineage.</li>
                    <li><strong>Genomic context for resistance genes.</strong> A <em>mecA</em> hit means more when you know it sits on SCCmec IVa inside an ST22-t032 background — a well-characterised UK/Ireland epidemic clone — vs. an isolated novel cassette.</li>
                    <li><strong>Horizontal transfer signals.</strong> When the same AMR gene appears across <em>unrelated</em> STs, plasmid or transposon spread is likely — a public health red flag.</li>
                </ul>
            </div>
        </div>

        <div class="alert-box alert-info">
            <i class="fas fa-mouse-pointer fa-2x"></i>
            <div>
                <h3>🖱️ How to Use This Tab</h3>
                <ul style="margin-top:8px;">
                    <li><strong>Search</strong> the box above the table — filters by sample ID, ST, spa, agr, capsule, SCCmec, or MRSA status in real time.</li>
                    <li><strong>Sort</strong> any column by clicking its header (⇅) — e.g. group all MRSA isolates together, or sort by MobileOG hits.</li>
                    <li><strong>Expand virulence genes</strong> — click the number in the Virulence column to reveal the full gene list inline, no page reload.</li>
                    <li><strong>Export</strong> as CSV with one click for downstream analysis in R, Python, or Excel.</li>
                    <li><strong>MobileOG Hits</strong> — total count of mobileOG-db protein family matches. Higher = more mobile-element machinery in the genome (plasmids, phages, IS elements). Detailed breakdown is on the <em>MGE</em> tab.</li>
                </ul>
            </div>
        </div>

        <input type="text" class="search-box" id="search-samples"
               onkeyup="searchTable('samples-table', 'search-samples')"
               placeholder="🔍 Search samples by ID, ST, spa, agr, capsule, SCCmec, MRSA status...">
        <div class="action-buttons">
            <button class="action-btn btn-primary" onclick="exportTableToCSV('samples-table', 'sample_overview.csv')"><i class="fas fa-download"></i> Export CSV</button>
            <button class="action-btn btn-light" onclick="document.getElementById('search-samples').value='';searchTable('samples-table','search-samples')"><i class="fas fa-sync"></i> Clear</button>
        </div>
        <div class="scrollable-table">
            <table id="samples-table" class="data-table">
                <thead><tr>
                    <th data-sort="string">Sample</th>
                    <th data-sort="string">MLST</th>
                    <th data-sort="string">spa Type</th>
                    <th data-sort="string">agr Type</th>
                    <th data-sort="string">Capsule Type</th>
                    <th data-sort="string">SCCmec (CGE)</th>
                    <th data-sort="string">SCCmec (RPet)</th>
                    <th data-sort="string">SCCmec Subtype</th>
                    <th data-sort="string">MRSA/MSSA</th>
                    <th data-sort="number">Virulence</th>
                    <th data-sort="number">MobileOG Hits</th>
                </tr></thead>
                <tbody>{rows}</tbody>
            </table>
        </div>'''

    # =========================================================================
    # FASTA QC
    # =========================================================================
    def _sec_qc(self, d):
        """FASTA quality control metrics per sample with fastANI acknowledgment."""
        qc_data = d.get('qc_data', {})
        if not qc_data:
            return self._alert('warning', 'fa-exclamation-circle',
                               '<h3>No QC Data Available</h3>'
                               '<p>The FASTA_QC_summary.html file was not found.</p>')
        all_metrics = set()
        for m in qc_data.values():
            all_metrics.update(m.keys())
        metric_list = sorted(all_metrics)

        credit = self._credit_bar('#17a2b8', '📊', 'FASTA QC Metrics',
            'Computed with <strong>Biopython</strong> — assembly statistics (contigs, N50, GC%, total length) from FASTA files.')

        fastani_credit = self._credit_bar('#17a2b8', '🧬',
            'Species Confirmation – fastANI',
            '<strong>fastANI</strong> developed by '
            '<a href="https://github.com/ParBLiSS/FastANI" target="_blank" style="color:#17a2b8;font-weight:bold;">ParBLiSS (Jain et al.)</a> — '
            'high-throughput Average Nucleotide Identity calculation for species-level identification.<br>'
            '<span style="font-size:0.9em;color:#6c757d;"><i class="fas fa-book-open"></i> '
            'Please cite: Jain C, Rodriguez-R LM, Phillippy AM, Konstantinidis KT, Aluru S. High throughput ANI analysis of 90K prokaryotic genomes reveals clear species boundaries. <em>Nat Commun</em>. 2018;9(1):5114. '
            '<a href="https://doi.org/10.1038/s41467-018-07641-9" target="_blank" style="color:#17a2b8;font-weight:bold;">🔗 DOI</a></span><br>'
            '<span style="font-size:0.9em;color:#6c757d;"><i class="fas fa-gratipay"></i> '
            'We are grateful to the developers for making this tool freely available for species-level identification.</span>')

        header_cells = ''.join(f'<th data-sort="number">{esc(m)}</th>' for m in metric_list)
        rows = ''
        for sample, metrics in sorted(qc_data.items()):
            cells = ''
            for m in metric_list:
                v = metrics.get(m, 'ND')
                if isinstance(v, float):
                    v = f"{v:,.0f}" if v > 1000 else f"{v:.2f}"
                cells += f'<td>{v}</td>'
            rows += f'<tr><td><strong>{esc(sample)}</strong></td>{cells}</tr>'

        return f'''
        {credit}
        {fastani_credit}
        <div class="alert-box alert-info">
            <i class="fas fa-chart-line fa-2x"></i>
            <div>
                <h3>📏 FASTA Quality Control</h3>
                <ul>
                    <li><strong>Number of contigs</strong> – lower is better; typical S. aureus: &lt;200 (good), &lt;100 (excellent).</li>
                    <li><strong>N50</strong> – higher is better; &gt;50 kb (good), &gt;100 kb (excellent).</li>
                    <li><strong>GC%</strong> – S. aureus is typically 32–33%.</li>
                    <li><strong>Total length</strong> – ~2.8 Mbp for a complete genome.</li>
                </ul>
            </div>
        </div>
        <input type="text" class="search-box" id="search-qc"
               onkeyup="searchTable('qc-table', 'search-qc')" placeholder="🔍 Search sample...">
        <div class="action-buttons">
            <button class="action-btn btn-primary" onclick="exportTableToCSV('qc-table', 'fasta_qc.csv')"><i class="fas fa-download"></i> Export QC Data</button>
        </div>
        <div class="master-scrollable-container">
            <table id="qc-table" class="data-table">
                <thead><tr><th data-sort="string">Sample</th>{header_cells}</tr></thead>
                <tbody>{rows}</tbody>
            </table>
        </div>'''

    # =========================================================================
    # MLST
    # =========================================================================
    def _sec_mlst(self, d):
        """MLST distribution and combinations with rich acknowledgements."""
        P = d['patterns']
        credit = self._credit_bar('#FF9800', '🧬',
            'MLST Typing – Acknowledgments &amp; Licensing',
            '<strong>MLST scheme</strong> powered by '
            '<a href="https://github.com/tseemann/mlst" target="_blank" style="color:#FF9800;font-weight:bold;">Prof. Torsten Seemann’s Perl scripts</a> '
            'and the '
            '<a href="https://pubmlst.org/" target="_blank" style="color:#FF9800;font-weight:bold;">PubMLST database</a> '
            '(Jolley et al., <em>Wellcome Open Res</em> 2018).<br>'
            '<span style="color:#856404;"><i class="fas fa-info-circle"></i> '
            '<strong>Note:</strong> Due to recent licensing changes and the new PubMLST data policy, '
            'the allele definitions used in this report are current as of <strong>2024</strong>. '
            'Check staphscope <code>--help</code> for instructions to update the MLST database.</span><br>'
            '<span style="font-size:0.9em;color:#6c757d;"><i class="fas fa-gratipay"></i> '
            'We are grateful to the PubMLST curators and Torsten Seemann for providing these '
            'essential open‑source tools.</span>')
        return f'''
        {credit}
        {self._alert('info', 'fa-code-branch',
            '<h3>🔬 MLST (Multi-Locus Sequence Typing)</h3>'
            '<p>MLST indexes internal fragments of seven housekeeping genes. Each unique combination defines a Sequence Type (ST) — the gold standard for global <em>S. aureus</em> epidemiology.</p>'
            f'<p><strong>{len(P.get("mlst_distribution", {}))} unique STs</strong> identified.</p>')}
        <h3>📊 ST Distribution</h3>
        {self._distribution_table('mlst-dist-table', P.get('mlst_distribution', Counter()), 'ST')}
        {self._combo_table('mlst-spa-table', 'ST – spa Combination', P.get('mlst_spa', {}))}
        {self._combo_table('mlst-scc-table', 'ST – SCCmec (CGE) Combination', P.get('mlst_sccmec_cge', {}))}
        {self._combo_table('mlst-sub-table', 'ST – SCCmec Subtype Combination', P.get('mlst_subtype', {}))}
        {self._combo_table('mlst-agr-table', 'ST – agr Combination', P.get('mlst_agr', {}))}
        {self._combo_table('mlst-cap-table', 'ST – Capsule Combination', P.get('mlst_capsule', {}))}'''

    # =========================================================================
    # SPA
    # =========================================================================
    def _sec_spa(self, d):
        """spa typing distribution and combinations with rich acknowledgements."""
        P = d['patterns']
        credit = self._credit_bar('#9C27B0', '🧬',
            'spa Typing – Acknowledgments',
            '<strong>spa typing</strong> powered by '
            '<a href="https://github.com/mjsull/spa_typing" target="_blank" style="color:#9C27B0;font-weight:bold;">original code by mjsull</a>, '
            'modified by <strong>JFSanchezHerrero</strong>, and the '
            '<a href="https://spa.ridom.de/" target="_blank" style="color:#9C27B0;font-weight:bold;">Ridom SpaServer database</a>.<br>'
            '<span style="font-size:0.9em;color:#6c757d;"><i class="fas fa-gratipay"></i> '
            'We are grateful to the developers and curators for maintaining this essential resource.</span>')
        return f'''
        {credit}
        {self._alert('info', 'fa-dna',
            '<h3>🧬 spa Typing – High-Resolution Outbreak Tracking</h3>'
            '<p>The <em>spa</em> gene encodes protein A; repeat region polymorphisms define spa types. This is one of the most discriminating single-locus typing methods for outbreak investigation.</p>'
            f'<p><strong>{len(P.get("spa_type_distribution", {}))} unique spa types</strong> identified.</p>')}
        <h3>📊 spa Type Distribution</h3>
        {self._distribution_table('spa-dist-table', P.get('spa_type_distribution', Counter()), 'spa Type')}
        {self._combo_table('spa-mlst-table', 'spa – ST Combination', P.get('mlst_spa', {}))}
        {self._combo_table('spa-scc-table', 'spa – SCCmec (CGE) Combination', P.get('spa_sccmec_cge', {}))}
        {self._combo_table('spa-sub-table', 'spa – SCCmec Subtype Combination', P.get('spa_subtype', {}))}
        {self._combo_table('spa-agr-table', 'spa – agr Combination', P.get('spa_agr', {}))}
        {self._combo_table('spa-cap-table', 'spa – Capsule Combination', P.get('spa_capsule', {}))}'''

    # =========================================================================
    # SCCMEC
    # =========================================================================
    def _sec_sccmec(self, d):
        """SCCmec CGE + RPet + Subtype sub-sections with rich acknowledgements for both callers."""
        P = d['patterns']

        cge_credit = self._credit_bar('#009688', '🛡️',
            'SCCmec Typing (CGE) – Acknowledgments',
            '<strong>SCCmecFinder</strong> developed by '
            '<a href="https://cge.cbs.dtu.dk/services/SCCmecFinder/" target="_blank" style="color:#009688;font-weight:bold;">the Center for Genomic Epidemiology (DTU)</a>. '
            'The database is curated by <strong>Anders Rhod Larsen</strong> — we are deeply grateful for his ongoing dedication.<br>'
            '<span style="font-size:0.9em;color:#6c757d;"><i class="fas fa-book-open"></i> '
            'Please cite: Kaya H, et al. SCCmecFinder, a Web‑Based Tool for Typing of Staphylococcal Cassette Chromosome <em>mec</em> in <em>Staphylococcus aureus</em> Using Whole‑Genome Sequence Data. <em>mSphere</em>. 2018;3(1):e00612-17. '
            '<a href="https://doi.org/10.1128/mSphere.00612-17" target="_blank" style="color:#009688;font-weight:bold;">🔗 DOI</a></span><br>'
            '<span style="font-size:0.9em;color:#6c757d;"><i class="fas fa-gratipay"></i> '
            'We thank the developers and curators for making this open‑source tool freely available.</span>')

        rpet_credit = self._credit_bar('#7c3aed', '🔬',
            'SCCmec Typing (RPet) – Acknowledgments',
            '<strong>sccmec</strong> (RPet caller) developed and maintained by '
            '<strong>Robert A. Petit III, PhD</strong> — '
            '<a href="https://github.com/rpetit3/sccmec" target="_blank" style="color:#7c3aed;font-weight:bold;">github.com/rpetit3/sccmec</a>.<br>'
            '<span style="font-size:0.9em;color:#6c757d;"><i class="fas fa-book-open"></i> '
            'Please cite: Petit RA III, Read TD. <em>Staphylococcus aureus</em> viewed from the perspective of 40,000+ genomes. <em>PeerJ</em>. 2018;6:e5261. '
            '<a href="https://doi.org/10.7717/peerj.5261" target="_blank" style="color:#7c3aed;font-weight:bold;">🔗 DOI</a> '
            '— this tool is a standalone, lightweight successor to the original Staphopia‑SCCmec module.</span><br>'
            '<span style="font-size:0.9em;color:#6c757d;"><i class="fas fa-gratipay"></i> '
            'We are deeply grateful to Robert Petit for his sustained contributions to open‑source <em>S. aureus</em> genomics.</span>')

        return f'''
        {cge_credit}
        {self._alert('info', 'fa-shield-alt',
            '<h3>🛡️ SCCmec – The MRSA Cassette</h3>'
            '<p>The staphylococcal cassette chromosome <em>mec</em> (SCCmec) carries <em>mecA</em> or <em>mecC</em>. '
            'Three sub-sections below: <strong>CGE caller</strong>, <strong>RPet caller</strong>, and the fine-grained '
            '<strong>Subtype</strong>. Use them side-by-side to spot caller disagreements.</p>')}

        <h3>📊 SCCmec Type Distribution — CGE caller</h3>
        {self._distribution_table('scc-cge-dist', P.get('sccmec_cge_distribution', Counter()), 'SCCmec (CGE)')}
        {self._combo_table('scc-cge-mlst', 'SCCmec (CGE) – ST', P.get('mlst_sccmec_cge', {}))}
        {self._combo_table('scc-cge-spa', 'SCCmec (CGE) – spa', P.get('spa_sccmec_cge', {}))}
        {self._combo_table('scc-cge-agr', 'SCCmec (CGE) – agr', P.get('agr_sccmec_cge', {}))}
        {self._combo_table('scc-cge-cap', 'SCCmec (CGE) – Capsule', P.get('capsule_sccmec_cge', {}))}

        <h3 style="margin-top:40px;">📊 SCCmec Type Distribution — RPet caller</h3>
        {rpet_credit}
        {self._distribution_table('scc-rpet-dist', P.get('sccmec_rpet_distribution', Counter()), 'SCCmec (RPet)')}
        {self._combo_table('scc-rpet-mlst', 'SCCmec (RPet) – ST', P.get('mlst_sccmec_rpet', {}))}
        {self._combo_table('scc-rpet-spa', 'SCCmec (RPet) – spa', P.get('spa_sccmec_rpet', {}))}
        {self._combo_table('scc-rpet-agr', 'SCCmec (RPet) – agr', P.get('agr_sccmec_rpet', {}))}
        {self._combo_table('scc-rpet-cap', 'SCCmec (RPet) – Capsule', P.get('capsule_sccmec_rpet', {}))}

        <h3 style="margin-top:40px;">📊 SCCmec Subtype Distribution</h3>
        {self._distribution_table('scc-sub-dist', P.get('sccmec_subtype_distribution', Counter()), 'SCCmec Subtype')}
        {self._combo_table('scc-sub-mlst', 'Subtype – ST', P.get('mlst_subtype', {}))}
        {self._combo_table('scc-sub-spa', 'Subtype – spa', P.get('spa_subtype', {}))}
        {self._combo_table('scc-sub-agr', 'Subtype – agr', P.get('agr_subtype', {}))}
        {self._combo_table('scc-sub-cap', 'Subtype – Capsule', P.get('capsule_subtype', {}))}'''

    # =========================================================================
    # CAPSULE
    # =========================================================================
    def _sec_capsule(self, d):
        """Capsule typing distribution and combinations."""
        P = d['patterns']
        return f'''
        {self._alert('info', 'fa-capsules',
            '<h3>💊 Capsule Typing – Serotype / Vaccine Marker</h3>'
            '<p><em>S. aureus</em> capsular polysaccharides (mainly type 5 and type 8) are major vaccine targets and affect phagocytosis evasion and immune recognition.</p>'
            f'<p><strong>{len(P.get("capsule_distribution", {}))} unique capsule types</strong> identified.</p>')}
        <h3>📊 Capsule Type Distribution</h3>
        {self._distribution_table('cap-dist-table', P.get('capsule_distribution', Counter()), 'Capsule Type')}
        {self._combo_table('cap-mlst-table', 'Capsule – ST', P.get('mlst_capsule', {}))}
        {self._combo_table('cap-spa-table', 'Capsule – spa', P.get('spa_capsule', {}))}
        {self._combo_table('cap-agr-table', 'Capsule – agr', P.get('agr_capsule', {}))}
        {self._combo_table('cap-scc-table', 'Capsule – SCCmec (CGE)', P.get('capsule_sccmec_cge', {}))}
        {self._combo_table('cap-sub-table', 'Capsule – SCCmec Subtype', P.get('capsule_subtype', {}))}'''
    # =========================================================================
    # MRSA
    # =========================================================================
    def _sec_mrsa(self, d):
        """Dedicated MRSA analysis with filtered combination tables."""
        P = d['patterns']
        samples_data = d['samples']
        mrsa_samples = [s for s, x in samples_data.items()
                        if 'MRSA' in x.get('typing', {}).get('MRSA_Status', '')]
        mrsa_status_dist = P.get('mrsa_status_distribution', Counter())
        combos = {
            'mrsa-mlst-spa':  P.get('mlst_spa', {}),
            'mrsa-mlst-sub':  P.get('mlst_subtype', {}),
            'mrsa-mlst-agr':  P.get('mlst_agr', {}),
            'mrsa-spa-agr':   P.get('spa_agr', {}),
            'mrsa-spa-cap':   P.get('spa_capsule', {}),
            'mrsa-mlst-cap':  P.get('mlst_capsule', {}),
        }
        html = f'''
        {self._alert('danger', 'fa-skull-crossbones',
            f'<h3>⚠️ MRSA – A Clinical Priority</h3>'
            f'<p><strong>{len(mrsa_samples)} MRSA samples</strong> identified.</p>')}
        <h3>📊 MRSA vs MSSA Distribution</h3>
        {self._distribution_table('mrsa-status-table', mrsa_status_dist, 'Status')}'''
        for tid, title, combo in [
            ('mrsa-mlst-spa-tbl', 'MRSA: ST – spa', combos['mrsa-mlst-spa']),
            ('mrsa-mlst-sub-tbl', 'MRSA: ST – SCCmec Subtype', combos['mrsa-mlst-sub']),
            ('mrsa-mlst-agr-tbl', 'MRSA: ST – agr', combos['mrsa-mlst-agr']),
            ('mrsa-spa-agr-tbl',  'MRSA: spa – agr', combos['mrsa-spa-agr']),
            ('mrsa-spa-cap-tbl',  'MRSA: spa – Capsule', combos['mrsa-spa-cap']),
            ('mrsa-mlst-cap-tbl', 'MRSA: ST – Capsule', combos['mrsa-mlst-cap']),
        ]:
            filtered = {k: v for k, v in combo.items()
                        if any(s in mrsa_samples for s in v)}
            if filtered:
                html += self._combo_table(tid, title, filtered)
        return html

    # =========================================================================
    # AGR
    # =========================================================================
    def _sec_agr(self, d):
        """agr typing distribution and combinations with rich acknowledgements."""
        P = d['patterns']
        credit = self._credit_bar('#8B5CF6', '🧬',
            'agr Typing – Acknowledgments',
            '<strong>AgrVATE</strong> typing tool developed by '
            '<a href="https://github.com/VishnuRaghuram94/AgrVATE" target="_blank" style="color:#8B5CF6;font-weight:bold;">Vishnu Raghuram (VishnuRaghuram94)</a> '
            'and maintained by <strong>Robert A. Petit III</strong> '
            '(<a href="https://github.com/rpetit3" target="_blank" style="color:#8B5CF6;">@rpetit3</a>).<br>'
            '<span style="font-size:0.9em;color:#6c757d;"><i class="fas fa-gratipay"></i> '
            'We are deeply grateful to the developers for making this open‑source tool available '
            'for accurate agr typing of <em>S. aureus</em> genomes.</span>')
        return f'''
        {credit}
        {self._alert('info', 'fa-dna',
            '<h3>🧬 agr Typing – Virulence Regulation</h3>'
            '<p>The accessory gene regulator (<em>agr</em>) system is a quorum-sensing circuit controlling virulence gene expression. Four types are recognised (I–IV); agr type influences toxin production, biofilm formation, and infection severity.</p>')}
        <h3>📊 agr Type Distribution</h3>
        {self._distribution_table('agr-dist-table', P.get('agr_type_distribution', Counter()), 'agr Type')}
        {self._combo_table('agr-mlst-table', 'agr – ST', P.get('mlst_agr', {}))}
        {self._combo_table('agr-spa-table', 'agr – spa', P.get('spa_agr', {}))}
        {self._combo_table('agr-scc-table', 'agr – SCCmec (CGE)', P.get('agr_sccmec_cge', {}))}
        {self._combo_table('agr-sub-table', 'agr – SCCmec Subtype', P.get('agr_subtype', {}))}
        {self._combo_table('agr-cap-table', 'agr – Capsule', P.get('agr_capsule', {}))}'''

    # =========================================================================
    # AMR
    # =========================================================================
    def _sec_amr(self, d):
        """AMR gene-centric tab with multi-database context, caveats, and role guide."""
        gene_centric = d['gene_centric']
        amr_dbs = gene_centric.get('amr_databases', {})
        total_samples = len(d['samples'])
        all_genes = []
        for genes in amr_dbs.values():
            all_genes.extend(genes)
        all_genes.sort(key=lambda x: x['count'], reverse=True)

        credit = self._credit_bar('#dc3545', '🧬',
            'AMR Detection Tools &amp; Databases',
            '<strong>ABRicate</strong> by '
            '<a href="https://github.com/tseemann/abricate" target="_blank" style="color:#dc3545;font-weight:bold;">Prof. Torsten Seemann</a> '
            '→ Powered by databases: '
            '<a href="https://card.mcmaster.ca/" target="_blank" style="color:#28a745;font-weight:600;">CARD</a>, '
            '<a href="https://www.mediterranee-infection.com/amr-databases/" target="_blank" style="color:#17a2b8;font-weight:600;">ARG-ANNOT</a>, '
            '<a href="https://www.ncbi.nlm.nih.gov/pathogens/antimicrobial-resistance/" target="_blank" style="color:#007bff;font-weight:600;">NCBI AMR</a>, '
            '<a href="https://megares.meglab.org/" target="_blank" style="color:#fd7e14;font-weight:600;">MEGARes</a>, '
            '<a href="https://genepi.food.dtu.dk/resfinder" target="_blank" style="color:#20c997;font-weight:600;">ResFinder</a>.<br>'
            '<strong>AMRFinderPlus</strong> by '
            '<a href="https://github.com/ncbi/amr" target="_blank" style="color:#dc3545;font-weight:bold;">NCBI</a> '
            '→ comprehensive resistance gene and point-mutation detection.<br>'
            '<span style="font-size:0.9em;color:#6c757d;"><i class="fas fa-gratipay"></i> '
            'We are deeply grateful to all developers and curators for their open‑source contributions.</span>')

        # ---- Educational box: why multiple databases ----
        why_multi_db = '''
        <div class="alert-box alert-info" style="border-left-color:#00695c; background:#e8f5e9; border-radius:8px; padding:18px 22px; margin:20px 0;">
            <div style="display:flex; gap:15px; align-items:flex-start;">
                <i class="fas fa-info-circle fa-2x" style="color:#00695c; margin-top:3px;"></i>
                <div>
                    <h4 style="margin:0 0 10px 0; color:#00695c; font-size:1.1em;">🔬 Why Do We Use Multiple Databases?</h4>
                    <p style="margin:6px 0; font-size:0.95em; line-height:1.5;">StaphScope screens every genome against <strong>five independent AMR databases</strong> plus <strong>AMRFinderPlus</strong> — because <strong>no single database is comprehensive</strong>. Each has unique strengths, biases, and update cadences.</p>
                    <p style="margin:10px 0 6px 0; font-size:0.95em; line-height:1.5;"><strong>⚠️ The problem with choosing one database:</strong></p>
                    <ul style="margin:6px 0 10px 20px; font-size:0.93em; line-height:1.6; color:#333;">
                        <li>Some researchers pick a favourite database (often CARD or ResFinder) and only report hits from that one.</li>
                        <li>This is <strong>fast but incomplete</strong> — a gene absent from CARD may still be present in MEGARes or AMRFinderPlus.</li>
                        <li>Conversely, single-database hits with weak support can be <strong>false positives</strong> that a second database would have flagged.</li>
                        <li>The result: <em>biased prevalence estimates</em> and <em>missed resistance signals</em> that matter clinically.</li>
                    </ul>
                    <p style="margin:10px 0 6px 0; font-size:0.95em; line-height:1.5;"><strong>✅ Our approach — report everything, provenance preserved:</strong></p>
                    <ul style="margin:6px 0 10px 20px; font-size:0.93em; line-height:1.6; color:#333;">
                        <li>All hits from all databases are kept <strong>separate and unmerged</strong> — the <strong>Database column</strong> tells you exactly which source found each gene.</li>
                        <li>You get the <strong>full picture</strong> and can filter/interpret yourself — no silent filtering, no cherry-picking.</li>
                        <li>Cross-database agreement becomes a <strong>confidence signal</strong> (see next box).</li>
                    </ul>
                    <p style="margin:6px 0 0 0; font-size:0.92em; line-height:1.5; background:#fff3cd; padding:8px 14px; border-radius:4px; border-left:3px solid #ffc107;">
                        <i class="fas fa-lightbulb" style="color:#856404;"></i>
                        <strong>Pro tip:</strong> Use the search box below to filter by database name (e.g. type <code>CARD</code> or <code>ResFinder</code>) to inspect only that database's hits. Group by typing to see which lineages carry which resistance signatures.
                    </p>
                </div>
            </div>
        </div>
        '''

        # ---- Confidence tiers box ----
        confidence_box = '''
        <div class="alert-box" style="border-left-color:#0891b2; background:#e0f2fe; border-radius:8px; padding:16px 20px; margin:20px 0;">
            <div style="display:flex; gap:15px; align-items:flex-start;">
                <i class="fas fa-layer-group fa-2x" style="color:#0891b2; margin-top:3px;"></i>
                <div>
                    <h4 style="margin:0 0 10px 0; color:#0891b2; font-size:1.05em;">🎯 Cross-Database Confidence Tiers</h4>
                    <p style="margin:6px 0; font-size:0.93em; line-height:1.5;">A gene detected by <strong>multiple databases</strong> is far more likely to be a true positive than a single-database hit.</p>
                    <ul style="margin:6px 0 0 20px; font-size:0.93em; line-height:1.6; color:#333;">
                        <li><span style="color:#16a34a; font-weight:bold;">🟢 High confidence</span> — gene found in <strong>3 or more databases</strong>.</li>
                        <li><span style="color:#f59e0b; font-weight:bold;">🟡 Moderate confidence</span> — gene found in <strong>2 databases</strong>.</li>
                        <li><span style="color:#dc2626; font-weight:bold;">🔴 Low confidence / investigate</span> — gene found in <strong>only 1 database</strong>. Could be a genuine novel or niche determinant, or a database-specific artefact.</li>
                    </ul>
                </div>
            </div>
        </div>
        '''

        # ---- Acquired vs intrinsic box ----
        acquired_intrinsic = '''
        <div class="alert-box" style="border-left-color:#6f42c1; background:#f3e8ff; border-radius:8px; padding:16px 20px; margin:20px 0;">
            <div style="display:flex; gap:15px; align-items:flex-start;">
                <i class="fas fa-dna fa-2x" style="color:#6f42c1; margin-top:3px;"></i>
                <div>
                    <h4 style="margin:0 0 10px 0; color:#6f42c1; font-size:1.05em;">🧬 Acquired vs Intrinsic Resistance — Both Matter</h4>
                    <p style="margin:6px 0; font-size:0.93em; line-height:1.5;">The AMR story is <strong>more than acquired genes</strong>. We report <strong>both</strong> because they answer different questions:</p>
                    <ul style="margin:6px 0 10px 20px; font-size:0.93em; line-height:1.6; color:#333;">
                        <li><strong style="color:#7c3aed;">Intrinsic genes</strong> — part of the species' baseline genome. Present in almost all <em>S. aureus</em>. Examples: <em>norA</em>, <em>norB</em>, <em>mepA</em>, <em>lmrS</em> (efflux pumps). They tell you what the bacterium is <em>naturally tolerant to</em> and set the floor for susceptibility.</li>
                        <li><strong style="color:#e11d48;">Acquired genes</strong> — gained by horizontal transfer (plasmids, transposons, SCC elements, phages). Examples: <em>mecA</em>, <em>mecC</em>, <em>erm</em>, <em>tetK</em>, <em>dfrG</em>, <em>aacA-aphD</em>. They tell you what <em>selective pressure</em> the population has been under, and predict clinical failure of specific drugs.</li>
                    </ul>
                    <p style="margin:6px 0 0 0; font-size:0.92em; line-height:1.5;"><i class="fas fa-lightbulb" style="color:#6f42c1;"></i> <strong>Why both:</strong> A population of <em>S. aureus</em> with only intrinsic efflux genes will still respond to standard therapy; a population with <em>mecA + ermC + tetK</em> is multidrug-resistant and needs alternative agents. Reporting only acquired genes hides the intrinsic baseline; reporting only intrinsic genes misses the acquired threat.</p>
                </div>
            </div>
        </div>
        '''

        # ---- Genotype ≠ phenotype caveat ----
        caveat = '''
        <div class="alert-box alert-warning" style="border-left-color:#ffc107;">
            <i class="fas fa-exclamation-triangle fa-2x"></i>
            <div>
                <h4 style="margin:0 0 8px 0; color:#856404;">⚠️ Gene Presence ≠ Phenotypic Resistance</h4>
                <p style="margin:6px 0; font-size:0.93em; line-height:1.6;">Detecting an AMR gene in a genome is <strong>necessary evidence</strong>, but not sufficient to declare phenotypic resistance. Several biological mechanisms can break the link:</p>
                <ul style="margin:6px 0 12px 20px; font-size:0.92em; line-height:1.6;">
                    <li><strong>Silent / truncated genes</strong> — <em>mecA</em> inside a non-functional SCCmec cassette, or an <em>erm</em> gene with a premature stop codon, produces no resistance.</li>
                    <li><strong>Expression regulation</strong> — <em>mecA</em> is inducible by <em>mecI/mecR1</em>; a mutation in the regulator can silence resistance.</li>
                    <li><strong>Mechanism matters</strong> — <em>erm</em> (rRNA methylation) and <em>msrA</em> (efflux) both produce MLS<sub>B</sub> resistance but have different spectrum and inducibility.</li>
                    <li><strong>Naming ambiguity</strong> — ResFinder appends <code>_1</code> to primary alleles (e.g. <code>mecA_1</code>); some databases use older synonyms for the same gene.</li>
                    <li><strong>Dose and route</strong> — a low-level efflux pump may be overcome by higher drug concentrations <em>in vivo</em> even if it fails <em>in vitro</em>.</li>
                </ul>
                <p style="margin:8px 0 0 0; font-size:0.92em; line-height:1.6; background:#fff8e1; padding:8px 12px; border-radius:4px; border-left:3px solid #f59e0b;">
                    <strong>Clinical bottom line:</strong> <strong>Antimicrobial Susceptibility Testing (AST)</strong> on the cultured isolate remains the gold standard for treatment decisions. Genomic AMR prediction is a powerful triage tool and surveillance instrument — not a replacement for AST.
                </p>
            </div>
        </div>
        '''

        # ---- Filter buttons ----
        filter_buttons = self._filter_buttons('amr-table', [
            ('mecA (MRSA)', 'mecA', 'btn-danger', 'fa-skull-crossbones'),
            ('mecC', 'mecC', 'btn-danger', 'fa-skull-crossbones'),
            ('bla (β-lactamase)', 'bla', 'btn-danger', 'fa-skull-crossbones'),
            ('vanA', 'vanA', 'btn-warning', 'fa-biohazard'),
            ('vanB', 'vanB', 'btn-warning', 'fa-biohazard'),
            ('erm (MLS_B)', 'erm', 'btn-info', 'fa-pills'),
            ('msrA', 'msrA', 'btn-info', 'fa-pills'),
            ('mphC', 'mphC', 'btn-info', 'fa-pills'),
            ('tet', 'tet', 'btn-secondary', 'fa-capsules'),
            ('aac', 'aac', 'btn-secondary', 'fa-syringe'),
            ('aph', 'aph', 'btn-secondary', 'fa-syringe'),
            ('ant', 'ant', 'btn-secondary', 'fa-syringe'),
            ('dfr (Trimethoprim)', 'dfr', 'btn-light', 'fa-tablets'),
            ('cat (Chloramphenicol)', 'cat', 'btn-light', 'fa-tablets'),
            ('fos (Fosfomycin)', 'fos', 'btn-light', 'fa-tablets'),
            ('pco (Copper)', 'pco', 'btn-light', 'fa-tablets'),
        ])

        # ---- Gene family roles ----
        info_block = self._gene_family_info('#F44336',
            'Role of each AMR gene family in <em>S. aureus</em>:', [
                ('mecA/mecC', 'Methicillin resistance (MRSA) — confers resistance to all β‑lactam antibiotics.'),
                ('bla (blaZ)', 'β‑lactamase — penicillin resistance.'),
                ('vanA/vanB', 'Vancomycin resistance — the last‑line antibiotic for MRSA.'),
                ('erm', 'Macrolide, lincosamide, streptogramin B resistance (MLS_B).'),
                ('msrA', 'Macrolide efflux pump.'),
                ('mphC', 'Macrolide phosphotransferase.'),
                ('tet', 'Tetracycline resistance (efflux or ribosomal protection).'),
                ('aac/aph/ant', 'Aminoglycoside-modifying enzymes (gentamicin, kanamycin, tobramycin).'),
                ('dfr', 'Trimethoprim resistance.'),
                ('cat', 'Chloramphenicol resistance.'),
                ('fos', 'Fosfomycin resistance.'),
                ('pco', 'Copper resistance — often linked to metal tolerance.'),
            ])

        # ---- Database roles ----
        db_roles = '''
        <div class="database-section" style="margin-top:30px;">
            <h3 style="color:#2c3e50; border-bottom:2px solid #3b82f6; padding-bottom:10px;">
                <i class="fas fa-database"></i> Roles &amp; Strengths of Each Database
            </h3>
            <p style="color:#666; margin-bottom:15px;">
                Each database is optimised for a different purpose. Knowing what each one is <em>good at</em> helps you interpret single-database hits.
            </p>
            <div style="display:grid; grid-template-columns:repeat(auto-fit, minmax(300px, 1fr)); gap:20px; margin:20px 0;">

                <div style="background:#f8f9fa; padding:16px; border-radius:10px; border-left:4px solid #28a745;">
                    <h4 style="color:#28a745; margin:0 0 8px 0;">🟢 CARD</h4>
                    <p style="font-size:0.9em; color:#333; margin:0;"><strong>Strengths:</strong> Expert-curated, strict SNP-based cutoffs, precise allele assignments. Reference standard for resistance mechanism annotations.<br><strong>Weaknesses:</strong> Strict thresholds may miss novel or divergent variants.<br><strong>Best for:</strong> Confident allele-level calls; mechanism of action.</p>
                </div>

                <div style="background:#f8f9fa; padding:16px; border-radius:10px; border-left:4px solid #17a2b8;">
                    <h4 style="color:#17a2b8; margin:0 0 8px 0;">🔵 ResFinder</h4>
                    <p style="font-size:0.9em; color:#333; margin:0;"><strong>Strengths:</strong> Highly sensitive; frequently updated; excellent for known acquired determinants.<br><strong>Weaknesses:</strong> Appends <code>_1</code> to primary alleles (<code>mecA_1</code> = <code>mecA</code>); occasionally splits alleles into subvariants that confuse reading.<br><strong>Best for:</strong> Broad sensitivity; newly described variants.</p>
                </div>

                <div style="background:#f8f9fa; padding:16px; border-radius:10px; border-left:4px solid #007bff;">
                    <h4 style="color:#007bff; margin:0 0 8px 0;">🔷 NCBI AMR</h4>
                    <p style="font-size:0.9em; color:#333; margin:0;"><strong>Strengths:</strong> Curated by NCBI; tightly linked to the Reference Gene Catalog used by AMRFinderPlus; strong provenance.<br><strong>Weaknesses:</strong> Overlaps substantially with AMRFinderPlus output.<br><strong>Best for:</strong> Cross-checking AMRFinderPlus calls.</p>
                </div>

                <div style="background:#f8f9fa; padding:16px; border-radius:10px; border-left:4px solid #fd7e14;">
                    <h4 style="color:#fd7e14; margin:0 0 8px 0;">🟠 MEGARes</h4>
                    <p style="font-size:0.9em; color:#333; margin:0;"><strong>Strengths:</strong> Hierarchy of gene families; includes biocide and metal resistance; designed for metagenomics but works well on isolates.<br><strong>Weaknesses:</strong> Broader families can produce "class-level" hits that obscure the exact allele.<br><strong>Best for:</strong> Environmental co-selection markers; family-level grouping.</p>
                </div>

                <div style="background:#f8f9fa; padding:16px; border-radius:10px; border-left:4px solid #6f42c1;">
                    <h4 style="color:#6f42c1; margin:0 0 8px 0;">🟣 ARG-ANNOT</h4>
                    <p style="font-size:0.9em; color:#333; margin:0;"><strong>Strengths:</strong> Historic, well-curated ARG catalogue; complements newer databases.<br><strong>Weaknesses:</strong> Updated less frequently than CARD or ResFinder; some entries superseded.<br><strong>Best for:</strong> Historical comparisons; confirming classic determinants.</p>
                </div>

                <div style="background:#f8f9fa; padding:16px; border-radius:10px; border-left:4px solid #dc3545;">
                    <h4 style="color:#dc3545; margin:0 0 8px 0;">🔴 AMRFinderPlus</h4>
                    <p style="font-size:0.9em; color:#333; margin:0;"><strong>Strengths:</strong> NCBI gold standard; includes <strong>point mutations</strong> (gyrA, parC, rpoB, 23S rRNA) alongside gene presence; excellent validation behind each entry.<br><strong>Weaknesses:</strong> Conservative on some efflux families; may not flag intrinsic genes by design.<br><strong>Best for:</strong> Clinical-grade gene calls; point-mutation resistance; the most defensible single source.</p>
                </div>

            </div>

            <div class="alert-box alert-info" style="border-left-color:#00695c; background:#e8f5e9;">
                <i class="fas fa-link fa-2x" style="color:#00695c;"></i>
                <div>
                    <strong>🧭 You decide — we give you the evidence, not the verdict.</strong>
                    <p style="margin-top:8px; font-size:0.95em;">StaphScope deliberately <strong>does not merge or prioritise</strong> hits across databases. Different research questions call for different choices:</p>
                    <ul style="margin:8px 0 0 20px; font-size:0.93em;">
                        <li><strong>Want conservative clinical calls?</strong> Filter to AMRFinderPlus hits only — they are the most validated.</li>
                        <li><strong>Want maximum sensitivity?</strong> Keep ResFinder + MEGARes + CARD together.</li>
                        <li><strong>Doing surveillance?</strong> Report all databases; use cross-database agreement as your confidence metric.</li>
                        <li><strong>Comparing to older literature?</strong> Check ARG-ANNOT — the historic standard.</li>
                    </ul>
                    <p style="margin-top:10px; font-size:0.92em; background:#fff3cd; padding:8px 14px; border-radius:4px; border-left:3px solid #ffc107;">
                        <i class="fas fa-lightbulb"></i> <strong>All hits are visible in the table above.</strong> Use the search box to filter by database (type <code>CARD</code>, <code>ResFinder</code>, etc.). Use the grouping dropdown to see which clonal backgrounds carry each hit. Export to CSV to merge or filter however your analysis requires — the raw provenance is preserved for you.
                    </p>
                </div>
            </div>
        </div>
        '''

        return f'''
        {credit}
        {self._alert('info', 'fa-biohazard',
            '<h3>🧬 AMR Genes – Gene-Centric View + Grouping</h3>'
            '<p>Each gene is shown with all genomes carrying it. Use the dropdown to reveal clone associations.</p>')}
        {why_multi_db}
        {confidence_box}
        {acquired_intrinsic}
        {caveat}
        {filter_buttons}
        {info_block}
        {self._gene_table('amr-table', 'All AMR Genes Across Databases',
                          all_genes, total_samples, self.analyzer.critical_amr_genes)}
        {self._database_cards(amr_dbs, {'amrfinder': 'AMRfinder'})}
        {db_roles}'''

    # =========================================================================
    # VIRULENCE
    # =========================================================================
    def _sec_virulence(self, d):
        """Virulence gene-centric tab with filter buttons, info box, and DB cards."""
        gene_centric = d['gene_centric']
        vir_dbs = gene_centric.get('virulence_databases', {})
        total_samples = len(d['samples'])
        all_genes = []
        for genes in vir_dbs.values():
            all_genes.extend(genes)
        all_genes.sort(key=lambda x: x['count'], reverse=True)

        credit = self._credit_bar('#E91E63', '🦠',
            'Virulence Detection Tools &amp; Databases',
            '<strong>ABRicate</strong> by '
            '<a href="https://github.com/tseemann/abricate" target="_blank" style="color:#E91E63;font-weight:bold;">Prof. Torsten Seemann</a> '
            '→ Powered by the '
            '<a href="http://www.mgc.ac.cn/VFs/" target="_blank" style="color:#E91E63;font-weight:bold;">VFDB (Virulence Factor Database)</a> '
            '— an open‑source resource for bacterial virulence factors.<br>'
            '<span style="font-size:0.9em;color:#6c757d;"><i class="fas fa-gratipay"></i> '
            'We thank the VFDB curators and Prof. Seemann for their invaluable contributions.</span>')

        filter_buttons = self._filter_buttons('vir-table', [
            ('PVL (luk)', 'luk', 'btn-danger', 'fa-skull'),
            ('TSST-1', 'tsst', 'btn-danger', 'fa-biohazard'),
            ('sea', 'sea', 'btn-warning', 'fa-filter'),
            ('seb', 'seb', 'btn-warning', 'fa-filter'),
            ('sec', 'sec', 'btn-warning', 'fa-filter'),
            ('sed', 'sed', 'btn-warning', 'fa-filter'),
            ('see', 'see', 'btn-warning', 'fa-filter'),
            ('seg', 'seg', 'btn-warning', 'fa-filter'),
            ('seh', 'seh', 'btn-warning', 'fa-filter'),
            ('sei', 'sei', 'btn-warning', 'fa-filter'),
            ('sej', 'sej', 'btn-warning', 'fa-filter'),
            ('sek', 'sek', 'btn-warning', 'fa-filter'),
            ('sel', 'sel', 'btn-warning', 'fa-filter'),
            ('eta', 'eta', 'btn-warning', 'fa-filter'),
            ('etb', 'etb', 'btn-warning', 'fa-filter'),
            ('hla', 'hla', 'btn-info', 'fa-filter'),
            ('hlb', 'hlb', 'btn-info', 'fa-filter'),
            ('hlg', 'hlg', 'btn-info', 'fa-filter'),
            ('hld', 'hld', 'btn-info', 'fa-filter'),
            ('ica (Biofilm)', 'ica', 'btn-secondary', 'fa-layer-group'),
            ('scn', 'scn', 'btn-secondary', 'fa-filter'),
            ('eap', 'eap', 'btn-secondary', 'fa-filter'),
            ('ads', 'ads', 'btn-secondary', 'fa-filter'),
            ('clfA/B', 'clf', 'btn-secondary', 'fa-filter'),
            ('fnbA/B', 'fnb', 'btn-secondary', 'fa-filter'),
            ('sdr (SdrC/D/E)', 'sdr', 'btn-secondary', 'fa-filter'),
            ('cap (Capsule)', 'cap', 'btn-secondary', 'fa-filter'),
            ('esa', 'esa', 'btn-secondary', 'fa-filter'),
            ('ess', 'ess', 'btn-secondary', 'fa-filter'),
            ('esx', 'esx', 'btn-secondary', 'fa-filter'),
            ('set (exotoxin-like)', 'set', 'btn-secondary', 'fa-filter'),
        ])

        info_block = self._gene_family_info('#E91E63',
            'Role of each virulence factor in <em>S. aureus</em>:', [
                ('PVL (lukF/S-PV)', 'Panton-Valentine leukocidin — leukocyte destruction; severe skin/soft tissue infection, necrotizing pneumonia.'),
                ('TSST-1 (tsst)', 'Toxic shock syndrome toxin-1 — causes TSS.'),
                ('Enterotoxins (sea–see, seg–seu)', 'Superantigens — food poisoning and toxic shock.'),
                ('Exfoliative toxins (eta, etb)', 'Staphylococcal scalded skin syndrome (SSSS).'),
                ('Hemolysins (hla, hlb, hlg, hld)', 'Lyse red blood cells; alpha-toxin (hla) is a major virulence factor.'),
                ('Biofilm (icaADBC)', 'Polysaccharide intercellular adhesin — device-related and chronic infections.'),
                ('Immune evasion (scn, eap, ads)', 'Inhibit complement, neutrophil chemotaxis, adenosine signalling.'),
                ('Adhesins (clfA/B, fnbA/B, sdrC/D/E)', 'Bind fibrinogen and fibronectin — promote attachment to host tissues and medical devices.'),
                ('Capsule (cap)', 'Polysaccharide capsule — protects against phagocytosis.'),
                ('ESAT-6 secretion (esa, ess, esx)', 'Type VII secretion system components — virulence and immune modulation.'),
                ('SET exotoxin-like proteins', 'Superantigen-like toxins modulating host immune responses.'),
            ])

        return f'''
        {credit}
        {self._alert('info', 'fa-virus',
            '<h3>🧬 Virulence Factors – Gene-Centric View</h3>'
            '<p>Use the grouping dropdown to see which clones carry specific virulence factors.</p>')}
        {filter_buttons}
        {info_block}
        {self._gene_table('vir-table', 'All Virulence Genes',
                          all_genes, total_samples, self.analyzer.critical_virulence_genes)}
        {self._database_cards(vir_dbs)}'''

    # =========================================================================
    # BACMET
    # =========================================================================
    def _sec_bacmet(self, d):
        """BACMET gene-centric tab with filter buttons, info box, and DB cards."""
        gene_centric = d['gene_centric']
        bacmet_dbs = gene_centric.get('bacmet_databases', {})
        total_samples = len(d['samples'])
        all_genes = []
        for genes in bacmet_dbs.values():
            all_genes.extend(genes)
        all_genes.sort(key=lambda x: x['count'], reverse=True)
        if not all_genes:
            return self._alert('warning', 'fa-flask',
                '<h3>No BACMET Data Available</h3>'
                '<p>Biocide and heavy metal resistance genes were not detected.</p>')

        credit = self._credit_bar('#FF5722', '🧪',
            'BACMET – Biocide &amp; Heavy Metal Resistance',
            '<strong>ABRicate</strong> by '
            '<a href="https://github.com/tseemann/abricate" target="_blank" style="color:#FF5722;font-weight:bold;">Prof. Torsten Seemann</a> '
            '→ Powered by the '
            '<a href="https://bacmet.biomedicine.gu.se/" target="_blank" style="color:#FF5722;font-weight:bold;">BacMet database</a> '
            '— an open‑source resource for antibacterial biocide and metal resistance genes.<br>'
            '<span style="font-size:0.9em;color:#6c757d;"><i class="fas fa-gratipay"></i> '
            'We thank the BacMet curators and Prof. Seemann for their invaluable contributions.</span>')

        filter_buttons = self._filter_buttons('bac-table', [
            ('qac (Quat. ammonium)', 'qac', 'btn-info', 'fa-filter'),
            ('cep (Chlorhexidine)', 'cep', 'btn-info', 'fa-filter'),
            ('form (Formaldehyde)', 'form', 'btn-info', 'fa-filter'),
            ('mer (Mercury)', 'mer', 'btn-warning', 'fa-filter'),
            ('ars (Arsenic)', 'ars', 'btn-warning', 'fa-filter'),
            ('cop (Copper)', 'cop', 'btn-warning', 'fa-filter'),
            ('sil (Silver)', 'sil', 'btn-warning', 'fa-filter'),
            ('cad (Cadmium)', 'cad', 'btn-warning', 'fa-filter'),
            ('znt (Zinc)', 'znt', 'btn-warning', 'fa-filter'),
            ('czc (Co-Zn-Cd)', 'czc', 'btn-warning', 'fa-filter'),
            ('chr (Chromate)', 'chr', 'btn-warning', 'fa-filter'),
            ('pbr (Lead)', 'pbr', 'btn-warning', 'fa-filter'),
            ('nik (Nickel)', 'nik', 'btn-warning', 'fa-filter'),
            ('soxR', 'soxR', 'btn-secondary', 'fa-filter'),
            ('cpxR', 'cpxR', 'btn-secondary', 'fa-filter'),
            ('baeR', 'baeR', 'btn-secondary', 'fa-filter'),
            ('emr', 'emr', 'btn-secondary', 'fa-filter'),
            ('norA', 'norA', 'btn-secondary', 'fa-filter'),
        ])

        info_block = self._gene_family_info('#FF5722',
            'Environmental co-selection — why these genes matter:', [
                ('qac family', 'Quaternary ammonium compounds — hospital disinfectants.'),
                ('cep', 'Chlorhexidine resistance (antiseptic).'),
                ('form', 'Formaldehyde resistance.'),
                ('mer', 'Mercury resistance — often on transposons and plasmids.'),
                ('ars', 'Arsenic resistance — common in environmental and clinical isolates.'),
                ('cop / sil', 'Copper and silver resistance — linked to metal-based antimicrobials.'),
                ('czc / cad / znt', 'Zinc, cadmium, cobalt efflux — co-selection with antibiotic resistance.'),
                ('chr', 'Chromate resistance.'),
                ('pbr', 'Lead resistance.'),
                ('nik', 'Nickel transport.'),
                ('soxR / cpxR / baeR', 'Stress response regulators that also upregulate multidrug efflux pumps.'),
                ('emr / sme / norA / mdeA', 'Multidrug efflux pumps that export biocides and antibiotics.'),
            ])

        return f'''
        {credit}
        {self._alert('info', 'fa-flask',
            '<h3>🧪 BACMET – Biocide / Heavy Metal Resistance</h3>'
            '<p>Environmental co-selection markers. Use the dropdown to see which clones carry them.</p>')}
        {filter_buttons}
        {info_block}
        {self._gene_table('bac-table', 'All BACMET Genes', all_genes, total_samples)}
        {self._database_cards(bacmet_dbs)}'''

    # =========================================================================
    # PLASMIDS
    # =========================================================================
    def _sec_plasmids(self, d):
        """Plasmid replicon tab with filter buttons, info box, and DB cards."""
        gene_centric = d['gene_centric']
        plas_dbs = gene_centric.get('plasmid_databases', {})
        total_samples = len(d['samples'])
        all_genes = []
        for genes in plas_dbs.values():
            all_genes.extend(genes)
        all_genes.sort(key=lambda x: x['count'], reverse=True)

        credit = self._credit_bar('#673AB7', '🔌',
            'Plasmid Replicon Typing',
            '<strong>ABRicate</strong> by '
            '<a href="https://github.com/tseemann/abricate" target="_blank" style="color:#673AB7;font-weight:bold;">Prof. Torsten Seemann</a> '
            '→ Powered by '
            '<a href="https://genepi.food.dtu.dk/PlasmidFinder/" target="_blank" style="color:#673AB7;font-weight:bold;">PlasmidFinder</a> '
            'from the '
            '<a href="https://genepi.food.dtu.dk/" target="_blank" style="color:#673AB7;font-weight:bold;">Center for Genomic Epidemiology (DTU)</a>.<br>'
            '<span style="font-size:0.9em;color:#6c757d;"><i class="fas fa-gratipay"></i> '
            'We thank the CGE team and Prof. Seemann for making these open‑source resources available.</span>')

        filter_buttons = self._filter_buttons('plasmid-table', [
            ('rep (Rolling-circle)', 'rep', 'btn-info', 'fa-filter'),
            ('pS (pS194 family)', 'pS', 'btn-info', 'fa-filter'),
            ('pT (pT181 family)', 'pT', 'btn-info', 'fa-filter'),
            ('pC (pC194 family)', 'pC', 'btn-info', 'fa-filter'),
            ('Inc (Inc groups)', 'Inc', 'btn-warning', 'fa-filter'),
            ('pUB (pUB110)', 'pUB', 'btn-warning', 'fa-filter'),
            ('pE (pE194)', 'pE', 'btn-warning', 'fa-filter'),
            ('pI (pI258)', 'pI', 'btn-warning', 'fa-filter'),
            ('pBL (pBL1)', 'pBL', 'btn-secondary', 'fa-filter'),
            ('pIP (pIP501)', 'pIP', 'btn-secondary', 'fa-filter'),
        ])

        info_block = self._gene_family_info('#673AB7',
            'Plasmid replicon families in <em>S. aureus</em>:', [
                ('Rolling-circle (rep)', 'Small, high-copy plasmids — antibiotic resistance.'),
                ('pS194 family (pS)', 'Small plasmids, often carry resistance determinants.'),
                ('pT181 family (pT)', 'Tetracycline resistance plasmids (e.g., tetK).'),
                ('pC194 family (pC)', 'Chloramphenicol resistance plasmids.'),
                ('pUB110', 'Kanamycin/neomycin resistance — widely used in S. aureus.'),
                ('pE194', 'Erythromycin resistance (ermC).'),
                ('pI258', 'Mercury resistance plasmid.'),
                ('pBL1', 'β-lactamase plasmid.'),
                ('pIP501', 'Multidrug resistance plasmid (streptomycin, chloramphenicol, erythromycin).'),
                ('Inc groups', 'Incompatibility groups common in Gram-negative plasmids.'),
            ])

        return f'''
        {credit}
        {self._alert('info', 'fa-plug',
            '<h3>🧬 Plasmid Replicons – Horizontal Gene Transfer Markers</h3>'
            '<p>Replicon families reveal which clones carry which plasmids.</p>')}
        {filter_buttons}
        {info_block}
        {self._gene_table('plasmid-table', 'Plasmid Replicons', all_genes, total_samples)}
        {self._database_cards(plas_dbs)}'''

    # =========================================================================
    # MUTATIONS
    # =========================================================================
    def _sec_mutation(self, d):
        """Point mutation tab with filter buttons, info box, and grouping."""
        mut = d.get('mutation_data', {})
        mutations = mut.get('mutations', [])
        total_samples = len(d['samples'])
        if not mutations:
            return self._alert('warning', 'fa-dna',
                '<h3>No Mutation Data Available</h3>'
                '<p>The mutation_summary.html file was not found.</p>')

        credit = self._credit_bar('#00BCD4', '🧬',
            'AMRFinderPlus – Point Mutations',
            '<strong>AMRFinderPlus</strong> by '
            '<a href="https://github.com/ncbi/amr" target="_blank" style="color:#00BCD4;font-weight:bold;">NCBI</a> '
            '— the most comprehensive database for antimicrobial resistance genes and point mutations.<br>'
            '<span style="font-size:0.9em;color:#6c757d;"><i class="fas fa-gratipay"></i> '
            'We thank the NCBI team for maintaining this essential open-resource for genomic surveillance.</span>')

        filter_buttons = self._filter_buttons('mutation-table', [
            ('Linezolid', 'LINEZOLID', 'btn-danger', 'fa-skull-crossbones'),
            ('Quinolone', 'QUINOLONE', 'btn-warning', 'fa-biohazard'),
            ('Rifampin', 'RIFAMPIN', 'btn-warning', 'fa-biohazard'),
            ('gyrA', 'gyrA', 'btn-info', 'fa-filter'),
            ('parC', 'parC', 'btn-info', 'fa-filter'),
            ('rpoB', 'rpoB', 'btn-info', 'fa-filter'),
            ('mprF (Daptomycin)', 'mprF', 'btn-info', 'fa-filter'),
            ('rplC (Linezolid)', 'rplC', 'btn-info', 'fa-filter'),
            ('rplD (Linezolid)', 'rplD', 'btn-info', 'fa-filter'),
            ('23S rRNA', '23S', 'btn-info', 'fa-filter'),
        ])

        info_block = self._gene_family_info('#00BCD4',
            'Clinical relevance of key mutations:', [
                ('23S rRNA (linezolid)', 'Mutations (e.g., G2576T, T2500A) confer linezolid resistance — last-line antibiotic for MRSA.'),
                ('gyrA / parC (quinolones)', 'QRDR mutations reduce susceptibility to fluoroquinolones (ciprofloxacin, levofloxacin).'),
                ('rpoB (rifampin)', 'High-level rifampin resistance — combination therapy consideration.'),
                ('mprF (daptomycin)', 'Daptomycin non-susceptibility.'),
                ('rplC / rplD (linezolid)', 'Ribosomal protein mutations also confer linezolid resistance.'),
                ('fusA (fusidic acid)', 'Fusidic acid resistance.'),
                ('mupA (mupirocin)', 'High-level mupirocin resistance.'),
            ])

        rows = ''
        for m in mutations:
            tags = ''.join(f'<span class="genome-tag">{esc(g)}</span>' for g in m['genomes'])
            rows += (f'<tr><td><strong>{esc(m["gene"])}</strong></td>'
                     f'<td>{esc(m["mutation"])}</td>'
                     f'<td>{esc(m["class"])}</td>'
                     f'<td>{esc(m["subclass"])}</td>'
                     f'<td><strong>{len(m["genomes"])}</strong> '
                     f'({len(m["genomes"]) / total_samples * 100:.1f}%)</td>'
                     f'<td><div class="genome-list">{tags}</div></td></tr>')
        return f'''
        {credit}
        {self._alert('info', 'fa-dna',
            '<h3>🧬 Point Mutations – Gene-Centric View</h3>'
            '<p>Mutations in gyrA, parC, rpoB, 23S rRNA, mprF, rplC, rplD confer resistance to key antibiotics.</p>')}
        {filter_buttons}
        {info_block}
        <h3>🔗 Group Genomes by Typing</h3>
        {self._grouping_dropdown('mutation-table')}
        <input type="text" class="search-box" id="search-mutation-table"
               onkeyup="searchTable('mutation-table','search-mutation-table')"
               placeholder="🔍 Search mutation by gene or mutation name...">
        <div class="master-scrollable-container">
            <table id="mutation-table" class="data-table">
                <thead><tr>
                    <th data-sort="string">Gene</th>
                    <th data-sort="string">Mutation</th>
                    <th data-sort="string">Class</th>
                    <th data-sort="string">Subclass</th>
                    <th data-sort="number">Count</th>
                    <th data-sort="string">Genomes (scrollable, groupable)</th>
                </tr></thead>
                <tbody>{rows}</tbody>
            </table>
        </div>'''

    # =========================================================================
    # MGE
    # =========================================================================
    def _sec_mge(self, d):
        """Mobile genetic element profile tab with educational context."""
        mge = d.get('mge_data', {})
        stats = mge.get('stats', {})
        per_sample = mge.get('per_sample', {})
        analysis = mge.get('analysis', {})
        note = mge.get('note', '')

        credit = self._credit_bar('#16A085', '🧬',
            'Mobile Genetic Elements – mobileOG-db',
            '<strong>mobileOG-db</strong> developed by '
            '<a href="https://github.com/clb21565/mobileOG-db" target="_blank" '
            'style="color:#16A085;font-weight:bold;">Brown, Mullet and colleagues</a> '
            '— a manually curated database of protein families mediating the life cycle of '
            'bacterial mobile genetic elements.<br>'
            '<span style="font-size:0.9em;color:#6c757d;"><i class="fas fa-book-open"></i> '
            'Please cite: Brown CL, Mullet J, Hindi F, Stoll JE, Gupta S, Choi M, Keenum I, '
            'Vikesland P, Pruden A, Zhang L. mobileOG-db: a Manually Curated Database of Protein '
            'Families Mediating the Life Cycle of Bacterial Mobile Genetic Elements. '
            '<em>Appl Environ Microbiol</em>. 2022;88(18):e00991-22. '
            '<a href="https://doi.org/10.1128/aem.00991-22" target="_blank" '
            'style="color:#16A085;font-weight:bold;">🔗 DOI</a></span><br>'
            '<span style="font-size:0.9em;color:#6c757d;"><i class="fas fa-gratipay"></i> '
            'We are grateful to the mobileOG-db curators and contributors for making this resource freely available.</span>')

        why_mge = '''
        <div class="alert-box alert-info" style="border-left-color:#16A085;">
            <i class="fas fa-dna fa-2x" style="color:#16A085;"></i>
            <div>
                <h3>🧬 Why Mobile Genetic Elements Matter for AMR</h3>
                <p>Mobile genetic elements (MGEs) are the <strong>physical vehicles of horizontal gene transfer</strong> — the mechanism by which bacteria acquire resistance and virulence genes from other bacteria. Without MGEs, most acquired resistance in <em>S. aureus</em> would not exist.</p>
                <ul style="margin-top:8px;">
                    <li><strong>Plasmids</strong> carry <em>blaZ</em>, <em>ermC</em>, <em>tetK</em>, and multidrug-resistance cassettes between lineages.</li>
                    <li><strong>SCCmec elements</strong> (a type of ICE) carry <em>mecA</em> — the defining marker of MRSA.</li>
                    <li><strong>Prophages</strong> carry virulence factors like PVL (<em>lukF/S-PV</em>), TSST-1 (<em>tsst</em>), and staphylokinase (<em>scn</em>).</li>
                    <li><strong>IS elements and transposons</strong> move resistance genes into new genomic contexts and can activate or disrupt neighbouring genes.</li>
                </ul>
                <p style="margin-top:8px;"><strong>Interpreting this tab:</strong> higher mobileOG hit counts suggest a genome that has accumulated more mobile-element machinery — potentially more permissive for acquiring new resistance genes.</p>
            </div>
        </div>'''

        category_cards = self._mge_category_cards()

        if not per_sample:
            return f'''
            {credit}
            {why_mge}
            {category_cards}
            {self._alert('warning', 'fa-mobile-alt',
                '<h3>No MGE Data Available</h3>'
                '<p>The MGE summary HTML was not found or could not be parsed.</p>')}'''

        cards = ''
        cards += self._stat_card(stats.get('Total Samples', len(per_sample)),
                                 'Total Samples', '#4CAF50', 'fa-vial')
        cards += self._stat_card(stats.get('With mobileOG hits', '—'),
                                 'With mobileOG hits', '#2196F3', 'fa-check')
        cards += self._stat_card(stats.get('Without mobileOG hits', '—'),
                                 'Without mobileOG hits', '#9E9E9E', 'fa-times')
        cards += self._stat_card(stats.get('Total mobileOG hits', analysis.get('grand_total', '—')),
                                 'Total mobileOG hits', '#FF9800', 'fa-database')
        cards += self._stat_card(stats.get('Key MGE-associated genes', '—'),
                                 'Key MGE-associated genes', '#8B5CF6', 'fa-star')
        cards += self._stat_card(stats.get('Runtime', '—'),
                                 'Runtime', '#E53935', 'fa-clock')

        cat_rows = ''
        for c in analysis.get('category_totals', []):
            cat_rows += (f'<tr><td><strong>{esc(c["category"])}</strong></td>'
                         f'<td>{c["total"]:,}</td>'
                         f'<td>{c["pct_of_all"]}%</td>'
                         f'<td>{c["mean"]}</td>'
                         f'<td>{c["max"]}</td>'
                         f'<td>{esc(c["max_sample"])}</td></tr>')

        all_cols = ['Sample']
        seen = set()
        for s in per_sample.values():
            for k in s.keys():
                if k not in seen:
                    seen.add(k)
                    all_cols.append(k)
        sample_header = ''.join(
            f'<th data-sort="{"string" if c == "Sample" else "number"}">{esc(c)}</th>'
            for c in all_cols)
        sample_rows = ''
        for sample, vals in sorted(per_sample.items()):
            cells = f'<td><strong>{esc(sample)}</strong></td>'
            for c in all_cols[1:]:
                cells += f'<td>{vals.get(c, 0)}</td>'
            sample_rows += f'<tr>{cells}</tr>'

        agg_html = self._build_mge_aggregate(per_sample)

        top_html = ''
        for col, label in [('Key MGE-assoc.', 'Key MGE-associated'),
                           ('IS-assoc.', 'IS-associated'),
                           ('Plasmid-assoc.', 'Plasmid-associated')]:
            top = sorted(per_sample.items(),
                         key=lambda x: x[1].get(col, 0), reverse=True)[:5]
            rows = ''.join(f'<tr><td><strong>{esc(s)}</strong></td>'
                           f'<td>{v.get(col, 0)}</td></tr>' for s, v in top)
            top_html += f'''<h4>Top 5 — {label}</h4>
            <table class="data-table">
                <thead><tr><th>Sample</th><th>{label} hits</th></tr></thead>
                <tbody>{rows}</tbody></table>'''

        note_html = ''
        if note:
            note_html = (f'<div class="alert-box alert-warning" style="font-size:.9em;">'
                         f'<strong>Interpretation note:</strong> {note}</div>')

        return f'''
        {credit}
        {self._alert('info', 'fa-mobile-alt',
            '<h3>🧬 Mobile Genetic Elements – mobileOG Profile</h3>'
            '<p>Per-sample counts of mobileOG-db protein family hits across multiple functional categories. '
            'mobileOG-db identifies protein families associated with mobile genetic element biology. '
            'Individual matches represent MGE-associated protein signatures and should not be interpreted '
            'as independent mobile genetic elements — multiple proteins from the same element may generate '
            'multiple hits, and some families (particularly replication, recombination, and repair proteins) '
            'also occur in the bacterial chromosome. Stronger element-level inference requires genomic context '
            'and co-localization of compatible MGE-associated genes.</p>')}
        <div class="stats-grid">{cards}</div>
        {note_html}
        {why_mge}

        <h3 style="margin-top:30px;">📚 Understanding the Categories</h3>
        <p style="color:#666;font-size:.92em;margin-bottom:10px;">
            Each column in the tables below corresponds to one functional category of the mobileOG-db.
            Hover any card to learn what the category means and why it matters clinically.
        </p>
        {category_cards}

        <h3 style="margin-top:30px;">📊 Category Totals (aggregated across all samples)</h3>
        <div class="scrollable-table">
            <table class="data-table">
                <thead><tr>
                    <th data-sort="string">Category</th>
                    <th data-sort="number">Total hits</th>
                    <th data-sort="number">% of all hits</th>
                    <th data-sort="number">Mean/sample</th>
                    <th data-sort="number">Max/sample</th>
                    <th data-sort="string">Sample with max</th>
                </tr></thead>
                <tbody>{cat_rows}</tbody>
            </table>
        </div>

        <h3 style="margin-top:30px;">🧬 Per-sample mobileOG Profile</h3>
        <input type="text" class="search-box" id="search-mge"
               onkeyup="searchTable('mge-table','search-mge')"
               placeholder="🔍 Search sample...">
        <div class="action-buttons">
            <button class="action-btn btn-primary" onclick="exportTableToCSV('mge-table', 'mge_profile.csv')">
                <i class="fas fa-download"></i> Export MGE CSV</button>
        </div>
        <div class="master-scrollable-container">
            <table id="mge-table" class="data-table">
                <thead><tr>{sample_header}</tr></thead>
                <tbody>{sample_rows}</tbody>
            </table>
        </div>

        <h3 style="margin-top:30px;">📈 Aggregate MGE Profile by Typing</h3>
        {agg_html}

        <h3 style="margin-top:30px;">🏆 Top-N Views</h3>
        <div style="display:grid;grid-template-columns:repeat(auto-fit,minmax(280px,1fr));gap:20px;">
            {top_html}
        </div>'''

    def _mge_category_cards(self) -> str:
        """Educational cards describing each mobileOG-db category."""
        cards = [
            ('Integr.', 'Integrase', '#ef4444', 'fa-cut',
             'Enzymes that catalyse the integration and excision of mobile elements into and out of the host chromosome. Includes serine and tyrosine recombinases.',
             'High integrase counts suggest active chromosomal integration of foreign DNA.'),
            ('Transfer', 'Transfer / Conjugation', '#f97316', 'fa-share-alt',
             'Type IV secretion systems (T4SS), relaxases, and coupling proteins that mediate the horizontal transfer of DNA between bacteria.',
             'The molecular machinery of horizontal gene transfer — the primary route by which resistance spreads.'),
            ('Stab.', 'Stability / Maintenance', '#eab308', 'fa-anchor',
             'Partitioning systems (parAB), toxin-antitoxin cassettes, and addiction modules that ensure plasmids and other MGEs are faithfully inherited.',
             'A plasmid with strong stability systems persists in a population even without antibiotic selection — a major AMR concern.'),
            ('Phage', 'Bacteriophage', '#22c55e', 'fa-virus',
             'Complete phage-related proteins including structural capsids, tails, and lysogeny modules. Bacteriophages are viruses that infect bacteria.',
             'Phages can carry virulence genes (e.g. PVL, TSST-1) and can mobilise chromosomal DNA via transduction.'),
            ('Repl.', 'Replication', '#14b8a6', 'fa-sync',
             'Replication initiators, primases, helicases, and single-stranded binding proteins that drive autonomous replication of plasmids and other MGEs.',
             'Only elements with a functional replication module can stably persist as extrachromosomal DNA.'),
            ('IS-assoc.', 'Insertion Sequences', '#0891b2', 'fa-exchange-alt',
             'Transposases and IS elements — the smallest autonomous MGEs. They move within and between genomes by a cut-and-paste or copy-and-paste mechanism.',
             'IS elements frequently disrupt genes or carry resistance genes into new contexts.'),
            ('ICE-assoc.', 'Integrative Conjugative Elements', '#3b82f6', 'fa-project-diagram',
             'Self-transmissible elements that integrate into the chromosome and can excise, circularise, and conjugate to a new host — including SCCmec-like elements.',
             'ICEs move large multi-gene cassettes between unrelated bacteria — a major vehicle for SCCmec and resistance islands.'),
            ('Plasmid-assoc.', 'Plasmid', '#8b5cf6', 'fa-plug',
             'Plasmid backbone, replication, partitioning, and conjugal transfer proteins. Plasmids are extrachromosomal, usually circular, self-replicating DNA molecules.',
             'Plasmids are the single most important vehicle for acquired antibiotic resistance in S. aureus — including blaZ, ermC, and tetK.'),
            ('Phage-assoc.', 'Phage-associated (lysogeny)', '#d946ef', 'fa-dna',
             'Prophage and lysogeny-associated proteins distinct from fully structured phage particles — including prophage integrases, repressors, and lysis modules.',
             'Prophage carriage often correlates with virulence gene acquisition (e.g. PVL prophage in ST121, ST152).'),
            ('Key MGE-assoc.', 'Key MGE signatures', '#dc2626', 'fa-star',
             'A curated subset of high-confidence mobileOG-db entries selected as strong indicators of mobile element biology.',
             'The cleanest single indicator of MGE activity in a genome — useful when comparing across large cohorts.'),
        ]
        html = ('<div style="display:grid;'
                'grid-template-columns:repeat(auto-fit,minmax(340px,1fr));'
                'gap:15px;margin:20px 0;">')
        for key, label, color, icon, definition, importance in cards:
            html += f'''
            <div style="background:#fff;padding:16px;border-radius:10px;
                        border-left:4px solid {color};
                        box-shadow:0 2px 8px rgba(0,0,0,.06);">
                <div style="display:flex;align-items:center;gap:10px;margin-bottom:8px;">
                    <i class="fas {icon}" style="color:{color};font-size:1.4em;"></i>
                    <div>
                        <strong style="color:{color};font-size:1.05em;">{label}</strong>
                        <div style="font-size:0.8em;color:#888;font-family:monospace;">{key}</div>
                    </div>
                </div>
                <p style="font-size:0.88em;color:#333;margin:0 0 8px 0;">{definition}</p>
                <p style="font-size:0.85em;color:#555;margin:0;padding-top:8px;
                          border-top:1px dashed #e0e0e0;">
                    <i class="fas fa-bullseye" style="color:{color};"></i>
                    <em>{importance}</em>
                </p>
            </div>'''
        html += '</div>'
        return html

    def _build_mge_aggregate(self, per_sample):
        """Server-side aggregation of MGE profiles by each typing field."""
        samples_data = getattr(self, '_current_samples_data', {})
        if not samples_data or not per_sample:
            return ''
        fields = [
            ('MLST', 'MLST'),
            ('spa_Type', 'spa'),
            ('agr_Type', 'agr'),
            ('capsule_type', 'Capsule'),
            ('SCCmec_CGE', 'SCCmec (CGE)'),
            ('SCCmec_RPet', 'SCCmec (RPet)'),
            ('SCCmec_Subtype', 'SCCmec Subtype'),
        ]
        sample_cols = list(next(iter(per_sample.values())).keys())

        blocks = ''
        for field_key, field_label in fields:
            groups = defaultdict(list)
            for sample, vals in per_sample.items():
                sd = samples_data.get(sample, {})
                typ = sd.get('typing', {}).get(field_key, 'Not Assigned')
                if typ in ('Not Assigned', 'ND', ''):
                    continue
                groups[typ].append(vals)
            if not groups:
                continue
            rows = ''
            for grp, records in sorted(groups.items(), key=lambda x: -len(x[1])):
                means = {}
                for c in sample_cols:
                    vals_ = [r.get(c, 0) for r in records
                             if isinstance(r.get(c, 0), (int, float))]
                    means[c] = round(sum(vals_) / len(vals_), 1) if vals_ else 0
                cells = ''.join(f'<td>{means[c]}</td>' for c in sample_cols)
                rows += (f'<tr><td><strong>{esc(grp)}</strong></td>'
                         f'<td>{len(records)}</td>{cells}</tr>')
            header = ''.join(f'<th>{esc(c)}</th>' for c in sample_cols)
            blocks += f'''
            <h4 style="margin-top:20px;">Mean MGE profile by {field_label}</h4>
            <div class="master-scrollable-container">
                <table class="data-table">
                    <thead><tr>
                        <th data-sort="string">{field_label}</th>
                        <th data-sort="number">N samples</th>
                        {header}
                    </tr></thead>
                    <tbody>{rows}</tbody>
                </table>
            </div>'''
        return blocks

    # =========================================================================
    # PATTERNS
    # =========================================================================
    def _sec_patterns(self, d):
        """Cross-genome pattern discovery: triple, four-way, high-risk, co-occurrence."""
        P = d['patterns']
        html = self._alert('info', 'fa-project-diagram',
            '<h3>🔍 Cross-Genome Pattern Discovery</h3>'
            '<p>Triple/four-way typing, gene co-occurrence, and high-risk combinations.</p>')
        html += self._combo_table('triple-table', 'Triple Typing (ST – spa – SCCmec CGE)',
                                  P.get('mlst_spa_sccmec_cge', {}))
        html += self._combo_table('four-way-table', 'Four-Way Typing (ST – spa – SCCmec CGE – agr)',
                                  P.get('mlst_spa_sccmec_agr', {}))
        high_risk = P.get('high_risk_combinations', [])
        if high_risk:
            rows = ''
            for c in high_risk:
                rows += (f'<tr><td><strong>{esc(c["sample"])}</strong></td>'
                         f'<td>{esc(c["mlst"])}</td><td>{esc(c["spa_type"])}</td>'
                         f'<td>{esc(c["sccmec_type"])}</td><td>{esc(c["agr_type"])}</td>'
                         f'<td>{", ".join(c["critical_amr_genes"])}</td>'
                         f'<td>{", ".join(c["critical_virulence_genes"])}</td></tr>')
            html += f'''<h3>⚠️ High-Risk Combinations</h3>
            <div class="master-scrollable-container"><table class="data-table">
                <thead><tr><th>Sample</th><th>MLST</th><th>spa</th><th>SCCmec</th><th>agr</th>
                <th>Critical AMR</th><th>Critical Virulence</th></tr></thead>
                <tbody>{rows}</tbody></table></div>'''
        cooc = P.get('gene_cooccurrence', {})
        if cooc:
            pairs = []
            for g1, partners in cooc.items():
                for g2, cnt in partners.items():
                    pairs.append((g1, g2, cnt))
            pairs.sort(key=lambda x: x[2], reverse=True)
            rows = ''.join(f'<tr><td>{esc(g1)}</td><td>{esc(g2)}</td><td>{cnt}</td></tr>'
                           for g1, g2, cnt in pairs[:500])
            html += f'''<h3>📈 Gene Co-occurrence (Top 500)</h3>
            <div class="master-scrollable-container"><table class="data-table">
                <thead><tr><th>Gene 1</th><th>Gene 2</th><th>Co-occurrence</th></tr></thead>
                <tbody>{rows}</tbody></table></div>'''
        return html

    # ------------------------------------------------------
    # Comparison of two samples
    #------------------------------------------------------

    def _sec_compare(self, d):
        """Pairwise + cluster comparison with colored visuals."""
        samples_data = d.get('samples', {})
        mutation_data = d.get('mutation_data', {})
        mutation_list = mutation_data.get('mutations', [])

        if len(samples_data) < 2:
            return self._alert('warning', 'fa-balance-scale',
                '<h3>Not Enough Samples</h3>'
                f'<p>Compare requires at least two samples. This dataset has '
                f'{len(samples_data)}.</p>')

        # --- Per-sample mutation map ---
        sample_mutations: Dict[str, List[str]] = defaultdict(list)
        for m in mutation_list:
            label = f"{m.get('gene', '')} {m.get('mutation', '')}".strip()
            for genome in m.get('genomes', []):
                sample_mutations[genome].append(label)

        # --- Build the JSON blob per sample ---
        compare_data: Dict[str, Dict[str, Any]] = {}

        def db_genes(sample: str, db: str) -> List[str]:
            return list(samples_data.get(sample, {})
                        .get('abricate_databases', {}).get(db, []))

        for sample, sd in samples_data.items():
            t = sd.get('typing', {})
            amr_genes = sorted(set(
                sd.get('amrfinder', {}).get('all_genes', [])
                + db_genes(sample, 'card')
                + db_genes(sample, 'resfinder')
                + db_genes(sample, 'argannot')
                + db_genes(sample, 'megares')
                + db_genes(sample, 'ncbi')
            ))
            virulence_genes = sorted(set(db_genes(sample, 'vfdb')))
            bacmet_genes = sorted(set(db_genes(sample, 'bacmet2')))
            plasmid_genes = sorted(set(db_genes(sample, 'plasmidfinder')))
            mutations = sorted(set(sample_mutations.get(sample, [])))

            compare_data[sample] = {
                'typing': {
                    'MLST':           t.get('MLST', 'Not Assigned'),
                    'spa':            t.get('spa_Type', 'Not Assigned'),
                    'SCCmec (CGE)':   t.get('SCCmec_CGE', 'Not Assigned'),
                    'SCCmec (RPet)':  t.get('SCCmec_RPet', 'Not Assigned'),
                    'SCCmec Subtype': t.get('SCCmec_Subtype', 'Not Assigned'),
                    'agr':            t.get('agr_Type', 'Not Assigned'),
                    'Capsule':        t.get('capsule_type', 'Not Assigned'),
                    'MRSA':           t.get('MRSA_Status', 'Not Assigned'),
                },
                'genes': {
                    'AMR':       amr_genes,
                    'Virulence': virulence_genes,
                    'BACMET':    bacmet_genes,
                    'Plasmids':  plasmid_genes,
                    'Mutations': mutations,
                },
            }

        sorted_samples = sorted(compare_data.keys())
        options_a = ''.join(
            f'<option value="{esc(s)}"{" selected" if i == 0 else ""}>{esc(s)}</option>'
            for i, s in enumerate(sorted_samples)
        )
        options_b = ''.join(
            f'<option value="{esc(s)}"{" selected" if i == 1 else ""}>{esc(s)}</option>'
            for i, s in enumerate(sorted_samples)
        )
        blob_json = json.dumps(compare_data, default=str, ensure_ascii=False)

        return f'''
        <div class="alert-box alert-info">
            <i class="fas fa-balance-scale fa-2x"></i>
            <div>
                <h3>⚖️ Compare &amp; Cluster Analysis</h3>
                <p>Two modes: <strong>Pair Compare</strong> for a deep side-by-side of two isolates, and
                <strong>Cluster Detection</strong> to find groups of highly similar samples across the whole dataset.</p>
                <ul style="margin-top:6px;">
                    <li><strong>Outbreak confirmation</strong> — do two isolates share MLST + spa + SCCmec + agr?</li>
                    <li><strong>Discordance analysis</strong> — what changed between an ancestor and a resistant descendant?</li>
                    <li><strong>Cluster detection</strong> — which samples form tight transmission clusters (≥95% similarity)?</li>
                </ul>
            </div>
        </div>

        <div class="mode-toggle">
            <button class="mode-btn active" data-mode="pair" onclick="setCompareMode('pair')">
                <i class="fas fa-balance-scale"></i> Pair Compare
            </button>
            <button class="mode-btn" data-mode="cluster" onclick="setCompareMode('cluster')">
                <i class="fas fa-project-diagram"></i> Cluster Detection
            </button>
        </div>

        <!-- ================= PAIR MODE ================= -->
        <div id="mode-pair" class="compare-mode active">
            <div class="compare-panel">
                <div class="compare-select-row">
                    <label for="compare-a">
                        <span>Sample A</span>
                        <select id="compare-a">{options_a}</select>
                    </label>
                    <button class="action-btn btn-primary" onclick="runCompare()">
                        <i class="fas fa-balance-scale"></i> Compare
                    </button>
                    <button class="action-btn btn-light" onclick="swapCompare()">
                        <i class="fas fa-exchange-alt"></i> Swap
                    </button>
                    <button class="action-btn btn-light" onclick="resetCompare()">
                        <i class="fas fa-sync"></i> Reset
                    </button>
                    <label for="compare-b">
                        <span>Sample B</span>
                        <select id="compare-b">{options_b}</select>
                    </label>
                </div>
                <div id="compare-verdict" style="display:none;"></div>
                <div id="compare-visuals"></div>
                <div id="compare-output"></div>
            </div>
        </div>

        <!-- ================= CLUSTER MODE ================= -->
        <div id="mode-cluster" class="compare-mode">
            <div class="compare-panel">
                <div class="cluster-controls">
                    <label>
                        <span>Similarity threshold:</span>
                        <select id="cluster-threshold" onchange="runCluster()">
                            <option value="95">Tight clusters (≥95%)</option>
                            <option value="90" selected>Moderate clusters (≥90%)</option>
                            <option value="85">Loose clusters (≥85%)</option>
                            <option value="80">Broad clusters (≥80%)</option>
                        </select>
                    </label>
                    <button class="action-btn btn-primary" onclick="runCluster()">
                        <i class="fas fa-project-diagram"></i> Detect Clusters
                    </button>
                    <span id="cluster-summary" class="results-counter"></span>
                </div>
                <div id="cluster-clusters"></div>
                <div id="cluster-heatmap"></div>
            </div>
        </div>

        <script>
        window.STAPHSCOPE_COMPARE_DATA = {blob_json};
        </script>'''

    # =========================================================================
    # AI GUIDE
    # =========================================================================
    def _sec_aiguide(self, d):
        """Long-form AI assistant guide with example questions and ethics."""
        return '''
        <div class="alert-box alert-info">
            <i class="fas fa-robot fa-2x"></i>
            <div>
                <h3>🤖 AI Assistant Guide – Unleash the Power of AI for Genomic Epidemiology</h3>
                <p>Use large language models (LLMs) like <strong>ChatGPT, Claude, or Gemini</strong> to interact with your <em>S. aureus</em> dataset — turning static reports into dynamic conversations.</p>
            </div>
        </div>
        <div style="margin:20px 0;">
            <div class="database-section">
                <h4><i class="fas fa-brain"></i> Why Use AI for Genomic Data Analysis?</h4>
                <ul>
                    <li><strong>Pattern recognition</strong> — spot epidemiological trends, clone associations, co-occurrence networks.</li>
                    <li><strong>Natural-language queries</strong> — ask in plain English, get instant answers without writing code.</li>
                    <li><strong>Hypothesis generation</strong> — uncover unexpected correlations that merit experimental follow-up.</li>
                    <li><strong>Literature synthesis</strong> — connect findings with published resistance mechanisms and clinical guidelines.</li>
                </ul>
                <p><span style="background:#fff3cd;padding:2px 8px;border-radius:4px;"><i class="fas fa-lightbulb"></i> <strong>Scientific note:</strong> AI is <em>pattern-finding</em>, not <em>causal-inferring</em>. Use it to suggest, then verify with wet-lab or clinical correlation.</span></p>
            </div>
            <div class="database-section">
                <h4><i class="fas fa-upload"></i> How to Feed This Report to AI</h4>
                <ol>
                    <li><strong>Upload the JSON file</strong> — <code>staphscope_ultimate_gene_centric_report.json</code> contains all structured data. Upload to ChatGPT (Advanced Data Analysis), Claude, or Gemini. <em>Best for quantitative queries.</em></li>
                    <li><strong>Upload the HTML report</strong> — modern AI tools parse HTML tables. Upload the <code>.html</code> file directly. <em>Great for visual context.</em></li>
                    <li><strong>Copy-paste specific tables</strong> — if you only need a quick insight, paste a table into chat. <em>Instant, no file upload needed.</em></li>
                </ol>
                <p style="margin-top:10px;background:#e8f5e9;padding:10px;border-radius:5px;">
                    <i class="fas fa-info-circle"></i> <strong>Pro tip:</strong> tell the AI: <em>"You are a bioinformatician analysing S. aureus genomes. The attached data contains typing, AMR, virulence, BACMET, mutation, and MGE information. Answer my questions with references to the data."</em>
                </p>
            </div>
            <div class="database-section">
                <h4><i class="fas fa-chart-line"></i> Scientifically Relevant Questions to Ask</h4>
                <div style="display:grid;grid-template-columns:1fr 1fr;gap:10px;">
                    <div style="background:#f8f9fa;padding:10px;border-radius:8px;">
                        <strong>🧬 Epidemiology &amp; Clonality</strong>
                        <ul style="margin-top:5px;font-size:.9em;">
                            <li>What are the most common MLST sequence types in this dataset?</li>
                            <li>Which spa types are dominant in MRSA vs MSSA?</li>
                            <li>Are there any ST-spa-SCCmec combinations with &gt;2 isolates?</li>
                            <li>What is the agr type distribution, and does it correlate with MRSA status?</li>
                            <li>Which clones carry the most resistance genes?</li>
                        </ul>
                    </div>
                    <div style="background:#f8f9fa;padding:10px;border-radius:8px;">
                        <strong>💊 Antimicrobial Resistance</strong>
                        <ul style="margin-top:5px;font-size:.9em;">
                            <li>How many samples carry mecA? What are their STs and spa types?</li>
                            <li>Are there any vanA/vanB positive samples? What is their SCCmec type?</li>
                            <li>Which AMR genes co-occur most frequently?</li>
                            <li>What is the distribution of tetracycline (tet) resistance genes?</li>
                            <li>Do any samples have combined β-lactam + macrolide resistance?</li>
                        </ul>
                    </div>
                    <div style="background:#f8f9fa;padding:10px;border-radius:8px;">
                        <strong>🦠 Virulence &amp; Toxins</strong>
                        <ul style="margin-top:5px;font-size:.9em;">
                            <li>Which samples carry PVL (lukF/S-PV)? Are they associated with specific STs or agr types?</li>
                            <li>List all samples with TSST-1 (tsst).</li>
                            <li>Which enterotoxin genes are most prevalent?</li>
                            <li>Is there a correlation between biofilm (ica) genes and MRSA?</li>
                            <li>Do any isolates carry both immune evasion and cytotoxin genes?</li>
                        </ul>
                    </div>
                    <div style="background:#f8f9fa;padding:10px;border-radius:8px;">
                        <strong>🧪 Mutations, MGE &amp; Biocides</strong>
                        <ul style="margin-top:5px;font-size:.9em;">
                            <li>What are the most frequent point mutations in gyrA or parC?</li>
                            <li>Are there any linezolid-related mutations (23S rRNA)?</li>
                            <li>Which samples carry qac genes (disinfectant resistance)?</li>
                            <li>Is mer (mercury) resistance linked to specific STs?</li>
                            <li>Which isolates have the highest mobileOG hit counts?</li>
                        </ul>
                    </div>
                </div>
            </div>
            <div class="database-section">
                <h4><i class="fas fa-balance-scale"></i> Scientific Rigour &amp; Ethical AI Use</h4>
                <ul>
                    <li><strong>AI is your co-pilot, not the pilot.</strong> Interpret AI-generated insights in the context of local epidemiology, clinical guidelines, and lab validation.</li>
                    <li><strong>Verify, verify, verify.</strong> Cross-check critical calls with primary literature, genome browsers, or secondary tools.</li>
                    <li><strong>No patient-identifiable data.</strong> Only upload aggregated, de-identified genomic data.</li>
                    <li><strong>Transparency in publications.</strong> Mention AI-assisted pattern discovery in methods, followed by manual curation.</li>
                    <li><strong>AI hallucination is real.</strong> If the AI confidently tells you that <em>mecA</em> is found in <em>E. coli</em> or that TSST-1 causes indigestion — <strong>don't believe it</strong>. Treat every AI statement as a hypothesis.</li>
                </ul>
            </div>
            <div class="database-section" style="background:#fff3cd;border-left:6px solid #ffc107;">
                <h4><i class="fas fa-smile-wink"></i> A (Mostly Serious) AI Survival Guide</h4>
                <ul>
                    <li><strong>If the AI says "I don't know"</strong> — trust it. It's being honest.</li>
                    <li><strong>If the AI says "It is widely known"</strong> — ask for a reference. It may have made it up.</li>
                    <li><strong>If the AI offers a completely novel evolutionary theory</strong> — check if your coffee is spiked.</li>
                    <li><strong>Remember:</strong> AI won't take your job — but a microbiologist who knows how to use AI might! Learn it, use it, and always keep a healthy dose of scepticism. 😉</li>
                </ul>
            </div>
        </div>'''

    # =========================================================================
    # CALL TO ACTION
    # =========================================================================
    def _sec_calltoaction(self, d):
        """Global AMR narrative and ESCAPE AMR project call to action."""
        return '''
        <div class="alert-box alert-info">
            <i class="fas fa-globe fa-2x"></i>
            <div>
                <h3>The Global Burden of AMR and Our Call to Action</h3>
                <p>Antimicrobial resistance (AMR) is one of the top global public health threats, with an estimated <strong>1.27 million direct deaths annually</strong>. <em>Staphylococcus aureus</em> is a major contributor — skin and soft tissue infections, bloodstream infections, pneumonia, endocarditis. Tracking AMR and virulence determinants is essential to inform treatment guidelines and infection control.</p>
                <p>We developed <strong>StaphScope</strong> to empower researchers and clinicians — especially in low-resource settings — to analyse their own sequencing data without extensive bioinformatics expertise.</p>
            </div>
        </div>
        <div style="background:#e8f5e9;padding:20px;border-radius:12px;margin:20px 0;">
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
            <p>These bacteria “escape” the effects of antibiotics — hence the name. But we believe the name is also a global call to action:</p>
            <div style="background:#fff3e0;padding:15px;border-radius:8px;margin:15px 0;">
                <p><strong>🔹 E</strong>veryone must join forces — researchers, clinicians, policymakers, and citizens.<br>
                <strong>🔹 S</strong>mart surveillance is our first line of defence. No more guessing — we need genomic data.<br>
                <strong>🔹 K</strong>nowledge must be shared openly. No paywalls, no closed silos.<br>
                <strong>🔹 A</strong>frica bears a heavy AMR burden, but African solutions are already emerging.<br>
                <strong>🔹 P</strong>revention is cheaper than cure. Let's stop resistant infections before they spread.<br>
                <strong>🔹 E</strong>very day we delay, more lives are at stake. The time to act is now, not tomorrow.</p>
            </div>
            <p><i class="fas fa-laugh-squint"></i> <strong>“We didn't choose the name ESKAPE because it sounds cool (though it does). We chose it because it reminds us every single day: we must ESCAPE the AMR crisis — together, urgently, and with the best science we have.”</strong><br>— Brown Beckley, lead developer (who secretly hopes this pun makes you smile, not roll your eyes 😉)</p>
        </div>
        <div style="text-align:center;margin:40px 0;">
            <i class="fas fa-star" style="font-size:3em;color:#ffc107;"></i>
            <h3>🤝 We Invite You to Contribute!</h3>
            <p><strong>If you find this tool useful, please:</strong></p>
            <div class="action-buttons" style="justify-content:center;">
                <a href="https://github.com/bbeckley-hub/staphscope-typing-tool" target="_blank" class="action-btn btn-primary" style="text-decoration:none;"><i class="fab fa-github"></i> ⭐ Star us on GitHub</a>
                <a href="mailto:brownbeckley94@gmail.com" class="action-btn btn-success" style="text-decoration:none;"><i class="fas fa-envelope"></i> Contact the team</a>
                <a href="https://github.com/bbeckley-hub/staphscope-typing-tool/issues" target="_blank" class="action-btn btn-warning" style="text-decoration:none;"><i class="fas fa-bug"></i> Report issues</a>
            </div>
            <p style="margin-top:20px;"><i class="fas fa-chalkboard-user"></i> <strong>We welcome collaborations</strong> to adapt this tool for other pathogens and to improve AMR surveillance in Africa and beyond.</p>
            <p><i class="fas fa-hand-holding-heart"></i> If you are a funder or organisation interested in supporting the <strong>ESCAPE AMR</strong> project, please reach out. Together we can build a free, open-source ecosystem for genomic surveillance of all ESKAPE pathogens.</p>
        </div>
        <div style="background:#f8f9fa;padding:15px;border-radius:8px;">
            <i class="fas fa-quote-left"></i> “AMR is a silent pandemic, but we have the tools to fight it — if we share them, if we teach each other, and if we act with urgency. Let's escape the era of untreatable infections.”<br>
            <strong>— The ESCAPE AMR Team, University of Ghana Medical School</strong>
        </div>'''

    # =========================================================================
    # CITATION
    # =========================================================================
    def _sec_citation(self, d):
        """Full citation accordion with clickable DOIs and per-citation colors."""
        palette = [
            '#e11d48', '#dc2626', '#ea580c', '#d97706', '#ca8a04', '#65a30d',
            '#16a34a', '#059669', '#0d9488', '#0891b2', '#0284c7', '#2563eb',
            '#4f46e5', '#7c3aed', '#9333ea', '#c026d3', '#db2777', '#be123c',
            '#f43f5e', '#0ea5e9', '#8b5cf6', '#f97316', '#10b981', '#facc15',
        ]

        main_citations = [
            ('StaphScope',
             'Beckley B, Amarh V. StaphScope: a species-optimized computational pipeline for rapid and accessible <em>Staphylococcus aureus</em> genotyping and surveillance. <em>BMC Genomics</em>. 2026;27:261.',
             'https://doi.org/10.1186/s12864-026-12609-x'),
        ]

        deps_citations = [
            ('MLST',
             'Seemann T. MLST: Scan contig files against PubMLST typing schemes. GitHub. 2018.',
             'https://github.com/tseemann/mlst'),
            ('PubMLST / BIGSdb',
             'Jolley KA, Bray JE, Maiden MCJ. Open-access bacterial population genomics: BIGSdb software, the PubMLST.org website and their applications. <em>Wellcome Open Res</em>. 2018;3:124.',
             'https://doi.org/10.12688/wellcomeopenres.14826.1'),
            ('spa typing',
             'Harmsen D, et al. Typing of methicillin-resistant <em>Staphylococcus aureus</em> in a university hospital setting by using novel software for spa repeat determination and database management. <em>J Clin Microbiol</em>. 2003;41(12):5442-8.',
             'https://doi.org/10.1128/JCM.41.12.5442-5448.2003'),
            ('SCCmecFinder',
             'Kaya H, et al. SCCmecFinder, a Web-Based Tool for Typing of Staphylococcal Cassette Chromosome <em>mec</em> in <em>Staphylococcus aureus</em> Using Whole-Genome Sequence Data. <em>mSphere</em>. 2018;3(1):e00612-17.',
             'https://doi.org/10.1128/mSphere.00612-17'),
            ('agrVATE',
             'Raghuram V, Alexander AM, Loo HQ, Petit RA 3rd, Goldberg JB, Read TD. Species-Wide Phylogenomics of the <em>Staphylococcus aureus</em> Agr Operon Revealed Convergent Evolution of Frameshift Mutations. <em>Microbiol Spectr</em>. 2022;10(1):e0133421.',
             'https://doi.org/10.1128/spectrum.01334-21'),
            ('Capsule Typing (cap5/cap8)',
             'Sau S, Bhasin N, Wann ER, Lee JC, Foster TJ, Lee CY. The <em>Staphylococcus aureus</em> allelic genetic loci for serotype 5 and 8 capsule expression contain the type-specific genes flanked by common genes. <em>Microbiology (Reading)</em>. 1997;143(Pt 7):2395-2405.',
             'https://doi.org/10.1099/00221287-143-7-2395'),
            ('fastANI',
             'Jain C, Rodriguez-R LM, Phillippy AM, Konstantinidis KT, Aluru S. High throughput ANI analysis of 90K prokaryotic genomes reveals clear species boundaries. <em>Nat Commun</em>. 2018;9(1):5114.',
             'https://doi.org/10.1038/s41467-018-07641-9'),
            ('AMRFinderPlus',
             'Feldgarden M, et al. AMRFinderPlus and the Reference Gene Catalog facilitate examination of the genomic links among antimicrobial resistance, stress response, and virulence. <em>Sci Rep</em>. 2021;11(1):12728.',
             'https://doi.org/10.1038/s41598-021-91456-0'),
            ('ABRicate',
             'Seemann T. ABRicate: mass screening of contigs for antibiotic resistance genes. GitHub. 2024.',
             'https://github.com/tseemann/abricate'),
            ('CARD',
             'McArthur AG, et al. The comprehensive antibiotic resistance database. <em>Antimicrob Agents Chemother</em>. 2013;57(7):3348-57.',
             'https://doi.org/10.1128/AAC.00419-13'),
            ('ResFinder',
             'Florensa AF, et al. ResFinder – an open online resource for identification of antimicrobial resistance genes in next-generation sequencing data and prediction of phenotypes from genotypes. <em>Microb Genom</em>. 2022;8(1):000748.',
             'https://doi.org/10.1099/mgen.0.000748'),
            ('VFDB',
             'Chen L, et al. VFDB 2012 update: toward the genetic diversity and molecular evolution of bacterial virulence factors. <em>Nucleic Acids Res</em>. 2012;40(Database issue):D641-5.',
             'https://doi.org/10.1093/nar/gkr989'),
            ('PlasmidFinder',
             'Carattoli A, et al. <em>In silico</em> detection and typing of plasmids using PlasmidFinder and plasmid multilocus sequence typing. <em>Antimicrob Agents Chemother</em>. 2014;58(7):3895-903.',
             'https://doi.org/10.1128/AAC.02412-14'),
            ('BacMet',
             'Pal C, et al. BacMet: antibacterial biocide and metal resistance genes database. <em>Nucleic Acids Res</em>. 2014;42(Database issue):D737-43.',
             'https://doi.org/10.1093/nar/gkt1252'),
            ('MEGARes',
             'Doster E, et al. MEGARes 2.0: a database for classification of antimicrobial drug, biocide and metal resistance determinants in metagenomic sequence data. <em>Nucleic Acids Res</em>. 2020;48(D1):D561-D569.',
             'https://doi.org/10.1093/nar/gkz1010'),
            ('ARG-ANNOT',
             'Gupta SK, et al. ARG-ANNOT, a new bioinformatic tool to discover antibiotic resistance genes in bacterial genomes. <em>Antimicrob Agents Chemother</em>. 2014;58(1):212-20.',
             'https://doi.org/10.1128/AAC.01310-13'),
            ('mobileOG-db',
             'Brown CL, Mullet J, Hindi F, Stoll JE, Gupta S, Choi M, Keenum I, Vikesland P, Pruden A, Zhang L. mobileOG-db: a Manually Curated Database of Protein Families Mediating the Life Cycle of Bacterial Mobile Genetic Elements. <em>Appl Environ Microbiol</em>. 2022;88(18):e00991-22.',
             'https://doi.org/10.1128/aem.00991-22'),
            ('Prodigal',
             'Hyatt D, et al. Prodigal: prokaryotic gene recognition and translation initiation site identification. <em>BMC Bioinformatics</em>. 2010;11:119.',
             'https://doi.org/10.1186/1471-2105-11-119'),
            ('DIAMOND',
             'Buchfink B, Xie C, Huson DH. Fast and sensitive protein alignment using DIAMOND. <em>Nat Methods</em>. 2015;12:59-60.',
             'https://doi.org/10.1038/nmeth.3176'),
            ('Biopython',
             'Cock PJ, et al. Biopython: freely available Python tools for computational molecular biology and bioinformatics. <em>Bioinformatics</em>. 2009;25(11):1422-3.',
             'https://doi.org/10.1093/bioinformatics/btp163'),
        ]

        def render_list(items, start_idx=0):
            html = ''
            for i, (name, text, url) in enumerate(items):
                color = palette[(start_idx + i) % len(palette)]
                plain = re.sub(r'<[^>]+>', '', f"{name} – {text}")
                copy_payload = plain.replace('"', '&quot;')
                link_btn = ''
                if url:
                    link_btn = (
                        f'<a class="citation-link" href="{url}" target="_blank" '
                        f'rel="noopener noreferrer" style="background:{color};">'
                        f'🔗 Open ↗</a>'
                    )
                html += (
                    f'<li class="citation-item" style="border-left:4px solid {color};">'
                    f'<div class="citation-body">'
                    f'<div class="citation-line">'
                    f'<strong class="citation-name" style="color:{color};">{name}</strong>'
                    f'<span class="citation-text"> – {text}</span>'
                    f'</div>'
                    f'<div class="citation-actions">'
                    f'{link_btn}'
                    f'<button class="copy-btn" data-citation="{copy_payload}">📋 Copy</button>'
                    f'</div></div></li>'
                )
            return html

        main_html = render_list(main_citations, start_idx=0)
        deps_html = render_list(deps_citations, start_idx=1)

        return f'''
    <style>
        .citation-item{{background:#fafbfc;border-radius:6px;margin-bottom:10px;padding:12px 14px;list-style:none;transition:box-shadow .2s,transform .2s;box-shadow:0 1px 3px rgba(0,0,0,.05);}}
        .citation-item:hover{{box-shadow:0 4px 12px rgba(0,0,0,.10);transform:translateX(2px);}}
        .citation-body{{display:flex;flex-direction:column;gap:8px;font-size:.92em;line-height:1.55;}}
        .citation-line{{display:block;}}
        .citation-name{{font-size:1em;font-weight:700;}}
        .citation-text{{color:#333;}}
        .citation-actions{{display:flex;gap:8px;flex-wrap:wrap;align-items:center;}}
        .citation-link{{display:inline-flex;align-items:center;gap:4px;padding:4px 14px;border-radius:16px;font-size:.82em;font-weight:600;color:white;text-decoration:none;transition:opacity .2s,transform .2s;}}
        .citation-link:hover{{opacity:.88;transform:translateY(-1px);}}
        .citation-actions .copy-btn{{background:#6b7280;color:white;border:none;padding:4px 14px;border-radius:16px;cursor:pointer;font-size:.82em;font-weight:600;transition:background .2s;}}
        .citation-actions .copy-btn:hover{{background:#4b5563;}}
    </style>
    <div class="alert-box alert-info">
        <i class="fas fa-quote-right fa-2x"></i>
        <div>
            <h3>📚 How to Cite StaphScope and Its Dependencies</h3>
            <p>If you use StaphScope in your research, please cite the main tool and the relevant third-party tools and databases. Each entry below is colour-coded and links directly to its DOI or source repository.</p>
        </div>
    </div>
    <div class="accordion">
        <div class="accordion-item">
            <div class="accordion-header">
                <span>📄 From the ESKAPE AMR Platform</span>
                <i class="fas fa-chevron-down"></i>
            </div>
            <div class="accordion-content" style="display:block;">
                <ul style="padding-left:0;margin:0;">{main_html}</ul>
            </div>
        </div>
        <div class="accordion-item">
            <div class="accordion-header">
                <span>🔧 Key Databases &amp; Methods</span>
                <i class="fas fa-chevron-down"></i>
            </div>
            <div class="accordion-content" style="display:block;">
                <ul style="padding-left:0;margin:0;">{deps_html}</ul>
            </div>
        </div>
    </div>
    <div class="alert-box alert-success" style="margin-top:20px;">
        <i class="fas fa-hand-peace"></i>
        <div>
            <strong>Suggested acknowledgement:</strong><br>
            "Genomic analysis was performed using StaphScope [Beckley &amp; Amarh, 2026], which integrates MLST [Seemann, 2018] using the PubMLST database [Jolley et al., 2018], ABRicate [Seemann, 2018], AMRFinderPlus [Feldgarden et al., 2021], SCCmecFinder [Kaya et al., 2018], agrVATE [Raghuram et al., 2022], and fastANI [Jain et al., 2018] for comprehensive <em>S. aureus</em> characterization. Capsule typing was performed using the cap5/cap8 locus reference [Sau et al., 1997]. Antimicrobial resistance genes were identified using the CARD [McArthur et al., 2013], ResFinder [Florensa et al., 2022], MEGARes [Doster et al., 2020], and ARG-ANNOT [Gupta et al., 2014] databases. For biocide and heavy metal resistance genes, BacMet [Pal et al., 2014] was used. Virulence and plasmid screening were performed with ABRicate using the VFDB [Chen et al., 2012] and PlasmidFinder [Carattoli et al., 2014] databases. Mutation detection was performed using AMRFinderPlus. Mobile genetic element profiling used mobileOG-db [Brown et al., 2022], Prodigal [Hyatt et al., 2010], and DIAMOND [Buchfink et al., 2015]. FASTA QC was performed using Biopython [Cock et al., 2009]."
        </div>
    </div>'''

    # =========================================================================
    # FUNDING
    # =========================================================================
    def _sec_funding(self, d):
        """Funding statement, support options, and ESKAPE AMR contribution invite."""
        return '''
        <div class="alert-box alert-info">
            <i class="fas fa-coffee fa-2x"></i>
            <div>
                <h3>☕ Funding &amp; Support – Keeping the Lights On (with code and caffeine)</h3>
                <p>StaphScope is an <strong>independent, unfunded project</strong> born out of passion for genomic surveillance and AMR research at the University of Ghana Medical School.</p>
                <p>No grants, no sponsors, no institutional backing — just a laptop, a lot of coffee, and a burning desire to help researchers fight antimicrobial resistance.</p>
            </div>
        </div>
        <div class="alert-box alert-warning">
            <i class="fas fa-heart fa-2x"></i>
            <div>
                <h3>💡 How You Can Help (Without Opening Your Wallet)</h3>
                <ul>
                    <li><strong>⭐ Star us on GitHub</strong> – It takes two seconds and makes us feel like rockstars.</li>
                    <li><strong>🐛 Report bugs</strong> – If something breaks, let us know. We'll fix it with joy.</li>
                    <li><strong>💡 Suggest features</strong> – Have an idea? We're all ears (and we actually implement them).</li>
                    <li><strong>🧬 Share your data</strong> – If you've used StaphScope and want to collaborate, we'd love to hear your story.</li>
                    <li><strong>📢 Spread the word</strong> – Tell your colleagues, tweet about it, or mention it in your next Zoom call.</li>
                    <li><strong>👋 Say hello!</strong> – Just drop an email to <strong>brownbeckley94@gmail.com</strong>. It makes our day.</li>
                </ul>
                <p><i class="fas fa-microbe"></i> <strong>Fun fact:</strong> This project runs on 100% volunteer tears, 0% grant money. But we're not bitter — we're just caffeinated.</p>
            </div>
        </div>
        <div class="alert-box alert-success">
            <i class="fas fa-hand-holding-heart"></i>
            <div>
                <h3>🤝 Contribute to the ESKAPE AMR Platform</h3>
                <p>We also maintain pipelines for other ESKAPE pathogens (AcinetoScope, Kleboscope, Pseudoscope, etc.). If you're a developer, bioinformatician, or just someone who loves clean code and bacteria, we welcome pull requests, issues, documentation improvements, and ideas for new databases.</p>
                <p>Visit our GitHub: <a href="https://github.com/bbeckley-hub" target="_blank">https://github.com/bbeckley-hub</a> — star, fork, and let's fight AMR together!</p>
                <p><strong>Brown Beckley</strong> — <i class="fas fa-envelope"></i> brownbeckley94@gmail.com</p>
                <p><i class="fas fa-laugh-beam"></i> <strong>P.S.</strong> If you ever meet Brown in person, buy him a coffee. He'll probably talk your ear off about SCCmec types, but it's worth it.</p>
            </div>
        </div>'''

    # =========================================================================
    # EXPORT
    # =========================================================================
    def _sec_export(self, d):
        """Export panel with one-click CSV downloads and full JSON access."""
        exports = [
            ('samples-table', 'Sample Overview CSV', 'fa-table', 'sample_overview.csv'),
            ('amr-table', 'AMR Genes CSV', 'fa-biohazard', 'amr_genes.csv'),
            ('vir-table', 'Virulence Genes CSV', 'fa-virus', 'virulence_genes.csv'),
            ('bac-table', 'BACMET Genes CSV', 'fa-flask', 'bacmet_genes.csv'),
            ('plasmid-table', 'Plasmid Replicons CSV', 'fa-plug', 'plasmid_replicons.csv'),
            ('mutation-table', 'Mutations CSV', 'fa-dna', 'mutations.csv'),
            ('qc-table', 'FASTA QC CSV', 'fa-chart-line', 'fasta_qc.csv'),
            ('mge-table', 'MGE Profile CSV', 'fa-mobile-alt', 'mge_profile.csv'),
        ]
        cards = ''
        for tid, label, icon, fname in exports:
            cards += f'''<div class="dashboard-card card-export"
                onclick="exportTableToCSV('{tid}', '{fname}')">
                <i class="fas {icon} fa-2x" style="color:#9E9E9E;"></i>
                <div class="card-label">{label}</div></div>'''
        cards += '''<div class="dashboard-card card-export"
            onclick="location.href='staphscope_ultimate_gene_centric_report.json'">
            <i class="fas fa-file-code fa-2x" style="color:#9E9E9E;"></i>
            <div class="card-label">Complete JSON Data</div></div>'''
        return f'''
        {self._alert('info', 'fa-download',
            '<h3>📥 Export Data</h3>'
            '<p>Download any table as CSV, or the full dataset as JSON.</p>')}
        <div style="display:grid;grid-template-columns:repeat(auto-fit,minmax(240px,1fr));gap:20px;margin-top:20px;">
            {cards}
        </div>'''


# =============================================================================
# ORCHESTRATOR
# =============================================================================
class StaphUltimateReporter:
    """Coordinates file discovery, parsing, analysis, and report generation."""

    def __init__(self, input_dir: Path):
        self.input_dir = Path(input_dir)
        self.output_dir = self.input_dir / "STAPHSCOPE_ULTIMATE_GENE_CENTRIC_REPORTS"
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.parser = StaphHTMLParser()
        self.analyzer = StaphDataAnalyzer()
        self.generator = StaphHTMLGenerator(self.analyzer)
        self.metadata = {
            "tool_name": "STAPHSCOPE Ultimate S. aureus Reporter",
            "version": "3.0.0",
            "author": "Brown Beckley <brownbeckley94@gmail.com>",
            "affiliation": "University of Ghana Medical School",
            "analysis_date": datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            "input_directory": str(self.input_dir),
        }

    def find_files(self) -> Dict[str, Optional[Path]]:
        """Discover all recognized input files in the input directory."""
        print("🔍 Searching for StaphScope reports...")
        found = {
            'master_tsv': None,
            'qc_html': None,
            'amrfinder_html': None,
            'abricate_htmls': [],
            'mutation_html': None,
            'mge_html': None,
        }
        master = self.input_dir / 'staphscope_comprehensive_report.tsv'
        if master.exists():
            found['master_tsv'] = master
            print(f"    ✓ Master TSV: {master.name}")

        for cand in ('FASTA_QC_summary.html', 'fasta_qc_summary.html'):
            p = self.input_dir / cand
            if p.exists():
                found['qc_html'] = p
                print(f"    ✓ QC HTML: {cand}")
                break

        for f in self.input_dir.glob('**/*.html'):
            if 'amrfinder_summary_report' in f.name.lower():
                found['amrfinder_html'] = f
                print(f"    ✓ AMRfinder: {f.name}")
                break

        for f in self.input_dir.glob('**/*.html'):
            fl = f.name.lower()
            if any(f'{db}_summary_report' in fl for db in self.parser.abricate_databases):
                found['abricate_htmls'].append(f)
                print(f"    ✓ ABRicate: {f.name}")

        for f in self.input_dir.glob('**/mutation_summary.html'):
            found['mutation_html'] = f
            print(f"    ✓ Mutation: {f.name}")
            break

        for f in self.input_dir.glob('**/*.html'):
            if 'mge' in f.name.lower() and 'summary' in f.name.lower():
                found['mge_html'] = f
                print(f"    ✓ MGE: {f.name}")
                break

        return found

    def integrate(self, files: Dict) -> Dict[str, Any]:
        """Run all parsers, merge results, and build derived tables."""
        print("\n🔗 Integrating data...")
        data = {
            'metadata': self.metadata,
            'samples': {},
            'patterns': {},
            'gene_centric': {},
            'qc_data': {},
            'mutation_data': {},
            'mge_data': {},
        }

        typing_data = {}
        if files['master_tsv']:
            typing_data = self.parser.load_master_tsv(files['master_tsv'])
        else:
            print("  ⚠️  Master TSV not found — typing columns will be 'Not Assigned'")

        if files['qc_html']:
            data['qc_data'] = self.parser.parse_qc_report(files['qc_html'])

        if files['mutation_html']:
            data['mutation_data'] = self.parser.parse_mutation_summary_html(files['mutation_html'])

        amr_by_sample, amr_freq = {}, {}
        if files['amrfinder_html']:
            amr_by_sample, amr_freq = self.parser.parse_amrfinder_report(files['amrfinder_html'])

        abricate_by_sample = defaultdict(dict)
        abricate_freq = {}
        for f in files['abricate_htmls']:
            db, by_sample, freq = self.parser.parse_abricate_report(f)
            if db != 'unknown':
                for s, g in by_sample.items():
                    abricate_by_sample[s][db] = g
                abricate_freq[db] = freq

        if files['mge_html']:
            mge_raw = self.parser.parse_mge_summary_html(files['mge_html'])
            mge_raw['analysis'] = self.analyzer.compute_mge_stats(mge_raw)
            data['mge_data'] = mge_raw

        all_samples = set(typing_data) | set(amr_by_sample) | set(abricate_by_sample) \
                      | set(data['qc_data']) | set(data['mge_data'].get('per_sample', {}))
        all_samples = sorted(all_samples)
        if not all_samples:
            print("❌ No samples found in any report")
            return {}

        print(f"📊 Found {len(all_samples)} unique samples")

        for sample in all_samples:
            t = typing_data.get(sample, {})
            typing = {
                'MLST':           t.get('MLST', 'Not Assigned'),
                'spa_Type':       t.get('spa_Type', 'Not Assigned'),
                'agr_Type':       t.get('agr_Type', 'Not Assigned'),
                'capsule_type':   t.get('capsule_type', 'Not Assigned'),
                'SCCmec_CGE':     t.get('SCCmec_CGE', 'Not Assigned'),
                'SCCmec_RPet':    t.get('SCCmec_RPet', 'Not Assigned'),
                'SCCmec_Subtype': t.get('SCCmec_Subtype', 'Not Assigned'),
                'MRSA_Status':    t.get('MRSA_Status', 'Not Assigned'),
            }
            mge_rec = data['mge_data'].get('per_sample', {}).get(sample, {})
            mge_hits = mge_rec.get('mobileOG Hits', 0) if isinstance(mge_rec, dict) else 0
            if not isinstance(mge_hits, (int, float)):
                mge_hits = 0

            data['samples'][sample] = {
                'typing': typing,
                'amrfinder': amr_by_sample.get(sample, {
                    'critical_genes': [], 'high_risk_genes': [], 'all_genes': []
                }),
                'abricate_databases': abricate_by_sample.get(sample, {}),
                'mge_hits': mge_hits,
            }

        data['gene_frequencies'] = {'amrfinder': amr_freq, 'abricate': abricate_freq}

        print("\n🧠 Building gene-centric + pattern tables...")
        data['gene_centric'] = self.analyzer.create_gene_centric_tables(data)
        data['patterns'] = self.analyzer.create_cross_genome_patterns(data)
        self.generator._current_samples_data = data['samples']
        return data

    def write_json(self, data: Dict[str, Any]) -> Path:
        """Write the full integrated dataset as JSON."""
        print("\n📝 Writing JSON report...")
        out = self.output_dir / "staphscope_ultimate_gene_centric_report.json"

        def serial(obj):
            if obj is None or isinstance(obj, (str, int, float, bool)):
                return obj
            if isinstance(obj, (list, tuple, set)):
                return [serial(x) for x in obj]
            if isinstance(obj, dict):
                return {str(k): serial(v) for k, v in obj.items()}
            if isinstance(obj, (Counter, defaultdict)):
                return {str(k): serial(v) for k, v in obj.items()}
            if isinstance(obj, Path):
                return str(obj)
            if hasattr(obj, 'isoformat'):
                return obj.isoformat()
            return str(obj)

        with open(out, 'w', encoding='utf-8') as f:
            json.dump(serial(data), f, indent=2, ensure_ascii=False)
        print(f"    ✅ JSON saved: {out}")
        return out

    def write_csvs(self, data: Dict[str, Any]):
        """Write all flat CSV exports next to the HTML report."""
        print("\n📊 Writing CSV reports...")

        rows = []
        for sample, s in data['samples'].items():
            t = s['typing']
            rows.append({
                'Sample': sample,
                'MLST': t['MLST'],
                'spa_Type': t['spa_Type'],
                'agr_Type': t['agr_Type'],
                'Capsule_Type': t['capsule_type'],
                'SCCmec_CGE': t['SCCmec_CGE'],
                'SCCmec_RPet': t['SCCmec_RPet'],
                'SCCmec_Subtype': t['SCCmec_Subtype'],
                'MRSA_Status': t['MRSA_Status'],
                'Virulence_Gene_Count': len(s.get('abricate_databases', {}).get('vfdb', [])),
                'MobileOG_Hits': s.get('mge_hits', 0),
            })
        pd.DataFrame(rows).to_csv(self.output_dir / "sample_overview.csv", index=False)

        gc = data.get('gene_centric', {})
        total = len(data['samples']) or 1
        for cat, fname in [('amr_databases', 'amr_genes.csv'),
                           ('virulence_databases', 'virulence_genes.csv'),
                           ('bacmet_databases', 'bacmet_genes.csv'),
                           ('plasmid_databases', 'plasmid_replicons.csv')]:
            out = []
            for db, genes in gc.get(cat, {}).items():
                for g in genes:
                    out.append({
                        'Gene': g['gene'],
                        'Database': g['database'],
                        'Count': g['count'],
                        'Percentage': f"{(g['count'] / total) * 100:.1f}%",
                        'Genomes': ';'.join(g.get('genomes', [])),
                    })
            if out:
                pd.DataFrame(out).to_csv(self.output_dir / fname, index=False)

        muts = data.get('mutation_data', {}).get('mutations', [])
        if muts:
            pd.DataFrame([{
                'Gene': m['gene'],
                'Mutation': m['mutation'],
                'Class': m['class'],
                'Subclass': m['subclass'],
                'Count': m['count'],
                'Genomes': ';'.join(m['genomes'])
            } for m in muts]).to_csv(self.output_dir / "mutations.csv", index=False)

        if data.get('qc_data'):
            pd.DataFrame([{'Sample': s, **m} for s, m in data['qc_data'].items()]) \
                .to_csv(self.output_dir / "fasta_qc.csv", index=False)

        mge_ps = data.get('mge_data', {}).get('per_sample', {})
        if mge_ps:
            pd.DataFrame([{'Sample': s, **v} for s, v in mge_ps.items()]) \
                .to_csv(self.output_dir / "mge_profile.csv", index=False)

        P = data['patterns']
        pat_rows = []
        for mlst, cnt in P.get('mlst_distribution', {}).items():
            pat_rows.append({'Pattern_Type': 'MLST_Distribution',
                             'Combination': mlst, 'Count': cnt})
        for key in ('mlst_spa_sccmec_cge', 'mlst_spa_sccmec_agr',
                    'mlst_spa_agr_capsule', 'mlst_spa_subtype_agr'):
            for combo, samples in P.get(key, {}).items():
                pat_rows.append({'Pattern_Type': key, 'Combination': combo,
                                 'Count': len(samples),
                                 'Samples': ';'.join(samples)})
        for c in P.get('high_risk_combinations', []):
            pat_rows.append({'Pattern_Type': 'High_Risk',
                             'Combination': c['sample'],
                             'Count': 1,
                             'Samples': c['sample']})
        if pat_rows:
            pd.DataFrame(pat_rows).to_csv(self.output_dir / "pattern_discovery.csv", index=False)
        print("    ✅ CSVs written")

    def run(self) -> bool:
        """Execute the full pipeline end-to-end."""
        print("=" * 80)
        print("🧬 STAPHSCOPE ULTIMATE S. AUREUS REPORTER v3.0.0")
        print("=" * 80)
        print(f"📁 Input:  {self.input_dir}")
        print(f"📁 Output: {self.output_dir}")

        files = self.find_files()
        if not any([files['master_tsv'], files['qc_html'], files['amrfinder_html'],
                    files['abricate_htmls'], files['mutation_html'], files['mge_html']]):
            print("❌ No recognized input files found")
            return False

        data = self.integrate(files)
        if not data:
            return False

        self.write_json(data)
        self.write_csvs(data)
        self.generator.generate_main_report(data, self.output_dir)

        n = len(data['samples'])
        mrsa = sum(1 for s in data['samples'].values()
                   if 'MRSA' in s['typing']['MRSA_Status'])
        mge_n = len(data['mge_data'].get('per_sample', {}))
        print("\n" + "=" * 80)
        print("✅ REPORT COMPLETE")
        print("=" * 80)
        print(f"   Samples:        {n}")
        print(f"   MRSA:           {mrsa}")
        print(f"   MGE-profiled:   {mge_n}")
        print(f"   Output dir:     {self.output_dir}")
        print("=" * 80)
        return True


def main():
    """Command-line entry point."""
    parser = argparse.ArgumentParser(
        description='STAPHSCOPE Ultimate S. aureus Reporter v3.0.0',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""Examples:\n  python staphscope_ultimate_reporter.py -i /path/to/reports\n\nAuthor: Brown Beckley <brownbeckley94@gmail.com>""")
    parser.add_argument('-i', '--input-dir', required=True,
                        help='Directory containing StaphScope reports')
    parser.add_argument('-o', '--output-dir', help='Custom output directory')
    args = parser.parse_args()

    inp = Path(args.input_dir)
    if not inp.exists():
        print(f"❌ Input directory not found: {inp}")
        sys.exit(1)

    reporter = StaphUltimateReporter(inp)
    if args.output_dir:
        reporter.output_dir = Path(args.output_dir)
        reporter.output_dir.mkdir(parents=True, exist_ok=True)

    if not reporter.run():
        sys.exit(1)


if __name__ == "__main__":
    main()
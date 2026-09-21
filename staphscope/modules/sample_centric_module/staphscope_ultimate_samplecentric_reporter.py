#!/usr/bin/env python3
"""
STAPHSCOPE ULTIMATE REPORTER - HYBRID GENE-CENTRIC & SAMPLE-CENTRIC
===================================================================
Version 2.0.0

Gene-centric for MLST/spa/SCCmec/Patterns.
Sample-centric interactive boxes for AMR, Virulence, BACMET, Plasmids, Mutations.
Lazy-rendered isolate boxes prevent browser crashes (SIGILL) on 1000+ sample datasets.

Highlights (v2.0.0):
- Single master TSV for all typing (MLST, spa, agr, capsule, SCCmec CGE/RPet/Subtype, MRSA)
- Lazy-loaded isolate boxes with "Show Details" toggle
- Multi-database education boxes, cross-DB confidence tiers, acquired vs intrinsic
- Genotype-phenotype caveat and per-DB role cards
- Rich acknowledgment bars with clickable DOIs on every typing tab
- fastANI credit in FASTA QC tab
- Color-coded capsule types (Type 5 green, Type 8 red) in Sample Overview
- Citation accordion with clickable DOIs and 24-color palette

Author: Brown Beckley <brownbeckley94@gmail.com>
Affiliation: University of Ghana Medical School
MIT
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


def esc(v) -> str:
    """Escape a value for safe HTML embedding."""
    if v is None:
        return ""
    return (str(v).replace("&", "&amp;").replace("<", "&lt;")
            .replace(">", "&gt;").replace('"', "&quot;"))


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
        """Strip file extension and directory path from a sample identifier."""
        sample = str(sample_id)
        for ext in ('.fna', '.fasta', '.fa', '.gb', '.gbk', '.gbff',
                    '.txt', '.tsv', '.csv'):
            if sample.endswith(ext):
                sample = sample[:-len(ext)]
        if '/' in sample or '\\' in sample:
            sample = Path(sample).name
        return sample.strip()

    # -------------------------------------------------------------------------
    # TSV LOADERS
    # -------------------------------------------------------------------------
    def load_typing_from_tsv(self, input_dir: Path) -> Dict[str, Dict]:
        """Load the master typing TSV (single source of truth for all typing)."""
        tsv_path = input_dir / 'staphscope_comprehensive_report.tsv'
        if not tsv_path.exists():
            return {}
        df = pd.read_csv(tsv_path, sep='\t', dtype=str).fillna('Not Assigned')
        typing = {}
        for _, row in df.iterrows():
            sample = str(row.get('Sample', '')).strip()
            if not sample or sample == 'Not Assigned':
                continue
            typing[sample] = {
                'MLST':            str(row.get('MLST', 'Not Assigned')).strip(),
                'spa_Type':        str(row.get('spa Type', 'Not Assigned')).strip(),
                'agr_Type':        str(row.get('agr Type', 'Not Assigned')).strip(),
                'capsule_type':    str(row.get('Capsule Type', 'Not Assigned')).strip(),
                'SCCmec_CGE':      str(row.get('SCCmec Type (CGE)', 'Not Assigned')).strip(),
                'SCCmec_RPet':     str(row.get('SCCmec Type (RPet)', 'Not Assigned')).strip(),
                'SCCmec_Subtype':  str(row.get('SCCmec Subtype', 'Not Assigned')).strip(),
                'MRSA_Status':     str(row.get('MRSA/MSSA Status', 'Not Assigned')).strip(),
            }
        print(f"  ✅ Loaded master typing TSV: {len(typing)} samples")
        return typing

    def load_amrfinder_from_tsv(self, input_dir: Path) -> Tuple[Dict[str, List[Dict]], Dict[str, int]]:
        """Load AMRFinderPlus per-sample gene details from TSV."""
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
        """Load per-database ABRicate per-sample details from TSVs."""
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
        """Load per-sample point mutations from TSV."""
        tsv_path = input_dir / 'mutation_summary.tsv'
        if not tsv_path.exists():
            return {}
        df = pd.read_csv(tsv_path, sep='\t')
        mutations_by_sample = defaultdict(list)
        for _, row in df.iterrows():
            sample = row['genome']

            def clean(v):
                return '' if pd.isna(v) else str(v)

            mutations_by_sample[sample].append({
                'gene':       clean(row.get('gene_symbol', '')),
                'mutation':   clean(row.get('element_name', '')),
                'class':      clean(row.get('class', '')),
                'subclass':   clean(row.get('subclass', '')),
                'contig':     clean(row.get('contig_id', '')),
                'start':      clean(row.get('start', '')),
                'stop':       clean(row.get('stop', '')),
                'strand':     clean(row.get('strand', '')),
                'coverage':   clean(row.get('coverage', '')),
                'identity':   clean(row.get('identity', '')),
                'accession':  clean(row.get('accession', '')),
            })
        return dict(mutations_by_sample)

    # -------------------------------------------------------------------------
    # HTML FALLBACK PARSERS
    # -------------------------------------------------------------------------
    def parse_html_table(self, html_content: str, table_index: int = 0) -> pd.DataFrame:
        """Parse the Nth HTML table in a string into a DataFrame."""
        try:
            soup = BeautifulSoup(html_content, 'html.parser')
            tables = soup.find_all('table')
            if not tables or table_index >= len(tables):
                return pd.DataFrame()
            table = tables[table_index]
            rows = table.find_all('tr')
            headers = [th.get_text().strip() for th in rows[0].find_all(['th', 'td'])]
            data = []
            for row in rows[1:]:
                cols = row.find_all(['td', 'th'])
                if cols:
                    row_data = [col.get_text().strip() for col in cols]
                    if len(row_data) == len(headers):
                        data.append(row_data)
            return pd.DataFrame(data, columns=headers) if data else pd.DataFrame()
        except Exception as e:
            print(f"  ⚠️ Table parsing error: {e}")
            return pd.DataFrame()

    def load_qc_from_html(self, file_path: Path) -> Dict[str, Dict]:
        """Parse FASTA QC HTML summary."""
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

    def parse_comprehensive_report(self, file_path: Path) -> Dict[str, Dict]:
        """Fallback: parse typing from comprehensive HTML report."""
        print(f"  🧬 Parsing Comprehensive HTML: {file_path.name}")
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
                return {}
            rows = typing_table.find_all('tr')
            if len(rows) < 2:
                return {}
            headers = [c.get_text().strip() for c in rows[0].find_all(['th', 'td'])]
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
            df.columns = [c.strip() for c in headers[:len(df.columns)]]
            df['normalized_sample'] = df[df.columns[0]].apply(self.normalize_sample_id)
            results = {}
            for _, row in df.iterrows():
                sample = row['normalized_sample']
                results[sample] = {
                    'MLST':           'Not Assigned',
                    'spa_Type':       'Not Assigned',
                    'agr_Type':       'Not Assigned',
                    'capsule_type':   'Not Assigned',
                    'SCCmec_CGE':     'Not Assigned',
                    'SCCmec_RPet':    'Not Assigned',
                    'SCCmec_Subtype': 'Not Assigned',
                    'MRSA_Status':    'Not Assigned',
                }
            print(f"    ✓ Found {len(results)} samples (typing fields empty)")
            return results
        except Exception as e:
            print(f"    ❌ Error parsing comprehensive HTML: {e}")
            return {}

    def parse_amrfinder_report(self, file_path: Path) -> Tuple[Dict[str, List], Dict[str, Dict]]:
        """Fallback: parse AMRFinderPlus HTML summary."""
        print(f"  🧬 Parsing AMRfinder HTML fallback: {file_path.name}")
        try:
            with open(file_path, 'r', encoding='utf-8') as f:
                html_content = f.read()
            soup = BeautifulSoup(html_content, 'html.parser')
            genes_by_genome = {}
            tables = soup.find_all('table')
            for table in tables:
                t = table.get_text()
                if 'Genome' in t and 'Critical Genes' in t:
                    df_genomes = pd.read_html(io.StringIO(str(table)))[0]
                    genome_col = next((c for c in df_genomes.columns
                                       if 'genome' in c.lower()), df_genomes.columns[0])
                    for _, row in df_genomes.iterrows():
                        sample = self.normalize_sample_id(row[genome_col])
                        genes_by_genome[sample] = {
                            'critical_genes': [], 'high_risk_genes': [], 'all_genes': []
                        }
                    break
            return genes_by_genome, {}
        except Exception as e:
            print(f"    ❌ Error parsing AMRfinder HTML: {e}")
            return {}, {}

    def parse_abricate_report(self, file_path: Path) -> Tuple[str, Dict[str, List], Dict[str, Dict]]:
        """Fallback: parse ABRicate HTML summary."""
        print(f"  🧬 Parsing ABRicate HTML fallback: {file_path.name}")
        try:
            with open(file_path, 'r', encoding='utf-8') as f:
                html_content = f.read()
            db_name = 'unknown'
            fname = file_path.name.lower()
            for db in self.abricate_databases:
                if db in fname:
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
            print(f"    ❌ Error parsing ABRicate HTML: {e}")
            return 'unknown', {}, {}


# -----------------------------------------------------------------------------
# DATA ANALYZER
# -----------------------------------------------------------------------------
class StaphDataAnalyzer:
    """Cross-genome patterns, gene-centric tables, MGE-like aggregations."""

    def __init__(self):
        self.critical_amr_genes = {
            'meca', 'mecc', 'vana', 'vanb', 'vanc',
            'erma', 'ermb', 'ermc', 'msra', 'mphc',
            'tetk', 'tetm', 'tetl'
        }
        self.high_priority_amr = [
            'mecA', 'mecC', 'vanA', 'vanB', 'ermA', 'ermB', 'ermC',
            'msrA', 'mphC', 'tetK', 'tetM', 'aacA-aphD',
            'ant(4\')-Ia', 'ant(6)-Ia', 'aph(3\')-IIIa', 'satA', 'dfrA', 'dfrG', 'cat'
        ]
        self.critical_virulence_genes = {
            'luks-pv', 'lukf-pv', 'tsst', 'sea', 'seb', 'sec', 'sed', 'see',
            'seg', 'seh', 'sei', 'sej', 'sek', 'sel', 'sem', 'sen', 'seo', 'sep',
            'seq', 'ser', 'seu', 'eta', 'etb', 'hla', 'hlb', 'hlg', 'hld',
        }
        self.high_priority_virulence = [
            'lukF-PV', 'lukS-PV', 'tsst', 'sea', 'seb', 'sec', 'sed', 'see',
            'seg', 'seh', 'sei', 'sej', 'sek', 'sel', 'sem', 'sen', 'seo', 'sep',
            'eta', 'etb', 'hla', 'hlb', 'hlg', 'hld'
        ]

    def create_gene_centric_tables(self, integrated_data: Dict[str, Any]) -> Dict[str, Any]:
        """Build gene-centric tables grouped by database category."""
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
                count = data if isinstance(data, (int, float)) else 0
                gene_list.append({
                    'gene': gene, 'database': 'AMRfinder',
                    'frequency': str(count), 'count': int(count),
                    'genomes': []
                })
            gene_centric['amr_databases']['amrfinder'] = sorted(
                gene_list, key=lambda x: x['count'], reverse=True)

        if 'abricate' in integrated_data.get('gene_frequencies', {}):
            for db_name, db_genes in integrated_data['gene_frequencies']['abricate'].items():
                gene_list = []
                for gene, data in db_genes.items():
                    count = data if isinstance(data, (int, float)) else 0
                    gene_list.append({
                        'gene': gene, 'database': db_name.upper(),
                        'frequency': str(count), 'count': int(count),
                        'genomes': []
                    })
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
        """Build combination tables and distributions for typing."""
        patterns = {
            'mlst_distribution': Counter(),
            'spa_type_distribution': Counter(),
            'agr_type_distribution': Counter(),
            'capsule_distribution': Counter(),
            'sccmec_cge_distribution': Counter(),
            'sccmec_rpet_distribution': Counter(),
            'sccmec_subtype_distribution': Counter(),
            'mrsa_status_distribution': Counter(),
            'mlst_spa_combinations': defaultdict(list),
            'mlst_sccmec_combinations': defaultdict(list),
            'spa_sccmec_combinations': defaultdict(list),
            'triple_combinations': defaultdict(list),
            'gene_cooccurrence': defaultdict(Counter),
            'high_risk_combinations': []
        }
        samples_data = integrated_data.get('samples', {})

        def ok(v):
            return v and v not in ('Not Assigned', 'ND', '', 'nan')

        for sample, data in samples_data.items():
            t = data.get('typing', {})
            mlst = t.get('MLST', 'Not Assigned')
            spa = t.get('spa_Type', 'Not Assigned')
            agr = t.get('agr_Type', 'Not Assigned')
            cap = t.get('capsule_type', 'Not Assigned')
            cge = t.get('SCCmec_CGE', 'Not Assigned')
            rpet = t.get('SCCmec_RPet', 'Not Assigned')
            sub = t.get('SCCmec_Subtype', 'Not Assigned')
            mrsa = t.get('MRSA_Status', 'Not Assigned')

            if ok(mlst): patterns['mlst_distribution'][mlst] += 1
            if ok(spa): patterns['spa_type_distribution'][spa] += 1
            if ok(agr): patterns['agr_type_distribution'][agr] += 1
            if ok(cap): patterns['capsule_distribution'][cap] += 1
            if ok(cge): patterns['sccmec_cge_distribution'][cge] += 1
            if ok(rpet): patterns['sccmec_rpet_distribution'][rpet] += 1
            if ok(sub): patterns['sccmec_subtype_distribution'][sub] += 1
            if ok(mrsa): patterns['mrsa_status_distribution'][mrsa] += 1

            if ok(mlst) and ok(spa):
                patterns['mlst_spa_combinations'][f"{mlst} - {spa}"].append(sample)
            if ok(mlst) and ok(cge):
                patterns['mlst_sccmec_combinations'][f"{mlst} - {cge}"].append(sample)
            if ok(spa) and ok(cge):
                patterns['spa_sccmec_combinations'][f"{spa} - {cge}"].append(sample)
            if ok(mlst) and ok(spa) and ok(cge):
                patterns['triple_combinations'][f"{mlst} - {spa} - {cge}"].append(sample)

            amr_genes = data.get('amrfinder', {}).get('all_genes', [])
            vir_genes = data.get('abricate_databases', {}).get('vfdb', [])
            critical_amr = [g for g in amr_genes
                            if any(c in str(g).lower() for c in self.critical_amr_genes)]
            critical_vir = [g for g in vir_genes
                            if any(c in str(g).lower() for c in self.critical_virulence_genes)]
            if critical_amr and critical_vir:
                patterns['high_risk_combinations'].append({
                    'sample': sample, 'mlst': mlst, 'spa_type': spa,
                    'sccmec_type': cge, 'mrsa_status': mrsa, 'agr_type': agr,
                    'critical_amr_genes': critical_amr,
                    'critical_virulence_genes': critical_vir,
                })

        # Convert defaultdict(list) to plain dict for JSON safety
        for key in ('mlst_spa_combinations', 'mlst_sccmec_combinations',
                    'spa_sccmec_combinations', 'triple_combinations'):
            patterns[key] = dict(patterns[key])
        return patterns


# -----------------------------------------------------------------------------
# HTML GENERATOR
# -----------------------------------------------------------------------------
class StaphHTMLGenerator:
    """Builds the interactive multi-tab HTML report with lazy-loaded boxes."""

    def __init__(self, data_analyzer: StaphDataAnalyzer):
        self.data_analyzer = data_analyzer
        self.tab_colors = {
            'summary': '#4CAF50', 'sample_overview': '#2196F3', 'qc': '#607D8B',
            'mlst': '#FF9800', 'spa': '#9C27B0', 'sccmec': '#009688',
            'mrsa': '#795548', 'agr': '#8B5CF6', 'amr': '#F44336',
            'virulence': '#E91E63', 'bacmet': '#FF5722', 'plasmids': '#673AB7',
            'mutation': '#00BCD4', 'patterns': '#3F51B5', 'aiguide': '#00BCD4',
            'citation': '#8BC34A', 'funding': '#FFC107', 'export': '#9E9E9E',
            'calltoaction': '#F472B6'
        }

    # -------------------------------------------------------------------------
    # Reusable HTML helpers
    # -------------------------------------------------------------------------
    def _credit_bar(self, color: str, icon: str, title: str, body: str) -> str:
        """Colored acknowledgement strip used at the top of tool-driven tabs."""
        return f'''
        <div class="scientific-note" style="background:linear-gradient(135deg,#f8f9fa 0%,#f0f4f8 100%);border-left:6px solid {color};margin-bottom:20px;padding:15px;border-radius:8px;">
            <div style="display:flex;align-items:center;gap:12px;flex-wrap:wrap;">
                <span style="font-size:1.4em;">{icon}</span>
                <div>
                    <strong style="font-size:1.1em;color:{color};">{title}</strong><br>
                    <span style="font-size:0.95em;color:#333;">{body}</span>
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
        """Colored stat card."""
        icon_html = (f'<i class="fas {icon} fa-2x" '
                     f'style="opacity:0.9;margin-bottom:8px;"></i>' if icon else '')
        return f'''
        <div class="stat-card" style="background:linear-gradient(135deg,{color} 0%,{color}dd 100%);">
            {icon_html}
            <div class="stat-value">{value}</div>
            <div class="stat-label">{label}</div>
        </div>'''

    def _filter_buttons(self, table_id: str, buttons: list) -> str:
        """Row of quick-filter buttons that populate the table's search box."""
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
        """Info box listing biological role of each gene family."""
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

    def _colorize_capsule_cell(self, value: str) -> str:
        """Color a capsule type value (Type 5 = green, Type 8 = red)."""
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

    def _capsule_badge(self, value: str) -> str:
        """Return a colored typing-badge span for a capsule type value."""
        v = str(value).strip()
        if v == 'Type 5':
            return '<span class="typing-badge" style="background:#d4edda;color:#155724;border-color:#b7e4c0;">Capsule: Type 5</span>'
        if v == 'Type 8':
            return '<span class="typing-badge" style="background:#f8d7da;color:#721c24;border-color:#f5c2c7;">Capsule: Type 8</span>'
        if v and v != 'Not Assigned':
            return f'<span class="typing-badge">Capsule: {esc(v)}</span>'
        return '<span class="typing-badge">Capsule: —</span>'
    # -------------------------------------------------------------------------
    # Multi-DB educational blocks (shared across AMR/Virulence/BACMET)
    # -------------------------------------------------------------------------
    def _multi_db_education(self) -> str:
        return '''
        <div class="alert-box" style="border-left-color:#00695c;background:#e8f5e9;border-radius:8px;padding:18px 22px;margin:20px 0;">
            <div style="display:flex;gap:15px;align-items:flex-start;">
                <i class="fas fa-info-circle fa-2x" style="color:#00695c;margin-top:3px;"></i>
                <div>
                    <h4 style="margin:0 0 10px 0;color:#00695c;font-size:1.1em;">🔬 Why Multiple Databases?</h4>
                    <p style="margin:6px 0;font-size:.95em;">StaphScope screens every genome against <strong>independent AMR databases</strong> because <strong>no single database is comprehensive</strong>. Each has unique strengths, biases, and update cadences.</p>
                    <p style="margin:10px 0 6px 0;font-size:.95em;"><strong>⚠️ The problem with choosing one database:</strong></p>
                    <ul style="margin:6px 0 10px 20px;font-size:.93em;">
                        <li>Some researchers pick a favourite DB (often CARD or ResFinder) and only report hits from that one.</li>
                        <li>This is <strong>fast but incomplete</strong> — a gene absent from CARD may still be present in MEGARes or AMRFinderPlus.</li>
                        <li>Single-DB hits with weak support can be <strong>false positives</strong> that a second database would have flagged.</li>
                        <li>Result: <em>biased prevalence estimates</em> and <em>missed resistance signals</em>.</li>
                    </ul>
                    <p style="margin:10px 0 6px 0;font-size:.95em;"><strong>✅ Our approach — report everything, provenance preserved:</strong></p>
                    <ul style="margin:6px 0 10px 20px;font-size:.93em;">
                        <li>All hits from all databases are kept <strong>separate and unmerged</strong> — the <strong>Database column</strong> tells you exactly which source found each gene.</li>
                        <li>You get the <strong>full picture</strong>; no silent filtering, no cherry-picking.</li>
                        <li>Cross-database agreement becomes a <strong>confidence signal</strong> (see next box).</li>
                    </ul>
                    <p style="margin:6px 0 0 0;font-size:.92em;background:#fff3cd;padding:8px 14px;border-radius:4px;border-left:3px solid #ffc107;">
                        <i class="fas fa-lightbulb" style="color:#856404;"></i>
                        <strong>Pro tip:</strong> Use the database dropdown below to inspect only one DB's hits, or group by typing to see which lineages carry which genes.
                    </p>
                </div>
            </div>
        </div>'''

    def _confidence_tiers(self) -> str:
        return '''
        <div class="alert-box" style="border-left-color:#0891b2;background:#e0f2fe;border-radius:8px;padding:16px 20px;margin:20px 0;">
            <div style="display:flex;gap:15px;align-items:flex-start;">
                <i class="fas fa-layer-group fa-2x" style="color:#0891b2;margin-top:3px;"></i>
                <div>
                    <h4 style="margin:0 0 10px 0;color:#0891b2;font-size:1.05em;">🎯 Cross-Database Confidence Tiers</h4>
                    <p style="margin:6px 0;font-size:.93em;">A gene detected by <strong>multiple databases</strong> is far more likely to be a true positive than a single-DB hit.</p>
                    <ul style="margin:6px 0 0 20px;font-size:.93em;">
                        <li><span style="color:#16a34a;font-weight:bold;">🟢 High confidence</span> — found in <strong>3 or more databases</strong>.</li>
                        <li><span style="color:#f59e0b;font-weight:bold;">🟡 Moderate confidence</span> — found in <strong>2 databases</strong>.</li>
                        <li><span style="color:#dc2626;font-weight:bold;">🔴 Low confidence / investigate</span> — found in <strong>only 1 database</strong>.</li>
                    </ul>
                </div>
            </div>
        </div>'''

    def _acquired_intrinsic(self, species: str = "S. aureus") -> str:
        return f'''
        <div class="alert-box" style="border-left-color:#6f42c1;background:#f3e8ff;border-radius:8px;padding:16px 20px;margin:20px 0;">
            <div style="display:flex;gap:15px;align-items:flex-start;">
                <i class="fas fa-dna fa-2x" style="color:#6f42c1;margin-top:3px;"></i>
                <div>
                    <h4 style="margin:0 0 10px 0;color:#6f42c1;font-size:1.05em;">🧬 Acquired vs Intrinsic Resistance — Both Matter</h4>
                    <p style="margin:6px 0;font-size:.93em;">The AMR story is <strong>more than acquired genes</strong>. We report <strong>both</strong> because they answer different questions:</p>
                    <ul style="margin:6px 0 10px 20px;font-size:.93em;">
                        <li><strong style="color:#7c3aed;">Intrinsic genes</strong> — baseline genome (e.g. <em>norA</em>, <em>mepA</em>, <em>lmrS</em> efflux pumps in {species}). Set the floor for susceptibility.</li>
                        <li><strong style="color:#e11d48;">Acquired genes</strong> — gained by horizontal transfer (<em>mecA</em>, <em>ermC</em>, <em>tetK</em>, <em>dfrG</em>). Predict clinical failure of specific drugs.</li>
                    </ul>
                    <p style="margin:6px 0 0 0;font-size:.92em;"><i class="fas fa-lightbulb" style="color:#6f42c1;"></i> <strong>Why both:</strong> only-acquired reports hide the intrinsic baseline; only-intrinsic reports miss the acquired threat.</p>
                </div>
            </div>
        </div>'''

    def _genotype_phenotype_caveat(self) -> str:
        return '''
        <div class="alert-box alert-warning" style="border-left-color:#ffc107;">
            <i class="fas fa-exclamation-triangle fa-2x"></i>
            <div>
                <h4 style="margin:0 0 8px 0;color:#856404;">⚠️ Gene Presence ≠ Phenotypic Resistance</h4>
                <p style="margin:6px 0;font-size:.93em;">Detecting an AMR gene is <strong>necessary evidence</strong> but not sufficient to declare phenotypic resistance. Several mechanisms can break the link:</p>
                <ul style="margin:6px 0 12px 20px;font-size:.92em;">
                    <li><strong>Silent / truncated genes</strong> — non-functional cassette.</li>
                    <li><strong>Expression regulation</strong> — inducible systems that may be off.</li>
                    <li><strong>Mechanism matters</strong> — <em>erm</em> (rRNA methylation) vs <em>msrA</em> (efflux) both give MLS<sub>B</sub> resistance but differ in spectrum.</li>
                    <li><strong>Naming ambiguity</strong> — ResFinder appends <code>_1</code> to primary alleles.</li>
                    <li><strong>Dose and route</strong> — low-level efflux may be overcome <em>in vivo</em>.</li>
                </ul>
                <p style="margin:8px 0 0 0;font-size:.92em;background:#fff8e1;padding:8px 12px;border-radius:4px;border-left:3px solid #f59e0b;">
                    <strong>Clinical bottom line:</strong> <strong>Antimicrobial Susceptibility Testing (AST)</strong> remains the gold standard. Genomic AMR prediction is a triage and surveillance tool — not a replacement.
                </p>
            </div>
        </div>'''

    def _db_roles_amr(self) -> str:
        return '''
        <div class="database-section" style="margin-top:30px;">
            <h3 style="color:#2c3e50;border-bottom:2px solid #3b82f6;padding-bottom:10px;">
                <i class="fas fa-database"></i> Roles &amp; Strengths of Each Database
            </h3>
            <p style="color:#666;margin-bottom:15px;">Each database is optimised for a different purpose.</p>
            <div style="display:grid;grid-template-columns:repeat(auto-fit,minmax(300px,1fr));gap:20px;margin:20px 0;">
                <div style="background:#f8f9fa;padding:16px;border-radius:10px;border-left:4px solid #28a745;">
                    <h4 style="color:#28a745;margin:0 0 8px 0;">🟢 CARD</h4>
                    <p style="font-size:.9em;margin:0;"><strong>Strengths:</strong> Expert-curated, strict SNP-based cutoffs, precise allele calls.<br><strong>Weaknesses:</strong> May miss novel variants.<br><strong>Best for:</strong> Confident allele-level calls.</p>
                </div>
                <div style="background:#f8f9fa;padding:16px;border-radius:10px;border-left:4px solid #17a2b8;">
                    <h4 style="color:#17a2b8;margin:0 0 8px 0;">🔵 ResFinder</h4>
                    <p style="font-size:.9em;margin:0;"><strong>Strengths:</strong> Highly sensitive, frequently updated.<br><strong>Weaknesses:</strong> Appends <code>_1</code> to primary alleles.<br><strong>Best for:</strong> Broad sensitivity.</p>
                </div>
                <div style="background:#f8f9fa;padding:16px;border-radius:10px;border-left:4px solid #007bff;">
                    <h4 style="color:#007bff;margin:0 0 8px 0;">🔷 NCBI AMR</h4>
                    <p style="font-size:.9em;margin:0;"><strong>Strengths:</strong> Curated, tightly linked to AMRFinderPlus's Reference Gene Catalog.<br><strong>Best for:</strong> Cross-checking AMRFinderPlus.</p>
                </div>
                <div style="background:#f8f9fa;padding:16px;border-radius:10px;border-left:4px solid #fd7e14;">
                    <h4 style="color:#fd7e14;margin:0 0 8px 0;">🟠 MEGARes</h4>
                    <p style="font-size:.9em;margin:0;"><strong>Strengths:</strong> Hierarchy of gene families; includes biocide + metal resistance.<br><strong>Best for:</strong> Environmental co-selection markers.</p>
                </div>
                <div style="background:#f8f9fa;padding:16px;border-radius:10px;border-left:4px solid #6f42c1;">
                    <h4 style="color:#6f42c1;margin:0 0 8px 0;">🟣 ARG-ANNOT</h4>
                    <p style="font-size:.9em;margin:0;"><strong>Strengths:</strong> Historic, well-curated ARG catalogue.<br><strong>Weaknesses:</strong> Updated less frequently.<br><strong>Best for:</strong> Historical comparisons.</p>
                </div>
                <div style="background:#f8f9fa;padding:16px;border-radius:10px;border-left:4px solid #dc3545;">
                    <h4 style="color:#dc3545;margin:0 0 8px 0;">🔴 AMRFinderPlus</h4>
                    <p style="font-size:.9em;margin:0;"><strong>Strengths:</strong> NCBI gold standard; includes <strong>point mutations</strong>.<br><strong>Best for:</strong> Clinical-grade calls; the most defensible single source.</p>
                </div>
            </div>
            <div class="alert-box alert-info" style="border-left-color:#00695c;background:#e8f5e9;">
                <i class="fas fa-link fa-2x" style="color:#00695c;"></i>
                <div>
                    <strong>🧭 You decide — we give you the evidence, not the verdict.</strong>
                    <p style="margin-top:8px;font-size:.95em;">StaphScope deliberately <strong>does not merge or prioritise</strong> hits across databases. Different questions call for different choices:</p>
                    <ul style="margin:8px 0 0 20px;font-size:.93em;">
                        <li><strong>Conservative clinical calls?</strong> Filter to AMRFinderPlus only.</li>
                        <li><strong>Maximum sensitivity?</strong> Keep ResFinder + MEGARes + CARD together.</li>
                        <li><strong>Surveillance?</strong> Report all databases; use cross-DB agreement as confidence.</li>
                        <li><strong>Historical comparison?</strong> Check ARG-ANNOT.</li>
                    </ul>
                    <p style="margin-top:10px;font-size:.92em;background:#fff3cd;padding:8px 14px;border-radius:4px;border-left:3px solid #ffc107;">
                        <i class="fas fa-lightbulb"></i> All hits are visible in the boxes above. Use the DB dropdown to filter. Export to CSV to merge however you need.
                    </p>
                </div>
            </div>
        </div>'''

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
        html = self._create_ultimate_html(
            metadata=metadata,
            samples_data=samples_data,
            patterns=patterns,
            gene_centric=gene_centric,
            integrated_data=integrated_data,
        )
        output_file = output_dir / "staphscope_ultimate_sample_centric_report.html"
        with open(output_file, 'w', encoding='utf-8') as f:
            f.write(html)
        print(f"    ✅ HTML report saved: {output_file}")
        return str(output_file)

    def _create_ultimate_html(self, **kwargs) -> str:
        samples_data = kwargs.get('samples_data', {})
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

        total_amr = sum(len(g) for g in kwargs['gene_centric'].get('amr_databases', {}).values())
        total_vir = sum(len(g) for g in kwargs['gene_centric'].get('virulence_databases', {}).values())
        total_bac = sum(len(g) for g in kwargs['gene_centric'].get('bacmet_databases', {}).values())

        TABS = [
            ('summary', 'Summary', 'chart-pie'),
            ('sample_overview', 'Sample Overview', 'list-alt'),
            ('qc', 'FASTA QC', 'chart-line'),
            ('mlst', 'MLST', 'code-branch'),
            ('spa', 'spa Typing', 'dna'),
            ('sccmec', 'SCCmec', 'shield-alt'),
            ('mrsa', 'MRSA', 'skull-crossbones'),
            ('agr', 'agr Typing', 'dna'),
            ('amr', 'AMR', 'biohazard'),
            ('virulence', 'Virulence', 'virus'),
            ('bacmet', 'BACMET', 'flask'),
            ('plasmids', 'Plasmids', 'plug'),
            ('mutation', 'Mutations', 'dna'),
            ('patterns', 'Patterns', 'project-diagram'),
            ('aiguide', 'AI Guide', 'robot'),
            ('calltoaction', 'Call to Action', 'globe'),
            ('citation', 'Citation', 'book'),
            ('funding', 'Funding', 'coffee'),
            ('export', 'Export', 'download'),
        ]

        nav_html = ''
        for i, (tid, title, icon) in enumerate(TABS):
            active = ' active' if i == 0 else ''
            nav_html += (f'<button class="tab-button {tid}{active}" '
                         f'onclick="switchTab(\'{tid}\')">'
                         f'<i class="fas fa-{icon}"></i> {title}</button>')

        section_methods = {
            'summary': self._generate_summary_section,
            'sample_overview': self._generate_sample_overview_section,
            'qc': self._generate_qc_section,
            'mlst': self._generate_mlst_section,
            'spa': self._generate_spa_section,
            'sccmec': self._generate_sccmec_section,
            'mrsa': self._generate_mrsa_section,
            'agr': self._generate_agr_section,
            'patterns': self._generate_pattern_discovery_section,
            'aiguide': self._generate_aiguide_section,
            'citation': self._generate_citation_section,
            'funding': self._generate_funding_section,
            'calltoaction': lambda kw: self._calltoaction_section(),
            'export': self._generate_export_section,
        }

        tabs_html = ''
        for i, (tid, title, icon) in enumerate(TABS):
            active = ' active' if i == 0 else ''
            color = self.tab_colors.get(tid, '#4CAF50')

            # Sample-centric tabs (AMR, Virulence, BACMET, Plasmids, Mutations) get lazy boxes
            if tid == 'amr':
                content = self._generate_sample_centric_boxes(
                    kwargs, 'amr', 'AMR',
                    ['amrfinder', 'resfinder', 'card', 'argannot', 'megares', 'ncbi'])
            elif tid == 'virulence':
                content = self._generate_sample_centric_boxes(
                    kwargs, 'virulence', 'Virulence', ['vfdb'])
            elif tid == 'bacmet':
                content = self._generate_sample_centric_boxes(
                    kwargs, 'bacmet', 'BACMET', ['bacmet2'])
            elif tid == 'plasmids':
                content = self._generate_sample_centric_boxes(
                    kwargs, 'plasmids', 'Plasmids', ['plasmidfinder'])
            elif tid == 'mutation':
                content = self._generate_mutation_boxes(kwargs)
            elif tid in section_methods:
                content = section_methods[tid](kwargs)
            else:
                content = ''

            tabs_html += f'''
            <div id="{tid}-tab" class="tab-content{active}">
                <h2 class="section-header {tid}-header" style="border-color:{color};">
                    <span><i class="fas fa-{icon}"></i> {title}</span>
                    <button class="print-section-btn" onclick="printSection('{tid}-tab')">
                        <i class="fas fa-print"></i> Print</button>
                </h2>
                {content}
            </div>'''

        dash_cards = f'''
            <div class="dashboard-card card-summary" onclick="switchTab('summary')">
                <div class="card-number">{len(kwargs['samples_data'])}</div>
                <div class="card-label">Total Samples</div></div>
            <div class="dashboard-card card-mlst" onclick="switchTab('mlst')">
                <div class="card-number">{len(kwargs['patterns'].get('mlst_distribution', {}))}</div>
                <div class="card-label">Unique STs</div></div>
            <div class="dashboard-card card-spa" onclick="switchTab('spa')">
                <div class="card-number">{len(kwargs['patterns'].get('spa_type_distribution', {}))}</div>
                <div class="card-label">spa Types</div></div>
            <div class="dashboard-card card-amr" onclick="switchTab('amr')">
                <div class="card-number">{total_amr}</div>
                <div class="card-label">AMR Genes</div></div>
            <div class="dashboard-card card-virulence" onclick="switchTab('virulence')">
                <div class="card-number">{total_vir}</div>
                <div class="card-label">Virulence Genes</div></div>
            <div class="dashboard-card card-patterns" onclick="switchTab('patterns')">
                <div class="card-number">{len(kwargs['patterns'].get('high_risk_combinations', []))}</div>
                <div class="card-label">High-Risk Combos</div></div>
            <div class="dashboard-card card-agr" onclick="switchTab('agr')">
                <div class="card-number">{len(kwargs['patterns'].get('agr_type_distribution', {}))}</div>
                <div class="card-label">agr Types</div></div>
        '''

        return f'''<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>STAPHSCOPE Ultimate S. aureus Report v2.0</title>
    <link rel="stylesheet" href="https://cdnjs.cloudflare.com/ajax/libs/font-awesome/6.4.0/css/all.min.css">
    {css}
    {js}
</head>
<body>
<div class="container">
    <div class="main-header">
        <h1><i class="fas fa-bacteria"></i> STAPHSCOPE Ultimate S. aureus Analysis Report</h1>
        <p>Hybrid Gene-Centric + Sample-Centric — Single-Source Typing, Lazy-Loaded Isolate Boxes</p>
        <div class="metadata-bar">
            <div class="metadata-item"><i class="fas fa-calendar"></i><span>Generated: {kwargs['metadata'].get('analysis_date', 'Unknown')}</span></div>
            <div class="metadata-item"><i class="fas fa-database"></i><span>Samples: {len(kwargs['samples_data'])}</span></div>
            <div class="metadata-item"><i class="fas fa-code-branch"></i><span>Tool: STAPHSCOPE Ultimate v2.0.0</span></div>
            <div class="metadata-item"><i class="fas fa-university"></i><span>University of Ghana Medical School</span></div>
        </div>
    </div>
    <div class="dashboard-grid">{dash_cards}</div>
    <div class="tab-navigation">{nav_html}</div>
    {tabs_html}
    <div class="footer">
        <h3>STAPHSCOPE Ultimate S. aureus Reporter v2.0.0</h3>
        <p>University of Ghana Medical School | Brown Beckley &lt;brownbeckley94@gmail.com&gt;</p>
        <p>Generated on {kwargs['metadata'].get('analysis_date', 'Unknown')}</p>
        <p>⭐ Please give a big STAR on GitHub if you found this useful!</p>
    </div>
</div>
</body>
</html>'''

    # -------------------------------------------------------------------------
    # CSS
    # -------------------------------------------------------------------------
    def _get_css(self) -> str:
        return """
        <style>
        :root {
            --summary-color: #4CAF50; --sample_overview-color: #2196F3;
            --qc-color: #607D8B; --mlst-color: #FF9800; --spa-color: #9C27B0;
            --sccmec-color: #009688; --mrsa-color: #795548; --amr-color: #F44336;
            --virulence-color: #E91E63; --bacmet-color: #FF5722;
            --plasmids-color: #673AB7; --mutation-color: #00BCD4;
            --patterns-color: #3F51B5; --aiguide-color: #00BCD4;
            --citation-color: #8BC34A; --funding-color: #FFC107;
            --export-color: #9E9E9E; --agr-color: #8B5CF6;
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
        .dashboard-card:hover { transform: translateY(-10px); box-shadow: 0 15px 30px rgba(0,0,0,0.2); }
        .card-summary { border-left-color: var(--summary-color); }
        .card-mlst { border-left-color: var(--mlst-color); }
        .card-spa { border-left-color: var(--spa-color); }
        .card-sccmec { border-left-color: var(--sccmec-color); }
        .card-mrsa { border-left-color: var(--mrsa-color); }
        .card-amr { border-left-color: var(--amr-color); }
        .card-virulence { border-left-color: var(--virulence-color); }
        .card-bacmet { border-left-color: var(--bacmet-color); }
        .card-plasmids { border-left-color: var(--plasmids-color); }
        .card-patterns { border-left-color: var(--patterns-color); }
        .card-agr { border-left-color: var(--agr-color); }
        .card-number { font-size: 3em; font-weight: bold; margin: 15px 0; background: linear-gradient(90deg, #006400, #228B22); -webkit-background-clip: text; -webkit-text-fill-color: transparent; }
        .card-label { font-size: 0.9em; color: #555; font-weight: 600; }
        .tab-navigation { display: flex; gap: 5px; margin-bottom: 20px; flex-wrap: wrap; background: white; padding: 15px; border-radius: 12px; box-shadow: 0 5px 20px rgba(0,0,0,0.1); position: sticky; top: 10px; z-index: 100; }
        .tab-button { padding: 12px 20px; background: #f5f5f5; border: none; border-radius: 8px; cursor: pointer; font-weight: 600; color: #666; transition: all 0.3s ease; display: flex; align-items: center; gap: 8px; font-size: 0.9em; }
        .tab-button.active { color: white; }
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
        .data-table td { padding: 12px; border-bottom: 1px solid #e0e0e0; vertical-align: top; }
        .data-table tr:hover { background: #f8f9fa; }
        .scrollable-table { max-height: none; overflow-y: auto; border: 1px solid #e0e0e0; border-radius: 8px; margin: 20px 0; width: 100%; }
        .master-scrollable-container { width: 100%; overflow-x: auto; border: 1px solid #e0e0e0; border-radius: 8px; margin: 20px 0; }
        .genome-list { display: flex; flex-wrap: wrap; gap: 5px; max-height: 200px; overflow-y: auto; padding: 5px; background: #f8f9fa; border-radius: 5px; }
        .genome-tag { display: inline-block; background: #e6ffe6; color: #006400; padding: 3px 10px; border-radius: 12px; font-size: 0.85em; border: 1px solid #b3ffb3; white-space: nowrap; margin: 2px; }
        .genome-tag.highlight { background-color: #ffff99 !important; color: #000 !important; border: 1px solid #ffc107; }
        .search-box { width: 100%; padding: 12px; margin-bottom: 20px; border: 2px solid #e0e0e0; border-radius: 8px; font-size: 1em; transition: all 0.3s ease; }
        .search-box:focus { outline: none; border-color: #006400; box-shadow: 0 0 0 3px rgba(0,100,0,0.1); }
        .badge { display: inline-block; padding: 5px 15px; border-radius: 20px; font-size: 0.85em; font-weight: 600; margin: 2px; }
        .badge-mrsa { background: #8B0000; color: white; }
        .badge-mssa { background: #4682B4; color: white; }
        .badge-critical { background: #DC143C; color: white; }
        .alert-box { padding: 20px; border-radius: 10px; margin: 20px 0; display: flex; align-items: flex-start; gap: 20px; border-left: 5px solid; }
        .alert-success { background: #d4edda; color: #155724; border-left-color: #28a745; }
        .alert-warning { background: #fff3cd; color: #856404; border-left-color: #ffc107; }
        .alert-danger { background: #f8d7da; color: #721c24; border-left-color: #dc3545; }
        .alert-info { background: #d1ecf1; color: #0c5460; border-left-color: #17a2b8; }
        .action-buttons { display: flex; gap: 10px; margin: 20px 0; flex-wrap: wrap; }
        .action-btn { padding: 10px 20px; border: none; border-radius: 8px; cursor: pointer; font-weight: 600; display: flex; align-items: center; gap: 8px; transition: all 0.3s ease; text-decoration: none; }
        .action-btn:hover { transform: translateY(-2px); box-shadow: 0 5px 15px rgba(0,0,0,0.2); }
        .btn-primary { background: #006400; color: white; }
        .btn-success { background: #28a745; color: white; }
        .btn-danger { background: #dc3545; color: white; }
        .btn-warning { background: #ffc107; color: black; }
        .btn-info { background: #17a2b8; color: white; }
        .btn-secondary { background: #6c757d; color: white; }
        .btn-light { background: #f8f9fa; color: #212529; border: 1px solid #dee2e6; }
        .database-section { margin: 30px 0; padding: 25px; border-radius: 12px; background: #f8f9fa; box-shadow: 0 3px 15px rgba(0,0,0,0.08); }
        .print-section-btn { background: #006400; color: white; border: none; border-radius: 5px; padding: 8px 15px; cursor: pointer; display: flex; align-items: center; gap: 5px; font-size: 0.9em; }
        .print-section-btn:hover { background: #228B22; }
        .footer { text-align: center; padding: 30px; color: white; margin-top: 40px; border-radius: 15px; background: linear-gradient(135deg, #2c3e50 0%, #34495e 100%); }
        .mrsa-highlight { background-color: #ffe6e6 !important; border-left: 3px solid #8B0000 !important; }
        .sort-icon { margin-left: 5px; font-size: 0.8em; opacity: 0.6; }
        .stat-card { color: white; padding: 18px; border-radius: 10px; text-align: center; box-shadow: 0 4px 15px rgba(0,0,0,0.15); transition: transform 0.2s; }
        .stat-card:hover { transform: translateY(-3px); }
        .stat-card .stat-value { font-size: 1.9em; font-weight: bold; margin-bottom: 4px; }
        .stat-card .stat-label { font-size: 0.85em; opacity: 0.95; text-transform: uppercase; letter-spacing: 0.5px; }
        .stats-grid { display: grid; grid-template-columns: repeat(auto-fit, minmax(180px, 1fr)); gap: 15px; margin: 20px 0; }
        .typing-badge { display: inline-block; padding: 3px 10px; border-radius: 12px; font-size: 0.8em; font-weight: 600; background: #e0e0e0; color: #333; border: 1px solid #ccc; }
        .typing-badge.badge-mrsa { background: #8B0000; color: white; border-color: #8B0000; }
        .typing-badge.badge-mssa { background: #4682B4; color: white; border-color: #4682B4; }
        .typing-badge.agr-I { background: #16a34a; color: white; border-color: #16a34a; }
        .typing-badge.agr-II { background: #2563eb; color: white; border-color: #2563eb; }
        .typing-badge.agr-III { background: #f59e0b; color: white; border-color: #f59e0b; }
        .typing-badge.agr-IV { background: #dc2626; color: white; border-color: #dc2626; }
        .typing-badge.agr-NA { background: #6b7280; color: white; border-color: #6b7280; }
        .accordion { margin: 20px 0; }
        .accordion-item { background: #f8f9fa; border: 1px solid #dee2e6; margin-bottom: 10px; border-radius: 8px; overflow: hidden; }
        .accordion-header { background: #e9ecef; padding: 12px 20px; cursor: pointer; font-weight: bold; color: #1e3a8a; display: flex; justify-content: space-between; align-items: center; }
        .accordion-header:hover { background: #dee2e6; }
        .accordion-content { padding: 15px 20px; border-top: 1px solid #dee2e6; background: white; }
        .copy-btn { background: #6b7280; color: white; border: none; padding: 4px 14px; border-radius: 16px; cursor: pointer; font-size: 0.82em; font-weight: 600; transition: background 0.2s; }
        .copy-btn:hover { background: #4b5563; }

        /* Sample-centric isolate boxes */
        .isolate-box { border: 1px solid #ddd; border-radius: 12px; margin-bottom: 20px; padding: 20px; background: #fafafa; box-shadow: 0 2px 8px rgba(0,0,0,0.06); transition: box-shadow 0.2s; }
        .isolate-box:hover { box-shadow: 0 4px 14px rgba(0,0,0,0.10); }
        .isolate-box .sample-header { display: flex; align-items: center; gap: 15px; flex-wrap: wrap; }
        .isolate-box .sample-header h3 { font-size: 1.4em; margin: 0; }
        .isolate-box .sample-header .total-badge { background: #006400; color: white; padding: 4px 16px; border-radius: 20px; font-weight: bold; font-size: 0.9em; }
        .isolate-box .sample-header .typing-info { display: flex; flex-wrap: wrap; gap: 8px; }
        .isolate-box .toggle-btn { background: #006400; color: white; border: none; padding: 8px 16px; border-radius: 8px; cursor: pointer; font-weight: 600; font-size: 0.85em; display: inline-flex; align-items: center; gap: 6px; margin-left: auto; transition: background 0.2s; }
        .isolate-box .toggle-btn:hover { background: #228B22; }
        .isolate-box .box-details { margin-top: 15px; }
        .database-table-wrapper { margin: 15px 0; overflow-x: auto; border: 1px solid #e0e0e0; border-radius: 8px; }
        .database-table-wrapper table { width: 100%; border-collapse: collapse; font-size: 0.85em; min-width: 800px; }
        .database-table-wrapper table th { background: #2c3e50; color: white; padding: 8px 12px; text-align: left; white-space: nowrap; }
        .database-table-wrapper table td { padding: 8px 12px; border-bottom: 1px solid #e0e0e0; white-space: nowrap; }
        .database-table-wrapper table tr:hover { background: #f1f1f1; }
        .db-title { font-weight: bold; color: #006400; margin: 10px 0 5px 0; font-size: 1.05em; border-left: 4px solid #006400; padding-left: 10px; }
        .filter-controls { display: flex; flex-wrap: wrap; gap: 10px; align-items: center; background: #f8f9fa; padding: 15px; border-radius: 8px; margin-bottom: 20px; }
        .filter-controls select { padding: 10px; border-radius: 8px; border: 2px solid #ddd; background: white; min-width: 150px; }
        .results-counter { font-size: 0.9em; color: #555; font-weight: 600; padding: 8px 12px; background: #f0f0f0; border-radius: 6px; white-space: nowrap; }

        @media print { body * { visibility: hidden; } .tab-content.active, .tab-content.active * { visibility: visible; } .tab-content.active { position: absolute; left: 0; top: 0; width: 100%; padding: 20px; box-shadow: none; border-radius: 0; } .print-section-btn, .tab-navigation, .dashboard-grid, .search-box, .action-buttons, .filter-controls, .toggle-btn { display: none !important; } .isolate-box .box-details { display: block !important; } .data-table { page-break-inside: auto; } .data-table tr { page-break-inside: avoid; } }
        @media (max-width: 768px) { body { min-width: auto; overflow-x: auto; } .container { padding: 10px; } .main-header h1 { font-size: 2em; } .tab-button { padding: 8px 12px; font-size: 0.8em; } .dashboard-grid { grid-template-columns: repeat(auto-fit, minmax(180px, 1fr)); } .data-table { font-size: 0.8em; } }
        </style>
        """

    # -------------------------------------------------------------------------
    # JavaScript
    # -------------------------------------------------------------------------
    def _get_js(self, typing_json: str) -> str:
        return f"""
        <script>
        var sampleTyping = {typing_json};
        window.STAPHSCOPE_BOX_DATA = window.STAPHSCOPE_BOX_DATA || {{}};

        function switchTab(tabName) {{
            document.querySelectorAll('.tab-content').forEach(t => t.classList.remove('active'));
            document.querySelectorAll('.tab-button').forEach(b => b.classList.remove('active'));
            var content = document.getElementById(tabName + '-tab');
            var button = document.querySelector('.tab-button.' + tabName);
            if (content) content.classList.add('active');
            if (button) button.classList.add('active');
            if (event && event.currentTarget) event.currentTarget.classList.add('active');
            window.location.hash = tabName;
        }}

        function escapeHtml(str) {{
            return String(str)
                .replace(/&/g, '&amp;').replace(/</g, '&lt;')
                .replace(/>/g, '&gt;').replace(/"/g, '&quot;');
        }}

        function searchTable(tableId, searchId) {{
            var input = document.getElementById(searchId);
            if (!input) return;
            var filter = input.value.toUpperCase();
            var table = document.getElementById(tableId);
            if (!table || !table.tBodies[0]) return;
            var rows = table.tBodies[0].rows;
            for (var i = 0; i < rows.length; i++) {{
                var cells = rows[i].getElementsByTagName('td');
                var found = false;
                for (var j = 0; j < cells.length; j++) {{
                    if (cells[j] && (cells[j].textContent || cells[j].innerText).toUpperCase().indexOf(filter) > -1) {{
                        found = true; break;
                    }}
                }}
                rows[i].style.display = found ? '' : 'none';
            }}
        }}

        function highlightGenome(tableId, searchId) {{
            var el = document.getElementById(searchId);
            if (!el) return;
            var filter = el.value.toUpperCase().trim();
            var table = document.getElementById(tableId);
            if (!table) return;
            table.querySelectorAll('.genome-tag').forEach(function(t) {{
                t.classList.remove('highlight');
                if (filter && t.textContent.toUpperCase().indexOf(filter) > -1) t.classList.add('highlight');
            }});
        }}

        // ---------- Lazy box rendering (SIGILL fix) ----------
        function toggleBoxDetails(tabId, sample) {{
            var box = document.querySelector('#' + tabId + '-tab .isolate-box[data-sample="' + CSS.escape(sample) + '"]');
            if (!box) return;
            var details = box.querySelector('.box-details');
            var btn = box.querySelector('.toggle-btn');
            var isOpen = details.style.display === 'block';
            if (isOpen) {{
                details.style.display = 'none';
                if (btn) btn.innerHTML = '<i class="fas fa-chevron-down"></i> Show Details';
                return;
            }}
            if (!details.innerHTML.trim()) {{
                details.innerHTML = (tabId === 'mutation')
                    ? buildMutationTable(sample)
                    : buildBoxContent(tabId, sample);
            }}
            details.style.display = 'block';
            if (btn) btn.innerHTML = '<i class="fas fa-chevron-up"></i> Hide Details';
        }}

        function buildBoxContent(tabId, sample) {{
            var root = window.STAPHSCOPE_BOX_DATA || {{}};
            var byTab = root[tabId] || {{}};
            var entry = byTab[sample];
            if (!entry) return '<p>No data available</p>';
            var html = '';
            if (entry.amrfinder && entry.amrfinder.length) {{
                html += buildDBTable('AMRfinder', entry.amrfinder, 'amrfinder');
            }}
            if (entry.abricate) {{
                for (var db in entry.abricate) {{
                    var genes = entry.abricate[db];
                    if (genes && genes.length) {{
                        html += buildDBTable(db.toUpperCase(), genes, db);
                    }}
                }}
            }}
            return html || '<p>No gene details</p>';
        }}

        function buildDBTable(dbName, genes, dbKey) {{
            if (!genes || !genes.length) return '';
            var keys = Object.keys(genes[0]);
            var priority = ['gene', 'product', 'coverage_percent', 'identity_percent', 'accession', 'contig', 'start', 'stop', 'class', 'subclass', 'scope', 'resistance'];
            var ordered = priority.filter(function(k) {{ return keys.indexOf(k) !== -1; }});
            keys.forEach(function(k) {{ if (ordered.indexOf(k) === -1) ordered.push(k); }});
            var html = '<div class="database-table-wrapper" data-db="' + dbKey + '">';
            html += '<div class="db-title">' + escapeHtml(dbName) + '</div>';
            html += '<table><thead><tr>';
            ordered.forEach(function(c) {{
                var display = c.replace(/_/g, ' ').replace(/\\b\\w/g, function(x) {{ return x.toUpperCase(); }});
                html += '<th>' + escapeHtml(display) + '</th>';
            }});
            html += '</tr></thead><tbody>';
            genes.forEach(function(g) {{
                html += '<tr>';
                ordered.forEach(function(c) {{
                    var v = (g[c] == null) ? '' : g[c];
                    html += '<td>' + escapeHtml(v) + '</td>';
                }});
                html += '</tr>';
            }});
            html += '</tbody></table></div>';
            return html;
        }}

        function buildMutationTable(sample) {{
            var root = window.STAPHSCOPE_BOX_DATA || {{}};
            var byTab = root['mutation'] || {{}};
            var entry = byTab[sample];
            if (!entry || !entry.mutations || !entry.mutations.length) return '<p>No mutations</p>';
            var cols = ['gene','mutation','class','subclass','contig','start','stop','strand','coverage','identity','accession'];
            var html = '<div class="database-table-wrapper" data-db="mutations">';
            html += '<div class="db-title">Mutations</div><table><thead><tr>';
            cols.forEach(function(c) {{
                var display = c.replace(/_/g, ' ').replace(/\\b\\w/g, function(x) {{ return x.toUpperCase(); }});
                html += '<th>' + display + '</th>';
            }});
            html += '</tr></thead><tbody>';
            entry.mutations.forEach(function(m) {{
                html += '<tr>';
                cols.forEach(function(c) {{
                    var v = (m[c] == null) ? '' : m[c];
                    html += '<td>' + escapeHtml(v) + '</td>';
                }});
                html += '</tr>';
            }});
            html += '</tbody></table></div>';
            return html;
        }}

        function expandAllBoxes(tabId) {{
            var boxes = document.querySelectorAll('#' + tabId + '-tab .isolate-box');
            var visibleBoxes = [];
            boxes.forEach(function(b) {{ if (b.style.display !== 'none') visibleBoxes.push(b); }});
            if (visibleBoxes.length > 100) {{
                if (!confirm('You have ' + visibleBoxes.length + ' visible samples. Expanding all may take a moment. Continue?')) return;
            }}
            var i = 0, batch = 20;
            function step() {{
                var end = Math.min(i + batch, visibleBoxes.length);
                for (; i < end; i++) {{
                    var box = visibleBoxes[i];
                    var sample = box.getAttribute('data-sample');
                    var details = box.querySelector('.box-details');
                    var btn = box.querySelector('.toggle-btn');
                    if (details && details.style.display !== 'block') {{
                        if (!details.innerHTML.trim()) {{
                            details.innerHTML = (tabId === 'mutation')
                                ? buildMutationTable(sample)
                                : buildBoxContent(tabId, sample);
                        }}
                        details.style.display = 'block';
                        if (btn) btn.innerHTML = '<i class="fas fa-chevron-up"></i> Hide Details';
                    }}
                }}
                if (i < visibleBoxes.length) setTimeout(step, 30);
            }}
            step();
        }}

        function collapseAllBoxes(tabId) {{
            document.querySelectorAll('#' + tabId + '-tab .isolate-box').forEach(function(box) {{
                var details = box.querySelector('.box-details');
                if (details) details.style.display = 'none';
                var btn = box.querySelector('.toggle-btn');
                if (btn) btn.innerHTML = '<i class="fas fa-chevron-down"></i> Show Details';
            }});
        }}

        function filterBoxes(tabId) {{
            var searchEl = document.getElementById('search-' + tabId);
            var dbEl = document.getElementById('dbFilter-' + tabId);
            var search = (searchEl ? searchEl.value : '').toUpperCase();
            var dbFilter = dbEl ? dbEl.value : 'all';
            var boxes = document.querySelectorAll('#' + tabId + '-tab .isolate-box');
            var visible = 0;
            boxes.forEach(function(box) {{
                var sample = box.getAttribute('data-sample') || '';
                var show = (!search || sample.toUpperCase().indexOf(search) !== -1);
                box.style.display = show ? '' : 'none';
                if (show) visible++;
                if (show) {{
                    var wrappers = box.querySelectorAll('.database-table-wrapper');
                    wrappers.forEach(function(w) {{
                        var dbName = w.getAttribute('data-db') || '';
                        var keep = (dbFilter === 'all' || dbName === dbFilter);
                        w.style.display = keep ? '' : 'none';
                        var title = w.previousElementSibling;
                        if (title && title.classList.contains('db-title')) {{
                            title.style.display = keep ? '' : 'none';
                        }}
                    }});
                }}
            }});
            var counter = document.getElementById('counter-' + tabId);
            if (counter) counter.textContent = visible + ' shown';
        }}

        function resetBoxFilters(tabId) {{
            var s = document.getElementById('search-' + tabId);
            var d = document.getElementById('dbFilter-' + tabId);
            if (s) s.value = '';
            if (d) d.value = 'all';
            filterBoxes(tabId);
        }}

        // ---------- General utilities ----------
        function sortTable(tableId, colIndex, type) {{
            var table = document.getElementById(tableId);
            if (!table || !table.tBodies[0]) return;
            var tbody = table.tBodies[0];
            var rows = Array.from(tbody.rows);
            var asc = table.getAttribute('data-sort-dir') !== 'asc';
            rows.sort(function(a, b) {{
                var av = a.cells[colIndex].innerText.trim();
                var bv = b.cells[colIndex].innerText.trim();
                if (type === 'number') {{
                    av = parseFloat(av.replace(/,/g, '')) || 0;
                    bv = parseFloat(bv.replace(/,/g, '')) || 0;
                    return asc ? av - bv : bv - av;
                }}
                return asc ? av.localeCompare(bv) : bv.localeCompare(av);
            }});
            tbody.append.apply(tbody, rows);
            table.setAttribute('data-sort-dir', asc ? 'asc' : 'desc');
        }}

        function printSection(sectionId) {{
            var content = document.getElementById(sectionId);
            if (!content) return;
            // Expand all boxes before printing
            content.querySelectorAll('.isolate-box').forEach(function(box) {{
                var tabId = box.getAttribute('data-tab');
                var sample = box.getAttribute('data-sample');
                var details = box.querySelector('.box-details');
                if (details && !details.innerHTML.trim()) {{
                    details.innerHTML = (tabId === 'mutation')
                        ? buildMutationTable(sample)
                        : buildBoxContent(tabId, sample);
                }}
            }});
            var w = window.open('', '_blank');
            var style = document.querySelector('style');
            w.document.write('<html><head><title>Print</title>');
            if (style) w.document.write('<style>' + style.textContent + '</style>');
            w.document.write('</head><body>' + content.innerHTML + '</body></html>');
            w.document.close();
            w.print();
        }}

        function exportTableToCSV(tableId, filename) {{
            var table = document.getElementById(tableId);
            if (!table) return;
            var rows = table.querySelectorAll('tr');
            var csv = [];
            for (var i = 0; i < rows.length; i++) {{
                var row = [], cols = rows[i].querySelectorAll('td, th');
                for (var j = 0; j < cols.length; j++) {{
                    row.push('"' + (cols[j].innerText || '').replace(/"/g, '""') + '"');
                }}
                csv.push(row.join(','));
            }}
            var blob = new Blob([csv.join('\\n')], {{ type: 'text/csv' }});
            var a = document.createElement('a');
            a.download = filename; a.href = URL.createObjectURL(blob);
            document.body.appendChild(a); a.click(); document.body.removeChild(a);
        }}

        document.addEventListener('DOMContentLoaded', function() {{
            var hash = window.location.hash.substring(1);
            var target = hash ? document.querySelector('.tab-button.' + hash) : document.querySelector('.tab-button');
            if (target) target.click();
            document.querySelectorAll('.data-table').forEach(function(table) {{
                var headers = table.querySelectorAll('th');
                headers.forEach(function(h, idx) {{
                    var type = h.getAttribute('data-sort') || 'string';
                    h.style.cursor = 'pointer';
                    h.addEventListener('click', function() {{ sortTable(table.id, idx, type); }});
                    var icon = document.createElement('span');
                    icon.className = 'sort-icon'; icon.innerHTML = '⇅';
                    h.appendChild(icon);
                }});
            }});
            document.querySelectorAll('.accordion-header').forEach(function(header) {{
                header.addEventListener('click', function() {{
                    var content = this.nextElementSibling;
                    content.style.display = content.style.display === 'block' ? 'none' : 'block';
                }});
            }});
            document.querySelectorAll('.copy-btn').forEach(function(b) {{
                b.addEventListener('click', function() {{
                    var c = this.getAttribute('data-citation') || '';
                    navigator.clipboard.writeText(c).then(() => {{
                        var t = this.innerHTML;
                        this.innerHTML = '✓ Copied!';
                        setTimeout(() => {{ this.innerHTML = t; }}, 2000);
                    }});
                }});
            }});
        }});
        </script>
        """

    # -------------------------------------------------------------------------
    # SUMMARY
    # -------------------------------------------------------------------------
    def _generate_summary_section(self, kwargs: Dict) -> str:
        samples_data = kwargs['samples_data']
        patterns = kwargs['patterns']
        gene_centric = kwargs['gene_centric']
        total = len(samples_data)
        total_amr = sum(len(g) for g in gene_centric.get('amr_databases', {}).values())
        total_vir = sum(len(g) for g in gene_centric.get('virulence_databases', {}).values())
        total_plasmids = sum(len(g) for g in gene_centric.get('plasmid_databases', {}).values())
        total_bacmet = sum(len(g) for g in gene_centric.get('bacmet_databases', {}).values())
        critical = len(patterns.get('high_risk_combinations', []))
        mrsa = sum(1 for s in samples_data.values()
                   if 'MRSA' in s.get('typing', {}).get('MRSA_Status', ''))
        mssa = sum(1 for s in samples_data.values()
                   if 'MSSA' in s.get('typing', {}).get('MRSA_Status', ''))
        agr_dist = patterns.get('agr_type_distribution', Counter())
        agr_str = ', '.join(f"{k}: {v}" for k, v in sorted(agr_dist.items())) or 'None'
        mutation_count = len(kwargs.get('integrated_data', {}).get('mutation_details', {}))

        return f'''
        <div class="alert-box alert-info">
            <i class="fas fa-info-circle fa-2x"></i>
            <div>
                <h3>📊 Hybrid Analysis Overview</h3>
                <p>This report analyses <strong>{total}</strong> <em>Staphylococcus aureus</em> genomes using a hybrid strategy:</p>
                <ul>
                    <li><strong>Gene-Centric</strong> for typing tabs (MLST, spa, SCCmec, MRSA) — each marker shown with all genomes that carry it.</li>
                    <li><strong>Sample-Centric, Lazy-Loaded</strong> for AMR, Virulence, BACMET, Plasmids, Mutations — each isolate in its own collapsed box; gene tables load on demand.</li>
                </ul>
                <p><strong>v2.0.0:</strong> single master TSV, lazy box rendering (browser-safe on 1000+ samples), enriched education boxes, clickable citation DOIs.</p>
            </div>
        </div>
        <div class="alert-box alert-success">
            <i class="fas fa-magic fa-2x"></i>
            <div>
                <h3>📘 How to Use the Isolate Boxes</h3>
                <ol>
                    <li>Open AMR, Virulence, BACMET, Plasmids, or Mutations.</li>
                    <li>Each isolate shows a <strong>header</strong> with sample name, total count, and typing badges.</li>
                    <li>Click <strong>Show Details</strong> inside a box to load its gene/mutation tables.</li>
                    <li>Filter by sample name (search) or by database (dropdown).</li>
                    <li>Use <strong>Expand All (visible)</strong> to load many boxes in batches — never crashes.</li>
                </ol>
            </div>
        </div>
        <h3><i class="fas fa-chart-bar"></i> Key Statistics</h3>
        <div class="scrollable-table"><table class="data-table">
            <thead><tr><th>Metric</th><th>Count</th><th>Details</th></tr></thead>
            <tbody>
                <tr><td>Total Samples Analysed</td><td><strong>{total}</strong></td><td>Complete genomic analysis</td></tr>
                <tr><td>MRSA Samples</td><td><span class="badge badge-mrsa">{mrsa}</span></td><td>Methicillin-resistant S. aureus</td></tr>
                <tr><td>MSSA Samples</td><td><span class="badge badge-mssa">{mssa}</span></td><td>Methicillin-sensitive S. aureus</td></tr>
                <tr><td>Unique MLST Types</td><td><strong>{len(patterns.get('mlst_distribution', {}))}</strong></td><td>Sequence types</td></tr>
                <tr><td>Unique spa Types</td><td><strong>{len(patterns.get('spa_type_distribution', {}))}</strong></td><td>Protein A typing</td></tr>
                <tr><td>Unique SCCmec (CGE)</td><td><strong>{len(patterns.get('sccmec_cge_distribution', {}))}</strong></td><td>SCCmec cassette (CGE)</td></tr>
                <tr><td>Unique SCCmec (RPet)</td><td><strong>{len(patterns.get('sccmec_rpet_distribution', {}))}</strong></td><td>SCCmec cassette (RPet)</td></tr>
                <tr><td>Unique SCCmec Subtypes</td><td><strong>{len(patterns.get('sccmec_subtype_distribution', {}))}</strong></td><td>Cassette subtypes</td></tr>
                <tr><td>Unique Capsule Types</td><td><strong>{len(patterns.get('capsule_distribution', {}))}</strong></td><td>Serotype / vaccine markers</td></tr>
                <tr><td>agr Types Detected</td><td><strong>{len(agr_dist)}</strong></td><td>Distribution: {agr_str}</td></tr>
                <tr><td>Samples with Mutations</td><td><strong>{mutation_count}</strong></td><td>Point mutations detected</td></tr>
                <tr><td>AMR Genes</td><td><strong>{total_amr}</strong></td><td>Across all AMR databases</td></tr>
                <tr><td>Virulence Genes</td><td><strong>{total_vir}</strong></td><td>From VFDB</td></tr>
                <tr><td>BACMET Genes</td><td><strong>{total_bacmet}</strong></td><td>Biocide and heavy-metal resistance</td></tr>
                <tr><td>Plasmid Replicons</td><td><strong>{total_plasmids}</strong></td><td>Plasmid families</td></tr>
                <tr><td>High-Risk AMR+Virulence</td><td><span class="badge badge-critical">{critical}</span></td><td>Samples with both critical classes</td></tr>
            </tbody>
        </table></div>'''

    # -------------------------------------------------------------------------
    # SAMPLE OVERVIEW (master TSV, capsule coloring)
    # -------------------------------------------------------------------------
    def _generate_sample_overview_section(self, kwargs: Dict) -> str:
        samples_data = kwargs['samples_data']
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

            row_class = 'class="mrsa-highlight"' if 'MRSA' in mrsa else ''
            status_badge = ('<span class="badge badge-mrsa">MRSA</span>' if 'MRSA' in mrsa
                            else ('<span class="badge badge-mssa">MSSA</span>' if 'MSSA' in mrsa
                                  else esc(mrsa)))
            agr_class = f"agr-{agr}" if agr in ('I', 'II', 'III', 'IV') else 'agr-NA'
            cap_html = self._colorize_capsule_cell(cap)

            vir_tags = ''.join(f'<span class="genome-tag">{esc(g)}</span>' for g in vir_genes)
            vir_details = (f'<details class="vir-details"><summary>{vir_count} gene(s)</summary>'
                           f'<div style="margin-top:6px;display:flex;flex-wrap:wrap;gap:4px;">'
                           f'{vir_tags or "<em>None</em>"}</div></details>')

            rows += f'''<tr {row_class}>
                <td><strong>{esc(sample)}</strong></td>
                <td>{esc(mlst)}</td>
                <td>{esc(spa)}</td>
                <td><span class="typing-badge {agr_class}">{esc(agr)}</span></td>
                <td>{cap_html}</td>
                <td>{esc(cge)}</td>
                <td>{esc(rpet)}</td>
                <td>{esc(sub)}</td>
                <td>{status_badge}</td>
                <td>{vir_details}</td>
            </tr>'''

        return f'''
        <div class="alert-box alert-info">
            <i class="fas fa-info-circle fa-2x"></i>
            <div>
                <h3>🧬 Sample Overview – The Population Snapshot</h3>
                <p>This is your <strong>master reference table</strong> — every isolate on one row, showing all the typing layers side by side. Use it to spot <strong>clonal clusters</strong>, track <strong>resistance-linked lineages</strong>, and flag <strong>outbreak candidates</strong> at a glance.</p>
                <ul style="margin-top:8px;">
                    <li><strong>MLST</strong> — global sequence type. The lingua franca for comparing your isolates to worldwide <em>S. aureus</em> populations.</li>
                    <li><strong>spa Type</strong> — fine-resolution outbreak tracking within a single ST.</li>
                    <li><strong>agr Type</strong> — quorum-sensing regulator (I–IV); influences virulence and biofilm behaviour.</li>
                    <li><strong>Capsule Type</strong> — serotype marker (Type 5, Type 8); informs vaccine coverage and immune-evasion potential.</li>
                    <li><strong>SCCmec (CGE / RPet / Subtype)</strong> — the <em>mec</em> cassette that defines MRSA lineage; two callers shown for cross-validation.</li>
                    <li><strong>MRSA / MSSA</strong> — the clinical bottom line for β-lactam therapy choice. MRSA rows are highlighted red.</li>
                    <li><strong>Virulence</strong> — click the count to expand the full VFDB gene list inline.</li>
                </ul>
                <p style="margin-top:8px;"><strong>Tip:</strong> Click any column header to sort. Combine <em>MLST + spa + SCCmec subtype</em> to identify probable transmission clusters.</p>
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
                </tr></thead>
                <tbody>{rows}</tbody>
            </table>
        </div>'''

    # -------------------------------------------------------------------------
    # FASTA QC (Biopython + fastANI credits)
    # -------------------------------------------------------------------------
    def _generate_qc_section(self, kwargs: Dict) -> str:
        qc_data = kwargs.get('integrated_data', {}).get('qc_data', {})
        if not qc_data:
            return '''
            <div class="alert-box alert-warning">
                <i class="fas fa-exclamation-circle fa-2x"></i>
                <div><h3>No QC Data Available</h3>
                <p>The FASTA_QC_summary.html file was not found or could not be parsed.</p></div>
            </div>'''
        all_metrics = set()
        for m in qc_data.values():
            all_metrics.update(m.keys())
        metric_list = sorted(all_metrics)

        credit = self._credit_bar('#17a2b8', '📊', 'FASTA QC Metrics',
            'Computed with <strong>Biopython</strong> — assembly statistics '
            '(contigs, N50, GC%, total length) from FASTA files.')

        fastani_credit = self._credit_bar('#17a2b8', '🧬',
            'Species Confirmation – fastANI',
            '<strong>fastANI</strong> developed by '
            '<a href="https://github.com/ParBLiSS/FastANI" target="_blank" '
            'style="color:#17a2b8;font-weight:bold;">ParBLiSS (Jain et al.)</a> — '
            'high-throughput Average Nucleotide Identity calculation for species-level identification.<br>'
            '<span style="font-size:0.9em;color:#6c757d;"><i class="fas fa-book-open"></i> '
            'Please cite: Jain C, Rodriguez-R LM, Phillippy AM, Konstantinidis KT, Aluru S. '
            'High throughput ANI analysis of 90K prokaryotic genomes reveals clear species boundaries. '
            '<em>Nat Commun</em>. 2018;9(1):5114. '
            '<a href="https://doi.org/10.1038/s41467-018-07641-9" target="_blank" '
            'style="color:#17a2b8;font-weight:bold;">🔗 DOI</a></span><br>'
            '<span style="font-size:0.9em;color:#6c757d;"><i class="fas fa-gratipay"></i> '
            'We are grateful to the developers for making this tool freely available.</span>')

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
                    <li><strong>Contigs</strong> – lower is better; typical S. aureus: &lt;200 (good), &lt;100 (excellent).</li>
                    <li><strong>N50</strong> – higher is better; &gt;50 kb (good), &gt;100 kb (excellent).</li>
                    <li><strong>GC%</strong> – S. aureus is typically 32–33%.</li>
                    <li><strong>Total length</strong> – ~2.8 Mbp for a complete genome.</li>
                    <li><strong>ANI ≥ 95%</strong> against an <em>S. aureus</em> reference confirms species identity (fastANI).</li>
                </ul>
            </div>
        </div>
        <input type="text" class="search-box" id="search-qc"
               onkeyup="searchTable('qc-table', 'search-qc')"
               placeholder="🔍 Search sample...">
        <div class="action-buttons">
            <button class="action-btn btn-primary" onclick="exportTableToCSV('qc-table', 'fasta_qc.csv')"><i class="fas fa-download"></i> Export QC Data</button>
        </div>
        <div class="master-scrollable-container">
            <table id="qc-table" class="data-table">
                <thead><tr><th data-sort="string">Sample</th>{header_cells}</tr></thead>
                <tbody>{rows}</tbody>
            </table>
        </div>'''

    # -------------------------------------------------------------------------
    # MLST / SPA / SCCMEC / MRSA / AGR — rich acknowledgements
    # -------------------------------------------------------------------------
    def _generate_mlst_section(self, kwargs: Dict) -> str:
        patterns = kwargs['patterns']
        mlst_dist = patterns.get('mlst_distribution', Counter())
        mlst_spa = patterns.get('mlst_spa_combinations', {})
        mlst_scc = patterns.get('mlst_sccmec_combinations', {})

        credit = self._credit_bar('#FF9800', '🧬',
            'MLST Typing – Acknowledgments &amp; Licensing',
            '<strong>MLST scheme</strong> powered by '
            '<a href="https://github.com/tseemann/mlst" target="_blank" style="color:#FF9800;font-weight:bold;">Prof. Torsten Seemann’s Perl scripts</a> '
            'and the <a href="https://pubmlst.org/" target="_blank" style="color:#FF9800;font-weight:bold;">PubMLST database</a> '
            '(Jolley et al., <em>Wellcome Open Res</em> 2018).<br>'
            '<span style="color:#856404;"><i class="fas fa-info-circle"></i> '
            '<strong>Note:</strong> Allele definitions current as of <strong>2024</strong>. Future updates may require manual downloads.</span><br>'
            '<span style="font-size:0.9em;color:#6c757d;"><i class="fas fa-gratipay"></i> '
            'We thank the PubMLST curators and Torsten Seemann for their invaluable work.</span>')

        rows = ''
        total = sum(mlst_dist.values())
        for mlst, count in mlst_dist.most_common():
            if mlst in ('ND', 'Not Assigned'):
                continue
            pct = (count / total * 100) if total else 0
            spas = [c.split(' - ')[1] for c in mlst_spa if c.startswith(f"{mlst} - ")]
            sccs = [c.split(' - ')[1] for c in mlst_scc if c.startswith(f"{mlst} - ")]
            rows += (f'<tr><td><strong>{esc(mlst)}</strong></td><td>{count}</td><td>{pct:.1f}%</td>'
                     f'<td>{esc(", ".join(sorted(set(spas))) or "ND")}</td>'
                     f'<td>{esc(", ".join(sorted(set(sccs))) or "ND")}</td></tr>')

        combo_rows = ''
        for combo, samples in sorted(mlst_spa.items(), key=lambda x: -len(x[1])):
            tags = ''.join(f'<span class="genome-tag">{esc(s)}</span>' for s in samples)
            combo_rows += (f'<tr><td><strong>{esc(combo)}</strong></td><td>{len(samples)}</td>'
                           f'<td><div class="genome-list">{tags}</div></td></tr>')

        return f'''
        {credit}
        {self._alert('info', 'fa-code-branch',
            '<h3>🔬 MLST (Multi-Locus Sequence Typing)</h3>'
            '<p>Seven housekeeping genes; each unique allele combination defines a Sequence Type (ST) — the gold standard for global <em>S. aureus</em> epidemiology.</p>'
            f'<p><strong>{len(mlst_dist)} unique STs</strong> identified.</p>')}
        <h3>📊 ST Distribution</h3>
        <div class="scrollable-table"><table class="data-table">
            <thead><tr><th>ST</th><th>Count</th><th>%</th><th>Associated spa Types</th><th>Associated SCCmec Types</th></tr></thead>
            <tbody>{rows}</tbody>
        </table></div>
        <h3>🔗 ST – spa Combinations</h3>
        <input type="text" class="search-box" id="search-mlst-spa" onkeyup="searchTable('mlst-spa-table','search-mlst-spa')" placeholder="🔍 Search ST-spa...">
        <div class="master-scrollable-container"><table id="mlst-spa-table" class="data-table">
            <thead><tr><th>ST-spa Combination</th><th>Count</th><th>Samples</th></tr></thead>
            <tbody>{combo_rows}</tbody>
        </table></div>'''

    def _generate_spa_section(self, kwargs: Dict) -> str:
        patterns = kwargs['patterns']
        spa_dist = patterns.get('spa_type_distribution', Counter())
        mlst_spa = patterns.get('mlst_spa_combinations', {})
        spa_scc = patterns.get('spa_sccmec_combinations', {})

        credit = self._credit_bar('#9C27B0', '🧬',
            'spa Typing – Acknowledgments',
            '<strong>spa typing</strong> powered by '
            '<a href="https://github.com/mjsull/spa_typing" target="_blank" style="color:#9C27B0;font-weight:bold;">original code by mjsull</a>, '
            'modified by <strong>JFSanchezHerrero</strong>, and the '
            '<a href="https://spa.ridom.de/" target="_blank" style="color:#9C27B0;font-weight:bold;">Ridom SpaServer database</a>.<br>'
            '<span style="font-size:0.9em;color:#6c757d;"><i class="fas fa-gratipay"></i> '
            'We thank the developers and curators for maintaining this essential resource.</span>')

        rows = ''
        total = sum(spa_dist.values())
        for spa, count in spa_dist.most_common():
            if spa in ('ND', 'Not Assigned'):
                continue
            pct = (count / total * 100) if total else 0
            sts = [c.split(' - ')[0] for c in mlst_spa if c.endswith(f" - {spa}")]
            rows += (f'<tr><td><strong>{esc(spa)}</strong></td><td>{count}</td><td>{pct:.1f}%</td>'
                     f'<td>{esc(", ".join(sorted(set(sts))) or "None")}</td></tr>')

        combo_rows = ''
        for combo, samples in sorted(spa_scc.items(), key=lambda x: -len(x[1])):
            tags = ''.join(f'<span class="genome-tag">{esc(s)}</span>' for s in samples)
            combo_rows += (f'<tr><td><strong>{esc(combo)}</strong></td><td>{len(samples)}</td>'
                           f'<td><div class="genome-list">{tags}</div></td></tr>')

        return f'''
        {credit}
        {self._alert('info', 'fa-dna',
            '<h3>🧬 spa Typing – High-Resolution Outbreak Tracking</h3>'
            '<p>The <em>spa</em> gene encodes protein A; repeat region polymorphisms define spa types.</p>'
            f'<p><strong>{len(spa_dist)} unique spa types</strong> identified.</p>')}
        <h3>📊 spa Type Distribution</h3>
        <div class="scrollable-table"><table class="data-table">
            <thead><tr><th>spa Type</th><th>Count</th><th>%</th><th>Common STs</th></tr></thead>
            <tbody>{rows}</tbody>
        </table></div>
        <h3>🔗 spa – SCCmec Combinations</h3>
        <input type="text" class="search-box" id="search-spa-scc" onkeyup="searchTable('spa-scc-table','search-spa-scc')" placeholder="🔍 Search spa-SCCmec...">
        <div class="master-scrollable-container"><table id="spa-scc-table" class="data-table">
            <thead><tr><th>spa-SCCmec Combination</th><th>Count</th><th>Samples</th></tr></thead>
            <tbody>{combo_rows}</tbody>
        </table></div>'''

    def _generate_sccmec_section(self, kwargs: Dict) -> str:
        patterns = kwargs['patterns']
        cge_dist = patterns.get('sccmec_cge_distribution', Counter())
        rpet_dist = patterns.get('sccmec_rpet_distribution', Counter())
        sub_dist = patterns.get('sccmec_subtype_distribution', Counter())

        cge_credit = self._credit_bar('#009688', '🛡️',
            'SCCmec Typing (CGE) – Acknowledgments',
            '<strong>SCCmecFinder</strong> by '
            '<a href="https://cge.cbs.dtu.dk/services/SCCmecFinder/" target="_blank" style="color:#009688;font-weight:bold;">Center for Genomic Epidemiology (DTU)</a>. '
            'Curated by <strong>Anders Rhod Larsen</strong>.<br>'
            '<span style="font-size:0.9em;color:#6c757d;"><i class="fas fa-book-open"></i> '
            'Cite: Kaya H, et al. <em>mSphere</em>. 2018;3(1):e00612-17. '
            '<a href="https://doi.org/10.1128/mSphere.00612-17" target="_blank" style="color:#009688;font-weight:bold;">🔗 DOI</a></span>')

        rpet_credit = self._credit_bar('#7c3aed', '🔬',
            'SCCmec Typing (RPet) – Acknowledgments',
            '<strong>sccmec</strong> developed by <strong>Robert A. Petit III, PhD</strong> — '
            '<a href="https://github.com/rpetit3/sccmec" target="_blank" style="color:#7c3aed;font-weight:bold;">github.com/rpetit3/sccmec</a>.<br>'
            '<span style="font-size:0.9em;color:#6c757d;"><i class="fas fa-book-open"></i> '
            'Cite: Petit RA III, Read TD. <em>PeerJ</em>. 2018;6:e5261. '
            '<a href="https://doi.org/10.7717/peerj.5261" target="_blank" style="color:#7c3aed;font-weight:bold;">🔗 DOI</a></span><br>'
            '<span style="font-size:0.9em;color:#6c757d;"><i class="fas fa-gratipay"></i> '
            'We thank Robert Petit for his sustained contributions to open-source <em>S. aureus</em> genomics.</span>')

        def dist_table(tid, dist, label):
            total = sum(dist.values())
            rows = ''
            for v, c in dist.most_common():
                if v in ('Not Assigned', 'ND', ''):
                    continue
                pct = (c / total * 100) if total else 0
                rows += f'<tr><td><strong>{esc(v)}</strong></td><td>{c}</td><td>{pct:.1f}%</td></tr>'
            return f'''
            <div class="scrollable-table"><table id="{tid}" class="data-table">
                <thead><tr><th>{label}</th><th>Count</th><th>%</th></tr></thead>
                <tbody>{rows}</tbody>
            </table></div>'''

        return f'''
        {cge_credit}
        {self._alert('info', 'fa-shield-alt',
            '<h3>🛡️ SCCmec – The MRSA Cassette</h3>'
            '<p>Three sub-sections: <strong>CGE caller</strong>, <strong>RPet caller</strong>, and the fine-grained <strong>Subtype</strong>.</p>')}
        <h3>📊 SCCmec Type Distribution — CGE caller</h3>
        {dist_table('scc-cge-dist', cge_dist, 'SCCmec (CGE)')}
        <h3 style="margin-top:30px;">📊 SCCmec Type Distribution — RPet caller</h3>
        {rpet_credit}
        {dist_table('scc-rpet-dist', rpet_dist, 'SCCmec (RPet)')}
        <h3 style="margin-top:30px;">📊 SCCmec Subtype Distribution</h3>
        {dist_table('scc-sub-dist', sub_dist, 'SCCmec Subtype')}'''

    def _generate_mrsa_section(self, kwargs: Dict) -> str:
        patterns = kwargs['patterns']
        samples_data = kwargs['samples_data']
        mrsa_status = patterns.get('mrsa_status_distribution', Counter())
        mrsa_samples = [s for s, d in samples_data.items()
                        if 'MRSA' in d.get('typing', {}).get('MRSA_Status', '')]

        mrsa_mlst_spa = defaultdict(list)
        mrsa_mlst_scc = defaultdict(list)
        for s in mrsa_samples:
            t = samples_data[s].get('typing', {})
            mlst = t.get('MLST', 'Not Assigned')
            spa = t.get('spa_Type', 'Not Assigned')
            cge = t.get('SCCmec_CGE', 'Not Assigned')
            if mlst != 'Not Assigned' and spa != 'Not Assigned':
                mrsa_mlst_spa[f"{mlst} - {spa}"].append(s)
            if mlst != 'Not Assigned' and cge != 'Not Assigned':
                mrsa_mlst_scc[f"{mlst} - {cge}"].append(s)

        def combo_block(tid, title, dict_):
            if not dict_:
                return ''
            rows = ''
            for combo, samples in sorted(dict_.items(), key=lambda x: -len(x[1])):
                tags = ''.join(f'<span class="genome-tag">{esc(s)}</span>' for s in samples)
                rows += (f'<tr><td><strong>{esc(combo)}</strong></td><td>{len(samples)}</td>'
                         f'<td><div class="genome-list">{tags}</div></td></tr>')
            return f'''
            <h3>🔗 {title}</h3>
            <input type="text" class="search-box" id="search-{tid}" onkeyup="searchTable('{tid}','search-{tid}')" placeholder="🔍 Search...">
            <div class="master-scrollable-container"><table id="{tid}" class="data-table">
                <thead><tr><th>{title}</th><th>Count</th><th>Samples</th></tr></thead>
                <tbody>{rows}</tbody>
            </table></div>'''

        status_rows = ''
        for status, count in mrsa_status.most_common():
            if status in ('Not Assigned', 'ND', ''):
                continue
            badge = ('<span class="badge badge-mrsa">MRSA</span>' if 'MRSA' in status
                     else '<span class="badge badge-mssa">MSSA</span>')
            status_rows += f'<tr><td>{badge}</td><td>{count}</td></tr>'

        return f'''
        {self._alert('danger', 'fa-skull-crossbones',
            f'<h3>⚠️ MRSA – A Clinical Priority</h3>'
            f'<p><strong>{len(mrsa_samples)} MRSA samples</strong> identified.</p>')}
        <h3>📊 MRSA vs MSSA</h3>
        <div class="scrollable-table"><table class="data-table">
            <thead><tr><th>Status</th><th>Count</th></tr></thead>
            <tbody>{status_rows}</tbody>
        </table></div>
        {combo_block('mrsa-mlst-spa', 'MRSA: ST – spa', mrsa_mlst_spa)}
        {combo_block('mrsa-mlst-scc', 'MRSA: ST – SCCmec (CGE)', mrsa_mlst_scc)}'''

    def _generate_agr_section(self, kwargs: Dict) -> str:
        patterns = kwargs['patterns']
        samples_data = kwargs['samples_data']
        agr_dist = patterns.get('agr_type_distribution', Counter())

        credit = self._credit_bar('#8B5CF6', '🧬',
            'agr Typing – Acknowledgments',
            '<strong>AgrVATE</strong> by '
            '<a href="https://github.com/VishnuRaghuram94/AgrVATE" target="_blank" style="color:#8B5CF6;font-weight:bold;">Vishnu Raghuram</a>, '
            'maintained by <strong>Robert A. Petit III</strong> '
            '(<a href="https://github.com/rpetit3" target="_blank" style="color:#8B5CF6;">@rpetit3</a>).<br>'
            '<span style="font-size:0.9em;color:#6c757d;"><i class="fas fa-book-open"></i> '
            'Cite: Raghuram V, et al. <em>Microbiol Spectr</em>. 2022;10(1):e0133421. '
            '<a href="https://doi.org/10.1128/spectrum.01334-21" target="_blank" style="color:#8B5CF6;font-weight:bold;">🔗 DOI</a></span><br>'
            '<span style="font-size:0.9em;color:#6c757d;"><i class="fas fa-gratipay"></i> '
            'We thank the developers for making this open-source tool available.</span>')

        dist_rows = ''
        total = sum(agr_dist.values())
        for typ in ['I', 'II', 'III', 'IV']:
            count = agr_dist.get(typ, 0)
            pct = (count / total * 100) if total else 0
            dist_rows += (f'<tr><td><span class="typing-badge agr-{typ}">{typ}</span></td>'
                          f'<td>{count}</td><td>{pct:.1f}%</td></tr>')

        samples_by_agr = defaultdict(list)
        for s, d in samples_data.items():
            a = d.get('typing', {}).get('agr_Type', 'Not Assigned')
            if a in ('I', 'II', 'III', 'IV'):
                samples_by_agr[a].append(s)

        sample_rows = ''
        for typ in ['I', 'II', 'III', 'IV']:
            samps = samples_by_agr.get(typ, [])
            if not samps:
                continue
            tags = ''.join(f'<span class="genome-tag">{esc(s)}</span>' for s in sorted(samps))
            sample_rows += (f'<tr><td><span class="typing-badge agr-{typ}">{typ}</span></td>'
                            f'<td>{len(samps)}</td>'
                            f'<td><div class="genome-list">{tags}</div></td></tr>')

        return f'''
        {credit}
        {self._alert('info', 'fa-dna',
            '<h3>🧬 agr Typing – Virulence Regulation</h3>'
            '<p>The accessory gene regulator (<em>agr</em>) system is a quorum-sensing circuit controlling virulence gene expression. Four types (I–IV).</p>')}
        <h3>📊 agr Type Distribution</h3>
        <div class="scrollable-table"><table class="data-table">
            <thead><tr><th>agr Type</th><th>Count</th><th>%</th></tr></thead>
            <tbody>{dist_rows}</tbody>
        </table></div>
        <h3>📋 Samples by agr Type</h3>
        <div class="master-scrollable-container"><table class="data-table">
            <thead><tr><th>agr Type</th><th>Count</th><th>Samples</th></tr></thead>
            <tbody>{sample_rows}</tbody>
        </table></div>'''

    # -------------------------------------------------------------------------
    # SAMPLE-CENTRIC BOXES (lazy rendering — SIGILL fix)
    # -------------------------------------------------------------------------
    def _generate_sample_centric_boxes(self, kwargs: Dict, tab_id: str, title: str, db_list: List[str]) -> str:
        """Lazy-rendered isolate boxes. Only headers are rendered server-side;
        gene tables are built on demand from an embedded JSON blob."""
        samples_data = kwargs.get('samples_data', {})
        amr_details = kwargs.get('integrated_data', {}).get('amrfinder_details', {})
        abricate_details = kwargs.get('integrated_data', {}).get('abricate_details', {})

        relevant_samples = []
        for sample in samples_data:
            has_data = False
            if 'amrfinder' in db_list and amr_details.get(sample):
                has_data = True
            for db in db_list:
                if db != 'amrfinder' and abricate_details.get(sample, {}).get(db):
                    has_data = True
            if has_data:
                relevant_samples.append(sample)
        relevant_samples.sort()

        # Educational blocks
        education_html = (
            self._multi_db_education() +
            self._confidence_tiers() +
            self._acquired_intrinsic() +
            self._genotype_phenotype_caveat()
        )

        if not relevant_samples:
            return f'''
            {education_html}
            <div class="alert-box alert-warning">
                <i class="fas fa-exclamation-circle fa-2x"></i>
                <div><h3>No {title} Data Available</h3>
                <p>No samples with {title} genes were found.</p></div>
            </div>'''

        # Build the JSON blob
        box_data = {}
        for sample in relevant_samples:
            entry = {'amrfinder': [], 'abricate': {}}
            if 'amrfinder' in db_list:
                entry['amrfinder'] = amr_details.get(sample, [])
            for db in db_list:
                if db != 'amrfinder':
                    entry['abricate'][db] = abricate_details.get(sample, {}).get(db, [])
            box_data[sample] = entry
        box_json = json.dumps(box_data, default=str, ensure_ascii=False)

        db_options = '<option value="all">All Databases</option>'
        for db in db_list:
            display = 'AMRfinder' if db == 'amrfinder' else db.upper()
            db_options += f'<option value="{db}">{display}</option>'

        html = f'''
        {education_html}
        <div class="alert-box alert-info">
            <i class="fas fa-info-circle fa-2x"></i>
            <div>
                <h3>🧬 {title} – Interactive Isolate Boxes (Lazy-Loaded)</h3>
                <p>Each box shows one isolate. Click <strong>Show Details</strong> to load its gene tables on demand — the report opens instantly even with thousands of samples.</p>
            </div>
        </div>
        <div class="filter-controls">
            <input type="text" class="search-box" id="search-{tab_id}" onkeyup="filterBoxes('{tab_id}')" placeholder="🔍 Search sample..." style="max-width:320px;">
            <select id="dbFilter-{tab_id}" onchange="filterBoxes('{tab_id}')">{db_options}</select>
            <button class="action-btn btn-info" onclick="expandAllBoxes('{tab_id}')"><i class="fas fa-expand-alt"></i> Expand All (visible)</button>
            <button class="action-btn btn-light" onclick="collapseAllBoxes('{tab_id}')"><i class="fas fa-compress-alt"></i> Collapse All</button>
            <button class="action-btn btn-success" onclick="resetBoxFilters('{tab_id}')"><i class="fas fa-sync"></i> Clear</button>
            <span class="results-counter" id="counter-{tab_id}">{len(relevant_samples)} shown</span>
        </div>
        <div id="box-container-{tab_id}">
        '''

        for sample in relevant_samples:
            sd = samples_data.get(sample, {})
            t = sd.get('typing', {})
            mlst = t.get('MLST', 'Not Assigned')
            spa = t.get('spa_Type', 'Not Assigned')
            cge = t.get('SCCmec_CGE', 'Not Assigned')
            mrsa = t.get('MRSA_Status', 'Not Assigned')
            agr_type = t.get('agr_Type', 'Not Assigned')

            mrsa_class = ('badge-mrsa' if 'MRSA' in mrsa
                          else 'badge-mssa' if 'MSSA' in mrsa else '')
            agr_class = f"agr-{agr_type}" if agr_type in ('I', 'II', 'III', 'IV') else 'agr-NA'

            total_genes = len(amr_details.get(sample, [])) if 'amrfinder' in db_list else 0
            for db in db_list:
                if db != 'amrfinder':
                    total_genes += len(abricate_details.get(sample, {}).get(db, []))
            cap_value = t.get('capsule_type', 'Not Assigned')

            html += f'''
            <div class="isolate-box" data-sample="{esc(sample)}" data-tab="{tab_id}">
                <div class="sample-header">
                    <h3><i class="fas fa-microbe"></i> {esc(sample)}</h3>
                    <span class="total-badge">Total Genes: {total_genes}</span>
                    <div class="typing-info">
                        <span class="typing-badge">ST: {esc(mlst)}</span>
                        <span class="typing-badge">spa: {esc(spa)}</span>
                        <span class="typing-badge">SCCmec: {esc(cge)}</span>
                        <span class="typing-badge {mrsa_class}">{esc(mrsa)}</span>
                        <span class="typing-badge {agr_class}">agr: {esc(agr_type)}</span>
                        {self._capsule_badge(cap_value)}
                    </div>
                    <button class="toggle-btn" onclick="toggleBoxDetails('{tab_id}', '{esc(sample)}')">
                        <i class="fas fa-chevron-down"></i> Show Details
                    </button>
                </div>
                <div class="box-details" style="display:none;"></div>
            </div>'''

        html += '</div>'
        html += f'''
        <script>
        window.STAPHSCOPE_BOX_DATA = window.STAPHSCOPE_BOX_DATA || {{}};
        window.STAPHSCOPE_BOX_DATA["{tab_id}"] = {box_json};
        </script>'''
        if tab_id == 'amr':
            html += self._db_roles_amr()
        return html

    # -------------------------------------------------------------------------
    # MUTATION BOXES (lazy rendering)
    # -------------------------------------------------------------------------
    def _generate_mutation_boxes(self, kwargs: Dict) -> str:
        """Sample-centric mutation view with lazy-rendered isolate boxes."""
        samples_data = kwargs.get('samples_data', {})
        integrated_data = kwargs.get('integrated_data', {})
        mutation_details = integrated_data.get('mutation_details', {})

        relevant_samples = [s for s in samples_data if mutation_details.get(s)]
        relevant_samples.sort()

        if not relevant_samples:
            return '''
            <div class="alert-box alert-warning">
                <i class="fas fa-exclamation-circle fa-2x"></i>
                <div><h3>No Mutation Data Available</h3>
                <p>No point mutations were found for any sample.</p></div>
            </div>'''

        box_data = {s: {'mutations': mutation_details.get(s, [])} for s in relevant_samples}
        box_json = json.dumps(box_data, default=str, ensure_ascii=False)

        credit = self._credit_bar('#00BCD4', '🧬',
            'AMRFinderPlus – Point Mutations',
            '<strong>AMRFinderPlus</strong> by '
            '<a href="https://github.com/ncbi/amr" target="_blank" '
            'style="color:#00BCD4;font-weight:bold;">NCBI</a> — the most comprehensive '
            'database for antimicrobial resistance genes and point mutations.<br>'
            '<span style="font-size:0.9em;color:#6c757d;"><i class="fas fa-book-open"></i> '
            'Cite: Feldgarden M, et al. <em>Sci Rep</em>. 2021;11(1):12728. '
            '<a href="https://doi.org/10.1038/s41598-021-91456-0" target="_blank" '
            'style="color:#00BCD4;font-weight:bold;">🔗 DOI</a></span><br>'
            '<span style="font-size:0.9em;color:#6c757d;"><i class="fas fa-gratipay"></i> '
            'We thank the NCBI team for maintaining this open resource.</span>')

        info_block = self._gene_family_info('#00BCD4',
            'Clinical relevance of key mutations:', [
                ('23S rRNA (linezolid)', 'Mutations (e.g., G2576T, T2500A) confer linezolid resistance.'),
                ('gyrA / parC (quinolones)', 'QRDR mutations reduce susceptibility to fluoroquinolones.'),
                ('rpoB (rifampin)', 'High-level rifampin resistance — combination therapy consideration.'),
                ('mprF (daptomycin)', 'Daptomycin non-susceptibility.'),
                ('rplC / rplD (linezolid)', 'Ribosomal protein mutations confer linezolid resistance.'),
                ('fusA (fusidic acid)', 'Fusidic acid resistance.'),
                ('mupA (mupirocin)', 'High-level mupirocin resistance.'),
            ])

        html = f'''
        {credit}
        {self._alert('info', 'fa-dna',
            '<h3>🧬 Point Mutations – Sample-Centric View (Lazy-Loaded)</h3>'
            '<p>Each box shows one isolate. Click <strong>Show Details</strong> to load its mutation table on demand.</p>')}
        {info_block}
        <div class="filter-controls">
            <input type="text" class="search-box" id="search-mutation"
                   onkeyup="filterBoxes('mutation')"
                   placeholder="🔍 Search sample..." style="max-width:320px;">
            <button class="action-btn btn-info" onclick="expandAllBoxes('mutation')">
                <i class="fas fa-expand-alt"></i> Expand All (visible)</button>
            <button class="action-btn btn-light" onclick="collapseAllBoxes('mutation')">
                <i class="fas fa-compress-alt"></i> Collapse All</button>
            <button class="action-btn btn-success" onclick="resetBoxFilters('mutation')">
                <i class="fas fa-sync"></i> Clear</button>
            <span class="results-counter" id="counter-mutation">{len(relevant_samples)} shown</span>
        </div>
        <div id="box-container-mutation">
        '''

        for sample in relevant_samples:
            sd = samples_data.get(sample, {})
            t = sd.get('typing', {})
            mlst = t.get('MLST', 'Not Assigned')
            spa = t.get('spa_Type', 'Not Assigned')
            cge = t.get('SCCmec_CGE', 'Not Assigned')
            mrsa = t.get('MRSA_Status', 'Not Assigned')
            agr_type = t.get('agr_Type', 'Not Assigned')

            mrsa_class = ('badge-mrsa' if 'MRSA' in mrsa
                          else 'badge-mssa' if 'MSSA' in mrsa else '')
            agr_class = f"agr-{agr_type}" if agr_type in ('I', 'II', 'III', 'IV') else 'agr-NA'

            n_mut = len(mutation_details.get(sample, []))
            cap_value = t.get('capsule_type', 'Not Assigned')

            html += f'''
            <div class="isolate-box" data-sample="{esc(sample)}" data-tab="mutation">
                <div class="sample-header">
                    <h3><i class="fas fa-microbe"></i> {esc(sample)}</h3>
                    <span class="total-badge">Total Mutations: {n_mut}</span>
                    <div class="typing-info">
                        <span class="typing-badge">ST: {esc(mlst)}</span>
                        <span class="typing-badge">spa: {esc(spa)}</span>
                        <span class="typing-badge">SCCmec: {esc(cge)}</span>
                        <span class="typing-badge {mrsa_class}">{esc(mrsa)}</span>
                        <span class="typing-badge {agr_class}">agr: {esc(agr_type)}</span>
                        {self._capsule_badge(cap_value)}
                    </div>
                    <button class="toggle-btn" onclick="toggleBoxDetails('mutation', '{esc(sample)}')">
                        <i class="fas fa-chevron-down"></i> Show Details
                    </button>
                </div>
                <div class="box-details" style="display:none;"></div>
            </div>'''

        html += '</div>'
        html += f'''
        <script>
        window.STAPHSCOPE_BOX_DATA = window.STAPHSCOPE_BOX_DATA || {{}};
        window.STAPHSCOPE_BOX_DATA["mutation"] = {box_json};
        </script>'''
        return html

    # -------------------------------------------------------------------------
    # PATTERN DISCOVERY
    # -------------------------------------------------------------------------
    def _generate_pattern_discovery_section(self, kwargs: Dict) -> str:
        patterns = kwargs['patterns']
        triple = patterns.get('triple_combinations', {})

        html = self._alert('info', 'fa-project-diagram',
            '<h3>🔍 Cross-Genome Pattern Discovery</h3>'
            '<p>Associations between typing results, gene co-occurrence, and high-risk combinations.</p>')

        # Triple typing
        triple_rows = ''
        for combo, samples in sorted(triple.items(), key=lambda x: -len(x[1])):
            tags = ''.join(f'<span class="genome-tag">{esc(s)}</span>' for s in samples)
            triple_rows += (f'<tr><td><strong>{esc(combo)}</strong></td><td>{len(samples)}</td>'
                            f'<td><div class="genome-list">{tags}</div></td></tr>')

        if triple_rows:
            html += f'''
            <h3>🔗 Triple Typing (ST – spa – SCCmec CGE)</h3>
            <input type="text" class="search-box" id="search-triple"
                   onkeyup="searchTable('triple-table','search-triple')"
                   placeholder="🔍 Search combination...">
            <div class="master-scrollable-container">
                <table id="triple-table" class="data-table">
                    <thead><tr><th>Combination</th><th>Count</th><th>Samples</th></tr></thead>
                    <tbody>{triple_rows}</tbody>
                </table>
            </div>'''

        # High-risk combos
        high_risk = patterns.get('high_risk_combinations', [])
        if high_risk:
            rows = ''
            for c in high_risk:
                rows += (f'<tr><td><strong>{esc(c["sample"])}</strong></td>'
                         f'<td>{esc(c["mlst"])}</td>'
                         f'<td>{esc(c["spa_type"])}</td>'
                         f'<td>{esc(c["sccmec_type"])}</td>'
                         f'<td>{esc(c["agr_type"])}</td>'
                         f'<td>{esc(", ".join(c["critical_amr_genes"]))}</td>'
                         f'<td>{esc(", ".join(c["critical_virulence_genes"]))}</td></tr>')
            html += f'''
            <h3>⚠️ High-Risk Combinations</h3>
            <div class="alert-box alert-danger">
                <i class="fas fa-radiation fa-2x"></i>
                <div><strong>{len(high_risk)} samples</strong> carry both critical AMR and virulence genes.</div>
            </div>
            <div class="master-scrollable-container">
                <table id="highrisk-table" class="data-table">
                    <thead><tr>
                        <th>Sample</th><th>MLST</th><th>spa</th><th>SCCmec</th>
                        <th>agr</th><th>Critical AMR</th><th>Critical Virulence</th>
                    </tr></thead>
                    <tbody>{rows}</tbody>
                </table>
            </div>'''

        # Gene co-occurrence (top 500)
        cooc = patterns.get('gene_cooccurrence', {})
        if cooc:
            pairs = []
            for g1, partners in cooc.items():
                for g2, cnt in partners.items():
                    pairs.append((g1, g2, cnt))
            pairs.sort(key=lambda x: -x[2])
            rows = ''.join(f'<tr><td>{esc(g1)}</td><td>{esc(g2)}</td><td>{cnt}</td></tr>'
                           for g1, g2, cnt in pairs[:500])
            html += f'''
            <h3>📈 Gene Co-occurrence (Top 500)</h3>
            <div class="master-scrollable-container">
                <table class="data-table">
                    <thead><tr><th>Gene 1</th><th>Gene 2</th><th>Co-occurrence</th></tr></thead>
                    <tbody>{rows}</tbody>
                </table>
            </div>'''
        return html

    # -------------------------------------------------------------------------
    # AI GUIDE
    # -------------------------------------------------------------------------
    def _generate_aiguide_section(self, kwargs: Dict) -> str:
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
                    <li><strong>Upload the JSON file</strong> — <code>staphscope_ultimate_sample_centric_report.json</code> contains all structured data. Best for precise quantitative queries.</li>
                    <li><strong>Upload the HTML report</strong> — modern AI tools parse HTML tables. Great for visual context.</li>
                    <li><strong>Copy-paste specific tables</strong> — quick insight without upload. Instant.</li>
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
                            <li>Does agr type correlate with MRSA status?</li>
                            <li>Which clones carry the most resistance genes?</li>
                        </ul>
                    </div>
                    <div style="background:#f8f9fa;padding:10px;border-radius:8px;">
                        <strong>💊 Antimicrobial Resistance</strong>
                        <ul style="margin-top:5px;font-size:.9em;">
                            <li>How many samples carry mecA? What are their STs and spa types?</li>
                            <li>Are there vanA/vanB positive samples?</li>
                            <li>Which AMR genes co-occur most frequently?</li>
                            <li>What is the distribution of tetracycline (tet) resistance?</li>
                            <li>Do any samples have combined β-lactam + macrolide resistance?</li>
                        </ul>
                    </div>
                    <div style="background:#f8f9fa;padding:10px;border-radius:8px;">
                        <strong>🦠 Virulence &amp; Toxins</strong>
                        <ul style="margin-top:5px;font-size:.9em;">
                            <li>Which samples carry PVL? Are they associated with specific STs or agr types?</li>
                            <li>List all samples with TSST-1.</li>
                            <li>Which enterotoxin genes are most prevalent?</li>
                            <li>Is there a correlation between biofilm (ica) genes and MRSA?</li>
                        </ul>
                    </div>
                    <div style="background:#f8f9fa;padding:10px;border-radius:8px;">
                        <strong>🧪 Mutations &amp; Biocides</strong>
                        <ul style="margin-top:5px;font-size:.9em;">
                            <li>What are the most frequent point mutations in gyrA or parC?</li>
                            <li>Are there any linezolid-related mutations (23S rRNA)?</li>
                            <li>Which samples carry qac genes (disinfectant resistance)?</li>
                            <li>Is mer resistance linked to specific STs?</li>
                        </ul>
                    </div>
                </div>
            </div>
            <div class="database-section">
                <h4><i class="fas fa-balance-scale"></i> Scientific Rigour &amp; Ethical AI Use</h4>
                <ul>
                    <li><strong>AI is your co-pilot, not the pilot.</strong> Interpret AI insights in context of local epidemiology and lab validation.</li>
                    <li><strong>Verify, verify, verify.</strong> Cross-check critical calls with primary literature or secondary tools.</li>
                    <li><strong>No patient-identifiable data.</strong> Only upload aggregated, de-identified genomic data.</li>
                    <li><strong>Transparency in publications.</strong> Mention AI-assisted pattern discovery in methods.</li>
                    <li><strong>AI hallucination is real.</strong> Treat every AI statement as a hypothesis, not a fact.</li>
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

    # -------------------------------------------------------------------------
    # CITATION (rich, clickable DOIs, 24-color palette)
    # -------------------------------------------------------------------------
    def _generate_citation_section(self, kwargs: Dict) -> str:
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
            ('MLST', 'Seemann T. MLST: Scan contig files against PubMLST typing schemes. GitHub. 2018.',
             'https://github.com/tseemann/mlst'),
            ('PubMLST / BIGSdb',
             'Jolley KA, Bray JE, Maiden MCJ. Open-access bacterial population genomics: BIGSdb software, the PubMLST.org website and their applications. <em>Wellcome Open Res</em>. 2018;3:124.',
             'https://doi.org/10.12688/wellcomeopenres.14826.1'),
            ('spa typing',
             'Harmsen D, et al. Typing of methicillin-resistant <em>Staphylococcus aureus</em> in a university hospital setting. <em>J Clin Microbiol</em>. 2003;41(12):5442-8.',
             'https://doi.org/10.1128/JCM.41.12.5442-5448.2003'),
            ('SCCmecFinder',
             'Kaya H, et al. SCCmecFinder, a Web-Based Tool for Typing of Staphylococcal Cassette Chromosome <em>mec</em> in <em>Staphylococcus aureus</em>. <em>mSphere</em>. 2018;3(1):e00612-17.',
             'https://doi.org/10.1128/mSphere.00612-17'),
            ('sccmec (RPet)',
             'Petit RA III, Read TD. <em>Staphylococcus aureus</em> viewed from the perspective of 40,000+ genomes. <em>PeerJ</em>. 2018;6:e5261.',
             'https://doi.org/10.7717/peerj.5261'),
            ('agrVATE',
             'Raghuram V, Alexander AM, Loo HQ, Petit RA 3rd, Goldberg JB, Read TD. Species-Wide Phylogenomics of the <em>Staphylococcus aureus</em> Agr Operon. <em>Microbiol Spectr</em>. 2022;10(1):e0133421.',
             'https://doi.org/10.1128/spectrum.01334-21'),
            ('Capsule Typing (cap5/cap8)',
             'Sau S, Bhasin N, Wann ER, Lee JC, Foster TJ, Lee CY. The <em>Staphylococcus aureus</em> allelic genetic loci for serotype 5 and 8 capsule expression. <em>Microbiology (Reading)</em>. 1997;143(Pt 7):2395-2405.',
             'https://doi.org/10.1099/00221287-143-7-2395'),
            ('fastANI',
             'Jain C, Rodriguez-R LM, Phillippy AM, Konstantinidis KT, Aluru S. High throughput ANI analysis of 90K prokaryotic genomes reveals clear species boundaries. <em>Nat Commun</em>. 2018;9(1):5114.',
             'https://doi.org/10.1038/s41467-018-07641-9'),
            ('AMRFinderPlus',
             'Feldgarden M, et al. AMRFinderPlus and the Reference Gene Catalog facilitate examination of the genomic links among antimicrobial resistance, stress response, and virulence. <em>Sci Rep</em>. 2021;11(1):12728.',
             'https://doi.org/10.1038/s41598-021-91456-0'),
            ('ABRicate', 'Seemann T. ABRicate: mass screening of contigs for antibiotic resistance genes. GitHub. 2024.',
             'https://github.com/tseemann/abricate'),
            ('CARD',
             'McArthur AG, et al. The comprehensive antibiotic resistance database. <em>Antimicrob Agents Chemother</em>. 2013;57(7):3348-57.',
             'https://doi.org/10.1128/AAC.00419-13'),
            ('ResFinder',
             'Florensa AF, et al. ResFinder – an open online resource for identification of antimicrobial resistance genes. <em>Microb Genom</em>. 2022;8(1):000748.',
             'https://doi.org/10.1099/mgen.0.000748'),
            ('VFDB',
             'Chen L, et al. VFDB 2012 update: toward the genetic diversity and molecular evolution of bacterial virulence factors. <em>Nucleic Acids Res</em>. 2012;40(Database issue):D641-5.',
             'https://doi.org/10.1093/nar/gkr989'),
            ('PlasmidFinder',
             'Carattoli A, et al. <em>In silico</em> detection and typing of plasmids using PlasmidFinder. <em>Antimicrob Agents Chemother</em>. 2014;58(7):3895-903.',
             'https://doi.org/10.1128/AAC.02412-14'),
            ('BacMet',
             'Pal C, et al. BacMet: antibacterial biocide and metal resistance genes database. <em>Nucleic Acids Res</em>. 2014;42(Database issue):D737-43.',
             'https://doi.org/10.1093/nar/gkt1252'),
            ('MEGARes',
             'Doster E, et al. MEGARes 2.0: a database for classification of antimicrobial drug, biocide and metal resistance determinants. <em>Nucleic Acids Res</em>. 2020;48(D1):D561-D569.',
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
        .citation-item {{ background:#fafbfc; border-radius:6px; margin-bottom:10px; padding:12px 14px; list-style:none; transition:box-shadow .2s,transform .2s; box-shadow:0 1px 3px rgba(0,0,0,.05); }}
        .citation-item:hover {{ box-shadow:0 4px 12px rgba(0,0,0,.10); transform:translateX(2px); }}
        .citation-body {{ display:flex; flex-direction:column; gap:8px; font-size:.92em; line-height:1.55; }}
        .citation-name {{ font-size:1em; font-weight:700; }}
        .citation-text {{ color:#333; }}
        .citation-actions {{ display:flex; gap:8px; flex-wrap:wrap; align-items:center; }}
        .citation-link {{ display:inline-flex; align-items:center; gap:4px; padding:4px 14px; border-radius:16px; font-size:.82em; font-weight:600; color:white; text-decoration:none; transition:opacity .2s,transform .2s; }}
        .citation-link:hover {{ opacity:.88; transform:translateY(-1px); }}
        .citation-actions .copy-btn {{ background:#6b7280; color:white; border:none; padding:4px 14px; border-radius:16px; cursor:pointer; font-size:.82em; font-weight:600; transition:background .2s; }}
        .citation-actions .copy-btn:hover {{ background:#4b5563; }}
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
                "Genomic analysis was performed using StaphScope [Beckley &amp; Amarh, 2026], which integrates MLST [Seemann, 2018] using the PubMLST database [Jolley et al., 2018], ABRicate [Seemann, 2018], AMRFinderPlus [Feldgarden et al., 2021], SCCmecFinder [Kaya et al., 2018], sccmec (RPet) [Petit &amp; Read, 2018], agrVATE [Raghuram et al., 2022], and fastANI [Jain et al., 2018] for comprehensive <em>S. aureus</em> characterization. Capsule typing used the cap5/cap8 locus reference [Sau et al., 1997]. Antimicrobial resistance genes were identified using the CARD [McArthur et al., 2013], ResFinder [Florensa et al., 2022], MEGARes [Doster et al., 2020], and ARG-ANNOT [Gupta et al., 2014] databases. For biocide and heavy metal resistance genes, BacMet [Pal et al., 2014] was used. Virulence and plasmid screening were performed with ABRicate using the VFDB [Chen et al., 2012] and PlasmidFinder [Carattoli et al., 2014] databases. Mutation detection was performed using AMRFinderPlus. Mobile genetic element profiling used mobileOG-db [Brown et al., 2022], Prodigal [Hyatt et al., 2010], and DIAMOND [Buchfink et al., 2015]. FASTA QC was performed using Biopython [Cock et al., 2009]."
            </div>
        </div>'''

    # -------------------------------------------------------------------------
    # FUNDING
    # -------------------------------------------------------------------------
    def _generate_funding_section(self, kwargs: Dict) -> str:
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

    # -------------------------------------------------------------------------
    # CALL TO ACTION
    # -------------------------------------------------------------------------
    def _calltoaction_section(self) -> str:
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
            <p>These bacteria "escape" the effects of antibiotics — hence the name. But we believe the name is also a global call to action:</p>
            <div style="background:#fff3e0;padding:15px;border-radius:8px;margin:15px 0;">
                <p><strong>🔹 E</strong>veryone must join forces — researchers, clinicians, policymakers, and citizens.<br>
                <strong>🔹 S</strong>mart surveillance is our first line of defence. No more guessing — we need genomic data.<br>
                <strong>🔹 K</strong>nowledge must be shared openly. No paywalls, no closed silos.<br>
                <strong>🔹 A</strong>frica bears a heavy AMR burden, but African solutions are already emerging.<br>
                <strong>🔹 P</strong>revention is cheaper than cure. Let's stop resistant infections before they spread.<br>
                <strong>🔹 E</strong>very day we delay, more lives are at stake. The time to act is now, not tomorrow.</p>
            </div>
            <p><i class="fas fa-laugh-squint"></i> <strong>"We didn't choose the name ESKAPE because it sounds cool (though it does). We chose it because it reminds us every single day: we must ESCAPE the AMR crisis — together, urgently, and with the best science we have."</strong><br>— Brown Beckley, lead developer (who secretly hopes this pun makes you smile, not roll your eyes 😉)</p>
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
            <i class="fas fa-quote-left"></i> "AMR is a silent pandemic, but we have the tools to fight it — if we share them, if we teach each other, and if we act with urgency. Let's escape the era of untreatable infections."<br>
            <strong>— The ESCAPE AMR Team, University of Ghana Medical School</strong>
        </div>'''

    # -------------------------------------------------------------------------
    # EXPORT
    # -------------------------------------------------------------------------
    def _generate_export_section(self, kwargs: Dict) -> str:
        return '''
        <div class="alert-box alert-info">
            <i class="fas fa-download fa-2x"></i>
            <div>
                <h3>📥 Export Data</h3>
                <p>Download tables as CSV using the buttons in each tab, or get the complete JSON data for AI-assisted analysis.</p>
            </div>
        </div>
        <div style="display:grid;grid-template-columns:repeat(auto-fit,minmax(280px,1fr));gap:20px;margin:30px 0;">
            <div class="dashboard-card card-export" onclick="exportTableToCSV('samples-table', 'sample_overview.csv')">
                <div style="font-size:2em;color:var(--export-color);"><i class="fas fa-table"></i></div>
                <div class="card-label">Sample Overview CSV</div></div>
            <div class="dashboard-card card-export" onclick="exportTableToCSV('qc-table', 'fasta_qc.csv')">
                <div style="font-size:2em;color:var(--export-color);"><i class="fas fa-chart-line"></i></div>
                <div class="card-label">FASTA QC CSV</div></div>
            <div class="dashboard-card card-export" onclick="location.href='staphscope_ultimate_sample_centric_report.json'">
                <div style="font-size:2em;color:var(--export-color);"><i class="fas fa-file-code"></i></div>
                <div class="card-label">Complete JSON Data</div></div>
        </div>
        <p style="color:#666;font-size:.92em;">
            <i class="fas fa-info-circle"></i>
            AMR, Virulence, BACMET, Plasmids, and Mutations tables are lazy-loaded. Click
            <strong>Show Details</strong> on each isolate box to reveal its tables; export them with
            the same <code>exportTableToCSV</code> function or copy from the JSON file.
        </p>'''


# -----------------------------------------------------------------------------
# ORCHESTRATOR
# -----------------------------------------------------------------------------
class StaphUltimateReporter:
    """Coordinates file discovery, parsing, integration, and report generation."""

    def __init__(self, input_dir: Path):
        self.input_dir = Path(input_dir)
        self.output_dir = self.input_dir / "STAPHSCOPE_ULTIMATE_SAMPLE_CENTRIC_REPORTS"
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.parser = StaphHTMLParser()
        self.analyzer = StaphDataAnalyzer()
        self.generator = StaphHTMLGenerator(self.analyzer)
        self.metadata = {
            "tool_name": "STAPHSCOPE Ultimate S. aureus Reporter",
            "version": "2.0.0",
            "author": "Brown Beckley <brownbeckley94@gmail.com>",
            "affiliation": "University of Ghana Medical School",
            "analysis_date": datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            "input_directory": str(self.input_dir),
        }

    def find_html_files(self) -> Dict[str, List[Path]]:
        """Discover QC HTML (fallback only; TSVs are the primary source)."""
        print("🔍 Searching for QC HTML fallback...")
        html_files = {'qc': []}
        for html_file in self.input_dir.glob("**/*.html"):
            fl = html_file.name.lower()
            if 'fasta_qc_summary' in fl or 'fasta_qc' in fl or 'qc_summary' in fl:
                html_files['qc'].append(html_file)
                print(f"    🎯 Found FASTA QC file: {html_file.name}")
        return html_files

    def integrate_all_data(self, html_files: Dict) -> Dict[str, Any]:
        """Load all sources (TSV primary), integrate per-sample data, build derived tables."""
        print("\n🔗 Integrating data (primary source: TSV files)…")
        integrated = {
            'metadata': self.metadata,
            'samples': {},
            'patterns': {},
            'gene_centric': {},
            'qc_data': {},
            'amrfinder_details': {},
            'abricate_details': {},
            'mutation_details': {},
        }

        # QC (HTML only)
        if html_files.get('qc'):
            integrated['qc_data'] = self.parser.load_qc_from_html(html_files['qc'][0])

        # Master typing TSV
        typing_data = self.parser.load_typing_from_tsv(self.input_dir)

        # Per-sample gene details
        amr_details, amr_freq = self.parser.load_amrfinder_from_tsv(self.input_dir)
        abricate_details, abricate_freq = self.parser.load_abricate_from_tsv(self.input_dir)

        # Per-sample mutations
        mutation_details = self.parser.load_mutations_from_tsv(self.input_dir)
        if mutation_details:
            integrated['mutation_details'] = mutation_details
            print(f"  ✅ Loaded mutations for {len(mutation_details)} samples")

        # Union of all sample IDs
        all_samples = (set(typing_data) | set(amr_details) | set(abricate_details)
                       | set(integrated['qc_data']) | set(mutation_details))
        all_samples = sorted(all_samples)
        if not all_samples:
            print("❌ No samples found in any source.")
            return {}

        print(f"📊 Found {len(all_samples)} unique samples")

        for sample in all_samples:
            typing = typing_data.get(sample, {
                'MLST': 'Not Assigned', 'spa_Type': 'Not Assigned',
                'agr_Type': 'Not Assigned', 'capsule_type': 'Not Assigned',
                'SCCmec_CGE': 'Not Assigned', 'SCCmec_RPet': 'Not Assigned',
                'SCCmec_Subtype': 'Not Assigned', 'MRSA_Status': 'Not Assigned',
            })

            amr_list = amr_details.get(sample, [])
            amr_gene_names = [d.get('gene', '') for d in amr_list if d.get('gene')]
            amr_info = {
                'critical_genes': [],
                'high_risk_genes': [],
                'all_genes': amr_gene_names,
            }
            for gene in amr_gene_names:
                gl = gene.lower()
                if gl in self.analyzer.critical_amr_genes:
                    amr_info['critical_genes'].append(gene)
                if gene in self.analyzer.high_priority_amr:
                    amr_info['high_risk_genes'].append(gene)

            abricate_info = {}
            for db, genes in abricate_details.get(sample, {}).items():
                abricate_info[db] = [g.get('gene', '') for g in genes if g.get('gene')]

            integrated['samples'][sample] = {
                'typing': typing,
                'amrfinder': amr_info,
                'abricate_databases': abricate_info,
            }

        integrated['amrfinder_details'] = amr_details
        integrated['abricate_details'] = abricate_details
        integrated['gene_frequencies'] = {
            'amrfinder': amr_freq,
            'abricate': abricate_freq,
        }

        print("\n🧠 Building gene-centric and pattern tables…")
        integrated['gene_centric'] = self.analyzer.create_gene_centric_tables(integrated)
        integrated['patterns'] = self.analyzer.create_cross_genome_patterns(integrated)
        return integrated

    def write_json(self, integrated_data: Dict[str, Any]) -> Path:
        """Write the full integrated dataset as JSON."""
        print("\n📝 Writing JSON report…")
        out = self.output_dir / "staphscope_ultimate_sample_centric_report.json"

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
            json.dump(serial(integrated_data), f, indent=2, ensure_ascii=False)
        print(f"    ✅ JSON saved: {out}")
        return out

    def write_csvs(self, integrated_data: Dict[str, Any]):
        """Write flat CSV exports next to the HTML report."""
        print("\n📊 Writing CSV reports…")

        # Sample overview
        rows = []
        for sample, data in integrated_data['samples'].items():
            t = data['typing']
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
                'Virulence_Gene_Count': len(data.get('abricate_databases', {}).get('vfdb', [])),
            })
        pd.DataFrame(rows).to_csv(self.output_dir / "sample_overview.csv", index=False)

        # Gene-centric CSVs
        gc = integrated_data.get('gene_centric', {})
        total = len(integrated_data['samples']) or 1
        for cat, fname in (
            ('amr_databases', 'amr_genes.csv'),
            ('virulence_databases', 'virulence_genes.csv'),
            ('bacmet_databases', 'bacmet_genes.csv'),
            ('plasmid_databases', 'plasmid_replicons.csv'),
        ):
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

        # Mutation CSV (per-sample flattened)
        muts = integrated_data.get('mutation_details', {})
        if muts:
            mut_rows = []
            for sample, mlist in muts.items():
                for m in mlist:
                    mut_rows.append({'Sample': sample, **m})
            if mut_rows:
                pd.DataFrame(mut_rows).to_csv(self.output_dir / "mutations.csv", index=False)

        # QC CSV
        if integrated_data.get('qc_data'):
            pd.DataFrame(
                [{'Sample': s, **m} for s, m in integrated_data['qc_data'].items()]
            ).to_csv(self.output_dir / "fasta_qc.csv", index=False)

        # Pattern discovery CSV
        P = integrated_data['patterns']
        pat_rows = []
        for k, v in P.get('mlst_distribution', {}).items():
            pat_rows.append({'Pattern_Type': 'MLST_Distribution',
                             'Combination': k, 'Count': v})
        for key in ('mlst_spa_combinations', 'mlst_sccmec_combinations',
                    'spa_sccmec_combinations', 'triple_combinations'):
            for combo, samples in P.get(key, {}).items():
                pat_rows.append({'Pattern_Type': key, 'Combination': combo,
                                 'Count': len(samples),
                                 'Samples': ';'.join(samples)})
        for c in P.get('high_risk_combinations', []):
            pat_rows.append({'Pattern_Type': 'High_Risk',
                             'Combination': c['sample'], 'Count': 1,
                             'Samples': c['sample']})
        if pat_rows:
            pd.DataFrame(pat_rows).to_csv(self.output_dir / "pattern_discovery.csv", index=False)

        print("    ✅ CSV reports written")

    def run(self) -> bool:
        """Execute the full pipeline end-to-end."""
        print("=" * 80)
        print("🧬 STAPHSCOPE ULTIMATE S. AUREUS REPORTER v2.0.0")
        print("   Hybrid: Gene-Centric for Typing + Lazy Sample-Centric for Genes")
        print("=" * 80)
        print(f"📁 Input:  {self.input_dir}")
        print(f"📁 Output: {self.output_dir}")

        html_files = self.find_html_files()
        integrated_data = self.integrate_all_data(html_files)
        if not integrated_data:
            return False

        print("\n" + "=" * 80)
        print("📊 GENERATING ULTIMATE STAPHSCOPE REPORTS")
        print("=" * 80)
        self.write_json(integrated_data)
        self.write_csvs(integrated_data)
        self.generator.generate_main_report(integrated_data, self.output_dir)

        n = len(integrated_data['samples'])
        mrsa = sum(1 for s in integrated_data['samples'].values()
                   if 'MRSA' in s['typing']['MRSA_Status'])
        n_agr = len(integrated_data['patterns'].get('agr_type_distribution', {}))
        n_mut = len(integrated_data.get('mutation_details', {}))

        print("\n" + "=" * 80)
        print("✅ REPORT COMPLETE")
        print("=" * 80)
        print(f"   Samples:              {n}")
        print(f"   MRSA:                 {mrsa}")
        print(f"   agr Types:            {n_agr}")
        print(f"   Samples with mutations: {n_mut}")
        print(f"   Output directory:     {self.output_dir}")
        print("=" * 80)
        return True


def main():
    """Command-line entry point."""
    parser = argparse.ArgumentParser(
        description='STAPHSCOPE Ultimate S. aureus Reporter v2.0.0',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""Examples:\n  python staphscope_ultimate_reporter.py -i /path/to/reports\n\nAuthor: Brown Beckley <brownbeckley94@gmail.com>""")
    parser.add_argument('-i', '--input-dir', required=True,
                        help='Directory containing StaphScope TSV summaries and QC HTML')
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
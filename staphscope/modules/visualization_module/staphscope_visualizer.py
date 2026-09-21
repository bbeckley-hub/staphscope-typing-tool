#!/usr/bin/env python3
"""
STAPHSCOPE VISUALIZER v2.0.0 - Unified Interactive + Static Report
===================================================================
Single-file visualization pipeline for S. aureus genomics.

Data source: the 8 files exported by the gene-centric module
  • staphscope_comprehensive_report.tsv   (typing master)
  • amr_genes.csv                          (AMR gene table)
  • virulence_genes.csv                    (VFDB gene table)
  • bacmet_genes.csv                       (BacMet gene table)
  • plasmid_replicons.csv                  (PlasmidFinder table)
  • mutations.csv                          (point-mutation table)
  • mge_profile.csv                        (mobileOG per-sample table)
  • fasta_qc.csv                           (FASTA QC + fastANI)

Produces:
  • staphscope_dashboard.html            — flagship interactive Plotly dashboard
  • PNG/ PDF/ SVG/                       — publication-quality static exports
  • DATA/                                — flat CSVs of every plot's underlying data
  • staphscope_visualization_report.txt  — plain-text summary
  • staphscope_visualizations_bundle.zip — one-click export bundle

Dashboard tabs:
  Overview · Typing · QC · AMR · Virulence · MGE · Resistance · Alerts · Story · Compare

Key features:
  • Auto-deduplicated gene tables (mecA == MECA == mec-A)
  • Auto-fitting axis labels (truncation + rotation + figure height)
  • Sunburst / Sankey / UpSet / Resistance matrix / Radial rings
  • QC box plots + ANI histogram with 95% species threshold line
  • Alert engine (severity-ranked, 5 rules)
  • Auto-narrative summary + Story Mode
  • Compare tool — side-by-side sample diff (deep, not shallow)
  • Cross-filtering — click a bar to reveal matching samples

Author: Brown Beckley <brownbeckley94@gmail.com>
Affiliation: University of Ghana Medical School
Version: 2.0.0
Date: 2026-09-20
MIT License
"""

import os
import re
import sys
import json
import argparse
import zipfile
from pathlib import Path
from datetime import datetime
from typing import Dict, List, Set, Tuple, Any, Optional
from collections import defaultdict, Counter
import warnings
warnings.filterwarnings('ignore')

import pandas as pd
import numpy as np
from bs4 import BeautifulSoup  # kept for future HTML fallbacks

# ------------------------------------------------------------------------------
# Matplotlib (static publication exports)
# ------------------------------------------------------------------------------
import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn as sns

plt.style.use('seaborn-v0_8-whitegrid')
mpl.rcParams['font.family'] = 'sans-serif'
mpl.rcParams['font.sans-serif'] = ['Arial', 'DejaVu Sans', 'Helvetica']
mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42
mpl.rcParams['savefig.dpi'] = 300
mpl.rcParams['savefig.bbox'] = 'tight'
mpl.rcParams['figure.max_open_warning'] = 50

# ------------------------------------------------------------------------------
# Plotly (interactive dashboard)
# ------------------------------------------------------------------------------
try:
    import plotly.graph_objects as go
    from plotly.subplots import make_subplots
    PLOTLY_AVAILABLE = True
except ImportError:
    PLOTLY_AVAILABLE = False


# ==============================================================================
# CONSTANTS
# ==============================================================================
TOOL_VERSION = "2.0.0"

PALETTE = {
    'mlst': '#FF9800', 'spa': '#9C27B0', 'sccmec': '#009688',
    'capsule': '#00ACC1', 'agr': '#8B5CF6', 'mrsa': '#DC143C',
    'mssa': '#4682B4', 'amr': '#F44336', 'virulence': '#E91E63',
    'bacmet': '#FF5722', 'plasmids': '#673AB7', 'mutation': '#00BCD4',
    'mge': '#16A085', 'primary': '#4CAF50', 'accent': '#7E22CE',
}

CHART_COLORS = [
    '#FF9800', '#9C27B0', '#009688', '#00ACC1', '#8B5CF6',
    '#F44336', '#E91E63', '#673AB7', '#16A085', '#4CAF50',
    '#0891B2', '#7C3AED', '#EC4899', '#14B8A6', '#F59E0B',
    '#8B5A2B', '#A0522D', '#CD5C5C', '#B8860B', '#2E8B57',
]

ALERT_RULES = {
    'critical_amr_virulence': {
        'severity': 'critical', 'icon': '☣️',
        'title': 'Critical AMR + Virulence Co-occurrence',
        'description': 'Carries both a last-resort resistance gene and a major virulence factor.',
    },
    'vancomycin_mrsa': {
        'severity': 'critical', 'icon': '💀',
        'title': 'Vancomycin-resistant MRSA',
        'description': 'MRSA isolate carries vanA or vanB — clinical red flag.',
    },
    'pvl_mrsa': {
        'severity': 'high', 'icon': '🦠',
        'title': 'PVL-positive MRSA',
        'description': 'PVL + MRSA — associated with severe skin and soft tissue infections.',
    },
    'mge_outlier': {
        'severity': 'medium', 'icon': '📱',
        'title': 'High MGE burden',
        'description': 'Sample is in the top 5% of mobileOG hit counts.',
    },
    'sccmec_disagreement': {
        'severity': 'medium', 'icon': '🔄',
        'title': 'SCCmec caller disagreement',
        'description': 'CGE and RPet callers disagree on the SCCmec type for this sample.',
    },
}

CRITICAL_AMR_GENES = {'meca', 'mecc', 'vana', 'vanb'}
CRITICAL_VIRULENCE_GENES = {'luks-pv', 'lukf-pv', 'tsst', 'tst'}


# ==============================================================================
# HELPERS
# ==============================================================================
_EXTENSIONS = ('.fna', '.fasta', '.fa', '.gb', '.gbk', '.gbff')


def normalize_sample(name: str) -> str:
    """Strip paths and FASTA extensions from a sample identifier."""
    s = str(name).strip()
    for ext in _EXTENSIONS:
        if s.endswith(ext):
            s = s[: -len(ext)]
    if '/' in s or '\\' in s:
        s = Path(s).name
    for ext in _EXTENSIONS:
        if s.endswith(ext):
            s = s[: -len(ext)]
    return s.strip()

def normalize_sccmec(value: str) -> str:
    """Normalize SCCmec type strings from different callers so they can be compared.

    CGE     : 'SCCmec_type_II(2A)'  → 'ii(2a)'
    RPet    : 'Type II(2A)'          → 'ii(2a)'
    Empty   : ''/'Not Assigned'/... → ''
    """
    s = str(value or '').strip()
    if not s or s in ('Not Assigned', 'ND', 'None', 'Not typed'):
        return ''
    # Strip known prefixes (case-insensitive)
    s = re.sub(r'^sccmec[_\s-]*type[_\s-]*', '', s, flags=re.I)
    s = re.sub(r'^type[_\s-]*', '', s, flags=re.I)
    # Collapse remaining whitespace / separators and lowercase for comparison
    return re.sub(r'[\s_]+', '', s).lower()


def esc(value: Any) -> str:
    """Escape a value for safe HTML embedding."""
    if value is None:
        return ""
    return (str(value).replace("&", "&amp;").replace("<", "&lt;")
            .replace(">", "&gt;").replace('"', "&quot;"))


def truncate_label(label: str, max_len: int = 35) -> str:
    """Truncate a label for axis rendering with an ellipsis."""
    s = str(label)
    if len(s) <= max_len:
        return s
    return s[: max_len - 1] + '…'


def auto_tickangle(n_categories: int) -> int:
    """Return a Plotly tickangle appropriate for the number of categories."""
    if n_categories <= 4:
        return 0
    if n_categories <= 8:
        return 25
    if n_categories <= 14:
        return 45
    if n_categories <= 22:
        return 60
    return 75


def auto_bar_height(n_categories: int, per_row: int = 22) -> int:
    """Return a Plotly figure height that scales with the number of categories."""
    return max(380, min(1400, 60 + per_row * n_categories))


# ==============================================================================
# GENE DEDUPLICATION
# ==============================================================================
# Keys and values are lowercase. Any gene name whose lowercased form matches a
# key will be merged into the corresponding canonical group.
GENE_ALIASES = {
    'meca': 'meca', 'mec a': 'meca', 'mec-a': 'meca',
    'mecc': 'mecc',
    'blaz': 'blaz', 'bla z': 'blaz', 'bla': 'blaz',
    'vana': 'vana', 'vanb': 'vanb', 'vanc': 'vanc',
    'erma': 'erma', 'ermb': 'ermb', 'ermc': 'ermc',
    'msra': 'msra', 'msrb': 'msrb',
    'mphc': 'mphc',
    'tetk': 'tetk', 'tetm': 'tetm', 'tetl': 'tetl',
    'tet38': 'tet38', 'tet(38)': 'tet38',
    'aaca-aphd': 'aaca-aphd', "aac(6')-aph(2'')": 'aaca-aphd',
    'aph3-iiia': 'aph3-iiia', "aph(3')-iiia": 'aph3-iiia',
    'ant4-ia': 'ant4-ia', "ant(4')-ia": 'ant4-ia',
    'ant6-ia': 'ant6-ia',
    'dfra': 'dfra', 'dfrg': 'dfrg',
    'cat': 'cat',
    'fosb': 'fosb',
    'pcob': 'pcob', 'pco': 'pcob',
    'luks-pv': 'luks-pv', 'lukf-pv': 'lukf-pv',
    'tsst': 'tsst', 'tst': 'tsst',
    'nora': 'nora',
    'mepa': 'mepa', 'lmrs': 'lmrs',
}


def canonical_gene_key(name: str) -> str:
    """Return a canonical key for gene deduplication (lowercase, alias-mapped)."""
    key = str(name).strip().lower()
    key = re.sub(r'\s+', ' ', key)
    return GENE_ALIASES.get(key, key)


def dedupe_gene_table(df: pd.DataFrame) -> pd.DataFrame:
    """Merge genes that differ only by case or alias.

    Combines the ``Genomes`` lists as a union; keeps the first-seen spelling
    for display. Joins databases if multiple.
    """
    if df.empty:
        return df
    df = df.copy()

    for col in ('Gene', 'Database', 'Count', 'Genomes'):
        if col not in df.columns:
            df[col] = '' if col != 'Count' else 0

    df['_key'] = df['Gene'].apply(canonical_gene_key)

    merged_rows = []
    for key, group in df.groupby('_key', sort=False):
        genomes: Set[str] = set()
        for g in group.get('Genomes', pd.Series(dtype=str)):
            for s in str(g).split(';'):
                s = s.strip()
                if s:
                    genomes.add(s)
        databases = sorted({str(d).strip() for d in group['Database'] if str(d).strip()})
        merged_rows.append({
            'Gene': group.iloc[0]['Gene'],
            'Database': ','.join(databases),
            'Count': len(genomes) if genomes else int(group['Count'].max()),
            'Genomes': ';'.join(sorted(genomes)),
        })

    out = pd.DataFrame(merged_rows)
    out = out.sort_values('Count', ascending=False).reset_index(drop=True)
    out['Percentage'] = (out['Count'] / out['Count'].max() * 100).round(1).astype(str) + '%'
    return out


# ==============================================================================
# DATA LOADER
# ==============================================================================
class StaphDataLoader:
    """Loads typing, gene tables, MGE, and QC data from the StaphScope module outputs."""

    GENE_CSVS = {
        'amr':       ('amr_genes.csv',         'AMR Genes'),
        'virulence': ('virulence_genes.csv',   'Virulence Genes'),
        'bacmet':    ('bacmet_genes.csv',      'BACMET Genes'),
        'plasmids':  ('plasmid_replicons.csv', 'Plasmid Replicons'),
    }
    MUTATION_CSV = 'mutations.csv'
    MGE_CSV = 'mge_profile.csv'
    QC_CSV = 'fasta_qc.csv'
    MASTER_TSV = 'staphscope_comprehensive_report.tsv'

    def __init__(self, input_dir: Path):
        self.input_dir = Path(input_dir)

    # ------------------------------------------------------------------
    # Typing (master TSV)
    # ------------------------------------------------------------------
    def load_typing(self) -> pd.DataFrame:
        path = self.input_dir / self.MASTER_TSV
        if not path.exists():
            print(f"  ⚠️ Master TSV not found: {path.name}")
            return pd.DataFrame()

        df = pd.read_csv(path, sep='\t', dtype=str).fillna('Not Assigned')
        rename = {
            'Sample': 'Sample',
            'MLST': 'MLST',
            'spa Type': 'spa_Type',
            'agr Type': 'agr_Type',
            'Capsule Type': 'Capsule_Type',
            'SCCmec Type (CGE)': 'SCCmec_CGE',
            'SCCmec Type (RPet)': 'SCCmec_RPet',
            'SCCmec Subtype': 'SCCmec_Subtype',
            'MRSA/MSSA Status': 'MRSA_Status',
        }
        df = df.rename(columns={k: v for k, v in rename.items() if k in df.columns})
        if 'Sample' not in df.columns:
            return pd.DataFrame()
        df['Sample'] = df['Sample'].apply(normalize_sample)
        if 'MLST' in df.columns:
            df['MLST'] = df['MLST'].astype(str).str.replace('ST', '', regex=False).str.strip()
        if 'MRSA_Status' in df.columns:
            df['MRSA_Status'] = df['MRSA_Status'].astype(str).str.upper()
        print(f"  ✅ Master TSV: {len(df)} samples")
        return df

    # ------------------------------------------------------------------
    # Gene tables (per-category CSV)
    # ------------------------------------------------------------------
    def load_gene_tables(self) -> Dict[str, pd.DataFrame]:
        out: Dict[str, pd.DataFrame] = {}
        for key, (fname, label) in self.GENE_CSVS.items():
            path = self.input_dir / fname
            if not path.exists():
                continue
            try:
                df = pd.read_csv(path, dtype=str).fillna('')
                df = dedupe_gene_table(df)
                out[key] = df
                print(f"  ✅ {label}: {len(df)} unique genes")
            except Exception as e:
                print(f"  ⚠️ {fname} load failed: {e}")
        return out

    def load_mutations(self) -> pd.DataFrame:
        path = self.input_dir / self.MUTATION_CSV
        if not path.exists():
            return pd.DataFrame()
        try:
            df = pd.read_csv(path, dtype=str).fillna('')
            df['Count'] = pd.to_numeric(df.get('Count', 0), errors='coerce').fillna(0).astype(int)
            df = df.sort_values('Count', ascending=False).reset_index(drop=True)
            print(f"  ✅ Mutations: {len(df)} entries")
            return df
        except Exception as e:
            print(f"  ⚠️ Mutations load failed: {e}")
            return pd.DataFrame()

    # ------------------------------------------------------------------
    # MGE (per-sample wide)
    # ------------------------------------------------------------------
    def load_mge(self) -> pd.DataFrame:
        path = self.input_dir / self.MGE_CSV
        if not path.exists():
            return pd.DataFrame()
        try:
            df = pd.read_csv(path, dtype=str).fillna('0')
            sample_col = next((c for c in df.columns if c.lower() == 'sample'), None)
            if not sample_col:
                return pd.DataFrame()
            df['Sample'] = df[sample_col].apply(normalize_sample)
            for col in df.columns:
                if col in ('Sample', sample_col):
                    continue
                df[col] = pd.to_numeric(
                    df[col].astype(str).str.replace(',', ''),
                    errors='coerce'
                ).fillna(0).astype(int)
            print(f"  ✅ MGE: {len(df)} samples")
            return df
        except Exception as e:
            print(f"  ⚠️ MGE load failed: {e}")
            return pd.DataFrame()

    # ------------------------------------------------------------------
    # FASTA QC + species
    # ------------------------------------------------------------------
    def load_qc(self) -> pd.DataFrame:
        path = self.input_dir / self.QC_CSV
        if not path.exists():
            return pd.DataFrame()
        try:
            df = pd.read_csv(path, dtype=str).fillna('')
            sample_col = next((c for c in df.columns if c.lower() == 'sample'), None)
            if sample_col:
                df['Sample'] = df[sample_col].apply(normalize_sample)
            numeric_keywords = ('N50', 'N75', 'N90', 'Length', 'Sequences',
                                'Content', 'bases', 'run', 'Homopolymer',
                                'Duplicate', 'Size', 'ANI')
            for col in df.columns:
                if col in ('Sample', sample_col):
                    continue
                if any(k.lower() in col.lower() for k in numeric_keywords):
                    df[col] = pd.to_numeric(
                        df[col].astype(str)
                              .str.replace('%', '').str.replace(',', '').str.strip(),
                        errors='coerce'
                    )
            print(f"  ✅ QC: {len(df)} samples")
            return df
        except Exception as e:
            print(f"  ⚠️ QC load failed: {e}")
            return pd.DataFrame()


# ==============================================================================
# GENE MATRIX
# ==============================================================================
def build_gene_matrix(gene_tables: Dict[str, pd.DataFrame]) -> pd.DataFrame:
    """Build a samples × genes 0/1 matrix from the CSVs' ``Genomes`` column."""
    rows: List[Tuple[str, str]] = []
    for cat, df in (gene_tables or {}).items():
        if df is None or df.empty:
            continue
        for _, r in df.iterrows():
            gene = str(r.get('Gene', '')).strip()
            if not gene:
                continue
            for s in str(r.get('Genomes', '')).split(';'):
                s = s.strip()
                if s:
                    rows.append((s, gene))
    if not rows:
        return pd.DataFrame()
    long = pd.DataFrame(rows, columns=['Sample', 'Gene']).drop_duplicates()
    return long.assign(Value=1).pivot_table(
        index='Sample', columns='Gene', values='Value', fill_value=0
    )


# ==============================================================================
# ALERT ENGINE
# ==============================================================================
class AlertEngine:
    """Evaluates alert rules against the dataset."""

    def __init__(self):
        self.alerts: List[Dict[str, Any]] = []

    def run(
        self,
        typing_df: pd.DataFrame,
        gene_tables: Dict[str, pd.DataFrame],
        mge_df: pd.DataFrame,
    ) -> List[Dict[str, Any]]:
        # Build per-sample gene sets
        sample_genes: Dict[str, Set[str]] = defaultdict(set)
        for df in (gene_tables or {}).values():
            if df is None or df.empty:
                continue
            for _, r in df.iterrows():
                gene = str(r.get('Gene', '')).strip().lower()
                for s in str(r.get('Genomes', '')).split(';'):
                    s = s.strip()
                    if s and gene:
                        sample_genes[s].add(gene)

        # Rule 1 — critical AMR + virulence
        for sample, genes in sample_genes.items():
            amr_hits = genes & CRITICAL_AMR_GENES
            vir_hits = genes & CRITICAL_VIRULENCE_GENES
            if amr_hits and vir_hits:
                self.alerts.append({
                    'id': 'critical_amr_virulence',
                    'sample': sample,
                    'detail': f"AMR: {', '.join(sorted(amr_hits))} · "
                              f"Virulence: {', '.join(sorted(vir_hits))}",
                })

        # Rule 2 — vancomycin-resistant MRSA
        if not typing_df.empty and 'MRSA_Status' in typing_df.columns:
            for _, row in typing_df.iterrows():
                if row.get('MRSA_Status') == 'MRSA':
                    s = row['Sample']
                    van = sample_genes[s] & {'vana', 'vanb'}
                    if van:
                        self.alerts.append({
                            'id': 'vancomycin_mrsa', 'sample': s,
                            'detail': f"van gene detected: {', '.join(sorted(van))}",
                        })

        # Rule 3 — PVL + MRSA
        if not typing_df.empty and 'MRSA_Status' in typing_df.columns:
            for _, row in typing_df.iterrows():
                if row.get('MRSA_Status') == 'MRSA':
                    s = row['Sample']
                    if sample_genes[s] & {'luks-pv', 'lukf-pv'}:
                        self.alerts.append({
                            'id': 'pvl_mrsa', 'sample': s,
                            'detail': 'MRSA + PVL detected',
                        })

        # Rule 4 — MGE outlier
        if mge_df is not None and not mge_df.empty and 'mobileOG Hits' in mge_df.columns:
            totals = dict(zip(mge_df['Sample'], mge_df['mobileOG Hits']))
            values = sorted(totals.values())
            if values:
                cutoff = values[int(len(values) * 0.95)]
                for s, t in totals.items():
                    if t >= cutoff and t > 0:
                        self.alerts.append({
                            'id': 'mge_outlier', 'sample': s,
                            'detail': f"{t} mobileOG hits (top 5%)",
                        })

        # Rule 5 — SCCmec caller disagreement (normalize caller-specific formatting)
        if (not typing_df.empty and 'SCCmec_CGE' in typing_df.columns
                and 'SCCmec_RPet' in typing_df.columns):
            for _, row in typing_df.iterrows():
                cge_raw = str(row.get('SCCmec_CGE', '')).strip()
                rpet_raw = str(row.get('SCCmec_RPet', '')).strip()
                cge_norm = normalize_sccmec(cge_raw)
                rpet_norm = normalize_sccmec(rpet_raw)
                # Skip if either caller didn't assign a type
                if not cge_norm or not rpet_norm:
                    continue
                # Only flag genuine disagreements (subtype-level differences included,
                # e.g. IIa vs IIb) — not formatting differences
                if cge_norm != rpet_norm:
                    self.alerts.append({
                        'id': 'sccmec_disagreement',
                        'sample': row['Sample'],
                        'detail': f"CGE: {cge_raw} · RPet: {rpet_raw}",
                    })

        order = {'critical': 0, 'high': 1, 'medium': 2}
        self.alerts.sort(key=lambda a: order.get(ALERT_RULES[a['id']]['severity'], 3))
        return self.alerts


# ==============================================================================
# INTERACTIVE DASHBOARD BUILDER
# ==============================================================================
class DashboardBuilder:
    """Builds the single-file Plotly dashboard."""

    def __init__(self, output_dir: Path):
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.charts: Dict[str, str] = {}
        self.alerts: List[Dict[str, Any]] = []
        self.narrative: List[str] = []
        self.story_chapters: List[Dict[str, Any]] = []
        self.meta: Dict[str, Any] = {}
        self.sample_table_html: str = ''

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------
    def build(
        self,
        typing_df: pd.DataFrame,
        gene_tables: Dict[str, pd.DataFrame],
        mge_df: pd.DataFrame,
        qc_df: pd.DataFrame,
    ) -> Optional[Path]:
        if typing_df.empty:
            print("  ⚠️ Dashboard skipped — no typing data")
            return None

        print("🎨 Building interactive dashboard...")

        gene_matrix = build_gene_matrix(gene_tables)

        # View-dependent charts
        for view in ('all', 'mrsa', 'mssa'):
            view_df = self._slice_view(typing_df, view)
            if view_df.empty and view != 'all':
                continue
            self._render_mrsa_donut(view, view_df)
            self._render_sunburst(view, view_df)
            self._render_sankey(view, view_df)
            self._render_bar(view, view_df, 'MLST', 'MLST Sequence Types')
            self._render_bar(view, view_df, 'spa_Type', 'spa Types')
            self._render_bar(view, view_df, 'SCCmec_CGE', 'SCCmec Types (CGE)')
            self._render_bar(view, view_df, 'SCCmec_Subtype', 'SCCmec Subtypes')
            self._render_bar(view, view_df, 'agr_Type', 'agr Types')
            self._render_bar(view, view_df, 'Capsule_Type', 'Capsule Types')

        # View-independent charts
        self._render_amr_cooccurrence(gene_matrix)
        self._render_amr_db_comparison(gene_tables)
        self._render_virulence_top(gene_tables)
        self._render_mge_breakdown(mge_df)
        self._render_qc_boxplots(qc_df)
        self._render_ani_chart(qc_df)
        self._render_upset(gene_matrix)
        self._render_resistance_matrix(gene_matrix)
        self._render_radial_rings(typing_df)

        self.sample_table_html = self._build_sample_table(typing_df)

        engine = AlertEngine()
        self.alerts = engine.run(typing_df, gene_tables, mge_df)
        self._narrate(typing_df, gene_tables, mge_df, qc_df)
        self._build_story(typing_df, gene_tables, gene_matrix, qc_df)

        self.meta = {
            'generated': datetime.now().strftime('%Y-%m-%d %H:%M:%S'),
            'n_samples': len(typing_df),
            'n_mrsa': int((typing_df.get('MRSA_Status', pd.Series()) == 'MRSA').sum()),
            'n_mssa': int((typing_df.get('MRSA_Status', pd.Series()) == 'MSSA').sum()),
            'n_alerts': len(self.alerts),
            'tool_version': TOOL_VERSION,
        }

        html = self._assemble_html(typing_df, qc_df)
        out = self.output_dir / 'staphscope_dashboard.html'
        out.write_text(html, encoding='utf-8')
        size_mb = out.stat().st_size / (1024 * 1024)
        print(f"  ✅ Dashboard: {out.name} ({size_mb:.1f} MB)")
        return out

    # ------------------------------------------------------------------
    # Helpers
    # ------------------------------------------------------------------
    @staticmethod
    def _slice_view(df: pd.DataFrame, view: str) -> pd.DataFrame:
        if view == 'all' or 'MRSA_Status' not in df.columns:
            return df
        return df[df['MRSA_Status'] == view.upper()]

    # ------------------------------------------------------------------
    # Chart renderers
    # ------------------------------------------------------------------
    def _store(self, key: str, fig, view: str = 'all'):
        html = fig.to_html(
            include_plotlyjs=False, full_html=False,
            div_id=f"chart-{key}-{view}",
            config={
                'displayModeBar': True, 'displaylogo': False,
                'modeBarButtonsToRemove': ['lasso2d', 'select2d', 'autoScale2d'],
                'toImageButtonOptions': {
                    'format': 'png',
                    'filename': f"staphscope_{key}_{view}",
                    'height': 900, 'width': 1400, 'scale': 2,
                },
            },
        )
        self.charts[f"{key}_{view}"] = html

    def _chart_block(self, key: str, title: str, subtitle: str = '') -> str:
        variants = []
        for v in ('all', 'mrsa', 'mssa'):
            chart_key = f"{key}_{v}"
            if chart_key not in self.charts:
                continue
            display = 'block' if v == 'all' else 'none'
            variants.append(
                f'<div class="chart-view" data-view="{v}" style="display:{display};">'
                f'{self.charts[chart_key]}</div>'
            )
        if not variants:
            return ''
        sub_html = f'<p class="chart-subtitle">{subtitle}</p>' if subtitle else ''
        return f'''<div class="chart-card">
            <h3 class="chart-title">{title}</h3>
            {sub_html}
            {''.join(variants)}
        </div>'''

    # ---------- MRSA donut ----------
    def _render_mrsa_donut(self, view: str, df: pd.DataFrame):
        if 'MRSA_Status' not in df.columns:
            return
        counts = df['MRSA_Status'].value_counts()
        counts = counts[counts.index.isin(['MRSA', 'MSSA'])]
        if counts.empty:
            return
        fig = go.Figure(go.Pie(
            labels=counts.index.tolist(),
            values=counts.values.tolist(),
            hole=0.55,
            marker=dict(
                colors=[PALETTE['mrsa'] if l == 'MRSA' else PALETTE['mssa']
                        for l in counts.index],
                line=dict(color='#0f172a', width=2),
            ),
            textinfo='label+percent+value',
            hovertemplate='<b>%{label}</b><br>%{value} samples (%{percent})<extra></extra>',
        ))
        fig.update_layout(
            template='plotly_dark', paper_bgcolor='rgba(0,0,0,0)',
            plot_bgcolor='rgba(0,0,0,0)', height=380, showlegend=False,
            margin=dict(t=10, l=10, r=10, b=10),
            annotations=[dict(
                text=f"<b>{counts.sum()}</b><br><span style='font-size:12px'>samples</span>",
                x=0.5, y=0.5, font=dict(size=22, color='#f1f5f9'), showarrow=False,
            )],
            font=dict(family='Segoe UI, sans-serif'),
        )
        self._store('mrsa_donut', fig, view)

    # ---------- Sunburst ----------
    def _render_sunburst(self, view: str, df: pd.DataFrame):
        if 'MRSA_Status' not in df.columns or 'SCCmec_CGE' not in df.columns:
            return
        df = df[df['MRSA_Status'].isin(['MRSA', 'MSSA'])].copy()
        if df.empty:
            return
        has_mlst = 'MLST' in df.columns

        labels: List[str] = []
        parents: List[str] = []
        ids: List[str] = []
        colors: List[str] = []

        def add(node_id: str, label: str, parent: str, color: str):
            if node_id in ids:
                return
            ids.append(node_id)
            labels.append(label)
            parents.append(parent)
            colors.append(color)

        for status in ('MRSA', 'MSSA'):
            sub = df[df['MRSA_Status'] == status]
            if sub.empty:
                continue
            add(f"status:{status}", status, '',
                PALETTE['mrsa' if status == 'MRSA' else 'mssa'])
            for scc in sub['SCCmec_CGE'].unique():
                if scc in ('Not Assigned', '', 'ND'):
                    continue
                scc_id = f"status:{status}|scc:{scc}"
                add(scc_id, scc, f"status:{status}", PALETTE['sccmec'])
                if not has_mlst:
                    continue
                sub2 = sub[sub['SCCmec_CGE'] == scc]
                for st in sub2['MLST'].unique():
                    if st in ('Not Assigned', '', 'ND'):
                        continue
                    st_id = f"{scc_id}|st:{st}"
                    add(st_id, f"ST{st}", scc_id, PALETTE['mlst'])

        values = []
        for node_id in ids:
            parts = node_id.split('|')
            mask = pd.Series(True, index=df.index)
            for part in parts:
                k, v = part.split(':', 1)
                if k == 'status':
                    mask &= (df['MRSA_Status'] == v)
                elif k == 'scc':
                    mask &= (df['SCCmec_CGE'] == v)
                elif k == 'st':
                    mask &= (df['MLST'] == v)
            values.append(int(mask.sum()))

        if not values:
            return
        fig = go.Figure(go.Sunburst(
            ids=ids, labels=labels, parents=parents, values=values,
            branchvalues='total',
            marker=dict(colors=colors, line=dict(color='#0f172a', width=2)),
            hovertemplate='<b>%{label}</b><br>Samples: %{value}<extra></extra>',
            textinfo='label+value', insidetextorientation='radial',
        ))
        fig.update_layout(
            template='plotly_dark', paper_bgcolor='rgba(0,0,0,0)',
            plot_bgcolor='rgba(0,0,0,0)', height=560,
            margin=dict(t=10, l=10, r=10, b=10),
            font=dict(family='Segoe UI, sans-serif', size=12),
        )
        self._store('sunburst', fig, view)

    # ---------- Sankey ----------
    def _render_sankey(self, view: str, df: pd.DataFrame):
        cols = ['MLST', 'SCCmec_CGE', 'agr_Type', 'Capsule_Type']
        if not all(c in df.columns for c in cols):
            return
        df = df.copy()
        for c in cols:
            df = df[~df[c].isin(['Not Assigned', '', 'ND', 'Unknown'])]
        if len(df) < 2:
            return

        nodes: List[str] = []
        node_idx: Dict[str, int] = {}
        src, tgt, val = [], [], []
        layer_palette = [PALETTE['mlst'], PALETTE['sccmec'],
                         PALETTE['agr'], PALETTE['capsule']]

        def nid(label: str, layer: int) -> int:
            key = f"{layer}:{label}"
            if key not in node_idx:
                node_idx[key] = len(nodes)
                nodes.append(label)
            return node_idx[key]

        for i in range(len(cols) - 1):
            a, b = cols[i], cols[i + 1]
            pairs = df.groupby([a, b]).size().reset_index(name='count')
            for _, row in pairs.iterrows():
                src.append(nid(str(row[a]), i))
                tgt.append(nid(str(row[b]), i + 1))
                val.append(int(row['count']))

        node_colors = []
        for i in range(len(nodes)):
            for key, idx in node_idx.items():
                if idx == i:
                    layer = int(key.split(':', 1)[0])
                    node_colors.append(layer_palette[layer % len(layer_palette)])
                    break

        fig = go.Figure(go.Sankey(
            node=dict(
                pad=20, thickness=18,
                line=dict(color='#0f172a', width=1),
                label=nodes, color=node_colors,
                hovertemplate='<b>%{label}</b><br>Flow: %{value}<extra></extra>',
            ),
            link=dict(
                source=src, target=tgt, value=val,
                color='rgba(148,163,184,0.25)',
                hovertemplate='%{source.label} → %{target.label}<br>%{value} samples<extra></extra>',
            ),
        ))
        fig.update_layout(
            template='plotly_dark', paper_bgcolor='rgba(0,0,0,0)',
            plot_bgcolor='rgba(0,0,0,0)', height=560,
            margin=dict(t=10, l=10, r=10, b=10),
            font=dict(family='Segoe UI, sans-serif', size=12),
        )
        self._store('sankey', fig, view)

    # ---------- Typing bar (auto labels) ----------
    def _render_bar(self, view: str, df: pd.DataFrame, column: str, title: str,
                    top_n: int = 20):
        if column not in df.columns:
            return
        counts = df[column].value_counts()
        counts = counts[~counts.index.isin(['Not Assigned', '', 'ND', 'Unknown'])]
        if counts.empty:
            return
        top = counts.head(top_n)
        n = len(top)

        max_label_len = 45 if n <= 8 else 35 if n <= 15 else 28
        labels = [truncate_label(str(i), max_label_len) for i in top.index]
        sample_lists = [
            ';'.join(df[df[column] == label]['Sample'].tolist())
            for label in top.index
        ]

        fig = go.Figure(go.Bar(
            x=top.values,
            y=labels,
            orientation='h',
            marker=dict(color=CHART_COLORS[:n], line=dict(color='#0f172a', width=1)),
            text=top.values, textposition='outside',
            customdata=sample_lists,
            hovertemplate='<b>%{y}</b><br>%{x} samples<extra></extra>',
        ))
        fig.update_layout(
            template='plotly_dark', paper_bgcolor='rgba(0,0,0,0)',
            plot_bgcolor='rgba(0,0,0,0)',
            height=auto_bar_height(n),
            showlegend=False,
            margin=dict(t=10, l=180, r=40, b=10),
            yaxis=dict(autorange='reversed', tickfont=dict(size=11), automargin=True),
            xaxis=dict(gridcolor='rgba(148,163,184,0.15)', title=None),
            font=dict(family='Segoe UI, sans-serif', size=11),
        )
        self._store(f'bar-{column}', fig, view)

    # ---------- AMR co-occurrence heatmap ----------
    def _render_amr_cooccurrence(self, gene_matrix: pd.DataFrame):
        if gene_matrix.empty:
            return
        top_genes = gene_matrix.sum(axis=0).sort_values(ascending=False).head(30).index.tolist()
        mat = gene_matrix[top_genes]
        mat = mat.loc[mat.sum(axis=1).sort_values(ascending=False).index]
        if mat.empty:
            return
        x_labels = [truncate_label(g, 24) for g in mat.columns]
        fig = go.Figure(go.Heatmap(
            z=mat.values,
            x=x_labels,
            y=mat.index.tolist(),
            colorscale=[[0.0, '#0f172a'], [1.0, PALETTE['amr']]],
            showscale=False, xgap=1, ygap=1,
            hovertemplate='<b>%{y}</b><br>Gene: %{x}<br>Present: %{z}<extra></extra>',
        ))
        fig.update_layout(
            template='plotly_dark', paper_bgcolor='rgba(0,0,0,0)',
            plot_bgcolor='rgba(0,0,0,0)',
            height=max(420, 20 * len(mat)),
            margin=dict(t=10, l=10, r=10, b=10),
            xaxis=dict(tickangle=auto_tickangle(len(x_labels)),
                       tickfont=dict(size=10), automargin=True),
            yaxis=dict(tickfont=dict(size=10), automargin=True),
            font=dict(family='Segoe UI, sans-serif'),
        )
        self._store('amr_cooccurrence', fig, 'all')

    # ---------- Category frequency comparison ----------
    def _render_amr_db_comparison(self, gene_tables: Dict[str, pd.DataFrame]):
        rows = []
        for cat, df in (gene_tables or {}).items():
            if df is None or df.empty:
                continue
            for _, r in df.iterrows():
                try:
                    rows.append({
                        'Category': cat.capitalize(),
                        'Count': int(r.get('Count', 0)),
                    })
                except Exception:
                    continue
        if not rows:
            return
        df = pd.DataFrame(rows)

        fig = go.Figure()
        for i, cat in enumerate(df['Category'].unique()):
            sub = df[df['Category'] == cat]['Count']
            fig.add_trace(go.Box(
                y=sub, name=cat,
                marker_color=CHART_COLORS[i % len(CHART_COLORS)],
                boxmean='sd',
                hovertemplate='<b>%{x}</b><br>Count: %{y}<extra></extra>',
            ))
        fig.update_layout(
            template='plotly_dark', paper_bgcolor='rgba(0,0,0,0)',
            plot_bgcolor='rgba(0,0,0,0)', height=440, showlegend=False,
            margin=dict(t=10, l=10, r=10, b=10),
            yaxis=dict(title='Number of genomes carrying gene',
                       gridcolor='rgba(148,163,184,0.15)'),
            xaxis=dict(tickangle=auto_tickangle(df['Category'].nunique())),
            font=dict(family='Segoe UI, sans-serif'),
        )
        self._store('amr_db_comparison', fig, 'all')

    # ---------- Virulence top genes ----------
    def _render_virulence_top(self, gene_tables: Dict[str, pd.DataFrame]):
        df = (gene_tables or {}).get('virulence')
        if df is None or df.empty:
            return
        top = df.nlargest(20, 'Count')
        n = len(top)
        labels = [truncate_label(g, 30) for g in top['Gene'].astype(str)]
        fig = go.Figure(go.Bar(
            x=top['Count'], y=labels,
            orientation='h',
            marker=dict(color=PALETTE['virulence']),
            text=top['Count'], textposition='outside',
            hovertemplate='<b>%{y}</b><br>%{x} genomes<extra></extra>',
        ))
        fig.update_layout(
            template='plotly_dark', paper_bgcolor='rgba(0,0,0,0)',
            plot_bgcolor='rgba(0,0,0,0)',
            height=auto_bar_height(n, 24),
            showlegend=False,
            margin=dict(t=10, l=180, r=40, b=10),
            yaxis=dict(autorange='reversed', tickfont=dict(size=10), automargin=True),
            xaxis=dict(gridcolor='rgba(148,163,184,0.15)', title=None),
            font=dict(family='Segoe UI, sans-serif'),
        )
        self._store('virulence_top', fig, 'all')

    # ---------- MGE breakdown ----------
    def _render_mge_breakdown(self, mge_df: pd.DataFrame):
        if mge_df is None or mge_df.empty:
            return
        skip = {'Sample', 'mobileOG Hits'}
        cats = [c for c in mge_df.columns if c not in skip]
        if not cats:
            return
        means = {c: float(mge_df[c].mean()) for c in cats}
        cats_sorted = sorted(means, key=lambda c: -means[c])
        n = len(cats_sorted)

        fig = go.Figure(go.Bar(
            x=cats_sorted, y=[means[c] for c in cats_sorted],
            marker=dict(color=CHART_COLORS[:n]),
            text=[f"{means[c]:.1f}" for c in cats_sorted], textposition='outside',
            hovertemplate='<b>%{x}</b><br>Mean hits/sample: %{y:.1f}<extra></extra>',
        ))
        fig.update_layout(
            template='plotly_dark', paper_bgcolor='rgba(0,0,0,0)',
            plot_bgcolor='rgba(0,0,0,0)', height=440, showlegend=False,
            margin=dict(t=10, l=10, r=10, b=80),
            xaxis=dict(tickangle=auto_tickangle(n), automargin=True),
            yaxis=dict(title='Mean mobileOG hits',
                       gridcolor='rgba(148,163,184,0.15)'),
            font=dict(family='Segoe UI, sans-serif'),
        )
        self._store('mge_breakdown', fig, 'all')

    # ---------- QC box plots ----------
    def _render_qc_boxplots(self, qc_df: pd.DataFrame):
        if qc_df is None or qc_df.empty:
            return

        metrics = []
        if 'N50' in qc_df.columns:
            metrics.append(('N50', 'N50 (bp)', '#FF9800'))
        if 'Total Length' in qc_df.columns:
            metrics.append(('Total Length', 'Assembly size (bp)', '#4CAF50'))
        if 'GC Content (%)' in qc_df.columns:
            metrics.append(('GC Content (%)', 'GC content (%)', '#009688'))
        if 'Total Sequences' in qc_df.columns:
            metrics.append(('Total Sequences', 'Number of contigs', '#7E22CE'))
        if 'Ambiguous Bases (%)' in qc_df.columns:
            metrics.append(('Ambiguous Bases (%)', 'Ambiguous bases (%)', '#F44336'))

        if not metrics:
            return

        n = len(metrics)
        cols = min(2, n)
        rows = (n + cols - 1) // cols
        fig = make_subplots(
            rows=rows, cols=cols,
            subplot_titles=[m[1] for m in metrics],
            vertical_spacing=0.15,
        )
        for i, (col, label, color) in enumerate(metrics):
            r = i // cols + 1
            c = i % cols + 1
            series = qc_df[col].dropna()
            if series.empty:
                continue
            fig.add_trace(
                go.Box(
                    y=series, name=label, marker_color=color,
                    boxmean='sd', boxpoints='all',
                    jitter=0.4, pointpos=0,
                    hovertemplate=f'<b>{label}</b><br>%{{y}}<extra></extra>',
                ),
                row=r, col=c,
            )
        fig.update_layout(
            template='plotly_dark', paper_bgcolor='rgba(0,0,0,0)',
            plot_bgcolor='rgba(0,0,0,0)',
            height=max(400, 300 * rows),
            showlegend=False,
            margin=dict(t=40, l=60, r=20, b=30),
            font=dict(family='Segoe UI, sans-serif'),
        )
        self._store('qc_boxplots', fig, 'all')

    # ---------- ANI chart ----------
    def _render_ani_chart(self, qc_df: pd.DataFrame):
        if qc_df is None or qc_df.empty or 'ANI (%)' not in qc_df.columns:
            return
        df = qc_df.dropna(subset=['ANI (%)'])
        if df.empty:
            return

        fig = go.Figure()
        fig.add_trace(go.Histogram(
            x=df['ANI (%)'],
            nbinsx=20,
            marker=dict(color=PALETTE['primary'],
                        line=dict(color='#0f172a', width=1)),
            hovertemplate='ANI %{x}<br>%{y} samples<extra></extra>',
            name='ANI',
        ))
        fig.add_vline(
            x=95, line_dash='dash', line_color=PALETTE['mrsa'],
            annotation_text='95% threshold (species boundary)',
            annotation_position='top right',
        )
        fig.update_layout(
            template='plotly_dark', paper_bgcolor='rgba(0,0,0,0)',
            plot_bgcolor='rgba(0,0,0,0)',
            height=400, showlegend=False,
            margin=dict(t=40, l=60, r=20, b=40),
            xaxis=dict(title='ANI (%)', gridcolor='rgba(148,163,184,0.15)'),
            yaxis=dict(title='Samples', gridcolor='rgba(148,163,184,0.15)'),
            font=dict(family='Segoe UI, sans-serif'),
        )
        self._store('ani_chart', fig, 'all')

    # ---------- UpSet-style combination matrix ----------
    def _render_upset(self, gene_matrix: pd.DataFrame):
        if gene_matrix.empty or gene_matrix.shape[1] < 2:
            return
        top_genes = gene_matrix.sum(axis=0).sort_values(ascending=False).head(10).index.tolist()
        mat = gene_matrix[top_genes]
        combos = mat.apply(lambda r: tuple(int(v) for v in r.values), axis=1)
        counts = Counter(combos)
        ranked = sorted(counts.items(), key=lambda kv: -kv[1])[:15]
        if not ranked:
            return
        labels = [' · '.join(top_genes[i] for i, b in enumerate(combo) if b)
                  or '(no genes)' for combo, _ in ranked]
        labels = [truncate_label(l, 60) for l in labels]
        values = [c for _, c in ranked]
        n = len(ranked)

        fig = go.Figure(go.Bar(
            x=values, y=labels, orientation='h',
            marker=dict(color=CHART_COLORS[:n]),
            text=values, textposition='outside',
            hovertemplate='<b>%{y}</b><br>%{x} samples<extra></extra>',
        ))
        fig.update_layout(
            template='plotly_dark', paper_bgcolor='rgba(0,0,0,0)',
            plot_bgcolor='rgba(0,0,0,0)',
            height=auto_bar_height(n, 26),
            showlegend=False,
            margin=dict(t=10, l=200, r=40, b=10),
            xaxis=dict(title='Samples with this combination',
                       gridcolor='rgba(148,163,184,0.15)'),
            yaxis=dict(autorange='reversed', tickfont=dict(size=10), automargin=True),
            font=dict(family='Segoe UI, sans-serif'),
        )
        self._store('upset', fig, 'all')

    # ---------- Resistance presence/absence matrix ----------
    def _render_resistance_matrix(self, gene_matrix: pd.DataFrame):
        if gene_matrix.empty:
            return
        top_genes = gene_matrix.sum(axis=0).sort_values(ascending=False).head(25).index.tolist()
        mat = gene_matrix[top_genes]
        if mat.empty:
            return
        try:
            from scipy.cluster.hierarchy import linkage, leaves_list
            from scipy.spatial.distance import pdist
            if len(mat) > 2:
                dists = pdist(mat.values, metric='jaccard')
                if np.any(dists):
                    Z = linkage(dists, method='average')
                    order = leaves_list(Z)
                    mat = mat.iloc[order]
        except Exception:
            pass

        x_labels = [truncate_label(g, 22) for g in mat.columns]
        fig = go.Figure(go.Heatmap(
            z=mat.values,
            x=x_labels,
            y=mat.index.tolist(),
            colorscale=[[0.0, '#1e293b'], [1.0, PALETTE['amr']]],
            showscale=False, xgap=1, ygap=1,
            hovertemplate='<b>%{y}</b><br>%{x}: %{z}<extra></extra>',
        ))
        fig.update_layout(
            template='plotly_dark', paper_bgcolor='rgba(0,0,0,0)',
            plot_bgcolor='rgba(0,0,0,0)',
            height=max(420, 20 * len(mat)),
            margin=dict(t=10, l=10, r=10, b=10),
            xaxis=dict(tickangle=auto_tickangle(len(x_labels)),
                       tickfont=dict(size=10), automargin=True),
            yaxis=dict(tickfont=dict(size=9), automargin=True),
            font=dict(family='Segoe UI, sans-serif'),
        )
        self._store('resistance_matrix', fig, 'all')

    # ---------- Radial typing rings ----------
    def _render_radial_rings(self, typing_df: pd.DataFrame):
        layers = [c for c in ('MLST', 'spa_Type', 'SCCmec_CGE', 'agr_Type',
                              'Capsule_Type', 'MRSA_Status')
                  if c in typing_df.columns]
        if not layers:
            return
        samples = typing_df['Sample'].tolist()
        theta = np.linspace(0, 360, len(samples), endpoint=False).tolist()

        fig = go.Figure()
        for layer in layers:
            values = typing_df[layer].astype(str).tolist()
            uniq = sorted({v for v in values if v not in ('Not Assigned', '', 'ND', 'Unknown')})
            cmap = {v: CHART_COLORS[i % len(CHART_COLORS)] for i, v in enumerate(uniq)}
            colors = [cmap.get(v, '#334155') for v in values]
            fig.add_trace(go.Barpolar(
                r=[1] * len(samples),
                theta=theta,
                width=[360 / len(samples) * 0.95] * len(samples),
                base=[layers.index(layer)] * len(samples),
                marker=dict(color=colors, line=dict(color='#0f172a', width=0.5)),
                name=layer,
                hovertemplate=('<b>%{customdata[0]}</b><br>'
                               + layer + ': %{customdata[1]}<extra></extra>'),
                customdata=list(zip(samples, values)),
            ))
        fig.update_layout(
            template='plotly_dark', paper_bgcolor='rgba(0,0,0,0)',
            plot_bgcolor='rgba(0,0,0,0)', height=560,
            margin=dict(t=20, l=20, r=20, b=20),
            polar=dict(
                radialaxis=dict(showticklabels=False, ticks='',
                                gridcolor='rgba(148,163,184,0.15)'),
                angularaxis=dict(showticklabels=False, ticks='',
                                 direction='clockwise'),
                bgcolor='rgba(0,0,0,0)',
            ),
            font=dict(family='Segoe UI, sans-serif'),
            showlegend=True,
            legend=dict(orientation='h', yanchor='bottom', y=-0.05,
                        x=0.5, xanchor='center'),
        )
        self._store('radial_rings', fig, 'all')

    # ------------------------------------------------------------------
    # Narrative, story
    # ------------------------------------------------------------------
    def _narrate(self, typing_df, gene_tables, mge_df, qc_df):
        n = len(typing_df)
        if n == 0:
            return
        mrsa = int((typing_df.get('MRSA_Status', pd.Series()) == 'MRSA').sum())
        pct = mrsa / n * 100

        top_st, top_n = '—', 0
        if 'MLST' in typing_df.columns:
            vc = typing_df['MLST'].value_counts()
            if not vc.empty:
                top_st, top_n = vc.idxmax(), int(vc.max())

        top_scc = '—'
        if 'SCCmec_CGE' in typing_df.columns:
            scc = typing_df['SCCmec_CGE'][~typing_df['SCCmec_CGE'].isin(
                ['Not Assigned', '', 'ND', 'Unknown'])]
            if not scc.empty:
                top_scc = scc.value_counts().idxmax()

        top_agr = '—'
        if 'agr_Type' in typing_df.columns:
            agr = typing_df['agr_Type'][~typing_df['agr_Type'].isin(
                ['Not Assigned', '', 'ND', 'Unknown'])]
            if not agr.empty:
                top_agr = agr.value_counts().idxmax()

        n_amr = len(gene_tables.get('amr', pd.DataFrame())) if gene_tables else 0
        n_vir = len(gene_tables.get('virulence', pd.DataFrame())) if gene_tables else 0

        self.narrative = [
            f"Analysed <b>{n}</b> S. aureus genomes.",
            f"<b>{mrsa}</b> MRSA ({pct:.1f}%) · <b>{n - mrsa}</b> MSSA.",
            f"Dominant lineage: <b>ST{top_st}</b> ({top_n} samples).",
            f"Top SCCmec: <b>{top_scc}</b> · Top agr: <b>{top_agr}</b>.",
            f"<b>{n_amr}</b> AMR genes · <b>{n_vir}</b> virulence genes detected.",
        ]

    def _build_story(self, typing_df, gene_tables, gene_matrix, qc_df):
        chapters = []
        n = len(typing_df)
        mrsa = int((typing_df.get('MRSA_Status', pd.Series()) == 'MRSA').sum())

        chapters.append({
            'title': 'Chapter 1 · The Cohort',
            'body': f"You are looking at <strong>{n}</strong> S. aureus genomes. "
                    f"{mrsa} of them are MRSA. That is your starting point.",
        })

        if 'MLST' in typing_df.columns and not typing_df['MLST'].empty:
            vc = typing_df['MLST'].value_counts().head(3)
            summary = ', '.join(f"ST{st} (n={c})" for st, c in vc.items())
            chapters.append({
                'title': 'Chapter 2 · The Dominant Lineages',
                'body': f"The top three sequence types are <strong>{summary}</strong>. "
                        f"These are the clones driving your dataset.",
            })

        if 'SCCmec_CGE' in typing_df.columns:
            scc = typing_df['SCCmec_CGE'][~typing_df['SCCmec_CGE'].isin(
                ['Not Assigned', '', 'ND', 'Unknown'])]
            if not scc.empty:
                top = scc.value_counts().idxmax()
                chapters.append({
                    'title': 'Chapter 3 · The Cassette Landscape',
                    'body': f"The dominant SCCmec cassette is <strong>{top}</strong>. "
                            f"That tells you which resistance lineage is most active here.",
                })

        vfdb = (gene_tables or {}).get('virulence')
        if vfdb is not None and not vfdb.empty:
            top_v = vfdb.nlargest(3, 'Count')
            names = ', '.join(top_v['Gene'].tolist())
            chapters.append({
                'title': 'Chapter 4 · The Virulence Arsenal',
                'body': f"The most prevalent virulence genes are <strong>{names}</strong>. "
                        f"These shape the clinical picture of infections from this cohort.",
            })

        if not gene_matrix.empty:
            n_genes = gene_matrix.shape[1]
            chapters.append({
                'title': 'Chapter 5 · The Resistance Landscape',
                'body': f"We detected <strong>{n_genes}</strong> distinct resistance genes "
                        f"across the cohort. Explore the Resistance tab to see which travel together.",
            })

        if qc_df is not None and not qc_df.empty and 'ANI (%)' in qc_df.columns:
            ani = qc_df['ANI (%)'].dropna()
            if not ani.empty:
                chapters.append({
                    'title': 'Chapter 6 · Species Confirmation',
                    'body': f"All genomes have ANI ≥ <strong>{ani.min():.2f}%</strong> against "
                            f"the S. aureus NCTC 8325 reference — species identity is confirmed.",
                })

        self.story_chapters = chapters

    # ------------------------------------------------------------------
    # Sample table
    # ------------------------------------------------------------------
    def _build_sample_table(self, typing_df: pd.DataFrame) -> str:
        cols = ['Sample'] + [c for c in
                             ('MLST', 'spa_Type', 'agr_Type', 'Capsule_Type',
                              'SCCmec_CGE', 'SCCmec_RPet', 'SCCmec_Subtype',
                              'MRSA_Status')
                             if c in typing_df.columns]
        rows = ''
        for _, r in typing_df.iterrows():
            cells = ''.join(f'<td>{esc(r.get(c, ""))}</td>' for c in cols)
            row_class = 'mrsa-row' if r.get('MRSA_Status') == 'MRSA' else ''
            rows += f'<tr class="{row_class}" data-sample="{esc(r["Sample"])}">{cells}</tr>'
        header = ''.join(f'<th>{esc(c)}</th>' for c in cols)
        return f'''<table class="sample-table" id="sample-table">
            <thead><tr>{header}</tr></thead>
            <tbody>{rows}</tbody>
        </table>'''

    # ------------------------------------------------------------------
    # HTML assembly
    # ------------------------------------------------------------------
    def _assemble_html(self, typing_df: pd.DataFrame, qc_df: pd.DataFrame) -> str:
        has_qc = qc_df is not None and not qc_df.empty
        return f'''<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>StaphScope Interactive Dashboard</title>
<script src="https://cdn.plot.ly/plotly-2.27.0.min.js" charset="utf-8"></script>
<style>{self._css()}</style>
</head>
<body>
{self._header()}
{self._narrative_block()}
{self._tabs(has_qc)}
<main class="dashboard-main">
{self._tab_overview()}
{self._tab_typing()}
{self._tab_qc() if has_qc else ''}
{self._tab_amr()}
{self._tab_virulence()}
{self._tab_mge()}
{self._tab_resistance()}
{self._tab_alerts()}
{self._tab_story()}
{self._tab_compare(typing_df)}
</main>
<footer class="dashboard-footer">
    <p><strong>StaphScope Interactive Dashboard v{TOOL_VERSION}</strong> · Generated {self.meta.get("generated","")}</p>
    <p>Author: Brown Beckley · University of Ghana Medical School · brownbeckley94@gmail.com</p>
    <p>⭐ Star us on GitHub if you find this useful!</p>
</footer>
<div id="filter-popup" class="filter-popup" style="display:none;">
    <div class="filter-popup-header">
        <span id="filter-title">Samples</span>
        <button id="filter-close" class="filter-close">×</button>
    </div>
    <div id="filter-body" class="filter-popup-body"></div>
</div>
<script>{self._js()}</script>
</body>
</html>'''

    def _header(self) -> str:
        m = self.meta
        return f'''<header class="dashboard-header">
    <div class="header-content">
        <h1>🧬 StaphScope Interactive Dashboard</h1>
        <p class="subtitle">Modern visual analytics for <em>Staphylococcus aureus</em> genomics</p>
        <div class="header-stats">
            <div class="header-stat"><div class="num">{m.get('n_samples',0)}</div><div class="lbl">Samples</div></div>
            <div class="header-stat"><div class="num">{m.get('n_mrsa',0)}</div><div class="lbl">MRSA</div></div>
            <div class="header-stat"><div class="num">{m.get('n_mssa',0)}</div><div class="lbl">MSSA</div></div>
        </div>
    </div>
</header>'''

    def _narrative_block(self) -> str:
        if not self.narrative:
            return ''
        items = ''.join(f'<li>{n}</li>' for n in self.narrative)
        return f'''<section class="narrative-bar">
    <h3>📖 What This Dataset Shows</h3>
    <ul>{items}</ul>
</section>'''

    def _tabs(self, has_qc: bool) -> str:
        tabs = [('overview', '📊 Overview'), ('typing', '🧬 Typing')]
        if has_qc:
            tabs.append(('qc', '📏 QC'))
        tabs += [
            ('amr', '💊 AMR'), ('virulence', '🦠 Virulence'),
            ('mge', '📱 MGE'), ('resistance', '🛡️ Resistance'),
            ('alerts', '🚨 Alerts'), ('story', '📖 Story'),
            ('compare', '⚖️ Compare'),
        ]
        buttons = ''.join(
            f'<button class="tab-btn{" active" if i==0 else ""}" data-tab="{t}">{l}</button>'
            for i, (t, l) in enumerate(tabs)
        )
        return f'''<nav class="tab-nav">{buttons}</nav>
<div class="view-filter">
    <label>Filter charts by status:</label>
    <button class="view-btn active" data-view="all">All samples</button>
    <button class="view-btn" data-view="mrsa">MRSA only</button>
    <button class="view-btn" data-view="mssa">MSSA only</button>
</div>'''

    def _tab_overview(self) -> str:
        return f'''<section class="tab-content active" data-tab-content="overview">
    <div class="grid-2">
        {self._chart_block('mrsa_donut', 'MRSA vs MSSA')}
        {self._chart_block('sunburst', 'Hierarchical Lineage — MRSA → SCCmec → ST')}
    </div>
    <div class="grid-full">
        {self._chart_block('sankey', 'Typing Flow — ST → SCCmec → agr → Capsule')}
    </div>
</section>'''

    def _tab_typing(self) -> str:
        return f'''<section class="tab-content" data-tab-content="typing">
    <div class="grid-2">
        {self._chart_block('bar-MLST', 'MLST Sequence Types')}
        {self._chart_block('bar-spa_Type', 'spa Types')}
        {self._chart_block('bar-SCCmec_CGE', 'SCCmec Types (CGE caller)')}
        {self._chart_block('bar-SCCmec_Subtype', 'SCCmec Subtypes')}
        {self._chart_block('bar-agr_Type', 'agr Types')}
        {self._chart_block('bar-Capsule_Type', 'Capsule Types')}
    </div>
</section>'''

    def _tab_qc(self) -> str:
        return f'''<section class="tab-content" data-tab-content="qc">
    <div class="grid-2">
        {self._chart_block('qc_boxplots', 'Assembly Quality Metrics')}
        {self._chart_block('ani_chart', 'Species Confirmation (fastANI)')}
    </div>
</section>'''

    def _tab_amr(self) -> str:
        return f'''<section class="tab-content" data-tab-content="amr">
    <div class="grid-full">
        {self._chart_block('amr_cooccurrence', 'AMR Gene Presence Matrix (Top 30 genes)')}
    </div>
    <div class="grid-full">
        {self._chart_block('amr_db_comparison', 'Gene Frequency Distribution by Category')}
    </div>
</section>'''

    def _tab_virulence(self) -> str:
        return f'''<section class="tab-content" data-tab-content="virulence">
    <div class="grid-full">
        {self._chart_block('virulence_top', 'Top Virulence Genes (VFDB)')}
    </div>
</section>'''

    def _tab_mge(self) -> str:
        return f'''<section class="tab-content" data-tab-content="mge">
    <div class="grid-full">
        {self._chart_block('mge_breakdown', 'Mean mobileOG Hits per Category')}
    </div>
    <div class="grid-full">
        {self._chart_block('radial_rings', 'Radial Typing Rings — one spoke per sample, one ring per typing layer')}
    </div>
</section>'''

    def _tab_resistance(self) -> str:
        return f'''<section class="tab-content" data-tab-content="resistance">
    <div class="grid-full">
        {self._chart_block('resistance_matrix', 'Resistance Presence / Absence Matrix (clustered)')}
    </div>
    <div class="grid-full">
        {self._chart_block('upset', 'Gene Combination Profiles (UpSet-style)')}
    </div>
</section>'''

    def _tab_alerts(self) -> str:
        return f'''<section class="tab-content" data-tab-content="alerts">
    {self._render_alerts()}
</section>'''

    def _render_alerts(self) -> str:
        if not self.alerts:
            return '<div class="no-alerts">✅ No alerts raised — the dataset is clean.</div>'
        groups = defaultdict(list)
        for a in self.alerts:
            groups[a['id']].append(a)
        html = '<div class="alert-panel">'
        for rule_id, items in groups.items():
            rule = ALERT_RULES[rule_id]
            lines = '<br>'.join(
                f'<b>{esc(a["sample"])}</b> — {esc(a["detail"])}'
                for a in items[:15]
            )
            if len(items) > 15:
                lines += f'<br><em>… and {len(items)-15} more</em>'
            html += f'''<div class="alert-card {rule['severity']}">
                <div class="alert-icon">{rule['icon']}</div>
                <div class="alert-body">
                    <h4>{rule['title']} <span class="severity-tag {rule['severity']}">{rule['severity']}</span></h4>
                    <p class="desc">{rule['description']}</p>
                    <div class="detail">{lines}</div>
                </div>
            </div>'''
        html += '</div>'
        return html

    def _tab_story(self) -> str:
        if not self.story_chapters:
            return '<section class="tab-content" data-tab-content="story"><p>No story available.</p></section>'
        cards = ''.join(f'''<div class="story-card">
            <div class="story-number">{i+1}</div>
            <div class="story-body">
                <h4>{esc(ch['title'])}</h4>
                <p>{ch['body']}</p>
            </div>
        </div>''' for i, ch in enumerate(self.story_chapters))
        return f'''<section class="tab-content" data-tab-content="story">
    <div class="story-container">{cards}</div>
</section>'''

    def _tab_compare(self, typing_df: pd.DataFrame) -> str:
        options = ''.join(
            f'<option value="{esc(s)}">{esc(s)}</option>'
            for s in sorted(typing_df['Sample'].tolist())
        ) if 'Sample' in typing_df.columns else ''
        return f'''<section class="tab-content" data-tab-content="compare">
    <div class="compare-controls">
        <label>Sample A:
            <select id="compare-a">{options}</select>
        </label>
        <label>Sample B:
            <select id="compare-b">{options}</select>
        </label>
        <button id="compare-run" class="btn-primary">Compare</button>
    </div>
    <div id="compare-output"></div>
    <details class="sample-table-details">
        <summary>📋 View all samples (click a row to filter other tabs)</summary>
        {self.sample_table_html}
    </details>
</section>'''

    # ------------------------------------------------------------------
    # CSS
    # ------------------------------------------------------------------
    def _css(self) -> str:
        return """
:root{
  --bg:#0f172a;--bg-card:#1e293b;--bg-hover:#334155;--bg-soft:#1a2538;
  --text:#f1f5f9;--text-muted:#94a3b8;--border:#334155;
  --primary:#4CAF50;--accent:#7E22CE;--mrsa:#DC143C;--mssa:#4682B4;
  --critical:#dc2626;--high:#ea580c;--medium:#eab308;
}
*{margin:0;padding:0;box-sizing:border-box}
html,body{background:var(--bg);color:var(--text);
  font-family:'Segoe UI',Tahoma,Geneva,Verdana,sans-serif;line-height:1.55}
body{padding-bottom:40px}

.dashboard-header{background:linear-gradient(135deg,#006400 0%,#228B22 45%,#7E22CE 100%);
  color:white;padding:38px 30px;border-radius:0 0 24px 24px;
  box-shadow:0 10px 40px rgba(0,0,0,.4);margin-bottom:24px}
.header-content{max-width:1600px;margin:0 auto}
.dashboard-header h1{font-size:2.4em;margin-bottom:6px;letter-spacing:-.5px}
.dashboard-header .subtitle{font-size:1.05em;opacity:.92;margin-bottom:20px}
.header-stats{display:flex;gap:16px;flex-wrap:wrap;margin-top:20px}
.header-stat{background:rgba(255,255,255,.15);padding:14px 22px;border-radius:12px;
  backdrop-filter:blur(10px);border:1px solid rgba(255,255,255,.15);min-width:120px}
.header-stat .num{font-size:1.9em;font-weight:700;line-height:1}
.header-stat .lbl{font-size:.78em;opacity:.88;text-transform:uppercase;letter-spacing:.5px;margin-top:4px}

.narrative-bar{max-width:1600px;margin:0 auto 24px;padding:20px 30px;
  background:var(--bg-card);border-radius:14px;border-left:5px solid var(--primary)}
.narrative-bar h3{font-size:1em;color:var(--primary);margin-bottom:10px;
  text-transform:uppercase;letter-spacing:.5px;font-weight:700}
.narrative-bar ul{list-style:none;padding:0;display:flex;flex-wrap:wrap;gap:16px 30px}
.narrative-bar li{font-size:1.02em}
.narrative-bar li b{color:var(--primary)}

.tab-nav{position:sticky;top:0;z-index:100;background:rgba(15,23,42,.95);
  backdrop-filter:blur(12px);padding:12px 30px;display:flex;gap:8px;flex-wrap:wrap;
  border-bottom:1px solid var(--border);margin-bottom:16px}
.tab-btn{padding:10px 20px;background:transparent;border:1px solid var(--border);
  color:var(--text-muted);border-radius:8px;cursor:pointer;font-weight:600;
  font-size:.9em;transition:all .2s}
.tab-btn:hover{background:var(--bg-hover);color:var(--text)}
.tab-btn.active{background:var(--primary);color:white;border-color:var(--primary);
  box-shadow:0 0 20px rgba(76,175,80,.35)}

.view-filter{max-width:1600px;margin:0 auto 24px;padding:0 30px;display:flex;
  gap:10px;align-items:center;flex-wrap:wrap}
.view-filter label{color:var(--text-muted);font-size:.9em;font-weight:600}
.view-btn{padding:8px 18px;background:var(--bg-card);border:1px solid var(--border);
  color:var(--text-muted);border-radius:20px;cursor:pointer;font-weight:600;
  font-size:.85em;transition:all .2s}
.view-btn:hover{color:var(--text)}
.view-btn.active{color:white}
.view-btn[data-view="all"].active{background:var(--primary);border-color:var(--primary)}
.view-btn[data-view="mrsa"].active{background:var(--mrsa);border-color:var(--mrsa)}
.view-btn[data-view="mssa"].active{background:var(--mssa);border-color:var(--mssa)}

.dashboard-main{max-width:1600px;margin:0 auto;padding:0 30px}
.tab-content{display:none;animation:fadeIn .3s}
.tab-content.active{display:block}
@keyframes fadeIn{from{opacity:0;transform:translateY(8px)}to{opacity:1;transform:translateY(0)}}

.grid-2{display:grid;grid-template-columns:repeat(auto-fit,minmax(520px,1fr));
  gap:20px;margin-bottom:20px}
.grid-full{margin-bottom:20px}
.chart-card{background:var(--bg-card);border:1px solid var(--border);
  border-radius:14px;padding:20px;transition:box-shadow .2s}
.chart-card:hover{box-shadow:0 8px 30px rgba(0,0,0,.3)}
.chart-title{font-size:1.05em;margin-bottom:6px;font-weight:600;
  display:flex;align-items:center;gap:8px}
.chart-title::before{content:'';width:4px;height:18px;background:var(--primary);
  border-radius:2px}
.chart-subtitle{font-size:.85em;color:var(--text-muted);margin-bottom:12px}

.alert-panel{display:flex;flex-direction:column;gap:12px}
.alert-card{background:var(--bg-card);border-radius:12px;padding:18px 22px;
  border-left:5px solid;display:flex;gap:16px;align-items:flex-start;
  transition:transform .15s}
.alert-card:hover{transform:translateX(4px)}
.alert-card.critical{border-left-color:var(--critical)}
.alert-card.high{border-left-color:var(--high)}
.alert-card.medium{border-left-color:var(--medium)}
.alert-icon{font-size:1.8em;line-height:1}
.alert-body{flex:1}
.alert-body h4{font-size:1.02em;margin-bottom:4px}
.alert-body .desc{font-size:.88em;color:var(--text-muted);margin-bottom:6px}
.alert-body .detail{font-size:.85em;font-family:'Consolas',monospace;
  background:var(--bg);padding:8px 12px;border-radius:6px;line-height:1.7}
.severity-tag{font-size:.72em;padding:3px 10px;border-radius:12px;
  text-transform:uppercase;letter-spacing:.5px;font-weight:700;margin-left:8px}
.severity-tag.critical{background:rgba(220,38,38,.2);color:#fca5a5}
.severity-tag.high{background:rgba(234,88,12,.2);color:#fdba74}
.severity-tag.medium{background:rgba(234,179,8,.2);color:#fde047}
.no-alerts{background:rgba(76,175,80,.1);border-left:5px solid var(--primary);
  padding:40px;border-radius:12px;text-align:center;font-size:1.1em;color:var(--primary)}

.story-container{display:flex;flex-direction:column;gap:16px}
.story-card{background:var(--bg-card);border:1px solid var(--border);
  border-radius:14px;padding:22px 26px;display:flex;gap:22px;align-items:flex-start}
.story-number{width:44px;height:44px;border-radius:50%;background:var(--primary);
  color:white;display:flex;align-items:center;justify-content:center;
  font-weight:800;font-size:1.2em;flex-shrink:0}
.story-body h4{margin-bottom:8px;font-size:1.1em}
.story-body p{color:var(--text-muted);font-size:.98em}
.story-body b{color:var(--primary)}

.compare-controls{background:var(--bg-card);border-radius:14px;padding:20px;
  display:flex;gap:16px;align-items:center;flex-wrap:wrap;margin-bottom:20px}
.compare-controls label{display:flex;gap:8px;align-items:center;
  color:var(--text-muted);font-weight:600;font-size:.92em}
.compare-controls select{background:var(--bg);color:var(--text);
  border:1px solid var(--border);border-radius:8px;padding:8px 12px;
  font-size:.92em;min-width:220px}
.btn-primary{background:var(--primary);color:white;border:none;
  padding:10px 22px;border-radius:8px;font-weight:700;cursor:pointer;
  transition:transform .2s}
.btn-primary:hover{transform:translateY(-2px)}

#compare-output{display:grid;grid-template-columns:1fr 1fr;gap:20px}
.compare-col{background:var(--bg-card);border-radius:14px;padding:20px;
  border:1px solid var(--border)}
.compare-col h4{color:var(--primary);margin-bottom:12px}
.compare-row{display:flex;justify-content:space-between;padding:8px 0;
  border-bottom:1px dashed var(--border);font-size:.92em}
.compare-row.diff{background:rgba(220,38,38,.08);border-left:3px solid var(--critical);
  padding-left:8px;margin-left:-8px}

.sample-table-details{margin-top:30px;background:var(--bg-card);
  border-radius:14px;padding:20px}
.sample-table-details summary{cursor:pointer;font-weight:700;
  color:var(--primary);font-size:1.02em;padding:6px 0}
.sample-table{width:100%;border-collapse:collapse;margin-top:16px;font-size:.88em}
.sample-table th{background:var(--bg);color:var(--text);padding:10px 12px;
  text-align:left;position:sticky;top:0;font-weight:600;border-bottom:2px solid var(--primary)}
.sample-table td{padding:8px 12px;border-bottom:1px solid var(--border)}
.sample-table tr:hover{background:var(--bg-hover)}
.sample-table tr.mrsa-row{background:rgba(220,38,38,.08)}
.sample-table tr.mrsa-row:hover{background:rgba(220,38,38,.15)}

.filter-popup{position:fixed;bottom:20px;right:20px;width:420px;max-height:60vh;
  background:var(--bg-card);border:1px solid var(--border);border-radius:12px;
  box-shadow:0 20px 60px rgba(0,0,0,.6);z-index:1000;overflow:hidden;
  display:flex;flex-direction:column}
.filter-popup-header{background:var(--primary);color:white;padding:12px 16px;
  display:flex;justify-content:space-between;align-items:center;font-weight:700}
.filter-close{background:transparent;border:none;color:white;font-size:1.5em;
  cursor:pointer;line-height:1}
.filter-popup-body{padding:14px 16px;overflow-y:auto;
  font-size:.85em;font-family:'Consolas',monospace;line-height:1.7}
.filter-popup-body span{display:inline-block;background:var(--bg);
  padding:3px 10px;border-radius:12px;margin:2px;border:1px solid var(--border)}

.dashboard-footer{max-width:1600px;margin:40px auto 0;padding:30px;
  text-align:center;color:var(--text-muted);font-size:.88em;
  border-top:1px solid var(--border)}
.dashboard-footer p{margin:4px 0}

@media (max-width:900px){
  .dashboard-header h1{font-size:1.7em}
  .grid-2{grid-template-columns:1fr}
  #compare-output{grid-template-columns:1fr}
  .dashboard-main,.tab-nav,.view-filter,.narrative-bar{padding-left:15px;padding-right:15px}
}
"""

    # ------------------------------------------------------------------
    # JS
    # ------------------------------------------------------------------
    def _js(self) -> str:
        return """
// ---- Tab switching ----
document.querySelectorAll('.tab-btn').forEach(btn => {
    btn.addEventListener('click', () => {
        const tab = btn.dataset.tab;
        document.querySelectorAll('.tab-btn').forEach(b => b.classList.remove('active'));
        document.querySelectorAll('.tab-content').forEach(c => c.classList.remove('active'));
        btn.classList.add('active');
        const content = document.querySelector(`[data-tab-content="${tab}"]`);
        if (content) content.classList.add('active');
        setTimeout(() => {
            document.querySelectorAll('.js-plotly-plot').forEach(div => {
                if (div.offsetWidth > 0) Plotly.Plots.resize(div);
            });
        }, 60);
    });
});

// ---- View filter ----
document.querySelectorAll('.view-btn').forEach(btn => {
    btn.addEventListener('click', () => {
        const view = btn.dataset.view;
        document.querySelectorAll('.view-btn').forEach(b => b.classList.remove('active'));
        btn.classList.add('active');
        document.querySelectorAll('.chart-view').forEach(v => {
            v.style.display = v.dataset.view === view ? 'block' : 'none';
        });
        setTimeout(() => {
            document.querySelectorAll('.chart-view:not([style*="display: none"]) .js-plotly-plot').forEach(div => {
                Plotly.Plots.resize(div);
            });
        }, 60);
    });
});

// ---- Cross-filter popup ----
const popup = document.getElementById('filter-popup');
const popupTitle = document.getElementById('filter-title');
const popupBody = document.getElementById('filter-body');
if (popup) {
    document.getElementById('filter-close').addEventListener('click', () => {
        popup.style.display = 'none';
    });
}

document.addEventListener('DOMContentLoaded', () => {
    setTimeout(() => {
        document.querySelectorAll('.js-plotly-plot').forEach(div => {
            if (typeof div.on !== 'function') return;
            div.on('plotly_click', (data) => {
                const pt = data.points[0];
                let samples = [];
                if (pt.customdata && typeof pt.customdata === 'string') {
                    samples = pt.customdata.split(';');
                } else if (Array.isArray(pt.customdata) && pt.customdata.length) {
                    samples = [String(pt.customdata[0])];
                }
                if (!samples.length) return;
                popupTitle.textContent = `${pt.label || pt.y || pt.x || 'Samples'} · ${samples.length} sample(s)`;
                popupBody.innerHTML = samples
                    .filter(s => s)
                    .map(s => `<span>${s}</span>`).join('');
                popup.style.display = 'flex';
            });
        });
    }, 500);
});

// ---- Compare tool ----
const compareBtn = document.getElementById('compare-run');
if (compareBtn) {
    compareBtn.addEventListener('click', () => {
        const a = document.getElementById('compare-a').value;
        const b = document.getElementById('compare-b').value;
        if (!a || !b) return;
        const rows = Array.from(document.querySelectorAll('#sample-table tr'));
        const rowA = rows.find(r => r.dataset.sample === a);
        const rowB = rows.find(r => r.dataset.sample === b);
        if (!rowA || !rowB) return;
        const headers = Array.from(rowA.parentElement.parentElement.querySelectorAll('thead th'))
            .map(th => th.textContent);
        const aVals = Array.from(rowA.querySelectorAll('td')).map(td => td.textContent);
        const bVals = Array.from(rowB.querySelectorAll('td')).map(td => td.textContent);

        let leftHtml = `<h4>${a}</h4>`;
        let rightHtml = `<h4>${b}</h4>`;
        let diffCount = 0;
        for (let i = 1; i < headers.length; i++) {
            const diff = aVals[i] !== bVals[i];
            if (diff) diffCount++;
            leftHtml += `<div class="compare-row ${diff ? 'diff' : ''}"><span>${headers[i]}</span><span>${aVals[i]}</span></div>`;
            rightHtml += `<div class="compare-row ${diff ? 'diff' : ''}"><span>${headers[i]}</span><span>${bVals[i]}</span></div>`;
        }
        document.getElementById('compare-output').innerHTML = `
            <div class="compare-col">${leftHtml}</div>
            <div class="compare-col">${rightHtml}</div>`;
        if (diffCount === 0) {
            document.getElementById('compare-output').insertAdjacentHTML('afterbegin',
                `<div style="grid-column:1/-1;padding:14px 20px;background:rgba(76,175,80,.15);border-left:4px solid var(--primary);border-radius:8px;margin-bottom:12px;">
                    ✅ These samples have identical typing profiles — possible transmission pair.
                </div>`);
        }
    });
}
"""


# ==============================================================================
# STATIC PLOTTER
# ==============================================================================
class StaticPlotter:
    """Publication-quality matplotlib figures."""

    def __init__(self, output_dir: Path):
        self.output_dir = Path(output_dir)
        self.subdirs = {
            'png': self.output_dir / 'PNG',
            'pdf': self.output_dir / 'PDF',
            'svg': self.output_dir / 'SVG',
            'data': self.output_dir / 'DATA',
        }
        for d in self.subdirs.values():
            d.mkdir(parents=True, exist_ok=True)

    def _save(self, fig, name: str):
        for fmt in ('png', 'pdf', 'svg'):
            try:
                fig.savefig(self.subdirs[fmt] / f"{name}.{fmt}",
                            dpi=300 if fmt == 'png' else None,
                            bbox_inches='tight')
            except Exception as e:
                print(f"    ⚠️ Save {name}.{fmt} failed: {e}")
        plt.close(fig)

    @staticmethod
    def _wrap_labels(labels, max_len=30):
        return [l if len(str(l)) <= max_len else str(l)[:max_len - 1] + '…' for l in labels]

    def typing_summary(self, typing_df: pd.DataFrame):
        if typing_df.empty:
            return
        fig, axes = plt.subplots(2, 2, figsize=(18, 14))

        # MRSA donut
        ax = axes[0, 0]
        if 'MRSA_Status' in typing_df.columns:
            counts = typing_df['MRSA_Status'].value_counts()
            counts = counts[counts.index.isin(['MRSA', 'MSSA'])]
            if not counts.empty:
                colors = [PALETTE['mrsa'] if l == 'MRSA' else PALETTE['mssa']
                          for l in counts.index]
                ax.pie(counts.values, labels=counts.index, colors=colors,
                       autopct='%1.1f%%', startangle=90,
                       wedgeprops=dict(width=0.4, edgecolor='white'))
                ax.set_title('MRSA vs MSSA', fontweight='bold', fontsize=13)
                pd.DataFrame({'Status': counts.index, 'Count': counts.values}) \
                    .to_csv(self.subdirs['data'] / 'mrsa_distribution.csv', index=False)

        # MLST
        ax = axes[0, 1]
        if 'MLST' in typing_df.columns:
            vc = typing_df['MLST'].value_counts().head(15)
            if not vc.empty:
                labels = self._wrap_labels([f"ST{s}" for s in vc.index])
                ax.barh(labels, vc.values, color=PALETTE['mlst'])
                ax.invert_yaxis()
                ax.set_xlabel('Samples')
                ax.set_title(f'MLST Top {len(vc)}', fontweight='bold', fontsize=13)
                ax.grid(alpha=0.3, axis='x')
                ax.tick_params(axis='y', labelsize=9)

        # SCCmec
        ax = axes[1, 0]
        if 'SCCmec_CGE' in typing_df.columns:
            scc = typing_df['SCCmec_CGE'][~typing_df['SCCmec_CGE'].isin(
                ['Not Assigned', '', 'ND', 'Unknown'])]
            if not scc.empty:
                vc = scc.value_counts().head(15)
                labels = self._wrap_labels(vc.index, max_len=28)
                ax.barh(labels, vc.values, color=PALETTE['sccmec'])
                ax.invert_yaxis()
                ax.set_xlabel('Samples')
                ax.set_title('SCCmec Types', fontweight='bold', fontsize=13)
                ax.grid(alpha=0.3, axis='x')
                ax.tick_params(axis='y', labelsize=9)

        # agr
        ax = axes[1, 1]
        if 'agr_Type' in typing_df.columns:
            agr = typing_df['agr_Type'][~typing_df['agr_Type'].isin(
                ['Not Assigned', '', 'ND', 'Unknown'])]
            if not agr.empty:
                vc = agr.value_counts()
                ax.bar(vc.index.astype(str), vc.values, color=PALETTE['agr'])
                ax.set_ylabel('Samples')
                ax.set_title('agr Types', fontweight='bold', fontsize=13)
                ax.grid(alpha=0.3, axis='y')

        plt.suptitle('S. aureus Typing Summary', fontsize=16, fontweight='bold')
        plt.tight_layout()
        self._save(fig, 'typing_summary')
        print("    ✓ typing_summary.[png/pdf/svg]")

    def gene_top_plot(self, gene_df: pd.DataFrame, title: str, filename: str, color: str):
        if gene_df is None or gene_df.empty:
            return
        top = gene_df.nlargest(20, 'Count')
        if top.empty:
            return
        labels = self._wrap_labels(top['Gene'].astype(str), max_len=40)
        height = max(6, 0.4 * len(top))
        fig, ax = plt.subplots(figsize=(11, height))
        ax.barh(labels, top['Count'], color=color)
        ax.invert_yaxis()
        ax.set_xlabel('Number of genomes')
        ax.set_title(title, fontweight='bold', fontsize=13)
        ax.grid(alpha=0.3, axis='x')
        ax.tick_params(axis='y', labelsize=9)
        plt.tight_layout()
        self._save(fig, filename)
        print(f"    ✓ {filename}.[png/pdf/svg]")

    def qc_boxplots(self, qc_df: pd.DataFrame):
        if qc_df is None or qc_df.empty:
            return
        metrics = [c for c in ('N50', 'Total Length', 'GC Content (%)',
                               'Total Sequences', 'Ambiguous Bases (%)')
                   if c in qc_df.columns]
        if not metrics:
            return
        n = len(metrics)
        cols = min(3, n)
        rows = (n + cols - 1) // cols
        fig, axes = plt.subplots(rows, cols, figsize=(6 * cols, 5 * rows))
        if n == 1:
            axes = np.array([axes])
        axes = axes.flatten()
        for i, col in enumerate(metrics):
            data = qc_df[col].dropna()
            if data.empty:
                continue
            ax = axes[i]
            ax.boxplot(data.values, vert=True, patch_artist=True,
                       boxprops=dict(facecolor=PALETTE['primary'], alpha=0.6),
                       medianprops=dict(color='black', linewidth=2))
            ax.scatter(np.random.normal(1, 0.04, len(data)), data.values,
                       alpha=0.6, s=20, color='#334155')
            ax.set_ylabel(col, fontsize=11)
            ax.set_title(col, fontweight='bold', fontsize=12)
            ax.grid(alpha=0.3, axis='y')
        for j in range(i + 1, len(axes)):
            axes[j].axis('off')
        plt.suptitle('Assembly QC Metrics', fontsize=15, fontweight='bold')
        plt.tight_layout()
        self._save(fig, 'qc_boxplots')
        print("    ✓ qc_boxplots.[png/pdf/svg]")


# ==============================================================================
# ORCHESTRATOR
# ==============================================================================
class StaphScopeVisualizer:
    """Coordinates parsing, dashboard, and static exports."""

    def __init__(self, input_dir: Path, output_dir: Optional[Path] = None):
        self.input_dir = Path(input_dir)
        self.output_dir = Path(output_dir) if output_dir else \
            self.input_dir / 'STAPHSCOPE_VISUALIZATIONS'
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.loader = StaphDataLoader(self.input_dir)
        self.plotter = StaticPlotter(self.output_dir)
        self.builder = DashboardBuilder(self.output_dir)

        self.typing_df = pd.DataFrame()
        self.gene_tables: Dict[str, pd.DataFrame] = {}
        self.mutations_df = pd.DataFrame()
        self.mge_df = pd.DataFrame()
        self.qc_df = pd.DataFrame()
        self.gene_matrix = pd.DataFrame()

    def run(self):
        print("=" * 70)
        print(f"🧬 STAPHSCOPE VISUALIZER v{TOOL_VERSION}")
        print("=" * 70)
        print(f"Input:  {self.input_dir}")
        print(f"Output: {self.output_dir}")
        print("-" * 70)

        start = datetime.now()

        print("\n📥 Loading data...")
        self.typing_df = self.loader.load_typing()
        self.gene_tables = self.loader.load_gene_tables()
        self.mutations_df = self.loader.load_mutations()
        self.mge_df = self.loader.load_mge()
        self.qc_df = self.loader.load_qc()

        if self.typing_df.empty and not self.gene_tables:
            print("❌ No data found. Nothing to visualize.")
            return

        self.gene_matrix = build_gene_matrix(self.gene_tables)

        print("\n📊 Generating static publication plots...")
        try:
            self.plotter.typing_summary(self.typing_df)
        except Exception as e:
            print(f"  ⚠️ typing_summary: {e}")
        try:
            self.plotter.gene_top_plot(self.gene_tables.get('amr'),
                                       'Top 20 AMR Genes', 'amr_top_genes',
                                       PALETTE['amr'])
        except Exception as e:
            print(f"  ⚠️ amr_top_genes: {e}")
        try:
            self.plotter.gene_top_plot(self.gene_tables.get('virulence'),
                                       'Top 20 Virulence Genes', 'virulence_top_genes',
                                       PALETTE['virulence'])
        except Exception as e:
            print(f"  ⚠️ virulence_top_genes: {e}")
        try:
            self.plotter.qc_boxplots(self.qc_df)
        except Exception as e:
            print(f"  ⚠️ qc_boxplots: {e}")

        if PLOTLY_AVAILABLE:
            print("\n🌐 Building interactive dashboard...")
            try:
                self.builder.build(
                    self.typing_df, self.gene_tables,
                    self.mge_df, self.qc_df)
            except Exception as e:
                print(f"  ⚠️ Dashboard failed: {e}")
                import traceback
                traceback.print_exc()
        else:
            print("\n⚠️ Plotly not available — skipping interactive dashboard")
            print("   Install with: pip install 'plotly>=5.0'")

        self._write_summary_report(datetime.now() - start)
        zip_path = self._build_export_bundle()
        if zip_path:
            print(f"\n📦 Export bundle: {zip_path.name}")

        print("\n" + "=" * 70)
        print("✅ VISUALIZATION PIPELINE COMPLETE")
        print("=" * 70)

    def _write_summary_report(self, duration):
        report = self.output_dir / 'staphscope_visualization_report.txt'
        with open(report, 'w', encoding='utf-8') as f:
            f.write("=" * 70 + "\n")
            f.write("STAPHSCOPE VISUALIZATION REPORT\n")
            f.write("=" * 70 + "\n\n")
            f.write(f"Generated: {datetime.now().isoformat()}\n")
            f.write(f"Duration: {duration}\n")
            f.write(f"Input:  {self.input_dir}\n")
            f.write(f"Output: {self.output_dir}\n\n")

            f.write("DATA SUMMARY\n" + "-" * 40 + "\n")
            f.write(f"Samples: {len(self.typing_df)}\n")
            if not self.typing_df.empty and 'MRSA_Status' in self.typing_df.columns:
                f.write(f"MRSA: {(self.typing_df['MRSA_Status'] == 'MRSA').sum()}\n")
            for cat, df in self.gene_tables.items():
                f.write(f"{cat.capitalize()} genes: {len(df)}\n")
            f.write(f"AMR matrix genes: {self.gene_matrix.shape[1] if not self.gene_matrix.empty else 0}\n")
            if not self.mge_df.empty:
                f.write(f"MGE samples: {len(self.mge_df)}\n")
            if not self.qc_df.empty:
                f.write(f"QC samples: {len(self.qc_df)}\n")

            if self.builder.alerts:
                f.write("\nALERTS\n" + "-" * 40 + "\n")
                for a in self.builder.alerts:
                    f.write(f"  [{ALERT_RULES[a['id']]['severity'].upper()}] "
                            f"{a['sample']}: {a['detail']}\n")

            f.write("\nOUTPUT FILES\n" + "-" * 40 + "\n")
            for d in ('PNG', 'PDF', 'SVG', 'DATA'):
                sub = self.output_dir / d
                if sub.exists():
                    files = list(sub.glob("*"))
                    f.write(f"\n{d} ({len(files)} files):\n")
                    for fp in sorted(files):
                        f.write(f"  • {fp.name}\n")

            dashboard = self.output_dir / 'staphscope_dashboard.html'
            if dashboard.exists():
                f.write("\nINTERACTIVE DASHBOARD\n" + "-" * 40 + "\n")
                f.write(f"  {dashboard.name}\n")
                f.write(f"  Open in browser: file://{dashboard.resolve()}\n")

        print(f"\n📋 Report: {report.name}")

    def _build_export_bundle(self) -> Optional[Path]:
        try:
            zip_path = self.output_dir.parent / 'staphscope_visualizations_bundle.zip'
            with zipfile.ZipFile(zip_path, 'w', zipfile.ZIP_DEFLATED) as zf:
                for file in self.output_dir.rglob('*'):
                    if file.is_file() and file != zip_path:
                        zf.write(file, file.relative_to(self.output_dir.parent))
            return zip_path
        except Exception as e:
            print(f"  ⚠️ Export bundle failed: {e}")
            return None


# ==============================================================================
# CLI
# ==============================================================================
def main():
    parser = argparse.ArgumentParser(
        description="STAPHSCOPE Visualizer — unified interactive + static report",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python staphscope_visualizer.py
  python staphscope_visualizer.py --input /path/to/results --output /path/to/vis
  python staphscope_visualizer.py --no-static      # dashboard only
  python staphscope_visualizer.py --no-dashboard   # static plots only
        """,
    )
    parser.add_argument('--input', '-i', type=str, default='.',
                        help='Directory containing StaphScope outputs (default: cwd)')
    parser.add_argument('--output', '-o', type=str, default=None,
                        help='Output directory (default: STAPHSCOPE_VISUALIZATIONS)')
    parser.add_argument('--no-static', action='store_true',
                        help='Skip static matplotlib exports')
    parser.add_argument('--no-dashboard', action='store_true',
                        help='Skip interactive dashboard')
    args = parser.parse_args()

    viz = StaphScopeVisualizer(
        input_dir=Path(args.input),
        output_dir=Path(args.output) if args.output else None,
    )
    if args.no_static:
        viz.plotter.typing_summary = lambda *a, **k: None
        viz.plotter.gene_top_plot = lambda *a, **k: None
        viz.plotter.qc_boxplots = lambda *a, **k: None
    if args.no_dashboard:
        viz.builder.build = lambda *a, **k: None

    try:
        viz.run()
    except Exception as e:
        print(f"\n❌ Error: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == '__main__':
    main()
#!/usr/bin/env python3
"""
StaphScope Comprehensive Report Generator
Creates unified report from MLST, spa, SCCmec (CGE + RPet), capsule, and agr typing results
Author: Brown Beckley <brownbeckley94@gmail.com>
Date: 2026-09-10 (Robust: any subset of input files is accepted)
Affiliation: University of Ghana Medical School - Department of Medical Biochemistry
Please kindly send a quick email for any technical issues or assistance
"""

import os
import sys
import json
import random
import argparse
import pandas as pd
from pathlib import Path
from collections import Counter
from datetime import datetime

# Science quotes for rotation
SCIENCE_QUOTES = [
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
    {"text": "Thank you for using STAPHSCOPE and don't forget to share with us any feature suggestions in mind.", "author": "Brown Beckley"},
    {"text": "Staphscope represents the convergence of genomic surveillance and clinical diagnostics, transforming raw sequences into actionable insights for infection control.", "author": "Brown Beckley"},
    {"text": "In the battle against antimicrobial resistance, tools like Staphscope are our eyes and ears, revealing the genetic blueprints of resistant pathogens.", "author": "Brown Beckley"},
    {"text": "Staphscope isn't just a tool; it's a comprehensive system that bridges the gap between sequencing data and public health action.", "author": "Brown Beckley"},
    {"text": "Through Staphscope, we turn the complexity of bacterial genomes into clear, interpretable reports, empowering clinicians and researchers alike.", "author": "Brown Beckley"},
    {"text": "Staphscope is a testament to the power of bioinformatics in the modern era, making advanced pathogen typing accessible to all.", "author": "Brown Beckley"}
]

# Values we treat as "no data"
NULL_VALUES = {'-', 'ND', '', 'UNKNOWN', 'Not Assigned', 'Not Assigned', 'None', 'nan', 'NaN'}

# All known input files (label → filename). Every entry is optional.
INPUT_FILES = {
    'mlst':        'mlst_summary.tsv',
    'spa':         'spa_summary.tsv',
    'sccmec_cge':  'staphscope_sccmec_cge_summary.tsv',
    'sccmec_rpet': 'staphscope_sccmec_rpet_summary.tsv',
    'capsule':     'staphscope_capsule_summary.tsv',
    'agr':         'agr_summary.tsv',
}


def normalize_sample_name(sample_name):
    """Remove FASTA extension if present for consistent matching across modules."""
    if not sample_name:
        return ""
    s = str(sample_name).strip()
    for ext in ('.fna', '.fasta', '.fa', '.fn'):
        if s.endswith(ext):
            return s[:-len(ext)]
    return s


def _clean(value, default="Not Assigned"):
    """Return a cleaned string value, defaulting when the value is empty / unknown."""
    if value is None:
        return default
    v = str(value).strip()
    if not v or v in NULL_VALUES:
        return default
    return v


def _safe_read_tsv(path):
    """Read a TSV file, returning an empty DataFrame on any error."""
    try:
        df = pd.read_csv(path, sep='\t', dtype=str)
        return df
    except Exception as e:
        print(f"   ✗ Error reading {path}: {e}")
        return pd.DataFrame()


# ---------------------------------------------------------------------------
# Loaders — all return {} if the file is missing or malformed
# ---------------------------------------------------------------------------
def load_mlst_data(mlst_file):
    print(f"📊 Loading MLST data from: {mlst_file}")
    if not Path(mlst_file).exists():
        print("   ⚠️  File not found — MLST column will be 'Not Assigned'")
        return {}
    df = _safe_read_tsv(mlst_file)
    if df.empty:
        return {}
    mlst_data = {}
    for _, row in df.iterrows():
        sample = normalize_sample_name(row.get('Sample', row.get('sample', '')))
        if not sample:
            continue
        st = _clean(row.get('ST'), "Not Assigned")
        mlst_data[sample] = f"ST{st}" if st != "Not Assigned" else "Not Assigned"
    print(f"   ✓ Loaded {len(mlst_data)} samples")
    return mlst_data


def load_spa_data(spa_file):
    print(f"📊 Loading spa data from: {spa_file}")
    if not Path(spa_file).exists():
        print("   ⚠️  File not found — spa column will be 'Not Assigned'")
        return {}
    df = _safe_read_tsv(spa_file)
    if df.empty:
        return {}
    spa_data = {}
    for _, row in df.iterrows():
        sample = normalize_sample_name(row.get('Sample', row.get('sample', '')))
        if not sample:
            continue
        spa_data[sample] = _clean(row.get('spa_Type'), "Not Assigned")
    print(f"   ✓ Loaded {len(spa_data)} samples")
    return spa_data


def load_sccmec_cge_data(sccmec_file):
    print(f"📊 Loading SCCmec (CGE) data from: {sccmec_file}")
    if not Path(sccmec_file).exists():
        print("   ⚠️  File not found — SCCmec (CGE) column will be 'Not Assigned'")
        return {}
    df = _safe_read_tsv(sccmec_file)
    if df.empty:
        return {}
    sccmec_data = {}
    for _, row in df.iterrows():
        sample = normalize_sample_name(row.get('Sample_Name', row.get('Sample', row.get('sample', ''))))
        if not sample:
            continue
        sccmec_data[sample] = {
            'sccmec_type_cge': _clean(row.get('SCCmec_Type'), "Not Assigned"),
            'mrsa_status_cge': _clean(row.get('MRSA_Status'), "Unknown"),
        }
    print(f"   ✓ Loaded {len(sccmec_data)} samples")
    return sccmec_data


def load_sccmec_rpet_data(rpet_file):
    print(f"📊 Loading SCCmec (RPet) data from: {rpet_file}")
    if not Path(rpet_file).exists():
        print("   ⚠️  File not found — SCCmec (RPet) / subtype columns will be 'Not Assigned'")
        return {}
    df = _safe_read_tsv(rpet_file)
    if df.empty:
        return {}
    rpet_data = {}
    for _, row in df.iterrows():
        sample = normalize_sample_name(row.get('sample', row.get('Sample', '')))
        if not sample:
            continue
        rpet_data[sample] = {
            'sccmec_type_rpet': _clean(row.get('sccmec_type'), "Not Assigned"),
            'sccmec_subtype': _clean(row.get('sccmec_subtype'), "Not Assigned"),
            'mrsa_status_rpet': _clean(row.get('mrsa_status'), "Unknown"),
        }
    print(f"   ✓ Loaded {len(rpet_data)} samples")
    return rpet_data


def load_capsule_data(capsule_file):
    print(f"📊 Loading capsule data from: {capsule_file}")
    if not Path(capsule_file).exists():
        print("   ⚠️  File not found — capsule column will be 'Not Assigned'")
        return {}
    df = _safe_read_tsv(capsule_file)
    if df.empty:
        return {}
    capsule_data = {}
    for _, row in df.iterrows():
        sample = normalize_sample_name(row.get('sample', row.get('Sample', '')))
        if not sample:
            continue
        capsule_data[sample] = _clean(row.get('cap_type'), "Not Assigned")
    print(f"   ✓ Loaded {len(capsule_data)} samples")
    return capsule_data


def load_agr_data(agr_file):
    print(f"📊 Loading agr data from: {agr_file}")
    if not Path(agr_file).exists():
        print("   ⚠️  File not found — agr column will be 'Not Assigned'")
        return {}
    df = _safe_read_tsv(agr_file)
    if df.empty:
        return {}
    agr_data = {}
    for _, row in df.iterrows():
        sample = normalize_sample_name(row.get('Sample', row.get('sample', '')))
        if not sample:
            continue
        agr_data[sample] = _clean(row.get('agr_Type'), "Not Assigned")
    print(f"   ✓ Loaded {len(agr_data)} samples")
    return agr_data


# ---------------------------------------------------------------------------
# Combine
# ---------------------------------------------------------------------------
def combine_data(mlst_data, spa_data, sccmec_cge_data, sccmec_rpet_data,
                 capsule_data, agr_data):
    """Merge all typing sources keyed by normalized sample name."""
    print("\n🔗 Combining data from all sources...")

    all_samples = set()
    for d in (mlst_data, spa_data, sccmec_cge_data, sccmec_rpet_data, capsule_data, agr_data):
        all_samples.update(d.keys())

    if not all_samples:
        print("   ⚠️  No samples found in any source")
        return []

    combined = []
    for sample in sorted(all_samples):
        cge = sccmec_cge_data.get(sample, {}) or {}
        rpet = sccmec_rpet_data.get(sample, {}) or {}

        # MRSA status: prefer RPet (has subtype context), fall back to CGE
        mrsa_status = rpet.get('mrsa_status_rpet') or cge.get('mrsa_status_cge') or "Unknown"
        if mrsa_status not in ("MRSA", "MSSA"):
            mrsa_status = "Unknown"

        combined.append({
            'sample': sample,
            'mlst': mlst_data.get(sample, "Not Assigned"),
            'spa_type': spa_data.get(sample, "Not Assigned"),
            'agr_type': agr_data.get(sample, "Not Assigned"),
            'capsule_type': capsule_data.get(sample, "Not Assigned"),
            'sccmec_type_cge': cge.get('sccmec_type_cge', "Not Assigned"),
            'sccmec_type_rpet': rpet.get('sccmec_type_rpet', "Not Assigned"),
            'sccmec_subtype': rpet.get('sccmec_subtype', "Not Assigned"),
            'mrsa_status': mrsa_status,
        })

    print(f"   ✓ Combined data for {len(combined)} samples")
    return combined


def compute_distributions(data):
    """Compute unique counts and frequency of each categorical field."""
    fields = {
        'agr_type': 'agr Type',
        'capsule_type': 'Capsule Type',
        'sccmec_type_cge': 'SCCmec Type (CGE)',
        'sccmec_type_rpet': 'SCCmec Type (RPet)',
        'sccmec_subtype': 'SCCmec Subtype',
        'mlst': 'MLST',
        'spa_type': 'spa Type',
    }
    distributions = {}
    for key, label in fields.items():
        counter = Counter()
        not_typed = 0
        for row in data:
            v = row.get(key, "Not Assigned")
            if v in ("Not Assigned", "Unknown", "Not Assigned", "None"):
                not_typed += 1
            else:
                counter[v] += 1
        ordered = sorted(counter.items(), key=lambda kv: (-kv[1], kv[0]))
        distributions[key] = {
            'label': label,
            'unique_count': len(ordered),
            'counts': [{'value': v, 'count': c} for v, c in ordered],
            'not_typed': not_typed,
        }
    return distributions


# ---------------------------------------------------------------------------
# Reports
# ---------------------------------------------------------------------------
TSV_COLUMNS = [
    ('sample', 'Sample'),
    ('mlst', 'MLST'),
    ('spa_type', 'spa Type'),
    ('agr_type', 'agr Type'),
    ('capsule_type', 'Capsule Type'),
    ('sccmec_type_cge', 'SCCmec Type (CGE)'),
    ('sccmec_type_rpet', 'SCCmec Type (RPet)'),
    ('sccmec_subtype', 'SCCmec Subtype'),
    ('mrsa_status', 'MRSA/MSSA Status'),
]


def generate_tsv_report(data, output_file):
    print(f"\n📄 Generating TSV report: {output_file}")
    try:
        df = pd.DataFrame(data)
        df = df[[k for k, _ in TSV_COLUMNS]]
        df.columns = [v for _, v in TSV_COLUMNS]
        df.to_csv(output_file, sep='\t', index=False)
        print("   ✓ TSV report generated successfully")
        return True
    except Exception as e:
        print(f"   ✗ Error generating TSV: {e}")
        return False


def generate_json_report(data, distributions, missing_sources, output_file):
    print(f"📄 Generating JSON report: {output_file}")
    try:
        total_samples = len(data)
        mrsa_count = sum(1 for d in data if d['mrsa_status'] == 'MRSA')
        mssa_count = sum(1 for d in data if d['mrsa_status'] == 'MSSA')
        unknown_count = total_samples - mrsa_count - mssa_count

        report_data = {
            'metadata': {
                'generated_date': datetime.now().isoformat(),
                'total_samples': total_samples,
                'report_type': 'StaphScope Comprehensive Report',
                'version': '2.1',
                'author': 'Brown Beckley',
                'email': 'brownbeckley94@gmail.com',
                'affiliation': 'University of Ghana Medical School - Department of Medical Biochemistry',
                'missing_sources': missing_sources,
            },
            'statistics': {
                'mrsa_count': mrsa_count,
                'mrsa_percentage': (mrsa_count / total_samples * 100) if total_samples else 0,
                'mssa_count': mssa_count,
                'mssa_percentage': (mssa_count / total_samples * 100) if total_samples else 0,
                'unknown_count': unknown_count,
                'samples_with_mlst': sum(1 for d in data if d['mlst'] != 'Not Assigned'),
                'samples_with_spa': sum(1 for d in data if d['spa_type'] != 'Not Assigned'),
                'samples_with_agr': sum(1 for d in data if d['agr_type'] != 'Not Assigned'),
                'samples_with_capsule': sum(1 for d in data if d['capsule_type'] != 'Not Assigned'),
                'samples_with_sccmec_cge': sum(1 for d in data if d['sccmec_type_cge'] != 'Not Assigned'),
                'samples_with_sccmec_rpet': sum(1 for d in data if d['sccmec_type_rpet'] != 'Not Assigned'),
                'samples_with_sccmec_subtype': sum(1 for d in data if d['sccmec_subtype'] != 'Not Assigned'),
                'unique_agr_types': distributions['agr_type']['unique_count'],
                'unique_capsule_types': distributions['capsule_type']['unique_count'],
                'unique_sccmec_types_cge': distributions['sccmec_type_cge']['unique_count'],
                'unique_sccmec_types_rpet': distributions['sccmec_type_rpet']['unique_count'],
                'unique_sccmec_subtypes': distributions['sccmec_subtype']['unique_count'],
            },
            'distributions': distributions,
            'samples': data,
        }

        with open(output_file, 'w', encoding='utf-8') as f:
            json.dump(report_data, f, indent=2, ensure_ascii=False)
        print("   ✓ JSON report generated successfully")
        return True
    except Exception as e:
        print(f"   ✗ Error generating JSON: {e}")
        return False


def _build_distribution_cards(distributions):
    palettes = {
        'agr_type': ('#8b5cf6', '#6d28d9'),
        'capsule_type': ('#14b8a6', '#0f766e'),
        'sccmec_type_cge': ('#f59e0b', '#b45309'),
        'sccmec_type_rpet': ('#ec4899', '#be185d'),
        'sccmec_subtype': ('#06b6d4', '#0e7490'),
        'mlst': ('#3b82f6', '#1e40af'),
        'spa_type': ('#10b981', '#047857'),
    }
    cards = []
    for key in ('agr_type', 'capsule_type', 'sccmec_type_cge',
                'sccmec_type_rpet', 'sccmec_subtype'):
        d = distributions.get(key)
        if not d:
            continue
        top, bottom = palettes.get(key, ('#64748b', '#334155'))
        total_typed = sum(c['count'] for c in d['counts'])
        max_count = max((c['count'] for c in d['counts']), default=1)

        if not d['counts']:
            bars_html = '<div class="dist-empty">No typed samples</div>'
        else:
            bars_html = ''
            for c in d['counts']:
                pct_of_typed = (c['count'] / total_typed * 100) if total_typed else 0
                bar_width = (c['count'] / max_count * 100) if max_count else 0
                bars_html += f'''
                    <div class="dist-row">
                        <div class="dist-label" title="{c['value']}">{c['value']}</div>
                        <div class="dist-bar-track">
                            <div class="dist-bar-fill"
                                 style="width: {bar_width:.1f}%;
                                        background: linear-gradient(90deg, {top} 0%, {bottom} 100%);">
                                <span class="dist-bar-count">{c['count']} ({pct_of_typed:.1f}%)</span>
                            </div>
                        </div>
                    </div>'''

        not_typed_html = ''
        if d['not_typed'] > 0:
            not_typed_html = f'<div class="dist-footnote">⚠️ Not Assigned / unknown: <strong>{d["not_typed"]}</strong></div>'

        cards.append(f'''
        <div class="dist-card">
            <div class="dist-header" style="background: linear-gradient(135deg, {top} 0%, {bottom} 100%);">
                <div class="dist-title">{d['label']}</div>
                <div class="dist-unique">
                    <span class="dist-unique-num">{d['unique_count']}</span>
                    <span class="dist-unique-lbl">unique</span>
                </div>
            </div>
            <div class="dist-body">
                {bars_html}
                {not_typed_html}
            </div>
        </div>''')

    return '<div class="dist-grid">' + ''.join(cards) + '</div>'


def generate_html_report(data, distributions, missing_sources, output_file):
    print(f"📄 Generating HTML report: {output_file}")
    try:
        random_quote = random.choice(SCIENCE_QUOTES)

        total_samples = len(data)
        mrsa_count = sum(1 for d in data if d['mrsa_status'] == 'MRSA')
        mssa_count = sum(1 for d in data if d['mrsa_status'] == 'MSSA')
        unknown_count = total_samples - mrsa_count - mssa_count

        table_rows = ""
        for row in data:
            status = row['mrsa_status']
            if status == 'MRSA':
                mrsa_class = 'mrsa-positive'
            elif status == 'MSSA':
                mrsa_class = 'mrsa-negative'
            else:
                mrsa_class = 'unknown-status'
            table_rows += f'''
                <tr data-mrsa="{status}">
                    <td class="col-sample"><strong>{row['sample']}</strong></td>
                    <td>{row['mlst']}</td>
                    <td>{row['spa_type']}</td>
                    <td>{row['agr_type']}</td>
                    <td>{row['capsule_type']}</td>
                    <td>{row['sccmec_type_cge']}</td>
                    <td>{row['sccmec_type_rpet']}</td>
                    <td>{row['sccmec_subtype']}</td>
                    <td class="{mrsa_class}">{status}</td>
                </tr>'''

        dist_cards_html = _build_distribution_cards(distributions)

        s = distributions
        unique_agr = s['agr_type']['unique_count']
        unique_cap = s['capsule_type']['unique_count']
        unique_scc_rpet = s['sccmec_type_rpet']['unique_count']
        unique_scc_sub = s['sccmec_subtype']['unique_count']

        # Missing-sources banner
        missing_banner = ''
        if missing_sources:
            items = ''.join(f'<li><code>{fname}</code> — <em>{label}</em></li>'
                            for label, fname in missing_sources)
            missing_banner = f'''
        <div class="missing-banner">
            <div class="missing-banner-title">⚠️ Some input files were not found — those columns are shown as <code>Not Assigned</code></div>
            <ul class="missing-banner-list">{items}</ul>
            <div class="missing-banner-hint">Run the corresponding module to produce them, then re-run this report.</div>
        </div>'''

        html_content = f'''<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>StaphScope Comprehensive Report</title>
    <style>
        * {{ margin: 0; padding: 0; box-sizing: border-box; }}

        body {{
            background: linear-gradient(135deg, #1e3c72 0%, #2a5298 50%, #7e22ce 100%);
            font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif;
            color: #ffffff;
            padding: 20px;
            min-height: 100vh;
        }}

        .container {{ max-width: 1500px; margin: 0 auto; }}

        .header {{ text-align: center; margin-bottom: 30px; }}

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

        .quote-text {{ font-size: 18px; font-style: italic; margin-bottom: 10px; }}
        .quote-author {{ font-size: 14px; color: #fbbf24; font-weight: bold; }}

        .report-section {{
            background: rgba(255, 255, 255, 0.97);
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

        .missing-banner {{
            background: #fef3c7;
            border-left: 5px solid #f59e0b;
            border-radius: 8px;
            padding: 14px 18px;
            margin-bottom: 20px;
            color: #78350f;
            font-size: 13px;
        }}

        .missing-banner-title {{ font-weight: 700; margin-bottom: 6px; }}
        .missing-banner-list {{ margin: 6px 0 6px 22px; }}
        .missing-banner-list code {{
            background: rgba(120, 53, 15, 0.15);
            padding: 1px 6px;
            border-radius: 4px;
            font-size: 12px;
        }}
        .missing-banner-hint {{ font-size: 12px; opacity: 0.8; margin-top: 4px; }}

        .stats-grid {{
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(180px, 1fr));
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
            transition: transform 0.2s ease;
        }}

        .stat-card:hover {{ transform: translateY(-3px); }}

        .stat-card.mrsa {{ background: linear-gradient(135deg, #ef4444 0%, #991b1b 100%); }}
        .stat-card.mssa {{ background: linear-gradient(135deg, #3b82f6 0%, #1e40af 100%); }}
        .stat-card.unknown {{ background: linear-gradient(135deg, #9ca3af 0%, #4b5563 100%); }}

        .stat-value {{ font-size: 24px; font-weight: bold; margin-bottom: 5px; }}
        .stat-label {{ font-size: 12px; opacity: 0.9; text-transform: uppercase; letter-spacing: 0.5px; }}

        .controls-bar {{
            display: flex;
            flex-wrap: wrap;
            gap: 12px;
            align-items: center;
            margin-bottom: 15px;
            padding: 15px;
            background: #f8fafc;
            border-radius: 8px;
            border: 1px solid #e2e8f0;
        }}

        .search-wrapper {{ position: relative; flex: 1; min-width: 240px; }}
        .search-wrapper input {{
            width: 100%;
            padding: 10px 14px 10px 38px;
            border: 2px solid #cbd5e1;
            border-radius: 8px;
            font-size: 14px;
            outline: none;
            transition: border-color 0.2s;
        }}
        .search-wrapper input:focus {{ border-color: #3b82f6; box-shadow: 0 0 0 3px rgba(59,130,246,0.15); }}
        .search-wrapper::before {{
            content: '🔍';
            position: absolute;
            left: 12px;
            top: 50%;
            transform: translateY(-50%);
            font-size: 14px;
        }}

        .filter-pills {{ display: flex; gap: 8px; flex-wrap: wrap; }}
        .filter-pill {{
            padding: 8px 16px;
            border-radius: 20px;
            border: 2px solid #cbd5e1;
            background: white;
            color: #475569;
            cursor: pointer;
            font-size: 13px;
            font-weight: 600;
            transition: all 0.2s;
        }}
        .filter-pill:hover {{ background: #f1f5f9; }}
        .filter-pill.active {{ background: #3b82f6; color: white; border-color: #3b82f6; }}
        .filter-pill.active.mrsa-filter {{ background: #dc2626; border-color: #dc2626; }}
        .filter-pill.active.mssa-filter {{ background: #2563eb; border-color: #2563eb; }}

        .results-counter {{
            font-size: 13px;
            color: #64748b;
            font-weight: 600;
            white-space: nowrap;
        }}

        .table-container {{
            overflow-x: auto;
            overflow-y: auto;
            max-height: 600px;
            border-radius: 8px;
            box-shadow: 0 2px 10px rgba(0,0,0,0.05);
        }}

        table {{ width: 100%; border-collapse: collapse; font-size: 13px; }}

        th {{
            background: linear-gradient(135deg, #3b82f6 0%, #1e40af 100%);
            color: white;
            padding: 12px 10px;
            text-align: left;
            position: sticky;
            top: 0;
            z-index: 2;
            cursor: pointer;
            user-select: none;
            white-space: nowrap;
        }}

        th::after {{
            content: '⇅';
            margin-left: 6px;
            opacity: 0.45;
            font-size: 11px;
        }}

        th.sorted-asc::after {{ content: '▲'; opacity: 1; }}
        th.sorted-desc::after {{ content: '▼'; opacity: 1; }}

        td {{ padding: 10px; border-bottom: 1px solid #e5e7eb; }}
        tr:nth-child(even) td {{ background-color: #f8fafc; }}
        tr:hover td {{ background-color: #e0f2fe; }}

        .col-sample {{ white-space: nowrap; }}

        .mrsa-positive {{ color: #dc2626; font-weight: bold; background: rgba(220,38,38,0.08); }}
        .mrsa-negative {{ color: #2563eb; font-weight: bold; background: rgba(37,99,235,0.08); }}
        .unknown-status {{ color: #6b7280; font-weight: bold; }}

        .legend {{
            display: flex;
            gap: 15px;
            align-items: center;
            flex-wrap: wrap;
            margin-top: 12px;
            font-size: 12px;
            color: #475569;
        }}

        .legend-item {{ display: flex; align-items: center; gap: 6px; }}
        .legend-swatch {{ width: 14px; height: 14px; border-radius: 3px; display: inline-block; }}

        .dist-grid {{
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(340px, 1fr));
            gap: 18px;
            margin-top: 10px;
        }}

        .dist-card {{
            background: white;
            border-radius: 10px;
            overflow: hidden;
            box-shadow: 0 4px 15px rgba(0,0,0,0.08);
            border: 1px solid #e2e8f0;
            transition: transform 0.2s, box-shadow 0.2s;
        }}

        .dist-card:hover {{ transform: translateY(-3px); box-shadow: 0 8px 25px rgba(0,0,0,0.12); }}

        .dist-header {{
            display: flex;
            justify-content: space-between;
            align-items: center;
            padding: 14px 18px;
            color: white;
        }}

        .dist-title {{ font-size: 15px; font-weight: 700; letter-spacing: 0.3px; }}

        .dist-unique {{
            display: flex;
            align-items: baseline;
            gap: 4px;
            background: rgba(255,255,255,0.22);
            padding: 4px 10px;
            border-radius: 20px;
            font-size: 11px;
        }}

        .dist-unique-num {{ font-size: 16px; font-weight: 800; }}
        .dist-unique-lbl {{ opacity: 0.9; }}

        .dist-body {{ padding: 15px 18px; }}

        .dist-row {{ margin-bottom: 10px; }}
        .dist-label {{
            font-size: 12px;
            font-weight: 600;
            color: #334155;
            margin-bottom: 4px;
            white-space: nowrap;
            overflow: hidden;
            text-overflow: ellipsis;
        }}

        .dist-bar-track {{ background: #e2e8f0; height: 22px; border-radius: 6px; overflow: hidden; }}
        .dist-bar-fill {{
            height: 100%;
            border-radius: 6px;
            transition: width 0.7s ease;
            display: flex;
            align-items: center;
            justify-content: flex-end;
            padding-right: 8px;
            min-width: 44px;
        }}
        .dist-bar-count {{ font-size: 11px; font-weight: 700; color: white; text-shadow: 0 1px 2px rgba(0,0,0,0.35); }}

        .dist-empty {{ color: #94a3b8; font-style: italic; font-size: 13px; text-align: center; padding: 12px 0; }}

        .dist-footnote {{
            margin-top: 10px;
            font-size: 11px;
            color: #b45309;
            background: #fef3c7;
            padding: 6px 10px;
            border-radius: 6px;
            border-left: 3px solid #f59e0b;
        }}

        .footer {{
            text-align: center;
            margin-top: 30px;
            padding: 20px;
            background: rgba(0, 0, 0, 0.3);
            border-radius: 10px;
            font-size: 14px;
            color: white;
        }}

        .timestamp {{ color: #fbbf24; font-weight: bold; }}

        @media (max-width: 768px) {{
            .ascii-art {{ font-size: 6px; }}
            table {{ font-size: 12px; }}
            .stats-grid {{ grid-template-columns: repeat(2, 1fr); }}
            .controls-bar {{ flex-direction: column; align-items: stretch; }}
        }}
    </style>
</head>
<body>
    <div class="container">
        <div class="header">
            <div class="ascii-container">
                <div class="ascii-art">███████╗████████╗ █████╗ ██████╗ ██╗  ██╗███████╗ ██████╗  ██████╗ ██████╗ ███████╗
██╔════╝╚══██╔══╝██╔══██╗██╔══██╗██║  ██║██╔════╝██╔════╝ ██╔═══██╗██╔══██╗██╔════╝
███████╗   ██║   ███████║██████╔╝███████║███████╗██║      ██║   ██║██████╔╝█████╗  
╚════██║   ██║   ██╔══██║██╔═══╝ ██╔══██║╚════██║██║      ██║   ██║██╔═══╝ ██╔══╝  
███████║   ██║   ██║  ██║██║     ██║  ██║███████║╚██████╗ ╚██████╔╝██║     ███████╗
╚══════╝   ╚═╝   ╚═╝  ╚═╝╚═╝     ╚═╝  ╚═╝╚══════╝ ╚═════╝  ╚═════╝ ╚═╝     ╚══════╝</div>
            </div>

            <div class="quote-container" id="quoteContainer">
                <div class="quote-text" id="quoteText">"{random_quote['text']}"</div>
                <div class="quote-author" id="quoteAuthor">— {random_quote['author']}</div>
            </div>
        </div>

        {missing_banner}

        <div class="report-section">
            <h2>📊 Comprehensive Typing Report</h2>

            <div class="stats-grid">
                <div class="stat-card">
                    <div class="stat-value">{total_samples}</div>
                    <div class="stat-label">Total Samples</div>
                </div>
                <div class="stat-card mrsa">
                    <div class="stat-value">{mrsa_count}</div>
                    <div class="stat-label">MRSA</div>
                </div>
                <div class="stat-card mssa">
                    <div class="stat-value">{mssa_count}</div>
                    <div class="stat-label">MSSA</div>
                </div>
                <div class="stat-card unknown">
                    <div class="stat-value">{unknown_count}</div>
                    <div class="stat-label">Unknown</div>
                </div>
                <div class="stat-card">
                    <div class="stat-value">{unique_agr}</div>
                    <div class="stat-label">Unique agr</div>
                </div>
                <div class="stat-card">
                    <div class="stat-value">{unique_cap}</div>
                    <div class="stat-label">Unique Capsules</div>
                </div>
                <div class="stat-card">
                    <div class="stat-value">{unique_scc_rpet}</div>
                    <div class="stat-label">Unique SCCmec (RPet)</div>
                </div>
                <div class="stat-card">
                    <div class="stat-value">{unique_scc_sub}</div>
                    <div class="stat-label">Unique Subtypes</div>
                </div>
            </div>

            <p><strong>Generated:</strong> {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}</p>
            <p><strong>Report includes:</strong> MLST, spa, agr, capsule, SCCmec (CGE + RPet) types and subtypes, MRSA/MSSA classification</p>
        </div>

        <div class="report-section">
            <h2>🧬 Sample Typing Data ({total_samples} samples)</h2>

            <div class="controls-bar">
                <div class="search-wrapper">
                    <input type="text" id="tableSearch" placeholder="Search by sample, MLST, spa, agr, capsule, SCCmec...">
                </div>
                <div class="filter-pills">
                    <button class="filter-pill active" data-filter="all">All</button>
                    <button class="filter-pill mrsa-filter" data-filter="MRSA">MRSA</button>
                    <button class="filter-pill mssa-filter" data-filter="MSSA">MSSA</button>
                    <button class="filter-pill" data-filter="Unknown">Unknown</button>
                </div>
                <span class="results-counter" id="resultsCounter">{total_samples} shown</span>
            </div>

            <div class="table-container">
                <table id="mainTable">
                    <thead>
                        <tr>
                            <th>Sample</th>
                            <th>MLST</th>
                            <th>spa Type</th>
                            <th>agr Type</th>
                            <th>Capsule Type</th>
                            <th>SCCmec Type (CGE)</th>
                            <th>SCCmec Type (RPet)</th>
                            <th>SCCmec Subtype</th>
                            <th>MRSA/MSSA Status</th>
                        </tr>
                    </thead>
                    <tbody>
                        {table_rows}
                    </tbody>
                </table>
            </div>

            <div class="legend">
                <div class="legend-item"><span class="legend-swatch" style="background:#dc2626;"></span>MRSA</div>
                <div class="legend-item"><span class="legend-swatch" style="background:#2563eb;"></span>MSSA</div>
                <div class="legend-item"><span class="legend-swatch" style="background:#6b7280;"></span>Unknown</div>
                <div class="legend-item" style="margin-left:auto;">💡 Click any column header to sort</div>
            </div>
        </div>

        <div class="report-section">
            <h2>📊 Type Distribution &amp; Frequency Analysis</h2>
            <p style="color:#64748b; margin-bottom: 15px; font-size: 13px;">
                Each card shows the number of unique values and frequency (as % of typed samples) for a given field.
            </p>
            {dist_cards_html}
        </div>

        <div class="footer">
            <p><strong>STAPHSCOPE</strong> - Comprehensive Analysis Report</p>
            <p class="timestamp">Generated: {datetime.now().strftime("%Y-%m-%d %H:%M:%S")}</p>
            <p><strong>Author:</strong> Brown Beckley | <strong>GitHub:</strong> bbeckley-hub</p>
            <p><strong>Email:</strong> brownbeckley94@gmail.com</p>
            <p><strong>Affiliation:</strong> University of Ghana Medical School - Department of Medical Biochemistry</p>
            <p><em>Report generated from MLST, spa, agr, capsule, and SCCmec (CGE + RPet) typing results.</em></p>
            <p><em>⭐ Star us on GitHub if you find this tool useful!</em></p>
            <p><em>Transforming fragmented genomic data into coherent biological narratives 🧬✨</em></p>
        </div>
    </div>

    <script>
        const quotes = {json.dumps(SCIENCE_QUOTES)};

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

        setInterval(displayQuote, 10000);

        const table = document.getElementById('mainTable');
        const tbody = table.querySelector('tbody');
        const allRows = Array.from(tbody.querySelectorAll('tr'));
        const searchInput = document.getElementById('tableSearch');
        const resultsCounter = document.getElementById('resultsCounter');
        const filterPills = document.querySelectorAll('.filter-pill');

        let activeFilter = 'all';

        function applyFilters() {{
            const term = (searchInput.value || '').toLowerCase().trim();
            let visible = 0;
            allRows.forEach(row => {{
                const text = row.textContent.toLowerCase();
                const mrsa = row.getAttribute('data-mrsa');
                const matchesText = !term || text.includes(term);
                const matchesFilter = activeFilter === 'all' || mrsa === activeFilter;
                if (matchesText && matchesFilter) {{
                    row.style.display = '';
                    visible++;
                }} else {{
                    row.style.display = 'none';
                }}
            }});
            resultsCounter.textContent = visible + ' shown';
        }}

        searchInput.addEventListener('input', applyFilters);

        filterPills.forEach(pill => {{
            pill.addEventListener('click', () => {{
                filterPills.forEach(p => p.classList.remove('active'));
                pill.classList.add('active');
                activeFilter = pill.getAttribute('data-filter');
                applyFilters();
            }});
        }});

        const headers = table.querySelectorAll('th');
        headers.forEach((th, idx) => {{
            th.addEventListener('click', () => {{
                const isAsc = th.classList.contains('sorted-asc');
                headers.forEach(h => h.classList.remove('sorted-asc', 'sorted-desc'));

                const sorted = allRows.slice().sort((a, b) => {{
                    const aText = a.children[idx].textContent.trim();
                    const bText = b.children[idx].textContent.trim();
                    const aNum = parseFloat(aText.replace(/[^0-9.\\-]/g, ''));
                    const bNum = parseFloat(bText.replace(/[^0-9.\\-]/g, ''));
                    const bothNum = !isNaN(aNum) && !isNaN(bNum);
                    if (bothNum) return isAsc ? bNum - aNum : aNum - bNum;
                    return isAsc ? bText.localeCompare(aText) : aText.localeCompare(bText);
                }});

                sorted.forEach(r => tbody.appendChild(r));
                th.classList.add(isAsc ? 'sorted-desc' : 'sorted-asc');
            }});
        }});
    </script>
</body>
</html>'''

        with open(output_file, 'w', encoding='utf-8') as f:
            f.write(html_content)

        print("   ✓ HTML report generated successfully")
        return True
    except Exception as e:
        print(f"   ✗ Error generating HTML: {e}")
        return False


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    parser = argparse.ArgumentParser(
        description='StaphScope Comprehensive Report Generator',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
All input files are OPTIONAL. The report is generated from whatever is
available in the target directory. Missing files simply produce 'Not Assigned'
values in the corresponding columns.

Recognized input files (in the target directory):
  • mlst_summary.tsv
  • spa_summary.tsv
  • staphscope_sccmec_cge_summary.tsv
  • staphscope_sccmec_rpet_summary.tsv
  • staphscope_capsule_summary.tsv
  • agr_summary.tsv
        """
    )
    parser.add_argument('--dir', '-d', default=None,
                        help='Directory containing the summary TSVs (default: current directory)')
    args = parser.parse_args()

    print("\n" + "="*60)
    print("STAPHSCOPE COMPREHENSIVE REPORT GENERATOR")
    print("="*60)

    target_dir = Path(args.dir).resolve() if args.dir else Path.cwd()
    if not target_dir.exists() or not target_dir.is_dir():
        print(f"❌ Error: directory does not exist: {target_dir}")
        sys.exit(1)

    print(f"📂 Working directory: {target_dir}")

    # Map label → existing path (or None)
    resolved = {key: (target_dir / fname) for key, fname in INPUT_FILES.items()}
    present = {k: p for k, p in resolved.items() if p.exists()}
    missing = {k: p for k, p in resolved.items() if not p.exists()}

    print("\n🔎 Input file check:")
    for key, fname in INPUT_FILES.items():
        mark = "✓" if key in present else "•"
        status = "found" if key in present else "missing (column will be 'Not Assigned')"
        print(f"   {mark} {fname:<45s} — {status}")

    if not present:
        print("\n❌ No recognized input files were found — nothing to report.")
        print("   Expected at least one of:")
        for fname in INPUT_FILES.values():
            print(f"     • {fname}")
        sys.exit(1)

    # Load whatever exists
    mlst_data       = load_mlst_data(present['mlst'])               if 'mlst' in present else {}
    spa_data        = load_spa_data(present['spa'])                 if 'spa' in present else {}
    sccmec_cge_data = load_sccmec_cge_data(present['sccmec_cge'])   if 'sccmec_cge' in present else {}
    sccmec_rpet_data= load_sccmec_rpet_data(present['sccmec_rpet']) if 'sccmec_rpet' in present else {}
    capsule_data    = load_capsule_data(present['capsule'])         if 'capsule' in present else {}
    agr_data        = load_agr_data(present['agr'])                 if 'agr' in present else {}

    combined_data = combine_data(
        mlst_data, spa_data, sccmec_cge_data, sccmec_rpet_data,
        capsule_data, agr_data
    )

    if not combined_data:
        print("\n❌ No samples were found in any of the available files.")
        sys.exit(1)

    distributions = compute_distributions(combined_data)

    # Prepare "missing sources" list for report metadata / banner
    missing_sources = [(key, INPUT_FILES[key]) for key in missing]

    base_name = "staphscope_comprehensive_report"
    tsv_success  = generate_tsv_report(combined_data, target_dir / f"{base_name}.tsv")
    json_success = generate_json_report(combined_data, distributions, missing_sources,
                                        target_dir / f"{base_name}.json")
    html_success = generate_html_report(combined_data, distributions, missing_sources,
                                        target_dir / f"{base_name}.html")

    print("\n" + "="*60)
    if tsv_success and json_success and html_success:
        print("✅ COMPREHENSIVE REPORT GENERATED SUCCESSFULLY!")
        print("="*60)
        print(f"📁 Reports saved in: {target_dir}")
        print(f"   • {base_name}.html (interactive HTML report)")
        print(f"   • {base_name}.json (JSON data + distributions)")
        print(f"   • {base_name}.tsv  (TSV data)")
        print(f"📊 Summary: {len(combined_data)} samples")

        mrsa_count = sum(1 for d in combined_data if d['mrsa_status'] == 'MRSA')
        mssa_count = sum(1 for d in combined_data if d['mrsa_status'] == 'MSSA')
        print(f"🦠 MRSA: {mrsa_count}, MSSA: {mssa_count}")
        print(f"🧬 Unique — agr: {distributions['agr_type']['unique_count']}, "
              f"capsule: {distributions['capsule_type']['unique_count']}, "
              f"SCCmec(RPet): {distributions['sccmec_type_rpet']['unique_count']}, "
              f"subtypes: {distributions['sccmec_subtype']['unique_count']}")
        if missing_sources:
            print("⚠️  Missing sources (columns shown as 'Not Assigned'):")
            for key, fname in missing_sources:
                print(f"     • {fname}")
    else:
        print("❌ Some reports failed to generate")

    print("="*60)


if __name__ == "__main__":
    main()
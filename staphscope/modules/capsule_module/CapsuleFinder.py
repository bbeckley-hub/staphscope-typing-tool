#!/usr/bin/env python3
"""
CapsuleFinder - Capsule Typing Module for StaphScope
Original Module/Logic Author: Riccardo Bollini (developer of StaphScan)
Integrated and extended by: Beckley Brown <brownbeckley94@gmail.com>
Date: 2026-09-17
Affiliation: University of Ghana Medical School - Department of Medical Biochemistry

Reverse BLAST against a capsule targets database to type S. aureus capsule
locus 5 or 8, assess operon completeness, and generate per-sample deep-dive
HTML/JSON reports plus top-level batch summaries.

Usage:
    python3 CapsuleFinder.py -i *.fna -d capsule_results -db_dir database
"""

import argparse
import io
import json
import re
import subprocess
import sys
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional, Set

import pandas as pd


# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

NOT_ASSIGNED = "Not Assigned"

OPERON_GENES = ["A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N", "O", "P"]

TYPE_LABELS = {
    "cap5": "Type 5",
    "cap8": "Type 8",
}


# ---------------------------------------------------------------------------
# Data containers
# ---------------------------------------------------------------------------

@dataclass
class BlastHit:
    """A single reverse-BLAST alignment record (database gene vs. assembly contig)."""
    qseqid: str
    sseqid: str
    pident: float
    length: int
    qlen: int
    evalue: float = 0.0
    bitscore: float = 0.0

    @property
    def coverage(self) -> float:
        """Percent of the database gene covered by this alignment."""
        if self.qlen <= 0:
            return 0.0
        return (self.length / self.qlen) * 100.0


@dataclass
class SampleResult:
    """All capsule typing results and supporting evidence for one sample."""
    sample: str
    input_file: str
    cap_type: str = NOT_ASSIGNED
    cap_completeness: str = NOT_ASSIGNED
    cap_genes: str = NOT_ASSIGNED
    cap5_score: int = 0
    cap8_score: int = 0
    found_genes: Set[str] = field(default_factory=set)
    gene_hits: Dict[str, BlastHit] = field(default_factory=dict)
    raw_hits: List[BlastHit] = field(default_factory=list)
    evidence: List[str] = field(default_factory=list)
    warnings: List[str] = field(default_factory=list)


# ---------------------------------------------------------------------------
# CapsuleTyper
# ---------------------------------------------------------------------------

class CapsuleTyper:
    """Reverse-BLAST capsule typing using a targets database of capsule locus genes."""

    def __init__(self, database_dir: Path, min_id: float = 90.0, min_cov: float = 80.0):
        """Initialise with a database directory containing a capsule targets FASTA file."""
        self.database_dir = Path(database_dir)
        self.min_id = min_id
        self.min_cov = min_cov
        self.db_fasta = self._locate_db()

    def _locate_db(self) -> Path:
        """Return the latest capsule targets FASTA file in the database directory."""
        matches = sorted(self.database_dir.glob("capsule_targets*.fasta"))
        if not matches:
            raise FileNotFoundError(
                f"Missing database: no file matching 'capsule_targets*.fasta' in {self.database_dir}"
            )
        return matches[-1]

    def check_db(self) -> bool:
        """Return True if the capsule database file exists."""
        return self.db_fasta.exists()

    def run(self, assembly_path: Path) -> SampleResult:
        """Run reverse BLAST and capsule typing on a single assembly."""
        sample_name = assembly_path.stem
        raw_hits = self._blast(assembly_path)

        result = SampleResult(
            sample=sample_name,
            input_file=str(assembly_path),
            raw_hits=raw_hits,
        )

        if not raw_hits:
            result.warnings.append("No capsule genes detected above thresholds.")
            result.evidence.append("No capsule locus genes found")
            return result

        best_per_gene: Dict[str, BlastHit] = {}
        for hit in raw_hits:
            gene = self._extract_gene_base(hit.qseqid)
            if gene not in best_per_gene or hit.coverage > best_per_gene[gene].coverage:
                best_per_gene[gene] = hit

        found_genes = set(best_per_gene.keys())
        result.gene_hits = best_per_gene
        result.found_genes = found_genes
        result.cap_genes = ";".join(sorted(found_genes))

        score_5 = sum(1 for g in found_genes if g.startswith("cap5"))
        score_8 = sum(1 for g in found_genes if g.startswith("cap8"))
        result.cap5_score = score_5
        result.cap8_score = score_8

        detected = self._call_type(score_5, score_8)
        if detected is None:
            result.warnings.append(
                f"Ambiguous capsule type: cap5 score = {score_5}, cap8 score = {score_8}"
            )
            result.evidence.append(f"Cap5 genes: {score_5}, Cap8 genes: {score_8}")
            return result

        result.cap_type = TYPE_LABELS[detected]
        result.cap_completeness = self._assess_completeness(detected, found_genes)

        if detected == "cap5":
            result.evidence.append(f"cap5 locus detected ({score_5} of 16 operon genes)")
        else:
            result.evidence.append(f"cap8 locus detected ({score_8} of 16 operon genes)")

        if score_5 > 0 and score_8 > 0:
            result.warnings.append(
                "Both cap5 and cap8 genes detected; possible contamination or mixed assembly."
            )

        return result

    @staticmethod
    def _extract_gene_base(qseqid: str) -> str:
        """
        Extract the canonical gene name from a database sequence ID.
        Matches cap5A / cap8K / cap5_A / cap5-A / cap5A1 / CAP8B and normalises
        to the form cap5[A-P] / cap8[A-P].
        """
        match = re.search(r"cap([58])[_\-]?([A-P])", str(qseqid), re.IGNORECASE)
        if match:
            return f"cap{match.group(1)}{match.group(2).upper()}"
        return str(qseqid)

    @staticmethod
    def _call_type(score_5: int, score_8: int) -> Optional[str]:
        """
        Decide cap5 vs cap8 from the number of distinct operon genes detected.
        Returns 'cap5', 'cap8', or None if the result is ambiguous.
        """
        if score_5 == 0 and score_8 == 0:
            return None
        if score_5 > score_8:
            return "cap5"
        if score_8 > score_5:
            return "cap8"
        return None

    @staticmethod
    def _assess_completeness(detected: str, found_genes: Set[str]) -> str:
        """Assess operon completeness for the detected capsule type."""
        expected = {f"{detected}{locus}" for locus in OPERON_GENES}
        found_count = len(expected & found_genes)
        total = len(expected)

        if found_count == total:
            return "Complete"
        if found_count >= total // 2:
            return "Partial"
        return "Incomplete"

    def _blast(self, assembly_path: Path) -> List[BlastHit]:
        """
        Run reverse blastn (database as query, assembly as subject) and return
        all hits passing identity and coverage thresholds.
        """
        cmd = [
            "blastn",
            "-task", "blastn",
            "-query", str(self.db_fasta),
            "-subject", str(assembly_path),
            "-outfmt", "6 qseqid sseqid pident length qlen evalue bitscore",
            "-max_target_seqs", "500",
        ]
        try:
            res = subprocess.run(cmd, capture_output=True, text=True)
            if not res.stdout or not res.stdout.strip():
                return []

            df = pd.read_csv(
                io.StringIO(res.stdout),
                sep="\t",
                names=["qseqid", "sseqid", "pident", "length", "qlen", "evalue", "bitscore"],
            )
            if df.empty:
                return []

            df["cov"] = (df["length"] / df["qlen"]) * 100
            df = df[(df["pident"] >= self.min_id) & (df["cov"] >= self.min_cov)]
            if df.empty:
                return []

            hits = [
                BlastHit(
                    qseqid=row["qseqid"],
                    sseqid=row["sseqid"],
                    pident=float(row["pident"]),
                    length=int(row["length"]),
                    qlen=int(row["qlen"]),
                    evalue=float(row["evalue"]),
                    bitscore=float(row["bitscore"]),
                )
                for _, row in df.iterrows()
            ]
            hits.sort(key=lambda h: (h.coverage, h.pident), reverse=True)
            return hits

        except Exception:
            return []


# ---------------------------------------------------------------------------
# Per-sample report writers
# ---------------------------------------------------------------------------

def write_sample_json(result: SampleResult, output_dir: Path) -> Path:
    """Write a per-sample JSON report with full evidence details."""
    output_dir.mkdir(parents=True, exist_ok=True)
    payload = {
        "sample": result.sample,
        "input_file": result.input_file,
        "cap_type": result.cap_type,
        "cap_completeness": result.cap_completeness,
        "cap_genes": result.cap_genes,
        "cap5_score": result.cap5_score,
        "cap8_score": result.cap8_score,
        "evidence": result.evidence,
        "warnings": result.warnings,
        "gene_hits": [
            {
                "gene": gene,
                "contig": hit.sseqid,
                "pident": round(hit.pident, 2),
                "length": hit.length,
                "gene_length": hit.qlen,
                "coverage": round(hit.coverage, 2),
                "bitscore": hit.bitscore,
            }
            for gene, hit in sorted(result.gene_hits.items())
        ],
        "raw_hits": [
            {
                "gene": h.qseqid,
                "contig": h.sseqid,
                "pident": round(h.pident, 2),
                "length": h.length,
                "gene_length": h.qlen,
                "coverage": round(h.coverage, 2),
                "bitscore": h.bitscore,
            }
            for h in result.raw_hits
        ],
    }
    path = output_dir / f"{result.sample}_capsule.json"
    with open(path, "w") as fh:
        json.dump(payload, fh, indent=2)
    return path


def _hit_rows_html(hits: List[BlastHit]) -> str:
    """Render BLAST hit rows as HTML table rows with coverage colouring."""
    if not hits:
        return "<tr><td colspan='7' style='text-align:center;color:#9ca3af;'>No hits above thresholds</td></tr>"
    rows = []
    for i, h in enumerate(hits):
        if i == 0:
            cls = "highlight-red"
        elif h.coverage >= 90.0:
            cls = "highlight-orange"
        elif h.coverage >= 70.0:
            cls = "highlight-yellow"
        else:
            cls = ""
        rows.append(
            f"<tr class='{cls}'>"
            f"<td>{i + 1}</td>"
            f"<td><strong>{h.qseqid}</strong></td>"
            f"<td>{h.sseqid}</td>"
            f"<td>{h.pident:.2f}</td>"
            f"<td>{h.length}</td>"
            f"<td>{h.qlen}</td>"
            f"<td>{h.coverage:.2f}</td>"
            f"</tr>"
        )
    return "\n".join(rows)


def write_sample_html(result: SampleResult, output_dir: Path) -> Path:
    """Write a per-sample deep-dive HTML report with all BLAST hits and evidence."""
    output_dir.mkdir(parents=True, exist_ok=True)

    type_slug = result.cap_type.lower().replace(" ", "-")
    badge_class = {
        "type-5": "status-type5",
        "type-8": "status-type8",
    }.get(type_slug, "status-unknown")

    evidence_items = "".join(f"<li>{e}</li>" for e in result.evidence)
    if not evidence_items:
        evidence_items = "<li>No evidence recorded</li>"

    warning_html = ""
    if result.warnings:
        warning_html = "<div class='warning-box'><strong>⚠️ Warnings:</strong><ul>" + \
            "".join(f"<li>{w}</li>" for w in result.warnings) + "</ul></div>"

    gene_cards = "".join(
        f"<div class='gene-card'>{g}</div>" for g in sorted(result.found_genes)
    )
    if not gene_cards:
        gene_cards = "<p style='color:#6b7280;'>No capsule genes detected above thresholds.</p>"

    summary_hits = list(result.gene_hits.values())
    summary_hits.sort(key=lambda h: (h.coverage, h.pident), reverse=True)
    summary_rows = _hit_rows_html(summary_hits)
    raw_rows = _hit_rows_html(result.raw_hits)

    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>STAPHSCOPE - {result.sample} Capsule Report</title>
<style>
* {{ margin:0; padding:0; box-sizing:border-box; }}
body {{
    background: linear-gradient(135deg, #1e3c72 0%, #2a5298 50%, #7e22ce 100%);
    font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif;
    color:#ffffff; padding:20px; min-height:100vh;
}}
.container {{ max-width:1400px; margin:0 auto; }}
.ascii-container {{
    background: rgba(0,0,0,0.7); padding:20px; border-radius:15px;
    margin-bottom:20px; box-shadow:0 8px 32px rgba(0,0,0,0.4);
    border:2px solid rgba(0,255,0,0.3);
}}
.ascii-art {{
    font-family:'Courier New', monospace; font-size:10px; line-height:1.1;
    white-space:pre; color:#00ff00; text-shadow:0 0 10px rgba(0,255,0,0.5);
    text-align:center; overflow-x:auto;
}}
.quote-container {{
    background: rgba(255,255,255,0.1); backdrop-filter:blur(10px);
    padding:20px; border-radius:10px; margin-bottom:30px; text-align:center;
    min-height:100px; display:flex; flex-direction:column;
    justify-content:center; box-shadow:0 4px 20px rgba(0,0,0,0.3);
    border:1px solid rgba(255,255,255,0.2); transition:opacity 0.5s ease-in-out;
}}
.quote-text {{ font-size:18px; font-style:italic; margin-bottom:10px; color:#ffffff; }}
.quote-author {{ font-size:14px; color:#fbbf24; font-weight:bold; }}
.report-section {{
    background: rgba(255,255,255,0.95); color:#1f2937;
    padding:25px; border-radius:10px; margin-bottom:20px;
    box-shadow:0 4px 15px rgba(0,0,0,0.2);
}}
.report-section h2 {{
    color:#1e3a8a; border-bottom:3px solid #3b82f6;
    padding-bottom:10px; margin-bottom:20px; font-size:24px;
}}
.report-section h3 {{
    color:#1e40af; margin-top:20px; margin-bottom:10px; font-size:18px;
}}
.metrics-grid {{
    display:grid; grid-template-columns:repeat(auto-fit, minmax(200px, 1fr));
    gap:20px; margin-top:15px;
}}
.metric-card {{
    background: linear-gradient(135deg, #8b5cf6 0%, #6d28d9 100%);
    color:white; padding:20px; border-radius:8px;
    box-shadow:0 4px 12px rgba(0,0,0,0.15); text-align:center;
}}
.metric-label {{
    font-size:13px; opacity:0.9; margin-bottom:5px;
    text-transform:uppercase; letter-spacing:0.5px;
}}
.metric-value {{ font-size:24px; font-weight:bold; }}
.status-badge {{
    display:inline-block; padding:8px 16px; border-radius:20px;
    font-weight:bold; font-size:16px; margin:10px 0;
}}
.status-type5 {{ background:#2563eb; color:white; }}
.status-type8 {{ background:#dc2626; color:white; }}
.status-unknown {{ background:#6b7280; color:white; }}
.detailed-table {{
    width:100%; border-collapse:collapse; margin-top:15px; font-size:13px;
}}
.detailed-table th {{
    background:#1e40af; color:white; padding:12px; text-align:left;
    font-weight:bold; position:sticky; top:0;
}}
.detailed-table td {{
    padding:10px; border-bottom:1px solid #e5e7eb; vertical-align:top;
}}
.detailed-table tr:hover {{ background:#f3f4f6; }}
.highlight-red {{ background-color:#fee2e2 !important; color:#dc2626 !important; font-weight:bold; }}
.highlight-orange {{ background-color:#ffedd5 !important; color:#ea580c !important; }}
.highlight-yellow {{ background-color:#fef3c7 !important; color:#d97706 !important; }}
.gene-grid {{
    display:grid; grid-template-columns:repeat(auto-fill, minmax(120px, 1fr));
    gap:10px; margin-top:10px;
}}
.gene-card {{
    background: linear-gradient(135deg, #3b82f6 0%, #1e40af 100%);
    color:white; padding:12px; border-radius:8px;
    text-align:center; font-weight:bold; font-size:13px;
}}
.warning-box {{
    background:#fef3c7; border-left:4px solid #f59e0b;
    padding:15px; margin:15px 0; border-radius:6px; color:#92400e;
}}
.footer {{
    text-align:center; margin-top:30px; padding:20px;
    background:rgba(0,0,0,0.3); border-radius:10px; font-size:14px;
}}
.timestamp {{ color:#fbbf24; font-weight:bold; }}
.authorship {{
    margin-top:15px; padding:15px; background:rgba(255,255,255,0.1);
    border-radius:8px; font-size:12px;
}}
@media (max-width:768px) {{
    .ascii-art {{ font-size:6px; }}
    .detailed-table {{ font-size:11px; }}
}}
</style>
</head>
<body>
<div class="container">
<div class="ascii-container">
<div class="ascii-art">███████╗████████╗ █████╗ ██████╗ ██╗  ██╗███████╗ ██████╗ ██████╗ ██████╗ ███████╗
██╔════╝╚══██╔══╝██╔══██╗██╔══██╗██║  ██║██╔════╝██╔════╝██╔═══██╗██╔══██╗██╔════╝
███████╗   ██║   ███████║██████╔╝███████║███████╗██║     ██║   ██║██████╔╝█████╗  
╚════██║   ██║   ██╔══██║██╔═══╝ ██╔══██║╚════██║██║     ██║   ██║██╔═══╝ ██╔══╝  
███████║   ██║   ██║  ██║██║     ██║  ██║███████║╚██████╗╚██████╔╝██║     ███████╗
╚══════╝   ╚═╝   ╚═╝  ╚═╝╚═╝     ╚═╝  ╚═╝╚══════╝ ╚═════╝ ╚═════╝ ╚═╝     ╚══════╝</div>
</div>

<div class="quote-container" id="quoteContainer">
    <div class="quote-text" id="quoteText"></div>
    <div class="quote-author" id="quoteAuthor"></div>
</div>

<div class="report-section">
<h2>📊 Sample Information</h2>
<div class="metrics-grid">
<div class="metric-card"><div class="metric-label">Sample</div><div class="metric-value" style="font-size:16px;">{result.sample}</div></div>
<div class="metric-card"><div class="metric-label">Capsule Type</div><div class="metric-value" style="font-size:18px;">{result.cap_type}</div></div>
<div class="metric-card"><div class="metric-label">Completeness</div><div class="metric-value" style="font-size:18px;">{result.cap_completeness}</div></div>
<div class="metric-card"><div class="metric-label">Genes Found</div><div class="metric-value">{len(result.found_genes)}</div></div>
</div>
</div>

<div class="report-section">
<h2>🧬 Capsule Classification</h2>
<div class="status-badge {badge_class}">{result.cap_type}</div>
<h3>Evidence</h3>
<ul style="margin-left:20px; line-height:1.8;">{evidence_items}</ul>
{warning_html}
</div>

<div class="report-section">
<h2>📈 Operon Scores</h2>
<div class="metrics-grid">
<div class="metric-card"><div class="metric-label">Cap5 Genes</div><div class="metric-value">{result.cap5_score}</div></div>
<div class="metric-card"><div class="metric-label">Cap8 Genes</div><div class="metric-value">{result.cap8_score}</div></div>
</div>
</div>

<div class="report-section">
<h2>🧪 Detected Gene Summary</h2>
<div class="gene-grid">{gene_cards}</div>
</div>

<div class="report-section">
<h2>🎯 Best Hit per Gene</h2>
<div style="max-height:500px; overflow-y:auto;">
<table class="detailed-table">
<thead>
<tr><th>Rank</th><th>Gene</th><th>Contig</th><th>% Identity</th><th>Alignment Length</th><th>Gene Length</th><th>Coverage %</th></tr>
</thead>
<tbody>
{summary_rows}
</tbody>
</table>
</div>
<p style="margin-top:10px; font-size:12px; color:#666;">
● Top hit in RED &nbsp;|&nbsp; <span style="color:#ea580c;">● Coverage ≥90% in ORANGE</span> &nbsp;|&nbsp; <span style="color:#d97706;">● Coverage ≥70% in YELLOW</span>
</p>
</div>

<div class="report-section">
<h2>🔍 All Raw BLAST Hits</h2>
<div style="max-height:500px; overflow-y:auto;">
<table class="detailed-table">
<thead>
<tr><th>Rank</th><th>Gene</th><th>Contig</th><th>% Identity</th><th>Alignment Length</th><th>Gene Length</th><th>Coverage %</th></tr>
</thead>
<tbody>
{raw_rows}
</tbody>
</table>
</div>
</div>

<div class="footer">
<p><strong>STAPHSCOPE</strong> - Capsule Typing</p>
<p class="timestamp">Generated: {time.strftime('%Y-%m-%d %H:%M:%S')}</p>
<div class="authorship">
<p><strong>Technical Support &amp; Inquiries:</strong></p>
<p>Original module author: Riccardo Bollini (StaphScan)</p>
<p>Integrated by: Brown Beckley</p>
<p>GitHub: <a href="https://github.com/bbeckley-hub" style="color:#fbbf24;">bbeckley-hub</a></p>
<p>Email: <a href="mailto:brownbeckley94@gmail.com" style="color:#fbbf24;">brownbeckley94@gmail.com</a></p>
<p>Affiliation: University of Ghana Medical School - Department of Medical Biochemistry</p>
</div>
</div>
</div>

<script>
const quotes = [
    {{ text: "The important thing is not to stop questioning. Curiosity has its own reason for existing.", author: "Albert Einstein" }},
    {{ text: "Science is not only a disciple of reason but also one of romance and passion.", author: "Stephen Hawking" }},
    {{ text: "Somewhere, something incredible is waiting to be known.", author: "Carl Sagan" }},
    {{ text: "The good thing about science is that it's true whether or not you believe in it.", author: "Neil deGrasse Tyson" }},
    {{ text: "In science, there are no shortcuts to truth.", author: "Karl Popper" }},
    {{ text: "Science knows no country, because knowledge belongs to humanity.", author: "Louis Pasteur" }},
    {{ text: "The science of today is the technology of tomorrow.", author: "Edward Teller" }},
    {{ text: "Nothing in life is to be feared, it is only to be understood.", author: "Marie Curie" }},
    {{ text: "Research is what I'm doing when I don't know what I'm doing.", author: "Wernher von Braun" }},
    {{ text: "The universe is not required to be in perfect harmony with human ambition.", author: "Carl Sagan" }}
];

const quoteContainer = document.getElementById('quoteContainer');
const quoteText = document.getElementById('quoteText');
const quoteAuthor = document.getElementById('quoteAuthor');

function displayQuote() {{
    quoteContainer.style.opacity = '0';
    setTimeout(() => {{
        const q = quotes[Math.floor(Math.random() * quotes.length)];
        quoteText.textContent = '"' + q.text + '"';
        quoteAuthor.textContent = '— ' + q.author;
        quoteContainer.style.opacity = '1';
    }}, 500);
}}

displayQuote();
setInterval(displayQuote, 10000);
</script>
</body>
</html>
"""
    path = output_dir / f"{result.sample}_capsule.html"
    with open(path, "w", encoding="utf-8") as fh:
        fh.write(html)
    return path


# ---------------------------------------------------------------------------
# Top-level summary writers
# ---------------------------------------------------------------------------

def write_summary_tsv(results: List[SampleResult], output_path: Path) -> None:
    """Write a batch TSV summary with one row per sample."""
    columns = [
        "sample", "cap_type", "cap_completeness",
        "cap5_score", "cap8_score", "cap_genes",
    ]
    with open(output_path, "w") as fh:
        fh.write("\t".join(columns) + "\n")
        for r in results:
            row = [
                r.sample,
                r.cap_type,
                r.cap_completeness,
                r.cap5_score,
                r.cap8_score,
                r.cap_genes,
            ]
            fh.write("\t".join(str(v) for v in row) + "\n")


def write_summary_json(results: List[SampleResult], output_path: Path) -> None:
    """Write a combined JSON summary for all samples."""
    payload = []
    for r in results:
        payload.append({
            "sample": r.sample,
            "cap_type": r.cap_type,
            "cap_completeness": r.cap_completeness,
            "cap_genes": r.cap_genes,
            "cap5_score": r.cap5_score,
            "cap8_score": r.cap8_score,
            "evidence": r.evidence,
            "warnings": r.warnings,
        })
    with open(output_path, "w") as fh:
        json.dump(payload, fh, indent=2)


def write_summary_html(results: List[SampleResult], output_path: Path, elapsed: float) -> None:
    """Write a StaphScope-styled batch summary HTML report."""
    total = len(results)
    type5_count = sum(1 for r in results if r.cap_type == "Type 5")
    type8_count = sum(1 for r in results if r.cap_type == "Type 8")
    not_assigned_count = total - type5_count - type8_count

    rows = []
    for r in results:
        badge_class = {
            "Type 5": "status-type5",
            "Type 8": "status-type8",
        }.get(r.cap_type, "status-unknown")
        rows.append(
            f"<tr>"
            f"<td><strong>{r.sample}</strong></td>"
            f"<td><span class='{badge_class}'>{r.cap_type}</span></td>"
            f"<td>{r.cap_completeness}</td>"
            f"<td>{r.cap5_score}</td>"
            f"<td>{r.cap8_score}</td>"
            f"<td class='genes'>{r.cap_genes}</td>"
            f"</tr>"
        )
    table_rows = "\n".join(rows)

    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>STAPHSCOPE - Capsule Batch Summary</title>
<style>
* {{ margin:0; padding:0; box-sizing:border-box; }}
body {{
    background: linear-gradient(135deg, #1e3c72 0%, #2a5298 50%, #7e22ce 100%);
    font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif;
    color:#ffffff; padding:20px; min-height:100vh;
}}
.container {{ max-width:1600px; margin:0 auto; }}
.ascii-container {{
    background: rgba(0,0,0,0.7); padding:20px; border-radius:15px;
    margin-bottom:20px; box-shadow:0 8px 32px rgba(0,0,0,0.4);
    border:2px solid rgba(0,255,0,0.3);
}}
.ascii-art {{
    font-family:'Courier New', monospace; font-size:10px; line-height:1.1;
    white-space:pre; color:#00ff00; text-shadow:0 0 10px rgba(0,255,0,0.5);
    text-align:center; overflow-x:auto;
}}
.quote-container {{
    background: rgba(255,255,255,0.1); backdrop-filter:blur(10px);
    padding:20px; border-radius:10px; margin-bottom:30px; text-align:center;
    min-height:100px; display:flex; flex-direction:column;
    justify-content:center; box-shadow:0 4px 20px rgba(0,0,0,0.3);
    border:1px solid rgba(255,255,255,0.2); transition:opacity 0.5s ease-in-out;
}}
.quote-text {{ font-size:18px; font-style:italic; margin-bottom:10px; }}
.quote-author {{ font-size:14px; color:#fbbf24; font-weight:bold; }}
.report-section {{
    background: rgba(255,255,255,0.95); color:#1f2937;
    padding:25px; border-radius:10px; margin-bottom:20px;
    box-shadow:0 4px 15px rgba(0,0,0,0.2);
}}
.report-section h2 {{
    color:#1e3a8a; border-bottom:3px solid #3b82f6;
    padding-bottom:10px; margin-bottom:20px; font-size:24px;
}}
.metrics-grid {{
    display:grid; grid-template-columns:repeat(auto-fit, minmax(180px, 1fr));
    gap:20px; margin-top:15px;
}}
.metric-card {{
    background: linear-gradient(135deg, #8b5cf6 0%, #6d28d9 100%);
    color:white; padding:20px; border-radius:8px;
    box-shadow:0 4px 12px rgba(0,0,0,0.15); text-align:center;
}}
.metric-label {{
    font-size:12px; opacity:0.9; margin-bottom:5px;
    text-transform:uppercase; letter-spacing:0.5px;
}}
.metric-value {{ font-size:28px; font-weight:bold; }}
.results-table {{
    width:100%; border-collapse:collapse; margin-top:15px; font-size:13px;
}}
.results-table th {{
    background:#1e40af; color:white; padding:12px; text-align:left;
    font-weight:bold; position:sticky; top:0;
}}
.results-table td {{
    padding:10px; border-bottom:1px solid #e5e7eb; vertical-align:top;
}}
.results-table tr:hover {{ background:#f3f4f6; }}
.status-type5 {{
    background:#2563eb; color:white; padding:4px 8px;
    border-radius:4px; font-weight:bold; font-size:11px;
}}
.status-type8 {{
    background:#dc2626; color:white; padding:4px 8px;
    border-radius:4px; font-weight:bold; font-size:11px;
}}
.status-unknown {{
    background:#6b7280; color:white; padding:4px 8px;
    border-radius:4px; font-weight:bold; font-size:11px;
}}
.genes {{ font-family:'Courier New', monospace; font-size:11px; color:#4b5563; }}
.footer {{
    text-align:center; margin-top:30px; padding:20px;
    background:rgba(0,0,0,0.3); border-radius:10px; font-size:14px;
}}
.timestamp {{ color:#fbbf24; font-weight:bold; }}
.authorship {{
    margin-top:15px; padding:15px; background:rgba(255,255,255,0.1);
    border-radius:8px; font-size:12px;
}}
@media (max-width:768px) {{
    .ascii-art {{ font-size:6px; }}
    .results-table {{ font-size:11px; }}
}}
</style>
</head>
<body>
<div class="container">
<div class="ascii-container">
<div class="ascii-art">███████╗████████╗ █████╗ ██████╗ ██╗  ██╗███████╗ ██████╗ ██████╗ ██████╗ ███████╗
██╔════╝╚══██╔══╝██╔══██╗██╔══██╗██║  ██║██╔════╝██╔════╝██╔═══██╗██╔══██╗██╔════╝
███████╗   ██║   ███████║██████╔╝███████║███████╗██║     ██║   ██║██████╔╝█████╗  
╚════██║   ██║   ██╔══██║██╔═══╝ ██╔══██║╚════██║██║     ██║   ██║██╔═══╝ ██╔══╝  
███████║   ██║   ██║  ██║██║     ██║  ██║███████║╚██████╗╚██████╔╝██║     ███████╗
╚══════╝   ╚═╝   ╚═╝  ╚═╝╚═╝     ╚═╝  ╚═╝╚══════╝ ╚═════╝ ╚═════╝ ╚═╝     ╚══════╝</div>
</div>
<div class="quote-container" id="quoteContainer">
<div class="quote-text" id="quoteText"></div>
<div class="quote-author" id="quoteAuthor"></div>
</div>

<div class="report-section">
<h2>📊 Capsule Batch Summary</h2>
<div class="metrics-grid">
<div class="metric-card"><div class="metric-label">Total Samples</div><div class="metric-value">{total}</div></div>
<div class="metric-card"><div class="metric-label">Type 5</div><div class="metric-value">{type5_count}</div></div>
<div class="metric-card"><div class="metric-label">Type 8</div><div class="metric-value">{type8_count}</div></div>
<div class="metric-card"><div class="metric-label">Not Assigned</div><div class="metric-value">{not_assigned_count}</div></div>
<div class="metric-card"><div class="metric-label">Runtime</div><div class="metric-value" style="font-size:20px;">{elapsed:.1f}s</div></div>
</div>
</div>

<div class="report-section">
<h2>🧬 Detailed Results</h2>
<div style="max-height:700px; overflow-y:auto;">
<table class="results-table">
<thead>
<tr>
<th>Sample</th>
<th>Capsule Type</th>
<th>Completeness</th>
<th>Cap5 Score</th>
<th>Cap8 Score</th>
<th>Detected Genes</th>
</tr>
</thead>
<tbody>
{table_rows}
</tbody>
</table>
</div>
</div>

<div class="footer">
<p><strong>STAPHSCOPE</strong> - Capsule Typing</p>
<p class="timestamp">Generated: {time.strftime('%Y-%m-%d %H:%M:%S')}</p>
<p style="margin-top:10px; font-size:12px;">Powered by reverse BLAST against capsule targets</p>
<div class="authorship">
<p><strong>Technical Support &amp; Inquiries:</strong></p>
<p>Original module author: Riccardo Bollini (StaphScan)</p>
<p>Integrated by: Brown Beckley</p>
<p>GitHub: <a href="https://github.com/bbeckley-hub" style="color:#fbbf24;">bbeckley-hub</a></p>
<p>Email: <a href="mailto:brownbeckley94@gmail.com" style="color:#fbbf24;">brownbeckley94@gmail.com</a></p>
<p>Affiliation: University of Ghana Medical School - Department of Medical Biochemistry</p>
</div>
</div>
</div>
<script>
const quotes = [
    {{ text: "The important thing is not to stop questioning. Curiosity has its own reason for existing.", author: "Albert Einstein" }},
    {{ text: "Science is not only a disciple of reason but also one of romance and passion.", author: "Stephen Hawking" }},
    {{ text: "Somewhere, something incredible is waiting to be known.", author: "Carl Sagan" }},
    {{ text: "The good thing about science is that it's true whether or not you believe in it.", author: "Neil deGrasse Tyson" }},
    {{ text: "In science, there are no shortcuts to truth.", author: "Karl Popper" }},
    {{ text: "Science knows no country, because knowledge belongs to humanity.", author: "Louis Pasteur" }},
    {{ text: "The science of today is the technology of tomorrow.", author: "Edward Teller" }},
    {{ text: "Nothing in life is to be feared, it is only to be understood.", author: "Marie Curie" }},
    {{ text: "Research is what I'm doing when I don't know what I'm doing.", author: "Wernher von Braun" }},
    {{ text: "The universe is not required to be in perfect harmony with human ambition.", author: "Carl Sagan" }}
];
const quoteContainer = document.getElementById('quoteContainer');
const quoteText = document.getElementById('quoteText');
const quoteAuthor = document.getElementById('quoteAuthor');
function displayQuote() {{
    quoteContainer.style.opacity = '0';
    setTimeout(() => {{
        const q = quotes[Math.floor(Math.random() * quotes.length)];
        quoteText.textContent = '"' + q.text + '"';
        quoteAuthor.textContent = '— ' + q.author;
        quoteContainer.style.opacity = '1';
    }}, 500);
}}
displayQuote();
setInterval(displayQuote, 10000);
</script>
</body>
</html>
"""
    with open(output_path, "w", encoding="utf-8") as fh:
        fh.write(html)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Reverse-BLAST capsule typing with per-sample and batch reporting",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "-i", "--input", nargs="+", required=True,
        help="One or more input FASTA files",
    )
    parser.add_argument(
        "-d", "--output-dir", required=True,
        help="Top-level output directory; per-sample folders and summary files are written here",
    )
    parser.add_argument(
        "-db_dir", "--database-dir", required=True,
        help="Directory containing the capsule_targets*.fasta database",
    )
    parser.add_argument(
        "-p", "--prefix", default="staphscope_capsule",
        help="Prefix for top-level summary files",
    )
    parser.add_argument(
        "--min-id", type=float, default=90.0,
        help="Minimum percent identity for a hit",
    )
    parser.add_argument(
        "--min-cov", type=float, default=80.0,
        help="Minimum gene coverage for a hit",
    )
    args = parser.parse_args()

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    start = time.time()
    try:
        typer = CapsuleTyper(
            Path(args.database_dir),
            min_id=args.min_id,
            min_cov=args.min_cov,
        )
    except FileNotFoundError as e:
        print(f"Error: {e}")
        sys.exit(1)

    if not typer.check_db():
        print(f"Error: capsule database missing in {args.database_dir}")
        sys.exit(1)

    results: List[SampleResult] = []
    for fasta in args.input:
        fasta_path = Path(fasta)
        if not fasta_path.exists():
            print(f"Skipping missing file: {fasta_path}")
            continue
        print(f"Processing {fasta_path.name} ...")
        result = typer.run(fasta_path)
        results.append(result)
        print(f"  {result.sample}: {result.cap_type} | {result.cap_completeness}")

        sample_out = output_dir / result.sample
        write_sample_json(result, sample_out)
        write_sample_html(result, sample_out)

    if not results:
        print("No samples processed.")
        sys.exit(1)

    write_summary_tsv(results, output_dir / f"{args.prefix}_summary.tsv")
    write_summary_json(results, output_dir / f"{args.prefix}_summary.json")
    write_summary_html(results, output_dir / f"{args.prefix}_summary.html", time.time() - start)

    print(f"\nResults written to {output_dir}")
    print(f"  {args.prefix}_summary.tsv")
    print(f"  {args.prefix}_summary.json")
    print(f"  {args.prefix}_summary.html")
    print(f"  <sample>/<sample>_capsule.html  (per-sample deep-dive)")
    print(f"  <sample>/<sample>_capsule.json  (per-sample JSON)")


if __name__ == "__main__":
    main()

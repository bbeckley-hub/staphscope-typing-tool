#!/usr/bin/env python3
"""
SCCmecFinder_RPet - Reference-based SCCmec Typing Module
Original Module/Logic Author: Robert Petit
Integrated and extended by: Beckley Brown <brownbeckley94@gmail.com>
Date: 2026-09-17
Affiliation: University of Ghana Medical School - Department of Medical Biochemistry

Combines BLAST-based SCCmec target typing and region-based subtyping with
per-sample deep-dive HTML/JSON reports and batch summaries.

Usage:
    python3 SCCmecFinder_RPet.py -i *.fna -d results -db_dir database
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
from typing import Any, Dict, List, Optional, Set, Tuple

import pandas as pd


# ---------------------------------------------------------------------------
# Type-name mapping (RPet gene-content names to CGE-style names)
# So you don't get confused about naming structure!!
# ---------------------------------------------------------------------------

RPET_TO_CGE = {
    "Type I(1B)": "SCCmec_type_I(1B)",
    "Type II(2A)": "SCCmec_type_II(2A)",
    "Type III(3A)": "SCCmec_type_III(3A)",
    "Type IV(2B)": "SCCmec_type_IV(2B)",
    "Type V (5C)": "SCCmec_type_V(5C2)",
    "Type VI(4B)": "SCCmec_type_VI(4B)",
    "Type VII (5C + IS12960D)": "SCCmec_type_VII(5C1)",
    "Type VIII(4A)": "SCCmec_type_VIII(4A)",
    "Type IX (1C)": "SCCmec_type_IX(1C2)",
    "Type X (A1B6-C)": "SCCmec_type_X(7C1)",
    "Type X (A1B6-B)": "SCCmec_type_X(7C1)",
    "Type XI (A1B3)": "SCCmec_type_XI(8E)",
    "Type XII (9C)": "SCCmec_type_XII(9C2)",
    "Type XIII (9A)": "SCCmec_type_XIII(9A)",
    "Type XIV (5A)": "SCCmec_type_XIV(5A)",
    "Type XV (A1B6-A)": "SCCmec_type_XV(7A)",
}

CCR_COMPLEX_NUMBERS = {
    "1": "1",
    "2": "2",
    "3": "3",
    "4": "4",
    "C1": "5",
    "C2": "9",
    "A1B6": "7",
    "A1B3": "8",
}

NOT_ASSIGNED = "Not Assigned"


# ---------------------------------------------------------------------------
# Type exclusions (mirrors the 'excludes' field of Robert's YAML schemas)
# ---------------------------------------------------------------------------

TYPE_EXCLUSIONS: Dict[str, Set[str]] = {
    "Type V (5C)": {"ccr Type 3", "mecI", "IS12960D"},
    "Type VII (5C + IS12960D)": {"ccr Type 3", "mecI"},
    "Type VIII(4A)": {"ccr Type 5"},
    "Type IX (1C)": {"IS1272"},
    "Type X (A1B6-C)": {"mecI"},
    "Type XII (9C)": {"mecI"},
    "Type XIV (5A)": {"ccr Type 3"},
}


# ---------------------------------------------------------------------------
# Data containers
# ---------------------------------------------------------------------------

@dataclass
class BlastHit:
    """A single BLAST alignment record with computed coverage."""
    sseqid: str
    pident: float
    length: int
    slen: int
    sstart: int
    send: int
    evalue: float = 0.0
    bitscore: float = 0.0

    @property
    def coverage(self) -> float:
        """Percent of the subject sequence covered by this alignment."""
        if self.slen <= 0:
            return 0.0
        return (self.length / self.slen) * 100.0


@dataclass
class SampleResult:
    """All typing results and supporting evidence for one sample."""
    sample: str
    input_file: str
    sccmec_type: str = NOT_ASSIGNED
    sccmec_subtype: str = NOT_ASSIGNED
    sccmec_genes: str = NOT_ASSIGNED
    cge_type: str = NOT_ASSIGNED
    mrsa_status: str = "MSSA"
    mrsa_evidence: List[str] = field(default_factory=list)
    mec_class: str = "None"
    ccr_complexes: List[str] = field(default_factory=list)
    found_genes: Set[str] = field(default_factory=set)
    target_hits: List[BlastHit] = field(default_factory=list)
    region_hits: List[BlastHit] = field(default_factory=list)
    warnings: List[str] = field(default_factory=list)


# ---------------------------------------------------------------------------
# SCCmecTyper
# ---------------------------------------------------------------------------

class SCCmecTyper:
    """Reference-based SCCmec typing using BLAST against targets and regions databases."""

    def __init__(self, database_dir: Path):
        """Initialise with a database directory containing targets and regions FASTA files."""
        self.database_dir = Path(database_dir)
        self.targets_fasta = self._locate_db("targets*.fasta")
        self.regions_fasta = self._locate_db("regions*.fasta")
        self.targets_db = self.database_dir / "sccmec_targets_db"
        self.regions_db = self.database_dir / "sccmec_regions_db"

    def _locate_db(self, pattern: str) -> Path:
        """Return the first file matching pattern in the database directory."""
        matches = list(self.database_dir.glob(pattern))
        if not matches:
            raise FileNotFoundError(
                f"Missing database: no file matching '{pattern}' in {self.database_dir}"
            )
        return matches[0]

    def check_db(self) -> bool:
        """Return True if all required database files and BLAST indices exist."""
        if not self.targets_fasta.exists():
            return False
        if not (self.database_dir / "sccmec_targets_db.nhr").exists():
            return False
        if not self.regions_fasta.exists():
            return False
        if not (self.database_dir / "sccmec_regions_db.nhr").exists():
            return False
        return True

    def run(self, assembly_path: Path) -> SampleResult:
        """Run target typing and region subtyping on a single assembly."""
        sample_name = assembly_path.stem

        target_hits = self._blast(assembly_path, self.targets_db, min_ident=90.0, min_cov=80.0)
        region_hits = self._blast(assembly_path, self.regions_db, min_ident=85.0, min_cov=83.0)

        result = SampleResult(
            sample=sample_name,
            input_file=str(assembly_path),
            target_hits=target_hits,
            region_hits=region_hits,
        )

        if not target_hits:
            result.warnings.append("No target genes detected; reporting as MSSA.")
            result.mrsa_status = "MSSA"
            result.mrsa_evidence.append("No mecA or mecC gene detected")
            result.sccmec_type = NOT_ASSIGNED
            result.sccmec_subtype = NOT_ASSIGNED
            result.sccmec_genes = NOT_ASSIGNED
            result.cge_type = NOT_ASSIGNED
            return result

        found_genes = self._extract_genes([h.sseqid for h in target_hits])
        mec_class = self._get_mec_class(found_genes)
        ccr_complexes = self._get_ccr_complex(found_genes)
        detected_aliases = self._detect_aliases(ccr_complexes, mec_class, found_genes)
        sccmec_type = self._assign_type(mec_class, ccr_complexes, found_genes, detected_aliases)
        sccmec_genes = ";".join(sorted(found_genes))
        cge_type = RPET_TO_CGE.get(sccmec_type, sccmec_type)

        base_type = self._parse_base_type(sccmec_type)
        sccmec_subtype = self._determine_subtype(region_hits, base_type)

        result.sccmec_type = sccmec_type
        result.sccmec_subtype = sccmec_subtype
        result.sccmec_genes = sccmec_genes
        result.cge_type = cge_type
        result.mec_class = mec_class
        result.ccr_complexes = ccr_complexes
        result.found_genes = found_genes

        if mec_class == "None":
            result.mrsa_status = "MSSA"
            result.mrsa_evidence.append("No mecA or mecC gene detected")
        else:
            result.mrsa_status = "MRSA"
            if "mecC" in found_genes:
                result.mrsa_evidence.append("mecC gene detected (MRSA determinant)")
            if "mecA" in found_genes:
                result.mrsa_evidence.append("mecA gene detected (MRSA determinant)")

        if result.sccmec_type.startswith("Composite"):
            result.warnings.append("Multiple SCCmec types detected; manual review recommended.")

        return result

    @staticmethod
    def _extract_genes(sseqids: List[str]) -> Set[str]:
        """Map raw BLAST subject identifiers to canonical SCCmec gene names."""
        genes: Set[str] = set()
        for raw_id in sseqids:
            lower = str(raw_id).lower()
            if "ccra1" in lower:
                genes.add("ccrA1")
            elif "ccrb1" in lower:
                genes.add("ccrB1")
            elif "ccra2" in lower:
                genes.add("ccrA2")
            elif "ccrb2" in lower:
                genes.add("ccrB2")
            elif "ccra3" in lower:
                genes.add("ccrA3")
            elif "ccrb3" in lower:
                genes.add("ccrB3")
            elif "ccra4" in lower:
                genes.add("ccrA4")
            elif "ccrb4" in lower:
                genes.add("ccrB4")
            elif "ccrb6" in lower:
                genes.add("ccrB6")
            elif "ccrc1" in lower:
                genes.add("ccrC1")
            elif "ccrc2" in lower:
                genes.add("ccrC2")
            elif "ccrc" in lower:
                genes.add("ccrC1")
            elif "meci" in lower:
                genes.add("mecI")
            elif "mecr1" in lower:
                genes.add("mecR1")
            elif "is1272" in lower:
                genes.add("IS1272")
            elif "is431" in lower:
                genes.add("IS431")
            elif "is12960d" in lower:
                genes.add("IS12960D")
            elif "mecc" in lower:
                genes.add("mecC")
            elif "meca" in lower:
                genes.add("mecA")
            elif "blaz" in lower:
                genes.add("blaZ")
        return genes

    @staticmethod
    def _get_mec_class(genes: Set[str]) -> str:
        """
        Determine mec complex class from detected genes.
        mecC is checked before mecI/mecR1 so that SCCmec XI (which carries both)
        is classified as class E rather than class A.
        """
        if "mecA" not in genes and "mecC" not in genes:
            return "None"
        if "mecC" in genes:
            return "E"
        if "mecI" in genes and "mecR1" in genes:
            return "A"
        if "IS1272" in genes:
            return "B"
        if "IS431" in genes:
            return "C"
        return "Unknown (mecA+)"

    @staticmethod
    def _get_ccr_complex(genes: Set[str]) -> List[str]:
        """Determine ccr gene complex classes from detected genes."""
        complexes: List[str] = []
        if "ccrA1" in genes and "ccrB1" in genes:
            complexes.append("1")
        if "ccrA2" in genes and "ccrB2" in genes:
            complexes.append("2")
        if "ccrA3" in genes and "ccrB3" in genes:
            complexes.append("3")
        if "ccrA4" in genes and "ccrB4" in genes:
            complexes.append("4")
        if "ccrC1" in genes:
            complexes.append("C1")
        if "ccrC2" in genes:
            complexes.append("C2")
        if "ccrA1" in genes and "ccrB6" in genes:
            complexes.append("A1B6")
        if "ccrA1" in genes and "ccrB3" in genes:
            complexes.append("A1B3")
        return complexes

    @staticmethod
    def _detect_aliases(
        ccr_complexes: List[str],
        mec_class: str,
        genes: Set[str],
    ) -> Set[str]:
        """
        Build the set of detected aliases and targets referenced by the exclusion rules.
        Returns ccr alias names, the mec class alias, and individual targets
        (mecI, IS12960D, IS1272) used as exclusion triggers in the YAML schema.
        """
        aliases: Set[str] = set()

        ccr_alias_map = {
            "1": "ccr Type 1",
            "2": "ccr Type 2",
            "3": "ccr Type 3",
            "4": "ccr Type 4",
            "C1": "ccr Type 5",
            "C2": "ccr Type 9",
            "A1B6": "ccr Type 7",
            "A1B3": "ccr Type 8",
        }
        for ccr in ccr_complexes:
            if ccr in ccr_alias_map:
                aliases.add(ccr_alias_map[ccr])

        if mec_class != "None":
            aliases.add(f"mec Class {mec_class}")

        for gene in ("mecI", "IS12960D", "IS1272"):
            if gene in genes:
                aliases.add(gene)

        return aliases

    @staticmethod
    def _assign_type(
        mec_class: str,
        ccr_list: List[str],
        all_genes: Set[str],
        detected_aliases: Set[str],
    ) -> str:
        """
        Assign SCCmec type from mec class and ccr complex combinations.
        Types whose exclusion rules are triggered by detected aliases or
        targets are discarded before the final result is assembled.
        """
        if mec_class == "None":
            return NOT_ASSIGNED
        if not ccr_list:
            return f"Orphan {mec_class} (No ccr)"

        types_found: List[str] = []
        for ccr in ccr_list:
            candidate: Optional[str] = None

            if ccr == "1" and mec_class == "B":
                candidate = "Type I(1B)"
            elif ccr == "2" and mec_class == "A":
                candidate = "Type II(2A)"
            elif ccr == "3" and mec_class == "A":
                candidate = "Type III(3A)"
            elif ccr == "2" and mec_class == "B":
                candidate = "Type IV(2B)"
            elif ccr == "4" and mec_class == "B":
                candidate = "Type VI(4B)"
            elif ccr == "4" and mec_class == "A":
                candidate = "Type VIII(4A)"
            elif ccr == "C1" and mec_class == "C":
                if "IS12960D" in all_genes:
                    candidate = "Type VII (5C + IS12960D)"
                else:
                    candidate = "Type V (5C)"
            elif ccr == "C1" and mec_class == "A":
                candidate = "Type XIV (5A)"
            elif ccr == "1" and mec_class == "C":
                candidate = "Type IX (1C)"
            elif ccr == "A1B6":
                if mec_class in ("C", "B"):
                    candidate = "Type X (A1B6-C)"
                elif mec_class == "A":
                    candidate = "Type XV (A1B6-A)"
            elif ccr == "A1B3":
                candidate = "Type XI (A1B3)"
            elif ccr == "C2" and mec_class == "C":
                candidate = "Type XII (9C)"
            elif ccr == "C2" and mec_class == "A":
                candidate = "Type XIII (9A)"
            else:
                candidate = f"Type ? ({ccr}-{mec_class})"

            if candidate is None:
                continue

            exclusions = TYPE_EXCLUSIONS.get(candidate, set())
            if exclusions and exclusions & detected_aliases:
                continue

            types_found.append(candidate)

        if not types_found:
            return NOT_ASSIGNED
        if len(types_found) == 1:
            return types_found[0]
        return f"Composite ({' + '.join(types_found)})"

    @staticmethod
    def _parse_base_type(sccmec_type: str) -> Optional[str]:
        """Extract roman numeral base type from a type string."""
        match = re.search(r"Type\s+([IVX]+)", sccmec_type)
        return match.group(1) if match else None

    def _determine_subtype(self, region_hits: List[BlastHit], base_type: Optional[str]) -> str:
        """Return the closest matching subtype from region BLAST hits."""
        if not region_hits or not base_type:
            return NOT_ASSIGNED

        subtype_candidates = []
        for hit in region_hits:
            name = str(hit.sseqid).split("|")[-1] if "|" in str(hit.sseqid) else str(hit.sseqid)
            if re.match(f"^{base_type}[a-z]*$", name):
                subtype_candidates.append((name, hit.coverage, hit.pident))

        if not subtype_candidates:
            return NOT_ASSIGNED

        subtype_candidates.sort(key=lambda x: (x[1], x[2]), reverse=True)
        return subtype_candidates[0][0]

    @staticmethod
    def _blast(assembly_path: Path, db_path: Path, min_ident: float, min_cov: float) -> List[BlastHit]:
        """
        Run blastn against a database and return filtered hits.
        Coverage is computed by merging overlapping alignments per reference
        sequence, making the result robust to fragmented assemblies.
        """
        cmd = [
            "blastn",
            "-query", str(assembly_path),
            "-db", str(db_path),
            "-outfmt", "6 sseqid pident length slen sstart send evalue bitscore",
            "-max_target_seqs", "100",
        ]
        try:
            res = subprocess.run(cmd, capture_output=True, text=True)
            if not res.stdout or not res.stdout.strip():
                return []

            df = pd.read_csv(
                io.StringIO(res.stdout),
                sep="\t",
                names=["sseqid", "pident", "length", "slen", "sstart", "send", "evalue", "bitscore"],
            )
            if df.empty:
                return []

            df = df[df["pident"] >= min_ident]
            if df.empty:
                return []

            merged_hits: Dict[str, Dict] = {}
            for sseqid, group in df.groupby("sseqid"):
                slen = int(group["slen"].iloc[0])
                intervals = [
                    sorted([int(row["sstart"]), int(row["send"])])
                    for _, row in group.iterrows()
                ]
                intervals.sort()

                merged: List[List[int]] = []
                for interval in intervals:
                    if not merged:
                        merged.append(interval)
                    elif interval[0] <= merged[-1][1] + 1:
                        merged[-1][1] = max(merged[-1][1], interval[1])
                    else:
                        merged.append(interval)

                cov_len = sum(end - start + 1 for start, end in merged)
                cov = (cov_len / slen) * 100
                weighted_ident = (group["pident"] * group["length"]).sum() / group["length"].sum()
                best = group.loc[group["bitscore"].idxmax()]

                merged_hits[sseqid] = {
                    "sseqid": sseqid,
                    "pident": weighted_ident,
                    "length": cov_len,
                    "slen": slen,
                    "sstart": int(best["sstart"]),
                    "send": int(best["send"]),
                    "evalue": float(best["evalue"]),
                    "bitscore": float(best["bitscore"]),
                    "cov": cov,
                }

            hits = [
                BlastHit(
                    sseqid=v["sseqid"],
                    pident=v["pident"],
                    length=v["length"],
                    slen=v["slen"],
                    sstart=v["sstart"],
                    send=v["send"],
                    evalue=v["evalue"],
                    bitscore=v["bitscore"],
                )
                for v in merged_hits.values()
                if v["cov"] >= min_cov
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
        "sccmec_type": result.sccmec_type,
        "sccmec_subtype": result.sccmec_subtype,
        "cge_equivalent_type": result.cge_type,
        "sccmec_genes": result.sccmec_genes,
        "mrsa_status": result.mrsa_status,
        "mrsa_evidence": result.mrsa_evidence,
        "mec_class": result.mec_class,
        "ccr_complexes": result.ccr_complexes,
        "found_genes": sorted(result.found_genes),
        "target_hits": [
            {
                "sseqid": h.sseqid,
                "pident": round(h.pident, 2),
                "length": h.length,
                "slen": h.slen,
                "coverage": round(h.coverage, 2),
                "bitscore": h.bitscore,
            }
            for h in result.target_hits
        ],
        "region_hits": [
            {
                "sseqid": h.sseqid,
                "pident": round(h.pident, 2),
                "length": h.length,
                "slen": h.slen,
                "coverage": round(h.coverage, 2),
                "bitscore": h.bitscore,
            }
            for h in result.region_hits
        ],
        "warnings": result.warnings,
    }
    path = output_dir / f"{result.sample}_sccmec_rpet.json"
    with open(path, "w") as fh:
        json.dump(payload, fh, indent=2)
    return path


def _hit_rows_html(hits: List[BlastHit]) -> str:
    """Render BLAST hit rows as HTML table rows with coverage colouring."""
    if not hits:
        return "<tr><td colspan='6' style='text-align:center;color:#9ca3af;'>No hits above thresholds</td></tr>"
    rows = []
    for i, h in enumerate(hits):
        if i == 0:
            cls = "highlight-red"
        elif h.coverage >= 70.0:
            cls = "highlight-orange"
        elif h.coverage >= 50.0:
            cls = "highlight-yellow"
        else:
            cls = ""
        rows.append(
            f"<tr class='{cls}'>"
            f"<td>{i + 1}</td>"
            f"<td><strong>{h.sseqid}</strong></td>"
            f"<td>{h.pident:.2f}</td>"
            f"<td>{h.length}</td>"
            f"<td>{h.slen}</td>"
            f"<td>{h.coverage:.2f}</td>"
            f"</tr>"
        )
    return "\n".join(rows)


def write_sample_html(result: SampleResult, output_dir: Path) -> Path:
    """Write a per-sample deep-dive HTML report with all BLAST hits and evidence."""
    output_dir.mkdir(parents=True, exist_ok=True)

    status_class = "status-mrsa" if result.mrsa_status == "MRSA" else "status-mssa"
    evidence_items = "".join(f"<li>{e}</li>" for e in result.mrsa_evidence)
    warning_html = ""
    if result.warnings:
        warning_html = "<div class='warning-box'><strong>⚠️ Warnings:</strong><ul>" + \
            "".join(f"<li>{w}</li>" for w in result.warnings) + "</ul></div>"

    gene_cards = "".join(
        f"<div class='gene-card'>{g}</div>" for g in sorted(result.found_genes)
    )
    if not gene_cards:
        gene_cards = "<p style='color:#6b7280;'>No genes detected above thresholds.</p>"

    ccr_display = ", ".join(result.ccr_complexes) if result.ccr_complexes else NOT_ASSIGNED
    target_rows = _hit_rows_html(result.target_hits)
    region_rows = _hit_rows_html(result.region_hits)

    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>STAPHSCOPE - {result.sample} SCCmec Report</title>
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
.status-mrsa {{ background:#dc2626; color:white; }}
.status-mssa {{ background:#16a34a; color:white; }}
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
    display:grid; grid-template-columns:repeat(auto-fill, minmax(150px, 1fr));
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
<div class="metric-card"><div class="metric-label">SCCmec Type</div><div class="metric-value" style="font-size:18px;">{result.sccmec_type}</div></div>
<div class="metric-card"><div class="metric-label">Subtype</div><div class="metric-value" style="font-size:18px;">{result.sccmec_subtype}</div></div>
<div class="metric-card"><div class="metric-label">Genes Found</div><div class="metric-value">{len(result.found_genes)}</div></div>
</div>
</div>

<div class="report-section">
<h2>🦠 MRSA/MSSA Classification</h2>
<div class="status-badge {status_class}">{result.mrsa_status}</div>
<h3>Evidence</h3>
<ul style="margin-left:20px; line-height:1.8;">{evidence_items}</ul>
{warning_html}
</div>

<div class="report-section">
<h2>🧬 SCCmec Typing</h2>
<p><strong>Type:</strong> {result.sccmec_type}</p>
<p><strong>Subtype:</strong> {result.sccmec_subtype}</p>
<p><strong>CGE Equivalent:</strong> {result.cge_type}</p>
<p><strong>MEC Class:</strong> {result.mec_class}</p>
<p><strong>CCR Complexes:</strong> {ccr_display}</p>
<p><strong>Detected Genes:</strong> {result.sccmec_genes}</p>
</div>

<div class="report-section">
<h2>🧪 Gene Detection Summary</h2>
<div class="gene-grid">{gene_cards}</div>
</div>

<div class="report-section">
<h2>🎯 Target BLAST Hits (mec/ccr genes)</h2>
<div style="max-height:500px; overflow-y:auto;">
<table class="detailed-table">
<thead>
<tr><th>Rank</th><th>Target</th><th>% Identity</th><th>Alignment Length</th><th>Subject Length</th><th>Coverage %</th></tr>
</thead>
<tbody>
{target_rows}
</tbody>
</table>
</div>
<p style="margin-top:10px; font-size:12px; color:#666;">
● Top hit in RED &nbsp;|&nbsp; <span style="color:#ea580c;">● Coverage ≥70% in ORANGE</span> &nbsp;|&nbsp; <span style="color:#d97706;">● Coverage ≥50% in YELLOW</span>
</p>
</div>

<div class="report-section">
<h2>🔍 Region BLAST Hits (subtyping)</h2>
<div style="max-height:500px; overflow-y:auto;">
<table class="detailed-table">
<thead>
<tr><th>Rank</th><th>Region</th><th>% Identity</th><th>Alignment Length</th><th>Subject Length</th><th>Coverage %</th></tr>
</thead>
<tbody>
{region_rows}
</tbody>
</table>
</div>
</div>

<div class="footer">
<p><strong>STAPHSCOPE</strong> - SCCmec Typing (RPet module)</p>
<p class="timestamp">Generated: {time.strftime('%Y-%m-%d %H:%M:%S')}</p>
<div class="authorship">
<p><strong>Technical Support &amp; Inquiries:</strong></p>
<p>Original module author: Robert Petit</p>
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
    path = output_dir / f"{result.sample}_sccmec_rpet.html"
    with open(path, "w", encoding="utf-8") as fh:
        fh.write(html)
    return path


# ---------------------------------------------------------------------------
# Top-level summary writers
# ---------------------------------------------------------------------------

def write_summary_tsv(results: List[SampleResult], output_path: Path) -> None:
    """Write a batch TSV summary with one row per sample."""
    columns = [
        "sample", "sccmec_type", "cge_equivalent", "sccmec_subtype",
        "mrsa_status", "mec_class", "ccr_complexes", "sccmec_genes",
    ]
    with open(output_path, "w") as fh:
        fh.write("\t".join(columns) + "\n")
        for r in results:
            row = [
                r.sample,
                r.sccmec_type,
                r.cge_type,
                r.sccmec_subtype,
                r.mrsa_status,
                r.mec_class,
                ",".join(r.ccr_complexes) if r.ccr_complexes else NOT_ASSIGNED,
                r.sccmec_genes,
            ]
            fh.write("\t".join(str(v) for v in row) + "\n")


def write_summary_json(results: List[SampleResult], output_path: Path) -> None:
    """Write a combined JSON summary for all samples."""
    payload = []
    for r in results:
        payload.append({
            "sample": r.sample,
            "sccmec_type": r.sccmec_type,
            "cge_equivalent_type": r.cge_type,
            "sccmec_subtype": r.sccmec_subtype,
            "mrsa_status": r.mrsa_status,
            "mrsa_evidence": r.mrsa_evidence,
            "mec_class": r.mec_class,
            "ccr_complexes": r.ccr_complexes,
            "sccmec_genes": r.sccmec_genes,
            "warnings": r.warnings,
        })
    with open(output_path, "w") as fh:
        json.dump(payload, fh, indent=2)


def write_summary_html(results: List[SampleResult], output_path: Path, elapsed: float) -> None:
    """Write a StaphScope-styled batch summary HTML report."""
    total = len(results)
    mrsa_count = sum(1 for r in results if r.mrsa_status == "MRSA")
    mssa_count = total - mrsa_count
    typed_count = sum(
        1 for r in results
        if r.sccmec_type not in (NOT_ASSIGNED, "-")
        and not r.sccmec_type.startswith("Orphan")
    )

    rows = []
    for r in results:
        status_class = "status-mrsa" if r.mrsa_status == "MRSA" else "status-mssa"
        ccr = ", ".join(r.ccr_complexes) if r.ccr_complexes else NOT_ASSIGNED
        rows.append(
            f"<tr>"
            f"<td><strong>{r.sample}</strong></td>"
            f"<td>{r.sccmec_type}</td>"
            f"<td>{r.cge_type}</td>"
            f"<td>{r.sccmec_subtype}</td>"
            f"<td><span class='{status_class}'>{r.mrsa_status}</span></td>"
            f"<td>{r.mec_class}</td>"
            f"<td>{ccr}</td>"
            f"</tr>"
        )
    table_rows = "\n".join(rows)

    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>STAPHSCOPE - SCCmec Batch Summary (RPet)</title>
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
.status-mrsa {{
    background:#dc2626; color:white; padding:4px 8px;
    border-radius:4px; font-weight:bold; font-size:11px;
}}
.status-mssa {{
    background:#16a34a; color:white; padding:4px 8px;
    border-radius:4px; font-weight:bold; font-size:11px;
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
<h2>📊 SCCmec Batch Summary (RPet module)</h2>
<div class="metrics-grid">
<div class="metric-card"><div class="metric-label">Total Samples</div><div class="metric-value">{total}</div></div>
<div class="metric-card"><div class="metric-label">MRSA</div><div class="metric-value">{mrsa_count}</div></div>
<div class="metric-card"><div class="metric-label">MSSA</div><div class="metric-value">{mssa_count}</div></div>
<div class="metric-card"><div class="metric-label">SCCmec Typed</div><div class="metric-value">{typed_count}</div></div>
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
<th>SCCmec Type</th>
<th>CGE Equivalent</th>
<th>Subtype</th>
<th>MRSA Status</th>
<th>MEC Class</th>
<th>CCR Complexes</th>
</tr>
</thead>
<tbody>
{table_rows}
</tbody>
</table>
</div>
</div>

<div class="footer">
<p><strong>STAPHSCOPE</strong> - SCCmec Typing (RPet module)</p>
<p class="timestamp">Generated: {time.strftime('%Y-%m-%d %H:%M:%S')}</p>
<p style="margin-top:10px; font-size:12px;">Powered by BLAST-based target and region typing</p>
<div class="authorship">
<p><strong>Technical Support &amp; Inquiries:</strong></p>
<p>Original module author: Robert Petit</p>
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
        description="Reference-based SCCmec typing with per-sample and batch reporting (RPet module)",
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
        help="Directory containing targets*.fasta, regions*.fasta and BLAST indices",
    )
    parser.add_argument(
        "-p", "--prefix", default="staphscope_sccmec_rpet",
        help="Prefix for top-level summary files",
    )
    args = parser.parse_args()

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    start = time.time()
    typer = SCCmecTyper(Path(args.database_dir))

    if not typer.check_db():
        print(f"Error: SCCmec databases missing or not indexed in {args.database_dir}")
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
        print(f"  {result.sample}: {result.sccmec_type} | {result.mrsa_status}")

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
    print(f"  <sample>/<sample>_sccmec_rpet.html  (per-sample deep-dive)")
    print(f"  <sample>/<sample>_sccmec_rpet.json  (per-sample JSON)")


if __name__ == "__main__":
    main()
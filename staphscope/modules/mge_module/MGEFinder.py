#!/usr/bin/env python3
"""
MGEFinder - Mobile Genetic Element Detection Module for StaphScope
Integrated by: Beckley Brown <brownbeckley94@gmail.com>
Date: 2026-09-17
Affiliation: University of Ghana Medical School - Department of Medical Biochemistry

Detects mobile genetic elements using Prodigal for ORF prediction and DIAMOND
against the mobileOG-db, then joins hits to the mobileOG-db metadata for
categorisation by functional class (integration/excision, transfer, phage,
stability/defense, replication/repair) and by source database class
(insertion sequences, integrative elements, plasmids, bacteriophages, multiple).

Supports automatic CPU and RAM detection with tiered scaling, plus two-level
parallelism: multiple genomes processed concurrently while DIAMOND threads are
split evenly across workers.

Per-sample deep-dive HTML/JSON reports plus top-level batch summaries.

Acknowledgements:
    - mobileOG-db team (Brown et al.) for the curated database
    - Doug Hyatt for Prodigal
    - Benjamin Buchfink for DIAMOND
    - James I. Mullet and Connor L. Brown for the mobileOG-pl parsing scripts

Database layout expected inside -db_dir:
    *.dmnd                            (DIAMOND database)
    mobileOG-db-beatrix-*-All.csv     (metadata)

Usage:
    python3 MGEFinder.py -i *.fna -d mge_results -db_dir database
    python3 MGEFinder.py -i *.fna -d mge_results -db_dir database --cpus 8
    python3 MGEFinder.py -i *.fna -d mge_results -db_dir mobileOG-db/beatrix-1-6_v1_all
"""

import argparse
import io
import json
import os
import re
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional, Set

import pandas as pd
import psutil


# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

NOT_ASSIGNED = "Not Assigned"

CATEGORY_ORDER = [
    "integration/excision",
    "transfer",
    "stability/transfer/defense",
    "phage",
    "replication/recombination/repair",
]

CATEGORY_LABELS = {
    "integration/excision": "Integration / Excision",
    "transfer": "Transfer (Conjugation / Competence)",
    "stability/transfer/defense": "Stability / Defense",
    "phage": "Phage",
    "replication/recombination/repair": "Replication / Repair",
}

SOURCE_CLASS_ORDER = [
    "Insertion sequences",
    "Integrative elements",
    "Plasmids",
    "Bacteriophages",
    "Multiple",
]

SOURCE_COLUMN_MAP = {
    "ISFinder": "Insertion sequences",
    "AICE": "Integrative elements",
    "ICE": "Integrative elements",
    "CIME": "Integrative elements",
    "IME": "Integrative elements",
    "immedb": "Integrative elements",
    "COMPASS": "Plasmids",
    "PlasmidRefSeq": "Plasmids",
    "ACLAME": "Multiple",
    "pVOG": "Bacteriophages",
    "GPD": "Bacteriophages",
}

KEY_GENES = {
    "mecA", "mecC", "ccrA", "ccrB", "ccrC",
    "tnp", "tnpA", "tnpB", "tnpC", "tnp1", "tnpM",
    "int", "intI", "xerC", "xerD",
    "traG", "traR", "traP", "virB", "virD", "virE",
    "mobA", "relaxase", "pre",
    "lukS", "lukF", "tst", "sea", "seb", "sec",
    "lytA", "ejh", "kilA", "ant", "ter", "mtp", "RINB",
    "ardA", "mazE", "mazF", "yefM", "yoeB",
}


# ---------------------------------------------------------------------------
# Data containers
# ---------------------------------------------------------------------------

@dataclass
class BlastHit:
    """A single DIAMOND alignment record with joined metadata and ORF context."""
    qseqid: str
    qtitle: str
    sseqid: str
    mobileog_id: str
    pident: float
    length: int
    qlen: int
    evalue: float = 0.0
    bitscore: float = 0.0
    orf_start: Optional[int] = None
    orf_stop: Optional[int] = None
    strand: str = ""
    partial: str = ""
    start_codon: str = ""
    rbs_motif: str = ""
    rbs_spacer: str = ""
    gc_content: Optional[float] = None
    gene_name: str = ""
    annotation: str = ""
    major_category: str = "other"
    minor_category: str = ""
    database: str = ""
    evidence: str = ""
    taxonomy: str = ""
    source_databases: List[str] = field(default_factory=list)
    source_classes: List[str] = field(default_factory=list)

    @property
    def coverage(self) -> float:
        """Percent of the query protein covered by this alignment."""
        if self.qlen <= 0:
            return 0.0
        return (self.length / self.qlen) * 100.0

    @property
    def specific_contig(self) -> str:
        """Contig identifier derived from the Prodigal ORF name."""
        if "_" in self.qseqid:
            return "_".join(self.qseqid.split("_")[:-1])
        return self.qseqid


@dataclass
class SampleResult:
    """All MGE detection results and supporting evidence for one sample."""
    sample: str
    input_file: str
    faa_file: str = ""
    total_hits: int = 0
    category_hits: Dict[str, List[BlastHit]] = field(default_factory=dict)
    key_gene_hits: List[BlastHit] = field(default_factory=list)
    raw_hits: List[BlastHit] = field(default_factory=list)
    contig_purity: List[Dict[str, Any]] = field(default_factory=list)
    evidence: List[str] = field(default_factory=list)
    warnings: List[str] = field(default_factory=list)

    @property
    def category_counts(self) -> Dict[str, int]:
        """Number of hits per functional category (ORF-level)."""
        counts: Dict[str, int] = {}
        for cat, hits in self.category_hits.items():
            counts[cat] = len(hits)
        return counts

    @property
    def source_class_counts(self) -> Dict[str, int]:
        """Number of hits per source-database class (multi-label)."""
        counts = {cls: 0 for cls in SOURCE_CLASS_ORDER}
        for h in self.raw_hits:
            for cls in (h.source_classes or ["Multiple"]):
                if cls in counts:
                    counts[cls] += 1
        return counts

    @property
    def elements_summary(self) -> str:
        """Semicolon-joined list of detected mobileOG gene names."""
        names: Set[str] = set()
        for hits in self.category_hits.values():
            for h in hits:
                if h.gene_name and h.gene_name != "NA":
                    names.add(h.gene_name)
        return ";".join(sorted(names)) if names else NOT_ASSIGNED


# ---------------------------------------------------------------------------
# MGETyper
# ---------------------------------------------------------------------------

class MGETyper:
    """Mobile genetic element detection via Prodigal + DIAMOND + mobileOG-db metadata."""

    def __init__(
        self,
        database_dir: Path,
        min_id: float = 80.0,
        min_cov: float = 60.0,
        max_evalue: float = 1e-20,
        cpus: Optional[int] = None,
    ):
        """Initialise with a database directory containing the DIAMOND db and metadata."""
        self.database_dir = Path(database_dir)
        self.min_id = min_id
        self.min_cov = min_cov
        self.max_evalue = max_evalue
        self.available_ram = self._get_available_ram()
        self.cpus = self._calculate_optimal_cpus(cpus)
        self.dmnd = self._locate_dmnd()
        self.metadata_csv = self._locate_metadata()
        self.metadata = self._load_metadata()

    # -- resource detection ----------------------------------------------

    @staticmethod
    def _get_available_ram() -> float:
        """Return available RAM in gigabytes."""
        try:
            return psutil.virtual_memory().available / (1024 ** 3)
        except Exception:
            return 8.0

    @staticmethod
    def _calculate_optimal_cpus(user_cpus: Optional[int] = None) -> int:
        """Return the number of CPU cores to use after applying tiered scaling."""
        if user_cpus is not None:
            return max(1, user_cpus)
        try:
            total = psutil.cpu_count(logical=False) or os.cpu_count() or 2
            if total <= 4:
                optimal = total
            elif total <= 8:
                optimal = total - 1
            elif total <= 16:
                optimal = max(8, total - 1)
            elif total <= 32:
                optimal = max(16, total - 3)
            else:
                optimal = min(32, int(total * 0.95))
            return max(1, min(optimal, total))
        except Exception:
            return os.cpu_count() or 4

    # -- database loading -------------------------------------------------

    def _locate_dmnd(self) -> Path:
        """Return the mobileOG DIAMOND database file."""
        matches = sorted(self.database_dir.glob("*.dmnd"))
        if not matches:
            raise FileNotFoundError(
                f"Missing database: no .dmnd file in {self.database_dir}"
            )
        return matches[0]

    def _locate_metadata(self) -> Path:
        """Return the mobileOG-db metadata CSV file."""
        matches = sorted(self.database_dir.glob("*metadata*.csv")) or \
                  sorted(self.database_dir.glob("mobileOG-db-*.csv"))
        if not matches:
            raise FileNotFoundError(
                f"Missing metadata: no metadata CSV in {self.database_dir}"
            )
        return matches[0]

    def _load_metadata(self) -> Dict[str, Dict[str, Any]]:
        """Load the mobileOG-db metadata CSV into a dict keyed by mobileOG ID."""
        try:
            df = pd.read_csv(self.metadata_csv, low_memory=False)
        except Exception:
            return {}

        available_cols = set(df.columns)
        source_cols = [c for c in SOURCE_COLUMN_MAP if c in available_cols]

        index: Dict[str, Dict[str, Any]] = {}
        for _, row in df.iterrows():
            mid = str(row.get("mobileOG Entry Name", "")).strip()
            if not mid:
                continue

            source_databases: List[str] = []
            source_classes: Set[str] = set()
            for col in source_cols:
                val = row.get(col)
                try:
                    if pd.notna(val) and float(val) > 0:
                        source_databases.append(col)
                        source_classes.add(SOURCE_COLUMN_MAP[col])
                except (ValueError, TypeError):
                    pass

            index[mid] = {
                "name": str(row.get("Name", "") or "").strip(),
                "annotation": str(row.get("Manual Annotation", "") or "").strip(),
                "major_category": str(row.get("Major mobileOG Category", "other") or "other").strip(),
                "minor_category": str(row.get("Minor mobileOG Categories", "") or "").strip(),
                "database": str(row.get("Database", "") or "").strip(),
                "evidence": str(row.get("Evidence", "") or "").strip(),
                "taxonomy": str(row.get("Taxonomy", "") or "").strip(),
                "source_databases": source_databases,
                "source_classes": sorted(source_classes),
            }
        return index

    def check_db(self) -> bool:
        """Return True if both the DIAMOND database and metadata are available."""
        return self.dmnd.exists() and bool(self.metadata)

    # -- main per-sample entry point --------------------------------------

    def run(
        self,
        assembly_path: Path,
        annotations_dir: Path,
        diamond_threads: Optional[int] = None,
    ) -> SampleResult:
        """Annotate (if needed), search, and type one assembly."""
        sample_name = assembly_path.stem
        threads = max(1, diamond_threads or self.cpus)

        faa_path = self._annotate(assembly_path, annotations_dir)
        if faa_path is None:
            result = SampleResult(sample=sample_name, input_file=str(assembly_path))
            result.warnings.append("ORF prediction failed; unable to run MGE detection.")
            return result

        raw_hits = self._diamond_search(faa_path, threads)
        result = SampleResult(
            sample=sample_name,
            input_file=str(assembly_path),
            faa_file=str(faa_path),
            raw_hits=raw_hits,
            total_hits=len(raw_hits),
        )

        if not raw_hits:
            result.warnings.append("No mobile genetic element proteins detected above thresholds.")
            result.evidence.append("No mobileOG-db matches found")
            return result

        for hit in raw_hits:
            result.category_hits.setdefault(hit.major_category, []).append(hit)
            if self._is_key_gene(hit):
                result.key_gene_hits.append(hit)

        result.contig_purity = self._compute_contig_purity(raw_hits)

        result.evidence.append(f"{result.total_hits} mobileOG-associated protein(s) matched")

        functional_counts = result.category_counts
        for cat in CATEGORY_ORDER:
            n = functional_counts.get(cat, 0)
            if n:
                result.evidence.append(f"{CATEGORY_LABELS[cat]}: {n}")

        source_counts = result.source_class_counts
        for cls in SOURCE_CLASS_ORDER:
            n = source_counts.get(cls, 0)
            if n:
                result.evidence.append(f"Source class {cls}: {n}")

        if result.key_gene_hits:
            key_names = sorted({h.gene_name for h in result.key_gene_hits if h.gene_name})
            if key_names:
                result.evidence.append(f"Key genes: {', '.join(key_names)}")

        if functional_counts.get("transfer", 0) >= 5:
            result.warnings.append(
                "Transfer-associated protein signatures detected; these may be consistent with "
                "conjugative machinery. Element-level assignment as an ICE or conjugative plasmid "
                "requires genomic context and co-localization of compatible genes."
            )

        return result

    # -- annotation -------------------------------------------------------

    @staticmethod
    def _annotate(assembly_path: Path, annotations_dir: Path) -> Optional[Path]:
        """Run Prodigal to predict ORFs, caching the .faa in annotations_dir."""
        annotations_dir.mkdir(parents=True, exist_ok=True)
        faa_path = annotations_dir / f"{assembly_path.stem}.faa"

        if faa_path.exists() and faa_path.stat().st_mtime >= assembly_path.stat().st_mtime:
            return faa_path

        cmd = [
            "prodigal",
            "-i", str(assembly_path),
            "-a", str(faa_path),
            "-p", "single",
            "-q",
            "-o", "/dev/null",
        ]
        try:
            subprocess.run(cmd, check=True, capture_output=True)
            return faa_path if faa_path.exists() else None
        except (subprocess.CalledProcessError, FileNotFoundError):
            return None

    # -- search -----------------------------------------------------------

    def _diamond_search(self, faa_path: Path, threads: int) -> List[BlastHit]:
        """Run DIAMOND blastp against the mobileOG-db and return joined hits."""
        cmd = [
            "diamond", "blastp",
            "--query", str(faa_path),
            "--db", str(self.dmnd),
            "--outfmt", "6", "qseqid", "qtitle", "sseqid", "pident", "length", "qlen", "evalue", "bitscore",
            "--max-target-seqs", "1",
            "--evalue", str(self.max_evalue),
            "--threads", str(threads),
            "--quiet",
        ]
        try:
            res = subprocess.run(cmd, capture_output=True, text=True)
            if not res.stdout or not res.stdout.strip():
                return []

            df = pd.read_csv(
                io.StringIO(res.stdout),
                sep="\t",
                names=["qseqid", "qtitle", "sseqid", "pident", "length", "qlen", "evalue", "bitscore"],
                dtype={"qtitle": str, "sseqid": str, "qseqid": str},
            )
            if df.empty:
                return []

            df["cov"] = (df["length"] / df["qlen"]) * 100
            df = df[(df["pident"] >= self.min_id) & (df["cov"] >= self.min_cov)]
            if df.empty:
                return []

            hits: List[BlastHit] = []
            for _, row in df.iterrows():
                mid = self._extract_mobileog_id(row["sseqid"])
                meta = self.metadata.get(mid, {})
                parsed = self._parse_qtitle(str(row["qtitle"]))

                hits.append(BlastHit(
                    qseqid=str(row["qseqid"]),
                    qtitle=str(row["qtitle"]),
                    sseqid=str(row["sseqid"]),
                    mobileog_id=mid,
                    pident=float(row["pident"]),
                    length=int(row["length"]),
                    qlen=int(row["qlen"]),
                    evalue=float(row["evalue"]),
                    bitscore=float(row["bitscore"]),
                    orf_start=parsed.get("orf_start"),
                    orf_stop=parsed.get("orf_stop"),
                    strand=parsed.get("strand", ""),
                    partial=parsed.get("partial", ""),
                    start_codon=parsed.get("start_codon", ""),
                    rbs_motif=parsed.get("rbs_motif", ""),
                    rbs_spacer=parsed.get("rbs_spacer", ""),
                    gc_content=parsed.get("gc_content"),
                    gene_name=meta.get("name", ""),
                    annotation=meta.get("annotation", ""),
                    major_category=meta.get("major_category", "other"),
                    minor_category=meta.get("minor_category", ""),
                    database=meta.get("database", ""),
                    evidence=meta.get("evidence", ""),
                    taxonomy=meta.get("taxonomy", ""),
                    source_databases=list(meta.get("source_databases", [])),
                    source_classes=list(meta.get("source_classes", [])),
                ))

            hits.sort(key=lambda h: (h.bitscore, h.pident), reverse=True)
            return hits

        except FileNotFoundError:
            return []
        except Exception:
            return []

    @staticmethod
    def _extract_mobileog_id(sseqid: str) -> str:
        """Extract the mobileOG ID (first pipe-delimited field) from a DIAMOND subject ID."""
        return str(sseqid).split("|")[0].strip()

    @staticmethod
    def _parse_qtitle(qtitle: str) -> Dict[str, Any]:
        """
        Parse a Prodigal ORF header line to recover genomic coordinates and ORF context.
        Expected format:
            {contig}_{orfnum} # {start} # {stop} # {strand} # ID={id};partial=...;start_type=...;...
        """
        result: Dict[str, Any] = {
            "orf_start": None,
            "orf_stop": None,
            "strand": "",
            "partial": "",
            "start_codon": "",
            "rbs_motif": "",
            "rbs_spacer": "",
            "gc_content": None,
        }

        parts = qtitle.split(" # ")
        if len(parts) < 5:
            return result

        try:
            result["orf_start"] = int(parts[1].strip())
            result["orf_stop"] = int(parts[2].strip())
            result["strand"] = "+" if parts[3].strip() == "1" else "-"
        except (ValueError, IndexError):
            return result

        for kv in parts[4].split(";"):
            if "=" not in kv:
                continue
            k, v = kv.split("=", 1)
            k = k.strip().lower()
            v = v.strip()
            if k == "partial":
                result["partial"] = v
            elif k == "start_type":
                result["start_codon"] = v
            elif k == "rbs_motif":
                result["rbs_motif"] = v
            elif k == "rbs_spacer":
                result["rbs_spacer"] = v
            elif k == "gc_cont":
                try:
                    result["gc_content"] = float(v)
                except ValueError:
                    pass

        return result

    @staticmethod
    def _compute_contig_purity(hits: List[BlastHit]) -> List[Dict[str, Any]]:
        """
        Compute per-contig source-class composition.
        Every (ORF hit x source class) pair contributes one count, so a hit
        assigned to multiple source classes is counted once per class.
        """
        if not hits:
            return []

        counts: Dict[str, Dict[str, int]] = {}
        for h in hits:
            contig = h.specific_contig
            if contig not in counts:
                counts[contig] = {cls: 0 for cls in SOURCE_CLASS_ORDER}
            classes = h.source_classes if h.source_classes else ["Multiple"]
            for cls in classes:
                if cls in counts[contig]:
                    counts[contig][cls] += 1

        rows: List[Dict[str, Any]] = []
        for contig, class_counts in counts.items():
            total = sum(class_counts.values())
            if total == 0:
                continue
            row: Dict[str, Any] = {"contig": contig}
            row.update(class_counts)
            row["Total"] = total
            for cls in SOURCE_CLASS_ORDER:
                row[f"{cls} %"] = round(class_counts[cls] / total * 100, 1)

            dominant = max(class_counts.items(), key=lambda kv: kv[1])
            row["Dominant"] = dominant[0] if dominant[1] > 0 else NOT_ASSIGNED
            rows.append(row)

        rows.sort(key=lambda r: r["Total"], reverse=True)
        return rows

    @staticmethod
    def _is_key_gene(hit: BlastHit) -> bool:
        """Return True if the hit's gene name is in the S. aureus key gene list."""
        if not hit.gene_name:
            return False
        lower = hit.gene_name.lower()
        return any(k.lower() in lower for k in KEY_GENES)


# ---------------------------------------------------------------------------
# Per-sample report writers
# ---------------------------------------------------------------------------

def write_sample_json(result: SampleResult, output_dir: Path) -> Path:
    """Write a per-sample JSON report with full evidence details."""
    output_dir.mkdir(parents=True, exist_ok=True)
    payload = {
        "sample": result.sample,
        "input_file": result.input_file,
        "faa_file": result.faa_file,
        "mobileog_hits": result.total_hits,
        "category_counts": result.category_counts,
        "source_class_counts": result.source_class_counts,
        "detected_elements": result.elements_summary,
        "evidence": result.evidence,
        "warnings": result.warnings,
        "contig_source_class_composition": result.contig_purity,
        "interpretation": {
            "hit_definition": "mobileOG-associated protein hit",
            "element_definition": "Individual protein hits do not represent independent mobile genetic elements",
            "element_level_inference": False,
            "caveats": [
                "Multiple proteins from the same mobile element may generate multiple hits.",
                "Some mobileOG-db families, particularly replication/recombination/repair proteins, also occur in the bacterial chromosome.",
                "A single protein may match multiple source databases and be counted in more than one source class.",
                "Stronger element-level inference requires genomic context, co-localization of compatible MGE-associated genes, and, where appropriate, complementary databases or analyses."
            ],
            "functional_category_note": "Functional categories are inherited from mobileOG-db and indicate association with MGE-related protein families. A protein classified as 'Transfer' is not thereby asserted to mediate horizontal transfer itself."
        },
        "key_mge_associated_genes": [
            {
                "query": h.qseqid,
                "mobileog_id": h.mobileog_id,
                "gene_name": h.gene_name,
                "major_category": h.major_category,
                "minor_category": h.minor_category,
                "source_classes": h.source_classes,
                "contig": h.specific_contig,
                "orf_start": h.orf_start,
                "orf_stop": h.orf_stop,
                "strand": h.strand,
                "pident": round(h.pident, 2),
                "coverage": round(h.coverage, 2),
                "annotation": h.annotation,
            }
            for h in result.key_gene_hits
        ],
        "all_hits": [
            {
                "query": h.qseqid,
                "mobileog_id": h.mobileog_id,
                "gene_name": h.gene_name,
                "major_category": h.major_category,
                "minor_category": h.minor_category,
                "source_classes": h.source_classes,
                "source_databases": h.source_databases,
                "contig": h.specific_contig,
                "orf_start": h.orf_start,
                "orf_stop": h.orf_stop,
                "strand": h.strand,
                "partial": h.partial,
                "start_codon": h.start_codon,
                "rbs_motif": h.rbs_motif,
                "rbs_spacer": h.rbs_spacer,
                "gc_content": h.gc_content,
                "pident": round(h.pident, 2),
                "length": h.length,
                "qlen": h.qlen,
                "coverage": round(h.coverage, 2),
                "bitscore": h.bitscore,
                "annotation": h.annotation,
                "database": h.database,
                "evidence": h.evidence,
            }
            for h in result.raw_hits
        ],
    }
    path = output_dir / f"{result.sample}_mge.json"
    with open(path, "w") as fh:
        json.dump(payload, fh, indent=2)
    return path


def _hit_rows_html(hits: List[BlastHit]) -> str:
    """Render MGE hit rows as HTML table rows with coverage colouring."""
    if not hits:
        return "<tr><td colspan='12' style='text-align:center;color:#9ca3af;'>No hits above thresholds</td></tr>"
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
        start = h.orf_start if h.orf_start is not None else "-"
        stop = h.orf_stop if h.orf_stop is not None else "-"
        strand = h.strand or "-"
        rows.append(
            f"<tr class='{cls}'>"
            f"<td>{i + 1}</td>"
            f"<td><strong>{h.mobileog_id}</strong></td>"
            f"<td>{h.gene_name or '-'}</td>"
            f"<td>{h.specific_contig}</td>"
            f"<td>{start}</td>"
            f"<td>{stop}</td>"
            f"<td>{strand}</td>"
            f"<td>{h.major_category}</td>"
            f"<td>{h.minor_category or '-'}</td>"
            f"<td>{h.pident:.2f}</td>"
            f"<td>{h.length}/{h.qlen}</td>"
            f"<td>{h.coverage:.2f}</td>"
            f"</tr>"
        )
    return "\n".join(rows)


def _category_cards_html(result: SampleResult) -> str:
    """Render per-category hit cards for the sample report."""
    blocks = []
    counts = result.category_counts
    for category in CATEGORY_ORDER:
        n = counts.get(category, 0)
        if not n:
            continue
        label = CATEGORY_LABELS[category]
        hits = result.category_hits.get(category, [])
        names = sorted({h.gene_name for h in hits if h.gene_name and h.gene_name != "NA"})
        cards = "".join(f"<div class='element-card'>{name}</div>" for name in names)
        blocks.append(
            f"<h3>{label} ({n})</h3>"
            f"<div class='element-grid'>{cards}</div>"
        )
    if not blocks:
        return "<p style='color:#6b7280;'>No mobile genetic element proteins detected above thresholds.</p>"
    return "\n".join(blocks)


def _key_genes_table_html(result: SampleResult) -> str:
    """Render a compact table of key-gene hits for the sample report."""
    if not result.key_gene_hits:
        return "<p style='color:#6b7280;'>No S. aureus key MGE-associated genes detected.</p>"
    rows = []
    for h in sorted(result.key_gene_hits, key=lambda x: (x.gene_name or "", x.specific_contig)):
        start = h.orf_start if h.orf_start is not None else "-"
        stop = h.orf_stop if h.orf_stop is not None else "-"
        strand = h.strand or "-"
        rows.append(
            f"<tr>"
            f"<td><strong>{h.gene_name}</strong></td>"
            f"<td>{h.major_category}</td>"
            f"<td>{h.mobileog_id}</td>"
            f"<td>{h.specific_contig}</td>"
            f"<td>{start}..{stop} ({strand})</td>"
            f"<td>{h.pident:.2f}</td>"
            f"<td>{h.coverage:.2f}</td>"
            f"</tr>"
        )
    return f"""<table class="detailed-table">
<thead>
<tr><th>Gene</th><th>Category</th><th>mobileOG ID</th><th>Contig</th><th>Position</th><th>% Identity</th><th>Coverage %</th></tr>
</thead>
<tbody>
{"".join(rows)}
</tbody>
</table>"""


def _contig_purity_table_html(result: SampleResult) -> str:
    """Render the per-contig source-class composition table."""
    if not result.contig_purity:
        return "<p style='color:#6b7280;'>No contig source-class composition data available.</p>"

    header_cells = "".join(f"<th>{cls}</th>" for cls in SOURCE_CLASS_ORDER)
    body_rows = []
    for row in result.contig_purity:
        cells = f"<td><strong>{row['contig']}</strong></td>"
        for cls in SOURCE_CLASS_ORDER:
            count = row.get(cls, 0)
            pct = row.get(f"{cls} %", 0.0)
            cells += f"<td>{count} <span style='color:#6b7280;font-size:11px;'>({pct:.1f}%)</span></td>"
        cells += f"<td><strong>{row['Total']}</strong></td>"
        cells += f"<td>{row['Dominant']}</td>"
        body_rows.append(f"<tr>{cells}</tr>")

    return f"""<table class="detailed-table">
<thead>
<tr>
<th>Contig</th>
{header_cells}
<th>Total Hits</th>
<th>Dominant Source Class</th>
</tr>
</thead>
<tbody>
{"".join(body_rows)}
</tbody>
</table>
<p style="margin-top:10px;font-size:12px;color:#666;">
Counts reflect associations with source databases, not reconstructed elements. Multi-label hits are counted once per matched source class; the total across classes may exceed the number of distinct ORFs.
</p>"""


def write_sample_html(result: SampleResult, output_dir: Path) -> Path:
    """Write a per-sample deep-dive HTML report with all hits and evidence."""
    output_dir.mkdir(parents=True, exist_ok=True)

    badge_class = "status-detected" if result.total_hits > 0 else "status-none"
    badge_label = "MGE-associated protein evidence detected" if result.total_hits > 0 else NOT_ASSIGNED

    evidence_items = "".join(f"<li>{e}</li>" for e in result.evidence)
    if not evidence_items:
        evidence_items = "<li>No evidence recorded</li>"

    warning_html = ""
    if result.warnings:
        warning_html = "<div class='warning-box'><strong>⚠️ Warnings:</strong><ul>" + \
            "".join(f"<li>{w}</li>" for w in result.warnings) + "</ul></div>"

    category_html = _category_cards_html(result)
    key_genes_html = _key_genes_table_html(result)
    contig_purity_html = _contig_purity_table_html(result)

    all_hits_sorted = sorted(result.raw_hits, key=lambda h: h.bitscore, reverse=True)
    all_rows = _hit_rows_html(all_hits_sorted)

    counts = result.category_counts
    source_counts = result.source_class_counts
    integration_count = counts.get("integration/excision", 0)
    transfer_count = counts.get("transfer", 0)
    stability_count = counts.get("stability/transfer/defense", 0)
    phage_count = counts.get("phage", 0)
    replication_count = counts.get("replication/recombination/repair", 0)

    is_count = source_counts.get("Insertion sequences", 0)
    integ_count = source_counts.get("Integrative elements", 0)
    plasmid_count = source_counts.get("Plasmids", 0)
    bacterio_count = source_counts.get("Bacteriophages", 0)
    multi_count = source_counts.get("Multiple", 0)

    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>STAPHSCOPE - {result.sample} MGE Report</title>
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
    display:grid; grid-template-columns:repeat(auto-fit, minmax(180px, 1fr));
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
.status-detected {{ background:#7c3aed; color:white; }}
.status-none {{ background:#6b7280; color:white; }}
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
.element-grid {{
    display:grid; grid-template-columns:repeat(auto-fill, minmax(140px, 1fr));
    gap:10px; margin-top:10px;
}}
.element-card {{
    background: linear-gradient(135deg, #3b82f6 0%, #1e40af 100%);
    color:white; padding:10px; border-radius:8px;
    text-align:center; font-weight:bold; font-size:12px;
    word-break:break-word;
}}
.warning-box {{
    background:#fef3c7; border-left:4px solid #f59e0b;
    padding:15px; margin:15px 0; border-radius:6px; color:#92400e;
}}
.interpretation-box {{
    background:#eff6ff; border-left:4px solid #3b82f6;
    padding:15px; margin:15px 0; border-radius:6px; color:#1e3a8a;
    font-size:0.95em; line-height:1.6;
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
<div class="metric-card"><div class="metric-label">mobileOG-associated hits</div><div class="metric-value">{result.total_hits}</div></div>
<div class="metric-card"><div class="metric-label">Key MGE-associated genes</div><div class="metric-value">{len(result.key_gene_hits)}</div></div>
<div class="metric-card"><div class="metric-label">Functional categories</div><div class="metric-value">{len(result.category_hits)}</div></div>
</div>
</div>

<div class="report-section">
<h2>🧬 MGE-associated Protein Evidence</h2>
<div class="status-badge {badge_class}">{badge_label}</div>
<div class="interpretation-box">
<strong>Interpretation note:</strong> mobileOG-db identifies protein families associated with
mobile genetic element biology. Individual matches represent MGE-associated protein signatures
and should not be interpreted as independent mobile genetic elements. Multiple proteins from
the same element may generate multiple hits, and some families &mdash; particularly replication,
recombination, and repair proteins &mdash; also occur in the bacterial chromosome. Stronger
element-level inference requires genomic context, co-localization of compatible MGE-associated
genes, and, where appropriate, complementary databases or analyses.
</div>
<h3>Evidence</h3>
<ul style="margin-left:20px; line-height:1.8;">{evidence_items}</ul>
{warning_html}
</div>

<div class="report-section">
<h2>📈 Functional Category Breakdown</h2>
<div class="metrics-grid">
<div class="metric-card"><div class="metric-label">Integration</div><div class="metric-value">{integration_count}</div></div>
<div class="metric-card"><div class="metric-label">Transfer</div><div class="metric-value">{transfer_count}</div></div>
<div class="metric-card"><div class="metric-label">Stability</div><div class="metric-value">{stability_count}</div></div>
<div class="metric-card"><div class="metric-label">Phage</div><div class="metric-value">{phage_count}</div></div>
<div class="metric-card"><div class="metric-label">Replication</div><div class="metric-value">{replication_count}</div></div>
</div>
<p style="margin-top:15px;font-size:0.9em;color:#6b7280;">
Functional categories are inherited from mobileOG-db and indicate association with MGE-related
protein families. A protein classified as &ldquo;Transfer&rdquo; is not thereby asserted to
mediate horizontal transfer itself.
</p>
</div>

<div class="report-section">
<h2>🧩 Source-class Association Breakdown</h2>
<div class="metrics-grid">
<div class="metric-card"><div class="metric-label">IS-associated</div><div class="metric-value">{is_count}</div></div>
<div class="metric-card"><div class="metric-label">ICE-associated</div><div class="metric-value">{integ_count}</div></div>
<div class="metric-card"><div class="metric-label">Plasmid-associated</div><div class="metric-value">{plasmid_count}</div></div>
<div class="metric-card"><div class="metric-label">Phage-associated</div><div class="metric-value">{bacterio_count}</div></div>
<div class="metric-card"><div class="metric-label">Multi-source</div><div class="metric-value">{multi_count}</div></div>
</div>
<p style="margin-top:15px;font-size:0.9em;color:#6b7280;">
Counts reflect associations with source databases, not reconstructed elements. A single protein
may match multiple source classes; the total across classes may exceed the number of distinct ORFs.
</p>
</div>

<div class="report-section">
<h2>🧩 Contig Source-class Composition</h2>
{contig_purity_html}
</div>

<div class="report-section">
<h2>🎯 Key MGE-associated Genes</h2>
{key_genes_html}
</div>

<div class="report-section">
<h2>🧪 Detected Protein Families by Functional Category</h2>
{category_html}
</div>

<div class="report-section">
<h2>🔍 All mobileOG-db Hits</h2>
<div style="max-height:600px; overflow-y:auto;">
<table class="detailed-table">
<thead>
<tr>
<th>Rank</th><th>mobileOG ID</th><th>Gene</th><th>Contig</th>
<th>Start</th><th>Stop</th><th>Strand</th>
<th>Major Category</th><th>Minor Category</th>
<th>% Identity</th><th>Align/Qlen</th><th>Coverage %</th>
</tr>
</thead>
<tbody>
{all_rows}
</tbody>
</table>
</div>
<p style="margin-top:10px; font-size:12px; color:#666;">
● Top hit in RED &nbsp;|&nbsp; <span style="color:#ea580c;">● Coverage ≥90% in ORANGE</span> &nbsp;|&nbsp; <span style="color:#d97706;">● Coverage ≥70% in YELLOW</span>
</p>
</div>

<div class="footer">
<p><strong>STAPHSCOPE</strong> - Mobile Genetic Element Detection</p>
<p class="timestamp">Generated: {time.strftime('%Y-%m-%d %H:%M:%S')}</p>
<div class="authorship">
<p><strong>Technical Support &amp; Inquiries:</strong></p>
<p>Integrated by: Brown Beckley</p>
<p>GitHub: <a href="https://github.com/bbeckley-hub" style="color:#fbbf24;">bbeckley-hub</a></p>
<p>Email: <a href="mailto:brownbeckley94@gmail.com" style="color:#fbbf24;">brownbeckley94@gmail.com</a></p>
<p>Affiliation: University of Ghana Medical School - Department of Medical Biochemistry</p>
<p style="margin-top:10px;font-size:11px;opacity:0.8;">
Powered by mobileOG-db (Brown et al.), Prodigal (Hyatt et al.), and DIAMOND (Buchfink et al.)<br>
Parsing logic adapted from mobileOG-pl (Mullet &amp; Brown)
</p>
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
    path = output_dir / f"{result.sample}_mge.html"
    with open(path, "w", encoding="utf-8") as fh:
        fh.write(html)
    return path


# ---------------------------------------------------------------------------
# Top-level summary writers
# ---------------------------------------------------------------------------

def write_summary_tsv(results: List[SampleResult], output_path: Path) -> None:
    """Write a batch TSV summary with one row per sample."""
    columns = [
        "sample", "total_hits", "key_genes",
        "functional_integration", "functional_transfer", "functional_stability",
        "functional_phage", "functional_replication",
        "source_IS", "source_integrative", "source_plasmid",
        "source_bacteriophage", "source_multiple",
        "detected_elements",
    ]
    with open(output_path, "w") as fh:
        fh.write("\t".join(columns) + "\n")
        for r in results:
            counts = r.category_counts
            source = r.source_class_counts
            row = [
                r.sample,
                r.total_hits,
                len(r.key_gene_hits),
                counts.get("integration/excision", 0),
                counts.get("transfer", 0),
                counts.get("stability/transfer/defense", 0),
                counts.get("phage", 0),
                counts.get("replication/recombination/repair", 0),
                source.get("Insertion sequences", 0),
                source.get("Integrative elements", 0),
                source.get("Plasmids", 0),
                source.get("Bacteriophages", 0),
                source.get("Multiple", 0),
                r.elements_summary,
            ]
            fh.write("\t".join(str(v) for v in row) + "\n")


def write_summary_json(results: List[SampleResult], output_path: Path) -> None:
    """Write a combined JSON summary for all samples."""
    payload = {
        "metadata": {
            "tool": "MGEFinder (StaphScope)",
            "database": "mobileOG-db",
            "interpretation": {
                "hit_definition": "mobileOG-associated protein hit",
                "element_definition": "Individual protein hits do not represent independent mobile genetic elements",
                "element_level_inference": False,
                "caveats": [
                    "Multiple proteins from the same mobile element may generate multiple hits.",
                    "Some mobileOG-db families, particularly replication/recombination/repair proteins, also occur in the bacterial chromosome.",
                    "A single protein may match multiple source databases and be counted in more than one source class.",
                    "Stronger element-level inference requires genomic context, co-localization of compatible MGE-associated genes, and, where appropriate, complementary databases or analyses."
                ],
                "functional_category_note": "Functional categories are inherited from mobileOG-db and indicate association with MGE-related protein families. A protein classified as 'Transfer' is not thereby asserted to mediate horizontal transfer itself."
            }
        },
        "samples": []
    }

    for r in results:
        payload["samples"].append({
            "sample": r.sample,
            "mobileog_hits": r.total_hits,
            "category_counts": r.category_counts,
            "source_class_counts": r.source_class_counts,
            "key_mge_associated_genes": sorted({h.gene_name for h in r.key_gene_hits if h.gene_name}),
            "detected_elements": r.elements_summary,
            "evidence": r.evidence,
            "warnings": r.warnings,
            "contig_source_class_composition": r.contig_purity,
        })

    with open(output_path, "w") as fh:
        json.dump(payload, fh, indent=2)


def write_summary_html(results: List[SampleResult], output_path: Path, elapsed: float) -> None:
    """Write a StaphScope-styled batch summary HTML report."""
    total = len(results)
    with_mges = sum(1 for r in results if r.total_hits > 0)
    without_mges = total - with_mges
    total_hits = sum(r.total_hits for r in results)
    total_key = sum(len(r.key_gene_hits) for r in results)

    rows = []
    for r in results:
        badge_class = "status-detected" if r.total_hits > 0 else "status-none"
        badge_label = str(r.total_hits) if r.total_hits > 0 else NOT_ASSIGNED
        c = r.category_counts
        s = r.source_class_counts
        rows.append(
            f"<tr>"
            f"<td><strong>{r.sample}</strong></td>"
            f"<td><span class='{badge_class}'>{badge_label}</span></td>"
            f"<td>{c.get('integration/excision', 0)}</td>"
            f"<td>{c.get('transfer', 0)}</td>"
            f"<td>{c.get('stability/transfer/defense', 0)}</td>"
            f"<td>{c.get('phage', 0)}</td>"
            f"<td>{c.get('replication/recombination/repair', 0)}</td>"
            f"<td>{s.get('Insertion sequences', 0)}</td>"
            f"<td>{s.get('Integrative elements', 0)}</td>"
            f"<td>{s.get('Plasmids', 0)}</td>"
            f"<td>{s.get('Bacteriophages', 0)}</td>"
            f"<td>{len(r.key_gene_hits)}</td>"
            f"</tr>"
        )
    table_rows = "\n".join(rows)

    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>STAPHSCOPE - MGE Batch Summary</title>
<style>
* {{ margin:0; padding:0; box-sizing:border-box; }}
body {{
    background: linear-gradient(135deg, #1e3c72 0%, #2a5298 50%, #7e22ce 100%);
    font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif;
    color:#ffffff; padding:20px; min-height:100vh;
}}
.container {{ max-width:1900px; margin:0 auto; }}
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
    width:100%; border-collapse:collapse; margin-top:15px; font-size:12px;
}}
.results-table th {{
    background:#1e40af; color:white; padding:10px; text-align:left;
    font-weight:bold; position:sticky; top:0;
}}
.results-table td {{
    padding:10px; border-bottom:1px solid #e5e7eb; vertical-align:top;
}}
.results-table tr:hover {{ background:#f3f4f6; }}
.status-detected {{
    background:#7c3aed; color:white; padding:4px 8px;
    border-radius:4px; font-weight:bold; font-size:11px;
}}
.status-none {{
    background:#6b7280; color:white; padding:4px 8px;
    border-radius:4px; font-weight:bold; font-size:11px;
}}
.interpretation-box {{
    background:#eff6ff; border-left:4px solid #3b82f6;
    padding:15px; margin:20px 0; border-radius:6px; color:#1e3a8a;
    font-size:0.95em; line-height:1.6;
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
<h2>📊 MGE-associated Protein Summary</h2>
<div class="metrics-grid">
<div class="metric-card"><div class="metric-label">Total Samples</div><div class="metric-value">{total}</div></div>
<div class="metric-card"><div class="metric-label">With mobileOG hits</div><div class="metric-value">{with_mges}</div></div>
<div class="metric-card"><div class="metric-label">Without mobileOG hits</div><div class="metric-value">{without_mges}</div></div>
<div class="metric-card"><div class="metric-label">Total mobileOG hits</div><div class="metric-value">{total_hits}</div></div>
<div class="metric-card"><div class="metric-label">Key MGE-associated genes</div><div class="metric-value">{total_key}</div></div>
<div class="metric-card"><div class="metric-label">Runtime</div><div class="metric-value" style="font-size:20px;">{elapsed:.1f}s</div></div>
</div>
<div class="interpretation-box">
<strong>Interpretation note:</strong> mobileOG-db identifies protein families associated with
mobile genetic element biology. Individual matches represent MGE-associated protein signatures
and should not be interpreted as independent mobile genetic elements. Multiple proteins from
the same element may generate multiple hits, and some families &mdash; particularly replication,
recombination, and repair proteins &mdash; also occur in the bacterial chromosome. Stronger
element-level inference requires genomic context and co-localization of compatible
MGE-associated genes.
</div>
</div>

<div class="report-section">
<h2>🧬 Detailed Results</h2>
<div style="max-height:700px; overflow-y:auto;">
<table class="results-table">
<thead>
<tr>
<th>Sample</th>
<th>mobileOG Hits</th>
<th>Integr.</th>
<th>Transfer</th>
<th>Stab.</th>
<th>Phage</th>
<th>Repl.</th>
<th>IS-assoc.</th>
<th>ICE-assoc.</th>
<th>Plasmid-assoc.</th>
<th>Phage-assoc.</th>
<th>Key MGE-assoc.</th>
</tr>
</thead>
<tbody>
{table_rows}
</tbody>
</table>
</div>
<p style="margin-top:15px;font-size:0.9em;color:#6b7280;">
Counts reflect associations with mobileOG-db families and source databases, not reconstructed
elements. A single protein may match multiple source classes; the total across classes may
exceed the number of distinct ORFs.
</p>
</div>

<div class="footer">
<p><strong>STAPHSCOPE</strong> - Mobile Genetic Element Detection</p>
<p class="timestamp">Generated: {time.strftime('%Y-%m-%d %H:%M:%S')}</p>
<p style="margin-top:10px; font-size:12px;">Powered by mobileOG-db, Prodigal, and DIAMOND</p>
<div class="authorship">
<p><strong>Technical Support &amp; Inquiries:</strong></p>
<p>Integrated by: Brown Beckley</p>
<p>GitHub: <a href="https://github.com/bbeckley-hub" style="color:#fbbf24;">bbeckley-hub</a></p>
<p>Email: <a href="mailto:brownbeckley94@gmail.com" style="color:#fbbf24;">brownbeckley94@gmail.com</a></p>
<p>Affiliation: University of Ghana Medical School - Department of Medical Biochemistry</p>
<p style="margin-top:10px;font-size:11px;opacity:0.8;">
Powered by mobileOG-db (Brown et al.), Prodigal (Hyatt et al.), and DIAMOND (Buchfink et al.)<br>
Parsing logic adapted from mobileOG-pl (Mullet &amp; Brown)
</p>
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
# Batch orchestration
# ---------------------------------------------------------------------------

def _resolve_batch_plan(total_cpus: int, n_jobs: int, ram_gb: float, user_workers: Optional[int] = None):
    """
    Compute (workers, threads_per_worker) so that workers * threads_per_worker
    never exceeds total_cpus. Uses RAM to cap concurrent workers conservatively.
    """
    workers = max(1, min(total_cpus, n_jobs, int(ram_gb / 1.5)))
    if user_workers is not None:
        workers = max(1, min(user_workers, total_cpus, n_jobs))
    threads_per_worker = max(1, total_cpus // workers)
    return workers, threads_per_worker


def _process_one(typer: MGETyper, fasta_path: Path, annotations_dir: Path, diamond_threads: int) -> SampleResult:
    """Worker function: run the typer on a single assembly."""
    return typer.run(fasta_path, annotations_dir, diamond_threads=diamond_threads)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Mobile genetic element detection via Prodigal + DIAMOND + mobileOG-db",
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
        help="Directory containing the mobileOG-db .dmnd and metadata CSV",
    )
    parser.add_argument(
        "-p", "--prefix", default="staphscope_mge",
        help="Prefix for top-level summary files",
    )
    parser.add_argument(
        "--annotations-dir", default=None,
        help="Directory for cached Prodigal .faa files (default: <output-dir>/annotations)",
    )
    parser.add_argument(
        "--min-id", type=float, default=80.0,
        help="Minimum percent identity for a DIAMOND hit",
    )
    parser.add_argument(
        "--min-cov", type=float, default=60.0,
        help="Minimum query coverage for a DIAMOND hit",
    )
    parser.add_argument(
        "--max-evalue", type=float, default=1e-20,
        help="Maximum e-value for a DIAMOND hit",
    )
    parser.add_argument(
        "--cpus", "-c", type=int, default=None,
        help="Total CPU core budget (default: auto-detect optimal for this system)",
    )
    parser.add_argument(
        "--workers", "-w", type=int, default=None,
        help="Override number of concurrent workers (default: auto-computed from CPU and RAM)",
    )
    args = parser.parse_args()

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    annotations_dir = Path(args.annotations_dir) if args.annotations_dir else output_dir / "annotations"
    annotations_dir.mkdir(parents=True, exist_ok=True)

    start = time.time()
    try:
        typer = MGETyper(
            Path(args.database_dir),
            min_id=args.min_id,
            min_cov=args.min_cov,
            max_evalue=args.max_evalue,
            cpus=args.cpus,
        )
    except FileNotFoundError as e:
        print(f"Error: {e}")
        sys.exit(1)

    if not typer.check_db():
        print(f"Error: mobileOG-db .dmnd or metadata missing in {args.database_dir}")
        sys.exit(1)

    print(f"Using DIAMOND database: {typer.dmnd.name}")
    print(f"Using metadata: {typer.metadata_csv.name}")
    print(f"Loaded {len(typer.metadata)} mobileOG entries")
    print(f"System resources: {typer.cpus} CPU cores | {typer.available_ram:.1f} GB available RAM")
    print()

    fasta_paths = [Path(f) for f in args.input if Path(f).exists()]
    missing = [f for f in args.input if not Path(f).exists()]
    for f in missing:
        print(f"Skipping missing file: {f}")

    if not fasta_paths:
        print("No samples processed.")
        sys.exit(1)

    workers, threads_per_worker = _resolve_batch_plan(
        total_cpus=typer.cpus,
        n_jobs=len(fasta_paths),
        ram_gb=typer.available_ram,
        user_workers=args.workers,
    )

    if len(fasta_paths) == 1 or workers == 1:
        print(f"Batch plan: sequential (1 worker × {threads_per_worker} thread(s) per worker)")
    else:
        print(f"Batch plan: {workers} concurrent worker(s) × {threads_per_worker} thread(s) per worker")
    print(f"Total thread budget: {workers * threads_per_worker} of {typer.cpus} available cores")
    print()

    results: List[SampleResult] = []

    if workers == 1 or len(fasta_paths) == 1:
        for fasta in fasta_paths:
            print(f"Processing {fasta.name} ...")
            result = typer.run(fasta, annotations_dir, diamond_threads=threads_per_worker)
            results.append(result)
            print(f"  {result.sample}: {result.total_hits} mobileOG hit(s) | {len(result.key_gene_hits)} key MGE-associated gene(s)")
            sample_out = output_dir / result.sample
            write_sample_json(result, sample_out)
            write_sample_html(result, sample_out)
    else:
        with ThreadPoolExecutor(max_workers=workers) as executor:
            futures = {
                executor.submit(_process_one, typer, fasta, annotations_dir, threads_per_worker): fasta
                for fasta in fasta_paths
            }
            for future in as_completed(futures):
                fasta = futures[future]
                try:
                    result = future.result()
                    results.append(result)
                    print(f"  {result.sample}: {result.total_hits} mobileOG hit(s) | {len(result.key_gene_hits)} key MGE-associated gene(s)")
                    sample_out = output_dir / result.sample
                    write_sample_json(result, sample_out)
                    write_sample_html(result, sample_out)
                except Exception as e:
                    print(f"  Failed: {fasta.name} - {e}")

    if not results:
        print("No samples processed successfully.")
        sys.exit(1)

    write_summary_tsv(results, output_dir / f"{args.prefix}_summary.tsv")
    write_summary_json(results, output_dir / f"{args.prefix}_summary.json")
    write_summary_html(results, output_dir / f"{args.prefix}_summary.html", time.time() - start)

    elapsed = time.time() - start
    print(f"\nResults written to {output_dir}")
    print(f"  {args.prefix}_summary.tsv")
    print(f"  {args.prefix}_summary.json")
    print(f"  {args.prefix}_summary.html")
    print(f"  <sample>/<sample>_mge.html  (per-sample deep-dive)")
    print(f"  <sample>/<sample>_mge.json  (per-sample JSON)")
    print(f"  annotations/                (cached Prodigal .faa files)")
    print(f"\nCompleted {len(results)} sample(s) in {elapsed:.1f} seconds")
    print(f"Processing mode: {'parallel (' + str(workers) + ' workers)' if workers > 1 else 'sequential'}")


if __name__ == "__main__":
    main()
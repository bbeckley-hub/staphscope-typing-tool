#!/usr/bin/env python3
"""
Staphscope core functionality
"""

from __future__ import annotations
import csv
import datetime
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Dict, List, Optional, Tuple

BANNER = (
    "\n=== Staphscope Typing Tool — v0.1.0 (2025-08-15) ===\n"
    "University of Ghana / K.N.U.S.T\n"
    "Author: Beckley Brown <brownbeckley94@gmail.com>\n"
)

FOOTER = (
    "\n------\n"
    "Done with MLST + spa + SCCmec typing. Enjoy your downstream analysis.\n"
    "Check GitHub for my Bioinformatics tools in the future: bbeckley-hub\n"
    "------\n"
)

# Default locations for SCCmecFinder
HOME = Path.home()
SCCMECFINDER_DIR = HOME / "sccmecfinder"
SCCMECFINDER_SCRIPT = SCCMECFINDER_DIR / "SCCmecFinder_v4.py"
SCCMEC_DB_DIR = SCCMECFINDER_DIR / "database"
SCCMEC_SCRIPT_DIR = SCCMECFINDER_DIR / "script_dir"

def log(msg: str):
    print(msg, flush=True)

def run_cmd(cmd: List[str], cwd: Optional[Path] = None, check: bool = True) -> subprocess.CompletedProcess:
    return subprocess.run(cmd, cwd=str(cwd) if cwd else None, capture_output=True, text=True, check=check)

def which_or_none(name: str) -> Optional[str]:
    return shutil.which(name)

def ensure_dir(p: Path):
    p.mkdir(parents=True, exist_ok=True)

def check_environment() -> Dict[str, Optional[str]]:
    tools = {
        "mlst": which_or_none("mlst"),
        "spatyper": which_or_none("spatyper") or which_or_none("spaTyper") or which_or_none("spaTyper.py"),
        "blastn": which_or_none("blastn"),
        "makeblastdb": which_or_none("makeblastdb"),
        "python3": which_or_none("python3")
    }
    return tools

def print_environment_report(tools: Dict[str, Optional[str]]):
    log("\n[Environment Check]")
    for k, v in tools.items():
        if v:
            log(f"  - {k:12s}: {v}")
        else:
            log(f"  - {k:12s}: NOT FOUND (please install / add to PATH)")
    
    log("\n[SCCmecFinder Check]")
    if SCCMECFINDER_SCRIPT.exists():
        log(f"  - SCCmecFinder_v4.py: Found at {SCCMECFINDER_SCRIPT}")
    else:
        log(f"  - SCCmecFinder_v4.py: NOT FOUND (expected at {SCCMECFINDER_SCRIPT})")
    
    if SCCMEC_DB_DIR.exists():
        log(f"  - SCCmec Database:   Found at {SCCMEC_DB_DIR}")
    else:
        log(f"  - SCCmec Database:   NOT FOUND (expected at {SCCMEC_DB_DIR})")
    
    if SCCMEC_SCRIPT_DIR.exists():
        log(f"  - SCCmec Script Dir: Found at {SCCMEC_SCRIPT_DIR}")
    else:
        log(f"  - SCCmec Script Dir: NOT FOUND (expected at {SCCMEC_SCRIPT_DIR})")

def try_update_databases(tools: Dict[str, Optional[str]]):
    """Attempt to update typing databases if tools support it."""
    if tools.get("mlst"):
        try:
            log("[Update] Running: mlst --update")
            run_cmd([tools["mlst"], "--update"], check=False)
        except Exception as e:
            log(f"[Update] mlst --update failed or skipped: {e}")

    # spaTyper update: not all builds support auto-update; add your command if available
    if tools.get("spatyper"):
        # Placeholder: if your spaTyper supports a DB update flag, add it here
        log("[Update] spaTyper auto-update: no standard command detected; skip (configure manually if needed).")

def run_mlst(sample: Path, tools: Dict[str, Optional[str]]) -> Tuple[str, str, Dict[str, str]]:
    """Run MLST and return (scheme, ST, loci dict)."""
    mlst_bin = tools.get("mlst")
    if not mlst_bin:
        return ("NA", "NA", {})
    try:
        cp = run_cmd([mlst_bin, str(sample)], check=False)
        if cp.returncode != 0:
            return ("NA", "NA", {})
        # Expected format: FILE\tscheme\tST\tabcZ(1)\t…
        line = cp.stdout.strip().splitlines()[-1]
        parts = line.split("\t")
        if len(parts) < 3:
            return ("NA", "NA", {})
        scheme, st = parts[1], parts[2]
        loci = {}
        for p in parts[3:]:
            if "(" in p and ")" in p:
                locus = p.split("(")[0]
                allele = p.split("(")[-1].rstrip(")")
                loci[locus] = allele
        return (scheme, st, loci)
    except Exception:
        return ("NA", "NA", {})

def parse_spatyper_results(tsv_path: Path) -> Tuple[str, Optional[float]]:
    """Parse spaTyper results TSV."""
    spa_type = "NA"
    score = None
    if not tsv_path.exists():
        return spa_type, score
    
    for encoding in ['utf-8', 'latin-1']:
        try:
            with tsv_path.open('r', encoding=encoding) as f:
                try:
                    reader = csv.DictReader(f, delimiter="\t")
                    for row in reader:
                        for key in ("spa_type", "spaType", "spa", "type", "Spa Type", "SpaType", "spa_type_result"):
                            if key in row and row[key]:
                                spa_type = row[key].strip()
                                if spa_type.lower().startswith("type:"):
                                    spa_type = spa_type.split(":")[1].strip()
                                if spa_type.lower().startswith("t"):
                                    spa_type = spa_type.lower()
                                break
                        if spa_type != "NA":
                            break
                        for key in ("score", "confidence", "identity", "Score", "Identity"):
                            if key in row and row[key]:
                                try:
                                    score = float(row[key])
                                except ValueError:
                                    pass
                    if spa_type != "NA":
                        return spa_type, score
                except csv.Error:
                    f.seek(0)
                    content = f.read()
                    match = re.search(r'\b[tT](\d{2,4})\b', content)
                    if match:
                        spa_type = f"t{match.group(1)}"
                        return spa_type, score
            break
        except UnicodeDecodeError:
            continue
            
    return spa_type, score

def run_spa(sample: Path, tools: Dict[str, Optional[str]], tmpdir: Path) -> Tuple[str, Optional[float]]:
    """Run spaTyper and return spa_type and score."""
    spa_bin = tools.get("spatyper")
    if not spa_bin:
        return ("NA", None)
    
    outdir = tmpdir / f"spa_{sample.stem}"
    ensure_dir(outdir)
    
    cmds = [
        [spa_bin, "-f", str(sample), "-o", str(outdir)],
        [spa_bin, str(sample), "-o", str(outdir)],
        [spa_bin, "-i", str(sample), "-o", str(outdir)],
        [spa_bin, str(sample)],
        [spa_bin, "-f", str(sample)]
    ]
    
    for cmd in cmds:
        try:
            log(f"  [spa] Trying command: {' '.join(cmd)}")
            cp = run_cmd(cmd, check=False)
            
            result_files = list(outdir.glob("*"))
            for result_file in result_files:
                if result_file.suffix in ['.tsv', '.txt', '.csv']:
                    spa_type, score = parse_spatyper_results(result_file)
                    if spa_type != "NA":
                        return spa_type, score
            
            output = cp.stdout + "\n" + cp.stderr
            match = re.search(r'\b[tT](\d{2,4})\b', output)
            if match:
                spa_type = f"t{match.group(1)}"
                return spa_type, None
                
            match = re.search(r'SPA type[:\s]*([tT]\d{2,4})\b', output)
            if match:
                spa_type = match.group(1).lower()
                return spa_type, None
        except Exception as e:
            log(f"  [spa] Command failed: {e}")
            continue

    return ("NA", None)

def run_sccmec_cge(sample: Path, tools: Dict[str, Optional[str]], tmpdir: Path) -> Tuple[str, List[Dict[str, str]]]:
    """Run CGE SCCmecFinder and return SCCmec type."""
    python3 = tools.get("python3")
    if not python3:
        return "NA", []
    
    if not SCCMECFINDER_SCRIPT.exists() or not SCCMEC_DB_DIR.exists():
        log(f"  [SCCmec] SCCmecFinder not found at {SCCMECFINDER_SCRIPT} or database not found at {SCCMEC_DB_DIR}")
        return "NA", []
    
    out_dir = tmpdir / "sccmec_results"
    ensure_dir(out_dir)
    
    db_input = out_dir / "db_input.fna"
    km_input = out_dir / "km_input.fna"
    shutil.copy(sample, db_input)
    shutil.copy(sample, km_input)
    
    cmd = [
        python3, str(SCCMECFINDER_SCRIPT),
        "-iDb", str(db_input),
        "-iKm", str(km_input),
        "-k", "90",
        "-l", "60",
        "-o", "SCCmecFinder_results.txt",
        "-d", str(out_dir),
        "-db_dir", str(SCCMEC_DB_DIR),
        "-sc_dir", str(SCCMEC_SCRIPT_DIR),
        "-db_choice", "reference"
    ]
    
    try:
        log(f"  [SCCmec] Running command: {' '.join(cmd)}")
        run_cmd(cmd, cwd=out_dir, check=True)
        
        result_file = out_dir / "SCCmecFinder_results.txt"
        if not result_file.exists():
            log(f"  [SCCmec] Result file not found: {result_file}")
            return "NA", []
        
        scc_type = "NA"
        content = result_file.read_text()
        
        gene_match = re.search(r'Prediction based on genes:\s*(.+)', content)
        if gene_match:
            scc_type = gene_match.group(1).strip()
        else:
            if "No SCCmec element was detected" in content:
                scc_type = "NA"
            else:
                homol_match = re.search(r'Prediction based on homology:\s*(.+)', content)
                if homol_match:
                    scc_type = homol_match.group(1).strip()
        
        if scc_type == "Not typable" or scc_type == "NA":
            scc_type = "NA"

        log(f"  [SCCmec] Final type: {scc_type}")
        return scc_type, []
    
    except subprocess.CalledProcessError as e:
        log(f"  [SCCmec] SCCmecFinder failed: {e}")
        log(f"  [SCCmec] Error output: {e.stderr}")
        return "NA", []
    except Exception as e:
        log(f"  [SCCmec] Error: {e}")
        return "NA", []

def process_sample(sample: Path, tools: Dict[str, Optional[str]], outdir: Path) -> Dict[str, object]:
    log(f"[Run] {sample.name}")
    ensure_dir(outdir)
    tmpdir = Path(tempfile.mkdtemp(prefix=f"staphscope_{sample.stem}_"))

    # MLST
    log("  Running MLST...")
    scheme, st, loci = run_mlst(sample, tools)

    # spa
    log("  Running spa typing...")
    spa_type, spa_score = run_spa(sample, tools, tmpdir)

    # SCCmec
    log("  Running SCCmec typing...")
    sccmec_type, scc_hits = run_sccmec_cge(sample, tools, tmpdir)

    # Clean up temporary directory
    try:
        shutil.rmtree(tmpdir)
    except Exception as e:
        log(f"  Warning: Failed to clean up temp directory: {e}")

    # Write per-sample JSON detail
    detail = {
        "sample": sample.name,
        "mlst": {"scheme": scheme, "ST": st, "loci": loci},
        "spa": {"spa_type": spa_type, "score": spa_score},
        "sccmec": {"type": sccmec_type, "hits": scc_hits},
    }
    with (outdir / f"{sample.stem}.staphscope.json").open("w") as f:
        json.dump(detail, f, indent=2)

    return {
        "sample": sample.name,
        "mlst_scheme": scheme,
        "mlst_ST": st,
        "spa_type": spa_type,
        "sccmec_type": sccmec_type,
    }

#!/usr/bin/env python3
"""
Staphscope core functionality — bundled MLST + spaTyper + SCCmecFinder
"""

from __future__ import annotations
import csv
import json
import shutil
import subprocess
import sys
import tempfile
import os
from pathlib import Path
from typing import Dict, List, Optional, Tuple

BANNER = "\n=== Staphscope Typing Tool — v0.2.0 (2025-08-20) ===\nUniversity of Ghana / K.N.U.S.T\nU.G.M.S - Department of Medical Biochemistry\nAuthor: Beckley Brown <brownbeckley94@gmail.com>\n"
FOOTER = """
------
Done with MLST + spa + SCCmec typing. Enjoy your downstream analysis.

If you use Staphscope Typing Tool in your research, please cite it as:

Brown, B. (2025). Staphscope Typing Tool: Unified MLST + spa + SCCmec typing for *Staphylococcus aureus*. 
K.N.U.S.T/U.G.M.S - Department of Medical Biochemistry, University of Ghana. 
GitHub repository: https://github.com/bbeckley-hub/staphscope-typing-tool
------
"""

# Base dir = staphscope package
BASE_DIR = Path(__file__).resolve().parent

# ===== Bundled Tools =====
MLST_BIN = BASE_DIR / "tools" / "mlst" / "bin" / "mlst"
SPATYPER_BIN = BASE_DIR / "tools" / "spatyper" / "spa_typing" / "main" / "spaTyper"
SCCMECFINDER_SCRIPT = BASE_DIR / "tools" / "sccmecfinder" / "SCCmecFinder_v4.py"
SCCMEC_DB_DIR = BASE_DIR / "tools" / "sccmecfinder" / "database"
SCCMEC_SCRIPT_DIR = BASE_DIR / "tools" / "sccmecfinder" / "script_dir"

# ===== Database Paths =====
SPA_REPEATS_FILE = BASE_DIR / "tools" / "spatyper" / "sparepeats.fasta"
SPA_TYPES_FILE = BASE_DIR / "tools" / "spatyper" / "spatypes.txt"

# ===== Utilities =====
def log(msg: str):
    print(msg, flush=True)

def run_cmd(cmd: List[str], cwd: Optional[Path] = None, env: Optional[Dict] = None) -> subprocess.CompletedProcess:
    return subprocess.run(cmd, cwd=str(cwd) if cwd else None, env=env, capture_output=True, text=True, check=False)

def ensure_dir(p: Path):
    p.mkdir(parents=True, exist_ok=True)

# ===== Environment =====
def check_environment() -> Dict[str, Optional[str]]:
    tools = {
        "mlst": str(MLST_BIN) if MLST_BIN.exists() else None,
        "spatyper": str(SPATYPER_BIN) if SPATYPER_BIN.exists() else None,
        "blastn": shutil.which("blastn"),
        "makeblastdb": shutil.which("makeblastdb"),
        "python3": shutil.which("python3")
    }
    return tools

def print_environment_report(tools: Dict[str, Optional[str]]):
    log("\n[Environment Check]")
    for k, v in tools.items():
        log(f"  - {k:12s}: {v or 'NOT FOUND'}")
    log("\n[SCCmecFinder Check]")
    log(f"  - SCCmecFinder_v4.py: {'Found' if SCCMECFINDER_SCRIPT.exists() else 'NOT FOUND'}")
    log(f"  - SCCmec Database:   {'Found' if SCCMEC_DB_DIR.exists() else 'NOT FOUND'}")
    log(f"  - SCCmec Script Dir: {'Found' if SCCMEC_SCRIPT_DIR.exists() else 'NOT FOUND'}")

# ===== MLST =====
def run_mlst(sample: Path, mlst_bin: Optional[str]) -> Tuple[str, str]:
    if not mlst_bin:
        return "NA", "NA"
    
    env = os.environ.copy()
    env['PERL5LIB'] = str(BASE_DIR / "tools" / "mlst" / "perl5")
    
    cp = run_cmd([mlst_bin, str(sample)], env=env)
    if cp.returncode != 0:
        return "NA", "NA"
    
    lines = cp.stdout.strip().splitlines()
    result_line = None
    for line in lines:
        if line.startswith(str(sample)):
            result_line = line
            break
    
    if not result_line:
        return "NA", "NA"
    
    parts = result_line.split("\t")
    if len(parts) < 3:
        return "NA", "NA"
    
    return parts[1], parts[2]

# ===== spaTyper =====
def run_spa(sample: Path, spatyper_bin: Optional[str], tmpdir: Path) -> str:
    if not spatyper_bin or not Path(spatyper_bin).exists() or not SPA_REPEATS_FILE.exists() or not SPA_TYPES_FILE.exists():
        return "NA"
    
    outfile = tmpdir / f"spa_result_{sample.stem}.txt"
    cmd = ["python3", spatyper_bin, "-f", str(sample), "--output", str(outfile)]
    
    env = os.environ.copy()
    env["PYTHONPATH"] = str(BASE_DIR / "tools" / "spatyper" / "spa_typing")
    
    cp = run_cmd(cmd, env=env)
    if cp.returncode != 0:
        return "NA"

    if outfile.exists():
        with outfile.open() as f:
            lines = [line.strip() for line in f if line.strip()]
            if len(lines) > 1 and "\t" in lines[1]:
                return lines[1].split("\t")[-1]
            elif len(lines) == 1 and "\t" in lines[0]:
                return lines[0].split("\t")[-1]
    return "NA"

# ===== SCCmecFinder =====
def run_sccmec_cge(sample: Path, python3: Optional[str], tmpdir: Path) -> str:
    if not python3 or not SCCMECFINDER_SCRIPT.exists():
        return "NA"
    
    outdir = tmpdir / "sccmec"
    ensure_dir(outdir)
    shutil.copy(sample, outdir / "db_input.fna")
    shutil.copy(sample, outdir / "km_input.fna")
    
    cmd = [
        python3, str(SCCMECFINDER_SCRIPT),
        "-iDb", str(outdir / "db_input.fna"),
        "-iKm", str(outdir / "km_input.fna"),
        "-k", "90", "-l", "60",
        "-o", "SCCmecFinder_results.txt",
        "-d", str(outdir),
        "-db_dir", str(SCCMEC_DB_DIR),
        "-sc_dir", str(SCCMEC_SCRIPT_DIR),
        "-db_choice", "reference"
    ]
    run_cmd(cmd, cwd=outdir)
    
    result_file = outdir / "SCCmecFinder_results.txt"
    if not result_file.exists():
        return "NA"
    
    import re
    match = re.search(r'Prediction based on genes:\s*(.+)', result_file.read_text())
    return match.group(1).strip() if match else "NA"

# ===== Process Sample =====
def process_sample(sample: Path, tools: Dict[str, Optional[str]], outdir: Path) -> Dict[str, str]:
    log(f"[Run] {sample.name}")
    ensure_dir(outdir)
    tmpdir = Path(tempfile.mkdtemp(prefix=f"staphscope_{sample.stem}_"))
    
    mlst_scheme, mlst_st = run_mlst(sample, tools.get("mlst"))
    spa_type = run_spa(sample, tools.get("spatyper"), tmpdir)
    sccmec_type = run_sccmec_cge(sample, tools.get("python3"), tmpdir)
    
    try:
        shutil.rmtree(tmpdir)
    except Exception:
        pass
    
    detail = {
        "sample": sample.name,
        "mlst_scheme": mlst_scheme,
        "mlst_ST": mlst_st,
        "spa_type": spa_type,
        "sccmec_type": sccmec_type
    }
    with (outdir / f"{sample.stem}.staphscope.json").open("w") as f:
        json.dump(detail, f, indent=2)
    
    return detail

# ===== Database Updates =====
def try_update_databases(tools: Dict[str, Optional[str]]):
    log("[Update] Database update requested, but this feature is not yet implemented.")

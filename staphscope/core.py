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

# ===== Utilities =====
def log(msg: str):
    print(msg, flush=True)

def run_cmd(cmd: List[str], cwd: Optional[Path] = None, env: Optional[Dict] = None) -> subprocess.CompletedProcess:
    return subprocess.run(cmd, cwd=str(cwd) if cwd else None, env=env, capture_output=True, text=True, check=False)

def ensure_dir(p: Path):
    p.mkdir(parents=True, exist_ok=True)

def get_bundled_tool_path(tool_name: str) -> Optional[Path]:
    """Dynamically get the path to bundled tools"""
    base_dir = Path(__file__).resolve().parent
    tool_paths = {
        'mlst': base_dir / "tools" / "mlst" / "bin" / "mlst",
        'spatyper': base_dir / "tools" / "spatyper" / "spa_typing" / "main" / "spaTyper",
        'sccmecfinder': base_dir / "tools" / "sccmecfinder" / "SCCmecFinder_v4.py",
        'sccmec_db': base_dir / "tools" / "sccmecfinder" / "database",
        'sccmec_script_dir': base_dir / "tools" / "sccmecfinder" / "script_dir",
        'spa_repeats': base_dir / "tools" / "spatyper" / "spa_typing" / "main" / "sparepeats.fasta",
        'spa_types': base_dir / "tools" / "spatyper" / "spa_typing" / "main" / "spatypes.txt",
    }
    path = tool_paths.get(tool_name)
    return path if path and path.exists() else None

# ===== Environment =====
def check_environment() -> Dict[str, Optional[str]]:
    tools = {
        "mlst": str(get_bundled_tool_path('mlst')) or shutil.which("mlst"),
        "spatyper": str(get_bundled_tool_path('spatyper')),
        "blastn": shutil.which("blastn"),
        "makeblastdb": shutil.which("makeblastdb"),
        "python": sys.executable
    }
    return tools

def print_environment_report(tools: Dict[str, Optional[str]]):
    log("\n[Environment Check]")
    for k, v in tools.items():
        log(f"  - {k:12s}: {v or 'NOT FOUND'}")
    
    log("\n[SCCmecFinder Check]")
    sccmec_finder = get_bundled_tool_path('sccmecfinder')
    sccmec_db = get_bundled_tool_path('sccmec_db')
    sccmec_script_dir = get_bundled_tool_path('sccmec_script_dir')
    log(f"  - SCCmecFinder_v4.py: {'Found' if sccmec_finder else 'NOT FOUND'}")
    log(f"  - SCCmec Database:   {'Found' if sccmec_db else 'NOT FOUND'}")
    log(f"  - SCCmec Script Dir: {'Found' if sccmec_script_dir else 'NOT FOUND'}")
    
    log("\n[spaTyper Check]")
    spa_repeats = get_bundled_tool_path('spa_repeats')
    spa_types = get_bundled_tool_path('spa_types')
    log(f"  - spa repeats: {'Found' if spa_repeats else 'NOT FOUND'}")
    log(f"  - spa types:   {'Found' if spa_types else 'NOT FOUND'}")

# ===== MLST =====
def run_mlst(sample: Path, mlst_bin: Optional[str]) -> Tuple[str, str]:
    if not mlst_bin or not os.access(mlst_bin, os.X_OK):
        log(f"[MLST Error] Binary not found or not executable: {mlst_bin}")
        return "NA", "NA"
    
    if not sample.exists():
        log(f"[MLST Error] Input file does not exist: {sample}")
        return "NA", "NA"
    
    # Set PERL5LIB environment variable
    env = os.environ.copy()
    perl5lib_path = Path(mlst_bin).parent.parent / "perl5"
    if perl5lib_path.exists():
        env['PERL5LIB'] = str(perl5lib_path)
    
    # Run from the directory of the mlst binary
    mlst_dir = Path(mlst_bin).parent
    cp = run_cmd([mlst_bin, str(sample)], cwd=mlst_dir, env=env)
    
    if cp.returncode != 0:
        log(f"[MLST Error] Command failed: {cp.stderr}")
        return "NA", "NA"
    
    lines = cp.stdout.strip().splitlines()
    for line in lines:
        # MLST outputs the full path, so we need to check if the line contains the sample name
        if str(sample.name) in line:
            parts = line.split("\t")
            if len(parts) >= 3:
                return parts[1], parts[2]
    
    log(f"[MLST Error] No valid output found for {sample.name}")
    log(f"[MLST Debug] Output was: {cp.stdout}")
    return "NA", "NA"
    
    lines = cp.stdout.strip().splitlines()
    for line in lines:
        if line.startswith(os.path.basename(str(sample))):
            parts = line.split("\t")
            if len(parts) >= 3:
                return parts[1], parts[2]
    
    log(f"[MLST Error] No valid output found for {sample.name}")
    return "NA", "NA"

# ===== spaTyper =====
def run_spa(sample: Path, spatyper_bin: Optional[str], tmpdir: Path) -> str:
    if not spatyper_bin:
        log("[spaTyper Error] Binary not found")
        return "NA"
    
    # The repeat and type files are in the same directory as the spaTyper script
    spa_repeats = Path(spatyper_bin).parent / "sparepeats.fasta"
    spa_types = Path(spatyper_bin).parent / "spatypes.txt"
    
    if not spa_repeats.exists() or not spa_types.exists():
        log("[spaTyper Error] Repeat or type files not found")
        return "NA"
    
    outfile = tmpdir / f"spa_result_{sample.stem}.txt"
    
    # Use current Python interpreter
    cmd = [sys.executable, spatyper_bin, "-f", str(sample), "--output", str(outfile),
           "-r", str(spa_repeats), "-o", str(spa_types)]
    
    # Set Python path for spa_typing modules
    env = os.environ.copy()
    spatyper_package_dir = Path(spatyper_bin).parent.parent  # This is the spa_typing directory
    env["PYTHONPATH"] = str(spatyper_package_dir)
    
    # Run from the directory containing the spaTyper script
    cp = run_cmd(cmd, cwd=Path(spatyper_bin).parent, env=env)
    if cp.returncode != 0:
        log(f"[spaTyper Error] Command failed: {cp.stderr}")
        return "NA"

    if outfile.exists():
        with outfile.open() as f:
            lines = [line.strip() for line in f if line.strip()]
            for line in lines:
                if "\t" in line and not line.startswith("Sequence name"):
                    parts = line.split("\t")
                    if len(parts) >= 3:
                        return parts[2]  # The type is the third column
    log(f"[spaTyper Error] No valid output found for {sample.name}")
    return "NA"

# ===== SCCmecFinder =====
def run_sccmec_cge(sample: Path, tmpdir: Path) -> str:
    sccmec_finder = get_bundled_tool_path('sccmecfinder')
    sccmec_db = get_bundled_tool_path('sccmec_db')
    sccmec_script_dir = get_bundled_tool_path('sccmec_script_dir')
    
    if not sccmec_finder or not sccmec_db or not sccmec_script_dir:
        log("[SCCmecFinder Error] Required files not found")
        return "NA"
    
    outdir = tmpdir / "sccmec"
    ensure_dir(outdir)
    
    # Copy input files
    shutil.copy(sample, outdir / "db_input.fna")
    shutil.copy(sample, outdir / "km_input.fna")
    
    # Use current Python interpreter
    cmd = [
        sys.executable, str(sccmec_finder),
        "-iDb", str(outdir / "db_input.fna"),
        "-iKm", str(outdir / "km_input.fna"),
        "-k", "90", "-l", "60",
        "-o", "SCCmecFinder_results.txt",
        "-d", str(outdir),
        "-db_dir", str(sccmec_db),
        "-sc_dir", str(sccmec_script_dir),
        "-db_choice", "reference"
    ]
    
    cp = run_cmd(cmd, cwd=outdir)
    if cp.returncode != 0:
        log(f"[SCCmecFinder Error] Command failed: {cp.stderr}")
        return "NA"
    
    result_file = outdir / "SCCmecFinder_results.txt"
    if not result_file.exists():
        log("[SCCmecFinder Error] Result file not created")
        return "NA"
    
    # Parse result
    import re
    content = result_file.read_text()
    match = re.search(r'Prediction based on genes:\s*(.+)', content)
    return match.group(1).strip() if match else "NA"

# ===== Process Sample =====
def process_sample(sample: Path, tools: Dict[str, Optional[str]], outdir: Path) -> Dict[str, str]:
    log(f"[Run] {sample.name}")
    ensure_dir(outdir)
    tmpdir = Path(tempfile.mkdtemp(prefix=f"staphscope_{sample.stem}_"))
    
    mlst_scheme, mlst_st = run_mlst(sample, tools.get("mlst"))
    spa_type = run_spa(sample, tools.get("spatyper"), tmpdir)
    sccmec_type = run_sccmec_cge(sample, tmpdir)
    
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

#!/usr/bin/env python3
"""
Staphscope CLI interface — unified MLST + spaTyper + SCCmec typing
"""

import os
import glob
import datetime
import argparse
import csv
import json
import sys
from pathlib import Path
from typing import Dict, List
from concurrent.futures import ThreadPoolExecutor, as_completed

from .core import (
    BANNER, FOOTER, check_environment, ensure_dir, log,
    print_environment_report, process_sample, try_update_databases
)

def parse_args(argv=None):
    p = argparse.ArgumentParser(
        description="Staphscope Typing Tool — unified MLST + spa + SCCmec typing "
                    "(bundled under staphscope/tools/)",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    p.add_argument("-i", "--inputs", nargs="+", required=False,
                   help="Input genome FASTA files (supports globs, e.g. '*.fasta').")
    p.add_argument("-o", "--outdir", type=lambda x: Path(x),
                   default=Path("staphscope_results"),
                   help="Output directory")
    p.add_argument("--threads", type=int, default=(os.cpu_count() or 2),
                   help="Number of parallel workers")
    p.add_argument("--check", action="store_true", help="Check environment and exit")
    p.add_argument("--update-db", action="store_true", help="Update databases before running")
    p.add_argument("--version", action="store_true", help="Print version banner and exit")
    return p.parse_args(argv)

def expand_inputs(patterns: List[str]) -> List[Path]:
    """Expand a list of input patterns into actual files."""
    files: List[Path] = []
    for pat in patterns:
        expanded = glob.glob(os.path.expanduser(pat))
        for f in expanded:
            p = Path(f)
            if p.exists() and p.is_file():
                files.append(p.resolve())
    return files

def main(argv=None):
    args = parse_args(argv)
    if args.version:
        print(BANNER)
        sys.exit(0)

    print(BANNER)
    tools = check_environment()
    print_environment_report(tools)

    if args.check:
        sys.exit(0)

    if args.update_db:
        try_update_databases(tools)
        sys.exit(0)

    if not args.inputs:
        log("\n[Error] No input files specified. Use -i to supply genomes.")
        sys.exit(2)

    samples = expand_inputs(args.inputs)
    if not samples:
        log("[Error] No input files matched.")
        sys.exit(2)

    ensure_dir(args.outdir)
    mlst_tsv = args.outdir / "mlst_results.tsv"
    spa_tsv = args.outdir / "spa_results.tsv"
    scc_tsv = args.outdir / "sccmec_results.tsv"
    combined_tsv = args.outdir / "staphscope_summary.tsv"

    # Write headers
    for f, headers in [
        (mlst_tsv, ["sample","mlst_scheme","mlst_ST"]),
        (spa_tsv, ["sample","spa_type"]),
        (scc_tsv, ["sample","sccmec_type"]),
        (combined_tsv, ["sample","mlst_scheme","mlst_ST","spa_type","sccmec_type"])
    ]:
        with f.open("w", newline="") as fh:
            csv.writer(fh, delimiter="\t").writerow(headers)

    results: List[Dict[str, str]] = []
    with ThreadPoolExecutor(max_workers=max(1, int(args.threads))) as ex:
        fut2sample = {ex.submit(process_sample, s, tools, args.outdir): s for s in samples}
        for fut in as_completed(fut2sample):
            sample = fut2sample[fut]
            try:
                res = fut.result()
                results.append(res)
                with mlst_tsv.open("a", newline="") as fh:
                    csv.writer(fh, delimiter="\t").writerow([res["sample"], res["mlst_scheme"], res["mlst_ST"]])
                with spa_tsv.open("a", newline="") as fh:
                    csv.writer(fh, delimiter="\t").writerow([res["sample"], res["spa_type"]])
                with scc_tsv.open("a", newline="") as fh:
                    csv.writer(fh, delimiter="\t").writerow([res["sample"], res["sccmec_type"]])
                with combined_tsv.open("a", newline="") as fh:
                    csv.writer(fh, delimiter="\t").writerow([
                        res["sample"], res["mlst_scheme"], res["mlst_ST"],
                        res["spa_type"], res["sccmec_type"]
                    ])
            except Exception as e:
                log(f"[Error] processing {sample}: {e}")

    meta = {
        "version": BANNER.strip(),
        "date": datetime.datetime.now().isoformat(),
        "author": "Beckley Brown",
        "department": "U.G.M.S - Department of Medical Biochemistry",
        "inputs": [str(s) for s in samples],
        "outdir": str(args.outdir.resolve()),
        "tools": {k:str(v) if v else "NA" for k,v in tools.items()}
    }
    with (args.outdir / "staphscope_run_meta.json").open("w") as f:
        json.dump(meta, f, indent=2)

    log(f"\n[Done] Wrote:\n  - {mlst_tsv}\n  - {spa_tsv}\n  - {scc_tsv}\n  - {combined_tsv}\n  - per-sample JSON in {args.outdir}")

    # Citation notice
    print("\nIf you use Staphscope Typing Tool in your research, please cite it as:\n")
    print("Brown, B. (2025). Staphscope Typing Tool: Unified MLST + spa + SCCmec typing for Staphylococcus aureus.")
    print("K.N.U.S.T/U.G.M.S - Department of Medical Biochemistry, University of Ghana.")
    print("GitHub repository: https://github.com/bbeckley-hub/staphscope-typing-tool\n")

    print(FOOTER)

if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        log("Interrupted.")
        sys.exit(130)
    except Exception as e:
        log(f"Fatal error: {e}")
        sys.exit(1)

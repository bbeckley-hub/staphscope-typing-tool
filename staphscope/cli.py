#!/usr/bin/env python3
"""
Staphscope CLI interface
"""
import datetime
import argparse
import csv
import json
import os
import sys
from pathlib import Path
from typing import Dict, List

from .core import (
    BANNER, FOOTER, check_environment, ensure_dir, log, 
    print_environment_report, process_sample, try_update_databases
)
from concurrent.futures import ThreadPoolExecutor, as_completed

def parse_args(argv=None):
    p = argparse.ArgumentParser(
        description="Staphscope Typing Tool — unified MLST + spa + SCCmec typing for Staphylococcus aureus",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("-i", "--inputs", nargs="+", required=False, help="Input genome files (FASTA/FASTQ). Accepts globs if shell-expansion is enabled.")
    p.add_argument("-o", "--outdir", type=lambda x: Path(x), default=Path("staphscope_results"), help="Output directory")
    p.add_argument("--threads", type=int, default=os.cpu_count() or 2, help="Number of parallel workers")
    p.add_argument("--check", action="store_true", help="Check environment and exit")
    p.add_argument("--update", action="store_true", help="Attempt to update underlying typing databases and exit")
    p.add_argument("--version", action="store_true", help="Print version banner and exit")
    return p.parse_args(argv)

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

    if args.update:
        try_update_databases(tools)
        sys.exit(0)

    if not args.inputs:
        log("\n[Error] No inputs specified. Use -i to supply one or more genome files.")
        sys.exit(2)

    ensure_dir(args.outdir)

    samples: List[Path] = []
    for pattern in args.inputs:
        if any(c in pattern for c in "*?["):
            for p in Path().glob(pattern):
                if p.is_file():
                    samples.append(p.resolve())
        else:
            p = Path(pattern)
            if p.exists() and p.is_file():
                samples.append(p.resolve())
    if not samples:
        log("[Error] No input files matched.")
        sys.exit(2)

    mlst_tsv = args.outdir / "mlst_results.tsv"
    spa_tsv = args.outdir / "spa_results.tsv"
    scc_tsv = args.outdir / "sccmec_results.tsv"
    combined_tsv = args.outdir / "staphscope_summary.tsv"

    with mlst_tsv.open("w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["sample", "mlst_scheme", "mlst_ST"])
    with spa_tsv.open("w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["sample", "spa_type"])
    with scc_tsv.open("w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["sample", "sccmec_type"])
    with combined_tsv.open("w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["sample", "mlst_scheme", "mlst_ST", "spa_type", "sccmec_type"])

    results: List[Dict[str, object]] = []
    with ThreadPoolExecutor(max_workers=max(1, int(args.threads))) as ex:
        fut2sample = {ex.submit(process_sample, s, tools, args.outdir): s for s in samples}
        for fut in as_completed(fut2sample):
            try:
                res = fut.result()
                results.append(res)
                with mlst_tsv.open("a", newline="") as f:
                    w = csv.writer(f, delimiter="\t")
                    w.writerow([res["sample"], res["mlst_scheme"], res["mlst_ST"]])
                with spa_tsv.open("a", newline="") as f:
                    w = csv.writer(f, delimiter="\t")
                    w.writerow([res["sample"], res["spa_type"]])
                with scc_tsv.open("a", newline="") as f:
                    w = csv.writer(f, delimiter="\t")
                    w.writerow([res["sample"], res["sccmec_type"]])
                with combined_tsv.open("a", newline="") as f:
                    w = csv.writer(f, delimiter="\t")
                    w.writerow([res["sample"], res["mlst_scheme"], res["mlst_ST"], res["spa_type"], res["sccmec_type"]])
            except Exception as e:
                log(f"Error processing sample: {e}")

    meta = {
        "version": "v0.1.0",
        "date": datetime.datetime.now().isoformat(),
        "author": "Beckley Brown",
        "affiliations": ["University of Ghana", "K.N.U.S.T"],
        "inputs": [str(s) for s in samples],
        "outdir": str(args.outdir.resolve()),
        "tools": tools,
    }
    with (args.outdir / "staphscope_run_meta.json").open("w") as f:
        json.dump(meta, f, indent=2)

    log(f"\n[Done] Wrote:\n  - {mlst_tsv}\n  - {spa_tsv}\n  - {scc_tsv}\n  - {combined_tsv}\n  - per-sample JSON in {args.outdir}")

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

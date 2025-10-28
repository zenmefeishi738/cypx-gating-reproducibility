#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
run_mmpbsa.py — One-file, reproducible wrapper for gmx_MMPBSA
=============================================================
- Generates the input file (mmpbsa_full.in) with your original settings
- Runs gmx_MMPBSA (PB with per-residue decomposition)
- Summarizes top-N residue contributors from the per-residue CSV
- Cross-platform (Linux/macOS/Windows), standard-library only

Defaults match your previous run:
  startframe=1, interval=25
  PB: istrng=0.150, indi=4.0, exdi=80.0, radiopt=0, inp=2
  Index groups: receptor=20, ligand=14
  Files: step5_1.tpr, prod_100_1000ns.xtc, topol.top, index.ndx
"""
from __future__ import annotations

import argparse
import csv
import datetime as _dt
import os
import shutil
import subprocess
import sys
from pathlib import Path
from typing import List, Tuple

def _echo(msg: str) -> None:
    ts = _dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    print(f"[{ts}] {msg}", flush=True)

def which(exe: str) -> str | None:
    return shutil.which(exe)

def write_mmpbsa_input(path: Path,
                       startframe: int,
                       interval: int,
                       istrng: float,
                       indi: float,
                       exdi: float,
                       radiopt: int,
                       inp: int,
                       idecomp: int = 1,
                       dec_verbose: int = 0) -> None:
    """Write mmpbsa_full.in with PB and per-residue decomposition blocks."""
    content = f"""&general
  startframe={startframe},
  interval={interval},
  keep_files=0,
  verbose=1
/
&pb
  istrng={istrng}, indi={indi}, exdi={exdi}, radiopt={radiopt}, inp={inp},
/
&decomp
  idecomp={idecomp},
  dec_verbose={dec_verbose}
/
"""
    path.write_text(content, encoding="utf-8")

def run_stream(cmd: List[str], log_file: Path) -> int:
    """Run a command and tee stdout/stderr to both console and a log file."""
    _echo(f"Running: {' '.join(cmd)}")
    log_file.parent.mkdir(parents=True, exist_ok=True)
    with subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True) as p, \
         log_file.open("w", encoding="utf-8", newline="") as lf:
        assert p.stdout is not None
        for line in p.stdout:
            sys.stdout.write(line)
            lf.write(line)
        return p.wait()

def summarize_topN(residue_csv: Path, out_csv: Path, topN: int = 15) -> Tuple[List[Tuple[str, float]], List[str]]:
    """Read per-residue CSV and write top-N most favorable (most negative) contributors."""
    if not residue_csv.exists():
        _echo(f"Per-residue CSV not found: {residue_csv}. Skip summary.")
        return [], []
    prefer_cols = ["ELE", "VDW", "PBSUR", "PBCAL", "PBSOL", "GBSUR", "GBCAL", "GBSOL"]
    rows: List[Tuple[str, float]] = []
    headers: List[str] = []
    with residue_csv.open(newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        headers = reader.fieldnames or []
        usable = [c for c in prefer_cols if c in headers]
        for r in reader:
            rid = r.get("Residue") or r.get("Residue name") or "UNK"
            s = 0.0
            for c in usable:
                try:
                    s += float(r.get(c, 0.0))
                except Exception:
                    pass
            rows.append((rid, s))
    rows.sort(key=lambda x: x[1])  # more negative -> stronger binding contribution
    top = rows[:topN]
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    with out_csv.open("w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["Residue", "dG_component_sum"])
        w.writerows(top)
    return top, headers

def build_argparser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="Run gmx_MMPBSA (PB, per-residue) and summarize top-N residues.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    p.add_argument("--gmx-mmpbsa", default=os.environ.get("GMX_MMPBSA", "gmx_MMPBSA"),
                   help="gmx_MMPBSA executable name or path")
    # Inputs
    p.add_argument("--tpr", default="step5_1.tpr", help="TPR file for the complex (e.g., production tpr)")
    p.add_argument("--xtc", default="prod_100_1000ns.xtc", help="trajectory (subsampled)")
    p.add_argument("--top", default="topol.top", help="system topology with ligand parameters")
    p.add_argument("--ndx", default="index.ndx", help="index file containing receptor/ligand groups")
    p.add_argument("--rec-grp", type=int, default=20, help="index group id for receptor")
    p.add_argument("--lig-grp", type=int, default=14, help="index group id for ligand")
    # Sampling & PB
    p.add_argument("--startframe", type=int, default=1, help="first frame index")
    p.add_argument("--interval", type=int, default=25, help="frame interval")
    p.add_argument("--istrng", type=float, default=0.150, help="ionic strength (PB)")
    p.add_argument("--indi", type=float, default=4.0, help="protein dielectric (PB)")
    p.add_argument("--exdi", type=float, default=80.0, help="solvent dielectric (PB)")
    p.add_argument("--radiopt", type=int, default=0, help="atomic radii set (PB)")
    p.add_argument("--inp", type=int, default=2, help="Poisson-type (PB)")
    # Outputs
    p.add_argument("--out-dat", default="MMPBSA_RESULTS.dat", help="overall results file")
    p.add_argument("--out-csv", default="MMPBSA_ENERGY_COMPONENTS.csv", help="per-residue CSV")
    p.add_argument("--topN", type=int, default=15, help="how many top residues to export")
    p.add_argument("--summary-csv", default="top_residues.csv", help="output CSV for top-N residues")
    # Misc
    p.add_argument("--workdir", default=".", help="working directory (where files are read/written)")
    p.add_argument("--input-name", default="mmpbsa_full.in", help="generated MMPBSA input filename")
    p.add_argument("--log", default="gmx_MMPBSA_run.log", help="log file capturing CLI output")
    return p

def main() -> int:
    args = build_argparser().parse_args()
    wd = Path(args.workdir).resolve()
    wd.mkdir(parents=True, exist_ok=True)

    # Resolve paths
    tpr = wd / args.tpr
    xtc = wd / args.xtc
    top = wd / args.top
    ndx = wd / args.ndx
    inp_path = wd / args.input_name
    log_path = wd / args.log
    out_dat = wd / args.out_dat
    out_csv = wd / args.out_csv
    sum_csv = wd / args.summary_csv

    # Sanity checks
    for f in [tpr, xtc, top, ndx]:
        if not f.exists():
            _echo(f"Missing file: {f}")
            return 2

    exe = which(args.gmx_mmpbsa)
    if not exe:
        _echo(f"Cannot find executable: {args.gmx_mmpbsa}. "
              f"Set --gmx-mmpbsa or export GMX_MMPBSA=/path/to/gmx_MMPBSA")
        return 3
    _echo(f"Using gmx_MMPBSA: {exe}")

    # Write input file
    write_mmpbsa_input(
        inp_path,
        startframe=args.startframe,
        interval=args.interval,
        istrng=args.istrng,
        indi=args.indi,
        exdi=args.exdi,
        radiopt=args.radiopt,
        inp=args.inp,
        idecomp=1,
        dec_verbose=0
    )
    _echo(f"Wrote input: {inp_path}")

    # Build command
    cmd = [
        exe, "-O", "-i", str(inp_path),
        "-cs", str(tpr),
        "-ct", str(xtc),
        "-cp", str(top),
        "-ci", str(ndx),
        "-cg", str(args.rec_grp), str(args.lig_grp),
        "-o", str(out_dat),
        "-eo", str(out_csv),
        "-nogui"
    ]

    # Run
    rc = run_stream(cmd, log_path)
    if rc != 0:
        _echo(f"gmx_MMPBSA exited with code {rc}. See log: {log_path}")
        return rc

    # Post-process
    top_rows, headers = summarize_topN(out_csv, sum_csv, args.topN)
    if top_rows:
        _echo(f"Top {len(top_rows)} residues written to: {sum_csv}")
        for rid, val in top_rows:
            print(f"{rid:>16s}  {val:12.4f}")
    else:
        _echo("No per-residue summary produced. (Missing CSV or empty.)")

    _echo("All done.")
    return 0

if __name__ == "__main__":
    sys.exit(main())

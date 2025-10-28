#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
occupancy_threshold_counts_closed.py

Plot grouped bar charts of H-bond counts above occupancy thresholds for four systems,
with a closed-frame style (top/right spines visible).

Differences from occupancy_threshold_counts.py:
- No numeric labels on top of bars by default (toggle with --annotate_counts)
- Closed frame: show top and right spines on the axes
- Keep integer y-axis ticks, legend, CSV export
"""

from pathlib import Path
import argparse
import sys
import math

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator


def load_occupancy_csv(path: Path) -> pd.DataFrame:
    """Read occupancy CSV and return a DataFrame with 'occupancy_percent' column."""
    if not path.exists():
        raise FileNotFoundError(str(path))

    df = pd.read_csv(path)
    if "occupancy_percent" in df.columns:
        out = df[["occupancy_percent"]].copy()
    elif "occupancy_fraction" in df.columns:
        out = df[["occupancy_fraction"]].copy()
        out["occupancy_percent"] = out["occupancy_fraction"] * 100.0
        out = out[["occupancy_percent"]]
    else:
        raise ValueError("CSV must contain 'occupancy_percent' or 'occupancy_fraction'.")
    return out.sort_values("occupancy_percent", ascending=False).reset_index(drop=True)


def boxify(ax, lw=1.2):
    """Ensure a closed frame with solid spines on all sides."""
    for side in ("top", "right", "bottom", "left"):
        ax.spines[side].set_visible(True)
        ax.spines[side].set_linewidth(lw)
    ax.tick_params(direction="out", width=lw)


def main():
    p = argparse.ArgumentParser(
        description="Plot counts above occupancy thresholds for 4 systems (closed frame)."
    )
    p.add_argument("--and_csv", default="cypx_AND_hbmap_prot_lig_AND_occupancy.csv")
    p.add_argument("--est_csv", default="cypx_est_hbmap_prot_lig_est_occupancy.csv")
    p.add_argument("--piy_csv", default="cypx_piy_hbmap_prot_lig_PIY_occupancy.csv")
    p.add_argument("--tes_csv", default="cypx_tes_hbmap_prot_lig_tes_occupancy.csv")
    p.add_argument("--thresholds", nargs="+", type=float, default=[1, 5, 10, 30],
                   help="occupancy thresholds in percent")
    p.add_argument("--out_png", default="fig_threshold_counts_box.png")
    p.add_argument("--out_pdf", default="fig_threshold_counts_box.pdf")
    p.add_argument("--out_csv", default="threshold_counts_table.csv")
    p.add_argument("--title", default="Counts above occupancy thresholds (100–1000 ns)")
    p.add_argument("--alpha", type=float, default=0.95, help="bar opacity (0–1)")
    p.add_argument("--annotate_counts", action="store_true",
                   help="if set, annotate bar tops with numeric counts")
    args = p.parse_args()

    in_files = {
        "AND": Path(args.and_csv),
        "EST": Path(args.est_csv),
        "PIY": Path(args.piy_csv),
        "TES": Path(args.tes_csv),
    }
    systems = ["PIY", "AND", "EST", "TES"]
    label_map = {
        "PIY": "CypX-PIY",
        "AND": "CypX-AND",
        "EST": "CypX-EST",
        "TES": "CypX-TES",
    }

    rows = []
    missing = []
    for sys_name in systems:
        f = in_files[sys_name]
        try:
            df = load_occupancy_csv(f)
        except Exception as e:
            missing.append(f"{sys_name}: {e}")
            continue
        for thr in args.thresholds:
            cnt = int((df["occupancy_percent"] >= float(thr)).sum())
            rows.append({"system": sys_name,
                         "threshold_percent": float(thr),
                         "count": cnt})

    if not rows:
        print("[ERROR] no usable data.", file=sys.stderr)
        if missing:
            print("\n".join(missing), file=sys.stderr)
        sys.exit(1)

    counts = pd.DataFrame(rows)
    counts["system"] = pd.Categorical(counts["system"], systems, ordered=True)
    counts.sort_values(["system", "threshold_percent"], inplace=True, ignore_index=True)
    counts.to_csv(args.out_csv, index=False)
    print(f"[OK] Saved counts table: {args.out_csv}")

    plt.rcParams.update({
        "axes.labelsize": 13,
        "xtick.labelsize": 13,
        "ytick.labelsize": 13,
        "figure.dpi": 300
    })

    thr = sorted([float(t) for t in args.thresholds])
    n_sys = len(systems)
    n_thr = len(thr)
    x = np.arange(n_sys)
    width = min(0.14, 0.7 / max(1, n_thr))
    colors = ["#4e79a7", "#f28e2c", "#59a14f", "#e15759"]

    fig, ax = plt.subplots(figsize=(9.5, 5.5))

    ymax = 0
    for j, t in enumerate(thr):
        sub = counts[counts["threshold_percent"] == t].set_index("system").reindex(systems)
        xo = x + (j - (n_thr - 1) / 2.0) * width
        vals = sub["count"].values.astype(int)

        bars = ax.bar(xo, vals, width=width,
                      color=colors[j % len(colors)], alpha=args.alpha, edgecolor="none")

        if args.annotate_counts:
            for xi, vi in zip(xo, vals):
                if vi > 0:
                    ax.text(xi, vi + 0.05, f"{vi:d}", ha="center", va="bottom", fontsize=12)

        ymax = max(ymax, vals.max() if len(vals) else 0)

    ax.set_xticks(x, [label_map[s] for s in systems])
    ax.set_ylabel("Count of H-bonds")
    ax.set_title(args.title, pad=10)
    ax.yaxis.set_major_locator(MaxNLocator(integer=True))
    ax.set_ylim(0, max(1, int(math.ceil(ymax)) + 1))

    # Legend
    legend_handles = [plt.Rectangle((0, 0), 1, 1, color=colors[i % len(colors)], alpha=args.alpha)
                      for i in range(n_thr)]
    legend_labels = [f">= {t:g}%" for t in thr]
    ax.legend(legend_handles, legend_labels, frameon=True, title=None)

    # Closed frame
    boxify(ax)

    fig.tight_layout()
    fig.savefig(args.out_png, bbox_inches="tight")
    fig.savefig(args.out_pdf, bbox_inches="tight")
    print(f"[OK] Saved figure: {args.out_png} / {args.out_pdf}")


if __name__ == "__main__":
    main()

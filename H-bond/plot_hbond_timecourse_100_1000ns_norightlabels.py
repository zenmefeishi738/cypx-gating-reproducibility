#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Plot protein–ligand H-bond count vs time for four CypX systems (2x2 grid).

This version fixes the x-axis to 100–1000 ns, keeps closed frames,
and removes any numbers/text on the **right side** of the right-column panels
(no right-side tick labels, no right-side y-axis label).
"""

from pathlib import Path
import argparse
import numpy as np
import matplotlib.pyplot as plt

# --------- Configs ---------
WINDOW_START_NS = 100.0
WINDOW_END_NS   = 1000.0
Y_MAX_FALLBACK  = 2.0        # If the max H-bond count < 2, set y upper bound to 2
LINEWIDTH       = 1.2        # Vertical line width
FIGSIZE         = (10, 6)
DPI             = 220

def read_xvg_two_cols(path: Path):
    """Read GROMACS .xvg with two numeric columns (time_ps, value).
    Skips lines starting with '@' or '#'. Returns numpy arrays (time_ps, value).
    """
    t, v = [], []
    with open(path, 'r', errors='ignore') as f:
        for line in f:
            line = line.strip()
            if not line or line[0] in ('@', '#'):
                continue
            parts = line.split()
            if len(parts) < 2:
                continue
            try:
                t.append(float(parts[0]))
                v.append(float(parts[1]))
            except ValueError:
                continue
    if len(t) == 0:
        return np.array([]), np.array([])
    return np.array(t), np.array(v)

def boxify(ax, lw=1.0):
    """Show all four spines with solid lines."""
    for side in ("top", "right", "bottom", "left"):
        ax.spines[side].set_visible(True)
        ax.spines[side].set_linewidth(lw)
    ax.tick_params(direction="out", width=lw)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--piy", default="cypx_piy_hbnum_prot_lig_PIY.xvg",
                    help="CypX-PIY hbnum .xvg")
    ap.add_argument("--and_", dest="and_", default="cypx_AND_hbnum_prot_lig_AND.xvg",
                    help="CypX-AND hbnum .xvg")
    ap.add_argument("--est", default="cypx_est_hbnum_prot_lig_est.xvg",
                    help="CypX-EST hbnum .xvg")
    ap.add_argument("--tes", default="cypx_tes_hbnum_prot_lig_tes.xvg",
                    help="CypX-TES hbnum .xvg")
    ap.add_argument("--out_png", default="fig_hbond_timecourse_100_1000ns_norightlabels.png")
    ap.add_argument("--out_pdf", default="fig_hbond_timecourse_100_1000ns_norightlabels.pdf")
    args = ap.parse_args()

    files = [
        ("CypX-PIY", Path(args.piy),  "#d62728"),  # red
        ("CypX-AND", Path(args.and_), "#17becf"),  # teal
        ("CypX-EST", Path(args.est),  "#9467bd"),  # purple
        ("CypX-TES", Path(args.tes),  "#2ca02c"),  # green
    ]

    series = []
    global_ymax = 0.0

    # Read series, convert ps->ns, and clip to 100–1000 ns
    for title, fp, color in files:
        if not fp.exists():
            print(f"[WARN] Missing file: {fp}")
            series.append((title, np.array([]), np.array([]), color))
            continue
        t_ps, y = read_xvg_two_cols(fp)
        t_ns = t_ps / 1000.0  # convert to ns
        mask = (t_ns >= WINDOW_START_NS) & (t_ns <= WINDOW_END_NS)
        t_ns = t_ns[mask]
        y = y[mask]
        series.append((title, t_ns, y, color))
        if y.size:
            global_ymax = max(global_ymax, float(np.max(y)))

    # Prepare figure (sharey=True for consistent scale)
    fig, axes = plt.subplots(2, 2, figsize=FIGSIZE, sharex=False, sharey=True)
    axes = axes.reshape(2, 2)

    y_top = max(Y_MAX_FALLBACK, global_ymax)
    if y_top <= 2.0:
        y_top = 2.0
    y_top *= 1.05  # headroom

    for r in range(2):
        for c in range(2):
            ax = axes[r, c]
            title, t_ns, y, color = series[r*2 + c]

            if t_ns.size and y.size:
                ax.vlines(t_ns, [0]*len(t_ns), y, color=color, lw=LINEWIDTH)
            else:
                ax.text(0.5, 0.5, "No data (100–1000 ns)", ha="center", va="center",
                        transform=ax.transAxes, fontsize=10, color="#666666")

            ax.set_title(title, y=1.02)
            ax.set_ylim(0, y_top)

            # x-axis fixed to 100–1000 ns
            ax.set_xlim(WINDOW_START_NS, WINDOW_END_NS)
            ax.set_xlabel("Time (ns)")
            boxify(ax)

    # Left column: standard y-axis with label
    for r in range(2):
        axes[r, 0].set_ylabel("Number of H-bonds")

    # Right column: KEEP closed frame but remove any right-side numbers/text
    for r in range(2):
        axr = axes[r, 1]
        # ensure only left ticks are shown
        axr.yaxis.set_ticks_position("left")   # ticks only on left side
        axr.tick_params(axis="y", labelright=False)  # no right-side labels
        # do not set a right-side y-label
        # keep the right spine visible via boxify()

    fig.tight_layout()
    fig.savefig(args.out_png, dpi=DPI, bbox_inches="tight")
    fig.savefig(args.out_pdf, dpi=DPI, bbox_inches="tight")
    print(f"[OK] Saved: {args.out_png}, {args.out_pdf}")

if __name__ == "__main__":
    main()

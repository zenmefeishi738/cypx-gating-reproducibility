# rmsd_plot_fixed_axes.py
# -*- coding: utf-8 -*-
"""
RMSD summary plot for CypX systems with:
- fixed axes window (so all figures can share identical inner plot size)
- larger overall figure size
- legend inside lower-right corner
- robust CSV parsing and unit handling
- pure ASCII labels (use mathtext for Angstrom to avoid encoding issues)

Author: <you>
"""

import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")  # safe for headless environments
import matplotlib.pyplot as plt


# ========= 1) Figure & style (copy these to your other figures for identical window) =========
FIG_SIZE = (7.2, 5.4)                # overall figure size (inches) - make all four figures use the same
AX_RECT  = (0.14, 0.16, 0.70, 0.72)  # (left, bottom, width, height) of the axes within the figure

TITLE_SIZE = 12
LABEL_SIZE = 11
TICK_SIZE  = 9.5
SPINE_LW   = 1.1
LINE_W     = 1.3


# ========= 2) Input files (edit paths/names as needed) =========
# key/name -> (filename, color)
FILES = {
    "CypX-APO": ("rmsd_cypx.csv", "black"),
    "CypX-PIY": ("rmsd_PIY.csv", "red"),
    "CypX-AND": ("rmsd_AND.csv", "blue"),
    "CypX-EST": ("rmsd_est.csv", "green"),
    "CypX-TES": ("rmsd_tes.csv", "purple"),
}
DATA_DIR = "."  # folder that contains the CSV files


# ========= 3) Time window and shading =========
T_EQUIL_NS = 100.0   # grey equilibration window (0-100 ns)
T_MAX_NS   = 1000.0  # x-axis upper limit


# ========= 4) Robust CSV reader (auto-detect columns and units) =========
def read_time_rmsd_ns_angstrom(path):
    """
    Returns time in ns (1D array) and RMSD in Angstrom (1D array).
    - Accepts time in ns or ps (auto-detected by header or magnitude)
    - Accepts RMSD in nm or Angstrom (auto-detected by header or magnitude)
    - Works with arbitrary column names (best-effort substring matching)
    """
    df = pd.read_csv(path, sep=None, engine="python")
    if df.shape[1] < 2:
        raise ValueError(f"Need at least 2 columns in {path}, got {df.shape[1]}")

    cols = [str(c).strip() for c in df.columns]
    low  = [c.lower() for c in cols]

    # time column
    time_hints = ["time", " ns", "(ns)", "time_ns", "time(ns)", "ps", "(ps)", "frame", "t"]
    t_idx = None
    for i, name in enumerate(low):
        if any(h in name for h in time_hints):
            t_idx = i
            break
    if t_idx is None:
        t_idx = 0  # fallback

    # rmsd column (not the time column)
    y_idx = 1 if df.shape[1] > 1 else 0
    for i, name in enumerate(low):
        if i == t_idx:
            continue
        if "rmsd" in name:
            y_idx = i
            break

    t = df.iloc[:, t_idx].astype(float).to_numpy()
    y = df.iloc[:, y_idx].astype(float).to_numpy()

    # ps -> ns (if header mentions ps, or values are very large)
    if "ps" in low[t_idx] or (np.nanmax(t) > 2e5):
        t = t / 1000.0

    # nm -> Angstrom (if header mentions nm, or values look like < 10)
    if ("nm" in low[y_idx]) or (np.nanmedian(y) < 10.0):
        y = y * 10.0

    return t, y


# ========= 5) Build the plot =========
fig = plt.figure(figsize=FIG_SIZE, dpi=600)
ax  = fig.add_axes(AX_RECT)  # FIXED inner axes window

missing = []

for label, (fname, color) in FILES.items():
    fpath = os.path.join(DATA_DIR, fname)
    if not os.path.exists(fpath):
        missing.append(fpath)
        continue

    t_ns, y_ang = read_time_rmsd_ns_angstrom(fpath)
    mask = (t_ns >= 0.0) & (t_ns <= T_MAX_NS)
    ax.plot(t_ns[mask], y_ang[mask], color=color, lw=LINE_W, label=label)

# shaded equilibration window and vertical line
ax.axvspan(0.0, T_EQUIL_NS, color="0.92", zorder=0)
ax.axvline(T_EQUIL_NS, ls="--", color="0.45", lw=1.1)

# axes labels, limits, and style (ASCII-safe Angstrom)
ax.set_xlim(0.0, T_MAX_NS)
ax.set_ylim(bottom=0)
ax.set_xlabel("Time (ns)", fontsize=LABEL_SIZE)
ax.set_ylabel(r"RMSD ($\AA$)", fontsize=LABEL_SIZE)
ax.set_title("RMSD of CypX systems", fontsize=TITLE_SIZE)

ax.tick_params(labelsize=TICK_SIZE, width=SPINE_LW*0.9, length=3)
for s in ax.spines.values():
    s.set_linewidth(SPINE_LW)

# legend INSIDE lower-right (opaque white box with thin border)
ax.legend(loc="lower right", fontsize=9,
          frameon=True, framealpha=0.90, facecolor="white",
          edgecolor="0.6", borderpad=0.5, handlelength=2)

# ========= 6) Save =========
out_png = os.path.join(DATA_DIR, "RMSD_fixed_axes_big_inlegend.png")
out_pdf = os.path.join(DATA_DIR, "RMSD_fixed_axes_big_inlegend.pdf")
fig.savefig(out_png)
fig.savefig(out_pdf)

print("Saved:", out_png, "|", out_pdf)
if missing:
    print("\nWarning: missing files:")
    for p in missing:
        print(" -", p)

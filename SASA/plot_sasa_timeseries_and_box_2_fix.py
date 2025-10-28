# plot_sasa_fixed_axes.py
# -*- coding: utf-8 -*-
"""
SASA visualization with FIXED inner axes window (match RMSD figure):
A) Time series (0–1000 ns) with 0–100 ns shaded + dashed line, legend inside LR
B) Boxplot for 100–1000 ns with Welch t vs APO (1 ns blocks) + BH-FDR and stars
"""

import re
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
from pathlib import Path
from scipy.stats import ttest_ind

# ================= Fixed drawing window (copy to all figures) =================
FIG_SIZE = (7.2, 5.4)                # overall figure size
AX_RECT  = (0.14, 0.16, 0.70, 0.72)  # (left, bottom, width, height) of axes

TITLE_SIZE = 12
LABEL_SIZE = 11
TICK_SIZE  = 9.5
SPINE_LW   = 1.1
LINE_W     = 1.2

# ================= Paths & data =================
WORKDIR = Path("/mnt/e/article_2/Rg_RASA_box")  # <-- 改成你的目录

# Time windows
T_SHOW_MIN, T_SHOW_MAX = 0.0, 1000.0
T_EQUIL = 100.0
T_BOX_MIN, T_BOX_MAX = 100.0, 1000.0

# Time-series cosmetics
SMOOTH_WIN = 5
DOWNSAMPLE = 1

# Boxplot zoom by quantiles (None for no zoom)
Y_ZOOM_QUANTILES = (0.01, 0.99)

# Outlier jitter (optional)
SHOW_JITTER = True
JITTER_POINT_SIZE = 12
JITTER_COLOR = "black"
JITTER_ALPHA = 0.30
JITTER_STD = 0.03

# Stats
BLOCK_NS = 1.0          # block length in ns
REF_HINT = "APO"        # reference keyword for statistics

INPUTS = [
    ("cypx_sasa.xvg",       "CypX-APO", "black"),
    ("cypx_piy_sasa.xvg",   "CypX-PIY", "red"),
    ("cypx_AND_sasa.xvg",   "CypX-AND", "blue"),
    ("cypx_est_sasa.xvg",   "CypX-EST", "green"),
    ("cypx_tes_sasa.xvg",   "CypX-TES", "purple"),
]

# ================= Helpers =================
def parse_units_from_xvg_header(text: str):
    t_unit = None; y_unit = None
    for line in text.splitlines():
        s = line.strip()
        if not s.startswith(("@", "#")):
            break
        if "@    xaxis" in s and "Time" in s:
            if "ps" in s.lower(): t_unit = "ps"
            if "ns" in s.lower(): t_unit = "ns"
        if "@    yaxis" in s and ("sasa" in s.lower() or "surface" in s.lower()):
            sl = s.lower()
            if "nm" in sl: y_unit = "nm2"
            if ("a^2" in sl) or ("ang" in sl) or ("å" in sl): y_unit = "a2"
    return (t_unit or "ns", y_unit or "nm2")

def read_xvg(path: Path):
    """Return (t_box, s_box, t_show, s_show) unified to ns and nm^2."""
    if not path.exists():
        return np.array([]), np.array([]), np.array([]), np.array([])
    text = path.read_text()
    t_unit, y_unit = parse_units_from_xvg_header(text)

    data = []
    for line in text.splitlines():
        s = line.strip()
        if not s or s.startswith(("#", "@")):
            continue
        parts = re.split(r"\s+", s)
        try:
            data.append([float(parts[0]), float(parts[1])])
        except Exception:
            continue
    if not data:
        return np.array([]), np.array([]), np.array([]), np.array([])

    arr = np.array(data, float)
    t, sasa = arr[:, 0], arr[:, 1]
    if t_unit == "ps": t = t / 1000.0      # ps -> ns
    if y_unit == "a2": sasa = sasa / 100.0 # Å^2 -> nm^2

    m_show = (t >= T_SHOW_MIN) & (t <= T_SHOW_MAX) & np.isfinite(sasa)
    m_box  = (t >= T_BOX_MIN ) & (t <= T_BOX_MAX ) & np.isfinite(sasa)
    return t[m_box], sasa[m_box], t[m_show], sasa[m_show]

def moving_average(y, win):
    if win <= 1 or y.size == 0:
        return y
    k = np.ones(win) / win
    return np.convolve(y, k, mode="same")

def compute_y_zoom(dlist, qpair):
    if qpair is None: return None
    allv = np.concatenate([d for d in dlist if d.size > 0]) if any(d.size for d in dlist) else np.array([])
    if allv.size == 0: return None
    qlo, qhi = np.quantile(allv, qpair)
    pad = 0.05 * (qhi - qlo) if qhi > qlo else 0.0
    return (qlo - pad, qhi + pad)

def tukey_whiskers(y):
    if y.size == 0: return np.nan, np.nan
    q1, q3 = np.percentile(y, [25, 75])
    iqr = q3 - q1
    lo_cap = q1 - 1.5 * iqr
    hi_cap = q3 + 1.5 * iqr
    lo = np.min(y[y >= lo_cap]) if np.any(y >= lo_cap) else np.min(y)
    hi = np.max(y[y <= hi_cap]) if np.any(y <= hi_cap) else np.max(y)
    return lo, hi

def block_average(t_ns, y, t0=T_BOX_MIN, block_ns=BLOCK_NS):
    if y.size == 0 or block_ns <= 0: return y
    idx = np.floor((t_ns - t0) / block_ns).astype(int)
    m = np.isfinite(y)
    idx, y = idx[m], y[m]
    if y.size == 0: return y
    maxbin = idx.max()
    sums   = np.bincount(idx, weights=y, minlength=maxbin + 1)
    counts = np.bincount(idx, minlength=maxbin + 1)
    with np.errstate(invalid="ignore"):
        means = sums / counts
    return means[np.isfinite(means)]

def fdr_bh(pvals):
    p = np.asarray(pvals, float); m = p.size
    order = np.argsort(p); ranked = p[order]
    adj = ranked * m / (np.arange(1, m + 1))
    adj = np.minimum.accumulate(adj[::-1])[::-1]
    out = np.ones_like(p); out[order] = np.minimum(adj, 1.0)
    return out

def star(q):
    return "***" if q < 1e-3 else ("**" if q < 1e-2 else ("*" if q < 0.05 else "ns"))

def fmt(x, spec):
    try:
        return "" if (x is None or (isinstance(x, float) and np.isnan(x))) else format(x, spec)
    except Exception:
        return ""

# ================= Read data =================
series = []  # (label, color, t_box, s_box, t_show, s_show_plot)
missing = []
for fn, lab, col in INPUTS:
    t_box, s_box, t_show, s_show = read_xvg(WORKDIR / fn)
    if t_box.size == 0 and t_show.size == 0:
        missing.append(fn);  continue
    if DOWNSAMPLE > 1:
        t_show = t_show[::DOWNSAMPLE]
        s_show = s_show[::DOWNSAMPLE]
    s_show_plot = moving_average(s_show, SMOOTH_WIN)
    series.append((lab, col, t_box, s_box, t_show, s_show_plot))
if missing:
    print("Missing/empty:", ", ".join(missing))
if not series:
    raise SystemExit("No usable data.")

# ================= Figure A: time series (fixed axes window) =================
figA = plt.figure(figsize=FIG_SIZE, dpi=600)
axA  = figA.add_axes(AX_RECT)

for lab, col, _, __, t_show, s_show in series:
    axA.plot(t_show, s_show, label=lab, color=col, linewidth=LINE_W)

axA.axvspan(T_SHOW_MIN, T_EQUIL, color="0.9", zorder=0)
axA.axvline(T_EQUIL, color="0.3", linestyle="--", linewidth=1.1)

axA.set_xlim(T_SHOW_MIN, T_SHOW_MAX)
axA.set_xlabel("Time (ns)", fontsize=LABEL_SIZE)
axA.set_ylabel(r"Total SASA (nm$^2$)", fontsize=LABEL_SIZE)
axA.set_title("SASA time series for CypX systems", fontsize=TITLE_SIZE)
axA.yaxis.set_major_locator(MaxNLocator(nbins=6))
axA.tick_params(labelsize=TICK_SIZE, width=SPINE_LW*0.9, length=3)
for s in axA.spines.values(): s.set_linewidth(SPINE_LW)

# legend inside lower-right
axA.legend(loc="lower right", fontsize=9,
           frameon=True, framealpha=0.90, facecolor="white",
           edgecolor="0.6", borderpad=0.5, handlelength=2)

out_ts = WORKDIR / "SASA_timeseries_fixed_axes.png"
figA.savefig(out_ts)
plt.close(figA)
print("Saved:", out_ts)

# ================= Figure B: boxplot + stats (fixed axes window) =================
labels = [lab for (lab, *_ ) in series]
data_box = [s_box for (_, _, _, s_box, _, _) in series]

# whisker limits
whisk_limits = [tukey_whiskers(y) for y in data_box]
w_lo = np.nanmin([lo for lo, _ in whisk_limits])
w_hi = np.nanmax([hi for _, hi in whisk_limits])

# zoom range
ylim_q = compute_y_zoom(data_box, Y_ZOOM_QUANTILES)
if ylim_q: ylo, yhi = ylim_q
else:      ylo, yhi = w_lo, w_hi
ylo = min(ylo, w_lo);  yhi = max(yhi, w_hi)

figB = plt.figure(figsize=FIG_SIZE, dpi=600)
axB  = figB.add_axes(AX_RECT)

bp = axB.boxplot(
    data_box,
    labels=labels,
    showfliers=True,
    whis=1.5,
    boxprops=dict(linewidth=1.2),
    medianprops=dict(color="orange", linewidth=1.4),
    whiskerprops=dict(linewidth=1.0),
    capprops=dict(linewidth=1.0),
    flierprops=dict(marker="o", markersize=5.5, markerfacecolor="none",
                    markeredgecolor="black", markeredgewidth=0.6, alpha=0.95)
)

# 关键修正：保持胡须与端帽在坐标轴内，避免越界
for artist in bp["whiskers"] + bp["caps"]:
    artist.set_clip_on(True)

# mean dots
for i, y in enumerate(data_box, start=1):
    axB.plot(i, np.mean(y), marker="o", markersize=6.0, color="black", zorder=3)

# optional jitter for outliers
if SHOW_JITTER:
    rng = np.random.default_rng(2025)
    for i, y in enumerate(data_box, start=1):
        q1, q3 = np.percentile(y, [25, 75]); iqr = q3 - q1
        lo_cap = q1 - 1.5 * iqr;  hi_cap = q3 + 1.5 * iqr
        m = (y < lo_cap) | (y > hi_cap)
        if np.any(m):
            xj = rng.normal(i, JITTER_STD, size=m.sum())
            axB.scatter(xj, y[m], s=JITTER_POINT_SIZE, alpha=JITTER_ALPHA,
                        linewidths=0, color=JITTER_COLOR, zorder=2)

axB.set_ylabel(r"Total SASA (nm$^2$)", fontsize=LABEL_SIZE)
axB.set_title("SASA distributions across CypX systems", fontsize=TITLE_SIZE)
axB.set_xticklabels(labels, rotation=20, ha="right")
axB.yaxis.set_major_locator(MaxNLocator(nbins=6))
axB.tick_params(labelsize=TICK_SIZE, width=SPINE_LW*0.9, length=3)
for s in axB.spines.values(): s.set_linewidth(SPINE_LW)

# ===== stats: Welch t vs APO on 1-ns blocks + FDR =====
series_blk = []
for (_, _, t_box, s_box, _, _) in series:
    series_blk.append(block_average(t_box, s_box, T_BOX_MIN, BLOCK_NS))

ref_idx = next((i for i, lab in enumerate(labels) if REF_HINT.lower() in lab.lower()), 0)
ref_blk = series_blk[ref_idx]

pvals = []
for i, y_blk in enumerate(series_blk):
    if i == ref_idx or len(y_blk) == 0:
        p = np.nan
    else:
        _, p = ttest_ind(ref_blk, y_blk, equal_var=False, nan_policy="omit")
    pvals.append(p)
pvals = np.array(pvals, float)
mask = np.isfinite(pvals)
qvals = pvals.copy();  qvals[mask] = fdr_bh(pvals[mask])

# 设定主 y 轴范围 + 底部轻微缓冲（关键修正）
axB.set_ylim(ylo, yhi)
ymin, ymax = axB.get_ylim()
yspan = max(ymax - ymin, 1.0)
axB.set_ylim(ymin - 0.02*yspan, ymax)  # 给底部 2% 的空间

# 写入显著性星号并给顶部留白
ymin, ymax = axB.get_ylim()
yspan = max(ymax - ymin, 1.0)
for i, q in enumerate(qvals, start=1):
    if not np.isfinite(q) or i-1 == ref_idx:  continue
    axB.text(i, ymax + 0.05*yspan, star(q), ha="center", va="bottom", fontsize=10)
axB.set_ylim(ymin, ymax + 0.12*yspan)

out_box = WORKDIR / "SASA_boxplot_fixed_axes.png"
figB.savefig(out_box)
plt.close(figB)
print("Saved:", out_box)

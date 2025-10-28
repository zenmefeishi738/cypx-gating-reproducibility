# plot_rg_box_fixed_axes.py
# -*- coding: utf-8 -*-
# Rg boxplots with fixed inner axes window to match RMSD figure

import re
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path
from scipy.stats import ttest_ind

# ========= 固定绘图区参数（与 RMSD 保持一致） =========
FIG_SIZE = (7.2, 5.4)                # 整图尺寸（英寸）
AX_RECT  = (0.14, 0.16, 0.70, 0.72)  # 轴在画布中的位置 (left, bottom, width, height)

TITLE_SIZE = 12
LABEL_SIZE = 11
TICK_SIZE  = 9.5
SPINE_LW   = 1.1

# ========= 其它设置 =========
WORKDIR    = Path("/mnt/e/article_2/Rg_RASA_box")
DISCARD_NS = 100.0
TMAX_NS    = 1000.0
BLOCK_NS   = 1.0        # 统计分块
REF_HINT   = "APO"
PLOT_RAW   = True       # True：绘原始帧；False：绘分块数据（更干净）

inputs = [
    ("cypx_gyrate.xvg",       "CypX-APO"),
    ("cypx_piy_gyrate.xvg",   "CypX-PIY"),
    ("cypx_AND_gyrate.xvg",   "CypX-AND"),
    ("cypx_est_gyrate.xvg",   "CypX-EST"),
    ("cypx_tes_gyrate.xvg",   "CypX-TES"),
]

def parse_units_from_xvg_header(text: str):
    t_unit = None; y_unit = None
    for line in text.splitlines():
        s = line.strip()
        if not s.startswith(("@", "#")):
            break
        if "@    xaxis" in s and "Time" in s:
            if "ps" in s.lower(): t_unit = "ps"
            if "ns" in s.lower(): t_unit = "ns"
        if "@    yaxis" in s and ("Rg" in s or "Gyration" in s):
            if "nm" in s.lower(): y_unit = "nm"
            if ("a)" in s.lower()) or ("å" in s.lower()) or ("ang" in s.lower()):
                y_unit = "A"
    return t_unit or "ns", y_unit or "nm"

def read_xvg_rg(path: Path):
    """Return t(ns), Rg(Å) for raw-windowed data and block-averaged data."""
    if not path.exists():
        return np.array([]), np.array([]), np.array([])
    text = path.read_text()
    t_unit, y_unit = parse_units_from_xvg_header(text)
    data = []
    for line in text.splitlines():
        s = line.strip()
        if (not s) or s.startswith(("#", "@")):
            continue
        parts = re.split(r"\s+", s)
        try:
            data.append([float(parts[0]), float(parts[1])])
        except Exception:
            continue
    if not data:
        return np.array([]), np.array([]), np.array([])
    arr = np.array(data, float)
    t = arr[:, 0]; rg = arr[:, 1]
    if t_unit == "ps": t = t / 1000.0          # ps → ns
    if y_unit.lower() == "nm": rg = rg * 10.0  # nm → Å

    m = (t >= DISCARD_NS) & (t <= TMAX_NS) & np.isfinite(rg)
    t, rg = t[m], rg[m]

    # block-averaged for stats
    if BLOCK_NS > 0:
        idx = np.floor((t - DISCARD_NS) / BLOCK_NS).astype(int)
        m2 = np.isfinite(rg)
        idx, rg2 = idx[m2], rg[m2]
        maxbin = idx.max() if rg2.size else -1
        if maxbin >= 0:
            sums = np.bincount(idx, weights=rg2, minlength=maxbin + 1)
            cnts = np.bincount(idx, minlength=maxbin + 1)
            with np.errstate(invalid="ignore"):
                blk = sums / cnts
            blk = blk[np.isfinite(blk)]
        else:
            blk = np.array([])
    else:
        blk = rg.copy()
    return t, rg, blk

def fdr_bh(pvals):
    p = np.asarray(pvals, float); m = p.size
    order = np.argsort(p); ranked = p[order]
    adj = ranked * m / (np.arange(1, m + 1))
    adj = np.minimum.accumulate(adj[::-1])[::-1]
    out = np.ones_like(p); out[order] = np.minimum(adj, 1.0)
    return out

def star(q):
    return "***" if q < 1e-3 else ("**" if q < 1e-2 else ("*" if q < 0.05 else "ns"))

# ===== 读入所有体系 =====
labels, series_raw, series_blk, missing = [], [], [], []
for fn, lab in inputs:
    t, rg_raw, rg_blk = read_xvg_rg(WORKDIR / fn)
    if rg_raw.size == 0:
        missing.append(fn);  continue
    labels.append(lab)
    series_raw.append(rg_raw)
    series_blk.append(rg_blk)

if missing:
    print("⚠ Missing or empty:", ", ".join(missing))
if not labels:
    raise SystemExit("No usable data.")

# 参考组（做显著性比较）
ref_idx = next((i for i, l in enumerate(labels) if REF_HINT.lower() in l.lower()), 0)
REF_LABEL = labels[ref_idx]
ref = series_blk[ref_idx]

# Welch t-test on block-averaged data, then FDR
pvals = []
for i, y in enumerate(series_blk):
    if i == ref_idx:
        p = np.nan
    else:
        _, p = ttest_ind(ref, y, equal_var=False, nan_policy="omit")
    pvals.append(p)
pvals = np.array(pvals, float)
mask = np.isfinite(pvals)
qvals = pvals.copy()
qvals[mask] = fdr_bh(pvals[mask])

# ========= 绘图：固定绘图区 =========
fig = plt.figure(figsize=FIG_SIZE, dpi=600)
ax  = fig.add_axes(AX_RECT)  # 关键：与 RMSD 同一内框大小

data_to_plot = series_raw if PLOT_RAW else series_blk

# 箱线图（细一点以免拥挤）
bp = ax.boxplot(
    data_to_plot,
    labels=labels,
    showfliers=True,
    whis=1.5,
    widths=0.45,
    patch_artist=True
)
# 盒子白底黑边
for box in bp['boxes']:
    box.set(facecolor='white', edgecolor='black', linewidth=1.0)
for med in bp['medians']:
    med.set(color='black', linewidth=1.0)

# 覆盖均值点
for i, y in enumerate(data_to_plot, start=1):
    ax.plot(i, np.mean(y), marker="o", markersize=3.0, color="black")

# 显著性星号（基于分块数据的 q 值）
ymax = max(np.max(y) for y in data_to_plot)
ymin = min(np.min(y) for y in data_to_plot)
yspan = max(ymax - ymin, 1.0)
for i, lab in enumerate(labels, start=1):
    if i - 1 == ref_idx:
        continue
    txt = star(qvals[i - 1]) if np.isfinite(qvals[i - 1]) else "ns"
    ax.text(i, ymax + 0.05 * yspan, txt, ha="center", va="bottom", fontsize=10)

# 坐标/样式（与 RMSD 对齐：坐标从 0 起；无网格；封闭边框）
ax.set_ylabel("Radius of gyration (Å)", fontsize=LABEL_SIZE)
ax.set_title("Rg distributions across CypX systems", fontsize=TITLE_SIZE)
ax.set_xticklabels(labels, rotation=20, ha="right")
ax.set_ylim(ymin - 0.05 * yspan, ymax + 0.12 * yspan)
ax.set_xlim(left=0.5, right=len(labels) + 0.5)  # 让箱线图左右有留白

ax.tick_params(labelsize=TICK_SIZE, width=SPINE_LW*0.9, length=3)
for s in ax.spines.values():
    s.set_linewidth(SPINE_LW)

# 保存（不使用 tight_layout，保持内框不变）
src = "raw" if PLOT_RAW else "block"
out_png = WORKDIR / f"Rg_boxplot_fixed_axes_{src}.png"
out_pdf = WORKDIR / f"Rg_boxplot_fixed_axes_{src}.pdf"
fig.savefig(out_png)
fig.savefig(out_pdf)
plt.close(fig)
print("Saved:", out_png, "|", out_pdf)

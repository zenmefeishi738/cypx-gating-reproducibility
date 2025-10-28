# plot_rmsf_fixed_axes_v2.py
# -*- coding: utf-8 -*-
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ---- 固定绘图区，与 RMSD 保持一致 ----
FIG_SIZE = (7.2, 5.4)
AX_RECT  = (0.14, 0.16, 0.70, 0.72)

TITLE_SIZE = 12
LABEL_SIZE = 11
TICK_SIZE  = 9.5
SPINE_LW   = 1.1
LINE_W     = 1.3

# 数据目录与文件
CANDIDATES = [r"E:\article_2\RMSF", "/mnt/e/article_2/RMSF", "."]
for cand in CANDIDATES:
    if os.path.isdir(cand):
        FOLDER = cand
        break

FILES = {
    "CypX-Apo": ("cypx_apo_rmsf.csv",       "black"),
    "CypX-PIY": ("cypx_3NC7_12A_rmsf.csv",  "red"),
    "CypX-AND": ("cypx_AND_rmsf.csv",       "blue"),
    "CypX-EST": ("cypx_est_rmsf.csv",       "green"),
    "CypX-TES": ("cypx_tes_rmsf.csv",       "purple"),
}

def read_rmsf(path):
    df = pd.read_csv(path, sep=None, engine="python")
    if df.shape[1] < 2:
        raise ValueError(f"Need at least 2 columns in {path}")
    cols = [str(c).strip() for c in df.columns]
    low  = [c.lower() for c in cols]

    # residue index
    res_hints = ["residue", "resid", "index", "residue index", "i"]
    x_idx = None
    for i, name in enumerate(low):
        if any(h in name for h in res_hints):
            x_idx = i; break
    if x_idx is None: x_idx = 0

    # rmsf value
    y_idx = 1 if df.shape[1] > 1 else 0
    for i, name in enumerate(low):
        if i == x_idx: continue
        if "rmsf" in name:
            y_idx = i; break

    x = df.iloc[:, x_idx].astype(float).to_numpy()
    y = df.iloc[:, y_idx].astype(float).to_numpy()

    # 单位 nm -> Å
    if ("nm" in low[y_idx]) or (np.nanmedian(y) < 10.0):
        y = y * 10.0
    return x, y

# ---------- 绘图（固定绘图区） ----------
fig = plt.figure(figsize=FIG_SIZE, dpi=600)
ax  = fig.add_axes(AX_RECT)

xmax = 0.0
missing = []
for label, (fname, color) in FILES.items():
    fpath = os.path.join(FOLDER, fname)
    if not os.path.exists(fpath):
        missing.append(fpath); continue
    x, y_ang = read_rmsf(fpath)
    ax.plot(x, y_ang, label=label, color=color, linewidth=LINE_W)
    if len(x) and np.nanmax(x) > xmax:
        xmax = float(np.nanmax(x))

# ---- 坐标从 0 开始（你要的微调）----
ax.set_xlim(left=0, right=xmax if xmax > 0 else None)  # x 轴从 0 起
ax.set_ylim(bottom=0)                                   # y 轴从 0 起

# 样式
ax.set_xlabel("Residue index", fontsize=LABEL_SIZE)
ax.set_ylabel(r"RMSF ($\AA$)", fontsize=LABEL_SIZE)
ax.set_title("RMSF of CypX Systems", fontsize=TITLE_SIZE)
ax.tick_params(labelsize=TICK_SIZE, width=SPINE_LW*0.9, length=3)
for s in ax.spines.values():
    s.set_linewidth(SPINE_LW)

# ---- 图例：图内左上角（你要的微调）----
ax.legend(loc="upper left", bbox_to_anchor=(0.02, 0.98),
          fontsize=9, frameon=True, framealpha=0.90,
          facecolor="white", edgecolor="0.6", borderpad=0.5,
          handlelength=2)

# 可选：高亮 mouth/cover（如需，取消注释）
# for a,b,c in [(55,90,'orange'), (172,189,'skyblue'),
#               (198,203,'skyblue'), (224,241,'skyblue')]:
#     ax.axvspan(a, b, color=c, alpha=0.15, lw=0)

# 保存
out_png = os.path.join(FOLDER, "RMSF_fixed_axes_big_ULlegend.png")
out_pdf = os.path.join(FOLDER, "RMSF_fixed_axes_big_ULlegend.pdf")
fig.savefig(out_png)
fig.savefig(out_pdf)
print("Saved:", out_png, "|", out_pdf)
if missing:
    print("Missing files:"); [print(" -", p) for p in missing]

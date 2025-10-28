# -*- coding: utf-8 -*-
"""
Plot PCA projection (PC1 vs PC2) for five systems.

- 从 WD 目录读取各体系的 2D 投影：每行形如 `idx  pc1  pc2`（或列≥3，取第2、3列）。
- 从 CSV 读取各体系 PC1/PC2 的 explained variance 百分比：
  需包含列：File, PC1 (%), PC2 (%), [PC1+PC2 (%)]（列名大小写/空格/括号不敏感）。
- 图例：
    * 标签中的下划线(_)替换为连字符(-)
    * 仅对 Apo 不附加 (PC1/PC2) 百分比
    * 放大图例圆点（markerscale）
"""

import re
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
from pathlib import Path
from os.path import basename, splitext

# ---------------------- I/O 与样式参数 ----------------------
WD = Path("/mnt/e/article_2/PCA")                 # 工作目录
CSV_PATH = WD / "pca_variance_summary.csv"        # 方差比例 CSV

FILES = [
    ("cypx_pca_proj_2d.dat",     "CypX_Apo", "black"),
    ("cypx_piy_pca_proj_2d.dat", "CypX_PIY", "red"),
    ("cypx_AND_pca_proj_2d.dat", "CypX_AND", "blue"),
    ("cypx_est_pca_proj_2d.dat", "CypX_EST", "green"),
    ("cypx_tes_pca_proj_2d.dat", "CypX_TES", "purple"),
]

DOWNSAMPLE = 5      # 帧降采样步长（>1 表示取每 N 帧）
POINT_SIZE = 3      # 散点大小
ALPHA = 0.25        # 散点透明度
LEGEND_MARKERSCALE = 2.5   # 图例圆点放大倍数，仅影响图例

# ---------------------- 工具函数 ----------------------
def load_three_col(path: Path) -> np.ndarray:
    """读取 2D 投影数据文件：取第 2、3 列为 PC1/PC2。"""
    xs, ys = [], []
    for line in path.read_text().splitlines():
        s = line.strip()
        if not s or s.startswith(("@", "#")):
            continue
        parts = re.split(r"\s+", s)
        if len(parts) < 3:
            continue
        try:
            xs.append(float(parts[1]))  # PC1
            ys.append(float(parts[2]))  # PC2
        except ValueError:
            continue
    if not xs:
        return np.empty((0, 2))
    return np.column_stack([xs, ys])

def normalize_colname(name: str) -> str:
    # 去掉空格、括号、百分号、下划线，转小写
    return re.sub(r"[\s_%()]+", "", str(name)).lower()

def normalize_file_token(s: str) -> str:
    """把 File 列内容规整成 key（忽略大小写/扩展名/路径）。"""
    base = splitext(basename(str(s)))[0].lower()
    for pref in ("cypx_", "cypx-", "cypx"):
        if base.startswith(pref):
            break
    return base

# File token -> 体系标签
TOKEN2LABEL = {
    "cypx": "CypX_Apo",
    "cypx_apo": "CypX_Apo",
    "cypx_piy": "CypX_PIY",
    "cypx_and": "CypX_AND",
    "cypx_est": "CypX_EST",
    "cypx_tes": "CypX_TES",
}

def read_variance_csv(csv_path: Path):
    """读取各体系 PC1/PC2 explained variance 百分比。"""
    df = pd.read_csv(csv_path)
    if df.empty:
        print("⚠️ CSV 为空")
        return {}, (None, None)

    df = df.rename(columns={c: normalize_colname(c) for c in df.columns})
    need = {"file", "pc1", "pc2"}
    if not need.issubset(df.columns):
        print("⚠️ 未找到期望列。CSV 列为：", list(df.columns))
        print("  需要列名（大小写/空白/括号不敏感）：'File', 'PC1 (%)', 'PC2 (%)'")
        return {}, (None, None)

    pct_map = {}
    for _, row in df.iterrows():
        token = normalize_file_token(row["file"])
        label = TOKEN2LABEL.get(token)
        if not label:
            for t, lab in TOKEN2LABEL.items():
                if t in token:
                    label = lab
                    break
        if not label:
            continue
        try:
            pc1 = float(row["pc1"])
            pc2 = float(row["pc2"])
        except Exception:
            continue
        pct_map[label] = (pc1, pc2)

    default_pair = pct_map.get("CypX_Apo")
    if default_pair is None and pct_map:
        default_pair = next(iter(pct_map.values()))
    return pct_map, default_pair if default_pair else (None, None)

# ---------------------- 读取数据 ----------------------
pct_map, default_pair = read_variance_csv(CSV_PATH)
default_pc1, default_pc2 = default_pair

series, missing = [], []
for fname, label, color in FILES:
    f = WD / fname
    if not f.exists():
        missing.append(fname + " (not found)")
        continue
    arr = load_three_col(f)
    if arr.size == 0:
        missing.append(fname + " (empty)")
        continue
    if DOWNSAMPLE > 1:
        arr = arr[::DOWNSAMPLE]
    series.append((label, color, arr))

# ---------------------- 绘图 ----------------------
plt.figure(figsize=(7.8, 6.6))

for label, color, arr in series:
    # 图例美化：下划线 -> 连字符
    pretty = label.replace("_", "-")

    # 仅非 Apo 附加 (PC1/PC2)
    if (label in pct_map) and (label != "CypX_Apo"):
        pc1, pc2 = pct_map[label]
        legend_label = f"{pretty} (PC1 {pc1:.1f}%, PC2 {pc2:.1f}%)"
    else:
        legend_label = pretty

    plt.scatter(arr[:, 0], arr[:, 1],
                s=POINT_SIZE, alpha=ALPHA, color=color, label=legend_label)

# 坐标轴：默认使用 Apo 的方差比例（若无则仅写 PC1/PC2）
if default_pc1 is not None and default_pc2 is not None:
    plt.xlabel(f"PC1 ({default_pc1:.1f}%)")
    plt.ylabel(f"PC2 ({default_pc2:.1f}%)")
else:
    plt.xlabel("PC1")
    plt.ylabel("PC2")

plt.title("PCA projection (PC1 vs PC2)")

plt.legend(
    frameon=False,
    ncol=2,
    scatterpoints=1,
    markerscale=LEGEND_MARKERSCALE,  # 放大图例圆点
    handletextpad=0.6,
    borderpad=0.2,
    labelspacing=0.6
)

plt.tight_layout()
out_png = WD / "PCA_PC1_PC2_scatter.png"
plt.savefig(out_png, dpi=350)
plt.close()
print("✅ Saved:", out_png)

if missing:
    print("⚠️ Missing/empty inputs:\n  - " + "\n  - ".join(missing))

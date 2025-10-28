# -*- coding: utf-8 -*-
"""
make_fes_mc.py —— Mouth–Cover 2D/3D Free Energy Surface (FES)

特点
- 自动匹配 mouth/cover 数据：优先识别 d_mouth_A_<TAG>.xvg / d_cover_A_<TAG>.xvg，
  也兼容 <TAG>_mouth(.csv/.xvg) + <TAG>_cover(.csv/.xvg)，以及单文件同时包含两列 mouth/cover 的情况。
- 计算 F = -kT ln P（kJ/mol），支持像素尺度高斯平滑、低计数过滤、范围裁剪、屋顶遮罩等。
- 缺省坐标：mouth 27–32 Å（X 轴）、cover 20–27 Å（Y 轴）；自由能上限 vmax=12 kJ/mol。

用法（WSL）
    cd /mnt/e/article_2/FES/mouth_cover
    python3 make_fes_mc.py --root . --out ./out_mc --mode 3d
或（关闭屋顶遮罩）
    python3 make_fes_mc.py --root . --out ./out_mc --mode 3d --no-mask-plateau
"""

import re
import argparse
from pathlib import Path
import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import cm

# -------------------- I/O 与匹配 --------------------

def _norm_colname(s: str) -> str:
    return re.sub(r"[\s_%()]+", "", str(s)).lower()

def _load_xvg_or_csv(path: Path) -> pd.DataFrame:
    """读取 CSV 或 XVG（去掉以@/#开头行），返回 DataFrame。"""
    if not path.exists():
        raise FileNotFoundError(path)

    if path.suffix.lower() == ".csv":
        return pd.read_csv(path)

    rows = []
    with path.open("r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            s = line.strip()
            if not s or s.startswith("@") or s.startswith("#"):
                continue
            parts = re.split(r"\s+", s)
            try:
                rows.append([float(x) for x in parts])
            except Exception:
                continue
    if not rows:
        return pd.DataFrame()
    maxlen = max(len(r) for r in rows)
    data = np.array([r + [np.nan]*(maxlen-len(r)) for r in rows], float)
    cols = [f"c{i}" for i in range(data.shape[1])]
    return pd.DataFrame(data, columns=cols)

def _find_value_col(df: pd.DataFrame) -> str:
    """选择“最后一个非时间列”作为数值列；若未识别到时间列则取最后一列。"""
    drop = []
    for c in df.columns:
        nc = _norm_colname(c)
        if "time" in nc or nc.endswith("ps") or nc.endswith("ns"):
            drop.append(c)
    cand = [c for c in df.columns if c not in drop]
    if not cand:
        cand = list(df.columns)
    return cand[-1]

def _find_col_by_keys(df: pd.DataFrame, keys) -> str | None:
    if isinstance(keys, str):
        keys = [keys]
    cmap = {_norm_colname(c): c for c in df.columns}
    # 完全匹配
    for k in keys:
        nk = _norm_colname(k)
        if nk in cmap:
            return cmap[nk]
    # 宽松包含
    for want in keys:
        wn = _norm_colname(want)
        for nk, orig in cmap.items():
            if wn in nk:
                return orig
    return None

def _load_singlefile_mc(path: Path):
    """从单文件读取 mouth/cover 两列。读取失败返回 (None, None)。"""
    df = _load_xvg_or_csv(path)
    if df.empty:
        return None, None

    mouth_keys = ["mouth_a", "mouth", "mouthdistance", "mouthang"]
    cover_keys = ["cover_a", "cover", "coverdistance", "coverdeg"]

    cm = _find_col_by_keys(df, mouth_keys)
    cc = _find_col_by_keys(df, cover_keys)

    # 如果没识别到列名，退化为“最后两列”
    if cm is None or cc is None:
        cand = list(df.columns)
        if len(cand) >= 2:
            cm, cc = cand[-2], cand[-1]
        else:
            return None, None

    try:
        m = pd.to_numeric(df[cm], errors="coerce").to_numpy()
        c = pd.to_numeric(df[cc], errors="coerce").to_numpy()
        mask = np.isfinite(m) & np.isfinite(c)
        m, c = m[mask], c[mask]
        n = min(len(m), len(c))
        return m[:n], c[:n]
    except Exception:
        return None, None

def _load_twofiles_mc(path_m: Path, path_c: Path):
    """从两文件分别读取 mouth/cover。读取失败返回 (None, None)。"""
    dm = _load_xvg_or_csv(path_m)
    dc = _load_xvg_or_csv(path_c)
    if dm.empty or dc.empty:
        return None, None
    try:
        cm = _find_value_col(dm)
        cc = _find_value_col(dc)
        m = pd.to_numeric(dm[cm], errors="coerce").to_numpy()
        c = pd.to_numeric(dc[cc], errors="coerce").to_numpy()
        mask = np.isfinite(m) & np.isfinite(c)
        m, c = m[mask], c[mask]
        n = min(len(m), len(c))
        return m[:n], c[:n]
    except Exception:
        return None, None

def _autodetect_mc(root: Path, tag: str):
    """
    按常见命名自动匹配 mouth/cover：
    - 优先：d_mouth_A_<TAG>.xvg + d_cover_A_<TAG>.xvg
    - 其次：<TAG>_mouth.* + <TAG>_cover.* 或 mouth_<TAG>.* + cover_<TAG>.*
    - 最后：<TAG>_mc.* / mc_<TAG>.* (单文件双列)
    """
    # 优先：d_mouth_A_TAG / d_cover_A_TAG
    pm = root / f"d_mouth_A_{tag}.xvg"
    pc = root / f"d_cover_A_{tag}.xvg"
    if pm.exists() and pc.exists():
        m, c = _load_twofiles_mc(pm, pc)
        if m is not None and c is not None and len(m) > 0:
            return m, c

    # 两文件：mouth + cover
    mouth_names = [
        f"{tag}_mouth.csv", f"mouth_{tag}.csv", f"{tag}_mouth.xvg", f"mouth_{tag}.xvg"
    ]
    cover_names = [
        f"{tag}_cover.csv", f"cover_{tag}.csv", f"{tag}_cover.xvg", f"cover_{tag}.xvg"
    ]
    for mn in mouth_names:
        pm = root / mn
        if not pm.exists():
            continue
        for cn in cover_names:
            pc = root / cn
            if not pc.exists():
                continue
            m, c = _load_twofiles_mc(pm, pc)
            if m is not None and c is not None and len(m) > 0:
                return m, c

    # 单文件双列
    for name in [f"{tag}_mc.csv", f"mc_{tag}.csv", f"{tag}_mouth_cover.csv",
                 f"{tag}_mc.xvg", f"mc_{tag}.xvg", f"{tag}_mouth_cover.xvg"]:
        p = root / name
        if p.exists():
            m, c = _load_singlefile_mc(p)
            if m is not None and c is not None and len(m) > 0:
                return m, c

    # 最后：扫描更宽松的 mouth/cover + tag 组合
    files = list(root.glob("*"))
    mouths = [p for p in files if re.search(fr"(mouth).*({tag})|({tag}).*(mouth)", p.name, re.I)]
    covers = [p for p in files if re.search(fr"(cover).*({tag})|({tag}).*(cover)", p.name, re.I)]
    for pm in mouths:
        for pc in covers:
            m, c = _load_twofiles_mc(pm, pc)
            if m is not None and c is not None and len(m) > 0:
                return m, c

    return None, None

# -------------------- 统计与 FES --------------------

def _gaussian_kernel_2d(sigma_pix: float):
    if sigma_pix <= 0:
        K = np.zeros((1, 1))
        K[0, 0] = 1.0
        return K
    r = int(3 * sigma_pix)
    x = np.arange(-r, r + 1)
    y = np.arange(-r, r + 1)
    X, Y = np.meshgrid(x, y, indexing="xy")
    K = np.exp(-(X**2 + Y**2) / (2.0 * sigma_pix**2))
    K /= K.sum()
    return K

def _smooth2d(H, sigma_pix: float):
    if sigma_pix <= 0:
        return H
    K = _gaussian_kernel_2d(sigma_pix)
    kx, ky = K.shape
    rx, ry = kx // 2, ky // 2
    Hp = np.pad(H, ((rx, rx), (ry, ry)), mode="constant", constant_values=0.0)
    out = np.zeros_like(H, dtype=float)
    for i in range(out.shape[0]):
        for j in range(out.shape[1]):
            patch = Hp[i:i+kx, j:j+ky]
            out[i, j] = np.sum(patch * K)
    return out

def make_fes(x, y, bins, xmin, xmax, ymin, ymax, tempK,
             count_min=3, sigma=1.0, vmax=12.0):
    """
    x=mouth, y=cover。计算二维直方、平滑、概率、自由能。
    返回 (Xcenters, Ycenters, F)；若样本太少返回 None。
    """
    # 只取范围内数据
    mask = (x >= xmin) & (x <= xmax) & (y >= ymin) & (y <= ymax)
    x = x[mask]; y = y[mask]
    if len(x) < 10:
        return None

    H, xe, ye = np.histogram2d(x, y, bins=bins,
                               range=[[xmin, xmax], [ymin, ymax]])
    H[H < count_min] = 0.0
    Hs = _smooth2d(H, float(sigma))

    total = np.sum(Hs)
    P = Hs / total if total > 0 else Hs

    kB = 0.0083144621  # kJ/(mol·K)
    with np.errstate(divide="ignore", invalid="ignore"):
        F = -kB * float(tempK) * np.log(P)

    # 以最小可见能量为零点；空值设为最大
    finite = np.isfinite(F)
    if not np.any(finite):
        return None
    F = F - np.nanmin(F[finite])
    F[~finite] = np.nanmax(F[finite])
    if vmax is not None:
        F = np.clip(F, 0, vmax)

    Xc = 0.5 * (xe[:-1] + xe[1:])
    Yc = 0.5 * (ye[:-1] + ye[1:])
    return Xc, Yc, F

# -------------------- 绘图 --------------------

def plot_fes_2d(X, Y, F, tag, out_png, xlab, ylab, vmax=12.0, mask_plateau=True):
    Xg, Yg = np.meshgrid(X, Y, indexing="xy")
    Z = F.T.copy()
    if mask_plateau:
        thr = 0.98 * vmax
        Z = np.where(Z >= thr, np.nan, Z)

    plt.figure(figsize=(7.2, 6.0))
    cf = plt.contourf(Xg, Yg, Z, levels=30, cmap="turbo")
    cbar = plt.colorbar(cf)
    cbar.set_label("Free Energy (kJ/mol)")
    cs = plt.contour(Xg, Yg, F.T, levels=10, colors="k", linewidths=0.3, alpha=0.5)
    plt.clabel(cs, inline=True, fontsize=7, fmt="%.1f")
    plt.xlabel(xlab); plt.ylabel(ylab)
    plt.title(f"FES: Mouth–Cover ({tag})")
    plt.tight_layout()
    plt.savefig(out_png, dpi=350)
    plt.close()

def plot_fes_3d(X, Y, F, tag, out_png, xlab, ylab, vmax=12.0, mask_plateau=True):
    from mpl_toolkits.mplot3d import Axes3D  # noqa: F401

    Xg, Yg = np.meshgrid(X, Y, indexing="xy")
    Z = F.T.copy()
    if mask_plateau:
        thr = 0.98 * vmax
        Z = np.where(Z >= thr, np.nan, Z)

    fig = plt.figure(figsize=(9.5, 8.2))
    ax = fig.add_subplot(111, projection="3d")

    surf = ax.plot_surface(Xg, Yg, Z, cmap=cm.turbo, linewidth=0, antialiased=True)
    # 投影等高线（使用原始 F，便于读数）
    ax.contour(Xg, Yg, F.T, zdir='z', offset=np.nanmin(F), cmap="turbo", levels=20, alpha=0.9)

    cb = fig.colorbar(surf, shrink=0.8, pad=0.12)
    cb.set_label("Free Energy (kJ/mol)")

    ax.set_xlabel(xlab)
    ax.set_ylabel(ylab)
    ax.set_zlabel("Free Energy (kJ/mol)")
    ax.set_title(f"FES: Mouth–Cover ({tag})")
    ax.view_init(elev=25, azim=-45)
    plt.tight_layout()
    plt.savefig(out_png, dpi=350)
    plt.close()

# -------------------- 主流程 --------------------

def main():
    ap = argparse.ArgumentParser(description="Build Mouth–Cover FES (2D/3D)")
    ap.add_argument("--root", type=str, required=True, help="数据根目录（含 CSV/XVG）")
    ap.add_argument("--out", type=str, required=True, help="输出目录")
    ap.add_argument("--mode", choices=["2d", "3d"], default="3d", help="绘图模式")
    ap.add_argument("--temp", type=float, default=300.0, help="温度 K")
    ap.add_argument("--bins", type=int, default=60, help="每轴网格数")
    ap.add_argument("--sigma", type=float, default=1.5, help="像素尺度高斯平滑 sigma")
    ap.add_argument("--count-min", type=int, default=1, help="计数阈值（小于则置零）")

    # 统一坐标范围（默认 mouth: 27–32 Å；cover: 20–27 Å）
    ap.add_argument("--x-min", type=float, default=27.0, help="mouth 最小 (Å)")
    ap.add_argument("--x-max", type=float, default=32.0, help="mouth 最大 (Å)")
    ap.add_argument("--y-min", type=float, default=20.0, help="cover 最小 (Å)")
    ap.add_argument("--y-max", type=float, default=27.0, help="cover 最大 (Å)")

    ap.add_argument("--vmax", type=float, default=12.0, help="自由能上限 (kJ/mol)")
    ap.add_argument("--systems", type=str, default="APO,PIY,AND,EST,TES",
                    help="体系列表，逗号分隔")
    ap.add_argument("--x-label", type=str, default="Mouth distance (Å)")
    ap.add_argument("--y-label", type=str, default="Cover distance (Å)")

    # 屋顶遮罩（默认开启；若要关闭，传 --no-mask-plateau）
    ap.add_argument("--no-mask-plateau", action="store_true",
                    help="禁用屋顶遮罩（大于 ~0.98*vmax 的区域不再置 NaN）")

    args = ap.parse_args()

    ROOT = Path(args.root).resolve()
    OUT = Path(args.out).resolve()
    OUT.mkdir(parents=True, exist_ok=True)

    tags = [t.strip() for t in args.systems.split(",") if t.strip()]

    print(f"[INFO] root={ROOT}")
    print(f"[INFO] out ={OUT}")
    print(f"[INFO] systems={tags}")
    print(f"[INFO] bins={args.bins}, sigma={args.sigma}, count_min={args.count_min}")
    print(f"[INFO] X∈[{args.x_min}, {args.x_max}], Y∈[{args.y_min}, {args.y_max}], vmax={args.vmax}")

    for tag in tags:
        m, c = _autodetect_mc(ROOT, tag)
        if m is None or c is None:
            print(f"[WARN] {tag}: 未找到 mouth/cover 数据，跳过。")
            continue

        res = make_fes(
            x=m, y=c,
            bins=int(args.bins),
            xmin=float(args.x_min), xmax=float(args.x_max),
            ymin=float(args.y_min), ymax=float(args.y_max),
            tempK=float(args.temp),
            count_min=int(args.count_min),
            sigma=float(args.sigma),
            vmax=float(args.vmax)
        )
        if res is None:
            print(f"[WARN] {tag}: 有效样本不足或全被过滤，跳过。")
            continue

        X, Y, F = res
        out_png = OUT / f"FES_Mouth_Cover_{tag}_{args.mode}.png"
        mask_plateau = (not args.no_mask_plateau)

        if args.mode == "2d":
            plot_fes_2d(X, Y, F, tag, out_png,
                        args.x_label, args.y_label,
                        vmax=args.vmax, mask_plateau=mask_plateau)
        else:
            plot_fes_3d(X, Y, F, tag, out_png,
                        args.x_label, args.y_label,
                        vmax=args.vmax, mask_plateau=mask_plateau)

        print(f"[OK] {tag} → {out_png}")

if __name__ == "__main__":
    main()

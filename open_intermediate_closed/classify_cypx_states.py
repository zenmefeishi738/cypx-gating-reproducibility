#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import argparse, re, numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

def read_xvg(path):
    xs, ys = [], []
    with open(path, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            if not line.strip() or line[0] in "#@":
                continue
            cols = re.split(r"\s+", line.strip())
            if len(cols) >= 2:
                xs.append(float(cols[0])); ys.append(float(cols[1]))
    return np.array(xs), np.array(ys)  # time, value

def read_two_series(mouth_xvg, cover_xvg):
    t1, m = read_xvg(mouth_xvg)
    t2, c = read_xvg(cover_xvg)
    # 对齐（以较短者为准）
    n = min(len(m), len(c))
    return t1[:n], m[:n], c[:n]

def kmeans3(data, max_iter=200, seed=0):
    """简单 k-means (k=3)"""
    rng = np.random.default_rng(seed)
    k = 3
    # 用分位点初始化更稳
    qs = np.quantile(data, [0.15, 0.50, 0.85], axis=0)
    centroids = qs.copy()
    for _ in range(max_iter):
        # assign
        d2 = ((data[:,None,:]-centroids[None,:,:])**2).sum(axis=2)
        lab = np.argmin(d2, axis=1)
        new = np.vstack([data[lab==i].mean(axis=0) if np.any(lab==i) else centroids[i] for i in range(k)])
        if np.allclose(new, centroids, atol=1e-6):
            break
        centroids = new
    return lab, centroids

def main():
    ap = argparse.ArgumentParser(description="用 mouth/cover 将轨迹分为 Open/Intermediate/Closed 三态")
    ap.add_argument("--mouth", required=True, help="d_mouth.xvg 或两列文件（Å）")
    ap.add_argument("--cover", required=True, help="d_cover.xvg 或两列文件（Å）")
    ap.add_argument("--out", default="state_assignments.csv")
    ap.add_argument("--scatter", default="states_scatter.png")
    ap.add_argument("--use_cuts", action="store_true", help="使用阈值法而非 k-means")
    ap.add_argument("--cuts", nargs=4, type=float,
                   metavar=("MOUTH_LOW","MOUTH_HIGH","COVER_LOW","COVER_HIGH"),
                   help="阈值法：Closed≤(low) & ≤(low)；Open≥(high) & ≥(high)；其余为Intermediate")
    return ap.parse_args()

if __name__ == "__main__":
    args = main()
    t_ps, mouth, cover = read_two_series(args.mouth, args.cover)   # Å, Å
    X = np.column_stack([mouth, cover])

    if args.use_cuts and args.cuts is not None:
        ml, mh, cl, ch = args.cuts
        lab = np.full(len(X), 1)  # 中间态=1
        lab[(X[:,0] <= ml) & (X[:,1] <= cl)] = 0  # Closed
        lab[(X[:,0] >= mh) & (X[:,1] >= ch)] = 2  # Open
        cent = np.array([X[lab==i].mean(axis=0) for i in (0,1,2)])
    else:
        lab_raw, cent = kmeans3(X)
        # 重命名: 按 mouth+cover 的均值从小到大 → Closed, Inter, Open
        order = np.argsort(cent.sum(axis=1))
        remap = {old:new for new, old in enumerate(order)}
        lab = np.array([remap[i] for i in lab_raw])
        cent = cent[order]

    state_names = np.array(["Closed","Intermediate","Open"])
    # 输出 CSV
    with open(args.out, "w") as f:
        f.write("time_ps,mouth_A,cover_A,state\n")
        for tt, m, c, s in zip(t_ps, mouth, cover, lab):
            f.write(f"{tt:.1f},{m:.3f},{c:.3f},{state_names[s]}\n")

    # 每个态挑“代表帧”（到该态质心最近的一个）
    reps = {}
    for i, name in enumerate(state_names):
        idx = np.where(lab==i)[0]
        if len(idx)==0: continue
        d2 = ((X[idx]-cent[i])**2).sum(axis=1)
        reps[name] = int(idx[np.argmin(d2)])

    # 写出各态时间列表（ps）
    for i, name in enumerate(state_names):
        idx = np.where(lab==i)[0]
        with open(f"{name.lower()}_times.txt","w") as f:
            for k in idx:
                f.write(f"{t_ps[k]:.1f}\n")
    with open("representatives.txt","w") as f:
        for name, k in reps.items():
            f.write(f"{name},{t_ps[k]:.1f},{mouth[k]:.3f},{cover[k]:.3f}\n")

    # 占比统计
    counts = np.array([(lab==i).sum() for i in range(3)])
    frac = counts / len(lab)
    print("[Occupancy] Closed/Inter/Open =", frac)

    # 画散点（颜色=三态，带质心）
    colors = np.array([[0.1,0.5,1.0,0.8],[0.3,0.8,0.3,0.8],[1.0,0.5,0.1,0.8]])
    plt.figure(figsize=(6,5),dpi=180)
    for i,name in enumerate(state_names):
        idx = lab==i
        plt.scatter(mouth[idx], cover[idx], s=5, c=[colors[i]], label=f"{name} ({frac[i]*100:.1f}%)")
        plt.plot(cent[i,0], cent[i,1], "k*", ms=12)
    plt.xlabel("BC–FG mouth distance (Å)")
    plt.ylabel("FG loop–heme Fe distance (Å)")
    plt.legend(frameon=False, loc="best")
    plt.tight_layout()
    plt.savefig(args.scatter, dpi=300)
    print(f"[OK] wrote {args.out}, representatives.txt, *_times.txt and {args.scatter}")

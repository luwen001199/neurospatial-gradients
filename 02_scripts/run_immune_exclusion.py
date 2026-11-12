"""主线2: 神经邻近区 CD8/T 细胞空间排斥指数。

定义
----
immune exclusion index (IEI) = log2[ mean(CD8 density z, ring 0-50um)
                                      / mean(CD8 density z, ring >400um) ]
负值 = 神经邻近区 CD8 更少 = 排斥。

三层统计
--------
A. 每样本: 秩空间偏 Spearman (CD8density vs distance, 残差化 log UMI / pct_mt)
   + Moran 谱随机化 (MSR) 零分布校准 (与基因级/读数级完全相同的置换框架);
   变体: 残差化额外加入 myCAFscore ("net of CAF" 的 CD8 排斥)。
B. 几何指数 IEI: 残差化 CD8 场的环带比值; 跨样本 Wilcoxon / sign test。
C. 介导分析: 每样本偏 Spearman (myCAFscore vs CD8density | UMI, MT),
   分近端 (<=100um) / 远端 (>200um) 两层; Fisher-z meta。
   外部锚定: GSE176078 NNLS 解卷积的 myCAF vs CD8_T spot 丰度偏相关 (如可用)。

输出
----
01_data/processed/immune_exclusion_per_sample.csv
03_results/tables/immune_exclusion_summary.csv

用法: python 02_scripts/run_immune_exclusion.py [--nsim 199]
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
from scipy import stats

sys.path.insert(0, r"D:/嗜神经浸润/00_pipeline")
sys.path.insert(0, r"D:/嗜神经浸润/02_scripts")

from nsgrad import config as CFG                      # noqa: E402
from run_msr_spatial import (                         # noqa: E402
    _rank_cols, build_knn_graph, moran_spectrum, ortho_cov,
)

DIR_H5 = CFG.DIR_PROC / "spot_level"
PROC = CFG.DIR_PROC
TAB = CFG.DIR_TAB
RING_ORDER = ["r1_0_50", "r2_50_100", "r3_100_200", "r4_200_400", "out_gt400"]
RING_D = {"r1_0_50": 25, "r2_50_100": 75, "r3_100_200": 150, "r4_200_400": 300,
          "out_gt400": 600}


def residualize(v, Q):
    """Gram-Schmidt 残差化: v - proj_Q(v)."""
    v = v.astype(float).copy()
    for j in range(Q.shape[1]):
        v -= Q[:, j] * (Q[:, j] @ v)
    return v


def msr_null(v, d, xy, Q, nsim, seed, bin_target=40):
    """对单个 spot 级场做 MSR: 返回 (rho_obs, mu, sd, u)."""
    n = len(v)
    R = v / (v.std() + 1e-12)
    rd = stats.rankdata(d)
    rd = (rd - rd.mean()) / (rd.std() + 1e-12)
    v_res = residualize(rd, Q)
    vv = float(v_res @ v_res)
    rtQ = R @ Q
    den = np.sqrt(max(n - float((rtQ ** 2).sum()), 1e-12) * vv)
    rho_obs = float(np.clip((R @ v_res) / den, -1, 1))

    B = build_knn_graph(xy)
    lam, V, S0 = moran_spectrum(B)
    order = np.argsort(lam)
    V_s = np.ascontiguousarray(V[:, order])
    lam_s = lam[order]
    del B, V, lam
    c = V_s.T @ R
    w = V_s.T @ v_res
    U = V_s.T @ Q

    nb = int(max(4, min(120, n // bin_target)))
    edges = np.linspace(0, n, nb + 1).astype(int)
    rng = np.random.default_rng(seed)
    rho_sim = np.empty(nsim)
    for r in range(nsim):
        cs = c.copy()
        for b in range(nb):
            s, e = edges[b], edges[b + 1]
            if e - s > 1:
                cs[s:e] = cs[s:e][rng.permutation(e - s)]
        cu = cs @ U
        denr = np.sqrt(max((cs @ cs) - float((cu ** 2).sum()), 1e-12) * vv)
        rho_sim[r] = np.clip((cs @ w) / denr, -1, 1)
    mu, sd = float(rho_sim.mean()), float(rho_sim.std(ddof=1))
    return rho_obs, mu, sd, (rho_obs - mu) / max(sd, 1e-4)


def process_sample(path: Path, nsim, seed):
    ad = sc.read_h5ad(path)
    obs = ad.obs
    keep = np.isfinite(obs["distance_to_nearest_nerve"].to_numpy(float))
    xy = np.column_stack([
        np.asarray(obs["x_um"].to_numpy(float)),
        np.asarray(obs["y_um"].to_numpy(float)),
    ])[keep]
    d = obs["distance_to_nearest_nerve"].to_numpy(float)[keep]
    ring = obs["ring_id"].astype(str).to_numpy()[keep]
    tc = obs["total_counts"].to_numpy(float)[keep]
    pm = obs["pct_mt"].to_numpy(float)[keep]
    cd8 = obs["CD8density"].to_numpy(float)[keep]
    caf = obs["myCAFscore"].to_numpy(float)[keep]
    tls = obs["TLSscore"].to_numpy(float)[keep] if "TLSscore" in obs else None
    n = int(keep.sum())

    # 秩空间 + 残差化
    rank = lambda a: stats.rankdata(a)
    Q2 = ortho_cov(np.log1p(tc), pm)                # UMI + MT
    Q3 = np.column_stack([Q2, rank(caf)])           # + CAF 含量
    # QR 正交化 Q3
    Q3, _ = np.linalg.qr(Q3)

    fields = {"CD8": rank(cd8)}
    if tls is not None:
        fields["TLS"] = rank(tls)
    out = {"sample": path.stem, "n_spot": n}
    msr_rows = []
    for nm, v in fields.items():
        v1 = residualize(v - v.mean(), Q2)
        v2 = residualize(v - v.mean(), Q3)          # net of CAF
        for tag, vv in [("", v1), ("_netcaf", v2)]:
            rho, mu, sd, u = msr_null(vv, d, xy, Q2, nsim, seed + hash(nm + tag) % 1000)
            out[f"rho{tag}_{nm}"] = rho
            out[f"u_msr{tag}_{nm}"] = u
        # 环带均值 (残差化场) -> IEI
        z = (v1 - v1.mean()) / (v1.std() + 1e-12)
        means = {}
        for rk in RING_ORDER:
            m = ring == rk
            if m.sum() >= 5:
                means[rk] = float(z[m].mean())
        out[f"ring_means_{nm}"] = means
        if "r1_0_50" in means and "out_gt400" in means:
            out[f"iei_{nm}"] = float(np.log2((means["r1_0_50"] + 0.25) /
                                             (means["out_gt400"] + 0.25)))
        # 单调性: 环带均值 vs 环中心距离
        if len(means) >= 3:
            ks = list(means)
            out[f"ring_trend_{nm}"] = float(stats.spearmanr(
                [RING_D[k] for k in ks], [means[k] for k in ks]).statistic)

    # 介导: myCAF vs CD8 偏相关 (| UMI, MT), 近端/远端分层
    for tag, msk in [("prox", (d <= 100)), ("dist", (d > 200))]:
        if msk.sum() >= 40:
            q = ortho_cov(np.log1p(tc[msk]), pm[msk])
            def pres(a, b):
                ra = residualize(rank(a) - rank(a).mean(), q)
                rb = residualize(rank(b) - rank(b).mean(), q)
                r = float(np.corrcoef(ra, rb)[0, 1])
                return r
            out[f"rho_mycaf_cd8_{tag}"] = pres(caf[msk], cd8[msk])
    del ad
    return out, msr_rows


def main(nsim: int):
    rows = []
    for i, p in enumerate(sorted(DIR_H5.glob("*.h5ad")), 1):
        try:
            o, _ = process_sample(p, nsim, i)
        except Exception as e:
            print(f"  [skip] {p.stem}: {type(e).__name__}: {e}")
            continue
        rows.append(o)
        print(f"  [{i}] {p.stem} IEI(CD8)={o.get('iei_CD8', float('nan')):+.2f} "
              f"rho={o.get('rho_CD8', float('nan')):+.3f} u={o.get('u_msr_CD8', float('nan')):+.1f}",
              flush=True)
    P = pd.DataFrame(rows)
    P.to_csv(PROC / "immune_exclusion_per_sample.csv", index=False)

    # ---- meta ----
    def meta(col):
        v = P[col].dropna()
        if len(v) < 10:
            return None
        t, p = stats.wilcoxon(v) if (v != 0).any() else (np.nan, np.nan)
        return dict(metric=col, n=len(v), median=float(v.median()),
                    mean=float(v.mean()), q1=float(v.quantile(.25)),
                    q3=float(v.quantile(.75)),
                    frac_neg=float((v < 0).mean()),
                    p_wilcoxon=float(p) if not np.isnan(p) else np.nan)

    metrics = [c for c in P.columns if c.startswith(("iei_", "rho", "u_msr", "ring_trend"))]
    Ms = [m for m in (meta(c) for c in metrics) if m]
    M = pd.DataFrame(Ms)
    # Stouffer meta of u_msr
    stou = []
    for c in [c for c in metrics if c.startswith("u_msr")]:
        u = P[c].dropna().to_numpy(float)
        if len(u) >= 10:
            stou.append(dict(metric=c, n=len(u),
                             stouffer_Z=float(u.sum() / np.sqrt(len(u))),
                             p_stouffer=float(2 * stats.norm.sf(abs(u.sum() / np.sqrt(len(u)))))))
    M = M.merge(pd.DataFrame(stou), on="metric", how="left") if stou else M
    M.to_csv(TAB / "immune_exclusion_summary.csv", index=False)
    print("\n", M.round(4).to_string(index=False))
    print(f"\nsaved: {TAB/'immune_exclusion_summary.csv'}")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--nsim", type=int, default=199)
    a = ap.parse_args()
    main(a.nsim)

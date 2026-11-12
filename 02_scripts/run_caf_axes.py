"""主线3: CAF 四轴 (myCAF / iCAF / apCAF / LRRC15+) 的空间分层与伪时间轨迹。

每 spot:
  1) 四轴 mean-z 评分 (主指标, 已在 run_axis_gradient 中经 MSR + 随机基线校准);
  2) NNLS 非负回归: min ||x_z - A w||, w >= 0, A 为四轴单位化基因面板向量
     (分解版, 用于轨迹图的连续构成)。

分析:
  a) 沿神经距离的环带分层 (NNLS 权重的 ring profile) —— 哪类贴神经、哪类外围;
  b) 每样本偏 Spearman (轴评分 vs 距离 | UMI, MT), Fisher-z meta + 样本级 bootstrap;
  c) diffusion pseudotime (神经 spot 为根, 复用 run_pseudotime 的 gauss 后端),
     轴评分 vs dpt 的偏 Spearman + 80% spot 重采样 x 100 bootstrap CI;
  d) 轨迹 bin 剖面 (dpt 十分位 x 各轴 + CD8/T) —— "神经 -> CAF 轴 -> 免疫读数"。

输出
----
01_data/processed/caf_axis_per_spot.csv     每 spot 四轴 mean-z + NNLS 权重 + dpt
01_data/processed/caf_axis_per_sample.csv   每样本梯度/轨迹 rho + bootstrap CI
01_data/processed/caf_axis_meta.csv         Fisher-z meta
01_data/processed/caf_axis_traj_bins.csv    轨迹 bin 剖面

用法: python 02_scripts/run_caf_axes.py [--nboot 100]
"""
from __future__ import annotations

import argparse
import importlib.util
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sparse
from scipy import stats
from scipy.optimize import nnls

sys.path.insert(0, r"D:/嗜神经浸润/00_pipeline")
sys.path.insert(0, r"D:/嗜神经浸润/02_scripts")

from nsgrad import config as CFG                      # noqa: E402
import niche_lib as NL                                # noqa: E402

DIR_H5 = CFG.DIR_PROC / "spot_level"
PROC = CFG.DIR_PROC
N_PCS = 15
RING_ORDER = ["r1_0_50", "r2_50_100", "r3_100_200", "r4_200_400", "out_gt400"]

AXES = {
    "myCAF": [g for g in NL.CELLTYPE_MARKERS["CAF_myofibroblastic"] if g != "LRRC15"],
    "iCAF": list(NL.CELLTYPE_MARKERS["CAF_inflammatory_dev"]),
    "apCAF": ["HLA-DRA", "HLA-DRB1", "HLA-DPA1", "HLA-DPB1", "CD74", "CIITA", "LAMP3", "CD70"],
    "LRRC15p": ["LRRC15"],
}
IMMUNE_READOUTS = {"CD8_T": list(NL.CELLTYPE_MARKERS["CD8_T"]),
                   "T_lymphocyte": list(NL.CELLTYPE_MARKERS["T_lymphocyte"])}


def _rank(a):
    return stats.rankdata(a)


def _ortho_cov(tc, pm):
    X = np.column_stack([np.log1p(tc), pm])
    X = X - X.mean(axis=0)
    Q, _ = np.linalg.qr(X)
    return Q


def _partial_rho(a, b, Q):
    ra, rb = _rank(a), _rank(b)
    ra = ra - ra.mean(); rb = rb - rb.mean()
    for j in range(Q.shape[1]):
        ra -= Q[:, j] * (Q[:, j] @ ra)
        rb -= Q[:, j] * (Q[:, j] @ rb)
    den = np.sqrt(float(ra @ ra) * float(rb @ rb))
    return float(np.clip((ra @ rb) / den, -1, 1)) if den > 1e-12 else np.nan


def process_sample(path: Path, nboot: int, seed: int):
    ad = sc.read_h5ad(path)
    obs = ad.obs
    keep = np.isfinite(obs["distance_to_nearest_nerve"].to_numpy(float))
    n = int(keep.sum())
    genes = np.asarray(ad.var_names)
    X = ad.X
    X = np.asarray(X.toarray() if sparse.issparse(X) else X)[keep]
    dist = obs["distance_to_nearest_nerve"].to_numpy(float)[keep]
    ring = obs["ring_id"].astype(str).to_numpy()[keep]
    tc = obs["total_counts"].to_numpy(float)[keep]
    pm = obs["pct_mt"].to_numpy(float)[keep]
    nerve = obs["nerve_spatial"].to_numpy(bool)[keep] if "nerve_spatial" in obs else np.zeros(n, bool)
    nsc = obs["nerve_score"].to_numpy(float)[keep] if "nerve_score" in obs else np.zeros(n)
    sample = path.stem

    # ---- 轴基因 z ----
    gidx = {g: i for i, g in enumerate(genes)}
    union = sorted({g for gs in AXES.values() for g in gs if g in gidx})
    lib = tc.copy(); lib[lib <= 0] = 1.0
    Xz = np.log1p(X[:, [gidx[g] for g in union]] / lib[:, None] * 1e4)
    Xz = (Xz - Xz.mean(axis=0)) / (Xz.std(axis=0) + 1e-12)
    uidx = {g: i for i, g in enumerate(union)}

    # mean-z 轴评分
    ax = {}
    for k, gs in AXES.items():
        cols = [uidx[g] for g in gs if g in uidx]
        ax[k] = Xz[:, cols].mean(axis=1) if cols else np.full(n, np.nan)

    # NNLS: A (G x 4), 列 = 轴面板单位向量
    A = np.zeros((len(union), len(AXES)))
    for k, (name, gs) in enumerate(AXES.items()):
        for g in gs:
            if g in uidx:
                A[uidx[g], k] = 1.0 / np.sqrt(max(1, len(gs)))
    W = np.zeros((n, len(AXES)))
    for i in range(n):
        try:
            W[i], _ = nnls(A, Xz[i])
        except Exception:
            W[i] = np.nan

    # ---- diffusion pseudotime (神经 spot 为根) ----
    try:
        ad2 = ad[keep].copy()
        ad2.X = sparse.csr_matrix(X)
        sc.pp.normalize_total(ad2, target_sum=1e4)
        sc.pp.log1p(ad2)
        try:
            sc.pp.highly_variable_genes(ad2, n_top_genes=2000, flavor="seurat")
            ad2 = ad2[:, ad2.var.highly_variable].copy()
        except Exception:
            pass
        sc.pp.pca(ad2, n_comps=min(N_PCS, ad2.n_vars - 1, n - 1),
                  svd_solver="arpack", random_state=seed)
        sc.pp.neighbors(ad2, n_neighbors=15, method="gauss", random_state=seed)
        sc.tl.diffmap(ad2, n_comps=min(15, ad2.n_vars - 1))
        if nerve.sum() >= 1:
            cand = np.where(nerve)[0]
            iroot = int(cand[np.argmax(nsc[cand])])
        else:
            iroot = int(np.argmin(dist))
        ad2.uns["iroot"] = iroot
        sc.tl.dpt(ad2)
        dpt = ad2.obs["dpt_pseudotime"].to_numpy(float)
    except Exception as e:
        print(f"  [dpt fail] {sample}: {type(e).__name__}")
        dpt = np.full(n, np.nan)
    del ad, ad2

    # ---- 统计 ----
    Q = _ortho_cov(tc, pm)
    row = {"sample": sample, "n_spot": n, "dpt_ok": int(np.isfinite(dpt).sum()),
           "iroot_nerve": bool(nerve.sum() >= 1)}
    for k in AXES:
        v = ax[k]
        ok = np.isfinite(v)
        row[f"rho_dist_{k}"] = _partial_rho(v[ok], dist[ok], Q[ok])
        row[f"n_meas_{k}"] = int(ok.sum())
        if np.isfinite(dpt).all() and ok.sum() > 100:
            row[f"rho_dpt_{k}"] = _partial_rho(v[ok], dpt[ok], Q[ok])
            # bootstrap: 80% spot 重采样
            rng = np.random.default_rng(seed + 31)
            bs = []
            ok_idx = np.where(ok)[0]
            for _ in range(nboot):
                sel = rng.choice(ok_idx, size=int(0.8 * len(ok_idx)), replace=False)
                bs.append(_partial_rho(v[sel], dpt[sel], Q[sel]))
            row[f"rho_dpt_{k}_lo"] = float(np.percentile(bs, 2.5))
            row[f"rho_dpt_{k}_hi"] = float(np.percentile(bs, 97.5))
    # 免疫读数 vs dpt (轨迹终点)
    for k, gs in IMMUNE_READOUTS.items():
        cols = [uidx[g] for g in gs if g in uidx]
        if cols:
            v = Xz[:, cols].mean(axis=1)
            if np.isfinite(dpt).all():
                row[f"rho_dpt_{k}"] = _partial_rho(v, dpt, Q)

    # ---- per-spot 输出 ----
    ps = pd.DataFrame({"sample": sample, "spot": np.arange(n), "dist": dist,
                       "ring": ring, "dpt": dpt})
    for k in AXES:
        ps[f"score_{k}"] = ax[k]
        ps[f"nnls_{k}"] = W[:, list(AXES).index(k)]
    for k, gs in IMMUNE_READOUTS.items():
        cols = [uidx[g] for g in gs if g in uidx]
        ps[f"score_{k}"] = Xz[:, cols].mean(axis=1) if cols else np.nan

    # ---- 轨迹 bin 剖面 (dpt 十分位) ----
    bins = []
    if np.isfinite(dpt).all():
        q = pd.qcut(dpt, 10, labels=False, duplicates="drop")
        for b in np.unique(q):
            m = q == b
            r = {"sample": sample, "bin": int(b), "n_spot": int(m.sum()),
                 "dpt_median": float(np.median(dpt[m])),
                 "dist_median": float(np.median(dist[m]))}
            for k in AXES:
                r[f"score_{k}"] = float(ax[k][m].mean()) if np.isfinite(ax[k][m]).any() else np.nan
                r[f"nnls_{k}"] = float(W[m, list(AXES).index(k)].mean())
            for k in IMMUNE_READOUTS:
                r[f"score_{k}"] = float(ps.loc[m, f"score_{k}"].mean())
            bins.append(r)
    return row, ps, pd.DataFrame(bins)


def main(nboot: int):
    rows, spots, bins = [], [], []
    for i, p in enumerate(sorted(DIR_H5.glob("*.h5ad")), 1):
        try:
            row, ps, bb = process_sample(p, nboot, i)
        except Exception as e:
            print(f"  [skip] {p.stem}: {type(e).__name__}: {e}")
            continue
        rows.append(row); spots.append(ps)
        if len(bb):
            bins.append(bb)
        print(f"  [{i}] {p.stem}: myCAF rho_dist={row.get('rho_dist_myCAF', float('nan')):+.3f} "
              f"rho_dpt={row.get('rho_dpt_myCAF', float('nan')):+.3f}", flush=True)

    R = pd.DataFrame(rows)
    R.to_csv(PROC / "caf_axis_per_sample.csv", index=False)
    pd.concat(spots, ignore_index=True).to_csv(PROC / "caf_axis_per_spot.csv", index=False)
    B = pd.concat(bins, ignore_index=True)
    B.to_csv(PROC / "caf_axis_traj_bins.csv", index=False)

    # ---- meta (Fisher-z) ----
    out = []
    for k in list(AXES) + list(IMMUNE_READOUTS):
        for pre in ("rho_dist_", "rho_dpt_"):
            col = pre + k
            v = R[col].dropna() if col in R else pd.Series(dtype=float)
            if len(v) < 10:
                continue
            z = np.arctanh(np.clip(v, -0.999, 0.999))
            t, p = stats.ttest_1samp(z, 0)
            # 样本级 bootstrap
            rng = np.random.default_rng(42)
            bm = []
            for _ in range(200):
                bs = rng.choice(z, size=len(z), replace=True)
                bm.append(np.tanh(bs.mean()))
            out.append(dict(axis=k, metric=pre.strip("_"), n=len(v),
                            rho_median=float(v.median()), rho_mean=float(np.tanh(z.mean())),
                            p_ttest=float(p),
                            boot_lo=float(np.percentile(bm, 2.5)),
                            boot_hi=float(np.percentile(bm, 97.5)),
                            frac_neg=float((v < 0).mean())))
    M = pd.DataFrame(out)
    M.to_csv(PROC / "caf_axis_meta.csv", index=False)
    print("\n", M.round(4).to_string(index=False))
    print(f"\nsaved: caf_axis_per_sample / per_spot / traj_bins / meta")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--nboot", type=int, default=100)
    a = ap.parse_args()
    main(a.nboot)

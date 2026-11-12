# -*- coding: utf-8 -*-
"""
方案 2-A 第三环：spot 水平细胞丰度与 NPP / 神经 spot 的共聚集与梯度。

分析（与 niche 读数同约定: 偏 Spearman 控制 log1p(total_counts) 与 pct_mt,
每样本一 rho, Fisher-z 合并, BH-FDR）:
  A. 细胞类型丰度 vs 到神经距离的梯度（负 rho = 神经近端富集）
  B. 神经 spot vs 非神经 spot 的各类型丰度差（Mann-Whitney per sample, 对数比合并）
  C. myCAF 丰度与两个 NPP 模块 spot 得分的共聚集（Spearman per sample, meta）
  D. 环带（0-50 … >400 µm）各类型丰度曲线

输出:
  01_data/processed/deconv_gradient_rho.csv     每样本各类型 rho
  01_data/processed/deconv_gradient_meta.csv    合并检验
  01_data/processed/deconv_nerve_vs_non.csv     神经 spot 丰度比较
  01_data/processed/deconv_module_coaggr.csv    模块-细胞共聚集
  01_data/processed/deconv_ring_profiles.csv    环带丰度
"""
from __future__ import annotations
import glob, time
from pathlib import Path
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse, stats

ROOT = Path(__file__).resolve().parents[1]
PROC = ROOT / "01_data/processed"

NEURAL_INTRINSIC = [
    "UCHL1", "TUBB3", "TUBA1A", "MAP2", "NEFL", "NEFM", "PLP1", "MPZ", "PMP22",
    "S100B", "NGFR", "SOX10", "SYP", "SCN9A", "SCN10A", "SCN7A", "PRX", "MBP",
    "TH", "DBH", "CHAT", "SLC18A3", "SLC5A7", "TRPV1", "CALCA", "CALCB", "VIP",
    "TAC1", "NPY", "SST", "PENK", "GAD1", "GAD2", "SLC6A1", "GAP43", "ELAVL4",
    "STMN2", "NRXN1", "SYT1",
]
PERINEURAL_STROMA = [
    "POSTN", "SPARC", "FN1", "CCN2", "THBS1", "THBS2", "COMP", "LRRC15",
    "COL1A1", "COL1A2", "COL3A1", "COL5A1", "COL5A2", "COL4A1", "COL4A2",
    "COL8A1", "COL10A1", "COL11A1", "COL12A1", "COL15A1", "COL18A1",
    "TAGLN", "ACTA2", "MYL9", "CNN1", "VIM", "TPM1", "FAP", "SULF1", "ASPN",
    "MMP1", "MMP2", "MMP9", "MMP11", "MMP13", "MMP14", "ADAM12", "LOX", "LOXL2",
    "TNC", "VCAN", "DCN", "LUM", "SERPINF1", "FSTL1", "CTHRC1", "SFRP1",
    "SERPINE1", "INHBA", "TGFB1", "PDGFRB", "PDGFRA", "RGS5", "MCAM",
    "CDH11", "NID1", "TGM2", "ANTXR1", "DKK3", "S100A4",
]
RING_ORDER = ["r1_0_50", "r2_50_100", "r3_100_200", "r4_200_400", "r5_gt400"]


def partial_spearman(x, y, tc, pm):
    """秩空间中控制 log1p(tc), pm 的偏 Spearman。"""
    rx = stats.rankdata(x); ry = stats.rankdata(y)
    rx = (rx - rx.mean()) / (rx.std() + 1e-12)
    ry = (ry - ry.mean()) / (ry.std() + 1e-12)
    cov = np.column_stack([np.log1p(tc) - np.log1p(tc).mean(), pm - pm.mean()])
    # Gram-Schmidt 正交化协变量
    Q = np.zeros_like(cov)
    for j in range(cov.shape[1]):
        v = cov[:, j].copy()
        for k in range(j):
            v -= Q[:, k] * (Q[:, k] @ v)
        nv = np.sqrt(v @ v)
        Q[:, j] = v / nv if nv > 1e-12 else 0.0
    def resid(a):
        out = a.copy().astype(float)
        for j in range(Q.shape[1]):
            out -= Q[:, j] * (Q[:, j] @ out)
        return out
    rxr, ryr = resid(rx), resid(ry)
    return float(rxr @ ryr / np.sqrt((rxr @ rxr) * (ryr @ ryr) + 1e-18))


def fisher_z_meta(rhos, n_eff=None):
    r = np.asarray(rhos, float)
    r = r[np.isfinite(r)]
    if len(r) == 0:
        return dict(k=0, rho_mean=np.nan, z=np.nan, p=np.nan)
    z = np.arctanh(np.clip(r, -0.999, 0.999))
    if n_eff is not None:
        w = np.asarray(n_eff, float)[np.isfinite(rhos)] - 3.0
    else:
        w = np.ones_like(z)
    zbar = (z * w).sum() / w.sum()
    se = 1.0 / np.sqrt(w.sum())
    from scipy.stats import norm
    p = 2 * norm.sf(abs(zbar / se))
    return dict(k=len(r), rho_mean=float(np.tanh(zbar)), z=float(zbar), p=float(p))


def bh(p):
    p = np.asarray(p, float)
    n = len(p); order = np.argsort(p)
    q = p[order] * n / (np.arange(n) + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    out = np.empty(n); out[order] = np.clip(q, 0, 1)
    return out


def main():
    t0 = time.time()
    sfx = "_nomodule" if "--no-module" in __import__("sys").argv else ""
    dv = pd.read_csv(PROC / f"deconv_spots{sfx}.csv")
    types = [c[5:] for c in dv.columns if c.startswith("frac_")]
    print(f"spots={len(dv)}, samples={dv['sample'].nunique()}, types={types}")

    # ---- A: 丰度-距离梯度（偏 Spearman）----
    rows, meta_rows = [], []
    for s, g in dv.groupby("sample"):
        keep = np.isfinite(g["distance_to_nearest_nerve"].to_numpy(float))
        g = g[keep]
        d = g["distance_to_nearest_nerve"].to_numpy(float)
        tc = g["total_counts"].to_numpy(float)
        pm = g["pct_mt"].to_numpy(float)
        for t in types:
            y = g[f"frac_{t}"].to_numpy(float)
            rho = partial_spearman(y, d, tc, pm)
            rows.append({"sample": s, "celltype": t, "rho": rho, "n": keep.sum()})
    grad = pd.DataFrame(rows)
    grad.to_csv(PROC / f"deconv_gradient_rho{sfx}.csv", index=False)

    for t in types:
        sub = grad[grad.celltype == t]
        m = fisher_z_meta(sub["rho"].to_numpy(), sub["n"].to_numpy())
        meta_rows.append({"celltype": t, **m})
    meta = pd.DataFrame(meta_rows)
    meta["fdr_bh"] = bh(meta["p"].to_numpy())
    meta = meta.sort_values("rho_mean")
    meta.to_csv(PROC / f"deconv_gradient_meta{sfx}.csv", index=False)
    print("\n==== A. 细胞丰度 ~ 距离(负=神经近端富集) ====")
    print(meta.round(3).to_string(index=False))

    # ---- B: 神经 spot vs 非神经 spot 丰度 ----
    b_rows = []
    for s, g in dv.groupby("sample"):
        g = g[g["nerve_spatial"].notna()]
        nv = g["nerve_spatial"].to_numpy(int) == 1
        if nv.sum() < 3 or (~nv).sum() < 3:
            continue
        for t in types:
            y = g[f"frac_{t}"].to_numpy(float)
            u, p = stats.mannwhitneyu(y[nv], y[~nv], alternative="two-sided")
            mn = float(y[nv].mean()); mo = float(y[~nv].mean())
            b_rows.append({
                "sample": s, "celltype": t, "n_nerve": int(nv.sum()),
                "mean_nerve": mn, "mean_other": mo,
                "log2fc": float(np.log2((mn + 1e-3) / (mo + 1e-3))),
                "p_mwu": p,
            })
    nerve = pd.DataFrame(b_rows)
    nerve.to_csv(PROC / f"deconv_nerve_vs_non{sfx}.csv", index=False)
    b_meta = []
    for t in types:
        sub = nerve[nerve.celltype == t].dropna()
        if len(sub) < 5:
            continue
        from scipy.stats import combine_pvalues
        stat, p = combine_pvalues(sub["p_mwu"].to_numpy(), method="stouffer")
        b_meta.append({"celltype": t, "k": len(sub),
                       "mean_log2fc": float(sub["log2fc"].mean()),
                       "median_log2fc": float(sub["log2fc"].median()),
                       "p_stouffer": p})
    b_meta = pd.DataFrame(b_meta)
    if len(b_meta):
        b_meta["fdr_bh"] = bh(b_meta["p_stouffer"].to_numpy())
    b_meta.to_csv(PROC / f"deconv_nerve_vs_non_meta{sfx}.csv", index=False)
    print("\n==== B. 神经 spot vs 其他 spot ====")
    print(b_meta.round(3).to_string(index=False))

    # ---- C: 模块得分 ~ myCAF(及其他类型) 共聚集 ----
    spot_scores = []
    for f in sorted(glob.glob(str(PROC / "spot_level/*.h5ad"))):
        ad = sc.read_h5ad(f)
        keep = ad.obs["in_tissue"].astype(bool).values if "in_tissue" in ad.obs else np.ones(ad.n_obs, bool)
        X = ad.X.toarray() if sparse.issparse(ad.X) else np.asarray(ad.X)
        X = X[keep]; obs = ad.obs.iloc[np.where(keep)[0]]
        lib = X.sum(axis=1); lib = np.where(lib <= 0, 1, lib)
        L = np.log1p(X / lib[:, None] * 1e4)
        gi = {g: j for j, g in enumerate(ad.var_names)}
        def score(genes):
            cols = [gi[g] for g in genes if g in gi]
            if len(cols) < max(3, 0.15 * len(genes)):
                return None
            Z = (L[:, cols] - L[:, cols].mean(0)) / (L[:, cols].std(0) + 1e-12)
            return Z.mean(axis=1)
        s_ni = score(NEURAL_INTRINSIC)
        s_ps = score(PERINEURAL_STROMA)
        if s_ni is None or s_ps is None:
            continue
        spot_scores.append(pd.DataFrame({
            "sample": ad.uns.get("sample_id", Path(f).stem),
            "barcode": obs.index, "neural_intrinsic": s_ni, "perineural_stroma": s_ps,
            "total_counts": obs["total_counts"].to_numpy(float),
            "pct_mt": obs["pct_mt"].to_numpy(float),
        }))
        del ad, X, L
    scores = pd.concat(spot_scores, ignore_index=True)
    print(f"\nmodule scores: {len(scores)} spots")

    c_rows = []
    dv_key = dv[["sample", "barcode"] + [f"frac_{t}" for t in types]]
    merged = scores.merge(dv_key, on=["sample", "barcode"], how="inner")
    print(f"merged: {len(merged)} spots")
    for s, g in merged.groupby("sample"):
        tc = g["total_counts"].to_numpy(float); pm = g["pct_mt"].to_numpy(float)
        for t in types:
            y = g[f"frac_{t}"].to_numpy(float)
            c_rows.append({"sample": s, "celltype": t,
                           "rho_neural": partial_spearman(y, g["neural_intrinsic"].to_numpy(float), tc, pm),
                           "rho_stroma": partial_spearman(y, g["perineural_stroma"].to_numpy(float), tc, pm),
                           "n": len(g)})
    coag = pd.DataFrame(c_rows)
    coag.to_csv(PROC / f"deconv_module_coaggr{sfx}.csv", index=False)
    c_meta = []
    for t in types:
        sub = coag[coag.celltype == t].dropna()
        m1 = fisher_z_meta(sub["rho_neural"].to_numpy(), sub["n"].to_numpy())
        m2 = fisher_z_meta(sub["rho_stroma"].to_numpy(), sub["n"].to_numpy())
        c_meta.append({"celltype": t, "k": m1["k"],
                       "rho_neural": m1["rho_mean"], "p_neural": m1["p"],
                       "rho_stroma": m2["rho_mean"], "p_stroma": m2["p"]})
    c_meta = pd.DataFrame(c_meta)
    c_meta["fdr_neural"] = bh(c_meta["p_neural"].to_numpy())
    c_meta["fdr_stroma"] = bh(c_meta["p_stroma"].to_numpy())
    c_meta.to_csv(PROC / f"deconv_module_coaggr_meta{sfx}.csv", index=False)
    print("\n==== C. 模块-细胞共聚集 (偏 Spearman, Fisher-z) ====")
    print(c_meta.round(3).to_string(index=False))

    # ---- D: 环带丰度 ----
    dv["ring_id"] = dv["ring_id"].astype(str)
    ring = dv[dv["ring_id"].isin(RING_ORDER)].groupby(["ring_id", "sample"])[
        [f"frac_{t}" for t in types]].mean().reset_index()
    ring.to_csv(PROC / f"deconv_ring_profiles{sfx}.csv", index=False)
    print(f"\nring profile samples: {ring['sample'].nunique()}")
    print(f"\ntotal {time.time()-t0:.0f}s")


if __name__ == "__main__":
    main()

# -*- coding: utf-8 -*-
"""
方案 2-A 第二环：把 GSE176078 scRNA 细胞类型参考锚定到 30 例 Visium spot。

算法（SPOTlight 式 seeded 解卷积的 Python 复刻）：
  1) 参考: GSE176078 各细胞类型的 pseudobulk mean CP10K（seeded profiles，不做自由 NMF）
  2) 基因面板: 每类型 top 标记基因 + NPP 两个模块全部基因
  3) 每 spot: 在 sqrt 方差稳定尺度上做 NNLS
     y_sqrt ≈ Σ_t w_t * profile_t_sqrt,  w >= 0
  4) 丰度 = w / Σw

输出:
  01_data/processed/deconv_spots.csv        每 spot 各类型丰度 + 空间/神经注释
  01_data/processed/deconv_reference.csv    参考面板(基因 × 类型)
  01_data/processed/deconv_qc.csv           每 spot 重构误差 + 样本级 QC
"""
from __future__ import annotations
import time, glob
from pathlib import Path
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse
from scipy.optimize import nnls

ROOT = Path(__file__).resolve().parents[1]
PROC = ROOT / "01_data/processed"
DIR_H5 = PROC / "spot_level"

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
MODULE_GENES = set(NEURAL_INTRINSIC) | set(PERINEURAL_STROMA)

MIN_CELLS_TYPE = 500      # 参考类型最低细胞数
TOP_MARKERS_PER_TYPE = 60
EXCLUDE_MODULE = "--no-module" in __import__("sys").argv   # 敏感性: 面板剔除 NPP 模块基因


def build_reference():
    long = pd.read_csv(PROC / "scrna_cell_origin_long.csv")
    piv_c = long.pivot_table(index="gene", columns="celltype", values="mean_cp10k")
    piv_f = long.pivot_table(index="gene", columns="celltype", values="frac_expr")
    ncells = long.drop_duplicates("celltype").set_index("celltype")["n_cells"]
    keep_types = [t for t in piv_c.columns if ncells.get(t, 0) >= MIN_CELLS_TYPE]
    piv_c, piv_f = piv_c[keep_types], piv_f[keep_types]
    # 标记基因: 每类型按 frac_expr × 相对特异性排序
    panel = set() if EXCLUDE_MODULE else set(MODULE_GENES)
    for t in keep_types:
        spec = piv_f[t] / (piv_f.max(axis=1) + 1e-9)
        score = (piv_f[t] * spec).sort_values(ascending=False)
        cand = [g for g in score.index
                if piv_f.loc[g, t] >= 0.10
                and not (EXCLUDE_MODULE and g in MODULE_GENES)]
        panel.update(cand[:TOP_MARKERS_PER_TYPE])
    panel = {g for g in panel if g in piv_c.index}   # 只保留 scRNA 可测基因
    ref = piv_c.loc[sorted(panel)]
    return ref, piv_f


def main():
    t0 = time.time()
    ref, piv_f = build_reference()
    print(f"reference: {ref.shape[0]} genes x {ref.shape[1]} types; "
          f"types = {list(ref.columns)}")
    sfx = "_nomodule" if EXCLUDE_MODULE else ""
    ref.to_csv(PROC / f"deconv_reference{sfx}.csv")

    R = np.sqrt(ref.values)              # genes x types (sqrt VST)
    gnames_ref = list(ref.index)

    spot_rows, qc_rows = [], []
    files = sorted(glob.glob(str(DIR_H5 / "*.h5ad")))
    print(f"{len(files)} Visium samples")
    for i, f in enumerate(files, 1):
        ad = sc.read_h5ad(f)
        obs = ad.obs
        keep = obs["in_tissue"].astype(bool).values if "in_tissue" in obs else np.ones(ad.n_obs, bool)
        X = ad.X.toarray() if sparse.issparse(ad.X) else np.asarray(ad.X)
        X = X[keep]
        obs = obs.iloc[np.where(keep)[0]]
        gnames = list(ad.var_names)
        lib = X.sum(axis=1)
        lib = np.where(lib <= 0, 1.0, lib)
        cp10k = X / lib[:, None] * 1e4
        # 基因对齐
        idx = {g: j for j, g in enumerate(gnames)}
        rows_g = [idx[g] for g in gnames_ref if g in idx]
        genes_used = [g for g in gnames_ref if g in idx]
        Y = np.sqrt(cp10k[:, rows_g])    # spots x panel (sqrt VST)
        A = R[[gnames_ref.index(g) for g in genes_used], :]

        W = np.zeros((Y.shape[0], A.shape[1]))
        err = np.zeros(Y.shape[0])
        for s in range(Y.shape[0]):
            w, rnorm = nnls(A, Y[s])
            W[s] = w
            err[s] = rnorm / (np.linalg.norm(Y[s]) + 1e-12)
        Wsum = W.sum(axis=1, keepdims=True)
        P = np.divide(W, Wsum, out=np.zeros_like(W), where=Wsum > 0)

        base = pd.DataFrame({
            "sample": ad.uns.get("sample_id", Path(f).stem),
            "barcode": obs.index,
            "x_um": obs["x_um"].to_numpy(float),
            "y_um": obs["y_um"].to_numpy(float),
            "distance_to_nearest_nerve": obs["distance_to_nearest_nerve"].to_numpy(float),
            "nerve_spatial": obs["nerve_spatial"].to_numpy(int) if "nerve_spatial" in obs else np.nan,
            "ring_id": obs["ring_id"].to_numpy(str) if "ring_id" in obs else np.nan,
            "total_counts": obs["total_counts"].to_numpy(float),
            "pct_mt": obs["pct_mt"].to_numpy(float),
        })
        frac = pd.DataFrame(P, columns=[f"frac_{t}" for t in ref.columns], index=base.index)
        out = pd.concat([base.reset_index(drop=True), frac.reset_index(drop=True)], axis=1)
        spot_rows.append(out)
        qc_rows.append({
            "sample": base["sample"].iloc[0],
            "n_spots": Y.shape[0],
            "panel_genes_used": len(genes_used),
            "median_rel_err": float(np.median(err)),
            "frac_spots_solved": float((Wsum.ravel() > 0).mean()),
        })
        print(f"  [{i}/{len(files)}] {base['sample'].iloc[0]}: {Y.shape[0]} spots, "
              f"rel_err={np.median(err):.3f}, {time.time()-t0:.0f}s")
        del ad, X, cp10k

    pd.concat(spot_rows, ignore_index=True).to_csv(PROC / f"deconv_spots{sfx}.csv", index=False)
    pd.DataFrame(qc_rows).to_csv(PROC / f"deconv_qc{sfx}.csv", index=False)
    print(f"done in {time.time()-t0:.0f}s")


if __name__ == "__main__":
    main()

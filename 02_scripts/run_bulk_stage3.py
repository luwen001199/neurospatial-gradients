# -*- coding: utf-8 -*-
"""Stage 3 —— 单细胞 myCAF 锚定（GSE176078, Wu et al. 2021, 26 例）。

原方案写法是「bulk 高低组 ↔ 外部 scRNA myCAF 比例」。两套样本集不可逐例对应，
故正确操作化为：在 scRNA 内部做 pseudobulk，检验「模块评分 ↔ myCAF 细胞占比」
—— 即签名是否能读出 myCAF 含量。这是签名有效性验证，不是组间比较。

输出: 03_results/tables/bulk_extrapolation/stage3_scrna_anchor.csv
"""
from __future__ import annotations

import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse, stats
from scipy.io import mmread

sys.path.insert(0, r"D:/嗜神经浸润/02_scripts")
import bulk_lib as BL  # noqa: E402

SRC = BL.ROOT / "01_data" / "external" / "gse176078"
OUT = BL.TAB

CAF_MY = "CAFs myCAF-like"
CAF_ICA = "CAFs MSC iCAF-like"


def load_sample(d: Path):
    # scipy fast_matrix_market 在 Windows 上不支持非 ASCII 路径 → 必须用文件句柄
    with open(d / "count_matrix_sparse.mtx", "rb") as fh:
        ad = mmread(fh)
    X = sparse.csc_matrix(ad)                                          # genes × cells
    genes = pd.read_csv(d / "count_matrix_genes.tsv", header=None)[0].astype(str).tolist()
    meta = pd.read_csv(d / "metadata.csv", index_col=0)
    meta.index = meta.index.astype(str)
    bar = pd.read_csv(d / "count_matrix_barcodes.tsv", header=None)[0].astype(str).tolist()
    if set(bar) != set(meta.index):
        keep = [b for b in bar if b in set(meta.index)]
        idx = {b: i for i, b in enumerate(bar)}
        X = X[:, [idx[b] for b in keep]]
        meta = meta.loc[keep]
        bar = keep
    else:
        meta = meta.loc[bar]
    return X, genes, meta


def main():
    mods = BL.load_modules()
    want = sorted({g for k, v in mods.items() for g in v})
    dirs = sorted([p for p in SRC.iterdir()
                   if p.is_dir() and (p / "count_matrix_sparse.mtx").exists()])
    print(f"样本 {len(dirs)} 个")

    rows = []
    t0 = time.time()
    for d in dirs:
        X, genes, meta = load_sample(d)
        gi = {g: i for i, g in enumerate(genes)}
        rows_idx = [gi[g] for g in want if g in gi]
        rows_g = [g for g in want if g in gi]
        total = float(X.sum())
        if total <= 0:
            continue
        sub = X[rows_idx, :]
        counts = np.asarray(sub.sum(axis=1)).ravel()
        cp10k = np.log1p(counts / total * 1e4)
        ct = meta["celltype_minor"].astype(str)
        n_cells = len(meta)
        n_my = int((ct == CAF_MY).sum())
        n_ica = int((ct == CAF_ICA).sum())
        n_caf = int(ct.str.startswith("CAFs").sum())
        sub_t = meta.get("subtype", pd.Series("NA", index=meta.index)).astype(str).iloc[0]
        rows.append({"sample": d.name, "subtype": sub_t, "n_cells": n_cells,
                     "n_myCAF": n_my, "n_iCAF": n_ica, "n_CAF": n_caf,
                     "frac_myCAF_of_all": n_my / n_cells,
                     "frac_myCAF_of_CAF": (n_my / n_caf) if n_caf else np.nan,
                     "frac_CAF_of_all": n_caf / n_cells,
                     **{g: v for g, v in zip(rows_g, cp10k)}})

    df = pd.DataFrame(rows).set_index("sample")
    print(f"读取完成 {time.time()-t0:.0f}s, 矩阵 {df.shape}")

    # pseudobulk 模块评分（只对能解析的基因）
    expr = df[[g for g in want if g in df.columns]].copy()
    sc = BL.score_all(expr, mods)
    out = pd.concat([df[["subtype", "n_cells", "frac_myCAF_of_all", "frac_myCAF_of_CAF",
                         "frac_CAF_of_all"]], sc], axis=1)

    print("\n模块评分 vs myCAF 占比 (Spearman)")
    stat_rows = []
    for m in ["Perineural_stroma", "Neural_intrinsic", "CAF_ECM", "LRRC15_myCAF", "NPP58"]:
        for tgt in ["frac_myCAF_of_all", "frac_myCAF_of_CAF", "frac_CAF_of_all"]:
            if m not in out or tgt not in out:
                continue
            d = out[[m, tgt]].dropna()
            if len(d) < 8:
                continue
            rho, p = stats.spearmanr(d[m], d[tgt])
            print(f"  {m:20s} ~ {tgt:20s} ρ={rho:+.3f}  P={p:.3g}  n={len(d)}")
            stat_rows.append({"module": m, "target": tgt, "n": len(d),
                              "rho": round(rho, 4), "p": p})

    # TNBC 亚组
    tn = out["subtype"].astype(str).str.upper().str.contains("TNBC")
    if tn.sum() >= 8:
        print(f"\nTNBC 亚组 n={int(tn.sum())}")
        for m in ["Perineural_stroma", "Neural_intrinsic", "CAF_ECM"]:
            d = out.loc[tn, [m, "frac_myCAF_of_all"]].dropna()
            if len(d) >= 8:
                rho, p = stats.spearmanr(d[m], d["frac_myCAF_of_all"])
                print(f"  {m:20s} ρ={rho:+.3f}  P={p:.3g}  n={len(d)}")
                stat_rows.append({"module": m, "target": "frac_myCAF_of_all_TNBC",
                                  "n": len(d), "rho": round(rho, 4), "p": p})

    out.to_csv(OUT / "stage3_scrna_per_sample.csv")
    pd.DataFrame(stat_rows).to_csv(OUT / "stage3_scrna_anchor.csv", index=False)
    print(f"\n输出 -> {OUT}")


if __name__ == "__main__":
    main()

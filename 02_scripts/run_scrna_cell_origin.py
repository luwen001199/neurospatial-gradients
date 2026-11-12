# -*- coding: utf-8 -*-
"""
方案 2-A 第一环：公共 TNBC scRNA (GSE176078, Wu et al. 2021) 锚定
NPP 两个模块基因（Neural_intrinsic 39 / Perineural_stroma 60）的细胞来源。

输出:
  01_data/processed/scrna_cell_origin_long.csv   gene x celltype 的 frac_expr / mean_cp10k
  01_data/processed/scrna_samples_meta.csv       每样本细胞构成
  01_data/processed/scrna_module_origin.csv      模块基因的来源归属表（主结果）
  05_reports/... 打印汇总
"""
from __future__ import annotations
import sys, time, json
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.io import mmread
from scipy import sparse

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "01_data/external/gse176078"
PROC = ROOT / "01_data/processed"

# ---- 模块基因（与 run_npp_split.py 严格一致）----
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

# 细胞类型层级: major + CAF 细分(用 celltype_minor)
CAF_SPLIT = {
    "CAFs myCAF-like": "CAF_myCAF",
    "CAFs MSC iCAF-like": "CAF_iCAF",
    "CAFs Transitioning": "CAF_transitioning",
}

def load_sample(d: Path):
    # 注意: scipy fast_matrix_market 在 Windows 上不支持非 ASCII 路径,
    # 必须用文件句柄绕过(本项目根目录含中文)。
    with open(d / "count_matrix_sparse.mtx", "rb") as fh:
        ad = mmread(fh)
    X = sparse.csc_matrix(ad)          # genes x cells
    genes = pd.read_csv(d / "count_matrix_genes.tsv", header=None)[0].astype(str).tolist()
    meta = pd.read_csv(d / "metadata.csv", index_col=0)
    meta.index = meta.index.astype(str)
    bar = pd.read_csv(d / "count_matrix_barcodes.tsv", header=None)[0].astype(str).tolist()
    # barcode 对齐
    if set(bar) != set(meta.index):
        keep = [b for b in bar if b in set(meta.index)]
        idx = {b: i for i, b in enumerate(bar)}
        X = X[:, [idx[b] for b in keep]]
        meta = meta.loc[keep]
        bar = keep
    else:
        meta = meta.loc[bar]
    return X, genes, meta

def refine_type(row):
    mj = str(row["celltype_major"])
    mn = str(row["celltype_minor"])
    if mj == "CAFs" and mn in CAF_SPLIT:
        return CAF_SPLIT[mn]
    return mj

def main():
    t0 = time.time()
    dirs = sorted([p for p in SRC.iterdir() if p.is_dir() and (p / "count_matrix_sparse.mtx").exists()])
    print(f"found {len(dirs)} samples: {[d.name for d in dirs]}")

    # 全基因级累计: type -> (gene -> count_sum, gene -> n_cells_expressing)
    store = {}          # type -> dict gene_symbol -> [sum_counts]
    store_cells = {}    # type -> int
    store_pos = {}      # type -> dict gene -> n_cells_expressing (全基因组)
    samples_meta = []

    for d in dirs:
        X, genes, meta = load_sample(d)
        types = meta.apply(refine_type, axis=1).values

        for t in np.unique(types):
            mask = types == t
            Xt = X[:, mask]
            s = np.asarray(Xt.sum(axis=1)).ravel()
            pos = np.asarray((Xt > 0).sum(axis=1)).ravel()
            n = int(mask.sum())
            if t not in store:
                store[t] = {}
                store_cells[t] = 0
                store_pos[t] = {}
            store_cells[t] += n
            for g, v in zip(genes, s):
                store[t][g] = store[t].get(g, 0.0) + float(v)
            for g, p in zip(genes, pos):
                store_pos[t][g] = store_pos[t].get(g, 0) + int(p)

        sub_meta = pd.Series(types).value_counts()
        samples_meta.append({
            "sample": d.name,
            "subtype": meta["subtype"].iloc[0],
            "n_cells": X.shape[1],
            "median_genes_per_cell": float(meta["nFeature_RNA"].median()),
        })
        print(f"  {d.name}: {X.shape[1]} cells, {len(sub_meta)} types, {time.time()-t0:.0f}s")
        del X

    # ---- 汇总 ----
    types = sorted(store.keys())
    all_genes = sorted(set().union(*[set(store[t].keys()) for t in types]))
    counts_df = pd.DataFrame({t: [store[t].get(g, 0.0) for g in all_genes] for t in types},
                             index=all_genes)
    cells_s = pd.Series(store_cells, name="n_cells").reindex(types)
    cp10k = counts_df.div(cells_s.values, axis=1) * 1e4
    frac = pd.DataFrame({t: [store_pos[t].get(g, 0) / store_cells[t] for g in all_genes] for t in types},
                        index=all_genes)
    frac.index.name = "gene"

    long_rows = []
    for g in all_genes:
        for t in types:
            long_rows.append({
                "gene": g, "celltype": t, "n_cells": int(cells_s[t]),
                "frac_expr": frac.loc[g, t], "mean_cp10k": cp10k.loc[g, t],
            })
    pd.DataFrame(long_rows).to_csv(PROC / "scrna_cell_origin_long.csv", index=False)

    # ---- 模块基因归属 ----
    out = []
    for gene, module in [(g, "Neural_intrinsic") for g in NEURAL_INTRINSIC] + \
                        [(g, "Perineural_stroma") for g in PERINEURAL_STROMA]:
        if gene not in counts_df.index:
            out.append({"gene": gene, "module": module, "status": "not_measured"})
            continue
        f = frac.loc[gene]; c = cp10k.loc[gene]
        top = f.sort_values(ascending=False)
        # 富集: 模块基因在该类型表达细胞占比 / 全部类型平均占比
        base = top.mean()
        out.append({
            "gene": gene, "module": module, "status": "measured",
            "top_frac_type": top.index[0], "top_frac": float(top.iloc[0]),
            "2nd_frac_type": top.index[1], "2nd_frac": float(top.iloc[1]) if len(top) > 1 else np.nan,
            "frac_myCAF": float(f.get("CAF_myCAF", np.nan)),
            "frac_iCAF": float(f.get("CAF_iCAF", np.nan)),
            "frac_CancerEpithelial": float(f.get("Cancer_Epithelial", np.nan)),
            "frac_Tcells": float(f.get("T-cells", np.nan)),
            "frac_Myeloid": float(f.get("Myeloid", np.nan)),
            "frac_Endothelial": float(f.get("Endothelial", np.nan)),
            "frac_PVL": float(f.get("PVL", np.nan)),
            "mean_cp10k_top": float(c[top.index[0]]),
            "enrichment_top": float(top.iloc[0] / (base + 1e-9)),
            "total_frac_expr": float(f.mean()),
        })
    mod = pd.DataFrame(out)
    mod.to_csv(PROC / "scrna_module_origin.csv", index=False)

    sm = pd.DataFrame(samples_meta)
    sm.to_csv(PROC / "scrna_samples_meta.csv", index=False)

    # ---- 控制台汇总 ----
    print("\n===== 细胞类型构成 =====")
    print(cells_s.to_string())
    print("\n===== Neural_intrinsic (39) =====")
    for r in mod[mod.module == "Neural_intrinsic"].itertuples():
        if r.status == "measured":
            print(f"  {r.gene:10s} top={r.top_frac_type:18s} frac={r.top_frac:.3f} "
                  f"(myCAF {r.frac_myCAF:.3f}, iCAF {r.frac_iCAF:.3f})")
        else:
            print(f"  {r.gene:10s} NOT MEASURED")
    print("\n===== Perineural_stroma (60) =====")
    for r in mod[mod.module == "Perineural_stroma"].itertuples():
        if r.status == "measured":
            print(f"  {r.gene:10s} top={r.top_frac_type:18s} frac={r.top_frac:.3f}")
        else:
            print(f"  {r.gene:10s} NOT MEASURED")
    print(f"\nTotal time: {time.time()-t0:.0f}s")

if __name__ == "__main__":
    main()

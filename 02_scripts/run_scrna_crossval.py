"""GSE75688 (Chung et al. 2017, Nat Commun) 交叉验证: NPP 模块基因的细胞来源归属.

第二锚定集, 独立于 GSE176078. 输入:
  01_data/external/gse75688/TPM.txt.gz          基因 x 样本(单细胞+混样) TPM
  01_data/external/gse75688/sample_info.txt.gz  每列样本的细胞类型注释
输出:
  01_data/processed/scrna_crossval_long.csv     gene x celltype -> frac_expr / mean_TPM
  01_data/processed/scrna_crossval_summary.csv  模块级别 A/B/C 归属复现对照
用法: python 02_scripts/run_scrna_crossval.py
"""
from __future__ import annotations

import gzip
import importlib.util
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path("D:/嗜神经浸润")
EXT = ROOT / "01_data/external/gse75688"
PROC = ROOT / "01_data/processed"

# 模块基因从主分析导入, 保证与 NPP 拆分严格一致
_spec = importlib.util.spec_from_file_location("npp_split", ROOT / "02_scripts/run_npp_split.py")
_npp = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_npp)
MODULES = _npp.MODULES
NEURAL_INTRINSIC = MODULES["Neural_intrinsic"]
PERINEURAL_STROMA = MODULES["Perineural_stroma"]

# 关键列名候选 (Chung 2017 sample_info: characteristics_ch1.9 -> cell type 等, 自动探测)
CELLTYPE_KEYS = ["cell_type", "celltype", "characteristics", "Cell Type", "type"]
SAMPLE_KEYS = ["sample", "Sample", "title", "ID", "sample_name"]


def load_sample_info(path: Path) -> pd.DataFrame:
    with gzip.open(path, "rt") as fh:
        info = pd.read_csv(fh, sep="\t", index_col=None)
    # 找细胞类型列
    ct_col = None
    for c in info.columns:
        lc = str(c).lower()
        if any(k in lc for k in ("cell", "type", "characteristics")):
            vals = info[c].astype(str).str.lower()
            if vals.nunique() < 40 and any(
                k in v for v in vals.unique() for k in
                ("cancer", "t cell", "t-cell", "fibro", "endo", "immune", "myeloid", "b cell", "normal")
            ):
                ct_col = c
                break
    if ct_col is None:
        # GSE75688 简化 sample_info: sample/type/index/index2/index3, 细胞类型在 index3
        _cands = []
        for c in info.columns:
            vals = info[c].astype(str)
            if set(vals.unique()) <= {"Tumor", "Bcell", "Tcell", "Myeloid", "Stromal", "Immune"}:
                _cands.append((info[c].nunique(), c))
        if _cands:  # 取分类最细的一列 (index3: Tumor/Bcell/Tcell/Myeloid/Stromal)
            ct_col = sorted(_cands)[-1][1]
    if ct_col is None:
        raise SystemExit(f"cannot find cell-type column in {path}; columns={info.columns.tolist()}")
    # 找样本名列 (与 TPM 列可对上)
    sam_col = None
    for c in info.columns:
        if info[c].astype(str).isin([p for p in _tpm_cols]).sum() > len(info) * 0.5:
            sam_col = c
            break
    if sam_col is None:  # 退化为行号/顺序对齐
        info = info.rename(columns={ct_col: "celltype"})
        info["sample"] = _tpm_cols[: len(info)]
    else:
        info = info.rename(columns={sam_col: "sample", ct_col: "celltype"})
    return info[["sample", "celltype"]]


if __name__ == "__main__":
    # 1) 读 TPM (基因 x 样本)
    with gzip.open(EXT / "TPM.txt.gz", "rt") as fh:
        tpm = pd.read_csv(fh, sep="\t", index_col=0)
    # 去掉 gene_id 行前缀: index 可能是 Ensembl, gene_name 列提供符号
    if "gene_name" in tpm.columns:
        tpm = tpm.set_index("gene_name")
    drop_meta = [c for c in ("gene_id", "gene_type") if c in tpm.columns]
    tpm = tpm.drop(columns=[c for c in drop_meta if c in tpm.columns])
    _tpm_cols = tpm.columns.tolist()
    tpm = tpm[~tpm.index.duplicated(keep="first")]
    # 只要单细胞列 (混样列名含 Pooled / Tumor / Normal)
    sc_cols = [c for c in tpm.columns
               if not any(k in c for k in ("Pooled", "_Tumor", "_Normal", "bulk"))]
    print(f"TPM matrix: {tpm.shape[0]} genes x {tpm.shape[1]} cols; single-cell cols: {len(sc_cols)}")

    # 2) 细胞类型注释
    info = load_sample_info(EXT / "sample_info.txt.gz")
    info["sample"] = info["sample"].astype(str)
    sc_cols = [c for c in sc_cols if c in set(info["sample"])]
    info = info[info["sample"].isin(sc_cols)].drop_duplicates("sample").set_index("sample")
    print("cell types:", info["celltype"].value_counts().to_dict())

    # 3) 每类型: 检出率 frac_expr + mean TPM
    X = tpm[sc_cols].reindex(columns=info.index)
    types = info["celltype"].values
    rows = []
    for t in pd.unique(types):
        cols = np.where(types == t)[0]
        Xt = X.iloc[:, cols]
        pos = (Xt > 0).sum(axis=1)
        rows.append(pd.DataFrame({
            "gene": X.index, "celltype": t, "n_cells": int(len(cols)),
            "frac_expr": (pos / max(len(cols), 1)).values,
            "mean_tpm": Xt.mean(axis=1).values,
        }))
    long = pd.concat(rows, ignore_index=True)
    long.to_csv(PROC / "scrna_crossval_long.csv", index=False)

    # 4) 模块归属复现
    gset = long.gene.unique()
    result = []
    for name, genes in [("Neural_intrinsic", NEURAL_INTRINSIC), ("Perineural_stroma", PERINEURAL_STROMA)]:
        have = [g for g in genes if g in gset]
        sub = long[long.gene.isin(have)]
        for g in have:
            sg = sub[sub.gene == g].set_index("celltype")["frac_expr"]
            top = sg.sort_values(ascending=False)
            result.append({
                "module": name, "gene": g,
                "top_type": top.index[0], "top_frac": float(top.iloc[0]),
                "2nd_type": top.index[1] if len(top) > 1 else "",
                "2nd_frac": float(top.iloc[1]) if len(top) > 1 else np.nan,
            })
    summ = pd.DataFrame(result)
    summ.to_csv(PROC / "scrna_crossval_summary.csv", index=False)

    # 5) 对照打印: 与 GSE176078 的主要结论是否一致
    for name in ("Neural_intrinsic", "Perineural_stroma"):
        sub = summ[summ.module == name]
        n = len(sub)
        print(f"\n== {name} (n={n} genes present in GSE75688)")
        print(sub.top_type.value_counts().head(8).to_string())
    print("\nsaved scrna_crossval_long.csv / scrna_crossval_summary.csv")

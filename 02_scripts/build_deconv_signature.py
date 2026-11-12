# -*- coding: utf-8 -*-
"""从 GSE176078 单细胞参考构建 bulk 反卷积签名矩阵（Stage 10 前置）。

做法等同 CIBERSORTx 的「自建签名矩阵」模式（custom signature matrix）：
按 celltype_minor 分箱 → 每箱算基因的平均 CP10K → 区室 × 基因矩阵。
与 LM22 的差别只在参考来源；混合模型与 NNLS 求解同构。

输出: 01_data/processed/deconv_signature_gse176078.csv
"""
from __future__ import annotations

import sys

import numpy as np
import pandas as pd

sys.path.insert(0, r"D:/嗜神经浸润/02_scripts")
import bulk_lib as BL  # noqa: E402

SRC = BL.ROOT / "01_data" / "external" / "gse176078"
OUT = BL.PROC / "deconv_signature_gse176078.csv"


def main():
    panel = pd.read_csv(BL.RAW / "brca_tcga" / "brca_tcga_panel.csv", index_col=0, nrows=1)
    genes = list(panel.columns)
    print(f"目标基因 {len(genes)}（来自 bulk 面板）")
    S = BL.build_signature(SRC, genes=genes, out_csv=None)
    print(f"签名矩阵 {S.shape}")
    print(S.round(2).to_string())
    # 记录每个区室与每个基因的有效覆盖
    cov = S.notna().sum(axis=1).to_dict()
    print("\n每区室可测基因数:", cov)
    S.to_csv(OUT)
    print(f"-> {OUT}")


if __name__ == "__main__":
    main()

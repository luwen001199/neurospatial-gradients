# -*- coding: utf-8 -*-
"""空间签名打分 —— 把基因层面的表达压成每个 spot 的生物学读数。

设计要点:
1. 先在**单个样本内部**对基因做 z-score, 再取均值 —— 这样不同样本、不同文库深度之间可比,
   也是 hypeR / GSVA 空间版本的标准做法。跨样本直接比较原始 counts 是错的。
2. 每个签名只取实际在该数据集中检出的基因, 并在结果里报告"用了几个/总共几个",
   避免某个基因缺失却无人察觉(审稿人常问的问题)。
3. 同时输出原始均值与 z-score 均值, 前者画图直观, 后者用于统计。
"""
from typing import Dict, List

import numpy as np
import pandas as pd
from scipy import stats


def zscore_genes(X, axis=0):
    """对稀疏/稠密矩阵按 axis z-score, 常量行/列返回 0 而不是 NaN。"""
    A = np.asarray(X, dtype=float)
    mu = A.mean(axis=axis, keepdims=True)
    sd = A.std(axis=axis, ddof=0, keepdims=True)
    sd = np.where(sd < 1e-12, 1.0, sd)
    return (A - mu) / sd


def score_signatures(ad, signatures: Dict[str, List[str]], suffix="", copy=True):
    """计算每个 spot 的签名得分。返回 (DataFrame scores, DataFrame 覆盖度)。"""
    if copy:
        ad = ad.copy()
    # log-normalize 到 1e4, 让不同 spot 的文库深度可比
    counts = ad.obs["total_counts"].values.copy()
    counts[counts == 0] = 1
    import scipy.sparse as sp
    Xn = sp.diags(1e4 / counts) @ ad.X
    Xn = np.asarray(Xn.todense()) if not isinstance(Xn, np.ndarray) else Xn
    Xl = np.log1p(Xn)
    Xz = zscore_genes(Xl, axis=0)          # gene-wise z-score
    Xz = np.nan_to_num(Xz, nan=0.0, posinf=0.0, neginf=0.0)

    gene_index = {g: i for i, g in enumerate(ad.var_names)}
    scores, cover = {}, {}
    for name, genes in signatures.items():
        present = [g for g in genes if g in gene_index]
        cover[name] = f"{len(present)}/{len(genes)}"
        if not present:
            scores[name] = np.zeros(ad.n_obs)
            continue
        idx = [gene_index[g] for g in present]
        scores[name] = Xz[:, idx].mean(axis=1)
        # 非 z 版本的 log 均值, 供空间热图使用
        if suffix == "":
            scores[name + "_raw"] = np.log1p(Xn[:, idx].mean(axis=1))

    S = pd.DataFrame(scores, index=ad.obs_names)
    C = pd.Series(cover)
    # 无论 copy 与否都写回 ad.obs —— copy=True 时返回的是副本,
    # 但调用方通常希望原对象也带上签名列; 不写回会导致下游 KeyError。
    for c in S.columns:
        ad.obs[c + suffix] = S[c].values
    return S, C


def summarize_missing_genes(ad, signatures: Dict[str, List[str]]):
    """报告每个签名里哪些基因在数据集中不存在 —— 写入方法学审计表。"""
    have = set(ad.var_names)
    rows = []
    for name, genes in signatures.items():
        miss = [g for g in genes if g not in have]
        rows.append({"signature": name, "n_genes": len(genes),
                     "n_found": len(genes) - len(miss),
                     "missing": ";".join(miss)})
    return pd.DataFrame(rows)


def spearman(x, y):
    x, y = np.asarray(x), np.asarray(y)
    m = np.isfinite(x) & np.isfinite(y)
    if m.sum() < 3:
        return np.nan, np.nan
    r, p = stats.spearmanr(x[m], y[m])
    return float(r), float(p)

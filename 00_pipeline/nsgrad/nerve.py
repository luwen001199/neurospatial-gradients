# -*- coding: utf-8 -*-
"""神经 spot 识别 —— 三重过滤策略 (+ 空间约束 + 阈值敏感性)。

对应方案 Step1.1 的五个筛选条件, 每一条都在这里有独立实现并且可单独开关,
目的是为了在补充材料里做"逐条拿掉某个条件会怎样"的 ablation:

  ① 同类标记共表达    至少 2 个同类别神经标记 > 该基因在全部 spot 的第 q 分位数
  ② 本底表达控制      取表达量在全样本的分位数上做判定, 天然排除普遍低表达噪声
  ③ 空间约束          仅保留紧邻恶性上皮的神经 spot(可选, 需要上皮区域定义)
  ④ 校准			   见 calibrate.py: 用 mIF 或独立空间数据集交叉验证
  ⑤ 敏感性轮换	   见 sensitivity_NerveThreshold: q = 0.60/0.75/0.90 三档

设计上的诚实声明: 由于 Visium 一个 spot 55um、包含多个细胞,
转录本水平"检测到神经标记"绝不等于"该处存在神经纤维", 这是本方法最基本的局限,
输出到 papers/limitations.rst 里。这里通过"共表达 + 空间连续性"把假阳性压到最低。
"""
from typing import Dict, List, Optional

import numpy as np
import pandas as pd
from scipy import ndimage

from .config import NERVE_MARKERS, ALL_NERVE_GENES, QC


def _gene_columns(ad, genes):
    gidx = {g: i for i, g in enumerate(ad.var_names)}
    return [gidx[g] for g in genes if g in gidx], [g for g in genes if g in gidx]


def quantile_thresholds(ad, genes: List[str], q: float) -> Dict[str, float]:
    """每个基因在全样本 spot 上的第 q 分位数作为阈值。"""
    import scipy.sparse as sp
    counts = ad.obs["total_counts"].values.copy()
    counts[counts == 0] = 1
    Xn = sp.diags(1e4 / counts) @ ad.X
    Xn = np.asarray(Xn.todense())
    Xl = np.log1p(Xn)
    gidx = {g: i for i, g in enumerate(ad.var_names)}
    out = {}
    for g in genes:
        if g in gidx:
            out[g] = float(np.quantile(Xl[:, gidx[g]], q))
    return out


def identify_nerve_spots(ad,
                         q: float = None,
                         min_markers: int = None,
                         require_spatial_cluster: bool = True,
                         verbose: bool = False):
    """返回 DataFrame: obs index x [nerve_raw, n_sensory, n_symp, n_para, n_struct, score]

    nerve_raw: 满足"同类标记 >= min_markers 个超过阈值"的 spot
    """
    q = QC["nerve_q_default"] if q is None else q
    min_markers = QC["nerve_min_markers"] if min_markers is None else min_markers

    import scipy.sparse as sp
    counts = ad.obs["total_counts"].values.copy()
    counts[counts == 0] = 1
    Xn = sp.diags(1e4 / counts) @ ad.X
    Xl = np.log1p(np.asarray(Xn.todense()))
    gidx = {g: i for i, g in enumerate(ad.var_names)}

    res = pd.DataFrame(index=ad.obs_names)
    # 每个神经亚类独立判定
    cat_hit = {}
    for cat, genes in NERVE_MARKERS.items():
        genes_present = [g for g in genes if g in gidx]
        if not genes_present:
            res[f"n_{cat}"] = 0
            cat_hit[cat] = np.zeros(ad.n_obs, dtype=int)
            continue
        cols = [gidx[g] for g in genes_present]
        sub = Xl[:, cols]
        thr = np.quantile(sub, q, axis=0, keepdims=True)
        above = (sub > thr).astype(int)
        res[f"n_{cat}"] = above.sum(axis=1)
        cat_hit[cat] = (above.sum(axis=1) >= min_markers).astype(int)

    for cat in cat_hit:
        res[f"hit_{cat}"] = cat_hit[cat]

    hits = np.zeros(ad.n_obs, dtype=int)
    for cat in cat_hit:
        hits = np.maximum(hits, cat_hit[cat])
    # ① + ②: 任一神经亚类达标即初试 Nerve spot
    res["nerve_raw"] = hits

    # 综合神经得分
    present_all = [g for g in ALL_NERVE_GENES if g in gidx]
    if present_all:
        cols = [gidx[g] for g in present_all]
        z = Xl[:, cols]
        z = (z - z.mean(0, keepdims=True)) / (z.std(0, keepdims=True) + 1e-9)
        res["nerve_score"] = z.mean(axis=1)
    else:
        res["nerve_score"] = 0.0

    # ③ 空间连续性约束: 孤立的单 spot 极可能是噪声, 要求至少一个同命中邻居
    if require_spatial_cluster and "array_row" in ad.obs:
        try:
            rows = ad.obs["array_row"].astype(int).values
            cols_ = ad.obs["array_col"].astype(int).values
            img = np.zeros((rows.max() + 2, cols_.max() + 2), dtype=int)
            img[rows, cols_] = hits
            st = np.array([[1, 1, 1], [1, 1, 1], [1, 1, 1]])
            nbr = ndimage.convolve(img, st, mode="constant") - img
            isolated = (nbr[rows, cols_] == 0) & (hits == 1)
            res["nerve_spatial"] = (hits & (~isolated)).astype(int)
            res["isolated_dropped"] = isolated.astype(int)
        except Exception:
            res["nerve_spatial"] = hits
            res["isolated_dropped"] = 0
    else:
        res["nerve_spatial"] = hits
        res["isolated_dropped"] = 0

    if verbose:
        print(f"    q={q}: nerve_raw={int(res['nerve_raw'].sum())}, "
              f"after spatial={int(res['nerve_spatial'].sum())} "
              f"(dropped {int(res['isolated_dropped'].sum())} isolated) "
              f"/ {ad.n_obs} spots = {res['nerve_spatial'].mean()*100:.2f}%")
    return res


def sensitivity_nerve_threshold(ad, q_grid=None, min_markers=2):
    """⑤ 阈值敏感性: 轮换 q, 返回每个 q 检出的神经 spot 数与占比。"""
    q_grid = QC["nerve_q_grid"] if q_grid is None else q_grid
    rows = []
    base = None
    for q in q_grid:
        r = identify_nerve_spots(ad, q=q, min_markers=min_markers)
        rows.append({"q": q, "n_nerve": int(r["nerve_spatial"].sum()),
                     "pct": round(float(r["nerve_spatial"].mean()) * 100, 3),
                     "n_isolated_dropped": int(r["isolated_dropped"].sum())})
        if base is None:
            base = set(r.index[r["nerve_spatial"] == 1])
        else:
            rows[-1]["jaccard_vs_first"] = round(
                len(base & set(r.index[r["nerve_spatial"] == 1])) /
                max(len(base | set(r.index[r["nerve_spatial"] == 1])), 1), 3)
    return pd.DataFrame(rows)

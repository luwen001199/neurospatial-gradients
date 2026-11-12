# -*- coding: utf-8 -*-
"""打破循环论证: 神经 spot 判定的跨标记集交叉验证。

问题背景
--------
主文用的 Nervesignal 读数里的基因(NGF/UCHL1/CALCA/CALCB/TRPV1/S100B/TUBB3/PLP1)
与判定神经 spot 所用的标记集**高度重叠**, 用它的 Cohen's d / AUC 去"验证"判定
本质是循环论证。真正有信息量的检验是: 用 A 亚类标记判定的 spot,
在**没有参与判定的 B 亚类标记**上是否也真的更高。

做法
----
1. 逐基因 AUC: 每个神经标记单独做"被判定 spot vs 未判定 spot"的判别能力。
   —— 回答"哪些基因在扛着这个判定"。
2. 跨亚类交叉验证: 对只由某一个亚类(sensory / sympathetic / parasympathetic /
   structural)触发的 spot, 计算其余三个亚类标记的平均 z-score 的 AUC。
   —— 若各亚类之间彼此独立地互相抬高, 说明检出的是真实存在的神经结构,
      而不是某一组基因的共表达假象。
3. 阴性对照: 随机打乱 spot 标签 200 次得到 AUC 的经验零分布,
   给出每个交叉验证 AUC 的经验 p 值(而不是只报一个孤零零的 AUC)。
"""
from __future__ import annotations
import sys
import traceback
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse as sp
from scipy import stats

ROOT = Path(r"D:/嗜神经浸润")
sys.path.insert(0, str(ROOT / "00_pipeline"))
from nsgrad import config as CFG      # noqa: E402
from nsgrad import io as nio         # noqa: E402
from nsgrad import nerve as nnerve   # noqa: E402

TAB = CFG.DIR_TAB
TAB.mkdir(parents=True, exist_ok=True)
N_PERM = 200
SEED = 42


def _auc(a, b):
    if len(a) < 3 or len(b) < 3:
        return np.nan
    u = stats.mannwhitneyu(a, b, alternative="two-sided")
    return float(u.statistic / (len(a) * len(b)))


def _perm_p(labels, score, n_perm=N_PERM, rng=None):
    """标签置换下的 AUC 零分布 -> 单侧经验 p (观测 AUC >= 零分布的比例)。"""
    rng = rng or np.random.default_rng(SEED)
    obs = _auc(score[labels == 1], score[labels == 0])
    if not np.isfinite(obs):
        return np.nan, np.nan
    null = []
    for _ in range(n_perm):
        perm = rng.permutation(labels)
        try:
            null.append(_auc(score[perm == 1], score[perm == 0]))
        except Exception:
            null.append(np.nan)
    null = np.array(null, float)
    null = null[np.isfinite(null)]
    if len(null) == 0:
        return obs, np.nan
    p = float((null >= obs).mean())
    return obs, p


def main():
    samples = []
    for d in sorted((CFG.DIR_RAW / "gse210616").glob("GSM*")):
        if not d.is_dir() or list(d.glob("*.part")):
            continue
        if not list(d.glob("*_filtered_feature_bc_matrix.h5")):
            continue
        if not list(d.glob("*tissue_positions*")):
            continue
        samples.append(d)
    print(f"[CV] {len(samples)} ready samples", flush=True)

    rng = np.random.default_rng(SEED)
    gene_rows, cat_rows, gene_rows_cat = [], [], []
    for i, d in enumerate(samples, 1):
        try:
            ad = nio.read_visium_sample(d, sample_id=d.name)
            if ad.n_obs < 50:
                continue
            nv = nnerve.identify_nerve_spots(ad, q=0.90, min_markers=2,
                                             require_spatial_cluster=True)
            called = (nv["nerve_spatial"] == 1).values
            if called.sum() < 5 or (~called).sum() < 20:
                print(f"    [{i}] {d.name} too few called spots ({int(called.sum())})",
                      flush=True)
                continue

            # --- 统一的 gene-wise z 表达矩阵 ---
            counts = ad.obs["total_counts"].values.copy()
            counts[counts == 0] = 1
            Xn = sp.diags(1e4 / counts) @ ad.X
            Xl = np.log1p(np.asarray(Xn.todense()))
            mu, sd = Xl.mean(0), Xl.std(0)
            sd = np.where(sd < 1e-12, 1.0, sd)
            Z = np.nan_to_num((Xl - mu) / sd)
            gidx = {g: j for j, g in enumerate(ad.var_names)}

            # 1) 逐基因 AUC
            for g in CFG.ALL_NERVE_GENES:
                if g not in gidx:
                    continue
                auc = _auc(Z[called, gidx[g]], Z[~called, gidx[g]])
                gene_rows.append(dict(sample=d.name, gene=g,
                                      mean_z=round(float(Z[:, gidx[g]].mean()), 3),
                                      auc=round(auc, 3) if np.isfinite(auc) else np.nan))
            # 2) 跨亚类交叉验证
            cats = list(CFG.NERVE_MARKERS.keys())
            hit = {c: nv[f"hit_{c}"].values.astype(int) for c in cats}
            for c in cats:
                only = hit[c].copy()
                for o in cats:
                    if o != c:
                        only = only * (1 - hit[o])
                only = (only == 1) & (called == 1)
                # 关键: 读数基因必须**完全不属于**触发调用的那个亚类。
                # UCHL1 同时出现在 sensory 与 structural 列表里, 不剔除的话
                # 交叉验证会被这一个基因污染, 变成半循环论证。
                call_genes = set(CFG.NERVE_MARKERS[c])
                others = [g for o in cats if o != c for g in CFG.NERVE_MARKERS[o]
                          if g not in call_genes]
                readout_name = ("markers_excluding_shared_with_call_panel")
                cols = [gidx[g] for g in others if g in gidx]
                if only.sum() < 5 or len(cols) == 0:
                    continue
                readout = Z[:, cols].mean(axis=1)
                lab = np.zeros(ad.n_obs, dtype=int)
                lab[only] = 1
                idx = np.where(lab == 1)[0]
                null_idx = np.where(lab == 0)[0]

                # 逐基因的跨亚类检验: 均值会被大量不表达基因稀释到 0.5,
                # 单基因 AUC 才能说明"到底哪几个基因互相印证"。
                for g in others:
                    if g not in gidx:
                        continue
                    ag = _auc(Z[idx, gidx[g]], Z[null_idx, gidx[g]])
                    if not np.isfinite(ag):
                        continue
                    try:
                        pg = float(stats.mannwhitneyu(Z[idx, gidx[g]],
                                                      Z[null_idx, gidx[g]],
                                                      alternative="two-sided").pvalue)
                    except Exception:
                        pg = np.nan
                    gene_rows_cat.append(dict(sample=d.name, call_category=c,
                                              gene=g, n_only=int(only.sum()),
                                              auc=round(ag, 3), p=round(pg, 5)))
                auc = _auc(readout[idx], readout[null_idx])
                # 标签置换零分布
                pooled = np.concatenate([idx, null_idx])
                null = []
                for _ in range(N_PERM):
                    pp = rng.permutation(pooled)
                    null.append(_auc(readout[pp[:len(idx)]], readout[pp[len(idx):]]))
                null = np.array([x for x in null if np.isfinite(x)], float)
                p = float((null >= auc).mean()) if len(null) else np.nan
                cat_rows.append(dict(sample=d.name, call_category=c,
                                     n_only=int(only.sum()),
                                     readout=readout_name,
                                     auc=round(auc, 3), p_perm=round(min(p, 1.0), 4)))
            print(f"    [{i}/{len(samples)}] {d.name} called={int(called.sum())}", flush=True)
        except Exception as e:
            print(f"    [{i}] {d.name} FAILED {type(e).__name__}: {str(e)[:90]}", flush=True)
            print(traceback.format_exc()[-250:], flush=True)

    G = pd.DataFrame(gene_rows)
    C = pd.DataFrame(cat_rows)
    GC = pd.DataFrame(gene_rows_cat)
    G.to_csv(TAB / "nerve_gene_level_auc.csv", index=False)
    C.to_csv(TAB / "nerve_crossvalidation_category.csv", index=False)
    GC.to_csv(TAB / "nerve_crossvalidation_gene_cross.csv", index=False)

    if len(G):
        print("\n[CV] 逐基因 AUC (跨样本中位数, 按 AUC 降序):")
        gs = (G.dropna(subset=["auc"]).groupby("gene")["auc"]
                .agg(["median", "count"]).sort_values("median", ascending=False))
        print(gs.round(3).to_string())
    if len(C):
        print("\n[CV] 跨亚类交叉验证汇总 (均值读数, 会被不表达基因稀释):")
        cs = C.dropna(subset=["auc"]).groupby("call_category")["auc"].agg(
            ["count", "median", "min", "max"]).round(3)
        print(cs.to_string())
        print(f"\n  p_perm<0.05 的比例: {(C['p_perm'] < 0.05).mean() * 100:.0f}%")
    if len(GC):
        print("\n[CV] 逐基因跨亚类验证 (同一基因不得属于触发调用的亚类):")
        gs = (GC.groupby(["call_category", "gene"])
                .agg(median_auc=("auc", "median"),
                     iqr_lo=("auc", lambda x: np.percentile(x, 25)),
                     iqr_hi=("auc", lambda x: np.percentile(x, 75)),
                     n_samples=("auc", "count"),
                     sig_samples=("p", lambda p: int(((p < 0.05)).sum())))
                .sort_values("median_auc", ascending=False))
        print(gs.round(3).to_string())
        gs.round(3).to_csv(TAB / "nerve_crossvalidation_gene_cross_summary.csv")
    print(f"\n  -> {TAB/'nerve_gene_level_auc.csv'}")
    print(f"  -> {TAB/'nerve_crossvalidation_category.csv'}")
    return G, C


if __name__ == "__main__":
    main()

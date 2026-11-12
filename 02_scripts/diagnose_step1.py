# -*- coding: utf-8 -*-
"""Step1 结果诊断 —— 在画任何 Figure 之前, 先确认参数到底能不能用。

要回答三个具体问题:
  Q1 拟合质量: R2 分布如何? 有多少比例的样本/指标是可信的?
  Q2 参数可辨识性: 有多少 b 撞到了边界 1/2000 或 1/25? 撞界 = 模型退化, 该样本该指标无效
  Q3 神经识别真实性: 被判定为 nerve 的 spot, 其 Nervesignal 是否真的显著高于其他 spot?
     这一步最关键 —— 如果标记集不能把 nerve spot 和非 nerve spot 区分开,
     那后面的梯度、分型、预后全是空中楼阁。用 Cohen's d 和 AUC 量化。
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(r"D:/嗜神经浸润")
sys.path.insert(0, str(ROOT / "00_pipeline"))
from nsgrad import config as CFG
from nsgrad import io as nio, signature as nsig, nerve as nnerve

PROC = CFG.DIR_PROC


def boundary_hits(gp, b_lo=1 / 2000, b_hi=1 / 25, eps=1e-6):
    b = gp["b"].astype(float)
    return int(((np.abs(b - b_lo) < eps) | (np.abs(b - b_hi) < eps)).sum())


def main():
    gp = pd.read_csv(PROC / "gradient_params.csv")
    ss = pd.read_csv(PROC / "sample_summary.csv")
    print("=" * 70)
    print("Q1 拟合质量 (R2)")
    print("=" * 70)
    print(gp.groupby("metric")["r2"].describe()[["count", "mean", "50%", "min", "max"]].round(3).to_string())
    for thr in (0.2, 0.3, 0.5):
        ok = gp[gp["r2"] >= thr]
        print(f"  R2>={thr}: {len(ok)}/{len(gp)} = {len(ok)/len(gp)*100:.0f}%")

    print()
    print("=" * 70)
    print("Q2 参数可辨识性 —— b 是否撞界")
    print("=" * 70)
    rows = []
    for m, sub in gp.groupby("metric"):
        n_lo = int((np.abs(sub["b"].astype(float) - 1 / 2000) < 1e-6).sum())
        n_hi = int((np.abs(sub["b"].astype(float) - 1 / 25) < 1e-6).sum())
        rows.append({"metric": m, "n": len(sub), "at_low_bound": n_lo,
                     "at_high_bound": n_hi,
                     "interpretable": len(sub) - n_lo - n_hi,
                     "pct_interpretable": round((len(sub) - n_lo - n_hi) / len(sub) * 100, 0)})
    print(pd.DataFrame(rows).to_string(index=False))

    print()
    print("=" * 70)
    print("Q3 神经 spot 识别真实性 —— nerve vs non-nerve 的 Nervesignal 差异")
    print("=" * 70)
    samples = []
    for d in sorted((CFG.DIR_RAW / "gse210616").glob("GSM*")):
        if list(d.glob("*.part")) or not list(d.glob("*_filtered_feature_bc_matrix.h5")):
            continue
        if not list(d.glob("*tissue_positions*")):
            continue
        samples.append(d)

    recs = []
    for d in samples:
        try:
            ad = nio.read_visium_sample(d, sample_id=d.name)
            S, _ = nsig.score_signatures(ad, CFG.GRADIENT_SIGNATURES)
            nv = nnerve.identify_nerve_spots(ad)
            mask = (nv["nerve_spatial"] == 1).values
            if mask.sum() < 5 or (~mask).sum() < 5:
                continue
            a = S["Nervesignal"].values[mask]
            b = S["Nervesignal"].values[~mask]
            sd = np.sqrt((a.var() * len(a) + b.var() * len(b)) / (len(a) + len(b) - 2))
            dcoh = (a.mean() - b.mean()) / sd if sd > 0 else np.nan
            try:
                u, p = stats.mannwhitneyu(a, b)
                auc = u / (len(a) * len(b))
            except Exception:
                p, auc = np.nan, np.nan
            recs.append({"sample": d.name, "pct_nerve": round(mask.mean() * 100, 1),
                         "d": round(float(dcoh), 2), "auc": round(float(auc), 3),
                         "p": p, "mean_nerve": round(float(a.mean()), 3),
                         "mean_other": round(float(b.mean()), 3)})
        except Exception as e:
            print(f"   {d.name} err {type(e).__name__}")

    R = pd.DataFrame(recs)
    print(R.to_string(index=False))
    if len(R):
        print(f"  中位数: pct_nerve={R['pct_nerve'].median():.1f}%  "
              f"Cohen d={R['d'].median():.2f}  AUC={R['auc'].median():.3f}")
        print(f"  d>0.8(大效应) 的样本比例: {(R['d']>0.8).mean()*100:.0f}%")
    R.to_csv(PROC / "nerve_validation.csv", index=False)
    return R


if __name__ == "__main__":
    main()

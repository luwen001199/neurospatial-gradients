# -*- coding: utf-8 -*-
"""Step1 主程序: 空间神经源点识别 + 环形梯度建模.

产出
----
01_data/processed/gradient_params.csv       每个样本 x 每个指标的 (a, b, c, r2, 1/b)
01_data/processed/ring_profiles.csv         每个样本每个环带的指标均值
01_data/processed/sample_summary.csv        样本级汇总(神经密度/ fewer的有效酒店 等)
01_data/processed/spot_level/{sample}.h5ad  带全部中间结果的 AnnData

刻意设计成"每样本独立 + 失败不计入":
  空间数据个体差异极大, 任何一个样本出问题都不应该中断整条流水线;
  失败的样本会写进 *_failed.csv 里对照, 而不是被静默丢掉。
"""
import sys
import traceback
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(r"D:/嗜神经浸润")
sys.path.insert(0, str(ROOT / "00_pipeline"))

from nsgrad import config as CFG
from nsgrad import io as nio
from nsgrad import signature as nsig
from nsgrad import nerve as nnerve
from nsgrad import gradient as ngrad

RAW = CFG.DIR_RAW / "gse210616"
PROC = CFG.DIR_PROC
SPOTDIR = PROC / "spot_level"
SPOTDIR.mkdir(parents=True, exist_ok=True)


def discover_samples():
    """直接扫描磁盘找出文件齐全的样本, 不依赖可能过期的 manifest。"""
    out = []
    for d in sorted(RAW.glob("GSM*")):
        if not d.is_dir():
            continue
        h5 = list(d.glob("*_filtered_feature_bc_matrix.h5"))
        pos = list(d.glob("*tissue_positions*"))
        sc = list(d.glob("*scalefactors*"))
        part = list(d.glob("*.part"))
        if h5 and pos and not part:
            out.append((d.name, d))
    return out


def process_one(sample_name, sample_dir, q=None, save_spot=True):
    ad = nio.read_visium_sample(sample_dir, sample_id=sample_name)
    if ad.n_obs < 50:
        raise ValueError(f"too few spots after QC: {ad.n_obs}")

    # --- 签名打分 ---
    S, cover = nsig.score_signatures(ad, CFG.GRADIENT_SIGNATURES)

    # --- 神经源点识别 ---
    if q is None:
        q = CFG.QC["nerve_q_default"]
    nv = nnerve.identify_nerve_spots(ad, q=q)
    mask = (nv["nerve_spatial"] == 1).values
    n_nerve = int(mask.sum())

    summary = {
        "sample": sample_name,
        "n_spots_qc": ad.n_obs,
        "n_genes": ad.n_vars,
        "median_umi": float(np.median(ad.obs["total_counts"])),
        "n_nerve_spots": n_nerve,
        "pct_nerve_spots": round(float(mask.mean()) * 100, 3),
        "n_nerve_isolated_dropped": int(nv["isolated_dropped"].sum()),
        "mean_nerve_score": float(nv["nerve_score"].mean()),
    }

    # --- 环形 + 梯度拟合 ---
    metrics = list(CFG.GRADIENT_SIGNATURES.keys())
    if n_nerve == 0:
        for m in metrics:
            summary[f"b_{m}"] = np.nan
        return summary, None, None, nv

    gdf, dist_field = ngrad.fit_sample_gradients(
        ad, {m: S[m].values for m in metrics}, mask,
        metrics=metrics, verbose=False)

    prof = ngrad.ring_profile(ad, {m: S[m].values for m in metrics},
                              dist_field, nerve_mask=mask)

    for _, r in gdf.iterrows():
        m = r["metric"]
        summary[f"a_{m}"] = r["a"]
        summary[f"b_{m}"] = r["b"]
        summary[f"c_{m}"] = r["c"]
        summary[f"r2_{m}"] = r["r2"]
        summary[f"decaypts_{m}"] = r.get("decay_points_um", np.nan)
        summary[f"fallback_{m}"] = bool(r.get("fallback_used", False))

    # --- 离散降级特征 ---
    disc, n_near, n_far = ngrad.discrete_ring_features(
        ad, {m: S[m].values for m in metrics}, mask, metrics=metrics)
    for k, v in disc.items():
        summary[k] = v
    summary["n_spot_near150"] = n_near
    summary["n_spot_far400"] = n_far

    prof.insert(0, "sample", sample_name)
    gdf.insert(0, "sample", sample_name)

    if save_spot:
        ad.obs = ad.obs.join(nv)
        ad.obs = ad.obs.join(dist_field)
        for c in S.columns:
            ad.obs[c] = S[c].values
        try:
            ad.write_h5ad(SPOTDIR / f"{sample_name}.h5ad")
        except Exception:
            pass
    return summary, gdf, prof, nv


def main():
    samples = discover_samples()
    print(f"[Step1] found {len(samples)} ready samples")
    if not samples:
        return
    summaries, allg, allp, fails = [], [], [], []
    for i, (name, d) in enumerate(samples, 1):
        try:
            s, g, p, _ = process_one(name, d)
            summaries.append(s)
            if g is not None:
                allg.append(g)
            if p is not None:
                allp.append(p)
            b = s.get("b_CollagenIndex", np.nan)
            print(f"  [{i}/{len(samples)}] {name} spots={s['n_spots_qc']:>5} "
                  f"nerve={s['n_nerve_spots']:>4} ({s['pct_nerve_spots']:>5.2f}%) "
                  f"b_coll={b if b == b else 'nan'}", flush=True)
        except Exception as e:
            fails.append({"sample": name, "error": f"{type(e).__name__}: {e}",
                          "tb": traceback.format_exc()[-400:]})
            print(f"  [{i}/{len(samples)}] {name} FAILED {type(e).__name__}: {str(e)[:80]}",
                  flush=True)

    pd.DataFrame(summaries).to_csv(PROC / "sample_summary.csv", index=False)
    if allg:
        pd.concat(allg, ignore_index=True).to_csv(PROC / "gradient_params.csv", index=False)
    if allp:
        pd.concat(allp, ignore_index=True).to_csv(PROC / "ring_profiles.csv", index=False)
    if fails:
        pd.DataFrame(fails).to_csv(PROC / "step1_failed.csv", index=False)

    print(f"\n[Step1] ok={len(summaries)} failed={len(fails)}")
    df = pd.DataFrame(summaries)
    if len(df):
        cols = [c for c in df.columns if c.startswith("b_")]
        print(df[["sample", "n_spots_qc", "n_nerve_spots", "pct_nerve_spots"] +
                 cols[:4]].to_string(index=False))
    print(f"[Step1] saved -> {PROC}")


if __name__ == "__main__":
    main()

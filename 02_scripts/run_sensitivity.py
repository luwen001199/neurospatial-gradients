# -*- coding: utf-8 -*-
"""补充材料的三套稳健性分析 —— 回答审稿人必然会问的三个问题。

S1 神经 spot 识别阈值是不是拍脑袋定的?
    -> q x min_markers x 空间约束 3x2x2 网格轮换, 看占比 / 效应量 / 下游梯度是否翻转

S2  consensus K=2 分型稳不稳?
    -> 特征留一 jackknife / 换聚类算法 / 样本自助法 / 换随机种子重跑 K 选择

S3  METABRIC 的 HR 是不是被少数样本撑起来的?
    -> 1000 次自助重抽的 HR 分布 + 百分位 CI + 未测量混杂 E-value

设计原则: 只做"诚实报告", 任何让结果变差的格子都照实输出, 不筛、不美化。
"""
from __future__ import annotations
import sys
import traceback
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import squareform
from sklearn.cluster import KMeans
from sklearn.mixture import GaussianMixture
from sklearn.metrics import adjusted_rand_score

warnings.filterwarnings("ignore")

ROOT = Path(r"D:/嗜神经浸润")
sys.path.insert(0, str(ROOT / "00_pipeline"))
from nsgrad import config as CFG          # noqa: E402
from nsgrad import io as nio             # noqa: E402
from nsgrad import signature as nsig     # noqa: E402
from nsgrad import nerve as nnerve       # noqa: E402
from nsgrad import gradient as ngrad     # noqa: E402
from nsgrad import surv as nsurv         # noqa: E402

PROC = CFG.DIR_PROC
TAB = CFG.DIR_TAB
TAB.mkdir(parents=True, exist_ok=True)

# 不参与神经 spot 判定的"独立神经读出" —— 用于打破循环论证:
# Nervesignal 里的基因与判定基因集有重叠, 用它算 AUC 有循环论证嫌疑;
# 下面这组 axon/neuron 基因既不出现在 NERVE_MARKERS 也不在 GRADIENT_SIGNATURES 里。
INDEPENDENT_NERVE_SET = ["SNAP25", "SYT1", "NEFH", "PRPH", "GFRA1", "NTRK1",
                         "RET", "TAC1", "ADCYAP1", "GAP43", "STMN2", "SCG2",
                         "VGF", "CHGB", "SLC17A6", "RTN1", "CHRNA7", "KIF5A"]

GRID_Q = [0.85, 0.90, 0.95]
GRID_MM = [2, 3]


# =========================================================================== #
# S1 神经阈值敏感性
# =========================================================================== #
def _effect(nv_scores, other_scores):
    a = np.asarray(nv_scores, float)
    b = np.asarray(other_scores, float)
    if len(a) < 5 or len(b) < 5:
        return np.nan, np.nan, np.nan
    sd = np.sqrt((a.var() * len(a) + b.var() * len(b)) / (len(a) + len(b) - 2))
    dcoh = (a.mean() - b.mean()) / sd if sd > 0 else np.nan
    try:
        u, p = stats.mannwhitneyu(a, b)
        auc = u / (len(a) * len(b))
    except Exception:
        p, auc = np.nan, np.nan
    return float(dcoh), float(auc), float(p)


def _indep_score(S, Indep):
    cols = [c for c in Indep.columns]
    return Indep[cols].mean(axis=1).values if cols else np.zeros(len(S))


def run_threshold_sensitivity():
    samples = []
    for d in sorted((CFG.DIR_RAW / "gse210616").glob("GSM*")):
        if not d.is_dir() or list(d.glob("*.part")):
            continue
        if not list(d.glob("*_filtered_feature_bc_matrix.h5")):
            continue
        if not list(d.glob("*tissue_positions*")):
            continue
        samples.append(d)
    print(f"[S1] {len(samples)} ready samples", flush=True)

    rows = []
    for i, d in enumerate(samples, 1):
        try:
            ad = nio.read_visium_sample(d, sample_id=d.name)
            if ad.n_obs < 50:
                continue
            S, _ = nsig.score_signatures(ad, CFG.GRADIENT_SIGNATURES)
            Indep, cov = nsig.score_signatures(ad, {"IndependentNerve": INDEPENDENT_NERVE_SET})
            print(f"    [{i}/{len(samples)}] {d.name} "
                  f"independent nerve set coverage {cov.get('IndependentNerve')}",
                  flush=True)
            for q in GRID_Q:
                for mm in GRID_MM:
                    for sp in (True, False):
                        nv = nnerve.identify_nerve_spots(ad, q=q, min_markers=mm,
                                                         require_spatial_cluster=sp)
                        mask = (nv["nerve_spatial"] == 1).values
                        rec = {"sample": d.name, "q": q, "min_markers": mm,
                               "spatial_constraint": sp,
                               "n_spots": int(ad.n_obs),
                               "n_nerve": int(mask.sum()),
                               "pct_nerve": round(float(mask.mean()) * 100, 3),
                               "n_isolated_dropped": int(nv["isolated_dropped"].sum())}
                        if 5 <= mask.sum() <= ad.n_obs - 5:
                            dm, am, pm = _effect(S["Nervesignal"].values[mask],
                                                 S["Nervesignal"].values[~mask])
                            di, ai, pi = _effect(Indep["IndependentNerve"].values[mask],
                                                 Indep["IndependentNerve"].values[~mask])
                            rec.update(cohen_d_main=round(dm, 3), auc_main=round(am, 3),
                                       p_main=float(f"{pm:.3g}"),
                                       cohen_d_indep=round(di, 3), auc_indep=round(ai, 3),
                                       p_indep=float(f"{pi:.3g}"))
                            dist_field, _ = ngrad.build_distance_fields(ad, mask)
                            dist = dist_field["distance_to_nearest_nerve"].values
                            for metric in ["Nervesignal", "CollagenIndex", "CD8density"]:
                                r = ngrad.fit_decay_aggregated(dist, S[metric].values)
                                rec[f"b_{metric}"] = r["b"]
                                rec[f"decay_um_{metric}"] = r.get("decay_points_um", np.nan)
                                rec[f"r2_{metric}"] = r["r2"]
                        else:
                            rec.update(cohen_d_main=np.nan, auc_main=np.nan, p_main=np.nan,
                                       cohen_d_indep=np.nan, auc_indep=np.nan, p_indep=np.nan)
                        rows.append(rec)
        except Exception as e:
            print(f"    [{i}/{len(samples)}] {d.name} FAILED {type(e).__name__}: {str(e)[:100]}",
                  flush=True)
            print(traceback.format_exc()[-300:], flush=True)
    out = pd.DataFrame(rows)
    out.to_csv(TAB / "sensitivity_nerve_threshold.csv", index=False)
    print("\n[S1] 阈值敏感性汇总 (按配置聚合中位数):")
    if len(out):
        agg = out.groupby(["q", "min_markers", "spatial_constraint"]).agg(
            n_samples=("sample", "nunique"),
            pct_nerve=("pct_nerve", "median"),
            cohen_d_main=("cohen_d_main", "median"),
            auc_main=("auc_main", "median"),
            cohen_d_indep=("cohen_d_indep", "median"),
            auc_indep=("auc_indep", "median"),
            decay_um_Nervesignal=("decay_um_Nervesignal", "median"),
            r2_Nervesignal=("r2_Nervesignal", "median"),
        ).round(3).reset_index()
        print(agg.to_string(index=False))
        agg.to_csv(TAB / "sensitivity_nerve_threshold_summary.csv", index=False)
    print(f"  -> {TAB/'sensitivity_nerve_threshold.csv'}")
    return out


# =========================================================================== #
# S2 分型稳健性
# =========================================================================== #
def _features(ss: pd.DataFrame):
    sys.path.insert(0, str(ROOT / "02_scripts"))
    from run_step2_subtyping import build_feature_matrix
    F = build_feature_matrix(ss)
    F = F.dropna(axis=1, how="all").fillna(F.median(numeric_only=True))
    return F


def run_clustering_robustness():
    ss = pd.read_csv(PROC / "sample_summary.csv")
    sub = pd.read_csv(PROC / "sample_subtypes.csv")
    F = _features(ss)
    common = [s for s in sub["sample"] if s in F.index]
    F = F.loc[common]
    ref = sub.set_index("sample").loc[common, "cluster"].values
    X = np.asarray(F.values, dtype=float)
    mu = X.mean(0)
    sd = X.std(0)
    sd = np.where(sd < 1e-12, 1.0, sd)
    Xz = (X - mu) / sd
    print(f"[S2] {len(common)} samples x {Xz.shape[1]} features; "
          f"reference clusters {np.bincount(ref)[1:]}")

    rows = []

    # (i) 特征留一 jackknife
    for j, col in enumerate(F.columns):
        Xj = np.delete(Xz, j, axis=1)
        try:
            lab = fcluster(linkage(Xj, method="ward"), t=2, criterion="maxclust")
            rows.append(dict(test="feature_jackknife", setting=col,
                             ari=round(float(adjusted_rand_score(ref, lab)), 3)))
        except Exception:
            rows.append(dict(test="feature_jackknife", setting=col, ari=np.nan))

    # (ii) 换聚类算法 / 距离
    algs = {
        "ward (reference)": lambda M: fcluster(linkage(M, method="ward"), t=2,
                                               criterion="maxclust"),
        "average linkage": lambda M: fcluster(linkage(M, method="average"), t=2,
                                              criterion="maxclust"),
        "complete linkage": lambda M: fcluster(linkage(M, method="complete"), t=2,
                                               criterion="maxclust"),
        "kmeans": lambda M: KMeans(2, n_init=50, random_state=42).fit_predict(M),
        "GMM(full)": lambda M: GaussianMixture(2, covariance_type="full",
                                               random_state=42, n_init=10)
                                   .fit_predict(M),
    }
    for name, fn in algs.items():
        try:
            lab = fn(Xz)
            rows.append(dict(test="algorithm", setting=name,
                             ari=round(float(adjusted_rand_score(ref, lab)), 3)))
        except Exception as e:
            rows.append(dict(test="algorithm", setting=name, ari=np.nan))

    # (iii) 样本自助法: 80% 子抽样重复聚类, 与全数据划分比 ARI
    rng = np.random.default_rng(42)
    n = len(Xz)
    aris, sizes = [], []
    for _ in range(1000):
        idx = np.sort(rng.choice(n, int(n * 0.8), replace=False))
        subX = Xz[idx]
        if len(np.unique(subX, axis=0)) < 2:
            continue
        try:
            lab = fcluster(linkage(subX, method="ward"), t=2, criterion="maxclust")
            cnt = np.bincount(lab)[1:]
            if min(cnt) < 3:          # 退化划分不算
                continue
            aris.append(adjusted_rand_score(ref[idx], lab))
            sizes.append(len(idx))
        except Exception:
            continue
    if aris:
        np.save(TAB / "sensitivity_clustering_bootstrap.npy", np.array(aris, float))
        rows.append(dict(test="sample_bootstrap", setting="80% x1000",
                         ari=round(float(np.median(aris)), 3)))
        print(f"[S2] bootstrap ARI median {np.median(aris):.3f} "
              f"(IQR {np.percentile(aris,25):.3f}-{np.percentile(aris,75):.3f}, "
              f"n={len(aris)})")

    # (iv) 换随机种子重跑 K 选择: K=2..4 的共识 silhouette 是否稳定地偏爱 K=2
    sys.path.insert(0, str(ROOT / "02_scripts"))
    from run_step2_subtyping import pick_k
    for seed in (1, 7, 42, 2024):
        try:
            res = pick_k(Xz, k_range=(2, 3, 4), seed=seed, n_boot=300)
            for k, v in res.items():
                rows.append(dict(test="k_selection_seed", setting=f"seed={seed},K={k}",
                                 ari=round(float(v["silhouette_consensus"]), 3)))
            best = max(res.items(), key=lambda kv: kv[1]["silhouette_consensus"])
            print(f"[S2] seed={seed}: best K={best[0]} "
                  f"sil={best[1]['silhouette_consensus']:.3f}")
        except Exception as e:
            print(f"[S2] seed={seed} failed {type(e).__name__}")

    out = pd.DataFrame(rows)
    out.to_csv(TAB / "sensitivity_clustering.csv", index=False)
    print(out.to_string(index=False))
    print(f"  -> {TAB/'sensitivity_clustering.csv'}")
    return out, aris


# =========================================================================== #
# S3 METABRIC 自助法 + E-value
# =========================================================================== #
def _evalue(rr):
    """VanderWeele-Ding E-value: 解释掉该关联所需的未测量混杂 RR 下限。"""
    if not np.isfinite(rr) or rr <= 0:
        return np.nan
    r = rr if rr >= 1 else 1.0 / rr
    return r + np.sqrt(r * (r - 1))


def run_metabric_robustness(n_boot=500, seed=42):
    d = pd.read_csv(PROC / "metabric_analysis_table.csv")
    print(f"[S3] METABRIC table {d.shape}")
    rng = np.random.default_rng(seed)
    rows = []
    boot_store = {}

    for outcome, tcol, ecol in [("RFS", "RFS_MONTHS", "RFS_EVENT"),
                                ("OS", "OS_MONTHS", "OS_EVENT")]:
        sub = d.dropna(subset=[tcol, ecol])
        sub = sub[sub[tcol] > 0]
        for var in ["NPPscore", "NPP_resid_stroma"]:
            covs = [var, "AGE", "nodes_pos", "ER_pos", "chemo", "hormone", "CAF_ECM"]
            s2 = sub.dropna(subset=covs)
            Z = pd.DataFrame({c: (s2[c] - s2[c].mean()) / s2[c].std(ddof=0)
                              for c in covs})
            keep = [c for c in covs if Z[c].std(ddof=0) > 1e-9]
            t = s2[tcol].values.astype(float)
            e = s2[ecol].values.astype(float)
            M = Z[keep].values
            point = [r for r in nsurv.cox(t, e, M, names=keep) if r.name == var][0]
            hrs = []
            for b in range(n_boot):
                idx = rng.integers(0, len(s2), len(s2))
                r = nsurv.cox(t[idx], e[idx], M[idx], names=keep)
                r = [x for x in r if x.name == var][0]
                if np.isfinite(r.hr):
                    hrs.append(r.hr)
                if len(hrs) and len(hrs) % 250 == 0:
                    pass
                if b >= n_boot:
                    break
            hrs = np.array(hrs, float)
            boot_store[f"{outcome}_{var}"] = hrs
            lo, hi = np.nanpercentile(hrs, [2.5, 97.5])
            rows.append(dict(
                outcome=outcome, variable=var, model="multivariable(z per SD)",
                n=int(len(s2)), events=int(e.sum()),
                hr=round(point.hr, 3), ci_lo=round(point.lo, 3), ci_hi=round(point.hi, 3),
                p=float(f"{point.p:.3g}"),
                boot_hr_median=round(float(np.nanmedian(hrs)), 3),
                boot_ci_lo=round(float(lo), 3), boot_ci_hi=round(float(hi), 3),
                boot_p_lt1=round(float(np.mean(hrs < 1)), 3),
                boot_n=int(len(hrs)),
                evalue_point=round(_evalue(point.hr), 2),
                evalue_ci=round(_evalue(point.hi if abs(point.hi - 1) < abs(point.lo - 1)
                                        else point.lo), 2)))
    out = pd.DataFrame(rows)
    out.to_csv(TAB / "sensitivity_metabric_bootstrap.csv", index=False)
    print(out.to_string(index=False))
    np.savez(TAB / "sensitivity_metabric_bootstrap_draws.npz", **boot_store)
    print(f"  -> {TAB/'sensitivity_metabric_bootstrap.csv'}")
    return out, boot_store


if __name__ == "__main__":
    which = sys.argv[1] if len(sys.argv) > 1 else "all"
    if which in ("all", "s1"):
        run_threshold_sensitivity()
    if which in ("all", "s2"):
        run_clustering_robustness()
    if which in ("all", "s3"):
        run_metabric_robustness()

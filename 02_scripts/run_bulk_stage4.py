# -*- coding: utf-8 -*-
"""Stage 4 + 6 —— 空间先验模块对 NAC 病理完全缓解 (pCR) 的预测力,
并用 510 组 size-matched 随机基因集做随机基线校准 (Stage 6)。

主队列: GSE163882 (222 NAC 前活检, 其中 TNBC 90 例, pCR 38)
复制  : GSE25065  (198 HER2-negative NAC, 有 pCR/RCB)

判定口径:
    OR per SD + 95% CI + Wald P + 交叉验证 AUC
    模块的 AUC / OR 在 510 组随机基因集中的百分位 (= 经验 P)
"""
from __future__ import annotations

import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy import stats
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import StratifiedKFold

sys.path.insert(0, r"D:/嗜神经浸润/02_scripts")
import bulk_lib as BL  # noqa: E402

warnings.filterwarnings("ignore")
OUT = BL.TAB
RAW = BL.RAW
SEED = 20251108
N_RANDOM = 510


def logit_or(x: np.ndarray, y: np.ndarray, cov: np.ndarray | None = None):
    """单变量/多变量 logistic -> (OR, lo, hi, p, n_event)。x 已标准化。"""
    X = x.reshape(-1, 1)
    names = ["x"]
    if cov is not None and cov.size:
        X = np.column_stack([X, cov])
        names += [f"c{i}" for i in range(cov.shape[1])]
    X = sm.add_constant(X, has_constant="add")
    try:
        fit = sm.Logit(y, X).fit(disp=0, maxiter=200)
    except Exception:
        return np.nan, np.nan, np.nan, np.nan, int(y.sum())
    b = fit.params[1]
    se = fit.bse[1]
    return float(np.exp(b)), float(np.exp(b - 1.96 * se)), float(np.exp(b + 1.96 * se)), \
        float(fit.pvalues[1]), int(y.sum())


def cv_auc(x: np.ndarray, y: np.ndarray, cov: np.ndarray | None = None, folds=5, reps=5):
    """重复分层 K 折 AUC（均值 ± SD）。"""
    X = x.reshape(-1, 1)
    if cov is not None and cov.size:
        X = np.column_stack([X, cov])
    aucs = []
    rng = np.random.default_rng(SEED)
    for _ in range(reps):
        skf = StratifiedKFold(n_splits=folds, shuffle=True, random_state=int(rng.integers(1e6)))
        pred = np.zeros(len(y))
        for tr, te in skf.split(X, y):
            m = LogisticRegression(max_iter=1000, C=1.0)
            m.fit(X[tr], y[tr])
            pred[te] = m.predict_proba(X[te])[:, 1]
        aucs.append(roc_auc_score(y, pred))
    return float(np.mean(aucs)), float(np.std(aucs))


def random_baseline(expr: pd.DataFrame, y: np.ndarray, size: int, n=N_RANDOM,
                    seed=SEED) -> np.ndarray:
    """n 组 size-matched 随机基因集的 AUC 分布。"""
    rng = np.random.default_rng(seed)
    genes = np.array(expr.columns)
    aucs = np.empty(n)
    for i in range(n):
        pick = rng.choice(genes, size=min(size, len(genes)), replace=False)
        s = BL.score_mean(expr, list(pick)).to_numpy()
        if np.isnan(s).all():
            aucs[i] = np.nan
            continue
        aucs[i] = roc_auc_score(y, s)
    return aucs


def run_cohort(name: str, expr: pd.DataFrame, y: np.ndarray, age: np.ndarray | None,
               mods: dict, do_random=True):
    rows = []
    sc = BL.score_all(expr, mods)
    print(f"\n[{name}] n={len(y)}, pCR={int(y.sum())} ({y.mean():.1%}), 模块 {sc.shape[1]}")
    print(f"  {'模块':22s} {'OR/SD':>7s} {'95%CI':>16s} {'P':>9s} {'AUC':>6s}")
    for m in ["NPP58", "Perineural_stroma", "Neural_intrinsic", "CAF_ECM",
              "CD8_cytotox", "LRRC15_myCAF", "TLS", "TGFb_core"]:
        if m not in sc:
            continue
        s = sc[m]
        ok = s.notna().to_numpy()
        if ok.sum() < 30:
            continue
        z = ((s - s.mean()) / s.std(ddof=0)).to_numpy()[ok]
        yy = y[ok]
        aa = None
        if age is not None:
            ag = age[ok]
            aa = np.where(np.isnan(ag), np.nanmean(ag), ag).reshape(-1, 1)
        orv, lo, hi, p, nev = logit_or(z, yy.astype(float, copy=False))
        orv2, lo2, hi2, p2, _ = logit_or(z, yy.astype(float, copy=False), aa) if aa is not None else (np.nan,) * 5
        auc, auc_sd = cv_auc(z, yy, aa)
        rows.append({"cohort": name, "module": m, "n": int(ok.sum()), "n_event": nev,
                     "OR_per_SD": round(orv, 3), "ci_lo": round(lo, 3), "ci_hi": round(hi, 3),
                     "p": p,
                     "OR_adj_age": round(orv2, 3) if not np.isnan(orv2) else np.nan,
                     "ci_lo_adj": round(lo2, 3) if not np.isnan(lo2) else np.nan,
                     "ci_hi_adj": round(hi2, 3) if not np.isnan(hi2) else np.nan,
                     "p_adj": p2,
                     "AUC_cv": round(auc, 3), "AUC_sd": round(auc_sd, 3)})
        print(f"  {m:22s} {orv:7.3f} [{lo:6.3f},{hi:6.3f}] {p:9.2e} {auc:6.3f}")

    # ---- Stage 6 随机基线
    if do_random:
        rnd_rows = []
        for m, genes in [("Perineural_stroma", mods["Perineural_stroma"]),
                         ("Neural_intrinsic", mods["Neural_intrinsic"]),
                         ("CAF_ECM", mods["CAF_ECM"]),
                         ("NPP58", mods["NPP_up"] + mods["NPP_down"])]:
            s = sc[m] if m in sc else BL.score_diff(expr, mods["NPP_up"], mods["NPP_down"])
            ok = s.notna().to_numpy()
            if ok.sum() < 30:
                continue
            yv = y[ok]
            obs = roc_auc_score(yv, s.to_numpy()[ok])
            dist = random_baseline(expr.loc[ok], yv, size=len(set(genes)))
            emp = float((np.abs(dist - 0.5) >= abs(obs - 0.5)).mean())
            rnd_rows.append({"cohort": name, "module": m, "n_genes": len(set(genes)),
                             "obs_auc": round(obs, 3),
                             "random_median": round(float(np.nanmedian(dist)), 3),
                             "random_p95": round(float(np.nanpercentile(dist, 95)), 3),
                             "empirical_p": emp, "n_random": int(np.isfinite(dist).sum())})
            print(f"  [Stage6] {m:22s} obs AUC={obs:.3f} vs random median "
                  f"{np.nanmedian(dist):.3f} (p95 {np.nanpercentile(dist,95):.3f}), emp P={emp:.3f}")
        return rows, rnd_rows
    return rows, []


def main():
    mods = BL.load_modules()
    all_rows, all_rnd = [], []

    # ---------------- GSE163882
    from run_bulk_stage1 import load_gse163882
    expr, clin = load_gse163882()
    print(f"[GSE163882] expr {expr.shape}, clin {clin.shape}, clin cols={list(clin.columns)[:12]}")
    if "response to nac" in clin:
        resp = clin["response to nac"].astype(str).str.strip()
        y = (resp == "pCR").to_numpy()
        age = pd.to_numeric(clin.get("age"), errors="coerce").to_numpy()
        er = clin.get("estrogen receptor status", pd.Series(index=clin.index, dtype=object)).astype(str)
        pr = clin.get("progesterone receptor status", pd.Series(index=clin.index, dtype=object)).astype(str)
        h2 = clin.get("her2 receptor status", pd.Series(index=clin.index, dtype=object)).astype(str)
        tn = (er.str.upper() == "N") & (pr.str.upper() == "N") & (h2.str.upper() == "N")
        ok_all = np.ones(len(y), bool)
        r, rr = run_cohort("GSE163882_all", expr.loc[ok_all], y[ok_all], age[ok_all], mods)
        all_rows += r; all_rnd += rr
        tnv = tn.to_numpy()
        r, rr = run_cohort("GSE163882_TNBC", expr.loc[tnv], y[tnv], age[tnv], mods)
        all_rows += r; all_rnd += rr

    # ---------------- GSE25065
    try:
        from run_bulk_stage1 import load_gse25065
        expr2, clin2 = load_gse25065()
        print(f"\n[GSE25065] expr {expr2.shape}, clin cols={list(clin2.columns)}")
        col = [c for c in clin2.columns if "pcr" in c.lower()]
        if col:
            resp = clin2[col[0]].astype(str).str.strip().str.upper()
            keep = resp.isin(["PCR", "RD"]).to_numpy()
            y2 = (resp == "PCR").to_numpy()
            age2 = None
            r, rr = run_cohort("GSE25065_all", expr2.loc[keep], y2[keep], age2, mods)
            all_rows += r; all_rnd += rr
    except Exception as ex:
        print(f"[skip] GSE25065: {ex}")

    pd.DataFrame(all_rows).to_csv(OUT / "stage4_pcr_association.csv", index=False)
    if all_rnd:
        pd.DataFrame(all_rnd).to_csv(OUT / "stage6_random_baseline.csv", index=False)
    print(f"\n输出 -> {OUT}")


if __name__ == "__main__":
    main()

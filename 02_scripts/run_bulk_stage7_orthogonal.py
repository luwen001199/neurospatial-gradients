# -*- coding: utf-8 -*-
"""Stage 7 —— PNI 模块与纤维化的正交化：扣除 desmoplasia 后还剩什么。

背景
----
Stage 1 已证实 Perineural_stroma 与通用 CAF/ECM 在 5 个 bulk 队列中 ρ=0.96–0.98，
即"基质总量"层面的信息被 desmoplasia 完全吸收。若要回答"PNI 相关基质是否有
独立于一般纤维化的生物学"，唯一严谨做法是把模块对 CAF/ECM 做秩空间正交化，
用残差（= 给定纤维化水平后剩余的那部分）重新跑免疫 / 耐药 / 生存三条链。

三级衰减口径（同一 target，三种校正）
    raw     : 模块原值 ~ target
    +CAF    : 控制 CAF/ECM（= 纤维化总量）
    +CAF+P  : 再控制肿瘤纯度代理（EPCAM/KRT 上皮评分）

输出
----
stage7_orthogonality.csv        模块可被纤维化解 %（R²）、残差 SD
stage7_attenuation.csv          三级衰减并排表（核心结果）
stage7_residual_immune.csv      残差 vs 免疫读数（含 FDR）
stage7_residual_pcr.csv         残差 vs pCR
stage7_residual_survival.csv    残差 vs OS（TCGA-BRCA）
"""
from __future__ import annotations

import sys
import warnings

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
SEED = 20251108

FOCUS = ["Perineural_stroma", "Neural_intrinsic", "NPP58", "LRRC15_myCAF",
         "CAF_ECM", "SPARC_POSTN_CORE", "COLLAGEN_MATURATION"]
IMMUNE = ["CD8_cytotox", "Treg", "M2_macro", "SPP1_macro", "TGFb_core",
          "CXCL12_axis", "LGALS1_axis", "IFN_gamma", "TLS", "MHC_II"]


# ------------------------------------------------------------------ 队列
def zs(s: pd.Series) -> pd.Series:
    return (s - s.mean()) / (s.std(ddof=0) or np.nan)


def build_cohorts() -> dict[str, tuple[pd.DataFrame, dict]]:
    from run_bulk_stage1 import load_tcga, load_gse163882, load_gse25065
    out = {}

    e, c = load_tcga()
    neg = lambda s: s.astype(str).str.lower().str.startswith(("neg", "n"))
    tn = neg(c["ER_Status_nature2012"]) & neg(c["PR_Status_nature2012"]) & \
        neg(c["HER2_Final_Status_nature2012"])
    meta = {
        "os_time": pd.to_numeric(c["OS_TIME"], errors="coerce"),
        "os_event": pd.to_numeric(c["OS_EVENT"], errors="coerce"),
        "age": pd.to_numeric(c["AGE"], errors="coerce"),
        "node": c["NODE"].astype(str),
        "stage": c["STAGE"].astype(str),
        "pam50": c["PAM50"].astype(str),
    }
    out["TCGA-BRCA"] = (e, meta)
    out["TCGA-TNBC"] = (e.loc[tn.fillna(False).values], {k: v.reindex(e.index)[tn.fillna(False).values]
                                                         for k, v in meta.items()})

    e, c2 = load_gse163882()
    resp = c2["response to nac"].astype(str).str.strip() if "response to nac" in c2 else pd.Series("", index=e.index)
    er = c2.get("estrogen receptor status", pd.Series("", index=e.index)).astype(str)
    pr = c2.get("progesterone receptor status", pd.Series("", index=e.index)).astype(str)
    h2 = c2.get("her2 receptor status", pd.Series("", index=e.index)).astype(str)
    tn2 = (er.str.upper() == "N") & (pr.str.upper() == "N") & (h2.str.upper() == "N")
    m2 = {"pcr": (resp == "pCR").astype(float).where(resp.isin(["pCR", "RD"])),
          "age": pd.to_numeric(c2.get("age"), errors="coerce"),
          "os_time": pd.Series(np.nan, index=e.index), "os_event": pd.Series(np.nan, index=e.index),
          "node": c2.get("node status", pd.Series("", index=e.index)).astype(str),
          "stage": pd.Series("", index=e.index), "pam50": pd.Series("", index=e.index)}
    out["GSE163882"] = (e, m2)
    out["GSE163882-TNBC"] = (e.loc[tn2.fillna(False).values],
                            {k: v[tn2.fillna(False).values] for k, v in m2.items()})

    e, c3 = load_gse25065()
    col = [x for x in c3.columns if "pcr" in x.lower()]
    pcr3 = pd.Series(np.nan, index=e.index)
    if col:
        r = c3[col[0]].astype(str).str.strip().str.upper()
        pcr3 = (r == "PCR").astype(float).where(r.isin(["PCR", "RD"]))
    m3 = {"pcr": pcr3, "age": pd.Series(np.nan, index=e.index),
          "os_time": pd.Series(np.nan, index=e.index), "os_event": pd.Series(np.nan, index=e.index),
          "node": pd.Series("", index=e.index), "stage": pd.Series("", index=e.index),
          "pam50": pd.Series("", index=e.index)}
    out["GSE25065"] = (e, m3)
    return out


def score_table(expr: pd.DataFrame, mods: dict) -> pd.DataFrame:
    sc = BL.score_all(expr, mods)
    for k, g in BL.PROGRAMS.items():
        s = BL.score_mean(expr, g)
        if s.notna().any():
            sc[k] = s
    sc["PurityProxy"] = BL.score_mean(expr, BL.PURITY_PROXY)
    return sc


# ------------------------------------------------------------------ 统计
def logit_or(x, y, cov: np.ndarray | None = None):
    X = np.asarray(x, float).reshape(-1, 1)
    if cov is not None and cov.size:
        X = np.column_stack([X, cov])
    X = sm.add_constant(X, has_constant="add")
    try:
        fit = sm.Logit(y, X).fit(disp=0, maxiter=300)
    except Exception:
        return np.nan, np.nan, np.nan, np.nan
    b, se = fit.params[1], fit.bse[1]
    return float(np.exp(b)), float(np.exp(b - 1.96 * se)), float(np.exp(b + 1.96 * se)), float(fit.pvalues[1])


def cv_auc(x, y, cov=None, folds=5, reps=5):
    X = np.asarray(x, float).reshape(-1, 1)
    if cov is not None and cov.size:
        X = np.column_stack([X, cov])
    rng = np.random.default_rng(SEED)
    aucs = []
    for _ in range(reps):
        skf = StratifiedKFold(folds, shuffle=True, random_state=int(rng.integers(1e6)))
        pred = np.zeros(len(y))
        for tr, te in skf.split(X, y):
            m = LogisticRegression(max_iter=1000).fit(X[tr], y[tr])
            pred[te] = m.predict_proba(X[te])[:, 1]
        aucs.append(roc_auc_score(y, pred))
    return float(np.mean(aucs))


def cox_hr(t: np.ndarray, d: np.ndarray, x: np.ndarray, cov: np.ndarray | None = None):
    """委托 bulk_lib 自实现（statsmodels PHReg 在 HR≈1 时 bse 退化）。"""
    return BL.cox_hr(t, d, x, cov)


# ------------------------------------------------------------------ 主流程
def main():
    mods = BL.load_modules()
    cohorts = build_cohorts()

    orth, atten, imm, pcr_rows, surv_rows = [], [], [], [], []

    for name, (expr, meta) in cohorts.items():
        if len(expr) < 30:
            print(f"[skip] {name} n={len(expr)}")
            continue
        sc = score_table(expr, mods)
        if "CAF_ECM" not in sc or "PurityProxy" not in sc:
            print(f"[skip] {name} 关键模块缺失")
            continue
        print(f"\n[{name}] n={len(sc)}")

        # ---- 1) 正交化：模块被纤维化解释的比例
        resid = {}
        for m in FOCUS:
            if m not in sc:
                continue
            if m == "CAF_ECM":
                continue
            r, r2 = BL.rank_resid(sc[m], sc["CAF_ECM"])
            resid[m] = r
            orth.append({"cohort": name, "module": m, "n": int(r.notna().sum()),
                         "R2_vs_CAFECM": round(r2, 4),
                         "pct_explained": round(100 * r2, 1),
                         "resid_sd": round(float(np.nanstd(r)), 3)})
        # 联合正交化（CAF/ECM + 纯度）：更保守口径
        resid2 = {}
        for m in FOCUS:
            if m not in sc or m == "CAF_ECM":
                continue
            r, r2 = BL.rank_resid(sc[m], sc[["CAF_ECM", "PurityProxy"]])
            resid2[m] = r

        # ---- 2) 三级衰减
        for m in ["Perineural_stroma", "Neural_intrinsic", "NPP58", "LRRC15_myCAF"]:
            if m not in sc:
                continue
            for t in IMMUNE:
                if t not in sc or t == m:
                    continue
                d = pd.DataFrame({"x": sc[m], "y": sc[t], "caf": sc["CAF_ECM"],
                                  "pur": sc["PurityProxy"]}).dropna()
                if len(d) < 40:
                    continue
                r0, p0 = stats.spearmanr(d["x"], d["y"])
                r1, p1, _ = BL.partial_spearman(d["x"], d["y"], d["caf"])
                r2v, p2, _ = BL.partial_spearman(d["x"], d["y"], d[["caf", "pur"]])
                atten.append({"cohort": name, "module": m, "target": t, "n": len(d),
                              "rho_raw": round(r0, 3), "p_raw": p0,
                              "rho_caf": round(r1, 3), "p_caf": p1,
                              "rho_caf_purity": round(r2v, 3), "p_caf_purity": p2,
                              "shrink_pct": round(100 * (1 - abs(r1) / abs(r0)), 1) if abs(r0) > 1e-6 else np.nan})
                imm.append({"cohort": name, "module_resid": m, "target": t, "n": len(d),
                            "rho_resid": round(r1, 3), "p_resid": p1})

        # ---- 3) pCR
        y = meta.get("pcr")
        if y is not None and y.notna().sum() >= 60 and y.sum() >= 15:
            ok = y.notna().to_numpy()
            yy = y.to_numpy()[ok].astype(float)
            age = meta.get("age")
            av = pd.to_numeric(age, errors="coerce").to_numpy()[ok] if age is not None else None
            cov = None
            if av is not None and np.isfinite(av).sum() > 0.8 * len(av):
                cov = np.where(np.isnan(av), np.nanmean(av), av).reshape(-1, 1)
            for m in FOCUS:
                if m not in sc:
                    continue
                for tag, s in (("raw", sc[m]), ("resid", resid.get(m))):
                    if s is None:
                        continue
                    sv = s.to_numpy()[ok]
                    if not np.isfinite(sv).all():
                        continue
                    z = (sv - sv.mean()) / sv.std(ddof=0)
                    orv, lo, hi, p = logit_or(z, yy, cov)
                    auc = cv_auc(z, yy, cov)
                    pcr_rows.append({"cohort": name, "module": m, "scale": tag,
                                     "n": int(ok.sum()), "n_event": int(yy.sum()),
                                     "OR_per_SD": round(orv, 3), "ci_lo": round(lo, 3),
                                     "ci_hi": round(hi, 3), "p": p, "AUC_cv": round(auc, 3)})
            print(f"  pCR n={int(ok.sum())}, event={int(yy.sum())}")

        # ---- 4) 生存
        t_ = pd.to_numeric(meta.get("os_time"), errors="coerce")
        d_ = pd.to_numeric(meta.get("os_event"), errors="coerce")
        if t_ is not None and t_.notna().sum() >= 200 and d_.sum() >= 40:
            ok = (t_.notna() & d_.notna() & (t_ > 0)).to_numpy()
            tv = t_.to_numpy()[ok]
            dv = d_.to_numpy()[ok].astype(int)
            age = pd.to_numeric(meta.get("age"), errors="coerce").to_numpy()[ok]
            node = meta.get("node").astype(str).to_numpy()[ok]
            npos = np.isin(node, ["N1", "N2", "N3"]).astype(float)
            cov = np.column_stack([np.where(np.isnan(age), np.nanmean(age), age), npos])
            for m in FOCUS:
                if m not in sc:
                    continue
                for tag, s in (("raw", sc[m]), ("resid", resid.get(m))):
                    if s is None:
                        continue
                    sv = s.to_numpy()[ok]
                    if not np.isfinite(sv).all():
                        continue
                    z = (sv - sv.mean()) / sv.std(ddof=0)
                    hr, lo, hi, p = cox_hr(tv, dv, z)
                    hr2, lo2, hi2, p2 = cox_hr(tv, dv, z, cov)
                    surv_rows.append({"cohort": name, "module": m, "scale": tag,
                                      "n": int(ok.sum()), "n_event": int(dv.sum()),
                                      "HR_per_SD": round(hr, 3), "ci_lo": round(lo, 3),
                                      "ci_hi": round(hi, 3), "p": p,
                                      "HR_adj": round(hr2, 3), "ci_lo_adj": round(lo2, 3),
                                      "ci_hi_adj": round(hi2, 3), "p_adj": p2})
            print(f"  OS n={int(ok.sum())}, events={int(dv.sum())}")

    # ---- 落盘 + FDR
    o = pd.DataFrame(orth)
    o.to_csv(OUT / "stage7_orthogonality.csv", index=False)
    a = pd.DataFrame(atten)
    a.to_csv(OUT / "stage7_attenuation.csv", index=False)
    i = pd.DataFrame(imm)
    if len(i):
        i["q_resid"] = BL.fdr_bh(i["p_resid"].to_numpy())
        i.to_csv(OUT / "stage7_residual_immune.csv", index=False)
    p_ = pd.DataFrame(pcr_rows)
    if len(p_):
        p_.to_csv(OUT / "stage7_residual_pcr.csv", index=False)
    s_ = pd.DataFrame(surv_rows)
    if len(s_):
        s_.to_csv(OUT / "stage7_residual_survival.csv", index=False)

    # ---- 控制台摘要
    print("\n" + "=" * 78)
    print("正交化：模块被 CAF/ECM 解释的比例")
    print(o.pivot_table(index="module", columns="cohort", values="pct_explained").round(1).to_string())

    if len(a):
        print("\n" + "=" * 78)
        print("三级衰减（PNI-stroma ~ 免疫排斥代理，跨队列）")
        sub = a[a["module"] == "Perineural_stroma"]
        piv = sub.pivot_table(index="target", columns="cohort",
                              values=["rho_raw", "rho_caf", "rho_caf_purity"])
        print(piv.round(2).to_string())

    if len(s_):
        print("\n" + "=" * 78)
        print("生存（HR per SD）")
        print(s_[["cohort", "module", "scale", "n", "n_event", "HR_per_SD", "p", "HR_adj", "p_adj"]]
              .to_string(index=False))
    print(f"\n输出 -> {OUT}")


if __name__ == "__main__":
    main()

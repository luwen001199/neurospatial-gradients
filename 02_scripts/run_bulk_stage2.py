# -*- coding: utf-8 -*-
"""Stage 2 —— 模块评分与免疫抑制表型的关联，并对通用基质负荷 (CAF/ECM) 做剥离。

设计要点
--------
* 每个关联都报两套: 原始 Spearman 与 控制 CAF/ECM 后的偏 Spearman。
  若控制后归零，说明该关联由 desmoplasia 共享，不是神经特异信号。
* 免疫读数一律为 mean-z 打分（本地可复现），不使用 CIBERSORTx / TIDE 在线 API；
  免疫排斥用本地 proxy EI = z(CAF_ECM) + z(M2) - z(CD8_cytotox)，并在方法中如实标注为 proxy。
* 输出: 03_results/tables/bulk_extrapolation/stage2_*.csv
"""
from __future__ import annotations

import sys
import warnings

import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, r"D:/嗜神经浸润/02_scripts")
import bulk_lib as BL  # noqa: E402

warnings.filterwarnings("ignore")
OUT = BL.TAB

FOCUS = ["NPP58", "Perineural_stroma", "Neural_intrinsic", "CAF_ECM", "LRRC15_myCAF"]
TARGETS = ["CD8_cytotox", "IFN_gamma", "Treg", "M2_macro", "SPP1_macro", "TGFb_core",
           "CXCL12_axis", "LGALS1_axis", "TLS", "MHC_II"]


def hedges_g(a: pd.Series, b: pd.Series) -> float:
    a, b = a.dropna(), b.dropna()
    if len(a) < 5 or len(b) < 5:
        return np.nan
    na, nb = len(a), len(b)
    sp = np.sqrt(((na - 1) * a.var(ddof=1) + (nb - 1) * b.var(ddof=1)) / (na + nb - 2))
    if sp == 0:
        return np.nan
    d = (a.mean() - b.mean()) / sp
    return float(d * (1 - 3 / (4 * (na + nb) - 9)))


def run(name: str, expr: pd.DataFrame, mods: dict):
    sc = BL.score_all(expr, mods)
    if "CAF_ECM" not in sc:
        print(f"[skip] {name}: 无 CAF_ECM 评分")
        return [], []
    # 免疫排斥 proxy
    if {"CAF_ECM", "M2_macro", "CD8_cytotox"} <= set(sc.columns):
        z = lambda s: (s - s.mean()) / s.std(ddof=0)
        sc["Excl_proxy"] = z(sc["CAF_ECM"]) + z(sc["M2_macro"]) - z(sc["CD8_cytotox"])

    rows, tert = [], []
    print(f"\n[{name}] n={len(sc)}")
    print(f"  {'模块':20s} {'读数':14s} {'ρ_raw':>7s} {'P':>9s} {'ρ|CAF':>7s} {'P':>9s} {'n':>5s}")
    for f in FOCUS + ["Excl_proxy"]:
        if f not in sc:
            continue
        for t in TARGETS + (["Excl_proxy"] if f != "Excl_proxy" else []):
            if t not in sc or t == f:
                continue
            cols = [f, t] + ([] if "CAF_ECM" in (f, t) else ["CAF_ECM"])
            d = sc[cols].dropna()
            if len(d) < 30:
                continue
            r_raw, p_raw = stats.spearmanr(d[f], d[t])
            if "CAF_ECM" in (f, t):
                # 控制变量与检验变量同一，偏相关无定义 → 只报原始相关
                r_par, p_par = np.nan, np.nan
            else:
                r_par, p_par, _ = BL.partial_spearman(d[f], d[t], d["CAF_ECM"])
            rows.append({"cohort": name, "module": f, "target": t, "n": len(d),
                         "rho_raw": round(r_raw, 4), "p_raw": p_raw,
                         "rho_partial": round(r_par, 4) if np.isfinite(r_par) else np.nan,
                         "p_partial": p_par})
            print(f"  {f:20s} {t:14s} {r_raw:+7.3f} {p_raw:9.2e} "
                  f"{r_par:+7.3f} {p_par:9.2e} {len(d):5d}")

        # 三分位效应量（分箱边界可能重复，逐项保护）
        try:
            q = pd.qcut(sc[f].rank(method="first"), 3, labels=["low", "mid", "high"])
        except Exception:
            q = None
        for t in TARGETS:
            if t not in sc or q is None:
                continue
            try:
                g = hedges_g(sc.loc[q == "high", t], sc.loc[q == "low", t])
                nH = int((q == "high").sum()); nL = int((q == "low").sum())
                if np.isfinite(g):
                    tert.append({"cohort": name, "module": f, "target": t,
                                 "hedges_g_high_vs_low": round(g, 3), "n_high": nH, "n_low": nL})
            except Exception:
                continue
    return rows, tert


def main():
    mods = BL.load_modules()
    all_rows, all_tert = [], []
    from run_bulk_stage1 import load_tcga, load_gse163882, load_gse25065

    # TCGA
    try:
        e, c = load_tcga()
        neg = lambda s: s.astype(str).str.lower().str.startswith(("neg", "n"))
        tn = neg(c["ER_Status_nature2012"]) & neg(c["PR_Status_nature2012"]) & \
             neg(c["HER2_Final_Status_nature2012"])
        r, t = run("TCGA-BRCA_all", e, mods); all_rows += r; all_tert += t
        r, t = run("TCGA-TNBC", e.loc[tn.fillna(False).values], mods); all_rows += r; all_tert += t
    except Exception as ex:
        import traceback; print(f"[skip] TCGA: {ex}"); traceback.print_exc()

    # GSE163882
    try:
        e, c = load_gse163882()
        tn = (c.get("estrogen receptor status", pd.Series(dtype=object)).astype(str).str.upper() == "N") & \
             (c.get("progesterone receptor status", pd.Series(dtype=object)).astype(str).str.upper() == "N") & \
             (c.get("her2 receptor status", pd.Series(dtype=object)).astype(str).str.upper() == "N")
        r, t = run("GSE163882_all", e, mods); all_rows += r; all_tert += t
        r, t = run("GSE163882_TNBC", e.loc[tn.reindex(e.index).fillna(False).values], mods)
        all_rows += r; all_tert += t
    except Exception as ex:
        import traceback; print(f"[skip] GSE163882: {ex}"); traceback.print_exc()

    # GSE25065
    try:
        e, c = load_gse25065()
        r, t = run("GSE25065_all", e, mods); all_rows += r; all_tert += t
    except Exception as ex:
        print(f"[skip] GSE25065: {ex}")

    pd.DataFrame(all_rows).to_csv(OUT / "stage2_immune_association.csv", index=False)
    pd.DataFrame(all_tert).to_csv(OUT / "stage2_tertile_effects.csv", index=False)
    print(f"\n输出 -> {OUT}")


if __name__ == "__main__":
    main()

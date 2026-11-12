# -*- coding: utf-8 -*-
"""Stage 12 —— 跨队列「随机基因集基线」正交性校准 + 跨癌种复现。

两个原稿没做、但审稿人会问的问题:

(1) 「neural 模块与 CAF/ECM 正交（ρ=0.040）」不能拿 0 当参照。任何一组
    共表达基因在这类数据里都会带一点正相关；正确的零假设是「同大小随机
    基因集」。本脚本给出每个队列中随机集的 ρ 分布，并报告观测值的 z 与
    经验 P。这同时解释了为什么 METABRIC 是 0.04 而 TCGA 是 0.20——
    如果两者都远低于各自随机基线，那「正交」是普遍性质，绝对值的差异
    只是基因检出率造成的。

(2) 嗜神经浸润在胰腺癌与头颈鳞癌中更常见，若 perineural field 是真实的
    生物学结构，它应在这些癌种里同样表现为「被 desmoplasia 主导」。

输出:
    03_results/tables/bulk_extrapolation/stage12_orthogonality_random.csv
    03_results/tables/bulk_extrapolation/stage12_pancancer.csv
"""
from __future__ import annotations

import sys

import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, r"D:/嗜神经浸润/02_scripts")
import bulk_lib as BL  # noqa: E402
from run_bulk_stage1 import load_gse163882, load_gse25065  # noqa: E402
from run_bulk_stage11_metabric import load_metabric  # noqa: E402

OUT = BL.TAB
RNG = np.random.default_rng(20261108)
N_RAND = 1000


# ------------------------------------------------------------------ 队列
def load_tcga(study: str = "brca_tcga"):
    expr = pd.read_csv(BL.RAW / study / f"{study}_panel.csv", index_col=0).astype(float)
    clin = pd.read_csv(BL.RAW / study / f"{study}_clin.csv", index_col=0).reindex(expr.index)
    t = pd.to_numeric(clin.get("OS_MONTHS"), errors="coerce")
    e = clin.get("OS_STATUS", pd.Series(index=expr.index, dtype=object))
    out = pd.DataFrame({"OS_MONTHS": t,
                        "OS_EVENT": e.astype(str).str.startswith("1").astype(float),
                        "AGE": pd.to_numeric(clin.get("AGE"), errors="coerce"),
                        "STAGE": clin.get("AJCC_PATHOLOGIC_TUMOR_STAGE")})
    out.loc[e.isna(), "OS_EVENT"] = np.nan
    return expr, out


def tnbc_mask(clin_index):
    xena = BL.RAW / "tcga_brca" / "BRCA_clinicalMatrix"
    if not xena.exists():
        return pd.Series(False, index=clin_index)
    c = pd.read_csv(xena, sep="\t", index_col=0, low_memory=False)
    neg = lambda s: s.astype(str).str.lower().str.startswith(("neg", "n"))
    m = (neg(c["ER_Status_nature2012"]) & neg(c["PR_Status_nature2012"])
         & neg(c["HER2_Final_Status_nature2012"]))
    return m.reindex(clin_index).fillna(False)


def build_cohorts() -> dict[str, tuple[pd.DataFrame, pd.DataFrame]]:
    co = {}
    e, c = load_tcga("brca_tcga")
    co["TCGA-BRCA"] = (e, c)
    tn = tnbc_mask(e.index)
    co["TCGA-TNBC"] = (e.loc[tn], c.loc[tn])
    try:
        e2, c2 = load_gse163882()
        co["GSE163882"] = (e2, pd.DataFrame(index=e2.index))
    except Exception as ex:
        print(f"[skip] GSE163882: {ex}")
    try:
        e3, c3 = load_gse25065()
        co["GSE25065"] = (e3, pd.DataFrame(index=e3.index))
    except Exception as ex:
        print(f"[skip] GSE25065: {ex}")
    try:
        e4, c4 = load_metabric()
        co["METABRIC"] = (e4, c4)
    except Exception as ex:
        print(f"[skip] METABRIC: {ex}")
    for st, tag in [("paad_tcga", "TCGA-PAAD"), ("hnsc_tcga", "TCGA-HNSC")]:
        p = BL.RAW / st / f"{st}_panel.csv"
        if p.exists():
            co[tag] = load_tcga(st)
        else:
            print(f"[wait] {tag} 面板未就绪")
    return co


# -------------------------------------------------- 随机集校准的正交性
def orthogonality_vs_random(expr: pd.DataFrame, mods: dict, name: str):
    caf = BL.score_mean(expr, mods["CAF_ECM"])
    genes_pool = list(expr.columns)
    rows = []
    for mod in ["Perineural_stroma", "Neural_intrinsic", "CD8_cytotox", "TGFb_core"]:
        gl = [g for g in mods[mod] if g in expr.columns]
        if len(gl) < 6:
            continue
        s = BL.score_mean(expr, gl)
        rho = stats.spearmanr(s, caf, nan_policy="omit")[0]
        pool = [g for g in genes_pool if g not in set(gl)]
        k = len(gl)
        null = np.empty(N_RAND)
        for i in range(N_RAND):
            sub = RNG.choice(pool, k, replace=False)
            null[i] = stats.spearmanr(BL.score_mean(expr, list(sub)), caf,
                                      nan_policy="omit")[0]
        mu, sd = float(np.nanmean(null)), float(np.nanstd(null))
        z = (rho - mu) / sd if sd > 0 else np.nan
        # 双侧经验 P（观测值在零分布中的位置）
        p_emp = float(min((null <= rho).mean(), (null >= rho).mean()) * 2)
        rows.append({"cohort": name, "module": mod, "n_genes": k, "n_sample": len(expr),
                     "rho_obs": round(float(rho), 4),
                     "null_median": round(mu, 4), "null_p2.5": round(float(np.percentile(null, 2.5)), 4),
                     "null_p97.5": round(float(np.percentile(null, 97.5)), 4),
                     "z_vs_null": round(float(z), 2), "p_empirical": p_emp})
        print(f"  {name:12s} {mod:20s} ρ={rho:+.3f} | 随机集 {mu:+.3f} "
              f"[{np.percentile(null,2.5):+.3f},{np.percentile(null,97.5):+.3f}] "
              f"z={z:+.1f} P={p_emp:.3g}")
    return pd.DataFrame(rows)


# -------------------------------------------------- 跨癌种主链
def pancancer_chain(expr: pd.DataFrame, clin: pd.DataFrame, mods: dict, name: str):
    logx = np.log1p(expr.clip(lower=0))
    sc = BL.score_all(logx, mods)
    if "CD8_cytotox" not in sc:
        sc["CD8_cytotox"] = BL.score_mean(logx, mods["CD8_cytotox"])
    resid, r2 = BL.rank_resid(sc["Perineural_stroma"], sc["CAF_ECM"])
    sc["Stroma_resid"] = resid
    r_tgf = p_tgf = np.nan
    if "TGFb_core" in sc:
        r_tgf, p_tgf, _ = BL.partial_spearman(sc["Stroma_resid"], sc["TGFb_core"], sc["CAF_ECM"])
    r_cd8, p_cd8, _ = BL.partial_spearman(sc["Stroma_resid"], sc["CD8_cytotox"], sc["CAF_ECM"])
    out = {"cohort": name, "n": len(expr), "r2_stroma_on_CAF": round(r2, 4),
           "rho_neural_CAF": round(float(stats.spearmanr(sc["Neural_intrinsic"],
                                                         sc["CAF_ECM"], nan_policy="omit")[0]), 3),
           "rho_resid_TGFb": round(float(r_tgf), 3), "p_resid_TGFb": p_tgf,
           "rho_resid_CD8": round(float(r_cd8), 3), "p_resid_CD8": p_cd8}
    if "OS_MONTHS" in clin and clin["OS_MONTHS"].notna().sum() > 80:
        d = pd.concat([sc[["Stroma_resid", "Perineural_stroma"]],
                       clin[["OS_MONTHS", "OS_EVENT"]]], axis=1).dropna()
        d = d[d["OS_MONTHS"] > 0]
        if len(d) > 80 and d["OS_EVENT"].sum() > 10:
            for key, col in [("HR_resid", "Stroma_resid"), ("HR_raw", "Perineural_stroma")]:
                z = (d[col] - d[col].mean()) / d[col].std(ddof=0)
                hr, lo, hi, p = BL.cox_hr(d["OS_MONTHS"].to_numpy(),
                                          d["OS_EVENT"].to_numpy(), z.to_numpy())
                out[key] = round(hr, 3); out[key + "_lo"] = round(lo, 3)
                out[key + "_hi"] = round(hi, 3); out[key + "_p"] = p
            out["n_os"] = len(d); out["events"] = int(d["OS_EVENT"].sum())
    print(f"  {name:12s} n={out['n']:4d} | CAF 解释 {r2*100:4.1f}% | "
          f"resid~TGF-β {r_tgf:+.3f} | resid~CD8 {r_cd8:+.3f}"
          + (f" | HR_resid {out['HR_resid']:.2f}" if "HR_resid" in out else ""))
    return out


def main():
    mods = BL.load_modules()
    co = build_cohorts()
    print(f"队列 {list(co)}")

    print("\n[1] 正交性 vs 随机基因集基线")
    orth = pd.concat([orthogonality_vs_random(e, mods, n) for n, (e, _) in co.items()],
                     ignore_index=True)
    orth["q"] = BL.fdr_bh(orth["p_empirical"].to_numpy())
    orth.to_csv(OUT / "stage12_orthogonality_random.csv", index=False)

    print("\n[2] 跨癌种主链")
    rows = [pancancer_chain(e, c, mods, n) for n, (e, c) in co.items()]
    pc = pd.DataFrame(rows)
    pc.to_csv(OUT / "stage12_pancancer.csv", index=False)
    print(f"\n输出 -> {OUT}")


if __name__ == "__main__":
    main()

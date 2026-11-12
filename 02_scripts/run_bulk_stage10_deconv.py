# -*- coding: utf-8 -*-
"""Stage 10 —— 免疫区室富集 + 免疫抑制链的正式中介分析。

原稿只在空间层看过 CD8；bulk 层从未做免疫组成/富集关联，也没有把
「myCAF 场 → TGF-β/巨噬 → CD8 抑制」量化成中介效应。本脚本补这两块。

反卷积可行性（已在 2026-11-08 实测，记录在此以免重复尝试）:
    256 基因靶向面板（ECM/神经为主）对 10 区室线性混合不可识别。
    - 线性 CP10K + 截距：免疫区室权重全 0，Fibroblast 0.56 / Epithelial 0.32；
    - sqrt 变换：同样退化；
    - 基因维 z 标度化：权重退化为均匀分布（R²=0.003）。
    结论：不做 CIBERSORTx 式比例反卷积，改用区室特异标记集富集分
    （xCell/Danaher 式），并用配对单细胞 pseudobulk 验证其能否读出真实占比。

输出:
    03_results/tables/bulk_extrapolation/stage10_score_validation.csv
    03_results/tables/bulk_extrapolation/stage10_deconv_composition.csv
    03_results/tables/bulk_extrapolation/stage10_mediation.csv
    03_results/tables/bulk_extrapolation/stage10_deconv_sample.csv
"""
from __future__ import annotations

import sys

import numpy as np
import pandas as pd
from scipy import sparse, stats
from scipy.io import mmread

sys.path.insert(0, r"D:/嗜神经浸润/02_scripts")
import bulk_lib as BL  # noqa: E402

OUT = BL.TAB
SRC = BL.ROOT / "01_data" / "external" / "gse176078"
N_BOOT = 2000
RNG = np.random.default_rng(20261108)


# ------------------------------------------------------------------ 数据
def load_brca():
    expr = pd.read_csv(BL.RAW / "brca_tcga" / "brca_tcga_panel.csv", index_col=0)
    clin = pd.read_csv(BL.RAW / "brca_tcga" / "brca_tcga_clin.csv", index_col=0)
    return expr.astype(float), clin.reindex(expr.index)


def tnbc_mask(clin: pd.DataFrame) -> pd.Series:
    xena = BL.RAW / "tcga_brca" / "BRCA_clinicalMatrix"
    if not xena.exists():
        return pd.Series(False, index=clin.index)
    c = pd.read_csv(xena, sep="\t", index_col=0, low_memory=False)
    neg = lambda s: s.astype(str).str.lower().str.startswith(("neg", "n"))
    m = (neg(c["ER_Status_nature2012"]) & neg(c["PR_Status_nature2012"])
         & neg(c["HER2_Final_Status_nature2012"]))
    return m.reindex(clin.index).fillna(False)


def scrna_validation(markers: dict[str, list[str]]) -> pd.DataFrame:
    """26 例单细胞 pseudobulk：区室富集分 vs 该区室真实细胞占比（Spearman）。"""
    genes = sorted({g for v in markers.values() for g in v})
    rows = []
    for d in sorted([p for p in SRC.iterdir()
                     if p.is_dir() and (p / "count_matrix_sparse.mtx").exists()]):
        with open(d / "count_matrix_sparse.mtx", "rb") as fh:
            X = sparse.csc_matrix(mmread(fh))
        gl = pd.read_csv(d / "count_matrix_genes.tsv", header=None)[0].astype(str).tolist()
        gi = {g: i for i, g in enumerate(gl)}
        keep = [(g, gi[g]) for g in genes if g in gi]
        meta = pd.read_csv(d / "metadata.csv", index_col=0)
        meta.index = meta.index.astype(str)
        bar = pd.read_csv(d / "count_matrix_barcodes.tsv", header=None)[0].astype(str).tolist()
        if len(bar) != X.shape[1]:
            continue
        meta = meta.reindex(bar)
        tot = np.asarray(X.sum(axis=0)).ravel().astype(float)
        tot[tot <= 0] = np.nan
        sub = np.asarray(X[[i for _, i in keep], :].todense(), float) / tot * 1e4
        rec = {"sample": d.name, "n_cells": len(bar)}
        for g, v in zip([g for g, _ in keep], np.nanmean(sub, axis=1)):
            rec[g] = np.log1p(v)
        ct = meta["celltype_minor"].astype(str).map(BL.MINOR_TO_COMPARTMENT)
        for comp in markers:
            rec[f"frac_{comp}"] = float((ct == comp).mean())
        rows.append(rec)
    df = pd.DataFrame(rows).set_index("sample")
    sc = BL.compartment_scores(df, pd.read_csv(BL.PROC / "deconv_signature_gse176078.csv",
                                               index_col=0))
    out = []
    for comp in sc.columns:
        f = f"frac_{comp}"
        if f not in df:
            continue
        d = pd.concat([sc[comp].rename("s"), df[f].rename("f")], axis=1).dropna()
        if len(d) < 8:
            continue
        rho, p = stats.spearmanr(d["s"], d["f"])
        out.append({"compartment": comp, "n_genes": len(markers.get(comp, [])),
                    "n_sample": len(d), "rho_score_vs_true_fraction": round(rho, 3),
                    "p": p})
    return pd.DataFrame(out)


# ------------------------------------------------------------------ 中介
def normal_score(s: pd.Series) -> np.ndarray:
    """秩 → 正态分数。用 (rank-0.5)/n 而非 rank/n，否则最大秩 ppf(1)=+inf。"""
    r = s.rank(method="average")
    return stats.norm.ppf((r - 0.5) / len(r))


def mediation(x: pd.Series, m: pd.Series, y: pd.Series, n_boot=N_BOOT):
    """product-of-coefficients 中介；三变量均转正态分数，OLS 取系数，bootstrap 取 ab CI。"""
    d = pd.concat([x.rename("x"), m.rename("m"), y.rename("y")], axis=1)
    d = d.replace([np.inf, -np.inf], np.nan).dropna()
    if len(d) < 40:
        return None
    X, M, Y = (normal_score(d[c]) for c in ("x", "m", "y"))
    n = len(d)

    def fit(xv, mv, yv):
        a = np.polyfit(xv, mv, 1)[0]
        A = np.column_stack([np.ones(len(xv)), xv, mv])
        coef, *_ = np.linalg.lstsq(A, yv, rcond=None)
        c = np.polyfit(xv, yv, 1)[0]
        return a, coef[2], c, coef[1]

    a, b, c, c_p = fit(X, M, Y)
    boots = np.empty(n_boot)
    for i in range(n_boot):
        idx = RNG.integers(0, n, n)
        ai, bi, _, _ = fit(X[idx], M[idx], Y[idx])
        boots[i] = ai * bi
    lo, hi = np.percentile(boots, [2.5, 97.5])
    return {"n": n, "a": a, "b": b, "c_total": c, "c_prime": c_p,
            "ab": a * b, "ab_lo": lo, "ab_hi": hi,
            "prop_mediated": (a * b) / c if abs(c) > 1e-9 else np.nan,
            "p_boot": float(2 * min((boots <= 0).mean(), (boots >= 0).mean()))}


# ------------------------------------------------------------------ 主流程
def run_cohort(name: str, expr: pd.DataFrame, mods: dict, S: pd.DataFrame):
    print(f"\n{'='*72}\n[{name}] n={len(expr)}")
    logx = np.log1p(expr.clip(lower=0))
    comp = BL.compartment_scores(logx, S)
    comps = list(comp.columns)
    print("  区室富集分: " + ", ".join(f"{c}({comp[c].median():+.2f})" for c in comps))

    sc = BL.score_all(logx, mods)
    sc["CD8_cytotox"] = BL.score_mean(logx, mods["CD8_cytotox"])
    resid, r2 = BL.rank_resid(sc["Perineural_stroma"], sc["CAF_ECM"])
    sc["Stroma_resid"] = resid
    print(f"  模块残差: 被 CAF/ECM 解释 {r2*100:.1f}% 秩方差")

    rows = []
    for c in comps:
        for mod in ["Perineural_stroma", "Stroma_resid", "CAF_ECM", "Neural_intrinsic"]:
            d = pd.concat([comp[c].rename("c"), sc[mod].rename("m")], axis=1).dropna()
            if len(d) < 30:
                continue
            rho, p = stats.spearmanr(d["c"], d["m"])
            q = pd.qcut(d["m"], 3, labels=["low", "mid", "high"], duplicates="drop")
            g = (d.loc[q == "high", "c"].mean() - d.loc[q == "low", "c"].mean()) / d["c"].std(ddof=0)
            rows.append({"cohort": name, "compartment": c, "module": mod, "n": len(d),
                         "rho": round(rho, 4), "p": p, "delta_high_minus_low_SD": round(g, 3)})
    comp_df = pd.DataFrame(rows)
    comp_df["q"] = BL.fdr_bh(comp_df["p"].to_numpy())
    print(f"\n  区室 ~ 模块（BH q<0.05，按 |ρ| 排序）")
    for mod in ["Perineural_stroma", "Stroma_resid"]:
        sub = comp_df[(comp_df.module == mod) & (comp_df.q < 0.05)]
        sub = sub.reindex(sub["rho"].abs().sort_values(ascending=False).index)
        for _, r in sub.head(8).iterrows():
            print(f"    {mod:18s} {r['compartment']:12s} ρ={r['rho']:+.3f} "
                  f"q={r['q']:.3g} Δ={r['delta_high_minus_low_SD']:+.2f} SD")

    med_rows = []
    for mv, label in [("TGFb_core", "TGF-β program"),
                      ("Macrophage", "macrophage enrichment"),
                      ("Treg", "Treg marker score"),
                      ("M2_macro", "M2 marker score"),
                      ("CD8_T", "CD8 T enrichment")]:
        M = comp[mv] if mv in comp.columns else sc.get(mv)
        if M is None:
            continue
        r = mediation(sc["Stroma_resid"], M, sc["CD8_cytotox"])
        if r is None:
            continue
        r.update({"cohort": name, "mediator": label,
                  "exposure": "perineural-stroma residual", "outcome": "CD8 cytotoxicity"})
        med_rows.append(r)
        print(f"  中介 {label:24s} a={r['a']:+.3f} b={r['b']:+.3f} c={r['c_total']:+.3f} "
              f"c'={r['c_prime']:+.3f} ab={r['ab']:+.3f} "
              f"[{r['ab_lo']:+.3f},{r['ab_hi']:+.3f}] prop={r['prop_mediated']:+.2f}")
    rev = mediation(sc["CD8_cytotox"], sc["TGFb_core"], sc["Stroma_resid"])
    if rev:
        rev.update({"cohort": name, "mediator": "TGF-β program",
                    "exposure": "CD8 cytotoxicity", "outcome": "perineural-stroma residual"})
        med_rows.append(rev)
        print(f"  反向   CD8 → TGF-β → 残差      ab={rev['ab']:+.3f} "
              f"[{rev['ab_lo']:+.3f},{rev['ab_hi']:+.3f}]")

    ps = pd.concat([comp, sc[["Perineural_stroma", "Stroma_resid", "CAF_ECM",
                              "Neural_intrinsic", "CD8_cytotox", "TGFb_core",
                              "M2_macro", "Treg"]]], axis=1)
    ps.insert(0, "cohort", name)
    return comp_df, pd.DataFrame(med_rows), ps


def main():
    mods = BL.load_modules()
    S = pd.read_csv(BL.PROC / "deconv_signature_gse176078.csv", index_col=0)
    markers = BL.specific_markers(S)
    print("区室特异标记数: " + ", ".join(f"{k}={len(v)}" for k, v in markers.items()))

    print("\n[验证] 单细胞 pseudobulk：富集分能否读出真实细胞占比（n=26 例）")
    val = scrna_validation(markers)
    print(val.to_string(index=False))
    val.to_csv(OUT / "stage10_score_validation.csv", index=False)

    expr, clin = load_brca()
    tn = tnbc_mask(clin)
    print(f"\nTCGA-BRCA n={len(expr)}; TNBC n={int(tn.sum())}")

    parts = [run_cohort(n, e, mods, S)
             for n, e in [("TCGA-BRCA", expr), ("TCGA-TNBC", expr.loc[tn])]]
    comp_df = pd.concat([p[0] for p in parts], ignore_index=True)
    med_df = pd.concat([p[1] for p in parts], ignore_index=True)
    ps_df = pd.concat([p[2] for p in parts])
    comp_df.to_csv(OUT / "stage10_deconv_composition.csv", index=False)
    med_df.to_csv(OUT / "stage10_mediation.csv", index=False)
    ps_df.to_csv(OUT / "stage10_deconv_sample.csv")
    print(f"\n输出 -> {OUT}")


if __name__ == "__main__":
    main()

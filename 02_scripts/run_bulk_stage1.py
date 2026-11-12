# -*- coding: utf-8 -*-
"""Stage 1 —— 空间先验模块向 bulk 的迁移与正交性复现。

判定标准（方案 Go/No-go）:
    perineural-stroma 模块 vs 通用 CAF/ECM 的 Spearman ρ ≥ 0.5 → 迁移成立
    neural-intrinsic vs CAF/ECM 的 ρ ≈ 0        → 正交性复现

输出: 03_results/tables/bulk_extrapolation/stage1_*.csv
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, r"D:/嗜神经浸润/02_scripts")
import bulk_lib as BL  # noqa: E402

RAW = BL.RAW
OUT = BL.TAB


# ------------------------------------------------------------------ 队列装载
def load_tcga() -> tuple[pd.DataFrame, pd.DataFrame]:
    """TCGA-BRCA: 表达走 cBioPortal per-gene 缓存, 临床走 Xena clinicalMatrix。"""
    expr = pd.read_csv(RAW / "tcga_brca" / "tcga_brca_target_genes.csv", index_col=0)
    expr.index.name = "sample"
    clin = pd.read_csv(RAW / "tcga_brca" / "tcga_brca_clin.csv", index_col=0, low_memory=False)
    clin = clin.reindex(expr.index)
    # 统一列名到 Stage 脚本其余部分期待的名字
    rename = {"ER": "ER_Status_nature2012", "PR": "PR_Status_nature2012",
              "HER2": "HER2_Final_Status_nature2012"}
    clin = clin.rename(columns={k: v for k, v in rename.items() if k in clin.columns})
    return expr.astype(float), clin


def load_gse163882() -> tuple[pd.DataFrame, pd.DataFrame]:
    """GSE163882: AmpliSeq 靶向 RNA-seq TPM（Ensembl 行 × BA 样本列）+ series matrix 标签。

    TPM 矩阵末列 'annotation' 存 gene symbol；样本列名为 BA00111 等院内编号，
    与 series matrix 的 !Sample_description 一致，据此对齐临床标签。
    """
    p = RAW / "gse163882" / "GSE163882_all.data.tpms_222Samples.csv.gz"
    df = pd.read_csv(p, index_col=0)
    sym = df.pop("annotation") if "annotation" in df.columns else None
    df = df.apply(pd.to_numeric, errors="coerce")
    if sym is not None:
        df.index = sym.astype(str).str.strip()
        df = df[~df.index.isin(["", "nan"])]
        df = df.groupby(level=0).mean()
    expr = df.T
    expr.index = [str(i).strip() for i in expr.index]

    mp = RAW / "gse163882" / "GSE163882_series_matrix.txt.gz"
    _, meta = BL.series_matrix_expr(mp)
    if "description" in meta.columns:
        meta = meta.set_index(meta["description"].astype(str).str.strip())
    meta.index = [str(i).strip() for i in meta.index]
    clin = meta.reindex(expr.index)
    return expr.astype(float), clin


def load_gse25065() -> tuple[pd.DataFrame, pd.DataFrame]:
    mp = RAW / "gse25065" / "GSE25065_series_matrix.txt.gz"
    expr_p, meta = BL.series_matrix_expr(mp)
    annot = BL.read_gpl_annot(RAW / "gse12276" / "GPL570.annot.gz")
    m = annot.set_index("probe")["symbol"].to_dict()
    hit = np.mean([p in m for p in expr_p.index])
    print(f"[GSE25065] probe 可映射率 (借 GPL570 注释, U133A ⊂ U133Plus2) = {hit:.1%}")
    expr_p = expr_p.loc[[p for p in expr_p.index if p in m]]
    expr_p.index = [m[p] for p in expr_p.index]
    expr = expr_p.groupby(level=0).mean().T
    clin = meta.reindex(expr.index)
    return expr.astype(float), clin


# ------------------------------------------------------------------ 主流程
def describe(name: str, expr: pd.DataFrame, mask: pd.Series | None) -> dict:
    e = expr if mask is None else expr.loc[mask.reindex(expr.index).fillna(False).values]
    return {"cohort": name, "n": len(e), "n_genes": e.shape[1]}


def corr_block(sc: pd.DataFrame, focus=("Perineural_stroma", "Neural_intrinsic", "NPP58"),
               ref=("CAF_ECM", "CD8_cytotox", "TLS")) -> pd.DataFrame:
    rows = []
    for f in focus:
        for r in ref:
            if f not in sc or r not in sc:
                continue
            d = sc[[f, r]].dropna()
            rho, p = stats.spearmanr(d[f], d[r])
            lo, hi = BL.boot_ci(d[f], d[r], n=1000)
            rows.append({"module": f, "reference": r, "n": len(d),
                         "rho": round(rho, 4), "p": p, "ci_lo": round(lo, 4), "ci_hi": round(hi, 4)})
    return pd.DataFrame(rows)


def partial_block(sc: pd.DataFrame, focus=("Perineural_stroma", "Neural_intrinsic", "NPP58"),
                  targets=("CD8_cytotox", "Treg", "M2_macro", "SPP1_macro", "TGFb_core",
                           "CXCL12_axis", "LGALS1_axis", "IFN_gamma", "TLS")) -> pd.DataFrame:
    rows = []
    if "CAF_ECM" not in sc:
        return pd.DataFrame()
    for f in focus:
        for t in targets:
            if f not in sc or t not in sc:
                continue
            raw_r, raw_p, n1 = BL.partial_spearman(sc[f], sc[t], pd.Series(0.0, index=sc.index))
            pr, pp, n2 = BL.partial_spearman(sc[f], sc[t], sc["CAF_ECM"])
            rows.append({"module": f, "target": t, "n": n2,
                         "rho_raw": round(raw_r, 4), "p_raw": raw_p,
                         "rho_partial_CAFECM": round(pr, 4), "p_partial": pp})
    return pd.DataFrame(rows)


def main():
    mods = BL.load_modules()
    cohorts = {}

    print("=" * 70)
    print("装载队列")
    try:
        e, c = load_tcga()
        cohorts["TCGA-BRCA"] = (e, c)
        print(f"  TCGA-BRCA: {e.shape[0]} 样本 × {e.shape[1]} 基因")
    except Exception as ex:
        print(f"  [skip] TCGA: {ex}")

    try:
        e, c = load_gse163882()
        cohorts["GSE163882"] = (e, c)
        print(f"  GSE163882: {e.shape[0]} 样本 × {e.shape[1]} 基因")
    except Exception as ex:
        print(f"  [skip] GSE163882: {ex}")

    try:
        e, c = load_gse25065()
        cohorts["GSE25065"] = (e, c)
        print(f"  GSE25065: {e.shape[0]} 样本 × {e.shape[1]} 基因")
    except Exception as ex:
        print(f"  [skip] GSE25065: {ex}")

    # ---- TNBC 子集
    subsets = {}
    if "TCGA-BRCA" in cohorts:
        e, c = cohorts["TCGA-BRCA"]
        def neg(s):
            return s.astype(str).str.lower().str.startswith(("neg", "n"))
        tn = neg(c["ER_Status_nature2012"]) & neg(c["PR_Status_nature2012"]) & \
             neg(c["HER2_Final_Status_nature2012"])
        subsets["TCGA-BRCA_all"] = (e, c, pd.Series(True, index=e.index))
        subsets["TCGA-TNBC"] = (e, c, tn)
        print(f"  TCGA TNBC n={int(tn.sum())} (总 n={len(e)})")

    if "GSE163882" in cohorts:
        e, c = cohorts["GSE163882"]
        tn = (c.get("estrogen receptor status") == "N") & \
             (c.get("progesterone receptor status") == "N") & \
             (c.get("her2 receptor status") == "N")
        subsets["GSE163882_all"] = (e, c, pd.Series(True, index=e.index))
        subsets["GSE163882-TNBC"] = (e, c, tn)
        print(f"  GSE163882 TNBC n={int(tn.sum())}")

    if "GSE25065" in cohorts:
        e, c = cohorts["GSE25065"]
        cols = {x.lower(): x for x in c.columns}
        er = c.get(cols.get("er_status_ihc", ""), pd.Series("", index=e.index))
        pr = c.get(cols.get("pr_status_ihc", ""), pd.Series("", index=e.index))
        h2 = c.get(cols.get("her2_status", ""), pd.Series("", index=e.index))
        tn = er.astype(str).str.upper().str.startswith("N") & \
             pr.astype(str).str.upper().str.startswith("N") & \
             h2.astype(str).str.upper().str.startswith("N")
        subsets["GSE25065_all"] = (e, c, pd.Series(True, index=e.index))
        subsets["GSE25065-TNBC"] = (e, c, tn)
        print(f"  GSE25065 TNBC n={int(tn.sum())} (总 n={len(e)})")

    # ---- 评分 + 相关
    all_corr, all_part, all_score = [], [], {}
    for name, (e, c, mask) in subsets.items():
        e2 = e.loc[mask.fillna(False).values]
        if len(e2) < 20:
            print(f"  [skip] {name}: n={len(e2)} 过小")
            continue
        sc = BL.score_all(e2, mods)
        all_score[name] = sc
        cb = corr_block(sc)
        cb.insert(0, "cohort", name)
        all_corr.append(cb)
        pb = partial_block(sc)
        if len(pb):
            pb.insert(0, "cohort", name)
            all_part.append(pb)
        print(f"  {name}: n={len(e2)}, 可用模块 {sc.shape[1]}")

    corr = pd.concat(all_corr, ignore_index=True)
    corr.to_csv(OUT / "stage1_module_correlations.csv", index=False)
    if all_part:
        part = pd.concat(all_part, ignore_index=True)
        part.to_csv(OUT / "stage1_partial_correlations.csv", index=False)

    # ---- Go/No-go
    print("\n" + "=" * 70)
    print("Stage 1 Go/No-go")
    for name, sc in all_score.items():
        if "Perineural_stroma" not in sc or "CAF_ECM" not in sc:
            continue
        r1 = stats.spearmanr(sc["Perineural_stroma"], sc["CAF_ECM"], nan_policy="omit")[0]
        r2 = stats.spearmanr(sc["Neural_intrinsic"], sc["CAF_ECM"], nan_policy="omit")[0]
        print(f"  {name:20s}  stroma~CAF/ECM ρ={r1:+.3f}  neural~CAF/ECM ρ={r2:+.3f}  "
              f"→ {'PASS' if r1 >= 0.5 else 'FAIL'}")

    # 汇总
    with open(OUT / "stage1_summary.md", "w", encoding="utf-8") as f:
        f.write("# Stage 1 — 模块迁移与正交性\n\n")
        f.write("## 队列规模\n\n")
        for name, sc in all_score.items():
            f.write(f"- {name}: n={len(sc)}\n")
        f.write("\n## 相关（Spearman, 1000× bootstrap 95% CI）\n\n")
        f.write("```\n" + corr.to_string(index=False) + "\n```\n")
        f.write("\n")
    print(f"\n输出 -> {OUT}")


if __name__ == "__main__":
    main()

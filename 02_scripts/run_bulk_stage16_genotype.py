# -*- coding: utf-8 -*-
"""Stage 16 — TCGA-BRCA 基因型桥: 模块 × 突变 landscape。

问题: bulk 主链全部停留在转录组表型; 基因组维度 (TMB / 驱动基因突变)
与 perineural stroma 模块的关系从未检验。

数据 (本地): brca_tcga_panel.csv (1093 样本 × 258 基因表达),
  brca_tcga_clin.csv (OS), brca_tcga_mutations.csv (982 样本:
  mutation_count + TP53/BRCA1/BRCA2/PIK3CA 状态, cBioPortal DETAILED fetch,
  MUTATION_COUNT 与逐条非同义计数 ρ=1.000 交叉验证)。

分析:
  1) TMB (log10 count+1) ~ 模块 Spearman, 校正纯度代理 + MKI67 的偏相关
  2) TP53 / PIK3CA / BRCA1 / BRCA2 突变 vs WT 的模块分布
     (Mann-Whitney + Cliff's delta + BH)
  3) 敏感性: stroma × TP53 对 OS 的交互 (Cox, z 标准化)

输出: stage16_tmb_correlation.csv / stage16_mutation_status.csv /
      stage16_tp53_os.csv
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
RAW = BL.RAW
MODULES = ["NPP58", "Perineural_stroma", "Neural_intrinsic", "CAF_ECM",
           "LRRC15_myCAF", "CD8_cytotox", "TGFb_core", "TLS", "IFN_gamma"]
PURITY = ["EPCAM", "KRT8", "KRT18", "KRT19", "KRT7", "CDH1", "MUC1"]


def zstd(s: pd.Series) -> pd.Series:
    return (s - s.mean()) / s.std(ddof=0)


def cliffs_delta(a: np.ndarray, b: np.ndarray) -> float:
    n1, n2 = len(a), len(b)
    if n1 == 0 or n2 == 0:
        return np.nan
    d = 0.0
    for x in a:
        d += (x > b).sum() - (x < b).sum()
    return d / (n1 * n2)


def main() -> None:
    expr = pd.read_csv(RAW / "brca_tcga" / "brca_tcga_panel.csv", index_col=0).astype(float)
    clin = pd.read_csv(RAW / "brca_tcga" / "brca_tcga_clin.csv", index_col=0)
    mut = pd.read_csv(RAW / "brca_tcga" / "brca_tcga_mutations.csv", index_col="sampleId")
    clin = clin.reindex(expr.index)
    mut = mut.reindex(expr.index)

    mods = BL.load_modules()
    sc = BL.score_all(expr, mods)
    purity = BL.score_mean(expr, [g for g in PURITY if g in expr.columns])
    mki67 = zstd(np.log1p(expr["MKI67"])) if "MKI67" in expr.columns else pd.Series(np.nan, index=expr.index)

    tmb = np.log10(mut["mutation_count"].astype(float).to_numpy() + 1)
    has_tmb = np.isfinite(tmb)
    print(f"[数据] 样本 {len(expr)}, 有突变数据 {has_tmb.sum()}, "
          f"TP53 mut {int(mut['TP53'].fillna(0).sum())}, "
          f"PIK3CA {int(mut['PIK3CA'].fillna(0).sum())}, "
          f"BRCA1/2 {int(mut['BRCA1'].fillna(0).sum())}/{int(mut['BRCA2'].fillna(0).sum())}")

    # ---- 1. TMB 相关
    rows = []
    for m in MODULES:
        if m not in sc:
            continue
        s = zstd(sc[m])
        ok = has_tmb & s.notna().to_numpy()
        rho, p = stats.spearmanr(tmb[ok], s.to_numpy()[ok])
        pr, pp, _ = BL.partial_spearman(pd.Series(tmb[ok]), s.to_numpy()[ok],
                                        pd.DataFrame({"purity": purity.to_numpy()[ok],
                                                      "mki67": mki67.to_numpy()[ok]}))
        rows.append({"module": m, "n": int(ok.sum()), "rho": round(rho, 4), "p": p,
                     "rho_partial_purity_mki67": round(pr, 4), "p_partial": pp})
    t1 = pd.DataFrame(rows)
    t1["q"] = BL.fdr_bh(t1["p"].to_numpy())
    t1.to_csv(OUT / "stage16_tmb_correlation.csv", index=False)
    print("\n[1] TMB 相关:")
    print(t1.to_string(index=False))

    # ---- 2. 突变状态 vs 模块
    rows = []
    for gene in ["TP53", "PIK3CA", "BRCA1", "BRCA2"]:
        g = mut[gene].fillna(0).to_numpy().astype(int)
        for m in MODULES:
            if m not in sc:
                continue
            s = zstd(sc[m]).to_numpy()
            ok = np.isfinite(s)
            a, b = s[ok & (g == 1)], s[ok & (g == 0)]
            if len(a) < 5:
                continue
            u, p = stats.mannwhitneyu(a, b, alternative="two-sided")
            rows.append({"gene": gene, "n_mut": len(a), "n_wt": len(b), "module": m,
                         "median_mut": round(float(np.median(a)), 3),
                         "median_wt": round(float(np.median(b)), 3),
                         "cliffs_delta": round(cliffs_delta(a, b), 3),
                         "p": p})
    t2 = pd.DataFrame(rows)
    t2["q"] = BL.fdr_bh(t2["p"].to_numpy())
    t2.to_csv(OUT / "stage16_mutation_status.csv", index=False)
    print("\n[2] 突变状态 (名义 P<0.05 的行):")
    print(t2[t2.p < 0.05].to_string(index=False))

    # BRCA1+2 合并
    g = ((mut["BRCA1"].fillna(0) + mut["BRCA2"].fillna(0)) > 0).to_numpy().astype(int)
    for m in ["Perineural_stroma", "CD8_cytotox", "TGFb_core", "LRRC15_myCAF"]:
        if m not in sc:
            continue
        s = zstd(sc[m]).to_numpy()
        ok = np.isfinite(s)
        a, b = s[ok & (g == 1)], s[ok & (g == 0)]
        if len(a) >= 5:
            u, p = stats.mannwhitneyu(a, b)
            print(f"  [BRCA1/2 合并 n={int(g.sum())}] {m:18s} Δ={cliffs_delta(a, b):+.3f} P={p:.3g}")

    # ---- 3. TP53 × stroma 对 OS 的交互
    osm = clin["OS_MONTHS"].astype(float).to_numpy()
    osd = clin["OS_STATUS"].astype(str).str.startswith("1:").to_numpy().astype(int)
    tp53 = mut["TP53"].fillna(0).to_numpy().astype(int)
    stroma = zstd(sc["Perineural_stroma"]).to_numpy()
    ok = np.isfinite(osm) & (osm > 0) & np.isfinite(stroma)
    inter = stroma * (tp53 - tp53.mean())
    iv = np.column_stack([stroma, tp53, inter])
    b, se, pv = BL.cox_breslow(osm[ok], osd[ok], iv)
    hr_i = float(np.exp(b[2]))
    t3 = pd.DataFrame([{"term": "stroma:TP53 interaction", "n": int(ok.sum()),
                        "n_event": int(osd[ok].sum()),
                        "HR_interaction": round(hr_i, 3),
                        "ci_lo": round(float(np.exp(b[2] - 1.96 * se[2])), 3),
                        "ci_hi": round(float(np.exp(b[2] + 1.96 * se[2])), 3),
                        "p": float(pv[2])}])
    t3.to_csv(OUT / "stage16_tp53_os.csv", index=False)
    print(f"\n[3] stroma × TP53 → OS 交互: HR={hr_i:.3f} P={pv[2]:.3g}")

    print(f"\n输出 -> {OUT}/stage16_*.csv")


if __name__ == "__main__":
    main()

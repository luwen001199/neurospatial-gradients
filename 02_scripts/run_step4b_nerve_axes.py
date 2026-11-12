# -*- coding: utf-8 -*-
"""Step4b —— METABRIC: 假设驱动的神经信号轴 + 基质校正后的 NPP。

背景: Step4a 显示空间推导的 NPP 打分在 bulk 层面对 RFS/OS 无显著预测力
(HR≈0.95, p=0.15)。一个显而易见的解释是: NPP 中富集了大量 ECM/CAF 基因,
在 bulk 组织里它的方差主要反映"通用促纤维增生程度", 而不是"神经近端性"。
因此这里做两件事:
  1) 报告 NPP 与通用 CAF/ECM 评分的相关性, 量化这种"被稀释"的程度;
  2) 用回归残差剥离通用基质成分, 检验残余的"神经特异"部分是否仍有预后信息;
  3) 同时检验方案预设的各条神经信号轴(CGRP / NGF-TrkA / GABA-BDNF / Semaphorin
     / 交感 / 副交感 / 感觉)在 bulk 中的预后方向, 这是假设驱动而非数据驱动的检验。
"""
from __future__ import annotations
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, r"D:/嗜神经浸润/00_pipeline")
from nsgrad import config as CFG          # noqa: E402
from nsgrad import metabric as mb         # noqa: E402
from nsgrad import bulk as nb             # noqa: E402
from nsgrad import surv as nsurv         # noqa: E402

PROC = CFG.DIR_PROC
DIR_RAW = CFG.DIR_RAW
TAB = CFG.DIR_TAB
TAB.mkdir(parents=True, exist_ok=True)

AXES = {
    "Sensory_nerve": ["NGF", "UCHL1", "TRPV1", "CALCA", "S100B", "SCN9A", "SCN10A"],
    "Sympathetic": ["TH", "DBH"],
    "Parasympathetic": ["SLC18A3", "CHAT"],
    "Structural_nerve": ["TUBB3", "MAP2", "NEFL", "NEFM", "PLP1", "MPZ", "SYP"],
    "CGRP_axis": ["CALCA", "CALCRL", "RAMP1", "CALCR"],
    "NGF_TrkA_axis": ["NGF", "NTRK1", "NGFR", "SORT1"],
    "GABA_BDNF_axis": ["GAD1", "GAD2", "ABAT", "SLC6A1", "GABBR1", "GABBR2",
                       "BDNF", "NTF3", "NTRK2"],
    "Semaphorin_axis": ["SEMA3A", "SEMA3C", "NRP1", "NRP2", "PLXNA1"],
    "CAF_ECM": ["COL1A1", "COL1A2", "COL3A1", "COL5A1", "COL11A1", "FN1", "POSTN",
                "ACTA2", "TAGLN", "FAP", "PDGFRB", "SPARC", "THBS2", "CTHRC1",
                "LOX", "INHBA"],
    "CD8_cytotox": ["CD8A", "CD8B", "GZMA", "GZMB", "PRF1", "NKG7", "IFNG"],
}


def z(s: pd.Series) -> pd.Series:
    return (s - s.mean()) / s.std(ddof=0)


def main():
    e1 = pd.read_csv(DIR_RAW / "metabric" / "metabric_mrna_targetgenes.csv", index_col=0)
    e2 = pd.read_csv(DIR_RAW / "metabric" / "metabric_mrna_npp.csv", index_col=0)
    expr = e1.join(e2, how="inner", rsuffix="_npp")
    print(f"[Step4b] 合并表达矩阵 {expr.shape}")

    import json
    sig = json.load(open(PROC / "npp_signature.json", encoding="utf-8"))
    npp = nb.signature_score(expr, sig["up"], sig["down"])

    S = pd.DataFrame({"NPPscore": npp})
    for k, genes in AXES.items():
        g = [x for x in genes if x in expr.columns]
        S[k] = nb.zscore_df(expr[g]).mean(axis=1) if g else np.nan
    S = S.dropna()
    print(f"[Step4b] 打分样本 n={len(S)}")

    # ---- 1) NPP 与通用基质评分的共线性 ----
    print("\n[Step4b] NPP 与各轴评分的 Spearman 相关:")
    corr_rows = []
    for c in S.columns:
        if c == "NPPscore":
            continue
        r = stats.spearmanr(S["NPPscore"], S[c])
        corr_rows.append(dict(pair="NPPscore", axis=c, rho=round(float(r.statistic), 3),
                              p=f"{r.pvalue:.2g}"))
        print(f"    NPPscore ~ {c:18s} rho={r.statistic:+.3f}  p={r.pvalue:.2g}")
    pd.DataFrame(corr_rows).to_csv(TAB / "step4b_npp_axis_correlation.csv", index=False)

    # ---- 2) 剥离通用基质后的残余 NPP ----
    X = S[["CAF_ECM"]].copy()
    X["CD8"] = S["CD8_cytotox"]
    Xc = (X - X.mean()) / X.std(ddof=0)
    Xc["intercept"] = 1.0
    beta, *_ = np.linalg.lstsq(Xc.values, z(S["NPPscore"]).values, rcond=None)
    resid = z(S["NPPscore"]).values - Xc.values @ beta
    S["NPP_resid_stroma"] = resid
    r2 = 1 - np.var(resid) / np.var(z(S["NPPscore"]).values)
    print(f"\n[Step4b] 通用 CAF/ECM + CD8 可解释 NPP 方差的 {r2*100:.1f}% "
          f"(残余方差 {100-r2*100:.1f}%)")

    # ---- 3) 合并临床 ----
    pat, smp = mb.load_clinical(DIR_RAW / "metabric")
    clin = mb.build_outcome_table(pat, smp)
    d = clin.join(S, how="inner")
    d = d[d["RFS_MONTHS"] > 0]
    # 注意: METABRIC 源数据里 ER_IHC 的阳性值拼写为 "POSITVE"(少一个 I), 两个拼写都要收
    d["ER_pos"] = (d.get("ER_IHC").astype(str).str.upper()
                   .isin(["POSITIVE", "POSITVE"])).astype(float)
    d["chemo"] = (d["CHEMOTHERAPY"].astype(str).str.upper() == "YES").astype(float)
    d["hormone"] = (d["HORMONE_THERAPY"].astype(str).str.upper() == "YES").astype(float)
    d["nodes_pos"] = d["NODES_POS"].fillna(0)
    print(f"[Step4b] 合并临床后 n={len(d)}")

    rows = []
    base_cov = ["AGE", "nodes_pos", "ER_pos", "chemo", "hormone"]
    for outcome, tcol, ecol in [("RFS", "RFS_MONTHS", "RFS_EVENT"),
                                ("OS", "OS_MONTHS", "OS_EVENT")]:
        sub = d.dropna(subset=[tcol, ecol])
        sub = sub[sub[tcol] > 0]
        print(f"\n===== {outcome} (n={len(sub)}, events={int(sub[ecol].sum())}) =====")
        out = {}
        for var in ["NPPscore", "NPP_resid_stroma"] + list(AXES.keys()):
            if var not in sub.columns:
                continue
            zz = z(sub[var]).values
            r = nsurv.cox(sub[tcol], sub[ecol], zz, names=[var])[0]
            out[var] = r
            # 多变量: 校正临床协变量 + 通用基质
            covs = [var] + [c for c in base_cov if c in sub.columns] + \
                   (["CAF_ECM"] if var != "CAF_ECM" else [])
            sub2 = sub.dropna(subset=covs)
            Zc = pd.DataFrame({c: z(sub2[c]) for c in covs})
            keep = [c for c in covs
                    if np.isfinite(Zc[c]).all() and Zc[c].std(ddof=0) > 1e-9]
            sub2 = sub2.loc[Zc.index]
            if len(keep) >= 2 and var in keep:
                r2m = nsurv.cox(sub2[tcol], sub2[ecol],
                                Zc[keep].values, names=keep)
                r2m = [r for r in r2m if r.name == var][0]
            else:
                r2m = nsurv.CoxResult(var, np.nan, np.nan, np.nan, np.nan,
                                      int(len(sub2)), int(sub2[ecol].sum()))
            rows.append(dict(cohort="METABRIC", outcome=outcome, model="univariate",
                             **vars(r)))
            rows.append(dict(cohort="METABRIC", outcome=outcome,
                             model="multivariable(+age,nodes,ER,chemo,hormone,CAF_ECM)",
                             **vars(r2m)))
            print(f"    {var:18s} uni HR={r.hr:.3f} ({r.lo:.2f}-{r.hi:.2f}) p={r.p:.3g} "
                  f"| adj HR={r2m.hr:.3f} ({r2m.lo:.2f}-{r2m.hi:.2f}) p={r2m.p:.3g}")

        # Basal 亚型内
        if "CLAUDIN_SUBTYPE" in sub.columns:
            for st in ["Basal", "LumA", "LumB", "Her2"]:
                ss = sub[sub["CLAUDIN_SUBTYPE"].astype(str) == st]
                if len(ss) < 40 or ss[ecol].sum() < 10:
                    continue
                for var in ["NPPscore", "NPP_resid_stroma", "CGRP_axis", "NGF_TrkA_axis"]:
                    if var not in ss.columns:
                        continue
                    r = nsurv.cox(ss[tcol], ss[ecol], z(ss[var]).values,
                                  names=[f"{var}[{st}]"])[0]
                    rows.append(dict(cohort="METABRIC", outcome=outcome,
                                     model=f"subtype_{st}", **vars(r)))
                    print(f"      [{st}] {var:16s} HR={r.hr:.3f} "
                          f"({r.lo:.2f}-{r.hi:.2f}) p={r.p:.3g}")

    df = pd.DataFrame(rows)
    df.to_csv(TAB / "step4b_metabric_nerve_axes.csv", index=False)
    print(f"\n  -> {TAB/'step4b_metabric_nerve_axes.csv'}")
    S.to_csv(PROC / "metabric_axis_scores.csv")
    print(f"  -> {PROC/'metabric_axis_scores.csv'}")

    # KM by NPP tertile (RFS) 数据
    d["NPP_tertile"] = pd.qcut(d["NPPscore"], 3, labels=["T1", "T2", "T3"])
    d.to_csv(PROC / "metabric_analysis_table.csv")
    return d, S


if __name__ == "__main__":
    main()

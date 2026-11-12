# -*- coding: utf-8 -*-
"""Step4a —— METABRIC 队列: 把空间推导的「神经近端程序」(NPP) 投影到 bulk,
检验其对 RFS / OS 的独立预后价值。

输出
----
01_data/processed/metabric_npp_score.csv
03_results/tables/step4_metabric_cox.csv
03_results/tables/step4_metabric_km.csv
"""
from __future__ import annotations
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, r"D:/嗜神经浸润/00_pipeline")
from nsgrad import config as CFG            # noqa: E402
from nsgrad import metabric as mb           # noqa: E402
from nsgrad import bulk as nb               # noqa: E402
from nsgrad import surv as nsurv           # noqa: E402

PROC = CFG.DIR_PROC
DIR_RAW = CFG.DIR_RAW
TAB = CFG.DIR_TAB
TAB.mkdir(parents=True, exist_ok=True)
SEED = 20251026
NPERM = 300


def load_signature():
    with open(PROC / "npp_signature.json", encoding="utf-8") as f:
        return json.load(f)


def build_expression(sig):
    genes = list(dict.fromkeys(sig["up"] + sig["down"]))
    print(f"[Step4a] NPP 签名基因 {len(genes)} 个 (up={len(sig['up'])}, down={len(sig['down'])})")
    cache = DIR_RAW / "metabric" / "metabric_mrna_npp.csv"
    df = nb.fetch_metabric_genes(genes, workers=8, cache=cache)
    print(f"    METABRIC 表达矩阵 {df.shape}")
    miss = [g for g in genes if g not in df.columns]
    if miss:
        print(f"    平台缺失 {len(miss)} 个基因: {miss[:10]}")
    return df


def main():
    sig = load_signature()
    expr = build_expression(sig)
    score = nb.signature_score(expr, sig["up"], sig["down"])
    score = score.dropna()

    pat, smp = mb.load_clinical(DIR_RAW / "metabric")
    clin = mb.build_outcome_table(pat, smp)
    d = clin.join(score, how="inner")
    print(f"\n[Step4a] 合并后 n={len(d)} (有表达+有临床)")

    d["NPP_z"] = (d["NPPscore"] - d["NPPscore"].mean()) / d["NPPscore"].std(ddof=0)
    d["NPP_tertile"] = pd.qcut(d["NPPscore"], 3, labels=["T1_low", "T2", "T3_high"])
    d["log2_size"] = np.log2(d["SIZE"].clip(lower=1))
    d["nodes_pos"] = d["NODES_POS"].fillna(0)
    # METABRIC 源数据 ER_IHC 阳性拼写为 "POSITVE", 两种写法都要收
    d["ER_pos"] = (d.get("ER_IHC").astype(str).str.upper()
                   .isin(["POSITIVE", "POSITVE"])).astype(float) \
        if "ER_IHC" in d else np.nan
    d["chemo"] = (d["CHEMOTHERAPY"].astype(str).str.upper() == "YES").astype(float) \
        if "CHEMOTHERAPY" in d else np.nan
    d["hormone"] = (d["HORMONE_THERAPY"].astype(str).str.upper() == "YES").astype(float) \
        if "HORMONE_THERAPY" in d else np.nan
    d["grade3"] = (d["GRADE"] == 3).astype(float)

    rows = []
    for outcome, tcol, ecol in [("RFS", "RFS_MONTHS", "RFS_EVENT"),
                                ("OS", "OS_MONTHS", "OS_EVENT")]:
        sub = d.dropna(subset=[tcol, ecol, "NPP_z"])
        sub = sub[sub[tcol] > 0]
        print(f"\n===== {outcome}: n={len(sub)} events={int(sub[ecol].sum())} =====")
        # 连续变量, 每 1 SD
        r = nsurv.cox(sub[tcol], sub[ecol], sub["NPP_z"].values, names=["NPPscore(per SD)"])
        nsurv.cox_table(r, f"{outcome} 单变量")
        rows += [dict(cohort="METABRIC", outcome=outcome, model="univariate",
                      **vars(r[0]))]
        # 多变量
        covs = ["NPP_z", "AGE", "log2_size", "grade3", "nodes_pos", "ER_pos",
                "chemo", "hormone"]
        covs = [c for c in covs if c in sub.columns and sub[c].notna().mean() > 0.5]
        sub2 = sub.dropna(subset=covs)
        # 必须剔除零方差列: 联合 Cox 里混入常数项会让整个模型求解失败
        Zc = pd.DataFrame({c: (sub2[c] - sub2[c].mean()) / sub2[c].std(ddof=0)
                           for c in covs})
        covs = [c for c in covs
                if np.isfinite(Zc[c]).all() and Zc[c].std(ddof=0) > 1e-9]
        r2 = nsurv.cox(sub2[tcol], sub2[ecol], Zc[covs].values, names=covs)
        nsurv.cox_table(r2, f"{outcome} 多变量 (n={len(sub2)})")
        for rr in r2:
            rows.append(dict(cohort="METABRIC", outcome=outcome, model="multivariable",
                             **vars(rr)))
        # 三分位 KM
        g = sub.groupby("NPP_tertile", observed=True)
        for lev, gg in g:
            pass
        lo = sub[sub["NPP_tertile"] == "T1_low"]; hi = sub[sub["NPP_tertile"] == "T3_high"]
        chi, p = nsurv.logrank(lo[tcol], lo[ecol], hi[tcol], hi[ecol])
        print(f"    {outcome} T3(高NPP) vs T1(低NPP) log-rank chi2={chi:.2f} p={p:.3g}")
        rows.append(dict(cohort="METABRIC", outcome=outcome, model="tertile_T3_vs_T1",
                         name="logrank", hr=np.nan, lo=np.nan, hi=np.nan, p=float(p),
                         n=int(len(lo) + len(hi)), events=int(lo[ecol].sum() + hi[ecol].sum()),
                         coef=np.nan, se=np.nan))

        # 亚型分层
        if "CLAUDIN_SUBTYPE" in sub.columns:
            for st in ["Basal", "LumA", "LumB", "Her2"]:
                ss = sub[sub["CLAUDIN_SUBTYPE"].astype(str) == st]
                if len(ss) < 40 or ss[ecol].sum() < 10:
                    continue
                rs = nsurv.cox(ss[tcol], ss[ecol], ss["NPP_z"].values,
                               names=[f"NPPscore[{st}]"])
                print(f"    [{st}] n={len(ss)} HR={rs[0].hr:.3f} "
                      f"({rs[0].lo:.2f}-{rs[0].hi:.2f}) p={rs[0].p:.3g}")
                rows.append(dict(cohort="METABRIC", outcome=outcome,
                                 model=f"subtype_{st}", **vars(rs[0])))

    # ---- 置换检验: 随机等大基因集能否得到同等强度的 HR ----
    print(f"\n[Step4a] 置换检验 ({NPERM} 次随机等大基因集) ...")
    rng = np.random.default_rng(SEED)
    all_genes = list(expr.columns)
    n_up, n_dn = len(sig["up"]), len(sig["down"])
    obs = {}
    for outcome, tcol, ecol in [("RFS", "RFS_MONTHS", "RFS_EVENT")]:
        sub = d.dropna(subset=[tcol, ecol]); sub = sub[sub[tcol] > 0]
        r = nsurv.cox(sub[tcol], sub[ecol], sub["NPP_z"].values, names=["x"])[0]
        obs[outcome] = r.hr
    null = {k: [] for k in obs}
    for i in range(NPERM):
        perm_up = list(rng.choice(all_genes, size=n_up, replace=False))
        rest = [g for g in all_genes if g not in set(perm_up)]
        perm_dn = list(rng.choice(rest, size=min(n_dn, len(rest)), replace=False)) if n_dn else []
        s = nb.signature_score(expr, perm_up, perm_dn).dropna()
        dd = d.join(s.rename("pz"), how="inner")
        dd["pz"] = (dd["pz"] - dd["pz"].mean()) / dd["pz"].std(ddof=0)
        for outcome, tcol, ecol in [("RFS", "RFS_MONTHS", "RFS_EVENT")]:
            sub = dd.dropna(subset=[tcol, ecol, "pz"]); sub = sub[sub[tcol] > 0]
            try:
                rr = nsurv.cox(sub[tcol], sub[ecol], sub["pz"].values, names=["x"])[0]
                null[outcome].append(rr.hr)
            except Exception:
                pass
    perm_rows = []
    for outcome, hrs in null.items():
        hrs = np.array([h for h in hrs if np.isfinite(h)])
        if len(hrs) < 20:
            continue
        o = obs[outcome]
        p_emp = float((np.abs(np.log(hrs)) >= abs(np.log(o))).mean())
        print(f"    {outcome}: 观测 HR={o:.3f}, 置换中位数 HR={np.median(hrs):.3f}, "
              f"经验 p={p_emp:.3f}")
        perm_rows.append(dict(outcome=outcome, observed_HR=round(o, 3),
                              null_median_HR=round(float(np.median(hrs)), 3),
                              null_2_5=round(float(np.percentile(hrs, 2.5)), 3),
                              null_97_5=round(float(np.percentile(hrs, 97.5)), 3),
                              n_perm=len(hrs), empirical_p=round(p_emp, 4)))

    out = pd.DataFrame(rows)
    out.to_csv(TAB / "step4_metabric_cox.csv", index=False)
    print(f"\n  -> {TAB/'step4_metabric_cox.csv'}")
    if perm_rows:
        pd.DataFrame(perm_rows).to_csv(TAB / "step4_metabric_permutation.csv", index=False)
        print(f"  -> {TAB/'step4_metabric_permutation.csv'}")

    d[["NPPscore", "NPP_z", "NPP_tertile", "RFS_MONTHS", "RFS_EVENT",
       "OS_MONTHS", "OS_EVENT", "AGE", "SIZE", "GRADE", "NODES_POS",
       "ER_IHC", "CHEMOTHERAPY", "CLAUDIN_SUBTYPE"]].to_csv(
        PROC / "metabric_npp_score.csv")
    print(f"  -> {PROC/'metabric_npp_score.csv'}")

    # KM 曲线数据
    km_rows = []
    for outcome, tcol, ecol in [("RFS", "RFS_MONTHS", "RFS_EVENT"),
                                ("OS", "OS_MONTHS", "OS_EVENT")]:
        sub = d.dropna(subset=[tcol, ecol]); sub = sub[sub[tcol] > 0]
        for lev in ["T1_low", "T2", "T3_high"]:
            gg = sub[sub["NPP_tertile"] == lev]
            tt, ss, nr, ne = nsurv.km(gg[tcol], gg[ecol])
            for a, b, c, e2 in zip(tt, ss, nr, ne):
                km_rows.append(dict(outcome=outcome, group=lev, t=a, S=b, n_risk=c, n_event=e2))
    pd.DataFrame(km_rows).to_csv(TAB / "step4_metabric_km.csv", index=False)
    print(f"  -> {TAB/'step4_metabric_km.csv'}")
    return d


if __name__ == "__main__":
    main()

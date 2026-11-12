# -*- coding: utf-8 -*-
"""NPP 签名的空间自相关稳健化修订 + METABRIC 复现。

目的
----
MSR 重算后, 若原 NPP 基因集大量无法在"保留空间自相关"的零分布下存活,
则该签名可能是空间光滑性的副产物。这里直接用 MSR 存活基因重建签名,
并在 METABRIC 中复现预后关联 —— 判断"空间层面的假阳性"是否传导到
bulk 层面的临床结论。

候选签名
--------
  A  NPP_original      原宽松集 (秩变换 bug 修复后)
  B  NPP_msr_fdr       MSR 零分布 FDR<0.05 且 |mean rho| >= 0.05
  C  NPP_msr_top50     MSR 零分布下 Stouffer Z 最极端的 up/down 各 50
                       (当 B 为空或过小时的替代方案)

输出
----
03_results/tables/npp_revision_metabric.csv
"""
from __future__ import annotations
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, r"D:/嗜神经浸润/00_pipeline")
from nsgrad import config as CFG  # noqa: E402
from nsgrad import metabric as mb  # noqa: E402
from nsgrad import bulk as nb  # noqa: E402
from nsgrad import surv as nsurv  # noqa: E402

PROC = CFG.DIR_PROC
TAB = CFG.DIR_TAB
RAW = CFG.DIR_RAW
TAB.mkdir(parents=True, exist_ok=True)

CAF = ["COL1A1", "COL1A2", "COL3A1", "COL5A1", "COL11A1", "FN1", "POSTN",
       "ACTA2", "TAGLN", "FAP", "PDGFRB", "SPARC", "THBS2", "CTHRC1",
       "LOX", "INHBA"]
CD8 = ["CD8A", "CD8B", "GZMA", "GZMB", "PRF1", "NKG7", "IFNG"]


def z(s: pd.Series) -> pd.Series:
    return (s - s.mean()) / s.std(ddof=0)


def load_expr():
    e1 = pd.read_csv(RAW / "metabric" / "metabric_mrna_targetgenes.csv", index_col=0)
    e2 = pd.read_csv(RAW / "metabric" / "metabric_mrna_npp.csv", index_col=0)
    expr = e1.join(e2, how="inner", rsuffix="_npp")
    print(f"[rev] 表达矩阵 {expr.shape}")
    return expr


def build_candidates(G: pd.DataFrame, sig_json: dict):
    up_o = [g for g in sig_json.get("up", [])]
    dn_o = [g for g in sig_json.get("down", [])]

    cands = {
        "NPP_original_fixed_rank": (up_o, dn_o),
    }
    # 秩变换 bug 修复前的签名(仅用于复现旧行为, 验证信号来源)
    pre = PROC / "npp_signature_PRE_RANKFIX.json"
    if pre.exists():
        old = json.load(open(pre, encoding="utf-8"))
        cands["NPP_pre_rankfix"] = (list(old.get("up", [])), list(old.get("down", [])))

    eff = G["mean_rho"].abs() >= 0.05
    m = G[(G["fdr_msr"] < 0.05) & eff]
    cands["NPP_msr_fdr"] = (list(m.index[m["mean_rho"] < 0]),
                            list(m.index[m["mean_rho"] > 0]))

    up_c = list(G[(G["mean_rho"] < 0)].sort_values("stouffer_Z").head(50).index)
    dn_c = list(G[(G["mean_rho"] > 0)].sort_values("stouffer_Z", ascending=False).head(50).index)
    cands["NPP_msr_top50"] = (up_c, dn_c)
    return cands


def evaluate(expr: pd.DataFrame, name: str, up: list, dn: list):
    up2 = [g for g in up if g in expr.columns]
    dn2 = [g for g in dn if g in expr.columns]
    if len(up2) < 5:
        print(f"[rev] {name}: 矩阵中仅剩 {len(up2)} 个 up 基因, 跳过")
        return None
    score = nb.signature_score(expr, up2, dn2)

    S = pd.DataFrame({"score": score})
    S["CAF_ECM"] = nb.zscore_df(expr[[g for g in CAF if g in expr.columns]]).mean(axis=1)
    S["CD8"] = nb.zscore_df(expr[[g for g in CD8 if g in expr.columns]]).mean(axis=1)
    S = S.dropna()
    rho_caf = float(stats.spearmanr(S["score"], S["CAF_ECM"]).statistic)

    Xc = pd.DataFrame({"CAF": z(S["CAF_ECM"]), "CD8": z(S["CD8"]), "int": 1.0})
    beta, *_ = np.linalg.lstsq(Xc.values, z(S["score"]).values, rcond=None)
    resid = z(S["score"]).values - Xc.values @ beta
    S["score_resid"] = resid
    r2 = float(1 - np.var(resid) / np.var(z(S["score"]).values))

    pat, smp = mb.load_clinical(RAW / "metabric")
    clin = mb.build_outcome_table(pat, smp)
    d = clin.join(S, how="inner")
    d = d[d["RFS_MONTHS"] > 0]
    d["ER_pos"] = (d.get("ER_IHC").astype(str).str.upper()
                   .isin(["POSITIVE", "POSITVE"])).astype(float)
    d["chemo"] = (d["CHEMOTHERAPY"].astype(str).str.upper() == "YES").astype(float)
    d["hormone"] = (d["HORMONE_THERAPY"].astype(str).str.upper() == "YES").astype(float)
    d["nodes_pos"] = d["NODES_POS"].fillna(0)

    covs = ["AGE", "nodes_pos", "ER_pos", "chemo", "hormone"]
    out = []
    for outcome, tcol, ecol in [("RFS", "RFS_MONTHS", "RFS_EVENT"),
                                ("OS", "OS_MONTHS", "OS_EVENT")]:
        if tcol not in d.columns:
            continue
        for var in ["score", "score_resid"]:
            sub = d.dropna(subset=[var, tcol, ecol] + covs)
            sub = sub[sub[tcol] > 0]
            if len(sub) < 100 or sub[ecol].sum() < 20:
                continue
            Z = pd.DataFrame({c: (sub[c] - sub[c].mean()) / sub[c].std(ddof=0)
                              for c in covs})
            Z[var] = (sub[var] - sub[var].mean()) / sub[var].std(ddof=0)
            names = covs + [var]
            try:
                res = nsurv.cox(sub[tcol].to_numpy(float),
                                sub[ecol].to_numpy(float),
                                Z[names].to_numpy(float), names=names)
            except Exception as e:  # noqa: BLE001
                print(f"   cox failed {outcome}/{var}: {e}")
                continue
            row = next((r for r in res if r.name == var), None)
            if row is None:
                continue
            out.append(dict(signature=name, n_genes_up=len(up2), n_genes_dn=len(dn2),
                            outcome=outcome, variable=var,
                            n=int(len(sub)), events=int(sub[ecol].sum()),
                            HR=round(float(row.hr), 4),
                            CI95_low=round(float(row.lo), 4),
                            CI95_high=round(float(row.hi), 4),
                            p=float(row.p),
                            rho_CAF_ECM=round(rho_caf, 3),
                            var_explained_by_stroma=round(r2, 3)))
    return out


def main():
    G = pd.read_csv(PROC / "msr_gene_level.csv", index_col=0)
    sig = json.load(open(PROC / "npp_signature.json", encoding="utf-8"))
    expr = load_expr()
    cands = build_candidates(G, sig)
    rows = []
    for name, (u, d) in cands.items():
        print(f"\n[rev] === {name}: up={len(u)} down={len(d)} ===")
        r = evaluate(expr, name, u, d)
        if r:
            rows.extend(r)
            for x in r:
                print(f"   {x['outcome']} {x['variable']:12s} HR={x['HR']:.3f} "
                      f"({x['CI95_low']:.3f}-{x['CI95_high']:.3f}) p={x['p']:.4g} "
                      f"n={x['n']} events={x['events']}")
    if rows:
        out = pd.DataFrame(rows)
        out.to_csv(TAB / "npp_revision_metabric.csv", index=False)
        print(f"\n  -> {TAB/'npp_revision_metabric.csv'}")


if __name__ == "__main__":
    main()

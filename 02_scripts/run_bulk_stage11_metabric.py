# -*- coding: utf-8 -*-
"""Stage 11 —— METABRIC 的交互项与「签名构建不稳定性」定量。

原稿 3.5 写了「basal-like shows no residual association」，但没有:
  (1) 正式检验 stroma module × CAF/ECM 的交互（模块效应是否随纤维化水平改变）；
  (2) 把「签名构建不稳定性」量化 —— 即换一批基因重建签名，ρ 会飘多少。
本脚本用「基因子集重采样」给出 ρ 的经验分布，作为 construction instability 的度量。

输出:
    03_results/tables/bulk_extrapolation/stage11_metabric_interaction.csv
    03_results/tables/bulk_extrapolation/stage11_construction_bootstrap.csv
"""
from __future__ import annotations

import json
import sys

import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, r"D:/嗜神经浸润/02_scripts")
import bulk_lib as BL  # noqa: E402

OUT = BL.TAB
META = BL.RAW / "metabric"
RNG = np.random.default_rng(20261108)
N_BOOT = 1000


def load_metabric():
    expr = pd.read_csv(META / "metabric_mrna_targetgenes.csv", index_col=0)
    expr.index = expr.index.astype(str)
    raw = json.loads((META / "clinical_patient.txt").read_text(encoding="utf-8"))
    rec = {}
    for r in raw:
        rec.setdefault(r["patientId"], {})[r["clinicalAttributeId"]] = r.get("value")
    clin = pd.DataFrame.from_dict(rec, orient="index")
    clin.index = clin.index.astype(str)
    out = pd.DataFrame({
        "OS_MONTHS": pd.to_numeric(clin.get("OS_MONTHS"), errors="coerce"),
        "OS_STATUS": clin.get("OS_STATUS"),
        "AGE": pd.to_numeric(clin.get("AGE_AT_DIAGNOSIS"), errors="coerce"),
        "SUBTYPE": clin.get("CLAUDIN_SUBTYPE"),
        "NODE": clin.get("LYMPH_NODES_EXAMINED_POSITIVE"),
        "GRADE": clin.get("GRADE"),
    })
    out = out.reindex(expr.index)
    out["OS_EVENT"] = out["OS_STATUS"].astype(str).str.startswith("1").astype(float)
    out.loc[out["OS_STATUS"].isna(), "OS_EVENT"] = np.nan
    return expr.astype(float), out


def interaction_cox(sc: pd.DataFrame, clin: pd.DataFrame, mod: str, cov: str):
    """Cox: h ~ mod_z + cov_z + mod_z:cov_z + age；报告交互项 HR 与 p。"""
    d = pd.concat([sc[[mod, cov]], clin[["OS_MONTHS", "OS_EVENT", "AGE"]]], axis=1).dropna()
    d = d[d["OS_MONTHS"] > 0]
    if len(d) < 100:
        return None
    z = lambda s: (s - s.mean()) / s.std(ddof=0)
    x1, x2 = z(d[mod]).to_numpy(), z(d[cov]).to_numpy()
    age = z(d["AGE"].fillna(d["AGE"].median())).to_numpy()
    X = np.column_stack([x1, x2, x1 * x2, age])
    b, se, p = BL.cox_breslow(d["OS_MONTHS"].to_numpy(), d["OS_EVENT"].to_numpy(), X)
    return {"n": len(d), "events": int(d["OS_EVENT"].sum()),
            "HR_main": float(np.exp(b[0])), "p_main": float(p[0]),
            "HR_cov": float(np.exp(b[1])), "p_cov": float(p[1]),
            "HR_interaction": float(np.exp(b[2])), "p_interaction": float(p[2])}


def construction_bootstrap(expr: pd.DataFrame, mods: dict, mod: str, cov: str,
                           frac: float = 0.8, n_boot: int = N_BOOT):
    """随机抽 frac 比例的模块基因重建模块，看 ρ(module, cov) 的漂移。"""
    genes = [g for g in mods[mod] if g in expr.columns]
    if len(genes) < 10:
        return None, np.nan
    caf = BL.score_mean(expr, mods[cov])
    full = BL.score_mean(expr, genes)
    rho_full = stats.spearmanr(full, caf, nan_policy="omit")[0]
    k = max(8, int(round(len(genes) * frac)))
    rhos = np.empty(n_boot)
    for i in range(n_boot):
        sub = RNG.choice(genes, k, replace=False)
        s = BL.score_mean(expr, list(sub))
        rhos[i] = stats.spearmanr(s, caf, nan_policy="omit")[0]
    # 随机基因集基线（同大小、非模块内基因）
    pool = [g for g in expr.columns if g not in set(genes)]
    rand = np.empty(n_boot)
    for i in range(n_boot):
        s = BL.score_mean(expr, list(RNG.choice(pool, k, replace=False)))
        rand[i] = stats.spearmanr(s, caf, nan_policy="omit")[0]
    return {"module": mod, "n_genes": len(genes), "subsample_size": k,
            "rho_full": round(float(rho_full), 4),
            "boot_median": round(float(np.median(rhos)), 4),
            "boot_p2.5": round(float(np.percentile(rhos, 2.5)), 4),
            "boot_p97.5": round(float(np.percentile(rhos, 97.5)), 4),
            "boot_range": round(float(rhos.max() - rhos.min()), 4),
            "random_median": round(float(np.median(rand)), 4),
            "random_p97.5": round(float(np.percentile(rand, 97.5)), 4)}, rhos


def main():
    mods = BL.load_modules()
    expr, clin = load_metabric()
    print(f"METABRIC n={len(expr)}; 基因 {expr.shape[1]}; "
          f"OS 可算 {clin['OS_MONTHS'].notna().sum()}, 事件 {int(clin['OS_EVENT'].sum())}")
    miss = [g for g in mods["Perineural_stroma"] if g not in expr.columns]
    print(f"stroma 模块缺 {len(miss)} 基因: {miss[:8]}")

    sc = BL.score_all(expr, mods)
    sc["CD8_cytotox"] = BL.score_mean(expr, mods["CD8_cytotox"])
    resid, r2 = BL.rank_resid(sc["Perineural_stroma"], sc["CAF_ECM"])
    sc["Stroma_resid"] = resid
    print(f"模块被 CAF/ECM 解释 {r2*100:.1f}% 秩方差")

    rows = []
    for mod in ["Perineural_stroma", "Stroma_resid", "NPP58", "Neural_intrinsic",
                "SPARC_POSTN_CORE", "COLLAGEN_MATURATION"]:
        if mod not in sc:
            continue
        for cov in ["CAF_ECM"]:
            r = interaction_cox(sc, clin, mod, cov)
            if r:
                r.update({"module": mod, "covariate": cov})
                rows.append(r)
    it = pd.DataFrame(rows)
    print("\n交互项（Cox，per SD，校正 age）")
    print(it[["module", "n", "events", "HR_main", "p_main", "HR_interaction",
              "p_interaction"]].to_string(index=False))
    it.to_csv(OUT / "stage11_metabric_interaction.csv", index=False)

    # 亚型分层（原稿 3.5 的说法需要数字支撑）
    print("\n亚型分层：模块残差 ~ OS（per SD）")
    sub_rows = []
    for st in sorted(clin["SUBTYPE"].dropna().unique()):
        mask = clin["SUBTYPE"] == st
        d = pd.concat([sc.loc[mask, "Stroma_resid"], clin.loc[mask, ["OS_MONTHS", "OS_EVENT"]]],
                      axis=1).dropna()
        d = d[d["OS_MONTHS"] > 0]
        if len(d) < 60 or d["OS_EVENT"].sum() < 15:
            continue
        z = (d["Stroma_resid"] - d["Stroma_resid"].mean()) / d["Stroma_resid"].std(ddof=0)
        hr, lo, hi, p = BL.cox_hr(d["OS_MONTHS"].to_numpy(), d["OS_EVENT"].to_numpy(), z.to_numpy())
        sub_rows.append({"subtype": st, "n": len(d), "events": int(d["OS_EVENT"].sum()),
                         "HR": round(hr, 3), "lo": round(lo, 3), "hi": round(hi, 3),
                         "p": p})
        print(f"  {st:16s} n={len(d):4d} ev={int(d['OS_EVENT'].sum()):3d} "
              f"HR={hr:.2f} [{lo:.2f},{hi:.2f}] p={p:.3g}")
    pd.DataFrame(sub_rows).to_csv(OUT / "stage11_metabric_subtype.csv", index=False)

    # 构建不稳定性
    print("\n签名构建 bootstrap（随机抽 80% 模块基因重建，1000 次）")
    brow, dists = [], {}
    for mod in ["Perineural_stroma", "Neural_intrinsic"]:
        r, dist = construction_bootstrap(expr, mods, mod, "CAF_ECM")
        if r:
            brow.append(r)
            dists[mod] = dist
            print(f"  {mod:20s} ρ_full={r['rho_full']:+.3f}  "
                  f"boot {r['boot_median']:+.3f} [{r['boot_p2.5']:+.3f},{r['boot_p97.5']:+.3f}]  "
                  f"随机集中位 {r['random_median']:+.3f}")
    cb = pd.DataFrame(brow)
    cb.to_csv(OUT / "stage11_construction_bootstrap.csv", index=False)
    if dists:
        pd.DataFrame(dists).to_csv(OUT / "stage11_construction_dist.csv", index=False)
    print(f"\n输出 -> {OUT}")


if __name__ == "__main__":
    main()

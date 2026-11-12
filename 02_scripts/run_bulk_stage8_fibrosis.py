# -*- coding: utf-8 -*-
"""Stage 8 —— 纤维化的「量」与「质」：哪一层携带预后/耐药信息。

Stage 7 显示 CAF/ECM 总量（"量"）对 OS 与 pCR 均无预测力，而扣除总量后剩下的
SPARC/POSTN 程序（"质"）同时预测更差 OS 与更低 pCR。本脚本把这条线验证到底：

    1. 递增 Cox（M0 单变量 → M1 +age/node → M2 +ER/Basal → M3 +CAF/ECM 总量）
       看 SPARC/POSTN 是否独立于临床分期、分子亚型与纤维化总量
    2. 三分位 KM + log-rank
    3. 单基因分解（SPARC / POSTN / THBS2 / TNC / CTHRC1 各自残差化）
    4. collagen 成熟程序 vs 胶原核心程序的对照（"质" vs "量"）

输出: stage8_cox_models.csv / stage8_km.csv / stage8_gene_decomposition.csv
"""
from __future__ import annotations

import sys
import warnings

import numpy as np
import pandas as pd

sys.path.insert(0, r"D:/嗜神经浸润/02_scripts")
import bulk_lib as BL  # noqa: E402
from run_bulk_stage7_orthogonal import build_cohorts, score_table, zs  # noqa: E402

warnings.filterwarnings("ignore")
OUT = BL.TAB


def num(s) -> np.ndarray:
    return pd.to_numeric(s, errors="coerce").to_numpy(float)


def km_curve(t, e, g, n_grid=60):
    """Kaplan-Meier 阶梯曲线（返回 times, surv）。"""
    t = np.asarray(t, float)
    e = np.asarray(e, float)
    g = np.asarray(g, bool)
    out = {}
    for lab, m in (("low", g), ("high", ~g)):
        tt, ee = t[m], e[m]
        times = np.sort(np.unique(tt[ee == 1]))
        s, n = 1.0, len(tt)
        xs, ys = [0.0], [1.0]
        for x in times:
            at = (tt >= x).sum()
            dd = ((tt == x) & (ee == 1)).sum()
            if at > 0:
                s *= 1 - dd / at
            xs.append(x)
            ys.append(s)
        out[lab] = (np.array(xs), np.array(ys), n, int(ee.sum()))
    return out


def main():
    mods = BL.load_modules()
    cohorts = build_cohorts()
    rows, km_rows, gene_rows = [], [], []

    for name in ["TCGA-BRCA", "TCGA-TNBC"]:
        if name not in cohorts:
            continue
        expr, meta = cohorts[name]
        sc = score_table(expr, mods)

        t = num(meta["os_time"])
        d = num(meta["os_event"])
        ok = np.isfinite(t) & np.isfinite(d) & (t > 0)
        if ok.sum() < 60 or np.nansum(d[ok]) < 15:
            print(f"[skip surv] {name}: n={int(ok.sum())}, event={int(np.nansum(d[ok]))}")
            continue

        tv, dv = t[ok], d[ok].astype(int)
        age = num(meta["age"])[ok]
        age = np.where(np.isnan(age), np.nanmean(age), age)
        node = np.isin(meta["node"].astype(str).to_numpy()[ok], ["N1", "N2", "N3"]).astype(float)
        er = np.where(meta.get("pam50") is not None,
                      meta["pam50"].astype(str).to_numpy()[ok] == "Basal", False).astype(float)

        print(f"\n[{name}] 生存 n={len(tv)}, events={int(dv.sum())}")

        targets = ["SPARC_POSTN_CORE", "Perineural_stroma", "CAF_ECM", "COLLAGEN_CORE",
                   "COLLAGEN_MATURATION", "MECHANOTRANSDUCTION", "EMT_CORE", "TGFb_PROGRAM"]
        for m in targets:
            if m not in sc:
                continue
            sv = sc[m].to_numpy()[ok]
            if not np.isfinite(sv).all():
                continue
            z = (sv - sv.mean()) / sv.std(ddof=0)
            # raw
            for mname, cov, tag in [
                ("M0_raw", None, "raw"),
                ("M1_raw", np.column_stack([age, node]), "raw"),
                ("M2_raw", np.column_stack([age, node, er]), "raw"),
            ]:
                hr, lo, hi, p = BL.cox_hr(tv, dv, z, cov)
                rows.append({"cohort": name, "module": m, "scale": tag, "model": mname,
                             "n": len(tv), "n_event": int(dv.sum()),
                             "HR": round(hr, 3), "lo": round(lo, 3), "hi": round(hi, 3), "p": p})
            # resid（对 CAF/ECM 正交化）
            r, r2 = BL.rank_resid(sc[m], sc["CAF_ECM"])
            rv = r.to_numpy()[ok]
            if not np.isfinite(rv).all():
                continue
            zr = (rv - rv.mean()) / rv.std(ddof=0)
            for mname, cov in [
                ("M0_resid", None),
                ("M1_resid", np.column_stack([age, node])),
                ("M2_resid", np.column_stack([age, node, er])),
                ("M3_resid", np.column_stack([age, node, er, (sc["CAF_ECM"].to_numpy()[ok] -
                                                              sc["CAF_ECM"].to_numpy()[ok].mean()) /
                                              sc["CAF_ECM"].to_numpy()[ok].std(ddof=0)])),
            ]:
                hr, lo, hi, p = BL.cox_hr(tv, dv, zr, cov)
                rows.append({"cohort": name, "module": m, "scale": "resid", "model": mname,
                             "n": len(tv), "n_event": int(dv.sum()),
                             "HR": round(hr, 3), "lo": round(lo, 3), "hi": round(hi, 3), "p": p})

        # ---- KM（SPARC_POSTN_CORE 残差三分位）
        if "SPARC_POSTN_CORE" in sc:
            r, _ = BL.rank_resid(sc["SPARC_POSTN_CORE"], sc["CAF_ECM"])
            rv = pd.Series(r.to_numpy()[ok], index=np.arange(ok.sum()))
            q = pd.qcut(rv, 3, labels=["T1", "T2", "T3"], duplicates="drop")
            hi = (q == "T3").to_numpy()
            lo = (q == "T1").to_numpy()
            chi, p = BL.logrank(tv[hi | lo], dv[hi | lo], hi[hi | lo].astype(int))
            print(f"  KM: SPARC_POSTN resid T3 vs T1 log-rank chi2={chi:.2f}, p={p:.4g}")
            km = km_curve(tv[hi | lo], dv[hi | lo], hi[hi | lo])
            for lab, (xs, ys, n_, ev) in km.items():
                for x, y in zip(xs, ys):
                    km_rows.append({"cohort": name, "group": lab, "time_days": x,
                                    "surv": y, "n": n_, "events": ev,
                                    "logrank_chi2": chi, "logrank_p": p})
            # 5 年生存率
            for lab, (xs, ys, n_, ev) in km.items():
                s5 = ys[np.searchsorted(xs, 1825, side="right") - 1]
                rows.append({"cohort": name, "module": "SPARC_POSTN_CORE_KM", "scale": f"resid_{lab}",
                             "model": "5y_survival", "n": n_, "n_event": ev,
                             "HR": round(float(s5), 3), "lo": np.nan, "hi": np.nan, "p": p})

        # ---- 单基因分解
        for g in ["SPARC", "POSTN", "THBS2", "TNC", "CTHRC1", "COL1A1", "LRRC15", "LOX", "PLOD2"]:
            if g not in expr.columns:
                continue
            v = expr[g].to_numpy(float)[ok]
            if not np.isfinite(v).all() or np.std(v) == 0:
                continue
            zg = (v - v.mean()) / v.std(ddof=0)
            hr, lo, hi, p = BL.cox_hr(tv, dv, zg)
            # 残差化（对 CAF_ECM）
            rr, _ = BL.rank_resid(expr[g], sc["CAF_ECM"])
            vv = rr.to_numpy()[ok]
            zr = (vv - vv.mean()) / vv.std(ddof=0)
            hr2, lo2, hi2, p2 = BL.cox_hr(tv, dv, zr, np.column_stack([age, node]))
            gene_rows.append({"cohort": name, "gene": g, "HR_raw": round(hr, 3),
                              "lo_raw": round(lo, 3), "hi_raw": round(hi, 3), "p_raw": p,
                              "HR_resid_adj": round(hr2, 3), "lo_resid": round(lo2, 3),
                              "hi_resid": round(hi2, 3), "p_resid_adj": p2})

    pd.DataFrame(rows).to_csv(OUT / "stage8_cox_models.csv", index=False)
    if km_rows:
        pd.DataFrame(km_rows).to_csv(OUT / "stage8_km.csv", index=False)
    if gene_rows:
        pd.DataFrame(gene_rows).to_csv(OUT / "stage8_gene_decomposition.csv", index=False)

    # ---- 摘要
    d = pd.DataFrame(rows)
    if len(d):
        print("\n" + "=" * 78)
        print("递增 Cox：SPARC/POSTN 是否独立于亚型与纤维化总量（HR per SD）")
        sub = d[d["module"].isin(["SPARC_POSTN_CORE", "CAF_ECM", "COLLAGEN_CORE",
                                  "COLLAGEN_MATURATION"])]
        print(sub[sub["model"] != "5y_survival"]
              [["cohort", "module", "scale", "model", "HR", "lo", "hi", "p"]]
              .to_string(index=False))
    if gene_rows:
        print("\n" + "=" * 78)
        print("单基因分解（TCGA-BRCA）")
        print(pd.DataFrame(gene_rows)[["cohort", "gene", "HR_raw", "p_raw",
                                       "HR_resid_adj", "p_resid_adj"]].to_string(index=False))
    print(f"\n输出 -> {OUT}")


if __name__ == "__main__":
    main()

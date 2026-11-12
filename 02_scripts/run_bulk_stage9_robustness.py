# -*- coding: utf-8 -*-
"""Stage 9 —— SPARC/POSTN 残差预后信号的稳健性检验。

Stage 8 的关键结果：扣除 CAF/ECM 总量后，SPARC_POSTN_CORE 残差预测更差 OS
（TCGA-BRCA，HR 1.30/SD，p=0.007，校正 age/node/亚型/总量后不变），而同批
7 个其它基质程序全部为 null。为避免"又一个假签名"，本脚本做四件事：

    1. 随机基因集基线：510 组 size-matched（8 基因）随机集，同法残差化后 Cox，
       给出观测 HR 的经验 P（与 Stage 6 的空间随机基线同一逻辑）
    2. Bootstrap 1000× 的 HR 95% CI（队列内重抽样）
    3. 加增殖（MKI67）校正——低增殖本身预后差，需排除
    4. 改进的分组 KM（上三分位 vs 其余）+ log-rank + 5 年生存差

输出: stage9_random_baseline.csv / stage9_robustness.csv / stage9_km_refined.csv
"""
from __future__ import annotations

import sys
import warnings

import numpy as np
import pandas as pd

sys.path.insert(0, r"D:/嗜神经浸润/02_scripts")
import bulk_lib as BL  # noqa: E402
from run_bulk_stage7_orthogonal import build_cohorts, score_table  # noqa: E402

warnings.filterwarnings("ignore")
OUT = BL.TAB
SEED = 20251108
N_RANDOM = 510


def surv_arrays(meta, ok_mask):
    t = pd.to_numeric(meta["os_time"], errors="coerce").to_numpy(float)
    d = pd.to_numeric(meta["os_event"], errors="coerce").to_numpy(float)
    ok = np.isfinite(t) & np.isfinite(d) & (t > 0)
    return t[ok], d[ok].astype(int), ok


def hr_of(score: pd.Series, caf: pd.Series, t, d, ok, cov=None):
    r, _ = BL.rank_resid(score, caf)
    rv = r.to_numpy(float)
    if not np.isfinite(rv).all():
        return np.nan
    rv = rv[ok]
    z = (rv - rv.mean()) / rv.std(ddof=0)
    hr, lo, hi, p = BL.cox_hr(t, d, z, cov)
    return hr, lo, hi, p


def km_groups(t, d, g, labels):
    """多组 KM。返回 [(label, times, surv, n, events)]。"""
    res = []
    for lab in labels:
        m = g == lab
        tt, ee = t[m], d[m]
        times = np.sort(np.unique(tt[ee == 1]))
        s = 1.0
        xs, ys = [0.0], [1.0]
        for x in times:
            at = (tt >= x).sum()
            dd = ((tt == x) & (ee == 1)).sum()
            if at > 0:
                s *= 1 - dd / at
            xs.append(x)
            ys.append(s)
        res.append((lab, np.array(xs), np.array(ys), int(m.sum()), int(ee.sum())))
    return res


def main():
    mods = BL.load_modules()
    cohorts = build_cohorts()
    rng = np.random.default_rng(SEED)
    rows, rand_rows, km_rows = [], [], []

    for name in ["TCGA-BRCA", "TCGA-TNBC"]:
        expr, meta = cohorts[name]
        sc = score_table(expr, mods)
        t, d, ok = surv_arrays(meta, None)
        if ok.sum() < 60:
            continue
        print(f"\n[{name}] n={len(t)}, events={int(d.sum())}")

        genes_all = np.array([g for g in expr.columns])
        prog = BL.PROGRAMS["SPARC_POSTN_CORE"]
        prog = [g for g in prog if g in expr.columns]
        print(f"  程序基因 {len(prog)}: {prog}")

        cov = None
        age = pd.to_numeric(meta["age"], errors="coerce").to_numpy(float)[ok]
        age = np.where(np.isnan(age), np.nanmean(age), age)
        node = np.isin(meta["node"].astype(str).to_numpy()[ok], ["N1", "N2", "N3"]).astype(float)
        basal = (meta["pam50"].astype(str).to_numpy()[ok] == "Basal").astype(float)
        mki = expr["MKI67"].to_numpy(float)[ok] if "MKI67" in expr else None
        cov = np.column_stack([age, node, basal])
        if mki is not None:
            mki = np.where(np.isnan(mki), np.nanmean(mki), mki)
            mki_z = (mki - mki.mean()) / (mki.std(ddof=0) or np.nan)
            cov4 = np.column_stack([age, node, basal, mki_z])

        # ---- 观测值
        s_prog = BL.score_mean(expr, prog)
        hr, lo, hi, p = hr_of(s_prog, sc["CAF_ECM"], t, d, ok)
        hr_a, lo_a, hi_a, p_a = hr_of(s_prog, sc["CAF_ECM"], t, d, ok, cov)
        rows.append({"cohort": name, "stat": "observed", "n": len(t), "n_event": int(d.sum()),
                     "HR": round(hr, 3), "lo": round(lo, 3), "hi": round(hi, 3), "p": p,
                     "HR_adj3": round(hr_a, 3), "lo_adj3": round(lo_a, 3),
                     "hi_adj3": round(hi_a, 3), "p_adj3": p_a})
        print(f"  SPARC/POSTN resid: HR={hr:.3f} [{lo:.3f},{hi:.3f}] p={p:.4g}; "
              f"adj(age/node/basal) HR={hr_a:.3f} p={p_a:.4g}")
        if mki is not None:
            hr_m, lo_m, hi_m, p_m = hr_of(s_prog, sc["CAF_ECM"], t, d, ok, cov4)
            rows.append({"cohort": name, "stat": "observed_adjust_MKI67", "n": len(t),
                         "n_event": int(d.sum()), "HR": round(hr_m, 3), "lo": round(lo_m, 3),
                         "hi": round(hi_m, 3), "p": p_m,
                         "HR_adj3": np.nan, "lo_adj3": np.nan, "hi_adj3": np.nan, "p_adj3": np.nan})
            print(f"  +MKI67 校正: HR={hr_m:.3f} [{lo_m:.3f},{hi_m:.3f}] p={p_m:.4g}")

        # ---- Bootstrap
        bts = []
        for _ in range(1000):
            idx = rng.integers(0, ok.sum(), ok.sum())
            sub = pd.DataFrame({"t": t[idx], "d": d[idx]})
            # 在子样本上重算（用原始 index 对齐 score）
            s_idx = np.where(ok)[0][idx]
            r, _ = BL.rank_resid(s_prog.iloc[s_idx].reset_index(drop=True),
                                 sc["CAF_ECM"].iloc[s_idx].reset_index(drop=True))
            rv = r.to_numpy(float)
            if not np.isfinite(rv).all():
                continue
            z = (rv - rv.mean()) / (rv.std(ddof=0) or np.nan)
            bts.append(BL.cox_hr(sub["t"].to_numpy(), sub["d"].to_numpy(), z)[0])
        bts = np.array([b for b in bts if np.isfinite(b)])
        rows.append({"cohort": name, "stat": "bootstrap_HR", "n": int(ok.sum()),
                     "n_event": int(d.sum()), "HR": round(float(np.median(bts)), 3),
                     "lo": round(float(np.percentile(bts, 2.5)), 3),
                     "hi": round(float(np.percentile(bts, 97.5)), 3),
                     "p": np.nan, "HR_adj3": np.nan, "lo_adj3": np.nan,
                     "hi_adj3": np.nan, "p_adj3": np.nan})
        print(f"  Bootstrap HR 中位 {np.median(bts):.3f} "
              f"[{np.percentile(bts,2.5):.3f},{np.percentile(bts,97.5):.3f}]")

        # ---- 随机基因集基线
        k = len(prog)
        rand = np.empty(N_RANDOM)
        for i in range(N_RANDOM):
            pick = list(rng.choice(genes_all, size=k, replace=False))
            try:
                rand[i] = hr_of(BL.score_mean(expr, pick), sc["CAF_ECM"], t, d, ok)[0]
            except Exception:
                rand[i] = np.nan
        rand = rand[np.isfinite(rand)]
        emp = float((np.abs(np.log(rand)) >= abs(np.log(hr))).mean()) if hr > 0 else np.nan
        rand_rows.append({"cohort": name, "n_genes": k, "obs_HR": round(hr, 3),
                          "rand_median_HR": round(float(np.median(rand)), 3),
                          "rand_p2.5": round(float(np.percentile(rand, 2.5)), 3),
                          "rand_p97.5": round(float(np.percentile(rand, 97.5)), 3),
                          "empirical_p": round(emp, 4) if np.isfinite(emp) else np.nan,
                          "n_random": int(len(rand))})
        print(f"  随机基线 (n={len(rand)}): 中位 HR={np.median(rand):.3f}, "
              f"95% 区间 [{np.percentile(rand,2.5):.3f},{np.percentile(rand,97.5):.3f}], "
              f"经验 P={emp:.4f}")

        # ---- KM：上三分位 vs 其余（z 化残差）
        r, _ = BL.rank_resid(s_prog, sc["CAF_ECM"])
        rv = r.to_numpy(float)[ok]
        q3 = np.quantile(rv, 2 / 3)
        grp = np.where(rv >= q3, "high", "rest")
        chi, plr = BL.logrank(t, d, (grp == "high").astype(int))
        print(f"  KM T3 vs rest: chi2={chi:.2f}, p={plr:.4g}")
        curves = km_groups(t, d, grp, ["high", "rest"])
        # 三组
        qs = np.quantile(rv, [1 / 3, 2 / 3])
        g3 = np.where(rv <= qs[0], "T1", np.where(rv >= qs[1], "T3", "T2"))
        for lab, xs, ys, n_, ev in curves:
            for x, y in zip(xs, ys):
                km_rows.append({"cohort": name, "panel": "T3_vs_rest", "group": lab,
                                "time_days": x, "surv": y, "n": n_, "events": ev,
                                "logrank_p": plr})
        for lab, xs, ys, n_, ev in km_groups(t, d, g3, ["T1", "T2", "T3"]):
            for x, y in zip(xs, ys):
                km_rows.append({"cohort": name, "panel": "tertile", "group": lab,
                                "time_days": x, "surv": y, "n": n_, "events": ev,
                                "logrank_p": np.nan})
        # 5 年生存率对比
        s5 = {}
        for lab, xs, ys, n_, ev in curves:
            s5[lab] = float(ys[max(0, np.searchsorted(xs, 1825, side="right") - 1)])
        print(f"  5 年 OS: high {s5.get('high', np.nan):.1%} vs rest {s5.get('rest', np.nan):.1%}")

    pd.DataFrame(rows).to_csv(OUT / "stage9_robustness.csv", index=False)
    pd.DataFrame(rand_rows).to_csv(OUT / "stage9_random_baseline.csv", index=False)
    if km_rows:
        pd.DataFrame(km_rows).to_csv(OUT / "stage9_km_refined.csv", index=False)
    print(f"\n输出 -> {OUT}")


if __name__ == "__main__":
    main()

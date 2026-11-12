# -*- coding: utf-8 -*-
"""Stage 15 — GSE25065 NAC 治疗语境的远端无复发生存 (DRFS)。

问题: 之前所有 OS/生存主链 (METABRIC, TCGA, GSE12276) 都是初治手术语境。
GSE25065 (Popovici 2010, n=198 HER2-negative NAC, taxane+anthracycline)
是唯一在本地、有 NAC 后 DRFS 随访的队列 —— Stage 4 只用了它的 pCR 标签。

新问题: 在 NAC 治疗语境下 (含残余病 RD 亚组), perineural stroma 模块是否
预测远端复发?

分析:
  1) DRFS ~ 模块 z-score Cox (univariate + adjusted age/AJCC/pCR)
     队列: ALL (n≈198), TNBC subset (ER-/PR-/HER2-), RD-only (残余病亚组)
  2) 510 组 size-matched 随机基因集同流程 Cox → 经验 P (双侧 |log HR|)
  3) 对照模块: CD8_cytotox, TGFb_core, LRRC15_myCAF, CAF_ECM, NPP58,
     Neural_intrinsic, TLS, IFN_gamma, M2_macro, Treg

输出: 03_results/tables/bulk_extrapolation/stage15_drfs_cox.csv
       stage15_random_baseline.csv
"""
from __future__ import annotations

import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, r"D:/嗜神经浸润/02_scripts")
import bulk_lib as BL  # noqa: E402
from run_bulk_stage1 import load_gse25065  # noqa: E402

warnings.filterwarnings("ignore")
OUT = BL.TAB
SEED = 20251109
N_RANDOM = 510
MODULES = ["NPP58", "Perineural_stroma", "Neural_intrinsic", "CAF_ECM",
           "LRRC15_myCAF", "CD8_cytotox", "TGFb_core", "TLS", "IFN_gamma",
           "M2_macro", "Treg"]
RND_MODULES = [("Perineural_stroma", "Perineural_stroma"),
               ("Neural_intrinsic", "Neural_intrinsic"),
               ("CAF_ECM", "CAF_ECM"),
               ("NPP58", "NPP_up+NPP_down")]


def zstd(s: pd.Series) -> pd.Series:
    return (s - s.mean()) / s.std(ddof=0)


def cox_random_baseline(expr: pd.DataFrame, t: np.ndarray, e: np.ndarray,
                        cov: np.ndarray | None, size: int) -> np.ndarray:
    """510 组 size-matched 随机集, 同流程 Cox → HR 分布。"""
    rng = np.random.default_rng(SEED)
    genes = expr.columns.to_list()
    hrs = np.full(N_RANDOM, np.nan)
    for i in range(N_RANDOM):
        pick = rng.choice(genes, size=min(size, len(genes)), replace=False)
        s = BL.score_mean(expr, list(pick))
        z = zstd(s).to_numpy()
        ok = np.isfinite(z)
        if ok.sum() < 40 or e[ok].sum() < 8:
            continue
        try:
            cc = cov[ok] if cov is not None else None
            hr, _, _, _ = BL.cox_hr(t[ok], e[ok], z[ok], cc)
            hrs[i] = np.log(hr) if np.isfinite(hr) and hr > 0 else np.nan
        except Exception:
            continue
    return hrs


def main() -> None:
    expr, clin = load_gse25065()
    print(f"[load] expr {expr.shape}, clin {clin.shape}")

    # ---- 临床字段
    ev = clin["drfs_1_event_0_censored"].astype(str).str.strip().astype(float)
    tt = clin["drfs_even_time_years"].astype(float)
    pcr = clin["pathologic_response_pcr_rd"].astype(str).str.strip().str.upper()
    age = clin["age_years"].astype(float)
    er = clin["er_status_ihc"].astype(str).str.strip()
    pr = clin["pr_status_ihc"].astype(str).str.strip()
    her2 = clin["her2_status"].astype(str).str.strip()
    stage = clin["clinical_ajcc_stage"].astype(str).str.strip()

    ok = np.isfinite(ev.to_numpy()) & np.isfinite(tt.to_numpy()) & (tt.to_numpy() > 0)
    print(f"[DRFS] 可评估 n={ok.sum()}, 事件={int(ev.to_numpy()[ok].sum())}, "
          f"随访中位 {np.nanmedian(tt.to_numpy()[ok]):.1f} 年")

    tnbc = ((er == "N") & (pr == "N") & (her2 == "N")).to_numpy()
    rd = (pcr == "RD").to_numpy()

    strata = {
        "ALL": ok,
        "TNBC_IHC": ok & tnbc,
        "RD_only": ok & rd,
        "pCR_only": ok & (pcr == "PCR").to_numpy(),
    }
    for k, m in strata.items():
        print(f"  {k:10s} n={int(m.sum())} events={int(ev.to_numpy()[m].sum())}")

    # ---- 协变量: age + stage(有序: I=1,II=2,III=3,IV=4) + pCR(0/1)
    # 注意 startswith("II") 会把 III* 误判为 2 —— 必须先数罗马字符
    def stage_ord(s: str) -> float:
        s = s.strip().rstrip("ABC")
        return {"I": 1, "II": 2, "III": 3, "IV": 4}.get(s, np.nan)

    st_ord = stage.map(stage_ord)
    cov_full = np.column_stack([
        age.to_numpy(float),
        st_ord.to_numpy(float),
        (pcr == "PCR").to_numpy(float),
    ])

    mods = BL.load_modules()
    sc = BL.score_all(expr, mods)
    print(f"[score] 模块分 {sc.shape}")

    rows, rnd_rows = [], []
    for name, mask in strata.items():
        if mask.sum() < 40:
            continue
        t = tt.to_numpy()[mask]
        e = ev.to_numpy()[mask]
        for m in MODULES:
            if m not in sc:
                continue
            s = sc[m]
            z = zstd(s).to_numpy()[mask]
            o = np.isfinite(z)
            if o.sum() < 40 or e[o].sum() < 6:
                continue
            hr, lo, hi, p = BL.cox_hr(t[o], e[o], z[o])
            cv = cov_full[mask][o]
            cvn = np.column_stack([np.where(np.isnan(cv[:, 0]), np.nanmean(cv[:, 0]), cv[:, 0]), cv[:, 1]])
            # 恒定协变量列 (RD/pCR 层内 pcr 恒定) 会共线退化 → 剔除方差为 0 的列
            keepcv = [j for j in range(cvn.shape[1]) if np.nanstd(cvn[:, j]) > 1e-8]
            cvn = cvn[:, keepcv]
            hr2, lo2, hi2, p2 = BL.cox_hr(t[o], e[o], z[o], cvn) if keepcv else (hr, lo, hi, p)
            rows.append({"cohort": name, "module": m, "n": int(o.sum()),
                         "n_event": int(e[o].sum()),
                         "HR_per_SD": round(hr, 3), "ci_lo": round(lo, 3),
                         "ci_hi": round(hi, 3), "p": p,
                         "HR_adj": round(hr2, 3), "ci_lo_adj": round(lo2, 3),
                         "ci_hi_adj": round(hi2, 3), "p_adj": p2})
            print(f"  [{name:9s}] {m:18s} HR={hr:5.3f} [{lo:.3f},{hi:.3f}] P={p:.3g} | adj P={p2:.3g}")

        # 随机基线 (ALL 与 RD_only 两个关键层)
        if name in ("ALL", "RD_only"):
            for m, key in RND_MODULES:
                genes = (mods["NPP_up"] + mods["NPP_down"]) if key == "NPP_up+NPP_down" else mods[m]
                genes = [g for g in set(genes) if g in expr.columns]
                s = sc[m]
                z = zstd(s).to_numpy()[mask]
                o = np.isfinite(z)
                if o.sum() < 40 or e[o].sum() < 6:
                    continue
                obs_hr, _, _, _ = BL.cox_hr(t[o], e[o], z[o])
                dist = cox_random_baseline(expr.loc[mask].loc[sc.index[mask][o]], t[o], e[o],
                                           None, size=len(genes))
                ll = np.log(obs_hr)
                emp = float((np.abs(dist[np.isfinite(dist)]) >= abs(ll)).mean()) if np.isfinite(ll) else np.nan
                rnd_rows.append({"cohort": name, "module": m, "n_genes": len(genes),
                                 "obs_HR": round(obs_hr, 3),
                                 "random_median_HR": round(float(np.exp(np.nanmedian(dist))), 3),
                                 "random_p97.5_HR": round(float(np.exp(np.nanpercentile(dist, 97.5))), 3),
                                 "empirical_p": emp,
                                 "n_random": int(np.isfinite(dist).sum())})
                print(f"  [Stage15-RND {name}] {m:18s} obs HR={obs_hr:.3f} "
                      f"vs random med {np.exp(np.nanmedian(dist)):.3f}, emp P={emp:.3f}")

    pd.DataFrame(rows).to_csv(OUT / "stage15_drfs_cox.csv", index=False)
    pd.DataFrame(rnd_rows).to_csv(OUT / "stage15_random_baseline.csv", index=False)
    print(f"\n输出 -> {OUT}/stage15_drfs_cox.csv, stage15_random_baseline.csv")


if __name__ == "__main__":
    main()

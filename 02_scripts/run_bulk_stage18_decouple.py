"""Stage 18 (bulk): CD8 效应/趋化解耦 + EMILIN1 反向对照.

用户第三轮方案第 2/3 条:
  "效应 T 细胞标志物在, 但迁移/趋化轴被 myCAF 压制" ——
    读数: CYT(GZMA/PRF1/GZMB) / IFN-gamma / 趋化轴(CXCL9/10/11, CCL5, CCL19/21, CXCL13) / 迁移受体(CCR7/CXCR3)
    检验: partial Spearman(读数, stroma 模块) 控 CAF/ECM 总纤维化 + 随机集基线
  EMILIN1 (TGF-beta 抑制性 ECM 组分) 反向对照: 预期与 CD8 正关联, 与 stroma 模块负/零.

队列: METABRIC / GSE163882 / GSE25065 (TCGA 靶向面板缺 CCL5/CCR7/EMILIN1, 不纳入).

输出: 03_results/tables/bulk_extrapolation/stage18_cd8_decoupling.csv
"""
from __future__ import annotations

import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

warnings.filterwarnings("ignore")
sys.path.insert(0, r"D:/嗜神经浸润/02_scripts")
import bulk_lib as BL                                    # noqa: E402
from run_bulk_stage1 import load_gse163882, load_gse25065  # noqa: E402

import importlib.util as _ilu
spec = _ilu.spec_from_file_location(
    "s11", Path(r"D:/嗜神经浸润/02_scripts/run_bulk_stage11_metabric.py"))
s11 = _ilu.module_from_spec(spec)
sys.modules["s11"] = s11
spec.loader.exec_module(s11)

RAW = Path(r"D:/嗜神经浸润/01_data/raw")
TAB = Path(r"D:/嗜神经浸润/03_results/tables/bulk_extrapolation")
SEED = 20261118
N_RANDOM = 500

# ------------------------------------------------------------------ 读数
READOUTS = {
    "CYT":            ["GZMA", "PRF1", "GZMB"],
    "CD8_cytotox":    BL.load_modules()["CD8_cytotox"],
    "IFN_gamma":      BL.load_modules()["IFN_gamma"],
    "Chemokine_axis": ["CXCL9", "CXCL10", "CXCL11", "CCL5", "CCL19", "CCL21", "CXCL13"],
    "CXCL9_10_11":    ["CXCL9", "CXCL10", "CXCL11"],
    "CCL_homeostatic": ["CCL19", "CCL21"],
    "T_migration":    ["CCR7", "CXCR3", "CXCR5"],
    "TLS":            BL.load_modules()["TLS"],
    "EMILIN1":        ["EMILIN1"],
    "LOX_path":       ["LOX", "LOXL2", "PLOD2"],
}
FOCUS_MODS = ["Perineural_stroma", "Neural_intrinsic", "CAF_ECM", "CD8_cytotox", "TLS"]


def cohorts() -> dict:
    out = {}
    out["GSE163882"], _ = load_gse163882()
    out["GSE25065"], _ = load_gse25065()
    out["METABRIC"], _ = s11.load_metabric()
    return out


def score(expr: pd.DataFrame, genes: list[str]) -> pd.Series:
    """小基因集(<4)兼容打分: 直接 mean-z, 不做 BL.score_mean 的 min_genes=4 过滤。"""
    avail = [g for g in genes if g in expr.columns]
    if not avail:
        return pd.Series(np.nan, index=expr.index)
    sub = expr[avail]
    z = (sub - sub.mean()) / (sub.std(ddof=0) + 1e-12)
    return z.mean(axis=1)


def random_rho(expr: pd.DataFrame, x: pd.Series, size: int,
               cov: pd.DataFrame, exclude: set[str]) -> np.ndarray:
    """size-matched 随机基因集与 x 的偏 Spearman 分布(控 cov)。

    池排除全部模块/程序/读数基因 —— METABRIC 本地仅 210 个靶向基因,
    其池固有偏倚(见脚本头注), 该队列随机基线仅作描述性参照。
    """
    rng = np.random.default_rng(SEED)
    pool = np.array([g for g in expr.columns if g not in exclude])
    rhos = np.full(N_RANDOM, np.nan)
    for i in range(N_RANDOM):
        pick = rng.choice(pool, size=min(size, len(pool)), replace=False)
        s = score(expr, list(pick))
        if np.isnan(s).all():
            continue
        r, p, n = BL.partial_spearman(s, x, cov)
        rhos[i] = r
    return rhos


def main() -> None:
    mods = BL.load_modules()
    caf = mods["CAF_ECM"]
    # 随机池排除所有模块/程序/读数基因
    exclude = set()
    for v in mods.values():
        exclude.update(v)
    for v in BL.PROGRAMS.values():
        exclude.update(v)
    for v in READOUTS.values():
        exclude.update(v)
    rows = []
    for cname, expr in cohorts().items():
        gset = set(expr.columns)
        stroma = score(expr, mods["Perineural_stroma"])
        cov = pd.DataFrame({"CAF_ECM": score(expr, caf)})
        biased_pool = expr.shape[1] <= 300   # METABRIC 210 基因靶向文件
        print(f"\n[{cname}] n={expr.shape[0]} genes={expr.shape[1]}"
              f"{' [pool-biased: target-gene file]' if biased_pool else ''}")
        for rname, glist in READOUTS.items():
            avail = [g for g in glist if g in gset]
            need = 1 if len(glist) == 1 else max(2, len(glist) // 2)
            if len(avail) < need:
                print(f"  {rname}: skipped (avail {len(avail)}/{len(glist)})")
                continue
            s = score(expr, avail)
            if np.isnan(s).all():
                continue
            r, p, n = BL.partial_spearman(s, stroma, cov)
            if biased_pool:
                rows.append(dict(
                    cohort=cname, readout=rname, n_genes_used=len(avail),
                    n_genes_total=len(glist), n=int(n),
                    rho_partial=float(r), p_partial=float(p),
                    random_median=np.nan, random_abs95=np.nan, p_empirical=np.nan))
                print(f"  {rname:16s} ρ={r:+.3f} p={p:.3g} (no random null)")
                continue
            rb = random_rho(expr, stroma, len(avail), cov, exclude)
            rb = rb[~np.isnan(rb)]
            emp_p = float((1 + int((np.abs(rb) >= abs(r)).sum())) / (1 + len(rb)))
            rows.append(dict(
                cohort=cname, readout=rname, n_genes_used=len(avail),
                n_genes_total=len(glist), n=int(n),
                rho_partial=float(r), p_partial=float(p),
                random_median=float(np.median(rb)), random_abs95=float(np.quantile(np.abs(rb), 0.95)),
                p_empirical=emp_p,
            ))
            print(f"  {rname:16s} ρ={r:+.3f} p={p:.3g} | rand med {np.median(rb):+.3f} empP={emp_p:.3f}")
        # EMILIN1 调制: EMILIN1 与 CD8/趋化轴的偏相关(控 CAF_ECM)
        if "EMILIN1" in gset:
            em = score(expr, ["EMILIN1"])
            for rname in ["CD8_cytotox", "Chemokine_axis", "CYT"]:
                glist = READOUTS[rname]
                avail = [g for g in glist if g in gset]
                if len(avail) < 2:
                    continue
                s = score(expr, avail)
                r, p, n = BL.partial_spearman(s, em, cov)
                rows.append(dict(
                    cohort=cname, readout=f"EMILIN1~{rname}", n_genes_used=1,
                    n_genes_total=1, n=int(n),
                    rho_partial=float(r), p_partial=float(p),
                    random_median=np.nan, random_abs95=np.nan, p_empirical=np.nan,
                ))
            # EMILIN1 vs stroma module
            r, p, n = BL.partial_spearman(em, stroma, cov)
            rows.append(dict(cohort=cname, readout="EMILIN1~Perineural_stroma",
                             n_genes_used=1, n_genes_total=1, n=int(n),
                             rho_partial=float(r), p_partial=float(p),
                             random_median=np.nan, random_abs95=np.nan, p_empirical=np.nan))

    out = pd.DataFrame(rows)
    out.to_csv(TAB / "stage18_cd8_decoupling.csv", index=False)
    print("\n=== cross-cohort consistency (rho_partial) ===")
    piv = out[~out.readout.str.contains("~")].pivot_table(index="readout", columns="cohort", values="rho_partial")
    print(piv.round(3).to_string())
    print("\nsaved:", TAB / "stage18_cd8_decoupling.csv")


if __name__ == "__main__":
    main()

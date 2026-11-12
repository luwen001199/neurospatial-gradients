"""Stage 19 (bulk): Schwann/卫星胶质评分 (SC score) 对照投影.

用户第三轮方案第 8 条: 把 Schwann 标记集投到 bulk 队列,
预期 SC score 与 stroma 模块高相关、与 neural-intrinsic 低相关 →
支持"bulk 层的『神经信号』很大一部分是 Schwann/基质侧, 不是神经元"。

SC signature (canonical Schwann/myelinating glia):
  SOX10, MPZ, MBP, PMP22, PLP1, NCAM1, PRNP, NEFL, NEFH, S100B
敏感性: 去除 SOX10 (乳腺中同时标记 basal 肿瘤细胞) 的 9 基因版本。

队列: GSE163882 / GSE25065 (全基因组阵列) + TCGA-BRCA 面板 (可用子集)。
      METABRIC 210 靶向基因文件缺 SOX10/MBP/PMP22/PRNP/NEFH/NCAM1, 不纳入。

输出: 03_results/tables/bulk_extrapolation/stage19_schwann.csv
"""
from __future__ import annotations

import sys
import warnings

import numpy as np
import pandas as pd
from scipy import stats

warnings.filterwarnings("ignore")
sys.path.insert(0, r"D:/嗜神经浸润/02_scripts")
import bulk_lib as BL                                    # noqa: E402
from run_bulk_stage1 import load_gse163882, load_gse25065  # noqa: E402

TAB = BL.TAB
SEED = 20261119
N_RANDOM = 500

SC_FULL = ["SOX10", "MPZ", "MBP", "PMP22", "PLP1", "NCAM1", "PRNP", "NEFL", "NEFH", "S100B"]
SC_NOSOX = [g for g in SC_FULL if g != "SOX10"]

NERVE_CALL = ["UCHL1", "TUBB3", "NEFH", "PRNP", "GAP43", "MPZ", "NCAM1"]  # 原稿神经spot定义基因
BASAL_KRT = ["KRT5", "KRT14", "KRT17"]


def score(expr, genes):
    avail = [g for g in genes if g in expr.columns]
    if not avail:
        return None, []
    sub = expr[avail]
    z = (sub - sub.mean()) / (sub.std(ddof=0) + 1e-12)
    return z.mean(axis=1), avail


def random_null(expr, target: pd.Series, size: int, exclude: set, cov=None):
    rng = np.random.default_rng(SEED)
    pool = np.array([g for g in expr.columns if g not in exclude])
    out = []
    for _ in range(N_RANDOM):
        pick = rng.choice(pool, size=min(size, len(pool)), replace=False)
        s, _ = score(expr, list(pick))
        if s is None or np.isnan(s).all():
            continue
        if cov is not None:
            r, p, n = BL.partial_spearman(s, target, cov)
        else:
            r = stats.spearmanr(s, target, nan_policy="omit")[0]
        out.append(r)
    return np.asarray(out)


def main() -> None:
    mods = BL.load_modules()
    exclude = set()
    for v in mods.values():
        exclude.update(v)
    for v in BL.PROGRAMS.values():
        exclude.update(v)
    exclude.update(SC_FULL + NERVE_CALL + BASAL_KRT)

    cohorts = {}
    cohorts["GSE163882"], _ = load_gse163882()
    cohorts["GSE25065"], _ = load_gse25065()
    expr_tcga = pd.read_csv(BL.RAW / "tcga_brca" / "tcga_brca_target_genes.csv", index_col=0).astype(float)
    cohorts["TCGA"] = expr_tcga

    rows = []
    for cname, expr in cohorts.items():
        sc_full, avail_f = score(expr, SC_FULL)
        sc_nosox, avail_n = score(expr, SC_NOSOX)
        stroma, _ = score(expr, mods["Perineural_stroma"])
        neural, _ = score(expr, mods["Neural_intrinsic"])
        caf, _ = score(expr, mods["CAF_ECM"])
        # SOX10 单基因 + basal 角蛋白对照 (SOX10 混杂探针)
        sox10 = None
        if "SOX10" in expr.columns:
            c = expr["SOX10"]
            sox10 = (c - c.mean()) / (c.std(ddof=0) + 1e-12)
        basal, avail_b = score(expr, BASAL_KRT)

        print(f"\n[{cname}] SC genes avail: full {len(avail_f)}/10, nosox {len(avail_n)}/9")
        for tag, sc, avail in [("SC_full", sc_full, avail_f), ("SC_nosox", sc_nosox, avail_n)]:
            if sc is None or np.isnan(sc).all():
                print(f"  {tag}: not measurable")
                continue
            r_st, p_st = stats.spearmanr(sc, stroma, nan_policy="omit")
            r_ne, p_ne = stats.spearmanr(sc, neural, nan_policy="omit")
            r_ca, p_ca = stats.spearmanr(sc, caf, nan_policy="omit")
            null = random_null(expr, stroma, len(avail), exclude)
            z = (r_st - null.mean()) / (null.std() + 1e-12)
            p_emp = float(min((null <= r_st).mean(), (null >= r_st).mean()) * 2)
            rows.append(dict(cohort=cname, signature=tag, n_genes=len(avail),
                             rho_vs_stroma=float(r_st), p_stroma=float(p_st),
                             rho_vs_neural=float(r_ne), p_neural=float(p_ne),
                             rho_vs_CAFECM=float(r_ca), p_CAFECM=float(p_ca),
                             null_median=float(np.median(null)), null_p975=float(np.percentile(null, 97.5)),
                             z_vs_null=float(z), p_empirical=p_emp))
            print(f"  {tag:9s} vs stroma ρ={r_st:+.3f} (P={p_st:.2g}) | vs neural ρ={r_ne:+.3f} "
                  f"| vs CAF/ECM ρ={r_ca:+.3f} | null {np.median(null):+.3f} z={z:+.1f} empP={p_emp:.3g}")
        if sox10 is not None:
            r, p = stats.spearmanr(sox10, stroma, nan_policy="omit")
            rb, _ = stats.spearmanr(sox10, basal) if basal is not None else (np.nan, np.nan)
            rows.append(dict(cohort=cname, signature="SOX10_single", n_genes=1,
                             rho_vs_stroma=float(r), p_stroma=float(p),
                             rho_vs_neural=np.nan, p_neural=np.nan,
                             rho_vs_CAFECM=np.nan, p_CAFECM=np.nan,
                             null_median=np.nan, null_p975=np.nan,
                             z_vs_null=np.nan, p_empirical=np.nan))
            print(f"  SOX10 single vs stroma ρ={r:+.3f} (P={p:.2g}) | vs basal-KRT ρ={rb:+.3f}")

    out = pd.DataFrame(rows)
    out.to_csv(TAB / "stage19_schwann.csv", index=False)
    print("\nsaved:", TAB / "stage19_schwann.csv")
    print(out.to_string(index=False))


if __name__ == "__main__":
    main()

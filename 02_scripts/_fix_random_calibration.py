"""修复 stage17/axis 随机基线校准 bug.

bug: R.groupby("name").agg(size=("rho","size")) 的 "size" 是**样本计数**而非基因集大小
     (named-agg 把 "size" 当聚合函数), 导致 bucket 按样本数分组、size-matching 失效,
     所有 p_emp 撞分辨率地板.

修复: 从随机集名 (RAND|s{size}|{idx}) 解析真实大小, 按基因集大小分桶,
     与观测读数按 |size 对数最近邻| 匹配.
重算: stage17_gradient_meta.csv (本 stage) 与 axis_meta_summary.csv (更正版).
"""
from __future__ import annotations

import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

TAB17 = Path(r"D:/嗜神经浸润/03_results/tables/spatial_tlscaf")
TABAX = Path(r"D:/嗜神经浸润/03_results/tables")


def bh(p):
    p = np.asarray(p, dtype=float)
    n = len(p)
    o = np.argsort(p)
    q = np.empty(n)
    prev = 1.0
    for i in range(n - 1, -1, -1):
        prev = min(prev, p[o[i]] * n / (i + 1))
        q[o[i]] = prev
    return np.minimum(q, 1.0)


def recalibrate(P: pd.DataFrame, members: dict, out_csv: Path,
                min_n: int = 10) -> pd.DataFrame:
    """P: per-sample 表; members: {name: gene-set size}."""
    is_rand = P.kind == "random"
    R = P[is_rand].copy()
    R["gsize"] = R.name.str.extract(r"RAND\|s(\d+)\|")[0].astype(float)
    Rg = R.groupby("gsize").agg(
        n_sets=("name", "nunique"),
    )
    # 每个随机集的 meta 统计量 (跨样本)
    rmeta = R.groupby("name").agg(
        gsize=("gsize", "first"),
        mean_rho=("rho", "mean"),
        stouffer=("u_msr", lambda u: float(np.sum(u) / np.sqrt(len(u)))),
    ).reset_index()
    rmeta["abs_rho"] = rmeta.mean_rho.abs()

    buckets = {int(s): g for s, g in rmeta.groupby("gsize")}
    sizes = np.array(sorted(buckets))
    print(f"[random] buckets: {[(s, len(g)) for s, g in buckets.items()]}")

    meta_rows = []
    for name, g in P[~is_rand].groupby("name", sort=True):
        n = g["rho"].notna().sum()
        if n < min_n:
            continue
        z = np.arctanh(np.clip(g["rho"], -0.999, 0.999))
        se = float(z.std(ddof=1) / np.sqrt(len(z)))
        rho_mean = float(np.tanh(z.mean()))
        ci = np.tanh([z.mean() - 1.96 * se, z.mean() + 1.96 * se])
        zz = np.asarray(g["u_msr"], dtype=float)
        stou = float(zz.sum() / np.sqrt(len(zz)))
        tstat, p_rho = stats.ttest_1samp(z, 0.0)
        meta_rows.append(dict(
            name=name, size=members.get(name, 1), n_samples=int(n),
            rho_mean=rho_mean, rho_lo=ci[0], rho_hi=ci[1], p_rho=p_rho,
            stouffer_Z=stou, p_stouffer=float(2 * stats.norm.sf(abs(stou))),
            mean_u=float(zz.mean()),
        ))
    M = pd.DataFrame(meta_rows)
    if M.empty:
        return M

    def nearest(s):
        return int(sizes[np.argmin(np.abs(np.log(sizes) - np.log(max(1.0, s))))])

    cal = []
    for r in M.itertuples():
        b = nearest(float(r.size))
        g = buckets[b]
        br = g.abs_rho.to_numpy(float)
        mr = g.mean_rho.to_numpy(float)
        ar = abs(float(r.rho_mean))
        zr = (ar - br.mean()) / (br.std(ddof=1) + 1e-12)
        sgn = np.sign(r.rho_mean) if r.rho_mean != 0 else 1.0
        same_dir = int((np.sign(mr) == np.sign(r.rho_mean)).sum())
        # 富集方向 one-sided: 与观测同向且幅度 >= |观测| 的随机集比例
        n_exceed_dir = int((mr * sgn >= ar).sum())
        cal.append(dict(
            bucket=b, n_random=int(len(g)),
            baseline_rho=float(br.mean()), excess_rho=ar - float(br.mean()),
            z_vs_random=float(zr),
            p_emp_one=float((1 + int((br >= ar).sum())) / (1 + len(br))),
            p_emp_dir=float((1 + n_exceed_dir) / (1 + len(g))),
            n_same_dir=same_dir,
        ))
    C = pd.DataFrame(cal)
    M = pd.concat([M.reset_index(drop=True), C], axis=1)
    M["fdr_rho"] = bh(M["p_rho"].to_numpy())
    M["fdr_emp_dir"] = bh(M["p_emp_dir"].to_numpy())
    M = M.sort_values("rho_mean")
    M.to_csv(out_csv, index=False)
    return M


# ---------------------------------------------------------------- stage17
sys.path.insert(0, r"D:/嗜神经浸润/02_scripts")
import run_bulk_stage17_tlscaf as S17  # noqa: E402

P17 = pd.read_csv(TAB17 / "stage17_gradient_per_sample.csv.gz")
members17 = {k: len(v) for k, v in S17.build_blocks([])["func"].items()}
# build_blocks([]) 的 func 基因集大小与跑批时一致(基因列表硬编码)
M17 = recalibrate(P17, members17, TAB17 / "stage17_gradient_meta.csv")
print(M17[["name", "size", "n_samples", "rho_mean", "rho_lo", "rho_hi",
           "baseline_rho", "p_emp_dir", "fdr_emp_dir"]].to_string(index=False))

# ---------------------------------------------------------------- axis (更正版)
PAX = pd.read_csv(r"D:/嗜神经浸润/01_data/processed/axis_per_sample.csv.gz")
# axis 的读数成员大小
import importlib
_npp_spec = importlib.util.spec_from_file_location(
    "npp_split", r"D:/嗜神经浸润/02_scripts/run_npp_split.py")
_npp = importlib.util.module_from_spec(_npp_spec)
_npp_spec.loader.exec_module(_npp)
PS = _npp.MODULES["Perineural_stroma"]
import niche_lib as NL  # noqa: E402
members_ax = {f"MOD|Perineural_stroma": len(PS)}
members_ax[f"MOD|Collagen_core"] = sum(1 for g in PS if g.startswith("COL"))
for g in ["LRRC15", "ACTA2", "POSTN", "TNC", "COL1A1", "SPARC", "FN1", "LOXL2"]:
    members_ax[f"GENE|{g}"] = 1
members_ax["AXIS|myCAF"] = 20
members_ax["AXIS|iCAF"] = 27
members_ax["AXIS|apCAF"] = 8
members_ax["AXIS|LRRC15p"] = 1
for ct in ["CD8_T", "T_lymphocyte", "Macrophage_monocyte"]:
    members_ax[f"CT|{ct}"] = len(NL.CELLTYPE_MARKERS[ct])
MAX = recalibrate(PAX, members_ax, TABAX / "axis_meta_summary.corrected.csv")
print(MAX[["name", "size", "n_samples", "rho_mean", "baseline_rho",
           "p_emp_dir", "fdr_emp_dir"]].to_string(index=False))

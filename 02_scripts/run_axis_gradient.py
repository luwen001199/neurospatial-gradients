"""主线1/2/3 的读数梯度扫描: stroma module / myCAF 标志基因 / CAF 四轴 / 免疫读数,
全部经偏 Spearman + MSR 随机化 + 大小匹配随机集基线三重校准。

复用 run_niche_spatial 的 process_sample (秩空间偏 Spearman, Gram-Schmidt 残差化
log UMI / pct_mt, Moran 谱随机化零分布) 与 make_random_sets (大小匹配基线)。

读数
----
stroma   : Perineural_stroma 模块 (60) + collagen 子模块 + 8 个 myCAF/ECM 标志基因
axis     : myCAF (20, 剔除 LRRC15) / iCAF (27) / apCAF (8) / LRRC15+ (1) 四轴
immune   : CD8_T / T_lymphocyte / Macrophage_monocyte 细胞类型读数
random   : 大小匹配随机基因集 (sizes 1,8,20,27,60 x 80 = 400)

输出
----
01_data/processed/axis_per_sample.csv.gz   每样本每读数 rho / u_msr / mu_sur / sd_sur
01_data/processed/axis_ring_profiles.csv   环剖面 (衰减曲线用)
03_results/tables/axis_meta_summary.csv    跨样本 meta + 双重校准 (MSR + 随机基线)

用法: python 02_scripts/run_axis_gradient.py [--nsim 199] [--only GSM...]
"""
from __future__ import annotations

import argparse
import importlib.util
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, r"D:/嗜神经浸润/00_pipeline")
sys.path.insert(0, r"D:/嗜神经浸润/02_scripts")

from nsgrad import config as CFG                      # noqa: E402
from run_niche_spatial import (                        # noqa: E402
    make_random_sets, process_sample, union_expressed_genes,
)
import niche_lib as NL                                # noqa: E402

DIR_H5 = CFG.DIR_PROC / "spot_level"
PROC = CFG.DIR_PROC
TAB = CFG.DIR_TAB
FDR = 0.05
RING_ORDER = ["r1_0_50", "r2_50_100", "r3_100_200", "r4_200_400", "out_gt400"]

# 模块基因与主分析严格一致
_spec = importlib.util.spec_from_file_location("npp_split", Path(__file__).parent / "run_npp_split.py")
_npp = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_npp)
PERINEURAL_STROMA = _npp.MODULES["Perineural_stroma"]

SINGLE_GENES = ["LRRC15", "ACTA2", "POSTN", "TNC", "COL1A1", "SPARC", "FN1", "LOXL2"]
COLLAGEN_CORE = [g for g in PERINEURAL_STROMA if g.startswith("COL")]
APCAF = ["HLA-DRA", "HLA-DRB1", "HLA-DPA1", "HLA-DPB1", "CD74", "CIITA", "LAMP3", "CD70"]


def build_blocks(pool: list[str]) -> dict:
    stroma = {
        "MOD|Perineural_stroma": list(PERINEURAL_STROMA),
        "MOD|Collagen_core": COLLAGEN_CORE,
    }
    for g in SINGLE_GENES:
        stroma[f"GENE|{g}"] = [g]
    axis = {
        "AXIS|myCAF": [g for g in NL.CELLTYPE_MARKERS["CAF_myofibroblastic"] if g != "LRRC15"],
        "AXIS|iCAF": list(NL.CELLTYPE_MARKERS["CAF_inflammatory_dev"]),
        "AXIS|apCAF": list(APCAF),
        "AXIS|LRRC15p": ["LRRC15"],
    }
    immune = {
        "CT|CD8_T": list(NL.CELLTYPE_MARKERS["CD8_T"]),
        "CT|T_lymphocyte": list(NL.CELLTYPE_MARKERS["T_lymphocyte"]),
        "CT|Macrophage_monocyte": list(NL.CELLTYPE_MARKERS["Macrophage_monocyte"]),
    }
    sizes = []
    for s, rep in [(1, 80), (8, 80), (20, 80), (27, 80), (60, 80)]:
        sizes.extend([s] * rep)
    rsets = make_random_sets(pool, sizes, 1, seed=7)
    return {"stroma": stroma, "axis": axis, "immune": immune, "random": rsets}


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


def main(nsim: int, only: str | None):
    files = sorted(DIR_H5.glob("*.h5ad"))
    if only:
        files = [f for f in files if only in f.stem]
    pool = sorted(union_expressed_genes(files))
    blocks = build_blocks(pool)
    print(f"[init] pool={len(pool)} readouts="
          f"{sum(len(blocks[k]) for k in ('stroma','axis','immune'))} random={len(blocks['random'])}")

    rows, rings = [], []
    t0 = time.time()
    for i, p in enumerate(files, 1):
        try:
            res = process_sample(p, nsim, 20261105, 40, blocks_global=blocks)
        except Exception as e:
            print(f"  [skip] {p.stem}: {type(e).__name__}: {e}")
            continue
        o, rg = res
        rows.append(o)
        if len(rg):
            rings.append(rg)
        print(f"  [{i}/{len(files)}] {p.stem}: K={len(o)} n={o['n_spot'].iloc[0]} "
              f"{time.time()-t0:.0f}s", flush=True)
    P = pd.concat(rows, ignore_index=True)
    P.to_csv(PROC / "axis_per_sample.csv.gz", index=False, compression="gzip")
    if rings:
        pd.concat(rings, ignore_index=True).to_csv(PROC / "axis_ring_profiles.csv", index=False)

    # ---------------- meta: Fisher-z / Stouffer / BH ----------------
    meta_rows = []
    for name, g in P[P.kind != "random"].groupby("name", sort=True):
        n = g["rho"].notna().sum()
        if n < 10:
            continue
        z = np.arctanh(np.clip(g["rho"], -0.999, 0.999))
        se = float(z.std(ddof=1) / np.sqrt(len(z)))    # 随机效应 SE
        rho_mean = float(np.tanh(z.mean()))
        ci = np.tanh([z.mean() - 1.96 * se, z.mean() + 1.96 * se])
        zz = np.asarray(g["u_msr"], dtype=float)
        stou = float(zz.sum() / np.sqrt(len(zz)))
        # 随机效应 meta: 对逐样本 Fisher-z 做单样本 t 检验
        tstat, p_rho = stats.ttest_1samp(z, 0.0)
        kind = name.split("|")[0]
        members = {"MOD": {"MOD|Perineural_stroma": len(PERINEURAL_STROMA),
                           "MOD|Collagen_core": len(COLLAGEN_CORE)},
                   "AXIS": {"AXIS|myCAF": 20, "AXIS|iCAF": 27, "AXIS|apCAF": 8, "AXIS|LRRC15p": 1},
                   "CT": {k: len(NL.CELLTYPE_MARKERS[k.split("|")[1]]) for k in
                          ["CT|CD8_T", "CT|T_lymphocyte", "CT|Macrophage_monocyte"]},
                   "GENE": {f"GENE|{g2}": 1 for g2 in SINGLE_GENES}}
        size = members.get(kind, {}).get(name, 1)
        meta_rows.append(dict(
            name=name, kind=kind, size=size, n_samples=int(n),
            rho_mean=rho_mean, rho_lo=ci[0], rho_hi=ci[1], p_rho=p_rho,
            stouffer_Z=stou, p_stouffer=float(2 * stats.norm.sf(abs(stou))),
            mean_u=float(zz.mean()),
        ))
    M = pd.DataFrame(meta_rows)
    if M.empty:
        print("[meta] no readout with n>=10 samples; skip meta stage")
        return

    # ---------------- 大小匹配随机基线 ----------------
    R = P[P.kind == "random"]
    R = R.groupby("name").agg(size=("rho", "size"), mean_rho=("rho", "mean"),
                              stouffer_Z=("u_msr", lambda u: float(
                                  np.sum(u) / np.sqrt(len(u)))))
    R = R.rename(columns={"size": "size_"}).reset_index()
    R["size"] = pd.to_numeric(R["size_"], errors="coerce")
    # 合并零分布: 桶内 z 标准化后合并 (分辨率 1/401)
    zr_all, zz_all = [], []
    for s, g in R.groupby("size"):
        br = np.abs(g["mean_rho"].to_numpy(float))
        bz = np.abs(g["stouffer_Z"].to_numpy(float))
        zr_all.extend(((np.abs(g["mean_rho"]) - br.mean()) / (br.std(ddof=1) + 1e-12)).tolist())
        zz_all.extend(((np.abs(g["stouffer_Z"]) - bz.mean()) / (bz.std(ddof=1) + 1e-12)).tolist())
    zr_all = np.sort(np.asarray(zr_all))
    zz_all = np.sort(np.asarray(zz_all))
    buckets = {int(s): g for s, g in R.groupby("size")}
    sizes = np.array(sorted(buckets))

    def nearest(s):
        return int(sizes[np.argmin(np.abs(np.log(sizes) - np.log(max(1.0, s))))])

    cal = []
    for r in M.itertuples():
        b = nearest(float(r.size))
        g = buckets[b]
        br = np.abs(g["mean_rho"].to_numpy(float))
        bz = np.abs(g["stouffer_Z"].to_numpy(float))
        ar, az = abs(float(r.rho_mean)), abs(float(r.stouffer_Z))
        zr = (ar - br.mean()) / (br.std(ddof=1) + 1e-12)
        zz = (az - bz.mean()) / (bz.std(ddof=1) + 1e-12)
        cal.append(dict(
            bucket=b, n_random=int(len(g)),
            baseline_rho=float(br.mean()), excess_rho=ar - float(br.mean()),
            z_vs_random=float(zr),
            p_emp_rho=float((1 + int((zr_all >= zr).sum())) / (1 + len(zr_all))),
            p_bucket=float((1 + int((br >= ar).sum())) / (1 + len(br))),
        ))
    C = pd.DataFrame(cal)
    M = pd.concat([M.reset_index(drop=True), C], axis=1)
    M["fdr_rho"] = bh(M["p_rho"].to_numpy())
    M["fdr_stouffer"] = bh(M["p_stouffer"].to_numpy())
    M["fdr_emp_rho"] = bh(M["p_emp_rho"].to_numpy())
    M["pass_msr"] = M["fdr_stouffer"] < FDR
    M["pass_random"] = M["fdr_emp_rho"] < FDR
    M["pass_both"] = M["pass_msr"] & M["pass_random"]
    M = M.sort_values("rho_mean")
    M.to_csv(TAB / "axis_meta_summary.csv", index=False)
    print(M[["name", "size", "n_samples", "rho_mean", "fdr_rho", "stouffer_Z",
             "fdr_stouffer", "baseline_rho", "p_emp_rho", "fdr_emp_rho", "pass_both"]]
          .round(4).to_string(index=False))
    print(f"\nsaved: {TAB/'axis_meta_summary.csv'}; "
          f"both-pass: {int(M['pass_both'].sum())}/{len(M)}")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--nsim", type=int, default=199)
    ap.add_argument("--only", type=str, default=None)
    a = ap.parse_args()
    main(a.nsim, a.only)

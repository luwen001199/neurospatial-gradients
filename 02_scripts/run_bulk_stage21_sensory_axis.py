"""Stage 21 (空间层): 与 Cell 2026 (Zhang SW et al., PMID 41650969) 的数据级对话.

Cell 2026 提出 NGF→感觉神经→CGRP→RAMP1+CAF→胶原屏障→免疫排斥 (TNBC 360例+小鼠).
本 stage 在原稿 43 例 Visium 上做两件此前未做的事:

Part A — 感觉神经信号轴的近神经距离梯度 (三重校准管线):
  CGRP_axis  (CALCA/CALCRL/RAMP1/CALCR)  — Cell 2026 核心轴的转录组代理
  CGRP_core  (CALCA/RAMP1)               — 配体+受体最小对
  NGF_TrkA   (NGF/NTRK1/NGFR/SORT1)      — 神经营养驱动轴
  Sensory_nf (TRPV1/SCN9A/SCN10A/S100B)  — 感觉神经特异标记 (无 NGF/UCHL1 混杂)
  单基因: CALCA RAMP1 CALCRL NGF NTRK1 NGFR TRPV1 SCN9A SCN10A S100B
  size-matched 随机集 [2,3,4,6,9]x60, 修正后 size-bucketed 校准 (同 Stage 17).

Part B — TLS 近/远功能比 (正式关闭"免疫聚集但无能"假说):
  TLS+ spots (TLS_ext 评分 > 样本中位) 内, 近神经(<=150um) vs 远(>=300um):
    CYT  效应轴  GZMA/GZMB/PRF1/NKG7/IFNG
    EXH  耗竭轴  PDCD1/CTLA4/LAG3/HAVCR2/TIGIT
    TGFB 响应轴  SERPINE1/TGFB1/TGFB2/TGFBR2/SMAD3
    以及 CYT-EXH 差 (效应/抑制比代理)
  逐样本 Mann-Whitney; meta 层单样本 t 检验.

输出: 03_results/tables/spatial_tlscaf/
  stage21_axis_gradient_meta.csv
  stage21_axis_gradient_per_sample.csv.gz
  stage21_tls_function_per_sample.csv
  stage21_tls_function_meta.csv
"""
from __future__ import annotations

import sys
import time
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
from scipy import stats

warnings.filterwarnings("ignore")
sys.path.insert(0, r"D:/嗜神经浸润/00_pipeline")
sys.path.insert(0, r"D:/嗜神经浸润/02_scripts")

from nsgrad import config as CFG                      # noqa: E402
from run_niche_spatial import (                        # noqa: E402
    make_random_sets, process_sample, union_expressed_genes,
)

DIR_H5 = CFG.DIR_PROC / "spot_level"
TAB = Path(r"D:/嗜神经浸润/03_results/tables/spatial_tlscaf")
TAB.mkdir(parents=True, exist_ok=True)

# ------------------------------------------------------------------ 读数定义
CGRP_AXIS = ["CALCA", "CALCRL", "RAMP1", "CALCR"]
CGRP_CORE = ["CALCA", "RAMP1"]
NGF_TRKA = ["NGF", "NTRK1", "NGFR", "SORT1"]
SENSORY_NF = ["TRPV1", "SCN9A", "SCN10A", "S100B"]
SINGLE_GENES = ["CALCA", "RAMP1", "CALCRL", "CALCR", "NGF", "NTRK1", "NGFR",
                "TRPV1", "SCN9A", "SCN10A", "S100B"]

TLS_EXT = ["CCL21", "CCL19", "CXCL13", "LAMP3", "MS4A1", "CD79A",
           "MEF2B", "POU2AF1", "BANK1"]
CYT = ["GZMA", "GZMB", "PRF1", "NKG7", "IFNG"]
EXH = ["PDCD1", "CTLA4", "LAG3", "HAVCR2", "TIGIT"]
TGFB = ["SERPINE1", "TGFB1", "TGFB2", "TGFBR2", "SMAD3"]


def build_blocks(pool: list[str]) -> dict:
    func = {
        "MOD|CGRP_axis": list(CGRP_AXIS),
        "MOD|CGRP_core": list(CGRP_CORE),
        "MOD|NGF_TrkA": list(NGF_TRKA),
        "MOD|Sensory_nf": list(SENSORY_NF),
    }
    for g in SINGLE_GENES:
        func[f"GENE|{g}"] = [g]
    sizes = []
    for s, rep in [(2, 60), (3, 60), (4, 60), (6, 60), (9, 60)]:
        sizes.extend([s] * rep)
    rsets = make_random_sets(pool, sizes, 1, seed=21)
    return {"func": func, "random": rsets}


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


# ------------------------------------------------------------------ Part B
def tls_function_sample(path: Path) -> dict | None:
    ad = sc.read_h5ad(path)
    obs = ad.obs
    dist = obs["distance_to_nearest_nerve"].to_numpy(dtype=float)
    keep = np.isfinite(dist)
    X = ad.X
    X = np.asarray(X.toarray() if hasattr(X, "toarray") else X)[keep]
    genes_all = np.asarray(ad.var_names)
    dist = dist[keep]
    tc = np.asarray(obs["total_counts"].to_numpy(dtype=float))[keep]
    pm = np.asarray(obs["pct_mt"].to_numpy(dtype=float))[keep]
    del ad

    from run_niche_spatial import MIN_SPOT_FRAC
    expressed = (X > 0).mean(axis=0)
    sel = expressed >= MIN_SPOT_FRAC
    X = X[:, sel].astype(np.float64)
    genes = genes_all[sel]

    lib = tc.copy(); lib[lib <= 0] = 1.0
    Xl = np.log1p(X / lib[:, None] * 1e4)
    Xl = (Xl - Xl.mean(axis=0)) / (Xl.std(axis=0) + 1e-12)

    gidx = {g: i for i, g in enumerate(genes)}

    def score(glist, min_frac=0.5):
        idx = [gidx[g] for g in glist if g in gidx]
        if len(idx) < max(2, int(np.ceil(len(glist) * min_frac))):
            return None, len(idx)
        return Xl[:, idx].mean(axis=1), len(idx)

    tls, n_tls = score(TLS_EXT)
    if tls is None:
        return None

    def resid(v):
        Z = np.column_stack([np.log(tc + 1.0), pm, np.ones(len(tc))])
        b, *_ = np.linalg.lstsq(Z, v, rcond=None)
        return v - Z @ b

    tls_r = resid(tls)
    hi = tls_r > np.median(tls_r)          # TLS+ spots (同 Stage 17 定义)
    near = (dist <= 150.0) & hi
    far = (dist >= 300.0) & hi
    out = {"sample": path.stem, "n_spot": int(len(dist)),
           "n_tls_genes": n_tls, "n_tls_hi": int(hi.sum()),
           "n_tls_near": int(near.sum()), "n_tls_far": int(far.sum())}

    panels = {"cyt": CYT, "exh": EXH, "tgfb": TGFB}
    scores = {}
    for k, glist in panels.items():
        v, n = score(glist)
        if v is not None:
            scores[k] = resid(v)
            out[f"n_{k}_genes"] = n
    if "cyt" in scores and "exh" in scores:
        scores["cyt_minus_exh"] = scores["cyt"] - scores["exh"]

    for k, v in scores.items():
        out[f"{k}_near"] = float(v[near].mean())
        out[f"{k}_far"] = float(v[far].mean())
        out[f"{k}_delta"] = float(v[near].mean() - v[far].mean())
        if near.sum() >= 10 and far.sum() >= 10:
            out[f"{k}_p_mwu"] = float(stats.mannwhitneyu(v[near], v[far]).pvalue)
    return out


def main(nsim: int = 199):
    files = sorted(DIR_H5.glob("*.h5ad"))
    pool = sorted(union_expressed_genes(files))
    blocks = build_blocks(pool)
    print(f"[init] pool={len(pool)} readouts={len(blocks['func'])} random={len(blocks['random'])}", flush=True)

    # ---------------- Part A: 轴梯度 ----------------
    rows = []
    t0 = time.time()
    for i, p in enumerate(files, 1):
        try:
            o, _ = process_sample(p, nsim, 20261111, 40, blocks_global=blocks)
        except Exception as e:
            print(f"  [skip] {p.stem}: {type(e).__name__}: {e}")
            continue
        rows.append(o)
        print(f"  [{i}/{len(files)}] {p.stem}: K={len(o)} {time.time()-t0:.0f}s", flush=True)
    P = pd.concat(rows, ignore_index=True)
    P.to_csv(TAB / "stage21_axis_gradient_per_sample.csv.gz", index=False, compression="gzip")

    meta_rows = []
    members = {k: len(v) for k, v in blocks["func"].items()}
    for name, g in P[P.kind != "random"].groupby("name", sort=True):
        n = g["rho"].notna().sum()
        if n < 10:
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

    R = P[P.kind == "random"]
    R = R.groupby("name").agg(size=("rho", "size"), mean_rho=("rho", "mean"),
                              stouffer_Z=("u_msr", lambda u: float(np.sum(u) / np.sqrt(len(u)))))
    R = R.rename(columns={"size": "size_"}).reset_index()
    R["size"] = pd.to_numeric(R["size_"], errors="coerce")
    zr_all = []
    for s, g in R.groupby("size"):
        br = np.abs(g["mean_rho"].to_numpy(float))
        zr_all.extend(((np.abs(g["mean_rho"]) - br.mean()) / (br.std(ddof=1) + 1e-12)).tolist())
    zr_all = np.sort(np.asarray(zr_all))
    buckets = {int(s): g for s, g in R.groupby("size")}
    sizes = np.array(sorted(buckets))

    def nearest(s):
        return int(sizes[np.argmin(np.abs(np.log(sizes) - np.log(max(1.0, s))))])

    cal = []
    for r in M.itertuples():
        b = nearest(float(r.size))
        g = buckets[b]
        br = np.abs(g["mean_rho"].to_numpy(float))
        ar = abs(float(r.rho_mean))
        zr = (ar - br.mean()) / (br.std(ddof=1) + 1e-12)
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
    M["fdr_emp_rho"] = bh(M["p_emp_rho"].to_numpy())
    M = M.sort_values("rho_mean")
    M.to_csv(TAB / "stage21_axis_gradient_meta.csv", index=False)
    print(M[["name", "size", "n_samples", "rho_mean", "rho_lo", "rho_hi",
             "p_emp_rho", "baseline_rho"]].to_string(index=False))

    # ---------------- Part B: TLS 近/远功能比 ----------------
    crows = []
    for p in files:
        try:
            r = tls_function_sample(p)
        except Exception as e:
            print(f"  [tls-skip] {p.stem}: {type(e).__name__}: {e}")
            continue
        if r:
            crows.append(r)
    CP = pd.DataFrame(crows)
    CP.to_csv(TAB / "stage21_tls_function_per_sample.csv", index=False)

    mrows = []
    for panel in ["cyt", "exh", "tgfb", "cyt_minus_exh"]:
        col = f"{panel}_delta"
        if col not in CP.columns:
            continue
        v = CP[col].dropna()
        if len(v) < 10:
            continue
        t, pv = stats.ttest_1samp(v, 0.0)
        mrows.append(dict(
            panel=panel, n=int(len(v)), mean=float(v.mean()),
            lo=float(v.mean() - 1.96 * v.std(ddof=1) / np.sqrt(len(v))),
            hi=float(v.mean() + 1.96 * v.std(ddof=1) / np.sqrt(len(v))),
            p=pv,
            n_med_near=int(np.median(CP["n_tls_near"])),
            n_med_far=int(np.median(CP["n_tls_far"])),
        ))
    CM = pd.DataFrame(mrows)
    CM.to_csv(TAB / "stage21_tls_function_meta.csv", index=False)
    print("\n=== TLS near-vs-far function ===")
    print(CM.to_string(index=False) if not CM.empty else "no measurable TLS function metric")


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--nsim", type=int, default=199)
    a = ap.parse_args()
    main(a.nsim)

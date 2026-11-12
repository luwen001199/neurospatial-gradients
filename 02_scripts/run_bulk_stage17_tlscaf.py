"""Stage 17 (空间层): TLS-CAF 共编程 / 胶原交联成熟度 / 缺氧-纤维化耦合.

用户第三轮方案的三条空间表型线, 全部走原稿三重校准管线
(秩偏 Spearman + Moran 谱随机化 MSR + size-matched 随机集基线).

Part A 梯度读数 (process_sample 复用):
  TLS 核心与扩展集            → 近神经富集已知(Z=-2.2), 此处给出基因集级数字
  CAF-TLS 通讯输出 (CXCL12/FAP/POSTN)
  胶原交联成熟度 (LOX/LOXL/PLOD/COL/FBN1/EMILIN1)
  缺氧 (CA9/VEGFA/SLC2A1/ANGPTL4/...)

Part B TLS-CAF 耦合检验:
  近神经 TLS 位点是否携带更强 CAF 输出 ——
  (i)  TLS+ spots 内 partial Spearman(CAFout, dist), 预期负 = 近神经 TLS 的 CAF 输出更高
  (ii) r(TLS, CAFout) 近带 vs 远带, Fisher-z 配对差

输出: 03_results/tables/spatial_tlscaf/
  stage17_gradient_meta.csv      梯度读数 meta (MSR + 随机基线双校准)
  stage17_gradient_per_sample.csv.gz
  stage17_coupling_per_sample.csv
  stage17_coupling_meta.csv
"""
from __future__ import annotations

import importlib.util
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
PROC = CFG.DIR_PROC
TAB = Path(r"D:/嗜神经浸润/03_results/tables/spatial_tlscaf")
TAB.mkdir(parents=True, exist_ok=True)

_spec = importlib.util.spec_from_file_location("npp_split", Path(__file__).parent / "run_npp_split.py")
_npp = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_npp)
PERINEURAL_STROMA = _npp.MODULES["Perineural_stroma"]

# ------------------------------------------------------------------ 读数定义
TLS_CORE = ["CCL21", "CCL19", "CXCL13", "LAMP3", "MS4A1", "CD79A"]
TLS_EXT = TLS_CORE + ["MEF2B", "POU2AF1", "BANK1"]
CAF_TLS_OUT = ["CXCL12", "FAP", "POSTN"]                     # CAF→TLS 通讯输出
COLL_MAT = ["LOX", "LOXL1", "LOXL2", "PLOD1", "PLOD2",
            "COL1A1", "COL3A1", "FBN1", "EMILIN1"]
LOX_ONLY = ["LOX", "LOXL1", "LOXL2", "PLOD1", "PLOD2"]
FIBRIL_ORG = ["FBN1", "EMILIN1", "EMILIN2", "LTBP1", "LTBP2"]
HYPOXIA = ["CA9", "VEGFA", "SLC2A1", "ANGPTL4", "ADM", "NDRG1",
           "BNIP3", "BNIP3L", "DDIT3", "PDK1"]
HYPOXIA_CORE = ["CA9", "VEGFA", "SLC2A1", "ANGPTL4"]
SINGLE_GENES = ["CXCL12", "FAP", "POSTN", "PLOD2", "EMILIN1", "FBN1",
                "CA9", "VEGFA", "SLC2A1", "CXCL13", "CCL21", "LOX", "COL3A1"]


def build_blocks(pool: list[str]) -> dict:
    func = {
        "MOD|TLS_core": list(TLS_CORE),
        "MOD|TLS_ext": list(TLS_EXT),
        "MOD|CAF_TLS_out": list(CAF_TLS_OUT),
        "MOD|Coll_mat": list(COLL_MAT),
        "MOD|LOX_only": list(LOX_ONLY),
        "MOD|Fibril_org": list(FIBRIL_ORG),
        "MOD|Hypoxia": list(HYPOXIA),
        "MOD|Hypoxia_core": list(HYPOXIA_CORE),
    }
    for g in SINGLE_GENES:
        func[f"GENE|{g}"] = [g]
    sizes = []
    for s, rep in [(3, 50), (4, 50), (6, 50), (9, 50), (10, 50), (12, 50)]:
        sizes.extend([s] * rep)
    rsets = make_random_sets(pool, sizes, 1, seed=17)
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
def _zs_genes(X_raw: np.ndarray, genes: np.ndarray, glist: list[str]):
    idx = [i for i, g in enumerate(genes) if g in glist]
    if len(idx) < max(2, len(glist) // 2):
        return None, len(idx)
    sub = X_raw[:, idx].astype(np.float64)
    return sub, len(idx)


def coupling_sample(path: Path) -> dict | None:
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

    expressed = (X > 0).mean(axis=0)
    from run_niche_spatial import MIN_SPOT_FRAC
    sel = expressed >= MIN_SPOT_FRAC
    X = X[:, sel].astype(np.float64)
    genes = genes_all[sel]

    lib = tc.copy(); lib[lib <= 0] = 1.0
    Xl = np.log1p(X / lib[:, None] * 1e4)
    Xl = (Xl - Xl.mean(axis=0)) / (Xl.std(axis=0) + 1e-12)

    gidx = {g: i for i, g in enumerate(genes)}

    def score(glist):
        idx = [gidx[g] for g in glist if g in gidx]
        if len(idx) < max(2, len(glist) // 2):
            return None, len(idx)
        return Xl[:, idx].mean(axis=1), len(idx)

    tls, n_tls = score(TLS_EXT)
    caf, n_caf = score(CAF_TLS_OUT)
    if tls is None or caf is None:
        return None

    # 控制变量: logUMI + pct_mt → 残差化
    def resid(v):
        Z = np.column_stack([np.log(tc + 1.0), pm, np.ones(len(tc))])
        b, *_ = np.linalg.lstsq(Z, v, rcond=None)
        return v - Z @ b

    tls_r, caf_r = resid(tls), resid(caf)
    out = {"sample": path.stem, "n_spot": int(len(dist)),
           "n_tls_genes": n_tls, "n_caf_genes": n_caf,
           "frac_tls_pos": float((tls > 0).mean()),
           "mean_tls": float(tls.mean()), "mean_cafout": float(caf.mean())}

    # (i) TLS+ spots 内 partial Spearman(CAFout, dist)
    hi = tls_r > np.median(tls_r)
    out["n_tls_hi"] = int(hi.sum())
    if hi.sum() >= 60:
        rho, p = stats.spearmanr(caf_r[hi], dist[hi])
        out["rho_cafout_dist_inTLS"] = float(rho)
        out["p_cafout_dist_inTLS"] = float(p)

    # TLS+ 且近/远的 CAFout 水平差 (配对样本内比较)
    near = dist <= 150.0
    far = dist >= 300.0
    if (hi & near).sum() >= 15 and (hi & far).sum() >= 15:
        out["cafout_tls_near"] = float(caf_r[hi & near].mean())
        out["cafout_tls_far"] = float(caf_r[hi & far].mean())
        u = stats.mannwhitneyu(caf_r[hi & near], caf_r[hi & far]).pvalue
        out["p_cafout_nearVfar_inTLS"] = float(u)

    # (ii) r(TLS, CAFout) 近带 vs 远带
    bands = {"near200": dist <= 200.0, "far200": dist > 200.0}
    for k, m in bands.items():
        if m.sum() >= 60:
            r = stats.spearmanr(tls_r[m], caf_r[m]).statistic
            out[f"r_tlscaf_{k}"] = float(r)
            out[f"n_{k}"] = int(m.sum())
    if "r_tlscaf_near200" in out and "r_tlscaf_far200" in out:
        z1 = np.arctanh(np.clip(out["r_tlscaf_near200"], -0.999, 0.999))
        z2 = np.arctanh(np.clip(out["r_tlscaf_far200"], -0.999, 0.999))
        n1, n2 = out["n_near200"], out["n_far200"]
        out["zdiff_nearVfar"] = (z1 - z2) / np.sqrt(1 / (n1 - 3) + 1 / (n2 - 3))
    return out


def main(nsim: int = 199):
    files = sorted(DIR_H5.glob("*.h5ad"))
    pool = sorted(union_expressed_genes(files))
    blocks = build_blocks(pool)
    print(f"[init] pool={len(pool)} readouts={len(blocks['func'])} random={len(blocks['random'])}", flush=True)

    rows, rings = [], []
    t0 = time.time()
    for i, p in enumerate(files, 1):
        try:
            o, rg = process_sample(p, nsim, 20261107, 40, blocks_global=blocks)
        except Exception as e:
            print(f"  [skip] {p.stem}: {type(e).__name__}: {e}")
            continue
        rows.append(o)
        if len(rg):
            rings.append(rg)
        print(f"  [{i}/{len(files)}] {p.stem}: K={len(o)} {time.time()-t0:.0f}s", flush=True)
    P = pd.concat(rows, ignore_index=True)
    P.to_csv(TAB / "stage17_gradient_per_sample.csv.gz", index=False, compression="gzip")
    if rings:
        pd.concat(rings, ignore_index=True).to_csv(TAB / "stage17_ring_profiles.csv", index=False)

    # ---------------- meta + 双重校准 ----------------
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
    M.to_csv(TAB / "stage17_gradient_meta.csv", index=False)
    print(M[["name", "size", "n_samples", "rho_mean", "rho_lo", "rho_hi",
             "p_emp_rho", "baseline_rho"]].to_string(index=False))

    # ---------------- Part B: TLS-CAF 耦合 ----------------
    crows = []
    for p in files:
        try:
            r = coupling_sample(p)
        except Exception as e:
            print(f"  [coup-skip] {p.stem}: {type(e).__name__}: {e}")
            continue
        if r:
            crows.append(r)
    CP = pd.DataFrame(crows)
    CP.to_csv(TAB / "stage17_coupling_per_sample.csv", index=False)

    mrows = []
    for col, label in [("rho_cafout_dist_inTLS", "CAFout~dist within TLS+ spots"),
                       ("cafout_tls_near", "CAFout level near-nerve TLS+"),
                       ("cafout_tls_far", "CAFout level far TLS+"),
                       ("zdiff_nearVfar", "Fisher-z(TLS-CAF r near vs far)")]:
        if col not in CP.columns:
            continue
        v = CP[col].dropna()
        if len(v) < 10:
            continue
        t, p = stats.ttest_1samp(v, 0.0)
        mrows.append(dict(metric=label, n=int(len(v)), mean=float(v.mean()),
                          lo=float(v.mean() - 1.96 * v.std(ddof=1) / np.sqrt(len(v))),
                          hi=float(v.mean() + 1.96 * v.std(ddof=1) / np.sqrt(len(v))),
                          p=float(p)))
    CM = pd.DataFrame(mrows)
    CM.to_csv(TAB / "stage17_coupling_meta.csv", index=False)
    print("\n=== TLS-CAF coupling ===")
    print(CM.to_string(index=False) if not CM.empty else "no measurable coupling metric")


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--nsim", type=int, default=199)
    a = ap.parse_args()
    main(a.nsim)

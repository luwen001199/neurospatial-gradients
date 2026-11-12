"""神经周围微环境的空间功能读数 —— 配体受体 / 通路 / 细胞类型 的距离梯度,
全部用 Moran 谱随机化 (MSR) 零分布校准。

复用 run_msr_spatial 的统计量定义: 秩空间偏 Spearman(剔除 rank(log UMI)、
rank(pct_mt)), 箱内置换系数域解析计算。区别在于被随机化的对象从"单基因场"
换成"功能读数场"(细胞类型丰度 / LR 共表达 / 通路活性), K ~ 400 而非 10419,
成本极低。

输出
----
01_data/processed/niche_per_sample.csv.gz   每样本每读数: rho, MSR 零分布均值/SD, u
01_data/processed/niche_ring_profiles.csv   每样本每读数: 各距离环的平均评分(用于画梯度曲线)
01_data/processed/niche_scores_meta.json    读数元数据(名称/类别/成员覆盖)
03_results/tables/niche_meta_summary.csv    跨样本 meta (Stouffer + 跨样本 t)
03_results/tables/niche_msr_*.csv           各类别的分表

用法
----
python run_niche_spatial.py scan  [--nsim 499] [--bin 40] [--only GSM...]
python run_niche_spatial.py meta
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
from scipy import stats

sys.path.insert(0, r"D:/嗜神经浸润/00_pipeline")
sys.path.insert(0, r"D:/嗜神经浸润/02_scripts")

from nsgrad import config as CFG                      # noqa: E402
from run_msr_spatial import (                         # noqa: E402
    _rank_cols, build_knn_graph, moran_spectrum, morans_i_from_coef,
    ortho_cov, MIN_SPOT_FRAC,
)
import niche_lib as NL                               # noqa: E402

DIR_H5 = CFG.DIR_PROC / "spot_level"
PROC = CFG.DIR_PROC
TAB = CFG.DIR_TAB
FDR = 0.05
RING_ORDER = ["r1_0_50", "r2_50_100", "r3_100_200", "r4_200_400", "out_gt400"]


# ------------------------------------------------------------------ 读数构建
def make_random_sets(gene_pool: list[str], sizes: list[int], n_per_size: int,
                     seed: int = 99):
    """大小匹配的随机基因集 —— 回答"该读数是否强于任何同等大小的基因集"。

    近神经区域在组织学上本就是不同的区域(基质/血管/边缘), 距离场 Moran's I
    ≈ 0.97 近乎确定性梯度, 因此**任何**带大尺度空间结构的基因集都会与之相关。
    随机基因集给出这条基线: 只有超出该基线的读数才能称为"神经周围特异"。
    """
    rng = np.random.default_rng(seed)
    pool = np.asarray(gene_pool)
    out = {}
    k = 0
    for s in sizes:
        for _ in range(n_per_size):
            g = rng.choice(pool, size=min(s, len(pool)), replace=False)
            out[f"RAND|s{s}|{k}"] = list(g)
            k += 1
    return out


def build_blocks(genes_avail: set[str] | None = None, random_sets=None):
    """细胞类型 + LR 对 + 通路 三类读数。

    genes_avail=None 时不过滤成员(仅做互斥去重), 由各样本自身的检出情况决定
    某个读数是否可算 —— 这样免疫细胞/低表达配体这类稀疏标记不会被"跨样本
    ≥10% spots"的阈值提前抹掉; 其可测性由 meta 阶段的 n_samples 如实反映。
    """
    ct_names, ct = NL.build_celltype_sets(genes_avail)
    lr = {}
    for lig, rec, axis, note in NL.LR_PAIRS:
        nm = f"LR|{lig}-{rec}|{axis}"
        lr[nm] = [lig, rec]
    pw_raw = NL.load_msigdb(collections=("HALLMARK", "REACTOME", "KEGG"))
    pw = pw_raw if genes_avail is None else {
        k: v for k, v in pw_raw.items()
        if sum(g in genes_avail for g in v) / len(v) >= 0.5}
    blocks = {"celltype": ct, "lr": lr, "pathway": pw}
    if random_sets:
        blocks["random"] = random_sets
    return blocks


# ---------------------------------------------------------------- per-sample
def process_sample(path: Path, nsim: int, seed: int, bin_target: int,
                   blocks_global: dict | None = None):
    ad = sc.read_h5ad(path)
    obs = ad.obs
    dist = obs["distance_to_nearest_nerve"].to_numpy(dtype=float)
    ring = obs["ring_id"].astype(str).to_numpy() if "ring_id" in obs else None
    keep = np.isfinite(dist)
    X_raw = ad.X
    X_raw = np.asarray(X_raw.toarray() if hasattr(X_raw, "toarray") else X_raw)[keep]
    xy = np.column_stack([
        np.asarray(obs["x_um"].to_numpy(dtype=float)),
        np.asarray(obs["y_um"].to_numpy(dtype=float)),
    ])[keep]
    d = dist[keep]
    ring = ring[keep] if ring is not None else None
    tc = np.asarray(obs["total_counts"].to_numpy(dtype=float))[keep]
    pm = np.asarray(obs["pct_mt"].to_numpy(dtype=float))[keep]
    genes_all = np.asarray(ad.var_names)
    del ad

    expressed = (X_raw > 0).mean(axis=0)
    sel = expressed >= MIN_SPOT_FRAC
    Xs = X_raw[:, sel].astype(np.float64)
    genes = genes_all[sel]

    lib = tc.copy()
    lib[lib <= 0] = 1.0
    Xl = np.log1p(Xs / lib[:, None] * 1e4)
    del Xs
    Xl = (Xl - Xl.mean(axis=0)) / (Xl.std(axis=0) + 1e-12)

    # ---- 读数矩阵 ----
    blocks = blocks_global if blocks_global is not None else \
        build_blocks(set(genes.tolist()))
    S, meta, W = NL.build_score_matrix(
        Xl, genes, blocks,
        min_frac_by_kind={"celltype": 0.5, "lr": 1.0, "pathway": 0.5,
                          "random": 0.5})
    del Xl, W
    gname = path.stem
    K = S.shape[1]

    # ---- 秩空间 ----
    R = _rank_cols(S)
    R -= R.mean(axis=0)
    R /= (R.std(axis=0) + 1e-12)
    rd = stats.rankdata(d)
    rd = (rd - rd.mean()) / (rd.std() + 1e-12)
    Q = ortho_cov(np.log1p(tc), pm)

    v_res = rd.copy()
    for j in range(Q.shape[1]):
        v_res -= Q[:, j] * (Q[:, j] @ v_res)
    vv = float(v_res @ v_res)
    rtQ = R.T @ Q
    den_obs = np.sqrt(np.maximum(R.shape[0] - (rtQ ** 2).sum(axis=1), 1e-12) * vv)
    rho_obs = np.clip((R.T @ v_res) / den_obs, -1.0, 1.0)

    # ---- Moran 谱 ----
    B = build_knn_graph(xy)
    lam, V, S0 = moran_spectrum(B)
    del B
    order = np.argsort(lam)
    lam_s = lam[order]
    V_s = np.ascontiguousarray(V[:, order])
    del V, lam
    n = R.shape[0]
    C = V_s.T @ R
    w = V_s.T @ v_res
    U = V_s.T @ Q
    moran_dist = float(morans_i_from_coef((V_s.T @ rd)[:, None], lam_s, S0, n)[0])
    moran_obs = morans_i_from_coef(C, lam_s, S0, n)
    del V_s

    nb = int(max(4, min(120, n // bin_target)))
    edges = np.linspace(0, n, nb + 1).astype(int)
    rng = np.random.default_rng(seed)
    rho_msr = np.empty((nsim, K), dtype=np.float64)
    moran_sur = np.empty((nsim, K), dtype=np.float64)
    for r in range(nsim):
        Cs = np.empty_like(C)
        for b in range(nb):
            s, e = edges[b], edges[b + 1]
            if e - s > 1:
                Cs[s:e] = C[s:e][rng.permutation(e - s)]
            else:
                Cs[s:e] = C[s:e]
        moran_sur[r] = morans_i_from_coef(Cs, lam_s, S0, n)
        num = w @ Cs
        Cu = Cs.T @ U
        den = np.sqrt(np.maximum((Cs * Cs).sum(axis=0) - (Cu ** 2).sum(axis=1), 1e-12) * vv)
        rho_msr[r] = np.clip(num / den, -1.0, 1.0)
    del C, Cs, Cu

    # ---- 环剖面(画梯度曲线用, 原始 z 空间) ----
    ring_rows = []
    if ring is not None:
        for rk in RING_ORDER:
            m = ring == rk
            if m.sum() >= 5:
                ring_rows.append(dict(sample=gname, ring=rk, n_spot=int(m.sum()),
                                      **{f"s{i}": float(S[m, i].mean()) for i in range(K)}))

    msr_mu = rho_msr.mean(axis=0)
    msr_sd = rho_msr.std(axis=0, ddof=1)
    sd_shrunk = np.maximum(msr_sd, 1e-4)          # 避免除零
    u = (rho_obs - msr_mu) / sd_shrunk

    out = pd.DataFrame({
        "sample": gname,
        "name": [m["name"] for m in meta],
        "kind": [m["kind"] for m in meta],
        "rho": rho_obs,
        "u_msr": u,
        "mu_sur": msr_mu,
        "sd_sur": msr_sd,
        "moran_obs": moran_obs,
        "moran_sur": moran_sur.mean(axis=0),
        "moran_dist": moran_dist,
        "n_spot": n,
    })
    ring_df = pd.DataFrame(ring_rows)
    if len(ring_df):
        ring_df = ring_df.rename(columns={f"s{i}": meta[i]["name"] for i in range(K)})
    return out, ring_df


# --------------------------------------------------------------------- scan
def union_expressed_genes(files, limit: int = 12):
    """读数成员的可用基因集。

    优先复用 Step3a/MSR 已确定的表达基因表(10,419 个, 与 NPP 完全相同的
    过滤标准), 保证细胞类型/LR/通路读数与基因级分析在同一基因宇宙内可比。
    """
    for cand in (PROC / "msr_gene_level.fine.csv", PROC / "msr_gene_level.csv"):
        if cand.exists():
            return set(pd.read_csv(cand, index_col=0).index.astype(str))
    import scanpy as _sc
    gs: set[str] = set()
    for p in files[:3]:
        try:
            ad = _sc.read_h5ad(p)
            X = np.asarray(ad.X.toarray())
            sel = (X > 0).mean(axis=0) >= MIN_SPOT_FRAC
            gs.update(np.asarray(ad.var_names)[sel].tolist())
            del ad, X
        except Exception:
            continue
    return gs


def stage_scan(nsim: int, seed: int, bin_target: int, only: str | None,
               mode: str = "main"):
    files = sorted(DIR_H5.glob("*.h5ad"))
    if only:
        files = [f for f in files if only in f.stem]
    # 先用各样本表达基因的交集确定统一的读数定义, 保证跨样本列名与成员一致
    if mode == "random":
        pool = sorted(union_expressed_genes(files))
        # 与真实读数的大小分布对齐
        sizes = []
        for s in (5, 10, 15, 20, 25, 30, 40, 50, 65, 80, 100, 130, 160, 200,
                  260, 330, 400):
            sizes.extend([s] * 30)
        rsets = make_random_sets(pool, sizes, 1, seed=seed)
        blocks_global = {"random": rsets}          # 只跑随机集, 节省算力
        print(f"  [init] random sets {len(rsets)} (pool {len(pool)})")
    else:
        blocks_global = build_blocks(None)
        print(f"  [init] celltype {len(blocks_global['celltype'])}, "
              f"lr {len(blocks_global['lr'])}, "
              f"pathway {len(blocks_global['pathway'])}")
    rows, rings = [], []
    t0 = time.time()
    for i, p in enumerate(files, 1):
        try:
            res = process_sample(p, nsim, seed, bin_target, blocks_global)
        except Exception as e:
            print(f"  [skip] {p.stem}: {type(e).__name__}: {e}")
            continue
        if res is None:
            continue
        o, rg = res
        rows.append(o)
        if len(rg):
            rings.append(rg)
        print(f"  [{i}/{len(files)}] {p.stem}: K={len(o)} n={o['n_spot'].iloc[0]} "
              f"{time.time()-t0:.0f}s", flush=True)
    P = pd.concat(rows, ignore_index=True)
    suffix = "" if mode == "main" else f".{mode}"
    P.to_csv(PROC / f"niche_per_sample{suffix}.csv.gz", index=False,
             compression="gzip")
    if rings and mode == "main":
        Rg = pd.concat(rings, ignore_index=True)
        Rg.to_csv(PROC / "niche_ring_profiles.csv", index=False)
    if mode == "main":
        meta = P[["name", "kind"]].drop_duplicates()
        meta.to_json(PROC / "niche_scores_meta.json", orient="records",
                     force_ascii=False)
    print(f"\n[niche] -> {PROC/f'niche_per_sample{suffix}.csv.gz'} ({len(P)} 行, "
          f"{P['sample'].nunique()} 样本, {P['name'].nunique()} 读数)")
    return P


# --------------------------------------------------------------------- meta
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


def stage_meta():
    P = pd.read_csv(PROC / "niche_per_sample.csv.gz")
    rows = []
    for (name, kind), g in P.groupby(["name", "kind"], sort=False):
        u = g["u_msr"].to_numpy(float)
        rho = g["rho"].to_numpy(float)
        ns = len(u)
        z = np.sum(u) / np.sqrt(ns)                       # Stouffer
        p_z = 2 * stats.norm.sf(abs(z))
        m, sd = u.mean(), u.std(ddof=1)
        t = m / (sd / np.sqrt(ns)) if sd > 0 else np.nan
        p_t = 2 * stats.t.sf(abs(t), df=ns - 1) if np.isfinite(t) else np.nan
        rows.append(dict(
            name=name, kind=kind, n_samples=ns,
            mean_rho=float(np.mean(rho)), median_rho=float(np.median(rho)),
            consistency=float((np.sign(rho) == np.sign(np.mean(rho))).mean()),
            stouffer_Z=float(z), p_msr=p_z,
            mean_u=float(m), t_between=float(t) if np.isfinite(t) else np.nan,
            p_msr_between=p_t,
            moran_obs=float(g["moran_obs"].mean()),
            moran_sur=float(g["moran_sur"].mean()),
            sd_sur=float(g["sd_sur"].mean()),
        ))
    Gm = pd.DataFrame(rows).set_index("name")
    Gm["fdr_msr"] = bh(Gm["p_msr"].to_numpy())
    Gm["fdr_msr_between"] = bh(np.nan_to_num(Gm["p_msr_between"].to_numpy(), nan=1.0))
    Gm = Gm.sort_values("stouffer_Z")
    Gm.to_csv(TAB / "niche_meta_summary.csv", index=True)

    for kind in ["celltype", "lr", "pathway"]:
        sub = Gm[Gm["kind"] == kind]
        if len(sub):
            sub.to_csv(TAB / f"niche_msr_{kind}.csv", index=True)

    print(f"[meta] 读数 {len(Gm)}  (celltype {int((Gm.kind=='celltype').sum())}, "
          f"lr {int((Gm.kind=='lr').sum())}, pathway {int((Gm.kind=='pathway').sum())})")
    for kind in ["celltype", "lr", "pathway"]:
        sub = Gm[Gm["kind"] == kind].copy()
        eff = sub["mean_rho"].abs() >= 0.03
        n_sig = int(((sub["fdr_msr"] < FDR) & eff).sum())
        print(f"\n--- {kind}: MSR FDR<0.05 且 |rho|>=0.03 共 {n_sig} / {len(sub)}")
        sig = sub[(sub["fdr_msr"] < FDR) & eff].sort_values("stouffer_Z")
        print(sig.head(15)[["n_samples", "mean_rho", "consistency", "stouffer_Z",
                            "fdr_msr", "fdr_msr_between"]].round(4).to_string())
        print("  反向 top5:")
        print(sig.tail(5)[["n_samples", "mean_rho", "consistency", "stouffer_Z",
                           "fdr_msr"]].round(4).to_string())
    return Gm


def stage_meta_random():
    """大小匹配随机基因集的零分布 —— 回答"多少随机集也能通过同一门槛"。

    距离场的 Moran's I ≈ 0.97, 近乎确定性大尺度梯度, 因此任何带空间结构的
    基因集都会与之相关。只有显著超出随机集基线的读数才能称为"神经周围特异"。
    """
    P = pd.read_csv(PROC / "niche_per_sample.random.csv.gz")
    if (P["kind"] != "random").all():
        raise SystemExit("[meta_random] 随机集结果为空 —— 请先跑 "
                         "`python run_niche_spatial.py scan --mode random`")
    rows = []
    for name, g in P.groupby("name"):
        u = g["u_msr"].to_numpy(float)
        rho = g["rho"].to_numpy(float)
        ns = len(u)
        z = float(u.sum() / np.sqrt(ns))
        size = int(str(name).split("|")[1][1:]) if "|" in str(name) else np.nan
        rows.append(dict(name=name, size=size, n_samples=ns,
                         mean_rho=float(rho.mean()), abs_rho=float(abs(rho.mean())),
                         stouffer_Z=z, p_msr=float(2 * stats.norm.sf(abs(z))),
                         mean_u=float(u.mean()), sd_u=float(u.std(ddof=1))))
    R = pd.DataFrame(rows)
    R.to_csv(TAB / "niche_random_null.csv", index=False)

    print(f"[meta_random] 随机集 {len(R)}  (size {R['size'].min()}~{R['size'].max()})")
    print("\n--- 按基因集大小的零分布 ---")
    print(R.groupby("size").agg(
        n_sets=("mean_rho", "size"),
        mean_abs_rho=("abs_rho", "mean"),
        sd_u=("sd_u", "mean")).round(4).to_string())
    q = [50, 75, 90, 95, 99]
    print("\n--- 随机集 |mean_rho| 分位 ---")
    print(dict(zip(q, np.round(np.percentile(R.abs_rho, q), 4))))
    print("--- 随机集 |Stouffer Z| 分位 ---")
    print(dict(zip(q, np.round(np.percentile(R.stouffer_Z.abs(), q), 3))))
    print(f"\n随机集通过率: p<0.05 = {(R.p_msr<0.05).mean():.3f}; "
          f"|rho|>=0.03 = {(R.abs_rho>=0.03).mean():.3f}; "
          f"两者同时 = {((R.p_msr<0.05)&(R.abs_rho>=0.03)).mean():.3f}")
    return R


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("stage", choices=["scan", "meta", "meta_random"])
    ap.add_argument("--nsim", type=int, default=499)
    ap.add_argument("--bin", type=int, default=40)
    ap.add_argument("--seed", type=int, default=7)
    ap.add_argument("--only", type=str, default=None)
    ap.add_argument("--mode", choices=["main", "random"], default="main")
    a = ap.parse_args()
    if a.stage == "scan":
        stage_scan(a.nsim, a.seed, a.bin, a.only, a.mode)
    elif a.stage == "meta":
        stage_meta()
    else:
        stage_meta_random()

"""方案五: 以神经 spot 为根的扩散伪时间 (DPT) —— 神经周围微环境的空间演化轴。

思路
----
Visium spot 的"状态"是连续的。若神经周围微环境存在有序的空间演化
(神经信号 -> CAF 活化 -> ECM 沉积 -> 肿瘤增殖), 那么在**表达空间**(而非
物理空间)上做扩散图、并指定神经 spot 为根, 得到的伪时间应当:
  (1) 与"到最近神经的距离"显著正相关 —— 轨迹与空间方向一致;
  (2) 各功能读数沿伪时间呈有序的先后峰值。

注意: 伪时间由表达构建, 距离完全不参与图的构建, 因此 (1) 是一个**非平凡**的
一致性检验 —— 若两者无关, 说明不存在统一的空间演化轴。

实现要点
--------
全程保持稀疏: normalize_total + log1p 在稀疏矩阵上完成, HVG 过滤后才子集化,
只对读数实际用到的基因稠密化, 避免 1376 x 36591 稠密阵导致的内存溢出。

输出
----
01_data/processed/pseudotime_per_sample.csv.gz  每样本: 伪时间-距离相关等
01_data/processed/pseudotime_profiles.csv       每样本 x 伪时间分箱的读数均值
03_results/tables/pseudotime_summary.csv        跨样本合并
03_results/tables/pseudotime_readout_order.csv  读数峰值顺序
"""
from __future__ import annotations

import gc
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse
from scipy import stats

sys.path.insert(0, r"D:/嗜神经浸润/00_pipeline")
sys.path.insert(0, r"D:/嗜神经浸润/02_scripts")

from nsgrad import config as CFG                      # noqa: E402
import niche_lib as NL                               # noqa: E402

DIR_H5 = CFG.DIR_PROC / "spot_level"
PROC = CFG.DIR_PROC
TAB = CFG.DIR_TAB
N_BINS = 10
N_PCS = 30

FOCUS_PATHWAYS = {
    "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
    "HALLMARK_TGF_BETA_SIGNALING",
    "HALLMARK_ANGIOGENESIS",
    "HALLMARK_HYPOXIA",
    "HALLMARK_G2M_CHECKPOINT",
    "HALLMARK_INTERFERON_GAMMA_RESPONSE",
    "HALLMARK_INFLAMMATORY_RESPONSE",
    "HALLMARK_COMPLEMENT",
    "REACTOME_EXTRACELLULAR_MATRIX_ORGANIZATION",
    "REACTOME_COLLAGEN_FORMATION",
    "REACTOME_AXON_GUIDANCE",
    "KEGG_FOCAL_ADHESION",
    "KEGG_ECM_RECEPTOR_INTERACTION",
}


def _partial_rho(a: np.ndarray, b: np.ndarray, Q: np.ndarray) -> float:
    """秩空间内剔除 Q(已正交化) 后的 Spearman."""
    ra = stats.rankdata(a)
    rb = stats.rankdata(b)
    ra = (ra - ra.mean()) / (ra.std() + 1e-12)
    rb = (rb - rb.mean()) / (rb.std() + 1e-12)
    for j in range(Q.shape[1]):
        ra -= Q[:, j] * (Q[:, j] @ ra)
        rb -= Q[:, j] * (Q[:, j] @ rb)
    return float(ra @ rb / np.sqrt((ra @ ra) * (rb @ rb) + 1e-30))


def _ortho_cov(tc: np.ndarray, pm: np.ndarray) -> np.ndarray:
    lib = np.where(tc <= 0, 1.0, tc)
    Q = np.column_stack([stats.rankdata(np.log1p(lib)), stats.rankdata(pm)])
    Q = (Q - Q.mean(axis=0)) / (Q.std(axis=0) + 1e-12)
    Qn = np.zeros_like(Q)
    for j in range(Q.shape[1]):
        v = Q[:, j].copy()
        for k in range(j):
            v -= Qn[:, k] * (Qn[:, k] @ v)
        nv = np.sqrt(v @ v)
        Qn[:, j] = v / nv if nv > 1e-12 else 0.0
    return Qn


_BLOCKS_CACHE: dict = {}


def get_blocks():
    """读数定义只构建一次(MSigDB 解析 ~1s, 30 个样本就浪费 30s)。"""
    if not _BLOCKS_CACHE:
        _, ct = NL.build_celltype_sets(None)
        pw_all = NL.load_msigdb(collections=("HALLMARK", "REACTOME", "KEGG"))
        pw_sel = {k: v for k, v in pw_all.items() if k in FOCUS_PATHWAYS}
        _BLOCKS_CACHE["blocks"] = {"celltype": ct, "pathway": pw_sel}
        needed: set[str] = set()
        for gl in list(ct.values()) + list(pw_sel.values()):
            needed.update(gl)
        _BLOCKS_CACHE["needed"] = needed
    return _BLOCKS_CACHE["blocks"], _BLOCKS_CACHE["needed"]


def sample_pseudotime(path: Path, seed: int = 5):
    ad = sc.read_h5ad(path)
    obs = ad.obs.copy()
    dist = obs["distance_to_nearest_nerve"].to_numpy(float)
    keep = np.isfinite(dist)
    ad = ad[keep].copy()
    dist = dist[keep]
    n = ad.n_obs
    if n < 200:
        return None

    # ---- 读数定义 (先确定用到的基因, 后面只稠密化这些列) ----
    blocks, needed = get_blocks()

    genes_all = np.asarray(ad.var_names)
    cols = np.array([i for i, g in enumerate(genes_all) if g in needed], dtype=int)

    # ---- 全程稀疏: normalize -> log1p ----
    if not sparse.issparse(ad.X):
        ad.X = sparse.csr_matrix(np.asarray(ad.X))
    sc.pp.normalize_total(ad, target_sum=1e4)
    sc.pp.log1p(ad)

    # ---- 读数矩阵: 只稠密化读数基因列 ----
    Xs = np.asarray(ad.X[:, cols].todense())
    sd = Xs.std(axis=0)
    m = sd > 1e-9
    Xs = Xs[:, m]
    gsub = genes_all[cols][m]
    Xs = (Xs - Xs.mean(axis=0)) / (Xs.std(axis=0) + 1e-12)
    S, meta, _ = NL.build_score_matrix(
        Xs, gsub, blocks,
        min_frac_by_kind={"celltype": 0.4, "pathway": 0.5})
    del Xs

    # ---- HVG -> PCA -> kNN -> diffmap -> DPT ----
    try:
        sc.pp.highly_variable_genes(ad, n_top_genes=2000, flavor="seurat")
        ad = ad[:, ad.var.highly_variable].copy()
    except Exception:
        pass
    sc.pp.pca(ad, n_comps=min(N_PCS, ad.n_vars - 1, n - 1),
              svd_solver="arpack", random_state=seed)
    # 注意: 默认的 umap/pynndescent 后端依赖 numba, 在 Windows 上会以
    # "LLVM ERROR: out of memory" 硬崩(与数据量无关)。改用 gauss 后端
    # (纯 sklearn + scipy, 确定性), 已在真实数据上验证。
    sc.pp.neighbors(ad, n_neighbors=15, method="gauss", random_state=seed)
    sc.tl.diffmap(ad, n_comps=min(15, ad.n_vars - 1))

    nerve = obs.loc[keep, "nerve_spatial"].to_numpy(bool) \
        if "nerve_spatial" in obs.columns else np.zeros(n, bool)
    nerve_score = obs.loc[keep, "nerve_score"].to_numpy(float) \
        if "nerve_score" in obs.columns else np.zeros(n)
    if nerve.sum() >= 1:
        cand = np.where(nerve)[0]
        iroot = int(cand[np.argmax(nerve_score[cand])])
    else:
        iroot = int(np.argmin(dist))
    ad.uns["iroot"] = iroot
    sc.tl.dpt(ad)
    dpt = ad.obs["dpt_pseudotime"].to_numpy(float)
    ok = np.isfinite(dpt)
    if ok.sum() < 100:
        return None

    tc = np.asarray(ad.obs["total_counts"].to_numpy(float))
    pm = np.asarray(ad.obs["pct_mt"].to_numpy(float))
    Qn = _ortho_cov(tc, pm)

    rho, p = stats.spearmanr(dpt[ok], dist[ok])
    rho_part = _partial_rho(dpt[ok], dist[ok], Qn[ok])

    # ---- 沿伪时间分箱的读数剖面 ----
    bins = pd.qcut(dpt[ok], N_BINS, labels=False, duplicates="drop")
    prof = []
    for b in range(N_BINS):
        m = np.zeros(n, bool)
        m[ok] = bins == b
        if m.sum() < 5:
            continue
        row = dict(sample=path.stem, bin=b, n_spot=int(m.sum()),
                   dpt_median=float(np.median(dpt[m])),
                   dist_median=float(np.median(dist[m])))
        for i, mm in enumerate(meta):
            row[mm["name"]] = float(S[m, i].mean())
        prof.append(row)

    # ---- 根 spot 是否落在近神经环内(可解释性) ----
    root_ring = str(obs.loc[keep, "ring_id"].iloc[iroot]) \
        if "ring_id" in obs.columns else ""

    del ad, S
    gc.collect()

    return dict(sample=path.stem, n_spot=int(ok.sum()),
                iroot_is_nerve=bool(nerve.sum() >= 1),
                root_ring=root_ring,
                rho_dpt_dist=float(rho), p_dpt_dist=float(p),
                rho_partial=float(rho_part)), \
        pd.DataFrame(prof)


def main():
    files = sorted(DIR_H5.glob("*.h5ad"))
    rows, profs = [], []
    t0 = time.time()
    for i, p in enumerate(files, 1):
        try:
            out = sample_pseudotime(p)
        except Exception as e:
            print(f"  [skip] {p.stem}: {type(e).__name__}: {e}", flush=True)
            gc.collect()
            continue
        if out is None:
            print(f"  [skip] {p.stem}: 不足", flush=True)
            continue
        r, pf = out
        rows.append(r)
        profs.append(pf)
        print(f"  [{i}/{len(files)}] {p.stem}: rho(dpt,dist)={r['rho_dpt_dist']:+.3f} "
              f"partial={r['rho_partial']:+.3f} root={r['root_ring']} "
              f"{time.time()-t0:.0f}s", flush=True)

    P = pd.DataFrame(rows)
    P.to_csv(PROC / "pseudotime_per_sample.csv.gz", index=False, compression="gzip")
    Pf = pd.concat(profs, ignore_index=True)
    Pf.to_csv(PROC / "pseudotime_profiles.csv", index=False)

    n = len(P)
    z = np.arctanh(np.clip(P["rho_dpt_dist"].to_numpy(float), -0.999, 0.999))
    zt = z.mean() / (z.std(ddof=1) / np.sqrt(n))
    zp = np.arctanh(np.clip(P["rho_partial"].to_numpy(float), -0.999, 0.999))
    zpt = zp.mean() / (zp.std(ddof=1) / np.sqrt(n))
    summ = dict(
        n_samples=n,
        median_rho_dpt_dist=float(np.median(P["rho_dpt_dist"])),
        mean_rho_dpt_dist=float(np.mean(P["rho_dpt_dist"])),
        frac_positive=float((P["rho_dpt_dist"] > 0).mean()),
        mean_rho_partial=float(np.mean(P["rho_partial"])),
        t_across_samples=float(zt),
        p_across_samples=float(2 * stats.t.sf(abs(zt), df=n - 1)),
        t_partial_across_samples=float(zpt),
        p_partial_across_samples=float(2 * stats.t.sf(abs(zpt), df=n - 1)),
        frac_root_is_nerve=float(P["iroot_is_nerve"].mean()),
    )
    summ["frac_root_in_r1"] = float((P["root_ring"] == "r1_0_50").mean())
    pd.DataFrame([summ]).to_csv(TAB / "pseudotime_summary.csv", index=False)
    print("\n=== 伪时间-距离一致性 ===")
    for k, v in summ.items():
        print(f"  {k}: {v}")

    # ---- 读数峰值顺序 ----
    num = Pf.drop(columns=["sample", "bin", "n_spot", "dpt_median", "dist_median"])
    m = num.groupby(Pf["bin"]).mean()
    order = []
    for c in m.columns:
        v = m[c].to_numpy(float)
        v = v - v.mean()
        order.append(dict(readout=c, peak_bin=int(np.argmax(v)),
                         amplitude=float(v.max() - v.min()),
                         rho_with_bin=float(stats.spearmanr(
                             np.arange(len(v)), v).statistic)))
    O = pd.DataFrame(order).sort_values("peak_bin")
    O.to_csv(TAB / "pseudotime_readout_order.csv", index=False)
    print("\n=== 读数沿伪时间的峰值顺序 (0=神经端, 9=远端) ===")
    print(O.to_string(index=False))
    return P, Pf


if __name__ == "__main__":
    main()

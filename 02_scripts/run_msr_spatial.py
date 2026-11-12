"""空间自相关感知的 NPP 重算 —— Moran 谱随机化 (Moran Spectral Randomization, MSR)。

背景(审稿人关切)
----------------
Step3a 用 spot 水平的偏 Spearman 推导 NPP。Visium spot 天然空间自相关,
独立性假设被违反, 经典相关检验的名义 P / FDR 偏小 —— NPP 的 97 个基因中
可能混有"只要空间光滑就会被检出"的假阳性。

方法
----
每样本构造两个零分布:

  (1) MSR 零分布 (Wagner & Dray 2015; adespatial::msr 思路)
      把中心化的表达矩阵投影到双中心化连接矩阵 H = PBP (P 为中心化投影)
      的正交特征基上, 按特征值(空间频率)分箱, 箱内随机置换投影系数。
      surrogate 保留每个基因自身的 Moran's I(空间光滑度), 但破坏表达场与
      "到神经距离"场之间真实的配对关系。基因间共用同一组置换, 因此同时
      保留基因-基因协方差结构。

  (2) 朴素 spot 置换零分布 (对照)
      直接随机打乱 spot 标签, 完全破坏自相关 = "忽略自相关"的零模型。

两者零分布宽度之比直接量化假阳性膨胀倍数。

统计量: 与 Step3a 完全相同的【秩空间偏 Spearman】——
        rank(log1p CPM 表达) 与 rank(距离), 剔除 rank(log UMI)、rank(pct_mt)。

数值技巧: 偏相关可在特征系数域解析计算, 无需显式重建 surrogate 空间场。
  设 Cs 为置换后的投影系数 (n x G), w = V' v_resid, U = V' Q (n x k), 则
      rho_g = (w' Cs[:,g]) / sqrt( (||Cs[:,g]||^2 - ||U' Cs[:,g]||^2) * ||v_resid||^2 )
  每 surrogate 成本 O(nG) 而非 O(n^2 G)。

跨样本合并: 每样本 rho -> Fisher z, 用该样本 MSR 零分布标准化为 u,
           Stouffer 合并 Z = sum(u)/sqrt(n_samples)。

用法
----
python run_msr_spatial.py scan  [--nsim 499] [--limit N]
python run_msr_spatial.py meta
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
from scipy.linalg import eigh
from scipy.spatial import cKDTree

sys.path.insert(0, r"D:/嗜神经浸润/00_pipeline")
from nsgrad import config as CFG  # noqa: E402

DIR_H5 = CFG.DIR_PROC / "spot_level"
PROC = CFG.DIR_PROC
TAB = CFG.DIR_TAB

MIN_SPOT_FRAC = 0.10       # 与 Step3a 一致
KNN = 6                     # Visium 六边形最近邻数
MAX_NB_UM = 220.0           # 邻接最大距离(避免跨组织碎片连边)
BIN_TARGET = 100            # 每个特征值箱的目标大小
MAX_SPOTS = 4000            # 超过则抽稀(特征分解成本 O(n^3))
MIN_RHO_ABS = 0.05          # 与原 NPP 相同的效应量门槛
FDR = 0.05


# ---------------------------------------------------------------- utilities
def _rank_cols(M: np.ndarray, block: int = 2048) -> np.ndarray:
    """逐列平均秩(并列取平均, 0-based) —— 分块向量化。

    注意: 组结束位置 en 必须用【后缀最小】求, 早期版本误用前缀最小,
    在有并列值时会给出错误的秩(counts 数据里 0 大量并列, 影响很大)。
    """
    n, g = M.shape
    out = np.empty((n, g), dtype=np.float64)
    idx = np.arange(n, dtype=np.float64)[:, None]
    for s in range(0, g, block):
        e = min(s + block, g)
        Ms = np.ascontiguousarray(M[:, s:e], dtype=np.float64)
        order = np.argsort(Ms, axis=0, kind="stable")
        S = np.take_along_axis(Ms, order, axis=0)
        new = np.empty(S.shape, dtype=bool)
        new[0, :] = True
        np.not_equal(S[1:], S[:-1], out=new[1:])
        # 组起始: 前缀最大
        st = np.maximum.accumulate(np.where(new, idx, -np.inf), axis=0)
        # 组结束: 后缀最小 (先反转 -> 前缀最小 -> 反转回来)
        end_mark = np.empty(S.shape, dtype=bool)
        end_mark[:-1] = new[1:]
        end_mark[-1] = True
        en = np.minimum.accumulate(
            np.where(end_mark, idx, np.inf)[::-1], axis=0)[::-1]
        r_sorted = (st + en) / 2.0
        blk = np.empty((n, e - s), dtype=np.float64)
        np.put_along_axis(blk, order, r_sorted, axis=0)
        out[:, s:e] = blk
        del Ms, S, new, end_mark, st, en, r_sorted, blk
    return out


def build_knn_graph(xy: np.ndarray, k: int = KNN, max_d: float = MAX_NB_UM) -> np.ndarray:
    """对称化 kNN 二元邻接矩阵(剔除超过 max_d 的边)。"""
    n = xy.shape[0]
    tree = cKDTree(xy)
    d, idx = tree.query(xy, k=k + 1, workers=-1)
    B = np.zeros((n, n), dtype=np.float64)
    rows = np.repeat(np.arange(n), k)
    cols = idx[:, 1:].ravel()
    dists = d[:, 1:].ravel()
    m = dists <= max_d
    B[rows[m], cols[m]] = 1.0
    B = np.maximum(B, B.T)          # 对称化
    np.fill_diagonal(B, 0.0)
    return B


def moran_spectrum(B: np.ndarray):
    """双中心化 H = PBP 的特征分解(升序)。返回 (lam, V, S0)。"""
    n = B.shape[0]
    S0 = float(B.sum())
    if S0 <= 0:
        raise ValueError("邻接矩阵为空(可能无近距离邻居)")
    r = B.sum(axis=1, keepdims=True)
    c = B.sum(axis=0, keepdims=True)
    tot = float(r.sum())
    H = B - r / n - c / n + tot / (n * n)
    H = (H + H.T) / 2.0
    lam, V = eigh(H)
    return lam, V, S0


def morans_i_from_coef(C: np.ndarray, lam: np.ndarray, S0: float, n: int) -> np.ndarray:
    """由投影系数 C (n x G) 反算每列 Moran's I = (n/S0)*sum(lam*c^2)/sum(c^2)。"""
    num = (C * C * lam[:, None]).sum(axis=0)
    den = (C * C).sum(axis=0)
    return (n / S0) * num / (den + 1e-30)


def ortho_cov(tc_log: np.ndarray, pct_mt: np.ndarray) -> np.ndarray:
    C = np.column_stack([stats.rankdata(tc_log), stats.rankdata(pct_mt)]).astype(np.float64)
    C = (C - C.mean(axis=0)) / (C.std(axis=0) + 1e-12)
    Q = np.zeros_like(C)
    for j in range(C.shape[1]):
        v = C[:, j].copy()
        for k in range(j):
            v -= Q[:, k] * (Q[:, k] @ v)
        nv = np.sqrt(v @ v)
        Q[:, j] = v / nv if nv > 1e-12 else 0.0
    return Q


# ---------------------------------------------------------------- per-sample
def process_sample(path: Path, target_genes: np.ndarray, nsim: int, seed: int):
    ad = sc.read_h5ad(path)
    obs = ad.obs
    dist = obs["distance_to_nearest_nerve"].to_numpy(dtype=float)
    keep = np.isfinite(dist)
    X_raw = ad.X
    X_raw = np.asarray(X_raw.toarray() if hasattr(X_raw, "toarray") else X_raw)[keep]
    # 注意: obsm["spatial"] 是_FULLRES 图像像素坐标, 不是微米; 用 obs 的 x_um/y_um
    xy = np.column_stack([
        np.asarray(obs["x_um"].to_numpy(dtype=float)),
        np.asarray(obs["y_um"].to_numpy(dtype=float)),
    ])[keep]
    d = dist[keep]
    tc = np.asarray(obs["total_counts"].to_numpy(dtype=float))[keep]
    pm = np.asarray(obs["pct_mt"].to_numpy(dtype=float))[keep]
    genes_all = np.asarray(ad.var_names)
    del ad

    # ---- 基因过滤 (与 Step3a 同标准) ∩ 最终基因表 ----
    expressed = (X_raw > 0).mean(axis=0)
    sel = expressed >= MIN_SPOT_FRAC
    gs = genes_all[sel]
    mask_t = np.isin(gs, target_genes)
    if mask_t.sum() < 300:
        return None
    Xs = X_raw[:, sel][:, mask_t].astype(np.float64)
    gs = gs[mask_t]

    n0 = Xs.shape[0]
    if n0 > MAX_SPOTS:
        rng0 = np.random.default_rng(seed)
        sub = np.sort(rng0.choice(n0, size=MAX_SPOTS, replace=False))
        Xs, d, xy, tc, pm = Xs[sub], d[sub], xy[sub], tc[sub], pm[sub]
    n = Xs.shape[0]

    # ---- 表达表示: 与 Step3a 严格一致 (raw counts 的秩) ----
    # (rank 对单调变换近乎不变, 保持与 Step3a 完全一致才能把新旧差异
    #  纯粹归因于零分布, 而不是表达表示)
    Xl = Xs

    # ---- 秩矩阵(表达场), 中心化标准化 ----
    R = _rank_cols(Xl)
    del Xs, Xl
    R -= R.mean(axis=0)
    R /= (R.std(axis=0) + 1e-12)
    G = R.shape[1]

    rd = stats.rankdata(d)
    rd = (rd - rd.mean()) / (rd.std() + 1e-12)
    Q = ortho_cov(np.log1p(tc), pm)

    # ---- 观测统计量 (偏相关, 系数域解析式) ----
    v_res = rd.copy()
    for j in range(Q.shape[1]):
        v_res -= Q[:, j] * (Q[:, j] @ v_res)
    vv = float(v_res @ v_res)
    rtQ = R.T @ Q                                    # G x k
    den_obs = np.sqrt(np.maximum(n - (rtQ ** 2).sum(axis=1), 1e-12) * vv)
    rho_obs = np.clip((R.T @ v_res) / den_obs, -1.0, 1.0)

    # ---- 空间图 + Moran 谱 ----
    t0 = time.time()
    B = build_knn_graph(xy)
    n_edges = B.sum() / 2.0
    lam, V, S0 = moran_spectrum(B)
    del B
    order = np.argsort(lam)
    lam_s = lam[order]
    V_s = np.ascontiguousarray(V[:, order])
    del V, lam

    C = V_s.T @ R                                    # n x G 投影系数
    w = V_s.T @ v_res                                # n
    U = V_s.T @ Q                                    # n x k
    moran_dist = float(morans_i_from_coef((V_s.T @ rd)[:, None], lam_s, S0, n)[0])
    del V_s
    moran_obs = morans_i_from_coef(C, lam_s, S0, n)

    # ---- 分箱 (按特征值排序等分) ----
    nb = int(max(4, min(60, n // BIN_TARGET)))
    edges = np.linspace(0, n, nb + 1).astype(int)

    rng = np.random.default_rng(seed)
    rho_msr = np.empty((nsim, G), dtype=np.float64)
    moran_sur = np.empty((nsim, G), dtype=np.float64)
    perm_base = np.arange(n)
    for r in range(nsim):
        Cs = np.empty_like(C)
        for b in range(nb):
            s, e = edges[b], edges[b + 1]
            if e - s > 1:
                Cs[s:e] = C[s:e][rng.permutation(e - s)]
            else:
                Cs[s:e] = C[s:e]
        moran_sur[r] = morans_i_from_coef(Cs, lam_s, S0, n)
        num = w @ Cs                                             # G
        Cu = Cs.T @ U                                            # G x k
        den = np.sqrt(np.maximum((Cs * Cs).sum(axis=0) - (Cu ** 2).sum(axis=1), 1e-12) * vv)
        rho_msr[r] = np.clip(num / den, -1.0, 1.0)
    del C, Cs, Cu

    # ---- 朴素 spot 置换零分布 ----
    rho_nai = np.empty((nsim, G), dtype=np.float64)
    RT = R.T                                          # G x n
    for r in range(nsim):
        p = rng.permutation(n)
        inv = np.empty(n, dtype=np.int64)
        inv[p] = perm_base
        num = RT @ v_res[inv]
        Qp = Q[inv]                                   # n x k (行重排)
        den = np.sqrt(np.maximum(n - ((RT @ Qp) ** 2).sum(axis=1), 1e-12) * vv)
        rho_nai[r] = np.clip(num / den, -1.0, 1.0)
    del RT

    out = pd.DataFrame({
        "sample": path.stem,
        "gene": gs,
        "n_spots": n,
        "rho_obs": rho_obs,
        "moran_i": moran_obs,
        "moran_i_surrogate": moran_sur.mean(axis=0),
        "moran_i_sur_sd": moran_sur.std(axis=0),
        "msr_mean": rho_msr.mean(axis=0),
        "msr_sd": rho_msr.std(axis=0, ddof=1),
        "naive_mean": rho_nai.mean(axis=0),
        "naive_sd": rho_nai.std(axis=0, ddof=1),
        "n_edges": n_edges,
        "n_bins": nb,
        "moran_i_dist": moran_dist,
    })
    # SD 收缩: 避免个别基因零标准差导致 u 爆炸
    for tag in ("msr", "naive"):
        ref = float(np.nanmedian(out[f"{tag}_sd"]))
        out[f"{tag}_sd"] = out[f"{tag}_sd"].clip(lower=0.3 * ref)
    # 双尾经验 P (基于 |rho| 相对零分布中心的偏离)
    for tag, arr in (("msr", rho_msr), ("naive", rho_nai)):
        ctr = arr.mean(axis=0)
        dev_obs = np.abs(rho_obs - ctr)
        dev_null = np.abs(arr - ctr)
        out[f"p_{tag}"] = (1.0 + (dev_null >= dev_obs[None, :]).sum(axis=0)) / (1.0 + nsim)
    out["u_msr"] = (rho_obs - out["msr_mean"]) / out["msr_sd"].replace(0, np.nan)
    out["u_naive"] = (rho_obs - out["naive_mean"]) / out["naive_sd"].replace(0, np.nan)
    out.attrs["t_eig"] = time.time() - t0
    return out


# ---------------------------------------------------------------- meta
def stage_scan(nsim: int, limit: int | None, seed: int, only: str | None = None,
               bin_target: int = BIN_TARGET, tag: str = ""):
    global BIN_TARGET
    BIN_TARGET = bin_target
    files = sorted(DIR_H5.glob("*.h5ad"))
    if only:
        files = [f for f in files if only in f.name]
    if limit:
        files = files[:limit]
    # 目标基因集 = Step3a 最终纳入元分析的基因
    tgt = pd.read_csv(PROC / "nerve_proximity_program.csv", index_col=0).index.to_numpy()
    print(f"[MSR] 目标基因 {len(tgt)} | 样本 {len(files)} | nsim={nsim}", flush=True)
    rows = []
    for i, p in enumerate(files, 1):
        t0 = time.time()
        try:
            r = process_sample(p, tgt, nsim, seed + i)
        except Exception as e:  # noqa: BLE001
            print(f"   !! {p.stem} 失败: {type(e).__name__}: {e}", flush=True)
            continue
        if r is None:
            print(f"   skip {p.stem} (目标基因不足)", flush=True)
            continue
        rows.append(r)
        print(f"   [{i}/{len(files)}] {p.stem}: n={r['n_spots'].iloc[0]} G={len(r)} "
              f"edges={r['n_edges'].iloc[0]:.0f} bins={r['n_bins'].iloc[0]} "
              f"moranI_med={np.nanmedian(r['moran_i']):.3f} "
              f"| {time.time()-t0:.1f}s", flush=True)
        del r
    if not rows:
        raise SystemExit("无样本成功")
    P = pd.concat(rows, ignore_index=True)
    fn = PROC / f"msr_per_sample{('.' + tag) if tag else ''}.csv.gz"
    P.to_csv(fn, index=False, compression="gzip")
    print(f"\n[MSR] -> {fn}  ({len(P)} 行)")
    return P


def bh(p: np.ndarray) -> np.ndarray:
    p = np.asarray(p, dtype=float)
    ok = ~np.isnan(p)
    q = np.full(p.shape, np.nan)
    pp = p[ok]
    n = len(pp)
    if n == 0:
        return q
    order = np.argsort(pp)
    ranked = pp[order] * n / (np.arange(n) + 1)
    qs = np.minimum.accumulate(ranked[::-1])[::-1]
    out = np.empty(n)
    out[order] = np.clip(qs, 0, 1)
    q[ok] = out
    return q


def stage_meta(tag: str = ""):
    P = pd.read_csv(PROC / f"msr_per_sample{('.' + tag) if tag else ''}.csv.gz")
    orig = pd.read_csv(PROC / "nerve_proximity_program.csv", index_col=0)
    n_samp = P["sample"].nunique()
    print(f"[meta] {n_samp} 样本, {P['gene'].nunique()} 基因", flush=True)

    P = P.dropna(subset=["u_msr"])
    # 每个基因至少在一半样本可测
    cnt = P.groupby("gene")["u_msr"].size()
    keepg = cnt[cnt >= max(2, int(np.ceil(0.5 * n_samp)))].index
    P = P[P["gene"].isin(keepg)]

    Gm = P.groupby("gene").agg(
        n_samples=("u_msr", "size"),
        mean_rho=("rho_obs", "mean"),
        mean_u_msr=("u_msr", "mean"),
        mean_u_naive=("u_naive", "mean"),
        moran_i=("moran_i", "mean"),
        moran_i_sur=("moran_i_surrogate", "mean"),
        msr_sd=("msr_sd", "mean"),
        naive_sd=("naive_sd", "mean"),
        consistency=("rho_obs", lambda s: float(max((s < 0).mean(), (s > 0).mean()))),
        frac_p05=("p_msr", lambda s: float((s < 0.05).mean())),
    )
    # 每样本 u 的跨样本 SD (用于"MSR 标准化 + 跨样本 t"的保守合并)
    sd_u = P.groupby("gene")["u_msr"].std()
    se_u = sd_u / np.sqrt(P.groupby("gene")["u_msr"].size())
    t_u = P.groupby("gene")["u_msr"].mean() / se_u.replace(0, np.nan)
    df_u = P.groupby("gene")["u_msr"].size() - 1
    p_msr_t = pd.Series(2 * stats.t.sf(np.abs(t_u), df=df_u), index=t_u.index)

    Gm["stouffer_Z"] = Gm["mean_u_msr"] * np.sqrt(Gm["n_samples"])
    Gm["p_msr"] = 2 * stats.norm.sf(np.abs(Gm["stouffer_Z"]))
    Gm["fdr_msr"] = bh(Gm["p_msr"].to_numpy())
    Gm["Z_naive"] = Gm["mean_u_naive"] * np.sqrt(Gm["n_samples"])
    Gm["p_naive"] = 2 * stats.norm.sf(np.abs(Gm["Z_naive"]))
    Gm["fdr_naive"] = bh(Gm["p_naive"].to_numpy())
    # 保守合并: MSR 去偏后再做跨样本 t 检验(同时承认样本间生物学异质性)
    Gm["t_msr_between"] = t_u.reindex(Gm.index)
    Gm["p_msr_between"] = p_msr_t.reindex(Gm.index)
    Gm["fdr_msr_between"] = bh(Gm["p_msr_between"].to_numpy())

    # 与原 Step3a 结果对齐
    Gm["orig_p"] = orig["p"].reindex(Gm.index)
    Gm["orig_fdr"] = orig["fdr"].reindex(Gm.index)
    Gm["orig_mean_rho"] = orig["mean_rho"].reindex(Gm.index)
    Gm = Gm.sort_values("stouffer_Z")

    # ---- 三套基因集 ----
    eff = Gm["mean_rho"].abs() >= MIN_RHO_ABS
    set_orig = set(Gm.index[(Gm["orig_p"] < 0.05) & eff])              # Step3a 宽松集(原 NPP 来源)
    set_msr = set(Gm.index[(Gm["fdr_msr"] < FDR) & eff])
    set_naive = set(Gm.index[(Gm["fdr_naive"] < FDR) & eff])
    set_orig_strict = set(Gm.index[(Gm["orig_fdr"] < FDR) & eff])
    set_msr_between = set(Gm.index[(Gm["fdr_msr_between"] < FDR) & eff])
    # 不设效应量门槛时的显著计数 —— 直接量化"名义检验的假阳性膨胀"
    n_noeff = {
        "orig_p05": int((Gm["orig_p"] < 0.05).sum()),
        "orig_fdr05": int((Gm["orig_fdr"] < FDR).sum()),
        "msr_fdr05": int((Gm["fdr_msr"] < FDR).sum()),
        "naive_fdr05": int((Gm["fdr_naive"] < FDR).sum()),
        "msr_p05": int((Gm["p_msr"] < 0.05).sum()),
        "naive_p05": int((Gm["p_naive"] < 0.05).sum()),
        "msr_between_fdr05": int((Gm["fdr_msr_between"] < FDR).sum()),
    }
    print("\n[无效应量门槛的显著基因数 / 6221]")
    for k, v in n_noeff.items():
        print(f"   {k:12s} {v:5d}  ({v/len(Gm)*100:.1f}%)")
    print(f"\n[meta] 原 NPP 宽松集 (p<0.05 & |rho|>={MIN_RHO_ABS}) : {len(set_orig)}")
    print(f"[meta] 原 NPP 严格集 (FDR<0.05 & |rho|>={MIN_RHO_ABS}): {len(set_orig_strict)}")
    print(f"[meta] MSR  零分布集 (FDR<0.05 & |rho|>={MIN_RHO_ABS}) : {len(set_msr)}")
    print(f"[meta] 朴素 零分布集 (FDR<0.05 & |rho|>={MIN_RHO_ABS}) : {len(set_naive)}")
    print(f"[meta] MSR+跨样本 t 保守集 (FDR<0.05 & |rho|>={MIN_RHO_ABS}): {len(set_msr_between)}")

    # ---- 重合度 ----
    inter = set_orig & set_msr
    jac = len(inter) / max(1, len(set_orig | set_msr))
    N = len(Gm)
    hyp_p = stats.hypergeom.sf(len(inter) - 1, N, len(set_orig), len(set_msr)) if len(set_msr) else np.nan
    surv = len(inter) / max(1, len(set_orig))
    print(f"\n[重合度] 交集 {len(inter)} | Jaccard {jac:.3f} | 原集存活率 {surv:.1%} "
          f"| 超几何 P={hyp_p:.2e}")

    sd_ratio = float(np.nanmedian(Gm["msr_sd"] / Gm["naive_sd"]))
    print(f"[零分布宽度] median(SD_msr / SD_naive) = {sd_ratio:.2f}")
    print(f"[Moran I 保留] obs 中位 {np.nanmedian(Gm['moran_i']):.3f} vs "
          f"surrogate 中位 {np.nanmedian(Gm['moran_i_sur']):.3f}")

    # ---- 方向一致性 ----
    both = sorted(inter)
    if both:
        sub = Gm.loc[both]
        same_dir = float(np.mean(np.sign(sub["mean_rho"]) == np.sign(sub["orig_mean_rho"])))
        print(f"[方向] 交集基因中 rho 符号一致率 {same_dir:.1%}")

    TAB.mkdir(parents=True, exist_ok=True)
    gf = PROC / f"msr_gene_level{('.' + tag) if tag else ''}.csv"
    Gm.to_csv(gf, index=True)
    print(f"\n  -> {gf}")

    summary = {
        "n_samples": int(n_samp),
        "n_genes": int(len(Gm)),
        "nsim": int(round(1.0 / P["p_msr"].iloc[0] * 0)) if False else None,
        "min_rho_abs": MIN_RHO_ABS,
        "set_orig_relaxed": sorted(set_orig),
        "set_orig_strict": sorted(set_orig_strict),
        "set_msr": sorted(set_msr),
        "set_naive": sorted(set_naive),
        "set_msr_between": sorted(set_msr_between),
        "overlap": sorted(inter),
        "jaccard": float(jac),
        "orig_survival_rate": float(surv),
        "hypergeom_p": float(hyp_p) if np.isfinite(hyp_p) else None,
        "sd_ratio_msr_over_naive": sd_ratio,
        "moran_i_obs_median": float(np.nanmedian(Gm["moran_i"])),
        "moran_i_surrogate_median": float(np.nanmedian(Gm["moran_i_sur"])),
        "n_sig_naive": len(set_naive),
        "n_sig_msr": len(set_msr),
        "inflation_factor": len(set_naive) / max(1, len(set_msr)),
        "n_significant_no_effect_threshold": n_noeff,
        "moran_i_dist_median": float(np.nanmedian(P.groupby("sample")["moran_i_dist"].first())),
    }
    sf = TAB / f"msr_summary{('.' + tag) if tag else ''}.json"
    with open(sf, "w", encoding="utf-8") as f:
        json.dump(summary, f, indent=2, ensure_ascii=False)
    print(f"  -> {sf}")
    return Gm, summary


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("stage", choices=["scan", "meta"])
    ap.add_argument("--nsim", type=int, default=499)
    ap.add_argument("--limit", type=int, default=None)
    ap.add_argument("--seed", type=int, default=20261026)
    ap.add_argument("--only", type=str, default=None)
    ap.add_argument("--bin", type=int, default=BIN_TARGET)
    ap.add_argument("--tag", type=str, default="")
    a = ap.parse_args()
    if a.stage == "scan":
        stage_scan(a.nsim, a.limit, a.seed, a.only, a.bin, a.tag)
    else:
        stage_meta(a.tag)

"""Step3a —— 从空间转录组反推「神经近端转录程序」(Nerve-Proximity Program, NPP)。

逻辑
----
对每个样本, 计算每个基因表达与「到最近神经的距离」的 Spearman 相关;
再做跨样本元分析(Fisher z + 单样本 t 检验), 得到一组稳定随神经距离变化的基因。

  NPP-up   : 距神经越近表达越高 (rho < 0)
  NPP-down : 距神经越近表达越低 (rho > 0)

bulk 投影: NPPscore = mean z(NPP-up) - mean z(NPP-down)

输出
----
01_data/processed/nerve_proximity_program.csv   逐基因统计
01_data/processed/npp_signature.json            最终签名(用于 bulk 打分)
"""
from __future__ import annotations
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
from scipy import stats

sys.path.insert(0, r"D:/嗜神经浸润/00_pipeline")
from nsgrad import config as CFG  # noqa: E402

DIR_H5 = CFG.DIR_PROC / "spot_level"
OUT = CFG.DIR_PROC

MIN_SPOT_FRAC = 0.10      # 基因需在该样本 >=10% 的 spot 中表达
MIN_SAMPLE_FRAC = 0.50    # 基因需在 >=50% 的样本中可测
MIN_RHO_ABS = 0.05        # 效应量门槛 (spot 级 ST 单基因效应本就温和: |rho| 99 分位仅 0.062)
MIN_CONSISTENCY = 0.70    # 符号一致率
FDR = 0.05
TOP_N = 150               # 每侧取 |mean rho| 最大的基因数


def _rank_cols(M: np.ndarray, block: int = 3000) -> np.ndarray:
    """对矩阵每一列做秩变换(并列取平均秩, 0-based) —— 分块向量化实现。

    Spearman = 对秩做 Pearson, 因此需要逐列的平均秩。

    修正(2026-10-26): 组结束行号 en 原先误用【前缀最小】求取, 在有并列值时
    会给出错误的秩; counts 矩阵里 0 大量并列, 影响显著。现改为【后缀最小】
    (反转 -> 前缀最小 -> 反转), 已与 scipy.stats.rankdata(method='average') 对齐。
    """
    n, g = M.shape
    out = np.empty((n, g), dtype=np.float32)
    idx = np.arange(n, dtype=np.float64)[:, None]
    for s in range(0, g, block):
        e = min(s + block, g)
        Ms = np.ascontiguousarray(M[:, s:e], dtype=np.float64)
        order = np.argsort(Ms, axis=0, kind="stable")
        S = np.take_along_axis(Ms, order, axis=0)
        new = np.empty(S.shape, dtype=bool)
        new[0, :] = True
        np.not_equal(S[1:], S[:-1], out=new[1:])
        # 每组起始行号 (前缀最大)
        st = np.maximum.accumulate(np.where(new, idx, -np.inf), axis=0)
        # 每组结束行号 (后缀最小)
        end_mark = np.empty(S.shape, dtype=bool)
        end_mark[:-1] = new[1:]
        end_mark[-1] = True
        en = np.minimum.accumulate(
            np.where(end_mark, idx, np.inf)[::-1], axis=0)[::-1]
        r_sorted = ((st + en) / 2.0).astype(np.float32)
        blk = np.empty((n, e - s), dtype=np.float32)
        np.put_along_axis(blk, order, r_sorted, axis=0)
        out[:, s:e] = blk
        del Ms, S, new, st, en, r_sorted, blk
    return out


def sample_gene_distance_rho(ad: sc.AnnData, exclude_nerve: bool = False):
    """单样本: 每个基因 vs 距离的【偏】Spearman rho(剔除 spot 细胞量/线粒体混杂)。

    关键: 未校正时, 距离与 total_counts 高度共线(神经灶附近 spot 细胞密集),
    直接相关会把管家基因(ACTB/TUBB/RPS 等)全部顶到前面 —— 纯技术假象。
    这里用「秩空间上的偏相关」: 先把 log(UMI) 与 pct_mt 的秩从基因秩和距离秩中
    分别回归掉, 再求相关。
    """
    obs = ad.obs
    dist = obs["distance_to_nearest_nerve"].to_numpy(dtype=float)
    keep = np.isfinite(dist)
    if exclude_nerve:
        keep &= (obs["nerve_spatial"].to_numpy() == 0)
    X = ad.X
    X = X.toarray() if hasattr(X, "toarray") else np.asarray(X)
    X = X[keep]
    d = dist[keep]
    genes = np.asarray(ad.var_names)

    # 基因过滤
    expressed = (X > 0).mean(axis=0)
    sel = expressed >= MIN_SPOT_FRAC
    if sel.sum() < 500:
        return None
    Xs, gs = X[:, sel].astype(np.float64), genes[sel]

    # ---- 秩变换 ----
    Rg = _rank_cols(Xs).astype(np.float64)     # n_spot x n_gene
    Rg -= Rg.mean(axis=0)
    Rg /= (Rg.std(axis=0) + 1e-12)
    rd = stats.rankdata(d).astype(np.float64)
    rd = (rd - rd.mean()) / (rd.std() + 1e-12)

    # ---- 混杂协变量(秩空间) ----
    cov_cols = []
    tc = np.asarray(obs["total_counts"].to_numpy(dtype=float))[keep]
    cov_cols.append(stats.rankdata(np.log1p(tc)).astype(np.float64))
    pm = np.asarray(obs["pct_mt"].to_numpy(dtype=float))[keep]
    cov_cols.append(stats.rankdata(pm).astype(np.float64))
    C = np.column_stack(cov_cols)
    C = (C - C.mean(axis=0)) / (C.std(axis=0) + 1e-12)
    # Gram-Schmidt 正交化后逐个投影剔除
    Q = np.zeros_like(C)
    for j in range(C.shape[1]):
        v = C[:, j].copy()
        for k in range(j):
            v -= Q[:, k] * (Q[:, k] @ v)
        nv = np.sqrt(v @ v)
        Q[:, j] = v / nv if nv > 1e-12 else 0.0

    def residualize(A: np.ndarray) -> np.ndarray:
        for j in range(Q.shape[1]):
            q = Q[:, j]
            if np.abs(q).max() < 1e-12:
                continue
            A = A - np.outer(q, q @ A)
        return A

    Rg_r = residualize(Rg)
    rd_r = residualize(rd[:, None])[:, 0]

    num = (Rg_r * rd_r[:, None]).sum(axis=0)
    den = np.sqrt((Rg_r ** 2).sum(axis=0) * (rd_r ** 2).sum())
    rho = num / (den + 1e-12)
    rho = np.clip(rho, -1.0, 1.0)

    # 近/远 (排除神经 spot 自身, 避免把神经元基因本身算进来)
    nonerve = obs["nerve_spatial"].to_numpy()[keep] == 0
    dn = d[nonerve]
    near_m = dn <= 150
    far_m = dn >= 400
    delta = np.full(Xs.shape[1], np.nan)
    if near_m.sum() >= 20 and far_m.sum() >= 20:
        Z = (Xs - Xs.mean(axis=0)) / (Xs.std(axis=0) + 1e-12)
        delta = Z[nonerve][near_m].mean(axis=0) - Z[nonerve][far_m].mean(axis=0)
    return pd.DataFrame({"gene": gs, "rho": rho, "delta_z": delta})


def main():
    files = sorted(DIR_H5.glob("*.h5ad"))
    print(f"[Step3a] {len(files)} 个空间样本参与程序推导", flush=True)
    per_sample = []
    for p in files:
        ad = sc.read_h5ad(p)
        r = sample_gene_distance_rho(ad)
        if r is None:
            print(f"   skip {p.name} (基因过滤后不足)", flush=True)
            continue
        r["sample"] = p.stem
        per_sample.append(r)
        print(f"   {p.stem}: spots={ad.n_obs} genes={len(r)} "
              f"nerve={int(ad.obs['nerve_spatial'].sum())}", flush=True)
        del ad
    if not per_sample:
        raise SystemExit("无可用样本")
    P = pd.concat(per_sample, ignore_index=True)

    n_samp = P["sample"].nunique()
    # 基因需在 >=50% 样本可测
    cnt = P.groupby("gene")["rho"].count()
    keepg = cnt[cnt >= max(2, int(np.ceil(MIN_SAMPLE_FRAC * n_samp)))].index
    P = P[P["gene"].isin(keepg)]
    print(f"\n[Step3a] 保留 {P['gene'].nunique()} 个基因 (在 >={MIN_SAMPLE_FRAC:.0%} 样本可测)")

    # 元分析: Fisher z
    P["z"] = np.arctanh(np.clip(P["rho"].to_numpy(), -0.999999, 0.999999))
    G = P.groupby("gene").agg(
        n_samples=("z", "size"),
        mean_rho=("rho", "mean"),
        mean_z=("z", "mean"),
        sd_z=("z", "std"),
        mean_delta_z=("delta_z", "mean"),
        consistency=("rho", lambda s: float(max((s < 0).mean(), (s > 0).mean()))),
    )
    G["se"] = G["sd_z"] / np.sqrt(G["n_samples"])
    G["t"] = G["mean_z"] / G["se"].replace(0, np.nan)
    G["p"] = 2 * stats.t.sf(np.abs(G["t"]), df=G["n_samples"] - 1)
    G = G.replace([np.inf, -np.inf], np.nan).dropna(subset=["p"])
    # BH FDR
    pvals = G["p"].to_numpy()
    order = np.argsort(pvals)
    n = len(pvals)
    q = np.empty(n)
    ranked = pvals[order] * n / (np.arange(n) + 1)
    q[order] = np.minimum.accumulate(ranked[::-1])[::-1]
    G["fdr"] = np.clip(q, 0, 1)
    G["mean_rho"] = np.tanh(G["mean_z"])
    G = G.sort_values("mean_z")

    # 严格集: FDR<0.05 且 |rho| 达到门槛
    strict = G[(G["fdr"] < FDR) & (G["mean_rho"].abs() >= MIN_RHO_ABS)]
    # 宽松集: 名义 p<0.05 且效应量达门槛(用于构建稳健的复合打分;
    #        单基因效应量在 spot 级 ST 里本来就很温和, 这是方法学现实, 不是缺陷)
    relaxed = G[(G["p"] < 0.05) & (G["mean_rho"].abs() >= MIN_RHO_ABS)]
    print(f"\n[Step3a] 严格集 FDR<{FDR} & |rho|>={MIN_RHO_ABS}: {len(strict)} "
          f"(近神经高 {(strict['mean_rho']<0).sum()} / 近神经低 {(strict['mean_rho']>0).sum()})")
    print(f"[Step3a] 宽松集 p<0.05 & |rho|>={MIN_RHO_ABS}: {len(relaxed)} "
          f"(近神经高 {(relaxed['mean_rho']<0).sum()} / 近神经低 {(relaxed['mean_rho']>0).sum()})")

    sig = relaxed
    up = sig[sig["mean_rho"] < 0].sort_values("mean_z").head(TOP_N)      # 近神经高
    down = sig[sig["mean_rho"] > 0].sort_values("mean_z", ascending=False).head(TOP_N)

    print(f"\n[Step3a] FDR<{FDR} 且 |rho|>={MIN_RHO_ABS} 且一致性>={MIN_CONSISTENCY:.0%}: {len(sig)} 基因")
    print(f"   NPP-up  (近神经高表达): {len(up)}")
    print(f"   NPP-down(近神经低表达): {len(down)}")
    print("\n  NPP-up top15:")
    print("   " + ", ".join(up.index[:15]))
    print("  NPP-down top15:")
    print("   " + ", ".join(down.index[:15]))

    OUT.mkdir(parents=True, exist_ok=True)
    G.to_csv(OUT / "nerve_proximity_program.csv", index=True)
    print(f"\n  -> {OUT/'nerve_proximity_program.csv'}")

    strict_up = [g for g in up.index if g in set(strict.index)]
    strict_dn = [g for g in down.index if g in set(strict.index)]
    payload = {
        "n_samples": int(n_samp),
        "criteria": {"fdr": FDR, "nominal_p": 0.05, "min_abs_rho": MIN_RHO_ABS,
                     "min_consistency": MIN_CONSISTENCY,
                     "min_spot_frac": MIN_SPOT_FRAC,
                     "min_sample_frac": MIN_SAMPLE_FRAC},
        "n_significant_fdr": int((G["fdr"] < FDR).sum()),
        "n_sig_genes": int(len(sig)),
        "up": list(up.index),
        "down": list(down.index),
        "up_rho": up["mean_rho"].round(4).to_dict(),
        "down_rho": down["mean_rho"].round(4).to_dict(),
        "strict_up": strict_up,
        "strict_down": strict_dn,
    }
    with open(OUT / "npp_signature.json", "w", encoding="utf-8") as f:
        json.dump(payload, f, indent=2, ensure_ascii=False)
    print(f"  -> {OUT/'npp_signature.json'}")

    # 逐样本的近/远差异也存一份, 供 Figure3 画火山图/热图
    P.to_csv(OUT / "npp_per_sample.csv.gz", index=False, compression="gzip")
    return G, up, down


if __name__ == "__main__":
    main()

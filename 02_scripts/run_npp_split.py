"""方案四: 把神经邻近程序(NPP)拆成"神经固有信号"与"神经诱导基质信号"。

科学问题
--------
NPP 在 bulk 中与通用 CAF/ECM 高度共线(rho=0.72), 导致预后信号随基因集边界摆动。
假设: NPP 其实是两类信号的混合 ——
  (A) Neural-intrinsic: 神经元/胶质/雪旺细胞固有的结构蛋白与递质合成酶
      (UCHL1/PGP9.5, TUBB3, MAP2, NEFL/NEFM, PLP1, MPZ, TH, DBH, CHAT ...)
  (B) Perineural-stroma: 神经诱导的基质/CAF 反应
      (POSTN, SPARC, FN1, 胶原, TAGLN, ACTA2, MMP ...)
若既往"PNI 签名"把两者混在一起, 则 bulk 层面的预后信号必然被(A)稀释、
被(B)的纤维化混杂主导 —— 这解释了领域内长期矛盾。

检验
----
空间(30 Visium): 两模块各自的神经距离梯度(MSR 校准), 以及【互相控制后】的
  独立梯度 —— 这是"两个信号可分离"的直接证据。
bulk(METABRIC 1980): 两模块的共线性对比、预后 Cox(单变量/互相调整/全调整)、
  bootstrap 稳定性。
外部(GSE12276 204): 两模块与脑转移竞争风险。

用法
----
python run_npp_split.py spatial
python run_npp_split.py bulk
python run_npp_split.py all
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

from nsgrad import config as CFG, surv as nsurv        # noqa: E402
from run_msr_spatial import (                          # noqa: E402
    _rank_cols, build_knn_graph, moran_spectrum,
    morans_i_from_coef, ortho_cov, MIN_SPOT_FRAC,
)

DIR_H5 = CFG.DIR_PROC / "spot_level"
PROC = CFG.DIR_PROC
TAB = CFG.DIR_TAB
META_FINE = PROC / "msr_gene_level.fine.csv"

# ------------------------------------------------------------------ 模块定义
# (A) 神经固有: 神经元/胶质/雪旺细胞结构蛋白 + 递质合成酶 + 递质肽
NEURAL_INTRINSIC = [
    "UCHL1", "TUBB3", "TUBA1A", "MAP2", "NEFL", "NEFM", "PLP1", "MPZ", "PMP22",
    "S100B", "NGFR", "SOX10", "SYP", "SCN9A", "SCN10A", "SCN7A", "PRX", "MBP",
    "TH", "DBH", "CHAT", "SLC18A3", "SLC5A7", "TRPV1", "CALCA", "CALCB", "VIP",
    "TAC1", "NPY", "SST", "PENK", "GAD1", "GAD2", "SLC6A1", "GAP43", "ELAVL4",
    "STMN2", "NRXN1", "SYT1",
]
# (B) 神经诱导基质: CAF/ECM (与 A 互斥)
PERINEURAL_STROMA = [
    "POSTN", "SPARC", "FN1", "CCN2", "THBS1", "THBS2", "COMP", "LRRC15",
    "COL1A1", "COL1A2", "COL3A1", "COL5A1", "COL5A2", "COL4A1", "COL4A2",
    "COL8A1", "COL10A1", "COL11A1", "COL12A1", "COL15A1", "COL18A1",
    "TAGLN", "ACTA2", "MYL9", "CNN1", "VIM", "TPM1", "FAP", "SULF1", "ASPN",
    "MMP1", "MMP2", "MMP9", "MMP11", "MMP13", "MMP14", "ADAM12", "LOX", "LOXL2",
    "TNC", "VCAN", "DCN", "LUM", "SERPINF1", "FSTL1", "CTHRC1", "SFRP1",
    "SERPINE1", "INHBA", "TGFB1", "PDGFRB", "PDGFRA", "RGS5", "MCAM",
    "CDH11", "NID1", "TGM2", "ANTXR1", "DKK3", "S100A4",
]
# 定义神经 spot 用的标记 —— 用于敏感性分析(剔除后重算, 排除循环论证)
NERVE_CALL_MARKERS = {
    "UCHL1", "TRPV1", "CALCA", "S100B", "SCN9A", "SCN10A", "NGF",
    "TH", "DBH", "SLC18A3", "CHAT", "TUBB3", "MAP2", "NEFL", "NEFM",
    "PLP1", "MPZ", "SYP",
}

MODULES = {
    "Neural_intrinsic": NEURAL_INTRINSIC,
    "Perineural_stroma": PERINEURAL_STROMA,
}


def msr_significant_genes(fdr=0.05, min_rho=0.05):
    G = pd.read_csv(META_FINE, index_col=0)
    m = G[(G["fdr_msr"] < fdr) & (G["mean_rho"].abs() >= min_rho)]
    return G, set(m.index)


# ------------------------------------------------------------- 空间: 梯度
def _module_scores(Xl, genes, mods, min_frac=0.30):
    idx = {g: i for i, g in enumerate(genes)}
    out = {}
    for nm, gl in mods.items():
        rows = [idx[g] for g in gl if g in idx]
        if len(rows) / max(1, len(gl)) < min_frac:
            continue
        out[nm] = Xl[:, rows].mean(axis=1)
    return out


def run_spatial(nsim=499, bin_target=40, seed=11):
    files = sorted(DIR_H5.glob("*.h5ad"))
    rows = []
    t0 = time.time()
    for i, p in enumerate(files, 1):
        ad = sc.read_h5ad(p)
        obs = ad.obs
        dist = obs["distance_to_nearest_nerve"].to_numpy(float)
        keep = np.isfinite(dist)
        X = np.asarray(ad.X.toarray())[keep]
        xy = np.column_stack([obs["x_um"].to_numpy(float),
                              obs["y_um"].to_numpy(float)])[keep]
        d = dist[keep]
        tc = np.asarray(obs["total_counts"].to_numpy(float))[keep]
        pm = np.asarray(obs["pct_mt"].to_numpy(float))[keep]
        gnames = np.asarray(ad.var_names)
        del ad
        # 模块基因用宽松检出标准: 神经固有基因极稀疏(神经 spot 仅占 ~2%),
        # 若套用"≥10% spots"的常规过滤会把 UCHL1/TUBB3/MAP2 全部抹掉。
        # 这里要求基因在 ≥5 个 spot 中检出, 且 logCPM 方差 > 0。
        lib = np.where(tc <= 0, 1.0, tc)
        Xl = np.log1p(X / lib[:, None] * 1e4)
        detected = (X > 0).sum(axis=0) >= 5
        keep_g = detected & (Xl.std(axis=0) > 1e-9)
        Xl = Xl[:, keep_g]
        gnames = gnames[keep_g]
        Xl = (Xl - Xl.mean(axis=0)) / (Xl.std(axis=0) + 1e-12)

        # 模块评分 + 敏感性版本(剔除神经 spot 定义标记)
        mods = dict(MODULES)
        mods_naive = {
            "Neural_intrinsic_noCallMarkers":
                [g for g in NEURAL_INTRINSIC if g not in NERVE_CALL_MARKERS],
            "Perineural_stroma": PERINEURAL_STROMA,
        }
        S = _module_scores(Xl, gnames, mods, min_frac=0.15)
        if "Neural_intrinsic" not in S:
            print(f"  [skip] {p.stem}: 神经模块覆盖不足")
            continue
        Sn = _module_scores(Xl, gnames, mods_naive, min_frac=0.15)
        del X

        # ---- 秩空间 + 协变量 ----
        M = np.column_stack([S["Neural_intrinsic"], S["Perineural_stroma"]])
        extra = []
        if "Neural_intrinsic_noCallMarkers" in Sn:
            extra = [Sn["Neural_intrinsic_noCallMarkers"]]
        Mall = np.column_stack([M] + extra) if extra else M
        names = ["Neural_intrinsic", "Perineural_stroma"] + \
                (["Neural_noCallMarkers"] if extra else [])
        n, K = Mall.shape
        R = _rank_cols(Mall)
        R -= R.mean(axis=0)
        R /= (R.std(axis=0) + 1e-12)
        rd = stats.rankdata(d)
        rd = (rd - rd.mean()) / (rd.std() + 1e-12)
        Q0 = ortho_cov(np.log1p(tc), pm)
        v_res = rd.copy()
        for j in range(Q0.shape[1]):
            v_res -= Q0[:, j] * (Q0[:, j] @ v_res)
        vv = float(v_res @ v_res)

        def partial_rho(Mx, cov=None):
            """Mx (n x k) 对 v_res 的偏相关, cov 为额外控制的 (n x m) 已中心化列。"""
            Ql = Q0 if cov is None else np.column_stack([Q0, cov])
            # Gram-Schmidt 正交化
            Qn = np.zeros_like(Ql)
            for j in range(Ql.shape[1]):
                v = Ql[:, j].astype(float).copy()
                for k2 in range(j):
                    v -= Qn[:, k2] * (Qn[:, k2] @ v)
                nv = np.sqrt(v @ v)
                Qn[:, j] = v / nv if nv > 1e-12 else 0.0
            vv2 = v_res.copy()
            for j in range(Qn.shape[1]):
                vv2 -= Qn[:, j] * (Qn[:, j] @ vv2)
            vv2n = float(vv2 @ vv2)
            Mx = Mx - Mx.mean(axis=0)
            Mxc = Mx.copy()
            for j in range(Qn.shape[1]):
                Mxc -= np.outer(Qn[:, j], Qn[:, j] @ Mx)
            num = Mxc.T @ v_res
            den = np.sqrt(np.maximum((Mxc * Mxc).sum(axis=0), 1e-12) * vv2n)
            return np.clip(num / den, -1, 1), Qn, vv2, vv2n

        # 观测: 自身梯度 / 控制对方后的梯度
        rho_self, _, _, _ = partial_rho(R)          # 全部 K 列(含 noCall)
        # 控制对方模块(秩空间)后再求自身梯度
        rho_adj = np.array([
            partial_rho(R[:, [0]], cov=R[:, [1]])[0][0],
            partial_rho(R[:, [1]], cov=R[:, [0]])[0][0],
        ])
        if extra:
            rho_nocall = float(rho_self[2])
        else:
            rho_nocall = np.nan

        # ---- MSR ----
        B = build_knn_graph(xy)
        lam, V, S0 = moran_spectrum(B)
        del B
        o = np.argsort(lam)
        lam_s = lam[o]
        V_s = np.ascontiguousarray(V[:, o])
        del V, lam
        C = V_s.T @ R
        w = V_s.T @ v_res
        U = V_s.T @ Q0
        del V_s
        nb = int(max(4, min(120, n // bin_target)))
        edges = np.linspace(0, n, nb + 1).astype(int)
        rng = np.random.default_rng(seed)
        rho_sur = np.empty((nsim, K), dtype=float)
        for r in range(nsim):
            Cs = np.empty_like(C)
            for b in range(nb):
                s, e = edges[b], edges[b + 1]
                if e - s > 1:
                    Cs[s:e] = C[s:e][rng.permutation(e - s)]
                else:
                    Cs[s:e] = C[s:e]
            num = w @ Cs
            Cu = Cs.T @ U
            den = np.sqrt(np.maximum((Cs * Cs).sum(axis=0) -
                                     (Cu ** 2).sum(axis=1), 1e-12) * vv)
            rho_sur[r] = np.clip(num / den, -1, 1)
        del C, Cs, Cu

        for j, nm in enumerate(names):
            mu = rho_sur[:, j].mean()
            sd = max(rho_sur[:, j].std(ddof=1), 1e-4)
            rows.append(dict(sample=p.stem, module=nm,
                             rho=rho_self[j], rho_adjusted=rho_adj[j] if j < 2 else np.nan,
                             u_msr=(rho_self[j] - mu) / sd,
                             sd_sur=sd, n_spot=n))
        print(f"  [{i}/{len(files)}] {p.stem} N={rho_self[0]:+.3f} "
              f"S={rho_self[1]:+.3f}  {time.time()-t0:.0f}s", flush=True)

    P = pd.DataFrame(rows)
    P.to_csv(PROC / "npp_split_spatial.csv", index=False)

    print("\n=== 模块梯度 (跨样本 meta) ===")
    out = []
    for mod, g in P.groupby("module"):
        u = g["u_msr"].dropna().to_numpy(float)
        ns = len(u)
        z = u.sum() / np.sqrt(ns) if ns else np.nan
        out.append(dict(module=mod, n_samples=ns,
                        mean_rho=g["rho"].mean(),
                        median_rho=g["rho"].median(),
                        consistency=(np.sign(g["rho"]) ==
                                     np.sign(g["rho"].mean())).mean(),
                        mean_rho_adjusted=g["rho_adjusted"].mean(),
                        stouffer_Z=z,
                        p_msr=2 * stats.norm.sf(abs(z)) if ns else np.nan))
    M = pd.DataFrame(out)
    M["fdr_msr"] = np.minimum(M["p_msr"] * len(M) / M["p_msr"].rank().values, 1.0)
    M.to_csv(TAB / "npp_split_spatial_meta.csv", index=False)
    print(M.round(4).to_string(index=False))
    return P, M


# ------------------------------------------------------------ bulk: METABRIC
def score_bulk(expr: pd.DataFrame, genes: list[str]) -> pd.Series:
    sub = expr[[g for g in genes if g in expr.columns]]
    Z = sub.sub(sub.mean(axis=0), axis=1).div(sub.std(axis=0).replace(0, np.nan), axis=1)
    return Z.mean(axis=1)


def run_bulk(n_boot=500, seed=42):
    expr = pd.read_csv("01_data/raw/metabric/metabric_mrna_targetgenes.csv", index_col=0)
    clin = pd.read_csv("01_data/processed/metabric_analysis_table.csv")
    print(f"[bulk] expr {expr.shape}, clin {clin.shape}")

    d = clin.set_index("PATIENT_ID").join(expr, how="inner") if "PATIENT_ID" in clin.columns \
        else clin.join(expr, how="inner")
    # 对齐索引
    if d.shape[0] == 0:
        d = clin.copy()
        for c in expr.columns:
            d[c] = expr[c].values
    scores = {}
    for nm, gl in MODULES.items():
        present = [g for g in gl if g in expr.columns]
        scores[nm] = score_bulk(expr, present)
        print(f"  {nm}: {len(present)}/{len(gl)} 基因可用于 bulk")
    Sc = pd.DataFrame(scores)
    d = d.join(Sc, how="left") if "Neural_intrinsic" not in d.columns else d

    # ---- 共线性 ----
    print("\n=== 模块共线性 (Spearman) ===")
    cor = []
    pairs = [("Neural_intrinsic", "Perineural_stroma"),
             ("Neural_intrinsic", "CAF_ECM"),
             ("Perineural_stroma", "CAF_ECM"),
             ("Neural_intrinsic", "NPPscore"),
             ("Perineural_stroma", "NPPscore")]
    for a, b in pairs:
        if a in d.columns and b in d.columns:
            m = d[[a, b]].dropna()
            r, p = stats.spearmanr(m[a], m[b])
            cor.append(dict(var1=a, var2=b, n=len(m), rho=r, p=p))
    C = pd.DataFrame(cor)
    print(C.round(4).to_string(index=False))
    C.to_csv(TAB / "npp_split_bulk_correlation.csv", index=False)

    # ---- Cox ----
    covs_base = ["AGE", "nodes_pos", "ER_pos", "chemo", "hormone"]
    rows = []
    for outcome, time_col, ev_col in [("RFS", "RFS_MONTHS", "RFS_EVENT"),
                                      ("OS", "OS_MONTHS", "OS_EVENT")]:
        for label, covs in [
            (f"{outcome}_N_unadj", ["Neural_intrinsic"]),
            (f"{outcome}_S_unadj", ["Perineural_stroma"]),
            (f"{outcome}_N_adj", ["Neural_intrinsic"] + covs_base),
            (f"{outcome}_S_adj", ["Perineural_stroma"] + covs_base),
            (f"{outcome}_mutual", ["Neural_intrinsic", "Perineural_stroma"] + covs_base),
            (f"{outcome}_mutual_CAF", ["Neural_intrinsic", "Perineural_stroma",
                                       "CAF_ECM"] + covs_base),
        ]:
            sub = d.dropna(subset=covs + [time_col, ev_col])
            sub = sub[sub[time_col] > 0]
            if len(sub) < 100:
                continue
            Z = pd.DataFrame({c: (sub[c] - sub[c].mean()) / sub[c].std(ddof=0)
                              for c in covs})
            t = sub[time_col].to_numpy(float)
            e = sub[ev_col].to_numpy(float)
            res = nsurv.cox(t, e, Z.to_numpy(float), names=covs)
            for r in res:
                rows.append(dict(outcome=outcome, model=label, variable=r.name,
                                 hr=r.hr, lo=r.lo, hi=r.hi, p=r.p,
                                 n=len(sub), events=int(e.sum())))
    R = pd.DataFrame(rows)
    R.to_csv(TAB / "npp_split_bulk_cox.csv", index=False)
    print("\n=== Cox ===")
    print(R[R.variable.isin(["Neural_intrinsic", "Perineural_stroma", "CAF_ECM"])]
          .round(4).to_string(index=False))

    # ---- bootstrap 稳定性 (mutual 模型) ----
    rng = np.random.default_rng(seed)
    draws = {}
    for outcome, time_col, ev_col in [("RFS", "RFS_MONTHS", "RFS_EVENT"),
                                      ("OS", "OS_MONTHS", "OS_EVENT")]:
        covs = ["Neural_intrinsic", "Perineural_stroma"] + covs_base
        sub = d.dropna(subset=covs + [time_col, ev_col])
        sub = sub[sub[time_col] > 0]
        Z0 = pd.DataFrame({c: (sub[c] - sub[c].mean()) / sub[c].std(ddof=0)
                           for c in covs})
        t = sub[time_col].to_numpy(float)
        e = sub[ev_col].to_numpy(float)
        acc = {c: [] for c in ["Neural_intrinsic", "Perineural_stroma"]}
        for b in range(n_boot):
            idx = rng.choice(len(sub), size=len(sub), replace=True)
            if e[idx].sum() < 10:
                continue
            try:
                res = nsurv.cox(t[idx], e[idx], Z0.to_numpy(float)[idx], names=covs)
            except Exception:
                continue
            for r in res:
                if r.name in acc:
                    acc[r.name].append(r.hr)
        for c, v in acc.items():
            v = np.asarray(v)
            if len(v) > 50:
                draws[f"{outcome}_{c}"] = v
                R.loc[len(R)] = dict(outcome=outcome, model=f"{outcome}_mutual_bootstrap",
                                     variable=c, hr=float(np.median(v)),
                                     lo=float(np.percentile(v, 2.5)),
                                     hi=float(np.percentile(v, 97.5)),
                                     p=float(np.mean(np.sign(v - 1) ==
                                                     np.sign(np.median(v) - 1))),
                                     n=len(sub), events=int(e.sum()))
    if draws:
        np.savez(TAB / "npp_split_bootstrap.npz", **draws)
    R.to_csv(TAB / "npp_split_bulk_cox.csv", index=False)
    print("\n=== bootstrap (mutual 模型, 500 次) ===")
    print(R[R.model.str.contains("bootstrap")].round(4).to_string(index=False))
    return d, R


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("stage", choices=["spatial", "bulk", "all"], default="all")
    a = ap.parse_args()
    if a.stage in ("spatial", "all"):
        run_spatial()
    if a.stage in ("bulk", "all"):
        run_bulk()

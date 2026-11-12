# -*- coding: utf-8 -*-
"""GSE213688 同平台外部复现 —— nsgrad 流程原样迁移。

目的
----
发现集 GSE210616(n=43, TNBC) 的结论能否在**独立 Visium 队列**上复现:
  1. 神经邻近 spot 是否存在, 占比量级是否一致 (发现集中位 2.3%)
  2. 固定 NPP(58 基因)迁移后是否仍随距离递减
  3. 基质模块 / 单个 myCAF 标记的 rho 方向与量级
  4. 衰减长度 1/b 的量级
  5. **去新推导**: 在本队列内独立重跑 NPP 推导, 与发现集 58 基因的重叠是否超随机

关键设计
--------
* 签名全部**固定**为发现集版本(不重新拟合), 这是真正的"迁移"检验;
* 输出目录独立 (01_data/processed/replication_gse213688/), 绝不覆盖主队列结果;
* 队列构成异质(CL-like / genomic-unstable TNBC / metaplastic), 因此同时报告
  "全队列" 与 "仅 TNBC 亚组" 两套元分析;
* 每样本独立, 失败不中断。
"""
from __future__ import annotations

import json
import sys
import traceback
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(r"D:/嗜神经浸润")
sys.path.insert(0, str(ROOT / "00_pipeline"))
sys.path.insert(0, str(ROOT / "02_scripts"))

from nsgrad import config as CFG  # noqa: E402
from nsgrad import io as nio  # noqa: E402
from nsgrad import nerve as nnerve  # noqa: E402
from nsgrad import gradient as ngrad  # noqa: E402
from nsgrad import signature as nsig  # noqa: E402
from run_step3a_nerve_program import sample_gene_distance_rho  # noqa: E402

RAW = CFG.DIR_RAW / "gse213688"
OUT = CFG.DIR_PROC / "replication_gse213688"
OUT.mkdir(parents=True, exist_ok=True)

# --------------------------------------------------------------- 固定签名
NPP = json.load(open(CFG.DIR_PROC / "npp_signature.json", encoding="utf-8"))
MOD = pd.read_csv(CFG.DIR_TAB / "supplementary" / "TableS8_module_definitions.csv")
MOD_GENES = {r["module"]: r["genes"].split(";") for _, r in MOD.iterrows()}

SIGNATURES: dict[str, list[str]] = dict(CFG.GRADIENT_SIGNATURES)
SIGNATURES.update({
    "NPP_total": NPP["up"] + NPP["down"],
    "NPP_up": NPP["up"],
    "NPP_down": NPP["down"],
    "NeuralModule": MOD_GENES["Neural_intrinsic"],
    "StromaModule": MOD_GENES["Perineural_stroma"],
})
for g in ["LRRC15", "SPARC", "POSTN", "ACTA2", "TNC", "COL1A1", "FN1", "LOXL2", "CCN2", "MMP13"]:
    SIGNATURES[f"g_{g}"] = [g]

GRAD_METRICS = ["NPP_total", "StromaModule", "NeuralModule", "CollagenIndex", "myCAFscore",
                "CD8density", "TLSscore", "g_LRRC15", "g_SPARC", "g_POSTN", "g_ACTA2", "g_TNC"]

# 发现集基准值(用于并排对照; 来自 n=43 主分析)
DISCOVERY = {
    "pct_nerve_median": 2.3,
    "NPP_total": -0.070,
    "StromaModule": -0.070,
    "NeuralModule": -0.122,
    "decay_len_um": 147.0,
}


def discover():
    out = []
    for d in sorted(RAW.glob("GSM*")):
        if not d.is_dir():
            continue
        h5 = list(d.glob("*filtered_feature_bc_matrix.h5"))
        pos = list(d.glob("*tissue_positions*"))
        part = list(d.glob("*.part"))
        if h5 and pos and not part:
            out.append((d.name, d))
    return out


def partial_spearman(y, d, tc, pm):
    """秩空间偏 Spearman: 从 y 秩与距离秩中剔除 log(UMI) 与 pct_mt 的秩。"""
    ok = np.isfinite(y) & np.isfinite(d)
    if ok.sum() < 50:
        return np.nan
    y, d, tc, pm = y[ok], d[ok], tc[ok], pm[ok]
    ry = stats.rankdata(y)
    rd = stats.rankdata(d)
    C = np.column_stack([stats.rankdata(np.log1p(tc)), stats.rankdata(pm)]).astype(float)
    def norm(v):
        v = v - v.mean()
        s = v.std()
        return v / s if s > 1e-12 else v * 0
    ry, rd = norm(ry), norm(rd)
    C = np.column_stack([norm(C[:, i]) for i in range(C.shape[1])])
    Q = np.zeros_like(C)
    for j in range(C.shape[1]):
        v = C[:, j].copy()
        for k in range(j):
            v -= Q[:, k] * (Q[:, k] @ v)
        nv = np.sqrt(v @ v)
        Q[:, j] = v / nv if nv > 1e-12 else 0.0
    def resid(A):
        for j in range(Q.shape[1]):
            if np.abs(Q[:, j]).max() < 1e-12:
                continue
            A = A - np.outer(Q[:, j], Q[:, j] @ A)
        return A
    ry_r = resid(ry[:, None])[:, 0]
    rd_r = resid(rd[:, None])[:, 0]
    den = np.sqrt((ry_r ** 2).sum() * (rd_r ** 2).sum())
    return float(np.clip((ry_r * rd_r).sum() / (den + 1e-12), -1, 1))


def process_one(name, d):
    ad = nio.read_visium_sample(d, sample_id=name)
    if ad.n_obs < 50:
        raise ValueError(f"too few spots: {ad.n_obs}")

    S, cover = nsig.score_signatures(ad, SIGNATURES)
    nv = nnerve.identify_nerve_spots(ad, q=CFG.QC["nerve_q_default"])
    mask = (nv["nerve_spatial"] == 1).values
    obs = ad.obs
    rec = {
        "sample": name,
        "n_spots_qc": int(ad.n_obs),
        "n_nerve_spots": int(mask.sum()),
        "pct_nerve_spots": 100.0 * float(mask.sum()) / max(int(ad.n_obs), 1),
    }
    if mask.sum() == 0:
        return rec, None, None, None

    gdf, dist_field = ngrad.fit_sample_gradients(
        ad, {m: S[m].values for m in GRAD_METRICS}, mask, metrics=GRAD_METRICS, verbose=False)
    for _, r in gdf.iterrows():
        m = r["metric"]
        rec[f"b_{m}"] = r["b"]
        rec[f"r2_{m}"] = r["r2"]
        rec[f"decaylen_{m}"] = r.get("decay_points_um", np.nan)

    prof = ngrad.ring_profile(ad, {m: S[m].values for m in GRAD_METRICS}, dist_field, nerve_mask=mask)
    prof.insert(0, "sample", name)
    gdf.insert(0, "sample", name)

    # 签名级偏相关(与发现集同一统计口径)
    dist = dist_field["distance_to_nearest_nerve"].values
    tc = obs["total_counts"].to_numpy(float)
    pm = obs["pct_mt"].to_numpy(float)
    for m in GRAD_METRICS:
        rec[f"rho_{m}"] = partial_spearman(np.asarray(S[m].values, float), dist, tc, pm)

    # 逐基因偏相关
    ad.obs = ad.obs.join(dist_field)
    ad.obs = ad.obs.join(nv)
    gene_df = sample_gene_distance_rho(ad)
    if gene_df is not None:
        gene_df.insert(0, "sample", name)
    return rec, gdf, prof, gene_df


def meta_rho(df, by):
    """按 rho 列做 Fisher-z 单样本 t 检验汇总。"""
    rows = []
    for key, sub in df.groupby(by):
        r = sub["rho"].dropna()
        n = len(r)
        if n < 3:
            rows.append({by: key, "n": n, "rho_mean": r.mean() if n else np.nan,
                         "t": np.nan, "p": np.nan})
            continue
        z = np.arctanh(np.clip(r, -0.999, 0.999))
        t, p = stats.ttest_1samp(z, 0)
        rows.append({by: key, "n": n, "rho_mean": float(np.tanh(z.mean())),
                     "rho_median": float(r.median()), "rho_sd": float(r.std()),
                     "t": float(t), "p": float(p),
                     "n_neg": int((r < 0).sum()), "n_pos": int((r > 0).sum())})
    out = pd.DataFrame(rows).sort_values("p")
    if len(out):
        out["fdr"] = np.minimum(1, out["p"] * len(out) / out["p"].rank())
    return out


def main():
    samples = discover()
    print(f"[repl] ready samples: {len(samples)}")
    if not samples:
        return
    recs, gdfs, profs, genes, fails = [], [], [], [], []
    for i, (name, d) in enumerate(samples, 1):
        try:
            r, g, p, gd = process_one(name, d)
            recs.append(r)
            if g is not None:
                gdfs.append(g); profs.append(p); genes.append(gd)
            print(f"  [{i}/{len(samples)}] {name} spots={r['n_spots_qc']} "
                  f"nerve={r['n_nerve_spots']} ({r['pct_nerve_spots']:.2f}%) "
                  f"rho_NPP={r.get('rho_NPP_total', float('nan')):.3f}", flush=True)
        except Exception as e:
            fails.append({"sample": name, "error": f"{type(e).__name__}: {e}",
                          "tb": traceback.format_exc()[-400:]})
            print(f"  [{i}/{len(samples)}] {name} FAILED {type(e).__name__}: {str(e)[:90]}", flush=True)

    pd.DataFrame(recs).to_csv(OUT / "sample_summary.csv", index=False)
    if gdfs:
        pd.concat(gdfs, ignore_index=True).to_csv(OUT / "gradient_params.csv", index=False)
    if profs:
        pd.concat(profs, ignore_index=True).to_csv(OUT / "ring_profiles.csv", index=False)
    if fails:
        pd.DataFrame(fails).to_csv(OUT / "failed.csv", index=False)

    if genes:
        G = pd.concat(genes, ignore_index=True)
        G.to_parquet(OUT / "gene_level_per_sample.parquet")
        # 每基因可测样本数与跨样本元分析
        meta = []
        for g, sub in G.groupby("gene"):
            r = sub["rho"].dropna()
            if len(r) < 3:
                continue
            z = np.arctanh(np.clip(r, -0.999, 0.999))
            t, p = stats.ttest_1samp(z, 0)
            meta.append({"gene": g, "n": len(r), "rho_mean": float(np.tanh(z.mean())),
                         "t": t, "p": p, "sign_consistency": float(max((r < 0).mean(), (r > 0).mean()))})
        M = pd.DataFrame(meta)
        M = M[M["n"] >= max(3, int(0.5 * len(genes)))]
        M["fdr"] = np.minimum(1, M["p"] * len(M) / M["p"].rank())
        M.sort_values("rho_mean").to_csv(OUT / "gene_meta.csv", index=False)

        # 去新推导的 NPP 与发现集的重叠
        sig = M[(M["p"] < 0.05) & (M["rho_mean"].abs() >= 0.05)]
        de_novo_up = set(sig[sig["rho_mean"] < 0].nlargest(150, "rho_mean", keep="all")["gene"]) if len(sig) else set()
        de_novo_up = set(sig[sig["rho_mean"] < 0]["gene"])
        de_novo_dn = set(sig[sig["rho_mean"] > 0]["gene"])
        disc = set(NPP["up"] + NPP["down"])
        hit = len(disc & (de_novo_up | de_novo_dn))
        p_hyper = stats.hypergeom.sf(hit - 1, len(M), len(de_novo_up | de_novo_dn), len(disc)) if hit else 1.0
        json.dump({
            "n_genes_tested": int(len(M)),
            "de_novo_n_up": len(de_novo_up), "de_novo_n_down": len(de_novo_dn),
            "overlap_with_discovery_58": hit,
            "hypergeom_p": float(p_hyper),
            "overlap_genes": sorted(disc & (de_novo_up | de_novo_dn)),
        }, open(OUT / "npp_overlap.json", "w", encoding="utf-8"), indent=1)

    # 签名级元分析
    R = pd.DataFrame(recs)
    rho_cols = [c for c in R.columns if c.startswith("rho_")]
    if rho_cols:
        long = R.melt(id_vars=["sample"], value_vars=rho_cols, var_name="metric", value_name="rho")
        long["metric"] = long["metric"].str.replace("^rho_", "", regex=True)
        long.dropna().to_csv(OUT / "signature_rho_per_sample.csv", index=False)
        meta_rho(long, "metric").to_csv(OUT / "signature_meta.csv", index=False)
    print(f"[repl] done -> {OUT}")


if __name__ == "__main__":
    main()

# -*- coding: utf-8 -*-
"""Figure 3 —— 神经近端转录程序 (NPP) 的表征与机制。

panels
------
3A 基因水平元分析火山图 (mean rho vs -log10 p)
3B Top-30 基因效应量条形图
3C NPP-up 基因的功能类别构成 + 超几何富集
3D 按神经距离的环带 NPP 剖面 (跨样本, 配对)
3E 示例样本的 NPP 空间分布图 + 神经源点
3F 每样本 NPP~距离 的 Spearman rho 箱线图 (Nervesignal 作阳性对照)
3G 各神经信号轴随距离的环带剖面
3H 两个空间亚型之间的 NPP 差异
"""
from __future__ import annotations
import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import stats

sys.path.insert(0, r"D:/嗜神经浸润/00_pipeline")
from nsgrad import config as CFG        # noqa: E402
from nsgrad import viz as nviz          # noqa: E402

PROC = CFG.DIR_PROC
FIG = CFG.DIR_FIG
FIG.mkdir(parents=True, exist_ok=True)
DPI = 300

AXES = {
    "Sensory_nerve": ["NGF", "UCHL1", "TRPV1", "CALCA", "S100B", "SCN9A", "SCN10A"],
    "Sympathetic": ["TH", "DBH", "CHGA", "SLC6A2"],
    "Parasympathetic": ["SLC18A3", "CHAT", "SLC5A7"],
    "CGRP_axis": ["CALCA", "CALCRL", "RAMP1", "CALCR", "RAMP2", "RAMP3"],
    "NGF_TrkA_axis": ["NGF", "NTRK1", "NGFR", "SORT1", "NTRK2", "NTRK3"],
    "GABA_BDNF_axis": ["GAD1", "GAD2", "ABAT", "SLC6A1", "GABBR1", "GABBR2", "BDNF", "NTF3"],
    "Semaphorin_axis": ["SEMA3A", "SEMA3C", "SEMA3D", "NRP1", "NRP2", "PLXNA1", "PLXND1"],
    "CAF_ECM": ["COL1A1", "COL1A2", "COL3A1", "COL5A1", "FN1", "POSTN", "ACTA2",
                "TAGLN", "FAP", "PDGFRB", "SPARC", "THBS2", "CTHRC1", "LOX"],
    "CD8_cytotox": ["CD8A", "CD8B", "GZMA", "GZMB", "PRF1", "NKG7", "IFNG"],
}
# 手工功能类别(用于 3C 的富集检验; 离线环境下无法跑 GO, 用 curated 列表代替并如实标注)
CATEGORIES = {
    "ECM / collagen": ["COL1A1", "COL1A2", "COL3A1", "COL4A1", "COL4A2", "COL5A1", "COL5A2",
                       "COL6A1", "COL6A2", "COL6A3", "COL11A1", "COL12A1", "COL14A1", "FN1",
                       "POSTN", "SPARC", "THBS1", "THBS2", "TNC", "VCAN", "LUM", "DCN",
                       "FBN1", "NID1", "NID2", "HSPG2", "LAMA4", "LAMB1", "LAMC1"],
    "CAF / myofibroblast": ["ACTA2", "TAGLN", "MYL9", "CNN1", "CALD1", "FAP", "PDGFRB",
                            "PDPN", "S100A4", "CTHRC1", "INHBA", "LRRC15", "SULF1",
                            "FSTL1", "CCN2", "TGFB1", "MMP2", "MMP9", "MMP11", "MMP14",
                            "LOX", "LOXL2", "SERPINE1", "TIMP1"],
    "Neuronal / glial": ["UCHL1", "TUBB3", "MAP2", "NEFL", "NEFM", "PLP1", "MPZ", "SYP",
                         "S100B", "SCN9A", "SCN10A", "TRPV1", "CALCA", "NGF", "MAP1B",
                         "TUBA1A", "TUBB2B", "STMN2", "UCHL1"],
    "Antigen presentation / IFN": ["HLA-A", "HLA-B", "HLA-C", "HLA-DRA", "HLA-DRB1",
                                   "B2M", "TAP1", "TAP2", "NLRC5", "PSMB8", "PSMB9",
                                   "STAT1", "IRF1", "IFI27", "IFI6", "ISG15", "GBP1"],
    "Ribosome / translation": [f"RPS{i}" for i in range(2, 30)] +
                              [f"RPL{i}" for i in range(3, 40)] + ["FAU", "EEF1A1", "RACK1"],
}


def load_sig():
    with open(PROC / "npp_signature.json", encoding="utf-8") as f:
        return json.load(f)


def ring_order(labels):
    """按环带下界排序; gt400 / out 类环永远排最后。"""
    import re as _re
    def key(x):
        s = str(x).lower()
        if "gt" in s or "out" in s or ">" in s:
            return 10 ** 6
        nums = _re.findall(r"\d+", s)
        return int(nums[1]) if len(nums) >= 2 else (int(nums[0]) if nums else 0)
    return sorted(set(labels), key=key)


def fig3a(G):
    fig, ax = plt.subplots(figsize=(6.6, 5.0))
    x = G["mean_rho"].to_numpy()
    y = -np.log10(np.clip(G["p"].to_numpy(), 1e-300, 1))
    sig = (G["p"].to_numpy() < 0.05) & (np.abs(x) >= 0.05)
    ax.scatter(x[~sig], y[~sig], s=6, c="#c9ced6", alpha=0.5, linewidths=0)
    ax.scatter(x[sig], y[sig], s=14, c="#c0392b", alpha=0.85, linewidths=0)
    top = G.reindex(np.argsort(np.abs(G["mean_rho"].to_numpy()))[::-1][:14])
    for g, r in top.iterrows():
        ax.annotate(g, (r["mean_rho"], -np.log10(max(r["p"], 1e-300))),
                    fontsize=7, xytext=(3, 3), textcoords="offset points")
    ax.axvline(0, color="#888", lw=0.7, ls="--")
    ax.axhline(-np.log10(0.05), color="#2c3e50", lw=0.7, ls=":")
    ax.axvline(-0.05, color="#7f8c8d", lw=0.6, ls=":")
    ax.axvline(0.05, color="#7f8c8d", lw=0.6, ls=":")
    ax.set_xlabel("Mean Spearman ρ (gene vs. distance-to-nearest-nerve)")
    ax.set_ylabel("−log10 P (meta-analysis across samples)")
    ax.set_title("A  Gene-level nerve-distance association", loc="left", fontweight="bold")
    ax.text(0.02, 0.97, "ρ < 0: enriched in the perineural niche",
            transform=ax.transAxes, fontsize=7.5, va="top", color="#555")
    fig.tight_layout()
    fig.savefig(FIG / "Fig3A_volcano.png", dpi=DPI)
    fig.savefig(FIG / "Fig3A_volcano.pdf")
    plt.close(fig)


def fig3b(G, sig):
    up = G.loc[[g for g in sig["up"] if g in G.index]].sort_values("mean_rho")
    dn = G.loc[[g for g in sig["down"] if g in G.index]].sort_values("mean_rho", ascending=False)
    show = pd.concat([up.head(22), dn.head(8)])
    fig, ax = plt.subplots(figsize=(5.4, 6.4))
    col = ["#c0392b" if v < 0 else "#2980b9" for v in show["mean_rho"]]
    ax.barh(range(len(show)), show["mean_rho"].to_numpy(), color=col)
    ax.set_yticks(range(len(show)))
    ax.set_yticklabels(show.index, fontsize=7.5)
    ax.invert_yaxis()
    ax.set_xlabel("Mean Spearman ρ vs. nerve distance")
    ax.set_title("B  Nerve-proximity program genes", loc="left", fontweight="bold")
    ax.axvline(0, color="#333", lw=0.8)
    fig.tight_layout()
    fig.savefig(FIG / "Fig3B_top_genes.png", dpi=DPI)
    fig.savefig(FIG / "Fig3B_top_genes.pdf")
    plt.close(fig)


def fig3c(G, sig):
    bg = set(G.index)
    rows = []
    for cat, genes in CATEGORIES.items():
        gs = set(genes) & bg
        hits = set(sig["up"]) & gs
        if len(gs) < 5:
            continue
        # 超几何: 从 |bg| 中抽 |up| 个, 命中该类别 k 个
        k = len(hits); K = len(gs); N = len(bg); n = len(sig["up"])
        p = stats.hypergeom.sf(k - 1, N, K, n)
        exp = n * K / N
        rows.append(dict(category=cat, k=k, K=K, expected=round(exp, 1),
                         fold=round(k / exp, 2) if exp > 0 else np.nan,
                         p=p, fdr=np.nan))
    if not rows:
        return None
    R = pd.DataFrame(rows)
    ps = R["p"].to_numpy()
    order = np.argsort(ps); n = len(ps)
    q = np.empty(n)
    q[order] = np.minimum.accumulate((ps[order] * n / (np.arange(n) + 1))[::-1])[::-1]
    R["fdr"] = np.clip(q, 0, 1)
    R = R.sort_values("fold", ascending=False)

    fig, ax = plt.subplots(figsize=(6.0, 3.6))
    colors = ["#c0392b" if f < 0.05 else "#95a5a6" for f in R["fdr"]]
    ax.bar(range(len(R)), R["fold"], color=colors)
    ax.set_xticks(range(len(R)))
    ax.set_xticklabels([c.replace(" / ", "\n") for c in R["category"]], fontsize=7.5)
    ax.axhline(1, color="#333", lw=0.8, ls="--")
    ax.set_ylabel("Observed / expected")
    ax.set_title("C  Functional composition of NPP-up genes", loc="left", fontweight="bold")
    for i, r in enumerate(R.itertuples()):
        ax.text(i, r.fold + 0.06, f"{r.k}/{r.K}\nFDR={r.fdr:.2g}", ha="center", fontsize=6.5)
    ax.set_ylim(0, max(R["fold"].max() * 1.45, 2))
    fig.tight_layout()
    fig.savefig(FIG / "Fig3C_category_enrichment.png", dpi=DPI)
    fig.savefig(FIG / "Fig3C_category_enrichment.pdf")
    plt.close(fig)
    return R


def compute_spot_scores():
    """对每个样本计算 spot 级 NPP 与各神经轴得分 + 环带剖面。"""
    sig = load_sig()
    files = sorted((PROC / "spot_level").glob("*.h5ad"))
    ring_rows, rho_rows, demo = [], [], None
    for p in files:
        ad = sc.read_h5ad(p)
        obs = ad.obs
        d = obs["distance_to_nearest_nerve"].to_numpy(float)
        X = ad.X
        X = X.toarray() if hasattr(X, "toarray") else np.asarray(X)
        genes = np.asarray(ad.var_names)
        gi = {g: i for i, g in enumerate(genes)}
        Z = (X - X.mean(0)) / (X.std(0) + 1e-9)     # spot 间 z-score
        up = [g for g in sig["up"] if g in gi]
        dn = [g for g in sig["down"] if g in gi]
        npp = Z[:, [gi[g] for g in up]].mean(1) - (Z[:, [gi[g] for g in dn]].mean(1)
                                                   if dn else 0.0)
        scores = {"NPP": npp, "Nervesignal": obs["Nervesignal"].to_numpy(float)}
        for k, gl in AXES.items():
            idx = [gi[g] for g in gl if g in gi]
            if idx:
                scores[k] = Z[:, idx].mean(1)
        S = pd.DataFrame(scores)
        S["ring"] = obs["ring_id"].to_numpy()
        S["dist"] = d
        S["sample"] = p.stem
        nonerve = obs["nerve_spatial"].to_numpy() == 0
        for k in S.columns[:-3]:
            try:
                r = stats.spearmanr(d[nonerve], S[k].to_numpy()[nonerve])
                rho_rows.append(dict(sample=p.stem, metric=k,
                                     rho=float(r.statistic), p=float(r.pvalue)))
            except Exception:
                pass
        for rlab, gg in S.groupby("ring"):
            row = dict(sample=p.stem, ring=rlab, n=len(gg),
                       mean_dist=float(np.median(gg["dist"])))
            for k in S.columns[:-3]:
                row[k] = float(gg[k].mean())
            ring_rows.append(row)
        if demo is None or (len(demo[1]) < ad.n_obs and
                            obs["nerve_spatial"].sum() > 20):
            demo = (p.stem, ad.obs.copy(), S.copy())
        del ad, X, Z
    return (pd.DataFrame(ring_rows), pd.DataFrame(rho_rows), demo, sig)


def fig3d(R):
    order = ring_order(R["ring"])
    fig, ax = plt.subplots(figsize=(5.6, 4.2))
    piv = R.pivot_table(index="sample", columns="ring", values="NPP")
    piv = piv[[c for c in order if c in piv.columns]]
    means = piv.mean(axis=0)
    xs = np.arange(len(means))
    for _, row in piv.iterrows():
        ax.plot(xs, row.to_numpy(), color="#bdc3c7", lw=0.7, alpha=0.8)
    ax.plot(xs, means.to_numpy(), color="#c0392b", lw=2.2, marker="o")
    ax.set_xticks(xs)
    ax.set_xticklabels([c.replace("r", "").replace("_", "–").replace("gt", ">") + " µm" for c in means.index], fontsize=7)
    ax.set_xlabel("Distance ring from nearest nerve source")
    ax.set_ylabel("NPP score (spot-level z)")
    ax.set_title("D  NPP declines with nerve distance", loc="left", fontweight="bold")
    ax.axhline(0, color="#888", lw=0.6, ls="--")
    # 配对检验: 最近环 vs 最远环
    a, b = piv.columns[0], piv.columns[-1]
    try:
        t = stats.ttest_rel(piv[a], piv[b])
        ax.text(0.98, 0.05, f"paired t, P={t.pvalue:.2g}", transform=ax.transAxes,
                ha="right", fontsize=7.5)
    except Exception:
        pass
    fig.tight_layout()
    fig.savefig(FIG / "Fig3D_ring_profile.png", dpi=DPI)
    fig.savefig(FIG / "Fig3D_ring_profile.pdf")
    plt.close(fig)


def fig3e(demo):
    name, obs, S = demo
    fig, ax = plt.subplots(figsize=(5.0, 4.6))
    sc_ = ax.scatter(obs["x_um"], obs["y_um"], c=S["NPP"], cmap="RdBu_r",
                     s=9, marker="s", linewidths=0)
    nv = obs["nerve_spatial"].to_numpy() == 1
    ax.scatter(obs["x_um"][nv], obs["y_um"][nv], facecolors="none",
               edgecolors="k", s=42, linewidths=0.7, label="nerve source")
    ax.invert_yaxis()
    ax.set_aspect("equal")
    ax.legend(fontsize=7, loc="upper right", frameon=True)
    cb = plt.colorbar(sc_, ax=ax, fraction=0.035, pad=0.02)
    cb.set_label("NPP score", fontsize=8)
    ax.set_xticks([]); ax.set_yticks([])
    ax.set_title(f"E  NPP spatial map ({name})", loc="left", fontweight="bold")
    fig.tight_layout()
    fig.savefig(FIG / "Fig3E_npp_map.png", dpi=DPI)
    fig.savefig(FIG / "Fig3E_npp_map.pdf")
    plt.close(fig)


def fig3f(rho):
    keep = [m for m in ["Nervesignal", "NPP"] + list(AXES.keys())
            if m in set(rho["metric"])]
    data = [rho.loc[rho["metric"] == m, "rho"].to_numpy() for m in keep]
    fig, ax = plt.subplots(figsize=(6.8, 4.0))
    bp = ax.boxplot(data, vert=True, patch_artist=True, widths=0.6, showfliers=False)
    for b in bp["boxes"]:
        b.set_facecolor("#dfe6e9"); b.set_edgecolor("#2d3436")
    for i, m in enumerate(keep, start=1):
        v = rho.loc[rho["metric"] == m, "rho"].to_numpy()
        rng = np.random.default_rng(i)
        ax.scatter(rng.normal(i, 0.06, len(v)), v, s=7, c="#2d3436", alpha=0.5, zorder=3)
    ax.axhline(0, color="#c0392b", lw=0.9, ls="--")
    ax.set_xticks(range(1, len(keep) + 1))
    ax.set_xticklabels(keep, rotation=35, ha="right", fontsize=8)
    ax.set_ylabel("Spearman ρ vs. nerve distance")
    ax.set_title("F  Per-sample gradient strength", loc="left", fontweight="bold")
    fig.tight_layout()
    fig.savefig(FIG / "Fig3F_gradient_rho_boxplot.png", dpi=DPI)
    fig.savefig(FIG / "Fig3F_gradient_rho_boxplot.pdf")
    plt.close(fig)


def fig3g(R):
    metrics = ["CGRP_axis", "NGF_TrkA_axis", "GABA_BDNF_axis",
               "Semaphorin_axis", "CAF_ECM", "CD8_cytotox"]
    metrics = [m for m in metrics if m in R.columns]
    order = ring_order(R["ring"])
    ncol, nrow = 3, 2
    fig, axes = plt.subplots(nrow, ncol, figsize=(10.5, 5.8), sharex=True)
    for i, m in enumerate(metrics):
        ax = axes.ravel()[i]
        piv = R.pivot_table(index="sample", columns="ring", values=m)
        piv = piv[[c for c in order if c in piv.columns]]
        z = (piv - piv.mean(axis=1).values[:, None])
        xs = np.arange(piv.shape[1])
        for _, row in z.iterrows():
            ax.plot(xs, row.to_numpy(), color="#bdc3c7", lw=0.6, alpha=0.8)
        ax.plot(xs, z.mean(axis=0).to_numpy(), color="#8e44ad", lw=2.0, marker="o", ms=3)
        ax.axhline(0, color="#888", lw=0.6, ls="--")
        ax.set_title(m.replace("_", " "), fontsize=9)
        if i >= 3:
            ax.set_xticks(xs)
            ax.set_xticklabels([c.replace("_", "\n") for c in piv.columns], fontsize=6.5)
    for j in range(len(metrics), ncol * nrow):
        axes.ravel()[j].axis("off")
    fig.suptitle("G  Neuro-signalling axes across distance rings", fontweight="bold", x=0.01,
                 ha="left", y=1.01)
    fig.tight_layout()
    fig.savefig(FIG / "Fig3G_axis_ring_profiles.png", dpi=DPI)
    fig.savefig(FIG / "Fig3G_axis_ring_profiles.pdf")
    plt.close(fig)


def fig3h():
    sub = pd.read_csv(PROC / "sample_subtypes.csv")
    ss = pd.read_csv(PROC / "sample_summary.csv")
    m = sub.merge(ss, left_on="sample", right_on="sample", how="left")
    fig, axes = plt.subplots(1, 2, figsize=(7.6, 3.7))
    feat = "pct_nerve_spots"
    if feat in m.columns and "subtype" in m.columns:
        groups = [m.loc[m["subtype"] == s, feat].dropna().to_numpy()
                  for s in sorted(m["subtype"].unique())]
        axes[0].boxplot(groups, tick_labels=sorted(m["subtype"].unique()),
                        patch_artist=True, widths=0.55)
        for b in axes[0].boxplots if False else []:
            pass
        for bb in axes[0].artists if hasattr(axes[0], "artists") else []:
            bb.set_facecolor("#dfe6e9")
        axes[0].set_ylabel("Nerve spot fraction (%)")
        try:
            g1, g2 = groups
            u = stats.mannwhitneyu(g1, g2)
            axes[0].set_title(f"Nerve density  (P={u.pvalue:.3g})", fontsize=9)
        except Exception:
            axes[0].set_title("Nerve density", fontsize=9)
    if "b_Nervesignal" in m.columns:
        b = m["b_Nervesignal"].astype(float)
        b = b[(b > 1 / 2000 * 1.001) & (b < 1 / 25 * 0.999)]
        g = [b[m.loc[b.index, "subtype"] == s].dropna().to_numpy()
             for s in sorted(m["subtype"].unique())]
        axes[1].boxplot(g, tick_labels=sorted(m["subtype"].unique()),
                        patch_artist=True, widths=0.55)
        for bb in axes[1].artists if hasattr(axes[1], "artists") else []:
            bb.set_facecolor("#dfe6e9")
        axes[1].set_ylabel("Decay constant b (1/µm)")
        axes[1].set_title("Nerve-signal decay constant", fontsize=9)
    fig.suptitle("H  Spatial subtypes", fontweight="bold", x=0.01, ha="left", y=1.02)
    fig.tight_layout()
    fig.savefig(FIG / "Fig3H_subtype_comparison.png", dpi=DPI)
    fig.savefig(FIG / "Fig3H_subtype_comparison.pdf")
    plt.close(fig)


def main():
    G = pd.read_csv(PROC / "nerve_proximity_program.csv", index_col=0)
    sig = load_sig()
    print("[Fig3] panels A-C")
    fig3a(G); print("  A")
    fig3b(G, sig); print("  B")
    R = fig3c(G, sig); print("  C")
    if R is not None:
        R.to_csv(CFG.DIR_TAB / "fig3c_category_enrichment.csv", index=False)
    print("[Fig3] computing spot-level scores ...")
    Rring, Rrho, demo, _ = compute_spot_scores()
    Rring.to_csv(PROC / "npp_ring_profiles.csv", index=False)
    Rrho.to_csv(PROC / "npp_gradient_rho.csv", index=False)
    print(f"  ring rows={len(Rring)} rho rows={len(Rrho)}")
    fig3d(Rring); print("  D")
    fig3e(demo); print("  E")
    fig3f(Rrho); print("  F")
    fig3g(Rring); print("  G")
    fig3h(); print("  H")
    print("[Fig3] done")


if __name__ == "__main__":
    main()

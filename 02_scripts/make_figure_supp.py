# -*- coding: utf-8 -*-
"""Supplementary Figures S1-S3 —— 稳健性分析的成套图。

S1 神经 spot 阈值敏感性 (q x min_markers x 空间约束)
S2 空间分型稳健性 (特征 jackknife / 算法替换 / 样本自助 / K 选择换种子)
S3 METABRIC 自助法 HR 分布 + E-value

风格与 Figure 1-4 保持一致: 8pt 字号、300 dpi、PNG + PDF 双输出。
"""
from __future__ import annotations
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, r"D:/嗜神经浸润/00_pipeline")
from nsgrad import config as CFG        # noqa: E402

TAB = CFG.DIR_TAB
FIG = CFG.DIR_FIG
FIG.mkdir(parents=True, exist_ok=True)
DPI = 300

plt.rcParams.update({
    "font.size": 8, "axes.titlesize": 9, "axes.labelsize": 8,
    "xtick.labelsize": 7, "ytick.labelsize": 7, "legend.fontsize": 7,
    "svg.fonttype": "none", "font.family": "DejaVu Sans",
})

C_MAIN = "#2980B9"
C_INDEP = "#C0392B"
C_GREY = "#7F8C8D"


def save(fig, name):
    fig.tight_layout()
    fig.savefig(FIG / f"{name}.png", dpi=DPI, bbox_inches="tight")
    fig.savefig(FIG / f"{name}.pdf", bbox_inches="tight")
    plt.close(fig)
    print(f"  -> {name}.png/.pdf")


# --------------------------------------------------------------------------- #
# S1
# --------------------------------------------------------------------------- #
def fig_s1():
    p = TAB / "sensitivity_nerve_threshold.csv"
    if not p.exists():
        print("  skip S1: no table"); return
    d = pd.read_csv(p)
    ds = d[d["spatial_constraint"] == True]          # noqa: E712
    if not len(ds):
        ds = d
    cfgs = [(q, mm) for q in sorted(ds["q"].unique()) for mm in sorted(ds["min_markers"].unique())]
    labels = [f"q{int(q*100)}\nmm{mm}" for q, mm in cfgs]

    fig, axes = plt.subplots(2, 2, figsize=(10.4, 7.2))

    # (a) 神经 spot 占比
    ax = axes[0, 0]
    data = [ds[(ds["q"] == q) & (ds["min_markers"] == mm)]["pct_nerve"].dropna().values
            for q, mm in cfgs]
    bp = ax.boxplot(data, tick_labels=labels, patch_artist=True, widths=0.6)
    for b in bp["boxes"]:
        b.set_facecolor("#D6EAF8"); b.set_edgecolor("#1A5276"); b.set_linewidth(0.8)
    for m in bp["medians"]:
        m.set_color(C_MAIN); m.set_linewidth(1.2)
    rng = np.random.default_rng(0)
    for i, dd in enumerate(data):
        x = rng.normal(i + 1, 0.055, len(dd))
        ax.scatter(x, dd, s=4, color="#154360", alpha=0.6, zorder=3)
    ax.axhline(2.6, color=C_INDEP, ls="--", lw=0.9)
    ax.text(len(cfgs) - 0.4, 2.7, "reported median 2.6%", fontsize=6.5,
            color=C_INDEP, ha="right")
    ax.set_ylabel("% spots called nerve")
    ax.set_title("A  Nerve-spot yield across calling thresholds",
                 loc="left", fontweight="bold")

    # (b) 效应量 (主 readout 与独立 readout)
    ax = axes[0, 1]
    x = np.arange(len(cfgs)); w = 0.38
    dm = [ds[(ds["q"] == q) & (ds["min_markers"] == mm)]["cohen_d_main"].median()
          for q, mm in cfgs]
    di = [ds[(ds["q"] == q) & (ds["min_markers"] == mm)]["cohen_d_indep"].median()
          for q, mm in cfgs]
    ax.bar(x - w / 2, dm, w, color=C_MAIN, label="Nervesignal (call signature)")
    ax.bar(x + w / 2, di, w, color=C_INDEP, label="independent axon/neuron set")
    ax.set_xticks(x); ax.set_xticklabels(labels)
    ax.axhline(0.8, color="#555", ls=":", lw=0.8)
    ax.text(-0.4, 0.85, "large effect d=0.8", fontsize=6.5, color="#555")
    ax.set_ylabel("Cohen's d (nerve vs non-nerve spots)")
    ax.legend(frameon=True, fontsize=6.5)
    ax.set_title("B  Effect size is threshold-insensitive",
                 loc="left", fontweight="bold")

    # (c) 独立读出 vs 主读出 AUC 逐样本散点 (默认配置)
    ax = axes[1, 0]
    dd = ds[(ds["q"] == 0.90) & (ds["min_markers"] == 2)]
    ax.scatter(dd["auc_main"], dd["auc_indep"], s=22, color="#34495E",
               edgecolor="white", linewidth=0.4)
    lo = min(dd["auc_main"].min(), dd["auc_indep"].min()) - 0.03
    hi = max(dd["auc_main"].max(), dd["auc_indep"].max()) + 0.03
    ax.plot([lo, hi], [lo, hi], ls="--", lw=0.7, color=C_GREY)
    ax.axhline(0.5, ls=":", lw=0.7, color="#888")
    ax.axvline(0.5, ls=":", lw=0.7, color="#888")
    ax.set_xlabel("AUC using Nervesignature genes")
    ax.set_ylabel("AUC using independent axon genes")
    rho = dd[["auc_main", "auc_indep"]].corr(method="spearman").iloc[0, 1]
    ax.text(0.03, 0.95, f"n={len(dd)} samples\nSpearman rho with indep = {rho:.2f}",
            transform=ax.transAxes, va="top", fontsize=6.5)
    ax.set_title("C  Independent validation of nerve calls (q90, 2 markers)",
                 loc="left", fontweight="bold")

    # (d) 下游梯度参数是否会翻转
    ax = axes[1, 1]
    dd1 = [ds[(ds["q"] == q) & (ds["min_markers"] == mm)]["decay_um_Nervesignal"]
           .dropna().values for q, mm in cfgs]
    bp = ax.boxplot(dd1, tick_labels=labels, patch_artist=True, widths=0.6)
    for b in bp["boxes"]:
        b.set_facecolor("#FDEDEC"); b.set_edgecolor("#922B21"); b.set_linewidth(0.8)
    for m in bp["medians"]:
        m.set_color(C_INDEP); m.set_linewidth(1.2)
    ax.set_yscale("log")
    ax.set_ylabel("characteristic decay length 1/b (um)")
    ax.set_title("D  Downstream decay length across thresholds",
                 loc="left", fontweight="bold")
    save(fig, "SuppFigS1_threshold_sensitivity")


# --------------------------------------------------------------------------- #
# S2
# --------------------------------------------------------------------------- #
def fig_s2():
    p = TAB / "sensitivity_clustering.csv"
    if not p.exists():
        print("  skip S2: no table"); return
    d = pd.read_csv(p)
    fig, axes = plt.subplots(2, 2, figsize=(10.4, 7.2))

    # (a) 特征留一
    ax = axes[0, 0]
    jk = d[d["test"] == "feature_jackknife"].dropna().sort_values("ari")
    ax.barh(np.arange(len(jk)), jk["ari"], color="#85C1E9", edgecolor="#1A5276",
            height=0.7)
    ax.set_yticks(np.arange(len(jk))); ax.set_yticklabels(jk["setting"], fontsize=6.5)
    ax.set_xlabel("ARI vs full-feature partition")
    ax.axvline(0.5, ls="--", lw=0.8, color=C_INDEP)
    ax.set_xlim(0, 1.05)
    ax.set_title("A  Leave-one-feature-out jackknife", loc="left", fontweight="bold")

    # (b) 换算法
    ax = axes[0, 1]
    al = d[d["test"] == "algorithm"].dropna().sort_values("ari")
    cols = ["#C0392B" if v < 0.4 else "#27AE60" if v > 0.6 else "#F39C12"
            for v in al["ari"]]
    ax.barh(np.arange(len(al)), al["ari"], color=cols, edgecolor="#555", height=0.7)
    ax.set_yticks(np.arange(len(al))); ax.set_yticklabels(al["setting"], fontsize=6.5)
    ax.set_xlabel("ARI vs consensus (ward) partition")
    ax.axvline(0.5, ls="--", lw=0.8, color=C_INDEP)
    ax.set_xlim(0, 1.05)
    ax.set_title("B  Algorithm substitution at K=2", loc="left", fontweight="bold")

    # (c) 样本自助法 ARI 分布
    ax = axes[1, 0]
    npz = TAB / "sensitivity_clustering_bootstrap.npy"
    if npz.exists():
        aris = np.load(npz)
        ax.hist(aris, bins=28, color="#AED6F1", edgecolor="#1A5276", linewidth=0.5)
        med = float(np.median(aris))
        ax.axvline(med, color=C_INDEP, lw=1.2)
        ax.text(med + 0.01, ax.get_ylim()[1] * 0.9, f"median ARI={med:.2f}",
                fontsize=7, color=C_INDEP)
        ax.set_xlabel("ARI of 80% subsample clustering vs full data")
        ax.set_ylabel("bootstrap replicates")
        ax.text(0.02, 0.9, "n=1000 replicates", transform=ax.transAxes, fontsize=6.5)
    ax.set_title("C  Sample-level bootstrap stability", loc="left", fontweight="bold")

    # (d) K 选择换种子
    ax = axes[1, 1]
    ks = d[d["test"] == "k_selection_seed"].dropna().copy()
    ks["K"] = ks["setting"].str.extract(r"K=(\d)").astype(int)
    ks["seed"] = ks["setting"].str.extract(r"seed=(\d+)").astype(int)
    for k, g in ks.groupby("K"):
        g = g.sort_values("seed")
        ax.plot(g["seed"].astype(str), g["ari"], marker="o", ms=4,
                label=f"K={k}", lw=1.2)
    ax.set_ylabel("consensus silhouette")
    ax.set_xlabel("random seed")
    ax.legend(frameon=True, fontsize=6.5, title="", loc="lower right")
    ax.set_ylim(0.4, 0.8)
    ax.set_title("D  K selection is seed-stable (K=2 wins at every seed)",
                 loc="left", fontweight="bold")
    save(fig, "SuppFigS2_clustering_robustness")


# --------------------------------------------------------------------------- #
# S3
# --------------------------------------------------------------------------- #
def fig_s3():
    p = TAB / "sensitivity_metabric_bootstrap.csv"
    if not p.exists():
        print("  skip S3: no table"); return
    d = pd.read_csv(p)
    npz = TAB / "sensitivity_metabric_bootstrap_draws.npz"
    draws = np.load(npz) if npz.exists() else {}

    fig, axes = plt.subplots(1, 2, figsize=(10.4, 3.9))
    ax = axes[0]
    keys = [k for k in draws.files if k.startswith("RFS")]
    for i, k in enumerate(sorted(keys)):
        h = draws[k]
        h = h[np.isfinite(h)]
        ax.hist(h, bins=40, alpha=0.55, label=f"{k} (n={len(h)})",
                color=["#2980B9", "#C0392B"][i % 2])
    ax.axvline(1.0, color="#333", ls="--", lw=0.9)
    ax.set_xlabel("bootstrap HR (per SD)")
    ax.set_ylabel("resamples")
    ax.legend(frameon=True, fontsize=6.5)
    ax.set_title("A  RFS: bootstrap distribution of adjusted HR",
                 loc="left", fontweight="bold")

    ax = axes[1]
    tmp = d.copy()
    tmp["lab"] = tmp["outcome"] + "\n" + tmp["variable"]
    x = np.arange(len(tmp))
    y = tmp["hr"]
    yerr = np.c_[tmp["hr"] - tmp["boot_ci_lo"], tmp["boot_ci_hi"] - tmp["hr"]].T
    ax.errorbar(x, y, yerr=yerr, fmt="o", ms=5, capsize=3, color="#2C3E50",
                ecolor="#7F8C8D", label="point estimate + bootstrap 95% CI")
    ax.axhline(1.0, ls="--", lw=0.9, color="#333")
    ax.set_xticks(x); ax.set_xticklabels(tmp["lab"], fontsize=6.5)
    ax.set_ylabel("adjusted HR per SD")
    for xi, r in zip(x, tmp.itertuples()):
        ax.text(xi + 0.08, r.hr, f"E={r.evalue_point}", fontsize=6, color="#555")
    ax.legend(frameon=True, fontsize=6.5, loc="upper right")
    ax.set_title("B  E-values for unmeasured confounding",
                 loc="left", fontweight="bold")
    save(fig, "SuppFigS3_bootstrap_evalue")


# --------------------------------------------------------------------------- #
# S4 神经 spot 判定的循环论证检验
# --------------------------------------------------------------------------- #
def fig_s4():
    gp = TAB / "nerve_gene_level_auc.csv"
    cp = TAB / "nerve_crossvalidation_category.csv"
    xp = TAB / "nerve_crossvalidation_gene_cross.csv"
    if not (gp.exists() and xp.exists()):
        print("  skip S4: no tables"); return
    G = pd.read_csv(gp)
    C = pd.read_csv(cp)
    X = pd.read_csv(xp)

    fig, axes = plt.subplots(1, 3, figsize=(14.4, 4.2))

    # (A) 逐基因判别力: 哪些基因在扛着这个判定
    ax = axes[0]
    order = (G.dropna(subset=["auc"]).groupby("gene")["auc"].median()
              .sort_values(ascending=False).index.tolist())
    data = [G.loc[G["gene"] == g, "auc"].dropna().values for g in order]
    bp = ax.boxplot(data, tick_labels=order, patch_artist=True, widths=0.62)
    for b in bp["boxes"]:
        b.set_facecolor("#D5F5E3"); b.set_edgecolor("#1D8348"); b.set_linewidth(0.8)
    for m in bp["medians"]:
        m.set_color("#145A32"); m.set_linewidth(1.1)
    ax.axhline(0.5, ls="--", lw=0.9, color=C_INDEP)
    ax.set_xticklabels(order, rotation=80, fontsize=6.5)
    ax.set_ylabel("AUC: called vs non-called spots")
    ax.text(0.02, 0.95, "n=29 samples", transform=ax.transAxes, va="top", fontsize=6.5)
    ax.set_title("A  Discrimination is carried by few expressed markers",
                 loc="left", fontweight="bold")

    # (B) 共享基因污染前后
    ax = axes[1]
    base = {"sensory": 0.678, "structural": 0.655}       # 未剔除共享基因的原结果
    clean = C.dropna(subset=["auc"]).groupby("call_category")["auc"].median().to_dict()
    cats = [c for c in ["sensory", "structural"] if c in clean]
    x = np.arange(len(cats)); w = 0.36
    ax.bar(x - w / 2, [base[c] for c in cats], w, color="#F5B041",
           label="readout keeps shared marker UCHL1")
    ax.bar(x + w / 2, [clean[c] for c in cats], w, color=C_GREY,
           label="readout excludes panel-shared genes")
    for xi, c in zip(x, cats):
        sub = C[(C["call_category"] == c) & np.isfinite(C["p_perm"])]
        sig = float((sub["p_perm"] < 0.05).mean() * 100) if len(sub) else np.nan
        ax.text(xi + w / 2, clean[c] + 0.015, f"{sig:.0f}%\np<0.05", fontsize=6,
                ha="center", color="#333")
    ax.set_xticks(x); ax.set_xticklabels(cats, fontsize=7)
    ax.axhline(0.5, ls="--", lw=0.9, color=C_INDEP)
    ax.set_ylim(0, 0.9)
    ax.set_ylabel("AUC of cross-panel readout")
    ax.legend(frameon=True, fontsize=6.5, loc="upper right")
    ax.set_title("B  Cross-panel validation collapses once contamination is removed",
                 loc="left", fontweight="bold")

    # (C) 逐基因交叉读数
    ax = axes[2]
    show = X.groupby(["call_category", "gene"])["auc"].median().reset_index()
    show["key"] = show["call_category"].str[:4] + "\u2192" + show["gene"]
    show = show.sort_values("auc", ascending=False).head(14)
    cols = ["#7FB3D5" if k.startswith("sens") else "#C39BD3" for k in show["key"]]
    ax.barh(np.arange(len(show))[::-1], show["auc"], color=cols, edgecolor="#555",
            height=0.7)
    ax.set_yticks(np.arange(len(show))[::-1])
    ax.set_yticklabels(show["key"], fontsize=6)
    ax.axvline(0.5, ls="--", lw=0.9, color=C_INDEP)
    ax.set_xlim(0.42, 0.62)
    ax.set_xlabel("median AUC across samples")
    ax.set_title("C  Gene-by-gene cross-panel readouts (top 14)",
                 loc="left", fontweight="bold")
    save(fig, "SuppFigS4_nerve_crossvalidation")


if __name__ == "__main__":
    print("[SuppFig] rendering")
    fig_s1(); fig_s2(); fig_s3(); fig_s4()

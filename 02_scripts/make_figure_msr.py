# -*- coding: utf-8 -*-
"""Supplementary Figure S5 —— 空间自相关感知的 NPP 重算 (Moran spectral randomization)。

面板
----
A  每样本表达场 Moran's I 分布 vs 距离场 Moran's I (自相关强度的量级对比)
B  每个基因的 MSR 零分布 SD vs 朴素置换零分布 SD (膨胀倍数)
C  四种判据下的显著基因数 (量化假阳性膨胀)
D  原 NPP 基因集与 MSR 基因集的重合度
E  MSR 标准化统计量 u 的观测分布 vs 理论 N(0,1) (整体校准诊断)
F  关键基因的效应量与两种零分布下的显著性对比

风格与 Fig1-4 / SuppFigS1-S4 一致: 8pt, 300 dpi, PNG + PDF。
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
from scipy import stats

sys.path.insert(0, r"D:/嗜神经浸润/00_pipeline")
from nsgrad import config as CFG  # noqa: E402

TAB = CFG.DIR_TAB
FIG = CFG.DIR_FIG
PROC = CFG.DIR_PROC
FIG.mkdir(parents=True, exist_ok=True)
DPI = 300

plt.rcParams.update({
    "font.size": 8, "axes.titlesize": 9, "axes.labelsize": 8,
    "xtick.labelsize": 7, "ytick.labelsize": 7, "legend.fontsize": 7,
    "figure.dpi": DPI, "savefig.dpi": DPI, "savefig.bbox": "tight",
    "font.family": "DejaVu Sans", "axes.linewidth": 0.6,
    "pdf.fonttype": 42, "ps.fonttype": 42,
})
C_MAIN, C_ALT, C_GREY, C_ACC = "#2C5F8A", "#B4534B", "#8A8A8A", "#4E8C6A"


def save(fig, name):
    fig.savefig(FIG / f"{name}.png", dpi=DPI, bbox_inches="tight")
    fig.savefig(FIG / f"{name}.pdf", bbox_inches="tight")
    plt.close(fig)
    print(f"   -> {name}.png/.pdf")


def panel_a(P: pd.DataFrame, ax):
    """每样本 Moran's I 箱线图 + 距离场 I。"""
    per = P.groupby("sample").agg(
        med=("moran_i", "median"), q25=("moran_i", lambda s: np.percentile(s, 25)),
        q75=("moran_i", lambda s: np.percentile(s, 75)),
        lo=("moran_i", lambda s: np.percentile(s, 5)),
        hi=("moran_i", lambda s: np.percentile(s, 95)),
        dI=("moran_i_dist", "first"), n=("n_spots", "first"))
    per = per.sort_values("med")
    x = np.arange(len(per))
    ax.bar(x, per["hi"] - per["lo"], bottom=per["lo"], width=0.7,
           color=C_MAIN, alpha=0.25, label="gene field 5–95%")
    ax.bar(x, per["q75"] - per["q25"], bottom=per["q25"], width=0.7,
           color=C_MAIN, alpha=0.65, label="gene field IQR")
    ax.plot(x, per["med"], "_", color="black", ms=9, label="gene field median")
    ax.plot(x, per["dI"], "o", color=C_ALT, ms=3.5, label="distance field I")
    ax.set_xticks(x[::3])
    ax.set_xticklabels([s.replace("GSM", "") for s in per.index[::3]], rotation=90)
    ax.set_ylabel("Moran's I")
    ax.set_xlabel("Visium sample (ordered by median gene-field I)")
    ax.set_title("A  Spatial autocorrelation: gene field vs nerve-distance field", loc="left")
    ax.legend(frameon=False, loc="upper left", fontsize=5.5, ncol=2,
              columnspacing=0.8, handlelength=1.2)
    ax.axhline(0, color="black", lw=0.4)
    return per


def panel_b(G: pd.DataFrame, ax, ratio):
    ax.scatter(G["naive_sd"], G["msr_sd"], s=3, color=C_MAIN, alpha=0.25,
               edgecolors="none")
    lim = [max(1e-4, float(np.nanpercentile(G["naive_sd"], 0.5))),
           float(np.nanpercentile(G["naive_sd"], 99.5))]
    ax.plot(lim, lim, "--", color="black", lw=0.7, label="y = x (no inflation)")
    for f, c, ls in ((ratio, C_ALT, "-"),):
        ax.plot(lim, [v * f for v in lim], color=c, lw=1.0, ls=ls,
                label=f"median ratio = {f:.2f}×")
    ax.set_xscale("log")
    ax.set_yscale("log")
    from matplotlib.ticker import LogLocator, FuncFormatter
    for axis in (ax.xaxis, ax.yaxis):
        axis.set_major_locator(LogLocator(base=10, numticks=4))
        axis.set_minor_locator(LogLocator(base=10, subs=(0.5,), numticks=4))
        axis.set_major_formatter(FuncFormatter(
            lambda v, _: f"{v:g}" if v >= 0.01 else f"{v:.2g}"))
    ax.set_xlabel("SD of null under naive spot permutation")
    ax.set_ylabel("SD of null under Moran spectral randomisation")
    ax.set_title("B  Null width: MSR vs naive (per gene, n = %d)" % len(G), loc="left")
    ax.legend(frameon=False, loc="upper left", fontsize=6)


def panel_c(n_noeff: dict, sets: dict, ax):
    labels = ["Naive null\n(ignores autocorr.)", "MSR null\n(Stouffer)",
              "MSR null +\nbetween-sample t", "Original\nt-test p<0.05",
              "Original\nBH-FDR<0.05"]
    vals = [n_noeff["naive_fdr05"], n_noeff["msr_fdr05"],
            n_noeff["msr_between_fdr05"], n_noeff["orig_p05"], n_noeff["orig_fdr05"]]
    cols = [C_ALT, C_MAIN, "#6B4E8C", C_GREY, C_GREY]
    ax.bar(np.arange(5), vals, color=cols, alpha=0.88, width=0.62)
    for i, v in enumerate(vals):
        ax.text(i, v, f"{v}", ha="center", va="bottom", fontsize=7)
    ax.set_xticks(range(5))
    ax.set_xticklabels(labels, fontsize=6)
    ax.set_ylabel("Number of genes called significant (of %d)" % n_noeff["n_universe"])
    ax.set_title("C  Genes called significant under five criteria (no effect-size filter)",
                 loc="left")
    ax.set_ylim(0, max(vals) * 1.28 if max(vals) > 0 else 1)
    infl = n_noeff["naive_fdr05"] / max(1, n_noeff["msr_fdr05"])
    txt = (f"Inflation from ignoring autocorrelation = {infl:.1f}×\n"
           f"With |mean rho| ≥ 0.05:  original {len(sets['orig'])}  |  "
           f"naive-null {len(sets['naive'])}  |  MSR-null {len(sets['msr'])}  |  "
           f"MSR+between-sample t {len(sets['msr_between'])}")
    ax.text(0.02, 0.94, txt, transform=ax.transAxes, fontsize=6, color="#333333",
            va="top")


def panel_d(sets: dict, overlap: list, jac: float, hyp_p, n_universe, ax):
    o, m, i = len(sets["orig"]), len(sets["msr"]), len(overlap)
    venn_like = {"Original NPP only": o - i, "Both": i, "MSR-null only": m - i}
    ax.barh([2], [o - i], color=C_GREY, alpha=0.85, height=0.5)
    ax.barh([1], [i], color=C_MAIN, height=0.5)
    ax.barh([0], [m - i], color=C_ALT, alpha=0.85, height=0.5)
    for y, k in zip([2, 1, 0], ["Original NPP only", "Both", "MSR-null only"]):
        ax.text(venn_like[k] + max(1, 0.02 * max(o, m)), y, str(venn_like[k]),
                va="center", fontsize=7)
    ax.set_yticks([0, 1, 2])
    ax.set_yticklabels(["MSR-null only", "Both", "Original NPP only"], fontsize=7)
    ax.set_xlabel("Number of genes")
    ax.set_title("D  Overlap of original NPP and MSR-validated gene sets", loc="left")
    surv = i / max(1, o)
    txt = (f"Jaccard = {jac:.3f}\n"
           f"Survival of original set = {surv:.1%}\n"
           f"Hypergeometric P = {hyp_p:.2e}\n"
           f"Under the most conservative criterion\n"
           f"(MSR null + between-sample t): 2 genes survive"
           if hyp_p is not None else
           f"Jaccard = {jac:.3f}\nSurvival of original set = {surv:.1%}")
    ax.text(0.98, 0.05, txt, transform=ax.transAxes, ha="right", va="bottom",
            fontsize=6.5,
            bbox=dict(boxstyle="round,pad=0.3", fc="#F5F5F5", ec="#CCCCCC", lw=0.5))
    ax.set_xlim(0, max(o, m, 1) * 1.25)


def panel_e(P: pd.DataFrame, ax):
    u = P["u_msr"].dropna().to_numpy()
    u = u[np.isfinite(u)]
    u = u[(u > -8) & (u < 8)]
    ax.hist(u, bins=80, density=True, color=C_MAIN, alpha=0.65,
            label=f"observed (n = {len(u):,})")
    xs = np.linspace(-4, 4, 400)
    ax.plot(xs, stats.norm.pdf(xs), color="black", lw=1.0, label="N(0, 1) reference")
    ax.axvline(0, color="black", lw=0.4)
    ax.set_xlabel("MSR-standardised statistic  u  (per gene per sample)")
    ax.set_ylabel("Density")
    ax.set_title("E  Calibration of the MSR null across all gene × sample pairs",
                 loc="left")
    frac = float(np.mean(np.abs(u) > 1.96))
    sd = float(np.std(u))
    ax.text(0.98, 0.95, f"SD = {sd:.2f}\n|u| > 1.96 in {frac:.1%} (expect 5.0%)",
            transform=ax.transAxes, ha="right", va="top", fontsize=6.5,
            bbox=dict(boxstyle="round,pad=0.3", fc="#F5F5F5", ec="#CCCCCC", lw=0.5))
    ax.legend(frameon=False, loc="upper left", fontsize=6)


def panel_f(G: pd.DataFrame, genes: list, ax):
    rows = [g for g in genes if g in G.index]
    if not rows:
        ax.axis("off")
        return
    sub = G.loc[rows]
    y = np.arange(len(sub))
    ax.barh(y + 0.2, -np.log10(sub["fdr_msr"].clip(lower=1e-8)), height=0.38,
            color=C_MAIN, alpha=0.9, label="−log10 FDR, MSR null (Stouffer)")
    ax.barh(y - 0.2, -np.log10(sub["fdr_msr_between"].clip(lower=1e-8)), height=0.38,
            color="#6B4E8C", alpha=0.85, label="−log10 FDR, MSR + between-sample t")
    for yy, (g, r) in zip(y, sub.iterrows()):
        ax.text(0.05, yy + 0.2, f"rho {r['mean_rho']:+.3f}", fontsize=5.5,
                va="center", color="#333333")
    ax.set_yticks(y)
    ax.set_yticklabels(sub.index, fontsize=7)
    ax.axvline(0, color="black", lw=0.4)
    ax.axvline(-np.log10(0.05), color=C_ALT, lw=0.7, ls="--")
    ax.text(-np.log10(0.05) + 0.15, len(sub) - 0.4, "FDR 0.05", fontsize=6,
            color=C_ALT, ha="left")
    ax.set_xlabel("−log10(FDR)")
    ax.set_title("F  Sentinel genes under the MSR null and the conservative null",
                 loc="left")
    ax.legend(frameon=False, loc="lower right", fontsize=6)


def main():
    P = pd.read_csv(PROC / "msr_per_sample.fine.csv.gz")
    G = pd.read_csv(PROC / "msr_gene_level.fine.csv", index_col=0)
    S = json.load(open(TAB / "msr_summary.fine.json", encoding="utf-8"))
    print(f"[S5] per-sample {P.shape}, gene-level {G.shape}")

    fig, axes = plt.subplots(2, 3, figsize=(13.2, 7.6))
    panel_a(P, axes[0, 0])
    panel_b(G, axes[0, 1], S["sd_ratio_msr_over_naive"])
    ne = dict(S["n_significant_no_effect_threshold"])
    ne["n_universe"] = len(G)
    panel_c(ne, {"orig": S["set_orig_relaxed"], "naive": S["set_naive"],
                 "msr": S["set_msr"], "msr_between": S["set_msr_between"]},
            axes[0, 2])
    panel_d({"orig": S["set_orig_relaxed"], "msr": S["set_msr"]},
            S["overlap"], S["jaccard"], S["hypergeom_p"], len(G), axes[1, 0])
    panel_e(P, axes[1, 1])
    sent = ["UCHL1", "TUBA1A", "TUBB3", "MMP13", "IFI27", "TPM1", "CCN2",
            "COL4A1", "VIM", "MMP9", "HLA-B", "SPARC", "POSTN", "LRRC15"]
    panel_f(G, sent, axes[1, 2])
    fig.tight_layout()
    save(fig, "SuppFigS5_spatial_autocorrelation")


if __name__ == "__main__":
    main()

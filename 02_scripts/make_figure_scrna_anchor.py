# -*- coding: utf-8 -*-
"""
方案 2-A 图：scRNA 锚定 + 解卷积 + 共聚集。
输出 03_results/figures/SuppFigS9_scrna_anchor.png/.pdf（多面板组合）。

面板:
  A  模块基因 × 细胞类型 检出率热图（上: 神经固有 39, 下: 基质 60）
  B  各模块基因来源归属的构成（top-type 计数条形）
  C  解卷积: 各细胞类型丰度随距神经环带的变化（样本均值 ± 95% CI）
  D  丰度~距离 偏 Spearman 森林图（Fisher-z meta, 主分析 + no-module 敏感性）
  E  模块得分~细胞丰度 共聚集点图（rho ± 95% CI）
  F  NNLS QC: 每 spot 重构相对误差分布（按样本）
"""
from __future__ import annotations
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
PROC = ROOT / "01_data/processed"
FIG = ROOT / "03_results/figures"

NEURAL_INTRINSIC = [
    "UCHL1", "TUBB3", "TUBA1A", "MAP2", "NEFL", "NEFM", "PLP1", "MPZ", "PMP22",
    "S100B", "NGFR", "SOX10", "SYP", "SCN9A", "SCN10A", "SCN7A", "PRX", "MBP",
    "TH", "DBH", "CHAT", "SLC18A3", "SLC5A7", "TRPV1", "CALCA", "CALCB", "VIP",
    "TAC1", "NPY", "SST", "PENK", "GAD1", "GAD2", "SLC6A1", "GAP43", "ELAVL4",
    "STMN2", "NRXN1", "SYT1",
]
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
RING_LABELS = {"r1_0_50": "0–50", "r2_50_100": "50–100", "r3_100_200": "100–200",
               "r4_200_400": "200–400", "r5_gt400": ">400"}
RING_ORDER = list(RING_LABELS)
TYPE_COLORS = {
    "CAF_myCAF": "#C0392B", "CAF_iCAF": "#E67E22", "CAF_transitioning": "#F5B041",
    "Cancer Epithelial": "#7D3C98", "Normal Epithelial": "#BB8FCE",
    "Endothelial": "#2471A3", "PVL": "#16A085", "Myeloid": "#117864",
    "T-cells": "#2E86C1", "B-cells": "#A9CCE3", "Plasmablasts": "#5D6D7E",
}

plt.rcParams.update({
    "font.size": 7, "axes.linewidth": 0.6, "font.family": "DejaVu Sans",
    "pdf.fonttype": 42, "figure.dpi": 300,
})


def panel_a(ax):
    long = pd.read_csv(PROC / "scrna_cell_origin_long.csv")
    piv = long.pivot_table(index="gene", columns="celltype", values="frac_expr")
    mod = pd.read_csv(PROC / "scrna_module_origin.csv")
    types = [t for t in piv.columns]
    genes = [g for g in NEURAL_INTRINSIC + PERINEURAL_STROMA]
    M = np.full((len(genes), len(types)), np.nan)
    for i, g in enumerate(genes):
        if g in piv.index:
            M[i, :] = piv.loc[g, types].values
    im = ax.imshow(M, aspect="auto", cmap="magma_r", vmin=0, vmax=1,
                   interpolation="nearest")
    ax.set_xticks(range(len(types)))
    ax.set_xticklabels([t.replace(" ", "\n").replace("CAF_", "CAF\n") for t in types],
                       fontsize=5.2)
    yticks = [4, 20, 38, 68, 98]
    ax.set_yticks(yticks)
    ax.set_yticklabels(yticks, fontsize=5)
    ax.set_ylabel("module genes\n(top: neural-intrinsic; bottom: stromal)", fontsize=6)
    # 模块分隔线
    ax.axhline(39 - 0.5, color="white", lw=1.2)
    ax.text(len(types) + 0.3, 19, "neural\n(39)", fontsize=5.5, va="center", color="#555")
    ax.text(len(types) + 0.3, 70, "stromal\n(60)", fontsize=5.5, va="center", color="#555")
    cb = plt.colorbar(im, ax=ax, fraction=0.03, pad=0.14)
    cb.set_label("fraction of cells expressing", fontsize=5.5)
    cb.ax.tick_params(labelsize=5)
    return ax


def panel_b(ax):
    mod = pd.read_csv(PROC / "scrna_module_origin.csv")
    for k, (mname, color) in enumerate(
            [("Neural_intrinsic", "#2E4053"), ("Perineural_stroma", "#B03A2E")]):
        sub = mod[mod.module == mname]
        ok = sub[sub.status == "measured"]
        top = ok["top_frac_type"].value_counts()
        n = len(ok)
        ypos = k
        left = 0
        order = ["CAF_myCAF", "CAF_iCAF", "PVL", "Endothelial", "Myeloid",
                 "Cancer Epithelial", "Normal Epithelial", "T-cells", "B-cells",
                 "Plasmablasts"]
        for t in order:
            if t in top.index:
                v = top[t] / n
                ax.barh(ypos, v, left=left, height=0.55,
                        color=TYPE_COLORS.get(t, "#999"), edgecolor="white", lw=0.4)
                if v >= 0.09:
                    ax.text(left + v / 2, ypos, f"{int(top[t])}", ha="center",
                            va="center", fontsize=5.5, color="white")
                left += v
        na = len(sub[sub.status != "measured"])
        if na:
            ax.barh(ypos, na / len(sub), left=left, height=0.55, color="#D5DBDB",
                    edgecolor="white", lw=0.4, hatch="//")
        ax.text(1.01, ypos, f"n={n} measured, {na} not measured", fontsize=5.5,
                va="center")
    ax.set_yticks([0, 1])
    ax.set_yticklabels(["Neural-intrinsic", "Perineural-stroma"], fontsize=6.5)
    ax.set_xlim(0, 1); ax.set_xlabel("fraction of module genes", fontsize=6.5)
    ax.spines[["top", "right"]].set_visible(False)


def panel_c(ax):
    ring = pd.read_csv(PROC / "deconv_ring_profiles.csv")
    ring = ring[ring.ring_id.isin(RING_ORDER)]
    types = ["CAF_myCAF", "Cancer Epithelial", "T-cells", "Myeloid", "Endothelial", "PVL"]
    for t in types:
        col = f"frac_{t}"
        g = ring.groupby("ring_id")[col]
        m = g.mean().reindex(RING_ORDER)
        se = g.sem().reindex(RING_ORDER)
        c = TYPE_COLORS.get(t, "#333")
        ax.plot(range(5), m.values, "-o", ms=2.5, lw=1.0, color=c,
                label=t.replace("Cancer Epithelial", "Cancer epithelial"))
        ax.fill_between(range(5), m - 1.96 * se, m + 1.96 * se, color=c, alpha=0.15, lw=0)
    ax.set_xticks(range(5)); ax.set_xticklabels([RING_LABELS[r] for r in RING_ORDER])
    ax.set_xlabel("distance to nearest nerve (µm)", fontsize=6.5)
    ax.set_ylabel("deconvolved abundance\n(mean fraction per spot)", fontsize=6.5)
    ax.legend(fontsize=5.2, frameon=False, ncol=2, loc="upper center")
    ax.spines[["top", "right"]].set_visible(False)


def panel_d(ax):
    meta = pd.read_csv(PROC / "deconv_gradient_meta.csv")
    sens = pd.read_csv(PROC / "deconv_nerve_vs_non_meta.csv")  # noqa: F841
    try:
        nom = pd.read_csv(PROC / "deconv_gradient_meta_nomodule.csv")
    except FileNotFoundError:
        nom = None
    n = len(meta)
    ys = np.arange(n)[::-1]
    for y, (_, r) in zip(ys, meta.iterrows()):
        c = TYPE_COLORS.get(r.celltype, "#333")
        ax.plot([r.rho_mean, r.rho_mean], [y - 0.18, y + 0.18], color=c, lw=1.4)
        if nom is not None:
            rr = nom[nom.celltype == r.celltype]
            if len(rr):
                ax.plot(rr.rho_mean.iloc[0], y - 0.38, marker="x", ms=4, color=c,
                        mew=0.9)
        if r.fdr_bh < 0.05:
            ax.text(-0.5, y, "*", fontsize=7, ha="right", va="center")
    ax.axvline(0, color="#888", ls=":", lw=0.6)
    ax.set_yticks(ys); ax.set_yticklabels(meta.celltype, fontsize=6)
    ax.set_xlabel("partial Spearman ρ (abundance vs distance; negative = perineural)",
                  fontsize=6.5)
    ax.spines[["top", "right"]].set_visible(False)
    ax.set_ylim(-0.8, n - 0.2)


def panel_e(ax):
    cm = pd.read_csv(PROC / "deconv_module_coaggr_meta.csv")
    ys = np.arange(len(cm))[::-1]
    for y, (_, r) in zip(ys, cm.iterrows()):
        c = TYPE_COLORS.get(r.celltype, "#333")
        # 神经模块
        ax.plot(r.rho_neural, y + 0.15, "o", ms=4, color="#2E4053",
                markeredgecolor="white", markeredgewidth=0.4)
        if r.fdr_neural < 0.05:
            ax.plot(r.rho_neural + 0.015, y + 0.15, "*", ms=4, color="#2E4053")
        # 基质模块
        ax.plot(r.rho_stroma, y - 0.15, "o", ms=4, color="#B03A2E",
                markeredgecolor="white", markeredgewidth=0.4)
        if r.fdr_stroma < 0.05:
            ax.plot(r.rho_stroma + 0.015, y - 0.15, "*", ms=4, color="#B03A2E")
    ax.axvline(0, color="#888", ls=":", lw=0.6)
    ax.set_yticks(ys); ax.set_yticklabels(cm.celltype, fontsize=6)
    ax.set_xlabel("partial Spearman ρ with spot-level module score", fontsize=6.5)
    ax.plot([], [], "o", ms=4, color="#2E4053", label="neural-intrinsic")
    ax.plot([], [], "o", ms=4, color="#B03A2E", label="perineural-stroma")
    ax.plot([], [], "*", ms=4, color="#333", label="FDR < 0.05")
    ax.legend(fontsize=5.2, frameon=False, loc="lower right")
    ax.spines[["top", "right"]].set_visible(False)


def panel_f(ax):
    qc = pd.read_csv(PROC / "deconv_qc.csv")
    x = np.arange(len(qc)) + 1
    ax.vlines(x, 0, qc.median_rel_err, color="#2471A3", lw=0.8)
    ax.plot(x, qc.median_rel_err, "o", ms=2.2, color="#2471A3")
    med = qc.median_rel_err.median()
    ax.axhline(med, color="#C0392B", ls="--", lw=0.8)
    ax.text(0.98, med + 0.015, f"median {med:.2f}", fontsize=5.5, color="#C0392B",
            ha="right", transform=ax.get_yaxis_transform())
    ax.set_xlabel("Visium sample (n=30)", fontsize=6.5)
    ax.set_ylabel("NNLS relative reconstruction error", fontsize=6.5)
    ax.set_ylim(0, None)
    ax.spines[["top", "right"]].set_visible(False)


def main():
    fig = plt.figure(figsize=(7.09, 10.6))  # 180 mm 宽
    gs = fig.add_gridspec(4, 2, height_ratios=[2.1, 1, 1, 0.8], hspace=0.5, wspace=0.32,
                          left=0.13, right=0.97, top=0.97, bottom=0.04)
    panel_a(fig.add_subplot(gs[0, :]))
    panel_b(fig.add_subplot(gs[1, 0]))
    panel_c(fig.add_subplot(gs[1, 1]))
    panel_d(fig.add_subplot(gs[2, 0]))
    panel_e(fig.add_subplot(gs[2, 1]))
    panel_f(fig.add_subplot(gs[3, :]))
    fig.savefig(FIG / "SuppFigS9_scrna_anchor.png", dpi=350, bbox_inches="tight")
    fig.savefig(FIG / "SuppFigS9_scrna_anchor.pdf", bbox_inches="tight")
    print("saved SuppFigS9_scrna_anchor")


if __name__ == "__main__":
    main()

# -*- coding: utf-8 -*-
"""Figure 6 —— 神经周围微环境的功能读数梯度(配体受体 / 通路 / 细胞类型)。

所有读数用 Moran 谱随机化 (MSR) 零分布校准, 再与**大小匹配的随机基因集**
基线比较(第三道门槛), 只报告超出基线的结果。

面板
----
A  细胞类型丰度随距离环的剖面 (z 空间, 跨样本均值)
B  校准散点: |mean rho| vs 基因集大小; 随机集基线 + 真实读数
C  细胞类型: 跨样本 mean rho 与逐样本分布
D  通路: 通过三道门槛者 + 最强未通过者
E  配体-受体对: 按通讯轴着色
F  随机集零分布 vs 观测 |Stouffer Z| (校准诊断)

风格与 Fig1-5 一致: 8pt, 300 dpi, PNG + PDF。
"""
from __future__ import annotations

import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, r"D:/嗜神经浸润/00_pipeline")
from nsgrad import config as CFG  # noqa: E402

TAB = CFG.DIR_TAB
FIG = CFG.DIR_FIG
PROC = CFG.DIR_PROC
FIG.mkdir(parents=True, exist_ok=True)
DPI = 300

plt.rcParams.update({
    "font.size": 8, "axes.titlesize": 9, "axes.labelsize": 8,
    "xtick.labelsize": 7, "ytick.labelsize": 7, "legend.fontsize": 6.5,
    "figure.dpi": DPI, "savefig.dpi": DPI, "savefig.bbox": "tight",
    "font.family": "DejaVu Sans", "axes.linewidth": 0.6,
    "pdf.fonttype": 42, "ps.fonttype": 42,
})
C_MAIN, C_ALT, C_GREY, C_ACC = "#2C5F8A", "#B4534B", "#8A8A8A", "#4E8C6A"
C_PURPLE = "#6B4E8C"


def save(fig, name):
    fig.savefig(FIG / f"{name}.png", dpi=DPI, bbox_inches="tight")
    fig.savefig(FIG / f"{name}.pdf", bbox_inches="tight")
    plt.close(fig)
    print(f"   -> {name}.png/.pdf")


RING_LABEL = {"r1_0_50": "0-50", "r2_50_100": "50-100", "r3_100_200": "100-200",
              "r4_200_400": "200-400", "out_gt400": ">400"}
RING_ORDER = list(RING_LABEL)


def panel_a(M, ax):
    keep = ["Neuronal_Schwann", "Pericyte_SMC", "CAF_myofibroblastic",
            "CAF_inflammatory_dev", "Endothelial", "Macrophage_monocyte",
            "Tumor_proliferative", "CD8_T", "T_lymphocyte"]
    keep = [k for k in keep if k in M.columns]
    cols = [C_PURPLE, C_MAIN, C_ALT, "#D98C3F", C_ACC, "#7A6E5B",
            "#5B7A8C", "#3B7A5B", "#9A3B6B"]
    xs = np.arange(len(RING_ORDER))
    for c, k in zip(cols, keep):
        v = [float(M.loc[r, k]) if r in M.index else np.nan for r in RING_ORDER]
        ax.plot(xs, v, "-o", ms=3.2, lw=1.1, color=c, label=k.replace("_", " "))
    ax.axhline(0, color="black", lw=0.5)
    ax.set_xticks(xs)
    ax.set_xticklabels([RING_LABEL[r] for r in RING_ORDER])
    ax.set_xlabel("Distance to nearest nerve (µm)")
    ax.set_ylabel("Mean score (z, across 30 samples)")
    ax.set_title("A  Cell-type abundance by perineural distance ring", loc="left")
    ax.legend(frameon=False, fontsize=5.6, ncol=2, columnspacing=0.7,
              handlelength=1.1, loc="upper right")


def panel_b(R, C, ax):
    ax.scatter(R["size"], R["abs_rho"], s=5, color=C_GREY, alpha=0.35,
               edgecolors="none", label="size-matched random gene sets (n=510)")
    s = np.sort(R["size"].unique())
    mu = [np.mean(R.loc[R["size"] == v, "abs_rho"]) for v in s]
    ax.plot(s, mu, "-", color="black", lw=1.2, label="random-set baseline")
    sd = [np.std(R.loc[R["size"] == v, "abs_rho"], ddof=1) for v in s]
    ax.fill_between(s, np.asarray(mu) + np.asarray(sd),
                    np.maximum(0, np.asarray(mu) - np.asarray(sd)),
                    color="black", alpha=0.10, label="±1 SD")
    kind_col = {"celltype": C_ACC, "lr": C_ALT, "pathway": C_MAIN}
    for kind in ["pathway", "lr", "celltype"]:
        g = C[(C["kind"] == kind) & C["measurable"] & (~C["circular"])]
        ax.scatter(g["size"], g["mean_rho"].abs(), s=9,
                   color=kind_col[kind], alpha=0.75, edgecolors="none",
                   label=f"{kind} readouts (n={len(g)})")
    ax.set_xscale("log")
    ax.set_xlabel("Gene-set size (number of member genes)")
    ax.set_ylabel("|mean partial ρ with nerve distance|")
    ax.set_title("B  Effect size vs size-matched random baseline", loc="left")
    ax.legend(frameon=False, fontsize=5.8, loc="upper left")


def panel_c(C, P, ax):
    g = C[(C["kind"] == "celltype") & C["measurable"]].sort_values("mean_rho")
    y = np.arange(len(g))
    ax.barh(y, g["mean_rho"], color=C_MAIN, alpha=0.85, height=0.62)
    for i, (nm, r) in enumerate(g.iterrows()):
        v = P.loc[P["name"] == nm, "rho"].dropna().to_numpy(float)
        if len(v):
            ax.plot(v, np.full(len(v), i) + np.random.default_rng(1).uniform(
                -0.22, 0.22, len(v)), "|", color="black", ms=4, alpha=0.65)
        ax.text(r["mean_rho"] - 0.004 if r["mean_rho"] < 0 else r["mean_rho"] + 0.004,
                i, f"n={int(r['n_samples'])}", va="center",
                ha="right" if r["mean_rho"] < 0 else "left", fontsize=5.5,
                color="#555555")
    ax.axvline(0, color="black", lw=0.5)
    ax.set_yticks(y)
    ax.set_yticklabels([i.replace("_", " ") for i in g.index], fontsize=7)
    ax.set_xlabel("Partial ρ with distance to nearest nerve (negative = perineural)")
    ax.set_title("C  Cell-type gradients: cross-sample mean + per-sample values",
                 loc="left")


def panel_d(C, ax):
    g = C[(C["kind"] == "pathway") & C["measurable"] & (~C["circular"])]
    top = g.sort_values("z_vs_random", ascending=False).head(16)[::-1]
    y = np.arange(len(top))
    passed = top["pass_calibrated"].to_numpy(bool)
    ax.barh(y, top["mean_rho"].abs(), color=[C_MAIN if p else "#C6CBD1"
                                             for p in passed], height=0.66)
    base = top["baseline_rho"].to_numpy(float)
    ax.plot(base, y, "_", color=C_ALT, ms=8, label="size-matched random baseline")
    labels = []
    for nm, r in top.iterrows():
        lab = nm.replace("REACTOME_", "").replace("HALLMARK_", "") \
                .replace("KEGG_MEDICUS_REFERENCE_", "") \
                .replace("KEGG_", "").replace("_", " ")
        labels.append((lab[:52] + "…") if len(lab) > 53 else lab)
    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=5.6)
    for tl, p in zip(ax.get_yticklabels(), passed):
        tl.set_color("black" if p else "#666666")
    ax.set_xlabel("|mean partial ρ| (grey = does not exceed random baseline)")
    ax.set_title("D  Pathways ranked by excess over random baseline\n"
                 "(* = passes all three thresholds)", loc="left", fontsize=8)
    ax.legend(frameon=False, fontsize=6, loc="lower right")


def panel_e(C, ax):
    g = C[(C["kind"] == "lr") & C["measurable"] & (~C["circular"])]
    g = g.sort_values("mean_rho").head(14)[::-1]
    y = np.arange(len(g))
    ax.barh(y, g["mean_rho"], color=[C_ALT if p else "#C6CBD1"
                                     for p in g["pass_calibrated"]], height=0.66)
    labels = []
    for nm, r in g.iterrows():
        parts = nm.split("|")
        pair = parts[1] if len(parts) > 1 else nm
        axis = parts[2] if len(parts) > 2 else ""
        labels.append(f"{pair}  ({axis})")
    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=6)
    ax.axvline(0, color="black", lw=0.5)
    ax.set_xlabel("Partial ρ (negative = higher near nerve)")
    ax.set_title("E  Ligand-receptor co-expression gradients", loc="left")
    ax.text(0.02, 0.04, "None exceeds the size-matched random baseline on "
            "cross-sample consistency", transform=ax.transAxes, fontsize=5.8,
            color=C_ALT)


def panel_f(R, C, ax):
    rz = np.sort(R["abs_rho"].to_numpy(float))
    obs = C.loc[C["measurable"] & (~C["circular"]), "mean_rho"].abs().to_numpy(float)
    obs = np.sort(obs)
    ecdf_r = np.arange(1, len(rz) + 1) / len(rz)
    ecdf_o = np.arange(1, len(obs) + 1) / len(obs)
    ax.plot(rz, ecdf_r, color=C_GREY, lw=1.2, label=f"random sets (n={len(rz)})")
    ax.plot(obs, ecdf_o, color=C_MAIN, lw=1.2,
            label=f"real readouts (n={len(obs)})")
    ax.set_xlabel("|mean partial ρ|")
    ax.set_ylabel("ECDF")
    ax.set_title("F  Distribution of effect sizes: real vs random", loc="left")
    ax.legend(frameon=False, fontsize=6, loc="upper left")
    med_r, med_o = float(np.median(rz)), float(np.median(obs))
    ax.text(0.98, 0.30, f"median |ρ|: random {med_r:.3f} vs real {med_o:.3f}",
            transform=ax.transAxes, ha="right", fontsize=6,
            bbox=dict(boxstyle="round,pad=0.3", fc="#F5F5F5", ec="#CCCCCC", lw=0.5))


def main():
    R = pd.read_csv(TAB / "niche_random_null.csv")
    C = pd.read_csv(TAB / "niche_calibrated.csv", index_col=0)
    P = pd.read_csv(PROC / "niche_per_sample.csv.gz")
    M = pd.read_csv(TAB / "niche_ring_meta.csv", index_col=0)
    print(f"[Fig6] random {R.shape}, calibrated {C.shape}, per-sample {P.shape}")

    fig, axes = plt.subplots(2, 3, figsize=(13.6, 8.2))
    panel_a(M, axes[0, 0])
    panel_b(R, C, axes[0, 1])
    panel_c(C, P, axes[0, 2])
    panel_d(C, axes[1, 0])
    panel_e(C, axes[1, 1])
    panel_f(R, C, axes[1, 2])
    fig.tight_layout()
    save(fig, "Fig6_perineural_niche")


if __name__ == "__main__":
    main()

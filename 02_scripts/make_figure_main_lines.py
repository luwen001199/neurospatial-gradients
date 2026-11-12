# -*- coding: utf-8 -*-
"""Figure 4 —— 三条主线: perineural stroma -> myCAF 空间梯度 / 免疫几何 / CAF 四轴分层。

面板
----
A  基质程序与 myCAF/ECM 标志基因的距离梯度森林图 (偏 Spearman meta,
   MSR + 大小匹配随机基线双重校准)
B  环带衰减曲线: Perineural-stroma 模块 / myCAF 轴 / LRRC15 / CD8 / TLS
C  免疫排斥几何: IEI 分布 (CD8 vs TLS) —— CD8 无排斥, TLS 近端富集
D  myCAF 与 CD8 的样本级正相关 (近端/远端分层)
E  CAF 四轴的环带分层 (NNLS 权重)
F  神经根伪时间轨迹: 四轴 + CD8/T 沿 dpt 十分位

风格与 Fig1-5 一致: 8pt, 300 dpi, PNG + PDF。
"""
from __future__ import annotations

import sys

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
C_ORANGE = "#C07A2B"

RING_LABEL = {"r1_0_50": "0-50", "r2_50_100": "50-100", "r3_100_200": "100-200",
              "r4_200_400": "200-400", "out_gt400": ">400"}
RING_MID = {"r1_0_50": 25, "r2_50_100": 75, "r3_100_200": 150, "r4_200_400": 300,
            "out_gt400": 600}
NICE = {"MOD|Perineural_stroma": "Perineural-stroma module",
        "MOD|Collagen_core": "Collagen core",
        "AXIS|myCAF": "myCAF axis", "AXIS|iCAF": "iCAF axis",
        "AXIS|apCAF": "apCAF axis", "AXIS|LRRC15p": "LRRC15$^{+}$",
        "CT|Macrophage_monocyte": "Macrophage",
        "GENE|SPARC": "SPARC", "GENE|POSTN": "POSTN", "GENE|ACTA2": "ACTA2",
        "GENE|FN1": "FN1", "GENE|LRRC15": "LRRC15", "GENE|TNC": "TNC",
        "GENE|COL1A1": "COL1A1", "GENE|LOXL2": "LOXL2"}


def save(fig, name):
    fig.savefig(FIG / f"{name}.png", dpi=DPI, bbox_inches="tight")
    fig.savefig(FIG / f"{name}.pdf", dpi=DPI, bbox_inches="tight")
    plt.close(fig)
    print(f"   -> {name}.png/.pdf")


# ------------------------------------------------------------------ 面板
def panel_a(ax):
    M = pd.read_csv(TAB / "axis_meta_summary.csv")
    M = M[M.kind.isin(["MOD", "AXIS", "GENE", "CT"])].copy()
    M["lab"] = M["name"].map(NICE)
    M = M.sort_values("rho_mean")
    y = np.arange(len(M))
    ax.errorbar(M["rho_mean"], y, xerr=[M["rho_mean"] - M["rho_lo"],
                                        M["rho_hi"] - M["rho_mean"]],
                fmt="o", ms=3.2, lw=0.7, color=C_MAIN, ecolor=C_MAIN,
                elinewidth=0.7, capsize=0)
    # 随机基线
    ax.scatter(M["baseline_rho"] * -1 * np.sign(M["rho_mean"]), y, marker="|",
               s=28, color=C_ALT, lw=1.0, zorder=3, label="size-matched random baseline")
    ax.axvline(0, color=C_GREY, lw=0.6, ls=":")
    for i, r in enumerate(M.itertuples()):
        if r.size == 1 and abs(r.rho_mean) < 0.04:
            continue
    ax.set_yticks(y)
    ax.set_yticklabels(M["lab"])
    ax.set_xlabel("partial Spearman $\\rho$ (readout vs distance)")
    ax.set_title("A  Perineural-stroma program and myCAF markers decay\n"
                 "with distance (MSR + size-matched calibration)", loc="left", fontsize=9)
    ax.legend(frameon=False, loc="lower right", fontsize=6.5)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)


def panel_b(ax):
    R = pd.read_csv(PROC / "axis_ring_profiles.csv")
    picks = ["MOD|Perineural_stroma", "AXIS|myCAF", "GENE|LRRC15", "CT|Macrophage_monocyte"]
    colors = [C_MAIN, C_ACC, C_ALT, C_GREY]
    labels = ["module (60)", "myCAF axis (20)", "LRRC15", "Macrophage"]
    for p, c, lb in zip(picks, colors, labels):
        g = R.groupby("ring")[p]
        sub = g.mean().reindex(RING_LABEL)
        se = g.sem().reindex(RING_LABEL)
        ks = [k for k in RING_LABEL if pd.notna(sub[k])]
        ax.errorbar([RING_MID[k] for k in ks], [sub[k] for k in ks],
                    yerr=[1.96 * se[k] for k in ks],
                    fmt="-o", ms=2.8, lw=1.0, color=c, label=lb, capsize=0)
    ax.axhline(0, color=C_GREY, lw=0.6, ls=":")
    ax.set_xscale("log")
    ax.set_xticks([25, 75, 150, 300, 600])
    ax.set_xticklabels(["25", "75", "150", "300", ">400"])
    ax.minorticks_off()
    ax.set_xlabel("distance to nearest nerve ($\\mu$m, ring midpoint)")
    ax.set_ylabel("ring mean $z$-score")
    ax.set_title("B  Ring-binned decay curves", loc="left", fontsize=9)
    ax.legend(frameon=False, fontsize=6.5)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)


def panel_c(ax):
    P = pd.read_csv(PROC / "immune_exclusion_per_sample.csv")
    for col, c, lb, dx in [("iei_CD8", C_MAIN, "CD8 density", -0.12),
                           ("iei_TLS", C_ALT, "TLS score", +0.12)]:
        v = P[col].dropna()
        rng = np.random.default_rng(3)
        x = rng.uniform(dx - 0.06, dx + 0.06, len(v))
        ax.scatter(x, v, s=9, color=c, alpha=0.75, edgecolors="white", lw=0.3)
        ax.hlines(v.median(), dx - 0.14, dx + 0.14, color=c, lw=1.4)
        m = stats.wilcoxon(v).pvalue
        ax.text(dx, 3.1, f"median {v.median():+.2f}\n$P$={m:.3f}", ha="center",
                fontsize=6.5, color=c)
    ax.axhline(0, color=C_GREY, lw=0.7, ls="--")
    ax.set_xticks([-0.12, 0.12])
    ax.set_xticklabels(["CD8 density", "TLS score"])
    ax.set_xlim(-0.35, 0.35)
    ax.set_ylim(-3.4, 4.0)
    ax.set_ylabel("immune exclusion index\nlog$_2$ (proximal / distal)")
    ax.set_title("C  No CD8 exclusion; TLS enriched near nerves", loc="left", fontsize=9)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)


def panel_d(ax):
    P = pd.read_csv(PROC / "immune_exclusion_per_sample.csv")
    v1 = P["rho_mycaf_cd8_prox"].dropna()
    v2 = P["rho_mycaf_cd8_dist"].dropna()
    for v, c, lb, dx in [(v1, C_ACC, "proximal ($\\leq$100 $\\mu$m)", -0.09),
                         (v2, C_MAIN, "distal (>200 $\\mu$m)", +0.09)]:
        rng = np.random.default_rng(5)
        x = rng.uniform(dx - 0.04, dx + 0.04, len(v))
        ax.scatter(x, v, s=9, color=c, alpha=0.75, edgecolors="white", lw=0.3)
        ax.hlines(v.median(), dx - 0.10, dx + 0.10, color=c, lw=1.4)
        t, p = stats.wilcoxon(v)
        ptxt = "$P$<0.0001" if p < 1e-4 else f"$P$={p:.4f}"
        ax.text(dx, v.min() - 0.012, f"median {v.median():+.3f}\n{ptxt}",
                ha="center", fontsize=6.5, color=c)
    ax.axhline(0, color=C_GREY, lw=0.7, ls="--")
    ax.set_xticks([-0.09, 0.09])
    ax.set_xticklabels(["proximal", "distal"])
    ax.set_xlim(-0.25, 0.25)
    ax.set_ylabel("partial $\\rho$ (myCAF vs CD8 | UMI, MT)")
    ax.set_title("D  myCAF co-localises with CD8, not exclusion", loc="left", fontsize=9)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)


def panel_e(ax):
    B = pd.read_csv(PROC / "caf_axis_traj_bins.csv")  # 未用; 环带用 per_spot
    PS = pd.read_csv(PROC / "caf_axis_per_spot.csv")
    axes = ["myCAF", "iCAF", "apCAF", "LRRC15p"]
    colors = [C_ACC, C_MAIN, C_PURPLE, C_ALT]
    for k, c in zip(axes, colors):
        sub = PS.groupby("ring")[f"nnls_{k}"].mean().reindex(RING_LABEL)
        ks = [r for r in RING_LABEL if pd.notna(sub[r])]
        ax.plot([RING_MID[r] for r in ks], [sub[r] for r in ks],
                "-o", ms=2.8, lw=1.0, color=c,
                label=k.replace("p", "$^{+}$") if k == "LRRC15p" else k)
    ax.set_xscale("log")
    ax.set_xticks([25, 75, 150, 300, 600])
    ax.set_xticklabels(["25", "75", "150", "300", ">400"])
    ax.minorticks_off()
    ax.set_xlabel("distance to nearest nerve ($\\mu$m)")
    ax.set_ylabel("NNLS axis weight (mean)")
    ax.set_title("E  CAF-axis layering along the perineural axis", loc="left", fontsize=9)
    ax.legend(frameon=False, fontsize=6.5, ncol=2)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)


def panel_f(ax):
    B = pd.read_csv(PROC / "caf_axis_traj_bins.csv")
    g = B.groupby("bin").agg(n=("n_spot", "sum"), dpt=("dpt_median", "mean"))
    axes = ["myCAF", "iCAF", "apCAF", "LRRC15p"]
    colors = [C_ACC, C_MAIN, C_PURPLE, C_ALT]
    for k, c in zip(axes, colors):
        m = B.groupby("bin")[f"score_{k}"].mean()
        ax.plot(g["dpt"], m.reindex(g.index), "-o", ms=2.8, lw=1.0, color=c,
                label=k.replace("p", "$^{+}$") if k == "LRRC15p" else k)
    m = B.groupby("bin")["score_CD8_T"].mean()
    ax.plot(g["dpt"], m.reindex(g.index), "--", lw=1.2, color=C_GREY, label="CD8_T")
    m = B.groupby("bin")["score_T_lymphocyte"].mean()
    ax.plot(g["dpt"], m.reindex(g.index), ":", lw=1.2, color=C_ORANGE, label="T lymph.")
    ax.axhline(0, color=C_GREY, lw=0.5, ls=":")
    ax.set_xlabel("diffusion pseudotime (nerve-rooted)")
    ax.set_ylabel("axis score (mean $z$)")
    ax.set_title("F  Nerve-rooted trajectory: stroma axes decline,\n"
                 "immune readouts track the stroma", loc="left", fontsize=9)
    ax.legend(frameon=False, fontsize=6.5, ncol=2)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)


def main():
    fig = plt.figure(figsize=(7.09, 8.6))     # 180 mm 宽
    gs = fig.add_gridspec(3, 2, height_ratios=[1.15, 1, 1], hspace=0.5, wspace=0.42,
                          left=0.17, right=0.97, top=0.95, bottom=0.05)
    panel_a(fig.add_subplot(gs[0, 0]))
    panel_b(fig.add_subplot(gs[0, 1]))
    panel_c(fig.add_subplot(gs[1, 0]))
    panel_d(fig.add_subplot(gs[1, 1]))
    panel_e(fig.add_subplot(gs[2, 0]))
    panel_f(fig.add_subplot(gs[2, 1]))
    save(fig, "Fig4_three_main_lines")


if __name__ == "__main__":
    main()

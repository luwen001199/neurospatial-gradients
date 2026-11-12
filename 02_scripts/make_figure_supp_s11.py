# -*- coding: utf-8 -*-
"""SuppFigS11 —— 空间先验模块在公共 bulk 队列中的外推（下游链：纤维化 / 耐药 / 免疫）。

panels
------
S10A 模块对 CAF/ECM 的秩正交化：被纤维化解释的方差百分比（5 队列）
S10B 三级衰减 slope chart：PNI-stroma ~ 免疫读数，raw → +CAF → +CAF+purity
S10C pCR 森林图：原始 vs 正交化残差（3 队列 × 2 模块）
S10D 8 个基质程序的 OS 森林：原始（空心）vs 残差（实心），TCGA-BRCA
S10E SPARC/POSTN 残差上三分位 KM（TCGA-BRCA）
S10F 随机基因集基线：观测 HR 与 510 组 size-matched 随机集分布
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
from nsgrad import config as CFG     # noqa: E402
from nsgrad import viz as nviz       # noqa: E402

TAB = CFG.ROOT / "03_results" / "tables" / "bulk_extrapolation"
FIG = CFG.DIR_FIG

C_OBS = "#C0392B"
C_ALT = "#2E86C1"
C_GREY = "#7F8C8D"
C_LIGHT = "#D5D8DC"

COHORT_SHORT = {
    "TCGA-BRCA": "TCGA", "TCGA-TNBC": "TCGA-TNBC",
    "GSE163882": "GSE163882", "GSE163882-TNBC": "GSE163882-TN",
    "GSE25065": "GSE25065",
}
MOD_LABEL = {
    "Perineural_stroma": "perineural-stroma (60)",
    "Neural_intrinsic": "neural-intrinsic (39)",
    "NPP58": "NPP (58)",
    "LRRC15_myCAF": "LRRC15$^{+}$ myCAF (8)",
    "CAF_ECM": "CAF/ECM (14)",
    "SPARC_POSTN_CORE": "SPARC/POSTN core (8)",
    "COLLAGEN_CORE": "collagen core (7)",
    "COLLAGEN_MATURATION": "collagen maturation (12)",
    "MECHANOTRANSDUCTION": "mechanotransduction (8)",
    "EMT_CORE": "EMT core (9)",
    "TGFb_PROGRAM": "TGF-$\\beta$ program (8)",
}
TGT_LABEL = {
    "CD8_cytotox": "CD8 cytotoxicity", "Treg": "Treg", "M2_macro": "M2 macrophage",
    "SPP1_macro": "SPP1$^{+}$ macrophage", "TGFb_core": "TGF-$\\beta$ core",
    "CXCL12_axis": "CXCL12 axis", "LGALS1_axis": "LGALS1 axis",
    "IFN_gamma": "IFN-$\\gamma$", "TLS": "TLS", "MHC_II": "MHC class II",
}


def panel_tag(ax, s):
    ax.text(-0.16, 1.06, s, transform=ax.transAxes, fontsize=9, weight="bold",
            va="bottom", ha="left")


# ---------------------------------------------------------------- S10A
def s10a(ax):
    d = pd.read_csv(TAB / "stage7_orthogonality.csv")
    order = ["Perineural_stroma", "LRRC15_myCAF", "SPARC_POSTN_CORE",
             "COLLAGEN_MATURATION", "Neural_intrinsic", "NPP58"]
    order = [m for m in order if m in set(d["module"])]
    rng = np.random.default_rng(7)
    for i, m in enumerate(order):
        v = d[d["module"] == m]["pct_explained"].to_numpy()
        jit = rng.uniform(-0.14, 0.14, len(v))
        ax.scatter(v, i + jit, s=13, facecolor="white", edgecolor=C_ALT,
                   linewidth=0.8, zorder=3)
        ax.plot([np.median(v)] * 2, [i - 0.26, i + 0.26], color=C_OBS, lw=1.6, zorder=4)
        ax.text(101, i, f"{np.median(v):.0f}%", va="center", fontsize=6.5, color="#333")
    ax.axvline(50, color=C_LIGHT, lw=0.9, ls="--", zorder=1)
    ax.set_yticks(range(len(order)))
    ax.set_yticklabels([MOD_LABEL.get(m, m) for m in order], fontsize=6.5)
    ax.invert_yaxis()
    ax.set_xlim(0, 112)
    ax.set_xticks([0, 25, 50, 75, 100])
    ax.set_xlabel("% of variance explained by CAF/ECM (rank space)", fontsize=7)
    ax.set_title("Modules are mostly desmoplasia", fontsize=8.5, loc="left")
    ax.text(0.5, 0.97, "dots = 5 bulk cohorts; bar = median",
            transform=ax.transAxes, ha="center", va="top", fontsize=5.6, color="#666")
    panel_tag(ax, "A")


# ---------------------------------------------------------------- S10B
def s10b(ax):
    d = pd.read_csv(TAB / "stage7_attenuation.csv")
    d = d[d["module"] == "Perineural_stroma"]
    tgts = ["TGFb_core", "M2_macro", "SPP1_macro", "MHC_II", "Treg",
            "TLS", "IFN_gamma", "CXCL12_axis", "CD8_cytotox"]
    tgts = [t for t in tgts if t in set(d["target"])]
    cols = ["rho_raw", "rho_caf", "rho_caf_purity"]
    xs = np.arange(len(cols))
    med = d.groupby("target")[cols].median()
    for t in tgts:
        sub = d[d["target"] == t]
        for _, r in sub.iterrows():
            ax.plot(xs, [r[c] for c in cols], color=C_LIGHT, lw=0.7, zorder=1)
    for t in tgts:
        y = [med.loc[t, c] for c in cols]
        lw = 2.0 if t in ("TGFb_core", "M2_macro") else 1.1
        col = C_OBS if t == "TGFb_core" else (C_ALT if t == "M2_macro" else C_GREY)
        ax.plot(xs, y, color=col, lw=lw, marker="o", ms=3.4, zorder=3,
                label=TGT_LABEL.get(t, t) if lw > 1.5 else None)
        pass
    # 右侧标签防重叠：按末值排序后强制最小间距
    ends = sorted([(med.loc[t, cols[-1]], t) for t in tgts])
    placed = []
    for yv, t in ends:
        yv2 = max(yv, placed[-1] + 0.045) if placed else yv
        placed.append(yv2)
        col, lw = C_GREY, 1.1
        if t == "TGFb_core":
            col, lw = C_OBS, 2.0
        elif t == "M2_macro":
            col, lw = C_ALT, 2.0
        ax.text(2.12, yv2, TGT_LABEL.get(t, t), fontsize=5.6, va="center",
                color=col if lw > 1.5 else "#555")
    ax.axhline(0, color="#999", lw=0.7, ls=":")
    ax.set_xticks(xs)
    ax.set_xticklabels(["raw", "+ CAF/ECM", "+ purity"], fontsize=6.5)
    ax.set_xlim(-0.25, 2.95)
    ax.set_ylabel("Spearman $\\rho$ with perineural-stroma", fontsize=7)
    ax.set_title("TGF-$\\beta$ survives fibrosis adjustment", fontsize=8.5, loc="left")
    panel_tag(ax, "B")


# ---------------------------------------------------------------- S10C
def s10c(ax):
    d = pd.read_csv(TAB / "stage7_residual_pcr.csv")
    d = d[d["module"].isin(["Perineural_stroma", "SPARC_POSTN_CORE"])]
    cohorts = ["GSE163882", "GSE163882-TNBC", "GSE25065"]
    rows = []
    for m in ["Perineural_stroma", "SPARC_POSTN_CORE"]:
        for c in cohorts:
            for sc in ["raw", "resid"]:
                sub = d[(d["module"] == m) & (d["cohort"] == c) & (d["scale"] == sc)]
                if len(sub):
                    rows.append(sub.iloc[0])
    y = 0
    yt, yl = [], []
    for r in rows:
        col = C_OBS if r["scale"] == "resid" else C_GREY
        mk = "s" if r["module"] == "Perineural_stroma" else "o"
        ax.plot([r["ci_lo"], r["ci_hi"]], [y, y], color=col, lw=0.9)
        ax.scatter([r["OR_per_SD"]], [y], s=16, marker=mk, color=col, zorder=3)
        p = r["p"]
        star = "***" if p < 0.001 else "**" if p < 0.01 else "*" if p < 0.05 else ""
        ax.text(3.05, y, f"{r['OR_per_SD']:.2f}{star}", fontsize=5.6, va="center", color=col)
        yt.append(y)
        yl.append(f"{COHORT_SHORT.get(r['cohort'], r['cohort'])} {r['scale']}")
        y += 1
        if r["module"] == "Perineural_stroma" and r["cohort"] == "GSE25065":
            y += 0.6
    ax.axvline(1, color="#555", lw=0.8, ls="--")
    ax.set_yticks(yt)
    ax.set_yticklabels(yl, fontsize=5.8)
    ax.invert_yaxis()
    ax.set_xscale("log")
    ax.set_xlim(0.45, 4.2)
    ax.set_xticks([0.5, 1, 2, 4])
    ax.set_xticklabels(["0.5", "1", "2", "4"])
    ax.minorticks_off()
    ax.set_xlabel("OR per SD for pCR  (squares = perineural-stroma, circles = SPARC/POSTN)",
                  fontsize=6.5)
    ax.set_title("Residual score: consistent direction, weak effect", fontsize=8.5, loc="left")
    panel_tag(ax, "C")


# ---------------------------------------------------------------- S10D
def s10d(ax):
    d = pd.read_csv(TAB / "stage8_cox_models.csv")
    d = d[(d["cohort"] == "TCGA-BRCA") & (d["model"].isin(["M0_raw", "M3_resid"]))]
    mods = ["SPARC_POSTN_CORE", "COLLAGEN_CORE", "COLLAGEN_MATURATION",
            "MECHANOTRANSDUCTION", "EMT_CORE", "TGFb_PROGRAM", "CAF_ECM",
            "Perineural_stroma"]
    for i, m in enumerate(mods):
        for sc, dy, col, mk in (("raw", -0.13, C_GREY, "o"), ("resid", 0.13, C_OBS, "s")):
            sub = d[(d["module"] == m) & (d["scale"] == sc) & (d["model"] == ("M0_raw" if sc == "raw" else "M3_resid"))]
            if not len(sub):
                continue
            r = sub.iloc[0]
            ax.plot([r["lo"], r["hi"]], [i + dy] * 2, color=col, lw=0.9)
            ax.scatter([r["HR"]], [i + dy], s=15, marker=mk, color=col, zorder=3,
                       facecolor=col if sc == "resid" else "white", edgecolor=col)
            if sc == "resid" and r["p"] < 0.05:
                ax.text(r["hi"] + 0.05, i + dy, f"p={r['p']:.3f}", fontsize=5.4, va="center", color=C_OBS)
    ax.axvline(1, color="#555", lw=0.8, ls="--")
    ax.set_yticks(range(len(mods)))
    ax.set_yticklabels([MOD_LABEL.get(m, m) for m in mods], fontsize=6.2)
    ax.invert_yaxis()
    ax.set_xscale("log")
    ax.set_xlim(0.55, 2.4)
    ax.set_xticks([0.6, 1.0, 1.5, 2.0])
    ax.set_xticklabels(["0.6", "1.0", "1.5", "2.0"])
    ax.set_xlabel("HR per SD for OS  (open circle = raw, filled square = adjusted residual)",
                  fontsize=6.5)
    ax.set_title("Only the SPARC/POSTN program predicts OS", fontsize=8.5, loc="left")
    panel_tag(ax, "D")


# ---------------------------------------------------------------- S10E
def s10e(ax):
    d = pd.read_csv(TAB / "stage9_km_refined.csv")
    d = d[(d["cohort"] == "TCGA-BRCA") & (d["panel"] == "T3_vs_rest")]
    for lab, col, ls in (("rest", C_GREY, "-"), ("high", C_OBS, "-")):
        s = d[d["group"] == lab].sort_values("time_days")
        n = int(s["n"].iloc[0])
        ev = int(s["events"].iloc[0])
        p = s["logrank_p"].iloc[0]
        ax.step(s["time_days"] / 365.25, s["surv"], where="post", color=col, lw=1.4,
                label=f"TCGA-BRCA {'upper tertile' if lab=='high' else 'lower 2/3'} (n={n}, d={ev})")
    ax.set_xlim(0, 15)
    ax.set_ylim(0.4, 1.02)
    ax.set_xlabel("Time (years)", fontsize=7)
    ax.set_ylabel("Overall survival", fontsize=7)
    ax.legend(fontsize=5.8, frameon=False, loc="lower left")
    ax.set_title(f"Log-rank P = {p:.2f}", fontsize=8.5, loc="left")
    panel_tag(ax, "E")


# ---------------------------------------------------------------- S10F
def s10f(ax):
    d = pd.read_csv(TAB / "stage9_random_baseline.csv")
    cohorts = ["TCGA-BRCA", "TCGA-TNBC"]
    for i, c in enumerate(cohorts):
        sub = d[d["cohort"] == c]
        if not len(sub):
            continue
        r = sub.iloc[0]
        ax.plot([r["rand_p2.5"], r["rand_p97.5"]], [i, i], color=C_LIGHT, lw=8,
                solid_capstyle="butt", zorder=2)
        ax.scatter([r["rand_median_HR"]], [i], s=26, marker="|", color=C_GREY, zorder=3)
        ax.scatter([r["obs_HR"]], [i], s=34, marker="D", color=C_OBS, zorder=4)
        ax.text(r["rand_p97.5"] + 0.04, i + 0.22,
                f"empirical P = {r['empirical_p']:.2f}", fontsize=6, color=C_OBS)
        ax.text(0.38, i + 0.30, f"n={int(r['n_random'])} random 8-gene sets",
                fontsize=5.5, color="#666")
    ax.axvline(1, color="#555", lw=0.8, ls="--")
    ax.set_yticks(range(len(cohorts)))
    ax.set_yticklabels([f"{c} (n={int(d[d['cohort']==c]['n_random'].iloc[0]) if len(d[d['cohort']==c]) else 0} sets)"
                        for c in cohorts], fontsize=6.2)
    ax.set_xlim(0.35, 2.1)
    ax.set_xlabel("HR per SD for OS after residualisation", fontsize=7)
    ax.set_title("Observed HR vs size-matched random sets", fontsize=8.5, loc="left")
    panel_tag(ax, "F")


def main():
    fig = plt.figure(figsize=(7.2, 8.4))
    gs = fig.add_gridspec(3, 2, hspace=0.62, wspace=0.42,
                          left=0.14, right=0.97, top=0.95, bottom=0.06)
    s10a(fig.add_subplot(gs[0, 0]))
    s10b(fig.add_subplot(gs[0, 1]))
    s10c(fig.add_subplot(gs[1, 0]))
    s10d(fig.add_subplot(gs[1, 1]))
    s10e(fig.add_subplot(gs[2, 0]))
    s10f(fig.add_subplot(gs[2, 1]))
    out = nviz.save(fig, "SuppFigS11_bulk_extrapolation")
    print(f"saved -> {out}")


if __name__ == "__main__":
    main()

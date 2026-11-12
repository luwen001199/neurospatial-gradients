# -*- coding: utf-8 -*-
"""SuppFigS12 —— 免疫区室、中介链与跨癌种复现（Stage 10–12）。

panels
------
S12A 随机基因集校准的正交性：stroma vs neural 在每个队列的 z（同一零假设下）
S12B 区室富集分的组成效度：与配对单细胞真实细胞占比的 Spearman（n=26）
S12C 中介分析森林：残差 → 中介 → CD8 细胞毒的间接效应 ab（两队列）
S12D 跨癌种主链：CAF/ECM 解释比例、残差~TGF-β、残差~CD8、残差对 OS 的 HR
"""
from __future__ import annotations

import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, r"D:/嗜神经浸润/00_pipeline")
from nsgrad import config as CFG     # noqa: E402
from nsgrad import viz as nviz       # noqa: E402

TAB = CFG.ROOT / "03_results" / "tables" / "bulk_extrapolation"

C_OBS = "#C0392B"
C_ALT = "#2E86C1"
C_GREY = "#7F8C8D"
C_LIGHT = "#D5D8DC"

SHORT = {"TCGA-BRCA": "TCGA-BRCA", "TCGA-TNBC": "TCGA-TNBC", "TCGA-PAAD": "TCGA-PAAD",
         "TCGA-HNSC": "TCGA-HNSC", "METABRIC": "METABRIC", "GSE163882": "GSE163882",
         "GSE25065": "GSE25065"}
MED_LABEL = {"TGF-β program": "TGF-$\\beta$ program",
             "macrophage enrichment": "macrophage",
             "M2 marker score": "M2 score", "Treg marker score": "Treg score",
             "CD8 T enrichment": "CD8 T"}


def tag(ax, s):
    ax.text(-0.20, 1.08, s, transform=ax.transAxes, fontsize=9, weight="bold",
            va="bottom", ha="left")


# ---------------------------------------------------------------- A
def s14a(ax):
    d = pd.read_csv(TAB / "stage12_orthogonality_random.csv")
    cohorts = ["TCGA-BRCA", "TCGA-TNBC", "TCGA-PAAD", "TCGA-HNSC",
               "METABRIC", "GSE163882", "GSE25065"]
    cohorts = [c for c in cohorts if c in set(d["cohort"])]
    h = 0.36
    for i, c in enumerate(cohorts):
        for j, (mod, col, lab, off) in enumerate([
                ("Perineural_stroma", C_OBS, "perineural-stroma", +h / 2),
                ("Neural_intrinsic", C_ALT, "neural-intrinsic", -h / 2)]):
            r = d[(d.cohort == c) & (d.module == mod)]
            if not len(r):
                continue
            z = float(r["z_vs_null"].iloc[0])
            ax.barh(i + off, z, height=h * 0.86, color=col, alpha=0.85,
                    edgecolor="white", linewidth=0.5)
            ax.text(z + (0.4 if z >= 0 else -0.4), i + off, f"{z:+.1f}",
                    va="center", ha="left" if z >= 0 else "right",
                    fontsize=5.8, color=col)
            if i == 0:
                ax.barh([], [], color=col, label=lab)
    ax.axvline(0, color="#333", lw=0.9)
    ax.set_yticks(range(len(cohorts)))
    ax.set_yticklabels([f"{SHORT[c]} (n={int(d[d.cohort==c]['n_sample'].iloc[0])})"
                        for c in cohorts], fontsize=6.4)
    ax.invert_yaxis()
    ax.set_xlim(-14.5, 12)
    ax.set_xlabel("z vs 1,000 size-matched random sets\n(negative = further from desmoplasia than chance)",
                  fontsize=6.4)
    ax.set_title("Orthogonality must be read against a null", fontsize=8.4, loc="left")
    ax.legend(fontsize=6.2, frameon=False, loc="lower left", bbox_to_anchor=(0.01, 0.02),
              handlelength=1.2, borderaxespad=0.1)
    tag(ax, "A")


# ---------------------------------------------------------------- B
def s14b(ax):
    d = pd.read_csv(TAB / "stage10_score_validation.csv")
    d = d.sort_values("rho_score_vs_true_fraction")
    y = np.arange(len(d))
    ax.barh(y, d["rho_score_vs_true_fraction"], height=0.62, color=C_ALT, alpha=0.85,
            edgecolor="white", linewidth=0.5)
    for i, (_, r) in enumerate(d.iterrows()):
        ax.text(r["rho_score_vs_true_fraction"] + 0.012, i,
                f"{r['rho_score_vs_true_fraction']:.2f} (n={int(r['n_genes'])})",
                va="center", fontsize=5.7, color="#333")
    ax.set_yticks(y)
    ax.set_yticklabels([c.replace("_", " ") for c in d["compartment"]], fontsize=6.4)
    ax.set_xlim(0, 1.0)
    ax.set_xlabel("Spearman ρ vs true cell fraction (pseudobulk, 26 patients; n = markers)",
                  fontsize=6.4)
    ax.set_title("Compartment scores track real composition", fontsize=8.4, loc="left")
    tag(ax, "B")


# ---------------------------------------------------------------- C
def s14c(ax):
    d = pd.read_csv(TAB / "stage10_mediation.csv")
    d = d[d.exposure.str.startswith("perineural")]
    meds = ["TGF-β program", "macrophage enrichment", "M2 marker score",
            "Treg marker score", "CD8 T enrichment"]
    meds = [m for m in meds if m in set(d["mediator"])]
    cohorts = ["TCGA-BRCA", "TCGA-TNBC"]
    cols = {cohorts[0]: C_OBS, cohorts[1]: C_ALT}
    for j, m in enumerate(meds):
        for k, co in enumerate(cohorts):
            r = d[(d.cohort == co) & (d.mediator == m)]
            if not len(r):
                continue
            r = r.iloc[0]
            y = j + (0.18 if k == 0 else -0.18)
            ax.plot([r["ab_lo"], r["ab_hi"]], [y, y], color=cols[co], lw=1.4,
                    solid_capstyle="round")
            ax.plot(r["ab"], y, "o", ms=4.2, color=cols[co],
                    mfc=cols[co] if k == 0 else "white", mew=1.2)
            ax.text(r["ab_hi"] + 0.008, y, f"{r['prop_mediated']*100:.0f}%",
                    va="center", fontsize=5.7, color=cols[co])
    ax.axvline(0, color="#333", lw=0.9)
    ax.set_yticks(range(len(meds)))
    ax.set_yticklabels([MED_LABEL.get(m, m) for m in meds], fontsize=6.4)
    ax.invert_yaxis()
    ax.set_xlim(-0.02, 0.52)
    ax.set_xlabel("indirect effect ab on CD8 cytotoxicity (95% CI)\n"
                  "filled = TCGA-BRCA, open = TCGA-TNBC; % = proportion mediated",
                  fontsize=6.4)
    ax.set_title("No suppression to mediate: all paths positive", fontsize=8.4, loc="left")
    tag(ax, "C")


# ---------------------------------------------------------------- D
def s14d(ax):
    d = pd.read_csv(TAB / "stage12_pancancer.csv")
    cohorts = ["TCGA-BRCA", "TCGA-TNBC", "GSE163882", "GSE25065",
               "METABRIC", "TCGA-PAAD", "TCGA-HNSC"]
    d = d.set_index("cohort").reindex([c for c in cohorts if c in set(d["cohort"])])
    cols = [("r2_stroma_on_CAF", "CAF/ECM\nrank $R^2$", 100),
            ("rho_resid_TGFb", "resid ~\nTGF-$\\beta$", 1),
            ("rho_resid_CD8", "resid ~\nCD8", 1),
            ("HR_resid", "resid vs OS\n(HR/SD)", 1)]
    M = np.full((len(d), len(cols)), np.nan)
    for j, (c, _, _) in enumerate(cols):
        if c in d:
            M[:, j] = pd.to_numeric(d[c], errors="coerce").to_numpy()
    ax.imshow(M.T, aspect="auto", cmap="RdBu_r", vmin=-1, vmax=1)
    for j in range(len(cols)):
        for i in range(len(d)):
            v = M[i, j]
            if not np.isfinite(v):
                ax.text(i, j, "n/a", ha="center", va="center", fontsize=5.4, color="#888")
                continue
            txt = f"{v*100:.0f}" if j == 0 else f"{v:+.2f}"
            ax.text(i, j, txt, ha="center", va="center", fontsize=6.0,
                    color="white" if abs(v) > 0.62 else "#222")
    ax.set_xticks(range(len(d)))
    ax.set_xticklabels([f"{SHORT.get(c, c)} ({int(d.loc[c,'n'])})" for c in d.index],
                       fontsize=5.6, rotation=32, ha="right", rotation_mode="anchor")
    ax.set_yticks(range(len(cols)))
    ax.set_yticklabels([c[1] for c in cols], fontsize=6.2)
    ax.set_title("Same three features in pancreatic and head–neck cancer",
                 fontsize=8.4, loc="left")
    tag(ax, "D")


def main():
    fig = plt.figure(figsize=(7.2, 6.6))
    gs = fig.add_gridspec(2, 2, hspace=0.52, wspace=0.52,
                          left=0.155, right=0.975, top=0.915, bottom=0.085)
    s14a(fig.add_subplot(gs[0, 0]))
    s14b(fig.add_subplot(gs[0, 1]))
    s14c(fig.add_subplot(gs[1, 0]))
    s14d(fig.add_subplot(gs[1, 1]))
    out = nviz.save(fig, "SuppFigS12_immune_pancancer")
    print(f"saved -> {out}")


if __name__ == "__main__":
    main()

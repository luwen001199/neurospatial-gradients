# -*- coding: utf-8 -*-
"""SuppFigS14 —— Xenium 单细胞分辨率神经邻近验证（Stage 14, GSE350571）。

panels
------
S14A 神经核点与分区规模: 每例核点数(log) + 近端/远端细胞数 (nerve-negative 例如实标注)
S14B 成纤维细胞模块评分 prox(≤150µm) vs dist(≥300µm) 哑铃图 (4 神经阳性例 × 4 模块)
S14C CD8 构成比 prox/dist + torus-shift 置换 P (方向不一致展示)
S14D 全细胞类型 prox/dist 构成比热图 (log2 ratio, 4 例 × 6 类型)
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

C_SIG = "#C0392B"
C_NS = "#2E86C1"
C_GREY = "#7F8C8D"
C_LIGHT = "#D5D8DC"

MODULES = ["Perineural_stroma", "LRRC15_myCAF", "CAF_ECM", "Neural_intrinsic"]
MOD_NICE = {"Perineural_stroma": "Perineural\nstroma", "LRRC15_myCAF": "LRRC15\nmyCAF",
            "CAF_ECM": "CAF/ECM", "Neural_intrinsic": "Neural\nintrinsic"}
CTYPES = ["myCAF", "fibroblast", "endothelial", "CD8_T", "macrophage", "B_cell"]
CT_NICE = {"myCAF": "myCAF", "fibroblast": "fibroblast", "endothelial": "endothelial",
           "CD8_T": "CD8 T", "macrophage": "macrophage", "B_cell": "B cell"}
NERVE_POS = ["P22", "P26", "P29", "P32"]


def tag(ax, s):
    ax.text(-0.22, 1.08, s, transform=ax.transAxes, fontsize=9, weight="bold",
            va="bottom", ha="left")


# ---------------------------------------------------------------- A
def s16a(ax):
    d = pd.read_csv(TAB / "stage14_xenium_nerve_summary.csv")
    d = d.set_index("sample")
    order = ["P17", "P22", "P24", "P26", "P29", "P32"]
    d = d.reindex(order)
    x = np.arange(len(d))
    npts = d["n_nerve_pts"].fillna(0).to_numpy(float)
    colors = [C_SIG if s in NERVE_POS else C_LIGHT for s in d.index]
    ax.bar(x, np.log10(npts + 1), color=colors, edgecolor="white", lw=0.5, width=0.62)
    for i, s in enumerate(d.index):
        v = npts[i]
        lab = f"{int(v)}" if v > 0 else "0"
        ax.text(i, np.log10(v + 1) + 0.06, lab, ha="center", fontsize=5.8,
                color="#333")
        if s not in NERVE_POS:
            ax.text(i, 0.12, "nerve-\nnegative", ha="center", fontsize=5.0,
                    color=C_GREY, rotation=90, va="bottom")
    ax.set_xticks(x)
    ax.set_xticklabels(d.index, fontsize=6.4)
    ax.set_ylabel("axonal-anchored nerve kernels\n(log$_{10}$, ≥1 shown)", fontsize=6.4)
    ax.set_title("Nerve kernels per Xenium section\n(MPZ-anchoring neutralised in P26)",
                 fontsize=8.0, loc="left")
    ax.tick_params(axis="y", labelsize=6.0)
    tag(ax, "A")


# ---------------------------------------------------------------- B
def s16b(ax):
    f = pd.read_csv(TAB / "stage14_xenium_fb_scores.csv")
    piv = f.pivot_table(index=["module", "sample"], columns="zone",
                        values="mean").reset_index()
    ylabels, y, pairs = [], [], []
    ypos = 0.0
    for m in MODULES:
        g = piv[piv.module == m].set_index("sample").reindex(NERVE_POS)
        for s in NERVE_POS:
            r = g.loc[s]
            pairs.append((ypos, r["prox"], r["dist"], m, s, int(
                f[(f.module == m) & (f['sample'] == s) & (f.zone == "prox")]["n"].iloc[0])))
            ylabels.append(f"{s}" if s == NERVE_POS[0] or True else "")
            ypos += 1
        ylabels[-1] = MOD_NICE[m]
        ypos += 0.9
    for (yy, prox, dist, m, s, n) in pairs:
        c = C_SIG if prox > dist else C_NS
        ax.plot([dist, prox], [yy, yy], color=C_LIGHT, lw=2.4, zorder=1)
        ax.plot(dist, yy, "o", color=C_GREY, ms=5, zorder=2)
        ax.plot(prox, yy, "o", color=c, ms=5.5, zorder=3)
        ax.text(1.88, yy, f"{n:,}", fontsize=4.8, color="#777", va="center")
    ax.set_yticks([p[0] for p in pairs])
    ax.set_yticklabels([p[3][:1] + p[4][-2:] if False else "" for p in pairs])
    # 模块标签放在组间隙
    for m in MODULES:
        g = piv[piv.module == m].set_index("sample").reindex(NERVE_POS)
        ys = [p[0] for p in pairs if p[3] == m]
        ax.text(-0.16, np.mean(ys), MOD_NICE[m], fontsize=6.0, ha="right",
                va="center", transform=ax.get_yaxis_transform())
        for k, s in enumerate(NERVE_POS):
            yy = [p[0] for p in pairs if p[3] == m and p[4] == s][0]
            ax.text(-0.035, yy, s, fontsize=5.2, ha="right", va="center",
                    color="#555", transform=ax.get_yaxis_transform())
    ax.set_xlabel("mean-z module score in fibroblasts", fontsize=6.4)
    ax.set_xlim(0.2, 1.95)
    ax.set_ylabel("")
    ax.set_yticklabels([])
    ax.set_title("Nerve-proximal (red) vs distal (grey) fibroblast programs\n"
                 "red/blue: proximal higher/lower; grey number = proximal n", fontsize=7.6, loc="left")
    ax.tick_params(axis="x", labelsize=6.0)
    tag(ax, "B")


# ---------------------------------------------------------------- C
def s16c(ax):
    d = pd.read_csv(TAB / "stage14_xenium_nerve_summary.csv")
    d = d[d['sample'].isin(NERVE_POS)].set_index("sample").reindex(NERVE_POS)
    x = np.arange(len(d))
    ratio = d["CD8_T_prox_dist_ratio"].to_numpy(float)
    for i, (r, s) in enumerate(zip(ratio, d.index)):
        p = d.loc[s, "CD8_T_permP"]
        c = C_SIG if p < 0.05 else C_NS
        ax.bar(i, r, color=c, alpha=0.85, width=0.6, edgecolor="white")
        ax.text(i, r + 0.03 * max(ratio), f"P={p:.2f}", ha="center", fontsize=5.6,
                color=c)
    ax.axhline(1.0, color="#333", lw=1.0, ls="--")
    ax.text(len(d) - 0.45, 1.03, "parity", fontsize=5.6, color="#333", ha="right")
    ax.set_xticks(x)
    ax.set_xticklabels(d.index, fontsize=6.4)
    ax.set_ylabel("CD8 T proximal/distal fraction", fontsize=6.4)
    ax.set_title("CD8 composition: directionally inconsistent\n"
                 "(P22 enriched, P29 depleted, P26/P32 null)", fontsize=8.0, loc="left")
    ax.tick_params(axis="y", labelsize=6.0)
    ax.set_ylim(0, max(ratio) * 1.25)
    tag(ax, "C")


# ---------------------------------------------------------------- D
def s16d(ax):
    d = pd.read_csv(TAB / "stage14_xenium_nerve_summary.csv")
    d = d[d['sample'].isin(NERVE_POS)].set_index("sample").reindex(NERVE_POS)
    M = np.full((len(CTYPES), len(d)), np.nan)
    for j, s in enumerate(d.index):
        for i, ct in enumerate(CTYPES):
            v = d.loc[s, f"{ct}_prox_dist_ratio"]
            if np.isfinite(v) and v > 0:
                M[i, j] = np.log2(v)
    im = ax.imshow(M, aspect="auto", cmap="RdBu_r", vmin=-1.2, vmax=1.2)
    for j in range(M.shape[1]):
        for i in range(M.shape[0]):
            v = M[i, j]
            if not np.isfinite(v):
                ax.text(j, i, "n/a", ha="center", va="center", fontsize=5.2,
                        color="#888")
                continue
            ax.text(j, i, f"{v:+.2f}", ha="center", va="center", fontsize=5.4,
                    color="white" if abs(v) > 0.7 else "#222")
    ax.set_xticks(range(M.shape[1]))
    ax.set_xticklabels(d.index, fontsize=6.4)
    ax.set_yticks(range(len(CTYPES)))
    ax.set_yticklabels([CT_NICE[c] for c in CTYPES], fontsize=6.2)
    ax.set_title("log$_2$ proximal/distal composition ratio\n"
                 "(fractions sum to 1; epithelium proximal-enriched)", fontsize=8.0,
                 loc="left")
    cb = plt.colorbar(im, ax=ax, fraction=0.035, pad=0.02)
    cb.ax.tick_params(labelsize=5.4)
    tag(ax, "D")


def main():
    fig = plt.figure(figsize=(7.2, 6.8))
    gs = fig.add_gridspec(2, 2, hspace=0.55, wspace=0.45,
                          left=0.30, right=0.965, top=0.915, bottom=0.075)
    s16a(fig.add_subplot(gs[0, 0]))
    s16b(fig.add_subplot(gs[0, 1]))
    s16c(fig.add_subplot(gs[1, 0]))
    s16d(fig.add_subplot(gs[1, 1]))
    out = nviz.save(fig, "SuppFigS14_xenium_validation")
    print(f"saved -> {out}")


if __name__ == "__main__":
    main()

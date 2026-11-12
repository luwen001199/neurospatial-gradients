# -*- coding: utf-8 -*-
"""SuppFigS13 —— 感觉神经/CGRP 轴近神经梯度 + TLS 近/远功能比（Stage 21，与 Cell 2026 对话）。

panels
------
S13A 轴读数距离梯度 forest: rho_mean + 95%CI, 灰菱形 = 该读数尺寸桶的随机基线;
      红 = 超基线 (p_bucket<0.05)
S13B TLS+ spots 内近(≤150µm)-远(≥300µm)功能比 forest: Δ mean + 95%CI, 0 参考线
S13C 样本异质性: Sensory_nf/S100B/CGRP_axis 的 per-sample rho 分布 (读出 CI 跨 0 的原因)
S13D TLS 功能 Δ 的 per-sample 分布 (与 B 互补)
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

TAB = CFG.ROOT / "03_results" / "tables" / "spatial_tlscaf"

C_SIG = "#C0392B"      # 超基线 / 显著
C_NS = "#2E86C1"       # 不超基线
C_GREY = "#7F8C8D"
C_ALT = "#8E44AD"
C_LIGHT = "#D5D8DC"
C_LIGHT = "#D5D8DC"

NICE = {
    "MOD|Sensory_nf": ("Sensory markers\n(TRPV1/SCN9A/SCN10A/S100B)", 4),
    "GENE|S100B": ("S100B", 27),
    "MOD|CGRP_axis": ("CGRP axis\n(CALCA/CALCRL/RAMP1/CALCR)", 18),
    "MOD|CGRP_core": ("CGRP ligand–receptor\ncore (CALCA/RAMP1)", 26),
    "GENE|CALCRL": ("CALCRL", 23),
    "GENE|TRPV1": ("TRPV1", 12),
    "GENE|RAMP1": ("RAMP1", 26),
}
PANEL_NICE = {
    "cyt": "Cytotoxicity",
    "exh": "Exhaustion",
    "tgfb": "TGF-β response",
    "cyt_minus_exh": "CYT − EXH",
}


def tag(ax, s):
    ax.text(-0.22, 1.08, s, transform=ax.transAxes, fontsize=9, weight="bold",
            va="bottom", ha="left")


# ---------------------------------------------------------------- A
def s15a(ax):
    d = pd.read_csv(TAB / "stage21_axis_gradient_meta.corrected.csv")
    d = d[d.name.isin(NICE)].copy()
    d["lab"] = d["name"].map(lambda k: NICE[k][0])
    d["nsamp"] = d["name"].map(lambda k: NICE[k][1])
    d = d.sort_values("rho_mean")
    y = np.arange(len(d))
    sig = d["p_bucket"] < 0.05
    for i, r in enumerate(d.itertuples()):
        c = C_SIG if sig.iloc[i] else C_NS
        ax.plot([r.rho_lo, r.rho_hi], [i, i], color=c, lw=1.6, alpha=0.9)
        ax.plot(r.rho_mean, i, "o", color=c, ms=5.5)
        ax.plot(-r.baseline_rho, i, "D", color=C_GREY, ms=4.5, alpha=0.9)
        ax.text(r.rho_hi - 0.004, i + 0.30, f"n={r.n_samples}",
                fontsize=5.4, color="#555", va="bottom", ha="right")
        if sig.iloc[i]:
            ax.text(r.rho_mean - 0.018, i, f"P={r.p_bucket:.3f}",
                    fontsize=5.4, color=C_SIG, va="center", ha="right")
    ax.axvline(0, color="#333", lw=0.9)
    ax.set_yticks(y)
    ax.set_yticklabels(d["lab"], fontsize=6.0)
    ax.invert_yaxis()
    ax.set_xlabel("partial Spearman ρ vs distance (negative = nerve-proximal)",
                  fontsize=6.4)
    ax.set_title("Sensory readouts exceed the null;\nCGRP axis weakly; LR core does not",
                 fontsize=7.6, loc="left")
    ax.plot([], [], "o", color=C_SIG, label="above null (P<0.05)")
    ax.plot([], [], "o", color=C_NS, label="not above null")
    ax.plot([], [], "D", color=C_GREY, label="random baseline (|ρ|)")
    ax.legend(fontsize=5.2, frameon=False, loc="upper right", bbox_to_anchor=(1.02, 0.16))
    ax.tick_params(axis="x", labelsize=6.0)
    tag(ax, "A")


# ---------------------------------------------------------------- B
def s15b(ax):
    d = pd.read_csv(TAB / "stage21_tls_function_meta.csv")
    d = d.set_index("panel").reindex(["cyt", "exh", "tgfb", "cyt_minus_exh"]).dropna(
        subset=["mean"]).reset_index()
    for i, r in enumerate(d.itertuples()):
        c = C_SIG if r.p < 0.05 else C_NS
        ax.plot([r.lo, r.hi], [i, i], color=c, lw=1.8)
        ax.plot(r.mean, i, "s", color=c, ms=6)
        ptxt = f"P={r.p:.2f}" if r.p >= 0.01 else f"P={r.p:.3f}"
        ax.text(0.135, i, f"{ptxt} (n={r.n})", fontsize=5.6, color=c,
                va="center", ha="right")
    ax.axvline(0, color="#333", lw=0.9)
    ax.axvspan(-0.02, 0.02, color=C_LIGHT, alpha=0.35, zorder=0)
    ax.set_yticks(range(len(d)))
    ax.set_yticklabels([PANEL_NICE[p] for p in d["panel"]], fontsize=6.0)
    ax.invert_yaxis()
    ax.set_xlabel("near-minus-far Δ within TLS$^+$ spots (residualised mean-z)",
                  fontsize=6.4)
    ax.set_xlim(-0.16, 0.14)
    ax.set_title("Near-nerve TLS are not anergic;\nexhaustion lower (opposite to anergy)",
                 fontsize=7.6, loc="left")
    ax.tick_params(axis="x", labelsize=6.0)
    tag(ax, "B")


# ---------------------------------------------------------------- C
def s15c(ax):
    P = pd.read_csv(TAB / "stage21_axis_gradient_per_sample.csv.gz")
    groups = [("MOD|Sensory_nf", "Sensory markers"), ("GENE|S100B", "S100B"),
              ("MOD|CGRP_axis", "CGRP axis")]
    rng = np.random.default_rng(20261111)
    for k, (name, lab) in enumerate(groups):
        g = P[(P.name == name) & (P.kind != "random")]["rho"].dropna().to_numpy(float)
        x = rng.normal(k, 0.055, size=len(g))
        c = [C_SIG, C_GREY, C_ALT][k]
        ax.scatter(x, g, s=9, color=c, alpha=0.65, edgecolor="white", lw=0.3)
        ax.hlines(np.median(g), k - 0.16, k + 0.16, color="#222", lw=1.4)
        ax.text(k + 0.22, np.median(g) - 0.03, f"med {np.median(g):+.3f}\n(n={len(g)})",
                fontsize=5.4, va="top", color="#333")
    ax.axhline(0, color="#333", lw=0.9)
    ax.set_xticks(range(len(groups)))
    ax.set_xticklabels([g[1] for g in groups], fontsize=6.4)
    ax.set_ylabel("per-sample ρ vs distance", fontsize=6.4)
    ax.set_title("Sample heterogeneity behind the CI\n(CGRP axis spans both signs)",
                 fontsize=8.0, loc="left")
    ax.tick_params(axis="y", labelsize=6.0)
    ax.set_xlim(-0.55, len(groups) - 0.28)
    tag(ax, "C")


# ---------------------------------------------------------------- D
def s15d(ax):
    CP = pd.read_csv(TAB / "stage21_tls_function_per_sample.csv")
    panels = [("cyt", "CYT", C_ALT), ("exh", "EXH", C_SIG),
              ("tgfb", "TGF-β resp.", C_GREY)]
    rng = np.random.default_rng(20261112)
    for k, (col, lab, c) in enumerate(panels):
        v = CP[f"{col}_delta"].dropna().to_numpy(float)
        x = rng.normal(k, 0.06, size=len(v))
        ax.scatter(x, v, s=9, color=c, alpha=0.65, edgecolor="white", lw=0.3)
        if len(v) >= 3:
            ax.hlines(np.median(v), k - 0.17, k + 0.17, color="#222", lw=1.4)
        ax.text(k + 0.24, 0, f"n={len(v)}", fontsize=5.6, va="center", color="#333")
    ax.axhline(0, color="#333", lw=0.9)
    ax.set_xticks(range(len(panels)))
    ax.set_xticklabels([p[1] for p in panels], fontsize=6.4)
    ax.set_ylabel("per-section Δ (near − far)", fontsize=6.4)
    ax.set_title("Per-section TLS functional deltas:\nno anergy signature anywhere",
                 fontsize=8.0, loc="left")
    ax.tick_params(axis="y", labelsize=6.0)
    ax.set_xlim(-0.55, len(panels) - 0.25)
    tag(ax, "D")


def main():
    fig = plt.figure(figsize=(8.2, 6.8))
    gs = fig.add_gridspec(2, 2, hspace=0.62, wspace=0.52,
                          left=0.235, right=0.965, top=0.915, bottom=0.075)
    s15a(fig.add_subplot(gs[0, 0]))
    s15b(fig.add_subplot(gs[0, 1]))
    s15c(fig.add_subplot(gs[1, 0]))
    s15d(fig.add_subplot(gs[1, 1]))
    out = nviz.save(fig, "SuppFigS13_sensory_axis_tls")
    print(f"saved -> {out}")


if __name__ == "__main__":
    main()

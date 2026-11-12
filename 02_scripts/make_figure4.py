# -*- coding: utf-8 -*-
"""Figure 4 —— bulk 队列临床验证 (METABRIC, n≈1980)。

panels
------
4A KM 曲线: RFS 按 NPP 三分位
4B KM 曲线: OS 按 NPP 三分位
4C 森林图: RFS 各神经轴 / NPP 的多变量校正 HR
4D 森林图: OS
4E NPP 与通用 CAF/ECM 评分的散点 (展示 bulk 中的"稀释"效应)
4F KM: 按基质校正后残余 NPP 三分位 (RFS)
4G 置换检验: 随机等大基因集的 HR 零分布 vs 观测值
4H 各 PAM50/Claudin 亚型内的 NPP 效应量
"""
from __future__ import annotations
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, r"D:/嗜神经浸润/00_pipeline")
from nsgrad import config as CFG        # noqa: E402
from nsgrad import surv as nsurv       # noqa: E402

PROC = CFG.DIR_PROC
TAB = CFG.DIR_TAB
FIG = CFG.DIR_FIG
FIG.mkdir(parents=True, exist_ok=True)
DPI = 300
COL = {"T1": "#2980b9", "T2": "#7f8c8d", "T3": "#c0392b",
       "T1_low": "#2980b9", "T3_high": "#c0392b"}


def _km_ax(ax, d, tcol, ecol, gcol, title):
    for lev in [x for x in ["T1_low", "T1", "T2", "T3", "T3_high"] if x in set(d[gcol])]:
        gg = d[d[gcol] == lev]
        t, s, nr, ne = nsurv.km(gg[tcol], gg[ecol])
        ax.step(t, s, where="post", lw=1.6, color=COL.get(lev, "#555"),
                label=f"{lev} (n={len(gg)})")
        cens = gg[gg[ecol] == 0]
        if len(cens):
            yy = np.interp(cens[tcol], t, s)
            ax.plot(cens[tcol], yy, "|", color=COL.get(lev, "#555"), ms=3, alpha=0.35)
    ax.set_ylim(0, 1.02)
    ax.set_xlabel("Time (months)")
    ax.set_ylabel("Survival probability")
    ax.set_title(title, loc="left", fontweight="bold", fontsize=10)
    ax.legend(fontsize=7, frameon=True, loc="lower left")
    ax.grid(alpha=0.2, lw=0.4)


def fig4ab(d):
    fig, axes = plt.subplots(1, 2, figsize=(10.4, 4.3))
    _km_ax(axes[0], d, "RFS_MONTHS", "RFS_EVENT", "NPP_tertile",
           "A  Relapse-free survival by NPP tertile")
    _km_ax(axes[1], d, "OS_MONTHS", "OS_EVENT", "NPP_tertile",
           "B  Overall survival by NPP tertile")
    for ax, (tc, ec) in zip(axes, [("RFS_MONTHS", "RFS_EVENT"),
                                   ("OS_MONTHS", "OS_EVENT")]):
        lo = d[d["NPP_tertile"] == "T1_low"]; hi = d[d["NPP_tertile"] == "T3_high"]
        chi, p = nsurv.logrank(lo[tc], lo[ec], hi[tc], hi[ec])
        ax.text(0.97, 0.30, f"log-rank P = {p:.3g}", transform=ax.transAxes,
                ha="right", fontsize=8)
    fig.tight_layout()
    fig.savefig(FIG / "Fig4AB_km_npp_tertiles.png", dpi=DPI)
    fig.savefig(FIG / "Fig4AB_km_npp_tertiles.pdf")
    plt.close(fig)


def forest(ax, rows, title, xlabel="Adjusted HR per SD (95% CI)"):
    rows = [r for r in rows if np.isfinite(r[1])]
    y = np.arange(len(rows))[::-1]
    for i, (nm, hr, lo, hi, p) in enumerate(rows):
        yy = y[i]
        c = "#c0392b" if (np.isfinite(p) and p < 0.05) else "#34495e"
        ax.plot([lo, hi], [yy, yy], color=c, lw=1.4)
        ax.plot(hr, yy, "s", color=c, ms=5)
    ax.axvline(1, color="#888", lw=0.8, ls="--")
    ax.set_yticks(y)
    ax.set_yticklabels([r[0] for r in rows], fontsize=8)
    ax.set_xlabel(xlabel)
    ax.set_title(title, loc="left", fontweight="bold", fontsize=10)
    ax.grid(axis="x", alpha=0.2, lw=0.4)
    for i, (nm, hr, lo, hi, p) in enumerate(rows):
        ptxt = f"{p:.3g}" if np.isfinite(p) else "NA"
        ax.text(ax.get_xlim()[1], y[i],
                f" {hr:.3f} ({lo:.2f}-{hi:.2f}) P={ptxt}",
                va="center", fontsize=6.8)


def fig4cd(tab):
    fig, axes = plt.subplots(1, 2, figsize=(11.6, 5.0))
    for ax, oc in zip(axes, ["RFS", "OS"]):
        sub = tab[(tab["outcome"] == oc) &
                  (tab["model"].str.startswith("multivariable"))]
        rows = [(r["name"], r["hr"], r["lo"], r["hi"], r["p"])
                for _, r in sub.iterrows()]
        forest(ax, rows, f"{'C' if oc == 'RFS' else 'D'}  "
                         f"{'Relapse-free' if oc == 'RFS' else 'Overall'} survival, "
                         "multivariable Cox")
        ax.set_xlim(0.7, 1.9)
    fig.tight_layout()
    fig.savefig(FIG / "Fig4CD_forest_nerve_axes.png", dpi=DPI)
    fig.savefig(FIG / "Fig4CD_forest_nerve_axes.pdf")
    plt.close(fig)


def fig4e(S):
    fig, ax = plt.subplots(figsize=(5.2, 4.4))
    ax.scatter(S["CAF_ECM"], S["NPPscore"], s=6, c="#34495e", alpha=0.18, linewidths=0)
    b, a = np.polyfit(S["CAF_ECM"], S["NPPscore"], 1)
    xs = np.linspace(S["CAF_ECM"].min(), S["CAF_ECM"].max(), 50)
    ax.plot(xs, a + b * xs, color="#c0392b", lw=1.8)
    r = stats.spearmanr(S["CAF_ECM"], S["NPPscore"])
    ax.set_xlabel("Generic CAF/ECM score (z)")
    ax.set_ylabel("NPP score (z)")
    ax.set_title("E  Bulk NPP is largely stromal", loc="left", fontweight="bold")
    ax.text(0.04, 0.95, f"Spearman ρ = {r.statistic:.3f}", transform=ax.transAxes,
            va="top", fontsize=9)
    fig.tight_layout()
    fig.savefig(FIG / "Fig4E_npp_vs_caf.png", dpi=DPI)
    fig.savefig(FIG / "Fig4E_npp_vs_caf.pdf")
    plt.close(fig)


def fig4f(d):
    fig, ax = plt.subplots(figsize=(5.4, 4.3))
    _km_ax(ax, d, "RFS_MONTHS", "RFS_EVENT", "NPPresid_tertile",
           "F  RFS by stroma-adjusted NPP tertile")
    lo = d[d["NPPresid_tertile"] == "T1_low"]; hi = d[d["NPPresid_tertile"] == "T3_high"]
    chi, p = nsurv.logrank(lo["RFS_MONTHS"], lo["RFS_EVENT"],
                           hi["RFS_MONTHS"], hi["RFS_EVENT"])
    ax.text(0.97, 0.30, f"log-rank P = {p:.3g}", transform=ax.transAxes,
            ha="right", fontsize=8)
    fig.tight_layout()
    fig.savefig(FIG / "Fig4F_km_npp_resid.png", dpi=DPI)
    fig.savefig(FIG / "Fig4F_km_npp_resid.pdf")
    plt.close(fig)


def fig4g(perm):
    if perm is None or len(perm) == 0:
        return
    r = perm.iloc[0]
    # 零分布无法从汇总统计重建, 这里画观测值与零区间的位置关系
    fig, ax = plt.subplots(figsize=(5.6, 2.8))
    xs = np.linspace(r["null_2_5"], r["null_97_5"], 200)
    sig = (r["null_97_5"] - r["null_2_5"]) / 4.0
    ax.plot(xs, np.exp(-0.5 * ((xs - r["null_median_HR"]) / max(sig, 1e-6)) ** 2),
            color="#7f8c8d", lw=1.6)
    ax.fill_between(xs, 0, np.exp(-0.5 * ((xs - r["null_median_HR"]) / max(sig, 1e-6)) ** 2),
                    color="#bdc3c7", alpha=0.5)
    ax.axvline(r["null_median_HR"], color="#7f8c8d", lw=1.0, ls="--",
               label=f"null median HR={r['null_median_HR']}")
    ax.axvline(r["observed_HR"], color="#c0392b", lw=2.0,
               label=f"observed HR={r['observed_HR']}")
    ax.set_yticks([])
    ax.set_xlabel("Hazard ratio (per SD)")
    ax.set_title(f"G  Permutation null ({int(r['n_perm'])} random gene sets), "
                 f"empirical P={r['empirical_p']}", loc="left",
                 fontweight="bold", fontsize=9)
    ax.legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(FIG / "Fig4G_permutation_null.png", dpi=DPI)
    fig.savefig(FIG / "Fig4G_permutation_null.pdf")
    plt.close(fig)


def fig4h(tab):
    sub = tab[tab["model"].str.startswith("subtype_") & (tab["outcome"] == "RFS")]
    sub = sub[sub["name"].str.contains("NPPscore|NPP_resid")]
    if len(sub) == 0:
        return
    fig, ax = plt.subplots(figsize=(6.0, 3.6))
    labs = [f"{r['name'].split('[')[1][:-1]}\n{r['name'].split('[')[0]}"
            for _, r in sub.iterrows()]
    x = np.arange(len(sub))
    cols = ["#c0392b" if "NPPscore" in n else "#2980b9" for n in sub["name"]]
    ax.errorbar(x, sub["hr"],
                yerr=[sub["hr"] - sub["lo"], sub["hi"] - sub["hr"]],
                fmt="none", ecolor=cols, lw=1.2)
    ax.scatter(x, sub["hr"], c=cols, s=28, zorder=3)
    ax.axhline(1, color="#888", lw=0.8, ls="--")
    ax.set_xticks(x); ax.set_xticklabels(labs, fontsize=7)
    ax.set_ylabel("HR per SD (within subtype)")
    ax.set_title("H  NPP effect within molecular subtypes (RFS)",
                 loc="left", fontweight="bold", fontsize=10)
    fig.tight_layout()
    fig.savefig(FIG / "Fig4H_subtype_forest.png", dpi=DPI)
    fig.savefig(FIG / "Fig4H_subtype_forest.pdf")
    plt.close(fig)


def main():
    d = pd.read_csv(PROC / "metabric_analysis_table.csv", index_col=0)
    S = pd.read_csv(PROC / "metabric_axis_scores.csv", index_col=0)
    tab = pd.read_csv(TAB / "step4b_metabric_nerve_axes.csv")
    perm_path = TAB / "step4_metabric_permutation.csv"
    perm = pd.read_csv(perm_path) if perm_path.exists() else None

    d["NPP_tertile"] = pd.qcut(d["NPPscore"], 3, labels=["T1_low", "T2", "T3_high"])
    d["NPPresid_tertile"] = pd.qcut(d["NPP_resid_stroma"], 3,
                                    labels=["T1_low", "T2", "T3_high"])
    fig4ab(d); print("  Fig4A/B")
    fig4cd(tab); print("  Fig4C/D")
    fig4e(S); print("  Fig4E")
    fig4f(d); print("  Fig4F")
    fig4g(perm); print("  Fig4G")
    fig4h(tab); print("  Fig4H")
    print("[Fig4] done")


if __name__ == "__main__":
    main()

# -*- coding: utf-8 -*-
"""Figure 5 —— 外部独立队列的器官特异性复发验证 (GSE12276, n=204)。

panels
------
5A NPP 三分位的脑转移 Aalen-Johansen 累积发生率曲线 (竞争事件=其他部位远处复发)
5B 全部神经轴/NPP 的 Fine-Gray 亚分布风险森林图 (cause1=脑转移)
5C 终点定义敏感性: 宽松(含 brain+other) vs 严格(仅 brain) 的 sHR 对比
5D 各复发部位的 NPP 与 CD8 评分分布
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
DPI = 300
COL = {"T1_low": "#2980b9", "T2": "#7f8c8d", "T3_high": "#c0392b"}

plt.rcParams.update({
    "font.size": 8, "axes.titlesize": 9, "axes.labelsize": 8,
    "xtick.labelsize": 7, "ytick.labelsize": 7, "legend.fontsize": 7,
    "svg.fonttype": "none", "font.family": "DejaVu Sans",
})


def fig5a(cif):
    fig, ax = plt.subplots(figsize=(4.6, 3.6))
    for lev, g in cif.groupby("tertile"):
        ax.step(g["time"], g["cif_brain"], where="post", lw=1.6,
                color=COL.get(lev, "#555"),
                label=f"{lev} (NPP)")
    ax.set_xlabel("Time from surgery to relapse (months)")
    ax.set_ylabel("Cumulative incidence of brain relapse")
    ax.set_ylim(0, None)
    ax.legend(frameon=True, title="", loc="upper left")
    ax.grid(alpha=0.2, lw=0.4)
    ax.set_title("A  Brain-relapse CIF by NPP tertile (GSE12276, n=204)",
                 loc="left", fontweight="bold")
    fig.tight_layout()
    fig.savefig(FIG / "Fig5A_cif_brain_npp.png", dpi=DPI, bbox_inches="tight")
    fig.savefig(FIG / "Fig5A_cif_brain_npp.pdf", bbox_inches="tight")
    plt.close(fig)
    print("  Fig5A")


def fig5b(res):
    r = res[res["cause"] == "brain_first_finegray"].dropna(subset=["hr"])
    r = r[~r["variable"].str.startswith("NPP_tertile")]
    order = ["NPPscore"] + [v for v in r["variable"] if v != "NPPscore"]
    r = r.set_index("variable").loc[order].reset_index()
    fig, ax = plt.subplots(figsize=(5.2, 4.4))
    y = np.arange(len(r))[::-1]
    ax.errorbar(r["hr"], y,
                xerr=np.c_[r["hr"] - r["lo"], r["hi"] - r["hr"]].T,
                fmt="o", ms=4.5, capsize=2.5, lw=1.1, color="#2C3E50",
                ecolor="#7F8C8D")
    ax.axvline(1.0, ls="--", lw=0.9, color="#C0392B")
    ax.set_yticks(y)
    ax.set_yticklabels(r["variable"], fontsize=7)
    ax.set_xscale("log")
    import matplotlib.ticker as mticker
    ax.xaxis.set_major_locator(mticker.FixedLocator([0.25, 0.5, 1, 2, 3]))
    ax.xaxis.set_major_formatter(mticker.FixedFormatter(["0.25", "0.5", "1", "2", "3"]))
    ax.xaxis.set_minor_locator(mticker.NullLocator())
    ax.set_xlim(0.2, 3.2)
    ax.set_xlabel("Fine-Gray subdistribution HR per SD (brain relapse)")
    for yi, (_, row) in zip(y, r.iterrows()):
        star = " *" if row["p"] < 0.05 else ""
        ax.text(row["hi"] * 1.06, yi, f"p={row['p']:.2g}{star}",
                va="center", fontsize=6, color="#555")
    ax.grid(alpha=0.2, lw=0.4, axis="x")
    ax.set_title("B  Competing-risk models, brain-first relapse",
                 loc="left", fontweight="bold")
    fig.tight_layout()
    fig.savefig(FIG / "Fig5B_forest_brainmet.png", dpi=DPI, bbox_inches="tight")
    fig.savefig(FIG / "Fig5B_forest_brainmet.pdf", bbox_inches="tight")
    plt.close(fig)
    print("  Fig5B")


def fig5c(res):
    r = res[res["variable"].isin(["NPPscore", "CD8_cytotox"]) &
            res["cause"].isin(["brain_first_finegray", "brain_only_strict_finegray"])]
    fig, ax = plt.subplots(figsize=(4.6, 3.2))
    lab = {"brain_first_finegray": "brain or brain+other (n=16)",
           "brain_only_strict_finegray": "brain only, strict (n=8)"}
    x = np.arange(len(lab))
    w = 0.36
    for i, var in enumerate(["NPPscore", "CD8_cytotox"]):
        s = r[r["variable"] == var].set_index("cause")
        hr = [s.loc[c, "hr"] for c in lab]
        lo = [s.loc[c, "lo"] for c in lab]
        hi = [s.loc[c, "hi"] for c in lab]
        ax.errorbar(x + (i - 0.5) * w, hr,
                    yerr=[np.array(hr) - np.array(lo), np.array(hi) - np.array(hr)],
                    fmt="o", ms=5, capsize=3, lw=1.2,
                    color=["#2980B9", "#C0392B"][i], label=var)
    ax.set_xticks(x)
    ax.set_xticklabels([f"{v}\n{s}" for v, s in
                        [("broad", "16 events"), ("strict", "8 events")]], fontsize=7)
    ax.set_xticklabels([f"broad definition\n(16 brain events)",
                        "strict brain-only\n(8 events)"], fontsize=7)
    ax.axhline(1.0, ls="--", lw=0.9, color="#333")
    ax.set_ylabel("subdistribution HR per SD")
    ax.legend(frameon=True, fontsize=6.5)
    ax.grid(alpha=0.2, lw=0.4, axis="y")
    ax.set_title("C  Endpoint-definition sensitivity",
                 loc="left", fontweight="bold")
    fig.tight_layout()
    fig.savefig(FIG / "Fig5C_sensitivity_endpoint.png", dpi=DPI, bbox_inches="tight")
    fig.savefig(FIG / "Fig5C_sensitivity_endpoint.pdf", bbox_inches="tight")
    plt.close(fig)
    print("  Fig5C")


def fig5d(d):
    site_order = ["brain", "brain+other", "other", "local"]
    site_lab = ["brain", "brain+other", "distant other", "local only"]
    fig, axes = plt.subplots(1, 2, figsize=(9.2, 3.4), sharey=False)
    for ax, var, ttl in zip(axes, ["NPPscore", "CD8_cytotox"],
                            ["NPP score", "CD8 cytotoxic score"]):
        data = [d.loc[d["site"] == s, var].dropna().values for s in site_order]
        bp = ax.boxplot(data, tick_labels=site_lab, patch_artist=True, widths=0.55)
        for b in bp["boxes"]:
            b.set_facecolor("#E8F6F3"); b.set_edgecolor("#148F77")
        for m in bp["medians"]:
            m.set_color("#0E6655"); m.set_linewidth(1.2)
        rng = np.random.default_rng(0)
        for i, dd in enumerate(data):
            ax.scatter(rng.normal(i + 1, 0.05, len(dd)), dd, s=5,
                       color="#117A65", alpha=0.55, zorder=3)
        ax.set_title(f"{ttl}", loc="left", fontweight="bold")
        ax.grid(alpha=0.2, lw=0.4, axis="y")
        ax.tick_params(axis="x", labelsize=6.5)
    axes[0].set_ylabel("z score")
    fig.tight_layout()
    fig.savefig(FIG / "Fig5D_scores_by_site.png", dpi=DPI, bbox_inches="tight")
    fig.savefig(FIG / "Fig5D_scores_by_site.pdf", bbox_inches="tight")
    plt.close(fig)
    print("  Fig5D")


def main():
    res = pd.read_csv(TAB / "step5_gse12276_brainmet.csv")
    cif = pd.read_csv(TAB / "step5_gse12276_cif_tertiles.csv")
    d = pd.read_csv(CFG.DIR_PROC / "gse12276_analysis_table.csv", index_col=0)
    fig5a(cif)
    fig5b(res)
    fig5c(res)
    fig5d(d)
    print("[Fig5] done")


if __name__ == "__main__":
    main()

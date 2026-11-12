# -*- coding: utf-8 -*-
"""Figure 7 —— 空间演化轴(扩散伪时间) 与 NPP 双模块拆分。

面板
----
A  每样本 ρ(伪时间, 神经距离): 伪时间只由表达构建, 与距离无关,
   却在 23/30 例中正相关 —— 空间演化轴的非平凡一致性检验
B  读数沿伪时间分箱的轨迹 (跨样本均值, 每条轨迹跨箱 z 标准化)
C  NPP 双模块的每样本空间梯度: 纯神经 vs 神经诱导基质 (配对散点)
D  bulk (METABRIC) 层面: 两模块与 CAF/ECM 及原 NPP 的相关
E  METABRIC 生存: 两模块的 HR (单变量 / 校正 / 互相校正 / +CAF)
F  bootstrap HR 分布 (250 次, RFS)

风格与 Fig1-6 一致: 8pt, 300 dpi, PNG + PDF。
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


def panel_a(P, S, ax):
    v = P["rho_dpt_dist"].to_numpy(float)
    x = np.arange(1, len(v) + 1)
    o = np.argsort(v)
    ax.axhline(0, color="black", lw=0.5)
    ax.bar(x, v[o], color=[C_MAIN if t > 0 else C_ALT for t in v[o]],
           alpha=0.85, width=0.72)
    ax.set_xticks(x)
    ax.set_xticklabels([str(i + 1) for i in range(len(v))], fontsize=5)
    ax.set_xlabel("Visium sample (ordered by ρ)")
    ax.set_ylabel("Spearman ρ (pseudotime, nerve distance)")
    ax.set_title("A  Expression-derived pseudotime vs physical distance",
                 loc="left")
    txt = (f"{int((v > 0).sum())}/{len(v)} samples positive\n"
           f"median ρ = {np.median(v):+.3f}\n"
           f"cross-sample t = {S['t_across_samples']:.2f}, "
           f"p = {S['p_across_samples']:.4f}\n"
           f"partial (UMI & pct_mt adjusted):\n"
           f"t = {S['t_partial_across_samples']:.2f}, "
           f"p = {S['p_partial_across_samples']:.4f}\n"
           f"root spot in 0-50 µm ring: {S['frac_root_in_r1']:.0%}")
    ax.text(0.98, 0.97, txt, transform=ax.transAxes, ha="right", va="top",
            fontsize=6,
            bbox=dict(boxstyle="round,pad=0.3", fc="#F5F5F5", ec="#CCCCCC", lw=0.5))


def panel_b(Pf, ax):
    sel = ["Neuronal_Schwann", "Pericyte_SMC", "CAF_myofibroblastic",
           "CAF_inflammatory_dev", "Endothelial", "Tumor_proliferative",
           "CD8_T", "T_lymphocyte", "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
           "HALLMARK_TGF_BETA_SIGNALING", "REACTOME_EXTRACELLULAR_MATRIX_ORGANIZATION",
           "HALLMARK_HYPOXIA"]
    cols = [C_PURPLE, C_MAIN, C_ALT, "#D98C3F", C_ACC, "#5B7A8C",
            "#3B7A5B", "#9A3B6B", "#8C6B4E", "#4E6B8C", "#6B8C4E", "#8C4E6B"]
    num = Pf.drop(columns=["sample", "bin", "n_spot", "dpt_median", "dist_median"])
    m = num.groupby(Pf["bin"]).mean()
    xs = m.index.to_numpy(float)
    for c, k in zip(cols, sel):
        if k not in m.columns:
            continue
        v = m[k].to_numpy(float)
        v = (v - v.mean()) / (v.std() + 1e-12)
        ax.plot(xs, v, "-o", ms=2.8, lw=1.1, color=c,
                label=k.replace("_", " ").replace("HALLMARK ", "")
                        .replace("REACTOME ", ""))
    ax.axhline(0, color="black", lw=0.5)
    ax.set_xlabel("Diffusion-pseudotime decile (0 = nerve-end root, 9 = distal)")
    ax.set_ylabel("Readout trajectory (z across bins)")
    ax.set_title("B  Readout ordering along the perineural trajectory", loc="left")
    ax.legend(frameon=False, fontsize=5.4, ncol=2, columnspacing=0.7,
              handlelength=1.1, loc="upper right")


def panel_c(ax):
    d = pd.read_csv(PROC / "npp_split_spatial.csv")
    g = d.pivot_table(index="sample", columns="module", values="rho")
    n = g["Neural_intrinsic"].to_numpy(float)
    s = g["Perineural_stroma"].to_numpy(float)
    ok = np.isfinite(n) & np.isfinite(s)
    n, s = n[ok], s[ok]
    ax.axhline(0, color="black", lw=0.4)
    ax.axvline(0, color="black", lw=0.4)
    ax.scatter(n, s, s=22, color=C_MAIN, alpha=0.8, edgecolors="white", lw=0.4)
    lim = [min(n.min(), s.min()) - 0.04, max(n.max(), s.max()) + 0.04]
    ax.plot(lim, lim, "--", color=C_GREY, lw=0.8, label="y = x")
    r, p = stats.spearmanr(n, s)
    ax.set_xlim(lim)
    ax.set_ylim(lim)
    ax.set_xlabel("Neural-intrinsic module ρ (per sample)")
    ax.set_ylabel("Perineural-stroma module ρ (per sample)")
    ax.set_title("C  Two modules, one direction: paired spatial gradients",
                 loc="left")
    ax.legend(frameon=False, fontsize=6, loc="upper left")
    n_sig = int(np.sum(n < 0))
    s_sig = int(np.sum(s < 0))
    ax.text(0.98, 0.04, f"Spearman ρ across samples = {r:+.2f}\n"
            f"negative in {n_sig}/30 (neural) vs {s_sig}/30 (stroma) samples",
            transform=ax.transAxes, ha="right", fontsize=6.5,
            bbox=dict(boxstyle="round,pad=0.3", fc="#F5F5F5", ec="#CCCCCC", lw=0.5))


def panel_d(B, ax):
    keep = [
        ("Neural vs Perineural-stroma", "Neural_intrinsic", "Perineural_stroma"),
        ("Neural vs CAF/ECM", "Neural_intrinsic", "CAF_ECM"),
        ("Stroma vs CAF/ECM", "Perineural_stroma", "CAF_ECM"),
        ("Neural vs original NPP", "Neural_intrinsic", "NPPscore"),
        ("Stroma vs original NPP", "Perineural_stroma", "NPPscore"),
    ]
    y, v, lab = [], [], []
    for i, (txt, a, b) in enumerate(keep):
        row = B[(B["var1"] == a) & (B["var2"] == b)]
        if not len(row):
            row = B[(B["var1"] == b) & (B["var2"] == a)]
        if not len(row):
            continue
        y.append(i)
        v.append(float(row["rho"].iloc[0]))
        lab.append(txt)
    cols = [C_GREY if abs(x) < 0.2 else C_MAIN for x in v]
    cols = [C_ALT if abs(x) > 0.9 else c for c, x in zip(cols, v)]
    ax.barh(y, v, color=cols, alpha=0.9, height=0.6)
    for i, x in zip(y, v):
        ax.text(x + (0.01 if x > 0 else -0.01), i, f"{x:+.3f}", va="center",
                ha="left" if x > 0 else "right", fontsize=6.5)
    ax.axvline(0, color="black", lw=0.5)
    ax.set_yticks(y)
    ax.set_yticklabels(lab, fontsize=7)
    ax.set_xlim(-0.1, 1.05)
    ax.set_xlabel("Spearman ρ in METABRIC (n = 1,976)")
    ax.set_title("D  Neural module orthogonal to, stroma module collinear with CAF/ECM",
                 loc="left", fontsize=8)


def panel_e(D, ax):
    keep = [("RFS", "RFS_N_unadj"), ("RFS", "RFS_S_unadj"),
            ("RFS", "RFS_N_adj"), ("RFS", "RFS_S_adj"),
            ("RFS", "RFS_mutual"), ("RFS", "RFS_mutual_CAF"),
            ("OS", "OS_N_unadj"), ("OS", "OS_S_unadj"),
            ("OS", "OS_N_adj"), ("OS", "OS_S_adj"),
            ("OS", "OS_mutual"), ("OS", "OS_mutual_CAF")]
    rows = []
    for oc, md in keep:
        r = D[(D["outcome"] == oc) & (D["model"] == md)]
        if len(r):
            rows.append((oc, md, float(r["hr"].iloc[0]),
                         float(r["lo"].iloc[0]), float(r["hi"].iloc[0]),
                         float(r["p"].iloc[0]),
                         r["variable"].iloc[0]))
    y = np.arange(len(rows))[::-1]
    ax.axvline(1, color="black", lw=0.6)
    for yy, (oc, md, hr, lo, hi, p, var) in zip(y, rows):
        c = C_MAIN if var == "Neural_intrinsic" else C_ALT
        ax.plot([lo, hi], [yy, yy], color=c, lw=1.4, solid_capstyle="butt")
        ax.plot(hr, yy, "o", ms=4, color=c)
    ax.set_yticks(y)
    ax.set_yticklabels([f"{oc} · {md.split('_', 1)[1]}" for oc, md, *_ in rows],
                       fontsize=6.2)
    ax.set_xlabel("Hazard ratio per SD increase (METABRIC)")
    ax.set_title("E  Prognostic models: neither module is a stable risk marker",
                 loc="left")
    ax.text(0.02, 0.02, "mutual_CAF = both modules + CAF/ECM\n"
            "(exploratory, not bootstrap-validated)", transform=ax.transAxes,
            ha="left", fontsize=5.6, color="#555555")


def panel_f(Z, ax):
    draws = Z["RFS_Neural_intrinsic"]
    draws_s = Z["RFS_Perineural_stroma"]
    parts = ax.violinplot([draws, draws_s], positions=[1, 2], widths=0.7,
                          showmedians=True)
    for pc, c in zip(parts["bodies"], [C_MAIN, C_ALT]):
        pc.set_facecolor(c)
        pc.set_alpha(0.55)
        pc.set_edgecolor("none")
    for k in ("cmedians", "cmins", "cmaxes", "cbars"):
        parts[k].set_color("black")
        parts[k].set_linewidth(0.8)
    ax.axhline(1, color="black", lw=0.6, ls="--")
    ax.set_xticks([1, 2])
    ax.set_xticklabels(["Neural-intrinsic", "Perineural-stroma"], fontsize=7)
    ax.set_ylabel("Bootstrap HR (RFS, 250 resamples)")
    ax.set_title("F  Bootstrap stability of the null prognostic result",
                 loc="left")
    for i, d in enumerate([draws, draws_s], 1):
        lo, hi = np.percentile(d, [2.5, 97.5])
        ax.text(i, hi + 0.004, f"95% CI {lo:.2f}-{hi:.2f}", ha="center",
                fontsize=6)


def main():
    P = pd.read_csv(PROC / "pseudotime_per_sample.csv.gz")
    Pf = pd.read_csv(PROC / "pseudotime_profiles.csv")
    S = pd.read_csv(TAB / "pseudotime_summary.csv").iloc[0].to_dict()
    Ps = pd.read_csv(TAB / "npp_split_spatial_meta.csv")
    B = pd.read_csv(TAB / "npp_split_bulk_correlation.csv")
    D = pd.read_csv(TAB / "npp_split_bulk_cox.csv")
    Z = np.load(TAB / "npp_split_bootstrap.npz", allow_pickle=True)
    print(f"[Fig7] pseudotime {P.shape}, profiles {Pf.shape}")

    fig, axes = plt.subplots(2, 3, figsize=(13.6, 8.2))
    panel_a(P, S, axes[0, 0])
    panel_b(Pf, axes[0, 1])
    panel_c(axes[0, 2])
    panel_d(B, axes[1, 0])
    panel_e(D, axes[1, 1])
    panel_f(Z, axes[1, 2])
    fig.tight_layout()
    save(fig, "Fig7_trajectory_and_modules")


if __name__ == "__main__":
    main()

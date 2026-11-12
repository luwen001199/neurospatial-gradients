# -*- coding: utf-8 -*-
"""Figure 1 —— 空间神经梯度建模 (Ng et al. style multi-panel).

panels:
  Fig1A  方法流程示意 (自己画, 不依赖数据)
  Fig1B  代表样本: Nervesignal 空间图 + 神经源点标注
  Fig1C  代表样本: 距离场 + 环形缓冲区
  Fig1D  代表样本: Collagen / myCAF / CD8 三张空间图
  Fig1E  代表样本: Nervesignal 拟合曲线(散点+分箱+拟合线+1/b)
  Fig1F  多样本 Nervesignal 环带衰减曲线叠加
  Fig1G  8 个指标的 b / 1/b / R2 分布
  Fig1H  神经识别有效性: nerve vs non-nerve 的 Nervesignal (效应量)
"""
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Circle, FancyArrowPatch
from scipy import stats

ROOT = Path(r"D:/嗜神经浸润")
sys.path.insert(0, str(ROOT / "00_pipeline"))

from nsgrad import config as CFG
from nsgrad import viz as nviz
from nsgrad import io as nio, signature as nsig, nerve as nnerve, gradient as ngrad

PROC = CFG.DIR_PROC
FIG = CFG.DIR_FIG
FIG.mkdir(parents=True, exist_ok=True)

METRICS = list(CFG.GRADIENT_SIGNATURES.keys())


def ready_samples():
    out = []
    for d in sorted((CFG.DIR_RAW / "gse210616").glob("GSM*")):
        if not d.is_dir() or list(d.glob("*.part")):
            continue
        if list(d.glob("*_filtered_feature_bc_matrix.h5")) and list(d.glob("*tissue_positions*")):
            out.append(d)
    return out


def load_one(d):
    ad = nio.read_visium_sample(d, sample_id=d.name)
    S, _ = nsig.score_signatures(ad, CFG.GRADIENT_SIGNATURES)
    nv = nnerve.identify_nerve_spots(ad)
    mask = (nv["nerve_spatial"] == 1).values
    distf, _ = ngrad.build_distance_fields(ad, mask)
    for c in S.columns:
        ad.obs[c] = S[c].values
    ad.obs["distance_to_nearest_nerve"] = distf["distance_to_nearest_nerve"].values
    return ad, S, nv, mask, distf


# ------------------------------------------------------------------ Fig1A
def fig1a():
    fig, ax = plt.subplots(figsize=(7.2, 2.0))
    ax.set_xlim(0, 100); ax.set_ylim(0, 26); ax.axis("off")

    def box(x, y, w, h, t, sub="", fc="#EBF5FB", ec="#2874A6"):
        ax.add_patch(plt.Rectangle((x, y), w, h, facecolor=fc, edgecolor=ec,
                                   lw=1.0, zorder=2))
        ax.text(x + w / 2, y + h / 2 + (0.9 if sub else 0), t, ha="center",
                va="center", fontsize=7, zorder=3, weight="bold")
        if sub:
            ax.text(x + w / 2, y + h / 2 - 1.1, sub, ha="center", va="center",
                    fontsize=5.2, zorder=3, color="#555")

    def arrow(x1, x2, y):
        ax.add_patch(FancyArrowPatch((x1, y), (x2, y), mutation_scale=9,
                                     color="#555", lw=1.0,
                                     arrowstyle="-|>", zorder=1))

    box(1, 13, 17, 8, "Visium ST", "43 TNBC samples")
    arrow(18, 20, 17)
    box(20, 13, 19, 8, "Gene signature", "z-score per gene,\nmean per spot")
    arrow(39, 42, 17)
    box(42, 13, 20, 8, "Nerve source", ">=2 co-expressed\nmarkers + spatial")
    arrow(62, 65, 17)
    box(65, 13, 17, 8, "Ring buffers", "0-50 / 50-100 /\n100-200 / 200-400 um")
    arrow(82, 84, 17)
    box(84, 13, 15, 8, "Exp decay fit", "y = a*e^(-bd) + c\n1/b = radius")

    # 下方: K=2 二分型 (探索性)
    box(26, 2.5, 22, 6.5, "Neural-Active", "high nerve fraction,\nsteep perineural decay", fc="#FDEDEC", ec="#C0392B")
    box(54, 2.5, 22, 6.5, "Neural-Cold", "low nerve fraction,\nweak gradient", fc="#EBF5FB", ec="#2E86C1")
    arrow(48, 41.5, 5.5); arrow(52, 63.5, 5.5)
    ax.text(50, 11.6, "consensus clustering, K = 2 (exploratory)", ha="center", fontsize=5.5, color="#666")

    # 衰减示意
    xs = np.linspace(0, 15, 100)
    axins = fig.add_axes([0.015, 0.10, 0.135, 0.36])
    for lbl, b, c in [("steep", 0.35, "#C0392B"), ("shallow", 0.08, "#2E86C1")]:
        axins.plot(xs, 1.8 * np.exp(-b * xs) + 0.15, color=c, lw=1.2, label=lbl)
    axins.axvline(1 / 0.35, color="#C0392B", ls=":", lw=0.8)
    axins.axvline(1 / 0.08, color="#2E86C1", ls=":", lw=0.8)
    axins.set_xticks([]); axins.set_yticks([])
    axins.set_xlabel("distance from nerve", fontsize=5)
    axins.set_ylabel("signal", fontsize=5)
    axins.legend(fontsize=4.6, frameon=False)
    for s in ("top", "right"):
        axins.spines[s].set_visible(False)
    return nviz.save(fig, "Fig1A_pipeline_schematic")


# ------------------------------------------------------------------ Fig1B/C/D/E
def fig1bcd(demo_ad, demo_S, demo_mask, demo_dist, demo_name):
    ad, S, mask, dist = demo_ad, demo_S, demo_mask, demo_dist
    # B: Nervesignal 空间 + 神经源点
    nviz.plot_spatial(ad, "Nervesignal",
                      title=f"Nerve signal (z) — {demo_name}",
                      nerve_mask=mask, size=16, name="Fig1B_nervesignal_map")
    # C: 距离场 + 环带
    nviz.plot_spatial(ad, "distance_to_nearest_nerve",
                      title="Distance to nearest nerve (um)",
                      cmap="plasma", cmap_div=False, nerve_mask=mask,
                      size=16, name="Fig1C_distance_field")
    # D: 三联图
    fig, axes = plt.subplots(1, 3, figsize=(9.6, 3.1))
    for ax, m in zip(axes, ["CollagenIndex", "myCAFscore", "CD8density"]):
        nviz.plot_spatial(ad, m, title=f"{m} (z)", nerve_mask=mask,
                          size=13, ax=ax, name=None)
    for a in axes:
        a.set_aspect("equal")
    nviz.save(fig, "Fig1D_collagen_myCAF_CD8_maps")
    # E: 拟合曲线
    for m in ["Nervesignal", "CollagenIndex", "CD8density"]:
        r = ngrad.fit_decay_aggregated(dist, S[m].values)
        nviz.plot_fit_single(dist, S[m].values, r, m, sample=demo_name,
                             name=f"Fig1E_fit_{m}")


# ------------------------------------------------------------------ Fig1F
def fig1f(all_profiles):
    p = pd.concat(all_profiles, ignore_index=True)
    for m in ["Nervesignal", "CollagenIndex", "CD8density"]:
        nviz.plot_decay_curves(p, m, name=f"Fig1F_decaycurves_{m}")


# ------------------------------------------------------------------ Fig1G
def fig1g(gp):
    import re
    piv_b = gp.pivot_table(index="sample", columns="metric", values="b")
    piv_r2 = gp.pivot_table(index="sample", columns="metric", values="r2")
    piv_dp = gp.pivot_table(index="sample", columns="metric", values="decay_points_um")
    order = ["Nervesignal"] + [c for c in METRICS if c != "Nervesignal"]
    piv_b = piv_b[order]; piv_r2 = piv_r2[order]; piv_dp = piv_dp[order]

    nviz.plot_heatmap(piv_r2, cmap="RdYlBu_r", vmin=0, vmax=1,
                      title="Fit quality R2 (sample x metric)",
                      xlabel="metric", ylabel="sample",
                      name="Fig1G1_r2_heatmap")

    fig, axes = plt.subplots(1, 3, figsize=(10.5, 3.0))
    for i, m in enumerate(["Nervesignal", "CollagenIndex", "CD8density"]):
        v = piv_b[m].dropna()
        v = v[(v > 1 / 2000 * 1.001) & (v < 1 / 25 * 0.999)]
        v = 1 / v
        v = v[v <= np.percentile(v, 98)] if len(v) > 3 else v
        if len(v) == 0:
            axes[i].text(0.5, 0.5, "no interpretable fits", ha="center",
                         transform=axes[i].transAxes, fontsize=7)
            axes[i].set_title(m); continue
        axes[i].hist(v, bins=min(14, max(4, len(v) // 2)), color="#7FB3D5",
                     edgecolor="white", lw=0.5)
        axes[i].axvline(v.median(), color="#C0392B", ls="--", lw=1.2,
                        label=f"median={v.median():.0f}um")
        axes[i].set_xlabel("effective radius 1/b (um)")
        axes[i].set_ylabel("n samples")
        axes[i].set_title(m, fontsize=9)
        axes[i].legend(fontsize=6, frameon=False)
        for s in ("top", "right"):
            axes[i].spines[s].set_visible(False)
    nviz.save(fig, "Fig1G2_decay_radius_distribution")

    # 可辨识性条形图
    rows = []
    for m in METRICS:
        sub = gp[gp["metric"] == m]
        b = sub["b"].astype(float)
        lo = int((np.abs(b - 1 / 2000) < 1e-6).sum())
        hi = int((np.abs(b - 1 / 25) < 1e-6).sum())
        rows.append({"metric": m, "interpretable": len(sub) - lo - hi,
                     "at_low_bound": lo, "at_high_bound": hi})
    D = pd.DataFrame(rows).set_index("metric").loc[order]
    fig, ax = plt.subplots(figsize=(4.2, 2.8))
    D.plot(kind="barh", stacked=True, ax=ax,
           color=["#27AE60", "#E74C3C", "#F39C12"], width=0.72)
    ax.set_xlabel("n samples"); ax.set_ylabel("")
    ax.set_title("Parameter identifiability", fontsize=9)
    ax.legend(["interpretable", "at lower bound", "at upper bound"],
              fontsize=5.5, frameon=False, loc="lower right")
    ax.tick_params(axis="y", labelsize=6)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    nviz.save(fig, "Fig1G3_identifiability")


# ------------------------------------------------------------------ Fig1H
def fig1h(prev, nrv):
    """H: (left) nerve-spot abundance across the 30 QC-passing samples;
    (right) per-sample Cohen's d versus AUC for the composite nerve score in
    called versus non-called spots (internal-consistency, NOT validation)."""
    fig, ax = plt.subplots(1, 2, figsize=(6.4, 2.7))
    v = prev["pct_nerve_spots"].astype(float).dropna()
    ax[0].hist(v, bins=12, color="#85C1E9", edgecolor="white")
    ax[0].axvline(v.median(), color="#C0392B", ls="--", lw=1.2,
                  label=f"median={v.median():.1f}% (n={len(v)})")
    ax[0].set_xlabel("nerve-associated spots (%)")
    ax[0].set_ylabel("n samples")
    ax[0].set_title("Nerve-spot abundance", fontsize=9)
    ax[0].legend(fontsize=6, frameon=False)
    for s in ("top", "right"):
        ax[0].spines[s].set_visible(False)

    x = nrv["d"].astype(float); y = nrv["auc"].astype(float)
    m = np.isfinite(x) & np.isfinite(y)
    ax[1].scatter(x[m], y[m], s=26, c="#2874A6", edgecolors="white", linewidths=0.4)
    ax[1].axvline(0.8, color="grey", ls=":", lw=0.8)
    ax[1].axhline(0.5, color="grey", ls=":", lw=0.8)
    ax[1].set_xlabel("Cohen's d (nerve vs non-nerve)")
    ax[1].set_ylabel("AUC")
    ax[1].set_title("Internal coherence of nerve calls", fontsize=9)
    for xx, yy, s in zip(x[m], y[m], nrv["sample"][m]):
        ax[1].annotate(str(s).replace("GSM", ""), (xx, yy), fontsize=3.6, alpha=0.6)
    for s in ("top", "right"):
        ax[1].spines[s].set_visible(False)
    nviz.save(fig, "Fig1H_nerve_call_validation")


def main():
    print("[Fig1] loading samples ...")
    samples = ready_samples()
    print(f"      {len(samples)} ready")
    gp = pd.read_csv(PROC / "gradient_params.csv")
    nrv = pd.read_csv(PROC / "nerve_validation.csv")
    prev = pd.read_csv(PROC / "sample_summary.csv")

    all_profiles = []
    demo = None
    for i, d in enumerate(samples):
        try:
            ad, S, nv, mask, distf = load_one(d)
        except Exception as e:
            print(f"      skip {d.name}: {type(e).__name__}")
            continue
        score = abs(mask.sum() - 150)   # 挑一个神经源点数量居中的样本做代表图
        if demo is None or score < demo[0]:
            demo = (score, ad, S, mask, distf, d.name)
        prof = ngrad.ring_profile(ad, {m: S[m].values for m in METRICS},
                                  distf, nerve_mask=mask)
        prof.insert(0, "sample", d.name)
        all_profiles.append(prof)

    fig1a()
    print("      Fig1A done")
    if demo:
        _, ad, S, mask, distf, name = demo
        print(f"      demo sample = {name}  nerve={mask.sum()}")
        fig1bcd(ad, S, mask, distf["distance_to_nearest_nerve"].values, name)
        print("      Fig1B-E done")
    if all_profiles:
        fig1f(all_profiles)
        print("      Fig1F done")
    fig1g(gp)
    print("      Fig1G done")
    fig1h(prev, nrv)
    print("      Fig1H done")
    print(f"[Fig1] all panels -> {FIG}")


if __name__ == "__main__":
    main()

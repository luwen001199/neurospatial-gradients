# -*- coding: utf-8 -*-
"""Figure 2 —— 空间神经-免疫-基质架构分型.

panels:
  Fig2A  PCA 样本散点 (按亚型着色) + 置信椭圆
  Fig2B  共识矩阵热图 (最优 K)
  Fig2C  K=2..4 的 PAC / silhouette 曲线
  Fig2D  特征 heatmap (样本 x 9 特征, 按亚型排序)
  Fig2E  关键特征按亚型的箱线图 (4 个: nerve_density / CD8_diff / collagen_diff / IFN_diff)
  Fig2F  每个亚型的代表样本空间图 (Nervesignal + CD8density 并排)
  Fig2G  亚型标记基因的空间分布 (TLSscore, Hypoxia)
"""
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Ellipse

ROOT = Path(r"D:/嗜神经浸润")
sys.path.insert(0, str(ROOT / "00_pipeline"))
from nsgrad import config as CFG
from nsgrad import viz as nviz
from nsgrad import io as nio, signature as nsig, nerve as nnerve, gradient as ngrad

PROC = CFG.DIR_PROC
FIG = CFG.DIR_FIG

PAL = {"Neural-Active": "#2E86C1", "Neural-Cold": "#7F8C8D",
       "Neural-Barriered": "#C0392B", "Neural-Inflamed": "#27AE60"}


def pca_2d(X):
    Xc = X - X.mean(0)
    U, S, Vt = np.linalg.svd(Xc, full_matrices=False)
    P = Xc @ Vt.T
    evr = S ** 2 / (S ** 2).sum()
    return P[:, :2], evr[:2], Vt.T[:, :2]


def fig2a(Fz, sub):
    P, evr, loadings = pca_2d(Fz.values)
    fig, ax = plt.subplots(figsize=(4.4, 3.6))
    for st, g in pd.DataFrame({"PC1": P[:, 0], "PC2": P[:, 1],
                               "st": sub.values}).groupby("st"):
        c = PAL.get(st, "#888")
        ax.scatter(g["PC1"], g["PC2"], c=c, s=52, alpha=0.9, label=st,
                   edgecolors="white", linewidths=0.7, zorder=3)
        # 95% 置信椭圆
        if len(g) >= 3:
            pts = g[["PC1", "PC2"]].values
            mu = pts.mean(0); cov = np.cov(pts.T)
            try:
                ev, evec = np.linalg.eigh(cov)
                order = np.argsort(ev)[::-1]
                ev, evec = ev[order], evec[:, order]
                ang = np.degrees(np.arctan2(evec[1, 0], evec[0, 0]))
                chi2 = 5.991
                w = np.sqrt(ev[0] * chi2); h = np.sqrt(ev[1] * chi2)
                ax.add_patch(Ellipse(mu, 2 * w, 2 * h, angle=ang, fill=False,
                                     edgecolor=c, lw=0.9, ls="--", alpha=0.7, zorder=2))
            except Exception:
                pass
    for i, s in enumerate(Fz.index):
        ax.annotate(str(s).replace("GSM", ""), (P[i, 0], P[i, 1]),
                    fontsize=4.6, alpha=0.7, xytext=(3, 3),
                    textcoords="offset points")
    ax.set_xlabel(f"PC1 ({evr[0]*100:.1f}%)")
    ax.set_ylabel(f"PC2 ({evr[1]*100:.1f}%)")
    ax.set_title("PCA of nerve-gradient features", fontsize=9)
    ax.legend(fontsize=6, frameon=False, loc="best")
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    nviz.save(fig, "Fig2A_pca_subtypes")


def fig2b(P, sub):
    fig, ax = plt.subplots(figsize=(4.0, 3.6))
    im = ax.imshow(P, cmap="Blues", vmin=0, vmax=1)
    cb = plt.colorbar(im, ax=ax, fraction=0.045)
    cb.ax.tick_params(labelsize=6); cb.outline.set_linewidth(0.5)
    order = np.argsort(pd.factorize(sub.values)[0])
    P2 = P[np.ix_(order, order)]
    ax.imshow(P2, cmap="Blues", vmin=0, vmax=1)
    # 画亚型边界
    bounds = np.where(np.diff(pd.factorize(sub.values)[0][order]) != 0)[0]
    for b in bounds:
        ax.axhline(b + 0.5, color="#C0392B", lw=1.0)
        ax.axvline(b + 0.5, color="#C0392B", lw=1.0)
    ax.set_title(f"Consensus matrix (K=2, 1000 subsamples)", fontsize=9)
    ax.set_xticks([]); ax.set_yticks([])
    nviz.save(fig, "Fig2B_consensus_matrix")


def fig2c():
    K = pd.read_csv(PROC / "consensus_k_selection.csv")
    fig, ax = plt.subplots(1, 2, figsize=(6.8, 2.6))
    ax[0].plot(K["k"], K["PAC"], "o-", color="#C0392B", lw=1.4)
    ax[0].set_xlabel("K"); ax[0].set_ylabel("PAC")
    ax[0].set_title("Proportion of ambiguous clustering", fontsize=9)
    ax[1].plot(K["k"], K["sil_consensus"], "s-", color="#2874A6", lw=1.4, label="consensus")
    ax[1].plot(K["k"], K["sil_data"], "^-", color="#7F8C8D", lw=1.2, label="raw features")
    ax[1].set_xlabel("K"); ax[1].set_ylabel("silhouette")
    ax[1].set_title("Clustering stability", fontsize=9)
    ax[1].legend(fontsize=6, frameon=False)
    for a in ax:
        for s in ("top", "right"):
            a.spines[s].set_visible(False)
        a.set_xticks(K["k"])
    nviz.save(fig, "Fig2C_k_selection")


def fig2d(Fz, sub):
    order = np.argsort(pd.factorize(sub.values)[0])
    fig, ax = plt.subplots(figsize=(6.2, 4.6))
    im = ax.imshow(Fz.values[order], aspect="auto", cmap="RdBu_r",
                   vmin=-2, vmax=2)
    cb = plt.colorbar(im, ax=ax, fraction=0.035)
    cb.ax.tick_params(labelsize=6); cb.outline.set_linewidth(0.5)
    ax.set_xticks(range(Fz.shape[1]))
    ax.set_yticks(range(Fz.shape[0]))
    ax.set_yticklabels([s.replace("GSM", "") for s in Fz.index[order]], fontsize=5.5)
    ax.set_xticklabels(Fz.columns, rotation=55, ha="right", fontsize=5.8)
    bounds = np.where(np.diff(pd.factorize(sub.values)[0][order]) != 0)[0]
    for b in bounds:
        ax.axhline(b + 0.5, color="black", lw=1.1)
    ax.set_title("Nerve-gradient features (z-scored, sorted by subtype)", fontsize=9)
    for s_ in ("top", "right"):
        ax.spines[s_].set_visible(False)
    nviz.save(fig, "Fig2D_feature_heatmap")


def fig2e(F, sub):
    keys = [("nerve_density_log", "log(nerve density)"),
            ("CD8_diff", "CD8 near-far diff (z)"),
            ("collagen_diff", "Collagen near-far diff (z)"),
            ("IFN_diff", "IFN near-far diff (z)")]
    from scipy import stats
    fig, axes = plt.subplots(1, 4, figsize=(11.5, 2.8))
    for i, (k, lab) in enumerate(keys):
        d = pd.DataFrame({"v": F[k].values, "st": sub.values})
        groups = sorted(d["st"].unique())
        data = [d.loc[d["st"] == g, "v"].values for g in groups]
        bp = axes[i].boxplot(data, tick_labels=groups, patch_artist=True, widths=0.55)
        for p_, g in zip(bp["boxes"], groups):
            p_.set_facecolor(PAL.get(g, "#95A5A6")); p_.set_alpha(0.8)
        rng = np.random.default_rng(0)
        for j, g in enumerate(groups):
            vv = d.loc[d["st"] == g, "v"].values
            axes[i].scatter(rng.normal(j + 1, 0.06, len(vv)), vv, s=6,
                            c="black", alpha=0.6, edgecolors="none")
        if len(groups) == 2:
            u, p = stats.mannwhitneyu(data[0], data[1])
            axes[i].text(0.5, 1.03, f"P={p:.3g}", transform=axes[i].transAxes,
                         ha="center", fontsize=6.5)
        axes[i].set_ylabel(lab, fontsize=7.5)
        axes[i].tick_params(labelsize=6)
        for s in ("top", "right"):
            axes[i].spines[s].set_visible(False)
    nviz.save(fig, "Fig2E_features_by_subtype")


def fig2f():
    subs = pd.read_csv(PROC / "sample_subtypes.csv")
    picks = {}
    for st, g in subs.groupby("subtype"):
        # 每个亚型挑一个有代表性的样本
        picks[st] = g["sample"].iloc[0]
    fig, axes = plt.subplots(2, len(picks), figsize=(3.1 * len(picks), 6.0))
    if len(picks) == 1:
        axes = axes.reshape(2, 1)
    for j, (st, sample) in enumerate(picks.items()):
        d = CFG.DIR_RAW / "gse210616" / sample
        try:
            ad = nio.read_visium_sample(d, sample_id=sample)
        except Exception:
            continue
        S, _ = nsig_score(ad)
        nv = nnerve.identify_nerve_spots(ad)
        mask = (nv["nerve_spatial"] == 1).values
        nviz.plot_spatial(ad, "Nervesignal", title=f"{st}\n{sample}", nerve_mask=mask,
                          size=12, ax=axes[0, j], name=None)
        nviz.plot_spatial(ad, "CD8density", title=f"CD8 density", size=12,
                          ax=axes[1, j], name=None)
    nviz.save(fig, "Fig2F_spatial_representative")


def nsig_score(ad):
    # copy=False: 直接写回 ad.obs, 否则下游 plot_spatial 找不到签名列
    return nsig.score_signatures(ad, CFG.GRADIENT_SIGNATURES, copy=False)


def fig2g():
    """两个亚型的 TLS / Hypoxia 空间对比"""
    subs = pd.read_csv(PROC / "sample_subtypes.csv")
    # 选每个亚型 nerve spot 最多的样本, 空间图更具信息量
    ss = pd.read_csv(PROC / "sample_summary.csv")
    merged = subs.merge(ss, on="sample")
    fig, axes = plt.subplots(2, 2, figsize=(6.4, 6.0))
    for i, st in enumerate(sorted(subs["subtype"].unique())):
        g = merged[merged["subtype"] == st].sort_values("n_nerve_spots", ascending=False)
        if g.empty:
            continue
        sample = g.iloc[0]["sample"]
        d = CFG.DIR_RAW / "gse210616" / sample
        try:
            ad = nio.read_visium_sample(d, sample_id=sample)
        except Exception:
            continue
        S, _ = nsig_score(ad)
        nv = nnerve.identify_nerve_spots(ad)
        mask = (nv["nerve_spatial"] == 1).values
        nviz.plot_spatial(ad, "TLSscore", title=f"{st}: TLS score", size=11,
                          ax=axes[i, 0], name=None)
        nviz.plot_spatial(ad, "Hypoxia", title=f"{st}: Hypoxia", size=11,
                          ax=axes[i, 1], name=None)
    nviz.save(fig, "Fig2G_TLS_hypoxia_by_subtype")


def main():
    subs = pd.read_csv(PROC / "sample_subtypes.csv")
    Fcols = [c for c in subs.columns if c not in ("sample", "cluster", "subtype")]
    F = subs.set_index("sample")[Fcols]
    sub = subs.set_index("sample")["subtype"]
    mu = F.mean(axis=0); sd = F.std(axis=0).replace(0, 1)
    Fz = (F - mu) / sd

    print("[Fig2] building panels ...")
    P = np.load(PROC / "consensus_matrix.npy")
    P = (P + P.T) / 2; np.fill_diagonal(P, 1.0)
    fig2a(Fz, sub)
    print("  Fig2A")
    fig2b(P, sub)
    print("  Fig2B")
    fig2c()
    print("  Fig2C")
    fig2d(Fz, sub)
    print("  Fig2D")
    fig2e(F, sub)
    print("  Fig2E")
    fig2g()
    print("  Fig2G")
    try:
        fig2f()
        print("  Fig2F")
    except Exception as e:
        print(f"  Fig2F skipped: {e}")
    print(f"[Fig2] done -> {FIG}")


if __name__ == "__main__":
    main()

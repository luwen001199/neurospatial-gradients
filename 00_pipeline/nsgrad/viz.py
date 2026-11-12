# -*- coding: utf-8 -*-
"""出版级绘图 —— 统一字号/色板/导出, 直接产出投稿可用 Figure panel。

约定:
  * 所有面板统一 dpi=300, 字号 8pt(标题 9pt), 与 Nature/Cell 系版面兼容
  * 每个 panel 单独保存为 png + pdf, 便于在 Illustrator 里重排组合
  * 空间散点图的坐标轴统一用微米, 且 X/Y 等比例(aspect='equal'),
    否则,spot 分布会被拉歪, 距离/形态的结论在视觉上失真
"""
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib import rcParams
from scipy.spatial import cKDTree

from .config import FIG, DIR_FIG

rcParams.update({
    "font.size": FIG["font_size"],
    "axes.titlesize": FIG["title_size"],
    "axes.labelsize": FIG["font_size"],
    "xtick.labelsize": FIG["font_size"] - 1,
    "ytick.labelsize": FIG["font_size"] - 1,
    "legend.fontsize": FIG["font_size"] - 1,
    "figure.dpi": FIG["dpi"],
    "savefig.dpi": FIG["dpi"],
    "savefig.bbox": "tight",
    "font.family": "DejaVu Sans",
    "axes.linewidth": 0.6,
    "pdf.fonttype": 42,     # 保证 PDF 里文字可编辑, 期刊要求
    "ps.fonttype": 42,
})


def save(fig, name, folder=DIR_FIG):
    folder = Path(folder)
    folder.mkdir(parents=True, exist_ok=True)
    fig.savefig(folder / f"{name}.png", dpi=FIG["dpi"], bbox_inches="tight")
    fig.savefig(folder / f"{name}.pdf", bbox_inches="tight")
    plt.close(fig)
    return folder / f"{name}.png"


# ----------------------------------------------------------------- 空间图
def plot_spatial(ad, color, title="", cmap=None, cmap_div=True,
                 nerve_mask=None, ax=None, size=14, alpha=0.9,
                 vmin=None, vmax=None, colorbar=True, name=None):
    """单张空间特征图。可选择叠加神经 spot 的轮廓十字标。"""
    cmap = cmap or (FIG["cmap_diverging"] if cmap_div else FIG["cmap_spatial"])
    if ax is None:
        _, ax = plt.subplots(figsize=(3.4, 3.2))
    x, y = ad.obs["x_um"].values, ad.obs["y_um"].values
    v = ad.obs[color].values.astype(float) if isinstance(color, str) else np.asarray(color, float)
    if vmin is None:
        vmin = np.nanpercentile(v, 2)
    if vmax is None:
        vmax = np.nanpercentile(v, 98)
    sc = ax.scatter(x, -y, c=v, cmap=cmap, s=size, alpha=alpha,
                    vmin=vmin, vmax=vmax, edgecolors="none", rasterized=True)
    if nerve_mask is not None and np.asarray(nerve_mask).sum() > 0:
        nm = np.asarray(nerve_mask)
        ax.scatter(x[nm], -y[nm], facecolors="none", edgecolors="black",
                   s=size * 2.2, linewidths=0.45, label="nerve")
    ax.set_aspect("equal")
    ax.axis("off")
    ax.set_title(title, fontsize=FIG["title_size"], pad=3)
    if colorbar:
        cb = plt.colorbar(sc, ax=ax, fraction=0.035, pad=0.02)
        cb.ax.tick_params(labelsize=FIG["font_size"] - 1.5)
        cb.outline.set_linewidth(0.5)
    if name:
        save(ax.get_figure(), name)
    return ax


def plot_decay_curves(prof_long, metric, hue="sample", ax=None, name=None,
                      max_samples=40, show_legend=False, jitter_x=True):
    """多样本的环带均值曲线叠加。横轴=平均距离, 纵轴=指标均值。"""
    if ax is None:
        _, ax = plt.subplots(figsize=(3.2, 2.8))
    sub = prof_long.dropna(subset=["mean_dist_um", f"{metric}_mean"])
    samples = sorted(sub[hue].unique())[:max_samples]
    cmap = plt.get_cmap("tab20")
    rng = np.random.default_rng(0)
    for i, s in enumerate(samples):
        d = sub[sub[hue] == s].sort_values("mean_dist_um")
        if len(d) < 2:
            continue
        xx = d["mean_dist_um"].values
        yy = d[f"{metric}_mean"].values
        if jitter_x:
            xx = xx + rng.normal(0, 1.5, len(xx))
        ax.plot(xx, yy, "-o", ms=2.2, lw=0.7, alpha=0.55,
                color=cmap(i % 20), rasterized=True)
    ax.axhline(0, color="grey", lw=0.5, ls=":")
    ax.set_xlabel("distance to nearest nerve (um)")
    ax.set_ylabel(f"{metric}\n(z-scored)")
    ax.set_title(f"{metric}: distance decay", fontsize=FIG["title_size"])
    ax.tick_params(width=0.5)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    if show_legend:
        ax.legend(fontsize=5, ncol=2, frameon=False)
    if name:
        save(ax.get_figure(), name)
    return ax


def plot_fit_single(dist, y, params, metric, sample="", ax=None, name=None):
    """单个样本的: 散点 + 聚合箱均值 + 拟合曲线。"""
    if ax is None:
        _, ax = plt.subplots(figsize=(3.4, 2.9))
    m = np.isfinite(dist) & np.isfinite(y)
    d, yy = dist[m], y[m]
    ax.scatter(d, yy, s=1.6, c="#BDC3C7", alpha=0.5, edgecolors="none",
               rasterized=True, label="spots")
    from .gradient import aggregate_distance_curve
    agg = aggregate_distance_curve(d, yy)
    if agg is not None:
        bx, by, bw = agg
        ax.plot(bx, by, "o-", color="#C0392B", ms=4, lw=1.2,
                label="binned mean", zorder=3)
    if params and np.isfinite(params.get("b", np.nan)):
        xs = np.linspace(0, np.nanpercentile(d, 99), 200)
        ys = params["a"] * np.exp(-params["b"] * xs) + params["c"]
        ax.plot(xs, ys, "-", color="#154360", lw=1.6, zorder=4,
                label=f"y={params['a']:.2f}·exp(-{params['b']:.4f}d){params['c']:+.2f}")
        ax.axvline(1.0 / params["b"], color="#154360", ls="--", lw=0.8, alpha=0.7)
        ax.text(1.0 / params["b"], np.nanpercentile(yy, 92),
                f" 1/b={1/params['b']:.0f}um", fontsize=6, color="#154360")
    ax.set_xlabel("distance to nearest nerve (um)")
    ax.set_ylabel(f"{metric} (z)")
    ax.set_title(f"{metric} — {sample}  $R^2$={params.get('r2', float('nan')):.2f}",
                 fontsize=FIG["title_size"])
    ax.legend(fontsize=5.5, frameon=False, loc="best")
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    if name:
        save(ax.get_figure(), name)
    return ax


def plot_param_distribution(df, col, group=None, ax=None, name=None, title=""):
    """参数分布: 直方图 + 按组着色的箱线图组合。"""
    if ax is None:
        _, ax = plt.subplots(figsize=(3.0, 2.6))
    v = pd.to_numeric(df[col], errors="coerce").dropna()
    if len(v) == 0:
        ax.text(0.5, 0.5, "no data", ha="center", transform=ax.transAxes)
        if name:
            save(ax.get_figure(), name)
        return ax
    lo, hi = float(np.percentile(v, [1, 99]))
    v = v[(v >= lo) & (v <= hi)]
    ax.hist(v, bins=min(18, max(6, len(v) // 3)), color="#7FB3D5",
            edgecolor="white", linewidth=0.5, alpha=0.9)
    ax.axvline(v.median(), color="#C0392B", lw=1.3, ls="--",
               label=f"median={v.median():.3g}")
    ax.set_xlabel(col)
    ax.set_ylabel("n samples")
    ax.set_title(title or col, fontsize=FIG["title_size"])
    ax.legend(fontsize=6, frameon=False)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    if name:
        save(ax.get_figure(), name)
    return ax


def plot_box_by_group(df, col, group, ax=None, name=None, palette=None, ylabel=None):
    palette = palette or FIG["palette_subtype"]
    if ax is None:
        _, ax = plt.subplots(figsize=(2.9, 2.8))
    d = df[[col, group]].dropna()
    if d.empty:
        ax.text(0.5, 0.5, "no data", ha="center", transform=ax.transAxes)
        if name:
            save(ax.get_figure(), name)
        return ax
    order = sorted(d[group].unique())
    data = [d.loc[d[group] == g, col].values for g in order]
    bp = ax.boxplot(data, labels=order, patch_artist=True, widths=0.62,
                    medianprops=dict(color="black", lw=1.2),
                    flierprops=dict(marker="o", ms=2, alpha=0.5))
    for patch, g in zip(bp["boxes"], order):
        patch.set_facecolor(palette.get(g, "#95A5A6"))
        patch.set_alpha(0.75)
        patch.set_linewidth(0.6)
    # 叠加散点, 让样本量一目了然
    rng = np.random.default_rng(1)
    for i, g in enumerate(order):
        vv = d.loc[d[group] == g, col].values
        ax.scatter(rng.normal(i + 1, 0.055, len(vv)), vv, s=4,
                   c="black", alpha=0.55, zorder=3, edgecolors="none")
    ax.set_xticklabels(order, rotation=32, ha="right", fontsize=6.5)
    ax.set_ylabel(ylabel or col)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    if name:
        save(ax.get_figure(), name)
    return ax


def plot_scatter_corr(x, y, labels=None, ax=None, name=None,
                      xlabel="", ylabel="", title="", annotate=False):
    from scipy import stats
    if ax is None:
        _, ax = plt.subplots(figsize=(3.0, 2.8))
    x = np.asarray(x, float); y = np.asarray(y, float)
    m = np.isfinite(x) & np.isfinite(y)
    ax.scatter(x[m], y[m], s=22, c="#2874A6", alpha=0.8, edgecolors="white",
               linewidths=0.4)
    if m.sum() >= 3:
        sl, ic, r, p, se = stats.linregress(x[m], y[m])
        xs = np.linspace(np.nanmin(x[m]), np.nanmax(x[m]), 50)
        ax.plot(xs, sl * xs + ic, color="#C0392B", lw=1.2)
        rs, ps = stats.spearmanr(x[m], y[m])
        ax.text(0.05, 0.94, f"Spearman r={rs:.2f}\np={ps:.2g}",
                transform=ax.transAxes, fontsize=6.5, va="top")
    if annotate and labels is not None:
        for xi, yi, lb in zip(x, y, labels):
            if np.isfinite(xi) and np.isfinite(yi):
                ax.annotate(str(lb)[:8], (xi, yi), fontsize=4.5, alpha=0.7)
    ax.set_xlabel(xlabel); ax.set_ylabel(ylabel)
    ax.set_title(title, fontsize=FIG["title_size"])
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    if name:
        save(ax.get_figure(), name)
    return ax


def plot_heatmap(mat, ax=None, name=None, cmap="RdBu_r", title="",
                 xlabel="", ylabel="", vmin=None, vmax=None, xticklabels=True):
    if ax is None:
        _, ax = plt.subplots(figsize=(4.2, 3.2))
    A = np.asarray(mat, float)
    if vmin is None:
        vmin = np.nanpercentile(A, 2)
    if vmax is None:
        vmax = np.nanpercentile(A, 98)
    im = ax.imshow(A, aspect="auto", cmap=cmap, vmin=vmin, vmax=vmax,
                   interpolation="nearest")
    cb = plt.colorbar(im, ax=ax, fraction=0.035, pad=0.02)
    cb.ax.tick_params(labelsize=6)
    cb.outline.set_linewidth(0.5)
    ax.set_title(title, fontsize=FIG["title_size"])
    ax.set_xlabel(xlabel); ax.set_ylabel(ylabel)
    if isinstance(mat, pd.DataFrame):
        ax.set_yticks(range(len(mat.index)))
        ax.set_yticklabels(mat.index, fontsize=6)
        if xticklabels:
            ax.set_xticks(range(len(mat.columns)))
            ax.set_xticklabels(mat.columns, rotation=60, ha="right", fontsize=5.5)
    ax.tick_params(width=0.4)
    if name:
        save(ax.get_figure(), name)
    return ax

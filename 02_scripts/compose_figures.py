# -*- coding: utf-8 -*-
"""把单面板 PNG 组合成投稿用的整图 (Figure 1 / 2 / 5)。

Fig3 / Fig4 已是组合图, 只做等比重采样到目标宽度。
布局规则: 每行若干面板等高排列, 行内宽度和 = 目标宽度; 字母标签按行内
从左到右、逐行向下的顺序添加 (面板内部已有字母的, 不重复添加)。
输出: 03_results/figures/composed/FigN_composite.png / .pdf
"""
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from PIL import Image

ROOT = Path(r"D:/嗜神经浸润")
FIGD = ROOT / "03_results" / "figures"
OUT = FIGD / "composed"
OUT.mkdir(exist_ok=True)

TARGET_W_IN = 180 / 25.4          # 180 mm
MAX_H_IN = 247 / 25.4             # 单页最大高度 (A4 减页边距)
GAP_IN = 0.04                     # 行/列间距
PAD_TOP_IN = 0.10                 # 每行顶部留给字母标签的空间
PAD_LEFT_IN = 0.09                # 左边距, 字母标签落在画布内 (修复贴边裁切)


def load(name):
    im = Image.open(FIGD / name)
    return np.asarray(im.convert("RGB")), im.size[0] / im.size[1]


def compose(rows, out_name, add_letters=True):
    """rows: list of (row_items, width_frac) 或直接 row_items。
    row_items: list of (filename, letter_or_None)。letter=None 表示不加标签。
    width_frac: 该行占目标宽度的比例 (默认 1.0)。"""
    norm = []
    for r in rows:
        if isinstance(r, tuple) and len(r) == 2 and isinstance(r[1], float):
            norm.append(r)
        else:
            norm.append((r, 1.0))
    # 1) 计算自然高度
    row_specs = []
    for row, wf in norm:
        aspects = [load(f)[1] for f, _ in row]
        h_nat = (TARGET_W_IN * wf) / sum(aspects)
        row_specs.append(h_nat)
    total_h = sum(row_specs) + len(rows) * (GAP_IN + PAD_TOP_IN) + PAD_TOP_IN
    scale = min(1.0, MAX_H_IN / total_h)
    content_w = TARGET_W_IN * scale
    padl = PAD_LEFT_IN * scale if add_letters else 0.0
    W = content_w + padl
    row_specs = [h * scale for h in row_specs]
    gap, pad = GAP_IN * scale, PAD_TOP_IN * scale

    fig_h = sum(row_specs) + len(rows) * (gap + pad) + pad  # 首行上方额外一段 pad, 字母不出画布上缘
    fig = plt.figure(figsize=(W, fig_h))
    y = fig_h - pad
    for (row, wf), rh in zip(norm, row_specs):
        y -= pad
        top = y
        y -= rh
        imgs = [load(f)[0] for f, _ in row]
        aspects = [load(f)[1] for f, _ in row]
        widths = [a * rh for a in aspects]
        total_w = sum(widths) + (len(row) - 1) * gap
        x = padl + (content_w - total_w) / 2
        for (fname, letter), img, w in zip(row, imgs, widths):
            ax = fig.add_axes([x / W, y / fig_h, w / W, rh / fig_h])
            ax.imshow(img, interpolation="lanczos", aspect="auto")
            ax.axis("off")
            if add_letters and letter:
                fig.text(x / W, top / fig_h + 0.002, letter, ha="left",
                         va="bottom", fontsize=9, fontweight="bold",
                         family="DejaVu Sans")
            x += w + gap
        y -= gap

    fig.savefig(OUT / f"{out_name}.png", dpi=300, facecolor="white")
    fig.savefig(OUT / f"{out_name}.pdf", facecolor="white")
    plt.close(fig)
    print(f"  {out_name}: {W*25.4:.0f} x {fig_h*25.4:.0f} mm")


F = "Fig1"
compose(
    [
        [("Fig1A_pipeline_schematic.png", "A")],
        [("Fig1B_nervesignal_map.png", "B"),
         ("Fig1C_distance_field.png", "C"),
         ("Fig1D_collagen_myCAF_CD8_maps.png", "D")],
        [("Fig1E_fit_Nervesignal.png", "E"),
         ("Fig1E_fit_CollagenIndex.png", None),
         ("Fig1E_fit_CD8density.png", None)],
        [("Fig1F_decaycurves_Nervesignal.png", "F"),
         ("Fig1F_decaycurves_CollagenIndex.png", None),
         ("Fig1F_decaycurves_CD8density.png", None),
         ("Fig1H_nerve_call_validation.png", "G")],
        [("Fig1G1_r2_heatmap.png", "H"),
         ("Fig1G2_decay_radius_distribution.png", None),
         ("Fig1G3_identifiability.png", None)],
    ],
    "Figure1_composite",
)

compose(
    [
        [("Fig2A_volcano.png", None), ("Fig2B_top_genes.png", None),
         ("Fig2C_category_enrichment.png", None)],
        [("Fig2D_ring_profile.png", None), ("Fig2E_npp_map.png", None),
         ("Fig2F_gradient_rho_boxplot.png", None)],
        [("Fig2G_axis_ring_profiles.png", None), ("Fig2H_subtype_comparison.png", None)],
    ],
    "Figure2_composite",
    add_letters=False,
)

compose(
    [
        [("Fig5AB_km_npp_tertiles.png", None)],
        [("Fig5CD_forest_nerve_axes.png", None)],
        ([("Fig5E_npp_vs_caf.png", None), ("Fig5F_km_npp_resid.png", None),
          ("Fig5G_permutation_null.png", None)], 1.0),
        ([("Fig5H_subtype_forest.png", None)], 0.5),
    ],
    "Figure5_composite",
    add_letters=False,
)

# Fig3 / Fig4 已经是组合图, 等比缩放到 180 mm 宽
for src, dst in [("Fig3_perineural_niche", "Figure3_composite"),
                 ("Fig4_trajectory_and_modules", "Figure4_composite")]:
    im = Image.open(FIGD / f"{src}.png").convert("RGB")
    w = int(180 / 25.4 * 300)
    h = int(im.size[1] * w / im.size[0])
    im.resize((w, h), Image.LANCZOS).save(OUT / f"{dst}.png", dpi=(300, 300))
    fig = plt.figure(figsize=(w / 300, h / 300))
    ax = fig.add_axes([0, 0, 1, 1]); ax.axis("off")
    ax.imshow(np.asarray(im), interpolation="lanczos", aspect="auto")
    fig.savefig(OUT / f"{dst}.pdf", facecolor="white")
    plt.close(fig)
    print(f"  {dst}: {w/300*25.4:.0f} x {h/300*25.4:.0f} mm")

print("done ->", OUT)

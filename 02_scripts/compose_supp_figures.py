# -*- coding: utf-8 -*-
"""把 SuppFigS7 (A-G) 与 SuppFigS8 (A-D) 的单面板 PNG 组合成整图。
布局规则同 compose_figures.py: 每行面板等高、行内宽度比例缩放,
行首面板加字母标签。输出 composed/SuppFigS7_composite.png/.pdf 等。
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

TARGET_W_IN = 180 / 25.4
MAX_H_IN = 247 / 25.4
GAP_IN = 0.04
PAD_TOP_IN = 0.10
PAD_LEFT_IN = 0.09   # 左边距, 字母标签落在画布内 (修复贴边裁切)


def load(name):
    im = Image.open(FIGD / name)
    return np.asarray(im.convert("RGB")), im.size[0] / im.size[1]


def compose(rows, out_name, add_letters=True):
    norm = []
    for r in rows:
        if isinstance(r, tuple) and len(r) == 2 and isinstance(r[1], float):
            norm.append(r)
        else:
            norm.append((r, 1.0))
    row_specs = []
    for row, wf in norm:
        aspects = [load(f)[1] for f, _ in row]
        row_specs.append((TARGET_W_IN * wf) / sum(aspects))
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


# S7: 空间亚型稳健性 (A PCA, B consensus, C k-selection, D feature heatmap,
#     E features-by-subtype 条带, F 空间代表图, G TLS/缺氧)
compose(
    [
        [("SuppFigS7A_pca_subtypes.png", "A"),
         ("SuppFigS7B_consensus_matrix.png", "B"),
         ("SuppFigS7C_k_selection.png", "C")],
        [("SuppFigS7D_feature_heatmap.png", "D"),
         ("SuppFigS7E_features_by_subtype.png", "E")],
        ([("SuppFigS7F_spatial_representative.png", "F"),
          ("SuppFigS7G_TLS_hypoxia_by_subtype.png", "G")], 0.72),
    ],
    "SuppFigS7_composite",
)

# S8: GSE12276 器官特异性复发 (A CIF, B forest, C 敏感性, D 按部位评分)
compose(
    [
        [("SuppFigS8A_cif_brain_npp.png", "A"),
         ("SuppFigS8B_forest_brainmet.png", "B"),
         ("SuppFigS8C_sensitivity_endpoint.png", "C")],
        ([("SuppFigS8D_scores_by_site.png", "D")], 0.82),
    ],
    "SuppFigS8_composite",
    add_letters=False,   # S8 各面板内部标题已带 A-D 前缀, 不重复加
)

print("done ->", OUT)

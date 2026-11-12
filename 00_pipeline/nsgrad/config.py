# -*- coding: utf-8 -*-
"""基因签名与全局参数 —— 全流水线唯一来源。

所有基因列表以 HGNC symbol 为准, 任何修改都会自动同步到下游每个模块,
避免'方法写一套、代码实现另一套'这类审稿人最爱的硬伤。
"""

# ----------------------------------------------------------------- 空间尺度
# Visium 标准: spot 直径 55 um, spot 中心间距 100 um
SPOT_DIAMETER_UM = 55.0
SPOT_PITCH_UM = 100.0

# 环形缓冲区边界 (um) —— 方案规定 50/100/200/400
RING_EDGES_UM = [0.0, 50.0, 100.0, 200.0, 400.0]
RING_LABELS = ["r1_0_50", "r2_50_100", "r3_100_200", "r4_200_400", "out_gt400"]

# ----------------------------------------------------------------- 神经标记
NERVE_MARKERS = {
    "sensory": ["NGF", "UCHL1", "TRPV1", "CALCA", "CALCB", "S100B", "SCN9A", "SCN10A"],
    "sympathetic": ["TH", "DBH", "CHGA", "SLC6A2"],
    "parasympathetic": ["SLC18A3", "CHAT", "SLC5A7"],
    "structural": ["TUBB3", "MAP2", "NEFL", "NEFM", "PLP1", "MPZ", "SYP", "UCHL1"],
}
ALL_NERVE_GENES = sorted({g for v in NERVE_MARKERS.values() for g in v})

# ----------------------------------------------------------------- 梯度读数
# 每个指标 = 一组基因的均值(z-scored gene-wise 后取均值), 用于拟合衰减曲线
GRADIENT_SIGNATURES = {
    "CollagenIndex": ["COL1A1", "COL1A2", "COL3A1", "COL5A1", "COL11A1", "FN1"],
    "myCAFscore":    ["ACTA2", "TAGLN", "MYL9", "FAP", "PDGFRB", "RAMP1", "CTHRC1"],
    "CD8density":    ["CD8A", "CD8B", "CD3D", "CD3E", "GZMA", "NKG7"],
    "TLSscore":      ["CCL21", "CCL19", "CXCL13", "LAMP3", "MS4A1", "CD79A"],
    "Nervesignal":   ["NGF", "UCHL1", "CALCA", "CALCB", "TRPV1", "S100B", "TUBB3", "PLP1"],
    "Hypoxia":       ["HIF1A", "SLC2A1", "ADM", "EGLN3", "CA9", "LDHA", "PGK1", "VEGFA"],
    "EMT":           ["TWIST1", "SNAI1", "SNAI2", "ZEB1", "ZEB2", "CDH2", "VIM"],
    "IFNsignature":  ["STAT1", "IRF1", "GBP1", "CXCL9", "CXCL10", "TAP1", "B2M"],
}

# ----------------------------------------------------------------- 分型特征
SUBTYPE_FEATURES = [
    "nerve_density",          # 神经 spot 密度
    "decay_b_median",         # 衰减常数中位数
    "immune_exclusion_idx",   # 免疫排斥指数 = 近神经 CD8 / 远神经 CD8
    "collagen_barrier_thick", # 胶原屏障厚度
    "nerve_prox_caf",         # 神经近邻 CAF 占比
    "nerve_prox_tcell",       # 神经近邻 T 细胞占比
]

# ----------------------------------------------------------------- 质控阈值
QC = {
    "min_umi_per_spot": 500,      # 低 UMI spot 过滤
    "min_genes_per_spot": 200,
    "mito_max_pct": 30.0,
    "nerve_q_default": 0.90,      # 方案要求默认 75% 分位数
    "nerve_q_grid": [0.85, 0.90, 0.95],   # 敏感性分析轮换阈值
    "nerve_min_markers": 2,       # 至少 2 个同类神经标记同时高表达
}

# ----------------------------------------------------------------- 绘图风格
FIG = {
    "dpi": 300,
    "font_size": 8,
    "title_size": 9,
    "palette_subtype": {
        "Neural-Barriered": "#C0392B",
        "Neural-Inflamed":  "#2E86C1",
        "Neural-Cold":      "#7F8C8D",
    },
    "cmap_spatial": "viridis",
    "cmap_diverging": "RdBu_r",
}

# ----------------------------------------------------------------- 目录
from pathlib import Path

ROOT = Path(r"D:/嗜神经浸润")
DIR_RAW = ROOT / "01_data" / "raw"
DIR_PROC = ROOT / "01_data" / "processed"
DIR_FIG = ROOT / "03_results" / "figures"
DIR_TAB = ROOT / "03_results" / "tables"
for _d in (DIR_PROC, DIR_FIG, DIR_TAB):
    _d.mkdir(parents=True, exist_ok=True)

"""神经周围微环境的"功能读数"基因集库 —— 供方案一/二/三/四共用。

设计原则
--------
1. **细胞类型标记为互斥(mutually exclusive)**: 每个 gene symbol 只归属于一个
   细胞类型(按 CELLTYPE_PRIORITY 顺序分配)。否则 myCAF 与 pericyte 共享 ACTA2、
   myCAF 与 myoepithelium 共享 MYH11 会让"丰度梯度"互相污染。
   来源: 乳腺 scRNA-seq 常规标记 (Wu 2021 Nat Genet; Pal 2021 EMBO J;
   Bartoschek 2018 Nat Commun; Costa 2018 Cancer Cell 的 myCAF/iCAF/dCAF 命名)。
2. **配体-受体对为 curated**: CellPhoneDB 全库在本环境不可下载(连接超时),
   因此采用文献 curated 的四条轴(神经->基质 / 神经->免疫 / 基质->免疫 /
   基质->肿瘤), 每条轴的成员都有明确的 PNI 或 CAF 生物学文献支撑。
   这是本分析的**先验知识来源**, 在论文中必须如此声明。
3. **通路来自真 MSigDB** (v2023.2.Hs), 按预定义关键词规则筛选, 不做人工挑选。
"""
from __future__ import annotations

import re
from pathlib import Path

import numpy as np

GMT = Path(r"D:/嗜神经浸润/01_data/raw/annotation/msigdb.v2023.2.Hs.symbols.gmt")

# ----------------------------------------------------------------- 细胞类型
# 顺序 = 互斥分配优先级(越靠前越"专有")
CELLTYPE_PRIORITY = [
    "Neuronal_Schwann", "Plasma_B", "B_lymphocyte", "Mast", "Neutrophil",
    "NK", "CD8_T", "DC", "Macrophage_monocyte", "T_lymphocyte",
    "Endothelial", "Pericyte_SMC", "Adipocyte",
    "CAF_myofibroblastic", "CAF_inflammatory_dev",
    "Tumor_proliferative", "Normal_epithelial",
]

CELLTYPE_MARKERS = {
    "Neuronal_Schwann": ["UCHL1", "TUBB3", "MAP2", "PLP1", "MPZ", "PMP22", "S100B",
                         "NGFR", "SOX10", "SCN7A", "CDH19", "MAL", "GAP43", "SYP",
                         "NEFL", "NEFM", "SCN9A", "SCN10A", "PRX", "MBP", "EGR2"],
    "Plasma_B": ["IGHG1", "IGHG2", "IGHG3", "IGHG4", "IGHA1", "MZB1", "JCHAIN",
                 "XBP1", "SDC1", "DERL3", "SSR4", "FKBP11", "IGLC1", "IGKC"],
    "B_lymphocyte": ["MS4A1", "CD79A", "CD79B", "CD19", "PAX5", "BANK1", "IGHM",
                     "IGHD", "CD37", "TNFRSF13C", "FCRL1", "CD22"],
    "Mast": ["TPSAB1", "TPSB2", "CPA3", "MS4A2", "CMA1", "HDC", "GATA2", "MS4A3"],
    "Neutrophil": ["FCGR3B", "CSF3R", "S100A8", "S100A9", "MNDA", "FPR1", "CXCR2",
                   "CEACAM8", "MMP8", "LTF", "CAMP", "DEFA1"],
    "NK": ["KLRD1", "KLRF1", "GNLY", "NCR1", "FGFBP2", "SPON2", "KLRC1", "KLRC2",
           "KIR2DL1", "NCR3"],
    "CD8_T": ["CD8A", "CD8B", "GZMA", "GZMB", "GZMK", "PRF1", "GZMH", "NKG7",
              "KLRG1", "CTSW"],
    "DC": ["CD1C", "FCER1A", "CLEC10A", "LAMP3", "CCR7", "CLEC9A", "CADM1", "CCL17",
           "IRF8", "ENHO", "CST3"],
    "Macrophage_monocyte": ["CD68", "CSF1R", "AIF1", "LYZ", "C1QA", "C1QB", "C1QC",
                            "TYROBP", "FCER1G", "CD163", "MRC1", "MSR1", "CD14",
                            "SIGLEC1", "ITGAM"],
    "T_lymphocyte": ["CD3D", "CD3E", "CD3G", "CD2", "TRAC", "TRBC2", "IL7R", "SKAP1",
                     "THEMIS", "LCK", "TRAT1", "ITK", "CD28"],
    "Endothelial": ["PECAM1", "VWF", "CDH5", "PLVAP", "ENG", "CLDN5", "FLT1", "EMCN",
                    "ESAM", "ROBO4", "SOX17", "EGFL7", "CD36"],
    "Pericyte_SMC": ["RGS5", "PDGFRB", "NOTCH3", "KCNJ8", "ABCC9", "DES", "CALD1",
                     "ACTA2", "MYH11", "CNN1", "TAGLN", "MYL9"],
    "Adipocyte": ["ADIPOQ", "FABP4", "PLIN1", "LEP", "CIDEC", "LPL", "PLIN4", "GPD1",
                  "ADH1B", "FABP5"],
    # myCAF: 收缩型/基质合成型 (Bartoschek 2018; Costa 2018 "myCAF")
    "CAF_myofibroblastic": ["POSTN", "COL11A1", "LRRC15", "MMP11", "COMP", "COL10A1",
                            "THBS2", "SULF1", "ASPN", "FAP", "COL8A1", "COL5A2",
                            "TNC", "FN1", "INHBA", "GREM1", "VCAN", "COL12A1",
                            "MXRA8", "SERPINE1", "TAGLN2"],
    # iCAF/dCAF: 炎性 + 发育样 (Costa 2018 "iCAF"; Bartoschek 2018 "dCAF")
    "CAF_inflammatory_dev": ["CXCL12", "CFD", "DCN", "LUM", "PDGFRA", "FBLN1", "FBLN2",
                             "IL6", "LIF", "CCL2", "CXCL1", "CXCL2", "SFRP2", "C7",
                             "IGF1", "PI16", "CD34", "SERPINF1", "GSN", "MFAP5",
                             "MGP", "GPC3", "FGF7", "FGF10", "WNT5A"],
    "Tumor_proliferative": ["MKI67", "TOP2A", "CCNB1", "PCNA", "TYMS", "CDK1", "BIRC5",
                            "UBE2C", "AURKA", "MYBL2", "RRM2", "CENPF", "CDC20",
                            "TYMS", "BUB1"],
    "Normal_epithelial": ["KRT5", "KRT14", "KRT15", "KRT17", "OXTR", "KRT6A", "TP63",
                          "KRT23", "SLPI"],
}


def build_celltype_sets(genes_available: set[str] | None = None):
    """返回 (names, dict name->list[str]), 已做互斥去重与可用性过滤。"""
    used: set[str] = set()
    out: dict[str, list[str]] = {}
    for ct in CELLTYPE_PRIORITY:
        gs = [g for g in CELLTYPE_MARKERS[ct] if g not in used]
        if genes_available is not None:
            gs = [g for g in gs if g in genes_available]
        if len(gs) >= 3:
            out[ct] = gs
            used.update(gs)
    return list(out.keys()), out


# ------------------------------------------------------------- 配体-受体对
# (ligand, receptor, axis, 简短生物学注释)
LR_PAIRS = [
    # --- 神经/神经元 -> 基质 ---
    ("NGF", "NTRK1", "neuro->stroma", "NGF/TrkA 促 CAF 活化与侵袭"),
    ("NGF", "NGFR", "neuro->stroma", "NGF/p75NTR 促侵袭与干性"),
    ("BDNF", "NTRK2", "neuro->stroma", "BDNF/TrkB"),
    ("NTF3", "NTRK3", "neuro->stroma", "NT-3/TrkC"),
    ("S100B", "AGER", "neuro->stroma", "S100B/RAGE 促 NF-kB 与纤维化"),
    ("NRG1", "ERBB3", "neuro->stroma", "Neuregulin-1/HER3"),
    ("NRG1", "ERBB4", "neuro->stroma", "Neuregulin-1/HER4"),
    ("SEMA3A", "NRP1", "neuro->stroma", "Semaphorin 3A 轴突导向/排斥"),
    ("SEMA3A", "PLXNA1", "neuro->stroma", "Sema3A/PlexinA1"),
    ("SEMA4D", "PLXNB1", "neuro->stroma", "Sema4D/PlexinB1"),
    ("SLIT2", "ROBO1", "neuro->stroma", "Slit2/Robo1 导向"),
    ("EFNB3", "EPHB1", "neuro->stroma", "Ephrin-B3/EphB1"),
    ("GAL", "GALR1", "neuro->stroma", "Galanin/GALR1"),
    ("NPY", "NPY1R", "neuro->stroma", "NPY/Y1 促血管与纤维化"),
    ("NPY", "NPY2R", "neuro->stroma", "NPY/Y2"),
    ("SST", "SSTR2", "neuro->stroma", "Somatostatin/SSTR2"),
    ("TAC1", "TACR1", "neuro->stroma", "P 物质/NK1 促基质活化"),
    ("ADM", "CALCRL", "neuro->stroma", "ADM/CLR"),
    ("PTN", "NRP1", "neuro->stroma", "Pleiotrophin/NRP1"),
    ("PSAP", "GPR37L1", "neuro->stroma", "Prosaposin/GPR37L1"),
    # --- 神经 -> 免疫 (神经-免疫轴) ---
    ("CALCA", "CALCRL", "neuro->immune", "CGRP/CLR 抑制 DC 成熟与 Th1"),
    ("CALCA", "RAMP1", "neuro->immune", "CGRP 受体组分 RAMP1"),
    ("VIP", "VIPR1", "neuro->immune", "VIP 诱导耐受性 DC / Treg"),
    ("VIP", "VIPR2", "neuro->immune", "VIP/VIPR2"),
    ("TAC1", "TACR1", "neuro->immune", "P 物质驱动神经源性炎症"),
    ("CHAT", "CHRNA7", "neuro->immune", "胆碱能抗炎通路"),
    ("NPY", "NPY1R", "neuro->immune", "NPY/Y1 抑制 NK 与 T 细胞"),
    ("ADRB2", "ARRB1", "neuro->immune", "beta2-肾上腺素受体信号"),
    ("CCL2", "CCR2", "neuro->immune", "单核细胞趋化"),
    ("CX3CL1", "CX3CR1", "neuro->immune", "Fractalkine 轴"),
    ("S100B", "AGER", "neuro->immune", "S100B/RAGE 促炎"),
    ("SEMA4D", "CD72", "neuro->immune", "Sema4D/CD72 免疫调节"),
    # --- 基质 -> 免疫 (排斥/屏障) ---
    ("POSTN", "ITGAV", "stroma->immune", "Periostin/integrin av 促排斥"),
    ("POSTN", "ITGB3", "stroma->immune", "Periostin/integrin beta3"),
    ("THBS1", "CD47", "stroma->immune", "TSP1/CD47 抑制免疫"),
    ("TGFB1", "TGFBR1", "stroma->immune", "TGF-beta 免疫抑制主轴"),
    ("TGFB1", "TGFBR2", "stroma->immune", "TGF-beta 受体 2"),
    ("COL1A1", "ITGA1", "stroma->immune", "胶原/整合素屏障"),
    ("COL1A1", "ITGA2", "stroma->immune", "胶原/整合素 alpha2"),
    ("COL1A1", "DDR1", "stroma->immune", "胶原/DDR1"),
    ("LAMA4", "ITGA6", "stroma->immune", "层粘连蛋白/整合素 a6"),
    ("FN1", "ITGA5", "stroma->immune", "纤连蛋白/整合素 a5"),
    ("SPP1", "CD44", "stroma->immune", "骨桥蛋白/CD44"),
    ("SPP1", "ITGAV", "stroma->immune", "骨桥蛋白/整合素 av"),
    ("LGALS9", "HAVCR2", "stroma->immune", "Galectin-9/TIM-3"),
    ("CD274", "PDCD1", "stroma->immune", "PD-L1/PD-1"),
    ("NT5E", "ADORA2A", "stroma->immune", "CD73-腺苷通路"),
    ("TNC", "TLR4", "stroma->immune", "Tenascin-C/TLR4"),
    ("CXCL12", "CXCR4", "stroma->immune", "CXCL12/CXCR4 免疫滞留"),
    ("CXCL12", "ACKR3", "stroma->immune", "CXCL12/CXCR7"),
    # --- 基质/CAF -> 肿瘤 ---
    ("HGF", "MET", "stroma->tumor", "HGF/c-Met"),
    ("FGF2", "FGFR1", "stroma->tumor", "FGF2/FGFR1"),
    ("FGF7", "FGFR2", "stroma->tumor", "FGF7/FGFR2"),
    ("IGF1", "IGF1R", "stroma->tumor", "IGF1/IGF1R"),
    ("PDGFA", "PDGFRA", "stroma->tumor", "PDGF-A/PDGFRa"),
    ("PDGFB", "PDGFRB", "stroma->tumor", "PDGF-B/PDGFRb"),
    ("WNT5A", "FZD2", "stroma->tumor", "WNT5A/FZD2"),
    ("WNT5A", "ROR2", "stroma->tumor", "WNT5A/ROR2"),
    ("IL6", "IL6R", "stroma->tumor", "IL6/IL6R"),
    ("TGFB1", "TGFBR1", "stroma->tumor", "TGF-beta 促 EMT"),
    ("LIF", "LIFR", "stroma->tumor", "LIF/LIFR"),
    ("GREM1", "BMP2", "stroma->tumor", "Gremlin-1 拮抗 BMP"),
    ("INHBA", "ACVR1B", "stroma->tumor", "Activin A/ALK4"),
    # --- 免疫 -> 基质 (反馈) ---
    ("IFNG", "IFNGR1", "immune->stroma", "IFNg 抗纤维化反馈"),
    ("TNF", "TNFRSF1A", "immune->stroma", "TNF/TNFR1"),
    ("IL1B", "IL1R1", "immune->stroma", "IL1b 促 CAF 炎性表型"),
    ("CCL2", "CCR2", "immune->stroma", "CCL2/CCR2 招募"),
]

LR_AXES = ["neuro->stroma", "neuro->immune", "stroma->immune",
           "stroma->tumor", "immune->stroma"]


# ------------------------------------------------------------------- 通路
PATHWAY_KEYWORDS = [
    "EXTRACELLULAR_MATRIX", "COLLAGEN", "FIBROBLAST", "FOCAL_ADHESION", "ADHESION",
    "TGF_BETA", "TGFB", "INTEGRIN", "PROTEOGLYCAN", "ELASTIN", "MATRISOME",
    "BASEMENT_MEMBRANE", "METALLOPROTEINASE", "MATRIX_METALLOPROTEINASE",
    "AXON", "NEURON", "NEUROGENESIS", "SYNAPSE", "NEUROTRANSMIT", "NEUROTROPH",
    "SCHWANN", "MYELIN", "CHOLINERGIC", "ADRENERGIC", "SEROTONIN", "GABA",
    "CHEMOKINE", "CYTOKINE", "INTERFERON", "ANTIGEN", "T_CELL", "B_CELL",
    "MACROPHAGE", "NEUTROPHIL", "COMPLEMENT", "IMMUNOGLOBULIN", "INFLAMMATORY",
    "INTERLEUKIN", "TOLL_LIKE", "IMMUNE", "LEUKOCYTE", "LYMPHOCYTE",
    "VASCULAR", "ANGIOGENESIS", "ENDOTHELIAL", "VASCULOGENESIS",
    "GLYCOLYSIS", "HYPOXIA", "OXIDATIVE_PHOSPHORYLATION", "CELL_CYCLE",
    "APOPTOSIS", "EPITHELIAL_MESENCHYMAL", "EMT", "WNT", "NOTCH", "MAPK",
    "PI3K", "NFKB", "JAK_STAT", "TNF", "MTOR", "HEDGEHOG", "HIPPO",
    "GPCR", "G_PROTEIN", "LIGAND", "RECEPTOR", "SIGNALING",
]
PATHWAY_COLLECTIONS = ("HALLMARK", "REACTOME", "KEGG", "WP", "GOBP", "GOCC", "PID")
PATHWAY_MIN, PATHWAY_MAX = 15, 400


def load_msigdb(keywords=PATHWAY_KEYWORDS, collections=PATHWAY_COLLECTIONS,
                min_size=PATHWAY_MIN, max_size=PATHWAY_MAX):
    """读 MSigDB gmt, 按关键词+集合大小筛选。返回 dict name->list[str]。

    名字格式为 "<COLLECTION>_<NAME>", 例如 "REACTOME_EXTRACELLULAR_MATRIX_ORGANIZATION"。
    """
    if not GMT.exists():
        raise SystemExit(f"缺少 MSigDB: {GMT}")
    pats = [re.compile(k, re.IGNORECASE) for k in keywords]
    out = {}
    with GMT.open(encoding="utf-8", errors="ignore") as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            name = f[0]
            coll = name.split("_")[0]
            if coll not in collections:
                continue
            if not any(p.search(name) for p in pats):
                continue
            genes = sorted({g.strip() for g in f[2:] if g.strip()})
            if not (min_size <= len(genes) <= max_size):
                continue
            out[name] = genes
    return out


def build_score_matrix(zg: np.ndarray, genes: np.ndarray, blocks: dict,
                       min_frac: float = 0.4, min_frac_by_kind: dict | None = None):
    """把 (细胞类型 / LR 对 / 通路) 三类读数统一成 spot x K 评分矩阵。

    zg    : n x G 基因级 z-score(跨 spot 标准化)
    genes : G 的 gene symbol
    blocks: {"celltype": {name: [genes]}, "lr": [(name, [lig, rec])], "pathway": {...}}

    评分 = 成员基因 z 的平均(线性) —— 线性保证系数域 MSR 可解析计算。
    成员缺失率 > (1-min_frac) 的读数丢弃。
    """
    idx = {g: i for i, g in enumerate(genes)}
    cols, names, kinds_meta = [], [], []
    kinds = [k for k in ("celltype", "lr", "pathway", "random") if k in blocks]
    kinds += [k for k in blocks if k not in kinds]
    for kind in kinds:
        mf = (min_frac_by_kind or {}).get(kind, min_frac)
        for nm, gl in blocks.get(kind, {}).items():
            rows = [idx[g] for g in gl if g in idx]
            if not rows or len(rows) / max(1, len(gl)) < mf:
                continue
            w = np.zeros(len(genes), dtype=np.float64)
            w[rows] = 1.0 / len(rows)
            cols.append(w)
            names.append(nm)
            kinds_meta.append(kind)
    if not cols:
        raise ValueError("没有任何可用读数")
    W = np.column_stack(cols)                 # G x K
    S = zg @ W                                # n x K
    meta = [{"name": n, "kind": k} for n, k in zip(names, kinds_meta)]
    return S, meta, W

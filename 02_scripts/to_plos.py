# -*- coding: utf-8 -*-
"""把 manuscript_main.md 转成 PLOS Computational Biology 投稿版 manuscript_PLOSCB.md。

变换:
1. PLOS 标题页 / Author Summary / 去 BIB 专属元素 (running title, Key Points, author profile)
2. 小节编号移除, 内部交叉引用改为名称引用
3. 图表引用改 PLOS S-notation (S13 Fig / S20 Table / S1 Text, Note 12)
4. 参考文献 66 条转 PLOS 格式 (作者 (年) 标题. 刊 卷: 页.), >6 作者截断
5. 句中粗体剥离 (保留段首结构粗体), 去 AI 口癖
"""
import re
from pathlib import Path

ROOT = Path(r"D:/嗜神经浸润")
SRC = ROOT / "04_manuscript/manuscript_main.md"
DST = ROOT / "04_manuscript/manuscript_PLOSCB.md"

md = SRC.read_text(encoding="utf-8")

# ---------------------------------------------------------------- 1. 摘要裁剪 (~300 词上限)
md = md.replace(
    "Whether that field can be quantified from spot-level transcriptomics, and whether it carries prognostic information beyond bulk RNA sequencing, is unknown.",
    "Whether that field can be quantified from spot-level transcriptomics is unknown.",
)

# ---------------------------------------------------------------- 2. 去 AI 口癖 (数字与统计量不动)
md = md.replace("they substantially measure desmoplasia", "they largely measure desmoplasia")
md = md.replace("The honest statement: at n = 43", "In plain terms: at n = 43")
md = md.replace("the whole range is the honest uncertainty", "the whole range is the construction uncertainty")
md = md.replace("the defensible statement is a nerve-proximal stromal niche with lymphocyte aggregation, not effector exclusion",
                "what the data support is a nerve-proximal stromal niche with lymphocyte aggregation, not effector exclusion")
md = md.replace("Fourth — decisively — the bulk association depends on signature construction:",
                "Fourth, the bulk association depends on signature construction:")
md = md.replace("Third — and most consequential for validity — adjacent Visium spots",
                "Third, and most consequential for validity, adjacent Visium spots")

# 3.3 独立成行的 "Three results emerged." (区别于 3.2 内嵌的)
md = re.sub(r"\nThree results emerged\.\n",
            "\nThree patterns emerged across the readout families.\n", md, count=1)

# ---------------------------------------------------------------- 3. 引用 S-notation / Fig 缩写
md = md.replace("Supplementary Figure S10C, S10D", "S10C and S10D Fig")
md = md.replace("Supplementary Figure S10E, S10F", "S10E and S10F Fig")
md = md.replace("Supplementary Figure S13, Note 13.28", "S13 Fig; S1 Text, Note 13.28")
md = md.replace("Supplementary Figure S14, Note 13.26", "S14 Fig; S1 Text, Note 13.26")
md = md.replace("Figures S11–S12", "S11 and S12 Figs")
md = re.sub(r"Supplementary Figure S(\d+)([A-F]?)", r"S\1\2 Fig", md)
md = re.sub(r"Supplementary Table S(\d+)", r"S\1 Table", md)
md = re.sub(r"Supplementary Note (\d+(?:\.\d+)?)", r"S1 Text, Note \1", md)
# 无前缀的 "(Table S15)" / "(Tables S4, S6, S7)"
md = re.sub(r"(?<!Supplementary )Tables (S\d)", r"\1s", md)   # 先占位, 防止 "Table S" 误吞
md = re.sub(r"(?<!Supplementary )Table (S\d)", r"\1", md)
md = md.replace("TablesS", "Tables S")
md = re.sub(r"\(S(\d+(?:–S\d+)?(?:, S\d+)*) (Tables?)\)", lambda m: "(%s %s)" % (m.group(1), m.group(2)), md)
# Note 残留 (如 "(Supplementary Notes 1–13" 在补充材料段, 该段整体替换, 此处兜底)
# 主图: Figure 5C, 5D -> Fig 5C and 5D
md = re.sub(r"Figure(s?) (\d)", r"Fig\1 \2", md)
md = re.sub(r"Fig (1[A-H]), (1[A-H](?:, 1[A-H])*)", lambda m: "Fig %s and %s" % (m.group(1), m.group(2).replace(",", " and")), md)

# ---------------------------------------------------------------- 4. 内部交叉引用改名称引用
md = md.replace("the partial-Spearman machinery of Methods 2.2–2.3", "the partial-Spearman machinery of Materials and methods")
md = md.replace("the gradient pipeline of Methods 2.2/2.4", "the gradient pipeline of Materials and methods")
md = md.replace("recomputed per sample as in Methods 2.5", "recomputed per sample as described in Materials and methods")
md = md.replace("against the MSR null (Methods 2.3;", "against the MSR null (Materials and methods;")
md = md.replace("re-derived the NPP as two modules (Methods 2.5)", "re-derived the NPP as two modules (Materials and methods)")
md = md.replace("the deconvolved T-cell gradient of Methods 2.8", "the deconvolved T-cell gradient of Materials and methods")
md = md.replace("(MSR plus size-matched random-gene-set baseline; Methods 2.9)", "(MSR plus size-matched random-gene-set baseline; see Materials and methods)")
md = md.replace("(Methods 2.3; Supplementary", "(Materials and methods; S")
md = md.replace("this adjustment is essential — otherwise housekeeping genes dominate the ranking because spot cellularity co-varies with nerve proximity; Supplementary Note 2)",
                "this adjustment is essential — otherwise housekeeping genes dominate the ranking because spot cellularity co-varies with nerve proximity; S1 Text, Note 2)")
md = md.replace("the stromal characterisation rests on Section 3.7", "the stromal characterisation rests on the three-line characterisation below")
md = md.replace("our result is positive on three calibrated lines (Section 3.7)", "our result is positive on three calibrated lines (Results, three-line characterisation)")
md = md.replace("the Section 3.7 myCAF anchoring", "the myCAF anchoring in Results")

# ---------------------------------------------------------------- 5. 小节标题去编号
md = re.sub(r"^## 1\. Introduction$", "## Introduction", md, flags=re.M)
md = re.sub(r"^## 2\. Materials and methods$", "## Materials and methods", md, flags=re.M)
md = re.sub(r"^## 5\. Conclusion$", "## Conclusion", md, flags=re.M)
md = re.sub(r"^### (\d+\.\d+) (.+)$", r"### \2", md, flags=re.M)

# ---------------------------------------------------------------- 6. 段首保留/句中剥离粗体
def strip_mid_bold(text: str) -> str:
    out_lines = []
    for line in text.split("\n"):
        if line.lstrip().startswith("#"):
            out_lines.append(line)
            continue
        # 段首粗体 (含编号列表首的 **X** 不在此列)
        m = re.match(r"^(\*\*[^*]+\*\*)", line)
        if m:
            head, rest = m.group(1), line[m.end():]
            rest = rest.replace("**", "")
            out_lines.append(head + rest)
        else:
            out_lines.append(line.replace("**", ""))
    return "\n".join(out_lines)

md = strip_mid_bold(md)

# ---------------------------------------------------------------- 7. 标题页 / Author Summary 替换
TITLE_PAGE = """# A calibrated spatial transcriptomics framework reveals a perineural stromal program in triple-negative breast cancer

**Authors:** [AUTHOR 1]¹, [AUTHOR 2]²*

**Affiliations:**

1. [DEPARTMENT], Hubei Provincial Maternal and Child Health Hospital (Hubei Women and Children's Hospital), Wuhan, Hubei, China
2. [DEPARTMENT], [INSTITUTION], [CITY], China

**Corresponding author:** [AUTHOR 2] ([corresponding@email]). ORCID: [0000-0000-0000-0000] ([AUTHOR 1]: [0000-0000-0000-0000])

**Keywords:** perineural invasion; spatial transcriptomics; triple-negative breast cancer; cancer-associated fibroblasts; Moran spectral randomisation; nerve–stroma axis

*One-sentence summary (submission form):* A two-stage calibration turns spot-level spatial transcriptomics into a distance field, revealing a ~150 µm myofibroblast-centred stromal niche around nerves in triple-negative breast cancer.

---

## Abstract
"""

# 原标题页 + Abstract 头整段替换
md = re.sub(
    r"^# .+?\n.*?\n---\n\n## Title page\n.*?\n---\n\n## Abstract\n",
    TITLE_PAGE,
    md, count=1, flags=re.S)

AUTHOR_SUMMARY = """
## Author Summary

Perineural invasion — the spread of cancer along nerves — is recorded in pathology reports as a simple yes-or-no observation, even though the tissue around a nerve changes gradually with distance. We asked whether this gradient can be measured from spatial transcriptomics data, in which thousands of spots of gene expression are sampled across a tumour section. Using 43 public triple-negative breast cancer samples, we located spots sitting next to nerves, measured how quickly each of about 10,000 genes fades with distance from those spots, and summarised the result as a 58-gene nerve-proximity program. A key technical point of our work is calibration: neighbouring spots are strongly correlated, so standard statistical tests produce many false discoveries — 47 times more in our data. We therefore added a second control, comparing every result against random gene sets of the same size. After calibration, the field around nerves in breast tumours turns out to be built mainly by myofibroblast-like cancer-associated fibroblasts and their matrix, over a radius of roughly 150 micrometres; immune cells are not measurably excluded. Nerve-proximity signatures used in the literature largely measure this generic scar-like tissue, which may explain why their prognostic value has been inconsistent. All data are public and the full analysis pipeline is open source.

---

"""

# Author Summary 插在 Introduction 前
md = md.replace("\n## Introduction\n", AUTHOR_SUMMARY + "## Introduction\n")

# 去掉 BIB Key Points 段
md = re.sub(r"\n## Key Points\n.*?\n---\n", "\n", md, count=1, flags=re.S)

# ---------------------------------------------------------------- 8. 背面事务段 (PLOS 措辞)
md = md.replace("""## Acknowledgements and Study funding

**Acknowledgements.** We thank the submitters of GSE210616, METABRIC and GSE12276 for making their data publicly available, and the developers of scanpy [60] and of the Moran spectral randomisation method [29] on which the calibration rests. [Acknowledgements for computational resources to be added.]

**Study funding.** [To be completed — list funders and grant numbers, or state "This work was supported by [funder] (grant [number])". If none: "This research received no specific grant from any funding agency in the public, commercial or not-for-profit sectors."]

**Competing interests.** The authors declare no competing interests.

**Ethics approval and consent to participate.** Not applicable. This study used exclusively publicly available, de-identified datasets released under the terms of NCBI GEO and cBioPortal; no human subjects, identifiable data or animals were involved.

---

## Data and code availability""", """## Acknowledgments

We thank the submitters of GSE210616, METABRIC and GSE12276 for making their data publicly available, and the developers of scanpy [60] and of the Moran spectral randomisation method [29] on which the calibration rests. [Acknowledgment for computational resources to be added.]

## Funding

[To be completed: "The author(s) received no specific funding for this work." or funder and grant number.]

## Competing interests

The authors have declared that no competing interests exist.

## Ethics statement

Not applicable. This study used exclusively publicly available, de-identified datasets released under the terms of NCBI GEO and cBioPortal; no human subjects, identifiable data, or animals were involved.

## Data and code availability""")

# ---------------------------------------------------------------- 9. Figure legends 标题 PLOS 化
md = md.replace("## Figure legends", "## Figure Legends")
md = re.sub(r"\*\*Figure (\d+)\.", r"**Fig \1.", md)

# ---------------------------------------------------------------- 10. Supporting information 段整体替换
SUPPORTING = """## Supporting information

**S1 Text.** Supplementary Notes 1–13: nerve-spot detection and threshold sensitivity; the housekeeping artefact and covariate adjustment; consensus-clustering robustness; cohort-candidate audit; cross-panel validation of nerve calls; bootstrap and E-value analysis; Moran spectral randomisation and the rank-transform correction; the three-threshold calibration of functional readouts; expression-derived diffusion pseudotime; module decomposition; external single-cell anchoring and deconvolution; the three-line characterisation of the perineural stroma; and external bulk-cohort extrapolation (TCGA-BRCA/TNBC, METABRIC, GSE163882, GSE25065, GSE91061, TCGA-PAAD, TCGA-HNSC, I-SPY2, GSE350571) including the size-bucket calibration erratum and the dialogue with the Cell 2026 sensory-neuron study.

**S1 Fig.** Nerve-spot detection threshold grid (see S15 Table).
**S2 Fig.** Spatial subtyping robustness (see S16 Table).
**S3 Fig.** METABRIC bootstrap hazard-ratio distributions and E-values (see S10 Table).
**S4 Fig.** Circularity audit of nerve calls.
**S5 Fig.** Spatial-autocorrelation calibration of the NPP (see S5 Table).
**S6 Fig.** Functional readouts of the perineural niche (see S4, S6, S7, S13, S17 Tables).
**S7 Fig.** Consensus clustering (see S16 Table).
**S8 Fig.** Organ-specific relapse in GSE12276 (see S14 Table).
**S9 Fig.** External single-cell anchoring of module genes and spot-level deconvolution (see S18, S19 Tables).
**S10 Fig.** Expression-derived trajectory and module decomposition (former main-figure content, now supporting the trajectory and decomposition results).
**S11 Fig.** External bulk-cohort extrapolation of the modules: fibrosis, immune, pCR/OS and ICB chains (see S22–S27 and S32 Tables).
**S12 Fig.** Immune compartments, mediation, pan-cancer replication, NAC-context DRFS, genotype bridge, spatial functional readouts and the I-SPY2 pCR split (see S28–S31 and S33–S35 Tables).
**S13 Fig.** Sensory-nerve/CGRP-axis proximity gradients and near-vs-far TLS functional ratios, the Cell 2026 dialogue (see S37 Table).
**S14 Fig.** Xenium single-cell-resolution nerve-proximity validation (see S36 Table).

**S1 Table.** Sample-level QC and nerve-spot summary.
**S2 Table.** Per-gene NPP statistics.
**S3 Table.** Signature-construction sensitivity.
**S4 Table.** Measurability and circularity audit of all 567 functional readouts.
**S5 Table.** MSR-calibrated per-gene statistics.
**S6 Table.** Pathway calibration (505 pathways).
**S7 Table.** Ligand–receptor audit (50 pairs).
**S8 Table.** Module definitions (39 + 60 genes) and per-sample spatial gradients.
**S9 Table.** Bulk Cox models in METABRIC.
**S10 Table.** Bootstrap hazard-ratio distributions and E-values.
**S11 Table.** Pseudotime per sample.
**S12 Table.** Prognostic incremental value.
**S13 Table.** Ring profiles of all functional readouts.
**S14 Table.** Organ-specific relapse analysis.
**S15 Table.** Nerve-spot threshold sensitivity.
**S16 Table.** Clustering robustness.
**S17 Table.** Size-matched random-gene-set null.
**S18 Table.** Single-cell cellular-origin assignment of module genes (GSE176078).
**S19 Table.** Spot-level deconvolution associations.
**S20 Table.** Calibrated gradients of the perineural-stroma program, CAF axes and myCAF/ECM markers (Line 1).
**S21 Table.** Immune-exclusion geometry and CAF-axis trajectory statistics (Lines 2–3).
**S22–S27 Tables.** Module migration, residualisation, immune attenuation, pCR and survival in external bulk cohorts.
**S28–S31 Tables.** Compartment enrichment and validation, mediation, METABRIC interaction and construction bootstrap, pan-cancer null-calibrated orthogonality.
**S32 Table.** ICB cross-indication exploration (GSE91061).
**S33 Table.** NAC-context DRFS (GSE25065) and the TCGA-BRCA genotype bridge.
**S34 Table.** Spatial TLS–CAF/collagen-maturation/hypoxia readouts, CD8 effector–chemokine decoupling, and Schwann/glial score projection (Notes 13.19–13.21).
**S35 Table.** I-SPY2 chemo-versus-IO pCR split (Note 13.24).
**S36 Table.** Xenium single-cell-resolution nerve-proximity validation (Note 13.26).
**S37 Table.** Sensory-nerve/CGRP-axis gradients and TLS functional ratios (Note 13.28).

Supporting figures (multi-panel composites) are in `03_results/figures/composed/`; supporting tables are provided as a single workbook, `03_results/tables/supplementary/Supplementary_Tables.xlsx`, one sheet per table; S1 Text is in `04_manuscript/supplementary/`."""

md = re.sub(
    r"\n## Supplementary material\n.*?\n\n---\n",
    "\n" + SUPPORTING + "\n\n---\n",
    md, count=1, flags=re.S)

# ---------------------------------------------------------------- 11. 内部提交清单 (PLOS 版)
md = re.sub(
    r"\n## Submission checklist \(internal — remove before submission\)\n.*$",
    """
## Submission checklist (internal — remove before submission)

- [ ] Fill placeholders: [AUTHOR 1], [AUTHOR 2], [DEPARTMENT], [INSTITUTION], [CITY], ORCIDs, corresponding email, funding sentence, GitHub owner, Zenodo DOI.
- [ ] Abstract ≤ 300 words; Author Summary plain-language (verify against current PLOS Comput Biol author guidelines).
- [ ] Alt text: verify each figure caption's alt text after any figure revision; enter at submission (required by PLOS).
- [ ] Figures: upload separately (TIFF/EPS/PDF per PLOS specs, `03_results/figures/composed/`); supporting information files uploaded as S-items; verify S-item numbering matches the Supporting information list.
- [ ] Suggest reviewers (PLOS form); provide repository link and Stride/ethics info (none required — public data).
- [ ] Confirm no related manuscripts under consideration elsewhere.
""", md, count=1, flags=re.S)

# ---------------------------------------------------------------- 12. 参考文献 PLOS 化
ref_start = md.find("\n## References\n")
body, refs = md[:ref_start], md[ref_start:]

def fmt_authors(authors: str) -> str:
    parts = [a.strip() for a in authors.split(", ") if a.strip()]
    if len(parts) > 6 and "et al" not in authors:
        return ", ".join(parts[:6]) + ", et al."
    return authors

def convert_ref(line: str) -> str:
    m = re.match(r"^(\d+)\. (.+?)\. (.+?)\. \*(.+?)\*\.? ?(\d{4});(.+?)\.?$", line)
    if m:
        n, auth, title, jour, year, vp = m.groups()
        auth = fmt_authors(auth)
        vp = vp.lstrip(";").strip()
        vp = re.sub(r"^(.+?):\s*(.+)$", r"\1: \2", vp)
        return f"{n}. {auth} ({year}) {title}. {jour} {vp}."
    # Cell 2026 带 doi 的在线优先
    m = re.match(r"^(\d+)\. (.+?)\. (.+?)\. \*(.+?)\*\. (\d{4})\. (doi:.+?)\.?$", line)
    if m:
        n, auth, title, jour, year, doi = m.groups()
        return f"{n}. {fmt_authors(auth)} ({year}) {title}. {jour}. {doi.replace('doi:', 'doi: ')}"
    # 书籍 A: "Authors. *Title*. ed. Publisher; year."
    m = re.match(r"^(\d+)\. (.+?)\. \*(.+?)\*\. (.*?)\.\s*(.+?);\s*(\d{4})\.?$", line)
    if m:
        n, auth, title, ed, pub, year = m.groups()
        return f"{n}. {fmt_authors(auth)} ({year}) {title}, {ed}. {pub}."
    # 书籍 B: "Authors. Title. City: Publisher; year."
    m = re.match(r"^(\d+)\. (.+?)\. (.+?)\. ([A-Z][^;]+?);\s*(\d{4})\.?$", line)
    if m:
        n, auth, title, pub, year = m.groups()
        return f"{n}. {fmt_authors(auth)} ({year}) {title}. {pub}."
    return line  # 未匹配, 保留原样待人工

out_lines, unfixed = [], []
for line in refs.split("\n"):
    if re.match(r"^\d+\. ", line):
        conv = convert_ref(line)
        if conv == line:
            unfixed.append(line)
        out_lines.append(conv)
    else:
        out_lines.append(line)
refs = "\n".join(out_lines)
md = body + refs

DST.write_text(md, encoding="utf-8")
print("written ->", DST)
print("未自动转换的参考文献 (%d):" % len(unfixed))
for u in unfixed:
    print("  ", u)

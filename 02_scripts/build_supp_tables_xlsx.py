# -*- coding: utf-8 -*-
"""把 03_results/tables/supplementary/ 下的 CSV 汇总为一个投稿用 xlsx 工作簿。

每个 Supplementary Table 一个 sheet, 首行为表题与来源说明, 数据从第 4 行开始。
输出: 03_results/tables/supplementary/Supplementary_Tables.xlsx
"""
import re
from pathlib import Path

import pandas as pd
from openpyxl.styles import Alignment, Font, PatternFill
from openpyxl.utils import get_column_letter

ROOT = Path(r"D:/嗜神经浸润")
SUPP = ROOT / "03_results/tables/supplementary"
OUT = SUPP / "Supplementary_Tables.xlsx"

TABLES = [
    ("S1", "Sample-level QC and nerve-spot summary (30 samples)",
     ["TableS1_sample_qc_and_nerve_spots.csv"]),
    ("S2", "Per-gene NPP statistics (10,419 genes)",
     ["TableS2_NPP_per_gene_statistics.csv"]),
    ("S3", "Signature-construction sensitivity: gene-set sizes, CAF/ECM correlations and hazard ratios for five constructions",
     ["TableS3_signature_construction_sensitivity.csv"]),
    ("S4", "Measurability and circularity audit of all 567 functional readouts "
           "(12 cell-type + 50 ligand-receptor + 505 pathway readouts; a further 5 marker panels could not be scored, giving 572 defined readouts in total)",
     ["TableS4_readout_measurability_and_circularity.csv"]),
    ("S5", "MSR-calibrated per-gene statistics (10,419 genes)",
     ["TableS5_MSR_per_gene_statistics.csv"]),
    ("S6", "Pathway calibration: 505 pathways under the two-stage (MSR + size-matched) calibration",
     ["TableS6_pathway_calibration.csv"]),
    ("S7", "Ligand-receptor audit: 50 pairs, measurability and calibration",
     ["TableS7_ligand_receptor_audit.csv"]),
    ("S8", "Neural-intrinsic (39 genes) and perineural-stroma (60 genes) module definitions",
     ["TableS8_module_definitions.csv"]),
    ("S8b", "Per-sample spatial gradients of the two modules",
     ["TableS8b_module_spatial_gradients.csv"]),
    ("S9", "Bulk Cox models in METABRIC (naive and residual NPP, axes, subtypes)",
     ["TableS9_bulk_cox_models_METABRIC.csv"]),
    ("S10", "Bootstrap hazard-ratio distributions and E-values (500 resamples x 4 models)",
     ["TableS10_bootstrap_and_evalue.csv"]),
    ("S11", "Pseudotime per sample: correlation with distance, root spot, covariate-adjusted",
     ["TableS11_pseudotime_per_sample.csv"]),
    ("S12", "Prognostic incremental value of spatial features over a clinical baseline",
     ["TableS12_prognostic_incremental_value.csv"]),
    ("S13", "Ring profiles of all functional readouts (5 rings x 567 readouts)",
     ["TableS13_readout_ring_profiles.csv"]),
    ("S14", "Organ-specific relapse analysis (GSE12276): Fine-Gray and Cox models",
     ["TableS14_organ_specific_relapse.csv"]),
    ("S15", "Nerve-spot threshold sensitivity (12 configurations)",
     ["TableS15_nerve_threshold_sensitivity.csv"]),
    ("S16", "Clustering robustness: jackknife, algorithm substitution, bootstrap",
     ["TableS16_clustering_robustness.csv"]),
    ("S17", "Size-matched random-gene-set null (510 sets x 30 samples)",
     ["TableS17_random_geneset_null.csv"]),
    ("S18", "Cellular origin of the 99 NPP module genes in GSE176078 scRNA (per cell type)",
     ["TableS18_module_cell_origin.csv"]),
    ("S19", "Spot-level deconvolution: gradients, nerve-spot contrast, module co-aggregation (primary and module-excluded panel)",
     ["TableS19_deconv_summary.csv"]),
    ("S20", "Line 1 - calibrated gradients of the perineural-stroma program, CAF axes and myCAF/ECM markers (MSR + size-matched baseline)",
     ["TableS20_axis_gradient_calibration.csv"]),
    ("S21", "Lines 2-3 - immune-exclusion geometry and CAF-axis trajectory statistics",
     ["TableS21_immune_geometry_caf_trajectory.csv"]),
    ("S22", "External bulk cohorts - module migration and neural-detection sensitivity (Note 13)",
     ["../bulk_extrapolation/stage1_module_correlations.csv",
      "../bulk_extrapolation/stage1_neural_detection_sensitivity.csv"]),
    ("S23", "External bulk cohorts - rank residualisation of modules against CAF/ECM and three-tier attenuation (Note 13)",
     ["../bulk_extrapolation/stage7_orthogonality.csv",
      "../bulk_extrapolation/stage7_attenuation.csv"]),
    ("S24", "External bulk cohorts - fibrosis-adjusted module residual vs immune readouts, BH-FDR (Note 13)",
     ["../bulk_extrapolation/stage7_residual_immune.csv"]),
    ("S25", "External bulk cohorts - pCR logistic models (raw and residualised) with random-gene-set baseline (Note 13)",
     ["../bulk_extrapolation/stage4_pcr_association.csv",
      "../bulk_extrapolation/stage7_residual_pcr.csv",
      "../bulk_extrapolation/stage6_random_baseline.csv"]),
    ("S26", "External bulk cohorts - incremental Cox models for eight matrix programs and single-gene decomposition (Note 13)",
     ["../bulk_extrapolation/stage8_cox_models.csv",
      "../bulk_extrapolation/stage8_gene_decomposition.csv"]),
    ("S27", "External bulk cohorts - robustness: bootstrap, random-gene-set baseline, KM (Note 13)",
     ["../bulk_extrapolation/stage9_robustness.csv",
      "../bulk_extrapolation/stage9_random_baseline.csv",
      "../bulk_extrapolation/stage9_km_refined.csv"]),
    ("S28", "External bulk cohorts - immune compartment enrichment: score validation on matched scRNA and composition vs modules (Note 13.9)",
     ["../bulk_extrapolation/stage10_score_validation.csv",
      "../bulk_extrapolation/stage10_deconv_composition.csv"]),
    ("S29", "External bulk cohorts - formal mediation of the stroma residual to CD8 cytotoxicity (Note 13.10)",
     ["../bulk_extrapolation/stage10_mediation.csv"]),
    ("S30", "METABRIC - module x CAF/ECM interaction, subtype-stratified OS, signature-construction bootstrap (Note 13.11)",
     ["../bulk_extrapolation/stage11_metabric_interaction.csv",
      "../bulk_extrapolation/stage11_metabric_subtype.csv",
      "../bulk_extrapolation/stage11_construction_bootstrap.csv"]),
    ("S31", "Pan-cancer - null-calibrated orthogonality and the fibrosis/TGF-beta/CD8/OS chain in seven cohorts (Note 13.12)",
     ["../bulk_extrapolation/stage12_orthogonality_random.csv",
      "../bulk_extrapolation/stage12_pancancer.csv"]),
    ("S32", "ICB cross-indication exploration - GSE91061 (Riaz 2017, melanoma, nivolumab +/- ipilimumab): module x RECIST response, random-set baseline, paired pre/on deltas (Note 13.14, exploratory)",
     ["../bulk_extrapolation/stage13_icb_logistic.csv",
      "../bulk_extrapolation/stage13_random_baseline.csv",
      "../bulk_extrapolation/stage13_paired_delta.csv"]),
    ("S33", "NAC-treatment-context DRFS in GSE25065 (module Cox per stratum with 510 random-set baseline) and the TCGA-BRCA genotype bridge (TMB, TP53/PIK3CA/BRCA1/BRCA2 status, TP53 interaction on OS) (Notes 13.16-13.17)",
     ["../bulk_extrapolation/stage15_drfs_cox.csv",
      "../bulk_extrapolation/stage15_random_baseline.csv",
      "../bulk_extrapolation/stage16_tmb_correlation.csv",
      "../bulk_extrapolation/stage16_mutation_status.csv",
      "../bulk_extrapolation/stage16_tp53_os.csv"]),
    ("S34", "Spatial functional readouts (TLS-CAF coupling, collagen maturation, hypoxia, size-bucketed null), CD8 effector/chemokine decoupling in bulk, and Schwann/glial score projection (Notes 13.19-13.21)",
     ["../spatial_tlscaf/stage17_gradient_meta.csv",
      "../spatial_tlscaf/stage17_coupling_per_sample.csv",
      "../spatial_tlscaf/stage17_coupling_meta.csv",
      "../bulk_extrapolation/stage18_cd8_decoupling.csv",
      "../bulk_extrapolation/stage19_schwann.csv"]),
    ("S35", "I-SPY2 chemo versus chemo+pembrolizumab pCR in TNBC: module logistic per arm, stroma x pembrolizumab interaction, random-set baseline (Note 13.24, exploratory cross-regimen)",
     ["../bulk_extrapolation/stage20_ispy2_pcr.csv"]),
    ("S36", "Xenium Prime 5K single-cell-resolution nerve-proximity validation (GSE350571, 6 FFPE breast tumours): axonal-anchored nerve kernels, cell-type proximal/distal composition with torus-shift permutation, fibroblast module scores proximal vs distal (Note 13.26, direction-level)",
     ["../bulk_extrapolation/stage14_xenium_nerve_summary.csv",
      "../bulk_extrapolation/stage14_xenium_fb_scores.csv"]),
    ("S37", "Sensory-nerve and CGRP-axis proximity gradients with size-matched random-set calibration, and near-vs-far TLS functional ratios (cytotoxicity, exhaustion, TGF-beta response) in TLS+ spots (Note 13.28, dialogue with Cell 2026)",
     ["../spatial_tlscaf/stage21_axis_gradient_meta.corrected.csv",
      "../spatial_tlscaf/stage21_tls_function_meta.csv"]),
]

HDR_FILL = PatternFill("solid", fgColor="DDEBF7")
TITLE_FONT = Font(bold=True, size=11)
SRC_FONT = Font(italic=True, size=9, color="555555")
HDR_FONT = Font(bold=True, size=10)


def sheet_for(wb, code, title, files):
    name = re.sub(r"[\\/*?:\[\]]", "_", f"Table_{code}")[:31]
    ws = wb.create_sheet(name)
    ws.cell(row=1, column=1, value=f"Supplementary Table {code}. {title}").font = TITLE_FONT
    ws.cell(row=2, column=1, value="Source: " + ", ".join(files)).font = SRC_FONT
    start = 4
    r = start
    for f in files:
        p = SUPP / f
        if not p.exists():
            ws.cell(row=r, column=1, value=f"[missing file: {f}]")
            r += 2
            continue
        df = pd.read_csv(p, nrows=20000)
        for j, col in enumerate(df.columns, start=1):
            c = ws.cell(row=r, column=j, value=str(col))
            c.font = HDR_FONT
            c.fill = HDR_FILL
            c.alignment = Alignment(wrap_text=True, vertical="center")
        for i, row in enumerate(df.itertuples(index=False), start=1):
            for j, v in enumerate(row, start=1):
                if isinstance(v, float) and (v != v):
                    ws.cell(row=r + i, column=j, value="")
                elif isinstance(v, (int, float)):
                    ws.cell(row=r + i, column=j, value=v)
                else:
                    ws.cell(row=r + i, column=j, value=str(v))
        # 列宽: 取前 200 行样本估算
        for j, col in enumerate(df.columns, start=1):
            try:
                sample = df[col].astype(str).head(200)
                w = max(len(str(col)) * 0.8, sample.str.len().max() if len(sample) else 8) + 2
            except Exception:
                w = 14
            ws.column_dimensions[get_column_letter(j)].width = min(max(w, 9), 42)
        ws.freeze_panes = ws.cell(row=r + 1, column=1)
        r += len(df) + 3
    return ws


def main():
    from openpyxl import Workbook
    wb = Workbook()
    ws = wb.active
    ws.title = "README"
    ws.column_dimensions["A"].width = 20
    ws.column_dimensions["B"].width = 110
    ws.cell(row=1, column=1, value="Supplementary Tables").font = Font(bold=True, size=13)
    ws.cell(row=2, column=1, value="Manuscript: A calibrated spatial transcriptomics framework reveals a "
                                   "perineural stromal program in triple-negative breast cancer")
    ws.cell(row=3, column=1, value="Note").font = Font(bold=True)
    ws.cell(row=3, column=2,
            value="Each sheet is one supplementary table. Numeric cells are rounded as produced by the pipeline; "
                  "readout counts follow the manuscript convention: 567 scoreable readouts = 12 cell-type panels "
                  "+ 50 ligand-receptor pairs + 505 pathways (a further 5 marker panels could not be scored).")
    for i, (code, title, files) in enumerate(TABLES, start=5):
        ws.cell(row=i, column=1, value=f"Table {code}").font = Font(bold=True)
        ws.cell(row=i, column=2, value=title)
    for code, title, files in TABLES:
        sheet_for(wb, code, title, files)
        print("  sheet", code, "ok")
    wb.save(OUT)
    print("saved ->", OUT)


if __name__ == "__main__":
    main()

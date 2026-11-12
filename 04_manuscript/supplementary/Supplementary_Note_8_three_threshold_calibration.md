# Supplementary Note 8 — The three-threshold calibration of functional readouts

## 8.1 Why a third threshold is needed

Readouts are scored as the mean of member-gene *z*-scores. Both the readout field and the distance field are spatially smooth, so the correlation between them increases with the number of genes averaged. We measured this directly with **510 random gene sets** (30 at each of 17 sizes from 5 to 400 genes) drawn from the same 10,503-gene expressed universe, scored and tested identically to the real readouts: the baseline |mean ρ| rises from 0.022 at 5 genes to 0.102 at 400 genes (0.030 at 10, 0.045 at 30, 0.069 at 100).

Under the MSR null alone, 96.7% of these random sets reach *P* < 0.05, and **88.2%** additionally reach |ρ| ≥ 0.03. A calibration that admits 88% of random gene sets is not a calibration. This is the reason for the third threshold.

Source table: `TableS17_random_geneset_null.csv`.

## 8.2 The three thresholds

1. **MSR calibration** (Note 7): the readout's own partial Spearman statistic standardised against its own MSR null, combined across samples by Stouffer's method and by a between-sample *t* test with BH correction.
2. **Measurability**: assessable in ≥10 of 43 samples. Readouts whose member genes are too sparsely detected produce unstable effect sizes and are excluded from inference.
3. **Non-circularity**: no membership overlap with the nerve-definition marker panel. Any readout containing those markers correlates with nerve distance by construction.

After thresholds 2 and 3, a readout must additionally **exceed its size-matched random baseline** on both effect size and between-sample consistency. Empirical *P* values come from the size-standardised, pooled residual distribution of the 510 random sets (resolution 1/511) with BH correction.

## 8.3 Measurability audit

Of 17 cell-type panels, five (*Adipocyte, B_lymphocyte, Mast, NK, Neuronal_Schwann*) yielded no scoreable readout because their markers were too sparsely detected (Adipocyte produced a gradient in a single sample). Among the twelve that returned a statistic, *DC* was assessable in 9 samples (<10, descriptive ring profiles only), while *CD8_T* (15 samples), *Neutrophil* (16) and *T_lymphocyte* (17) cleared the measurability threshold. Of 50 ligand–receptor pairs, 34 were measurable in ≥10 samples. Of 507 pathways, 475 were measurable.

Source table: `TableS4_readout_measurability_and_circularity.csv`.

## 8.4 Ligand–receptor measurability in detail

The canonical neurotrophic and axon-guidance axes could not be assessed: *NGF–NTRK1*, *NGF–NGFR*, *SEMA3A–NRP1* and *BDNF–NTRK2* had too few detected spots for a stable score; *S100B–AGER* was measurable in 7 samples and *CALCA–RAMP1* in 2. The large apparent effects in these pairs (ρ up to −0.27) are exactly the small-*n* artefact the measurability threshold exists to remove. Among measurable pairs the strongest gradients were stromal: *COL1A1–ITGA2*, *PDGFB–PDGFRB*, *TNC–TLR4*, *TGFB1–TGFBR1/2* (|ρ| 0.03–0.05, several reaching MSR-calibrated FDR < 0.05) but none exceeded the random baseline on both criteria (empirical FDR ≥ 0.13 on between-sample consistency).

Source table: `TableS7_ligand_receptor_audit.csv`.

## 8.5 What passes at the full cohort (n = 43)

**No pathway passed all thresholds jointly.** Eleven pathway readouts exceeded the size-matched baseline on effect size alone (empirical FDR = 0.024, the 1/511 resolution floor): ECM interactions, EPH/ephrin and EPHB forward signalling, β-catenin-independent Wnt, apoptosis, NOTCH4 (and its negative regulation), Hedgehog ligand biogenesis and on-state, semaphorin-3A axon repulsion, NCAM neurite outgrowth, MHC class I folding/assembly with cross-presentation, and IFN-γ signalling (|ρ| 0.061–0.119 against baselines of 0.032–0.088). On between-sample consistency the same readouts sat at empirical FDR 0.07–0.08, just above the joint threshold. No cell-type and no ligand–receptor readout passed on either criterion jointly.

The reason the pathway layer tightens with cohort size is visible in the random null itself: the pooled size-standardised |Z| distribution of the 510 random sets sharpens as *n* grows (median |Z| 4.38 at n = 30 → 5.43 at n = 43), so the empirical-FDR bar for "cross-sample consistency" rises mechanically with *n*. We therefore read the 0/507 joint result not as "no biology" but as the honest statement that at n = 43 no pathway *category* is distinguishable from a size-matched random gene set at joint FDR < 0.05; the stromal characterisation of the paper rests on the module- and axis-level analyses of Section 3.7, calibrated against the same machinery. Effect sizes are reported as enrichment of *machinery*, not of activity.

Source table: `TableS6_pathway_calibration.csv`. Figure: Supplementary Figure S6, Figure 3.

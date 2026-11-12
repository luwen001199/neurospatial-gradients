# Supplementary Note 13 — External bulk-cohort extrapolation: fibrosis, immune and chemotherapy-response chains

This note documents the external validation of the spatially derived modules in five public bulk cohorts, executed after the discovery analyses were frozen. It answers one question: **what survives when the perineural-stroma module is carried out of the spatial context and confronted with fibrosis quantification, immune readouts, pathological complete response (pCR) and overall survival?** The headline results are stated in Discussion Section 4.3 of the manuscript; this note gives the full numbers, the cohort audit and the failure modes.

All code is in `02_scripts/run_bulk_stage1.py` (module migration), `run_bulk_stage2.py` (immune associations), `run_bulk_stage4.py` (pCR), `run_bulk_stage7_orthogonal.py` (residualisation and three-tier attenuation), `run_bulk_stage8_fibrosis.py` (matrix quantity vs quality), `run_bulk_stage9_robustness.py` (bootstrap, random-gene-set baseline, KM); figure `02_scripts/make_figure_supp_s11.py` (Supplementary Figure S11). Result tables: `03_results/tables/bulk_extrapolation/`.

## N13.1 Cohort audit

| cohort | platform | n | TNBC subset | labels used |
|---|---|---|---|---|
| TCGA-BRCA | RNA-seq (RSEM), cBioPortal per-gene API | 1,100 | 122 (ER/PR/HER2 negative) | OS (761 evaluable, 96 events), PAM50, age, node |
| GSE163882 | AmpliSeq targeted RNA-seq, TPM, NAC baseline biopsies | 222 | 90 (38 pCR) | pCR, age, receptor status |
| GSE25065 | Affymetrix U133A, HER2-negative NAC | 198 (182 labelled) | 66 | pCR (42 events) |

Module scores are mean-*z* of member genes (`bulk_lib.score_mean` / `score_diff`); the 58-gene NPP, 39-gene neural-intrinsic and 60-gene perineural-stroma definitions are those of Notes 10/12, unchanged. The curated CAF/ECM reference is the 14-gene panel of Methods 2.x (COL1A1/2, COL3A1, COL5A1, FN1, POSTN, ACTA2, TAGLN, FAP, PDGFRB, SPARC, THBS2, CTHRC1, LOX). Probe-to-symbol mapping for GSE25065 borrows the GPL570 annotation (94.9% of U133A probes mapped); TCGA expression was fetched gene-wise from cBioPortal (`brca_tcga_rna_seq_v2_mrna`) rather than from a single monolithic matrix, with 219 of 226 target genes recovered.

## N13.2 Stage 1 — migration succeeds, orthogonality does not

Spearman correlations (bootstrap 95% CI, 1,000 replicates; full table `stage1_module_correlations.csv`):

| cohort | n | stroma ~ CAF/ECM | neural ~ CAF/ECM |
|---|---|---|---|
| TCGA-BRCA | 1,100 | 0.973 | 0.198 |
| TCGA-TNBC | 122 | 0.964 | 0.252 |
| GSE163882 | 222 | 0.973 | 0.406 |
| GSE163882-TNBC | 90 | 0.970 | 0.569 |
| GSE25065 | 198 | 0.972 | 0.124 |
| *METABRIC (manuscript, Note 10)* | *1,976* | *0.975* | ***0.040*** |

The stroma module ≡ desmoplasia result replicates in every cohort (ρ = 0.96–0.98). The neural-intrinsic module's **absolute** correlation with CAF/ECM does not replicate (0.04 in METABRIC, 0.12–0.57 elsewhere), and the inflation tracks gene detection: stratifying GSE163882 by detection rate, the high-detection subset shows the *strongest* correlation (≥ 50%: ρ = 0.488; < 50%: ρ = 0.117; `stage1_neural_detection_sensitivity.csv`).

An absolute ρ cannot be interpreted against zero, because any co-expressed gene set carries some correlation with a stromal score. We therefore re-expressed the claim against the correct null (Stage 12, N13.12): 1,000 size-matched random gene sets drawn from the same panel, in each cohort. Against that null the neural module is **more** orthogonal to CAF/ECM than chance in five of seven cohorts (z = −2.0 to −11.1; empirical *P* = 0.044 to < 0.001) and indistinguishable from random in two (GSE163882 z = +1.5, *P* = 0.14; GSE25065 z = +0.4, *P* = 0.65), while the stroma module sits 6.5–10.0 SD **above** the same null in all seven. The defensible statement, adopted in Results 3.5, is therefore: the neural module is orthogonal to desmoplasia *relative to a size-matched null in most cohorts*, with a platform-dependent absolute magnitude; the stroma module is desmoplasia.

## N13.3 Stage 7 — residualisation: what remains after fibrosis is partialled out

Because the stroma module shares ≥ 93% of its rank variance with CAF/ECM (`stage7_orthogonality.csv`: R² 0.930–0.946 across cohorts; LRRC15⁺ myCAF panel 0.935–0.961; SPARC/POSTN core 0.818–0.864; collagen maturation 0.51–0.68), all downstream tests use (i) raw scores, (ii) partial Spearman controlling CAF/ECM, (iii) partial Spearman controlling CAF/ECM plus an epithelial purity proxy (mean *z* of EPCAM, KRT8/18/19/7, CDH1, MUC1), and (iv) rank-space residuals of the module against CAF/ECM (`stage7_attenuation.csv`).

Perineural-stroma module vs immune readouts, ρ at three adjustment tiers (median across five cohorts; full table `stage7_attenuation.csv`):

| target | raw | + CAF/ECM | + CAF/ECM + purity |
|---|---|---|---|
| TGF-β core | +0.67 | +0.37 | +0.37 |
| M2 macrophage | +0.43 | +0.27 | +0.25 |
| SPP1⁺ macrophage | +0.43 | +0.27 | +0.25 |
| MHC class II | +0.26 | +0.23 | +0.23 |
| Treg | +0.16 | +0.20 | +0.15 |
| TLS | +0.10 | +0.21 | +0.19 |
| IFN-γ | +0.14 | +0.20 | +0.16 |
| CD8 cytotoxicity | +0.09 | +0.22 | +0.19 |

Findings: (1) the purity proxy is not a confounder — adding it changes nothing beyond CAF/ECM; (2) TGF-β signalling retains 50–80% of its raw association after full fibrosis adjustment and is the strongest surviving companion in every cohort (residual ρ 0.29–0.57, BH *q* < 1.5 × 10⁻³ in all five; 36 of 40 module–target pairs at *q* < 0.05, `stage7_residual_immune.csv`); (3) **CD8 cytotoxicity is not negatively associated with the stroma module at any adjustment tier** — the weak calibrated CD8 decline near nerves observed spatially (Note 12, Z = −3.14) is a local geometric phenomenon that leaves no imprint on bulk CD8 abundance, which instead correlates *positively* with the fibrosis-adjusted stroma residual.

## N13.4 Stage 4 — pCR under neoadjuvant chemotherapy

Logistic regression, OR per SD, adjusted for age where available; 5-fold cross-validated AUC (`stage4_pcr_association.csv`, `stage7_residual_pcr.csv`). Neither the NPP, the stroma module, the neural module, CAF/ECM nor the LRRC15⁺ myCAF panel predicts pCR in any cohort (all raw OR 0.74–1.07, *P* > 0.10; AUC 0.42–0.63). The 510-set size-matched random baseline (Stage 6, `stage6_random_baseline.csv`) confirms the null: observed AUCs 0.45–0.51 sit at or below the random median (0.51–0.57), empirical *P* 0.42–0.98. The only replicated positive predictor is CD8 cytotoxicity (GSE163882 OR 1.69 [1.03, 2.80], *P* = 0.037, AUC 0.67; TNBC subset 1.73 [1.04, 2.86], *P* = 0.033), consistent with the literature.

Residualised (fibrosis-adjusted) scores do not rescue the prediction: the perineural-stroma residual is directionally consistent across all three cohorts (OR 1.32 / 1.54 / 1.16, all *P* > 0.05) and the SPARC/POSTN-core residual is significant in GSE163882 (1.37 [1.02, 1.85], *P* = 0.038; TNBC 1.68 [1.05, 2.69], *P* = 0.030) but reverses in GSE25065 (0.80, *P* = 0.22). We report this as failure to replicate rather than as a discovery (Supplementary Figure S11C).

## N13.5 Stage 8 — matrix quantity versus matrix quality (overall survival)

Eight matrix programs (collagen core, collagen maturation, SPARC/POSTN core, mechanotransduction, EMT core, TGF-β program, CAF/ECM total, perineural-stroma module) were tested against OS in TCGA-BRCA (n = 761, 96 events; Cox, Breslow ties, Newton–Raphson implemented in `bulk_lib.cox_breslow` after statsmodels PHReg produced degenerate standard errors at HR ≈ 1) with incremental adjustment M0 (univariate) → M1 (+ age, node) → M2 (+ Basal) → M3 (+ CAF/ECM total) (`stage8_cox_models.csv`):

| program | M0 raw HR | M3 residual HR | P (M3) |
|---|---|---|---|
| SPARC/POSTN core | 1.117 | **1.297 (1.073, 1.567)** | 0.0072 |
| CAF/ECM total | 1.027 | 0.976 | 0.91 |
| Collagen core | 0.984 | 0.897 | 0.28 |
| Collagen maturation | 1.028 | 1.080 | 0.52 |
| Mechanotransduction | 1.002 | 0.931 | 0.55 |
| EMT core | 1.046 | 1.091 | 0.43 |
| TGF-β program | 1.114 | 1.111 | 0.35 |
| Perineural-stroma module | 1.030 | 1.036 | 0.70 |

In TCGA-TNBC (117 evaluable, 19 events) the SPARC/POSTN residual is directionally consistent (HR 1.51–1.53 after adjustment) but not significant (*P* = 0.052–0.053) and should be read only as concordant. Single-gene decomposition shows the signal is program-level, not gene-level: no individual member (SPARC, POSTN, THBS2, TNC, CTHRC1, COL1A1, LRRC15, LOX, PLOD2) reaches significance after residualisation and adjustment (`stage8_gene_decomposition.csv`).

## N13.6 Stage 9 — robustness of the SPARC/POSTN residual signal

(`stage9_robustness.csv`, `stage9_random_baseline.csv`)

1. **Incremental adjustment.** HR 1.249 (M0) → 1.301 (M1: age/node) → 1.301 (M2: + Basal) → 1.297 (M3: + CAF/ECM total) → 1.281 (+ MKI67) — the association is not explained by proliferation, subtype or stromal quantity.
2. **Bootstrap.** 1,000 resamples: median HR 1.259, 95% CI 1.061–1.511 (TCGA-BRCA); 1.277 (0.809–2.108) in TCGA-TNBC.
3. **Random-gene-set baseline.** 510 size-matched random 8-gene sets, residualised identically: random median HR 0.92 (95% range 0.75–1.10) in TCGA-BRCA; the observed 1.25 exceeds |log HR| of 92% of random sets (empirical *P* = 0.080). In TCGA-TNBC the empirical *P* is 0.38. The signal therefore does **not** clear the dual-calibration standard applied to the spatial analyses, and the manuscript labels it exploratory.
4. **Grouped survival.** Upper-tertile vs lower-two-thirds KM: log-rank *P* = 0.37 (TCGA-BRCA), 5-year OS 77.0% vs 81.3%; in TCGA-TNBC 50.4% vs 73.3% (*P* = 0.19, 19 events).

## N13.7 Interpretation and limits

- The bulk shadow of the perineural-stroma module is **desmoplasia quantity plus a TGF-β/M2-polarised myCAF program**; it is not T-cell exclusion.
- Matrix *quantity* carries no prognostic information in any program tested; the only hint of prognostic structure lies in the SPARC/POSTN core *relative to* quantity, and it fails the random-set calibration.
- These cohorts contain no spatial nerve annotation, so this note tests the *bulk footprint* of the modules, not the perineural field itself; no claim in the manuscript's Results 3.1–3.7 depends on it.
- The stage-1 orthogonality failure modifies one sentence of the original manuscript narrative (neural-intrinsic ρ = 0.040 is now stated as METABRIC-specific), which we consider a self-correction enabled by the audit rather than a change of conclusion.

## N13.8 Data and code availability for this note

TCGA-BRCA via cBioPortal API (`brca_tcga`, `brca_tcga_rna_seq_v2_mrna`) and the Xena BRCA clinical matrix; GSE163882 via GEO series matrix and supplementary TPM file; GSE25065 via GEO series matrix with GPL570 annotation. All fetched artifacts cached under `01_data/raw/{tcga_brca,gse163882,gse25065}`; fetch script `02_scripts/_fetch_tcga_cbioportal.py`. R was not used; all statistics are NumPy/SciPy/statsmodels/scikit-learn implementations.


## N13.9 Stage 10 — immune compartments: why deconvolution was abandoned, and what replaced it

(`stage10_score_validation.csv`, `stage10_deconv_composition.csv`)

**Deconvolution is not identifiable from this panel.** We built a 10-compartment signature matrix (CD8 T, CD4 T, NK, B/plasma, macrophage, monocyte, DC, endothelial, fibroblast, epithelial) directly from GSE176078 single-cell reference (100,064 cells, `build_deconv_signature.py`, CP10K means per `celltype_minor`), which is the same construction as a CIBERSORTx custom signature matrix. NNLS on the 256 genes available in the bulk panel failed under three different scalings: (i) linear CP10K with intercept — all immune weights zero, fibroblast 0.56 / epithelial 0.32; (ii) square-root transform — identical degeneracy; (iii) per-gene *z*-standardisation — weights collapse to a near-uniform distribution (R² = 0.003). The reason is structural: the panel is a targeted ECM/neural panel in which 57 of 256 genes are fibroblast-specific and 5–14 genes carry each immune lineage, so the design matrix is rank-deficient for the mixture problem. CIBERSORTx itself is not runnable here (token + web upload) and xCell requires R.

**What replaced it.** xCell/Danaher-style enrichment scores: for each compartment we kept genes with CP10K ≥ 0.5 in that compartment and ≥ 2-fold over the next-highest compartment (117 genes total), then scored bulk samples as the mean *z* of those genes. This does not estimate proportions, but its ability to read composition is testable: in the 26 matched single-cell samples (pseudobulk), each compartment score correlates with the true cell fraction at ρ = 0.79–0.92 (DC 0.79, NK 0.81, fibroblast 0.83, CD4 0.84, CD8 0.84, epithelial 0.88, macrophage 0.91, endothelial 0.92, B/plasma 0.92; all *P* < 2 × 10⁻⁶). The scores therefore track composition closely enough for association testing, and we label them enrichment, not proportion.

**Composition of the perineural-stroma gradient (TCGA-BRCA, n = 1,093; TNBC n = 122).** Against the raw stroma module: fibroblast ρ = +0.989 (Δ high vs low = +2.12 SD), endothelial +0.575, macrophage +0.404, CD4 T +0.197, NK +0.157, epithelial +0.154, DC +0.131, B/plasma +0.118 (all BH *q* < 10⁻³). The module is thus not "stroma replacing immune cells"; it co-varies positively with every leukocyte compartment except that the epithelial fraction does not fall. This is the compositional form of the same message as N13.3: the field is a leukocyte-rich, matrix-rich stroma.

After residualising the module against CAF/ECM the leukocyte associations **strengthen rather than invert** (endothelial +0.460, DC +0.369, macrophage +0.355, CD4 +0.331, NK +0.328, CD8 +0.281, B +0.273; fibroblast falls to +0.132), i.e. the fibrosis-independent part of the module is the part most associated with immune infiltration. TNBC reproduces this ordering (DC +0.442, macrophage +0.392, endothelial +0.368, CD8 +0.312).

## N13.10 Stage 10 — formal mediation: where the CD8 association comes from

(`stage10_mediation.csv`; 2,000 bootstrap replicates, product-of-coefficients, all variables rank-normal-scored)

Exposure = perineural-stroma residual (against CAF/ECM); outcome = CD8 cytotoxicity; total effect *c* = +0.306 (TCGA-BRCA) and +0.361 (TCGA-TNBC) — note the sign: **positive**.

| mediator | *a* | *b* | *ab* [95% CI] | proportion mediated |
|---|---|---|---|---|
| TGF-β program | +0.315 | +0.254 | +0.080 [0.060, 0.101] | 26% |
| macrophage enrichment | +0.347 | +0.564 | +0.195 [0.164, 0.230] | 64% |
| M2 marker score | +0.340 | +0.631 | +0.215 [0.179, 0.251] | 70% |
| Treg marker score | +0.310 | +0.828 | +0.257 [0.214, 0.298] | 84% |
| CD8 T enrichment | +0.289 | +0.930 | +0.269 [0.218, 0.315] | 88% |

TNBC gives the same ordering with wider intervals (macrophage 72%, M2 67%, Treg 61%, CD8 T 75%) and a non-significant TGF-β path (9%, CI −0.02 to 0.10). The reverse-direction model (CD8 → TGF-β → residual) is also non-zero (+0.079 [0.060, 0.100]), so the mediation is not directional evidence — bulk data cannot order these variables in time. What the analysis does establish is arithmetic: 60–88% of the residual's CD8 association is shared with myeloid/lymphoid infiltration, leaving a direct path *c*′ = +0.04 to +0.11. There is no suppression to mediate.

## N13.11 Stage 11 — METABRIC: interaction, subtype strata, and how unstable the signature really is

(`stage11_metabric_interaction.csv`, `stage11_metabric_subtype.csv`, `stage11_construction_bootstrap.csv`)

1. **Interaction.** Cox with module × CAF/ECM (both per SD, adjusted for age; n = 1,979, 1,143 events): interaction HR 0.966 (*P* = 0.12) for the stroma module, 1.051 (*P* = 0.32) for its residual, 1.004 (*P* = 0.88) for the NPP, 1.027 (*P* = 0.27) for the neural module. The module's (null) survival association does not depend on desmoplasia level — the two are the same variable, so there is no separable effect to modulate.
2. **Subtype strata.** Residual score vs OS: Basal 0.95 [0.79, 1.15] (*P* = 0.61), Her2 0.88 [0.74, 1.05] (*P* = 0.16), LumA 0.91 [0.83, 1.00] (*P* = 0.061), LumB 0.95 [0.85, 1.06] (*P* = 0.35), Normal 0.91 [0.73, 1.14] (*P* = 0.40), claudin-low 0.91 [0.76, 1.10] (*P* = 0.32). The manuscript's "basal-like shows no residual association" is thus quantified and extended: no subtype shows one.
3. **Construction instability, measured properly.** Rebuilding the module from a random 80% of its genes (1,000 draws) gives ρ with CAF/ECM of 0.972 [0.955, 0.982] for the perineural-stroma module — i.e. the "signature ≡ desmoplasia" result is *not* an artefact of which 60 genes were chosen; any 48 of them reproduce it. For the neural module the 80% rebuilds give 0.036 [−0.050, 0.119] against a random-gene-set median of 0.352: here the instability is in the opposite direction — the module is stable and it is *below* chance, which is the null-calibrated orthogonality statement of N13.2.

## N13.12 Stage 12 — pan-cancer replication and null-calibrated orthogonality

(`stage12_orthogonality_random.csv`, `stage12_pancancer.csv`; 1,000 size-matched random sets per module per cohort, drawn from the same targeted panel)

**Null calibration** (ρ observed vs random-set distribution):

| cohort | n | stroma ρ (z) | neural ρ (z, *P*) | CD8 ρ (z) |
|---|---|---|---|---|
| TCGA-BRCA | 1,093 | 0.972 (+9.4) | 0.189 (−5.1, < 0.001) | +0.042 (−2.3) |
| TCGA-TNBC | 122 | 0.964 (+7.1) | 0.252 (−3.9, < 0.001) | −0.139 (−2.7) |
| TCGA-PAAD | 178 | 0.978 (+7.4) | 0.136 (−11.1, < 0.001) | +0.261 (−1.6) |
| TCGA-HNSC | 520 | 0.956 (+8.5) | 0.338 (−3.5, < 0.001) | −0.074 (−3.3) |
| METABRIC | 1,980 | 0.976 (+10.0) | 0.041 (−2.0, 0.044) | −0.018 (−1.1) |
| GSE163882 | 222 | 0.973 (+8.5) | 0.406 (+1.5, 0.14) | +0.018 (−1.1) |
| GSE25065 | 198 | 0.972 (+6.5) | 0.124 (+0.4, 0.65) | −0.071 (−0.6) |

Note that the null median itself varies 0.03–0.77 across cohorts and modules, because the panel composition and platform differ; this is precisely why a fixed ρ threshold cannot be used and why the manuscript reports cohort-internal z scores.

**Pan-cancer chain.** The three-part structure reproduces in tumours where perineural invasion is clinically common:

| cohort | n | CAF/ECM explains (rank R²) | residual ~ TGF-β | residual ~ CD8 | residual vs OS (HR per SD) |
|---|---|---|---|---|---|
| TCGA-BRCA | 1,093 | 95.8% | +0.456 | +0.302 | 0.95 |
| TCGA-TNBC | 122 | 95.4% | +0.379 | +0.373 | 0.70 |
| GSE163882 | 222 | 96.2% | +0.576 | +0.394 | — |
| GSE25065 | 198 | 94.6% | +0.254 | +0.120 | — |
| METABRIC | 1,980 | 95.4% | n/a | +0.378 | 0.92 |
| TCGA-PAAD | 178 | 96.0% | +0.364 | +0.382 | 1.02 |
| TCGA-HNSC | 520 | 96.2% | +0.431 | +0.116 | 0.96 |

So "the bulk proxy of the perineural field is desmoplasia plus a TGF-β-rich, immune-infiltrated myCAF stroma, with no prognostic content" holds in pancreatic adenocarcinoma and head-and-neck squamous carcinoma as well as in breast cancer, and in none of the seven cohorts does the fibrosis-adjusted residual predict survival adversely.

## N13.13 Code and data added for N13.9–N13.12

`02_scripts/build_deconv_signature.py` (single-cell signature matrix), `_fetch_cbio_panel.py` (generic cBioPortal panel fetch, `brca_tcga`/`paad_tcga`/`hnsc_tcga`), `run_bulk_stage10_deconv.py`, `run_bulk_stage11_metabric.py`, `run_bulk_stage12_pancancer.py`; new cached matrices under `01_data/raw/{brca_tcga,paad_tcga,hnsc_tcga}` and `01_data/processed/deconv_signature_gse176078.csv`. CIBERSORTx, xCell, TIDE and IMvigor210 were **not** used: the first needs a token and web upload, the second and fourth need R (unavailable), the third accepts only interactive submission. Every immune statement in this note is therefore a local, reproducible implementation and is labelled as such.

## N13.14 Stage 13 — ICB response, cross-indication (exploratory, not part of the main conclusions)

**Design and boundary.** GSE91061 (Riaz et al., *Cell* 2017, PMID 29033130) is the only public ICB cohort with a full expression matrix obtainable in this environment: 109 RNA-seq samples (51 pre-treatment, 58 on-treatment) from 65 melanoma patients treated with nivolumab ± ipilimumab, with RECIST categories (PD / SD / PR-CR) but no treatment-arm label and no PFS/OS in the GEO record. This is therefore a cross-indication, hypothesis-generating analysis. It is reported separately and does not enter the main result set or any multiplicity-corrected inference. The planned "TNBC-like" basal-keratin subset was not viable: basal keratins are essentially unexpressed in melanoma (interquartile width of the keratin score = 0.10), so all 49 evaluable pre-treatment samples were analysed together.

**Results** (`stage13_icb_logistic.csv`, `stage13_random_baseline.csv`, `stage13_paired_delta.csv`; logistic per-SD OR, median-split Fisher exact, 510 size-matched random sets of 60 genes, 42 patients with paired pre/on samples):

| module | OR (95% CI) | AUC | *P* (logistic) | random-set OR median | empirical *P* (AUC) |
|---|---|---|---|---|---|
| Perineural stroma | 1.06 (0.37–3.03) | 0.52 | 0.91 | 0.55 | 0.61 |
| Neural-intrinsic | 0.40 (0.03–6.38) | 0.42 | 0.52 | 0.55 | 0.22 |
| NPP58 | 1.53 (0.46–5.15) | 0.60 | 0.49 | 0.55 | 0.06 |
| CAF/ECM | 1.09 (0.45–2.66) | 0.52 | 0.85 | 0.55 | 0.57 |
| LRRC15 myCAF | 1.07 (0.45–2.55) | 0.52 | 0.88 | 0.55 | 0.57 |
| **CD8 cytotoxicity** | **1.62 (0.77–3.41)** | **0.62** | 0.20 | 0.55 | **0.016** |
| TGF-β core | 1.78 (0.48–6.58) | 0.59 | 0.38 | 0.55 | 0.06 |

Ten of 49 evaluable patients responded. No stroma- or neural-lineage module predicts response, and none differs from the random-set null (all empirical *P* ≥ 0.22); the direction is null in the paired pre→on analysis as well (42 pairs, all *P* > 0.09). The only read-out that separates from the null is CD8 cytotoxicity (AUC 0.62, empirical *P* = 0.016), consistent with melanoma ICB biology and with the chemotherapy cohort, where CD8 cytotoxicity was likewise the only replicable positive predictor. Treg reached nominal logistic significance (OR 2.31, *P* = 0.047) that does not survive Fisher exact or BH correction and is treated as noise.

**Interpretation.** The perineural stroma module carries no information about checkpoint-blockade outcome in this cross-indication setting, mirroring its null behaviour under neoadjuvant chemotherapy (N13.4). Because the cohort is melanoma, arm labels and survival data are unavailable, and n = 49, this observation is directional only; it is consistent with, but does not test, the absence of an immune-exclusion phenotype at the bulk level (N13.9, N13.10).

## N13.15 Code and data added for N13.14

`02_scripts/run_bulk_stage13_icb.py`; cached data `01_data/raw/gse91061/` (FPKM matrix, GEO sample metadata, symbol→Entrez map). Random-set scoring used the same 17,949-gene detected pool as the fixed modules. IMvigor210 could not be used (R package distribution only, no R runtime); this is declared rather than substituted.

## N13.16 Stage 15 — disease-free survival under neoadjuvant chemotherapy (GSE25065)

**Rationale.** Every survival analysis so far (METABRIC, TCGA-BRCA, GSE12276) sits in a primary-surgery, treatment-unstratified context. GSE25065 (Popovici et al. 2010) is the only local cohort with survival follow-up in a *neoadjuvant* context: 198 HER2-negative patients treated with taxane–anthracycline NAC, distant-relapse-free survival (DRFS) with 45 events over a median 3.2 years, plus pCR/RD, age and AJCC stage. Stage 4 used only its pCR label; this stage asks whether the modules predict relapse once the treatment context is fixed.

**Design.** Cox DRFS ~ module z-score, unadjusted and adjusted for age, AJCC stage (I–IV ordinal) and pCR, in four strata: ALL (n = 198), TNBC by IHC (n = 64), RD-only (residual disease, n = 140) and pCR-only (n = 42, 3 events — reported but uninformative). 510 size-matched random gene sets per module re-run through the same Cox pipeline (ALL and RD strata) as the null-calibrated baseline. A stage-ordinal mapping bug (roman-numeral prefix matching collapsed IIIA–C into stage II, leaving a near-constant covariate that produced complete separation) was caught and fixed before interpretation; the corrected stage distribution is II n = 107, III n = 89.

**Results** (`stage15_drfs_cox.csv`, `stage15_random_baseline.csv`):

| stratum | module | HR/SD (95% CI) | adj. *P* | random-set HR median | empirical *P* |
|---|---|---|---|---|---|
| ALL | Perineural stroma | 1.12 (0.84–1.47) | 0.33 | 1.04 | 0.57 |
| ALL | Neural-intrinsic | 1.19 (0.88–1.60) | 0.24 | 1.05 | 0.39 |
| ALL | CAF/ECM | 1.09 (0.82–1.45) | 0.38 | 1.02 | 0.65 |
| RD-only | Perineural stroma | 1.17 (0.87–1.58) | 0.19 | 1.08 | 0.47 |
| RD-only | CD8 cytotoxicity | **1.52 (1.10–2.08)** | **0.016** | — | — |
| ALL | CD8 cytotoxicity | 1.33 (1.00–1.77) | 0.050 | — | — |

No stroma- or neural-lineage module predicts DRFS in any stratum, and none separates from the random-set null (empirical *P* 0.33–0.65). The only signal is again CD8 cytotoxicity, here in the **adverse** direction: pre-treatment CD8 cytotoxic score predicts *worse* distant relapse, strongest in the residual-disease stratum (HR 1.52, adjusted *P* = 0.016). This direction agrees with the GSE12276 CD8 brain-relapse hazard (HR 1.43) and with the clinical literature on stromal TILs in post-NAC residual disease (Dieci et al.), while pre-NAC CD8 simultaneously predicts pCR benefit — an apparently paradoxical but internally consistent pattern: CD8 predicts *response*, and among non-responders marks residual aggressive disease. TNBC-IHC stratum: all modules null (HR 0.84–0.99).

## N13.17 Stage 16 — the genotype bridge (TCGA-BRCA mutations × modules)

**Rationale.** The bulk chain so far lives entirely in transcriptomic phenotype space. This stage crosses the modules with the genomic axis — tumour mutational burden and driver-gene status — in the 974 TCGA-BRCA samples with both expression (258-gene panel, Stage 7 cache) and mutation data. Mutation tables were fetched from cBioPortal (detailed projection, all 60,864 mutation records); the per-sample non-silent count from the fetch agrees perfectly with the clinical MUTATION_COUNT attribute (Spearman ρ = 1.000, n = 974), validating the pipeline.

**TMB.** Stroma modules show small negative raw correlations with log10 TMB (LRRC15 myCAF ρ = −0.11, TGF-β core ρ = −0.14, q < 0.001) that **vanish after adjustment for an epithelial purity proxy and MKI67** (|ρ| ≤ 0.045, all *P* > 0.16): the raw association is a tumour-content artefact. The only partial correlations that survive are CD8 cytotoxicity (ρ = +0.10, *P* = 0.0017) and IFN-γ (ρ = +0.10, *P* = 0.0018) — the canonical TMB → immune-infiltration axis, replicating known biology and serving as the positive control for this stage.

**Driver-gene status.** TP53-mutant tumours (n = 299) have lower LRRC15 myCAF (Cliff's Δ = −0.11, q = 0.028) and, expectedly, lower immune read-outs (CD8 Δ = +0.23 in WT, q = 3 × 10⁻⁸; IFN-γ q = 2 × 10⁻¹¹). PIK3CA-mutant tumours (n = 317) have *higher* stroma modules (perineural stroma Δ = +0.18, CAF/ECM Δ = +0.22, LRRC15 Δ = +0.24, all q < 2 × 10⁻⁴) and lower neural-intrinsic (Δ = −0.13, q = 0.002); however the 258-gene panel contains no ESR1/ERBB2/PGR, so **subtype confounding cannot be excluded** and this association is reported as descriptive. BRCA1/2 mutation carriers (combined n = 27) show the expected HRD immune phenotype (CD8 cytotoxicity Δ = +0.32, *P* = 0.005); their perineural-stroma elevation (Δ = +0.22, *P* = 0.056) is nominal only.

**Interaction.** Perineural stroma × TP53 on overall survival: HR 1.02 (95% CI 0.74–1.42; *P* = 0.87). No genotype stratifies the (already null) stroma–survival relationship.

**Interpretation.** The genomic axis does not rescue the module: TMB associations are purity artefacts, the TP53 interaction is null, and the remaining genotype associations either replicate known biology (TMB → CD8/IFN; HRD → CD8) or are subtype-confounded (PIK3CA). The genotype bridge thus closes on the same note as the phenotype bridge — the only cross-cohort replicable signals are immune read-outs, not the perineural stroma module.

## N13.18 Code and data added for N13.16–N13.17

`02_scripts/run_bulk_stage15_gse25065_drfs.py`, `02_scripts/run_bulk_stage16_genotype.py`, `02_scripts/_fetch_tcga_mutations.py`; cached data `01_data/raw/brca_tcga/brca_tcga_mutations.csv` (982 samples × mutation counts + TP53/BRCA1/BRCA2/PIK3CA status) and stage15/16 tables in `03_results/tables/bulk_extrapolation/`.

## N13.19 Stage 17 — spatial functional readouts: TLS–CAF output, collagen maturation and hypoxia

**Rationale.** Three hypothesised functional phenotypes of the perineural niche were tested at the spatial level with the manuscript's dual-calibration pipeline (rank-space partial Spearman vs distance, MSR randomisation, size-matched random gene sets, 300 sets at sizes 3–12 × 50): (i) near-nerve tertiary-lymphoid structures carrying stronger CAF output ("TLS–CAF co-programming"); (ii) collagen-crosslinking maturation (LOX/LOXL/PLOD/COL/FBN1/EMILIN1, 9 genes) declining more steeply than total collagen; (iii) a hypoxia–fibrosis coupling (10-gene signature CA9/VEGFA/SLC2A1/ANGPTL4/ADM/NDRG1/BNIP3/BNIP3L/DDIT3/PDK1).

**Results** (`stage17_gradient_meta.csv`, corrected size-bucketed null, resolution 1/51 per bucket). (i) The CAF–TLS output triple (CXCL12/FAP/POSTN) shows a weak negative raw gradient (ρ = −0.031) that does **not** exceed its size-matched null (empirical P = 0.098); TLS-marker fields are measurable in only 16/43 samples and, in raw distance terms, sit slightly *away* from nerve sources (TLS core +0.037, P_dir = 0.020) — the near-nerve TLS enrichment reported in the main text emerges only after CAF-content calibration, i.e. TLS are near nerve *relative to their stromal setting*, not in absolute distance. A direct coupling test (TLS-high spots; CAF output vs distance within those spots; near-vs-far levels) rejects the co-programming hypothesis: within TLS+ spots CAF output *increases* with distance (mean ρ = +0.06, P = 0.005, n = 16 samples), far TLS+ spots carry higher CAF output than near ones (+0.071 vs −0.006), and the near-vs-far Fisher-z difference is null (P = 0.83). (ii) Collagen-maturation machinery does not outperform total collagen (Coll_mat ρ = −0.026 vs null P = 0.49; LOX-only ρ = −0.010, P = 0.90; EMILIN1 −0.023, FBN1 −0.022, both below the size-1 baseline's tail) — the perineural collagen field is not preferentially crosslink-matured. (iii) The one direction-matched exceedance is **hypoxia, depleted near nerves**: the 10-gene hypoxia signature increases with distance (+0.043; core 4 genes +0.041; CA9 +0.027, SLC2A1 +0.031, VEGFA +0.031; each P_dir at the 1/51 floor, family FDR 0.051) — no size-matched random set reproduced a positive gradient of this magnitude. The nerve-proximal ~150 µm niche is therefore not a hypoxic fibrotic core; its stromal density is a vascularised, perineural-cuff-like architecture rather than the hypoxia–ECM coupling described in generic TNBC fibrosis.

## N13.20 Stage 18 — no CD8 effector/chemokine decoupling in bulk

**Rationale.** A hypothesised immunosuppressive phenotype — effector markers retained while the migration/chemokine axis is suppressed — was tested by partial Spearman of each readout against the stroma module controlling the CAF/ECM score, with 500 size-matched random gene sets (pool excluding all module/program/readout genes) in GSE163882 (n = 222), GSE25065 (n = 198) and METABRIC (n = 1,980; random null not applicable — the local 210-gene target file is pool-biased, partials reported descriptively).

**Results** (`stage18_cd8_decoupling.csv`). The decoupling hypothesis is rejected: cytotoxicity (CYT, CD8 cytotoxicity, IFN-γ) and the chemokine/migration axis (CXCL9/10/11, CCL5, CCL19/21, CXCL13; CCR7/CXCR3) travel **together**, all positively associated with the stroma module in every cohort (ρ_partial +0.10 to +0.39). The chemokine axis exceeds its size-matched random null in GSE25065 (+0.256, empirical P = 0.020; CXCL9/10/11 +0.237, P = 0.046), and nothing in GSE163882 exceeds a null whose median itself sits at +0.30 (a library-composition shared factor of the AmpliSeq platform — itself a caution against uncalibrated partial correlations). EMILIN1, the proposed TGF-β-antagonistic reverse control, does not replicate (GSE163882 −0.16 nominal, not beyond null; GSE25065 +0.18, P_emp = 0.14). The bulk immune shadow of the stroma module is therefore uniformly *aggregative* — effectors and chemokines together — the same picture as the spatial layer's TLS aggregation rather than exclusion, and the opposite of the hypothesised chemokine suppression.

## N13.21 Stage 19 — the bulk "neural" signal is Schwann/glial, not neuronal

**Rationale.** A canonical Schwann/myelinating-glia score (SOX10, MPZ, MBP, PMP22, PLP1, NCAM1, PRNP, NEFL, NEFH, S100B; 9-gene SOX10-excluded sensitivity variant) was projected onto the three full-transcriptome bulk cohorts with 500 size-matched random sets (pool excluding all module/nerve/SC genes).

**Results** (`stage19_schwann.csv`). The SC score correlates far more strongly with the neural-intrinsic module (ρ = 0.52 GSE25065, 0.57 TCGA, 0.77 GSE163882) than with the stroma module (0.25–0.53, exceeding the size-matched null only in TCGA, z = +2.6 to +2.9, empirical P ≤ 0.002; GSE163882/GSE25065 null-adjusted P ≈ 0.13). Removing SOX10 does not attenuate either correlation. SOX10 alone tracks basal keratins (ρ = 0.62–0.74), confirming that single-gene nerve markers in bulk breast data are subtype-confounded. Together with the spatial finding that nerve-centred ligand–receptor signalling is unmeasurable at spot resolution and the scRNA observation that Schwann cells are absent from dissociated references, the bulk-level "nerve signal" decomposes into a Schwann/glial component (carried by the neural-intrinsic module) and a matrix component (carried by the stroma module) — neither is a neuronal-intrinsic signal.

## N13.22 Calibration erratum for the axis-level random null (Note 12 and Table S20)

The size-matched random-gene-set calibration of the axis-gradient analysis (Note 12 / Table S20) standardised random sets into buckets keyed by **sample count** rather than gene-set size (a named-aggregation bug), which collapsed the null and put every empirical P at the resolution floor. The corrected, size-bucketed null (Note 12 table, regenerated Table S20) changes the status of two readouts: the 60-gene perineural-stroma module (FDR_emp 0.145) and the 20-gene myCAF axis (0.238) no longer exceed their size-matched baselines, while single matrisome markers (SPARC, TNC, ACTA2, POSTN, LRRC15), the collagen core, the apCAF axis and (descriptively) the away-from-nerve CD8/T panels retain direction-matched exceedance. The primary 511-set calibration of the functional readouts (Methods 2.9, Note 8) parses gene-set size correctly and is unaffected. The manuscript text (Section 3.7, Conclusion) has been amended accordingly; the substantive interpretation — perineural-specific excess carried by individual matrisome genes, set-level averaging dominated by the generic smoothness of the distance field — is the spatial counterpart of the bulk collinearity result and strengthens the paper's central methodological claim.

## N13.23 Code and data added for N13.19–N13.21

`02_scripts/run_bulk_stage17_tlscaf.py` (spatial functional readouts + TLS–CAF coupling), `run_bulk_stage18_decouple.py`, `run_bulk_stage19_schwann.py`, `_fix_random_calibration.py` (erratum re-computation for Note 12/Table S20); outputs in `03_results/tables/spatial_tlscaf/` (stage17_*) and `03_results/tables/bulk_extrapolation/` (stage18_*, stage19_*), plus `03_results/tables/axis_meta_summary.corrected.csv`.

## N13.24 Stage 20 — I-SPY2 chemo versus chemo+pembrolizumab pCR (GSE194040)

**Rationale and data.** The I-SPY2-990 mRNA resource (GSE194040) is the only public breast-cancer cohort pairing pre-treatment expression with pCR across randomised regimens: 988 samples / 987 patients, 14 arms, including paclitaxel-only (the control backbone; n = 179) and paclitaxel + pembrolizumab (n = 69). Clinical fields (arm, HR, HER2, pCR) were parsed from the family.soft records; the gene-level matrix (19,134 genes × 988) aligned to patient IDs at 99.8%. Analyses: module z-scores → pCR logistic per arm within the TNBC (HR−/HER2−) stratum; stroma × pembrolizumab interaction on the two-arm TNBC subset; 510 size-matched random gene sets for the stroma module in the TNBC-all stratum. I-SPY2 arms are **parallel randomisations, not head-to-head**; this stage is hypothesis-generating, cross-regimen, and excluded from primary conclusions.

**Results** (`stage20_ispy2_pcr.csv`). TNBC all arms (n = 362, 142 pCR): the stroma module is null (OR 0.90, 0.73–1.12; AUC 0.47 at the random-set median, empirical P = 0.998); the neural-intrinsic module shows a nominal protective direction (OR 0.75, 0.60–0.93) not carried through random-set calibration; CD8 cytotoxicity is again the only positive predictor (OR 1.51, 1.21–1.89; AUC 0.61). By arm (TNBC): paclitaxel-only (n = 85, 17 pCR) — stroma OR 1.48 (0.86–2.56, P = 0.16), the "fibrotic resistance" direction without significance; paclitaxel + pembrolizumab (n = 29, 19 pCR) — stroma OR 0.42 (0.17–1.07, P = 0.07), CAF/ECM OR 0.35 (0.12–0.98), CD8 OR 3.33 (1.08–10.3, AUC 0.74). The stroma × pembrolizumab interaction on the two-arm TNBC subset (n = 114) is nominally significant: OR 0.28 (0.10–0.84, P = 0.023) — i.e. the stroma-associated pCR disadvantage in the chemo-only direction is attenuated and reverses under checkpoint blockade. We flag three limits: 29 TNBC patients in the pembrolizumab arm; the counterintuitive CAF/ECM direction in that arm; and the absence of arm-level random-set calibration. Direction only, hypothesis-generating: a fibrotic stroma field that disadvantages chemotherapy response may not disadvantage — and by this interaction estimate may favour — immunotherapy response, consistent with the aggregative (TLS-rich) immune characterisation of Sections 4.2/4.3.

## N13.25 Code and data added for N13.24

`02_scripts/run_bulk_stage20_ispy2.py`; cached data `01_data/raw/gse194040/` (gene-level matrix 49.3 MB, family.soft clinical, clinical_full.csv); output `03_results/tables/bulk_extrapolation/stage20_ispy2_pcr.csv`.

## N13.26 Stage 14 — Xenium single-cell-resolution nerve-proximity validation (GSE350571)

**Rationale and data.** The main-text limitation that the neural-intrinsic module can be located only in situ (Limitation 2) called for single-cell-resolution platforms. GSE350571 provides six FFPE primary breast tumours on the Xenium Prime 5K platform (5,001 predesigned targets plus a 100-gene TLS/immune add-on; chemistry v2; onboard analysis 3.3.0.1) with cell-level counts, cell centroids and full transcript coordinates (0.43–0.67 million cells per section). Blocks were selected to be TLS-rich (selection bias declared); one section per patient; the cohort is not TNBC-enriched.

**Nerve definition (calibrated on this cohort).** Nerve bundles were defined as density kernels of axonal transcripts — ≥5 transcripts of ≥2 of *NEFH, TUBB3, GAP43* within 30 µm — each confirmed by myelin (≥2 *MPZ* transcripts within 50 µm). *PRNP* was excluded as diffuse background (P29: 148,396 transcripts, 43.8% in dense neighbourhoods), *UCHL1* for myoepithelial confounding, *SOX10* (basal-keratin proxy; Section 3.7) and *NCAM1* (NK-cell CD56). *MPZ* itself proved tumour-widely expressed in two sections (P26: 440,612 transcripts, 95.5% in dense neighbourhoods — an *MPZ*-expressing tumour, not myelin), which the axonal anchoring neutralises: with it P26 contributes 59 kernel points instead of a spurious 121,687. Nerve structures were detectable in four sections (P22: 327 kernel points; P26: 59; P29: 2,326; P32: 22); P17 and P24 sit at the axonal-density floor (4 and 0 kernel points) and are reported as nerve-negative sections. Composition contrasts are calibrated by torus-shift permutation of the nerve point cloud (999 per section).

**Results** (`stage14_xenium_nerve_summary.csv`, `stage14_xenium_fb_scores.csv`; Supplementary Figure S14). Fibroblasts within 150 µm of nerve kernels versus ≥300 µm away (identical marker-rule cell typing, mean-z module scores on log1p CP10K): the perineural-stroma module is higher proximally in three of four nerve-positive sections — P29 +0.211 (n = 13,221 vs 35,111), P32 +0.140 (529 vs 64,649), P22 +0.083 (3,251 vs 25,945) — and null in P26 (−0.008; proximal n = 721). The LRRC15-myCAF score is proximally higher in all four (+0.246, +0.215, +0.042, +0.032) and CAF/ECM in three (+0.269, +0.245, +0.056; P26 −0.100). CD8 composition ratios are directionally inconsistent (P22 1.36, permutation *P* = 0.042, enriched; P29 0.51, *P* = 0.046, depleted; P26 1.08 and P32 0.72, null) — reproducing the main-text picture of no consistent geometric CD8 exclusion. Two interpretive guards: proximal zones are tumour-enriched (P29 epithelial fraction 3.2–5.1× distal), so non-epithelial composition fractions drop mechanically (fractions sum to 1) and inference rests on the within-fibroblast module-score contrasts, which are immune to that coupling; and the sections are TLS-rich-selected, single per patient, and not TNBC-specific.

**Interpretation.** At single-cell resolution in an independent cohort, the nerve-proximal fibroblast carries a higher perineural-stroma and LRRC15-myCAF program in three of four evaluable sections — a directional, n = 4 replication of Section 3.7 Line 1 that requires no deconvolution reference (cell-type assignment is marker-rule-based). No meta-statistic is applied at four sections; the result is direction-level external validation, closing the loop the main text left open at Limitation 2.

## N13.27 Code and data added for N13.26

`02_scripts/run_bulk_stage14_xenium.py`, `02_scripts/_dl_gse350571.sh`, `02_scripts/_chunkdl_url.sh`; cached data `01_data/raw/gse350571/` (six-sample cells/h5/transcripts, gene panels, GSM manifest); outputs `03_results/tables/bulk_extrapolation/stage14_xenium_nerve_summary.csv` and `stage14_xenium_fb_scores.csv`. Technical notes: the GEO `*.parquet.gz` supplements are gzip-wrapped parquet and were stream-decompressed before predicate-pushdown reading; the chunked downloader was hardened (strict part-name concatenation, integrity-gated skip) after a retry-tail contamination incident (P17 h5 and P24/P26 transcripts re-fetched and byte-verified).

## N13.28 Stage 21 — sensory-nerve and CGRP-axis proximity gradients, and the near-nerve TLS functional ratio (dialogue with Cell 2026)

**Rationale.** A recent TNBC study (Zhang S-W et al., *Cell* 2026, PMID 41650969; 360 clinical samples plus mouse models) demonstrated causally that sensory nerves drive immune exclusion through NGF→neuron→CGRP→RAMP1⁺-CAF→cAMP/PKA/CREB1 collagen deposition, with rimegepant (CGRP-receptor antagonist) restoring infiltration. Two questions follow for our spot-level framework: (i) is the sensory-nerve/CGRP signal measurable in the ~150 µm perineural field at all; (ii) does the near-nerve TLS aggregation (main text Section 4.2, Z = −2.10 to −2.20) carry a functional deficit — the "aggregated but anergic" phenotype.

**Part A — axis gradients** (`stage21_axis_gradient_meta.corrected.csv`; 43 sections, dual calibration as in Stage 17; the initial meta suffered the same named-aggregation bucket bug as N13.22 and was recomputed from `stage21_axis_gradient_per_sample.csv.gz` by extracting gene-set size from the `RAND|s{size}|` prefix — script `_fix_stage21_calibration.py`; Supplementary Figure S13). Read-outs: CGRP axis (CALCA/CALCRL/RAMP1/CALCR), CGRP ligand–receptor core (CALCA/RAMP1), NGF–TrkA axis (NGF/NTRK1/NGFR/SORT1), sensory-nerve markers (TRPV1/SCN9A/SCN10A/S100B) and single genes. Findings: the sensory-marker set is nerve-proximal at ρ = −0.072 [−0.115, −0.028] (n = 10 samples passing expression thresholds; size-bucket empirical *P* = 0.016, direction-matched *P* = 0.016); *S100B* alone replicates it at ρ = −0.061 [−0.097, −0.025], n = 27, above the size-1 bucket null — the strongest non-matrix nerve-proximal read-out measured in this manuscript. The CGRP-axis proxy is weakly nerve-proximal (ρ = −0.043 [−0.102, 0.016], n = 18, above the size-4 null *P* = 0.007 global / 0.016 bucket, but sample-heterogeneous); the CALCA–RAMP1 two-gene core does not exceed its null (ρ = −0.026, bucket *P* = 0.18), consistent with Limitation 4 (the pair is assessable as a ligand–receptor test in only 2 of 43 samples). The NGF–TrkA axis, and NGF, NTRK1, NGFR, CALCA, CALCR, SCN9A, SCN10A individually, fall below the 10-sample measurability floor. The bulk counterpart (METABRIC, Step 4b) gives CGRP-axis ~ NPP ρ = +0.158 (P = 1.8 × 10⁻¹²) with no survival content.

**Part B — TLS near-vs-far functional ratio** (`stage21_tls_function_meta.csv`, `stage21_tls_function_per_sample.csv`; TLS⁺ spots defined identically to Stage 17 — TLS-extended score above the within-sample median; near ≤ 150 µm, far ≥ 300 µm; per-sample Mann–Whitney, one-sample *t* across sections; median 77 near / 378 far TLS⁺ spots). The "aggregated but anergic" hypothesis is rejected: cytotoxicity (GZMA/GZMB/PRF1/NKG7/IFNG) near-minus-far Δ = −0.038 [−0.102, 0.027], *P* = 0.28; exhaustion (PDCD1/CTLA4/LAG3/HAVCR2/TIGIT) Δ = −0.058 [−0.100, −0.016], *P* = 0.024 — **depleted**, i.e. the direction opposite to anergy; TGF-β-response (SERPINE1/TGFB1/TGFB2/TGFBR2/SMAD3) Δ = −0.023, *P* = 0.34; the cytotoxicity-minus-exhaustion contrast Δ = +0.020, *P* = 0.48 (n = 10–12 sections per panel, 12 for cytotoxicity). Near-nerve TLS are not functionally impaired relative to distant TLS at spot resolution.

**Interpretation.** Three statements bound what Visium resolution can add to the Cell 2026 mechanism. (1) The macroscopic half is measurable and directionally consistent: the field around nerves is sensory-marker-positive and carries a weak CGRP-axis transcript proxy, co-located with the myCAF output field of Section 3.7. (2) The causal chain's two ends are below resolution: the NGF driver is unmeasurable, and the CALCA–RAMP1 pair cannot be tested as receptor–ligand co-localisation (Limitation 4). (3) The immune end differs by scale: in PNI-stratified TNBC the Cell 2026 phenotype is tumour-level immune exclusion with TLS loss, whereas at field level the same ~150 µm niche aggregates TLS whose effector output is intact (this stage). A two-scale reading — neural CGRP signalling remodels matrix below spot resolution, while the ≥55 µm field-level consequence is myCAF output with lymphocyte aggregation rather than geometric exclusion — keeps both observations, and their divergence (tumour-level exclusion vs field-level aggregation) is itself a testable question for co-registered Visium HD or Xenium with neuronal and CGRP-pathway markers.

## N13.29 Code and data added for N13.28

`02_scripts/run_bulk_stage21_sensory_axis.py` (Part A gradient blocks through the Stage 17 dual-calibration pipeline plus Part B TLS functional panels), `02_scripts/_fix_stage21_calibration.py` (size-bucket recomputation; the in-script meta retained the N13.22 named-aggregation defect and was superseded); outputs `03_results/tables/spatial_tlscaf/stage21_axis_gradient_per_sample.csv.gz`, `stage21_axis_gradient_meta.corrected.csv`, `stage21_tls_function_per_sample.csv`, `stage21_tls_function_meta.csv`.

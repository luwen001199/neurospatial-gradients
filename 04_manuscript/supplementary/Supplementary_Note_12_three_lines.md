# Supplementary Note 12 — Three-line characterisation of the perineural stroma

This note documents the full statistical detail behind Results Section 3.7 (manuscript) and the three-line analysis of the perineural stroma: (1) the calibrated myCAF/ECM gradient, (2) immune-exclusion geometry, and (3) CAF-axis layering along the nerve-rooted trajectory. All analyses run in `02_scripts/run_axis_gradient.py`, `run_immune_exclusion.py` and `run_caf_axes.py`; outputs are `01_data/processed/axis_*.csv`, `immune_exclusion_*.csv`, `caf_axis_*.csv` and `03_results/tables/axis_meta_summary.csv`, `immune_exclusion_summary.csv`. All values below are from the full 43-sample cohort.

## N1. Line 1 — readouts and dual calibration

Readouts (mean *z* of member genes per spot):

| readout | size | members |
|---|---|---|
| Perineural-stroma module | 60 | module definition of Note 10 (unchanged) |
| Collagen core | 13 | COL* genes within the module |
| myCAF axis | 20 | Note-6 myofibroblastic-CAF panel minus *LRRC15* |
| iCAF axis | 27 | inflammatory/dev-CAF panel |
| apCAF axis | 8 | *HLA-DRA, HLA-DRB1, HLA-DPA1, HLA-DPB1, CD74, CIITA, LAMP3, CD70* |
| LRRC15⁺ axis | 1 | *LRRC15* |
| single markers | 1 each | *LRRC15, ACTA2, POSTN, TNC, COL1A1, SPARC, FN1, LOXL2* |
| Macrophage panel | 15 | Note-6 macrophage/monocyte panel |
| CD8 T panel | 10 | descriptive only (measurable in 16/43 samples) |
| Pan-T panel | 13 | descriptive only (measurable in 17/43 samples) |

Statistics per sample: rank-space partial Spearman against distance-to-nearest-nerve, controlling for log-total-counts and percent-mitochondrial (Gram–Schmidt, identical to Methods 2.2). Nulls: (i) Moran spectral randomisation, 199 draws, 40-spot coefficient bins (identical to Note 7); (ii) 400 size-matched random gene sets at sizes {1, 8, 20, 27, 60} × 80, drawn from the same expressed-gene universe as the readout pipeline. Cross-sample inference: random-effects *t* test on Fisher-*z* coefficients; Stouffer combination of MSR-calibrated *u*; BH at 0.05 on both families; empirical *P* direction-matched within the correct size bucket (resolution 1/81 per bucket). *(Calibration correction: an earlier implementation standardised random sets into buckets keyed by sample count rather than gene-set size, collapsing the size-matched null and putting every empirical P at the resolution floor; the table below is from the corrected, size-bucketed null, and the same correction was applied to Supplementary Table S20.)*

Full results (43 samples; LRRC15 34/43, TNC 38/43 measurable):

| readout | ρ mean (95% CI) | P(RE-t) | Stouffer Z | P(Stouffer) | baseline ρ | P_emp(dir) | FDR_emp | pass MSR | pass random | pass both |
|---|---|---|---|---|---|---|---|---|---|---|
| apCAF axis | −0.071 (−0.119, −0.021) | 0.008 | −4.94 | 7.9 × 10⁻⁷ | 0.024 | 0.012 | 0.030 | yes | yes | yes |
| Perineural-stroma module | −0.070 (−0.125, −0.014) | 0.018 | −4.37 | 1.2 × 10⁻⁵ | 0.059 | 0.111 | 0.145 | yes | **no** | no |
| Collagen core | −0.057 (−0.114, −0.000) | 0.056 | −3.38 | 7.2 × 10⁻⁴ | 0.039 | 0.025 | 0.042 | yes | yes | yes |
| iCAF axis | −0.055 (−0.098, −0.012) | 0.016 | −3.84 | 1.2 × 10⁻⁴ | 0.044 | 0.049 | 0.070 | yes | (borderline) | no |
| TNC (n = 38) | −0.054 (−0.083, −0.026) | 0.0006 | −5.34 | 9.5 × 10⁻⁸ | 0.012 | 0.012 | 0.030 | yes | yes | yes |
| SPARC | −0.052 (−0.106, +0.002) | 0.066 | −3.36 | 7.9 × 10⁻⁴ | 0.012 | 0.012 | 0.030 | yes | yes | yes |
| myCAF axis | −0.047 (−0.097, +0.003) | 0.072 | −3.23 | 1.3 × 10⁻³ | 0.039 | 0.210 | 0.238 | yes | **no** | no |
| ACTA2 | −0.038 (−0.077, +0.002) | 0.068 | −3.14 | 1.7 × 10⁻³ | 0.012 | 0.012 | 0.030 | yes | yes | yes |
| POSTN | −0.037 (−0.081, +0.008) | 0.119 | −2.77 | 5.6 × 10⁻³ | 0.012 | 0.012 | 0.030 | yes | yes | yes |
| Macrophage panel | −0.033 (−0.087, +0.021) | 0.234 | −3.00 | 2.7 × 10⁻³ | 0.039 | 0.716 | 0.716 | yes | no | no |
| LRRC15 (n = 34) | −0.030 (−0.064, +0.004) | 0.088 | −2.39 | 0.017 | 0.012 | 0.025 | 0.042 | yes | yes | yes |
| COL1A1 | −0.026 (−0.083, +0.031) | 0.375 | −1.42 | 0.16 | 0.012 | 0.037 | 0.057 | no | (borderline) | no |
| FN1 | −0.020 (−0.071, +0.031) | 0.451 | −1.61 | 0.11 | 0.012 | 0.148 | 0.180 | no | no | no |
| LOXL2 | −0.016 (−0.043, +0.011) | 0.243 | −1.97 | 0.048 | 0.012 | 0.235 | 0.249 | no | no | no |
| CD8 T (n = 16, descriptive) | +0.010 | 0.730 | +0.30 | 0.76 | 0.024 | 0.012† | 0.030† | no | (yes, away) | no |
| Pan-T (n = 17, descriptive) | +0.015 | 0.639 | +0.43 | 0.67 | 0.039 | 0.012† | 0.030† | no | (yes, away) | no |

† Direction-matched exceedance in the *positive* (away-from-nerve) direction: no size-matched random set of 8–13 genes produced a positive gradient of this magnitude, a descriptive observation consistent with the raw ring profiles.

The corrected size-bucketed null changes one earlier claim. Single markers (*TNC, SPARC, ACTA2, POSTN, LRRC15*), the collagen core and the apCAF axis exceed their size-matched baselines (FDR_emp 0.030–0.042), and the iCAF axis is borderline (0.070). However, the **60-gene module and the 20-gene myCAF axis do not** (FDR_emp 0.145 and 0.238): random 60- and 20-gene sets average member *z*-fields into progressively smoother fields whose expected |ρ| grows with set size (baseline 0.012 at 1 gene → 0.059 at 60), so a 60-gene gradient of −0.070 is exactly what any equally large smooth set produces. The perineural-*specific* excess is therefore carried by the individual matrisome genes, not by set-level averaging — the spatial counterpart of the bulk finding that the module is ρ = 0.97 collinear with generic desmoplasia. Both statements are two faces of the same phenomenon, and we report the module-level gradient as MSR-validated but not size-matched-null-exceeding. Effect sizes are uniformly smaller than in the 30-sample interim (e.g. module −0.112 → −0.070, SPARC −0.083 → −0.052, LRRC15 −0.056 → −0.030), the expected direction as cohort size grows.

Module gradient under the module-excluded deconvolution panel (sensitivity; see N4): deconvolved myCAF abundance gradient ρ = −0.078 (module-excluded) versus −0.062 (primary), both Fisher-*z* meta *P* < 10⁻⁴⁷ — the gradient is not a reference-panel artefact.

Cross-reference to the single-cell atlas (Note 11): 40 of 60 measurable perineural-stroma module genes peak in the myCAF-like compartment (COL1A2 99%, COL1A1 98%, COL3A1 98%, SPARC 98%, LUM 98%, DCN 95%, POSTN 88% of myCAF-like cells expressing); *ACTA2* is shared with pericytes (90% PVL vs 65% myCAF-like); neural-intrinsic genes are dispersed (TUBA1A 75% and PMP22 71% myCAF-like, UCHL1 20% myCAF-like, TUBB3 10% Cancer Epithelial). LRRC15 is expressed by 11.5% of myCAF-like cells — a sparse but specific subset, consistent with its role as a tumour-boundary activation marker rather than a household CAF marker. Directional cross-check in GSE75688: 48 of 59 measurable module genes peak in stromal cells.

## N2. Line 2 — immune-exclusion geometry

Definitions per sample (n = 43):
- Residualised fields: rank( CD8-density ) and rank( TLS-score ), Gram–Schmidt residualised against log-total-counts and percent-mitochondrial (`rho_CD8`, `rho_TLS`); a second variant additionally removes rank( myCAF-score ) (`*_netcaf`).
- Immune exclusion index IEI = log₂[(mean residualised field in the 0–50 µm ring + 0.25)/(mean in the >400 µm ring + 0.25)]; negative = proximal depletion. Requires ≥5 spots per ring (35/43 samples informative).
- Perineural trend calibrated against the MSR null of Note 7 (199 draws/sample); cross-sample Stouffer.

| metric | n | median | IQR | frac negative | Wilcoxon P | Stouffer Z | P |
|---|---|---|---|---|---|---|---|
| ρ (CD8 vs distance) | 43 | −0.022 | −0.073…+0.045 | 0.60 | 0.16 | — | — |
| u_MSR (CD8) | 43 | −0.52 | — | 0.60 | 0.058 | **−3.14** | **0.0017** |
| u_MSR (CD8, net of CAF) | 43 | −0.48 | — | 0.58 | 0.075 | **−2.98** | **0.0028** |
| IEI (CD8) | 35 | +0.22 | −0.58…+1.23 | 0.40 | 0.30 | — | — |
| ρ (TLS vs distance) | 43 | −0.022 | −0.071…+0.063 | 0.60 | 0.41 | — | — |
| u_MSR (TLS) | 43 | −0.52 | — | 0.60 | 0.21 | **−2.20** | **0.028** |
| u_MSR (TLS, net of CAF) | 43 | −0.44 | — | 0.60 | 0.24 | **−2.10** | **0.036** |
| IEI (TLS) | 35 | **+0.77** | −0.16…+1.20 | 0.26 | **0.024** | — | — |
| ρ (myCAF vs CD8, ≤100 µm) | 37 | +0.135 | +0.036…+0.204 | 0.14 | **2.9 × 10⁻⁵** | — | — |
| ρ (myCAF vs CD8, >200 µm) | 43 | +0.102 | +0.032…+0.165 | 0.16 | **2.0 × 10⁻⁶** | — | — |

Interpretation. At the full cohort the calibrated CD8 marker field shows a weak but significant perineural decline (Stouffer Z = −3.14, P = 0.0017; net of CAF Z = −2.98), an order of magnitude weaker than the stromal gradients of N1 and invisible to the geometric exclusion index, which remains null (IEI +0.22, P = 0.30). The TLS field is reproducibly nerve-proximal (Z = −2.20; net of CAF Z = −2.10; IEI +0.77, P = 0.024) and myCAF co-localises positively with CD8 density in both compartments. We read the weak CD8 decline as a TLS-composition effect on the marker field rather than effector exclusion. The deconvolution-based T-cell gradient of Note 11 (ρ = +0.018, i.e. relative T-cell sparsity near nerves) describes cell-type *composition*; the TLS/CD8 marker-field statistics describe lymphocyte *aggregation* — different contrasts that we report without adjudicating, flagging the discrepancy as resolution- and panel-dependent.

## N3. Line 3 — CAF-axis layering and nerve-rooted trajectory

Two complementary decompositions. First, seeded-NNLS deconvolution (Methods 2.8; per-spot cell-type abundance fractions) by ring:

| ring | 0–50 µm | 50–100 | 100–200 | 200–400 | >400 µm |
|---|---|---|---|---|---|
| deconvolved myCAF fraction | 0.277 | 0.279 | 0.266 | 0.251 | — |
| deconvolved Myeloid fraction | 0.134 | 0.133 | 0.132 | 0.127 | — |
| deconvolved T-cell fraction | 0.033 | 0.032 | 0.034 | 0.036 | — |

Second, NNLS decomposition onto the four CAF-axis profiles (`run_caf_axes.py`; Figure 4E shows mean weights per ring):

| axis | 0–50 µm | 50–100 | 100–200 | 200–400 | >400 µm |
|---|---|---|---|---|---|
| myCAF | 1.051 | 1.034 | 0.960 | 0.827 | 0.819 |
| iCAF | 0.651 | 0.665 | 0.625 | 0.563 | 0.569 |
| apCAF | 0.721 | 0.737 | 0.695 | 0.673 | 0.536 |
| LRRC15⁺ | 0.331 | 0.330 | 0.333 | 0.297 | 0.325 |

Both decompositions agree on the ordering: myCAF carries the largest nerve-proximal weight/fraction, and the remaining axes are lower and flatter.

Nerve-rooted pseudotime (Methods 2.5 pipeline, recomputed per sample; root = highest-scoring nerve spot, in the 0–50 µm ring in 43/43 informative samples). Axis score vs pseudotime, per-sample partial Spearman (| UMI, MT), cross-sample random-effects *t*:

| axis | median ρ (dpt) | mean ρ | P | bootstrap 95% CI | ρ (distance) | P (distance) |
|---|---|---|---|---|---|---|
| myCAF | −0.279 | −0.223 | 0.0092 | −0.355…−0.083 | −0.043 | 0.15 |
| apCAF | −0.250 | −0.204 | 0.0066 | −0.320…−0.075 | −0.083 | 0.013 |
| iCAF | −0.205 | −0.172 | 0.019 | −0.289…−0.036 | −0.044 | 0.10 |
| LRRC15⁺ | −0.035 | −0.063 | 0.099 | −0.128…+0.005 | +0.010 | 0.26 |

Stability: per-sample 80% spot resampling (100 draws) gives percentile intervals excluding zero in the pooled direction for all mature-CAF axes; sample-level bootstrap of the cross-sample mean (200 draws) is tabulated above. Negative ρ = the axis is highest at the nerve-rooted (pseudotime-0) end. The trajectory, which aggregates expression structure beyond a single distance readout, yields stronger and more significant ordering than the distance field for myCAF and iCAF, whereas the apCAF distance gradient (P = 0.013) is comparable.

## N4. Caveats

1. **Panel dependence of spot-level contrasts.** The nerve-spot enrichment of myCAF abundance flips sign under the module-excluded deconvolution panel (mean log₂FC −0.018 → +0.071) while the distance gradient is stable and, if anything, strengthens (ρ = −0.062 → −0.078; both Fisher-*z* meta *P* < 10⁻⁴⁷). Proximal *enrichment* contrasts are therefore panel-dependent; *gradient* claims are not. Only the latter are used inferentially.
2. **TLS score** is a composite marker field (Note 6), not a histological annotation; the nerve-proximal TLS enrichment should be verified against mIF of matched sections.
3. **apCAF markers** (HLA-II/CD74) are shared with antigen-presenting myeloid cells; the apCAF axis should be read as an antigen-presentation axis, not a fibroblast-specific identity.
4. The exclusion index is a geometric summary of a residualised field; its null distribution is calibrated through the MSR machinery on the underlying field, not on the ratio itself.
5. **Cohort-size dependence.** All Line 1 effect sizes shrank from the 30-sample interim to the full 43-sample cohort (N1); calibrated gradients should be quoted with full-cohort values (manuscript Limitation 8).

# Supplementary Note 10 — Decomposition of the NPP into neural-intrinsic and perineural-stroma modules

## 10.1 Module definitions

The program was decomposed by curation into a **neural-intrinsic** module (39 genes: neuronal, glial and Schwann structural proteins, transmitter-synthesis enzymes, transmitter receptors and neuropeptides — *UCHL1, TUBB3, TUBA1A, MAP2, NEFL, NEFM, PLP1, MPZ, TH, DBH, CHAT, TRPV1, CALCA, CALCB, VIP* and others) and a **perineural-stroma** module (60 genes: matrix, collagen, CAF-activation and remodelling genes — *POSTN, SPARC, FN1, collagens, TAGLN, ACTA2, MMPs, LRRC15, THBS1* and others). Full membership: `TableS8_module_definitions.csv`.

## 10.2 Spatial gradients (43 samples)

Both modules decline towards nerve sources (mean partial ρ: neural-intrinsic −0.122, consistent in 41/43 samples; perineural-stroma −0.068, consistent in 28/43). After each is partialled out of the other in rank space, both retain a significant gradient (−0.106 and −0.054) — two separable signals sharing one spatial direction rather than one signal measured twice.

## 10.3 Marker-removal sensitivity

The neural-intrinsic module retains its gradient after **every gene used to define nerve spots is removed from it** (ρ = −0.054, consistent in 37/43 samples; MSR Stouffer Z = −6.71, *P* = 1.9 × 10⁻¹¹). Its perineural enrichment is therefore not an artefact of the spot-calling markers. This is the single most important control in the paper: it is what separates "nerve-associated transcription" from "a restatement of how we called the spots".

Source table: `TableS8b_module_spatial_gradients.csv`.

## 10.4 Bulk correlations (METABRIC, n = 1,976)

| Pair | ρ | *P* |
|---|---|---|
| Neural-intrinsic ~ Perineural-stroma | 0.050 | 0.026 |
| Neural-intrinsic ~ CAF/ECM | 0.040 | 0.075 |
| Perineural-stroma ~ CAF/ECM | 0.975 | <10⁻³⁰⁰ |
| Neural-intrinsic ~ original NPP | −0.101 | 6.6 × 10⁻⁶ |
| Perineural-stroma ~ original NPP | 0.230 | 4.4 × 10⁻²⁵ |

The two modules are nearly uncorrelated with each other; the neural-intrinsic module is essentially orthogonal to generic desmoplasia (ρ = 0.040, not significant), whereas the perineural-stroma module is almost a duplicate of it (ρ = 0.975). The 58-gene NPP sits between them and is dominated by the stromal half.

## 10.5 Survival

Per-SD hazard ratios: neural-intrinsic RFS 1.03 (0.97–1.09), stroma 1.00 (0.93–1.07); after mutual adjustment 1.05 and 1.01; 500-replicate bootstrap distributions centred on 1.0 for both (95% CI 0.99–1.18 and 0.94–1.08). Neither module is prognostic. In an exploratory model containing both modules **and** the generic CAF/ECM score, the stroma module flipped protective (RFS HR = 0.62, *P* = 0.006; OS HR = 0.73, *P* = 0.037) — a suppression effect consistent with its near-collinearity with desmoplasia. It is reported as exploratory because it was not bootstrap-validated.

## 10.6 Incremental prognostic value

A ridge- and lasso-penalised Cox model combining spatial-gradient and cell-composition features added at most ΔC-index = +0.008 (cross-validated) over a clinical baseline (age, nodal status, ER status, systemic therapy, CD8 cytotoxic score) for RFS and −0.000 for OS. The only variable with a stable independent signal was cytotoxic T-cell content (RFS HR = 0.88, *P* = 0.0019). Source table: `TableS12_prognostic_incremental_value.csv`.

## 10.7 Conclusion

Published PNI signatures that mix neural-intrinsic with perineural-stroma genes inherit the desmoplastic confound in proportion to their stromal content. Splitting the two explains why bulk-level hazard ratios for such signatures have been unstable across studies: the hazard ratio of the mixture depends on the mixing ratio, which is an arbitrary product of gene-set construction.

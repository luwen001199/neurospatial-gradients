# Supplementary Note 7 — Moran spectral randomisation: algorithm, binning sensitivity, and the rank-transform correction

## 7.1 Algorithm

For each sample we build a binary symmetric 6-nearest-neighbour graph on spot centroids (edge distance ≤ 220 µm), double-centre it (H = PHP) and eigendecompose. The centred, rank-transformed expression matrix is projected onto the eigenvectors; projection coefficients are permuted at random **within eigenvalue bins** (bins of ≈40 eigenvalues, ordered by eigenvalue), with the same permutation applied across genes so that the gene–gene covariance structure is preserved as well; the permuted coefficients yield the surrogate statistic in closed form without reconstructing the spatial field. With *C* the permuted coefficient matrix, *w* = V′r (r the residualised distance ranks) and *U* = V′Q (Q the orthogonalised covariates),

ρ = w′C_g / √[(‖C_g‖² − ‖U′C_g‖²)·‖r‖²].

Each surrogate retains the gene's Moran's *I* (median retention 96% of the observed value) while randomising its spatial phase. 499 surrogates per gene per sample. Reference implementation: Wagner & Dray (2015), `adespatial::msr` [29].

## 7.2 Eigenvalue-binning sensitivity

Coarser bins preserve autocorrelation less faithfully but randomise more freely. With bins of ≈100 eigenvalues, Moran's-*I* retention falls to 87% (from 97% at bins of 40) and the null becomes permissive: 56 of 56 program genes "survive" — a null that rejects nothing is not a calibration. With bins of 40, 39 of 54 up-program genes survive (hypergeometric *P* = 9.2 × 10⁻⁷⁷ against 69 MSR-significant genes in the 10,503-gene universe). We report the finer binning as primary because it is the configuration that actually preserves the target autocorrelation, and we note the sensitivity explicitly because the choice is consequential.

## 7.3 The rank-transform correction and its consequences

An early implementation of the average-rank transform mishandled ties: the cumulative offset for tied groups was accumulated in the wrong direction relative to the sorted order. Because Visium counts contain a very large number of exact zeros, tie handling is not an edge case here — it is the modal case. The corrected transform was verified element-wise against `scipy.stats.rankdata` on every sample.

Effects of the correction:
- The testable gene pool changed from 6,221 to 10,503 genes (43-sample cohort).
- The METABRIC association that had appeared protective under the pre-correction pipeline (OS residual HR = 0.910, *P* = 0.0013) disappeared; under the corrected pipeline the same quantity is adverse (OS residual HR = 1.139, *P* = 8 × 10⁻⁶). We attribute the earlier protective estimate to the implementation error, not to biology, and we report both in `TableS3_signature_construction_sensitivity.csv`.
- The bootstrap/E-value table (Note 6) was recomputed accordingly.

## 7.4 What MSR does and does not correct

MSR preserves each gene's own spatial autocorrelation while breaking its pairing with the distance field. It therefore corrects the *autocorrelation* inflation. It does **not** correct a second, independent inflation: because the distance field is near-deterministically smooth (median Moran's *I* = 0.97), any large gene set averaged into a score produces a smooth field that correlates with it. Quantified in Note 8, this is why a size-matched random-gene-set baseline is required as a second calibration stage.

Source tables: `TableS5_MSR_per_gene_statistics.csv`, `TableS2_NPP_per_gene_statistics.csv`. Figure: Supplementary Figure S5.

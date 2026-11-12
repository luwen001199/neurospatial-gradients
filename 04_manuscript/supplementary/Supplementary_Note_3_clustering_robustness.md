# Supplementary Note 3 — Robustness of the sample-level spatial partition

**Scope.** The two-architecture partition of samples (Neural-Active n = 15, Neural-Cold n = 28) is reported in the main text as exploratory. This note documents the four robustness workstreams that support that qualifier.

**K selection.** Consensus clustering (80% subsampling, 1,000 iterations, Ward linkage on *z*-scaled features) was run for K = 2–4. K = 2 maximised the consensus-matrix silhouette (0.72) against K = 3 (0.55) and K = 4 (0.46); PAC at K = 2 was 0.61. Re-running under four random seeds (300 subsamples each; three completed) gave K = 2 silhouettes 0.719–0.736, so the *choice of K* is seed-stable.

**Partition stability.** Three stress tests were applied. (i) 80% subsample reclustering (1,000 replicates) reproduced the reference labels with median adjusted Rand index 0.66 (IQR 0.33–0.88) — moderate, not good. (ii) Feature-wise jackknife: ARI 0.32–0.91 across the nine spatial features, the myCAF-gradient and decay-constant features being the least stable inputs and the nerve-signal and CD8 features the most stable. (iii) Algorithm substitution at K = 2: complete linkage ARI 0.91 and K-means 0.82, but average linkage −0.02 and full-covariance Gaussian mixture −0.05 — agreement with the reference partition is algorithm-dependent.

**Interpretation.** K = 2 is consistently selected, but the labels themselves are moderately reproducible in subsampling and highly sensitive to the clustering algorithm. With 43 samples we do not consider the partition a finding. It is used nowhere in the inferential chain: the NPP, the MSR calibration, the functional readouts, the pseudotime and the module decomposition are all computed per sample and combined by meta-analysis, never by cluster.

Source table: `TableS16_clustering_robustness.csv`. Figures: Supplementary Figures S2 and S7.

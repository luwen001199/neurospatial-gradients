# Supplementary Note 6 — Bootstrap and E-value analysis of the bulk survival associations

**Purpose.** The main text declines to claim a prognostic association for the NPP. This note shows that the decision is not driven by a single model fit but survives resampling, and quantifies how much unmeasured confounding would be needed to explain the residual estimates.

**Design.** Four Cox models — RFS and OS × naive NPP and stroma-adjusted (residual) NPP — fitted with the multivariable covariate set used in the main text (age, nodal status, ER status, chemotherapy, hormone therapy and generic CAF/ECM content), all exposures *z*-scored so that hazard ratios are per SD. 500 nonparametric bootstrap resamples per model; percentile confidence intervals and the fraction of resamples with HR < 1. VanderWeele–Ding E-values computed on the point estimate and on the confidence limit closest to the null.

Source table: `TableS10_bootstrap_and_evalue.csv`. Figure: Supplementary Figure S3.

**Results (revised).** All four estimates are adverse and attenuated:

| Outcome | Exposure | HR (95% CI) | *P* | Bootstrap HR (95% CI) | E-value (CI limit) |
|---|---|---|---|---|---|
| RFS | NPP | 1.063 (0.976–1.157) | 0.161 | 1.058 (0.977–1.160) | 1.32 (1.18) |
| RFS | residual NPP | 1.064 (0.991–1.143) | 0.086 | 1.061 (0.992–1.144) | 1.33 (1.10) |
| OS | NPP | 1.058 (0.985–1.136) | 0.124 | 1.053 (0.982–1.137) | 1.30 (1.14) |
| OS | residual NPP | 1.056 (0.995–1.121) | 0.071 | 1.056 (0.992–1.125) | 1.30 (1.07) |

Bootstrap intervals match the model-based intervals closely, so the attenuation is not a small-sample artefact of a single fit. The residual-NPP estimates sit at *P* = 0.086 and 0.071 — borderline, not significant — and the E-values of 1.30–1.33 (confidence limit 1.07–1.18) mean that an unmeasured confounder associated with both the program and the outcome by a risk ratio of about 1.3 would suffice to explain the point estimate, and one of about 1.1 would suffice to explain the confidence limit. These are small effects, well within the range of routine clinical confounding.

**Note on a superseded version.** An earlier run of this analysis reported protective estimates (RFS HR 0.890, OS HR 0.871) with E-values 1.4–1.6. Those numbers were computed against an NPP score that predated the rank-transform correction described in Note 7 and are superseded. The table above is the current, corrected version.

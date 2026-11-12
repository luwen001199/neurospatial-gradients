# Supplementary Note 4 — Audit of cohort candidates that could not support the analysis plan

The original study design included two analyses that were abandoned for lack of suitable public data. This note records what was found, so that the negative result is auditable rather than merely asserted.

**Immune-checkpoint blockade response.** Requirement: TNBC, bulk or spatial expression, documented response to anti-PD-1/PD-L1.
- *GSE169246* — single-cell RNA-seq of sorted CD45⁺ cells. Cell-type restriction removes the stromal and epithelial compartments that carry the program, and the series carries no response annotation.
- *GSE319641* (NeoTRIP) — baseline and on-treatment tumour RNA-seq, 401 samples. Released fields are patient identifier, timepoint and batch only. No treatment-arm allocation, no pCR status. Response cannot be reconstructed from the deposited metadata.
No public cohort satisfying the requirement was identified. The immunotherapy question is therefore not addressed in this paper.

**Spatial replication.** At the time of analysis no second public breast Visium series with the sample count needed for a perineural distance analysis was available for download within the project's constraints; this is listed as a limitation and as the first priority for follow-up work.

**Consequence.** Two analyses in the original plan — immunotherapy prediction and independent spatial replication — are reported as not performed, with the accession-level reasons above, rather than being silently dropped.

# Supplementary Note 11 — External single-cell anchoring and spot-level deconvolution

## 11.1 Purpose and design

The Visium platform integrates ~10–30 cells per spot, so the cellular provenance of the two NPP modules cannot be read off the spatial data alone. We therefore anchored the modules to an external, publicly available single-cell atlas: **GSE176078** (Wu et al., *Nat Genet* 2021), comprising 26 tumour samples (TNBC, ER-positive/HER2-negative and HER2-positive) plus adjacent normal tissue, ~100,000 cells, with published annotations at major-type resolution and CAF-subtype resolution (myCAF-like s4/s5, MSC/iCAF-like s1/s2, transitioning s3). A secondary Smart-seq2 cohort (**GSE75688**, Chung et al. 2017) served as a directional cross-check for stromal-gene origin.

Three questions were asked, in this order:

1. **Which cell types express the module genes?** (gene-level, in scRNA)
2. **Which cell populations occupy the perineural niche in the Visium sections?** (spot-level seeded-NNLS deconvolution)
3. **Do the deconvolved populations co-aggregate with the module scores and with nerve spots?** (spot-level association statistics identical to the main pipeline)

## 11.2 Gene-level origin (GSE176078)

Per-sample matrices (29,733 genes; MatrixMarket) were parsed, cells assigned to types (major types; CAFs split into myCAF-like / MSC-iCAF-like / transitioning), and aggregated into pseudobulk mean counts-per-10k and per-gene expression fractions per type.

**Perineural-stroma module (60 genes; 59 measurable).** The dominant sources are myCAF-like cells: *SPARC* (98% of myCAF-like cells expressing), *COL1A1* (98%), *COL1A2* (99%), *COL3A1* (98%), *LUM* (98%), *POSTN* (88%), *CTHRC1* (72%), *MMP11* (55%), *VCAN* (48%), *ASPN* (46%), *FAP* (35%), *LRRC15* (12%). Contractile/activation genes are shared with perivascular-like (PVL) cells (*ACTA2* 90% PVL vs 65% myCAF-like, *MYL9, TAGLN, RGS5, MCAM*), *DCN* marks MSC/iCAF-like cells (95%), basement-membrane collagens are co-expressed by endothelium (*COL4A1, COL4A2*), and *VIM*/*S100A4*/*MMP9*/*TGFB1* have substantial myeloid expression. **40 of 59 measurable stromal genes have their maximum expression fraction in myCAF-like cells** (exact counts in Supplementary Table S18).

**Neural-intrinsic module (39 genes; 38 measurable).** The dissociated tumour-cell atlas contains no neuronal, Schwann or glial population. 29 of 38 measurable genes stay below 10% expressing fraction in every annotated type, and 19 below 3%. The handful with higher values map to non-neural sources — *TUBA1A* (75%) and *PMP22* (71%) myCAF-like, *MBP* (33%) myeloid, *UCHL1* (20%) myCAF-like, *S100B* (18%) normal epithelial — consistent with low-level ambient/broad-expression cross-mapping rather than a specific neuronal source; one gene is not reliably measurable on the 10x 3′ platform. This is the expected pattern for transcripts that in tissue are confined to axons and Schwann sheaths: the cells that carry them are rare, fragile and largely lost in tumour dissociation, so their signal is recoverable **in situ** (spatially, next to nerve structures) but not from a dissociated tumour-cell reference. We therefore report the neural-intrinsic module as *in-situ recoverable, single-cell unassignable*, rather than forcing a spurious attribution.

## 11.3 Spot-level deconvolution (seeded NNLS)

A reference panel (union of the top 60 markers per type plus measurable module genes; square-root variance-stabilised pseudobulk profiles) was regressed per spot by non-negative least squares; weights normalised to sum to one give per-spot abundance fractions. Median relative reconstruction error across the 43 Visium samples was 0.408 (IQR 0.392–0.439; 0.397 for the module-excluded panel; Supplementary Figure S9F).

**Primary specification results (all 43 samples).** myCAF-like abundance decreases with distance to the nearest nerve (partial Spearman ρ = −0.062, Fisher-*z* *P* = 3.0 × 10⁻⁴⁸, BH *q* < 10⁻⁴⁶), together with myeloid (−0.075) and endothelial (−0.041) fractions, while cancer-epithelial, T-cell and B-cell abundances increase with distance (ρ = +0.015, +0.018, +0.015, all *q* < 10⁻³) — i.e., nerve spots sit in a myCAF/myeloid-rich, tumour- and lymphocyte-poor microenvironment. Directly contrasting nerve spots with all other spots gives a consistent depletion of tumour-epithelial, B- and T-cell fractions (mean log₂FC −0.51, −0.69, −0.37; Stouffer *P* < 10⁻⁶) and enrichment of myeloid and normal-epithelial fractions (+0.18, +0.17); the myCAF contrast itself is panel-dependent (log₂FC −0.018 primary vs +0.071 module-excluded) and is therefore not used inferentially (see Note 12, caveat N4.1).

**Module co-aggregation.** The perineural-stroma module score is strongly positively associated with deconvolved myCAF abundance (ρ = +0.887, *q* < 10⁻²⁰⁰) and negatively with epithelial and normal-epithelial fractions (−0.432 and −0.540); the neural-intrinsic module is weakly positively associated with myCAF abundance (ρ = +0.138, *q* < 10⁻¹⁰⁰) — consistent with a two-component niche in which nerve-proximal stromal output is myCAF-dominated.

**Circularity control.** Because the primary panel includes module genes, all associations were re-estimated with a **module-excluded panel** (all 99 module genes removed from the reference; Table S19, "module-excluded" rows). myCAF abundance retained — and slightly strengthened — its negative distance gradient (ρ = −0.078), the stroma-module co-aggregation remained strong (ρ = 0.623), and the neural-module association was unchanged (ρ = +0.146). The perineural myCAF gradient is therefore not a restatement of the module definition.

## 11.4 What this anchoring does and does not establish

Established: the spatial perineural program is carried, at the stromal arm, by myCAF-like cells (with PVL and endothelial contributions to specific gene classes); the neural-intrinsic arm is not expressible by any dissociated tumour-cell type and must reflect nerve-associated cells recoverable only in situ.

Not established: sample-matched cell identity (the reference cohort is external); Schwann-cell biology (no Schwann population exists in either public breast scRNA reference, so the neural module's in-situ source cannot be named at single-cell resolution); any statement about ligand–receptor signalling (see Limitation 4).

# NeuroSpatialGradients

A calibrated spatial-transcriptomics framework for measuring the perineural
microenvironment of triple-negative breast cancer.

Companion code for the manuscript *"A calibrated spatial transcriptomics
framework reveals a perineural stromal program in triple-negative breast
cancer"*.

Archived on Zenodo: [10.5281/zenodo.24165453](https://zenodo.org/records/24165453)
(release v1.0.0).

## What it does

1. **Nerve-associated spot calling** from 10x Visium data, using class-wise
   neuronal marker panels with a spatial-cluster constraint, plus a full
   threshold-sensitivity grid.
2. **Distance field and ring modelling** — every spot gets a continuous
   distance to the nearest nerve source; microenvironmental signatures are
   fitted with ring-binned single-exponential decay models.
3. **De novo program derivation** — per-sample partial Spearman gradients
   (adjusted for spot cellularity) meta-analysed by Fisher-*z* and a
   between-sample *t* test.
4. **Moran spectral randomisation (MSR) calibration** — a null that preserves
   each gene's own spatial autocorrelation while breaking its pairing with the
   distance field. Closed-form surrogate statistics in the coefficient domain;
   ~15 CPU-minutes per sample for 10,419 genes × 499 surrogates.
5. **Size-matched random-gene-set calibration** — because a near-deterministically
   smooth distance field correlates with *any* sufficiently large gene set, a
   second threshold is required: 510 random sets at 17 matched sizes.
6. **Functional readouts** — cell-type panels, curated ligand–receptor pairs
   and MSigDB pathways, all passed through the three-threshold calibration
   (MSR, measurability, non-circularity, size-matched baseline).
7. **Expression-derived diffusion pseudotime** rooted at nerve spots — built
   without any spatial information; its agreement with physical distance is a
   non-trivial consistency check.
8. **Module decomposition** — neural-intrinsic versus perineural-stroma
   modules, with mutual adjustment and a marker-removal sensitivity module.
9. **Bulk projection** — METABRIC and GSE12276 survival analyses (Cox,
   Fine–Gray competing risks, bootstrap, E-values, penalised incremental-value
   models).

## Repository layout

```
00_pipeline/nsgrad/     Core library (io, signature, nerve, gradient, surv,
                        metabric, bulk, viz, config)
02_scripts/             Analysis scripts in numbered execution order
03_results/             Tables and figures produced by the pipeline
04_manuscript/          Manuscript source and supplementary notes
05_reports/             Analysis reports
Snakefile               Snakemake workflow (end-to-end)
environment.yml         Pinned conda environment
Dockerfile              Reproducible container
RUN.md                  Runtime and memory budget per step
```

## Quick start

```bash
# Option A: conda
conda env create -f environment.yml
conda activate nsgrad
snakemake --cores 4

# Option B: Docker
docker build -t nsgrad .
docker run -v $(pwd):/work nsgrad snakemake --cores 4
```

All input data are public: GSE210616 (per-sample downloads), METABRIC via
cBioPortal, GSE12276, MSigDB v2023.2. The workflow downloads nothing
automatically; see `RUN.md` for the exact files and URLs expected under
`01_data/raw/`.

## Reproducibility

- Every figure and table regenerates from the raw downloads.
- Random seeds are fixed in `config.py`.
- The rank transform is verified element-wise against
  `scipy.stats.rankdata` on every sample (see `02_scripts/run_step3a_nerve_program.py`).
- Known platform note: `sc.pp.neighbors(..., method="gauss")` is used instead
  of the `umap`/pynndescent default, which aborts on some Windows/numba
  combinations.

## Citation

See `CITATION.cff`. The manuscript's supplementary notes (1–10) document the
calibration logic, the rank-transform correction and its consequences, and the
measurability audits in full.

## License

MIT — see `LICENSE`.

# Runtime and memory budget

Hardware used for the manuscript: Windows 11, 8 physical cores, 32 GB RAM,
Python 3.13. All steps are single-machine; no cluster is required.

## Spatial cohort (GSE210616, 30 samples, ~1,400 spots each)

| Step | Script | Wall time | Peak RSS | Notes |
|---|---|---|---|---|
| QC + nerve calling | `run_step1_gradient.py` | ~6 min | <1 GB | all 12 threshold configs |
| Gene-gradient scan | `run_step3a_nerve_program.py` | ~4 min | <1 GB | 10,419 genes × 30 samples |
| MSR calibration | `run_msr_spatial.py` | **~15 CPU-min per sample** (7.5 h serial / ~1 h at 8 cores) | ~1.5 GB | 499 surrogates × 10,419 genes; closed-form, no field reconstruction |
| Functional readouts scan | `run_niche_spatial.py scan` | ~2 min | <1 GB | 567 readouts |
| Functional readouts meta | `run_niche_spatial.py meta` | <1 min | <1 GB | |
| Random-gene-set scan | `run_niche_spatial.py scan --mode random --nsim 499` | ~2 min | <1 GB | 510 random sets |
| Random null calibration | `run_niche_spatial.py meta_random` + `check_niche_random_calibration.py` | <2 min | <1 GB | |
| Pseudotime | `run_pseudotime.py` | ~30 s per sample | ~1.2 GB | `method="gauss"` neighbours; the `umap`/pynndescent backend aborts on some Windows/numba combos |
| NPP split (spatial) | `run_npp_split.py spatial` | ~3 min | <1 GB | |
| Cross-panel validation | `run_nerve_crossvalidation.py` | ~2 min | <1 GB | |

## Bulk cohorts

| Step | Script | Wall time | Peak RSS |
|---|---|---|---|
| METABRIC Cox + KM | `run_step4_metabric.py` | ~1 min | <0.5 GB |
| Nerve axes + residual | `run_step4b_nerve_axes.py` | ~1 min | <0.5 GB |
| Signature revision | `run_npp_revision.py` | ~2 min | <1 GB |
| GSE12276 Fine–Gray | `run_step5_gse12276.py` | ~2 min | <0.5 GB |
| Bootstrap + E-values (500 reps × 4 models) | `run_sensitivity.py` | ~8 min | <0.5 GB |
| NPP split (bulk + bootstrap 250) | `run_npp_split.py bulk` | ~4 min | <0.5 GB |
| Incremental-value models | `run_prognostic_model.py` | ~2 min | <0.5 GB |

## Figures

Every `make_figure*.py` runs in well under 2 minutes and <1 GB.

## End-to-end

Serial single-machine runtime is dominated by the MSR scan
(~1 h at 8 cores for the spatial cohort). Everything else combined is
under 1 hour.

## Memory notes

- The MSR scan is O(n_genes × n_eigenvectors) in the coefficient domain and
  never reconstructs a full surrogate expression field; peak RSS stays around
  1.5 GB.
- Pseudotime densifies only the readout columns and the 2,000 HVGs, never the
  full 36,601-gene matrix.
- If you rerun the MSR scan with finer eigenvalue bins, expect wall time to
  grow roughly linearly in the number of bins.

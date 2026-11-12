# Snakemake workflow — NeuroSpatialGradients, end to end
# Run: snakemake --cores 4
# All paths are relative to the project root.

configfile: "00_pipeline/nsgrad/config.py" if False else None

rule all:
    input:
        # figures
        "03_results/figures/Fig1B_nervesignal_map.png",
        "03_results/figures/Fig2A_volcano.png",
        "03_results/figures/Fig3_perineural_niche.png",
        "03_results/figures/Fig4_trajectory_and_modules.png",
        "03_results/figures/Fig5AB_km_npp_tertiles.png",
        "03_results/figures/SuppFigS5_spatial_autocorrelation.png",
        "03_results/figures/SuppFigS6_functional_readouts.png",
        # key tables
        "03_results/tables/supplementary/TableS5_MSR_per_gene_statistics.csv",
        "03_results/tables/supplementary/TableS4_readout_measurability_and_circularity.csv",
        "03_results/tables/supplementary/TableS10_bootstrap_and_evalue.csv",


rule step1_gradient:
    input:  scripts=["02_scripts/run_step1_gradient.py"]
    output: "01_data/processed/gradient_params.csv",
            "01_data/processed/ring_profiles.csv"
    shell:  "python 02_scripts/run_step1_gradient.py"

rule step2_subtyping:
    input:  "01_data/processed/gradient_params.csv"
    output: "03_results/tables/consensus_k_selection.csv",
            "03_results/tables/sensitivity_clustering.csv"
    shell:  "python 02_scripts/run_step2_subtyping.py"

rule step3a_nerve_program:
    input:  "01_data/processed/gradient_params.csv"
    output: "01_data/processed/nerve_proximity_program.csv",
            "01_data/processed/npp_signature.json",
            "01_data/processed/npp_per_sample.csv.gz"
    shell:  "python 02_scripts/run_step3a_nerve_program.py"

rule step3b_msr:
    input:  "01_data/processed/npp_signature.json"
    output: "01_data/processed/msr_gene_level.fine.csv",
            "03_results/tables/msr_summary.fine.json"
    shell:  "python 02_scripts/run_msr_spatial.py"

rule step3c_niche:
    input:  "01_data/processed/msr_gene_level.fine.csv"
    output: "01_data/processed/niche_per_sample.csv.gz",
            "03_results/tables/niche_meta_summary.csv"
    shell:  "python 02_scripts/run_niche_spatial.py scan && "
            "python 02_scripts/run_niche_spatial.py meta && "
            "python 02_scripts/run_niche_spatial.py scan --mode random --nsim 499 && "
            "python 02_scripts/run_niche_spatial.py meta_random"

rule step3d_calibration:
    input:  "03_results/tables/niche_meta_summary.csv"
    output: "03_results/tables/supplementary/TableS4_readout_measurability_and_circularity.csv"
    shell:  "python 02_scripts/check_niche_circularity.py && "
            "python 02_scripts/check_niche_random_calibration.py"

rule step4_metabric:
    input:  "01_data/processed/npp_signature.json"
    output: "03_results/tables/step4_metabric_cox.csv"
    shell:  "python 02_scripts/run_step4_metabric.py"

rule step4b_nerve_axes:
    input:  "03_results/tables/step4_metabric_cox.csv"
    output: "03_results/tables/step4b_metabric_nerve_axes.csv"
    shell:  "python 02_scripts/run_step4b_nerve_axes.py"

rule step4c_revision:
    input:  "03_results/tables/step4b_metabric_nerve_axes.csv"
    output: "03_results/tables/npp_revision_metabric.csv"
    shell:  "python 02_scripts/run_npp_revision.py"

rule step5_gse12276:
    input:  "03_results/tables/npp_revision_metabric.csv"
    output: "03_results/tables/step5_gse12276_brainmet.csv"
    shell:  "python 02_scripts/run_step5_gse12276.py"

rule sensitivity:
    input:  "03_results/tables/step5_gse12276_brainmet.csv"
    output: "03_results/tables/sensitivity_metabric_bootstrap.csv"
    shell:  "python 02_scripts/run_sensitivity.py"

rule pseudotime:
    input:  "01_data/processed/niche_per_sample.csv.gz"
    output: "03_results/tables/pseudotime_summary.csv"
    shell:  "python 02_scripts/run_pseudotime.py"

rule npp_split:
    input:  "03_results/tables/pseudotime_summary.csv"
    output: "03_results/tables/npp_split_spatial_meta.csv",
            "03_results/tables/npp_split_bulk_cox.csv"
    shell:  "python 02_scripts/run_npp_split.py spatial && "
            "python 02_scripts/run_npp_split.py bulk"

rule prognostic_model:
    input:  "03_results/tables/npp_split_bulk_cox.csv"
    output: "03_results/tables/prognostic_incremental_value.csv"
    shell:  "python 02_scripts/run_prognostic_model.py"

# ---- figures ----
rule figure1:
    input:  "01_data/processed/gradient_params.csv"
    output: "03_results/figures/Fig1B_nervesignal_map.png"
    shell:  "python 02_scripts/make_figure1.py"

rule figure2:
    input:  "01_data/processed/npp_signature.json"
    output: "03_results/figures/Fig2A_volcano.png"
    shell:  "python 02_scripts/make_figure3.py"

rule figure3:
    input:  "03_results/tables/niche_meta_summary.csv"
    output: "03_results/figures/Fig3_perineural_niche.png"
    shell:  "python 02_scripts/make_figure_niche.py"

rule figure4:
    input:  "03_results/tables/pseudotime_summary.csv"
    output: "03_results/figures/Fig4_trajectory_and_modules.png"
    shell:  "python 02_scripts/make_figure_trajectory.py"

rule figure5:
    input:  "03_results/tables/step4_metabric_cox.csv"
    output: "03_results/figures/Fig5AB_km_npp_tertiles.png"
    shell:  "python 02_scripts/make_figure4.py"

rule figure_supp:
    input:  "03_results/tables/msr_summary.fine.json"
    output: "03_results/figures/SuppFigS5_spatial_autocorrelation.png",
            "03_results/figures/SuppFigS6_functional_readouts.png"
    shell:  "python 02_scripts/make_figure_msr.py && "
            "python 02_scripts/make_figure_supp_s6.py"

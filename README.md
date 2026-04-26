# YELLOW RoUTIne / rUTI Minimal Analysis Repository

This repository is intentionally trimmed to the **current canonical analysis workflow** anchored by `RUN_COMPLETE_ANALYSIS.sh`.

It keeps only:
- the numbered pipeline scripts required for the current run order,
- helper code required by those scripts,
- minimal configuration/reference files,
- documentation needed to run and understand the workflow.

It excludes historical UTIGA/poster/conference material, legacy/proposed scripts, and generated outputs.

## Canonical pipeline entrypoint

```bash
bash RUN_COMPLETE_ANALYSIS.sh
```

## Script order used by the pipeline

### Phase 0
- `00a_load_clean_clinical.R`
- `00b_classify_episodes.R`
- `00c_plot_clinical_summary.R`

### Phase 1
- `12a_wgs_qc.R`
- `12b_core_snp.R`
- `12c_panaroo.R`
- `13_visualise_panaroo_selection.R`
- `02_gene_presence_analysis.R`
- `06_MLST.R`

### Phase 1b
- `03_plotting.R`
- `04_gene_breakdown.R`
- `05_gene_overview_plots.R`
- `07_explore_MLST.R`
- `08_core_vs_plasmid.R`
- `09_inc_plasmid_network.R`
- `10_replicon_heatmap.R`

### Phase 2
- `11_compare_strains.R`
- `14_genotype_phenotype_model.R`
- `17_lineage_analysis.R`

### Phase 3
- `15_longitudinal_patterns.R`
- `16_within_host_evolution.R`
- `18_annotate_variants.R`
- `20_variant_annotation_deep.R`
- `19_host_context.R`
- `21_publication_figures.R`

### Phase 4
- `22_vf_build_analysis_dataset.R`
- `23_vf_cross_sectional.R`
- `24_vf_longitudinal_dynamics.R`
- `25_vf_lineage_vf_interaction.R`

## Required software

### Conda / system tools
Install tools used by the pipeline, including:
- `abricate`
- `mlst`
- `parsnp`
- `panaroo`
- `snp-dists`
- `prokka`
- `mummer`

Use your environment manager of choice (conda/mamba recommended).

### R runtime
R plus packages used by the numbered scripts and helper files (e.g. tidyverse/data.table/ggplot2 stack and modeling/network packages where referenced).

## Inputs not committed

Large/sensitive/raw data are intentionally not committed. You are expected to provide local inputs such as:
- assembly FASTA/annotation files,
- raw/cleaned clinical input tables,
- tool-specific reference inputs.

See script path constants in `00_config.R` and the cleanup manifest for scope.

## Generated outputs (not committed)

Pipeline outputs are expected under directories such as:
- `results/`
- `plots/`
- `logs/`
- tool output folders (e.g. panaroo/parsnp/mlst/abricate outputs)

These are ignored via `.gitignore` to keep the repository analysis-only and reviewable.

## Cleanup manifest

Repository keep/remove rationale and cleanup rules are documented in:
- `REPO_CLEANUP_MANIFEST.md`

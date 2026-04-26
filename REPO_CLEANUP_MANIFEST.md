# REPO_CLEANUP_MANIFEST

Date: 2026-04-26
Branch: `work`
Anchor workflow: `RUN_COMPLETE_ANALYSIS.sh`

## Scope
This repository is trimmed to the minimal YELLOW RoUTIne / rUTI analysis workflow.

## Kept (required)

### 1) Canonical runner
- `RUN_COMPLETE_ANALYSIS.sh`

### 2) Canonical numbered scripts called by runner
- `00a_load_clean_clinical.R`
- `00b_classify_episodes.R`
- `00c_plot_clinical_summary.R`
- `12a_wgs_qc.R`
- `12b_core_snp.R`
- `12c_panaroo.R`
- `13_visualise_panaroo_selection.R`
- `02_gene_presence_analysis.R`
- `06_MLST.R`
- `03_plotting.R`
- `04_gene_breakdown.R`
- `05_gene_overview_plots.R`
- `07_explore_MLST.R`
- `08_core_vs_plasmid.R`
- `09_inc_plasmid_network.R`
- `10_replicon_heatmap.R`
- `11_compare_strains.R`
- `14_genotype_phenotype_model.R`
- `17_lineage_analysis.R`
- `15_longitudinal_patterns.R`
- `16_within_host_evolution.R`
- `18_annotate_variants.R`
- `20_variant_annotation_deep.R`
- `19_host_context.R`
- `21_publication_figures.R`
- `22_vf_build_analysis_dataset.R`
- `23_vf_cross_sectional.R`
- `24_vf_longitudinal_dynamics.R`
- `25_vf_lineage_vf_interaction.R`

### 3) Directly sourced helpers
- `00_config.R`
- `11_compare_strains_helpers.R`
- `R/clinical_helpers.R`
- `R/plot_helpers.R`
- `R/wgs_helpers.R`

### 4) Required workflow docs/config
- `README.md`
- `.gitignore`
- `REPO_CLEANUP_MANIFEST.md`
- Environment files: `environment-rutis-core.yml`, `env-wgs.yml`, `env-annot.yml`

### 5) Required small input file(s)
- `assembly_metadata.csv` (used by canonical scripts)

## Removed
- UTIGA/poster/conference/slides material
- generated outputs (`plots/`, result artifacts)
- legacy/proposed/debug/backup/one-off scripts
- non-canonical helper/doc files not required to run the canonical pipeline

## Validation checklist
- `bash -n RUN_COMPLETE_ANALYSIS.sh`
- Confirm every `Rscript` target in `RUN_COMPLETE_ANALYSIS.sh` exists
- Confirm `source(...)` paths in canonical scripts exist
- Parse-check canonical `.R` scripts with `Rscript -e "parse(file='...')"` when `Rscript` is available

## Sensitive data policy
If private/raw clinical or sequencing data is found, keep it out of the working tree, ignore it in Git, and note that old history may still contain prior commits.

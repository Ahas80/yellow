# REPO_CLEANUP_MANIFEST

Date: 2026-04-25
Branch: `work`
Anchor workflow: `RUN_COMPLETE_ANALYSIS.sh`

## Cleanup goal
Retain only files needed for the current YELLOW RoUTIne / rUTI canonical analysis pipeline and remove UTIGA/poster/legacy/generated clutter.

## Canonical scripts called directly by `RUN_COMPLETE_ANALYSIS.sh`

- Phase 0: `00a_load_clean_clinical.R`, `00b_classify_episodes.R`, `00c_plot_clinical_summary.R`
- Phase 1: `12a_wgs_qc.R`, `12b_core_snp.R`, `12c_panaroo.R`, `13_visualise_panaroo_selection.R`, `02_gene_presence_analysis.R`, `06_MLST.R`
- Phase 1b: `03_plotting.R`, `04_gene_breakdown.R`, `05_gene_overview_plots.R`, `07_explore_MLST.R`, `08_core_vs_plasmid.R`, `09_inc_plasmid_network.R`, `10_replicon_heatmap.R`
- Phase 2: `11_compare_strains.R`, `14_genotype_phenotype_model.R`, `17_lineage_analysis.R`
- Phase 3: `15_longitudinal_patterns.R`, `16_within_host_evolution.R`, `18_annotate_variants.R`, `20_variant_annotation_deep.R`, `19_host_context.R`, `21_publication_figures.R`
- Phase 4: `22_vf_build_analysis_dataset.R`, `23_vf_cross_sectional.R`, `24_vf_longitudinal_dynamics.R`, `25_vf_lineage_vf_interaction.R`

## Keep / remove policy

### Keep
1. `RUN_COMPLETE_ANALYSIS.sh`
2. Numbered canonical pipeline scripts listed above
3. Direct helper files sourced by kept scripts (e.g. selected files in `R/`)
4. Required environment and user documentation (`README.md`, `.gitignore`, env YAMLs)
5. Small static lookup/reference files required as inputs

### Remove
1. Generated outputs (`plots/`, `results/`, `logs/`, caches/intermediates)
2. Poster/conference/slides and presentation artifacts
3. UTIGA-specific material
4. Legacy/proposed/backup/debug/one-off scripts not in canonical pipeline
5. Duplicate drafts superseded by kept numbered scripts

### Needs decision (if encountered)
Any file not directly called/sourced/read by canonical scripts and not clearly generated output should be marked as uncertain in review notes instead of silently retained.

## Safety checks planned
- `bash -n RUN_COMPLETE_ANALYSIS.sh`
- Confirm every `Rscript` target in `RUN_COMPLETE_ANALYSIS.sh` exists
- Parse-check kept `.R` scripts (`Rscript -e "parse(file='...')"`) if R is available
- Verify `source("...")` references resolve to existing files

## Notes on sensitive/raw data
If private or raw sequencing/clinical data appears during cleanup, remove from working tree, block via `.gitignore`, and flag that it may still exist in Git history.

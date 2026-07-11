# Backup Branch Manifest

## Backup metadata

- Backup created: 2026-05-15 18:10:16 CEST
- Current Git branch before backup: `codex-source-only-github-update`
- Source commit inspected before backup: `d14d0e5`
- New backup branch: `backup/core-pipeline-20260515`
- Backup style: source-only orphan branch built from an isolated temporary Git index, so local generated files were not deleted or overwritten.

## Included files

### Numbered R pipeline scripts

- `00_config.R`
- `00_input_snapshot.R`
- `00_make_assembly_metadata.r`
- `00a_load_clean_clinical.R`
- `00b_classify_episodes.R`
- `00c_plot_clinical_summary.R`
- `00d_derive_plot_timepoints.R`
- `02_gene_presence_analysis.R`
- `03_plotting.R`
- `04_gene_breakdown.R`
- `05_gene_overview_plots.R`
- `06_MLST.R`
- `07_explore_MLST.R`
- `08_core_vs_plasmid.R`
- `09_inc_plasmid_network.R`
- `10_replicon_heatmap.R`
- `11_compare_strains.R`
- `11_compare_strains_helpers.R`
- `12a_wgs_qc.R`
- `12b_core_snp.R`
- `12c_panaroo.R`
- `13_visualise_panaroo_selection.R`
- `14_genotype_phenotype_model.R`
- `15_longitudinal_patterns.R`
- `16_within_host_evolution.R`
- `17_lineage_analysis.R`
- `18_annotate_variants.R`
- `19_host_context.R`
- `20_variant_annotation_deep.R`
- `21_publication_figures.R`
- `22_vf_build_analysis_dataset.R`
- `23_vf_cross_sectional.R`
- `24_vf_longitudinal_dynamics.R`
- `25_vf_lineage_vf_interaction.R`
- `26_vf_define_gene_modules.R`
- `27_vf_score_framework.R`
- `28_vf_transition_case_studies.R`
- `29_vf_amr_combined_profile.R`
- `30_vf_project_summary_tables.R`
- `31_audit_fasta_usage.R`
- `31_audit_uti_denominator_drop.R`

### Required sourced R helpers

- `R/clinical_helpers.R`
- `R/pipeline_qc_helpers.R`
- `R/plot_helpers.R`
- `R/wgs_helpers.R`

### Runner, documentation, and environment files

- `.gitignore`
- `BACKUP_BRANCH_MANIFEST.md`
- `FOLDER_MAP.md`
- `README.md`
- `RUN_COMPLETE_ANALYSIS.sh`
- `env-annot.yml`
- `env-wgs.yml`
- `environment-rutis-core.yml`

## Intentionally excluded files and folders

- Raw sequencing and assembly data: `ont-yellow-routine-fastas/`, `data/`, `assemblies.list`, `*.fasta`, `*.fa`, `*.fna`, `*.fastq`, `*.fq`, `*.fastq.gz`, `*.fq.gz`
- Genomics tool outputs: `*.bam`, `*.sam`, `*.vcf`, `*.vcf.gz`, `*.gff`, `*.gff3`, `*.gbk`, `*.gbff`
- Generated outputs: `results/`, `plots/`, `logs/`, `Rplots.pdf`
- Archived and legacy generated material: `archive/`, `legacy/unused_modules/`, `not-assemblies/`
- Binary papers/slides and poster outputs: `docs/papers/`, `docs/slides/`, `UTIGA/`
- Local or temporary files: `.DS_Store`, `.Rhistory`, `.RData`, `.Rproj.user/`, `tmp/`, `temp/`, `cache/`, `.snakemake/`
- Generated data tables that are required at runtime but should be rebuilt or supplied locally: `assembly_metadata.csv`, `status_map.csv`, `status_map_full.csv`, `mlst_all.tsv`, `vf_pa_all.csv`
- Convenience or optional non-numbered scripts not required by the current runner, including `RUN_IN_TERMINAL.sh`, `scripts/wgs_audit.sh`, `scripts/12d_wgs_badsize_rescue_screen.R`, `scripts/32_compare_primary_vs_rescue_vf.R`, `scripts/audit_vf_ready_exclusions_standalone.R`, and `scripts/regenerate_missing_gffs.R`

## Sensitive or unusually large files found

- `OVERVIEW E.coli batch 1-6 - CLEAN_RC.xlsx` is a required local clinical input but was left uncommitted because it may contain participant-level clinical data.
- Large generated genomics outputs were found under `results/`, including FASTA, GFF, VCF, Panaroo, Parsnp, and strain-comparison artifacts.
- Many generated figure files were found under `plots/`.
- Archived generated files and installers were found under `archive/cleanup_2026-05-14/` and were intentionally excluded.

## Required local files that may be missing after restore

This branch is source-only. To run the full workflow after cloning or restoring it, provide the local input data and generated intermediates expected by `00_config.R` and `RUN_COMPLETE_ANALYSIS.sh`, especially:

- `OVERVIEW E.coli batch 1-6 - CLEAN_RC.xlsx`
- raw assembly/sequence files under the local data paths documented in the project
- generated intermediate tables in `results/` when resuming instead of rebuilding

## How to restore or use this branch

1. Fetch and switch to the branch: `git fetch origin backup/core-pipeline-20260515 && git switch backup/core-pipeline-20260515`
2. Place the required local input data back into the documented local paths.
3. Create or activate the conda/R environment from the included environment YAML files.
4. Run `bash RUN_COMPLETE_ANALYSIS.sh` from the repository root.

This branch is intended as a safe GitHub backup of the code and workflow surface, not a complete archive of local biological data or generated results.

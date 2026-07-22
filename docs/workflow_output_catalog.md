# Workflow Output, Figure, and Statistics Catalog

This companion to `workflow_flowchart.md` catalogs what the runner scripts
produce and what statistical or visual evidence each step adds. The order follows
`RUN_COMPLETE_ANALYSIS.sh`.

Use this file when you need to answer practical questions such as:

- Which script produced this plot or table?
- Which scripts run formal or exploratory statistical tests?
- Which outputs are denominator/QC checks rather than biological findings?
- Which figures belong to the final manuscript-style figure pack?

## Source of Truth

- Run order: `RUN_COMPLETE_ANALYSIS.sh`.
- Script behavior: headers, `source(...)` calls, read/write calls, model/test
  calls, and plotting calls in each runner script.
- Figure metadata where available:
  - `results/vf/vf_figure_index.csv`
  - `results/vf/vf_visualisation_audit.csv`
  - `results/vf/uti_not_uti_diagnostic_figure_metadata.csv`
  - `results/statistical_sensitivity/statistical_sensitivity_figure_metadata.csv`
  - `results/final_figures/final_figure_manifest.csv`
  - `results/summary/table_16_uti_not_uti_diagnostic_figures.csv`

## Important Interpretation Guardrails

- The primary clinical contrast is `UTI` vs `Not_UTI`. Older ASB-vs-UTI outputs
  are legacy unless explicitly archived or labelled as old.
- The VF-ready UTI denominator is sparse. Treat broad screens as exploratory
  unless a script explicitly models participant structure and still carries the
  sparse-data caveat.
- Repeated isolates from the same participant are not independent unless the
  method explicitly accounts for participant clustering or collapses within
  participant.
- Sequence type is lineage context. ST agreement does not prove same-strain
  persistence; SNP evidence is primary where available.
- Display-only Uricult/poster timepoints are for visualization and explanation,
  not statistical covariates.
- VF gene presence is genomic detection, not expression, activity, or causality.

## Phase 0: Clinical Data Foundation

### `00a_load_clean_clinical.R`

- Role: raw clinical merge and harmonization.
- Reads: clinical batch CSVs configured in `00_config.R`.
- Main work: standardizes column names and types, detects signs-and-symptoms
  fields, normalizes participant and timepoint identifiers, and preserves missing
  batch warnings.
- Tables/reports:
  - `results/clinical/uti_detected_column_map.csv`
  - `results/clinical/intermediate/clinical_merged.rds`
- Statistics: none.
- Images: none.
- Downstream: `00b_classify_episodes.R` classifies episodes from the merged RDS.

### `00b_classify_episodes.R`

- Role: primary UTI vs Not_UTI episode classification.
- Reads: `results/clinical/intermediate/clinical_merged.rds`.
- Main work: applies catheter-aware symptom rules, urine culture support,
  CFU-threshold logic, manual curation fields, duplicate handling, and
  legacy-status comparison.
- Tables/reports:
  - `results/clinical/status_map.csv`
  - `results/clinical/status_map_legacy_comparison.csv`
  - `results/clinical/uti_binary_classification_audit.csv`
  - `results/clinical/uti_reclassification_movement_table.csv`
  - `results/clinical/uti_symptom_rule_audit.csv`
  - `results/clinical/uti_cfu_threshold_audit.csv`
  - UTI attrition outputs written through shared QC helpers.
- Statistics: rule counts and audit counts only.
- Images: none.
- Downstream: every clinical, WGS, VF, modelling, longitudinal, and validation
  step uses `status_map.csv`.

### `00d_derive_plot_timepoints.R`

- Role: display-only Uricult/event timepoint placement.
- Reads: `status_map.csv` and raw batch collection dates.
- Main work: recovers collection dates, orders routine visits, places Uricult
  rows between routine visits where possible, and records placement confidence.
- Tables/reports:
  - `results/clinical/status_map_with_poster_tp.csv`
  - `results/clinical/unplaced_uricult_rows.csv`
  - `results/clinical/display_only_uricult_rows.csv`
  - Optional `UTIGA/data/status_map_with_poster_tp.csv` if that folder exists.
- Statistics: none.
- Images: none.
- Guardrail: poster/display timepoints explain chronology visually only.

### `00c_plot_clinical_summary.R`

- Role: clinical denominator visualization and QA.
- Reads: `status_map.csv` and poster timepoint map if present.
- Main work: summarizes status counts, assembly availability by status,
  reclassification, Not_UTI subgroups, symptom-rule provenance, CFU-threshold
  provenance, and episode inspection rows.
- Tables/reports:
  - `results/clinical/assembly_metrics_by_status.csv`
  - `results/clinical/waterfall_counts.csv`
  - `results/clinical/uti_reclassification_plot_table.csv`
  - `results/clinical/not_uti_subgroup_by_batch_event.csv`
  - `results/clinical/uti_symptom_rule_plot_table.csv`
  - `results/clinical/uti_cfu_threshold_plot_table.csv`
  - `results/clinical/uti_episode_inspection_table.csv`
- Statistics: descriptive counts only.
- Images:
  - `plots/clinical/trajectories_heatmap.png`
  - `plots/clinical/transitions_alluvial_or_heatmap.png`
  - `plots/clinical/assembly_contigs_boxplot.png`
  - `plots/clinical/waterfall_counts.png`
  - `plots/clinical/uti_reclassification_heatmap.png`
  - `plots/clinical/not_uti_subgroup_by_batch_event.png`
  - `plots/clinical/uti_symptom_rule_provenance.png`
  - `plots/clinical/uti_cfu_threshold_provenance.png`

## Phase 1: WGS Processing

### `12a_wgs_qc.R`

- Role: assembly QC and canonical assembly selection.
- Reads: `assembly_metadata.csv`, FASTA assemblies, and `status_map.csv` when
  available.
- Main work: computes total length, contig count, N50, GC percent, read errors,
  QC pass/fail status, manual-curation exclusions, and one selected canonical
  assembly per participant-timepoint.
- Tables/reports:
  - `results/wgs/qc_summary.csv`
  - `results/qc/canonical_assembly_selection.csv`
  - `results/qc/wgs_qc_manual_curation_excluded_rows.csv`
  - `results/qc/qc_selection_bias_by_status.csv`
  - `results/qc/qc_selection_bias_report.txt`
  - Denominator summaries written by shared QC helpers.
- Statistics: Fisher exact QC-selection bias check by primary status when the
  status table is valid.
- Images:
  - `plots/wgs/wgs_qc_n50_vs_contigs.png`
- Guardrail: genomic comparisons are conditional on QC-passing WGS availability.

### `12b_core_snp.R`

- Role: core SNP alignment and strain-distance generation.
- Reads: canonical assembly selection or `results/wgs/qc_summary.csv`, selected
  FASTA assemblies, and external tools `parsnp`, `harvesttools`, and `snp-dists`.
- Main work: fingerprints current FASTA inputs, writes an input manifest and
  staleness report, stages FASTAs, runs Parsnp, converts GGR/XMFA to FASTA when
  needed, computes pairwise SNP distances, builds a neighbor-joining tree, and
  classifies pairs as Same, Related, or Different by SNP thresholds.
- Tables/reports:
  - `results/wgs/core/core_snp_input_manifest.csv`
  - `results/wgs/core/core_snp_staleness_report.txt`
  - `results/wgs/core/core_snp_input_manifest.hash`
  - `results/wgs/core/snp_dists.tsv`
  - `results/wgs/core/strain_pairs.csv`
- Other artifacts:
  - `results/wgs/core/core_genome.tree`
  - `results/wgs/core/parsnp_out/`
- Statistics: pairwise SNP distances and rule-based strain calls, not an
  association test.
- Images: none.

### `12c_panaroo.R`

- Role: GFF preparation and Panaroo pangenome generation.
- Reads: canonical QC PASS assemblies, FASTA/GFF inventories, Prokka folders,
  and Panaroo.
- Main work: builds assembly/GFF inventory, regenerates missing GFFs, refuses
  stale metadata/QC inventories, fingerprints Panaroo inputs, runs or skips
  Panaroo based on manifest freshness, and records completeness.
- Tables/reports:
  - `results/wgs/pan/panaroo_input_manifest.csv`
  - `results/wgs/pan/missing_gffs.csv`
  - `results/wgs/pan/regenerate_missing_gffs_summary.csv`
  - `results/wgs/pan/panaroo_staleness_report.txt`
  - `results/wgs/pan/panaroo_input_manifest.hash`
  - `results/wgs/pan/gene_presence_absence.csv`
- Statistics: none.
- Images: none.
- Guardrail: if the selected assembly set changed, Panaroo must be refreshed
  before interpreting VF or accessory-gene outputs.

### `13_visualise_panaroo_selection.R`

- Role: WGS/Panaroo selection visibility.
- Reads: QC summary, canonical selection, Panaroo manifest, assembly metadata,
  and status map.
- Main work: summarizes selected, missing, excluded, and available assemblies by
  participant and timepoint, then checks selection bias by status.
- Tables/reports:
  - `results/wgs/panaroo_selection_summary.csv`
  - `results/wgs/panaroo_selection_detailed.csv`
  - `results/wgs/panaroo_timepoint_selection_summary.csv`
  - `results/qc/qc_selection_bias_by_status.csv`
  - `results/qc/qc_selection_bias_report.txt` when Fisher testing can be
    computed.
  - Denominator summary rows written by shared QC helpers.
- Statistics: Fisher exact selection-bias check by primary status when valid.
- Images:
  - `plots/wgs/panaroo_selection_matrix.png`
  - `plots/wgs/panaroo_selection_overview.png`
  - `plots/wgs/panaroo_selection_samples_per_timepoint.png`
  - `plots/wgs/panaroo_timepoint_selection_bar.png`
  - `plots/wgs/panaroo_timepoint_selection_tile.png`
- Guardrail: this explains selection; it does not change inclusion rules.

### `02_gene_presence_analysis.R`

- Role: VF gene hit parsing and canonical presence/absence matrix.
- Reads: Panaroo or VF detection outputs, ABRicate VFDB-style hits when present,
  selected assembly information, and optional selection arguments.
- Main work: applies identity/coverage filters, removes manually curated rows,
  creates long VF hit tables, constructs isolate-by-gene binary matrices, and
  summarizes gene prevalence.
- Tables/reports:
  - `results/qc/vf_abricate_manual_curation_excluded_rows.csv`
  - `results/vf/vf_hits_all.rds`
  - `results/vf/vf_hits_all.csv`
  - `results/vf/vf_pa_all.csv`
  - `results/vf/stats_gene_level.csv`
  - Core and variable gene lists.
- Statistics: descriptive gene prevalence and burden only.
- Images:
  - `plots/vf/core_bar_top25_all.png`
  - `plots/vf/core_histogram_all.png`
- Guardrail: VF presence/absence does not measure expression or virulence
  activity.

### `06_MLST.R`

- Role: active provider-preferred MLST integration.
- Reads: provider/RIVM SeqSphere MLST source, local MLST provenance if needed,
  and MLST helper scripts.
- Main work: runs or reuses local MLST provenance, compares provider vs local
  calls, builds provider-preferred ST outputs, and verifies downstream scripts
  use the active MLST source.
- Tables/reports:
  - `results/mlst/mlst_provider_preferred.csv`
  - `results/mlst/mlst_provider_preferred_all.csv`
  - `results/mlst/mlst_provider_source_audit.csv`
  - `results/mlst/mlst_all.tsv`
  - `results/mlst/mlst_with_meta.csv`
  - Source-usage verification outputs from helper scripts.
- Statistics: source counts and provenance audits only.
- Images: none.
- Guardrail: ST is lineage context, not same-strain proof.

## Phase 1b: Additional Plots and Exploration

### `03_plotting.R` optional

- Role: legacy exploratory plotting, skipped unless
  `RUN_LEGACY_EXPLORATORY_PLOTS=1`.
- Reads: `vf_pa_all.csv`, gene stats, status map, MLST, pairwise metrics,
  plasmid/nitrate files when present.
- Main work: writes older descriptive stats, persistence plots, legacy status
  plots, phylogeny/timeline/network views, and nitrate heatmap if data exist.
- Tables/reports:
  - `results/stats_sample_level.csv`
  - `results/stats_participant_level.csv`
  - `results/stats_cohort_level.csv`
  - `results/stats_descriptive.md`
  - `results/persistence_summary.csv`
  - `results/stats_vf_burden_wilcox.txt`
- Statistics: descriptive summaries and optional Wilcoxon text output; legacy
  Fisher/volcano material is not the current primary analysis.
- Images:
  - `plots/core_bar_top25_all.png`
  - `plots/core_histogram_all.png`
  - `plots/richness_by_timepoint.png`
  - `plots/richness_trajectories_numeric.png`
  - `plots/persistence/upset_P*.png`
  - `plots/upset_genes_by_status.png`
  - `plots/legacy/old_asb_uti_outputs/volcano_UTI_vs_ASB_legacy.png`
  - `plots/phylogeny/core_tree_phenotype.png`
  - `plots/phylogeny/core_tree_base.png`
  - `plots/epidemiology/st_distribution_stacked.png`
  - `plots/genomics/snp_distance_violin.png`
  - `plots/timelines/swimmer_plot_top20.png`
  - `plots/epidemiology/transmission_network.png`
  - `plots/genomics/virulence_nitrate_heatmap.png`
- Guardrail: optional and legacy. Current canonical UTI vs Not_UTI VF figures
  come from scripts `23` through `35`.

### `04_gene_breakdown.R`

- Role: focused gene annotation and category summaries.
- Reads: `vf_pa_all.csv`, `gene_map.csv` if present, status map, and VF helper
  definitions.
- Main work: creates/updates gene map, annotates genes, counts per-sample
  categories, extracts nitrate genes, and optionally runs focus-gene models.
- Tables/reports:
  - `results/vf/gene_map.csv`
  - `results/vf/annotated_gene_table.csv`
  - `results/vf/per_sample_category_counts.csv`
  - `results/vf/nitrate_presence_matrix.csv`
  - `results/vf/diff_focus_genes_UTI_vs_Not_UTI_glmm.csv`
- Statistics: optional focus-gene GLMM with participant random intercept,
  GLM fallback, and BH-adjusted p-values.
- Images:
  - `plots/vf/nitrate_upset.png`
- Guardrail: focus-gene modelling is narrow and should not replace the main
  genotype-phenotype model in `14`.

### `05_gene_overview_plots.R`

- Role: descriptive VF gene overview.
- Reads: `vf_pa_all.csv`.
- Main work: identifies core and variable genes, ranks top prevalent genes, and
  plots the variable-gene matrix.
- Tables/reports: no new primary tables.
- Statistics: descriptive prevalence only.
- Images:
  - `plots/vf/gene_prevalence_bar.png`
  - `plots/vf/variable_gene_heatmap.png`
  - `plots/vf/variable_gene_heatmap.pdf`

### `07_explore_MLST.R`

- Role: exploratory MLST distribution.
- Reads: active MLST output and assembly metadata.
- Main work: counts ST frequencies, joins metadata when possible, and checks
  duplicate metadata rows.
- Tables/reports:
  - `results/mlst/ST_frequencies.csv`
  - `results/debug/meta_duplicates.csv` when duplicates are detected.
  - Exploratory MLST-with-metadata output configured by `00_config.R`.
- Statistics: descriptive ST counts only.
- Images:
  - `plots/mlst/top20_STs.png`
  - `plots/mlst/top20_STs.pdf`

### `08_core_vs_plasmid.R`

- Role: ST and plasmid/pMLST association exploration.
- Reads: provider-preferred MLST, canonical assembly selection, external `mlst`,
  and PlasmidFinder/pMLST tools.
- Main work: writes chromosomal ST frequencies, checks available plasmid MLST
  schemes, runs pMLST where possible, falls back to PlasmidFinder replicons, and
  tests ST-plasmid associations.
- Tables/reports:
  - `results/mlst/ST_core_freq.csv`
  - `results/mlst/pmlst_logs/mlst_list.txt`
  - `results/mlst/pMLST_hits_long.csv`
  - `results/mlst/plasmid_types_per_isolate.csv`
  - `results/mlst/plasmid_replicons_long.csv`
  - `results/mlst/plasmid_replicons_wide.csv`
  - `results/mlst/ST_plasmid_associations.csv`
- Statistics: Fisher-style ST-plasmid association tests where valid.
- Images: none in this script.

### `09_inc_plasmid_network.R`

- Role: PlasmidFinder replicon matrix and network plots.
- Reads: canonical assembly selection and provider-preferred MLST.
- Main work: runs ABRicate PlasmidFinder, caches tabular outputs, builds
  replicon long and presence/absence tables, and plots replicon networks.
- Tables/reports:
  - `results/plasmids/plasmidfinder_input_manifest.csv`
  - `results/plasmids/plasmidfinder_hits_long.csv`
  - `results/plasmids/plasmidfinder_presence_absence.csv`
  - `results/plasmids/abricate_cache/`
- Statistics: descriptive co-occurrence only.
- Images:
  - `plots/plasmids/replicon_cooccurrence.pdf`
  - `plots/plasmids/ST_vs_replicon_network.pdf`
- Downstream: `10`, `11`, `14`, `29`, and final summaries can use replicon
  presence/absence.

### `10_replicon_heatmap.R`

- Role: plasmid replicon heatmap.
- Reads: plasmid replicon matrix and MLST/status annotations.
- Main work: orders isolates and replicons, writes ordering debug files, and
  renders heatmaps.
- Tables/reports:
  - `logs/debug/replicons_kept_order.txt`
  - `logs/debug/isolates_order.txt`
- Statistics: none.
- Images:
  - `plots/plasmids/replicon_heatmap.png`
  - `plots/plasmids/replicon_heatmap.pdf`

## Phase 2: Comparative Genomics

### `11_compare_strains.R --participants ALL`

- Role: within-participant strain comparison.
- Reads: status map, SNP distances, VF matrix, plasmid/replicon data, MLST,
  assembly metadata, and helper functions.
- Main work: builds within-participant pairs, calculates ANI/SNP/VF/plasmid
  similarities, applies SNP-threshold same-strain rules, and creates summary
  plots.
- Tables/reports:
  - `results/strain_compare/pairwise_metrics.csv`
  - `results/strain_compare/summary_counts.csv`
  - `results/strain_compare/summary_by_participant.csv`
  - `results/strain_compare/stats_within_vs_between.csv`
  - `results/strain_compare/stats_by_status.csv`
  - `results/strain_compare/README.txt`
- Statistics: Kruskal-Wallis summaries for selected metrics when enough data
  exist.
- Images:
  - `plots/strain_compare/heatmap_vf_jaccard.png`
  - `plots/strain_compare/heatmap_inc_jaccard.png`
  - `plots/strain_compare/identity_vs_snps_scatter.png`
  - `plots/strain_compare/snp_distance_violin.png`
  - `plots/strain_compare/network_same_strain.png`
  - `plots/strain_compare/timeline_by_participant.png`

### `22_vf_build_analysis_dataset.R`

- Role: canonical VF-ready analysis dataset.
- Reads: `vf_pa_all.csv`, status map, QC/selection files, MLST, manual curation
  fields, and bridge inputs.
- Main work: validates anchor files, anti-joins VF and status rows, bridges
  Uricult clinical episodes to UTI-N WGS rows, filters primary/manual-curated
  rows, checks duplicates, adds ST and VF burden, and writes dataset diagnostics.
- Tables/reports:
  - `results/qc/mlst_duplicate_participant_timepoint_st_conflicts.csv`
  - `results/vf/vf_without_status_rows.csv`
  - `results/vf/status_without_vf_rows.csv`
  - `results/vf/vf_status_unmatched_rows.csv`
  - `results/qc/uricult_bridge_audit.csv`
  - `results/qc/uricult_bridge_sensitivity_alternatives.csv`
  - `results/qc/uricult_bridge_key_conflicts.csv`
  - `results/qc/vf_ready_manual_curation_excluded_rows.csv`
  - Duplicate-key diagnostics for VF-ready rows when present.
  - `results/vf/vf_gene_annotation_gap_report.csv`
  - `results/vf/vf_gene_annotation_gap_report.txt`
  - `results/vf/uti_vf_episode_inspection.csv`
  - `results/vf/vf_analysis_ready.csv`
  - `results/vf/vf_binary_uti_ready.csv`
  - `results/vf/vf_dataset_diagnostics.txt`
- Statistics: denominator and QC counts only.
- Images: none.
- Guardrail: downstream VF scripts should use this file instead of rejoining raw
  clinical and VF inputs.

### `14_genotype_phenotype_model.R`

- Role: main genotype-phenotype association modelling.
- Reads: status map, VF features, plasmid replicons, active MLST, assembly
  metadata, and `vf_binary_uti_ready.csv` when available.
- Main work: builds a model denominator, filters features by prevalence, runs
  Fisher screening, fits participant-aware GLMMs, falls back to GLM as needed,
  adjusts p-values with BH FDR, writes sparse/separation warnings, and optionally
  writes random forest importance.
- Tables/reports:
  - `results/models/model_dataset_denominator.csv`
  - `results/models/gwas_univariable_stats.csv`
  - `results/models/gwas_multivariable_glmm.csv`
  - `results/models/model_interpretation_warnings.txt`
  - Primary binary model warnings file configured in `00_config.R`
  - `results/models/rf_variable_importance.csv` when the optional RF branch runs.
- Statistics: Fisher exact feature screening, GLMM/GLM odds ratios, confidence
  intervals, BH FDR, convergence flags, sparse/separation flags.
- Images:
  - `plots/volcano_plot_UTI_vs_Not_UTI.png`
  - `plots/forest_plot_top_hits.png`
  - `plots/vf/vf_gene_screening_vs_model_evidence.png`
- Guardrail: this is the main inferential association layer, but sparse UTI
  counts and lineage/event structure still limit claims.

### `17_lineage_analysis.R`

- Role: ST risk and lineage context.
- Reads: status map, active MLST, and VF-ready data if available.
- Main work: counts STs by primary status, computes UTI proportions, performs
  ST-specific Fisher tests, adjusts FDR, and plots lineage risk.
- Tables/reports:
  - `results/lineage/st_risk_profile.csv`
- Statistics: Fisher exact tests with BH FDR where denominators allow.
- Images:
  - `results/lineage/st_risk_plot.png`
- Guardrail: lineage summaries are sparse and diagnostic, not causal.

## Phase 3: Longitudinal and Mechanism

### `15_longitudinal_patterns.R`

- Role: participant timelines and phenotype-switch candidates.
- Reads: pairwise strain metrics, clinical status, VF/MLST context.
- Main work: orders episodes, labels transitions, groups same-strain evidence,
  identifies status-switch candidates, and summarizes Not_UTI to UTI transitions.
- Tables/reports:
  - `results/longitudinal/participant_timelines.csv`
  - `results/longitudinal/transitions.csv`
  - `results/longitudinal/phenotype_switch_candidates.csv`
  - `results/longitudinal/strain_persistence_stats.csv`
  - `results/longitudinal/not_uti_uti_transition_summary.csv`
- Statistics: descriptive transition and persistence counts.
- Images:
  - `results/longitudinal/swimmer_plot.png`

### `16_within_host_evolution.R`

- Role: selected switch-candidate genome comparisons.
- Reads: `phenotype_switch_candidates.csv`, selected FASTA assemblies, and
  strain comparison helpers.
- Main work: finds paired assemblies, runs NUCmer/MUMmer comparisons, counts
  SNPs and gene changes where possible, and writes evolution events.
- Tables/reports:
  - `results/longitudinal/evolution_events.csv`
  - `results/longitudinal/evolution_summary.txt`
  - Raw NUCmer `.snps` files for compared pairs.
- Statistics: none.
- Images: none.

### `18_annotate_variants.R`

- Role: raw SNP/indel parsing.
- Reads: `.snps` files from `16_within_host_evolution.R`.
- Main work: parses positions, reference/alternate bases, pair labels, and
  variant types into a tidy table.
- Tables/reports:
  - `results/longitudinal/annotated_snps.csv`
- Statistics: none.
- Images: none.

### `20_variant_annotation_deep.R`

- Role: gene-level SNP annotation.
- Reads: `annotated_snps.csv`, Prokka/Panaroo GFF files, and assembly metadata.
- Main work: finds available GFFs, parses gene coordinates, maps variants to
  genes/products, and reports missing coverage.
- Tables/reports:
  - `results/longitudinal/variant_annotation_detailed.csv`
- Statistics: none.
- Images: none.

### `19_host_context.R`

- Role: host/clinical context for switch candidates.
- Reads: `phenotype_switch_candidates.csv` and `clinical_merged.rds`.
- Main work: joins catheter, symptoms, antibiotics, and other clinical metadata
  onto transition pairs.
- Tables/reports:
  - `results/longitudinal/host_context_table.csv`
- Statistics: descriptive table only.
- Images: none.

### `21_publication_figures.R`

- Role: early publication-style longitudinal figures.
- Reads: timelines, switch candidates, variant annotations, and host context.
- Main work: renders a swimmer plot and mutation map.
- Tables/reports: none.
- Statistics: none.
- Images:
  - `plots/publication/Fig1_Swimmer_Plot.png`
  - `plots/publication/Fig2_Mutation_Map.png`
- Guardrail: final current figure pack is produced later by `35`.

## Phase 4: VF Deep Analysis and Final Reporting

### `23_vf_cross_sectional.R`

- Role: cross-sectional VF comparison by primary status.
- Reads: `vf_analysis_ready.csv` and `gene_map.csv`.
- Main work: computes VF burden, gene prevalence, category burden,
  depth-stratified summaries, participant summaries, paired-resident views, and
  exploratory feature/category tests.
- Tables/reports:
  - `results/vf/vf_burden_by_status.csv`
  - `results/vf/vf_gene_prevalence_by_status.csv`
  - `results/vf/vf_fisher_exploratory.csv`
  - `results/vf/vf_gene_prevalence_tests.csv`
  - `results/vf/vf_burden_by_status_stratified.csv`
  - `results/vf/vf_gene_prevalence_stratified.csv`
  - `results/vf/vf_enrichment_stratified.csv`
  - `results/vf/vf_category_burden_by_status.csv`
  - `results/vf/vf_category_enrichment.csv`
  - `results/vf/vf_category_burden_plot_summary.csv`
  - `results/vf/vf_cross_sectional_summary.txt`
- Statistics: exploratory Fisher tests with BH q-values; category
  Fisher/Wilcoxon summaries in companion CSVs.
- Images:
  - `plots/vf/vf_burden_by_status.png`
  - `plots/vf/vf_burden_boxplot.png`
  - `plots/vf/vf_burden_participant_summary.png`
  - `plots/vf/vf_burden_paired_uti_not_uti.png`
  - `plots/vf/vf_top_gene_prevalence.png`
  - `plots/vf/vf_gene_prevalence_difference_uti_not_uti.png`
  - `plots/vf/vf_gene_prevalence_heatmap.png`
  - `plots/vf/vf_category_burden_by_status.png`
  - `plots/vf/vf_category_barplot.png`

### `24_vf_longitudinal_dynamics.R`

- Role: VF stability and change across consecutive within-resident isolates.
- Reads: `vf_analysis_ready.csv`, `gene_map.csv`, and pairwise strain metrics.
- Main work: orders transitions, computes Jaccard similarity, counts gene
  gains/losses, attaches SNP/ST strain context, stratifies transition summaries,
  and writes longitudinal interpretation text.
- Tables/reports:
  - `results/vf/vf_longitudinal_transitions.csv`
  - `results/vf/vf_transition_summary_by_type.csv`
  - `results/vf/vf_transitions_stratified.csv`
  - `results/vf/vf_transition_summary_stratified.csv`
  - `results/vf/vf_same_strain_vf_stability_summary.csv`
  - `results/vf/vf_replacement_vf_change_summary.csv`
  - `results/vf/vf_strain_context_by_transition_summary.csv`
  - `results/vf/vf_same_strain_by_ST_summary.csv`
  - `results/vf/vf_longitudinal_summary.txt`
- Statistics: descriptive similarity and gain/loss summaries.
- Images:
  - `plots/vf/vf_same_strain_jaccard_by_transition.png`
  - `plots/vf/vf_same_strain_gene_gain_loss.png`
  - `plots/vf/vf_same_strain_by_ST.png`
  - `plots/vf/vf_jaccard_by_strain_context.png`
  - `plots/vf/vf_replacement_vs_same_strain_vf_change.png`
  - `plots/vf/vf_jaccard_by_transition.png`
  - `plots/vf/vf_within_host_jaccard_distribution.png`
  - `plots/vf/vf_jaccard_by_days_between_samples.png`
  - `plots/vf/vf_jaccard_same_vs_different_st.png`
  - `plots/vf/vf_gene_gain_loss_consecutive_pairs.png`

### `25_vf_lineage_vf_interaction.R`

- Role: lineage confounding and denominator diagnostics.
- Reads: `vf_analysis_ready.csv`, gene map, status/QC context, MLST/ST fields,
  and bridge diagnostics.
- Main work: summarizes VF burden by ST and by ST/status, tests lineage/status
  composition where possible, audits batch/event/timepoint structure, visualizes
  QC selection and Uricult bridge denominators, and writes a confounding summary.
- Tables/reports:
  - `results/vf/vf_burden_by_st.csv`
  - `results/vf/vf_burden_by_st_and_status.csv`
  - `results/vf/vf_lineage_confounding_summary.txt`
- Statistics: Kruskal-Wallis, Wilcoxon, and simulated Fisher diagnostics where
  possible. Sparse within-ST UTI counts limit inference.
- Images:
  - `plots/vf/vf_burden_by_st.png`
  - `plots/vf/vf_burden_by_top_st.png`
  - `plots/vf/vf_burden_st_x_status.png`
  - `plots/vf/vf_st_composition_by_status.png`
  - `plots/vf/vf_batch_by_status.png`
  - `plots/vf/vf_event_type_by_status.png`
  - `plots/vf/vf_status_timepoint_event_tile.png`
  - `plots/vf/vf_qc_selection_by_status.png`
  - `plots/vf/vf_denominator_flow.png`
  - `plots/vf/vf_uricult_join_diagnostic.png`

### `26_vf_define_gene_modules.R`

- Role: VF gene-to-module curation and module matrix.
- Reads: `vf_analysis_ready.csv`, `vf_pa_all.csv`, and `gene_map.csv`.
- Main work: assigns genes to modules with confidence levels, writes unassigned
  review rows, builds module presence/counts per episode, documents category
  definitions, and audits annotation gaps.
- Tables/reports:
  - `results/vf/vf_module_qc_report.txt`
  - `results/vf/vf_gene_annotation_gap_report.csv`
  - `results/vf/vf_gene_annotation_gap_report.txt`
  - `results/vf/gene_module_map.csv`
  - `results/vf/vf_module_presence_by_episode.csv`
  - `results/vf/vf_module_summary.csv`
  - `results/vf/vf_module_assignment_audit.csv`
  - `results/vf/gene_module_unassigned_review.csv`
  - `results/vf/vf_category_definitions.csv`
  - `results/vf/vf_module_definition_notes.md`
- Statistics: descriptive module counts and prevalence.
- Images:
  - `plots/vf/module_gene_counts.png`
  - `plots/vf/module_prevalence_by_status.png`
  - `plots/vf/vf_category_composition_by_status.png`
  - `plots/vf/vf_module_assignment_confidence.png`

### `27_vf_score_framework.R`

- Role: supplementary VF marker/system endpoint and ordination framework.
- Reads: `vf_analysis_ready.csv`, gene module map, and module presence tables.
- Main work: calculates ExPEC-like marker groups, the ExPEC-like `>=2 of 5`
  classifier, UPEC system counts/fractions, descriptive VF burden columns, and
  ordination views; summarizes by primary status and ST; runs exploratory
  endpoint tests; computes endpoint correlations; and writes PCA/PCoA
  coordinates when possible.
- Tables/reports:
  - `results/vf/vf_score_table.csv`
  - `results/vf/vf_expec_marker_definitions.csv`
  - `results/vf/vf_expec_marker_summary_by_status.csv`
  - `results/vf/vf_expec_marker_tests.csv`
  - `results/vf/vf_score_endpoint_catalog.csv`
  - `results/vf/vf_module_score_table.csv`
  - `results/vf/vf_upec_score_components.csv`
  - `results/vf/vf_score_summary_by_status.csv`
  - `results/vf/vf_score_summary_by_ST.csv`
  - `results/vf/vf_score_summary_by_status_within_ST.csv`
  - `results/vf/vf_score_tests_exploratory.csv`
  - `results/vf/vf_score_correlations.csv`
  - `results/vf/vf_pca_coordinates.csv`
  - `results/vf/vf_pca_loadings.csv`
  - `results/vf/vf_pcoa_jaccard_coordinates.csv`
  - `results/vf/vf_score_framework_summary.txt`
- Statistics: exploratory Wilcoxon endpoint tests with BH q-values, Fisher
  marker-group tests with BH q-values, Spearman endpoint correlations, PCA, and
  Jaccard PCoA.
- Images:
  - `plots/vf/vf_expec_marker_prevalence_by_status.png`
  - `plots/vf/vf_scores_by_status.png`
  - `plots/vf/vf_score_effect_summary_uti_not_uti.png`
  - `plots/vf/vf_scores_by_ST.png`
  - `plots/vf/vf_score_correlation_heatmap.png`
  - `plots/vf/vf_pca_status.png`
  - `plots/vf/vf_pca_ST.png`
  - `plots/vf/vf_pcoa_jaccard_status.png`
  - `plots/vf/vf_pcoa_jaccard_ST.png`
- Guardrail: because no validated UTI-specific VF score exists for this cohort,
  composite measures are supplementary; primary interpretation should use
  marker-group and module-level patterns with lineage and longitudinal context.

### `28_vf_transition_case_studies.R`

- Role: clinical-first transition case studies.
- Reads: ordered clinical status map, VF-ready data, module outputs, endpoint
  outputs, longitudinal outputs, and strain metrics.
- Main work: builds a transition index, keeps transitions with missing WGS/VF
  endpoints visible, computes gene/module/endpoint changes, classifies strain
  context, writes case notes, and plots transition evidence.
- Tables/reports:
  - `results/qc/vf_transition_duplicate_episode_lookup.csv`
  - `results/vf/vf_transition_case_index.csv`
  - `results/vf/vf_transition_case_summary.csv`
  - `results/vf/vf_transition_gene_changes.csv`
  - `results/vf/vf_transition_module_changes.csv`
  - `results/vf/vf_transition_score_changes.csv`
  - `results/vf/vf_transition_strain_context.csv`
  - `results/vf/vf_transition_case_notes.csv`
  - `results/vf/vf_transition_case_study_summary.txt`
- Statistics: descriptive transition deltas only.
- Images:
  - `plots/vf/vf_transition_case_timeline.png`
  - `plots/vf/vf_transition_score_slopeplot.png`
  - `plots/vf/vf_transition_module_change_heatmap.png`
  - `plots/vf/vf_transition_gene_gain_loss_tile.png`
  - `plots/vf/vf_not_uti_uti_transition_strain_context.png`
  - `plots/vf/vf_not_uti_uti_transition_case_classes.png`
  - `plots/vf/vf_transition_snp_vs_vf_jaccard.png`

### `29_vf_amr_combined_profile.R`

- Role: authoritative genomic-AMR analysis and VF/plasmid integration.
- Reads: the exact 532-assembly Longcycler manifest, matched Prokka
  FAA/FFN/FNA/GFF annotations, the 371-pair canonical transition table,
  VF-ready endpoints and plasmid/replicon tables.
- Main work: validates sequence-equivalent annotations; runs SHA-bound
  ABRicate-ResFinder at 80/80, AMRFinderPlus 4.2.7 with organism
  `Escherichia`, and ResFinder 4.7.2 with PointFinder; harmonizes determinants;
  audits caller discrepancies without voting; builds episode/resident
  prevalence and adjacent-pair profiles; performs the resident-bootstrap
  longitudinal comparison; and joins AMR to VF/plasmid context.
- Tables/reports:
  - `results/amr/provenance/input_manifest.csv`
  - `results/amr/provenance/run_manifest.csv`
  - `results/amr/provenance/tool_database_versions.csv`
  - `results/amr/provenance/published_output_manifest.csv`
  - `results/amr/harmonized_determinants_long.csv`
  - `results/amr/episode_amr_profiles.csv`
  - `results/amr/resident_amr_profiles.csv`
  - `results/amr/caller_concordance_discrepancies.csv`
  - `results/amr/caller_coverage_summary.csv`
  - `results/amr/caller_determinant_count_summary.csv`
  - `results/amr/gene_prevalence_episode_resident.csv`
  - `results/amr/class_prevalence_episode_resident.csv`
  - `results/amr/mutation_prevalence_episode_resident.csv`
  - `results/amr/resfinder_predicted_phenotypes_genomic_not_ast.csv`
  - `results/amr/adjacent_pair_amr_profiles_371.csv`
  - `results/amr/not_uti_to_uti_amr_profiles_9.csv`
  - `results/amr/longitudinal_resident_bootstrap_inference.csv`
  - `results/amr/longitudinal_summary_by_snp_context.csv`
  - `results/amr/validation_checks.csv`
  - `results/amr/interpretation_report.md`
  - `results/amr/RUN_COMPLETE.txt`
  - `results/vf_amr/vf_amr_input_availability_report.txt`
  - `results/vf_amr/vf_plasmid_combined_profile.csv`
  - `results/vf_amr/vf_amr_combined_profile_table.csv`
  - `results/vf_amr/replicon_summary_by_status.csv`
  - `results/vf_amr/vf_amr_score_summary_by_status.csv`
  - `results/vf_amr/replicon_summary_by_ST.csv`
  - `results/vf_amr/vf_amr_score_summary_by_ST.csv`
  - `results/vf_amr/vf_amr_profile_groups.csv`
  - `results/vf_amr/vf_plasmid_correlation.csv`
- Statistics: primary informative acquired-gene gain/loss across 371 adjacent
  pairs, plus acquired-gene Jaccard similarity, class change and mutation
  change. The ≤25 versus >25 SNP comparison uses 10,000 resident-cluster
  bootstrap replicates (seed 20260712) and adjusts for days between samples.
  The nine Not_UTI-to-UTI transitions are descriptive only.
- Images:
  - `plots/amr/most_prevalent_informative_acquired_genes.png`
  - `plots/amr/amr_profile_stability_by_direct_snp_context.png`
  - `plots/amr/caller_concordance_by_determinant_class.png`
  - `plots/vf_amr/replicon_burden_by_status.png`
  - `plots/vf_amr/vf_vs_replicon_scatter.png`
  - `plots/vf_amr/replicon_heatmap_top_STs.png`
  - `plots/vf_amr/vf_plasmid_analysis_scope.png`
- Guardrail: AMRFinderPlus defines the primary profile; ResFinder/PointFinder is
  complementary and ABRicate is the legacy comparison. `mdf(A)` remains raw
  and in sensitivity/QC fields but is excluded from primary burden, gain/loss
  and Jaccard calculations. Outputs are genomic predictions, not phenotypic
  AST, and AMR remains supplementary to RQ01-RQ10.

### `32_uti_not_uti_diagnostic_stats.R`

- Role: UTI vs Not_UTI denominator and sparse-count diagnostic layer.
- Reads: status map, VF-ready data, score/model outputs, exclusion reports,
  quarantine reports, and transition outputs.
- Main work: validates clinical and VF counts, writes denominator flow and
  attrition tables, audits near-miss rows, bootstraps supplementary endpoint effects, performs
  exploratory feature Fisher tests, runs leave-one-UTI sensitivity, calculates
  paired resident and transition score deltas, checks duplicate cultures, and
  writes figure metadata.
- Tables/reports:
  - `results/audit/uti_not_uti_denominator_flow.csv`
  - `results/audit/uti_not_uti_near_miss_rows.csv`
  - `results/vf/uti_not_uti_attrition_summary.csv`
  - `results/vf/uti_not_uti_bootstrap_effects.csv`
  - `results/vf/uti_not_uti_feature_fisher_exploratory.csv`
  - `results/vf/uti_not_uti_leave_one_uti_out.csv`
  - `results/vf/uti_not_uti_paired_participant_deltas.csv`
  - `results/vf/uti_not_uti_paired_participant_tests.csv`
  - `results/vf/uti_not_uti_transition_score_tests.csv`
  - `results/audit/duplicate_culture_qc_31036.csv`
  - `results/vf/uti_not_uti_power_precision_context.csv`
  - `results/vf/uti_not_uti_test_interpretation_table.csv`
  - `results/vf/uti_not_uti_diagnostic_summary.csv`
  - `results/vf/uti_not_uti_diagnostic_figure_metadata.csv`
- Statistics: participant bootstrap confidence intervals, exploratory Fisher
  feature tests, leave-one-UTI stability, paired resident sign/signed-rank
  summaries, transition sign/signed-rank summaries, and sparse expected-count
  context.
- Images:
  - `plots/clinical/uti_not_uti_clinical_rule_flow.png`
  - `plots/vf/uti_not_uti_denominator_waterfall.png`
  - `plots/clinical/uti_not_uti_near_miss_evidence_heatmap.png`
  - `plots/vf/uti_not_uti_bootstrap_effect_forest.png`
  - `plots/vf/uti_not_uti_leave_one_uti_out_stability.png`
  - `plots/vf/uti_not_uti_paired_participant_slopeplot.png`
  - `plots/vf/uti_not_uti_transition_delta_forest.png`
  - `plots/vf/uti_not_uti_sparse_power_precision.png`
  - `plots/vf/duplicate_culture_qc_31036.png`

### `33_mechanism_first_addon.R`

- Role: mechanism-first synthesis for Not_UTI to UTI transitions.
- Reads: transition case outputs, VF-ready data, status maps, module/score
  changes, strain context, host context, variant annotations, Panaroo accessory
  data, plasmid context, and script 29's validated episode/focused AMR profiles.
- Main work: builds an evidence casebook, attaches host/strain/VF/module/accessory
  evidence, consumes the nine validated focused AMR profiles, summarizes transition
  mechanisms, and writes validation checks.
- Tables/reports:
  - `results/mechanism/accessory_gene_transition_changes.csv`
  - `results/mechanism/not_uti_to_uti_casebook.csv`
  - `results/mechanism/transition_mechanism_summary.csv`
  - `results/mechanism/host_context_transition_summary.csv`
  - `results/mechanism/near_miss_sensitivity_summary.csv`
  - `results/mechanism/near_miss_sensitivity_summary.md`
  - `results/mechanism/not_uti_to_uti_casebook.md`
  - `results/mechanism/mechanism_validation_checks.csv`
- Statistics: descriptive mechanism classification and validation only.
- Images:
  - `plots/mechanism/not_uti_to_uti_case_matrix.png`
  - `plots/mechanism/strain_replacement_vs_stability.png`
  - `plots/mechanism/host_context_transition_heatmap.png`
- Guardrail: evidence buckets organize uncertainty; they do not prove causality.

### `34_robustness_first_addon.R`

- Role: robustness and claim-safety diagnostics.
- Reads: status map, VF-ready data, score tables, model outputs, QC attrition
  files, module tables, and transition outputs.
- Main work: summarizes denominator robustness, QC attrition, expanded near-miss
  sensitivity, model stability flags, leave-one-UTI sensitivity, bootstrap score
  robustness, sparse power/precision, and claim-safety matrix.
- Tables/reports:
  - `results/robustness/denominator_robustness_summary.csv`
  - `results/robustness/qc_attrition_robustness.csv`
  - `results/robustness/near_miss_expanded_score_summary.csv`
  - `results/robustness/near_miss_expanded_score_contrasts.csv`
  - `results/robustness/near_miss_expanded_module_fisher.csv`
  - `results/robustness/model_stability_summary.csv`
  - `results/robustness/leave_one_uti_sensitivity_summary.csv`
  - `results/robustness/bootstrap_score_robustness.csv`
  - `results/robustness/power_precision_summary.csv`
  - `results/robustness/near_miss_expanded_top_score_contrasts.csv`
  - `results/robustness/top_glmm_robustness_flags.csv`
  - `results/robustness/robustness_claim_matrix.csv`
  - `results/robustness/robustness_validation_checks.csv`
  - `results/robustness/robustness_summary.md`
- Statistics: near-miss score contrasts, module Fisher sensitivity, model
  stability flags, leave-one-UTI sensitivity, bootstrap robustness, and
  sparse-count precision context.
- Images:
  - `plots/robustness/qc_retention_by_status.png`
  - `plots/robustness/near_miss_score_shift.png`
  - `plots/robustness/model_stability_flags.png`

### `36_statistical_sensitivity_addon.R`

- Role: targeted statistical sensitivity.
- Reads: status map, VF-ready data, VF endpoint table, module outputs, transition
  outputs, model/lineage context, and PCoA coordinates.
- Main work: collapses supplementary endpoint tests by participant, calculates paired binary
  feature deltas, tests transition module gain/loss enrichment, fits endpoint-level
  GLMM/GLM sensitivity models, renders sensitivity figures, and writes validation.
- Tables/reports:
  - `results/statistical_sensitivity/participant_collapsed_score_values.csv`
  - `results/statistical_sensitivity/participant_collapsed_score_summary.csv`
  - `results/statistical_sensitivity/participant_collapsed_score_tests.csv`
  - `results/statistical_sensitivity/paired_binary_feature_catalog.csv`
  - `results/statistical_sensitivity/paired_binary_feature_deltas.csv`
  - `results/statistical_sensitivity/paired_binary_feature_sensitivity.csv`
  - `results/statistical_sensitivity/transition_module_gain_loss_enrichment.csv`
  - `results/statistical_sensitivity/not_uti_to_uti_module_change_matrix.csv`
  - `results/statistical_sensitivity/score_glmm_sensitivity.csv`
  - `results/statistical_sensitivity/statistical_sensitivity_figure_metadata.csv`
  - `results/statistical_sensitivity/statistical_sensitivity_summary.md`
  - `results/statistical_sensitivity/statistical_sensitivity_validation_checks.csv`
- Statistics: Wilcoxon, paired sign/signed-rank tests, bootstrap confidence
  intervals, Fisher odds ratios, GLMM/GLM odds ratios, and BH-adjusted q-values.
- Images:
  - `plots/statistical_sensitivity/lineage_confounding_panel.png`
  - `plots/statistical_sensitivity/lineage_confounding_panel.pdf`
  - `plots/statistical_sensitivity/paired_resident_expec_marker_slopeplot.png`
  - `plots/statistical_sensitivity/paired_resident_expec_marker_slopeplot.pdf`
  - `plots/statistical_sensitivity/not_uti_to_uti_module_gain_loss_heatmap.png`
  - `plots/statistical_sensitivity/not_uti_to_uti_module_gain_loss_heatmap.pdf`
  - `plots/statistical_sensitivity/vf_module_pcoa_primary_status.png`
  - `plots/statistical_sensitivity/vf_module_pcoa_primary_status.pdf`
- Guardrail: deliberately targeted, not a broad discovery layer.

### `30_vf_project_summary_tables.R`

- Role: final summary table and figure-index builder.
- Reads: clinical, WGS, MLST, VF, model, module, score, longitudinal,
  transition, lineage, diagnostic, mechanism, robustness, and statistical
  sensitivity outputs.
- Main work: validates denominators and freshness, writes numbered summary
  tables, creates a visualization audit and figure index, writes a key-results
  Markdown summary, and bundles summary tables.
- Tables/reports:
  - `results/summary/table_01_cohort_episode_flow.csv`
  - `results/summary/table_02_clinical_status_counts.csv`
  - `results/summary/table_03_wgs_qc_summary.csv`
  - `results/summary/table_04_mlst_st_summary.csv`
  - `results/summary/table_05_vf_gene_category_summary.csv`
  - `results/summary/table_06_vf_module_summary.csv`
  - `results/summary/table_07_vf_score_summary_by_status.csv`
  - `results/summary/table_08_vf_score_summary_by_ST.csv`
  - `results/summary/table_09_longitudinal_vf_stability_summary.csv`
  - `results/summary/table_09_longitudinal_vf_stability.csv`
  - `results/summary/table_10_not_uti_uti_transition_cases.csv`
  - `results/summary/table_10_all_transition_context.csv`
  - `results/summary/table_11_lineage_context_summary.csv`
  - `results/summary/table_12_missing_data_audit.csv`
  - `results/summary/table_13_genomic_amr_summary.csv`
  - `results/summary/table_14_uti_not_uti_diagnostics.csv`
  - `results/summary/table_15_uti_not_uti_test_interpretation.csv`
  - `results/summary/table_16_uti_not_uti_diagnostic_figures.csv`
  - `results/summary/table_17_participant_collapsed_score_tests.csv`
  - `results/summary/table_18_paired_binary_feature_sensitivity.csv`
  - `results/summary/table_19_transition_module_gain_loss_enrichment.csv`
  - `results/summary/table_20_score_glmm_sensitivity.csv`
  - `results/summary/table_21_statistical_sensitivity_validation.csv`
  - `results/vf/vf_figure_index.csv`
  - `results/vf/vf_visualisation_audit.csv`
  - `results/summary/summary_qc_log.txt`
  - `results/summary/final_key_results_summary.md`
  - `results/summary/final_summary_tables.rds`
  - `results/summary/final_summary_tables.xlsx` when `writexl` is available.
- Statistics: no new tests. It summarizes and labels upstream evidence.
- Images: no new images.

### `35_final_figure_pack.R`

- Role: final manuscript-style figure pack.
- Reads: validation checks, status map, VF-ready data, mechanism casebook,
  robustness outputs, diagnostic outputs, statistical sensitivity outputs,
  variant tables, and summary tables.
- Main work: validates required inputs, renders main and supplementary figures,
  writes a manifest, and writes captions.
- Tables/reports:
  - `results/final_figures/final_figure_validation_checks.csv`
  - `results/final_figures/final_figure_manifest.csv`
  - `results/final_figures/final_figure_captions.md`
- Statistics: no new tests. Figures combine validated upstream diagnostics.
- Main figures:
  - `plots/final_figures/primary_denominator_and_uncertainty.png` and `.pdf`
  - `plots/final_figures/not_uti_to_uti_mechanism_casebook.png` and `.pdf`
  - `plots/final_figures/strain_stability_and_host_context.png` and `.pdf`
  - `plots/final_figures/global_vf_signal_and_robustness.png` and `.pdf`
- Supplementary figures:
  - `plots/final_figures/transition_mechanisms_by_transition_type.png` and `.pdf`
  - `plots/final_figures/accessory_plasmid_amr_changes.png` and `.pdf`
  - `plots/final_figures/near_miss_and_sparse_precision.png` and `.pdf`
  - `plots/final_figures/leave_one_uti_stability.png` and `.pdf`
  - `plots/final_figures/prioritised_variant_map.png` and `.pdf`
  - `plots/final_figures/lineage_confounding_diagnostic.png` and `.pdf`
  - `plots/final_figures/paired_resident_expec_marker.png` and `.pdf`
  - `plots/final_figures/not_uti_to_uti_module_gain_loss.png` and `.pdf`
  - `plots/final_figures/vf_module_pcoa_primary_status.png` and `.pdf`

### `scripts/archive_legacy_asb_uti_outputs.R`

- Role: archive stale generated ASB-vs-UTI outputs after current UTI-vs-Not_UTI
  outputs are produced.
- Reads: `results/` and `plots/`.
- Main work: finds generated filenames that still advertise old ASB-vs-UTI
  contrasts, moves them to legacy archive locations, and writes a manifest.
- Tables/reports:
  - `results/qc/legacy_status_archive_manifest.csv`
- Statistics: none.
- Images: none generated.
- Guardrail: archival reduces confusion; archived outputs are not current
  evidence.

### `scripts/verify_uti_not_uti_alignment.R`

- Role: final post-run alignment gate.
- Reads: status maps, VF-ready dataset, summary outputs, quarantine reports,
  transition labels, and current results/plots directories.
- Main work: checks required primary status fields, expected denominators,
  duplicates/manual curation, poster-map freshness, VF-ready alignment,
  quarantine rows, transition labels, and absence of current generated
  ASB-vs-UTI outputs.
- Tables/reports:
  - `results/qc/uti_not_uti_alignment_checks.csv`
  - `results/qc/uti_not_uti_alignment_checks.txt`
- Statistics: none.
- Images: none.
- Guardrail: validation gate only. It does not transform data.

## Statistical Methods Catalog

| Method or output | Scripts | Interpretation |
|:--|:--|:--|
| Rule/count summaries | `00b`, `00c`, `12a`, `22`, `30` | Denominator and QC explanation. |
| Fisher exact or Fisher-style tests | `12a`, `13`, `08`, `14`, `17`, `23`, `32`, `34`, `36` | Exact-style tests are useful with sparse counts, but repeated residents and lineage can still confound results. |
| Wilcoxon tests | `03`, `23`, `25`, `27`, `32`, `36` | Mostly exploratory or sensitivity summaries unless explicitly collapsed/paired. |
| Kruskal-Wallis summaries | `11`, `25` | Descriptive group comparisons, not causal inference. |
| GLMM with participant random intercept | `04`, `14`, `36` | Participant-aware model layer; still limited by sparse UTI counts and separation risks. |
| GLM fallback | `04`, `14`, `36` | Used when mixed models fail or are requested as simpler sensitivity. |
| BH FDR adjustment | `04`, `14`, `17`, `23`, `27`, `32`, `36` | Controls false-discovery rate within the tested family, not across all exploratory analyses in the project. |
| Spearman correlation | `27`, `29` | Describes score/plasmid or score/score relationships; not prediction. |
| PCA and Jaccard PCoA | `27`, `36`, `35` | Descriptive ordination; not adjusted for repeated residents or lineage. |
| Participant bootstrap | `32`, `34`, `35` | Descriptive uncertainty under sparse UTI denominator. |
| Leave-one-UTI sensitivity | `32`, `34`, `35` | Tests dependence on individual UTI rows; instability is expected with sparse UTI counts. |
| Paired sign/signed-rank summaries | `32`, `36` | Resident-paired sensitivity, only among residents observed in both statuses. |

## Figure Manifest Crosswalk

The most authoritative figure-level metadata are generated near the end of the
pipeline:

- `results/vf/vf_figure_index.csv`: VF analysis figure registry with script,
  denominator, statistical test, evidence label, output path, and limitation.
- `results/vf/vf_visualisation_audit.csv`: visualization audit and whether
  expected figure files exist.
- `results/vf/uti_not_uti_diagnostic_figure_metadata.csv`: denominator and
  sparse-count diagnostic figure metadata from script `32`.
- `results/statistical_sensitivity/statistical_sensitivity_figure_metadata.csv`:
  statistical sensitivity figure metadata from script `36`.
- `results/final_figures/final_figure_manifest.csv`: final main and
  supplementary figure pack manifest from script `35`.

When explaining a figure to a collaborator, prefer the manifest fields over
memory: figure ID, source script, source inputs, denominator, statistical test,
evidence type, and interpretation limitations.

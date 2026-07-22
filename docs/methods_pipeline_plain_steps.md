# Plain-English Bullet Methods From the Full Pipeline

This document explains what happens in the analysis pipeline in simple bullet points. It is meant to help you understand the methods step by step. It is not a polished manuscript methods section.

## Sources Used

- Main run order: `RUN_COMPLETE_ANALYSIS.sh`
- Current pipeline explanations: `FOLDER_MAP.md`, `docs/LECTURER_README.md`, `docs/pipeline_architecture.md`, `docs/workflow_flowchart.md`, and `docs/workflow_output_catalog.md`
- Superseded drafts and older output guidance are retained under `docs/legacy_asb_uti_docs/superseded_2025_2026/` for historical reference only.
- Shared setup and helper files: `00_config.R`, `R/clinical_helpers.R`, `R/wgs_helpers.R`, `R/pipeline_qc_helpers.R`, `R/plot_helpers.R`, `11_compare_strains_helpers.R`
- WGS prerequisite files: `00_make_assembly_metadata.r`, `assembly_metadata.csv`
- Clinical runner scripts: `00a_load_clean_clinical.R`, `00b_classify_episodes.R`, `00d_derive_plot_timepoints.R`, `00c_plot_clinical_summary.R`
- WGS runner scripts: `12a_wgs_qc.R`, `12b_core_snp.R`, `12c_panaroo.R`, `13_visualise_panaroo_selection.R`, `02_gene_presence_analysis.R`, `06_MLST.R`
- Additional descriptive scripts: optional `03_plotting.R`, `04_gene_breakdown.R`, `05_gene_overview_plots.R`, `07_explore_MLST.R`, `08_core_vs_plasmid.R`, `09_inc_plasmid_network.R`, `10_replicon_heatmap.R`
- Comparative and longitudinal scripts: `11_compare_strains.R`, `22_vf_build_analysis_dataset.R`, `14_genotype_phenotype_model.R`, `17_lineage_analysis.R`, `15_longitudinal_patterns.R`, `16_within_host_evolution.R`, `18_annotate_variants.R`, `20_variant_annotation_deep.R`, `19_host_context.R`, `21_publication_figures.R`
- VF deep-analysis scripts: `23_vf_cross_sectional.R`, `24_vf_longitudinal_dynamics.R`, `25_vf_lineage_vf_interaction.R`, `26_vf_define_gene_modules.R`, `27_vf_score_framework.R`, `28_vf_transition_case_studies.R`, `29_vf_amr_combined_profile.R`, `32_uti_not_uti_diagnostic_stats.R`, `33_mechanism_first_addon.R`, `34_robustness_first_addon.R`, `36_statistical_sensitivity_addon.R`, `30_vf_project_summary_tables.R`, `35_final_figure_pack.R`
- Final hygiene and validation scripts: `scripts/archive_legacy_asb_uti_outputs.R`, `scripts/verify_uti_not_uti_alignment.R`

## Big Picture

- The clinical data are cleaned first.
- Each episode is classified using the current primary outcome: `UTI` or `Not_UTI`.
- The bacterial genome data are quality-checked, compared, and converted into gene, virulence factor, plasmid, and sequence type tables.
- Isolates from the same participant are compared over time.
- The pipeline asks whether a participant kept the same strain, acquired a new strain, or had a clinical switch while carrying a similar strain.
- Virulence factor patterns are then checked cross-sectionally, longitudinally, by lineage, by gene module, and with robustness/sensitivity checks.
- Final tables, final figures, and validation checks are produced at the end.

## Phase 0: Clinical Data Foundation

### Step 1: Set up shared project rules

- Script/tool: `00_config.R` plus helper files
- What goes in:
  - No clinical episode table directly.
  - Project paths, constants, thresholds, helper functions, and shared output folders.
- What happens:
  - The pipeline defines where inputs, outputs, plots, and logs live.
  - It also defines shared interpretation rules, including the same-strain SNP threshold.
- What comes out:
  - Shared settings used by later scripts.
- Why it matters:
  - Every script needs to use the same folders, status definitions, colors, and thresholds.

### Step 2: Load and harmonize raw clinical files

- Script/tool: `00a_load_clean_clinical.R`
- What goes in:
  - Raw clinical batch files from `data/inputs/`.
- What happens:
  - Batch files are read and combined.
  - Column names and participant/timepoint labels are standardized.
  - Symptom columns are detected and mapped for later rule-based classification.
- What comes out:
  - `results/clinical/intermediate/clinical_merged.rds`
  - `results/clinical/uti_detected_column_map.csv`
- Why it matters:
  - Before classifying infection status, all clinical data need to be in one consistent table.

### Step 3: Classify episodes as UTI or Not_UTI

- Script/tool: `00b_classify_episodes.R`
- What goes in:
  - `clinical_merged.rds`
  - Clinical helper rules from `R/clinical_helpers.R`
- What happens:
  - The script applies catheter-aware symptom rules.
  - It checks culture support using the >=10^3 CFU/mL threshold.
  - It creates the current primary status: `UTI` versus `Not_UTI`.
  - It keeps older ASB/UTI/Negative labels only as legacy comparison fields.
- What comes out:
  - `results/clinical/status_map.csv`
  - Classification audit files, reclassification tables, symptom-rule audit files, and CFU-threshold audit files.
- Why it matters:
  - `status_map.csv` is the master clinical file used by nearly every downstream script.

### Step 4: Create display-only timepoints for event samples

- Script/tool: `00d_derive_plot_timepoints.R`
- What goes in:
  - `status_map.csv`
  - Raw collection-date information when available.
- What happens:
  - Uricult/event samples are placed into an approximate visual order.
  - The script records whether placement was possible and how confident it was.
- What comes out:
  - `results/clinical/status_map_with_poster_tp.csv`
  - `results/clinical/unplaced_uricult_rows.csv`
  - `results/clinical/display_only_uricult_rows.csv`
- Why it matters:
  - This helps plots show event samples in a readable order.
- Guardrail:
  - These poster/display timepoints are for visualization only, not statistical modelling.

### Step 5: Plot and audit the clinical denominator

- Script/tool: `00c_plot_clinical_summary.R`
- What goes in:
  - `status_map.csv`
  - Display timepoint map when available.
- What happens:
  - The script summarizes UTI and Not_UTI counts.
  - It shows how episodes move from legacy labels to the current primary labels.
  - It visualizes symptom-rule, culture-rule, subgroup, and assembly-availability context.
- What comes out:
  - Clinical plots in `plots/clinical/`
  - Clinical summary and inspection tables in `results/clinical/`
- Why it matters:
  - It checks what denominator the genomic analysis will use.

## Phase 1: Whole-Genome Sequencing Processing

### Step 6: Check assembly quality and choose canonical assemblies

- Script/tool: `12a_wgs_qc.R`
- What goes in:
  - `assembly_metadata.csv`
  - Assembly FASTA files
  - `status_map.csv` when available.
- What happens:
  - The script checks assembly quality, including contig count, total length, N50, GC content, missing files, and manual exclusions.
  - It selects the canonical assembly to represent each participant-timepoint.
- What comes out:
  - `results/wgs/qc_summary.csv`
  - `results/qc/canonical_assembly_selection.csv`
  - QC and selection-bias reports.
- Why it matters:
  - The rest of the genome pipeline should use the same selected assembly set.

### Step 7: Build core SNP distances

- Script/tool: `12b_core_snp.R`
- What goes in:
  - Canonical selected assemblies from `12a_wgs_qc.R`.
- What happens:
  - FASTA inputs are staged and fingerprinted.
  - Parsnp creates a core-genome alignment.
  - `snp-dists` calculates pairwise SNP distances.
  - Strain-distance outputs and freshness reports are written.
- What comes out:
  - `results/wgs/core/snp_dists.tsv`
  - `results/wgs/core/strain_pairs.csv`
  - Core alignment/tree artifacts and staleness reports.
- Why it matters:
  - Core SNP distance is the main evidence used to judge same-strain persistence.

### Step 8: Build the pangenome input and run Panaroo

- Script/tool: `12c_panaroo.R`
- What goes in:
  - Canonical QC-passing assemblies.
  - GFF annotation files, or FASTA files that can be annotated.
- What happens:
  - Missing or empty GFFs are regenerated when possible.
  - The script writes a manifest of input GFFs.
  - Panaroo is run or skipped depending on whether outputs are fresh.
- What comes out:
  - Panaroo outputs under `results/wgs/pan/`
  - GFF manifests, missing-GFF reports, and staleness reports.
- Why it matters:
  - The pangenome tells the pipeline which genes are present across isolates.

### Step 9: Visualize which assemblies entered WGS/Panaroo

- Script/tool: `13_visualise_panaroo_selection.R`
- What goes in:
  - QC summary files.
  - Canonical selection files.
  - Panaroo manifests.
  - Clinical status map.
- What happens:
  - The script shows which participant-timepoints were selected, missing, or excluded.
  - It checks whether WGS/Panaroo selection is imbalanced by clinical status when possible.
- What comes out:
  - Panaroo selection plots in `plots/wgs/`
  - Selection diagnostics in `results/wgs/` and `results/qc/`
- Why it matters:
  - This explains why some clinical rows appear or do not appear in genomic analysis.

### Step 10: Detect virulence factor genes and build a binary matrix

- Script/tool: `02_gene_presence_analysis.R`
- What goes in:
  - ABRicate/VFDB-style hit files or generated gene-detection outputs.
  - Selected assembly information.
- What happens:
  - VF hits are filtered by identity and coverage rules.
  - The script builds an isolate-by-gene matrix.
  - Each gene is marked present or absent for each isolate.
- What comes out:
  - `results/vf/vf_hits_all.csv`
  - `results/vf/vf_hits_all.rds`
  - `results/vf/vf_pa_all.csv`
  - `results/vf/stats_gene_level.csv`
- Why it matters:
  - This is the base virulence factor table used by the later VF analyses.
- Guardrail:
  - Gene presence does not prove gene expression, activity, virulence, or causality.

### Step 11: Add active sequence type information

- Script/tool: `06_MLST.R`
- What goes in:
  - Provider/RIVM MLST calls when available.
  - Local MLST provenance when needed.
  - Assembly metadata.
- What happens:
  - ST labels are normalized.
  - Provider-preferred and local calls are reconciled.
  - Source/provenance checks are written.
- What comes out:
  - `results/mlst/mlst_all.tsv`
  - `results/mlst/mlst_provider_preferred.csv`
  - MLST source audit files.
- Why it matters:
  - ST gives lineage context, such as whether isolates belong to the same broad bacterial lineage.
- Guardrail:
  - Same ST does not prove same strain. SNP evidence is still primary.

## Phase 1b: Additional Descriptive Genomics and Exploration

### Step 12: Optionally run legacy exploratory plots

- Script/tool: `03_plotting.R`
- What goes in:
  - VF matrix, status map, MLST, pairwise metrics, and optional plasmid/nitrate files.
- What happens:
  - Older exploratory plots and descriptive summaries are made only if `RUN_LEGACY_EXPLORATORY_PLOTS=1`.
- What comes out:
  - Legacy exploratory plots and descriptive statistics.
- Why it matters:
  - This can help with review, but it is not the current primary UTI versus Not_UTI analysis.

### Step 13: Create focused gene and category summaries

- Script/tool: `04_gene_breakdown.R`
- What goes in:
  - `vf_pa_all.csv`
  - Existing or generated gene-map information.
- What happens:
  - VF genes are mapped to categories.
  - Per-sample gene category counts are created.
  - Optional focused gene tests are written when possible.
- What comes out:
  - `results/vf/gene_map.csv`
  - `results/vf/annotated_gene_table.csv`
  - `results/vf/per_sample_category_counts.csv`
- Why it matters:
  - Later scripts need gene categories and gene metadata to interpret the VF matrix.

### Step 14: Make descriptive gene overview plots

- Script/tool: `05_gene_overview_plots.R`
- What goes in:
  - `vf_pa_all.csv`
- What happens:
  - Common and variable genes are summarized.
  - Prevalence bars and heatmaps are produced.
- What comes out:
  - `plots/vf/gene_prevalence_bar.png`
  - `plots/vf/variable_gene_heatmap.png`
- Why it matters:
  - This gives a quick view of the gene landscape before statistical testing.

### Step 15: Explore ST distribution

- Script/tool: `07_explore_MLST.R`
- What goes in:
  - Active MLST output.
  - Assembly metadata when available.
- What happens:
  - STs are counted and checked across samples.
  - Common STs are plotted.
- What comes out:
  - `results/mlst/ST_frequencies.csv`
  - `plots/mlst/top20_STs.png`
- Why it matters:
  - It shows which bacterial lineages dominate the cohort.

### Step 16: Explore plasmid and replicon context

- Scripts/tools: `08_core_vs_plasmid.R`, `09_inc_plasmid_network.R`, `10_replicon_heatmap.R`
- What goes in:
  - Assemblies.
  - Active MLST output.
  - PlasmidFinder/pMLST outputs or inputs.
- What happens:
  - Plasmid replicons are detected or loaded.
  - Replicon presence/absence tables are built.
  - Plasmid co-occurrence networks and heatmaps are generated.
  - ST-plasmid association checks are run when valid.
- What comes out:
  - `results/plasmids/plasmidfinder_presence_absence.csv`
  - `results/plasmids/plasmidfinder_hits_long.csv`
  - `results/mlst/ST_plasmid_associations.csv`
  - Plasmid network and heatmap plots.
- Why it matters:
  - Plasmids are accessory-genome context and can help interpret gene gain/loss.
- Guardrail:
  - Plasmid similarity does not override SNP-based strain evidence.

## Phase 2: Comparative Genomics

### Step 17: Compare isolates within participants

- Script/tool: `11_compare_strains.R --participants ALL`
- What goes in:
  - `status_map.csv`
  - Core SNP distances.
  - `vf_pa_all.csv`
  - Plasmid/replicon data.
  - MLST and assembly metadata.
- What happens:
  - Isolate pairs are built, especially within each participant.
  - The script calculates SNP distance, VF similarity, plasmid similarity, and lineage context.
  - The main same-strain rule is:
    - 0-25 core-genome SNPs = strong same-strain evidence.
    - More than 25 SNPs = above the same-strain SNP threshold.
    - Missing SNPs = missing SNP evidence.
  - ST is reported separately as `Same ST`, `Different ST`, or `Missing ST evidence`.
  - Same ST alone does not promote a pair to same-strain status when SNP evidence is missing.
- What comes out:
  - `results/strain_compare/pairwise_metrics.csv`
  - `results/strain_compare/summary_counts.csv`
  - `results/strain_compare/summary_by_participant.csv`
  - Strain comparison plots.
- Why it matters:
  - This is the main step for deciding whether a participant likely kept the same strain or acquired a replacement strain.

### Step 18: Build the canonical VF-ready dataset

- Script/tool: `22_vf_build_analysis_dataset.R`
- What goes in:
  - `vf_pa_all.csv`
  - `status_map.csv`
  - QC/selection files.
  - Active MLST.
  - Manual curation fields.
- What happens:
  - VF rows are joined to the primary clinical status.
  - The script checks unmatched rows, duplicate keys, Uricult bridging, manual exclusions, and gene annotation gaps.
  - It adds ST and VF burden columns.
- What comes out:
  - `results/vf/vf_analysis_ready.csv`
  - `results/vf/vf_binary_uti_ready.csv`
  - VF dataset diagnostics and QC reports.
- Why it matters:
  - Later VF scripts should use this file rather than rejoining raw clinical and VF tables independently.

### Step 19: Model gene/status associations

- Script/tool: `14_genotype_phenotype_model.R`
- What goes in:
  - VF-ready data.
  - Status map.
  - VF features.
  - Plasmid replicons.
  - MLST and assembly metadata.
- What happens:
  - Candidate features are screened.
  - Fisher tests are used for exploratory screening.
  - GLMMs with participant random intercepts are fitted when possible.
  - GLM fallback models are used when mixed models fail.
  - P-values are adjusted using Benjamini-Hochberg FDR.
- What comes out:
  - `results/models/gwas_univariable_stats.csv`
  - `results/models/gwas_multivariable_glmm.csv`
  - Model warnings and plots.
- Why it matters:
  - This is the main participant-aware gene/status association layer.
- Guardrail:
  - Sparse UTI counts and model separation risks mean results still need cautious interpretation.

### Step 20: Check lineage risk

- Script/tool: `17_lineage_analysis.R`
- What goes in:
  - Active MLST data.
  - Status map or VF-ready status data.
- What happens:
  - STs are counted by UTI and Not_UTI status.
  - UTI proportions are calculated for STs.
  - Fisher tests with FDR adjustment are run where denominators allow.
- What comes out:
  - `results/lineage/st_risk_profile.csv`
  - `results/lineage/st_risk_plot.png`
- Why it matters:
  - This checks whether apparent UTI patterns might be linked to bacterial lineage.

## Phase 3: Longitudinal and Mechanism Analysis

### Step 21: Reconstruct participant timelines

- Script/tool: `15_longitudinal_patterns.R`
- What goes in:
  - `results/strain_compare/pairwise_metrics.csv`
  - Clinical status data.
  - VF/MLST context.
- What happens:
  - Episodes are ordered within each participant.
  - Same-strain graph edges are used to group isolates into strain IDs.
  - Clinical transitions are labelled.
  - Candidate phenotype switches are identified, especially Not_UTI -> UTI.
- What comes out:
  - `results/longitudinal/participant_timelines.csv`
  - `results/longitudinal/transitions.csv`
  - `results/longitudinal/phenotype_switch_candidates.csv`
  - `results/longitudinal/swimmer_plot.png`
- Why it matters:
  - This turns separate samples into patient-level stories over time.

### Step 22: Compare genomes in switch candidates

- Script/tool: `16_within_host_evolution.R`
- What goes in:
  - `phenotype_switch_candidates.csv`
  - Selected FASTA assemblies.
- What happens:
  - Paired assemblies are compared for selected switch candidates.
  - NUCmer/MUMmer is used to identify SNPs and indels.
  - Gene gain/loss evidence is summarized where possible.
- What comes out:
  - `results/longitudinal/evolution_events.csv`
  - `results/longitudinal/evolution_summary.txt`
  - Raw NUCmer `.snps` files.
- Why it matters:
  - This asks whether a clinical switch is accompanied by bacterial genomic change.

### Step 23: Parse raw variants

- Script/tool: `18_annotate_variants.R`
- What goes in:
  - Raw `.snps` files from NUCmer.
- What happens:
  - Variant positions, bases, pair labels, and variant types are parsed into a tidy table.
- What comes out:
  - `results/longitudinal/annotated_snps.csv`
- Why it matters:
  - Raw variant files are hard to read directly, so this turns them into structured data.

### Step 24: Map variants to genes

- Script/tool: `20_variant_annotation_deep.R`
- What goes in:
  - `annotated_snps.csv`
  - Prokka/Panaroo GFF annotation files.
  - Assembly metadata.
- What happens:
  - Variants are mapped to genes or intergenic regions.
  - Gene names and product annotations are added where available.
  - Missing GFF coverage is reported.
- What comes out:
  - `results/longitudinal/variant_annotation_detailed.csv`
- Why it matters:
  - This helps move from "there was a mutation" to "this mutation may affect this gene or region."

### Step 25: Add host and clinical context

- Script/tool: `19_host_context.R`
- What goes in:
  - `phenotype_switch_candidates.csv`
  - `clinical_merged.rds`
- What happens:
  - Catheter status, symptoms, antibiotics, and other clinical metadata are joined onto switch pairs.
- What comes out:
  - `results/longitudinal/host_context_table.csv`
- Why it matters:
  - Clinical switches may reflect host or care context, not only bacterial genetic change.

### Step 26: Make early publication-style longitudinal figures

- Script/tool: `21_publication_figures.R`
- What goes in:
  - Timelines.
  - Switch candidates.
  - Variant annotations.
  - Host context.
- What happens:
  - Early swimmer and mutation-map figures are rendered.
- What comes out:
  - `plots/publication/Fig1_Swimmer_Plot.png`
  - `plots/publication/Fig2_Mutation_Map.png`
- Why it matters:
  - These figures help communicate longitudinal strain and mechanism findings.
- Guardrail:
  - The current final figure pack is produced later by `35_final_figure_pack.R`.

## Phase 4: VF Deep Analysis and Final Reporting

### Step 27: Compare VF profiles by UTI status

- Script/tool: `23_vf_cross_sectional.R`
- What goes in:
  - `vf_analysis_ready.csv`
  - `gene_map.csv`
- What happens:
  - VF burden, gene prevalence, and category burden are compared between UTI and Not_UTI.
  - Exploratory Fisher/Wilcoxon-style summaries are written where appropriate.
- What comes out:
  - `results/vf/vf_burden_by_status.csv`
  - `results/vf/vf_gene_prevalence_by_status.csv`
  - `results/vf/vf_fisher_exploratory.csv`
  - Cross-sectional VF plots.
- Why it matters:
  - This gives the main descriptive view of whether VF content differs by clinical status.
- Guardrail:
  - These screens are exploratory because repeated episodes from the same participant are not independent.

### Step 28: Measure VF changes over time

- Script/tool: `24_vf_longitudinal_dynamics.R`
- What goes in:
  - `vf_analysis_ready.csv`
  - `gene_map.csv`
  - `pairwise_metrics.csv`
- What happens:
  - Consecutive within-participant transitions are ordered.
  - VF Jaccard similarity is calculated.
  - Gene gains and losses are counted.
  - Results are stratified by transition type, strain context, ST, and timing.
- What comes out:
  - `results/vf/vf_longitudinal_transitions.csv`
  - `results/vf/vf_transition_summary_by_type.csv`
  - Same-strain and replacement VF stability summaries.
  - Longitudinal VF plots.
- Why it matters:
  - This asks whether VF profiles stay stable or change within participants over time.

### Step 29: Check lineage confounding

- Script/tool: `25_vf_lineage_vf_interaction.R`
- What goes in:
  - `vf_analysis_ready.csv`
  - Gene map and MLST/ST context.
  - QC and bridge diagnostics.
- What happens:
  - VF burden is summarized by ST and by ST/status.
  - Batch, event type, timepoint, QC selection, and Uricult bridging are checked.
  - Sparse within-ST UTI counts are flagged.
- What comes out:
  - `results/vf/vf_burden_by_st.csv`
  - `results/vf/vf_burden_by_st_and_status.csv`
  - `results/vf/vf_lineage_confounding_summary.txt`
  - Lineage and denominator diagnostic plots.
- Why it matters:
  - It helps separate possible VF effects from lineage and denominator structure.

### Step 30: Group genes into VF modules

- Script/tool: `26_vf_define_gene_modules.R`
- What goes in:
  - `vf_analysis_ready.csv`
  - `vf_pa_all.csv`
  - `gene_map.csv`
- What happens:
  - VF genes are assigned to interpretable modules.
  - Module presence and counts are calculated per episode.
  - Unassigned or uncertain genes are saved for review.
- What comes out:
  - `results/vf/gene_module_map.csv`
  - `results/vf/vf_module_presence_by_episode.csv`
  - `results/vf/vf_module_summary.csv`
  - Module QC reports and plots.
- Why it matters:
  - Modules make the gene matrix easier to interpret biologically.

### Step 31: Calculate supplementary VF marker/system endpoints and ordination views

- Script/tool: `27_vf_score_framework.R`
- What goes in:
  - `vf_analysis_ready.csv`
  - Gene module map.
  - Module presence table.
- What happens:
  - ExPEC-like marker groups, the ExPEC-like `>=2 of 5` classifier, UPEC system counts/fractions, descriptive VF burden columns, PCA, and Jaccard PCoA views are calculated where possible.
  - Supplementary endpoints are summarized by UTI status and ST.
  - Exploratory endpoint tests, marker-group tests, and correlations are written.
- What comes out:
  - `results/vf/vf_score_table.csv`
  - `results/vf/vf_expec_marker_definitions.csv`
  - `results/vf/vf_expec_marker_summary_by_status.csv`
  - `results/vf/vf_score_endpoint_catalog.csv`
  - `results/vf/vf_module_score_table.csv`
  - `results/vf/vf_score_tests_exploratory.csv`
  - PCA/PCoA tables and endpoint plots.
- Why it matters:
  - The framework separates literature-aligned marker groups, module-level UPEC-associated systems, and descriptive gene burden.
- Guardrail:
  - Because no validated UTI-specific VF score exists for this cohort, composite endpoints are supplementary and not validated clinical predictors.

### Step 32: Build transition case studies

- Script/tool: `28_vf_transition_case_studies.R`
- What goes in:
  - Ordered status map.
  - VF-ready data.
  - Module and supplementary endpoint outputs.
  - Longitudinal outputs.
  - Strain metrics.
- What happens:
  - Every clinical transition is indexed.
  - The script records whether WGS/VF/module/endpoint evidence is available.
  - Gene, module, and supplementary endpoint changes are calculated for transitions.
  - Strain context is attached.
- What comes out:
  - `results/vf/vf_transition_case_index.csv`
  - `results/vf/vf_transition_case_summary.csv`
  - `results/vf/vf_transition_gene_changes.csv`
  - `results/vf/vf_transition_module_changes.csv`
  - `results/vf/vf_transition_score_changes.csv`
  - Transition plots.
- Why it matters:
  - This creates a clinical-first view of what changed around important transitions.

### Step 33: Run genomic AMR and integrate VF/plasmid context

- Script/tool: `29_vf_amr_combined_profile.R`
- What goes in:
  - VF-ready data.
  - VF endpoint tables.
  - Plasmid/replicon tables.
  - The exact 532 selected Longcycler FASTAs and matched Prokka annotations.
  - The canonical 371 adjacent-pair table.
- What happens:
  - SHA-bound ABRicate-ResFinder, AMRFinderPlus 4.2.7 and
    ResFinder/PointFinder 4.7.2 calls are run or reused from valid caches.
  - AMRFinderPlus uses Prokka annotation mode, organism `Escherichia`, and the
    pinned project-local database recorded under `data/amr_runtime/databases`.
  - AMRFinderPlus acquired genes and known mutations define the primary profile.
  - Calls are harmonized and discrepancies are audited without a voting rule.
  - Episode/resident prevalence and all 371 adjacent-pair profiles are built.
  - The ≤25 versus >25 SNP comparison uses 10,000 resident bootstraps and
    time-between-samples adjustment.
  - The nine Not_UTI-to-UTI transitions are described without regression.
  - VF and plasmid context are joined to the validated episode AMR profiles.
- What comes out:
  - `results/amr/episode_amr_profiles.csv`
  - `results/amr/resident_amr_profiles.csv`
  - `results/amr/harmonized_determinants_long.csv`
  - `results/amr/caller_coverage_summary.csv`
  - `results/amr/gene_prevalence_episode_resident.csv`
  - `results/amr/class_prevalence_episode_resident.csv`
  - `results/amr/mutation_prevalence_episode_resident.csv`
  - `results/amr/resfinder_predicted_phenotypes_genomic_not_ast.csv`
  - `results/amr/adjacent_pair_amr_profiles_371.csv`
  - `results/amr/not_uti_to_uti_amr_profiles_9.csv`
  - `results/amr/longitudinal_resident_bootstrap_inference.csv`
  - `results/amr/validation_checks.csv`
  - `results/amr/interpretation_report.md`
  - `results/summary/table_13_genomic_amr_summary.csv`
  - `results/vf_amr/vf_amr_input_availability_report.txt`
  - `results/vf_amr/vf_plasmid_combined_profile.csv`
  - `results/vf_amr/vf_amr_combined_profile_table.csv`
  - Plasmid/AMR-context plots.
- Why it matters:
  - This adds genomic resistance mechanism and longitudinal context to RQ01-RQ10
    without creating another research question.
- Guardrail:
  - `mdf(A)` is retained raw and in sensitivity metrics but excluded from
    primary acquired-gene burden, gain/loss and Jaccard calculations.
  - These are genomic predictions, not phenotypic AST.

### Step 34: Diagnose the UTI versus Not_UTI denominator

- Script/tool: `32_uti_not_uti_diagnostic_stats.R`
- What goes in:
  - Status map.
  - VF-ready data.
  - Endpoint/model outputs.
  - Exclusion, quarantine, and transition outputs.
- What happens:
  - Clinical and VF counts are validated.
  - Denominator flow and attrition are summarized.
  - Near-miss rows are audited.
  - Sparse-count uncertainty, leave-one-UTI sensitivity, paired participant deltas, and transition deltas are calculated.
- What comes out:
  - `results/audit/uti_not_uti_denominator_flow.csv`
  - `results/vf/uti_not_uti_attrition_summary.csv`
  - `results/vf/uti_not_uti_bootstrap_effects.csv`
  - `results/vf/uti_not_uti_leave_one_uti_out.csv`
  - Diagnostic plots and figure metadata.
- Why it matters:
  - It explains how much evidence the UTI versus Not_UTI contrast really has.

### Step 35: Build a mechanism-first casebook

- Script/tool: `33_mechanism_first_addon.R`
- What goes in:
  - Transition case outputs.
  - VF-ready data.
  - Status maps.
  - Module/endpoint changes.
  - Strain context.
  - Host context.
  - Variant annotations.
  - Panaroo accessory gene data.
  - Plasmid context and validated script-29 AMR profiles.
- What happens:
  - Not_UTI -> UTI transitions are organized into evidence buckets.
  - Host, strain, VF, module, accessory-gene, plasmid, variant, and validated genomic-AMR evidence are combined.
  - A readable casebook is written.
- What comes out:
  - `results/mechanism/not_uti_to_uti_casebook.csv`
  - `results/mechanism/not_uti_to_uti_casebook.md`
  - `results/mechanism/transition_mechanism_summary.csv`
  - Mechanism plots and validation checks.
- Why it matters:
  - This gives a structured explanation of possible mechanisms without overclaiming causality.

### Step 36: Run robustness checks

- Script/tool: `34_robustness_first_addon.R`
- What goes in:
  - Status map.
  - VF-ready data.
  - Endpoint tables.
  - Model outputs.
  - QC attrition files.
  - Module and transition outputs.
- What happens:
  - The script checks denominator robustness, near-miss sensitivity, model stability, leave-one-UTI behavior, bootstrap robustness, and sparse power/precision.
- What comes out:
  - `results/robustness/robustness_summary.md`
  - `results/robustness/robustness_claim_matrix.csv`
  - Robustness validation files and plots.
- Why it matters:
  - It helps decide which claims are stable enough to report and which need caution.

### Step 37: Run targeted statistical sensitivity checks

- Script/tool: `36_statistical_sensitivity_addon.R`
- What goes in:
  - Status map.
  - VF-ready data.
  - VF endpoint table.
  - Module outputs.
  - Transition outputs.
  - Model and lineage context.
- What happens:
  - Participant-collapsed supplementary endpoint tests are run.
  - Paired binary feature changes are assessed.
  - Transition module gain/loss enrichment is tested.
  - Endpoint-level GLMM/GLM sensitivity models are fitted where possible.
  - BH-adjusted q-values are written for new test families.
- What comes out:
  - Files under `results/statistical_sensitivity/`
  - Plots under `plots/statistical_sensitivity/`
- Why it matters:
  - It provides a more cautious statistical layer for the sparse UTI denominator.

### Step 38: Build final summary tables

- Script/tool: `30_vf_project_summary_tables.R`
- What goes in:
  - Clinical, WGS, MLST, VF, model, module, score, longitudinal, transition, lineage, diagnostic, mechanism, robustness, and statistical sensitivity outputs.
- What happens:
  - The script gathers upstream results.
  - It validates denominators and freshness.
  - It writes numbered summary tables.
  - It builds a figure index and visualization audit.
  - It writes a final key-results summary.
- What comes out:
  - `results/summary/table_01_*.csv` and later numbered tables.
  - `results/summary/final_key_results_summary.md`
  - `results/summary/final_summary_tables.rds`
  - `results/vf/vf_figure_index.csv`
  - `results/vf/vf_visualisation_audit.csv`
- Why it matters:
  - This creates the main reporting package from all earlier outputs.

### Step 39: Render the final figure pack

- Script/tool: `35_final_figure_pack.R`
- What goes in:
  - Validation checks.
  - Status map.
  - VF-ready data.
  - Mechanism casebook.
  - Robustness outputs.
  - Diagnostic outputs.
  - Statistical sensitivity outputs.
  - Variant tables.
  - Summary tables.
- What happens:
  - Required inputs are checked.
  - Four main figures and supplementary figures are rendered.
  - Captions and a figure manifest are written.
- What comes out:
  - `results/final_figures/final_figure_manifest.csv`
  - `results/final_figures/final_figure_captions.md`
  - Final PNG/PDF figures under `plots/final_figures/`
- Why it matters:
  - This is the final manuscript-style visual output layer.

### Step 40: Archive legacy outputs and verify final alignment

- Scripts/tools: `scripts/archive_legacy_asb_uti_outputs.R`, `scripts/verify_uti_not_uti_alignment.R`
- What goes in:
  - Current `results/` and `plots/` folders.
  - Status maps.
  - VF-ready dataset.
  - Summary outputs.
  - Quarantine and transition checks.
- What happens:
  - Old generated files that still advertise ASB-vs-UTI contrasts are moved to legacy archive locations.
  - The final status/VF alignment is checked.
  - Required primary UTI versus Not_UTI fields, denominators, curation rules, transition labels, and stale-output risks are checked.
- What comes out:
  - `results/qc/legacy_status_archive_manifest.csv`
  - `results/qc/uti_not_uti_alignment_checks.csv`
  - `results/qc/uti_not_uti_alignment_checks.txt`
- Why it matters:
  - This reduces the chance of accidentally using old ASB-vs-UTI files as current evidence.

## Simple End-to-End Flow

- Start with raw clinical data.
- Merge and clean the clinical data.
- Classify episodes as `UTI` or `Not_UTI`.
- Create display-only timepoint labels for event samples.
- Check assembly quality.
- Select canonical assemblies.
- Build core SNP distances.
- Run or refresh Panaroo pangenome outputs.
- Detect virulence factor genes.
- Build the VF presence/absence matrix.
- Add sequence type and plasmid context.
- Compare isolates from the same participant.
- Use 0-25 core-genome SNPs as strong same-strain evidence.
- Treat more than 25 SNPs as above the same-strain threshold.
- Treat missing SNPs as missing evidence.
- Use ST only as lineage context, not same-strain proof.
- Build the canonical VF-ready dataset.
- Run genotype-phenotype models and lineage checks.
- Build participant timelines.
- Identify Not_UTI -> UTI and other clinical transitions.
- Compare genomes for switch candidates.
- Parse and annotate variants.
- Add host context.
- Compare VF burden, genes, modules, and scores by status.
- Check longitudinal VF stability and gene gain/loss.
- Check lineage, denominator, robustness, and statistical sensitivity.
- Build mechanism casebooks.
- Produce final summary tables.
- Render final figures.
- Archive stale legacy outputs and verify current UTI-vs-Not_UTI alignment.

## Key Files to Remember

- `results/clinical/status_map.csv`
  - Master clinical status table.
- `results/clinical/status_map_with_poster_tp.csv`
  - Display-only timepoint map for visual ordering.
- `results/wgs/qc_summary.csv`
  - Assembly QC summary.
- `results/wgs/core/snp_dists.tsv`
  - Core SNP distance table.
- `results/vf/vf_pa_all.csv`
  - VF gene presence/absence matrix.
- `results/vf/vf_analysis_ready.csv`
  - Canonical VF plus clinical status analysis dataset.
- `results/vf/vf_binary_uti_ready.csv`
  - Model-ready VF/status dataset.
- `results/mlst/mlst_all.tsv`
  - Active MLST sequence type table.
- `results/plasmids/plasmidfinder_presence_absence.csv`
  - Plasmid replicon presence/absence table.
- `results/strain_compare/pairwise_metrics.csv`
  - Pairwise isolate comparison and same-strain context table.
- `results/longitudinal/participant_timelines.csv`
  - Participant timeline table.
- `results/longitudinal/phenotype_switch_candidates.csv`
  - Candidate clinical switches for deeper review.
- `results/longitudinal/variant_annotation_detailed.csv`
  - Variant-to-gene annotation table.
- `results/models/gwas_multivariable_glmm.csv`
  - Main participant-aware genotype-phenotype model output.
- `results/mechanism/not_uti_to_uti_casebook.md`
  - Mechanism-first transition casebook.
- `results/summary/final_key_results_summary.md`
  - Final key-results summary.
- `results/final_figures/final_figure_manifest.csv`
  - Final figure manifest.

## Main Guardrails

- The current primary clinical comparison is `UTI` versus `Not_UTI`.
- Legacy ASB/UTI/Negative fields are retained for comparison, not as the primary outcome.
- The same-strain threshold is 0-25 core-genome SNPs.
- Same ST does not prove same strain.
- Missing SNP evidence remains missing evidence.
- VF gene presence is genomic detection, not expression or causality.
- Display-only Uricult/poster timepoints are for visual explanation, not modelling.
- Many VF status tests are exploratory because UTI counts are sparse and repeated participant episodes are not independent.
- Participant-aware GLMMs and sensitivity analyses are used where possible, but sparse-data caveats still apply.
- Script 29 requires complete genomic-AMR profiles for all 532 assemblies;
  genomic determinants and predicted phenotypes are not measured susceptibility.

## One-Sentence Version

- The pipeline cleans and classifies clinical episodes, processes bacterial genomes, detects virulence genes, plasmids, and lineages, compares repeated isolates using SNP evidence, reconstructs timelines, investigates Not_UTI -> UTI switches, checks VF and mechanism evidence with robustness/sensitivity layers, and produces validated final tables and figures.

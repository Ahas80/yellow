#!/usr/bin/env Rscript
# ==============================================================================
# 30_vf_project_summary_tables.R
# ==============================================================================
#
# GOAL:
#   Generate final thesis/manuscript-ready summary tables for the VF-focused
#   YELLOW RoUTIne/rUTI analysis. This script collects, validates, summarises,
#   and exports key clinical, WGS, MLST, VF, module, score, longitudinal,
#   transition, lineage, and optional VF+AMR/plasmid outputs.
#
# METHOD:
#   1. Load canonical pipeline outputs and optional extension outputs.
#   2. Validate denominators, freshness, status labels, and missingness.
#   3. Build final numbered CSV tables plus optional XLSX/RDS bundles.
#   4. Write a markdown key-results summary using current generated outputs.
#
# OUTPUT:
#   - results/summary/table_01_cohort_episode_flow.csv
#   - results/summary/table_02_clinical_status_counts.csv
#   - results/summary/table_03_wgs_qc_summary.csv
#   - results/summary/table_04_mlst_st_summary.csv
#   - results/summary/table_05_vf_gene_category_summary.csv
#   - results/summary/table_06_vf_module_summary.csv
#   - results/summary/table_07_vf_score_summary_by_status.csv
#   - results/summary/table_08_vf_score_summary_by_ST.csv
#   - results/summary/table_09_longitudinal_vf_stability_summary.csv
#   - results/summary/table_10_not_uti_uti_transition_cases.csv
#   - results/summary/table_11_lineage_context_summary.csv
#   - results/summary/table_12_missing_data_audit.csv
#   - results/summary/table_13_optional_vf_amr_summary.csv
#   - results/summary/table_14_uti_not_uti_diagnostics.csv
#   - results/summary/final_key_results_summary.md
#   - results/summary/final_summary_tables.xlsx [optional]
#   - results/summary/final_summary_tables.rds  [optional]
#   - results/summary/summary_qc_log.txt
#
# ==============================================================================

source("00_config.R")
suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
  library(stringr)
  library(tibble)
})

msg("Starting 30_vf_project_summary_tables.R")

DIR_SUMMARY <- file.path(DIR_RESULTS, "summary")
ensure_dir(DIR_SUMMARY)

# ==============================================================================
# 1. HELPERS
# ==============================================================================
safe_read <- function(path, ...) {
  if (file.exists(path)) {
    msg("  Loading %s", basename(path))
    read_csv(path, show_col_types = FALSE, ...)
  } else {
    msg("  MISSING: %s", path)
    NULL
  }
}

mtime_chr <- function(path) {
  if (!file.exists(path)) return(NA_character_)
  format(file.info(path)$mtime, "%Y-%m-%d %H:%M")
}

is_stale <- function(target, upstream) {
  file.exists(target) && file.exists(upstream) &&
    file.info(target)$mtime < file.info(upstream)$mtime
}

vf_ready_matches_pa <- function(vf_ready_df, vf_pa_df) {
  if (is.null(vf_ready_df) || is.null(vf_pa_df)) return(FALSE)
  ready_genes <- canonical_vf_gene_cols(names(vf_ready_df), required = FALSE)
  pa_genes <- canonical_vf_gene_cols(names(vf_pa_df), required = FALSE)
  nrow(vf_ready_df) == nrow(vf_pa_df) && setequal(ready_genes, pa_genes)
}

normalise_tp_label <- function(x) {
  normalise_timepoint_preserve_events(x)
}

standardise_episode_keys <- function(df, tp_col = NULL, genomics = FALSE) {
  if (is.null(df)) return(NULL)
  out <- df %>% mutate(Participant_id = as.character(Participant_id))
  if ("UTI_Status" %in% names(out) || "Infection_Status" %in% names(out)) {
    out <- prefer_primary_uti_status(out)
  }
  if (!is.null(tp_col) && tp_col %in% names(out)) {
    out <- out %>% mutate(tp_lab = normalise_tp_label(.data[[tp_col]]))
  } else if ("tp_lab" %in% names(out)) {
    out <- out %>% mutate(tp_lab = normalise_tp_label(tp_lab))
  }
  out <- apply_manual_sample_curation(out, context = "30_summary_tables") %>%
    { if (genomics) filter_primary_genomics(.) else filter_primary_analysis(.) }
  out
}

status_counts <- function(df, source_file = NA_character_) {
  if (is.null(df) || !"Infection_Status" %in% names(df)) {
    return(tibble(n_ASB = NA_integer_, n_UTI = NA_integer_, n_Not_UTI = NA_integer_,
                  n_Negative = NA_integer_, n_other_status = NA_integer_))
  }
  status <- if ("UTI_Status" %in% names(df)) df$UTI_Status else df$Infection_Status
  legacy <- if ("Infection_Status_legacy" %in% names(df)) df$Infection_Status_legacy else df$Infection_Status
  tibble(
    n_ASB = sum(legacy == "ASB", na.rm = TRUE),
    n_UTI = sum(status == "UTI", na.rm = TRUE),
    n_Not_UTI = sum(status == "Not_UTI", na.rm = TRUE),
    n_Negative = sum(legacy == "Negative", na.rm = TRUE),
    n_other_status = sum(is.na(status) | !status %in% c("UTI", "Not_UTI"), na.rm = TRUE)
  )
}

flow_row <- function(stage_order, stage_name, df, source_file, reason = NA_character_) {
  sc <- status_counts(df)
  tibble(
    stage_order = stage_order,
    stage_name = stage_name,
    n_episodes = if (!is.null(df)) nrow(df) else NA_integer_,
    n_participants = if (!is.null(df) && "Participant_id" %in% names(df)) n_distinct(df$Participant_id) else NA_integer_,
    n_ASB = sc$n_ASB,
    n_UTI = sc$n_UTI,
    n_Not_UTI = sc$n_Not_UTI,
    n_Negative = sc$n_Negative,
    n_other_status = sc$n_other_status,
    reason_for_drop_from_previous = reason,
    source_file = source_file
  )
}

missing_row <- function(category, df, source_file, impact, likely_reason = NA_character_) {
  sc <- status_counts(df)
  tibble(
    missingness_category = category,
    n_episodes = if (!is.null(df)) nrow(df) else NA_integer_,
    n_participants = if (!is.null(df) && "Participant_id" %in% names(df)) n_distinct(df$Participant_id) else NA_integer_,
    n_ASB = sc$n_ASB,
    n_UTI = sc$n_UTI,
    n_Not_UTI = sc$n_Not_UTI,
    n_Negative = sc$n_Negative,
    n_other_status = sc$n_other_status,
    likely_reason = likely_reason,
    impact_on_analysis = impact,
    source_file = source_file
  )
}

score_summary_from_table <- function(df, group_col, score_cols) {
  if (is.null(df) || !group_col %in% names(df)) return(tibble())
  score_cols <- intersect(score_cols, names(df))
  if (length(score_cols) == 0) return(tibble())
  df %>%
    select(Participant_id, all_of(group_col), all_of(score_cols)) %>%
    filter(!is.na(.data[[group_col]])) %>%
    pivot_longer(cols = all_of(score_cols), names_to = "score_name", values_to = "value") %>%
    group_by(.data[[group_col]], score_name) %>%
    summarise(
      n_episodes = sum(!is.na(value)),
      n_participants = n_distinct(Participant_id[!is.na(value)]),
      median = ifelse(all(is.na(value)), NA_real_, median(value, na.rm = TRUE)),
      q25 = ifelse(all(is.na(value)), NA_real_, quantile(value, 0.25, na.rm = TRUE)),
      q75 = ifelse(all(is.na(value)), NA_real_, quantile(value, 0.75, na.rm = TRUE)),
      mean = ifelse(all(is.na(value)), NA_real_, round(mean(value, na.rm = TRUE), 2)),
      sd = ifelse(all(is.na(value)), NA_real_, round(sd(value, na.rm = TRUE), 2)),
      min = ifelse(all(is.na(value)), NA_real_, min(value, na.rm = TRUE)),
      max = ifelse(all(is.na(value)), NA_real_, max(value, na.rm = TRUE)),
      .groups = "drop"
    )
}

# ==============================================================================
# 2. LOAD INPUTS
# ==============================================================================
path_status <- FILE_STATUS_MAP
path_status_poster <- FILE_STATUS_MAP_POSTER
path_vf_ready <- FILE_VF_READY
path_vf_pa <- FILE_VF_PA
path_canonical_selection <- file.path(DIR_QC, "canonical_assembly_selection.csv")
path_gene_map <- file.path(DIR_VF, "gene_map.csv")
path_qc_wgs <- file.path(DIR_WGS, "qc_summary.csv")
path_mlst_freq <- file.path(DIR_MLST, "ST_frequencies.csv")
path_mlst_meta <- FILE_MLST_CANONICAL
path_mod_summary <- file.path(DIR_VF, "vf_module_summary.csv")
path_mod_map <- file.path(DIR_VF, "gene_module_map.csv")
path_scores <- file.path(DIR_VF, "vf_score_table.csv")
path_scores_status <- file.path(DIR_VF, "vf_score_summary_by_status.csv")
path_scores_st <- file.path(DIR_VF, "vf_score_summary_by_ST.csv")
path_score_catalog <- file.path(DIR_VF, "vf_score_endpoint_catalog.csv")
path_expec_marker_defs <- file.path(DIR_VF, "vf_expec_marker_definitions.csv")
path_expec_marker_summary <- file.path(DIR_VF, "vf_expec_marker_summary_by_status.csv")
path_expec_marker_tests <- file.path(DIR_VF, "vf_expec_marker_tests.csv")
path_vf_trans <- file.path(DIR_VF, "vf_longitudinal_transitions.csv")
path_vf_same_strain_summary <- file.path(DIR_VF, "vf_same_strain_vf_stability_summary.csv")
path_vf_replacement_summary <- file.path(DIR_VF, "vf_replacement_vf_change_summary.csv")
path_vf_strain_context_summary <- file.path(DIR_VF, "vf_strain_context_by_transition_summary.csv")
path_vf_same_strain_by_ST <- file.path(DIR_VF, "vf_same_strain_by_ST_summary.csv")
path_case_index <- file.path(DIR_VF, "vf_transition_case_index.csv")
path_case_summary <- file.path(DIR_VF, "vf_transition_case_summary.csv")
path_strain_ctx <- file.path(DIR_VF, "vf_transition_strain_context.csv")
path_vf_amr_combined <- file.path(DIR_RESULTS, "vf_amr", "vf_amr_combined_profile_table.csv")
path_vf_plasmid <- file.path(DIR_RESULTS, "vf_amr", "vf_plasmid_combined_profile.csv")
path_vf_amr_status <- file.path(DIR_RESULTS, "vf_amr", "vf_amr_score_summary_by_status.csv")
path_denominator <- file.path(DIR_QC, "pipeline_denominator_summary.csv")
path_metadata_manifest <- file.path(DIR_QC, "metadata_input_manifest.csv")
path_metadata_crosswalk <- file.path(DIR_QC, "metadata_input_crosswalk_audit.csv")
path_fasta_audit <- file.path(DIR_QC, "fasta_usage_audit.csv")
path_uti_attrition <- file.path(DIR_QC, "uti_attrition_episode_level.csv")
path_core_stale <- file.path(DIR_WGS, "core", "core_snp_staleness_report.txt")
path_panaroo_stale <- file.path(DIR_WGS, "pan", "panaroo_staleness_report.txt")
path_panaroo_manifest <- file.path(DIR_WGS, "pan", "panaroo_input_manifest.csv")
path_model_warnings <- file.path(DIR_MODELS, "model_interpretation_warnings.txt")
path_qc_bias <- file.path(DIR_QC, "qc_selection_bias_by_status.csv")
path_vf_gap <- file.path(DIR_VF, "vf_gene_annotation_gap_report.csv")
path_model_univ <- file.path(DIR_MODELS, "gwas_univariable_stats.csv")
path_model_glmm <- file.path(DIR_MODELS, "gwas_multivariable_glmm.csv")
path_model_denom <- file.path(DIR_MODELS, "model_dataset_denominator.csv")
path_score_tests <- file.path(DIR_VF, "vf_score_tests_exploratory.csv")
path_mod_assign <- file.path(DIR_VF, "vf_module_assignment_audit.csv")
path_vf_amr_report <- file.path(DIR_RESULTS, "vf_amr", "vf_amr_input_availability_report.txt")
path_diag_summary <- file.path(DIR_VF, "uti_not_uti_diagnostic_summary.csv")
path_diag_bootstrap <- file.path(DIR_VF, "uti_not_uti_bootstrap_effects.csv")
path_diag_feature_fisher <- file.path(DIR_VF, "uti_not_uti_feature_fisher_exploratory.csv")
path_diag_leave_one <- file.path(DIR_VF, "uti_not_uti_leave_one_uti_out.csv")
path_diag_power <- file.path(DIR_VF, "uti_not_uti_power_precision_context.csv")
path_diag_paired_deltas <- file.path(DIR_VF, "uti_not_uti_paired_participant_deltas.csv")
path_diag_paired_tests <- file.path(DIR_VF, "uti_not_uti_paired_participant_tests.csv")
path_diag_transition_tests <- file.path(DIR_VF, "uti_not_uti_transition_score_tests.csv")
path_diag_interpretation <- file.path(DIR_VF, "uti_not_uti_test_interpretation_table.csv")
path_diag_fig_meta <- file.path(DIR_VF, "uti_not_uti_diagnostic_figure_metadata.csv")
path_diag_attrition <- file.path(DIR_VF, "uti_not_uti_attrition_summary.csv")
path_diag_flow <- file.path(DIR_RESULTS, "audit", "uti_not_uti_denominator_flow.csv")
path_diag_near_miss <- file.path(DIR_RESULTS, "audit", "uti_not_uti_near_miss_rows.csv")
path_diag_duplicate_qc <- file.path(DIR_RESULTS, "audit", "duplicate_culture_qc_31036.csv")
path_stat_score_values <- file.path(DIR_RESULTS, "statistical_sensitivity", "participant_collapsed_score_values.csv")
path_stat_score_tests <- file.path(DIR_RESULTS, "statistical_sensitivity", "participant_collapsed_score_tests.csv")
path_stat_feature_catalog <- file.path(DIR_RESULTS, "statistical_sensitivity", "paired_binary_feature_catalog.csv")
path_stat_feature_deltas <- file.path(DIR_RESULTS, "statistical_sensitivity", "paired_binary_feature_deltas.csv")
path_stat_feature_tests <- file.path(DIR_RESULTS, "statistical_sensitivity", "paired_binary_feature_sensitivity.csv")
path_stat_transition_module <- file.path(DIR_RESULTS, "statistical_sensitivity", "transition_module_gain_loss_enrichment.csv")
path_stat_glmm <- file.path(DIR_RESULTS, "statistical_sensitivity", "score_glmm_sensitivity.csv")
path_stat_validation <- file.path(DIR_RESULTS, "statistical_sensitivity", "statistical_sensitivity_validation_checks.csv")
path_stat_fig_meta <- file.path(DIR_RESULTS, "statistical_sensitivity", "statistical_sensitivity_figure_metadata.csv")

status_map <- standardise_episode_keys(safe_read(path_status), "Timepoint")
status_poster_path_selected <- tryCatch(
  select_primary_status_map(prefer_poster = TRUE, require_fresh = TRUE,
                            caller = "30_vf_project_summary_tables.R", quiet = TRUE),
  error = function(e) NA_character_
)
status_poster <- if (!is.na(status_poster_path_selected) && identical(status_poster_path_selected, path_status_poster)) {
  standardise_episode_keys(safe_read(path_status_poster), "Timepoint")
} else {
  msg("  Skipping stale or non-primary status_map_with_poster_tp.csv in summary tables.")
  NULL
}
poster_status_unsuitable <- !file.exists(path_status_poster) ||
  is_stale(path_status_poster, path_status) ||
  !status_map_has_primary_fields(path_status_poster)
canonical_selection <- safe_read(path_canonical_selection)
if (is.null(canonical_selection)) stop("Missing canonical assembly selection: ", path_canonical_selection)
assembler_col <- intersect(c("assembler", "Assembler"), names(canonical_selection))[1]
if (is.na(assembler_col) ||
    !all(c("selected_canonical", "QC_PASS", "Participant_id", "tp_lab") %in% names(canonical_selection))) {
  stop("Canonical assembly selection lacks Longcycler-primary selection fields: ", path_canonical_selection)
}
active_longcycler_keys <- canonical_selection %>%
  mutate(
    Participant_id = as.character(.data$Participant_id),
    tp_lab = normalise_tp_label(.data$tp_lab),
    selected_canonical = as_pipeline_bool(.data$selected_canonical),
    QC_PASS = as_pipeline_bool(.data$QC_PASS),
    active_assembler = str_to_lower(as.character(.data[[assembler_col]]))
  ) %>%
  filter(.data$selected_canonical %in% TRUE,
         .data$QC_PASS %in% TRUE,
         .data$active_assembler == "longcycler") %>%
  distinct(.data$Participant_id, .data$tp_lab)
if (nrow(active_longcycler_keys) == 0) {
  stop("Canonical manifest contains no selected QC-passing Longcycler episode keys.")
}
restrict_to_active_longcycler <- function(df) {
  if (is.null(df) || !all(c("Participant_id", "tp_lab") %in% names(df))) return(df)
  df %>% semi_join(active_longcycler_keys, by = c("Participant_id", "tp_lab"))
}

vf_pa <- standardise_episode_keys(safe_read(path_vf_pa), genomics = TRUE) %>%
  restrict_to_active_longcycler()
vf_ready <- standardise_episode_keys(safe_read(path_vf_ready), genomics = TRUE) %>%
  restrict_to_active_longcycler()
gene_map <- safe_read(path_gene_map)
qc_wgs <- safe_read(path_qc_wgs)
mlst_freq <- safe_read(path_mlst_freq)
mlst_meta <- standardise_episode_keys(safe_read(path_mlst_meta), "Timepoint", genomics = TRUE) %>%
  restrict_to_active_longcycler()
mod_summary <- safe_read(path_mod_summary)
mod_map <- safe_read(path_mod_map)
score_table <- standardise_episode_keys(safe_read(path_scores), genomics = TRUE) %>%
  restrict_to_active_longcycler()
scores_status <- safe_read(path_scores_status)
scores_st <- safe_read(path_scores_st)
score_catalog <- safe_read(path_score_catalog)
expec_marker_defs <- safe_read(path_expec_marker_defs)
expec_marker_summary <- safe_read(path_expec_marker_summary)
expec_marker_tests <- safe_read(path_expec_marker_tests)
vf_trans <- standardise_episode_keys(safe_read(path_vf_trans))
vf_same_strain_summary <- safe_read(path_vf_same_strain_summary)
vf_replacement_summary <- safe_read(path_vf_replacement_summary)
vf_strain_context_summary <- safe_read(path_vf_strain_context_summary)
vf_same_strain_by_ST <- safe_read(path_vf_same_strain_by_ST)
case_index <- standardise_episode_keys(safe_read(path_case_index))
case_summary <- standardise_episode_keys(safe_read(path_case_summary))
normalise_transition_flag <- function(df) {
  if (is.null(df)) return(NULL)
  if (!"is_not_uti_to_uti" %in% names(df)) {
    df$is_not_uti_to_uti <- if (all(c("from_status", "to_status") %in% names(df))) {
      df$from_status == "Not_UTI" & df$to_status == "UTI"
    } else {
      FALSE
    }
  }
  df
}
case_index <- normalise_transition_flag(case_index)
case_summary <- normalise_transition_flag(case_summary)
strain_ctx <- safe_read(path_strain_ctx)
vf_amr_combined <- standardise_episode_keys(safe_read(path_vf_amr_combined), genomics = TRUE)
vf_plasmid <- standardise_episode_keys(safe_read(path_vf_plasmid), genomics = TRUE)
vf_amr_status <- safe_read(path_vf_amr_status)
denominator_summary <- safe_read(path_denominator)
metadata_manifest <- safe_read(path_metadata_manifest)
metadata_crosswalk <- safe_read(path_metadata_crosswalk)
fasta_audit <- safe_read(path_fasta_audit)
uti_attrition <- safe_read(path_uti_attrition)
panaroo_manifest <- safe_read(path_panaroo_manifest)
qc_bias <- safe_read(path_qc_bias)
vf_gap <- safe_read(path_vf_gap)
model_denom <- standardise_episode_keys(safe_read(path_model_denom)) %>%
  restrict_to_active_longcycler()
diag_summary <- safe_read(path_diag_summary)
diag_bootstrap <- safe_read(path_diag_bootstrap)
diag_feature_fisher <- safe_read(path_diag_feature_fisher)
diag_leave_one <- safe_read(path_diag_leave_one)
diag_power <- safe_read(path_diag_power)
diag_paired_deltas <- safe_read(path_diag_paired_deltas)
diag_paired_tests <- safe_read(path_diag_paired_tests)
diag_transition_tests <- safe_read(path_diag_transition_tests)
diag_interpretation <- safe_read(path_diag_interpretation)
diag_fig_meta <- safe_read(path_diag_fig_meta)
diag_attrition <- safe_read(path_diag_attrition)
diag_flow <- safe_read(path_diag_flow)
diag_near_miss <- safe_read(path_diag_near_miss)
diag_duplicate_qc <- safe_read(path_diag_duplicate_qc)
stat_score_values <- safe_read(path_stat_score_values)
stat_score_tests <- safe_read(path_stat_score_tests)
stat_feature_catalog <- safe_read(path_stat_feature_catalog)
stat_feature_deltas <- safe_read(path_stat_feature_deltas)
stat_feature_tests <- safe_read(path_stat_feature_tests)
stat_transition_module <- safe_read(path_stat_transition_module)
stat_glmm <- safe_read(path_stat_glmm)
stat_validation <- safe_read(path_stat_validation)
stat_fig_meta <- safe_read(path_stat_fig_meta)

if (is.null(status_map)) stop("Missing required primary UTI status map: ", path_status)
if (is.null(vf_pa)) stop("Missing required raw VF P/A matrix: ", path_vf_pa)
if (is.null(vf_ready)) stop("Missing required canonical VF-ready table: ", path_vf_ready)
if (nrow(vf_pa) != nrow(active_longcycler_keys)) {
  stop("VF P/A matrix does not contain every active selected QC-pass Longcycler episode key.")
}
if (nrow(vf_ready) != nrow(active_longcycler_keys)) {
  stop("VF-ready table does not contain every active selected QC-pass Longcycler episode key.")
}

# ==============================================================================
# 3. TABLE 01: COHORT AND EPISODE FLOW
# ==============================================================================
uti_not_uti_model <- vf_ready %>% filter(Infection_Status %in% c("Not_UTI", "UTI"))
longitudinal_subset <- if (!is.null(vf_trans)) {
  vf_ready %>% semi_join(vf_trans %>% select(Participant_id), by = "Participant_id")
} else NULL
not_uti_uti_cases <- if (!is.null(case_index)) {
  case_index %>% filter(is_not_uti_to_uti %in% TRUE)
} else NULL
not_uti_uti_wgs_cases <- if (!is.null(not_uti_uti_cases)) not_uti_uti_cases %>% filter(has_vf_pair %in% TRUE) else NULL
vf_amr_complete <- if (!is.null(vf_amr_combined)) {
  vf_amr_combined %>% filter(amr_data_available %in% TRUE)
} else NULL

t01 <- bind_rows(
  flow_row(1, "All-batch clinical episodes", status_map, path_status, "Canonical clinical denominator"),
  flow_row(2, "Ordered/poster clinical episodes", status_poster, path_status_poster,
           "Fresh ordered file used only for timelines when available"),
  flow_row(3, "Raw VF presence/absence rows", vf_pa, path_vf_pa, "Sequenced/VFDB-called episode rows"),
  flow_row(4, "Active Longcycler-primary VF-ready episodes", vf_ready, path_vf_ready,
           "Selected QC-pass Longcycler VF rows joined to primary UTI status/ST by script 22"),
  flow_row(5, "UTI/Not_UTI VF modelling subset", uti_not_uti_model, path_vf_ready,
           "Rows with missing primary UTI status excluded"),
  flow_row(6, "Repeated-measures VF longitudinal subset", longitudinal_subset, path_vf_trans,
           "Participants represented in VF transition output"),
  flow_row(7, "Clinical Not_UTI->UTI transition cases", not_uti_uti_cases, path_case_index,
           "All ordered clinical Not_UTI->UTI transitions indexed by script 28"),
  flow_row(8, "WGS/VF-linked Not_UTI->UTI transition cases", not_uti_uti_wgs_cases, path_case_index,
           "Both transition endpoints have VF-ready rows"),
  flow_row(9, "Optional true VF+AMR complete cases", vf_amr_complete, path_vf_amr_combined,
           "Only populated if dedicated AMR screening outputs exist")
)
write_csv(t01, file.path(DIR_SUMMARY, "table_01_cohort_episode_flow.csv"))

# ==============================================================================
# 4. TABLE 02: CLINICAL STATUS COUNTS BY DATASET LAYER
# ==============================================================================
status_layer <- function(df, layer, source_file) {
  if (is.null(df) || !"Infection_Status" %in% names(df)) return(tibble())
  df %>%
    group_by(Infection_Status) %>%
    summarise(
      n_episodes = n(),
      n_participants = n_distinct(Participant_id),
      pct_episodes = round(100 * n_episodes / nrow(df), 1),
      .groups = "drop"
    ) %>%
    mutate(dataset_layer = layer, source_file = source_file) %>%
    select(dataset_layer, Infection_Status, n_episodes, n_participants, pct_episodes, source_file)
}

t02 <- bind_rows(
  status_layer(status_map, "all_batch_clinical", path_status),
  status_layer(status_poster, "ordered_poster_clinical", path_status_poster),
  status_layer(vf_ready, "vf_ready", path_vf_ready),
  status_layer(uti_not_uti_model, "uti_not_uti_vf_subset", path_vf_ready)
)
write_csv(t02, file.path(DIR_SUMMARY, "table_02_clinical_status_counts.csv"))

# ==============================================================================
# 5. TABLE 03: WGS QC SUMMARY
# ==============================================================================
if (!is.null(qc_wgs)) {
  metric_cols <- intersect(c("total_bp", "n_contigs", "N50", "GC_pct"), names(qc_wgs))
  qc_group <- qc_wgs %>%
    mutate(qc_group = case_when(
      "QC_PASS" %in% names(.) & QC_PASS %in% TRUE ~ "pass",
      "QC_PASS" %in% names(.) & QC_PASS %in% FALSE ~ "fail",
      TRUE ~ "all"
    ))
  metric_rows <- bind_rows(lapply(c("all", sort(unique(qc_group$qc_group))), function(grp) {
    df <- if (grp == "all") qc_group else qc_group %>% filter(qc_group == grp)
    bind_rows(lapply(metric_cols, function(metric) {
      vals <- df[[metric]]
      tibble(
        metric = metric,
        group = grp,
        n = sum(!is.na(vals)),
        median = median(vals, na.rm = TRUE),
        q25 = quantile(vals, 0.25, na.rm = TRUE),
        q75 = quantile(vals, 0.75, na.rm = TRUE),
        min = min(vals, na.rm = TRUE),
        max = max(vals, na.rm = TRUE),
        source_file = path_qc_wgs
      )
    }))
  }))
  count_rows <- tibble(
    metric = c("n_assemblies", "n_qc_pass", "n_qc_fail"),
    group = "all",
    n = c(nrow(qc_wgs),
          if ("QC_PASS" %in% names(qc_wgs)) sum(qc_wgs$QC_PASS %in% TRUE, na.rm = TRUE) else NA_integer_,
          if ("QC_PASS" %in% names(qc_wgs)) sum(qc_wgs$QC_PASS %in% FALSE, na.rm = TRUE) else NA_integer_),
    median = NA_real_, q25 = NA_real_, q75 = NA_real_, min = NA_real_, max = NA_real_,
    source_file = path_qc_wgs
  )
  t03 <- bind_rows(count_rows, metric_rows)
} else {
  t03 <- tibble(metric = "wgs_qc_unavailable", group = "all", n = NA_integer_,
                median = NA_real_, q25 = NA_real_, q75 = NA_real_,
                min = NA_real_, max = NA_real_, source_file = path_qc_wgs)
}
write_csv(t03, file.path(DIR_SUMMARY, "table_03_wgs_qc_summary.csv"))

# ==============================================================================
# 6. TABLE 04: MLST/ST SUMMARY
# ==============================================================================
if ("ST" %in% names(vf_ready)) {
  t04 <- vf_ready %>%
    mutate(ST_source = if ("ST_source" %in% names(.)) ST_source else NA_character_) %>%
    filter(!is.na(ST), ST != "-") %>%
    group_by(ST) %>%
    summarise(
      n_episodes_or_isolates = n(),
      n_participants = n_distinct(Participant_id),
      ST_sources = paste(sort(unique(na.omit(ST_source))), collapse = ";"),
      pct = round(100 * n() / nrow(vf_ready), 1),
      n_legacy_ASB = if ("Infection_Status_legacy" %in% names(.)) sum(Infection_Status_legacy == "ASB", na.rm = TRUE) else NA_integer_,
      n_UTI = sum(Infection_Status == "UTI", na.rm = TRUE),
      n_Not_UTI = sum(Infection_Status == "Not_UTI", na.rm = TRUE),
      n_legacy_Negative = if ("Infection_Status_legacy" %in% names(.)) sum(Infection_Status_legacy == "Negative", na.rm = TRUE) else NA_integer_,
      median_vf_burden = median(vf_count_total, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(desc(n_episodes_or_isolates)) %>%
    mutate(rank = row_number(), source_file = path_vf_ready)
} else if (!is.null(mlst_freq)) {
  t04 <- mlst_freq %>% mutate(source_file = path_mlst_freq)
} else {
  t04 <- tibble(ST = character(), source_file = path_mlst_freq)
}
write_csv(t04, file.path(DIR_SUMMARY, "table_04_mlst_st_summary.csv"))

# ==============================================================================
# 7. TABLE 05: VF GENE AND CATEGORY SUMMARY
# ==============================================================================
if (!is.null(gene_map) && !is.null(vf_ready)) {
  gene_cols <- canonical_vf_gene_cols(names(vf_ready), required = FALSE)
  gene_summary <- tibble(Gene = gene_cols) %>%
    left_join(gene_map %>% select(Gene, Category), by = "Gene") %>%
    mutate(Category = coalesce(Category, "Unassigned")) %>%
    rowwise() %>%
    mutate(
      summary_level = "gene",
      n_episodes_present = sum(vf_ready[[Gene]], na.rm = TRUE),
      pct_episodes_present = round(100 * n_episodes_present / nrow(vf_ready), 1),
      n_Not_UTI_present = sum(vf_ready[[Gene]][vf_ready$Infection_Status == "Not_UTI"], na.rm = TRUE),
      n_legacy_ASB_present = if ("Infection_Status_legacy" %in% names(vf_ready)) sum(vf_ready[[Gene]][vf_ready$Infection_Status_legacy == "ASB"], na.rm = TRUE) else NA_integer_,
      n_UTI_present = sum(vf_ready[[Gene]][vf_ready$Infection_Status == "UTI"], na.rm = TRUE),
      n_legacy_Negative_present = if ("Infection_Status_legacy" %in% names(vf_ready)) sum(vf_ready[[Gene]][vf_ready$Infection_Status_legacy == "Negative"], na.rm = TRUE) else NA_integer_,
      assignment_status = ifelse(Category == "Unassigned", "unassigned", "assigned")
    ) %>%
    ungroup() %>%
    select(summary_level, Gene, Category, n_episodes_present, pct_episodes_present,
           n_Not_UTI_present, n_UTI_present, n_legacy_ASB_present,
           n_legacy_Negative_present, assignment_status)

  category_summary <- gene_summary %>%
    group_by(Category) %>%
    summarise(
      summary_level = "category",
      Gene = NA_character_,
      n_episodes_present = NA_integer_,
      pct_episodes_present = NA_real_,
      n_Not_UTI_present = sum(n_Not_UTI_present > 0, na.rm = TRUE),
      n_UTI_present = sum(n_UTI_present > 0, na.rm = TRUE),
      n_legacy_ASB_present = sum(n_legacy_ASB_present > 0, na.rm = TRUE),
      n_legacy_Negative_present = sum(n_legacy_Negative_present > 0, na.rm = TRUE),
      assignment_status = ifelse(first(Category) == "Unassigned", "unassigned", "assigned"),
      n_genes_in_category = n(),
      .groups = "drop"
    ) %>%
    select(names(gene_summary), n_genes_in_category)

  t05 <- bind_rows(
    gene_summary %>% mutate(n_genes_in_category = NA_integer_),
    category_summary
  ) %>%
    mutate(source_file = paste(path_vf_ready, path_gene_map, sep = "; "))
} else {
  t05 <- tibble(summary_level = character(), source_file = paste(path_vf_ready, path_gene_map, sep = "; "))
}
write_csv(t05, file.path(DIR_SUMMARY, "table_05_vf_gene_category_summary.csv"))

# ==============================================================================
# 8. TABLES 06-08: MODULE AND SCORE SUMMARIES
# ==============================================================================
t06 <- if (!is.null(mod_summary)) mod_summary else tibble(module_id = character())
write_csv(t06, file.path(DIR_SUMMARY, "table_06_vf_module_summary.csv"))
if (!is.null(expec_marker_defs)) {
  write_csv(expec_marker_defs, file.path(DIR_SUMMARY, "table_07a_expec_marker_definitions.csv"))
}
if (!is.null(expec_marker_summary)) {
  write_csv(expec_marker_summary, file.path(DIR_SUMMARY, "table_07b_expec_marker_summary_by_status.csv"))
}
if (!is.null(expec_marker_tests)) {
  write_csv(expec_marker_tests, file.path(DIR_SUMMARY, "table_07c_expec_marker_tests.csv"))
}

score_cols <- c("expec_marker_count", "upec_system_count", "upec_system_fraction",
                "total_vf_count_all", "total_vf_count_curated",
                "total_vf_count_unassigned", "low_confidence_count")
t07 <- if (!is.null(scores_status)) scores_status else score_summary_from_table(score_table, "Infection_Status", score_cols)
t08 <- if (!is.null(scores_st)) scores_st else score_summary_from_table(score_table, "ST", score_cols)
write_csv(t07, file.path(DIR_SUMMARY, "table_07_vf_score_summary_by_status.csv"))
write_csv(t08, file.path(DIR_SUMMARY, "table_08_vf_score_summary_by_ST.csv"))

# ==============================================================================
# 9. TABLE 09: LONGITUDINAL VF STABILITY SUMMARY
# ==============================================================================
if (!is.null(vf_trans)) {
  vf_trans_context <- vf_trans
  if (!"snp_strain_context" %in% names(vf_trans_context)) {
    vf_trans_context$snp_strain_context <- if ("SNPs" %in% names(vf_trans_context)) {
      case_when(
        is.na(vf_trans_context$SNPs) ~ "Missing SNP evidence",
        vf_trans_context$SNPs <= strain_snp_threshold() ~ "Strong same strain",
        TRUE ~ "Above same-strain SNP threshold"
      )
    } else {
      "Missing SNP evidence"
    }
  }
  if (!"st_lineage_context" %in% names(vf_trans_context)) {
    if ("same_ST" %in% names(vf_trans_context)) {
      vf_trans_context$st_lineage_context <- case_when(
        vf_trans_context$same_ST %in% TRUE ~ "Same ST",
        vf_trans_context$same_ST %in% FALSE ~ "Different ST",
        TRUE ~ "Missing ST evidence"
      )
    } else {
      vf_trans_context$st_lineage_context <- "Missing ST evidence"
    }
  }
  if (!"pair_interpretation" %in% names(vf_trans_context)) {
    vf_trans_context$pair_interpretation <- if ("same_strain_evidence" %in% names(vf_trans_context)) {
      as.character(vf_trans_context$same_strain_evidence)
    } else {
      "Missing strain metrics"
    }
  }
  vf_trans_context$snp_strain_context[is.na(vf_trans_context$snp_strain_context)] <- "Missing SNP evidence"
  vf_trans_context$st_lineage_context[is.na(vf_trans_context$st_lineage_context)] <- "Missing ST evidence"
  vf_trans_context$pair_interpretation[is.na(vf_trans_context$pair_interpretation)] <- "Missing strain metrics"
  if (!"same_strain_snp_threshold" %in% names(vf_trans_context)) {
    vf_trans_context$same_strain_snp_threshold <- strain_snp_threshold()
  }

  t09 <- vf_trans_context %>%
    group_by(snp_strain_context, pair_interpretation, st_lineage_context, transition_type) %>%
    summarise(
      n_transitions = n(),
      n_participants = n_distinct(Participant_id),
      median_vf_jaccard = median(jaccard_similarity, na.rm = TRUE),
      q25_vf_jaccard = ifelse(all(is.na(jaccard_similarity)), NA_real_, quantile(jaccard_similarity, 0.25, na.rm = TRUE)),
      q75_vf_jaccard = ifelse(all(is.na(jaccard_similarity)), NA_real_, quantile(jaccard_similarity, 0.75, na.rm = TRUE)),
      n_zero_gene_change = sum(n_gained == 0 & n_lost == 0, na.rm = TRUE),
      n_any_gene_change = sum(n_gained > 0 | n_lost > 0, na.rm = TRUE),
      median_genes_gained = median(n_gained, na.rm = TRUE),
      median_genes_lost = median(n_lost, na.rm = TRUE),
      same_strain_snp_threshold = suppressWarnings(max(as.numeric(same_strain_snp_threshold), na.rm = TRUE)),
      strong_same_strain_snp_range = sprintf("0-%d SNPs", strain_snp_threshold()),
      above_same_strain_snp_rule = sprintf(">%d SNPs = above same-strain SNP threshold", strain_snp_threshold()),
      st_lineage_rule = "ST is secondary lineage context: Same ST, Different ST, or Missing ST evidence",
      replacement_likely_rule = "Different ST or pairwise Different when SNPs do not support same strain",
      missing_strain_metrics_rule = "Missing SNP evidence with no usable ST/pairwise context or missing WGS/VF endpoint",
      source_file = path_vf_trans,
      .groups = "drop"
    ) %>%
    mutate(
      same_strain_snp_threshold = ifelse(is.infinite(same_strain_snp_threshold), strain_snp_threshold(), same_strain_snp_threshold),
      snp_strain_context = factor(
        snp_strain_context,
        levels = c("Strong same strain", "Above same-strain SNP threshold",
                   "Missing SNP evidence")
      ),
      pair_interpretation = factor(
        pair_interpretation,
        levels = c("Strong same strain",
                   "Conflict: SNP same-strain but ST differs",
                   "Same lineage, not same strain by SNP",
                   "ST-consistent, SNP missing",
                   "Above same-strain SNP threshold",
                   "Missing SNP evidence",
                   "Replacement likely",
                   "Missing strain metrics")
      ),
      st_lineage_context = factor(
        st_lineage_context,
        levels = c("Same ST", "Different ST", "Missing ST evidence")
      )
    ) %>%
    arrange(snp_strain_context, pair_interpretation, st_lineage_context, desc(n_transitions), transition_type) %>%
    mutate(
      snp_strain_context = as.character(snp_strain_context),
      pair_interpretation = as.character(pair_interpretation),
      st_lineage_context = as.character(st_lineage_context)
    )
} else {
  t09 <- tibble(snp_strain_context = character(), pair_interpretation = character(),
                st_lineage_context = character(), transition_type = character(),
                source_file = path_vf_trans)
}
write_csv(t09, file.path(DIR_SUMMARY, "table_09_longitudinal_vf_stability_summary.csv"))
write_csv(t09, file.path(DIR_SUMMARY, "table_09_longitudinal_vf_stability.csv")) # backward-compatible alias

# ==============================================================================
# 10. TABLE 10: Not_UTI->UTI TRANSITION CASES
# ==============================================================================
if (!is.null(case_summary)) {
  t10 <- case_summary %>%
    filter(is_not_uti_to_uti %in% TRUE |
             str_detect(transition_type, "^Not_UTI->UTI$")) %>%
    select(any_of(c("case_id", "Participant_id", "from_tp", "to_tp", "from_status", "to_status",
                    "transition_type", "has_vf_pair", "has_module_pair", "has_score_pair",
                    "is_uricult_transition", "timing_caveat", "ST_from", "ST_to",
                    "same_ST", "SNPs", "AvgIdentity", "Pairwise_Classification",
                    "Pairwise_RuleUsed", "snp_strain_context", "st_lineage_context",
                    "pair_interpretation", "same_strain_evidence", "strain_context_level",
                    "same_strain_snp_threshold", "vf_jaccard",
                    "module_jaccard", "n_vf_genes_gained", "n_vf_genes_lost",
                    "n_modules_gained", "n_modules_lost",
                    "delta_expec_marker_count", "delta_upec_system_count",
                    "delta_upec_system_fraction", "delta_upec_module_score",
                    "case_class", "missing_data_note", "interpretation_short")))
} else if (!is.null(case_index)) {
  t10 <- case_index %>% filter(is_not_uti_to_uti %in% TRUE)
} else {
  t10 <- tibble(case_id = character())
}
write_csv(t10, file.path(DIR_SUMMARY, "table_10_not_uti_uti_transition_cases.csv"))
if (!is.null(case_summary)) {
  write_csv(case_summary, file.path(DIR_SUMMARY, "table_10_all_transition_context.csv"))
}

# ==============================================================================
# 11. TABLE 11: LINEAGE CONTEXT SUMMARY
# ==============================================================================
if ("ST" %in% names(vf_ready)) {
  top_st_summary <- vf_ready %>%
    filter(!is.na(ST), ST != "-") %>%
    group_by(ST) %>%
    summarise(
      n_episodes = n(),
      n_participants = n_distinct(Participant_id),
      n_Not_UTI = sum(Infection_Status == "Not_UTI", na.rm = TRUE),
      n_UTI = sum(Infection_Status == "UTI", na.rm = TRUE),
      median_vf_burden = median(vf_count_total, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    filter(n_episodes >= 5) %>%
    arrange(desc(n_episodes)) %>%
    mutate(
      question = "Are VF burden and primary UTI status distributed across common STs?",
      test_or_summary = "Descriptive common-ST summary",
      result_value = sprintf("median VF burden %.1f; Not_UTI=%d; UTI=%d", median_vf_burden, n_Not_UTI, n_UTI),
      p_value = NA_real_,
      interpretation = "Lineage context should be considered before interpreting VF-status differences",
      source_file = path_vf_ready
    ) %>%
    select(question, ST, test_or_summary, n_episodes, n_participants,
           n_Not_UTI, n_UTI, result_value, p_value, interpretation, source_file)
  t11 <- top_st_summary
} else {
  t11 <- tibble(question = "ST unavailable", source_file = path_vf_ready)
}
write_csv(t11, file.path(DIR_SUMMARY, "table_11_lineage_context_summary.csv"))

# ==============================================================================
# 12. TABLE 12: MISSING-DATA AUDIT
# ==============================================================================
status_for_join <- status_map %>%
  mutate(Participant_id = as.character(Participant_id),
         Episode_ID = as.character(Episode_ID),
         tp_lab = normalise_timepoint_preserve_events(tp_lab)) %>%
  select(Participant_id, Episode_ID, tp_lab, Infection_Status, Batch)
vf_keys <- vf_ready %>%
  mutate(Participant_id = as.character(Participant_id),
         Episode_ID = as.character(Episode_ID),
         tp_lab = normalise_timepoint_preserve_events(tp_lab)) %>%
  select(Participant_id, Episode_ID, tp_lab, Infection_Status, Batch, any_of("uricult_bridge_applied"))
vf_pa_keys <- vf_pa %>% select(Participant_id, tp_lab)

clinical_no_vf <- status_for_join %>%
  anti_join(
    vf_keys %>% filter(!is.na(Infection_Status), !is.na(Episode_ID)),
    by = c("Participant_id", "Episode_ID")
  )
raw_vf_no_status <- vf_ready %>%
  filter(is.na(Infection_Status)) %>%
  transmute(Participant_id, tp_lab, Infection_Status = NA_character_)
vf_missing_st <- vf_ready %>% filter(!"ST" %in% names(.) | is.na(ST) | ST == "-")
vf_missing_status <- vf_ready %>% filter(is.na(Infection_Status))
transition_missing_wgs <- if (!is.null(case_index)) {
  case_index %>%
    filter(is_not_uti_to_uti %in% TRUE, !(has_vf_pair %in% TRUE)) %>%
    transmute(Participant_id, tp_lab = to_tp, Infection_Status = To_Status)
} else NULL
transition_missing_snp <- if (!is.null(case_summary)) {
  case_summary %>%
    filter(has_vf_pair %in% TRUE, is.na(SNPs)) %>%
    transmute(Participant_id, tp_lab = to_tp, Infection_Status = to_status)
} else NULL

t12 <- bind_rows(
  missing_row("Clinical episodes without VF-ready row", clinical_no_vf, paste(path_status, path_vf_ready, sep = "; "),
              "Excluded from VF/module/score analyses", "No WGS/VFDB row or unmatched episode key"),
  missing_row("Raw VF rows without primary UTI status", raw_vf_no_status, paste(path_vf_pa, path_status, sep = "; "),
              "Clinical interpretation unavailable for these VF rows", "Unmatched participant/timepoint"),
  missing_row("VF-ready episodes missing ST", vf_missing_st, path_vf_ready,
              "Lineage-aware summaries incomplete", "Missing/non-typable MLST or join issue"),
  missing_row("VF-ready episodes missing primary UTI status", vf_missing_status, path_vf_ready,
              "Excluded from status-stratified VF analyses", "Clinical join missing"),
  missing_row("Not_UTI->UTI transitions missing WGS/VF endpoint", transition_missing_wgs, path_case_index,
              "Cannot calculate VF/module/score transition changes", "One or both endpoints not in VF-ready data"),
  missing_row("WGS-linked transition cases missing SNP distance", transition_missing_snp, path_case_summary,
              "Same-strain interpretation weaker", "No pairwise/evolution SNP metric"),
  tibble(
    missingness_category = "Dedicated true AMR screening output",
    n_episodes = NA_integer_,
    n_participants = NA_integer_,
    n_ASB = NA_integer_, n_UTI = NA_integer_, n_Negative = NA_integer_, n_other_status = NA_integer_,
    likely_reason = if (dir.exists(file.path(DIR_RESULTS, "amr"))) "AMR directory exists but parser/status unknown" else "No results/amr directory detected",
    impact_on_analysis = "Script 29 remains VF+plasmid/mobile-context only; no true VF+AMR claims",
    source_file = file.path(DIR_RESULTS, "amr")
  ),
  tibble(
    missingness_category = "primary status_map_with_poster_tp freshness",
    n_episodes = if (!is.null(status_poster)) nrow(status_poster) else NA_integer_,
    n_participants = if (!is.null(status_poster)) n_distinct(status_poster$Participant_id) else NA_integer_,
    n_ASB = status_counts(status_poster)$n_ASB,
    n_UTI = status_counts(status_poster)$n_UTI,
    n_Negative = status_counts(status_poster)$n_Negative,
    n_other_status = status_counts(status_poster)$n_other_status,
    likely_reason = if (!poster_status_unsuitable) "Fresh primary poster status map" else "Missing, older than status_map.csv, or lacks UTI_Status/UTI_binary",
    impact_on_analysis = "Fresh file needed for best Uricult/date-aware transition ordering",
    source_file = path_status_poster
  ),
  tibble(
    missingness_category = "vf_analysis_ready freshness",
    n_episodes = nrow(vf_ready),
    n_participants = n_distinct(vf_ready$Participant_id),
    n_ASB = status_counts(vf_ready)$n_ASB,
    n_UTI = status_counts(vf_ready)$n_UTI,
    n_Negative = status_counts(vf_ready)$n_Negative,
    n_other_status = status_counts(vf_ready)$n_other_status,
    likely_reason = if (is_stale(path_vf_ready, path_vf_pa) && !vf_ready_matches_pa(vf_ready, vf_pa)) {
      "Older than vf_pa_all.csv and row/gene content differs"
    } else if (is_stale(path_vf_ready, path_vf_pa)) {
      "Timestamp older than vf_pa_all.csv, but row/gene content matches"
    } else {
      "Fresh relative to vf_pa_all.csv"
    },
    impact_on_analysis = if (is_stale(path_vf_ready, path_vf_pa) && !vf_ready_matches_pa(vf_ready, vf_pa)) {
      "Re-run script 22 before trusting scripts 26-30"
    } else {
      "Canonical VF denominator current"
    },
    source_file = path_vf_ready
  )
)
write_csv(t12, file.path(DIR_SUMMARY, "table_12_missing_data_audit.csv"))

# ==============================================================================
# 13. TABLE 13: OPTIONAL VF+AMR SUMMARY
# ==============================================================================
if (!is.null(vf_amr_status)) {
  note_vec <- if ("note" %in% names(vf_amr_status)) vf_amr_status$note else rep("", nrow(vf_amr_status))
  t13 <- vf_amr_status %>%
    mutate(analysis_available = any(!str_detect(note_vec, "No dedicated AMR"), na.rm = TRUE),
           source_file = path_vf_amr_status)
} else if (!is.null(vf_amr_combined)) {
  t13 <- vf_amr_combined %>%
    summarise(
      analysis_available = any(amr_data_available %in% TRUE, na.rm = TRUE),
      metric = "vf_amr_combined_profile_table",
      group = "all",
      n = n(),
      summary_value = ifelse(any(amr_data_available %in% TRUE), "True AMR data present", "No true AMR data integrated"),
      interpretation = "VF+plasmid/mobile-context table only unless amr_data_available is TRUE",
      source_file = path_vf_amr_combined
    )
} else {
  t13 <- tibble(
    analysis_available = FALSE,
    metric = "true_vf_amr_analysis",
    group = "all",
    n = NA_integer_,
    summary_value = "No script 29 combined table detected",
    interpretation = "Optional VF+AMR analysis was not performed",
    source_file = path_vf_amr_combined
  )
}
write_csv(t13, file.path(DIR_SUMMARY, "table_13_optional_vf_amr_summary.csv"))
write_csv(t13, file.path(DIR_SUMMARY, "table_13_optional_vf_plasmid_summary.csv")) # backward-compatible alias

# ==============================================================================
# 14. TABLE 14: PRIMARY UTI/NOT_UTI DIAGNOSTIC LAYER
# ==============================================================================
t14 <- if (!is.null(diag_summary)) {
  diag_summary %>%
    mutate(source_file = path_diag_summary)
} else {
  tibble(
    metric = "uti_not_uti_diagnostic_summary",
    value = NA_real_,
    source = path_diag_summary,
    interpretation = "Diagnostic layer unavailable; run 32_uti_not_uti_diagnostic_stats.R.",
    source_file = path_diag_summary
  )
}
write_csv(t14, file.path(DIR_SUMMARY, "table_14_uti_not_uti_diagnostics.csv"))

if (!is.null(diag_interpretation)) {
  write_csv(diag_interpretation, file.path(DIR_SUMMARY, "table_15_uti_not_uti_test_interpretation.csv"))
}
if (!is.null(diag_fig_meta)) {
  write_csv(diag_fig_meta, file.path(DIR_SUMMARY, "table_16_uti_not_uti_diagnostic_figures.csv"))
}

# ==============================================================================
# 14b. TABLE 17: TARGETED STATISTICAL SENSITIVITY ADD-ON
# ==============================================================================
if (!is.null(stat_score_tests)) {
  write_csv(stat_score_tests, file.path(DIR_SUMMARY, "table_17_participant_collapsed_score_tests.csv"))
}
if (!is.null(stat_feature_tests)) {
  write_csv(stat_feature_tests, file.path(DIR_SUMMARY, "table_18_paired_binary_feature_sensitivity.csv"))
}
if (!is.null(stat_transition_module)) {
  write_csv(stat_transition_module, file.path(DIR_SUMMARY, "table_19_transition_module_gain_loss_enrichment.csv"))
}
if (!is.null(stat_glmm)) {
  write_csv(stat_glmm, file.path(DIR_SUMMARY, "table_20_score_glmm_sensitivity.csv"))
}
if (!is.null(stat_validation)) {
  write_csv(stat_validation, file.path(DIR_SUMMARY, "table_21_statistical_sensitivity_validation.csv"))
}

# ==============================================================================
# 15. VF FIGURE INDEX AND VISUALISATION AUDIT
# ==============================================================================
vf_denominator_label <- sprintf(
  "Active selected QC-pass Longcycler VF-ready isolates n=%d; participants n=%d; primary UTI n=%d; primary Not_UTI n=%d; missing/other primary status n=%d",
  nrow(vf_ready),
  n_distinct(vf_ready$Participant_id),
  status_counts(vf_ready)$n_UTI,
  status_counts(vf_ready)$n_Not_UTI,
  status_counts(vf_ready)$n_other_status
)

clinical_denominator_label <- sprintf(
  "Clinical episodes n=%d; participants n=%d; primary UTI n=%d; primary Not_UTI n=%d; missing/other primary status n=%d",
  nrow(status_map),
  n_distinct(status_map$Participant_id),
  status_counts(status_map)$n_UTI,
  status_counts(status_map)$n_Not_UTI,
  status_counts(status_map)$n_other_status
)

transition_denominator_label <- if (!is.null(vf_trans)) {
  sprintf("Consecutive within-resident VF comparisons n=%d; participants n=%d",
          nrow(vf_trans), n_distinct(vf_trans$Participant_id))
} else {
  "Consecutive within-resident VF comparisons unavailable"
}

model_denominator_label <- if (!is.null(model_denom)) {
  sprintf("UTI/Not_UTI model rows n=%d; participants n=%d; UTI n=%d; Not_UTI n=%d",
          nrow(model_denom), n_distinct(model_denom$Participant_id),
          status_counts(model_denom)$n_UTI, status_counts(model_denom)$n_Not_UTI)
} else {
  "UTI/Not_UTI model denominator unavailable; see results/models/model_dataset_denominator.csv"
}

amr_scope_denominator_label <- if (!is.null(vf_amr_combined)) {
  sprintf("VF+plasmid combined rows n=%d; participants n=%d; true AMR rows n=%d",
          nrow(vf_amr_combined), n_distinct(vf_amr_combined$Participant_id),
          if ("amr_data_available" %in% names(vf_amr_combined)) sum(vf_amr_combined$amr_data_available %in% TRUE, na.rm = TRUE) else 0)
} else {
  "VF+plasmid combined table unavailable"
}

figure_index <- tribble(
  ~module, ~script, ~figure_id, ~output_file, ~title, ~input_data, ~denominator, ~level_of_analysis, ~statistical_test, ~exploratory_or_confirmatory, ~interpretation_limitations,
  "Burden", "23_vf_cross_sectional.R", "vf_burden_by_status", "plots/vf/vf_burden_by_status.png", "Distribution of E. coli virulence factor burden by primary UTI status", path_vf_ready, vf_denominator_label, "isolate-level", "None shown", "Descriptive", "Repeated isolates, small UTI denominator, and ST/lineage/event-type confounding limit causal interpretation.",
  "Burden", "23_vf_cross_sectional.R", "vf_burden_participant_summary", "plots/vf/vf_burden_participant_summary.png", "Participant-level summary of E. coli VF burden by primary UTI status", path_vf_ready, vf_denominator_label, "participant-status summary", "None shown", "Descriptive", "Participant summaries reduce repeated-isolate weighting but do not model confounding.",
  "Burden", "23_vf_cross_sectional.R", "vf_burden_paired_uti_not_uti", "plots/vf/vf_burden_paired_uti_not_uti.png", "Paired participant-level VF burden in residents with UTI and Not_UTI isolates", path_vf_ready, vf_denominator_label, "paired participant-level UTI/Not_UTI comparison", "None shown", "Descriptive/exploratory", "Restricted to residents with both UTI and Not_UTI VF-ready isolates; UTI n is small.",
  "Gene prevalence", "23_vf_cross_sectional.R", "vf_top_gene_prevalence", "plots/vf/vf_top_gene_prevalence.png", "Most prevalent virulence factor genes among VF/WGS-linked E. coli isolates", path_vf_ready, vf_denominator_label, "isolate-level gene prevalence", "None shown", "Descriptive", "Prevalence ranking is not an association test and does not imply expression or function.",
  "Gene prevalence", "23_vf_cross_sectional.R", "vf_gene_prevalence_difference_uti_not_uti", "plots/vf/vf_gene_prevalence_difference_uti_not_uti.png", "Virulence factor genes with largest UTI/Not_UTI prevalence differences", path_vf_ready, vf_denominator_label, "isolate-level UTI/Not_UTI gene prevalence difference", "Fisher tests in companion CSV; BH q-values", "Exploratory", "Repeated isolates are not modelled; rare genes give unstable contrasts; lineage may confound differences.",
  "Gene prevalence", "23_vf_cross_sectional.R", "vf_gene_prevalence_heatmap", "plots/vf/vf_gene_prevalence_heatmap.png", "VF gene prevalence heatmap across primary UTI status", path_vf_ready, vf_denominator_label, "isolate-level gene prevalence heatmap", "None shown", "Descriptive", "Percentages are denominator-sensitive and should be interpreted with status counts.",
  "Model evidence bridge", "14_genotype_phenotype_model.R", "vf_gene_screening_vs_model_evidence", "plots/vf/vf_gene_screening_vs_model_evidence.png", "Exploratory VF gene screening versus participant-aware model evidence", paste(path_model_univ, path_model_glmm, sep = "; "), model_denominator_label, "gene/feature-level association screen", "Fisher screening plus GLMM/GLM model OR, CI, and FDR", "Exploratory/model diagnostic", "Nominal Fisher hits do not imply robust corrected or participant-aware model support; sparse fits and lineage/event structure limit inference.",
  "Category profiles", "23_vf_cross_sectional.R", "vf_category_burden_by_status", "plots/vf/vf_category_burden_by_status.png", "Virulence factor category burden across primary UTI status", path_vf_ready, vf_denominator_label, "isolate-level category burden", "Category Fisher/Wilcoxon tests in CSV; BH q-values", "Exploratory", "Categories are descriptive, not validated causal virulence scores.",
  "Category profiles", "26_vf_define_gene_modules.R", "module_gene_counts", "plots/vf/module_gene_counts.png", "Virulence factor genes per curated biological module", paste(path_gene_map, path_vf_ready, sep = "; "), vf_denominator_label, "gene-to-module curation", "None", "Descriptive", "Module definitions are analysis curation units, not validated UTI scores.",
  "Category profiles", "26_vf_define_gene_modules.R", "module_prevalence_by_status", "plots/vf/module_prevalence_by_status.png", "Virulence factor module prevalence across primary UTI status", path_vf_ready, vf_denominator_label, "isolate-level module presence", "None shown", "Descriptive/exploratory", "Repeated measures and lineage structure are not modelled.",
  "Category profiles", "26_vf_define_gene_modules.R", "vf_category_composition_by_status", "plots/vf/vf_category_composition_by_status.png", "Virulence factor category profiles across primary UTI status", path_vf_ready, vf_denominator_label, "isolate-level category composition", "None shown", "Descriptive", "Category composition does not imply causal virulence.",
  "Annotation confidence", "26_vf_define_gene_modules.R", "vf_module_assignment_confidence", "plots/vf/vf_module_assignment_confidence.png", "VF module assignment confidence and gene-map coverage", paste(path_mod_assign, path_vf_gap, sep = "; "), "Gene/module assignment audit rows; see vf_gene_annotation_gap_report.csv for current VF-matrix gap counts", "gene-to-module curation diagnostic", "None shown", "Diagnostic/descriptive", "Low-confidence, unassigned, and gene-map-absent VFDB hits should not be interpreted as validated virulence modules.",
  "Longitudinal stability", "24_vf_longitudinal_dynamics.R", "vf_same_strain_jaccard_by_transition", "plots/vf/vf_same_strain_jaccard_by_transition.png", "VF similarity across repeated isolates with strong same-strain evidence", paste(path_vf_trans, path_vf_same_strain_summary, sep = "; "), transition_denominator_label, "consecutive same-strain within-resident isolate-pair comparison", "None shown", "Descriptive", sprintf("Same-strain calls use <=%d SNPs; ST is secondary lineage context.", strain_snp_threshold()),
  "Longitudinal stability", "24_vf_longitudinal_dynamics.R", "vf_same_strain_gene_gain_loss", "plots/vf/vf_same_strain_gene_gain_loss.png", "VF gene gains and losses within strong same-strain repeated isolates", paste(path_vf_trans, path_vf_same_strain_summary, sep = "; "), transition_denominator_label, "same-strain within-resident gene gain/loss summary", "None shown", "Descriptive", "Gene gain/loss within low-SNP pairs can still reflect assembly or annotation differences and should be read as descriptive.",
  "Longitudinal stability", "24_vf_longitudinal_dynamics.R", "vf_jaccard_by_strain_context", "plots/vf/vf_jaccard_by_strain_context.png", "VF similarity by SNP context and secondary pair interpretation", paste(path_vf_trans, path_vf_strain_context_summary, sep = "; "), transition_denominator_label, "consecutive within-resident isolate-pair comparison stratified by SNP context and pair interpretation", "None shown", "Descriptive/diagnostic", sprintf("Primary same-strain interpretation uses <=%d SNPs; ST is secondary lineage context and must not be merged into same-strain evidence.", strain_snp_threshold()),
  "Longitudinal stability", "24_vf_longitudinal_dynamics.R", "vf_replacement_vs_same_strain_vf_change", "plots/vf/vf_replacement_vs_same_strain_vf_change.png", "VF change is separated for SNP-defined same-strain persistence versus replacement", paste(path_vf_trans, path_vf_same_strain_summary, path_vf_replacement_summary, sep = "; "), transition_denominator_label, "consecutive within-resident isolate-pair comparison", "None shown", "Descriptive/diagnostic", "Replacement-like pairs are interpreted separately from SNP-defined within-strain VF stability/change.",
  "Longitudinal stability", "24_vf_longitudinal_dynamics.R", "vf_same_strain_by_ST", "plots/vf/vf_same_strain_by_ST.png", "Strong same-strain VF stability by sequence type", paste(path_vf_trans, path_vf_same_strain_by_ST, sep = "; "), transition_denominator_label, "same-strain pair summary grouped by ST", "None shown", "Secondary lineage diagnostic", "ST summaries are used after SNP-defined same-strain classification to check lineage structure and sparse-ST confounding.",
  "Longitudinal stability", "24_vf_longitudinal_dynamics.R", "vf_jaccard_by_transition", "plots/vf/vf_jaccard_by_transition.png", "Within-resident virulence factor similarity across repeated E. coli isolates", path_vf_ready, transition_denominator_label, "consecutive within-resident isolate-pair comparison", "None shown", "Descriptive", "Intervals vary; transition-specific UTI comparisons are sparse.",
  "Longitudinal stability", "24_vf_longitudinal_dynamics.R", "vf_within_host_jaccard_distribution", "plots/vf/vf_within_host_jaccard_distribution.png", "Distribution of within-resident VF profile similarity", path_vf_ready, transition_denominator_label, "consecutive within-resident isolate-pair comparison", "None shown", "Descriptive", "High similarity supports stability but does not prove same-strain persistence without SNP/ANI context.",
  "Longitudinal stability", "24_vf_longitudinal_dynamics.R", "vf_jaccard_by_days_between_samples", "plots/vf/vf_jaccard_by_days_between_samples.png", "VF similarity versus time between repeated isolates", path_vf_ready, transition_denominator_label, "dated consecutive within-resident isolate-pair comparison", "None shown", "Descriptive", "Only pairs with parseable dates are shown.",
  "Longitudinal stability", "24_vf_longitudinal_dynamics.R", "vf_jaccard_same_vs_different_st", "plots/vf/vf_jaccard_same_vs_different_st.png", "Secondary lineage diagnostic: VF similarity by sequence-type consistency", path_vf_ready, transition_denominator_label, "consecutive pair stratified by ST agreement", "None shown", "Secondary lineage diagnostic", sprintf("ST agreement is reviewed after SNP-defined same-strain classification (<=%d SNPs) and does not alone prove same strain.", strain_snp_threshold()),
  "Longitudinal stability", "24_vf_longitudinal_dynamics.R", "vf_gene_gain_loss_consecutive_pairs", "plots/vf/vf_gene_gain_loss_consecutive_pairs.png", "VF gene gains and losses between repeated E. coli isolates", path_vf_ready, transition_denominator_label, "consecutive within-resident gene gain/loss summary", "None shown", "Descriptive", "Gain/loss can reflect replacement, assembly/calling differences, or true gene-content change.",
  "Transition case studies", "28_vf_transition_case_studies.R", "vf_transition_case_timeline", "plots/vf/vf_transition_case_timeline.png", "Clinical timelines for residents with Not_UTI-to-UTI transitions", paste(path_status, path_status_poster, path_vf_ready, sep = "; "), clinical_denominator_label, "clinical episode timeline", "None shown", "Descriptive", "Uricult ordering uses dates where available; fallback/display ordering must be cited.",
  "Transition case studies", "28_vf_transition_case_studies.R", "vf_transition_score_slopeplot", "plots/vf/vf_transition_score_slopeplot.png", "Supplementary VF endpoint changes across Not_UTI-to-UTI transitions", path_vf_ready, "WGS/VF-linked Not_UTI-to-UTI transition cases", "transition case-study pair", "None shown", "Descriptive", "Changes may reflect lineage replacement or technical differences and do not imply causality.",
  "Transition case studies", "28_vf_transition_case_studies.R", "vf_transition_module_change_heatmap", "plots/vf/vf_transition_module_change_heatmap.png", "VF module changes in Not_UTI-to-UTI transition cases", path_vf_ready, "WGS/VF-linked Not_UTI-to-UTI transition cases", "transition case-study pair", "None shown", "Descriptive", "Modules are descriptive and should be interpreted after SNP same-strain classification; ST is secondary lineage context.",
  "Transition case studies", "28_vf_transition_case_studies.R", "vf_transition_gene_gain_loss_tile", "plots/vf/vf_transition_gene_gain_loss_tile.png", "VF gene gains and losses in Not_UTI-to-UTI transition cases", path_vf_ready, "WGS/VF-linked Not_UTI-to-UTI transition cases", "transition case-study pair", "None shown", "Descriptive", "Gene gain/loss should be interpreted after SNP same-strain classification; ST is secondary lineage context.",
  "Transition case studies", "28_vf_transition_case_studies.R", "vf_transition_snp_vs_vf_jaccard", "plots/vf/vf_transition_snp_vs_vf_jaccard.png", "SNP distance versus VF similarity in Not_UTI-to-UTI transition cases", path_vf_ready, "Not_UTI-to-UTI cases with WGS/VF, SNP distance, and VF Jaccard", "transition case-study pair", "None shown", "Descriptive/diagnostic", "Low SNP distance plus high VF similarity supports persistence but does not establish VF causality for UTI.",
  "Transition case studies", "28_vf_transition_case_studies.R", "vf_not_uti_uti_transition_strain_context", "plots/vf/vf_not_uti_uti_transition_strain_context.png", "Not_UTI-to-UTI VF changes require SNP-first and ST-secondary interpretation", paste(path_case_summary, path_strain_ctx, sep = "; "), "Not_UTI-to-UTI cases with WGS/VF endpoints, SNP distance, and VF burden change", "transition case-study strain-context diagnostic", "None shown", "Descriptive/diagnostic", "VF burden changes should be read through SNP-defined same-strain evidence first; ST explains lineage context but does not prove same strain.",
  "Transition case studies", "28_vf_transition_case_studies.R", "vf_not_uti_uti_transition_case_classes", "plots/vf/vf_not_uti_uti_transition_case_classes.png", "Not_UTI-to-UTI transition case classes and missing genomic endpoints", path_case_summary, "All clinical Not_UTI-to-UTI transition cases, including missing WGS/VF endpoints", "transition case classification", "None shown", "Descriptive/diagnostic", "Missing WGS/VF endpoints remain part of the clinical transition denominator and should not be hidden.",
  "Confounding checks", "25_vf_lineage_vf_interaction.R", "vf_burden_by_st", "plots/vf/vf_burden_by_st.png", "Virulence factor burden varies by E. coli sequence type", path_vf_ready, vf_denominator_label, "isolate-level ST diagnostic", "Kruskal-Wallis in summary text", "Exploratory/diagnostic", "Repeated isolates are not modelled; sparse STs are filtered.",
  "Confounding checks", "25_vf_lineage_vf_interaction.R", "vf_burden_st_x_status", "plots/vf/vf_burden_st_x_status.png", "UTI/Not_UTI VF burden contrasts within E. coli sequence types", path_vf_ready, vf_denominator_label, "isolate-level within-ST diagnostic", "Wilcoxon in summary text when possible", "Exploratory/diagnostic", "Within-ST UTI counts are often too small for stable inference.",
  "Confounding checks", "25_vf_lineage_vf_interaction.R", "vf_st_composition_by_status", "plots/vf/vf_st_composition_by_status.png", "Primary UTI status distribution across E. coli sequence types", path_vf_ready, vf_denominator_label, "isolate-level ST composition diagnostic", "Fisher simulated in summary text", "Exploratory/diagnostic", "ST distribution differences can confound naive VF-status associations.",
  "Confounding checks", "25_vf_lineage_vf_interaction.R", "vf_batch_by_status", "plots/vf/vf_batch_by_status.png", "Batch structure across VF-ready primary UTI status groups", path_vf_ready, vf_denominator_label, "isolate-level batch diagnostic", "None shown", "Diagnostic", "Batch imbalance can affect interpretation of associations.",
  "Confounding checks", "25_vf_lineage_vf_interaction.R", "vf_event_type_by_status", "plots/vf/vf_event_type_by_status.png", "Sampling context across VF-ready primary UTI status groups", path_vf_ready, vf_denominator_label, "isolate-level event-type diagnostic", "None shown", "Diagnostic", "Routine and event-driven samples are not interchangeable denominators.",
  "Confounding checks", "25_vf_lineage_vf_interaction.R", "vf_status_timepoint_event_tile", "plots/vf/vf_status_timepoint_event_tile.png", "Primary UTI status, timepoint, and event context in VF-ready isolates", path_vf_ready, vf_denominator_label, "isolate-level timepoint and event-context diagnostic", "None shown", "Diagnostic", "UTI-labelled VF/WGS rows are structurally concentrated in event-driven labels and should not be read as routine numeric timepoints.",
  "Confounding checks", "25_vf_lineage_vf_interaction.R", "vf_qc_selection_by_status", "plots/vf/vf_qc_selection_by_status.png", "WGS/QC selection structure by primary UTI status", path_qc_bias, "QC/selection rows by primary UTI status from qc_selection_bias_by_status.csv", "QC and selection-bias diagnostic", "None shown", "Diagnostic", "Unequal QC or selection inclusion across primary UTI status can alter VF-ready denominators.",
  "Confounding checks", "25_vf_lineage_vf_interaction.R", "vf_denominator_flow", "plots/vf/vf_denominator_flow.png", "Clinical-to-genomic denominator flow for VF analysis", paste(path_status, path_vf_pa, path_vf_ready, sep = "; "), paste(clinical_denominator_label, vf_denominator_label, sep = " -> "), "denominator-flow diagnostic", "None", "Diagnostic", "Do not hide attrition from clinical episodes to VF-ready analyses.",
  "Confounding checks", "25_vf_lineage_vf_interaction.R", "vf_uricult_join_diagnostic", "plots/vf/vf_uricult_join_diagnostic.png", "Uricult clinical events require harmonisation with UTI-labelled WGS isolates", paste(path_status, path_vf_ready, sep = "; "), "Clinical UTI versus VF-ready UTI bridge counts", "join/denominator diagnostic", "None", "Diagnostic", "Uricult clinical labels and UTI-N WGS labels require explicit bridge assumptions.",
  "Score framework", "27_vf_score_framework.R", "vf_expec_marker_prevalence_by_status", "plots/vf/vf_expec_marker_prevalence_by_status.png", "ExPEC-like marker prevalence by primary UTI status", paste(path_expec_marker_summary, path_expec_marker_tests, sep = "; "), vf_denominator_label, "isolate-level marker prevalence", "Fisher exact with BH q-values", "Supplementary/exploratory", "Marker groups are literature-aligned but not validated UTI predictors in this cohort; repeated residents and lineage are not modelled.",
  "Score framework", "27_vf_score_framework.R", "vf_scores_by_status", "plots/vf/vf_scores_by_status.png", "Supplementary VF marker, system, and burden summaries by primary UTI status", path_vf_ready, vf_denominator_label, "isolate-level supplementary endpoint summary", "Wilcoxon in results/vf/vf_score_tests_exploratory.csv", "Supplementary/exploratory", "Composite endpoints are not validated predictors; repeated isolates and lineage are not modelled.",
  "Score framework", "27_vf_score_framework.R", "vf_score_effect_summary_uti_not_uti", "plots/vf/vf_score_effect_summary_uti_not_uti.png", "Supplementary UTI/Not_UTI VF endpoint differences", path_score_tests, vf_denominator_label, "isolate-level supplementary endpoint effect summary", "Wilcoxon with BH q-values", "Supplementary/exploratory", "Median endpoint differences and BH q-values are descriptive because participant clustering, ST/lineage, timepoint, batch, and event type are not modelled here.",
  "Score framework", "27_vf_score_framework.R", "vf_scores_by_ST", "plots/vf/vf_scores_by_ST.png", "Supplementary VF marker/system summaries by E. coli sequence type", path_vf_ready, vf_denominator_label, "isolate-level lineage diagnostic", "None shown", "Diagnostic/descriptive", "VF marker/system endpoints are expected to be lineage structured.",
  "Score framework", "27_vf_score_framework.R", "vf_score_correlation_heatmap", "plots/vf/vf_score_correlation_heatmap.png", "Spearman correlations among supplementary VF endpoints", path_vf_ready, vf_denominator_label, "endpoint correlation", "Spearman correlation", "Descriptive", "Correlations reflect shared marker/module components and do not indicate clinical prediction.",
  "Score framework", "27_vf_score_framework.R", "vf_pca_status", "plots/vf/vf_pca_status.png", "Exploratory PCA of VF module profiles by primary UTI status", path_vf_ready, vf_denominator_label, "ordination", "PCA", "Exploratory/descriptive", "Ordination is not adjusted for repeated residents or lineage.",
  "Score framework", "27_vf_score_framework.R", "vf_pca_ST", "plots/vf/vf_pca_ST.png", "Exploratory PCA of VF module profiles by sequence type", path_vf_ready, vf_denominator_label, "ordination", "PCA", "Diagnostic/descriptive", "Lineage clustering should be considered as confounding.",
  "Score framework", "27_vf_score_framework.R", "vf_pcoa_jaccard_status", "plots/vf/vf_pcoa_jaccard_status.png", "Exploratory Jaccard PCoA of VF module profiles by primary UTI status", path_vf_ready, vf_denominator_label, "ordination", "Jaccard PCoA", "Exploratory/descriptive", "Ordination is descriptive and denominator-sensitive.",
  "Score framework", "27_vf_score_framework.R", "vf_pcoa_jaccard_ST", "plots/vf/vf_pcoa_jaccard_ST.png", "Exploratory Jaccard PCoA of VF module profiles by sequence type", path_vf_ready, vf_denominator_label, "ordination", "Jaccard PCoA", "Diagnostic/descriptive", "Lineage clustering should be considered as confounding.",
  "Optional VF+plasmid", "29_vf_amr_combined_profile.R", "vf_plasmid_analysis_scope", "plots/vf_amr/vf_plasmid_analysis_scope.png", "Scope of VF, plasmid, and AMR data integration", paste(path_vf_amr_report, path_vf_amr_combined, sep = "; "), amr_scope_denominator_label, "input availability and analysis-scope diagnostic", "None shown", "Diagnostic", "Script 29 is VF+plasmid/mobile-context only unless dedicated AMR screening rows are present; plasmid summaries are not true AMR analysis."
) %>%
  mutate(
    output_exists = file.exists(file.path(DIR_ROOT, output_file)),
    note = ifelse(output_exists, "Generated or expected from current pipeline", "Optional/skipped if required inputs were absent")
  )

if (!is.null(diag_fig_meta) && nrow(diag_fig_meta) > 0) {
  diagnostic_title <- c(
    uti_not_uti_clinical_rule_flow = "Primary UTI / Not_UTI decision-flow diagnostic",
    uti_not_uti_denominator_waterfall = "Clinical-to-VF/model denominator diagnostic",
    uti_not_uti_near_miss_evidence_heatmap = "Near-miss legacy UTI evidence heatmap",
    uti_not_uti_bootstrap_effect_forest = "Participant-bootstrap VF endpoint effect sizes",
    uti_not_uti_leave_one_uti_out_stability = "Leave-one-UTI-out stability diagnostic",
    uti_not_uti_paired_participant_slopeplot = "Within-resident paired UTI / Not_UTI VF endpoint slopes",
    uti_not_uti_transition_delta_forest = "Transition endpoint-change diagnostic",
    uti_not_uti_sparse_power_precision = "Sparse-UTI precision context",
    duplicate_culture_qc_31036 = "Duplicate-culture QC for 31036 UTI-1 versus UTI-2"
  )
  diag_rows <- diag_fig_meta %>%
    mutate(
      output_file = if_else(
        str_starts(.data$file_path, fixed(DIR_ROOT)),
        substring(.data$file_path, nchar(DIR_ROOT) + 2L),
        .data$file_path
      )
    ) %>%
    transmute(
      module = "Primary UTI/Not_UTI diagnostics",
      script = "32_uti_not_uti_diagnostic_stats.R",
      figure_id,
      output_file,
      title = recode(.data$figure_id, !!!diagnostic_title, .default = .data$figure_id),
      input_data = paste(path_status, path_vf_ready, path_model_denom, path_diag_summary, sep = "; "),
      denominator = sprintf(
        "Clinical primary n=583 (UTI 18, Not_UTI 565); active selected QC-pass Longcycler VF/model primary n=%d (UTI %d, Not_UTI %d, missing status %d)",
        nrow(vf_ready), status_counts(vf_ready)$n_UTI,
        status_counts(vf_ready)$n_Not_UTI, status_counts(vf_ready)$n_other_status
      ),
      level_of_analysis = .data$evidence_type,
      statistical_test = case_when(
        str_detect(.data$figure_id, "bootstrap") ~ "Participant bootstrap confidence intervals",
        str_detect(.data$figure_id, "leave_one") ~ "Leave-one-UTI-out stability",
        str_detect(.data$figure_id, "paired") ~ "Participant-paired deltas with sign/signed-rank summaries in CSV",
        str_detect(.data$figure_id, "transition") ~ "Transition sign/signed-rank summaries in CSV",
        str_detect(.data$figure_id, "power") ~ "Expected Fisher-test sparse-count context",
        TRUE ~ "None shown"
      ),
      exploratory_or_confirmatory = case_when(
        str_detect(.data$evidence_type, regex("bootstrap|paired|transition|context", ignore_case = TRUE)) ~
          "Exploratory/diagnostic",
        TRUE ~ "Descriptive/diagnostic"
      ),
      interpretation_limitations,
      output_exists = file.exists(.data$file_path),
      note = ifelse(.data$output_exists, "Generated by script 32", "Run script 32 to regenerate")
    )
  figure_index <- bind_rows(figure_index, diag_rows)
}

visualisation_audit <- figure_index %>%
  transmute(
    module,
    script_name = script,
    approximate_section_or_line = case_when(
      script == "23_vf_cross_sectional.R" ~ "Section 7, VF visualisation modules 01-03",
      script == "14_genotype_phenotype_model.R" ~ "Model plotting section, evidence bridge",
      script == "24_vf_longitudinal_dynamics.R" ~ "Section 5, longitudinal stability plots",
      script == "25_vf_lineage_vf_interaction.R" ~ "Section 5, confounding diagnostics",
      script == "26_vf_define_gene_modules.R" ~ "Section 11, module/category plots",
      script == "27_vf_score_framework.R" ~ "Section 10, score/ordination plots",
      script == "28_vf_transition_case_studies.R" ~ "Section 10, transition figures",
      script == "29_vf_amr_combined_profile.R" ~ "Section 8, VF+plasmid plots",
      script == "32_uti_not_uti_diagnostic_stats.R" ~ "Primary UTI/Not_UTI diagnostic layer",
      TRUE ~ "See script"
    ),
    plot_object_name = figure_id,
    ggsave_output_path = output_file,
    input_dataset = input_data,
    level_of_analysis,
    current_title = title,
    current_subtitle = "See ggplot subtitle/caption in current script",
    current_axis_labels = "Scientifically labelled in current script",
    current_legend_labels = "Scientifically labelled in current script where legend is present",
    current_caption = "Caption includes input data, denominator, analysis level, repeated-measures caveat, exploratory status, UTI denominator, and lineage/batch/event caveats where relevant",
    denominator_shown_or_missing = "Shown in subtitle/caption and/or figure index",
    statistical_annotation_shown_or_missing = statistical_test,
    repeated_measure_acknowledged = str_detect(interpretation_limitations, regex("repeated|resident", ignore_case = TRUE)) | module %in% c("Longitudinal stability", "Transition case studies"),
    suitable_for_thesis_manuscript_poster_use = output_exists,
    classification = case_when(
      figure_id %in% c("vf_not_uti_uti_transition_strain_context", "vf_not_uti_uti_transition_case_classes") ~
        "ADD companion plot / MARK diagnostic",
      module == "Primary UTI/Not_UTI diagnostics" ~ "MARK diagnostic",
      module %in% c("Confounding checks") ~ "ADD companion plot / MARK diagnostic",
      module %in% c("Annotation confidence", "Optional VF+plasmid") ~ "ADD companion plot / MARK diagnostic",
      str_detect(exploratory_or_confirmatory, regex("Exploratory", ignore_case = TRUE)) ~ "MARK exploratory",
      TRUE ~ "AMEND labels/caption only"
    ),
    problem_found = case_when(
      module == "Confounding checks" ~ "Needed explicit denominator, ST, batch, event-type, and Uricult-bridge diagnostic framing.",
      module == "Model evidence bridge" ~ "Nominal Fisher screening hits needed a direct visual comparison with corrected participant-aware model evidence.",
      module == "Annotation confidence" ~ "Module plots needed an explicit annotation-confidence and gene-map coverage diagnostic.",
      module == "Optional VF+plasmid" ~ "VF+plasmid outputs needed a scope figure preventing interpretation as true AMR integration.",
      module == "Primary UTI/Not_UTI diagnostics" ~ "Needed explicit sparse-count, denominator, exclusion, duplicate, and near-miss diagnostics under the primary UTI rule.",
      module == "Longitudinal stability" ~ "Older logic risked omitting T5/T6 and UTI-N event labels and needed clearer comparison scope.",
      module == "Gene prevalence" ~ "Gene selection and multiple-testing/exploratory status needed explicit wording.",
      module == "Burden" ~ "Burden definition and repeated-measures/UTI denominator caveats needed to be visible.",
      module == "Category profiles" ~ "Category/module definitions could be read as causal scores without clearer captions.",
      TRUE ~ "Needed clearer exploratory/descriptive interpretation."
    ),
    amendment_made = "Updated titles/subtitles/captions and, where relevant, added companion plots or tables using vf_analysis_ready.csv.",
    why_it_matters = interpretation_limitations,
    output_exists,
    note
  )

write_csv(figure_index, file.path(DIR_VF, "vf_figure_index.csv"))
write_csv(visualisation_audit, file.path(DIR_VF, "vf_visualisation_audit.csv"))

# ==============================================================================
# 15. FINAL KEY RESULTS MARKDOWN + QC LOG
# ==============================================================================
qc_log <- character()
qa <- function(...) qc_log <<- c(qc_log, sprintf(...))
qa("=== SUMMARY QC LOG ===")
qa("Generated: %s", format(Sys.time()))
qa("status_map.csv modified: %s", mtime_chr(path_status))
qa("status_map_with_poster_tp.csv modified: %s", mtime_chr(path_status_poster))
qa("vf_pa_all.csv modified: %s", mtime_chr(path_vf_pa))
qa("vf_analysis_ready.csv modified: %s", mtime_chr(path_vf_ready))
qa("uti_not_uti_diagnostic_summary.csv modified: %s", mtime_chr(path_diag_summary))
qa("status_map_with_poster_tp unsuitable for primary summaries: %s", poster_status_unsuitable)
qa("vf_analysis_ready stale by timestamp: %s", is_stale(path_vf_ready, path_vf_pa))
qa("vf_analysis_ready row/gene content matches vf_pa_all: %s", vf_ready_matches_pa(vf_ready, vf_pa))
if (!is.null(diag_summary)) {
  diag_counts <- diag_summary %>%
    filter(.data$metric %in% c(
      "primary_clinical_uti", "primary_clinical_not_uti",
      "primary_vf_model_uti", "primary_vf_model_not_uti",
      "primary_vf_model_missing_status"
    )) %>%
    transmute(label = paste0(.data$metric, "=", .data$value)) %>%
    pull(.data$label)
  qa("script_32_diagnostic_counts: %s", paste(diag_counts, collapse = "; "))
}
if (file.exists(path_core_stale)) {
  core_status <- grep("^Status:", readLines(path_core_stale, warn = FALSE), value = TRUE)
  qa("core_snp_status: %s", paste(core_status, collapse = " | "))
}
if (file.exists(path_panaroo_stale)) {
  pan_status <- grep("^Status:", readLines(path_panaroo_stale, warn = FALSE), value = TRUE)
  qa("panaroo_status: %s", paste(pan_status, collapse = " | "))
}
if (file.exists(path_model_warnings)) {
  qa("model_interpretation_warnings:")
  for (ln in readLines(path_model_warnings, warn = FALSE)) qa("  %s", ln)
}

legacy_paths <- c(
  file.path(DIR_ROOT, "KEY_FINDINGS.md"),
  file.path(DIR_ROOT, "FINAL_SUMMARY.md"),
  file.path(DIR_ROOT, "docs", "FINAL_SUMMARY.md"),
  file.path(DIR_ROOT, "docs", "research_outcomes.md"),
  file.path(DIR_RESULTS, "status_map.csv"),
  file.path(DIR_VF, "vf_analysis_ready_OLD.csv"),
  file.path(DIR_VF, "vf_burden_by_status_OLD.csv"),
  file.path(DIR_VF, "vf_longitudinal_transitions_OLD.csv")
)
legacy_found <- legacy_paths[file.exists(legacy_paths)]
if (length(legacy_found) > 0) {
  qa("Legacy/stale-looking files detected; do not use as truth unless manually verified:")
  for (p in legacy_found) qa("  - %s", p)
}
writeLines(qc_log, file.path(DIR_SUMMARY, "summary_qc_log.txt"))

md <- character()
ma <- function(...) md <<- c(md, sprintf(...))

vf_gene_cols <- canonical_vf_gene_cols(names(vf_ready), required = FALSE)
raw_vf_gene_cols <- canonical_vf_gene_cols(names(vf_pa), required = FALSE)

ma("# YELLOW RoUTIne / rUTI VF-Focused Key Results Summary")
ma("")
ma("**Generated:** %s", format(Sys.time(), "%Y-%m-%d %H:%M"))
ma("")
ma("## Cohort Availability")
ma("- Clinical episodes: **%d** (%d participants)", nrow(status_map), n_distinct(status_map$Participant_id))
ma("- Clinical primary status counts: UTI **%d**, Not_UTI **%d**, Other/unknown primary status **%d**",
   status_counts(status_map)$n_UTI, status_counts(status_map)$n_Not_UTI,
   status_counts(status_map)$n_other_status)
if (!is.null(status_poster)) {
  ma("- Ordered/poster clinical episodes: **%d** (%d participants)", nrow(status_poster), n_distinct(status_poster$Participant_id))
}
ma("")
ma("## VF/WGS-Ready Dataset")
ma("- Raw VF matrix rows: **%d** (%d participants)", nrow(vf_pa), n_distinct(vf_pa$Participant_id))
ma("- Raw VF gene columns: **%d**", length(raw_vf_gene_cols))
ma("- Active selected QC-pass Longcycler VF-ready episodes: **%d** (%d participants)", nrow(vf_ready), n_distinct(vf_ready$Participant_id))
ma("- VF-ready primary status counts: UTI **%d**, Not_UTI **%d**",
   status_counts(vf_ready)$n_UTI, status_counts(vf_ready)$n_Not_UTI)
ma("- VF-ready gene columns: **%d**", length(vf_gene_cols))
if ("ST" %in% names(vf_ready)) {
  ma("- Distinct STs in VF-ready dataset: **%d**", n_distinct(vf_ready$ST[!is.na(vf_ready$ST) & vf_ready$ST != "-"]))
}
ma("- Median total VF burden: **%.0f** (IQR %.0f-%.0f)",
   median(vf_ready$vf_count_total, na.rm = TRUE),
   quantile(vf_ready$vf_count_total, 0.25, na.rm = TRUE),
   quantile(vf_ready$vf_count_total, 0.75, na.rm = TRUE))
ma("")
ma("## Primary UTI/Not_UTI Diagnostic Layer")
if (!is.null(diag_summary)) {
  diag_metric <- function(metric) {
    val <- diag_summary %>% filter(.data$metric == !!metric) %>% pull(.data$value)
    if (length(val) == 0) NA_real_ else val[[1]]
  }
  ma("- Clinical primary denominator confirmed by script 32: **%s** rows = UTI **%s**, Not_UTI **%s**.",
     diag_metric("primary_clinical_total"),
     diag_metric("primary_clinical_uti"),
     diag_metric("primary_clinical_not_uti"))
  ma("- VF/model primary denominator confirmed by script 32: **%s** rows = UTI **%s**, Not_UTI **%s**, missing status **%s**.",
     diag_metric("primary_vf_model_total"),
     diag_metric("primary_vf_model_uti"),
     diag_metric("primary_vf_model_not_uti"),
     diag_metric("primary_vf_model_missing_status"))
  ma("- Near-miss rows under the current rule: **%s**; manual primary exclusions: **%s**; quarantined failed/not-expected FASTA rows: **%s**.",
     diag_metric("near_miss_rows"),
     diag_metric("manual_primary_exclusions"),
     diag_metric("quarantined_failed_fastas"))
  ma("- Script 32 writes bootstrap supplementary endpoint effects, Fisher feature screens, leave-one-UTI-out stability, paired resident deltas, transition diagnostics, sparse-count precision context, and duplicate-culture QC metadata.")
  ma("- These diagnostics are descriptive/exploratory. They support interpretation under the current UTI rule; they do not relabel near-miss, duplicate, unknown, or quarantined rows.")
} else {
  ma("- Diagnostic layer unavailable. Run `Rscript 32_uti_not_uti_diagnostic_stats.R` before using final summary diagnostics.")
}
ma("")
ma("## VF Modules And Supplementary Marker/System Endpoints")
if (!is.null(mod_summary)) {
  ma("- VF modules defined: **%d**", nrow(mod_summary))
  ma("- UPEC-candidate modules: **%d**", sum(mod_summary$upec_score_candidate, na.rm = TRUE))
}
if (!is.null(expec_marker_defs)) {
  ma("- ExPEC-like marker framework written to `table_07a_expec_marker_definitions.csv`; kpsM is treated as a proxy marker unless subtype is confirmed.")
}
if (!is.null(expec_marker_summary)) {
  expec_like_rows <- expec_marker_summary %>%
    filter(.data$marker == "expec_like") %>%
    select(Infection_Status, pct_present, n_present, n_episodes)
  if (nrow(expec_like_rows) > 0) {
    ma("- ExPEC-like prevalence by primary status is available in `table_07b_expec_marker_summary_by_status.csv`: %s.",
       paste(sprintf("%s %.1f%% (%d/%d)",
                     expec_like_rows$Infection_Status,
                     expec_like_rows$pct_present,
                     expec_like_rows$n_present,
                     expec_like_rows$n_episodes), collapse = "; "))
  }
}
if (nrow(t07) > 0) {
  uti_scores <- t07 %>% filter(Infection_Status == "UTI") %>% select(score_name, median)
  ma("- Supplementary endpoint summaries by status written to `table_07_vf_score_summary_by_status.csv` (%d rows)", nrow(t07))
  if (nrow(uti_scores) > 0) {
    ma("- UTI endpoint medians available for: %s", paste(uti_scores$score_name, collapse = ", "))
  }
}
if (!is.null(stat_score_tests)) {
  ma("- Targeted participant-collapsed supplementary endpoint sensitivity is available for **%d** prespecified endpoints in `table_17_participant_collapsed_score_tests.csv`.", nrow(stat_score_tests))
}
if (!is.null(stat_glmm)) {
  ma("- Prespecified endpoint-level GLMM sensitivity is available in `table_20_score_glmm_sensitivity.csv`; collapsed ST uses ST131, ST141, Other typed ST, and Missing ST only.")
}
ma("")
ma("## Longitudinal Same-Strain VF Stability")
ma("- Same-strain-first interpretation uses **<=%d SNPs** as the strong same-strain threshold from prior YELLOW study work.", strain_snp_threshold())
ma("- SNP context is separated from ST: **0-%d SNPs = Strong same strain**, **>%d SNPs = Above same-strain SNP threshold**, and missing SNPs remain **Missing SNP evidence**.", strain_snp_threshold(), strain_snp_threshold())
ma("- ST context is secondary lineage information only: Same ST, Different ST, or Missing ST evidence. Same ST does not prove same strain.")
if (nrow(t09) > 0) {
  strong_transitions <- t09 %>%
    filter(snp_strain_context == "Strong same strain") %>%
    summarise(n = sum(n_transitions, na.rm = TRUE), .groups = "drop") %>%
    pull(n)
  ma("- Table 09 now presents SNP context first, then combined pair interpretation, then ST lineage context; SNP-defined strong same-strain transition rows total **%d**.", strong_transitions %||% 0L)
}
if (!is.null(vf_same_strain_by_ST) && nrow(vf_same_strain_by_ST) > 0) {
  ma("- Same-strain-by-ST summary is written as a secondary lineage diagnostic with **%d** ST rows.", nrow(vf_same_strain_by_ST))
}
ma("")
ma("## Not_UTI->UTI Transition Cases")
if (!is.null(not_uti_uti_cases)) {
  ma("- Clinical Not_UTI->UTI transitions indexed: **%d**", nrow(not_uti_uti_cases))
  ma("- WGS/VF-linked Not_UTI->UTI transitions: **%d**", nrow(not_uti_uti_wgs_cases))
  ma("- Not_UTI->UTI transitions missing WGS/VF endpoint: **%d**", nrow(not_uti_uti_cases) - nrow(not_uti_uti_wgs_cases))
}
if (nrow(t10) > 0 && "case_class" %in% names(t10)) {
  stable_n <- sum(str_detect(t10$case_class, "stable VF/module"), na.rm = TRUE)
  ma("- Same-strain stable VF/module Not_UTI->UTI cases: **%d**", stable_n)
}
if (!is.null(stat_transition_module)) {
  ma("- Transition-level module gain/loss enrichment is available in `table_19_transition_module_gain_loss_enrichment.csv`; tests are exploratory because transitions can repeat within residents.")
}
if (!is.null(stat_feature_tests)) {
  ma("- Paired binary feature sensitivity is available in `table_18_paired_binary_feature_sensitivity.csv`; features are restricted to prespecified/top Fisher, GLMM, and module candidates.")
}
ma("")
ma("## Optional VF+AMR/Plasmid")
if (!is.null(vf_amr_combined) && "amr_data_available" %in% names(vf_amr_combined)) {
  if (any(vf_amr_combined$amr_data_available %in% TRUE, na.rm = TRUE)) {
    ma("- True AMR data were integrated for at least one episode.")
  } else {
    ma("- No dedicated AMR database screening output was integrated; script 29 produced VF+plasmid/mobile-context summaries only.")
  }
} else {
  ma("- Optional VF+AMR combined profile table was not available.")
}
ma("")
ma("## QC And Caveats")
if (!is.null(denominator_summary)) {
  ma("- Denominator summary available in `results/qc/pipeline_denominator_summary.csv` (%d rows).", nrow(denominator_summary))
}
ma("- VF figure index and visualisation audit written to `results/vf/vf_figure_index.csv` and `results/vf/vf_visualisation_audit.csv`.")
if (!is.null(uti_attrition)) {
  ma("- UTI attrition audit: clinical UTI **%d**, retained in VF-ready UTI set **%d**.",
     nrow(uti_attrition), sum(uti_attrition$in_vf_analysis_ready %in% TRUE, na.rm = TRUE))
}
if (!is.null(fasta_audit)) {
  ma("- FASTA audit separates canonical candidate FASTAs from legacy/quarantine/cache files; review `results/qc/fasta_usage_audit_summary.txt` before final WGS claims.")
}
if (file.exists(path_core_stale)) {
  core_status <- grep("^Status:", readLines(path_core_stale, warn = FALSE), value = TRUE)
  ma("- Core SNP status: %s", paste(core_status, collapse = "; "))
}
if (file.exists(path_panaroo_stale)) {
  pan_status <- grep("^Status:", readLines(path_panaroo_stale, warn = FALSE), value = TRUE)
  ma("- Panaroo status: %s", paste(pan_status, collapse = "; "))
}
if (file.exists(path_model_warnings)) {
  ma("- Model interpretation warnings are written to `results/models/model_interpretation_warnings.txt`; association claims should be labelled exploratory/underpowered where indicated.")
}
if (!is.null(stat_validation)) {
  ma("- Statistical sensitivity add-on validation written to `results/statistical_sensitivity/statistical_sensitivity_validation_checks.csv`; all new p-value families include BH-adjusted q-values.")
}
if (!is.null(vf_gap)) {
  ma("- VF annotation gap report available in `results/vf/vf_gene_annotation_gap_report.csv`; unassigned genes are separated from curated burden and supplementary marker/system endpoints.")
}
if (is_stale(path_vf_ready, path_vf_pa) && !vf_ready_matches_pa(vf_ready, vf_pa)) {
  ma("- **WARNING:** `vf_analysis_ready.csv` is older than `vf_pa_all.csv` and row/gene content differs; rerun script 22 before trusting final VF extension outputs.")
} else if (is_stale(path_vf_ready, path_vf_pa)) {
  ma("- `vf_analysis_ready.csv` has an older timestamp than `vf_pa_all.csv`, but row count and VF gene set match; this was treated as current.")
}
if (poster_status_unsuitable) {
  ma("- **WARNING:** `status_map_with_poster_tp.csv` is missing, older than `status_map.csv`, or lacks primary UTI columns; rerun script 00d before transition interpretation.")
}
if (length(legacy_found) > 0) {
  ma("- Legacy/stale-looking summary files were detected and are listed in `summary_qc_log.txt`; use current `results/summary/` outputs as truth.")
}
ma("- UTI sample size in the VF-linked dataset is limited; supplementary endpoint tests are descriptive/exploratory and underpowered for definitive UTI association.")
ma("- Repeated measures, ST/lineage structure, remaining missing/ambiguous WGS/VF endpoints, and Uricult bridge assumptions limit causal interpretation.")
ma("- VF presence does not imply expression or functional activity.")
ma("- Do not interpret VFDB-derived or plasmid-only summaries as true AMR analysis.")
ma("- Pangenome/core-SNP and same-strain conclusions require current input hashes; stale or incomplete reports mean those outputs require rerun after GFF/core-SNP refresh.")

writeLines(md, file.path(DIR_SUMMARY, "final_key_results_summary.md"))

# ==============================================================================
# 16. OPTIONAL BUNDLES
# ==============================================================================
all_tables <- list(
  table_01_cohort_episode_flow = t01,
  table_02_clinical_status_counts = t02,
  table_03_wgs_qc_summary = t03,
  table_04_mlst_st_summary = t04,
  table_05_vf_gene_category_summary = t05,
  table_06_vf_module_summary = t06,
  table_07_vf_score_summary_by_status = t07,
  table_07a_expec_marker_definitions = expec_marker_defs %||% tibble(),
  table_07b_expec_marker_summary_by_status = expec_marker_summary %||% tibble(),
  table_07c_expec_marker_tests = expec_marker_tests %||% tibble(),
  table_08_vf_score_summary_by_ST = t08,
  table_09_longitudinal_vf_stability_summary = t09,
  vf_same_strain_vf_stability_summary = vf_same_strain_summary %||% tibble(),
  vf_replacement_vf_change_summary = vf_replacement_summary %||% tibble(),
  vf_strain_context_by_transition_summary = vf_strain_context_summary %||% tibble(),
  vf_same_strain_by_ST_summary = vf_same_strain_by_ST %||% tibble(),
  table_10_not_uti_uti_transition_cases = t10,
  table_11_lineage_context_summary = t11,
  table_12_missing_data_audit = t12,
  table_13_optional_vf_amr_summary = t13,
  table_14_uti_not_uti_diagnostics = t14,
  table_15_uti_not_uti_test_interpretation = diag_interpretation %||% tibble(),
  table_16_uti_not_uti_diagnostic_figures = diag_fig_meta %||% tibble(),
  table_17_participant_collapsed_score_tests = stat_score_tests %||% tibble(),
  table_18_paired_binary_feature_sensitivity = stat_feature_tests %||% tibble(),
  table_19_transition_module_gain_loss_enrichment = stat_transition_module %||% tibble(),
  table_20_score_glmm_sensitivity = stat_glmm %||% tibble(),
  table_21_statistical_sensitivity_validation = stat_validation %||% tibble(),
  stat_participant_collapsed_score_values = stat_score_values %||% tibble(),
  stat_paired_binary_feature_catalog = stat_feature_catalog %||% tibble(),
  stat_paired_binary_feature_deltas = stat_feature_deltas %||% tibble(),
  stat_figure_metadata = stat_fig_meta %||% tibble(),
  diag_uti_not_uti_bootstrap_effects = diag_bootstrap %||% tibble(),
  diag_uti_not_uti_feature_fisher = diag_feature_fisher %||% tibble(),
  diag_uti_not_uti_leave_one_uti_out = diag_leave_one %||% tibble(),
  diag_uti_not_uti_power_precision = diag_power %||% tibble(),
  diag_uti_not_uti_paired_deltas = diag_paired_deltas %||% tibble(),
  diag_uti_not_uti_paired_tests = diag_paired_tests %||% tibble(),
  diag_uti_not_uti_transition_tests = diag_transition_tests %||% tibble(),
  diag_uti_not_uti_attrition = diag_attrition %||% tibble(),
  diag_uti_not_uti_denominator_flow = diag_flow %||% tibble(),
  diag_uti_not_uti_near_miss_rows = diag_near_miss %||% tibble(),
  diag_duplicate_culture_qc_31036 = diag_duplicate_qc %||% tibble(),
  qa_denominator_summary = denominator_summary %||% tibble(),
  qa_metadata_input_manifest = metadata_manifest %||% tibble(),
  qa_metadata_crosswalk = metadata_crosswalk %||% tibble(),
  qa_fasta_audit = fasta_audit %||% tibble(),
  qa_uti_attrition = uti_attrition %||% tibble(),
  qa_panaroo_manifest = panaroo_manifest %||% tibble(),
  qa_qc_selection_bias = qc_bias %||% tibble(),
  qa_vf_annotation_gap = vf_gap %||% tibble(),
  vf_figure_index = figure_index,
  vf_visualisation_audit = visualisation_audit
)
saveRDS(all_tables, file.path(DIR_SUMMARY, "final_summary_tables.rds"))
msg("RDS bundle saved")

tryCatch({
  if (requireNamespace("writexl", quietly = TRUE)) {
    workbook_tables <- all_tables
    writexl::write_xlsx(workbook_tables, file.path(DIR_SUMMARY, "final_summary_tables.xlsx"))
    msg("Excel workbook saved")
  } else {
    msg("writexl not available; Excel export skipped")
  }
}, error = function(e) msg("Excel export skipped: %s", e$message))

msg("✓ 30_vf_project_summary_tables.R complete.")

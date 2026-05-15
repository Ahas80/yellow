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
#   - results/summary/table_10_asb_uti_transition_cases.csv
#   - results/summary/table_11_lineage_context_summary.csv
#   - results/summary/table_12_missing_data_audit.csv
#   - results/summary/table_13_optional_vf_amr_summary.csv
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

standardise_episode_keys <- function(df, tp_col = NULL) {
  if (is.null(df)) return(NULL)
  out <- df %>% mutate(Participant_id = as.character(Participant_id))
  if (!is.null(tp_col) && tp_col %in% names(out)) {
    out <- out %>% mutate(tp_lab = normalise_tp_label(.data[[tp_col]]))
  } else if ("tp_lab" %in% names(out)) {
    out <- out %>% mutate(tp_lab = normalise_tp_label(tp_lab))
  }
  out
}

status_counts <- function(df, source_file = NA_character_) {
  if (is.null(df) || !"Infection_Status" %in% names(df)) {
    return(tibble(n_ASB = NA_integer_, n_UTI = NA_integer_,
                  n_Negative = NA_integer_, n_other_status = NA_integer_))
  }
  tibble(
    n_ASB = sum(df$Infection_Status == "ASB", na.rm = TRUE),
    n_UTI = sum(df$Infection_Status == "UTI", na.rm = TRUE),
    n_Negative = sum(df$Infection_Status == "Negative", na.rm = TRUE),
    n_other_status = sum(!df$Infection_Status %in% c("ASB", "UTI", "Negative"), na.rm = TRUE)
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
path_gene_map <- file.path(DIR_VF, "gene_map.csv")
path_qc_wgs <- file.path(DIR_WGS, "qc_summary.csv")
path_mlst_freq <- file.path(DIR_MLST, "ST_frequencies.csv")
path_mlst_meta <- FILE_MLST_CANONICAL
path_mod_summary <- file.path(DIR_VF, "vf_module_summary.csv")
path_mod_map <- file.path(DIR_VF, "gene_module_map.csv")
path_scores <- file.path(DIR_VF, "vf_score_table.csv")
path_scores_status <- file.path(DIR_VF, "vf_score_summary_by_status.csv")
path_scores_st <- file.path(DIR_VF, "vf_score_summary_by_ST.csv")
path_vf_trans <- file.path(DIR_VF, "vf_longitudinal_transitions.csv")
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

status_map <- standardise_episode_keys(safe_read(path_status), "Timepoint")
status_poster <- standardise_episode_keys(safe_read(path_status_poster), "Timepoint")
vf_pa <- standardise_episode_keys(safe_read(path_vf_pa))
vf_ready <- standardise_episode_keys(safe_read(path_vf_ready))
gene_map <- safe_read(path_gene_map)
qc_wgs <- safe_read(path_qc_wgs)
mlst_freq <- safe_read(path_mlst_freq)
mlst_meta <- standardise_episode_keys(safe_read(path_mlst_meta), "Timepoint")
mod_summary <- safe_read(path_mod_summary)
mod_map <- safe_read(path_mod_map)
score_table <- standardise_episode_keys(safe_read(path_scores))
scores_status <- safe_read(path_scores_status)
scores_st <- safe_read(path_scores_st)
vf_trans <- standardise_episode_keys(safe_read(path_vf_trans))
case_index <- standardise_episode_keys(safe_read(path_case_index))
case_summary <- standardise_episode_keys(safe_read(path_case_summary))
strain_ctx <- safe_read(path_strain_ctx)
vf_amr_combined <- standardise_episode_keys(safe_read(path_vf_amr_combined))
vf_plasmid <- standardise_episode_keys(safe_read(path_vf_plasmid))
vf_amr_status <- safe_read(path_vf_amr_status)
denominator_summary <- safe_read(path_denominator)
metadata_manifest <- safe_read(path_metadata_manifest)
metadata_crosswalk <- safe_read(path_metadata_crosswalk)
fasta_audit <- safe_read(path_fasta_audit)
uti_attrition <- safe_read(path_uti_attrition)
panaroo_manifest <- safe_read(path_panaroo_manifest)
qc_bias <- safe_read(path_qc_bias)
vf_gap <- safe_read(path_vf_gap)
model_denom <- standardise_episode_keys(safe_read(path_model_denom))

if (is.null(status_map)) stop("Missing required clinical status map: ", path_status)
if (is.null(vf_pa)) stop("Missing required raw VF P/A matrix: ", path_vf_pa)
if (is.null(vf_ready)) stop("Missing required canonical VF-ready table: ", path_vf_ready)

# ==============================================================================
# 3. TABLE 01: COHORT AND EPISODE FLOW
# ==============================================================================
asb_uti_model <- vf_ready %>% filter(Infection_Status %in% c("ASB", "UTI"))
longitudinal_subset <- if (!is.null(vf_trans)) {
  vf_ready %>% semi_join(vf_trans %>% select(Participant_id), by = "Participant_id")
} else NULL
asb_uti_cases <- if (!is.null(case_index)) case_index %>% filter(is_asb_to_uti %in% TRUE) else NULL
asb_uti_wgs_cases <- if (!is.null(asb_uti_cases)) asb_uti_cases %>% filter(has_vf_pair %in% TRUE) else NULL
vf_amr_complete <- if (!is.null(vf_amr_combined)) {
  vf_amr_combined %>% filter(amr_data_available %in% TRUE)
} else NULL

t01 <- bind_rows(
  flow_row(1, "All-batch clinical episodes", status_map, path_status, "Canonical clinical denominator"),
  flow_row(2, "Ordered/poster clinical episodes", status_poster, path_status_poster,
           "Fresh ordered file used only for timelines when available"),
  flow_row(3, "Raw VF presence/absence rows", vf_pa, path_vf_pa, "Sequenced/VFDB-called episode rows"),
  flow_row(4, "Canonical VF-ready episodes", vf_ready, path_vf_ready,
           "Raw VF rows joined to clinical status/ST by script 22"),
  flow_row(5, "ASB/UTI VF modelling subset", asb_uti_model, path_vf_ready,
           "Negative and unknown statuses excluded"),
  flow_row(6, "Repeated-measures VF longitudinal subset", longitudinal_subset, path_vf_trans,
           "Participants represented in VF transition output"),
  flow_row(7, "Clinical ASB->UTI transition cases", asb_uti_cases, path_case_index,
           "All ordered clinical ASB->UTI transitions indexed by script 28"),
  flow_row(8, "WGS/VF-linked ASB->UTI transition cases", asb_uti_wgs_cases, path_case_index,
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
  status_layer(asb_uti_model, "asb_uti_vf_subset", path_vf_ready)
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
    filter(!is.na(ST), ST != "-") %>%
    group_by(ST) %>%
    summarise(
      n_episodes_or_isolates = n(),
      n_participants = n_distinct(Participant_id),
      pct = round(100 * n() / nrow(vf_ready), 1),
      n_ASB = sum(Infection_Status == "ASB", na.rm = TRUE),
      n_UTI = sum(Infection_Status == "UTI", na.rm = TRUE),
      n_Negative = sum(Infection_Status == "Negative", na.rm = TRUE),
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
      n_ASB_present = sum(vf_ready[[Gene]][vf_ready$Infection_Status == "ASB"], na.rm = TRUE),
      n_UTI_present = sum(vf_ready[[Gene]][vf_ready$Infection_Status == "UTI"], na.rm = TRUE),
      n_Negative_present = sum(vf_ready[[Gene]][vf_ready$Infection_Status == "Negative"], na.rm = TRUE),
      assignment_status = ifelse(Category == "Unassigned", "unassigned", "assigned")
    ) %>%
    ungroup() %>%
    select(summary_level, Gene, Category, n_episodes_present, pct_episodes_present,
           n_ASB_present, n_UTI_present, n_Negative_present, assignment_status)

  category_summary <- gene_summary %>%
    group_by(Category) %>%
    summarise(
      summary_level = "category",
      Gene = NA_character_,
      n_episodes_present = NA_integer_,
      pct_episodes_present = NA_real_,
      n_ASB_present = sum(n_ASB_present > 0, na.rm = TRUE),
      n_UTI_present = sum(n_UTI_present > 0, na.rm = TRUE),
      n_Negative_present = sum(n_Negative_present > 0, na.rm = TRUE),
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

score_cols <- c("total_vf_count_all", "total_vf_count_curated",
                "total_vf_count_upec_candidate", "total_vf_count_unassigned",
                "low_confidence_count", "vf_count_total", "vf_count_curated", "vf_count_unassigned",
                "upec_gene_score_unweighted", "upec_gene_score",
                "module_count_present", "n_modules_present",
                "upec_module_score_unweighted", "n_upec_modules_present",
                "upec_score_fraction")
t07 <- if (!is.null(scores_status)) scores_status else score_summary_from_table(score_table, "Infection_Status", score_cols)
t08 <- if (!is.null(scores_st)) scores_st else score_summary_from_table(score_table, "ST", score_cols)
write_csv(t07, file.path(DIR_SUMMARY, "table_07_vf_score_summary_by_status.csv"))
write_csv(t08, file.path(DIR_SUMMARY, "table_08_vf_score_summary_by_ST.csv"))

# ==============================================================================
# 9. TABLE 09: LONGITUDINAL VF STABILITY SUMMARY
# ==============================================================================
if (!is.null(vf_trans)) {
  t09 <- vf_trans %>%
    group_by(transition_type) %>%
    summarise(
      n_transitions = n(),
      n_participants = n_distinct(Participant_id),
      median_vf_jaccard = median(jaccard_similarity, na.rm = TRUE),
      q25_vf_jaccard = quantile(jaccard_similarity, 0.25, na.rm = TRUE),
      q75_vf_jaccard = quantile(jaccard_similarity, 0.75, na.rm = TRUE),
      n_zero_gene_change = sum(n_gained == 0 & n_lost == 0, na.rm = TRUE),
      n_any_gene_change = sum(n_gained > 0 | n_lost > 0, na.rm = TRUE),
      median_genes_gained = median(n_gained, na.rm = TRUE),
      median_genes_lost = median(n_lost, na.rm = TRUE),
      source_file = path_vf_trans,
      .groups = "drop"
    ) %>%
    arrange(desc(n_transitions))
} else {
  t09 <- tibble(transition_type = character(), source_file = path_vf_trans)
}
write_csv(t09, file.path(DIR_SUMMARY, "table_09_longitudinal_vf_stability_summary.csv"))
write_csv(t09, file.path(DIR_SUMMARY, "table_09_longitudinal_vf_stability.csv")) # backward-compatible alias

# ==============================================================================
# 10. TABLE 10: ASB->UTI TRANSITION CASES
# ==============================================================================
if (!is.null(case_summary)) {
  t10 <- case_summary %>%
    filter(is_asb_to_uti %in% TRUE | str_detect(transition_type, "ASB.*UTI")) %>%
    select(any_of(c("case_id", "Participant_id", "from_tp", "to_tp", "from_status", "to_status",
                    "transition_type", "has_vf_pair", "has_module_pair", "has_score_pair",
                    "is_uricult_transition", "timing_caveat", "ST_from", "ST_to",
                    "same_ST", "SNPs", "same_strain_evidence", "vf_jaccard",
                    "module_jaccard", "n_vf_genes_gained", "n_vf_genes_lost",
                    "n_modules_gained", "n_modules_lost", "delta_upec_module_score",
                    "case_class", "missing_data_note", "interpretation_short")))
} else if (!is.null(case_index)) {
  t10 <- case_index %>% filter(is_asb_to_uti %in% TRUE)
} else {
  t10 <- tibble(case_id = character())
}
write_csv(t10, file.path(DIR_SUMMARY, "table_10_asb_uti_transition_cases.csv"))

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
      n_ASB = sum(Infection_Status == "ASB", na.rm = TRUE),
      n_UTI = sum(Infection_Status == "UTI", na.rm = TRUE),
      median_vf_burden = median(vf_count_total, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    filter(n_episodes >= 5) %>%
    arrange(desc(n_episodes)) %>%
    mutate(
      question = "Are VF burden and clinical status distributed across common STs?",
      test_or_summary = "Descriptive common-ST summary",
      result_value = sprintf("median VF burden %.1f; ASB=%d; UTI=%d", median_vf_burden, n_ASB, n_UTI),
      p_value = NA_real_,
      interpretation = "Lineage context should be considered before interpreting VF-status differences",
      source_file = path_vf_ready
    ) %>%
    select(question, ST, test_or_summary, n_episodes, n_participants,
           n_ASB, n_UTI, result_value, p_value, interpretation, source_file)
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
    filter(is_asb_to_uti %in% TRUE, !(has_vf_pair %in% TRUE)) %>%
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
  missing_row("Raw VF rows without clinical status", raw_vf_no_status, paste(path_vf_pa, path_status, sep = "; "),
              "Clinical interpretation unavailable for these VF rows", "Unmatched participant/timepoint"),
  missing_row("VF-ready episodes missing ST", vf_missing_st, path_vf_ready,
              "Lineage-aware summaries incomplete", "Missing/non-typable MLST or join issue"),
  missing_row("VF-ready episodes missing Infection_Status", vf_missing_status, path_vf_ready,
              "Excluded from status-stratified VF analyses", "Clinical join missing"),
  missing_row("ASB->UTI transitions missing WGS/VF endpoint", transition_missing_wgs, path_case_index,
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
    missingness_category = "status_map_with_poster_tp freshness",
    n_episodes = if (!is.null(status_poster)) nrow(status_poster) else NA_integer_,
    n_participants = if (!is.null(status_poster)) n_distinct(status_poster$Participant_id) else NA_integer_,
    n_ASB = status_counts(status_poster)$n_ASB,
    n_UTI = status_counts(status_poster)$n_UTI,
    n_Negative = status_counts(status_poster)$n_Negative,
    n_other_status = status_counts(status_poster)$n_other_status,
    likely_reason = if (file.exists(path_status_poster) && !is_stale(path_status_poster, path_status)) "Fresh" else "Missing or older than status_map.csv",
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
# 14. VF FIGURE INDEX AND VISUALISATION AUDIT
# ==============================================================================
vf_denominator_label <- sprintf(
  "VF-ready isolates n=%d; participants n=%d; ASB n=%d; UTI n=%d; Negative n=%d; missing/other status n=%d",
  nrow(vf_ready),
  n_distinct(vf_ready$Participant_id),
  status_counts(vf_ready)$n_ASB,
  status_counts(vf_ready)$n_UTI,
  status_counts(vf_ready)$n_Negative,
  sum(is.na(vf_ready$Infection_Status) | !vf_ready$Infection_Status %in% c("ASB", "UTI", "Negative"), na.rm = TRUE)
)

clinical_denominator_label <- sprintf(
  "Clinical episodes n=%d; participants n=%d; ASB n=%d; UTI n=%d; Negative n=%d",
  nrow(status_map),
  n_distinct(status_map$Participant_id),
  status_counts(status_map)$n_ASB,
  status_counts(status_map)$n_UTI,
  status_counts(status_map)$n_Negative
)

transition_denominator_label <- if (!is.null(vf_trans)) {
  sprintf("Consecutive within-resident VF comparisons n=%d; participants n=%d",
          nrow(vf_trans), n_distinct(vf_trans$Participant_id))
} else {
  "Consecutive within-resident VF comparisons unavailable"
}

model_denominator_label <- if (!is.null(model_denom)) {
  sprintf("ASB/UTI model rows n=%d; participants n=%d; ASB n=%d; UTI n=%d",
          nrow(model_denom), n_distinct(model_denom$Participant_id),
          status_counts(model_denom)$n_ASB, status_counts(model_denom)$n_UTI)
} else {
  "ASB/UTI model denominator unavailable; see results/models/model_dataset_denominator.csv"
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
  "Burden", "23_vf_cross_sectional.R", "vf_burden_by_status", "plots/vf/vf_burden_by_status.png", "Distribution of E. coli virulence factor burden by clinical status", path_vf_ready, vf_denominator_label, "isolate-level", "None shown", "Descriptive", "Repeated isolates, small UTI denominator, and ST/lineage/event-type confounding limit causal interpretation.",
  "Burden", "23_vf_cross_sectional.R", "vf_burden_participant_summary", "plots/vf/vf_burden_participant_summary.png", "Participant-level summary of E. coli VF burden by clinical status", path_vf_ready, vf_denominator_label, "participant-status summary", "None shown", "Descriptive", "Participant summaries reduce repeated-isolate weighting but do not model confounding.",
  "Burden", "23_vf_cross_sectional.R", "vf_burden_paired_asb_uti", "plots/vf/vf_burden_paired_asb_uti.png", "Paired participant-level VF burden in residents with ASB and UTI isolates", path_vf_ready, vf_denominator_label, "paired participant-level ASB-UTI comparison", "None shown", "Descriptive/exploratory", "Restricted to residents with both ASB and UTI VF-ready isolates; UTI n is small.",
  "Gene prevalence", "23_vf_cross_sectional.R", "vf_top_gene_prevalence", "plots/vf/vf_top_gene_prevalence.png", "Most prevalent virulence factor genes among VF/WGS-linked E. coli isolates", path_vf_ready, vf_denominator_label, "isolate-level gene prevalence", "None shown", "Descriptive", "Prevalence ranking is not an association test and does not imply expression or function.",
  "Gene prevalence", "23_vf_cross_sectional.R", "vf_gene_prevalence_difference_asb_uti", "plots/vf/vf_gene_prevalence_difference_asb_uti.png", "Virulence factor genes with largest ASB-UTI prevalence differences", path_vf_ready, vf_denominator_label, "isolate-level ASB-UTI gene prevalence difference", "Fisher tests in companion CSV; BH q-values", "Exploratory", "Repeated isolates are not modelled; rare genes give unstable contrasts; lineage may confound differences.",
  "Gene prevalence", "23_vf_cross_sectional.R", "vf_gene_prevalence_heatmap", "plots/vf/vf_gene_prevalence_heatmap.png", "VF gene prevalence heatmap across clinical states", path_vf_ready, vf_denominator_label, "isolate-level gene prevalence heatmap", "None shown", "Descriptive", "Percentages are denominator-sensitive and should be interpreted with status counts.",
  "Model evidence bridge", "14_genotype_phenotype_model.R", "vf_gene_screening_vs_model_evidence", "plots/vf/vf_gene_screening_vs_model_evidence.png", "Exploratory VF gene screening versus participant-aware model evidence", paste(path_model_univ, path_model_glmm, sep = "; "), model_denominator_label, "gene/feature-level association screen", "Fisher screening plus GLMM/GLM model OR, CI, and FDR", "Exploratory/model diagnostic", "Nominal Fisher hits do not imply robust corrected or participant-aware model support; sparse fits and lineage/event structure limit inference.",
  "Category profiles", "23_vf_cross_sectional.R", "vf_category_burden_by_status", "plots/vf/vf_category_burden_by_status.png", "Virulence factor category burden across clinical states", path_vf_ready, vf_denominator_label, "isolate-level category burden", "Category Fisher/Wilcoxon tests in CSV; BH q-values", "Exploratory", "Categories are descriptive, not validated causal virulence scores.",
  "Category profiles", "26_vf_define_gene_modules.R", "module_gene_counts", "plots/vf/module_gene_counts.png", "Virulence factor genes per curated biological module", paste(path_gene_map, path_vf_ready, sep = "; "), vf_denominator_label, "gene-to-module curation", "None", "Descriptive", "Module definitions are analysis curation units, not validated UTI scores.",
  "Category profiles", "26_vf_define_gene_modules.R", "module_prevalence_by_status", "plots/vf/module_prevalence_by_status.png", "Virulence factor module prevalence across clinical states", path_vf_ready, vf_denominator_label, "isolate-level module presence", "None shown", "Descriptive/exploratory", "Repeated measures and lineage structure are not modelled.",
  "Category profiles", "26_vf_define_gene_modules.R", "vf_category_composition_by_status", "plots/vf/vf_category_composition_by_status.png", "Virulence factor category profiles across clinical states", path_vf_ready, vf_denominator_label, "isolate-level category composition", "None shown", "Descriptive", "Category composition does not imply causal virulence.",
  "Annotation confidence", "26_vf_define_gene_modules.R", "vf_module_assignment_confidence", "plots/vf/vf_module_assignment_confidence.png", "VF module assignment confidence and gene-map coverage", paste(path_mod_assign, path_vf_gap, sep = "; "), "Gene/module assignment audit rows; see vf_gene_annotation_gap_report.csv for current VF-matrix gap counts", "gene-to-module curation diagnostic", "None shown", "Diagnostic/descriptive", "Low-confidence, unassigned, and gene-map-absent VFDB hits should not be interpreted as validated virulence modules.",
  "Longitudinal stability", "24_vf_longitudinal_dynamics.R", "vf_jaccard_by_transition", "plots/vf/vf_jaccard_by_transition.png", "Within-resident virulence factor similarity across repeated E. coli isolates", path_vf_ready, transition_denominator_label, "consecutive within-resident isolate-pair comparison", "None shown", "Descriptive", "Intervals vary; transition-specific UTI comparisons are sparse.",
  "Longitudinal stability", "24_vf_longitudinal_dynamics.R", "vf_within_host_jaccard_distribution", "plots/vf/vf_within_host_jaccard_distribution.png", "Distribution of within-resident VF profile similarity", path_vf_ready, transition_denominator_label, "consecutive within-resident isolate-pair comparison", "None shown", "Descriptive", "High similarity supports stability but does not prove same-strain persistence without SNP/ANI context.",
  "Longitudinal stability", "24_vf_longitudinal_dynamics.R", "vf_jaccard_by_days_between_samples", "plots/vf/vf_jaccard_by_days_between_samples.png", "VF similarity versus time between repeated isolates", path_vf_ready, transition_denominator_label, "dated consecutive within-resident isolate-pair comparison", "None shown", "Descriptive", "Only pairs with parseable dates are shown.",
  "Longitudinal stability", "24_vf_longitudinal_dynamics.R", "vf_jaccard_same_vs_different_st", "plots/vf/vf_jaccard_same_vs_different_st.png", "Within-resident VF similarity by sequence-type consistency", path_vf_ready, transition_denominator_label, "consecutive pair stratified by ST agreement", "None shown", "Diagnostic/descriptive", "ST agreement is a lineage diagnostic and does not alone prove same strain.",
  "Longitudinal stability", "24_vf_longitudinal_dynamics.R", "vf_gene_gain_loss_consecutive_pairs", "plots/vf/vf_gene_gain_loss_consecutive_pairs.png", "VF gene gains and losses between repeated E. coli isolates", path_vf_ready, transition_denominator_label, "consecutive within-resident gene gain/loss summary", "None shown", "Descriptive", "Gain/loss can reflect replacement, assembly/calling differences, or true gene-content change.",
  "Transition case studies", "28_vf_transition_case_studies.R", "vf_transition_case_timeline", "plots/vf/vf_transition_case_timeline.png", "Clinical timelines for residents with ASB-to-UTI transitions", paste(path_status, path_status_poster, path_vf_ready, sep = "; "), clinical_denominator_label, "clinical episode timeline", "None shown", "Descriptive", "Uricult ordering uses dates where available; fallback/display ordering must be cited.",
  "Transition case studies", "28_vf_transition_case_studies.R", "vf_transition_score_slopeplot", "plots/vf/vf_transition_score_slopeplot.png", "Virulence factor score changes across ASB-to-UTI transitions", path_vf_ready, "WGS/VF-linked ASB-to-UTI transition cases", "transition case-study pair", "None shown", "Descriptive", "Changes may reflect lineage replacement or technical differences and do not imply causality.",
  "Transition case studies", "28_vf_transition_case_studies.R", "vf_transition_module_change_heatmap", "plots/vf/vf_transition_module_change_heatmap.png", "VF module changes in ASB-to-UTI transition cases", path_vf_ready, "WGS/VF-linked ASB-to-UTI transition cases", "transition case-study pair", "None shown", "Descriptive", "Modules are descriptive and should be interpreted with ST/SNP evidence.",
  "Transition case studies", "28_vf_transition_case_studies.R", "vf_transition_gene_gain_loss_tile", "plots/vf/vf_transition_gene_gain_loss_tile.png", "VF gene gains and losses in ASB-to-UTI transition cases", path_vf_ready, "WGS/VF-linked ASB-to-UTI transition cases", "transition case-study pair", "None shown", "Descriptive", "Gene gain/loss should be interpreted with lineage and strain context.",
  "Transition case studies", "28_vf_transition_case_studies.R", "vf_transition_snp_vs_vf_jaccard", "plots/vf/vf_transition_snp_vs_vf_jaccard.png", "SNP distance versus VF similarity in ASB-to-UTI transition cases", path_vf_ready, "ASB-to-UTI cases with WGS/VF, SNP distance, and VF Jaccard", "transition case-study pair", "None shown", "Descriptive/diagnostic", "Low SNP distance plus high VF similarity supports persistence but does not establish VF causality for UTI.",
  "Transition case studies", "28_vf_transition_case_studies.R", "vf_asb_uti_transition_strain_context", "plots/vf/vf_asb_uti_transition_strain_context.png", "ASB-to-UTI VF changes require strain-context interpretation", paste(path_case_summary, path_strain_ctx, sep = "; "), "ASB-to-UTI cases with WGS/VF endpoints, SNP distance, and VF burden change", "transition case-study strain-context diagnostic", "None shown", "Descriptive/diagnostic", "VF burden changes across ASB-to-UTI cases may reflect strain replacement, technical differences, or true gene-content change; they do not establish UTI causality.",
  "Transition case studies", "28_vf_transition_case_studies.R", "vf_asb_uti_transition_case_classes", "plots/vf/vf_asb_uti_transition_case_classes.png", "ASB-to-UTI transition case classes and missing genomic endpoints", path_case_summary, "All clinical ASB-to-UTI transition cases, including missing WGS/VF endpoints", "transition case classification", "None shown", "Descriptive/diagnostic", "Missing WGS/VF endpoints remain part of the clinical transition denominator and should not be hidden.",
  "Confounding checks", "25_vf_lineage_vf_interaction.R", "vf_burden_by_st", "plots/vf/vf_burden_by_st.png", "Virulence factor burden varies by E. coli sequence type", path_vf_ready, vf_denominator_label, "isolate-level ST diagnostic", "Kruskal-Wallis in summary text", "Exploratory/diagnostic", "Repeated isolates are not modelled; sparse STs are filtered.",
  "Confounding checks", "25_vf_lineage_vf_interaction.R", "vf_burden_st_x_status", "plots/vf/vf_burden_st_x_status.png", "ASB-UTI VF burden contrasts within E. coli sequence types", path_vf_ready, vf_denominator_label, "isolate-level within-ST diagnostic", "Wilcoxon in summary text when possible", "Exploratory/diagnostic", "Within-ST UTI counts are often too small for stable inference.",
  "Confounding checks", "25_vf_lineage_vf_interaction.R", "vf_st_composition_by_status", "plots/vf/vf_st_composition_by_status.png", "Clinical status distribution across E. coli sequence types", path_vf_ready, vf_denominator_label, "isolate-level ST composition diagnostic", "Fisher simulated in summary text", "Exploratory/diagnostic", "ST distribution differences can confound naive VF-status associations.",
  "Confounding checks", "25_vf_lineage_vf_interaction.R", "vf_batch_by_status", "plots/vf/vf_batch_by_status.png", "Batch structure across VF-ready clinical states", path_vf_ready, vf_denominator_label, "isolate-level batch diagnostic", "None shown", "Diagnostic", "Batch imbalance can affect interpretation of associations.",
  "Confounding checks", "25_vf_lineage_vf_interaction.R", "vf_event_type_by_status", "plots/vf/vf_event_type_by_status.png", "Sampling context across VF-ready clinical states", path_vf_ready, vf_denominator_label, "isolate-level event-type diagnostic", "None shown", "Diagnostic", "Routine and event-driven samples are not interchangeable denominators.",
  "Confounding checks", "25_vf_lineage_vf_interaction.R", "vf_status_timepoint_event_tile", "plots/vf/vf_status_timepoint_event_tile.png", "Clinical status, timepoint, and event context in VF-ready isolates", path_vf_ready, vf_denominator_label, "isolate-level timepoint and event-context diagnostic", "None shown", "Diagnostic", "UTI-labelled VF/WGS rows are structurally concentrated in event-driven labels and should not be read as routine numeric timepoints.",
  "Confounding checks", "25_vf_lineage_vf_interaction.R", "vf_qc_selection_by_status", "plots/vf/vf_qc_selection_by_status.png", "WGS/QC selection structure by clinical status", path_qc_bias, "QC/selection rows by clinical status from qc_selection_bias_by_status.csv", "QC and selection-bias diagnostic", "None shown", "Diagnostic", "Unequal QC or selection inclusion across clinical status can alter VF-ready denominators.",
  "Confounding checks", "25_vf_lineage_vf_interaction.R", "vf_denominator_flow", "plots/vf/vf_denominator_flow.png", "Clinical-to-genomic denominator flow for VF analysis", paste(path_status, path_vf_pa, path_vf_ready, sep = "; "), paste(clinical_denominator_label, vf_denominator_label, sep = " -> "), "denominator-flow diagnostic", "None", "Diagnostic", "Do not hide attrition from clinical episodes to VF-ready analyses.",
  "Confounding checks", "25_vf_lineage_vf_interaction.R", "vf_uricult_join_diagnostic", "plots/vf/vf_uricult_join_diagnostic.png", "Uricult clinical events require harmonisation with UTI-labelled WGS isolates", paste(path_status, path_vf_ready, sep = "; "), "Clinical UTI versus VF-ready UTI bridge counts", "join/denominator diagnostic", "None", "Diagnostic", "Uricult clinical labels and UTI-N WGS labels require explicit bridge assumptions.",
  "Score framework", "27_vf_score_framework.R", "vf_scores_by_status", "plots/vf/vf_scores_by_status.png", "Exploratory VF score distributions by clinical status", path_vf_ready, vf_denominator_label, "isolate-level score summary", "Wilcoxon in results/vf/vf_score_tests_exploratory.csv", "Exploratory", "Scores are not validated predictors; repeated isolates and lineage are not modelled.",
  "Score framework", "27_vf_score_framework.R", "vf_score_effect_summary_asb_uti", "plots/vf/vf_score_effect_summary_asb_uti.png", "Exploratory ASB-UTI VF score differences", path_score_tests, vf_denominator_label, "isolate-level score effect summary", "Wilcoxon with BH q-values", "Exploratory", "Median score differences and BH q-values are descriptive because participant clustering, ST/lineage, timepoint, batch, and event type are not modelled here.",
  "Score framework", "27_vf_score_framework.R", "vf_scores_by_ST", "plots/vf/vf_scores_by_ST.png", "Exploratory VF score distributions by E. coli sequence type", path_vf_ready, vf_denominator_label, "isolate-level lineage diagnostic", "None shown", "Diagnostic/descriptive", "VF scores are expected to be lineage structured.",
  "Score framework", "27_vf_score_framework.R", "vf_score_correlation_heatmap", "plots/vf/vf_score_correlation_heatmap.png", "Spearman correlations among exploratory VF scores", path_vf_ready, vf_denominator_label, "score-score correlation", "Spearman correlation", "Descriptive", "Correlations reflect shared score components and do not indicate clinical prediction.",
  "Score framework", "27_vf_score_framework.R", "vf_pca_status", "plots/vf/vf_pca_status.png", "Exploratory PCA of VF module profiles by clinical status", path_vf_ready, vf_denominator_label, "ordination", "PCA", "Exploratory/descriptive", "Ordination is not adjusted for repeated residents or lineage.",
  "Score framework", "27_vf_score_framework.R", "vf_pca_ST", "plots/vf/vf_pca_ST.png", "Exploratory PCA of VF module profiles by sequence type", path_vf_ready, vf_denominator_label, "ordination", "PCA", "Diagnostic/descriptive", "Lineage clustering should be considered as confounding.",
  "Score framework", "27_vf_score_framework.R", "vf_pcoa_jaccard_status", "plots/vf/vf_pcoa_jaccard_status.png", "Exploratory Jaccard PCoA of VF module profiles by clinical status", path_vf_ready, vf_denominator_label, "ordination", "Jaccard PCoA", "Exploratory/descriptive", "Ordination is descriptive and denominator-sensitive.",
  "Score framework", "27_vf_score_framework.R", "vf_pcoa_jaccard_ST", "plots/vf/vf_pcoa_jaccard_ST.png", "Exploratory Jaccard PCoA of VF module profiles by sequence type", path_vf_ready, vf_denominator_label, "ordination", "Jaccard PCoA", "Diagnostic/descriptive", "Lineage clustering should be considered as confounding.",
  "Optional VF+plasmid", "29_vf_amr_combined_profile.R", "vf_plasmid_analysis_scope", "plots/vf_amr/vf_plasmid_analysis_scope.png", "Scope of VF, plasmid, and AMR data integration", paste(path_vf_amr_report, path_vf_amr_combined, sep = "; "), amr_scope_denominator_label, "input availability and analysis-scope diagnostic", "None shown", "Diagnostic", "Script 29 is VF+plasmid/mobile-context only unless dedicated AMR screening rows are present; plasmid summaries are not true AMR analysis."
) %>%
  mutate(
    output_exists = file.exists(file.path(DIR_ROOT, output_file)),
    note = ifelse(output_exists, "Generated or expected from current pipeline", "Optional/skipped if required inputs were absent")
  )

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
      figure_id %in% c("vf_asb_uti_transition_strain_context", "vf_asb_uti_transition_case_classes") ~
        "ADD companion plot / MARK diagnostic",
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
qa("status_map_with_poster_tp stale: %s", is_stale(path_status_poster, path_status))
qa("vf_analysis_ready stale by timestamp: %s", is_stale(path_vf_ready, path_vf_pa))
qa("vf_analysis_ready row/gene content matches vf_pa_all: %s", vf_ready_matches_pa(vf_ready, vf_pa))
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
ma("- Clinical status counts: ASB **%d**, UTI **%d**, Negative **%d**, Other/unknown **%d**",
   status_counts(status_map)$n_ASB, status_counts(status_map)$n_UTI,
   status_counts(status_map)$n_Negative, status_counts(status_map)$n_other_status)
if (!is.null(status_poster)) {
  ma("- Ordered/poster clinical episodes: **%d** (%d participants)", nrow(status_poster), n_distinct(status_poster$Participant_id))
}
ma("")
ma("## VF/WGS-Ready Dataset")
ma("- Raw VF matrix rows: **%d** (%d participants)", nrow(vf_pa), n_distinct(vf_pa$Participant_id))
ma("- Raw VF gene columns: **%d**", length(raw_vf_gene_cols))
ma("- Canonical VF-ready episodes: **%d** (%d participants)", nrow(vf_ready), n_distinct(vf_ready$Participant_id))
ma("- VF-ready status counts: ASB **%d**, UTI **%d**, Negative **%d**",
   status_counts(vf_ready)$n_ASB, status_counts(vf_ready)$n_UTI, status_counts(vf_ready)$n_Negative)
ma("- VF-ready gene columns: **%d**", length(vf_gene_cols))
if ("ST" %in% names(vf_ready)) {
  ma("- Distinct STs in VF-ready dataset: **%d**", n_distinct(vf_ready$ST[!is.na(vf_ready$ST) & vf_ready$ST != "-"]))
}
ma("- Median total VF burden: **%.0f** (IQR %.0f-%.0f)",
   median(vf_ready$vf_count_total, na.rm = TRUE),
   quantile(vf_ready$vf_count_total, 0.25, na.rm = TRUE),
   quantile(vf_ready$vf_count_total, 0.75, na.rm = TRUE))
ma("")
ma("## VF Modules And Scores")
if (!is.null(mod_summary)) {
  ma("- VF modules defined: **%d**", nrow(mod_summary))
  ma("- UPEC-candidate modules: **%d**", sum(mod_summary$upec_score_candidate, na.rm = TRUE))
}
if (nrow(t07) > 0) {
  uti_scores <- t07 %>% filter(Infection_Status == "UTI") %>% select(score_name, median)
  ma("- Score summaries by status written to `table_07_vf_score_summary_by_status.csv` (%d rows)", nrow(t07))
  if (nrow(uti_scores) > 0) {
    ma("- UTI score medians available for: %s", paste(uti_scores$score_name, collapse = ", "))
  }
}
ma("")
ma("## ASB->UTI Transition Cases")
if (!is.null(asb_uti_cases)) {
  ma("- Clinical ASB->UTI transitions indexed: **%d**", nrow(asb_uti_cases))
  ma("- WGS/VF-linked ASB->UTI transitions: **%d**", nrow(asb_uti_wgs_cases))
  ma("- ASB->UTI transitions missing WGS/VF endpoint: **%d**", nrow(asb_uti_cases) - nrow(asb_uti_wgs_cases))
}
if (nrow(t10) > 0 && "case_class" %in% names(t10)) {
  stable_n <- sum(str_detect(t10$case_class, "stable VF/module"), na.rm = TRUE)
  ma("- Same-strain stable VF/module ASB->UTI cases: **%d**", stable_n)
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
if (!is.null(vf_gap)) {
  ma("- VF annotation gap report available in `results/vf/vf_gene_annotation_gap_report.csv`; unassigned genes are separated from curated/UPEC-candidate counts.")
}
if (is_stale(path_vf_ready, path_vf_pa) && !vf_ready_matches_pa(vf_ready, vf_pa)) {
  ma("- **WARNING:** `vf_analysis_ready.csv` is older than `vf_pa_all.csv` and row/gene content differs; rerun script 22 before trusting final VF extension outputs.")
} else if (is_stale(path_vf_ready, path_vf_pa)) {
  ma("- `vf_analysis_ready.csv` has an older timestamp than `vf_pa_all.csv`, but row count and VF gene set match; this was treated as current.")
}
if (is_stale(path_status_poster, path_status)) {
  ma("- **WARNING:** `status_map_with_poster_tp.csv` is missing or older than `status_map.csv`; rerun script 00d before transition interpretation.")
}
if (length(legacy_found) > 0) {
  ma("- Legacy/stale-looking summary files were detected and are listed in `summary_qc_log.txt`; use current `results/summary/` outputs as truth.")
}
ma("- UTI sample size in the VF-linked dataset is limited; score tests are descriptive/exploratory and underpowered for definitive UTI association.")
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
  table_08_vf_score_summary_by_ST = t08,
  table_09_longitudinal_vf_stability_summary = t09,
  table_10_asb_uti_transition_cases = t10,
  table_11_lineage_context_summary = t11,
  table_12_missing_data_audit = t12,
  table_13_optional_vf_amr_summary = t13,
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

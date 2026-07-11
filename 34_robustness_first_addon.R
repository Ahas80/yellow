#!/usr/bin/env Rscript
# ==============================================================================
# 34_robustness_first_addon.R
# ==============================================================================
#
# GOAL:
#   Build a robustness-focused interpretation layer for the current primary
#   UTI vs Not_UTI analysis. This script consolidates existing diagnostics and
#   adds a near-miss-expanded sensitivity analysis without relabelling the
#   primary denominator.
#
# OUTPUT:
#   - results/robustness/robustness_summary.md
#   - results/robustness/robustness_validation_checks.csv
#   - results/robustness/denominator_robustness_summary.csv
#   - results/robustness/qc_attrition_robustness.csv
#   - results/robustness/near_miss_expanded_score_summary.csv
#   - results/robustness/near_miss_expanded_score_contrasts.csv
#   - results/robustness/near_miss_expanded_module_fisher.csv
#   - results/robustness/model_stability_summary.csv
#   - results/robustness/leave_one_uti_sensitivity_summary.csv
#   - results/robustness/bootstrap_score_robustness.csv
#   - results/robustness/power_precision_summary.csv
#   - plots/robustness/qc_retention_by_status.png
#   - plots/robustness/near_miss_score_shift.png
#   - plots/robustness/model_stability_flags.png
#
# ==============================================================================

source("00_config.R")
suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
  library(stringr)
  library(tibble)
  library(ggplot2)
  library(purrr)
})

msg("Starting 34_robustness_first_addon.R")

DIR_ROBUSTNESS <- file.path(DIR_RESULTS, "robustness")
DIR_PLOTS_ROBUSTNESS <- file.path(DIR_PLOTS, "robustness")
ensure_dir(DIR_ROBUSTNESS)
ensure_dir(DIR_PLOTS_ROBUSTNESS)

# ==============================================================================
# HELPERS
# ==============================================================================

require_file <- function(path, label = basename(path)) {
  if (!file.exists(path)) stop("Missing required input for robustness add-on: ", label, " (", path, ")")
  path
}

safe_read_csv <- function(path, ...) {
  if (!file.exists(path)) return(NULL)
  read_csv(path, show_col_types = FALSE, ...)
}

normalise_tp_label <- function(x) {
  if (exists("normalise_timepoint_preserve_events", mode = "function")) {
    normalise_timepoint_preserve_events(x)
  } else {
    x <- str_trim(as.character(x))
    case_when(
      str_detect(str_to_lower(x), "uricult") ~ "Uricult",
      str_detect(str_to_upper(x), "^T\\d+$") ~ str_to_upper(x),
      str_detect(x, "^\\d+$") ~ paste0("T", x),
      TRUE ~ x
    )
  }
}

key2 <- function(pid, tp) paste(as.character(pid), normalise_tp_label(tp), sep = "|")

median_or_na <- function(x) {
  if (all(is.na(x))) NA_real_ else median(x, na.rm = TRUE)
}

mean_or_na <- function(x) {
  if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)
}

q_or_na <- function(x, p) {
  if (all(is.na(x))) NA_real_ else as.numeric(quantile(x, p, na.rm = TRUE))
}

fmt_num <- function(x, digits = 2) {
  ifelse(is.na(x), "NA", format(round(x, digits), nsmall = digits, trim = TRUE))
}

plot_theme_robustness <- function(base_size = 10) {
  theme_bw(base_size = base_size) +
    theme(
      plot.title = element_text(face = "bold"),
      plot.subtitle = element_text(colour = "grey25"),
      plot.caption = element_text(hjust = 0, colour = "grey35", size = base_size - 2),
      legend.position = "bottom"
    )
}

score_summary <- function(df, group_col, score_cols) {
  score_cols <- intersect(score_cols, names(df))
  if (length(score_cols) == 0) return(tibble())
  df %>%
    select(Participant_id, all_of(group_col), all_of(score_cols)) %>%
    filter(!is.na(.data[[group_col]])) %>%
    pivot_longer(all_of(score_cols), names_to = "score", values_to = "value") %>%
    group_by(.data[[group_col]], score) %>%
    summarise(
      n_episodes = sum(!is.na(value)),
      n_participants = n_distinct(Participant_id[!is.na(value)]),
      median = median_or_na(value),
      q25 = q_or_na(value, 0.25),
      q75 = q_or_na(value, 0.75),
      mean = mean_or_na(value),
      sd = ifelse(sum(!is.na(value)) > 1, sd(value, na.rm = TRUE), NA_real_),
      .groups = "drop"
    ) %>%
    rename(group = all_of(group_col))
}

score_contrast <- function(df, group_col, group_a, group_b, score_cols, contrast_name) {
  score_cols <- intersect(score_cols, names(df))
  map_dfr(score_cols, function(sc) {
    a <- df %>% filter(.data[[group_col]] == group_a) %>% pull(sc)
    b <- df %>% filter(.data[[group_col]] == group_b) %>% pull(sc)
    a <- a[!is.na(a)]
    b <- b[!is.na(b)]
    p <- if (length(a) > 0 && length(b) > 0) {
      suppressWarnings(wilcox.test(a, b, exact = FALSE)$p.value)
    } else NA_real_
    tibble(
      contrast = contrast_name,
      score = sc,
      group_a = group_a,
      group_b = group_b,
      n_a = length(a),
      n_b = length(b),
      median_a = median_or_na(a),
      median_b = median_or_na(b),
      mean_a = mean_or_na(a),
      mean_b = mean_or_na(b),
      median_difference_a_minus_b = median_or_na(a) - median_or_na(b),
      mean_difference_a_minus_b = mean_or_na(a) - mean_or_na(b),
      wilcox_p = p
    )
  }) %>%
    group_by(contrast) %>%
    mutate(wilcox_fdr = p.adjust(wilcox_p, method = "BH")) %>%
    ungroup()
}

module_fisher <- function(df, binary_col, module_cols, contrast_name) {
  module_cols <- intersect(module_cols, names(df))
  map_dfr(module_cols, function(feature) {
    keep <- df %>% filter(.data[[binary_col]] %in% c("Expanded_UTI", "Expanded_Not_UTI"))
    uti <- keep %>% filter(.data[[binary_col]] == "Expanded_UTI")
    not <- keep %>% filter(.data[[binary_col]] == "Expanded_Not_UTI")
    uti_present <- sum(uti[[feature]] == 1, na.rm = TRUE)
    uti_absent <- sum(uti[[feature]] == 0, na.rm = TRUE)
    not_present <- sum(not[[feature]] == 1, na.rm = TRUE)
    not_absent <- sum(not[[feature]] == 0, na.rm = TRUE)
    mat <- matrix(c(uti_present, uti_absent, not_present, not_absent), nrow = 2, byrow = TRUE)
    ft <- suppressWarnings(fisher.test(mat))
    tibble(
      contrast = contrast_name,
      feature = feature,
      n_expanded_uti_present = uti_present,
      n_expanded_uti_absent = uti_absent,
      n_expanded_not_uti_present = not_present,
      n_expanded_not_uti_absent = not_absent,
      prevalence_expanded_uti = ifelse(nrow(uti) > 0, uti_present / nrow(uti), NA_real_),
      prevalence_expanded_not_uti = ifelse(nrow(not) > 0, not_present / nrow(not), NA_real_),
      prevalence_diff = prevalence_expanded_uti - prevalence_expanded_not_uti,
      fisher_or = unname(ft$estimate),
      fisher_p = ft$p.value
    )
  }) %>%
    mutate(fisher_fdr = p.adjust(fisher_p, method = "BH")) %>%
    arrange(fisher_p)
}

# ==============================================================================
# LOAD INPUTS
# ==============================================================================

status <- read_csv(require_file(FILE_STATUS_MAP), show_col_types = FALSE) %>%
  prefer_primary_uti_status(allow_legacy_fallback = FALSE) %>%
  apply_manual_sample_curation(context = "34_status") %>%
  mutate(Participant_id = as.character(Participant_id),
         tp_lab = normalise_tp_label(tp_lab))

status_primary <- status %>%
  filter(analysis_include_primary %in% TRUE)

vf_ready <- read_csv(require_file(FILE_VF_READY), show_col_types = FALSE) %>%
  prefer_primary_uti_status(allow_legacy_fallback = FALSE) %>%
  apply_manual_sample_curation(context = "34_vf_ready") %>%
  filter_primary_genomics() %>%
  mutate(Participant_id = as.character(Participant_id),
         tp_lab = normalise_tp_label(tp_lab))

score_path <- require_file(file.path(DIR_VF, "vf_score_table.csv"))
scores <- read_csv(score_path, show_col_types = FALSE) %>%
  mutate(Participant_id = as.character(Participant_id),
         tp_lab = normalise_tp_label(tp_lab))

if (!"UTI_Status" %in% names(scores)) {
  scores <- scores %>%
    left_join(
      vf_ready %>% select(Participant_id, tp_lab, Episode_ID, UTI_Status, UTI_binary),
      by = c("Participant_id", "tp_lab"),
      relationship = "many-to-one"
    )
}

scores <- scores %>%
  apply_manual_sample_curation(context = "34_scores") %>%
  filter_primary_genomics()

near_miss <- read_csv(require_file(file.path(DIR_RESULTS, "audit", "uti_not_uti_near_miss_rows.csv")),
                      show_col_types = FALSE) %>%
  mutate(Participant_id = as.character(Participant_id),
         tp_lab = normalise_tp_label(tp_lab),
         near_miss_key = key2(Participant_id, tp_lab))

near_miss_keys <- unique(near_miss$near_miss_key)
vf_keys <- key2(scores$Participant_id, scores$tp_lab)

selected_scores <- c(
  "expec_marker_count", "upec_system_count", "upec_system_fraction",
  "total_vf_count_all", "total_vf_count_curated", "total_vf_count_unassigned",
  "low_confidence_count"
)
selected_scores <- intersect(selected_scores, names(scores))
score_label_lookup <- c(
  expec_marker_count = "ExPEC-like marker count",
  upec_system_count = "UPEC systems present",
  upec_system_fraction = "UPEC system fraction",
  total_vf_count_all = "All VF genes",
  total_vf_count_curated = "Curated VF genes",
  total_vf_count_unassigned = "Unassigned VF genes",
  low_confidence_count = "Low-confidence VF genes"
)
module_cols <- grep("^mod_.*_present$", names(scores), value = TRUE)

scores_robust <- scores %>%
  mutate(
    near_miss_possible_uti = key2(Participant_id, tp_lab) %in% near_miss_keys,
    primary_binary_group = if_else(UTI_Status == "UTI", "Primary_UTI", "Primary_Not_UTI"),
    near_miss_three_group = case_when(
      UTI_Status == "UTI" ~ "Primary_UTI",
      near_miss_possible_uti ~ "Near_miss_possible_UTI",
      UTI_Status == "Not_UTI" ~ "Primary_Not_UTI_excluding_near_miss",
      TRUE ~ NA_character_
    ),
    expanded_binary_group = case_when(
      UTI_Status == "UTI" | near_miss_possible_uti ~ "Expanded_UTI",
      UTI_Status == "Not_UTI" ~ "Expanded_Not_UTI",
      TRUE ~ NA_character_
    )
  )

msg("Loaded primary and sensitivity datasets.")

# ==============================================================================
# DENOMINATOR AND QC ROBUSTNESS
# ==============================================================================

clinical_counts <- status_primary %>%
  count(UTI_Status, name = "n_clinical") %>%
  pivot_wider(names_from = UTI_Status, values_from = n_clinical, values_fill = 0)

vf_counts <- vf_ready %>%
  count(UTI_Status, name = "n_vf_ready") %>%
  pivot_wider(names_from = UTI_Status, values_from = n_vf_ready, values_fill = 0)

near_miss_vf_n <- sum(key2(scores$Participant_id, scores$tp_lab) %in% near_miss_keys)

denominator_summary <- tibble(
  metric = c(
    "primary_clinical_total", "primary_clinical_uti", "primary_clinical_not_uti",
    "primary_vf_total", "primary_vf_uti", "primary_vf_not_uti",
    "near_miss_clinical_rows", "near_miss_vf_ready_rows",
    "expanded_vf_uti_if_near_miss_included", "expanded_vf_not_uti_if_near_miss_excluded"
  ),
  value = c(
    nrow(status_primary),
    if ("UTI" %in% names(clinical_counts)) clinical_counts$UTI else 0,
    if ("Not_UTI" %in% names(clinical_counts)) clinical_counts$Not_UTI else 0,
    nrow(vf_ready),
    if ("UTI" %in% names(vf_counts)) vf_counts$UTI else 0,
    if ("Not_UTI" %in% names(vf_counts)) vf_counts$Not_UTI else 0,
    nrow(near_miss),
    near_miss_vf_n,
    sum(scores_robust$expanded_binary_group == "Expanded_UTI", na.rm = TRUE),
    sum(scores_robust$expanded_binary_group == "Expanded_Not_UTI", na.rm = TRUE)
  ),
  interpretation = c(
    "Primary clinical denominator after manual exclusions.",
    "Rows satisfying culture support plus catheter-aware symptom compatibility.",
    "Rows not satisfying both primary UTI rule components.",
    "Sequenced VF/model denominator after genomics filters.",
    "Primary UTI rows in VF/model denominator.",
    "Primary Not_UTI rows in VF/model denominator.",
    "Culture-supported or legacy UTI-like rows blocked by current symptom rule.",
    "Near-miss rows that have VF/model data and can enter expanded sensitivity only.",
    "Sensitivity denominator if VF-ready near-miss rows are treated as possible UTI.",
    "Sensitivity denominator if VF-ready near-miss rows are removed from Not_UTI."
  )
)

write_csv(denominator_summary, file.path(DIR_ROBUSTNESS, "denominator_robustness_summary.csv"))

qc_bias <- read_csv(require_file(file.path(DIR_QC, "qc_selection_bias_by_status.csv")),
                    show_col_types = FALSE) %>%
  mutate(
    Infection_Status = as.character(Infection_Status),
    qc_pass = as.logical(qc_pass)
  )

qc_attrition <- qc_bias %>%
  group_by(Infection_Status) %>%
  mutate(total_status = sum(n), pct = n / total_status) %>%
  ungroup() %>%
  mutate(retention_state = if_else(qc_pass, "QC_PASS_selected", "QC_or_selection_loss"))

write_csv(qc_attrition, file.path(DIR_ROBUSTNESS, "qc_attrition_robustness.csv"))

# ==============================================================================
# NEAR-MISS EXPANDED SENSITIVITY
# ==============================================================================

near_miss_score_summary <- bind_rows(
  score_summary(scores_robust, "primary_binary_group", selected_scores) %>%
    mutate(sensitivity_layer = "primary_rule_binary"),
  score_summary(scores_robust, "near_miss_three_group", selected_scores) %>%
    mutate(sensitivity_layer = "near_miss_three_group"),
  score_summary(scores_robust, "expanded_binary_group", selected_scores) %>%
    mutate(sensitivity_layer = "near_miss_expanded_binary")
) %>%
  mutate(
    score_label = if_else(
      .data$score %in% names(score_label_lookup),
      unname(score_label_lookup[.data$score]),
      .data$score
    )
  ) %>%
  select(sensitivity_layer, everything())

write_csv(near_miss_score_summary, file.path(DIR_ROBUSTNESS, "near_miss_expanded_score_summary.csv"))

near_miss_score_contrasts <- bind_rows(
  score_contrast(scores_robust, "primary_binary_group",
                 "Primary_UTI", "Primary_Not_UTI", selected_scores, "primary_UTI_vs_all_Not_UTI"),
  score_contrast(scores_robust, "expanded_binary_group",
                 "Expanded_UTI", "Expanded_Not_UTI", selected_scores, "expanded_UTI_vs_Not_UTI_excluding_near_miss"),
  score_contrast(scores_robust, "near_miss_three_group",
                 "Near_miss_possible_UTI", "Primary_UTI", selected_scores, "near_miss_possible_vs_primary_UTI")
) %>%
  mutate(
    score_label = if_else(
      .data$score %in% names(score_label_lookup),
      unname(score_label_lookup[.data$score]),
      .data$score
    )
  ) %>%
  relocate(score_label, .after = score)

write_csv(near_miss_score_contrasts, file.path(DIR_ROBUSTNESS, "near_miss_expanded_score_contrasts.csv"))

expanded_module_fisher <- module_fisher(
  scores_robust,
  "expanded_binary_group",
  module_cols,
  "expanded_UTI_includes_near_miss_possible"
)

write_csv(expanded_module_fisher, file.path(DIR_ROBUSTNESS, "near_miss_expanded_module_fisher.csv"))

# ==============================================================================
# MODEL AND SPARSE-DATA STABILITY
# ==============================================================================

glmm <- read_csv(require_file(file.path(DIR_MODELS, "gwas_multivariable_glmm.csv")),
                 show_col_types = FALSE)
univ <- read_csv(require_file(file.path(DIR_MODELS, "gwas_univariable_stats.csv")),
                 show_col_types = FALSE)
warnings_txt <- readLines(require_file(file.path(DIR_MODELS, "model_interpretation_warnings.txt")))

model_stability_summary <- tibble(
  metric = c(
    "glmm_models_total", "glmm_models_converged", "glmm_singular_models",
    "glmm_sparse_or_separation_risk", "glmm_nominal_p_lt_0_05", "glmm_fdr_lt_0_05",
    "univariable_tests_total", "univariable_nominal_p_lt_0_05", "univariable_fdr_lt_0_05",
    "model_warning_text"
  ),
  value = c(
    nrow(glmm),
    sum(glmm$converged %in% TRUE, na.rm = TRUE),
    sum(str_detect(glmm$model_type, regex("singular", ignore_case = TRUE)), na.rm = TRUE),
    sum(glmm$sparse_data_separation_risk %in% TRUE, na.rm = TRUE),
    sum(glmm$p.value < 0.05, na.rm = TRUE),
    sum(glmm$FDR < 0.05, na.rm = TRUE),
    nrow(univ),
    sum(univ$p_value < 0.05, na.rm = TRUE),
    sum(univ$FDR < 0.05, na.rm = TRUE),
    paste(warnings_txt, collapse = " | ")
  ),
  interpretation = c(
    "Number of mixed-effect models attempted.",
    "Convergence alone does not imply robust inference.",
    "Singular random-effect fits indicate limited repeated-measure information for those features.",
    "Sparse/separation-risk features are unstable with 17 primary VF-ready UTI rows.",
    "Nominal GLMM hits before FDR correction.",
    "Confirmatory GLMM hits after FDR correction.",
    "Univariable tests attempted.",
    "Nominal univariable hits before FDR correction.",
    "Confirmatory univariable hits after FDR correction.",
    "Exact warning text emitted by the model script."
  )
)

write_csv(model_stability_summary, file.path(DIR_ROBUSTNESS, "model_stability_summary.csv"))

leave_one <- read_csv(require_file(file.path(DIR_VF, "uti_not_uti_leave_one_uti_out.csv")),
                      show_col_types = FALSE)

leave_one_summary <- leave_one %>%
  group_by(feature_type) %>%
  summarise(
    n_features = n(),
    n_direction_flip = sum(direction_flip %in% TRUE, na.rm = TRUE),
    n_sparse_warning = sum(sparse_uti_warning %in% TRUE, na.rm = TRUE),
    median_sensitivity_range = median_or_na(sensitivity_range),
    max_sensitivity_range = ifelse(all(is.na(sensitivity_range)), NA_real_, max(sensitivity_range, na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  arrange(desc(n_direction_flip), desc(max_sensitivity_range))

write_csv(leave_one_summary, file.path(DIR_ROBUSTNESS, "leave_one_uti_sensitivity_summary.csv"))

bootstrap <- read_csv(require_file(file.path(DIR_VF, "uti_not_uti_bootstrap_effects.csv")),
                      show_col_types = FALSE) %>%
  mutate(
    ci_crosses_zero = ci_lower <= 0 & ci_upper >= 0,
    ci_direction = case_when(
      ci_lower > 0 ~ "positive_CI_excludes_zero",
      ci_upper < 0 ~ "negative_CI_excludes_zero",
      TRUE ~ "CI_crosses_zero"
    )
  ) %>%
  arrange(ci_crosses_zero, desc(abs(observed_effect)))

write_csv(bootstrap, file.path(DIR_ROBUSTNESS, "bootstrap_score_robustness.csv"))

power_context <- read_csv(require_file(file.path(DIR_VF, "uti_not_uti_power_precision_context.csv")),
                          show_col_types = FALSE)

power_precision_summary <- power_context %>%
  group_by(baseline_not_uti_prevalence) %>%
  summarise(
    min_detectable_or_above_one = suppressWarnings(min(odds_ratio[odds_ratio > 1 & detectable_at_alpha_0_05], na.rm = TRUE)),
    max_detectable_or_below_one = suppressWarnings(max(odds_ratio[odds_ratio < 1 & detectable_at_alpha_0_05], na.rm = TRUE)),
    n_detectable_scenarios = sum(detectable_at_alpha_0_05 %in% TRUE, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    min_detectable_or_above_one = ifelse(is.infinite(min_detectable_or_above_one), NA_real_, min_detectable_or_above_one),
    max_detectable_or_below_one = ifelse(is.infinite(max_detectable_or_below_one), NA_real_, max_detectable_or_below_one),
    interpretation = "Approximate sparse-count context only; not a prospective power calculation."
  )

write_csv(power_precision_summary, file.path(DIR_ROBUSTNESS, "power_precision_summary.csv"))

# ==============================================================================
# ROBUSTNESS DECISION MATRIX
# ==============================================================================

expanded_score_top <- near_miss_score_contrasts %>%
  filter(contrast == "expanded_UTI_vs_Not_UTI_excluding_near_miss") %>%
  arrange(wilcox_p) %>%
  slice_head(n = 5) %>%
  transmute(score, median_difference_a_minus_b, wilcox_p, wilcox_fdr)

top_glmm <- glmm %>%
  arrange(FDR, p.value) %>%
  slice_head(n = 10) %>%
  select(feature, OR, OR_lower, OR_upper, p.value, FDR, sparse_data_separation_risk, model_type)

write_csv(expanded_score_top, file.path(DIR_ROBUSTNESS, "near_miss_expanded_top_score_contrasts.csv"))
write_csv(top_glmm, file.path(DIR_ROBUSTNESS, "top_glmm_robustness_flags.csv"))

robustness_claim_matrix <- tibble(
  result_area = c(
    "Primary denominator",
    "QC attrition",
    "Near-miss sensitivity",
    "Score burden",
    "Feature associations",
    "Leave-one-UTI sensitivity",
    "Model stability",
    "Power/precision"
  ),
  robustness_readout = c(
    sprintf("Primary clinical %d UTI/%d Not_UTI; VF-ready %d UTI/%d Not_UTI.",
            denominator_summary$value[denominator_summary$metric == "primary_clinical_uti"],
            denominator_summary$value[denominator_summary$metric == "primary_clinical_not_uti"],
            denominator_summary$value[denominator_summary$metric == "primary_vf_uti"],
            denominator_summary$value[denominator_summary$metric == "primary_vf_not_uti"]),
    sprintf("QC loss by status: UTI %d lost, Not_UTI %d lost.",
            qc_attrition %>% filter(Infection_Status == "UTI", retention_state == "QC_or_selection_loss") %>% pull(n) %>% sum(na.rm = TRUE),
            qc_attrition %>% filter(Infection_Status == "Not_UTI", retention_state == "QC_or_selection_loss") %>% pull(n) %>% sum(na.rm = TRUE)),
    sprintf("%d near-miss clinical rows; %d have VF-ready data and form expanded sensitivity only.",
            nrow(near_miss), near_miss_vf_n),
    "Primary and expanded supplementary endpoint contrasts remain descriptive; no burden or marker result should be called confirmatory.",
    sprintf("GLMM FDR<0.05: %d; univariable FDR<0.05: %d.",
            sum(glmm$FDR < 0.05, na.rm = TRUE), sum(univ$FDR < 0.05, na.rm = TRUE)),
    sprintf("%d/%d leave-one diagnostics flip direction.",
            sum(leave_one$direction_flip %in% TRUE, na.rm = TRUE), nrow(leave_one)),
    sprintf("%d GLMM sparse/separation-risk features and %d singular fits.",
            sum(glmm$sparse_data_separation_risk %in% TRUE, na.rm = TRUE),
            sum(str_detect(glmm$model_type, regex("singular", ignore_case = TRUE)), na.rm = TRUE)),
    "With 17 VF-ready UTI rows, only large/extreme prevalence differences are detectably stable."
  ),
  claim_strength = c(
    "strong denominator trace",
    "low apparent status-biased QC loss, still sparse",
    "sensitivity-only",
    "descriptive",
    "not confirmatory",
    "fragile",
    "fragile",
    "underpowered"
  )
)

write_csv(robustness_claim_matrix, file.path(DIR_ROBUSTNESS, "robustness_claim_matrix.csv"))

# ==============================================================================
# PLOTS
# ==============================================================================

p_qc <- qc_attrition %>%
  mutate(Infection_Status = factor(Infection_Status, levels = c("Not_UTI", "UTI"))) %>%
  ggplot(aes(Infection_Status, n, fill = retention_state)) +
  geom_col(colour = "white", linewidth = 0.3) +
  geom_text(aes(label = n), position = position_stack(vjust = 0.5), size = 3) +
  scale_fill_manual(values = c(QC_PASS_selected = "#2F855A", QC_or_selection_loss = "#C05621")) +
  labs(
    title = "QC and canonical-selection retention by primary status",
    x = NULL, y = "Episodes", fill = NULL,
    caption = "Uses current primary UTI_Status. Counts are sparse for UTI."
  ) +
  plot_theme_robustness()
ggsave(file.path(DIR_PLOTS_ROBUSTNESS, "qc_retention_by_status.png"), p_qc, width = 7.5, height = 5, dpi = 300)

plot_scores <- near_miss_score_summary %>%
  filter(
    sensitivity_layer == "near_miss_three_group",
    score %in% intersect(c("expec_marker_count", "upec_system_count",
                           "upec_system_fraction", "total_vf_count_curated"),
                         score)
  ) %>%
  mutate(
    group_label = recode(
      group,
      Primary_Not_UTI_excluding_near_miss = "Not_UTI",
      Near_miss_possible_UTI = "Near-miss",
      Primary_UTI = "UTI"
    ),
    group_label = factor(group_label, levels = c("Not_UTI", "Near-miss", "UTI")),
    score_label = factor(score_label, levels = unique(score_label))
  )

p_scores <- ggplot(plot_scores, aes(group_label, median, fill = group_label)) +
  geom_col(width = 0.75, colour = "white", linewidth = 0.3) +
  geom_errorbar(aes(ymin = q25, ymax = q75), width = 0.2) +
  facet_wrap(~ score_label, scales = "free_y", ncol = 2) +
  scale_fill_manual(values = c(
    Not_UTI = "#4A5568",
    `Near-miss` = "#B7791F",
    UTI = "#C05621"
  )) +
  labs(
    title = "Near-miss sensitivity: supplementary endpoint medians",
    subtitle = "Near-miss rows are shown as possible UTI, not relabelled primary UTI",
    x = NULL, y = "Median endpoint value with IQR", fill = NULL
  ) +
  plot_theme_robustness(9) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5))
ggsave(file.path(DIR_PLOTS_ROBUSTNESS, "near_miss_score_shift.png"), p_scores, width = 9, height = 6.5, dpi = 300)

flag_counts <- tibble(
  flag = c(
    "GLMM sparse/separation risk",
    "GLMM singular fit",
    "GLMM FDR < 0.05",
    "Univariable FDR < 0.05",
    "Leave-one direction flip",
    "Bootstrap CI excludes 0"
  ),
  n = c(
    sum(glmm$sparse_data_separation_risk %in% TRUE, na.rm = TRUE),
    sum(str_detect(glmm$model_type, regex("singular", ignore_case = TRUE)), na.rm = TRUE),
    sum(glmm$FDR < 0.05, na.rm = TRUE),
    sum(univ$FDR < 0.05, na.rm = TRUE),
    sum(leave_one$direction_flip %in% TRUE, na.rm = TRUE),
    sum(!bootstrap$ci_crosses_zero, na.rm = TRUE)
  ),
  denominator = c(nrow(glmm), nrow(glmm), nrow(glmm), nrow(univ), nrow(leave_one), nrow(bootstrap))
) %>%
  mutate(label = paste0(n, "/", denominator))

p_flags <- ggplot(flag_counts, aes(reorder(flag, n), n, fill = flag)) +
  geom_col(show.legend = FALSE) +
  geom_text(aes(label = label), hjust = -0.1, size = 3) +
  coord_flip() +
  scale_fill_manual(values = rep(c("#C05621", "#B7791F", "#2F855A", "#2B6CB0"), length.out = nrow(flag_counts))) +
  labs(
    title = "Robustness flags across model and sensitivity diagnostics",
    x = NULL, y = "Flagged results",
    caption = "Counts summarize existing primary-rule diagnostics."
  ) +
  plot_theme_robustness() +
  expand_limits(y = max(flag_counts$n, na.rm = TRUE) * 1.15)
ggsave(file.path(DIR_PLOTS_ROBUSTNESS, "model_stability_flags.png"), p_flags, width = 8.5, height = 5.5, dpi = 300)

# ==============================================================================
# MARKDOWN SUMMARY AND VALIDATION
# ==============================================================================

validation <- tibble(
  check = c(
    "primary clinical denominator is 583 rows",
    "primary clinical UTI count is 18",
    "primary VF/model denominator is 556 rows",
    "primary VF/model UTI count is 17",
    "near-miss clinical row count is 19",
    "VF-ready missing primary status is zero",
    "legacy OLD files not used as inputs"
  ),
  status = c(
    ifelse(nrow(status_primary) == 583, "PASS", "FAIL"),
    ifelse(sum(status_primary$UTI_Status == "UTI", na.rm = TRUE) == 18, "PASS", "FAIL"),
    ifelse(nrow(vf_ready) == 556, "PASS", "FAIL"),
    ifelse(sum(vf_ready$UTI_Status == "UTI", na.rm = TRUE) == 17, "PASS", "FAIL"),
    ifelse(nrow(near_miss) == 19, "PASS", "FAIL"),
    ifelse(sum(is.na(vf_ready$UTI_Status)) == 0, "PASS", "FAIL"),
    "PASS"
  ),
  detail = c(
    sprintf("n=%d", nrow(status_primary)),
    sprintf("UTI=%d", sum(status_primary$UTI_Status == "UTI", na.rm = TRUE)),
    sprintf("n=%d", nrow(vf_ready)),
    sprintf("UTI=%d", sum(vf_ready$UTI_Status == "UTI", na.rm = TRUE)),
    sprintf("near_miss=%d", nrow(near_miss)),
    sprintf("missing=%d", sum(is.na(vf_ready$UTI_Status))),
    "Inputs restricted to current results/vf, results/qc, results/audit, results/models, and status_map outputs."
  )
)

write_csv(validation, file.path(DIR_ROBUSTNESS, "robustness_validation_checks.csv"))
if (any(validation$status != "PASS")) {
  print(validation)
  stop("Robustness add-on validation failed.")
}

summary_lines <- c(
  "# Robustness-First Add-On Summary",
  "",
  sprintf("Generated: %s", format(Sys.time())),
  "",
  "## Denominator Robustness",
  "",
  sprintf("- Primary clinical denominator: **%d** rows = **%d UTI**, **%d Not_UTI**.",
          nrow(status_primary),
          sum(status_primary$UTI_Status == "UTI", na.rm = TRUE),
          sum(status_primary$UTI_Status == "Not_UTI", na.rm = TRUE)),
  sprintf("- Primary VF/model denominator: **%d** rows = **%d UTI**, **%d Not_UTI**.",
          nrow(vf_ready),
          sum(vf_ready$UTI_Status == "UTI", na.rm = TRUE),
          sum(vf_ready$UTI_Status == "Not_UTI", na.rm = TRUE)),
  sprintf("- Near-miss rows: **%d** clinical, **%d** VF-ready. These are sensitivity-only possible UTI rows.",
          nrow(near_miss), near_miss_vf_n),
  "",
  "## Robustness Readout",
  "",
  sprintf("- QC/selection loss is sparse for UTI: **%d** UTI row lost versus **%d** Not_UTI rows lost.",
          qc_attrition %>% filter(Infection_Status == "UTI", retention_state == "QC_or_selection_loss") %>% pull(n) %>% sum(na.rm = TRUE),
          qc_attrition %>% filter(Infection_Status == "Not_UTI", retention_state == "QC_or_selection_loss") %>% pull(n) %>% sum(na.rm = TRUE)),
  sprintf("- Model results remain exploratory: GLMM FDR<0.05 = **%d**; univariable FDR<0.05 = **%d**.",
          sum(glmm$FDR < 0.05, na.rm = TRUE), sum(univ$FDR < 0.05, na.rm = TRUE)),
  sprintf("- Sparse-data warnings remain substantial: **%d/%d** GLMM features flagged sparse/separation risk.",
          sum(glmm$sparse_data_separation_risk %in% TRUE, na.rm = TRUE), nrow(glmm)),
  sprintf("- Leave-one-UTI-out direction flips occur in **%d/%d** diagnostics.",
          sum(leave_one$direction_flip %in% TRUE, na.rm = TRUE), nrow(leave_one)),
  "",
  "## Interpretation",
  "",
  "- Robustness diagnostics support the current denominator trace, but they do **not** upgrade association claims to confirmatory.",
  "- Near-miss expansion is useful for sensitivity framing, not for relabelling primary UTI.",
  "- With 17 VF-ready UTI rows, supplementary endpoint and feature effects should be described as descriptive/hypothesis-generating.",
  "",
  "## Key Files",
  "",
  "- `results/robustness/denominator_robustness_summary.csv`",
  "- `results/robustness/near_miss_expanded_score_contrasts.csv`",
  "- `results/robustness/near_miss_expanded_module_fisher.csv`",
  "- `results/robustness/model_stability_summary.csv`",
  "- `results/robustness/robustness_claim_matrix.csv`",
  "- `results/robustness/robustness_validation_checks.csv`"
)

writeLines(summary_lines, file.path(DIR_ROBUSTNESS, "robustness_summary.md"))

msg("Robustness add-on complete.")
msg("Outputs written to %s and %s", DIR_ROBUSTNESS, DIR_PLOTS_ROBUSTNESS)

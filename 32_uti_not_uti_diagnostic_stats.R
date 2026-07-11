#!/usr/bin/env Rscript

# ==============================================================================
# 32_uti_not_uti_diagnostic_stats.R
# ==============================================================================
#
# Goal:
#   Add a conservative diagnostic layer for the primary UTI vs Not_UTI
#   denominator. The script does not change labels or inclusion rules. It
#   explains the current counts, quantifies sparse-count uncertainty, and writes
#   interpretation-ready figures/tables for the final summary.
# ==============================================================================

source("00_config.R")

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
  library(stringr)
  library(purrr)
  library(ggplot2)
})

set.seed(20260522)

DIR_AUDIT <- file.path(DIR_RESULTS, "audit")
ensure_dir(DIR_AUDIT)
ensure_dir(DIR_VF)
ensure_dir(DIR_PLOTS_CLINICAL)
ensure_dir(DIR_PLOTS_VF)

msg("Starting 32_uti_not_uti_diagnostic_stats.R")

safe_read <- function(path, required = FALSE) {
  if (!file.exists(path)) {
    if (required) stop("Missing required input: ", path)
    return(tibble())
  }
  readr::read_csv(path, show_col_types = FALSE)
}

safe_num <- function(x) suppressWarnings(as.numeric(x))

collapse_nonmissing <- function(x) {
  x <- unique(as.character(x))
  x <- x[!is.na(x) & nzchar(x) & x != "."]
  if (length(x) == 0) NA_character_ else paste(x, collapse = "; ")
}

clean_keyed <- function(df, context, genomics = FALSE, primary = TRUE) {
  if (nrow(df) == 0) return(df)
  df <- tibble::as_tibble(df)
  if (!"tp_lab" %in% names(df) && "Timepoint" %in% names(df)) {
    df$tp_lab <- normalise_timepoint_preserve_events(df$Timepoint)
  }
  if (!"tp_lab" %in% names(df)) {
    stop("Input in ", context, " has no tp_lab or Timepoint column.")
  }
  df <- df %>%
    mutate(
      Participant_id = as.character(.data$Participant_id),
      tp_lab = as.character(.data$tp_lab)
    ) %>%
    prefer_primary_uti_status() %>%
    apply_manual_sample_curation(context = context)

  if (genomics) {
    df <- filter_primary_genomics(df)
  } else if (primary) {
    df <- filter_primary_analysis(df)
  }
  df
}

save_plot <- function(plot, path, width = 9, height = 5.5, dpi = 300) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  ggplot2::ggsave(path, plot, width = width, height = height, dpi = dpi, bg = "white")
  invisible(path)
}

empty_plot <- function(title, subtitle = "No eligible rows under the current primary analysis filters.") {
  ggplot() +
    annotate("text", x = 0, y = 0, label = subtitle, size = 4.2, colour = "grey30") +
    labs(title = title) +
    theme_void(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, margin = margin(b = 12))
    )
}

exact_p <- function(expr) {
  tryCatch(expr, error = function(e) NA_real_, warning = function(w) suppressWarnings(expr))
}

participant_without_leading_zero <- function(x) {
  out <- str_replace(as.character(x), "^0+", "")
  ifelse(is.na(out) | !nzchar(out), as.character(x), out)
}

# ------------------------------------------------------------------------------
# Inputs and primary denominators
# ------------------------------------------------------------------------------

status_all <- safe_read(FILE_STATUS_MAP, required = TRUE) %>%
  mutate(
    Participant_id = as.character(.data$Participant_id),
    tp_lab = as.character(.data$tp_lab)
  ) %>%
  prefer_primary_uti_status(allow_legacy_fallback = FALSE) %>%
  apply_manual_sample_curation(context = "32_status_all")

status_primary <- filter_primary_analysis(status_all)

vf_all <- safe_read(FILE_VF_READY, required = TRUE) %>%
  mutate(
    Participant_id = as.character(.data$Participant_id),
    tp_lab = as.character(.data$tp_lab)
  ) %>%
  prefer_primary_uti_status() %>%
  apply_manual_sample_curation(context = "32_vf_ready_all")

vf <- filter_primary_genomics(vf_all)

score_file <- file.path(DIR_VF, "vf_score_table.csv")
score_raw <- safe_read(score_file)
score_df <- if (nrow(score_raw) > 0) {
  clean_keyed(score_raw, context = "32_score_table", genomics = TRUE)
} else {
  vf
}

model_denom <- safe_read(file.path(DIR_MODELS, "model_dataset_denominator.csv")) %>%
  { if (nrow(.) > 0) clean_keyed(., context = "32_model_denominator", genomics = TRUE) else . }

manual_exclusions <- safe_read(file.path(DIR_AUDIT, "primary_clinical_manual_exclusions.csv"))
vf_manual_exclusions <- safe_read(file.path(DIR_AUDIT, "primary_vf_manual_exclusions.csv"))
quarantined_fastas <- safe_read(file.path(DIR_AUDIT, "quarantined_failed_or_not_expected_fastas.csv"))
transition_score_changes <- safe_read(file.path(DIR_VF, "vf_transition_score_changes.csv"))

clinical_counts <- status_primary %>%
  count(UTI_Status, name = "n") %>%
  tidyr::complete(UTI_Status = c("UTI", "Not_UTI"), fill = list(n = 0L))

vf_counts <- vf %>%
  mutate(UTI_Status = if_else(is.na(.data$UTI_Status), "<missing>", .data$UTI_Status)) %>%
  count(UTI_Status, name = "n") %>%
  tidyr::complete(UTI_Status = c("UTI", "Not_UTI", "<missing>"), fill = list(n = 0L))

expected_clinical <- c(UTI = 18L, Not_UTI = 565L)
expected_vf <- c(UTI = 17L, Not_UTI = 539L, `<missing>` = 0L)

clinical_named <- setNames(clinical_counts$n, clinical_counts$UTI_Status)
vf_named <- setNames(vf_counts$n, vf_counts$UTI_Status)
if (!all(clinical_named[names(expected_clinical)] == expected_clinical)) {
  stop("Primary clinical counts changed unexpectedly. Expected 18 UTI and 565 Not_UTI.")
}
if (!all(vf_named[names(expected_vf)] == expected_vf)) {
  stop("Primary VF counts changed unexpectedly. Expected 17 UTI, 539 Not_UTI, and 0 missing.")
}

score_candidates <- c(
  "expec_marker_count", "upec_system_count", "upec_system_fraction",
  "total_vf_count_all", "total_vf_count_curated",
  "total_vf_count_unassigned", "low_confidence_count"
)
score_cols <- intersect(score_candidates, names(score_df))
score_cols <- score_cols[vapply(score_df[score_cols], is.numeric, logical(1))]
score_cols <- unique(score_cols)
if (length(score_cols) == 0) {
  stop("No numeric supplementary VF endpoint columns were found for diagnostics.")
}
score_label_lookup <- c(
  expec_marker_count = "ExPEC-like marker count",
  upec_system_count = "UPEC systems present",
  upec_system_fraction = "UPEC system fraction",
  total_vf_count_all = "All VFDB genes (descriptive burden)",
  total_vf_count_curated = "Curated VF genes (descriptive burden)",
  total_vf_count_unassigned = "Unassigned VFDB genes (QC burden)",
  low_confidence_count = "Low-confidence VF genes (QC burden)"
)

score_small <- score_df %>%
  select(any_of(c("Participant_id", "tp_lab", "Episode_ID", "UTI_Status", score_cols))) %>%
  mutate(across(all_of(score_cols), safe_num))

# ------------------------------------------------------------------------------
# Clinical rule decision flow and near-miss table
# ------------------------------------------------------------------------------

clinical_flow <- tibble::tribble(
  ~stage_order, ~stage, ~n, ~stage_group, ~interpretation,
  1L, "Raw status-map rows", nrow(status_all), "source",
  "All classified rows retained for audit before primary manual exclusions.",
  2L, "Primary clinical rows", nrow(status_primary), "primary",
  "Rows with analysis_include_primary == TRUE.",
  3L, "Culture supports UTI", sum(status_primary$culture_supports_uti %in% TRUE, na.rm = TRUE), "rule",
  "Primary rows where the culture side of the UTI rule is met.",
  4L, "Symptom compatible", sum(status_primary$symptom_compatible_uti %in% TRUE, na.rm = TRUE), "rule",
  "Primary rows where the catheter-aware symptom side of the UTI rule is met.",
  5L, "Primary UTI", sum(status_primary$UTI_Status == "UTI", na.rm = TRUE), "final",
  "Rows satisfying both culture support and symptom compatibility.",
  6L, "Primary Not_UTI", sum(status_primary$UTI_Status == "Not_UTI", na.rm = TRUE), "final",
  "All primary rows not satisfying both UTI rule components.",
  7L, "Primary VF/model rows", nrow(vf), "genomics",
  "VF-ready rows with analysis_include_primary == TRUE and genomics_expected_include == TRUE.",
  8L, "Missing VF status", sum(is.na(vf$UTI_Status)), "guardrail",
  "VF/model rows without a primary clinical status after matching."
)
write_csv(clinical_flow, file.path(DIR_AUDIT, "uti_not_uti_denominator_flow.csv"))

flow_plot <- clinical_flow %>%
  mutate(stage = factor(.data$stage, levels = .data$stage)) %>%
  ggplot(aes(x = stage, y = n, fill = stage_group)) +
  geom_col(width = 0.68, colour = "grey25", linewidth = 0.2) +
  geom_text(aes(label = n), vjust = -0.35, size = 3.4) +
  scale_fill_manual(
    values = c(
      source = "#6B7280", primary = "#2563EB", rule = "#0F766E",
      final = "#B45309", genomics = "#7C3AED", guardrail = "#BE123C"
    ),
    guide = "none"
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(
    title = "Primary UTI / Not_UTI Decision Flow",
    subtitle = "Counts are descriptive checkpoints; rows are not relabelled by this diagnostic script.",
    x = NULL,
    y = "Rows"
  ) +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))
save_plot(flow_plot, file.path(DIR_PLOTS_CLINICAL, "uti_not_uti_clinical_rule_flow.png"), width = 10, height = 5.8)

near_miss <- status_primary %>%
  filter(
    .data$UTI_Status == "Not_UTI",
    .data$Infection_Status_legacy == "UTI"
  ) %>%
  mutate(
    near_miss_reason = case_when(
      .data$culture_supports_uti %in% TRUE & .data$symptom_compatible_uti %in% FALSE ~
        "Culture supports UTI, symptom rule not met",
      .data$culture_supports_uti %in% TRUE & is.na(.data$symptom_compatible_uti) ~
        "Culture supports UTI, symptom rule unknown",
      TRUE ~ "Legacy UTI but primary UTI rule not met"
    )
  ) %>%
  arrange(.data$near_miss_reason, .data$Participant_id, .data$tp_lab) %>%
  select(any_of(c(
    "Participant_id", "Timepoint", "tp_lab", "Collection_Date", "Batch",
    "UTI_Label", "UTI_Status", "Not_UTI_subgroup", "Infection_Status_legacy",
    "near_miss_reason", "Urine_collection_method", "catheter_rule",
    "culture_supports_uti", "cfu_raw", "cfu_raw_parsed", "cfu_threshold_source",
    "symptom_compatible_uti", "symptom_rule_met", "local_urinary_symptom_any",
    "systemic_symptom_any", "dysuria_present", "urgency_present",
    "frequency_present", "incontinence_present", "pus_present",
    "flankpain_present", "suprapubic_pain_present", "fever_present",
    "rigors_present", "delirium_present", "UTI_classification_reason",
    "Episode_ID"
  )))
write_csv(near_miss, file.path(DIR_AUDIT, "uti_not_uti_near_miss_rows.csv"))

evidence_cols <- intersect(
  c(
    "culture_supports_uti", "symptom_compatible_uti",
    "local_urinary_symptom_any", "systemic_symptom_any", "dysuria_present",
    "urgency_present", "frequency_present", "incontinence_present",
    "pus_present", "flankpain_present", "suprapubic_pain_present",
    "fever_present", "rigors_present", "delirium_present"
  ),
  names(near_miss)
)

if (nrow(near_miss) > 0 && length(evidence_cols) > 0) {
  near_miss_heat <- near_miss %>%
    mutate(row_label = paste(.data$Participant_id, .data$tp_lab)) %>%
    select(row_label, near_miss_reason, all_of(evidence_cols)) %>%
    pivot_longer(all_of(evidence_cols), names_to = "evidence_field", values_to = "value") %>%
    mutate(
      evidence_value = case_when(
        is.na(.data$value) ~ "Unknown",
        as.character(.data$value) %in% c("TRUE", "1", "true", "Yes", "yes") ~ "Present/true",
        as.character(.data$value) %in% c("FALSE", "0", "false", "No", "no") ~ "Absent/false",
        TRUE ~ as.character(.data$value)
      ),
      row_label = factor(.data$row_label, levels = rev(unique(.data$row_label))),
      evidence_field = factor(.data$evidence_field, levels = evidence_cols)
    )

  heat_plot <- ggplot(near_miss_heat, aes(x = evidence_field, y = row_label, fill = evidence_value)) +
    geom_tile(colour = "white", linewidth = 0.25) +
    scale_fill_manual(
      values = c("Present/true" = "#0F766E", "Absent/false" = "#E5E7EB", "Unknown" = "#F59E0B"),
      na.value = "#F3F4F6"
    ) +
    labs(
      title = "Near-Miss Legacy UTI Rows Now Classified Not_UTI",
      subtitle = "Rows keep their current labels; the plot shows which rule component prevents primary UTI status.",
      x = NULL,
      y = "Participant / timepoint",
      fill = "Evidence"
    ) +
    theme_bw(base_size = 10) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.grid = element_blank())
} else {
  heat_plot <- empty_plot("Near-Miss Legacy UTI Rows Now Classified Not_UTI")
}
save_plot(heat_plot, file.path(DIR_PLOTS_CLINICAL, "uti_not_uti_near_miss_evidence_heatmap.png"), width = 12, height = 7)

# ------------------------------------------------------------------------------
# Clinical to VF/model attrition
# ------------------------------------------------------------------------------

clinical_without_vf <- status_primary %>%
  anti_join(vf %>% distinct(.data$Participant_id, .data$tp_lab), by = c("Participant_id", "tp_lab"))

attrition <- tibble::tribble(
  ~stage_order, ~stage, ~n, ~category, ~interpretation,
  1L, "Primary clinical rows", nrow(status_primary), "starting denominator",
  "Clinical primary denominator after manual primary exclusions.",
  2L, "Primary clinical UTI", sum(status_primary$UTI_Status == "UTI", na.rm = TRUE), "clinical subset",
  "Clinical rows meeting the primary culture plus symptom rule.",
  3L, "Primary clinical Not_UTI", sum(status_primary$UTI_Status == "Not_UTI", na.rm = TRUE), "clinical subset",
  "Clinical rows not meeting both primary UTI rule components.",
  4L, "Clinical rows without primary VF/model row", nrow(clinical_without_vf), "loss category",
  "Primary clinical rows absent from the VF/model denominator after genomics filters.",
  5L, "Manual primary exclusions", nrow(manual_exclusions), "loss category",
  "Rows retained for audit but excluded from primary analyses.",
  6L, "Quarantined failed/not-expected FASTA rows", nrow(quarantined_fastas), "loss category",
  "Rows removed from active genomics expected-denominator checks.",
  7L, "Primary VF/model rows", nrow(vf), "final denominator",
  "Rows used in primary VF/model analyses.",
  8L, "Model denominator rows", if (nrow(model_denom) > 0) nrow(model_denom) else NA_integer_, "final denominator",
  "Rows recorded by 14_genotype_phenotype_model.R model denominator output."
)
write_csv(attrition, file.path(DIR_VF, "uti_not_uti_attrition_summary.csv"))

attrition_plot <- attrition %>%
  mutate(stage = factor(.data$stage, levels = .data$stage)) %>%
  ggplot(aes(stage, n, fill = category)) +
  geom_col(width = 0.68, colour = "grey25", linewidth = 0.2) +
  geom_text(aes(label = if_else(is.na(n), "NA", as.character(n))), vjust = -0.35, size = 3.3) +
  scale_fill_manual(
    values = c(
      "starting denominator" = "#2563EB", "clinical subset" = "#0F766E",
      "loss category" = "#B45309", "final denominator" = "#7C3AED"
    )
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(
    title = "Clinical to VF/Model Denominator Context",
    subtitle = "Loss categories are listed explicitly so failed/quarantined and duplicate/unknown rows do not pollute active denominators.",
    x = NULL,
    y = "Rows",
    fill = NULL
  ) +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1), legend.position = "bottom")
save_plot(attrition_plot, file.path(DIR_PLOTS_VF, "uti_not_uti_denominator_waterfall.png"), width = 11, height = 6)

# ------------------------------------------------------------------------------
# Bootstrap supplementary endpoint effects
# ------------------------------------------------------------------------------

score_effects <- function(df, cols) {
  bind_rows(lapply(cols, function(col) {
    uti_vals <- df %>% filter(.data$UTI_Status == "UTI") %>% pull(all_of(col)) %>% safe_num()
    not_vals <- df %>% filter(.data$UTI_Status == "Not_UTI") %>% pull(all_of(col)) %>% safe_num()
    tibble(
      score = col,
      n_uti = sum(!is.na(uti_vals)),
      n_not_uti = sum(!is.na(not_vals)),
      mean_uti = mean(uti_vals, na.rm = TRUE),
      mean_not_uti = mean(not_vals, na.rm = TRUE),
      median_uti = median(uti_vals, na.rm = TRUE),
      median_not_uti = median(not_vals, na.rm = TRUE),
      mean_diff = mean_uti - mean_not_uti,
      median_diff = median_uti - median_not_uti
    )
  }))
}

observed_score_effects <- score_effects(score_small, score_cols)

bootstrap_participants <- unique(score_small$Participant_id)
bootstrap_B <- 1000L
boot_effects <- bind_rows(lapply(seq_len(bootstrap_B), function(i) {
  sampled <- sample(bootstrap_participants, length(bootstrap_participants), replace = TRUE)
  boot_df <- bind_rows(lapply(seq_along(sampled), function(j) {
    score_small %>%
      filter(.data$Participant_id == sampled[[j]]) %>%
      mutate(.bootstrap_participant_draw = j)
  }))
  score_effects(boot_df, score_cols) %>% mutate(bootstrap_id = i)
}))

bootstrap_effects <- boot_effects %>%
  select(bootstrap_id, score, mean_diff, median_diff) %>%
  pivot_longer(c(mean_diff, median_diff), names_to = "statistic", values_to = "bootstrap_effect") %>%
  group_by(.data$score, .data$statistic) %>%
  summarise(
    ci_lower = quantile(.data$bootstrap_effect, 0.025, na.rm = TRUE, names = FALSE),
    ci_upper = quantile(.data$bootstrap_effect, 0.975, na.rm = TRUE, names = FALSE),
    bootstrap_sd = sd(.data$bootstrap_effect, na.rm = TRUE),
    bootstrap_n = sum(!is.na(.data$bootstrap_effect)),
    .groups = "drop"
  ) %>%
  left_join(
    observed_score_effects %>%
      select(score, n_uti, n_not_uti, mean_diff, median_diff) %>%
      pivot_longer(c(mean_diff, median_diff), names_to = "statistic", values_to = "observed_effect"),
    by = c("score", "statistic")
  ) %>%
  mutate(
    score_label = score_label_lookup[.data$score],
    statistic = recode(.data$statistic, mean_diff = "mean difference", median_diff = "median difference"),
    interpretation = "Participant bootstrap CI for UTI - Not_UTI supplementary endpoint differences; descriptive because UTI n is sparse and residents can repeat."
  ) %>%
  arrange(.data$statistic, desc(abs(.data$observed_effect)))
write_csv(bootstrap_effects, file.path(DIR_VF, "uti_not_uti_bootstrap_effects.csv"))

forest_data <- bootstrap_effects %>%
  filter(.data$statistic == "median difference") %>%
  arrange(desc(abs(.data$observed_effect))) %>%
  slice_head(n = 12) %>%
  mutate(score_label = factor(.data$score_label, levels = rev(.data$score_label)))

forest_plot <- ggplot(forest_data, aes(x = observed_effect, y = score_label)) +
  geom_vline(xintercept = 0, colour = "grey55", linewidth = 0.4) +
  geom_errorbarh(aes(xmin = ci_lower, xmax = ci_upper), height = 0.22, colour = "#2563EB") +
  geom_point(size = 2.2, colour = "#0F766E") +
  labs(
    title = "Participant-Bootstrap VF Endpoint Effects",
    subtitle = "Effect is UTI - Not_UTI. Intervals show bootstrap uncertainty, not confirmatory inference.",
    x = "Median difference",
    y = NULL
  ) +
  theme_bw(base_size = 11)
save_plot(forest_plot, file.path(DIR_PLOTS_VF, "uti_not_uti_bootstrap_effect_forest.png"), width = 8.5, height = 6)

# ------------------------------------------------------------------------------
# Fisher exact screening and leave-one-UTI-out stability
# ------------------------------------------------------------------------------

binary_feature_effect <- function(df, feature) {
  value <- suppressWarnings(as.integer(df[[feature]]))
  value[is.na(value)] <- 0L
  status <- df$UTI_Status
  a <- sum(status == "UTI" & value == 1L, na.rm = TRUE)
  b <- sum(status == "UTI" & value == 0L, na.rm = TRUE)
  c <- sum(status == "Not_UTI" & value == 1L, na.rm = TRUE)
  d <- sum(status == "Not_UTI" & value == 0L, na.rm = TRUE)
  mat <- matrix(c(a, b, c, d), nrow = 2, byrow = TRUE)
  fisher <- tryCatch(fisher.test(mat), error = function(e) NULL)
  or_ha <- ((a + 0.5) / (b + 0.5)) / ((c + 0.5) / (d + 0.5))
  tibble(
    feature = feature,
    n_uti_present = a,
    n_uti_absent = b,
    n_not_uti_present = c,
    n_not_uti_absent = d,
    prevalence_uti = a / max(1, a + b),
    prevalence_not_uti = c / max(1, c + d),
    prevalence_diff = prevalence_uti - prevalence_not_uti,
    log2_or_ha = log2(or_ha),
    fisher_p = if (is.null(fisher)) NA_real_ else fisher$p.value,
    fisher_or = if (is.null(fisher)) NA_real_ else unname(fisher$estimate)
  )
}

module_cols <- names(score_df)[str_detect(names(score_df), "^mod_.*_present$")]
module_cols <- module_cols[vapply(score_df[module_cols], is.numeric, logical(1))]
gene_cols <- canonical_vf_gene_cols(names(vf), required = FALSE)

fisher_existing <- safe_read(file.path(DIR_VF, "vf_fisher_exploratory.csv"))
top_existing_genes <- if (nrow(fisher_existing) > 0 && "gene" %in% names(fisher_existing)) {
  fisher_existing %>%
    filter(.data$gene %in% gene_cols) %>%
    arrange(.data$p_value) %>%
    pull(.data$gene) %>%
    unique() %>%
    head(25)
} else {
  character()
}

if (length(top_existing_genes) < 25 && length(gene_cols) > 0) {
  gene_prevalence <- vf %>%
    summarise(across(all_of(gene_cols), ~ mean(safe_num(.x) > 0, na.rm = TRUE))) %>%
    pivot_longer(everything(), names_to = "gene", values_to = "prevalence") %>%
    arrange(desc(.data$prevalence)) %>%
    pull(.data$gene)
  top_existing_genes <- unique(c(top_existing_genes, head(gene_prevalence, 25)))
}

feature_cols <- unique(c(head(top_existing_genes, 25), module_cols))
feature_cols <- intersect(feature_cols, names(score_df))
feature_small <- score_df %>%
  select(any_of(c("Participant_id", "tp_lab", "Episode_ID", "UTI_Status", feature_cols))) %>%
  mutate(across(all_of(feature_cols), ~ as.integer(safe_num(.x) > 0)))

feature_fisher <- if (length(feature_cols) > 0) {
  bind_rows(lapply(feature_cols, function(feature) {
    binary_feature_effect(feature_small, feature) %>%
      mutate(
        feature_type = if_else(str_detect(feature, "^mod_"), "module", "gene"),
        note = "Exploratory Fisher exact test; repeated residents and lineage are not modelled."
      )
  })) %>%
    arrange(.data$fisher_p, desc(abs(.data$prevalence_diff)))
} else {
  tibble()
}
write_csv(feature_fisher, file.path(DIR_VF, "uti_not_uti_feature_fisher_exploratory.csv"))

uti_rows <- score_small %>%
  filter(.data$UTI_Status == "UTI") %>%
  mutate(uti_row_id = paste(.data$Participant_id, .data$tp_lab, sep = "__")) %>%
  select(Participant_id, tp_lab, uti_row_id)

leave_one_score <- bind_rows(lapply(score_cols, function(score) {
  full <- observed_score_effects %>% filter(.data$score == !!score)
  by_removed <- bind_rows(lapply(seq_len(nrow(uti_rows)), function(i) {
    removed <- uti_rows[i, , drop = FALSE]
    tmp <- score_small %>%
      filter(!(.data$Participant_id == removed$Participant_id & .data$tp_lab == removed$tp_lab))
    eff <- score_effects(tmp, score)
    tibble(
      removed_uti = removed$uti_row_id,
      statistic = c("mean difference", "median difference"),
      effect_without_one_uti = c(eff$mean_diff, eff$median_diff)
    )
  }))
  bind_rows(
    tibble(
      feature_type = "score",
      feature = score,
      statistic = "mean difference",
      full_effect = full$mean_diff,
      n_uti = full$n_uti,
      n_not_uti = full$n_not_uti
    ),
    tibble(
      feature_type = "score",
      feature = score,
      statistic = "median difference",
      full_effect = full$median_diff,
      n_uti = full$n_uti,
      n_not_uti = full$n_not_uti
    )
  ) %>%
    left_join(
      by_removed %>%
        group_by(.data$statistic) %>%
        summarise(
          min_without_one_uti = min(.data$effect_without_one_uti, na.rm = TRUE),
          max_without_one_uti = max(.data$effect_without_one_uti, na.rm = TRUE),
          removed_uti_most_lower = .data$removed_uti[which.min(.data$effect_without_one_uti)][1],
          removed_uti_most_upper = .data$removed_uti[which.max(.data$effect_without_one_uti)][1],
          .groups = "drop"
        ),
      by = "statistic"
    )
}))

top_feature_cols <- if (nrow(feature_fisher) > 0) {
  feature_fisher %>%
    arrange(.data$fisher_p, desc(abs(.data$prevalence_diff))) %>%
    pull(.data$feature) %>%
    unique() %>%
    head(30)
} else {
  character()
}

leave_one_feature <- if (length(top_feature_cols) > 0) {
  bind_rows(lapply(top_feature_cols, function(feature) {
    full <- binary_feature_effect(feature_small, feature)
    by_removed <- bind_rows(lapply(seq_len(nrow(uti_rows)), function(i) {
      removed <- uti_rows[i, , drop = FALSE]
      tmp <- feature_small %>%
        filter(!(.data$Participant_id == removed$Participant_id & .data$tp_lab == removed$tp_lab))
      eff <- binary_feature_effect(tmp, feature)
      tibble(
        removed_uti = removed$uti_row_id,
        statistic = c("prevalence difference", "log2 odds ratio"),
        effect_without_one_uti = c(eff$prevalence_diff, eff$log2_or_ha)
      )
    }))
    tibble(
      feature_type = if_else(str_detect(feature, "^mod_"), "module", "gene"),
      feature = feature,
      statistic = c("prevalence difference", "log2 odds ratio"),
      full_effect = c(full$prevalence_diff, full$log2_or_ha),
      n_uti = full$n_uti_present + full$n_uti_absent,
      n_not_uti = full$n_not_uti_present + full$n_not_uti_absent
    ) %>%
      left_join(
        by_removed %>%
          group_by(.data$statistic) %>%
          summarise(
            min_without_one_uti = min(.data$effect_without_one_uti, na.rm = TRUE),
            max_without_one_uti = max(.data$effect_without_one_uti, na.rm = TRUE),
            removed_uti_most_lower = .data$removed_uti[which.min(.data$effect_without_one_uti)][1],
            removed_uti_most_upper = .data$removed_uti[which.max(.data$effect_without_one_uti)][1],
            .groups = "drop"
          ),
        by = "statistic"
      )
  }))
} else {
  tibble()
}

leave_one <- bind_rows(leave_one_score, leave_one_feature) %>%
  mutate(
    sensitivity_range = .data$max_without_one_uti - .data$min_without_one_uti,
    direction_flip = sign(.data$min_without_one_uti) != sign(.data$max_without_one_uti),
    sparse_uti_warning = .data$n_uti < 20,
    interpretation = "Leave-one-UTI-out diagnostic; large ranges or direction flips mean a result may be driven by one UTI row."
  ) %>%
  arrange(desc(.data$direction_flip), desc(.data$sensitivity_range))
write_csv(leave_one, file.path(DIR_VF, "uti_not_uti_leave_one_uti_out.csv"))

stability_plot_data <- leave_one %>%
  filter(
    (.data$feature_type == "score" & .data$statistic == "median difference") |
      (.data$feature_type != "score" & .data$statistic == "prevalence difference")
  ) %>%
  arrange(desc(.data$sensitivity_range)) %>%
  slice_head(n = 20) %>%
  mutate(feature = factor(.data$feature, levels = rev(.data$feature)))

stability_plot <- if (nrow(stability_plot_data) > 0) {
  ggplot(stability_plot_data, aes(y = feature, x = full_effect, colour = feature_type)) +
    geom_vline(xintercept = 0, colour = "grey60", linewidth = 0.35) +
    geom_errorbarh(aes(xmin = min_without_one_uti, xmax = max_without_one_uti), height = 0.22) +
    geom_point(size = 2.2) +
    facet_wrap(~statistic, scales = "free_x") +
    labs(
      title = "Leave-One-UTI-Out Stability",
      subtitle = "Error bars show the range after removing each UTI row one at a time.",
      x = "Effect after UTI removal range",
      y = NULL,
      colour = "Feature type"
    ) +
    theme_bw(base_size = 11)
} else {
  empty_plot("Leave-One-UTI-Out Stability")
}
save_plot(stability_plot, file.path(DIR_PLOTS_VF, "uti_not_uti_leave_one_uti_out_stability.png"), width = 10, height = 7)

# ------------------------------------------------------------------------------
# Paired resident and transition diagnostics
# ------------------------------------------------------------------------------

paired_rows <- bind_rows(lapply(score_cols, function(score) {
  score_small %>%
    group_by(.data$Participant_id, .data$UTI_Status) %>%
    summarise(value = mean(.data[[score]], na.rm = TRUE), .groups = "drop") %>%
    filter(.data$UTI_Status %in% c("UTI", "Not_UTI")) %>%
    pivot_wider(names_from = "UTI_Status", values_from = "value") %>%
    filter(!is.na(.data$UTI), !is.na(.data$Not_UTI)) %>%
    mutate(
      score = score,
      delta_uti_minus_not_uti = .data$UTI - .data$Not_UTI
    ) %>%
    select(Participant_id, score, Not_UTI, UTI, delta_uti_minus_not_uti)
}))
write_csv(paired_rows, file.path(DIR_VF, "uti_not_uti_paired_participant_deltas.csv"))

paired_tests <- paired_rows %>%
  group_by(.data$score) %>%
  summarise(
    n_paired_participants = n(),
    median_delta = median(.data$delta_uti_minus_not_uti, na.rm = TRUE),
    mean_delta = mean(.data$delta_uti_minus_not_uti, na.rm = TRUE),
    n_positive = sum(.data$delta_uti_minus_not_uti > 0, na.rm = TRUE),
    n_negative = sum(.data$delta_uti_minus_not_uti < 0, na.rm = TRUE),
    n_nonzero = n_positive + n_negative,
    sign_test_p = if_else(
      n_nonzero > 0,
      as.numeric(binom.test(n_positive, n_nonzero, p = 0.5)$p.value),
      NA_real_
    ),
    signed_rank_p = if_else(
      n_paired_participants >= 3,
      exact_p(wilcox.test(.data$delta_uti_minus_not_uti, mu = 0, exact = FALSE)$p.value),
      NA_real_
    ),
    interpretation = "Participant-paired descriptive comparison for residents with both UTI and Not_UTI VF rows.",
    .groups = "drop"
  ) %>%
  arrange(.data$sign_test_p)
write_csv(paired_tests, file.path(DIR_VF, "uti_not_uti_paired_participant_tests.csv"))

paired_plot_scores <- intersect(
  c("expec_marker_count", "upec_system_count", "upec_system_fraction", "total_vf_count_curated"),
  unique(paired_rows$score)
)
if (length(paired_plot_scores) == 0) paired_plot_scores <- head(unique(paired_rows$score), 4)

paired_plot <- if (nrow(paired_rows) > 0 && length(paired_plot_scores) > 0) {
  paired_rows %>%
    filter(.data$score %in% paired_plot_scores) %>%
    pivot_longer(c(Not_UTI, UTI), names_to = "status", values_to = "value") %>%
    mutate(
      status = factor(.data$status, levels = c("Not_UTI", "UTI")),
      score_label = factor(score_label_lookup[.data$score],
                           levels = score_label_lookup[paired_plot_scores])
    ) %>%
    ggplot(aes(x = status, y = value, group = Participant_id)) +
    geom_line(colour = "grey50", alpha = 0.65) +
    geom_point(aes(colour = status), size = 2) +
    facet_wrap(~score_label, scales = "free_y") +
    scale_colour_manual(values = c(Not_UTI = "#2563EB", UTI = "#B45309")) +
    labs(
      title = "Within-Resident UTI / Not_UTI VF Endpoint Pairs",
      subtitle = "Only residents with both statuses are shown; repeated rows are averaged within status.",
      x = NULL,
      y = "Participant-level mean endpoint value",
      colour = NULL
    ) +
    theme_bw(base_size = 11) +
    theme(legend.position = "bottom")
} else {
  empty_plot("Within-Resident UTI / Not_UTI VF Endpoint Pairs")
}
save_plot(paired_plot, file.path(DIR_PLOTS_VF, "uti_not_uti_paired_participant_slopeplot.png"), width = 10, height = 6)

transition_tests <- if (nrow(transition_score_changes) > 0 &&
                        all(c("score_name", "delta") %in% names(transition_score_changes))) {
  transition_score_changes %>%
    mutate(delta = safe_num(.data$delta)) %>%
    group_by(.data$score_name) %>%
    summarise(
      n_transitions = sum(!is.na(.data$delta)),
      median_delta = median(.data$delta, na.rm = TRUE),
      mean_delta = mean(.data$delta, na.rm = TRUE),
      q25_delta = quantile(.data$delta, 0.25, na.rm = TRUE, names = FALSE),
      q75_delta = quantile(.data$delta, 0.75, na.rm = TRUE, names = FALSE),
      n_positive = sum(.data$delta > 0, na.rm = TRUE),
      n_negative = sum(.data$delta < 0, na.rm = TRUE),
      n_nonzero = n_positive + n_negative,
      sign_test_p = if_else(
        n_nonzero > 0,
        as.numeric(binom.test(n_positive, n_nonzero, p = 0.5)$p.value),
        NA_real_
      ),
      signed_rank_p = if_else(
        n_transitions >= 3,
        exact_p(wilcox.test(.data$delta, mu = 0, exact = FALSE)$p.value),
        NA_real_
      ),
      interpretation = "Transition-level score-change diagnostic; not a relabelling rule and not independent when residents repeat.",
      .groups = "drop"
    ) %>%
    arrange(.data$sign_test_p)
} else {
  tibble()
}
write_csv(transition_tests, file.path(DIR_VF, "uti_not_uti_transition_score_tests.csv"))

transition_plot <- if (nrow(transition_tests) > 0) {
  transition_tests %>%
    arrange(desc(abs(.data$median_delta))) %>%
    slice_head(n = 14) %>%
    mutate(score_name = factor(.data$score_name, levels = rev(.data$score_name))) %>%
    ggplot(aes(x = median_delta, y = score_name)) +
    geom_vline(xintercept = 0, colour = "grey60", linewidth = 0.35) +
    geom_errorbarh(aes(xmin = q25_delta, xmax = q75_delta), height = 0.22, colour = "#2563EB") +
    geom_point(colour = "#B45309", size = 2.2) +
    labs(
      title = "Not_UTI to UTI Transition Score Changes",
      subtitle = "Point is median delta; bar is the interquartile range across transition cases.",
      x = "Score delta",
      y = NULL
    ) +
    theme_bw(base_size = 11)
} else {
  empty_plot("Not_UTI to UTI Transition Score Changes")
}
save_plot(transition_plot, file.path(DIR_PLOTS_VF, "uti_not_uti_transition_delta_forest.png"), width = 9, height = 6)

# ------------------------------------------------------------------------------
# Duplicate culture QC: 31036 UTI-1 vs UTI-2
# ------------------------------------------------------------------------------

status_31036 <- status_all %>%
  mutate(participant_no_zero = participant_without_leading_zero(.data$Participant_id)) %>%
  filter(.data$participant_no_zero == "31036", .data$tp_lab %in% c("UTI-1", "UTI-2")) %>%
  select(any_of(c(
    "Participant_id", "Timepoint", "tp_lab", "Collection_Date", "Batch",
    "UTI_Label", "UTI_Status", "Not_UTI_subgroup", "Infection_Status_legacy",
    "analysis_include_primary", "analysis_exclusion_reason", "duplicate_role",
    "duplicate_of_participant_id", "duplicate_of_tp_lab",
    "allow_secondary_duplicate_qc", "duplicate_use_note", "manual_curation_note",
    "culture_supports_uti", "symptom_compatible_uti", "Episode_ID"
  )))

score_31036 <- score_raw %>%
  { if (nrow(.) > 0) . else tibble() } %>%
  mutate(
    Participant_id = as.character(.data$Participant_id),
    tp_lab = as.character(.data$tp_lab),
    participant_no_zero = participant_without_leading_zero(.data$Participant_id)
  ) %>%
  filter(.data$participant_no_zero == "31036", .data$tp_lab %in% c("UTI-1", "UTI-2")) %>%
  select(any_of(c(
    "Participant_id", "tp_lab", "ST", "vf_count_total",
    "total_vf_count_curated", "expec_marker_count",
    "upec_system_count", "upec_system_fraction"
  )))

duplicate_qc <- status_31036 %>%
  full_join(score_31036, by = c("Participant_id", "tp_lab")) %>%
  mutate(
    duplicate_qc_role = case_when(
      .data$tp_lab == "UTI-1" ~ "primary duplicate representative",
      .data$tp_lab == "UTI-2" ~ "secondary duplicate, excluded from primary",
      TRUE ~ "not targeted"
    ),
    duplicate_qc_use = case_when(
      .data$tp_lab == "UTI-2" ~
        "Use only in explicit duplicate-culture or longitudinal QC, with reuse stated.",
      .data$tp_lab == "UTI-1" ~
        "Use as the primary representative if one duplicate UTI sample is required.",
      TRUE ~ NA_character_
    ),
    present_in_primary_vf_score_table = !is.na(.data$vf_count_total) | !is.na(.data$total_vf_count_curated)
  ) %>%
  arrange(.data$tp_lab)
write_csv(duplicate_qc, file.path(DIR_AUDIT, "duplicate_culture_qc_31036.csv"))

dup_plot <- if (nrow(duplicate_qc) > 0) {
  duplicate_qc %>%
    mutate(
      score_value = coalesce(safe_num(.data$total_vf_count_curated), safe_num(.data$vf_count_total), 0),
      label = paste0(
        .data$UTI_Status,
        "\nprimary=", .data$analysis_include_primary,
        "\nVF endpoint=", if_else(.data$present_in_primary_vf_score_table, "yes", "no")
      )
    ) %>%
    ggplot(aes(x = tp_lab, y = score_value, fill = duplicate_qc_role)) +
    geom_col(width = 0.58, colour = "grey25", linewidth = 0.2) +
    geom_text(aes(label = label), vjust = -0.15, size = 3.3, lineheight = 0.9) +
    scale_fill_manual(values = c(
      "primary duplicate representative" = "#0F766E",
      "secondary duplicate, excluded from primary" = "#B45309",
      "not targeted" = "#6B7280"
    )) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.3))) +
    labs(
      title = "Duplicate Culture QC: 31036 UTI-1 vs UTI-2",
      subtitle = "UTI-2 is retained only for explicit duplicate QC and is not a primary analysis row.",
      x = NULL,
      y = "Curated VF burden if present",
      fill = NULL
    ) +
    theme_bw(base_size = 11) +
    theme(legend.position = "bottom")
} else {
  empty_plot("Duplicate Culture QC: 31036 UTI-1 vs UTI-2")
}
save_plot(dup_plot, file.path(DIR_PLOTS_VF, "duplicate_culture_qc_31036.png"), width = 8.5, height = 5.5)

# ------------------------------------------------------------------------------
# Sparse-count power/precision context
# ------------------------------------------------------------------------------

vf_n_uti <- sum(vf$UTI_Status == "UTI", na.rm = TRUE)
vf_n_not_uti <- sum(vf$UTI_Status == "Not_UTI", na.rm = TRUE)

power_precision <- expand_grid(
  baseline_not_uti_prevalence = seq(0.05, 0.95, by = 0.05),
  odds_ratio = c(0.20, 0.33, 0.50, 0.75, 1.00, 1.50, 2.00, 3.00, 5.00, 8.00)
) %>%
  mutate(
    implied_uti_prevalence = odds_ratio * baseline_not_uti_prevalence /
      (1 - baseline_not_uti_prevalence + odds_ratio * baseline_not_uti_prevalence),
    not_uti_present_expected = round(vf_n_not_uti * baseline_not_uti_prevalence),
    not_uti_absent_expected = vf_n_not_uti - .data$not_uti_present_expected,
    uti_present_expected = round(vf_n_uti * .data$implied_uti_prevalence),
    uti_absent_expected = vf_n_uti - .data$uti_present_expected,
    fisher_p_expected = pmap_dbl(
      list(.data$uti_present_expected, .data$uti_absent_expected,
           .data$not_uti_present_expected, .data$not_uti_absent_expected),
      function(a, b, c, d) {
        fisher.test(matrix(c(a, b, c, d), nrow = 2, byrow = TRUE))$p.value
      }
    ),
    detectable_at_alpha_0_05 = .data$fisher_p_expected < 0.05,
    interpretation = "Approximate sparse-count context using expected 2x2 counts, not a formal prospective power analysis."
  )
write_csv(power_precision, file.path(DIR_VF, "uti_not_uti_power_precision_context.csv"))

power_plot <- power_precision %>%
  mutate(odds_ratio = factor(.data$odds_ratio, levels = sort(unique(.data$odds_ratio)))) %>%
  ggplot(aes(x = baseline_not_uti_prevalence, y = odds_ratio, fill = -log10(pmax(fisher_p_expected, 1e-12)))) +
  geom_tile(colour = "white", linewidth = 0.25) +
  geom_text(aes(label = if_else(detectable_at_alpha_0_05, "*", "")), colour = "white", size = 4) +
  scale_fill_gradient(low = "#E5E7EB", high = "#7C2D12", name = "-log10 p") +
  scale_x_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(
    title = "Sparse-UTI Precision Context",
    subtitle = sprintf("Expected Fisher-test behaviour with n=%d UTI and n=%d Not_UTI. Stars mark p < 0.05.", vf_n_uti, vf_n_not_uti),
    x = "Assumed Not_UTI feature prevalence",
    y = "Assumed odds ratio"
  ) +
  theme_bw(base_size = 11) +
  theme(panel.grid = element_blank())
save_plot(power_plot, file.path(DIR_PLOTS_VF, "uti_not_uti_sparse_power_precision.png"), width = 9, height = 5.8)

# ------------------------------------------------------------------------------
# Interpretation and summary tables
# ------------------------------------------------------------------------------

interpretation_table <- tibble::tribble(
  ~diagnostic, ~evidence_type, ~primary_use, ~what_it_can_support, ~what_it_cannot_support,
  "Clinical rule decision flow", "descriptive", "denominator explanation",
  "Explains why rows are UTI, Not_UTI, excluded, or VF/model eligible.",
  "Does not change clinical labels or infer biology.",
  "Near-miss evidence heatmap", "descriptive", "rule audit",
  "Identifies culture-supported rows blocked by symptom incompatibility or unknown symptom rule.",
  "Cannot promote near-miss rows to UTI without upstream symptom/collection correction.",
  "Bootstrap VF endpoint effects", "exploratory participant bootstrap", "effect size and uncertainty",
  "Shows UTI - Not_UTI supplementary endpoint differences with resident-aware resampling.",
  "Does not provide confirmatory p-values or remove lineage/repeated-measure confounding.",
  "Fisher exact feature screen", "exploratory exact test", "feature prioritisation",
  "Highlights genes/modules with large sparse-count prevalence differences.",
  "Does not model repeated residents, ST, or multiple-testing robustness by itself.",
  "Leave-one-UTI-out stability", "sensitivity diagnostic", "sparse UTI robustness",
  "Shows whether endpoint/gene/module effects depend on a single UTI row.",
  "Does not prove a stable association if all UTIs are sparse or lineage-clustered.",
  "Paired resident endpoint deltas", "paired descriptive test", "within-resident contrast",
  "Uses residents who have both UTI and Not_UTI rows to reduce between-resident confounding.",
  "Does not represent residents who have only one status.",
  "Transition endpoint-change tests", "transition diagnostic", "longitudinal context",
    "Summarises supplementary endpoint changes around Not_UTI to UTI transitions.",
  "Does not relabel transition rows or prove causality.",
  "Duplicate culture QC 31036", "quarantine/QC", "duplicate handling traceability",
  "Documents that UTI-1 is primary and UTI-2 is secondary duplicate QC only.",
  "Does not include UTI-2 in primary denominator.",
  "Sparse-count precision context", "approximate context", "interpretation calibration",
  "Shows why many feature effects are unstable with 17 UTI VF rows.",
  "Not a replacement for a formal sample-size or prospective power analysis."
)
write_csv(interpretation_table, file.path(DIR_VF, "uti_not_uti_test_interpretation_table.csv"))

diagnostic_summary <- tibble::tribble(
  ~metric, ~value, ~source, ~interpretation,
  "primary_clinical_total", nrow(status_primary), FILE_STATUS_MAP,
  "Clinical primary denominator after manual exclusions.",
  "primary_clinical_uti", sum(status_primary$UTI_Status == "UTI", na.rm = TRUE), FILE_STATUS_MAP,
  "Clinical rows meeting culture support plus symptom compatibility.",
  "primary_clinical_not_uti", sum(status_primary$UTI_Status == "Not_UTI", na.rm = TRUE), FILE_STATUS_MAP,
  "Clinical rows not meeting both primary UTI rule components.",
  "primary_vf_model_total", nrow(vf), FILE_VF_READY,
  "VF/model denominator after primary and genomics filters.",
  "primary_vf_model_uti", vf_n_uti, FILE_VF_READY,
  "Sequenced VF/model rows with primary UTI status.",
  "primary_vf_model_not_uti", vf_n_not_uti, FILE_VF_READY,
  "Sequenced VF/model rows with primary Not_UTI status.",
  "primary_vf_model_missing_status", sum(is.na(vf$UTI_Status)), FILE_VF_READY,
  "Must remain zero for the current cleaned inputs.",
  "near_miss_rows", nrow(near_miss), file.path(DIR_AUDIT, "uti_not_uti_near_miss_rows.csv"),
  "Legacy UTI/culture-supported rows now primary Not_UTI.",
  "manual_primary_exclusions", nrow(manual_exclusions), file.path(DIR_AUDIT, "primary_clinical_manual_exclusions.csv"),
  "Unknown participant and duplicate culture rows excluded from primary analysis.",
  "quarantined_failed_fastas", nrow(quarantined_fastas), file.path(DIR_AUDIT, "quarantined_failed_or_not_expected_fastas.csv"),
  "Failed/not-expected FASTA rows outside active genomics expected denominator.",
  "bootstrap_scores_tested", length(score_cols), file.path(DIR_VF, "uti_not_uti_bootstrap_effects.csv"),
  "Numeric VF burden/score columns with participant-bootstrap diagnostics.",
  "feature_fisher_tests", nrow(feature_fisher), file.path(DIR_VF, "uti_not_uti_feature_fisher_exploratory.csv"),
  "Exploratory gene/module Fisher exact screens.",
  "leave_one_diagnostics", nrow(leave_one), file.path(DIR_VF, "uti_not_uti_leave_one_uti_out.csv"),
  "Score/gene/module stability rows after leave-one-UTI-out analysis.",
  "paired_participant_tests", nrow(paired_tests), file.path(DIR_VF, "uti_not_uti_paired_participant_tests.csv"),
  "Scores with resident-paired UTI vs Not_UTI deltas.",
  "transition_score_tests", nrow(transition_tests), file.path(DIR_VF, "uti_not_uti_transition_score_tests.csv"),
  "Transition score-change sign/signed-rank summaries."
)
write_csv(diagnostic_summary, file.path(DIR_VF, "uti_not_uti_diagnostic_summary.csv"))

figure_metadata <- tibble::tribble(
  ~figure_id, ~file_path, ~evidence_type, ~interpretation_limitations,
  "uti_not_uti_clinical_rule_flow", file.path(DIR_PLOTS_CLINICAL, "uti_not_uti_clinical_rule_flow.png"),
  "descriptive", "Denominator/rule explanation only.",
  "uti_not_uti_denominator_waterfall", file.path(DIR_PLOTS_VF, "uti_not_uti_denominator_waterfall.png"),
  "descriptive", "Attrition context, not a statistical test.",
  "uti_not_uti_near_miss_evidence_heatmap", file.path(DIR_PLOTS_CLINICAL, "uti_not_uti_near_miss_evidence_heatmap.png"),
  "descriptive", "Near-miss rows are not relabelled.",
  "uti_not_uti_bootstrap_effect_forest", file.path(DIR_PLOTS_VF, "uti_not_uti_bootstrap_effect_forest.png"),
  "exploratory participant bootstrap", "UTI n is sparse; intervals are descriptive.",
  "uti_not_uti_leave_one_uti_out_stability", file.path(DIR_PLOTS_VF, "uti_not_uti_leave_one_uti_out_stability.png"),
  "sensitivity diagnostic", "Large instability implies single-UTI dependence.",
  "uti_not_uti_paired_participant_slopeplot", file.path(DIR_PLOTS_VF, "uti_not_uti_paired_participant_slopeplot.png"),
  "paired descriptive", "Only includes residents with both statuses.",
  "uti_not_uti_transition_delta_forest", file.path(DIR_PLOTS_VF, "uti_not_uti_transition_delta_forest.png"),
  "transition diagnostic", "Transition rows are not independent causal events.",
  "uti_not_uti_sparse_power_precision", file.path(DIR_PLOTS_VF, "uti_not_uti_sparse_power_precision.png"),
  "approximate context", "Expected-count context, not formal power.",
  "duplicate_culture_qc_31036", file.path(DIR_PLOTS_VF, "duplicate_culture_qc_31036.png"),
  "duplicate QC", "UTI-2 is excluded from primary denominators."
) %>%
  mutate(output_exists = file.exists(.data$file_path))
write_csv(figure_metadata, file.path(DIR_VF, "uti_not_uti_diagnostic_figure_metadata.csv"))

msg("Wrote UTI/Not_UTI diagnostic outputs.")
msg("Clinical primary: %d UTI, %d Not_UTI.", expected_clinical[["UTI"]], expected_clinical[["Not_UTI"]])
msg("VF/model primary: %d UTI, %d Not_UTI, %d missing status.",
    expected_vf[["UTI"]], expected_vf[["Not_UTI"]], expected_vf[["<missing>"]])

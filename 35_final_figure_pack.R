#!/usr/bin/env Rscript

# ==============================================================================
# 35_final_figure_pack.R
# ==============================================================================
#
# Goal:
#   Render a manuscript-ready final figure pack from the current validated
#   primary UTI/Not_UTI, mechanism, and robustness outputs.
#
# Guardrails:
#   - This script does not reclassify episodes or rerun upstream analyses.
#   - It does not use legacy ASB-vs-UTI or archived OLD generated outputs.
#   - Figures are descriptive or sensitivity-focused, not confirmatory.
# ==============================================================================

source("00_config.R")
source("R/plot_helpers.R")

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
  library(scales)
  library(patchwork)
  library(ggrepel)
})

msg("Starting 35_final_figure_pack.R")

DIR_FINAL_RESULTS <- file.path(DIR_RESULTS, "final_figures")
DIR_FINAL_PLOTS <- file.path(DIR_PLOTS, "final_figures")
DIR_MECHANISM <- file.path(DIR_RESULTS, "mechanism")
DIR_ROBUSTNESS <- file.path(DIR_RESULTS, "robustness")
DIR_STAT <- file.path(DIR_RESULTS, "statistical_sensitivity")
DIR_AUDIT <- file.path(DIR_RESULTS, "audit")
ensure_dir(DIR_FINAL_RESULTS)
ensure_dir(DIR_FINAL_PLOTS)

required_paths <- c(
  status_map = FILE_STATUS_MAP,
  vf_ready = FILE_VF_READY,
  casebook = file.path(DIR_MECHANISM, "not_uti_to_uti_casebook.csv"),
  mechanism_summary = file.path(DIR_MECHANISM, "transition_mechanism_summary.csv"),
  mechanism_validation = file.path(DIR_MECHANISM, "mechanism_validation_checks.csv"),
  denominator = file.path(DIR_ROBUSTNESS, "denominator_robustness_summary.csv"),
  robustness_validation = file.path(DIR_ROBUSTNESS, "robustness_validation_checks.csv"),
  score_contrasts = file.path(DIR_ROBUSTNESS, "near_miss_expanded_score_contrasts.csv"),
  model_stability = file.path(DIR_ROBUSTNESS, "model_stability_summary.csv"),
  leave_one_summary = file.path(DIR_ROBUSTNESS, "leave_one_uti_sensitivity_summary.csv"),
  bootstrap = file.path(DIR_ROBUSTNESS, "bootstrap_score_robustness.csv"),
  near_miss = file.path(DIR_AUDIT, "uti_not_uti_near_miss_rows.csv"),
  power_precision = file.path(DIR_VF, "uti_not_uti_power_precision_context.csv"),
  leave_one_detail = file.path(DIR_VF, "uti_not_uti_leave_one_uti_out.csv"),
  variants = file.path(DIR_RESULTS, "longitudinal", "variant_annotation_detailed.csv"),
  stat_validation = file.path(DIR_STAT, "statistical_sensitivity_validation_checks.csv"),
  stat_score_values = file.path(DIR_STAT, "participant_collapsed_score_values.csv"),
  stat_module_matrix = file.path(DIR_STAT, "not_uti_to_uti_module_change_matrix.csv"),
  stat_transition_module = file.path(DIR_STAT, "transition_module_gain_loss_enrichment.csv"),
  stat_pcoa = file.path(DIR_VF, "vf_pcoa_jaccard_coordinates.csv")
)

missing_paths <- required_paths[!file.exists(required_paths)]
if (length(missing_paths) > 0) {
  stop("Missing required final-figure input(s): ", paste(missing_paths, collapse = ", "))
}
if (any(str_detect(required_paths, regex("_OLD|old_asb|legacy", ignore_case = TRUE)))) {
  stop("Final figure inputs include a legacy or OLD path; refusing to render.")
}

read_current <- function(path) read_csv(path, show_col_types = FALSE)

mechanism_validation <- read_current(required_paths[["mechanism_validation"]])
robustness_validation <- read_current(required_paths[["robustness_validation"]])
stat_validation <- read_current(required_paths[["stat_validation"]])
if (any(mechanism_validation$status != "PASS")) {
  stop("Mechanism validation includes a non-PASS check; rerun or repair the mechanism add-on first.")
}
if (any(robustness_validation$status != "PASS")) {
  stop("Robustness validation includes a non-PASS check; rerun or repair the robustness add-on first.")
}
if (any(stat_validation$status != "PASS")) {
  stop("Statistical sensitivity validation includes a non-PASS check; rerun or repair script 36 first.")
}

status <- read_current(required_paths[["status_map"]]) %>%
  mutate(
    Participant_id = as.character(.data$Participant_id),
    tp_lab = as.character(.data$tp_lab)
  ) %>%
  filter(.data$analysis_include_primary %in% TRUE, .data$UTI_Status %in% c("UTI", "Not_UTI"))

vf_ready <- read_current(required_paths[["vf_ready"]]) %>%
  mutate(
    Participant_id = as.character(.data$Participant_id),
    tp_lab = as.character(.data$tp_lab)
  ) %>%
  filter(.data$UTI_Status %in% c("UTI", "Not_UTI"))

casebook <- read_current(required_paths[["casebook"]]) %>%
  mutate(Participant_id = as.character(.data$Participant_id))
if (!"snp_strain_context" %in% names(casebook)) {
  casebook <- casebook %>%
    mutate(
      snp_strain_context = case_when(
        is.na(.data$SNPs) ~ "Missing SNP evidence",
        .data$SNPs <= strain_snp_threshold() ~ "Strong same strain",
        TRUE ~ "Above same-strain SNP threshold"
      )
    )
}
if (!"st_lineage_context" %in% names(casebook)) {
  casebook <- casebook %>%
    mutate(
      st_lineage_context = case_when(
        .data$same_ST %in% TRUE ~ "Same ST",
        .data$same_ST %in% FALSE ~ "Different ST",
        TRUE ~ "Missing ST evidence"
      )
    )
}
if (!"pair_interpretation" %in% names(casebook)) {
  casebook <- casebook %>% mutate(pair_interpretation = .data$same_strain_evidence)
}
mechanism_summary <- read_current(required_paths[["mechanism_summary"]])
denominator <- read_current(required_paths[["denominator"]])
score_contrasts <- read_current(required_paths[["score_contrasts"]])
model_stability <- read_current(required_paths[["model_stability"]])
leave_one_summary <- read_current(required_paths[["leave_one_summary"]])
bootstrap <- read_current(required_paths[["bootstrap"]])
near_miss <- read_current(required_paths[["near_miss"]]) %>%
  mutate(Participant_id = as.character(.data$Participant_id))
power_precision <- read_current(required_paths[["power_precision"]])
leave_one_detail <- read_current(required_paths[["leave_one_detail"]])
variants <- read_current(required_paths[["variants"]]) %>%
  mutate(Participant_id = as.character(.data$Participant_id))
stat_score_values <- read_current(required_paths[["stat_score_values"]]) %>%
  mutate(
    Participant_id = as.character(.data$Participant_id),
    UTI_Status = factor(.data$UTI_Status, levels = c("Not_UTI", "UTI"))
  )
stat_module_matrix <- read_current(required_paths[["stat_module_matrix"]]) %>%
  mutate(Participant_id = as.character(.data$Participant_id))
stat_transition_module <- read_current(required_paths[["stat_transition_module"]])
stat_pcoa <- read_current(required_paths[["stat_pcoa"]]) %>%
  prefer_primary_uti_status(allow_legacy_fallback = FALSE) %>%
  mutate(
    Participant_id = as.character(.data$Participant_id),
    UTI_Status = factor(.data$UTI_Status, levels = c("Not_UTI", "UTI"))
  )

metric_num <- function(df, metric_name) {
  out <- df %>% filter(.data$metric == metric_name) %>% pull(.data$value)
  if (length(out) != 1) stop("Expected one denominator metric: ", metric_name)
  suppressWarnings(as.numeric(out))
}

clinical_total <- metric_num(denominator, "primary_clinical_total")
clinical_uti <- metric_num(denominator, "primary_clinical_uti")
clinical_not <- metric_num(denominator, "primary_clinical_not_uti")
vf_total <- metric_num(denominator, "primary_vf_total")
vf_uti <- metric_num(denominator, "primary_vf_uti")
vf_not <- metric_num(denominator, "primary_vf_not_uti")
near_miss_clinical <- metric_num(denominator, "near_miss_clinical_rows")
near_miss_vf <- metric_num(denominator, "near_miss_vf_ready_rows")
expanded_uti <- metric_num(denominator, "expanded_vf_uti_if_near_miss_included")
expanded_not <- metric_num(denominator, "expanded_vf_not_uti_if_near_miss_excluded")

vf_keys <- vf_ready %>% distinct(.data$Participant_id, .data$tp_lab)
primary_retention <- status %>%
  left_join(vf_keys %>% mutate(retained_for_vf = TRUE),
            by = c("Participant_id", "tp_lab")) %>%
  mutate(
    retained_for_vf = coalesce(.data$retained_for_vf, FALSE),
    retention_state = if_else(.data$retained_for_vf, "VF/model retained", "No VF/model endpoint")
  ) %>%
  count(.data$UTI_Status, .data$retention_state, name = "n") %>%
  complete(
    UTI_Status = c("Not_UTI", "UTI"),
    retention_state = c("VF/model retained", "No VF/model endpoint"),
    fill = list(n = 0L)
  )

check_row <- function(check, ok, detail) {
  tibble(check = check, status = if_else(ok, "PASS", "FAIL"), detail = as.character(detail))
}

validation <- bind_rows(
  check_row("all upstream mechanism checks passed", all(mechanism_validation$status == "PASS"),
            sprintf("%d mechanism checks imported", nrow(mechanism_validation))),
  check_row("all upstream robustness checks passed", all(robustness_validation$status == "PASS"),
            sprintf("%d robustness checks imported", nrow(robustness_validation))),
  check_row("primary clinical denominator is 583 with 18 UTI",
            clinical_total == 583 && clinical_uti == 18 && clinical_not == 565,
            sprintf("n=%d; UTI=%d; Not_UTI=%d", clinical_total, clinical_uti, clinical_not)),
  check_row("primary VF/model denominator is 556 with 17 UTI",
            vf_total == 556 && vf_uti == 17 && vf_not == 539,
            sprintf("n=%d; UTI=%d; Not_UTI=%d", vf_total, vf_uti, vf_not)),
  check_row("near-miss denominator matches current sensitivity outputs",
            near_miss_clinical == 19 && near_miss_vf == 18 &&
              expanded_uti == 35 && expanded_not == 521,
            sprintf("clinical=%d; VF-ready=%d; expanded=%d/%d",
                    near_miss_clinical, near_miss_vf, expanded_uti, expanded_not)),
  check_row("casebook includes 11 clinical Not_UTI-to-UTI transitions",
            nrow(casebook) == 11,
            sprintf("n=%d", nrow(casebook))),
  check_row("casebook includes 10 WGS/VF-linked and 1 missing endpoint",
            sum(casebook$has_vf_pair %in% TRUE) == 10 &&
              sum(!(casebook$has_vf_pair %in% TRUE)) == 1,
            sprintf("linked=%d; missing=%d",
                    sum(casebook$has_vf_pair %in% TRUE),
                    sum(!(casebook$has_vf_pair %in% TRUE)))),
  check_row("casebook includes no Uricult-linked transition",
            sum(casebook$is_uricult_transition %in% TRUE, na.rm = TRUE) == 0,
            sprintf("uricult=%d", sum(casebook$is_uricult_transition %in% TRUE, na.rm = TRUE))),
  check_row("primary-key retention uses primary denominator directly",
            sum(primary_retention$n[primary_retention$UTI_Status == "Not_UTI" &
                                      primary_retention$retention_state == "VF/model retained"]) == 539 &&
              sum(primary_retention$n[primary_retention$UTI_Status == "Not_UTI" &
                                      primary_retention$retention_state == "No VF/model endpoint"]) == 26 &&
              sum(primary_retention$n[primary_retention$UTI_Status == "UTI" &
                                      primary_retention$retention_state == "VF/model retained"]) == 17 &&
              sum(primary_retention$n[primary_retention$UTI_Status == "UTI" &
                                      primary_retention$retention_state == "No VF/model endpoint"]) == 1,
            "Direct primary-key join: Not_UTI 539 retained/26 missing; UTI 17 retained/1 missing."),
  check_row("no legacy or OLD source path used",
            !any(str_detect(required_paths, regex("_OLD|old_asb|legacy", ignore_case = TRUE))),
            "Inputs are restricted to current primary-status outputs.")
)

write_csv(validation, file.path(DIR_FINAL_RESULTS, "final_figure_validation_checks.csv"))
if (any(validation$status != "PASS")) {
  stop("Final figure validation failed; see results/final_figures/final_figure_validation_checks.csv.")
}

amr_report_path <- file.path(DIR_MECHANISM, "amr_screen_report.txt")
amr_complete <- file.exists(amr_report_path) &&
  any(str_detect(readLines(amr_report_path, warn = FALSE), fixed("Status: COMPLETE")))

mechanism_levels <- c(
  "same_strain_stable_profile",
  "same_strain_genomic_change",
  "strain_replacement",
  "uncertain",
  "missing_wgs_endpoint"
)
mechanism_labels <- c(
  same_strain_stable_profile = "Same strain, stable profile",
  same_strain_genomic_change = "Same strain, profile change",
  strain_replacement = "Strain replacement",
  uncertain = "Uncertain",
  missing_wgs_endpoint = "Missing genomic endpoint"
)
mechanism_cols <- c(
  same_strain_stable_profile = "#2E8B57",
  same_strain_genomic_change = "#E69F00",
  strain_replacement = "#C44E19",
  uncertain = "#718096",
  missing_wgs_endpoint = "#BDBDBD"
)

final_theme <- function(base_size = 10) {
  theme_bw(base_size = base_size) +
    theme(
      plot.title = element_text(face = "bold", size = rel(1.12)),
      plot.subtitle = element_text(colour = "grey25", size = rel(0.92)),
      plot.caption = element_text(colour = "grey35", hjust = 0, size = rel(0.78)),
      legend.position = "bottom",
      legend.title = element_text(face = "bold"),
      strip.background = element_rect(fill = "#F3F4F6", colour = "grey75"),
      strip.text = element_text(face = "bold"),
      panel.grid.minor = element_blank(),
      axis.title = element_text(face = "bold")
    )
}

human_case <- function(participant, from_tp, to_tp) {
  paste0(participant, ": ", from_tp, " to ", to_tp)
}

normalise_st_label <- function(x) {
  x <- str_trim(as.character(x))
  unknown <- c("", "-", "ST-", "NA", "N/A", "UNKNOWN", "UNK", "NT",
               "NON-TYPABLE", "NONTYPABLE", "NOT TYPED")
  x[str_to_upper(x) %in% unknown] <- NA_character_
  x
}

st_group4 <- function(x) {
  st <- normalise_st_label(x)
  case_when(
    st == "131" ~ "ST131",
    st == "141" ~ "ST141",
    is.na(st) ~ "Missing ST",
    TRUE ~ "Other typed ST"
  )
}

casebook <- casebook %>%
  mutate(
    mechanism_bucket = factor(.data$mechanism_bucket, levels = mechanism_levels),
    mechanism_label = factor(mechanism_labels[as.character(.data$mechanism_bucket)],
                             levels = mechanism_labels[mechanism_levels]),
    case_label = human_case(.data$Participant_id, .data$from_tp, .data$to_tp)
  ) %>%
  arrange(.data$mechanism_bucket, is.na(.data$SNPs), .data$SNPs, .data$Participant_id)

case_order <- rev(casebook$case_label)
casebook <- casebook %>% mutate(case_label = factor(.data$case_label, levels = case_order))

manifest <- tibble()
register_figure <- function(figure_id, figure_class, plot, width, height, caption,
                            source_inputs, evidence_type, limitations, include_amr = FALSE) {
  png_path <- file.path(DIR_FINAL_PLOTS, paste0(figure_id, ".png"))
  pdf_path <- file.path(DIR_FINAL_PLOTS, paste0(figure_id, ".pdf"))
  ggsave(png_path, plot, width = width, height = height, dpi = 300, bg = "white")
  ggsave(pdf_path, plot, width = width, height = height, device = "pdf", bg = "white")
  manifest <<- bind_rows(
    manifest,
    tibble(
      figure_id = figure_id,
      figure_class = figure_class,
      png_path = png_path,
      pdf_path = pdf_path,
      caption = caption,
      source_inputs = paste(source_inputs, collapse = "; "),
      evidence_type = evidence_type,
      interpretation_limitations = limitations,
      amr_screen_included = include_amr && amr_complete
    )
  )
}

# ------------------------------------------------------------------------------
# Main Figure 1: Denominators and sensitivity framing
# ------------------------------------------------------------------------------
clinical_plot_data <- tibble(
  Status = factor(c("Not_UTI", "UTI"), levels = c("Not_UTI", "UTI")),
  n = c(clinical_not, clinical_uti),
  cohort = "Clinical episodes"
)
p1a <- ggplot(clinical_plot_data, aes(.data$cohort, .data$n, fill = .data$Status)) +
  geom_col(width = 0.62, colour = "white") +
  geom_text(aes(label = .data$n), position = position_stack(vjust = 0.5),
            colour = "white", fontface = "bold", size = 3.6) +
  scale_fill_uti_status(drop = FALSE) +
  labs(
    title = "Primary clinical status",
    subtitle = "583 episodes: 18 UTI, 565 Not_UTI",
    x = NULL, y = "Episodes", fill = "Primary status"
  ) +
  final_theme(10) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

retention_plot_data <- primary_retention %>%
  mutate(
    UTI_Status = factor(.data$UTI_Status, levels = c("Not_UTI", "UTI")),
    retention_state = factor(.data$retention_state,
                             levels = c("VF/model retained", "No VF/model endpoint"))
  )
p1b <- ggplot(retention_plot_data, aes(.data$n, .data$retention_state, fill = .data$retention_state)) +
  geom_col(width = 0.6) +
  geom_text(aes(label = if_else(.data$n > 0, as.character(.data$n), "")),
            hjust = -0.2, size = 3.3, fontface = "bold") +
  scale_fill_manual(values = c("VF/model retained" = "#2E8B57",
                               "No VF/model endpoint" = "#D9A441")) +
  facet_grid(rows = vars(UTI_Status), scales = "free_y") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(
    title = "Clinical-to-VF retention",
    subtitle = "556 retained: 17 UTI, 539 Not_UTI",
    x = "Episodes", y = NULL, fill = NULL
  ) +
  final_theme(10) +
  theme(axis.text.y = element_text(size = 8))

sensitivity_data <- tibble(
  analysis = factor(rep(c("Primary rule", "Possible-UTI\nsensitivity"), each = 2),
                    levels = c("Primary rule", "Possible-UTI\nsensitivity")),
  segment = factor(c("Reference Not_UTI", "UTI", "Reference Not_UTI", "Possible UTI"),
                   levels = c("Reference Not_UTI", "UTI", "Possible UTI")),
  n = c(vf_not, vf_uti, expanded_not, expanded_uti)
)
p1c <- ggplot(sensitivity_data, aes(.data$analysis, .data$n, fill = .data$segment)) +
  geom_col(width = 0.62, colour = "white") +
  geom_text(aes(label = .data$n), position = position_stack(vjust = 0.5),
            size = 3.5, fontface = "bold") +
  scale_fill_manual(values = c("Reference Not_UTI" = uti_status_cols[["Not_UTI"]],
                               "UTI" = uti_status_cols[["UTI"]],
                               "Possible UTI" = "#E69F00")) +
  labs(
    title = "Near-miss sensitivity",
    subtitle = "18 VF-ready near-miss rows are sensitivity only",
    x = NULL, y = "VF/model episodes", fill = NULL
  ) +
  final_theme(10)

main_fig_1 <- (p1a | p1b | p1c) +
  plot_annotation(
    title = "Sparse UTI denominators define the interpretation boundary",
    subtitle = "Primary analyses use UTI_Status; near-miss rows are never relabelled in the main analysis.",
    caption = "Descriptive denominator audit. Retention is derived by joining primary clinical episode keys directly to VF-ready keys.",
    tag_levels = "A",
    theme = theme(plot.title = element_text(face = "bold", size = 16),
                  plot.subtitle = element_text(size = 11),
                  plot.caption = element_text(size = 9, colour = "grey35"))
  )
register_figure(
  "primary_denominator_and_uncertainty", "Main", main_fig_1, 15, 5.8,
  "Primary status and VF-ready denominators, with near-miss sensitivity framing. Primary clinical analyses include 18 UTI and 565 Not_UTI episodes; genomic analyses include 17 UTI and 539 Not_UTI episodes.",
  required_paths[c("status_map", "vf_ready", "denominator")],
  "Descriptive denominator audit", "Near-miss rows are sensitivity-only and are not primary UTI events."
)

# ------------------------------------------------------------------------------
# Main Figure 2: Mechanism casebook
# ------------------------------------------------------------------------------
case_counts <- casebook %>%
  count(.data$mechanism_bucket, name = "n") %>%
  complete(mechanism_bucket = factor(mechanism_levels, levels = mechanism_levels),
           fill = list(n = 0L)) %>%
  mutate(mechanism_label = factor(mechanism_labels[as.character(.data$mechanism_bucket)],
                                  levels = rev(mechanism_labels[mechanism_levels])))

p2a <- ggplot(case_counts, aes(.data$n, .data$mechanism_label, fill = .data$mechanism_bucket)) +
  geom_col(width = 0.65) +
  geom_text(aes(label = .data$n), hjust = -0.35, fontface = "bold", size = 3.6) +
  scale_fill_manual(values = mechanism_cols, guide = "none") +
  scale_x_continuous(limits = c(0, max(case_counts$n) + 1), breaks = pretty_breaks()) +
  labs(title = "Mechanism categories", subtitle = "11 clinical Not_UTI to UTI transitions",
       x = "Transitions", y = NULL) +
  final_theme(10)

case_tiles <- bind_rows(
  casebook %>% transmute(case_label, feature = "WGS/VF endpoint",
                         value = if_else(.data$has_vf_pair %in% TRUE, "Available", "Missing")),
  casebook %>% transmute(case_label, feature = "SNP context",
                         value = case_when(
                           .data$snp_strain_context == "Strong same strain" ~ "Strong same strain",
                           .data$snp_strain_context == "Above same-strain SNP threshold" ~ "Above SNP threshold",
                           .data$snp_strain_context == "Missing SNP evidence" ~ "Missing SNP",
                           TRUE ~ "Uncertain")),
  casebook %>% transmute(case_label, feature = "ST context",
                         value = case_when(
                           .data$st_lineage_context == "Same ST" ~ "Same ST",
                           .data$st_lineage_context == "Different ST" ~ "Different ST",
                           .data$st_lineage_context == "Missing ST evidence" ~ "Missing ST",
                           TRUE ~ "Uncertain")),
  casebook %>% transmute(case_label, feature = "VF/module profile",
                         value = case_when(
                           !(.data$has_vf_pair %in% TRUE) ~ "Missing",
                           coalesce(.data$n_vf_genes_gained, 0) + coalesce(.data$n_vf_genes_lost, 0) +
                             coalesce(.data$n_modules_gained, 0) + coalesce(.data$n_modules_lost, 0) > 0 ~ "Changed",
                           TRUE ~ "Stable")),
  casebook %>% transmute(case_label, feature = "Plasmid/AMR profile",
                         value = case_when(
                           !(.data$has_vf_pair %in% TRUE) ~ "Missing",
                           (coalesce(.data$n_plasmid_gained, 0) + coalesce(.data$n_plasmid_lost, 0) +
                              if (amr_complete) {
                                coalesce(.data$n_amr_gained, 0) + coalesce(.data$n_amr_lost, 0)
                              } else {
                                0
                              }) > 0 ~ "Changed",
                           TRUE ~ "Stable")),
  casebook %>% transmute(case_label, feature = "Symptom state",
                         value = if_else(.data$symptom_state_change == "symptoms_emerged",
                                         "Symptoms emerged", "No recorded emergence")),
  casebook %>% transmute(case_label, feature = "Collection/catheter",
                         value = if_else(coalesce(.data$collection_method_changed, FALSE) |
                                           coalesce(.data$catheter_rule_changed, FALSE),
                                         "Changed context", "No context shift"))
) %>%
  mutate(feature = factor(.data$feature,
                          levels = c("WGS/VF endpoint", "SNP context", "ST context", "VF/module profile",
                                     "Plasmid/AMR profile", "Symptom state", "Collection/catheter")))

tile_cols <- c(
  "Available" = "#4C9F70", "Missing" = "#D0D5DD",
  "Strong same strain" = mechanism_cols[["same_strain_stable_profile"]],
  "Above SNP threshold" = "#D69E2E", "Missing SNP" = "#9CA3AF",
  "Same ST" = "#009E73", "Different ST" = "#D55E00", "Missing ST" = "#9CA3AF",
  "Uncertain" = mechanism_cols[["uncertain"]],
  "Stable" = "#A7D7C5", "Changed" = "#E69F00",
  "Symptoms emerged" = "#6B5B95", "No recorded emergence" = "#D0D5DD",
  "Changed context" = "#4C78A8", "No context shift" = "#E5E7EB"
)
p2b <- ggplot(case_tiles, aes(.data$feature, .data$case_label, fill = .data$value)) +
  geom_tile(colour = "white", linewidth = 0.7) +
  scale_fill_manual(values = tile_cols, name = "Observed state") +
  labs(
    title = "Casewise evidence matrix",
    subtitle = "Ordered by mechanism category and SNP evidence",
    x = NULL, y = NULL
  ) +
  final_theme(9) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1),
        legend.position = "bottom")

main_fig_2 <- p2a + p2b +
  plot_layout(widths = c(0.34, 0.66), guides = "collect") +
  plot_annotation(
    title = "Not_UTI-to-UTI transitions show heterogeneous bacterial context",
    subtitle = "Four events show same-strain stable profiles; three are consistent with strain replacement.",
    caption = "Descriptive mechanism classification; all 11 clinical transitions are retained, including one with a missing genomic endpoint.",
    tag_levels = "A",
    theme = theme(plot.title = element_text(face = "bold", size = 16),
                  plot.subtitle = element_text(size = 11),
                  plot.caption = element_text(size = 9, colour = "grey35"))
  ) &
  theme(legend.position = "bottom")
register_figure(
  "not_uti_to_uti_mechanism_casebook", "Main", main_fig_2, 16, 8.2,
  "Mechanism classification and evidence matrix for the 11 clinical Not_UTI-to-UTI transitions. The casebook includes 10 WGS/VF-linked transitions and one missing genomic endpoint.",
  c(required_paths[c("casebook", "mechanism_validation")], amr_report_path),
  "Descriptive mechanism casebook", "Buckets organise evidence; they do not prove a biological mechanism."
  , include_amr = TRUE
)

# ------------------------------------------------------------------------------
# Main Figure 3: Strain/VF evidence and host context
# ------------------------------------------------------------------------------
wgs_cases <- casebook %>% filter(.data$has_vf_pair %in% TRUE, !is.na(.data$SNPs), !is.na(.data$vf_jaccard))
p3a <- ggplot(wgs_cases, aes(.data$SNPs, .data$vf_jaccard, colour = .data$mechanism_bucket,
                             label = .data$Participant_id)) +
  geom_point(size = 3.2) +
  ggrepel::geom_text_repel(size = 3, show.legend = FALSE, max.overlaps = Inf,
                           box.padding = 0.35, point.padding = 0.25) +
  scale_x_log10(labels = label_number(big.mark = ",")) +
  scale_colour_manual(values = mechanism_cols, labels = mechanism_labels, name = "Mechanism category") +
  coord_cartesian(ylim = c(0.45, 1.04)) +
  labs(
    title = "Strain distance and VF similarity",
    subtitle = "10 WGS/VF-linked Not_UTI to UTI cases",
    x = "SNP distance (log scale)", y = "VF Jaccard similarity"
  ) +
  final_theme(10)

host_tiles <- bind_rows(
  casebook %>% transmute(case_label, feature = "Symptoms emerged",
                         value = .data$symptom_state_change == "symptoms_emerged"),
  casebook %>% transmute(case_label, feature = "Culture support persisted",
                         value = .data$culture_state_change == "culture_support_persisted"),
  casebook %>% transmute(case_label, feature = "Collection method changed",
                         value = coalesce(.data$collection_method_changed, FALSE)),
  casebook %>% transmute(case_label, feature = "Catheter rule changed",
                         value = coalesce(.data$catheter_rule_changed, FALSE))
) %>%
  mutate(feature = factor(.data$feature,
                          levels = c("Symptoms emerged", "Culture support persisted",
                                     "Collection method changed", "Catheter rule changed",
                                     "Days between samples")),
         display = if_else(.data$value, "Yes", "No"))
interval_tiles <- casebook %>%
  transmute(case_label, feature = factor("Days between samples",
                                         levels = levels(host_tiles$feature)),
            display = paste0(.data$days_between, " d"))

p3b <- ggplot(host_tiles, aes(.data$feature, .data$case_label, fill = .data$display)) +
  geom_tile(colour = "white", linewidth = 0.7) +
  geom_text(aes(label = .data$display), size = 2.8) +
  geom_tile(data = interval_tiles, aes(.data$feature, .data$case_label),
            inherit.aes = FALSE, fill = "#F1F3F5", colour = "white", linewidth = 0.7) +
  geom_text(data = interval_tiles, aes(.data$feature, .data$case_label, label = .data$display),
            inherit.aes = FALSE, size = 2.8) +
  scale_fill_manual(values = c("Yes" = "#6B5B95", "No" = "#E5E7EB"), name = NULL) +
  labs(
    title = "Clinical context around transition",
    subtitle = "UTI emerges on pre-existing culture support in all 11 cases",
    x = NULL, y = NULL
  ) +
  final_theme(9) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1))

main_fig_3 <- p3a + p3b +
  plot_layout(widths = c(0.42, 0.58), guides = "collect") +
  plot_annotation(
    title = "Stable bacterial profiles coexist with changing clinical state",
    subtitle = "Same-strain, VF-stable cases motivate host-state and regulatory hypotheses rather than simple VF acquisition.",
    caption = "SNP distance, VF similarity, and host-context observations are descriptive; UTI case counts are sparse.",
    tag_levels = "A",
    theme = theme(plot.title = element_text(face = "bold", size = 16),
                  plot.subtitle = element_text(size = 11),
                  plot.caption = element_text(size = 9, colour = "grey35"))
  ) &
  theme(legend.position = "bottom")
register_figure(
  "strain_stability_and_host_context", "Main", main_fig_3, 16, 8,
  "Strain distance and VF similarity alongside clinical-context changes among Not_UTI-to-UTI transitions.",
  required_paths[c("casebook")],
  "Descriptive transition context", "Low SNP distance and stable VF profiles do not distinguish host state from unmeasured bacterial regulation."
)

# ------------------------------------------------------------------------------
# Main Figure 4: Global VF signal and robustness
# ------------------------------------------------------------------------------
forest_labels <- c(
  expec_marker_count = "ExPEC-like markers",
  upec_system_count = "UPEC systems",
  total_vf_count_all = "All VF genes",
  total_vf_count_curated = "Curated VF genes"
)
forest <- bootstrap %>%
  filter(.data$score %in% names(forest_labels), .data$statistic == "mean difference") %>%
  mutate(score_label = factor(forest_labels[.data$score],
                              levels = rev(unname(forest_labels))))
p4a <- ggplot(forest, aes(.data$observed_effect, .data$score_label)) +
  geom_vline(xintercept = 0, colour = "grey35", linewidth = 0.5) +
  geom_errorbarh(aes(xmin = .data$ci_lower, xmax = .data$ci_upper),
                 height = 0.14, colour = "#4C78A8", linewidth = 0.7) +
  geom_point(size = 2.8, colour = "#0F766E") +
  labs(
    title = "Participant-bootstrap effects",
    subtitle = "Mean difference: UTI minus Not_UTI",
    x = "Difference in endpoint value", y = NULL
  ) +
  final_theme(9)

get_stability <- function(metric) {
  suppressWarnings(as.numeric(model_stability$value[model_stability$metric == metric][1]))
}
flag_data <- tibble(
  flag = factor(c("GLMM sparse/separation risk", "Leave-one direction flip",
                  "GLMM singular fit", "GLMM FDR < 0.05", "Univariable FDR < 0.05"),
                levels = rev(c("GLMM sparse/separation risk", "Leave-one direction flip",
                               "GLMM singular fit", "GLMM FDR < 0.05", "Univariable FDR < 0.05"))),
  n = c(get_stability("glmm_sparse_or_separation_risk"),
        sum(leave_one_summary$n_direction_flip),
        get_stability("glmm_singular_models"),
        get_stability("glmm_fdr_lt_0_05"),
        get_stability("univariable_fdr_lt_0_05")),
  denom = c(get_stability("glmm_models_total"),
            sum(leave_one_summary$n_features),
            get_stability("glmm_models_total"),
            get_stability("glmm_models_total"),
            get_stability("univariable_tests_total"))
) %>%
  mutate(label = paste0(.data$n, "/", .data$denom))
p4b <- ggplot(flag_data, aes(.data$n, .data$flag)) +
  geom_col(fill = "#C45A24", width = 0.63) +
  geom_text(aes(label = .data$label), hjust = -0.15, size = 3.1) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.18))) +
  labs(title = "Stability flags", subtitle = "Sparse data dominate model interpretation",
       x = "Flagged results", y = NULL) +
  final_theme(9)

contrast_labels <- c(
  expec_marker_count = "ExPEC-like markers",
  upec_system_count = "UPEC systems",
  total_vf_count_curated = "Curated VF genes",
  total_vf_count_all = "All VF genes"
)
contrast_plot <- score_contrasts %>%
  filter(.data$score %in% names(contrast_labels),
         .data$contrast %in% c("primary_UTI_vs_all_Not_UTI",
                               "expanded_UTI_vs_Not_UTI_excluding_near_miss")) %>%
  mutate(
    score_label = factor(contrast_labels[.data$score],
                         levels = rev(unname(contrast_labels))),
    sensitivity = recode(.data$contrast,
                         primary_UTI_vs_all_Not_UTI = "Primary rule (17 UTI)",
                         expanded_UTI_vs_Not_UTI_excluding_near_miss = "Possible-UTI sensitivity (35)")
  )
p4c <- ggplot(contrast_plot, aes(.data$median_difference_a_minus_b, .data$score_label,
                                 colour = .data$sensitivity)) +
  geom_vline(xintercept = 0, colour = "grey35", linewidth = 0.5) +
  geom_point(size = 3, position = position_dodge(width = 0.42)) +
  scale_colour_manual(values = c("Primary rule (17 UTI)" = uti_status_cols[["UTI"]],
                                 "Possible-UTI sensitivity (35)" = "#E69F00")) +
  labs(
    title = "Near-miss sensitivity",
    subtitle = "Median difference versus reference Not_UTI",
    x = "Difference in endpoint value", y = NULL, colour = NULL
  ) +
  final_theme(9)

main_fig_4 <- (p4a | p4b) / p4c +
  plot_layout(heights = c(0.55, 0.45), guides = "collect") +
  plot_annotation(
    title = "Global VF associations remain exploratory and unstable",
    subtitle = "No GLMM or univariable feature association remains significant after FDR correction.",
    caption = "Participant bootstrap, model diagnostics, and near-miss sensitivity are interpretive robustness checks, not confirmatory inference.",
    tag_levels = "A",
    theme = theme(plot.title = element_text(face = "bold", size = 16),
                  plot.subtitle = element_text(size = 11),
                  plot.caption = element_text(size = 9, colour = "grey35"))
  ) &
  theme(legend.position = "bottom")
register_figure(
  "global_vf_signal_and_robustness", "Main", main_fig_4, 15, 9,
  "Participant-bootstrap VF supplementary endpoint effects, model stability flags, and near-miss sensitivity contrasts in the VF/model cohort.",
  required_paths[c("bootstrap", "model_stability", "leave_one_summary", "score_contrasts")],
  "Exploratory robustness diagnostics", "No adjusted association is confirmatory; only 17 VF-ready primary UTI episodes are available."
)

# ------------------------------------------------------------------------------
# Supplementary Figure S1: Mechanisms by transition type
# ------------------------------------------------------------------------------
transition_types <- c("Not_UTI->Not_UTI", "Not_UTI->UTI", "UTI->Not_UTI")
s1_data <- mechanism_summary %>%
  filter(.data$transition_type %in% transition_types) %>%
  complete(transition_type = transition_types,
           mechanism_bucket = mechanism_levels,
           fill = list(n_transitions = 0)) %>%
  mutate(
    transition_label = factor(recode(.data$transition_type,
                                     "Not_UTI->Not_UTI" = "Not_UTI to Not_UTI",
                                     "Not_UTI->UTI" = "Not_UTI to UTI",
                                     "UTI->Not_UTI" = "UTI to Not_UTI"),
                              levels = c("Not_UTI to Not_UTI", "Not_UTI to UTI", "UTI to Not_UTI")),
    mechanism_bucket = factor(.data$mechanism_bucket, levels = mechanism_levels)
  ) %>%
  group_by(.data$transition_label) %>%
  mutate(total = sum(.data$n_transitions),
         fraction = .data$n_transitions / .data$total,
         n_label = if_else(.data$fraction >= 0.07 | .data$transition_label == "Not_UTI to UTI",
                           as.character(.data$n_transitions), ""))
s1 <- ggplot(s1_data, aes(.data$transition_label, .data$fraction, fill = .data$mechanism_bucket)) +
  geom_col(width = 0.65, colour = "white") +
  geom_text(aes(label = .data$n_label), position = position_stack(vjust = 0.5),
            size = 3.4, colour = "white", fontface = "bold") +
  scale_fill_manual(values = mechanism_cols, labels = mechanism_labels, name = "Mechanism category",
                    guide = guide_legend(nrow = 2, byrow = TRUE)) +
  scale_y_continuous(labels = percent_format()) +
  labs(
    title = "Transition mechanism composition",
    subtitle = "Mechanism buckets are descriptive and transition denominators differ markedly",
    x = NULL, y = "Proportion of transition type",
    caption = "Not_UTI-to-UTI n=11; UTI-to-Not_UTI n=14; Not_UTI-to-Not_UTI n=391."
  ) +
  final_theme(11) +
  theme(legend.text = element_text(size = 9))
register_figure(
  "transition_mechanisms_by_transition_type", "Supplementary", s1, 9, 6.4,
  "Descriptive composition of mechanism categories across within-resident transition types.",
  required_paths[c("mechanism_summary")],
  "Descriptive comparison", "Unequal transition denominators prevent inferential comparison of proportions."
)

# ------------------------------------------------------------------------------
# Supplementary Figure S2: Accessory, plasmid, and AMR changes
# ------------------------------------------------------------------------------
s2_accessory <- casebook %>%
  transmute(
    case_label,
    gained = coalesce(.data$n_accessory_genes_gained, 0),
    lost = coalesce(.data$n_accessory_genes_lost, 0)
  ) %>%
  pivot_longer(c("gained", "lost"), names_to = "direction", values_to = "n") %>%
  mutate(
    signed_log = if_else(.data$direction == "gained", log10(.data$n + 1), -log10(.data$n + 1)),
    label = if_else(.data$n == 0, "", as.character(.data$n))
  )
p_s2a <- ggplot(s2_accessory, aes(.data$signed_log, .data$case_label, fill = .data$direction)) +
  geom_vline(xintercept = 0, colour = "grey45") +
  geom_col(width = 0.67) +
  geom_text(aes(label = .data$label),
            hjust = ifelse(s2_accessory$direction == "gained", -0.1, 1.1),
            size = 2.8) +
  scale_fill_manual(values = c(gained = "#2E8B57", lost = "#C44E19"),
                    labels = c(gained = "Gained", lost = "Lost"), name = "Accessory genes") +
  scale_x_continuous(
    breaks = c(-log10(1001), -log10(101), -log10(11), 0,
               log10(11), log10(101), log10(1001)),
    labels = c("-1000", "-100", "-10", "0", "+10", "+100", "+1000")
  ) +
  labs(title = "Accessory-gene changes", subtitle = "Signed log scale; count labels show raw Panaroo changes",
       x = "Genes lost / gained", y = NULL) +
  final_theme(9)

s2_mobile <- bind_rows(
  casebook %>% transmute(case_label, feature = "Plasmid gained", n = coalesce(.data$n_plasmid_gained, 0)),
  casebook %>% transmute(case_label, feature = "Plasmid lost", n = coalesce(.data$n_plasmid_lost, 0)),
  casebook %>% transmute(case_label, feature = "AMR gained", n = if (amr_complete) coalesce(.data$n_amr_gained, 0) else NA_real_),
  casebook %>% transmute(case_label, feature = "AMR lost", n = if (amr_complete) coalesce(.data$n_amr_lost, 0) else NA_real_)
) %>%
  mutate(
    feature = factor(.data$feature, levels = c("Plasmid gained", "Plasmid lost", "AMR gained", "AMR lost")),
    category = case_when(is.na(.data$n) ~ "Unavailable", .data$n > 0 ~ "Change detected", TRUE ~ "No change"),
    label = if_else(is.na(.data$n), "NA", as.character(.data$n))
  )
p_s2b <- ggplot(s2_mobile, aes(.data$feature, .data$case_label, fill = .data$category)) +
  geom_tile(colour = "white", linewidth = 0.7) +
  geom_text(aes(label = .data$label), size = 2.8) +
  scale_fill_manual(values = c("Change detected" = "#E69F00", "No change" = "#E5E7EB",
                               "Unavailable" = "#BDBDBD"), name = NULL) +
  labs(title = "Plasmid and AMR changes",
       subtitle = if_else(amr_complete, "AMR screened by cached ABRicate ResFinder output",
                          "AMR output unavailable; AMR tiles intentionally blank"),
       x = NULL, y = NULL) +
  final_theme(9) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1))
s2 <- p_s2a + p_s2b +
  plot_layout(widths = c(0.55, 0.45), guides = "collect") +
  plot_annotation(
    title = "Exploratory accessory and mobile-resistance context",
    caption = "Accessory-gene and AMR changes are exploratory leads and do not establish a UTI mechanism.",
    tag_levels = "A",
    theme = theme(plot.title = element_text(face = "bold", size = 15),
                  plot.caption = element_text(size = 9, colour = "grey35"))
  ) &
  theme(legend.position = "bottom")
register_figure(
  "accessory_plasmid_amr_changes", "Supplementary", s2, 15, 8.4,
  "Exploratory Panaroo accessory-gene, plasmid-replicon, and completed ResFinder AMR changes for Not_UTI-to-UTI cases.",
  c(required_paths[["casebook"]], amr_report_path),
  "Exploratory genome screen context", "Accessory calls may include annotation or presence-absence noise; AMR changes are not VF evidence.",
  include_amr = TRUE
)

# ------------------------------------------------------------------------------
# Supplementary Figure S3: Near-miss rule and sparse precision
# ------------------------------------------------------------------------------
near_features <- bind_rows(
  near_miss %>% transmute(case_label = human_case(.data$Participant_id, .data$tp_lab, ""),
                          feature = "Culture supports UTI", value = .data$culture_supports_uti),
  near_miss %>% transmute(case_label = human_case(.data$Participant_id, .data$tp_lab, ""),
                          feature = "Local urinary symptom", value = .data$local_urinary_symptom_any),
  near_miss %>% transmute(case_label = human_case(.data$Participant_id, .data$tp_lab, ""),
                          feature = "Systemic symptom", value = .data$systemic_symptom_any),
  near_miss %>% transmute(case_label = human_case(.data$Participant_id, .data$tp_lab, ""),
                          feature = "Symptom rule met", value = .data$symptom_compatible_uti),
  near_miss %>% transmute(case_label = human_case(.data$Participant_id, .data$tp_lab, ""),
                          feature = "Indwelling catheter rule", value = .data$catheter_rule == "B_indwelling")
) %>%
  mutate(
    feature = factor(.data$feature, levels = c("Culture supports UTI", "Local urinary symptom",
                                               "Systemic symptom", "Symptom rule met",
                                               "Indwelling catheter rule")),
    case_label = factor(.data$case_label, levels = rev(unique(.data$case_label))),
    display = if_else(coalesce(.data$value, FALSE), "Yes", "No")
  )
p_s3a <- ggplot(near_features, aes(.data$feature, .data$case_label, fill = .data$display)) +
  geom_tile(colour = "white", linewidth = 0.55) +
  scale_fill_manual(values = c(Yes = "#D55E00", No = "#E5E7EB"), name = NULL) +
  labs(title = "Near-miss evidence audit", subtitle = "19 clinical rows remain Not_UTI under the primary rule",
       x = NULL, y = NULL) +
  final_theme(8.5) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1))

p_s3b <- power_precision %>%
  mutate(
    odds_ratio = factor(.data$odds_ratio, levels = sort(unique(.data$odds_ratio))),
    baseline_label = paste0(round(.data$baseline_not_uti_prevalence * 100), "%"),
    detectable = if_else(.data$detectable_at_alpha_0_05, "Detectable", "Not detectable")
  ) %>%
  ggplot(aes(.data$odds_ratio, .data$baseline_label, fill = .data$detectable)) +
  geom_tile(colour = "white", linewidth = 0.55) +
  scale_fill_manual(values = c("Detectable" = "#2E8B57", "Not detectable" = "#E5E7EB")) +
  labs(
    title = "Sparse-count precision context",
    subtitle = "Expected-count illustration with 17 UTI VF-ready episodes",
    x = "Assumed odds ratio", y = "Not_UTI feature prevalence", fill = NULL
  ) +
  final_theme(9)
s3 <- p_s3a + p_s3b +
  plot_layout(widths = c(0.56, 0.44), guides = "collect") +
  plot_annotation(
    title = "Endpoint sensitivity and sparse-count limits",
    caption = "Near-miss rows are not relabelled; precision tiles are approximate interpretation context, not formal prospective power.",
    tag_levels = "A",
    theme = theme(plot.title = element_text(face = "bold", size = 15),
                  plot.caption = element_text(size = 9, colour = "grey35"))
  ) &
  theme(legend.position = "bottom")
register_figure(
  "near_miss_and_sparse_precision", "Supplementary", s3, 15, 8.7,
  "Current-rule near-miss audit and approximate sparse-count precision context.",
  required_paths[c("near_miss", "power_precision")],
  "Sensitivity diagnostic", "Near-miss episodes do not enter the primary UTI group; precision context is approximate."
)

# ------------------------------------------------------------------------------
# Supplementary Figure S4: Leave-one-UTI stability
# ------------------------------------------------------------------------------
score_features <- c("expec_marker_count", "upec_system_count", "total_vf_count_curated")
score_names <- c(
  expec_marker_count = "ExPEC-like markers",
  upec_system_count = "UPEC systems",
  total_vf_count_curated = "Curated VF genes"
)
s4_scores <- leave_one_detail %>%
  filter(.data$feature_type == "score", .data$statistic == "mean difference",
         .data$feature %in% score_features) %>%
  mutate(label = factor(score_names[.data$feature], levels = rev(unname(score_names))))
p_s4a <- ggplot(s4_scores, aes(.data$full_effect, .data$label)) +
  geom_vline(xintercept = 0, colour = "grey35") +
  geom_segment(aes(x = .data$min_without_one_uti, xend = .data$max_without_one_uti,
                   yend = .data$label), linewidth = 1.2, colour = "#4C78A8") +
  geom_point(size = 3, colour = "#0F766E") +
  labs(title = "Endpoint stability", subtitle = "Range after removing one UTI episode",
       x = "UTI minus Not_UTI mean difference", y = NULL) +
  final_theme(10)

s4_modules <- leave_one_detail %>%
  filter(.data$feature_type == "module", .data$statistic == "log2 odds ratio") %>%
  arrange(desc(.data$direction_flip), desc(.data$sensitivity_range)) %>%
  slice_head(n = 8) %>%
  mutate(
    label = .data$feature %>%
      str_remove("^mod_") %>%
      str_remove("_present$") %>%
      str_replace_all("_", " ") %>%
      str_to_sentence(),
    label = factor(.data$label, levels = rev(.data$label)),
    flip = if_else(.data$direction_flip, "Direction flip", "No direction flip")
  )
p_s4b <- ggplot(s4_modules, aes(.data$full_effect, .data$label, colour = .data$flip)) +
  geom_vline(xintercept = 0, colour = "grey35") +
  geom_segment(aes(x = .data$min_without_one_uti, xend = .data$max_without_one_uti,
                   yend = .data$label), linewidth = 1.1) +
  geom_point(size = 2.8) +
  scale_colour_manual(values = c("Direction flip" = "#C44E19", "No direction flip" = "#4C78A8")) +
  labs(title = "Most sensitive modules", subtitle = "Log2 odds-ratio range after removing one UTI",
       x = "Log2 odds ratio", y = NULL, colour = NULL) +
  final_theme(9)
s4 <- p_s4a + p_s4b +
  plot_layout(widths = c(0.42, 0.58), guides = "collect") +
  plot_annotation(
    title = "Leave-one-UTI-out stability",
    subtitle = "Twelve of 98 diagnostic effects change direction when one UTI episode is removed.",
    caption = "A direction flip indicates dependence on individual sparse UTI observations, not a corrected association result.",
    tag_levels = "A",
    theme = theme(plot.title = element_text(face = "bold", size = 15),
                  plot.subtitle = element_text(size = 10.5),
                  plot.caption = element_text(size = 9, colour = "grey35"))
  ) &
  theme(legend.position = "bottom")
register_figure(
  "leave_one_uti_stability", "Supplementary", s4, 15, 7,
  "Leave-one-UTI-out sensitivity for representative supplementary VF endpoints and the most unstable module results.",
  required_paths[c("leave_one_detail", "leave_one_summary")],
  "Sparse-case sensitivity diagnostic", "Instability is expected with 17 UTI cases and does not constitute association testing."
)

# ------------------------------------------------------------------------------
# Supplementary Figure S5: Prioritised variant map
# ------------------------------------------------------------------------------
variant_cases <- casebook %>%
  filter(.data$mechanism_bucket %in% c("same_strain_stable_profile",
                                       "same_strain_genomic_change", "uncertain")) %>%
  select("Participant_id", "from_tp", "to_tp", "case_label", "mechanism_bucket")
s5_variants <- variants %>%
  inner_join(variant_cases,
             by = c("Participant_id" = "Participant_id",
                    "From_Time" = "from_tp", "To_Time" = "to_tp")) %>%
  mutate(
    label = if_else(!is.na(.data$Gene) & nzchar(.data$Gene), .data$Gene, NA_character_),
    case_label = factor(.data$case_label, levels = rev(unique(as.character(variant_cases$case_label))))
  )
if (nrow(s5_variants) == 0) stop("No prioritised same-strain/uncertain detailed variants are available to plot.")
max_variant_pos <- max(s5_variants$Pos_Ref, na.rm = TRUE)
s5_labels <- s5_variants %>%
  filter(!is.na(.data$label)) %>%
  group_by(.data$case_label, .data$label) %>%
  slice_min(.data$Pos_Ref, n = 1, with_ties = FALSE) %>%
  ungroup()
s5_counts <- s5_variants %>%
  count(.data$case_label, name = "n_snps") %>%
  mutate(Pos_Ref = max_variant_pos * 1.02,
         label = paste0("n = ", .data$n_snps))
s5 <- ggplot(s5_variants, aes(.data$Pos_Ref, .data$case_label, colour = .data$mechanism_bucket)) +
  geom_segment(aes(x = 0, xend = max_variant_pos * 1.01,
                   yend = .data$case_label), colour = "grey88", linewidth = 1.2) +
  geom_point(size = 2.1, alpha = 0.85) +
  ggrepel::geom_text_repel(data = s5_labels, aes(label = .data$label),
                           size = 2.7, min.segment.length = 0, box.padding = 0.25,
                           point.padding = 0.15, max.overlaps = Inf, show.legend = FALSE) +
  geom_text(data = s5_counts, aes(x = .data$Pos_Ref, y = .data$case_label, label = .data$label),
            colour = "grey25",
            hjust = 0, size = 3, inherit.aes = FALSE) +
  scale_colour_manual(values = mechanism_cols, labels = mechanism_labels, name = "Mechanism category") +
  scale_x_continuous(labels = unit_format(unit = "Mb", scale = 1e-6),
                     expand = expansion(mult = c(0.01, 0.13))) +
  labs(
    title = "Prioritised within-strain variant map",
    subtitle = "Detailed annotations shown only for same-strain or uncertain Not_UTI-to-UTI cases with available GFF annotation",
    x = "Reference genome position", y = NULL,
    caption = "Gene labels identify annotated variant loci; unlabeled points are retained in counts. This is descriptive candidate prioritisation."
  ) +
  final_theme(10) +
  theme(legend.position = "bottom")
register_figure(
  "prioritised_variant_map", "Supplementary", s5, 13, 6.8,
  "Prioritised detailed-annotation SNP map for same-strain or uncertain Not_UTI-to-UTI cases with available annotations.",
  required_paths[c("casebook", "variants")],
  "Descriptive variant prioritisation", "Only annotated candidate comparisons are shown; absence from this panel is not absence of genomic change."
)

# ------------------------------------------------------------------------------
# Supplementary Figure S6: Lineage confounding diagnostic
# ------------------------------------------------------------------------------
status_st_data <- vf_ready %>%
  mutate(
    ST_norm = normalise_st_label(.data$ST),
    ST_plot = if_else(is.na(.data$ST_norm), "Missing ST", paste0("ST", .data$ST_norm)),
    UTI_Status = factor(.data$UTI_Status, levels = c("Not_UTI", "UTI"))
  )
top_st_groups <- status_st_data %>%
  count(.data$ST_plot, sort = TRUE) %>%
  filter(.data$ST_plot != "Missing ST") %>%
  slice_head(n = 8) %>%
  pull(.data$ST_plot)
st_comp <- status_st_data %>%
  mutate(
    ST_group_plot = case_when(
      .data$ST_plot %in% top_st_groups ~ .data$ST_plot,
      .data$ST_plot == "Missing ST" ~ "Missing ST",
      TRUE ~ "Other STs"
    )
  ) %>%
  count(.data$UTI_Status, .data$ST_group_plot, name = "n") %>%
  group_by(.data$UTI_Status) %>%
  mutate(prop = .data$n / sum(.data$n)) %>%
  ungroup()
base_st_cols <- c(
  "ST131" = "#4C78A8", "ST141" = "#F58518", "ST69" = "#54A24B", "ST73" = "#B279A2",
  "ST12" = "#72B7B2", "ST95" = "#E45756", "ST127" = "#9D755D", "ST10" = "#EECA3B",
  "Other STs" = "#BDBDBD", "Missing ST" = "#E5E7EB"
)
extra_st_groups <- setdiff(unique(st_comp$ST_group_plot), names(base_st_cols))
extra_st_cols <- if (length(extra_st_groups) > 0) {
  setNames(grDevices::hcl.colors(length(extra_st_groups), palette = "Dark 3"), extra_st_groups)
} else character()
st_cols <- c(base_st_cols, extra_st_cols)[unique(st_comp$ST_group_plot)]

p_s6a <- ggplot(st_comp, aes(.data$UTI_Status, .data$prop, fill = .data$ST_group_plot)) +
  geom_col(width = 0.66, colour = "white", linewidth = 0.25) +
  scale_y_continuous(labels = percent_format()) +
  scale_fill_manual(values = st_cols, name = "Sequence type") +
  labs(
    title = "Sequence-type composition",
    x = "Primary UTI status",
    y = "Within-status proportion of VF-ready isolates"
  ) +
  final_theme(9)

sts_both <- status_st_data %>%
  filter(!is.na(.data$ST_norm)) %>%
  count(.data$ST_norm, .data$UTI_Status) %>%
  group_by(.data$ST_norm) %>%
  filter(n_distinct(.data$UTI_Status) == 2) %>%
  pull(.data$ST_norm) %>%
  unique()

p_s6b <- ggplot(status_st_data %>% filter(.data$ST_norm %in% sts_both),
                aes(.data$UTI_Status, .data$total_vf_count_all, fill = .data$UTI_Status)) +
  geom_boxplot(outlier.shape = NA, width = 0.5, alpha = 0.75) +
  geom_jitter(width = 0.12, height = 0, alpha = 0.45, size = 1.3) +
  facet_wrap(~ paste0("ST", ST_norm), scales = "free_x") +
  scale_fill_uti_status(drop = FALSE) +
  labs(
    title = "Within-ST VF burden",
    x = "Primary UTI status",
    y = "Detected VF genes per isolate"
  ) +
  final_theme(9) +
  theme(legend.position = "none")

s6 <- p_s6a + p_s6b +
  plot_layout(widths = c(0.45, 0.55), guides = "collect") +
  plot_annotation(
    title = "Lineage confounding diagnostic",
    subtitle = "VF burden differs strongly by ST, but ST composition was not significantly different by primary status in this sparse UTI set.",
    caption = "Exploratory diagnostic; UTI n=17 and most STs contain too few UTI isolates for within-ST inference.",
    tag_levels = "A",
    theme = theme(plot.title = element_text(face = "bold", size = 15),
                  plot.subtitle = element_text(size = 10.5),
                  plot.caption = element_text(size = 9, colour = "grey35"))
  ) &
  theme(legend.position = "bottom")
register_figure(
  "lineage_confounding_diagnostic", "Supplementary", s6, 13, 6.5,
  "Lineage confounding diagnostic combining sequence-type composition by primary status with within-ST VF burden for STs containing both UTI and Not_UTI isolates.",
  c(required_paths[c("vf_ready", "stat_validation")]),
  "Exploratory lineage diagnostic", "ST structure is sparse for UTI episodes; within-ST contrasts are underpowered."
)

# ------------------------------------------------------------------------------
# Supplementary Figure S7: Paired resident ExPEC-like marker slopeplot
# ------------------------------------------------------------------------------
s7_data <- stat_score_values %>%
  filter(.data$score == "expec_marker_count") %>%
  group_by(.data$Participant_id) %>%
  filter(all(c("Not_UTI", "UTI") %in% as.character(.data$UTI_Status))) %>%
  ungroup() %>%
  mutate(UTI_Status = factor(.data$UTI_Status, levels = c("Not_UTI", "UTI")))

s7 <- ggplot(s7_data, aes(.data$UTI_Status, .data$participant_status_median, group = .data$Participant_id)) +
  geom_line(colour = "#B8B8B8", linewidth = 0.45, alpha = 0.85) +
  geom_point(aes(colour = .data$UTI_Status), size = 2.2) +
  scale_colour_uti_status(drop = FALSE) +
  labs(
    title = "Paired resident ExPEC-like marker count",
    subtitle = "Each line is one resident; this visual reduces between-resident confounding but includes only residents observed in both states.",
    x = "Primary status within resident",
    y = "Participant-median ExPEC-like marker count",
    colour = "Primary status",
    caption = "Paired resident-level descriptive sensitivity; residents observed in only one state are not represented."
  ) +
  final_theme(10)
register_figure(
  "paired_resident_expec_marker", "Supplementary", s7, 7.5, 5.8,
  "Each line is one resident; this visual reduces between-resident confounding but includes only residents observed in both states.",
  required_paths[c("stat_score_values")],
  "Paired descriptive sensitivity", "Only residents with both statuses are shown; this does not represent residents observed in one state only."
)

# ------------------------------------------------------------------------------
# Supplementary Figure S8: Not_UTI-to-UTI module gain/loss heatmap
# ------------------------------------------------------------------------------
s8_plot_data <- stat_module_matrix %>%
  mutate(
    case_label = factor(.data$case_label, levels = unique(.data$case_label)),
    module_label = factor(.data$module_label, levels = rev(unique(.data$module_label))),
    change_display = factor(.data$change_display,
                            levels = c("Gained", "Lost", "Stable present", "Stable absent"))
  )

s8 <- ggplot(s8_plot_data, aes(.data$case_label, .data$module_label, fill = .data$change_display)) +
  geom_tile(colour = "white", linewidth = 0.25) +
  scale_fill_manual(
    values = c("Gained" = "#009E73", "Lost" = "#CC79A7",
               "Stable present" = "#0072B2", "Stable absent" = "#F2F2F2"),
    name = "Module change"
  ) +
  labs(
    title = "Not_UTI-to-UTI module gain/loss heatmap",
    subtitle = "Module changes across WGS-linked Not_UTI-to-UTI transitions; descriptive, not causal.",
    x = "Transition case",
    y = "VF module",
    caption = "Transitions can repeat within residents; stable/present and stable/absent states are shown to keep absence explicit."
  ) +
  final_theme(8.5) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
register_figure(
  "not_uti_to_uti_module_gain_loss", "Supplementary", s8, 12.5, 9.5,
  "Module gain/loss heatmap across WGS-linked Not_UTI-to-UTI transitions.",
  required_paths[c("stat_module_matrix", "stat_transition_module")],
  "Descriptive transition sensitivity", "Module changes are presence/absence observations and do not prove causality."
)

# ------------------------------------------------------------------------------
# Supplementary Figure S9: VF module PCoA
# ------------------------------------------------------------------------------
stat_pcoa <- stat_pcoa %>%
  mutate(
    ST_group = factor(st_group4(.data$ST), levels = c("ST131", "ST141", "Other typed ST", "Missing ST"))
  ) %>%
  filter(.data$UTI_Status %in% c("Not_UTI", "UTI"))
axis1_lab <- sprintf("PCoA1 (%.1f%% variation)", unique(stat_pcoa$var_Axis1)[1])
axis2_lab <- sprintf("PCoA2 (%.1f%% variation)", unique(stat_pcoa$var_Axis2)[1])
s9 <- ggplot(stat_pcoa, aes(.data$Axis1, .data$Axis2, colour = .data$UTI_Status, shape = .data$ST_group)) +
  geom_point(size = 2.2, alpha = 0.78) +
  scale_colour_uti_status(drop = FALSE) +
  scale_shape_manual(values = c("ST131" = 16, "ST141" = 17, "Other typed ST" = 1, "Missing ST" = 4)) +
  labs(
    title = "Global VF module PCoA",
    subtitle = "Global VF module profiles show whether UTI episodes separate from Not_UTI episodes; interpretation is limited by sparse UTI counts and lineage structure.",
    x = axis1_lab,
    y = axis2_lab,
    colour = "Primary status",
    shape = "Sequence type group",
    caption = "Jaccard PCoA from module presence/absence; descriptive and unadjusted for repeated residents or lineage."
  ) +
  final_theme(10)
register_figure(
  "vf_module_pcoa_primary_status", "Supplementary", s9, 8.2, 6.2,
  "Global VF module profiles show whether UTI episodes separate from Not_UTI episodes; interpretation is limited by sparse UTI counts and lineage structure.",
  required_paths[c("stat_pcoa", "stat_validation")],
  "Exploratory global profile visualisation", "Separation patterns are descriptive and are not adjusted for repeated residents or lineage."
)

write_csv(manifest, file.path(DIR_FINAL_RESULTS, "final_figure_manifest.csv"))

caption_lines <- c(
  "# Final Figure Captions",
  "",
  sprintf("Generated: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  "",
  "All figures use the current primary `UTI_Status` definition. Analyses are descriptive or sensitivity-focused because the VF/model denominator includes only 17 primary UTI episodes.",
  ""
)
for (i in seq_len(nrow(manifest))) {
  caption_lines <- c(
    caption_lines,
    paste0("## ", manifest$figure_id[i]),
    "",
    manifest$caption[i],
    "",
    paste0("Evidence type: ", manifest$evidence_type[i], "."),
    paste0("Interpretation limitation: ", manifest$interpretation_limitations[i]),
    ""
  )
}
caption_lines <- c(
  caption_lines,
  "## Reproducibility Note",
  "",
  "Run `Rscript 35_final_figure_pack.R` after the mechanism and robustness add-ons have passed validation. Existing publication figures are not overwritten.",
  "",
  "The clinical-to-VF retention panel is calculated directly from primary episode keys joined to `vf_analysis_ready.csv`; it is not taken from a pre-exclusion QC summary table."
)
writeLines(caption_lines, file.path(DIR_FINAL_RESULTS, "final_figure_captions.md"))

msg("Rendered %d final figures as PNG and PDF.", nrow(manifest))
msg("Outputs written to %s and %s", DIR_FINAL_RESULTS, DIR_FINAL_PLOTS)
msg("35_final_figure_pack.R complete.")

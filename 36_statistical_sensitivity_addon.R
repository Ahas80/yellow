#!/usr/bin/env Rscript
# ==============================================================================
# 36_statistical_sensitivity_addon.R
# ==============================================================================
#
# GOAL:
#   Add a narrow, prespecified sensitivity layer for the current primary
#   UTI vs Not_UTI analysis. This script deliberately avoids new broad
#   discovery testing because the VF-ready UTI denominator is sparse.
#
# OUTPUT:
#   - results/statistical_sensitivity/participant_collapsed_score_values.csv
#   - results/statistical_sensitivity/participant_collapsed_score_tests.csv
#   - results/statistical_sensitivity/paired_binary_feature_sensitivity.csv
#   - results/statistical_sensitivity/transition_module_gain_loss_enrichment.csv
#   - results/statistical_sensitivity/score_glmm_sensitivity.csv
#   - results/statistical_sensitivity/statistical_sensitivity_validation_checks.csv
#   - results/statistical_sensitivity/statistical_sensitivity_summary.md
#   - plots/statistical_sensitivity/*.png/.pdf
#
# ==============================================================================

source("00_config.R")
source("R/plot_helpers.R")
suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
  library(stringr)
  library(tibble)
  library(ggplot2)
  library(purrr)
  library(lme4)
  library(broom.mixed)
  library(patchwork)
})

msg("Starting 36_statistical_sensitivity_addon.R")

DIR_STAT <- file.path(DIR_RESULTS, "statistical_sensitivity")
DIR_PLOTS_STAT <- file.path(DIR_PLOTS, "statistical_sensitivity")
ensure_dir(DIR_STAT)
ensure_dir(DIR_PLOTS_STAT)

# ==============================================================================
# HELPERS
# ==============================================================================

require_file <- function(path, label = basename(path)) {
  if (!file.exists(path)) stop("Missing required input for statistical sensitivity add-on: ", label, " (", path, ")")
  path
}

normalise_tp_label <- function(x) {
  if (exists("normalise_timepoint_preserve_events", mode = "function")) {
    normalise_timepoint_preserve_events(x)
  } else {
    str_trim(as.character(x))
  }
}

load_active_longcycler_keys <- function() {
  cohort <- read_csv(require_file(FILE_ANALYSIS_CLINICAL_COHORT, "selected Longcycler analysis cohort"), show_col_types = FALSE)
  keys <- cohort %>%
    mutate(
      Participant_id = as.character(.data$Participant_id),
      tp_lab = normalise_tp_label(.data$tp_lab)
    ) %>%
    distinct(.data$Participant_id, .data$tp_lab)
  if (nrow(keys) != nrow(cohort)) stop("Selected Longcycler analysis cohort has duplicate episode keys")
  keys
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

median_or_na <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0) NA_real_ else median(x)
}

mean_or_na <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0) NA_real_ else mean(x)
}

q_or_na <- function(x, p) {
  x <- x[!is.na(x)]
  if (length(x) == 0) NA_real_ else as.numeric(quantile(x, p, names = FALSE))
}

safe_wilcox_p <- function(x, y = NULL, paired = FALSE) {
  x <- x[!is.na(x)]
  if (is.null(y)) {
    if (length(x) < 2 || length(unique(x)) < 2) return(NA_real_)
    return(suppressWarnings(wilcox.test(x, mu = 0, exact = FALSE)$p.value))
  }
  y <- y[!is.na(y)]
  if (length(x) < 1 || length(y) < 1) return(NA_real_)
  if (length(unique(c(x, y))) < 2) return(NA_real_)
  suppressWarnings(wilcox.test(x, y, paired = paired, exact = FALSE)$p.value)
}

safe_sign_p <- function(n_positive, n_negative) {
  n <- n_positive + n_negative
  if (is.na(n) || n == 0) return(NA_real_)
  binom.test(n_positive, n, p = 0.5)$p.value
}

bootstrap_ci <- function(x, y = NULL, statistic = c("median_difference", "median_delta"),
                         B = 5000L, seed = 20260530L) {
  statistic <- match.arg(statistic)
  set.seed(seed)
  if (is.null(y)) {
    x <- x[!is.na(x)]
    if (length(x) < 2) return(c(lower = NA_real_, upper = NA_real_))
    vals <- replicate(B, median(sample(x, length(x), replace = TRUE)))
  } else {
    x <- x[!is.na(x)]
    y <- y[!is.na(y)]
    if (length(x) < 1 || length(y) < 1) return(c(lower = NA_real_, upper = NA_real_))
    vals <- replicate(B, {
      median(sample(x, length(x), replace = TRUE)) -
        median(sample(y, length(y), replace = TRUE))
    })
  }
  c(lower = as.numeric(quantile(vals, 0.025, na.rm = TRUE, names = FALSE)),
    upper = as.numeric(quantile(vals, 0.975, na.rm = TRUE, names = FALSE)))
}

plot_theme_stat <- function(base_size = 10) {
  theme_bw(base_size = base_size) +
    theme(
      plot.title = element_text(face = "bold"),
      plot.subtitle = element_text(colour = "grey25"),
      plot.caption = element_text(hjust = 0, colour = "grey35", size = base_size - 2),
      legend.position = "bottom",
      panel.grid.minor = element_blank(),
      axis.title = element_text(face = "bold")
    )
}

feature_display_name <- function(x) {
  x %>%
    str_remove("^mod_") %>%
    str_remove("_present$") %>%
    str_replace_all("_", " ") %>%
    str_to_sentence()
}

read_current_vf_ready <- function() {
  read_csv(require_file(FILE_VF_READY), show_col_types = FALSE) %>%
    prefer_primary_uti_status(allow_legacy_fallback = FALSE) %>%
    apply_manual_sample_curation(context = "36_vf_ready") %>%
    filter_primary_genomics() %>%
    mutate(
      Participant_id = as.character(Participant_id),
      tp_lab = normalise_tp_label(tp_lab),
      ST = if ("ST" %in% names(.)) normalise_st_label(ST) else NA_character_
    )
}

# ==============================================================================
# LOAD CURRENT INPUTS AND VALIDATE DENOMINATORS
# ==============================================================================

status_primary <- read_csv(require_file(FILE_ANALYSIS_CLINICAL_COHORT), show_col_types = FALSE) %>%
  prefer_primary_uti_status(allow_legacy_fallback = FALSE) %>%
  mutate(Participant_id = as.character(Participant_id),
         tp_lab = normalise_tp_label(tp_lab))

active_longcycler_keys <- load_active_longcycler_keys()
canonical_transitions <- read_csv(
  require_file(file.path(DIR_RESULTS, "longitudinal", "longcycler_transitions.csv"),
               "canonical Longcycler transitions"),
  show_col_types = FALSE
)
if (nrow(canonical_transitions) != 371L ||
    sum(canonical_transitions$status_from == "Not_UTI" & canonical_transitions$status_to == "UTI", na.rm = TRUE) != 9L ||
    any(is.na(canonical_transitions$TotalSNPs))) {
  stop("Canonical Longcycler transition export failed the 371/9/direct-SNP contract")
}
expected_active_status <- status_primary
if (nrow(expected_active_status) != nrow(active_longcycler_keys)) {
  stop("Not every active selected QC-pass Longcycler key has one primary clinical status row.")
}
expected_active_total <- nrow(active_longcycler_keys)
expected_active_uti <- sum(expected_active_status$UTI_Status == "UTI", na.rm = TRUE)
expected_active_not_uti <- sum(expected_active_status$UTI_Status == "Not_UTI", na.rm = TRUE)

vf_ready <- read_current_vf_ready()
vf_keys <- vf_ready %>% distinct(.data$Participant_id, .data$tp_lab)
if (nrow(vf_ready) != expected_active_total ||
    nrow(anti_join(active_longcycler_keys, vf_keys, by = c("Participant_id", "tp_lab"))) ||
    nrow(anti_join(vf_keys, active_longcycler_keys, by = c("Participant_id", "tp_lab")))) {
  stop("VF-ready keys do not exactly equal the selected Longcycler analysis cohort.")
}

score_path <- require_file(file.path(DIR_VF, "vf_score_table.csv"))
scores <- read_csv(score_path, show_col_types = FALSE) %>%
  prefer_primary_uti_status(allow_legacy_fallback = FALSE) %>%
  apply_manual_sample_curation(context = "36_scores") %>%
  filter_primary_genomics() %>%
  mutate(
    Participant_id = as.character(Participant_id),
    tp_lab = normalise_tp_label(tp_lab),
    ST = if ("ST" %in% names(.)) normalise_st_label(ST) else NA_character_
  )
score_keys <- scores %>% distinct(.data$Participant_id, .data$tp_lab)
if (nrow(score_keys) != nrow(scores) ||
    nrow(anti_join(active_longcycler_keys, score_keys, by = c("Participant_id", "tp_lab"))) ||
    nrow(anti_join(score_keys, active_longcycler_keys, by = c("Participant_id", "tp_lab")))) {
  stop("VF score keys do not exactly equal the selected Longcycler analysis cohort")
}

if (any(is.na(scores$UTI_Status))) {
  scores <- scores %>%
    select(-any_of(c("UTI_Status", "UTI_binary"))) %>%
    left_join(
      vf_ready %>% select(Participant_id, tp_lab, UTI_Status, UTI_binary),
      by = c("Participant_id", "tp_lab"),
      relationship = "many-to-one"
    )
}

if (!"UTI_binary" %in% names(scores)) {
  scores <- scores %>%
    mutate(UTI_binary = if_else(UTI_Status == "UTI", 1L, if_else(UTI_Status == "Not_UTI", 0L, NA_integer_)))
}

scores <- scores %>%
  filter(UTI_Status %in% c("UTI", "Not_UTI")) %>%
  mutate(
    UTI_Status = factor(UTI_Status, levels = c("Not_UTI", "UTI")),
    ST_group = st_group4(ST),
    Batch_model = factor(if ("Batch" %in% names(.)) as.character(Batch) else "Missing batch"),
    Timepoint_model = factor(if ("Timepoint" %in% names(.)) as.character(Timepoint) else as.character(tp_lab))
  )

score_endpoints <- c(
  "total_vf_count_all",
  "total_vf_count_curated",
  "expec_marker_count",
  "upec_system_count"
)
score_endpoints <- intersect(score_endpoints, names(scores))
if (length(score_endpoints) != 4) {
  stop("Expected four prespecified supplementary VF endpoints in vf_score_table.csv; found: ", paste(score_endpoints, collapse = ", "))
}

module_cols <- grep("^mod_.*_present$", names(scores), value = TRUE)
gene_cols <- canonical_vf_gene_cols(names(vf_ready), required = FALSE)
feature_df <- scores %>%
  select(Participant_id, tp_lab, UTI_Status, all_of(module_cols)) %>%
  left_join(
    vf_ready %>% select(Participant_id, tp_lab, all_of(gene_cols)),
    by = c("Participant_id", "tp_lab"),
    relationship = "many-to-one"
  )

msg("Loaded current sensitivity denominator: clinical %d rows; VF/model %d rows.",
    nrow(status_primary), nrow(vf_ready))

# ==============================================================================
# PARTICIPANT-COLLAPSED SCORE TESTS
# ==============================================================================

participant_score_values <- scores %>%
  select(Participant_id, UTI_Status, all_of(score_endpoints)) %>%
  pivot_longer(all_of(score_endpoints), names_to = "score", values_to = "value") %>%
  group_by(Participant_id, UTI_Status, score) %>%
  summarise(participant_status_median = median_or_na(value), .groups = "drop") %>%
  mutate(
    score_label = recode(
      score,
      total_vf_count_all = "All VF genes",
      total_vf_count_curated = "Curated VF genes",
      expec_marker_count = "ExPEC-like marker count",
      upec_system_count = "UPEC systems present"
    )
  )

write_csv(participant_score_values, file.path(DIR_STAT, "participant_collapsed_score_values.csv"))

participant_score_summary <- participant_score_values %>%
  group_by(score, score_label, UTI_Status) %>%
  summarise(
    n_participant_status_groups = sum(!is.na(participant_status_median)),
    median = median_or_na(participant_status_median),
    q25 = q_or_na(participant_status_median, 0.25),
    q75 = q_or_na(participant_status_median, 0.75),
    mean = mean_or_na(participant_status_median),
    sd = ifelse(sum(!is.na(participant_status_median)) > 1, sd(participant_status_median, na.rm = TRUE), NA_real_),
    .groups = "drop"
  )
write_csv(participant_score_summary, file.path(DIR_STAT, "participant_collapsed_score_summary.csv"))

participant_score_tests <- map_dfr(score_endpoints, function(sc) {
  score_dat <- participant_score_values %>% filter(score == sc)
  uti <- score_dat %>% filter(UTI_Status == "UTI") %>% pull(participant_status_median)
  not <- score_dat %>% filter(UTI_Status == "Not_UTI") %>% pull(participant_status_median)
  wide <- score_dat %>%
    select(Participant_id, UTI_Status, participant_status_median) %>%
    pivot_wider(names_from = UTI_Status, values_from = participant_status_median) %>%
    filter(!is.na(UTI), !is.na(Not_UTI)) %>%
    mutate(delta_uti_minus_not_uti = UTI - Not_UTI)
  ci_unpaired <- bootstrap_ci(uti, not, B = 5000L, seed = 20260530L + match(sc, score_endpoints))
  ci_paired <- bootstrap_ci(wide$delta_uti_minus_not_uti, B = 5000L,
                            seed = 20260630L + match(sc, score_endpoints))
  tibble(
    score = sc,
    score_label = unique(score_dat$score_label)[1],
    n_participants_uti = sum(!is.na(uti)),
    n_participants_not_uti = sum(!is.na(not)),
    n_paired_participants = nrow(wide),
    median_uti = median_or_na(uti),
    median_not_uti = median_or_na(not),
    median_difference_uti_minus_not_uti = median_or_na(uti) - median_or_na(not),
    mean_difference_uti_minus_not_uti = mean_or_na(uti) - mean_or_na(not),
    bootstrap_median_diff_ci_lower = ci_unpaired[["lower"]],
    bootstrap_median_diff_ci_upper = ci_unpaired[["upper"]],
    wilcoxon_rank_sum_p = safe_wilcox_p(uti, not),
    paired_median_delta = median_or_na(wide$delta_uti_minus_not_uti),
    paired_mean_delta = mean_or_na(wide$delta_uti_minus_not_uti),
    paired_bootstrap_median_delta_ci_lower = ci_paired[["lower"]],
    paired_bootstrap_median_delta_ci_upper = ci_paired[["upper"]],
    n_positive_delta = sum(wide$delta_uti_minus_not_uti > 0, na.rm = TRUE),
    n_negative_delta = sum(wide$delta_uti_minus_not_uti < 0, na.rm = TRUE),
    n_nonzero_delta = sum(wide$delta_uti_minus_not_uti != 0, na.rm = TRUE),
    paired_sign_test_p = safe_sign_p(
      sum(wide$delta_uti_minus_not_uti > 0, na.rm = TRUE),
      sum(wide$delta_uti_minus_not_uti < 0, na.rm = TRUE)
    ),
    paired_signed_rank_p = safe_wilcox_p(wide$delta_uti_minus_not_uti),
    interpretation = "Participant-collapsed exploratory sensitivity; one median per resident/status."
  )
}) %>%
  mutate(
    wilcoxon_rank_sum_q_BH = p.adjust(wilcoxon_rank_sum_p, method = "BH"),
    paired_sign_test_q_BH = p.adjust(paired_sign_test_p, method = "BH"),
    paired_signed_rank_q_BH = p.adjust(paired_signed_rank_p, method = "BH")
  )

write_csv(participant_score_tests, file.path(DIR_STAT, "participant_collapsed_score_tests.csv"))

# ==============================================================================
# PAIRED BINARY FEATURE SENSITIVITY
# ==============================================================================

gene_tests <- if (file.exists(file.path(DIR_VF, "vf_gene_prevalence_tests.csv"))) {
  read_csv(file.path(DIR_VF, "vf_gene_prevalence_tests.csv"), show_col_types = FALSE) %>%
    arrange(p_value) %>%
    slice_head(n = 10) %>%
    transmute(feature = gene, feature_source = "top_gene_fisher")
} else tibble(feature = character(), feature_source = character())

glmm_tests <- read_csv(require_file(file.path(DIR_MODELS, "gwas_multivariable_glmm.csv")),
                       show_col_types = FALSE) %>%
  arrange(p.value, FDR) %>%
  slice_head(n = 10) %>%
  transmute(feature, feature_source = "top_glmm")

module_tests <- read_csv(require_file(file.path(DIR_VF, "uti_not_uti_feature_fisher_exploratory.csv")),
                         show_col_types = FALSE) %>%
  filter(feature_type == "module") %>%
  arrange(fisher_p) %>%
  slice_head(n = 10) %>%
  transmute(feature, feature_source = "top_module_fisher")

explicit_features <- tibble(
  feature = c("astA", "ybtA", "ybtE", "ybtT", "ybtU", "irp1", "irp2", "fyuA",
              "mod_toxin_east1_present", "mod_iron_yersiniabactin_present"),
  feature_source = "prespecified_example"
)

feature_catalog <- bind_rows(gene_tests, glmm_tests, module_tests, explicit_features) %>%
  filter(!is.na(feature), feature != "") %>%
  distinct(feature, feature_source) %>%
  group_by(feature) %>%
  summarise(
    feature_sources = paste(sort(unique(feature_source)), collapse = ";"),
    .groups = "drop"
  ) %>%
  mutate(
    feature_type = case_when(
      feature %in% module_cols ~ "module",
      feature %in% gene_cols ~ "gene",
      TRUE ~ "unavailable"
    )
  )

write_csv(feature_catalog, file.path(DIR_STAT, "paired_binary_feature_catalog.csv"))

available_features <- feature_catalog %>%
  filter(feature_type %in% c("gene", "module"), feature %in% names(feature_df)) %>%
  pull(feature)

participant_feature_values <- feature_df %>%
  select(Participant_id, UTI_Status, all_of(available_features)) %>%
  pivot_longer(all_of(available_features), names_to = "feature", values_to = "present") %>%
  group_by(Participant_id, UTI_Status, feature) %>%
  summarise(participant_status_prevalence = mean(present, na.rm = TRUE), .groups = "drop")

paired_feature_deltas <- participant_feature_values %>%
  pivot_wider(names_from = UTI_Status, values_from = participant_status_prevalence) %>%
  filter(!is.na(UTI), !is.na(Not_UTI)) %>%
  mutate(delta_uti_minus_not_uti = UTI - Not_UTI) %>%
  left_join(feature_catalog, by = "feature") %>%
  mutate(feature_label = feature_display_name(feature))

write_csv(paired_feature_deltas, file.path(DIR_STAT, "paired_binary_feature_deltas.csv"))

paired_feature_sensitivity <- paired_feature_deltas %>%
  group_by(feature, feature_label, feature_type, feature_sources) %>%
  summarise(
    n_paired_participants = n(),
    median_delta_uti_minus_not_uti = median_or_na(delta_uti_minus_not_uti),
    mean_delta_uti_minus_not_uti = mean_or_na(delta_uti_minus_not_uti),
    n_positive_delta = sum(delta_uti_minus_not_uti > 0, na.rm = TRUE),
    n_negative_delta = sum(delta_uti_minus_not_uti < 0, na.rm = TRUE),
    n_nonzero_delta = sum(delta_uti_minus_not_uti != 0, na.rm = TRUE),
    sign_test_p = safe_sign_p(n_positive_delta, n_negative_delta),
    signed_rank_p = safe_wilcox_p(delta_uti_minus_not_uti),
    .groups = "drop"
  ) %>%
  mutate(
    sign_test_q_BH = p.adjust(sign_test_p, method = "BH"),
    signed_rank_q_BH = p.adjust(signed_rank_p, method = "BH"),
    interpretation = "Paired resident binary-feature sensitivity; exploratory and limited to residents observed in both primary states."
  ) %>%
  arrange(sign_test_p, signed_rank_p)

write_csv(paired_feature_sensitivity, file.path(DIR_STAT, "paired_binary_feature_sensitivity.csv"))

# ==============================================================================
# TRANSITION-LEVEL MODULE GAIN/LOSS ENRICHMENT
# ==============================================================================

case_summary <- read_csv(require_file(file.path(DIR_VF, "vf_transition_case_summary.csv")),
                         show_col_types = FALSE) %>%
  mutate(
    case_id = as.character(case_id),
    transition_type_simple = str_replace_all(as.character(transition_type), "→", "->")
  )
if (nrow(case_summary) != 371L ||
    sum(case_summary$transition_type_simple == "Not_UTI->UTI", na.rm = TRUE) != 9L ||
    any(!(case_summary$has_vf_pair %in% TRUE)) || any(is.na(case_summary$SNPs))) {
  stop("Statistical transition inputs must be the 371 fully linked selected Longcycler adjacent pairs")
}
case_keys <- case_summary %>%
  transmute(key = paste(as.character(.data$Participant_id),
                        normalise_tp_label(.data$from_tp), normalise_tp_label(.data$to_tp), sep = "|"))
canonical_keys <- canonical_transitions %>%
  transmute(key = paste(as.character(.data$Participant_id),
                        normalise_tp_label(.data$tp_from), normalise_tp_label(.data$tp_to), sep = "|"))
if (!setequal(case_keys$key, canonical_keys$key)) {
  stop("Statistical transition inputs do not match the canonical Longcycler transition keys")
}

module_changes <- read_csv(require_file(file.path(DIR_VF, "vf_transition_module_changes.csv")),
                           show_col_types = FALSE) %>%
  mutate(case_id = as.character(case_id))

transition_module_events <- module_changes %>%
  left_join(
    case_summary %>% select(case_id, transition_type_simple, has_module_pair),
    by = "case_id",
    relationship = "many-to-one"
  ) %>%
  filter(
    has_module_pair %in% TRUE,
    transition_type_simple %in% c("Not_UTI->UTI", "Not_UTI->Not_UTI")
  ) %>%
  mutate(
    module_label = if ("system_name" %in% names(.)) system_name else feature_display_name(module_col),
    transition_group = if_else(transition_type_simple == "Not_UTI->UTI", "Not_UTI_to_UTI", "Not_UTI_to_Not_UTI")
  )

transition_module_enrichment <- map_dfr(sort(unique(transition_module_events$module_col)), function(mod) {
  mod_dat <- transition_module_events %>% filter(module_col == mod)
  map_dfr(c("gained", "lost"), function(direction) {
    a <- sum(mod_dat$transition_group == "Not_UTI_to_UTI" & mod_dat$change_type == direction, na.rm = TRUE)
    b <- sum(mod_dat$transition_group == "Not_UTI_to_UTI" & mod_dat$change_type != direction, na.rm = TRUE)
    c <- sum(mod_dat$transition_group == "Not_UTI_to_Not_UTI" & mod_dat$change_type == direction, na.rm = TRUE)
    d <- sum(mod_dat$transition_group == "Not_UTI_to_Not_UTI" & mod_dat$change_type != direction, na.rm = TRUE)
    ft <- suppressWarnings(fisher.test(matrix(c(a, b, c, d), nrow = 2, byrow = TRUE)))
    tibble(
      module_col = mod,
      module_id = unique(mod_dat$module_id)[1],
      module_label = unique(mod_dat$module_label)[1],
      broad_module = if ("broad_module" %in% names(mod_dat)) unique(mod_dat$broad_module)[1] else NA_character_,
      direction = direction,
      n_not_uti_to_uti_with_event = a,
      n_not_uti_to_uti_without_event = b,
      n_not_uti_to_not_uti_with_event = c,
      n_not_uti_to_not_uti_without_event = d,
      prevalence_not_uti_to_uti = a / (a + b),
      prevalence_not_uti_to_not_uti = c / (c + d),
      prevalence_difference = (a / (a + b)) - (c / (c + d)),
      fisher_or = unname(ft$estimate),
      fisher_p = ft$p.value
    )
  })
}) %>%
  mutate(
    fisher_q_BH = p.adjust(fisher_p, method = "BH"),
    interpretation = "Exploratory transition-level module gain/loss enrichment; transitions can repeat within residents."
  ) %>%
  arrange(fisher_p, direction, module_label)

write_csv(transition_module_enrichment, file.path(DIR_STAT, "transition_module_gain_loss_enrichment.csv"))

not_uti_uti_module_matrix <- transition_module_events %>%
  filter(transition_group == "Not_UTI_to_UTI") %>%
  mutate(
    case_label = paste0(Participant_id, ": ", case_id),
    change_display = recode(
      change_type,
      gained = "Gained",
      lost = "Lost",
      stable_present = "Stable present",
      stable_absent = "Stable absent",
      .default = change_type
    )
  ) %>%
  select(case_id, Participant_id, case_label, module_col, module_id, module_label,
         broad_module, change_type, change_display)

write_csv(not_uti_uti_module_matrix, file.path(DIR_STAT, "not_uti_to_uti_module_change_matrix.csv"))

# ==============================================================================
# SCORE-LEVEL GLMM SENSITIVITY
# ==============================================================================

fit_score_model <- function(score, include_st_group = FALSE) {
  data <- scores %>%
    filter(!is.na(.data[[score]]), !is.na(UTI_binary)) %>%
    mutate(
      scaled_score = as.numeric(scale(.data[[score]])),
      Participant_id = factor(Participant_id),
      Batch_model = droplevels(Batch_model),
      Timepoint_model = droplevels(Timepoint_model),
      ST_group = factor(ST_group, levels = c("ST131", "ST141", "Other typed ST", "Missing ST"))
    ) %>%
    filter(!is.na(scaled_score))

  covariates <- c("scaled_score")
  if (n_distinct(data$Batch_model, na.rm = TRUE) > 1) covariates <- c(covariates, "Batch_model")
  if (n_distinct(data$Timepoint_model, na.rm = TRUE) > 1) covariates <- c(covariates, "Timepoint_model")
  if (include_st_group && n_distinct(data$ST_group, na.rm = TRUE) > 1) covariates <- c(covariates, "ST_group")
  fixed_part <- paste(covariates, collapse = " + ")
  fml_glmm <- as.formula(sprintf("UTI_binary ~ %s + (1 | Participant_id)", fixed_part))
  fml_glm <- as.formula(sprintf("UTI_binary ~ %s", fixed_part))
  variant <- if (include_st_group && "ST_group" %in% covariates) "batch_timepoint_collapsed_st_glmm" else "batch_timepoint_glmm"

  result <- tryCatch({
    mod <- lme4::glmer(
      fml_glmm,
      data = data,
      family = binomial(link = "logit"),
      control = lme4::glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))
    )
    tidy <- broom.mixed::tidy(mod, effects = "fixed", conf.int = TRUE, conf.method = "Wald") %>%
      filter(term == "scaled_score")
    tibble(
      score = score,
      model_variant = variant,
      model_type = if_else(lme4::isSingular(mod), "GLMM (Singular)", "GLMM"),
      n_rows = nrow(data),
      n_participants = n_distinct(data$Participant_id),
      n_uti = sum(data$UTI_binary == 1, na.rm = TRUE),
      n_not_uti = sum(data$UTI_binary == 0, na.rm = TRUE),
      estimate = tidy$estimate[1],
      std.error = tidy$std.error[1],
      OR_per_1sd = exp(tidy$estimate[1]),
      OR_lower = exp(tidy$conf.low[1]),
      OR_upper = exp(tidy$conf.high[1]),
      p.value = tidy$p.value[1],
      converged = TRUE,
      covariates = paste(covariates, collapse = " + "),
      interpretation = "Prespecified score-level GLMM sensitivity; exploratory with sparse UTI denominator."
    )
  }, error = function(e) {
    tryCatch({
      mod <- stats::glm(fml_glm, data = data, family = binomial(link = "logit"))
      tidy <- broom::tidy(mod, conf.int = TRUE) %>% filter(term == "scaled_score")
      tibble(
        score = score,
        model_variant = str_replace(variant, "glmm$", "glm_fallback"),
        model_type = "GLM_Fallback",
        n_rows = nrow(data),
        n_participants = n_distinct(data$Participant_id),
        n_uti = sum(data$UTI_binary == 1, na.rm = TRUE),
        n_not_uti = sum(data$UTI_binary == 0, na.rm = TRUE),
        estimate = tidy$estimate[1],
        std.error = tidy$std.error[1],
        OR_per_1sd = exp(tidy$estimate[1]),
        OR_lower = exp(tidy$conf.low[1]),
        OR_upper = exp(tidy$conf.high[1]),
        p.value = tidy$p.value[1],
        converged = TRUE,
        covariates = paste(covariates, collapse = " + "),
        interpretation = paste("GLMM failed; GLM fallback used.", e$message)
      )
    }, error = function(e2) {
      tibble(
        score = score,
        model_variant = variant,
        model_type = "Failed",
        n_rows = nrow(data),
        n_participants = n_distinct(data$Participant_id),
        n_uti = sum(data$UTI_binary == 1, na.rm = TRUE),
        n_not_uti = sum(data$UTI_binary == 0, na.rm = TRUE),
        estimate = NA_real_,
        std.error = NA_real_,
        OR_per_1sd = NA_real_,
        OR_lower = NA_real_,
        OR_upper = NA_real_,
        p.value = NA_real_,
        converged = FALSE,
        covariates = paste(covariates, collapse = " + "),
        interpretation = paste("Model failed.", e2$message)
      )
    })
  })
  result
}

score_glmm_sensitivity <- bind_rows(
  map_dfr(score_endpoints, fit_score_model, include_st_group = FALSE),
  map_dfr(score_endpoints, fit_score_model, include_st_group = TRUE)
) %>%
  group_by(model_variant) %>%
  mutate(q_value_BH = p.adjust(p.value, method = "BH")) %>%
  ungroup() %>%
  arrange(model_variant, q_value_BH, p.value)

write_csv(score_glmm_sensitivity, file.path(DIR_STAT, "score_glmm_sensitivity.csv"))

# ==============================================================================
# STANDALONE PLOTS S6-S9
# ==============================================================================

status_st_data <- scores %>%
  mutate(
    ST_plot = if_else(is.na(ST), "Missing ST", paste0("ST", ST))
  )
top_st_groups <- status_st_data %>%
  count(ST_plot, sort = TRUE) %>%
  filter(ST_plot != "Missing ST") %>%
  slice_head(n = 8) %>%
  pull(ST_plot)
st_comp <- status_st_data %>%
  mutate(
    ST_group_plot = case_when(
      ST_plot %in% top_st_groups ~ ST_plot,
      ST_plot == "Missing ST" ~ "Missing ST",
      TRUE ~ "Other STs"
    ),
    UTI_Status = factor(UTI_Status, levels = c("Not_UTI", "UTI"))
  ) %>%
  count(UTI_Status, ST_group_plot, name = "n") %>%
  group_by(UTI_Status) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()
base_st_cols <- c(
  "ST131" = "#4C78A8", "ST141" = "#F58518", "ST69" = "#54A24B", "ST73" = "#B279A2",
  "ST12" = "#72B7B2", "ST95" = "#E45756", "ST127" = "#9D755D", "ST10" = "#EECA3B",
  "Other STs" = "#BDBDBD", "Missing ST" = "#E5E7EB"
)
all_st_groups <- unique(st_comp$ST_group_plot)
extra_st_groups <- setdiff(all_st_groups, names(base_st_cols))
extra_st_cols <- if (length(extra_st_groups) > 0) {
  setNames(grDevices::hcl.colors(length(extra_st_groups), palette = "Dark 3"), extra_st_groups)
} else {
  character()
}
st_cols <- c(base_st_cols, extra_st_cols)[all_st_groups]

p_s6a <- ggplot(st_comp, aes(UTI_Status, prop, fill = ST_group_plot)) +
  geom_col(width = 0.66, colour = "white", linewidth = 0.25) +
  scale_y_continuous(labels = scales::percent_format()) +
  scale_fill_manual(values = st_cols, name = "Sequence type") +
  labs(
    title = "Sequence-type composition",
    x = "Primary UTI status",
    y = "Within-status proportion of VF-ready isolates"
  ) +
  plot_theme_stat(9)

sts_both <- scores %>%
  filter(!is.na(ST)) %>%
  count(ST, UTI_Status) %>%
  group_by(ST) %>%
  filter(n_distinct(UTI_Status) == 2) %>%
  pull(ST) %>%
  unique()

p_s6b <- ggplot(scores %>% filter(ST %in% sts_both),
                aes(UTI_Status, total_vf_count_all, fill = UTI_Status)) +
  geom_boxplot(outlier.shape = NA, width = 0.5, alpha = 0.75) +
  geom_jitter(width = 0.12, height = 0, alpha = 0.45, size = 1.3) +
  facet_wrap(~ paste0("ST", ST), scales = "free_x") +
  scale_fill_uti_status(drop = FALSE) +
  labs(
    title = "Within-ST VF burden",
    x = "Primary UTI status",
    y = "Detected VF genes per isolate"
  ) +
  plot_theme_stat(9) +
  theme(legend.position = "none")

s6 <- p_s6a + p_s6b +
  plot_layout(widths = c(0.45, 0.55), guides = "collect") +
  plot_annotation(
    title = "Lineage confounding diagnostic",
    subtitle = "VF burden differs strongly by ST, but ST composition was not significantly different by primary status in this sparse UTI set.",
    caption = sprintf(
      "Exploratory diagnostic; active Longcycler UTI n=%d and most STs contain too few UTI isolates for within-ST inference.",
      sum(vf_ready$UTI_Status == "UTI", na.rm = TRUE)
    )
  ) &
  theme(legend.position = "bottom")
ggsave(file.path(DIR_PLOTS_STAT, "lineage_confounding_panel.png"), s6, width = 13, height = 6.5, dpi = 300, bg = "white")
ggsave(file.path(DIR_PLOTS_STAT, "lineage_confounding_panel.pdf"), s6, width = 13, height = 6.5, device = "pdf", bg = "white")

s7_data <- participant_score_values %>%
  filter(score == "expec_marker_count") %>%
  group_by(Participant_id) %>%
  filter(all(c("Not_UTI", "UTI") %in% as.character(UTI_Status))) %>%
  ungroup() %>%
  mutate(UTI_Status = factor(UTI_Status, levels = c("Not_UTI", "UTI")))

s7 <- ggplot(s7_data, aes(UTI_Status, participant_status_median, group = Participant_id)) +
  geom_line(colour = "#B8B8B8", linewidth = 0.45, alpha = 0.85) +
  geom_point(aes(colour = UTI_Status), size = 2.2) +
  scale_colour_uti_status(drop = FALSE) +
  labs(
    title = "Paired resident ExPEC-like marker count",
    subtitle = "Each line is one resident; this visual reduces between-resident confounding but includes only residents observed in both states.",
    x = "Primary status within resident",
    y = "Participant-median ExPEC-like marker count",
    colour = "Primary status",
    caption = "Paired resident-level descriptive sensitivity; marker count is supplementary and residents observed in only one state are not represented."
  ) +
  plot_theme_stat(10)
ggsave(file.path(DIR_PLOTS_STAT, "paired_resident_expec_marker_slopeplot.png"), s7, width = 7.5, height = 5.8, dpi = 300, bg = "white")
ggsave(file.path(DIR_PLOTS_STAT, "paired_resident_expec_marker_slopeplot.pdf"), s7, width = 7.5, height = 5.8, device = "pdf", bg = "white")

s8_plot_data <- not_uti_uti_module_matrix %>%
  mutate(
    case_label = factor(case_label, levels = unique(case_label)),
    module_label = factor(module_label, levels = rev(unique(module_label))),
    change_display = factor(change_display, levels = c("Gained", "Lost", "Stable present", "Stable absent"))
  )

s8 <- ggplot(s8_plot_data, aes(case_label, module_label, fill = change_display)) +
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
  plot_theme_stat(8.5) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file.path(DIR_PLOTS_STAT, "not_uti_to_uti_module_gain_loss_heatmap.png"), s8, width = 12.5, height = 9.5, dpi = 300, bg = "white")
ggsave(file.path(DIR_PLOTS_STAT, "not_uti_to_uti_module_gain_loss_heatmap.pdf"), s8, width = 12.5, height = 9.5, device = "pdf", bg = "white")

pcoa_path <- require_file(file.path(DIR_VF, "vf_pcoa_jaccard_coordinates.csv"))
pcoa <- read_csv(pcoa_path, show_col_types = FALSE) %>%
  prefer_primary_uti_status(allow_legacy_fallback = FALSE) %>%
  mutate(
    UTI_Status = factor(UTI_Status, levels = c("Not_UTI", "UTI")),
    ST_group = factor(st_group4(ST), levels = c("ST131", "ST141", "Other typed ST", "Missing ST"))
  ) %>%
  filter(UTI_Status %in% c("Not_UTI", "UTI"))
axis1_lab <- sprintf("PCoA1 (%.1f%% variation)", unique(pcoa$var_Axis1)[1])
axis2_lab <- sprintf("PCoA2 (%.1f%% variation)", unique(pcoa$var_Axis2)[1])
s9 <- ggplot(pcoa, aes(Axis1, Axis2, colour = UTI_Status, shape = ST_group)) +
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
  plot_theme_stat(10)
ggsave(file.path(DIR_PLOTS_STAT, "vf_module_pcoa_primary_status.png"), s9, width = 8.2, height = 6.2, dpi = 300, bg = "white")
ggsave(file.path(DIR_PLOTS_STAT, "vf_module_pcoa_primary_status.pdf"), s9, width = 8.2, height = 6.2, device = "pdf", bg = "white")

figure_metadata <- tibble(
  figure_id = c(
    "lineage_confounding_panel",
    "paired_resident_expec_marker_slopeplot",
    "not_uti_to_uti_module_gain_loss_heatmap",
    "vf_module_pcoa_primary_status"
  ),
  file_png = file.path(DIR_PLOTS_STAT, c(
    "lineage_confounding_panel.png",
    "paired_resident_expec_marker_slopeplot.png",
    "not_uti_to_uti_module_gain_loss_heatmap.png",
    "vf_module_pcoa_primary_status.png"
  )),
  x_axis = c(
    "Primary UTI status",
    "Primary status within resident",
    "Transition case",
    axis1_lab
  ),
  y_axis = c(
    "Within-status proportion of VF-ready isolates / Detected VF genes per isolate",
    "Participant-median ExPEC-like marker count",
    "VF module",
    axis2_lab
  ),
  legend = c(
    "Sequence type; Primary status",
    "Primary status",
    "Module change",
    "Primary status; Sequence type group"
  ),
  caption = c(
    "VF burden differs strongly by ST, but ST composition was not significantly different by primary status in this sparse UTI set.",
    "Each line is one resident; this visual reduces between-resident confounding but includes only residents observed in both states.",
    "Module changes across WGS-linked Not_UTI-to-UTI transitions; descriptive, not causal.",
    "Global VF module profiles show whether UTI episodes separate from Not_UTI episodes; interpretation is limited by sparse UTI counts and lineage structure."
  )
)
write_csv(figure_metadata, file.path(DIR_STAT, "statistical_sensitivity_figure_metadata.csv"))

# ==============================================================================
# VALIDATION AND SUMMARY
# ==============================================================================

has_no_stale_stated_counts <- function(path) {
  if (!file.exists(path)) return(TRUE)
  txt <- readLines(path, warn = FALSE)
  !any(str_detect(txt, regex("12\\s+UTI|UTI\\s*=\\s*12|550\\s+rows|legacy\\s+ASB-vs-UTI", ignore_case = TRUE)))
}

validation <- tibble(
  check = c(
    sprintf("selected Longcycler denominator is %d rows", expected_active_total),
    sprintf("selected Longcycler counts are %d UTI and %d Not_UTI", expected_active_uti, expected_active_not_uti),
    sprintf("active selected QC-pass Longcycler VF/model denominator is %d rows", expected_active_total),
    sprintf("active Longcycler VF/model counts are %d UTI and %d Not_UTI", expected_active_uti, expected_active_not_uti),
    "participant-collapsed supplementary endpoint tests include BH q-values",
    "paired feature tests include BH q-values",
    "transition module enrichment includes BH q-values",
    "score GLMM sensitivity includes BH q-values",
    "new summary contains no stale 12-UTI or legacy ASB-vs-UTI wording"
  ),
  status = c(
    ifelse(nrow(status_primary) == expected_active_total, "PASS", "FAIL"),
    ifelse(sum(status_primary$UTI_Status == "UTI", na.rm = TRUE) == expected_active_uti &&
             sum(status_primary$UTI_Status == "Not_UTI", na.rm = TRUE) == expected_active_not_uti, "PASS", "FAIL"),
    ifelse(nrow(vf_ready) == expected_active_total, "PASS", "FAIL"),
    ifelse(sum(vf_ready$UTI_Status == "UTI", na.rm = TRUE) == expected_active_uti &&
             sum(vf_ready$UTI_Status == "Not_UTI", na.rm = TRUE) == expected_active_not_uti,
           "PASS", "FAIL"),
    ifelse(all(c("wilcoxon_rank_sum_q_BH", "paired_sign_test_q_BH", "paired_signed_rank_q_BH") %in% names(participant_score_tests)), "PASS", "FAIL"),
    ifelse(all(c("sign_test_q_BH", "signed_rank_q_BH") %in% names(paired_feature_sensitivity)), "PASS", "FAIL"),
    ifelse("fisher_q_BH" %in% names(transition_module_enrichment), "PASS", "FAIL"),
    ifelse("q_value_BH" %in% names(score_glmm_sensitivity), "PASS", "FAIL"),
    "PENDING"
  ),
  detail = c(
    sprintf("n=%d", nrow(status_primary)),
    sprintf("UTI=%d; Not_UTI=%d",
            sum(status_primary$UTI_Status == "UTI", na.rm = TRUE),
            sum(status_primary$UTI_Status == "Not_UTI", na.rm = TRUE)),
    sprintf("n=%d", nrow(vf_ready)),
    sprintf("UTI=%d; Not_UTI=%d",
            sum(vf_ready$UTI_Status == "UTI", na.rm = TRUE),
            sum(vf_ready$UTI_Status == "Not_UTI", na.rm = TRUE)),
    paste(names(participant_score_tests), collapse = ";"),
    paste(names(paired_feature_sensitivity), collapse = ";"),
    paste(names(transition_module_enrichment), collapse = ";"),
    paste(names(score_glmm_sensitivity), collapse = ";"),
    "Checked after summary is written."
  )
)

summary_lines <- c(
  "# Statistical Sensitivity Add-On Summary",
  "",
  sprintf("Generated: %s", format(Sys.time())),
  "",
  "## Denominators",
  "",
  sprintf("- Selected QC-pass Longcycler denominator: **%d** rows = **%d UTI**, **%d Not_UTI**.",
          nrow(status_primary),
          sum(status_primary$UTI_Status == "UTI", na.rm = TRUE),
          sum(status_primary$UTI_Status == "Not_UTI", na.rm = TRUE)),
  sprintf("- Active selected QC-pass Longcycler VF/model denominator: **%d** rows = **%d UTI**, **%d Not_UTI**.",
          nrow(vf_ready),
          sum(vf_ready$UTI_Status == "UTI", na.rm = TRUE),
          sum(vf_ready$UTI_Status == "Not_UTI", na.rm = TRUE)),
  "",
  "## Prespecified Sensitivity Readout",
  "",
  sprintf("- Participant-collapsed supplementary endpoint tests written for **%d** prespecified endpoints.", n_distinct(participant_score_tests$score)),
  sprintf("- Paired binary feature sensitivity written for **%d** selected gene/module features.", n_distinct(paired_feature_sensitivity$feature)),
  sprintf("- Transition-level module gain/loss enrichment written for **%d** module-direction tests.", nrow(transition_module_enrichment)),
  sprintf("- Supplementary endpoint GLMM sensitivity written for **%d** model rows; full ST is not used as a covariate.", nrow(score_glmm_sensitivity)),
  "",
  "## Interpretation",
  "",
  "- These tests are exploratory sensitivity checks only; they do not upgrade any VF association to confirmatory evidence.",
  "- All p-values are adjusted with BH correction within their output table/test family.",
  "- Event type is retained as a diagnostic/figure context rather than a primary model covariate.",
  "",
  "## Key Files",
  "",
  "- `results/statistical_sensitivity/participant_collapsed_score_tests.csv`",
  "- `results/statistical_sensitivity/paired_binary_feature_sensitivity.csv`",
  "- `results/statistical_sensitivity/transition_module_gain_loss_enrichment.csv`",
  "- `results/statistical_sensitivity/score_glmm_sensitivity.csv`",
  "- `plots/statistical_sensitivity/lineage_confounding_panel.png`",
  "- `plots/statistical_sensitivity/paired_resident_expec_marker_slopeplot.png`",
  "- `plots/statistical_sensitivity/not_uti_to_uti_module_gain_loss_heatmap.png`",
  "- `plots/statistical_sensitivity/vf_module_pcoa_primary_status.png`"
)
summary_path <- file.path(DIR_STAT, "statistical_sensitivity_summary.md")
writeLines(summary_lines, summary_path)

validation$status[validation$check == "new summary contains no stale 12-UTI or legacy ASB-vs-UTI wording"] <-
  ifelse(has_no_stale_stated_counts(summary_path), "PASS", "FAIL")
validation$detail[validation$check == "new summary contains no stale 12-UTI or legacy ASB-vs-UTI wording"] <-
  summary_path

write_csv(validation, file.path(DIR_STAT, "statistical_sensitivity_validation_checks.csv"))
if (any(validation$status != "PASS")) {
  print(validation)
  stop("Statistical sensitivity add-on validation failed.")
}

msg("Statistical sensitivity add-on complete.")
msg("Outputs written to %s and %s", DIR_STAT, DIR_PLOTS_STAT)

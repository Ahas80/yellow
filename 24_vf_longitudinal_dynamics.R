#!/usr/bin/env Rscript
# ==============================================================================
# 24_vf_longitudinal_dynamics.R
# ==============================================================================
#
# GOAL:
#   Assess within-individual VF profile stability across ALL consecutive
#   timepoint pairs.  Answers: "How stable are VF profiles within the same
#   patient over time?"
#
# WHY THIS SCRIPT EXISTS:
#   This is a longitudinal cohort — each participant contributes multiple
#   urine samples over months/years.  Understanding whether VF profiles
#   are stable within a host is critical:
#     - If VF profiles are highly stable,  it suggests the same bacterial
#       strain persists and VF content is a lineage trait, not dynamically
#       acquired.
#     - If VF profiles change frequently, it could mean VF gain/loss events
#       or strain replacement, which has clinical implications.
#
#   The existing script 16_within_host_evolution.R only examines VF changes
#   for phenotype-switch pairs. This script uses the primary UTI/Not_UTI status.
#   This is a much narrower question.  This script examines ALL consecutive
#   timepoint pairs regardless of status change, giving a comprehensive
#   picture of VF dynamics.
#
# METHODOLOGY:
#   For each consecutive pair of timepoints within a participant:
#     1. Extract the set of VF genes present at time t and time t+1
#     2. Compute Jaccard similarity = |intersection| / |union|
#        (1.0 = identical VF profile, 0.0 = completely different)
#     3. Identify specific genes gained and lost
#     4. Label the transition type (e.g., Not_UTI→Not_UTI, Not_UTI→UTI)
#
# INPUTS:
#   - results/vf/vf_analysis_ready.csv   (from 22_vf_build_analysis_dataset.R)
#   - results/vf/gene_map.csv            (for gene metadata)
#   - results/strain_compare/pairwise_metrics.csv
#     (same-strain/replacement evidence from 11_compare_strains.R)
#
# OUTPUTS (all in results/vf/):
#   - vf_longitudinal_transitions.csv      One row per consecutive pair
#                                          (full cohort)
#   - vf_transition_summary_by_type.csv    Aggregated by transition type
#   - vf_transitions_stratified.csv        All pairs, by depth (≥2/≥3/≥4)
#   - vf_transition_summary_stratified.csv Summaries, by depth
#   - vf_same_strain_vf_stability_summary.csv
#   - vf_replacement_vf_change_summary.csv
#   - vf_strain_context_by_transition_summary.csv
#   - vf_same_strain_by_ST_summary.csv
#   - vf_longitudinal_summary.txt          Human-readable summary
#
# PLOTS (in plots/vf/):
#   - vf_jaccard_by_transition.png         Jaccard boxplot by transition type
#
# TIMEPOINT ORDERING:
#   T0 < T1 < T2 < T3 < T4 < Uricult
#   Uricult (event-triggered) samples are placed last.  This is a convention
#   for ordering; the actual temporal spacing between timepoints varies.
#
# SUPERSEDES:
#   - Section 4 (within-individual dynamics) of compute_vf_abstract_stats.R
#   - Longitudinal sections of compute_vf_stratified_by_depth.R
# ==============================================================================

source("00_config.R")
source("R/plot_helpers.R")  # Primary UTI status colours
suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
})

msg("Starting 24_vf_longitudinal_dynamics.R")

# ==============================================================================
# VF visualisation shared helpers
# ==============================================================================

STATUS_LEVELS <- c("UTI", "Not_UTI", "Unknown")

status_for_plot <- function(x) {
  x <- as.character(x)
  x[is.na(x) | x == ""] <- "Unknown"
  x[!x %in% STATUS_LEVELS] <- "Unknown"
  factor(x, levels = STATUS_LEVELS)
}

status_count_text <- function(data) {
  if (!"Infection_Status" %in% names(data)) return("status denominator not applicable")
  counts <- data %>%
    mutate(.status = status_for_plot(Infection_Status)) %>%
    count(.status, name = "n") %>%
    filter(n > 0) %>%
    mutate(label = paste0(as.character(.status), " n=", n))
  paste(counts$label, collapse = "; ")
}

vf_caption <- function(input_file, data, analysis_unit,
                       p_value_note = "No inferential p-value is shown; longitudinal summaries are descriptive.",
                       extra_note = NULL) {
  n_uti <- if ("status_to" %in% names(data)) {
    sum(data$status_from == "UTI" | data$status_to == "UTI", na.rm = TRUE)
  } else if ("Infection_Status" %in% names(data)) {
    sum(data$Infection_Status == "UTI", na.rm = TRUE)
  } else {
    NA_integer_
  }
  paste(
    sprintf("Data: %s.", input_file),
    sprintf("Denominator: n=%d within-resident comparisons from %d participants.",
            nrow(data), n_distinct(data$Participant_id)),
    sprintf("Level of analysis: %s.", analysis_unit),
    "Jaccard similarity = 1 indicates identical binary VF presence/absence profiles.",
    "Comparisons are consecutive observed isolates within Participant_id; time intervals can vary.",
    p_value_note,
    sprintf("Comparisons involving UTI are limited (n=%d), so transition-specific patterns are underpowered.", n_uti),
    sprintf("SNP context: 0-%d SNPs = strong same strain; >%d SNPs = above same-strain SNP threshold; missing SNPs remain missing SNP evidence. ST is reported separately as lineage context.",
            strain_snp_threshold(), strain_snp_threshold()),
    extra_note %||% "",
    sep = " "
  ) %>% str_squish()
}

plot_theme_vf <- function(base_size = 11) {
  theme_bw(base_size = base_size) +
    theme(
      plot.caption = element_text(hjust = 0, size = base_size - 3, colour = "grey35"),
      plot.subtitle = element_text(colour = "grey25"),
      legend.position = "bottom"
    )
}

parse_date_safe <- function(x) {
  x <- as.character(x)
  out <- suppressWarnings(as.Date(x, format = "%d/%m/%Y"))
  missing <- is.na(out)
  if (any(missing)) out[missing] <- suppressWarnings(as.Date(x[missing], format = "%Y-%m-%d"))
  out
}

fallback_time_order <- function(tp) {
  tp <- as.character(tp)
  upper <- str_to_upper(tp)
  case_when(
    str_detect(upper, "^T\\d+$") ~ suppressWarnings(as.numeric(str_remove(upper, "^T"))),
    str_detect(upper, "^UTI-\\d+$") ~ 100 + suppressWarnings(as.numeric(str_remove(upper, "^UTI-"))),
    str_detect(tp, regex("uricult", ignore_case = TRUE)) ~ 99,
    TRUE ~ NA_real_
  )
}

normalise_st_label <- function(x) {
  x <- as.character(x)
  x <- str_trim(x)
  unknown <- c("", "-", "ST-", "NA", "N/A", "UNKNOWN", "UNK", "NT",
               "NON-TYPABLE", "NONTYPABLE", "NOT TYPED")
  x[str_to_upper(x) %in% unknown] <- NA_character_
  x
}

stop_if_stale <- function(target, upstream, target_label, upstream_label) {
  if (!file.exists(target) || !file.exists(upstream)) return(invisible(FALSE))
  target_mtime <- file.info(target)$mtime
  upstream_mtime <- file.info(upstream)$mtime
  if (!is.na(target_mtime) && !is.na(upstream_mtime) && target_mtime < upstream_mtime) {
    if (basename(target) == "vf_analysis_ready.csv" && basename(upstream) == "vf_pa_all.csv") {
      target_probe <- read_csv(target, show_col_types = FALSE, n_max = Inf)
      upstream_probe <- read_csv(upstream, show_col_types = FALSE, n_max = Inf)
      target_genes <- canonical_vf_gene_cols(names(target_probe), vf_pa_file = upstream)
      upstream_genes <- canonical_vf_gene_cols(names(upstream_probe), vf_pa_file = upstream)
      if (nrow(target_probe) == nrow(upstream_probe) && setequal(target_genes, upstream_genes)) {
        msg("WARNING: %s timestamp is older than %s, but row count and gene set match; continuing.",
            target_label, upstream_label)
        return(invisible(FALSE))
      }
    }
    stop(sprintf(
      "%s is older than %s. Re-run 22_vf_build_analysis_dataset.R before script 24 so longitudinal plots use the current VF matrix.\n  %s: %s\n  %s: %s",
      target_label, upstream_label, target_label, format(target_mtime),
      upstream_label, format(upstream_mtime)
    ))
  }
}

# ==============================================================================
# 1. LOAD DATA
# ==============================================================================

ready_file <- FILE_VF_READY
if (!file.exists(ready_file)) stop("Missing ", ready_file, ". Run 22_vf_build_analysis_dataset.R first.")
stop_if_stale(ready_file, FILE_VF_PA, "vf_analysis_ready.csv", "vf_pa_all.csv")
vf_ready <- read_csv(ready_file, show_col_types = FALSE) %>%
  prefer_primary_uti_status() %>%
  apply_manual_sample_curation(context = "24_vf_ready") %>%
  filter_primary_genomics() %>%
  mutate(Participant_id = as.character(Participant_id),
         ST = if ("ST" %in% names(.)) normalise_st_label(ST) else NA_character_)

gene_map <- read_csv(file.path(DIR_VF, "gene_map.csv"), show_col_types = FALSE) %>%
  mutate(Gene = as.character(Gene),
         Category = coalesce(as.character(Category), "Unassigned"))

# Individual VF gene columns are defined by the canonical P/A matrix from
# 02_gene_presence_analysis.R. Do not infer genes by excluding metadata here:
# vf_analysis_ready.csv also contains clinical provenance and score columns.
cat_cols  <- grep("^cat_", names(vf_ready), value = TRUE)
gene_cols <- canonical_vf_gene_cols(names(vf_ready))

msg("Loaded: %d rows, %d gene columns", nrow(vf_ready), length(gene_cols))
msg("VF-ready denominator by status: %s", status_count_text(vf_ready))

pairwise_file <- file.path(DIR_STRAIN, "pairwise_metrics.csv")
if (!file.exists(pairwise_file)) {
  stop("Missing ", pairwise_file, ". Run 11_compare_strains.R before script 24 so VF longitudinal dynamics are interpreted same-strain-first.")
}
pairwise <- read_csv(pairwise_file, show_col_types = FALSE) %>%
  prepare_pairwise_for_strain_context()
msg("Loaded pairwise strain metrics: %d rows; same-strain SNP threshold <=%d",
    nrow(pairwise), strain_snp_threshold())

# ==============================================================================
# 2. TRANSITION-BUILDING FUNCTION
# ==============================================================================
# This function takes a filtered dataset and builds all consecutive-pair
# transitions for participants with ≥2 timepoints.

build_transitions <- function(data, gene_cols, cohort_label = "all") {

  # Filter to episodes with primary UTI status, then sort by true collection date
  # where available. Display-only fallback order is used only when dates are
  # unavailable and is explicitly flagged in the output.
  long_data <- data %>%
    filter(!is.na(Infection_Status)) %>%
    mutate(
      Collection_Date_parsed = if ("Collection_Date" %in% names(.)) parse_date_safe(Collection_Date) else as.Date(NA),
      fallback_order = fallback_time_order(tp_lab)
    ) %>%
    group_by(Participant_id) %>%
    mutate(
      first_collection_date = if (all(is.na(Collection_Date_parsed))) as.Date(NA) else min(Collection_Date_parsed, na.rm = TRUE),
      date_order = as.numeric(Collection_Date_parsed - first_collection_date),
      time_order = coalesce(date_order, fallback_order),
      time_order_source = case_when(
        !is.na(date_order) ~ "Collection_Date",
        !is.na(fallback_order) & str_detect(tp_lab, regex("uricult|^UTI-", ignore_case = TRUE)) ~
          "Fallback_event_label_order",
        !is.na(fallback_order) ~ "Fallback_routine_timepoint_order",
        TRUE ~ "Unavailable"
      )
    ) %>%
    ungroup() %>%
    filter(!is.na(time_order)) %>%
    arrange(Participant_id, time_order, tp_lab)

  # Keep only participants with ≥2 timepoints (need at least a pair)
  multi_tp <- long_data %>%
    group_by(Participant_id) %>%
    filter(n() >= 2) %>%
    ungroup()

  if (nrow(multi_tp) == 0) {
    msg("  No participants with ≥2 timepoints in cohort '%s'.  Skipping.", cohort_label)
    return(list(transitions = tibble(), summary = tibble()))
  }

  # Build consecutive pairs for each participant.
  # For a participant with timepoints T0, T1, T2 we get pairs:
  #   (T0→T1) and (T1→T2)
  transitions <- list()
  for (pid in unique(multi_tp$Participant_id)) {
    pid_data <- multi_tp %>%
      filter(Participant_id == pid) %>%
      arrange(time_order, tp_lab)

    if (nrow(pid_data) < 2) next

    for (i in 1:(nrow(pid_data) - 1)) {
      rf <- pid_data[i, ]      # "from" row (earlier timepoint)
      rt <- pid_data[i + 1, ]  # "to" row   (later timepoint)

      # Extract the sets of VF genes present at each timepoint
      vf_from <- gene_cols[which(rf[, gene_cols] > 0)]
      vf_to   <- gene_cols[which(rt[, gene_cols] > 0)]

      # Compute set operations:
      #   gained = genes present at t+1 but NOT at t (potentially acquired)
      #   lost   = genes present at t but NOT at t+1 (potentially deleted/replaced)
      #   stable = genes present at BOTH timepoints
      gained  <- setdiff(vf_to, vf_from)
      lost    <- setdiff(vf_from, vf_to)
      stable  <- intersect(vf_from, vf_to)
      union_s <- union(vf_from, vf_to)

      # Jaccard similarity = |A ∩ B| / |A ∪ B|
      # Ranges from 0 (no overlap) to 1 (identical sets).
      # This is the primary measure of VF profile stability.
      jac <- if (length(union_s) == 0) NA_real_ else length(stable) / length(union_s)

      st_from <- if ("ST" %in% names(rf)) normalise_st_label(rf$ST) else NA_character_
      st_to <- if ("ST" %in% names(rt)) normalise_st_label(rt$ST) else NA_character_
      strain_ctx <- lookup_strain_context(
        pairwise = pairwise,
        pid = pid,
        tp_from = as.character(rf$tp_lab),
        tp_to = as.character(rt$tp_lab),
        ST_from = st_from,
        ST_to = st_to,
        has_vf_pair = TRUE
      )
      date_from <- rf$Collection_Date_parsed
      date_to <- rt$Collection_Date_parsed
      days_between <- if (!is.na(date_from) && !is.na(date_to)) {
        as.numeric(date_to - date_from)
      } else {
        NA_real_
      }

      transitions[[length(transitions) + 1]] <- tibble(
        Participant_id  = pid,
        tp_from         = as.character(rf$tp_lab),
        tp_to           = as.character(rt$tp_lab),
        status_from     = rf$Infection_Status,
        status_to       = rt$Infection_Status,
        transition_type = paste0(rf$Infection_Status, "\u2192", rt$Infection_Status),
        event_type_from = if ("Event_type" %in% names(rf)) as.character(rf$Event_type) else NA_character_,
        event_type_to   = if ("Event_type" %in% names(rt)) as.character(rt$Event_type) else NA_character_,
        collection_date_from = as.character(date_from),
        collection_date_to   = as.character(date_to),
        days_between_samples = days_between,
        time_order_from = rf$time_order,
        time_order_to = rt$time_order,
        time_order_source = paste(rf$time_order_source, rt$time_order_source, sep = " -> "),
        comparison_scope = "consecutive observed within-resident pair",
        ST_from = strain_ctx$ST_from,
        ST_to = strain_ctx$ST_to,
        same_ST = strain_ctx$same_ST,
        SNPs = strain_ctx$SNPs,
        AvgIdentity = strain_ctx$AvgIdentity,
        Pairwise_Classification = strain_ctx$Pairwise_Classification,
        Pairwise_RuleUsed = strain_ctx$Pairwise_RuleUsed,
        snp_strain_context = as.character(strain_ctx$snp_strain_context),
        st_lineage_context = as.character(strain_ctx$st_lineage_context),
        pair_interpretation = as.character(strain_ctx$pair_interpretation),
        same_strain_evidence = strain_ctx$same_strain_evidence,
        strain_context_level = as.character(strain_ctx$strain_context_level),
        replacement_flag = strain_ctx$replacement_flag,
        strain_context_note = strain_ctx$strain_context_note,
        same_strain_snp_threshold = strain_ctx$same_strain_snp_threshold,
        vf_count_from   = rf$vf_count_total,
        vf_count_to     = rt$vf_count_total,
        n_gained        = length(gained),
        n_lost          = length(lost),
        n_stable        = length(stable),
        n_from          = length(vf_from),
        n_to            = length(vf_to),
        jaccard_similarity = round(jac, 3),
        any_vf_change   = length(gained) > 0 | length(lost) > 0,
        genes_gained    = paste(sort(gained), collapse = ";"),
        genes_lost      = paste(sort(lost),   collapse = ";")
      )
    }
  }

  trans_df <- bind_rows(transitions) %>%
    mutate(cohort = cohort_label)

  # Aggregate by transition type (e.g., Not_UTI->Not_UTI, Not_UTI->UTI, etc.)
  # This summary is useful for comparing VF stability across different
  # clinical trajectory types.
  trans_summary <- trans_df %>%
    group_by(transition_type) %>%
    summarise(
      n_transitions    = n(),
      n_participants   = n_distinct(Participant_id),
      n_with_days_between = sum(!is.na(days_between_samples)),
      median_days_between = ifelse(all(is.na(days_between_samples)), NA_real_, median(days_between_samples, na.rm = TRUE)),
      n_snp_strong_same_strain = sum(snp_strain_context == "Strong same strain", na.rm = TRUE),
      n_snp_above_same_strain_threshold = sum(snp_strain_context == "Above same-strain SNP threshold", na.rm = TRUE),
      n_snp_missing = sum(snp_strain_context == "Missing SNP evidence", na.rm = TRUE),
      n_pair_replacement_likely = sum(pair_interpretation == "Replacement likely", na.rm = TRUE),
      n_pair_same_lineage_not_same_strain = sum(pair_interpretation == "Same lineage, not same strain by SNP", na.rm = TRUE),
      n_pair_st_consistent_snp_missing = sum(pair_interpretation == "ST-consistent, SNP missing", na.rm = TRUE),
      n_pair_missing_strain_metrics = sum(pair_interpretation == "Missing strain metrics", na.rm = TRUE),
      n_same_ST = sum(same_ST %in% TRUE, na.rm = TRUE),
      n_different_ST = sum(same_ST %in% FALSE, na.rm = TRUE),
      median_gained    = median(n_gained),
      mean_gained      = round(mean(n_gained), 1),
      median_lost      = median(n_lost),
      mean_lost        = round(mean(n_lost), 1),
      median_stable    = median(n_stable),
      median_jaccard   = round(median(jaccard_similarity, na.rm = TRUE), 3),
      mean_jaccard     = round(mean(jaccard_similarity, na.rm = TRUE), 3),
      pct_no_change    = round(mean(!any_vf_change) * 100, 1),
      .groups          = "drop"
    ) %>%
    mutate(cohort = cohort_label) %>%
    arrange(desc(n_transitions))

  list(transitions = trans_df, summary = trans_summary)
}

# ==============================================================================
# 3. RUN FOR FULL COHORT + STRATIFIED BY TIMEPOINT DEPTH
# ==============================================================================
# Same stratification logic as 23_: analyses are repeated for participants
# with ≥2, ≥3, and ≥4 timepoints to assess sensitivity to follow-up depth.

cohorts <- list(
  "all"    = vf_ready,
  ">=2 tp" = vf_ready %>% filter(!is.na(n_timepoints) & n_timepoints >= 2),
  ">=3 tp" = vf_ready %>% filter(!is.na(n_timepoints) & n_timepoints >= 3),
  ">=4 tp" = vf_ready %>% filter(!is.na(n_timepoints) & n_timepoints >= 4)
)

results <- lapply(names(cohorts), function(label) {
  msg("  Building transitions for cohort: %s (%d rows)", label, nrow(cohorts[[label]]))
  build_transitions(cohorts[[label]], gene_cols, label)
})
names(results) <- names(cohorts)

# ==============================================================================
# 4. WRITE OUTPUTS
# ==============================================================================

msg("Writing outputs...")

# Full-cohort outputs
write_csv(results[["all"]]$transitions, file.path(DIR_VF, "vf_longitudinal_transitions.csv"))
write_csv(results[["all"]]$summary,     file.path(DIR_VF, "vf_transition_summary_by_type.csv"))

# Stratified outputs (all depth levels stacked in single CSVs, with
# a "cohort" column to distinguish them)
all_trans <- bind_rows(lapply(results, `[[`, "transitions"))
all_summ  <- bind_rows(lapply(results, `[[`, "summary"))
full_trans <- results[["all"]]$transitions

write_csv(all_trans, file.path(DIR_VF, "vf_transitions_stratified.csv"))
write_csv(all_summ,  file.path(DIR_VF, "vf_transition_summary_stratified.csv"))

vf_stability_summary <- function(data, group_cols) {
  summary_cols <- c(
    "n_transitions", "n_participants", "median_vf_jaccard", "q25_vf_jaccard",
    "q75_vf_jaccard", "pct_no_vf_gene_change", "median_genes_gained",
    "median_genes_lost", "median_snp", "n_same_ST", "n_different_ST",
    "same_strain_snp_threshold", "strong_same_strain_snp_range",
    "above_same_strain_snp_rule", "st_lineage_rule", "replacement_likely_rule",
    "missing_strain_metrics_rule"
  )
  if (nrow(data) == 0) {
    empty_groups <- setNames(rep(list(character()), length(group_cols)), group_cols)
    numeric_summary <- setNames(rep(list(numeric()), 12), summary_cols[1:12])
    character_summary <- setNames(rep(list(character()), length(summary_cols) - 12), summary_cols[-seq_len(12)])
    return(as_tibble(c(empty_groups, numeric_summary, character_summary)))
  }
  data %>%
    group_by(across(all_of(group_cols))) %>%
    summarise(
      n_transitions = n(),
      n_participants = n_distinct(Participant_id),
      median_vf_jaccard = median(jaccard_similarity, na.rm = TRUE),
      q25_vf_jaccard = ifelse(all(is.na(jaccard_similarity)), NA_real_, quantile(jaccard_similarity, 0.25, na.rm = TRUE)),
      q75_vf_jaccard = ifelse(all(is.na(jaccard_similarity)), NA_real_, quantile(jaccard_similarity, 0.75, na.rm = TRUE)),
      pct_no_vf_gene_change = round(mean(!any_vf_change, na.rm = TRUE) * 100, 1),
      median_genes_gained = median(n_gained, na.rm = TRUE),
      median_genes_lost = median(n_lost, na.rm = TRUE),
      median_snp = ifelse(all(is.na(SNPs)), NA_real_, median(SNPs, na.rm = TRUE)),
      n_same_ST = sum(same_ST %in% TRUE, na.rm = TRUE),
      n_different_ST = sum(same_ST %in% FALSE, na.rm = TRUE),
      same_strain_snp_threshold = strain_snp_threshold(),
      strong_same_strain_snp_range = sprintf("0-%d SNPs", strain_snp_threshold()),
      above_same_strain_snp_rule = sprintf(">%d SNPs = above same-strain SNP threshold", strain_snp_threshold()),
      st_lineage_rule = "ST is reported separately as Same ST, Different ST, or Missing ST evidence",
      replacement_likely_rule = "Different ST or pairwise Different when SNPs do not support same strain",
      missing_strain_metrics_rule = "Missing SNP evidence with no usable ST/pairwise context or missing WGS/VF endpoint",
      .groups = "drop"
    ) %>%
    arrange(desc(n_transitions))
}

same_strain_trans <- full_trans %>%
  filter(snp_strain_context == "Strong same strain")
replacement_trans <- full_trans %>%
  filter(pair_interpretation == "Replacement likely")

write_csv(
  vf_stability_summary(same_strain_trans, c("transition_type")),
  file.path(DIR_VF, "vf_same_strain_vf_stability_summary.csv")
)
write_csv(
  vf_stability_summary(replacement_trans, c("transition_type")),
  file.path(DIR_VF, "vf_replacement_vf_change_summary.csv")
)
write_csv(
  vf_stability_summary(full_trans, c("transition_type", "snp_strain_context", "pair_interpretation", "st_lineage_context")),
  file.path(DIR_VF, "vf_strain_context_by_transition_summary.csv")
)
write_csv(
  same_strain_trans %>%
    mutate(ST_context = ifelse(!is.na(ST_from), paste0("ST", ST_from), "Missing/non-typable ST")) %>%
    vf_stability_summary(c("ST_context", "transition_type")),
  file.path(DIR_VF, "vf_same_strain_by_ST_summary.csv")
)

# ==============================================================================
# 5. PLOT: JACCARD SIMILARITY BY TRANSITION TYPE
# ==============================================================================
# This figure shows how stable VF profiles are within each type of clinical
# transition (e.g., Not_UTI→Not_UTI vs Not_UTI→UTI). High Jaccard values (close to 1)
# indicate the VF profile barely changes between consecutive observed isolates.

ensure_dir(DIR_PLOTS_VF)

if (nrow(full_trans) > 0) {
  full_trans <- full_trans %>%
    mutate(
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
    )
  type_counts <- full_trans %>% count(transition_type, sort = TRUE)
  common_types <- type_counts %>% filter(n >= 3) %>% pull(transition_type)

  if (nrow(same_strain_trans) > 0) {
    same_plot <- same_strain_trans %>%
      mutate(plot_type = ifelse(transition_type %in% common_types,
                                transition_type, "Other"),
             plot_type = factor(plot_type, levels = c(common_types, "Other")))
    p_same_jac <- ggplot(same_plot,
                         aes(x = reorder(plot_type, -jaccard_similarity, FUN = median),
                             y = jaccard_similarity)) +
      geom_boxplot(outlier.shape = NA, width = 0.5, fill = "#BFE3D0") +
      geom_jitter(width = 0.15, alpha = 0.55, size = 1.5, colour = "#2F855A") +
      geom_hline(yintercept = 1.0, linetype = "dashed", colour = "grey50") +
      labs(
        title = "VF stability within strong same-strain repeated isolates",
        subtitle = sprintf("Same strain defined first by SNP distance <=%d; ST is secondary context",
                           strain_snp_threshold()),
        x = "Clinical transition type",
        y = "Jaccard similarity (1 = identical VF profile)",
        caption = vf_caption(
          ready_file, same_plot, "same-strain consecutive within-resident isolate-pair comparison",
          extra_note = sprintf("Only pairs with strong same-strain evidence under the <=%d SNP threshold are shown.",
                               strain_snp_threshold())
        )
      ) +
      plot_theme_vf(base_size = 11) +
      theme(axis.text.x = element_text(angle = 30, hjust = 1))

    ggsave(file.path(DIR_PLOTS_VF, "vf_same_strain_jaccard_by_transition.png"),
           p_same_jac, width = 8.5, height = 5.8, dpi = 300)

    same_gain_loss <- same_strain_trans %>%
      select(Participant_id, transition_type, n_gained, n_lost) %>%
      pivot_longer(cols = c(n_gained, n_lost), names_to = "change_type", values_to = "n_genes") %>%
      mutate(change_type = recode(change_type,
                                  "n_gained" = "VF genes gained",
                                  "n_lost" = "VF genes lost"))
    p_same_gain_loss <- ggplot(same_gain_loss,
                               aes(x = change_type, y = n_genes, fill = change_type)) +
      geom_boxplot(outlier.shape = NA, alpha = 0.65, width = 0.55) +
      geom_jitter(width = 0.14, alpha = 0.35, size = 1.1) +
      facet_wrap(~transition_type) +
      scale_fill_manual(values = c("VF genes gained" = "#D55E00",
                                   "VF genes lost" = "#0072B2")) +
      labs(
        title = "VF gene gains and losses within strong same-strain pairs",
        subtitle = sprintf("Restricted to repeated isolates with SNP distance <=%d where SNPs are available",
                           strain_snp_threshold()),
        x = NULL,
        y = "Number of VF genes",
        caption = vf_caption(
          ready_file, same_strain_trans, "same-strain consecutive within-resident gene gain/loss summary",
          extra_note = "This same-strain-first view separates within-strain VF stability from strain replacement."
        )
      ) +
      plot_theme_vf(base_size = 10) +
      theme(axis.text.x = element_text(angle = 20, hjust = 1),
            legend.position = "none")

    ggsave(file.path(DIR_PLOTS_VF, "vf_same_strain_gene_gain_loss.png"),
           p_same_gain_loss, width = 10, height = 7, dpi = 300)

    same_st_plot <- same_strain_trans %>%
      mutate(ST_context = ifelse(!is.na(ST_from), paste0("ST", ST_from), "Missing/non-typable ST")) %>%
      group_by(ST_context) %>%
      mutate(n_ST_context = n()) %>%
      ungroup() %>%
      filter(n_ST_context >= 2)
    if (nrow(same_st_plot) > 0) {
      p_same_st <- ggplot(same_st_plot,
                          aes(x = reorder(ST_context, -jaccard_similarity, FUN = median),
                              y = jaccard_similarity)) +
        geom_boxplot(outlier.shape = NA, fill = "#E6F2EF", width = 0.55) +
        geom_jitter(width = 0.12, alpha = 0.55, size = 1.3, colour = "#2C7A7B") +
        labs(
          title = "Same-strain VF stability by sequence type",
          subtitle = "ST is shown after SNP-based same-strain filtering",
          x = "Sequence type",
          y = "Jaccard similarity",
          caption = vf_caption(
            ready_file, same_st_plot, "same-strain pair comparison stratified by ST",
            extra_note = "Only ST groups with at least two strong same-strain pairs are shown."
          )
        ) +
        plot_theme_vf(base_size = 10) +
        theme(axis.text.x = element_text(angle = 35, hjust = 1))

      ggsave(file.path(DIR_PLOTS_VF, "vf_same_strain_by_ST.png"),
             p_same_st, width = 8, height = 5.4, dpi = 300)
    } else {
      msg("Skipping same-strain by-ST plot: fewer than two same-strain pairs per ST.")
    }
  } else {
    msg("Skipping same-strain VF plots: no strong same-strain transitions under <=%d SNP threshold.", strain_snp_threshold())
  }

  p_ctx_jac <- ggplot(full_trans,
                      aes(x = pair_interpretation, y = jaccard_similarity,
                          fill = pair_interpretation)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.7, width = 0.58) +
    geom_jitter(width = 0.15, alpha = 0.4, size = 1.2) +
    scale_fill_manual(values = c("Strong same strain" = "#2F855A",
                                 "Conflict: SNP same-strain but ST differs" = "#6B46C1",
                                 "Same lineage, not same strain by SNP" = "#B7791F",
                                 "ST-consistent, SNP missing" = "#6B7280",
                                 "Above same-strain SNP threshold" = "#D69E2E",
                                 "Missing SNP evidence" = "grey70",
                                 "Replacement likely" = "#C05621",
                                 "Missing strain metrics" = "grey65")) +
    labs(
      title = "VF similarity by SNP-defined strain context and secondary lineage interpretation",
      subtitle = sprintf("SNP context uses <=%d for strong same strain; ST is shown separately as lineage context",
                         strain_snp_threshold()),
      x = "Pair interpretation",
      y = "Jaccard similarity",
      caption = vf_caption(
        ready_file, full_trans, "consecutive pair comparison by strain-context tier",
        extra_note = "SNP-defined same-strain calls drive the VF stability analysis; ST agreement is secondary lineage context and does not prove same strain."
      )
    ) +
    plot_theme_vf(base_size = 10) +
    theme(axis.text.x = element_text(angle = 25, hjust = 1),
          legend.position = "none")

  ggsave(file.path(DIR_PLOTS_VF, "vf_jaccard_by_strain_context.png"),
         p_ctx_jac, width = 8.5, height = 5.6, dpi = 300)

  change_context <- full_trans %>%
    group_by(pair_interpretation) %>%
    summarise(
      n_transitions = n(),
      pct_any_vf_change = mean(any_vf_change, na.rm = TRUE),
      median_total_gene_changes = median(n_gained + n_lost, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    filter(pair_interpretation %in% c("Strong same strain", "Replacement likely"))
  if (nrow(change_context) > 0) {
    p_change_ctx <- ggplot(change_context,
                           aes(x = pair_interpretation, y = pct_any_vf_change,
                               fill = pair_interpretation)) +
      geom_col(width = 0.62) +
      geom_text(aes(label = sprintf("%d pairs", n_transitions)), vjust = -0.35, size = 3.4) +
      scale_y_continuous(labels = scales::percent, expand = expansion(mult = c(0, 0.16))) +
      scale_fill_manual(values = c("Strong same strain" = "#2F855A",
                                   "Replacement likely" = "#C05621")) +
      labs(
        title = "VF change is separated into same-strain and replacement contexts",
        subtitle = "Proportion of repeated-isolate pairs with any VF gene gain or loss",
        x = NULL,
        y = "Pairs with any VF gene change",
        caption = vf_caption(
          ready_file, full_trans, "same-strain versus replacement VF-change diagnostic",
          extra_note = "This guards against interpreting replacement-driven profile differences as within-strain VF evolution."
        )
      ) +
      plot_theme_vf(base_size = 10) +
      theme(legend.position = "none")

    ggsave(file.path(DIR_PLOTS_VF, "vf_replacement_vs_same_strain_vf_change.png"),
           p_change_ctx, width = 7.5, height = 5.3, dpi = 300)
  }

  # Group rare transition types (those with < 3 occurrences) into "Other"
  # to keep the plot readable
  plot_trans <- full_trans %>%
    mutate(plot_type = ifelse(transition_type %in% common_types,
                              transition_type, "Other"),
           plot_type = factor(plot_type, levels = c(common_types, "Other")))

  p_jac <- ggplot(plot_trans, aes(x = reorder(plot_type, -jaccard_similarity, FUN = median),
                                   y = jaccard_similarity)) +
    geom_boxplot(outlier.shape = NA, width = 0.5, fill = "grey90") +
    geom_jitter(width = 0.15, alpha = 0.5, size = 1.5, colour = "#0072B2") +
    geom_hline(yintercept = 1.0, linetype = "dashed", colour = "grey50") +
    labs(
      title = "Within-resident virulence factor similarity across repeated E. coli isolates",
      subtitle = "Consecutive observed within-resident pairs; Jaccard similarity from binary VF gene presence/absence",
      x = "Clinical transition type",
      y = "Jaccard similarity (1 = identical VF profile)",
      caption = vf_caption(
        ready_file, full_trans, "consecutive within-resident isolate-pair comparison",
        extra_note = "Rare transition types with fewer than three comparisons are grouped as Other for readability."
      )
    ) +
    plot_theme_vf(base_size = 11) +
    theme(axis.text.x = element_text(angle = 30, hjust = 1))

  ggsave(file.path(DIR_PLOTS_VF, "vf_jaccard_by_transition.png"), p_jac,
         width = 8.5, height = 5.8, dpi = 300)

  p_dist <- ggplot(full_trans, aes(x = jaccard_similarity)) +
    geom_histogram(binwidth = 0.02, boundary = 0, fill = "#0072B2",
                   colour = "white", alpha = 0.85) +
    geom_vline(xintercept = 1, linetype = "dashed", colour = "grey45") +
    labs(
      title = "Distribution of within-resident VF profile similarity",
      subtitle = sprintf("%d consecutive comparisons from %d residents; median Jaccard %.3f",
                         nrow(full_trans), n_distinct(full_trans$Participant_id),
                         median(full_trans$jaccard_similarity, na.rm = TRUE)),
      x = "Jaccard similarity",
      y = "Consecutive isolate-pair comparisons",
      caption = vf_caption(
        ready_file, full_trans, "consecutive within-resident isolate-pair comparison",
        extra_note = "This plot summarizes all transition types together and does not imply causal VF stability."
      )
    ) +
    plot_theme_vf(base_size = 11)

  ggsave(file.path(DIR_PLOTS_VF, "vf_within_host_jaccard_distribution.png"),
         p_dist, width = 7.5, height = 5.4, dpi = 300)

  days_plot <- full_trans %>% filter(!is.na(days_between_samples))
  if (nrow(days_plot) > 0) {
    p_days <- ggplot(days_plot, aes(x = days_between_samples, y = jaccard_similarity)) +
      geom_point(aes(colour = transition_type), alpha = 0.65, size = 1.8) +
      geom_smooth(method = "loess", se = FALSE, colour = "grey30", linewidth = 0.6) +
      labs(
        title = "VF similarity versus time between repeated isolates",
        subtitle = "Days are calculated from Collection_Date where both dates are available",
        x = "Days between samples",
        y = "Jaccard similarity",
        colour = "Transition type",
        caption = vf_caption(
          ready_file, days_plot, "dated consecutive within-resident isolate-pair comparison",
          extra_note = "Pairs lacking usable Collection_Date values are excluded from this days-based diagnostic only."
        )
      ) +
      plot_theme_vf(base_size = 11)

    ggsave(file.path(DIR_PLOTS_VF, "vf_jaccard_by_days_between_samples.png"),
           p_days, width = 8, height = 5.6, dpi = 300)
  } else {
    msg("Skipping Jaccard-by-days plot: no pairs have two parseable Collection_Date values.")
  }

  st_plot <- full_trans %>%
    filter(st_lineage_context != "Missing ST evidence")
  if (nrow(st_plot) > 0 && n_distinct(st_plot$st_lineage_context) >= 1) {
    p_st <- ggplot(st_plot, aes(x = st_lineage_context, y = jaccard_similarity,
                                fill = st_lineage_context)) +
      geom_boxplot(outlier.shape = NA, alpha = 0.6, width = 0.55) +
      geom_jitter(width = 0.12, alpha = 0.45, size = 1.3) +
      scale_fill_manual(values = c("Same ST" = "#009E73", "Different ST" = "#D55E00")) +
      labs(
        title = "Secondary lineage diagnostic: VF similarity by ST context",
        subtitle = "ST is shown after SNP-defined same-strain interpretation and does not prove same strain",
        x = NULL,
        y = "Jaccard similarity",
        caption = vf_caption(
          ready_file, st_plot, "consecutive within-resident isolate-pair comparison stratified by ST",
          extra_note = "Use this only as lineage/confounding context after reviewing SNP-defined same-strain calls."
        )
      ) +
      plot_theme_vf(base_size = 11) +
      theme(legend.position = "none")

    ggsave(file.path(DIR_PLOTS_VF, "vf_jaccard_same_vs_different_st.png"),
           p_st, width = 6.5, height = 5.2, dpi = 300)
  } else {
    msg("Skipping secondary ST lineage-context Jaccard plot: insufficient ST-labelled pairs.")
  }

  gain_loss <- full_trans %>%
    select(Participant_id, transition_type, n_gained, n_lost) %>%
    pivot_longer(cols = c(n_gained, n_lost), names_to = "change_type", values_to = "n_genes") %>%
    mutate(change_type = recode(change_type,
                                "n_gained" = "VF genes gained",
                                "n_lost" = "VF genes lost"))

  p_gain_loss <- ggplot(gain_loss,
                        aes(x = change_type, y = n_genes, fill = change_type)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.65, width = 0.55) +
    geom_jitter(width = 0.14, alpha = 0.35, size = 1.1) +
    facet_wrap(~transition_type) +
    scale_fill_manual(values = c("VF genes gained" = "#D55E00",
                                 "VF genes lost" = "#0072B2")) +
    labs(
      title = "VF gene gains and losses between repeated E. coli isolates",
      subtitle = "Consecutive observed within-resident pairs; transition-specific panels are descriptive",
      x = NULL,
      y = "Number of VF genes",
      caption = vf_caption(
        ready_file, full_trans, "consecutive within-resident gene gain/loss summary",
        extra_note = "Gene gains/losses can reflect strain replacement, assembly/VF calling differences, or true gene-content change."
      )
    ) +
    plot_theme_vf(base_size = 10) +
    theme(axis.text.x = element_text(angle = 20, hjust = 1),
          legend.position = "none")

  ggsave(file.path(DIR_PLOTS_VF, "vf_gene_gain_loss_consecutive_pairs.png"),
         p_gain_loss, width = 10, height = 7, dpi = 300)
}

# ==============================================================================
# 6. HUMAN-READABLE SUMMARY
# ==============================================================================
# Plain-text summary of longitudinal findings, suitable for reference
# when drafting the abstract or results section.

summary_lines <- character()
sl <- function(...) summary_lines <<- c(summary_lines, sprintf(...))

sl("=== VF LONGITUDINAL DYNAMICS SUMMARY ===")
sl("Generated: %s", format(Sys.time()))
sl("")

if (nrow(full_trans) > 0) {
  sl("Total transitions (full cohort): %d from %d participants",
     nrow(full_trans), n_distinct(full_trans$Participant_id))
  sl("Median Jaccard: %.3f", median(full_trans$jaccard_similarity, na.rm = TRUE))
  sl("Mean Jaccard: %.3f",   mean(full_trans$jaccard_similarity, na.rm = TRUE))
  sl("%% with zero VF change: %.1f%%", mean(!full_trans$any_vf_change) * 100)
  sl("Median genes gained: %.0f", median(full_trans$n_gained))
  sl("Median genes lost: %.0f",   median(full_trans$n_lost))
  sl("")

  # Special focus on Not_UTI→UTI transitions: are specific genes
  # gained when a participant transitions from Not_UTI to UTI?
  a2u <- full_trans %>% filter(transition_type == "Not_UTI\u2192UTI")
  sl("Not_UTI→UTI transitions: %d", nrow(a2u))
  if (nrow(a2u) > 0) {
    sl("  Median Jaccard: %.3f", median(a2u$jaccard_similarity, na.rm = TRUE))
    sl("  Median gained: %.0f, lost: %.0f", median(a2u$n_gained), median(a2u$n_lost))
    if (any(a2u$n_gained > 0)) {
      top_gained <- a2u %>%
        filter(genes_gained != "") %>%
        pull(genes_gained) %>%
        strsplit(";") %>% unlist() %>%
        table() %>% sort(decreasing = TRUE)
      sl("  Most commonly gained: %s",
         paste(head(names(top_gained), 5), collapse = ", "))
    }
  }

  sl("")
  sl("Transition type breakdown:")
  s <- results[["all"]]$summary
  for (i in seq_len(nrow(s))) {
    r <- s[i, ]
    sl("  %s: n=%d, median_jaccard=%.3f, %%noΔ=%.0f%%",
       r$transition_type, r$n_transitions, r$median_jaccard, r$pct_no_change)
  }
} else {
  sl("No transitions computed (no participants with ≥2 timepoints).")
}

writeLines(summary_lines, file.path(DIR_VF, "vf_longitudinal_summary.txt"))
msg("Summary written to %s", file.path(DIR_VF, "vf_longitudinal_summary.txt"))

msg("✓ 24_vf_longitudinal_dynamics.R complete.")

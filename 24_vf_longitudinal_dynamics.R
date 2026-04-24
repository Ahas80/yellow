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
#   for phenotype-switch pairs (same strain, status goes ASB→UTI or vice versa).
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
#     4. Label the transition type (e.g., ASB→ASB, ASB→UTI)
#
# INPUTS:
#   - results/vf/vf_analysis_ready.csv   (from 22_vf_build_analysis_dataset.R)
#   - results/vf/gene_map.csv            (for gene metadata)
#
# OUTPUTS (all in results/vf/):
#   - vf_longitudinal_transitions.csv      One row per consecutive pair
#                                          (full cohort)
#   - vf_transition_summary_by_type.csv    Aggregated by transition type
#   - vf_transitions_stratified.csv        All pairs, by depth (≥2/≥3/≥4)
#   - vf_transition_summary_stratified.csv Summaries, by depth
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
source("R/plot_helpers.R")  # Canonical infection status colours
suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
})

msg("Starting 24_vf_longitudinal_dynamics.R")

# ==============================================================================
# 1. LOAD DATA
# ==============================================================================

ready_file <- file.path(DIR_VF, "vf_analysis_ready.csv")
if (!file.exists(ready_file)) stop("Missing ", ready_file, ". Run 22_vf_build_analysis_dataset.R first.")
vf_ready <- read_csv(ready_file, show_col_types = FALSE) %>%
  mutate(Participant_id = as.character(Participant_id))

gene_map <- read_csv(file.path(DIR_VF, "gene_map.csv"), show_col_types = FALSE) %>%
  mutate(Gene = as.character(Gene),
         Category = coalesce(as.character(Category), "Unassigned"))

# Identify gene columns (exclude metadata and category counts)
skip_cols <- c("Participant_id", "tp_lab", "Infection_Status", "Batch", "ST",
               "vf_count_total", "is_ecoli", "n_timepoints")
cat_cols  <- grep("^cat_", names(vf_ready), value = TRUE)
gene_cols <- setdiff(names(vf_ready), c(skip_cols, cat_cols))

# Define the intended temporal ordering of timepoints.
# T0 is the baseline visit, T1–T4 are follow-ups, and Uricult events
# (triggered by suspected UTI) are placed last.
tp_order <- c("T0", "T1", "T2", "T3", "T4", "Uricult")

msg("Loaded: %d rows, %d gene columns", nrow(vf_ready), length(gene_cols))

# ==============================================================================
# 2. TRANSITION-BUILDING FUNCTION
# ==============================================================================
# This function takes a filtered dataset and builds all consecutive-pair
# transitions for participants with ≥2 timepoints.

build_transitions <- function(data, gene_cols, tp_order, cohort_label = "all") {

  # Filter to episodes with clinical status that fall within our known
  # timepoint ordering, then sort by participant and timepoint rank.
  long_data <- data %>%
    filter(!is.na(Infection_Status), tp_lab %in% tp_order) %>%
    mutate(tp_rank = match(tp_lab, tp_order)) %>%
    arrange(Participant_id, tp_rank)

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
      arrange(tp_rank)

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

      transitions[[length(transitions) + 1]] <- tibble(
        Participant_id  = pid,
        tp_from         = as.character(rf$tp_lab),
        tp_to           = as.character(rt$tp_lab),
        status_from     = rf$Infection_Status,
        status_to       = rt$Infection_Status,
        transition_type = paste0(rf$Infection_Status, "\u2192", rt$Infection_Status),
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

  # Aggregate by transition type (e.g., ASB→ASB, ASB→UTI, etc.)
  # This summary is useful for comparing VF stability across different
  # clinical trajectory types.
  trans_summary <- trans_df %>%
    group_by(transition_type) %>%
    summarise(
      n_transitions    = n(),
      n_participants   = n_distinct(Participant_id),
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
  build_transitions(cohorts[[label]], gene_cols, tp_order, label)
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

write_csv(all_trans, file.path(DIR_VF, "vf_transitions_stratified.csv"))
write_csv(all_summ,  file.path(DIR_VF, "vf_transition_summary_stratified.csv"))

# ==============================================================================
# 5. PLOT: JACCARD SIMILARITY BY TRANSITION TYPE
# ==============================================================================
# This figure shows how stable VF profiles are within each type of clinical
# transition (e.g., ASB→ASB vs ASB→UTI).  High Jaccard values (close to 1)
# indicate the VF profile barely changes between timepoints.

ensure_dir(DIR_PLOTS_VF)

full_trans <- results[["all"]]$transitions
if (nrow(full_trans) > 0) {
  # Group rare transition types (those with < 3 occurrences) into "Other"
  # to keep the plot readable
  type_counts <- full_trans %>% count(transition_type, sort = TRUE)
  common_types <- type_counts %>% filter(n >= 3) %>% pull(transition_type)
  plot_trans <- full_trans %>%
    mutate(plot_type = ifelse(transition_type %in% common_types,
                              transition_type, "Other"))

  p_jac <- ggplot(plot_trans, aes(x = reorder(plot_type, -jaccard_similarity, FUN = median),
                                   y = jaccard_similarity)) +
    geom_boxplot(outlier.shape = NA, width = 0.5, fill = "grey90") +
    geom_jitter(width = 0.15, alpha = 0.5, size = 1.5, colour = "#0072B2") +
    geom_hline(yintercept = 1.0, linetype = "dashed", colour = "grey50") +
    labs(title = "VF Profile Stability Across Consecutive Timepoints",
         subtitle = sprintf("%d transitions from %d participants",
                            nrow(full_trans), n_distinct(full_trans$Participant_id)),
         x = "Transition Type",
         y = "Jaccard Similarity (1 = identical VF profile)") +
    theme_minimal(base_size = 11) +
    theme(axis.text.x = element_text(angle = 30, hjust = 1))

  ggsave(file.path(DIR_PLOTS_VF, "vf_jaccard_by_transition.png"), p_jac,
         width = 8, height = 5, dpi = 300)
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
  sl("Median genes gained: %d", median(full_trans$n_gained))
  sl("Median genes lost: %d",   median(full_trans$n_lost))
  sl("")

  # Special focus on ASB→UTI transitions: are specific genes
  # gained when a participant transitions from ASB to UTI?
  a2u <- full_trans %>% filter(transition_type == "ASB\u2192UTI")
  sl("ASB→UTI transitions: %d", nrow(a2u))
  if (nrow(a2u) > 0) {
    sl("  Median Jaccard: %.3f", median(a2u$jaccard_similarity, na.rm = TRUE))
    sl("  Median gained: %d, lost: %d", median(a2u$n_gained), median(a2u$n_lost))
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

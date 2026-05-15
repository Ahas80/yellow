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
# VF visualisation shared helpers
# ==============================================================================

STATUS_LEVELS <- c("ASB", "UTI", "Negative", "Culture-positive/S&S unknown", "Unknown")

status_for_plot <- function(x) {
  x <- as.character(x)
  x[is.na(x) | x == ""] <- "Unknown"
  x[!x %in% STATUS_LEVELS] <- "Culture-positive/S&S unknown"
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
    "Same-ST and different-ST comparisons are diagnostic for lineage/strain replacement, not causal tests.",
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

# ==============================================================================
# 2. TRANSITION-BUILDING FUNCTION
# ==============================================================================
# This function takes a filtered dataset and builds all consecutive-pair
# transitions for participants with ≥2 timepoints.

build_transitions <- function(data, gene_cols, cohort_label = "all") {

  # Filter to episodes with clinical status, then sort by true collection date
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
        ST_from = st_from,
        ST_to = st_to,
        same_ST = ifelse(!is.na(st_from) & !is.na(st_to), st_from == st_to, NA),
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
      n_with_days_between = sum(!is.na(days_between_samples)),
      median_days_between = ifelse(all(is.na(days_between_samples)), NA_real_, median(days_between_samples, na.rm = TRUE)),
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

write_csv(all_trans, file.path(DIR_VF, "vf_transitions_stratified.csv"))
write_csv(all_summ,  file.path(DIR_VF, "vf_transition_summary_stratified.csv"))

# ==============================================================================
# 5. PLOT: JACCARD SIMILARITY BY TRANSITION TYPE
# ==============================================================================
# This figure shows how stable VF profiles are within each type of clinical
# transition (e.g., ASB→ASB vs ASB→UTI). High Jaccard values (close to 1)
# indicate the VF profile barely changes between consecutive observed isolates.

ensure_dir(DIR_PLOTS_VF)

full_trans <- results[["all"]]$transitions
if (nrow(full_trans) > 0) {
  # Group rare transition types (those with < 3 occurrences) into "Other"
  # to keep the plot readable
  type_counts <- full_trans %>% count(transition_type, sort = TRUE)
  common_types <- type_counts %>% filter(n >= 3) %>% pull(transition_type)
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
    filter(!is.na(same_ST)) %>%
    mutate(ST_comparison = ifelse(same_ST, "Same ST", "Different ST"))
  if (nrow(st_plot) > 0 && n_distinct(st_plot$ST_comparison) >= 1) {
    p_st <- ggplot(st_plot, aes(x = ST_comparison, y = jaccard_similarity,
                                fill = ST_comparison)) +
      geom_boxplot(outlier.shape = NA, alpha = 0.6, width = 0.55) +
      geom_jitter(width = 0.12, alpha = 0.45, size = 1.3) +
      scale_fill_manual(values = c("Same ST" = "#009E73", "Different ST" = "#D55E00")) +
      labs(
        title = "Within-resident VF similarity by sequence-type consistency",
        subtitle = "Same-ST comparisons help distinguish persistent lineage from replacement-like changes",
        x = NULL,
        y = "Jaccard similarity",
        caption = vf_caption(
          ready_file, st_plot, "consecutive within-resident isolate-pair comparison stratified by ST",
          extra_note = "ST agreement is a lineage diagnostic; it does not prove the same strain without SNP/ANI context."
        )
      ) +
      plot_theme_vf(base_size = 11) +
      theme(legend.position = "none")

    ggsave(file.path(DIR_PLOTS_VF, "vf_jaccard_same_vs_different_st.png"),
           p_st, width = 6.5, height = 5.2, dpi = 300)
  } else {
    msg("Skipping same-ST vs different-ST Jaccard plot: insufficient ST-labelled pairs.")
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

#!/usr/bin/env Rscript
# ==============================================================================
# poster_01_reconcile_numbers.R
# Purpose: Reconcile abstract claims against current data files.
# Produces a line-by-line verification report.
# Run from: UTIGA/scripts/
# ==============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
  library(stringr)
  library(glue)
})

# --- PATHS (relative to UTIGA/) ---
DATA_DIR    <- file.path(dirname(getwd()), "data")  # UTIGA/data
OUT_DIR     <- file.path(dirname(getwd()), "reference")
status_file <- file.path(DATA_DIR, "status_map.csv")
mlst_file   <- file.path(DATA_DIR, "mlst_with_meta.csv")
class_file  <- file.path(DATA_DIR, "class_inputs_full.csv")

cat("=== POSTER NUMBER RECONCILIATION ===\n\n")

# --- 1. Load Data ---
status_df <- read_csv(status_file, show_col_types = FALSE) %>%
  mutate(Participant_id = as.character(Participant_id))

mlst_df <- read_csv(mlst_file, show_col_types = FALSE) %>%
  mutate(Participant_id = as.character(Participant_id))

class_inputs <- read_csv(class_file, show_col_types = FALSE) %>%
  mutate(Participant_id = as.character(Participant_id)) %>%
  select(Participant_id, tp_lab, isolate_ID)

cat(glue("Status map rows: {nrow(status_df)}\n"))
cat(glue("MLST rows: {nrow(mlst_df)}\n\n"))

# --- 2. Count participants by timepoint depth ---
tp_counts <- status_df %>%
  filter(!is.na(Infection_Status)) %>%
  group_by(Participant_id) %>%
  summarise(n_tp = n_distinct(Timepoint), .groups = "drop")

for (K in c(2, 3, 4)) {
  n <- sum(tp_counts$n_tp >= K)
  cat(glue("Participants with >= {K} timepoints: {n}\n"))
}

# --- 3. Build linked dataset (same as get_stratified_stats.R) ---
clean_df <- status_df %>%
  filter(!is.na(Infection_Status)) %>%
  group_by(Participant_id) %>%
  mutate(n_tp = n_distinct(Timepoint)) %>%
  ungroup() %>%
  filter(cfu_recorded_any == TRUE)

calc_poster_stats <- function(df, min_tps, label) {
  cohort <- df %>% filter(n_tp >= min_tps)
  pids <- unique(cohort$Participant_id)
  n_pids <- length(pids)
  n_episodes <- nrow(cohort)

  # Link to MLST via class_inputs
  cohort_isolates <- cohort %>%
    left_join(class_inputs, by = c("Participant_id", "Timepoint" = "tp_lab"),
              relationship = "many-to-many")

  cohort_mlst <- cohort_isolates %>%
    left_join(mlst_df %>% select(Isolate_ID, ST, Participant_id, Timepoint),
              by = c("isolate_ID" = "Isolate_ID"),
              relationship = "many-to-many",
              suffix = c("", ".mlst")) %>%
    filter(!is.na(ST)) %>%
    distinct(Participant_id, Timepoint, .keep_all = TRUE) %>%
    mutate(
      tp_num_val = case_when(
        grepl("T[0-9]+", Timepoint) ~ as.numeric(gsub("T", "", Timepoint)),
        Timepoint == "Uricult" ~ 99,
        TRUE ~ 999
      )
    )

  # Total isolates / episodes with ST
  n_with_st <- nrow(cohort_mlst)

  # Typeable vs untypeable
  n_untypeable <- sum(cohort_mlst$ST == "-", na.rm = TRUE)
  n_typeable <- n_with_st - n_untypeable

  # Distinct STs (excluding "-")
  typeable_sts <- cohort_mlst %>% filter(ST != "-") %>% pull(ST)
  n_distinct_sts <- n_distinct(typeable_sts)

  # Top STs (typeable only)
  st_tab <- table(typeable_sts)
  st_df <- as.data.frame(st_tab) %>%
    rename(ST = typeable_sts, n = Freq) %>%
    arrange(desc(n)) %>%
    mutate(pct = round(n / sum(n) * 100, 1))

  # UTI episodes
  uti_mlst <- cohort_mlst %>% filter(Infection_Status == "UTI", ST != "-")
  n_uti_typeable <- nrow(uti_mlst)
  uti_st_tab <- as.data.frame(table(uti_mlst$ST)) %>%
    rename(ST = Var1, n = Freq) %>%
    arrange(desc(n))

  # --- Pair-level stability (matching abstract methodology) ---
  st_pairs <- cohort_mlst %>%
    filter(ST != "-") %>%
    distinct(Participant_id, tp_num_val, ST, Infection_Status) %>%
    group_by(Participant_id) %>%
    arrange(tp_num_val) %>%
    mutate(Next_ST = lead(ST)) %>%
    filter(!is.na(Next_ST)) %>%
    ungroup()

  n_pairs <- nrow(st_pairs)
  n_same <- sum(st_pairs$ST == st_pairs$Next_ST)
  pct_same <- if (n_pairs > 0) round(n_same / n_pairs * 100, 1) else NA

  # Participant-level stability
  part_stability <- st_pairs %>%
    group_by(Participant_id) %>%
    summarise(any_switch = any(ST != Next_ST), .groups = "drop")
  n_switchers <- sum(part_stability$any_switch)
  pct_switch <- if (nrow(part_stability) > 0) round(n_switchers / nrow(part_stability) * 100, 1) else NA

  # UTI residents with stability info (for >=3 TP)
  uti_pids <- cohort_mlst %>%
    filter(Infection_Status == "UTI") %>%
    pull(Participant_id) %>%
    unique()

  uti_stability <- st_pairs %>%
    filter(Participant_id %in% uti_pids) %>%
    group_by(Participant_id) %>%
    summarise(any_switch = any(ST != Next_ST), .groups = "drop")

  n_uti_stable <- sum(!uti_stability$any_switch)
  n_uti_switched <- sum(uti_stability$any_switch)

  # --- OUTPUT ---
  cat(glue("\n\n{'='} {label} {'='}\n"))
  cat(glue("  Clinical participants: {n_pids}\n"))
  cat(glue("  Episodes with CFU: {n_episodes}\n"))
  cat(glue("  Episodes with ST linked: {n_with_st}\n"))
  cat(glue("  Untypeable ('-'): {n_untypeable} ({round(n_untypeable/n_with_st*100,1)}%)\n"))
  cat(glue("  Typeable isolates: {n_typeable}\n"))
  cat(glue("  Distinct STs (typeable): {n_distinct_sts}\n"))
  cat("\n  Top 5 STs:\n")
  print(head(st_df, 5))
  cat(glue("\n  UTI typeable isolates: {n_uti_typeable}\n"))
  cat("  UTI ST distribution:\n")
  print(head(uti_st_tab, 5))
  cat(glue("\n  PAIR-LEVEL STABILITY:\n"))
  cat(glue("    Consecutive pairs: {n_pairs}\n"))
  cat(glue("    Same ST: {n_same} ({pct_same}%)\n"))
  cat(glue("  PARTICIPANT-LEVEL:\n"))
  cat(glue("    Participants with >=2 linked: {nrow(part_stability)}\n"))
  cat(glue("    Switchers: {n_switchers} ({pct_switch}%)\n"))
  cat(glue("  UTI RESIDENTS (with ST pairs):\n"))
  cat(glue("    Total: {nrow(uti_stability)}\n"))
  cat(glue("    Stable: {n_uti_stable}\n"))
  cat(glue("    Switched: {n_uti_switched}\n"))
}

calc_poster_stats(clean_df, 2, "COHORT >= 2 TIMEPOINTS")
calc_poster_stats(clean_df, 3, "COHORT >= 3 TIMEPOINTS")
calc_poster_stats(clean_df, 4, "COHORT >= 4 TIMEPOINTS")

cat("\n\n=== ABSTRACT COMPARISON ===\n")
cat("Abstract claims (submitted):\n")
cat("  >= 2 TP: 83 participants, 181 isolates, 38 STs, 31 untypeable (17.1%)\n")
cat("  Top STs: ST43 (11.6%), ST4 (7.7%), ST6 (6.6%), ST1 (5.5%)\n")
cat("  UTI typeable: 33, ST6=8, ST10=5\n")
cat("  Stability: 80.2% (77/96 pairs)\n")
cat("  >= 3 TP: 54 participants, 119 isolates\n")
cat("  Stability: 81.8% (54/66 pairs)\n")
cat("  Switching: ~20%\n")
cat("  UTI >=3 TP: 12 residents, 8 stable, 4 switched\n")
cat("\nCompare these against the numbers printed above.\n")
cat("Discrepancies may be due to:\n")
cat("  1) Different filtering (cfu_recorded_any vs all episodes)\n")
cat("  2) Different counting (episodes vs isolates vs linked pairs)\n")
cat("  3) Whether '-' STs count as switches\n")

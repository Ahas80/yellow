#!/usr/bin/env Rscript
# ==============================================================================
# 00d_derive_plot_timepoints.R
# ------------------------------------------------------------------------------
# Role: [Data Prep / Plotting-Only] — Derive poster-safe half-step timepoint
#       labels for Uricult events.
#
# ⚠️  WARNING: The variables created by this script are POSTER / DISPLAY
#     APPROXIMATIONS ONLY.  They are NOT to be used as statistical covariates
#     or modelling variables.  For modelling, derive Days_since_T0 (continuous
#     elapsed time) from Collection_Date.
#
# Inputs:
#   - results/clinical/status_map.csv          (from 00b_classify_episodes.R)
#   - data/inputs/batch1.csv, batch2.csv, batch3.csv  (for Collection_Date)
#
# Outputs:
#   - results/clinical/status_map_with_poster_tp.csv
#     (copy also placed in UTIGA/data/ if that directory exists)
#
# Does NOT modify:
#   - Original Timepoint column
#   - Original status_map.csv
#
# Usage:
#   Rscript 00d_derive_plot_timepoints.R
# ==============================================================================

source("00_config.R")
source(here::here("R", "clinical_helpers.R"))

suppressPackageStartupMessages({
    library(dplyr)
    library(readr)
    library(stringr)
    library(lubridate)
})

msg("Starting 00d_derive_plot_timepoints.R")

# ==============================================================================
# 1. Load Status Map
# ==============================================================================
status_file <- file.path(DIR_CLINICAL, "status_map.csv")
if (!file.exists(status_file)) {
    stop("status_map.csv not found. Run 00b_classify_episodes.R first.")
}
status_map <- read_csv(status_file, show_col_types = FALSE) %>%
    mutate(Participant_id = as.character(Participant_id))

msg("Loaded %d episodes from status_map.csv", nrow(status_map))

# ==============================================================================
# 2. Recover Collection_Date from Batch CSVs
# ==============================================================================
# status_map.csv does not carry Collection_Date (lost in group_by collapse).
# We recover it from the raw batch files.
read_batch_dates <- function(fname) {
    p <- file.path("data", "inputs", fname)
    if (!file.exists(p)) {
        warning("Batch file not found: ", p)
        return(tibble(Participant_id = character(),
                      Timepoint = character(),
                      Collection_Date_raw = character()))
    }
    read_csv(p, show_col_types = FALSE) %>%
        transmute(
            Participant_id = as.character(Participant_id),
            Timepoint      = canon_tp(Timepoint),
            Collection_Date_raw = as.character(Collection_Date)
        )
}

batch_dates <- bind_rows(
    read_batch_dates("batch1.csv"),
    read_batch_dates("batch2.csv"),
    read_batch_dates("batch3.csv")
) %>%
    mutate(
        date_parsed = lubridate::dmy(Collection_Date_raw)
    )

# Take the earliest date per (Participant_id, Timepoint) episode
episode_dates <- batch_dates %>%
    filter(!is.na(date_parsed)) %>%
    group_by(Participant_id, Timepoint) %>%
    summarise(
        Collection_Date     = min(Collection_Date_raw, na.rm = TRUE),
        Collection_Date_dmy = min(date_parsed, na.rm = TRUE),
        .groups = "drop"
    )

msg("Recovered dates for %d episodes from batch CSVs", nrow(episode_dates))

# Join dates into status_map
status_aug <- status_map %>%
    left_join(episode_dates, by = c("Participant_id", "Timepoint"))

n_with_date <- sum(!is.na(status_aug$Collection_Date_dmy))
msg("Joined dates: %d/%d episodes have Collection_Date", n_with_date, nrow(status_aug))

# ==============================================================================
# 3. Build Routine-Visit Lookup Per Participant
# ==============================================================================
# Extract numeric timepoint index for routine visits (T0=0, T1=1, …)
routine_visits <- status_aug %>%
    filter(grepl("^T[0-9]+$", Timepoint), !is.na(Collection_Date_dmy)) %>%
    mutate(tp_int = as.integer(sub("^T", "", Timepoint))) %>%
    select(Participant_id, tp_int, tp_date = Collection_Date_dmy) %>%
    distinct() %>%
    arrange(Participant_id, tp_int)

msg("Built routine-visit lookup: %d participant-timepoint pairs", nrow(routine_visits))

# ==============================================================================
# 4. Assign Poster Half-Step Labels
# ==============================================================================
# Rule:
#   For routine visits → integer (T0=0, T1=1, etc.)
#   For Uricults:
#     A) Both flanking routine visits known:
#        - Consecutive flanking (TX → T(X+1)): assign TX.5
#        - Skipped visits (TX → T(X+N)): assign based on timing
#          · ≤90 days after prior → T{prior}.5
#          · >90 days after prior → T{prior + floor(days/90)}.5
#          · Confidence = "High"
#     B) Only prior routine visit known, ≤90 days away:
#        - Assign T{prior}.5, Confidence = "Moderate"
#     C) Only next routine visit known, ≤90 days away, next ≥ T1:
#        - Assign T{next-1}.5, Confidence = "Moderate"
#     D) All other cases (>90d single-sided, unlinkable, no routine visits):
#        - Assign NA, Confidence = "Excluded"

MAX_SINGLE_SIDE_DAYS <- 90

assign_poster_tp <- function(df, routine_lkp) {
    # Pre-allocate output columns
    n <- nrow(df)
    plot_label <- character(n)
    plot_num   <- numeric(n)
    confidence <- character(n)

    for (i in seq_len(n)) {
        tp  <- df$Timepoint[i]
        pid <- df$Participant_id[i]
        dt  <- df$Collection_Date_dmy[i]

        # --- Routine visit: straightforward integer ---
        if (grepl("^T[0-9]+$", tp)) {
            tp_int <- as.integer(sub("^T", "", tp))
            plot_label[i] <- tp
            plot_num[i]   <- tp_int
            confidence[i] <- "Scheduled"
            next
        }

        # --- Not a Uricult: leave NA ---
        if (!grepl("(?i)uricult", tp)) {
            plot_label[i] <- tp
            plot_num[i]   <- NA_real_
            confidence[i] <- "Excluded"
            next
        }

        # --- Uricult: attempt placement ---
        # Check for unlinkable participant IDs
        pid_lower <- tolower(pid)
        if (pid_lower %in% c("still to be linked", "niet te koppelen", "") ||
            is.na(pid) || is.na(dt)) {
            plot_label[i] <- "Uricult_unplaced"
            plot_num[i]   <- NA_real_
            confidence[i] <- "Excluded"
            next
        }

        # Get this participant's routine visits
        rv <- routine_lkp %>% filter(Participant_id == pid) %>% arrange(tp_int)

        if (nrow(rv) == 0) {
            # No routine visits at all
            plot_label[i] <- "Uricult_unplaced"
            plot_num[i]   <- NA_real_
            confidence[i] <- "Excluded"
            next
        }

        # Find flanking visits
        before <- rv %>% filter(tp_date <= dt)
        after  <- rv %>% filter(tp_date >  dt)

        has_before <- nrow(before) > 0
        has_after  <- nrow(after) > 0

        if (has_before) {
            prior_row <- before %>% slice_tail(n = 1)
            prior_int  <- prior_row$tp_int
            prior_date <- prior_row$tp_date
            days_after_prior <- as.numeric(difftime(dt, prior_date, units = "days"))
        }

        if (has_after) {
            next_row <- after %>% slice_head(n = 1)
            next_int  <- next_row$tp_int
            next_date <- next_row$tp_date
            days_before_next <- as.numeric(difftime(next_date, dt, units = "days"))
        }

        # --- Case A: Both flanking visits known ---
        if (has_before && has_after) {
            if (next_int - prior_int == 1) {
                # Consecutive (e.g. T0 → T1)
                half <- prior_int + 0.5
                plot_label[i] <- paste0("T", prior_int, ".5")
                plot_num[i]   <- half
                confidence[i] <- "High"
            } else {
                # Skipped visits (e.g. T0 → T2)
                if (days_after_prior <= MAX_SINGLE_SIDE_DAYS) {
                    half <- prior_int + 0.5
                    plot_label[i] <- paste0("T", prior_int, ".5")
                    plot_num[i]   <- half
                    confidence[i] <- "High"
                } else {
                    # Place in a later sub-interval
                    intervals_past <- floor(days_after_prior / MAX_SINGLE_SIDE_DAYS)
                    est_tp <- min(prior_int + intervals_past, next_int - 1)
                    half <- est_tp + 0.5
                    plot_label[i] <- paste0("T", est_tp, ".5")
                    plot_num[i]   <- half
                    confidence[i] <- "High"
                }
            }
            next
        }

        # --- Case B: Prior only ---
        if (has_before && !has_after) {
            if (days_after_prior <= MAX_SINGLE_SIDE_DAYS) {
                half <- prior_int + 0.5
                plot_label[i] <- paste0("T", prior_int, ".5")
                plot_num[i]   <- half
                confidence[i] <- "Moderate"
            } else {
                plot_label[i] <- "Uricult_unplaced"
                plot_num[i]   <- NA_real_
                confidence[i] <- "Excluded"
            }
            next
        }

        # --- Case C: Next only ---
        if (!has_before && has_after) {
            if (days_before_next <= MAX_SINGLE_SIDE_DAYS && next_int >= 1) {
                half <- next_int - 0.5
                plot_label[i] <- paste0("T", next_int - 1L, ".5")
                plot_num[i]   <- half
                confidence[i] <- "Moderate"
            } else {
                plot_label[i] <- "Uricult_unplaced"
                plot_num[i]   <- NA_real_
                confidence[i] <- "Excluded"
            }
            next
        }

        # Fallback (should not reach here)
        plot_label[i] <- "Uricult_unplaced"
        plot_num[i]   <- NA_real_
        confidence[i] <- "Excluded"
    }

    df %>% mutate(
        Plot_TP_Label_Poster = plot_label,
        Plot_TP_Num_Poster   = plot_num,
        Placement_Confidence = confidence
    )
}

status_out <- assign_poster_tp(status_aug, routine_visits)

# ==============================================================================
# 5. Report Summary
# ==============================================================================
uricult_rows <- status_out %>% filter(grepl("(?i)uricult", Timepoint, perl = TRUE))

msg("\n=== URICULT HALF-STEP ASSIGNMENT SUMMARY ===")
msg("Total Uricult rows: %d", nrow(uricult_rows))
msg("")

# Counts by label
label_counts <- uricult_rows %>%
    count(Plot_TP_Label_Poster, Placement_Confidence, name = "n") %>%
    arrange(Plot_TP_Num_Poster = NA)
msg("By label:")
for (j in seq_len(nrow(label_counts))) {
    msg("  %-20s %-12s  n=%d",
        label_counts$Plot_TP_Label_Poster[j],
        paste0("[", label_counts$Placement_Confidence[j], "]"),
        label_counts$n[j])
}

# Placed vs excluded
n_placed   <- sum(!is.na(uricult_rows$Plot_TP_Num_Poster))
n_excluded <- sum(is.na(uricult_rows$Plot_TP_Num_Poster))
msg("")
msg("Placed:   %d Uricult rows (%d%%)", n_placed, round(n_placed / nrow(uricult_rows) * 100))
msg("Excluded: %d Uricult rows (%d%%)", n_excluded, round(n_excluded / nrow(uricult_rows) * 100))

# Confidence breakdown
conf_counts <- uricult_rows %>% count(Placement_Confidence, name = "n")
msg("")
msg("By confidence:")
for (j in seq_len(nrow(conf_counts))) {
    msg("  %-12s  n=%d", conf_counts$Placement_Confidence[j], conf_counts$n[j])
}

# Sanity: routine visits should be integers
routine_check <- status_out %>%
    filter(grepl("^T[0-9]+$", Timepoint)) %>%
    mutate(is_int = Plot_TP_Num_Poster == round(Plot_TP_Num_Poster))
n_bad <- sum(!routine_check$is_int, na.rm = TRUE)
if (n_bad > 0) {
    warning("SANITY FAILED: ", n_bad, " routine visits have non-integer Plot_TP_Num_Poster")
} else {
    msg("✓ All routine visits have integer Plot_TP_Num_Poster")
}

# ==============================================================================
# 6. Write Outputs
# ==============================================================================
# Select final columns: everything from status_map + new poster columns + date
out_cols <- c(
    names(status_map),
    "Collection_Date",
    "Plot_TP_Label_Poster",
    "Plot_TP_Num_Poster",
    "Placement_Confidence"
)

status_final <- status_out %>% select(any_of(out_cols))

out_path <- file.path(DIR_CLINICAL, "status_map_with_poster_tp.csv")
write_csv(status_final, out_path)
msg("Wrote %s (%d rows, %d cols)", out_path, nrow(status_final), ncol(status_final))

# Also copy to UTIGA/data/ if that directory exists (for poster scripts)
utiga_data <- file.path(DIR_ROOT, "UTIGA", "data")
if (dir.exists(utiga_data)) {
    utiga_path <- file.path(utiga_data, "status_map_with_poster_tp.csv")
    write_csv(status_final, utiga_path)
    msg("Also copied to %s", utiga_path)
}

msg("✓ Done. Original status_map.csv is UNCHANGED.")
msg("")
msg("⚠️  REMINDER: Plot_TP_Label_Poster and Plot_TP_Num_Poster are")
msg("   POSTER / DISPLAY APPROXIMATIONS ONLY.")
msg("   Do NOT use them as statistical covariates.")

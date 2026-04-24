# Load libraries
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(tidyr))

# Paths
status_path <- "results/class_inputs_full.csv"
pairwise_path <- "results/pairwise_stats.csv"
mlst_path <- "results/mlst/mlst_with_meta.csv"

# 1. Load Data
status_df <- read_csv(status_path, show_col_types = FALSE) %>% rename(isolate_id = isolate_ID, pid = Participant_id)
mlst_df <- read_csv(mlst_path, show_col_types = FALSE) %>% rename(isolate_id = Isolate_ID)
pairwise_df <- read_csv(pairwise_path, show_col_types = FALSE)

# 2. Strict Deduplication & Joining
# Status may have duplicates. Group by unique Episode (PID+Timepoint+Isolate)
status_clean <- status_df %>%
    filter(!is.na(isolate_id)) %>%
    distinct(pid, tp_num, Infection_Status, isolate_id, .keep_all = TRUE)

# Join MLST (Inner Join creates the "Complete" set)
complete_df <- status_clean %>%
    inner_join(mlst_df %>% select(isolate_id, ST, scheme, everything()), by = "isolate_id") %>%
    # Deduplicate again just in case MLST had duplicates
    distinct(pid, tp_num, Infection_Status, isolate_id, .keep_all = TRUE)

# 3. Corrected Population Counts (Unique Episodes)
n_participants <- n_distinct(complete_df$pid)
episodes <- complete_df %>% distinct(pid, tp_num, Infection_Status)
n_asb <- sum(episodes$Infection_Status == "ASB")
n_uti <- sum(episodes$Infection_Status == "UTI")

cat("\n--- CORRECTED POPULATION COUNTS ---\n")
cat("Participants:", n_participants, "\n")
cat("ASB Episodes:", n_asb, "\n")
cat("UTI Episodes:", n_uti, "\n")

# 4. Re-check Overlap (ASB+UTI or Neg+UTI)
# Need to include Negative status for the user's "transitions" request, if present in complete_df
# But standard analysis is ASB vs UTI.
# Let's check PID overlap in complete_df
pids_complete <- complete_df %>%
    group_by(pid) %>%
    summarize(
        has_asb = any(Infection_Status == "ASB"),
        has_uti = any(Infection_Status == "UTI"),
        has_neg = any(Infection_Status == "Negative")
    )

overlap_asb_uti <- pids_complete %>%
    filter(has_asb & has_uti) %>%
    pull(pid)
overlap_neg_uti <- pids_complete %>%
    filter(has_neg & has_uti) %>%
    pull(pid)

cat("\n--- GENOMICALLY COMPLETE TRANSITIONS ---\n")
cat("Participants with ASB+UTI (Complete):", length(overlap_asb_uti), "\n")
if (length(overlap_asb_uti) > 0) print(overlap_asb_uti)
cat("Participants with Neg+UTI (Complete):", length(overlap_neg_uti), "\n")
if (length(overlap_neg_uti) > 0) print(overlap_neg_uti)

# 5. Extract Pairs for these overlaps
target_pids <- unique(c(overlap_asb_uti, overlap_neg_uti))

if (length(target_pids) > 0) {
    # Logic: For each PID, get all valid isolates, calculate similarity.
    # Use pairwise_stats if available, else MLST.
    # We'll just list the isolates for manual verification script output first.

    for (p in target_pids) {
        cat("\nParticipant", p, "Isolates:\n")
        iso <- complete_df %>%
            filter(pid == p) %>%
            select(tp_num, Infection_Status, isolate_id, ST)
        print(iso)

        # Check SNP distance in pairwise
        # Filter pairwise for this PID
        p_pairs <- pairwise_df %>% filter(Participant_id == p)
        if (nrow(p_pairs) > 0) {
            cat("  SNP Distances found:\n")
            # Print relevant columns
            print(p_pairs %>% select(Isolate_ID, path_B, TotalSnpCnt) %>% head())
        }
    }
}

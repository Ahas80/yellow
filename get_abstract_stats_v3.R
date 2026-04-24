# Script: get_abstract_stats_v3.R (Corrected Source)
# Goal: Generate robust, verifiable statistics for the R- UTI Abstract
# Features: Deduplication, Cross-File Verification, Timepoint Filtering

suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(stringr))

# --- USER CONFIGURATION (CHANGE THESE SETTINGS) ---

# 1. TIMEPOINT FILTERING
# How many "Complete" episodes (Clinical + Isolate + Genotype) must a participant have to be included?
# - Set to 1 to include almost everyone who has data.
# - Set to 2 or more to find longitudinal participants only.
MIN_COMPLETE_EPISODES <- 1

# How many total clinical timepoints (including empty/negative ones) must they have?
# - Set to 0 to ignore this check.
MIN_CLINICAL_TIMEPOINTS <- 0

# 2. STATUS INCLUSION (URICULT / NEGATIVE)
# Which clinical statuses do you want to include in the "Complete" counts?
# - "ASB", "UTI": Standard infections.
# - "Negative": Usually screening samples or Uricult negatives.
# - Remove "Negative" from this list if you want to EXCLUDE Uricult results from the final stats.
INCLUDED_STATUSES <- c("ASB", "UTI", "Negative")
# Example to EXCLUDE Uricult:
# INCLUDED_STATUSES <- c("ASB", "UTI")

# Paths
# status_map_path REMOVED - Found to be incomplete subset.
class_inputs_path <- "results/class_inputs_full.csv" # The Comprehensive Source
mlst_path <- "results/mlst/mlst_with_meta.csv" # Genotyped Truth
pairwise_path <- "results/pairwise_stats.csv" # SNP Distances

# --- 1. DATA LOADING ---
message("Loading Data...")
class_inputs <- read_csv(class_inputs_path, show_col_types = FALSE) %>%
    mutate(Participant_id = as.character(Participant_id), tp_num = as.character(tp_num))

mlst_df <- read_csv(mlst_path, show_col_types = FALSE) %>%
    mutate(Participant_id = as.character(Participant_id))

pairwise_df <- read_csv(pairwise_path, show_col_types = FALSE) %>%
    mutate(Participant_id = as.character(Participant_id))

# --- 2. DATA CLEANING & MERGING (THE FUNNEL) ---

# A. Clinical/Sequencing Base (class_inputs)
# This file contains one row per isolate (or per episode if culture negative).
# We treat this as the base.
base_df <- class_inputs %>%
    rename(pid = Participant_id, isolate_id = isolate_ID, status_clinical = Infection_Status) %>%
    # Deduplicate: Unique Episode + Isolate Combo
    distinct(pid, tp_num, isolate_id, .keep_all = TRUE)

# B. Genomic Data (MLST)
gen_data <- mlst_df %>%
    rename(isolate_id = Isolate_ID) %>%
    select(isolate_id, ST, scheme) %>%
    distinct(isolate_id, .keep_all = TRUE) # Unique Genotype per Isolate

# --- 3. BUILD MASTER TABLE ---

master_df <- base_df %>%
    # Join Genomic Info
    left_join(gen_data, by = "isolate_id") %>%
    mutate(
        has_clinical = TRUE,
        has_isolate = !is.na(isolate_id), # Some rows have NA isolate_id (Cult Neg)
        has_genome = !is.na(ST),
        is_complete_episode = has_clinical & has_isolate & has_genome
    )

# --- 4. VERIFICATION & TRANSITION CHECK ---

message("\n--- TRANSITION VERIFICATION (AUDIT) ---")
check_transitions <- function(df, pids, note) {
    sub <- df %>%
        filter(pid %in% pids) %>%
        select(pid, tp_lab, status_clinical, isolate_id, ST, has_genome)
    cat("\nCheck:", note, "\n")
    print(sub)
}
check_transitions(master_df, "122006", "Pt 122006")
check_transitions(master_df, "120003", "Pt 120003")

# --- 5. FILTERING ---

pid_stats <- master_df %>%
    group_by(pid) %>%
    summarize(
        n_clinical_tp = n_distinct(tp_num),
        n_complete_episodes = n_distinct(paste(tp_num, isolate_id)[is_complete_episode]),
        has_asb_complete = any(status_clinical == "ASB" & is_complete_episode),
        has_uti_complete = any(status_clinical == "UTI" & is_complete_episode)
    )

valid_pids <- pid_stats %>%
    filter(n_complete_episodes >= MIN_COMPLETE_EPISODES) %>%
    filter(n_clinical_tp >= MIN_CLINICAL_TIMEPOINTS) %>%
    pull(pid)

final_df <- master_df %>%
    filter(pid %in% valid_pids) %>%
    filter(is_complete_episode) %>%
    filter(status_clinical %in% INCLUDED_STATUSES)

# --- 6. STATISTICS GENERATION ---

message("\n--- FINAL STATISTICS ---")
message(paste("Filter: Min Complete Episodes =", MIN_COMPLETE_EPISODES))

# Population N
n_final_pids <- n_distinct(final_df$pid)
cat("Total Participants (Filtered):", n_final_pids, "\n")

# Episodes Details
# Deduplicate unique episodes (Ignore Isolate dups for count unless different status)
unique_episodes <- final_df %>% distinct(pid, tp_num, status_clinical)
counts_status <- table(unique_episodes$status_clinical)
print(counts_status)

# Top STs
cat("\n--- TOP STs ---\n")
top_st <- final_df %>%
    group_by(status_clinical, ST) %>%
    tally() %>%
    arrange(status_clinical, desc(n)) %>%
    group_by(status_clinical) %>%
    slice_head(n = 5)
print(top_st)

# --- 7. EXPORT VERIFIABLE TABLE ---
verify_table <- pid_stats %>%
    filter(pid %in% valid_pids) %>%
    select(pid, n_clinical_tp, n_complete_episodes) %>%
    arrange(desc(n_complete_episodes))

write_csv(verify_table, "results/verification_table_v3.csv")
message("\nVerification table written to 'results/verification_table_v3.csv'")

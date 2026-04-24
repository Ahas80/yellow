
# Script: verify_clinical_stats.R
# Goal: Attempt to replicate the "Clinical Abstract" numbers from the User's Image
# Image Claims:
# - N = 71 participants (with >= 2 timepoints)
# - Episodes = 236
# - ASB = 187 (79.2%)
# - UTI = 17 (7.2%)
# - Negative = 32 (13.6%)
# - Transitions: 14 ASB -> UTI

suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(tidyr))

# Load Clinical Data
# Trying status_map_full.csv as it usually represents the full clinical picture
status_df <- read_csv("results/status_map_full.csv", show_col_types = FALSE) %>%
  mutate(Participant_id = as.character(Participant_id))

# Clean Statuses to match Abstract categories
# Map to: ASB, UTI, Negative
# Note: Check exact spelling in file
clean_df <- status_df %>%
  mutate(
    Status_Simple = case_when(
      Infection_Status == "ASB" ~ "ASB",
      Infection_Status == "UTI" ~ "UTI",
      Infection_Status %in% c("Negative", "Cult_Neg", "Culture_Negative") ~ "Negative",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(Status_Simple)) 

# Filter: >= 2 Timepoints per Participant
valid_pids <- clean_df %>%
  group_by(Participant_id) %>%
  summarize(n_tp = n_distinct(Timepoint)) %>%
  filter(n_tp >= 2) %>%
  pull(Participant_id)

cohort_df <- clean_df %>%
  filter(Participant_id %in% valid_pids)

# --- STATS CHECK ---

message("--- REPLICATION ATTEMPT ---")
cat("Participants:", n_distinct(cohort_df$Participant_id), "(Image: 71)\n")
cat("Episodes:", nrow(cohort_df), "(Image: 236)\n")

cat("\nStatus Counts:\n")
print(table(cohort_df$Status_Simple))
# Image: ASB 187, UTI 17, Neg 32

message("\n--- TRANSITION ANALYSIS ---")
# Need to sort by timepoint to calculate transitions
# Assuming Timepoint column issortable or we can extract numeric
# T0, T1, T2... Uricult?
# Hard to sort mixed strings. Let's try basic alphanumeric sort or extraction.

trans_df <- cohort_df %>%
  mutate(
    tp_num = case_when(
      grepl("T[0-9]+", Timepoint) ~ as.numeric(gsub("T", "", Timepoint)),
      Timepoint == "Uricult" ~ 99,
      TRUE ~ 999
    )
  ) %>%
  group_by(Participant_id) %>%
  arrange(tp_num) %>%
  mutate(
    Next_Status = lead(Status_Simple)
  ) %>%
  filter(!is.na(Next_Status))

cat("Total Transitions:", nrow(trans_df), "(Image: 165)\n")

cat("\nTransitions from ASB:\n")
asb_trans <- trans_df %>% filter(Status_Simple == "ASB")
print(table(asb_trans$Next_Status))
# Image: 14 ASB -> UTI


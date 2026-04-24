suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(readr))

# 1. Load Data
status_df <- read_csv("status_map.csv", show_col_types = FALSE) %>%
    mutate(Participant_id = as.character(Participant_id))

# Recreate Clean DF (Verified Logic)
clean_df <- status_df %>%
    filter(!is.na(Infection_Status)) %>%
    group_by(Participant_id) %>%
    mutate(n_tp = n_distinct(Timepoint)) %>%
    ungroup() %>%
    filter(cfu_recorded_any == TRUE)

cat("Clean DF (CFU+) Rows:", nrow(clean_df), "\n")
cat("Clean DF >= 2 TPs Rows:", nrow(clean_df %>% filter(n_tp >= 2)), "\n")

# 2. Check Deduplication Logic
raw_inputs <- read_csv("results/class_inputs_full.csv", show_col_types = FALSE) %>%
    mutate(Participant_id = as.character(Participant_id))

dedup_inputs <- raw_inputs %>%
    group_by(Participant_id, tp_lab) %>%
    slice(1) %>%
    ungroup()

cat("Dedup Inputs Rows:", nrow(dedup_inputs), "\n")
cat("Duplicated Keys in Dedup?:", nrow(dedup_inputs %>% group_by(Participant_id, tp_lab) %>% filter(n() > 1)), "\n")

# 3. Simulate Join
cohort <- clean_df %>% filter(n_tp >= 2)
joined <- cohort %>%
    left_join(dedup_inputs, by = c("Participant_id", "Timepoint" = "tp_lab"))

cat("Joined Rows:", nrow(joined), "\n")

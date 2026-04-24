# Load libraries
library(dplyr)
library(readr)

df <- read_csv("results/clinical/status_map.csv", show_col_types = FALSE)

# 1. Identify the 19 valid UTIs
utis <- df %>%
    filter(Infection_Status == "UTI") %>%
    filter(Participant_id != "Still to be linked")

# 2. Identify the 2 participants likely excluded (by timepoint count)
excluded_ids <- c("120014", "122005")

# 3. Inspect these 2 specifically
cat("=== Inspecting Excluded Candidates ===\n")
excluded_rows <- utis %>% filter(Participant_id %in% excluded_ids)
print(excluded_rows %>% select(Participant_id, Timepoint, cfu_recorded_any, raw_CFU_examples, raw_BEO_examples, culture_pos_epi, Sx_present_any, cfu_ge_1e5_any) %>% print(width = Inf))

# 4. Confirm they have NO other timepoints in the full dataset
cat("\n=== Checking Timepoint Counts ===\n")
tp_counts <- df %>%
    filter(Participant_id %in% excluded_ids) %>%
    group_by(Participant_id) %>%
    summarise(n_timepoints = n(), timepoints = paste(Timepoint, collapse = ", "))
print(tp_counts)

# 5. Check if any *other* UTIs are weird
cat("\n=== Checking Other UTIs ===\n")
other_utis <- utis %>% filter(!Participant_id %in% excluded_ids)
# print brief summary
print(other_utis %>% select(Participant_id, cfu_recorded_any, culture_pos_epi) %>% summary())

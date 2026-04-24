# Load libraries
library(dplyr)
library(readr)

cat("Reading results/status_map_full.csv...\n")
if (!file.exists("results/status_map_full.csv")) {
    cat("File not found! using results/clinical/status_map.csv instead as fallback.\n")
    df <- read_csv("results/clinical/status_map.csv", show_col_types = FALSE)
} else {
    df <- read_csv("results/status_map_full.csv", show_col_types = FALSE)
}

# 1. Raw UTI count
utis <- df %>% filter(Infection_Status == "UTI")
cat("Raw UTI count:", nrow(utis), "\n")

# 2. Replicate verify_clinical_stats.R logic
valid_pids <- df %>%
    group_by(Participant_id) %>%
    summarize(n_tp = n_distinct(Timepoint)) %>%
    filter(n_tp >= 2) %>%
    pull(Participant_id)

cat("Number of participants with >= 2 timepoints:", length(valid_pids), "\n")

# 3. Filter UTIs by valid_pids
utis_filtered <- utis %>% filter(Participant_id %in% valid_pids)
cat("UTI count after >= 2 timepoints filter:", nrow(utis_filtered), "\n")

# 4. Identify the dropped UTIs
dropped_utis <- utis %>% filter(!Participant_id %in% valid_pids)
cat("Dropped UTIs count:", nrow(dropped_utis), "\n")

if (nrow(dropped_utis) > 0) {
    cat("\n=== Dropped UTIs Details ===\n")
    # Check relevant columns
    print(dropped_utis %>%
        select(Participant_id, Timepoint, Infection_Status) %>%
        print(width = Inf))

    # Check if they have CFU data
    # Note: verify_clinical_stats.R doesn't seem to check CFU, just timepoints.
    # But let's check if the columns exist in this file
    cols <- colnames(dropped_utis)
    cfu_cols <- cols[grepl("cfu|beoord|raw", cols, ignore.case = TRUE)]

    if (length(cfu_cols) > 0) {
        print(dropped_utis %>% select(Participant_id, any_of(cfu_cols)))
    } else {
        cat("CFU columns not found in status_map_full.csv\n")
    }
}

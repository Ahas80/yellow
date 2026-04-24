library(tidyverse)

# Load data
df_status <- read_csv("/Users/Aamir/Desktop/rUTIs/results/clinical/status_map.csv", show_col_types = FALSE)
df_mlst <- read_csv("/Users/Aamir/Desktop/rUTIs/results/mlst/mlst_with_meta.csv", show_col_types = FALSE)

# Clean Status Map
if ("Participant_id" %in% names(df_status) && !"participant_id" %in% names(df_status)) {
    df_status <- df_status %>% rename(participant_id = Participant_id)
}
# Clean MLST like the main script
if ("Participant_id" %in% names(df_mlst) && !"participant_id" %in% names(df_mlst)) {
    df_mlst <- df_mlst %>% rename(participant_id = Participant_id)
}
df_mlst_clean <- df_mlst %>%
    mutate(
        is_complete = if_else(str_detect(tolower(as.character(if ("mlst_complete" %in% names(.)) mlst_complete else "FALSE")), "true|yes"), 1, 0),
        n_loci = if ("n_loci_typed" %in% names(.)) n_loci_typed else 7
    ) %>%
    arrange(participant_id, Timepoint, desc(is_complete), desc(n_loci)) %>%
    distinct(participant_id, Timepoint, .keep_all = TRUE) # Keep best MLST per episode

# Check math for K=2
K <- 2
cat(sprintf("\n--- K=%d Analysis ---\n", K))

eligible_ids <- df_status %>%
    group_by(participant_id) %>%
    summarise(n_tp = n_distinct(Timepoint)) %>%
    filter(n_tp >= K) %>%
    pull(participant_id)

df_cohort <- df_status %>%
    filter(participant_id %in% eligible_ids)

cat("Clinical Cohort:\n")
cat(sprintf("Participants: %d\n", n_distinct(df_cohort$participant_id)))
cat(sprintf("Episodes: %d\n", nrow(df_cohort)))

# Join MLST
df_merged <- inner_join(df_cohort, df_mlst_clean, by = c("participant_id", "Timepoint"))

cat("\nMerged (Typed) Cohort:\n")
n_merged_part <- n_distinct(df_merged$participant_id)
n_merged_ep <- nrow(df_merged)
cat(sprintf("Participants with >=1 Typed Isolate: %d\n", n_merged_part))
cat(sprintf("Total Typed Episodes: %d\n", n_merged_ep))

# Pairs Analysis
stable_check <- df_merged %>%
    group_by(participant_id) %>%
    mutate(n_typed = n()) %>%
    filter(n_typed >= 2) %>%
    arrange(Timepoint) # Use Timepoint or derived num

n_participants_with_pairs <- n_distinct(stable_check$participant_id)
n_episodes_in_pairs_cohort <- nrow(stable_check)

# Manual Pair Count
# Pairs = Sum(n_i - 1) for each participant
pair_counts <- stable_check %>%
    group_by(participant_id) %>%
    summarise(n = n(), pairs = n - 1)

total_pairs <- sum(pair_counts$pairs)

cat("\nPair Analysis:\n")
cat(sprintf("Participants contributing pairs (>=2 typed): %d\n", n_participants_with_pairs))
cat(sprintf("Episodes in these participants: %d\n", n_episodes_in_pairs_cohort))
cat(sprintf("Calculated Pairs (Sum n_i - 1): %d\n", total_pairs))

# Check difference
cat(sprintf("Theoretical Max Pairs (Total Typed - Total Part): %d - %d = %d\n", n_merged_ep, n_merged_part, n_merged_ep - n_merged_part))

# K=3 Analysis
K <- 3
cat(sprintf("\n--- K=%d Analysis ---\n", K))

eligible_ids_3 <- df_status %>%
    group_by(participant_id) %>%
    summarise(n_tp = n_distinct(Timepoint)) %>%
    filter(n_tp >= K) %>%
    pull(participant_id)

df_cohort_3 <- df_status %>%
    filter(participant_id %in% eligible_ids_3)

df_merged_3 <- inner_join(df_cohort_3, df_mlst_clean, by = c("participant_id", "Timepoint"))

cat(sprintf("Participants (Typed): %d\n", n_distinct(df_merged_3$participant_id)))
cat(sprintf("Episodes (Typed): %d\n", nrow(df_merged_3)))

stable_check_3 <- df_merged_3 %>%
    group_by(participant_id) %>%
    filter(n() >= 2)

pair_counts_3 <- stable_check_3 %>%
    group_by(participant_id) %>%
    summarise(n = n(), pairs = n - 1)

total_pairs_3 <- sum(pair_counts_3$pairs)
cat(sprintf("Calculated Pairs: %d\n", total_pairs_3))

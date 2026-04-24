#!/usr/bin/env Rscript

# identify_st_minus_individuals.R
# Identify which participants have the 31 episodes with "ST-" (untypeable)

suppressPackageStartupMessages({
    library(tidyverse)
    library(knitr)
})

# Load data
cat("Loading data...\n")
df_status <- read_csv("results/clinical/status_map.csv", show_col_types = FALSE) %>%
    mutate(Participant_id = as.character(Participant_id))

df_mlst <- read_csv("results/mlst/mlst_with_meta.csv", show_col_types = FALSE) %>%
    mutate(Participant_id = as.character(Participant_id))

# Strategy: Deduplicate MLST per episode (Participant_id + Timepoint)
# Prefer complete MLST (no new allele, no ambiguity) -> higher n_loci -> first row
ALLELE_COLS <- c("dinb", "icda", "pabb", "polb", "putp", "trpa", "trpb", "uida")
actual_allele_cols <- intersect(names(df_mlst), ALLELE_COLS)

df_mlst_clean <- df_mlst %>%
    mutate(
        ST_clean = suppressWarnings(as.numeric(as.character(ST))),
        has_new = replace_na(has_new_allele, FALSE),
        is_ambig = replace_na(ambiguous_call, FALSE),
        mlst_complete = !has_new & !is_ambig & !is.na(ST_clean),
        n_loci_typed = rowSums(!is.na(select(., all_of(actual_allele_cols))))
    ) %>%
    group_by(Participant_id, Timepoint) %>%
    arrange(desc(mlst_complete), desc(n_loci_typed), Isolate_ID) %>%
    mutate(rank = row_number()) %>%
    ungroup() %>%
    filter(rank == 1)

# Join MLST to status_map
df_linked <- df_status %>%
    left_join(df_mlst_clean %>% select(Participant_id, Timepoint, ST, ST_clean, Isolate_ID, has_new, is_ambig),
        by = c("Participant_id", "Timepoint")
    )

# Filter to ASB/UTI episodes only (as per the abstract analysis)
df_asb_uti <- df_linked %>%
    filter(Infection_Status %in% c("ASB", "UTI"))

cat("\n=== TOTAL EPISODES ===\n")
cat(glue::glue("Total ASB/UTI episodes: {nrow(df_asb_uti)}\n"))
cat(glue::glue("Episodes with ST data: {sum(!is.na(df_asb_uti$ST))}\n"))
cat(glue::glue("Episodes with missing ST: {sum(is.na(df_asb_uti$ST))}\n\n"))

# Identify episodes with ST = "-" or equivalent (untypeable)
df_st_minus <- df_asb_uti %>%
    filter(ST == "-" | is.na(ST_clean))

cat("=== UNTYPEABLE ST EPISODES ===\n")
cat(glue::glue("Total episodes with untypeable ST (ST = '-' or no numeric ST): {nrow(df_st_minus)}\n\n"))

# Identify unique participants
unique_participants <- df_st_minus %>%
    distinct(Participant_id) %>%
    pull(Participant_id) %>%
    sort()

cat("=== UNIQUE PARTICIPANTS WITH ST- EPISODES ===\n")
cat(glue::glue("Number of unique participants: {length(unique_participants)}\n"))
cat("Participant IDs:\n")
cat(paste(unique_participants, collapse = ", "), "\n\n")

# Summary by participant: how many ST- episodes does each have?
participant_summary <- df_st_minus %>%
    group_by(Participant_id) %>%
    summarise(
        n_st_minus_episodes = n(),
        timepoints = paste(Timepoint, collapse = ", "),
        infection_statuses = paste(Infection_Status, collapse = ", "),
        .groups = "drop"
    ) %>%
    arrange(desc(n_st_minus_episodes))

cat("=== BREAKDOWN BY PARTICIPANT ===\n")
print(kable(participant_summary))

# Detailed list of all ST- episodes
cat("\n=== DETAILED LIST OF ALL ST- EPISODES ===\n")
detailed_list <- df_st_minus %>%
    select(Participant_id, Timepoint, Infection_Status, ST, Isolate_ID, has_new, is_ambig) %>%
    arrange(Participant_id, Timepoint)

print(kable(detailed_list))

# Check how many participants have >=2 timepoints (the primary cohort from abstract)
ppt_counts <- df_status %>%
    group_by(Participant_id) %>%
    summarise(n_timepoints = n_distinct(Timepoint))

eligible_ids_k2 <- ppt_counts %>%
    filter(n_timepoints >= 2) %>%
    pull(Participant_id)

# Of the ST- participants, how many are in the >=2 timepoints cohort?
st_minus_in_k2 <- unique_participants[unique_participants %in% eligible_ids_k2]

cat("\n=== CONTEXT: PRIMARY COHORT (>=2 TIMEPOINTS) ===\n")
cat(glue::glue("Total participants with >=2 timepoints: {length(eligible_ids_k2)}\n"))
cat(glue::glue("ST- participants in >=2 timepoints cohort: {length(st_minus_in_k2)}\n"))
cat("These participants:\n")
cat(paste(st_minus_in_k2, collapse = ", "), "\n\n")

# Save output to a text file
output_file <- "ST_minus_individuals_report.txt"
sink(output_file)

cat("========================================\n")
cat("REPORT: Individuals with ST- Episodes\n")
cat("========================================\n")
cat(glue::glue("Generated: {Sys.time()}\n\n"))

cat("=== SUMMARY ===\n")
cat(glue::glue("Total ASB/UTI episodes analyzed: {nrow(df_asb_uti)}\n"))
cat(glue::glue("Episodes with untypeable ST (ST-): {nrow(df_st_minus)}\n"))
cat(glue::glue("Percentage untypeable: {round(nrow(df_st_minus)/nrow(df_asb_uti)*100, 1)}%\n\n"))

cat(glue::glue("Unique participants with ST- episodes: {length(unique_participants)}\n"))
cat("Participant IDs with ST- episodes:\n")
cat(paste(unique_participants, collapse = ", "), "\n\n")

cat("=== BREAKDOWN BY PARTICIPANT ===\n")
print(kable(participant_summary))

cat("\n=== DETAILED EPISODE LIST ===\n")
print(kable(detailed_list))

cat("\n=== COHORT CONTEXT ===\n")
cat(glue::glue("Participants in >=2 timepoints cohort (primary analysis): {length(eligible_ids_k2)}\n"))
cat(glue::glue("ST- participants in primary cohort: {length(st_minus_in_k2)}\n"))

sink()

cat(glue::glue("\n✓ Report saved to: {output_file}\n"))
cat("Done!\n")

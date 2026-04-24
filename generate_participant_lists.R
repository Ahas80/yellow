#!/usr/bin/env Rscript

# generate_participant_lists.R
# Create formatted lists of participants with ST- and no MLST data

suppressPackageStartupMessages({
    library(tidyverse)
    library(knitr)
})

# Load data
df_status <- read_csv("results/clinical/status_map.csv", show_col_types = FALSE) %>%
    mutate(Participant_id = as.character(Participant_id))

df_mlst <- read_csv("results/mlst/mlst_with_meta.csv", show_col_types = FALSE) %>%
    mutate(Participant_id = as.character(Participant_id))

# Deduplicate MLST
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

# Join
df_linked <- df_status %>%
    left_join(df_mlst_clean %>% select(Participant_id, Timepoint, ST, ST_clean, Isolate_ID),
        by = c("Participant_id", "Timepoint")
    )

# Filter to ASB/UTI
df_asb_uti <- df_linked %>%
    filter(Infection_Status %in% c("ASB", "UTI"))

# === LIST 1: Participants with ST = "-" ===
cat("==============================================\n")
cat("LIST 1: Participants with ST = '-' Episodes\n")
cat("==============================================\n\n")

st_minus_episodes <- df_asb_uti %>%
    filter(ST == "-" & !is.na(ST))

st_minus_participants <- st_minus_episodes %>%
    group_by(Participant_id) %>%
    summarise(
        n_episodes = n(),
        timepoints = paste(sort(Timepoint), collapse = ", "),
        statuses = paste(Infection_Status, collapse = ", "),
        .groups = "drop"
    ) %>%
    arrange(desc(n_episodes), Participant_id)

cat(glue::glue("Total participants: {nrow(st_minus_participants)}\n"))
cat(glue::glue("Total episodes: {nrow(st_minus_episodes)}\n\n"))

print(kable(st_minus_participants))

cat("\n\nParticipant IDs (comma-separated):\n")
cat(paste(sort(st_minus_participants$Participant_id), collapse = ", "))
cat("\n\n")

# === LIST 2: Participants with No MLST data ===
cat("==============================================\n")
cat("LIST 2: Participants with No MLST Data\n")
cat("==============================================\n\n")

no_mlst_episodes <- df_asb_uti %>%
    filter(is.na(ST))

no_mlst_participants <- no_mlst_episodes %>%
    group_by(Participant_id) %>%
    summarise(
        n_episodes = n(),
        timepoints = paste(sort(Timepoint), collapse = ", "),
        statuses = paste(Infection_Status, collapse = ", "),
        .groups = "drop"
    ) %>%
    arrange(desc(n_episodes), Participant_id)

cat(glue::glue("Total participants: {nrow(no_mlst_participants)}\n"))
cat(glue::glue("Total episodes: {nrow(no_mlst_episodes)}\n\n"))

print(kable(no_mlst_participants))

cat("\n\nParticipant IDs (comma-separated):\n")
cat(paste(sort(no_mlst_participants$Participant_id), collapse = ", "))
cat("\n\n")

# === COMBINED SUMMARY ===
cat("==============================================\n")
cat("COMBINED SUMMARY\n")
cat("==============================================\n\n")

all_untypeable <- unique(c(st_minus_participants$Participant_id, no_mlst_participants$Participant_id))
only_st_minus <- setdiff(st_minus_participants$Participant_id, no_mlst_participants$Participant_id)
only_no_mlst <- setdiff(no_mlst_participants$Participant_id, st_minus_participants$Participant_id)
both_types <- intersect(st_minus_participants$Participant_id, no_mlst_participants$Participant_id)

cat(glue::glue("Total unique participants with untypeable episodes: {length(all_untypeable)}\n\n"))

cat(glue::glue("Breakdown:\n"))
cat(glue::glue("  - Only ST='-' (no NA episodes): {length(only_st_minus)}\n"))
cat(glue::glue("  - Only no MLST data (no ST='-' episodes): {length(only_no_mlst)}\n"))
cat(glue::glue("  - Both types (have ST='-' AND NA episodes): {length(both_types)}\n\n"))

if (length(both_types) > 0) {
    cat("Participants with BOTH types:\n")
    cat(paste(sort(both_types), collapse = ", "))
    cat("\n\n")
}

# Save to file
sink("participant_lists_formatted.txt")

cat("==============================================\n")
cat("COMPLETE PARTICIPANT LISTS\n")
cat("Generated: ", as.character(Sys.time()), "\n")
cat("==============================================\n\n")

cat("LIST 1: Participants with ST = '-' Episodes\n")
cat("--------------------------------------------\n")
cat(glue::glue("Total: {nrow(st_minus_participants)} participants, {nrow(st_minus_episodes)} episodes\n\n"))
print(kable(st_minus_participants))
cat("\n\nIDs: ", paste(sort(st_minus_participants$Participant_id), collapse = ", "), "\n\n\n")

cat("LIST 2: Participants with No MLST Data\n")
cat("---------------------------------------\n")
cat(glue::glue("Total: {nrow(no_mlst_participants)} participants, {nrow(no_mlst_episodes)} episodes\n\n"))
print(kable(no_mlst_participants))
cat("\n\nIDs: ", paste(sort(no_mlst_participants$Participant_id), collapse = ", "), "\n\n\n")

cat("COMBINED SUMMARY\n")
cat("----------------\n")
cat(glue::glue("Total unique: {length(all_untypeable)}\n"))
cat(glue::glue("Only ST='-': {length(only_st_minus)}\n"))
cat(glue::glue("Only no data: {length(only_no_mlst)}\n"))
cat(glue::glue("Both types: {length(both_types)}\n"))
if (length(both_types) > 0) {
    cat("\nBoth types: ", paste(sort(both_types), collapse = ", "), "\n")
}

sink()

cat("✓ Report saved to: participant_lists_formatted.txt\n")
cat("Done!\n")

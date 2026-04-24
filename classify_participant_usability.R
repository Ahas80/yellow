#!/usr/bin/env Rscript

# classify_participant_usability.R
# Determine which participants are unusable for ST analysis and why

suppressPackageStartupMessages({
    library(tidyverse)
    library(knitr)
    library(glue)
})

cat("===================================================\n")
cat("PARTICIPANT USABILITY CLASSIFICATION\n")
cat("===================================================\n\n")

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

# Filter to ASB/UTI only
df_asb_uti <- df_linked %>%
    filter(Infection_Status %in% c("ASB", "UTI"))

# Classify each episode
df_classified <- df_asb_uti %>%
    mutate(
        episode_status = case_when(
            is.na(ST) ~ "No data",
            ST == "-" ~ "Ambiguous",
            !is.na(ST_clean) ~ "Usable",
            TRUE ~ "Unknown"
        )
    )

# Summarize by participant
participant_summary <- df_classified %>%
    group_by(Participant_id) %>%
    summarise(
        n_total = n(),
        n_usable = sum(episode_status == "Usable"),
        n_ambiguous = sum(episode_status == "Ambiguous"),
        n_no_data = sum(episode_status == "No data"),
        pct_usable = round(n_usable / n_total * 100, 1),
        timepoints = paste(sort(Timepoint), collapse = ", "),
        .groups = "drop"
    ) %>%
    mutate(
        usability_category = case_when(
            n_usable == 0 & n_total > 0 ~ "COMPLETELY UNUSABLE",
            n_usable == n_total ~ "FULLY USABLE",
            n_usable > 0 & n_usable < n_total ~ "PARTIALLY USABLE",
            TRUE ~ "CHECK"
        ),
        unusable_type = case_when(
            n_usable == 0 & n_ambiguous > 0 & n_no_data == 0 ~ "All ambiguous",
            n_usable == 0 & n_no_data > 0 & n_ambiguous == 0 ~ "All no data",
            n_usable == 0 & n_ambiguous > 0 & n_no_data > 0 ~ "Mixed (ambiguous + no data)",
            n_usable > 0 & (n_ambiguous > 0 | n_no_data > 0) ~ "Partially missing",
            n_usable == n_total ~ "None - all usable",
            TRUE ~ "N/A"
        )
    ) %>%
    arrange(usability_category, desc(n_total), Participant_id)

# === CATEGORY 1: Completely unusable ===
completely_unusable <- participant_summary %>%
    filter(usability_category == "COMPLETELY UNUSABLE")

cat("CATEGORY 1: COMPLETELY UNUSABLE\n")
cat("================================\n")
cat(glue("Total: {nrow(completely_unusable)} participants\n"))
cat("These participants have NO usable ST data for ANY episode\n\n")

if (nrow(completely_unusable) > 0) {
    cat("Breakdown by type:\n")
    type_summary <- completely_unusable %>%
        count(unusable_type, name = "n_participants") %>%
        arrange(desc(n_participants))
    print(kable(type_summary))

    cat("\n\nComplete list:\n")
    print(kable(completely_unusable %>%
        select(Participant_id, n_total, n_ambiguous, n_no_data, unusable_type, timepoints)))

    cat("\n\nParticipant IDs (comma-separated):\n")
    cat(paste(sort(completely_unusable$Participant_id), collapse = ", "))
    cat("\n\n")
}

# === CATEGORY 2: Partially unusable ===
partially_unusable <- participant_summary %>%
    filter(usability_category == "PARTIALLY USABLE")

cat("\nCATEGORY 2: PARTIALLY USABLE\n")
cat("============================\n")
cat(glue("Total: {nrow(partially_unusable)} participants\n"))
cat("These participants have SOME usable ST data, but not all episodes\n\n")

if (nrow(partially_unusable) > 0) {
    cat("Summary by usability percentage:\n")
    pct_bins <- partially_unusable %>%
        mutate(pct_bin = cut(pct_usable,
            breaks = c(0, 25, 50, 75, 100),
            labels = c("1-25%", ">25-50%", ">50-75%", ">75-99%")
        )) %>%
        count(pct_bin, name = "n_participants")
    print(kable(pct_bins))

    cat("\n\nMost affected (lowest % usable):\n")
    print(kable(head(partially_unusable %>%
        select(Participant_id, n_total, n_usable, n_ambiguous, n_no_data, pct_usable, timepoints), 10)))

    cat("\n\nAll partially usable participant IDs:\n")
    cat(paste(sort(partially_unusable$Participant_id), collapse = ", "))
    cat("\n\n")
}

# === CATEGORY 3: Fully usable ===
fully_usable <- participant_summary %>%
    filter(usability_category == "FULLY USABLE")

cat("\nCATEGORY 3: FULLY USABLE\n")
cat("========================\n")
cat(glue("Total: {nrow(fully_usable)} participants\n"))
cat("These participants have 100% usable ST data for all episodes\n\n")

cat(glue("(Not listing all {nrow(fully_usable)} - these are good for analysis)\n\n"))

# === IMPACT ANALYSIS ===
cat("\nIMPACT ANALYSIS\n")
cat("===============\n\n")

total_ppts <- nrow(participant_summary)
total_episodes <- sum(participant_summary$n_total)
usable_episodes <- sum(participant_summary$n_usable)
unusable_episodes <- total_episodes - usable_episodes

cat("Overall cohort impact:\n")
cat(glue("  Total participants (ASB/UTI): {total_ppts}\n"))
cat(glue("  Total episodes: {total_episodes}\n"))
cat(glue("  Usable episodes: {usable_episodes} ({round(usable_episodes/total_episodes*100,1)}%)\n"))
cat(glue("  Unusable episodes: {unusable_episodes} ({round(unusable_episodes/total_episodes*100,1)}%)\n\n"))

cat("Participant-level impact:\n")
cat(glue("  Completely unusable: {nrow(completely_unusable)} ({round(nrow(completely_unusable)/total_ppts*100,1)}%)\n"))
cat(glue("  Partially usable: {nrow(partially_unusable)} ({round(nrow(partially_unusable)/total_ppts*100,1)}%)\n"))
cat(glue("  Fully usable: {nrow(fully_usable)} ({round(nrow(fully_usable)/total_ppts*100,1)}%)\n\n"))

# Episode-level detail by participant category
cat("Episodes by participant category:\n")
episodes_by_cat <- participant_summary %>%
    group_by(usability_category) %>%
    summarise(
        n_participants = n(),
        total_episodes = sum(n_total),
        usable_episodes = sum(n_usable),
        pct_usable = round(usable_episodes / total_episodes * 100, 1)
    )
print(kable(episodes_by_cat))

cat("\n")

# === RECOMMENDATIONS ===
cat("\nRECOMMENDATIONS\n")
cat("===============\n\n")

cat("1. EXCLUDE FROM ST ANALYSIS:\n")
cat(glue("   {nrow(completely_unusable)} participants with 0% usable data\n"))
cat(glue("   Impact: Lose {sum(completely_unusable$n_total)} episodes\n\n"))

cat("2. USE WITH CAUTION:\n")
cat(glue("   {nrow(partially_unusable)} participants with partial data\n"))
cat("   Depending on analysis:\n")
cat("   - Longitudinal ST stability: May be usable if >=2 typed episodes\n")
cat("   - ST distribution: Can use available episodes\n")
cat("   - Participant-level stats: May introduce bias\n\n")

cat("3. FULLY USABLE:\n")
cat(glue("   {nrow(fully_usable)} participants (100% typed)\n"))
cat("   Safe for all ST analyses\n\n")

# Save report
sink("participant_usability_classification.txt")

cat("===================================================\n")
cat("PARTICIPANT USABILITY CLASSIFICATION\n")
cat(glue("Generated: {Sys.time()}\n"))
cat("===================================================\n\n")

cat("CATEGORY 1: COMPLETELY UNUSABLE\n")
cat(glue("{nrow(completely_unusable)} participants with NO usable ST data\n\n"))
if (nrow(completely_unusable) > 0) {
    print(kable(completely_unusable %>%
        select(Participant_id, n_total, n_ambiguous, n_no_data, unusable_type, timepoints)))
    cat("\nIDs: ", paste(sort(completely_unusable$Participant_id), collapse = ", "), "\n\n")
}

cat("\nCATEGORY 2: PARTIALLY USABLE\n")
cat(glue("{nrow(partially_unusable)} participants with SOME usable ST data\n\n"))
if (nrow(partially_unusable) > 0) {
    print(kable(partially_unusable %>%
        select(Participant_id, n_total, n_usable, pct_usable, unusable_type, timepoints)))
    cat("\nIDs: ", paste(sort(partially_unusable$Participant_id), collapse = ", "), "\n\n")
}

cat("\nCATEGORY 3: FULLY USABLE\n")
cat(glue("{nrow(fully_usable)} participants with 100% usable ST data\n"))
cat("(All episodes successfully typed)\n\n")

cat("\nIMPACT SUMMARY:\n")
print(kable(episodes_by_cat))

sink()

cat("✓ Report saved to: participant_usability_classification.txt\n")
cat("Done!\n")

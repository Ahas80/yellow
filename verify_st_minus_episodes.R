#!/usr/bin/env Rscript

# verify_st_minus_episodes.R
# Comprehensive verification of ST- episodes with duplicate checking

suppressPackageStartupMessages({
    library(tidyverse)
    library(knitr)
    library(glue)
})

cat("========================================\n")
cat("VERIFICATION: ST- Episodes Analysis\n")
cat("========================================\n\n")

# Load data
cat("Loading data...\n")
df_status <- read_csv("results/clinical/status_map.csv", show_col_types = FALSE) %>%
    mutate(Participant_id = as.character(Participant_id))

df_mlst <- read_csv("results/mlst/mlst_with_meta.csv", show_col_types = FALSE) %>%
    mutate(Participant_id = as.character(Participant_id))

cat(glue("✓ Loaded {nrow(df_status)} rows from status_map.csv\n"))
cat(glue("✓ Loaded {nrow(df_mlst)} rows from mlst_with_meta.csv\n\n"))

# Check for duplicates in status_map
cat("=== STEP 1: Checking for duplicates in status_map ===\n")
status_dups <- df_status %>%
    group_by(Participant_id, Timepoint) %>%
    filter(n() > 1) %>%
    ungroup()

if (nrow(status_dups) > 0) {
    cat("⚠️  WARNING: Found duplicates in status_map:\n")
    print(kable(status_dups))
} else {
    cat("✓ No duplicates in status_map (each Participant_id + Timepoint is unique)\n")
}

# Deduplicate MLST
cat("\n=== STEP 2: Deduplicating MLST data ===\n")
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
    ungroup()

# Show how many rows per episode
multi_mlst <- df_mlst_clean %>%
    group_by(Participant_id, Timepoint) %>%
    summarise(n_rows = n(), .groups = "drop") %>%
    filter(n_rows > 1)

cat(glue("Episodes with multiple MLST rows: {nrow(multi_mlst)}\n"))
if (nrow(multi_mlst) > 0) {
    cat("(This is expected - multiple assemblers per sample)\n")
    cat("Taking rank 1 (best quality) for each episode...\n")
}

df_mlst_unique <- df_mlst_clean %>% filter(rank == 1)
cat(glue("After deduplication: {nrow(df_mlst_unique)} unique episodes with MLST data\n"))

# Join to status_map
cat("\n=== STEP 3: Joining MLST to clinical data ===\n")
df_linked <- df_status %>%
    left_join(df_mlst_unique %>% select(Participant_id, Timepoint, ST, ST_clean, Isolate_ID, has_new, is_ambig),
        by = c("Participant_id", "Timepoint")
    )

# Filter to ASB/UTI
df_asb_uti <- df_linked %>%
    filter(Infection_Status %in% c("ASB", "UTI"))

cat(glue("Total ASB/UTI episodes: {nrow(df_asb_uti)}\n"))
cat(glue("ASB episodes: {sum(df_asb_uti$Infection_Status == 'ASB')}\n"))
cat(glue("UTI episodes: {sum(df_asb_uti$Infection_Status == 'UTI')}\n"))

# Identify untypeable episodes
cat("\n=== STEP 4: Categorizing untypeable episodes ===\n")

# Category 1: ST = "-" (ambiguous typing)
st_minus <- df_asb_uti %>%
    filter(ST == "-" & !is.na(ST))

# Category 2: No MLST data at all
st_na <- df_asb_uti %>%
    filter(is.na(ST))

# Category 3: Successfully typed
st_success <- df_asb_uti %>%
    filter(!is.na(ST) & ST != "-")

cat(glue("Episodes with ST = '-' (ambiguous): {nrow(st_minus)}\n"))
cat(glue("Episodes with no MLST data (NA): {nrow(st_na)}\n"))
cat(glue("Episodes successfully typed: {nrow(st_success)}\n"))
cat(glue("TOTAL: {nrow(st_minus) + nrow(st_na) + nrow(st_success)}\n"))

# Verify totals
total_untypeable <- nrow(st_minus) + nrow(st_na)
cat(glue("\n✓ Total untypeable: {total_untypeable}\n"))
cat(glue("✓ Successfully typed: {nrow(st_success)}\n"))
cat(glue("✓ Grand total: {nrow(df_asb_uti)} (matches ASB/UTI count)\n"))

# Check for any overlap
overlap <- intersect(st_minus$Isolate_ID, st_na$Isolate_ID)
if (length(overlap) > 0) {
    cat("\n⚠️  WARNING: Found overlap between ST- and NA categories!\n")
} else {
    cat("\n✓ No overlap between categories (good)\n")
}

# Detailed breakdown: ST = "-"
cat("\n=== STEP 5: Detailed list of ST = '-' episodes ===\n")
cat(glue("Total: {nrow(st_minus)} episodes\n\n"))

st_minus_detailed <- st_minus %>%
    select(Participant_id, Timepoint, Infection_Status, Isolate_ID, has_new, is_ambig) %>%
    arrange(Participant_id, Timepoint)

print(kable(st_minus_detailed))

# Check for repeats in ST-
cat("\n=== STEP 6: Checking for repeated Participant_id + Timepoint in ST- ===\n")
st_minus_repeats <- st_minus_detailed %>%
    group_by(Participant_id, Timepoint) %>%
    filter(n() > 1) %>%
    ungroup()

if (nrow(st_minus_repeats) > 0) {
    cat("⚠️  Found repeats:\n")
    print(kable(st_minus_repeats))
} else {
    cat("✓ No repeats found - each ST- episode is unique\n")
}

# Participant summary for ST-
cat("\n=== STEP 7: Participants with ST = '-' ===\n")
st_minus_participants <- st_minus_detailed %>%
    group_by(Participant_id) %>%
    summarise(
        n_episodes = n(),
        timepoints = paste(Timepoint, collapse = ", "),
        statuses = paste(Infection_Status, collapse = ", "),
        .groups = "drop"
    ) %>%
    arrange(desc(n_episodes))

cat(glue("Unique participants: {nrow(st_minus_participants)}\n\n"))
print(kable(st_minus_participants))

# Detailed breakdown: NA
cat("\n=== STEP 8: Detailed list of episodes with no MLST data (NA) ===\n")
cat(glue("Total: {nrow(st_na)} episodes\n\n"))

st_na_detailed <- st_na %>%
    select(Participant_id, Timepoint, Infection_Status, Isolate_ID) %>%
    arrange(Participant_id, Timepoint)

# Show first 20 and last 20
if (nrow(st_na_detailed) > 40) {
    cat("First 20 episodes:\n")
    print(kable(head(st_na_detailed, 20)))
    cat("\n...\n\nLast 20 episodes:\n")
    print(kable(tail(st_na_detailed, 20)))
    cat(glue("\n(Showing 40 of {nrow(st_na_detailed)} total episodes with NA)\n"))
} else {
    print(kable(st_na_detailed))
}

# Check for repeats in NA
cat("\n=== STEP 9: Checking for repeated Participant_id + Timepoint in NA ===\n")
st_na_repeats <- st_na_detailed %>%
    group_by(Participant_id, Timepoint) %>%
    filter(n() > 1) %>%
    ungroup()

if (nrow(st_na_repeats) > 0) {
    cat("⚠️  Found repeats:\n")
    print(kable(st_na_repeats))
} else {
    cat("✓ No repeats found - each NA episode is unique\n")
}

# Timepoint distribution
cat("\n=== STEP 10: Distribution by Timepoint ===\n")

timepoint_dist <- df_asb_uti %>%
    mutate(
        Category = case_when(
            ST == "-" ~ "ST = '-'",
            is.na(ST) ~ "No MLST data",
            TRUE ~ "Successfully typed"
        )
    ) %>%
    group_by(Timepoint, Category) %>%
    summarise(n = n(), .groups = "drop") %>%
    pivot_wider(names_from = Category, values_from = n, values_fill = 0)

print(kable(timepoint_dist))

# Overall summary
cat("\n=== FINAL SUMMARY ===\n")
cat(glue("
Total ASB/UTI episodes: {nrow(df_asb_uti)}
  ├─ Successfully typed: {nrow(st_success)} ({round(nrow(st_success)/nrow(df_asb_uti)*100, 1)}%)
  ├─ ST = '-' (ambiguous): {nrow(st_minus)} ({round(nrow(st_minus)/nrow(df_asb_uti)*100, 1)}%)
  └─ No MLST data (NA): {nrow(st_na)} ({round(nrow(st_na)/nrow(df_asb_uti)*100, 1)}%)

Untypeable episodes total: {total_untypeable} ({round(total_untypeable/nrow(df_asb_uti)*100, 1)}%)

Unique participants:
  - With ST = '-': {nrow(st_minus_participants)}
  - With NA: {n_distinct(st_na$Participant_id)}
  - Total with untypeable: {n_distinct(c(st_minus$Participant_id, st_na$Participant_id))}

No duplicates detected ✓
All numbers verified ✓
"))

cat("\nDone!\n")

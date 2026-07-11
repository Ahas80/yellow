#!/usr/bin/env Rscript
# ==============================================================================
# 00a_load_clean_clinical.R
# ==============================================================================
#
# GOAL:
#   Load and merge clinical data batches into a single harmonised dataset.
#   This is the very first step in the pipeline — all subsequent clinical
#   classification depends on this merge.
#
# WHY THIS SCRIPT EXISTS:
#   The YELLOW RoUTIne study collected clinical data in recruitment batches.
#   Each batch has slightly different column naming conventions and data
#   formatting.  This script reconciles those differences so that downstream
#   scripts can work with a consistent dataset.
#
#   Signs & Symptoms (S&S) columns are deliberately separated from metadata
#   because they require special boolean parsing (Dutch-language values like
#   "ja"/"nee", numeric codes, etc.).
#
# BATCH HANDLING (v2 — Batches 1–6):
#   - The list of clinical batch CSVs is defined centrally in 00_config.R
#     (BATCH_DEFINITIONS / BATCH_CLINICAL_CSVS).
#   - Batches 4–6 may not yet have clinical CSV files.  The script loads
#     only those CSVs that exist on disk and warns (not crashes) about missing ones.
#   - Batch 4–6 clinical CSVs are exported from the master overview workbook
#     into data/inputs/ and then listed in BATCH_DEFINITIONS in 00_config.R.
#
# INPUTS:
#   - data/inputs/batch*.csv  (clinical exports — as many as exist)
#
# OUTPUTS:
#   - results/clinical/intermediate/clinical_merged.rds
#     (contains metadata, symptom data, and S&S table as a named list)
#
# DOWNSTREAM:
#   → 00b_classify_episodes.R reads clinical_merged.rds to assign
#     primary UTI_Status (UTI/Not_UTI) plus legacy Infection_Status.
# ==============================================================================

source("00_config.R")
source(here::here("R", "clinical_helpers.R"))

suppressPackageStartupMessages({
    library(dplyr)
    library(readr)
    library(stringr)
    library(tidyr)
})

msg("Starting 00a_load_clean_clinical.R")

# Ensure intermediate directory exists
DIR_INTERMEDIATE <- file.path(DIR_CLINICAL, "intermediate")
ensure_dir(DIR_INTERMEDIATE)

# ==============================================================================
# 1. Load Batches — dynamically from BATCH_DEFINITIONS
# ==============================================================================
# read_batch: reads a single batch CSV, trying data/inputs/ first, then root.
# Returns NULL (not a crash) if the file does not exist.
read_batch <- function(fname) {
    if (is.na(fname)) return(NULL)
    p1 <- file.path("data", "inputs", fname)
    if (file.exists(p1)) {
        return(read_csv(p1, show_col_types = FALSE))
    }
    if (file.exists(fname)) {
        return(read_csv(fname, show_col_types = FALSE))
    }
    # File does not exist — return NULL rather than crashing.
    NULL
}

# Load all defined batches.  Keep track of which ones were found.
batch_data <- list()
batch_ids_loaded <- integer(0)
batch_ids_missing <- integer(0)

for (bdef in BATCH_DEFINITIONS) {
    csv_name <- bdef$clinical_csv
    bid <- bdef$batch_id
    df <- read_batch(csv_name)
    if (!is.null(df)) {
        batch_data[[as.character(bid)]] <- df
        batch_ids_loaded <- c(batch_ids_loaded, bid)
    } else {
        batch_ids_missing <- c(batch_ids_missing, bid)
    }
}

if (length(batch_ids_loaded) == 0) {
    stop("No clinical batch CSV files found. Cannot proceed without clinical data.")
}

msg("Loaded %d clinical batch CSV(s): %s",
    length(batch_ids_loaded), paste(batch_ids_loaded, collapse = ", "))
if (length(batch_ids_missing) > 0) {
    msg("⚠ Missing clinical CSV(s) for batch(es): %s",
        paste(batch_ids_missing, collapse = ", "))
}

# ==============================================================================
# 2. Clean and Merge Metadata (excluding S&S columns)
# ==============================================================================
drop_SnS_cols <- function(df) {
    df %>%
        select(
            -matches(
                paste(
                    "S&S|sns.*status|signs.*symptoms.*status|sympt|dysur|urg|urge|fever|koorts|pijn|pain|burn",
                    "|(^|[_ -])sx[_ -]?present($|[_ -])",
                    sep = ""
                ),
                ignore.case = TRUE
            )
        )
}

# Clean each batch and harmonize column types
cleaned_list <- lapply(names(batch_data), function(bid_str) {
    df <- batch_data[[bid_str]]
    df_clean <- drop_SnS_cols(df)
    # Coerce all columns to character for safe binding
    df_clean <- df_clean %>% mutate(across(everything(), as.character))
    df_clean$Batch <- as.integer(bid_str)
    df_clean
})

# Find common columns and bind
all_cols <- lapply(cleaned_list, names)
common_cols <- Reduce(intersect, all_cols)
msg("Common columns across %d batches: %d", length(cleaned_list), length(common_cols))

metadata <- bind_rows(cleaned_list)

# Standardize ID & Timepoint
standardize_id_tp <- function(df) {
    nm <- names(df)
    id_col <- nm[grepl("^participant[_ ]?id$", nm, ignore.case = TRUE)][1]
    tp_col <- nm[grepl("^time\\s*point$|^timepoint$", nm, ignore.case = TRUE)][1]
    if (is.na(id_col) || is.na(tp_col)) stop("Could not find Participant_id / Timepoint columns in metadata.")
    df %>%
        dplyr::rename(Participant_id = all_of(id_col), Timepoint = all_of(tp_col)) %>%
        dplyr::mutate(
            Participant_id = as.character(Participant_id),
            Timepoint      = canon_tp(Timepoint)
        )
}

metadata <- metadata %>% standardize_id_tp()

# [STAT] Validate Participant_ID uniqueness across batches
# Check if same Participant_id appears in multiple batches with conflicting data
batch_conflicts <- metadata %>%
    group_by(Participant_id) %>%
    summarise(
        n_batches = n_distinct(Batch),
        batches = paste(sort(unique(Batch)), collapse = ","),
        .groups = "drop"
    ) %>%
    filter(n_batches > 1)

if (nrow(batch_conflicts) > 0) {
    msg("⚠ Warning: %d participants appear in multiple batches:", nrow(batch_conflicts))
    print(batch_conflicts)
    msg("This is expected if the same individuals were sampled across batches.")
    msg("If batches represent DIFFERENT cohorts, this is a critical error.")
} else {
    msg("✓ Participant IDs are unique to batches (no cross-batch participants)")
}

msg("Merged metadata: %d rows from %d batches", nrow(metadata), length(batch_ids_loaded))

# ==============================================================================
# 3. Recover S&S Columns (for detailed symptom logic)
# ==============================================================================
symptom_cols2 <- function(df) {
    nm <- names(df)
    id_col <- nm[grepl("^participant[_ ]?id$", nm, ignore.case = TRUE)][1]
    tp_col <- nm[grepl("^time\\s*point$|^timepoint$", nm, ignore.case = TRUE)][1]
    if (is.na(id_col) || is.na(tp_col)) {
        return(NULL)
    }

    df_std <- df %>%
        rename(Participant_id = all_of(id_col), Timepoint = all_of(tp_col)) %>%
        mutate(
            Participant_id = as.character(Participant_id),
            Timepoint = canon_tp(Timepoint)
        )

    detected <- detect_symptom_columns(df_std)
    symptom_cols <- unique(detected$raw_column)
    status_like <- names(df_std)[grepl(
        "status|sns.*status|signs.*symptoms.*status|^no\\s+s\\s*&\\s*s$|^geen\\s+s\\s*&\\s*s$",
        names(df_std),
        ignore.case = TRUE
    )]
    symptom_cols <- setdiff(symptom_cols, status_like)

    if (length(symptom_cols) == 0) {
        return(df_std %>% select(Participant_id, Timepoint))
    }

    df_std %>%
        select(Participant_id, Timepoint, all_of(symptom_cols)) %>%
        mutate(across(-c(Participant_id, Timepoint), as.character))
}

sym_list <- lapply(batch_data, symptom_cols2)
sym_all_raw <- bind_rows(sym_list[!vapply(sym_list, is.null, logical(1))])

detected_symptom_columns <- bind_rows(lapply(names(batch_data), function(bid_str) {
    df <- batch_data[[bid_str]]
    detected <- detect_symptom_columns(df)
    if (nrow(detected) == 0) {
        tibble(
            Batch = as.integer(bid_str),
            raw_column = NA_character_,
            symptom = NA_character_,
            used_for_uti_rule = NA,
            rule_role = "no_usable_symptom_columns_detected"
        )
    } else {
        detected %>% mutate(Batch = as.integer(bid_str), .before = raw_column)
    }
}))
write_csv(detected_symptom_columns, FILE_UTI_DETECTED_COLUMN_MAP)
msg("Wrote detected S&S column map to %s", FILE_UTI_DETECTED_COLUMN_MAP)

# ==============================================================================
# 4. Extract SnS Status (Fallback for B1/B2)
# ==============================================================================
empty_sns_frame <- function() {
    tibble(
        Participant_id = character(), Timepoint = character(),
        SnS_status = integer(), has_SnS = logical(), Batch = integer()
    )
}

extract_sns_status <- function(df, batch_id) {
    nm <- names(df)
    id_col <- nm[grepl("^participant[_ ]?id$", nm, ignore.case = TRUE)][1]
    tp_col <- nm[grepl("^time\\s*point$|^timepoint$", nm, ignore.case = TRUE)][1]
    if (is.na(id_col) || is.na(tp_col)) {
        return(empty_sns_frame())
    }

    df_std <- df %>%
        rename(Participant_id = all_of(id_col), Timepoint = all_of(tp_col)) %>%
        mutate(
            Participant_id = as.character(Participant_id),
            Timepoint = canon_tp(Timepoint)
        )

    status_col <- nm[grepl("s\\s*&\\s*s.*status|sns.*status|signs.*symptoms.*status", nm, ignore.case = TRUE)][1]
    if (is.na(status_col)) {
        return(empty_sns_frame())
    }

    df_std %>%
        transmute(
            Participant_id, Timepoint,
            SnS_status = suppressWarnings(as.integer(readr::parse_number(.data[[status_col]]))),
            has_SnS = !is.na(SnS_status),
            Batch = batch_id
        ) %>%
        filter(!is.na(SnS_status)) %>%
        group_by(Participant_id, Timepoint, Batch) %>%
        summarise(
            SnS_status = max(SnS_status, na.rm = TRUE),
            has_SnS = any(has_SnS, na.rm = TRUE),
            .groups = "drop"
        )
}

# Extract S&S status from batches 1 and 2 (those that have this fallback)
sns_parts <- list()
for (bid_str in names(batch_data)) {
    bid <- as.integer(bid_str)
    if (bid %in% c(1L, 2L)) {
        sns_parts[[bid_str]] <- extract_sns_status(batch_data[[bid_str]], bid)
    }
}
sns_12 <- bind_rows(sns_parts)
if (ncol(sns_12) == 0) sns_12 <- empty_sns_frame()

# ==============================================================================
# 5. Save Intermediate Data
# ==============================================================================
saveRDS(list(
    metadata = metadata,
    sym_all_raw = sym_all_raw,
    sns_12 = sns_12
), file.path(DIR_INTERMEDIATE, "clinical_merged.rds"))

msg("Saved intermediate data to %s", file.path(DIR_INTERMEDIATE, "clinical_merged.rds"))

#!/usr/bin/env Rscript
# ==============================================================================
# 00a_load_clean_clinical.R
# ------------------------------------------------------------------------------
# Role: [Data Prep] - Load and merge raw clinical batches.
#
# Inputs:
#   - data/inputs/batch1.csv (and batch2, batch3...)
#
# Outputs:
#   - results/clinical/intermediate/clinical_merged.rds
#
# Usage:
#   Rscript 00a_load_clean_clinical.R
#
# Biological/Statistical purpose:
#   - Harmonizes clinical data from multiple collection batches to create a
#     single, consistent dataset for downstream analysis.
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

# 1. Load Batches
# ------------------------------------------------------------------------------
read_batch <- function(fname) {
    # Try data/inputs/ first, then root
    p1 <- file.path("data", "inputs", fname)
    if (file.exists(p1)) {
        return(read_csv(p1, show_col_types = FALSE))
    }
    if (file.exists(fname)) {
        return(read_csv(fname, show_col_types = FALSE))
    }
    stop("Could not find input file: ", fname)
}

batch1 <- read_batch("batch1.csv")
batch2 <- read_batch("batch2.csv")
batch3 <- read_batch("batch3.csv")

msg("Loaded batches: B1=%d, B2=%d, B3=%d rows", nrow(batch1), nrow(batch2), nrow(batch3))

# 2. Clean and Merge Metadata (excluding S&S columns)
# ------------------------------------------------------------------------------
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

b1_clean <- drop_SnS_cols(batch1)
b2_clean <- drop_SnS_cols(batch2)
b3_clean <- drop_SnS_cols(batch3)

common_cols <- Reduce(intersect, list(names(b1_clean), names(b2_clean), names(b3_clean)))

metadata <- bind_rows(
    mutate(b1_clean, across(all_of(common_cols), as.character), Batch = 1L),
    mutate(b2_clean, across(all_of(common_cols), as.character), Batch = 2L),
    mutate(b3_clean, across(all_of(common_cols), as.character), Batch = 3L)
)

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

msg("Merged metadata: %d rows", nrow(metadata))

# 3. Recover S&S Columns (for detailed symptom logic)
# ------------------------------------------------------------------------------
symptom_cols2 <- function(df) {
    nm <- names(df)
    id_col <- nm[grepl("^participant[_ ]?id$", nm, ignore.case = TRUE)][1]
    tp_col <- nm[grepl("^time\\s*point$|^timepoint$", nm, ignore.case = TRUE)][1]
    if (is.na(id_col) || is.na(tp_col)) {
        return(NULL)
    }

    df %>%
        rename(Participant_id = all_of(id_col), Timepoint = all_of(tp_col)) %>%
        mutate(
            Participant_id = as.character(Participant_id),
            Timepoint = canon_tp(Timepoint)
        ) %>%
        select(
            Participant_id, Timepoint,
            matches("sympt|dysur|urg|urge|fever|koorts|pijn|pain|burn", ignore.case = TRUE),
            -matches("status|sns.*status|signs.*symptoms.*status", ignore.case = TRUE)
        ) %>%
        mutate(across(-c(Participant_id, Timepoint), as.character))
}

sym_b1 <- symptom_cols2(batch1)
sym_b2 <- symptom_cols2(batch2)
sym_b3 <- symptom_cols2(batch3)
sym_list <- list(sym_b1, sym_b2, sym_b3)
sym_all_raw <- bind_rows(sym_list[!vapply(sym_list, is.null, logical(1))])

# 4. Extract SnS Status (Fallback for B1/B2)
# ------------------------------------------------------------------------------
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

sns_12 <- bind_rows(
    extract_sns_status(batch1, 1L),
    extract_sns_status(batch2, 2L)
)
if (ncol(sns_12) == 0) sns_12 <- empty_sns_frame()

# 5. Save Intermediate Data
# ------------------------------------------------------------------------------
saveRDS(list(
    metadata = metadata,
    sym_all_raw = sym_all_raw,
    sns_12 = sns_12
), file.path(DIR_INTERMEDIATE, "clinical_merged.rds"))

msg("Saved intermediate data to %s", file.path(DIR_INTERMEDIATE, "clinical_merged.rds"))

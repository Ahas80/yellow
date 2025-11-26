#!/usr/bin/env Rscript
# ==============================================================================
# 00_make_assembly_metadata.r
# ------------------------------------------------------------------------------
# GOAL:
#   Generate the master metadata file (assembly_metadata.csv) linking isolates
#   to participants and timepoints.
#
# METHOD:
#   1. Load MLST results (provenance) or scan FASTA directory.
#   2. Normalize Isolate IDs.
#   3. Map to Participant/Timepoint using mapping file or VF hits.
#   4. Write master metadata CSV.
#
# INPUTS:
#   - 00_config.R
#   - results/mlst/mlst_all.tsv (optional)
#   - ont-yellow-routine-fastas/ (raw FASTAs)
#
# OUTPUTS:
#   - assembly_metadata.csv
#   - results/assembly_metadata_TODO_fill.csv (if missing metadata)
# ==============================================================================

# 1. Load Configuration & Libraries
source("00_config.R")
suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
  library(tidyr)
  library(fs)
})

# 2. Helper Functions
# Normalise isolate IDs to ensure consistent joining
norm_iso_id <- function(x) {
  x <- as.character(x)
  x <- basename(x)
  x <- tools::file_path_sans_ext(x)
  # Keep "24110099601-1" pattern if found
  m <- stringr::str_match(x, "(\\d+-\\d+)")
  ifelse(!is.na(m[, 2]), m[, 2], x)
}

# 3. Strategy A: Start from MLST provenance if available
source_tbl <- NULL
if (file.exists(FILE_MLST_ALL)) {
  source_tbl <- readr::read_tsv(FILE_MLST_ALL, show_col_types = FALSE)
}

if (!is.null(source_tbl)) {
  # Try to derive Isolate_ID from common fields
  if (!"Isolate_ID" %in% names(source_tbl)) {
    if ("file_name" %in% names(source_tbl)) {
      source_tbl <- source_tbl %>%
        mutate(Isolate_ID = tools::file_path_sans_ext(basename(file_name)))
    } else if ("full_path" %in% names(source_tbl)) {
      source_tbl <- source_tbl %>%
        mutate(Isolate_ID = tools::file_path_sans_ext(basename(full_path)))
    }
  }
}

# 4. Strategy B: Scan FASTA directory
fa_files <- if (dir.exists(DIR_FASTAS)) {
  fs::dir_ls(DIR_FASTAS, regexp = "\\.(fa|fna|fasta)(\\.gz)?$", recurse = FALSE)
} else {
  character(0)
}

scan_tbl <- if (length(fa_files)) {
  tibble(full_path = fa_files) %>%
    mutate(
      file_name = basename(full_path),
      Isolate_ID = tools::file_path_sans_ext(file_name)
    )
} else {
  NULL
}

# 5. Merge Strategies
cand <- list(source_tbl, scan_tbl) %>%
  purrr::keep(~ !is.null(.x)) %>%
  purrr::reduce(dplyr::bind_rows, .init = NULL)

if (is.null(cand)) {
  stop("Could not find MLST outputs or FASTA files to seed assembly metadata.")
}

cand <- cand %>%
  mutate(Isolate_ID = norm_iso_id(Isolate_ID)) %>%
  distinct(Isolate_ID, .keep_all = TRUE)

# 6. Get Participant_id + Timepoint mapping
pid_tp <- NULL
map_file <- file.path(DIR_RESULTS, "assembly_pid_timepoint_map.csv")

if (file.exists(map_file)) {
  pid_tp <- readr::read_csv(map_file, show_col_types = FALSE) %>%
    mutate(
      Isolate_ID     = norm_iso_id(Isolate_ID),
      Participant_id = as.character(Participant_id),
      Timepoint      = as.character(Timepoint)
    )
} else if (file.exists(FILE_VF_HITS)) {
  # Fallback: vf_hits_all.rds
  hits <- readRDS(FILE_VF_HITS)
  link_col <- intersect(c("Isolate_ID", "file_name", "full_path"), names(hits))[1]

  if (!is.na(link_col)) {
    pid_tp <- hits %>%
      mutate(.link = norm_iso_id(.data[[link_col]])) %>%
      transmute(
        Isolate_ID     = .link,
        Participant_id = as.character(Participant_id),
        Timepoint      = as.character(if ("Timepoint" %in% names(.)) Timepoint else NA_character_)
      ) %>%
      distinct()
  }
}

# 7. Combine and Finalize
if (is.null(pid_tp)) {
  # No mapping found, initialize empty columns
  meta <- cand %>%
    mutate(Participant_id = NA_character_, Timepoint = NA_character_)
} else {
  meta <- cand %>%
    left_join(pid_tp, by = "Isolate_ID")
}

meta <- meta %>%
  # Ensure required columns exist
  mutate(
    Participant_id = if ("Participant_id" %in% names(.)) Participant_id else NA_character_,
    Timepoint      = if ("Timepoint" %in% names(.)) Timepoint else NA_character_
  ) %>%
  select(
    Isolate_ID,
    Participant_id,
    Timepoint,
    file_name = any_of("file_name"),
    full_path = any_of("full_path")
  ) %>%
  arrange(Participant_id, Timepoint, Isolate_ID)

# 8. Write Outputs
# Write to root (legacy compatibility) and results (cleaner)
readr::write_csv(meta, FILE_METADATA)
readr::write_csv(meta, file.path(DIR_RESULTS, "assembly_metadata.csv"))

# Check for missing data
todo <- meta %>% filter(is.na(Participant_id) | is.na(Timepoint))
if (nrow(todo)) {
  todo_file <- file.path(DIR_RESULTS, "assembly_metadata_TODO_fill.csv")
  readr::write_csv(todo, todo_file)
  message("⚠ Rows needing manual fill: ", todo_file, " (", nrow(todo), ")")
}

message("✓ Metadata generated: ", FILE_METADATA)

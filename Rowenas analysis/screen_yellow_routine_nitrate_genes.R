#!/usr/bin/env Rscript

# ==============================================================================
# USER SETTING
# ==============================================================================
# Set this to the YELLOW project folder on the computer running the script.
# Example:
# WORKING_DIRECTORY <- "D:/Batch 1-6 Ecoli FULL sequence YELLOW Study"
#
# The script then searches inside that folder for:
#   - raw_fasta/ containing longcycler FASTA files
#   - reference gene FASTAs such as CsgD.fasta, NarG.fasta, ...
#   - OVERVIEW E.coli batch 1-6 - CLEAN_RC.xlsx
#   - YELLOW dataset stage Rowena Studie meetmomenten long format_1.sav
# ==============================================================================
WORKING_DIRECTORY <- ""

# ==============================================================================
# OVERVIEW OF THE WORKFLOW
# ==============================================================================
# This script does four main things:
#
# 1. Find the project inputs from WORKING_DIRECTORY.
#    It does not rebuild metadata. It reads the existing overview workbook,
#    the existing SPSS nitrite file, the existing reference genes, and the
#    existing longcycler FASTA files.
#
# 2. Parse the longcycler FASTA filenames into isolate IDs.
#    The file names often contain sequencing run/barcode text that is not in
#    the metadata workbook. The parser removes that technical prefix so the
#    FASTA can be matched back to the overview workbook's isolate_id/Isolaat id.
#
# 3. Screen each longcycler assembly for the selected nitrate/nitrite genes.
#    This intentionally follows the older script logic and uses
#    Biostrings::matchPattern(..., max.mismatch = 60), rather than BLAST.
#
# 4. Join the gene presence table to dipstick nitrite results and run the
#    first-pass descriptive/Fisher/BH analysis.
#    Fisher exact tests are used because many selected genes are very common,
#    so absent groups can be too small for stable logistic regression.
# ==============================================================================

# Load packages quietly so that the script output focuses on useful progress
# messages and any real errors.
suppressPackageStartupMessages({
  library(Biostrings)
  library(dplyr)
  library(haven)
  library(readr)
  library(readxl)
  library(stringr)
  library(tidyr)
})

# Optional command-line settings. These are useful for testing without editing
# the script itself, for example:
#   Rscript screen_yellow_routine_nitrate_genes.R --max-isolates=5
#   Rscript screen_yellow_routine_nitrate_genes.R --working-directory="D:/..."
args <- commandArgs(trailingOnly = TRUE)

# Helper for reading arguments written as --name=value.
get_arg <- function(name, default = NULL) {
  hit <- args[str_detect(args, paste0("^", fixed(name), "="))]
  if (length(hit) == 0) {
    return(default)
  }
  str_remove(hit[[1]], paste0("^", fixed(name), "="))
}

max_isolates <- suppressWarnings(as.integer(get_arg("--max-isolates", NA_character_)))
max_mismatch <- suppressWarnings(as.integer(get_arg("--max-mismatch", "60")))
working_directory <- get_arg("--working-directory", WORKING_DIRECTORY)

# These are the genes that will be screened. The script expects one reference
# FASTA per gene, usually named CsgD.fasta, NarG.fasta, etc.
expected_gene_names <- c(
  "CsgD", "NarG", "NarH", "NarI", "NarJ", "NarL",
  "NarP", "NarQ", "NarV", "NarX", "NarY", "NarZ"
)
expected_reference_files <- paste0(expected_gene_names, ".fasta")

# Small message wrapper so progress output is consistent.
message_line <- function(...) {
  message(paste0(...))
}

# Used when a required file is missing. If there are several missing files,
# this prints them as a readable list.
stop_with_files <- function(message_text, files = character()) {
  if (length(files) > 0) {
    stop(paste0(message_text, "\n", paste(files, collapse = "\n")), call. = FALSE)
  }
  stop(message_text, call. = FALSE)
}

# Find the folder where this script lives. This is only used as a fallback
# search location for reference FASTAs if WORKING_DIRECTORY does not contain
# them directly.
get_script_dir <- function() {
  cmd <- commandArgs(trailingOnly = FALSE)
  file_arg <- cmd[startsWith(cmd, "--file=")]
  if (length(file_arg) > 0) {
    return(dirname(normalizePath(
      sub("^--file=", "", file_arg[[1]]),
      winslash = "/",
      mustWork = FALSE
    )))
  }
  normalizePath(getwd(), winslash = "/", mustWork = TRUE)
}

# Return a path plus its parent folders. This lets the script auto-detect the
# project folder when WORKING_DIRECTORY is left blank and the script is launched
# from somewhere inside the project.
parent_dirs <- function(path, max_depth = 8) {
  if (!dir.exists(path)) {
    return(character())
  }
  path <- normalizePath(path, winslash = "/", mustWork = TRUE)
  out <- path
  for (i in seq_len(max_depth)) {
    parent <- dirname(tail(out, 1))
    if (identical(parent, tail(out, 1))) {
      break
    }
    out <- c(out, parent)
  }
  unique(out)
}

# Score candidate folders based on the files/folders we expect in the YELLOW
# project. The highest-scoring folder is treated as the project root if the user
# has not set WORKING_DIRECTORY.
score_project_root <- function(path) {
  score <- 0L
  if (dir.exists(file.path(path, "raw_fasta"))) {
    score <- score + 5L
  }
  if (file.exists(file.path(path, "OVERVIEW E.coli batch 1-6 - CLEAN_RC.xlsx"))) {
    score <- score + 4L
  }
  if (file.exists(file.path(path, "YELLOW dataset stage Rowena Studie meetmomenten long format_1.sav"))) {
    score <- score + 4L
  }
  if (dir.exists(file.path(path, "reference_genes"))) {
    score <- score + 2L
  }
  score <- score + sum(file.exists(file.path(path, expected_reference_files)))
  score
}

# Choose the project root. If WORKING_DIRECTORY is provided, it wins. If it is
# blank, the script tries to infer the project root from the current/script
# folder and nearby parent folders.
detect_project_root <- function(working_directory) {
  if (!is.null(working_directory) && nzchar(working_directory)) {
    if (!dir.exists(working_directory)) {
      stop("WORKING_DIRECTORY does not exist: ", working_directory, call. = FALSE)
    }
    return(normalizePath(working_directory, winslash = "/", mustWork = TRUE))
  }

  script_dir <- get_script_dir()
  candidates <- unique(c(
    script_dir,
    getwd(),
    parent_dirs(script_dir),
    parent_dirs(getwd())
  ))
  candidates <- candidates[dir.exists(candidates)]

  scores <- vapply(candidates, score_project_root, integer(1))
  if (length(scores) == 0 || max(scores) == 0) {
    return(normalizePath(getwd(), winslash = "/", mustWork = TRUE))
  }
  normalizePath(candidates[which.max(scores)], winslash = "/", mustWork = TRUE)
}

# Find an exact file name anywhere inside the project folder. This is used for
# the overview workbook and the SPSS nitrite file. Matching is case-insensitive
# but otherwise exact, so it will not accidentally pick the wrong file.
find_named_file <- function(root, file_name, required = TRUE) {
  direct <- file.path(root, file_name)
  if (file.exists(direct)) {
    return(normalizePath(direct, winslash = "/", mustWork = TRUE))
  }

  found <- list.files(root, recursive = TRUE, full.names = TRUE)
  found <- found[tolower(basename(found)) == tolower(file_name)]
  found <- found[file.exists(found)]
  if (length(found) > 0) {
    return(normalizePath(found[[1]], winslash = "/", mustWork = TRUE))
  }

  if (required) {
    stop("Could not find required file inside WORKING_DIRECTORY: ", file_name, call. = FALSE)
  }
  NA_character_
}

# Locate the raw_fasta folder. The user said all FASTAs are in raw_fasta, so
# the script searches for that folder and then only uses longcycler FASTAs from it.
find_raw_fasta_dir <- function(root) {
  direct <- file.path(root, "raw_fasta")
  if (dir.exists(direct)) {
    return(normalizePath(direct, winslash = "/", mustWork = TRUE))
  }

  found <- list.files(root, pattern = "^raw_fasta$", recursive = TRUE, full.names = TRUE)
  found <- found[dir.exists(found)]
  if (length(found) > 0) {
    return(normalizePath(found[[1]], winslash = "/", mustWork = TRUE))
  }

  stop("Could not find a raw_fasta folder inside WORKING_DIRECTORY.", call. = FALSE)
}

# Locate the selected gene reference FASTAs. It first checks reference_genes/
# and references/ if present, then the project folder, then the script folder.
# The script stops early if any expected gene reference is missing.
find_reference_files <- function(root) {
  search_roots <- unique(c(
    file.path(root, "reference_genes"),
    file.path(root, "references"),
    root,
    get_script_dir()
  ))
  search_roots <- search_roots[dir.exists(search_roots)]

  out <- tibble(gene = expected_gene_names, reference_file = NA_character_)
  for (i in seq_along(expected_gene_names)) {
    wanted <- expected_reference_files[[i]]
    hits <- character()
    for (search_root in search_roots) {
      direct <- file.path(search_root, wanted)
      if (file.exists(direct)) {
        hits <- c(hits, direct)
      }
      nested <- list.files(search_root, recursive = TRUE, full.names = TRUE)
      hits <- c(hits, nested[tolower(basename(nested)) == tolower(wanted)])
    }
    hits <- unique(hits[file.exists(hits)])
    if (length(hits) > 0) {
      out$reference_file[[i]] <- normalizePath(hits[[1]], winslash = "/", mustWork = TRUE)
    }
  }

  missing <- out %>% filter(is.na(.data$reference_file))
  if (nrow(missing) > 0) {
    stop_with_files(
      "Could not find these required reference FASTA files:",
      missing$gene
    )
  }
  out
}

standardise_key <- function(x) {
  x %>%
    as.character() %>%
    str_replace_all("[^A-Za-z0-9]+", "_") %>%
    str_replace_all("_+", "_") %>%
    str_replace_all("^_|_$", "") %>%
    str_to_lower()
}

# Find the first column name that matches one of the allowed spellings.
# This is useful because the same column might be written as "Isolaat id",
# "isolate_id", "Isolate ID", or another close variant.
first_existing_col <- function(names_vector, candidates) {
  keys <- standardise_key(names_vector)
  candidate_keys <- standardise_key(candidates)
  idx <- match(candidate_keys, keys)
  idx <- idx[!is.na(idx)]
  if (length(idx) == 0) {
    return(NA_character_)
  }
  names_vector[[idx[[1]]]]
}

# Read the existing overview workbook and standardise only the columns this
# analysis needs. This does not rebuild metadata or create a new workbook.
read_overview_metadata <- function(path) {
  # The workbook can contain title rows before the actual column headers, so
  # first read a preview and detect the most likely header row.
  preview <- readxl::read_excel(
    path,
    col_names = FALSE,
    n_max = 40,
    .name_repair = "minimal"
  )

  # Header rows are scored by how many expected column names they contain.
  header_score <- function(row_values) {
    keys <- standardise_key(row_values)
    sum(keys %in% c("isolaat_id", "isolate_id", "isolateid", "isolate", "participant_id", "timepoint", "meetmoment"))
  }

  scores <- apply(as.data.frame(preview), 1, header_score)
  header_row <- which.max(scores)
  if (length(header_row) == 0 || scores[[header_row]] < 2) {
    stop(
      "Could not detect a metadata header row in the overview workbook. ",
      "Expected columns such as Isolaat id, participant_id and timepoint/meetmoment.",
      call. = FALSE
    )
  }

  # Re-read the workbook using the detected header row.
  metadata_raw <- readxl::read_excel(
    path,
    skip = header_row - 1,
    .name_repair = "unique_quiet"
  )

  # Detect required columns using multiple acceptable names.
  participant_col <- first_existing_col(
    names(metadata_raw),
    c("participant_id", "Participant id", "participant", "deelnemer_id", "studienummer")
  )
  isolate_col <- first_existing_col(
    names(metadata_raw),
    c("Isolaat id", "isolaat_id", "isolate_id", "isolate id", "Isolate_ID", "Isolate ID", "isolateid", "isolate")
  )
  timepoint_col <- first_existing_col(
    names(metadata_raw),
    c("timepoint", "Timepoint", "meetmoment", "Meetmoment", "tp_lab")
  )

  if (is.na(participant_col) || is.na(isolate_col) || is.na(timepoint_col)) {
    stop(
      "Overview workbook was found, but key columns were not detected. ",
      "Need participant_id, Isolaat id and timepoint/meetmoment.",
      call. = FALSE
    )
  }

  # Optional columns are retained if present, but filled with NA if absent.
  # This keeps the script robust to slightly different workbook versions.
  optional <- list(
    Clinical_Beoord = c("Clinical_Beoord", "clinical_beoord", "klinische_beoordeling"),
    Clinical_CFU_Count = c("Clinical_CFU_Count", "clinical_cfu_count", "cfu"),
    Urine_collection_method = c("Urine_collection_method", "urine_collection_method"),
    UTI_Label = c("UTI_Label", "uti_label"),
    Batch = c("Batch", "batch")
  )

  out <- metadata_raw %>%
    mutate(
      Participant_id = str_trim(as.character(.data[[participant_col]])),
      Isolate_ID = str_trim(as.character(.data[[isolate_col]])),
      Timepoint = str_trim(as.character(.data[[timepoint_col]]))
    )

  for (new_name in names(optional)) {
    old_name <- first_existing_col(names(metadata_raw), optional[[new_name]])
    out[[new_name]] <- if (is.na(old_name)) NA_character_ else as.character(metadata_raw[[old_name]])
  }

  # Convert the timepoint into tp_lab, the label used for joining to nitrite.
  # Routine SPSS-style timepoints 1-7 map to T0-T6.
  out %>%
    filter(!is.na(.data$Isolate_ID), nzchar(.data$Isolate_ID)) %>%
    mutate(
      tp_lab = case_when(
        str_detect(.data$Timepoint, regex("^T[0-6]$", ignore_case = TRUE)) ~ str_to_upper(.data$Timepoint),
        str_detect(.data$Timepoint, regex("^UTI", ignore_case = TRUE)) ~ str_replace_all(.data$Timepoint, "_", "-"),
        str_detect(.data$Timepoint, "^[1-7]$") ~ paste0("T", as.integer(.data$Timepoint) - 1L),
        str_detect(.data$Timepoint, "^[0-6]$") ~ paste0("T", .data$Timepoint),
        TRUE ~ .data$Timepoint
      ),
      Event_type = if_else(
        str_detect(.data$tp_lab, regex("^UTI", ignore_case = TRUE)),
        "UTI_event",
        "routine"
      )
    ) %>%
    select(
      Participant_id,
      Timepoint,
      tp_lab,
      Event_type,
      Isolate_ID,
      Batch,
      UTI_Label,
      Clinical_Beoord,
      Clinical_CFU_Count,
      Urine_collection_method,
      everything()
    ) %>%
    distinct(.data$Isolate_ID, .keep_all = TRUE)
}

# Convert a longcycler FASTA filename into the isolate ID that should match the
# overview workbook. For example:
#   APR1_barcode01_ABC-1_longcycler.fasta -> ABC-1
#   PR2-barcode4-ABC_longcycler.fa        -> ABC
parse_longcycler_isolate_id <- function(file_name) {
  file_name %>%
    str_remove(regex("[.](fasta|fa|fna)$", ignore_case = TRUE)) %>%
    str_remove(regex("^A?PR[0-9]+[-_ ]*barcode[-_ ]*[0-9]+[-_ ]*", ignore_case = TRUE)) %>%
    str_remove(regex("[-_ ]*longcycler$", ignore_case = TRUE)) %>%
    str_remove(regex("[-_ ]*longcycler[-_ ]*assembly$", ignore_case = TRUE)) %>%
    str_replace_all("_", "-") %>%
    str_replace_all("--+", "-") %>%
    str_replace_all("^-|-$", "") %>%
    str_trim()
}

# More forgiving isolate ID used only for joining. It removes punctuation and
# casing differences, so "ABC-1", "ABC_1", and "abc 1" can still match.
normalise_isolate_id_for_join <- function(x) {
  x %>%
    as.character() %>%
    str_to_lower() %>%
    str_remove(regex("[.](fasta|fa|fna)$", ignore_case = TRUE)) %>%
    str_remove(regex("^a?pr[0-9]+[-_ ]*barcode[-_ ]*[0-9]+[-_ ]*", ignore_case = TRUE)) %>%
    str_remove_all(regex("longcycler|flye|assembly", ignore_case = TRUE)) %>%
    str_replace_all("[^a-z0-9]", "")
}

# Some files may include a terminal suffix such as "-1" while another file uses
# the base ID only. This helper lets the join try both forms.
drop_terminal_isolate_suffix <- function(x) {
  x %>%
    str_remove(regex("[-_][0-9]+$", ignore_case = TRUE)) %>%
    str_replace_all("_", "-") %>%
    str_trim()
}

# List only longcycler FASTAs from raw_fasta. Other assemblers are deliberately
# ignored because this requested workflow is longcycler-only.
list_longcycler_fastas <- function(raw_fasta_dir) {
  fasta_files <- list.files(
    raw_fasta_dir,
    pattern = "[.](fasta|fa|fna)$",
    recursive = TRUE,
    full.names = TRUE,
    ignore.case = TRUE
  )
  fasta_files <- fasta_files[str_detect(basename(fasta_files), regex("longcycler", ignore_case = TRUE))]
  fasta_files <- normalizePath(fasta_files, winslash = "/", mustWork = TRUE)
  raw_prefix <- paste0(normalizePath(raw_fasta_dir, winslash = "/", mustWork = TRUE), "/")

  # The result is also an audit table: it records the original filename,
  # parsed isolate ID, sequencing run, sequencing barcode, and full path.
  tibble(
    full_path = fasta_files,
    relative_path = if_else(
      startsWith(fasta_files, raw_prefix),
      substring(fasta_files, nchar(raw_prefix) + 1L),
      basename(fasta_files)
    ),
    file_name = basename(fasta_files)
  ) %>%
    mutate(
      Isolate_ID_raw = parse_longcycler_isolate_id(.data$file_name),
      Isolate_ID_no_suffix = drop_terminal_isolate_suffix(.data$Isolate_ID_raw),
      Isolate_ID_join_key = normalise_isolate_id_for_join(.data$Isolate_ID_raw),
      Isolate_ID_no_suffix_join_key = normalise_isolate_id_for_join(.data$Isolate_ID_no_suffix),
      sequencing_run = str_extract(.data$file_name, regex("^A?PR[0-9]+", ignore_case = TRUE)),
      sequencing_barcode = str_match(.data$file_name, regex("barcode[-_ ]*0*([0-9]+)", ignore_case = TRUE))[, 2],
      run_barcode = if_else(
        !is.na(.data$sequencing_run) & !is.na(.data$sequencing_barcode),
        paste0(.data$sequencing_run, "_barcode", .data$sequencing_barcode),
        NA_character_
      ),
      Assembler = "longcycler"
    ) %>%
    arrange(.data$file_name)
}

# Match parsed FASTA isolate IDs to rows in the overview workbook.
#
# Matching is tried in stages from strict to forgiving:
#   1. exact_isolate_id
#   2. unique_id_without_terminal_suffix
#   3. normalised_isolate_id
#   4. normalised_id_without_terminal_suffix
#
# The selected_yellow_routine_fastas.csv output includes match_type so she can
# see exactly how each FASTA matched back to the workbook.
match_fastas_to_metadata <- function(longcycler_tbl, overview_tbl) {
  overview_lookup <- overview_tbl %>%
    mutate(
      Isolate_ID_overview = .data$Isolate_ID,
      Isolate_ID_no_suffix = drop_terminal_isolate_suffix(.data$Isolate_ID),
      Isolate_ID_join_key = normalise_isolate_id_for_join(.data$Isolate_ID),
      Isolate_ID_no_suffix_join_key = normalise_isolate_id_for_join(.data$Isolate_ID_no_suffix)
    )

  overview_columns <- names(overview_tbl)

  # Only allow forgiving matches if the overview key is unique. This prevents
  # one FASTA from being joined to multiple possible metadata rows.
  make_lookup <- function(key_column) {
    overview_lookup %>%
      filter(!is.na(.data[[key_column]]), nzchar(.data[[key_column]])) %>%
      add_count(.data[[key_column]], name = "overview_key_n") %>%
      filter(.data$overview_key_n == 1L) %>%
      select(any_of(unique(c(key_column, overview_columns, "Isolate_ID_overview"))))
  }

  # Run one matching stage, then pass only the still-unmatched FASTAs to the
  # next stage.
  run_match_stage <- function(remaining, left_key, right_key, match_type) {
    by_keys <- setNames(right_key, left_key)
    joined <- remaining %>%
      left_join(make_lookup(right_key), by = by_keys)

    list(
      matched = joined %>%
        filter(!is.na(.data$Participant_id)) %>%
        mutate(match_type = match_type),
      remaining = joined %>%
        filter(is.na(.data$Participant_id)) %>%
        select(all_of(names(remaining)))
    )
  }

  remaining <- longcycler_tbl %>%
    mutate(.fasta_row_id = row_number())

  # Order matters: exact matching is preferred before any forgiving matching.
  stages <- list(
    c("Isolate_ID_raw", "Isolate_ID_overview", "exact_isolate_id"),
    c("Isolate_ID_no_suffix", "Isolate_ID_no_suffix", "unique_id_without_terminal_suffix"),
    c("Isolate_ID_join_key", "Isolate_ID_join_key", "normalised_isolate_id"),
    c("Isolate_ID_no_suffix_join_key", "Isolate_ID_no_suffix_join_key", "normalised_id_without_terminal_suffix")
  )

  matched_list <- list()
  for (stage in stages) {
    stage_result <- run_match_stage(
      remaining = remaining,
      left_key = stage[[1]],
      right_key = stage[[2]],
      match_type = stage[[3]]
    )
    matched_list <- c(matched_list, list(stage_result$matched))
    remaining <- stage_result$remaining
  }

  bind_rows(matched_list, remaining %>% mutate(match_type = NA_character_)) %>%
    arrange(.data$.fasta_row_id) %>%
    select(-.fasta_row_id) %>%
    mutate(
      Isolate_ID = coalesce(.data$Isolate_ID, .data$Isolate_ID_raw),
      fasta_matched_overview = !is.na(.data$Participant_id)
    )
}

# Read the first sequence from a reference FASTA file. The original workflow
# assumes one reference sequence per gene FASTA, so only the first sequence is used.
read_first_reference_sequence <- function(path) {
  ref <- readDNAStringSet(path)
  if (length(ref) == 0) {
    stop("Reference FASTA has no sequences: ", path, call. = FALSE)
  }
  as.character(ref[[1]])
}

# Search one assembled genome for one reference gene sequence.
#
# This is the old-script approach:
#   - read all contigs from the assembly FASTA
#   - paste them together into one genome string
#   - search for the reference gene using matchPattern()
#   - allow max_mismatch differences, default 60
#
# If no match is found, NA is returned and the gene is treated as absent.
extract_gene_from_genome <- function(genome_file, gene_sequence, max_mismatch = 60) {
  genome <- readDNAStringSet(genome_file)
  if (length(genome) == 0) {
    return(NA_character_)
  }

  genome_sequence <- paste(as.character(genome), collapse = "")
  matches <- matchPattern(gene_sequence, genome_sequence, max.mismatch = max_mismatch)
  if (length(matches) == 0) {
    return(NA_character_)
  }

  as.character(matches[[1]])
}

# Run the gene screen for every matched longcycler FASTA and every selected
# reference gene. This produces a long table with one row per isolate-gene pair.
screen_gene_presence <- function(fasta_tbl, reference_tbl, max_mismatch) {
  # Preload reference sequences once so they are not reread for every isolate.
  reference_tbl <- reference_tbl %>%
    mutate(reference_sequence = vapply(.data$reference_file, read_first_reference_sequence, character(1)))

  out <- vector("list", nrow(fasta_tbl) * nrow(reference_tbl))
  k <- 1L

  for (i in seq_len(nrow(fasta_tbl))) {
    if (i %% 25 == 0 || i == 1 || i == nrow(fasta_tbl)) {
      message_line("Screening isolate ", i, " of ", nrow(fasta_tbl), ": ", fasta_tbl$file_name[[i]])
    }

    for (j in seq_len(nrow(reference_tbl))) {
      # If one isolate/gene search fails for a file-reading reason, record NA
      # and continue rather than losing the whole run.
      sequence <- tryCatch(
        extract_gene_from_genome(
          genome_file = fasta_tbl$full_path[[i]],
          gene_sequence = reference_tbl$reference_sequence[[j]],
          max_mismatch = max_mismatch
        ),
        error = function(e) {
          warning(
            "Could not screen ", basename(fasta_tbl$full_path[[i]]),
            " for ", reference_tbl$gene[[j]], ": ", conditionMessage(e),
            call. = FALSE
          )
          NA_character_
        }
      )

      # Store both the gene presence result and the extracted sequence. The
      # sequence is needed for allele calling; sequence_found becomes gene_*
      # presence/absence later.
      out[[k]] <- tibble(
        Isolate_ID = fasta_tbl$Isolate_ID[[i]],
        Isolate_ID_raw = fasta_tbl$Isolate_ID_raw[[i]],
        Participant_id = fasta_tbl$Participant_id[[i]],
        Timepoint = fasta_tbl$Timepoint[[i]],
        tp_lab = fasta_tbl$tp_lab[[i]],
        Event_type = fasta_tbl$Event_type[[i]],
        sequencing_run = fasta_tbl$sequencing_run[[i]],
        sequencing_barcode = fasta_tbl$sequencing_barcode[[i]],
        run_barcode = fasta_tbl$run_barcode[[i]],
        file_name = fasta_tbl$file_name[[i]],
        full_path = fasta_tbl$full_path[[i]],
        gene = reference_tbl$gene[[j]],
        reference_file = basename(reference_tbl$reference_file[[j]]),
        max_mismatch = max_mismatch,
        sequence_found = !is.na(sequence),
        sequence = sequence
      )
      k <- k + 1L
    }
  }

  bind_rows(out)
}

# Assign allele names within each gene. For example, if NarG has three unique
# extracted sequences, they become NarG_Allele1, NarG_Allele2, NarG_Allele3.
assign_gene_alleles <- function(extraction_tbl) {
  present <- extraction_tbl %>%
    filter(.data$sequence_found)

  if (nrow(present) == 0) {
    return(tibble())
  }

  present %>%
    group_by(gene) %>%
    group_modify(function(.x, .y) {
      unique_sequences <- sort(unique(.x$sequence))
      allele_map <- tibble(
        sequence = unique_sequences,
        allele = paste0(.y$gene[[1]], "_Allele", seq_along(unique_sequences))
      )
      .x %>%
        left_join(allele_map, by = "sequence")
    }) %>%
    ungroup()
}

# Convert the long isolate-gene table into one row per isolate with binary gene
# columns: gene_CsgD, gene_NarG, etc. These are the columns tested against
# dipstick nitrite status.
build_presence_wide <- function(extraction_tbl) {
  extraction_tbl %>%
    mutate(
      gene_col = paste0("gene_", .data$gene),
      value = as.integer(.data$sequence_found)
    ) %>%
    select(
      Isolate_ID,
      Isolate_ID_raw,
      Participant_id,
      Timepoint,
      tp_lab,
      Event_type,
      sequencing_run,
      sequencing_barcode,
      run_barcode,
      file_name,
      gene_col,
      value
    ) %>%
    distinct() %>%
    pivot_wider(
      names_from = gene_col,
      values_from = value,
      values_fill = 0
    )
}

# Convert per-gene allele calls into one row per isolate. The genotype string is
# a compact concatenation of all selected gene allele calls, using "missing" for
# genes without an extracted sequence.
build_profile_wide <- function(allele_tbl) {
  if (nrow(allele_tbl) == 0) {
    return(tibble())
  }

  profile <- allele_tbl %>%
    select(Isolate_ID, gene, allele) %>%
    distinct() %>%
    pivot_wider(names_from = gene, values_from = allele)

  gene_cols <- intersect(expected_gene_names, names(profile))
  if (length(gene_cols) == 0) {
    profile$genotype <- NA_character_
  } else {
    profile$genotype <- apply(
      profile[, gene_cols, drop = FALSE],
      1,
      function(x) paste(replace_na(as.character(x), "missing"), collapse = "_")
    )
  }
  profile
}

# Convert raw nitrite values into the labels used for analysis.
# Positive can be SPSS code 1 or text such as Positief/positive/pos.
# Negative can be SPSS code 2 or text such as Negatief/negative/neg.
# Anything else is treated as missing and excluded from Fisher tests.
normalise_nitrite_value <- function(x) {
  x_chr <- as.character(x)
  x_num <- suppressWarnings(as.numeric(x_chr))
  x_low <- str_to_lower(str_trim(x_chr))

  case_when(
    !is.na(x_num) & x_num == 1 ~ "positive",
    !is.na(x_num) & x_num == 2 ~ "negative",
    str_detect(x_low, "positief|positive|pos") ~ "positive",
    str_detect(x_low, "negatief|negative|neg") ~ "negative",
    TRUE ~ "missing"
  )
}

# Keep the human-readable label from an SPSS labelled vector when available.
# This is useful in audit outputs because code 1/2 alone is less readable.
nitrite_label_text <- function(x) {
  if (inherits(x, "haven_labelled") || inherits(x, "labelled")) {
    return(as.character(haven::as_factor(x, levels = "labels")))
  }
  as.character(x)
}

# Read the routine nitrite SPSS .sav file and convert its timepoint to tp_lab.
# Routine timepoints 1-7 are mapped to T0-T6 so they can join to the FASTA rows.
read_routine_nitrite <- function(path) {
  raw <- haven::read_sav(path)

  participant_col <- first_existing_col(names(raw), c("participant_id", "Participant id"))
  timepoint_col <- first_existing_col(names(raw), c("timepoint", "Timepoint", "meetmoment"))
  nitrite_col <- first_existing_col(names(raw), c("urinestick_nitriet", "dipstick_nitrite", "nitrite"))

  if (is.na(participant_col) || is.na(timepoint_col) || is.na(nitrite_col)) {
    stop(
      "Nitrite file was found, but required columns were not detected. ",
      "Need participant_id, timepoint and urinestick_nitriet.",
      call. = FALSE
    )
  }

  raw %>%
    transmute(
      Participant_id = str_trim(as.character(.data[[participant_col]])),
      spss_timepoint = as.character(.data[[timepoint_col]]),
      spss_timepoint_integer = suppressWarnings(as.integer(as.character(.data[[timepoint_col]]))),
      tp_lab = case_when(
        !is.na(.data$spss_timepoint_integer) &
          .data$spss_timepoint_integer >= 1L &
          .data$spss_timepoint_integer <= 7L ~ paste0("T", .data$spss_timepoint_integer - 1L),
        str_detect(.data$spss_timepoint, regex("^T[0-6]$", ignore_case = TRUE)) ~ str_to_upper(.data$spss_timepoint),
        TRUE ~ .data$spss_timepoint
      ),
      urinestick_nitriet_code = as.character(.data[[nitrite_col]]),
      urinestick_nitriet_label = nitrite_label_text(.data[[nitrite_col]]),
      nitrite_binary = normalise_nitrite_value(.data[[nitrite_col]])
    ) %>%
    distinct(Participant_id, tp_lab, .keep_all = TRUE)
}

# Optional UTI nitrite support. If the UTI file is not present, this returns an
# empty tibble and the routine analysis continues normally.
read_optional_uti_nitrite <- function(path) {
  if (is.na(path) || !file.exists(path)) {
    return(tibble())
  }

  raw <- haven::read_sav(path)
  participant_col <- first_existing_col(names(raw), c("participant_id", "Participant id"))
  timepoint_col <- first_existing_col(names(raw), c("timepoint", "Timepoint", "meetmoment"))
  nitrite_col <- first_existing_col(names(raw), c("urinestick_nitriet", "dipstick_nitrite", "nitrite"))

  if (is.na(participant_col) || is.na(timepoint_col) || is.na(nitrite_col)) {
    warning("Optional UTI nitrite file was found, but required columns were not detected. Skipping it.", call. = FALSE)
    return(tibble())
  }

  raw %>%
    transmute(
      Participant_id = str_trim(as.character(.data[[participant_col]])),
      spss_timepoint = as.character(.data[[timepoint_col]]),
      spss_timepoint_integer = suppressWarnings(as.integer(as.character(.data[[timepoint_col]]))),
      tp_lab = case_when(
        !is.na(.data$spss_timepoint_integer) ~ paste0("UTI-", .data$spss_timepoint_integer),
        str_detect(.data$spss_timepoint, regex("^UTI", ignore_case = TRUE)) ~ str_replace_all(.data$spss_timepoint, "_", "-"),
        TRUE ~ paste0("UTI-", .data$spss_timepoint)
      ),
      urinestick_nitriet_code = as.character(.data[[nitrite_col]]),
      urinestick_nitriet_label = nitrite_label_text(.data[[nitrite_col]]),
      nitrite_binary = normalise_nitrite_value(.data[[nitrite_col]])
    ) %>%
    distinct(Participant_id, tp_lab, .keep_all = TRUE)
}

# After joining nitrite to the FASTA/gene table, label each row so unmatched
# nitrite results can be audited instead of silently disappearing.
add_nitrite_status <- function(joined_tbl) {
  joined_tbl %>%
    mutate(
      nitrite_join_status = case_when(
        .data$nitrite_binary %in% c("positive", "negative") ~ "matched_binary_result",
        .data$nitrite_binary == "missing" ~ "matched_missing_or_unclear_result",
        TRUE ~ "no_matching_nitrite_row"
      ),
      nitrite_positive = case_when(
        .data$nitrite_binary == "positive" ~ TRUE,
        .data$nitrite_binary == "negative" ~ FALSE,
        TRUE ~ NA
      )
    )
}

# Compare one binary gene/module feature against dipstick nitrite status.
# The output contains the 2x2 table counts, percentages, and Fisher p-value.
summarise_feature_vs_nitrite <- function(analysis_rows, feature_col) {
  tmp <- analysis_rows %>%
    mutate(
      feature_present = as.integer(coalesce(as.integer(.data[[feature_col]]), 0L)),
      nitrite_binary = factor(.data$nitrite_binary, levels = c("negative", "positive"))
    )

  present_positive <- sum(tmp$feature_present == 1 & tmp$nitrite_binary == "positive", na.rm = TRUE)
  present_negative <- sum(tmp$feature_present == 1 & tmp$nitrite_binary == "negative", na.rm = TRUE)
  absent_positive <- sum(tmp$feature_present == 0 & tmp$nitrite_binary == "positive", na.rm = TRUE)
  absent_negative <- sum(tmp$feature_present == 0 & tmp$nitrite_binary == "negative", na.rm = TRUE)

  # This creates the 2x2 table tested by fisher.test():
  #             nitrite positive   nitrite negative
  # present     present_positive   present_negative
  # absent      absent_positive    absent_negative
  contingency <- matrix(
    c(present_positive, present_negative, absent_positive, absent_negative),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(
      feature = c("present", "absent"),
      nitrite = c("positive", "negative")
    )
  )

  fisher_p <- tryCatch(
    fisher.test(contingency)$p.value,
    error = function(e) NA_real_
  )

  tibble(
    feature = feature_col,
    present_n = present_positive + present_negative,
    absent_n = absent_positive + absent_negative,
    present_positive = present_positive,
    present_negative = present_negative,
    absent_positive = absent_positive,
    absent_negative = absent_negative,
    present_positive_percent = if_else(
      present_positive + present_negative > 0,
      100 * present_positive / (present_positive + present_negative),
      NA_real_
    ),
    absent_positive_percent = if_else(
      absent_positive + absent_negative > 0,
      100 * absent_positive / (absent_positive + absent_negative),
      NA_real_
    ),
    fisher_p_value = fisher_p
  )
}

# Build a long descriptive table with denominators and percentages for every
# selected gene. This is often the most useful table when genes are nearly
# universal and formal tests have sparse absent groups.
build_descriptive_counts <- function(analysis_rows, feature_cols) {
  if (length(feature_cols) == 0 || nrow(analysis_rows) == 0) {
    return(tibble())
  }

  analysis_rows %>%
    select(nitrite_binary, all_of(feature_cols)) %>%
    pivot_longer(
      cols = all_of(feature_cols),
      names_to = "feature",
      values_to = "present"
    ) %>%
    mutate(
      present = if_else(coalesce(as.integer(.data$present), 0L) == 1L, "present", "absent")
    ) %>%
    count(feature, present, nitrite_binary, name = "n") %>%
    group_by(feature, present) %>%
    mutate(
      denominator = sum(.data$n),
      percent = if_else(.data$denominator > 0, 100 * .data$n / .data$denominator, NA_real_)
    ) %>%
    ungroup()
}

# Write a short plain-English summary file. This is not used by the code, but is
# useful for quickly interpreting the Fisher/BH output.
write_summary_text <- function(path, fisher_tbl, analysis_rows) {
  n_analysis <- nrow(analysis_rows)
  nitrite_counts <- analysis_rows %>%
    count(nitrite_binary, name = "n") %>%
    arrange(nitrite_binary)

  best <- fisher_tbl %>%
    filter(!is.na(.data$fisher_p_value)) %>%
    arrange(.data$fisher_p_value) %>%
    slice_head(n = 1)

  bh_hits <- fisher_tbl %>%
    filter(!is.na(.data$bh_q_value), .data$bh_q_value < 0.05)

  lines <- c(
    "Nitrite/gene association summary",
    "",
    paste0("Analyzable rows with binary nitrite results: ", n_analysis),
    paste0(
      "Nitrite counts: ",
      paste0(nitrite_counts$nitrite_binary, "=", nitrite_counts$n, collapse = ", ")
    ),
    "",
    "Fisher exact tests were used because most selected nitrate/nitrite genes are expected to be very common and absent groups can be very small.",
    "Benjamini-Hochberg correction was used to account for testing multiple genes.",
    "Logistic regression is not the preferred first-pass method here because sparse absent groups can make model estimates unstable.",
    "These results should be interpreted as exploratory/descriptive because some participants contribute repeated timepoints.",
    ""
  )

  if (nrow(best) > 0) {
    lines <- c(
      lines,
      paste0(
        "Top unadjusted feature: ", best$feature[[1]],
        " (Fisher p=", signif(best$fisher_p_value[[1]], 4),
        ", BH q=", signif(best$bh_q_value[[1]], 4), ")."
      )
    )
  }

  if (nrow(bh_hits) == 0) {
    lines <- c(
      lines,
      "No BH-adjusted q-values were below 0.05, so there is no clear evidence that the selected genes are individually associated with dipstick nitrite status in this first-pass analysis."
    )
  } else {
    lines <- c(
      lines,
      paste0(
        "Features with BH-adjusted q-values below 0.05: ",
        paste(bh_hits$feature, collapse = ", "),
        "."
      )
    )
  }

  writeLines(lines, path)
}

# ==============================================================================
# MAIN WORKFLOW
# ==============================================================================

# 1. Work out where the project is and where outputs should be written.
project_root <- detect_project_root(working_directory)
output_dir <- file.path(project_root, "outputs", "nitrate_blast")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

message_line("Using project root: ", project_root)
message_line("Writing outputs to: ", output_dir)

# 2. Find the required input files/folders inside WORKING_DIRECTORY.
#    The metadata workbook is read as-is. It is not rebuilt.
raw_fasta_dir <- find_raw_fasta_dir(project_root)
reference_tbl <- find_reference_files(project_root)
overview_path <- find_named_file(project_root, "OVERVIEW E.coli batch 1-6 - CLEAN_RC.xlsx", required = TRUE)
routine_nitrite_path <- find_named_file(
  project_root,
  "YELLOW dataset stage Rowena Studie meetmomenten long format_1.sav",
  required = FALSE
)
if (is.na(routine_nitrite_path)) {
  # Keep a fallback for the older filename that had a space before "_1.sav".
  routine_nitrite_path <- find_named_file(
    project_root,
    "YELLOW dataset stage Rowena Studie meetmomenten long format _1.sav",
    required = FALSE
  )
}
if (is.na(routine_nitrite_path)) {
  stop(
    "Could not find the routine nitrite SPSS file inside WORKING_DIRECTORY: ",
    "YELLOW dataset stage Rowena Studie meetmomenten long format_1.sav",
    call. = FALSE
  )
}
optional_uti_nitrite_path <- find_named_file(
  project_root,
  "YELLOW dataset stage Rowena NEW dipslide long format.sav",
  required = FALSE
)

message_line("raw_fasta folder: ", raw_fasta_dir)
message_line("Overview workbook: ", overview_path)
message_line("Routine nitrite file: ", routine_nitrite_path)
if (!is.na(optional_uti_nitrite_path)) {
  message_line("Optional UTI nitrite file: ", optional_uti_nitrite_path)
}

# 3. List longcycler FASTAs and parse their isolate IDs from filenames.
longcycler_fastas <- list_longcycler_fastas(raw_fasta_dir)
if (nrow(longcycler_fastas) == 0) {
  stop("No longcycler FASTA files were found in raw_fasta.", call. = FALSE)
}

# For testing, --max-isolates lets her run only the first few FASTAs quickly.
if (!is.na(max_isolates) && max_isolates > 0) {
  longcycler_fastas <- longcycler_fastas %>% slice_head(n = max_isolates)
  message_line("Limiting run to --max-isolates=", max_isolates)
}

# 4. Read the overview workbook and match FASTAs to overview isolate IDs.
overview_metadata <- read_overview_metadata(overview_path)
matched_fastas <- match_fastas_to_metadata(longcycler_fastas, overview_metadata)

# Audit files are deliberately written before stopping. If matching fails, these
# files show exactly what isolate IDs were parsed and which metadata rows did
# not match.
write_csv(longcycler_fastas, file.path(output_dir, "audit_longcycler_fastas_found.csv"), na = "")
write_csv(
  longcycler_fastas %>% select(file_name, Isolate_ID_raw, Isolate_ID_no_suffix, everything()),
  file.path(output_dir, "audit_longcycler_parsed_isolate_ids.csv"),
  na = ""
)
write_csv(
  matched_fastas %>% filter(!.data$fasta_matched_overview),
  file.path(output_dir, "audit_longcycler_fastas_unmatched_to_overview.csv"),
  na = ""
)
write_csv(
  overview_metadata %>%
    anti_join(matched_fastas %>% filter(.data$fasta_matched_overview), by = "Isolate_ID"),
  file.path(output_dir, "audit_overview_rows_without_longcycler_fasta.csv"),
  na = ""
)

matched_for_screening <- matched_fastas %>%
  filter(.data$fasta_matched_overview) %>%
  arrange(.data$Isolate_ID, .data$file_name)

write_csv(matched_for_screening, file.path(output_dir, "selected_yellow_routine_fastas.csv"), na = "")

if (nrow(matched_for_screening) == 0) {
  stop(
    "Longcycler FASTAs were found, but none matched overview Isolaat id. ",
    "Check audit_longcycler_parsed_isolate_ids.csv and audit_longcycler_fastas_unmatched_to_overview.csv.",
    call. = FALSE
  )
}

message_line("Longcycler FASTAs found: ", nrow(longcycler_fastas))
message_line("Longcycler FASTAs matched to overview: ", nrow(matched_for_screening))

# 5. Screen matched longcycler FASTAs against the selected reference genes.
#    This is where the gene_CsgD/gene_NarG/etc. columns ultimately come from.
extraction_tbl <- screen_gene_presence(matched_for_screening, reference_tbl, max_mismatch)
allele_tbl <- assign_gene_alleles(extraction_tbl)
presence_wide <- build_presence_wide(extraction_tbl)
profile_wide <- build_profile_wide(allele_tbl)

# Ensure every expected gene column exists, even if a gene was absent from every
# isolate. This prevents downstream Fisher/BH code from failing with missing
# gene_ columns.
for (gene_col in paste0("gene_", expected_gene_names)) {
  if (!gene_col %in% names(presence_wide)) {
    presence_wide[[gene_col]] <- 0L
  }
}

gene_cols <- paste0("gene_", expected_gene_names)

# 6. Combine metadata, FASTA information, gene presence/absence, and genotype
#    profiles into the main trend-ready table.
trend_ready <- matched_for_screening %>%
  select(
    Isolate_ID,
    Isolate_ID_raw,
    Participant_id,
    Timepoint,
    tp_lab,
    Event_type,
    Batch,
    UTI_Label,
    Clinical_Beoord,
    Clinical_CFU_Count,
    Urine_collection_method,
    sequencing_run,
    sequencing_barcode,
    run_barcode,
    Assembler,
    file_name,
    full_path,
    match_type
  ) %>%
  left_join(
    presence_wide %>%
      select(Isolate_ID, all_of(gene_cols)),
    by = "Isolate_ID"
  ) %>%
  left_join(profile_wide, by = "Isolate_ID") %>%
  mutate(
    # Summary features used for descriptive biological interpretation.
    genes_screened = length(expected_gene_names),
    genes_present = rowSums(across(all_of(gene_cols)), na.rm = TRUE),
    genes_absent = .data$genes_screened - .data$genes_present,
    all_screened_genes_present = .data$genes_present == .data$genes_screened,
    narA_complete = .data$gene_NarG == 1L & .data$gene_NarH == 1L & .data$gene_NarI == 1L & .data$gene_NarJ == 1L,
    narZ_complete_available_refs = .data$gene_NarZ == 1L & .data$gene_NarY == 1L & .data$gene_NarV == 1L,
    nar_regulators_complete = .data$gene_NarL == 1L & .data$gene_NarP == 1L & .data$gene_NarQ == 1L & .data$gene_NarX == 1L,
    csgD_present = .data$gene_CsgD == 1L,
    nitrate_profile_01 = paste0(
      "narA=", as.integer(.data$narA_complete),
      ";narZ=", as.integer(.data$narZ_complete_available_refs),
      ";reg=", as.integer(.data$nar_regulators_complete),
      ";csgD=", as.integer(.data$csgD_present)
    )
  )

# Gene-level counts and absence checklist help with QC. For example, if a gene
# is absent in almost everyone, she can spot that before interpreting Fisher tests.
gene_counts <- extraction_tbl %>%
  count(gene, sequence_found, name = "n") %>%
  group_by(gene) %>%
  mutate(total_isolates = sum(.data$n), percent = 100 * .data$n / .data$total_isolates) %>%
  ungroup()

absence_checklist <- extraction_tbl %>%
  filter(!.data$sequence_found) %>%
  select(
    Isolate_ID,
    Participant_id,
    tp_lab,
    Event_type,
    file_name,
    gene,
    reference_file
  )

# Full extraction table: one row for every isolate-gene pair, including genes
# where no sequence was found.
write_csv(extraction_tbl, file.path(output_dir, "sequence_extraction_output_ref_vs_fasta_all_batches.csv"), na = "")

# Cleaned extraction table: only rows where an actual sequence was extracted.
# This is easier to inspect manually because absent genes are removed.
write_csv(
  extraction_tbl %>%
    filter(.data$sequence_found, !is.na(.data$sequence), nzchar(.data$sequence)),
  file.path(output_dir, "sequence_extraction_output_ref_vs_fasta_all_batches_clean.csv"),
  na = ""
)

write_csv(allele_tbl, file.path(output_dir, "Allele_typing_output_allele_profile_per_gene_all_batches.csv"), na = "")
write_csv(profile_wide, file.path(output_dir, "Allele_typing_output_allele_profile_per_isolate_all_batches.csv"), na = "")
write_csv(
  profile_wide %>% count(genotype, name = "n_isolates") %>% arrange(desc(n_isolates)),
  file.path(output_dir, "Genotype_counts_all_batches.csv"),
  na = ""
)
write_csv(
  trend_ready %>% select(Isolate_ID, Participant_id, tp_lab, genotype, everything()),
  file.path(output_dir, "Genotype_isolate_combination_overview_all_batches.csv"),
  na = ""
)
write_csv(extraction_tbl, file.path(output_dir, "nitrate_gene_presence_long.csv"), na = "")
write_csv(presence_wide, file.path(output_dir, "nitrate_gene_presence_by_isolate.csv"), na = "")
write_csv(gene_counts, file.path(output_dir, "nitrate_gene_presence_counts.csv"), na = "")
write_csv(absence_checklist, file.path(output_dir, "nitrate_gene_absence_checklist.csv"), na = "")
write_csv(trend_ready, file.path(output_dir, "nitrate_trend_ready_needs_nitrite.csv"), na = "")

# 7. Read nitrite results and join them by Participant_id + tp_lab.
routine_nitrite <- read_routine_nitrite(routine_nitrite_path)
uti_nitrite <- read_optional_uti_nitrite(optional_uti_nitrite_path)
nitrite_clean <- bind_rows(routine_nitrite, uti_nitrite) %>%
  distinct(Participant_id, tp_lab, .keep_all = TRUE)

# Duplicate nitrite keys would make a join ambiguous, so they are written as an
# audit table even though the script keeps the first distinct row.
duplicate_nitrite_keys <- bind_rows(routine_nitrite, uti_nitrite) %>%
  count(Participant_id, tp_lab, name = "n") %>%
  filter(.data$n > 1)
write_csv(duplicate_nitrite_keys, file.path(output_dir, "audit_duplicate_nitrite_keys.csv"), na = "")

nitrite_joined <- trend_ready %>%
  left_join(nitrite_clean, by = c("Participant_id", "tp_lab")) %>%
  add_nitrite_status()

write_csv(nitrite_joined, file.path(output_dir, "nitrate_barcode_nitrite_joined.csv"), na = "")
write_csv(
  nitrite_joined %>% filter(.data$nitrite_join_status == "no_matching_nitrite_row"),
  file.path(output_dir, "nitrate_barcode_nitrite_unmatched.csv"),
  na = ""
)
write_csv(nitrite_joined, file.path(output_dir, "nitrate_gene_nitrite_joined_for_fisher.csv"), na = "")
write_csv(
  nitrite_joined %>%
    count(nitrite_join_status, Event_type, name = "n"),
  file.path(output_dir, "audit_nitrite_join_status.csv"),
  na = ""
)

# Only rows with binary positive/negative nitrite results are used for Fisher/BH.
# Missing and unmatched rows stay in audit outputs but are excluded from tests.
analysis_rows <- nitrite_joined %>%
  filter(.data$nitrite_binary %in% c("positive", "negative"))

if (length(gene_cols) == 0 || !all(gene_cols %in% names(analysis_rows))) {
  stop("No gene_ columns found for Fisher/BH analysis after gene screening.", call. = FALSE)
}

# 8. Test each selected gene for association with nitrite status.
#    p-values are Benjamini-Hochberg adjusted because multiple genes are tested.
gene_fisher <- bind_rows(lapply(gene_cols, function(feature) {
  summarise_feature_vs_nitrite(analysis_rows, feature)
})) %>%
  mutate(
    bh_q_value = p.adjust(.data$fisher_p_value, method = "BH"),
    tested_rows = nrow(analysis_rows)
  ) %>%
  arrange(.data$bh_q_value, .data$fisher_p_value, .data$feature)

# Module/grouped features are useful biological summaries. These are exploratory
# and should be interpreted alongside the individual gene results.
module_cols <- c(
  "all_screened_genes_present",
  "narA_complete",
  "narZ_complete_available_refs",
  "nar_regulators_complete",
  "csgD_present"
)
module_fisher <- bind_rows(lapply(module_cols, function(feature) {
  summarise_feature_vs_nitrite(
    analysis_rows %>% mutate(across(all_of(feature), as.integer)),
    feature
  )
})) %>%
  mutate(
    bh_q_value = p.adjust(.data$fisher_p_value, method = "BH"),
    tested_rows = nrow(analysis_rows)
  ) %>%
  arrange(.data$bh_q_value, .data$fisher_p_value, .data$feature)

descriptive_counts <- build_descriptive_counts(analysis_rows, gene_cols)

# 9. Write the association outputs.
write_csv(gene_fisher, file.path(output_dir, "gene_vs_nitrite_fisher_bh.csv"), na = "")
write_csv(module_fisher, file.path(output_dir, "module_vs_nitrite_fisher_bh.csv"), na = "")
write_csv(descriptive_counts, file.path(output_dir, "gene_vs_nitrite_descriptive_counts.csv"), na = "")

# 10. Write barcode/run/profile descriptive summaries. These are descriptive QC
#     summaries, not formal statistics, because run-barcode IDs are usually
#     singleton technical identifiers.
write_csv(
  analysis_rows %>%
    group_by(sequencing_run) %>%
    summarise(
      n = n(),
      positive_n = sum(.data$nitrite_binary == "positive"),
      negative_n = sum(.data$nitrite_binary == "negative"),
      positive_percent = 100 * .data$positive_n / .data$n,
      .groups = "drop"
    ),
  file.path(output_dir, "nitrite_by_sequencing_run.csv"),
  na = ""
)
write_csv(
  analysis_rows %>%
    group_by(sequencing_barcode) %>%
    summarise(
      n = n(),
      positive_n = sum(.data$nitrite_binary == "positive"),
      negative_n = sum(.data$nitrite_binary == "negative"),
      positive_percent = 100 * .data$positive_n / .data$n,
      .groups = "drop"
    ),
  file.path(output_dir, "nitrite_by_barcode_number.csv"),
  na = ""
)
write_csv(
  analysis_rows %>%
    group_by(run_barcode) %>%
    summarise(
      n = n(),
      positive_n = sum(.data$nitrite_binary == "positive"),
      negative_n = sum(.data$nitrite_binary == "negative"),
      positive_percent = 100 * .data$positive_n / .data$n,
      .groups = "drop"
    ),
  file.path(output_dir, "nitrite_by_run_barcode.csv"),
  na = ""
)
write_csv(
  analysis_rows %>%
    group_by(nitrate_profile_01) %>%
    summarise(
      n = n(),
      positive_n = sum(.data$nitrite_binary == "positive"),
      negative_n = sum(.data$nitrite_binary == "negative"),
      positive_percent = 100 * .data$positive_n / .data$n,
      .groups = "drop"
    ),
  file.path(output_dir, "nitrite_by_nitrate_gene_profile.csv"),
  na = ""
)

write_summary_text(
  file.path(output_dir, "nitrite_gene_association_summary.txt"),
  gene_fisher,
  analysis_rows
)

# 11. Final QC summary. This is the quickest file to inspect after running the
#     script because it reports how many FASTAs were found, how many matched the
#     overview workbook, and how many nitrite rows were analyzable.
qc_summary <- tibble(
  metric = c(
    "longcycler_fastas_found",
    "longcycler_fastas_matched_to_overview",
    "genes_screened",
    "extraction_rows",
    "isolates_with_presence_rows",
    "nitrite_joined_rows",
    "binary_nitrite_rows",
    "nitrite_positive_rows",
    "nitrite_negative_rows",
    "unmatched_nitrite_rows",
    "unmatched_nitrite_uti_event_rows"
  ),
  value = c(
    nrow(longcycler_fastas),
    nrow(matched_for_screening),
    length(expected_gene_names),
    nrow(extraction_tbl),
    nrow(presence_wide),
    nrow(nitrite_joined),
    nrow(analysis_rows),
    sum(analysis_rows$nitrite_binary == "positive"),
    sum(analysis_rows$nitrite_binary == "negative"),
    sum(nitrite_joined$nitrite_join_status == "no_matching_nitrite_row"),
    sum(nitrite_joined$nitrite_join_status == "no_matching_nitrite_row" & nitrite_joined$Event_type == "UTI_event")
  )
)
write_csv(qc_summary, file.path(output_dir, "nitrate_nitrite_qc_summary.csv"), na = "")

# Console summary for a quick sanity check without opening the CSVs.
message_line("")
message_line("Done.")
message_line("Longcycler FASTAs found: ", nrow(longcycler_fastas))
message_line("Longcycler FASTAs matched to overview: ", nrow(matched_for_screening))
message_line("Presence matrix rows: ", nrow(presence_wide))
message_line("Binary nitrite rows: ", nrow(analysis_rows))
message_line("Nitrite positive/negative: ",
             sum(analysis_rows$nitrite_binary == "positive"), "/",
             sum(analysis_rows$nitrite_binary == "negative"))
message_line("Main Fisher/BH output: ", file.path(output_dir, "gene_vs_nitrite_fisher_bh.csv"))

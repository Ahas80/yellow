#!/usr/bin/env Rscript
# ==============================================================================
# 31_audit_fasta_usage.R
# ==============================================================================
#
# GOAL:
#   Audit whether every FASTA assembly on disk is represented in the key
#   metadata and downstream WGS/VF outputs used by the rUTIs pipeline.
#
# METHOD:
#   1. Recursively scan DIR_FASTAS for FASTA files.
#   2. Normalise FASTA names to isolate IDs using conservative filename rules.
#   3. Compare the FASTA scan against assembly_metadata.csv, assemblies.list,
#      WGS QC, MLST, and VF outputs when those files are available.
#   4. Report missing, unused, stale, duplicated, and assembler-conflict cases.
#
# INPUTS:
#   - 00_config.R
#   - DIR_FASTAS
#   - FILE_METADATA
#   - FILE_ASSEMBLIES
#   - results/wgs/qc_summary.csv or equivalent QC locations
#   - FILE_MLST_CANONICAL
#   - FILE_VF_HITS
#   - FILE_VF_PA
#
# OUTPUTS:
#   - results/qc/fasta_usage_audit.csv
#   - results/qc/fasta_usage_audit_summary.txt
#   - results/qc/fasta_duplicate_or_assembler_conflicts.csv
#   - results/qc/fasta_batch_usage_summary.csv
#
# USAGE:
#   Rscript 31_audit_fasta_usage.R
#
# Biological/statistical purpose:
#   Ensures that each sequenced E. coli assembly, including newly added batches
#   4-6, is visible to QC, MLST, VF screening, and related downstream analyses.
#   This helps prevent silent sample loss and flags multiple assemblies per
#   isolate before biological interpretation.
# ==============================================================================

# ------------------------------------------------------------------------------
# 1. Configuration and libraries
# ------------------------------------------------------------------------------
source("00_config.R")

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
  library(tibble)
  library(tidyr)
  library(purrr)
  library(fs)
})

# Set to FALSE if you want the audit to write reports and warn only, even when
# unused FASTAs or stale metadata references are found.
STRICT_FAIL <- TRUE

ensure_dir(DIR_QC)

OUT_AUDIT <- file.path(DIR_QC, "fasta_usage_audit.csv")
OUT_SUMMARY <- file.path(DIR_QC, "fasta_usage_audit_summary.txt")
OUT_CONFLICTS <- file.path(DIR_QC, "fasta_duplicate_or_assembler_conflicts.csv")
OUT_BATCH <- file.path(DIR_QC, "fasta_batch_usage_summary.csv")
OUT_GFF_PANAROO <- file.path(DIR_QC, "gff_panaroo_audit.csv")

log_audit <- function(...) {
  message(format(Sys.time(), "[%H:%M:%S] "), sprintf(...))
}

`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

# ------------------------------------------------------------------------------
# 2. Helper functions
# ------------------------------------------------------------------------------

# Strip FASTA extensions, assembler labels, and common assembly suffixes, then
# try to recover the laboratory isolate ID. The function is intentionally
# conservative: if no isolate-like token is found, it returns the cleaned stem.
normalise_isolate_id <- function(x) {
  if (length(x) == 0) return(character())

  vapply(x, function(xx) {
    if (is.na(xx) || !nzchar(trimws(as.character(xx)))) return(NA_character_)

    stem <- basename(as.character(xx))
    stem <- str_remove(stem, regex("\\.(fasta|fna|fa)(\\.gz)?$", ignore_case = TRUE))
    stem <- str_replace_all(stem, "\\s+", "_")
    stem <- str_replace_all(stem, "\\.+", "_")
    stem <- str_replace_all(stem, "_+", "_")

    # Remove repeated trailing assembly/assembler descriptors while preserving
    # true isolate suffixes such as "-1" in 24110099601-1.
    suffix_pattern <- regex(
      "([_-](flye|long[-_]?cycler|unicycler|assembly|assemblies|contigs?|scaffolds?|consensus|polished|medaka|racon|final|draft|complete|circulari[sz]ed|chromosome|plasmid))+$",
      ignore_case = TRUE
    )
    for (i in seq_len(4)) {
      new_stem <- str_remove(stem, suffix_pattern)
      if (identical(new_stem, stem)) break
      stem <- new_stem
    }
    stem <- str_replace_all(stem, "_+$", "")

    # Common naming pattern:
    #   PR0048_barcode52_2521C199301-1-flye.fasta
    barcode_match <- str_match(
      stem,
      regex("barcode[0-9]+[_-]+([A-Za-z0-9]+(?:-[0-9]+)?)", ignore_case = TRUE)
    )[, 2]
    if (!is.na(barcode_match)) return(str_to_upper(barcode_match))

    # Isolate IDs with the batch/year-style C token, e.g. 2521C199301-1.
    c_match <- str_extract(stem, regex("\\d{4}C\\d{4,}(?:-\\d+)?", ignore_case = TRUE))
    if (!is.na(c_match)) return(str_to_upper(c_match))

    # Pure numeric isolate IDs, e.g. 24110099601-1.
    numeric_match <- str_extract(stem, "\\d{8,}(?:-\\d+)?")
    if (!is.na(numeric_match)) return(str_to_upper(numeric_match))

    # Short fallback for unusual barcode entries such as PR0049_barcode91_35.
    short_barcode_match <- str_match(
      stem,
      regex("barcode[0-9]+[_-]+([^_/-]+(?:-[0-9]+)?)", ignore_case = TRUE)
    )[, 2]
    if (!is.na(short_barcode_match)) return(str_to_upper(short_barcode_match))

    str_to_upper(stem)
  }, character(1), USE.NAMES = FALSE)
}

fasta_extension <- function(x) {
  ext <- str_extract(basename(x), regex("\\.(fasta|fna|fa)(\\.gz)?$", ignore_case = TRUE))
  tolower(ext)
}

fasta_stem <- function(x) {
  str_remove(basename(x), regex("\\.(fasta|fna|fa)(\\.gz)?$", ignore_case = TRUE))
}

detect_probable_batch <- function(x) {
  if (length(x) == 0) return(character())

  vapply(x, function(xx) {
    s <- str_to_lower(as.character(xx))
    if (str_detect(s, "batch\\s*4\\s*[-_]\\s*6|batch4-6|batch4_6")) return("Batch 4-6")
    m <- str_match(s, "batch\\s*[_-]?\\s*([0-9]+)")[, 2]
    if (!is.na(m)) return(paste0("Batch ", m))
    "unknown"
  }, character(1), USE.NAMES = FALSE)
}

detect_assembler <- function(x) {
  if (length(x) == 0) return(character())

  case_when(
    str_detect(x, regex("long[-_]?cycler", ignore_case = TRUE)) ~ "LongCycler",
    str_detect(x, regex("flye", ignore_case = TRUE)) ~ "Flye",
    str_detect(x, regex("unicycler", ignore_case = TRUE)) ~ "Unicycler",
    TRUE ~ "unknown"
  )
}

normalise_path_safe <- function(x) {
  if (length(x) == 0) return(character())
  vapply(x, function(xx) {
    if (is.na(xx) || !nzchar(trimws(as.character(xx)))) return(NA_character_)
    normalizePath(as.character(xx), winslash = "/", mustWork = FALSE)
  }, character(1), USE.NAMES = FALSE)
}

is_absolute_path <- function(x) {
  str_detect(x, "^/|^[A-Za-z]:")
}

path_variants <- function(x) {
  x <- as.character(x)
  x <- x[!is.na(x) & nzchar(trimws(x))]
  if (length(x) == 0) return(character())

  variants <- purrr::map(x, function(xx) {
    xx <- trimws(xx)
    if (is_absolute_path(xx)) {
      xx
    } else {
      c(xx, file.path(DIR_ROOT, xx), file.path(DIR_FASTAS, xx))
    }
  }) %>%
    unlist(use.names = FALSE) %>%
    unique()

  normalise_path_safe(variants)
}

first_nonempty_value <- function(dat, cols) {
  cols <- cols[cols %in% names(dat)]
  values <- rep(NA_character_, nrow(dat))
  used_cols <- rep(NA_character_, nrow(dat))

  for (col in cols) {
    candidate <- trimws(as.character(dat[[col]]))
    candidate[candidate == ""] <- NA_character_
    take <- is.na(values) & !is.na(candidate)
    values[take] <- candidate[take]
    used_cols[take] <- col
  }

  tibble(value = values, column = used_cols)
}

build_source_evidence <- function(dat, source_name, id_cols = character(), file_cols = character()) {
  if (is.null(dat) || nrow(dat) == 0) {
    return(list(
      source = source_name,
      n_rows = 0L,
      id_cols = character(),
      file_cols = character(),
      id_values = character(),
      file_values = character(),
      ids = character(),
      file_names = character(),
      path_variants = character(),
      row_index = integer(),
      row_id = character(),
      row_file_value = character(),
      row_id_value = character(),
      has_id_identifiers = FALSE,
      has_file_identifiers = FALSE
    ))
  }

  id_cols <- id_cols[id_cols %in% names(dat)]
  file_cols <- file_cols[file_cols %in% names(dat)]

  id_first <- first_nonempty_value(dat, id_cols)
  file_first <- first_nonempty_value(dat, file_cols)

  id_values <- id_first$value[!is.na(id_first$value)]
  file_values <- file_first$value[!is.na(file_first$value)]

  ids <- unique(c(normalise_isolate_id(id_values), normalise_isolate_id(file_values)))
  ids <- ids[!is.na(ids) & nzchar(ids)]

  file_names <- unique(basename(file_values))
  file_names <- file_names[!is.na(file_names) & nzchar(file_names)]

  list(
    source = source_name,
    n_rows = nrow(dat),
    id_cols = id_cols,
    file_cols = file_cols,
    id_values = id_values,
    file_values = file_values,
    ids = ids,
    file_names = file_names,
    path_variants = unique(path_variants(file_values)),
    row_index = seq_len(nrow(dat)),
    row_id = normalise_isolate_id(ifelse(!is.na(id_first$value), id_first$value, file_first$value)),
    row_file_value = file_first$value,
    row_id_value = id_first$value,
    has_id_identifiers = length(id_values) > 0,
    has_file_identifiers = length(file_values) > 0
  )
}

mark_scan_presence <- function(scan_tbl, evidence, allow_id_fallback = FALSE) {
  if (is.null(evidence)) return(rep(NA, nrow(scan_tbl)))

  file_match <- rep(FALSE, nrow(scan_tbl))
  id_match <- rep(FALSE, nrow(scan_tbl))

  if (isTRUE(evidence$has_file_identifiers)) {
    file_match <- scan_tbl$full_path_norm %in% evidence$path_variants |
      scan_tbl$file_name %in% evidence$file_names
  }

  if (isTRUE(evidence$has_id_identifiers)) {
    id_match <- scan_tbl$normalised_isolate_id %in% evidence$ids
  }

  if (isTRUE(evidence$has_file_identifiers) && isTRUE(allow_id_fallback)) {
    return(file_match | id_match)
  }

  if (isTRUE(evidence$has_file_identifiers)) {
    return(file_match)
  }

  if (isTRUE(evidence$has_id_identifiers)) {
    return(id_match)
  }

  rep(NA, nrow(scan_tbl))
}

source_rows_without_fasta <- function(evidence, scan_tbl) {
  if (is.null(evidence) || evidence$n_rows == 0) {
    return(tibble(source = character(), source_row = integer(), source_value = character(),
                  normalised_isolate_id = character()))
  }

  scan_ids <- unique(scan_tbl$normalised_isolate_id)
  scan_names <- unique(scan_tbl$file_name)
  scan_paths <- unique(scan_tbl$full_path_norm)

  row_file_paths <- purrr::map(evidence$row_file_value, path_variants)
  row_file_match <- purrr::map_lgl(seq_along(row_file_paths), function(i) {
    vals <- row_file_paths[[i]]
    file_value <- evidence$row_file_value[i]
    if (length(vals) > 0 && any(vals %in% scan_paths)) return(TRUE)
    if (!is.na(file_value) && basename(file_value) %in% scan_names) return(TRUE)
    FALSE
  })

  row_id_match <- !is.na(evidence$row_id) & evidence$row_id %in% scan_ids
  matched <- row_file_match | row_id_match

  tibble(
    source = evidence$source,
    source_row = evidence$row_index,
    source_value = dplyr::coalesce(evidence$row_file_value, evidence$row_id_value),
    normalised_isolate_id = evidence$row_id
  ) %>%
    filter(!matched)
}

safe_read_csv <- function(path) {
  tryCatch(readr::read_csv(path, show_col_types = FALSE, progress = FALSE),
           error = function(e) {
             warning("Could not read CSV: ", path, " (", e$message, ")")
             NULL
           })
}

safe_read_tsv <- function(path) {
  tryCatch(readr::read_tsv(path, show_col_types = FALSE, progress = FALSE),
           error = function(e) {
             warning("Could not read TSV: ", path, " (", e$message, ")")
             NULL
           })
}

safe_read_delim_auto <- function(path) {
  if (grepl("\\.tsv$", path, ignore.case = TRUE)) {
    safe_read_tsv(path)
  } else {
    safe_read_csv(path)
  }
}

safe_read_rds <- function(path) {
  tryCatch(readRDS(path),
           error = function(e) {
             warning("Could not read RDS: ", path, " (", e$message, ")")
             NULL
           })
}

count_true <- function(x) {
  sum(x %in% TRUE, na.rm = TRUE)
}

count_false_if_available <- function(x) {
  if (all(is.na(x))) return(NA_integer_)
  sum(!(x %in% TRUE), na.rm = TRUE)
}

collapse_unique <- function(x) {
  x <- unique(as.character(x[!is.na(x) & nzchar(as.character(x))]))
  if (length(x) == 0) return(NA_character_)
  paste(x, collapse = "; ")
}

format_missing_from <- function(row) {
  checks <- c(
    metadata = row[["in_metadata"]],
    assemblies_list = row[["in_assemblies_list"]],
    wgs_qc = row[["in_wgs_qc"]],
    mlst = row[["in_mlst"]],
    vf_hits = row[["in_vf_hits"]],
    vf_pa = row[["in_vf_pa"]]
  )

  missing <- names(checks)[!is.na(checks) & !(checks %in% TRUE)]
  if (length(missing) == 0) "" else paste(missing, collapse = "; ")
}

# ------------------------------------------------------------------------------
# 3. FASTA discovery
# ------------------------------------------------------------------------------
log_audit("Scanning FASTA directory: %s", DIR_FASTAS)

fasta_scan_raw <- discover_project_fastas(DIR_FASTAS, recursive = TRUE, include_excluded = TRUE)

fasta_scan <- fasta_scan_raw %>%
  mutate(
    full_path_norm = normalise_path_safe(full_path),
    normalised_isolate_id = Isolate_ID,
    probable_batch = detect_probable_batch(full_path),
    assembler_guess = detect_assembler(full_path)
  ) %>%
  arrange(normalised_isolate_id, file_name, full_path)

log_audit("Found %d FASTA files representing %d normalised isolate IDs.",
          nrow(fasta_scan), n_distinct(fasta_scan$normalised_isolate_id))
log_audit("Canonical candidate input FASTAs: %d; excluded legacy/cache/nested FASTAs: %d.",
          sum(fasta_scan$include_in_metadata), sum(!fasta_scan$include_in_metadata))

# ------------------------------------------------------------------------------
# 4. assembly_metadata.csv
# ------------------------------------------------------------------------------
metadata_df <- NULL
metadata_evidence <- NULL
metadata_missing_rows <- tibble(source = character(), source_row = integer(),
                                source_value = character(), normalised_isolate_id = character())
metadata_duplicate_ids <- tibble(normalised_isolate_id = character(), n_rows = integer())

metadata_id_cols <- c("Isolate_ID", "isolate_ID", "isolate_id", "IsolateID", "SampleID", "Sample_ID")
metadata_file_cols <- c("full_path", "file_name", "assembly", "fasta", "path", "FILE", "filename")

if (file.exists(FILE_METADATA)) {
  log_audit("Reading metadata: %s", FILE_METADATA)
  metadata_df <- safe_read_csv(FILE_METADATA)
  metadata_evidence <- build_source_evidence(metadata_df, "assembly_metadata", metadata_id_cols, metadata_file_cols)
  metadata_missing_rows <- source_rows_without_fasta(metadata_evidence, fasta_scan)

  metadata_duplicate_ids <- tibble(normalised_isolate_id = metadata_evidence$row_id) %>%
    filter(!is.na(normalised_isolate_id)) %>%
    dplyr::count(normalised_isolate_id, name = "n_rows") %>%
    filter(n_rows > 1) %>%
    arrange(desc(n_rows), normalised_isolate_id)

  if ("Batch" %in% names(metadata_df)) {
    metadata_batch_by_id <- tibble(
      normalised_isolate_id = metadata_evidence$row_id,
      metadata_batch = as.character(metadata_df$Batch)
    ) %>%
      filter(!is.na(normalised_isolate_id), !is.na(metadata_batch), nzchar(metadata_batch)) %>%
      group_by(normalised_isolate_id) %>%
      summarise(metadata_batch = first(metadata_batch), .groups = "drop")

    metadata_batch_by_file <- tibble(
      file_name = basename(metadata_evidence$row_file_value),
      metadata_batch_file = as.character(metadata_df$Batch)
    ) %>%
      filter(!is.na(file_name), nzchar(file_name),
             !is.na(metadata_batch_file), nzchar(metadata_batch_file)) %>%
      group_by(file_name) %>%
      summarise(metadata_batch_file = first(metadata_batch_file), .groups = "drop")

    fasta_scan <- fasta_scan %>%
      left_join(metadata_batch_by_id, by = "normalised_isolate_id") %>%
      left_join(metadata_batch_by_file, by = "file_name") %>%
      mutate(
        probable_batch = if_else(
          probable_batch == "unknown" & !is.na(coalesce(metadata_batch, metadata_batch_file)),
          paste0("Batch ", coalesce(metadata_batch, metadata_batch_file)),
          probable_batch
        )
      ) %>%
      select(-metadata_batch, -metadata_batch_file)
  }
} else {
  warning("Metadata file not available: ", FILE_METADATA)
}

# ------------------------------------------------------------------------------
# 5. assemblies.list
# ------------------------------------------------------------------------------
assemblies_evidence <- NULL
assemblies_list_df <- tibble(
  listed_path = character(),
  listed_path_exists = logical(),
  duplicate_entry = logical(),
  normalised_isolate_id = character()
)

if (file.exists(FILE_ASSEMBLIES)) {
  log_audit("Reading assemblies list: %s", FILE_ASSEMBLIES)
  list_lines <- readr::read_lines(FILE_ASSEMBLIES, progress = FALSE)
  list_lines <- trimws(list_lines)
  list_lines <- list_lines[nzchar(list_lines)]

  assemblies_list_df <- tibble(listed_path = list_lines) %>%
    mutate(
      file_name = basename(listed_path),
      normalised_isolate_id = normalise_isolate_id(file_name),
      listed_path_exists = purrr::map_lgl(listed_path, ~ any(file.exists(c(
        .x,
        file.path(DIR_ROOT, .x),
        file.path(DIR_FASTAS, .x)
      )))),
      duplicate_entry = duplicated(listed_path) | duplicated(listed_path, fromLast = TRUE)
    )

  assemblies_evidence <- build_source_evidence(
    assemblies_list_df,
    "assemblies.list",
    id_cols = character(),
    file_cols = c("listed_path")
  )
} else {
  warning("Assemblies list not available: ", FILE_ASSEMBLIES)
}

# ------------------------------------------------------------------------------
# 6. WGS QC output
# ------------------------------------------------------------------------------
qc_candidates <- unique(c(
  file.path(DIR_WGS, "qc_summary.csv"),
  file.path(DIR_WGS_QC, "qc_summary.csv"),
  file.path(DIR_RESULTS, "qc", "qc_summary.csv")
))
qc_file <- qc_candidates[file.exists(qc_candidates)][1] %||% NA_character_
qc_df <- NULL
qc_evidence <- NULL
qc_missing_rows <- tibble(source = character(), source_row = integer(),
                          source_value = character(), normalised_isolate_id = character())

qc_id_cols <- c("Isolate_ID", "isolate_ID", "isolate_id", "IsolateID", "SampleID", "Sample_ID")
qc_file_cols <- c("full_path", "file_name", "assembly", "fasta", "path", "FILE", "filename")

if (!is.na(qc_file)) {
  log_audit("Reading WGS QC: %s", qc_file)
  qc_df <- safe_read_csv(qc_file)
  qc_evidence <- build_source_evidence(qc_df, "wgs_qc", qc_id_cols, qc_file_cols)
  qc_missing_rows <- source_rows_without_fasta(qc_evidence, fasta_scan)
} else {
  warning("No WGS QC summary found in expected locations.")
}

# ------------------------------------------------------------------------------
# 7. MLST output
# ------------------------------------------------------------------------------
mlst_df <- NULL
mlst_evidence <- NULL
mlst_missing_rows <- tibble(source = character(), source_row = integer(),
                            source_value = character(), normalised_isolate_id = character())
mlst_duplicate_ids <- tibble(normalised_isolate_id = character(), n_records = integer())

mlst_id_cols <- c("Isolate_ID", "isolate_ID", "isolate_id", "IsolateID", "SampleID", "Sample_ID")
mlst_file_cols <- c("full_path", "file_name", "FILE", "file", "filename", "assembly", "fasta", "path", "#FILE")

if (file.exists(FILE_MLST_CANONICAL)) {
  log_audit("Reading active provider-preferred MLST: %s", FILE_MLST_CANONICAL)
  mlst_df <- safe_read_delim_auto(FILE_MLST_CANONICAL)
  mlst_evidence <- build_source_evidence(mlst_df, "mlst", mlst_id_cols, mlst_file_cols)
  mlst_missing_rows <- source_rows_without_fasta(mlst_evidence, fasta_scan)

  mlst_duplicate_ids <- tibble(normalised_isolate_id = mlst_evidence$row_id) %>%
    filter(!is.na(normalised_isolate_id)) %>%
    dplyr::count(normalised_isolate_id, name = "n_records") %>%
    filter(n_records > 1) %>%
    arrange(desc(n_records), normalised_isolate_id)
} else {
  warning("Active provider-preferred MLST output not available: ", FILE_MLST_CANONICAL)
}

# ------------------------------------------------------------------------------
# 8. VF outputs
# ------------------------------------------------------------------------------
vf_hits_df <- NULL
vf_hits_evidence <- NULL
vf_hits_missing_rows <- tibble(source = character(), source_row = integer(),
                               source_value = character(), normalised_isolate_id = character())
vf_hits_note <- "not available"

vf_id_cols <- c("Isolate_ID", "isolate_ID", "isolate_id", "IsolateID", "SampleID", "Sample_ID")
vf_file_cols <- c("full_path", "file_name", "FILE", "file", "filename", "assembly", "fasta", "path", "#FILE")

if (file.exists(FILE_VF_HITS)) {
  log_audit("Reading VF hits: %s", FILE_VF_HITS)
  vf_hits_df <- safe_read_rds(FILE_VF_HITS)
  if (!is.null(vf_hits_df)) {
    vf_hits_df <- as_tibble(vf_hits_df)
    vf_hits_evidence <- build_source_evidence(vf_hits_df, "vf_hits", vf_id_cols, vf_file_cols)
    vf_hits_missing_rows <- source_rows_without_fasta(vf_hits_evidence, fasta_scan)

    if (!vf_hits_evidence$has_file_identifiers && !vf_hits_evidence$has_id_identifiers &&
        all(c("Participant_id", "tp_lab") %in% names(vf_hits_df))) {
      vf_hits_note <- "participant-timepoint only; direct FASTA tracing not possible"
    } else if (!vf_hits_evidence$has_file_identifiers && !vf_hits_evidence$has_id_identifiers) {
      vf_hits_note <- "no isolate or file column detected"
    } else {
      vf_hits_note <- "FASTA-level or isolate-level identifiers detected"
    }
  }
} else {
  warning("VF hits output not available: ", FILE_VF_HITS)
}

vf_pa_df <- NULL
vf_pa_evidence <- NULL
vf_pa_missing_rows <- tibble(source = character(), source_row = integer(),
                             source_value = character(), normalised_isolate_id = character())
vf_pa_note <- "not available"
vf_pa_gene_cols <- character()

if (file.exists(FILE_VF_PA)) {
  log_audit("Reading VF presence/absence: %s", FILE_VF_PA)
  vf_pa_df <- safe_read_csv(FILE_VF_PA)
  if (!is.null(vf_pa_df)) {
    vf_pa_evidence <- build_source_evidence(vf_pa_df, "vf_pa", vf_id_cols, vf_file_cols)
    vf_pa_missing_rows <- source_rows_without_fasta(vf_pa_evidence, fasta_scan)
    vf_pa_gene_cols <- canonical_vf_gene_cols(names(vf_pa_df), vf_pa_file = FILE_VF_PA, required = FALSE)

    if (!vf_pa_evidence$has_file_identifiers && !vf_pa_evidence$has_id_identifiers) {
      vf_pa_note <- "participant-timepoint level; cannot alone prove every FASTA was screened"
    } else {
      vf_pa_note <- "FASTA-level or isolate-level identifiers detected"
    }
  }
} else {
  warning("VF presence/absence output not available: ", FILE_VF_PA)
}

# ------------------------------------------------------------------------------
# 9. Duplicate and assembler-conflict checks
# ------------------------------------------------------------------------------
id_conflicts <- fasta_scan %>%
  group_by(normalised_isolate_id) %>%
  summarise(
    conflict_key = first(normalised_isolate_id),
    n_fastas = n(),
    n_assemblers = n_distinct(assembler_guess),
    n_directories = n_distinct(directory),
    assemblers = collapse_unique(assembler_guess),
    directories = collapse_unique(directory),
    file_names = collapse_unique(file_name),
    full_paths = collapse_unique(full_path),
    has_flye = any(assembler_guess == "Flye"),
    has_longcycler = any(assembler_guess == "LongCycler"),
    conflict_type = paste(c(
      if (n() > 1) "MULTIPLE_FASTAS_SAME_ISOLATE_ID",
      if (any(assembler_guess == "Flye") && any(assembler_guess == "LongCycler")) "FLYE_AND_LONGCYCLER_PRESENT",
      if (n_distinct(directory) > 1) "SAME_ISOLATE_IN_MULTIPLE_DIRECTORIES"
    ), collapse = "; "),
    .groups = "drop"
  ) %>%
  filter(nzchar(conflict_type)) %>%
  mutate(conflict_level = "normalised_isolate_id")

file_name_conflicts <- fasta_scan %>%
  group_by(file_name) %>%
  summarise(
    conflict_key = first(file_name),
    normalised_isolate_id = collapse_unique(normalised_isolate_id),
    n_fastas = n(),
    n_assemblers = n_distinct(assembler_guess),
    n_directories = n_distinct(directory),
    assemblers = collapse_unique(assembler_guess),
    directories = collapse_unique(directory),
    file_names = collapse_unique(file_name),
    full_paths = collapse_unique(full_path),
    has_flye = any(assembler_guess == "Flye"),
    has_longcycler = any(assembler_guess == "LongCycler"),
    conflict_type = if_else(n_distinct(full_path) > 1, "EXACT_FILE_NAME_IN_MULTIPLE_FOLDERS", ""),
    .groups = "drop"
  ) %>%
  filter(nzchar(conflict_type)) %>%
  mutate(conflict_level = "file_name")

conflicts <- bind_rows(
  id_conflicts %>%
    select(conflict_level, conflict_key, normalised_isolate_id, everything()),
  file_name_conflicts %>%
    select(conflict_level, conflict_key, normalised_isolate_id, everything())
) %>%
  distinct(conflict_level, conflict_key, conflict_type, .keep_all = TRUE) %>%
  arrange(conflict_level, conflict_key)

conflicted_ids <- unique(id_conflicts$normalised_isolate_id)

# ------------------------------------------------------------------------------
# 10. Main per-FASTA audit table
# ------------------------------------------------------------------------------
audit_tbl <- fasta_scan
audit_tbl$in_metadata <- mark_scan_presence(audit_tbl, metadata_evidence, allow_id_fallback = TRUE)
audit_tbl$in_assemblies_list <- mark_scan_presence(audit_tbl, assemblies_evidence)
audit_tbl$in_wgs_qc <- mark_scan_presence(audit_tbl, qc_evidence)
audit_tbl$in_mlst <- mark_scan_presence(audit_tbl, mlst_evidence)
audit_tbl$in_vf_hits <- mark_scan_presence(audit_tbl, vf_hits_evidence)
audit_tbl$in_vf_pa <- mark_scan_presence(audit_tbl, vf_pa_evidence)

available_downstream_cols <- c("in_assemblies_list", "in_wgs_qc", "in_mlst", "in_vf_hits", "in_vf_pa")
available_downstream_cols <- available_downstream_cols[
  !vapply(audit_tbl[available_downstream_cols], function(x) all(is.na(x)), logical(1))
]

if (length(available_downstream_cols) > 0) {
  downstream_true <- as.data.frame(lapply(audit_tbl[available_downstream_cols], function(x) x %in% TRUE))
  audit_tbl$any_downstream_evidence <- rowSums(downstream_true) > 0
  audit_tbl$n_available_downstream_checks <- length(available_downstream_cols)
  audit_tbl$n_missing_available_downstream <- rowSums(!downstream_true)
} else {
  audit_tbl$any_downstream_evidence <- FALSE
  audit_tbl$n_available_downstream_checks <- 0L
  audit_tbl$n_missing_available_downstream <- 0L
}

audit_tbl$missing_from <- purrr::pmap_chr(
  audit_tbl %>% select(in_metadata, in_assemblies_list, in_wgs_qc, in_mlst, in_vf_hits, in_vf_pa),
  function(in_metadata, in_assemblies_list, in_wgs_qc, in_mlst, in_vf_hits, in_vf_pa) {
    checks <- c(
      metadata = in_metadata,
      assemblies_list = in_assemblies_list,
      wgs_qc = in_wgs_qc,
      mlst = in_mlst,
      vf_hits = in_vf_hits,
      vf_pa = in_vf_pa
    )
    missing <- names(checks)[!is.na(checks) & !(checks %in% TRUE)]
    if (length(missing) == 0) "" else paste(missing, collapse = "; ")
  }
)

audit_tbl <- audit_tbl %>%
  mutate(
    extra_fasta_category = case_when(
      !include_in_metadata ~ fasta_class,
      in_metadata %in% FALSE & normalised_isolate_id %in% conflicted_ids ~ "naming_normalisation_conflict",
      in_metadata %in% FALSE & !any_downstream_evidence ~ "unexpected_unlinked",
      in_metadata %in% FALSE & any_downstream_evidence ~ "downstream_without_metadata",
      normalised_isolate_id %in% conflicted_ids &
        str_detect(conflicts$conflict_type[match(normalised_isolate_id, conflicts$normalised_isolate_id)] %||% "", "FLYE_AND_LONGCYCLER_PRESENT") ~ "legitimate_assembler_alternative",
      TRUE ~ "expected_metadata_linked"
    ),
    audit_status = case_when(
      !include_in_metadata ~ "EXCLUDED_FROM_METADATA_UNIVERSE",
      normalised_isolate_id %in% conflicted_ids ~ "DUPLICATE_OR_CONFLICT",
      in_metadata %in% FALSE & any_downstream_evidence ~ "DOWNSTREAM_WITHOUT_METADATA",
      in_metadata %in% TRUE & !any_downstream_evidence ~ "METADATA_ONLY_NO_DOWNSTREAM",
      in_metadata %in% FALSE & !any_downstream_evidence ~ "NOT_USED_ANYWHERE",
      in_metadata %in% TRUE & any_downstream_evidence & n_missing_available_downstream == 0 ~ "OK_USED",
      in_metadata %in% TRUE & any_downstream_evidence & n_missing_available_downstream > 0 ~ "PARTIAL_USED",
      TRUE ~ "UNKNOWN_CHECK_MANUALLY"
    )
  ) %>%
  select(
    normalised_isolate_id,
    file_name,
    full_path,
    probable_batch,
    assembler_guess,
    in_metadata,
    in_assemblies_list,
    in_wgs_qc,
    in_mlst,
    in_vf_hits,
    in_vf_pa,
    any_downstream_evidence,
    missing_from,
    audit_status,
    extra_fasta_category,
    include_in_metadata,
    fasta_class,
    exclusion_reason,
    file_stem,
    extension,
    directory
  ) %>%
  arrange(audit_status, probable_batch, normalised_isolate_id, file_name)

# ------------------------------------------------------------------------------
# 11. Batch-specific summary
# ------------------------------------------------------------------------------
batch_summary <- audit_tbl %>%
  group_by(probable_batch) %>%
  summarise(
    n_fastas_on_disk = n(),
    n_in_metadata = count_true(in_metadata),
    n_in_assemblies_list = count_true(in_assemblies_list),
    n_in_qc = count_true(in_wgs_qc),
    n_in_mlst = count_true(in_mlst),
    n_in_vf_hits = count_true(in_vf_hits),
    n_missing_metadata = count_false_if_available(in_metadata),
    n_missing_qc = count_false_if_available(in_wgs_qc),
    n_missing_mlst = count_false_if_available(in_mlst),
    n_missing_vf_hits = count_false_if_available(in_vf_hits),
    .groups = "drop"
  ) %>%
  arrange(probable_batch)

# ------------------------------------------------------------------------------
# 12. GFF/Panaroo QA summary
# ------------------------------------------------------------------------------
gff_inventory <- build_assembly_gff_inventory(write_outputs = TRUE)
gff_summary <- gff_inventory$summary
summary_value <- function(metric, default = NA_character_) {
  hit <- gff_summary$value[gff_summary$metric == metric]
  if (length(hit) == 0 || is.na(hit[1])) default else hit[1]
}

regen_gffs_file <- file.path(DIR_WGS, "pan", "regenerate_missing_gffs_summary.csv")
panaroo_status_file <- file.path(DIR_WGS, "pan", "panaroo_staleness_report.txt")

regen_gffs <- if (file.exists(regen_gffs_file)) safe_read_csv(regen_gffs_file) else NULL
panaroo_status_lines <- if (file.exists(panaroo_status_file)) {
  grep("^Status:", readLines(panaroo_status_file, warn = FALSE), value = TRUE)
} else {
  character()
}

metadata_gff_expected <- as.integer(summary_value("metadata_linked_fastas"))
metadata_gff_available <- as.integer(summary_value("metadata_linked_gffs_available"))
metadata_gff_missing <- as.integer(summary_value("metadata_linked_gffs_missing_warning_only"))
expected_missing_fastas <- as.integer(summary_value("expected_metadata_rows_missing_fasta"))
unexpected_candidate_fastas <- as.integer(summary_value("unexpected_unlinked_candidate_fastas"))
new_candidate_fastas <- as.integer(summary_value("new_candidate_fastas_since_metadata_scan"))
removed_candidate_fastas <- as.integer(summary_value("removed_candidate_fastas_since_metadata_scan"))
gff_expected <- as.integer(summary_value("panaroo_eligible_assemblies"))
gff_available <- as.integer(summary_value("panaroo_gffs_available"))
missing_gff_count <- as.integer(summary_value("panaroo_gffs_missing_required"))
prokka_failure_count <- if (!is.null(regen_gffs) && "gff_created" %in% names(regen_gffs)) {
  sum(!(regen_gffs$gff_created %in% TRUE), na.rm = TRUE)
} else {
  NA_integer_
}
panaroo_status_text <- if (length(panaroo_status_lines) > 0) {
  paste(panaroo_status_lines, collapse = " | ")
} else {
  "not available"
}
panaroo_green <- any(grepl("^Status: GREEN", panaroo_status_lines))
gff_complete <- !is.na(gff_expected) && gff_expected > 0 &&
  !is.na(gff_available) && !is.na(missing_gff_count) &&
  gff_available == gff_expected && missing_gff_count == 0
metadata_inventory_stale <- isTRUE(gff_inventory$metadata_stale)
qc_selection_stale <- isTRUE(gff_inventory$qc_stale)
gff_panaroo_passed <- gff_complete && panaroo_green && !metadata_inventory_stale && !qc_selection_stale

fmt_qa <- function(x) {
  if (length(x) == 0 || is.na(x)) "not available" else as.character(x)
}

gff_panaroo_audit <- tibble(
  check = c(
    "metadata_linked_fastas",
    "metadata_linked_gffs_available",
    "metadata_linked_gffs_missing_warning_only",
    "expected_metadata_rows_missing_fasta",
    "unexpected_unlinked_candidate_fastas",
    "new_candidate_fastas_since_metadata_scan",
    "removed_candidate_fastas_since_metadata_scan",
    "panaroo_eligible_assemblies",
    "panaroo_gffs_available",
    "panaroo_gffs_missing_required",
    "prokka_regeneration_failures",
    "metadata_inventory_stale",
    "qc_selection_stale",
    "panaroo_status",
    "gff_panaroo_passed"
  ),
  value = c(
    fmt_qa(metadata_gff_expected),
    fmt_qa(metadata_gff_available),
    fmt_qa(metadata_gff_missing),
    fmt_qa(expected_missing_fastas),
    fmt_qa(unexpected_candidate_fastas),
    fmt_qa(new_candidate_fastas),
    fmt_qa(removed_candidate_fastas),
    fmt_qa(gff_expected),
    fmt_qa(gff_available),
    fmt_qa(missing_gff_count),
    fmt_qa(prokka_failure_count),
    ifelse(metadata_inventory_stale, "YES", "NO"),
    ifelse(qc_selection_stale, "YES", "NO"),
    panaroo_status_text,
    ifelse(gff_panaroo_passed, "YES", "NO")
  )
)

# ------------------------------------------------------------------------------
# 13. Write machine-readable outputs
# ------------------------------------------------------------------------------
readr::write_csv(audit_tbl, OUT_AUDIT)
readr::write_csv(conflicts, OUT_CONFLICTS)
readr::write_csv(batch_summary, OUT_BATCH)
readr::write_csv(gff_panaroo_audit, OUT_GFF_PANAROO)

# ------------------------------------------------------------------------------
# 14. Human-readable summary report
# ------------------------------------------------------------------------------
unused_tbl <- audit_tbl %>% filter(audit_status == "NOT_USED_ANYWHERE")
unused_candidate_tbl <- audit_tbl %>% filter(include_in_metadata, audit_status == "NOT_USED_ANYWHERE")
missing_metadata_tbl <- audit_tbl %>% filter(in_metadata %in% FALSE)
missing_candidate_metadata_tbl <- audit_tbl %>% filter(include_in_metadata, in_metadata %in% FALSE)
missing_mlst_tbl <- audit_tbl %>% filter(in_mlst %in% FALSE)
missing_vf_hits_tbl <- audit_tbl %>% filter(in_vf_hits %in% FALSE)
problem_tbl <- audit_tbl %>%
  filter(audit_status != "OK_USED" | missing_from != "") %>%
  arrange(audit_status, desc(!is.na(missing_from) & nzchar(missing_from)), probable_batch, normalised_isolate_id) %>%
  slice_head(n = 50)

availability <- tibble(
  check = c("assembly_metadata", "assemblies.list", "wgs_qc", "mlst", "vf_hits", "vf_pa"),
  path = c(FILE_METADATA, FILE_ASSEMBLIES, qc_file %||% NA_character_, FILE_MLST_CANONICAL, FILE_VF_HITS, FILE_VF_PA),
  available = c(
    file.exists(FILE_METADATA),
    file.exists(FILE_ASSEMBLIES),
    !is.na(qc_file),
    file.exists(FILE_MLST_CANONICAL),
    file.exists(FILE_VF_HITS),
    file.exists(FILE_VF_PA)
  ),
  note = c(
    if (is.null(metadata_evidence)) "not available" else paste("identifier columns:", paste(c(metadata_evidence$id_cols, metadata_evidence$file_cols), collapse = ", ")),
    if (is.null(assemblies_evidence)) "not available" else paste0("entries: ", nrow(assemblies_list_df), "; stale paths: ", sum(!assemblies_list_df$listed_path_exists)),
    if (is.null(qc_evidence)) "not available" else paste("identifier columns:", paste(c(qc_evidence$id_cols, qc_evidence$file_cols), collapse = ", ")),
    if (is.null(mlst_evidence)) "not available" else paste("identifier columns:", paste(c(mlst_evidence$id_cols, mlst_evidence$file_cols), collapse = ", ")),
    vf_hits_note,
    vf_pa_note
  )
)

batch_report <- capture.output(print(batch_summary, n = Inf))
gff_panaroo_report <- capture.output(print(gff_panaroo_audit, n = Inf, width = Inf))
status_report <- capture.output(print(audit_tbl %>% dplyr::count(audit_status, sort = TRUE), n = Inf))
availability_report <- capture.output(print(availability, n = Inf, width = Inf))
problem_report <- capture.output(print(
  problem_tbl %>%
    select(normalised_isolate_id, file_name, probable_batch, assembler_guess, missing_from, audit_status),
  n = 50,
  width = Inf
))
metadata_missing_report <- capture.output(print(metadata_missing_rows %>% slice_head(n = 50), n = 50, width = Inf))
conflict_report <- capture.output(print(conflicts %>% slice_head(n = 50), n = 50, width = Inf))

report_lines <- c(
  "FASTA Usage Audit Summary",
  "=========================",
  paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  paste0("FASTA directory scanned: ", DIR_FASTAS),
  "",
  "Headline counts",
  "---------------",
  paste0("Total FASTAs found: ", nrow(audit_tbl)),
  paste0("Canonical candidate FASTAs: ", sum(audit_tbl$include_in_metadata)),
  paste0("Excluded legacy/cache/nested FASTAs: ", sum(!audit_tbl$include_in_metadata)),
  paste0("Total unique normalised isolate IDs found: ", n_distinct(audit_tbl$normalised_isolate_id)),
  paste0("Represented in assembly_metadata.csv: ", count_true(audit_tbl$in_metadata)),
  paste0("Represented in assemblies.list: ", count_true(audit_tbl$in_assemblies_list)),
  paste0("Represented in WGS QC: ", count_true(audit_tbl$in_wgs_qc)),
  paste0("Represented in MLST: ", count_true(audit_tbl$in_mlst)),
  paste0("Represented in VF hits: ", count_true(audit_tbl$in_vf_hits)),
  paste0("Represented in VF PA where detectable: ", count_true(audit_tbl$in_vf_pa)),
  paste0("Completely unused FASTAs: ", nrow(unused_tbl)),
  paste0("Completely unused candidate FASTAs: ", nrow(unused_candidate_tbl)),
  paste0("FASTAs missing from metadata: ", nrow(missing_metadata_tbl)),
  paste0("Candidate FASTAs missing from metadata: ", nrow(missing_candidate_metadata_tbl)),
  paste0("FASTAs missing from MLST: ", nrow(missing_mlst_tbl)),
  paste0("FASTAs missing from VF hits: ", nrow(missing_vf_hits_tbl)),
  paste0("Duplicate/conflict rows: ", nrow(conflicts)),
  paste0("Metadata rows pointing to missing FASTAs: ", nrow(metadata_missing_rows)),
  paste0("Duplicated isolate IDs in metadata: ", nrow(metadata_duplicate_ids)),
  paste0("Duplicated MLST records per isolate: ", nrow(mlst_duplicate_ids)),
  "",
  "GFF/Panaroo QA",
  "--------------",
  gff_panaroo_report,
  "",
  "Check availability",
  "------------------",
  availability_report,
  "",
  "Audit status counts",
  "-------------------",
  status_report,
  "",
  "Batch-level summary",
  "-------------------",
  batch_report,
  "",
  "Top 50 problematic FASTAs",
  "-------------------------",
  problem_report,
  "",
  "Metadata rows pointing to missing FASTAs (top 50)",
  "------------------------------------------------",
  metadata_missing_report,
  "",
  "Duplicate or assembler conflicts (top 50)",
  "-----------------------------------------",
  conflict_report,
  "",
  "Interpretation",
  "--------------",
  "These files are on disk but appear unused: rows with audit_status NOT_USED_ANYWHERE.",
  "These metadata rows point to missing FASTAs: rows listed above under metadata missing FASTAs.",
  "These isolates have multiple assemblies and may need explicit selection logic: rows in fasta_duplicate_or_assembler_conflicts.csv.",
  "Flye versus LongCycler/long-cycler duplicates are reported only; this script does not delete, choose, or rename assemblies.",
  "VF presence/absence is participant-timepoint level, so it may not be sufficient to prove FASTA-level usage unless vf_hits_all.rds contains isolate/file identifiers.",
  "GFF verification is two-tier: missing GFFs for all metadata-linked FASTAs are warnings, but missing Panaroo-eligible GFFs block final analysis.",
  "",
  "Output files",
  "------------",
  paste0("Main audit table: ", OUT_AUDIT),
  paste0("Batch summary: ", OUT_BATCH),
  paste0("Duplicate/conflict table: ", OUT_CONFLICTS),
  paste0("GFF/Panaroo audit table: ", OUT_GFF_PANAROO),
  paste0("All-metadata GFF inventory: ", gff_inventory$paths$all_metadata),
  paste0("Panaroo GFF inventory: ", gff_inventory$paths$panaroo),
  paste0("GFF inventory summary: ", gff_inventory$paths$summary),
  paste0("Text report: ", OUT_SUMMARY)
)

readr::write_lines(report_lines, OUT_SUMMARY)

# ------------------------------------------------------------------------------
# 15. Console summary and exit status
# ------------------------------------------------------------------------------
serious_issue_count <- nrow(unused_candidate_tbl) + nrow(metadata_missing_rows) +
  ifelse(gff_panaroo_passed, 0L, 1L)
audit_passed <- serious_issue_count == 0

cat("\nFASTA usage audit complete\n")
cat("--------------------------\n")
cat("Total FASTAs: ", nrow(audit_tbl), "\n", sep = "")
cat("Canonical candidate FASTAs: ", sum(audit_tbl$include_in_metadata), "\n", sep = "")
cat("Excluded legacy/cache/nested FASTAs: ", sum(!audit_tbl$include_in_metadata), "\n", sep = "")
cat("Unused FASTAs: ", nrow(unused_tbl), "\n", sep = "")
cat("Unused candidate FASTAs: ", nrow(unused_candidate_tbl), "\n", sep = "")
cat("Missing metadata: ", nrow(missing_metadata_tbl), "\n", sep = "")
cat("Candidate FASTAs missing metadata: ", nrow(missing_candidate_metadata_tbl), "\n", sep = "")
cat("Metadata rows pointing to missing FASTAs: ", nrow(metadata_missing_rows), "\n", sep = "")
cat("Duplicate/conflict rows: ", nrow(conflicts), "\n", sep = "")
cat("Metadata-linked GFFs available: ", fmt_qa(metadata_gff_available), "/", fmt_qa(metadata_gff_expected), "\n", sep = "")
cat("Metadata-linked missing GFFs (warning only): ", fmt_qa(metadata_gff_missing), "\n", sep = "")
cat("Panaroo GFFs available: ", fmt_qa(gff_available), "/", fmt_qa(gff_expected), "\n", sep = "")
cat("Panaroo missing GFFs (blocking): ", fmt_qa(missing_gff_count), "\n", sep = "")
cat("Prokka regeneration failures: ", fmt_qa(prokka_failure_count), "\n", sep = "")
cat("Metadata inventory stale: ", ifelse(metadata_inventory_stale, "YES", "NO"), "\n", sep = "")
cat("QC selection stale: ", ifelse(qc_selection_stale, "YES", "NO"), "\n", sep = "")
cat("Panaroo status: ", panaroo_status_text, "\n", sep = "")
cat("Output audit: ", OUT_AUDIT, "\n", sep = "")
cat("Output summary: ", OUT_SUMMARY, "\n", sep = "")
cat("Output conflicts: ", OUT_CONFLICTS, "\n", sep = "")
cat("Output batch summary: ", OUT_BATCH, "\n", sep = "")
cat("Output GFF/Panaroo audit: ", OUT_GFF_PANAROO, "\n", sep = "")
cat("Audit passed: ", ifelse(audit_passed, "YES", "NO"), "\n", sep = "")

append_denominator_summary(
  audit_tbl,
  "31_audit_fasta_usage.R",
  "fasta_usage_audit",
  "assembly",
  OUT_AUDIT,
  "Uses the same canonical FASTA discovery helper as script 00; excluded FASTAs are classified rather than mixed into the metadata universe"
)

if (!audit_passed && isTRUE(STRICT_FAIL)) {
  quit(save = "no", status = 1)
}

quit(save = "no", status = 0)

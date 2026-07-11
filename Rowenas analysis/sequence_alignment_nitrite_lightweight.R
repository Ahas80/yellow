#!/usr/bin/env Rscript

# Lightweight one-file workflow:
# 1) screen Yellow routine FASTAs for local nitrate/nitrite reference genes
# 2) summarize gene presence profiles by isolate
# 3) join SPSS nitrite results
# 4) write barcode/nitrite trend CSVs

required_packages <- c("Biostrings", "dplyr", "haven", "readr", "readxl", "stringr", "tidyr")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0) {
  stop("Missing required R packages: ", paste(missing_packages, collapse = ", "))
}

suppressPackageStartupMessages({
  library(Biostrings)
  library(dplyr)
  library(haven)
  library(readr)
  library(stringr)
  library(tidyr)
})

args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(name, default = NULL) {
  prefix <- paste0(name, "=")
  hit <- args[startsWith(args, prefix)]
  if (length(hit) == 0) return(default)
  sub(paste0("^", prefix), "", hit[[1]])
}

get_script_dir <- function() {
  cmd <- commandArgs(trailingOnly = FALSE)
  file_arg <- cmd[startsWith(cmd, "--file=")]
  if (length(file_arg) > 0) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg[[1]]), winslash = "/", mustWork = FALSE)))
  }
  normalizePath(getwd(), winslash = "/", mustWork = TRUE)
}

normalise_pid <- function(x) {
  x_chr <- trimws(as.character(x))
  is_intish <- str_detect(x_chr, "^[0-9]+([.]0+)?$")
  x_num <- suppressWarnings(as.integer(as.numeric(x_chr)))
  ifelse(is_intish & !is.na(x_num), as.character(x_num), x_chr)
}

as_flag <- function(x) {
  x %in% c(TRUE, "TRUE", "true", "True", "1", 1L, 1)
}

assert_required_cols <- function(df, required, context) {
  missing <- setdiff(required, names(df))
  if (length(missing) > 0) {
    stop(context, " is missing required columns: ", paste(missing, collapse = ", "))
  }
}

drop_missing_paths <- function(x) {
  unique(x[!is.na(x) & nzchar(x)])
}

existing_dirs <- function(paths) {
  paths <- drop_missing_paths(paths)
  paths <- paths[dir.exists(paths)]
  unique(normalizePath(paths, winslash = "/", mustWork = TRUE))
}

parent_dirs <- function(path, max_depth = 8) {
  if (!dir.exists(path)) return(character())
  path <- normalizePath(path, winslash = "/", mustWork = TRUE)
  out <- path
  for (i in seq_len(max_depth)) {
    parent <- dirname(tail(out, 1))
    if (identical(parent, tail(out, 1))) break
    out <- c(out, parent)
  }
  unique(out)
}

project_root_score <- function(path) {
  score <- 0L
  if (file.exists(file.path(path, "assembly_metadata.csv"))) score <- score + 4L
  if (dir.exists(file.path(path, "results")) &&
      file.exists(file.path(path, "results", "assembly_metadata.csv"))) {
    score <- score + 2L
  }
  if (file.exists(file.path(path, "OVERVIEW E.coli batch 1-6 - CLEAN_RC.xlsx"))) {
    score <- score + 2L
  }
  if (length(list.files(path, pattern = "[.]Rproj$", full.names = TRUE)) > 0) {
    score <- score + 2L
  }
  if (dir.exists(file.path(path, "Rowenas analysis"))) score <- score + 1L
  if (file.exists(file.path(path, "YELLOW dataset stage Rowena Studie meetmomenten long format_1.sav"))) {
    score <- score + 1L
  }
  if (file.exists(file.path(path, "SPSS_nitrite_results"))) score <- score + 1L
  score
}

find_project_root <- function(script_dir) {
  env_root <- Sys.getenv("YELLOW_PROJECT_ROOT", unset = "")
  windows_roots <- c(
    "D:/Batch 1-6 E.coli FULL sequence YELLOW Study",
    "D:/Batch 1-6 E.coli FULL sequence YELLOW study",
    "D:/Batch 1-6 Ecoli FULL sequence YELLOW Study",
    "D:/Batch 1-6 Ecoli FULL sequence YELLOW study"
  )
  seeds <- existing_dirs(c(env_root, windows_roots, script_dir, getwd()))
  candidates <- unique(unlist(lapply(seeds, parent_dirs), use.names = FALSE))
  candidates <- candidates[dir.exists(candidates)]
  scores <- vapply(candidates, project_root_score, integer(1))
  if (any(scores > 0)) {
    return(candidates[order(scores, decreasing = TRUE)][1])
  }
  stop(
    "Could not find the YELLOW project root. Put this script somewhere inside ",
    "the project folder, or set Sys.setenv(YELLOW_PROJECT_ROOT = 'path/to/project')."
  )
}

find_project_file <- function(project_root, file_name, label, extra_dirs = character()) {
  direct_dirs <- existing_dirs(c(
    extra_dirs,
    project_root,
    file.path(project_root, "Rowenas analysis"),
    file.path(project_root, "metadata"),
    file.path(project_root, "results"),
    file.path(project_root, "data")
  ))
  direct_hits <- file.path(direct_dirs, file_name)
  direct_hits <- direct_hits[file.exists(direct_hits)]
  if (length(direct_hits) > 0) {
    return(normalizePath(direct_hits[[1]], winslash = "/", mustWork = TRUE))
  }

  recursive_hits <- list.files(
    project_root,
    pattern = paste0("^", file_name, "$"),
    recursive = TRUE,
    full.names = TRUE,
    ignore.case = FALSE
  )
  recursive_hits <- recursive_hits[file.exists(recursive_hits)]
  if (length(recursive_hits) > 0) {
    return(normalizePath(recursive_hits[[1]], winslash = "/", mustWork = TRUE))
  }

  stop("Cannot find ", label, " named '", file_name, "' under: ", project_root)
}

find_optional_project_file <- function(project_root, file_name, extra_dirs = character()) {
  direct_dirs <- existing_dirs(c(
    extra_dirs,
    project_root,
    file.path(project_root, "Rowenas analysis"),
    file.path(project_root, "metadata"),
    file.path(project_root, "results"),
    file.path(project_root, "data")
  ))
  direct_hits <- file.path(direct_dirs, file_name)
  direct_hits <- direct_hits[file.exists(direct_hits)]
  if (length(direct_hits) > 0) {
    return(normalizePath(direct_hits[[1]], winslash = "/", mustWork = TRUE))
  }

  recursive_hits <- list.files(
    project_root,
    pattern = paste0("^", file_name, "$"),
    recursive = TRUE,
    full.names = TRUE,
    ignore.case = FALSE
  )
  recursive_hits <- recursive_hits[file.exists(recursive_hits)]
  if (length(recursive_hits) > 0) {
    return(normalizePath(recursive_hits[[1]], winslash = "/", mustWork = TRUE))
  }

  NA_character_
}

find_reference_dir <- function(project_root, script_dir, expected_files) {
  direct_dirs <- existing_dirs(c(
    script_dir,
    getwd(),
    file.path(project_root, "Rowenas analysis"),
    file.path(project_root, "reference_genes"),
    file.path(project_root, "references"),
    file.path(project_root, "metadata"),
    file.path(project_root, "data")
  ))
  score_dir <- function(path) {
    sum(file.exists(file.path(path, expected_files)))
  }
  direct_scores <- vapply(direct_dirs, score_dir, integer(1))
  if (length(direct_scores) > 0 && max(direct_scores) > 0) {
    return(direct_dirs[which.max(direct_scores)])
  }

  recursive_hits <- list.files(
    project_root,
    pattern = "[.]fasta$",
    recursive = TRUE,
    full.names = TRUE,
    ignore.case = TRUE
  )
  recursive_hits <- recursive_hits[basename(recursive_hits) %in% expected_files]
  if (length(recursive_hits) == 0) {
    stop("Cannot find the selected nitrate/nitrite reference FASTAs under: ", project_root)
  }

  hit_dirs <- unique(dirname(recursive_hits))
  hit_scores <- vapply(hit_dirs, score_dir, integer(1))
  hit_dirs[which.max(hit_scores)]
}

normalise_overview_timepoint <- function(x) {
  x_chr <- trimws(as.character(x))
  x_chr[x_chr %in% c("", "NA", "NaN")] <- NA_character_
  numeric_tp <- suppressWarnings(as.integer(as.numeric(x_chr)))
  case_when(
    is.na(x_chr) ~ NA_character_,
    str_detect(x_chr, regex("^T[0-9]+$", ignore_case = TRUE)) ~ str_to_upper(x_chr),
    str_detect(x_chr, regex("^UTI[-_ ]*[0-9]+$", ignore_case = TRUE)) ~
      str_replace_all(str_to_upper(x_chr), "[_ ]+", "-"),
    str_detect(x_chr, "^[0-9]+([.]0+)?$") ~ paste0("T", numeric_tp - 1L),
    TRUE ~ x_chr
  )
}

read_overview_metadata <- function(overview_file) {
  if (!requireNamespace("readxl", quietly = TRUE)) {
    stop(
      "Package 'readxl' is required to read the overview workbook. ",
      "Install it once with install.packages('readxl')."
    )
  }

  sheets <- readxl::excel_sheets(overview_file)
  sheet <- if ("Batches overzicht" %in% sheets) "Batches overzicht" else sheets[[1]]
  raw <- readxl::read_excel(overview_file, sheet = sheet, skip = 3)

  rename_if_present <- function(df, old, new) {
    if (old %in% names(df) && !new %in% names(df)) {
      names(df)[names(df) == old] <- new
    }
    df
  }

  raw <- raw %>%
    rename_if_present("Isolaat id", "Isolate_ID") %>%
    rename_if_present("Meetmoment", "Timepoint") %>%
    rename_if_present("Organisme", "Organism") %>%
    rename_if_present("Beoord", "Clinical_Beoord") %>%
    rename_if_present("Kiemgetal", "Clinical_CFU_Count") %>%
    rename_if_present("Urine opvang methode", "Urine_collection_method")

  for (col in c(
    "Batch", "Participant_id", "Isolate_ID", "Timepoint", "Organism",
    "Clinical_Beoord", "Clinical_CFU_Count", "Urine_collection_method"
  )) {
    if (!col %in% names(raw)) raw[[col]] <- NA_character_
  }

  raw %>%
    transmute(
      Batch = suppressWarnings(as.integer(.data$Batch)),
      Participant_id = trimws(as.character(.data$Participant_id)),
      Isolate_ID = trimws(as.character(.data$Isolate_ID)),
      Timepoint = trimws(as.character(.data$Timepoint)),
      tp_lab = normalise_overview_timepoint(.data$Timepoint),
      Event_type = ifelse(str_detect(.data$tp_lab, regex("^UTI", ignore_case = TRUE)),
                          "UTI_event", "routine"),
      UTI_Label = ifelse(.data$Event_type == "UTI_event", .data$tp_lab, NA_character_),
      Clinical_Beoord = as.character(.data$Clinical_Beoord),
      Clinical_CFU_Count = as.character(.data$Clinical_CFU_Count),
      Clinical_Organism = as.character(.data$Organism),
      Urine_collection_method = as.character(.data$Urine_collection_method),
      overview_row_number = row_number()
    ) %>%
    filter(!is.na(.data$Isolate_ID), .data$Isolate_ID != "")
}

discover_project_assembly_fastas <- function(project_root, expected_reference_files) {
  preferred_roots <- existing_dirs(c(
    file.path(project_root, "ont-yellow-routine-fastas"),
    file.path(project_root, "assemblies"),
    file.path(project_root, "fastas"),
    file.path(project_root, "data", "assemblies"),
    file.path(project_root, "data", "fastas")
  ))
  search_roots <- if (length(preferred_roots) > 0) preferred_roots else project_root

  fasta_files <- unique(unlist(lapply(search_roots, function(root) {
    list.files(
      root,
      pattern = "[.](fasta|fa|fna)$",
      recursive = TRUE,
      full.names = TRUE,
      ignore.case = TRUE
    )
  }), use.names = FALSE))

  if (length(fasta_files) == 0) {
    stop("No assembly FASTA/FA/FNA files found under: ", paste(search_roots, collapse = ", "))
  }

  fasta_files <- normalizePath(fasta_files, winslash = "/", mustWork = TRUE)
  project_norm <- normalizePath(project_root, winslash = "/", mustWork = TRUE)
  rel_path <- ifelse(
    startsWith(fasta_files, paste0(project_norm, "/")),
    substring(fasta_files, nchar(project_norm) + 2L),
    fasta_files
  )

  tibble(full_path = fasta_files, relative_path = rel_path, file_name = basename(fasta_files)) %>%
    filter(
      !.data$file_name %in% expected_reference_files,
      .data$file_name != "reference_genes_combined.fasta",
      str_detect(.data$file_name, regex("^PR[0-9]+_barcode[0-9]+_", ignore_case = TRUE)),
      str_detect(.data$file_name, regex("(flye|longcycler)", ignore_case = TRUE)),
      !str_detect(.data$relative_path, regex("(^|/)outputs/", ignore_case = TRUE)),
      !str_detect(.data$relative_path, regex("(^|/)results/prokka/", ignore_case = TRUE)),
      !str_detect(.data$relative_path, regex("(^|/)not-assemblies/", ignore_case = TRUE)),
      !str_detect(.data$relative_path, regex("(^|/)legacy/", ignore_case = TRUE))
    ) %>%
    mutate(
      extension = tolower(tools::file_ext(.data$file_name)),
      Assembly_Base_ID = str_remove(.data$file_name, regex("[.](fasta|fa|fna)$", ignore_case = TRUE)),
      Assembly_ID = paste0(.data$Assembly_Base_ID, "__", .data$extension),
      Assembler = case_when(
        str_detect(.data$Assembly_Base_ID, regex("longcycler", ignore_case = TRUE)) ~ "longcycler",
        str_detect(.data$Assembly_Base_ID, regex("flye", ignore_case = TRUE)) ~ "flye",
        TRUE ~ NA_character_
      ),
      Isolate_ID_raw = .data$Assembly_Base_ID %>%
        str_remove(regex("^PR[0-9]+_barcode[0-9]+_", ignore_case = TRUE)) %>%
        str_remove(regex("[-_](flye|longcycler)$", ignore_case = TRUE))
    )
}

summarise_fasta_metrics <- function(path) {
  if (is.na(path) || !file.exists(path)) {
    return(tibble(num_contigs = NA_integer_, total_bases = NA_real_))
  }
  tryCatch({
    seqs <- readDNAStringSet(path)
    tibble(num_contigs = length(seqs), total_bases = as.numeric(sum(width(seqs))))
  }, error = function(e) tibble(num_contigs = NA_integer_, total_bases = NA_real_))
}

build_assembly_metadata_from_overview_and_fastas <- function(
  project_root,
  overview_file,
  expected_reference_files
) {
  if (is.na(overview_file) || !file.exists(overview_file)) {
    stop(
      "Cannot build metadata by FASTA discovery because the overview workbook ",
      "was not found."
    )
  }

  message("Building analysis table from overview workbook plus project FASTA discovery.")
  overview <- read_overview_metadata(overview_file)
  fasta_tbl <- discover_project_assembly_fastas(project_root, expected_reference_files)

  suffix_lookup <- overview %>%
    mutate(file_id_without_suffix = str_remove(.data$Isolate_ID, "-[0-9]+$")) %>%
    filter(.data$file_id_without_suffix != .data$Isolate_ID) %>%
    group_by(.data$file_id_without_suffix) %>%
    summarise(
      matched_isolate_id = first(.data$Isolate_ID),
      n_expected_ids = n_distinct(.data$Isolate_ID),
      .groups = "drop"
    ) %>%
    filter(.data$n_expected_ids == 1)

  fasta_tbl <- fasta_tbl %>%
    left_join(suffix_lookup, by = c("Isolate_ID_raw" = "file_id_without_suffix")) %>%
    mutate(
      Isolate_ID = case_when(
        .data$Isolate_ID_raw %in% overview$Isolate_ID ~ .data$Isolate_ID_raw,
        !str_detect(.data$Isolate_ID_raw, "-[0-9]+$") & !is.na(.data$matched_isolate_id) ~
          .data$matched_isolate_id,
        TRUE ~ .data$Isolate_ID_raw
      )
    )

  metrics <- bind_rows(lapply(fasta_tbl$full_path, summarise_fasta_metrics))
  fasta_tbl <- bind_cols(fasta_tbl, metrics)

  metadata <- overview %>%
    left_join(
      fasta_tbl %>%
        select(
          Isolate_ID, Assembly_ID, Assembly_Base_ID, file_name, full_path,
          relative_path, Assembler, Isolate_ID_raw, num_contigs, total_bases
        ),
      by = "Isolate_ID"
    ) %>%
    mutate(
      file_exists = !is.na(.data$full_path) & file.exists(.data$full_path),
      usable_fasta = .data$file_exists,
      fasta_path = .data$full_path,
      assembler = .data$Assembler,
      found = .data$file_exists,
      analysis_include_primary = .data$file_exists,
      metadata_source_status = ifelse(.data$file_exists, "overview_plus_discovered_fasta", "overview_no_fasta"),
      Assembly_ID = ifelse(is.na(.data$Assembly_ID), paste0("missing__", .data$Isolate_ID), .data$Assembly_ID)
    )

  curation_file <- find_optional_project_file(
    project_root,
    "manual_sample_curation.csv",
    extra_dirs = c(file.path(project_root, "data"), file.path(project_root, "metadata"))
  )
  if (!is.na(curation_file)) {
    curation <- read_csv(curation_file, show_col_types = FALSE) %>%
      transmute(
        Isolate_ID = trimws(as.character(.data$Isolate_ID)),
        exclude_primary = as_flag(.data$exclude_primary),
        exclude_from_genomics_expected = as_flag(.data$exclude_from_genomics_expected),
        analysis_exclusion_reason = as.character(.data$exclude_reason),
        duplicate_role = as.character(.data$duplicate_role),
        duplicate_of_participant_id = as.character(.data$duplicate_of_participant_id),
        duplicate_of_tp_lab = as.character(.data$duplicate_of_tp_lab),
        allow_secondary_duplicate_qc = as_flag(.data$allow_secondary_duplicate_qc),
        duplicate_use_note = as.character(.data$duplicate_use_note),
        manual_curation_note = as.character(.data$curation_note),
        manual_curation_source = as.character(.data$curation_source)
      ) %>%
      filter(!is.na(.data$Isolate_ID), .data$Isolate_ID != "")

    metadata <- metadata %>%
      left_join(curation, by = "Isolate_ID") %>%
      mutate(
        exclude_primary = coalesce(.data$exclude_primary, FALSE),
        exclude_from_genomics_expected = coalesce(.data$exclude_from_genomics_expected, FALSE),
        analysis_include_primary = .data$file_exists & !.data$exclude_primary,
        analysis_exclusion_reason = ifelse(
          .data$exclude_primary,
          .data$analysis_exclusion_reason,
          NA_character_
        ),
        genomics_expected_include = !.data$exclude_from_genomics_expected,
        duplicate_role = coalesce(.data$duplicate_role, "not_duplicate"),
        allow_secondary_duplicate_qc = coalesce(.data$allow_secondary_duplicate_qc, FALSE),
        manual_curation_applied = .data$exclude_primary | .data$exclude_from_genomics_expected
      )
  } else {
    metadata <- metadata %>%
      mutate(
        analysis_exclusion_reason = NA_character_,
        duplicate_role = "not_duplicate",
        duplicate_of_participant_id = NA_character_,
        duplicate_of_tp_lab = NA_character_,
        allow_secondary_duplicate_qc = FALSE,
        duplicate_use_note = NA_character_,
        genomics_expected_include = TRUE,
        manual_curation_applied = FALSE,
        manual_curation_note = NA_character_,
        manual_curation_source = NA_character_
      )
  }

  unmatched_fastas <- fasta_tbl %>%
    filter(!.data$Isolate_ID %in% overview$Isolate_ID)
  if (nrow(unmatched_fastas) > 0) {
    warning(
      nrow(unmatched_fastas),
      " discovered FASTA file(s) did not match an overview Isolate_ID and will be excluded."
    )
  }

  metadata
}

pick_named_column <- function(df, candidates) {
  exact_hit <- candidates[candidates %in% names(df)]
  if (length(exact_hit) > 0) return(exact_hit[[1]])

  normalise_name <- function(x) str_replace_all(str_to_lower(x), "[^a-z0-9]+", "")
  names_norm <- normalise_name(names(df))
  cand_norm <- normalise_name(candidates)
  hit_index <- match(cand_norm, names_norm, nomatch = 0L)
  hit_index <- hit_index[hit_index > 0L]
  if (length(hit_index) > 0) return(names(df)[hit_index[[1]]])

  NA_character_
}

read_nitrite_source <- function(path) {
  ext <- tolower(tools::file_ext(path))
  raw <- if (ext == "sav") {
    haven::read_sav(path)
  } else {
    readRDS(path)
  }

  if (!is.data.frame(raw)) {
    stop("Nitrite source is not a data frame/tibble: ", path)
  }

  pid_col <- pick_named_column(raw, c("participant_id", "Participant_id", "participantid"))
  tp_col <- pick_named_column(raw, c("timepoint", "Timepoint", "meetmoment", "tp_lab"))
  nitrite_col <- pick_named_column(raw, c(
    "urinestick_nitriet", "urinestick_nitrite", "nitriet", "nitrite"
  ))

  missing <- c(
    if (is.na(pid_col)) "participant_id" else character(),
    if (is.na(tp_col)) "timepoint" else character(),
    if (is.na(nitrite_col)) "urinestick_nitriet" else character()
  )
  if (length(missing) > 0) {
    stop(
      "Nitrite source is missing required field(s): ",
      paste(missing, collapse = ", "),
      ". Available columns are: ",
      paste(names(raw), collapse = ", ")
    )
  }

  raw %>%
    transmute(
      participant_id = .data[[pid_col]],
      timepoint = .data[[tp_col]],
      urinestick_nitriet = .data[[nitrite_col]]
    )
}

safe_median <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0) return(NA_real_)
  median(x)
}

summarise_nitrite <- function(df, ...) {
  df %>%
    group_by(...) %>%
    summarise(
      total_rows = n(),
      matched_binary_n = sum(.data$nitrite_binary %in% c("positive", "negative"), na.rm = TRUE),
      positive_n = sum(.data$nitrite_binary == "positive", na.rm = TRUE),
      negative_n = sum(.data$nitrite_binary == "negative", na.rm = TRUE),
      missing_or_unmatched_n = .data$total_rows - .data$matched_binary_n,
      positive_pct_among_matched = ifelse(
        .data$matched_binary_n > 0,
        round(100 * .data$positive_n / .data$matched_binary_n, 1),
        NA_real_
      ),
      distinct_participants = n_distinct(.data$Participant_id),
      distinct_isolates = n_distinct(.data$Isolate_ID),
      .groups = "drop"
    ) %>%
    arrange(desc(.data$matched_binary_n), desc(.data$positive_pct_among_matched))
}

max_isolates <- suppressWarnings(as.integer(get_arg("--max-isolates", NA_character_)))
identity_threshold <- as.numeric(get_arg("--min-pident", "80"))
coverage_threshold <- as.numeric(get_arg("--min-qcov", "80"))

if (is.na(identity_threshold) || is.na(coverage_threshold)) {
  stop("--min-pident and --min-qcov must be numeric.")
}

expected_reference_files <- paste0(
  c("CsgD", "NarG", "NarH", "NarI", "NarJ", "NarL",
    "NarP", "NarQ", "NarV", "NarX", "NarY", "NarZ"),
  ".fasta"
)

script_dir <- get_script_dir()
project_root <- find_project_root(script_dir)
overview_metadata_file <- find_project_file(
  project_root,
  "OVERVIEW E.coli batch 1-6 - CLEAN_RC.xlsx",
  "overview metadata workbook",
  extra_dirs = c(script_dir, getwd())
)
nitrite_sav_file <- find_optional_project_file(
  project_root,
  "YELLOW dataset stage Rowena Studie meetmomenten long format_1.sav",
  extra_dirs = c(script_dir, getwd())
)
nitrite_rds_file <- find_optional_project_file(
  project_root,
  "SPSS_nitrite_results",
  extra_dirs = c(script_dir, getwd())
)
nitrite_file <- if (!is.na(nitrite_sav_file)) nitrite_sav_file else nitrite_rds_file
if (is.na(nitrite_file)) {
  stop(
    "Cannot find nitrite source. Expected the SPSS file named ",
    "'YELLOW dataset stage Rowena Studie meetmomenten long format_1.sav' ",
    "under the project folder."
  )
}
reference_dir <- find_reference_dir(project_root, script_dir, expected_reference_files)
out_dir <- file.path(script_dir, "outputs", "nitrate_blast")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(nitrite_file)) {
  stop("Cannot find nitrite source file: ", nitrite_file)
}
if (Sys.which("blastn") == "") {
  stop("Cannot find blastn on PATH. Install BLAST+ or add blastn to PATH.")
}

message("Using project root: ", project_root)
message("Using overview workbook for isolate metadata: ", overview_metadata_file)
message("Discovering assembly FASTA paths under project root.")
message("Using reference FASTA folder: ", reference_dir)
message("Using nitrite source: ", nitrite_file)
message("Writing outputs to: ", out_dir)

message("Loading overview metadata and discovering FASTA files...")
assembly_metadata <- build_assembly_metadata_from_overview_and_fastas(
  project_root,
  overview_metadata_file,
  expected_reference_files
)

required_metadata_cols <- c(
  "Assembly_ID", "Isolate_ID", "Participant_id", "tp_lab", "Timepoint",
  "Event_type", "Assembler", "file_name", "full_path", "file_exists",
  "usable_fasta", "analysis_include_primary", "num_contigs", "total_bases",
  "UTI_Label", "Clinical_Beoord", "Clinical_CFU_Count",
  "Urine_collection_method"
)
assert_required_cols(assembly_metadata, required_metadata_cols, overview_metadata_file)

selected_fastas <- assembly_metadata %>%
  mutate(
    extension = tolower(tools::file_ext(.data$file_name)),
    assembler_rank = case_when(
      .data$Assembler == "longcycler" & .data$extension == "fasta" ~ 1L,
      .data$Assembler == "flye" & .data$extension == "fasta" ~ 2L,
      .data$Assembler == "longcycler" ~ 3L,
      .data$Assembler == "flye" ~ 4L,
      TRUE ~ 5L
    ),
    num_contigs_rank = ifelse(is.na(.data$num_contigs), Inf, .data$num_contigs),
    total_bases_rank = ifelse(is.na(.data$total_bases), -Inf, .data$total_bases),
    sequencing_run = str_extract(.data$file_name, "^PR[0-9]+"),
    sequencing_barcode = str_match(.data$file_name, "barcode([0-9]+)")[, 2]
  ) %>%
  filter(
    as_flag(.data$file_exists),
    as_flag(.data$usable_fasta),
    as_flag(.data$analysis_include_primary),
    .data$extension %in% c("fasta", "fna", "fa")
  ) %>%
  arrange(
    .data$Isolate_ID,
    .data$assembler_rank,
    .data$num_contigs_rank,
    desc(.data$total_bases_rank)
  ) %>%
  group_by(.data$Isolate_ID) %>%
  slice(1) %>%
  ungroup()

if (!is.na(max_isolates)) {
  selected_fastas <- selected_fastas %>% slice_head(n = max_isolates)
}

if (nrow(selected_fastas) == 0) {
  stop("No usable FASTA assemblies selected from assembly metadata.")
}

message("Selected ", nrow(selected_fastas), " FASTA assemblies for screening.")
write_csv(selected_fastas, file.path(out_dir, "selected_yellow_routine_fastas.csv"))

reference_files <- file.path(reference_dir, expected_reference_files)
missing_reference_files <- reference_files[!file.exists(reference_files)]
if (length(missing_reference_files) > 0) {
  stop(
    "Missing expected reference FASTA files: ",
    paste(basename(missing_reference_files), collapse = ", "),
    ". Looked in: ",
    reference_dir
  )
}

message("Combining ", length(reference_files), " reference genes...")
reference_sets <- lapply(reference_files, function(path) {
  x <- readDNAStringSet(path)
  if (length(x) == 0) stop("Reference FASTA is empty: ", path)
  x <- x[1]
  names(x) <- tools::file_path_sans_ext(basename(path))
  x
})
references <- do.call(c, reference_sets)

reference_lengths <- tibble(
  gene = names(references),
  gene_lower = str_to_lower(names(references)),
  qlen_reference = as.integer(width(references))
)

reference_fasta <- file.path(out_dir, "reference_genes_combined.fasta")
writeXStringSet(references, reference_fasta)
write_csv(reference_lengths, file.path(out_dir, "reference_gene_lengths.csv"))

blast_outfmt_fields <- c(
  "qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
  "qstart", "qend", "sstart", "send", "evalue", "bitscore",
  "qlen", "slen"
)
blast_columns <- c(
  "gene", "contig", "pident", "alignment_length", "mismatch", "gapopen",
  "qstart", "qend", "sstart", "send", "evalue", "bitscore",
  "qlen", "slen"
)

run_blast_one <- function(row_index) {
  row <- selected_fastas[row_index, ]
  out_file <- tempfile(fileext = ".blast.tsv")
  err_file <- tempfile(fileext = ".blast.err")
  outfmt <- paste("6", paste(blast_outfmt_fields, collapse = " "))

  status <- system2(
    "blastn",
    args = c(
      "-query", shQuote(reference_fasta),
      "-subject", shQuote(row$full_path),
      "-task", "blastn",
      "-dust", "no",
      "-soft_masking", "false",
      "-evalue", "1e-10",
      "-outfmt", shQuote(outfmt),
      "-out", shQuote(out_file)
    ),
    stdout = TRUE,
    stderr = err_file
  )

  if (!is.null(attr(status, "status")) && attr(status, "status") != 0) {
    warning("blastn failed for ", row$file_name, ": ",
            paste(readLines(err_file, warn = FALSE), collapse = " "))
    return(tibble())
  }
  if (!file.exists(out_file) || file.size(out_file) == 0) {
    return(tibble())
  }

  read_tsv(out_file, col_names = blast_columns, show_col_types = FALSE) %>%
    mutate(
      Isolate_ID = row$Isolate_ID,
      Participant_id = row$Participant_id,
      tp_lab = row$tp_lab,
      Timepoint = row$Timepoint,
      Event_type = row$Event_type,
      Assembly_ID = row$Assembly_ID,
      Assembler = row$Assembler,
      file_name = row$file_name,
      full_path = row$full_path,
      sequencing_run = row$sequencing_run,
      sequencing_barcode = row$sequencing_barcode,
      UTI_Label = row$UTI_Label,
      Clinical_Beoord = row$Clinical_Beoord,
      Clinical_CFU_Count = row$Clinical_CFU_Count,
      Urine_collection_method = row$Urine_collection_method
    )
}

message("Running BLAST screens...")
raw_hits <- vector("list", nrow(selected_fastas))
for (i in seq_len(nrow(selected_fastas))) {
  if (i == 1 || i %% 25 == 0 || i == nrow(selected_fastas)) {
    message("  ", i, "/", nrow(selected_fastas), " assemblies")
  }
  raw_hits[[i]] <- run_blast_one(i)
}

raw_hits <- bind_rows(raw_hits)
if (nrow(raw_hits) == 0) {
  raw_hits <- tibble(
    gene = character(), contig = character(), pident = numeric(),
    alignment_length = integer(), mismatch = integer(), gapopen = integer(),
    qstart = integer(), qend = integer(), sstart = integer(), send = integer(),
    evalue = numeric(), bitscore = numeric(), qlen = integer(), slen = integer(),
    Isolate_ID = character(), Participant_id = character(), tp_lab = character(),
    Timepoint = character(), Event_type = character(), Assembly_ID = character(),
    Assembler = character(), file_name = character(), full_path = character(),
    sequencing_run = character(), sequencing_barcode = character(),
    UTI_Label = character(), Clinical_Beoord = character(),
    Clinical_CFU_Count = character(), Urine_collection_method = character()
  )
}

raw_hits <- raw_hits %>%
  mutate(
    gene_lower = str_to_lower(.data$gene),
    qcov = 100 * .data$alignment_length / .data$qlen,
    strand = ifelse(.data$sstart <= .data$send, "plus", "minus"),
    present = .data$pident >= identity_threshold & .data$qcov >= coverage_threshold
  )

write_csv(raw_hits, file.path(out_dir, "nitrate_gene_blast_hits_raw.csv"))

all_isolate_gene_pairs <- tidyr::crossing(
  selected_fastas %>%
    select(
      Isolate_ID, Participant_id, tp_lab, Timepoint, Event_type,
      Assembly_ID, Assembler, file_name, full_path, sequencing_run,
      sequencing_barcode, UTI_Label, Clinical_Beoord, Clinical_CFU_Count,
      Urine_collection_method
    ),
  reference_lengths
)

best_hits <- raw_hits %>%
  group_by(.data$Isolate_ID, .data$gene) %>%
  arrange(
    desc(.data$present),
    desc(.data$qcov),
    desc(.data$pident),
    desc(.data$bitscore),
    .data$evalue,
    .by_group = TRUE
  ) %>%
  slice(1) %>%
  ungroup()

presence_long <- all_isolate_gene_pairs %>%
  left_join(
    best_hits %>%
      select(
        Isolate_ID, gene, contig, pident, alignment_length, mismatch, gapopen,
        qstart, qend, sstart, send, evalue, bitscore, qlen, slen, qcov,
        strand, present
      ),
    by = c("Isolate_ID", "gene")
  ) %>%
  mutate(
    present = replace_na(.data$present, FALSE),
    pident = ifelse(.data$present, .data$pident, NA_real_),
    qcov = ifelse(.data$present, .data$qcov, NA_real_)
  )

write_csv(best_hits, file.path(out_dir, "nitrate_gene_best_hits.csv"))

absence_checklist <- all_isolate_gene_pairs %>%
  left_join(
    best_hits %>%
      select(
        Isolate_ID, gene, contig, pident, alignment_length, mismatch, gapopen,
        qstart, qend, sstart, send, evalue, bitscore, qlen, slen, qcov,
        strand, present
      ),
    by = c("Isolate_ID", "gene")
  ) %>%
  mutate(
    no_blast_hit = is.na(.data$contig),
    present = replace_na(.data$present, FALSE)
  ) %>%
  filter(!.data$present) %>%
  arrange(.data$Isolate_ID, .data$gene)

write_csv(
  absence_checklist,
  file.path(out_dir, "nitrate_gene_absence_or_below_threshold_checklist.csv")
)
write_csv(presence_long, file.path(out_dir, "nitrate_gene_presence_long.csv"))

presence_wide <- presence_long %>%
  mutate(value = as.integer(.data$present)) %>%
  select(
    Isolate_ID, Participant_id, tp_lab, Timepoint, Event_type,
    Assembly_ID, Assembler, file_name, sequencing_run, sequencing_barcode,
    UTI_Label, Clinical_Beoord, Clinical_CFU_Count, Urine_collection_method,
    gene, value
  ) %>%
  pivot_wider(names_from = gene, values_from = value,
              names_prefix = "gene_", values_fill = 0)

write_csv(presence_wide, file.path(out_dir, "nitrate_gene_presence_by_isolate.csv"))

module_summary <- presence_long %>%
  group_by(
    .data$Isolate_ID, .data$Participant_id, .data$tp_lab, .data$Timepoint,
    .data$Event_type, .data$Assembly_ID, .data$Assembler, .data$file_name,
    .data$sequencing_run, .data$sequencing_barcode, .data$UTI_Label,
    .data$Clinical_Beoord, .data$Clinical_CFU_Count,
    .data$Urine_collection_method
  ) %>%
  summarise(
    genes_present = sum(.data$present),
    genes_screened = n(),
    narA_complete = all(c("narg", "narh", "nari", "narj") %in%
                          .data$gene_lower[.data$present]),
    narZ_complete_available_refs = all(c("narz", "nary", "narv") %in%
                                         .data$gene_lower[.data$present]),
    nar_regulators_complete = all(c("narx", "narq", "narl", "narp") %in%
                                    .data$gene_lower[.data$present]),
    csgD_present = "csgd" %in% .data$gene_lower[.data$present],
    nitrate_profile_01 = paste0(as.integer(.data$present), collapse = ""),
    genes_absent = paste(.data$gene[!.data$present], collapse = ";"),
    .groups = "drop"
  )

write_csv(module_summary, file.path(out_dir, "nitrate_operon_summary_by_isolate.csv"))

trend_ready <- module_summary %>%
  left_join(presence_wide, by = c(
    "Isolate_ID", "Participant_id", "tp_lab", "Timepoint", "Event_type",
    "Assembly_ID", "Assembler", "file_name", "sequencing_run",
    "sequencing_barcode", "UTI_Label", "Clinical_Beoord",
    "Clinical_CFU_Count", "Urine_collection_method"
  )) %>%
  mutate(
    dipstick_nitrite = NA_character_,
    nitrite_join_note = paste(
      "Dipstick nitrite is joined by Participant_id + tp_lab in this script.",
      "UTI-* event rows are reported separately when no routine nitrite row exists."
    )
  )

write_csv(trend_ready, file.path(out_dir, "nitrate_trend_ready_needs_nitrite.csv"))

gene_counts <- presence_long %>%
  group_by(.data$gene) %>%
  summarise(
    present_n = sum(.data$present),
    screened_n = n(),
    present_pct = round(100 * .data$present_n / .data$screened_n, 1),
    median_pident = round(safe_median(.data$pident), 2),
    median_qcov = round(safe_median(.data$qcov), 2),
    .groups = "drop"
  ) %>%
  arrange(.data$gene)

write_csv(gene_counts, file.path(out_dir, "nitrate_gene_presence_counts.csv"))

message("Reading nitrite source...")
nitrite_raw <- read_nitrite_source(nitrite_file)

nitrite_clean <- nitrite_raw %>%
  transmute(
    Participant_id = normalise_pid(.data$participant_id),
    spss_timepoint = as.character(.data$timepoint),
    spss_timepoint_integer = suppressWarnings(as.integer(.data$timepoint)),
    tp_lab = ifelse(
      !is.na(.data$spss_timepoint_integer),
      paste0("T", .data$spss_timepoint_integer - 1L),
      NA_character_
    ),
    urinestick_nitriet_code = suppressWarnings(as.numeric(.data$urinestick_nitriet)),
    urinestick_nitriet_label = as.character(haven::as_factor(.data$urinestick_nitriet)),
    nitrite_binary = case_when(
      .data$urinestick_nitriet_code == 1 |
        .data$urinestick_nitriet_label == "Positief" ~ "positive",
      .data$urinestick_nitriet_code == 2 |
        .data$urinestick_nitriet_label == "Negatief" ~ "negative",
      is.na(.data$urinestick_nitriet_code) |
        is.na(.data$urinestick_nitriet_label) |
        .data$urinestick_nitriet_label %in% c(
          "Not done", "Asked but unknown", "Not asked",
          "Not applicable", "Measurement failed"
        ) ~ "missing",
      TRUE ~ "missing"
    )
  )

duplicate_keys <- nitrite_clean %>%
  count(.data$Participant_id, .data$tp_lab, name = "n") %>%
  filter(.data$n > 1)

if (nrow(duplicate_keys) > 0) {
  write_csv(duplicate_keys, file.path(out_dir, "nitrite_duplicate_participant_timepoint_keys.csv"))
  stop(
    "Nitrite source has duplicate Participant_id + tp_lab keys. ",
    "See outputs/nitrate_blast/nitrite_duplicate_participant_timepoint_keys.csv"
  )
}

base_prepped <- trend_ready %>%
  select(-any_of(c("dipstick_nitrite", "nitrite_join_note"))) %>%
  mutate(
    Participant_id = normalise_pid(.data$Participant_id),
    sequencing_run = as.character(.data$sequencing_run),
    sequencing_barcode = str_pad(as.character(.data$sequencing_barcode), width = 2, pad = "0"),
    run_barcode = paste0(.data$sequencing_run, "_barcode", .data$sequencing_barcode),
    all_12_screened_genes_present = .data$genes_present == .data$genes_screened
  )

joined <- base_prepped %>%
  left_join(
    nitrite_clean,
    by = c("Participant_id", "tp_lab"),
    relationship = "many-to-one"
  ) %>%
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
  ) %>%
  relocate(
    run_barcode, nitrite_binary, nitrite_positive,
    nitrite_join_status, urinestick_nitriet_code,
    urinestick_nitriet_label, spss_timepoint,
    spss_timepoint_integer,
    .after = sequencing_barcode
  )

unmatched <- joined %>%
  filter(.data$nitrite_join_status != "matched_binary_result") %>%
  arrange(.data$Event_type, .data$Participant_id, .data$tp_lab, .data$Isolate_ID)

by_run <- summarise_nitrite(joined, .data$sequencing_run)

by_barcode_number <- summarise_nitrite(joined, .data$sequencing_barcode) %>%
  rename(barcode_number = sequencing_barcode)

by_run_barcode <- summarise_nitrite(
  joined,
  .data$sequencing_run,
  .data$sequencing_barcode,
  .data$run_barcode
)

by_profile <- summarise_nitrite(
  joined,
  .data$nitrate_profile_01,
  .data$genes_absent,
  .data$genes_present,
  .data$genes_screened,
  .data$all_12_screened_genes_present,
  .data$narA_complete,
  .data$narZ_complete_available_refs,
  .data$nar_regulators_complete
)

join_qc <- tibble(
  metric = c(
    "selected_isolates",
    "reference_genes",
    "presence_long_rows",
    "base_rows",
    "nitrite_source_rows",
    "nitrite_unique_participant_timepoint_keys",
    "joined_rows",
    "matched_binary_rows",
    "positive_rows",
    "negative_rows",
    "unmatched_or_missing_rows",
    "unmatched_uti_event_rows",
    "unmatched_routine_rows"
  ),
  value = c(
    nrow(selected_fastas),
    nrow(reference_lengths),
    nrow(presence_long),
    nrow(trend_ready),
    nrow(nitrite_clean),
    n_distinct(paste(nitrite_clean$Participant_id, nitrite_clean$tp_lab, sep = "__")),
    nrow(joined),
    sum(joined$nitrite_binary %in% c("positive", "negative"), na.rm = TRUE),
    sum(joined$nitrite_binary == "positive", na.rm = TRUE),
    sum(joined$nitrite_binary == "negative", na.rm = TRUE),
    nrow(unmatched),
    sum(unmatched$Event_type == "UTI_event", na.rm = TRUE),
    sum(unmatched$Event_type == "Routine", na.rm = TRUE)
  )
)

write_csv(joined, file.path(out_dir, "nitrate_barcode_nitrite_joined.csv"))
write_csv(unmatched, file.path(out_dir, "nitrate_barcode_nitrite_unmatched.csv"))
write_csv(by_run, file.path(out_dir, "nitrite_by_sequencing_run.csv"))
write_csv(by_barcode_number, file.path(out_dir, "nitrite_by_barcode_number.csv"))
write_csv(by_run_barcode, file.path(out_dir, "nitrite_by_run_barcode.csv"))
write_csv(by_profile, file.path(out_dir, "nitrite_by_nitrate_gene_profile.csv"))
write_csv(join_qc, file.path(out_dir, "nitrite_barcode_join_qc.csv"))

summary_file <- file.path(out_dir, "nitrate_blast_summary.txt")
writeLines(c(
  "Yellow routine FASTA nitrate/nitrite gene BLAST + nitrite barcode workflow",
  paste("Run date:", Sys.time()),
  paste("Assemblies screened:", nrow(selected_fastas)),
  paste("Reference genes screened:", nrow(reference_lengths)),
  paste("Presence threshold: pident >=", identity_threshold,
        "and query coverage >=", coverage_threshold),
  paste("Joined rows:", nrow(joined)),
  paste("Matched binary nitrite rows:",
        sum(joined$nitrite_binary %in% c("positive", "negative"), na.rm = TRUE)),
  paste("Nitrite positive rows:", sum(joined$nitrite_binary == "positive", na.rm = TRUE)),
  paste("Nitrite negative rows:", sum(joined$nitrite_binary == "negative", na.rm = TRUE)),
  paste("Unmatched/missing rows:", nrow(unmatched)),
  "",
  "Key outputs:",
  "- selected_yellow_routine_fastas.csv",
  "- nitrate_gene_blast_hits_raw.csv",
  "- nitrate_gene_best_hits.csv",
  "- nitrate_gene_absence_or_below_threshold_checklist.csv",
  "- nitrate_gene_presence_long.csv",
  "- nitrate_gene_presence_by_isolate.csv",
  "- nitrate_operon_summary_by_isolate.csv",
  "- nitrate_trend_ready_needs_nitrite.csv",
  "- nitrate_gene_presence_counts.csv",
  "- nitrate_barcode_nitrite_joined.csv",
  "- nitrate_barcode_nitrite_unmatched.csv",
  "- nitrite_by_sequencing_run.csv",
  "- nitrite_by_barcode_number.csv",
  "- nitrite_by_run_barcode.csv",
  "- nitrite_by_nitrate_gene_profile.csv",
  "- nitrite_barcode_join_qc.csv"
), summary_file)

message("Done. Outputs written to: ", out_dir)
print(join_qc)

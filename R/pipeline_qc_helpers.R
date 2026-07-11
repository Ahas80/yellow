# ==============================================================================
# R/pipeline_qc_helpers.R
# ------------------------------------------------------------------------------
# Shared source-level QA helpers for the rUTI/YELLOW RoUTIne pipeline.
# These helpers are intentionally small and dependency-light so they can be
# sourced by early clinical scripts and WGS/VF scripts without changing their
# original purpose.
# ==============================================================================

`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

safe_chr <- function(x) {
  out <- as.character(x)
  out[is.na(out)] <- NA_character_
  out
}

sanitize_key_part <- function(x) {
  x <- safe_chr(x)
  x <- trimws(x)
  x[is.na(x) | x == ""] <- "NA"
  gsub("[^A-Za-z0-9._-]+", "_", x)
}

normalise_timepoint_preserve_events <- function(x) {
  x <- trimws(safe_chr(x))
  upper <- toupper(x)
  num <- suppressWarnings(as.integer(stringr::str_extract(x, "[0-9]+")))
  dplyr::case_when(
    is.na(x) | x == "" ~ NA_character_,
    stringr::str_detect(x, stringr::regex("uricult", ignore_case = TRUE)) ~ "Uricult",
    stringr::str_detect(upper, "^UTI[-_A-Za-z0-9?]*$") ~
      stringr::str_replace_all(upper, "_", "-"),
    stringr::str_detect(upper, "^T[0-9]+(\\.[0-9]+)?$") ~
      paste0("T", stringr::str_extract(upper, "[0-9]+(\\.[0-9]+)?")),
    stringr::str_detect(x, "^[0-9]+$") ~ paste0("T", x),
    !is.na(num) & stringr::str_detect(x, stringr::regex("^(timepoint|time point)\\s*[0-9]+$", ignore_case = TRUE)) ~
      paste0("T", num),
    TRUE ~ x
  )
}

episode_event_type <- function(timepoint) {
  tp <- normalise_timepoint_preserve_events(timepoint)
  dplyr::case_when(
    is.na(tp) ~ "Unknown",
    stringr::str_detect(tp, stringr::regex("uricult", ignore_case = TRUE)) ~ "Uricult",
    stringr::str_detect(tp, "^UTI") ~ "UTI_event",
    stringr::str_detect(tp, "^T[0-9]+(\\.[0-9]+)?$") ~ "Routine",
    TRUE ~ "Other"
  )
}

build_episode_id <- function(df,
                             participant_col = "Participant_id",
                             timepoint_col = NULL,
                             event_col = NULL,
                             date_col = NULL) {
  if (!participant_col %in% names(df)) {
    stop("build_episode_id: missing participant column: ", participant_col)
  }
  if (is.null(timepoint_col)) {
    timepoint_col <- if ("Timepoint" %in% names(df)) "Timepoint" else if ("tp_lab" %in% names(df)) "tp_lab" else NULL
  }
  if (is.null(timepoint_col) || !timepoint_col %in% names(df)) {
    stop("build_episode_id: missing timepoint column")
  }

  tp <- normalise_timepoint_preserve_events(df[[timepoint_col]])
  event <- if (!is.null(event_col) && event_col %in% names(df)) {
    safe_chr(df[[event_col]])
  } else {
    episode_event_type(tp)
  }
  date_part <- if (!is.null(date_col) && date_col %in% names(df)) {
    safe_chr(df[[date_col]])
  } else if ("Collection_Date" %in% names(df)) {
    safe_chr(df$Collection_Date)
  } else {
    NA_character_
  }

  base <- paste(
    "EP",
    sanitize_key_part(df[[participant_col]]),
    sanitize_key_part(tp),
    sanitize_key_part(event),
    ifelse(is.na(date_part) | date_part == "", "nodate", sanitize_key_part(date_part)),
    sep = "__"
  )

  # Only append an occurrence suffix when the candidate ID is not unique. This
  # avoids silently merging true duplicate clinical/event rows while keeping IDs
  # stable for the common one-row-per-event case.
  dup_idx <- ave(seq_along(base), base, FUN = seq_along)
  dup_n <- ave(seq_along(base), base, FUN = length)
  ifelse(dup_n > 1, paste0(base, "__occ", dup_idx), base)
}

normalise_isolate_id <- function(x) {
  if (length(x) == 0) return(character())
  vapply(x, function(xx) {
    if (is.na(xx) || !nzchar(trimws(as.character(xx)))) return(NA_character_)
    stem <- basename(as.character(xx))
    stem <- stringr::str_remove(stem, stringr::regex("\\.(fasta|fna|fa)(\\.gz)?$", ignore_case = TRUE))
    stem <- stringr::str_replace_all(stem, "\\s+", "_")
    stem <- stringr::str_replace_all(stem, "\\.+", "_")
    stem <- stringr::str_replace_all(stem, "_+", "_")
    suffix_pattern <- stringr::regex(
      "([_-](flye|long[-_]?cycler|unicycler|assembly|assemblies|contigs?|scaffolds?|consensus|polished|medaka|racon|final|draft|complete|circulari[sz]ed|chromosome|plasmid))+$",
      ignore_case = TRUE
    )
    for (i in seq_len(4)) {
      new_stem <- stringr::str_remove(stem, suffix_pattern)
      if (identical(new_stem, stem)) break
      stem <- new_stem
    }
    stem <- stringr::str_replace_all(stem, "_+$", "")

    barcode_match <- stringr::str_match(
      stem,
      stringr::regex("barcode[0-9]+[_-]+([A-Za-z0-9]+(?:-[0-9]+)?)", ignore_case = TRUE)
    )[, 2]
    if (!is.na(barcode_match)) return(stringr::str_to_upper(barcode_match))

    c_match <- stringr::str_extract(stem, stringr::regex("\\d{4}C\\d{4,}(?:-\\d+)?", ignore_case = TRUE))
    if (!is.na(c_match)) return(stringr::str_to_upper(c_match))

    numeric_match <- stringr::str_extract(stem, "\\d{8,}(?:-\\d+)?")
    if (!is.na(numeric_match)) return(stringr::str_to_upper(numeric_match))

    short_barcode_match <- stringr::str_match(
      stem,
      stringr::regex("barcode[0-9]+[_-]+([^_/-]+(?:-[0-9]+)?)", ignore_case = TRUE)
    )[, 2]
    if (!is.na(short_barcode_match)) return(stringr::str_to_upper(short_barcode_match))

    stringr::str_to_upper(stem)
  }, character(1), USE.NAMES = FALSE)
}

detect_assembler <- function(x) {
  x <- safe_chr(x)
  dplyr::case_when(
    stringr::str_detect(x, stringr::regex("long[-_]?cycler", ignore_case = TRUE)) ~ "longcycler",
    stringr::str_detect(x, stringr::regex("flye", ignore_case = TRUE)) ~ "flye",
    stringr::str_detect(x, stringr::regex("unicycler", ignore_case = TRUE)) ~ "unicycler",
    TRUE ~ "unknown"
  )
}

fasta_extension <- function(x) {
  tolower(stringr::str_extract(basename(x), stringr::regex("\\.(fasta|fna|fa)(\\.gz)?$", ignore_case = TRUE)))
}

discover_project_fastas <- function(fasta_dir = DIR_FASTAS,
                                    recursive = TRUE,
                                    include_excluded = TRUE) {
  allowed_regex <- "\\.(fasta|fa|fna)(\\.gz)?$"
  if (!dir.exists(fasta_dir)) {
    return(tibble::tibble(
      full_path = character(), relative_path = character(), file_name = character(),
      directory = character(), extension = character(), file_stem = character(),
      Isolate_ID = character(), Assembler = character(), include_in_metadata = logical(),
      fasta_class = character(), exclusion_reason = character()
    ))
  }

  paths <- list.files(fasta_dir, pattern = allowed_regex, recursive = recursive,
                      full.names = TRUE, ignore.case = TRUE)
  if (length(paths) == 0) {
    return(tibble::tibble(
      full_path = character(), relative_path = character(), file_name = character(),
      directory = character(), extension = character(), file_stem = character(),
      Isolate_ID = character(), Assembler = character(), include_in_metadata = logical(),
      fasta_class = character(), exclusion_reason = character()
    ))
  }

  root_norm <- normalizePath(fasta_dir, winslash = "/", mustWork = FALSE)
  full_norm <- normalizePath(paths, winslash = "/", mustWork = FALSE)
  rel <- sub(paste0("^", gsub("([\\W])", "\\\\\\1", root_norm), "/?"), "", full_norm)
  dir_rel <- dirname(rel)
  dir_rel[dir_rel == "."] <- ""

  excluded_dir <- stringr::str_detect(
    rel,
    stringr::regex("(^|/)(_|\\.)?(quarantine|legacy|cache|tmp|temp|old|archive)(/|$)", ignore_case = TRUE)
  )
  output_like_name <- stringr::str_detect(
    basename(full_norm),
    stringr::regex("^(core_alignment|combined_DNA|combined_protein|multi_|parsnp|consensus_)", ignore_case = TRUE)
  )
  nested_unexpected <- dir_rel != "" & !excluded_dir

  include <- !excluded_dir & !output_like_name & !nested_unexpected
  class <- dplyr::case_when(
    excluded_dir ~ "legacy_or_quarantine",
    output_like_name ~ "derived_output_or_cache",
    nested_unexpected ~ "outside_expected_input_directory",
    TRUE ~ "candidate_input_fasta"
  )
  reason <- dplyr::case_when(
    excluded_dir ~ "Path is under quarantine/legacy/cache/temp/archive directory",
    output_like_name ~ "Filename looks like derived core/alignment/combined output",
    nested_unexpected ~ "Nested FASTA outside expected top-level FASTA input",
    TRUE ~ NA_character_
  )

  out <- tibble::tibble(
    full_path = full_norm,
    relative_path = rel,
    file_name = basename(full_norm),
    directory = dirname(full_norm),
    extension = fasta_extension(full_norm),
    file_stem = stringr::str_remove(basename(full_norm), stringr::regex("\\.(fasta|fna|fa)(\\.gz)?$", ignore_case = TRUE)),
    Isolate_ID = normalise_isolate_id(full_norm),
    Assembler = detect_assembler(full_norm),
    include_in_metadata = include,
    fasta_class = class,
    exclusion_reason = reason
  )

  if (!include_excluded) out <- dplyr::filter(out, include_in_metadata)
  out
}

assert_unique_keys <- function(df, keys, context = "table", out_path = NULL, stop_on_duplicate = FALSE) {
  missing <- setdiff(keys, names(df))
  if (length(missing) > 0) stop("assert_unique_keys: missing keys in ", context, ": ", paste(missing, collapse = ", "))
  dupes <- df |>
    dplyr::count(dplyr::across(dplyr::all_of(keys)), name = "n") |>
    dplyr::filter(.data$n > 1)
  if (nrow(dupes) > 0 && !is.null(out_path)) {
    dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
    rows <- df |> dplyr::semi_join(dupes |> dplyr::select(dplyr::all_of(keys)), by = keys)
    readr::write_csv(rows, out_path)
  }
  if (nrow(dupes) > 0 && stop_on_duplicate) {
    stop(context, " has ", nrow(dupes), " duplicated key combination(s): ", paste(keys, collapse = " + "))
  }
  invisible(dupes)
}

hash_input_manifest <- function(df, cols = NULL) {
  if (!is.null(cols)) df <- df[, intersect(cols, names(df)), drop = FALSE]
  df <- as.data.frame(df, stringsAsFactors = FALSE)
  if (nrow(df) > 0) {
    df[] <- lapply(df, function(x) ifelse(is.na(x), "<NA>", as.character(x)))
    ord_cols <- names(df)
    df <- df[do.call(order, df[ord_cols]), , drop = FALSE]
  }
  tmp <- tempfile("manifest_hash_")
  on.exit(unlink(tmp), add = TRUE)
  writeLines(c(paste(names(df), collapse = "\t"), apply(df, 1, paste, collapse = "\t")), tmp)
  unname(tools::md5sum(tmp))
}

write_input_manifest <- function(paths, out_file, role = "input") {
  paths <- unique(paths[!is.na(paths) & nzchar(paths)])
  manifest <- tibble::tibble(
    role = role,
    path = normalizePath(paths, winslash = "/", mustWork = FALSE),
    exists = file.exists(paths),
    file_size = ifelse(file.exists(paths), file.size(paths), NA_real_),
    modified_time = ifelse(file.exists(paths), as.character(file.info(paths)$mtime), NA_character_),
    md5 = ifelse(file.exists(paths), unname(tools::md5sum(paths)), NA_character_)
  )
  dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)
  readr::write_csv(manifest, out_file)
  manifest
}

find_prokka_bin <- function(bin = "prokka") {
  candidates <- unique(c(
    file.path(Sys.getenv("CONDA_PREFIX"), "bin", bin),
    file.path(Sys.getenv("HOME"), "miniforge_x86/envs/yellow-wgs-x86/bin", bin),
    file.path(Sys.getenv("HOME"), "miniconda3/envs/asm-snp-x86/bin", bin),
    file.path(Sys.getenv("HOME"), "miniconda3/envs/wgs/bin", bin),
    file.path(Sys.getenv("HOME"), "miniconda3/envs/panaroo/bin", bin),
    bin
  ))
  hits <- Sys.which(candidates)
  hits <- hits[hits != ""]
  if (length(hits) == 0) return("")
  unname(hits[1])
}

set_prokka_java_home <- function() {
  conda_prefix <- Sys.getenv("CONDA_PREFIX", "")
  java_candidates <- c(
    file.path(conda_prefix, "lib", "jvm"),
    file.path(Sys.getenv("HOME"), "miniforge_x86/envs/yellow-wgs-x86/lib/jvm"),
    file.path(Sys.getenv("HOME"), "miniconda3/envs/asm-snp-x86/lib/jvm"),
    file.path(Sys.getenv("HOME"), "miniconda3/envs/wgs/lib/jvm"),
    file.path(Sys.getenv("HOME"), "miniconda3/envs/panaroo/lib/jvm")
  )
  java_home <- java_candidates[dir.exists(java_candidates)][1]
  if (!is.na(java_home)) Sys.setenv(JAVA_HOME = java_home)
  invisible(java_home)
}

collect_pipeline_gffs <- function(gff_dirs = c(DIR_PROKKA, DIR_PROKKA_SLIM)) {
  gff_dirs <- unique(gff_dirs[!is.na(gff_dirs) & nzchar(gff_dirs)])
  existing_dirs <- gff_dirs[dir.exists(gff_dirs)]
  if (length(existing_dirs) == 0) return(character())

  all_gffs <- unlist(lapply(existing_dirs, function(d) {
    list.files(d, pattern = "\\.gff$", recursive = TRUE, full.names = TRUE, ignore.case = TRUE)
  }), use.names = FALSE)
  unique(normalizePath(all_gffs, winslash = "/", mustWork = FALSE))
}

choose_pipeline_gff <- function(assembly_id, assembly_base_id, all_gffs) {
  if (length(all_gffs) == 0) return(NA_character_)
  lookup_ids <- unique(stats::na.omit(c(as.character(assembly_id), as.character(assembly_base_id))))
  lookup_ids <- lookup_ids[nzchar(lookup_ids)]
  if (length(lookup_ids) == 0) return(NA_character_)

  matches <- all_gffs[vapply(all_gffs, function(path) {
    any(vapply(lookup_ids, function(id) {
      stringr::str_detect(basename(path), stringr::fixed(id)) ||
        stringr::str_detect(dirname(path), stringr::fixed(id))
    }, logical(1)))
  }, logical(1))]
  if (length(matches) == 0) return(NA_character_)

  pref <- matches[stringr::str_detect(matches, stringr::regex("prefixed|slim|min270", ignore_case = TRUE))]
  if (length(pref) > 0) return(sort(pref)[1])
  sort(matches)[1]
}

attach_pipeline_gff_paths <- function(expected,
                                      gff_dirs = c(DIR_PROKKA, DIR_PROKKA_SLIM)) {
  if (nrow(expected) == 0) {
    expected$gff_path <- character()
    expected$gff_available <- logical()
    expected$gff_file_size <- numeric()
    expected$gff_modified_time <- character()
    return(expected)
  }

  all_gffs <- collect_pipeline_gffs(gff_dirs)
  gff_path <- rep(NA_character_, nrow(expected))
  if (length(all_gffs) > 0) {
    gff_index <- tibble::tibble(
      gff_path = all_gffs,
      basename_id = stringr::str_remove(basename(all_gffs), stringr::regex("\\.gff$", ignore_case = TRUE)),
      dir_id = basename(dirname(all_gffs)),
      preferred = stringr::str_detect(all_gffs, stringr::regex("prefixed|slim|min270", ignore_case = TRUE))
    ) |>
      dplyr::mutate(
        basename_clean_id = stringr::str_remove(.data$basename_id, stringr::regex("(\\.(min270|prefixed|slim))+$", ignore_case = TRUE))
      )
    gff_keys <- dplyr::bind_rows(
      gff_index |> dplyr::transmute(key = .data$basename_id, .data$gff_path, .data$preferred),
      gff_index |> dplyr::transmute(key = .data$basename_clean_id, .data$gff_path, .data$preferred),
      gff_index |> dplyr::transmute(key = .data$dir_id, .data$gff_path, .data$preferred)
    ) |>
      dplyr::filter(!is.na(.data$key), nzchar(.data$key)) |>
      dplyr::arrange(dplyr::desc(.data$preferred), .data$gff_path) |>
      dplyr::distinct(.data$key, .keep_all = TRUE)

    base_match <- match(as.character(expected$Assembly_Base_ID), gff_keys$key)
    id_match <- match(as.character(expected$Assembly_ID), gff_keys$key)
    take_base <- !is.na(base_match)
    gff_path[take_base] <- gff_keys$gff_path[base_match[take_base]]
    take_id <- is.na(gff_path) & !is.na(id_match)
    gff_path[take_id] <- gff_keys$gff_path[id_match[take_id]]

    fallback <- if (identical(Sys.getenv("GFF_ALLOW_SUBSTRING_MATCH", "0"), "1")) which(is.na(gff_path)) else integer()
    if (length(fallback) > 0) {
      gff_path[fallback] <- mapply(
        choose_pipeline_gff,
        expected$Assembly_ID[fallback],
        expected$Assembly_Base_ID[fallback],
        MoreArgs = list(all_gffs = all_gffs),
        USE.NAMES = FALSE
      )
    }
  }

  gff_exists <- !is.na(gff_path) & file.exists(gff_path)
  gff_size <- rep(NA_real_, length(gff_path))
  gff_mtime <- rep(NA_character_, length(gff_path))
  if (any(gff_exists)) {
    gff_size[gff_exists] <- file.size(gff_path[gff_exists])
    gff_mtime[gff_exists] <- as.character(file.info(gff_path[gff_exists])$mtime)
  }

  expected$gff_path <- gff_path
  expected$gff_available <- gff_exists & !is.na(gff_size) & gff_size > 0
  expected$gff_file_size <- gff_size
  expected$gff_modified_time <- gff_mtime
  expected
}

write_pipeline_gff_manifest <- function(manifest, manifest_file, missing_report) {
  dir.create(dirname(manifest_file), recursive = TRUE, showWarnings = FALSE)
  dir.create(dirname(missing_report), recursive = TRUE, showWarnings = FALSE)
  readr::write_csv(manifest, manifest_file)

  missing <- manifest |>
    dplyr::filter(!(.data$gff_available %in% TRUE)) |>
    dplyr::select(dplyr::any_of(c(
      "Assembly_ID", "Assembly_Base_ID", "Participant_id", "tp_lab",
      "Assembler", "fasta_path"
    )))
  readr::write_csv(missing, missing_report)
  invisible(missing)
}

as_pipeline_bool <- function(x, default = FALSE) {
  if (is.null(x)) return(rep(default, 0))
  if (is.logical(x)) return(ifelse(is.na(x), default, x))
  y <- tolower(trimws(as.character(x)))
  out <- y %in% c("true", "t", "1", "yes", "y")
  out[is.na(y) | y == ""] <- default
  out
}

curation_default_columns <- function(df) {
  df <- tibble::as_tibble(df)
  n <- nrow(df)
  if (!"analysis_include_primary" %in% names(df)) df$analysis_include_primary <- rep(TRUE, n)
  df$analysis_include_primary <- as_pipeline_bool(df$analysis_include_primary, default = TRUE)
  if (!"analysis_exclusion_reason" %in% names(df)) df$analysis_exclusion_reason <- rep(NA_character_, n)
  if (!"duplicate_role" %in% names(df)) df$duplicate_role <- rep("not_duplicate", n)
  df$duplicate_role <- safe_chr(df$duplicate_role)
  df$duplicate_role[is.na(df$duplicate_role) | !nzchar(trimws(df$duplicate_role))] <- "not_duplicate"
  if (!"duplicate_of_participant_id" %in% names(df)) df$duplicate_of_participant_id <- rep(NA_character_, n)
  if (!"duplicate_of_tp_lab" %in% names(df)) df$duplicate_of_tp_lab <- rep(NA_character_, n)
  if (!"duplicate_use_note" %in% names(df)) df$duplicate_use_note <- rep(NA_character_, n)
  if (!"allow_secondary_duplicate_qc" %in% names(df)) df$allow_secondary_duplicate_qc <- rep(FALSE, n)
  df$allow_secondary_duplicate_qc <- as_pipeline_bool(df$allow_secondary_duplicate_qc, default = FALSE)
  if (!"genomics_expected_include" %in% names(df)) df$genomics_expected_include <- rep(TRUE, n)
  df$genomics_expected_include <- as_pipeline_bool(df$genomics_expected_include, default = TRUE)
  if (!"genomics_exclusion_reason" %in% names(df)) df$genomics_exclusion_reason <- rep(NA_character_, n)
  if (!"manual_curation_applied" %in% names(df)) df$manual_curation_applied <- rep(FALSE, n)
  df$manual_curation_applied <- as_pipeline_bool(df$manual_curation_applied, default = FALSE)
  if (!"manual_curation_note" %in% names(df)) df$manual_curation_note <- rep(NA_character_, n)
  if (!"manual_curation_source" %in% names(df)) df$manual_curation_source <- rep(NA_character_, n)
  df
}

read_manual_sample_curation <- function(path = if (exists("FILE_MANUAL_SAMPLE_CURATION", inherits = TRUE)) FILE_MANUAL_SAMPLE_CURATION else file.path("data", "manual_sample_curation.csv")) {
  required <- c(
    "Participant_id", "tp_lab", "Timepoint", "Isolate_ID", "UTI_Label",
    "exclude_primary", "exclude_from_genomics_expected", "exclude_reason",
    "duplicate_role", "duplicate_of_participant_id", "duplicate_of_tp_lab",
    "allow_secondary_duplicate_qc", "duplicate_use_note", "curation_note",
    "curation_source"
  )
  if (!file.exists(path)) {
    out <- tibble::as_tibble(stats::setNames(rep(list(character()), length(required)), required))
    out$exclude_primary <- logical()
    out$exclude_from_genomics_expected <- logical()
    out$allow_secondary_duplicate_qc <- logical()
    return(out)
  }

  cur <- readr::read_csv(path, show_col_types = FALSE, na = c("", "NA", "N/A")) |>
    tibble::as_tibble()
  missing <- setdiff(required, names(cur))
  for (col in missing) cur[[col]] <- NA
  cur <- cur |>
    dplyr::mutate(
      dplyr::across(dplyr::any_of(c(
        "Participant_id", "tp_lab", "Timepoint", "Isolate_ID", "UTI_Label",
        "exclude_reason", "duplicate_role", "duplicate_of_participant_id",
        "duplicate_of_tp_lab", "duplicate_use_note", "curation_note",
        "curation_source"
      )), safe_chr),
      Participant_id = stringr::str_trim(.data$Participant_id),
      tp_lab = dplyr::if_else(
        is.na(.data$tp_lab) | !nzchar(stringr::str_trim(.data$tp_lab)),
        normalise_timepoint_preserve_events(.data$Timepoint),
        normalise_timepoint_preserve_events(.data$tp_lab)
      ),
      Isolate_ID = normalise_isolate_id(.data$Isolate_ID),
      UTI_Label = stringr::str_trim(.data$UTI_Label),
      exclude_primary = as_pipeline_bool(.data$exclude_primary, default = FALSE),
      exclude_from_genomics_expected = as_pipeline_bool(.data$exclude_from_genomics_expected, default = FALSE),
      allow_secondary_duplicate_qc = as_pipeline_bool(.data$allow_secondary_duplicate_qc, default = FALSE),
      duplicate_role = dplyr::if_else(
        is.na(.data$duplicate_role) | !nzchar(stringr::str_trim(.data$duplicate_role)),
        "not_duplicate",
        .data$duplicate_role
      )
    )
  cur$.curation_row_id <- seq_len(nrow(cur))
  cur
}

apply_manual_sample_curation <- function(df,
                                         context = "unspecified",
                                         path = if (exists("FILE_MANUAL_SAMPLE_CURATION", inherits = TRUE)) FILE_MANUAL_SAMPLE_CURATION else file.path("data", "manual_sample_curation.csv"),
                                         write_audit = FALSE,
                                         audit_path = if (exists("FILE_SAMPLE_CURATION_AUDIT", inherits = TRUE)) FILE_SAMPLE_CURATION_AUDIT else file.path("results", "qc", "manual_sample_curation_applied.csv")) {
  df <- curation_default_columns(df)
  if (nrow(df) == 0) return(df)

  cur <- read_manual_sample_curation(path)
  if (nrow(cur) == 0) return(df)
  if (!all(c("Participant_id", "tp_lab") %in% names(df))) {
    if ("Participant_id" %in% names(df) && "Timepoint" %in% names(df)) {
      df$tp_lab <- normalise_timepoint_preserve_events(df$Timepoint)
    } else {
      return(df)
    }
  }

  isolate_col <- first_existing_col(df, c("Isolate_ID", "isolate_ID"))
  label_col <- first_existing_col(df, c("UTI_Label", "UTI label", "UTI_label"))
  keyed <- df |>
    dplyr::mutate(
      .curation_row_index = dplyr::row_number(),
      .curation_participant = stringr::str_trim(safe_chr(.data$Participant_id)),
      .curation_tp = normalise_timepoint_preserve_events(.data$tp_lab),
      .curation_isolate = if (!is.null(isolate_col)) normalise_isolate_id(.data[[isolate_col]]) else NA_character_,
      .curation_uti_label = if (!is.null(label_col)) stringr::str_trim(safe_chr(.data[[label_col]])) else NA_character_
    )
  cur_keyed <- cur |>
    dplyr::mutate(
      .curation_participant = stringr::str_trim(safe_chr(.data$Participant_id)),
      .curation_tp = normalise_timepoint_preserve_events(.data$tp_lab),
      .curation_isolate = normalise_isolate_id(.data$Isolate_ID),
      .curation_uti_label = stringr::str_trim(safe_chr(.data$UTI_Label)),
      .curation_has_isolate = !is.na(.data$.curation_isolate) & nzchar(.data$.curation_isolate),
      .curation_has_label = !is.na(.data$.curation_uti_label) & nzchar(.data$.curation_uti_label)
    )

  candidates <- keyed |>
    dplyr::select(dplyr::all_of(c(
      ".curation_row_index", ".curation_participant", ".curation_tp",
      ".curation_isolate", ".curation_uti_label"
    ))) |>
    dplyr::inner_join(
      cur_keyed |>
        dplyr::select(dplyr::all_of(c(
          ".curation_row_id", ".curation_participant", ".curation_tp",
          ".curation_isolate", ".curation_uti_label",
          ".curation_has_isolate", ".curation_has_label",
          "exclude_primary", "exclude_from_genomics_expected",
          "exclude_reason", "duplicate_role",
          "duplicate_of_participant_id", "duplicate_of_tp_lab",
          "allow_secondary_duplicate_qc", "duplicate_use_note",
          "curation_note", "curation_source"
        ))) |>
        dplyr::rename(.cur_isolate = .curation_isolate,
                      .cur_uti_label = .curation_uti_label),
      by = c(".curation_participant", ".curation_tp"),
      relationship = "many-to-many"
    ) |>
    dplyr::mutate(
      .df_has_isolate = !is.na(.data$.curation_isolate) & nzchar(.data$.curation_isolate),
      .df_has_label = !is.na(.data$.curation_uti_label) & nzchar(.data$.curation_uti_label),
      .isolate_match = .data$.df_has_isolate & .data$.curation_isolate == .data$.cur_isolate,
      .label_match = .data$.df_has_label & .data$.curation_uti_label == .data$.cur_uti_label,
      .match_ok = dplyr::case_when(
        .data$.curation_has_isolate & .data$.curation_has_label ~ .data$.isolate_match | .data$.label_match,
        .data$.curation_has_isolate & !.data$.curation_has_label ~ .data$.isolate_match,
        !.data$.curation_has_isolate & .data$.curation_has_label ~ .data$.label_match,
        TRUE ~ TRUE
      )
    ) |>
    dplyr::filter(.data$.match_ok)

  if (nrow(candidates) == 0) return(df)

  dup_matches <- candidates |>
    dplyr::count(.data$.curation_row_index, name = "n") |>
    dplyr::filter(.data$n > 1)
  if (nrow(dup_matches) > 0) {
    stop("Manual curation produced multiple matches in ", context,
         ". Refine data/manual_sample_curation.csv before continuing.")
  }

  match <- candidates |>
    dplyr::arrange(.data$.curation_row_index) |>
    dplyr::distinct(.data$.curation_row_index, .keep_all = TRUE)
  idx <- match$.curation_row_index

  df$manual_curation_applied[idx] <- TRUE
  df$analysis_include_primary[idx] <- df$analysis_include_primary[idx] & !match$exclude_primary
  df$genomics_expected_include[idx] <- df$genomics_expected_include[idx] & !match$exclude_from_genomics_expected

  primary_reason <- match$exclude_reason[match$exclude_primary]
  primary_idx <- idx[match$exclude_primary]
  if (length(primary_idx) > 0) df$analysis_exclusion_reason[primary_idx] <- primary_reason

  genomics_reason <- match$exclude_reason[match$exclude_from_genomics_expected]
  genomics_idx <- idx[match$exclude_from_genomics_expected]
  if (length(genomics_idx) > 0) df$genomics_exclusion_reason[genomics_idx] <- genomics_reason

  replacement <- function(values, current, default_blank = NA_character_) {
    out <- current
    take <- !is.na(values) & nzchar(stringr::str_trim(values))
    out[idx[take]] <- values[take]
    out[is.na(out) | !nzchar(stringr::str_trim(out))] <- default_blank
    out
  }
  df$duplicate_role <- replacement(match$duplicate_role, df$duplicate_role, "not_duplicate")
  df$duplicate_of_participant_id <- replacement(match$duplicate_of_participant_id, df$duplicate_of_participant_id)
  df$duplicate_of_tp_lab <- replacement(match$duplicate_of_tp_lab, df$duplicate_of_tp_lab)
  df$duplicate_use_note <- replacement(match$duplicate_use_note, df$duplicate_use_note)
  df$manual_curation_note <- replacement(match$curation_note, df$manual_curation_note)
  df$manual_curation_source <- replacement(match$curation_source, df$manual_curation_source)
  df$allow_secondary_duplicate_qc[idx] <- df$allow_secondary_duplicate_qc[idx] | match$allow_secondary_duplicate_qc

  if (write_audit && length(idx) > 0) {
    dir.create(dirname(audit_path), recursive = TRUE, showWarnings = FALSE)
    audit <- df[idx, , drop = FALSE] |>
      dplyr::mutate(curation_context = context, .before = 1) |>
      dplyr::select(dplyr::any_of(c(
        "curation_context", "Participant_id", "Timepoint", "tp_lab", "Isolate_ID",
        "UTI_Label", "UTI_Status", "Infection_Status", "analysis_include_primary",
        "analysis_exclusion_reason", "genomics_expected_include",
        "genomics_exclusion_reason", "duplicate_role",
        "duplicate_of_participant_id", "duplicate_of_tp_lab",
        "allow_secondary_duplicate_qc", "duplicate_use_note",
        "manual_curation_note", "manual_curation_source"
      ))) |>
      dplyr::mutate(dplyr::across(dplyr::everything(), as.character))
    existing <- if (file.exists(audit_path)) {
      readr::read_csv(audit_path, show_col_types = FALSE) |>
        dplyr::mutate(dplyr::across(dplyr::everything(), as.character))
    } else {
      tibble::tibble()
    }
    if ("curation_context" %in% names(existing)) existing <- dplyr::filter(existing, .data$curation_context != .env$context)
    readr::write_csv(dplyr::bind_rows(existing, audit), audit_path)
  }

  df
}

filter_primary_analysis <- function(df) {
  df <- curation_default_columns(df)
  dplyr::filter(df, .data$analysis_include_primary %in% TRUE)
}

filter_primary_genomics <- function(df) {
  df <- curation_default_columns(df)
  dplyr::filter(df, .data$analysis_include_primary %in% TRUE,
                .data$genomics_expected_include %in% TRUE)
}

prefer_primary_uti_status <- function(df, status_col = "Infection_Status", allow_legacy_fallback = FALSE) {
  if (is.null(df)) return(df)
  df <- tibble::as_tibble(df)

  if (!"UTI_Status" %in% names(df)) {
    if (status_col %in% names(df)) {
      status_vals <- unique(stats::na.omit(safe_chr(df[[status_col]])))
      if (all(status_vals %in% c("UTI", "Not_UTI"))) {
        df$UTI_Status <- safe_chr(df[[status_col]])
      } else if (allow_legacy_fallback) {
        df$UTI_Status <- ifelse(safe_chr(df[[status_col]]) == "UTI", "UTI", "Not_UTI")
      } else {
        df$UTI_Status <- NA_character_
      }
    } else {
      df$UTI_Status <- NA_character_
    }
  }
  df$UTI_Status <- safe_chr(df$UTI_Status)
  df$UTI_Status[is.na(df$UTI_Status) | !nzchar(df$UTI_Status)] <- NA_character_

  if (!"UTI_binary" %in% names(df)) {
    df$UTI_binary <- ifelse(df$UTI_Status == "UTI", 1L, ifelse(df$UTI_Status == "Not_UTI", 0L, NA_integer_))
  }
  if (!"Not_UTI_subgroup" %in% names(df)) {
    df$Not_UTI_subgroup <- ifelse(df$UTI_Status == "Not_UTI", "unknown_or_indeterminate", NA_character_)
  }
  if (!"UTI_definition_version" %in% names(df)) {
    df$UTI_definition_version <- if (exists("UTI_DEFINITION_VERSION", inherits = TRUE)) UTI_DEFINITION_VERSION else "legacy_fallback"
  }

  if (status_col %in% names(df)) {
    if (!"Infection_Status_legacy" %in% names(df)) df$Infection_Status_legacy <- safe_chr(df[[status_col]])
    if (!"Infection_Status_old" %in% names(df)) df$Infection_Status_old <- df$Infection_Status_legacy
  } else {
    df[[status_col]] <- NA_character_
  }

  df$Primary_Status <- df$UTI_Status
  df[[status_col]] <- df$UTI_Status
  df
}

status_map_has_primary_fields <- function(path) {
  if (!file.exists(path)) return(FALSE)
  cols <- tryCatch(
    names(readr::read_csv(path, n_max = 0, show_col_types = FALSE)),
    error = function(e) character()
  )
  all(c("UTI_Status", "UTI_binary") %in% cols)
}

status_map_is_stale <- function(target, upstream = FILE_STATUS_MAP) {
  file.exists(target) && file.exists(upstream) &&
    file.info(target)$mtime < file.info(upstream)$mtime
}

select_primary_status_map <- function(prefer_poster = TRUE,
                                      require_fresh = TRUE,
                                      caller = NULL,
                                      quiet = FALSE) {
  if (!file.exists(FILE_STATUS_MAP)) {
    stop("Missing primary status map: ", FILE_STATUS_MAP, ". Run 00b_classify_episodes.R first.")
  }
  if (!status_map_has_primary_fields(FILE_STATUS_MAP)) {
    stop("Primary status_map.csv lacks UTI_Status/UTI_binary. Rerun 00b_classify_episodes.R.")
  }

  use_poster <- prefer_poster && file.exists(FILE_STATUS_MAP_POSTER)
  poster_reason <- NULL
  if (use_poster && require_fresh && status_map_is_stale(FILE_STATUS_MAP_POSTER, FILE_STATUS_MAP)) {
    use_poster <- FALSE
    poster_reason <- "older than status_map.csv"
  }
  if (use_poster && !status_map_has_primary_fields(FILE_STATUS_MAP_POSTER)) {
    use_poster <- FALSE
    poster_reason <- "missing UTI_Status/UTI_binary"
  }

  selected <- if (use_poster) FILE_STATUS_MAP_POSTER else FILE_STATUS_MAP
  if (!quiet && prefer_poster && file.exists(FILE_STATUS_MAP_POSTER) && identical(selected, FILE_STATUS_MAP)) {
    prefix <- if (!is.null(caller) && nzchar(caller)) paste0(caller, ": ") else ""
    message(prefix, "using status_map.csv for primary UTI status; status_map_with_poster_tp.csv is ",
            poster_reason %||% "not suitable for primary status")
  }
  selected
}

read_primary_status_map <- function(prefer_poster = TRUE,
                                    require_fresh = TRUE,
                                    caller = NULL,
                                    standardise_tp = TRUE) {
  path <- select_primary_status_map(prefer_poster = prefer_poster,
                                    require_fresh = require_fresh,
                                    caller = caller)
  df <- readr::read_csv(path, show_col_types = FALSE) %>%
    prefer_primary_uti_status(allow_legacy_fallback = FALSE)
  if (standardise_tp) {
    if ("tp_lab" %in% names(df)) {
      df$tp_lab <- normalise_timepoint_preserve_events(df$tp_lab)
    } else if ("Timepoint" %in% names(df)) {
      df$tp_lab <- normalise_timepoint_preserve_events(df$Timepoint)
    }
  }
  attr(df, "source_file") <- path
  df
}

first_existing_col <- function(df, cols) {
  hit <- intersect(cols, names(df))
  if (length(hit) == 0) return(NULL)
  hit[1]
}

standardise_gff_assembly_table <- function(df, source_label = "assembly_table") {
  if (is.null(df) || nrow(df) == 0) {
    return(tibble::tibble(
      Assembly_ID = character(),
      Assembly_Base_ID = character(),
      Isolate_ID = character(),
      Participant_id = character(),
      tp_lab = character(),
      Timepoint = character(),
      Assembler = character(),
      file_name = character(),
      full_path = character(),
      fasta_path = character(),
      file_exists = logical(),
      found = logical(),
      usable_fasta = logical(),
      QC_PASS = logical(),
      selected_canonical = logical(),
      analysis_include_primary = logical(),
      analysis_exclusion_reason = character(),
      genomics_expected_include = logical(),
      genomics_exclusion_reason = character(),
      duplicate_role = character(),
      allow_secondary_duplicate_qc = logical(),
      source_label = character()
    ))
  }

  df <- tibble::as_tibble(df)
  n <- nrow(df)

  file_col <- first_existing_col(df, c("full_path", "fasta_path", "assembly", "path"))
  name_col <- first_existing_col(df, c("file_name", "filename", "FILE", "file"))

  full_path <- if (!is.null(file_col)) as.character(df[[file_col]]) else rep(NA_character_, n)
  full_path[is.na(full_path) | !nzchar(trimws(full_path))] <- NA_character_

  file_name <- if (!is.null(name_col)) as.character(df[[name_col]]) else basename(full_path)
  file_name[is.na(file_name) | !nzchar(trimws(file_name))] <- basename(full_path[is.na(file_name) | !nzchar(trimws(file_name))])
  file_name[is.na(file_name) | file_name == "NA"] <- NA_character_

  fill_from_name <- is.na(full_path) & !is.na(file_name) & nzchar(file_name)
  full_path[fill_from_name] <- file.path(DIR_FASTAS, file_name[fill_from_name])
  full_path_norm <- ifelse(
    is.na(full_path),
    NA_character_,
    normalizePath(full_path, winslash = "/", mustWork = FALSE)
  )

  current_file_exists <- !is.na(full_path_norm) & file.exists(full_path_norm)
  stem <- stringr::str_remove(basename(full_path_norm), stringr::regex("\\.(fasta|fa|fna)(\\.gz)?$", ignore_case = TRUE))
  stem[is.na(stem) | stem == "NA"] <- NA_character_
  ext <- fasta_extension(full_path_norm)
  ext_clean <- stringr::str_replace_all(ext, "\\.", "_")
  ext_clean[is.na(ext_clean) | !nzchar(ext_clean)] <- "fasta"

  assembly_base <- if ("Assembly_Base_ID" %in% names(df)) as.character(df$Assembly_Base_ID) else stem
  assembly_base[is.na(assembly_base) | !nzchar(trimws(assembly_base))] <- stem[is.na(assembly_base) | !nzchar(trimws(assembly_base))]

  assembly_id <- if ("Assembly_ID" %in% names(df)) as.character(df$Assembly_ID) else paste0(assembly_base, "__", ext_clean)
  assembly_id[is.na(assembly_id) | !nzchar(trimws(assembly_id))] <- paste0(assembly_base[is.na(assembly_id) | !nzchar(trimws(assembly_id))], "__", ext_clean[is.na(assembly_id) | !nzchar(trimws(assembly_id))])

  isolate_id <- if ("Isolate_ID" %in% names(df)) as.character(df$Isolate_ID) else normalise_isolate_id(ifelse(!is.na(file_name), file_name, full_path_norm))
  participant_id <- if ("Participant_id" %in% names(df)) as.character(df$Participant_id) else rep(NA_character_, n)
  timepoint <- if ("Timepoint" %in% names(df)) as.character(df$Timepoint) else if ("tp_lab" %in% names(df)) as.character(df$tp_lab) else rep(NA_character_, n)
  tp_lab <- if ("tp_lab" %in% names(df)) normalise_timepoint_preserve_events(df$tp_lab) else normalise_timepoint_preserve_events(timepoint)
  assembler <- if ("Assembler" %in% names(df)) as.character(df$Assembler) else if ("assembler" %in% names(df)) as.character(df$assembler) else detect_assembler(ifelse(!is.na(file_name), file_name, full_path_norm))

  found <- if ("found" %in% names(df)) as_pipeline_bool(df$found, default = !is.na(full_path_norm)) else !is.na(full_path_norm)
  usable <- if ("usable_fasta" %in% names(df)) as_pipeline_bool(df$usable_fasta, default = current_file_exists) else current_file_exists
  qc_pass <- if ("QC_PASS" %in% names(df)) as_pipeline_bool(df$QC_PASS, default = FALSE) else rep(FALSE, n)
  selected <- if ("selected_canonical" %in% names(df)) as_pipeline_bool(df$selected_canonical, default = FALSE) else rep(FALSE, n)
  analysis_include_primary <- if ("analysis_include_primary" %in% names(df)) as_pipeline_bool(df$analysis_include_primary, default = TRUE) else rep(TRUE, n)
  analysis_exclusion_reason <- if ("analysis_exclusion_reason" %in% names(df)) as.character(df$analysis_exclusion_reason) else rep(NA_character_, n)
  genomics_expected_include <- if ("genomics_expected_include" %in% names(df)) as_pipeline_bool(df$genomics_expected_include, default = TRUE) else rep(TRUE, n)
  genomics_exclusion_reason <- if ("genomics_exclusion_reason" %in% names(df)) as.character(df$genomics_exclusion_reason) else rep(NA_character_, n)
  duplicate_role <- if ("duplicate_role" %in% names(df)) as.character(df$duplicate_role) else rep("not_duplicate", n)
  duplicate_role[is.na(duplicate_role) | !nzchar(trimws(duplicate_role))] <- "not_duplicate"
  allow_secondary_duplicate_qc <- if ("allow_secondary_duplicate_qc" %in% names(df)) as_pipeline_bool(df$allow_secondary_duplicate_qc, default = FALSE) else rep(FALSE, n)

  tibble::tibble(
    Assembly_ID = assembly_id,
    Assembly_Base_ID = assembly_base,
    Isolate_ID = isolate_id,
    Participant_id = participant_id,
    tp_lab = tp_lab,
    Timepoint = timepoint,
    Assembler = assembler,
    file_name = file_name,
    full_path = full_path_norm,
    fasta_path = full_path_norm,
    file_exists = current_file_exists,
    found = found,
    usable_fasta = usable & current_file_exists,
    QC_PASS = qc_pass,
    selected_canonical = selected,
    analysis_include_primary = analysis_include_primary,
    analysis_exclusion_reason = analysis_exclusion_reason,
    genomics_expected_include = genomics_expected_include,
    genomics_exclusion_reason = genomics_exclusion_reason,
    duplicate_role = duplicate_role,
    allow_secondary_duplicate_qc = allow_secondary_duplicate_qc,
    source_label = source_label
  )
}

build_assembly_gff_inventory <- function(metadata_file = FILE_METADATA,
                                         canonical_file = file.path(DIR_QC, "canonical_assembly_selection.csv"),
                                         fasta_dir = DIR_FASTAS,
                                         gff_dirs = c(DIR_PROKKA, DIR_PROKKA_SLIM),
                                         metadata_scan_file = file.path(DIR_QC, "metadata_fasta_discovery_manifest.csv"),
                                         out_all_metadata = file.path(DIR_QC, "gff_inventory_all_metadata_fastas.csv"),
                                         out_panaroo = file.path(DIR_WGS, "pan", "panaroo_gff_inventory.csv"),
                                         out_summary = file.path(DIR_QC, "gff_inventory_summary.csv"),
                                         write_outputs = TRUE) {
  if (!file.exists(metadata_file)) {
    stop("Assembly metadata is missing: ", metadata_file, "\nRun 00_make_assembly_metadata.r first.")
  }

  current_scan <- discover_project_fastas(fasta_dir, recursive = TRUE, include_excluded = TRUE)
  current_candidates <- current_scan |>
    dplyr::filter(.data$include_in_metadata) |>
    dplyr::mutate(full_path_norm = normalizePath(.data$full_path, winslash = "/", mustWork = FALSE))
  current_candidate_paths <- unique(current_candidates$full_path_norm)

  metadata_scan <- if (file.exists(metadata_scan_file)) {
    readr::read_csv(metadata_scan_file, show_col_types = FALSE)
  } else {
    tibble::tibble(full_path = character(), include_in_metadata = logical())
  }
  scan_candidates <- metadata_scan
  if ("include_in_metadata" %in% names(scan_candidates)) {
    scan_candidates <- dplyr::filter(scan_candidates, .data$include_in_metadata)
  }
  snapshot_candidate_paths <- if ("full_path" %in% names(scan_candidates)) {
    unique(normalizePath(scan_candidates$full_path, winslash = "/", mustWork = FALSE))
  } else {
    character()
  }

  new_candidate_paths <- setdiff(current_candidate_paths, snapshot_candidate_paths)
  removed_candidate_paths <- setdiff(snapshot_candidate_paths, current_candidate_paths)

  metadata_raw <- readr::read_csv(metadata_file, show_col_types = FALSE)
  metadata_std <- standardise_gff_assembly_table(metadata_raw, "assembly_metadata")
  metadata_paths <- unique(metadata_std$full_path[!is.na(metadata_std$full_path) & nzchar(metadata_std$full_path)])

  metadata_linked <- metadata_std |>
    dplyr::filter(!is.na(.data$full_path), nzchar(.data$full_path), .data$found %in% TRUE) |>
    dplyr::mutate(
      inventory_tier = "metadata_linked_fasta",
      gff_requirement = "warning_only",
      metadata_or_scan_status = dplyr::case_when(
        !.data$file_exists ~ "metadata_path_missing_on_disk",
        .data$full_path %in% current_candidate_paths ~ "metadata_linked_current_fasta",
        TRUE ~ "metadata_path_outside_current_candidate_scan"
      )
    )

  metadata_missing_expected <- metadata_std |>
    dplyr::filter(!(.data$found %in% TRUE) | is.na(.data$full_path) | !nzchar(.data$full_path)) |>
    dplyr::mutate(
      inventory_tier = "expected_metadata_row_missing_fasta",
      gff_requirement = "not_applicable",
      metadata_or_scan_status = "expected_metadata_row_without_fasta"
    )

  unexpected_candidates <- current_candidates |>
    dplyr::filter(!(.data$full_path_norm %in% metadata_paths)) |>
    dplyr::transmute(
      Assembly_ID = paste0(.data$file_stem, "__", stringr::str_replace_all(.data$extension, "\\.", "_")),
      Assembly_Base_ID = .data$file_stem,
      Isolate_ID = .data$Isolate_ID,
      Participant_id = NA_character_,
      tp_lab = NA_character_,
      Timepoint = NA_character_,
      Assembler = .data$Assembler,
      file_name = .data$file_name,
      full_path = .data$full_path_norm,
      fasta_path = .data$full_path_norm,
      file_exists = file.exists(.data$full_path_norm),
      found = TRUE,
      usable_fasta = file.exists(.data$full_path_norm),
      QC_PASS = FALSE,
      selected_canonical = FALSE,
      source_label = "current_fasta_scan",
      inventory_tier = "unexpected_unlinked_candidate_fasta",
      gff_requirement = "not_applicable",
      metadata_or_scan_status = dplyr::case_when(
        .data$full_path_norm %in% new_candidate_paths ~ "new_since_metadata_scan",
        TRUE ~ "unexpected_unlinked_candidate_fasta"
      )
    )

  all_inventory <- dplyr::bind_rows(
    metadata_linked,
    metadata_missing_expected,
    unexpected_candidates
  ) |>
    attach_pipeline_gff_paths(gff_dirs = gff_dirs) |>
    dplyr::mutate(
      gff_missing = !(.data$gff_available %in% TRUE),
      active_expected_for_canonical_qc = .data$found %in% TRUE &
        .data$file_exists %in% TRUE &
        .data$analysis_include_primary %in% TRUE &
        .data$genomics_expected_include %in% TRUE,
      final_analysis_blocker = .data$metadata_or_scan_status %in% c(
        "metadata_path_missing_on_disk",
        "new_since_metadata_scan"
      ) & .data$active_expected_for_canonical_qc
    )

  canonical_missing <- !file.exists(canonical_file)
  canonical_raw <- if (!canonical_missing) {
    readr::read_csv(canonical_file, show_col_types = FALSE)
  } else {
    tibble::tibble()
  }
  canonical_std <- standardise_gff_assembly_table(canonical_raw, "canonical_assembly_selection")
  panaroo_inventory <- canonical_std |>
    dplyr::filter(.data$selected_canonical %in% TRUE, .data$QC_PASS %in% TRUE, .data$file_exists) |>
    dplyr::mutate(
      inventory_tier = "panaroo_eligible_canonical_qc_pass",
      gff_requirement = "required_for_panaroo",
      metadata_or_scan_status = "panaroo_input"
    ) |>
    attach_pipeline_gff_paths(gff_dirs = gff_dirs) |>
    dplyr::mutate(
      gff_missing = !(.data$gff_available %in% TRUE),
      final_analysis_blocker = .data$gff_missing
    ) |>
    dplyr::arrange(.data$Participant_id, .data$tp_lab, .data$Assembly_ID)

  metadata_existing_ids <- metadata_linked |>
    dplyr::filter(
      .data$file_exists,
      .data$found %in% TRUE,
      .data$analysis_include_primary %in% TRUE,
      .data$genomics_expected_include %in% TRUE
    ) |>
    dplyr::pull(.data$Assembly_ID) |>
    unique()
  canonical_all_ids <- canonical_std |>
    dplyr::filter(.data$file_exists) |>
    dplyr::pull(.data$Assembly_ID) |>
    unique()

  qc_missing_metadata_ids <- setdiff(metadata_existing_ids, canonical_all_ids)
  qc_extra_ids <- setdiff(canonical_all_ids, metadata_existing_ids)
  canonical_older_than_metadata <- !canonical_missing &&
    file.exists(metadata_file) &&
    file.info(canonical_file)$mtime < file.info(metadata_file)$mtime
  active_metadata_missing_on_disk <- metadata_linked |>
    dplyr::filter(
      .data$metadata_or_scan_status == "metadata_path_missing_on_disk",
      .data$analysis_include_primary %in% TRUE,
      .data$genomics_expected_include %in% TRUE
    )

  metadata_stale_messages <- c(
    if (length(new_candidate_paths) > 0) sprintf("%d new candidate FASTA(s) were added since the last metadata scan.", length(new_candidate_paths)),
    if (length(removed_candidate_paths) > 0) sprintf("%d candidate FASTA(s) from the last metadata scan are no longer present.", length(removed_candidate_paths)),
    if (nrow(active_metadata_missing_on_disk) > 0) {
      sprintf("%d active metadata-linked FASTA path(s) are missing on disk.", nrow(active_metadata_missing_on_disk))
    }
  )
  qc_stale_messages <- c(
    if (canonical_missing) sprintf("Canonical assembly selection is missing: %s", canonical_file),
    if (canonical_older_than_metadata) "Canonical assembly selection is older than assembly_metadata.csv.",
    if (length(qc_missing_metadata_ids) > 0) sprintf("%d metadata-linked FASTA(s) are absent from canonical_assembly_selection.csv.", length(qc_missing_metadata_ids)),
    if (length(qc_extra_ids) > 0) sprintf("%d canonical selection row(s) are no longer present in assembly_metadata.csv.", length(qc_extra_ids))
  )

  metadata_stale <- length(metadata_stale_messages) > 0
  qc_stale <- length(qc_stale_messages) > 0
  metadata_inventory_rows <- all_inventory |>
    dplyr::filter(.data$inventory_tier == "metadata_linked_fasta")
  metadata_gffs_available <- sum(metadata_inventory_rows$gff_available %in% TRUE, na.rm = TRUE)
  metadata_gffs_missing <- sum(!(metadata_inventory_rows$gff_available %in% TRUE), na.rm = TRUE)
  panaroo_gffs_available <- sum(panaroo_inventory$gff_available %in% TRUE, na.rm = TRUE)
  panaroo_gffs_missing <- sum(!(panaroo_inventory$gff_available %in% TRUE), na.rm = TRUE)

  summary <- tibble::tibble(
    metric = c(
      "current_candidate_fastas",
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
      "metadata_inventory_stale",
      "qc_selection_stale"
    ),
    value = as.character(c(
      length(current_candidate_paths),
      nrow(metadata_linked),
      metadata_gffs_available,
      metadata_gffs_missing,
      nrow(metadata_missing_expected),
      nrow(unexpected_candidates),
      length(new_candidate_paths),
      length(removed_candidate_paths),
      nrow(panaroo_inventory),
      panaroo_gffs_available,
      panaroo_gffs_missing,
      metadata_stale,
      qc_stale
    )),
    severity = c(
      "info", "info", "info", "warning", "warning", "warning",
      ifelse(length(new_candidate_paths) > 0, "error", "info"),
      ifelse(length(removed_candidate_paths) > 0, "error", "info"),
      "info", "info",
      ifelse(panaroo_gffs_missing > 0, "error", "info"),
      ifelse(metadata_stale, "error", "info"),
      ifelse(qc_stale, "error", "info")
    )
  )

  if (write_outputs) {
    dir.create(dirname(out_all_metadata), recursive = TRUE, showWarnings = FALSE)
    dir.create(dirname(out_panaroo), recursive = TRUE, showWarnings = FALSE)
    dir.create(dirname(out_summary), recursive = TRUE, showWarnings = FALSE)
    readr::write_csv(all_inventory, out_all_metadata)
    readr::write_csv(panaroo_inventory, out_panaroo)
    readr::write_csv(summary, out_summary)
  }

  list(
    all_inventory = all_inventory,
    metadata_linked_inventory = metadata_linked,
    panaroo_inventory = panaroo_inventory,
    summary = summary,
    metadata_stale = metadata_stale,
    qc_stale = qc_stale,
    metadata_stale_messages = metadata_stale_messages,
    qc_stale_messages = qc_stale_messages,
    paths = list(
      all_metadata = out_all_metadata,
      panaroo = out_panaroo,
      summary = out_summary,
      metadata_scan = metadata_scan_file,
      metadata = metadata_file,
      canonical = canonical_file
    )
  )
}

prepare_gff_regen_jobs <- function(missing_df,
                                   output_dir = DIR_PROKKA,
                                   log_dir = file.path(DIR_LOGS, "prokka_missing_gffs")) {
  if (is.null(missing_df)) missing_df <- tibble::tibble()
  for (col in c("Assembly_Base_ID", "fasta_path", "Assembler", "Participant_id", "tp_lab")) {
    if (!col %in% names(missing_df)) missing_df[[col]] <- NA_character_
  }

  missing_df |>
    dplyr::mutate(
      Assembly_Base_ID = as.character(.data$Assembly_Base_ID),
      fasta_path = as.character(.data$fasta_path),
      out_dir = file.path(output_dir, .data$Assembly_Base_ID),
      gff_path = file.path(.data$out_dir, paste0(.data$Assembly_Base_ID, ".gff")),
      log_path = file.path(log_dir, paste0(.data$Assembly_Base_ID, ".log"))
    ) |>
    dplyr::distinct(.data$Assembly_Base_ID, .keep_all = TRUE) |>
    dplyr::filter(!is.na(.data$Assembly_Base_ID), nzchar(.data$Assembly_Base_ID),
                  !is.na(.data$fasta_path), nzchar(.data$fasta_path))
}

log_pipeline_gff <- function(logger, fmt, ...) {
  txt <- sprintf(fmt, ...)
  if (is.null(logger)) {
    message(txt)
  } else {
    logger(txt)
  }
}

empty_gff_regen_summary <- function() {
  tibble::tibble(
    Assembly_Base_ID = character(),
    Participant_id = character(),
    tp_lab = character(),
    Assembler = character(),
    status = integer(),
    gff_created = logical(),
    gff_path = character(),
    log_path = character(),
    reason = character()
  )
}

regenerate_missing_gffs <- function(missing_df,
                                    output_dir = DIR_PROKKA,
                                    log_dir = file.path(DIR_LOGS, "prokka_missing_gffs"),
                                    summary_file = file.path(DIR_WGS, "pan", "regenerate_missing_gffs_summary.csv"),
                                    missing_fastas_file = file.path(DIR_WGS, "pan", "missing_gffs_no_fasta.csv"),
                                    cores = NULL,
                                    dry_run = identical(Sys.getenv("GFF_REGEN_DRY_RUN", "0"), "1"),
                                    limit = NULL,
                                    prokka_bin = NULL,
                                    logger = NULL) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(dirname(summary_file), recursive = TRUE, showWarnings = FALSE)

  jobs <- prepare_gff_regen_jobs(missing_df, output_dir = output_dir, log_dir = log_dir)

  if (is.null(limit)) {
    limit_raw <- Sys.getenv("GFF_REGEN_LIMIT", "0")
    limit <- suppressWarnings(as.integer(limit_raw))
  }

  missing_fastas <- jobs |>
    dplyr::filter(!file.exists(.data$fasta_path))
  if (nrow(missing_fastas) > 0) {
    readr::write_csv(missing_fastas, missing_fastas_file)
    log_pipeline_gff(
      logger,
      "WARNING: %d manifest row(s) point to missing FASTA files; wrote %s",
      nrow(missing_fastas),
      missing_fastas_file
    )
  }

  missing_fasta_results <- if (nrow(missing_fastas) > 0) {
    missing_fastas |>
      dplyr::transmute(
        Assembly_Base_ID = .data$Assembly_Base_ID,
        Participant_id = as.character(.data$Participant_id),
        tp_lab = as.character(.data$tp_lab),
        Assembler = as.character(.data$Assembler),
        status = NA_integer_,
        gff_created = FALSE,
        gff_path = .data$gff_path,
        log_path = .data$log_path,
        reason = "missing_fasta"
      )
  } else {
    empty_gff_regen_summary()
  }

  jobs <- jobs |>
    dplyr::filter(file.exists(.data$fasta_path)) |>
    dplyr::filter(!file.exists(.data$gff_path) | file.info(.data$gff_path)$size <= 0)

  if (!is.na(limit) && limit > 0) jobs <- utils::head(jobs, limit)

  log_pipeline_gff(logger, "Prokka jobs to run: %d", nrow(jobs))
  log_pipeline_gff(logger, "Output directory: %s", output_dir)
  log_pipeline_gff(logger, "Logs: %s", log_dir)

  if (nrow(jobs) == 0) {
    summary <- missing_fasta_results
    readr::write_csv(summary, summary_file)
    log_pipeline_gff(logger, "No missing GFFs need Prokka regeneration.")
    return(list(
      summary = summary,
      jobs = jobs,
      missing_fastas = missing_fastas,
      dry_run = dry_run,
      status = ifelse(nrow(missing_fasta_results) > 0, 1L, 0L),
      summary_file = summary_file
    ))
  }

  if (is.null(prokka_bin) || identical(prokka_bin, "")) {
    prokka_bin <- find_prokka_bin("prokka")
  }
  if (identical(prokka_bin, "")) {
    stop("Could not find prokka on PATH. Activate the yellow-wgs-x86 environment first.")
  }
  set_prokka_java_home()
  log_pipeline_gff(logger, "Prokka: %s", prokka_bin)

  if (isTRUE(dry_run)) {
    dry_summary <- jobs |>
      dplyr::transmute(
        Assembly_Base_ID = .data$Assembly_Base_ID,
        Participant_id = as.character(.data$Participant_id),
        tp_lab = as.character(.data$tp_lab),
        Assembler = as.character(.data$Assembler),
        status = NA_integer_,
        gff_created = file.exists(.data$gff_path) & file.info(.data$gff_path)$size > 0,
        gff_path = .data$gff_path,
        log_path = .data$log_path,
        reason = "dry_run"
      )
    summary <- dplyr::bind_rows(missing_fasta_results, dry_summary)
    readr::write_csv(summary, summary_file)
    log_pipeline_gff(logger, "Dry run requested; not running Prokka.")
    return(list(
      summary = summary,
      jobs = jobs,
      missing_fastas = missing_fastas,
      dry_run = TRUE,
      status = 0L,
      summary_file = summary_file
    ))
  }

  default_cores <- if (exists("CORES_USE")) min(CORES_USE, 6L) else 4L
  if (is.null(cores)) {
    cores_raw <- Sys.getenv("GFF_REGEN_CORES", unset = as.character(default_cores))
    if (is.na(cores_raw) || identical(cores_raw, "")) cores_raw <- as.character(default_cores)
    cores <- suppressWarnings(as.integer(cores_raw))
  }
  if (is.na(cores) || cores < 1) cores <- 1L
  cores <- min(cores, nrow(jobs))
  log_pipeline_gff(logger, "Running Prokka with %d parallel job(s), 1 CPU per job.", cores)

  run_one <- function(i) {
    row <- jobs[i, ]
    args <- c(
      "--outdir", row$out_dir,
      "--prefix", row$Assembly_Base_ID,
      "--force",
      "--compliant",
      "--centre", "X",
      "--cpus", "1",
      "--kingdom", "Bacteria",
      "--genus", "Escherichia",
      "--species", "coli",
      row$fasta_path
    )

    status <- suppressWarnings(system2(prokka_bin, args = args, stdout = row$log_path, stderr = row$log_path))
    status <- as.integer(if (is.null(status)) 0L else status)
    ok <- file.exists(row$gff_path) && file.info(row$gff_path)$size > 0
    reason <- if (ok) {
      "created"
    } else if (!identical(status, 0L)) {
      paste0("prokka_exit_", status)
    } else {
      "gff_missing_or_empty"
    }

    tibble::tibble(
      Assembly_Base_ID = row$Assembly_Base_ID,
      Participant_id = as.character(row$Participant_id),
      tp_lab = as.character(row$tp_lab),
      Assembler = as.character(row$Assembler),
      status = status,
      gff_created = ok,
      gff_path = row$gff_path,
      log_path = row$log_path,
      reason = reason
    )
  }

  res <- dplyr::bind_rows(parallel::mclapply(seq_len(nrow(jobs)), run_one, mc.cores = cores))
  summary <- dplyr::bind_rows(missing_fasta_results, res)
  readr::write_csv(summary, summary_file)

  log_pipeline_gff(logger, "Finished. GFFs created: %d/%d", sum(res$gff_created), nrow(res))
  log_pipeline_gff(logger, "Summary: %s", summary_file)

  if (any(!summary$gff_created)) {
    failed_file <- sub("\\.csv$", "_failed.csv", summary_file)
    readr::write_csv(summary |> dplyr::filter(!.data$gff_created), failed_file)
    log_pipeline_gff(logger, "Failures: %s", failed_file)
  }

  list(
    summary = summary,
    jobs = jobs,
    missing_fastas = missing_fastas,
    dry_run = FALSE,
    status = ifelse(any(!summary$gff_created), 1L, 0L),
    summary_file = summary_file
  )
}

canonical_vf_gene_cols <- function(vf_ready_names = NULL,
                                   vf_pa_file = FILE_VF_PA,
                                   required = TRUE) {
  if (!file.exists(vf_pa_file)) {
    if (required) stop("Missing canonical VF P/A matrix: ", vf_pa_file)
    return(character())
  }

  pa_names <- names(readr::read_csv(vf_pa_file, show_col_types = FALSE, n_max = 0))
  pa_meta_cols <- c("Participant_id", "tp_lab", "Episode_ID", "Event_type",
                    "Assembly_ID", "Isolate_ID")
  gene_cols <- setdiff(pa_names, pa_meta_cols)

  if (!is.null(vf_ready_names)) {
    missing_from_ready <- setdiff(gene_cols, vf_ready_names)
    if (length(missing_from_ready) > 0) {
      stop(
        "VF-ready dataset is missing ", length(missing_from_ready),
        " canonical VF gene column(s): ",
        paste(head(missing_from_ready, 10), collapse = ", "),
        if (length(missing_from_ready) > 10) ", ..." else ""
      )
    }
    gene_cols <- intersect(gene_cols, vf_ready_names)
  }

  gene_cols
}

append_denominator_summary <- function(df,
                                       script,
                                       stage,
                                       analysis_unit,
                                       input_file = NA_character_,
                                       notes = NA_character_,
                                       out_file = file.path(DIR_QC, "pipeline_denominator_summary.csv")) {
  dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)
  if (is.null(df)) df <- tibble::tibble()

  status_col <- if ("UTI_Status" %in% names(df)) "UTI_Status" else if ("Infection_Status" %in% names(df)) "Infection_Status" else if ("Status_Simple" %in% names(df)) "Status_Simple" else NA_character_
  tp <- if ("tp_lab" %in% names(df)) {
    normalise_timepoint_preserve_events(df$tp_lab)
  } else if ("Timepoint" %in% names(df)) {
    normalise_timepoint_preserve_events(df$Timepoint)
  } else {
    rep(NA_character_, nrow(df))
  }
  pid <- if ("Participant_id" %in% names(df)) safe_chr(df$Participant_id) else rep(NA_character_, nrow(df))
  assembly_cols <- intersect(c("Assembly_ID", "full_path", "fasta_path", "file_name"), names(df))
  n_assemblies <- if (length(assembly_cols) > 0 && nrow(df) > 0) {
    nrow(dplyr::distinct(df, dplyr::across(dplyr::all_of(assembly_cols[1]))))
  } else NA_integer_
  statuses <- if (!is.na(status_col)) safe_chr(df[[status_col]]) else rep(NA_character_, nrow(df))
  legacy_statuses <- if ("Infection_Status_old" %in% names(df)) {
    safe_chr(df$Infection_Status_old)
  } else if ("Infection_Status_legacy" %in% names(df)) {
    safe_chr(df$Infection_Status_legacy)
  } else if ("Infection_Status" %in% names(df)) {
    safe_chr(df$Infection_Status)
  } else {
    rep(NA_character_, nrow(df))
  }
  not_uti_subgroup <- if ("Not_UTI_subgroup" %in% names(df)) safe_chr(df$Not_UTI_subgroup) else rep(NA_character_, nrow(df))

  row <- tibble::tibble(
    timestamp = as.character(Sys.time()),
    script = script,
    stage = stage,
    analysis_unit = analysis_unit,
    input_file = paste(input_file, collapse = "; "),
    n_rows = nrow(df),
    n_unique_participants = dplyr::n_distinct(pid[!is.na(pid)]),
    n_unique_participant_timepoints = dplyr::n_distinct(paste(pid, tp, sep = "|")[!is.na(pid) & !is.na(tp)]),
    n_assemblies = n_assemblies,
    primary_status_col = status_col,
    n_ASB = sum(legacy_statuses == "ASB", na.rm = TRUE),
    n_UTI = sum(statuses == "UTI", na.rm = TRUE),
    n_Not_UTI = sum(statuses == "Not_UTI", na.rm = TRUE),
    n_legacy_UTI = sum(legacy_statuses == "UTI", na.rm = TRUE),
    n_Negative = sum(legacy_statuses == "Negative", na.rm = TRUE),
    n_bacteriuria_not_UTI = sum(not_uti_subgroup == "bacteriuria_not_UTI", na.rm = TRUE),
    n_culture_negative_or_below_threshold = sum(not_uti_subgroup == "culture_negative_or_below_threshold", na.rm = TRUE),
    n_unknown_or_indeterminate = sum(not_uti_subgroup == "unknown_or_indeterminate", na.rm = TRUE),
    n_missing_status = if (!is.na(status_col)) sum(is.na(statuses) | statuses == "", na.rm = TRUE) else NA_integer_,
    notes = notes
  )

  existing <- if (file.exists(out_file)) readr::read_csv(out_file, show_col_types = FALSE) else tibble::tibble()
  if ("timestamp" %in% names(existing)) existing$timestamp <- as.character(existing$timestamp)
  if (nrow(existing) > 0 && all(c("script", "stage", "analysis_unit") %in% names(existing))) {
    existing <- existing |>
      dplyr::filter(!(.data$script == .env$script & .data$stage == .env$stage & .data$analysis_unit == .env$analysis_unit))
  }
  if (nrow(existing) == 0) existing <- tibble::tibble()
  readr::write_csv(dplyr::bind_rows(existing, row), out_file)
  invisible(row)
}

select_canonical_assemblies <- function(qc_df) {
  if (!"QC_PASS" %in% names(qc_df)) stop("select_canonical_assemblies requires QC_PASS column")
  if (!"full_path" %in% names(qc_df)) qc_df$full_path <- NA_character_
  if (!"file_name" %in% names(qc_df)) qc_df$file_name <- basename(qc_df$full_path)
  if (!"Participant_id" %in% names(qc_df)) qc_df$Participant_id <- NA_character_
  if (!"Timepoint" %in% names(qc_df) && !"tp_lab" %in% names(qc_df)) qc_df$Timepoint <- NA_character_
  if (!"N50" %in% names(qc_df)) qc_df$N50 <- NA_real_
  if (!"n_contigs" %in% names(qc_df)) qc_df$n_contigs <- NA_real_
  out <- qc_df |>
    dplyr::mutate(
      Participant_id = safe_chr(.data$Participant_id),
      tp_lab = if ("tp_lab" %in% names(qc_df)) normalise_timepoint_preserve_events(.data$tp_lab) else normalise_timepoint_preserve_events(.data$Timepoint),
      Assembly_ID = if ("Assembly_ID" %in% names(qc_df)) .data$Assembly_ID else tools::file_path_sans_ext(basename(ifelse(!is.na(.data$full_path), .data$full_path, .data$file_name))),
      Assembler = if ("Assembler" %in% names(qc_df)) .data$Assembler else if ("assembler" %in% names(qc_df)) .data$assembler else detect_assembler(ifelse(!is.na(.data$file_name), .data$file_name, .data$full_path)),
      file_exists = !is.na(.data$full_path) & file.exists(.data$full_path),
      usable_fasta = .data$file_exists & .data$QC_PASS,
      qc_priority = ifelse(.data$QC_PASS, 0L, 1L),
      assembler_priority = dplyr::case_when(
        .data$Assembler == "longcycler" ~ 1L,
        .data$Assembler == "flye" ~ 2L,
        TRUE ~ 3L
      ),
      n50_priority = -1 * suppressWarnings(as.numeric(.data$N50)),
      contig_priority = suppressWarnings(as.numeric(.data$n_contigs))
    ) |>
    dplyr::group_by(.data$Participant_id, .data$tp_lab) |>
    dplyr::arrange(.data$qc_priority, .data$assembler_priority, .data$n50_priority, .data$contig_priority, .data$Assembly_ID, .by_group = TRUE) |>
    dplyr::mutate(
      canonical_rank = dplyr::row_number(),
      selected_canonical = .data$canonical_rank == 1L & .data$QC_PASS,
      canonical_reason = dplyr::case_when(
        .data$selected_canonical & .data$Assembler == "longcycler" ~ "QC PASS; preferred longcycler",
        .data$selected_canonical & .data$Assembler == "flye" ~ "QC PASS; flye selected because preferred assembly unavailable/failing",
        .data$selected_canonical ~ "QC PASS; selected by deterministic fallback",
        any(.data$QC_PASS, na.rm = TRUE) ~ "Alternative assembly not selected",
        TRUE ~ "No QC PASS assembly for participant-timepoint"
      )
    ) |>
    dplyr::ungroup()
  out
}

write_uti_attrition_outputs <- function(status_file = file.path(DIR_CLINICAL, "status_map.csv"),
                                        out_episode = file.path(DIR_QC, "uti_attrition_episode_level.csv"),
                                        out_stage = file.path(DIR_QC, "uti_attrition_by_pipeline_stage.csv"),
                                        out_md = file.path(DIR_RESULTS, "summary", "uti_attrition_summary.md")) {
  if (!file.exists(status_file)) return(invisible(NULL))
  dir.create(dirname(out_episode), recursive = TRUE, showWarnings = FALSE)
  dir.create(dirname(out_md), recursive = TRUE, showWarnings = FALSE)

  status <- readr::read_csv(status_file, show_col_types = FALSE)
  status$Participant_id <- safe_chr(status$Participant_id)
  status$tp_lab <- if ("tp_lab" %in% names(status)) {
    normalise_timepoint_preserve_events(status$tp_lab)
  } else {
    normalise_timepoint_preserve_events(status$Timepoint)
  }
  if (!"Event_type" %in% names(status)) status$Event_type <- episode_event_type(status$tp_lab)
  if (!"Episode_ID" %in% names(status)) {
    status$Episode_ID <- build_episode_id(
      status,
      timepoint_col = if ("Timepoint" %in% names(status)) "Timepoint" else "tp_lab",
      event_col = "Event_type",
      date_col = if ("Collection_Date" %in% names(status)) "Collection_Date" else NULL
    )
  }
  if (!"Collection_Date" %in% names(status)) status$Collection_Date <- NA_character_
  if (!"UTI_Status" %in% names(status)) {
    stop("write_uti_attrition_outputs: status map lacks UTI_Status; refusing legacy status fallback.")
  }
  uti <- status |> dplyr::filter(.data$UTI_Status == "UTI")
  if (nrow(uti) == 0) return(invisible(NULL))

  key <- function(df, tp_col = "tp_lab") paste(safe_chr(df$Participant_id), normalise_timepoint_preserve_events(df[[tp_col]]), sep = "|")
  uti$key <- key(uti)
  uti$episode_key <- paste(safe_chr(uti$Participant_id), safe_chr(uti$Episode_ID), sep = "|")

  bridge <- if (file.exists(file.path(DIR_QC, "uricult_bridge_audit.csv"))) {
    readr::read_csv(file.path(DIR_QC, "uricult_bridge_audit.csv"), show_col_types = FALSE) |>
      dplyr::filter(.data$selected %in% TRUE) |>
      dplyr::transmute(
        episode_key = paste(safe_chr(.data$Participant_id_clinical), safe_chr(.data$Episode_ID_clinical), sep = "|"),
        bridged_wgs_key = paste(safe_chr(.data$Participant_id_clinical), normalise_timepoint_preserve_events(.data$mapped_tp_lab), sep = "|"),
        bridged_wgs_tp_lab = normalise_timepoint_preserve_events(.data$mapped_tp_lab)
      ) |>
      dplyr::distinct(.data$episode_key, .keep_all = TRUE)
  } else {
    tibble::tibble(episode_key = character(), bridged_wgs_key = character(), bridged_wgs_tp_lab = character())
  }
  bridge_match <- match(uti$episode_key, bridge$episode_key)
  uti$linked_via_uricult_bridge <- !is.na(bridge_match)
  uti$wgs_key <- ifelse(uti$linked_via_uricult_bridge, bridge$bridged_wgs_key[bridge_match], uti$key)
  uti$wgs_tp_lab <- ifelse(uti$linked_via_uricult_bridge, bridge$bridged_wgs_tp_lab[bridge_match], uti$tp_lab)

  meta <- if (file.exists(FILE_METADATA)) {
    m <- readr::read_csv(FILE_METADATA, show_col_types = FALSE)
    m$Participant_id <- safe_chr(m$Participant_id)
    m$tp_lab <- if ("tp_lab" %in% names(m)) normalise_timepoint_preserve_events(m$tp_lab) else normalise_timepoint_preserve_events(m$Timepoint)
    if (!"full_path" %in% names(m) && "fasta_path" %in% names(m)) m$full_path <- m$fasta_path
    m$key <- paste(m$Participant_id, m$tp_lab, sep = "|")
    m$has_usable <- if ("usable_fasta" %in% names(m)) m$usable_fasta else (!is.na(m$full_path) & file.exists(m$full_path))
    m
  } else tibble::tibble(key = character(), has_usable = logical())
  canonical <- if (file.exists(file.path(DIR_QC, "canonical_assembly_selection.csv"))) {
    readr::read_csv(file.path(DIR_QC, "canonical_assembly_selection.csv"), show_col_types = FALSE) |>
      dplyr::mutate(Participant_id = safe_chr(.data$Participant_id),
                    tp_lab = normalise_timepoint_preserve_events(.data$tp_lab),
                    key = paste(.data$Participant_id, .data$tp_lab, sep = "|"))
  } else tibble::tibble(key = character(), selected_canonical = logical(), QC_PASS = logical())
  qc <- if (file.exists(file.path(DIR_WGS, "qc_summary.csv"))) {
    q <- readr::read_csv(file.path(DIR_WGS, "qc_summary.csv"), show_col_types = FALSE)
    q$Participant_id <- safe_chr(q$Participant_id)
    q$tp_lab <- if ("tp_lab" %in% names(q)) normalise_timepoint_preserve_events(q$tp_lab) else normalise_timepoint_preserve_events(q$Timepoint)
    q$key <- paste(q$Participant_id, q$tp_lab, sep = "|")
    q
  } else tibble::tibble(key = character(), QC_PASS = logical())
  vf_pa <- if (file.exists(FILE_VF_PA)) {
    readr::read_csv(FILE_VF_PA, show_col_types = FALSE) |>
      dplyr::mutate(Participant_id = safe_chr(.data$Participant_id),
                    tp_lab = normalise_timepoint_preserve_events(.data$tp_lab),
                    key = paste(.data$Participant_id, .data$tp_lab, sep = "|"))
  } else tibble::tibble(key = character())
  vf_ready <- if (file.exists(file.path(DIR_VF, "vf_analysis_ready.csv"))) {
    vr <- readr::read_csv(file.path(DIR_VF, "vf_analysis_ready.csv"), show_col_types = FALSE)
    if (!"Episode_ID" %in% names(vr)) vr$Episode_ID <- NA_character_
    vr |>
      dplyr::mutate(Participant_id = safe_chr(.data$Participant_id),
                    tp_lab = normalise_timepoint_preserve_events(.data$tp_lab),
                    key = paste(.data$Participant_id, .data$tp_lab, sep = "|"),
                    episode_key = paste(.data$Participant_id, safe_chr(.data$Episode_ID), sep = "|"))
  } else tibble::tibble(key = character(), Infection_Status = character(), UTI_Status = character())
  if (!"UTI_Status" %in% names(vf_ready)) {
    vf_ready$UTI_Status <- NA_character_
  }
  mlst <- if (file.exists(FILE_MLST_CANONICAL)) {
    ml <- readr::read_csv(FILE_MLST_CANONICAL, show_col_types = FALSE)
    ml$Participant_id <- safe_chr(ml$Participant_id)
    ml$tp_lab <- if ("tp_lab" %in% names(ml)) normalise_timepoint_preserve_events(ml$tp_lab) else normalise_timepoint_preserve_events(ml$Timepoint)
    ml$key <- paste(ml$Participant_id, ml$tp_lab, sep = "|")
    ml
  } else tibble::tibble(key = character())
  model_ds <- if (file.exists(file.path(DIR_MODELS, "model_dataset_denominator.csv"))) {
    readr::read_csv(file.path(DIR_MODELS, "model_dataset_denominator.csv"), show_col_types = FALSE)
  } else tibble::tibble()
  model_keys <- if (all(c("Participant_id", "Timepoint") %in% names(model_ds))) {
    paste(safe_chr(model_ds$Participant_id), normalise_timepoint_preserve_events(model_ds$Timepoint), sep = "|")
  } else character()
  model_episode_keys <- if (all(c("Participant_id", "Episode_ID") %in% names(model_ds))) {
    paste(safe_chr(model_ds$Participant_id), safe_chr(model_ds$Episode_ID), sep = "|")
  } else character()
  trans <- if (file.exists(file.path(DIR_VF, "vf_transition_case_index.csv"))) {
    readr::read_csv(file.path(DIR_VF, "vf_transition_case_index.csv"), show_col_types = FALSE)
  } else tibble::tibble()
  trans_keys <- unique(c(
    if (all(c("Participant_id", "from_tp") %in% names(trans))) paste(safe_chr(trans$Participant_id), normalise_timepoint_preserve_events(trans$from_tp), sep = "|") else character(),
    if (all(c("Participant_id", "to_tp") %in% names(trans))) paste(safe_chr(trans$Participant_id), normalise_timepoint_preserve_events(trans$to_tp), sep = "|") else character()
  ))
  trans_episode_keys <- unique(c(
    if (all(c("Participant_id", "from_episode_id") %in% names(trans))) paste(safe_chr(trans$Participant_id), safe_chr(trans$from_episode_id), sep = "|") else character(),
    if (all(c("Participant_id", "to_episode_id") %in% names(trans))) paste(safe_chr(trans$Participant_id), safe_chr(trans$to_episode_id), sep = "|") else character()
  ))

  episode <- uti |>
    dplyr::transmute(
      Participant_id, tp_lab, wgs_tp_lab, Episode_ID, Collection_Date,
      Event_type, Uricult_flag = .data$Event_type == "Uricult",
      linked_via_uricult_bridge,
      in_assembly_metadata = .data$wgs_key %in% unique(meta$key),
      has_usable_fasta = .data$wgs_key %in% unique(meta$key[meta$has_usable %in% TRUE]),
      has_selected_canonical = .data$wgs_key %in% unique(canonical$key[canonical$selected_canonical %in% TRUE]),
      qc_present = .data$wgs_key %in% unique(qc$key),
      qc_pass = .data$wgs_key %in% unique(qc$key[qc$QC_PASS %in% TRUE]),
      in_vf_pa_all = .data$wgs_key %in% unique(vf_pa$key),
      in_vf_analysis_ready = .data$wgs_key %in% unique(vf_ready$key[vf_ready$UTI_Status == "UTI"]) |
        .data$episode_key %in% unique(vf_ready$episode_key[vf_ready$UTI_Status == "UTI"]),
      in_mlst = .data$wgs_key %in% unique(mlst$key),
      in_model_14_dataset = .data$wgs_key %in% model_keys | .data$episode_key %in% model_episode_keys,
      in_transition_case_studies = .data$wgs_key %in% trans_keys | .data$episode_key %in% trans_episode_keys
    ) |>
    dplyr::mutate(
      reason_lost_primary = dplyr::case_when(
        !.data$in_assembly_metadata ~ "No assembly metadata/WGS link",
        !.data$has_usable_fasta ~ "No usable FASTA",
        !.data$qc_present ~ "Missing WGS QC row",
        !.data$qc_pass ~ "QC fail",
        !.data$in_vf_pa_all ~ "Missing VF result or key mismatch",
        !.data$in_vf_analysis_ready ~ "Missing VF-ready UTI row",
        TRUE ~ "Retained in VF-ready UTI set"
      ),
      reason_lost_detailed = dplyr::case_when(
        .data$linked_via_uricult_bridge & .data$in_vf_analysis_ready ~ "Uricult UTI episode linked through audited clinical Episode_ID to WGS UTI-N row",
        .data$Uricult_flag & !.data$in_assembly_metadata ~ "Uricult UTI episode has no matching isolate/assembly metadata; cannot infer WGS without trusted mapping",
        !.data$in_vf_pa_all & .data$qc_pass ~ "QC-passing WGS exists but no VF P/A row under the preserved episode key; rerun script 02/22 after key repair",
        !.data$in_vf_analysis_ready & .data$in_vf_pa_all ~ "VF P/A exists but status join did not retain UTI label; check Episode_ID/key join",
        TRUE ~ .data$reason_lost_primary
      )
    )

  readr::write_csv(episode, out_episode)
  stage <- tibble::tibble(
    stage = c("clinical_UTI", "in_assembly_metadata", "has_usable_fasta", "qc_present", "qc_pass",
              "in_vf_pa_all", "in_vf_analysis_ready", "in_mlst", "in_model_14_dataset", "in_transition_case_studies"),
    n_UTI = c(nrow(episode), sum(episode$in_assembly_metadata), sum(episode$has_usable_fasta),
              sum(episode$qc_present), sum(episode$qc_pass), sum(episode$in_vf_pa_all),
              sum(episode$in_vf_analysis_ready), sum(episode$in_mlst),
              sum(episode$in_model_14_dataset), sum(episode$in_transition_case_studies))
  )
  readr::write_csv(stage, out_stage)
  md <- c(
    "# UTI Attrition Summary",
    "",
    sprintf("Generated: %s", format(Sys.time())),
    "",
    sprintf("- Clinical UTI episodes: **%d**", nrow(episode)),
    sprintf("- UTI with assembly metadata: **%d**", sum(episode$in_assembly_metadata)),
    sprintf("- UTI with QC PASS assembly: **%d**", sum(episode$qc_pass)),
    sprintf("- UTI retained in VF-ready dataset: **%d**", sum(episode$in_vf_analysis_ready)),
    "",
    "## Primary Loss Reasons",
    "",
    paste0("- ", names(table(episode$reason_lost_primary)), ": ", as.integer(table(episode$reason_lost_primary)))
  )
  writeLines(md, out_md)
  invisible(episode)
}

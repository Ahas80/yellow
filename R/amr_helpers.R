# ==============================================================================
# Genomic AMR helpers for 29_vf_amr_combined_profile.R
# ==============================================================================
#
# This module intentionally contains no top-level analysis. Script 29 is the
# numbered owner and calls run_genomic_amr_analysis(). The public pure helpers
# are also used by tests.
# ==============================================================================

AMR_SCHEMA_VERSION <- "genomic_amr_v2"
AMR_EXPECTED_EPISODES <- 532L
AMR_EXPECTED_RESIDENTS <- 161L
AMR_EXPECTED_ADJACENT <- 371L
AMR_EXPECTED_FOCUSED <- 9L
AMR_BOOTSTRAP_SEED <- 20260712L
AMR_AMRFINDER_VERSION <- "4.2.7"
AMR_AMRFINDER_DB_VERSION <- "2026-05-15.1"
AMR_RESFINDER_VERSION <- "4.7.2"
AMR_RESFINDER_DB_COMMIT <- "eecf0aa207594fe6d51badf808473de62b28cb06"
AMR_POINTFINDER_DB_COMMIT <- "44ce624a806c6d2b70f7e39841a5f9cb4d9010aa"

amr_normalise_symbol <- function(x) {
  x <- trimws(as.character(x))
  x[is.na(x)] <- ""
  x <- sub("_[0-9]+$", "", x)
  x <- sub("^mdfA$", "mdf(A)", x, ignore.case = TRUE)
  x <- sub("^tetA$", "tet(A)", x, ignore.case = TRUE)
  x <- sub("^tetB$", "tet(B)", x, ignore.case = TRUE)
  x <- sub("^tetC$", "tet(C)", x, ignore.case = TRUE)
  x <- sub("^tetD$", "tet(D)", x, ignore.case = TRUE)
  x
}

amr_normalise_mutation <- function(symbol) {
  symbol <- trimws(as.character(symbol))
  symbol[is.na(symbol)] <- ""
  symbol <- sub("^([^_:[:space:]]+)_([A-Za-z*][0-9]+[A-Za-z*])$",
                "\\1:\\2", symbol)
  symbol
}

amr_gene_family <- function(x) {
  x <- amr_normalise_symbol(x)
  out <- x
  out <- sub("^(bla[A-Za-z]+)[-_]?[0-9].*$", "\\1", out)
  out <- sub("^(dfrA)[0-9]+.*$", "\\1", out, ignore.case = TRUE)
  out <- sub("^(sul)[0-9]+.*$", "\\1", out, ignore.case = TRUE)
  out <- sub("^((?:aac|aad|aph|ant)\\([^)]*\\)).*$", "\\1", out, ignore.case = TRUE)
  out <- sub("^((?:qnr)[A-Za-z]+)[0-9]*$", "\\1", out, ignore.case = TRUE)
  out
}

amr_is_mdfa <- function(x) {
  tolower(gsub("[()_ -]", "", amr_normalise_symbol(x))) == "mdfa"
}

amr_is_background <- function(x) {
  amr_is_mdfa(x)
}

amr_infer_class <- function(symbol, supplied = NA_character_) {
  supplied <- trimws(as.character(supplied))
  symbol <- tolower(amr_normalise_symbol(symbol))
  supplied_lower <- tolower(supplied)
  supplied_class <- dplyr::case_when(
    grepl("beta.?lact|penicillin|cephalo|carbapenem|aztreonam",
          supplied_lower) ~ "Beta-lactams",
    grepl("aminoglycos|aminocyclitol|gentamicin|tobramycin|amikacin|streptomycin",
          supplied_lower) ~ "Aminoglycosides",
    grepl("tetracycl", supplied_lower) ~ "Tetracyclines",
    grepl("sulfon|sulphon", supplied_lower) ~ "Sulfonamides",
    grepl("trimethoprim", supplied_lower) ~ "Trimethoprim",
    grepl("quinolone|fluoroquinolone|ciprofloxacin|nalidixic",
          supplied_lower) ~ "Fluoroquinolones",
    grepl("phenicol|chloramphenicol|florfenicol", supplied_lower) ~ "Phenicols",
    grepl("macrolide|lincosamide|streptogramin|azithromycin",
          supplied_lower) ~ "Macrolides/lincosamides/streptogramins",
    grepl("fosfomycin", supplied_lower) ~ "Fosfomycin",
    grepl("fosmidomycin", supplied_lower) ~ "Fosmidomycin",
    grepl("polymyxin|colistin", supplied_lower) ~ "Polymyxins",
    grepl("rifamp", supplied_lower) ~ "Rifamycins",
    grepl("glycopeptide|vancomycin|teicoplanin", supplied_lower) ~ "Glycopeptides",
    grepl("multidrug|efflux", supplied_lower) ~ "Multidrug efflux",
    TRUE ~ NA_character_
  )
  inferred <- dplyr::case_when(
    grepl("^bla", symbol) ~ "Beta-lactams",
    grepl("^(aac|aad|aph|ant)", symbol) ~ "Aminoglycosides",
    grepl("^tet", symbol) ~ "Tetracyclines",
    grepl("^sul", symbol) ~ "Sulfonamides",
    grepl("^dfr", symbol) ~ "Trimethoprim",
    grepl("^(qnr|gyr[ab]|par[ce])", symbol) ~ "Fluoroquinolones",
    grepl("^(cat|cml|flo)", symbol) ~ "Phenicols",
    grepl("^(mph|erm|mef)", symbol) ~ "Macrolides",
    grepl("^fos", symbol) ~ "Fosfomycin",
    grepl("^mcr", symbol) ~ "Polymyxins",
    grepl("^mdf", symbol) ~ "Multidrug efflux",
    TRUE ~ "Other/unspecified"
  )
  ifelse(!is.na(supplied_class), supplied_class, inferred)
}

amr_jaccard <- function(a, b, empty_empty = 1) {
  a <- sort(unique(as.character(a[!is.na(a) & nzchar(a)])))
  b <- sort(unique(as.character(b[!is.na(b) & nzchar(b)])))
  union_set <- union(a, b)
  if (!length(union_set)) return(as.numeric(empty_empty))
  length(intersect(a, b)) / length(union_set)
}

amr_cache_key <- function(fasta_sha256, annotation_sha256 = NA_character_,
                          tool, tool_version, database_fingerprint, parameters) {
  payload <- list(
    schema = AMR_SCHEMA_VERSION,
    fasta_sha256 = tolower(as.character(fasta_sha256)),
    annotation_sha256 = tolower(as.character(annotation_sha256)),
    tool = as.character(tool),
    tool_version = as.character(tool_version),
    database_fingerprint = as.character(database_fingerprint),
    parameters = parameters
  )
  substr(digest::digest(payload, algo = "sha256", serialize = TRUE), 1L, 24L)
}

amr_sequence_content_hash <- function(path) {
  dna <- Biostrings::readDNAStringSet(path, format = "fasta")
  sequence_hashes <- vapply(
    as.character(dna),
    function(sequence) digest::digest(toupper(sequence), algo = "sha256", serialize = FALSE),
    character(1)
  )
  digest::digest(sort(unname(sequence_hashes)), algo = "sha256", serialize = TRUE)
}

amr_pick_column <- function(df, candidates, default = NA_character_) {
  hit <- candidates[candidates %in% names(df)]
  if (!length(hit)) return(rep(default, nrow(df)))
  as.character(df[[hit[[1L]]]])
}

amr_safe_number <- function(x) {
  suppressWarnings(readr::parse_number(as.character(x)))
}

amr_collapse <- function(x) {
  x <- sort(unique(as.character(x[!is.na(x) & nzchar(x)])))
  paste(x, collapse = ";")
}

amr_split <- function(x) {
  if (!length(x) || is.na(x) || !nzchar(x)) return(character())
  unique(trimws(strsplit(as.character(x), ";", fixed = TRUE)[[1L]]))
}

amr_quote_args <- function(args) {
  vapply(as.character(args), shQuote, character(1), type = "sh")
}

amr_command_version <- function(command, args) {
  out <- tryCatch(
    suppressWarnings(system2(
      command, amr_quote_args(args), stdout = TRUE, stderr = TRUE
    )),
    error = function(e) structure(character(), status = 127L)
  )
  status <- attr(out, "status")
  if (is.null(status)) status <- 0L
  if (!identical(as.integer(status), 0L)) return("")
  trimws(paste(out, collapse = "\n"))
}

amr_git_head <- function(path) {
  amr_command_version("git", c("-C", path, "rev-parse", "HEAD"))
}

amr_sha256_file <- function(path) {
  unname(digest::digest(path, algo = "sha256", file = TRUE))
}

amr_atomic_write_csv <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  temporary <- tempfile(pattern = paste0(".", basename(path), "."),
                        tmpdir = dirname(path), fileext = ".tmp")
  on.exit(unlink(temporary), add = TRUE)
  readr::write_csv(x, temporary, na = "")
  if (!file.rename(temporary, path)) stop("Could not publish ", path, call. = FALSE)
  invisible(path)
}

amr_atomic_write_lines <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  temporary <- tempfile(pattern = paste0(".", basename(path), "."),
                        tmpdir = dirname(path), fileext = ".tmp")
  on.exit(unlink(temporary), add = TRUE)
  writeLines(as.character(x), temporary, useBytes = TRUE)
  if (!file.rename(temporary, path)) stop("Could not publish ", path, call. = FALSE)
  invisible(path)
}

amr_directory_fingerprint <- function(path, include = NULL) {
  if (!dir.exists(path)) return("")
  files <- if (is.null(include)) {
    list.files(path, recursive = TRUE, full.names = TRUE, include.dirs = FALSE)
  } else {
    file.path(path, include)
  }
  files <- sort(files[file.exists(files)])
  if (!length(files)) return("")
  manifest <- tibble::tibble(
    relative_path = sub(
      paste0("^", gsub("([][{}()+*^$|\\\\?.])", "\\\\\\1",
                       normalizePath(path, winslash = "/", mustWork = TRUE)),
             "/?"),
      "",
      normalizePath(files, winslash = "/", mustWork = TRUE)
    ),
    bytes = unname(file.info(files)$size),
    sha256 = vapply(files, amr_sha256_file, character(1))
  )
  digest::digest(manifest, algo = "sha256", serialize = TRUE)
}

amr_resolve_runtime <- function(root) {
  prefix <- Sys.getenv("AMR_RUNTIME_PREFIX", file.path(root, "data", "amr_runtime", "env"))
  db_root <- Sys.getenv("AMR_DATABASE_ROOT", file.path(root, "data", "amr_runtime", "databases"))
  abricate <- Sys.getenv("ABRICATE_BIN", Sys.which("abricate"))
  amrfinder <- Sys.getenv("AMRFINDER_BIN", file.path(prefix, "bin", "amrfinder"))
  python <- Sys.getenv("RESFINDER_PYTHON", file.path(prefix, "bin", "python"))
  kma <- Sys.getenv("KMA_BIN", file.path(prefix, "bin", "kma"))
  res_db <- Sys.getenv("RESFINDER_DB", file.path(db_root, "resfinder_db"))
  point_db <- Sys.getenv("POINTFINDER_DB", file.path(db_root, "pointfinder_db"))
  amrfinder_db_version <- Sys.getenv(
    "AMRFINDER_DB_VERSION", AMR_AMRFINDER_DB_VERSION
  )
  amrfinder_db_root <- Sys.getenv(
    "AMRFINDER_DB_ROOT", file.path(db_root, "amrfinderplus")
  )
  amrfinder_db <- Sys.getenv(
    "AMRFINDER_DB", file.path(amrfinder_db_root, amrfinder_db_version)
  )
  required <- c(abricate = abricate, amrfinder = amrfinder, python = python,
                kma = kma, amrfinder_db = amrfinder_db,
                resfinder_db = res_db, pointfinder_db = point_db)
  missing <- names(required)[!nzchar(required) | !file.exists(required)]
  if (length(missing)) {
    stop(
      "Pinned genomic-AMR runtime is incomplete (missing: ",
      paste(missing, collapse = ", "),
      "). Run scripts/setup_amr_runtime.sh, then rerun script 29.",
      call. = FALSE
    )
  }
  list(
    prefix = prefix, abricate = abricate, amrfinder = amrfinder,
    python = python, kma = kma, amrfinder_db = amrfinder_db,
    amrfinder_db_root = amrfinder_db_root,
    amrfinder_db_version = amrfinder_db_version,
    resfinder_db = res_db, pointfinder_db = point_db
  )
}

amr_runtime_provenance <- function(runtime) {
  abricate_version <- amr_command_version(runtime$abricate, "--version")
  amrfinder_version <- amr_command_version(runtime$amrfinder, "--version")
  resfinder_version <- amr_command_version(
    runtime$python,
    c("-c", "import importlib.metadata as m; print(m.version('resfinder'))")
  )
  abricate_help <- amr_command_version(runtime$abricate, "--help")
  datadir_line <- grep(
    "--datadir", strsplit(abricate_help, "\n", fixed = TRUE)[[1L]],
    value = TRUE
  )
  datadir <- if (length(datadir_line)) {
    sub(".*Databases folder \\[([^]]+)\\].*", "\\1", datadir_line[[1L]])
  } else {
    ""
  }
  if (!dir.exists(datadir)) {
    candidates <- Sys.glob(
      file.path(dirname(dirname(runtime$abricate)), "Cellar", "abricate",
                "*", "libexec", "db")
    )
    candidates <- candidates[dir.exists(candidates)]
    datadir <- if (length(candidates)) candidates[[1L]] else ""
  }
  abr_sequence <- file.path(datadir, "resfinder", "sequences")
  if (!file.exists(abr_sequence)) {
    stop("Could not fingerprint the ABRicate ResFinder database: ", abr_sequence,
         call. = FALSE)
  }
  amrfinder_db <- amr_command_version(
    runtime$amrfinder, c("-d", runtime$amrfinder_db, "-V")
  )
  amrfinder_db_hash <- amr_directory_fingerprint(
    runtime$amrfinder_db,
    include = c(
      "database_format_version.txt", "AMRProt.fa", "AMR_CDS.fa",
      "AMRProt-mutation.tsv", "fam.tsv", "taxgroup.tsv"
    )
  )
  resfinder_commit <- amr_git_head(runtime$resfinder_db)
  pointfinder_commit <- amr_git_head(runtime$pointfinder_db)
  tibble::tibble(
    component = c(
      "abricate", "abricate_resfinder_database", "amrfinderplus",
      "amrfinderplus_database", "resfinder", "resfinder_database",
      "pointfinder_database", "kma"
    ),
    version_or_commit = c(
      abricate_version,
      amr_sha256_file(abr_sequence),
      amrfinder_version,
      paste(
        paste0("version=", runtime$amrfinder_db_version),
        paste0("fingerprint=", amrfinder_db_hash),
        gsub("[\r\n]+", " | ", amrfinder_db),
        sep = ";"
      ),
      resfinder_version,
      resfinder_commit,
      pointfinder_commit,
      amr_command_version(runtime$kma, "-v")
    ),
    path = c(
      runtime$abricate, abr_sequence, runtime$amrfinder, runtime$amrfinder_db,
      runtime$python, runtime$resfinder_db, runtime$pointfinder_db, runtime$kma
    )
  )
}

amr_validate_manifest <- function(root) {
  assembly_path <- file.path(root, "results", "qc", "analysis_assembly_manifest.csv")
  annotation_path <- file.path(root, "results", "wgs", "pan", "panaroo_input_manifest.csv")
  if (!file.exists(assembly_path) || !file.exists(annotation_path)) {
    stop("Selected assembly or Panaroo/Prokka manifest is missing.", call. = FALSE)
  }
  assembly <- readr::read_csv(assembly_path, show_col_types = FALSE) |>
    dplyr::mutate(
      Participant_id = as.character(.data$Participant_id),
      tp_lab = normalise_timepoint_preserve_events(.data$tp_lab),
      fasta_path = normalizePath(.data$full_path, winslash = "/", mustWork = TRUE),
      fasta_sha256 = tolower(as.character(.data$fasta_sha256))
    ) |>
    dplyr::select(
      dplyr::any_of(c("Participant_id", "tp_lab", "Episode_ID", "Event_type",
                      "UTI_Status", "Assembly_ID", "Assembly_Base_ID", "Isolate_ID")),
      "fasta_path", "fasta_sha256"
    )
  annotation <- readr::read_csv(annotation_path, show_col_types = FALSE) |>
    dplyr::mutate(
      Participant_id = as.character(.data$Participant_id),
      tp_lab = normalise_timepoint_preserve_events(.data$tp_lab),
      fasta_path = normalizePath(.data$fasta_path, winslash = "/", mustWork = TRUE),
      fasta_sha256 = tolower(as.character(.data$fasta_sha256)),
      gff_path = normalizePath(.data$gff_path, winslash = "/", mustWork = TRUE),
      gff_sha256 = tolower(as.character(.data$gff_sha256))
    ) |>
    dplyr::select(
      "Participant_id", "tp_lab", "fasta_path", "fasta_sha256",
      "gff_path", "gff_sha256"
    )
  manifest <- assembly |>
    dplyr::inner_join(
      annotation,
      by = c("Participant_id", "tp_lab", "fasta_path", "fasta_sha256"),
      relationship = "one-to-one"
    )
  if (nrow(manifest) != AMR_EXPECTED_EPISODES ||
      dplyr::n_distinct(manifest$Participant_id) != AMR_EXPECTED_RESIDENTS ||
      anyDuplicated(paste(manifest$Participant_id, manifest$tp_lab, sep = "|"))) {
    stop("AMR input contract is not exactly 532 unique episodes from 161 residents.",
         call. = FALSE)
  }
  observed_fasta <- vapply(manifest$fasta_path, amr_sha256_file, character(1))
  observed_gff <- vapply(manifest$gff_path, amr_sha256_file, character(1))
  if (!all(tolower(observed_fasta) == manifest$fasta_sha256) ||
      !all(tolower(observed_gff) == manifest$gff_sha256)) {
    stop("Selected FASTA or Prokka GFF SHA-256 provenance changed.", call. = FALSE)
  }
  stem <- sub("\\.gff$", "", manifest$gff_path, ignore.case = TRUE)
  # Panaroo may use a filtered/prefixed derivative of the matched Prokka GFF.
  # FAA/FFN/FNA remain in the same Prokka bundle under the original stem.
  annotation_stem <- sub(
    "(\\.min270)?(\\.prefixed)?$", "", stem, ignore.case = TRUE
  )
  manifest <- manifest |>
    dplyr::mutate(
      amrfinder_gff_path = paste0(annotation_stem, ".gff"),
      faa_path = paste0(annotation_stem, ".faa"),
      ffn_path = paste0(annotation_stem, ".ffn"),
      fna_path = paste0(annotation_stem, ".fna")
    )
  annotation_files <- unlist(
    manifest[c("amrfinder_gff_path", "faa_path", "ffn_path", "fna_path")],
                             use.names = FALSE)
  if (!all(file.exists(annotation_files))) {
    stop("One or more selected Prokka FAA/FFN/FNA files are missing.", call. = FALSE)
  }
  manifest <- manifest |>
    dplyr::mutate(
      amrfinder_gff_sha256 = vapply(
        .data$amrfinder_gff_path, amr_sha256_file, character(1)
      ),
      faa_sha256 = vapply(.data$faa_path, amr_sha256_file, character(1)),
      ffn_sha256 = vapply(.data$ffn_path, amr_sha256_file, character(1)),
      fna_sha256 = vapply(.data$fna_path, amr_sha256_file, character(1)),
      selected_sequence_hash = vapply(.data$fasta_path, amr_sequence_content_hash,
                                      character(1)),
      prokka_sequence_hash = vapply(.data$fna_path, amr_sequence_content_hash,
                                    character(1)),
      annotation_sequence_equivalent =
        .data$selected_sequence_hash == .data$prokka_sequence_hash,
      annotation_bundle_sha256 = vapply(
        seq_len(dplyr::n()),
        function(i) digest::digest(
          c(.data$gff_sha256[[i]], .data$faa_sha256[[i]],
            .data$amrfinder_gff_sha256[[i]], .data$ffn_sha256[[i]],
            .data$fna_sha256[[i]]),
          algo = "sha256", serialize = TRUE
        ),
        character(1)
      )
    )
  if (!all(manifest$annotation_sequence_equivalent)) {
    stop("Selected FASTA and matched Prokka FNA are not sequence-equivalent for ",
         sum(!manifest$annotation_sequence_equivalent), " episode(s).", call. = FALSE)
  }
  manifest
}

amr_run_process <- function(command, args, stdout, stderr,
                            wd = NULL, env = character()) {
  dir.create(dirname(stdout), recursive = TRUE, showWarnings = FALSE)
  dir.create(dirname(stderr), recursive = TRUE, showWarnings = FALSE)
  started_at <- Sys.time()
  old_wd <- getwd()
  if (!is.null(wd)) {
    setwd(wd)
    on.exit(setwd(old_wd), add = TRUE)
  }
  error_message <- ""
  status <- tryCatch(
    suppressWarnings(system2(
      command, amr_quote_args(args), stdout = stdout, stderr = stderr,
      env = env, wait = TRUE
    )),
    error = function(e) {
      error_message <<- conditionMessage(e)
      127L
    }
  )
  if (is.null(status)) status <- 0L
  status <- as.integer(status)
  completed_at <- Sys.time()
  if (nzchar(error_message)) {
    cat(error_message, "\n", file = stderr, append = TRUE)
  }
  if (status != 0L) {
    stop("AMR caller failed (status ", status, "): ", command, " ",
         paste(args, collapse = " "), ". See ", stderr, call. = FALSE)
  }
  list(
    exit_status = status,
    started_at = format(started_at, "%Y-%m-%dT%H:%M:%S%z"),
    completed_at = format(completed_at, "%Y-%m-%dT%H:%M:%S%z"),
    command = paste(c(shQuote(command), amr_quote_args(args)), collapse = " ")
  )
}

amr_marker_values <- function(path) {
  if (!file.exists(path)) return(list())
  lines <- readLines(path, warn = FALSE)
  key_value <- lines[grepl("=", lines, fixed = TRUE)]
  if (!length(key_value)) return(list())
  keys <- sub("=.*$", "", key_value)
  values <- sub("^[^=]*=", "", key_value)
  stats::setNames(as.list(values), keys)
}

amr_cache_complete <- function(marker, cache_key, required_paths) {
  values <- amr_marker_values(marker)
  file.exists(marker) && length(required_paths) > 0L &&
    all(file.exists(required_paths)) &&
    identical(values$schema, AMR_SCHEMA_VERSION) &&
    identical(values$cache_key, cache_key) &&
    identical(values$exit_status, "0")
}

amr_write_completion_marker <- function(path, cache_key, run_result,
                                        caller, output_path) {
  amr_atomic_write_lines(
    c(
      paste0("schema=", AMR_SCHEMA_VERSION),
      paste0("cache_key=", cache_key),
      paste0("caller=", caller),
      paste0("exit_status=", run_result$exit_status),
      paste0("started_at=", run_result$started_at),
      paste0("completed_at=", run_result$completed_at),
      paste0("output_path=", output_path),
      paste0("command=", run_result$command)
    ),
    path
  )
}

amr_publish_file <- function(temporary, final) {
  if (!file.exists(temporary)) {
    stop("Caller did not produce expected file: ", temporary, call. = FALSE)
  }
  if (file.exists(final)) unlink(final)
  if (!file.rename(temporary, final)) {
    stop("Could not publish caller output: ", final, call. = FALSE)
  }
  invisible(final)
}

amr_run_callers <- function(manifest, runtime, provenance, output_root,
                            workers = 4L) {
  raw_root <- file.path(output_root, "raw")
  dirs <- file.path(raw_root, c("abricate_resfinder", "amrfinderplus", "resfinder"))
  vapply(dirs, dir.create, logical(1), recursive = TRUE, showWarnings = FALSE)

  version_for <- function(component) {
    provenance$version_or_commit[match(component, provenance$component)]
  }
  abr_params <- list(database = "resfinder", min_identity = 80, min_coverage = 80)
  amrfinder_threads <- suppressWarnings(as.integer(Sys.getenv(
    "AMRFINDER_THREADS_PER_CALL", "1"
  )))
  if (is.na(amrfinder_threads) || amrfinder_threads < 1L) {
    amrfinder_threads <- 1L
  }
  amrf_params <- list(
    organism = "Escherichia", mode = "nucleotide+protein+gff",
    annotation_format = "prokka",
    database_version = runtime$amrfinder_db_version,
    threads = amrfinder_threads
  )
  rf_params <- list(
    species = "Escherichia coli", acquired = TRUE, pointfinder = TRUE,
    min_identity = 0.80, min_coverage = 0.80
  )
  make_record <- function(row, caller, cache_key, output_path, marker,
                          stderr_path, cache_status) {
    values <- amr_marker_values(marker)
    marker_value <- function(name, default = NA_character_) {
      value <- values[[name]]
      if (is.null(value) || !length(value)) default else as.character(value)
    }
    tibble::tibble(
      Participant_id = row$Participant_id,
      tp_lab = row$tp_lab,
      Assembly_ID = row$Assembly_ID,
      caller = caller,
      cache_key = cache_key,
      output_path = output_path,
      stderr_path = stderr_path,
      completion_marker = marker,
      exit_status = suppressWarnings(as.integer(marker_value("exit_status"))),
      started_at = marker_value("started_at"),
      completed_at = marker_value("completed_at"),
      command = marker_value("command"),
      status = ifelse(identical(values$exit_status, "0"), "complete", "failed"),
      cache_status = cache_status
    )
  }
  worker_one <- function(i) {
    row <- manifest[i, , drop = FALSE]
    temporary_paths <- character()
    on.exit({
      for (path in temporary_paths) {
        if (file.exists(path)) unlink(path, recursive = TRUE, force = TRUE)
      }
    }, add = TRUE)
    identifier <- if (
      "Assembly_Base_ID" %in% names(row) &&
        !is.na(row$Assembly_Base_ID[[1L]]) &&
        nzchar(row$Assembly_Base_ID[[1L]])
    ) {
      row$Assembly_Base_ID[[1L]]
    } else {
      row$Assembly_ID[[1L]]
    }
    safe_id <- gsub("[^A-Za-z0-9_.-]", "_", identifier)
    records <- list()

    abr_key <- amr_cache_key(
      row$fasta_sha256, tool = "ABRicate-ResFinder",
      tool_version = version_for("abricate"),
      database_fingerprint = version_for("abricate_resfinder_database"),
      parameters = abr_params
    )
    abr_out <- file.path(dirs[[1L]], paste0(safe_id, "__", abr_key, ".tsv"))
    abr_err <- sub("\\.tsv$", ".stderr.txt", abr_out)
    abr_marker <- sub("\\.tsv$", ".complete.txt", abr_out)
    abr_reused <- amr_cache_complete(
      abr_marker, abr_key, c(abr_out, abr_err)
    )
    if (!abr_reused) {
      abr_tmp <- tempfile(
        paste0(".", safe_id, "__", abr_key, "."),
        tmpdir = dirs[[1L]], fileext = ".tsv.tmp"
      )
      abr_err_tmp <- paste0(abr_tmp, ".stderr.txt")
      temporary_paths <- c(temporary_paths, abr_tmp, abr_err_tmp)
      result <- amr_run_process(
        runtime$abricate,
        c("--quiet", "--threads", "1", "--db", "resfinder",
          "--minid", "80", "--mincov", "80", row$fasta_path),
        abr_tmp, abr_err_tmp
      )
      amr_publish_file(abr_tmp, abr_out)
      amr_publish_file(abr_err_tmp, abr_err)
      amr_write_completion_marker(
        abr_marker, abr_key, result, "abricate_resfinder", abr_out
      )
    }
    records[[1L]] <- make_record(
      row, "abricate_resfinder", abr_key, abr_out, abr_marker, abr_err,
      ifelse(abr_reused, "reused", "generated")
    )

    amrf_key <- amr_cache_key(
      row$fasta_sha256, row$annotation_bundle_sha256,
      tool = "AMRFinderPlus", tool_version = version_for("amrfinderplus"),
      database_fingerprint = version_for("amrfinderplus_database"),
      parameters = amrf_params
    )
    amrf_out <- file.path(dirs[[2L]], paste0(safe_id, "__", amrf_key, ".tsv"))
    amrf_err <- sub("\\.tsv$", ".stderr.txt", amrf_out)
    amrf_marker <- sub("\\.tsv$", ".complete.txt", amrf_out)
    amrf_reused <- amr_cache_complete(
      amrf_marker, amrf_key, c(amrf_out, amrf_err)
    )
    if (!amrf_reused) {
      amrf_tmp <- tempfile(
        paste0(".", safe_id, "__", amrf_key, "."),
        tmpdir = dirs[[2L]], fileext = ".tsv.tmp"
      )
      amrf_err_tmp <- paste0(amrf_tmp, ".stderr.txt")
      temporary_paths <- c(temporary_paths, amrf_tmp, amrf_err_tmp)
      result <- amr_run_process(
        runtime$amrfinder,
        c(
          "-n", row$fna_path,
          "-p", row$faa_path,
          "-g", row$amrfinder_gff_path,
          "-a", "prokka",
          "-O", "Escherichia",
          "-d", runtime$amrfinder_db,
          "--threads", as.character(amrfinder_threads),
          "--name", row$Assembly_ID
        ),
        amrf_tmp, amrf_err_tmp
      )
      amr_publish_file(amrf_tmp, amrf_out)
      amr_publish_file(amrf_err_tmp, amrf_err)
      amr_write_completion_marker(
        amrf_marker, amrf_key, result, "amrfinderplus", amrf_out
      )
    }
    records[[2L]] <- make_record(
      row, "amrfinderplus", amrf_key, amrf_out, amrf_marker, amrf_err,
      ifelse(amrf_reused, "reused", "generated")
    )

    rf_key <- amr_cache_key(
      row$fasta_sha256, tool = "ResFinder+PointFinder",
      tool_version = version_for("resfinder"),
      database_fingerprint = paste(
        version_for("resfinder_database"),
        version_for("pointfinder_database"), sep = "+"
      ),
      parameters = rf_params
    )
    rf_dir <- file.path(dirs[[3L]], paste0(safe_id, "__", rf_key))
    rf_marker <- file.path(rf_dir, "RUN_COMPLETE.txt")
    rf_err <- file.path(rf_dir, "stderr.txt")
    rf_stdout <- file.path(rf_dir, "stdout.txt")
    rf_reused <- amr_cache_complete(
      rf_marker, rf_key,
      c(
        rf_err, rf_stdout,
        file.path(rf_dir, "ResFinder_results_tab.txt"),
        file.path(rf_dir, "PointFinder_results.txt"),
        file.path(rf_dir, "pheno_table_escherichia_coli.txt")
      )
    )
    if (!rf_reused) {
      rf_tmp <- tempfile(
        paste0(".", safe_id, "__", rf_key, "."),
        tmpdir = dirs[[3L]]
      )
      dir.create(rf_tmp, recursive = TRUE, showWarnings = FALSE)
      temporary_paths <- c(temporary_paths, rf_tmp)
      result <- amr_run_process(
        runtime$python,
        c(
          "-m", "resfinder", "-ifa", row$fasta_path, "-o", rf_tmp,
          "-s", "Escherichia coli", "--acquired", "--point",
          "-db_res", runtime$resfinder_db, "-db_point", runtime$pointfinder_db,
          "-t", "0.80", "-l", "0.80"
        ),
        file.path(rf_tmp, "stdout.txt"), file.path(rf_tmp, "stderr.txt"),
        env = c(
          paste0("PATH=", dirname(runtime$kma), ":", Sys.getenv("PATH"))
        )
      )
      amr_write_completion_marker(
        file.path(rf_tmp, "RUN_COMPLETE.txt"), rf_key, result,
        "resfinder_pointfinder", rf_dir
      )
      if (dir.exists(rf_dir)) unlink(rf_dir, recursive = TRUE, force = TRUE)
      if (!file.rename(rf_tmp, rf_dir)) {
        stop("Could not publish ResFinder output directory: ", rf_dir)
      }
    }
    records[[3L]] <- make_record(
      row, "resfinder_pointfinder", rf_key, rf_dir, rf_marker, rf_err,
      ifelse(rf_reused, "reused", "generated")
    )
    message(sprintf(
      "[genomic AMR %d/%d] completed all callers for %s",
      i, nrow(manifest), identifier
    ))
    dplyr::bind_rows(records)
  }
  indices <- seq_len(nrow(manifest))
  workers <- max(1L, min(as.integer(workers), length(indices)))
  runs <- if (.Platform$OS.type == "unix" && workers > 1L) {
    parallel::mclapply(indices, worker_one, mc.cores = workers,
                       mc.preschedule = FALSE)
  } else {
    lapply(indices, worker_one)
  }
  failed <- which(vapply(runs, inherits, logical(1), what = "try-error"))
  if (length(failed)) {
    stop(
      "AMR caller worker(s) failed for manifest row(s): ",
      paste(failed, collapse = ", "),
      ". First error: ", as.character(runs[[failed[[1L]]]]),
      call. = FALSE
    )
  }
  dplyr::bind_rows(runs)
}

amr_empty_harmonized <- function() {
  tibble::tibble(
    Participant_id = character(), tp_lab = character(),
    Assembly_ID = character(), caller = character(),
    determinant_type = character(), raw_symbol = character(),
    normalized_symbol = character(), gene_family = character(),
    drug_class = character(), subclass = character(), identity = double(),
    coverage = double(), contig = character(), start = double(), end = double(),
    accession = character(), product = character(), scope = character(),
    method = character(), background_flag = logical(),
    evidence_role = character()
  )
}

amr_read_abricate <- function(path, meta) {
  if (!file.exists(path) || file.size(path) == 0L) return(amr_empty_harmonized())
  x <- suppressWarnings(readr::read_tsv(
    path, show_col_types = FALSE, col_types = readr::cols(.default = "c"),
    name_repair = "minimal"
  ))
  if (!nrow(x) || !"GENE" %in% names(x)) return(amr_empty_harmonized())
  symbol <- amr_normalise_symbol(x$GENE)
  tibble::tibble(
    Participant_id = meta$Participant_id, tp_lab = meta$tp_lab,
    Assembly_ID = meta$Assembly_ID, caller = "abricate_resfinder",
    determinant_type = "acquired_gene", raw_symbol = as.character(x$GENE),
    normalized_symbol = symbol, gene_family = amr_gene_family(symbol),
    drug_class = amr_infer_class(
      symbol,
      amr_pick_column(x, c("RESISTANCE", "PHENOTYPE"), NA_character_)
    ),
    subclass = amr_pick_column(x, c("PRODUCT"), NA_character_),
    identity = amr_safe_number(amr_pick_column(x, c("%IDENTITY", "IDENTITY"))),
    coverage = amr_safe_number(amr_pick_column(x, c("%COVERAGE", "COVERAGE"))),
    contig = amr_pick_column(x, c("SEQUENCE", "CONTIG")),
    start = amr_safe_number(amr_pick_column(x, c("START"))),
    end = amr_safe_number(amr_pick_column(x, c("END"))),
    accession = amr_pick_column(x, c("ACCESSION")),
    product = amr_pick_column(x, c("PRODUCT")),
    scope = NA_character_, method = "BLAST nucleotide",
    background_flag = amr_is_background(symbol),
    evidence_role = "legacy comparison"
  ) |>
    dplyr::filter(
      (is.na(.data$identity) | .data$identity >= 80),
      (is.na(.data$coverage) | .data$coverage >= 80)
    ) |>
    # The historical comparison table is determinant-level rather than
    # alignment-level. ABRicate can emit the same gene, identity and coverage
    # more than once when equivalent alignments occur on repeated sequence.
    # Preserve every alignment in the SHA-bound raw TSV, but deterministically
    # retain one representative here so the clean legacy 80/80 regression
    # anchor remains exactly reproducible.
    dplyr::arrange(.data$contig, .data$start, .data$end) |>
    dplyr::distinct(
      .data$Participant_id, .data$tp_lab, .data$Assembly_ID,
      .data$caller, .data$determinant_type, .data$raw_symbol,
      .data$normalized_symbol, .data$gene_family, .data$drug_class,
      .data$subclass, .data$identity, .data$coverage, .data$accession,
      .data$product, .data$scope, .data$method, .data$background_flag,
      .data$evidence_role,
      .keep_all = TRUE
    )
}

amr_read_amrfinder <- function(path, meta) {
  if (!file.exists(path) || file.size(path) == 0L) return(amr_empty_harmonized())
  x <- suppressWarnings(readr::read_tsv(
    path, show_col_types = FALSE, col_types = readr::cols(.default = "c"),
    name_repair = "minimal"
  ))
  if (!nrow(x)) return(amr_empty_harmonized())
  type <- amr_pick_column(x, c("Type", "Element type"))
  subtype <- amr_pick_column(x, c("Subtype", "Element subtype"))
  method <- amr_pick_column(x, c("Method"))
  is_mutation <- grepl("POINT|MUTATION", paste(subtype, method),
                       ignore.case = TRUE)
  element_symbol <- amr_pick_column(
    x, c("Element symbol", "Gene symbol", "Gene", "Symbol")
  )
  element_name <- amr_pick_column(
    x, c(
      "Element name", "Closest reference name", "Name of closest sequence",
      "Protein name", "Mutation"
    )
  )
  element_symbol[is.na(element_symbol)] <- ""
  element_name[is.na(element_name)] <- ""
  raw <- ifelse(nzchar(element_symbol), element_symbol, element_name)
  symbol <- ifelse(
    is_mutation,
    amr_normalise_mutation(raw),
    amr_normalise_symbol(raw)
  )
  gene <- ifelse(
    is_mutation,
    sub("[:_].*$", "", element_symbol),
    element_symbol
  )
  tibble::tibble(
    Participant_id = meta$Participant_id, tp_lab = meta$tp_lab,
    Assembly_ID = meta$Assembly_ID, caller = "amrfinderplus",
    determinant_type = ifelse(is_mutation, "point_mutation", "acquired_gene"),
    raw_symbol = raw, normalized_symbol = symbol,
    gene_family = ifelse(is_mutation, amr_normalise_symbol(gene),
                         amr_gene_family(symbol)),
    drug_class = amr_infer_class(
      ifelse(nzchar(element_symbol), element_symbol, symbol),
      amr_pick_column(x, c("Class"))
    ),
    subclass = amr_pick_column(x, c("Subclass")),
    identity = amr_safe_number(amr_pick_column(
      x, c("% Identity to reference", "% Identity to reference sequence",
           "% Identity")
    )),
    coverage = amr_safe_number(amr_pick_column(
      x, c("% Coverage of reference", "% Coverage of reference sequence",
           "% Coverage")
    )),
    contig = amr_pick_column(x, c("Contig id", "Contig", "Sequence name")),
    start = amr_safe_number(amr_pick_column(x, c("Start", "Start coordinate"))),
    end = amr_safe_number(amr_pick_column(x, c("Stop", "End", "Stop coordinate"))),
    accession = amr_pick_column(
      x, c("Closest reference accession", "Accession of closest sequence",
           "Accession")
    ),
    product = element_name, scope = amr_pick_column(x, c("Scope")),
    method = method, background_flag = amr_is_background(symbol),
    evidence_role = "primary caller"
  ) |>
    dplyr::filter(
      nzchar(.data$normalized_symbol),
      is.na(type) | !nzchar(type) | toupper(type) == "AMR"
    ) |>
    dplyr::distinct()
}

amr_find_result_file <- function(directory, patterns) {
  files <- list.files(directory, recursive = TRUE, full.names = TRUE,
                      include.dirs = FALSE)
  for (pattern in patterns) {
    hit <- files[grepl(pattern, basename(files), ignore.case = TRUE)]
    if (length(hit)) return(hit[[1L]])
  }
  ""
}

amr_read_resfinder_gene_table <- function(directory, meta) {
  path <- amr_find_result_file(
    directory,
    c("^ResFinder_results_tab\\.txt$", "resfinder.*results.*\\.tsv$",
      "resfinder.*results.*\\.txt$")
  )
  if (!nzchar(path) || file.size(path) == 0L) return(amr_empty_harmonized())
  x <- suppressWarnings(readr::read_tsv(
    path, show_col_types = FALSE, col_types = readr::cols(.default = "c"),
    name_repair = "minimal"
  ))
  if (!nrow(x)) return(amr_empty_harmonized())
  raw <- amr_pick_column(
    x, c("Resistance gene", "Gene", "GENE", "Template", "Gene name")
  )
  symbol <- amr_normalise_symbol(raw)
  identity <- amr_safe_number(amr_pick_column(
    x, c("Identity", "%Identity", "%IDENTITY")
  ))
  identity[identity <= 1] <- identity[identity <= 1] * 100
  coverage <- amr_safe_number(amr_pick_column(
    x, c("Coverage", "%Coverage", "%COVERAGE")
  ))
  coverage[coverage <= 1] <- coverage[coverage <= 1] * 100
  position <- amr_pick_column(x, c("Position in contig"))
  position_start <- amr_safe_number(position)
  position_end <- suppressWarnings(readr::parse_number(
    sub("^.*[.][.]|^.*-", "", position)
  ))
  tibble::tibble(
    Participant_id = meta$Participant_id, tp_lab = meta$tp_lab,
    Assembly_ID = meta$Assembly_ID, caller = "resfinder",
    determinant_type = "acquired_gene", raw_symbol = raw,
    normalized_symbol = symbol, gene_family = amr_gene_family(symbol),
    drug_class = amr_infer_class(
      symbol,
      amr_pick_column(x, c("Phenotype", "Resistance", "Drug class"))
    ),
    subclass = amr_pick_column(x, c("Phenotype", "Resistance")),
    identity = identity, coverage = coverage,
    contig = amr_pick_column(x, c("Contig", "Contig name")),
    start = dplyr::coalesce(
      amr_safe_number(amr_pick_column(x, c("Start"))), position_start
    ),
    end = dplyr::coalesce(
      amr_safe_number(amr_pick_column(x, c("End"))), position_end
    ),
    accession = amr_pick_column(x, c("Accession no.", "Accession", "Template")),
    product = amr_pick_column(x, c("Phenotype", "Resistance")),
    scope = NA_character_, method = "ResFinder acquired-gene model",
    background_flag = amr_is_background(symbol),
    evidence_role = "complementary evidence"
  ) |>
    dplyr::filter(
      nzchar(.data$normalized_symbol),
      is.na(.data$identity) | .data$identity >= 80,
      is.na(.data$coverage) | .data$coverage >= 80
    ) |>
    dplyr::distinct()
}

amr_read_pointfinder_table <- function(directory, meta) {
  path <- amr_find_result_file(
    directory,
    c("^PointFinder_results\\.txt$", "pointfinder.*results.*\\.tsv$",
      "pointfinder.*results.*\\.txt$")
  )
  if (!nzchar(path) || file.size(path) == 0L) return(amr_empty_harmonized())
  x <- suppressWarnings(readr::read_tsv(
    path, show_col_types = FALSE, col_types = readr::cols(.default = "c"),
    name_repair = "minimal"
  ))
  if (!nrow(x)) return(amr_empty_harmonized())
  mutation <- amr_pick_column(x, c("Mutation"))
  gene <- amr_pick_column(x, c("Gene", "GENE", "Template"))
  gene[is.na(gene)] <- ""
  mutation[is.na(mutation)] <- ""
  gene <- ifelse(
    nzchar(gene), gene,
    sub("^([A-Za-z0-9-]+).*$", "\\1", mutation)
  )
  amino_change <- amr_pick_column(x, c("Amino acid change"))
  nucleotide_change <- amr_pick_column(x, c("Nucleotide change"))
  amino_change[is.na(amino_change)] <- ""
  nucleotide_change[is.na(nucleotide_change)] <- ""
  change <- ifelse(
    nzchar(amino_change), amino_change,
    ifelse(nzchar(nucleotide_change), nucleotide_change, mutation)
  )
  raw <- ifelse(
    nzchar(gene) & nzchar(change), paste0(gene, ":", change),
    ifelse(nzchar(mutation), mutation, gene)
  )
  tibble::tibble(
    Participant_id = meta$Participant_id, tp_lab = meta$tp_lab,
    Assembly_ID = meta$Assembly_ID, caller = "pointfinder",
    determinant_type = "point_mutation", raw_symbol = raw,
    normalized_symbol = amr_normalise_symbol(raw),
    gene_family = amr_normalise_symbol(gene),
    drug_class = amr_infer_class(
      gene, amr_pick_column(x, c("Resistance", "Phenotype"))
    ),
    subclass = amr_pick_column(x, c("Resistance", "Phenotype")),
    identity = amr_safe_number(amr_pick_column(x, c("Identity", "%Identity"))),
    coverage = NA_real_, contig = amr_pick_column(x, c("Contig")),
    start = amr_safe_number(amr_pick_column(x, c("Position", "Start"))),
    end = NA_real_, accession = amr_pick_column(x, c("Template", "Accession")),
    product = amr_pick_column(x, c("Resistance", "Phenotype")),
    scope = NA_character_, method = "PointFinder mutation model",
    background_flag = FALSE, evidence_role = "complementary evidence"
  ) |>
    dplyr::filter(nzchar(.data$normalized_symbol)) |>
    dplyr::distinct()
}

amr_read_predicted_phenotype <- function(directory, meta) {
  path <- amr_find_result_file(
    directory,
    c(
      "^pheno_table_escherichia_coli\\.txt$",
      "^pheno_table\\.txt$",
      "phenotype.*\\.tsv$", "phenotype.*\\.txt$"
    )
  )
  if (!nzchar(path) || file.size(path) == 0L) return(tibble::tibble())
  lines <- readLines(path, warn = FALSE)
  header_index <- grep(
    "^#[[:space:]]*Antimicrobial[[:space:]]*\t", lines
  )
  if (!length(header_index)) return(tibble::tibble())
  header_index <- header_index[[1L]]
  payload <- c(
    sub("^#[[:space:]]*", "", lines[[header_index]]),
    lines[seq.int(header_index + 1L, length(lines))]
  )
  payload <- payload[nzchar(payload) & !grepl("^#", payload)]
  x <- suppressWarnings(tibble::as_tibble(utils::read.delim(
    text = paste(payload, collapse = "\n"),
    sep = "\t", quote = "", comment.char = "",
    check.names = FALSE, stringsAsFactors = FALSE
  )))
  if (!nrow(x)) return(tibble::tibble())
  tibble::tibble(
    Participant_id = meta$Participant_id,
    tp_lab = meta$tp_lab,
    Assembly_ID = meta$Assembly_ID,
    antimicrobial = amr_pick_column(x, c("Antimicrobial")),
    drug_class = amr_pick_column(x, c("Class")),
    wgs_predicted_phenotype = amr_pick_column(
      x, c("WGS-predicted phenotype", "Predicted phenotype")
    ),
    match_level = suppressWarnings(as.integer(amr_pick_column(x, c("Match")))),
    genetic_background = amr_pick_column(x, c("Genetic background")),
    interpretation_scope = "genomic prediction—not phenotypic AST",
    source_file = path
  )
}

amr_parse_callers <- function(run_manifest) {
  rows <- split(run_manifest, seq_len(nrow(run_manifest)))
  harmonized <- lapply(rows, function(meta) {
    switch(
      meta$caller[[1L]],
      abricate_resfinder = amr_read_abricate(meta$output_path[[1L]], meta),
      amrfinderplus = amr_read_amrfinder(meta$output_path[[1L]], meta),
      resfinder_pointfinder = dplyr::bind_rows(
        amr_read_resfinder_gene_table(meta$output_path[[1L]], meta),
        amr_read_pointfinder_table(meta$output_path[[1L]], meta)
      ),
      amr_empty_harmonized()
    )
  })
  long <- dplyr::bind_rows(harmonized) |>
    dplyr::mutate(
      drug_class = dplyr::if_else(
        is.na(.data$drug_class) | !nzchar(.data$drug_class),
        amr_infer_class(.data$normalized_symbol), .data$drug_class
      )
    ) |>
    dplyr::distinct()
  phenotype_manifest <- run_manifest |>
    dplyr::filter(.data$caller == "resfinder_pointfinder")
  phenotype_rows <- lapply(
    split(phenotype_manifest, seq_len(nrow(phenotype_manifest))),
    function(meta) {
      amr_read_predicted_phenotype(meta$output_path[[1L]], meta)
    }
  )
  list(long = long, phenotype = dplyr::bind_rows(phenotype_rows))
}

amr_build_episode_profiles <- function(manifest, long) {
  primary <- long |>
    dplyr::filter(.data$caller == "amrfinderplus")
  acquired <- primary |>
    dplyr::filter(.data$determinant_type == "acquired_gene")
  informative <- acquired |>
    dplyr::filter(!.data$background_flag)
  mutations <- primary |>
    dplyr::filter(.data$determinant_type == "point_mutation")
  any_mdfa <- long |>
    dplyr::filter(.data$background_flag) |>
    dplyr::distinct(.data$Participant_id, .data$tp_lab) |>
    dplyr::mutate(mdfA_detected = TRUE)

  summarize_sets <- function(x, value, name) {
    x |>
      dplyr::group_by(.data$Participant_id, .data$tp_lab) |>
      dplyr::summarise("{name}" := amr_collapse({{ value }}), .groups = "drop")
  }
  genes <- summarize_sets(informative, .data$normalized_symbol,
                          "informative_acquired_genes")
  genes_with_background <- summarize_sets(acquired, .data$normalized_symbol,
                                          "primary_acquired_genes_including_background")
  classes <- summarize_sets(informative, .data$drug_class,
                            "informative_acquired_classes")
  mutation_sets <- summarize_sets(mutations, .data$normalized_symbol,
                                  "primary_known_mutations")
  profiles <- manifest |>
    dplyr::select(
      dplyr::any_of(c("Participant_id", "tp_lab", "Episode_ID", "Event_type",
                      "UTI_Status", "Assembly_ID", "Assembly_Base_ID", "Isolate_ID")),
      .data$fasta_path, .data$fasta_sha256
    ) |>
    dplyr::left_join(genes, by = c("Participant_id", "tp_lab")) |>
    dplyr::left_join(genes_with_background, by = c("Participant_id", "tp_lab")) |>
    dplyr::left_join(classes, by = c("Participant_id", "tp_lab")) |>
    dplyr::left_join(mutation_sets, by = c("Participant_id", "tp_lab")) |>
    dplyr::left_join(any_mdfa, by = c("Participant_id", "tp_lab")) |>
    dplyr::mutate(
      dplyr::across(
        c("informative_acquired_genes",
          "primary_acquired_genes_including_background",
          "informative_acquired_classes", "primary_known_mutations"),
        ~ tidyr::replace_na(as.character(.x), "")
      ),
      mdfA_detected = tidyr::replace_na(.data$mdfA_detected, FALSE),
      acquired_genes_sensitivity_including_mdfA = dplyr::if_else(
        .data$mdfA_detected &
          !grepl("(^|;)mdf\\(A\\)($|;)",
                 .data$primary_acquired_genes_including_background),
        dplyr::if_else(
          nzchar(.data$primary_acquired_genes_including_background),
          paste0(.data$primary_acquired_genes_including_background, ";mdf(A)"),
          "mdf(A)"
        ),
        .data$primary_acquired_genes_including_background
      ),
      amr_gene_count_informative = lengths(lapply(
        .data$informative_acquired_genes, amr_split
      )),
      amr_gene_count_including_mdfA = lengths(lapply(
        .data$acquired_genes_sensitivity_including_mdfA, amr_split
      )),
      amr_class_count = lengths(lapply(
        .data$informative_acquired_classes, amr_split
      )),
      amr_mutation_count = lengths(lapply(
        .data$primary_known_mutations, amr_split
      )),
      any_informative_acquired_amr = .data$amr_gene_count_informative > 0L,
      any_acquired_amr_including_mdfA = .data$amr_gene_count_including_mdfA > 0L,
      primary_profile_caller = "AMRFinderPlus acquired genes and known resistance mutations",
      interpretation_scope = "genomic determinants—not phenotypic AST"
    )
  profiles
}

amr_build_concordance <- function(long) {
  acquired <- long |>
    dplyr::filter(.data$determinant_type == "acquired_gene") |>
    dplyr::mutate(
      caller_group = dplyr::case_when(
        .data$caller == "amrfinderplus" ~ "AMRFinderPlus",
        .data$caller == "resfinder" ~ "ResFinder",
        .data$caller == "abricate_resfinder" ~ "ABRicate-ResFinder",
        TRUE ~ .data$caller
      ),
      determinant_key = .data$gene_family
    ) |>
    dplyr::distinct(
      .data$Participant_id, .data$tp_lab, .data$determinant_key,
      .data$drug_class, .data$caller_group
    ) |>
    dplyr::mutate(present = TRUE) |>
    tidyr::pivot_wider(
      names_from = .data$caller_group, values_from = .data$present,
      values_fill = FALSE, names_prefix = "called_by_"
    )
  needed <- c(
    "called_by_AMRFinderPlus", "called_by_ResFinder",
    "called_by_ABRicate-ResFinder"
  )
  for (column in setdiff(needed, names(acquired))) acquired[[column]] <- FALSE
  acquired |>
    dplyr::mutate(
      caller_count = rowSums(dplyr::across(dplyr::all_of(needed))),
      agreement_class = dplyr::case_when(
        .data$called_by_AMRFinderPlus &
          .data$called_by_ResFinder &
          .data[["called_by_ABRicate-ResFinder"]] ~ "all three callers",
        .data$called_by_AMRFinderPlus & .data$called_by_ResFinder ~
          "primary + full ResFinder",
        .data$called_by_AMRFinderPlus &
          .data[["called_by_ABRicate-ResFinder"]] ~
          "primary + legacy ABRicate",
        .data$called_by_AMRFinderPlus ~ "primary only",
        .data$called_by_ResFinder &
          .data[["called_by_ABRicate-ResFinder"]] ~
          "ResFinder + legacy (not promoted)",
        .data$called_by_ResFinder ~ "full ResFinder only (not promoted)",
        TRUE ~ "legacy ABRicate only (not promoted)"
      ),
      primary_profile_inclusion = .data$called_by_AMRFinderPlus,
      discrepancy_rule =
        "Caller agreement is audited descriptively; no two-of-three voting rule is used."
    ) |>
    dplyr::arrange(
      .data$Participant_id, .data$tp_lab, .data$determinant_key
    )
}

amr_build_prevalence <- function(profiles, long) {
  episode_gene <- long |>
    dplyr::filter(
      .data$caller == "amrfinderplus",
      .data$determinant_type == "acquired_gene",
      !.data$background_flag
    ) |>
    dplyr::distinct(
      .data$Participant_id, .data$tp_lab, gene = .data$normalized_symbol
    )
  episode_class <- long |>
    dplyr::filter(
      .data$caller == "amrfinderplus",
      .data$determinant_type == "acquired_gene",
      !.data$background_flag
    ) |>
    dplyr::distinct(
      .data$Participant_id, .data$tp_lab, resistance_class = .data$drug_class
    )
  episode_mutation <- long |>
    dplyr::filter(
      .data$caller == "amrfinderplus",
      .data$determinant_type == "point_mutation"
    ) |>
    dplyr::distinct(
      .data$Participant_id, .data$tp_lab,
      mutation = .data$normalized_symbol
    )
  gene_episode <- episode_gene |>
    dplyr::count(.data$gene, name = "episodes_positive") |>
    dplyr::mutate(
      denominator_episodes = nrow(profiles),
      episode_prevalence = .data$episodes_positive / .data$denominator_episodes
    )
  gene_resident <- episode_gene |>
    dplyr::distinct(.data$Participant_id, .data$gene) |>
    dplyr::count(.data$gene, name = "residents_positive") |>
    dplyr::mutate(
      denominator_residents = dplyr::n_distinct(profiles$Participant_id),
      resident_prevalence = .data$residents_positive / .data$denominator_residents
    )
  class_episode <- episode_class |>
    dplyr::count(.data$resistance_class, name = "episodes_positive") |>
    dplyr::mutate(
      denominator_episodes = nrow(profiles),
      episode_prevalence = .data$episodes_positive / .data$denominator_episodes
    )
  class_resident <- episode_class |>
    dplyr::distinct(.data$Participant_id, .data$resistance_class) |>
    dplyr::count(.data$resistance_class, name = "residents_positive") |>
    dplyr::mutate(
      denominator_residents = dplyr::n_distinct(profiles$Participant_id),
      resident_prevalence = .data$residents_positive / .data$denominator_residents
    )
  mutation_episode <- episode_mutation |>
    dplyr::count(.data$mutation, name = "episodes_positive") |>
    dplyr::mutate(
      denominator_episodes = nrow(profiles),
      episode_prevalence = .data$episodes_positive / .data$denominator_episodes
    )
  mutation_resident <- episode_mutation |>
    dplyr::distinct(.data$Participant_id, .data$mutation) |>
    dplyr::count(.data$mutation, name = "residents_positive") |>
    dplyr::mutate(
      denominator_residents = dplyr::n_distinct(profiles$Participant_id),
      resident_prevalence = .data$residents_positive / .data$denominator_residents
    )
  list(
    gene = dplyr::full_join(gene_episode, gene_resident, by = "gene") |>
      dplyr::arrange(dplyr::desc(.data$episode_prevalence), .data$gene),
    class = dplyr::full_join(
      class_episode, class_resident, by = "resistance_class"
    ) |>
      dplyr::arrange(
        dplyr::desc(.data$episode_prevalence), .data$resistance_class
      ),
    mutation = dplyr::full_join(
      mutation_episode, mutation_resident, by = "mutation"
    ) |>
      dplyr::arrange(
        dplyr::desc(.data$episode_prevalence), .data$mutation
      )
  )
}

amr_set_change <- function(from, to) {
  from <- amr_split(from)
  to <- amr_split(to)
  list(
    gained = setdiff(to, from), lost = setdiff(from, to),
    jaccard = amr_jaccard(from, to), changed = !setequal(from, to)
  )
}

amr_build_transitions <- function(root, profiles) {
  path <- file.path(root, "results", "longitudinal", "longcycler_transitions.csv")
  if (!file.exists(path)) stop("Canonical adjacent transition table is missing.")
  transitions <- readr::read_csv(path, show_col_types = FALSE) |>
    dplyr::mutate(
      Participant_id = as.character(.data$Participant_id),
      tp_from = normalise_timepoint_preserve_events(.data$tp_from),
      tp_to = normalise_timepoint_preserve_events(.data$tp_to)
    )
  if (nrow(transitions) != AMR_EXPECTED_ADJACENT ||
      anyDuplicated(transitions$pair_key)) {
    stop("AMR longitudinal universe is not exactly 371 unique adjacent pairs.")
  }
  endpoint <- profiles |>
    dplyr::select(
      .data$Participant_id, .data$tp_lab,
      .data$informative_acquired_genes,
      .data$acquired_genes_sensitivity_including_mdfA,
      .data$informative_acquired_classes,
      .data$primary_known_mutations, .data$mdfA_detected
    )
  from <- endpoint |>
    dplyr::rename_with(~ paste0(.x, "_from"), -.data$Participant_id) |>
    dplyr::rename(tp_from = .data$tp_lab_from)
  to <- endpoint |>
    dplyr::rename_with(~ paste0(.x, "_to"), -.data$Participant_id) |>
    dplyr::rename(tp_to = .data$tp_lab_to)
  transitions <- transitions |>
    dplyr::left_join(
      from, by = c("Participant_id", "tp_from"), relationship = "many-to-one"
    ) |>
    dplyr::left_join(
      to, by = c("Participant_id", "tp_to"), relationship = "many-to-one"
    )
  required_profiles <- c(
    "informative_acquired_genes_from", "informative_acquired_genes_to",
    "informative_acquired_classes_from", "informative_acquired_classes_to",
    "primary_known_mutations_from", "primary_known_mutations_to"
  )
  if (anyNA(transitions[required_profiles])) {
    stop("At least one adjacent-pair endpoint lacks an AMR profile.", call. = FALSE)
  }
  changes <- lapply(seq_len(nrow(transitions)), function(i) {
    informative <- amr_set_change(
      transitions$informative_acquired_genes_from[[i]],
      transitions$informative_acquired_genes_to[[i]]
    )
    sensitivity <- amr_set_change(
      transitions$acquired_genes_sensitivity_including_mdfA_from[[i]],
      transitions$acquired_genes_sensitivity_including_mdfA_to[[i]]
    )
    classes <- amr_set_change(
      transitions$informative_acquired_classes_from[[i]],
      transitions$informative_acquired_classes_to[[i]]
    )
    mutations <- amr_set_change(
      transitions$primary_known_mutations_from[[i]],
      transitions$primary_known_mutations_to[[i]]
    )
    tibble::tibble(
      informative_genes_gained = amr_collapse(informative$gained),
      informative_genes_lost = amr_collapse(informative$lost),
      n_informative_genes_gained = length(informative$gained),
      n_informative_genes_lost = length(informative$lost),
      any_informative_acquired_gene_gain_or_loss = informative$changed,
      acquired_gene_jaccard = informative$jaccard,
      acquired_gene_jaccard_including_mdfA = sensitivity$jaccard,
      classes_gained = amr_collapse(classes$gained),
      classes_lost = amr_collapse(classes$lost),
      any_class_change = classes$changed,
      mutations_gained = amr_collapse(mutations$gained),
      mutations_lost = amr_collapse(mutations$lost),
      any_mutation_profile_change = mutations$changed,
      mdfA_threshold_change_qc_flag =
        transitions$mdfA_detected_from[[i]] != transitions$mdfA_detected_to[[i]]
    )
  })
  dplyr::bind_cols(transitions, dplyr::bind_rows(changes)) |>
    dplyr::mutate(
      direct_snp_context = dplyr::if_else(
        .data$TotalSNPs <= 25, "≤25 SNP", ">25 SNP"
      ),
      focused_strain_context = dplyr::case_when(
        .data$TotalSNPs <= 25 ~ "same-strain context",
        !is.na(.data$ST_from) & !is.na(.data$ST_to) &
          as.character(.data$ST_from) != as.character(.data$ST_to) ~
          "replacement context",
        TRUE ~ "uncertain context"
      ),
      interpretation_scope =
        "Genomic AMR context only; AMR profiles do not define strain identity."
    )
}

amr_cluster_resample <- function(df) {
  ids <- unique(as.character(df$Participant_id))
  sampled <- sample(ids, length(ids), replace = TRUE)
  dplyr::bind_rows(lapply(seq_along(sampled), function(i) {
    x <- df[as.character(df$Participant_id) == sampled[[i]], , drop = FALSE]
    x$.bootstrap_cluster <- i
    x
  }))
}

amr_binary_estimand <- function(df) {
  z <- df |>
    dplyr::transmute(
      outcome = as.integer(.data$any_informative_acquired_gene_gain_or_loss),
      close25 = as.integer(.data$TotalSNPs <= 25),
      days = as.numeric(.data$days_between_samples)
    ) |>
    dplyr::filter(
      !is.na(.data$outcome), !is.na(.data$close25), is.finite(.data$days)
    )
  if (nrow(z) < 20L || dplyr::n_distinct(z$outcome) != 2L ||
      dplyr::n_distinct(z$close25) != 2L) {
    return(c(odds_ratio = NA_real_, adjusted_risk_difference = NA_real_))
  }
  fit <- suppressWarnings(stats::glm(
    outcome ~ close25 + splines::ns(days, df = 3),
    family = stats::binomial(), data = z,
    control = stats::glm.control(maxit = 100)
  ))
  p1 <- z
  p0 <- z
  p1$close25 <- 1L
  p0$close25 <- 0L
  c(
    odds_ratio = exp(unname(stats::coef(fit)[["close25"]])),
    adjusted_risk_difference = mean(
      stats::predict(fit, p1, type = "response") -
        stats::predict(fit, p0, type = "response"),
      na.rm = TRUE
    )
  )
}

amr_continuous_estimand <- function(df) {
  z <- df |>
    dplyr::transmute(
      outcome = as.numeric(.data$acquired_gene_jaccard),
      close25 = as.integer(.data$TotalSNPs <= 25),
      days = as.numeric(.data$days_between_samples)
    ) |>
    dplyr::filter(
      is.finite(.data$outcome), !is.na(.data$close25), is.finite(.data$days)
    )
  if (nrow(z) < 20L || dplyr::n_distinct(z$close25) != 2L) return(NA_real_)
  fit <- stats::lm(outcome ~ close25 + splines::ns(days, df = 3), data = z)
  p1 <- z
  p0 <- z
  p1$close25 <- 1L
  p0$close25 <- 0L
  mean(stats::predict(fit, p1) - stats::predict(fit, p0), na.rm = TRUE)
}

amr_bootstrap_inference <- function(transitions, reps = 10000L,
                                    seed = AMR_BOOTSTRAP_SEED) {
  reps <- as.integer(reps)
  if (is.na(reps) || reps < 100L) {
    stop("AMR_BOOTSTRAP_REPS must be at least 100.", call. = FALSE)
  }
  point_binary <- amr_binary_estimand(transitions)
  point_jaccard <- amr_continuous_estimand(transitions)
  set.seed(seed)
  draws <- replicate(reps, {
    z <- amr_cluster_resample(transitions)
    c(amr_binary_estimand(z), jaccard_difference = amr_continuous_estimand(z))
  })
  draws <- t(draws)
  interval <- function(values, estimate, estimand) {
    good <- values[is.finite(values)]
    tibble::tibble(
      estimand = estimand, estimate = estimate,
      ci_lower = if (length(good)) unname(stats::quantile(good, 0.025)) else NA_real_,
      ci_upper = if (length(good)) unname(stats::quantile(good, 0.975)) else NA_real_,
      bootstrap_reps_requested = reps,
      bootstrap_reps_valid = length(good)
    )
  }
  dplyr::bind_rows(
    interval(draws[, "odds_ratio"], point_binary[["odds_ratio"]],
             "odds ratio: ≤25 versus >25 SNP"),
    interval(
      draws[, "adjusted_risk_difference"],
      point_binary[["adjusted_risk_difference"]],
      "adjusted risk difference: ≤25 minus >25 SNP"
    ),
    interval(
      draws[, "jaccard_difference"], point_jaccard,
      "adjusted mean Jaccard difference: ≤25 minus >25 SNP"
    )
  ) |>
    dplyr::mutate(
      n_pairs = nrow(transitions),
      n_residents = dplyr::n_distinct(transitions$Participant_id),
      seed = seed,
      model =
        "outcome ~ SNP≤25 + natural spline(days between samples, df=3)",
      multiplicity =
        "Prespecified profile-level outcomes; no gene-level hypothesis tests."
    )
}

amr_add_validation <- function(checks, check, observed, required, pass,
                               severity = "critical", detail = "") {
  checks[[length(checks) + 1L]] <- tibble::tibble(
    check = check, observed = as.character(observed),
    required = as.character(required), pass = isTRUE(pass),
    severity = severity, detail = detail
  )
  checks
}

amr_build_validation <- function(manifest, run_manifest, long, profiles,
                                 transitions, focused, phenotype, provenance) {
  checks <- list()
  add <- function(check, observed, required, pass,
                  severity = "critical", detail = "") {
    checks <<- amr_add_validation(
      checks, check, observed, required, pass, severity, detail
    )
  }
  keys <- paste(profiles$Participant_id, profiles$tp_lab, sep = "|")
  add("episode profiles", nrow(profiles), AMR_EXPECTED_EPISODES,
      nrow(profiles) == AMR_EXPECTED_EPISODES)
  add("unique episode keys", sum(duplicated(keys)), 0L, !anyDuplicated(keys))
  add("resident denominator", dplyr::n_distinct(profiles$Participant_id),
      AMR_EXPECTED_RESIDENTS,
      dplyr::n_distinct(profiles$Participant_id) == AMR_EXPECTED_RESIDENTS)
  add(
    "descriptive UTI-event profile denominator",
    sum(profiles$Event_type == "UTI_event", na.rm = TRUE), 32L,
    sum(profiles$Event_type == "UTI_event", na.rm = TRUE) == 32L,
    detail = "Descriptive only; not an additional UTI-vs-Not_UTI endpoint."
  )
  add("sequence-equivalent Prokka annotations",
      sum(manifest$annotation_sequence_equivalent), AMR_EXPECTED_EPISODES,
      all(manifest$annotation_sequence_equivalent))
  caller_counts <- run_manifest |>
    dplyr::filter(.data$status == "complete") |>
    dplyr::count(.data$caller)
  for (caller in c(
    "abricate_resfinder", "amrfinderplus", "resfinder_pointfinder"
  )) {
    observed <- caller_counts$n[match(caller, caller_counts$caller)]
    if (is.na(observed)) observed <- 0L
    add(paste(caller, "successful outputs"), observed, AMR_EXPECTED_EPISODES,
        observed == AMR_EXPECTED_EPISODES)
  }
  add(
    "cache completion markers",
    sum(vapply(
      seq_len(nrow(run_manifest)),
      function(i) amr_cache_complete(
        run_manifest$completion_marker[[i]],
        run_manifest$cache_key[[i]],
        c(run_manifest$output_path[[i]], run_manifest$stderr_path[[i]])
      ),
      logical(1)
    )),
    nrow(run_manifest),
    all(vapply(
      seq_len(nrow(run_manifest)),
      function(i) amr_cache_complete(
        run_manifest$completion_marker[[i]],
        run_manifest$cache_key[[i]],
        c(run_manifest$output_path[[i]], run_manifest$stderr_path[[i]])
      ),
      logical(1)
    ))
  )
  legacy <- long |>
    dplyr::filter(.data$caller == "abricate_resfinder")
  add("legacy ABRicate hit rows", nrow(legacy), 1401L, nrow(legacy) == 1401L,
      detail = "Regression anchor from the clean ≥80% identity/coverage screen.")
  legacy_positive <- legacy |>
    dplyr::distinct(.data$Participant_id, .data$tp_lab) |>
    nrow()
  add("legacy ABRicate hit-positive episodes", legacy_positive, 531L,
      legacy_positive == 531L)
  legacy_mdfa <- legacy |>
    dplyr::filter(.data$background_flag) |>
    dplyr::distinct(.data$Participant_id, .data$tp_lab) |>
    nrow()
  add("legacy mdf(A)-positive episodes", legacy_mdfa, 530L,
      legacy_mdfa == 530L)
  add("adjacent-pair AMR profiles", nrow(transitions), AMR_EXPECTED_ADJACENT,
      nrow(transitions) == AMR_EXPECTED_ADJACENT)
  add("focused Not_UTI→UTI AMR profiles", nrow(focused), AMR_EXPECTED_FOCUSED,
      nrow(focused) == AMR_EXPECTED_FOCUSED)
  add("allele suffix normalization", amr_normalise_symbol("blaTEM-1_1"),
      "blaTEM-1", identical(amr_normalise_symbol("blaTEM-1_1"), "blaTEM-1"))
  add("mdf(A) background flag", amr_is_background("mdfA"), TRUE,
      isTRUE(amr_is_background("mdfA")))
  add("mutation classification available",
      sum(long$determinant_type == "point_mutation"), ">=0 parsed calls",
      all(long$determinant_type %in% c("acquired_gene", "point_mutation")))
  phenotype_episodes <- phenotype |>
    dplyr::distinct(.data$Participant_id, .data$tp_lab) |>
    nrow()
  add(
    "ResFinder genomic-prediction profiles",
    phenotype_episodes, AMR_EXPECTED_EPISODES,
    phenotype_episodes == AMR_EXPECTED_EPISODES,
    detail = "Predicted phenotypes are genomic predictions, not phenotypic AST."
  )
  distorted <- amr_jaccard(c("mdf(A)", "blaTEM"), c("mdf(A)", "tet(A)"))
  informative <- amr_jaccard(c("blaTEM"), c("tet(A)"))
  add("mdf(A) Jaccard distortion unit test",
      paste(distorted, informative, sep = " / "), "1/3 including; 0 excluding",
      isTRUE(all.equal(distorted, 1 / 3)) &&
        isTRUE(all.equal(informative, 0)))
  version <- provenance$version_or_commit[
    match("amrfinderplus", provenance$component)
  ]
  add("AMRFinderPlus pinned version", version, AMR_AMRFINDER_VERSION,
      grepl(gsub("\\.", "\\\\.", AMR_AMRFINDER_VERSION), version))
  version <- provenance$version_or_commit[
    match("amrfinderplus_database", provenance$component)
  ]
  add(
    "AMRFinderPlus pinned project-local database", version,
    AMR_AMRFINDER_DB_VERSION,
    grepl(AMR_AMRFINDER_DB_VERSION, version, fixed = TRUE) &&
      startsWith(
        normalizePath(
          provenance$path[
            match("amrfinderplus_database", provenance$component)
          ],
          winslash = "/", mustWork = TRUE
        ),
        normalizePath(
          file.path(getwd(), "data", "amr_runtime", "databases"),
          winslash = "/", mustWork = TRUE
        )
      )
  )
  version <- provenance$version_or_commit[
    match("resfinder", provenance$component)
  ]
  add("ResFinder pinned version", version, AMR_RESFINDER_VERSION,
      grepl(gsub("\\.", "\\\\.", AMR_RESFINDER_VERSION), version))
  version <- provenance$version_or_commit[
    match("resfinder_database", provenance$component)
  ]
  add("ResFinder pinned database commit", version, AMR_RESFINDER_DB_COMMIT,
      identical(version, AMR_RESFINDER_DB_COMMIT))
  version <- provenance$version_or_commit[
    match("pointfinder_database", provenance$component)
  ]
  add("PointFinder pinned database commit", version, AMR_POINTFINDER_DB_COMMIT,
      identical(version, AMR_POINTFINDER_DB_COMMIT))
  add("deterministic cache output identities",
      sum(duplicated(run_manifest[c(
        "Participant_id", "tp_lab", "caller", "output_path"
      )])), 0L,
      !anyDuplicated(run_manifest[c(
        "Participant_id", "tp_lab", "caller", "output_path"
      )]),
      detail = "Each key binds input hashes, versions, database fingerprint and parameters.")
  dplyr::bind_rows(checks)
}

amr_publish_plots <- function(gene_prevalence, transitions, concordance,
                              plot_root) {
  dir.create(plot_root, recursive = TRUE, showWarnings = FALSE)
  transitions <- transitions |>
    dplyr::mutate(
      direct_snp_context = factor(
        .data$direct_snp_context, levels = c("≤25 SNP", ">25 SNP")
      )
    )
  top <- gene_prevalence |>
    dplyr::slice_max(.data$episode_prevalence, n = 20, with_ties = FALSE) |>
    dplyr::mutate(gene = stats::reorder(.data$gene, .data$episode_prevalence))
  p1 <- ggplot2::ggplot(
    top, ggplot2::aes(.data$episode_prevalence, .data$gene)
  ) +
    ggplot2::geom_col(fill = "#0072B2", width = 0.72) +
    ggplot2::geom_text(
      ggplot2::aes(label = paste0(.data$episodes_positive, "/532")),
      hjust = -0.08, size = 3
    ) +
    ggplot2::scale_x_continuous(labels = scales::percent_format(accuracy = 1)) +
    ggplot2::coord_cartesian(
      xlim = c(0, max(c(0.05, top$episode_prevalence * 1.2), na.rm = TRUE)),
      clip = "off"
    ) +
    ggplot2::labs(
      title = "Most prevalent informative acquired AMR genes",
      subtitle = "Primary AMRFinderPlus profile; mdf(A) excluded",
      x = "Episode prevalence (n=532)", y = NULL,
      caption =
        "Genomic determinant prevalence—not measured susceptibility. Descriptive only."
    ) +
    theme_ruti_publication()
  context_counts <- transitions |>
    dplyr::count(.data$direct_snp_context, name = "n_pairs") |>
    dplyr::mutate(
      context_label = paste0(.data$direct_snp_context, "\n(n=",
                             .data$n_pairs, " pairs)")
    )
  context_labels <- stats::setNames(
    context_counts$context_label, context_counts$direct_snp_context
  )
  p2 <- ggplot2::ggplot(
    transitions,
    ggplot2::aes(.data$direct_snp_context, .data$acquired_gene_jaccard,
                 fill = .data$direct_snp_context)
  ) +
    ggplot2::geom_boxplot(width = 0.55, outlier.shape = NA, alpha = 0.75) +
    ggplot2::geom_jitter(width = 0.14, alpha = 0.25, size = 0.9) +
    ggplot2::scale_fill_manual(
      values = c("≤25 SNP" = "#0072B2", ">25 SNP" = "#D55E00")
    ) +
    ggplot2::scale_x_discrete(labels = context_labels) +
    ggplot2::coord_cartesian(ylim = c(0, 1)) +
    ggplot2::labs(
      title = "AMR-profile stability by direct SNP context",
      subtitle = "371 adjacent pairs; informative acquired genes, excluding mdf(A)",
      x = "Direct pairwise SNP context", y = "Acquired-gene Jaccard similarity",
      caption =
        "Pair-level display. Resident-bootstrap estimates adjusted for time are tabulated."
    ) +
    theme_ruti_publication(legend_position = "none")
  concordance_summary <- concordance |>
    tidyr::pivot_longer(
      dplyr::starts_with("called_by_"), names_to = "caller",
      values_to = "detected"
    ) |>
    dplyr::filter(.data$detected) |>
    dplyr::mutate(
      caller = sub("^called_by_", "", .data$caller),
      drug_class = dplyr::if_else(
        is.na(.data$drug_class) | !nzchar(.data$drug_class),
        "Other/unspecified", .data$drug_class
      )
    ) |>
    dplyr::count(.data$drug_class, .data$caller, name = "episode_determinants")
  p3 <- ggplot2::ggplot(
    concordance_summary,
    ggplot2::aes(.data$caller, .data$drug_class, fill = .data$episode_determinants)
  ) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.25) +
    ggplot2::scale_fill_viridis_c(option = "C") +
    ggplot2::labs(
      title = "Caller concordance by determinant class",
      subtitle = paste0(
        "532 episodes per caller; caller-specific detections audited, ",
        "no two-of-three vote"
      ),
      x = NULL, y = NULL, fill = "Episode-\ndeterminants",
      caption =
        "Differences reflect caller databases and models; AMRFinderPlus defines the primary profile."
    ) +
    theme_ruti_publication()
  plots <- list(
    most_prevalent_informative_acquired_genes = p1,
    amr_profile_stability_by_direct_snp_context = p2,
    caller_concordance_by_determinant_class = p3
  )
  paths <- file.path(plot_root, paste0(names(plots), ".png"))
  for (i in seq_along(plots)) {
    ggplot2::ggsave(
      paths[[i]], plots[[i]], width = if (i == 3L) 9 else 7.5,
      height = if (i == 3L) 6.5 else 5.2, dpi = 300, bg = "white"
    )
  }
  paths
}

amr_write_interpretation <- function(path, profiles, transitions, focused,
                                     inference, validation) {
  close <- transitions |>
    dplyr::filter(.data$TotalSNPs <= 25)
  distant <- transitions |>
    dplyr::filter(.data$TotalSNPs > 25)
  lines <- c(
    "# Genomic AMR interpretation",
    "",
    sprintf(
      "Script 29 produced primary genomic-AMR profiles for %d episodes from %d residents.",
      nrow(profiles), dplyr::n_distinct(profiles$Participant_id)
    ),
    sprintf(
      "%d episodes carried at least one informative AMRFinderPlus acquired-gene determinant; mdf(A) was tracked separately in %d episodes.",
      sum(profiles$any_informative_acquired_amr),
      sum(profiles$mdfA_detected)
    ),
    sprintf(
      "Across %d adjacent pairs, informative acquired-gene gain/loss occurred in %d/%d ≤25-SNP pairs and %d/%d >25-SNP pairs.",
      nrow(transitions),
      sum(close$any_informative_acquired_gene_gain_or_loss), nrow(close),
      sum(distant$any_informative_acquired_gene_gain_or_loss), nrow(distant)
    ),
    sprintf(
      "The focused table contains all %d Not_UTI→UTI transitions; it is descriptive and no regression was fitted.",
      nrow(focused)
    ),
    "",
    "## How to read the callers",
    "",
    "- AMRFinderPlus acquired genes and known resistance mutations define the primary episode profile.",
    "- ResFinder/PointFinder provides complementary acquired-gene, mutation and predicted-phenotype evidence.",
    "- ABRicate–ResFinder is retained as a SHA-bound legacy comparison.",
    "- Caller agreement is an audit, not a two-of-three consensus vote.",
    "",
    "## mdf(A)",
    "",
    "mdf(A) is retained in raw results and sensitivity metrics, but excluded from primary burden, any-acquired-AMR, gain/loss and Jaccard calculations because its near-ubiquity makes it an uninformative background chromosomal efflux marker. An apparent longitudinal mdf(A) change is flagged for sequence/threshold QC.",
    "",
    "## Interpretation boundaries",
    "",
    "- These data support genomic AMR determinants and predicted phenotypes, not measured susceptibility.",
    "- No phenotypic AST or antibiotic-exposure data are available; clinical susceptibility and treatment-effect claims are prohibited.",
    "- AMR profiles add mechanism and longitudinal context but cannot define strain identity, transmission or UTI causation.",
    "- UTI-event profiles are descriptive only; no additional UTI-versus-Not_UTI endpoint was tested.",
    "",
    sprintf(
      "Validation: %d/%d checks passed. Adjusted inference rows: %d.",
      sum(validation$pass), nrow(validation), nrow(inference)
    )
  )
  writeLines(lines, path)
}

run_genomic_amr_analysis <- function(root = DIR_ROOT,
                                     output_root = file.path(root, "results", "amr"),
                                     plot_root = file.path(root, "plots", "amr")) {
  required_packages <- c(
    "Biostrings", "digest", "dplyr", "ggplot2", "purrr",
    "readr", "scales", "tibble", "tidyr"
  )
  absent <- required_packages[
    !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
  ]
  if (length(absent)) {
    stop("Missing R packages for genomic AMR: ", paste(absent, collapse = ", "))
  }
  dir.create(output_root, recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(output_root, "provenance"), recursive = TRUE,
             showWarnings = FALSE)

  manifest <- amr_validate_manifest(root)
  runtime <- amr_resolve_runtime(root)
  provenance <- amr_runtime_provenance(runtime)
  if (!grepl("4\\.2\\.7", provenance$version_or_commit[
    provenance$component == "amrfinderplus"
  ])) {
    stop("AMRFinderPlus must be pinned to 4.2.7.", call. = FALSE)
  }
  if (!grepl("4\\.7\\.2", provenance$version_or_commit[
    provenance$component == "resfinder"
  ])) {
    stop("ResFinder must be pinned to 4.7.2.", call. = FALSE)
  }
  if (!identical(
    provenance$version_or_commit[
      provenance$component == "resfinder_database"
    ],
    AMR_RESFINDER_DB_COMMIT
  )) {
    stop("ResFinder database checkout is not at the pinned commit.", call. = FALSE)
  }
  if (!identical(
    provenance$version_or_commit[
      provenance$component == "pointfinder_database"
    ],
    AMR_POINTFINDER_DB_COMMIT
  )) {
    stop("PointFinder database checkout is not at the pinned commit.", call. = FALSE)
  }
  if (!grepl(
    AMR_AMRFINDER_DB_VERSION,
    provenance$version_or_commit[
      provenance$component == "amrfinderplus_database"
    ],
    fixed = TRUE
  )) {
    stop("AMRFinderPlus database is not at the pinned version.", call. = FALSE)
  }
  workers <- suppressWarnings(as.integer(Sys.getenv(
    "AMR_CALLER_WORKERS", as.character(min(4L, parallel::detectCores()))
  )))
  if (is.na(workers) || workers < 1L) workers <- 1L
  run_manifest <- amr_run_callers(
    manifest, runtime, provenance, output_root, workers
  )
  parsed <- amr_parse_callers(run_manifest)
  long <- parsed$long
  profiles <- amr_build_episode_profiles(manifest, long)
  concordance <- amr_build_concordance(long)
  prevalence <- amr_build_prevalence(profiles, long)
  union_sets <- function(x) {
    amr_collapse(unique(unlist(lapply(x, amr_split), use.names = FALSE)))
  }
  resident_profiles <- profiles |>
    dplyr::group_by(.data$Participant_id) |>
    dplyr::summarise(
      n_episodes = dplyr::n(),
      informative_acquired_genes = union_sets(
        .data$informative_acquired_genes
      ),
      informative_acquired_classes = union_sets(
        .data$informative_acquired_classes
      ),
      primary_known_mutations = union_sets(.data$primary_known_mutations),
      mdfA_detected = any(.data$mdfA_detected),
      any_informative_acquired_amr = any(
        .data$any_informative_acquired_amr
      ),
      .groups = "drop"
    ) |>
    dplyr::mutate(
      amr_gene_count_informative = lengths(lapply(
        .data$informative_acquired_genes, amr_split
      )),
      amr_class_count = lengths(lapply(
        .data$informative_acquired_classes, amr_split
      )),
      amr_mutation_count = lengths(lapply(
        .data$primary_known_mutations, amr_split
      )),
      interpretation_scope = "genomic determinants—not phenotypic AST"
    )
  transitions <- amr_build_transitions(root, profiles)
  focused <- transitions |>
    dplyr::filter(.data$status_from == "Not_UTI", .data$status_to == "UTI")
  reps <- suppressWarnings(as.integer(Sys.getenv("AMR_BOOTSTRAP_REPS", "10000")))
  inference <- amr_bootstrap_inference(
    transitions, reps = reps, seed = AMR_BOOTSTRAP_SEED
  )
  validation <- amr_build_validation(
    manifest, run_manifest, long, profiles, transitions, focused,
    parsed$phenotype, provenance
  )

  amr_atomic_write_csv(
    manifest, file.path(output_root, "provenance", "input_manifest.csv")
  )
  amr_atomic_write_csv(
    provenance, file.path(output_root, "provenance", "tool_database_versions.csv")
  )
  amr_atomic_write_csv(
    run_manifest, file.path(output_root, "provenance", "run_manifest.csv")
  )
  amr_atomic_write_csv(
    long, file.path(output_root, "harmonized_determinants_long.csv")
  )
  amr_atomic_write_csv(
    profiles, file.path(output_root, "episode_amr_profiles.csv")
  )
  amr_atomic_write_csv(
    resident_profiles, file.path(output_root, "resident_amr_profiles.csv")
  )
  amr_atomic_write_csv(
    concordance, file.path(output_root, "caller_concordance_discrepancies.csv")
  )
  amr_atomic_write_csv(
    prevalence$gene, file.path(output_root, "gene_prevalence_episode_resident.csv")
  )
  amr_atomic_write_csv(
    prevalence$class, file.path(output_root, "class_prevalence_episode_resident.csv")
  )
  amr_atomic_write_csv(
    prevalence$mutation,
    file.path(output_root, "mutation_prevalence_episode_resident.csv")
  )
  amr_atomic_write_csv(
    parsed$phenotype,
    file.path(output_root, "resfinder_predicted_phenotypes_genomic_not_ast.csv")
  )
  amr_atomic_write_csv(
    transitions, file.path(output_root, "adjacent_pair_amr_profiles_371.csv")
  )
  amr_atomic_write_csv(
    focused, file.path(output_root, "not_uti_to_uti_amr_profiles_9.csv")
  )
  amr_atomic_write_csv(
    inference, file.path(output_root, "longitudinal_resident_bootstrap_inference.csv")
  )
  descriptive <- transitions |>
    dplyr::group_by(.data$direct_snp_context) |>
    dplyr::summarise(
      n_pairs = dplyr::n(),
      n_residents = dplyr::n_distinct(.data$Participant_id),
      any_gain_loss_n = sum(.data$any_informative_acquired_gene_gain_or_loss),
      any_gain_loss_pct =
        100 * mean(.data$any_informative_acquired_gene_gain_or_loss),
      median_jaccard = stats::median(.data$acquired_gene_jaccard),
      q25_jaccard = stats::quantile(.data$acquired_gene_jaccard, 0.25),
      q75_jaccard = stats::quantile(.data$acquired_gene_jaccard, 0.75),
      class_change_n = sum(.data$any_class_change),
      mutation_change_n = sum(.data$any_mutation_profile_change),
      .groups = "drop"
    )
  amr_atomic_write_csv(
    descriptive, file.path(output_root, "longitudinal_summary_by_snp_context.csv")
  )
  uti_events <- profiles |>
    dplyr::filter(.data$Event_type == "UTI_event") |>
    dplyr::mutate(
      analysis_role =
        "Descriptive UTI-event genomic AMR profile; no formal UTI-vs-Not_UTI test."
    )
  amr_atomic_write_csv(
    uti_events, file.path(output_root, "uti_event_amr_profiles_descriptive.csv")
  )
  amr_atomic_write_csv(
    validation, file.path(output_root, "validation_checks.csv")
  )
  plot_paths <- amr_publish_plots(
    prevalence$gene, transitions, concordance, plot_root
  )
  amr_write_interpretation(
    file.path(output_root, "interpretation_report.md"),
    profiles, transitions, focused, inference, validation
  )
  raw_summary <- long |>
    dplyr::count(.data$caller, .data$determinant_type, name = "n_calls")
  amr_atomic_write_csv(
    raw_summary, file.path(output_root, "caller_determinant_count_summary.csv")
  )
  caller_coverage <- run_manifest |>
    dplyr::group_by(.data$caller) |>
    dplyr::summarise(
      n_expected = AMR_EXPECTED_EPISODES,
      n_complete = sum(.data$status == "complete" & .data$exit_status == 0L),
      n_generated = sum(.data$cache_status == "generated"),
      n_reused = sum(.data$cache_status == "reused"),
      all_completion_markers_present = all(file.exists(.data$completion_marker)),
      .groups = "drop"
    ) |>
    dplyr::mutate(
      interpretation_scope = dplyr::case_when(
        .data$caller == "amrfinderplus" ~
          "Primary genomic determinant caller",
        .data$caller == "resfinder_pointfinder" ~
          "Complementary genomic prediction—not phenotypic AST",
        TRUE ~ "Legacy comparison"
      )
    )
  amr_atomic_write_csv(
    caller_coverage, file.path(output_root, "caller_coverage_summary.csv")
  )
  if (any(!validation$pass & validation$severity == "critical")) {
    failed <- validation$check[!validation$pass &
                                 validation$severity == "critical"]
    stop(
      "Genomic AMR validation failed: ", paste(failed, collapse = "; "),
      ". See results/amr/validation_checks.csv.", call. = FALSE
    )
  }
  published <- c(
    file.path(output_root, "provenance", c(
      "input_manifest.csv", "tool_database_versions.csv", "run_manifest.csv"
    )),
    file.path(output_root, c(
      "harmonized_determinants_long.csv",
      "episode_amr_profiles.csv", "resident_amr_profiles.csv",
      "caller_concordance_discrepancies.csv",
      "gene_prevalence_episode_resident.csv",
      "class_prevalence_episode_resident.csv",
      "mutation_prevalence_episode_resident.csv",
      "resfinder_predicted_phenotypes_genomic_not_ast.csv",
      "adjacent_pair_amr_profiles_371.csv",
      "not_uti_to_uti_amr_profiles_9.csv",
      "longitudinal_resident_bootstrap_inference.csv",
      "longitudinal_summary_by_snp_context.csv",
      "uti_event_amr_profiles_descriptive.csv",
      "validation_checks.csv", "interpretation_report.md",
      "caller_determinant_count_summary.csv", "caller_coverage_summary.csv"
    )),
    plot_paths
  )
  published <- unique(published[file.exists(published)])
  output_manifest <- tibble::tibble(
    path = normalizePath(published, winslash = "/", mustWork = TRUE),
    bytes = unname(file.info(published)$size),
    sha256 = vapply(published, amr_sha256_file, character(1)),
    generated_by = "29_vf_amr_combined_profile.R",
    schema = AMR_SCHEMA_VERSION
  )
  amr_atomic_write_csv(
    output_manifest,
    file.path(output_root, "provenance", "published_output_manifest.csv")
  )
  amr_atomic_write_lines(
    c(
      paste0("schema=", AMR_SCHEMA_VERSION),
      "status=complete",
      paste0("episodes=", nrow(profiles)),
      paste0("residents=", nrow(resident_profiles)),
      paste0("adjacent_pairs=", nrow(transitions)),
      paste0("focused_transitions=", nrow(focused)),
      paste0("completed_at=", format(
        Sys.time(), "%Y-%m-%dT%H:%M:%S%z"
      ))
    ),
    file.path(output_root, "RUN_COMPLETE.txt")
  )
  list(
    manifest = manifest, run_manifest = run_manifest, provenance = provenance,
    long = long, profiles = profiles, residents = resident_profiles,
    concordance = concordance,
    transitions = transitions, focused = focused, inference = inference,
    validation = validation, plots = plot_paths
  )
}

#!/usr/bin/env Rscript
# =============================================================
# 11_compare_strains_helpers.R
# Shared helpers for strain comparison across participants/timepoints
# Yellow RoUTIne – production-grade utilities
# =============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
  library(stringr)
  library(purrr)
  library(tibble)
})
if (!exists("FILE_VF_PA")) source("00_config.R")
if (!exists("prefer_primary_uti_status")) source("R/pipeline_qc_helpers.R")

# ------------- timepoint normalization ---------------------------------------
tp_norm <- function(x) {
  tp_chr <- as.character(x)
  is_uricult <- stringr::str_detect(tp_chr, stringr::regex("uricult", ignore_case = TRUE))
  tp_num <- suppressWarnings(as.integer(stringr::str_extract(tp_chr, "\\d+")))
  tp_lab <- dplyr::case_when(
    is_uricult ~ "Uricult",
    !is.na(tp_num) ~ paste0("T", tp_num),
    TRUE ~ "Unscheduled"
  )
  tp_levels <- c(
    paste0("T", sort(unique(tp_num[!is.na(tp_num)]))),
    "Uricult", "Unscheduled"
  )
  tibble::tibble(tp_lab = factor(tp_lab, levels = tp_levels), tp_num = tp_num)
}

# ------------- IO helpers ----------------------------------------------------
safe_dir_create <- function(path) {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  invisible(path)
}

safe_write_csv <- function(x, path, ...) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  tmp <- paste0(path, ".tmp")
  readr::write_csv(x, tmp, ...)
  if (file.exists(path)) unlink(path)
  file.rename(tmp, path)
  invisible(path)
}

archive_legacy_strain_outputs <- function(outdir) {
  pairwise_file <- file.path(outdir, "pairwise_metrics.csv")
  if (!file.exists(pairwise_file)) return(invisible(NULL))
  header <- tryCatch(
    names(readr::read_csv(pairwise_file, n_max = 0, show_col_types = FALSE)),
    error = function(e) character()
  )
  provenance_cols <- c(
    "Fasta_SHA256_A", "Fasta_SHA256_B", "dnadiff_cache_signature",
    "dnadiff_sidecar_path"
  )
  if (all(provenance_cols %in% header)) return(invisible(NULL))

  stamp <- format(Sys.time(), "%Y%m%dT%H%M%SZ", tz = "UTC")
  archive_dir <- file.path(outdir, "archive", paste0("pre_sha256_cache_", stamp))
  suffix <- 0L
  while (dir.exists(archive_dir)) {
    suffix <- suffix + 1L
    archive_dir <- file.path(outdir, "archive", paste0("pre_sha256_cache_", stamp, "_", suffix))
  }
  safe_dir_create(archive_dir)
  candidates <- c(
    "pairwise_metrics.csv", "summary_counts.csv", "summary_by_participant.csv",
    "stats_within_vs_between.csv", "stats_by_status.csv", "README.txt"
  )
  existing <- candidates[file.exists(file.path(outdir, candidates))]
  copied <- file.copy(
    file.path(outdir, existing), file.path(archive_dir, existing),
    overwrite = FALSE, copy.date = TRUE
  )
  if (length(existing) && !all(copied)) {
    stop("Could not archive all pre-provenance strain-comparison outputs.")
  }
  legacy_cache <- file.path(outdir, "nucmer_cache")
  manifest <- c(
    "Pre-SHA256 strain-comparison archive",
    paste0("Archived at UTC: ", format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")),
    paste0("Copied files: ", if (length(existing)) paste(existing, collapse = ", ") else "none"),
    paste0("Legacy key-only cache preserved in place: ", normalizePath(legacy_cache, winslash = "/", mustWork = FALSE)),
    "The legacy cache is intentionally not read by the SHA-256 provenance workflow."
  )
  writeLines(manifest, file.path(archive_dir, "ARCHIVE_MANIFEST.txt"))
  invisible(archive_dir)
}

timestamp_msg <- function(...) {
  cat(format(Sys.time(), "[%Y-%m-%d %H:%M:%S] "))
  message(...)
}

usable_fasta_path <- function(path) {
  path <- as.character(path)
  !is.na(path) & nzchar(path) & path != "NA" & file.exists(path)
}

normalise_existing_path <- function(path) {
  path <- as.character(path)
  if (length(path) != 1L || !usable_fasta_path(path)) {
    stop("Expected one existing FASTA path, received: ", paste(path, collapse = ", "))
  }
  normalizePath(path, winslash = "/", mustWork = TRUE)
}

fasta_fingerprint <- function(path) {
  if (!requireNamespace("digest", quietly = TRUE)) {
    stop("Package 'digest' is required for SHA-256 FASTA provenance.")
  }
  path_norm <- normalise_existing_path(path)
  info <- file.info(path_norm)
  list(
    path = path_norm,
    sha256 = unname(digest::digest(path_norm, algo = "sha256", file = TRUE)),
    size_bytes = unname(as.numeric(info$size)),
    mtime_utc = format(info$mtime, "%Y-%m-%dT%H:%M:%OS6Z", tz = "UTC")
  )
}

add_fasta_fingerprints <- function(assemblies) {
  if (!nrow(assemblies)) return(assemblies)
  fps <- lapply(assemblies$full_path, fasta_fingerprint)
  assemblies %>%
    mutate(
      full_path = vapply(fps, `[[`, character(1), "path"),
      full_path_norm = full_path,
      fasta_sha256 = vapply(fps, `[[`, character(1), "sha256"),
      fasta_size_bytes = vapply(fps, `[[`, numeric(1), "size_bytes"),
      fasta_mtime_utc = vapply(fps, `[[`, character(1), "mtime_utc")
    )
}

# ------------- tool detection ------------------------------------------------
has_tool <- function(bin) nzchar(Sys.which(bin))

detect_tools <- function() {
  list(
    nucmer = has_tool("nucmer"),
    dnadiff = has_tool("dnadiff"),
    `delta-filter` = has_tool("delta-filter"),
    mash = has_tool("mash")
  )
}

# ------------- load core tables ---------------------------------------------
load_core_tables <- function() {
  # Assemblies: the validated analysis manifest is the sole authority. It
  # contains selected, QC-passing Longcycler FASTAs only.
  asm <- load_analysis_assemblies(FILE_ANALYSIS_ASSEMBLY_MANIFEST, require_files = TRUE)
  required_asm <- c(
    "Participant_id", "tp_lab", "Isolate_ID", "assembler", "QC_PASS",
    "selected_canonical"
  )
  missing_asm <- setdiff(required_asm, names(asm))
  if (length(missing_asm)) {
    stop("Canonical assembly selection lacks columns: ", paste(missing_asm, collapse = ", "))
  }
  if (!"full_path" %in% names(asm) && !"fasta_path" %in% names(asm)) {
    stop("Canonical assembly selection must contain full_path or fasta_path.")
  }
  if (!"full_path" %in% names(asm)) asm$full_path <- asm$fasta_path
  asm <- asm %>%
    mutate(
      Participant_id = as.character(Participant_id),
      tp_lab = as.character(tp_lab),
      full_path = as.character(full_path),
      selected_canonical = as_pipeline_bool(selected_canonical),
      QC_PASS = as_pipeline_bool(QC_PASS)
    ) %>%
    filter(selected_canonical %in% TRUE, QC_PASS %in% TRUE)
  assert_analysis_assembly_manifest(
    asm,
    context = "strain-comparison assembly manifest",
    require_selected = TRUE,
    require_qc = TRUE,
    require_files = TRUE,
    require_unique_episode = TRUE
  )
  if (!nrow(asm)) stop("Canonical selection contains no selected QC-passing assemblies.")
  missing_fasta <- !usable_fasta_path(asm$full_path)
  if (any(missing_fasta)) {
    stop(
      sum(missing_fasta), " selected canonical FASTA(s) are missing or unusable:\n",
      paste(asm$full_path[missing_fasta], collapse = "\n")
    )
  }
  duplicate_episode <- asm %>% count(Participant_id, tp_lab, name = "n") %>% filter(n != 1L)
  if (nrow(duplicate_episode)) {
    stop(
      "Canonical selection must contain exactly one selected QC-passing FASTA per participant/timepoint. Problem keys:\n",
      paste0(duplicate_episode$Participant_id, "__", duplicate_episode$tp_lab, " (n=", duplicate_episode$n, ")", collapse = "\n")
    )
  }
  asm <- add_fasta_fingerprints(asm) %>% mutate(assembly_available = TRUE)

  # status map (optional; primary analyses use results/clinical/status_map.csv)
  status_map <- NULL
  status_candidates <- c(FILE_STATUS_MAP)
  status_file <- status_candidates[file.exists(status_candidates)][1]
  if (!is.na(status_file)) {
    status_map <- readr::read_csv(status_file, show_col_types = FALSE)
    if (!"tp_lab" %in% names(status_map)) {
      if (!"Timepoint" %in% names(status_map)) stop("status_map.csv must contain 'tp_lab' or 'Timepoint'")
      status_map <- dplyr::bind_cols(status_map, tp_norm(status_map$Timepoint))
    }
    status_map <- status_map %>%
      prefer_primary_uti_status() %>%
      select(
        Participant_id, tp_lab, Infection_Status, UTI_Status, UTI_binary,
        Not_UTI_subgroup, Infection_Status_legacy, Infection_Status_old,
        UTI_definition_version
      ) %>%
      distinct()
  }

  # MLST: active Longcycler-derived provider-preferred chromosomal STs only.
  mlst_file <- FILE_MLST_CANONICAL
  if (!file.exists(mlst_file)) {
    stop("No provider-preferred MLST file found - run 06_MLST.R, scripts/compare_mlst_sources.R, and scripts/integrate_provider_mlst.R")
  }
  mlst <- if (grepl("\\.tsv$", mlst_file, ignore.case = TRUE)) {
    readr::read_tsv(mlst_file, show_col_types = FALSE)
  } else {
    readr::read_csv(mlst_file, show_col_types = FALSE)
  }
  st_col <- names(mlst)[tolower(names(mlst)) == "st"]
  if (!length(st_col)) {
    alt <- names(mlst)[grepl("(?i)^sequence[_ ]?type$", names(mlst))]
    if (length(alt)) st_col <- alt[1]
  }
  if (!length(st_col)) stop("Could not find ST column in provider-preferred MLST file")
  if (!"ST" %in% names(mlst)) mlst$ST <- mlst[[st_col[1]]]
  if (!"full_path" %in% names(mlst)) stop(mlst_file, " lacks full_path; Longcycler provenance cannot be verified.")
  mlst <- mlst %>%
    mutate(full_path = normalizePath(full_path, winslash = "/", mustWork = FALSE)) %>%
    filter(full_path %in% asm$full_path)
  if ("ST_source" %in% names(mlst) && "provider_assembler" %in% names(mlst)) {
    bad_provider <- mlst$ST_source == "provider_qc95" &
      (is.na(mlst$provider_assembler) | tolower(mlst$provider_assembler) != ANALYSIS_ASSEMBLER)
    if (any(bad_provider, na.rm = TRUE)) {
      stop("MLST provider provenance includes a non-Longcycler call. Rerun 06_MLST.R.")
    }
  }
  # keep only relevant columns
  mlst <- mlst %>% select(Isolate_ID = any_of(c("Isolate_ID", "isolate_id")), ST, any_of(c("ST_source", "ST_provider", "ST_local")), everything())

  # pMLST (wide)
  pmlst_file <- "results/mlst/plasmid_types_per_isolate.csv"
  pmlst_wide <- NULL
  if (file.exists(pmlst_file)) {
    pmlst_wide <- readr::read_csv(pmlst_file, show_col_types = FALSE)
    # ensure Isolate_ID column case
    if (!"Isolate_ID" %in% names(pmlst_wide)) {
      iso <- names(pmlst_wide)[grepl("(?i)^isolate[_ ]?id$", names(pmlst_wide))][1]
      if (!is.na(iso)) names(pmlst_wide)[names(pmlst_wide) == iso] <- "Isolate_ID"
    }
  }

  # VF presence/absence
  if (!file.exists(FILE_VF_PA)) {
    stop(FILE_VF_PA, " not found – run 02_gene_presence_analysis.R")
  }
  vf_pa <- readr::read_csv(FILE_VF_PA, show_col_types = FALSE)
  if (!"tp_lab" %in% names(vf_pa)) {
    if (!"Timepoint" %in% names(vf_pa)) stop("vf_pa_all.csv must have 'tp_lab' or 'Timepoint'")
    vf_pa <- dplyr::bind_cols(vf_pa, tp_norm(vf_pa$Timepoint))
  }
  vf_pa <- vf_pa %>% mutate(SampleKey = paste(Participant_id, as.character(tp_lab), sep = "__"))

  # replicon presence/absence (plasmidfinder)
  inc_pa <- NULL
  inc_candidates <- c(
    file.path("results", "plasmids", "plasmidfinder_presence_absence.csv"),
    file.path("results", "plasmidfinder_presence_absence.csv")
  )
  inc_file <- inc_candidates[file.exists(inc_candidates)][1]
  if (!is.na(inc_file)) {
    inc_pa <- readr::read_csv(inc_file, show_col_types = FALSE)
    # ensure Isolate_ID column exists
    if (!"Isolate_ID" %in% names(inc_pa)) {
      iso <- names(inc_pa)[grepl("(?i)^isolate[_ ]?id$", names(inc_pa))][1]
      if (!is.na(iso)) names(inc_pa)[names(inc_pa) == iso] <- "Isolate_ID"
    }
  }

  # Assembly-first predicted plasmid profiles from numbered script 09b.
  # These are kept separate from PlasmidFinder replicon-marker profiles.
  mob_file <- file.path(
    "results", "plasmids", "mob_suite", "episode_plasmid_profiles.csv"
  )
  mob_marker <- file.path(
    "results", "plasmids", "mob_suite", "RUN_COMPLETE.txt"
  )
  mob_profiles <- NULL
  if (file.exists(mob_file) && file.exists(mob_marker)) {
    mob_profiles <- readr::read_csv(mob_file, show_col_types = FALSE)
    if (!"Isolate_ID" %in% names(mob_profiles)) {
      stop("MOB episode profiles lack Isolate_ID.")
    }
  }

  list(
    assemblies = asm,
    status_map = status_map,
    mlst = mlst,
    vf_pa = vf_pa,
    inc_pa = inc_pa,
    pmlst_wide = pmlst_wide,
    mob_profiles = mob_profiles
  )
}

# ------------- sample resolution --------------------------------------------
resolve_sample <- function(Participant_id, tp_lab, assemblies, prefer_assembler = NULL) {
  required <- c("Participant_id", "tp_lab", "full_path", "selected_canonical", "QC_PASS")
  missing <- setdiff(required, names(assemblies))
  if (length(missing)) stop("Assembly table lacks columns: ", paste(missing, collapse = ", "))
  df <- assemblies %>% filter(
    as.character(Participant_id) == as.character(!!Participant_id),
    as.character(tp_lab) == as.character(!!tp_lab),
    selected_canonical %in% TRUE,
    QC_PASS %in% TRUE,
    usable_fasta_path(full_path)
  )
  if (!nrow(df)) {
    return(NULL)
  }
  if (nrow(df) != 1L) {
    stop(
      "Expected exactly one selected QC-passing canonical assembly for ",
      Participant_id, "__", tp_lab, "; found ", nrow(df), "."
    )
  }
  assert_analysis_assembly_manifest(
    df,
    context = paste0("strain-comparison endpoint ", Participant_id, "__", tp_lab),
    require_selected = TRUE,
    require_qc = TRUE,
    require_files = TRUE,
    require_unique_episode = TRUE
  )
  df
}

# ------------- pair helpers --------------------------------------------------
make_sample_key <- function(pid, tp) paste(pid, as.character(tp), sep = "__")

# ------------- similarity metrics -------------------------------------------
parse_dnadiff_report <- function(report_path) {
  if (!file.exists(report_path)) {
    return(tibble(AvgIdentity = NA_real_, TotalSNPs = NA_real_))
  }
  L <- readLines(report_path, warn = FALSE)
  get_num <- function(key) {
    m <- grep(key, L, value = TRUE)
    if (!length(m)) {
      return(NA_real_)
    }
    as.numeric(stringr::str_extract(m[1], "\\d+\\.?\\d*"))
  }
  tibble(AvgIdentity = get_num("AvgIdentity"), TotalSNPs = get_num("TotalSNPs"))
}

dnadiff_cache_schema <- "dnadiff_sha256_cache_v1"

sanitize_cache_key <- function(key) {
  gsub("[^A-Za-z0-9_-]+", "_", as.character(key))
}

dnadiff_tool_version <- local({
  cached <- NULL
  function(bin = Sys.which("dnadiff")) {
    if (!is.null(cached)) return(cached)
    if (!nzchar(bin)) return(NA_character_)
    out <- suppressWarnings(tryCatch(
      system2(bin, "--version", stdout = TRUE, stderr = TRUE),
      error = function(e) character()
    ))
    cached <<- if (length(out)) paste(trimws(out[nzchar(trimws(out))]), collapse = " | ") else NA_character_
    cached
  }
})

dnadiff_cache_spec <- function(a_fasta, b_fasta, cache_dir, key, a_fingerprint = NULL, b_fingerprint = NULL) {
  if (!requireNamespace("digest", quietly = TRUE)) {
    stop("Package 'digest' is required for dnadiff cache signatures.")
  }
  a <- if (is.null(a_fingerprint)) fasta_fingerprint(a_fasta) else a_fingerprint
  b <- if (is.null(b_fingerprint)) fasta_fingerprint(b_fasta) else b_fingerprint
  # dnadiff is invoked in A/B order, so the signature deliberately preserves
  # input order.  The sidecar records both roles explicitly.
  signature_payload <- paste(
    dnadiff_cache_schema, as.character(key),
    a$path, a$sha256, b$path, b$sha256,
    sep = "\n"
  )
  signature <- digest::digest(signature_payload, algo = "sha256", serialize = FALSE)
  prefix <- file.path(
    cache_dir,
    paste0(sanitize_cache_key(key), "__", substr(signature, 1L, 24L))
  )
  list(
    schema_version = dnadiff_cache_schema,
    pair_key = as.character(key),
    cache_signature = signature,
    prefix = prefix,
    report_path = paste0(prefix, ".report"),
    sidecar_path = paste0(prefix, ".provenance.json"),
    log_path = paste0(prefix, ".dnadiff.log"),
    input_a = a,
    input_b = b
  )
}

read_dnadiff_sidecar <- function(path) {
  if (!file.exists(path) || !requireNamespace("jsonlite", quietly = TRUE)) return(NULL)
  tryCatch(jsonlite::read_json(path, simplifyVector = TRUE), error = function(e) NULL)
}

dnadiff_cache_is_valid <- function(spec) {
  if (!file.exists(spec$report_path) || !file.exists(spec$sidecar_path)) return(FALSE)
  sidecar <- read_dnadiff_sidecar(spec$sidecar_path)
  if (is.null(sidecar)) return(FALSE)
  required <- c("schema_version", "pair_key", "cache_signature", "input_a", "input_b", "report")
  if (!all(required %in% names(sidecar))) return(FALSE)
  report_hash <- tryCatch(fasta_fingerprint(spec$report_path)$sha256, error = function(e) NA_character_)
  identical(as.character(sidecar$schema_version), spec$schema_version) &&
    identical(as.character(sidecar$pair_key), spec$pair_key) &&
    identical(as.character(sidecar$cache_signature), spec$cache_signature) &&
    identical(as.character(sidecar$input_a$path), spec$input_a$path) &&
    identical(as.character(sidecar$input_a$sha256), spec$input_a$sha256) &&
    identical(as.character(sidecar$input_b$path), spec$input_b$path) &&
    identical(as.character(sidecar$input_b$sha256), spec$input_b$sha256) &&
    identical(as.character(sidecar$report$sha256), report_hash)
}

write_dnadiff_sidecar <- function(spec, command, version) {
  if (!requireNamespace("jsonlite", quietly = TRUE)) {
    stop("Package 'jsonlite' is required for dnadiff cache provenance.")
  }
  report_fp <- fasta_fingerprint(spec$report_path)
  payload <- list(
    schema_version = spec$schema_version,
    pair_key = spec$pair_key,
    cache_signature = spec$cache_signature,
    input_a = spec$input_a,
    input_b = spec$input_b,
    report = report_fp,
    dnadiff = list(
      executable = normalise_existing_path(Sys.which("dnadiff")),
      version = version,
      command = command,
      generated_at_utc = format(Sys.time(), "%Y-%m-%dT%H:%M:%OS6Z", tz = "UTC")
    )
  )
  tmp <- paste0(spec$sidecar_path, ".tmp-", Sys.getpid())
  jsonlite::write_json(payload, tmp, auto_unbox = TRUE, pretty = TRUE, na = "null")
  if (file.exists(spec$sidecar_path)) unlink(spec$sidecar_path)
  if (!file.rename(tmp, spec$sidecar_path)) stop("Could not atomically write ", spec$sidecar_path)
  invisible(payload)
}

empty_dnadiff_result <- function(cache_status = "unavailable") {
  tibble(
    AvgIdentity = NA_real_, TotalSNPs = NA_real_,
    dnadiff_report_path = NA_character_, dnadiff_sidecar_path = NA_character_,
    dnadiff_report_sha256 = NA_character_, dnadiff_cache_signature = NA_character_,
    dnadiff_cache_status = cache_status, dnadiff_version = NA_character_
  )
}

run_dnadiff <- function(a_fasta, b_fasta, cache_dir, key, a_fingerprint = NULL, b_fingerprint = NULL) {
  safe_dir_create(cache_dir)
  if (!has_tool("dnadiff") || !usable_fasta_path(a_fasta) || !usable_fasta_path(b_fasta)) {
    return(empty_dnadiff_result())
  }
  spec <- dnadiff_cache_spec(
    a_fasta, b_fasta, cache_dir, key,
    a_fingerprint = a_fingerprint, b_fingerprint = b_fingerprint
  )
  cache_valid <- dnadiff_cache_is_valid(spec)
  version <- dnadiff_tool_version()
  cache_status <- "reused"
  if (!cache_valid) {
    cache_status <- "generated"
    args <- c("-p", shQuote(spec$prefix), shQuote(spec$input_a$path), shQuote(spec$input_b$path))
    command <- paste(c("dnadiff", args), collapse = " ")
    status <- suppressWarnings(system2(
      "dnadiff", args,
      stdout = spec$log_path, stderr = spec$log_path
    ))
    if (status != 0) {
      return(empty_dnadiff_result("failed"))
    }
    if (!file.exists(spec$report_path)) return(empty_dnadiff_result("failed_missing_report"))
    write_dnadiff_sidecar(spec, command = command, version = version)
  }
  parsed <- parse_dnadiff_report(spec$report_path)
  report_sha256 <- fasta_fingerprint(spec$report_path)$sha256
  parsed %>%
    mutate(
      dnadiff_report_path = normalizePath(spec$report_path, winslash = "/", mustWork = TRUE),
      dnadiff_sidecar_path = normalizePath(spec$sidecar_path, winslash = "/", mustWork = TRUE),
      dnadiff_report_sha256 = report_sha256,
      dnadiff_cache_signature = spec$cache_signature,
      dnadiff_cache_status = cache_status,
      dnadiff_version = version
    )
}

mash_distance <- function(a_fasta, b_fasta) {
  if (!has_tool("mash") || !usable_fasta_path(a_fasta) || !usable_fasta_path(b_fasta)) {
    return(NA_real_)
  }
  out <- suppressWarnings(tryCatch(
    system2("mash", c("dist", a_fasta, b_fasta), stdout = TRUE, stderr = TRUE),
    error = function(e) character()
  ))
  if (!length(out)) {
    return(NA_real_)
  }
  # typical format: a\tb\tdist\tvar\tshared-hashes
  result_line <- out[vapply(strsplit(out, "\t", fixed = TRUE), function(parts) {
    length(parts) >= 3 && !is.na(suppressWarnings(as.numeric(parts[3])))
  }, logical(1))]
  if (!length(result_line)) {
    return(NA_real_)
  }
  parts <- strsplit(result_line[1], "\t", fixed = TRUE)[[1]]
  if (length(parts) < 3) {
    return(NA_real_)
  }
  as.numeric(parts[3])
}

booleanize_df <- function(df) {
  df[] <- lapply(df, function(col) {
    if (is.logical(col)) {
      return(as.integer(col))
    }
    if (is.numeric(col)) {
      return(as.integer(col > 0))
    }
    # character: treat "1"/"0" as ints, else 0/1
    suppressWarnings(as.integer(col)) %>% tidyr::replace_na(0L)
  })
  as.data.frame(df)
}

jaccard_from_wide <- function(df_wide, row_key_cols = c("Participant_id", "tp_lab"), sampleA, sampleB, cols_regex = NULL) {
  stopifnot(all(row_key_cols %in% names(df_wide)))
  dat <- df_wide
  if (!is.null(cols_regex)) {
    keep <- grep(cols_regex, names(dat), value = TRUE)
    dat <- dat[, c(row_key_cols, keep), drop = FALSE]
  }
  # ensure one row per sample key
  dat$SampleKey <- make_sample_key(dat[[row_key_cols[1]]], dat[[row_key_cols[2]]])
  A <- dat %>% filter(SampleKey == sampleA)
  B <- dat %>% filter(SampleKey == sampleB)
  if (!nrow(A) || !nrow(B)) {
    return(list(jaccard = NA_real_, n_int = NA_integer_, n_union = NA_integer_))
  }
  mat <- rbind(A, B) %>%
    select(-all_of(c(row_key_cols, "SampleKey"))) %>%
    booleanize_df() %>%
    as.matrix()
  if (!ncol(mat)) {
    return(list(jaccard = NA_real_, n_int = 0L, n_union = 0L))
  }
  a <- as.integer(mat[1, ] > 0)
  b <- as.integer(mat[2, ] > 0)
  inter <- sum(a & b)
  un <- sum(a | b)
  list(jaccard = if (un == 0) NA_real_ else inter / un, n_int = inter, n_union = un)
}

empty_profile_metrics <- function() {
  list(
    jaccard = NA_real_,
    n_intersection = NA_integer_,
    n_union = NA_integer_,
    both_empty = NA,
    available = FALSE,
    gains = NA_character_,
    losses = NA_character_,
    n_gains = NA_integer_,
    n_losses = NA_integer_
  )
}

set_profile_metrics <- function(a, b, available_a = TRUE, available_b = TRUE) {
  if (!isTRUE(available_a) || !isTRUE(available_b)) return(empty_profile_metrics())
  a <- sort(unique(as.character(a[!is.na(a) & nzchar(a)])))
  b <- sort(unique(as.character(b[!is.na(b) & nzchar(b)])))
  union_set <- union(a, b)
  intersection_set <- intersect(a, b)
  gains <- setdiff(b, a)
  losses <- setdiff(a, b)
  both_empty <- !length(a) && !length(b)
  list(
    jaccard = if (both_empty) 1 else length(intersection_set) / length(union_set),
    n_intersection = length(intersection_set),
    n_union = length(union_set),
    both_empty = both_empty,
    available = TRUE,
    gains = paste(gains, collapse = ";"),
    losses = paste(losses, collapse = ";"),
    n_gains = length(gains),
    n_losses = length(losses)
  )
}

binary_profile_metrics <- function(pa, isolate_a, isolate_b, feature_cols = NULL) {
  if (is.null(pa) || !"Isolate_ID" %in% names(pa)) return(empty_profile_metrics())
  row_a <- pa[match(isolate_a, pa$Isolate_ID), , drop = FALSE]
  row_b <- pa[match(isolate_b, pa$Isolate_ID), , drop = FALSE]
  if (!nrow(row_a) || !nrow(row_b) ||
      is.na(match(isolate_a, pa$Isolate_ID)) ||
      is.na(match(isolate_b, pa$Isolate_ID))) {
    return(empty_profile_metrics())
  }
  if (is.null(feature_cols)) feature_cols <- setdiff(names(pa), "Isolate_ID")
  feature_cols <- intersect(feature_cols, names(pa))
  if (!length(feature_cols)) return(set_profile_metrics(character(), character()))
  a_values <- suppressWarnings(as.numeric(unlist(
    row_a[1, feature_cols, drop = FALSE], use.names = FALSE
  )))
  b_values <- suppressWarnings(as.numeric(unlist(
    row_b[1, feature_cols, drop = FALSE], use.names = FALSE
  )))
  if (any(is.na(a_values)) || any(is.na(b_values))) return(empty_profile_metrics())
  set_profile_metrics(
    feature_cols[a_values > 0],
    feature_cols[b_values > 0]
  )
}

delimited_profile_metrics <- function(
    profiles, isolate_a, isolate_b, value_col,
    status_col = "mob_call_status", delimiter = ";") {
  if (is.null(profiles) ||
      !all(c("Isolate_ID", value_col, status_col) %in% names(profiles))) {
    return(empty_profile_metrics())
  }
  idx_a <- match(isolate_a, profiles$Isolate_ID)
  idx_b <- match(isolate_b, profiles$Isolate_ID)
  if (is.na(idx_a) || is.na(idx_b)) return(empty_profile_metrics())
  status_a <- identical(as.character(profiles[[status_col]][idx_a]), "complete")
  status_b <- identical(as.character(profiles[[status_col]][idx_b]), "complete")
  parse_values <- function(x) {
    if (is.na(x) || !nzchar(x)) return(character())
    trimws(strsplit(as.character(x), delimiter, fixed = TRUE)[[1L]])
  }
  set_profile_metrics(
    parse_values(profiles[[value_col]][idx_a]),
    parse_values(profiles[[value_col]][idx_b]),
    status_a, status_b
  )
}

# ------------- classification -------------------------------------------------
# [EPI] Strain Persistence Classification Rules
# This function defines what constitutes "strain persistence" (Same) versus
# "strain replacement" (Different) or "strain evolution/relatedness" (Related).
# These definitions are critical because they dictate which pairs are analyzed
# in the deep-dive phenotype-switch analysis (Script 16).
classify_pair <- function(metrics, thresholds = list(id = 99.9, snps = strain_snp_threshold(), vf = 0.9, inc = 0.8, vf_rel = 0.7, inc_rel = 0.7)) {
  ST_equal <- isTRUE(metrics$ST_equal)
  has_id <- !is.na(metrics$AvgIdentity) && !is.na(metrics$TotalSNPs)
  mean_acc <- mean(c(metrics$VF_Jaccard, metrics$Inc_Jaccard), na.rm = TRUE)

  # Helper: treat NA as passing (1.0) for Jaccard (implies empty union -> identical in emptiness)
  pass_j <- function(val, thresh) is.na(val) || val >= thresh

  # Same strain criteria
  same_rule <- ST_equal && (
    (!has_id || (metrics$AvgIdentity >= thresholds$id && metrics$TotalSNPs <= thresholds$snps))
  ) && pass_j(metrics$VF_Jaccard, thresholds$vf) && pass_j(metrics$Inc_Jaccard, thresholds$inc)

  if (same_rule) {
    return(list(Classification = "Same", RuleUsed = "ST + ID/SNP + accessory"))
  }

  # Related criteria
  related_rule <- ST_equal && (
    (has_id && metrics$AvgIdentity >= (thresholds$id - 0.9)) || (!is.na(metrics$MashDistance) && metrics$MashDistance <= 0.02) || (!is.na(mean_acc) && mean_acc >= 0.7)
  )
  if (related_rule) {
    return(list(Classification = "Related", RuleUsed = "ST + partial concordance"))
  }

  # Fallbacks when identity unavailable
  if (!has_id && ST_equal && !is.na(mean_acc)) {
    if (mean_acc >= 0.9) {
      return(list(Classification = "Same", RuleUsed = "ST + accessory >=0.9 (no ID)"))
    }
    if (mean_acc >= 0.7) {
      return(list(Classification = "Related", RuleUsed = "ST + accessory >=0.7 (no ID)"))
    }
  }
  list(Classification = "Different", RuleUsed = "Default")
}

# ------------- plotting helpers ---------------------------------------------
ggsave_safe <- function(filename, plot, width = 7, height = 5, dpi = 300) {
  dir.create(dirname(filename), recursive = TRUE, showWarnings = FALSE)
  width <- min(width, 49)
  height <- min(height, 49)
  if (grepl("\\.png$", filename, ignore.case = TRUE) && requireNamespace("ragg", quietly = TRUE)) {
    ggplot2::ggsave(filename, plot, device = ragg::agg_png, width = width, height = height, dpi = dpi, units = "in")
  } else {
    ggplot2::ggsave(filename, plot, width = width, height = height, dpi = dpi, units = "in")
  }
}

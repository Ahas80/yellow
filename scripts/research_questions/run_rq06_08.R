#!/usr/bin/env Rscript

# ==============================================================================
# RQ06-RQ08: VF stability, event-sample VF burden, and ST proxy performance
# ==============================================================================
# This module is isolated from superseded output generations.
# It rebuilds VFDB calls for the authoritative 532-row Longcycler manifest,
# binds every cache entry to the FASTA SHA-256 and VFDB/tool settings, and only
# uses pairwise rows whose endpoint paths and SHA-256 values match that manifest.
#
# Outputs are written exclusively below results/research_questions/.
# ==============================================================================

options(stringsAsFactors = FALSE, warn = 1)

required_packages <- c(
  "dplyr", "readr", "tidyr", "purrr", "stringr", "tibble",
  "ggplot2", "digest"
)
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages)) {
  stop("Missing required R package(s): ", paste(missing_packages, collapse = ", "))
}

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
  library(purrr)
  library(stringr)
  library(tibble)
  library(ggplot2)
})

EXPECTED_EPISODES <- 532L
EXPECTED_PARTICIPANTS <- 161L
EXPECTED_ALL_PAIRS <- 893L
EXPECTED_ADJACENT_PAIRS <- 371L
EXPECTED_EVENT_ROWS <- 32L
EXPECTED_EVENT_RESIDENTS <- 29L
EXPECTED_EVENT_UTI <- 15L
EXPECTED_EVENT_NOT_UTI <- 17L
EXPECTED_PAIRED_RESIDENTS <- 12L
VFDB_MIN_ID <- 80
VFDB_MIN_COV <- 80
PROVIDER_MLST_MIN_GOOD_TARGETS <- 95
BOOT_REPS <- suppressWarnings(as.integer(Sys.getenv("RQ_BOOTSTRAP_REPS", "10000")))
if (is.na(BOOT_REPS) || BOOT_REPS < 100L) {
  stop("RQ_BOOTSTRAP_REPS must be at least 100.", call. = FALSE)
}
RQ_SEED <- suppressWarnings(as.integer(Sys.getenv("RQ_SEED", "20260712")))
if (is.na(RQ_SEED)) RQ_SEED <- 20260712L
RQ_WORKERS <- suppressWarnings(as.integer(Sys.getenv("RQ_WORKERS", "8")))
if (is.na(RQ_WORKERS) || RQ_WORKERS < 1L) RQ_WORKERS <- 1L
detected_cores <- suppressWarnings(parallel::detectCores(logical = FALSE))
if (length(detected_cores) != 1L || is.na(detected_cores) || detected_cores < 2L) {
  detected_cores <- suppressWarnings(parallel::detectCores(logical = TRUE))
}
if (length(detected_cores) != 1L || is.na(detected_cores) || detected_cores < 2L) detected_cores <- 2L
RQ_WORKERS <- min(RQ_WORKERS, max(1L, detected_cores - 1L), 8L)

root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
if (!file.exists(file.path(root, "00_config.R"))) {
  stop("Run this script from the rUTIs project root; 00_config.R was not found.")
}

dir_rq <- file.path(root, "results", "research_questions")
dir_inputs <- file.path(dir_rq, "_inputs")
dir_cache <- file.path(dir_inputs, "vfdb_cache_sha256_v1")
dir_rq06 <- file.path(dir_rq, "RQ06")
dir_rq07 <- file.path(dir_rq, "RQ07")
dir_rq08 <- file.path(dir_rq, "RQ08")
for (d in c(dir_rq, dir_inputs, dir_cache, dir_rq06, dir_rq07, dir_rq08)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

path_manifest <- file.path(root, "results", "qc", "analysis_assembly_manifest.csv")
path_status <- file.path(root, "results", "clinical", "status_map.csv")
path_pairwise <- file.path(root, "results", "strain_compare", "pairwise_metrics.csv")
path_mlst <- file.path(root, "results", "mlst", "mlst_provider_preferred.csv")
path_module_map <- file.path(root, "results", "vf", "gene_module_map.csv")

messagef <- function(...) message(format(Sys.time(), "[%H:%M:%S] "), sprintf(...))

atomic_replace <- function(tmp, path) {
  if (file.exists(path) && !file.remove(path)) stop("Could not replace output: ", path)
  if (!file.rename(tmp, path)) stop("Could not atomically rename output: ", path)
  invisible(path)
}

atomic_write_csv <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  tmp <- tempfile(pattern = paste0(basename(path), "."), tmpdir = dirname(path), fileext = ".tmp")
  readr::write_csv(x, tmp, na = "")
  atomic_replace(tmp, path)
}

atomic_write_lines <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  tmp <- tempfile(pattern = paste0(basename(path), "."), tmpdir = dirname(path), fileext = ".tmp")
  writeLines(x, tmp, useBytes = TRUE)
  atomic_replace(tmp, path)
}

atomic_save_rds <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  tmp <- tempfile(pattern = paste0(basename(path), "."), tmpdir = dirname(path), fileext = ".tmp")
  saveRDS(x, tmp)
  atomic_replace(tmp, path)
}

atomic_ggsave <- function(plot, path, width, height, dpi = 300) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  ext <- tools::file_ext(path)
  tmp <- tempfile(pattern = paste0(tools::file_path_sans_ext(basename(path)), "."),
                  tmpdir = dirname(path), fileext = paste0(".tmp.", ext))
  ggplot2::ggsave(tmp, plot = plot, width = width, height = height, dpi = dpi, bg = "white")
  atomic_replace(tmp, path)
}

as_bool <- function(x) {
  if (is.logical(x)) return(replace(x, is.na(x), FALSE))
  tolower(trimws(as.character(x))) %in% c("true", "t", "1", "yes", "y")
}

norm_path <- function(x, must_work = FALSE) {
  normalizePath(as.character(x), winslash = "/", mustWork = must_work)
}

episode_key <- function(pid, tp) paste(as.character(pid), as.character(tp), sep = "||")

unordered_pair_id <- function(key_a, key_b) {
  lo <- ifelse(key_a <= key_b, key_a, key_b)
  hi <- ifelse(key_a <= key_b, key_b, key_a)
  paste(lo, hi, sep = "~~")
}

parse_date_multi <- function(x) {
  x <- as.character(x)
  out <- as.Date(rep(NA_character_, length(x)))
  for (fmt in c("%d/%m/%Y", "%d-%m-%Y", "%Y-%m-%d", "%Y/%m/%d")) {
    idx <- is.na(out) & !is.na(x) & nzchar(trimws(x))
    out[idx] <- as.Date(x[idx], format = fmt)
  }
  out
}

sha256_file <- function(path) {
  digest::digest(file = path, algo = "sha256", serialize = FALSE)
}

file_sha256_or_na <- function(path) {
  if (!file.exists(path)) return(NA_character_)
  sha256_file(path)
}

checks <- tibble(
  check_id = character(), status = character(), expected = character(),
  observed = character(), detail = character()
)

record_check <- function(check_id, status, expected, observed, detail = "") {
  checks <<- bind_rows(checks, tibble(
    check_id = check_id,
    status = status,
    expected = as.character(expected),
    observed = as.character(observed),
    detail = as.character(detail)
  ))
  invisible(identical(status, "PASS"))
}

write_run_state <- function(state, reason = "") {
  atomic_write_csv(checks, file.path(dir_inputs, "provenance_checks.csv"))
  atomic_write_csv(
    tibble(
      state = state,
      reason = reason,
      generated_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
      boot_reps = BOOT_REPS,
      rq_workers = RQ_WORKERS,
      script_sha256 = file_sha256_or_na(file.path(root, "scripts", "research_questions", "run_rq06_08.R"))
    ),
    file.path(dir_inputs, "execution_status.csv")
  )
}

block <- function(check_id, expected, observed, detail) {
  record_check(check_id, "FAIL", expected, observed, detail)
  write_run_state("BLOCKED", paste(check_id, detail, sep = ": "))
  stop("BLOCKED [", check_id, "]: ", detail, call. = FALSE)
}

require_columns <- function(df, columns, label) {
  missing <- setdiff(columns, names(df))
  if (length(missing)) {
    block(paste0(label, "_columns"), paste(columns, collapse = ";"),
          paste(names(df), collapse = ";"),
          paste("Missing required column(s):", paste(missing, collapse = ", ")))
  }
  record_check(paste0(label, "_columns"), "PASS", paste(columns, collapse = ";"),
               paste(columns, collapse = ";"), "Required columns present")
}

safe_ratio <- function(num, den) ifelse(den > 0, num / den, NA_real_)

auc_rank <- function(outcome, score) {
  keep <- is.finite(score) & !is.na(outcome)
  y <- as.integer(outcome[keep])
  s <- score[keep]
  n1 <- sum(y == 1L)
  n0 <- sum(y == 0L)
  if (!n1 || !n0) return(NA_real_)
  (sum(rank(s, ties.method = "average")[y == 1L]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

cluster_resample <- function(df, id_col = "Participant_id") {
  ids <- unique(as.character(df[[id_col]]))
  sampled <- sample(ids, length(ids), replace = TRUE)
  bind_rows(lapply(seq_along(sampled), function(i) {
    z <- df[as.character(df[[id_col]]) == sampled[i], , drop = FALSE]
    z$.bootstrap_cluster <- i
    z
  }))
}

bootstrap_stat <- function(df, statistic, reps = BOOT_REPS, seed = RQ_SEED) {
  set.seed(seed)
  out <- replicate(reps, {
    z <- cluster_resample(df)
    tryCatch(statistic(z), error = function(e) NA_real_)
  })
  as.numeric(out)
}

bootstrap_interval <- function(draws, estimate, conf = 0.95) {
  good <- draws[is.finite(draws)]
  alpha <- (1 - conf) / 2
  tibble(
    estimate = estimate,
    ci_lower = if (length(good)) unname(quantile(good, alpha, na.rm = TRUE)) else NA_real_,
    ci_upper = if (length(good)) unname(quantile(good, 1 - alpha, na.rm = TRUE)) else NA_real_,
    bootstrap_reps_requested = length(draws),
    bootstrap_reps_valid = length(good)
  )
}

messagef("Loading and validating authoritative Longcycler manifest")
for (p in c(path_manifest, path_status, path_pairwise, path_mlst, path_module_map)) {
  if (!file.exists(p)) block("required_input_exists", "all required inputs", p, "Required input is absent")
}

manifest <- read_csv(path_manifest, show_col_types = FALSE, progress = FALSE)
require_columns(
  manifest,
  c("Participant_id", "tp_lab", "Assembly_ID", "Isolate_ID", "full_path",
    "assembler", "selected_canonical", "QC_PASS", "Collection_Date", "Event_type"),
  "analysis_manifest"
)
manifest <- manifest %>%
  mutate(
    Participant_id = as.character(Participant_id),
    tp_lab = as.character(tp_lab),
    episode_key = episode_key(Participant_id, tp_lab),
    full_path = norm_path(full_path, must_work = FALSE),
    assembler = str_to_lower(as.character(assembler)),
    selected_canonical = as_bool(selected_canonical),
    QC_PASS = as_bool(QC_PASS)
  )

if (nrow(manifest) != EXPECTED_EPISODES) {
  block("manifest_n", EXPECTED_EPISODES, nrow(manifest), "Authoritative manifest denominator changed")
}
record_check("manifest_n", "PASS", EXPECTED_EPISODES, nrow(manifest), "Exact authoritative denominator")
if (n_distinct(manifest$Participant_id) != EXPECTED_PARTICIPANTS) {
  block("manifest_participants", EXPECTED_PARTICIPANTS, n_distinct(manifest$Participant_id),
        "Authoritative resident denominator changed")
}
record_check("manifest_participants", "PASS", EXPECTED_PARTICIPANTS,
             n_distinct(manifest$Participant_id), "Exact resident denominator")
if (anyDuplicated(manifest$episode_key)) {
  block("manifest_unique_keys", "532 unique keys", n_distinct(manifest$episode_key),
        "Duplicate participant-timepoint keys in authoritative manifest")
}
record_check("manifest_unique_keys", "PASS", EXPECTED_EPISODES,
             n_distinct(manifest$episode_key), "One row per participant-timepoint")
bad_manifest <- manifest$assembler != "longcycler" | !manifest$selected_canonical |
  !manifest$QC_PASS | !file.exists(manifest$full_path)
if (any(bad_manifest)) {
  block("manifest_contract", "all Longcycler, selected, QC-pass, existing",
        sum(bad_manifest), "Manifest violates the primary-analysis assembly contract")
}
record_check("manifest_contract", "PASS", "0 violations", 0,
             "All FASTAs are selected QC-passing Longcycler files")

messagef("Calculating SHA-256 for %d selected FASTAs using %d workers", nrow(manifest), RQ_WORKERS)
sha_values <- if (.Platform$OS.type == "unix" && RQ_WORKERS > 1L) {
  unlist(parallel::mclapply(manifest$full_path, sha256_file, mc.cores = RQ_WORKERS), use.names = FALSE)
} else {
  vapply(manifest$full_path, sha256_file, character(1))
}
manifest$fasta_sha256 <- sha_values
if (any(is.na(manifest$fasta_sha256) | nchar(manifest$fasta_sha256) != 64L)) {
  block("manifest_sha256", "532 valid SHA-256 values",
        sum(!is.na(manifest$fasta_sha256) & nchar(manifest$fasta_sha256) == 64L),
        "Could not fingerprint every selected FASTA")
}
record_check("manifest_sha256", "PASS", EXPECTED_EPISODES,
             sum(nchar(manifest$fasta_sha256) == 64L), "Every selected FASTA content-bound")

status <- read_csv(path_status, show_col_types = FALSE, progress = FALSE)
require_columns(
  status,
  c("Participant_id", "tp_lab", "UTI_Status", "analysis_include_primary",
    "Collection_Date", "Event_type", "Batch"),
  "status_map"
)
status <- status %>%
  mutate(
    Participant_id = as.character(Participant_id),
    tp_lab = as.character(tp_lab),
    episode_key = episode_key(Participant_id, tp_lab),
    analysis_include_primary = as_bool(analysis_include_primary)
  ) %>%
  filter(analysis_include_primary)
if (anyDuplicated(status$episode_key)) {
  block("status_unique_keys", "unique primary keys", n_distinct(status$episode_key),
        "Primary status map has duplicate participant-timepoint keys")
}
status_idx <- match(manifest$episode_key, status$episode_key)
if (anyNA(status_idx)) {
  block("manifest_status_key_match", "all 532 keys", sum(!is.na(status_idx)),
        "At least one authoritative assembly key lacks a primary clinical status")
}
if (any(as.character(manifest$Collection_Date) != as.character(status$Collection_Date[status_idx])) ||
    any(as.character(manifest$Event_type) != as.character(status$Event_type[status_idx]))) {
  block("manifest_status_metadata_match", "identical date and event type", "mismatch",
        "Clinical and assembly metadata disagree for an authoritative episode key")
}
record_check("manifest_status_key_match", "PASS", EXPECTED_EPISODES,
             sum(!is.na(status_idx)), "All authoritative keys have one primary status")
record_check("manifest_status_metadata_match", "PASS", "identical", "identical",
             "Collection date and event type agree for every key")

episode <- manifest %>%
  select(Participant_id, tp_lab, episode_key, Assembly_ID, Isolate_ID, full_path,
         fasta_sha256, assembler, Collection_Date, Event_type) %>%
  mutate(
    UTI_Status = as.character(status$UTI_Status[status_idx]),
    Batch = as.character(status$Batch[status_idx]),
    collection_date = parse_date_multi(Collection_Date)
  )
if (any(!episode$UTI_Status %in% c("UTI", "Not_UTI")) || anyNA(episode$collection_date)) {
  block("episode_analysis_fields", "binary status and parseable date for all 532", "invalid values",
        "A selected episode lacks a binary primary status or parseable collection date")
}

# ----------------------------------------------------------------------------
# Clean SHA-bound VFDB rebuild
# ----------------------------------------------------------------------------
abricate_bin <- Sys.which("abricate")
if (!nzchar(abricate_bin)) {
  block("abricate_available", "abricate in PATH", "missing",
        "A clean 532-isolate VFDB rebuild is required; stale 556 tables are forbidden")
}
abricate_version <- paste(system2(abricate_bin, "--version", stdout = TRUE, stderr = TRUE), collapse = " | ")
db_list <- system2(abricate_bin, "--list", stdout = TRUE, stderr = TRUE)
vfdb_line <- db_list[str_detect(db_list, "^vfdb\\s")]
if (!length(vfdb_line)) {
  block("vfdb_available", "vfdb listed by abricate", "missing",
        "ABRicate VFDB database is unavailable")
}
db_list_sha256 <- digest::digest(paste(db_list, collapse = "\n"), algo = "sha256", serialize = FALSE)
cache_context <- digest::digest(
  paste(abricate_version, vfdb_line, db_list_sha256, VFDB_MIN_ID, VFDB_MIN_COV, sep = "||"),
  algo = "sha256", serialize = FALSE
)
record_check("abricate_available", "PASS", "available", abricate_bin, abricate_version)
record_check("vfdb_available", "PASS", "listed", vfdb_line[1],
             paste0("database-list SHA-256=", db_list_sha256))

sanitize_file_part <- function(x) str_replace_all(as.character(x), "[^A-Za-z0-9._-]+", "_")

run_vfdb_one <- function(i) {
  row <- manifest[i, , drop = FALSE]
  stem <- paste0(
    sanitize_file_part(row$Assembly_ID), ".",
    substr(row$fasta_sha256, 1, 16), ".",
    substr(cache_context, 1, 12), ".id", VFDB_MIN_ID, ".cov", VFDB_MIN_COV
  )
  tab_path <- file.path(dir_cache, paste0(stem, ".tsv"))
  meta_path <- file.path(dir_cache, paste0(stem, ".meta.rds"))
  err_path <- file.path(dir_cache, paste0(stem, ".stderr.txt"))
  reused <- FALSE

  if (file.exists(tab_path) && file.exists(meta_path)) {
    meta <- tryCatch(readRDS(meta_path), error = function(e) NULL)
    valid <- !is.null(meta) && identical(meta$fasta_sha256, row$fasta_sha256) &&
      identical(meta$cache_context, cache_context) && identical(meta$full_path, row$full_path) &&
      identical(as.integer(meta$exit_status), 0L)
    if (valid) reused <- TRUE
  }

  if (!reused) {
    tmp_tab <- tempfile(pattern = paste0(stem, "."), tmpdir = dir_cache, fileext = ".tmp.tsv")
    tmp_err <- tempfile(pattern = paste0(stem, "."), tmpdir = dir_cache, fileext = ".tmp.err")
    exit_status <- suppressWarnings(system2(
      abricate_bin,
      c("--quiet", "--db", "vfdb", "--mincov", as.character(VFDB_MIN_COV),
        "--minid", as.character(VFDB_MIN_ID), row$full_path),
      stdout = tmp_tab, stderr = tmp_err
    ))
    exit_status <- as.integer(if (is.null(exit_status)) 0L else exit_status)
    if (file.exists(tmp_err) && file.size(tmp_err) > 0) {
      if (file.exists(err_path)) file.remove(err_path)
      file.rename(tmp_err, err_path)
    } else if (file.exists(tmp_err)) {
      file.remove(tmp_err)
    }
    if (file.exists(tab_path)) file.remove(tab_path)
    if (!file.rename(tmp_tab, tab_path)) exit_status <- 99L
    meta <- list(
      full_path = row$full_path,
      fasta_sha256 = row$fasta_sha256,
      cache_context = cache_context,
      exit_status = exit_status,
      abricate_version = abricate_version,
      vfdb_line = vfdb_line[1],
      generated_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")
    )
    atomic_save_rds(meta, meta_path)
  } else {
    exit_status <- 0L
  }

  n_hits <- if (exit_status == 0L && file.exists(tab_path)) {
    max(0L, length(readLines(tab_path, warn = FALSE)) - 1L)
  } else {
    NA_integer_
  }
  tibble(
    episode_key = row$episode_key,
    Participant_id = row$Participant_id,
    tp_lab = row$tp_lab,
    Assembly_ID = row$Assembly_ID,
    full_path = row$full_path,
    fasta_sha256 = row$fasta_sha256,
    cache_context = cache_context,
    cache_file = tab_path,
    cache_reused = reused,
    exit_status = exit_status,
    call_status = case_when(
      exit_status != 0L ~ "error",
      n_hits == 0L ~ "success_zero_hits",
      TRUE ~ "success_hits"
    ),
    n_hits = n_hits
  )
}

messagef("Running/reusing clean VFDB calls for %d FASTAs", nrow(manifest))
vfdb_run_list <- if (.Platform$OS.type == "unix" && RQ_WORKERS > 1L) {
  parallel::mclapply(seq_len(nrow(manifest)), run_vfdb_one, mc.cores = RQ_WORKERS)
} else {
  lapply(seq_len(nrow(manifest)), run_vfdb_one)
}
vfdb_run <- bind_rows(vfdb_run_list)

# ABRicate can fail transiently under a large parallel launch. Retry only failed
# calls once, sequentially, while preserving both attempt statuses for audit.
retry_keys <- vfdb_run %>%
  filter(.data$call_status == "error") %>%
  pull(.data$episode_key)
if (length(retry_keys)) {
  messagef("Retrying %d transient VFDB call(s) sequentially", length(retry_keys))
  retry_idx <- match(retry_keys, manifest$episode_key)
  retry_run <- bind_rows(lapply(retry_idx, run_vfdb_one))
  atomic_write_csv(
    bind_rows(
      vfdb_run %>%
        filter(.data$episode_key %in% retry_keys) %>%
        mutate(attempt = "parallel_initial"),
      retry_run %>% mutate(attempt = "sequential_retry")
    ),
    file.path(dir_inputs, "vfdb_retry_diagnostics.csv")
  )
  vfdb_run <- bind_rows(
    vfdb_run %>% filter(!.data$episode_key %in% retry_keys),
    retry_run
  ) %>% arrange(match(.data$episode_key, manifest$episode_key))
}
atomic_write_csv(vfdb_run, file.path(dir_inputs, "vfdb_run_manifest_532.csv"))
if (nrow(vfdb_run) != EXPECTED_EPISODES || any(vfdb_run$call_status == "error")) {
  block("vfdb_complete", "532 successful explicit hit/zero-hit results",
        paste0(nrow(vfdb_run), " rows; ", sum(vfdb_run$call_status == "error"), " errors"),
        "Clean VFDB screening did not complete for every selected FASTA")
}
record_check("vfdb_complete", "PASS", EXPECTED_EPISODES, nrow(vfdb_run),
             paste0(sum(vfdb_run$cache_reused), " SHA-valid caches reused; ",
                    sum(vfdb_run$call_status == "success_zero_hits"), " explicit zero-hit rows"))

read_vfdb_tab <- function(tab_path, episode_key_value, pid, tp, assembly_id, fasta_sha) {
  if (!file.exists(tab_path) || length(readLines(tab_path, n = 2L, warn = FALSE)) <= 1L) {
    return(tibble())
  }
  x <- suppressWarnings(read_tsv(tab_path, show_col_types = FALSE, progress = FALSE,
                                 col_types = cols(.default = col_character())))
  if (!nrow(x)) return(tibble())
  if (!"GENE" %in% names(x)) stop("VFDB output lacks GENE: ", tab_path)
  x %>%
    mutate(
      episode_key = episode_key_value,
      Participant_id = as.character(pid),
      tp_lab = as.character(tp),
      Assembly_ID = as.character(assembly_id),
      fasta_sha256 = fasta_sha,
      identity = suppressWarnings(parse_number(.data[["%IDENTITY"]])),
      coverage = suppressWarnings(parse_number(.data[["%COVERAGE"]]))
    )
}

vf_hits <- pmap_dfr(
  vfdb_run %>% transmute(
    tab_path = cache_file,
    episode_key_value = episode_key,
    pid = Participant_id,
    tp = tp_lab,
    assembly_id = Assembly_ID,
    fasta_sha = fasta_sha256
  ),
  read_vfdb_tab
)
if (nrow(vf_hits) && any(is.na(vf_hits$GENE) | !nzchar(trimws(vf_hits$GENE)))) {
  block("vfdb_gene_names", "all nonmissing", sum(!is.na(vf_hits$GENE)),
        "Clean VFDB output contains an empty gene name")
}
if (nrow(vf_hits) && "#FILE" %in% names(vf_hits)) {
  expected_path <- manifest$full_path[match(vf_hits$episode_key, manifest$episode_key)]
  hit_paths <- norm_path(vf_hits[["#FILE"]], must_work = FALSE)
  if (any(hit_paths != expected_path)) {
    block("vfdb_hit_path_identity", "all hit paths equal selected FASTA", "mismatch",
          "At least one VFDB hit was not generated from its selected FASTA")
  }
}
record_check("vfdb_hit_path_identity", "PASS", "all selected paths", "all selected paths",
             "Hit-level FASTA identity reconciled")

vf_hits_out <- vf_hits %>%
  select(episode_key, Participant_id, tp_lab, Assembly_ID, fasta_sha256,
         any_of(c("GENE", "ACCESSION", "PRODUCT", "SEQUENCE", "START", "END",
                  "STRAND", "%IDENTITY", "%COVERAGE", "identity", "coverage")))
atomic_write_csv(vf_hits_out, file.path(dir_inputs, "vfdb_hits_long_532.csv"))

gene_names <- sort(unique(as.character(vf_hits$GENE)))
vf_pa <- vf_hits %>%
  distinct(episode_key, GENE) %>%
  mutate(present = 1L) %>%
  pivot_wider(names_from = GENE, values_from = present, values_fill = 0L) %>%
  right_join(episode %>% select(episode_key), by = "episode_key")
for (g in gene_names) {
  if (!g %in% names(vf_pa)) vf_pa[[g]] <- 0L
  vf_pa[[g]][is.na(vf_pa[[g]])] <- 0L
}
vf_pa <- vf_pa %>% arrange(match(episode_key, episode$episode_key))
if (nrow(vf_pa) != EXPECTED_EPISODES || !setequal(vf_pa$episode_key, episode$episode_key)) {
  block("vf_matrix_keys", "exact 532 episode keys", nrow(vf_pa),
        "Fresh VF matrix does not equal the authoritative episode-key set")
}
record_check("vf_matrix_keys", "PASS", EXPECTED_EPISODES, nrow(vf_pa),
             paste(length(gene_names), "fresh VFDB gene columns"))

module_map <- read_csv(path_module_map, show_col_types = FALSE, progress = FALSE)
require_columns(module_map,
                c("Gene", "module_id", "assignment_confidence", "primary_assignment"),
                "gene_module_map")
module_map <- module_map %>%
  mutate(primary_assignment = as_bool(primary_assignment))
curated_genes <- module_map %>%
  filter(primary_assignment,
         assignment_confidence %in% c("High", "Moderate"),
         module_id != "unassigned") %>%
  pull(Gene) %>%
  unique() %>%
  intersect(gene_names)
if (!length(curated_genes)) {
  block("curated_gene_definition", "at least one fresh mapped gene", 0,
        "No fresh VFDB gene matches the frozen curated-gene definition")
}
record_check("curated_gene_definition", "PASS", "non-empty", length(curated_genes),
             paste0("gene_module_map SHA-256=", sha256_file(path_module_map)))

marker_present <- function(df, candidates) {
  candidates <- intersect(candidates, names(df))
  if (!length(candidates)) return(rep(0L, nrow(df)))
  as.integer(rowSums(as.data.frame(df[, candidates, drop = FALSE]), na.rm = TRUE) > 0)
}
afa_dra_genes <- unique(c(
  grep("^(afaB|afaC)", gene_names, value = TRUE),
  intersect(c("draB", "draC"), gene_names)
))
vf_metrics <- vf_pa %>%
  transmute(
    episode_key,
    total_vf_count_all = if (length(gene_names)) rowSums(as.data.frame(vf_pa[, gene_names, drop = FALSE])) else 0L,
    total_vf_count_curated = rowSums(as.data.frame(vf_pa[, curated_genes, drop = FALSE])),
    expec_pap = marker_present(vf_pa, c("papA", "papC")),
    expec_sfa_foc = marker_present(vf_pa, c("sfaD", "sfaE", "focD", "focE")),
    expec_afa_dra = marker_present(vf_pa, afa_dra_genes),
    expec_iutA = marker_present(vf_pa, "iutA"),
    expec_kpsM = marker_present(vf_pa, "kpsM")
  ) %>%
  mutate(five_marker_expec_count = expec_pap + expec_sfa_foc + expec_afa_dra + expec_iutA + expec_kpsM)

episode_vf <- episode %>%
  left_join(vf_metrics, by = "episode_key", relationship = "one-to-one")
if (anyNA(episode_vf$total_vf_count_curated) || nrow(episode_vf) != EXPECTED_EPISODES) {
  block("episode_vf_complete", "532 complete metric rows", sum(!is.na(episode_vf$total_vf_count_curated)),
        "Fresh episode VF metrics are incomplete")
}
record_check("episode_vf_complete", "PASS", EXPECTED_EPISODES, nrow(episode_vf),
             "Fresh curated burden and frozen five-marker count derived")

atomic_write_csv(
  vf_pa %>% left_join(episode %>% select(episode_key, Participant_id, tp_lab, Assembly_ID, fasta_sha256),
                      by = "episode_key", relationship = "one-to-one") %>%
    select(Participant_id, tp_lab, episode_key, Assembly_ID, fasta_sha256, all_of(gene_names)),
  file.path(dir_inputs, "vf_presence_absence_532.csv")
)
atomic_write_csv(episode_vf, file.path(dir_inputs, "episode_vf_metrics_532.csv"))

# ----------------------------------------------------------------------------
# Exact direct-pair universe and provider-only ST layer
# ----------------------------------------------------------------------------
make_expected_pairs <- function(df) {
  bind_rows(lapply(split(df, df$Participant_id), function(z) {
    if (nrow(z) < 2L) return(tibble())
    cmb <- combn(seq_len(nrow(z)), 2L)
    tibble(
      Participant_id = z$Participant_id[cmb[1, ]],
      episode_key_A_expected = z$episode_key[cmb[1, ]],
      episode_key_B_expected = z$episode_key[cmb[2, ]],
      expected_pair_id = unordered_pair_id(z$episode_key[cmb[1, ]], z$episode_key[cmb[2, ]])
    )
  }))
}

expected_pairs <- make_expected_pairs(episode)
if (nrow(expected_pairs) != EXPECTED_ALL_PAIRS || anyDuplicated(expected_pairs$expected_pair_id)) {
  block("expected_pair_universe", EXPECTED_ALL_PAIRS, nrow(expected_pairs),
        "The authoritative manifest no longer generates the approved direct-pair universe")
}
record_check("expected_pair_universe", "PASS", EXPECTED_ALL_PAIRS, nrow(expected_pairs),
             "All unordered within-resident combinations")

pairwise <- read_csv(path_pairwise, show_col_types = FALSE, progress = FALSE)
require_columns(
  pairwise,
  c("Participant_id_A", "Timepoint_A", "Participant_id_B", "Timepoint_B",
    "Fasta_path_A", "Fasta_path_B", "Fasta_SHA256_A", "Fasta_SHA256_B",
    "TotalSNPs", "MashDistance", "VF_Jaccard",
    "Replicon_Jaccard", "Replicon_Both_Empty",
    "Replicon_Profile_Available", "MOB_Cluster_Jaccard",
    "MOB_Cluster_Both_Empty", "MOB_Profile_Available",
    "Predicted_Plasmid_Count_A", "Predicted_Plasmid_Count_B",
    "MOB_High_Confidence_Profile_Both", "strict_same_strain",
    "pair_interpretation"),
  "pairwise_metrics"
)
pairwise <- pairwise %>%
  mutate(
    key_A = episode_key(Participant_id_A, Timepoint_A),
    key_B = episode_key(Participant_id_B, Timepoint_B),
    pair_id = unordered_pair_id(key_A, key_B),
    path_A = norm_path(Fasta_path_A, must_work = FALSE),
    path_B = norm_path(Fasta_path_B, must_work = FALSE)
  ) %>%
  filter(key_A %in% episode$episode_key, key_B %in% episode$episode_key)

if (nrow(pairwise) != EXPECTED_ALL_PAIRS || anyDuplicated(pairwise$pair_id) ||
    !setequal(pairwise$pair_id, expected_pairs$expected_pair_id)) {
  block("pairwise_key_universe", EXPECTED_ALL_PAIRS, nrow(pairwise),
        "Current pairwise rows do not equal the 893 expected unordered pairs")
}
idx_a <- match(pairwise$key_A, episode$episode_key)
idx_b <- match(pairwise$key_B, episode$episode_key)
endpoint_ok <- pairwise$path_A == episode$full_path[idx_a] &
  pairwise$path_B == episode$full_path[idx_b] &
  tolower(pairwise$Fasta_SHA256_A) == tolower(episode$fasta_sha256[idx_a]) &
  tolower(pairwise$Fasta_SHA256_B) == tolower(episode$fasta_sha256[idx_b])
if (any(!endpoint_ok) || anyNA(pairwise$TotalSNPs)) {
  block("pairwise_endpoint_identity", "893 path/SHA-matched pairs with SNPs",
        paste(sum(endpoint_ok), "identity matches;", sum(!is.na(pairwise$TotalSNPs)), "SNP rows"),
        "Stale or mismatched pairwise evidence is forbidden")
}
record_check("pairwise_key_universe", "PASS", EXPECTED_ALL_PAIRS, nrow(pairwise),
             "Exact expected pair IDs")
record_check("pairwise_endpoint_identity", "PASS", EXPECTED_ALL_PAIRS, sum(endpoint_ok),
             "Every endpoint path and SHA-256 equals the selected FASTA")

mlst <- read_csv(path_mlst, show_col_types = FALSE, progress = FALSE)
require_columns(mlst, c("full_path", "ST_source", "ST_provider"), "provider_mlst")
mlst <- mlst %>%
  mutate(
    full_path = norm_path(full_path, must_work = FALSE),
    ST_source = as.character(ST_source),
    ST_provider = as.character(ST_provider)
  ) %>%
  filter(full_path %in% episode$full_path)
if (anyDuplicated(mlst$full_path)) {
  block("provider_mlst_unique_path", "one MLST row per selected path", n_distinct(mlst$full_path),
        "Provider MLST contains duplicate selected FASTA paths")
}
provider_map <- episode %>%
  select(episode_key, full_path) %>%
  left_join(mlst %>% select(full_path, ST_source, ST_provider),
            by = "full_path", relationship = "one-to-one") %>%
  mutate(
    provider_ST = if_else(
      ST_source == "provider_qc95" & !is.na(ST_provider) & nzchar(ST_provider),
      ST_provider,
      NA_character_
    )
  )
record_check("provider_mlst_unique_path", "PASS", "unique selected paths", nrow(mlst),
             paste(sum(!is.na(provider_map$provider_ST)), "provider-only typed episodes"))

pair_dat <- pairwise %>%
  transmute(
    Participant_id = as.character(Participant_id_A),
    key_A, key_B, pair_id,
    TotalSNPs = as.numeric(TotalSNPs),
    MashDistance = as.numeric(MashDistance),
    pairwise_VF_Jaccard = as.numeric(VF_Jaccard),
    Replicon_Jaccard = as.numeric(Replicon_Jaccard),
    Replicon_Both_Empty = as_bool(Replicon_Both_Empty),
    Replicon_Profile_Available = as_bool(Replicon_Profile_Available),
    MOB_Cluster_Jaccard = as.numeric(MOB_Cluster_Jaccard),
    MOB_Cluster_Both_Empty = as_bool(MOB_Cluster_Both_Empty),
    MOB_Profile_Available = as_bool(MOB_Profile_Available),
    Predicted_Plasmid_Count_A = as.integer(Predicted_Plasmid_Count_A),
    Predicted_Plasmid_Count_B = as.integer(Predicted_Plasmid_Count_B),
    MOB_High_Confidence_Profile_Both =
      as_bool(MOB_High_Confidence_Profile_Both),
    strict_same_strain = as_bool(strict_same_strain),
    pair_interpretation = as.character(pair_interpretation)
  ) %>%
  mutate(
    provider_ST_A = provider_map$provider_ST[match(key_A, provider_map$episode_key)],
    provider_ST_B = provider_map$provider_ST[match(key_B, provider_map$episode_key)],
    provider_ST_both = !is.na(provider_ST_A) & !is.na(provider_ST_B),
    same_provider_ST = if_else(provider_ST_both, provider_ST_A == provider_ST_B, NA)
  )

atomic_write_csv(
  pair_dat %>%
    mutate(
      fasta_sha256_A = episode$fasta_sha256[match(key_A, episode$episode_key)],
      fasta_sha256_B = episode$fasta_sha256[match(key_B, episode$episode_key)]
    ),
  file.path(dir_inputs, "direct_pair_metrics_893.csv")
)

# The analysis sections are below. Their first operation independently rebuilds
# the 371 adjacent-pair set from dates and joins fresh VF calls by episode key.

# ----------------------------------------------------------------------------
# Shared pair-level VF calculations and the exact adjacent-pair universe
# ----------------------------------------------------------------------------
fresh_pair_metrics <- function(pair_table) {
  ia <- match(pair_table$key_A, vf_pa$episode_key)
  ib <- match(pair_table$key_B, vf_pa$episode_key)
  if (anyNA(ia) || anyNA(ib)) {
    block("fresh_pair_vf_keys", "all pair endpoints in fresh VF matrix", "missing",
          "A direct-pair endpoint is absent from the clean 532 VF matrix")
  }
  ma <- as.matrix(vf_pa[ia, gene_names, drop = FALSE])
  mb <- as.matrix(vf_pa[ib, gene_names, drop = FALSE])
  storage.mode(ma) <- "integer"
  storage.mode(mb) <- "integer"
  intersection_n <- rowSums(ma == 1L & mb == 1L)
  union_n <- rowSums(ma == 1L | mb == 1L)
  pair_table %>%
    mutate(
      fresh_vf_any_difference = rowSums(ma != mb) > 0L,
      fresh_vf_jaccard = if_else(union_n > 0L, intersection_n / union_n, NA_real_),
      fresh_vf_union_n = as.integer(union_n),
      fresh_vf_intersection_n = as.integer(intersection_n)
    )
}

pair_dat <- fresh_pair_metrics(pair_dat)

build_adjacent_pairs <- function(df) {
  bind_rows(lapply(split(df, df$Participant_id), function(z) {
    z <- z %>% arrange(collection_date, tp_lab, episode_key)
    if (nrow(z) < 2L) return(tibble())
    tibble(
      Participant_id = z$Participant_id[-nrow(z)],
      key_from = z$episode_key[-nrow(z)],
      key_to = z$episode_key[-1L],
      date_from = z$collection_date[-nrow(z)],
      date_to = z$collection_date[-1L],
      days_between = as.numeric(z$collection_date[-1L] - z$collection_date[-nrow(z)]),
      status_from = z$UTI_Status[-nrow(z)],
      status_to = z$UTI_Status[-1L],
      event_type_from = z$Event_type[-nrow(z)],
      event_type_to = z$Event_type[-1L],
      curated_from = z$total_vf_count_curated[-nrow(z)],
      curated_to = z$total_vf_count_curated[-1L],
      pair_id = unordered_pair_id(z$episode_key[-nrow(z)], z$episode_key[-1L])
    )
  }))
}

adjacent <- build_adjacent_pairs(episode_vf)
if (nrow(adjacent) != EXPECTED_ADJACENT_PAIRS || anyDuplicated(adjacent$pair_id) ||
    any(adjacent$days_between <= 0 | is.na(adjacent$days_between))) {
  block("adjacent_pair_universe", EXPECTED_ADJACENT_PAIRS, nrow(adjacent),
        "Date-ordered adjacent pairs are duplicated, non-positive, missing, or changed in number")
}
adjacent <- adjacent %>%
  left_join(pair_dat, by = c("Participant_id", "pair_id"), relationship = "one-to-one")
if (anyNA(adjacent$TotalSNPs) || nrow(adjacent) != EXPECTED_ADJACENT_PAIRS) {
  block("adjacent_direct_pair_join", "371 direct SNP rows", sum(!is.na(adjacent$TotalSNPs)),
        "Adjacent clinical pairs did not reconcile to direct path/SHA-validated pairs")
}
record_check("adjacent_pair_universe", "PASS", EXPECTED_ADJACENT_PAIRS, nrow(adjacent),
             paste(n_distinct(adjacent$Participant_id), "residents; all intervals positive"))
record_check("adjacent_direct_pair_join", "PASS", EXPECTED_ADJACENT_PAIRS,
             sum(!is.na(adjacent$TotalSNPs)), "Every adjacent pair has direct SNP evidence")

# Reorient the fresh VF matrix from arbitrary pairwise A/B order to chronological
# from/to order, then retain pair-level difference lists without testing genes.
idx_from <- match(adjacent$key_from, vf_pa$episode_key)
idx_to <- match(adjacent$key_to, vf_pa$episode_key)
mat_from <- as.matrix(vf_pa[idx_from, gene_names, drop = FALSE])
mat_to <- as.matrix(vf_pa[idx_to, gene_names, drop = FALSE])
storage.mode(mat_from) <- "integer"
storage.mode(mat_to) <- "integer"
genes_from_only <- lapply(seq_len(nrow(adjacent)), function(i) {
  gene_names[mat_from[i, ] == 1L & mat_to[i, ] == 0L]
})
genes_to_only <- lapply(seq_len(nrow(adjacent)), function(i) {
  gene_names[mat_from[i, ] == 0L & mat_to[i, ] == 1L]
})
adjacent <- adjacent %>%
  mutate(
    any_vf_difference = rowSums(mat_from != mat_to) > 0L,
    vf_jaccard = {
      un <- rowSums(mat_from == 1L | mat_to == 1L)
      int <- rowSums(mat_from == 1L & mat_to == 1L)
      ifelse(un > 0L, int / un, NA_real_)
    },
    abs_curated_count_change = abs(curated_to - curated_from),
    n_genes_from_only = lengths(genes_from_only),
    n_genes_to_only = lengths(genes_to_only),
    genes_from_only = vapply(genes_from_only, paste, collapse = ";", character(1)),
    genes_to_only = vapply(genes_to_only, paste, collapse = ";", character(1)),
    close25 = TotalSNPs <= 25,
    routine_routine = event_type_from == "Routine" & event_type_to == "Routine"
  )

atomic_write_csv(adjacent, file.path(dir_inputs, "adjacent_pair_metrics_371.csv"))

bootstrap_vector <- function(df, statistic, metric_names, reps = BOOT_REPS, seed = RQ_SEED) {
  set.seed(seed)
  values <- replicate(reps, {
    z <- cluster_resample(df)
    ans <- tryCatch(as.numeric(statistic(z)), error = function(e) rep(NA_real_, length(metric_names)))
    if (length(ans) != length(metric_names)) rep(NA_real_, length(metric_names)) else ans
  })
  values <- t(values)
  colnames(values) <- metric_names
  values
}

model_binary_standardized <- function(df, exposure) {
  z <- df %>%
    mutate(.exposure = as.integer(.data[[exposure]])) %>%
    filter(!is.na(any_vf_difference), !is.na(.exposure), is.finite(days_between))
  if (nrow(z) < 20L || n_distinct(z$.exposure) != 2L || n_distinct(z$any_vf_difference) != 2L) {
    return(c(log_or = NA_real_, adjusted_rd = NA_real_))
  }
  fit <- suppressWarnings(glm(
    any_vf_difference ~ .exposure + splines::ns(days_between, df = 3),
    family = binomial(), data = z, control = glm.control(maxit = 100)
  ))
  nd1 <- z
  nd0 <- z
  nd1$.exposure <- 1L
  nd0$.exposure <- 0L
  p1 <- suppressWarnings(predict(fit, newdata = nd1, type = "response"))
  p0 <- suppressWarnings(predict(fit, newdata = nd0, type = "response"))
  c(log_or = unname(coef(fit)[[".exposure"]]), adjusted_rd = mean(p1 - p0, na.rm = TRUE))
}

binary_model_inference <- function(df, exposure, label, seed_offset = 0L) {
  exposure_name <- exposure
  point <- model_binary_standardized(df, exposure_name)
  draws <- bootstrap_vector(
    df,
    function(z) model_binary_standardized(z, exposure_name),
    c("log_or", "adjusted_rd"),
    seed = RQ_SEED + seed_offset
  )
  bind_rows(
    bootstrap_interval(exp(draws[, "log_or"]), exp(point[["log_or"]])) %>%
      mutate(estimand = "odds_ratio", .before = 1),
    bootstrap_interval(draws[, "adjusted_rd"], point[["adjusted_rd"]]) %>%
      mutate(estimand = "adjusted_risk_difference", .before = 1)
  ) %>%
    mutate(
      analysis = label,
      exposure = exposure_name,
      n_pairs = nrow(df),
      n_residents = n_distinct(df$Participant_id),
      n_exposed = sum(df[[exposure_name]] %in% TRUE, na.rm = TRUE),
      outcome_events = sum(df$any_vf_difference %in% TRUE, na.rm = TRUE),
      model = "any binary VF difference ~ exposure + natural spline(days, df=3)",
      .before = 1
    )
}

model_continuous_standardized <- function(df, outcome, exposure = "close25") {
  z <- df %>%
    mutate(.exposure = as.integer(.data[[exposure]]), .outcome = as.numeric(.data[[outcome]])) %>%
    filter(is.finite(.outcome), !is.na(.exposure), is.finite(days_between))
  if (nrow(z) < 20L || n_distinct(z$.exposure) != 2L) return(NA_real_)
  fit <- lm(.outcome ~ .exposure + splines::ns(days_between, df = 3), data = z)
  nd1 <- z
  nd0 <- z
  nd1$.exposure <- 1L
  nd0$.exposure <- 0L
  mean(predict(fit, nd1) - predict(fit, nd0), na.rm = TRUE)
}

continuous_model_inference <- function(df, outcome, label, seed_offset = 0L) {
  point <- model_continuous_standardized(df, outcome)
  draws <- bootstrap_stat(
    df,
    function(z) model_continuous_standardized(z, outcome),
    seed = RQ_SEED + seed_offset
  )
  bootstrap_interval(draws, point) %>%
    mutate(
      analysis = label,
      outcome = outcome,
      estimand = "adjusted_mean_difference_close_minus_not_close",
      n_pairs = nrow(df),
      n_residents = n_distinct(df$Participant_id),
      model = paste0(outcome, " ~ SNP<=25 + natural spline(days, df=3)"),
      .before = 1
    )
}

# ==============================================================================
# RQ06: Does direct SNP closeness predict VF stability in adjacent isolates?
# ==============================================================================
messagef("RQ06: fitting resident-bootstrap adjacent-pair models")

rq06_primary <- binary_model_inference(adjacent, "close25", "primary_SNP_le_25", 600L)

rq06_threshold <- map_dfr(c(10L, 25L, 50L), function(threshold) {
  z <- adjacent %>% mutate(.threshold_exposure = TotalSNPs <= threshold)
  binary_model_inference(
    z, ".threshold_exposure", paste0("threshold_SNP_le_", threshold), 610L + threshold
  )
})

rq06_provider <- adjacent %>% filter(provider_ST_both) %>%
  mutate(provider_same_ST_exposure = same_provider_ST)
rq06_provider_inf <- binary_model_inference(
  rq06_provider, "provider_same_ST_exposure", "provider_only_same_ST_proxy", 700L
)

rq06_routine <- adjacent %>% filter(routine_routine)
rq06_routine_inf <- binary_model_inference(
  rq06_routine, "close25", "routine_to_routine_only_SNP_le_25", 710L
)

rq06_inference <- bind_rows(rq06_primary, rq06_threshold, rq06_provider_inf, rq06_routine_inf)
atomic_write_csv(rq06_inference, file.path(dir_rq06, "rq06_resident_bootstrap_inference.csv"))

rq06_secondary <- bind_rows(
  continuous_model_inference(adjacent, "vf_jaccard", "secondary_vf_jaccard", 720L),
  continuous_model_inference(adjacent, "abs_curated_count_change", "secondary_abs_curated_count_change", 730L)
)
atomic_write_csv(rq06_secondary, file.path(dir_rq06, "rq06_secondary_adjusted_contrasts.csv"))

rq06_groups <- adjacent %>%
  mutate(snp_context = if_else(close25, "SNP <=25", "SNP >25")) %>%
  group_by(snp_context) %>%
  summarise(
    n_pairs = n(),
    n_residents = n_distinct(Participant_id),
    median_days = median(days_between),
    any_vf_difference_n = sum(any_vf_difference),
    any_vf_difference_pct = 100 * mean(any_vf_difference),
    median_vf_jaccard = median(vf_jaccard, na.rm = TRUE),
    q25_vf_jaccard = quantile(vf_jaccard, 0.25, na.rm = TRUE),
    q75_vf_jaccard = quantile(vf_jaccard, 0.75, na.rm = TRUE),
    median_abs_curated_count_change = median(abs_curated_count_change, na.rm = TRUE),
    .groups = "drop"
  )
atomic_write_csv(rq06_groups, file.path(dir_rq06, "rq06_descriptive_summary.csv"))

rq06_pair_changes <- adjacent %>%
  select(Participant_id, key_from, key_to, date_from, date_to, days_between,
         status_from, status_to, event_type_from, event_type_to, TotalSNPs,
         close25, provider_ST_A, provider_ST_B, same_provider_ST,
         any_vf_difference, vf_jaccard, abs_curated_count_change,
         n_genes_from_only, n_genes_to_only, genes_from_only, genes_to_only)
atomic_write_csv(rq06_pair_changes, file.path(dir_rq06, "rq06_adjacent_pair_vf_changes.csv"))

gene_difference_counts <- bind_rows(
  tibble(
    pair_id = adjacent$pair_id,
    Participant_id = adjacent$Participant_id,
    snp_context = if_else(adjacent$close25, "SNP<=25", "SNP>25"),
    direction = "from_only",
    gene = genes_from_only
  ) %>% unnest_longer(gene),
  tibble(
    pair_id = adjacent$pair_id,
    Participant_id = adjacent$Participant_id,
    snp_context = if_else(adjacent$close25, "SNP<=25", "SNP>25"),
    direction = "to_only",
    gene = genes_to_only
  ) %>% unnest_longer(gene)
) %>%
  filter(!is.na(gene), nzchar(gene)) %>%
  count(snp_context, direction, gene, name = "n_pairs") %>%
  arrange(snp_context, direction, desc(n_pairs), gene) %>%
  mutate(note = "Descriptive endpoint-specific difference count; no gene-level hypothesis test performed.")
atomic_write_csv(gene_difference_counts, file.path(dir_rq06, "rq06_gene_difference_counts_descriptive.csv"))

p_rq06 <- ggplot(adjacent %>% mutate(snp_context = if_else(close25, "SNP <=25", "SNP >25")),
                  aes(snp_context, vf_jaccard, fill = snp_context)) +
  geom_boxplot(width = 0.55, outlier.shape = NA, alpha = 0.75) +
  geom_jitter(width = 0.15, alpha = 0.28, size = 1) +
  scale_fill_manual(values = c("SNP <=25" = "#0072B2", "SNP >25" = "#D55E00")) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(
    title = "VF profile similarity in adjacent within-resident isolates",
    subtitle = "Direct SNP distance; fresh SHA-bound VFDB calls from 532 Longcycler assemblies",
    x = NULL, y = "VF Jaccard similarity",
    caption = "Descriptive pair-level plot. Resident-bootstrap adjusted inference is reported separately."
  ) +
  theme_bw(base_size = 11) + theme(legend.position = "none")
atomic_ggsave(p_rq06, file.path(dir_rq06, "rq06_vf_jaccard_by_snp_context.png"), 7, 5)

prediction_grid <- expand_grid(
  days_between = seq(min(adjacent$days_between), max(adjacent$days_between), length.out = 120),
  close25 = c(FALSE, TRUE)
)
prediction_fit <- suppressWarnings(glm(
  any_vf_difference ~ close25 + splines::ns(days_between, df = 3),
  family = binomial(), data = adjacent
))
prediction_grid$predicted_probability <- predict(prediction_fit, prediction_grid, type = "response")
p_rq06_days <- ggplot(prediction_grid,
                       aes(days_between, predicted_probability, colour = close25)) +
  geom_line(linewidth = 1) +
  scale_colour_manual(values = c("FALSE" = "#D55E00", "TRUE" = "#0072B2"),
                      labels = c("SNP >25", "SNP <=25"), name = "Direct SNP context") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(
    title = "Adjusted probability of any binary VF difference",
    subtitle = "Natural spline for days between adjacent samples (df=3)",
    x = "Days between samples", y = "Standardized model probability",
    caption = "Curves are descriptive model estimates; uncertainty is resident-bootstrap and tabulated."
  ) +
  theme_bw(base_size = 11) + theme(legend.position = "bottom")
atomic_ggsave(p_rq06_days, file.path(dir_rq06, "rq06_adjusted_probability_by_days.png"), 7.5, 5)

# ==============================================================================
# RQ07: VF endpoints within UTI-event samples
# ==============================================================================
messagef("RQ07: analysing UTI-event samples and prespecified sensitivities")

event_samples <- episode_vf %>%
  filter(Event_type == "UTI_event") %>%
  arrange(Participant_id, collection_date, tp_lab)
event_anchor <- c(
  rows = nrow(event_samples),
  residents = n_distinct(event_samples$Participant_id),
  uti = sum(event_samples$UTI_Status == "UTI"),
  not_uti = sum(event_samples$UTI_Status == "Not_UTI")
)
event_expected <- c(
  rows = EXPECTED_EVENT_ROWS, residents = EXPECTED_EVENT_RESIDENTS,
  uti = EXPECTED_EVENT_UTI, not_uti = EXPECTED_EVENT_NOT_UTI
)
if (!identical(as.integer(event_anchor), as.integer(event_expected))) {
  block("rq07_event_anchor", paste(event_expected, collapse = ";"),
        paste(event_anchor, collapse = ";"),
        "The approved UTI-event sample denominator changed")
}
record_check("rq07_event_anchor", "PASS", paste(event_expected, collapse = ";"),
             paste(event_anchor, collapse = ";"), "Exact event-sample denominator")

median_difference <- function(df, endpoint) {
  x <- as.numeric(df[[endpoint]][df$UTI_Status == "UTI"])
  y <- as.numeric(df[[endpoint]][df$UTI_Status == "Not_UTI"])
  if (!length(x) || !length(y)) return(NA_real_)
  median(x, na.rm = TRUE) - median(y, na.rm = TRUE)
}

median_difference_inference <- function(df, endpoint, analysis, seed_offset = 0L) {
  point <- median_difference(df, endpoint)
  draws <- bootstrap_stat(
    df,
    function(z) median_difference(z, endpoint),
    seed = RQ_SEED + seed_offset
  )
  valid <- draws[is.finite(draws)]
  tail_p <- if (length(valid)) {
    min(1, 2 * min((sum(valid <= 0) + 1) / (length(valid) + 1),
                   (sum(valid >= 0) + 1) / (length(valid) + 1)))
  } else {
    NA_real_
  }
  bootstrap_interval(draws, point) %>%
    mutate(
      analysis = analysis,
      endpoint = endpoint,
      estimand = "median_UTI_minus_median_Not_UTI",
      n_rows = nrow(df),
      n_residents = n_distinct(df$Participant_id),
      n_uti = sum(df$UTI_Status == "UTI"),
      n_not_uti = sum(df$UTI_Status == "Not_UTI"),
      cluster_bootstrap_tail_p_exploratory = tail_p,
      .before = 1
    )
}

rq07_primary <- bind_rows(
  median_difference_inference(event_samples, "total_vf_count_curated",
                              "all_UTI_event_samples", 800L),
  median_difference_inference(event_samples, "five_marker_expec_count",
                              "all_UTI_event_samples", 810L)
) %>%
  mutate(
    p_holm_two_frozen_endpoints = p.adjust(cluster_bootstrap_tail_p_exploratory, method = "holm"),
    p_value_note = "Optional two-sided resident-cluster bootstrap tail areas; Holm correction covers only the two frozen endpoints."
  )
atomic_write_csv(rq07_primary, file.path(dir_rq07, "rq07_event_sample_inference.csv"))

earliest_event <- event_samples %>%
  group_by(Participant_id) %>%
  arrange(collection_date, tp_lab, .by_group = TRUE) %>%
  slice_head(n = 1L) %>%
  ungroup()
if (nrow(earliest_event) != EXPECTED_EVENT_RESIDENTS) {
  block("rq07_earliest_event", EXPECTED_EVENT_RESIDENTS, nrow(earliest_event),
        "Earliest-event sensitivity does not contain one row per event resident")
}
rq07_earliest <- bind_rows(
  median_difference_inference(earliest_event, "total_vf_count_curated", "earliest_event_per_resident", 820L),
  median_difference_inference(earliest_event, "five_marker_expec_count", "earliest_event_per_resident", 830L)
)
atomic_write_csv(rq07_earliest, file.path(dir_rq07, "rq07_earliest_event_sensitivity.csv"))

both_status_residents <- episode_vf %>%
  group_by(Participant_id) %>%
  summarise(n_status = n_distinct(UTI_Status), .groups = "drop") %>%
  filter(n_status == 2L) %>%
  pull(Participant_id)
if (length(both_status_residents) != EXPECTED_PAIRED_RESIDENTS) {
  block("rq07_paired_residents", EXPECTED_PAIRED_RESIDENTS, length(both_status_residents),
        "The approved within-resident paired sensitivity denominator changed")
}

nearest_pairs <- bind_rows(lapply(both_status_residents, function(pid) {
  z <- episode_vf %>% filter(Participant_id == pid)
  u <- z %>% filter(UTI_Status == "UTI")
  n <- z %>% filter(UTI_Status == "Not_UTI")
  grid <- tidyr::crossing(iu = seq_len(nrow(u)), inot = seq_len(nrow(n))) %>%
    mutate(
      absolute_days = abs(as.numeric(u$collection_date[iu] - n$collection_date[inot])),
      uti_key = u$episode_key[iu],
      not_uti_key = n$episode_key[inot]
    ) %>%
    arrange(absolute_days, uti_key, not_uti_key) %>%
    slice_head(n = 1L)
  tibble(
    Participant_id = pid,
    uti_key = grid$uti_key,
    not_uti_key = grid$not_uti_key,
    absolute_days = grid$absolute_days,
    curated_delta_uti_minus_not = u$total_vf_count_curated[grid$iu] - n$total_vf_count_curated[grid$inot],
    expec_delta_uti_minus_not = u$five_marker_expec_count[grid$iu] - n$five_marker_expec_count[grid$inot]
  )
}))
if (nrow(nearest_pairs) != EXPECTED_PAIRED_RESIDENTS || anyDuplicated(nearest_pairs$Participant_id)) {
  block("rq07_nearest_pair_rows", EXPECTED_PAIRED_RESIDENTS, nrow(nearest_pairs),
        "Nearest paired sensitivity did not produce one pair per eligible resident")
}
record_check("rq07_paired_residents", "PASS", EXPECTED_PAIRED_RESIDENTS, nrow(nearest_pairs),
             "One nearest UTI/Not_UTI pair per resident")

paired_delta_inference <- function(df, delta_col, endpoint, seed_offset) {
  statistic <- function(z) median(as.numeric(z[[delta_col]]), na.rm = TRUE)
  point <- statistic(df)
  draws <- bootstrap_stat(df, statistic, seed = RQ_SEED + seed_offset)
  bootstrap_interval(draws, point) %>%
    mutate(
      analysis = "nearest_within_resident_pair",
      endpoint = endpoint,
      estimand = "median_within_resident_UTI_minus_Not_UTI",
      n_pairs = nrow(df),
      n_positive = sum(df[[delta_col]] > 0),
      n_zero = sum(df[[delta_col]] == 0),
      n_negative = sum(df[[delta_col]] < 0),
      .before = 1
    )
}
rq07_paired <- bind_rows(
  paired_delta_inference(nearest_pairs, "curated_delta_uti_minus_not",
                         "total_vf_count_curated", 840L),
  paired_delta_inference(nearest_pairs, "expec_delta_uti_minus_not",
                         "five_marker_expec_count", 850L)
)
atomic_write_csv(nearest_pairs, file.path(dir_rq07, "rq07_nearest_paired_resident_values.csv"))
atomic_write_csv(rq07_paired, file.path(dir_rq07, "rq07_nearest_paired_sensitivity.csv"))
atomic_write_csv(
  event_samples %>%
    select(Participant_id, tp_lab, episode_key, collection_date, UTI_Status,
           total_vf_count_curated, five_marker_expec_count,
           expec_pap, expec_sfa_foc, expec_afa_dra, expec_iutA, expec_kpsM),
  file.path(dir_rq07, "rq07_uti_event_samples.csv")
)

p_rq07 <- ggplot(event_samples,
                  aes(UTI_Status, total_vf_count_curated, fill = UTI_Status)) +
  geom_boxplot(width = 0.5, outlier.shape = NA, alpha = 0.72) +
  geom_jitter(width = 0.12, size = 1.8, alpha = 0.68) +
  scale_fill_manual(values = c("Not_UTI" = "#0072B2", "UTI" = "#D55E00")) +
  labs(
    title = "Curated VF burden within UTI-event samples",
    subtitle = sprintf("%d samples from %d residents; fresh SHA-bound VFDB calls",
                       nrow(event_samples), n_distinct(event_samples$Participant_id)),
    x = "Primary status", y = "Curated detected VF genes",
    caption = "The event-sample comparison is descriptive; resident-cluster bootstrap intervals are tabulated."
  ) +
  theme_bw(base_size = 11) + theme(legend.position = "none")
atomic_ggsave(p_rq07, file.path(dir_rq07, "rq07_curated_vf_burden_event_samples.png"), 6.5, 5)

paired_plot <- nearest_pairs %>%
  select(Participant_id, uti_key, not_uti_key) %>%
  pivot_longer(c(uti_key, not_uti_key), names_to = "status_key", values_to = "episode_key") %>%
  mutate(UTI_Status = if_else(status_key == "uti_key", "UTI", "Not_UTI")) %>%
  left_join(episode_vf %>% select(episode_key, total_vf_count_curated),
            by = "episode_key", relationship = "many-to-one")
p_rq07_pair <- ggplot(paired_plot,
                       aes(UTI_Status, total_vf_count_curated, group = Participant_id)) +
  geom_line(alpha = 0.45, colour = "grey45") +
  geom_point(aes(colour = UTI_Status), size = 2.2) +
  scale_colour_manual(values = c("Not_UTI" = "#0072B2", "UTI" = "#D55E00")) +
  labs(
    title = "Nearest within-resident UTI/Not_UTI sensitivity",
    subtitle = sprintf("One nearest dated pair for each of %d residents", nrow(nearest_pairs)),
    x = NULL, y = "Curated detected VF genes",
    caption = "Lines identify paired observations but resident identifiers are not displayed."
  ) +
  theme_bw(base_size = 11) + theme(legend.position = "none")
atomic_ggsave(p_rq07_pair, file.path(dir_rq07, "rq07_nearest_paired_curated_burden.png"), 6.5, 5)

# ==============================================================================
# RQ08: Provider ST as a proxy for direct <=25-SNP closeness
# ==============================================================================
messagef("RQ08: evaluating provider-only ST and continuous genomic proxies")

provider_pairs <- pair_dat %>% filter(provider_ST_both)
if (!nrow(provider_pairs) || anyNA(provider_pairs$same_provider_ST)) {
  block("rq08_provider_pair_set", "non-empty complete provider-only pair set",
        nrow(provider_pairs), "No complete provider-only ST pair universe is available")
}
record_check("rq08_provider_pair_set", "PASS", "provider-only both endpoints",
             nrow(provider_pairs), "Local fallback and missing ST calls excluded")

diagnostic_metrics <- function(df, threshold) {
  actual <- df$TotalSNPs <= threshold
  test <- as.logical(df$same_provider_ST)
  keep <- !is.na(actual) & !is.na(test)
  actual <- actual[keep]
  test <- test[keep]
  tp <- sum(test & actual)
  fp <- sum(test & !actual)
  tn <- sum(!test & !actual)
  fn <- sum(!test & actual)
  sensitivity <- safe_ratio(tp, tp + fn)
  specificity <- safe_ratio(tn, tn + fp)
  c(
    sensitivity = sensitivity,
    specificity = specificity,
    ppv = safe_ratio(tp, tp + fp),
    npv = safe_ratio(tn, tn + fn),
    balanced_accuracy = mean(c(sensitivity, specificity), na.rm = TRUE)
  )
}

diagnostic_inference <- function(df, threshold, scope, seed_offset) {
  metric_names <- c("sensitivity", "specificity", "ppv", "npv", "balanced_accuracy")
  point <- diagnostic_metrics(df, threshold)
  draws <- bootstrap_vector(
    df,
    function(z) diagnostic_metrics(z, threshold),
    metric_names,
    seed = RQ_SEED + seed_offset
  )
  map_dfr(metric_names, function(metric) {
    bootstrap_interval(draws[, metric], point[[metric]]) %>%
      mutate(metric = metric, .before = 1)
  }) %>%
    mutate(
      scope = scope,
      reference_snp_threshold = threshold,
      n_pairs = nrow(df),
      n_residents = n_distinct(df$Participant_id),
      n_reference_close = sum(df$TotalSNPs <= threshold),
      n_same_provider_ST = sum(df$same_provider_ST),
      test_definition = "same provider SeqSphere ST; both endpoints provider_qc95 Longcycler",
      .before = 1
    )
}

adjacent_pair_ids <- adjacent$pair_id
rq08_diagnostic <- bind_rows(
  map_dfr(c(10L, 25L, 50L), function(threshold) {
    diagnostic_inference(provider_pairs, threshold, "all_within_resident_pairs",
                         900L + threshold)
  }),
  map_dfr(c(10L, 25L, 50L), function(threshold) {
    diagnostic_inference(provider_pairs %>% filter(pair_id %in% adjacent_pair_ids),
                         threshold, "adjacent_pairs_only", 950L + threshold)
  })
)
atomic_write_csv(rq08_diagnostic, file.path(dir_rq08, "rq08_provider_st_diagnostic_metrics.csv"))

auc_inference <- function(df, score_col, direction, label, seed_offset) {
  z <- df %>%
    mutate(.score = as.numeric(.data[[score_col]]), close25 = TotalSNPs <= 25) %>%
    filter(is.finite(.score))
  if (direction == "lower_is_close") z$.score <- -z$.score
  point <- auc_rank(z$close25, z$.score)
  draws <- bootstrap_stat(
    z,
    function(b) auc_rank(b$close25, b$.score),
    seed = RQ_SEED + seed_offset
  )
  bootstrap_interval(draws, point) %>%
    mutate(
      predictor = label,
      score_column = score_col,
      direction = direction,
      n_pairs = nrow(z),
      n_residents = n_distinct(z$Participant_id),
      reference = "direct SNP <=25",
      .before = 1
    )
}

rq08_auc <- bind_rows(
  auc_inference(pair_dat, "MashDistance", "lower_is_close", "Mash distance", 1000L),
  auc_inference(pair_dat, "fresh_vf_jaccard", "higher_is_close", "Fresh VF Jaccard", 1010L)
)

fit_loor_composite <- function(df, reference_threshold = 25L) {
  z <- df %>%
    mutate(
      close_reference = as.integer(TotalSNPs <= reference_threshold),
      same_st_num = as.integer(same_provider_ST)
    ) %>%
    filter(provider_ST_both, is.finite(MashDistance), is.finite(fresh_vf_jaccard))
  ids <- unique(z$Participant_id)
  predictions <- bind_rows(lapply(ids, function(held_out) {
    train <- z %>% filter(Participant_id != held_out)
    test <- z %>% filter(Participant_id == held_out)
    mash_mean <- mean(train$MashDistance)
    mash_sd <- sd(train$MashDistance)
    vf_mean <- mean(train$fresh_vf_jaccard)
    vf_sd <- sd(train$fresh_vf_jaccard)
    if (!is.finite(mash_sd) || mash_sd == 0) mash_sd <- 1
    if (!is.finite(vf_sd) || vf_sd == 0) vf_sd <- 1
    train <- train %>%
      mutate(mash_z = (MashDistance - mash_mean) / mash_sd,
             vf_z = (fresh_vf_jaccard - vf_mean) / vf_sd)
    test <- test %>%
      mutate(mash_z = (MashDistance - mash_mean) / mash_sd,
             vf_z = (fresh_vf_jaccard - vf_mean) / vf_sd)
    fit <- tryCatch(
      suppressWarnings(glm(close_reference ~ same_st_num + mash_z + vf_z,
                           data = train, family = binomial(),
                           control = glm.control(maxit = 100))),
      error = function(e) NULL
    )
    pred <- if (is.null(fit)) rep(NA_real_, nrow(test)) else
      suppressWarnings(as.numeric(predict(fit, newdata = test, type = "response")))
    test %>%
      transmute(Participant_id, pair_id, reference_threshold,
                outcome = close_reference, predicted_probability = pred)
  }))
  predictions
}

composite_specs <- tribble(
  ~scope, ~threshold,
  "all_within_resident_pairs", 10L,
  "all_within_resident_pairs", 25L,
  "all_within_resident_pairs", 50L,
  "adjacent_pairs_only", 25L
)
composite_predictions <- pmap_dfr(composite_specs, function(scope, threshold) {
  z <- provider_pairs
  if (scope == "adjacent_pairs_only") z <- z %>% filter(pair_id %in% adjacent_pair_ids)
  fit_loor_composite(z, threshold) %>% mutate(scope = scope)
})
atomic_write_csv(composite_predictions, file.path(dir_rq08, "rq08_loor_composite_predictions.csv"))

composite_summary <- composite_predictions %>%
  group_by(scope, reference_threshold) %>%
  group_modify(~ {
    z <- .x %>% filter(is.finite(predicted_probability))
    point_auc <- auc_rank(z$outcome, z$predicted_probability)
    draws <- bootstrap_stat(
      z,
      function(b) auc_rank(b$outcome, b$predicted_probability),
      seed = RQ_SEED + 1100L + unique(.y$reference_threshold)
    )
    bootstrap_interval(draws, point_auc) %>%
      mutate(
        n_pairs = nrow(z),
        n_residents = n_distinct(z$Participant_id),
        brier_score = mean((z$predicted_probability - z$outcome)^2),
        sensitivity_at_0_5 = safe_ratio(sum(z$predicted_probability >= 0.5 & z$outcome == 1),
                                        sum(z$outcome == 1)),
        specificity_at_0_5 = safe_ratio(sum(z$predicted_probability < 0.5 & z$outcome == 0),
                                        sum(z$outcome == 0)),
        model = "leave-one-resident-out logistic: same provider ST + scaled Mash distance + scaled fresh VF Jaccard"
      )
  }) %>%
  ungroup() %>%
  rename(auc = estimate, auc_ci_lower = ci_lower, auc_ci_upper = ci_upper)
atomic_write_csv(composite_summary, file.path(dir_rq08, "rq08_loor_composite_summary.csv"))

composite_primary <- composite_predictions %>%
  filter(scope == "all_within_resident_pairs", reference_threshold == 25L,
         is.finite(predicted_probability))
rq08_auc <- bind_rows(
  rq08_auc,
  bootstrap_interval(
    bootstrap_stat(composite_primary,
                   function(z) auc_rank(z$outcome, z$predicted_probability),
                   seed = RQ_SEED + 1110L),
    auc_rank(composite_primary$outcome, composite_primary$predicted_probability)
  ) %>%
    mutate(
      predictor = "LOOR composite logistic",
      score_column = "predicted_probability",
      direction = "higher_is_close",
      n_pairs = nrow(composite_primary),
      n_residents = n_distinct(composite_primary$Participant_id),
      reference = "direct SNP <=25",
      .before = 1
    )
)
atomic_write_csv(rq08_auc, file.path(dir_rq08, "rq08_roc_auc_resident_bootstrap.csv"))

roc_coordinates <- function(outcome, score, predictor, points = 201L) {
  keep <- !is.na(outcome) & is.finite(score)
  y <- as.integer(outcome[keep])
  s <- score[keep]
  thresholds <- unique(as.numeric(quantile(s, probs = seq(0, 1, length.out = points),
                                           na.rm = TRUE, names = FALSE)))
  bind_rows(lapply(sort(thresholds, decreasing = TRUE), function(thr) {
    pred <- s >= thr
    tibble(
      predictor = predictor,
      threshold = thr,
      sensitivity = safe_ratio(sum(pred & y == 1L), sum(y == 1L)),
      false_positive_rate = 1 - safe_ratio(sum(!pred & y == 0L), sum(y == 0L))
    )
  }))
}

roc_coords <- bind_rows(
  roc_coordinates(pair_dat$TotalSNPs <= 25, -pair_dat$MashDistance, "Mash distance"),
  roc_coordinates(pair_dat$TotalSNPs <= 25, pair_dat$fresh_vf_jaccard, "Fresh VF Jaccard"),
  roc_coordinates(composite_primary$outcome, composite_primary$predicted_probability,
                  "LOOR composite logistic")
)
atomic_write_csv(roc_coords, file.path(dir_rq08, "rq08_roc_coordinates.csv"))

p_rq08 <- ggplot(roc_coords, aes(false_positive_rate, sensitivity, colour = predictor)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey65") +
  geom_path(linewidth = 0.9) +
  coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
  labs(
    title = "Proxy discrimination for direct within-resident SNP closeness",
    subtitle = "Operational reference: direct assembly SNP distance <=25",
    x = "False-positive rate", y = "Sensitivity", colour = "Predictor",
    caption = "Composite predictions are leave-one-resident-out; AUC uncertainty uses resident bootstrap."
  ) +
  theme_bw(base_size = 11) + theme(legend.position = "bottom")
atomic_ggsave(p_rq08, file.path(dir_rq08, "rq08_roc_curves.png"), 7, 6)

metric_plot <- rq08_diagnostic %>%
  filter(scope == "all_within_resident_pairs", reference_snp_threshold == 25L)
p_rq08_metrics <- ggplot(metric_plot, aes(metric, estimate)) +
  geom_col(fill = "#4C78A8", width = 0.68) +
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0.18) +
  coord_cartesian(ylim = c(0, 1)) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(
    title = "Provider-ST performance against direct SNP <=25",
    subtitle = sprintf("%d within-resident pairs with provider ST at both endpoints",
                       unique(metric_plot$n_pairs)),
    x = NULL, y = "Performance",
    caption = "Bars are point estimates; intervals use resident-cluster bootstrap."
  ) +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 25, hjust = 1))
atomic_ggsave(p_rq08_metrics, file.path(dir_rq08, "rq08_provider_st_performance.png"), 7, 5)

# Prespecified predicted-plasmid extension. The sourced module uses only the
# already validated 532/893/371 universes and blocks on partial script-09b/29
# outputs.
source(file.path(
  root, "scripts", "research_questions",
  "plasmid_mechanism_addon_rq06_08.R"
))

# ----------------------------------------------------------------------------
# Reproducibility manifests and concise question summaries
# ----------------------------------------------------------------------------
input_manifest <- tibble(
  role = c("analysis_assembly_manifest", "clinical_status_map", "pairwise_metrics",
           "provider_mlst", "gene_module_map", "mob_episode_profiles",
           "episode_mechanism_profiles", "analysis_script",
           "plasmid_addon_script"),
  path = c(path_manifest, path_status, path_pairwise, path_mlst, path_module_map,
           path_mob_profiles, path_mechanism_profiles,
           file.path(root, "scripts", "research_questions", "run_rq06_08.R"),
           file.path(
             root, "scripts", "research_questions",
             "plasmid_mechanism_addon_rq06_08.R"
           ))
) %>%
  mutate(
    exists = file.exists(path),
    file_size = ifelse(exists, file.size(path), NA_real_),
    modified_time = ifelse(exists, as.character(file.info(path)$mtime), NA_character_),
    sha256 = vapply(path, file_sha256_or_na, character(1))
  )
atomic_write_csv(input_manifest, file.path(dir_inputs, "source_input_manifest.csv"))

method_manifest <- tibble(
  parameter = c(
    "expected_episodes", "expected_all_pairs", "expected_adjacent_pairs",
    "vfdb_min_identity_pct", "vfdb_min_coverage_pct", "same_strain_snp_threshold",
    "bootstrap_unit", "bootstrap_reps", "random_seed",
    "provider_mlst_min_good_targets_pct", "provider_st_policy",
    "abricate_version", "vfdb_database_record", "abricate_database_list_sha256",
    "vfdb_cache_context_sha256", "plasmid_interpretation",
    "plasmid_composite_selection_rule"
  ),
  value = c(
    EXPECTED_EPISODES, EXPECTED_ALL_PAIRS, EXPECTED_ADJACENT_PAIRS,
    VFDB_MIN_ID, VFDB_MIN_COV, 25,
    "resident/Participant_id", BOOT_REPS, RQ_SEED,
    PROVIDER_MLST_MIN_GOOD_TARGETS,
    "provider_qc95 call key/path-linked to the selected Longcycler episode; local fallback excluded",
    abricate_version, vfdb_line[1], db_list_sha256, cache_context,
    "assembly-based predicted plasmid context; not confirmed circularity, transfer or transmission",
    "add plasmid metric only if resident-bootstrap paired delta-AUC 95% interval is above zero"
  )
)
atomic_write_csv(method_manifest, file.path(dir_inputs, "method_manifest.csv"))

atomic_write_lines(c(
  "# RQ06 — Adjacent-pair VF stability",
  "",
  sprintf("- Analysis universe: %d date-ordered adjacent pairs from %d residents.",
          nrow(adjacent), n_distinct(adjacent$Participant_id)),
  "- Primary exposure: direct assembly SNP distance <=25.",
  "- Primary outcome: any binary fresh-VFDB profile difference.",
  "- Inference: resident-cluster bootstrap GLM with a natural spline for days between samples.",
  "- Plasmid extension: any corrected replicon-profile change and any MOB primary-cluster change; 10/25/50-SNP, both-empty and high-confidence sensitivities.",
  "- Gene difference tables are descriptive only; no broad gene-level tests are performed."
), file.path(dir_rq06, "README.md"))

atomic_write_lines(c(
  "# RQ07 — VF endpoints within UTI-event samples",
  "",
  sprintf("- Primary event-sample universe: %d samples from %d residents (%d UTI, %d Not_UTI).",
          nrow(event_samples), n_distinct(event_samples$Participant_id),
          sum(event_samples$UTI_Status == "UTI"), sum(event_samples$UTI_Status == "Not_UTI")),
  "- Frozen endpoints: curated VF burden and five-marker ExPEC count.",
  "- Exploratory plasmid endpoints: predicted plasmid count and plasmid-binned informative VF/AMR burden only.",
  "- Inference: resident-cluster bootstrap median differences; two optional tail-area p-values use Holm correction.",
  "- Sensitivities: earliest UTI-event per resident and nearest within-resident UTI/Not_UTI pair."
), file.path(dir_rq07, "README.md"))

atomic_write_lines(c(
  "# RQ08 — Provider-ST proxy performance",
  "",
  sprintf("- Direct pair universe: %d path/SHA-validated within-resident pairs.", nrow(pair_dat)),
  sprintf("- Provider-only evaluable pairs: %d; local fallback labels are excluded.", nrow(provider_pairs)),
  "- Operational reference: direct assembly SNP distance <=25.",
  "- Metrics: sensitivity, specificity, PPV, NPV, balanced accuracy, and resident-bootstrap intervals.",
  "- Continuous proxies: Mash distance and fresh VF Jaccard. Composite predictions use leave-one-resident-out logistic regression.",
  "- Plasmid proxies: corrected replicon and MOB-cluster Jaccard; a plasmid metric enters the selected composite only if its paired resident-bootstrap delta-AUC interval supports improvement.",
  "- Sensitivities: SNP cutoffs 10/25/50 and adjacent-pair restriction."
), file.path(dir_rq08, "README.md"))

write_analysis_status <- function(path, rq, detail) {
  atomic_write_csv(
    tibble(
      research_question = rq,
      status = "complete",
      reason = "analysis_completed",
      detail = detail,
      run_timestamp_utc = format(Sys.time(), tz = "UTC", usetz = TRUE)
    ),
    file.path(path, "analysis_status.csv")
  )
}

write_analysis_status(
  dir_rq06,
  "RQ06",
  sprintf("%d adjacent Longcycler-linked direct pairs", nrow(adjacent))
)
write_analysis_status(
  dir_rq07,
  "RQ07",
  sprintf("%d Longcycler-linked UTI-event samples from %d residents",
          nrow(event_samples), n_distinct(event_samples$Participant_id))
)
write_analysis_status(
  dir_rq08,
  "RQ08",
  sprintf("%d Longcycler-linked direct pairs; %d provider-ST evaluable pairs",
          nrow(pair_dat), nrow(provider_pairs))
)

record_check("rq06_outputs", "PASS", "models/tables/plots", "written",
             "No gene-level hypothesis testing")
record_check("rq07_outputs", "PASS", "32/29 event anchor and paired sensitivities", "written",
             "Only two frozen endpoints")
record_check("rq08_outputs", "PASS", "893 pairs and provider-only metrics", "written",
             "Direct SNP reference; local fallback excluded")
write_run_state("COMPLETE", "RQ06-RQ08 completed with all provenance checks passing")

messagef("RQ06-RQ08 complete. Outputs: %s", dir_rq)

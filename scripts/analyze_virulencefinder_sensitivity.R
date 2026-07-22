#!/usr/bin/env Rscript

# =============================================================================
# CGE VirulenceFinder sensitivity analysis for RQ06-RQ08
# =============================================================================
# Scientific and operational settings are read from the commented file
# config/virulencefinder_sensitivity.toml through the batch runner's validated
# effective_config.json. Change settings there, not in this script.
#
# This script never writes inside results/research_questions/. The completed
# ABRicate/VFDB analysis remains primary; every output here is explicitly a CGE
# VirulenceFinder sensitivity or concordance result.
# =============================================================================

options(stringsAsFactors = FALSE, warn = 1)

required_packages <- c("dplyr", "readr", "tidyr", "purrr", "jsonlite", "digest", "tibble")
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
  library(tibble)
})

root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
if (!file.exists(file.path(root, "config", "virulencefinder_sensitivity.toml"))) {
  stop("Run this script from the rUTIs project root.", call. = FALSE)
}

output <- file.path(root, "results", "virulencefinder_cge_3_2_1")
out_rq <- file.path(output, "rq_sensitivity")
out_rq06 <- file.path(out_rq, "RQ06")
out_rq07 <- file.path(out_rq, "RQ07")
out_rq08 <- file.path(out_rq, "RQ08")
for (directory in c(output, out_rq, out_rq06, out_rq07, out_rq08)) {
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
}

final_marker <- file.path(output, "RUN_COMPLETE.txt")
if (file.exists(final_marker)) file.remove(final_marker)

timestamp_utc <- function() format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")
sha256_file <- function(path) digest::digest(file = path, algo = "sha256", serialize = FALSE)
as_bool <- function(x) tolower(trimws(as.character(x))) %in% c("true", "t", "1", "yes", "y")
safe_ratio <- function(n, d) ifelse(d > 0, n / d, NA_real_)

atomic_replace <- function(tmp, path) {
  if (file.exists(path) && !file.remove(path)) stop("Could not replace output: ", path)
  if (!file.rename(tmp, path)) stop("Could not publish output: ", path)
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

checks <- tibble(
  check_id = character(), status = character(), expected = character(),
  observed = character(), detail = character()
)

check <- function(id, condition, expected, observed, detail = "") {
  checks <<- bind_rows(checks, tibble(
    check_id = id,
    status = if (isTRUE(condition)) "PASS" else "FAIL",
    expected = as.character(expected), observed = as.character(observed), detail = detail
  ))
  if (!isTRUE(condition)) {
    atomic_write_csv(checks, file.path(out_rq, "acceptance_checks.csv"))
    stop("BLOCKED [", id, "]: ", detail, call. = FALSE)
  }
  invisible(TRUE)
}

require_columns <- function(data, columns, label) {
  missing <- setdiff(columns, names(data))
  check(paste0(label, "_columns"), !length(missing), paste(columns, collapse = ";"),
        paste(names(data), collapse = ";"),
        if (length(missing)) paste("Missing:", paste(missing, collapse = ", ")) else "Required columns present")
}

tree_fingerprint <- function(directory) {
  paths <- sort(list.files(directory, recursive = TRUE, full.names = TRUE, all.files = TRUE,
                           include.dirs = FALSE, no.. = TRUE))
  records <- vapply(paths, function(path) {
    paste(sub(paste0("^", directory, "/?"), "", path), file.info(path)$size, sha256_file(path), sep = "|")
  }, character(1))
  digest::digest(paste(records, collapse = "\n"), algo = "sha256", serialize = FALSE)
}

auc_rank <- function(outcome, score) {
  keep <- !is.na(outcome) & is.finite(score)
  y <- as.integer(outcome[keep])
  score <- as.numeric(score[keep])
  n1 <- sum(y == 1L)
  n0 <- sum(y == 0L)
  if (!n1 || !n0) return(NA_real_)
  (sum(rank(score, ties.method = "average")[y == 1L]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

bootstrap_interval <- function(draws, estimate) {
  good <- as.numeric(draws)[is.finite(as.numeric(draws))]
  alpha <- (1 - CONFIDENCE_LEVEL) / 2
  tibble(
    estimate = estimate,
    ci_lower = if (length(good)) unname(quantile(good, alpha, na.rm = TRUE)) else NA_real_,
    ci_upper = if (length(good)) unname(quantile(good, 1 - alpha, na.rm = TRUE)) else NA_real_,
    bootstrap_reps_requested = length(draws), bootstrap_reps_valid = length(good)
  )
}

# Generates the complete resident-resampling index matrix from the fixed seed
# before optional forking. Results therefore do not depend on worker scheduling.
bootstrap_apply <- function(data, statistic, reps, seed, workers) {
  ids <- unique(as.character(data$Participant_id))
  row_groups <- split(seq_len(nrow(data)), factor(as.character(data$Participant_id), levels = ids))
  set.seed(seed)
  sampled <- matrix(sample.int(length(ids), length(ids) * reps, replace = TRUE), nrow = length(ids))
  one <- function(index) {
    groups <- row_groups[sampled[, index]]
    rows <- unlist(groups, use.names = FALSE)
    z <- data[rows, , drop = FALSE]
    z$.bootstrap_cluster <- rep(seq_along(groups), lengths(groups))
    tryCatch(as.numeric(statistic(z)), error = function(error) NA_real_)
  }
  if (workers > 1L && .Platform$OS.type != "windows") {
    values <- parallel::mclapply(seq_len(reps), one, mc.cores = workers,
                                 mc.preschedule = TRUE, mc.set.seed = FALSE)
  } else {
    values <- lapply(seq_len(reps), one)
  }
  width <- max(vapply(values, length, integer(1)))
  matrix(unlist(lapply(values, function(value) {
    if (length(value) == width) value else rep(NA_real_, width)
  })), nrow = width)[, seq_len(reps), drop = FALSE]
}

parse_date_multi <- function(x) {
  x <- as.character(x)
  out <- as.Date(rep(NA_character_, length(x)))
  for (format in c("%d/%m/%Y", "%d-%m-%Y", "%Y-%m-%d", "%Y/%m/%d")) {
    use <- is.na(out) & !is.na(x) & nzchar(trimws(x))
    out[use] <- as.Date(x[use], format = format)
  }
  out
}

effective_path <- file.path(output, "effective_config.json")
check("effective_config_exists", file.exists(effective_path), "present", effective_path,
      "The Python preflight must complete first")
effective <- jsonlite::fromJSON(effective_path, simplifyVector = TRUE)
cfg <- effective$config
cohort <- cfg$cohort
analysis <- cfg$analysis
profiles <- names(cfg$profiles)
expected_profiles <- c("web_default_id90_cov60", "matched_id80_cov80")
check("profile_contract", setequal(profiles, expected_profiles), paste(expected_profiles, collapse = ";"),
      paste(profiles, collapse = ";"), "Both prespecified profiles are required")

BOOT_REPS <- as.integer(analysis$bootstrap_reps)
RQ_SEED <- as.integer(analysis$seed)
SPLINE_DF <- as.integer(analysis$spline_degrees_freedom)
SNP_THRESHOLDS <- as.integer(analysis$snp_thresholds)
PRIMARY_SNP <- as.integer(analysis$primary_snp_threshold)
CONFIDENCE_LEVEL <- as.numeric(analysis$confidence_level)
BOOT_WORKERS <- as.integer(if (is.null(analysis$bootstrap_workers)) 1L else analysis$bootstrap_workers)
BOOT_WORKERS <- max(1L, min(BOOT_WORKERS, 8L))
check("final_bootstrap_contract", BOOT_REPS == 10000L && RQ_SEED == 20260712L && CONFIDENCE_LEVEL == 0.95,
      "10000 reps; seed 20260712; confidence 0.95",
      paste(BOOT_REPS, RQ_SEED, CONFIDENCE_LEVEL, sep = ";"),
      "Final sensitivity release uses the prespecified resampling contract")

primary_dir <- cfg$paths$primary_release_dir
primary_before <- tree_fingerprint(primary_dir)
check("primary_run_marker", file.exists(file.path(primary_dir, "RUN_COMPLETE.txt")) &&
        grepl("PASS", paste(readLines(file.path(primary_dir, "RUN_COMPLETE.txt"), warn = FALSE), collapse = " ")),
      "PASS marker", "present", "Primary ABRicate/VFDB release remains authoritative")
primary_checks <- read_csv(file.path(primary_dir, "final_contract_checks.csv"), show_col_types = FALSE)
check("primary_18_contract_checks", nrow(primary_checks) == 18L && all(as_bool(primary_checks$pass)),
      "18 PASS", paste(nrow(primary_checks), sum(as_bool(primary_checks$pass)), sep = ";"),
      "Primary denominator contract must be intact")
primary_status <- read_csv(file.path(primary_dir, "final_question_status.csv"), show_col_types = FALSE)
check("primary_rq_status", setequal(primary_status$research_question[primary_status$status == "complete"],
                                     sprintf("RQ%02d", 1:10)),
      "RQ01-RQ10 complete", paste(primary_status$research_question[primary_status$status == "complete"], collapse = ";"),
      "No sensitivity result is generated over an incomplete primary release")

batch_marker <- file.path(output, "BATCH_COMPLETE.txt")
check("batch_complete_marker", file.exists(batch_marker), "present", batch_marker,
      "All VirulenceFinder jobs must complete before RQ integration")
run_manifest <- read_csv(file.path(output, "run_manifest.csv"), show_col_types = FALSE)
require_columns(run_manifest, c("profile", "Assembly_ID", "episode_key", "fasta_sha256", "status"), "run_manifest")
expected_jobs <- as.integer(cohort$expected_assemblies) * length(expected_profiles)
check("run_manifest_jobs", nrow(run_manifest) == expected_jobs && all(run_manifest$status == "success"),
      paste(expected_jobs, "successful"), paste(nrow(run_manifest), sum(run_manifest$status == "success"), sep = ";"),
      "Successful zero-hit jobs are valid successes")
run_by_profile <- run_manifest %>% count(profile, name = "n")
check("run_manifest_profile_counts", all(run_by_profile$n == as.integer(cohort$expected_assemblies)) &&
        setequal(run_by_profile$profile, expected_profiles), "532 per profile",
      paste(run_by_profile$profile, run_by_profile$n, sep = "=", collapse = ";"), "No profile may be incomplete")

manifest <- read_csv(cfg$paths$manifest, show_col_types = FALSE, progress = FALSE) %>%
  mutate(
    Participant_id = as.character(Participant_id), tp_lab = as.character(tp_lab),
    episode_key = paste(Participant_id, tp_lab, sep = "||"),
    fasta_sha256 = tolower(fasta_sha256),
    collection_date = parse_date_multi(Collection_Date)
  )
check("manifest_anchor", nrow(manifest) == as.integer(cohort$expected_assemblies) &&
        n_distinct(manifest$Participant_id) == as.integer(cohort$expected_residents),
      "532 genomes; 161 residents", paste(nrow(manifest), n_distinct(manifest$Participant_id), sep = ";"),
      "Authoritative Longcycler denominator")

metrics <- read_csv(file.path(output, "episode_metrics.csv"), show_col_types = FALSE, progress = FALSE) %>%
  mutate(Participant_id = as.character(Participant_id), fasta_sha256 = tolower(fasta_sha256))
require_columns(metrics, c("profile", "Participant_id", "episode_key", "fasta_sha256",
                           "distinct_cge_gene_family_count", "CGE_available_ExPEC_group_count"), "episode_metrics")
check("episode_metrics_rows", nrow(metrics) == expected_jobs, expected_jobs, nrow(metrics),
      "One metric row per assembly and profile")

pa <- list()
metadata_columns <- c("Participant_id", "tp_lab", "episode_key", "Assembly_ID", "fasta_sha256",
                      "UTI_Status", "Event_type", "full_path")
for (profile in expected_profiles) {
  path <- file.path(output, paste0("presence_absence_", profile, ".csv"))
  z <- read_csv(path, show_col_types = FALSE, progress = FALSE) %>%
    mutate(Participant_id = as.character(Participant_id), fasta_sha256 = tolower(fasta_sha256))
  genes <- setdiff(names(z), metadata_columns)
  anchor <- manifest %>% select(episode_key, manifest_sha = fasta_sha256)
  validated <- z %>% left_join(anchor, by = "episode_key", relationship = "many-to-one")
  check(paste0("pa_", profile), nrow(z) == as.integer(cohort$expected_assemblies) &&
          length(genes) == as.integer(cohort$expected_combined_gene_families) &&
          !anyDuplicated(z$episode_key) && all(z$fasta_sha256 == validated$manifest_sha) &&
          all(as.matrix(z[, genes, drop = FALSE]) %in% c(0, 1)),
        "532 rows; 681 families; binary; SHA-bound",
        paste(nrow(z), length(genes), sum(z$fasta_sha256 == validated$manifest_sha), sep = ";"),
        "Pinned stx2 overlap gives 681 unique gene-family columns")
  pa[[profile]] <- list(data = z, genes = genes)
}

primary_episode <- read_csv(cfg$paths$primary_episode_metrics, show_col_types = FALSE, progress = FALSE) %>%
  mutate(Participant_id = as.character(Participant_id), fasta_sha256 = tolower(fasta_sha256),
         collection_date = as.Date(collection_date))
adjacent_path <- cfg$paths$adjacent_pairs
adjacent_base <- read_csv(adjacent_path, show_col_types = FALSE, progress = FALSE) %>%
  mutate(Participant_id = as.character(Participant_id))
direct <- read_csv(cfg$paths$direct_pairs, show_col_types = FALSE, progress = FALSE) %>%
  mutate(Participant_id = as.character(Participant_id),
         fasta_sha256_A = tolower(fasta_sha256_A), fasta_sha256_B = tolower(fasta_sha256_B))
require_columns(direct, c("Participant_id", "key_A", "key_B", "pair_id", "TotalSNPs", "MashDistance",
                          "provider_ST_both", "same_provider_ST", "fasta_sha256_A", "fasta_sha256_B"), "direct_pairs")

manifest_sha <- setNames(manifest$fasta_sha256, manifest$episode_key)
direct_hash_ok <- direct$fasta_sha256_A == manifest_sha[direct$key_A] &
  direct$fasta_sha256_B == manifest_sha[direct$key_B]
check("direct_pair_anchor", nrow(direct) == as.integer(cohort$expected_direct_pairs) &&
        !anyDuplicated(direct$pair_id) && all(direct_hash_ok),
      "893 unique SHA-bound pairs", paste(nrow(direct), sum(direct_hash_ok), sep = ";"),
      "CGE Jaccard is evaluated against the approved direct SNP endpoints")
check("adjacent_pair_anchor", nrow(adjacent_base) == as.integer(cohort$expected_adjacent_pairs) &&
        n_distinct(adjacent_base$Participant_id) == as.integer(cohort$expected_adjacent_residents) &&
        all(adjacent_base$pair_id %in% direct$pair_id),
      "371 pairs; 139 residents", paste(nrow(adjacent_base), n_distinct(adjacent_base$Participant_id), sep = ";"),
      "Exact approved adjacent-pair set")

resident_key <- manifest %>% distinct(Participant_id) %>% arrange(Participant_id) %>%
  mutate(resident_label = sprintf("RES%03d", row_number()))
episode_key_map <- manifest %>% arrange(Participant_id, collection_date, tp_lab) %>%
  mutate(episode_label = sprintf("EPI%04d", row_number())) %>%
  select(episode_key, episode_label)

pair_metrics <- function(pair_table, key_a, key_b, profile) {
  z <- pa[[profile]]$data
  genes <- pa[[profile]]$genes
  ia <- match(pair_table[[key_a]], z$episode_key)
  ib <- match(pair_table[[key_b]], z$episode_key)
  check(paste0("pair_endpoint_matrix_", profile, "_", key_a), !anyNA(ia) && !anyNA(ib),
        "all endpoints present", sum(!is.na(ia) & !is.na(ib)), "Every pair endpoint must occur in the 532-row matrix")
  a <- as.matrix(z[ia, genes, drop = FALSE])
  b <- as.matrix(z[ib, genes, drop = FALSE])
  storage.mode(a) <- "integer"; storage.mode(b) <- "integer"
  intersection <- rowSums(a == 1L & b == 1L)
  union <- rowSums(a == 1L | b == 1L)
  from_only <- lapply(seq_len(nrow(a)), function(i) genes[a[i, ] == 1L & b[i, ] == 0L])
  to_only <- lapply(seq_len(nrow(a)), function(i) genes[a[i, ] == 0L & b[i, ] == 1L])
  pair_table %>% mutate(
    any_gene_family_difference = rowSums(a != b) > 0L,
    cge_jaccard = ifelse(union > 0L, intersection / union, 1),
    cge_intersection_n = intersection, cge_union_n = union,
    absolute_gene_family_count_change = abs(rowSums(a) - rowSums(b)),
    n_families_from_only = lengths(from_only), n_families_to_only = lengths(to_only),
    families_from_only = vapply(from_only, paste, collapse = ";", character(1)),
    families_to_only = vapply(to_only, paste, collapse = ";", character(1))
  )
}

model_binary <- function(data, exposure) {
  z <- data %>% mutate(.exposure = as.integer(.data[[exposure]])) %>%
    filter(!is.na(any_gene_family_difference), !is.na(.exposure), is.finite(days_between))
  if (nrow(z) < 20L || n_distinct(z$.exposure) != 2L ||
      n_distinct(z$any_gene_family_difference) != 2L) return(c(log_or = NA_real_, adjusted_rd = NA_real_))
  formula <- as.formula(paste0("any_gene_family_difference ~ .exposure + splines::ns(days_between, df=", SPLINE_DF, ")"))
  fit <- suppressWarnings(glm(formula, family = binomial(), data = z, control = glm.control(maxit = 100)))
  nd1 <- z; nd0 <- z; nd1$.exposure <- 1L; nd0$.exposure <- 0L
  c(log_or = unname(coef(fit)[[".exposure"]]),
    adjusted_rd = mean(predict(fit, nd1, type = "response") - predict(fit, nd0, type = "response"), na.rm = TRUE))
}

binary_inference <- function(data, exposure, analysis_label, profile, seed_offset) {
  point <- model_binary(data, exposure)
  draws <- bootstrap_apply(data, function(z) model_binary(z, exposure), BOOT_REPS,
                           RQ_SEED + seed_offset, BOOT_WORKERS)
  n_exposed_value <- sum(data[[exposure]] %in% TRUE, na.rm = TRUE)
  bind_rows(
    bootstrap_interval(exp(draws[1, ]), exp(point[[1]])) %>% mutate(estimand = "odds_ratio"),
    bootstrap_interval(draws[2, ], point[[2]]) %>% mutate(estimand = "adjusted_risk_difference")
  ) %>% mutate(
    profile = .env$profile, analysis = .env$analysis_label, exposure = .env$exposure,
    n_pairs = nrow(data), n_residents = n_distinct(data$Participant_id),
    n_exposed = n_exposed_value,
    outcome_events = sum(data$any_gene_family_difference %in% TRUE, na.rm = TRUE),
    model = paste0("any CGE gene-family difference ~ exposure + natural spline(days, df=", SPLINE_DF, ")"),
    .before = 1
  )
}

model_continuous <- function(data, outcome) {
  z <- data %>% mutate(.exposure = as.integer(TotalSNPs <= PRIMARY_SNP),
                       .outcome = as.numeric(.data[[outcome]])) %>%
    filter(is.finite(.outcome), is.finite(days_between))
  if (nrow(z) < 20L || n_distinct(z$.exposure) != 2L) return(NA_real_)
  formula <- as.formula(paste0(".outcome ~ .exposure + splines::ns(days_between, df=", SPLINE_DF, ")"))
  fit <- lm(formula, data = z)
  nd1 <- z; nd0 <- z; nd1$.exposure <- 1L; nd0$.exposure <- 0L
  mean(predict(fit, nd1) - predict(fit, nd0), na.rm = TRUE)
}

continuous_inference <- function(data, outcome, profile, seed_offset) {
  point <- model_continuous(data, outcome)
  draws <- bootstrap_apply(data, function(z) model_continuous(z, outcome), BOOT_REPS,
                           RQ_SEED + seed_offset, BOOT_WORKERS)
  bootstrap_interval(draws[1, ], point) %>% mutate(
    profile = profile, outcome = outcome,
    estimand = "adjusted_mean_difference_SNP_close_minus_not_close",
    n_pairs = nrow(data), n_residents = n_distinct(data$Participant_id),
    model = paste0(outcome, " ~ SNP<=", PRIMARY_SNP, " + natural spline(days, df=", SPLINE_DF, ")"),
    .before = 1
  )
}

# -----------------------------------------------------------------------------
# RQ06: adjacent-pair CGE gene-family stability
# -----------------------------------------------------------------------------
rq06_inference <- list(); rq06_continuous <- list(); rq06_pairs <- list(); rq06_gene_changes <- list()
for (profile_index in seq_along(expected_profiles)) {
  profile <- expected_profiles[[profile_index]]
  message("RQ06 ", profile, ": pair metrics and resident bootstraps")
  z <- pair_metrics(adjacent_base, "key_from", "key_to", profile) %>%
    mutate(close_primary = TotalSNPs <= PRIMARY_SNP,
           routine_routine = event_type_from == "Routine" & event_type_to == "Routine")
  rq06_inference[[length(rq06_inference) + 1L]] <- binary_inference(
    z, "close_primary", paste0("primary_SNP_le_", PRIMARY_SNP), profile, 600L + 1000L * profile_index)
  for (threshold in SNP_THRESHOLDS) {
    threshold_name <- paste0("close_", threshold)
    z[[threshold_name]] <- z$TotalSNPs <= threshold
    rq06_inference[[length(rq06_inference) + 1L]] <- binary_inference(
      z, threshold_name, paste0("threshold_SNP_le_", threshold), profile,
      610L + threshold + 1000L * profile_index)
  }
  rq06_inference[[length(rq06_inference) + 1L]] <- binary_inference(
    filter(z, routine_routine), "close_primary", "routine_to_routine_only", profile,
    710L + 1000L * profile_index)
  rq06_continuous[[length(rq06_continuous) + 1L]] <- continuous_inference(
    z, "cge_jaccard", profile, 720L + 1000L * profile_index)
  rq06_continuous[[length(rq06_continuous) + 1L]] <- continuous_inference(
    z, "absolute_gene_family_count_change", profile, 730L + 1000L * profile_index)

  deidentified <- z %>%
    left_join(resident_key, by = "Participant_id", relationship = "many-to-one") %>%
    left_join(rename(episode_key_map, key_from = episode_key, episode_from_label = episode_label), by = "key_from") %>%
    left_join(rename(episode_key_map, key_to = episode_key, episode_to_label = episode_label), by = "key_to") %>%
    transmute(profile, resident_label, episode_from_label, episode_to_label,
              date_from, date_to, days_between, status_from, status_to,
              event_type_from, event_type_to, TotalSNPs, close_primary,
              any_gene_family_difference, cge_jaccard,
              absolute_gene_family_count_change, n_families_from_only,
              n_families_to_only, families_from_only, families_to_only)
  rq06_pairs[[length(rq06_pairs) + 1L]] <- deidentified
  changes <- bind_rows(
    z %>% transmute(direction = "from_only", context = if_else(close_primary, "SNP_close", "SNP_not_close"),
                    gene_family = strsplit(families_from_only, ";", fixed = TRUE)),
    z %>% transmute(direction = "to_only", context = if_else(close_primary, "SNP_close", "SNP_not_close"),
                    gene_family = strsplit(families_to_only, ";", fixed = TRUE))
  ) %>% unnest_longer(gene_family) %>% filter(nzchar(gene_family)) %>%
    count(direction, context, gene_family, name = "n_pairs") %>%
    mutate(profile = profile, note = "Descriptive only; no gene-level hypothesis test.", .before = 1)
  rq06_gene_changes[[length(rq06_gene_changes) + 1L]] <- changes
}
rq06_inference <- bind_rows(rq06_inference)
rq06_continuous <- bind_rows(rq06_continuous)
rq06_pairs <- bind_rows(rq06_pairs)
rq06_gene_changes <- bind_rows(rq06_gene_changes)
atomic_write_csv(rq06_inference, file.path(out_rq06, "rq06_resident_bootstrap_inference.csv"))
atomic_write_csv(rq06_continuous, file.path(out_rq06, "rq06_secondary_adjusted_contrasts.csv"))
atomic_write_csv(rq06_pairs, file.path(out_rq06, "rq06_adjacent_pair_changes_deidentified.csv"))
atomic_write_csv(rq06_gene_changes, file.path(out_rq06, "rq06_gene_family_changes_descriptive.csv"))

rq06_descriptive <- rq06_pairs %>% mutate(snp_context = if_else(close_primary, "SNP_close", "SNP_not_close")) %>%
  group_by(profile, snp_context) %>% summarise(
    n_pairs = n(), n_residents = n_distinct(resident_label),
    any_difference_n = sum(any_gene_family_difference),
    any_difference_pct = 100 * mean(any_gene_family_difference),
    median_jaccard = median(cge_jaccard),
    median_absolute_count_change = median(absolute_gene_family_count_change), .groups = "drop")
atomic_write_csv(rq06_descriptive, file.path(out_rq06, "rq06_descriptive_summary.csv"))

# -----------------------------------------------------------------------------
# RQ07: event-sample burden analogue and explicitly limited ExPEC proxy
# -----------------------------------------------------------------------------
event_primary <- primary_episode %>% filter(Event_type == "UTI_event")
event_anchor <- c(nrow(event_primary), n_distinct(event_primary$Participant_id),
                  sum(event_primary$UTI_Status == "UTI"), sum(event_primary$UTI_Status == "Not_UTI"))
expected_event <- c(as.integer(cohort$expected_event_samples), as.integer(cohort$expected_event_residents),
                    as.integer(cohort$expected_event_uti), as.integer(cohort$expected_event_not_uti))
check("rq07_event_anchor", identical(as.integer(event_anchor), expected_event),
      paste(expected_event, collapse = ";"), paste(event_anchor, collapse = ";"),
      "Exact approved 32-sample UTI-event anchor")

median_difference <- function(data, endpoint) {
  median(data[[endpoint]][data$UTI_Status == "UTI"], na.rm = TRUE) -
    median(data[[endpoint]][data$UTI_Status == "Not_UTI"], na.rm = TRUE)
}

median_inference <- function(data, endpoint, analysis_label, profile, seed_offset) {
  point <- median_difference(data, endpoint)
  draws <- bootstrap_apply(data, function(z) median_difference(z, endpoint), BOOT_REPS,
                           RQ_SEED + seed_offset, BOOT_WORKERS)
  bootstrap_interval(draws[1, ], point) %>% mutate(
    profile = profile, analysis = analysis_label, endpoint = endpoint,
    estimand = "median_UTI_minus_median_Not_UTI", n_rows = nrow(data),
    n_residents = n_distinct(data$Participant_id), n_uti = sum(data$UTI_Status == "UTI"),
    n_not_uti = sum(data$UTI_Status == "Not_UTI"),
    inference_note = "Robustness estimate with resident-bootstrap CI; no new confirmatory p-value.", .before = 1)
}

paired_inference <- function(data, endpoint, profile, seed_offset) {
  delta <- paste0(endpoint, "_delta")
  statistic <- function(z) median(z[[delta]], na.rm = TRUE)
  point <- statistic(data)
  draws <- bootstrap_apply(data, statistic, BOOT_REPS, RQ_SEED + seed_offset, BOOT_WORKERS)
  bootstrap_interval(draws[1, ], point) %>% mutate(
    profile = profile, analysis = "nearest_within_resident_pair", endpoint = endpoint,
    estimand = "median_within_resident_UTI_minus_Not_UTI", n_pairs = nrow(data),
    n_positive = sum(data[[delta]] > 0), n_zero = sum(data[[delta]] == 0),
    n_negative = sum(data[[delta]] < 0),
    inference_note = "Robustness estimate with resident-bootstrap CI; no new confirmatory p-value.", .before = 1)
}

rq07_results <- list(); rq07_events <- list(); rq07_nearest <- list(); proxy_agreement <- list()
for (profile_index in seq_along(expected_profiles)) {
  profile <- expected_profiles[[profile_index]]
  pm <- metrics %>% filter(.data$profile == .env$profile) %>%
    select(episode_key, distinct_cge_gene_family_count, CGE_available_ExPEC_group_count,
           starts_with("cge_"), stx1_present, stx2_present)
  full <- primary_episode %>% left_join(pm, by = "episode_key", relationship = "one-to-one")
  check(paste0("rq07_metric_join_", profile), nrow(full) == as.integer(cohort$expected_assemblies) &&
          !anyNA(full$distinct_cge_gene_family_count), "532 complete", sum(!is.na(full$distinct_cge_gene_family_count)),
        "CGE metrics must join one-to-one to approved episodes")
  event <- full %>% filter(Event_type == "UTI_event") %>% arrange(Participant_id, collection_date, tp_lab)
  endpoints <- c("distinct_cge_gene_family_count", "CGE_available_ExPEC_group_count")
  for (endpoint_index in seq_along(endpoints)) {
    endpoint <- endpoints[[endpoint_index]]
    rq07_results[[length(rq07_results) + 1L]] <- median_inference(
      event, endpoint, "all_UTI_event_samples", profile,
      800L + 10L * endpoint_index + 1000L * profile_index)
  }
  earliest <- event %>% group_by(Participant_id) %>% arrange(collection_date, tp_lab, .by_group = TRUE) %>%
    slice_head(n = 1L) %>% ungroup()
  check(paste0("rq07_earliest_", profile), nrow(earliest) == as.integer(cohort$expected_event_residents),
        cohort$expected_event_residents, nrow(earliest), "One earliest event sample per event resident")
  for (endpoint_index in seq_along(endpoints)) {
    endpoint <- endpoints[[endpoint_index]]
    rq07_results[[length(rq07_results) + 1L]] <- median_inference(
      earliest, endpoint, "earliest_event_per_resident", profile,
      820L + 10L * endpoint_index + 1000L * profile_index)
  }
  paired_ids <- full %>% group_by(Participant_id) %>%
    summarise(n_status = n_distinct(UTI_Status), .groups = "drop") %>%
    filter(n_status == 2L) %>% pull(Participant_id)
  check(paste0("rq07_paired_residents_", profile), length(paired_ids) == as.integer(cohort$expected_paired_residents),
        cohort$expected_paired_residents, length(paired_ids), "Fixed nearest UTI/Not_UTI sensitivity")
  nearest <- bind_rows(lapply(paired_ids, function(pid) {
    z <- filter(full, Participant_id == pid)
    uti <- filter(z, UTI_Status == "UTI")
    not <- filter(z, UTI_Status == "Not_UTI")
    candidates <- crossing(iu = seq_len(nrow(uti)), inot = seq_len(nrow(not))) %>%
      mutate(days = abs(as.numeric(uti$collection_date[iu] - not$collection_date[inot])),
             uti_key = uti$episode_key[iu], not_key = not$episode_key[inot]) %>%
      arrange(days, uti_key, not_key) %>% slice_head(n = 1L)
    tibble(
      Participant_id = pid, uti_key = candidates$uti_key, not_key = candidates$not_key,
      absolute_days = candidates$days,
      distinct_cge_gene_family_count_delta =
        uti$distinct_cge_gene_family_count[candidates$iu] - not$distinct_cge_gene_family_count[candidates$inot],
      CGE_available_ExPEC_group_count_delta =
        uti$CGE_available_ExPEC_group_count[candidates$iu] - not$CGE_available_ExPEC_group_count[candidates$inot]
    )
  }))
  for (endpoint_index in seq_along(endpoints)) {
    endpoint <- endpoints[[endpoint_index]]
    rq07_results[[length(rq07_results) + 1L]] <- paired_inference(
      nearest, endpoint, profile, 840L + 10L * endpoint_index + 1000L * profile_index)
  }
  rq07_nearest[[length(rq07_nearest) + 1L]] <- nearest %>%
    left_join(resident_key, by = "Participant_id") %>%
    left_join(rename(episode_key_map, uti_key = episode_key, uti_episode_label = episode_label), by = "uti_key") %>%
    left_join(rename(episode_key_map, not_key = episode_key, not_uti_episode_label = episode_label), by = "not_key") %>%
    transmute(profile, resident_label, uti_episode_label, not_uti_episode_label, absolute_days,
              distinct_cge_gene_family_count_delta, CGE_available_ExPEC_group_count_delta)
  rq07_events[[length(rq07_events) + 1L]] <- event %>%
    left_join(resident_key, by = "Participant_id") %>% left_join(episode_key_map, by = "episode_key") %>%
    transmute(profile, resident_label, episode_label, collection_date, UTI_Status,
              distinct_cge_gene_family_count, CGE_available_ExPEC_group_count,
              cge_pap_group, cge_sfa_foc_group, cge_afa_dra_group,
              cge_aerobactin_group, cge_capsule_group, stx1_present, stx2_present)
  for (scope_name in c("all_532", "event_32")) {
    scope_data <- if (scope_name == "all_532") full else event
    difference <- scope_data$CGE_available_ExPEC_group_count - scope_data$five_marker_expec_count
    proxy_agreement[[length(proxy_agreement) + 1L]] <- tibble(
      profile = profile, scope = scope_name, n = nrow(scope_data),
      exact_count_agreement_n = sum(difference == 0),
      exact_count_agreement_pct = 100 * mean(difference == 0),
      mean_absolute_difference = mean(abs(difference)), median_absolute_difference = median(abs(difference)),
      spearman_correlation = suppressWarnings(cor(scope_data$CGE_available_ExPEC_group_count,
                                                   scope_data$five_marker_expec_count,
                                                   method = "spearman", use = "complete.obs")),
      note = "CGE proxy is not the primary five-marker ExPEC count and is not numerically interchangeable."
    )
  }
}
rq07_results <- bind_rows(rq07_results)
atomic_write_csv(rq07_results, file.path(out_rq07, "rq07_resident_bootstrap_robustness_estimates.csv"))
atomic_write_csv(bind_rows(rq07_events), file.path(out_rq07, "rq07_event_samples_deidentified.csv"))
atomic_write_csv(bind_rows(rq07_nearest), file.path(out_rq07, "rq07_nearest_pairs_deidentified.csv"))
atomic_write_csv(bind_rows(proxy_agreement), file.path(out_rq07, "rq07_cge_proxy_vs_primary_five_marker.csv"))

# -----------------------------------------------------------------------------
# RQ08: CGE Jaccard as a continuous relatedness proxy and LOOR composite
# -----------------------------------------------------------------------------
auc_inference <- function(data, score, label, profile, seed_offset) {
  z <- data %>% mutate(close = TotalSNPs <= PRIMARY_SNP, .score = as.numeric(.data[[score]])) %>%
    filter(is.finite(.score))
  point <- auc_rank(z$close, z$.score)
  draws <- bootstrap_apply(z, function(b) auc_rank(b$close, b$.score), BOOT_REPS,
                           RQ_SEED + seed_offset, BOOT_WORKERS)
  bootstrap_interval(draws[1, ], point) %>% mutate(
    profile = profile, predictor = label, score_column = score,
    direction = "higher_is_close", n_pairs = nrow(z), n_residents = n_distinct(z$Participant_id),
    reference = paste0("direct SNP <=", PRIMARY_SNP), .before = 1)
}

fit_loor <- function(data, profile) {
  z <- data %>% mutate(close_reference = as.integer(TotalSNPs <= PRIMARY_SNP),
                       same_st_num = as.integer(same_provider_ST)) %>%
    filter(as_bool(provider_ST_both), is.finite(MashDistance), is.finite(cge_jaccard))
  ids <- unique(z$Participant_id)
  bind_rows(lapply(ids, function(held_out) {
    train <- filter(z, Participant_id != held_out)
    test <- filter(z, Participant_id == held_out)
    mash_mean <- mean(train$MashDistance); mash_sd <- sd(train$MashDistance)
    cge_mean <- mean(train$cge_jaccard); cge_sd <- sd(train$cge_jaccard)
    if (!is.finite(mash_sd) || mash_sd == 0) mash_sd <- 1
    if (!is.finite(cge_sd) || cge_sd == 0) cge_sd <- 1
    train <- train %>% mutate(mash_z = (MashDistance - mash_mean) / mash_sd,
                              cge_z = (cge_jaccard - cge_mean) / cge_sd)
    test <- test %>% mutate(mash_z = (MashDistance - mash_mean) / mash_sd,
                            cge_z = (cge_jaccard - cge_mean) / cge_sd)
    fit <- tryCatch(suppressWarnings(glm(close_reference ~ same_st_num + mash_z + cge_z,
                                         family = binomial(), data = train,
                                         control = glm.control(maxit = 100))), error = function(error) NULL)
    prediction <- if (is.null(fit)) rep(NA_real_, nrow(test)) else
      suppressWarnings(as.numeric(predict(fit, newdata = test, type = "response")))
    test %>% transmute(Participant_id, pair_id, profile = profile,
                       outcome = close_reference, predicted_probability = prediction)
  }))
}

rq08_auc <- list(); rq08_composite <- list(); rq08_predictions <- list()
for (profile_index in seq_along(expected_profiles)) {
  profile <- expected_profiles[[profile_index]]
  message("RQ08 ", profile, ": direct-pair Jaccard and LOOR model")
  z <- pair_metrics(direct, "key_A", "key_B", profile)
  rq08_auc[[length(rq08_auc) + 1L]] <- auc_inference(
    z, "cge_jaccard", "CGE gene-family Jaccard", profile, 1000L + 1000L * profile_index)
  predictions <- fit_loor(z, profile)
  valid <- filter(predictions, is.finite(predicted_probability))
  point <- auc_rank(valid$outcome, valid$predicted_probability)
  draws <- bootstrap_apply(valid, function(b) auc_rank(b$outcome, b$predicted_probability), BOOT_REPS,
                           RQ_SEED + 1100L + 1000L * profile_index, BOOT_WORKERS)
  rq08_composite[[length(rq08_composite) + 1L]] <- bootstrap_interval(draws[1, ], point) %>%
    mutate(profile = profile, scope = "all_within_resident_pairs", reference_threshold = PRIMARY_SNP,
           n_pairs = nrow(valid), n_residents = n_distinct(valid$Participant_id),
           brier_score = mean((valid$predicted_probability - valid$outcome)^2),
           sensitivity_at_0_5 = safe_ratio(sum(valid$predicted_probability >= 0.5 & valid$outcome == 1),
                                           sum(valid$outcome == 1)),
           specificity_at_0_5 = safe_ratio(sum(valid$predicted_probability < 0.5 & valid$outcome == 0),
                                           sum(valid$outcome == 0)),
           model = "LOOR logistic: provider ST + scaled Mash distance + scaled CGE Jaccard", .before = 1) %>%
    rename(auc = estimate, auc_ci_lower = ci_lower, auc_ci_upper = ci_upper)
  rq08_predictions[[length(rq08_predictions) + 1L]] <- valid %>%
    left_join(resident_key, by = "Participant_id") %>%
    transmute(profile, resident_label, pair_label = sprintf("PAIR%04d", match(pair_id, sort(unique(direct$pair_id)))),
              outcome, predicted_probability)
}
rq08_auc <- bind_rows(rq08_auc)
rq08_composite <- bind_rows(rq08_composite)
atomic_write_csv(rq08_auc, file.path(out_rq08, "rq08_cge_jaccard_auc_resident_bootstrap.csv"))
atomic_write_csv(rq08_composite, file.path(out_rq08, "rq08_loor_composite_summary.csv"))
atomic_write_csv(bind_rows(rq08_predictions), file.path(out_rq08, "rq08_loor_predictions_deidentified.csv"))

primary_auc <- read_csv(file.path(primary_dir, "RQ08", "rq08_roc_auc_resident_bootstrap.csv"), show_col_types = FALSE)
primary_composite <- read_csv(file.path(primary_dir, "RQ08", "rq08_loor_composite_summary.csv"), show_col_types = FALSE) %>%
  filter(scope == "all_within_resident_pairs", reference_threshold == PRIMARY_SNP)
primary_vf_auc <- primary_auc %>% filter(predictor == "Fresh VF Jaccard") %>% slice_head(n = 1)
rq08_comparison <- bind_rows(
  rq08_auc %>% transmute(profile, comparison = "continuous_jaccard_auc",
                         cge_estimate = estimate, cge_ci_lower = ci_lower, cge_ci_upper = ci_upper,
                         primary_vfdb_estimate = primary_vf_auc$estimate,
                         cge_minus_primary = estimate - primary_vf_auc$estimate),
  rq08_composite %>% transmute(profile, comparison = "loor_composite_auc",
                               cge_estimate = auc, cge_ci_lower = auc_ci_lower, cge_ci_upper = auc_ci_upper,
                               primary_vfdb_estimate = primary_composite$auc,
                               cge_minus_primary = auc - primary_composite$auc)
) %>% mutate(note = "Provider-ST diagnostic sensitivity, specificity, PPV, NPV and balanced accuracy are unchanged.")
atomic_write_csv(rq08_comparison, file.path(out_rq08, "rq08_comparison_with_primary_vfdb.csv"))

# -----------------------------------------------------------------------------
# Cross-RQ interpretation and final release gates
# -----------------------------------------------------------------------------
primary_rq06_bin <- read_csv(file.path(primary_dir, "RQ06", "rq06_resident_bootstrap_inference.csv"), show_col_types = FALSE) %>%
  filter(analysis == paste0("primary_SNP_le_", PRIMARY_SNP), estimand == "odds_ratio") %>% slice_head(n = 1)
primary_rq06_jac <- read_csv(file.path(primary_dir, "RQ06", "rq06_secondary_adjusted_contrasts.csv"), show_col_types = FALSE) %>%
  filter(outcome == "vf_jaccard") %>% slice_head(n = 1)
rq06_comparison <- bind_rows(
  rq06_inference %>% filter(analysis == paste0("primary_SNP_le_", PRIMARY_SNP), estimand == "odds_ratio") %>%
    transmute(profile, endpoint = "any_profile_difference_odds_ratio",
              cge_estimate = estimate, cge_ci_lower = ci_lower, cge_ci_upper = ci_upper,
              primary_vfdb_estimate = primary_rq06_bin$estimate,
              same_direction = (estimate < 1) == (primary_rq06_bin$estimate < 1)),
  rq06_continuous %>% filter(outcome == "cge_jaccard") %>%
    transmute(profile, endpoint = "adjusted_jaccard_difference",
              cge_estimate = estimate, cge_ci_lower = ci_lower, cge_ci_upper = ci_upper,
              primary_vfdb_estimate = primary_rq06_jac$estimate,
              same_direction = sign(estimate) == sign(primary_rq06_jac$estimate))
)
atomic_write_csv(rq06_comparison, file.path(out_rq06, "rq06_comparison_with_primary_vfdb.csv"))

primary_rq07 <- read_csv(file.path(primary_dir, "RQ07", "rq07_event_sample_inference.csv"), show_col_types = FALSE)
rq07_comparison <- rq07_results %>% filter(analysis == "all_UTI_event_samples") %>%
  mutate(primary_endpoint = if_else(endpoint == "distinct_cge_gene_family_count",
                                    "total_vf_count_curated", "five_marker_expec_count")) %>%
  left_join(primary_rq07 %>% select(primary_endpoint = endpoint,
                                    primary_estimate = estimate, primary_ci_lower = ci_lower,
                                    primary_ci_upper = ci_upper), by = "primary_endpoint") %>%
  transmute(profile, cge_endpoint = endpoint, primary_endpoint,
            cge_estimate = estimate, cge_ci_lower = ci_lower, cge_ci_upper = ci_upper,
            primary_estimate, primary_ci_lower, primary_ci_upper,
            same_direction = sign(cge_estimate) == sign(primary_estimate),
            note = "CGE total-family burden and limited ExPEC proxy are database-specific robustness endpoints, not replacements.")
atomic_write_csv(rq07_comparison, file.path(out_rq07, "rq07_comparison_with_primary_vfdb.csv"))

interpretation <- bind_rows(
  rq06_comparison %>% group_by(profile) %>% summarise(
    research_question = "RQ06", direction_agreement = all(same_direction),
    summary = if_else(direction_agreement,
                      "CGE stability results agree in direction with the primary VFDB analysis.",
                      "At least one CGE stability contrast differs in direction from the primary VFDB analysis."), .groups = "drop"),
  rq07_comparison %>% group_by(profile) %>% summarise(
    research_question = "RQ07", direction_agreement = all(same_direction),
    summary = if_else(direction_agreement,
                      "CGE event-sample robustness estimates agree in direction with the corresponding primary VFDB endpoints.",
                      "At least one CGE event-sample robustness estimate differs in direction; database definitions limit direct comparison."),
    .groups = "drop"),
  rq08_comparison %>% group_by(profile) %>% summarise(
    research_question = "RQ08", direction_agreement = all(cge_estimate > 0.5),
    summary = if_else(direction_agreement,
                      "CGE Jaccard and the CGE composite both retain out-of-resident discrimination above chance.",
                      "At least one CGE relatedness proxy does not retain discrimination above chance."), .groups = "drop")
) %>% mutate(
  interpretation_impact = if_else(direction_agreement, "primary interpretation not made less robust",
                                  "primary interpretation requires caution in this sensitivity layer"),
  analysis_role = "CGE VirulenceFinder sensitivity analysis; primary ABRicate/VFDB results unchanged"
)
atomic_write_csv(interpretation, file.path(out_rq, "interpretation_summary.csv"))

primary_after <- tree_fingerprint(primary_dir)
check("primary_tree_unchanged", identical(primary_before, primary_after), primary_before, primary_after,
      "No file inside the active primary research-question output tree may change")
check("rq06_final_denominator", nrow(rq06_pairs) == 2L * as.integer(cohort$expected_adjacent_pairs),
      2L * as.integer(cohort$expected_adjacent_pairs), nrow(rq06_pairs), "371 pairs per profile")
check("rq07_final_denominator", nrow(bind_rows(rq07_events)) == 2L * as.integer(cohort$expected_event_samples),
      2L * as.integer(cohort$expected_event_samples), nrow(bind_rows(rq07_events)), "32 event samples per profile")
check("rq08_final_denominator", all(rq08_auc$n_pairs == as.integer(cohort$expected_direct_pairs)),
      paste(cohort$expected_direct_pairs, "per profile"), paste(rq08_auc$n_pairs, collapse = ";"), "893 direct pairs per profile")
check("bootstrap_outputs", all(rq06_inference$bootstrap_reps_requested == BOOT_REPS) &&
        all(rq06_continuous$bootstrap_reps_requested == BOOT_REPS) &&
        all(rq07_results$bootstrap_reps_requested == BOOT_REPS) &&
        all(rq08_auc$bootstrap_reps_requested == BOOT_REPS) &&
        all(rq08_composite$bootstrap_reps_requested == BOOT_REPS),
      BOOT_REPS, "all final tables", "All uncertainty uses the prespecified resident bootstrap")

atomic_write_csv(checks, file.path(out_rq, "acceptance_checks.csv"))
check_hashes <- tibble(
  file = c("config", "R_analysis_script", "run_manifest", "hits_long", "episode_metrics", "rq06_inference",
           "rq07_inference", "rq08_auc", "primary_release_tree"),
  sha256 = c(sha256_file(file.path(root, "config", "virulencefinder_sensitivity.toml")),
             sha256_file(file.path(root, "scripts", "analyze_virulencefinder_sensitivity.R")),
             sha256_file(file.path(output, "run_manifest.csv")),
             sha256_file(file.path(output, "hits_long.csv")),
             sha256_file(file.path(output, "episode_metrics.csv")),
             sha256_file(file.path(out_rq06, "rq06_resident_bootstrap_inference.csv")),
             sha256_file(file.path(out_rq07, "rq07_resident_bootstrap_robustness_estimates.csv")),
             sha256_file(file.path(out_rq08, "rq08_cge_jaccard_auc_resident_bootstrap.csv")),
             primary_after)
)
atomic_write_csv(check_hashes, file.path(out_rq, "release_hashes.csv"))

report_lines <- c(
  "# CGE VirulenceFinder sensitivity analysis: final report",
  "",
  paste0("Generated: ", timestamp_utc()),
  "",
  "This is a separate CGE VirulenceFinder sensitivity analysis for RQ06-RQ08. The primary ABRicate/VFDB 80/80 release remains unchanged and authoritative.",
  "",
  "## Scientific contract",
  "",
  paste0("- 532 QC-passed canonical Longcycler assemblies from 161 residents."),
  paste0("- VirulenceFinder 3.2.1; database 2.1.0 commit ", cfg$software$database_commit, "; BLAST 2.17.0+."),
  "- Profiles: web-default identity 90%/coverage 60%, and threshold-matched identity 80%/coverage 80%.",
  "- Presence/absence is defined over 681 unique gene families (680 virulence_ecoli plus stx1/stx2, with stx2 shared across the selected databases).",
  paste0("- Resident bootstrap: ", BOOT_REPS, " replicates; seed ", RQ_SEED, "."),
  "",
  "## Interpretation",
  "",
  paste0("- ", interpretation$profile, " / ", interpretation$research_question, ": ",
         interpretation$summary, " Impact: ", interpretation$interpretation_impact, "."),
  "",
  "Raw total burdens across CGE and VFDB are not interpreted as acquisition or loss because the feature universes differ. The CGE_available_ExPEC_group_count is explicitly a limited proxy, not the primary five-marker ExPEC count. No new confirmatory RQ07 p-values were introduced.",
  "",
  "All acceptance checks passed; the primary research-question tree fingerprint was identical before and after this analysis."
)
atomic_write_lines(report_lines, file.path(out_rq, "FINAL_REPORT.md"))

atomic_write_csv(tibble(
  state = "COMPLETE", generated_at_utc = timestamp_utc(),
  assemblies = as.integer(cohort$expected_assemblies), jobs = expected_jobs,
  profiles = length(expected_profiles), bootstrap_reps = BOOT_REPS,
  seed = RQ_SEED, primary_release_modified = FALSE
), file.path(out_rq, "RUN_STATUS.csv"))

atomic_write_lines(c(
  "VirulenceFinder sensitivity release: PASS",
  paste0("Generated: ", timestamp_utc()),
  paste0("Assemblies: ", cohort$expected_assemblies),
  paste0("Profiles: ", length(expected_profiles)),
  paste0("Jobs: ", expected_jobs),
  paste0("Gene-family universe: ", cohort$expected_combined_gene_families),
  paste0("RQ06 adjacent pairs/profile: ", cohort$expected_adjacent_pairs),
  paste0("RQ07 event samples/profile: ", cohort$expected_event_samples),
  paste0("RQ08 direct pairs/profile: ", cohort$expected_direct_pairs),
  paste0("Bootstrap replicates: ", BOOT_REPS),
  paste0("Seed: ", RQ_SEED),
  "Primary ABRicate/VFDB release modified: NO",
  "All acceptance checks: PASS"
), final_marker)

message("CGE VirulenceFinder RQ06-RQ08 sensitivity release: PASS")

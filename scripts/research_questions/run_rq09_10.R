#!/usr/bin/env Rscript

# RQ09-RQ10 release analyses for the exact selected Longcycler cohort.

options(stringsAsFactors = FALSE, warn = 1)

suppressPackageStartupMessages({
  library(digest)
  library(dplyr)
  library(jsonlite)
  library(purrr)
  library(readr)
  library(stringr)
  library(tidyr)
})

EXPECTED_EPISODES <- 532L
EXPECTED_RESIDENTS <- 161L
EXPECTED_UTI <- 16L
EXPECTED_NOT_UTI <- 516L
EXPECTED_ALL_PAIRS <- 893L
EXPECTED_ADJACENT_PAIRS <- 371L
EXPECTED_ADJACENT_RESIDENTS <- 139L
EXPECTED_ADJACENT_LE25 <- 140L
EXPECTED_NOT_UTI_TO_UTI <- 9L
EXPECTED_NOT_UTI_TO_UTI_LE25 <- 5L

project_root <- normalizePath(
  if (file.exists("00_config.R")) "." else file.path("..", ".."),
  winslash = "/", mustWork = TRUE
)
setwd(project_root)

out_root <- file.path("results", "research_questions")
input_root <- file.path(out_root, "_inputs")
rq_dirs <- setNames(file.path(out_root, c("RQ09", "RQ10")), c("RQ09", "RQ10"))
dir.create(input_root, recursive = TRUE, showWarnings = FALSE)
walk(rq_dirs, dir.create, recursive = TRUE, showWarnings = FALSE)

run_timestamp <- format(Sys.time(), tz = "UTC", usetz = TRUE)
workers <- suppressWarnings(as.integer(Sys.getenv("RQ_WORKERS", "8")))
n_perm <- suppressWarnings(as.integer(Sys.getenv("RQ_PERMUTATIONS", "10000")))
n_boot <- suppressWarnings(as.integer(Sys.getenv("RQ_BOOTSTRAP_REPS", "10000")))
seed <- suppressWarnings(as.integer(Sys.getenv("RQ_SEED", "20260712")))
if (is.na(workers) || workers < 1L) workers <- 1L
if (is.na(n_perm) || n_perm < 100L) stop("RQ_PERMUTATIONS must be at least 100.")
if (is.na(n_boot) || n_boot < 100L) stop("RQ_BOOTSTRAP_REPS must be at least 100.")
if (is.na(seed)) seed <- 20260712L
allow_reduced_reps <- tolower(Sys.getenv(
  "RQ_ALLOW_REDUCED_REPS", "0"
)) %in% c("1", "true", "yes")
if (
  !allow_reduced_reps &&
    (n_perm != 10000L || n_boot != 10000L)
) {
  stop(
    "Final RQ09/RQ10 release requires exactly 10,000 permutations and ",
    "10,000 bootstrap replicates. Set RQ_ALLOW_REDUCED_REPS=1 only for ",
    "non-release smoke tests.",
    call. = FALSE
  )
}

forbidden_release_token <- paste0("fl", "ye")

assert_release_clean <- function(x, context) {
  values <- as.character(unlist(x, recursive = TRUE, use.names = TRUE))
  values <- values[!is.na(values)]
  if (any(str_detect(str_to_lower(values), fixed(forbidden_release_token)))) {
    stop("Refusing to publish excluded-assembler content in ", context, call. = FALSE)
  }
  invisible(TRUE)
}

atomic_write_csv <- function(x, path) {
  assert_release_clean(c(names(x), unlist(x, recursive = TRUE, use.names = FALSE)), path)
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  tmp <- tempfile(pattern = paste0(".", basename(path), "."), tmpdir = dirname(path))
  readr::write_csv(x, tmp, na = "")
  if (!file.rename(tmp, path)) {
    file.copy(tmp, path, overwrite = TRUE)
    unlink(tmp)
  }
  invisible(path)
}

atomic_write_lines <- function(x, path) {
  assert_release_clean(x, path)
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  tmp <- tempfile(pattern = paste0(".", basename(path), "."), tmpdir = dirname(path))
  writeLines(x, tmp, useBytes = TRUE)
  if (!file.rename(tmp, path)) {
    file.copy(tmp, path, overwrite = TRUE)
    unlink(tmp)
  }
  invisible(path)
}

atomic_write_json <- function(x, path) {
  assert_release_clean(x, path)
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  tmp <- tempfile(pattern = paste0(".", basename(path), "."), tmpdir = dirname(path))
  jsonlite::write_json(x, tmp, auto_unbox = TRUE, pretty = TRUE, null = "null", na = "null")
  if (!file.rename(tmp, path)) {
    file.copy(tmp, path, overwrite = TRUE)
    unlink(tmp)
  }
  invisible(path)
}

write_status <- function(rq, status, reason, detail = NA_character_) {
  atomic_write_csv(tibble(
    research_question = rq,
    status = status,
    reason = reason,
    detail = detail,
    run_timestamp_utc = run_timestamp
  ), file.path(rq_dirs[[rq]], "analysis_status.csv"))
}

require_columns <- function(x, cols, context) {
  missing <- setdiff(cols, names(x))
  if (length(missing)) stop(context, " lacks: ", paste(missing, collapse = ", "), call. = FALSE)
}

as_flag <- function(x) {
  if (is.logical(x)) return(replace(x, is.na(x), FALSE))
  str_to_lower(str_trim(as.character(x))) %in% c("true", "t", "1", "yes", "y")
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
  for (fmt in c("%Y-%m-%d", "%d/%m/%Y", "%d-%m-%Y", "%Y/%m/%d")) {
    take <- is.na(out) & !is.na(x) & nzchar(str_trim(x))
    out[take] <- as.Date(x[take], format = fmt)
  }
  out
}

sha256_file <- function(path) unname(digest::digest(path, algo = "sha256", file = TRUE))

parallel_lapply <- function(x, fun, mc.cores = workers) {
  if (.Platform$OS.type != "unix" || mc.cores <= 1L) return(lapply(x, fun))
  parallel::mclapply(x, fun, mc.cores = mc.cores, mc.preschedule = TRUE)
}

tool_version <- function(tool, args = "--version") {
  exe <- Sys.which(tool)
  if (exe == "") return(NA_character_)
  paste(suppressWarnings(system2(exe, args, stdout = TRUE, stderr = TRUE)), collapse = " ")
}

run_system2 <- function(command, args, stdout, stderr) {
  status <- suppressWarnings(system2(command, args, stdout = stdout, stderr = stderr))
  if (is.null(status)) 0L else as.integer(status)
}

required_inputs <- c(
  cohort = file.path("results", "clinical", "analysis_cohort_longcycler.csv"),
  manifest = file.path("results", "qc", "analysis_assembly_manifest.csv"),
  provider_mlst = file.path("results", "mlst", "mlst_provider_preferred.csv"),
  transitions = file.path("results", "longitudinal", "longcycler_transitions.csv"),
  direct_pairs = file.path(input_root, "direct_pair_metrics_893.csv")
)
missing_inputs <- required_inputs[!file.exists(required_inputs)]
if (length(missing_inputs)) stop("Missing required input(s): ", paste(missing_inputs, collapse = ", "))

cohort <- read_csv(required_inputs[["cohort"]], show_col_types = FALSE) %>%
  mutate(
    Participant_id = as.character(.data$Participant_id),
    tp_lab = as.character(.data$tp_lab),
    Assembly_ID = as.character(.data$Assembly_ID),
    Isolate_ID = as.character(.data$Isolate_ID),
    episode_key = episode_key(.data$Participant_id, .data$tp_lab),
    full_path = norm_path(.data$full_path),
    analysis_assembler = str_to_lower(coalesce(
      if ("assembler" %in% names(.)) as.character(.data$assembler) else NA_character_,
      if ("Assembler" %in% names(.)) as.character(.data$Assembler) else NA_character_
    )),
    QC_PASS = as_flag(.data$QC_PASS),
    selected_canonical = as_flag(.data$selected_canonical),
    Collection_Date_parsed = parse_date_multi(.data$Collection_Date)
  )

cohort_checks <- tibble(
  check = c(
    "selected_episode_rows", "selected_residents", "operational_uti_rows",
    "operational_not_uti_rows", "unique_episode_keys", "selected_assembler_only",
    "all_qc_pass", "all_selected", "all_fasta_exist"
  ),
  observed = c(
    nrow(cohort), n_distinct(cohort$Participant_id), sum(cohort$UTI_Status == "UTI", na.rm = TRUE),
    sum(cohort$UTI_Status == "Not_UTI", na.rm = TRUE), n_distinct(cohort$episode_key),
    sum(cohort$analysis_assembler == "longcycler", na.rm = TRUE), sum(cohort$QC_PASS, na.rm = TRUE),
    sum(cohort$selected_canonical, na.rm = TRUE), sum(file.exists(cohort$full_path))
  ),
  expected = c(
    EXPECTED_EPISODES, EXPECTED_RESIDENTS, EXPECTED_UTI, EXPECTED_NOT_UTI,
    EXPECTED_EPISODES, EXPECTED_EPISODES, EXPECTED_EPISODES, EXPECTED_EPISODES,
    EXPECTED_EPISODES
  )
) %>% mutate(passed = .data$observed == .data$expected)
if (anyNA(cohort$Collection_Date_parsed)) {
  cohort_checks <- bind_rows(cohort_checks, tibble(
    check = "all_collection_dates_parse", observed = sum(!is.na(cohort$Collection_Date_parsed)),
    expected = EXPECTED_EPISODES, passed = FALSE
  ))
}
atomic_write_csv(cohort_checks, file.path(input_root, "rq09_10_selected_cohort_checks.csv"))
if (!all(cohort_checks$passed)) stop("Selected-cohort invariants failed.", call. = FALSE)

manifest <- read_csv(required_inputs[["manifest"]], show_col_types = FALSE) %>%
  mutate(
    Participant_id = as.character(.data$Participant_id),
    tp_lab = as.character(.data$tp_lab),
    episode_key = episode_key(.data$Participant_id, .data$tp_lab),
    full_path = norm_path(.data$full_path)
  )
if (nrow(manifest) != EXPECTED_EPISODES || anyDuplicated(manifest$episode_key) ||
    !setequal(manifest$episode_key, cohort$episode_key)) {
  stop("Analysis manifest does not exactly match the selected clinical cohort.", call. = FALSE)
}
manifest_match <- match(cohort$episode_key, manifest$episode_key)
if (anyNA(manifest_match) || any(cohort$full_path != manifest$full_path[manifest_match])) {
  stop("Selected cohort and manifest disagree on exact FASTA paths.", call. = FALSE)
}

message("Hashing ", EXPECTED_EPISODES, " selected FASTAs with ", workers, " worker(s).")
cohort$fasta_sha256 <- unlist(parallel_lapply(cohort$full_path, sha256_file), use.names = FALSE)
if (any(nchar(cohort$fasta_sha256) != 64L)) stop("Could not hash every selected FASTA.")
atomic_write_csv(
  cohort %>% select(
    episode_key, Assembly_ID, Isolate_ID, Participant_id,
    tp_lab, UTI_Status, Event_type, full_path, fasta_sha256
  ),
  file.path(input_root, "rq09_10_selected_fasta_sha256.csv")
)

direct_pairs <- read_csv(required_inputs[["direct_pairs"]], show_col_types = FALSE)
require_columns(
  direct_pairs,
  c("Participant_id", "key_A", "key_B", "pair_id", "fasta_sha256_A", "fasta_sha256_B"),
  "direct_pair_metrics_893"
)
expected_pairs <- bind_rows(lapply(split(cohort, cohort$Participant_id), function(z) {
  if (nrow(z) < 2L) return(tibble(pair_id = character()))
  cmb <- combn(z$episode_key, 2L)
  tibble(pair_id = unordered_pair_id(cmb[1, ], cmb[2, ]))
}))
if (nrow(expected_pairs) != EXPECTED_ALL_PAIRS || nrow(direct_pairs) != EXPECTED_ALL_PAIRS ||
    anyDuplicated(direct_pairs$pair_id) || !setequal(direct_pairs$pair_id, expected_pairs$pair_id)) {
  stop("Direct-pair table is not the exact 893-pair selected-cohort universe.", call. = FALSE)
}
idx_a <- match(direct_pairs$key_A, cohort$episode_key)
idx_b <- match(direct_pairs$key_B, cohort$episode_key)
if (anyNA(idx_a) || anyNA(idx_b) ||
    any(str_to_lower(direct_pairs$fasta_sha256_A) != str_to_lower(cohort$fasta_sha256[idx_a])) ||
    any(str_to_lower(direct_pairs$fasta_sha256_B) != str_to_lower(cohort$fasta_sha256[idx_b]))) {
  stop("Direct-pair endpoints are not hash-bound to the selected cohort.", call. = FALSE)
}

provider_all <- read_csv(required_inputs[["provider_mlst"]], show_col_types = FALSE)
require_columns(provider_all, c("full_path", "ST_source", "ST_provider"), "provider MLST")
provider_subset <- provider_all %>%
  transmute(
    full_path = norm_path(.data$full_path),
    ST_source = as.character(.data$ST_source),
    ST_provider = as.character(.data$ST_provider)
  ) %>%
  filter(.data$full_path %in% cohort$full_path)
if (anyDuplicated(provider_subset$full_path)) {
  stop("Provider MLST has duplicate rows for a selected FASTA path.", call. = FALSE)
}
provider_st <- cohort %>%
  select(
    episode_key, Assembly_ID, Participant_id, tp_lab,
    Event_type, Collection_Date_parsed, full_path
  ) %>%
  left_join(provider_subset, by = "full_path", relationship = "one-to-one") %>%
  filter(
    .data$ST_source == "provider_qc95", !is.na(.data$ST_provider),
    nzchar(str_trim(.data$ST_provider))
  ) %>%
  mutate(ST = str_trim(.data$ST_provider)) %>%
  select(-ST_provider)
if (!nrow(provider_st)) stop("No key-linked provider ST calls were available.", call. = FALSE)
atomic_write_csv(
  provider_st %>% select(
    episode_key, Assembly_ID, Participant_id, tp_lab,
    Event_type, full_path, ST, ST_source
  ),
  file.path(input_root, "provider_st_keylinked_selected.csv")
)

resident_key <- cohort %>% distinct(.data$Participant_id) %>% arrange(.data$Participant_id) %>%
  mutate(Resident_Label = sprintf("Resident_%03d", row_number()))
episode_label_key <- cohort %>% distinct(.data$episode_key) %>% arrange(.data$episode_key) %>%
  mutate(Episode_Label = sprintf("Episode_%03d", row_number()))
genome_key <- cohort %>% distinct(.data$Assembly_ID) %>% arrange(.data$Assembly_ID) %>%
  mutate(Genome_Label = sprintf("Genome_%03d", row_number()))
atomic_write_csv(resident_key, file.path(input_root, "private_resident_label_key.csv"))
atomic_write_csv(episode_label_key, file.path(input_root, "private_episode_label_key.csv"))
atomic_write_csv(genome_key, file.path(input_root, "private_genome_label_key.csv"))

valid_mash_sidecar <- function(sidecar_path, signature, output_path, sketch_path) {
  if (!file.exists(sidecar_path) || !file.exists(output_path) || !file.exists(sketch_path)) return(FALSE)
  side <- tryCatch(read_json(sidecar_path, simplifyVector = TRUE), error = function(e) NULL)
  if (is.null(side)) return(FALSE)
  identical(as.character(side$signature), signature) &&
    identical(as.character(side$output_sha256), sha256_file(output_path)) &&
    isTRUE(as.integer(side$exit_code) == 0L)
}

weighted_median <- function(x, w) {
  keep <- is.finite(x) & is.finite(w) & w > 0
  if (!any(keep)) return(NA_real_)
  x <- x[keep]
  w <- w[keep]
  ord <- order(x)
  x <- x[ord]
  w <- w[ord]
  x[[which(cumsum(w) >= sum(w) / 2)[[1]]]]
}

run_rq09 <- function() {
  rq <- "RQ09"
  out_dir <- rq_dirs[[rq]]
  mash <- Sys.which("mash")
  if (mash == "") {
    write_status(rq, "blocked", "mash_not_available", "No result was fabricated.")
    return(FALSE)
  }
  mash_version <- tool_version("mash")

  eligible_st <- provider_st %>%
    count(.data$ST, name = "n_genomes") %>%
    left_join(
      provider_st %>% distinct(.data$ST, .data$Participant_id) %>% count(.data$ST, name = "n_residents"),
      by = "ST"
    ) %>%
    filter(.data$n_genomes >= 5L, .data$n_residents >= 2L) %>%
    arrange(desc(.data$n_genomes), .data$ST)
  atomic_write_csv(eligible_st, file.path(out_dir, "eligible_provider_st.csv"))
  if (!nrow(eligible_st)) {
    write_status(rq, "blocked", "no_provider_st_meets_size_rule")
    return(FALSE)
  }

  mash_input <- cohort %>% arrange(.data$Assembly_ID)
  signature <- digest(
    paste(
      "rq09_mash_selected_v2", mash_version, "k=21", "s=1000",
      paste(mash_input$Assembly_ID, mash_input$fasta_sha256, sep = "=", collapse = "\n"),
      sep = "\n"
    ),
    algo = "sha256", serialize = FALSE
  )
  mash_dir <- file.path(input_root, "rq09_mash", signature)
  dir.create(mash_dir, recursive = TRUE, showWarnings = FALSE)
  list_path <- file.path(mash_dir, "selected_fasta_paths.txt")
  sketch_prefix <- file.path(mash_dir, "selected_532")
  sketch_path <- paste0(sketch_prefix, ".msh")
  dist_path <- file.path(mash_dir, "all_vs_all_dist.tsv")
  sidecar_path <- file.path(mash_dir, "provenance.json")
  atomic_write_lines(mash_input$full_path, list_path)

  cache_ok <- valid_mash_sidecar(sidecar_path, signature, dist_path, sketch_path)
  if (!cache_ok) {
    sketch_tmp_prefix <- file.path(mash_dir, paste0("selected_532.tmp.", Sys.getpid()))
    sketch_tmp <- paste0(sketch_tmp_prefix, ".msh")
    sketch_err <- file.path(mash_dir, "sketch.stderr.txt")
    sketch_status <- run_system2(
      mash,
      c("sketch", "-p", as.character(workers), "-o", sketch_tmp_prefix, "-l", list_path),
      stdout = file.path(mash_dir, "sketch.stdout.txt"), stderr = sketch_err
    )
    if (sketch_status != 0L || !file.exists(sketch_tmp)) {
      write_status(rq, "blocked", "mash_sketch_failed")
      return(FALSE)
    }
    file.copy(sketch_tmp, sketch_path, overwrite = TRUE)
    unlink(sketch_tmp)

    dist_tmp <- paste0(dist_path, ".tmp.", Sys.getpid())
    dist_err <- file.path(mash_dir, "dist.stderr.txt")
    dist_status <- run_system2(
      mash,
      c("dist", "-p", as.character(workers), sketch_path, sketch_path),
      stdout = dist_tmp, stderr = dist_err
    )
    if (dist_status != 0L || !file.exists(dist_tmp) || file.size(dist_tmp) == 0L) {
      write_status(rq, "blocked", "mash_dist_failed")
      return(FALSE)
    }
    file.copy(dist_tmp, dist_path, overwrite = TRUE)
    unlink(dist_tmp)
    atomic_write_json(list(
      signature = signature,
      created_utc = run_timestamp,
      mash_executable = mash,
      mash_version = mash_version,
      parameters = list(k = 21L, sketch_size = 1000L),
      n_genomes = nrow(mash_input),
      manifest_sha256 = sha256_file(required_inputs[["manifest"]]),
      output_sha256 = sha256_file(dist_path),
      exit_code = 0L
    ), sidecar_path)
  }

  mash_raw <- read_tsv(
    dist_path,
    col_names = c("reference", "query", "mash_distance", "p_value", "shared_hashes"),
    col_types = cols(
      reference = col_character(), query = col_character(), mash_distance = col_double(),
      p_value = col_double(), shared_hashes = col_character()
    ),
    show_col_types = FALSE, progress = FALSE
  ) %>% mutate(reference = norm_path(.data$reference), query = norm_path(.data$query))
  path_index <- setNames(seq_len(nrow(mash_input)), mash_input$full_path)
  ref_i <- unname(path_index[mash_raw$reference])
  qry_i <- unname(path_index[mash_raw$query])
  if (anyNA(ref_i) || anyNA(qry_i)) stop("Mash output contains an unselected endpoint path.")
  keep <- ref_i < qry_i
  mash_unique <- mash_raw[keep, ] %>% transmute(
    Assembly_ID_A = mash_input$Assembly_ID[ref_i[keep]],
    Assembly_ID_B = mash_input$Assembly_ID[qry_i[keep]],
    .data$mash_distance, .data$p_value, .data$shared_hashes
  )

  st_meta <- provider_st %>% select(
    Assembly_ID, Participant_id, Collection_Date_parsed, ST
  )
  pairs <- mash_unique %>%
    inner_join(
      st_meta %>% rename(
        Assembly_ID_A = Assembly_ID, Participant_id_A = Participant_id,
        date_A = Collection_Date_parsed, ST_A = ST
      ), by = "Assembly_ID_A"
    ) %>%
    inner_join(
      st_meta %>% rename(
        Assembly_ID_B = Assembly_ID, Participant_id_B = Participant_id,
        date_B = Collection_Date_parsed, ST_B = ST
      ), by = "Assembly_ID_B"
    ) %>%
    filter(.data$ST_A == .data$ST_B, .data$ST_A %in% eligible_st$ST) %>%
    transmute(
      ST = .data$ST_A, .data$Assembly_ID_A, .data$Assembly_ID_B,
      .data$Participant_id_A, .data$Participant_id_B,
      days_between = abs(as.numeric(.data$date_B - .data$date_A)),
      same_resident = .data$Participant_id_A == .data$Participant_id_B,
      .data$mash_distance, .data$p_value, .data$shared_hashes
    )

  estimate_effect <- function(dat) {
    by_st <- dat %>%
      group_by(.data$ST, .data$same_resident) %>%
      summarise(n_pairs = n(), median_mash = median(.data$mash_distance), .groups = "drop") %>%
      pivot_wider(
        names_from = same_resident, values_from = c(n_pairs, median_mash),
        names_glue = "{.value}_{ifelse(same_resident, 'within', 'between')}"
      ) %>%
      filter(!is.na(.data$median_mash_within), !is.na(.data$median_mash_between)) %>%
      mutate(within_minus_between = .data$median_mash_within - .data$median_mash_between)
    list(
      by_st = by_st,
      estimate = if (nrow(by_st)) mean(by_st$within_minus_between) else NA_real_,
      n_st = nrow(by_st)
    )
  }

  primary_pairs <- pairs %>% filter(!is.na(.data$days_between), .data$days_between <= 180)
  primary <- estimate_effect(primary_pairs)
  if (!is.finite(primary$estimate) || primary$n_st == 0L) {
    write_status(rq, "blocked", "no_st_has_both_pair_types_within_180_days")
    return(FALSE)
  }
  effect_st <- primary$by_st$ST
  pair_split <- split(filter(primary_pairs, .data$ST %in% effect_st), ~ ST)
  genome_split <- split(filter(provider_st, .data$ST %in% effect_st), ~ ST)

  set.seed(seed + 9L)
  perm_values <- replicate(n_perm, {
    differences <- map_dbl(effect_st, function(st) {
      p <- pair_split[[st]]
      g <- genome_split[[st]]
      labels <- sample(g$Participant_id, nrow(g), replace = FALSE)
      names(labels) <- g$Assembly_ID
      within <- labels[p$Assembly_ID_A] == labels[p$Assembly_ID_B]
      if (!any(within) || all(within)) return(NA_real_)
      median(p$mash_distance[within]) - median(p$mash_distance[!within])
    })
    if (all(is.na(differences))) NA_real_ else mean(differences, na.rm = TRUE)
  })
  permutation_p <- (1 + sum(abs(perm_values) >= abs(primary$estimate), na.rm = TRUE)) /
    (1 + sum(is.finite(perm_values)))

  bootstrap_residents <- sort(unique(c(
    primary_pairs$Participant_id_A, primary_pairs$Participant_id_B
  )))
  set.seed(seed + 90L)
  bootstrap_tbl <- map_dfr(seq_len(n_boot), function(iteration) {
    # Positive Exp(1) resident multipliers avoid the zero-cell failures of the
    # superseded resident-frequency resample while retaining resident-level
    # clustering. A between-resident pair receives the product of its two
    # endpoint multipliers; a within-resident pair receives its resident's
    # multiplier.
    resident_weight <- setNames(
      rexp(length(bootstrap_residents), rate = 1),
      bootstrap_residents
    )
    differences <- map_dbl(effect_st, function(st) {
      p <- pair_split[[st]]
      w_a <- resident_weight[p$Participant_id_A]
      w_b <- resident_weight[p$Participant_id_B]
      weight <- ifelse(p$same_resident, w_a, w_a * w_b)
      weighted_median(p$mash_distance[p$same_resident], weight[p$same_resident]) -
        weighted_median(p$mash_distance[!p$same_resident], weight[!p$same_resident])
    })
    valid <- sum(is.finite(differences))
    tibble(
      iteration = iteration,
      bootstrap_method = "Exp(1) resident-weight multiplier",
      valid_st = valid,
      all_primary_st_valid = valid == length(effect_st),
      equal_weight_difference = if (valid == length(effect_st)) mean(differences) else NA_real_
    )
  })
  boot_values <- bootstrap_tbl$equal_weight_difference
  valid_bootstrap_replicates <- sum(is.finite(boot_values))
  required_valid_bootstrap <- if (n_boot == 10000L) {
    9990L
  } else {
    ceiling(0.999 * n_boot)
  }
  if (valid_bootstrap_replicates < required_valid_bootstrap) {
    stop(
      "RQ09 exponential resident-weight multiplier bootstrap produced only ",
      valid_bootstrap_replicates, "/", n_boot, " valid replicates; required ",
      required_valid_bootstrap, ".",
      call. = FALSE
    )
  }

  primary_display <- primary_pairs %>%
    left_join(genome_key %>% rename(Assembly_ID_A = Assembly_ID, Genome_Label_A = Genome_Label), by = "Assembly_ID_A") %>%
    left_join(genome_key %>% rename(Assembly_ID_B = Assembly_ID, Genome_Label_B = Genome_Label), by = "Assembly_ID_B") %>%
    left_join(resident_key %>% rename(Participant_id_A = Participant_id, Resident_Label_A = Resident_Label), by = "Participant_id_A") %>%
    left_join(resident_key %>% rename(Participant_id_B = Participant_id, Resident_Label_B = Resident_Label), by = "Participant_id_B") %>%
    transmute(
      Pair_Label = sprintf("RQ09_Pair_%05d", row_number()), .data$ST,
      .data$Genome_Label_A, .data$Genome_Label_B, .data$Resident_Label_A,
      .data$Resident_Label_B, .data$days_between, .data$same_resident,
      .data$mash_distance, .data$p_value, .data$shared_hashes
    )
  atomic_write_csv(primary_display, file.path(out_dir, "same_st_pairs_within_180_days.csv"))
  atomic_write_csv(primary$by_st, file.path(out_dir, "st_equal_weight_effects_180_days.csv"))
  atomic_write_csv(
    tibble(iteration = seq_len(n_perm), equal_weight_difference = perm_values),
    file.path(out_dir, "resident_label_permutation_distribution.csv")
  )
  atomic_write_csv(
    bootstrap_tbl,
    file.path(
      out_dir,
      "resident_exponential_multiplier_bootstrap_distribution.csv"
    )
  )
  atomic_write_csv(tibble(
    estimand = "Mean across provider STs of within-resident median Mash distance minus between-resident median Mash distance",
    selected_genomes = nrow(cohort),
    provider_typed_genomes = nrow(provider_st),
    eligible_st = nrow(eligible_st),
    st_in_estimand = primary$n_st,
    pairs_within_180_days = nrow(primary_pairs),
    observed_difference = primary$estimate,
    primary_inference = "resident-label permutation",
    confidence_interval_method =
      "Exp(1) resident-weight multiplier bootstrap",
    bootstrap_ci_low = unname(quantile(boot_values, 0.025, na.rm = TRUE)),
    bootstrap_ci_high = unname(quantile(boot_values, 0.975, na.rm = TRUE)),
    bootstrap_valid_replicates = valid_bootstrap_replicates,
    bootstrap_failed_replicates = sum(!is.finite(boot_values)),
    bootstrap_minimum_valid_required = required_valid_bootstrap,
    bootstrap_validity_gate_passed =
      valid_bootstrap_replicates >= required_valid_bootstrap,
    permutation_p_two_sided = permutation_p,
    permutations = n_perm,
    bootstrap_replicates = n_boot
  ), file.path(out_dir, "primary_result.csv"))
  atomic_write_csv(tibble(
    input = c("selected_clinical_cohort", "selected_assembly_manifest", "key_linked_provider_mlst", "mash_distance_cache"),
    path = c(required_inputs[["cohort"]], required_inputs[["manifest"]], required_inputs[["provider_mlst"]], dist_path),
    sha256 = c(
      sha256_file(required_inputs[["cohort"]]), sha256_file(required_inputs[["manifest"]]),
      sha256_file(required_inputs[["provider_mlst"]]), sha256_file(dist_path)
    ),
    role = c("exact_532_episode_analysis_cohort", "selected_qc_passing_fastas", "exact_selected_path_subset", "sha_bound_all_vs_all_mash")
  ), file.path(out_dir, "input_provenance.csv"))
  atomic_write_lines(c(
    "# RQ09 — same-provider-ST Mash distance",
    "",
    "Status: complete.",
    "",
    "The analysis is restricted to the exact 532 selected episodes from 161 residents.",
    "Provider calls are accepted only by exact selected FASTA path; no assembler field is required or published.",
    "STs require at least five typed genomes and two residents. The primary window is 180 days.",
    "Negative estimates mean lower Mash distance within residents. The resident-label",
    "permutation is primary; the confidence interval uses an Exp(1) resident-weight",
    "multiplier bootstrap. Final inference uses 10,000 replicates and requires at least",
    "9,990 valid bootstrap estimates. Clinical UTI_Status is the frozen",
    "operational phenotype and is not used as a causal exposure in this question."
  ), file.path(out_dir, "README.md"))
  write_status(rq, "complete", "analysis_completed", sprintf("%d selected genomes; Mash %s", nrow(cohort), mash_version))
  TRUE
}

validate_dnadiff_sidecar <- function(row) {
  side_path <- as.character(row$dnadiff_sidecar_path[[1]])
  report_path <- as.character(row$dnadiff_report_path[[1]])
  if (!file.exists(side_path) || !file.exists(report_path)) return(FALSE)
  side <- tryCatch(read_json(side_path, simplifyVector = TRUE), error = function(e) NULL)
  if (is.null(side) || is.null(side$input_a) || is.null(side$input_b) || is.null(side$report)) return(FALSE)
  live_report_sha <- sha256_file(report_path)
  identical(as.character(side$cache_signature), as.character(row$dnadiff_cache_signature[[1]])) &&
    identical(norm_path(side$input_a$path), norm_path(row$Fasta_path_A[[1]])) &&
    identical(norm_path(side$input_b$path), norm_path(row$Fasta_path_B[[1]])) &&
    identical(str_to_lower(as.character(side$input_a$sha256)), str_to_lower(as.character(row$Fasta_SHA256_A[[1]]))) &&
    identical(str_to_lower(as.character(side$input_b$sha256)), str_to_lower(as.character(row$Fasta_SHA256_B[[1]]))) &&
    identical(norm_path(side$report$path), norm_path(report_path)) &&
    identical(str_to_lower(as.character(side$report$sha256)), str_to_lower(live_report_sha)) &&
    identical(str_to_lower(as.character(row$dnadiff_report_sha256[[1]])), str_to_lower(live_report_sha))
}

run_rq10 <- function() {
  rq <- "RQ10"
  out_dir <- rq_dirs[[rq]]
  transitions <- read_csv(required_inputs[["transitions"]], show_col_types = FALSE)
  require_columns(transitions, c(
    "Participant_id", "tp_from", "tp_to", "status_from", "status_to",
    "days_between_samples", "TotalSNPs", "fasta_path_from", "fasta_path_to",
    "Fasta_path_A", "Fasta_path_B", "Fasta_SHA256_A", "Fasta_SHA256_B",
    "dnadiff_report_path", "dnadiff_sidecar_path", "dnadiff_report_sha256",
    "dnadiff_cache_signature", "dnadiff_cache_status", "dnadiff_version"
  ), "selected adjacent transitions")
  transitions <- transitions %>% mutate(
    Participant_id = as.character(.data$Participant_id),
    tp_from = as.character(.data$tp_from),
    tp_to = as.character(.data$tp_to),
    key_from = episode_key(.data$Participant_id, .data$tp_from),
    key_to = episode_key(.data$Participant_id, .data$tp_to),
    fasta_path_from = norm_path(.data$fasta_path_from),
    fasta_path_to = norm_path(.data$fasta_path_to),
    Fasta_path_A = norm_path(.data$Fasta_path_A),
    Fasta_path_B = norm_path(.data$Fasta_path_B),
    days_between_samples = as.numeric(.data$days_between_samples),
    TotalSNPs = as.numeric(.data$TotalSNPs)
  )

  idx_from <- match(transitions$key_from, cohort$episode_key)
  idx_to <- match(transitions$key_to, cohort$episode_key)
  endpoint_valid <- !is.na(idx_from) & !is.na(idx_to) &
    transitions$fasta_path_from == cohort$full_path[idx_from] &
    transitions$fasta_path_to == cohort$full_path[idx_to]
  hash_idx_a <- match(transitions$Fasta_path_A, cohort$full_path)
  hash_idx_b <- match(transitions$Fasta_path_B, cohort$full_path)
  hash_valid <- !is.na(hash_idx_a) & !is.na(hash_idx_b) &
    !is.na(transitions$Fasta_SHA256_A) & !is.na(transitions$Fasta_SHA256_B) &
    str_to_lower(transitions$Fasta_SHA256_A) == str_to_lower(cohort$fasta_sha256[hash_idx_a]) &
    str_to_lower(transitions$Fasta_SHA256_B) == str_to_lower(cohort$fasta_sha256[hash_idx_b])
  cache_status_valid <- transitions$dnadiff_cache_status %in% c("fresh", "generated", "reused")
  sidecar_valid <- vapply(
    seq_len(nrow(transitions)),
    function(i) validate_dnadiff_sidecar(transitions[i, , drop = FALSE]),
    logical(1)
  )
  transition_checks <- tibble(
    check = c(
      "selected_adjacent_pairs", "selected_adjacent_residents", "selected_endpoint_paths",
      "selected_endpoint_hashes", "validated_cache_status", "validated_sidecar_and_report",
      "snp_values_present", "pairs_at_operational_25_snp_boundary",
      "not_uti_to_uti_pairs", "not_uti_to_uti_pairs_at_operational_boundary",
      "all_direct_within_resident_pairs"
    ),
    observed = c(
      nrow(transitions), n_distinct(transitions$Participant_id), sum(endpoint_valid),
      sum(hash_valid, na.rm = TRUE), sum(cache_status_valid, na.rm = TRUE), sum(sidecar_valid), sum(!is.na(transitions$TotalSNPs)),
      sum(transitions$TotalSNPs <= 25, na.rm = TRUE),
      sum(transitions$status_from == "Not_UTI" & transitions$status_to == "UTI"),
      sum(transitions$status_from == "Not_UTI" & transitions$status_to == "UTI" & transitions$TotalSNPs <= 25, na.rm = TRUE),
      nrow(direct_pairs)
    ),
    expected = c(
      EXPECTED_ADJACENT_PAIRS, EXPECTED_ADJACENT_RESIDENTS, EXPECTED_ADJACENT_PAIRS,
      EXPECTED_ADJACENT_PAIRS, EXPECTED_ADJACENT_PAIRS, EXPECTED_ADJACENT_PAIRS,
      EXPECTED_ADJACENT_PAIRS, EXPECTED_ADJACENT_LE25, EXPECTED_NOT_UTI_TO_UTI,
      EXPECTED_NOT_UTI_TO_UTI_LE25, EXPECTED_ALL_PAIRS
    )
  ) %>% mutate(passed = .data$observed == .data$expected)
  atomic_write_csv(transition_checks, file.path(out_dir, "input_checks.csv"))
  if (!all(transition_checks$passed)) {
    write_status(rq, "blocked", "selected_transition_provenance_failed")
    return(FALSE)
  }

  st_key <- provider_st %>% select(
    episode_key, ST, Assembly_ID, Event_type
  )
  adjacent_all <- transitions %>%
    left_join(
      st_key %>% rename(
        key_from = episode_key, ST_from_provider = ST,
        Assembly_ID_from_provider = Assembly_ID, event_type_from = Event_type
      ), by = "key_from", relationship = "many-to-one"
    ) %>%
    left_join(
      st_key %>% rename(
        key_to = episode_key, ST_to_provider = ST,
        Assembly_ID_to_provider = Assembly_ID, event_type_to = Event_type
      ), by = "key_to", relationship = "many-to-one"
    ) %>%
    mutate(
      provider_ST_from_available = !is.na(.data$ST_from_provider),
      provider_ST_to_available = !is.na(.data$ST_to_provider),
      provider_ST_both_available = .data$provider_ST_from_available & .data$provider_ST_to_available,
      ST_changed = .data$ST_from_provider != .data$ST_to_provider,
      routine_both = .data$event_type_from == "Routine" & .data$event_type_to == "Routine",
      within_365_days = .data$days_between_samples <= 365,
      ST_relation = if_else(.data$ST_changed, "ST_changed", "ST_same"),
      SNP_relation = if_else(.data$TotalSNPs <= 25, "SNP_le_25", "SNP_gt_25")
    )
  atomic_write_csv(bind_rows(
    tibble(stage = "selected_adjacent_transitions", n_pairs = nrow(adjacent_all)),
    tibble(stage = "provider_ST_available_at_from_endpoint", n_pairs = sum(adjacent_all$provider_ST_from_available)),
    tibble(stage = "provider_ST_available_at_to_endpoint", n_pairs = sum(adjacent_all$provider_ST_to_available)),
    tibble(stage = "provider_ST_available_at_both_endpoints", n_pairs = sum(adjacent_all$provider_ST_both_available)),
    tibble(stage = "excluded_missing_provider_ST_at_either_endpoint", n_pairs = sum(!adjacent_all$provider_ST_both_available))
  ), file.path(out_dir, "provider_st_endpoint_attrition.csv"))
  adjacent <- adjacent_all %>% filter(.data$provider_ST_both_available)
  if (!nrow(adjacent) || nrow(adjacent) > EXPECTED_ADJACENT_PAIRS) {
    write_status(rq, "blocked", "unexpected_provider_complete_case_denominator")
    return(FALSE)
  }

  fit_spline <- function(data) {
    data <- filter(data, is.finite(.data$days_between_samples), !is.na(.data$ST_changed))
    if (nrow(data) < 20L || n_distinct(data$ST_changed) < 2L ||
        n_distinct(data$days_between_samples) < 5L) return(NULL)
    suppressWarnings(tryCatch(
      glm(ST_changed ~ splines::ns(days_between_samples, df = 3), data = data, family = binomial()),
      error = function(e) NULL
    ))
  }
  analyses <- list(
    primary = adjacent,
    routine_both = filter(adjacent, .data$routine_both),
    within_365_days = filter(adjacent, .data$within_365_days)
  )
  models <- map(analyses, fit_spline)
  primary_model <- models[["primary"]]
  if (is.null(primary_model)) {
    write_status(rq, "blocked", "primary_time_model_not_estimable")
    return(FALSE)
  }
  coefficient_table <- imap_dfr(models, function(model, analysis) {
    dat <- analyses[[analysis]]
    if (is.null(model)) return(tibble(
      analysis = analysis, n_pairs = nrow(dat), n_residents = n_distinct(dat$Participant_id),
      term = NA_character_, estimate_log_odds = NA_real_, standard_error = NA_real_,
      z_value = NA_real_, p_value = NA_real_, model_status = "not_estimable"
    ))
    sm <- coef(summary(model))
    tibble(
      analysis = analysis, n_pairs = nrow(dat), n_residents = n_distinct(dat$Participant_id),
      term = rownames(sm), estimate_log_odds = sm[, "Estimate"], standard_error = sm[, "Std. Error"],
      z_value = sm[, "z value"], p_value = sm[, "Pr(>|z|)"],
      model_status = if_else(model$converged, "converged", "not_converged")
    )
  })

  time_grid <- tibble(days_between_samples = seq(
    min(adjacent$days_between_samples), max(adjacent$days_between_samples), length.out = 61L
  ))
  primary_curve <- time_grid %>% mutate(
    predicted_turnover_probability = as.numeric(predict(primary_model, newdata = time_grid, type = "response"))
  )
  residents <- unique(adjacent$Participant_id)
  set.seed(seed + 10L)
  bootstrap_samples <- replicate(n_boot, sample(residents, length(residents), replace = TRUE), simplify = FALSE)
  bootstrap_results <- parallel_lapply(seq_len(n_boot), function(iteration) {
    sampled <- bootstrap_samples[[iteration]]
    boot_data <- map2_dfr(sampled, seq_along(sampled), function(pid, draw) {
      filter(adjacent, .data$Participant_id == pid) %>% mutate(bootstrap_cluster = draw)
    })
    model <- fit_spline(boot_data)
    list(
      proportion = mean(boot_data$ST_changed),
      curve = if (is.null(model)) rep(NA_real_, nrow(time_grid)) else
        as.numeric(predict(model, newdata = time_grid, type = "response"))
    )
  })
  bootstrap_proportions <- tibble(
    iteration = seq_len(n_boot),
    turnover_proportion = map_dbl(bootstrap_results, "proportion")
  )
  curve_matrix <- do.call(rbind, map(bootstrap_results, "curve"))
  curve_summary <- primary_curve %>% mutate(
    bootstrap_median = apply(curve_matrix, 2, median, na.rm = TRUE),
    bootstrap_ci_low = apply(curve_matrix, 2, quantile, probs = 0.025, na.rm = TRUE),
    bootstrap_ci_high = apply(curve_matrix, 2, quantile, probs = 0.975, na.rm = TRUE),
    valid_bootstrap_models = colSums(is.finite(curve_matrix))
  )

  adjacent_display <- adjacent %>%
    left_join(resident_key, by = "Participant_id") %>%
    left_join(episode_label_key %>% rename(key_from = episode_key, Episode_Label_From = Episode_Label), by = "key_from") %>%
    left_join(episode_label_key %>% rename(key_to = episode_key, Episode_Label_To = Episode_Label), by = "key_to") %>%
    transmute(
      Transition_Label = sprintf("RQ10_Transition_%03d", row_number()), .data$Resident_Label,
      .data$Episode_Label_From, .data$Episode_Label_To, ST_From = .data$ST_from_provider,
      ST_To = .data$ST_to_provider, .data$ST_changed, .data$event_type_from,
      .data$event_type_to, .data$routine_both, .data$days_between_samples,
      .data$within_365_days, .data$TotalSNPs, .data$SNP_relation, .data$ST_relation,
      .data$dnadiff_cache_status, .data$dnadiff_version
    )
  atomic_write_csv(adjacent_display, file.path(out_dir, "provider_only_adjacent_st_pairs.csv"))
  atomic_write_csv(coefficient_table, file.path(out_dir, "spline_model_coefficients.csv"))
  atomic_write_csv(imap_dfr(analyses, function(dat, analysis) tibble(
    analysis = analysis, n_pairs = nrow(dat), n_residents = n_distinct(dat$Participant_id),
    n_st_changed = sum(dat$ST_changed), turnover_proportion = mean(dat$ST_changed),
    median_days = median(dat$days_between_samples), max_days = max(dat$days_between_samples),
    model_estimable = !is.null(models[[analysis]])
  )), file.path(out_dir, "turnover_sensitivity_summary.csv"))
  atomic_write_csv(bootstrap_proportions, file.path(out_dir, "resident_bootstrap_turnover_proportion.csv"))
  atomic_write_csv(curve_summary, file.path(out_dir, "resident_bootstrap_time_curve.csv"))
  atomic_write_csv(tibble(
    estimand = "Proportion of adjacent selected-cohort pairs with different provider-qc95 STs",
    complete_case_pairs = nrow(adjacent), residents = n_distinct(adjacent$Participant_id),
    st_changes = sum(adjacent$ST_changed), observed_turnover_proportion = mean(adjacent$ST_changed),
    bootstrap_ci_low = quantile(bootstrap_proportions$turnover_proportion, 0.025, na.rm = TRUE),
    bootstrap_ci_high = quantile(bootstrap_proportions$turnover_proportion, 0.975, na.rm = TRUE),
    requested_bootstrap_replicates = n_boot,
    valid_bootstrap_replicates = sum(is.finite(bootstrap_proportions$turnover_proportion)),
    point_spline_model_converged = isTRUE(primary_model$converged)
  ), file.path(out_dir, "primary_turnover_estimate.csv"))
  atomic_write_csv(
    adjacent %>% count(.data$ST_relation, .data$SNP_relation, name = "n_transitions") %>%
      group_by(.data$ST_relation) %>% mutate(row_proportion = .data$n_transitions / sum(.data$n_transitions)) %>% ungroup(),
    file.path(out_dir, "st_by_snp_cross_tab.csv")
  )
  atomic_write_csv(tibble(
    input = c("selected_clinical_cohort", "selected_assembly_manifest", "key_linked_provider_mlst", "validated_adjacent_transitions", "direct_pair_universe"),
    path = unname(required_inputs[c("cohort", "manifest", "provider_mlst", "transitions", "direct_pairs")]),
    sha256 = map_chr(unname(required_inputs[c("cohort", "manifest", "provider_mlst", "transitions", "direct_pairs")]), sha256_file),
    role = c(
      "exact_532_episode_analysis_cohort", "selected_qc_passing_fastas", "exact_selected_path_subset",
      "371_sha_and_sidecar_validated_adjacent_pairs", "893_hash_bound_within_resident_pairs"
    )
  ), file.path(out_dir, "input_provenance.csv"))
  atomic_write_lines(c(
    "# RQ10 — adjacent provider-ST turnover",
    "",
    "Status: complete.",
    "",
    "The denominator is the exact 371 adjacent pairs rebuilt within the selected 532-episode cohort.",
    "All 893 direct within-resident pairs are checked as the complete longitudinal pair universe.",
    "Every adjacent endpoint path and FASTA hash must match the selected cohort. DNAdiff evidence",
    "is accepted only when its report hash, cache signature, endpoint hashes, and JSON sidecar agree.",
    "Validated fresh, generated, and reused cache states are allowed; no legacy path is reconstructed.",
    "Provider STs are linked only by exact selected FASTA path. Missing calls are complete-case exclusions",
    "and do not bridge visits. The 25-SNP rule and UTI_Status are frozen operational definitions.",
    "These descriptive models do not establish transmission, treatment effects, or causality."
  ), file.path(out_dir, "README.md"))
  write_status(
    rq, "complete", "analysis_completed",
    sprintf("%d validated adjacent pairs; %d provider-complete pairs", nrow(transitions), nrow(adjacent))
  )
  TRUE
}

run_safely <- function(rq, fun) {
  tryCatch(
    fun(),
    error = function(e) {
      write_status(rq, "blocked", "unexpected_runtime_error", conditionMessage(e))
      atomic_write_lines(c(
        paste0("# ", rq, " blocked"), "",
        "The analysis stopped rather than suppressing an error or using stale output.",
        paste0("Error: ", conditionMessage(e))
      ), file.path(rq_dirs[[rq]], "README.md"))
      message(rq, " blocked: ", conditionMessage(e))
      FALSE
    }
  )
}

results <- c(
  RQ09 = run_safely("RQ09", run_rq09),
  RQ10 = run_safely("RQ10", run_rq10)
)
atomic_write_csv(tibble(
  research_question = names(results),
  completed_without_block = as.logical(results),
  run_timestamp_utc = run_timestamp,
  selected_episode_rows = nrow(cohort),
  selected_residents = n_distinct(cohort$Participant_id),
  operational_uti_rows = sum(cohort$UTI_Status == "UTI", na.rm = TRUE),
  operational_not_uti_rows = sum(cohort$UTI_Status == "Not_UTI", na.rm = TRUE),
  direct_pair_rows = nrow(direct_pairs),
  provider_typed_rows = nrow(provider_st)
), file.path(input_root, "RQ09_10_run_summary.csv"))
atomic_write_lines(capture.output(sessionInfo()), file.path(input_root, "RQ09_10_sessionInfo.txt"))
message("RQ09-RQ10 runner finished: ", paste(names(results), results, sep = "=", collapse = ", "))

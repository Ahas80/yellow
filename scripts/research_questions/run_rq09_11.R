#!/usr/bin/env Rscript

# RQ09-RQ11 publication-oriented research-question analyses.
#
# This runner is intentionally isolated from the numbered pipeline. It reads the
# current 532-genome Longcycler analysis manifest, never promotes the superseded
# 556-genome mixed-assembler denominator, and writes only beneath
# results/research_questions/{_inputs,RQ09,RQ10,RQ11}.

options(stringsAsFactors = FALSE, warn = 1)

suppressPackageStartupMessages({
  library(digest)
  library(dplyr)
  library(ggplot2)
  library(jsonlite)
  library(purrr)
  library(readr)
  library(stringr)
  library(tidyr)
})

`%||%` <- function(x, y) if (is.null(x) || length(x) == 0L || is.na(x[[1]])) y else x

project_root <- normalizePath(if (file.exists("00_config.R")) "." else file.path("..", ".."),
                              winslash = "/", mustWork = TRUE)
setwd(project_root)

out_root <- file.path("results", "research_questions")
input_root <- file.path(out_root, "_inputs")
rq_dirs <- setNames(file.path(out_root, c("RQ09", "RQ10", "RQ11")), c("RQ09", "RQ10", "RQ11"))
dir.create(input_root, recursive = TRUE, showWarnings = FALSE)
walk(rq_dirs, dir.create, recursive = TRUE, showWarnings = FALSE)
legacy_root_summary <- file.path(out_root, "RQ09_11_run_summary.csv")
if (file.exists(legacy_root_summary)) unlink(legacy_root_summary)

run_timestamp <- format(Sys.time(), tz = "UTC", usetz = TRUE)
workers <- max(1L, as.integer(Sys.getenv("RQ_WORKERS", "8")))
n_perm <- max(1L, as.integer(Sys.getenv("RQ_PERMUTATIONS", "10000")))
n_boot <- max(1L, as.integer(Sys.getenv("RQ_BOOTSTRAP_REPS", "10000")))
seed <- as.integer(Sys.getenv("RQ_SEED", "20260712"))

atomic_write_csv <- function(x, path) {
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
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  tmp <- tempfile(pattern = paste0(".", basename(path), "."), tmpdir = dirname(path))
  writeLines(x, tmp, useBytes = TRUE)
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

sha256_file <- function(path) unname(digest::digest(path, algo = "sha256", file = TRUE))
norm_path <- function(path) normalizePath(path, winslash = "/", mustWork = FALSE)

required_inputs <- c(
  manifest = file.path("results", "qc", "analysis_assembly_manifest.csv"),
  candidates = file.path("results", "qc", "canonical_assembly_selection.csv"),
  provider_mlst = file.path("results", "mlst", "mlst_provider_preferred.csv"),
  longcycler_transitions = file.path("results", "sensitivity", "longcycler_only", "longcycler_transitions.csv")
)
missing_inputs <- required_inputs[!file.exists(required_inputs)]
if (length(missing_inputs)) stop("Missing required input(s): ", paste(missing_inputs, collapse = ", "))

manifest <- read_csv(required_inputs[["manifest"]], show_col_types = FALSE) %>%
  mutate(
    Participant_id = as.character(Participant_id),
    tp_lab = as.character(tp_lab),
    Assembly_ID = as.character(Assembly_ID),
    Isolate_ID = as.character(Isolate_ID),
    full_path = norm_path(full_path),
    Collection_Date_parsed = as.Date(Collection_Date)
  )

manifest_checks <- tibble(
  check = c("row_count_532", "unique_assembly_id", "longcycler_only", "all_qc_pass", "all_fasta_exist"),
  passed = c(
    nrow(manifest) == 532L,
    n_distinct(manifest$Assembly_ID) == 532L,
    all(tolower(manifest$Assembler) == "longcycler"),
    all(manifest$QC_PASS %in% TRUE),
    all(file.exists(manifest$full_path))
  ),
  observed = c(
    as.character(nrow(manifest)),
    as.character(n_distinct(manifest$Assembly_ID)),
    paste(sort(unique(manifest$Assembler)), collapse = ";"),
    as.character(sum(manifest$QC_PASS %in% TRUE)),
    as.character(sum(file.exists(manifest$full_path)))
  )
)
atomic_write_csv(manifest_checks, file.path(input_root, "current_532_manifest_checks.csv"))
if (!all(manifest_checks$passed)) stop("Current analysis manifest failed its 532-genome invariants.")

# Stable display labels keep participant/isolate keys out of the public RQ
# folders. Reproducibility mappings stay in the explicitly internal _inputs
# area. Labels are deterministic for this frozen manifest.
resident_key <- manifest %>%
  distinct(Participant_id) %>%
  arrange(Participant_id) %>%
  mutate(Resident_Label = sprintf("Resident_%03d", row_number()))
episode_key <- manifest %>%
  distinct(Participant_id, tp_lab) %>%
  arrange(Participant_id, tp_lab) %>%
  mutate(Episode_Label = sprintf("Episode_%03d", row_number()))
assembly_key <- manifest %>%
  distinct(Assembly_ID) %>%
  arrange(Assembly_ID) %>%
  mutate(Genome_Label = sprintf("Genome_%03d", row_number()))
atomic_write_csv(resident_key, file.path(input_root, "private_resident_label_key.csv"))
atomic_write_csv(episode_key, file.path(input_root, "private_episode_label_key.csv"))
atomic_write_csv(assembly_key, file.path(input_root, "private_genome_label_key.csv"))

message("Draft initialized: current 532-genome manifest validated. RQ09-RQ11 implementation follows.")

# -----------------------------------------------------------------------------
# Provider STs: the active provider table is allowed only after exact subsetting
# to the current manifest and path/assembler validation. Its extra 24 rows are
# never carried into an estimand.
# -----------------------------------------------------------------------------

provider_all <- read_csv(required_inputs[["provider_mlst"]], show_col_types = FALSE) %>%
  mutate(
    Assembly_ID = as.character(Assembly_ID),
    ST = as.character(ST),
    full_path = norm_path(full_path)
  )

provider_required <- c(
  "Assembly_ID", "ST", "ST_source", "ST_provider", "provider_ST_called",
  "provider_file", "provider_assembler", "Assembler", "full_path"
)
if (length(setdiff(provider_required, names(provider_all)))) {
  stop("Provider MLST table lacks: ", paste(setdiff(provider_required, names(provider_all)), collapse = ", "))
}

provider_current <- manifest %>%
  select(
    Assembly_ID, Isolate_ID, Participant_id, tp_lab, Event_type,
    Collection_Date_parsed, manifest_path = full_path
  ) %>%
  left_join(
    provider_all %>%
      select(
        Assembly_ID, ST, ST_source, ST_provider, provider_ST_called,
        provider_file, provider_assembler, mlst_assembler = Assembler,
        mlst_path = full_path
      ),
    by = "Assembly_ID"
  ) %>%
  mutate(
    path_exact = !is.na(mlst_path) & manifest_path == mlst_path,
    provider_provenance_valid =
      path_exact &
      tolower(mlst_assembler) == "longcycler" &
      ST_source == "provider_qc95" &
      provider_ST_called %in% TRUE &
      !is.na(ST) & ST != ""
  )

provider_checks <- tibble(
  check = c(
    "current_manifest_rows_joined_once", "all_paths_exact",
    "no_non_longcycler_active_rows", "provider_qc95_calls"
  ),
  passed = c(
    nrow(provider_current) == 532L,
    all(provider_current$path_exact),
    all(tolower(provider_current$mlst_assembler) == "longcycler"),
    sum(provider_current$provider_provenance_valid) == 509L
  ),
  observed = c(
    as.character(nrow(provider_current)),
    as.character(sum(provider_current$path_exact)),
    paste(sort(unique(provider_current$mlst_assembler)), collapse = ";"),
    as.character(sum(provider_current$provider_provenance_valid))
  )
)
atomic_write_csv(provider_checks, file.path(input_root, "provider_st_current_checks.csv"))
if (!all(provider_checks$passed)) stop("Current provider-ST provenance checks failed.")

provider_st <- provider_current %>%
  filter(provider_provenance_valid) %>%
  select(
    Assembly_ID, Isolate_ID, Participant_id, tp_lab, Event_type,
    Collection_Date_parsed, full_path = manifest_path, ST,
    ST_source, ST_provider, provider_file, provider_assembler
  )
atomic_write_csv(provider_st, file.path(input_root, "provider_st_current_509.csv"))

parallel_lapply <- function(x, fun, mc.cores = workers) {
  if (.Platform$OS.type == "windows" || mc.cores <= 1L) return(lapply(x, fun))
  parallel::mclapply(x, fun, mc.cores = mc.cores, mc.preschedule = TRUE)
}

hash_paths <- function(paths) {
  paths <- unique(norm_path(paths))
  stopifnot(all(file.exists(paths)))
  values <- unlist(parallel_lapply(paths, sha256_file), use.names = FALSE)
  tibble(full_path = paths, fasta_sha256 = values, fasta_size_bytes = file.size(paths))
}

manifest_hashes <- hash_paths(manifest$full_path)
current_manifest_hashed <- manifest %>%
  left_join(manifest_hashes, by = "full_path")
atomic_write_csv(
  current_manifest_hashed %>%
    select(
      Assembly_ID, Isolate_ID, Participant_id, tp_lab, Event_type,
      Collection_Date, Assembler, full_path, fasta_sha256, fasta_size_bytes
    ),
  file.path(input_root, "current_532_fasta_manifest_sha256.csv")
)

tool_version <- function(tool, args = "--version") {
  exe <- Sys.which(tool)
  if (exe == "") return(NA_character_)
  out <- suppressWarnings(system2(exe, args, stdout = TRUE, stderr = TRUE))
  paste(out, collapse = " ")
}

run_system2 <- function(command, args, stdout, stderr) {
  status <- suppressWarnings(system2(command, args, stdout = stdout, stderr = stderr))
  if (is.null(status)) 0L else as.integer(status)
}

valid_sidecar <- function(sidecar_path, expected, output_path) {
  if (!file.exists(sidecar_path) || !file.exists(output_path)) return(FALSE)
  side <- tryCatch(jsonlite::read_json(sidecar_path, simplifyVector = TRUE), error = function(e) NULL)
  if (is.null(side)) return(FALSE)
  isTRUE(identical(as.character(side$signature), as.character(expected))) &&
    isTRUE(side$exit_code == 0L) &&
    isTRUE(identical(as.character(side$output_sha256), sha256_file(output_path)))
}

# -----------------------------------------------------------------------------
# RQ09. Are same-provider-ST genomes more similar within than between residents?
# -----------------------------------------------------------------------------

run_rq09 <- function() {
  rq <- "RQ09"
  out_dir <- rq_dirs[[rq]]
  old_bootstrap <- file.path(out_dir, "resident_block_bootstrap_distribution.csv")
  if (file.exists(old_bootstrap)) unlink(old_bootstrap)
  mash <- Sys.which("mash")
  if (mash == "") {
    write_status(rq, "blocked", "mash_not_available",
                 "No Mash calculation or inferential output was fabricated.")
    atomic_write_csv(tibble(tool = "mash", available = FALSE), file.path(out_dir, "tool_check.csv"))
    return(invisible(FALSE))
  }

  mash_version <- tool_version("mash")
  eligible_st <- provider_st %>%
    count(ST, name = "n_genomes") %>%
    left_join(provider_st %>% distinct(ST, Participant_id) %>% count(ST, name = "n_residents"), by = "ST") %>%
    filter(n_genomes >= 5L, n_residents >= 2L) %>%
    arrange(desc(n_genomes), ST)
  atomic_write_csv(eligible_st, file.path(out_dir, "eligible_provider_st.csv"))
  if (!nrow(eligible_st)) {
    write_status(rq, "blocked", "no_provider_st_meets_size_rule")
    return(invisible(FALSE))
  }

  mash_input <- current_manifest_hashed %>% arrange(Assembly_ID)
  signature <- digest::digest(
    paste(
      "rq09_mash_v1", mash_version, "k=21", "s=1000",
      paste(mash_input$Assembly_ID, mash_input$fasta_sha256, sep = "=", collapse = "\n"),
      sep = "\n"
    ),
    algo = "sha256", serialize = FALSE
  )
  mash_dir <- file.path(input_root, "rq09_mash", signature)
  dir.create(mash_dir, recursive = TRUE, showWarnings = FALSE)
  list_path <- file.path(mash_dir, "fasta_paths.txt")
  sketch_prefix <- file.path(mash_dir, "current_532")
  sketch_path <- paste0(sketch_prefix, ".msh")
  dist_path <- file.path(mash_dir, "all_vs_all_dist.tsv")
  sidecar_path <- file.path(mash_dir, "provenance.json")
  atomic_write_lines(mash_input$full_path, list_path)

  cache_ok <- valid_sidecar(sidecar_path, signature, dist_path) && file.exists(sketch_path)
  if (!cache_ok) {
    sketch_tmp_prefix <- file.path(mash_dir, paste0("current_532.tmp.", Sys.getpid()))
    sketch_tmp <- paste0(sketch_tmp_prefix, ".msh")
    sketch_err <- file.path(mash_dir, "sketch.stderr.txt")
    sketch_status <- run_system2(
      mash,
      c("sketch", "-p", as.character(workers), "-o", sketch_tmp_prefix, "-l", list_path),
      stdout = file.path(mash_dir, "sketch.stdout.txt"), stderr = sketch_err
    )
    if (sketch_status != 0L || !file.exists(sketch_tmp)) {
      write_status(rq, "blocked", "mash_sketch_failed",
                   paste(readLines(sketch_err, warn = FALSE), collapse = " | "))
      return(invisible(FALSE))
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
      write_status(rq, "blocked", "mash_dist_failed",
                   paste(readLines(dist_err, warn = FALSE), collapse = " | "))
      return(invisible(FALSE))
    }
    file.copy(dist_tmp, dist_path, overwrite = TRUE)
    unlink(dist_tmp)
    jsonlite::write_json(list(
      signature = signature,
      created_utc = run_timestamp,
      mash_executable = mash,
      mash_version = mash_version,
      parameters = list(k = 21L, sketch_size = 1000L),
      n_genomes = nrow(mash_input),
      manifest_sha256 = sha256_file(required_inputs[["manifest"]]),
      output_sha256 = sha256_file(dist_path),
      exit_code = 0L
    ), sidecar_path, auto_unbox = TRUE, pretty = TRUE)
  }

  mash_raw <- read_tsv(
    dist_path,
    col_names = c("reference", "query", "mash_distance", "p_value", "shared_hashes"),
    col_types = cols(
      reference = col_character(), query = col_character(), mash_distance = col_double(),
      p_value = col_double(), shared_hashes = col_character()
    ),
    show_col_types = FALSE, progress = FALSE
  ) %>%
    mutate(reference = norm_path(reference), query = norm_path(query))

  path_index <- setNames(seq_len(nrow(mash_input)), mash_input$full_path)
  ref_i <- unname(path_index[mash_raw$reference])
  qry_i <- unname(path_index[mash_raw$query])
  if (anyNA(ref_i) || anyNA(qry_i)) {
    write_status(rq, "blocked", "mash_output_path_mismatch",
                 sprintf("Unmatched reference/query rows: %d/%d", sum(is.na(ref_i)), sum(is.na(qry_i))))
    return(invisible(FALSE))
  }

  keep <- ref_i < qry_i
  ref_i <- ref_i[keep]
  qry_i <- qry_i[keep]
  mash_unique <- mash_raw[keep, ] %>%
    transmute(
      Assembly_ID_A = mash_input$Assembly_ID[ref_i],
      Assembly_ID_B = mash_input$Assembly_ID[qry_i],
      mash_distance, p_value, shared_hashes
    )

  st_meta <- provider_st %>%
    select(Assembly_ID, Participant_id, Collection_Date_parsed, ST)
  pairs <- mash_unique %>%
    inner_join(st_meta %>% rename_with(~ paste0(.x, "_A"), -Assembly_ID) %>%
                 rename(Assembly_ID_A = Assembly_ID), by = "Assembly_ID_A") %>%
    inner_join(st_meta %>% rename_with(~ paste0(.x, "_B"), -Assembly_ID) %>%
                 rename(Assembly_ID_B = Assembly_ID), by = "Assembly_ID_B") %>%
    filter(ST_A == ST_B, ST_A %in% eligible_st$ST) %>%
    transmute(
      ST = ST_A,
      Assembly_ID_A, Assembly_ID_B,
      Participant_id_A, Participant_id_B,
      date_A = Collection_Date_parsed_A,
      date_B = Collection_Date_parsed_B,
      days_between = abs(as.numeric(date_B - date_A)),
      same_resident = Participant_id_A == Participant_id_B,
      mash_distance, p_value, shared_hashes,
      resident_pair = if_else(
        Participant_id_A <= Participant_id_B,
        paste(Participant_id_A, Participant_id_B, sep = "||"),
        paste(Participant_id_B, Participant_id_A, sep = "||")
      )
    )

  estimate_effect <- function(pair_data) {
    st_effects <- pair_data %>%
      group_by(ST, same_resident) %>%
      summarise(n_pairs = n(), median_mash = median(mash_distance), .groups = "drop") %>%
      pivot_wider(
        names_from = same_resident,
        values_from = c(n_pairs, median_mash),
        names_glue = "{.value}_{ifelse(same_resident, 'within', 'between')}"
      ) %>%
      filter(!is.na(median_mash_within), !is.na(median_mash_between)) %>%
      mutate(within_minus_between = median_mash_within - median_mash_between)
    list(
      by_st = st_effects,
      estimate = if (nrow(st_effects)) mean(st_effects$within_minus_between) else NA_real_,
      n_st = nrow(st_effects)
    )
  }

  primary_pairs <- pairs %>% filter(!is.na(days_between), days_between <= 180)
  primary <- estimate_effect(primary_pairs)
  primary_pairs_display <- primary_pairs %>%
    left_join(assembly_key %>% rename(Assembly_ID_A = Assembly_ID, Genome_Label_A = Genome_Label), by = "Assembly_ID_A") %>%
    left_join(assembly_key %>% rename(Assembly_ID_B = Assembly_ID, Genome_Label_B = Genome_Label), by = "Assembly_ID_B") %>%
    left_join(resident_key %>% rename(Participant_id_A = Participant_id, Resident_Label_A = Resident_Label), by = "Participant_id_A") %>%
    left_join(resident_key %>% rename(Participant_id_B = Participant_id, Resident_Label_B = Resident_Label), by = "Participant_id_B") %>%
    transmute(
      Pair_Label = sprintf("RQ09_Pair_%05d", row_number()),
      ST, Genome_Label_A, Genome_Label_B, Resident_Label_A, Resident_Label_B,
      days_between, same_resident,
      mash_distance, p_value, shared_hashes
    )
  atomic_write_csv(primary_pairs_display, file.path(out_dir, "same_st_pairs_within_180_days.csv"))
  atomic_write_csv(primary$by_st, file.path(out_dir, "st_equal_weight_effects_180_days.csv"))
  if (!is.finite(primary$estimate) || primary$n_st == 0L) {
    write_status(rq, "blocked", "no_st_has_both_within_and_between_resident_pairs_within_180_days")
    return(invisible(FALSE))
  }

  effect_st <- primary$by_st$ST
  pair_split <- split(filter(primary_pairs, ST %in% effect_st), ~ ST)
  genome_split <- split(filter(provider_st, ST %in% effect_st), ~ ST)
  set.seed(seed + 9L)
  permutation_one <- function() {
    values <- map_dbl(effect_st, function(st) {
      p <- pair_split[[st]]
      g <- genome_split[[st]]
      labels <- sample(g$Participant_id, length(g$Participant_id), replace = FALSE)
      names(labels) <- g$Assembly_ID
      within <- labels[p$Assembly_ID_A] == labels[p$Assembly_ID_B]
      if (!any(within) || all(within)) return(NA_real_)
      median(p$mash_distance[within]) - median(p$mash_distance[!within])
    })
    if (all(is.na(values))) NA_real_ else mean(values, na.rm = TRUE)
  }
  perm_values <- replicate(n_perm, permutation_one())
  perm_tbl <- tibble(iteration = seq_len(n_perm), equal_weight_difference = perm_values)
  atomic_write_csv(perm_tbl, file.path(out_dir, "resident_label_permutation_distribution.csv"))
  permutation_p <- (1 + sum(abs(perm_values) >= abs(primary$estimate), na.rm = TRUE)) /
    (1 + sum(is.finite(perm_values)))

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
  bootstrap_residents <- sort(unique(c(primary_pairs$Participant_id_A, primary_pairs$Participant_id_B)))
  set.seed(seed + 90L)
  bootstrap_one <- function(iteration) {
    sampled <- sample(bootstrap_residents, length(bootstrap_residents), replace = TRUE)
    multiplicity <- table(factor(sampled, levels = bootstrap_residents))
    diffs <- map_dbl(effect_st, function(st) {
      p <- pair_split[[st]]
      m_a <- as.numeric(multiplicity[match(p$Participant_id_A, bootstrap_residents)])
      m_b <- as.numeric(multiplicity[match(p$Participant_id_B, bootstrap_residents)])
      pair_weight <- ifelse(p$same_resident, m_a, m_a * m_b)
      within_median <- weighted_median(p$mash_distance[p$same_resident], pair_weight[p$same_resident])
      between_median <- weighted_median(p$mash_distance[!p$same_resident], pair_weight[!p$same_resident])
      within_median - between_median
    })
    valid_st <- sum(is.finite(diffs))
    tibble(
      iteration = iteration,
      valid_st = valid_st,
      all_primary_st_valid = valid_st == length(effect_st),
      equal_weight_difference = if (valid_st == length(effect_st)) mean(diffs) else NA_real_
    )
  }
  bootstrap_tbl <- map_dfr(seq_len(n_boot), bootstrap_one)
  boot_values <- bootstrap_tbl$equal_weight_difference
  atomic_write_csv(
    bootstrap_tbl,
    file.path(out_dir, "resident_frequency_bootstrap_distribution.csv")
  )

  dominant_st <- eligible_st$ST[[which.max(eligible_st$n_genomes)]]
  sensitivity_specs <- tribble(
    ~analysis, ~window_days, ~excluded_st,
    "primary_180_days", 180, NA_character_,
    "window_90_days", 90, NA_character_,
    "window_365_days", 365, NA_character_,
    "leave_dominant_st_out_180_days", 180, dominant_st
  )
  sensitivity <- pmap_dfr(sensitivity_specs, function(analysis, window_days, excluded_st) {
    dat <- pairs %>% filter(!is.na(days_between), days_between <= window_days)
    if (!is.na(excluded_st)) dat <- filter(dat, ST != excluded_st)
    est <- estimate_effect(dat)
    tibble(
      analysis = analysis,
      window_days = window_days,
      excluded_st = excluded_st,
      n_pairs = nrow(dat),
      n_within_resident_pairs = sum(dat$same_resident),
      n_between_resident_pairs = sum(!dat$same_resident),
      n_st_equal_weighted = est$n_st,
      equal_weight_within_minus_between_mash = est$estimate
    )
  })
  atomic_write_csv(sensitivity, file.path(out_dir, "sensitivity_analyses.csv"))

  primary_summary <- tibble(
    estimand = "Mean across provider STs of (within-resident median Mash distance - between-resident median Mash distance)",
    n_current_genomes = nrow(manifest),
    n_provider_st_genomes = nrow(provider_st),
    n_eligible_st = nrow(eligible_st),
    n_st_in_estimand = primary$n_st,
    n_pairs_within_180_days = nrow(primary_pairs),
    observed_difference = primary$estimate,
    bootstrap_ci_low = unname(quantile(boot_values, 0.025, na.rm = TRUE)),
    bootstrap_ci_high = unname(quantile(boot_values, 0.975, na.rm = TRUE)),
    bootstrap_valid_replicates = sum(is.finite(boot_values)),
    bootstrap_failed_replicates = sum(!is.finite(boot_values)),
    permutation_p_two_sided = permutation_p,
    permutations = n_perm,
    bootstrap_replicates = n_boot
  )
  atomic_write_csv(primary_summary, file.path(out_dir, "primary_result.csv"))
  atomic_write_csv(tibble(
    input = c("analysis_manifest", "provider_mlst", "mash_distance_cache"),
    path = c(required_inputs[["manifest"]], required_inputs[["provider_mlst"]], dist_path),
    sha256 = c(sha256_file(required_inputs[["manifest"]]),
               sha256_file(required_inputs[["provider_mlst"]]), sha256_file(dist_path)),
    role = c("current_532_genomes", "exact_Assembly_ID_subset_provider_ST", "fresh_SHA_bound_all_vs_all_Mash")
  ), file.path(out_dir, "input_provenance.csv"))
  atomic_write_lines(c(
    "# RQ09 — same-provider-ST Mash distance",
    "",
    "Status: complete.",
    "",
    "The denominator is the current 532-genome Longcycler manifest. Only exact, path-validated",
    "provider_qc95 ST calls are used. STs require at least five genomes and two residents.",
    "The primary window is 180 days. The estimand gives each ST equal weight. Negative values",
    "mean lower Mash distance within residents. The null is generated by 10,000 resident-label",
    "permutations within ST. Uncertainty resamples residents as clusters; within-resident pair",
    "weights are the resident draw multiplicity and between-resident weights are the product of",
    "the two multiplicities. Failed replicates lacking both pair types for every primary ST are explicit.",
    "Sensitivity analyses use 90 and 365 days and remove the most frequent eligible ST."
  ), file.path(out_dir, "README.md"))
  write_status(rq, "complete", "analysis_completed", sprintf("Mash %s; %d current genomes", mash_version, nrow(manifest)))
  invisible(TRUE)
}

# -----------------------------------------------------------------------------
# RQ10. How often do provider STs turn over between adjacent current-genome
# episodes, and does turnover vary with time between samples?
# -----------------------------------------------------------------------------

run_rq10 <- function() {
  rq <- "RQ10"
  out_dir <- rq_dirs[[rq]]
  transitions <- read_csv(required_inputs[["longcycler_transitions"]], show_col_types = FALSE) %>%
    mutate(
      Participant_id = as.character(Participant_id),
      tp_from = as.character(tp_from),
      tp_to = as.character(tp_to),
      fasta_path_from = norm_path(fasta_path_from),
      fasta_path_to = norm_path(fasta_path_to)
    )

  current_paths <- manifest$full_path
  transition_checks <- tibble(
    check = c(
      "current_longcycler_adjacency_371", "all_from_paths_current",
      "all_to_paths_current", "all_dnadiff_cache_sha_bound",
      "all_dnadiff_sidecars_exist", "all_snp_values_present"
    ),
    passed = c(
      nrow(transitions) == 371L,
      all(transitions$fasta_path_from %in% current_paths),
      all(transitions$fasta_path_to %in% current_paths),
      all(transitions$dnadiff_cache_status %in% c("fresh", "generated")),
      all(file.exists(transitions$dnadiff_sidecar_path)),
      all(!is.na(transitions$TotalSNPs))
    ),
    observed = c(
      as.character(nrow(transitions)),
      as.character(sum(transitions$fasta_path_from %in% current_paths)),
      as.character(sum(transitions$fasta_path_to %in% current_paths)),
      paste(names(table(transitions$dnadiff_cache_status)), as.integer(table(transitions$dnadiff_cache_status)), collapse = ";"),
      as.character(sum(file.exists(transitions$dnadiff_sidecar_path))),
      as.character(sum(!is.na(transitions$TotalSNPs)))
    )
  )
  atomic_write_csv(transition_checks, file.path(out_dir, "input_checks.csv"))
  if (!all(transition_checks$passed)) {
    write_status(rq, "blocked", "current_longcycler_transition_provenance_failed")
    return(invisible(FALSE))
  }

  st_key <- provider_st %>%
    select(Participant_id, tp_lab, ST, Assembly_ID, Event_type) %>%
    distinct()
  if (nrow(st_key) != n_distinct(paste(st_key$Participant_id, st_key$tp_lab, sep = "||"))) {
    write_status(rq, "blocked", "provider_st_key_not_unique")
    return(invisible(FALSE))
  }

  adjacent_all <- transitions %>%
    left_join(
      st_key %>% rename(tp_from = tp_lab, ST_from_provider = ST,
                        Assembly_ID_from_provider = Assembly_ID, event_type_from = Event_type),
      by = c("Participant_id", "tp_from")
    ) %>%
    left_join(
      st_key %>% rename(tp_to = tp_lab, ST_to_provider = ST,
                        Assembly_ID_to_provider = Assembly_ID, event_type_to = Event_type),
      by = c("Participant_id", "tp_to")
    ) %>%
    mutate(
      provider_ST_from_available = !is.na(ST_from_provider),
      provider_ST_to_available = !is.na(ST_to_provider),
      provider_ST_both_available = provider_ST_from_available & provider_ST_to_available,
      ST_changed = ST_from_provider != ST_to_provider,
      days_between_samples = as.numeric(days_between_samples),
      routine_both = event_type_from == "Routine" & event_type_to == "Routine",
      within_365_days = days_between_samples <= 365,
      ST_relation = if_else(ST_changed, "ST_changed", "ST_same"),
      SNP_relation = case_when(
        is.na(TotalSNPs) ~ "SNP_missing",
        TotalSNPs <= 25 ~ "SNP_le_25",
        TRUE ~ "SNP_gt_25"
      )
    )
  provider_attrition <- bind_rows(
    tibble(stage = "current_Longcycler_adjacent_transitions", n_pairs = nrow(adjacent_all)),
    tibble(stage = "provider_ST_available_at_from_endpoint", n_pairs = sum(adjacent_all$provider_ST_from_available)),
    tibble(stage = "provider_ST_available_at_to_endpoint", n_pairs = sum(adjacent_all$provider_ST_to_available)),
    tibble(stage = "provider_ST_available_at_both_endpoints", n_pairs = sum(adjacent_all$provider_ST_both_available)),
    tibble(stage = "excluded_missing_provider_ST_at_either_endpoint", n_pairs = sum(!adjacent_all$provider_ST_both_available))
  )
  atomic_write_csv(provider_attrition, file.path(out_dir, "provider_st_endpoint_attrition.csv"))
  adjacent <- adjacent_all %>% filter(provider_ST_both_available)

  if (!nrow(adjacent) || nrow(adjacent) > 355L) {
    write_status(rq, "blocked", "unexpected_provider_only_adjacency_denominator",
                 sprintf("Observed %d; required 1-355.", nrow(adjacent)))
    return(invisible(FALSE))
  }
  adjacent_display <- adjacent %>%
    left_join(resident_key, by = "Participant_id") %>%
    left_join(
      episode_key %>% rename(tp_from = tp_lab, Episode_Label_From = Episode_Label),
      by = c("Participant_id", "tp_from")
    ) %>%
    left_join(
      episode_key %>% rename(tp_to = tp_lab, Episode_Label_To = Episode_Label),
      by = c("Participant_id", "tp_to")
    ) %>%
    transmute(
      Transition_Label = sprintf("RQ10_Transition_%03d", row_number()),
      Resident_Label, Episode_Label_From, Episode_Label_To,
      ST_From = ST_from_provider, ST_To = ST_to_provider, ST_changed,
      event_type_from, event_type_to, routine_both,
      days_between_samples, within_365_days,
      TotalSNPs, SNP_relation, ST_relation,
      AvgIdentity, MashDistance, Classification,
      dnadiff_cache_status, dnadiff_version
    )
  atomic_write_csv(adjacent_display, file.path(out_dir, "provider_only_adjacent_st_pairs.csv"))

  fit_spline <- function(data) {
    data <- filter(data, is.finite(days_between_samples), !is.na(ST_changed))
    if (nrow(data) < 20L || n_distinct(data$ST_changed) < 2L || n_distinct(data$days_between_samples) < 5L) {
      return(NULL)
    }
    suppressWarnings(tryCatch(
      glm(ST_changed ~ splines::ns(days_between_samples, df = 3),
          data = data, family = binomial()),
      error = function(e) NULL
    ))
  }

  model_coefficients <- function(model, analysis, n, residents) {
    if (is.null(model)) {
      return(tibble(
        analysis = analysis, n_pairs = n, n_residents = residents,
        term = NA_character_, estimate_log_odds = NA_real_, standard_error = NA_real_,
        z_value = NA_real_, p_value = NA_real_, model_status = "not_estimable"
      ))
    }
    sm <- coef(summary(model))
    tibble(
      analysis = analysis,
      n_pairs = n,
      n_residents = residents,
      term = rownames(sm),
      estimate_log_odds = sm[, "Estimate"],
      standard_error = sm[, "Std. Error"],
      z_value = sm[, "z value"],
      p_value = sm[, "Pr(>|z|)"],
      model_status = if_else(model$converged, "converged", "not_converged")
    )
  }

  analyses <- list(
    primary = adjacent,
    routine_both = filter(adjacent, routine_both),
    within_365_days = filter(adjacent, within_365_days)
  )
  models <- map(analyses, fit_spline)
  coefficient_table <- imap_dfr(models, function(model, nm) {
    dat <- analyses[[nm]]
    model_coefficients(model, nm, nrow(dat), n_distinct(dat$Participant_id))
  })
  atomic_write_csv(coefficient_table, file.path(out_dir, "spline_model_coefficients.csv"))

  sensitivity_summary <- imap_dfr(analyses, function(dat, nm) {
    tibble(
      analysis = nm,
      n_pairs = nrow(dat),
      n_residents = n_distinct(dat$Participant_id),
      n_st_changed = sum(dat$ST_changed),
      turnover_proportion = mean(dat$ST_changed),
      median_days = median(dat$days_between_samples),
      max_days = max(dat$days_between_samples),
      model_estimable = !is.null(models[[nm]])
    )
  })
  atomic_write_csv(sensitivity_summary, file.path(out_dir, "turnover_sensitivity_summary.csv"))

  primary_model <- models[["primary"]]
  if (is.null(primary_model)) {
    write_status(rq, "blocked", "primary_time_model_not_estimable")
    return(invisible(FALSE))
  }

  time_grid <- tibble(days_between_samples = seq(
    min(adjacent$days_between_samples), max(adjacent$days_between_samples), length.out = 61L
  ))
  primary_curve <- time_grid %>%
    mutate(predicted_turnover_probability = as.numeric(predict(primary_model, newdata = time_grid, type = "response")))

  residents <- unique(adjacent$Participant_id)
  set.seed(seed + 10L)
  bootstrap_one <- function(iteration) {
    sampled <- sample(residents, length(residents), replace = TRUE)
    boot_data <- map2_dfr(sampled, seq_along(sampled), function(pid, draw) {
      filter(adjacent, Participant_id == pid) %>% mutate(bootstrap_cluster = draw)
    })
    model <- fit_spline(boot_data)
    list(
      iteration = iteration,
      proportion = mean(boot_data$ST_changed),
      curve = if (is.null(model)) rep(NA_real_, nrow(time_grid)) else
        as.numeric(predict(model, newdata = time_grid, type = "response"))
    )
  }
  bootstrap_results <- parallel_lapply(seq_len(n_boot), bootstrap_one)
  bootstrap_proportions <- tibble(
    iteration = seq_len(n_boot),
    turnover_proportion = map_dbl(bootstrap_results, "proportion")
  )
  curve_matrix <- do.call(rbind, map(bootstrap_results, "curve"))
  curve_summary <- primary_curve %>%
    mutate(
      bootstrap_median = apply(curve_matrix, 2, median, na.rm = TRUE),
      bootstrap_ci_low = apply(curve_matrix, 2, quantile, probs = 0.025, na.rm = TRUE),
      bootstrap_ci_high = apply(curve_matrix, 2, quantile, probs = 0.975, na.rm = TRUE),
      valid_bootstrap_models = colSums(is.finite(curve_matrix))
    )
  valid_proportion_boot <- is.finite(bootstrap_proportions$turnover_proportion)
  primary_turnover <- tibble(
    estimand = "Proportion of adjacent current-genome pairs with different provider-qc95 STs",
    complete_case_pairs = nrow(adjacent),
    residents = n_distinct(adjacent$Participant_id),
    st_changes = sum(adjacent$ST_changed),
    observed_turnover_proportion = mean(adjacent$ST_changed),
    bootstrap_ci_low = quantile(bootstrap_proportions$turnover_proportion, 0.025, na.rm = TRUE),
    bootstrap_ci_high = quantile(bootstrap_proportions$turnover_proportion, 0.975, na.rm = TRUE),
    requested_bootstrap_replicates = n_boot,
    valid_bootstrap_replicates = sum(valid_proportion_boot),
    failed_bootstrap_replicates = sum(!valid_proportion_boot),
    point_spline_model_converged = isTRUE(primary_model$converged)
  )
  atomic_write_csv(bootstrap_proportions, file.path(out_dir, "resident_bootstrap_turnover_proportion.csv"))
  atomic_write_csv(curve_summary, file.path(out_dir, "resident_bootstrap_time_curve.csv"))
  atomic_write_csv(primary_turnover, file.path(out_dir, "primary_turnover_estimate.csv"))

  st_snp <- adjacent %>%
    count(ST_relation, SNP_relation, name = "n_transitions") %>%
    group_by(ST_relation) %>%
    mutate(row_proportion = n_transitions / sum(n_transitions)) %>%
    ungroup()
  atomic_write_csv(st_snp, file.path(out_dir, "st_by_snp_cross_tab.csv"))

  atomic_write_csv(tibble(
    input = c("analysis_manifest", "provider_mlst", "longcycler_transitions"),
    path = unname(required_inputs[c("manifest", "provider_mlst", "longcycler_transitions")]),
    sha256 = map_chr(unname(required_inputs[c("manifest", "provider_mlst", "longcycler_transitions")]), sha256_file),
    role = c("current_532_genomes", "exact_provider_qc95_calls", "current_SHA_provenanced_adjacent_dnadiff_pairs")
  ), file.path(out_dir, "input_provenance.csv"))
  atomic_write_lines(c(
    "# RQ10 — adjacent provider-ST turnover",
    "",
    "Status: complete.",
    "",
    "Adjacency is inherited only from the verified 371-transition Longcycler table. The",
    "analysis then requires provider_qc95 ST calls at both endpoints; missing ST calls do",
    "not cause visits to be bridged. The binary outcome is exact provider-ST change.",
    "A logistic model uses a three-degree-of-freedom natural spline for days between samples.",
    "Uncertainty for the turnover proportion and time curve resamples residents as clusters.",
    "Sensitivities require two routine endpoints or a gap of at most 365 days. The ST-by-SNP",
    "table is descriptive and uses the operational 25-SNP boundary."
  ), file.path(out_dir, "README.md"))
  write_status(rq, "complete", "analysis_completed",
               sprintf("%d provider-only adjacent pairs from %d residents", nrow(adjacent), n_distinct(adjacent$Participant_id)))
  invisible(TRUE)
}

# -----------------------------------------------------------------------------
# RQ11. Paired assembler measurement agreement. This is a methods/sensitivity
# question, not an AMR analysis. Calls are rerun with identical settings in an
# isolated cache whose key contains the exact FASTA SHA-256 and tool/database
# signature. Existing basename-only caches are never trusted.
# -----------------------------------------------------------------------------

write_json_atomic <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  tmp <- tempfile(pattern = paste0(".", basename(path), "."), tmpdir = dirname(path))
  jsonlite::write_json(x, tmp, auto_unbox = TRUE, pretty = TRUE, null = "null", na = "null")
  if (!file.rename(tmp, path)) {
    file.copy(tmp, path, overwrite = TRUE)
    unlink(tmp)
  }
  invisible(path)
}

abricate_database_table <- function(executable) {
  raw <- suppressWarnings(system2(executable, "--list", stdout = TRUE, stderr = TRUE))
  tryCatch(
    read_tsv(I(paste(raw, collapse = "\n")), show_col_types = FALSE, progress = FALSE),
    error = function(e) tibble()
  )
}

count_abricate_rows <- function(path) {
  if (!file.exists(path) || file.size(path) == 0L) return(NA_integer_)
  lines <- readLines(path, warn = FALSE)
  sum(nzchar(lines) & !startsWith(lines, "#"))
}

count_mlst_rows <- function(path) {
  if (!file.exists(path) || file.size(path) == 0L) return(0L)
  sum(nzchar(readLines(path, warn = FALSE)))
}

run_sha_bound_call <- function(job, tools, signatures, cache_root) {
  layer <- job$layer[[1]]
  input_sha <- job$fasta_sha256[[1]]
  input_path <- job$full_path[[1]]
  signature <- signatures[[layer]]
  exe <- tools[[layer]]
  layer_dir <- file.path(cache_root, layer, digest::digest(signature, algo = "sha256", serialize = FALSE))
  dir.create(layer_dir, recursive = TRUE, showWarnings = FALSE)
  extension <- if (layer == "mlst_ecoli") "csv" else "tsv"
  output_path <- file.path(layer_dir, paste0(input_sha, ".", extension))
  sidecar_path <- file.path(layer_dir, paste0(input_sha, ".provenance.json"))
  stderr_path <- file.path(layer_dir, paste0(input_sha, ".stderr.txt"))

  cached <- FALSE
  call_status <- "error"
  exit_code <- NA_integer_
  call_rows <- NA_integer_
  if (exe != "" && file.exists(sidecar_path) && file.exists(output_path)) {
    side <- tryCatch(read_json(sidecar_path, simplifyVector = TRUE), error = function(e) NULL)
    cached <- !is.null(side) &&
      identical(as.character(side$signature), as.character(signature)) &&
      identical(as.character(side$input_sha256), as.character(input_sha)) &&
      isTRUE(side$exit_code == 0L) &&
      identical(as.character(side$output_sha256), sha256_file(output_path))
    if (cached) {
      call_status <- as.character(side$call_status)
      exit_code <- as.integer(side$exit_code)
      call_rows <- as.integer(side$call_rows)
    }
  }

  if (!cached) {
    if (exe == "") {
      call_status <- "tool_missing"
      exit_code <- 127L
      call_rows <- NA_integer_
    } else {
      tmp_output <- tempfile(pattern = paste0(".", input_sha, "."), tmpdir = layer_dir)
      tmp_stderr <- tempfile(pattern = paste0(".", input_sha, ".stderr."), tmpdir = layer_dir)
      args <- switch(
        layer,
        mlst_ecoli = c("--quiet", "--scheme", "ecoli", "--csv", input_path),
        vfdb = c("--quiet", "--db", "vfdb", "--mincov", "80", "--minid", "80", input_path),
        plasmidfinder = c("--quiet", "--db", "plasmidfinder", "--mincov", "80", "--minid", "80", input_path),
        stop("Unknown call layer: ", layer)
      )
      exit_code <- run_system2(exe, args, stdout = tmp_output, stderr = tmp_stderr)
      call_rows <- if (layer == "mlst_ecoli") count_mlst_rows(tmp_output) else count_abricate_rows(tmp_output)
      if (exit_code != 0L || is.na(call_rows)) {
        call_status <- "error"
      } else if (layer == "mlst_ecoli" && call_rows == 0L) {
        call_status <- "no_call"
      } else if (layer != "mlst_ecoli" && call_rows == 0L) {
        call_status <- "zero_hit"
      } else {
        call_status <- "success"
      }
      if (file.exists(tmp_output)) file.copy(tmp_output, output_path, overwrite = TRUE)
      if (file.exists(tmp_stderr)) file.copy(tmp_stderr, stderr_path, overwrite = TRUE)
      unlink(c(tmp_output, tmp_stderr))
      output_sha <- if (file.exists(output_path)) sha256_file(output_path) else NA_character_
      write_json_atomic(list(
        cache_schema = "rq11_sha_bound_call_v1",
        signature = signature,
        layer = layer,
        input_path = input_path,
        input_sha256 = input_sha,
        output_path = output_path,
        output_sha256 = output_sha,
        call_status = call_status,
        call_rows = call_rows,
        exit_code = exit_code,
        created_utc = format(Sys.time(), tz = "UTC", usetz = TRUE)
      ), sidecar_path)
    }
  }

  tibble(
    layer = layer,
    input_sha256 = input_sha,
    representative_input_path = input_path,
    tool_signature = signature,
    output_path = output_path,
    sidecar_path = sidecar_path,
    call_status = call_status,
    call_rows = call_rows,
    exit_code = exit_code,
    cache_reused = cached
  )
}

parse_mlst_st <- function(path, status) {
  if (status != "success" || !file.exists(path)) return(NA_character_)
  x <- tryCatch(
    utils::read.csv(path, header = FALSE, colClasses = "character", check.names = FALSE),
    error = function(e) NULL
  )
  if (is.null(x) || nrow(x) != 1L || ncol(x) < 3L) return(NA_character_)
  as.character(x[[3]][[1]])
}

parse_abricate_features <- function(path, status) {
  if (!status %in% c("success", "zero_hit") || !file.exists(path)) return(NULL)
  if (status == "zero_hit") return(character())
  x <- tryCatch(
    read_tsv(path, show_col_types = FALSE, progress = FALSE, name_repair = "minimal"),
    error = function(e) NULL
  )
  if (is.null(x) || !"GENE" %in% names(x)) return(NULL)
  sort(unique(na.omit(as.character(x$GENE))))
}

set_jaccard <- function(a, b) {
  if (is.null(a) || is.null(b)) return(NA_real_)
  union_n <- length(union(a, b))
  if (union_n == 0L) return(1)
  length(intersect(a, b)) / union_n
}

exact_mcnemar <- function(long_sets, flye_sets, layer) {
  features <- sort(unique(c(unlist(long_sets, use.names = FALSE), unlist(flye_sets, use.names = FALSE))))
  if (!length(features)) return(tibble())
  out <- map_dfr(features, function(feature) {
    long_present <- map_lgl(long_sets, ~ feature %in% .x)
    flye_present <- map_lgl(flye_sets, ~ feature %in% .x)
    long_only <- sum(long_present & !flye_present)
    flye_only <- sum(!long_present & flye_present)
    discordant <- long_only + flye_only
    p <- if (discordant == 0L) 1 else
      stats::binom.test(long_only, discordant, p = 0.5, alternative = "two.sided")$p.value
    tibble(
      layer = layer,
      feature = feature,
      n_pairs = length(long_sets),
      longcycler_present = sum(long_present),
      flye_present = sum(flye_present),
      longcycler_only = long_only,
      flye_only = flye_only,
      discordant_pairs = discordant,
      exact_mcnemar_p = p
    )
  })
  out %>% mutate(BH_FDR = p.adjust(exact_mcnemar_p, method = "BH"))
}

run_rq11 <- function() {
  rq <- "RQ11"
  out_dir <- rq_dirs[[rq]]
  curated_map_path <- file.path("results", "vf", "gene_module_map.csv")
  if (!file.exists(curated_map_path)) {
    write_status(rq, "blocked", "curated_vf_definition_missing", curated_map_path)
    return(invisible(FALSE))
  }
  curated_map_hash <- sha256_file(curated_map_path)
  curated_map <- read_csv(curated_map_path, show_col_types = FALSE) %>%
    mutate(
      primary_assignment = tolower(as.character(primary_assignment)) %in% c("true", "t", "1", "yes"),
      mapping_sha256 = curated_map_hash,
      mapping_definition_version = "primary_assignment__high_or_moderate_confidence__non_unassigned_v1"
    )
  curated_required <- c("Gene", "module_id", "assignment_confidence", "primary_assignment")
  if (length(setdiff(curated_required, names(curated_map)))) {
    write_status(rq, "blocked", "curated_vf_definition_schema_invalid",
                 paste(setdiff(curated_required, names(curated_map)), collapse = ","))
    return(invisible(FALSE))
  }
  curated_gene_definition <- curated_map %>%
    filter(
      primary_assignment,
      assignment_confidence %in% c("High", "Moderate"),
      module_id != "unassigned"
    ) %>%
    distinct(Gene, module_id, assignment_confidence, mapping_sha256, mapping_definition_version)
  curated_genes <- curated_gene_definition$Gene
  if (!length(curated_genes)) {
    write_status(rq, "blocked", "curated_vf_definition_empty")
    return(invisible(FALSE))
  }
  atomic_write_csv(curated_gene_definition, file.path(out_dir, "curated_vf_definition.csv"))

  candidates <- read_csv(required_inputs[["candidates"]], show_col_types = FALSE) %>%
    mutate(
      Isolate_ID = as.character(Isolate_ID),
      Assembly_ID = as.character(Assembly_ID),
      Assembler = tolower(Assembler),
      full_path = norm_path(full_path),
      extension = tolower(tools::file_ext(full_path))
    )

  flye_candidates <- candidates %>%
    filter(
      Isolate_ID %in% manifest$Isolate_ID,
      Assembler == "flye", QC_PASS %in% TRUE,
      file.exists(full_path)
    ) %>%
    group_by(Isolate_ID) %>%
    arrange(
      if_else(extension == "fasta", 0L, 1L),
      desc(N50), n_contigs, Assembly_ID, .by_group = TRUE
    ) %>%
    mutate(
      n_qcpass_flye_candidates = n(),
      rq11_flye_rank = row_number(),
      rq11_selection_reason = if_else(
        rq11_flye_rank == 1L,
        "deterministic QC-pass Flye selection: prefer .fasta, higher N50, fewer contigs, Assembly_ID",
        "not selected for paired assembler comparison"
      )
    ) %>%
    ungroup()
  # Raw candidate paths/keys are reproducibility inputs, not public RQ output.
  atomic_write_csv(flye_candidates, file.path(input_root, "rq11_private_flye_candidate_selection_audit.csv"))
  atomic_write_csv(
    flye_candidates %>%
      count(extension, n_qcpass_flye_candidates, selected = rq11_flye_rank == 1L, name = "n_candidate_rows"),
    file.path(out_dir, "flye_candidate_selection_summary.csv")
  )
  old_public_audit <- file.path(out_dir, "flye_candidate_selection_audit.csv")
  if (file.exists(old_public_audit)) unlink(old_public_audit)
  flye_selected <- filter(flye_candidates, rq11_flye_rank == 1L)

  long_selected <- manifest %>%
    transmute(
      Isolate_ID, Participant_id, tp_lab,
      Assembly_ID, Assembler = "longcycler", full_path,
      N50, n_contigs, QC_PASS
    )
  pair_keys <- long_selected %>%
    inner_join(
      flye_selected %>%
        transmute(
          Isolate_ID,
          flye_Assembly_ID = Assembly_ID,
          flye_full_path = full_path,
          flye_N50 = N50,
          flye_n_contigs = n_contigs,
          flye_QC_PASS = QC_PASS,
          n_qcpass_flye_candidates
        ),
      by = "Isolate_ID"
    )
  if (nrow(pair_keys) != 527L || n_distinct(pair_keys$Isolate_ID) != 527L) {
    write_status(rq, "blocked", "unexpected_paired_qcpass_denominator",
                 sprintf("Observed %d rows and %d isolate keys; expected 527.",
                         nrow(pair_keys), n_distinct(pair_keys$Isolate_ID)))
    return(invisible(FALSE))
  }
  rq11_pair_key <- pair_keys %>%
    distinct(Isolate_ID, Participant_id, tp_lab) %>%
    arrange(Isolate_ID) %>%
    mutate(Paired_Isolate_Label = sprintf("Paired_Isolate_%03d", row_number())) %>%
    left_join(resident_key, by = "Participant_id")
  atomic_write_csv(rq11_pair_key, file.path(input_root, "rq11_private_paired_isolate_label_key.csv"))

  pair_long <- bind_rows(
    pair_keys %>% transmute(
      Isolate_ID, Participant_id, tp_lab, Assembler = "longcycler",
      Assembly_ID, full_path
    ),
    pair_keys %>% transmute(
      Isolate_ID, Participant_id, tp_lab, Assembler = "flye",
      Assembly_ID = flye_Assembly_ID, full_path = flye_full_path
    )
  )
  pair_hashes <- hash_paths(pair_long$full_path)
  pair_long <- pair_long %>% left_join(pair_hashes, by = "full_path")
  atomic_write_csv(pair_long, file.path(input_root, "rq11_paired_527_assemblies_sha256.csv"))

  paired_hashed <- pair_long %>%
    select(Isolate_ID, Participant_id, tp_lab, Assembler, Assembly_ID, full_path, fasta_sha256) %>%
    pivot_wider(
      names_from = Assembler,
      values_from = c(Assembly_ID, full_path, fasta_sha256),
      names_glue = "{.value}_{Assembler}"
    ) %>%
    left_join(
      rq11_pair_key %>% select(Isolate_ID, Paired_Isolate_Label, Resident_Label),
      by = "Isolate_ID"
    )

  mlst_exe <- Sys.which("mlst")
  abricate_exe <- Sys.which("abricate")
  db_table <- if (abricate_exe == "") tibble() else abricate_database_table(abricate_exe)
  db_vf <- if (nrow(db_table)) filter(db_table, DATABASE == "vfdb") else tibble()
  db_plasmid <- if (nrow(db_table)) filter(db_table, DATABASE == "plasmidfinder") else tibble()

  tools <- c(
    mlst_ecoli = unname(mlst_exe),
    vfdb = unname(abricate_exe),
    plasmidfinder = unname(abricate_exe)
  )
  versions <- c(
    mlst_ecoli = tool_version("mlst"),
    vfdb = tool_version("abricate"),
    plasmidfinder = tool_version("abricate")
  )
  db_signature <- function(row, db) {
    if (!nrow(row)) return(paste0(db, ":MISSING"))
    paste(row$DATABASE[[1]], row$SEQUENCES[[1]], row$DBTYPE[[1]], row$DATE[[1]], sep = ":")
  }
  signatures <- c(
    mlst_ecoli = paste("mlst", versions[["mlst_ecoli"]], "scheme=ecoli", sep = "|"),
    vfdb = paste("abricate", versions[["vfdb"]], db_signature(db_vf, "vfdb"), "mincov=80", "minid=80", sep = "|"),
    plasmidfinder = paste("abricate", versions[["plasmidfinder"]],
                          db_signature(db_plasmid, "plasmidfinder"), "mincov=80", "minid=80", sep = "|")
  )
  tool_checks <- tibble(
    layer = names(tools),
    executable = unname(tools),
    available = unname(tools) != "" & c(TRUE, nrow(db_vf) == 1L, nrow(db_plasmid) == 1L),
    version = unname(versions),
    signature = unname(signatures)
  )
  atomic_write_csv(tool_checks, file.path(out_dir, "tool_and_database_checks.csv"))

  # A missing tool/database creates explicit error rows and a blocked layer; it
  # never triggers use of the old 556-genome caches.
  unavailable <- tool_checks$layer[!tool_checks$available]
  if (length(unavailable)) tools[unavailable] <- ""
  unique_inputs <- pair_long %>%
    distinct(fasta_sha256, .keep_all = TRUE) %>%
    select(fasta_sha256, full_path)
  jobs <- tidyr::crossing(layer = names(tools), unique_inputs)
  cache_root <- file.path(input_root, "rq11_sha_bound_call_cache")
  dir.create(cache_root, recursive = TRUE, showWarnings = FALSE)
  message(sprintf("RQ11: evaluating %d unique SHA-bound tool calls with %d workers.", nrow(jobs), workers))
  job_list <- split(jobs, seq_len(nrow(jobs)))
  call_results <- parallel_lapply(job_list, function(job) {
    tryCatch(
      run_sha_bound_call(job, tools = tools, signatures = signatures, cache_root = cache_root),
      error = function(e) tibble(
        layer = as.character(job$layer[[1]]),
        input_sha256 = as.character(job$fasta_sha256[[1]]),
        representative_input_path = as.character(job$full_path[[1]]),
        tool_signature = signatures[[as.character(job$layer[[1]])]],
        output_path = NA_character_, sidecar_path = NA_character_,
        call_status = "runner_error", call_rows = NA_integer_, exit_code = 1L,
        cache_reused = FALSE, runner_error = conditionMessage(e)
      )
    )
  })
  unique_call_status <- bind_rows(call_results)
  atomic_write_csv(unique_call_status, file.path(input_root, "rq11_unique_call_status.csv"))

  assembly_call_status <- pair_long %>%
    tidyr::crossing(layer = names(tools)) %>%
    left_join(
      unique_call_status,
      by = c("layer", "fasta_sha256" = "input_sha256")
    ) %>%
    mutate(
      provenance_valid = case_when(
        layer == "mlst_ecoli" ~ exit_code == 0L & call_status %in% c("success", "no_call"),
        TRUE ~ exit_code == 0L & call_status %in% c("success", "zero_hit")
      )
    )
  atomic_write_csv(assembly_call_status, file.path(input_root, "rq11_private_assembly_call_status.csv"))
  old_public_call_status <- file.path(out_dir, "assembly_call_status.csv")
  if (file.exists(old_public_call_status)) unlink(old_public_call_status)

  layer_status <- assembly_call_status %>%
    group_by(layer) %>%
    summarise(
      expected_assembly_calls = n(),
      provenance_valid_assembly_calls = sum(provenance_valid),
      success_calls = sum(call_status == "success"),
      zero_hit_calls = sum(call_status == "zero_hit"),
      no_call_rows = sum(call_status == "no_call"),
      error_or_missing_calls = sum(!provenance_valid),
      cache_reused_calls = sum(cache_reused),
      .groups = "drop"
    ) %>%
    mutate(
      layer_analysis_status = if_else(
        provenance_valid_assembly_calls == expected_assembly_calls,
        "provenance_complete", "blocked_incomplete_provenance"
      ),
      blocked_reason = if_else(
        layer_analysis_status == "provenance_complete", NA_character_,
        "One or more exact FASTA/tool/database calls lacked a validated SHA-bound output; no stale cache substitution is allowed."
      )
    )

  status_lookup <- split(unique_call_status, unique_call_status$layer)
  status_vector <- function(layer, field) {
    x <- status_lookup[[layer]]
    setNames(x[[field]], x$input_sha256)
  }
  mlst_status <- status_vector("mlst_ecoli", "call_status")
  mlst_outputs <- status_vector("mlst_ecoli", "output_path")
  mlst_st <- setNames(
    map2_chr(mlst_outputs, mlst_status, parse_mlst_st),
    names(mlst_outputs)
  )
  feature_maps <- map(c("vfdb", "plasmidfinder"), function(layer) {
    statuses <- status_vector(layer, "call_status")
    outputs <- status_vector(layer, "output_path")
    sets <- map2(outputs, statuses, parse_abricate_features)
    names(sets) <- names(outputs)
    sets
  })
  names(feature_maps) <- c("vfdb", "plasmidfinder")

  mlst_pairs <- paired_hashed %>%
    transmute(
      Paired_Isolate_Label, Resident_Label,
      mlst_status_longcycler = unname(mlst_status[fasta_sha256_longcycler]),
      mlst_status_flye = unname(mlst_status[fasta_sha256_flye]),
      ST_longcycler = unname(mlst_st[fasta_sha256_longcycler]),
      ST_flye = unname(mlst_st[fasta_sha256_flye])
    ) %>%
    mutate(
      both_st_callable =
        mlst_status_longcycler == "success" & mlst_status_flye == "success" &
        !is.na(ST_longcycler) & !is.na(ST_flye) &
        !ST_longcycler %in% c("", "-") & !ST_flye %in% c("", "-"),
      exact_ST_agreement = if_else(both_st_callable, ST_longcycler == ST_flye, NA)
    )
  atomic_write_csv(mlst_pairs, file.path(out_dir, "paired_same_scheme_mlst.csv"))

  set_for <- function(map, sha) unname(map[sha])
  metric_layer <- function(layer) {
    sets <- feature_maps[[layer]]
    map_dfr(seq_len(nrow(paired_hashed)), function(i) {
      long_set <- sets[[paired_hashed$fasta_sha256_longcycler[[i]]]]
      flye_set <- sets[[paired_hashed$fasta_sha256_flye[[i]]]]
      tibble(
        Paired_Isolate_Label = paired_hashed$Paired_Isolate_Label[[i]],
        Resident_Label = paired_hashed$Resident_Label[[i]],
        layer = layer,
        call_pair_valid = !is.null(long_set) & !is.null(flye_set),
        n_longcycler = if (is.null(long_set)) NA_integer_ else length(long_set),
        n_flye = if (is.null(flye_set)) NA_integer_ else length(flye_set),
        count_difference_flye_minus_longcycler =
          if (is.null(long_set) || is.null(flye_set)) NA_integer_ else length(flye_set) - length(long_set),
        jaccard = set_jaccard(long_set, flye_set)
      )
    })
  }
  feature_metrics <- bind_rows(metric_layer("vfdb"), metric_layer("plasmidfinder"))
  atomic_write_csv(feature_metrics, file.path(out_dir, "paired_vf_replicon_metrics.csv"))

  paired_valid_counts <- feature_metrics %>%
    group_by(layer) %>%
    summarise(provenance_valid_pairs = sum(call_pair_valid), .groups = "drop")
  layer_status <- layer_status %>%
    left_join(paired_valid_counts, by = "layer") %>%
    mutate(
      provenance_valid_pairs = if_else(layer == "mlst_ecoli", 527L, provenance_valid_pairs),
      analysis_denominator = case_when(
        layer == "mlst_ecoli" ~ sum(mlst_pairs$both_st_callable),
        TRUE ~ provenance_valid_pairs
      )
    )
  atomic_write_csv(layer_status, file.path(out_dir, "layer_status_and_denominators.csv"))

  complete_feature_layers <- layer_status %>%
    filter(layer %in% c("vfdb", "plasmidfinder"), layer_analysis_status == "provenance_complete") %>%
    pull(layer)
  mcnemar_results <- map_dfr(complete_feature_layers, function(layer) {
    sets <- feature_maps[[layer]]
    long_sets <- unname(sets[paired_hashed$fasta_sha256_longcycler])
    flye_sets <- unname(sets[paired_hashed$fasta_sha256_flye])
    exact_mcnemar(long_sets, flye_sets, layer)
  })
  atomic_write_csv(mcnemar_results, file.path(out_dir, "feature_mcnemar_exact_BH.csv"))

  bland_altman_pairs <- feature_metrics %>%
    filter(call_pair_valid) %>%
    transmute(
      Paired_Isolate_Label, Resident_Label, layer,
      mean_feature_count = (n_longcycler + n_flye) / 2,
      difference_flye_minus_longcycler = count_difference_flye_minus_longcycler
    )
  bland_altman_summary <- bland_altman_pairs %>%
    group_by(layer) %>%
    summarise(
      n_pairs = n(),
      mean_bias_flye_minus_longcycler = mean(difference_flye_minus_longcycler),
      sd_difference = sd(difference_flye_minus_longcycler),
      lower_limit_agreement = mean_bias_flye_minus_longcycler - 1.96 * sd_difference,
      upper_limit_agreement = mean_bias_flye_minus_longcycler + 1.96 * sd_difference,
      median_difference = median(difference_flye_minus_longcycler),
      .groups = "drop"
    )
  atomic_write_csv(bland_altman_pairs, file.path(out_dir, "bland_altman_pairs.csv"))
  atomic_write_csv(bland_altman_summary, file.path(out_dir, "bland_altman_summary.csv"))

  agreement_summary <- bind_rows(
    tibble(
      endpoint = "same_scheme_ecoli_ST_exact_agreement",
      n_pairs = sum(mlst_pairs$both_st_callable),
      estimate = mean(mlst_pairs$exact_ST_agreement, na.rm = TRUE)
    ),
    feature_metrics %>%
      filter(call_pair_valid) %>%
      group_by(layer) %>%
      summarise(
        endpoint = paste0(layer, "_mean_jaccard"),
        n_pairs = n(), estimate = mean(jaccard), .groups = "drop"
      ) %>% select(endpoint, n_pairs, estimate)
  )

  set.seed(seed + 11L)
  bootstrap_indices <- replicate(n_boot, sample.int(nrow(paired_hashed), nrow(paired_hashed), replace = TRUE), simplify = FALSE)
  vf_wide <- feature_metrics %>% filter(layer == "vfdb")
  plasmid_wide <- feature_metrics %>% filter(layer == "plasmidfinder")
  bootstrap <- map_dfr(seq_along(bootstrap_indices), function(iteration) {
    idx <- bootstrap_indices[[iteration]]
    tibble(
      iteration = iteration,
      ST_exact_agreement = mean(mlst_pairs$exact_ST_agreement[idx], na.rm = TRUE),
      VF_mean_jaccard = mean(vf_wide$jaccard[idx], na.rm = TRUE),
      VF_mean_count_difference = mean(vf_wide$count_difference_flye_minus_longcycler[idx], na.rm = TRUE),
      replicon_mean_jaccard = mean(plasmid_wide$jaccard[idx], na.rm = TRUE),
      replicon_mean_count_difference = mean(plasmid_wide$count_difference_flye_minus_longcycler[idx], na.rm = TRUE)
    )
  })
  atomic_write_csv(bootstrap, file.path(out_dir, "isolate_bootstrap_distribution.csv"))
  bootstrap_summary <- bootstrap %>%
    pivot_longer(-iteration, names_to = "endpoint", values_to = "estimate") %>%
    group_by(endpoint) %>%
    summarise(
      bootstrap_median = median(estimate, na.rm = TRUE),
      bootstrap_ci_low = quantile(estimate, 0.025, na.rm = TRUE),
      bootstrap_ci_high = quantile(estimate, 0.975, na.rm = TRUE),
      valid_replicates = sum(is.finite(estimate)),
      .groups = "drop"
    )
  agreement_summary <- agreement_summary %>%
    left_join(
      bootstrap_summary %>%
        mutate(endpoint = recode(
          endpoint,
          ST_exact_agreement = "same_scheme_ecoli_ST_exact_agreement",
          VF_mean_jaccard = "vfdb_mean_jaccard",
          replicon_mean_jaccard = "plasmidfinder_mean_jaccard"
        )),
      by = "endpoint"
    )
  atomic_write_csv(bootstrap_summary, file.path(out_dir, "isolate_bootstrap_summary.csv"))
  atomic_write_csv(agreement_summary, file.path(out_dir, "agreement_summary.csv"))

  atomic_write_csv(tibble(
    input = c("analysis_manifest", "canonical_candidate_table", "paired_fasta_manifest"),
    path = c(required_inputs[["manifest"]], required_inputs[["candidates"]],
             file.path(input_root, "rq11_paired_527_assemblies_sha256.csv")),
    sha256 = c(sha256_file(required_inputs[["manifest"]]),
               sha256_file(required_inputs[["candidates"]]),
               sha256_file(file.path(input_root, "rq11_paired_527_assemblies_sha256.csv"))),
    role = c("current_532_Longcycler_keys", "QC-pass_Flye_counterpart_discovery", "527_exact_paired_FASTA_SHA256_keys")
  ), file.path(out_dir, "input_provenance.csv"))

  all_layers_complete <- all(layer_status$layer_analysis_status == "provenance_complete")
  atomic_write_lines(c(
    "# RQ11 — paired Longcycler/Flye measurement agreement",
    "",
    paste0("Status: ", if (all_layers_complete) "complete." else "partially blocked; see layer_status_and_denominators.csv."),
    "",
    "The 527 pairs are the current Longcycler analysis keys with a deterministic QC-pass Flye mate.",
    "Every comparison is based on a fresh or validated SHA-256-bound isolated cache. Calls use",
    "the same local E. coli MLST scheme, VFDB at 80% identity/80% coverage, and PlasmidFinder",
    "at 80% identity/80% coverage for both assemblers. Zero-hit outputs are valid and explicit;",
    "errors or missing tools block that layer instead of invoking the old 556-genome tables.",
    "Outputs include exact ST agreement, VF/replicon Jaccard and count differences, isolate-level",
    "bootstrap intervals, exact McNemar tests with BH correction, and Bland-Altman limits.",
    "AMR is deliberately excluded because RQ11 specifies no AMR comparison."
  ), file.path(out_dir, "README.md"))
  write_status(
    rq,
    if (all_layers_complete) "complete" else "partially_blocked",
    if (all_layers_complete) "all_SHA_bound_layers_complete" else "one_or_more_layers_incomplete",
    sprintf("527 paired isolates; %d unique SHA-bound calls evaluated", nrow(jobs))
  )
  invisible(all_layers_complete)
}

run_safely <- function(rq, fun) {
  tryCatch(
    fun(),
    error = function(e) {
      write_status(rq, "blocked", "unexpected_runtime_error", conditionMessage(e))
      atomic_write_lines(c(
        paste0("# ", rq, " blocked"), "",
        "The runner stopped this research question rather than suppressing an error or using stale output.",
        "", paste0("Error: ", conditionMessage(e))
      ), file.path(rq_dirs[[rq]], "README.md"))
      message(rq, " blocked: ", conditionMessage(e))
      FALSE
    }
  )
}

results <- c(
  RQ09 = run_safely("RQ09", run_rq09),
  RQ10 = run_safely("RQ10", run_rq10),
  RQ11 = run_safely("RQ11", run_rq11)
)

atomic_write_csv(tibble(
  research_question = names(results),
  completed_without_block = as.logical(results),
  run_timestamp_utc = run_timestamp,
  current_manifest_rows = nrow(manifest),
  provider_qc95_rows = nrow(provider_st)
), file.path(input_root, "RQ09_11_run_summary.csv"))

atomic_write_lines(capture.output(sessionInfo()), file.path(input_root, "RQ09_11_sessionInfo.txt"))
message("RQ09-RQ11 runner finished: ", paste(names(results), results, sep = "=", collapse = ", "))

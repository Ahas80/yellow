#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
})

project_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
if (!file.exists(file.path(project_root, "00_config.R"))) {
  stop("Run this script from the rUTIs project root.", call. = FALSE)
}

out_root <- file.path(project_root, "results", "research_questions")
dir.create(out_root, recursive = TRUE, showWarnings = FALSE)

atomic_write_csv <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  tmp <- tempfile(pattern = paste0(".", basename(path), "."), tmpdir = dirname(path), fileext = ".csv")
  readr::write_csv(x, tmp, na = "")
  if (!file.rename(tmp, path)) stop("Could not atomically publish ", path, call. = FALSE)
  invisible(path)
}

atomic_write_lines <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  tmp <- tempfile(pattern = paste0(".", basename(path), "."), tmpdir = dirname(path), fileext = ".txt")
  writeLines(x, tmp, useBytes = TRUE)
  if (!file.rename(tmp, path)) stop("Could not atomically publish ", path, call. = FALSE)
  invisible(path)
}

release_marker <- file.path(out_root, "RUN_COMPLETE.txt")
if (file.exists(release_marker) && !file.remove(release_marker)) {
  stop("Could not invalidate the previous release marker.", call. = FALSE)
}

# A release rerun must publish an exact inventory, not overwrite a subset of a
# previous run. Remove every active question directory and every generated
# shared-input artifact before running the modules. Only the two SHA-validated
# caches are retained; their sidecars are revalidated by their owning modules.
remove_generated_paths <- function(paths, label) {
  path_exists <- function(path) {
    link_target <- Sys.readlink(path)
    file.exists(path) || dir.exists(path) || (!is.na(link_target) && nzchar(link_target))
  }
  paths <- unique(paths[vapply(paths, path_exists, logical(1))])
  paths <- paths[order(nchar(paths), decreasing = TRUE)]
  if (!length(paths)) return(invisible(character()))
  status <- vapply(paths, unlink, integer(1), recursive = TRUE, force = TRUE)
  remaining <- paths[vapply(paths, path_exists, logical(1))]
  if (any(status != 0L) || length(remaining)) {
    stop(
      "Could not clear ", label, ": ",
      paste(unique(c(paths[status != 0L], remaining)), collapse = ", "),
      call. = FALSE
    )
  }
  invisible(paths)
}

active_question_dirs <- file.path(out_root, sprintf("RQ%02d", 1:10))
remove_generated_paths(active_question_dirs, "prior RQ01-RQ10 output directories")

input_root <- file.path(out_root, "_inputs")
preserved_cache_names <- c("vfdb_cache_sha256_v1", "rq09_mash")
prior_shared_inputs <- if (dir.exists(input_root)) {
  entries <- list.files(input_root, full.names = TRUE, all.files = TRUE, no.. = TRUE)
  entries[!basename(entries) %in% preserved_cache_names]
} else {
  character()
}
remove_generated_paths(prior_shared_inputs, "prior shared RQ input artifacts")
remove_generated_paths(
  file.path(out_root, c("final_contract_checks.csv", "final_question_status.csv")),
  "prior final RQ checks"
)

# Retired question artifacts must not coexist with the active RQ01-RQ10 tree,
# including caches or neutral-depth files carrying an old runner name.
retired_paths <- list.files(
  out_root,
  pattern = "(rq11|rq09_11)",
  recursive = TRUE,
  full.names = TRUE,
  all.files = TRUE,
  include.dirs = TRUE,
  ignore.case = TRUE,
  no.. = TRUE
)
remove_generated_paths(retired_paths, "retired RQ11/RQ09_11 artifacts")
remaining_retired <- list.files(
  out_root,
  pattern = "(rq11|rq09_11)",
  recursive = TRUE,
  full.names = TRUE,
  all.files = TRUE,
  include.dirs = TRUE,
  ignore.case = TRUE,
  no.. = TRUE
)
if (length(remaining_retired)) {
  stop("Retired RQ11/RQ09_11 artifacts remain after cleanup.", call. = FALSE)
}

# This is the release runner, not a development entry point.  Refuse a caller's
# reduced-resampling override and pass one common contract to every subprocess.
requested_boot <- Sys.getenv("RQ_BOOTSTRAP_REPS", unset = "")
requested_perm <- Sys.getenv("RQ_PERMUTATIONS", unset = "")
if (nzchar(requested_boot) && requested_boot != "10000") {
  stop("run_all.R requires RQ_BOOTSTRAP_REPS=10000 for publication outputs.", call. = FALSE)
}
if (nzchar(requested_perm) && requested_perm != "10000") {
  stop("run_all.R requires RQ_PERMUTATIONS=10000 for publication outputs.", call. = FALSE)
}
Sys.setenv(
  RQ_BOOTSTRAP_REPS = "10000",
  RQ_PERMUTATIONS = "10000",
  RQ_SEED = "20260712"
)

scripts <- file.path(
  project_root,
  "scripts",
  "research_questions",
  c("run_rq01_05.R", "run_rq06_08.R", "run_rq09_10.R")
)
missing_scripts <- scripts[!file.exists(scripts)]
if (length(missing_scripts)) {
  stop("Missing research-question script(s): ", paste(missing_scripts, collapse = ", "), call. = FALSE)
}

run_one <- function(path) {
  message("Running ", basename(path), " ...")
  status <- system2(file.path(R.home("bin"), "Rscript"), shQuote(path))
  if (!identical(status, 0L)) {
    stop("Research-question script failed: ", path, call. = FALSE)
  }
}

invisible(lapply(scripts, run_one))

# Final source-of-truth contract checks. These are deliberately repeated after
# every module so a partially stale research directory cannot be called final.
status <- read_csv(file.path(project_root, "results", "clinical", "status_map.csv"), show_col_types = FALSE) %>%
  filter(.data$analysis_include_primary %in% TRUE)
manifest <- read_csv(file.path(project_root, "results", "qc", "analysis_assembly_manifest.csv"), show_col_types = FALSE)
transitions <- read_csv(
  file.path(project_root, "results", "longitudinal", "longcycler_transitions.csv"),
  show_col_types = FALSE
)

episode_context <- status %>%
  transmute(
    Participant_id = as.character(.data$Participant_id),
    tp_lab = as.character(.data$tp_lab),
    Event_type = as.character(.data$Event_type),
    UTI_Status = as.character(.data$UTI_Status)
  )
manifest_context <- manifest %>%
  transmute(
    Participant_id = as.character(.data$Participant_id),
    tp_lab = as.character(.data$tp_lab)
  ) %>%
  left_join(episode_context, by = c("Participant_id", "tp_lab"), relationship = "one-to-one")
if (anyNA(manifest_context$UTI_Status) || anyNA(manifest_context$Event_type)) {
  stop("At least one selected assembly lacks authoritative clinical context.", call. = FALSE)
}

transition_context <- transitions %>%
  mutate(
    Participant_id = as.character(.data$Participant_id),
    tp_from = as.character(.data$tp_from),
    tp_to = as.character(.data$tp_to)
  ) %>%
  left_join(
    episode_context %>% select(Participant_id, tp_from = tp_lab, event_from = Event_type),
    by = c("Participant_id", "tp_from"), relationship = "many-to-one"
  ) %>%
  left_join(
    episode_context %>% select(Participant_id, tp_to = tp_lab, event_to = Event_type),
    by = c("Participant_id", "tp_to"), relationship = "many-to-one"
  ) %>%
  mutate(event_involved = .data$event_from == "UTI_event" | .data$event_to == "UTI_event")

all_pair_path <- file.path(out_root, "_inputs", "direct_pair_metrics_893.csv")
if (!file.exists(all_pair_path)) {
  stop("RQ06--RQ10 did not publish the required exact-denominator pair table.", call. = FALSE)
}
all_pairs <- read_csv(all_pair_path, show_col_types = FALSE)

checks <- tibble::tribble(
  ~check, ~observed, ~expected,
  "source_attrition_clinical_episodes", nrow(status), 583,
  "source_attrition_clinical_residents", dplyr::n_distinct(status$Participant_id), 166,
  "selected_longcycler_genomes", nrow(manifest), 532,
  "selected_longcycler_residents", dplyr::n_distinct(manifest$Participant_id), 161,
  "selected_operational_uti_genomes", sum(manifest_context$UTI_Status == "UTI"), 16,
  "selected_operational_not_uti_genomes", sum(manifest_context$UTI_Status == "Not_UTI"), 516,
  "adjacent_longcycler_pairs", nrow(transitions), 371,
  "adjacent_pair_residents", dplyr::n_distinct(transitions$Participant_id), 139,
  "adjacent_pairs_le25", sum(transitions$TotalSNPs <= 25, na.rm = TRUE), 140,
  "not_uti_to_uti_pairs", sum(transitions$status_from == "Not_UTI" & transitions$status_to == "UTI"), 9,
  "not_uti_to_uti_pairs_le25", sum(transitions$status_from == "Not_UTI" & transitions$status_to == "UTI" & transitions$TotalSNPs <= 25), 5,
  "routine_to_routine_intervals", sum(transition_context$event_from == "Routine" & transition_context$event_to == "Routine"), 322,
  "event_involved_intervals", sum(transition_context$event_involved), 49,
  "all_direct_within_resident_pairs", nrow(all_pairs), 893,
  "selected_uti_event_genomes", sum(manifest_context$Event_type == "UTI_event"), 32,
  "selected_uti_event_residents", dplyr::n_distinct(manifest_context$Participant_id[manifest_context$Event_type == "UTI_event"]), 29,
  "selected_uti_event_operational_uti", sum(manifest_context$Event_type == "UTI_event" & manifest_context$UTI_Status == "UTI"), 15,
  "selected_uti_event_operational_not_uti", sum(manifest_context$Event_type == "UTI_event" & manifest_context$UTI_Status == "Not_UTI"), 17
) %>%
  mutate(pass = .data$observed == .data$expected)

atomic_write_csv(checks, file.path(out_root, "final_contract_checks.csv"))
if (!all(checks$pass)) {
  stop("One or more final research-question denominator contracts failed.", call. = FALSE)
}

# A caught/structured downstream error is still a failed release.  Every
# question must explicitly report completion before this runner may emit PASS.
status_files <- file.path(out_root, sprintf("RQ%02d", 1:10), "analysis_status.csv")
if (any(!file.exists(status_files))) {
  stop("Missing per-question analysis status: ",
       paste(status_files[!file.exists(status_files)], collapse = ", "), call. = FALSE)
}
question_status <- bind_rows(lapply(status_files, function(path) {
  read_csv(path, show_col_types = FALSE) %>% mutate(status_file = path)
}))
if (nrow(question_status) != 10L || any(question_status$status != "complete")) {
  atomic_write_csv(question_status, file.path(out_root, "final_question_status.csv"))
  stop("At least one research question is not complete; no release PASS was written.", call. = FALSE)
}
atomic_write_csv(question_status, file.path(out_root, "final_question_status.csv"))

rq06_execution <- read_csv(file.path(out_root, "_inputs", "execution_status.csv"), show_col_types = FALSE)
if (nrow(rq06_execution) != 1L || rq06_execution$state[[1]] != "COMPLETE" ||
    rq06_execution$boot_reps[[1]] != 10000L) {
  stop("RQ06--RQ08 execution did not complete with 10,000 bootstrap replicates.", call. = FALSE)
}

provenance_01_05 <- lapply(sprintf("RQ%02d", 1:5), function(rq) {
  x <- read_csv(file.path(out_root, rq, "provenance.csv"), show_col_types = FALSE)
  as.integer(x$value[x$field == "bootstrap_reps"][[1]])
})
if (any(unlist(provenance_01_05) != 10000L)) {
  stop("RQ01--RQ05 provenance does not record 10,000 bootstrap replicates.", call. = FALSE)
}

# Input manifests may include human-readable exclusion notes, but no declared
# input path is allowed to point into the excluded Rowena tree.
input_manifests <- list.files(
  out_root,
  pattern = "(input_manifest|input_provenance|source_input_manifest).*\\.csv$",
  recursive = TRUE,
  full.names = TRUE,
  ignore.case = TRUE
)
path_violations <- character()
for (f in input_manifests) {
  x <- read_csv(f, show_col_types = FALSE)
  path_cols <- grep("path|file|input", names(x), value = TRUE, ignore.case = TRUE)
  if (!length(path_cols)) next
  vals <- unlist(x[path_cols], use.names = FALSE)
  vals <- vals[!is.na(vals)]
  if (any(grepl("Rowenas analysis", vals, fixed = TRUE))) path_violations <- c(path_violations, f)
}
if (length(path_violations)) {
  stop("Excluded Rowena input detected in: ", paste(unique(path_violations), collapse = ", "), call. = FALSE)
}


# Publication-facing de-identified tables may contain generated labels and the
# prespecified scientific measurements, but never raw project key columns.
deidentified_tables <- list.files(
  out_root,
  pattern = "(deidentified|case_matrix|case_table).*\\.csv$",
  recursive = TRUE,
  full.names = TRUE,
  ignore.case = TRUE
)
raw_key_names <- c("Participant_id", "tp_lab", "Episode_ID", "Isolate_ID", "Assembly_ID")
privacy_violations <- character()
for (f in deidentified_tables) {
  x <- read_csv(f, show_col_types = FALSE)
  if (any(names(x) %in% raw_key_names)) privacy_violations <- c(privacy_violations, f)
}
if (length(privacy_violations)) {
  stop("Raw project keys found in de-identified output(s): ",
       paste(unique(privacy_violations), collapse = ", "), call. = FALSE)
}

# The active RQ release may contain only selected Longcycler-linked provenance.
release_files <- list.files(out_root, recursive = TRUE, full.names = TRUE,
                            all.files = TRUE, no.. = TRUE)
forbidden_path <- release_files[grepl("flye", release_files, ignore.case = TRUE)]
if (length(forbidden_path)) {
  stop("Retired assembler path detected in active RQ release.", call. = FALSE)
}
text_files <- release_files[grepl("\\.(csv|tsv|txt|md|json|ndjson)$", release_files, ignore.case = TRUE)]
forbidden_content_pattern <- "(^|[^[:alnum:]])flye([^[:alnum:]]|$)"
forbidden_text <- text_files[vapply(text_files, function(path) {
  any(grepl(
    forbidden_content_pattern,
    readLines(path, warn = FALSE),
    ignore.case = TRUE
  ))
}, logical(1))]
if (length(forbidden_text)) {
  stop("Retired assembler content detected in active RQ release.", call. = FALSE)
}

atomic_write_lines(
  c(
    "Research-question analysis runner: PASS",
    sprintf("Generated: %s", format(Sys.time(), tz = "Europe/Amsterdam", usetz = TRUE)),
    "All RQ01-RQ10 modules completed successfully.",
    "All release runs used seed 20260712 and 10,000 resamples/permutations where prespecified.",
    "Final denominator contracts passed.",
    "No declared input path points into Rowenas analysis/."
  ),
  release_marker
)

message("Research-question analysis runner complete.")

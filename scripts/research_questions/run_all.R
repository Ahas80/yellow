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
  c("run_rq01_05.R", "run_rq06_08.R", "run_rq09_11.R")
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
  file.path(project_root, "results", "sensitivity", "longcycler_only", "longcycler_transitions.csv"),
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
    episode_context %>% select(.data$Participant_id, tp_from = .data$tp_lab, event_from = .data$Event_type),
    by = c("Participant_id", "tp_from"), relationship = "many-to-one"
  ) %>%
  left_join(
    episode_context %>% select(.data$Participant_id, tp_to = .data$tp_lab, event_to = .data$Event_type),
    by = c("Participant_id", "tp_to"), relationship = "many-to-one"
  ) %>%
  mutate(event_involved = .data$event_from == "UTI_event" | .data$event_to == "UTI_event")

all_pair_path <- file.path(out_root, "_inputs", "direct_pair_metrics_893.csv")
rq11_pair_path <- file.path(out_root, "_inputs", "rq11_paired_527_assemblies_sha256.csv")
if (!file.exists(all_pair_path) || !file.exists(rq11_pair_path)) {
  stop("RQ06--RQ11 did not publish their required exact-denominator input tables.", call. = FALSE)
}
all_pairs <- read_csv(all_pair_path, show_col_types = FALSE)
rq11_pairs <- read_csv(rq11_pair_path, show_col_types = FALSE)
rq11_episode_keys <- paste(rq11_pairs$Participant_id, rq11_pairs$tp_lab, sep = "||")

checks <- tibble::tribble(
  ~check, ~observed, ~expected,
  "eligible_clinical_episodes", nrow(status), 583,
  "eligible_clinical_residents", dplyr::n_distinct(status$Participant_id), 166,
  "operational_uti_episodes", sum(status$UTI_Status == "UTI"), 18,
  "operational_not_uti_episodes", sum(status$UTI_Status == "Not_UTI"), 565,
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
  "selected_uti_event_operational_not_uti", sum(manifest_context$Event_type == "UTI_event" & manifest_context$UTI_Status == "Not_UTI"), 17,
  "rq11_paired_assembly_rows", nrow(rq11_pairs), 1054,
  "rq11_paired_episode_keys", dplyr::n_distinct(rq11_episode_keys), 527,
  "rq11_longcycler_rows", sum(tolower(rq11_pairs$Assembler) == "longcycler"), 527,
  "rq11_flye_rows", sum(tolower(rq11_pairs$Assembler) == "flye"), 527
) %>%
  mutate(pass = .data$observed == .data$expected)

atomic_write_csv(checks, file.path(out_root, "final_contract_checks.csv"))
if (!all(checks$pass)) {
  stop("One or more final research-question denominator contracts failed.", call. = FALSE)
}

# A caught/structured downstream error is still a failed release.  Every
# question must explicitly report completion before this runner may emit PASS.
status_files <- file.path(out_root, sprintf("RQ%02d", 1:11), "analysis_status.csv")
if (any(!file.exists(status_files))) {
  stop("Missing per-question analysis status: ",
       paste(status_files[!file.exists(status_files)], collapse = ", "), call. = FALSE)
}
question_status <- bind_rows(lapply(status_files, function(path) {
  read_csv(path, show_col_types = FALSE) %>% mutate(status_file = path)
}))
if (nrow(question_status) != 11L || any(question_status$status != "complete")) {
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

atomic_write_lines(
  c(
    "Research-question analysis runner: PASS",
    sprintf("Generated: %s", format(Sys.time(), tz = "Europe/Amsterdam", usetz = TRUE)),
    "All RQ01--RQ11 modules completed successfully.",
    "All release runs used seed 20260712 and 10,000 resamples/permutations where prespecified.",
    "Final denominator contracts passed.",
    "No declared input path points into Rowenas analysis/."
  ),
  release_marker
)

message("Research-question analysis runner complete.")

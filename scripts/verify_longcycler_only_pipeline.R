#!/usr/bin/env Rscript
# Verify that active genomic outputs are derived exclusively from the selected,
# QC-passing Longcycler analysis manifest. Flye may remain in candidate/QC audit
# tables, but it must never appear in an active input or endpoint.

source("00_config.R")
suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
  library(tidyr)
})

args <- commandArgs(trailingOnly = TRUE)
stage_arg <- args[str_detect(args, "^--stage=")]
if (!length(stage_arg) && any(args == "--stage")) {
  i <- match("--stage", args)
  stage <- if (i < length(args)) args[[i + 1L]] else "final"
} else {
  stage <- if (length(stage_arg)) str_remove(stage_arg[[1]], "^--stage=") else "final"
}
stage <- tolower(stage)
if (!stage %in% c("upstream", "final")) stop("--stage must be upstream or final")

out_csv <- file.path(DIR_QC, "longcycler_only_pipeline_verification.csv")
out_txt <- file.path(DIR_QC, "longcycler_only_pipeline_verification.txt")
checks <- list()
add_check <- function(component, file, metric, observed, requirement, pass) {
  checks[[length(checks) + 1L]] <<- tibble(
    checked_at = as.character(Sys.time()),
    stage = stage,
    component = component,
    file = file,
    metric = metric,
    observed = as.character(observed),
    requirement = requirement,
    pass = isTRUE(pass)
  )
}

normalise_path <- function(x) normalizePath(as.character(x), winslash = "/", mustWork = FALSE)
episode_keys <- function(df) paste(as.character(df$Participant_id), normalise_timepoint_preserve_events(df$tp_lab), sep = "|")
require_file <- function(path, component) {
  if (!file.exists(path)) stop(component, " is missing: ", path)
  path
}

manifest_path <- require_file(FILE_ANALYSIS_ASSEMBLY_MANIFEST, "Longcycler analysis manifest")
manifest <- load_analysis_assemblies(manifest_path, require_files = TRUE) %>%
  mutate(full_path = normalise_path(full_path))
manifest_keys <- episode_keys(manifest)
expected_pairs <- manifest %>%
  count(Participant_id, name = "n") %>%
  summarise(n_pairs = sum(n * (n - 1) / 2)) %>%
  pull(n_pairs)
expected_transitions <- nrow(manifest) - n_distinct(manifest$Participant_id)
add_check("analysis_manifest", manifest_path, "rows", nrow(manifest), ">0", nrow(manifest) > 0)
add_check("analysis_manifest", manifest_path, "non_longcycler_rows", sum(normalise_assembler_column(manifest) != ANALYSIS_ASSEMBLER), "0", all(normalise_assembler_column(manifest) == ANALYSIS_ASSEMBLER))
add_check("analysis_manifest", manifest_path, "duplicate_episode_keys", sum(duplicated(manifest_keys)), "0", !anyDuplicated(manifest_keys))

canonical_path <- require_file(FILE_CANONICAL_ASSEMBLY_SELECTION, "Canonical candidate/QC table")
canonical <- read_csv(canonical_path, show_col_types = FALSE) %>%
  mutate(
    selected_canonical = as_pipeline_bool(selected_canonical),
    QC_PASS = as_pipeline_bool(QC_PASS)
  )
selected <- canonical %>% filter(selected_canonical, QC_PASS)
add_check("canonical_selection", canonical_path, "selected_rows", nrow(selected), as.character(nrow(manifest)), nrow(selected) == nrow(manifest))
add_check("canonical_selection", canonical_path, "selected_non_longcycler_rows", sum(normalise_assembler_column(selected) != ANALYSIS_ASSEMBLER), "0", all(normalise_assembler_column(selected) == ANALYSIS_ASSEMBLER))

vf_path <- require_file(FILE_VF_PA, "VF presence/absence matrix")
vf <- read_csv(vf_path, show_col_types = FALSE)
vf_keys <- episode_keys(vf)
add_check("vf", vf_path, "rows", nrow(vf), as.character(nrow(manifest)), nrow(vf) == nrow(manifest))
add_check("vf", vf_path, "episode_key_set_matches_manifest", setequal(vf_keys, manifest_keys), "TRUE", setequal(vf_keys, manifest_keys))

mlst_path <- require_file(FILE_MLST_CANONICAL, "Active MLST table")
mlst <- read_csv(mlst_path, show_col_types = FALSE)
if (!"full_path" %in% names(mlst)) stop(mlst_path, " lacks full_path")
mlst_paths <- normalise_path(mlst$full_path)
bad_provider <- if (all(c("ST_source", "provider_assembler") %in% names(mlst))) {
  mlst$ST_source == "provider_qc95" &
    (is.na(mlst$provider_assembler) | tolower(mlst$provider_assembler) != ANALYSIS_ASSEMBLER)
} else {
  rep(TRUE, nrow(mlst))
}
add_check("mlst", mlst_path, "rows", nrow(mlst), as.character(nrow(manifest)), nrow(mlst) == nrow(manifest))
add_check("mlst", mlst_path, "paths_outside_manifest", sum(!(mlst_paths %in% manifest$full_path)), "0", all(mlst_paths %in% manifest$full_path))
add_check("mlst", mlst_path, "non_longcycler_provider_rows", sum(bad_provider, na.rm = TRUE), "0", !any(bad_provider, na.rm = TRUE))

manifest_specs <- tribble(
  ~component, ~path,
  "core_snp", file.path(DIR_WGS, "core", "core_snp_input_manifest.csv"),
  "panaroo", file.path(DIR_WGS, "pan", "panaroo_input_manifest.csv")
)
for (i in seq_len(nrow(manifest_specs))) {
  path <- require_file(manifest_specs$path[[i]], paste0(manifest_specs$component[[i]], " manifest"))
  x <- read_csv(path, show_col_types = FALSE)
  x_assembler <- normalise_assembler_column(x)
  add_check(manifest_specs$component[[i]], path, "rows", nrow(x), as.character(nrow(manifest)), nrow(x) == nrow(manifest))
  add_check(manifest_specs$component[[i]], path, "non_longcycler_rows", sum(is.na(x_assembler) | x_assembler != ANALYSIS_ASSEMBLER), "0", all(!is.na(x_assembler) & x_assembler == ANALYSIS_ASSEMBLER))
}

if (stage == "final") {
  ready_path <- require_file(FILE_VF_READY, "VF analysis-ready table")
  ready <- read_csv(ready_path, show_col_types = FALSE)
  ready_keys <- episode_keys(ready)
  add_check("vf_ready", ready_path, "rows", nrow(ready), as.character(nrow(manifest)), nrow(ready) == nrow(manifest))
  add_check("vf_ready", ready_path, "episode_key_set_matches_manifest", setequal(ready_keys, manifest_keys), "TRUE", setequal(ready_keys, manifest_keys))

  pair_path <- require_file(file.path(DIR_STRAIN, "pairwise_metrics.csv"), "Pairwise strain metrics")
  pairwise <- read_csv(pair_path, show_col_types = FALSE)
  endpoint_longcycler <- tolower(pairwise$Assembler_A) == ANALYSIS_ASSEMBLER & tolower(pairwise$Assembler_B) == ANALYSIS_ASSEMBLER
  endpoint_paths <- c(normalise_path(pairwise$Fasta_path_A), normalise_path(pairwise$Fasta_path_B))
  add_check("pairwise", pair_path, "rows", nrow(pairwise), as.character(expected_pairs), nrow(pairwise) == expected_pairs)
  add_check("pairwise", pair_path, "non_longcycler_endpoint_rows", sum(!endpoint_longcycler), "0", all(endpoint_longcycler))
  add_check("pairwise", pair_path, "endpoint_paths_outside_manifest", sum(!(endpoint_paths %in% manifest$full_path)), "0", all(endpoint_paths %in% manifest$full_path))

  transition_path <- require_file(file.path(DIR_VF, "vf_longitudinal_transitions.csv"), "Longitudinal transition table")
  transitions <- read_csv(transition_path, show_col_types = FALSE)
  if ("cohort" %in% names(transitions)) transitions <- transitions %>% filter(cohort == "all")
  add_check("longitudinal", transition_path, "rows", nrow(transitions), as.character(expected_transitions), nrow(transitions) == expected_transitions)
}

result <- bind_rows(checks)
write_csv(result, out_csv)
failed <- result %>% filter(!pass)
writeLines(
  c(
    "Longcycler-only pipeline verification",
    paste0("Generated: ", Sys.time()),
    paste0("Stage: ", stage),
    paste0("Analysis manifest rows: ", nrow(manifest)),
    paste0("Participants: ", n_distinct(manifest$Participant_id)),
    paste0("Expected within-participant pairs: ", expected_pairs),
    paste0("Expected adjacent transitions: ", expected_transitions),
    paste0("Checks passed: ", sum(result$pass), "/", nrow(result)),
    if (nrow(failed)) paste0("FAILED: ", failed$component, " / ", failed$metric, " observed=", failed$observed, " required=", failed$requirement) else "PASS: no active Flye input or endpoint detected."
  ),
  out_txt
)
if (nrow(failed)) stop("Longcycler-only verification failed. See ", out_csv)
message("Longcycler-only verification passed: ", nrow(result), " checks. Report: ", out_csv)

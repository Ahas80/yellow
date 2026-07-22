#!/usr/bin/env Rscript

# Release gate for the selected, QC-passing Longcycler cohort. The verifier is
# intentionally independent of analysis summaries: it derives every expected
# denominator from the selected manifest and then checks the fixed release
# anchors requested for this run.

source("00_config.R")
suppressPackageStartupMessages({
  library(digest)
  library(dplyr)
  library(readr)
  library(stringr)
  library(tibble)
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
if (!stage %in% c("upstream", "final")) stop("--stage must be upstream or final", call. = FALSE)

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

out_csv <- file.path(DIR_QC, "longcycler_only_pipeline_verification.csv")
out_txt <- file.path(DIR_QC, "longcycler_only_pipeline_verification.txt")
checks <- list()

add_check <- function(component, file, metric, observed, requirement, pass) {
  checks[[length(checks) + 1L]] <<- tibble(
    checked_at = as.character(Sys.time()),
    stage = stage,
    component = component,
    file = as.character(file),
    metric = metric,
    observed = as.character(observed),
    requirement = as.character(requirement),
    pass = isTRUE(pass)
  )
}

normalise_path <- function(x, must_work = FALSE) {
  x <- as.character(x)
  out <- rep(NA_character_, length(x))
  ok <- !is.na(x) & nzchar(x)
  out[ok] <- normalizePath(x[ok], winslash = "/", mustWork = must_work)
  out
}

episode_keys <- function(df, pid = "Participant_id", tp = NULL) {
  if (is.null(tp)) {
    tp <- intersect(c("tp_lab", "Timepoint", "timepoint"), names(df))[[1]]
  }
  paste(as.character(df[[pid]]), normalise_timepoint_preserve_events(df[[tp]]), sep = "|")
}

require_file <- function(path, component) {
  if (!file.exists(path)) stop(component, " is missing: ", path, call. = FALSE)
  path
}

require_columns <- function(df, columns, component) {
  missing <- setdiff(columns, names(df))
  if (length(missing)) {
    stop(component, " lacks required columns: ", paste(missing, collapse = ", "), call. = FALSE)
  }
}

hashes_match <- function(paths, expected) {
  paths <- normalise_path(paths, must_work = TRUE)
  observed <- vapply(paths, digest, character(1), algo = "sha256", file = TRUE)
  !is.na(expected) & nzchar(expected) & tolower(observed) == tolower(as.character(expected))
}

manifest_path <- require_file(FILE_ANALYSIS_ASSEMBLY_MANIFEST, "Selected analysis manifest")
manifest <- load_analysis_assemblies(manifest_path, require_files = TRUE) %>%
  mutate(
    Participant_id = as.character(Participant_id),
    tp_lab = normalise_timepoint_preserve_events(tp_lab),
    full_path = normalise_path(full_path, must_work = TRUE)
  )
require_columns(manifest, c("Participant_id", "tp_lab", "full_path"), "Selected analysis manifest")
if (!"fasta_sha256" %in% names(manifest)) {
  stop("Selected analysis manifest lacks fasta_sha256 provenance.", call. = FALSE)
}

manifest_keys <- episode_keys(manifest)
manifest_paths <- manifest$full_path
manifest_path_hash <- paste(manifest$full_path, tolower(manifest$fasta_sha256), sep = "\n")
expected_pairs <- manifest %>%
  count(Participant_id, name = "n") %>%
  summarise(n_pairs = sum(n * (n - 1L) / 2L)) %>%
  pull(n_pairs)
expected_transitions <- nrow(manifest) - n_distinct(manifest$Participant_id)

add_check("analysis_manifest", manifest_path, "rows", nrow(manifest), EXPECTED_EPISODES, nrow(manifest) == EXPECTED_EPISODES)
add_check("analysis_manifest", manifest_path, "residents", n_distinct(manifest$Participant_id), EXPECTED_RESIDENTS, n_distinct(manifest$Participant_id) == EXPECTED_RESIDENTS)
add_check("analysis_manifest", manifest_path, "duplicate_episode_keys", sum(duplicated(manifest_keys)), 0L, !anyDuplicated(manifest_keys))
add_check("analysis_manifest", manifest_path, "duplicate_paths", sum(duplicated(manifest_paths)), 0L, !anyDuplicated(manifest_paths))
manifest_assembler <- normalise_assembler_column(manifest)
add_check("analysis_manifest", manifest_path, "rows_outside_selected_assembler", sum(is.na(manifest_assembler) | manifest_assembler != ANALYSIS_ASSEMBLER), 0L, all(!is.na(manifest_assembler) & manifest_assembler == ANALYSIS_ASSEMBLER))
add_check("analysis_manifest", manifest_path, "fasta_hashes_match", sum(hashes_match(manifest$full_path, manifest$fasta_sha256)), EXPECTED_EPISODES, all(hashes_match(manifest$full_path, manifest$fasta_sha256)))
add_check("analysis_manifest", manifest_path, "derived_all_pairs", expected_pairs, EXPECTED_ALL_PAIRS, expected_pairs == EXPECTED_ALL_PAIRS)
add_check("analysis_manifest", manifest_path, "derived_adjacent_pairs", expected_transitions, EXPECTED_ADJACENT_PAIRS, expected_transitions == EXPECTED_ADJACENT_PAIRS)

metadata_path <- require_file(FILE_ASSEMBLY_META_ALL, "Assembly metadata")
metadata <- read_csv(metadata_path, show_col_types = FALSE)
metadata_assembler <- normalise_assembler_column(metadata)
add_check("assembly_metadata", metadata_path, "rows_outside_selected_assembler", sum(is.na(metadata_assembler) | metadata_assembler != ANALYSIS_ASSEMBLER), 0L, all(!is.na(metadata_assembler) & metadata_assembler == ANALYSIS_ASSEMBLER))

canonical_path <- require_file(FILE_CANONICAL_ASSEMBLY_SELECTION, "Candidate/QC table")
canonical <- read_csv(canonical_path, show_col_types = FALSE) %>%
  mutate(
    selected_canonical = as_pipeline_bool(selected_canonical),
    QC_PASS = as_pipeline_bool(QC_PASS)
  )
canonical_assembler <- normalise_assembler_column(canonical)
selected <- canonical %>%
  filter(selected_canonical, QC_PASS) %>%
  mutate(full_path = normalise_path(full_path, must_work = TRUE))
selected_keys <- episode_keys(selected)
add_check("canonical_selection", canonical_path, "all_rows_outside_selected_assembler", sum(is.na(canonical_assembler) | canonical_assembler != ANALYSIS_ASSEMBLER), 0L, all(!is.na(canonical_assembler) & canonical_assembler == ANALYSIS_ASSEMBLER))
add_check("canonical_selection", canonical_path, "selected_rows", nrow(selected), EXPECTED_EPISODES, nrow(selected) == EXPECTED_EPISODES)
add_check("canonical_selection", canonical_path, "selected_keys_match_manifest", setequal(selected_keys, manifest_keys), TRUE, setequal(selected_keys, manifest_keys))
add_check("canonical_selection", canonical_path, "selected_paths_match_manifest", setequal(selected$full_path, manifest_paths), TRUE, setequal(selected$full_path, manifest_paths))

cohort_path <- require_file(FILE_ANALYSIS_CLINICAL_COHORT, "Selected clinical analysis cohort")
cohort <- read_csv(cohort_path, show_col_types = FALSE) %>%
  mutate(
    Participant_id = as.character(Participant_id),
    tp_lab = normalise_timepoint_preserve_events(tp_lab)
  )
cohort_mtime <- file.info(cohort_path)$mtime
status_column <- intersect(c("uti_status", "UTI_Status", "status", "Status"), names(cohort))
if (!length(status_column)) stop("Selected clinical analysis cohort lacks a UTI status column.", call. = FALSE)
status_column <- status_column[[1]]
cohort_keys <- episode_keys(cohort)
cohort_status <- as.character(cohort[[status_column]])
add_check("analysis_cohort", cohort_path, "rows", nrow(cohort), EXPECTED_EPISODES, nrow(cohort) == EXPECTED_EPISODES)
add_check("analysis_cohort", cohort_path, "residents", n_distinct(cohort$Participant_id), EXPECTED_RESIDENTS, n_distinct(cohort$Participant_id) == EXPECTED_RESIDENTS)
add_check("analysis_cohort", cohort_path, "episode_keys_match_manifest", setequal(cohort_keys, manifest_keys), TRUE, setequal(cohort_keys, manifest_keys))
add_check("analysis_cohort", cohort_path, "UTI_rows", sum(cohort_status == "UTI", na.rm = TRUE), EXPECTED_UTI, sum(cohort_status == "UTI", na.rm = TRUE) == EXPECTED_UTI)
add_check("analysis_cohort", cohort_path, "Not_UTI_rows", sum(cohort_status == "Not_UTI", na.rm = TRUE), EXPECTED_NOT_UTI, sum(cohort_status == "Not_UTI", na.rm = TRUE) == EXPECTED_NOT_UTI)
if (all(c("fasta_path", "fasta_sha256") %in% names(cohort))) {
  cohort_path_hash <- paste(normalise_path(cohort$fasta_path, must_work = TRUE), tolower(cohort$fasta_sha256), sep = "\n")
  add_check("analysis_cohort", cohort_path, "path_hash_set_matches_manifest", setequal(cohort_path_hash, manifest_path_hash), TRUE, setequal(cohort_path_hash, manifest_path_hash))
} else {
  add_check("analysis_cohort", cohort_path, "path_hash_provenance_columns", "missing", "fasta_path and fasta_sha256", FALSE)
}

check_episode_table <- function(component, path, require_path_hash = FALSE) {
  path <- require_file(path, component)
  x <- read_csv(path, show_col_types = FALSE)
  keys <- episode_keys(x)
  add_check(component, path, "rows", nrow(x), EXPECTED_EPISODES, nrow(x) == EXPECTED_EPISODES)
  add_check(component, path, "episode_key_set_matches_manifest", setequal(keys, manifest_keys), TRUE, setequal(keys, manifest_keys))
  if (require_path_hash) {
    path_col <- intersect(c("full_path", "fasta_path"), names(x))
    sha_col <- intersect(c("fasta_sha256", "Fasta_SHA256"), names(x))
    if (!length(path_col) || !length(sha_col)) {
      add_check(component, path, "path_hash_provenance_columns", "missing", "path and SHA-256", FALSE)
    } else {
      observed <- paste(normalise_path(x[[path_col[[1]]]], must_work = TRUE), tolower(x[[sha_col[[1]]]]), sep = "\n")
      add_check(component, path, "path_hash_set_matches_manifest", setequal(observed, manifest_path_hash), TRUE, setequal(observed, manifest_path_hash))
    }
  }
  invisible(x)
}

vf <- check_episode_table("vf", FILE_VF_PA)
mlst <- check_episode_table("mlst", FILE_MLST_CANONICAL, require_path_hash = TRUE)
if (!"ST_source" %in% names(mlst)) {
  add_check("mlst", FILE_MLST_CANONICAL, "source_provenance", "missing", "ST_source present", FALSE)
} else {
  allowed_sources <- c(
    "provider_qc95",
    "local_fallback_provider_conflict",
    "local_fallback_provider_missing",
    "missing_provider_conflict",
    "missing"
  )
  add_check("mlst", FILE_MLST_CANONICAL, "unrecognised_source_rows", sum(!mlst$ST_source %in% allowed_sources), 0L, all(mlst$ST_source %in% allowed_sources))
}

manifest_specs <- tribble(
  ~component, ~path, ~path_candidates, ~sha_candidates,
  "plasmidfinder", file.path(DIR_PLASMIDS, "plasmidfinder_input_manifest.csv"), list(c("fasta_path", "full_path")), list(c("fasta_sha256", "Fasta_SHA256")),
  "core_snp", file.path(DIR_WGS_CORE, "core_snp_input_manifest.csv"), list(c("fasta_path", "full_path")), list(c("fasta_sha256", "Fasta_SHA256")),
  "panaroo", file.path(DIR_WGS_PAN, "panaroo_input_manifest.csv"), list(c("fasta_path", "full_path")), list(c("fasta_sha256", "Fasta_SHA256"))
)
panaroo_expected_samples <- character()

for (i in seq_len(nrow(manifest_specs))) {
  component <- manifest_specs$component[[i]]
  path <- require_file(manifest_specs$path[[i]], paste0(component, " input manifest"))
  x <- read_csv(path, show_col_types = FALSE)
  # tribble() stores the candidate vectors as one nested list per row. Flatten
  # that cell before matching, otherwise intersect() compares a list object to
  # character column names and falsely reports valid provenance as missing.
  path_candidates <- unlist(manifest_specs$path_candidates[[i]], use.names = FALSE)
  sha_candidates <- unlist(manifest_specs$sha_candidates[[i]], use.names = FALSE)
  path_col <- intersect(path_candidates, names(x))
  sha_col <- intersect(sha_candidates, names(x))
  assembler <- normalise_assembler_column(x)
  add_check(component, path, "rows", nrow(x), EXPECTED_EPISODES, nrow(x) == EXPECTED_EPISODES)
  add_check(component, path, "rows_outside_selected_assembler", sum(is.na(assembler) | assembler != ANALYSIS_ASSEMBLER), 0L, all(!is.na(assembler) & assembler == ANALYSIS_ASSEMBLER))
  if (!length(path_col) || !length(sha_col)) {
    add_check(component, path, "path_hash_provenance_columns", "missing", "path and SHA-256", FALSE)
  } else {
    observed <- paste(normalise_path(x[[path_col[[1]]]], must_work = TRUE), tolower(x[[sha_col[[1]]]]), sep = "\n")
    add_check(component, path, "path_hash_set_matches_manifest", setequal(observed, manifest_path_hash), TRUE, setequal(observed, manifest_path_hash))
    add_check(component, path, "file_hashes_match", sum(hashes_match(x[[path_col[[1]]]], x[[sha_col[[1]]]])), EXPECTED_EPISODES, all(hashes_match(x[[path_col[[1]]]], x[[sha_col[[1]]]])))
  }
  if (component == "panaroo") {
    require_columns(x, c("Assembly_Base_ID", "gff_path", "gff_sha256", "gff_available"), "Panaroo input manifest")
    gff_ok <- as_pipeline_bool(x$gff_available) & file.exists(x$gff_path) & hashes_match(x$gff_path, x$gff_sha256)
    add_check(component, path, "valid_exact_gffs", sum(gff_ok), EXPECTED_EPISODES, all(gff_ok))
    panaroo_expected_samples <- sub(
      "\\.gff$",
      "",
      basename(as.character(x$gff_path)),
      ignore.case = TRUE
    )
  }
}

# The mechanism analysis consumes Panaroo's Roary-compatible table rather than
# the primary table.  Its sample columns therefore belong in the release
# contract too; otherwise an old mixed-cohort compatibility export could be
# combined with a current manifest.
panaroo_roary_path <- require_file(
  file.path(DIR_WGS_PAN, "gene_presence_absence_roary.csv"),
  "Panaroo Roary-compatible presence/absence table"
)
panaroo_metadata_columns <- c(
  "Gene", "Non-unique Gene name", "Annotation", "No. isolates",
  "No. sequences", "Avg sequences per isolate", "Genome Fragment",
  "Order within Fragment", "Accessory Fragment",
  "Accessory Order with Fragment", "QC", "Min group size nuc",
  "Max group size nuc", "Avg group size nuc"
)
panaroo_roary_columns <- names(read_csv(
  panaroo_roary_path,
  n_max = 0,
  show_col_types = FALSE
))
panaroo_roary_samples <- setdiff(panaroo_roary_columns, panaroo_metadata_columns)
add_check(
  "panaroo",
  panaroo_roary_path,
  "roary_sample_columns",
  length(panaroo_roary_samples),
  EXPECTED_EPISODES,
  setequal(panaroo_roary_samples, panaroo_expected_samples)
)
add_check(
  "panaroo",
  panaroo_roary_path,
  "roary_sample_columns_outside_manifest",
  sum(!panaroo_roary_samples %in% panaroo_expected_samples),
  0L,
  all(panaroo_roary_samples %in% panaroo_expected_samples)
)

# Canonical PlasmidFinder outputs are exact GENE-label profiles with an
# explicit successful/no-hit status for every selected assembly.
pf_marker <- require_file(
  file.path(DIR_PLASMIDS, "PLASMIDFINDER_RUN_COMPLETE.txt"),
  "PlasmidFinder completion marker"
)
pf_run_path <- require_file(
  file.path(DIR_PLASMIDS, "plasmidfinder_run_manifest.csv"),
  "PlasmidFinder run manifest"
)
pf_run <- read_csv(pf_run_path, show_col_types = FALSE)
require_columns(
  pf_run,
  c("Isolate_ID", "fasta_path", "fasta_sha256", "call_status", "no_hit"),
  "PlasmidFinder run manifest"
)
pf_path_hash <- paste(
  normalise_path(pf_run$fasta_path, must_work = TRUE),
  tolower(pf_run$fasta_sha256),
  sep = "\n"
)
add_check("plasmidfinder", pf_run_path, "completed_calls",
          sum(pf_run$call_status == "complete"), EXPECTED_EPISODES,
          nrow(pf_run) == EXPECTED_EPISODES &&
            all(pf_run$call_status == "complete"))
add_check("plasmidfinder", pf_run_path, "valid_no_hit_calls",
          sum(as_pipeline_bool(pf_run$no_hit)), 110L,
          sum(as_pipeline_bool(pf_run$no_hit)) == 110L)
add_check("plasmidfinder", pf_run_path, "path_hash_set_matches_manifest",
          setequal(pf_path_hash, manifest_path_hash), TRUE,
          setequal(pf_path_hash, manifest_path_hash))

pf_pa_path <- require_file(
  file.path(DIR_PLASMIDS, "plasmidfinder_presence_absence.csv"),
  "PlasmidFinder GENE-label matrix"
)
pf_pa <- read_csv(pf_pa_path, show_col_types = FALSE)
pf_catalog_path <- require_file(
  file.path(DIR_PLASMIDS, "plasmidfinder_replicon_catalog.csv"),
  "PlasmidFinder replicon catalog"
)
pf_catalog <- read_csv(pf_catalog_path, show_col_types = FALSE)
pf_features <- setdiff(names(pf_pa), "Isolate_ID")
add_check("plasmidfinder", pf_pa_path, "episode_rows",
          nrow(pf_pa), EXPECTED_EPISODES,
          nrow(pf_pa) == EXPECTED_EPISODES &&
            setequal(pf_pa$Isolate_ID, manifest$Isolate_ID))
add_check("plasmidfinder", pf_pa_path, "GENE_features",
          length(pf_features), 42L,
          length(pf_features) == 42L &&
            all(pf_features %in% pf_catalog$GENE) &&
            !any(pf_features %in% pf_catalog$accession))
ap001918 <- pf_catalog %>%
  filter(accession == "AP001918") %>%
  distinct(GENE) %>%
  pull(GENE)
add_check("plasmidfinder", pf_catalog_path, "AP001918_separate_GENE_labels",
          paste(sort(ap001918), collapse = ";"),
          "IncFIA_1;IncFIB(AP001918)_1;IncFIC(FII)_1",
          all(c(
            "IncFIA_1", "IncFIB(AP001918)_1", "IncFIC(FII)_1"
          ) %in% ap001918))

# MOB-suite is an assembly-only prediction layer. Every selected input contig
# must be retained as predicted plasmid, predicted chromosome, or unassigned.
mob_root <- file.path(DIR_PLASMIDS, "mob_suite")
mob_marker <- require_file(
  file.path(mob_root, "RUN_COMPLETE.txt"),
  "MOB-suite completion marker"
)
mob_status_path <- require_file(
  file.path(mob_root, "sample_status.csv"),
  "MOB-suite sample status"
)
mob_status <- read_csv(mob_status_path, show_col_types = FALSE)
mob_contigs_path <- require_file(
  file.path(mob_root, "contig_assignments.csv"),
  "MOB-suite contig assignments"
)
mob_contigs <- read_csv(mob_contigs_path, show_col_types = FALSE)
mob_profiles_path <- require_file(
  file.path(mob_root, "episode_plasmid_profiles.csv"),
  "MOB-suite episode profiles"
)
mob_profiles <- read_csv(mob_profiles_path, show_col_types = FALSE)
require_columns(
  mob_contigs,
  c("Isolate_ID", "contig_id", "molecule_type", "fasta_path", "fasta_sha256"),
  "MOB contig assignments"
)
allowed_molecules <- c(
  "predicted_plasmid", "predicted_chromosome", "unassigned"
)
expected_input_contigs <- sum(as.integer(manifest$n_contigs))
add_check("mob_suite", mob_status_path, "completed_samples",
          sum(mob_status$status == "complete"), EXPECTED_EPISODES,
          nrow(mob_status) == EXPECTED_EPISODES &&
            all(mob_status$status == "complete") &&
            setequal(mob_status$Isolate_ID, manifest$Isolate_ID))
add_check("mob_suite", mob_profiles_path, "episode_profiles",
          nrow(mob_profiles), EXPECTED_EPISODES,
          nrow(mob_profiles) == EXPECTED_EPISODES &&
            setequal(mob_profiles$Isolate_ID, manifest$Isolate_ID))
add_check("mob_suite", mob_contigs_path, "input_contigs_accounted",
          nrow(mob_contigs), expected_input_contigs,
          nrow(mob_contigs) == expected_input_contigs &&
            !anyDuplicated(mob_contigs[c("Isolate_ID", "contig_id")]) &&
            all(mob_contigs$molecule_type %in% allowed_molecules))
mob_marker_text <- paste(readLines(mob_marker, warn = FALSE), collapse = "\n")
add_check("mob_suite", mob_marker, "interpretation_scope",
          mob_marker_text,
          "assembly-based prediction; no transfer/transmission claim",
          str_detect(mob_marker_text, fixed("assembly-based predicted")) &&
            str_detect(mob_marker_text, fixed("no circularity, transfer, transmission")))

if (stage == "final") {
  ready <- check_episode_table("vf_ready", FILE_VF_READY)

  # The binary VF-hit table is consumed directly by downstream models.  CSV
  # endpoint checks alone cannot prove that a stale or mixed-assembler RDS was
  # not retained, so validate its endpoint/path contract independently.
  vf_hits_path <- require_file(FILE_VF_HITS, "VF hit RDS")
  vf_hits <- readRDS(vf_hits_path)
  if (!is.data.frame(vf_hits)) {
    stop("VF hit RDS must contain a data frame: ", vf_hits_path, call. = FALSE)
  }
  require_columns(vf_hits, c("Participant_id", "tp_lab", "full_path"), "VF hit RDS")
  vf_hit_keys <- episode_keys(vf_hits)
  vf_hit_paths <- normalise_path(vf_hits$full_path, must_work = TRUE)
  vf_hit_assembler <- normalise_assembler_column(vf_hits)
  add_check("vf_hits_rds", vf_hits_path, "rows", nrow(vf_hits), ">0", nrow(vf_hits) > 0L)
  add_check("vf_hits_rds", vf_hits_path, "endpoint_keys_outside_manifest", sum(!vf_hit_keys %in% manifest_keys), 0L, all(vf_hit_keys %in% manifest_keys))
  add_check("vf_hits_rds", vf_hits_path, "paths_outside_manifest", sum(!vf_hit_paths %in% manifest_paths), 0L, all(vf_hit_paths %in% manifest_paths))
  add_check("vf_hits_rds", vf_hits_path, "rows_outside_selected_assembler", sum(is.na(vf_hit_assembler) | vf_hit_assembler != ANALYSIS_ASSEMBLER), 0L, all(!is.na(vf_hit_assembler) & vf_hit_assembler == ANALYSIS_ASSEMBLER))
  vf_hits_mtime <- file.info(vf_hits_path)$mtime
  add_check("vf_hits_rds", vf_hits_path, "not_older_than_selected_cohort", as.character(vf_hits_mtime), paste0(">=", cohort_mtime), !is.na(vf_hits_mtime) && vf_hits_mtime >= cohort_mtime)

  count_audit_path <- require_file(
    file.path(DIR_RESULTS, "audit", "uti_status_count_explanation.md"),
    "Selected-cohort UTI/Not_UTI count audit"
  )
  count_audit_mtime <- file.info(count_audit_path)$mtime
  add_check(
    "count_audit",
    count_audit_path,
    "not_older_than_selected_cohort",
    as.character(count_audit_mtime),
    paste0(">=", cohort_mtime),
    !is.na(count_audit_mtime) && count_audit_mtime >= cohort_mtime
  )

  pair_path <- require_file(file.path(DIR_STRAIN, "pairwise_metrics.csv"), "Pairwise strain metrics")
  pairwise <- read_csv(pair_path, show_col_types = FALSE)
  require_columns(pairwise, c("Participant_id_A", "Participant_id_B", "Timepoint_A", "Timepoint_B", "Fasta_path_A", "Fasta_path_B", "Fasta_SHA256_A", "Fasta_SHA256_B", "TotalSNPs"), "Pairwise strain metrics")
  pair_path_hash <- c(
    paste(normalise_path(pairwise$Fasta_path_A, must_work = TRUE), tolower(pairwise$Fasta_SHA256_A), sep = "\n"),
    paste(normalise_path(pairwise$Fasta_path_B, must_work = TRUE), tolower(pairwise$Fasta_SHA256_B), sep = "\n")
  )
  pair_keys_a <- paste(pairwise$Participant_id_A, normalise_timepoint_preserve_events(pairwise$Timepoint_A), sep = "|")
  pair_keys_b <- paste(pairwise$Participant_id_B, normalise_timepoint_preserve_events(pairwise$Timepoint_B), sep = "|")
  add_check("pairwise", pair_path, "rows", nrow(pairwise), EXPECTED_ALL_PAIRS, nrow(pairwise) == EXPECTED_ALL_PAIRS)
  add_check("pairwise", pair_path, "cross_resident_rows", sum(pairwise$Participant_id_A != pairwise$Participant_id_B), 0L, all(pairwise$Participant_id_A == pairwise$Participant_id_B))
  add_check("pairwise", pair_path, "endpoint_keys_outside_manifest", sum(!c(pair_keys_a, pair_keys_b) %in% manifest_keys), 0L, all(c(pair_keys_a, pair_keys_b) %in% manifest_keys))
  add_check("pairwise", pair_path, "endpoint_path_hashes_outside_manifest", sum(!pair_path_hash %in% manifest_path_hash), 0L, all(pair_path_hash %in% manifest_path_hash))
  add_check("pairwise", pair_path, "missing_SNP_values", sum(is.na(pairwise$TotalSNPs)), 0L, !anyNA(pairwise$TotalSNPs))

  transition_path <- require_file(file.path(DIR_RESULTS, "longitudinal", "longcycler_transitions.csv"), "Canonical adjacent transitions")
  transitions <- read_csv(transition_path, show_col_types = FALSE)
  require_columns(transitions, c("Participant_id", "status_from", "status_to", "TotalSNPs"), "Canonical adjacent transitions")
  add_check("longitudinal", transition_path, "rows", nrow(transitions), EXPECTED_ADJACENT_PAIRS, nrow(transitions) == EXPECTED_ADJACENT_PAIRS)
  add_check("longitudinal", transition_path, "residents", n_distinct(transitions$Participant_id), EXPECTED_ADJACENT_RESIDENTS, n_distinct(transitions$Participant_id) == EXPECTED_ADJACENT_RESIDENTS)
  add_check("longitudinal", transition_path, "pairs_at_or_below_25", sum(transitions$TotalSNPs <= SAME_STRAIN_SNP_THRESHOLD, na.rm = TRUE), EXPECTED_ADJACENT_LE25, sum(transitions$TotalSNPs <= SAME_STRAIN_SNP_THRESHOLD, na.rm = TRUE) == EXPECTED_ADJACENT_LE25)
  n_to_uti <- sum(transitions$status_from == "Not_UTI" & transitions$status_to == "UTI")
  n_to_uti_le25 <- sum(transitions$status_from == "Not_UTI" & transitions$status_to == "UTI" & transitions$TotalSNPs <= SAME_STRAIN_SNP_THRESHOLD)
  add_check("longitudinal", transition_path, "Not_UTI_to_UTI_pairs", n_to_uti, EXPECTED_NOT_UTI_TO_UTI, n_to_uti == EXPECTED_NOT_UTI_TO_UTI)
  add_check("longitudinal", transition_path, "Not_UTI_to_UTI_pairs_at_or_below_25", n_to_uti_le25, EXPECTED_NOT_UTI_TO_UTI_LE25, n_to_uti_le25 == EXPECTED_NOT_UTI_TO_UTI_LE25)

  amr_profile_path <- require_file(
    file.path(DIR_RESULTS, "amr", "episode_amr_profiles.csv"),
    "Script-29 genomic AMR episode profiles"
  )
  amr_profiles <- read_csv(amr_profile_path, show_col_types = FALSE)
  require_columns(
    amr_profiles,
    c(
      "Participant_id", "tp_lab", "fasta_path", "fasta_sha256",
      "amr_gene_count_informative", "mdfA_detected",
      "any_informative_acquired_amr", "interpretation_scope"
    ),
    "Script-29 genomic AMR episode profiles"
  )
  amr_keys <- episode_keys(amr_profiles)
  amr_path_hash <- paste(
    normalise_path(amr_profiles$fasta_path, must_work = TRUE),
    tolower(amr_profiles$fasta_sha256), sep = "\n"
  )
  add_check("genomic_amr", amr_profile_path, "rows", nrow(amr_profiles),
            EXPECTED_EPISODES, nrow(amr_profiles) == EXPECTED_EPISODES)
  add_check(
    "genomic_amr", amr_profile_path, "episode_key_set_matches_manifest",
    setequal(amr_keys, manifest_keys), TRUE, setequal(amr_keys, manifest_keys)
  )
  add_check(
    "genomic_amr", amr_profile_path, "path_hash_set_matches_manifest",
    setequal(amr_path_hash, manifest_path_hash), TRUE,
    setequal(amr_path_hash, manifest_path_hash)
  )
  add_check(
    "genomic_amr", amr_profile_path, "genomic_not_AST_scope",
    unique(amr_profiles$interpretation_scope),
    "genomic determinants—not phenotypic AST",
    all(str_detect(
      amr_profiles$interpretation_scope,
      regex("not phenotypic AST", ignore_case = TRUE)
    ))
  )
  amr_resident_path <- require_file(
    file.path(DIR_RESULTS, "amr", "resident_amr_profiles.csv"),
    "Script-29 genomic AMR resident profiles"
  )
  amr_residents <- read_csv(amr_resident_path, show_col_types = FALSE)
  require_columns(
    amr_residents,
    c(
      "Participant_id", "n_episodes", "amr_gene_count_informative",
      "amr_mutation_count", "mdfA_detected",
      "any_informative_acquired_amr", "interpretation_scope"
    ),
    "Script-29 genomic AMR resident profiles"
  )
  add_check(
    "genomic_amr", amr_resident_path, "resident_profiles",
    nrow(amr_residents), EXPECTED_RESIDENTS,
    nrow(amr_residents) == EXPECTED_RESIDENTS &&
      n_distinct(amr_residents$Participant_id) == EXPECTED_RESIDENTS
  )

  amr_run_path <- require_file(
    file.path(DIR_RESULTS, "amr", "provenance", "run_manifest.csv"),
    "Script-29 AMR caller run manifest"
  )
  amr_runs <- read_csv(amr_run_path, show_col_types = FALSE)
  require_columns(
    amr_runs,
    c(
      "Participant_id", "tp_lab", "caller", "cache_key", "status",
      "exit_status", "completion_marker", "stderr_path", "command"
    ),
    "Script-29 AMR caller run manifest"
  )
  amr_run_counts <- amr_runs %>%
    filter(status == "complete") %>%
    count(caller)
  required_amr_callers <- c(
    "abricate_resfinder", "amrfinderplus", "resfinder_pointfinder"
  )
  caller_complete <- setNames(amr_run_counts$n, amr_run_counts$caller)
  caller_complete <- caller_complete[required_amr_callers]
  caller_complete[is.na(caller_complete)] <- 0L
  add_check(
    "genomic_amr", amr_run_path, "complete_outputs_by_required_caller",
    paste(names(caller_complete), caller_complete, sep = "=", collapse = ";"),
    paste(required_amr_callers, EXPECTED_EPISODES, sep = "=", collapse = ";"),
    all(caller_complete == EXPECTED_EPISODES)
  )
  add_check(
    "genomic_amr", amr_run_path, "successful_completion_markers",
    sum(file.exists(amr_runs$completion_marker) & amr_runs$exit_status == 0L),
    nrow(amr_runs),
    all(file.exists(amr_runs$completion_marker)) &&
      all(amr_runs$exit_status == 0L)
  )
  amr_coverage_path <- require_file(
    file.path(DIR_RESULTS, "amr", "caller_coverage_summary.csv"),
    "Script-29 AMR caller coverage summary"
  )
  amr_coverage <- read_csv(amr_coverage_path, show_col_types = FALSE)
  add_check(
    "genomic_amr", amr_coverage_path, "caller_coverage_rows",
    paste(amr_coverage$caller, amr_coverage$n_complete, sep = "=", collapse = ";"),
    paste(required_amr_callers, EXPECTED_EPISODES, sep = "=", collapse = ";"),
    nrow(amr_coverage) == length(required_amr_callers) &&
      setequal(amr_coverage$caller, required_amr_callers) &&
      all(amr_coverage$n_complete == EXPECTED_EPISODES)
  )

  amr_transition_path <- require_file(
    file.path(DIR_RESULTS, "amr", "adjacent_pair_amr_profiles_371.csv"),
    "Script-29 adjacent-pair AMR profiles"
  )
  amr_transitions <- read_csv(amr_transition_path, show_col_types = FALSE)
  add_check(
    "genomic_amr", amr_transition_path, "adjacent_pair_profiles",
    nrow(amr_transitions), EXPECTED_ADJACENT_PAIRS,
    nrow(amr_transitions) == EXPECTED_ADJACENT_PAIRS
  )
  focused_amr_path <- require_file(
    file.path(DIR_RESULTS, "amr", "not_uti_to_uti_amr_profiles_9.csv"),
    "Script-29 focused AMR profiles"
  )
  focused_amr <- read_csv(focused_amr_path, show_col_types = FALSE)
  add_check(
    "genomic_amr", focused_amr_path, "Not_UTI_to_UTI_profiles",
    nrow(focused_amr), EXPECTED_NOT_UTI_TO_UTI,
    nrow(focused_amr) == EXPECTED_NOT_UTI_TO_UTI
  )
  amr_inference_path <- require_file(
    file.path(
      DIR_RESULTS, "amr",
      "longitudinal_resident_bootstrap_inference.csv"
    ),
    "Script-29 AMR resident-cluster bootstrap inference"
  )
  amr_inference <- read_csv(amr_inference_path, show_col_types = FALSE)
  require_columns(
    amr_inference,
    c(
      "estimand", "estimate", "ci_lower", "ci_upper",
      "bootstrap_reps_requested", "seed", "n_pairs"
    ),
    "Script-29 AMR resident-cluster bootstrap inference"
  )
  add_check(
    "genomic_amr", amr_inference_path, "bootstrap_contract",
    sprintf(
      "rows=%d; reps=%s; seed=%s; pairs=%s",
      nrow(amr_inference),
      paste(unique(amr_inference$bootstrap_reps_requested), collapse = ";"),
      paste(unique(amr_inference$seed), collapse = ";"),
      paste(unique(amr_inference$n_pairs), collapse = ";")
    ),
    "3 rows; 10000 replicates; seed 20260712; 371 pairs",
    nrow(amr_inference) == 3L &&
      all(amr_inference$bootstrap_reps_requested == 10000L) &&
      all(amr_inference$seed == 20260712L) &&
      all(amr_inference$n_pairs == EXPECTED_ADJACENT_PAIRS)
  )
  amr_phenotype_path <- require_file(
    file.path(
      DIR_RESULTS, "amr",
      "resfinder_predicted_phenotypes_genomic_not_ast.csv"
    ),
    "Script-29 ResFinder genomic-prediction table"
  )
  amr_phenotype <- read_csv(amr_phenotype_path, show_col_types = FALSE)
  add_check(
    "genomic_amr", amr_phenotype_path, "genomic_prediction_episode_profiles",
    n_distinct(paste(amr_phenotype$Participant_id, amr_phenotype$tp_lab)),
    EXPECTED_EPISODES,
    n_distinct(paste(
      amr_phenotype$Participant_id, amr_phenotype$tp_lab
    )) == EXPECTED_EPISODES &&
      all(str_detect(
        amr_phenotype$interpretation_scope,
        regex("not phenotypic AST", ignore_case = TRUE)
      ))
  )
  amr_summary_path <- require_file(
    file.path(DIR_RESULTS, "summary", "table_13_genomic_amr_summary.csv"),
    "Table 13 genomic AMR result summary"
  )
  amr_summary <- read_csv(amr_summary_path, show_col_types = FALSE)
  required_amr_metrics <- c(
    "any_informative_acquired_gene",
    "mdfA_detected_background",
    "known_resistance_mutation",
    "adjacent_informative_acquired_gene_gain_or_loss"
  )
  add_check(
    "genomic_amr", amr_summary_path, "result_summary_metrics",
    sum(required_amr_metrics %in% amr_summary$metric),
    length(required_amr_metrics),
    all(required_amr_metrics %in% amr_summary$metric)
  )
  amr_plot_paths <- file.path(
    DIR_PLOTS, "amr",
    paste0(
      c(
        "most_prevalent_informative_acquired_genes",
        "amr_profile_stability_by_direct_snp_context",
        "caller_concordance_by_determinant_class"
      ),
      ".png"
    )
  )
  add_check(
    "genomic_amr", paste(amr_plot_paths, collapse = ";"),
    "supplementary_plot_files", sum(file.exists(amr_plot_paths)),
    length(amr_plot_paths), all(file.exists(amr_plot_paths))
  )
  stale_report_paths <- file.path(
    DIR_RESULTS, "summary",
    c("report_key_warnings.csv", "report_claim_safety_table.csv")
  )
  stale_report_text <- paste(
    unlist(lapply(
      stale_report_paths[file.exists(stale_report_paths)],
      readLines, warn = FALSE
    )),
    collapse = "\n"
  )
  add_check(
    "genomic_amr", paste(stale_report_paths, collapse = ";"),
    "no_stale_unavailable_AMR_claim",
    str_detect(
      stale_report_text,
      regex(
        "no true AMR|true AMR unavailable|AMR (screening )?unavailable",
        ignore_case = TRUE
      )
    ),
    FALSE,
    !str_detect(
      stale_report_text,
      regex(
        "no true AMR|true AMR unavailable|AMR (screening )?unavailable",
        ignore_case = TRUE
      )
    )
  )
  amr_validation_path <- require_file(
    file.path(DIR_RESULTS, "amr", "validation_checks.csv"),
    "Script-29 AMR validation"
  )
  amr_validation <- read_csv(amr_validation_path, show_col_types = FALSE)
  add_check(
    "genomic_amr", amr_validation_path, "critical_validation_failures",
    sum(!amr_validation$pass & amr_validation$severity == "critical"),
    0L, !any(!amr_validation$pass & amr_validation$severity == "critical")
  )
  amr_marker_path <- require_file(
    file.path(DIR_RESULTS, "amr", "RUN_COMPLETE.txt"),
    "Script-29 AMR completion marker"
  )
  amr_marker_text <- paste(readLines(amr_marker_path, warn = FALSE),
                           collapse = "\n")
  add_check(
    "genomic_amr", amr_marker_path, "completion_marker_contract",
    amr_marker_text,
    "status=complete; episodes=532; residents=161; adjacent_pairs=371; focused_transitions=9",
    all(str_detect(
      amr_marker_text,
      fixed(c(
        "status=complete", "episodes=532", "residents=161",
        "adjacent_pairs=371", "focused_transitions=9"
      ))
    ))
  )
  amr_output_manifest_path <- require_file(
    file.path(
      DIR_RESULTS, "amr", "provenance", "published_output_manifest.csv"
    ),
    "Script-29 AMR published-output manifest"
  )
  amr_output_manifest <- read_csv(
    amr_output_manifest_path, show_col_types = FALSE
  )
  amr_manifest_exists <- file.exists(amr_output_manifest$path)
  amr_manifest_hashes <- rep(NA_character_, nrow(amr_output_manifest))
  amr_manifest_hashes[amr_manifest_exists] <- vapply(
    amr_output_manifest$path[amr_manifest_exists],
    digest, character(1), algo = "sha256", file = TRUE
  )
  add_check(
    "genomic_amr", amr_output_manifest_path,
    "published_outputs_exist_and_match_sha256",
    sum(
      amr_manifest_exists &
        tolower(amr_manifest_hashes) ==
        tolower(amr_output_manifest$sha256),
      na.rm = TRUE
    ),
    nrow(amr_output_manifest),
    all(amr_manifest_exists) &&
      all(
        tolower(amr_manifest_hashes) ==
          tolower(amr_output_manifest$sha256)
      )
  )

  casebook_path <- require_file(file.path(DIR_RESULTS, "mechanism", "not_uti_to_uti_casebook.csv"), "Mechanism casebook")
  casebook <- read_csv(casebook_path, show_col_types = FALSE)
  add_check("casebook", casebook_path, "rows", nrow(casebook), EXPECTED_NOT_UTI_TO_UTI, nrow(casebook) == EXPECTED_NOT_UTI_TO_UTI)
  linkage_column <- intersect(c("has_vf_pair", "wgs_linked"), names(casebook))
  if (length(linkage_column)) {
    linked <- as_pipeline_bool(casebook[[linkage_column[[1]]]])
    add_check("casebook", casebook_path, "linked_rows", sum(linked, na.rm = TRUE), EXPECTED_NOT_UTI_TO_UTI, all(linked))
    add_check("casebook", casebook_path, "missing_linkage_rows", sum(!linked | is.na(linked)), 0L, all(linked))
  } else {
    add_check("casebook", casebook_path, "linkage_column", "missing", "has_vf_pair or wgs_linked", FALSE)
  }

  rq_root <- file.path(DIR_RESULTS, "research_questions")
  rq_marker <- require_file(file.path(rq_root, "RUN_COMPLETE.txt"), "RQ01-RQ10 completion marker")
  rq_status_expected <- file.path(rq_root, sprintf("RQ%02d", 1:10), "analysis_status.csv")
  rq_status <- rq_status_expected[file.exists(rq_status_expected)]
  retired_paths <- list.files(
    rq_root,
    pattern = "(rq11|rq09_11)",
    recursive = TRUE,
    full.names = TRUE,
    include.dirs = TRUE,
    all.files = TRUE,
    ignore.case = TRUE,
    no.. = TRUE
  )
  add_check("research_questions", rq_marker, "numbered_status_files", length(rq_status), 10L, length(rq_status) == 10L)
  add_check("research_questions", rq_root, "retired_question_paths", length(retired_paths), 0L, length(retired_paths) == 0L)
  rq_marker_text <- paste(readLines(rq_marker, warn = FALSE), collapse = "\n")
  rq_marker_scope_ok <- str_detect(rq_marker_text, fixed("RQ01-RQ10")) ||
    str_detect(rq_marker_text, fixed("RQ01--RQ10"))
  add_check("research_questions", rq_marker, "marker_scope", rq_marker_scope_ok, TRUE, rq_marker_scope_ok)

  registry_path <- require_file(
    file.path(DIR_RESULTS, "pipeline", "longcycler_release_claim_registry.json"),
    "Longcycler release claim registry"
  )
  registry <- jsonlite::fromJSON(registry_path, simplifyVector = TRUE)
  add_check("claim_registry", registry_path, "analytical_episodes", registry$analytical_cohort$episodes, EXPECTED_EPISODES, identical(as.integer(registry$analytical_cohort$episodes), EXPECTED_EPISODES))
  add_check("claim_registry", registry_path, "analytical_residents", registry$analytical_cohort$residents, EXPECTED_RESIDENTS, identical(as.integer(registry$analytical_cohort$residents), EXPECTED_RESIDENTS))
  add_check("claim_registry", registry_path, "operational_UTI", registry$analytical_cohort$operational_UTI, EXPECTED_UTI, identical(as.integer(registry$analytical_cohort$operational_UTI), EXPECTED_UTI))
  add_check("claim_registry", registry_path, "operational_Not_UTI", registry$analytical_cohort$operational_Not_UTI, EXPECTED_NOT_UTI, identical(as.integer(registry$analytical_cohort$operational_Not_UTI), EXPECTED_NOT_UTI))
  add_check("claim_registry", registry_path, "vfdb_min_identity_pct", registry$method_contract$vfdb$min_identity_pct, 80L, identical(as.integer(registry$method_contract$vfdb$min_identity_pct), 80L))
  add_check("claim_registry", registry_path, "vfdb_min_coverage_pct", registry$method_contract$vfdb$min_coverage_pct, 80L, identical(as.integer(registry$method_contract$vfdb$min_coverage_pct), 80L))
  add_check(
    "claim_registry", registry_path, "genomic_amr_primary_caller",
    registry$method_contract$genomic_amr$primary_caller,
    "AMRFinderPlus 4.2.7",
    str_detect(
      registry$method_contract$genomic_amr$primary_caller,
      fixed("AMRFinderPlus 4.2.7")
    )
  )
  add_check(
    "claim_registry", registry_path, "genomic_amr_interpretation_scope",
    registry$method_contract$genomic_amr$interpretation,
    "not phenotypic AST",
    str_detect(
      registry$method_contract$genomic_amr$interpretation,
      regex("not phenotypic AST", ignore_case = TRUE)
    )
  )
  add_check("claim_registry", registry_path, "provider_mlst_min_good_targets_pct", registry$method_contract$mlst$provider_min_good_targets_pct, 95L, identical(as.integer(registry$method_contract$mlst$provider_min_good_targets_pct), 95L))
  add_check("claim_registry", registry_path, "direct_pair_snp_threshold", registry$method_contract$direct_pair_evidence$operational_snp_threshold, SAME_STRAIN_SNP_THRESHOLD, identical(as.integer(registry$method_contract$direct_pair_evidence$operational_snp_threshold), SAME_STRAIN_SNP_THRESHOLD))
  add_check("claim_registry", registry_path, "research_question_count", registry$research_questions$count, 10L, identical(as.integer(registry$research_questions$count), 10L))

  registry_sources <- registry$sources
  required_source_roles <- c(
    "cohort", "manifest", "vf", "mlst", "pairs", "transitions",
    "panaroo_manifest", "panaroo_roary_compatibility", "casebook",
    "near_miss", "rq_marker", "vf_method_manifest",
    "amr_episode_profiles", "amr_resident_profiles",
    "amr_adjacent_profiles", "amr_focused_profiles",
    "amr_longitudinal_inference", "amr_caller_coverage",
    "amr_result_summary", "amr_validation", "amr_provenance",
    "plasmidfinder_run", "plasmidfinder_database",
    "plasmidfinder_gene_matrix", "mob_suite_marker", "mob_suite_run",
    "mob_contig_assignments", "mob_predicted_plasmids",
    "mob_episode_profiles", "plasmid_gene_locations",
    "plasmid_episode_mechanism_profiles", "plasmid_adjacent_metrics",
    "plasmid_focused_metrics", "plasmid_location_validation",
    "plasmid_result_summary",
    "project_configuration", "assembly_qc_implementation",
    "provider_mlst_implementation", "vfdb_implementation",
    "genomic_amr_implementation",
    "plasmidfinder_implementation", "mob_reconstruction_implementation",
    "plasmid_localization_implementation",
    "direct_pair_implementation", "core_genome_implementation",
    "pangenome_implementation"
  )
  if (!is.data.frame(registry_sources) ||
      !all(c("role", "path", "sha256") %in% names(registry_sources))) {
    stop("Claim registry sources are not a role/path/SHA-256 table.", call. = FALSE)
  }
  registry_source_paths_exist <- file.exists(registry_sources$path)
  registry_source_hashes <- rep(NA_character_, nrow(registry_sources))
  registry_source_hashes[registry_source_paths_exist] <- vapply(
    registry_sources$path[registry_source_paths_exist],
    digest,
    character(1),
    algo = "sha256",
    file = TRUE
  )
  registry_source_hash_match <- registry_source_paths_exist &
    !is.na(registry_sources$sha256) &
    tolower(registry_source_hashes) == tolower(as.character(registry_sources$sha256))
  add_check("claim_registry", registry_path, "required_source_roles", sum(required_source_roles %in% registry_sources$role), length(required_source_roles), all(required_source_roles %in% registry_sources$role))
  add_check("claim_registry", registry_path, "missing_source_files", sum(!registry_source_paths_exist), 0L, all(registry_source_paths_exist))
  add_check("claim_registry", registry_path, "source_sha256_mismatches", sum(!registry_source_hash_match), 0L, all(registry_source_hash_match))

  registry_plots <- as.character(registry$plot_files)
  registry_plot_exists <- length(registry_plots) > 0L && all(file.exists(registry_plots))
  registry_plot_historical <- str_detect(
    registry_plots,
    regex("/(?:legacy|archive[^/]*|backup[^/]*)/", ignore_case = TRUE)
  ) | str_detect(
    basename(registry_plots),
    regex("(?:ASB.*UTI|UTI.*ASB)", ignore_case = TRUE)
  )
  registry_plot_mtime <- file.info(registry_plots)$mtime
  registry_plot_fresh <- length(registry_plots) > 0L &&
    all(!is.na(registry_plot_mtime) & registry_plot_mtime >= cohort_mtime)
  add_check("claim_registry", registry_path, "current_plot_files", length(registry_plots), ">0 existing files", registry_plot_exists)
  add_check("claim_registry", registry_path, "historical_plot_paths", sum(registry_plot_historical), 0L, !any(registry_plot_historical))
  add_check("claim_registry", registry_path, "plot_files_not_older_than_selected_cohort", sum(!is.na(registry_plot_mtime) & registry_plot_mtime >= cohort_mtime), length(registry_plots), registry_plot_fresh)

  generated_roots <- c(DIR_RESULTS, DIR_PLOTS, DIR_LOGS)
  generated_paths <- unlist(lapply(generated_roots[dir.exists(generated_roots)], function(root) {
    list.files(root, recursive = TRUE, full.names = TRUE, include.dirs = TRUE, all.files = TRUE, no.. = TRUE)
  }), use.names = FALSE)
  bad_path_count <- sum(str_detect(generated_paths, regex("flye", ignore_case = TRUE)))
  add_check("release_hygiene", paste(generated_roots, collapse = ";"), "forbidden_generated_paths", bad_path_count, 0L, bad_path_count == 0L)

  text_globs <- c("*.csv", "*.tsv", "*.txt", "*.md", "*.json", "*.jsonl", "*.ndjson", "*.log", "*.yaml", "*.yml", "*.html", "*.xml", "*.out", "*.err")
  # Long sequence fields can contain the coincidental amino-acid motif FLYE.
  # Match the retired assembler only as a delimited value/path component.
  forbidden_content_pattern <- "(^|[^[:alnum:]])flye([^[:alnum:]]|$)"
  # Results, plots and logs are largely gitignored. The final content gate must
  # scan ignored generated files too or a neutral-named stale artifact can
  # evade the Longcycler-only release check.
  rg_args <- c(
    "-l", "-i", "--hidden", "--no-ignore", "--no-messages",
    "--glob", shQuote("!longcycler_only_pipeline_verification.*")
  )
  for (glob in text_globs) rg_args <- c(rg_args, "--glob", shQuote(glob))
  rg_args <- c(
    rg_args,
    shQuote(forbidden_content_pattern),
    generated_roots[dir.exists(generated_roots)]
  )
  rg_bin <- unname(Sys.which("rg"))
  if (!nzchar(rg_bin)) {
    candidates <- c(
      "/Applications/ChatGPT.app/Contents/Resources/rg",
      "/opt/homebrew/bin/rg",
      "/usr/local/bin/rg"
    )
    rg_bin <- candidates[file.exists(candidates)][1]
  }
  if (is.na(rg_bin) || !nzchar(rg_bin)) stop("ripgrep is required for release content verification.")
  content_hits <- suppressWarnings(system2(rg_bin, rg_args, stdout = TRUE, stderr = FALSE))
  rg_status <- attr(content_hits, "status")
  if (identical(rg_status, 1L)) {
    content_hits <- character()
  } else if (!is.null(rg_status) && !identical(rg_status, 0L)) {
    stop("Release-content ripgrep scan failed with exit status ", rg_status, ".")
  }
  content_hits <- unique(content_hits[nzchar(content_hits)])
  content_hits <- normalizePath(content_hits, winslash = "/", mustWork = FALSE)

  # Repository-wide audit ledgers intentionally retain names of retired
  # artifacts as historical provenance.  Exempt only those ledgers, not the
  # current final manifests/captions/checks or reference-aware analysis tables
  # in the same directory.
  audit_provenance_dir <- file.path(DIR_RESULTS, "figure_audit")
  audit_provenance_files <- c(
    file.path(audit_provenance_dir, "artifact_census.csv"),
    file.path(audit_provenance_dir, "figure_inventory.csv"),
    file.path(audit_provenance_dir, "validation_results.txt"),
    if (dir.exists(audit_provenance_dir)) {
      list.files(
        audit_provenance_dir,
        pattern = "^baseline_run_.*\\.txt$",
        full.names = TRUE,
        recursive = FALSE
      )
    } else {
      character()
    }
  )
  audit_provenance_paths <- unique(normalizePath(
    audit_provenance_files,
    winslash = "/",
    mustWork = FALSE
  ))
  exempted_audit_provenance <- intersect(content_hits, audit_provenance_paths)
  content_hits <- setdiff(content_hits, audit_provenance_paths)
  add_check(
    "release_hygiene",
    paste(audit_provenance_files, collapse = ";"),
    "exempted_audit_provenance_files",
    length(exempted_audit_provenance),
    "informational; exact file allowlist only",
    TRUE
  )
  add_check("release_hygiene", paste(generated_roots, collapse = ";"), "forbidden_text_content_files", length(content_hits), 0L, length(content_hits) == 0L)

  # These were historical active-looking result locations, not explicit
  # archives.  They are intentionally absent from the Longcycler-only release.
  retired_active_roots <- c(
    file.path(DIR_RESULTS, "sensitivity"),
    file.path(DIR_RESULTS, "intermediate")
  )
  retired_active_present <- file.exists(retired_active_roots) | dir.exists(retired_active_roots)
  add_check(
    "release_hygiene",
    paste(retired_active_roots, collapse = ";"),
    "retired_active_result_roots",
    sum(retired_active_present),
    0L,
    !any(retired_active_present)
  )
}

result <- bind_rows(checks)
dir.create(DIR_QC, recursive = TRUE, showWarnings = FALSE)
write_csv(result, out_csv)
failed <- result %>% filter(!pass)
summary_lines <- c(
  "Longcycler-only pipeline verification",
  paste0("Generated: ", Sys.time()),
  paste0("Stage: ", stage),
  paste0("Selected analysis episodes: ", nrow(manifest)),
  paste0("Selected residents: ", n_distinct(manifest$Participant_id)),
  paste0("Derived within-resident pairs: ", expected_pairs),
  paste0("Derived adjacent transitions: ", expected_transitions),
  paste0("Checks passed: ", sum(result$pass), "/", nrow(result)),
  if (nrow(failed)) {
    paste0("FAILED: ", failed$component, " / ", failed$metric, " observed=", failed$observed, " required=", failed$requirement)
  } else {
    "PASS: active release uses only selected Longcycler inputs and endpoints."
  }
)
writeLines(summary_lines, out_txt)

if (nrow(failed)) stop("Longcycler-only verification failed. See ", out_csv, call. = FALSE)
message("Longcycler-only verification passed: ", nrow(result), " checks. Report: ", out_csv)

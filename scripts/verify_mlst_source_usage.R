#!/usr/bin/env Rscript
# ==============================================================================
# verify_mlst_source_usage.R
# ------------------------------------------------------------------------------
# Guardrail that active downstream analyses consume provider-preferred MLST, not
# the local mlst outputs except in explicit provenance/audit scripts.
# ==============================================================================

source("00_config.R")

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
  library(tibble)
})

if (!identical(FILE_MLST_CANONICAL, FILE_MLST_PROVIDER_PREFERRED)) {
  stop("FILE_MLST_CANONICAL is not provider-preferred: ", FILE_MLST_CANONICAL)
}

if (!file.exists(FILE_MLST_CANONICAL)) {
  stop("Provider-preferred MLST file is missing: ", FILE_MLST_CANONICAL)
}

canonical_file <- FILE_ANALYSIS_ASSEMBLY_MANIFEST
if (!file.exists(canonical_file)) {
  stop("Canonical assembly selection is mandatory: ", canonical_file)
}

canonical <- load_analysis_assemblies(canonical_file, require_files = TRUE)
required_canonical <- c("Participant_id", "tp_lab", "Isolate_ID", "selected_canonical", "QC_PASS")
missing_canonical <- setdiff(required_canonical, names(canonical))
if (length(missing_canonical) > 0) {
  stop("Canonical assembly selection lacks required column(s): ", paste(missing_canonical, collapse = ", "))
}
if (!"full_path" %in% names(canonical) && "fasta_path" %in% names(canonical)) canonical$full_path <- canonical$fasta_path
if (!"full_path" %in% names(canonical)) stop("Canonical assembly selection lacks full_path/fasta_path.")
if (!"file_name" %in% names(canonical)) canonical$file_name <- basename(canonical$full_path)

canonical <- canonical %>%
  mutate(
    selected_canonical = as_pipeline_bool(selected_canonical),
    QC_PASS = as_pipeline_bool(QC_PASS),
    active_assembler = str_to_lower(coalesce(
      if ("assembler" %in% names(.)) as.character(assembler) else NA_character_,
      if ("Assembler" %in% names(.)) as.character(Assembler) else NA_character_,
      detect_assembler(coalesce(as.character(full_path), as.character(file_name)))
    )),
    full_path = normalizePath(as.character(full_path), winslash = "/", mustWork = TRUE),
    fasta_sha256 = vapply(full_path, digest::digest, character(1), algo = "sha256", file = TRUE)
  ) %>%
  filter(selected_canonical %in% TRUE, QC_PASS %in% TRUE)

if (nrow(canonical) == 0) stop("Canonical selection contains no selected QC-passing assemblies.")
if (any(is.na(canonical$active_assembler) | canonical$active_assembler != "longcycler")) {
  stop("Selected canonical MLST denominator contains non-Longcycler or unknown assembler provenance.")
}
if (anyDuplicated(canonical[c("Participant_id", "tp_lab")])) {
  stop("Selected Longcycler canonical manifest has duplicate participant-timepoint rows.")
}

mlst <- read_csv(FILE_MLST_CANONICAL, show_col_types = FALSE, progress = FALSE)
required_cols <- c(
  "ST", "ST_source", "ST_provider", "ST_local", "provider_PercGoodTargets",
  "provider_file", "provider_batch_match", "provider_assembler", "full_path", "fasta_sha256",
  "provider_has_classic_7_loci", "ST_numeric_comparable_to_local",
  "ST_provider_below_qc95", "provider_below_qc95_PercGoodTargets",
  "local_mlst_complete", "local_ambiguous_call"
)
missing_cols <- setdiff(required_cols, names(mlst))
if (length(missing_cols) > 0) {
  stop("Provider-preferred MLST file lacks required provenance column(s): ", paste(missing_cols, collapse = ", "))
}

if (nrow(mlst) != nrow(canonical)) {
  stop(
    "Provider-preferred MLST denominator is ", nrow(mlst),
    "; current selected Longcycler manifest has ", nrow(canonical), "."
  )
}
if (!"file_name" %in% names(mlst)) mlst$file_name <- basename(mlst$full_path)

mlst <- mlst %>%
  mutate(
    full_path = normalizePath(as.character(full_path), winslash = "/", mustWork = FALSE),
    active_assembler = str_to_lower(coalesce(
      if ("assembler" %in% names(.)) as.character(assembler) else NA_character_,
      if ("Assembler" %in% names(.)) as.character(Assembler) else NA_character_,
      if ("local_assembler" %in% names(.)) as.character(local_assembler) else NA_character_,
      detect_assembler(coalesce(as.character(full_path), as.character(file_name)))
    )),
    provider_assembler = str_to_lower(str_squish(as.character(provider_assembler)))
  )

if (!setequal(mlst$full_path, canonical$full_path)) {
  stop("Provider-preferred MLST FASTA paths do not exactly match the selected Longcycler manifest.")
}
if (!setequal(
  paste(mlst$full_path, mlst$fasta_sha256, sep = "\n"),
  paste(canonical$full_path, canonical$fasta_sha256, sep = "\n")
)) stop("Provider-preferred MLST FASTA path/SHA-256 pairs do not exactly match the selected Longcycler manifest.")
if (any(is.na(mlst$active_assembler) | mlst$active_assembler != "longcycler")) {
  stop("Provider-preferred MLST contains non-Longcycler or missing active assembly provenance.")
}

allowed_sources <- c(
  "provider_qc95", "local_fallback_provider_missing",
  "local_fallback_provider_conflict", "missing", "missing_provider_conflict"
)
unexpected_sources <- setdiff(unique(na.omit(mlst$ST_source)), allowed_sources)
if (length(unexpected_sources) > 0) {
  stop("Unexpected active ST_source value(s): ", paste(unexpected_sources, collapse = ", "))
}
if (any(is.na(mlst$ST_source) | !nzchar(as.character(mlst$ST_source)))) {
  stop("Provider-preferred MLST contains missing ST_source provenance.")
}

provider_primary <- mlst %>% filter(ST_source == "provider_qc95")
if (nrow(provider_primary) > 0 && any(provider_primary$provider_assembler != "longcycler")) {
  stop("Provider-primary MLST is not tied to the selected Longcycler manifest.")
}
if (nrow(provider_primary) > 0 && any(!(provider_primary$provider_has_classic_7_loci %in% TRUE))) {
  stop("Provider-primary MLST lacks classic seven-locus scheme evidence.")
}
if (nrow(provider_primary) > 0 && any(
  is.na(provider_primary$ST_provider) |
    !nzchar(as.character(provider_primary$ST_provider)) |
    is.na(provider_primary$provider_PercGoodTargets) |
    provider_primary$provider_PercGoodTargets < 95
)) {
  stop("Provider-primary MLST contains a missing ST or a call below the QC95 threshold.")
}
dual_usable <- provider_primary %>% filter(!is.na(.data$ST_local), nzchar(as.character(.data$ST_local)), .data$ST_local != "-")
if (nrow(dual_usable) > 0 && any(as.character(dual_usable$ST_provider) != as.character(dual_usable$ST_local))) {
  stop("Provider/local classic seven-locus ST discordance is present in the active MLST layer.")
}

local_fallback <- mlst %>% filter(str_starts(ST_source, "local_fallback"))
if (nrow(local_fallback) > 0 && any(is.na(local_fallback$ST_local))) {
  stop("A local-fallback MLST row lacks its Longcycler local ST provenance.")
}
if (nrow(local_fallback) > 0 && any(
  !(local_fallback$local_mlst_complete %in% TRUE) |
    local_fallback$local_ambiguous_call %in% TRUE
)) {
  stop("A local-fallback MLST row is incomplete or ambiguous at the classic seven loci.")
}
fallback_below <- local_fallback %>%
  filter(!is.na(.data$ST_provider_below_qc95), nzchar(as.character(.data$ST_provider_below_qc95)))
if (nrow(fallback_below) > 0 && any(
  as.character(fallback_below$ST_local) != as.character(fallback_below$ST_provider_below_qc95)
)) {
  stop("A local fallback disagrees with its retained below-QC95 provider ST evidence.")
}
missing_rows <- mlst %>% filter(.data$ST_source %in% c("missing", "missing_provider_conflict"))
if (nrow(missing_rows) > 0 && any(!is.na(missing_rows$ST) & nzchar(as.character(missing_rows$ST)))) {
  stop("An MLST row labelled missing contains an accepted active ST.")
}
called <- mlst %>% filter(!is.na(.data$ST), nzchar(as.character(.data$ST)), .data$ST != "-")
if (nrow(called) > 0 && any(!(called$ST_numeric_comparable_to_local %in% TRUE))) {
  stop("A usable active MLST call is not certified as classic seven-locus comparable.")
}

source_counts <- mlst %>% count(ST_source, name = "n")
provider_n <- source_counts$n[match("provider_qc95", source_counts$ST_source)]
provider_n <- ifelse(is.na(provider_n), 0L, provider_n)
if (sum(source_counts$n) != nrow(canonical)) {
  stop("MLST source assignments do not sum to the current Longcycler manifest denominator.")
}

scan_roots <- c(".", "R", "scripts")
scan_files <- unique(unlist(lapply(scan_roots, function(root) {
  if (!dir.exists(root)) return(character())
  list.files(root, pattern = "\\.(R|r)$", recursive = TRUE, full.names = TRUE)
})))

allow_files <- normalizePath(c(
  "00_config.R",
  "06_MLST.R",
  "scripts/compare_mlst_sources.R",
  "R/provider_mlst_integration.R",
  "scripts/verify_mlst_source_usage.R"
), winslash = "/", mustWork = FALSE)

scan_files <- scan_files[!str_detect(scan_files, "/legacy/|^legacy/|/archive/|^archive/")]
scan_files_norm <- normalizePath(scan_files, winslash = "/", mustWork = FALSE)
scan_files <- scan_files[!scan_files_norm %in% allow_files]

forbidden_read_patterns <- c(
  "read[a-zA-Z_:.]*\\s*\\([^\\n)]*FILE_MLST_ALL",
  "read[a-zA-Z_:.]*\\s*\\([^\\n)]*FILE_MLST_ISOLATE_EXPLORATORY",
  "read[a-zA-Z_:.]*\\s*\\([^\\n)]*mlst_all\\.tsv",
  "read[a-zA-Z_:.]*\\s*\\([^\\n)]*mlst_with_meta\\.csv",
  "read_if_exists\\s*\\([^\\n)]*mlst_with_meta\\.csv"
)

violations <- list()
for (path in scan_files) {
  lines <- readLines(path, warn = FALSE)
  active_lines <- lines[!str_detect(str_trim(lines), "^#")]
  hit <- Reduce(`|`, lapply(forbidden_read_patterns, function(pattern) {
    str_detect(active_lines, regex(pattern, ignore_case = TRUE))
  }))
  if (any(hit, na.rm = TRUE)) {
    violations[[path]] <- tibble(
      file = path,
      line = which(!str_detect(str_trim(lines), "^#"))[which(hit)],
      text = active_lines[hit]
    )
  }
}

violations <- bind_rows(violations)
if (nrow(violations) > 0) {
  write_csv(violations, file.path(DIR_QC, "mlst_source_usage_violations.csv"))
  stop(
    "Found active script(s) reading local MLST directly. See ",
    file.path(DIR_QC, "mlst_source_usage_violations.csv")
  )
}

write_csv(source_counts, file.path(DIR_QC, "mlst_active_source_counts.csv"))
msg(
  "Longcycler-only MLST source usage verified: canonical=%d, provider_qc95=%d, local_fallback=%d.",
  nrow(canonical), provider_n, nrow(local_fallback)
)

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

mlst <- read_csv(FILE_MLST_CANONICAL, show_col_types = FALSE, progress = FALSE)
required_cols <- c(
  "ST", "ST_source", "ST_provider", "ST_local", "provider_PercGoodTargets",
  "provider_file", "provider_batch_match", "provider_assembler"
)
missing_cols <- setdiff(required_cols, names(mlst))
if (length(missing_cols) > 0) {
  stop("Provider-preferred MLST file lacks required provenance column(s): ", paste(missing_cols, collapse = ", "))
}

if (nrow(mlst) != 556L) {
  stop("Provider-preferred MLST denominator is ", nrow(mlst), "; expected 556 canonical isolates.")
}

source_counts <- mlst %>% count(ST_source, name = "n")
provider_n <- source_counts$n[match("provider_qc95", source_counts$ST_source)]
provider_n <- ifelse(is.na(provider_n), 0L, provider_n)
if (provider_n < 527L) {
  stop("Provider QC95 ST coverage is ", provider_n, "; expected at least 527.")
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
msg("MLST source usage verified. Provider-preferred ST is active; provider_qc95 rows: %d.", provider_n)

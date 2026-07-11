# ==============================================================================
# 00_input_snapshot.R
# ------------------------------------------------------------------------------
# GOAL:
#   Create debug snapshots (head/dims) of key input files.
#
# METHOD:
#   1. Read first N rows/cols of specified CSV/TSV files.
#   2. Write dimensions and head to debug directory.
#
# INPUTS:
#   - 00_config.R
#   - results/vf_pa_all.csv
#   - results/mlst/mlst_provider_preferred.csv
#
# OUTPUTS:
#   - results/within_person/debug/input_snapshot/
# ==============================================================================

source("00_config.R")
suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
})

DIR_DEBUG_SNAP <- file.path(DIR_RESULTS, "within_person", "debug", "input_snapshot")
ensure_dir(DIR_DEBUG_SNAP)
D <- DIR_DEBUG_SNAP

dump <- function(p, n = 200, m = 80) {
  if (!file.exists(p)) {
    return(invisible())
  }
  x <- tryCatch(
    readr::read_delim(p, delim = if (grepl("\\.tsv$", p)) "\t" else ",", show_col_types = FALSE),
    error = function(e) NULL
  )
  if (is.null(x)) {
    return(invisible())
  }
  readr::write_csv(tibble(file = p, nrow = nrow(x), ncol = ncol(x)), file.path(D, "dims.csv"), append = TRUE)
  readr::write_lines(paste(names(x), collapse = ","), file.path(D, paste0(basename(p), ".names.txt")))
  # head rows and first m cols to keep tiny
  cols <- seq_len(min(ncol(x), m))
  readr::write_csv(head(x[, cols, drop = FALSE], n), file.path(D, paste0(basename(p), ".head.csv")))
}

dump(FILE_VF_PA)
dump("results/mlst/mlst_matrix.csv")
dump(FILE_MLST_CANONICAL)
dump(FILE_MLST_PROVIDER_PREFERRED_ALL)
dump("results/plasmidfinder_presence_absence.csv")
dump(FILE_STATUS_MAP)

# list nucmer reports (paths only)
rp <- list.files("results/nucmer", pattern = "run_dd.report$", recursive = TRUE, full.names = TRUE)
readr::write_lines(rp, file.path(D, "nucmer_reports.txt"))
cat("Wrote previews under:", normalizePath(D), "\n")

#!/usr/bin/env Rscript
# Regenerate Prokka GFFs for the current Panaroo missing-GFF manifest.
#
# Normal full-pipeline runs do this automatically inside 12c_panaroo.R. This
# wrapper remains available for manual recovery or dry-run inspection.
#
# Usage:
#   GFF_REGEN_CORES=6 Rscript scripts/regenerate_missing_gffs.R
#   GFF_REGEN_DRY_RUN=1 Rscript scripts/regenerate_missing_gffs.R

source("00_config.R")

suppressPackageStartupMessages({
  library(readr)
})

log_line <- function(...) {
  cat(sprintf("[%s] %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), sprintf(...)))
  flush.console()
}

missing_file <- file.path(DIR_WGS, "pan", "missing_gffs.csv")
if (!file.exists(missing_file)) {
  stop("Missing-GFF manifest not found: ", missing_file, "\nRun Rscript 12c_panaroo.R first.")
}

missing_df <- read_csv(missing_file, show_col_types = FALSE)

log_line("Missing-GFF manifest: %s", missing_file)
log_line("Rows in manifest: %d", nrow(missing_df))

result <- regenerate_missing_gffs(
  missing_df,
  summary_file = file.path(DIR_WGS, "pan", "regenerate_missing_gffs_summary.csv"),
  logger = function(x) log_line("%s", x)
)

if (isTRUE(result$dry_run) && nrow(result$jobs) > 0) {
  print(
    result$jobs |>
      dplyr::select(dplyr::any_of(c(
        "Assembly_Base_ID", "Assembler", "Participant_id", "tp_lab",
        "fasta_path", "gff_path", "log_path"
      ))) |>
      utils::head(20),
    n = 20
  )
}

log_line("Summary: %s", result$summary_file)

if (!identical(as.integer(result$status), 0L)) {
  quit(save = "no", status = 1)
}

quit(save = "no", status = 0)

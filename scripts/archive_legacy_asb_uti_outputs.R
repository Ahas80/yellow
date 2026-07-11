#!/usr/bin/env Rscript
# ==============================================================================
# Archive stale ASB-vs-UTI generated outputs
# ------------------------------------------------------------------------------
# The primary pipeline contrast is UTI vs Not_UTI. Old generated files whose
# filenames still advertise ASB-vs-UTI are moved out of current results/plots
# folders so they cannot be mistaken for current figures or tables.
# ==============================================================================

source("00_config.R")

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
  library(tibble)
})

msg("Starting scripts/archive_legacy_asb_uti_outputs.R")

legacy_root_results <- file.path(DIR_RESULTS, "legacy", "old_asb_uti_outputs")
legacy_root_plots <- file.path(DIR_PLOTS, "legacy", "old_asb_uti_outputs")
ensure_dir(legacy_root_results)
ensure_dir(legacy_root_plots)

legacy_name <- regex(
  "(UTI_vs_ASB|ASB_vs_UTI|asb_uti|vf_asb_uti|table_10_asb_uti|diff_.*UTI_vs_ASB|permanova_UTI_vs_ASB)",
  ignore_case = TRUE
)

collect_candidates <- function(root) {
  if (!dir.exists(root)) return(character())
  all_files <- list.files(root, recursive = TRUE, full.names = TRUE, all.files = FALSE, no.. = TRUE)
  all_files <- all_files[file.exists(all_files) & !dir.exists(all_files)]
  all_files <- all_files[str_detect(basename(all_files), legacy_name)]
  all_files[!str_detect(all_files, regex("/legacy/", ignore_case = TRUE))]
}

move_one <- function(path, root, archive_root) {
  rel <- sub(paste0("^", gsub("([\\W])", "\\\\\\1", normalizePath(root, winslash = "/", mustWork = FALSE)), "/?"),
             "", normalizePath(path, winslash = "/", mustWork = FALSE))
  dest <- file.path(archive_root, rel)
  ensure_dir(dirname(dest))
  ok <- suppressWarnings(file.rename(path, dest))
  if (!ok) {
    ok <- file.copy(path, dest, overwrite = TRUE)
    if (ok) unlink(path)
  }
  tibble(source_path = path, archived_path = dest, archived = ok)
}

result_files <- collect_candidates(DIR_RESULTS)
plot_files <- collect_candidates(DIR_PLOTS)

manifest <- bind_rows(
  bind_rows(lapply(result_files, move_one, root = DIR_RESULTS, archive_root = legacy_root_results)),
  bind_rows(lapply(plot_files, move_one, root = DIR_PLOTS, archive_root = legacy_root_plots))
) %>%
  mutate(archived_at = format(Sys.time()))

if (nrow(manifest) == 0) {
  manifest <- tibble(
    source_path = character(),
    archived_path = character(),
    archived = logical(),
    archived_at = character()
  )
}

manifest_path <- file.path(DIR_QC, "legacy_status_archive_manifest.csv")
write_csv(manifest, manifest_path)
msg("Archived %d stale ASB-vs-UTI generated output(s); manifest: %s",
    sum(manifest$archived %in% TRUE), manifest_path)

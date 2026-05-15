#!/usr/bin/env Rscript
# ==============================================================================
# 12c_panaroo.R
# ------------------------------------------------------------------------------
# Prepare Prokka GFF annotations for the current canonical QC PASS assembly set,
# then run Panaroo when the pangenome outputs are missing, stale, or explicitly
# forced. Existing valid GFFs are reused; only missing/empty GFFs are generated.
# ==============================================================================

source("00_config.R")
source("R/wgs_helpers.R")

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
})

DIR_PAN <- file.path(DIR_WGS, "pan")
ensure_dir(DIR_PAN)

log_info("Starting 12c_panaroo.R")

output_file <- file.path(DIR_PAN, "gene_presence_absence.csv")
manifest_file <- file.path(DIR_PAN, "panaroo_input_manifest.csv")
hash_file <- file.path(DIR_PAN, "panaroo_input_manifest.hash")
missing_report <- file.path(DIR_PAN, "missing_gffs.csv")
stale_report <- file.path(DIR_PAN, "panaroo_staleness_report.txt")
regen_summary_file <- file.path(DIR_PAN, "regenerate_missing_gffs_summary.csv")
force_panaroo <- identical(Sys.getenv("FORCE_RERUN_PANAROO", "0"), "1")
dry_run_requested <- identical(Sys.getenv("GFF_REGEN_DRY_RUN", "0"), "1")
gff_dirs <- c(DIR_PROKKA, DIR_PROKKA_SLIM)
input_source <- "dynamic assembly/GFF inventory from DIR_FASTAS, assembly_metadata.csv, and canonical_assembly_selection.csv"
current_inventory <- NULL

write_panaroo_report <- function(manifest,
                                 current_hash,
                                 previous_hash,
                                 status_line,
                                 extra = character()) {
  if (is.null(manifest)) manifest <- tibble()
  gff_available <- if ("gff_available" %in% names(manifest)) sum(manifest$gff_available %in% TRUE) else 0L
  gff_missing <- if ("gff_available" %in% names(manifest)) sum(!(manifest$gff_available %in% TRUE)) else NA_integer_
  coverage <- if (nrow(manifest) > 0) 100 * gff_available / nrow(manifest) else 0
  inventory_lines <- if (!is.null(current_inventory)) {
    summary_lookup <- stats::setNames(current_inventory$summary$value, current_inventory$summary$metric)
    c(
      sprintf("All metadata-linked FASTAs: %s", summary_lookup[["metadata_linked_fastas"]] %||% "not available"),
      sprintf("All metadata-linked GFFs available: %s", summary_lookup[["metadata_linked_gffs_available"]] %||% "not available"),
      sprintf("All metadata-linked GFFs missing (warning only): %s", summary_lookup[["metadata_linked_gffs_missing_warning_only"]] %||% "not available"),
      sprintf("Unexpected/unlinked candidate FASTAs: %s", summary_lookup[["unexpected_unlinked_candidate_fastas"]] %||% "not available"),
      sprintf("New candidate FASTAs since metadata scan: %s", summary_lookup[["new_candidate_fastas_since_metadata_scan"]] %||% "not available"),
      sprintf("GFF inventory summary: %s", current_inventory$paths$summary),
      sprintf("All-metadata GFF inventory: %s", current_inventory$paths$all_metadata),
      sprintf("Panaroo GFF inventory: %s", current_inventory$paths$panaroo)
    )
  } else {
    character()
  }

  writeLines(
    c(
      "Panaroo input staleness/completeness report",
      sprintf("Generated: %s", format(Sys.time())),
      sprintf("Input source: %s", input_source),
      sprintf("Expected canonical QC PASS assemblies: %d", nrow(manifest)),
      sprintf("GFFs available: %d", gff_available),
      sprintf("Missing GFFs: %d", gff_missing),
      sprintf("GFF coverage: %.1f%%", coverage),
      sprintf("Current manifest hash: %s", current_hash %||% "<not calculated>"),
      sprintf("Previous manifest hash: %s", ifelse(is.na(previous_hash), "<none>", previous_hash)),
      sprintf("FORCE_RERUN_PANAROO: %s", Sys.getenv("FORCE_RERUN_PANAROO", "0")),
      sprintf("GFF_REGEN_CORES: %s", {
        x <- Sys.getenv("GFF_REGEN_CORES", unset = "auto")
        ifelse(identical(x, ""), "auto", x)
      }),
      sprintf("GFF regeneration summary: %s", regen_summary_file),
      inventory_lines,
      status_line,
      extra
    ),
    stale_report
  )
}

refresh_manifest <- function() {
  inventory <- build_assembly_gff_inventory(gff_dirs = gff_dirs, write_outputs = TRUE)
  current_inventory <<- inventory
  manifest <- inventory$panaroo_inventory
  write_pipeline_gff_manifest(manifest, manifest_file, missing_report)
  manifest
}

manifest <- refresh_manifest()
log_info("Metadata-linked FASTAs in current inventory: ", nrow(current_inventory$metadata_linked_inventory))
log_info("Panaroo-eligible canonical QC PASS assemblies: ", nrow(manifest))

current_hash <- hash_input_manifest(manifest)
previous_hash <- if (file.exists(hash_file)) readLines(hash_file, warn = FALSE)[1] else NA_character_

if (current_inventory$metadata_stale || current_inventory$qc_stale) {
  stale_messages <- c(current_inventory$metadata_stale_messages, current_inventory$qc_stale_messages)
  write_panaroo_report(
    manifest,
    current_hash,
    previous_hash,
    "Status: RED - FASTA metadata/QC inventory is stale relative to current inputs.",
    c(
      stale_messages,
      "Next step: rerun 00_make_assembly_metadata.r and 12a_wgs_qc.R before 12c_panaroo.R."
    )
  )
  stop(
    "FASTA metadata/QC inventory is stale. ",
    paste(stale_messages, collapse = " "),
    " Rerun 00_make_assembly_metadata.r and 12a_wgs_qc.R before 12c_panaroo.R."
  )
}

missing_before <- manifest %>% filter(!(gff_available %in% TRUE))
log_info("GFFs available before regeneration: ", sum(manifest$gff_available %in% TRUE), "/", nrow(manifest))

# Always refresh the regeneration summary so downstream QA reports the current
# state, even when there is nothing to do.
regen_result <- regenerate_missing_gffs(
  missing_before,
  summary_file = regen_summary_file,
  logger = function(x) log_info(x)
)

if (dry_run_requested) {
  manifest <- refresh_manifest()
  current_hash <- hash_input_manifest(manifest)
  previous_hash <- if (file.exists(hash_file)) readLines(hash_file, warn = FALSE)[1] else NA_character_
  write_panaroo_report(
    manifest,
    current_hash,
    previous_hash,
    "Status: DRY-RUN - GFF inventory was rebuilt but Prokka/Panaroo were not run.",
    "Next step: rerun without GFF_REGEN_DRY_RUN=1 to generate missing GFFs and/or refresh Panaroo."
  )
  log_warn("GFF_REGEN_DRY_RUN=1; not running Panaroo.")
  quit(save = "no", status = 0)
}

manifest <- refresh_manifest()
missing_after <- manifest %>% filter(!(gff_available %in% TRUE))
log_info("GFFs available after regeneration: ", sum(manifest$gff_available %in% TRUE), "/", nrow(manifest))

current_hash <- hash_input_manifest(manifest)
previous_hash <- if (file.exists(hash_file)) readLines(hash_file, warn = FALSE)[1] else NA_character_
outputs_exist <- file.exists(output_file)
hash_matches <- outputs_exist && identical(current_hash, previous_hash)

append_denominator_summary(
  manifest,
  "12c_panaroo.R",
  "panaroo_input_manifest",
  "participant_timepoint",
  manifest_file,
  sprintf(
    "GFF coverage %.1f%%; missing GFFs %d; output is current only when the Panaroo manifest hash matches",
    ifelse(nrow(manifest) > 0, 100 * sum(manifest$gff_available %in% TRUE) / nrow(manifest), 0),
    nrow(missing_after)
  )
)

if (nrow(missing_after) > 0) {
  write_panaroo_report(
    manifest,
    current_hash,
    previous_hash,
    "Status: RED - Panaroo input incomplete because required GFFs are still missing.",
    sprintf("Next step: review %s and %s.", missing_report, regen_summary_file)
  )
  stop(
    "Panaroo input incomplete: ", nrow(missing_after), " required GFF(s) are still missing. ",
    "Review ", missing_report, " and ", regen_summary_file, "."
  )
}

if (hash_matches && !force_panaroo) {
  write_panaroo_report(
    manifest,
    current_hash,
    previous_hash,
    "Status: GREEN - Panaroo output matches the current input manifest."
  )
  log_info("Panaroo output exists and input manifest hash matches. Skipping run.")
  quit(save = "no", status = 0)
}

if (sum(manifest$gff_available %in% TRUE) < 2) {
  stop("Not enough GFFs found for pangenome analysis.")
}

panaroo_bin <- find_prokka_bin("panaroo")
if (identical(panaroo_bin, "")) stop("Panaroo is required for this module.")

valid_gffs <- manifest %>% filter(gff_available) %>% pull(gff_path)
run_reason <- if (!outputs_exist) {
  "Panaroo output is missing."
} else if (force_panaroo && hash_matches) {
  "FORCE_RERUN_PANAROO=1 requested a refresh despite a matching manifest."
} else {
  "Panaroo output is stale relative to the current GFF input manifest."
}

write_panaroo_report(
  manifest,
  current_hash,
  previous_hash,
  "Status: RUNNING - Panaroo is being generated/refreshed for current inputs.",
  run_reason
)

log_info("Running Panaroo on ", length(valid_gffs), " GFF files. ", run_reason)
run_started <- Sys.time()
args <- c(
  "-i", valid_gffs,
  "-o", DIR_PAN,
  "--clean-mode", "strict",
  "--remove-invalid-genes",
  "-t", Sys.getenv("PANAROO_THREADS", "1")
)

res <- system2(panaroo_bin, args = args)
res <- as.integer(if (is.null(res)) 0L else res)
if (res != 0) {
  write_panaroo_report(
    manifest,
    current_hash,
    previous_hash,
    sprintf("Status: RED - Panaroo failed with exit code %d.", res),
    "Next step: inspect the Panaroo console/log output above."
  )
  log_error("Panaroo failed with exit code ", res)
  stop("Panaroo execution failed.")
}

if (!file.exists(output_file)) {
  write_panaroo_report(
    manifest,
    current_hash,
    previous_hash,
    "Status: RED - Panaroo completed but gene_presence_absence.csv is absent.",
    sprintf("Expected output: %s", output_file)
  )
  stop("Panaroo completed but expected output is absent: ", output_file)
}

if (file.info(output_file)$mtime < run_started) {
  write_panaroo_report(
    manifest,
    current_hash,
    previous_hash,
    "Status: RED - Panaroo completed but gene_presence_absence.csv was not refreshed.",
    sprintf("Expected refreshed output: %s", output_file)
  )
  stop("Panaroo completed but output was not refreshed: ", output_file)
}

writeLines(current_hash, hash_file)
write_panaroo_report(
  manifest,
  current_hash,
  previous_hash,
  "Status: GREEN - Panaroo output was generated/refreshed for current inputs."
)

log_info("12c_panaroo.R complete.")

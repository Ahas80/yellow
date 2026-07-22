#!/usr/bin/env Rscript
# ==============================================================================
# 12c_panaroo.R
# ------------------------------------------------------------------------------
# Prepare Prokka GFF annotations for the validated Longcycler-only analysis
# manifest, then run Panaroo when required pangenome outputs are missing, stale,
# or explicitly forced. Existing valid GFFs are reused; only missing/empty GFFs
# are generated.
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
command_log_file <- file.path(DIR_PAN, "panaroo_command.log")
force_panaroo <- identical(Sys.getenv("FORCE_RERUN_PANAROO", "0"), "1")
dry_run_requested <- identical(Sys.getenv("GFF_REGEN_DRY_RUN", "0"), "1")
resume_output_validation <- identical(
  Sys.getenv("RESUME_PANAROO_OUTPUT_VALIDATION", "0"),
  "1"
)
# Preserve the report that initiated an interrupted Panaroo run before the
# inventory refresh below rewrites any generated manifests.  Resume validation
# uses this recorded start time, not the mtime of an identical regenerated CSV.
resume_run_report <- if (resume_output_validation && file.exists(stale_report)) {
  readLines(stale_report, warn = FALSE)
} else {
  character()
}
gff_dirs <- c(DIR_PROKKA, DIR_PROKKA_SLIM)
input_source <- FILE_ANALYSIS_ASSEMBLY_MANIFEST
current_inventory <- NULL

required_outputs <- file.path(
  DIR_PAN,
  c(
    "gene_presence_absence.csv",
    "gene_presence_absence_roary.csv",
    "gene_presence_absence.Rtab",
    "pan_genome_reference.fa",
    "final_graph.gml",
    "gene_data.csv",
    "summary_statistics.txt",
    "panaroo_command.log"
  )
)

files_complete <- function(paths) {
  all(file.exists(paths) & !is.na(file.size(paths)) & file.size(paths) > 0)
}

panaroo_hash_cols <- c(
  "Assembly_ID", "Participant_id", "tp_lab", "Assembler",
  "selection_policy_version", "fasta_path", "fasta_sha256",
  "gff_path", "gff_sha256"
)

panaroo_input_hash <- function(manifest) hash_input_manifest(manifest, panaroo_hash_cols)

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
      sprintf("Required outputs complete: %s", files_complete(required_outputs)),
      sprintf("Required outputs: %s", paste(required_outputs, collapse = "; ")),
      sprintf("FORCE_RERUN_PANAROO: %s", Sys.getenv("FORCE_RERUN_PANAROO", "0")),
      sprintf(
        "RESUME_PANAROO_OUTPUT_VALIDATION: %s",
        Sys.getenv("RESUME_PANAROO_OUTPUT_VALIDATION", "0")
      ),
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
  manifest <- inventory$panaroo_inventory |>
    dplyr::mutate(
      fasta_path = normalizePath(.data$fasta_path, winslash = "/", mustWork = FALSE),
      gff_path = ifelse(
        is.na(.data$gff_path),
        NA_character_,
        normalizePath(.data$gff_path, winslash = "/", mustWork = FALSE)
      )
    ) |>
    add_file_content_sha256("fasta_path", "fasta_sha256") |>
    add_file_content_sha256("gff_path", "gff_sha256")
  if (any(is.na(manifest$fasta_sha256) | nchar(manifest$fasta_sha256) != 64L)) {
    stop("Failed to create SHA-256 provenance for every Panaroo FASTA input.")
  }
  bad_gff_hash <- manifest$gff_available %in% TRUE &
    (is.na(manifest$gff_sha256) | nchar(manifest$gff_sha256) != 64L)
  if (any(bad_gff_hash)) stop("Failed to create SHA-256 provenance for available Panaroo GFF inputs.")
  write_pipeline_gff_manifest(manifest, manifest_file, missing_report)
  manifest
}

manifest <- refresh_manifest()
log_info("Metadata-linked FASTAs in current inventory: ", nrow(current_inventory$metadata_linked_inventory))
log_info("Panaroo-eligible Longcycler analysis-manifest assemblies: ", nrow(manifest))

current_hash <- panaroo_input_hash(manifest)
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
  current_hash <- panaroo_input_hash(manifest)
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

current_hash <- panaroo_input_hash(manifest)
previous_hash <- if (file.exists(hash_file)) readLines(hash_file, warn = FALSE)[1] else NA_character_
outputs_complete <- files_complete(required_outputs)
hash_matches <- outputs_complete && identical(current_hash, previous_hash)
valid_gffs <- manifest %>% filter(gff_available) %>% pull(gff_path)

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

if (resume_output_validation) {
  if (force_panaroo) {
    stop(
      "RESUME_PANAROO_OUTPUT_VALIDATION=1 cannot be combined with FORCE_RERUN_PANAROO=1."
    )
  }
  if (!outputs_complete) {
    missing_outputs <- required_outputs[
      !file.exists(required_outputs) |
        is.na(file.size(required_outputs)) |
        file.size(required_outputs) <= 0
    ]
    stop(
      "Cannot resume Panaroo output validation because required outputs are missing/empty: ",
      paste(missing_outputs, collapse = ", ")
    )
  }

  command_lines <- readLines(command_log_file, warn = FALSE)
  logged_hash_lines <- grep(
    "^Input manifest SHA-256: [[:xdigit:]]{64}$",
    command_lines,
    value = TRUE
  )
  if (length(logged_hash_lines) != 1L) {
    stop("Panaroo command log does not contain exactly one SHA-256 input-manifest line.")
  }
  logged_hash <- sub("^Input manifest SHA-256: ", "", logged_hash_lines)
  if (!identical(tolower(logged_hash), tolower(current_hash))) {
    stop(
      "Panaroo command-log input hash does not match the current exact FASTA/GFF manifest."
    )
  }
  logged_exit_status <- grep("^Exit status: [0-9]+$", command_lines, value = TRUE)
  if (length(logged_exit_status) != 1L || !identical(logged_exit_status, "Exit status: 0")) {
    stop("Panaroo command log does not prove a successful exit status of zero.")
  }

  command_line_index <- match(logged_hash_lines, command_lines) + 2L
  if (length(command_line_index) != 1L || command_line_index > length(command_lines)) {
    stop("Panaroo command line is absent from the saved command log.")
  }
  command_tokens <- scan(
    text = command_lines[[command_line_index]],
    what = character(),
    quiet = TRUE
  )
  input_flag <- match("-i", command_tokens)
  output_flag <- match("-o", command_tokens)
  if (
    is.na(input_flag) || is.na(output_flag) ||
      input_flag >= output_flag - 1L || output_flag >= length(command_tokens)
  ) {
    stop("Saved Panaroo command does not contain a valid -i ... -o argument block.")
  }
  logged_gffs <- normalizePath(
    command_tokens[seq.int(input_flag + 1L, output_flag - 1L)],
    winslash = "/",
    mustWork = FALSE
  )
  expected_gffs <- normalizePath(valid_gffs, winslash = "/", mustWork = FALSE)
  logged_output_dir <- normalizePath(
    command_tokens[[output_flag + 1L]],
    winslash = "/",
    mustWork = FALSE
  )
  if (!identical(logged_gffs, expected_gffs)) {
    stop(
      "Saved Panaroo command GFF inputs do not exactly match the current ordered 532-GFF manifest."
    )
  }
  if (!identical(logged_output_dir, normalizePath(DIR_PAN, winslash = "/", mustWork = FALSE))) {
    stop("Saved Panaroo command output directory does not match the active Panaroo directory.")
  }
  if (
    !all(c("--clean-mode", "strict", "--remove-invalid-genes") %in% command_tokens)
  ) {
    stop("Saved Panaroo command does not match the required strict cleaning policy.")
  }

  running_status <- grep("^Status: RUNNING - Panaroo", resume_run_report, value = TRUE)
  run_started_lines <- grep(
    "^Generated: [0-9]{4}-[0-9]{2}-[0-9]{2} [0-9]{2}:[0-9]{2}:[0-9]{2}$",
    resume_run_report,
    value = TRUE
  )
  if (length(running_status) != 1L || length(run_started_lines) != 1L) {
    stop(
      "The preserved Panaroo report does not prove exactly one interrupted RUNNING invocation."
    )
  }
  recorded_run_start <- as.POSIXct(
    sub("^Generated: ", "", run_started_lines),
    format = "%Y-%m-%d %H:%M:%S",
    tz = Sys.timezone()
  )
  output_mtimes <- file.info(required_outputs)$mtime
  if (
    is.na(recorded_run_start) ||
      any(is.na(output_mtimes) | output_mtimes < recorded_run_start - 2)
  ) {
    stop("One or more Panaroo outputs predate the preserved start of the saved run.")
  }

  roary_metadata_columns <- c(
    "Gene", "Non-unique Gene name", "Annotation", "No. isolates",
    "No. sequences", "Avg sequences per isolate", "Genome Fragment",
    "Order within Fragment", "Accessory Fragment",
    "Accessory Order with Fragment", "QC", "Min group size nuc",
    "Max group size nuc", "Avg group size nuc"
  )
  roary_columns <- names(read_csv(
    file.path(DIR_PAN, "gene_presence_absence_roary.csv"),
    n_max = 0,
    show_col_types = FALSE,
    name_repair = "minimal"
  ))
  roary_samples <- setdiff(roary_columns, roary_metadata_columns)
  expected_samples <- sub("\\.gff$", "", basename(expected_gffs), ignore.case = TRUE)
  if (
    length(roary_samples) != length(expected_samples) ||
      anyDuplicated(roary_samples) || anyDuplicated(expected_samples) ||
      !setequal(roary_samples, expected_samples)
  ) {
    stop(
      "Panaroo Roary-compatible sample columns do not exactly match the current 532-GFF manifest."
    )
  }

  writeLines(current_hash, hash_file)
  write_panaroo_report(
    manifest,
    current_hash,
    previous_hash,
    paste0(
      "Status: GREEN - Existing Panaroo outputs were resumed after exact command-log, ",
      "SHA-256 manifest, output-freshness, and 532-sample-column validation."
    ),
    c(
      sprintf("Validated saved GFF inputs: %d", length(logged_gffs)),
      sprintf("Validated Roary sample columns: %d", length(roary_samples)),
      sprintf("Preserved Panaroo run start: %s", format(recorded_run_start))
    )
  )
  log_info(
    "Validated and certified existing Panaroo outputs for ",
    length(logged_gffs),
    " exact Longcycler GFF inputs; Panaroo was not rerun."
  )
  quit(save = "no", status = 0)
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

run_reason <- if (!outputs_complete) {
  "One or more required Panaroo outputs are missing or empty."
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
# Remove the previous completion marker and required outputs before invoking
# Panaroo. A failed/interrupted run must never leave a reusable GREEN state.
unlink(c(hash_file, required_outputs), force = TRUE)
args <- c(
  "-i", valid_gffs,
  "-o", DIR_PAN,
  "--clean-mode", "strict",
  "--remove-invalid-genes",
  "-t", Sys.getenv("PANAROO_THREADS", "1")
)

panaroo_capture <- tempfile("panaroo-command-")
on.exit(unlink(panaroo_capture), add = TRUE)
res <- system2(panaroo_bin, args = args, stdout = panaroo_capture, stderr = panaroo_capture)
res <- as.integer(if (is.null(res)) 0L else res)
writeLines(
  c(
    sprintf("Generated: %s", format(Sys.time())),
    sprintf("Input manifest SHA-256: %s", current_hash),
    sprintf("Exit status: %d", res),
    paste(c(shQuote(panaroo_bin), vapply(args, shQuote, character(1))), collapse = " "),
    "",
    readLines(panaroo_capture, warn = FALSE)
  ),
  command_log_file
)
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

if (!files_complete(required_outputs)) {
  missing_outputs <- required_outputs[
    !file.exists(required_outputs) | is.na(file.size(required_outputs)) | file.size(required_outputs) <= 0
  ]
  write_panaroo_report(
    manifest,
    current_hash,
    previous_hash,
    "Status: RED - Panaroo completed but one or more required outputs are absent or empty.",
    sprintf("Missing/empty required outputs: %s", paste(missing_outputs, collapse = "; "))
  )
  stop("Panaroo completed without all required outputs: ", paste(missing_outputs, collapse = ", "))
}

output_mtimes <- file.info(required_outputs)$mtime
# Required outputs were unlinked immediately before the run; the small
# tolerance avoids false failures on filesystems with coarse mtime resolution.
refresh_cutoff <- run_started - 2
if (any(is.na(output_mtimes) | output_mtimes < refresh_cutoff)) {
  stale_outputs <- required_outputs[is.na(output_mtimes) | output_mtimes < refresh_cutoff]
  write_panaroo_report(
    manifest,
    current_hash,
    previous_hash,
    "Status: RED - Panaroo completed but one or more required outputs were not refreshed.",
    sprintf("Not refreshed: %s", paste(stale_outputs, collapse = "; "))
  )
  stop("Panaroo completed but required outputs were not refreshed: ", paste(stale_outputs, collapse = ", "))
}

writeLines(current_hash, hash_file)
write_panaroo_report(
  manifest,
  current_hash,
  previous_hash,
  "Status: GREEN - Panaroo output was generated/refreshed for current inputs."
)

log_info("12c_panaroo.R complete.")

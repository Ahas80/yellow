#!/usr/bin/env Rscript

# Remove superseded generated artifacts that cannot be part of the active
# Longcycler-only release. Raw FASTAs, provider source files, and source code are
# deliberately outside the allowed roots and are never touched.

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
  library(tibble)
})

if (!file.exists("00_config.R")) {
  stop("Run this script from the rUTIs project root.", call. = FALSE)
}

args <- commandArgs(trailingOnly = TRUE)
apply_cleanup <- "--apply" %in% args
generated_roots <- c("results", "plots", "logs")
forbidden_pattern <- regex("flye", ignore_case = TRUE)
# Content cleanup must match the assembler token as a delimited value/path
# component.  A raw substring scan also matches ordinary amino-acid sequence
# motifs (for example FLYE) inside valid Panaroo intermediates.
forbidden_content_pattern <- "(^|[^[:alnum:]])flye([^[:alnum:]]|$)"

# These files are provenance ledgers, not active scientific inputs.  They may
# intentionally name retired artifacts while documenting the repository-wide
# audit.  Keep the exemption file-specific: current final-figure manifests,
# captions, validation checks and reference-aware variant tables remain inside
# the normal forbidden-content gate.
audit_provenance_dir <- file.path("results", "figure_audit")
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

path_exists <- function(path) {
  link_target <- Sys.readlink(path)
  file.exists(path) || dir.exists(path) || (!is.na(link_target) && nzchar(link_target))
}

all_paths <- unlist(lapply(generated_roots, function(root) {
  if (!dir.exists(root)) return(character())
  list.files(root, recursive = TRUE, full.names = TRUE, include.dirs = TRUE,
             all.files = TRUE, no.. = TRUE)
}), use.names = FALSE)

path_targets <- all_paths[str_detect(all_paths, forbidden_pattern)]
relative_target <- !startsWith(path_targets, "/")
path_targets[relative_target] <- file.path(
  normalizePath(".", winslash = "/", mustWork = TRUE),
  path_targets[relative_target]
)
path_targets <- normalizePath(path_targets, winslash = "/", mustWork = FALSE)

# A neutral-named table or log can still retain a superseded input path. Remove
# such text artifacts as well so no stale provenance survives into the release.
text_globs <- c(
  "*.csv", "*.tsv", "*.txt", "*.md", "*.json", "*.jsonl", "*.ndjson",
  "*.log", "*.yaml", "*.yml", "*.html", "*.xml", "*.out", "*.err"
)
# Generated roots are intentionally gitignored, so a normal ripgrep scan can
# miss exactly the stale artifacts this release sweep is meant to remove.
rg_args <- c("-l", "-i", "--hidden", "--no-ignore", "--no-messages")
# system2() passes arguments through a shell on this platform. Quote wildcard
# globs explicitly so the project root cannot expand them into positional file
# arguments before ripgrep receives them.
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
if (is.na(rg_bin) || !nzchar(rg_bin)) stop("ripgrep is required for generated-content cleanup.")
content_targets <- suppressWarnings(system2(rg_bin, rg_args, stdout = TRUE, stderr = FALSE))
rg_status <- attr(content_targets, "status")
if (identical(rg_status, 1L)) {
  content_targets <- character()
} else if (!is.null(rg_status) && !identical(rg_status, 0L)) {
  stop("Generated-content ripgrep scan failed with exit status ", rg_status, ".")
}
content_targets <- content_targets[nzchar(content_targets)]
content_targets <- normalizePath(content_targets, winslash = "/", mustWork = FALSE)
exempted_audit_provenance <- intersect(content_targets, audit_provenance_paths)
content_targets <- setdiff(content_targets, audit_provenance_paths)

# Remove superseded, inactive result trees that this complete run does not
# regenerate.  The former bad-size rescue tree is an older mixed-assembler,
# pre-operational-status sensitivity analysis (and also contains binary/image
# artifacts that cannot be found by the text-content scan).  The release now
# has one canonical transition table under results/longitudinal and the active
# clinical intermediate lives under results/clinical/intermediate.
retired_generated_paths <- c(
  file.path("results", "sensitivity"),
  file.path("results", "intermediate")
)
retired_generated_paths <- retired_generated_paths[
  vapply(retired_generated_paths, path_exists, logical(1))
]
retired_generated_paths <- normalizePath(retired_generated_paths, winslash = "/", mustWork = FALSE)

targets <- c(path_targets, content_targets, retired_generated_paths)
targets <- unique(normalizePath(targets, winslash = "/", mustWork = FALSE))

# Never silently delete a Panaroo artifact created for the currently selected
# cohort.  A hit here means either the content matcher regressed or the live
# Panaroo output contains forbidden provenance and the pipeline must stop for
# review.  Stale pre-run Panaroo files remain cleanable because they predate the
# selected manifest.
selected_manifest <- file.path("results", "qc", "analysis_assembly_manifest.csv")
current_panaroo_names <- c(
  "combined_DNA_CDS.fasta",
  "combined_protein_CDS.fasta",
  "combined_protein_cdhit_out.txt",
  "combined_protein_cdhit_out.txt.clstr",
  "gene_data.csv",
  "gene_presence_absence.csv",
  "gene_presence_absence_roary.csv",
  "gene_presence_absence.Rtab",
  "pan_genome_reference.fa",
  "final_graph.gml",
  "summary_statistics.txt",
  "panaroo_command.log"
)
current_panaroo_paths <- normalizePath(
  file.path("results", "wgs", "pan", current_panaroo_names),
  winslash = "/",
  mustWork = FALSE
)
if (file.exists(selected_manifest)) {
  manifest_mtime <- file.info(selected_manifest)$mtime
  current_panaroo_paths <- current_panaroo_paths[
    file.exists(current_panaroo_paths) &
      !is.na(file.info(current_panaroo_paths)$mtime) &
      file.info(current_panaroo_paths)$mtime >= manifest_mtime
  ]
  protected_hits <- intersect(targets, current_panaroo_paths)
  if (length(protected_hits)) {
    stop(
      "Cleanup matched current manifest-bound Panaroo artifact(s): ",
      paste(protected_hits, collapse = "; "),
      ". Refusing to delete live/current pangenome output.",
      call. = FALSE
    )
  }
}

# Remove the deepest paths first so parent-directory removal is deterministic.
targets <- targets[order(nchar(targets), decreasing = TRUE)]
target_files <- targets[file.exists(targets) & !dir.exists(targets)]
bytes <- if (length(target_files)) sum(file.size(target_files), na.rm = TRUE) else 0

summary <- tibble(
  generated_root = generated_roots,
  matching_paths = vapply(generated_roots, function(root) {
    sum(startsWith(path_targets, normalizePath(root, winslash = "/", mustWork = FALSE)))
  }, integer(1)),
  matching_text_files = vapply(generated_roots, function(root) {
    sum(startsWith(content_targets, normalizePath(root, winslash = "/", mustWork = FALSE)))
  }, integer(1)),
  exempted_audit_provenance_files = vapply(generated_roots, function(root) {
    sum(startsWith(exempted_audit_provenance, normalizePath(root, winslash = "/", mustWork = FALSE)))
  }, integer(1)),
  retired_duplicate_paths = vapply(generated_roots, function(root) {
    sum(startsWith(retired_generated_paths, normalizePath(root, winslash = "/", mustWork = FALSE)))
  }, integer(1)),
  apply_requested = apply_cleanup
)

if (!apply_cleanup) {
  print(summary)
  if (length(exempted_audit_provenance)) {
    message("Preserving audit-only provenance file(s): ",
            paste(exempted_audit_provenance, collapse = "; "))
  }
  message(sprintf("Dry run: %d generated path(s), %.2f GiB. Re-run with --apply.",
                  length(targets), bytes / 1024^3))
  quit(save = "no", status = 0)
}

if (length(targets)) {
  failed <- targets[!vapply(targets, function(path) {
    if (!path_exists(path)) return(TRUE)
    unlink(path, recursive = TRUE, force = TRUE)
    !path_exists(path)
  }, logical(1))]
  if (length(failed)) {
    stop("Could not remove ", length(failed), " superseded generated path(s).",
         call. = FALSE)
  }
}

dir.create(file.path("results", "pipeline"), recursive = TRUE, showWarnings = FALSE)
write_csv(
  summary %>% mutate(removed_bytes_total = bytes),
  file.path("results", "pipeline", "release_cleanup_summary.csv")
)

# A prior completion marker must never survive into a new run.
markers <- c(
  file.path("results", "pipeline", "RUN_COMPLETE.txt"),
  file.path("results", "research_questions", "RUN_COMPLETE.txt"),
  file.path("results", "pipeline", "longcycler_release_claim_registry.json"),
  file.path("results", "pipeline", "DELIVERABLES_COMPLETE.txt")
)
active_markers <- markers[vapply(markers, path_exists, logical(1))]
marker_status <- vapply(active_markers, unlink, integer(1), force = TRUE)
if (any(marker_status != 0L) || any(vapply(active_markers, path_exists, logical(1)))) {
  stop("Could not invalidate one or more prior release markers.", call. = FALSE)
}

message(sprintf("Release cleanup complete: %d generated path(s), %.2f GiB.",
                length(targets), bytes / 1024^3))

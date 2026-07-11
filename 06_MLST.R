#!/usr/bin/env Rscript
# ==============================================================================
# 06_MLST.R
# ==============================================================================
#
# GOAL:
#   Build the active MLST layer for the project using RIVM/provider SeqSphere
#   MLST as the primary chromosomal ST source.
#
# WHY THIS SCRIPT EXISTS:
#   Downstream analyses should have one canonical MLST entrypoint. The old local
#   PubMLST/mlst runner has been moved to scripts/run_local_mlst_deprecated.R and
#   is retained only for provenance, sensitivity checks, and labelled fallback.
#
# ACTIVE OUTPUTS:
#   - results/mlst/mlst_provider_preferred.csv       active canonical ST table
#   - results/mlst/mlst_provider_preferred_all.csv   active isolate/assembly ST table
#   - results/mlst/mlst_provider_source_audit.csv    provider/local provenance audit
#
# LOCAL PROVENANCE OUTPUTS:
#   - results/mlst/mlst_all.tsv
#   - results/mlst/mlst_with_meta.csv
#
# Biological/Statistical purpose:
#   - Uses provider/RIVM MLST at PercGoodTargets >= 95 as the authoritative ST.
#   - Uses local MLST only as an explicit labelled fallback when provider QC95 ST
#     is missing.
#   - Preserves ST_source/ST_provider/ST_local so analyses never silently mix
#     local and provider nomenclature.
# ==============================================================================

source("00_config.R")

suppressPackageStartupMessages({
  library(readr)
})

run_pipeline_script <- function(path, label) {
  if (!file.exists(path)) {
    stop("Missing required MLST helper script: ", path)
  }
  msg("%s", label)
  sys.source(path, envir = new.env(parent = globalenv()))
}

local_runner <- file.path("scripts", "run_local_mlst_deprecated.R")
source_audit <- file.path("scripts", "compare_mlst_sources.R")
source_integrator <- file.path("R", "provider_mlst_integration.R")
source_guardrail <- file.path("scripts", "verify_mlst_source_usage.R")

msg("Starting active RIVM/provider MLST pipeline.")
msg("Active MLST target: %s", FILE_MLST_CANONICAL)

if (!file.exists(FILE_MLST_LOCAL_CANONICAL)) {
  msg(
    "Local MLST provenance is missing (%s); running deprecated local mlst runner once for fallback evidence.",
    FILE_MLST_LOCAL_CANONICAL
  )
  run_pipeline_script(local_runner, "Running deprecated local MLST provenance step.")
} else {
  msg("Using existing local MLST provenance: %s", FILE_MLST_LOCAL_CANONICAL)
}

if (!file.exists(FILE_MLST_LOCAL_CANONICAL)) {
  stop("Local MLST provenance is still missing after attempted generation: ", FILE_MLST_LOCAL_CANONICAL)
}

run_pipeline_script(source_audit, "Refreshing provider-vs-local MLST source audit.")
run_pipeline_script(source_integrator, "Building active RIVM/provider-preferred MLST outputs.")
run_pipeline_script(source_guardrail, "Verifying active scripts use provider-preferred MLST.")

if (!file.exists(FILE_MLST_PROVIDER_PREFERRED)) {
  stop("Provider-preferred MLST output was not created: ", FILE_MLST_PROVIDER_PREFERRED)
}

active_mlst <- read_csv(FILE_MLST_PROVIDER_PREFERRED, show_col_types = FALSE, progress = FALSE)
source_counts <- if ("ST_source" %in% names(active_mlst)) {
  table(active_mlst$ST_source, useNA = "ifany")
} else {
  stop("Provider-preferred MLST output lacks ST_source: ", FILE_MLST_PROVIDER_PREFERRED)
}

msg("Active RIVM/provider MLST complete.")
msg("ST source counts: %s", paste(names(source_counts), as.integer(source_counts), sep = "=", collapse = "; "))
msg("Use %s for downstream analyses.", FILE_MLST_PROVIDER_PREFERRED)

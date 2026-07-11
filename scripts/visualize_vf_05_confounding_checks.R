#!/usr/bin/env Rscript
# ==============================================================================
# Deprecated compatibility wrapper: VF visualisation module 05 (confounding)
# ==============================================================================

message("DEPRECATED: scripts/visualize_vf_05_confounding_checks.R is a compatibility wrapper.")
message("Canonical VF lineage/ST, batch, event-type, and denominator diagnostics now live in 25_vf_lineage_vf_interaction.R.")
message("The canonical script uses results/vf/vf_analysis_ready.csv; this wrapper does not read stale results/vf/vf_analysis_table.csv.")

source("25_vf_lineage_vf_interaction.R", chdir = TRUE)

#!/usr/bin/env Rscript
# ==============================================================================
# Deprecated compatibility wrapper: VF visualisation module 01 (burden)
# ==============================================================================

message("DEPRECATED: scripts/visualize_vf_01_burden.R is a compatibility wrapper.")
message("Canonical VF burden plots now live in 23_vf_cross_sectional.R and use results/vf/vf_analysis_ready.csv.")
message("This wrapper does not read stale results/vf/vf_analysis_table.csv.")

source("23_vf_cross_sectional.R", chdir = TRUE)

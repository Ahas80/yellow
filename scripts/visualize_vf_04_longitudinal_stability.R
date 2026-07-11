#!/usr/bin/env Rscript
# ==============================================================================
# Deprecated compatibility wrapper: VF visualisation module 04 (longitudinal)
# ==============================================================================

message("DEPRECATED: scripts/visualize_vf_04_longitudinal_stability.R is a compatibility wrapper.")
message("Canonical within-resident VF stability plots now live in 24_vf_longitudinal_dynamics.R and use results/vf/vf_analysis_ready.csv.")
message("This wrapper does not read stale results/vf/vf_analysis_table.csv.")

source("24_vf_longitudinal_dynamics.R", chdir = TRUE)

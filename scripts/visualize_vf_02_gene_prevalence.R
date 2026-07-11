#!/usr/bin/env Rscript
# ==============================================================================
# Deprecated compatibility wrapper: VF visualisation module 02 (gene prevalence)
# ==============================================================================

message("DEPRECATED: scripts/visualize_vf_02_gene_prevalence.R is a compatibility wrapper.")
message("Canonical VF gene prevalence plots now live in 23_vf_cross_sectional.R and use results/vf/vf_analysis_ready.csv.")
message("This wrapper does not read stale results/vf/vf_analysis_table.csv.")

source("23_vf_cross_sectional.R", chdir = TRUE)

#!/usr/bin/env Rscript
# ==============================================================================
# Deprecated compatibility wrapper: VF visualisation module 03 (category profiles)
# ==============================================================================

message("DEPRECATED: scripts/visualize_vf_03_category_profiles.R is a compatibility wrapper.")
message("Canonical VF category/module plots now live in 23_vf_cross_sectional.R and 26_vf_define_gene_modules.R.")
message("Both canonical scripts use results/vf/vf_analysis_ready.csv; this wrapper does not read stale results/vf/vf_analysis_table.csv.")

source("26_vf_define_gene_modules.R", chdir = TRUE)

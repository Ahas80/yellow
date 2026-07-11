#!/usr/bin/env Rscript
# ==============================================================================
# 12d_plasmid_analysis.R
# ------------------------------------------------------------------------------
# GOAL:
#   Analyze plasmid content and structure.
#   1. Load PlasmidFinder results.
#   2. (Optional) Cluster plasmids using Mash if plasmid-only FASTAs exist.
#   3. Generate summary of plasmid distribution.
#
# METHOD:
#   - Load 'results/mlst/plasmid_replicons_long.csv'.
#   - Aggregate by participant/timepoint.
#   - Identify shared plasmids.
#
# INPUTS:
#   - results/mlst/plasmid_replicons_long.csv
#   - results/wgs/qc_summary.csv
#
# OUTPUTS:
#   - results/wgs/plasmids/plasmid_summary.csv
# ==============================================================================

# 1. Load Configuration & Libraries
source("00_config.R")
source("R/wgs_helpers.R")

suppressPackageStartupMessages({
    library(dplyr)
    library(readr)
    library(tidyr)
})

# 2. Setup
# ------------------------------------------------------------------------------
DIR_PLASMIDS <- file.path(DIR_RESULTS, "wgs", "plasmids")
ensure_dir(DIR_PLASMIDS)

log_info("Starting 12d_plasmid_analysis.R")

FILE_PF <- file.path(DIR_RESULTS, "mlst", "plasmid_replicons_long.csv")

if (!file.exists(FILE_PF)) {
    log_warn("PlasmidFinder results not found (", FILE_PF, "). Skipping plasmid analysis.")
    quit(save = "no", status = 0)
}

# 3. Load Data
# ------------------------------------------------------------------------------
pf_df <- read_csv(FILE_PF, show_col_types = FALSE)
qc_df <- load_qc_summary()

# Filter for valid genomes
valid_ids <- qc_df %>%
    filter(QC_PASS) %>%
    pull(file_name) %>%
    tools::file_path_sans_ext()
pf_df <- pf_df %>% filter(isolate_id %in% valid_ids)

log_info("Loaded plasmid hits for ", n_distinct(pf_df$isolate_id), " valid genomes.")

# 4. Analysis
# ------------------------------------------------------------------------------
# Summary by Replicon
replicon_stats <- pf_df %>%
    count(replicon, sort = TRUE) %>%
    mutate(pct = n / n_distinct(pf_df$isolate_id) * 100)

write_csv(replicon_stats, file.path(DIR_PLASMIDS, "replicon_stats.csv"))

# Summary by Isolate (Wide)
pf_wide <- pf_df %>%
    distinct(isolate_id, replicon) %>%
    mutate(val = 1) %>%
    pivot_wider(names_from = replicon, values_from = val, values_fill = 0)

write_csv(pf_wide, file.path(DIR_PLASMIDS, "plasmid_matrix.csv"))

# 5. Advanced: Mash Clustering (Placeholder)
# ------------------------------------------------------------------------------
# If we had extracted plasmid FASTAs, we would run 'mash dist' here.
# For now, we rely on replicon profiles.

log_info("12d_plasmid_analysis.R complete.")

#!/usr/bin/env Rscript
# ==============================================================================
# 12c_panaroo.R
# ------------------------------------------------------------------------------
# Role: [Pangenome] - Perform Pangenome Analysis using Panaroo.
#
# Inputs:
#   - results/wgs/qc_summary.csv
#   - results/prokka/*.gff
#
# Outputs:
#   - results/wgs/pan/gene_presence_absence.csv
#   - results/wgs/pan/struct_presence_absence.Rtab
#
# Usage:
#   Rscript 12c_panaroo.R
#
# Biological/Statistical purpose:
#   - Defines the pangenome (core + accessory genes) of the cohort.
#   - Generates gene presence/absence matrices for downstream association testing.
# ==============================================================================

# 1. Load Configuration & Libraries
source("00_config.R")
source("R/wgs_helpers.R")

suppressPackageStartupMessages({
    library(dplyr)
    library(readr)
    library(fs)
    library(stringr)
})

# 2. Setup
# ------------------------------------------------------------------------------
DIR_PAN <- file.path(DIR_WGS, "pan")
ensure_dir(DIR_PAN)

# Assuming GFFs are in results/prokka/ (standard pipeline)
# Adjust if they are elsewhere
DIR_GFFS <- file.path(DIR_RESULTS, "prokka")

log_info("Starting 12c_panaroo.R")

# Check tools
has_panaroo <- check_wgs_tool("panaroo")
if (!has_panaroo) stop("Panaroo is required for this module.")

# 3. Load QC Data & Map GFFs
# ------------------------------------------------------------------------------
qc_df <- load_qc_summary()
valid_ids <- qc_df %>%
    filter(QC_PASS) %>%
    pull(file_name) %>%
    tools::file_path_sans_ext()

log_info("Looking for GFFs for ", length(valid_ids), " valid genomes...")

# Find all GFFs
all_gffs <- dir_ls(DIR_GFFS, recurse = TRUE, glob = "*.gff")

# Match GFFs to valid IDs
# Assuming GFF filename contains the ID or is in a folder named ID
# Heuristic: Check if ID is in the path
valid_gffs <- c()
for (id in valid_ids) {
    # Strict match: filename starts with ID or folder is ID
    matches <- all_gffs[str_detect(path_file(all_gffs), fixed(id))]
    if (length(matches) > 0) {
        valid_gffs <- c(valid_gffs, matches[1]) # Take first match
    } else {
        log_warn("Missing GFF for ID: ", id)
    }
}

if (length(valid_gffs) < 2) {
    stop("Not enough GFFs found for pangenome analysis.")
}

log_info("Found ", length(valid_gffs), " matching GFF files.")

# 4. Run Panaroo
# ------------------------------------------------------------------------------
# --clean-mode strict: remove contamination
# --remove-invalid-genes: robust parsing
# -a core: core genome alignment (optional, we did Parsnp)
# --threads

cmd_panaroo <- paste(
    "panaroo",
    "-i", paste(valid_gffs, collapse = " "),
    "-o", DIR_PAN,
    "--clean-mode", "strict",
    "--remove-invalid-genes",
    "-t", CORES_USE,
    "--verbose"
)

# Check if already run (Panaroo fails if output dir exists and is not empty)
if (file.exists(file.path(DIR_PAN, "gene_presence_absence.csv"))) {
    log_warn("Panaroo output exists. Skipping run. Delete ", DIR_PAN, " to force rerun.")
} else {
    log_info("Running Panaroo (this may take a while)...")
    # Panaroo can fail if output dir exists, so we might need to clean it carefully
    # But we want to keep the dir structure. Panaroo usually complains if *files* exist.

    res <- system(cmd_panaroo)

    if (res != 0) {
        log_error("Panaroo failed with exit code ", res)
        stop("Panaroo execution failed.")
    }

    log_info("Panaroo completed successfully.")
}

log_info("12c_panaroo.R complete.")

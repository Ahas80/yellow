#!/usr/bin/env Rscript
# ==============================================================================
# repair_prokka.R
# ------------------------------------------------------------------------------
# Purpose:
#   Identify samples where Prokka failed (missing .gff files) and re-run Prokka
#   with the '--compliant' flag to handle long contig IDs.
#
# Usage:
#   Rscript repair_prokka.R
# ==============================================================================

source("00_config.R")
source("R/wgs_helpers.R")
suppressPackageStartupMessages({
    library(dplyr)
    library(fs)
    library(stringr)
    library(parallel)
})

# 1. Setup
# ------------------------------------------------------------------------------
log_info("Starting Prokka repair...")

if (!dir.exists(DIR_FASTAS)) {
    stop("FASTA directory not found: ", DIR_FASTAS)
}

ensure_dir(DIR_PROKKA)

# Set PATH to include conda env bin (for java, tbl2asn, etc.)
conda_bin <- "/Users/Aamir/miniforge_x86/envs/yellow-wgs-x86/bin"
java_home <- "/Users/Aamir/miniforge_x86/envs/yellow-wgs-x86/lib/jvm"

Sys.setenv(PATH = paste(conda_bin, Sys.getenv("PATH"), sep = ":"))
Sys.setenv(JAVA_HOME = java_home)


# 2. Identify Missing GFFs
# ------------------------------------------------------------------------------
# Get all input FASTA files
fasta_files <- dir_ls(DIR_FASTAS, glob = "*.fasta")
sample_ids <- path_ext_remove(path_file(fasta_files))

log_info("Found ", length(sample_ids), " input FASTA files.")

# Check for corresponding GFFs
missing_samples <- c()

for (id in sample_ids) {
    # Expected output directory and GFF file
    out_dir <- file.path(DIR_PROKKA, id)
    gff_file <- file.path(out_dir, paste0(id, ".gff"))

    if (!file.exists(gff_file)) {
        missing_samples <- c(missing_samples, id)
    }
}

if (length(missing_samples) == 0) {
    log_info("All samples have valid GFF files. No repair needed.")
    quit(save = "no")
}

log_info("Found ", length(missing_samples), " samples with missing GFF files.")
log_info("First 5 missing: ", paste(head(missing_samples, 5), collapse = ", "))

# 3. Define Repair Function
# ------------------------------------------------------------------------------
run_prokka_compliant <- function(sample_id, fasta_path, out_dir) {
    # Ensure output directory is clean (Prokka fails if dir exists)
    if (dir.exists(out_dir)) {
        fs::dir_delete(out_dir)
    }

    # Construct Prokka command
    # --compliant: Force Genbank/Sequin compliant names (fixes long ID issue)
    # --centre X: Required by --compliant
    # --force: Overwrite output
    cmd <- paste(
        "/Users/Aamir/miniforge_x86/envs/yellow-wgs-x86/bin/perl",
        "/Users/Aamir/miniforge_x86/envs/yellow-wgs-x86/bin/prokka",
        "--outdir", out_dir,
        "--prefix", sample_id,
        "--compliant",
        "--centre", "X",
        "--cpus", "1", # Use 1 cpu per job, parallelize via R
        fasta_path
    )

    res <- system(cmd, ignore.stdout = TRUE, ignore.stderr = TRUE)

    if (res == 0) {
        return(TRUE)
    } else {
        return(FALSE)
    }
}

# 4. Execute Repair
# ------------------------------------------------------------------------------
log_info("Re-running Prokka for missing samples (using ", CORES_USE, " cores)...")

# Prepare arguments for parallel execution
jobs <- data.frame(
    id = missing_samples,
    fasta = file.path(DIR_FASTAS, paste0(missing_samples, ".fasta")),
    out = file.path(DIR_PROKKA, missing_samples),
    stringsAsFactors = FALSE
)

# Run in parallel
results <- mclapply(1:nrow(jobs), function(i) {
    id <- jobs$id[i]
    fasta <- jobs$fasta[i]
    out <- jobs$out[i]

    msg("Processing %s (%d/%d)...", id, i, nrow(jobs))
    success <- run_prokka_compliant(id, fasta, out)

    if (success) {
        return(id)
    } else {
        return(NULL)
    }
}, mc.cores = CORES_USE)

# 5. Summary
# ------------------------------------------------------------------------------
success_ids <- unlist(results)
failed_count <- length(missing_samples) - length(success_ids)

log_info("Repair complete.")
log_info("Successfully repaired: ", length(success_ids))
log_info("Failed: ", failed_count)

if (failed_count > 0) {
    log_warn("Some samples failed to process. Check logs.")
} else {
    log_info("All missing GFFs have been generated.")
}

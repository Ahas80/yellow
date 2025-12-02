#!/usr/bin/env Rscript
# ==============================================================================
# 12b_core_snp.R
# ------------------------------------------------------------------------------
# Role: [Phylogeny] - Perform core genome SNP calling and alignment.
#
# Inputs:
#   - results/wgs/qc_summary.csv
#   - data/assemblies/*.fasta
#
# Outputs:
#   - results/wgs/core/core.aln.fasta
#   - results/wgs/core/snp_dists.tsv
#   - results/wgs/core/strain_pairs.csv
#   - results/wgs/core/core_genome.tree
#
# Usage:
#   Rscript 12b_core_snp.R
#
# Biological/Statistical purpose:
#   - Constructs a core genome alignment to infer phylogenetic relationships.
#   - Calculates pairwise SNP distances to identify "Same" vs "Different" strains.
# ==============================================================================

# 1. Load Configuration & Libraries
source("00_config.R")
source("R/wgs_helpers.R")

suppressPackageStartupMessages({
    library(dplyr)
    library(readr)
    library(fs)
    library(ape)
    library(tidyr)
})

# 2. Setup
# ------------------------------------------------------------------------------
DIR_CORE <- file.path(DIR_WGS, "core")
ensure_dir(DIR_CORE)

log_info("Starting 12b_core_snp.R")

# Check if outputs already exist
tree_file <- file.path(DIR_CORE, "core_genome.tree")
dist_file <- file.path(DIR_CORE, "snp_dists.tsv")

if (file.exists(tree_file) && file.exists(dist_file)) {
    log_info("Core genome tree and SNP distances already exist at %s", DIR_CORE)
    log_info("Skipping Parsnp run to save time. Delete 'results/wgs/core' to force re-run.")
    quit(save = "no", status = 0)
}

# Check tools
has_parsnp <- check_wgs_tool("parsnp")
has_snp_dists <- check_wgs_tool("snp-dists")

if (!has_parsnp) stop("Parsnp is required for this module.")

# 3. Load QC Data
# ------------------------------------------------------------------------------
qc_df <- load_qc_summary()
valid_genomes <- qc_df %>%
    filter(QC_PASS) %>%
    pull(full_path)

if (length(valid_genomes) < 2) {
    stop("Not enough valid genomes for core SNP analysis (Need >= 2, found ", length(valid_genomes), ")")
}

log_info("Processing ", length(valid_genomes), " QC-passed genomes.")

# 4. Prepare Input for Parsnp
# ------------------------------------------------------------------------------
# Parsnp takes a directory of FASTAs or a file list.
# We'll create a temporary directory with symlinks to ensure clean input.
temp_fasta_dir <- file.path(DIR_CORE, "temp_fastas")
if (dir_exists(temp_fasta_dir)) dir_delete(temp_fasta_dir)
dir_create(temp_fasta_dir)

log_info("Staging FASTAs in ", temp_fasta_dir)
# Symlink files
for (f in valid_genomes) {
    link_name <- file.path(temp_fasta_dir, fs::path_file(f))
    if (!file.exists(link_name)) fs::link_create(f, link_name)
}

# 5. Run Parsnp
# ------------------------------------------------------------------------------
# -c: core genome alignment
# -r !: random reference
# -p: threads
# -o: output dir
# -n: skip tree construction (avoids RAxML "Too few species" error on identical seqs)
parsnp_out <- file.path(DIR_CORE, "parsnp_out")
if (dir_exists(parsnp_out)) dir_delete(parsnp_out) # Parsnp requires clean dir

cmd_parsnp <- paste(
    "parsnp",
    "-c",
    "-r !",
    "-d", temp_fasta_dir,
    "-o", parsnp_out,
    "-p", CORES_USE,
    "-x", # Filter duplicates (helps with "Too few species")
    "--verbose"
)

log_info("Running Parsnp...")
res <- system(cmd_parsnp)

# Check if alignment exists even if Parsnp failed (e.g. at tree step)
xmfa_file <- file.path(parsnp_out, "parsnp.xmfa")

if (res != 0) {
    if (file.exists(xmfa_file)) {
        log_warn("Parsnp exited with error (likely RAxML failure), but alignment exists. Proceeding...")
    } else {
        log_error("Parsnp failed with exit code ", res)
        stop("Parsnp execution failed and no alignment found.")
    }
} else {
    log_info("Parsnp completed successfully.")
}

# 5b. Convert XMFA/GGR to FASTA (if needed)
# ------------------------------------------------------------------------------
aln_file <- file.path(parsnp_out, "parsnp.fasta")
ggr_file <- file.path(parsnp_out, "parsnp.ggr")

if (!file.exists(aln_file) && file.exists(ggr_file)) {
    log_info("Converting GGR to FASTA...")
    # harvesttools -i <ggr> -M <fasta>
    cmd_conv <- paste("harvesttools -i", ggr_file, "-M", aln_file)
    system(cmd_conv)
}

# 6. Run snp-dists & Build Tree
# ------------------------------------------------------------------------------
dist_file <- file.path(DIR_CORE, "snp_dists.tsv")
tree_file <- file.path(DIR_CORE, "core_genome.tree")

if (file.exists(aln_file) && has_snp_dists) {
    log_info("Running snp-dists...")
    cmd_dists <- paste("snp-dists", aln_file, ">", dist_file)
    system(cmd_dists)
    log_info("SNP distances written to ", dist_file)

    # 6b. Build Tree (Robust Fix for "Too Few Species")
    # ------------------------------------------------------------------------------
    # Parsnp's RAxML often fails on identical sequences. We build a robust NJ tree here.
    log_info("Building Neighbor-Joining tree from SNP distances...")

    # Load distances as matrix
    dists_df <- read_tsv(dist_file, show_col_types = FALSE)
    dist_mat <- as.matrix(dists_df[, -1])
    rownames(dist_mat) <- colnames(dist_mat)

    # Build NJ tree
    tree <- ape::nj(as.dist(dist_mat))

    # Save tree
    ape::write.tree(tree, file = tree_file)
    log_info("Tree saved to ", tree_file)

    # 7. Pairwise Classification
    # ------------------------------------------------------------------------------
    log_info("Classifying strain pairs...")

    # Load distances (already loaded as dists_df, but keeping flow)
    dists <- dists_df

    # Convert to long format
    # Columns are SampleIDs, rows are SampleIDs (first col is id)
    # Rename first col to 'A'
    colnames(dists)[1] <- "A"

    pairs_long <- dists %>%
        pivot_longer(-A, names_to = "B", values_to = "snps") %>%
        filter(A < B) # Unique pairs only

    # Load Genome Sizes for normalization (optional, but good for SNPs/Mb)
    # We'll use a fixed threshold for now based on standard E. coli (5Mb)
    # Heuristic:
    #   Same Strain: <= 20 SNPs (approx 4 SNPs/Mb)
    #   Related:     <= 1000 SNPs
    #   Different:   > 1000 SNPs

    pairs_classified <- pairs_long %>%
        mutate(
            call = case_when(
                snps <= 20 ~ "Same",
                snps <= 1000 ~ "Related",
                TRUE ~ "Different"
            )
        )

    # Save
    pairs_file <- file.path(DIR_CORE, "strain_pairs.csv")
    write_csv(pairs_classified, pairs_file)
    log_info("Classified ", nrow(pairs_classified), " pairs. Saved to ", pairs_file)

    # Summary
    log_info("Pair Summary:")
    print(table(pairs_classified$call))
} else {
    log_warn("Skipping snp-dists (Alignment missing or tool not found).")
}

# Cleanup
# dir_delete(temp_fasta_dir) # Optional: keep for debugging

log_info("12b_core_snp.R complete.")

#!/usr/bin/env Rscript
# ==============================================================================
# 12b_core_snp.R
# ==============================================================================
#
# GOAL:
#   Perform core genome SNP calling and alignment to quantify genetic distance
#   between isolates.  SNP distances are used by script 11 (strain comparison)
#   and script 15 (longitudinal patterns) to distinguish "same strain" from
#   "different strain" within the same participant.
#
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

tree_file <- file.path(DIR_CORE, "core_genome.tree")
dist_file <- file.path(DIR_CORE, "snp_dists.tsv")
manifest_file <- file.path(DIR_CORE, "core_snp_input_manifest.csv")
hash_file <- file.path(DIR_CORE, "core_snp_input_manifest.hash")
stale_report <- file.path(DIR_CORE, "core_snp_staleness_report.txt")
force_core <- identical(Sys.getenv("FORCE_RERUN_CORE_SNP", "0"), "1")

# 3. Load current selected assembly inputs and fingerprint them before deciding
# whether an existing core-SNP directory is current.
canonical_file <- file.path(DIR_QC, "canonical_assembly_selection.csv")
if (file.exists(canonical_file)) {
    qc_df <- read_csv(canonical_file, show_col_types = FALSE)
    selected_df <- qc_df %>%
        filter(selected_canonical %in% TRUE, QC_PASS %in% TRUE)
    input_source <- canonical_file
} else {
    qc_df <- load_qc_summary()
    selected_df <- qc_df %>% filter(QC_PASS %in% TRUE)
    input_source <- file.path(DIR_WGS, "qc_summary.csv")
    log_warn("Canonical assembly selection not found; using all QC PASS assemblies for fingerprinting.")
}

if (!"tp_lab" %in% names(selected_df) && "Timepoint" %in% names(selected_df)) {
    selected_df$tp_lab <- normalise_timepoint_preserve_events(selected_df$Timepoint)
}
if (!"Assembly_ID" %in% names(selected_df)) {
    selected_df$Assembly_ID <- tools::file_path_sans_ext(basename(selected_df$full_path))
}
if (!"Assembler" %in% names(selected_df)) {
    selected_df$Assembler <- if ("assembler" %in% names(selected_df)) selected_df$assembler else detect_assembler(selected_df$full_path)
}

manifest <- selected_df %>%
    transmute(
        Assembly_ID,
        Participant_id = as.character(Participant_id),
        tp_lab = normalise_timepoint_preserve_events(tp_lab),
        Assembler = Assembler,
        full_path = normalizePath(full_path, winslash = "/", mustWork = FALSE),
        file_size = ifelse(file.exists(full_path), file.size(full_path), NA_real_),
        modified_time = ifelse(file.exists(full_path), as.character(file.info(full_path)$mtime), NA_character_)
    ) %>%
    arrange(Participant_id, tp_lab, Assembly_ID)

write_csv(manifest, manifest_file)
current_hash <- hash_input_manifest(manifest)
previous_hash <- if (file.exists(hash_file)) readLines(hash_file, warn = FALSE)[1] else NA_character_
outputs_exist <- file.exists(tree_file) && file.exists(dist_file)

writeLines(
    c(
        "Core SNP input staleness report",
        sprintf("Generated: %s", format(Sys.time())),
        sprintf("Input source: %s", input_source),
        sprintf("Selected assemblies in current manifest: %d", nrow(manifest)),
        sprintf("Current manifest hash: %s", current_hash),
        sprintf("Previous manifest hash: %s", ifelse(is.na(previous_hash), "<none>", previous_hash)),
        sprintf("Outputs exist: %s", outputs_exist),
        sprintf("FORCE_RERUN_CORE_SNP: %s", Sys.getenv("FORCE_RERUN_CORE_SNP", "0")),
        if (outputs_exist && identical(current_hash, previous_hash)) "Status: GREEN - outputs match current input manifest."
        else if (outputs_exist) "Status: RED - Core SNP outputs are stale relative to current assembly inputs."
        else "Status: MISSING - core SNP outputs are absent."
    ),
    stale_report
)

append_denominator_summary(
    manifest,
    "12b_core_snp.R",
    "core_snp_input_manifest",
    "participant_timepoint",
    manifest_file,
    "Fingerprint of canonical selected assemblies used to decide whether Parsnp outputs are current"
)

if (outputs_exist && identical(current_hash, previous_hash)) {
    log_info("Core genome tree and SNP distances already exist and input manifest hash matches.")
    log_info("Skipping Parsnp run.")
    quit(save = "no", status = 0)
}

if (outputs_exist && !force_core) {
    log_warn("Core SNP outputs are stale relative to current assembly inputs.")
    log_warn("Not rerunning Parsnp because FORCE_RERUN_CORE_SNP is not 1.")
    log_warn("Review ", stale_report, " and rerun with FORCE_RERUN_CORE_SNP=1 when ready.")
    quit(save = "no", status = 0)
}

# Check tools
has_parsnp <- check_wgs_tool("parsnp")
has_snp_dists <- check_wgs_tool("snp-dists")

if (!has_parsnp) stop("Parsnp is required for this module.")

valid_genomes <- manifest %>%
    filter(!is.na(full_path), file.exists(full_path)) %>%
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
    #   Same Strain: <= configured YELLOW study SNP threshold
    #   Related:     <= 1000 SNPs
    #   Different:   > 1000 SNPs
    same_strain_snp_threshold <- strain_snp_threshold()

    pairs_classified <- pairs_long %>%
        mutate(
            call = case_when(
                snps <= same_strain_snp_threshold ~ "Same",
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

if (file.exists(tree_file) && file.exists(dist_file)) {
    writeLines(current_hash, hash_file)
    writeLines(
        c(
            "Core SNP input staleness report",
            sprintf("Generated: %s", format(Sys.time())),
            sprintf("Input source: %s", input_source),
            sprintf("Selected assemblies in current manifest: %d", nrow(manifest)),
            sprintf("Current manifest hash: %s", current_hash),
            sprintf("Previous manifest hash: %s", ifelse(is.na(previous_hash), "<none>", previous_hash)),
            "Status: GREEN - outputs were generated/refreshed for the current input manifest."
        ),
        stale_report
    )
}

# Cleanup
# dir_delete(temp_fasta_dir) # Optional: keep for debugging

log_info("12b_core_snp.R complete.")

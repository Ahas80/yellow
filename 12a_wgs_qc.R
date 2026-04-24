#!/usr/bin/env Rscript
# ==============================================================================
# 12a_wgs_qc.R
# ==============================================================================
#
# GOAL:
#   Perform quality control on all genome assemblies: N50, contig count,
#   total bases, GC content.  Identifies poor-quality assemblies that should
#   be flagged or excluded from downstream comparisons.
#
# ------------------------------------------------------------------------------
# Role: [QC] - Perform Quality Control (QC) on all assemblies.
#
# Inputs:
#   - assembly_metadata.csv
#   - data/assemblies/*.fasta
#
# Outputs:
#   - results/wgs/qc_summary.csv
#   - plots/wgs/wgs_qc_n50_vs_contigs.png
#
# Usage:
#   Rscript 12a_wgs_qc.R
#
# Biological/Statistical purpose:
#   - Filters out low-quality assemblies (contamination, fragmentation) to ensure
#     downstream analyses (SNP, Pangenome) are robust.
# ==============================================================================

# 1. Load Configuration & Libraries
source("00_config.R")
source("R/wgs_helpers.R")

suppressPackageStartupMessages({
    library(dplyr)
    library(readr)
    library(tidyr)
    library(purrr)
    library(furrr)
    library(seqinr) # For FASTA reading
})

# 2. Setup
# ------------------------------------------------------------------------------
QC_CONFIG <- get_qc_config()
DIR_WGS_OUT <- DIR_WGS
ensure_dir(DIR_WGS_OUT)

log_info("Starting 12a_wgs_qc.R")
log_info("QC Thresholds:")
log_info("  Max Contigs: ", QC_CONFIG$MAX_CONTIGS)
log_info("  Min N50:     ", QC_CONFIG$MIN_N50)
log_info("  Size Range:  ", QC_CONFIG$MIN_GENOME_SIZE, " - ", QC_CONFIG$MAX_GENOME_SIZE)

# 3. Load Metadata
# ------------------------------------------------------------------------------
if (!file.exists(FILE_METADATA)) stop("Missing ", FILE_METADATA)
meta_df <- read_csv(FILE_METADATA, show_col_types = FALSE)

# Ensure full_path
if (!"full_path" %in% names(meta_df)) {
    meta_df$full_path <- file.path(DIR_FASTAS, meta_df$file_name)
}

# Filter missing files
meta_df <- meta_df %>% filter(file.exists(full_path))
log_info("Found ", nrow(meta_df), " assemblies to process.")

# 4. Compute Metrics
# ------------------------------------------------------------------------------
compute_assembly_stats <- function(fasta_path) {
    tryCatch(
        {
            seqs <- seqinr::read.fasta(fasta_path, seqtype = "DNA", as.string = TRUE, forceDNAtolower = FALSE)

            # Calculate lengths
            lens <- vapply(seqs, nchar, integer(1))
            total_len <- sum(lens)
            n_contigs <- length(lens)

            # N50
            lens_sorted <- sort(lens, decreasing = TRUE)
            csum <- cumsum(lens_sorted)
            n50_idx <- which(csum >= total_len / 2)[1]
            n50 <- lens_sorted[n50_idx]

            # GC Content
            # Concatenate all sequences for global GC
            all_seq <- paste(unlist(seqs), collapse = "")
            gc_count <- stringr::str_count(all_seq, "[GCgc]")
            gc_pct <- (gc_count / total_len) * 100

            tibble(
                total_bp = total_len,
                n_contigs = n_contigs,
                N50 = n50,
                GC_pct = gc_pct,
                error = NA_character_
            )
        },
        error = function(e) {
            tibble(
                total_bp = NA_real_,
                n_contigs = NA_integer_,
                N50 = NA_real_,
                GC_pct = NA_real_,
                error = as.character(e$message)
            )
        }
    )
}

# Run in parallel
plan(multisession, workers = CORES_USE)
log_info("Computing metrics using ", CORES_USE, " cores...")

stats_res <- meta_df %>%
    mutate(stats = future_map(full_path, compute_assembly_stats, .progress = TRUE)) %>%
    unnest(stats)

plan(sequential)

# 5. Apply Filters
# ------------------------------------------------------------------------------
qc_res <- stats_res %>%
    mutate(
        QC_PASS = TRUE,
        QC_REASON = ""
    ) %>%
    rowwise() %>%
    mutate(
        # Check Contigs
        QC_PASS = ifelse(n_contigs > QC_CONFIG$MAX_CONTIGS, FALSE, QC_PASS),
        QC_REASON = ifelse(n_contigs > QC_CONFIG$MAX_CONTIGS, paste0(QC_REASON, "HighContigs;"), QC_REASON),

        # Check N50
        QC_PASS = ifelse(N50 < QC_CONFIG$MIN_N50, FALSE, QC_PASS),
        QC_REASON = ifelse(N50 < QC_CONFIG$MIN_N50, paste0(QC_REASON, "LowN50;"), QC_REASON),

        # Check Size
        QC_PASS = ifelse(total_bp < QC_CONFIG$MIN_GENOME_SIZE | total_bp > QC_CONFIG$MAX_GENOME_SIZE, FALSE, QC_PASS),
        QC_REASON = ifelse(total_bp < QC_CONFIG$MIN_GENOME_SIZE | total_bp > QC_CONFIG$MAX_GENOME_SIZE, paste0(QC_REASON, "BadSize;"), QC_REASON),

        # Check Errors
        QC_PASS = ifelse(!is.na(error), FALSE, QC_PASS),
        QC_REASON = ifelse(!is.na(error), paste0(QC_REASON, "ReadError;"), QC_REASON)
    ) %>%
    ungroup() %>%
    mutate(QC_REASON = ifelse(QC_PASS, "PASS", QC_REASON))

# 6. Summary & Output
# ------------------------------------------------------------------------------
pass_count <- sum(qc_res$QC_PASS)
fail_count <- sum(!qc_res$QC_PASS)

log_info("QC Complete.")
log_info("  PASS: ", pass_count)
log_info("  FAIL: ", fail_count)

outfile <- file.path(DIR_WGS_OUT, "qc_summary.csv")
write_csv(qc_res, outfile)
log_info("Written QC summary to: ", outfile)

# Optional: Plot QC metrics
if (pass_count > 0) {
    # N50 vs Contigs Plot
    library(ggplot2)
    g <- ggplot(qc_res, aes(x = n_contigs, y = N50, color = QC_PASS)) +
        geom_point(alpha = 0.6) +
        scale_y_log10() +
        scale_x_log10() +
        geom_vline(xintercept = QC_CONFIG$MAX_CONTIGS, linetype = "dashed") +
        geom_hline(yintercept = QC_CONFIG$MIN_N50, linetype = "dashed") +
        labs(
            title = "Assembly QC: N50 vs Contig Count",
            subtitle = paste("Pass:", pass_count, "| Fail:", fail_count)
        ) +
        theme_minimal()

    ggsave(file.path(DIR_PLOTS_WGS, "wgs_qc_n50_vs_contigs.png"), g, width = 7, height = 5)
}

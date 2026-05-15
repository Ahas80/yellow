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
meta_df <- meta_df %>% filter(!is.na(full_path) & file.exists(full_path))
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

# Canonical assembly selection for participant-timepoint analyses.  QC remains
# assembly-level, but downstream biological episode analyses should use only the
# selected row per Participant_id x tp_lab.
canonical_selection <- select_canonical_assemblies(qc_res)
canonical_file <- file.path(DIR_QC, "canonical_assembly_selection.csv")
write_csv(canonical_selection, canonical_file)
log_info("Written canonical assembly selection to: ", canonical_file)

append_denominator_summary(
    qc_res,
    "12a_wgs_qc.R",
    "assembly_qc",
    "assembly",
    FILE_METADATA,
    "Assembly-level QC; assembler alternatives remain explicit rows"
)
append_denominator_summary(
    canonical_selection %>% filter(selected_canonical),
    "12a_wgs_qc.R",
    "canonical_selected_assemblies",
    "participant_timepoint",
    canonical_file,
    "One QC PASS selected assembly per Participant_id x tp_lab; longcycler preferred over flye when both pass"
)

# QC selection bias by infection status.  This uses exact-style Fisher testing
# because UTI counts are sparse and the chi-square approximation can be invalid.
status_file <- FILE_STATUS_MAP
if (file.exists(status_file)) {
    status <- read_csv(status_file, show_col_types = FALSE) %>%
        mutate(
            Participant_id = as.character(Participant_id),
            tp_lab = if ("tp_lab" %in% names(.)) normalise_timepoint_preserve_events(tp_lab) else normalise_timepoint_preserve_events(Timepoint)
        )
    status_dupes <- assert_unique_keys(
        status,
        c("Participant_id", "tp_lab"),
        context = "status_map for QC selection bias",
        out_path = file.path(DIR_QC, "status_map_duplicate_episode_keys.csv")
    )
    if (nrow(status_dupes) == 0) {
        qc_episode <- canonical_selection %>%
            group_by(Participant_id, tp_lab) %>%
            summarise(
                has_wgs = any(file_exists, na.rm = TRUE),
                qc_pass = any(QC_PASS, na.rm = TRUE),
                selected_canonical = any(selected_canonical, na.rm = TRUE),
                .groups = "drop"
            )

        bias_df <- status %>%
            left_join(qc_episode, by = c("Participant_id", "tp_lab"), relationship = "one-to-one") %>%
            mutate(
                has_wgs = coalesce(has_wgs, FALSE),
                qc_pass = coalesce(qc_pass, FALSE),
                selected_canonical = coalesce(selected_canonical, FALSE)
            )

        bias_summary <- bias_df %>%
            dplyr::count(Infection_Status, qc_pass, name = "n") %>%
            group_by(Infection_Status) %>%
            mutate(status_total = sum(n), pct = n / status_total) %>%
            ungroup()
        write_csv(bias_summary, file.path(DIR_QC, "qc_selection_bias_by_status.csv"))

        tbl <- xtabs(n ~ qc_pass + Infection_Status, data = bias_summary)
        fisher_res <- tryCatch(fisher.test(tbl), error = function(e) NULL)
        worst_loss <- bias_summary %>%
            filter(qc_pass %in% FALSE) %>%
            arrange(desc(pct), desc(n)) %>%
            slice_head(n = 1)
        report <- c(
            "QC selection bias by infection status",
            sprintf("Generated: %s", format(Sys.time())),
            "",
            capture.output(print(tbl)),
            "",
            if (!is.null(fisher_res)) sprintf("Fisher exact p-value: %.5g", fisher_res$p.value) else "Fisher exact test could not be computed.",
            if (nrow(worst_loss) > 0) sprintf("Largest QC loss proportion: %s (%.1f%%, n=%d)",
                                             worst_loss$Infection_Status[1], 100 * worst_loss$pct[1], worst_loss$n[1]) else "No QC losses by status.",
            "Interpretation: if QC pass differs by status, genomic comparisons may be selection-biased and should be described as conditional on available QC-passing WGS."
        )
        writeLines(report, file.path(DIR_QC, "qc_selection_bias_report.txt"))
    } else {
        writeLines(
            c(
                "QC selection bias by infection status",
                sprintf("Generated: %s", format(Sys.time())),
                "",
                "RED: status_map has duplicated Participant_id + tp_lab keys. Bias analysis was not collapsed silently.",
                "Review results/qc/status_map_duplicate_episode_keys.csv."
            ),
            file.path(DIR_QC, "qc_selection_bias_report.txt")
        )
    }
}

write_uti_attrition_outputs()

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

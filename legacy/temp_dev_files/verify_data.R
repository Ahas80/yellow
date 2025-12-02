#!/usr/bin/env Rscript
source("00_config.R")
library(tidyverse)

cat("=== Phase 4: Data Verification ===\n\n")

# 1. Check qc_summary.csv Participant_id correctness
cat("1. Checking QC Summary Participant IDs...\n")
qc_file <- "results/wgs/wgs_qc_summary.csv"
if (file.exists(qc_file)) {
    qc <- read_csv(qc_file, show_col_types = FALSE)
    cat("   - Total rows:", nrow(qc), "\n")
    cat("   - Participant_id column present:", "Participant_id" %in% names(qc), "\n")
    if ("Participant_id" %in% names(qc)) {
        # Check for PR#### batch IDs (should NOT be present)
        batch_ids <- qc %>% filter(str_detect(Participant_id, "^PR[0-9]+$"))
        cat("   - Rows with batch IDs (PR####):", nrow(batch_ids), "(should be 0)\n")
        # Check for correct IDs (5-6 digit numbers)
        correct_ids <- qc %>% filter(str_detect(Participant_id, "^[0-9]{5,6}$"))
        cat("   - Rows with correct IDs (5-6 digits):", nrow(correct_ids), "\n")
        cat("   - Sample Participant IDs:", paste(head(qc$Participant_id, 5), collapse = ", "), "\n")
    }
} else {
    cat("   - QC summary file not found at", qc_file, "\n")
}

# 2. Check strain_pairs.csv within-host logic
cat("\n2. Checking strain_pairs.csv within-host logic...\n")
pairs_file <- "results/wgs/core/strain_pairs.csv"
if (file.exists(pairs_file)) {
    pairs <- read_csv(pairs_file, show_col_types = FALSE)
    cat("   - Total pairs:", nrow(pairs), "\n")
    cat("   - Column names:", paste(names(pairs), collapse = ", "), "\n")
    # Extract Participant IDs from Isolate IDs
    pairs <- pairs %>%
        mutate(
            Pid_A = str_extract(A, "^[A-Z0-9]+"),
            Pid_B = str_extract(B, "^[A-Z0-9]+"),
            within_host = Pid_A == Pid_B
        )
    cat("   - Within-host pairs:", sum(pairs$within_host, na.rm = TRUE), "\n")
    cat("   - Between-host pairs:", sum(!pairs$within_host, na.rm = TRUE), "\n")
    cat("   - Same/Related/Different breakdown:\n")
    print(table(pairs$call))
} else {
    cat("   - strain_pairs.csv not found\n")
}

# 3. Check pairwise_metrics.csv Participant_id grouping
cat("\n3. Checking pairwise_metrics.csv Participant_id grouping...\n")
metrics_file <- "results/strain_compare/pairwise_metrics.csv"
if (file.exists(metrics_file)) {
    metrics <- read_csv(metrics_file, show_col_types = FALSE)
    cat("   - Total rows:", nrow(metrics), "\n")
    cat("   - Participant_id_A present:", "Participant_id_A" %in% names(metrics), "\n")
    cat("   - Participant_id_B present:", "Participant_id_B" %in% names(metrics), "\n")
    if ("Participant_id_A" %in% names(metrics)) {
        # Check for batch IDs
        batch_A <- metrics %>% filter(str_detect(Participant_id_A, "^PR[0-9]+$"))
        batch_B <- metrics %>% filter(str_detect(Participant_id_B, "^PR[0-9]+$"))
        cat("   - Rows with batch ID in A:", nrow(batch_A), "(should be 0)\n")
        cat("   - Rows with batch ID in B:", nrow(batch_B), "(should be 0)\n")
        # Check within-participant pairs
        cat("   - Within-participant TRUE:", sum(metrics$within_participant, na.rm = TRUE), "\n")
        cat("   - Within-participant FALSE:", sum(!metrics$within_participant, na.rm = TRUE), "\n")
        cat("   - Sample Participant pairs:\n")
        print(head(metrics[, c("Participant_id_A", "Participant_id_B", "within_participant")], 5))
    }
} else {
    cat("   - pairwise_metrics.csv not found\n")
}

# 4. Check GWAS model convergence
cat("\n4. Checking GWAS model convergence...\n")
gwas_file <- "results/models/gwas_multivariable_glmm.csv"
if (file.exists(gwas_file)) {
    gwas <- read_csv(gwas_file, show_col_types = FALSE)
    cat("   - Total genes tested:", nrow(gwas), "\n")
    cat("   - Column names:", paste(names(gwas), collapse = ", "), "\n")
    cat("   - Significant genes (p<0.05):", sum(gwas$p_value < 0.05, na.rm = TRUE), "\n")
    if ("p_value" %in% names(gwas)) {
        sig <- gwas %>%
            filter(p_value < 0.05) %>%
            arrange(p_value)
        cat("   - Top significant genes:\n")
        print(sig[, c("gene", "p_value")], n = 10)
    }
} else {
    cat("   - GWAS results file not found\n")
}

cat("\n=== Verification Complete ===\n")

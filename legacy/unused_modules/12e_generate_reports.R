#!/usr/bin/env Rscript
# ==============================================================================
# 12e_generate_reports.R
# ------------------------------------------------------------------------------
# Role: [Reporting] - Generate per-participant WGS reports.
#
# Inputs:
#   - results/wgs/qc_summary.csv
#   - results/wgs/core/strain_pairs.csv
#   - results/wgs/pan/gene_presence_absence.csv
#
# Outputs:
#   - results/reports/Participant_X_WGS_Report.pdf
#
# Usage:
#   Rscript 12e_generate_reports.R
#
# Biological/Statistical purpose:
#   - Summarizes WGS findings for each participant to facilitate clinical interpretation.
#   - Integrates QC, Phylogeny, and Pangenome results.
# ==============================================================================

# 1. Load Configuration & Libraries
source("00_config.R")
source("R/wgs_helpers.R")

suppressPackageStartupMessages({
    library(dplyr)
    library(readr)
    library(rmarkdown)
    library(ggplot2)
})

# 2. Setup
# ------------------------------------------------------------------------------
DIR_REPORTS <- DIR_REPORTS
ensure_dir(DIR_REPORTS)

log_info("Starting 12e_generate_reports.R")

# 3. Load Data
# ------------------------------------------------------------------------------
qc_df <- load_qc_summary()
valid_pids <- qc_df %>%
    filter(QC_PASS) %>%
    pull(Participant_id) %>%
    unique()

log_info("Generating reports for ", length(valid_pids), " participants.")

# 4. Define Report Template (Inline for simplicity, or external file)
# ------------------------------------------------------------------------------
# We'll create a simple Rmd file on the fly
template_path <- file.path("R", "wgs_report_template.Rmd")

if (!file.exists(template_path)) {
    log_info("Creating report template at ", template_path)
    cat('---
title: "WGS Report: Participant {{pid}}"
date: "`r Sys.Date()`"
output: pdf_document
params:
  pid: "NA"
  qc_df: "NA"
  pairs_path: "NA"
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = FALSE, warning = FALSE, message = FALSE)
library(dplyr)
library(readr)
library(stringr)
library(knitr)
```

## 1. Assembly Quality Control

```{r qc_table}
p_qc <- params$qc_df %>% filter(Participant_id == params$pid) %>% select(Timepoint, n_contigs, N50, total_bp, QC_PASS)
kable(p_qc, caption = "Assembly Metrics")
```

## 2. Core Genome SNPs

```{r snp_pairs}
if(file.exists(params$pairs_path)) {
  pairs <- read_csv(params$pairs_path, show_col_types = FALSE)

  # Filter for pairs involving this participant (assuming SampleID contains PID)
  # Robust check: A or B contains the PID string
  my_pairs <- pairs %>%
    filter(str_detect(A, params$pid) | str_detect(B, params$pid)) %>%
    arrange(snps)

  if(nrow(my_pairs) > 0) {
    kable(head(my_pairs, 10), caption = "Top Related Strains (SNP distance)")
  } else {
    cat("No close relatives found involving this participant.")
  }
} else {
  cat("Strain comparison data not available (run 12b_core_snp.R).")
}
```

## 3. Pangenome & Plasmids

(Placeholder for gene content summary)

', file = template_path)
}

# 5. Generate Reports
# ------------------------------------------------------------------------------
pairs_csv <- file.path(DIR_WGS, "core", "strain_pairs.csv")

for (pid in valid_pids) {
    outfile <- file.path(DIR_REPORTS, paste0("Participant_", pid, "_WGS_Report.pdf"))

    # Skip if exists (optional)
    # if (file.exists(outfile)) next

    tryCatch(
        {
            render(
                input = template_path,
                output_file = outfile,
                params = list(
                    pid = pid,
                    qc_df = qc_df,
                    pairs_path = pairs_csv
                ),
                quiet = TRUE
            )
            log_info("Generated report for ", pid)
        },
        error = function(e) {
            log_error("Failed to generate report for ", pid, ": ", e$message)
        }
    )
}

log_info("12e_generate_reports.R complete.")

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
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = FALSE, warning = FALSE, message = FALSE)
library(dplyr)
library(ggplot2)
library(knitr)
```

## 1. Assembly Quality Control

```{r qc_table}
p_qc <- params$qc_df %>% filter(Participant_id == params$pid) %>% select(Timepoint, n_contigs, N50, total_bp, QC_PASS)
kable(p_qc, caption = "Assembly Metrics")
```

## 2. Core Genome SNPs

(Placeholder for SNP tree or matrix)
```{r snp_pairs}
pairs_file <- file.path("results/wgs/core/strain_pairs.csv")
if(file.exists(pairs_file)) {
  pairs <- read_csv(pairs_file, show_col_types = FALSE)
  my_pairs <- pairs %>% filter(A == params$pid | B == params$pid)
  if(nrow(my_pairs) > 0) {
    kable(head(my_pairs, 10), caption = "Top Related Strains (SNP distance)")
  } else {
    cat("No close relatives found.")
  }
}
```

## 3. Pangenome & Plasmids

(Placeholder for gene content summary)

', file = template_path)
}

# 5. Generate Reports
# ------------------------------------------------------------------------------
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
                    qc_df = qc_df
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

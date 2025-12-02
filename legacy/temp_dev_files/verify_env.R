#!/usr/bin/env Rscript
# verify_env.R
# Comprehensive check of all R dependencies for the rUTIs pipeline

cat("========================================\n")
cat("Verifying R Environment for rUTIs Pipeline\n")
cat("========================================\n\n")

# List of all packages used across all scripts
required_packages <- c(
    # Core
    "here", "optparse", "tidyverse", "data.table", "readr", "dplyr", "tidyr",
    "stringr", "ggplot2", "tibble", "forcats", "purrr", "vroom", "fs", "glue",

    # Statistics & Modeling
    "lme4", "broom.mixed", "future", "future.apply", "furrr",

    # Visualization
    "pheatmap", "ComplexUpset", "gridExtra", "scales", "ggrepel", "ggraph",
    "patchwork", "ragg", "RColorBrewer", "knitr", "rmarkdown",

    # Genomics/Bioinformatics (Bioconductor)
    "Biostrings", "ggtree", "ape", "ComplexHeatmap", "seqinr",

    # System/Utils
    "jsonlite", "lubridate", "processx", "testthat"
)

missing <- c()
loaded <- c()

for (pkg in required_packages) {
    if (requireNamespace(pkg, quietly = TRUE)) {
        cat(sprintf("✓ %-20s : Installed\n", pkg))
        loaded <- c(loaded, pkg)
    } else {
        cat(sprintf("❌ %-20s : MISSING\n", pkg))
        missing <- c(missing, pkg)
    }
}

cat("\n========================================\n")
if (length(missing) > 0) {
    cat("⚠️  MISSING PACKAGES FOUND:\n")
    cat(paste(missing, collapse = ", "), "\n\n")
    cat("Please run the following to fix:\n")

    # Generate conda install command for missing packages
    # Map some R packages to their conda names if different
    conda_pkgs <- paste0("r-", tolower(missing))
    # Fix known naming differences
    conda_pkgs <- gsub("r-biostrings", "bioconductor-biostrings", conda_pkgs)
    conda_pkgs <- gsub("r-ggtree", "bioconductor-ggtree", conda_pkgs)
    conda_pkgs <- gsub("r-complexheatmap", "bioconductor-complexheatmap", conda_pkgs)
    conda_pkgs <- gsub("r-ape", "r-ape", conda_pkgs) # usually r-ape
    conda_pkgs <- gsub("r-complexupset", "r-complexupset", conda_pkgs)

    cat("conda install -c conda-forge -c bioconda", paste(conda_pkgs, collapse = " "), "\n")
    quit(status = 1)
} else {
    cat("✅ All dependencies are installed and loadable.\n")
    cat("   You are ready to run the pipeline!\n")
    quit(status = 0)
}

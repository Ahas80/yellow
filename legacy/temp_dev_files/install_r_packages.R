#!/usr/bin/env Rscript
# Install required R packages for rUTIs pipeline

cat("Installing required R packages...\n")
cat("This may take a few minutes...\n\n")

# CRAN packages
cran_packages <- c(
    "here", # Path management
    "optparse", # CLI arguments
    "tidyverse", # Data manipulation (includes dplyr, ggplot2, tidyr, readr, purrr, tibble, stringr, forcats)
    "data.table", # Fast data handling
    "vroom", # Fast file reading
    "fs", # File system operations
    "glue", # String interpolation
    "jsonlite", # JSON handling
    "lubridate", # Date handling
    "pheatmap", # Heatmaps
    "ComplexUpset", # UpSet plots
    "RColorBrewer", # Color palettes
    "gridExtra", # Multiple plots
    "patchwork", # Plot composition
    "ggrepel", # Label repulsion
    "ggraph", # Network visualization
    "ragg", # High-quality graphics device
    "scales", # Plot scales
    "igraph", # Network analysis
    "future", # Parallel processing core
    "future.apply", # Parallel apply
    "furrr", # Parallel mapping
    "lme4", # Mixed models
    "broom.mixed", # Model tidying
    "seqinr", # Sequence analysis
    "knitr", # Dynamic report generation
    "rmarkdown", # Markdown rendering
    "processx", # External process control
    "testthat", # Unit testing
    "randomForest", # Random Forest models
    "broom", # Model tidying
    "BiocManager" # For Bioconductor packages
)

# Bioconductor packages
bioc_packages <- c(
    "Biostrings", # Sequence manipulation
    "ComplexHeatmap", # Advanced heatmaps
    "ggtree", # Tree visualization
    "ape" # Phylogenetics
)

cat("=== STEP 1/2: Installing CRAN packages ===\n")

# Check which packages are already installed
installed <- installed.packages()[, "Package"]
to_install <- cran_packages[!cran_packages %in% installed]

if (length(to_install) > 0) {
    cat("Installing", length(to_install), "CRAN packages:\n")
    cat(paste(to_install, collapse = ", "), "\n\n")

    install.packages(to_install,
        repos = "https://cloud.r-project.org",
        dependencies = TRUE,
        Ncpus = 4
    )
    cat("✓ CRAN packages installed!\n\n")
} else {
    cat("✓ All CRAN packages already installed!\n\n")
}

cat("=== STEP 2/2: Installing Bioconductor packages ===\n")

# Ensure BiocManager is installed
if (!"BiocManager" %in% installed.packages()[, "Package"]) {
    install.packages("BiocManager", repos = "https://cloud.r-project.org")
}

# Check which Bioconductor packages need installation
installed <- installed.packages()[, "Package"]
bioc_to_install <- bioc_packages[!bioc_packages %in% installed]

if (length(bioc_to_install) > 0) {
    cat("Installing", length(bioc_to_install), "Bioconductor packages:\n")
    cat(paste(bioc_to_install, collapse = ", "), "\n\n")

    BiocManager::install(bioc_to_install, update = FALSE, ask = FALSE)
    cat("✓ Bioconductor packages installed!\n\n")
} else {
    cat("✓ All Bioconductor packages already installed!\n\n")
}

cat("\n✅ ALL R PACKAGES INSTALLATION COMPLETE!\n")
cat("You can now run: bash RUN_COMPLETE_ANALYSIS.sh\n")

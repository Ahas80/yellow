#!/usr/bin/env Rscript
# ==============================================================================
# 10_replicon_heatmap.R
# ------------------------------------------------------------------------------
# Role: [Descriptive] - Generate a presence/absence heatmap of plasmid replicons.
#
# Inputs:
#   - results/plasmids/plasmidfinder_presence_absence.csv
#   - results/mlst/mlst_all.tsv
#
# Outputs:
#   - plots/plasmids/replicon_heatmap.pdf
#   - plots/plasmids/replicon_heatmap.png
#
# Usage:
#   Rscript 10_replicon_heatmap.R
#
# Biological/Statistical purpose:
#   - Visualizes the plasmid content landscape across the cohort.
#   - Highlights associations between specific plasmids and STs.
# ==============================================================================

# 1. Load Configuration & Libraries
source("00_config.R")
suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
  library(RColorBrewer)
  library(optparse)
})

# 2. CLI Options
# ------------------------------------------------------------------------------
option_list <- list(
  make_option(c("-i", "--input"),
    type = "character", default = file.path(DIR_PLASMIDS, "plasmidfinder_presence_absence.csv"),
    help = "Input plasmid presence/absence CSV [default: %default]"
  ),
  make_option(c("-m", "--mlst"),
    type = "character", default = file.path(DIR_MLST, "mlst_all.tsv"),
    help = "Input MLST TSV file [default: %default]"
  ),
  make_option(c("-o", "--outdir"),
    type = "character", default = DIR_PLOTS_PLASMIDS,
    help = "Output directory for plots [default: %default]"
  )
)

opt <- parse_args(OptionParser(option_list = option_list))

# 3. Setup & Validation
# ------------------------------------------------------------------------------
ensure_dir(opt$outdir)
ensure_dir(file.path(DIR_LOGS, "debug"))

FILE_MAT <- opt$input
FILE_MLST <- opt$mlst
OUT_PDF <- file.path(opt$outdir, "replicon_heatmap.pdf")
OUT_PNG <- file.path(opt$outdir, "replicon_heatmap.png")

# Check for heatmap packages
use_complex <- requireNamespace("ComplexHeatmap", quietly = TRUE)
use_pheat <- requireNamespace("pheatmap", quietly = TRUE)

if (!use_complex && !use_pheat) {
  stop("Neither ComplexHeatmap nor pheatmap is installed. Please install one.")
}

msg("Starting 10_replicon_heatmap.R")
msg("Input Matrix: %s", FILE_MAT)
msg("Input MLST:   %s", FILE_MLST)
msg("Output Dir:   %s", opt$outdir)

# 4. Load Data
# ------------------------------------------------------------------------------
if (!file.exists(FILE_MAT)) stop("Missing input matrix: ", FILE_MAT)
mat_df <- read_csv(FILE_MAT, show_col_types = FALSE) %>% as.data.frame()

if (nrow(mat_df) == 0) stop("Empty matrix file.")

# Handle Isolate_ID
if (!"Isolate_ID" %in% names(mat_df)) {
  # Try first column if Isolate_ID missing
  rownames(mat_df) <- mat_df[[1]]
  mat_df <- mat_df[, -1]
} else {
  rownames(mat_df) <- mat_df$Isolate_ID
  mat_df$Isolate_ID <- NULL
}

mat <- as.matrix(mat_df)

if (anyNA(mat)) {
  warning("NAs in matrix, coercing to 0.")
  mat[is.na(mat)] <- 0
}

msg("Matrix dimensions: %d isolates x %d replicons", nrow(mat), ncol(mat))

# 5. Annotations (ST)
# ------------------------------------------------------------------------------
st_vec <- rep("Unknown", nrow(mat))
names(st_vec) <- rownames(mat)

if (file.exists(FILE_MLST)) {
  mlst <- read_tsv(FILE_MLST, show_col_types = FALSE)
  if ("ST" %in% names(mlst)) {
    # Ensure Isolate_ID matching
    match_idx <- match(rownames(mat), mlst$Isolate_ID)
    st_vals <- mlst$ST[match_idx]
    st_vec[!is.na(st_vals)] <- st_vals[!is.na(st_vals)]
  } else {
    warning("MLST file found but no 'ST' column.")
  }
} else {
  warning("MLST file not found: ", FILE_MLST)
}

# Group rare STs
st_counts <- sort(table(st_vec), decreasing = TRUE)
keep_st <- names(st_counts)[seq_len(min(20, length(st_counts)))]
st_vec_plot <- ifelse(st_vec %in% keep_st, st_vec, "Other_ST")

# Colors
palette_fun <- colorRampPalette(brewer.pal(8, "Set3"))
st_levels <- sort(unique(st_vec_plot))
st_colours <- setNames(palette_fun(length(st_levels)), st_levels)

# 6. Filter & Order
# ------------------------------------------------------------------------------
rep_prev <- colSums(mat, na.rm = TRUE)
keep_cols <- which(rep_prev >= 5)

if (length(keep_cols) == 0) {
  msg("No replicons in >=5 isolates, keeping all.")
  heat <- mat
} else {
  heat <- mat[, keep_cols, drop = FALSE]
}

col_order <- order(colSums(heat, na.rm = TRUE), decreasing = TRUE)
row_order <- order(rowSums(heat, na.rm = TRUE), decreasing = TRUE)

# Provenance
writeLines(colnames(heat)[col_order], file.path(DIR_LOGS, "debug", "replicons_kept_order.txt"))
writeLines(rownames(heat)[row_order], file.path(DIR_LOGS, "debug", "isolates_order.txt"))

# 7. Plot
# ------------------------------------------------------------------------------
msg("Generating plots...")

if (use_complex) {
  library(ComplexHeatmap)
  ha_row <- rowAnnotation(ST = st_vec_plot, col = list(ST = st_colours), show_annotation_name = FALSE)

  pdf(OUT_PDF, width = 10, height = 8)
  draw(Heatmap(heat,
    name = "Inc\npresent", col = c("white", "black"),
    show_row_names = FALSE, cluster_rows = TRUE, cluster_columns = TRUE,
    row_order = row_order, column_order = col_order,
    column_title = "Plasmid Replicon Presence/Absence Heatmap",
    left_annotation = ha_row
  ))
  dev.off()
} else {
  library(pheatmap)
  ann <- data.frame(ST = st_vec_plot, row.names = rownames(heat))

  pdf(OUT_PDF, width = 10, height = 8)
  pheatmap(heat,
    color = c("white", "black"), cluster_rows = TRUE, cluster_cols = TRUE,
    annotation_row = ann, annotation_colors = list(ST = st_colours),
    show_rownames = FALSE, border_color = NA,
    main = "Plasmid Replicon Presence/Absence Heatmap"
  )
  dev.off()
}

# PNG version
if (file.exists(OUT_PDF) && requireNamespace("ragg", quietly = TRUE)) {
  ragg::agg_png(OUT_PNG, width = 1600, height = 1200, res = 150)
  op <- par(mar = c(4, 4, 2, 1))
  image(t(heat[row_order, col_order, drop = FALSE]), axes = FALSE, col = c("white", "black"))
  par(op)
  dev.off()
}

msg("✓ Plasmid heatmap generated at %s", OUT_PDF)

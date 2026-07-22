#!/usr/bin/env Rscript
# ==============================================================================
# 10_replicon_heatmap.R
# ==============================================================================
#
# GOAL:
#   Generate a clustered presence/absence heatmap of plasmid replicons across
#   all isolates, annotated by ST.  This visualises which plasmid types are
#   shared across lineages and which are lineage-specific.
#
# The input columns are exact PlasmidFinder GENE labels. ACCESSION values are
# metadata only and must never be used as feature keys.
#
# ------------------------------------------------------------------------------
# Role: [Descriptive] - Generate a presence/absence heatmap of plasmid replicons.
#
# Inputs:
#   - results/plasmids/plasmidfinder_presence_absence.csv
#   - results/mlst/mlst_provider_preferred.csv
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
    type = "character", default = FILE_MLST_CANONICAL,
    help = "Input active provider-preferred MLST file [default: %default]"
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
selection <- load_analysis_assemblies(FILE_ANALYSIS_ASSEMBLY_MANIFEST, require_files = TRUE)
expected_isolates <- as.character(selection$Isolate_ID)

if (!"Isolate_ID" %in% names(mat_df)) {
  stop("Plasmid presence/absence matrix lacks Isolate_ID; positional fallback is not allowed.")
}
matrix_isolates <- as.character(mat_df$Isolate_ID)
if (nrow(mat_df) != 532L || anyDuplicated(matrix_isolates) ||
    !setequal(matrix_isolates, expected_isolates)) {
  stop("Plasmid heatmap input does not exactly match the 532 selected Longcycler assemblies.")
}
rownames(mat_df) <- matrix_isolates
mat_df$Isolate_ID <- NULL

catalog_path <- file.path(DIR_PLASMIDS, "plasmidfinder_replicon_catalog.csv")
if (!file.exists(catalog_path)) {
  stop("Missing canonical PlasmidFinder GENE/accession catalog: ", catalog_path)
}
catalog <- read_csv(catalog_path, show_col_types = FALSE)
if (!all(c("GENE", "accession") %in% names(catalog))) {
  stop("PlasmidFinder catalog lacks GENE/accession columns.")
}
if (any(names(mat_df) %in% unique(catalog$accession)) ||
    !all(names(mat_df) %in% unique(catalog$GENE))) {
  stop("Heatmap features must be exact PlasmidFinder GENE labels, never accession IDs.")
}
ap001918_genes <- catalog %>%
  filter(accession == "AP001918") %>%
  distinct(GENE) %>%
  pull(GENE)
if (!all(c("IncFIB(AP001918)_1", "IncFIC(FII)_1", "IncFIA_1") %in% ap001918_genes)) {
  stop("AP001918 gene-level feature separation failed.")
}

if (!all(vapply(mat_df, is.numeric, logical(1)))) {
  bad <- names(mat_df)[!vapply(mat_df, is.numeric, logical(1))]
  stop(
    "Replicon matrix contains non-numeric column(s): ",
    paste(bad, collapse = ", ")
  )
}
mat <- as.matrix(mat_df)
if (ncol(mat) == 0L) stop("Replicon matrix contains no replicon columns.")
invalid_binary <- !is.na(mat) & !mat %in% c(0, 1)
if (any(invalid_binary)) {
  bad_cols <- unique(colnames(mat)[col(mat)[invalid_binary]])
  stop(
    "Replicon matrix contains non-binary non-missing values in column(s): ",
    paste(bad_cols, collapse = ", ")
  )
}

msg("Matrix dimensions: %d isolates x %d replicons", nrow(mat), ncol(mat))
msg("Unavailable replicon calls retained as unavailable: %d", sum(is.na(mat)))

# 5. Annotations (ST)
# ------------------------------------------------------------------------------
st_vec <- rep("Unknown", nrow(mat))
names(st_vec) <- rownames(mat)

if (!file.exists(FILE_MLST)) stop("Active MLST file is missing: ", FILE_MLST)
mlst <- if (grepl("\\.tsv$", FILE_MLST, ignore.case = TRUE)) {
  read_tsv(FILE_MLST, show_col_types = FALSE)
} else {
  read_csv(FILE_MLST, show_col_types = FALSE)
}
if (!all(c("Isolate_ID", "ST") %in% names(mlst))) {
  stop("Active MLST file lacks Isolate_ID or ST.")
}
mlst_isolates <- as.character(mlst$Isolate_ID)
if (nrow(mlst) != 532L || anyDuplicated(mlst_isolates) ||
    !setequal(mlst_isolates, expected_isolates)) {
  stop("Active MLST rows do not exactly match the 532 selected Longcycler assemblies.")
}
match_idx <- match(rownames(mat), mlst_isolates)
st_vals <- mlst$ST[match_idx]
st_vec[!is.na(st_vals)] <- st_vals[!is.na(st_vals)]

# Group rare STs
st_counts <- sort(table(st_vec), decreasing = TRUE)
keep_st <- names(st_counts)[seq_len(min(20, length(st_counts)))]
st_vec_plot <- ifelse(st_vec %in% keep_st, st_vec, "Other_ST")
names(st_vec_plot) <- names(st_vec)

# Colors
palette_fun <- colorRampPalette(brewer.pal(8, "Set3"))
st_levels <- unique(c(keep_st, if (any(st_vec_plot == "Other_ST")) "Other_ST"))
st_levels <- st_levels[st_levels %in% unique(st_vec_plot)]
st_colours <- setNames(palette_fun(length(st_levels)), st_levels)
if (any(!st_vec_plot %in% names(st_colours))) {
  stop("ST annotation palette does not cover every plotted ST group.")
}

# 6. Filter & Order
# ------------------------------------------------------------------------------
rep_present <- colSums(mat == 1, na.rm = TRUE)
rep_available <- colSums(!is.na(mat))
rep_missing <- nrow(mat) - rep_available
keep_cols <- which(rep_present >= 5)

if (length(keep_cols) == 0) {
  msg("No replicons in >=5 isolates, keeping all.")
  heat <- mat
} else {
  heat <- mat[, keep_cols, drop = FALSE]
}

# Deterministic biological/descriptive ordering. Rows are grouped by the ST
# annotation, then by observed replicon burden; columns are ordered by observed
# prevalence and name. Default Euclidean clustering is inappropriate for binary
# presence/absence data and would also mishandle unavailable calls.
heat_present <- colSums(heat == 1, na.rm = TRUE)
col_order <- order(-heat_present, colnames(heat))
row_burden <- rowSums(heat == 1, na.rm = TRUE)
row_missing <- rowSums(is.na(heat))
row_order <- order(
  match(st_vec_plot, st_levels),
  -row_burden,
  row_missing,
  rownames(heat)
)
heat <- heat[row_order, col_order, drop = FALSE]

# Provenance
writeLines(colnames(heat), file.path(DIR_LOGS, "debug", "replicons_kept_order.txt"))
writeLines(rownames(heat), file.path(DIR_LOGS, "debug", "isolates_order.txt"))
write_csv(
  tibble(
    replicon = colnames(mat),
    n_present = unname(rep_present),
    n_available = unname(rep_available),
    n_unavailable = unname(rep_missing),
    retained = colnames(mat) %in% colnames(heat)
  ),
  file.path(DIR_LOGS, "debug", "replicon_heatmap_call_summary.csv")
)

# 7. Plot
# ------------------------------------------------------------------------------
msg("Generating plots...")

replicon_state_colours <- c(
  "Absent" = "#F2F2F2",
  "Present" = "#1A1A1A",
  "Unavailable" = "#7A7A7A"
)
plot_title <- paste0(
  "PlasmidFinder replicon-marker presence/absence\n",
  "Rows grouped by ST; exact GENE labels ordered by prevalence; grey = unavailable"
)

render_device <- function(path, device = c("pdf", "png"), draw_fun) {
  device <- match.arg(device)
  if (identical(device, "pdf")) {
    grDevices::pdf(path, width = 10, height = 8, bg = "white")
  } else {
    grDevices::png(
      path, width = 10, height = 8, units = "in", res = 300,
      bg = "white"
    )
  }
  device_open <- TRUE
  on.exit({
    if (device_open) grDevices::dev.off()
  }, add = TRUE)
  draw_fun()
  grDevices::dev.off()
  device_open <- FALSE
  invisible(path)
}

if (use_complex) {
  heat_state <- matrix(
    "Unavailable", nrow = nrow(heat), ncol = ncol(heat),
    dimnames = dimnames(heat)
  )
  heat_state[!is.na(heat) & heat == 0] <- "Absent"
  heat_state[!is.na(heat) & heat == 1] <- "Present"
  st_annotation <- factor(st_vec_plot[rownames(heat)], levels = st_levels)
  if (anyNA(st_annotation)) stop("Ordered ST annotation contains unmapped values.")

  draw_heatmap <- function() {
    ha_row <- ComplexHeatmap::rowAnnotation(
      ST = st_annotation,
      col = list(ST = st_colours),
      show_annotation_name = FALSE
    )
    ComplexHeatmap::draw(ComplexHeatmap::Heatmap(
      heat_state,
      name = "Replicon state",
      col = replicon_state_colours,
      show_row_names = FALSE,
      cluster_rows = FALSE,
      cluster_columns = FALSE,
      column_title = plot_title,
      left_annotation = ha_row,
      heatmap_legend_param = list(
        at = names(replicon_state_colours),
        labels = names(replicon_state_colours)
      )
    ))
  }
  render_device(OUT_PDF, "pdf", draw_heatmap)
  render_device(OUT_PNG, "png", draw_heatmap)
} else {
  heat_plot <- heat
  heat_plot[is.na(heat_plot)] <- 2
  ann <- data.frame(
    ST = factor(st_vec_plot[rownames(heat)], levels = st_levels),
    row.names = rownames(heat)
  )
  if (anyNA(ann$ST)) stop("Ordered ST annotation contains unmapped values.")
  draw_heatmap <- function() {
    pheatmap::pheatmap(
      heat_plot,
      color = unname(replicon_state_colours),
      breaks = c(-0.5, 0.5, 1.5, 2.5),
      legend_breaks = c(0, 1, 2),
      legend_labels = names(replicon_state_colours),
      cluster_rows = FALSE,
      cluster_cols = FALSE,
      annotation_row = ann,
      annotation_colors = list(ST = st_colours),
      show_rownames = FALSE,
      border_color = NA,
      main = plot_title
    )
  }
  render_device(OUT_PDF, "pdf", draw_heatmap)
  render_device(OUT_PNG, "png", draw_heatmap)
}

heat_outputs <- c(OUT_PDF, OUT_PNG)
if (any(!file.exists(heat_outputs)) || any(file.info(heat_outputs)$size <= 0)) {
  stop("One or more replicon heatmap outputs are missing or empty.")
}

msg("✓ Plasmid heatmaps generated at %s and %s", OUT_PDF, OUT_PNG)

#!/usr/bin/env Rscript
# ==============================================================================
# 05_gene_overview_plots.R
# ==============================================================================
#
# GOAL:
#   Visualise VF gene prevalence across the cohort and classify genes as
#   "core" (present in >95% of isolates) or "variable" (present in <95%).
#   Produces heatmaps and bar charts for thesis/poster figures.
#
# WHY THIS SCRIPT EXISTS:
#   Understanding which VF genes are universally carried (core) vs sporadically
#   present (variable) is biologically important.  Core genes (e.g., fimH)
#   are likely essential for colonisation regardless of clinical outcome.
#   Variable genes are the candidates most likely to differ between UTI and
#   Not_UTI — they are the features worth testing in the GLMM (script 14).
#
# INPUTS:
#   - results/vf/vf_pa_all.csv           (from 02_gene_presence_analysis.R)
#
# OUTPUTS:
#   - results/vf/core_gene_list.txt      (genes present in >95% of isolates)
#   - results/vf/variable_gene_list.txt  (genes present in <95%)
#   - plots/vf/gene_prevalence_bar.png
#   - plots/vf/variable_gene_heatmap.png, .pdf
# ==============================================================================
#
# Usage:
#   Rscript 05_gene_overview_plots.R
#
# Biological/Statistical purpose:
#   - Distinguishes "core" VF genes (high prevalence) from "variable" ones.
#   - Visualizes the distribution of variable genes across the cohort.
# ==============================================================================

# 1. Load Configuration & Libraries
source("00_config.R")
suppressPackageStartupMessages({
  library(tidyverse)
  library(pheatmap)
})

# 2. Load Data
# ------------------------------------------------------------------------------
if (!file.exists(FILE_VF_PA)) stop("Missing ", FILE_VF_PA)
vf <- read_csv(FILE_VF_PA, show_col_types = FALSE)

ensure_dir(DIR_PLOTS_VF)

# 3. Timepoint Normalization
# ------------------------------------------------------------------------------
tp_norm <- function(x) {
  x <- as.character(x)
  is_uricult <- str_detect(x, regex("uricult", ignore_case = TRUE))
  tp_num <- suppressWarnings(as.integer(str_extract(x, "\\d+")))
  tp_lab <- case_when(
    is_uricult ~ "Uricult",
    !is.na(tp_num) ~ paste0("T", tp_num),
    TRUE ~ "Unscheduled"
  )
  tp_levels <- c(
    paste0("T", sort(unique(tp_num[!is.na(tp_num)]))),
    "Uricult", "Unscheduled"
  )
  tibble(tp_lab = factor(tp_lab, levels = tp_levels))
}

if (!"tp_lab" %in% names(vf)) {
  if ("Timepoint" %in% names(vf)) {
    vf <- vf |> bind_cols(tp_norm(vf$Timepoint))
  } else {
    stop("vf_pa_all.csv lacks both 'tp_lab' and 'Timepoint' columns")
  }
}

# Script 02 now keeps provenance columns such as Episode_ID in vf_pa_all.
# Only numeric binary columns are VF genes; metadata must never enter gene
# prevalence or heatmap calculations.
meta_cols <- intersect(c(
  "Participant_id", "tp_lab", "Timepoint", "Episode_ID", "Assembly_ID",
  "Assembly_Base_ID", "Isolate_ID", "Event_type", "Batch", "Collection_Date",
  "UTI_Label", "Urine_collection_method", "Infection_Status", "Status_Confidence_epi",
  "Sx_source_epi", "ST", "file_name", "full_path", "fasta_path"
), names(vf))
gene_cols <- canonical_vf_gene_cols(names(vf), vf_pa_file = FILE_VF_PA)
if (length(gene_cols) == 0) {
  stop("No numeric VF gene columns found in ", FILE_VF_PA)
}

# 4. Tidy Matrix
# ------------------------------------------------------------------------------
# Blank or missing calls mean unavailable, not biological absence. Non-missing
# gene cells must already be explicit binary presence/absence calls.
vf[gene_cols] <- vf[gene_cols] |>
  mutate(across(everything(), ~ {
    raw <- trimws(as.character(.x))
    raw_missing <- is.na(.x) | !nzchar(raw)
    value <- suppressWarnings(as.numeric(raw))
    value[raw_missing] <- NA_real_
    invalid_parse <- !raw_missing & is.na(value)
    invalid_binary <- !is.na(value) & !value %in% c(0, 1)
    if (any(invalid_parse) || any(invalid_binary)) {
      stop(
        "VF gene column '", cur_column(),
        "' contains non-binary or unparseable non-missing values."
      )
    }
    value
  }))

collapse_binary_calls <- function(x) {
  if (all(is.na(x))) return(NA_real_)
  if (any(x == 1, na.rm = TRUE)) return(1)
  # Absence is supported only when every contributing call is available and 0.
  if (anyNA(x)) return(NA_real_)
  0
}

vf <- vf |>
  mutate(
    Participant_id = as.character(Participant_id),
    tp_lab = as.character(tp_lab)
  )

# Collapse exact duplicate metadata rows conservatively: presence can be
# supported by any available call, whereas absence requires all calls to be 0.
vf <- vf |>
  group_by(across(all_of(meta_cols))) |>
  summarise(across(all_of(gene_cols), collapse_binary_calls), .groups = "drop") |>
  arrange(across(all_of(meta_cols)))

# Matrix for Heatmap
if (anyNA(vf$Participant_id) || anyNA(vf$tp_lab)) {
  stop("VF heatmap rows require non-missing Participant_id and tp_lab.")
}
row_ids <- paste(vf$Participant_id, vf$tp_lab, sep = "_")
if (anyDuplicated(row_ids)) {
  stop(
    "VF heatmap row labels are not unique at participant-timepoint level; ",
    "refusing to merge or silently relabel repeated observations."
  )
}
mat <- vf |>
  unite(row_id, Participant_id, tp_lab, sep = "_", remove = FALSE) |>
  column_to_rownames("row_id") |>
  select(all_of(gene_cols)) |>
  as.matrix()

# 5. Core vs Variable Genes
# ------------------------------------------------------------------------------
n_iso <- nrow(mat)
n_genes <- ncol(mat)
message(sprintf("Matrix dimensions: %d isolates x %d genes", n_iso, n_genes))

# [BIO] NOTE: "Core" here means high-prevalence VF genes (≥95% of isolates)
# This is NOT the same as the bacterial core genome (housekeeping genes)
# VF genes from VFDB are mostly ACCESSORY (pathotype-specific)
# True E. coli core genome: ~2000-3000 genes; VFDB captures dozens of VFs

# Calculate prevalence per gene
gene_present <- colSums(mat == 1, na.rm = TRUE)
gene_absent <- colSums(mat == 0, na.rm = TRUE)
gene_available <- colSums(!is.na(mat))
gene_missing <- n_iso - gene_available
if (any(gene_available == 0)) {
  stop(
    "VF matrix contains gene(s) with no available calls: ",
    paste(names(gene_available)[gene_available == 0], collapse = ", ")
  )
}
gene_prevalence <- gene_present / gene_available
message("Unavailable VF gene calls retained as unavailable: ", sum(gene_missing))

# [BIO] Use 95% threshold instead of 100% (accounts for assembly/annotation gaps)
core_threshold <- 0.95
core_genes <- names(which(gene_prevalence >= core_threshold))

# Variable here means at least one observed presence and one observed absence.
# This preserves the historical presence/absence-variation definition; it is
# not the complement of the >=95% high-prevalence list.
variable_genes <- names(which(gene_present > 0 & gene_absent > 0))

write_lines(core_genes, file.path(DIR_VF, "core_gene_list.txt"))
write_lines(variable_genes, file.path(DIR_VF, "variable_gene_list.txt"))

message("Core genes (>=95%): ", length(core_genes))
message("Variable genes (observed present and absent): ", length(variable_genes))
message("Absent genes among available calls (0%): ", length(which(gene_present == 0)))

if (length(variable_genes) == 0 && n_genes > 0) {
  message("⚠ No variable genes found. This implies all genes are either high-prevalence (>=95%) or absent.")
  message("  Check if the input matrix is binary (0/1) and correctly parsed.")
  # Diagnostic: print first few rows/cols
  print(mat[1:min(5, n_iso), 1:min(5, n_genes)])
}

# 6. Prevalence Bar Plot
# ------------------------------------------------------------------------------
gene_prev <- tibble(
  gene = colnames(mat),
  n_present = unname(gene_present),
  n_available = unname(gene_available),
  n_missing = unname(gene_missing),
  prevalence = unname(gene_prevalence)
) |>
  arrange(desc(prevalence))

topN <- 40
prev_plot <- ggplot(
  slice_max(gene_prev, prevalence, n = topN),
  aes(reorder(gene, prevalence), prevalence)
) +
  geom_col(fill = "steelblue") +
  coord_flip() +
  scale_y_continuous(labels = scales::percent) +
  labs(
    x = NULL,
    y = "Available isolate calls with gene detected",
    title = "Most prevalent virulence factor genes among VF/WGS-linked isolates",
    subtitle = "Selection criterion: top 40 genes by isolate-level prevalence in the canonical VF matrix",
    caption = sprintf(
      "Data: %s. Cohort: %d VF/WGS-linked E. coli isolates; each gene's prevalence denominator excludes explicitly unavailable calls. Missing is not treated as absence. This overview is descriptive and is not a UTI-vs-Not_UTI association test.",
      FILE_VF_PA, nrow(mat)
    )
  ) +
  theme_bw(base_size = 10) +
  theme(plot.caption = element_text(hjust = 0, size = 7, colour = "grey35"))

ggsave(file.path(DIR_PLOTS_VF, "gene_prevalence_bar.png"), prev_plot, width = 6, height = 8, dpi = 300)

# 7. Variable Gene Heatmap
# ------------------------------------------------------------------------------
if (length(variable_genes) > 0) {
  var_mat <- mat[, variable_genes, drop = FALSE]
  # A deterministic biological/descriptive order is more interpretable than
  # pheatmap's default Euclidean clustering of binary calls. Order genes by
  # observed prevalence, then name; keep participant/timepoint row order.
  var_order <- gene_prev |>
    filter(gene %in% variable_genes) |>
    arrange(desc(prevalence), gene) |>
    pull(gene)
  var_mat <- var_mat[, var_order, drop = FALSE]
  var_mat_plot <- var_mat
  var_mat_plot[is.na(var_mat_plot)] <- 2

  ann_row <- vf |>
    unite(row_id, Participant_id, tp_lab, sep = "_", remove = FALSE) |>
    mutate(Timepoint = factor(tp_lab)) |>
    select(row_id, Timepoint) |>
    column_to_rownames("row_id")
  ann_row <- ann_row[rownames(var_mat_plot), , drop = FALSE]
  if (!identical(rownames(ann_row), rownames(var_mat_plot))) {
    stop("VF heatmap row annotation does not exactly align with the matrix.")
  }

  row_labels <- rownames(var_mat_plot)
  heat_width <- max(9, 0.08 * ncol(var_mat_plot) + 3.5)
  heat_height <- max(10, 0.05 * nrow(var_mat_plot) + 2)
  gene_status_colours <- c("#F2F2F2", "#005AA0", "#6A6A6A")
  gene_status_breaks <- c(-0.5, 0.5, 1.5, 2.5)
  gene_status_legend_breaks <- c(0, 1, 2)
  gene_status_legend_labels <- c("Absent", "Present", "Unavailable")
  heat_title <- paste0(
    "Variable VF gene presence/absence across VF/WGS-linked isolates\n",
    "Genes ordered by observed prevalence; rows ordered by participant and timepoint"
  )

  heat_file_png <- file.path(DIR_PLOTS_VF, "variable_gene_heatmap.png")
  heat_file_pdf <- file.path(DIR_PLOTS_VF, "variable_gene_heatmap.pdf")

  # PNG
  pheatmap(var_mat_plot,
    cluster_rows   = FALSE,
    cluster_cols   = FALSE,
    show_rownames  = TRUE,
    labels_row     = row_labels,
    fontsize       = 8,
    fontsize_row   = 3,
    fontsize_col   = 5,
    annotation_row = ann_row,
    color          = gene_status_colours,
    breaks         = gene_status_breaks,
    legend_breaks  = gene_status_legend_breaks,
    legend_labels  = gene_status_legend_labels,
    main           = heat_title,
    filename       = heat_file_png,
    width          = heat_width,
    height         = heat_height
  )

  # PDF
  pheatmap(var_mat_plot,
    cluster_rows   = FALSE,
    cluster_cols   = FALSE,
    show_rownames  = TRUE,
    labels_row     = row_labels,
    fontsize       = 8,
    fontsize_row   = 3,
    fontsize_col   = 6,
    color          = gene_status_colours,
    breaks         = gene_status_breaks,
    legend_breaks  = gene_status_legend_breaks,
    legend_labels  = gene_status_legend_labels,
    main           = heat_title,
    annotation_row = ann_row,
    filename       = heat_file_pdf,
    width          = heat_width,
    height         = heat_height
  )

  heat_outputs <- c(heat_file_png, heat_file_pdf)
  if (any(!file.exists(heat_outputs)) || any(file.info(heat_outputs)$size <= 0)) {
    stop("One or more VF heatmap outputs are missing or empty.")
  }

  message("✓ Heatmaps generated.")
} else {
  message("⚠ No variable genes found, skipping heatmap.")
}

message("✓ Overview plots complete.")

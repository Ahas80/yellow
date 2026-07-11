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
vf[gene_cols] <- vf[gene_cols] |>
  mutate(across(everything(), ~ replace_na(as.numeric(.x), 0)))

vf <- vf |>
  mutate(
    Participant_id = as.character(Participant_id),
    tp_lab = as.character(tp_lab)
  )

# Maximize per participant-timepoint (handle duplicates/assemblers)
vf <- vf |>
  group_by(across(all_of(meta_cols))) |>
  summarise(across(all_of(gene_cols), max), .groups = "drop") |>
  arrange(across(all_of(meta_cols)))

# Matrix for Heatmap
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
gene_sums <- colSums(mat)

# [BIO] Use 95% threshold instead of 100% (accounts for assembly/annotation gaps)
core_threshold <- 0.95
core_genes <- names(which(gene_sums >= core_threshold * n_iso))

# If strict 100% yields too few, one might relax this (e.g. >= 99%),
# but for now we just report it.
# Variable genes = present in at least one but not all
variable_genes <- names(which(gene_sums > 0 & gene_sums < n_iso))

write_lines(core_genes, file.path(DIR_VF, "core_gene_list.txt"))
write_lines(variable_genes, file.path(DIR_VF, "variable_gene_list.txt"))

message("Core genes (>=95%): ", length(core_genes))
message("Variable genes (<95% & >0%): ", length(variable_genes))
message("Absent genes (0%): ", length(which(gene_sums == 0)))

if (length(variable_genes) == 0 && n_genes > 0) {
  message("⚠ No variable genes found. This implies all genes are either high-prevalence (>=95%) or absent.")
  message("  Check if the input matrix is binary (0/1) and correctly parsed.")
  # Diagnostic: print first few rows/cols
  print(mat[1:min(5, n_iso), 1:min(5, n_genes)])
}

# 6. Prevalence Bar Plot
# ------------------------------------------------------------------------------
gene_prev <- colSums(mat) |>
  enframe(name = "gene", value = "n_iso") |>
  mutate(prevalence = n_iso / nrow(mat)) |>
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
    y = "VF/WGS-linked isolates with gene detected",
    title = "Most prevalent virulence factor genes among VF/WGS-linked isolates",
    subtitle = "Selection criterion: top 40 genes by isolate-level prevalence in the canonical VF matrix",
    caption = sprintf(
      "Data: %s. Denominator: %d VF/WGS-linked E. coli isolates. This overview is descriptive and is not a UTI-vs-Not_UTI association test.",
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
  var_mat <- var_mat[, colSums(var_mat) > 0, drop = FALSE]

  ann_row <- vf |>
    unite(row_id, Participant_id, tp_lab, sep = "_", remove = FALSE) |>
    mutate(Timepoint = factor(tp_lab)) |>
    select(row_id, Timepoint) |>
    column_to_rownames("row_id")

  row_labels <- rownames(var_mat)
  heat_width <- max(9, 0.08 * ncol(var_mat) + 3.5)
  heat_height <- max(10, 0.05 * nrow(var_mat) + 2)
  gene_status_colours <- c("#f2f2f2", "#005aa0")
  gene_status_breaks <- c(-0.01, 0.5, 1.01)
  gene_status_legend_breaks <- c(0, 1)
  gene_status_legend_labels <- c("Absent (light grey)", "Present (blue)")
  heat_title <- "Variable VF gene presence/absence across VF/WGS-linked isolates\nLight grey = absent gene; blue = present gene"

  heat_file_png <- file.path(DIR_PLOTS_VF, "variable_gene_heatmap.png")
  heat_file_pdf <- file.path(DIR_PLOTS_VF, "variable_gene_heatmap.pdf")

  # PNG
  pheatmap(var_mat,
    cluster_rows   = FALSE,
    cluster_cols   = TRUE,
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
  pdf(heat_file_pdf, width = heat_width, height = heat_height)
  pheatmap(var_mat,
    cluster_rows   = FALSE,
    cluster_cols   = TRUE,
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
    annotation_row = ann_row
  )
  dev.off()

  message("✓ Heatmaps generated.")
} else {
  message("⚠ No variable genes found, skipping heatmap.")
}

message("✓ Overview plots complete.")

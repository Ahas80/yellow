#!/usr/bin/env Rscript
# ==============================================================================
# 05_gene_overview_plots.R
# ------------------------------------------------------------------------------
# Role: [Descriptive] - Visualize gene prevalence and variability.
#
# Inputs:
#   - results/vf/vf_pa_all.csv
#
# Outputs:
#   - results/vf/core_gene_list.txt
#   - results/vf/variable_gene_list.txt
#   - plots/vf/gene_prevalence_bar.png
#   - plots/vf/variable_gene_heatmap.png
#   - plots/vf/variable_gene_heatmap.pdf
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

meta_cols <- c("Participant_id", "tp_lab")
gene_cols <- setdiff(names(vf), meta_cols)

# 4. Tidy Matrix
# ------------------------------------------------------------------------------
vf[gene_cols] <- vf[gene_cols] |>
  mutate(across(everything(), ~ replace_na(as.numeric(.x), 0)))

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

message("Core genes (100%): ", length(core_genes))
message("Variable genes (<100% & >0%): ", length(variable_genes))
message("Absent genes (0%): ", length(which(gene_sums == 0)))

if (length(variable_genes) == 0 && n_genes > 0) {
  message("⚠ No variable genes found. This implies all genes are either Core (100%) or Absent (0%).")
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
    x = NULL, y = "Number of Isolates",
    title = "Top 40 Most Prevalent Virulence Factor Genes"
  ) +
  theme_minimal(base_size = 10)

ggsave(file.path(DIR_PLOTS_VF, "gene_prevalence_bar.png"), prev_plot, width = 6, height = 8, dpi = 300)

# 7. Variable Gene Heatmap
# ------------------------------------------------------------------------------
if (length(variable_genes) > 0) {
  var_mat <- mat[, variable_genes, drop = FALSE]
  var_mat <- var_mat[, colSums(var_mat) > 0, drop = FALSE]

  ann_row <- vf |>
    unite(row_id, Participant_id, tp_lab, sep = "_", remove = FALSE) |>
    select(row_id, Participant_id, Timepoint = tp_lab) |>
    column_to_rownames("row_id")

  heat_file_png <- file.path(DIR_PLOTS_VF, "variable_gene_heatmap.png")
  heat_file_pdf <- file.path(DIR_PLOTS_VF, "variable_gene_heatmap.pdf")

  # PNG
  pheatmap(var_mat,
    cluster_rows   = FALSE,
    cluster_cols   = TRUE,
    show_rownames  = FALSE,
    fontsize_col   = 5,
    annotation_row = ann_row,
    main           = "Variable Virulence Gene Presence/Absence",
    filename       = heat_file_png
  )

  # PDF
  pdf(heat_file_pdf, width = 0.18 * ncol(var_mat) + 2, height = 8)
  pheatmap(var_mat,
    cluster_rows   = FALSE,
    cluster_cols   = TRUE,
    show_rownames  = FALSE,
    fontsize_col   = 6,
    main           = "Variable Virulence Gene Presence/Absence",
    annotation_row = ann_row
  )
  dev.off()

  message("✓ Heatmaps generated.")
} else {
  message("⚠ No variable genes found, skipping heatmap.")
}

message("✓ Overview plots complete.")

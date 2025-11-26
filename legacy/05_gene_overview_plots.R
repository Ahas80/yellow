#!/usr/bin/env Rscript
# =============================================================
#  05_gene_overview_plots.R
#  • bar-plot of gene prevalence
#  • heat-map of variable genes (all labels saved to PNG & PDF)
# =============================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(pheatmap)
})

# ---------------------------------------------------------------- paths -------
csv_in   <- "results/vf_pa_all.csv"
plot_dir <- "results/plots"
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------- 1 · load -------
vf <- read_csv(csv_in, show_col_types = FALSE)

meta_cols  <- c("Participant_id", "Timepoint")
gene_cols  <- setdiff(names(vf), meta_cols)

# ----------------------------------------------------- 2 · tidy matrix -------
vf[gene_cols] <- vf[gene_cols] |>
  mutate(across(everything(), ~ replace_na(as.numeric(.x), 0)))

vf <- vf |>
  group_by(across(all_of(meta_cols))) |>
  summarise(across(all_of(gene_cols), max), .groups = "drop") |>
  arrange(across(all_of(meta_cols)))

# strictly-numeric presence/absence matrix
mat <- vf |>
  unite(row_id, Participant_id, Timepoint, sep = "_", remove = FALSE) |>
  column_to_rownames("row_id") |>
  select(all_of(gene_cols)) |>
  as.matrix()

# ------------------------------------------------ 3 · core vs variable --------
core_genes     <- names(which(colSums(mat) == nrow(mat)))
write_lines(core_genes, "results/core_gene_list.txt")

variable_genes <- setdiff(colnames(mat), core_genes)
write_lines(variable_genes, "results/variable_gene_list.txt")


cat("\n===== Core genes (present in *every* isolate) =====\n",
    paste(core_genes, collapse = ", "), "\n\n")
cat("===== Variable genes (present in < 100 % isolates) =====\n",
    paste(variable_genes, collapse = ", "), "\n\n")

# -------------------------------------- 4 · prevalence bar-plot --------------
gene_prev <- colSums(mat) |>
  enframe(name = "gene", value = "n_iso") |>
  mutate(prevalence = n_iso / nrow(mat)) |>
  arrange(desc(prevalence))

topN <- 40
prev_plot <- ggplot(slice_max(gene_prev, prevalence, n = topN),
                    aes(reorder(gene, prevalence), prevalence)) +
  geom_col(fill = "steelblue") +
  coord_flip() +
  scale_y_continuous(labels = scales::percent) +
  labs(x = NULL, y = "Isolates",
       title = paste("Top", topN, "genes by prevalence")) +
  theme_minimal(base_size = 10)

bar_file <- file.path(plot_dir, "gene_prevalence_bar.png")
ggsave(bar_file, prev_plot, width = 6, height = 8, dpi = 300)

# ------------------------------ 5 · variable-gene heat-map --------------------
var_mat <- mat[, variable_genes, drop = FALSE]
var_mat <- var_mat[, colSums(var_mat) > 0, drop = FALSE]   # safety

ann_row <- vf |>
  unite(row_id, Participant_id, Timepoint, sep = "_", remove = FALSE) |>
  select(row_id, Participant_id, Timepoint) |>
  column_to_rownames("row_id")

heat_file_png <- file.path(plot_dir, "variable_gene_heatmap.png")
heat_file_pdf <- file.path(plot_dir, "variable_gene_heatmap.pdf")  # wide PDF

# ---- PNG (bitmap) -----------------------------------------------------------
pheatmap(var_mat,
         cluster_rows   = FALSE,
         cluster_cols   = TRUE,
         show_rownames  = FALSE,
         fontsize_col   = 5,
         annotation_row = ann_row,
         filename       = heat_file_png)

# ---- PDF (vector: every label stays crisp) ----------------------------------
pdf(heat_file_pdf,
    width  = 0.18 * ncol(var_mat) + 2,   # ≈5 mm / gene + margins
    height = 8)
pheatmap(var_mat,
         cluster_rows   = FALSE,
         cluster_cols   = TRUE,
         show_rownames  = FALSE,
         fontsize_col   = 6,
         annotation_row = ann_row)
dev.off()

# --------------------------- 6 · optional on-screen preview -------------------
if (interactive()) {
  message("Showing heat-map on screen …")
  pheatmap(var_mat,
           cluster_rows   = FALSE,
           cluster_cols   = TRUE,
           show_rownames  = FALSE,
           fontsize_col   = 5,
           annotation_row = ann_row,
           main           = "Variable genes across participants & time-points")
}

# ------------------------------------------------------------------ done ------
message("\n✓  Outputs written:\n  • ", bar_file,
        "\n  • ", heat_file_png,
        "\n  • ", heat_file_pdf, "\n")

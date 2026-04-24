#!/usr/bin/env Rscript
# ==============================================================================
# visualize_vf_03_category_profiles.R — Figure 3: VF Category Composition
# ==============================================================================

suppressPackageStartupMessages({
    library(dplyr)
    library(readr)
    library(ggplot2)
    library(tidyr)
})

STATUS_COLORS <- c(ASB = "#4A90D9", UTI = "#D94A4A", Negative = "#888888")
STATUS_ORDER <- c("ASB", "UTI", "Negative")

CAT_COLORS <- c(
    "Adhesion/Fimbriae" = "#E69F00", "Iron acquisition" = "#56B4E9",
    "Toxins" = "#CC79A7", "Capsule/Surface" = "#009E73",
    "Invasion/Evasion" = "#F0E442", "Unassigned" = "#999999"
)

outdir <- "plots/vf"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
sumdir <- "results/vf/fig_summaries"
dir.create(sumdir, recursive = TRUE, showWarnings = FALSE)

# ---- Load ----
df <- read_csv("results/vf/vf_analysis_table.csv", show_col_types = FALSE) %>%
    filter(!is.na(Status)) %>%
    mutate(Status = factor(Status, levels = STATUS_ORDER))

cat_cols <- grep("^cat_", names(df), value = TRUE)
cat_cols <- cat_cols[!grepl("Unassigned_matrix", cat_cols)] # avoid duplicate

cat("Figure 3: Category Profiles | n =", nrow(df), "episodes\n")

# ---- Panel A: Stacked bar (mean category counts by status) ----
cat_long <- df %>%
    select(Participant_id, tp_lab, Status, all_of(cat_cols)) %>%
    pivot_longer(cols = all_of(cat_cols), names_to = "Category", values_to = "count") %>%
    mutate(
        Category = gsub("^cat_", "", Category),
        Category = gsub("_", "/", Category)
    )

cat_means <- cat_long %>%
    group_by(Status, Category) %>%
    summarise(mean_count = mean(count, na.rm = TRUE), .groups = "drop")

p_a <- ggplot(cat_means, aes(x = Status, y = mean_count, fill = Category)) +
    geom_col(position = "stack", width = 0.6, alpha = 0.85) +
    scale_fill_manual(values = CAT_COLORS) +
    labs(
        title = "VF Category Composition by Clinical Status",
        subtitle = "Mean gene count per category per episode",
        x = NULL, y = "Mean VF Gene Count", fill = "VF Category"
    ) +
    theme_minimal(base_size = 12) +
    theme(plot.subtitle = element_text(size = 9, color = "grey40"))

# ---- Panel B: Per-category boxplot (faceted) ----
p_b <- ggplot(cat_long, aes(x = Status, y = count, fill = Status)) +
    geom_boxplot(alpha = 0.5, outlier.size = 0.5) +
    facet_wrap(~Category, scales = "free_y", ncol = 3) +
    scale_fill_manual(values = STATUS_COLORS) +
    labs(
        title = "VF Gene Counts by Category and Status",
        subtitle = "Each panel is one functional category",
        x = NULL, y = "Gene Count"
    ) +
    theme_minimal(base_size = 10) +
    theme(
        legend.position = "none",
        strip.text = element_text(face = "bold", size = 9),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        plot.subtitle = element_text(size = 9, color = "grey40")
    )

# ---- Panel C: Heatmap (median by category × status) ----
cat_medians <- cat_long %>%
    group_by(Status, Category) %>%
    summarise(median_count = median(count, na.rm = TRUE), .groups = "drop")

p_c <- ggplot(cat_medians, aes(x = Status, y = Category, fill = median_count)) +
    geom_tile(color = "white", linewidth = 1) +
    geom_text(aes(label = round(median_count, 1)), size = 4) +
    scale_fill_gradient(low = "white", high = "#2166AC") +
    labs(
        title = "Median VF Gene Count by Category × Status",
        x = NULL, y = NULL, fill = "Median\nCount"
    ) +
    theme_minimal(base_size = 12) +
    theme(axis.text.y = element_text(size = 10))

# ---- Save ----
ggsave(file.path(outdir, "fig3a_category_stacked.png"), p_a, width = 6, height = 5, dpi = 300)
ggsave(file.path(outdir, "fig3b_category_boxplots.png"), p_b, width = 9, height = 6, dpi = 300)
ggsave(file.path(outdir, "fig3c_category_heatmap.png"), p_c, width = 6, height = 4, dpi = 300)
ggsave(file.path(outdir, "fig3a_category_stacked.pdf"), p_a, width = 6, height = 5)

write_csv(cat_means, file.path(sumdir, "fig3_category_means.csv"))
cat("Saved Figure 3 panels to", outdir, "\n")

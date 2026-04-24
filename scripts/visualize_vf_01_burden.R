#!/usr/bin/env Rscript
# ==============================================================================
# visualize_vf_01_burden.R — Figure 1: VF Burden by Status
# ==============================================================================

suppressPackageStartupMessages({
    library(dplyr)
    library(readr)
    library(ggplot2)
    library(tidyr)
})

# ---- Config ----
STATUS_COLORS <- c(ASB = "#4A90D9", UTI = "#D94A4A", Negative = "#888888")
STATUS_ORDER <- c("ASB", "UTI", "Negative")

outdir <- "plots/vf"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
sumdir <- "results/vf/fig_summaries"
dir.create(sumdir, recursive = TRUE, showWarnings = FALSE)

# ---- Load ----
df <- read_csv("results/vf/vf_analysis_table.csv", show_col_types = FALSE) %>%
    filter(!is.na(Status)) %>%
    mutate(Status = factor(Status, levels = STATUS_ORDER))

cat("Figure 1: VF Burden | n =", nrow(df), "episodes,", n_distinct(df$Participant_id), "participants\n")

# ---- Panel A: Episode-level raincloud ----
p_a <- ggplot(df, aes(x = Status, y = vf_count, fill = Status)) +
    geom_violin(alpha = 0.3, width = 0.8, color = NA) +
    geom_boxplot(width = 0.15, outlier.shape = NA, alpha = 0.6) +
    geom_jitter(width = 0.12, alpha = 0.35, size = 1.2, aes(color = Status)) +
    scale_fill_manual(values = STATUS_COLORS) +
    scale_color_manual(values = STATUS_COLORS) +
    labs(
        title = "VF Burden by Clinical Status",
        subtitle = "Episode-level (note: same participant can contribute multiple episodes)",
        x = NULL, y = "VF Gene Count (out of 164)"
    ) +
    annotate("text",
        x = 1:3,
        y = rep(max(df$vf_count) + 5, 3),
        label = paste0("n=", table(df$Status)[STATUS_ORDER]),
        size = 3.2, fontface = "italic"
    ) +
    theme_minimal(base_size = 12) +
    theme(
        legend.position = "none",
        plot.subtitle = element_text(size = 9, color = "grey40")
    )

# ---- Panel B: Participant-aggregated ----
pt_mean <- df %>%
    group_by(Participant_id, Status) %>%
    summarise(mean_vf = mean(vf_count), .groups = "drop")

# Identify participants with BOTH ASB and UTI
both <- pt_mean %>%
    filter(Status %in% c("ASB", "UTI")) %>%
    group_by(Participant_id) %>%
    filter(n() > 1) %>%
    ungroup()

p_b <- ggplot(pt_mean, aes(x = Status, y = mean_vf)) +
    geom_boxplot(aes(fill = Status), alpha = 0.3, width = 0.4, outlier.shape = NA) +
    geom_jitter(aes(color = Status), width = 0.12, alpha = 0.5, size = 1.5) +
    {
        if (nrow(both) > 0) {
            geom_line(
                data = both, aes(group = Participant_id),
                alpha = 0.4, color = "grey30", linewidth = 0.4
            )
        }
    } +
    scale_fill_manual(values = STATUS_COLORS) +
    scale_color_manual(values = STATUS_COLORS) +
    labs(
        title = "Participant-Aggregated VF Burden",
        subtitle = "One mean per participant per status (lines connect paired participants)",
        x = NULL, y = "Mean VF Gene Count"
    ) +
    annotate("text",
        x = 1:3,
        y = rep(max(pt_mean$mean_vf, na.rm = TRUE) + 5, 3),
        label = paste0("n=", table(pt_mean$Status)[STATUS_ORDER]),
        size = 3.2, fontface = "italic"
    ) +
    theme_minimal(base_size = 12) +
    theme(
        legend.position = "none",
        plot.subtitle = element_text(size = 9, color = "grey40")
    )

# ---- Panel C: Faceted by timepoint ----
p_c <- ggplot(df, aes(x = Status, y = vf_count, fill = Status)) +
    geom_boxplot(alpha = 0.5, width = 0.6, outlier.size = 1) +
    geom_jitter(width = 0.1, alpha = 0.3, size = 0.8, aes(color = Status)) +
    facet_wrap(~tp_lab, nrow = 1) +
    scale_fill_manual(values = STATUS_COLORS) +
    scale_color_manual(values = STATUS_COLORS) +
    labs(
        title = "VF Burden by Status and Timepoint",
        subtitle = "Note: All UTI episodes are at Uricult",
        x = NULL, y = "VF Gene Count"
    ) +
    theme_minimal(base_size = 11) +
    theme(
        legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        plot.subtitle = element_text(size = 9, color = "grey40")
    )

# ---- Save ----
ggsave(file.path(outdir, "fig1a_burden_raincloud.png"), p_a, width = 6, height = 5, dpi = 300)
ggsave(file.path(outdir, "fig1b_burden_participant_agg.png"), p_b, width = 6, height = 5, dpi = 300)
ggsave(file.path(outdir, "fig1c_burden_by_timepoint.png"), p_c, width = 10, height = 5, dpi = 300)
ggsave(file.path(outdir, "fig1a_burden_raincloud.pdf"), p_a, width = 6, height = 5)

cat("Saved Figure 1 panels to", outdir, "\n")

# ---- Summary CSV ----
burden_summary <- df %>%
    group_by(Status) %>%
    summarise(
        n_episodes = n(), n_participants = n_distinct(Participant_id),
        mean_vf = round(mean(vf_count), 1), sd_vf = round(sd(vf_count), 1),
        median_vf = median(vf_count),
        q25 = quantile(vf_count, 0.25), q75 = quantile(vf_count, 0.75),
        .groups = "drop"
    )
write_csv(burden_summary, file.path(sumdir, "fig1_burden_summary.csv"))
cat("Summary written to", file.path(sumdir, "fig1_burden_summary.csv"), "\n")

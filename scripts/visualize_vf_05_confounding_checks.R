#!/usr/bin/env Rscript
# ==============================================================================
# visualize_vf_05_confounding_checks.R — Figure 5: Timepoint Confounding
# ==============================================================================

suppressPackageStartupMessages({
    library(dplyr)
    library(readr)
    library(ggplot2)
    library(tidyr)
})

STATUS_COLORS <- c(ASB = "#4A90D9", UTI = "#D94A4A", Negative = "#888888")
STATUS_ORDER <- c("ASB", "UTI", "Negative")

outdir <- "plots/vf"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
sumdir <- "results/vf/fig_summaries"
dir.create(sumdir, recursive = TRUE, showWarnings = FALSE)

# ---- Load ----
df <- read_csv("results/vf/vf_analysis_table.csv", show_col_types = FALSE) %>%
    filter(!is.na(Status)) %>%
    mutate(
        Status = factor(Status, levels = STATUS_ORDER),
        tp_lab = factor(tp_lab, levels = c("T0", "T1", "T2", "Uricult"))
    )

cat("Figure 5: Confounding Checks | n =", nrow(df), "\n")

# ---- Panel A: Balloon/tile plot Status × Timepoint ----
cross <- df %>%
    count(Status, tp_lab, .drop = FALSE) %>%
    mutate(label = ifelse(n == 0, "", as.character(n)))

p_a <- ggplot(cross, aes(x = tp_lab, y = Status)) +
    geom_tile(aes(fill = n), color = "white", linewidth = 1.5) +
    geom_text(aes(label = label), size = 5, fontface = "bold") +
    scale_fill_gradient(low = "white", high = "#2166AC", limits = c(0, NA)) +
    labs(
        title = "Status × Timepoint Cross-Tabulation",
        subtitle = "All 16 UTI episodes cluster at Uricult — a structural confound",
        x = "Timepoint", y = "Clinical Status", fill = "n episodes"
    ) +
    theme_minimal(base_size = 13) +
    theme(
        panel.grid = element_blank(),
        plot.subtitle = element_text(size = 10, color = "grey40")
    )

# ---- Panel B: Scheduled visits only (T0+T1+T2): ASB vs Negative ----
scheduled <- df %>% filter(tp_lab %in% c("T0", "T1", "T2"))

p_b <- ggplot(scheduled, aes(x = Status, y = vf_count, fill = Status)) +
    geom_violin(alpha = 0.3, width = 0.7, color = NA) +
    geom_boxplot(width = 0.15, outlier.shape = NA, alpha = 0.6) +
    geom_jitter(width = 0.1, alpha = 0.3, size = 1, aes(color = Status)) +
    scale_fill_manual(values = STATUS_COLORS) +
    scale_color_manual(values = STATUS_COLORS) +
    labs(
        title = "VF Burden at Scheduled Visits (T0/T1/T2 Only)",
        subtitle = sprintf(
            "Excludes Uricult entirely | n_ASB=%d, n_Neg=%d (no UTI at scheduled visits)",
            sum(scheduled$Status == "ASB"), sum(scheduled$Status == "Negative")
        ),
        x = NULL, y = "VF Gene Count"
    ) +
    annotate("text",
        x = c(1, 3), y = rep(max(scheduled$vf_count) + 3, 2),
        label = c(
            paste0("n=", sum(scheduled$Status == "ASB")),
            paste0("n=", sum(scheduled$Status == "Negative"))
        ),
        size = 3.5, fontface = "italic"
    ) +
    theme_minimal(base_size = 12) +
    theme(
        legend.position = "none",
        plot.subtitle = element_text(size = 9, color = "grey40")
    )

# ---- Panel C: Uricult only ----
uricult <- df %>% filter(tp_lab == "Uricult")

p_c <- ggplot(uricult, aes(x = Status, y = vf_count, fill = Status)) +
    geom_boxplot(alpha = 0.5, width = 0.5, outlier.shape = NA) +
    geom_jitter(width = 0.1, alpha = 0.5, size = 2, aes(color = Status)) +
    scale_fill_manual(values = STATUS_COLORS) +
    scale_color_manual(values = STATUS_COLORS) +
    labs(
        title = "VF Burden at Uricult Only",
        subtitle = sprintf(
            "n_UTI=%d, n_Neg=%d | Constant timepoint eliminates temporal confound",
            sum(uricult$Status == "UTI"), sum(uricult$Status == "Negative")
        ),
        x = NULL, y = "VF Gene Count"
    ) +
    theme_minimal(base_size = 12) +
    theme(
        legend.position = "none",
        plot.subtitle = element_text(size = 9, color = "grey40")
    )

# ---- Panel D: Within-person: participants with both scheduled AND Uricult ----
paired <- df %>%
    mutate(visit_type = ifelse(tp_lab == "Uricult", "Uricult", "Scheduled")) %>%
    group_by(Participant_id) %>%
    filter("Uricult" %in% visit_type & "Scheduled" %in% visit_type) %>%
    ungroup()

if (n_distinct(paired$Participant_id) >= 2) {
    paired_agg <- paired %>%
        group_by(Participant_id, visit_type) %>%
        summarise(
            mean_vf = mean(vf_count), status_at = paste(unique(Status), collapse = "/"),
            .groups = "drop"
        )

    p_d <- ggplot(paired_agg, aes(x = visit_type, y = mean_vf)) +
        geom_line(aes(group = Participant_id), alpha = 0.3, color = "grey50") +
        geom_point(aes(color = status_at), size = 2, alpha = 0.7) +
        labs(
            title = "Within-Person: Scheduled vs Uricult VF Burden",
            subtitle = sprintf(
                "%d participants with both scheduled and Uricult data",
                n_distinct(paired_agg$Participant_id)
            ),
            x = "Visit Type", y = "Mean VF Count", color = "Status at Visit"
        ) +
        theme_minimal(base_size = 12) +
        theme(plot.subtitle = element_text(size = 9, color = "grey40"))
} else {
    p_d <- ggplot() +
        annotate("text",
            x = 0.5, y = 0.5,
            label = "Insufficient paired data"
        ) +
        theme_void()
}

# ---- Save ----
ggsave(file.path(outdir, "fig5a_status_timepoint_tile.png"), p_a, width = 7, height = 4, dpi = 300)
ggsave(file.path(outdir, "fig5b_burden_scheduled_only.png"), p_b, width = 6, height = 5, dpi = 300)
ggsave(file.path(outdir, "fig5c_burden_uricult_only.png"), p_c, width = 5, height = 5, dpi = 300)
ggsave(file.path(outdir, "fig5d_within_person_paired.png"), p_d, width = 7, height = 5, dpi = 300)
ggsave(file.path(outdir, "fig5a_status_timepoint_tile.pdf"), p_a, width = 7, height = 4)

# Summary
write_csv(cross, file.path(sumdir, "fig5_status_timepoint_cross.csv"))
cat("Saved Figure 5 panels to", outdir, "\n")

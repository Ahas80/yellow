#!/usr/bin/env Rscript
# ==============================================================================
# visualize_vf_04_longitudinal_stability.R — Figure 4: Within-Person VF Dynamics
# ==============================================================================

suppressPackageStartupMessages({
    library(dplyr)
    library(readr)
    library(ggplot2)
    library(tidyr)
    library(stringr)
})

TRANS_COLORS <- c(
    "ASB→ASB" = "#4A90D9", "ASB→UTI" = "#D94A4A",
    "ASB→Negative" = "#A0A0A0", "Negative→ASB" = "#7EC8E3",
    "Negative→Negative" = "#C0C0C0", "Negative→UTI" = "#FFB347"
)

outdir <- "plots/vf"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
sumdir <- "results/vf/fig_summaries"
dir.create(sumdir, recursive = TRUE, showWarnings = FALSE)

# ---- Load transitions ----
trans <- read_csv("results/vf/vf_longitudinal_transitions.csv", show_col_types = FALSE)

cat("--- DIAGNOSTICS BEFORE FIXING ---\n")
cat("Rows:", nrow(trans), "\n")
cat("Unique transition types:\n")
print(table(trans$transition_type, useNA = "always"))

# Fix NA genes_gained/lost
trans <- trans %>%
    mutate(
        genes_gained = ifelse(is.na(genes_gained), "", genes_gained),
        genes_lost = ifelse(is.na(genes_lost), "", genes_lost)
    )

# Convert to factor directly without string manipulation, assuming "→" is already UTF-8
trans <- trans %>%
    mutate(transition_type = factor(transition_type,
        levels = c(
            "ASB→ASB", "ASB→UTI", "ASB→Negative",
            "Negative→ASB", "Negative→Negative", "Negative→UTI"
        )
    ))

cat("\n--- DIAGNOSTICS AFTER FIXING ---\n")
cat("Unique transition types:\n")
print(table(trans$transition_type, useNA = "always"))

a2u_check <- trans %>% filter(transition_type == "ASB→UTI")
cat("\nNumber of ASB→UTI rows:", nrow(a2u_check), "\n")
cat("Number of ASB→UTI rows with non-empty genes_gained:", sum(a2u_check$genes_gained != ""), "\n")
if (nrow(a2u_check) > 0) {
    gained_list <- a2u_check$genes_gained[a2u_check$genes_gained != ""]
    gained_genes <- unlist(strsplit(gained_list, ";"))
    cat("Top 10 gained genes in ASB→UTI:\n")
    print(head(sort(table(gained_genes), decreasing = TRUE), 10))
}

gmap <- read_csv("results/vf/gene_map.csv", show_col_types = FALSE)

cat(
    "\nFigure 4: Longitudinal | n =", nrow(trans), "transitions,",
    n_distinct(trans$Participant_id), "participants\n"
)

# ---- Panel A: Jaccard distribution histogram by transition type ----
p_a <- ggplot(trans, aes(x = jaccard_similarity, fill = transition_type)) +
    geom_histogram(binwidth = 0.02, alpha = 0.7, position = "identity", color = "white", linewidth = 0.2) +
    facet_wrap(~transition_type, ncol = 3, scales = "free_y") +
    scale_fill_manual(values = TRANS_COLORS) +
    labs(
        title = "VF Profile Stability: Jaccard Similarity Distribution",
        subtitle = sprintf(
            "96 transitions | %.1f%% show no VF change (Jaccard = 1.0)",
            mean(trans$jaccard_similarity == 1, na.rm = TRUE) * 100
        ),
        x = "Jaccard Similarity", y = "Count"
    ) +
    theme_minimal(base_size = 11) +
    theme(
        legend.position = "none", strip.text = element_text(face = "bold", size = 9),
        plot.subtitle = element_text(size = 9, color = "grey40")
    )

# ---- Panel B: Slopegraph VF count from → to ----
slope_data <- trans %>%
    select(Participant_id, transition_type, vf_count_from, vf_count_to) %>%
    pivot_longer(
        cols = c(vf_count_from, vf_count_to),
        names_to = "time", values_to = "vf_count"
    ) %>%
    mutate(time = factor(ifelse(time == "vf_count_from", "From", "To"), levels = c("From", "To")))

p_b <- ggplot(slope_data, aes(x = time, y = vf_count, group = Participant_id)) +
    geom_line(aes(color = transition_type), alpha = 0.4, linewidth = 0.5) +
    geom_point(aes(color = transition_type), size = 1.2, alpha = 0.6) +
    facet_wrap(~transition_type, ncol = 3) +
    scale_color_manual(values = TRANS_COLORS) +
    labs(
        title = "VF Count Change per Transition",
        subtitle = "Each line = one consecutive-timepoint pair",
        x = NULL, y = "VF Gene Count"
    ) +
    theme_minimal(base_size = 11) +
    theme(
        legend.position = "none", strip.text = element_text(face = "bold", size = 9),
        plot.subtitle = element_text(size = 9, color = "grey40")
    )

# ---- Panel C: Top gained genes in ASB→UTI ----
a2u <- trans %>% filter(transition_type == "ASB→UTI", genes_gained != "")
if (nrow(a2u) > 0) {
    gained_counts <- a2u %>%
        pull(genes_gained) %>%
        strsplit(";") %>%
        unlist() %>%
        as.data.frame() %>%
        setNames("gene") %>%
        count(gene, sort = TRUE) %>%
        head(15) %>%
        left_join(gmap %>% select(Gene, Category), by = c("gene" = "Gene")) %>%
        mutate(
            Category = coalesce(Category, "Unassigned"),
            gene = factor(gene, levels = rev(gene))
        )

    p_c <- ggplot(gained_counts, aes(x = n, y = gene, fill = Category)) +
        geom_col(alpha = 0.8, width = 0.6) +
        geom_text(aes(label = paste0(n, "/", nrow(a2u_check))), hjust = -0.2, size = 3) +
        labs(
            title = "Most Commonly Gained VF Genes in ASB→UTI Transitions",
            subtitle = sprintf("%d / %d ASB→UTI transitions had gene gains", nrow(a2u), nrow(a2u_check)),
            x = "Number of Transitions", y = NULL, fill = "VF Category"
        ) +
        theme_minimal(base_size = 11) +
        theme(plot.subtitle = element_text(size = 9, color = "grey40"))
} else {
    p_c <- ggplot() +
        annotate("text", x = 0.5, y = 0.5, label = "No VF gains in ASB→UTI") +
        theme_void()
}

# ---- Panel D: Heatmap of VF changes per ASB→UTI transition ----
a2u_all <- trans %>% filter(transition_type == "ASB→UTI")
if (nrow(a2u_all) > 0) {
    # Parse gains and losses into a matrix
    all_genes_changed <- unique(c(
        unlist(strsplit(a2u_all$genes_gained, ";")),
        unlist(strsplit(a2u_all$genes_lost, ";"))
    ))
    all_genes_changed <- all_genes_changed[all_genes_changed != ""]

    if (length(all_genes_changed) > 0) {
        heatmap_data <- lapply(seq_len(nrow(a2u_all)), function(i) {
            row <- a2u_all[i, ]
            gained <- strsplit(row$genes_gained, ";")[[1]]
            lost <- strsplit(row$genes_lost, ";")[[1]]
            vals <- setNames(rep(0, length(all_genes_changed)), all_genes_changed)
            vals[gained[gained != ""]] <- 1
            vals[lost[lost != ""]] <- -1
            tibble(participant = row$Participant_id, transition_idx = i, gene = names(vals), change = vals)
        })
        hm_df <- bind_rows(heatmap_data) %>%
            mutate(change_label = case_when(change == 1 ~ "Gained", change == -1 ~ "Lost", TRUE ~ "No change")) %>%
            mutate(participant_idx = paste0(participant, "_", transition_idx))

        # Filter to most variable genes
        gene_var <- hm_df %>%
            group_by(gene) %>%
            summarise(n_changes = sum(change != 0), .groups = "drop") %>%
            filter(n_changes > 0) %>%
            arrange(desc(n_changes))

        hm_plot <- hm_df %>%
            filter(gene %in% head(gene_var$gene, 25)) %>%
            mutate(gene = factor(gene, levels = rev(gene_var$gene[gene_var$gene %in% gene])))

        p_d <- ggplot(hm_plot, aes(x = participant_idx, y = gene, fill = change_label)) +
            geom_tile(color = "white", linewidth = 0.5) +
            scale_fill_manual(values = c("Gained" = "#D94A4A", "Lost" = "#4A90D9", "No change" = "#F0F0F0")) +
            labs(
                title = "VF Gene Changes in ASB→UTI Transitions",
                subtitle = "Each column = one participant's ASB→UTI pair",
                x = "Participant Transition", y = NULL, fill = "Change"
            ) +
            theme_minimal(base_size = 10) +
            theme(
                axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
                plot.subtitle = element_text(size = 9, color = "grey40")
            )
    } else {
        p_d <- ggplot() +
            annotate("text", x = 0.5, y = 0.5, label = "No gene changes") +
            theme_void()
    }
} else {
    p_d <- ggplot() +
        annotate("text", x = 0.5, y = 0.5, label = "No ASB→UTI transitions") +
        theme_void()
}

# ---- Save ----
ggsave(file.path(outdir, "fig4a_jaccard_histograms.png"), p_a, width = 9, height = 6, dpi = 300)
ggsave(file.path(outdir, "fig4b_vf_slopegraph.png"), p_b, width = 9, height = 6, dpi = 300)
ggsave(file.path(outdir, "fig4c_gained_genes_a2u.png"), p_c, width = 7, height = 5, dpi = 300)
ggsave(file.path(outdir, "fig4d_changes_heatmap_a2u.png"), p_d, width = 8, height = 6, dpi = 300)
ggsave(file.path(outdir, "fig4a_jaccard_histograms.pdf"), p_a, width = 9, height = 6)

# Summary
trans_summary <- trans %>%
    group_by(transition_type) %>%
    summarise(
        n = n(), median_jac = median(jaccard_similarity, na.rm = TRUE),
        pct_no_change = mean(!any_vf_change, na.rm = TRUE) * 100,
        mean_gained = mean(n_gained, na.rm = TRUE), mean_lost = mean(n_lost, na.rm = TRUE),
        .groups = "drop"
    )
write_csv(trans_summary, file.path(sumdir, "fig4_transition_summary.csv"))
cat("\nSaved Figure 4 panels to", outdir, "\n")

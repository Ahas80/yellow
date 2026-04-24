#!/usr/bin/env Rscript
# ==============================================================================
# visualize_vf_02_gene_prevalence.R — Figure 2: Gene Prevalence + Enrichment
# ==============================================================================

suppressPackageStartupMessages({
    library(dplyr)
    library(readr)
    library(ggplot2)
    library(tidyr)
    library(stringr)
})

STATUS_COLORS <- c(ASB = "#4A90D9", UTI = "#D94A4A", Negative = "#888888")

outdir <- "plots/vf"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
sumdir <- "results/vf/fig_summaries"
dir.create(sumdir, recursive = TRUE, showWarnings = FALSE)

# ---- Load ----
df <- read_csv("results/vf/vf_analysis_table.csv", show_col_types = FALSE) %>%
    filter(!is.na(Status))

gmap <- read_csv("results/vf/gene_map.csv", show_col_types = FALSE)

gene_cols <- setdiff(
    names(read_csv("results/vf/vf_pa_all.csv", show_col_types = FALSE, n_max = 0)),
    c("Participant_id", "tp_lab")
)

asb_uti <- df %>% filter(Status %in% c("ASB", "UTI"))
n_asb <- sum(asb_uti$Status == "ASB")
n_uti <- sum(asb_uti$Status == "UTI")

cat("Figure 2: Gene Prevalence | ASB n=", n_asb, ", UTI n=", n_uti, "\n")

# ---- Compute prevalence + Fisher tests ----
results <- lapply(gene_cols, function(g) {
    asb_n <- sum(asb_uti[[g]][asb_uti$Status == "ASB"] > 0, na.rm = TRUE)
    uti_n <- sum(asb_uti[[g]][asb_uti$Status == "UTI"] > 0, na.rm = TRUE)
    asb_p <- asb_n / n_asb * 100
    uti_p <- uti_n / n_uti * 100
    delta <- uti_p - asb_p

    # Fisher test
    tab <- table(gene = asb_uti[[g]] > 0, uti = asb_uti$Status == "UTI")
    ft_or <- NA_real_
    ft_p <- NA_real_
    ci_lo <- NA_real_
    ci_hi <- NA_real_
    if (nrow(tab) >= 2 && ncol(tab) >= 2) {
        ft <- tryCatch(fisher.test(tab), error = function(e) NULL)
        if (!is.null(ft)) {
            ft_or <- ft$estimate
            ft_p <- ft$p.value
            ci_lo <- ft$conf.int[1]
            ci_hi <- ft$conf.int[2]
        }
    }

    tibble(
        gene = g, ASB_n = asb_n, UTI_n = uti_n,
        ASB_pct = round(asb_p, 1), UTI_pct = round(uti_p, 1),
        delta = round(delta, 1), OR = round(ft_or, 3),
        p = round(ft_p, 4), CI_lo = round(ci_lo, 3), CI_hi = round(ci_hi, 3)
    )
})

res <- bind_rows(results) %>%
    filter(!is.na(OR)) %>%
    mutate(p_adj = round(p.adjust(p, method = "BH"), 4)) %>%
    left_join(gmap %>% select(Gene, Category), by = c("gene" = "Gene")) %>%
    mutate(Category = coalesce(Category, "Unassigned"))

write_csv(res, file.path(sumdir, "fig2_gene_prevalence_enrichment.csv"))

# ---- Panel A: Lollipop of top 20 by |Δ| ----
top20 <- res %>%
    slice_max(abs(delta), n = 20) %>%
    mutate(
        gene = factor(gene, levels = gene[order(delta)]),
        direction = ifelse(delta > 0, "Higher in UTI", "Higher in ASB")
    )

p_a <- ggplot(top20, aes(x = delta, y = gene, fill = direction)) +
    geom_col(alpha = 0.7, width = 0.6) +
    geom_vline(xintercept = 0, linewidth = 0.5) +
    geom_text(
        aes(
            label = sprintf("%d/%d vs %d/%d", UTI_n, n_uti, ASB_n, n_asb),
            x = ifelse(delta > 0, delta + 1, delta - 1)
        ),
        hjust = ifelse(top20$delta > 0, 0, 1), size = 2.3
    ) +
    scale_fill_manual(values = c("Higher in UTI" = "#D94A4A", "Higher in ASB" = "#4A90D9")) +
    labs(
        title = "Top 20 VF Genes by Prevalence Difference (UTI − ASB)",
        subtitle = sprintf("n_ASB = %d, n_UTI = %d | All tests exploratory (repeated measures)", n_asb, n_uti),
        x = "Δ Prevalence (UTI − ASB, percentage points)", y = NULL, fill = NULL
    ) +
    theme_minimal(base_size = 11) +
    theme(
        legend.position = "bottom",
        plot.subtitle = element_text(size = 9, color = "grey40"),
        axis.text.y = element_text(size = 8)
    )

# ---- Panel B: Forest plot of ORs ----
top15 <- res %>%
    filter(!is.na(OR), is.finite(OR), OR > 0) %>%
    slice_min(p, n = 15) %>%
    mutate(
        gene = factor(gene, levels = gene[order(OR)]),
        sig_label = case_when(
            p_adj < 0.05 ~ "*",
            p < 0.05 ~ "†",
            TRUE ~ ""
        )
    )

p_b <- ggplot(top15, aes(x = OR, y = gene)) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "grey50") +
    geom_errorbarh(aes(xmin = CI_lo, xmax = pmin(CI_hi, 50)),
        height = 0.25, color = "grey40"
    ) +
    geom_point(aes(color = Category), size = 2.5) +
    geom_text(aes(label = sprintf("p=%.3f%s", p, sig_label)),
        hjust = -0.2, size = 2.5, nudge_x = 0.1
    ) +
    scale_x_log10() +
    labs(
        title = "Fisher's Exact Test: VF Gene Enrichment in UTI vs ASB",
        subtitle = "Exploratory, unadjusted | † nominal p<0.05; * BH-adjusted p<0.05",
        x = "Odds Ratio (log scale)", y = NULL, color = "VF Category"
    ) +
    theme_minimal(base_size = 11) +
    theme(
        legend.position = "right",
        plot.subtitle = element_text(size = 9, color = "grey40")
    )

# ---- Save ----
ggsave(file.path(outdir, "fig2a_gene_prevalence_lollipop.png"), p_a, width = 8, height = 7, dpi = 300)
ggsave(file.path(outdir, "fig2b_enrichment_forest.png"), p_b, width = 8, height = 6, dpi = 300)
ggsave(file.path(outdir, "fig2a_gene_prevalence_lollipop.pdf"), p_a, width = 8, height = 7)
ggsave(file.path(outdir, "fig2b_enrichment_forest.pdf"), p_b, width = 8, height = 6)

cat("Saved Figure 2 panels to", outdir, "\n")

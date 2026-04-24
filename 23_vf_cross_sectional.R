#!/usr/bin/env Rscript
# ==============================================================================
# 23_vf_cross_sectional.R
# ==============================================================================
#
# GOAL:
#   Cross-sectional comparison of virulence factor (VF) profiles between
#   clinical states: ASB vs UTI.  Produces burden summaries, gene-level
#   prevalence tables, category-level enrichment, and exploratory Fisher
#   exact tests.
#
# WHY THIS SCRIPT EXISTS:
#   The project's central question is whether bacteria causing UTI carry
#   a different VF arsenal from those causing ASB in elderly nursing home
#   residents.  This script provides the core descriptive and exploratory
#   statistical tables needed to answer that question.
#
#   Previous iterations of this analysis were scattered across 03_plotting.R
#   (quick boxplot), 04_gene_breakdown.R (focus-gene GLMM), and the
#   un-numbered compute_vf_abstract_stats.R (monolithic).  This script
#   consolidates the cross-sectional logic into a single, well-structured
#   numbered pipeline step.
#
# IMPORTANT STATISTICAL CAVEAT:
#   Fisher exact tests in this script are EXPLORATORY.  Our cohort has
#   repeated measures (the same participant contributes multiple episodes),
#   which violates the independence assumption of Fisher's test.
#   Formal inference accounting for within-participant correlation is done
#   via GLMM in 14_genotype_phenotype_model.R.  The Fisher tests here are
#   useful screening tools only.
#
# INPUTS:
#   - results/vf/vf_analysis_ready.csv     (from 22_vf_build_analysis_dataset.R)
#   - results/vf/gene_map.csv              (Gene → Category mapping)
#
# OUTPUTS (all in results/vf/):
#   - vf_burden_by_status.csv              Full-cohort burden summary
#   - vf_burden_by_status_stratified.csv   Burden by depth (≥2/≥3/≥4 timepoints)
#   - vf_category_burden_by_status.csv     Category-level burden
#   - vf_gene_prevalence_by_status.csv     Per-gene prevalence in each status
#   - vf_gene_prevalence_stratified.csv    Prevalence by depth
#   - vf_fisher_exploratory.csv            Fisher exact tests (full cohort)
#   - vf_enrichment_stratified.csv         Fisher tests by depth
#   - vf_category_enrichment.csv           Category-level Fisher tests
#   - vf_cross_sectional_summary.txt       Human-readable abstract-ready summary
#
# PLOTS (in plots/vf/):
#   - vf_burden_boxplot.png                Boxplot of VF burden (ASB vs UTI)
#   - vf_category_barplot.png              Median category counts by status
#
# DEPTH STRATIFICATION:
#   All analyses are automatically run for four cohort definitions:
#     "all"    — all episodes regardless of follow-up depth
#     "≥2 tp"  — only participants with ≥2 clinical timepoints
#     "≥3 tp"  — only participants with ≥3 timepoints
#     "≥4 tp"  — only participants with ≥4 timepoints
#   This prevents the need for separate scripts for different depth cutoffs.
#
# SUPERSEDES:
#   - Sections 3–4 of compute_vf_abstract_stats.R (burden, prevalence, Fisher)
#   - Sections C1–C3 of compute_vf_stratified_by_depth.R
# ==============================================================================

source("00_config.R")
source("R/plot_helpers.R")  # Canonical infection status colours
suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
  library(stringr)
  library(purrr)    # for compact()
  library(ggplot2)
})

msg("Starting 23_vf_cross_sectional.R")

# ==============================================================================
# 1. LOAD DATA
# ==============================================================================
# We load the canonical analysis-ready dataset produced by 22_.
# This contains one row per episode with VF gene P/A, status, ST, and burden.

ready_file <- file.path(DIR_VF, "vf_analysis_ready.csv")
if (!file.exists(ready_file)) stop("Missing ", ready_file, ". Run 22_vf_build_analysis_dataset.R first.")
vf_ready <- read_csv(ready_file, show_col_types = FALSE) %>%
  mutate(Participant_id = as.character(Participant_id))

gene_map_file <- file.path(DIR_VF, "gene_map.csv")
gene_map <- read_csv(gene_map_file, show_col_types = FALSE) %>%
  mutate(Gene = as.character(Gene),
         Category = coalesce(as.character(Category), "Unassigned"))

# Identify column groups:
#   skip_cols  = metadata columns (not VF genes)
#   cat_cols   = category-level count columns (e.g., cat_Iron_acquisition)
#   gene_cols  = individual VF gene columns (binary 0/1)
skip_cols <- c("Participant_id", "tp_lab", "Infection_Status", "Batch", "ST",
               "vf_count_total", "is_ecoli", "n_timepoints")
cat_cols  <- grep("^cat_", names(vf_ready), value = TRUE)
gene_cols <- setdiff(names(vf_ready), c(skip_cols, cat_cols))

msg("Loaded: %d rows, %d gene columns, %d participants",
    nrow(vf_ready), length(gene_cols), n_distinct(vf_ready$Participant_id))

# ==============================================================================
# 2. HELPER: COMPUTE ALL CROSS-SECTIONAL METRICS FOR ONE COHORT
# ==============================================================================
# This function encapsulates the full cross-sectional analysis so it can be
# called once per depth stratum without code duplication.

compute_cross_sectional <- function(data, cohort_label, gene_cols, gene_map) {

  # Filter to episodes that have a clinical status assignment
  with_status <- data %>% filter(!is.na(Infection_Status))

  # ------------------------------------------------------------------
  # C1: VF BURDEN BY STATUS
  # ------------------------------------------------------------------
  # For each clinical status (ASB, UTI, Negative), compute descriptive
  # statistics for total VF gene count.  These are the core numbers
  # reported in the abstract (e.g., "median VF burden was 80 in ASB
  # vs 80.5 in UTI").
  burden <- with_status %>%
    group_by(Infection_Status) %>%
    summarise(
      n_episodes     = n(),
      n_participants = n_distinct(Participant_id),
      mean_vf   = round(mean(vf_count_total), 1),
      sd_vf     = round(sd(vf_count_total), 1),
      median_vf = median(vf_count_total),
      q25_vf    = quantile(vf_count_total, 0.25),
      q75_vf    = quantile(vf_count_total, 0.75),
      min_vf    = min(vf_count_total),
      max_vf    = max(vf_count_total),
      .groups   = "drop"
    ) %>%
    mutate(cohort = cohort_label)

  # ------------------------------------------------------------------
  # C2: GENE-LEVEL PREVALENCE
  # ------------------------------------------------------------------
  # For each VF gene, compute the percentage of episodes where it is
  # present, separately for ASB, UTI, and Negative.
  # The delta (UTI% - ASB%) highlights genes that are differentially
  # prevalent between symptomatic and asymptomatic infections.
  asb_uti <- with_status %>% filter(Infection_Status %in% c("ASB", "UTI"))
  n_asb <- sum(asb_uti$Infection_Status == "ASB")
  n_uti <- sum(asb_uti$Infection_Status == "UTI")

  gene_prev <- lapply(gene_cols, function(g) {
    df <- with_status %>%
      group_by(Infection_Status) %>%
      summarise(n_present = sum(.data[[g]] > 0, na.rm = TRUE),
                n_total = n(), .groups = "drop") %>%
      mutate(pct = round(n_present / n_total * 100, 1))

    # Helper to safely extract a value for a given status
    get_val <- function(status, col) {
      r <- df %>% filter(Infection_Status == status)
      if (nrow(r) == 0) return(NA)
      r[[col]]
    }

    tibble(
      gene = g,
      ASB_n = get_val("ASB", "n_present"), ASB_N = get_val("ASB", "n_total"),
      ASB_pct = get_val("ASB", "pct"),
      UTI_n = get_val("UTI", "n_present"), UTI_N = get_val("UTI", "n_total"),
      UTI_pct = get_val("UTI", "pct"),
      Neg_n = get_val("Negative", "n_present"), Neg_N = get_val("Negative", "n_total"),
      Neg_pct = get_val("Negative", "pct")
    )
  }) %>%
    bind_rows() %>%
    mutate(delta_UTI_minus_ASB = round(UTI_pct - ASB_pct, 1)) %>%
    left_join(gene_map %>% select(Gene, Category), by = c("gene" = "Gene")) %>%
    mutate(Category = coalesce(Category, "Unassigned"),
           cohort = cohort_label) %>%
    arrange(desc(abs(delta_UTI_minus_ASB)))

  # ------------------------------------------------------------------
  # C3: EXPLORATORY FISHER EXACT TESTS (UTI vs ASB)
  # ------------------------------------------------------------------
  # For each gene, test whether its prevalence differs significantly
  # between UTI and ASB episodes using Fisher's exact test.
  #
  # ⚠ These tests assume independent observations.  Since participants
  # contribute multiple episodes, this assumption is VIOLATED.  These
  # results are screening tools only.  Formal inference uses GLMM in
  # script 14_genotype_phenotype_model.R.
  enrichment <- NULL
  if (n_uti >= 2) {
    enrichment <- lapply(gene_cols, function(g) {
      tab <- table(gene_present = asb_uti[[g]] > 0,
                   is_uti       = asb_uti$Infection_Status == "UTI")
      # Need a proper 2×2 table (gene must have variation)
      if (nrow(tab) < 2 || ncol(tab) < 2) return(NULL)
      ft <- tryCatch(fisher.test(tab), error = function(e) NULL)
      if (is.null(ft)) return(NULL)
      tibble(gene = g,
             OR = round(ft$estimate, 3),
             p_value = round(ft$p.value, 4),
             CI_lower = round(ft$conf.int[1], 3),
             CI_upper = round(ft$conf.int[2], 3))
    }) %>%
      bind_rows()

    if (nrow(enrichment) > 0) {
      # Apply Benjamini-Hochberg correction for multiple testing
      enrichment <- enrichment %>%
        mutate(p_adj_BH = round(p.adjust(p_value, method = "BH"), 4),
               note = "EXPLORATORY: repeated measures violate independence") %>%
        left_join(gene_map %>% select(Gene, Category), by = c("gene" = "Gene")) %>%
        mutate(Category = coalesce(Category, "Unassigned"),
               cohort = cohort_label) %>%
        arrange(p_value)
    }
  }

  list(burden = burden, gene_prev = gene_prev, enrichment = enrichment)
}

# ==============================================================================
# 3. RUN FOR FULL COHORT + STRATIFIED BY TIMEPOINT DEPTH
# ==============================================================================
# We run the same analysis four times: once for all participants, and once
# each for participants with ≥2, ≥3, and ≥4 timepoints.
# Stratification matters because participants with more timepoints are
# likely the ones with persistent bacteriuria — restricting to them gives
# a more "longitudinally informative" subset.

cohorts <- list(
  "all"    = vf_ready,
  ">=2 tp" = vf_ready %>% filter(!is.na(n_timepoints) & n_timepoints >= 2),
  ">=3 tp" = vf_ready %>% filter(!is.na(n_timepoints) & n_timepoints >= 3),
  ">=4 tp" = vf_ready %>% filter(!is.na(n_timepoints) & n_timepoints >= 4)
)

results <- lapply(names(cohorts), function(label) {
  msg("  Computing cross-sectional for cohort: %s (%d rows)", label, nrow(cohorts[[label]]))
  compute_cross_sectional(cohorts[[label]], label, gene_cols, gene_map)
})
names(results) <- names(cohorts)

# ==============================================================================
# 4. CATEGORY-LEVEL ENRICHMENT (full cohort only)
# ==============================================================================
# In addition to testing individual genes, we test whether VF *categories*
# (e.g., "Iron acquisition" as a group) are enriched in UTI vs ASB.
# This increases statistical power by grouping functionally related genes,
# and is more biologically interpretable.

cat_cols_ready <- grep("^cat_", names(vf_ready), value = TRUE)
asb_uti_full   <- vf_ready %>%
  filter(!is.na(Infection_Status), Infection_Status %in% c("ASB", "UTI"))

cat_enrichment <- NULL
if (length(cat_cols_ready) > 0 && sum(asb_uti_full$Infection_Status == "UTI") >= 2) {
  cat_enrichment <- lapply(cat_cols_ready, function(cc) {
    # Binarise: does this episode have ≥1 gene from this category?
    asb_uti_full$cat_present <- asb_uti_full[[cc]] > 0
    tab <- table(asb_uti_full$cat_present, asb_uti_full$Infection_Status == "UTI")
    if (nrow(tab) < 2 || ncol(tab) < 2) return(NULL)
    ft <- tryCatch(fisher.test(tab), error = function(e) NULL)
    if (is.null(ft)) return(NULL)

    # Also run a Wilcoxon rank-sum test on the raw counts (non-binarised)
    # to test whether the NUMBER of genes in this category differs by status.
    wt <- tryCatch(
      wilcox.test(as.formula(paste0("`", cc, "` ~ Infection_Status")),
                  data = asb_uti_full),
      error = function(e) NULL
    )

    tibble(
      category = cc,
      fisher_OR    = round(ft$estimate, 3),
      fisher_p     = round(ft$p.value, 4),
      wilcox_p     = if (!is.null(wt)) round(wt$p.value, 4) else NA_real_,
      median_ASB   = median(asb_uti_full[[cc]][asb_uti_full$Infection_Status == "ASB"]),
      median_UTI   = median(asb_uti_full[[cc]][asb_uti_full$Infection_Status == "UTI"])
    )
  }) %>% bind_rows()

  if (!is.null(cat_enrichment) && nrow(cat_enrichment) > 0) {
    cat_enrichment <- cat_enrichment %>%
      mutate(fisher_p_adj_BH = round(p.adjust(fisher_p, method = "BH"), 4),
             note = "EXPLORATORY: repeated measures not accounted for") %>%
      arrange(fisher_p)
  }
}

# ==============================================================================
# 5. CATEGORY-LEVEL BURDEN TABLE
# ==============================================================================
# Descriptive summary: for each clinical status, what is the median and mean
# count of genes in each VF category?  This feeds directly into category
# barplots and manuscript Table 2-style summaries.

cat_burden <- vf_ready %>%
  filter(!is.na(Infection_Status)) %>%
  group_by(Infection_Status) %>%
  summarise(
    across(all_of(cat_cols_ready),
           list(median = ~median(., na.rm = TRUE),
                mean   = ~round(mean(., na.rm = TRUE), 1)),
           .names = "{.col}__{.fn}"),
    .groups = "drop"
  )

# ==============================================================================
# 6. WRITE OUTPUTS
# ==============================================================================

msg("Writing outputs...")

# Full-cohort outputs
write_csv(results[["all"]]$burden,    file.path(DIR_VF, "vf_burden_by_status.csv"))
write_csv(results[["all"]]$gene_prev, file.path(DIR_VF, "vf_gene_prevalence_by_status.csv"))
if (!is.null(results[["all"]]$enrichment) && nrow(results[["all"]]$enrichment) > 0) {
  write_csv(results[["all"]]$enrichment, file.path(DIR_VF, "vf_fisher_exploratory.csv"))
}

# Stratified outputs (all depth levels combined into single CSVs)
all_burden   <- bind_rows(lapply(results, `[[`, "burden"))
all_prev     <- bind_rows(lapply(results, `[[`, "gene_prev"))
all_enriched <- bind_rows(compact(lapply(results, `[[`, "enrichment")))

write_csv(all_burden,   file.path(DIR_VF, "vf_burden_by_status_stratified.csv"))
write_csv(all_prev,     file.path(DIR_VF, "vf_gene_prevalence_stratified.csv"))
if (nrow(all_enriched) > 0) {
  write_csv(all_enriched, file.path(DIR_VF, "vf_enrichment_stratified.csv"))
}

# Category-level outputs
write_csv(cat_burden, file.path(DIR_VF, "vf_category_burden_by_status.csv"))
if (!is.null(cat_enrichment) && nrow(cat_enrichment) > 0) {
  write_csv(cat_enrichment, file.path(DIR_VF, "vf_category_enrichment.csv"))
}

# ==============================================================================
# 7. PLOTS
# ==============================================================================

ensure_dir(DIR_PLOTS_VF)

# PLOT 1: VF Burden Boxplot (ASB vs UTI)
#   This is the "headline" figure: side-by-side box/jitter plot showing
#   that total VF burden is similar between ASB and UTI.
plot_data <- vf_ready %>%
  filter(!is.na(Infection_Status), Infection_Status %in% c("ASB", "UTI"))

if (nrow(plot_data) > 0) {
  p_burden <- ggplot(plot_data, aes(x = Infection_Status, y = vf_count_total,
                                     fill = Infection_Status)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.7, width = 0.5) +
    geom_jitter(width = 0.15, alpha = 0.4, size = 1.5) +
    scale_fill_infection() +
    labs(title = "VF Gene Burden by Clinical Status",
         subtitle = sprintf("ASB n=%d, UTI n=%d | Wilcoxon exploratory",
                            sum(plot_data$Infection_Status == "ASB"),
                            sum(plot_data$Infection_Status == "UTI")),
         x = "Episode Type", y = "Total VF Genes Detected") +
    theme_minimal(base_size = 12) +
    theme(legend.position = "none")
  ggsave(file.path(DIR_PLOTS_VF, "vf_burden_boxplot.png"), p_burden,
         width = 5, height = 5, dpi = 300)
}

# PLOT 2: Category Barplot
#   Shows median count of genes per VF category, side by side for ASB and UTI.
#   Useful for identifying which *functional groups* might differ.
if (length(cat_cols_ready) > 0) {
  cat_long <- vf_ready %>%
    filter(!is.na(Infection_Status), Infection_Status %in% c("ASB", "UTI")) %>%
    select(Participant_id, Infection_Status, all_of(cat_cols_ready)) %>%
    pivot_longer(cols = all_of(cat_cols_ready),
                 names_to = "Category", values_to = "count") %>%
    mutate(Category = gsub("^cat_", "", Category),
           Category = gsub("_", "/", Category))

  cat_summary <- cat_long %>%
    group_by(Category, Infection_Status) %>%
    summarise(median_count = median(count, na.rm = TRUE), .groups = "drop")

  p_cat <- ggplot(cat_summary,
                  aes(x = reorder(Category, -median_count),
                      y = median_count, fill = Infection_Status)) +
    geom_col(position = "dodge", width = 0.6) +
    scale_fill_infection() +
    labs(title = "Median VF Category Counts by Status",
         x = NULL, y = "Median Gene Count", fill = "Episode Type") +
    theme_minimal(base_size = 11) +
    theme(axis.text.x = element_text(angle = 35, hjust = 1)) +
    coord_cartesian(clip = "off")
  ggsave(file.path(DIR_PLOTS_VF, "vf_category_barplot.png"), p_cat,
         width = 7, height = 5, dpi = 300)
}

# ==============================================================================
# 8. HUMAN-READABLE SUMMARY
# ==============================================================================
# This text file is designed to be copy-pasted or referenced directly when
# writing the abstract or results section.  It summarises the key numbers
# in plain language.

summary_lines <- character()
sl <- function(...) summary_lines <<- c(summary_lines, sprintf(...))

sl("=== VF CROSS-SECTIONAL SUMMARY ===")
sl("Generated: %s", format(Sys.time()))
sl("")

b <- results[["all"]]$burden
for (i in seq_len(nrow(b))) {
  r <- b[i, ]
  sl("%s: n=%d (%d participants), median=%.0f (IQR %.0f–%.0f), mean=%.1f±%.1f",
     r$Infection_Status, r$n_episodes, r$n_participants,
     r$median_vf, r$q25_vf, r$q75_vf, r$mean_vf, r$sd_vf)
}

enr <- results[["all"]]$enrichment
if (!is.null(enr) && nrow(enr) > 0) {
  sl("")
  sl("EXPLORATORY Fisher exact (UTI vs ASB, full cohort):")
  sl("  Genes with p < 0.05: %d", sum(enr$p_value < 0.05, na.rm = TRUE))
  sl("  Genes with p_adj(BH) < 0.05: %d", sum(enr$p_adj_BH < 0.05, na.rm = TRUE))
  top5 <- head(enr, 5)
  for (i in seq_len(nrow(top5))) {
    r <- top5[i, ]
    sl("  %s (%s): OR=%.2f, p=%.4f, p_adj=%.4f",
       r$gene, r$Category, r$OR, r$p_value, r$p_adj_BH)
  }
  sl("")
  sl("⚠  WARNING: Fisher tests assume independent observations.")
  sl("   Episodes from the same participant are NOT independent.")
  sl("   Use 14_genotype_phenotype_model.R (GLMM) for formal inference.")
}

writeLines(summary_lines, file.path(DIR_VF, "vf_cross_sectional_summary.txt"))
msg("Summary written to %s", file.path(DIR_VF, "vf_cross_sectional_summary.txt"))

msg("✓ 23_vf_cross_sectional.R complete.")

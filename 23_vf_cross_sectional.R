#!/usr/bin/env Rscript
# ==============================================================================
# 23_vf_cross_sectional.R
# ==============================================================================
#
# GOAL:
#   Cross-sectional comparison of virulence factor (VF) profiles between
#   primary UTI status groups: UTI vs Not_UTI. Produces burden summaries, gene-level
#   prevalence tables, category-level enrichment, and exploratory Fisher
#   exact tests.
#
# WHY THIS SCRIPT EXISTS:
#   The project's central question is whether bacteria causing UTI carry
#   a different VF arsenal from all Not_UTI episodes in elderly nursing home
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
#   - vf_burden_boxplot.png                Boxplot of VF burden (UTI vs Not_UTI)
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
source("R/plot_helpers.R")  # Primary UTI status colours
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
# VF visualisation shared helpers
# ==============================================================================

STATUS_LEVELS <- c("UTI", "Not_UTI", "Unknown")

status_for_plot <- function(x) {
  x <- as.character(x)
  x[is.na(x) | x == ""] <- "Unknown"
  x[!x %in% STATUS_LEVELS] <- "Unknown"
  factor(x, levels = STATUS_LEVELS)
}

status_count_text <- function(data) {
  counts <- data %>%
    mutate(.status = status_for_plot(Infection_Status)) %>%
    count(.status, name = "n") %>%
    filter(n > 0) %>%
    mutate(label = paste0(as.character(.status), " n=", n))
  paste(counts$label, collapse = "; ")
}

vf_caption <- function(input_file, data, analysis_unit,
                       p_value_note = "P-values are descriptive/exploratory unless derived from mixed-effects models with Participant_id as a random effect.",
                       multiple_testing_note = NULL,
                       extra_note = NULL) {
  repeated <- if ("Participant_id" %in% names(data)) {
    any(table(data$Participant_id) > 1)
  } else {
    FALSE
  }
  n_uti <- if ("Infection_Status" %in% names(data)) {
    sum(data$Infection_Status == "UTI", na.rm = TRUE)
  } else {
    NA_integer_
  }
  paste(
    sprintf("Data: %s.", input_file),
    sprintf("Denominator: n=%d VF/WGS-linked E. coli isolates from %d participants (%s).",
            nrow(data), n_distinct(data$Participant_id), status_count_text(data)),
    sprintf("Level of analysis: %s.", analysis_unit),
    if (repeated) "Residents may contribute repeated isolates; isolate-level tests can be pseudoreplicated." else "Each participant contributes one plotted value.",
    p_value_note,
    multiple_testing_note %||% "",
    sprintf("UTI denominator is small (n=%d), so UTI-vs-Not_UTI contrasts are underpowered.", n_uti),
    "Not_UTI is heterogeneous; subgroup/legacy outputs should be used for sensitivity interpretation.",
    "ST/lineage, batch, timepoint, and event-driven sampling may confound apparent VF-status associations.",
    extra_note %||% "",
    sep = " "
  ) %>% str_squish()
}

plot_theme_vf <- function(base_size = 11) {
  theme_bw(base_size = base_size) +
    theme(
      plot.caption = element_text(hjust = 0, size = base_size - 3, colour = "grey35"),
      plot.subtitle = element_text(colour = "grey25"),
      legend.position = "bottom"
    )
}

stop_if_stale <- function(target, upstream, target_label, upstream_label) {
  if (!file.exists(target) || !file.exists(upstream)) return(invisible(FALSE))
  target_mtime <- file.info(target)$mtime
  upstream_mtime <- file.info(upstream)$mtime
  if (!is.na(target_mtime) && !is.na(upstream_mtime) && target_mtime < upstream_mtime) {
    if (basename(target) == "vf_analysis_ready.csv" && basename(upstream) == "vf_pa_all.csv") {
      target_probe <- read_csv(target, show_col_types = FALSE, n_max = Inf)
      upstream_probe <- read_csv(upstream, show_col_types = FALSE, n_max = Inf)
      target_genes <- canonical_vf_gene_cols(names(target_probe), vf_pa_file = upstream)
      upstream_genes <- canonical_vf_gene_cols(names(upstream_probe), vf_pa_file = upstream)
      if (nrow(target_probe) == nrow(upstream_probe) && setequal(target_genes, upstream_genes)) {
        msg("WARNING: %s timestamp is older than %s, but row count and gene set match; continuing.",
            target_label, upstream_label)
        return(invisible(FALSE))
      }
    }
    stop(sprintf(
      "%s is older than %s. Re-run 22_vf_build_analysis_dataset.R before script 23 so plots use the current VF matrix.\n  %s: %s\n  %s: %s",
      target_label, upstream_label, target_label, format(target_mtime),
      upstream_label, format(upstream_mtime)
    ))
  }
}

format_category_label <- function(x) {
  recode(
    gsub("^cat_", "", x),
    "Adhesion_Fimbriae" = "Adhesion / fimbriae",
    "Iron_acquisition" = "Iron acquisition",
    "Toxins" = "Toxins",
    "Capsule_Surface" = "Capsule / surface structures",
    "Invasion_Evasion" = "Invasion / immune evasion",
    "Unassigned" = "Unassigned VFDB hits",
    "Unassigned_matrix" = "Unassigned VFDB hits not in gene map",
    .default = gsub("_", " ", gsub("^cat_", "", x))
  )
}

# ==============================================================================
# 1. LOAD DATA
# ==============================================================================
# We load the canonical analysis-ready dataset produced by 22_.
# This contains one row per episode with VF gene P/A, status, ST, and burden.

ready_file <- FILE_VF_READY
if (!file.exists(ready_file)) stop("Missing ", ready_file, ". Run 22_vf_build_analysis_dataset.R first.")
stop_if_stale(ready_file, FILE_VF_PA, "vf_analysis_ready.csv", "vf_pa_all.csv")
vf_ready <- read_csv(ready_file, show_col_types = FALSE) %>%
  prefer_primary_uti_status() %>%
  apply_manual_sample_curation(context = "23_vf_ready") %>%
  filter_primary_genomics() %>%
  mutate(Participant_id = as.character(Participant_id),
         Infection_Status_plot = status_for_plot(Infection_Status))

gene_map_file <- file.path(DIR_VF, "gene_map.csv")
gene_map <- read_csv(gene_map_file, show_col_types = FALSE) %>%
  mutate(Gene = as.character(Gene),
         Category = coalesce(as.character(Category), "Unassigned"))

# Individual VF gene columns are defined by the canonical P/A matrix from
# 02_gene_presence_analysis.R. Do not infer genes by excluding metadata here:
# vf_analysis_ready.csv also contains clinical provenance and score columns.
cat_cols  <- grep("^cat_", names(vf_ready), value = TRUE)
gene_cols <- canonical_vf_gene_cols(names(vf_ready))

msg("Loaded: %d rows, %d gene columns, %d participants",
    nrow(vf_ready), length(gene_cols), n_distinct(vf_ready$Participant_id))
msg("VF-ready denominator by status: %s", status_count_text(vf_ready))
msg("Repeated-measures warning: %d/%d participants contribute >1 VF-ready isolate.",
    sum(table(vf_ready$Participant_id) > 1), n_distinct(vf_ready$Participant_id))

# ==============================================================================
# 2. HELPER: COMPUTE ALL CROSS-SECTIONAL METRICS FOR ONE COHORT
# ==============================================================================
# This function encapsulates the full cross-sectional analysis so it can be
# called once per depth stratum without code duplication.

compute_cross_sectional <- function(data, cohort_label, gene_cols, gene_map) {

  # Filter to episodes that have a primary UTI status assignment
  with_status <- data %>% filter(!is.na(Infection_Status))

  # ------------------------------------------------------------------
  # C1: VF BURDEN BY STATUS
  # ------------------------------------------------------------------
  # For each primary UTI status (UTI, Not_UTI), compute descriptive
  # statistics for total VF gene count.  These are the core numbers
  # reported in the abstract (e.g., "median VF burden was 80 in Not_UTI
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
  # present, separately for UTI and Not_UTI.
  # The delta (UTI% - Not_UTI%) highlights genes that are differentially
  # prevalent in the new primary clinical contrast.
  uti_not_uti <- with_status %>% filter(Infection_Status %in% c("Not_UTI", "UTI"))
  n_not_uti <- sum(uti_not_uti$Infection_Status == "Not_UTI")
  n_uti <- sum(uti_not_uti$Infection_Status == "UTI")

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
      Not_UTI_n = get_val("Not_UTI", "n_present"), Not_UTI_N = get_val("Not_UTI", "n_total"),
      Not_UTI_pct = get_val("Not_UTI", "pct"),
      UTI_n = get_val("UTI", "n_present"), UTI_N = get_val("UTI", "n_total"),
      UTI_pct = get_val("UTI", "pct")
    )
  }) %>%
    bind_rows() %>%
    mutate(delta_UTI_minus_Not_UTI = round(UTI_pct - Not_UTI_pct, 1)) %>%
    left_join(gene_map %>% select(Gene, Category), by = c("gene" = "Gene")) %>%
    mutate(Category = coalesce(Category, "Unassigned"),
           cohort = cohort_label) %>%
    arrange(desc(abs(delta_UTI_minus_Not_UTI)))

  # ------------------------------------------------------------------
  # C3: EXPLORATORY FISHER EXACT TESTS (UTI vs Not_UTI)
  # ------------------------------------------------------------------
  # For each gene, test whether its prevalence differs significantly
  # between UTI and Not_UTI episodes using Fisher's exact test.
  #
  # ⚠ These tests assume independent observations.  Since participants
  # contribute multiple episodes, this assumption is VIOLATED.  These
  # results are screening tools only.  Formal inference uses GLMM in
  # script 14_genotype_phenotype_model.R.
  enrichment <- NULL
  if (n_uti >= 2) {
    enrichment <- lapply(gene_cols, function(g) {
      tab <- table(gene_present = uti_not_uti[[g]] > 0,
                   is_uti       = uti_not_uti$Infection_Status == "UTI")
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
# (e.g., "Iron acquisition" as a group) are enriched in UTI vs Not_UTI.
# This increases statistical power by grouping functionally related genes,
# and is more biologically interpretable.

cat_cols_ready <- grep("^cat_", names(vf_ready), value = TRUE)
uti_not_uti_full   <- vf_ready %>%
  filter(!is.na(Infection_Status), Infection_Status %in% c("Not_UTI", "UTI"))

cat_enrichment <- NULL
if (length(cat_cols_ready) > 0 && sum(uti_not_uti_full$Infection_Status == "UTI") >= 2) {
  cat_enrichment <- lapply(cat_cols_ready, function(cc) {
    # Binarise: does this episode have ≥1 gene from this category?
    uti_not_uti_full$cat_present <- uti_not_uti_full[[cc]] > 0
    tab <- table(uti_not_uti_full$cat_present, uti_not_uti_full$Infection_Status == "UTI")
    if (nrow(tab) < 2 || ncol(tab) < 2) return(NULL)
    ft <- tryCatch(fisher.test(tab), error = function(e) NULL)
    if (is.null(ft)) return(NULL)

    # Also run a Wilcoxon rank-sum test on the raw counts (non-binarised)
    # to test whether the NUMBER of genes in this category differs by status.
    wt <- tryCatch(
      wilcox.test(as.formula(paste0("`", cc, "` ~ Infection_Status")),
                  data = uti_not_uti_full),
      error = function(e) NULL
    )

    tibble(
      category = cc,
      fisher_OR    = round(ft$estimate, 3),
      fisher_p     = round(ft$p.value, 4),
      wilcox_p     = if (!is.null(wt)) round(wt$p.value, 4) else NA_real_,
      median_Not_UTI = median(uti_not_uti_full[[cc]][uti_not_uti_full$Infection_Status == "Not_UTI"]),
      median_UTI   = median(uti_not_uti_full[[cc]][uti_not_uti_full$Infection_Status == "UTI"])
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
# Descriptive summary: for each primary UTI status, what is the median and mean
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
  gene_tests <- results[["all"]]$enrichment %>%
    left_join(
      results[["all"]]$gene_prev %>%
        select(gene, Not_UTI_n, Not_UTI_N, UTI_n, UTI_N),
      by = "gene"
    ) %>%
    transmute(
      gene,
      Category,
      test = "Fisher exact test, UTI vs Not_UTI",
      OR,
      CI_lower,
      CI_upper,
      p_value,
      q_value_BH = p_adj_BH,
      Not_UTI_n,
      Not_UTI_N,
      UTI_n,
      UTI_N,
      exploratory_or_confirmatory = "Exploratory",
      interpretation_limitations = "Repeated resident isolates are not accounted for; rare genes produce unstable odds ratios; ST/lineage and event-driven sampling may confound UTI-vs-Not_UTI contrasts. Not_UTI is heterogeneous.",
      input_data = ready_file
    )
  write_csv(gene_tests, file.path(DIR_VF, "vf_gene_prevalence_tests.csv"))
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

# =============================================================================
# VF visualisation module 01: Burden
# =============================================================================

plot_data <- vf_ready %>%
  filter(!is.na(Infection_Status), Infection_Status %in% c("UTI", "Not_UTI")) %>%
  mutate(Infection_Status = status_for_plot(Infection_Status))

if (nrow(plot_data) > 0) {
  msg("Plotting VF burden by status: %s", status_count_text(plot_data))
  y_annot <- max(plot_data$vf_count_total, na.rm = TRUE) + 4
  status_n <- plot_data %>%
    count(Infection_Status, name = "n") %>%
    mutate(y = y_annot, label = paste0("n=", n))

  p_burden <- ggplot(plot_data, aes(x = Infection_Status, y = vf_count_total,
                                    fill = Infection_Status)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.65, width = 0.55) +
    geom_jitter(aes(colour = Infection_Status), width = 0.15, alpha = 0.35, size = 1.4) +
    geom_text(data = status_n, aes(x = Infection_Status, y = y, label = label),
              inherit.aes = FALSE, size = 3.2, fontface = "italic") +
    scale_fill_uti_status() +
    scale_colour_uti_status() +
    scale_y_continuous(expand = expansion(mult = c(0.04, 0.12))) +
    labs(
      title = "Distribution of E. coli virulence factor burden by primary UTI status",
      subtitle = "VF burden = number of detected virulence factor genes per VF/WGS-linked isolate",
      x = "Primary UTI status",
      y = "Detected VF genes per isolate",
      caption = vf_caption(
        ready_file, plot_data, "isolate-level",
        p_value_note = "No inferential p-value is shown on this plot; isolate-level UTI-vs-Not_UTI tests in this script are exploratory.",
        extra_note = "Not_UTI subgroups are retained in the dataset for sensitivity interpretation."
      )
    ) +
    plot_theme_vf(base_size = 12) +
    theme(legend.position = "none")

  ggsave(file.path(DIR_PLOTS_VF, "vf_burden_by_status.png"), p_burden,
         width = 7, height = 5.8, dpi = 300)
  ggsave(file.path(DIR_PLOTS_VF, "vf_burden_boxplot.png"), p_burden,
         width = 7, height = 5.8, dpi = 300)

  participant_burden <- plot_data %>%
    group_by(Participant_id, Infection_Status) %>%
    summarise(mean_vf_burden = mean(vf_count_total, na.rm = TRUE),
              n_isolates = n(), .groups = "drop")

  p_participant <- ggplot(participant_burden,
                          aes(x = Infection_Status, y = mean_vf_burden,
                              fill = Infection_Status)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.55, width = 0.5) +
    geom_jitter(aes(size = n_isolates, colour = Infection_Status),
                width = 0.15, alpha = 0.45) +
    scale_fill_uti_status() +
    scale_colour_uti_status() +
    scale_size_continuous(range = c(1.2, 3.5), breaks = c(1, 2, 4, 6, 8),
                          name = "Isolates in\nparticipant-status mean") +
    labs(
      title = "Participant-level summary of E. coli VF burden by primary UTI status",
      subtitle = "Each point is one participant-status mean, reducing repeated-isolate weighting",
      x = "Primary UTI status",
      y = "Mean detected VF genes per participant-status",
      caption = vf_caption(
        ready_file, participant_burden, "participant-status summary",
        p_value_note = "No p-value is shown; this is a descriptive companion to isolate-level burden plots.",
        extra_note = "Lines are not shown because this plot summarises all available participant-status combinations."
      )
    ) +
    plot_theme_vf(base_size = 11)

  ggsave(file.path(DIR_PLOTS_VF, "vf_burden_participant_summary.png"),
         p_participant, width = 7, height = 5.6, dpi = 300)

  paired_burden <- participant_burden %>%
    filter(Infection_Status %in% c("Not_UTI", "UTI")) %>%
    group_by(Participant_id) %>%
    filter(n_distinct(Infection_Status) == 2) %>%
    ungroup()

  if (n_distinct(paired_burden$Participant_id) >= 2) {
    msg("Plotting paired UTI-Not_UTI burden for %d participants.",
        n_distinct(paired_burden$Participant_id))
    p_paired <- ggplot(paired_burden,
                       aes(x = Infection_Status, y = mean_vf_burden,
                           group = Participant_id)) +
      geom_line(alpha = 0.35, colour = "grey40", linewidth = 0.4) +
      geom_point(aes(colour = Infection_Status), size = 2.2, alpha = 0.8) +
      scale_colour_uti_status() +
      labs(
        title = "Paired participant-level VF burden in residents with UTI and Not_UTI isolates",
        subtitle = sprintf("Each line connects Not_UTI and UTI participant-status means for n=%d residents",
                           n_distinct(paired_burden$Participant_id)),
        x = "Primary UTI status",
        y = "Mean detected VF genes per participant-status",
        caption = vf_caption(
          ready_file, paired_burden, "paired participant-level UTI-Not_UTI comparison",
          p_value_note = "No p-value is shown; this companion plot is descriptive and does not model within-resident clustering.",
          extra_note = "The UTI denominator is restricted to residents who also contributed Not_UTI isolates."
        )
      ) +
      plot_theme_vf(base_size = 11) +
      theme(legend.position = "none")

    ggsave(file.path(DIR_PLOTS_VF, "vf_burden_paired_uti_not_uti.png"),
           p_paired, width = 6, height = 5.2, dpi = 300)
  } else {
    msg("Skipping paired UTI-Not_UTI burden plot: fewer than two paired participants.")
  }
}

# =============================================================================
# VF visualisation module 02: Gene prevalence
# =============================================================================

gene_prev_all <- results[["all"]]$gene_prev
if (nrow(gene_prev_all) > 0) {
  gene_prev_plot <- gene_prev_all %>%
    mutate(
      overall_n = coalesce(Not_UTI_n, 0) + coalesce(UTI_n, 0),
      overall_N = coalesce(Not_UTI_N, 0) + coalesce(UTI_N, 0),
      overall_pct = if_else(overall_N > 0, round(100 * overall_n / overall_N, 1), NA_real_),
      Category = coalesce(Category, "Unassigned")
    )

  top_prev <- gene_prev_plot %>%
    slice_max(overall_pct, n = 40, with_ties = FALSE) %>%
    mutate(gene = factor(gene, levels = rev(gene)))

  p_top_gene <- ggplot(top_prev, aes(x = overall_pct, y = gene, fill = Category)) +
    geom_col(width = 0.7, alpha = 0.85) +
    geom_text(aes(label = paste0(overall_n, "/", overall_N)), hjust = -0.1, size = 2.4) +
    scale_x_continuous(labels = function(x) paste0(x, "%"),
                       expand = expansion(mult = c(0, 0.18))) +
    labs(
      title = "Most prevalent virulence factor genes among VF/WGS-linked E. coli isolates",
      subtitle = "Selection criterion: top 40 genes by prevalence across UTI and Not_UTI VF-ready isolates",
      x = "Isolates with gene detected",
      y = NULL,
      fill = "VF category",
      caption = vf_caption(
        ready_file, plot_data, "isolate-level gene prevalence",
        p_value_note = "No p-values are shown in this prevalence ranking.",
        extra_note = "Gene presence/absence is derived from the canonical VF matrix."
      )
    ) +
    plot_theme_vf(base_size = 10)

  ggsave(file.path(DIR_PLOTS_VF, "vf_top_gene_prevalence.png"),
         p_top_gene, width = 8, height = 9, dpi = 300)

  diff_top <- gene_prev_plot %>%
    filter(!is.na(delta_UTI_minus_Not_UTI), !is.na(Not_UTI_N), !is.na(UTI_N), UTI_N > 0, Not_UTI_N > 0) %>%
    slice_max(abs(delta_UTI_minus_Not_UTI), n = 30, with_ties = FALSE) %>%
    mutate(
      gene = factor(gene, levels = gene[order(delta_UTI_minus_Not_UTI)]),
      direction = if_else(delta_UTI_minus_Not_UTI >= 0, "Higher in UTI", "Higher in Not_UTI")
    )

  if (nrow(diff_top) > 0) {
    p_diff <- ggplot(diff_top, aes(x = delta_UTI_minus_Not_UTI, y = gene, fill = direction)) +
      geom_col(width = 0.65, alpha = 0.85) +
      geom_vline(xintercept = 0, linewidth = 0.4, colour = "grey45") +
      geom_text(
        aes(label = sprintf("UTI %d/%d; Not_UTI %d/%d", UTI_n, UTI_N, Not_UTI_n, Not_UTI_N),
            x = ifelse(delta_UTI_minus_Not_UTI >= 0,
                       delta_UTI_minus_Not_UTI + 1.5,
                       delta_UTI_minus_Not_UTI - 1.5)),
        hjust = ifelse(diff_top$delta_UTI_minus_Not_UTI >= 0, 0, 1),
        size = 2.25
      ) +
      scale_fill_manual(values = c("Higher in UTI" = "#D55E00", "Higher in Not_UTI" = "#0072B2")) +
      scale_x_continuous(labels = function(x) paste0(x, " pp"),
                         expand = expansion(mult = c(0.2, 0.2))) +
      labs(
        title = "Virulence factor genes with largest UTI-Not_UTI prevalence differences",
        subtitle = "Selection criterion: top 30 genes by absolute UTI minus Not_UTI prevalence difference",
        x = "Prevalence difference (UTI minus Not_UTI, percentage points)",
        y = NULL,
        fill = NULL,
        caption = vf_caption(
          ready_file, uti_not_uti_full, "isolate-level UTI-Not_UTI gene prevalence difference",
          p_value_note = "Differences are descriptive; gene-level Fisher tests are exploratory because residents may contribute repeated isolates.",
          multiple_testing_note = "Multiple testing is addressed in results/vf/vf_gene_prevalence_tests.csv using BH/FDR q-values.",
          extra_note = "Rare genes can produce unstable contrasts and should be interpreted with gene counts."
        )
      ) +
      plot_theme_vf(base_size = 10)

    ggsave(file.path(DIR_PLOTS_VF, "vf_gene_prevalence_difference_uti_not_uti.png"),
           p_diff, width = 9, height = 8, dpi = 300)
  }

  heat_genes <- gene_prev_plot %>%
    slice_max(overall_pct, n = 50, with_ties = FALSE) %>%
    pull(gene)
  heat_df <- gene_prev_plot %>%
    filter(gene %in% heat_genes) %>%
    select(gene, Category, Not_UTI_pct, UTI_pct) %>%
    pivot_longer(cols = c(Not_UTI_pct, UTI_pct),
                 names_to = "Infection_Status", values_to = "prevalence_pct") %>%
    mutate(
      Infection_Status = recode(Infection_Status,
                                "Not_UTI_pct" = "Not_UTI", "UTI_pct" = "UTI"),
      Infection_Status = factor(Infection_Status, levels = c("UTI", "Not_UTI")),
      gene = factor(gene, levels = rev(heat_genes))
    )

  p_heat <- ggplot(heat_df, aes(x = Infection_Status, y = gene, fill = prevalence_pct)) +
    geom_tile(colour = "white", linewidth = 0.2) +
    geom_text(aes(label = ifelse(is.na(prevalence_pct), "", sprintf("%.0f", prevalence_pct))),
              size = 2.1) +
    scale_fill_gradient(low = "white", high = "#0072B2", na.value = "grey90",
                        labels = function(x) paste0(x, "%")) +
    labs(
      title = "VF gene prevalence heatmap across primary UTI status",
      subtitle = "Selection criterion: top 50 genes by overall prevalence in VF-ready isolates",
      x = "Primary UTI status",
      y = NULL,
      fill = "Prevalence",
      caption = vf_caption(
        ready_file, plot_data, "isolate-level gene prevalence heatmap",
        p_value_note = "No p-values are shown; heatmap is descriptive.",
        extra_note = "Percentages are calculated within each primary UTI status denominator."
      )
    ) +
    plot_theme_vf(base_size = 9)

  ggsave(file.path(DIR_PLOTS_VF, "vf_gene_prevalence_heatmap.png"),
         p_heat, width = 6.5, height = 11, dpi = 300)
}

# =============================================================================
# VF visualisation module 03: Category burden companion plot
# =============================================================================

if (length(cat_cols_ready) > 0) {
  cat_long <- vf_ready %>%
    filter(!is.na(Infection_Status), Infection_Status %in% c("UTI", "Not_UTI")) %>%
    mutate(Infection_Status = status_for_plot(Infection_Status)) %>%
    select(Participant_id, Infection_Status, all_of(cat_cols_ready)) %>%
    pivot_longer(cols = all_of(cat_cols_ready),
                 names_to = "Category", values_to = "count") %>%
    mutate(Category = format_category_label(Category))

  cat_summary <- cat_long %>%
    group_by(Category, Infection_Status) %>%
    summarise(median_count = median(count, na.rm = TRUE),
              mean_count = mean(count, na.rm = TRUE),
              .groups = "drop")

  p_cat <- ggplot(cat_long,
                  aes(x = Infection_Status, y = count, fill = Infection_Status)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.55, width = 0.55) +
    geom_jitter(aes(colour = Infection_Status), width = 0.12, alpha = 0.25, size = 0.8) +
    facet_wrap(~Category, scales = "free_y", ncol = 3) +
    scale_fill_uti_status() +
    scale_colour_uti_status() +
    labs(
      title = "Virulence factor category burden across primary UTI status",
      subtitle = "Category burden = number of detected genes per isolate within each curated VF category",
      x = "Primary UTI status",
      y = "Detected genes in category",
      caption = vf_caption(
        ready_file, vf_ready, "isolate-level category burden",
        p_value_note = "Category-level Fisher/Wilcoxon tests in result tables are exploratory and not adjusted for repeated resident samples.",
        multiple_testing_note = "Category-level p-values are BH-adjusted in results/vf/vf_category_enrichment.csv.",
        extra_note = "Categories are descriptive biological groupings and should not be interpreted as validated causal virulence scores."
      )
    ) +
    plot_theme_vf(base_size = 10) +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 35, hjust = 1))

  ggsave(file.path(DIR_PLOTS_VF, "vf_category_burden_by_status.png"), p_cat,
         width = 9, height = 6.5, dpi = 300)
  ggsave(file.path(DIR_PLOTS_VF, "vf_category_barplot.png"), p_cat,
         width = 9, height = 6.5, dpi = 300)

  write_csv(cat_summary, file.path(DIR_VF, "vf_category_burden_plot_summary.csv"))
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
  sl("EXPLORATORY Fisher exact (UTI vs Not_UTI, full cohort):")
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

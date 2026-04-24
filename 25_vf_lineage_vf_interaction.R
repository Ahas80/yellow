#!/usr/bin/env Rscript
# ==============================================================================
# 25_vf_lineage_vf_interaction.R
# ==============================================================================
#
# GOAL:
#   Determine whether observed VF burden differences between ASB and UTI
#   are driven by clinical status or by the underlying bacterial lineage (ST).
#
# WHY THIS SCRIPT IS CRITICAL:
#   This is the study's most important confounding check.  Consider this
#   scenario:
#     - ST131 naturally carries 90 VF genes.
#     - ST73 naturally carries 70 VF genes.
#     - If ST131 happens to cause more UTIs in our cohort, then a naive
#       comparison would show "UTI episodes have more VFs" — but this is
#       entirely a lineage artefact, not a genuine VF→UTI association.
#
#   Without this analysis, any VF-status association claim is suspect.
#
#   The script addresses three sequential questions:
#     Q1: Does VF burden vary across STs?
#         (Kruskal-Wallis test)
#     Q2: Within each major ST, is there a VF burden difference between
#         ASB and UTI?
#         (Within-ST Wilcoxon tests — small sample sizes expected)
#     Q3: Is ST composition different between ASB and UTI?
#         (Fisher exact on ST×Status contingency table)
#
#   Interpretation:
#     Q1=YES + Q3=YES → Lineage IS a confounder.  Add ST as a covariate
#                        in 14_genotype_phenotype_model.R.
#     Q1=YES + Q3=NO  → STs carry different VFs but aren't differentially
#                        distributed.  Less confounding concern.
#     Q1=NO           → Lineage is unlikely to confound VF-status results.
#
# RELATIONSHIP TO EXISTING SCRIPTS:
#   - 17_lineage_analysis.R computes UTI *risk* per ST (what proportion of
#     each ST's episodes are UTI?).  It does NOT examine VF burden per ST.
#     This script fills that gap.
#   - 14_genotype_phenotype_model.R runs GLMM on individual VF genes.
#     If this script detects confounding, the recommendation is to add ST
#     as a covariate there, NOT to build a separate modelling script.
#
# INPUTS:
#   - results/vf/vf_analysis_ready.csv   (from 22_vf_build_analysis_dataset.R)
#   - results/vf/gene_map.csv            (for category-level medians per ST)
#
# OUTPUTS (all in results/vf/):
#   - vf_burden_by_st.csv                  VF burden per ST (top STs only)
#   - vf_burden_by_st_and_status.csv       VF burden per ST × Status
#   - vf_lineage_confounding_summary.txt   Human-readable confounding report
#
# PLOTS (in plots/vf/):
#   - vf_burden_by_st.png                  Boxplot of VF burden per top STs
#   - vf_burden_st_x_status.png            Faceted: ASB vs UTI within each ST
# ==============================================================================

source("00_config.R")
source("R/plot_helpers.R")  # Canonical infection status colours
suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
})

msg("Starting 25_vf_lineage_vf_interaction.R")

# ==============================================================================
# 1. LOAD DATA
# ==============================================================================

ready_file <- file.path(DIR_VF, "vf_analysis_ready.csv")
if (!file.exists(ready_file)) stop("Missing ", ready_file, ". Run 22_vf_build_analysis_dataset.R first.")
vf_ready <- read_csv(ready_file, show_col_types = FALSE) %>%
  mutate(Participant_id = as.character(Participant_id),
         ST = as.character(ST))

gene_map <- read_csv(file.path(DIR_VF, "gene_map.csv"), show_col_types = FALSE) %>%
  mutate(Gene = as.character(Gene),
         Category = coalesce(as.character(Category), "Unassigned"))

cat_cols <- grep("^cat_", names(vf_ready), value = TRUE)

msg("Loaded: %d rows, %d with ST, %d distinct STs",
    nrow(vf_ready), sum(!is.na(vf_ready$ST)), n_distinct(vf_ready$ST, na.rm = TRUE))

# ==============================================================================
# 2. CHECK ST DATA AVAILABILITY
# ==============================================================================
# If MLST data was not joined in 22_ (e.g., mlst_with_meta.csv was missing),
# this script cannot run.  Fail gracefully with a clear message.

has_st <- vf_ready %>% filter(!is.na(ST), !is.na(Infection_Status))

if (nrow(has_st) < 10) {
  msg("WARNING: Only %d episodes have both ST and Infection_Status.", nrow(has_st))
  msg("Lineage confounding analysis will be limited.")
}

if (nrow(has_st) == 0) {
  msg("ERROR: No episodes with ST data.  Cannot run lineage analysis.")
  msg("Ensure 22_vf_build_analysis_dataset.R successfully joined MLST data.")
  msg("Writing empty outputs and exiting.")

  writeLines("No ST data available for lineage–VF analysis.",
             file.path(DIR_VF, "vf_lineage_confounding_summary.txt"))
  quit(save = "no")
}

# ==============================================================================
# 3. VF BURDEN BY SEQUENCE TYPE
# ==============================================================================
# We focus on STs with sufficient data (≥5 episodes) to make meaningful
# comparisons.  With <5 episodes, a single outlier dominates the summary.

st_counts <- has_st %>%
  count(ST, name = "n_episodes") %>%
  arrange(desc(n_episodes))

min_episodes <- 5
top_sts <- st_counts %>% filter(n_episodes >= min_episodes) %>% pull(ST)

msg("STs with ≥%d episodes: %d (out of %d total)",
    min_episodes, length(top_sts), nrow(st_counts))

# If no ST has ≥5 episodes, relax the threshold to ≥3 to still produce output
if (length(top_sts) == 0) {
  msg("WARNING: No ST has ≥%d episodes.  Relaxing threshold to ≥3.", min_episodes)
  min_episodes <- 3
  top_sts <- st_counts %>% filter(n_episodes >= min_episodes) %>% pull(ST)
}

# Compute VF burden descriptive statistics for each top ST.
# This table directly answers Q1: do different STs carry different VF loads?
burden_by_st <- has_st %>%
  filter(ST %in% top_sts) %>%
  group_by(ST) %>%
  summarise(
    n_episodes     = n(),
    n_participants = n_distinct(Participant_id),
    mean_vf   = round(mean(vf_count_total), 1),
    sd_vf     = round(sd(vf_count_total), 1),
    median_vf = median(vf_count_total),
    q25_vf    = quantile(vf_count_total, 0.25),
    q75_vf    = quantile(vf_count_total, 0.75),
    .groups   = "drop"
  ) %>%
  arrange(desc(median_vf))

# Also add category-level medians per ST (e.g., how many Iron Acquisition
# genes does ST131 typically carry vs ST73?)
for (cc in cat_cols) {
  cat_med <- has_st %>%
    filter(ST %in% top_sts) %>%
    group_by(ST) %>%
    summarise(!!paste0(cc, "_median") := median(.data[[cc]], na.rm = TRUE),
              .groups = "drop")
  burden_by_st <- burden_by_st %>% left_join(cat_med, by = "ST")
}

write_csv(burden_by_st, file.path(DIR_VF, "vf_burden_by_st.csv"))

# Compute VF burden per ST × Status (ASB/UTI only).
# This answers Q2: within the same ST, does VF burden differ by status?
asb_uti_st <- has_st %>%
  filter(ST %in% top_sts, Infection_Status %in% c("ASB", "UTI"))

burden_st_status <- asb_uti_st %>%
  group_by(ST, Infection_Status) %>%
  summarise(
    n_episodes = n(),
    mean_vf    = round(mean(vf_count_total), 1),
    median_vf  = median(vf_count_total),
    .groups    = "drop"
  ) %>%
  arrange(ST, Infection_Status)

write_csv(burden_st_status, file.path(DIR_VF, "vf_burden_by_st_and_status.csv"))

# ==============================================================================
# 4. STATISTICAL TESTS FOR CONFOUNDING
# ==============================================================================
# We apply the three-question framework described in the header.

summary_lines <- character()
sl <- function(...) summary_lines <<- c(summary_lines, sprintf(...))

sl("=== VF–LINEAGE CONFOUNDING ASSESSMENT ===")
sl("Generated: %s", format(Sys.time()))
sl("")

# ------------------------------------------------------------------
# Q1: Does VF burden differ significantly across STs?
# ------------------------------------------------------------------
# Kruskal-Wallis is a non-parametric test for differences in medians
# across >2 groups (distribution-free alternative to one-way ANOVA).
# If significant, it means STs carry genuinely different VF loads,
# which is a prerequisite for lineage confounding.
if (length(top_sts) >= 2) {
  kw_data <- has_st %>% filter(ST %in% top_sts)
  kw <- kruskal.test(vf_count_total ~ factor(ST), data = kw_data)
  sl("Q1: Does VF burden differ across STs?")
  sl("  Kruskal–Wallis chi-sq = %.2f, df = %d, p = %.4f",
     kw$statistic, kw$parameter, kw$p.value)
  if (kw$p.value < 0.05) {
    sl("  → YES: Significant VF burden variation across STs.")
    sl("  → This means lineage is a potential confounder.")
  } else {
    sl("  → NO: VF burden does not significantly vary across STs.")
  }
  sl("")
}

# ------------------------------------------------------------------
# Q2: Within each major ST, does VF burden differ between ASB and UTI?
# ------------------------------------------------------------------
# This is the key within-lineage test.  If VF burden is the SAME for
# ASB and UTI episodes of the same ST, then VF differences observed
# in the overall cohort are likely driven by ST composition, not by
# VF content per se.
#
# NOTE: Small sample sizes within STs are expected.  Many STs will have
# too few UTI episodes for a valid test.  This is a known limitation.
sl("Q2: Within each top ST, does VF burden differ between ASB and UTI?")
for (st in top_sts) {
  st_data <- asb_uti_st %>% filter(ST == st)
  n_asb <- sum(st_data$Infection_Status == "ASB")
  n_uti <- sum(st_data$Infection_Status == "UTI")

  if (n_asb >= 2 && n_uti >= 2) {
    # Wilcoxon rank-sum (Mann-Whitney U) test: non-parametric comparison
    # of two independent groups.
    wt <- tryCatch(
      wilcox.test(vf_count_total ~ Infection_Status, data = st_data),
      error = function(e) NULL
    )
    if (!is.null(wt)) {
      sl("  ST%s (ASB=%d, UTI=%d): Wilcoxon p = %.4f, median ASB=%.0f, UTI=%.0f",
         st, n_asb, n_uti, wt$p.value,
         median(st_data$vf_count_total[st_data$Infection_Status == "ASB"]),
         median(st_data$vf_count_total[st_data$Infection_Status == "UTI"]))
    }
  } else {
    sl("  ST%s (ASB=%d, UTI=%d): Too few in one group for within-ST test.", st, n_asb, n_uti)
  }
}

# ------------------------------------------------------------------
# Q3: Does ST composition differ between ASB and UTI?
# ------------------------------------------------------------------
# This tests whether certain STs are over-represented in UTI vs ASB.
# If yes, and Q1 is also yes, then lineage is a LIKELY CONFOUNDER:
# UTI episodes might have higher VF counts simply because they tend
# to harbour VF-heavy STs, not because VFs cause UTI.
sl("")
sl("Q3: Does ST composition differ between ASB and UTI?")
st_status_tab <- asb_uti_st %>%
  filter(ST %in% top_sts) %>%
  count(ST, Infection_Status) %>%
  pivot_wider(names_from = Infection_Status, values_from = n, values_fill = 0)

if (nrow(st_status_tab) >= 2 && ncol(st_status_tab) >= 3) {
  mat <- as.matrix(st_status_tab %>% select(-ST))
  # Use simulated p-values because the contingency table may be sparse
  ft <- tryCatch(fisher.test(mat, simulate.p.value = TRUE, B = 10000),
                 error = function(e) NULL)
  if (!is.null(ft)) {
    sl("  Fisher exact (simulated): p = %.4f", ft$p.value)
    if (ft$p.value < 0.05) {
      sl("  → YES: ST composition differs between ASB and UTI.")
      sl("  → Combined with Q1, this means lineage is a LIKELY CONFOUNDER.")
    } else {
      sl("  → NO: ST composition is not significantly different.")
    }
  }
}

# ------------------------------------------------------------------
# Interpretation guide
# ------------------------------------------------------------------
sl("")
sl("INTERPRETATION GUIDANCE:")
sl("  If Q1=YES and Q3=YES → VF–status differences may be ST-driven artefacts.")
sl("  → Add ST as covariate in 14_genotype_phenotype_model.R.")
sl("  If Q1=YES and Q3=NO  → STs carry different VFs but aren't over-represented")
sl("    in either status.  Less confounding concern.")
sl("  If Q1=NO  → Lineage is unlikely to confound VF–status results.")

# ==============================================================================
# 5. PLOTS
# ==============================================================================

ensure_dir(DIR_PLOTS_VF)

# PLOT 1: VF burden boxplot by ST
#   Shows the range of VF gene counts for each major ST.
#   Large differences here (Q1=YES) mean that lineage choice alone
#   substantially determines VF load.
if (length(top_sts) >= 2) {
  plot_st <- has_st %>%
    filter(ST %in% top_sts) %>%
    mutate(ST_label = paste0("ST", ST))

  p_st <- ggplot(plot_st, aes(x = reorder(ST_label, -vf_count_total, FUN = median),
                                y = vf_count_total)) +
    geom_boxplot(outlier.shape = NA, width = 0.5, fill = "grey90") +
    geom_jitter(width = 0.15, alpha = 0.5, size = 1.5) +
    labs(title = "VF Gene Burden by Sequence Type",
         subtitle = sprintf("STs with ≥%d episodes", min_episodes),
         x = "Sequence Type", y = "Total VF Genes") +
    theme_minimal(base_size = 11) +
    theme(axis.text.x = element_text(angle = 30, hjust = 1))

  ggsave(file.path(DIR_PLOTS_VF, "vf_burden_by_st.png"), p_st,
         width = max(6, length(top_sts) * 0.8), height = 5, dpi = 300)
}

# PLOT 2: VF burden by ST × Status (ASB vs UTI, faceted)
#   This is the visual complement to the Q2 tests above.
#   For each ST that has both ASB and UTI episodes, show side-by-side
#   boxplots.  If boxes overlap heavily within each ST, the VF-status
#   association is likely an artefact of ST composition.
if (nrow(asb_uti_st) > 0 && length(top_sts) >= 2) {
  plot_st_status <- asb_uti_st %>%
    mutate(ST_label = paste0("ST", ST))

  # Only facet STs that have at least 1 episode of each status
  sts_both <- plot_st_status %>%
    count(ST_label, Infection_Status) %>%
    group_by(ST_label) %>%
    filter(n() >= 2) %>%
    pull(ST_label) %>% unique()

  if (length(sts_both) >= 1) {
    p_stx <- ggplot(plot_st_status %>% filter(ST_label %in% sts_both),
                    aes(x = Infection_Status, y = vf_count_total,
                        fill = Infection_Status)) +
      geom_boxplot(outlier.shape = NA, width = 0.5, alpha = 0.7) +
      geom_jitter(width = 0.15, alpha = 0.4, size = 1.2) +
      scale_fill_infection() +
      facet_wrap(~ST_label, scales = "free_x") +
      labs(title = "VF Burden: ASB vs UTI Within Each Sequence Type",
           subtitle = "STs with both ASB and UTI episodes",
           x = NULL, y = "Total VF Genes") +
      theme_minimal(base_size = 10) +
      theme(legend.position = "none",
            strip.text = element_text(face = "bold"))

    ggsave(file.path(DIR_PLOTS_VF, "vf_burden_st_x_status.png"), p_stx,
           width = max(6, length(sts_both) * 2.5),
           height = max(4, ceiling(length(sts_both) / 3) * 3),
           dpi = 300)
  }
}

# ==============================================================================
# 6. WRITE SUMMARY
# ==============================================================================

writeLines(summary_lines, file.path(DIR_VF, "vf_lineage_confounding_summary.txt"))
msg("Summary written to %s", file.path(DIR_VF, "vf_lineage_confounding_summary.txt"))

# Also print to console
cat(paste(summary_lines, collapse = "\n"), "\n")

msg("✓ 25_vf_lineage_vf_interaction.R complete.")

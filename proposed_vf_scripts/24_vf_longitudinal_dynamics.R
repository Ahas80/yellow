#!/usr/bin/env Rscript
# ==============================================================================
# 24_vf_longitudinal_dynamics.R
# ------------------------------------------------------------------------------
# Role: [VF Longitudinal] - Quantify VF stability and change across consecutive
# episodes within participants.
#
# Inputs:
#   - results/vf/vf_transition_dataset.csv
#   - results/vf/vf_episode_dataset.csv
#
# Outputs:
#   - results/vf/vf_transition_summary.csv
#   - results/vf/vf_zero_change_summary.csv
#   - results/vf/vf_asb_to_uti_candidates.csv
#   - plots/vf/vf_jaccard_by_transition_type.png
#
# Notes:
#   - This script is descriptive by design.
#   - It should stay separate from cross-sectional screening to avoid mixing
#     episode-level and transition-level inference.
# ==============================================================================

source("00_config.R")
source("R/plot_helpers.R")

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(tidyr)
  library(forcats)
})

msg <- function(...) message(sprintf(...))

in_trans <- file.path(DIR_RESULTS, "vf", "vf_transition_dataset.csv")
in_episode <- file.path(DIR_RESULTS, "vf", "vf_episode_dataset.csv")
if (!file.exists(in_trans)) stop("Run 22_vf_build_analysis_dataset.R first.")
if (!file.exists(in_episode)) stop("Run 22_vf_build_analysis_dataset.R first.")

tr <- read_csv(in_trans, show_col_types = FALSE)
episode <- read_csv(in_episode, show_col_types = FALSE)

# 1. Transition summaries
summ <- tr %>%
  group_by(transition_type) %>%
  summarise(
    n = n(),
    median_jaccard = median(jaccard_similarity, na.rm = TRUE),
    mean_jaccard = mean(jaccard_similarity, na.rm = TRUE),
    zero_change_n = sum(zero_vf_change, na.rm = TRUE),
    zero_change_pct = 100 * mean(zero_vf_change, na.rm = TRUE),
    mean_gained = mean(gained_genes, na.rm = TRUE),
    mean_lost = mean(lost_genes, na.rm = TRUE),
    median_gained = median(gained_genes, na.rm = TRUE),
    median_lost = median(lost_genes, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(n))
write_csv(summ, file.path(DIR_RESULTS, "vf", "vf_transition_summary.csv"))

zero_tbl <- tr %>%
  summarise(
    transitions = n(),
    zero_change_n = sum(zero_vf_change, na.rm = TRUE),
    zero_change_pct = 100 * mean(zero_vf_change, na.rm = TRUE),
    median_jaccard = median(jaccard_similarity, na.rm = TRUE)
  )
write_csv(zero_tbl, file.path(DIR_RESULTS, "vf", "vf_zero_change_summary.csv"))

asb_to_uti <- tr %>%
  filter(Status_A == "ASB", Status_B == "UTI") %>%
  arrange(desc(gained_genes), desc(lost_genes), jaccard_similarity)
write_csv(asb_to_uti, file.path(DIR_RESULTS, "vf", "vf_asb_to_uti_candidates.csv"))

# 2. Plot
p <- ggplot(tr, aes(x = fct_reorder(transition_type, jaccard_similarity, .fun = median, na.rm = TRUE),
                    y = jaccard_similarity)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.75) +
  geom_jitter(width = 0.18, alpha = 0.45, size = 1.5) +
  coord_flip() +
  labs(
    title = "VF Jaccard similarity across consecutive transitions",
    x = "Transition type",
    y = "Jaccard similarity"
  ) +
  theme_minimal(base_size = 12)

ggsave(file.path(DIR_PLOTS, "vf", "vf_jaccard_by_transition_type.png"), p, width = 7.5, height = 5.5, dpi = 300)
msg("Done.")

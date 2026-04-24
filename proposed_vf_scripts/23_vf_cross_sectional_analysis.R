#!/usr/bin/env Rscript
# ==============================================================================
# 23_vf_cross_sectional_analysis.R
# ------------------------------------------------------------------------------
# Role: [VF Cross-sectional] - Compare VF burden and individual VF prevalence
# across clinical states, with clear separation between exploratory and more
# defensible models.
#
# Inputs:
#   - results/vf/vf_episode_dataset.csv
#   - results/vf/vf_gene_long_dataset.csv
#
# Outputs:
#   - results/vf/vf_burden_by_status.csv
#   - results/vf/vf_gene_prevalence_asb_vs_uti.csv
#   - results/vf/vf_category_burden_by_status.csv
#   - results/vf/vf_status_timepoint_confounding_table.csv
#   - plots/vf/vf_burden_by_status_boxplot.png
#
# Notes:
#   - Fisher tests are retained as exploratory screening.
#   - Effect sizes and BH correction are mandatory.
#   - Status-timepoint confounding is explicitly reported because UTI episodes
#     may be structurally concentrated at Uricult.
# ==============================================================================

source("00_config.R")
source("R/plot_helpers.R")

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
  library(ggplot2)
  library(purrr)
  library(forcats)
})

msg <- function(...) message(sprintf(...))

in_episode <- file.path(DIR_RESULTS, "vf", "vf_episode_dataset.csv")
in_gene    <- file.path(DIR_RESULTS, "vf", "vf_gene_long_dataset.csv")
if (!file.exists(in_episode)) stop("Run 22_vf_build_analysis_dataset.R first.")
if (!file.exists(in_gene)) stop("Run 22_vf_build_analysis_dataset.R first.")

episode <- read_csv(in_episode, show_col_types = FALSE)
gene_long <- read_csv(in_gene, show_col_types = FALSE)

dir.create(file.path(DIR_RESULTS, "vf"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(DIR_PLOTS, "vf"), recursive = TRUE, showWarnings = FALSE)

main <- episode %>% filter(Infection_Status %in% c("ASB", "UTI"))

# 1. Burden summary
burden_tbl <- main %>%
  group_by(Infection_Status) %>%
  summarise(
    n = n(),
    participants = n_distinct(Participant_id),
    median_vf = median(VF_Burden, na.rm = TRUE),
    iqr_low = quantile(VF_Burden, 0.25, na.rm = TRUE),
    iqr_high = quantile(VF_Burden, 0.75, na.rm = TRUE),
    mean_vf = mean(VF_Burden, na.rm = TRUE),
    sd_vf = sd(VF_Burden, na.rm = TRUE),
    min_vf = min(VF_Burden, na.rm = TRUE),
    max_vf = max(VF_Burden, na.rm = TRUE),
    .groups = "drop"
  )
write_csv(burden_tbl, file.path(DIR_RESULTS, "vf", "vf_burden_by_status.csv"))

# 2. Gene-level prevalence ASB vs UTI
fisher_gene <- gene_long %>%
  filter(Infection_Status %in% c("ASB", "UTI")) %>%
  group_by(gene) %>%
  group_modify(~ {
    d <- .x
    a <- sum(d$present[d$Infection_Status == "UTI"] == 1, na.rm = TRUE)
    b <- sum(d$present[d$Infection_Status == "UTI"] == 0, na.rm = TRUE)
    c <- sum(d$present[d$Infection_Status == "ASB"] == 1, na.rm = TRUE)
    e <- sum(d$present[d$Infection_Status == "ASB"] == 0, na.rm = TRUE)
    m <- matrix(c(a, b, c, e), nrow = 2, byrow = TRUE)
    ft <- fisher.test(m)
    tibble(
      uti_present = a,
      uti_absent = b,
      asb_present = c,
      asb_absent = e,
      prev_uti = a / (a + b),
      prev_asb = c / (c + e),
      prev_diff = (a / (a + b)) - (c / (c + e)),
      odds_ratio = unname(ft$estimate),
      p_value = ft$p.value
    )
  }) %>%
  ungroup() %>%
  mutate(p_adj_bh = p.adjust(p_value, method = "BH")) %>%
  arrange(p_value)
write_csv(fisher_gene, file.path(DIR_RESULTS, "vf", "vf_gene_prevalence_asb_vs_uti.csv"))

# 3. Category burden summaries if available
cat_cols <- grep("^cat_", names(main), value = TRUE)
if (length(cat_cols) > 0) {
  cat_tbl <- main %>%
    group_by(Infection_Status) %>%
    summarise(across(all_of(cat_cols), list(mean = ~ mean(.x, na.rm = TRUE), median = ~ median(.x, na.rm = TRUE))), .groups = "drop")
  write_csv(cat_tbl, file.path(DIR_RESULTS, "vf", "vf_category_burden_by_status.csv"))
}

# 4. Explicit status-timepoint confounding table
conf_tbl <- episode %>%
  count(Timepoint, Infection_Status, name = "n") %>%
  tidyr::pivot_wider(names_from = Infection_Status, values_from = n, values_fill = 0)
write_csv(conf_tbl, file.path(DIR_RESULTS, "vf", "vf_status_timepoint_confounding_table.csv"))

# 5. Plot
p <- ggplot(main, aes(x = Infection_Status, y = VF_Burden, fill = Infection_Status)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.8) +
  geom_jitter(width = 0.15, alpha = 0.5, size = 1.7) +
  labs(
    title = "Virulence-factor burden by clinical status",
    subtitle = "Exploratory cross-sectional comparison",
    x = NULL,
    y = "Distinct VF genes detected"
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "none")

ggsave(file.path(DIR_PLOTS, "vf", "vf_burden_by_status_boxplot.png"), p, width = 7, height = 5, dpi = 300)
msg("Done.")

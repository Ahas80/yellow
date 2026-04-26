#!/usr/bin/env Rscript
# ==============================================================================
# 17_lineage_analysis.R
# ==============================================================================
#
# GOAL:
#   Compute UTI risk per Sequence Type: for each ST with sufficient data,
#   what proportion of its episodes are UTI?  This identifies STs that are
#   disproportionately associated with symptomatic infection in this cohort.
#
# NOTE:
#   Script 25_vf_lineage_vf_interaction.R extends this by linking VF burden
#   to ST and testing for lineage confounding of VF–status associations.
#
# ------------------------------------------------------------------------------
# Role: [Analysis] - Priority 3: Lineage-Specific Virulence.
#
# Inputs:
#   - results/clinical/status_map.csv
#   - results/mlst/mlst_with_meta.csv (or mlst_all.tsv)
#   - results/models/gwas_multivariable_glmm.csv (to compare hits)
#
# Outputs:
#   - results/lineage/st_risk_profile.csv
#   - results/lineage/st_risk_plot.png
#
# Purpose:
#   - Determine if specific STs (e.g. ST131) are inherently more "symptomatic"
#     than others in this cohort.
#   - Risk stratification of clones.
#
# KEY DESIGN DECISIONS:
#   - Computes per-ST UTI proportions only for analysable ST strata to avoid
#     unstable rates from single-observation lineages.
#
# POSITION IN PIPELINE:
#   - Phase 2 comparative genomics lineage-risk summary, complements 14 and 25.
#
# NOTES / LIMITATIONS:
#   - ST-risk estimates are confounded by repeated measures and small per-ST n.
# ==============================================================================

source("00_config.R")
source("R/plot_helpers.R")
suppressPackageStartupMessages({
    library(dplyr)
    library(readr)
    library(ggplot2)
    library(tidyr)
    library(scales)
})

msg("Starting 17_lineage_analysis.R")

# 1. Load Data
# ------------------------------------------------------------------------------
# Status
status <- read_csv(file.path(DIR_CLINICAL, "status_map.csv"), show_col_types = FALSE) %>%
    filter(Infection_Status %in% c("UTI", "ASB"))

# MLST
# Try mlst_with_meta first, else build it
mlst_file <- file.path(DIR_MLST, "mlst_with_meta.csv")
if (file.exists(mlst_file)) {
    mlst <- read_csv(mlst_file, show_col_types = FALSE)
} else {
    # Fallback to mlst_all.tsv and join with metadata
    mlst_raw <- read_tsv(file.path(DIR_MLST, "mlst_all.tsv"), show_col_types = FALSE)
    meta <- read_csv("assembly_metadata.csv", show_col_types = FALSE)
    mlst <- mlst_raw %>%
        select(Isolate_ID, ST) %>%
        inner_join(meta, by = c("Isolate_ID" = "file_name")) # Adjust join if needed
}

# 2. Harmonize & Join
# ------------------------------------------------------------------------------
# We need to link ST to Episode (Participant + Timepoint)
# MLST is per Isolate. Status is per Episode.
# We'll take the ST of the isolate associated with the episode.

# Ensure columns match for join
mlst_clean <- mlst %>%
    select(Participant_id, Timepoint, ST) %>%
    distinct() %>%
    mutate(
        Participant_id = as.character(Participant_id),
        Timepoint = as.character(Timepoint)
    )

status_clean <- status %>%
    mutate(
        Participant_id = as.character(Participant_id),
        Timepoint = as.character(Timepoint)
    )

data_merged <- status_clean %>%
    inner_join(mlst_clean, by = c("Participant_id", "Timepoint"), relationship = "many-to-many")

msg("Linked %d episodes to ST data", nrow(data_merged))

# 3. Calculate Risk per ST
# ------------------------------------------------------------------------------
# Filter for common STs (e.g. n >= 5 episodes)
st_counts <- data_merged %>%
    count(ST) %>%
    arrange(desc(n))

top_sts <- st_counts %>%
    filter(n >= 5) %>%
    pull(ST)

msg("Analyzing %d STs with >= 5 episodes", length(top_sts))

st_risk <- data_merged %>%
    filter(ST %in% top_sts) %>%
    group_by(ST) %>%
    summarise(
        n_total = n(),
        n_UTI = sum(Infection_Status == "UTI"),
        n_ASB = sum(Infection_Status == "ASB"),
        Risk_UTI = n_UTI / n_total,
        .groups = "drop"
    ) %>%
    arrange(desc(Risk_UTI))

# 4. Statistical Test (Fisher's Exact per ST)
# ------------------------------------------------------------------------------
# [STAT] NOTE ON INDEPENDENCE:
# We compare each ST against "All Other STs" using Fisher's Exact test.
# Because the same participant can contribute multiple episodes (and multiple
# isolates of the same ST), this test violates the assumption of independent
# observations. It should be interpreted as an exploratory screen for
# over-represented STs, not as formal inference.
calc_p <- function(st_target, df) {
    # Contingency Table
    #       UTI  ASB
    # Target  a    b
    # Other   c    d

    a <- sum(df$ST == st_target & df$Infection_Status == "UTI")
    b <- sum(df$ST == st_target & df$Infection_Status == "ASB")
    c <- sum(df$ST != st_target & df$Infection_Status == "UTI")
    d <- sum(df$ST != st_target & df$Infection_Status == "ASB")

    mat <- matrix(c(a, c, b, d), nrow = 2)
    test <- fisher.test(mat)

    tibble(
        ST = st_target,
        OR = test$estimate,
        p_value = test$p.value,
        CI_low = test$conf.int[1],
        CI_high = test$conf.int[2]
    )
}

stats_res <- purrr::map_dfr(top_sts, calc_p, df = data_merged)

final_res <- st_risk %>%
    left_join(stats_res, by = "ST") %>%
    mutate(FDR = p.adjust(p_value, method = "BH")) %>%
    arrange(p_value)

print(final_res)

# 5. Save & Plot
# ------------------------------------------------------------------------------
out_dir <- file.path(DIR_RESULTS, "lineage")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

write_csv(final_res, file.path(out_dir, "st_risk_profile.csv"))

# Plot
p <- ggplot(final_res, aes(x = reorder(ST, Risk_UTI), y = Risk_UTI, fill = Risk_UTI)) +
    geom_col(color = "black", width = 0.7) +
    geom_errorbar(
        aes(
            ymin = pmax(0, Risk_UTI - 1.96 * sqrt(Risk_UTI * (1 - Risk_UTI) / n_total)),
            ymax = pmin(1, Risk_UTI + 1.96 * sqrt(Risk_UTI * (1 - Risk_UTI) / n_total))
        ),
        width = 0.2, alpha = 0.5
    ) +
    geom_text(aes(label = sprintf("n=%d", n_total)), vjust = -0.5, size = 6, fontface = "bold") +
    scale_y_continuous(labels = scales::percent, limits = c(0, 1.1)) +
    scale_fill_gradient(low = infection_cols[["ASB"]], high = infection_cols[["UTI"]]) +
    labs(
        title = "UTI Risk by Sequence Type",
        subtitle = "Proportion of episodes that are Symptomatic UTI (vs ASB)",
        x = "Sequence Type",
        y = "UTI Risk (%)"
    ) +
    theme_minimal(base_size = 20) +
    theme(
      plot.title = element_text(size = 26, face = "bold"),
      plot.subtitle = element_text(size = 18, color = "grey40"),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 16, face = "bold"),
      axis.text.y = element_text(size = 16, face = "bold"),
      plot.margin = margin(15, 15, 15, 15)
    )

ggsave(file.path(out_dir, "st_risk_plot.png"), p, width = 12, height = 8, dpi = 600)

msg("Saved results to %s", out_dir)

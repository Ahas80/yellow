#!/usr/bin/env Rscript
# ==============================================================================
# 17_lineage_analysis.R
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
# Compare each ST against "All Other STs"
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
    geom_text(aes(label = sprintf("n=%d", n_total)), vjust = -0.5, size = 3) +
    scale_y_continuous(labels = scales::percent, limits = c(0, 1.1)) +
    scale_fill_gradient(low = infection_cols[["ASB"]], high = infection_cols[["UTI"]]) +
    labs(
        title = "UTI Risk by Sequence Type",
        subtitle = "Proportion of episodes that are Symptomatic UTI (vs ASB)",
        x = "Sequence Type",
        y = "UTI Risk (%)"
    ) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(out_dir, "st_risk_plot.png"), p, width = 8, height = 6)

msg("Saved results to %s", out_dir)

#!/usr/bin/env Rscript
# ==============================================================================
# 17_lineage_analysis.R
# ==============================================================================
#
# GOAL:
#   Compute UTI risk per Sequence Type: for each ST with sufficient data,
#   what proportion of its episodes are UTI?  This identifies STs that are
#   disproportionately associated with catheter-aware UTI episodes in this cohort.
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
#   - results/mlst/mlst_provider_preferred.csv
#   - results/models/gwas_multivariable_glmm.csv (to compare hits)
#
# Outputs:
#   - results/lineage/st_risk_profile.csv
#   - results/lineage/st_risk_plot.png
#
# Purpose:
#   - Determine if specific STs (e.g. ST131) are disproportionately represented
#     among UTI episodes in this cohort.
#   - Risk stratification of clones.
# ==============================================================================

source("00_config.R")
source("R/pipeline_qc_helpers.R")
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
status <- read_csv(FILE_STATUS_MAP, show_col_types = FALSE) %>%
    prefer_primary_uti_status() %>%
    filter(UTI_Status %in% c("UTI", "Not_UTI"))

# MLST
mlst_file <- FILE_MLST_CANONICAL
if (file.exists(mlst_file)) {
    mlst <- read_csv(mlst_file, show_col_types = FALSE)
} else {
    stop("Provider-preferred MLST missing: ", mlst_file, ". Run scripts/integrate_provider_mlst.R.")
}

# 2. Harmonize & Join
# ------------------------------------------------------------------------------
# We need to link ST to Episode (Participant + Timepoint)
# MLST is per Isolate. Status is per Episode.
# We'll take the ST of the isolate associated with the episode.

# Ensure columns match for join
mlst_clean <- mlst %>%
    mutate(Timepoint = if ("tp_lab" %in% names(.)) tp_lab else Timepoint) %>%
    select(any_of(c("Participant_id", "Timepoint", "ST", "ST_source", "ST_provider", "ST_local"))) %>%
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

if (file.exists(FILE_VF_READY)) {
    vf_ready_lineage <- read_csv(FILE_VF_READY, show_col_types = FALSE) %>%
        prefer_primary_uti_status() %>%
        mutate(
            Participant_id = as.character(Participant_id),
            Timepoint = normalise_timepoint_preserve_events(tp_lab),
            Infection_Status = as.character(UTI_Status),
            ST = as.character(ST)
        ) %>%
        filter(Infection_Status %in% c("UTI", "Not_UTI")) %>%
        filter(!is.na(ST), ST != "") %>%
        select(any_of(c("Participant_id", "Timepoint", "Episode_ID", "Infection_Status", "ST", "ST_source", "ST_provider", "ST_local", "uricult_bridge_applied"))) %>%
        distinct()

    if (nrow(vf_ready_lineage) > 0) {
        data_merged <- vf_ready_lineage
        msg(
            "Using canonical vf_analysis_ready.csv for lineage risk: %d episodes (%d UTI, %d Not_UTI; %d Uricult-bridged)",
            nrow(data_merged),
            sum(data_merged$Infection_Status == "UTI", na.rm = TRUE),
            sum(data_merged$Infection_Status == "Not_UTI", na.rm = TRUE),
            if ("uricult_bridge_applied" %in% names(data_merged)) sum(data_merged$uricult_bridge_applied %in% TRUE, na.rm = TRUE) else 0L
        )
    }
}

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
        n_Not_UTI = sum(Infection_Status == "Not_UTI"),
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
    #       UTI  Not_UTI
    # Target  a    b
    # Other   c    d

    a <- sum(df$ST == st_target & df$Infection_Status == "UTI")
    b <- sum(df$ST == st_target & df$Infection_Status == "Not_UTI")
    c <- sum(df$ST != st_target & df$Infection_Status == "UTI")
    d <- sum(df$ST != st_target & df$Infection_Status == "Not_UTI")

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
    scale_fill_gradient(low = uti_status_cols[["Not_UTI"]], high = uti_status_cols[["UTI"]]) +
    labs(
        title = "UTI Risk by Sequence Type",
        subtitle = "Proportion of VF/WGS-linked episodes classified as UTI versus Not_UTI",
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

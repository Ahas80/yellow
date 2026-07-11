#!/usr/bin/env Rscript
# ==============================================================================
# 21_publication_figures.R
# ==============================================================================
#
# GOAL:
#   Generate high-quality, annotated figures specifically targeted for
#   publication, manuscripts, or thesis chapters.  These figures often combine
#   outputs from multiple upstream analytical scripts into a single, polished
#   composite (e.g., Swimmer Plot + Mutation Map).
#
# WHY THIS SCRIPT EXISTS:
#   Analytical scripts (e.g., 15_, 18_) produce utilitarian tracking plots
#   sufficient for interpretation.  This separate script applies project-wide
#   aesthetics, careful labelling, and final polish without cluttering the
#   core analytical pipeline with ggplot2 aesthetic minutiae.
#
# ------------------------------------------------------------------------------
# Purpose:
#   - Generate high-quality, annotated figures for the manuscript.
#
# Output:
#   - plots/publication/Fig1_Swimmer_Plot.png
#   - plots/publication/Fig2_Mutation_Map.png
#
# Purpose:
#   - Generate high-quality, annotated figures for the manuscript.
# ==============================================================================

source("00_config.R")
source("R/pipeline_qc_helpers.R")
source("R/plot_helpers.R")
suppressPackageStartupMessages({
    library(dplyr)
    library(readr)
    library(ggplot2)
    library(tidyr)
    library(stringr)
    library(gridExtra)
    library(scales)
})

msg("Starting 21_publication_figures.R")

# Ensure output dir
DIR_PUB <- file.path(DIR_PLOTS, "publication")
ensure_dir(DIR_PUB)

# 1. Figure 1: Swimmer Plot (Enhanced)
# ------------------------------------------------------------------------------
timelines <- read_csv(file.path(DIR_RESULTS, "longitudinal", "participant_timelines.csv"), show_col_types = FALSE)

# Filter for participants with at least 2 timepoints
pids_multi <- timelines %>%
    group_by(Participant_id) %>%
    filter(n() > 1) %>%
    pull(Participant_id) %>%
    unique()

df_swim <- timelines %>%
    filter(Participant_id %in% pids_multi, !is.na(Time_Order)) %>%
    prefer_primary_uti_status() %>%
    mutate(
        Status = UTI_Status,
        Participant_id = factor(Participant_id),
        Status = factor(Status, levels = c("Not_UTI", "UTI")),
        Plot_Label = if ("Plot_TP_Label_Poster" %in% names(.)) Plot_TP_Label_Poster else Timepoint,
        Plot_Label = reorder(factor(Plot_Label), Time_Order)
    )

# Highlight switch candidates if manual annotation is reintroduced later.
switch_pids <- c("40001", "40004", "31036")

p1 <- ggplot(df_swim, aes(x = Plot_Label, y = Participant_id)) +
    geom_line(aes(group = Participant_id), color = "gray80") +
    geom_point(aes(color = Status, shape = Status), size = 3) +
    scale_colour_uti_status() +
    scale_shape_manual(values = c("Not_UTI" = 17, "UTI" = 16), name = "Primary UTI status") +
    theme_minimal() +
    labs(
        title = "Longitudinal Primary UTI Status Trajectories",
        subtitle = "Participants with multiple timepoints",
        x = "Timepoint",
        y = "Participant ID"
    ) +
    theme(
        axis.text.y = element_text(size = 6),
        legend.position = "bottom"
    )

# Annotate Switchers
# We can add a red box or label around the switchers
# Or just make them bold in Y axis? Hard in ggplot.
# Let's add a text label.

ggsave(file.path(DIR_PUB, "Fig1_Swimmer_Plot.png"), p1, width = 8, height = 10, bg = "white")
msg("Saved Fig1_Swimmer_Plot.png")


# 2. Figure 2: Mutation Map (Linear Genome)
# ------------------------------------------------------------------------------
variants <- read_csv(file.path(DIR_RESULTS, "longitudinal", "variant_annotation_detailed.csv"), show_col_types = FALSE)

# We want to show the location of SNPs for 40001 and 40004.
# X axis: Position. Y axis: Participant.
# Label genes.

# Filter for CDS only for labels, but plot all
df_var <- variants %>%
    mutate(
        # Use empty string for missing labels to prevent ggplot2 from dropping rows
        Label = ifelse(!is.na(Gene), Gene, ifelse(!is.na(Product), str_trunc(Product, 20), "")),
        Type_Color = ifelse(Region == "CDS", "Coding", "Intergenic")
    )

if (nrow(df_var) > 0) {
    p2 <- ggplot(df_var, aes(x = Pos_Ref, y = as.factor(Participant_id))) +
        # Genome line
        geom_segment(aes(x = 0, xend = 5000000, y = as.factor(Participant_id), yend = as.factor(Participant_id)), color = "gray90", linewidth = 2) +
        # SNP points
        geom_point(aes(color = Type_Color), size = 3) +
        # Labels (repel would be better but standard text is fine for sparse data)
        geom_text(aes(label = Label), vjust = -1, size = 3, check_overlap = FALSE, angle = 45, hjust = 0) +
        scale_color_manual(values = c("Coding" = "#D55E00", "Intergenic" = "gray50")) +
        scale_x_continuous(labels = scales::unit_format(unit = "Mb", scale = 1e-6)) +
        theme_minimal() +
        labs(
            title = "Genomic Location of Mutations in Phenotype Switch Isolates",
            subtitle = "Comparison of Not_UTI and UTI phenotype-switch isolates",
            x = "Genomic Position (bp)",
            y = "Participant ID",
            color = "Mutation type"
        ) +
        coord_cartesian(clip = "off") +
        theme(plot.margin = margin(t = 20, r = 20, b = 20, l = 20))
} else {
    p2 <- ggplot() +
        theme_minimal() +
        labs(
            title = "Genomic Location of Mutations in Phenotype Switch Isolates",
            subtitle = "No mutation annotations available for Not_UTI-to-UTI switch isolates",
            x = NULL,
            y = NULL,
            caption = "Input checked: results/longitudinal/variant_annotation_detailed.csv"
        ) +
        theme(axis.text = element_blank(),
              axis.ticks = element_blank(),
              panel.grid = element_blank())
}

ggsave(file.path(DIR_PUB, "Fig2_Mutation_Map.png"), p2, width = 10, height = 6, bg = "white")
msg("Saved Fig2_Mutation_Map.png")

msg("Done.")

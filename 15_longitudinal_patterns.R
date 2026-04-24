#!/usr/bin/env Rscript
# ==============================================================================
# 15_longitudinal_patterns.R
# ==============================================================================
#
# GOAL:
#   Construct longitudinal timelines for each participant: assign strain IDs
#   based on pairwise similarity, track strain persistence duration, and
#   identify "phenotype-switch" events where the same strain changes clinical
#   status (e.g., ASB at T0 → UTI at T2).  These switch events are the
#   candidates for deep within-host evolution analysis in script 16.
#
# ------------------------------------------------------------------------------
# Role: [Analysis] - Construct longitudinal timelines and identify transitions.
#
# Inputs:
#   - results/clinical/status_map.csv
#   - results/strain_compare/pairwise_metrics.csv
#   - results/clinical/intermediate/clinical_merged.rds (for dates if avail)
#
# Outputs:
#   - results/longitudinal/participant_timelines.csv
#   - results/longitudinal/transitions.csv
#   - results/longitudinal/swimmer_plot.png
#
# Purpose:
#   - Assigns global "Strain IDs" to isolates based on pairwise "Same" clusters.
#   - Tracks persistence duration.
#   - Identifies specific "Phenotype Switch" events (Same Strain, ASB <-> UTI)
#     for downstream evolutionary analysis.
# ==============================================================================

source("00_config.R")
source("R/plot_helpers.R")
suppressPackageStartupMessages({
    library(dplyr)
    library(readr)
    library(tidyr)
    library(stringr)
    library(igraph)
    library(ggplot2)
    library(lubridate)
    library(forcats)
})

msg("Starting 15_longitudinal_patterns.R")

# 1. Load Data
# ------------------------------------------------------------------------------
# Clinical Status
status_file <- file.path(DIR_CLINICAL, "status_map_with_poster_tp.csv")
if (!file.exists(status_file)) status_file <- file.path(DIR_CLINICAL, "status_map.csv")
status_map <- read_csv(status_file, show_col_types = FALSE)
msg("Loaded %d episodes from %s", nrow(status_map), basename(status_file))

# Pairwise Metrics
pairwise_file <- file.path(DIR_STRAIN, "pairwise_metrics.csv")
if (!file.exists(pairwise_file)) stop("Missing pairwise_metrics.csv. Run 11_compare_strains.R first.")
pairwise <- read_csv(pairwise_file, show_col_types = FALSE)
msg("Loaded %d pairwise comparisons", nrow(pairwise))

# 2. Assign Global Strain IDs (Graph Clustering)
# ------------------------------------------------------------------------------
msg("Assigning Strain IDs based on 'Same' classification...")

# [EPI] Why Graph Clustering?
# The pairwise comparisons in Script 11 only tell us if Isolate A == Isolate B.
# To track a strain over time (A == B, B == C), we treat each "Same" classification
# as an edge in an undirected graph. The connected components of this graph
# represent a single contiguous strain lineage (Strain_ID) persisting within the host.

# Filter for "Same" edges
edges <- pairwise %>%
    filter(Classification == "Same") %>%
    select(from = SampleKey_A, to = SampleKey_B)

# Get all unique samples from status_map to ensure singletons are included
# We need to construct SampleKey for status_map entries to match pairwise
# SampleKey format: Participant_id__Timepoint (as used in 11_compare_strains.R)
status_map <- status_map %>%
    mutate(SampleKey = paste0(Participant_id, "__", Timepoint))

all_nodes <- unique(status_map$SampleKey)

# Build graph
g <- graph_from_data_frame(edges, directed = FALSE, vertices = data.frame(name = all_nodes))

# Find clusters
comps <- components(g)
membership <- comps$membership

# Create mapping table
strain_map <- tibble(
    SampleKey = names(membership),
    Cluster_ID = membership
) %>%
    mutate(Strain_ID = paste0("Strain_", Cluster_ID))

# Join back to status_map
timeline <- status_map %>%
    left_join(strain_map, by = "SampleKey", relationship = "many-to-many") %>%
    select(-Cluster_ID)

msg("Assigned %d unique Strain IDs across %d episodes", n_distinct(timeline$Strain_ID, na.rm = TRUE), nrow(timeline))

# [Transparency Check] Report episodes with missing Strain IDs (QC Failures/No Seq)
missing_strain <- sum(is.na(timeline$Strain_ID))
if (missing_strain > 0) {
    msg(
        "NOTE: %d episodes (%.1f%%) have no assigned Strain ID (Sequencing failed or unavailable). These are retained in the timeline as 'No_Seq' but excluded from strain comparisons.",
        missing_strain,
        missing_strain / nrow(timeline) * 100
    )
}

# 3. Order Timelines & Calculate Transitions
# ------------------------------------------------------------------------------
# We need to order timepoints.
# Logic: T0 < T1 < T2 < ... < Uricult (treat Uricult as distinct or last?)
# Let's try to extract numeric time or use a factor level.

parse_time_order <- function(tp) {
    case_when(
        tp == "T0" ~ 0,
        tp == "T1" ~ 1,
        tp == "T2" ~ 2,
        tp == "T3" ~ 3,
        tp == "T4" ~ 4,
        str_detect(tp, "(?i)uricult") ~ 99, # Put Uricult last or treat as separate
        TRUE ~ 999
    )
}

timeline <- timeline %>%
    mutate(
        Time_Order = if ("Plot_TP_Num_Poster" %in% names(.)) {
            Plot_TP_Num_Poster
        } else {
            parse_time_order(Timepoint)
        },
        # Clean Infection Status for plotting
        Status_Simple = case_when(
            Infection_Status == "UTI" ~ "UTI",
            Infection_Status == "ASB" ~ "ASB",
            Infection_Status == "Negative" ~ "Negative",
            TRUE ~ "Other"
        )
    ) %>%
    # WARNING: Uricults without a Plot_TP_Num_Poster are excluded from ordered analysis
    filter(!is.na(Time_Order)) %>%
    arrange(Participant_id, Time_Order)

# Identify Transitions
transitions <- timeline %>%
    group_by(Participant_id) %>%
    mutate(
        Prev_Strain = lag(Strain_ID),
        Prev_Status = lag(Status_Simple),
        Prev_Time = lag(Timepoint),

        # Transition Types
        Is_Consecutive = !is.na(Prev_Strain), # Is there a previous point?
        Same_Strain = Strain_ID == Prev_Strain,
        Status_Change = Status_Simple != Prev_Status,
        Event_Type = case_when(
            !Is_Consecutive ~ "Initial",
            Same_Strain & Status_Change ~ "Phenotype_Switch",
            Same_Strain & !Status_Change ~ "Persistence",
            !Same_Strain ~ "Strain_Replacement",
            TRUE ~ "Unknown"
        )
    ) %>%
    ungroup()

# 4. Extract "Phenotype Switch" Candidates (Priority 1)
# ------------------------------------------------------------------------------
switch_events <- transitions %>%
    filter(Event_Type == "Phenotype_Switch") %>%
    select(Participant_id,
        From_Time = Prev_Time, To_Time = Timepoint,
        From_Status = Prev_Status, To_Status = Status_Simple, Strain_ID
    )

msg("Found %d 'Phenotype Switch' events (Same Strain, Different Status)", nrow(switch_events))
if (nrow(switch_events) > 0) {
    print(head(switch_events))
}

# 5. Calculate Persistence Stats (Priority 2)
# ------------------------------------------------------------------------------
# For each Strain_ID in a participant, how many timepoints?
persistence_stats <- timeline %>%
    group_by(Participant_id, Strain_ID) %>%
    summarise(
        n_timepoints = n(),
        timepoints = paste(Timepoint, collapse = "->"),
        statuses = paste(Status_Simple, collapse = "->"),
        has_UTI = any(Status_Simple == "UTI"),
        has_ASB = any(Status_Simple == "ASB"),
        .groups = "drop"
    ) %>%
    arrange(desc(n_timepoints))

# 6. Save Outputs
# ------------------------------------------------------------------------------
out_dir <- file.path(DIR_RESULTS, "longitudinal")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

write_csv(timeline, file.path(out_dir, "participant_timelines.csv"))
write_csv(transitions, file.path(out_dir, "transitions.csv"))
write_csv(switch_events, file.path(out_dir, "phenotype_switch_candidates.csv"))
write_csv(persistence_stats, file.path(out_dir, "strain_persistence_stats.csv"))

msg("Saved tables to %s", out_dir)

# 7. Visualization: Swimmer Plot
# ------------------------------------------------------------------------------
# Filter to participants with > 1 timepoint for cleaner plot
plot_data <- timeline %>%
    filter(!is.na(Time_Order)) %>%
    group_by(Participant_id) %>%
    filter(n() > 1) %>%
    ungroup()

if (nrow(plot_data) > 0) {
    # Ensure Timepoint is ordered correctly
    plot_data <- plot_data %>%
        mutate(
            Plot_Label = if ("Plot_TP_Label_Poster" %in% names(.)) Plot_TP_Label_Poster else Timepoint,
            Plot_Label = fct_reorder(Plot_Label, Time_Order)
        )

    p <- ggplot(plot_data, aes(x = Plot_Label, y = Participant_id)) +
        # Line connecting points
        geom_line(aes(group = Participant_id), color = "grey80") +
        # Points colored by Status, shaped by Strain (if feasible, or just Status)
        geom_point(aes(color = Status_Simple, shape = Status_Simple), size = 3) +
        # Add Strain ID labels (optional, might be cluttered)
        # geom_text(aes(label = Strain_ID), vjust = -0.5, size = 2) +
        scale_colour_infection() +
        # Use filled shapes so colors are visible (not hollow circles)
        scale_shape_manual(
            name = "Infection Status",
            values = c("UTI" = 16, "ASB" = 17, "Negative" = 15), # 16=filled circle, 17=filled triangle, 15=filled square
            breaks = c("UTI", "ASB", "Negative")
        ) +
        theme_minimal() +
        labs(
            title = "Longitudinal Infection Status (Participants > 1 Timepoint)",
            subtitle = "Participants with >1 timepoint",
            x = "Timepoint",
            y = "Participant",
            color = "Infection Status",
            shape = "Infection Status"
        ) +
        theme(axis.text.y = element_text(size = 6))

    ggsave(file.path(out_dir, "swimmer_plot.png"), p, width = 10, height = 12, dpi = 300)
    msg("Saved swimmer_plot.png")
}

msg("Analysis Complete.")

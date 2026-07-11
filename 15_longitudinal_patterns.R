#!/usr/bin/env Rscript
# ==============================================================================
# 15_longitudinal_patterns.R
# ==============================================================================
#
# GOAL:
#   Construct longitudinal timelines for each participant: assign strain IDs
#   based on pairwise similarity, track strain persistence duration, and
#   identify "phenotype-switch" events where the same strain changes primary
#   primary UTI status (e.g., Not_UTI at T0 -> UTI at T2). These switch events are the
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
#   - Identifies specific "Phenotype Switch" events (Same Strain, Not_UTI <-> UTI)
#     for downstream evolutionary analysis.
# ==============================================================================

source("00_config.R")
source("R/pipeline_qc_helpers.R")
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
# Clinical Status. Prefer the poster-timepoint file only when it is fresh and
# already carries the primary UTI_Status columns; otherwise use status_map.csv so
# stale legacy ASB/UTI/Negative labels cannot re-enter timelines.
status_map <- read_primary_status_map(
    prefer_poster = TRUE,
    require_fresh = TRUE,
    caller = "15_longitudinal_patterns.R"
)
status_file <- attr(status_map, "source_file") %||% FILE_STATUS_MAP
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
    transmute(
        from = as.character(SampleKey_A),
        to = as.character(SampleKey_B)
    ) %>%
    filter(!is.na(from), !is.na(to), nzchar(from), nzchar(to), from != to) %>%
    distinct()

# Get all unique samples from status_map to ensure singletons are included
# We need to construct SampleKey for status_map entries to match pairwise
# SampleKey format: Participant_id__Timepoint (as used in 11_compare_strains.R)
status_map <- status_map %>%
    mutate(
        tp_lab = if ("tp_lab" %in% names(.)) normalise_timepoint_preserve_events(tp_lab) else normalise_timepoint_preserve_events(Timepoint),
        SampleKey = paste0(Participant_id, "__", tp_lab)
    )

timeline_nodes <- unique(status_map$SampleKey)
edge_nodes <- unique(c(edges$from, edges$to))
edge_only_nodes <- setdiff(edge_nodes, timeline_nodes)

if (length(edge_only_nodes) > 0) {
    msg(
        "NOTE: %d same-strain graph node(s) occur in pairwise_metrics.csv but not in the clinical timeline. Including them as graph-only vertices.",
        length(edge_only_nodes)
    )
    msg("First graph-only nodes: %s", paste(head(edge_only_nodes, 10), collapse = ", "))
}

all_nodes <- unique(c(timeline_nodes, edge_nodes))

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
    left_join(strain_map, by = "SampleKey", relationship = "many-to-one") %>%
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
    num <- suppressWarnings(readr::parse_number(tp))
    case_when(
        str_detect(tp, regex("^T[0-9]+$", ignore_case = TRUE)) ~ as.numeric(num),
        str_detect(tp, regex("uricult", ignore_case = TRUE)) ~ NA_real_,
        TRUE ~ NA_real_
    )
}

timeline <- timeline %>%
    mutate(
        Collection_Date_parsed = if ("Collection_Date" %in% names(.)) suppressWarnings(lubridate::dmy(Collection_Date)) else as.Date(NA),
        Time_Order_Date = ifelse(!is.na(Collection_Date_parsed),
                                 as.numeric(Collection_Date_parsed - ave(Collection_Date_parsed, Participant_id, FUN = function(x) min(x, na.rm = TRUE))),
                                 NA_real_),
        Time_Order_Fallback = parse_time_order(tp_lab),
        Time_Order = coalesce(Time_Order_Date, Time_Order_Fallback),
        Time_Order_Source = case_when(
            !is.na(Time_Order_Date) ~ "Collection_Date",
            !is.na(Time_Order_Fallback) ~ "Routine_timepoint_number",
            TRUE ~ "Unavailable"
        ),
        # Clean primary UTI status for plotting
        Status_Simple = case_when(
            UTI_Status == "UTI" ~ "UTI",
            UTI_Status == "Not_UTI" ~ "Not_UTI",
            TRUE ~ "Other"
        )
    ) %>%
    # Poster half-step labels are display-only; they are not used as statistical covariates here.
    filter(!is.na(Time_Order)) %>%
    arrange(Participant_id, Time_Order, tp_lab)

# Identify Transitions
transitions <- timeline %>%
    group_by(Participant_id) %>%
    mutate(
        Prev_Strain = lag(Strain_ID),
        Prev_Status = lag(Status_Simple),
        Prev_Time = lag(tp_lab),

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
        From_Time = Prev_Time, To_Time = tp_lab,
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
        has_Not_UTI = any(Status_Simple == "Not_UTI"),
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

not_uti_uti_transition_summary <- transitions %>%
    summarise(
        clinical_Not_UTI_to_UTI = sum(Prev_Status == "Not_UTI" & Status_Simple == "UTI", na.rm = TRUE),
        same_strain_Not_UTI_to_UTI = sum(Prev_Status == "Not_UTI" & Status_Simple == "UTI" & Same_Strain %in% TRUE, na.rm = TRUE),
        uricult_related_Not_UTI_to_UTI = sum(Prev_Status == "Not_UTI" & Status_Simple == "UTI" &
                                             (str_detect(Prev_Time, regex("uricult", ignore_case = TRUE)) |
                                                str_detect(tp_lab, regex("uricult", ignore_case = TRUE))), na.rm = TRUE)
    )
write_csv(not_uti_uti_transition_summary, file.path(out_dir, "not_uti_uti_transition_summary.csv"))

append_denominator_summary(
    timeline,
    "15_longitudinal_patterns.R",
    "participant_timelines",
    "clinical_episode",
    file.path(out_dir, "participant_timelines.csv"),
    "Clinical longitudinal ordering uses Collection_Date where available; poster half-step labels are display-only"
)

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
        scale_colour_uti_status(na.value = "grey75") +
        # Use filled shapes so colors are visible (not hollow circles)
        scale_shape_manual(
            name = "Primary UTI status",
            values = c("UTI" = 16, "Not_UTI" = 17, "Other" = 15),
            breaks = c("UTI", "Not_UTI", "Other")
        ) +
        theme_minimal() +
        labs(
            title = "Longitudinal Primary UTI Status (Participants > 1 Timepoint)",
            subtitle = "Participants with >1 timepoint",
            x = "Timepoint",
            y = "Participant",
            color = "Primary UTI status",
            shape = "Primary UTI status"
        ) +
        theme(axis.text.y = element_text(size = 6))

    ggsave(file.path(out_dir, "swimmer_plot.png"), p, width = 10, height = 12, dpi = 300)
    msg("Saved swimmer_plot.png")
}

msg("Analysis Complete.")

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
#   - results/clinical/analysis_cohort_longcycler.csv
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
# The analytical timeline is the selected QC-pass Longcycler cohort. The broader
# clinical status map is source attrition/QC only and must not enter this script.
if (!file.exists(FILE_ANALYSIS_CLINICAL_COHORT)) {
    stop("Missing ", FILE_ANALYSIS_CLINICAL_COHORT,
         ". Run the canonical assembly/QC cohort build before script 15.")
}
status_map <- readr::read_csv(FILE_ANALYSIS_CLINICAL_COHORT, show_col_types = FALSE) %>%
    prefer_primary_uti_status(allow_legacy_fallback = FALSE) %>%
    mutate(
        Participant_id = as.character(.data$Participant_id),
        tp_lab = normalise_timepoint_preserve_events(.data$tp_lab)
    )
assert_unique_keys(
    status_map, c("Participant_id", "tp_lab"),
    context = "15 Longcycler analysis cohort",
    out_path = file.path(DIR_QC, "15_duplicate_analysis_cohort_keys.csv")
)
msg("Loaded %d selected Longcycler episodes from %s",
    nrow(status_map), basename(FILE_ANALYSIS_CLINICAL_COHORT))

# Pairwise Metrics
pairwise_file <- file.path(DIR_STRAIN, "pairwise_metrics.csv")
if (!file.exists(pairwise_file)) stop("Missing pairwise_metrics.csv. Run 11_compare_strains.R first.")
pairwise <- read_csv(pairwise_file, show_col_types = FALSE) %>%
    mutate(
        Participant_id_A = as.character(.data$Participant_id_A),
        Participant_id_B = as.character(.data$Participant_id_B),
        Timepoint_A = normalise_timepoint_preserve_events(.data$Timepoint_A),
        Timepoint_B = normalise_timepoint_preserve_events(.data$Timepoint_B),
        SampleKey_A = paste0(.data$Participant_id_A, "__", .data$Timepoint_A),
        SampleKey_B = paste0(.data$Participant_id_B, "__", .data$Timepoint_B),
        pair_key = if_else(
            .data$SampleKey_A <= .data$SampleKey_B,
            paste(.data$SampleKey_A, .data$SampleKey_B, sep = "||"),
            paste(.data$SampleKey_B, .data$SampleKey_A, sep = "||")
        )
    )
if (anyDuplicated(pairwise$pair_key)) stop("pairwise_metrics.csv has duplicate unordered endpoint pairs")
if (any(tolower(pairwise$Assembler_A) != ANALYSIS_ASSEMBLER |
        tolower(pairwise$Assembler_B) != ANALYSIS_ASSEMBLER, na.rm = TRUE)) {
    stop("pairwise_metrics.csv contains an endpoint outside the selected Longcycler cohort")
}
msg("Loaded %d pairwise comparisons", nrow(pairwise))

canonical_transition_file <- file.path(DIR_RESULTS, "longitudinal", "longcycler_transitions.csv")
if (!file.exists(canonical_transition_file)) {
    stop("Missing canonical selected-Longcycler transition export: ", canonical_transition_file)
}
canonical_transitions <- read_csv(canonical_transition_file, show_col_types = FALSE) %>%
    mutate(
        Participant_id = as.character(.data$Participant_id),
        tp_from = normalise_timepoint_preserve_events(.data$tp_from),
        tp_to = normalise_timepoint_preserve_events(.data$tp_to),
        from_sample_key = paste0(.data$Participant_id, "__", .data$tp_from),
        to_sample_key = paste0(.data$Participant_id, "__", .data$tp_to),
        pair_key = if_else(
            .data$from_sample_key <= .data$to_sample_key,
            paste(.data$from_sample_key, .data$to_sample_key, sep = "||"),
            paste(.data$to_sample_key, .data$from_sample_key, sep = "||")
        )
    )
if (nrow(canonical_transitions) != 371L || anyDuplicated(canonical_transitions$pair_key)) {
    stop("Canonical selected-Longcycler transition export must contain 371 unique adjacent pairs")
}
canonical_transition_join_cols <- c(
    "pair_key", "TotalSNPs", "AvgIdentity", "strict_same_strain",
    "legacy_accessory_composite_classification",
    "legacy_accessory_composite_rule",
    "snp_strain_context", "st_lineage_context", "pair_interpretation",
    "Assembler_A", "Assembler_B", "Fasta_SHA256_A", "Fasta_SHA256_B"
)
missing_transition_join_cols <- setdiff(
    canonical_transition_join_cols,
    names(canonical_transitions)
)
if (length(missing_transition_join_cols)) {
    stop(
        "Canonical selected-Longcycler transition export lacks required join column(s): ",
        paste(missing_transition_join_cols, collapse = ", ")
    )
}

# 2. Assign Global Strain IDs (Graph Clustering)
# ------------------------------------------------------------------------------
msg("Assigning strain-context components from direct SNP-supported edges...")

# Graph components provide lineage context only. Every adjacent transition below
# is classified from its own direct pairwise SNP record; component membership is
# never substituted for missing direct evidence and never defines Same_Strain.

# Only direct pairwise SNP evidence defines a same-strain edge. The historical
# accessory composite remains available as a non-canonical compatibility field.
edges <- pairwise %>%
    filter(strict_same_strain %in% TRUE) %>%
    transmute(
        from = as.character(SampleKey_A),
        to = as.character(SampleKey_B)
    ) %>%
    filter(!is.na(from), !is.na(to), nzchar(from), nzchar(to), from != to) %>%
    distinct()

status_map <- status_map %>%
    mutate(
        SampleKey = paste0(Participant_id, "__", tp_lab)
    )

timeline_nodes <- unique(status_map$SampleKey)
edges <- edges %>% filter(.data$from %in% timeline_nodes, .data$to %in% timeline_nodes)
all_nodes <- timeline_nodes

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

msg("Assigned %d graph-component context IDs across %d selected Longcycler episodes",
    n_distinct(timeline$Strain_ID, na.rm = TRUE), nrow(timeline))

# 3. Order Timelines & Calculate Transitions
# ------------------------------------------------------------------------------
# We need to order timepoints.
# Logic: T0 < T1 < T2 < ... < Uricult (treat Uricult as distinct or last?)
# Let's try to extract numeric time or use a factor level.

parse_time_order <- function(tp) {
    num <- suppressWarnings(readr::parse_number(tp))
    case_when(
        str_detect(tp, regex("^T[0-9]+$", ignore_case = TRUE)) ~ as.numeric(num),
        str_detect(tp, regex("uricult", ignore_case = TRUE)) ~ 99,
        str_detect(tp, regex("^UTI-[0-9]+$", ignore_case = TRUE)) ~ 100 + as.numeric(num),
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
if (nrow(timeline) != nrow(status_map)) {
    stop("Every selected Longcycler episode must have a usable longitudinal order")
}

# Identify adjacent selected-episode transitions and attach the direct pair
# record in orientation-independent form.
transitions <- timeline %>%
    group_by(Participant_id) %>%
    mutate(
        Prev_Strain = lag(Strain_ID),
        Prev_Status = lag(Status_Simple),
        Prev_Time = lag(tp_lab),
        Is_Consecutive = !is.na(Prev_Time)
    ) %>%
    ungroup() %>%
    filter(.data$Is_Consecutive) %>%
    mutate(
        from_sample_key = paste0(.data$Participant_id, "__", .data$Prev_Time),
        to_sample_key = paste0(.data$Participant_id, "__", .data$tp_lab),
        pair_key = if_else(
            .data$from_sample_key <= .data$to_sample_key,
            paste(.data$from_sample_key, .data$to_sample_key, sep = "||"),
            paste(.data$to_sample_key, .data$from_sample_key, sep = "||")
        )
    ) %>%
    left_join(
        canonical_transitions %>%
            select(all_of(canonical_transition_join_cols)),
        by = "pair_key", relationship = "many-to-one"
    ) %>%
    mutate(
        Direct_Pair_Evidence = !is.na(.data$TotalSNPs),
        Same_Strain = .data$Direct_Pair_Evidence & .data$TotalSNPs <= strain_snp_threshold(),
        Same_Graph_Component = .data$Strain_ID == .data$Prev_Strain,
        Status_Change = .data$Status_Simple != .data$Prev_Status,
        Event_Type = case_when(
            !.data$Direct_Pair_Evidence ~ "Missing_Direct_Pair_Evidence",
            .data$Same_Strain & .data$Status_Change ~ "Phenotype_Switch",
            .data$Same_Strain & !.data$Status_Change ~ "Persistence",
            !.data$Same_Strain ~ "Strain_Replacement",
            TRUE ~ "Unknown"
        )
    )

expected_transitions <- timeline %>%
    count(.data$Participant_id, name = "n_episodes") %>%
    summarise(n = sum(pmax(.data$n_episodes - 1L, 0L))) %>%
    pull(.data$n)
if (nrow(transitions) != expected_transitions) {
    stop("Adjacent transition count is inconsistent with the selected Longcycler timeline")
}
if (nrow(transitions) != 371L ||
    sum(transitions$Prev_Status == "Not_UTI" & transitions$Status_Simple == "UTI", na.rm = TRUE) != 9L) {
    stop("Selected Longcycler transition contract failed: expected 371 adjacent pairs and 9 Not_UTI-to-UTI pairs")
}
if (any(!transitions$Direct_Pair_Evidence)) {
    stop("One or more selected Longcycler adjacent transitions lack direct pairwise SNP evidence")
}
if (any(tolower(transitions$Assembler_A) != ANALYSIS_ASSEMBLER |
        tolower(transitions$Assembler_B) != ANALYSIS_ASSEMBLER)) {
    stop("A selected timeline transition contains a non-Longcycler pairwise endpoint")
}

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
    "selected_longcycler_episode",
    file.path(out_dir, "participant_timelines.csv"),
    "Exact selected QC-pass Longcycler cohort; direct pair evidence defines adjacent transition strain context"
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

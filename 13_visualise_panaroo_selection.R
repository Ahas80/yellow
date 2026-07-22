#!/usr/bin/env Rscript
# ==============================================================================
# 13_visualise_panaroo_selection.R
# ==============================================================================
#
# GOAL:
#   Visualise which assemblies were selected for Panaroo input and check
#   whether the selection introduces any bias (e.g., over-representing
#   certain batches, statuses, or assemblers).  This is a QC step to
#   ensure pangenome analysis is not biased by input selection.
#
# ------------------------------------------------------------------------------
# Role: [Visualization] - Visualize Panaroo input selection and QC bias.
#
# Inputs:
#   - results/wgs/qc_summary.csv
#   - results/status_map.csv
#
# Outputs:
#   - results/wgs/panaroo_selection_summary.csv
#   - results/wgs/panaroo_selection_detailed.csv
#   - results/wgs/panaroo_timepoint_selection_summary.csv
#   - plots/wgs/panaroo_selection_matrix.png
#   - plots/wgs/panaroo_selection_overview.png
#   - plots/wgs/panaroo_timepoint_selection_bar.png
#   - plots/wgs/panaroo_timepoint_selection_tile.png
#   - plots/wgs/panaroo_selection_samples_per_timepoint.png
#
# Usage:
#   Rscript 13_visualise_panaroo_selection.R
#
# Biological/Statistical purpose:
#   - Visualizes which samples passed QC and were selected for pangenome analysis.
#   - Checks for selection bias (e.g., are UTI samples more likely to fail QC?).
# ==============================================================================

# 1. Load Configuration & Libraries
source("00_config.R")
source("R/wgs_helpers.R")
source("R/pipeline_qc_helpers.R")
source("R/plot_helpers.R")

suppressPackageStartupMessages({
    library(tidyverse)
    library(fs)
    library(ggplot2)
})

# The live/older complete-analysis runner may already have been parsed before a
# post-Panaroo cleanup hook was added.  Enforce the same boundary here, before
# any downstream result is read or written.  A cleanup summary older than the
# newly selected cohort proves that only the start-of-run sweep has happened.
cleanup_script <- file.path("scripts", "prepare_longcycler_release.R")
cleanup_summary <- file.path(DIR_RESULTS, "pipeline", "release_cleanup_summary.csv")
selected_cohort <- FILE_ANALYSIS_CLINICAL_COHORT
if (!file.exists(selected_cohort)) {
    stop("Selected Longcycler cohort is missing before downstream cleanup: ", selected_cohort)
}
selected_cohort_mtime <- file.info(selected_cohort)$mtime
cleanup_summary_mtime <- file.info(cleanup_summary)$mtime
cleanup_is_current <- file.exists(cleanup_summary) &&
    !is.na(cleanup_summary_mtime) &&
    !is.na(selected_cohort_mtime) &&
    cleanup_summary_mtime >= selected_cohort_mtime
if (!cleanup_is_current) {
    if (!file.exists(cleanup_script)) {
        stop("Pre-downstream release cleanup script is missing: ", cleanup_script)
    }
    message("Applying required post-Panaroo Longcycler-only generated-output sweep...")
    cleanup_status <- system2(
        file.path(R.home("bin"), "Rscript"),
        c(cleanup_script, "--apply")
    )
    if (!identical(as.integer(cleanup_status), 0L)) {
        stop("Post-Panaroo Longcycler-only cleanup failed with exit status ", cleanup_status, ".")
    }
}
cleanup_summary_mtime <- file.info(cleanup_summary)$mtime
cleanup_is_current <- file.exists(cleanup_summary) &&
    !is.na(cleanup_summary_mtime) &&
    !is.na(selected_cohort_mtime) &&
    cleanup_summary_mtime >= selected_cohort_mtime
if (!cleanup_is_current) {
    stop("Post-Panaroo cleanup did not publish a current release-cleanup summary.")
}

canon_tp <- function(x) {
    normalise_timepoint_preserve_events(x)
}

format_qc_reason <- function(qc_pass, qc_reason) {
    reason <- str_trim(as.character(qc_reason))
    reason <- str_remove(reason, ";+$")
    reason <- str_replace_all(reason, c(
        "HighContigs" = "high contig count",
        "LowN50" = "low N50",
        "BadSize" = "assembly size outside QC range",
        "ReadError" = "assembly read error"
    ))
    reason <- str_replace_all(reason, ";+", " + ")
    case_when(
        qc_pass %in% TRUE ~ "QC pass",
        is.na(qc_pass) ~ "QC status unavailable",
        is.na(reason) | !nzchar(reason) ~ "QC failure: reason unavailable",
        TRUE ~ paste0("QC failure: ", reason)
    )
}

# 2. Load Centralized QC Data (from 12a_wgs_qc.R)
# ------------------------------------------------------------------------------
qc_file <- file.path(DIR_WGS, "qc_summary.csv")

if (!file.exists(qc_file)) {
    stop("QC summary not found. Please run 12a_wgs_qc.R first.")
}

# Load QC summary
qc_df <- read_csv(qc_file, show_col_types = FALSE)

# Prepare 'df' for plotting
# We map QC_PASS to "Selected" and QC_REASON to "Reason"
df <- qc_df %>%
    mutate(
        Assembly_ID = if ("Assembly_ID" %in% names(.)) Assembly_ID else tools::file_path_sans_ext(basename(full_path)),
        tp_lab = if ("tp_lab" %in% names(.)) normalise_timepoint_preserve_events(tp_lab) else normalise_timepoint_preserve_events(Timepoint),
        Timepoint = tp_lab,
        Selected = QC_PASS,
        Reason = ifelse(QC_PASS, "QC PASS", QC_REASON)
    )

canonical_file <- file.path(DIR_QC, "canonical_assembly_selection.csv")
if (!file.exists(canonical_file)) stop("Canonical Longcycler selection is missing: ", canonical_file)
canonical <- read_csv(canonical_file, show_col_types = FALSE) %>%
    select(any_of(c("Assembly_ID", "selected_canonical", "canonical_reason")))
df <- df %>%
    left_join(canonical, by = "Assembly_ID", relationship = "one-to-one") %>%
    mutate(selected_canonical = coalesce(as_pipeline_bool(selected_canonical), FALSE))

pan_manifest_file <- file.path(DIR_WGS, "pan", "panaroo_input_manifest.csv")
if (!file.exists(pan_manifest_file)) stop("Strict Panaroo input manifest is missing: ", pan_manifest_file)
pan_manifest <- read_csv(pan_manifest_file, show_col_types = FALSE) %>%
    select(any_of(c("Assembly_ID", "gff_available", "gff_path")))
df <- df %>%
    left_join(pan_manifest, by = "Assembly_ID", relationship = "one-to-one") %>%
    mutate(
        gff_available = coalesce(as_pipeline_bool(gff_available), FALSE),
        included_in_current_panaroo = selected_canonical & gff_available
    )

# 3. Summary
selected_count <- sum(df$selected_canonical, na.rm = TRUE)
included_count <- sum(df$included_in_current_panaroo, na.rm = TRUE)
message(sprintf("Canonical selected %d out of %d assemblies; %d have GFFs and are eligible for current Panaroo input.",
                selected_count, nrow(df), included_count))

# 4. Export List (Optional, as 12a already did it, but keeping for compatibility)
out_csv <- file.path(DIR_WGS, "panaroo_selection_summary.csv")
write_csv(df, out_csv)
message(sprintf("Summary saved to %s", out_csv))

# Order timepoints logically for plotting
df <- df %>%
    mutate(
        tp_num = suppressWarnings(as.integer(str_extract(Timepoint, "\\d+"))),
        tp_num = ifelse(is.na(tp_num), 999, tp_num), # Put non-numeric at end
        Timepoint = fct_reorder(Timepoint, tp_num)
    )

# Plot 3: Participant vs Timepoint Status
# We want to see for each person, what they have. Reader-facing plot labels are
# kept separate so compatibility tables retain the raw QC_REASON field.
participant_levels <- df %>%
    distinct(Participant_id) %>%
    mutate(
        participant_text = as.character(Participant_id),
        participant_number = suppressWarnings(as.numeric(participant_text))
    ) %>%
    arrange(is.na(participant_number), participant_number, participant_text) %>%
    pull(participant_text)
panaroo_matrix_df <- df %>%
    mutate(
        Reason_plot = format_qc_reason(QC_PASS, QC_REASON),
        Participant_plot = factor(as.character(Participant_id),
                                  levels = rev(participant_levels))
    )
reason_levels <- unique(panaroo_matrix_df$Reason_plot)
reason_levels <- c(
    intersect("QC pass", reason_levels),
    sort(setdiff(reason_levels, "QC pass"))
)
reason_colours <- setNames(
    ifelse(
        reason_levels == "QC pass", "#0072B2",
        ifelse(reason_levels == "QC status unavailable", "#6A6A6A", "#D55E00")
    ),
    reason_levels
)
assert_ruti_scale_levels(
    panaroo_matrix_df$Reason_plot, reason_colours,
    context = "Panaroo QC-selection matrix fill",
    allow_na = FALSE,
    require_all_palette_levels = TRUE
)

p3 <- ggplot(panaroo_matrix_df, aes(x = Timepoint, y = Participant_plot, fill = Reason_plot)) +
    geom_tile(color = "white", linewidth = 0.2) +
    scale_fill_manual(values = reason_colours, breaks = reason_levels) +
    labs(
        title = "Panaroo Sample Selection Matrix",
        subtitle = "Blue = QC pass; vermillion = QC failure (reader-facing reason)",
        x = "Timepoint",
        y = "Participant ID",
        fill = "Assembly QC outcome"
    ) +
    theme_minimal() +
    theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text.y = element_text(size = 6),
        panel.grid = element_blank()
    )

plot_file3 <- file.path(DIR_PLOTS_WGS, "panaroo_selection_matrix.png")
ggsave(plot_file3, p3, width = 10, height = 12, bg = "white", dpi = 300) # Taller for many participants
message(sprintf("Matrix Plot saved to %s", plot_file3))

# Export detailed CSV with metadata
out_csv_detailed <- file.path(DIR_WGS, "panaroo_selection_detailed.csv")
write_csv(df, out_csv_detailed)
message(sprintf("Detailed summary saved to %s", out_csv_detailed))

# Plot 4: Overview - How many individuals have X valid isolates?
# Group by participant and count selected Longcycler assemblies.
overview_df <- df %>%
    group_by(Participant_id) %>%
    summarise(
        n_valid = sum(Selected),
        valid_timepoints = paste(sort(Timepoint[Selected]), collapse = ",")
    ) %>%
    ungroup()

# Count how many participants have 0, 1, 2, ... valid isolates
# n_valid counts selected Longcycler episode assemblies.
overview_counts <- overview_df %>%
    dplyr::count(n_valid) %>%
    mutate(label = paste0(n, " participants"))

p4 <- ggplot(overview_counts, aes(x = factor(n_valid), y = n)) +
    geom_col(fill = "dodgerblue", width = 0.7) +
    geom_text(aes(label = n), vjust = -0.5) +
    labs(
        title = "Distribution of Selected Longcycler Episodes per Participant",
        subtitle = "Selected, QC-passing Longcycler episode assemblies",
        x = "Number of Selected Episodes per Participant",
        y = "Number of Participants"
    ) +
    theme_minimal() +
    theme(
        panel.grid.major.x = element_blank(),
        plot.title = element_text(face = "bold")
    )

plot_file4 <- file.path(DIR_PLOTS_WGS, "panaroo_selection_overview.png")
ggsave(plot_file4, p4, width = 8, height = 6, bg = "white", dpi = 300)
message(sprintf("Overview Plot saved to %s", plot_file4))

# ---- 4 · Timepoint-level selection summary ----------------------------------

message("\nSummarising Panaroo selection at Participant x Timepoint level ...")

# Define output directories for this section
OUTDIR <- DIR_WGS
PLOTDIR <- DIR_PLOTS_WGS

# 4.1 Use the Longcycler candidate/QC table with exact current metadata.
assembly_tp <- df %>%
    mutate(selected_for_panaroo = included_in_current_panaroo)

# 4.4 Collapse to Participant x Timepoint ----
# Here we decide which "timepoints" (participant x timepoint) survive QC.
timepoint_df <- assembly_tp %>%
    group_by(Participant_id, Timepoint) %>%
    summarise(
        n_isolates_total = n(),
        n_isolates_selected = sum(selected_for_panaroo, na.rm = TRUE),
        status = case_when(
            n_isolates_selected == 0L ~ "not included (no GFF-backed canonical assembly)",
            n_isolates_selected == n_isolates_total ~ "included (selected Longcycler)",
            TRUE ~ "included (selected Longcycler/GFF subset)"
        ),
        .groups = "drop"
    )

# 4.5 Save table for inspection in Excel/R ----
timepoint_csv <- file.path(OUTDIR, "panaroo_timepoint_selection_summary.csv")
write_csv(timepoint_df, timepoint_csv)
message("Wrote timepoint selection summary to: ", timepoint_csv)

# 4.6 Barplot: how many timepoints are kept vs dropped per Timepoint ----

# Make an explicit, stable ordering of timepoints (T0, T1, T2, ... then anything else)
tp_order <- timepoint_df %>%
    reframe(Timepoint = sort(unique(Timepoint))) %>%
    pull(Timepoint)

tp_bar <- timepoint_df %>%
    mutate(
        Timepoint = factor(Timepoint, levels = tp_order),
        status_simple = if_else(
            n_isolates_selected > 0L,
            "≥1 isolate kept for Panaroo",
            "0 isolates included (timepoint absent)"
        )
    ) %>%
    dplyr::count(Timepoint, status_simple, name = "n_timepoints")

tp_bar_colours <- c(
    "≥1 isolate kept for Panaroo" = "dodgerblue",
    "0 isolates included (timepoint absent)" = "firebrick"
)
assert_ruti_scale_levels(
    tp_bar$status_simple, tp_bar_colours,
    context = "Panaroo timepoint bar fill", allow_na = FALSE
)

p_tp_bar <- ggplot(tp_bar, aes(x = Timepoint, y = n_timepoints, fill = status_simple)) +
    geom_col(position = "stack", colour = "black", linewidth = 0.2) +
    labs(
        title = "Panaroo Selection Outcome by Timepoint",
        subtitle = "Each bar = study timepoint; height = number of participant x timepoint combinations",
        x = "Study timepoint (per participant)",
        y = "Number of participant x timepoint combinations",
        fill = "Panaroo outcome"
    ) +
    scale_fill_manual(values = tp_bar_colours) +
    theme_bw(base_size = 10) +
    theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "right"
    )

bar_path <- file.path(PLOTDIR, "panaroo_timepoint_selection_bar.png")
ggsave(bar_path, p_tp_bar, width = 160, height = 110, units = "mm", dpi = 300)
message("Saved barplot of timepoint selection to: ", bar_path)

# 4.7 Tile plot: Participant x Timepoint grid ----
# Each row is one participant; we suppress labels but make the axis meaning explicit.

tp_tile <- timepoint_df %>%
    mutate(
        Participant_id = factor(Participant_id),
        Timepoint = factor(Timepoint, levels = tp_order)
    )

tp_tile_colours <- c(
    "included (selected Longcycler)" = "dodgerblue",
    "included (selected Longcycler/GFF subset)" = "skyblue",
    "not included (no GFF-backed canonical assembly)" = "firebrick"
)
assert_ruti_scale_levels(
    tp_tile$status, tp_tile_colours,
    context = "Panaroo timepoint tile fill", allow_na = FALSE
)

p_tp_tile <- ggplot(tp_tile, aes(x = Timepoint, y = Participant_id, fill = status)) +
    geom_tile(colour = "grey80") +
    scale_y_discrete(guide = "none") +
    scale_fill_manual(values = tp_tile_colours) +
    labs(
        title = "Panaroo Selection Status Heatmap",
        subtitle = "Each tile = one participant x timepoint; colour shows Panaroo QC outcome",
        x = "Study timepoint",
        y = "Participants (one row per participant)",
        fill = "Panaroo QC status"
    ) +
    theme_bw(base_size = 8) +
    theme(
        axis.text.x = element_text(angle = 45, hjust = 1)
    )

tile_path <- file.path(PLOTDIR, "panaroo_timepoint_selection_tile.png")
ggsave(tile_path, p_tp_tile, width = 140, height = 180, units = "mm", dpi = 300)
message("Saved tile plot of timepoint selection to: ", tile_path)

# 4.8 Zoomed-out Sample Utilization Plot (Samples per Timepoint) ----
# User request: "all the timepoints vs the number of samples that we can utilize vs the ones we cannot utilize"

sample_util_df <- df %>%
    mutate(
        Timepoint = factor(Timepoint, levels = tp_order),
        Status = if_else(included_in_current_panaroo, "Included in current Panaroo input", "Not included in current Panaroo input")
    ) %>%
    dplyr::count(Timepoint, Status)

sample_util_colours <- c(
    "Included in current Panaroo input" = "dodgerblue",
    "Not included in current Panaroo input" = "grey50"
)
assert_ruti_scale_levels(
    sample_util_df$Status, sample_util_colours,
    context = "Panaroo sample-utilisation fill", allow_na = FALSE
)

p5 <- ggplot(sample_util_df, aes(x = Timepoint, y = n, fill = Status)) +
    geom_col(position = "stack", width = 0.7) +
    geom_text(aes(label = n), position = position_stack(vjust = 0.5), size = 3, color = "white") +
    scale_fill_manual(values = sample_util_colours) +
    labs(
        title = "Sample Utilization for Panaroo by Timepoint",
        subtitle = sprintf("All %d Longcycler candidate/QC rows; %d selected with exact GFFs", nrow(df), included_count),
        x = "Timepoint",
        y = "Number of Samples",
        fill = "Status"
    ) +
    theme_minimal() +
    theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.major.x = element_blank()
    )

plot_file5 <- file.path(PLOTDIR, "panaroo_selection_samples_per_timepoint.png")
ggsave(plot_file5, p5, width = 10, height = 6, bg = "white", dpi = 300)
message(sprintf("Sample Utilization Plot saved to %s", plot_file5))

# 5. QC Bias Analysis (Primary UTI Status)
# ------------------------------------------------------------------------------
# [STAT] Check if QC failure rate differs by primary UTI status
# (e.g., are UTI samples more likely to fail QC than Not_UTI?)

status_map_file <- file.path(DIR_CLINICAL, "status_map.csv")

if (file.exists(status_map_file)) {
    message("\nRunning QC Bias Analysis (QC Pass vs primary UTI status)...")
    status_map <- read_csv(status_map_file, show_col_types = FALSE) %>%
        prefer_primary_uti_status() %>%
        apply_manual_sample_curation(context = "panaroo_qc_attrition_context") %>%
        filter_primary_analysis() %>%
        mutate(
            Participant_id = as.character(Participant_id),
            tp_clean = if ("tp_lab" %in% names(.)) canon_tp(tp_lab) else canon_tp(Timepoint)
        ) %>%
        select(Participant_id, tp_clean, Infection_Status) %>%
        distinct()
    if (nrow(status_map) != 583L || n_distinct(status_map$Participant_id) != 166L) {
        stop("QC attrition context must be the labelled 583-episode, 166-resident clinical source cohort.")
    }

    status_dupes <- assert_unique_keys(
        status_map,
        c("Participant_id", "tp_clean"),
        context = "status_map for Panaroo QC bias",
        out_path = file.path(DIR_QC, "status_map_duplicate_episode_keys.csv")
    )

    qc_episode <- df %>%
        mutate(Participant_id = as.character(Participant_id), tp_clean = canon_tp(Timepoint)) %>%
        group_by(Participant_id, tp_clean) %>%
        summarise(
            qc_pass = any(QC_PASS, na.rm = TRUE),
            selected_canonical = any(selected_canonical, na.rm = TRUE),
            gff_available = any(gff_available, na.rm = TRUE),
            included_in_current_panaroo = any(included_in_current_panaroo, na.rm = TRUE),
            .groups = "drop"
        )

    qc_status <- if (nrow(status_dupes) == 0) {
        status_map %>%
            left_join(qc_episode, by = c("Participant_id", "tp_clean"), relationship = "one-to-one") %>%
            mutate(
                qc_pass = coalesce(qc_pass, FALSE),
                selected_canonical = coalesce(selected_canonical, FALSE),
                gff_available = coalesce(gff_available, FALSE),
                included_in_current_panaroo = coalesce(included_in_current_panaroo, FALSE)
            )
    } else {
        tibble()
    }

    if (nrow(qc_status) > 0) {
        tbl <- table(qc_status$qc_pass, qc_status$Infection_Status)

        message("QC Pass/Fail by primary UTI status:")
        print(tbl)

        bias_summary <- qc_status %>%
            dplyr::count(Infection_Status, qc_pass, selected_canonical, gff_available, included_in_current_panaroo, name = "n") %>%
            mutate(dataset_layer = "source_clinical_attrition_qc_context", .before = 1)
        write_csv(bias_summary, file.path(DIR_QC, "qc_selection_bias_by_status.csv"))

        if (all(dim(tbl) > 1)) {
            test_res <- fisher.test(tbl)
            message("\nFisher exact test for QC bias:")
            print(test_res)

            loss <- qc_status %>%
                group_by(Infection_Status) %>%
                summarise(n = n(), n_qc_fail = sum(!qc_pass), pct_qc_fail = n_qc_fail / n, .groups = "drop") %>%
                arrange(desc(pct_qc_fail))
            writeLines(
                c(
                    "QC selection bias by primary UTI status",
                    sprintf("Generated: %s", format(Sys.time())),
                    "Dataset layer: full clinical source retained only for attrition/QC context (583 episodes).",
                    "",
                    capture.output(print(tbl)),
                    "",
                    sprintf("Fisher exact p-value: %.5g", test_res$p.value),
                    sprintf("Largest QC loss: %s (%.1f%%, n=%d)",
                            loss$Infection_Status[1], 100 * loss$pct_qc_fail[1], loss$n_qc_fail[1]),
                    "Interpretation: sparse primary-status counts require Fisher/exact-style checks. Differences in QC pass by status may bias genomic comparisons."
                ),
                file.path(DIR_QC, "qc_selection_bias_report.txt")
            )
        }
    }
}

append_denominator_summary(
    df,
    "13_visualise_panaroo_selection.R",
    "panaroo_selection_visualisation",
    "assembly",
    qc_file,
    "Distinguishes QC PASS, canonical selected, GFF available, and included in current Panaroo input"
)

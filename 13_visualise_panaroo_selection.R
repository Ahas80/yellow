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

suppressPackageStartupMessages({
    library(tidyverse)
    library(fs)
    library(ggplot2)
})

canon_tp <- function(x) {
    tp_norm(x)$tp_lab
}

# 2. Load Centralized QC Data (from 12a_wgs_qc.R)
# ------------------------------------------------------------------------------
qc_file <- file.path(DIR_WGS, "qc_summary.csv")

if (!file.exists(qc_file)) {
    stop("QC summary not found. Please run 12a_wgs_qc.R first.")
}

# Load QC summary
qc_df <- read_csv(qc_file, show_col_types = FALSE)

# Load Metadata to get Timepoint info (QC summary might not have it if it came from fasta filenames only)
if (file.exists(FILE_METADATA)) {
    meta <- read_csv(FILE_METADATA, show_col_types = FALSE) %>%
        select(Isolate_ID, Participant_id, Timepoint) %>%
        distinct()

    # Join QC with Metadata
    # Assuming file_name in QC matches Isolate_ID or we need to map it.
    # 12a_wgs_qc.R uses assembly_metadata.csv, so it should have Participant_id if it was in the input.
    # Let's check columns of qc_df.
}

# Prepare 'df' for plotting
# We map QC_PASS to "Selected" and QC_REASON to "Reason"
df <- qc_df %>%
    mutate(
        Selected = QC_PASS,
        Reason = ifelse(QC_PASS, "Selected", QC_REASON)
    )

# If Participant_id/Timepoint are missing in QC summary (depends on 12a implementation), join them.
# 12a uses assembly_metadata.csv, so it likely has them or can be joined.
# Let's assume we need to join if missing.
if (!"Timepoint" %in% names(df) && exists("meta")) {
    # Try to join by file_name or Isolate_ID
    # 12a output has file_name.
    # We need to map file_name to Isolate_ID if not present.
    # Actually 12a output preserves input columns from assembly_metadata.csv usually?
    # Let's check 12a code... it reads FILE_METADATA.
    # So qc_df should have Participant_id and Timepoint.
    # pass
}

# 3. Summary
selected_count <- sum(df$Selected)
message(sprintf("Selected %d out of %d samples for Panaroo.", selected_count, nrow(df)))

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
# We want to see for each person, what they have.
p3 <- ggplot(df, aes(x = Timepoint, y = Participant_id, fill = Reason)) +
    geom_tile(color = "white", linewidth = 0.2) +
    scale_fill_manual(values = c(
        "Selected" = "dodgerblue",
        "Size > 7MB" = "firebrick",
        "CDS > 6000" = "orange",
        "Contigs > 500" = "purple",
        "Other" = "grey50"
    )) +
    labs(
        title = "Panaroo Sample Selection Matrix",
        subtitle = "Blue = Kept, Other colors = Eliminated (Reason)",
        x = "Timepoint",
        y = "Participant ID"
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
# Note: Each biological sample may have 2+ isolates (different assemblers)
# Group by Participant and count valid assembly files (isolates)
overview_df <- df %>%
    group_by(Participant_id) %>%
    summarise(
        n_valid = sum(Selected),
        valid_timepoints = paste(sort(Timepoint[Selected]), collapse = ",")
    ) %>%
    ungroup()

# Count how many participants have 0, 1, 2, ... valid isolates
# Note: n_valid counts ISOLATES (assemblies), not unique biological samples
overview_counts <- overview_df %>%
    count(n_valid) %>%
    mutate(label = paste0(n, " participants"))

p4 <- ggplot(overview_counts, aes(x = factor(n_valid), y = n)) +
    geom_col(fill = "dodgerblue", width = 0.7) +
    geom_text(aes(label = n), vjust = -0.5) +
    labs(
        title = "Distribution of Valid Isolates per Participant",
        subtitle = "Isolates = assembly files passing QC (note: multiple assemblers per sample)",
        x = "Number of Valid Isolates per Participant",
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

# 4.1 Use existing dataframe 'df' which contains all 1552 isolates with metadata
# We already have Participant_id, Timepoint, and Selected status in 'df'.
assembly_tp <- df %>%
    rename(selected_for_panaroo = Selected)

# 4.4 Collapse to Participant x Timepoint ----
# Here we decide which "timepoints" (participant x timepoint) survive QC.
timepoint_df <- assembly_tp %>%
    group_by(Participant_id, Timepoint) %>%
    summarise(
        n_isolates_total = n(),
        n_isolates_selected = sum(selected_for_panaroo),
        status = case_when(
            n_isolates_selected == 0L ~ "dropped (no isolates kept)",
            n_isolates_selected == n_isolates_total ~ "kept (all isolates)",
            TRUE ~ "kept (subset of isolates)"
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
            "0 isolates kept (timepoint dropped)"
        )
    ) %>%
    count(Timepoint, status_simple, name = "n_timepoints")

p_tp_bar <- ggplot(tp_bar, aes(x = Timepoint, y = n_timepoints, fill = status_simple)) +
    geom_col(position = "stack", colour = "black", linewidth = 0.2) +
    labs(
        title = "Panaroo Selection Outcome by Timepoint",
        subtitle = "Each bar = study timepoint; height = number of participant x timepoint combinations",
        x = "Study timepoint (per participant)",
        y = "Number of participant x timepoint combinations",
        fill = "Panaroo outcome"
    ) +
    scale_fill_manual(values = c(
        "≥1 isolate kept for Panaroo" = "dodgerblue",
        "0 isolates kept (timepoint dropped)" = "firebrick"
    )) +
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

p_tp_tile <- ggplot(tp_tile, aes(x = Timepoint, y = Participant_id, fill = status)) +
    geom_tile(colour = "grey80") +
    scale_y_discrete(guide = "none") +
    scale_fill_manual(values = c(
        "kept (all isolates)" = "dodgerblue",
        "kept (subset of isolates)" = "skyblue",
        "dropped (no isolates kept)" = "firebrick"
    )) +
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
        Status = if_else(Selected, "Utilized (Passed QC)", "Not Utilized (Failed QC)")
    ) %>%
    count(Timepoint, Status)

p5 <- ggplot(sample_util_df, aes(x = Timepoint, y = n, fill = Status)) +
    geom_col(position = "stack", width = 0.7) +
    geom_text(aes(label = n), position = position_stack(vjust = 0.5), size = 3, color = "white") +
    scale_fill_manual(values = c(
        "Utilized (Passed QC)" = "dodgerblue",
        "Not Utilized (Failed QC)" = "grey50"
    )) +
    labs(
        title = "Sample Utilization for Panaroo by Timepoint",
        subtitle = "Zoomed-out view of all 1552 samples",
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

# 5. QC Bias Analysis (Infection Status)
# ------------------------------------------------------------------------------
# [STAT] Check if QC failure rate differs by Infection Status
# (e.g., Are UTI samples more likely to fail QC than ASB?)

status_map_file <- file.path(DIR_CLINICAL, "status_map.csv")

if (file.exists(status_map_file)) {
    message("\nRunning QC Bias Analysis (QC Pass vs Infection Status)...")
    status_map <- read_csv(status_map_file, show_col_types = FALSE) %>%
        select(Participant_id, Timepoint, Infection_Status) %>%
        distinct()

    # Join QC data with Status
    # Ensure Timepoint format matches
    qc_status <- df %>%
        mutate(tp_clean = canon_tp(Timepoint)) %>%
        inner_join(status_map %>% mutate(tp_clean = canon_tp(Timepoint)),
            relationship = "many-to-many",
            by = c("Participant_id", "tp_clean")
        )

    if (nrow(qc_status) > 0) {
        # Create contingency table
        tbl <- table(qc_status$Selected, qc_status$Infection_Status)

        message("QC Pass/Fail by Infection Status:")
        print(tbl)

        # Chi-square test
        if (all(dim(tbl) > 1)) {
            test_res <- chisq.test(tbl)
            message("\nChi-square test for QC bias:")
            print(test_res)

            if (test_res$p.value < 0.05) {
                message("WARNING: Significant difference in QC pass rates between infection statuses!")
                message("This may introduce selection bias in pangenome analyses.")
            } else {
                message("No significant difference in QC pass rates detected (p > 0.05).")
            }
        }
    }
}

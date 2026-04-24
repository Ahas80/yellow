#!/usr/bin/env Rscript
# ==============================================================================
# UTIGA/scripts/poster_01_all_swimmer.R
# ------------------------------------------------------------------------------
# Generation of a swimmer plot for ALL individuals, including singles,
# using the exact poster-safe Uricult 0.5x plotting labels.
#
# Output:
#   - UTIGA/poster_figures/fig0_all_swimmer.png
#   - UTIGA/poster_figures/fig0_all_swimmer.pdf
# ==============================================================================

# Ensure paths correctly resolve relative to the standard rUTIs working directory
# We assume this script is run from `rUTIs/UTIGA/scripts/` or `rUTIs/`.
# Sourcing from central config handling
tryCatch({
    source(file.path("..", "..", "00_config.R"))
    source(file.path("..", "..", "R", "plot_helpers.R"))
}, error = function(e) {
    source("00_config.R")
    source("R/plot_helpers.R")
})

suppressPackageStartupMessages({
    library(dplyr)
    library(readr)
    library(ggplot2)
    library(forcats)
})

msg("Starting compilation of comprehensive all-individuals poster swimmer plot...")

DATA_DIR <- file.path(DIR_ROOT, "results", "clinical")
OUT_DIR  <- file.path(DIR_ROOT, "UTIGA", "poster_figures")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# 1. Load data
status_file <- file.path(DATA_DIR, "status_map_with_poster_tp.csv")
if (!file.exists(status_file)) {
    status_file <- file.path(DATA_DIR, "status_map.csv")
    warning("status_map_with_poster_tp.csv not found. Falling back to status_map.csv. You might be missing standard .5 labels")
}
df_all <- read_csv(status_file, show_col_types = FALSE)

# 2. Add plotting bounds & safety checks
df_swim <- df_all %>%
    mutate(
        Participant_id = factor(Participant_id),
        Status = factor(Infection_Status, levels = c("Negative", "ASB", "UTI")),
        Time_Order = if ("Plot_TP_Num_Poster" %in% names(.)) Plot_TP_Num_Poster else as.numeric(gsub("T", "", Timepoint)),
        Plot_Label = if ("Plot_TP_Label_Poster" %in% names(.)) Plot_TP_Label_Poster else Timepoint
    ) %>%
    filter(!is.na(Time_Order)) %>%
    # Ensure factor order
    mutate(Plot_Label = reorder(factor(Plot_Label), Time_Order))

msg(sprintf("Plotting %d valid episode/events across %d unique participants.", nrow(df_swim), length(unique(df_swim$Participant_id))))

# 3. Custom Shape logic (Differentiate Uricult visually from standard visits)
# We will use shape to distinguish Uricults from scheduled boundaries
df_swim <- df_swim %>%
    mutate(
        Visit_Type = ifelse(grepl("uricult|\\.5", Plot_Label, ignore.case=TRUE), "Unscheduled (Uricult)", "Scheduled Routine")
    )

# 4. Generate Plot
p <- ggplot(df_swim, aes(x = Plot_Label, y = fct_rev(Participant_id))) +
    # Reference line for longitudinal span per patient (Only applies if they have > 1 visit, handled silently)
    geom_line(aes(group = Participant_id), color = "gray85", linewidth = 0.5) +
    # Overlay the true observations
    geom_point(aes(color = Status, shape = Visit_Type), size = 3.5, alpha = 0.9) +
    
    # Customise aesthetics matching standard palette
    scale_color_manual(values = c("Negative" = "#43a2ca", "ASB" = "#fec44f", "UTI" = "#de2d26", "Other" = "grey50")) +
    scale_shape_manual(values = c("Scheduled Routine" = 16, "Unscheduled (Uricult)" = 18)) +
    
    theme_minimal(base_size = 14) +
    labs(
        title = "Longitudinal Cohort Routine and Uricult Trajectories",
        subtitle = "Chronological arrangement mapping all participants including single-point evaluations",
        x = "Timeline Intervals (Derived)",
        y = "Participant ID",
        color = "Clinical Status",
        shape = "Event Source"
    ) +
    theme(
        axis.text.x = element_text(angle = 45, hjust = 1, face="bold"),
        axis.text.y = element_text(size = 7),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        legend.position = "right"
    )

# 5. Output rendering
png_path <- file.path(OUT_DIR, "fig0_all_swimmer.png")
pdf_path <- file.path(OUT_DIR, "fig0_all_swimmer.pdf")

# We use an expanded height because plotting >100 individuals dynamically requires vertical density
tall_height <- max(10, length(unique(df_swim$Participant_id)) * 0.15)

ggsave(png_path, p, width = 11, height = tall_height, bg = "white", dpi = 300)
ggsave(pdf_path, p, width = 11, height = tall_height, bg = "white")

msg(sprintf("Saved successfully to %s", png_path))
msg("Script Complete.")

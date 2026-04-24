#!/usr/bin/env Rscript
# ==============================================================================
# poster_02_flow_diagram.R
# Purpose: Generate CONSORT-style flow diagram for the poster.
# Run from: UTIGA/scripts/
# Output:   UTIGA/poster_figures/fig1_cohort_flow.png
# ==============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(grid)
  library(glue)
})

# --- PATHS ---
DATA_DIR <- file.path(dirname(getwd()), "data")
OUT_DIR  <- file.path(dirname(getwd()), "poster_figures")
dir.create(OUT_DIR, showWarnings = FALSE)

status_df <- read_csv(file.path(DATA_DIR, "status_map.csv"), show_col_types = FALSE) %>%
  mutate(Participant_id = as.character(Participant_id))
mlst_df <- read_csv(file.path(DATA_DIR, "mlst_with_meta.csv"), show_col_types = FALSE) %>%
  mutate(Participant_id = as.character(Participant_id))

# --- Compute Counts ---
# Total unique participants in status_map
all_pids <- unique(status_df$Participant_id)
n_total_pids <- length(all_pids)
n_total_episodes <- nrow(status_df)

# Participants with >=2 timepoints
tp_counts <- status_df %>%
  filter(!is.na(Infection_Status)) %>%
  group_by(Participant_id) %>%
  summarise(n_tp = n_distinct(Timepoint), .groups = "drop")

pids_2tp <- tp_counts %>% filter(n_tp >= 2) %>% pull(Participant_id)
n_pids_2tp <- length(pids_2tp)
n_episodes_2tp <- status_df %>% filter(Participant_id %in% pids_2tp) %>% nrow()

# ASB/UTI episodes
asb_uti_2tp <- status_df %>%
  filter(Participant_id %in% pids_2tp, Infection_Status %in% c("ASB", "UTI"))
n_asb_uti <- nrow(asb_uti_2tp)

# With ST linked (excluding untypeable)
mlst_linked <- mlst_df %>%
  filter(Participant_id %in% pids_2tp, !is.na(ST))
n_with_st <- n_distinct(
  paste0(mlst_linked$Participant_id, "__", mlst_linked$Timepoint)
)

typeable <- mlst_linked %>% filter(ST != "-")
n_typeable <- n_distinct(
  paste0(typeable$Participant_id, "__", typeable$Timepoint)
)

n_untypeable_episodes <- n_with_st - n_typeable

# >=3 TP subset
pids_3tp <- tp_counts %>% filter(n_tp >= 3) %>% pull(Participant_id)
n_pids_3tp <- length(pids_3tp)

# --- Draw Flow Diagram ---
# Using ggplot2 with annotate() for a clean box-and-arrow diagram

box_w <- 6.5
box_h <- 1.4
arrow_color <- "grey40"
box_fill <- "#E8F0FE"
box_border <- "#3367D6"
highlight_fill <- "#FFF3E0"
highlight_border <- "#E65100"

flow_data <- data.frame(
  x = c(0, 0, 0, 0, 0, -8, 8),
  y = c(10, 8, 6, 4, 2, 4, 6),
  label = c(
    glue("YELLOW RoUTIne Cohort\n{n_total_pids} participants, {n_total_episodes} episodes"),
    glue("≥2 Timepoints\n{n_pids_2tp} participants, {n_episodes_2tp} episodes"),
    glue("ASB or UTI episodes\n{n_asb_uti} episodes"),
    glue("ST linked (WGS + MLST)\n{n_with_st} episodes"),
    glue("Typeable isolates (analysis set)\n{n_typeable} episodes"),
    glue("Untypeable ('-')\n{n_untypeable_episodes} episodes excluded"),
    glue("Excluded: <2 timepoints\n{n_total_pids - n_pids_2tp} participants")
  ),
  is_exclusion = c(FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE)
)

p <- ggplot() +
  # Main boxes
  annotate("rect",
    xmin = flow_data$x[1:5] - box_w/2, xmax = flow_data$x[1:5] + box_w/2,
    ymin = flow_data$y[1:5] - box_h/2, ymax = flow_data$y[1:5] + box_h/2,
    fill = box_fill, color = box_border, linewidth = 0.8
  ) +
  # Exclusion boxes
  annotate("rect",
    xmin = flow_data$x[6:7] - box_w/2, xmax = flow_data$x[6:7] + box_w/2,
    ymin = flow_data$y[6:7] - box_h/2, ymax = flow_data$y[6:7] + box_h/2,
    fill = highlight_fill, color = highlight_border, linewidth = 0.6, linetype = "dashed"
  ) +
  # Labels
  annotate("text", x = flow_data$x, y = flow_data$y, label = flow_data$label,
    size = 5.5, fontface = "plain", lineheight = 1.1
  ) +
  # Main arrows (vertical)
  annotate("segment",
    x = 0, xend = 0,
    y = c(10 - box_h/2, 8 - box_h/2, 6 - box_h/2, 4 - box_h/2),
    yend = c(8 + box_h/2, 6 + box_h/2, 4 + box_h/2, 2 + box_h/2),
    arrow = arrow(length = unit(0.2, "cm"), type = "closed"),
    color = arrow_color, linewidth = 0.8
  ) +
  # Side arrows to exclusion boxes
  annotate("segment",
    x = c(0, 0), xend = c(-8 + box_w/2, 8 - box_w/2),
    y = c(5, 7), yend = c(4, 6),
    arrow = arrow(length = unit(0.15, "cm"), type = "closed"),
    color = highlight_border, linewidth = 0.6, linetype = "dashed"
  ) +
  # Subtitle annotation
  annotate("text", x = 0, y = 0.5,
    label = glue("Sensitivity analysis: ≥3 timepoints sub-cohort = {n_pids_3tp} participants"),
    size = 4.5, fontface = "italic", color = "grey40"
  ) +
  coord_cartesian(xlim = c(-12, 12), ylim = c(-0.2, 11)) +
  theme_void() +
  labs(title = "Study Flow Diagram") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 24, face = "bold", margin = margin(b = 10)),
    plot.margin = margin(10, 10, 10, 10)
  )

ggsave(file.path(OUT_DIR, "fig1_cohort_flow.png"), p, width = 12, height = 8, dpi = 600, bg = "white")
cat("Saved: poster_figures/fig1_cohort_flow.png\n")

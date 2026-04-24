#!/usr/bin/env Rscript
# ==============================================================================
# poster_05_trajectories.R
# Purpose: Selected participant ST trajectory plot for the poster.
#          Shows 10-15 illustrative participants with ≥3 timepoints,
#          including stable, switching, and UTI cases.
# Run from: UTIGA/scripts/
# Output:   UTIGA/poster_figures/fig5_trajectories.png
# ==============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
  library(ggplot2)
  library(glue)
  library(forcats)
  library(stringr)
})

# --- PATHS ---
DATA_DIR <- file.path(dirname(getwd()), "data")
OUT_DIR  <- file.path(dirname(getwd()), "poster_figures")
dir.create(OUT_DIR, showWarnings = FALSE)

status_file <- file.path(DATA_DIR, "status_map_with_poster_tp.csv")
if (!file.exists(status_file)) status_file <- file.path(DATA_DIR, "status_map.csv")
status_df <- read_csv(status_file, show_col_types = FALSE) %>%
  mutate(Participant_id = as.character(Participant_id))
mlst_df <- read_csv(file.path(DATA_DIR, "mlst_with_meta.csv"), show_col_types = FALSE) %>%
  mutate(Participant_id = as.character(Participant_id))
class_inputs <- read_csv(file.path(DATA_DIR, "class_inputs_full.csv"), show_col_types = FALSE) %>%
  mutate(Participant_id = as.character(Participant_id)) %>%
  select(Participant_id, tp_lab, isolate_ID)

# --- Filter to ≥3-TP cohort and link STs ---
tp_counts <- status_df %>%
  filter(!is.na(Infection_Status)) %>%
  group_by(Participant_id) %>%
  summarise(n_tp = n_distinct(Timepoint), .groups = "drop")

pids_3tp <- tp_counts %>% filter(n_tp >= 3) %>% pull(Participant_id)

cohort <- status_df %>%
  filter(Participant_id %in% pids_3tp,
         Infection_Status %in% c("ASB", "UTI")) %>%
  mutate(tp_num = if("Plot_TP_Num_Poster" %in% names(.)) {
    Plot_TP_Num_Poster
  } else {
    case_when(
      grepl("T[0-9]+", Timepoint) ~ as.numeric(gsub("T", "", Timepoint)),
      Timepoint == "Uricult" ~ 5, TRUE ~ 999
    )
  }) %>%
  filter(!is.na(tp_num))

cohort_mlst <- cohort %>%
  left_join(class_inputs, by = c("Participant_id", "Timepoint" = "tp_lab"),
            relationship = "many-to-many") %>%
  left_join(mlst_df %>% select(Isolate_ID, ST),
            by = c("isolate_ID" = "Isolate_ID"),
            relationship = "many-to-many") %>%
  distinct(Participant_id, Timepoint, .keep_all = TRUE) %>%
  mutate(ST_display = ifelse(is.na(ST) | ST == "-", "Untypeable", paste0("ST", ST)))

# --- Classify participants ---
pid_summary <- cohort_mlst %>%
  filter(!is.na(ST), ST != "-") %>%
  group_by(Participant_id) %>%
  summarise(
    n_episodes = n(),
    n_sts = n_distinct(ST),
    has_uti = any(Infection_Status == "UTI"),
    is_stable = n_distinct(ST) == 1,
    .groups = "drop"
  ) %>%
  filter(n_episodes >= 2)

# Select illustrative cases (plotting ALL UTI residents for completeness)
uti_switchers   <- pid_summary %>% filter(has_uti, !is_stable)
uti_stable      <- pid_summary %>% filter(has_uti, is_stable)

# Optionally include a few ASB cases for contrast (or set n=0 to drop them)
asb_switchers   <- pid_summary %>% filter(!has_uti, !is_stable) %>% slice_head(n = 2)
asb_stable      <- pid_summary %>% filter(!has_uti, is_stable, n_episodes >= 3) %>% slice_head(n = 2)

selected_pids <- unique(c(
  uti_switchers$Participant_id,
  uti_stable$Participant_id,
  asb_switchers$Participant_id,
  asb_stable$Participant_id
))

cat(glue("Selected {length(selected_pids)} participants for trajectory plot\n"))
cat(glue("  UTI + switched: {nrow(uti_switchers)}\n"))
cat(glue("  UTI + stable:   {nrow(uti_stable)}\n"))
cat(glue("  ASB + switched: {nrow(asb_switchers)}\n"))
cat(glue("  ASB + stable:   {nrow(asb_stable)}\n"))

# --- Prepare plot data ---
plot_df <- cohort_mlst %>%
  filter(Participant_id %in% selected_pids) %>%
  arrange(Participant_id, tp_num)

# Clean timepoint labels
plot_df <- plot_df %>%
  mutate(
    Timepoint_label = ifelse(Timepoint == "Uricult", "UTI\nepisode", Timepoint),
    PID_label = paste0("P", Participant_id)
  )

# Add category annotation
plot_df <- plot_df %>%
  left_join(pid_summary %>% select(Participant_id, has_uti, is_stable),
            by = "Participant_id") %>%
  mutate(
    Category = case_when(
      has_uti & !is_stable ~ "UTI + ST switch",
      has_uti & is_stable  ~ "UTI + Stable ST",
      !has_uti & !is_stable ~ "ASB + ST switch",
      TRUE ~ "ASB + Stable ST"
    ),
    Category = factor(Category, levels = c("UTI + ST switch", "UTI + Stable ST",
                                            "ASB + ST switch", "ASB + Stable ST"))
  )

# Order participants: by category then PID
pid_order <- plot_df %>%
  distinct(PID_label, Category) %>%
  arrange(Category, PID_label) %>%
  pull(PID_label)

plot_df$PID_label <- factor(plot_df$PID_label, levels = rev(pid_order))

# --- Color palette for top STs ---
top_sts_all <- plot_df %>%
  filter(ST_display != "Untypeable") %>%
  count(ST_display, sort = TRUE) %>%
  pull(ST_display)

n_colors <- length(top_sts_all)
qual_palette <- c("#2196F3", "#FF5722", "#4CAF50", "#9C27B0", "#FFC107",
                  "#00BCD4", "#E91E63", "#795548", "#607D8B", "#FF9800",
                  "#3F51B5", "#CDDC39", "#009688", "#F44336", "#673AB7")

st_pal <- setNames(qual_palette[seq_len(min(n_colors, length(qual_palette)))],
                   top_sts_all[seq_len(min(n_colors, length(qual_palette)))])
st_pal["Untypeable"] <- "#BDBDBD"

# --- Plot ---
p <- ggplot(plot_df, aes(x = tp_num, y = PID_label)) +
  # Connect points per participant
  geom_line(aes(group = Participant_id), color = "grey75", linewidth = 1.5) +
  # Points: color = ST, shape = infection status
  geom_point(aes(color = ST_display, shape = Infection_Status), size = 8, stroke = 1.5) +
  scale_color_manual(values = st_pal, name = "Sequence Type") +
  scale_shape_manual(
    values = c("ASB" = 17, "UTI" = 16),  # triangle = ASB, circle = UTI
    name = "Infection Status"
  ) +
  scale_x_continuous(
    breaks = c(0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5),
    labels = c("T0", "T0.5", "T1", "T1.5", "T2", "T2.5", "T3", "T3.5", "T4", "T4.5", "UTI\nepisode"),
    minor_breaks = NULL
  ) +
  facet_grid(Category ~ ., scales = "free_y", space = "free_y",
             switch = "y") +
  labs(
    title = expression(paste("Longitudinal ST Trajectories in Selected ", italic("E. coli"), " Carriers")),
    subtitle = "≥3-timepoint sub-cohort · ▲ = ASB, ● = UTI episode",
    x = "Timepoint",
    y = NULL
  ) +
  theme_minimal(base_size = 18) +
  theme(
    plot.title = element_text(face = "bold", size = 22, hjust = 0.5),
    plot.subtitle = element_text(color = "grey40", size = 16, hjust = 0.5),
    strip.text.y.left = element_text(angle = 0, face = "bold", size = 12, hjust = 1),
    strip.placement = "outside",
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "right",
    legend.box = "vertical",
    axis.text.y = element_text(size = 12, face = "bold"),
    plot.margin = margin(10, 10, 10, 10)
  ) +
  guides(
    color = guide_legend(order = 1, override.aes = list(size = 6)),
    shape = guide_legend(order = 2, override.aes = list(size = 6))
  )

height <- max(10, 0.6 * length(selected_pids) + 4)
ggsave(file.path(OUT_DIR, "fig5_trajectories.png"), p,
       width = 15, height = height, dpi = 600, bg = "white")
cat(glue("Saved: poster_figures/fig5_trajectories.png ({length(selected_pids)} participants)\n"))

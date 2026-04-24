#!/usr/bin/env Rscript
# ==============================================================================
# poster_06_asb_vs_uti.R
# Purpose: Poster-quality ASB vs UTI ST comparison figure.
# Run from: UTIGA/scripts/
# Output:   UTIGA/poster_figures/fig4_st_by_status.png
# ==============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
  library(ggplot2)
  library(glue)
  library(forcats)
})

# --- PATHS ---
DATA_DIR <- file.path(dirname(getwd()), "data")
OUT_DIR  <- file.path(dirname(getwd()), "poster_figures")
dir.create(OUT_DIR, showWarnings = FALSE)

status_df <- read_csv(file.path(DATA_DIR, "status_map.csv"), show_col_types = FALSE) %>%
  mutate(Participant_id = as.character(Participant_id))
mlst_df <- read_csv(file.path(DATA_DIR, "mlst_with_meta.csv"), show_col_types = FALSE) %>%
  mutate(Participant_id = as.character(Participant_id))
class_inputs <- read_csv(file.path(DATA_DIR, "class_inputs_full.csv"), show_col_types = FALSE) %>%
  mutate(Participant_id = as.character(Participant_id)) %>%
  select(Participant_id, tp_lab, isolate_ID)

# --- Filter to >=2 TP, ASB/UTI, typeable ---
tp_counts <- status_df %>%
  filter(!is.na(Infection_Status)) %>%
  group_by(Participant_id) %>%
  summarise(n_tp = n_distinct(Timepoint), .groups = "drop")

pids_2tp <- tp_counts %>% filter(n_tp >= 2) %>% pull(Participant_id)

cohort <- status_df %>%
  filter(Participant_id %in% pids_2tp,
         Infection_Status %in% c("ASB", "UTI"),
         cfu_recorded_any == TRUE)

cohort_mlst <- cohort %>%
  left_join(class_inputs, by = c("Participant_id", "Timepoint" = "tp_lab"),
            relationship = "many-to-many") %>%
  left_join(mlst_df %>% select(Isolate_ID, ST),
            by = c("isolate_ID" = "Isolate_ID"),
            relationship = "many-to-many") %>%
  filter(!is.na(ST), ST != "-") %>%
  distinct(Participant_id, Timepoint, .keep_all = TRUE)

# --- Identify top STs across both groups ---
overall_top <- cohort_mlst %>%
  count(ST, sort = TRUE) %>%
  slice_head(n = 8) %>%
  pull(ST)

plot_df <- cohort_mlst %>%
  mutate(ST_group = ifelse(ST %in% overall_top, paste0("ST", ST), "Other")) %>%
  count(Infection_Status, ST_group) %>%
  group_by(Infection_Status) %>%
  mutate(pct = n / sum(n) * 100) %>%
  ungroup()

# Order by overall frequency
st_order <- plot_df %>%
  group_by(ST_group) %>%
  summarise(total = sum(n), .groups = "drop") %>%
  arrange(desc(total)) %>%
  pull(ST_group)

# "Other" always last
st_order <- c(setdiff(st_order, "Other"), "Other")
plot_df$ST_group <- factor(plot_df$ST_group, levels = rev(st_order))

# Counts per status
n_asb <- sum(plot_df$n[plot_df$Infection_Status == "ASB"])
n_uti <- sum(plot_df$n[plot_df$Infection_Status == "UTI"])

# --- Color palette ---
st_colors <- c(
  "ST43" = "#2196F3", "ST4"  = "#4CAF50", "ST6"  = "#FF5722",
  "ST1"  = "#9C27B0", "ST22" = "#FFC107", "ST2"  = "#FF9800",
  "ST87" = "#795548", "ST10" = "#E91E63", "ST33" = "#CDDC39",
  "ST3"  = "#00BCD4", "Other" = "#E0E0E0"
)

# --- Plot: grouped bar chart ---
p <- ggplot(plot_df, aes(x = ST_group, y = pct, fill = Infection_Status)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.65,
           color = "white", linewidth = 0.3) +
  geom_text(aes(label = paste0(round(pct, 0), "%")),
            position = position_dodge(width = 0.7), vjust = -0.4,
            size = 5.5, fontface = "bold") +
  scale_fill_manual(
    values = c("ASB" = "#42A5F5", "UTI" = "#EF5350"),
    name = "Episode Type",
    labels = c(glue("ASB (n={n_asb})"), glue("UTI (n={n_uti})"))
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15)),
                     labels = function(x) paste0(x, "%")) +
  coord_flip() +
  labs(
    title = expression(paste(italic("E. coli"), " ST Distribution: ASB vs UTI")),
    subtitle = "≥2-timepoint cohort · Typeable isolates only",
    x = NULL, y = "Proportion of Episodes (%)"
  ) +
  theme_minimal(base_size = 20) +
  theme(
    plot.title = element_text(face = "bold", size = 24, hjust = 0.5),
    plot.subtitle = element_text(color = "grey40", size = 16, hjust = 0.5),
    legend.position = "top",
    legend.text = element_text(size = 16),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.y = element_text(face = "bold", size = 16),
    plot.margin = margin(15, 15, 15, 15)
  )

ggsave(file.path(OUT_DIR, "fig4_st_by_status.png"), p,
       width = 12, height = 8, dpi = 600, bg = "white")
cat("Saved: poster_figures/fig4_st_by_status.png\n")

# --- Print summary stats for captions ---
cat("\n=== ASB vs UTI Summary ===\n")
cat(glue("ASB episodes (typeable): {n_asb}\n"))
cat(glue("UTI episodes (typeable): {n_uti}\n"))
uti_top <- cohort_mlst %>%
  filter(Infection_Status == "UTI") %>%
  count(ST, sort = TRUE)
cat("Top STs in UTI:\n")
print(uti_top)

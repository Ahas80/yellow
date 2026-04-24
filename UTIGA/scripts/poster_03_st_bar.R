#!/usr/bin/env Rscript
# ==============================================================================
# poster_03_st_bar.R
# Purpose: Top ST distribution bar chart — filtered to ≥2-TP cohort,
#          typeable only, with percentage labels.
# Run from: UTIGA/scripts/
# Output:   UTIGA/poster_figures/fig2_st_distribution.png
# ==============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
  library(ggplot2)
  library(glue)
  library(forcats)
  library(scales)
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

# --- Filter to >=2 TP cohort ---
tp_counts <- status_df %>%
  filter(!is.na(Infection_Status)) %>%
  group_by(Participant_id) %>%
  summarise(n_tp = n_distinct(Timepoint), .groups = "drop")

pids_2tp <- tp_counts %>% filter(n_tp >= 2) %>% pull(Participant_id)

cohort <- status_df %>%
  filter(Participant_id %in% pids_2tp,
         Infection_Status %in% c("ASB", "UTI"),
         cfu_recorded_any == TRUE)

# Link to MLST
cohort_linked <- cohort %>%
  left_join(class_inputs, by = c("Participant_id", "Timepoint" = "tp_lab"),
            relationship = "many-to-many") %>%
  left_join(mlst_df %>% select(Isolate_ID, ST),
            by = c("isolate_ID" = "Isolate_ID"),
            relationship = "many-to-many") %>%
  filter(!is.na(ST), ST != "-") %>%
  distinct(Participant_id, Timepoint, .keep_all = TRUE)

cat(glue("Typeable isolates in ≥2-TP cohort: {nrow(cohort_linked)}\n"))
cat(glue("Distinct STs: {n_distinct(cohort_linked$ST)}\n"))

# --- ST frequency table ---
st_freq <- cohort_linked %>%
  count(ST, sort = TRUE) %>%
  mutate(pct = n / sum(n) * 100,
         label = paste0(round(pct, 1), "%"))

# Top 12 STs
n_show <- min(12, nrow(st_freq))
top_sts <- st_freq %>% slice_head(n = n_show)
other_n <- sum(st_freq$n) - sum(top_sts$n)

if (other_n > 0) {
  top_sts <- bind_rows(
    top_sts,
    tibble(ST = "Other", n = other_n,
           pct = other_n / sum(st_freq$n) * 100,
           label = paste0(round(other_n / sum(st_freq$n) * 100, 1), "%"))
  )
}

# --- Color palette ---
# Use a curated palette for known global lineages
st_colors <- c(
  "43"  = "#2196F3", "4"   = "#4CAF50", "6"   = "#FF5722",
  "1"   = "#9C27B0", "3"   = "#00BCD4", "22"  = "#FFC107",
  "10"  = "#E91E63", "87"  = "#795548", "35"  = "#607D8B",
  "2"   = "#FF9800", "658" = "#3F51B5", "33"  = "#CDDC39",
  "Other" = "#BDBDBD"
)

top_sts$fill <- ifelse(top_sts$ST %in% names(st_colors),
                       st_colors[top_sts$ST], "#90A4AE")

# --- Plot ---
p <- ggplot(top_sts, aes(x = reorder(paste0("ST", ST), n), y = n)) +
  geom_col(aes(fill = ST), width = 0.7, show.legend = FALSE) +
  geom_text(aes(label = glue("{n} ({label})")),
            hjust = -0.08, size = 5.5, fontface = "plain") +
  scale_fill_manual(values = setNames(top_sts$fill, top_sts$ST)) +
  coord_flip(clip = "off") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.25))) +
  labs(
    title = expression(paste("Top ", italic("E. coli"), " Sequence Types")),
    subtitle = glue("≥2-timepoint cohort · {sum(st_freq$n)} typeable isolates · {n_distinct(cohort_linked$ST)} distinct STs"),
    x = NULL,
    y = "Number of Isolates"
  ) +
  theme_minimal(base_size = 20) +
  theme(
    plot.title = element_text(face = "bold", size = 26),
    plot.subtitle = element_text(color = "grey40", size = 16),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.y = element_text(face = "bold", size = 16),
    plot.margin = margin(15, 60, 15, 15)
  )

ggsave(file.path(OUT_DIR, "fig2_st_distribution.png"), p,
       width = 12, height = 8, dpi = 600, bg = "white")
cat("Saved: poster_figures/fig2_st_distribution.png\n")

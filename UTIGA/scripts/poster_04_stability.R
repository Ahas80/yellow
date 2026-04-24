#!/usr/bin/env Rscript
# ==============================================================================
# poster_04_stability.R
# Purpose: Within-host ST stability comparison figure.
#          Side-by-side for ≥2-TP and ≥3-TP cohorts.
#          Uses PAIR-LEVEL counting (matching abstract methodology).
# Run from: UTIGA/scripts/
# Output:   UTIGA/poster_figures/fig3_stability.png
# ==============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
  library(ggplot2)
  library(glue)
  library(scales)
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

# --- Helper: compute pair-level stability for a given TP threshold ---
compute_stability <- function(status_df, mlst_df, class_inputs, min_tps) {
  tp_counts <- status_df %>%
    filter(!is.na(Infection_Status)) %>%
    group_by(Participant_id) %>%
    summarise(n_tp = n_distinct(Timepoint), .groups = "drop")

  pids <- tp_counts %>% filter(n_tp >= min_tps) %>% pull(Participant_id)

  cohort <- status_df %>%
    filter(Participant_id %in% pids,
           Infection_Status %in% c("ASB", "UTI"),
           cfu_recorded_any == TRUE)

  cohort_mlst <- cohort %>%
    left_join(class_inputs, by = c("Participant_id", "Timepoint" = "tp_lab"),
              relationship = "many-to-many") %>%
    left_join(mlst_df %>% select(Isolate_ID, ST),
              by = c("isolate_ID" = "Isolate_ID"),
              relationship = "many-to-many") %>%
    filter(!is.na(ST), ST != "-") %>%
    distinct(Participant_id, Timepoint, .keep_all = TRUE) %>%
    mutate(tp_num = if("Plot_TP_Num_Poster" %in% names(.)) {
      Plot_TP_Num_Poster
    } else {
      case_when(
        grepl("T[0-9]+", Timepoint) ~ as.numeric(gsub("T", "", Timepoint)),
        Timepoint == "Uricult" ~ 99, TRUE ~ 999
      )
    }) %>%
    filter(!is.na(tp_num))

  pairs <- cohort_mlst %>%
    distinct(Participant_id, tp_num, ST) %>%
    group_by(Participant_id) %>%
    arrange(tp_num) %>%
    mutate(Next_ST = lead(ST)) %>%
    filter(!is.na(Next_ST)) %>%
    ungroup() %>%
    mutate(Same = ST == Next_ST)

  n_pairs <- nrow(pairs)
  n_same  <- sum(pairs$Same)
  n_diff  <- n_pairs - n_same

  tibble(
    Cohort = glue("≥{min_tps} Timepoints\n(n={length(pids)} participants)"),
    Same_ST = n_same,
    Different_ST = n_diff,
    Total_Pairs = n_pairs,
    Pct_Same = round(n_same / n_pairs * 100, 1)
  )
}

# --- Compute for both cohorts ---
stab_2 <- compute_stability(status_df, mlst_df, class_inputs, 2)
stab_3 <- compute_stability(status_df, mlst_df, class_inputs, 3)

stab_all <- bind_rows(stab_2, stab_3) %>%
  pivot_longer(cols = c(Same_ST, Different_ST),
               names_to = "Category", values_to = "Count") %>%
  mutate(
    Category = factor(Category, levels = c("Different_ST", "Same_ST"),
                      labels = c("ST switched", "Same ST maintained")),
    Pct = Count / Total_Pairs * 100
  )

cat("Stability Data:\n")
print(stab_all)

# --- Plot ---
p <- ggplot(stab_all, aes(x = Cohort, y = Count, fill = Category)) +
  geom_col(width = 0.55, color = "white", linewidth = 0.5) +
  geom_text(aes(label = glue("{Count}\n({round(Pct, 1)}%)")),
            position = position_stack(vjust = 0.5),
            size = 6.5, fontface = "bold", color = "white") +
  geom_text(data = bind_rows(stab_2, stab_3),
            aes(x = Cohort, y = Total_Pairs, fill = NULL,
                label = glue("{Pct_Same}% stable")),
            vjust = -0.5, size = 7, fontface = "bold", color = "#2E7D32") +
  scale_fill_manual(
    values = c("Same ST maintained" = "#2E7D32", "ST switched" = "#E65100"),
    name = NULL
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(
    title = expression(paste("Within-Host ", italic("E. coli"), " ST Stability")),
    subtitle = "Consecutive isolate pairs from ASB/UTI episodes (typeable only)",
    x = NULL,
    y = "Number of Consecutive Isolate Pairs"
  ) +
  theme_minimal(base_size = 20) +
  theme(
    plot.title = element_text(face = "bold", size = 26, hjust = 0.5),
    plot.subtitle = element_text(color = "grey40", size = 16, hjust = 0.5),
    legend.position = "top",
    legend.text = element_text(size = 18),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(size = 18, face = "bold"),
    plot.margin = margin(15, 15, 15, 15)
  )

ggsave(file.path(OUT_DIR, "fig3_stability.png"), p,
       width = 11, height = 8, dpi = 600, bg = "white")
cat("Saved: poster_figures/fig3_stability.png\n")

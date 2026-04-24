#!/usr/bin/env Rscript
# ==============================================================================
# 07_explore_MLST.R
# ==============================================================================
#
# GOAL:
#   Summarise and visualise MLST results: ST frequency distributions, top-20
#   ST bar chart, and merged metadata for downstream lineage analyses.
#   Produces mlst_with_meta.csv which is the authoritative MLST+metadata join
#   used by scripts 14, 17, 22, and 25.
#
# ------------------------------------------------------------------------------
# Role: [Descriptive] - Explore and summarize MLST data.
#
# Inputs:
#   - results/mlst/mlst_all.tsv
#   - assembly_metadata.csv
#
# Outputs:
#   - results/mlst/ST_frequencies.csv
#   - results/mlst/mlst_with_meta.csv
#   - plots/mlst/top20_STs.pdf
#   - plots/mlst/top20_STs.png
#
# Usage:
#   Rscript 07_explore_MLST.R
#
# Biological/Statistical purpose:
#   - Visualizes the distribution of Sequence Types (STs) in the cohort.
#   - Identifies dominant lineages (e.g., ST131, ST73).
# ==============================================================================

# 1. Load Configuration & Libraries
source("00_config.R")
suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(fs)
  library(stringr)
  library(tidyr)
  library(scales)
  library(ragg)
})

# 2. Configuration
# ------------------------------------------------------------------------------
FILE_MLST <- FILE_MLST_ALL
DIR_DEBUG <- file.path(DIR_LOGS, "debug")
ensure_dir(DIR_DEBUG)
ensure_dir(DIR_PLOTS_MLST)

# 3. Load Data
# ------------------------------------------------------------------------------
if (!file.exists(FILE_MLST)) stop("Missing ", FILE_MLST)
mlst <- read_tsv(FILE_MLST, show_col_types = FALSE)

# Ensure ST column
st_col <- names(mlst)[tolower(names(mlst)) == "st"]
if (!length(st_col)) stop("No ST column found in MLST results.")
mlst$ST <- as.character(mlst[[st_col[1]]])

# Load Metadata
if (file.exists(FILE_METADATA)) {
  meta <- read_csv(FILE_METADATA, show_col_types = FALSE)
} else {
  warning("Metadata not found, skipping joins.")
  meta <- NULL
}

# 4. Collapse to One Row per Isolate
# ------------------------------------------------------------------------------
# Prefer complete MLST -> more loci typed -> first row
mlst_in <- mlst %>%
  mutate(
    n_loci_typed  = coalesce(as.numeric(n_loci_typed), -Inf),
    mlst_complete = coalesce(as.logical(mlst_complete), FALSE)
  ) %>%
  group_by(Isolate_ID) %>%
  arrange(desc(mlst_complete), desc(n_loci_typed), .by_group = TRUE) %>%
  slice_head(n = 1) %>%
  ungroup() %>%
  select(-any_of(c("n_loci_typed", "mlst_complete")))

# 5. ST Frequencies
# ------------------------------------------------------------------------------
st_freq <- mlst_in %>%
  filter(!is.na(ST), ST != "") %>%
  count(ST, sort = TRUE) %>%
  mutate(pct = percent(n / sum(n)))

write_csv(st_freq, file.path(DIR_MLST, "ST_frequencies.csv"))

# 6. Plot Top 20 STs
# ------------------------------------------------------------------------------
top20 <- st_freq %>% slice_max(n, n = 20, with_ties = FALSE)

p <- ggplot(top20, aes(x = reorder(factor(ST), n), y = n)) +
  geom_col(fill = "steelblue") +
  coord_flip() +
  labs(
    x = "Sequence type (ST)", y = "Number of Isolates",
    title = "Top 20 Most Frequent E. coli Sequence Types"
  ) +
  theme_minimal(base_size = 12)

ggsave(file.path(DIR_PLOTS_MLST, "top20_STs.pdf"), p, width = 7, height = 5)
ggsave(file.path(DIR_PLOTS_MLST, "top20_STs.png"), p, width = 7, height = 5, dpi = 300)

# 7. Join Metadata (if available)
# ------------------------------------------------------------------------------
if (!is.null(meta)) {
  iso_col <- names(meta)[grepl("(?i)^isolate[_ ]?id$", names(meta))][1]

  if (!is.na(iso_col) && "Isolate_ID" %in% names(mlst_in)) {
    # Diagnostic: Check duplicates
    dup_check <- meta %>%
      group_by(.data[[iso_col]]) %>%
      count() %>%
      filter(n > 1)
    if (nrow(dup_check) > 0) {
      write_csv(dup_check, file.path(DIR_DEBUG, "meta_duplicates.csv"))
      warning("Metadata has duplicate Isolate_IDs. See results/debug/meta_duplicates.csv")
    }

    # Join
    mlst_with_meta <- mlst_in %>%
      left_join(meta, by = setNames(iso_col, "Isolate_ID"), suffix = c("", "_meta"))

    write_csv(mlst_with_meta, file.path(DIR_MLST, "mlst_with_meta.csv"))
    msg("Joined MLST with metadata.")
  }
}

msg("✓ MLST exploration complete.")

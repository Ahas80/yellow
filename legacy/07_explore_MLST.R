#!/usr/bin/env Rscript
# =============================================================
# 07_explore_mlst.R  –  quick descriptive stats & QC
# Requires the output of 06_MLST.R in results/mlst/
# =============================================================

## ---- 0 · libraries ---------------------------------------------------------
suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(fs)
  library(stringr)
  library(tidyr)
  library(scales)
})

## ---- 1 · file paths --------------------------------------------------------
# EDIT HERE if you moved files
mlst_file    <- "results/mlst/mlst_all.tsv"
meta_file    <- "assembly_metadata.csv"
plot_file    <- "results/mlst/top20_STs.pdf"
out_freq_csv <- "results/mlst/ST_frequencies.csv"

## ---- 2 · load data ---------------------------------------------------------
mlst <- read_tsv(mlst_file,  show_col_types = FALSE)
meta <- read_csv(meta_file,  show_col_types = FALSE)

## ---- 3 · ST frequency table -----------------------------------------------
st_freq <- mlst %>%
  count(ST, sort = TRUE) %>%
  mutate(pct = percent(n / sum(n)))

write_csv(st_freq, out_freq_csv)
message("✓ wrote ST frequency table → ", out_freq_csv)

## ---- 4 · bar-plot of top 20 STs -------------------------------------------
top20 <- st_freq %>% slice_max(n, n = 20, with_ties = FALSE)

ggplot(top20, aes(x = reorder(ST, n), y = n)) +
  geom_col() +
  coord_flip() +
  labs(x = "Sequence type (ST)", y = "Assemblies",
       title = "Top 20 E. coli STs – Yellow RoUTIne") +
  theme_classic(base_size = 12)

ggsave(plot_file, width = 7, height = 5)
message("✓ saved plot → ", plot_file)

## ---- 5 · join ST back onto clinical meta

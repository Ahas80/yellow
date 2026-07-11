#!/usr/bin/env Rscript

# generate_non_typable_report.R
# Extracts active MLST non-typable/missing isolates and generates comparison tables.

suppressPackageStartupMessages({
    library(tidyverse)
})
source("00_config.R")

cat("Loading MLST data...\n")
df_mlst <- read_csv(FILE_MLST_CANONICAL, show_col_types = FALSE)

# "Non-typable" STs are encoded as missing-like values after provider promotion.
cat("Filtering for non-typable STs...\n")
missing_tokens <- c("", "-", "?", "NA", "N/A", "nan", "missing", "unknown")
df_untypable <- df_mlst %>%
    filter(is.na(ST) | stringr::str_squish(stringr::str_to_upper(as.character(ST))) %in% stringr::str_to_upper(missing_tokens))

cat(glue::glue("Found {nrow(df_untypable)} non-typable isolates out of {nrow(df_mlst)} total.\n"))

# Minimal Table for the colleague
minimal_table <- df_untypable %>%
    select(any_of(c("Isolate_ID", "Participant_id", "Timepoint", "tp_lab", "ST", "ST_source", "ST_provider", "ST_local", "Collection_Date"))) %>%
    arrange(Participant_id, Timepoint)

# Extended Table with QC and Assembly Metrics
ALLELE_COLS <- c("dinb", "icda", "pabb", "polb", "putp", "trpa", "trpb", "uida")
actual_allele_cols <- intersect(names(df_untypable), ALLELE_COLS)
if (!"n_loci_typed" %in% names(df_untypable)) df_untypable$n_loci_typed <- NA_real_

extended_table <- df_untypable %>%
    mutate(n_loci_typed = if (length(actual_allele_cols) > 0) rowSums(!is.na(select(., all_of(actual_allele_cols)))) else n_loci_typed) %>%
    select(any_of(c(
        "Isolate_ID", "Participant_id", "Timepoint", "tp_lab", "ST", "ST_source",
        "ST_provider", "ST_local", "provider_PercGoodTargets", "Collection_Date",
        "has_new_allele", "ambiguous_call", "n_loci_typed", "num_contigs", "total_bases"
    ))) %>%
    arrange(Participant_id, Timepoint)

write_csv(minimal_table, "results/mlst/non_typable_isolates_minimal.csv")
write_csv(extended_table, "results/mlst/non_typable_isolates_extended.csv")

cat("Generated: results/mlst/non_typable_isolates_minimal.csv\n")
cat("Generated: results/mlst/non_typable_isolates_extended.csv\n")
cat("Done.\n")

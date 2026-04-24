#!/usr/bin/env Rscript

# generate_non_typable_report.R
# Extracts MLST non-typable ("-") isolates and generates comparison tables.

suppressPackageStartupMessages({
    library(tidyverse)
})

cat("Loading MLST data...\n")
df_mlst <- read_csv("results/mlst/mlst_with_meta.csv", show_col_types = FALSE)

# "Non-typable" STs are encoded as "-" (or literal NA)
cat("Filtering for non-typable STs...\n")
df_untypable <- df_mlst %>% filter(ST == "-" | is.na(ST))

cat(glue::glue("Found {nrow(df_untypable)} non-typable isolates out of {nrow(df_mlst)} total.\n"))

# Minimal Table for the colleague
minimal_table <- df_untypable %>%
    select(Isolate_ID, Participant_id, Timepoint, ST, Collection_Date) %>%
    arrange(Participant_id, Timepoint)

# Extended Table with QC and Assembly Metrics
ALLELE_COLS <- c("dinb", "icda", "pabb", "polb", "putp", "trpa", "trpb", "uida")
actual_allele_cols <- intersect(names(df_untypable), ALLELE_COLS)

extended_table <- df_untypable %>%
    mutate(n_loci_typed = rowSums(!is.na(select(., all_of(actual_allele_cols))))) %>%
    select(Isolate_ID, Participant_id, Timepoint, ST, Collection_Date,
           has_new_allele, ambiguous_call, n_loci_typed, num_contigs, total_bases) %>%
    arrange(Participant_id, Timepoint)

write_csv(minimal_table, "results/mlst/non_typable_isolates_minimal.csv")
write_csv(extended_table, "results/mlst/non_typable_isolates_extended.csv")

cat("Generated: results/mlst/non_typable_isolates_minimal.csv\n")
cat("Generated: results/mlst/non_typable_isolates_extended.csv\n")
cat("Done.\n")

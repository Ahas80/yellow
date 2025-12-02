source("00_config.R")
library(tidyverse)
library(fs)

dist_file <- file.path(DIR_WGS_CORE, "mash_dists.tsv")
pairs_file <- file.path(DIR_WGS_CORE, "strain_pairs.csv")

if (!file.exists(dist_file)) stop("mash_dists.tsv not found")

msg("Loading Mash distances...")
# Mash output: Ref, Query, Dist, P-val, Matching-Hashes
dists <- read_tsv(dist_file, col_names = c("A_path", "B_path", "dist", "pval", "hashes"), show_col_types = FALSE)

msg("Processing pairs...")
pairs_classified <- dists %>%
    mutate(
        A = path_file(A_path) %>% str_remove("\\.fasta$"),
        B = path_file(B_path) %>% str_remove("\\.fasta$")
    ) %>%
    filter(A < B) %>% # Unique pairs
    mutate(
        est_snps = dist * 5000000, # Approx 5Mb genome
        call = case_when(
            est_snps <= 20 ~ "Same",
            est_snps <= 1000 ~ "Related",
            TRUE ~ "Different"
        )
    ) %>%
    select(A, B, snps = est_snps, call)

write_csv(pairs_classified, pairs_file)
msg("Saved %d pairs to %s", nrow(pairs_classified), pairs_file)

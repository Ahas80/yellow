source("00_config.R")
library(tidyverse)

dist_file <- file.path(DIR_WGS_CORE, "snp_dists.tsv")
pairs_file <- file.path(DIR_WGS_CORE, "strain_pairs.csv")

if (!file.exists(dist_file)) stop("snp_dists.tsv not found")

msg("Loading distances...")
dists <- read_tsv(dist_file, show_col_types = FALSE)
colnames(dists)[1] <- "A"

msg("Processing pairs...")
pairs_long <- dists %>%
    pivot_longer(-A, names_to = "B", values_to = "snps") %>%
    filter(A < B)

pairs_classified <- pairs_long %>%
    mutate(
        call = case_when(
            snps <= 20 ~ "Same",
            snps <= 1000 ~ "Related",
            TRUE ~ "Different"
        )
    )

write_csv(pairs_classified, pairs_file)
msg("Saved %d pairs to %s", nrow(pairs_classified), pairs_file)

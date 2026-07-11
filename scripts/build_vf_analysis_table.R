#!/usr/bin/env Rscript
# ==============================================================================
# build_vf_analysis_table.R
# ==============================================================================
# Builds a clean single analysis-ready table and runs a data audit.
# Outputs: results/vf/vf_analysis_table.csv
# ==============================================================================

suppressPackageStartupMessages({
    library(dplyr)
    library(readr)
    library(tidyr)
    library(stringr)
})
source("00_config.R")
source("R/pipeline_qc_helpers.R")

cat("=== build_vf_analysis_table.R ===\n")

# ---- Load ----
vf_pa <- read_csv(FILE_VF_PA, show_col_types = FALSE) %>%
    mutate(Participant_id = as.character(Participant_id), tp_lab = as.character(tp_lab))
status <- read_csv(FILE_STATUS_MAP, show_col_types = FALSE) %>%
    prefer_primary_uti_status() %>%
    mutate(Participant_id = as.character(Participant_id))
gmap <- read_csv(file.path(DIR_VF, "gene_map.csv"), show_col_types = FALSE)

gene_cols <- canonical_vf_gene_cols(names(vf_pa))

# ---- Normalize timepoint keys ----
status <- status %>%
    mutate(tp_lab = case_when(
        str_detect(Timepoint, regex("uricult", TRUE)) ~ "Uricult",
        !is.na(suppressWarnings(as.integer(str_extract(Timepoint, "\\d+")))) ~
            paste0("T", str_extract(Timepoint, "\\d+")),
        TRUE ~ as.character(Timepoint)
    )) %>%
    group_by(Participant_id, tp_lab) %>%
    slice(1) %>%
    ungroup()

# ---- Join ----
df <- vf_pa %>%
    left_join(status %>% select(Participant_id, tp_lab, Infection_Status, Batch),
        by = c("Participant_id", "tp_lab")
    ) %>%
    mutate(Status = Infection_Status)

# ---- VF burden + categories ----
df$vf_count <- rowSums(df[, gene_cols, drop = FALSE], na.rm = TRUE)

cat_map <- split(gmap$Gene, gmap$Category)
for (cat_name in names(cat_map)) {
    col <- paste0("cat_", gsub("[/ ]", "_", cat_name))
    matching <- intersect(cat_map[[cat_name]], gene_cols)
    df[[col]] <- if (length(matching)) rowSums(df[, matching, drop = FALSE], na.rm = TRUE) else 0L
}

# Participant-level tp counts
tp_counts <- df %>%
    filter(!is.na(Status)) %>%
    group_by(Participant_id) %>%
    summarise(n_tp = n_distinct(tp_lab), .groups = "drop")
df <- df %>% left_join(tp_counts, by = "Participant_id")

# ---- iro operon score ----
iro_genes <- intersect(c("iroB", "iroC", "iroD", "iroE", "iroN"), gene_cols)
df$iro_score <- rowSums(df[, iro_genes, drop = FALSE], na.rm = TRUE)
df$iro_present <- df$iro_score > 0

# ---- Write ----
write_csv(df, "results/vf/vf_analysis_table.csv")

# ---- AUDIT ----
cat("\n=== DATA AUDIT ===\n\n")
cat("Episodes:", nrow(df), "\n")
cat("Participants:", n_distinct(df$Participant_id), "\n")
cat("Gene cols:", length(gene_cols), "\n")
cat("Status breakdown:\n")
print(table(df$Status, useNA = "ifany"))
cat("\nStatus × Timepoint:\n")
print(table(df$Status, df$tp_lab))
cat("\nParticipants by n_tp:\n")
print(table(tp_counts$n_tp))
cat("\nVF burden summary:\n")
print(summary(df$vf_count))
cat("\niro+ prevalence:\n")
print(table(df$Status, df$iro_present, dnn = c("Status", "iro+")))
cat("\nDone.\n")

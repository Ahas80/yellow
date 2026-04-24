#!/usr/bin/env Rscript
# ==============================================================================
# 25_vf_lineage_context.R
# ------------------------------------------------------------------------------
# Role: [VF Context] - Place VF burden and candidate genes in lineage context.
#
# Inputs:
#   - results/vf/vf_episode_dataset.csv
#   - results/vf/vf_gene_prevalence_asb_vs_uti.csv (optional, from step 23)
#
# Outputs:
#   - results/vf/vf_burden_by_st.csv
#   - results/vf/vf_status_within_st.csv
#   - results/vf/vf_candidate_genes_by_st.csv
#
# Why this matters:
#   - A strong apparent VF/status signal may simply reflect clonal structure.
#   - This script quantifies whether burden/candidate genes are dominated by ST.
# ==============================================================================

source("00_config.R")

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
})

msg <- function(...) message(sprintf(...))

in_episode <- file.path(DIR_RESULTS, "vf", "vf_episode_dataset.csv")
in_hits <- file.path(DIR_RESULTS, "vf", "vf_gene_prevalence_asb_vs_uti.csv")
if (!file.exists(in_episode)) stop("Run 22_vf_build_analysis_dataset.R first.")

episode <- read_csv(in_episode, show_col_types = FALSE)

st_col <- c("ST", "SequenceType")[c("ST", "SequenceType") %in% names(episode)][1]
if (is.na(st_col)) stop("No ST column found in vf_episode_dataset.csv")

episode <- episode %>% mutate(ST_resolved = as.character(.data[[st_col]]))

burden_by_st <- episode %>%
  filter(!is.na(ST_resolved), ST_resolved != "") %>%
  group_by(ST_resolved) %>%
  summarise(
    n = n(),
    participants = n_distinct(Participant_id),
    median_vf_burden = median(VF_Burden, na.rm = TRUE),
    mean_vf_burden = mean(VF_Burden, na.rm = TRUE),
    uti_n = sum(Infection_Status == "UTI", na.rm = TRUE),
    asb_n = sum(Infection_Status == "ASB", na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(n))
write_csv(burden_by_st, file.path(DIR_RESULTS, "vf", "vf_burden_by_st.csv"))

status_within_st <- episode %>%
  filter(!is.na(ST_resolved), Infection_Status %in% c("ASB", "UTI")) %>%
  count(ST_resolved, Infection_Status, name = "n") %>%
  tidyr::pivot_wider(names_from = Infection_Status, values_from = n, values_fill = 0) %>%
  mutate(total = ASB + UTI,
         uti_prop = ifelse(total == 0, NA_real_, UTI / total)) %>%
  arrange(desc(total), desc(uti_prop))
write_csv(status_within_st, file.path(DIR_RESULTS, "vf", "vf_status_within_st.csv"))

if (file.exists(in_hits)) {
  hits <- read_csv(in_hits, show_col_types = FALSE)
  top_genes <- hits %>%
    arrange(p_value) %>%
    slice_head(n = 25) %>%
    pull(gene)

  gene_cols <- intersect(top_genes, names(episode))
  if (length(gene_cols) > 0) {
    cand_by_st <- episode %>%
      filter(!is.na(ST_resolved), ST_resolved != "") %>%
      group_by(ST_resolved) %>%
      summarise(across(all_of(gene_cols), ~ mean(.x, na.rm = TRUE)), .groups = "drop")
    write_csv(cand_by_st, file.path(DIR_RESULTS, "vf", "vf_candidate_genes_by_st.csv"))
  }
}

msg("Done.")

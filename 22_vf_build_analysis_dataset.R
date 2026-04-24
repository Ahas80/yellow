#!/usr/bin/env Rscript
# ==============================================================================
# 22_vf_build_analysis_dataset.R
# ==============================================================================
#
# GOAL:
#   Build the SINGLE canonical VF analysis-ready episode-level dataset by
#   joining four upstream data sources:
#     1. VF presence/absence matrix  (from 02_gene_presence_analysis.R)
#     2. Clinical episode status     (from 00b_classify_episodes.R)
#     3. Gene → Category mapping     (from 04_gene_breakdown.R)
#     4. MLST Sequence Types         (from 06_MLST.R)
#
# WHY THIS SCRIPT EXISTS:
#   Before this script, every downstream VF analysis independently re-derived
#   the same join of VF data + clinical status.  Different scripts used
#   slightly different join logic, leading to subtle inconsistencies in
#   denominators and status assignments.  This script centralises that join
#   into a single, validated artifact.
#
#   All downstream VF scripts (23_, 24_, 25_) depend ONLY on the output
#   of this script.  They never re-join raw inputs themselves.
#
# INPUTS:
#   - results/vf/vf_pa_all.csv             Participant × timepoint gene matrix
#                                           (binary 0/1 per VF gene)
#   - results/clinical/status_map.csv       Clinical episode classification
#                                           (ASB / UTI / Negative per episode)
#   - results/vf/gene_map.csv              Gene-to-category mapping
#                                           (e.g., fimH → Adhesion/Fimbriae)
#   - results/mlst/mlst_with_meta.csv       MLST typing results (optional)
#                                           (Sequence Type per isolate)
#
# OUTPUTS:
#   - results/vf/vf_analysis_ready.csv      The canonical dataset.  Each row
#                                           is one episode (Participant × timepoint)
#                                           with gene P/A, status, ST, burden counts,
#                                           and category-level counts.
#   - results/vf/vf_dataset_diagnostics.txt Human-readable log of merge quality:
#                                           duplicates, unmatched rows, etc.
#
# KEY DESIGN DECISIONS:
#   1. ALL genes are retained, including those labelled "Unassigned" in gene_map.
#      Rationale: dropping ~50% of detected genes would undercount VF burden.
#      Category summaries still show "Unassigned" as its own category.
#   2. Union-based gene calling (present in EITHER Flye OR Longcycler assembly)
#      is inherited from upstream 02_gene_presence_analysis.R. We do NOT change
#      that logic here — consistency with the established pipeline is paramount.
#   3. No statistical testing or plotting belongs in this script.  This is
#      purely a data preparation step.
#
# SUPERSEDES:
#   The dataset-building section (lines 1–200) of compute_vf_abstract_stats.R
# ==============================================================================

source("00_config.R")
suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
  library(stringr)
})

msg("Starting 22_vf_build_analysis_dataset.R")

# ==============================================================================
# 1. LOAD ANCHOR FILES
# ==============================================================================
# Each anchor file is produced by a specific upstream script.  If any is
# missing, the pipeline has not been run in order.

# --- Anchor 1: VF Presence/Absence Matrix ---
#   Produced by 02_gene_presence_analysis.R
#   Structure: rows = participant × timepoint, columns = VF genes (0/1)
#   This uses the union of hits across assemblers: a gene is called "present"
#   if it was detected in EITHER the Flye or Longcycler assembly for that
#   participant at that timepoint.
if (!file.exists(FILE_VF_PA)) stop("Missing ", FILE_VF_PA, ". Run 02_gene_presence_analysis.R first.")
vf_pa <- read_csv(FILE_VF_PA, show_col_types = FALSE) %>%
  mutate(Participant_id = as.character(Participant_id),
         tp_lab         = as.character(tp_lab))

# Separate metadata columns from gene columns.
# The first two columns are always Participant_id and tp_lab.
# Everything else is a VF gene (binary).
meta_cols <- c("Participant_id", "tp_lab")
gene_cols <- setdiff(names(vf_pa), meta_cols)

msg("Loaded VF P/A matrix: %d rows × %d gene columns, %d participants",
    nrow(vf_pa), length(gene_cols), n_distinct(vf_pa$Participant_id))

# --- Anchor 2: Clinical Status Map ---
#   Produced by 00b_classify_episodes.R
#   Maps each participant × timepoint to an Infection_Status:
#     ASB      = Asymptomatic Bacteriuria (culture-positive, no symptoms)
#     UTI      = Urinary Tract Infection (culture-positive + symptoms)
#     Negative = Culture-negative episode
#   Also contains Batch (recruitment batch) which is useful as a covariate.
status_file <- file.path(DIR_CLINICAL, "status_map.csv")
if (!file.exists(status_file)) stop("Missing ", status_file, ". Run 00b_classify_episodes.R first.")
status_map <- read_csv(status_file, show_col_types = FALSE) %>%
  mutate(Participant_id = as.character(Participant_id))

# The VF matrix uses "tp_lab" (T0, T1, T2, Uricult) while the status map uses
# "Timepoint" (which may have various formats like "0", "1", "Uricult").
# This helper normalises both to the same format so the join works correctly.
normalize_tp <- function(x) {
  x <- as.character(x)
  is_uricult <- str_detect(x, regex("uricult", ignore_case = TRUE))
  tp_num     <- suppressWarnings(as.integer(str_extract(x, "\\d+")))
  case_when(
    is_uricult  ~ "Uricult",
    !is.na(tp_num) ~ paste0("T", tp_num),
    TRUE ~ x
  )
}

status_map <- status_map %>% mutate(tp_lab = normalize_tp(Timepoint))

msg("Loaded status_map: %d rows. Status breakdown:", nrow(status_map))
print(table(status_map$Infection_Status, useNA = "ifany"))

# --- Anchor 3: Gene Category Map ---
#   Produced by 04_gene_breakdown.R
#   Maps each VF gene name to a functional category (e.g., "Iron acquisition",
#   "Adhesion/Fimbriae", "Toxins", "Capsule/Surface", "Invasion/Evasion").
#   Genes not mapped to a known category get "Unassigned" — these are still
#   included in all analyses per our design decision.
gene_map_file <- file.path(DIR_VF, "gene_map.csv")
if (!file.exists(gene_map_file)) stop("Missing ", gene_map_file, ". Run 04_gene_breakdown.R first.")
gene_map <- read_csv(gene_map_file, show_col_types = FALSE) %>%
  mutate(Gene = as.character(Gene),
         Category = coalesce(as.character(Category), "Unassigned"))

msg("Loaded gene_map: %d genes across %d categories",
    nrow(gene_map), n_distinct(gene_map$Category))

# --- Anchor 4: MLST Sequence Types (Optional) ---
#   Produced by 06_MLST.R
#   Provides the Sequence Type (ST) for each isolate.  ST is critical for
#   the lineage confounding analysis in 25_vf_lineage_vf_interaction.R:
#   we need to know whether VF differences between ASB and UTI are genuine
#   or simply reflect different STs carrying different VF arsenals.
#   If this file is missing, the pipeline continues but the ST column is NA.
mlst_file <- file.path(DIR_MLST, "mlst_with_meta.csv")
mlst_available <- file.exists(mlst_file)
if (mlst_available) {
  mlst <- read_csv(mlst_file, show_col_types = FALSE) %>%
    mutate(Participant_id = as.character(Participant_id),
           tp_lab = normalize_tp(Timepoint),
           ST = as.character(ST)) %>%
    select(Participant_id, tp_lab, ST) %>%
    # Some participants may have multiple assemblies per timepoint (e.g., from
    # different assemblers).  Keep only one ST per participant × timepoint.
    distinct(Participant_id, tp_lab, .keep_all = TRUE)
  msg("Loaded MLST: %d rows", nrow(mlst))
} else {
  msg("WARNING: MLST file not found (%s). ST column will be NA.", mlst_file)
  mlst <- NULL
}

# ==============================================================================
# 2. MERGE DIAGNOSTICS
# ==============================================================================
# Before merging, we check for data quality issues that could silently corrupt
# the analysis.  The diagnostics log is written to a separate file so it can
# be reviewed without re-running the script.

diag_lines <- character()
log_diag <- function(...) {
  line <- sprintf(...)
  diag_lines <<- c(diag_lines, line)
  cat(line, "\n")
}

log_diag("=== MERGE DIAGNOSTICS ===")
log_diag("Timestamp: %s", format(Sys.time()))
log_diag("")

# CHECK 1: Duplicate join keys
#   If the same (Participant_id, tp_lab) appears more than once in either
#   table, the left_join will produce unintended row multiplication.
vf_dupes <- vf_pa %>% count(Participant_id, tp_lab) %>% filter(n > 1)
status_dupes <- status_map %>% count(Participant_id, tp_lab) %>% filter(n > 1)

log_diag("Duplicate keys in vf_pa: %d", nrow(vf_dupes))
log_diag("Duplicate keys in status_map: %d", nrow(status_dupes))

# If status_map has duplicates (e.g., from multiple Uricult events at the same
# timepoint label), deduplicate by keeping the first row per key.
if (nrow(status_dupes) > 0) {
  log_diag("  Deduplicating status_map (taking first row per key)")
  status_map <- status_map %>%
    group_by(Participant_id, tp_lab) %>%
    slice(1) %>%
    ungroup()
}

# CHECK 2: Anti-join diagnostics
#   Which VF rows have no matching clinical status? (Should be zero if
#   every sequenced isolate has been clinically classified.)
#   Which status rows have no VF data? (Expected: episodes where bacteria
#   were not sequenced, or culture-negative episodes.)
vf_not_in_status <- vf_pa %>% anti_join(status_map, by = c("Participant_id", "tp_lab"))
status_not_in_vf <- status_map %>% anti_join(vf_pa, by = c("Participant_id", "tp_lab"))

log_diag("VF rows NOT in status_map: %d", nrow(vf_not_in_status))
if (nrow(vf_not_in_status) > 0) {
  log_diag("  Unmatched participants: %s",
           paste(unique(vf_not_in_status$Participant_id), collapse = ", "))
}
log_diag("Status rows NOT in VF: %d (episodes without sequenced VF data)", nrow(status_not_in_vf))

# ==============================================================================
# 3. BUILD ANALYSIS-READY DATASET
# ==============================================================================
# This is the core of the script: merge all inputs into a single table where
# each row = one episode and columns include VF genes, clinical status, ST,
# total VF burden, and category-level VF counts.

msg("Building analysis-ready dataset...")

# STEP 3a: Join VF P/A with clinical status
#   Use left_join to keep all VF rows (we have sequence data for these).
#   Episodes without clinical status will get NA for Infection_Status —
#   this is expected and logged above.
vf_ready <- vf_pa %>%
  left_join(
    status_map %>% select(Participant_id, tp_lab, Infection_Status, Batch),
    by = c("Participant_id", "tp_lab")
  )

# STEP 3b: Add MLST Sequence Type if available
#   ST is needed by 25_vf_lineage_vf_interaction.R to check whether VF
#   differences are driven by lineage rather than clinical status.
if (!is.null(mlst)) {
  vf_ready <- vf_ready %>%
    left_join(mlst, by = c("Participant_id", "tp_lab"))
  log_diag("MLST joined: %d rows with ST, %d without",
           sum(!is.na(vf_ready$ST)), sum(is.na(vf_ready$ST)))
} else {
  vf_ready$ST <- NA_character_
}

# STEP 3c: Compute total VF burden per episode
#   This is the count of all VF genes detected (including Unassigned).
#   It is the primary dependent variable for cross-sectional comparisons.
vf_ready <- vf_ready %>%
  mutate(vf_count_total = rowSums(across(all_of(gene_cols)), na.rm = TRUE))

# STEP 3d: Compute category-level VF counts
#   For each functional category in gene_map (e.g., "Adhesion/Fimbriae",
#   "Iron acquisition"), count how many genes from that category are present
#   in each episode.  This enables category-level enrichment testing in 23_.
categories <- unique(gene_map$Category)
for (cat in categories) {
  cat_genes <- gene_map %>% filter(Category == cat) %>% pull(Gene)
  matching <- intersect(cat_genes, gene_cols)
  col_name <- paste0("cat_", gsub("[/ ]", "_", cat))
  if (length(matching) > 0) {
    vf_ready[[col_name]] <- rowSums(vf_ready[, matching, drop = FALSE], na.rm = TRUE)
  } else {
    vf_ready[[col_name]] <- 0L
  }
}

# STEP 3e: Handle genes in the P/A matrix that are NOT in gene_map
#   These are VF genes detected by Abricate but not yet assigned to a
#   functional category.  We count them separately so they are not lost.
genes_not_in_map <- setdiff(gene_cols, gene_map$Gene)
if (length(genes_not_in_map) > 0) {
  vf_ready$cat_Unassigned_matrix <- rowSums(
    vf_ready[, genes_not_in_map, drop = FALSE], na.rm = TRUE
  )
  log_diag("Genes in matrix but NOT in gene_map: %d (counted in cat_Unassigned_matrix)",
           length(genes_not_in_map))
}

# STEP 3f: Compute timepoint depth per participant
#   This counts how many distinct timepoints each participant contributed
#   (among episodes with clinical status).  Used by 23_ and 24_ to stratify
#   analyses by ≥2, ≥3, ≥4 timepoints — important because participants with
#   more timepoints provide more information about longitudinal stability.
tp_depth <- vf_ready %>%
  filter(!is.na(Infection_Status)) %>%
  group_by(Participant_id) %>%
  summarise(n_timepoints = n_distinct(tp_lab), .groups = "drop")

vf_ready <- vf_ready %>%
  left_join(tp_depth, by = "Participant_id")

# ==============================================================================
# 4. WRITE OUTPUTS
# ==============================================================================

out_file <- file.path(DIR_VF, "vf_analysis_ready.csv")
write_csv(vf_ready, out_file)

# Log a comprehensive summary of what was produced
log_diag("")
log_diag("=== OUTPUT SUMMARY ===")
log_diag("File: %s", out_file)
log_diag("Rows: %d", nrow(vf_ready))
log_diag("Participants: %d", n_distinct(vf_ready$Participant_id))
log_diag("Gene columns: %d", length(gene_cols))
log_diag("With Infection_Status: %d", sum(!is.na(vf_ready$Infection_Status)))
log_diag("Without Infection_Status: %d", sum(is.na(vf_ready$Infection_Status)))
log_diag("ASB: %d", sum(vf_ready$Infection_Status == "ASB", na.rm = TRUE))
log_diag("UTI: %d", sum(vf_ready$Infection_Status == "UTI", na.rm = TRUE))
log_diag("Negative: %d", sum(vf_ready$Infection_Status == "Negative", na.rm = TRUE))
log_diag("Mean vf_count_total: %.1f", mean(vf_ready$vf_count_total))
log_diag("Median vf_count_total: %.0f", median(vf_ready$vf_count_total))
if (mlst_available) {
  log_diag("With ST: %d, Without ST: %d, Distinct STs: %d",
           sum(!is.na(vf_ready$ST)), sum(is.na(vf_ready$ST)),
           n_distinct(vf_ready$ST, na.rm = TRUE))
}
log_diag("Timepoint depth range: %d–%d",
         min(vf_ready$n_timepoints, na.rm = TRUE),
         max(vf_ready$n_timepoints, na.rm = TRUE))

# Write diagnostics to a separate text file for review
diag_file <- file.path(DIR_VF, "vf_dataset_diagnostics.txt")
writeLines(diag_lines, diag_file)
msg("Diagnostics written to %s", diag_file)

msg("✓ 22_vf_build_analysis_dataset.R complete.")

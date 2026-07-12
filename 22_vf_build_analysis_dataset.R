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
#   the same join of VF data + primary UTI status. Different scripts used
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
#                                           (primary UTI / Not_UTI plus legacy status per episode)
#   - results/vf/gene_map.csv              Gene-to-category mapping
#                                           (e.g., fimH → Adhesion/Fimbriae)
#   - results/mlst/mlst_provider_preferred.csv  Active provider-preferred MLST typing results (optional)
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
#   2. Every row is derived from one selected, QC-passing Longcycler assembly.
#      Flye candidates are retained only in upstream QC/audit tables and cannot
#      enter this active analysis dataset as a fallback.
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
  library(optparse)
})

option_list <- list(
  make_option(c("--vf_pa"),
    type = "character", default = FILE_VF_PA,
    help = "VF presence/absence matrix to join [default: %default]"
  ),
  make_option(c("--out"),
    type = "character", default = FILE_VF_READY,
    help = "Analysis-ready output CSV [default: %default]"
  ),
  make_option(c("--diagnostics_out"),
    type = "character", default = NA_character_,
    help = "Diagnostics text output [default: <dirname(out)>/vf_dataset_diagnostics.txt]"
  ),
  make_option(c("--selection_file"),
    type = "character", default = FILE_ANALYSIS_ASSEMBLY_MANIFEST,
    help = "Assembly selection CSV used to filter MLST rows [default: %default]"
  ),
  make_option(c("--selection_column"),
    type = "character", default = "selected_canonical",
    help = "Logical column in --selection_file used to filter MLST rows [default: %default]"
  ),
  make_option(c("--qc_out_dir"),
    type = "character", default = NA_character_,
    help = "Directory for diagnostic bridge/unmatched tables [default: canonical DIR_QC or <dirname(out)>/qc]"
  )
)
opt <- parse_args(OptionParser(option_list = option_list))

out_file <- opt$out
out_dir <- dirname(out_file)
ensure_dir(out_dir)
is_primary_run <- identical(normalizePath(out_file, winslash = "/", mustWork = FALSE),
                            normalizePath(FILE_VF_READY, winslash = "/", mustWork = FALSE))
vf_diag_dir <- if (is_primary_run) DIR_VF else out_dir
qc_diag_dir <- if (!is.na(opt$qc_out_dir) && nzchar(opt$qc_out_dir)) opt$qc_out_dir else if (is_primary_run) DIR_QC else file.path(out_dir, "qc")
ensure_dir(vf_diag_dir)
ensure_dir(qc_diag_dir)
diag_file <- if (!is.na(opt$diagnostics_out) && nzchar(opt$diagnostics_out)) opt$diagnostics_out else file.path(vf_diag_dir, "vf_dataset_diagnostics.txt")
out_stem <- tools::file_path_sans_ext(basename(out_file))
stage_suffix <- if (is_primary_run) "" else paste0("_", stringr::str_replace(out_stem, "^vf_analysis_ready_?", ""))

msg("Starting 22_vf_build_analysis_dataset.R")
msg("VF P/A input: %s", opt$vf_pa)
msg("VF-ready output: %s", out_file)

# ==============================================================================
# 1. LOAD ANCHOR FILES
# ==============================================================================
# Each anchor file is produced by a specific upstream script.  If any is
# missing, the pipeline has not been run in order.

# --- Anchor 1: VF Presence/Absence Matrix ---
#   Produced by 02_gene_presence_analysis.R
#   Structure: rows = participant × timepoint, columns = VF genes (0/1)
#   This contains hits from one selected QC-passing Longcycler assembly per
#   participant-timepoint.
if (!file.exists(opt$vf_pa)) stop("Missing ", opt$vf_pa, ". Run 02_gene_presence_analysis.R first.")
selection <- load_analysis_assemblies(opt$selection_file, require_files = TRUE)
selection <- selection %>%
  mutate(
    Participant_id = as.character(Participant_id),
    tp_lab = normalise_timepoint_preserve_events(tp_lab),
    full_path = normalizePath(full_path, winslash = "/", mustWork = FALSE)
  )
vf_pa <- read_csv(opt$vf_pa, show_col_types = FALSE) %>%
  mutate(Participant_id = as.character(Participant_id),
         tp_lab = normalise_timepoint_preserve_events(tp_lab),
         Event_type = episode_event_type(tp_lab))

vf_keys <- vf_pa %>% distinct(Participant_id, tp_lab)
selection_keys <- selection %>% distinct(Participant_id, tp_lab)
if (nrow(anti_join(vf_keys, selection_keys, by = c("Participant_id", "tp_lab"))) > 0 ||
    nrow(anti_join(selection_keys, vf_keys, by = c("Participant_id", "tp_lab"))) > 0) {
  stop("VF presence/absence rows do not exactly match the Longcycler-only analysis manifest. Rerun script 02.")
}

# Separate metadata columns from gene columns.
# Common metadata columns are not VF genes.  Everything else is a binary VF
# gene column.
gene_cols <- canonical_vf_gene_cols(names(vf_pa))

msg("Loaded VF P/A matrix: %d rows × %d gene columns, %d participants",
    nrow(vf_pa), length(gene_cols), n_distinct(vf_pa$Participant_id))

# --- Anchor 2: Clinical Status Map ---
#   Produced by 00b_classify_episodes.R
#   Maps each participant × timepoint to the primary UTI_Status:
#     UTI      = catheter-aware S&S + culture support at >=10^3 CFU/mL
#     Not_UTI  = all non-UTI episodes, with Not_UTI_subgroup retained
#   Infection_Status is normalised to the primary UTI/Not_UTI status here.
#   Legacy ASB / UTI / Negative labels are retained only in explicit legacy
#   columns for audit and sensitivity checks.
#   Also contains Batch (recruitment batch) which is useful as a covariate.
status_file <- FILE_STATUS_MAP
if (!file.exists(status_file)) stop("Missing ", status_file, ". Run 00b_classify_episodes.R first.")
status_map <- read_csv(status_file, show_col_types = FALSE) %>%
  mutate(Participant_id = as.character(Participant_id)) %>%
  apply_manual_sample_curation(context = "vf_ready_status_map")

status_map <- status_map %>%
  mutate(
    tp_lab = if ("tp_lab" %in% names(.)) normalise_timepoint_preserve_events(tp_lab) else normalise_timepoint_preserve_events(Timepoint),
    Event_type = if ("Event_type" %in% names(.)) as.character(Event_type) else episode_event_type(tp_lab),
    Collection_Date = if ("Collection_Date" %in% names(.)) as.character(Collection_Date) else NA_character_
  )
if (!"UTI_Status" %in% names(status_map)) {
  stop("status_map lacks UTI_Status. Rerun 00b_classify_episodes.R; refusing legacy ASB/UTI/Negative fallback for primary VF dataset.")
}
if (!"UTI_binary" %in% names(status_map)) status_map$UTI_binary <- as.integer(status_map$UTI_Status == "UTI")
if (!"Infection_Status_legacy" %in% names(status_map) && "Infection_Status" %in% names(status_map)) {
  status_map$Infection_Status_legacy <- status_map$Infection_Status
}
if (!"Infection_Status_old" %in% names(status_map) && "Infection_Status_legacy" %in% names(status_map)) {
  status_map$Infection_Status_old <- status_map$Infection_Status_legacy
}
status_map <- prefer_primary_uti_status(status_map, allow_legacy_fallback = FALSE)
if (!"Episode_ID" %in% names(status_map)) {
  status_map$Episode_ID <- build_episode_id(status_map, timepoint_col = "tp_lab",
                                            event_col = "Event_type",
                                            date_col = "Collection_Date")
}

msg("Loaded status_map: %d rows. Primary status breakdown:", nrow(status_map))
print(table(status_map$UTI_Status, useNA = "ifany"))
msg("Loaded status_map primary-included breakdown:")
print(table(filter_primary_analysis(status_map)$UTI_Status, useNA = "ifany"))
msg("Legacy status breakdown:")
print(table(status_map$Infection_Status_legacy, useNA = "ifany"))

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
if (!"Subcategory" %in% names(gene_map)) gene_map$Subcategory <- NA_character_

msg("Loaded gene_map: %d genes across %d categories",
    nrow(gene_map), n_distinct(gene_map$Category))

# --- Anchor 4: MLST Sequence Types (Optional) ---
#   Produced by 06_MLST.R
#   Provides the Sequence Type (ST) for each isolate.  ST is critical for
#   the lineage confounding analysis in 25_vf_lineage_vf_interaction.R:
#   we need to know whether VF differences between UTI and Not_UTI are genuine
#   or simply reflect different STs carrying different VF arsenals.
#   If this file is missing, the pipeline continues but the ST column is NA.
mlst_file <- FILE_MLST_CANONICAL
mlst_available <- file.exists(mlst_file)
if (mlst_available) {
  mlst <- read_csv(mlst_file, show_col_types = FALSE) %>%
    mutate(Participant_id = as.character(Participant_id),
           tp_lab = if ("tp_lab" %in% names(.)) normalise_timepoint_preserve_events(tp_lab) else normalise_timepoint_preserve_events(Timepoint),
           ST = as.character(ST)) %>%
    select(any_of(c(
      "Participant_id", "tp_lab", "ST", "ST_source", "ST_provider", "ST_local",
      "provider_PercGoodTargets", "provider_file", "provider_batch_match",
      "provider_assembler", "full_path", "file_name", "assembler", "Assembler"
    )))

  if (!"full_path" %in% names(mlst)) stop(mlst_file, " lacks full_path; cannot verify Longcycler provenance.")
  canonical_paths <- selection$full_path
  mlst <- mlst %>%
    mutate(full_path = normalizePath(full_path, winslash = "/", mustWork = FALSE)) %>%
    filter(full_path %in% canonical_paths)
  if ("ST_source" %in% names(mlst) && "provider_assembler" %in% names(mlst)) {
    bad_provider <- mlst$ST_source == "provider_qc95" &
      (is.na(mlst$provider_assembler) | tolower(mlst$provider_assembler) != ANALYSIS_ASSEMBLER)
    if (any(bad_provider, na.rm = TRUE)) {
      stop("Active MLST contains provider calls that are not Longcycler-derived. Rerun 06_MLST.R.")
    }
  }

  mlst_conflicts <- mlst %>%
    group_by(Participant_id, tp_lab) %>%
    summarise(n_ST = n_distinct(ST[!is.na(ST)]), ST_values = paste(sort(unique(na.omit(ST))), collapse = ";"), .groups = "drop") %>%
    filter(n_ST > 1)
  if (nrow(mlst_conflicts) > 0) {
    write_csv(mlst_conflicts, file.path(DIR_QC, "mlst_duplicate_participant_timepoint_st_conflicts.csv"))
    msg("WARNING: %d participant-timepoints have conflicting ST calls; ST set to NA for those keys.", nrow(mlst_conflicts))
  }
  mlst <- mlst %>%
    mutate(
      ST_source = if ("ST_source" %in% names(.)) as.character(ST_source) else NA_character_,
      ST_provider = if ("ST_provider" %in% names(.)) as.character(ST_provider) else NA_character_,
      ST_local = if ("ST_local" %in% names(.)) as.character(ST_local) else NA_character_,
      provider_PercGoodTargets = if ("provider_PercGoodTargets" %in% names(.)) as.numeric(provider_PercGoodTargets) else NA_real_,
      provider_file = if ("provider_file" %in% names(.)) as.character(provider_file) else NA_character_,
      provider_batch_match = if ("provider_batch_match" %in% names(.)) as.character(provider_batch_match) else NA_character_,
      provider_assembler = if ("provider_assembler" %in% names(.)) as.character(provider_assembler) else NA_character_
    ) %>%
    group_by(Participant_id, tp_lab) %>%
    summarise(
      ST = if (n_distinct(ST[!is.na(ST)]) == 1) first(na.omit(ST)) else NA_character_,
      ST_source = { vals <- na.omit(ST_source); if (length(vals)) vals[1] else NA_character_ },
      ST_provider = { vals <- na.omit(ST_provider); if (length(vals)) vals[1] else NA_character_ },
      ST_local = { vals <- na.omit(ST_local); if (length(vals)) vals[1] else NA_character_ },
      provider_PercGoodTargets = { vals <- na.omit(provider_PercGoodTargets); if (length(vals)) vals[1] else NA_real_ },
      provider_file = { vals <- na.omit(provider_file); if (length(vals)) vals[1] else NA_character_ },
      provider_batch_match = { vals <- na.omit(provider_batch_match); if (length(vals)) vals[1] else NA_character_ },
      provider_assembler = { vals <- na.omit(provider_assembler); if (length(vals)) vals[1] else NA_character_ },
      .groups = "drop"
    )
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
vf_dupes <- vf_pa %>% dplyr::count(Participant_id, tp_lab) %>% filter(n > 1)
status_dupes <- status_map %>% dplyr::count(Participant_id, tp_lab) %>% filter(n > 1)

log_diag("Duplicate keys in vf_pa: %d", nrow(vf_dupes))
log_diag("Duplicate keys in status_map: %d", nrow(status_dupes))

if (nrow(status_dupes) > 0) {
  dup_file <- file.path(DIR_QC, "status_map_duplicate_episode_keys.csv")
  status_map %>%
    semi_join(status_dupes %>% select(Participant_id, tp_lab), by = c("Participant_id", "tp_lab")) %>%
    write_csv(dup_file)
  log_diag("  RED: duplicate status_map keys were written to %s", dup_file)
  stop("status_map has duplicate Participant_id + tp_lab keys. This script will not silently deduplicate; repair Episode_ID/timepoint keys upstream.")
}

# CHECK 2: Anti-join diagnostics
#   Which VF rows have no matching primary UTI status? (Should be zero if
#   every sequenced isolate has been clinically classified.)
#   Which status rows have no VF data? (Expected: episodes where bacteria
#   were not sequenced, or culture-negative episodes.)
vf_not_in_status <- vf_pa %>% anti_join(status_map, by = c("Participant_id", "tp_lab"))
status_not_in_vf <- filter_primary_analysis(status_map) %>% anti_join(vf_pa, by = c("Participant_id", "tp_lab"))
write_csv(vf_not_in_status, file.path(vf_diag_dir, "vf_without_status_rows.csv"))
write_csv(status_not_in_vf, file.path(vf_diag_dir, "status_without_vf_rows.csv"))
write_csv(
  bind_rows(
    vf_not_in_status %>% mutate(unmatched_source = "vf_without_status"),
    status_not_in_vf %>% mutate(unmatched_source = "status_without_vf")
  ),
  file.path(vf_diag_dir, "vf_status_unmatched_rows.csv")
)

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
# each row = one episode and columns include VF genes, primary UTI status, ST,
# total VF burden, and category-level VF counts.

msg("Building analysis-ready dataset...")

# ==============================================================================
# STEP 3a-pre: Bridge Uricult clinical episodes to their UTI-N WGS rows
# ==============================================================================
# Clinical Uricult episodes are classified with tp_lab = "Uricult" in
# status_map.  Their WGS data, however, appears in vf_pa under tp_lab =
# "UTI-1", "UTI-2", etc. (the overview spreadsheet labels them as UTI-N
# meetmomenten).  Without this bridge the left_join below silently drops
# every Uricult UTI.
#
# POLICY (primary analysis):
#   One representative UTI-N VF row per clinical Uricult episode.
#   Selection priority:
#     1. Date-matched candidates  →  lowest UTI-N number.
#     2. Participant-only match   →  lowest UTI-N number.
#   All alternative (non-selected) UTI-N rows are written to a sensitivity
#   audit table so they can be used in supplementary analyses.
#   A single clinical Uricult episode must NOT inflate multiple rows in the
#   primary UTI-vs-Not_UTI comparison.
# ==============================================================================

status_join_cols <- c("Participant_id", "tp_lab", "Episode_ID", "Collection_Date",
                      "UTI_Status", "UTI_binary", "Not_UTI_subgroup",
                      "Infection_Status", "Infection_Status_legacy",
                      "Infection_Status_old", "UTI_definition_version",
                      "UTI_classification_confidence", "UTI_classification_reason",
                      "Batch", "Status_Confidence_epi", "Sx_source_epi",
                      "UTI_Label", "Urine_collection_method",
                      "urine_collection_method_raw", "urine_collection_method_norm",
                      "catheter_rule", "symptom_compatible_uti", "symptom_rule_met",
                      "local_urinary_symptom_any", "systemic_symptom_any",
                      "flankpain_present", "dysuria_present", "urgency_present",
                      "frequency_present", "incontinence_present", "pus_present",
                      "fever_present", "rigors_present", "delirium_present",
                      "suprapubic_pain_present", "other_sxs_present",
                      "cfu_raw", "cfu_raw_parsed", "cfu_ge_1e3", "cfu_ge_1e4",
                      "cfu_ge_1e5", "culture_supports_uti",
                      "cfu_threshold_used_for_uti", "cfu_threshold_source",
                      "beoord_cat",
                      "analysis_include_primary", "analysis_exclusion_reason",
                      "duplicate_role", "duplicate_of_participant_id",
                      "duplicate_of_tp_lab", "allow_secondary_duplicate_qc",
                      "duplicate_use_note", "genomics_expected_include",
                      "genomics_exclusion_reason", "manual_curation_applied",
                      "manual_curation_note", "manual_curation_source")
status_for_join <- status_map %>% select(any_of(status_join_cols))

# --- Uricult clinical rows (exclude unsafe/unlinked participants) ---
uricult_clinical <- status_for_join %>%
  filter(tp_lab == "Uricult") %>%
  filter(analysis_include_primary %in% TRUE) %>%
  filter(!Participant_id %in% c("Still to be linked", "Niet te koppelen", "UNKNOWN"))

# --- UTI-N rows in vf_pa that have NO direct clinical match ---
vf_uti_n_unmatched <- vf_pa %>%
  filter(str_detect(tp_lab, "^UTI-\\d+$")) %>%
  anti_join(status_for_join, by = c("Participant_id", "tp_lab")) %>%
  select(Participant_id, tp_lab, Episode_ID_wgs = Episode_ID) %>%
  mutate(
    uti_n_num  = as.integer(str_extract(tp_lab, "\\d+$")),
    wgs_date_raw = vapply(
      strsplit(Episode_ID_wgs, "__", fixed = TRUE),
      function(x) if (length(x) >= 5) x[5] else NA_character_,
      character(1)
    ),
    wgs_date = as.Date(gsub("_", "/", wgs_date_raw), format = "%d/%m/%Y")
  )

if (nrow(uricult_clinical) > 0 && nrow(vf_uti_n_unmatched) > 0) {

  # Parse semicolon-separated clinical Collection_Dates into a list column
  uricult_with_dates <- uricult_clinical %>%
    mutate(
      .row = row_number(),
      clinical_dates = lapply(strsplit(coalesce(Collection_Date, ""), ";"), function(d) {
        as.Date(trimws(d), format = "%d/%m/%Y")
      })
    )

  # Candidate bridge: every Uricult × unmatched-UTI-N pair for same participant
  bridge_candidates <- uricult_with_dates %>%
    inner_join(vf_uti_n_unmatched, by = "Participant_id",
               relationship = "many-to-many") %>%
    mutate(
      date_match = mapply(function(cd, wd) {
        if (is.na(wd) || length(cd) == 0 || all(is.na(cd))) return(FALSE)
        wd %in% cd
      }, clinical_dates, wgs_date),
      match_basis = if_else(date_match, "participant_and_date", "participant_only")
    )

  # --- Full bridge audit (every candidate, selected or not) ---
  bridge_audit <- bridge_candidates %>%
    arrange(Participant_id, Episode_ID, uti_n_num) %>%
    transmute(
      Participant_id_clinical  = Participant_id,
      Episode_ID_clinical      = Episode_ID,
      tp_lab_clinical          = tp_lab.x,
      UTI_Status,
      UTI_binary,
      Not_UTI_subgroup,
      Infection_Status,
      Infection_Status_legacy,
      UTI_Label,
      Collection_Date_clinical = Collection_Date,
      mapped_tp_lab            = tp_lab.y,
      mapped_uti_n_num         = uti_n_num,
      Episode_ID_wgs,
      wgs_date                 = as.character(wgs_date),
      date_match,
      match_basis
    )

  # --- Select one representative per clinical Uricult episode ---
  bridge_selected <- bridge_candidates %>%
    group_by(Participant_id, Episode_ID) %>%
    arrange(desc(date_match), uti_n_num) %>%
    slice_head(n = 1) %>%
    ungroup() %>%
    mutate(
      selection_reason = if_else(
        date_match,
        "selected: date-matched, lowest UTI-N number",
        "selected: participant-only match, lowest UTI-N number"
      )
    )

  # Mark selected/not-selected in audit
  sel_key <- paste0(bridge_selected$Participant_id, "||",
                    bridge_selected$Episode_ID, "||",
                    bridge_selected$tp_lab.y)
  bridge_audit <- bridge_audit %>%
    mutate(
      .key = paste0(Participant_id_clinical, "||",
                    Episode_ID_clinical, "||", mapped_tp_lab),
      selected         = .key %in% sel_key,
      selection_reason = if_else(
        selected,
        bridge_selected$selection_reason[match(.key, sel_key)],
        "not selected: alternative UTI-N for same clinical episode"
      )
    ) %>%
    select(-.key)

  write_csv(bridge_audit, file.path(qc_diag_dir, "uricult_bridge_audit.csv"))

  # Sensitivity table: non-selected alternatives
  bridge_sensitivity <- bridge_audit %>%
    filter(!selected) %>%
    mutate(reason_excluded =
      "Alternative UTI-N for multi-isolate Uricult episode; available for sensitivity analysis")
  write_csv(bridge_sensitivity,
            file.path(qc_diag_dir, "uricult_bridge_sensitivity_alternatives.csv"))

  msg("Uricult bridge audit: %d candidates, %d selected (%d date-matched, %d participant-only).",
      nrow(bridge_audit), sum(bridge_audit$selected),
      sum(bridge_selected$date_match),
      sum(!bridge_selected$date_match))

  # --- Build expanded status map ---
  bridged_rows <- bridge_selected %>%
    mutate(tp_lab = tp_lab.y) %>%
    select(any_of(status_join_cols)) %>%
    mutate(uricult_bridge_applied = TRUE)

  status_map_expanded <- bind_rows(
    status_for_join %>% mutate(uricult_bridge_applied = FALSE),
    bridged_rows
  )

  # Resolve any duplicate Participant_id + tp_lab (bridged row wins)
  dup_keys <- status_map_expanded %>% count(Participant_id, tp_lab) %>% filter(n > 1)
  if (nrow(dup_keys) > 0) {
    msg("WARNING: Bridge created %d duplicate Participant_id + tp_lab keys; keeping bridged row.",
        nrow(dup_keys))
    write_csv(dup_keys, file.path(qc_diag_dir, "uricult_bridge_key_conflicts.csv"))
    status_map_expanded <- status_map_expanded %>%
      group_by(Participant_id, tp_lab) %>%
      arrange(desc(uricult_bridge_applied)) %>%
      slice_head(n = 1) %>%
      ungroup()
  }
} else {
  msg("No Uricult bridge rows needed (uricult=%d, unmatched UTI-N=%d).",
      nrow(uricult_clinical), nrow(vf_uti_n_unmatched))
  empty_bridge_audit <- tibble(
    Participant_id_clinical = character(),
    Episode_ID_clinical = character(),
    tp_lab_clinical = character(),
    UTI_Status = character(),
    UTI_binary = integer(),
    Not_UTI_subgroup = character(),
    Infection_Status = character(),
    Infection_Status_legacy = character(),
    UTI_Label = character(),
    Collection_Date_clinical = character(),
    mapped_tp_lab = character(),
    mapped_uti_n_num = integer(),
    Episode_ID_wgs = character(),
    wgs_date = character(),
    date_match = logical(),
    match_basis = character(),
    selected = logical(),
    selection_reason = character()
  )
  write_csv(empty_bridge_audit, file.path(qc_diag_dir, "uricult_bridge_audit.csv"))
  write_csv(
    empty_bridge_audit %>%
      mutate(reason_excluded = character()),
    file.path(qc_diag_dir, "uricult_bridge_sensitivity_alternatives.csv")
  )
  status_map_expanded <- status_for_join %>%
    mutate(uricult_bridge_applied = FALSE)
}

# ==============================================================================
# STEP 3a: Join VF P/A with primary UTI status (using expanded status map)
# ==============================================================================
#   left_join keeps all VF rows. Episodes without primary UTI status get NA
#   for UTI_Status / Infection_Status — this is expected and logged above.
vf_ready <- vf_pa %>%
  left_join(
    status_map_expanded,
    by = c("Participant_id", "tp_lab"),
    relationship = "many-to-one"
  )

if (all(c("Episode_ID.x", "Episode_ID.y") %in% names(vf_ready))) {
  vf_ready <- vf_ready %>%
    mutate(Episode_ID = coalesce(Episode_ID.y, Episode_ID.x)) %>%
    select(-Episode_ID.x, -Episode_ID.y)
}

vf_ready <- curation_default_columns(vf_ready)
vf_ready_excluded <- vf_ready %>%
  filter(!(analysis_include_primary %in% TRUE) | !(genomics_expected_include %in% TRUE))
if (nrow(vf_ready_excluded) > 0) {
  write_csv(vf_ready_excluded, file.path(qc_diag_dir, "vf_ready_manual_curation_excluded_rows.csv"))
  msg("Manual curation excludes %d VF-ready row(s) from the primary VF/model denominator.", nrow(vf_ready_excluded))
}
vf_ready <- vf_ready %>%
  filter(analysis_include_primary %in% TRUE, genomics_expected_include %in% TRUE)

# ==============================================================================
# STEP 3a-post: Duplicate safety checks
# ==============================================================================
# Check 1: unique Participant_id + tp_lab in vf_ready
vf_key_dupes <- vf_ready %>% count(Participant_id, tp_lab) %>% filter(n > 1)
if (nrow(vf_key_dupes) > 0) {
  write_csv(
    vf_ready %>% semi_join(vf_key_dupes, by = c("Participant_id", "tp_lab")),
    file.path(qc_diag_dir, "vf_ready_duplicate_participant_tp_keys.csv")
  )
  msg("WARNING: %d duplicate Participant_id + tp_lab in vf_ready.", nrow(vf_key_dupes))
}
# Check 2: unique clinical Episode_ID among status-stratified rows
vf_ep_dupes <- vf_ready %>%
  filter(!is.na(UTI_Status), !is.na(Episode_ID)) %>%
  count(Episode_ID) %>% filter(n > 1)
if (nrow(vf_ep_dupes) > 0) {
  write_csv(
    vf_ready %>% semi_join(vf_ep_dupes, by = "Episode_ID"),
    file.path(qc_diag_dir, "vf_ready_duplicate_clinical_episode_ids.csv")
  )
  msg("WARNING: %d clinical Episode_IDs duplicated in status-stratified vf_ready.", nrow(vf_ep_dupes))
}

# Post-bridge summary
n_uti_vf <- sum(vf_ready$UTI_Status == "UTI", na.rm = TRUE)
n_not_uti_vf <- sum(vf_ready$UTI_Status == "Not_UTI", na.rm = TRUE)
n_bridged_vf <- sum(vf_ready$uricult_bridge_applied %in% TRUE, na.rm = TRUE)
n_bridged_uti_vf <- sum(vf_ready$UTI_Status == "UTI" & vf_ready$uricult_bridge_applied %in% TRUE, na.rm = TRUE)
log_diag(
  "Post-bridge primary status rows: %d UTI, %d Not_UTI (UTI via Uricult bridge: %d; total bridged rows: %d)",
  n_uti_vf, n_not_uti_vf, n_bridged_uti_vf, n_bridged_vf
)
if (n_uti_vf < 20) {
  msg("WARNING: Primary UTI count in VF-ready data is %d (<20). Downstream association models are exploratory and may be unstable; this does not by itself indicate Uricult bridge failure.", n_uti_vf)
}

# STEP 3b: Add MLST Sequence Type if available
#   ST is needed by 25_vf_lineage_vf_interaction.R to check whether VF
#   differences are driven by lineage rather than primary UTI status.
if (!is.null(mlst)) {
  vf_ready <- vf_ready %>%
    left_join(mlst, by = c("Participant_id", "tp_lab"), relationship = "many-to-one")
  log_diag("MLST joined: %d rows with ST, %d without",
           sum(!is.na(vf_ready$ST)), sum(is.na(vf_ready$ST)))
} else {
  vf_ready$ST <- NA_character_
  vf_ready$ST_source <- NA_character_
  vf_ready$ST_provider <- NA_character_
  vf_ready$ST_local <- NA_character_
  vf_ready$provider_PercGoodTargets <- NA_real_
}

# STEP 3c: Compute total VF burden per episode
#   This is the count of all VF genes detected (including Unassigned).
#   It is the primary dependent variable for cross-sectional comparisons.
genes_not_in_map <- setdiff(gene_cols, gene_map$Gene)
curated_genes <- intersect(gene_cols, gene_map$Gene)
upec_pattern <- regex("adhesion|fimbr|iron|toxin|capsule|surface|invasion|evasion|serum|protectin|autotransporter|siderophore", ignore_case = TRUE)
upec_candidate_genes <- gene_map %>%
  filter(str_detect(Category, upec_pattern) | str_detect(coalesce(Subcategory, ""), upec_pattern)) %>%
  pull(Gene) %>%
  intersect(gene_cols)

vf_ready <- vf_ready %>%
  mutate(
    vf_count_total = rowSums(across(all_of(gene_cols)), na.rm = TRUE),
    total_vf_count_all = vf_count_total,
    total_vf_count_curated = if (length(curated_genes) > 0) rowSums(across(all_of(curated_genes)), na.rm = TRUE) else 0L,
    total_vf_count_upec_candidate = if (length(upec_candidate_genes) > 0) rowSums(across(all_of(upec_candidate_genes)), na.rm = TRUE) else 0L,
    total_vf_count_unassigned = if (length(genes_not_in_map) > 0) rowSums(across(all_of(genes_not_in_map)), na.rm = TRUE) else 0L,
    low_confidence_count = 0L
  )

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
if (length(genes_not_in_map) > 0) {
  vf_ready$cat_Unassigned_matrix <- rowSums(
    vf_ready[, genes_not_in_map, drop = FALSE], na.rm = TRUE
  )
  log_diag("Genes in matrix but NOT in gene_map: %d (counted in cat_Unassigned_matrix)",
           length(genes_not_in_map))
}

gap_report <- tibble(
  metric = c(
    "vf_gene_columns",
    "genes_in_gene_map",
    "genes_in_matrix_not_gene_map",
    "curated_genes_in_matrix",
    "upec_candidate_genes_in_matrix",
    "unassigned_fraction_of_matrix"
  ),
  value = c(
    length(gene_cols),
    nrow(gene_map),
    length(genes_not_in_map),
    length(curated_genes),
    length(upec_candidate_genes),
    round(length(genes_not_in_map) / max(1, length(gene_cols)), 4)
  ),
  notes = c(
    "All VF matrix columns retained",
    "Curated mapping file loaded by script 22",
    "Detected VF columns without curated category assignment",
    "Matrix genes present in gene_map",
    "Heuristic UPEC-candidate categories used for descriptive burden scores",
    "High values mean total burden is heavily influenced by unassigned genes"
  )
)
write_csv(gap_report, file.path(vf_diag_dir, "vf_gene_annotation_gap_report.csv"))
writeLines(
  c(
    "VF gene annotation gap report",
    sprintf("Generated: %s", format(Sys.time())),
    sprintf("VF matrix gene columns: %d", length(gene_cols)),
    sprintf("gene_map genes: %d", nrow(gene_map)),
    sprintf("Genes in matrix but not gene_map: %d", length(genes_not_in_map)),
    sprintf("Unassigned fraction of matrix: %.1f%%", 100 * length(genes_not_in_map) / max(1, length(gene_cols))),
    "Interpretation: unassigned/low-confidence genes are retained for transparency but should not be interpreted as curated UPEC-relevant features unless separately labelled."
  ),
  file.path(vf_diag_dir, "vf_gene_annotation_gap_report.txt")
)

# STEP 3f: Compute timepoint depth per participant
#   This counts how many distinct timepoints each participant contributed
#   (among episodes with primary UTI status). Used by 23_ and 24_ to stratify
#   analyses by ≥2, ≥3, ≥4 timepoints — important because participants with
#   more timepoints provide more information about longitudinal stability.
tp_depth <- vf_ready %>%
  filter(!is.na(UTI_Status)) %>%
  group_by(Participant_id) %>%
  summarise(n_timepoints = n_distinct(tp_lab), .groups = "drop")

vf_ready <- vf_ready %>%
  left_join(tp_depth, by = "Participant_id")

uti_vf_inspection <- vf_ready %>%
  filter(UTI_Status == "UTI") %>%
  select(any_of(c(
    "Participant_id", "tp_lab", "Episode_ID", "Timepoint", "Batch",
    "UTI_definition_version", "UTI_Status", "UTI_binary",
    "Not_UTI_subgroup", "uricult_bridge_applied", "mapped_tp_lab",
    "catheter_rule", "symptom_rule_met", "symptom_compatible_uti",
    "culture_supports_uti", "cfu_raw", "cfu_raw_parsed",
    "cfu_ge_1e3", "cfu_ge_1e5", "cfu_threshold_source",
    "UTI_classification_confidence", "UTI_classification_reason",
    "vf_count_total", "ST", "n_timepoints"
  ))) %>%
  arrange(Participant_id, tp_lab)
write_csv(uti_vf_inspection, file.path(DIR_VF, "uti_vf_episode_inspection.csv"))

# ==============================================================================
# 4. WRITE OUTPUTS
# ==============================================================================

vf_ready <- vf_ready %>%
  mutate(
    Infection_Status_legacy = if ("Infection_Status_legacy" %in% names(.)) as.character(Infection_Status_legacy) else as.character(Infection_Status),
    Infection_Status_old = if ("Infection_Status_old" %in% names(.)) as.character(Infection_Status_old) else Infection_Status_legacy,
    Infection_Status = as.character(UTI_Status),
    Primary_Status = as.character(UTI_Status)
  )

write_csv(vf_ready, out_file)

if (is_primary_run) {
  vf_binary_ready <- vf_ready %>%
    filter(!is.na(UTI_binary), UTI_Status %in% c("UTI", "Not_UTI"))
  write_csv(vf_binary_ready, FILE_VF_BINARY_UTI_READY)
  msg("Saved binary UTI-ready VF dataset (%d rows) to %s", nrow(vf_binary_ready), FILE_VF_BINARY_UTI_READY)
}

# Log a comprehensive summary of what was produced
log_diag("")
log_diag("=== OUTPUT SUMMARY ===")
log_diag("File: %s", out_file)
log_diag("Rows: %d", nrow(vf_ready))
log_diag("Participants: %d", n_distinct(vf_ready$Participant_id))
log_diag("Gene columns: %d", length(gene_cols))
log_diag("With UTI_Status: %d", sum(!is.na(vf_ready$UTI_Status)))
log_diag("Without UTI_Status: %d", sum(is.na(vf_ready$UTI_Status)))
log_diag("Primary UTI: %d", sum(vf_ready$UTI_Status == "UTI", na.rm = TRUE))
log_diag("Primary Not_UTI: %d", sum(vf_ready$UTI_Status == "Not_UTI", na.rm = TRUE))
legacy_status_vec <- if ("Infection_Status_legacy" %in% names(vf_ready)) {
  vf_ready$Infection_Status_legacy
} else {
  vf_ready$Infection_Status
}
log_diag("Legacy ASB: %d", sum(legacy_status_vec == "ASB", na.rm = TRUE))
log_diag("Legacy UTI: %d", sum(legacy_status_vec == "UTI", na.rm = TRUE))
log_diag("Legacy Negative: %d", sum(legacy_status_vec == "Negative", na.rm = TRUE))
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
writeLines(diag_lines, diag_file)
msg("Diagnostics written to %s", diag_file)

append_denominator_summary(
  vf_pa,
  "22_vf_build_analysis_dataset.R",
  paste0("vf_pa_input", stage_suffix),
  "participant_timepoint",
  opt$vf_pa,
  "VF P/A rows after preserved-event timepoint normalisation"
)
append_denominator_summary(
  vf_ready,
  "22_vf_build_analysis_dataset.R",
  paste0("vf_analysis_ready", stage_suffix),
  "participant_timepoint",
  out_file,
  "Canonical VF-ready dataset; no status_map first-row deduplication"
)
if (is_primary_run) {
  write_uti_attrition_outputs()
}

msg("✓ 22_vf_build_analysis_dataset.R complete.")

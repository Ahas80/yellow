#!/usr/bin/env Rscript

source("00_config.R")

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
  library(tidyr)
})

out_dir <- "results/audit"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

overview_workbook <- "OVERVIEW E.coli batch 1-6 - CLEAN_RC.xlsx"

collapse_nonmissing <- function(x) {
  x <- unique(as.character(x))
  x <- x[!is.na(x) & nzchar(x) & x != "."]
  if (length(x) == 0) {
    NA_character_
  } else {
    paste(x, collapse = "; ")
  }
}

is_true <- function(x) {
  tolower(as.character(x)) %in% c("true", "t", "1", "yes")
}

md_table <- function(df) {
  if (nrow(df) == 0) {
    return("_None._")
  }

  df <- as.data.frame(df, stringsAsFactors = FALSE)
  df[] <- lapply(df, function(x) ifelse(is.na(x), "", as.character(x)))

  header <- paste0("| ", paste(names(df), collapse = " | "), " |")
  sep <- paste0("| ", paste(rep("---", ncol(df)), collapse = " | "), " |")
  rows <- apply(df, 1, function(x) paste0("| ", paste(x, collapse = " | "), " |"))

  paste(c(header, sep, rows), collapse = "\n")
}

diag_summary_path <- file.path("results", "vf", "uti_not_uti_diagnostic_summary.csv")
diag_interpretation_path <- file.path("results", "vf", "uti_not_uti_test_interpretation_table.csv")
diag_fig_meta_path <- file.path("results", "vf", "uti_not_uti_diagnostic_figure_metadata.csv")
diag_summary <- if (file.exists(diag_summary_path)) {
  read_csv(diag_summary_path, show_col_types = FALSE)
} else {
  tibble()
}
diag_interpretation <- if (file.exists(diag_interpretation_path)) {
  read_csv(diag_interpretation_path, show_col_types = FALSE)
} else {
  tibble()
}
diag_fig_meta <- if (file.exists(diag_fig_meta_path)) {
  read_csv(diag_fig_meta_path, show_col_types = FALSE)
} else {
  tibble()
}

status_all <- read_csv("results/clinical/status_map.csv", show_col_types = FALSE) %>%
  mutate(
    Participant_id = as.character(Participant_id),
    tp_lab = as.character(tp_lab)
  ) %>%
  apply_manual_sample_curation(context = "uti_count_audit_status")
status <- read_csv(FILE_ANALYSIS_CLINICAL_COHORT, show_col_types = FALSE) %>%
  prefer_primary_uti_status(allow_legacy_fallback = FALSE) %>%
  mutate(
    Participant_id = as.character(.data$Participant_id),
    tp_lab = normalise_timepoint_preserve_events(.data$tp_lab)
  )
active_longcycler_keys <- status %>%
  distinct(.data$Participant_id, .data$tp_lab)
if (nrow(active_longcycler_keys) != nrow(status)) {
  stop("Selected Longcycler analysis cohort contains duplicate episode keys.")
}

vf_all <- read_csv("results/vf/vf_analysis_ready.csv", show_col_types = FALSE) %>%
  mutate(
    Participant_id = as.character(Participant_id),
    tp_lab = as.character(tp_lab)
  ) %>%
  apply_manual_sample_curation(context = "uti_count_audit_vf")
vf <- filter_primary_genomics(vf_all) %>%
  mutate(tp_lab = normalise_timepoint_preserve_events(.data$tp_lab))
vf_keys <- vf %>% distinct(.data$Participant_id, .data$tp_lab)
if (nrow(vf) != nrow(active_longcycler_keys) ||
    nrow(anti_join(active_longcycler_keys, vf_keys, by = c("Participant_id", "tp_lab"))) ||
    nrow(anti_join(vf_keys, active_longcycler_keys, by = c("Participant_id", "tp_lab")))) {
  stop("VF-ready keys do not exactly equal the selected Longcycler analysis cohort.")
}

bridge <- read_csv("results/qc/uricult_bridge_audit.csv", show_col_types = FALSE) %>%
  mutate(
    Participant_id = as.character(Participant_id_clinical),
    tp_lab = as.character(mapped_tp_lab)
  )

meta <- read_csv("assembly_metadata.csv", show_col_types = FALSE) %>%
  mutate(
    Participant_id = as.character(Participant_id),
    tp_lab = as.character(tp_lab)
  ) %>%
  apply_manual_sample_curation(context = "uti_count_audit_metadata")
meta_all <- meta
meta <- filter_primary_genomics(meta)

vf_pa <- read_csv("results/vf/vf_pa_all.csv", show_col_types = FALSE) %>%
  mutate(
    Participant_id = as.character(Participant_id),
    tp_lab = as.character(tp_lab)
  ) %>%
  apply_manual_sample_curation(context = "uti_count_audit_vf_pa") %>%
  filter_primary_genomics() %>%
  mutate(tp_lab = normalise_timepoint_preserve_events(.data$tp_lab))
vf_pa_keys <- vf_pa %>% distinct(.data$Participant_id, .data$tp_lab)
if (nrow(vf_pa_keys) != nrow(vf_pa) ||
    nrow(anti_join(active_longcycler_keys, vf_pa_keys, by = c("Participant_id", "tp_lab"))) ||
    nrow(anti_join(vf_pa_keys, active_longcycler_keys, by = c("Participant_id", "tp_lab")))) {
  stop("VF presence/absence keys do not exactly equal the selected Longcycler analysis cohort.")
}

clinical_excluded_primary <- status_all %>%
  filter(!(analysis_include_primary %in% TRUE)) %>%
  arrange(Participant_id, tp_lab) %>%
  select(any_of(c(
    "Participant_id", "Timepoint", "tp_lab", "Collection_Date", "Batch",
    "UTI_Label", "UTI_Status", "Not_UTI_subgroup", "Infection_Status_legacy",
    "analysis_include_primary", "analysis_exclusion_reason",
    "duplicate_role", "duplicate_of_participant_id", "duplicate_of_tp_lab",
    "allow_secondary_duplicate_qc", "duplicate_use_note", "manual_curation_note"
  )))
write_csv(clinical_excluded_primary, file.path(out_dir, "primary_clinical_manual_exclusions.csv"))

vf_excluded_in_ready <- vf_all %>%
  filter(!(analysis_include_primary %in% TRUE) | !(genomics_expected_include %in% TRUE)) %>%
  arrange(Participant_id, tp_lab) %>%
  select(any_of(c(
    "Participant_id", "tp_lab", "Episode_ID", "Collection_Date", "Batch",
    "UTI_Label", "UTI_Status", "Not_UTI_subgroup", "ST",
    "analysis_include_primary", "analysis_exclusion_reason",
    "genomics_expected_include", "genomics_exclusion_reason",
    "duplicate_role", "duplicate_of_participant_id", "duplicate_of_tp_lab",
    "allow_secondary_duplicate_qc", "duplicate_use_note", "manual_curation_note"
  )))

vf_excluded_upstream <- meta_all %>%
  filter(!(analysis_include_primary %in% TRUE)) %>%
  group_by(Participant_id, tp_lab) %>%
  summarise(
    Isolate_ID = collapse_nonmissing(Isolate_ID),
    n_metadata_rows = n(),
    n_found_fastas = sum(found %in% TRUE, na.rm = TRUE),
    analysis_include_primary = all(analysis_include_primary %in% TRUE),
    analysis_exclusion_reason = collapse_nonmissing(analysis_exclusion_reason),
    genomics_expected_include = all(genomics_expected_include %in% TRUE),
    genomics_exclusion_reason = collapse_nonmissing(genomics_exclusion_reason),
    duplicate_role = collapse_nonmissing(duplicate_role),
    duplicate_of_participant_id = collapse_nonmissing(duplicate_of_participant_id),
    duplicate_of_tp_lab = collapse_nonmissing(duplicate_of_tp_lab),
    allow_secondary_duplicate_qc = any(allow_secondary_duplicate_qc %in% TRUE, na.rm = TRUE),
    duplicate_use_note = collapse_nonmissing(duplicate_use_note),
    manual_curation_note = collapse_nonmissing(manual_curation_note),
    .groups = "drop"
  ) %>%
  left_join(
    status_all %>%
      select(any_of(c(
        "Participant_id", "tp_lab", "Episode_ID", "Collection_Date", "Batch",
        "UTI_Label", "UTI_Status", "Not_UTI_subgroup", "ST"
      ))),
    by = c("Participant_id", "tp_lab")
  ) %>%
  mutate(exclusion_layer = "excluded_upstream_before_vf_ready")

vf_excluded_primary <- bind_rows(
  vf_excluded_in_ready %>%
    mutate(exclusion_layer = "present_in_vf_ready_but_excluded") %>%
    mutate(across(everything(), as.character)),
  vf_excluded_upstream %>%
    mutate(across(everything(), as.character))
) %>%
  arrange(Participant_id, tp_lab) %>%
  select(any_of(c(
    "exclusion_layer", "Participant_id", "tp_lab", "Episode_ID",
    "Collection_Date", "Batch", "UTI_Label", "UTI_Status",
    "Not_UTI_subgroup", "ST", "Isolate_ID", "n_metadata_rows",
    "n_found_fastas", "analysis_include_primary",
    "analysis_exclusion_reason", "genomics_expected_include",
    "genomics_exclusion_reason", "duplicate_role",
    "duplicate_of_participant_id", "duplicate_of_tp_lab",
    "allow_secondary_duplicate_qc", "duplicate_use_note",
    "manual_curation_note"
  )))
write_csv(vf_excluded_primary, file.path(out_dir, "primary_vf_manual_exclusions.csv"))

quarantined_fastas <- meta_all %>%
  filter(!(genomics_expected_include %in% TRUE)) %>%
  arrange(Batch, Participant_id, tp_lab, Isolate_ID) %>%
  select(any_of(c(
    "Participant_id", "tp_lab", "Timepoint", "Isolate_ID", "Batch",
    "Collection_Date", "UTI_Label", "Clinical_CFU_Count", "Clinical_Beoord",
    "Urine_collection_method", "found", "file_exists", "usable_fasta",
    "genomics_expected_include", "genomics_exclusion_reason",
    "manual_curation_note"
  )))
write_csv(quarantined_fastas, file.path(out_dir, "quarantined_failed_or_not_expected_fastas.csv"))

clinical_counts <- status %>%
  count(UTI_Status, name = "n")

vf_counts <- vf %>%
  mutate(UTI_Status = if_else(is.na(UTI_Status), "<missing>", UTI_Status)) %>%
  count(UTI_Status, name = "n")

legacy_cross <- status %>%
  count(Infection_Status_legacy, UTI_Status, name = "n") %>%
  arrange(Infection_Status_legacy, UTI_Status)

vf_legacy_cross <- vf %>%
  mutate(
    UTI_Status = if_else(is.na(UTI_Status), "<missing>", UTI_Status),
    Infection_Status_legacy = if_else(
      is.na(Infection_Status_legacy),
      "<missing>",
      Infection_Status_legacy
    )
  ) %>%
  count(Infection_Status_legacy, UTI_Status, name = "n") %>%
  arrange(Infection_Status_legacy, UTI_Status)

primary_uti <- status %>%
  filter(UTI_Status == "UTI") %>%
  arrange(Participant_id, tp_lab) %>%
  select(any_of(c(
    "Participant_id", "Timepoint", "tp_lab", "Collection_Date", "Batch",
    "UTI_Label", "Infection_Status_legacy", "UTI_Status", "UTI_binary",
    "Urine_collection_method", "urine_collection_method_norm", "catheter_rule",
    "culture_supports_uti", "cfu_raw", "cfu_raw_parsed", "cfu_ge_1e3",
    "cfu_ge_1e5", "symptom_compatible_uti", "symptom_rule_met",
    "local_urinary_symptom_any", "systemic_symptom_any", "dysuria_present",
    "urgency_present", "frequency_present", "incontinence_present",
    "pus_present", "flankpain_present", "suprapubic_pain_present",
    "fever_present", "rigors_present", "delirium_present",
    "UTI_classification_reason", "Episode_ID"
  )))

write_csv(primary_uti, file.path(out_dir, "primary_uti_rows.csv"))

participant_100011_status <- status %>%
  filter(Participant_id == "100011") %>%
  mutate(
    Not_UTI_subgroup = coalesce(Not_UTI_subgroup, ""),
    status_reason = case_when(
      UTI_Status == "UTI" ~ "Culture support and compatible symptoms",
      culture_supports_uti == TRUE & symptom_compatible_uti == FALSE ~
        "Culture support, but symptoms not compatible",
      culture_supports_uti == TRUE & is.na(symptom_compatible_uti) ~
        "Culture support, but symptom rule unknown",
      culture_supports_uti == FALSE ~
        "Culture does not meet current threshold",
      TRUE ~ "Indeterminate under current rule"
    )
  ) %>%
  mutate(
    timepoint_group = case_when(
      str_detect(tp_lab, "^T[0-9]+$") ~ 0L,
      str_detect(tp_lab, "^UTI-[0-9]+$") ~ 1L,
      TRUE ~ 2L
    ),
    timepoint_number = suppressWarnings(as.integer(str_extract(tp_lab, "[0-9]+$")))
  ) %>%
  arrange(timepoint_group, timepoint_number, tp_lab) %>%
  select(any_of(c(
    "Participant_id", "tp_lab", "Collection_Date", "UTI_Label",
    "UTI_Status", "Not_UTI_subgroup", "UTI_binary",
    "Infection_Status_legacy", "Urine_collection_method", "catheter_rule",
    "culture_supports_uti", "cfu_threshold_source",
    "symptom_compatible_uti", "local_urinary_symptom_any",
    "systemic_symptom_any", "status_reason"
  )))

write_csv(participant_100011_status, file.path(out_dir, "participant_100011_status_audit.csv"))

legacy_reclass <- status %>%
  filter(Infection_Status_legacy == "UTI", UTI_Status == "Not_UTI") %>%
  mutate(
    exclusion_bucket = case_when(
      symptom_compatible_uti == FALSE ~ "Culture supports UTI, but symptom rule not met",
      is.na(symptom_compatible_uti) ~ "Culture supports UTI, but symptom rule unknown",
      TRUE ~ "Other"
    )
  ) %>%
  arrange(exclusion_bucket, Participant_id, tp_lab) %>%
  select(any_of(c(
    "Participant_id", "Timepoint", "tp_lab", "Collection_Date", "Batch",
    "UTI_Label", "Infection_Status_legacy", "UTI_Status", "Not_UTI_subgroup",
    "exclusion_bucket", "Urine_collection_method", "urine_collection_method_norm",
    "catheter_rule", "culture_supports_uti", "cfu_raw", "cfu_raw_parsed",
    "cfu_ge_1e3", "cfu_ge_1e5", "symptom_compatible_uti",
    "symptom_rule_met", "local_urinary_symptom_any", "systemic_symptom_any",
    "dysuria_present", "urgency_present", "frequency_present",
    "incontinence_present", "pus_present", "flankpain_present",
    "suprapubic_pain_present", "fever_present", "rigors_present",
    "delirium_present", "UTI_classification_reason", "Episode_ID"
  )))

write_csv(legacy_reclass, file.path(out_dir, "legacy_uti_reclassified_not_uti.csv"))

legacy_reclass_counts <- legacy_reclass %>%
  count(exclusion_bucket, name = "n")

criteria_mismatch <- status %>%
  filter(
    culture_supports_uti == TRUE,
    symptom_compatible_uti == TRUE,
    UTI_Status != "UTI"
  ) %>%
  select(any_of(c(
    "Participant_id", "Timepoint", "tp_lab", "UTI_Status",
    "Infection_Status_legacy", "culture_supports_uti",
    "symptom_compatible_uti", "symptom_rule_met",
    "UTI_classification_reason", "Episode_ID"
  )))

write_csv(
  criteria_mismatch,
  file.path(out_dir, "status_map_rows_meeting_uti_criteria_but_not_labelled_uti.csv")
)

bridge_selected <- bridge %>%
  filter(is_true(selected)) %>%
  transmute(
    bridge_clinical_episode = Episode_ID_clinical,
    selected_counterpart_tp_lab = mapped_tp_lab,
    selected_counterpart_wgs_episode = Episode_ID_wgs,
    selected_counterpart_status = UTI_Status,
    selected_counterpart_legacy_status = Infection_Status_legacy
  )

bridge_alt <- bridge %>%
  filter(!is_true(selected)) %>%
  transmute(
    Participant_id,
    tp_lab = mapped_tp_lab,
    bridge_clinical_episode = Episode_ID_clinical,
    bridge_clinical_tp_lab = tp_lab_clinical,
    bridge_clinical_status = UTI_Status,
    bridge_clinical_legacy_status = Infection_Status_legacy,
    bridge_not_uti_subgroup = Not_UTI_subgroup,
    bridge_collection_date = Collection_Date_clinical,
    bridge_uti_label = UTI_Label,
    bridge_wgs_episode = Episode_ID_wgs,
    bridge_wgs_date = wgs_date,
    bridge_date_match = date_match,
    bridge_match_basis = match_basis,
    bridge_selection_reason = selection_reason
  ) %>%
  left_join(bridge_selected, by = "bridge_clinical_episode")

meta_sum <- meta %>%
  group_by(Participant_id, tp_lab) %>%
  summarise(
    n_metadata_assemblies = n_distinct(Assembly_ID),
    assembly_base_ids = collapse_nonmissing(Assembly_Base_ID),
    isolate_ids = collapse_nonmissing(Isolate_ID),
    metadata_collection_date = collapse_nonmissing(Collection_Date),
    metadata_batch = collapse_nonmissing(Batch),
    metadata_uti_label = collapse_nonmissing(UTI_Label),
    metadata_clinical_cfu_count = collapse_nonmissing(Clinical_CFU_Count),
    metadata_clinical_beoord = collapse_nonmissing(Clinical_Beoord),
    metadata_clinical_organism = collapse_nonmissing(Clinical_Organism),
    metadata_urine_collection_method = collapse_nonmissing(Urine_collection_method),
    metadata_population = collapse_nonmissing(Population),
    metadata_obj = collapse_nonmissing(Obj),
    .groups = "drop"
  )

vf_missing <- vf %>%
  filter(is.na(UTI_Status)) %>%
  select(any_of(c(
    "Participant_id", "tp_lab", "Episode_ID", "Event_type",
    "Collection_Date", "Batch", "ST", "Infection_Status",
    "Infection_Status_legacy", "UTI_Label", "uricult_bridge_applied",
    "vf_count_total", "total_vf_count_curated", "Primary_Status"
  ))) %>%
  left_join(bridge_alt, by = c("Participant_id", "tp_lab")) %>%
  left_join(meta_sum, by = c("Participant_id", "tp_lab")) %>%
  mutate(
    missing_status_category = if_else(
      !is.na(bridge_clinical_episode),
      "alternative_unselected_uricult_wgs_isolate",
      "wgs_only_uti_event_not_in_status_map"
    ),
    can_receive_primary_status_now = case_when(
      missing_status_category == "alternative_unselected_uricult_wgs_isolate" ~
        "No for primary denominator; it is an alternate WGS isolate for a clinical episode already represented by the selected mapped isolate.",
      TRUE ~
        "No under current criteria; no matching clinical status-map episode/symptom evidence is present."
    ),
    possible_unlabelled_uti_flag = case_when(
      missing_status_category == "wgs_only_uti_event_not_in_status_map" ~
        "Yes: candidate WGS-only UTI event; needs clinical episode/symptom evidence before primary labelling.",
      TRUE ~
        "No: status belongs to the bridged clinical episode, but this row is deliberately not selected for the primary denominator."
    )
  ) %>%
  arrange(missing_status_category, Participant_id, tp_lab)

write_csv(vf_missing, file.path(out_dir, "vf_ready_missing_status_rows.csv"))

missing_counts <- vf_missing %>%
  count(missing_status_category, name = "n")

wgs_only <- vf_missing %>%
  filter(missing_status_category == "wgs_only_uti_event_not_in_status_map")

write_csv(wgs_only, file.path(out_dir, "wgs_only_possible_unlabelled_uti_rows.csv"))

wgs_only_status_context <- status %>%
  semi_join(wgs_only %>% distinct(Participant_id), by = "Participant_id") %>%
  arrange(Participant_id, tp_lab) %>%
  select(any_of(c(
    "Participant_id", "Timepoint", "tp_lab", "Collection_Date", "Batch",
    "UTI_Label", "Infection_Status_legacy", "UTI_Status",
    "Not_UTI_subgroup", "culture_supports_uti", "symptom_compatible_uti",
    "symptom_rule_met", "UTI_classification_reason", "Episode_ID"
  )))

write_csv(
  wgs_only_status_context,
  file.path(out_dir, "clinical_status_rows_for_wgs_only_participants.csv")
)

status_keys <- status %>%
  distinct(Participant_id, tp_lab)

direct_unmatched <- vf_pa %>%
  anti_join(status_keys, by = c("Participant_id", "tp_lab")) %>%
  filter(str_detect(tp_lab, "^UTI-")) %>%
  select(Participant_id, tp_lab, Episode_ID) %>%
  left_join(
    bridge %>%
      transmute(
        Participant_id,
        tp_lab = mapped_tp_lab,
        bridge_selected = is_true(selected),
        bridge_clinical_episode = Episode_ID_clinical,
        bridge_clinical_status = UTI_Status,
        bridge_selection_reason = selection_reason
      ),
    by = c("Participant_id", "tp_lab")
  ) %>%
  mutate(
    resolution = case_when(
      bridge_selected == TRUE ~ "selected_by_uricult_bridge_and_inherits_clinical_status",
      bridge_selected == FALSE ~ "alternate_uricult_bridge_candidate_not_selected",
      TRUE ~ "no_clinical_status_or_uricult_bridge_match"
    )
  ) %>%
  arrange(resolution, Participant_id, tp_lab)

write_csv(
  direct_unmatched,
  file.path(out_dir, "vf_uti_n_direct_unmatched_bridge_resolution.csv")
)

direct_unmatched_counts <- direct_unmatched %>%
  count(resolution, name = "n")

beoord_fallback_rows <- status %>%
  filter(str_detect(coalesce(cfu_threshold_source, ""), "^beoord_plus[123]_fallback_1e[345]$")) %>%
  arrange(cfu_threshold_source, UTI_Status, Participant_id, tp_lab) %>%
  select(any_of(c(
    "Participant_id", "Timepoint", "tp_lab", "Collection_Date", "Batch",
    "UTI_Label", "Infection_Status_legacy", "UTI_Status", "Not_UTI_subgroup",
    "beoord_cat", "cfu_raw", "cfu_threshold_source", "culture_supports_uti",
    "symptom_compatible_uti", "symptom_rule_met", "UTI_classification_reason",
    "Episode_ID"
  )))
write_csv(beoord_fallback_rows, file.path(out_dir, "beoord_fallback_primary_culture_support_rows.csv"))

beoord_fallback_counts <- beoord_fallback_rows %>%
  count(cfu_threshold_source, UTI_Status, Not_UTI_subgroup, name = "n") %>%
  arrange(cfu_threshold_source, UTI_Status, Not_UTI_subgroup)

clinical_total <- nrow(status)
clinical_uti <- sum(status$UTI_Status == "UTI", na.rm = TRUE)
clinical_not_uti <- sum(status$UTI_Status == "Not_UTI", na.rm = TRUE)
clinical_missing <- sum(is.na(status$UTI_Status))
clinical_raw_total <- nrow(status_all)
clinical_excluded_n <- nrow(clinical_excluded_primary)
vf_total <- nrow(vf)
vf_uti <- sum(vf$UTI_Status == "UTI", na.rm = TRUE)
vf_not_uti <- sum(vf$UTI_Status == "Not_UTI", na.rm = TRUE)
vf_missing_n <- sum(is.na(vf$UTI_Status))
vf_raw_total <- nrow(vf_all)
vf_excluded_n <- nrow(vf_excluded_primary)
quarantined_fasta_n <- nrow(quarantined_fastas)
legacy_uti_total <- sum(status$Infection_Status_legacy == "UTI", na.rm = TRUE)
criteria_mismatch_n <- nrow(criteria_mismatch)
alt_missing_n <- sum(
  vf_missing$missing_status_category == "alternative_unselected_uricult_wgs_isolate",
  na.rm = TRUE
)
wgs_only_n <- nrow(wgs_only)
direct_unmatched_n <- nrow(direct_unmatched)

missing_status_text <- if (vf_missing_n == 0) {
  c(
    "No VF rows currently have missing `UTI_Status`. Every VF row now has a direct `Participant_id + tp_lab` match in the clinical status map, so no Uricult bridge rows are needed for the current inputs."
  )
} else {
  c(
    sprintf(
      "%d VF rows currently have missing `UTI_Status`. These are WGS/VF rows whose `Participant_id + tp_lab` key did not receive a clinical status after direct matching and Uricult bridging.",
      vf_missing_n
    ),
    sprintf(
      "%d missing rows are alternate UTI-N WGS isolates for Uricult clinical episodes where another mapped UTI-N isolate was selected as the primary representative.",
      alt_missing_n
    ),
    sprintf(
      "%d missing rows are WGS-only `UTI_event` rows with UTI-like metadata but no matching clinical episode/status-map row.",
      wgs_only_n
    )
  )
}

unlabelled_text <- if (wgs_only_n == 0) {
  "No WGS-only possible unlabelled UTI rows remain in `vf_analysis_ready.csv`."
} else {
  sprintf(
    "These %d rows are the cleanest candidates to manually curate as possible missing clinical UTI episodes.",
    wgs_only_n
  )
}

bridge_text <- if (direct_unmatched_n == 0) {
  "No UTI-N VF rows are currently unmatched to the clinical status map. The clean workbook inputs make the previous bridge unnecessary."
} else {
  sprintf(
    "Before the Uricult bridge, %d UTI-N VF rows did not directly match the clinical status map.",
    direct_unmatched_n
  )
}

count_logic_reference <- tibble::tribble(
  ~count_name, ~row_unit, ~source_table, ~logic_followed,
  "primary clinical_status_map total",
  "selected genomic episode",
  "results/clinical/analysis_cohort_longcycler.csv",
  "Count the exact selected QC-pass Longcycler participant-timepoint keys after one-to-one clinical status linkage.",
  "primary clinical UTI",
  "clinical episode",
  "results/clinical/analysis_cohort_longcycler.csv",
  "Within selected Longcycler rows, count UTI_Status == 'UTI'. This requires culture_supports_uti == TRUE and symptom_compatible_uti == TRUE.",
  "primary clinical Not_UTI",
  "clinical episode",
  "results/clinical/analysis_cohort_longcycler.csv",
  "Within selected Longcycler rows, count UTI_Status == 'Not_UTI'.",
  "primary VF-ready total",
  "sequenced VF participant-timepoint",
  "results/vf/vf_analysis_ready.csv",
  "Count every row in vf_pa_all after joining primary clinical status by Participant_id + tp_lab, then keep analysis_include_primary == TRUE and genomics_expected_include == TRUE.",
  "primary VF-ready UTI",
  "sequenced VF participant-timepoint",
  "results/vf/vf_analysis_ready.csv",
  "Within primary-included/genomics-expected VF rows, count joined UTI_Status == 'UTI'. These are sequenced rows whose clinical episode meets the primary UTI rule.",
  "primary VF-ready Not_UTI",
  "sequenced VF participant-timepoint",
  "results/vf/vf_analysis_ready.csv",
  "Within primary-included/genomics-expected VF rows, count joined UTI_Status == 'Not_UTI'.",
  "primary VF-ready missing status",
  "sequenced VF participant-timepoint",
  "results/vf/vf_analysis_ready.csv",
  "Within primary-included/genomics-expected VF rows, count rows where UTI_Status is NA after joining by Participant_id + tp_lab. A missing row means no matching clinical status-map key, not a clinical row that failed UTI criteria.",
  "legacy UTI",
  "clinical episode",
  "results/clinical/status_map.csv",
  "Count rows where Infection_Status_legacy == 'UTI'. This is retained for audit and follows the older culture-positive plus broad S&S/population logic, not the current primary UTI rule.",
  "legacy UTI reclassified to Not_UTI",
  "clinical episode",
  "results/clinical/status_map.csv",
  "Count rows where Infection_Status_legacy == 'UTI' and UTI_Status == 'Not_UTI'. These are culture-supported but fail or lack the current catheter-aware symptom requirement.",
  "possible unlabelled UTI candidates",
  "sequenced VF participant-timepoint",
  "results/audit/wgs_only_possible_unlabelled_uti_rows.csv",
  "Count missing-status VF rows whose UTI-N WGS event has no matching clinical status-map row. Current count is zero after repopulating batch inputs from the clean workbooks.",
  "manual primary exclusions",
  "clinical/VF row",
  "data/manual_sample_curation.csv",
  "Rows with exclude_primary == TRUE are retained in source/audit tables but removed from primary clinical, VF, and model denominators.",
  "quarantined failed/not-expected FASTA rows",
  "expected isolate",
  "data/manual_sample_curation.csv",
  "Rows with exclude_from_genomics_expected == TRUE are retained in metadata quarantine reports but removed from active genomics expected-denominator checks."
)
write_csv(count_logic_reference, file.path(out_dir, "uti_count_logic_reference.csv"))

report <- c(
  "# UTI / Not_UTI Count Audit",
  "",
  "Generated from the selected cohort in `results/clinical/analysis_cohort_longcycler.csv` and its matching VF outputs. The broader status map and assembly metadata are used only for source attrition/QC.",
  "",
  "## Bottom line",
  "",
  sprintf(
    "- The selected QC-pass Longcycler denominator has %d episodes: %d `UTI`, %d `Not_UTI`, and %d missing primary statuses.",
    clinical_total, clinical_uti, clinical_not_uti, clinical_missing
  ),
  sprintf(
    "- The raw/audit clinical status map still has %d rows; %d are excluded from primary analyses but retained for traceability.",
    clinical_raw_total, clinical_excluded_n
  ),
  sprintf(
    "- The primary VF analysis-ready table has %d WGS/VF rows after primary/genomics filters: %d `UTI`, %d `Not_UTI`, and %d missing status rows.",
    vf_total, vf_uti, vf_not_uti, vf_missing_n
  ),
  sprintf(
    "- The final raw/audit VF table has %d rows after upstream filtering; %d clinical/genomics rows are listed as excluded from VF/model construction.",
    vf_raw_total, vf_excluded_n
  ),
  sprintf(
    "- The failed/not-expected FASTA quarantine has %d expected isolate row(s) outside active genomics completeness checks.",
    quarantined_fasta_n
  ),
  sprintf(
    "- There are %d clinical status-map rows that meet the current primary UTI criterion (`culture_supports_uti == TRUE` and `symptom_compatible_uti == TRUE`) but are labelled `Not_UTI`.",
    criteria_mismatch_n
  ),
  "",
  "## Current counts",
  "",
  "Selected Longcycler analysis cohort:",
  md_table(clinical_counts),
  "",
  "Primary VF analysis-ready table:",
  md_table(vf_counts),
  "",
  "Manual primary clinical exclusions:",
  md_table(clinical_excluded_primary),
  "",
  "Manual primary VF/model exclusions:",
  md_table(vf_excluded_primary),
  "",
  "Quarantined failed/not-expected FASTA rows:",
  md_table(quarantined_fastas),
  "",
  "Legacy clinical label to primary status:",
  md_table(legacy_cross),
  "",
  "VF legacy label to primary status:",
  md_table(vf_legacy_cross),
  "",
  "## Exact logic behind the counts",
  "",
  "### Input rows",
  "",
  "- `data/inputs/batch1.csv` to `batch6.csv` are regenerated from the workbook inputs before classification.",
  sprintf("- Batches 1 to 3 use the three named sequencing workbooks, but `%s` is treated as the canonical source for corrected participant IDs, UTI-N timepoints, collection dates, isolate IDs, organism, CFU/growth fields, urine collection method, and symptom columns.", overview_workbook),
  "- For batches 1 to 3, `Population` and `UWI#` are supplemented from the individual workbook rows when available, because the clean RC overview does not always retain those fields.",
  sprintf("- Batches 4 to 6 are populated directly from `%s`.", overview_workbook),
  "",
  "### Clinical episode rows",
  "",
  sprintf(
    "- The analytical denominator is `nrow(analysis_cohort_longcycler.csv)`: currently %d rows.",
    clinical_total
  ),
  "- `00b_classify_episodes.R` groups cleaned clinical rows into one row per `Participant_id + Timepoint` / `tp_lab` episode.",
  "- It computes `UTI_Status` first, then applies `data/manual_sample_curation.csv` to set primary inclusion flags.",
  "- The count `primary clinical UTI` is `sum(UTI_Status == 'UTI')` within the selected Longcycler cohort.",
  "- The count `primary clinical Not_UTI` is `sum(UTI_Status == 'Not_UTI')` within the selected Longcycler cohort.",
  "- Excluded rows keep their computed clinical status in the audit table; they are not relabelled.",
  "",
  "### Primary UTI rule",
  "",
  "A row is primary `UTI` only if both sides are true:",
  "",
  "```r",
  "culture_supports_uti == TRUE & symptom_compatible_uti == TRUE",
  "```",
  "",
  "Culture support is derived as follows:",
  "",
  "- If a usable CFU value is present, the current primary threshold is `>= 1e3 CFU/mL` using the parsed lower bound.",
  "- If CFU is missing, the `Beoord` column can provide a culture-support fallback using `+ = 10^3`, `++ = 10^4`, and `+++ = 10^5`.",
  "- Because the primary threshold is `>= 1e3 CFU/mL`, any parsed `Beoord` category of `+`, `++`, or `+++` supports the culture side of the current UTI rule.",
  "- Ambiguous or missing culture evidence remains unknown.",
  "",
  "Symptom compatibility is catheter-aware:",
  "",
  "- For non-indwelling urine collection (`spontaan_geloosd` or `inco`), the row needs a local urinary symptom: dysuria, urgency, frequency, incontinence, or pus. Flank pain only counts if paired with a systemic symptom.",
  "- For indwelling catheter rows, the row needs a systemic symptom: fever, rigors/chills, or delirium.",
  "- If the urine collection method is unknown, the symptom rule is `Unknown`, so the row cannot become primary UTI under the current criteria.",
  "",
  "Final status logic:",
  "",
  "- `UTI_Status = 'UTI'` when culture support and symptom compatibility are both true.",
  "- Every other classified clinical episode is `UTI_Status = 'Not_UTI'`.",
  "- `Not_UTI_subgroup` explains why: usually `bacteriuria_not_UTI`, `culture_negative_or_below_threshold`, or `unknown_or_indeterminate`.",
  "",
  "### Beoord fallback impact",
  "",
  "The `Beoord` fallback is now interpreted as `+ = 10^3`, `++ = 10^4`, and `+++ = 10^5` when CFU text is missing.",
  "",
  md_table(beoord_fallback_counts),
  "",
  sprintf(
    "There are %d current clinical rows where culture support comes from `Beoord` rather than a CFU string. None of these rows are primary UTI unless the symptom rule is also compatible.",
    nrow(beoord_fallback_rows)
  ),
  "",
  "The row-level table is written to `beoord_fallback_primary_culture_support_rows.csv`.",
  "",
  "### VF-ready rows",
  "",
  sprintf(
    "- The active VF denominator is the exact key-matched selected Longcycler cohort: currently %d rows.",
    vf_total
  ),
  "- Each row is a sequenced VF participant-timepoint from `vf_pa_all.csv`.",
  "- `22_vf_build_analysis_dataset.R` joins clinical status to VF rows by `Participant_id + tp_lab`.",
  "- `Primary VF-ready UTI` is `sum(vf_analysis_ready$UTI_Status == 'UTI')` within active primary/genomics rows.",
  "- `Primary VF-ready Not_UTI` is `sum(vf_analysis_ready$UTI_Status == 'Not_UTI')` within active primary/genomics rows.",
  "- `Primary VF-ready missing status` is `sum(is.na(vf_analysis_ready$UTI_Status))` within active primary/genomics rows. Missing means the VF row has no matching clinical status key. It is not a third clinical category.",
  "",
  "### Statistical diagnostics added for interpretation",
  "",
  "The new active diagnostic layer is generated by `32_uti_not_uti_diagnostic_stats.R`. It does not change the UTI definition and does not relabel any near-miss, duplicate, unknown, or quarantined row.",
  "",
  if (nrow(diag_summary) == 0) {
    "`results/vf/uti_not_uti_diagnostic_summary.csv` was not found. Run `Rscript 32_uti_not_uti_diagnostic_stats.R` to regenerate the diagnostic layer."
  } else {
    md_table(diag_summary %>% select(metric, value, interpretation))
  },
  "",
  "Diagnostic interpretation guide:",
  "",
  if (nrow(diag_interpretation) == 0) {
    "_Diagnostic interpretation table not found._"
  } else {
    md_table(diag_interpretation %>% select(
      diagnostic, evidence_type, primary_use, what_it_can_support, what_it_cannot_support
    ))
  },
  "",
  "Diagnostic figures generated by script 32:",
  "",
  if (nrow(diag_fig_meta) == 0) {
    "_Diagnostic figure metadata not found._"
  } else {
    md_table(diag_fig_meta %>% transmute(
      figure_id,
      file_path,
      evidence_type,
      output_exists,
      interpretation_limitations
    ))
  },
  "",
  "Key reading of these diagnostics:",
  "",
  sprintf(
    "- The decision-flow and denominator plots use the selected Longcycler analytical denominator `%d = %d UTI + %d Not_UTI`; the broader `%d`-row status map appears only as source attrition/QC.",
    vf_total, vf_uti, vf_not_uti, clinical_raw_total
  ),
  "- The near-miss table/heatmap identifies rows that look clinically suspicious under legacy logic but remain `Not_UTI` because the current symptom rule is not met or is unknown.",
  "- The bootstrap, Fisher, paired, transition, and leave-one-UTI-out outputs are conservative interpretation tools for sparse UTI counts. They are useful for effect sizes, uncertainty, and stability, but not for definitive causal claims.",
  "- The duplicate-culture QC output keeps `31036 UTI-2` outside primary denominators and documents it only as a secondary duplicate comparison row.",
  "",
  "### Manual curation/exclusion policy",
  "",
  "- `UNKNOWN / UTI-? / 2517C049401-1 / UTI label 89003` is excluded from all primary analyses because the participant identity is unknown.",
  "- `31036 / UTI-2 / 24240213401-1 / UTI label 39009` is excluded from primary analyses as a duplicate dipslide culture for the same UTI suspicion; `31036 UTI-1` remains the primary representative.",
  "- Secondary duplicate rows may be used only in explicit duplicate-culture or longitudinal QC outputs, with the reuse stated.",
  "- The four failed/not-expected FASTA rows are excluded from active genomics expected-denominator checks and written to `quarantined_failed_or_not_expected_fastas.csv`.",
  "",
  "### Legacy rows",
  "",
  "- `Infection_Status_legacy` is kept only for audit. It preserves the older ASB / UTI / Negative style label.",
  sprintf(
    "- The legacy UTI count is currently %d: `sum(Infection_Status_legacy == 'UTI')`.",
    legacy_uti_total
  ),
  sprintf(
    "- The legacy UTI rows reclassified to primary Not_UTI are currently %d: `Infection_Status_legacy == 'UTI' & UTI_Status == 'Not_UTI'`.",
    nrow(legacy_reclass)
  ),
  "- Those rows stay Not_UTI unless their upstream symptom or urine-collection evidence changes enough to satisfy the current primary rule.",
  "",
  "### Participant 100011 check",
  "",
  "Participant `100011` is included as a concrete row-level check because it previously exposed a workbook/keying difference (`UTI-5` was absent before the `_CLEAN_RC` refresh).",
  "",
  md_table(participant_100011_status %>% select(
    tp_lab, Collection_Date, UTI_Label, UTI_Status, Not_UTI_subgroup,
    UTI_binary, Infection_Status_legacy, Urine_collection_method,
    culture_supports_uti, symptom_compatible_uti, status_reason
  )),
  "",
  "`100011 UTI-4` and `100011 UTI-5` are current primary `UTI` rows because both `culture_supports_uti` and `symptom_compatible_uti` are true. The other `100011` rows are current `Not_UTI` because culture support is present but the symptom rule is not compatible.",
  "",
  "### Machine-readable logic reference",
  "",
  "The count-by-count logic is also written to `uti_count_logic_reference.csv`.",
  "",
  sprintf("## Why there are only %d clinical UTIs", clinical_uti),
  "",
  sprintf(
    "The current primary UTI definition requires both culture support and catheter-aware symptom compatibility. The %d clinical UTI rows are exactly the rows satisfying both pieces of evidence. The old/legacy label has %d clinical `UTI` rows, but %d of those are now `Not_UTI` because the culture evidence is present while the symptom criterion is either not met or unknown.",
    clinical_uti, legacy_uti_total, nrow(legacy_reclass)
  ),
  "",
  md_table(legacy_reclass_counts),
  "",
  sprintf(
    "The full %d-row table is written to `legacy_uti_reclassified_not_uti.csv`. These are the main clinically suspicious rows if you expected more UTIs, but they should remain `Not_UTI` under the current criteria unless the source symptom data are corrected upstream.",
    nrow(legacy_reclass)
  ),
  "",
  sprintf("## Why %d VF rows have missing status", vf_missing_n),
  "",
  md_table(missing_counts),
  "",
  missing_status_text,
  "",
  "## Possible unlabelled UTI candidates",
  "",
  if (nrow(wgs_only) == 0) {
    "_None._"
  } else {
    md_table(wgs_only %>% select(
      Participant_id, tp_lab, Episode_ID, ST, metadata_collection_date,
      metadata_uti_label, metadata_clinical_cfu_count,
      metadata_clinical_beoord, metadata_urine_collection_method,
      metadata_population
    ))
  },
  "",
  unlabelled_text,
  "",
  "## Direct unmatched UTI-N bridge resolution",
  "",
  bridge_text,
  "",
  md_table(direct_unmatched_counts),
  "",
  "## Output files",
  "",
  sprintf("- `primary_uti_rows.csv`: the %d rows currently classified as primary UTI.", clinical_uti),
  "- `primary_clinical_manual_exclusions.csv`: clinical rows retained for audit but removed from primary analyses.",
  "- `primary_vf_manual_exclusions.csv`: VF/model rows retained for audit but removed from primary analyses.",
  "- `quarantined_failed_or_not_expected_fastas.csv`: failed/not-expected FASTA rows removed from active genomics expected-denominator checks.",
  "- `participant_100011_status_audit.csv`: all current status-map rows for participant 100011, including the rule fields behind each status.",
  sprintf("- `legacy_uti_reclassified_not_uti.csv`: the %d old UTI rows that are now Not_UTI, with the reason bucket.", nrow(legacy_reclass)),
  "- `beoord_fallback_primary_culture_support_rows.csv`: rows where `Beoord` supplies the primary culture-support evidence.",
  sprintf("- `vf_ready_missing_status_rows.csv`: the %d missing VF rows with category, metadata, and whether they can receive a primary status now.", vf_missing_n),
  sprintf("- `wgs_only_possible_unlabelled_uti_rows.csv`: the %d WGS-only possible unlabelled UTI candidates.", wgs_only_n),
  "- `clinical_status_rows_for_wgs_only_participants.csv`: existing status-map rows for WGS-only candidates, when present.",
  sprintf("- `vf_uti_n_direct_unmatched_bridge_resolution.csv`: all %d initially unmatched UTI-N VF rows and how the bridge handled them.", direct_unmatched_n),
  "- `status_map_rows_meeting_uti_criteria_but_not_labelled_uti.csv`: should be empty; written as a guardrail check.",
  "- `results/vf/uti_not_uti_diagnostic_summary.csv`: primary count and diagnostic-layer summary from script 32.",
  "- `results/audit/uti_not_uti_near_miss_rows.csv`: legacy/culture-supported near-miss rows that remain Not_UTI under the current rule.",
  "- `results/vf/uti_not_uti_bootstrap_effects.csv`: participant-bootstrap UTI minus Not_UTI VF score effect sizes and confidence intervals.",
  "- `results/vf/uti_not_uti_feature_fisher_exploratory.csv`: exploratory Fisher exact gene/module screens.",
  "- `results/vf/uti_not_uti_leave_one_uti_out.csv`: leave-one-UTI-out sparse-count stability diagnostics.",
  "- `results/vf/uti_not_uti_paired_participant_deltas.csv`: within-resident UTI versus Not_UTI score deltas where both statuses exist.",
  "- `results/vf/uti_not_uti_transition_score_tests.csv`: Not_UTI-to-UTI transition score-change summaries.",
  sprintf("- `results/vf/uti_not_uti_power_precision_context.csv`: sparse UTI precision context for n=%d active Longcycler UTI VF rows.", vf_uti),
  "- `results/audit/duplicate_culture_qc_31036.csv`: explicit duplicate-culture QC table for `31036 UTI-1` and `31036 UTI-2`."
)

write_lines(report, file.path(out_dir, "uti_status_count_explanation.md"))

message("Wrote audit report and tables to ", normalizePath(out_dir))
message("Clinical counts:")
print(clinical_counts)
message("VF counts:")
print(vf_counts)
message("Manual primary clinical exclusions:")
print(clinical_excluded_primary)
message("Quarantined failed/not-expected FASTA rows:")
print(quarantined_fastas)
message("Missing status categories:")
print(missing_counts)
message("Legacy UTI reclassification counts:")
print(legacy_reclass_counts)

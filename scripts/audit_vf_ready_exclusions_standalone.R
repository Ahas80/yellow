#!/usr/bin/env Rscript
# ==============================================================================
# audit_vf_ready_exclusions_standalone.R
# ------------------------------------------------------------------------------
# Explain which clinical episodes and FASTA assemblies are not usable in the
# canonical VF-ready dataset, and why.
#
# Standalone diagnostic only:
#   - not part of RUN_COMPLETE_ANALYSIS.sh
#   - not a numbered pipeline step
#   - does not change VF calls or analysis denominators
#
# Run from the project root:
#   Rscript scripts/audit_vf_ready_exclusions_standalone.R
#
# Outputs are written under results/qc/ for manual double-checking.
# ==============================================================================

source("00_config.R")
source("R/pipeline_qc_helpers.R")

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
  library(tidyr)
})

ensure_dir(DIR_QC)

OUT_FLOW <- file.path(DIR_QC, "vf_ready_clinical_episode_flow.csv")
OUT_EXCLUDED_EPISODES <- file.path(DIR_QC, "vf_ready_excluded_clinical_episodes.csv")
OUT_EXCLUDED_FASTAS <- file.path(DIR_QC, "vf_ready_excluded_fastas.csv")
OUT_UNLINKED_FASTAS <- file.path(DIR_QC, "unlinked_candidate_fastas_reuse_assessment.csv")
OUT_REPORT <- file.path(DIR_QC, "vf_ready_exclusion_audit.md")

collapse_values <- function(x) {
  x <- unique(as.character(x[!is.na(x) & nzchar(as.character(x))]))
  if (length(x) == 0) NA_character_ else paste(x, collapse = "; ")
}

as_bool <- function(x) {
  if (is.logical(x)) return(x)
  x <- trimws(tolower(as.character(x)))
  case_when(
    x %in% c("true", "t", "1", "yes", "y") ~ TRUE,
    x %in% c("false", "f", "0", "no", "n") ~ FALSE,
    TRUE ~ NA
  )
}

status <- read_csv(FILE_STATUS_MAP, show_col_types = FALSE) %>%
  prefer_primary_uti_status() %>%
  mutate(
    Participant_id = as.character(Participant_id),
    tp_lab = normalise_timepoint_preserve_events(tp_lab),
    Event_type = episode_event_type(tp_lab)
  )

vf_ready <- read_csv(FILE_VF_READY, show_col_types = FALSE) %>%
  prefer_primary_uti_status() %>%
  mutate(
    Participant_id = as.character(Participant_id),
    tp_lab = normalise_timepoint_preserve_events(tp_lab),
    Episode_ID = as.character(Episode_ID)
  )

canonical <- read_csv(file.path(DIR_QC, "canonical_assembly_selection.csv"), show_col_types = FALSE) %>%
  mutate(
    Participant_id = as.character(Participant_id),
    tp_lab = normalise_timepoint_preserve_events(tp_lab),
    file_exists = as_bool(file_exists),
    usable_fasta = as_bool(usable_fasta),
    QC_PASS = as_bool(QC_PASS),
    selected_canonical = as_bool(selected_canonical)
  )

bridge_file <- file.path(DIR_QC, "uricult_bridge_audit.csv")
bridge_selected <- if (file.exists(bridge_file)) {
  read_csv(bridge_file, show_col_types = FALSE) %>%
    filter(selected %in% TRUE) %>%
    transmute(
      Episode_ID = as.character(Episode_ID_clinical),
      mapped_wgs_tp_lab = normalise_timepoint_preserve_events(mapped_tp_lab),
      linked_via_uricult_bridge = TRUE,
      bridge_basis = match_basis,
      bridge_selection_reason = selection_reason
    )
} else {
  tibble(
    Episode_ID = character(),
    mapped_wgs_tp_lab = character(),
    linked_via_uricult_bridge = logical(),
    bridge_basis = character(),
    bridge_selection_reason = character()
  )
}

ready_episodes <- vf_ready %>%
  filter(!is.na(Infection_Status), !is.na(Episode_ID)) %>%
  transmute(
    Episode_ID,
    ready_tp_lab = tp_lab,
    ready_status = Infection_Status,
    retained_in_vf_ready = TRUE
  ) %>%
  distinct(Episode_ID, .keep_all = TRUE)

assembly_summary <- canonical %>%
  group_by(Participant_id, tp_lab) %>%
  summarise(
    n_assembly_rows = n(),
    n_file_exists = sum(file_exists %in% TRUE, na.rm = TRUE),
    n_qc_pass = sum(QC_PASS %in% TRUE, na.rm = TRUE),
    n_selected_canonical = sum(selected_canonical %in% TRUE, na.rm = TRUE),
    all_fastas = collapse_values(file_name),
    selected_fastas = collapse_values(file_name[selected_canonical %in% TRUE]),
    qc_fail_fastas = collapse_values(file_name[!(QC_PASS %in% TRUE)]),
    qc_fail_reasons = collapse_values(QC_REASON[!(QC_PASS %in% TRUE)]),
    selected_canonical_reason = collapse_values(canonical_reason[selected_canonical %in% TRUE]),
    .groups = "drop"
  )

episode_flow <- status %>%
  left_join(ready_episodes, by = "Episode_ID") %>%
  left_join(bridge_selected, by = "Episode_ID") %>%
  mutate(
    retained_in_vf_ready = retained_in_vf_ready %in% TRUE,
    linked_via_uricult_bridge = linked_via_uricult_bridge %in% TRUE,
    wgs_tp_lab = coalesce(mapped_wgs_tp_lab, tp_lab)
  ) %>%
  left_join(assembly_summary, by = c("Participant_id", "wgs_tp_lab" = "tp_lab")) %>%
  mutate(
    across(
      c(n_assembly_rows, n_file_exists, n_qc_pass, n_selected_canonical),
      ~ tidyr::replace_na(.x, 0L)
    ),
    exclusion_reason = case_when(
      retained_in_vf_ready ~ "Retained in VF-ready dataset",
      Participant_id %in% c("Still to be linked", "Niet te koppelen") ~
        "Untrusted participant identifier; no safe Participant_id/timepoint WGS mapping",
      n_assembly_rows == 0 ~
        "No assembly metadata/WGS link for this clinical episode",
      n_file_exists == 0 ~
        "Expected assembly metadata row exists, but FASTA file is absent on disk",
      n_qc_pass == 0 ~
        paste0("FASTA exists, but no candidate assembly passes WGS QC", if_else(!is.na(qc_fail_reasons), paste0(" (", qc_fail_reasons, ")"), "")),
      n_selected_canonical == 0 ~
        "QC-passing assembly exists, but no selected canonical assembly was produced",
      TRUE ~ "Not retained despite WGS evidence; investigate key mismatch"
    ),
    reuse_assessment = case_when(
      retained_in_vf_ready ~
        "Already usable in the canonical VF-ready dataset",
      str_detect(exclusion_reason, "Untrusted participant") ~
        "Potentially recoverable only if source records can confidently resolve the participant/timepoint identity",
      str_detect(exclusion_reason, "No assembly metadata") ~
        "Potentially recoverable if a trusted clinical-to-WGS/isolate mapping can be reconstructed",
      str_detect(exclusion_reason, "absent on disk") ~
        "Recoverable if the missing FASTA is located/restored, then metadata/QC/VF steps are rerun",
      str_detect(exclusion_reason, "no candidate assembly passes WGS QC") ~
        "Not suitable for the canonical denominator as-is; possible only after reassembly/re-QC or as labelled sensitivity analysis",
      str_detect(exclusion_reason, "no selected canonical") ~
        "Potentially recoverable by reviewing canonical selection and rerunning VF steps",
      TRUE ~ "Manual investigation required"
    )
  ) %>%
  arrange(Infection_Status, Participant_id, tp_lab)

excluded_episodes <- episode_flow %>%
  filter(!retained_in_vf_ready) %>%
  select(
    Participant_id, tp_lab, wgs_tp_lab, Event_type, Infection_Status,
    Episode_ID, Collection_Date, Batch, UTI_Label, Urine_collection_method,
    n_assembly_rows, n_file_exists, n_qc_pass, n_selected_canonical,
    all_fastas, qc_fail_fastas, qc_fail_reasons,
    exclusion_reason, reuse_assessment
  )

excluded_fasta_keys <- excluded_episodes %>%
  select(
    Participant_id, clinical_tp_lab = tp_lab, wgs_tp_lab, Infection_Status,
    clinical_event_type = Event_type,
    clinical_episode_id = Episode_ID, exclusion_reason, reuse_assessment
  )

excluded_fastas <- canonical %>%
  inner_join(excluded_fasta_keys, by = c("Participant_id", "tp_lab" = "wgs_tp_lab")) %>%
  transmute(
    Participant_id,
    clinical_tp_lab,
    wgs_tp_lab = tp_lab,
    Infection_Status,
    Event_type = clinical_event_type,
    Episode_ID = clinical_episode_id,
    Isolate_ID,
    Assembly_ID,
    Assembler,
    file_name,
    full_path,
    file_exists,
    QC_PASS,
    QC_REASON,
    selected_canonical,
    canonical_reason,
    exclusion_reason,
    reuse_assessment
  ) %>%
  arrange(Infection_Status, Participant_id, clinical_tp_lab, Assembly_ID)

unlinked <- if (file.exists(file.path(DIR_QC, "00_unexpected_assemblies.csv"))) {
  read_csv(file.path(DIR_QC, "00_unexpected_assemblies.csv"), show_col_types = FALSE)
} else if (file.exists(file.path(DIR_QC, "unlinked_unexpected_fastas.csv"))) {
  read_csv(file.path(DIR_QC, "unlinked_unexpected_fastas.csv"), show_col_types = FALSE) %>%
    filter(fasta_class == "candidate_input_fasta", fasta_source_status == "unexpected_unlinked")
} else {
  tibble()
}

unlinked_assessment <- unlinked %>%
  mutate(
    exclusion_reason = "FASTA is present on disk but does not map to the overview expected isolate universe or batch clinical CSVs",
    reuse_assessment = "Potentially recoverable only if a trusted lab crosswalk can map this FASTA to a real participant/timepoint/isolate; otherwise exclude from biological episode analyses"
  ) %>%
  select(any_of(c(
    "Isolate_ID", "Assembler", "file_name", "full_path", "in_overview", "in_batch_csv",
    "fasta_class", "fasta_source_status", "exclusion_reason", "reuse_assessment"
  ))) %>%
  arrange(Isolate_ID, Assembler, file_name)

write_csv(episode_flow, OUT_FLOW)
write_csv(excluded_episodes, OUT_EXCLUDED_EPISODES)
write_csv(excluded_fastas, OUT_EXCLUDED_FASTAS)
write_csv(unlinked_assessment, OUT_UNLINKED_FASTAS)

flow_summary <- episode_flow %>%
  count(Infection_Status, retained_in_vf_ready, name = "n") %>%
  arrange(Infection_Status, desc(retained_in_vf_ready))

reason_summary <- excluded_episodes %>%
  count(Infection_Status, Event_type, exclusion_reason, name = "n") %>%
  arrange(Infection_Status, Event_type, exclusion_reason)

report <- c(
  "# VF-ready FASTA/episode exclusion audit",
  "",
  sprintf("Generated: %s", format(Sys.time())),
  "",
  "## Current clinical episode flow",
  "",
  paste(capture.output(print(flow_summary, n = Inf)), collapse = "\n"),
  "",
  "## Excluded clinical episodes by reason",
  "",
  paste(capture.output(print(reason_summary, n = Inf, width = Inf)), collapse = "\n"),
  "",
  "## Output files",
  "",
  sprintf("- Clinical episode flow: `%s`", OUT_FLOW),
  sprintf("- Excluded clinical episodes: `%s`", OUT_EXCLUDED_EPISODES),
  sprintf("- FASTA rows attached to excluded episodes: `%s`", OUT_EXCLUDED_FASTAS),
  sprintf("- Unlinked candidate FASTA reuse assessment: `%s`", OUT_UNLINKED_FASTAS),
  "",
  "## Interpretation",
  "",
  "The main canonical VF-ready denominator excludes episodes when no trusted clinical-to-WGS link exists, when a linked FASTA fails WGS QC, or when the participant identifier is explicitly unresolved. FASTAs that are merely Flye/LongCycler alternatives are not biological losses; the canonical selector keeps one QC-passing assembly per participant-timepoint for denominator stability."
)
writeLines(report, OUT_REPORT)

msg("Wrote clinical episode flow: %s", OUT_FLOW)
msg("Wrote excluded clinical episodes: %s", OUT_EXCLUDED_EPISODES)
msg("Wrote excluded FASTAs: %s", OUT_EXCLUDED_FASTAS)
msg("Wrote unlinked FASTA assessment: %s", OUT_UNLINKED_FASTAS)
msg("Wrote markdown report: %s", OUT_REPORT)

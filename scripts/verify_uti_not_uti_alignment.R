#!/usr/bin/env Rscript
# ==============================================================================
# Verify UTI vs Not_UTI alignment
# ------------------------------------------------------------------------------
# This is a lightweight post-run gate for the primary clinical-status redesign.
# It checks that current outputs use UTI_Status/UTI_binary and that stale
# ASB-vs-UTI generated outputs are not left in primary result/plot folders.
# ==============================================================================

source("00_config.R")

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
  library(tibble)
})

msg("Starting scripts/verify_uti_not_uti_alignment.R")

checks <- tibble(check = character(), status = character(), detail = character())
add_check <- function(check, ok, detail = "") {
  checks <<- bind_rows(checks, tibble(
    check = check,
    status = if (isTRUE(ok)) "PASS" else "FAIL",
    detail = detail
  ))
}

read_optional <- function(path) {
  if (!file.exists(path)) return(NULL)
  read_csv(path, show_col_types = FALSE)
}

status <- read_optional(FILE_STATUS_MAP)
poster <- read_optional(FILE_STATUS_MAP_POSTER)
analysis_cohort <- read_optional(FILE_ANALYSIS_CLINICAL_COHORT)
vf_ready <- read_optional(FILE_VF_READY)
table02 <- read_optional(file.path(DIR_RESULTS, "summary", "table_02_clinical_status_counts.csv"))
table10 <- read_optional(file.path(DIR_RESULTS, "summary", "table_10_not_uti_uti_transition_cases.csv"))
transition_index <- read_optional(file.path(DIR_VF, "vf_transition_case_index.csv"))
canonical_transitions <- read_optional(file.path(DIR_RESULTS, "longitudinal", "longcycler_transitions.csv"))

active_longcycler_keys <- NULL
add_check("selected Longcycler analysis cohort exists", !is.null(analysis_cohort), FILE_ANALYSIS_CLINICAL_COHORT)
if (!is.null(analysis_cohort)) {
  analysis_cohort <- analysis_cohort %>%
    mutate(Participant_id = as.character(.data$Participant_id),
           tp_lab = normalise_timepoint_preserve_events(.data$tp_lab))
  active_longcycler_keys <- analysis_cohort %>% distinct(.data$Participant_id, .data$tp_lab)
  add_check("selected Longcycler cohort keys are unique",
            nrow(active_longcycler_keys) == nrow(analysis_cohort),
            sprintf("rows=%d unique_keys=%d", nrow(analysis_cohort), nrow(active_longcycler_keys)))
  add_check("selected Longcycler cohort contract is 532 rows, 161 residents, 16 UTI, 516 Not_UTI",
            nrow(analysis_cohort) == 532L && n_distinct(analysis_cohort$Participant_id) == 161L &&
              sum(analysis_cohort$UTI_Status == "UTI", na.rm = TRUE) == 16L &&
              sum(analysis_cohort$UTI_Status == "Not_UTI", na.rm = TRUE) == 516L,
            sprintf("n=%d residents=%d UTI=%d Not_UTI=%d", nrow(analysis_cohort),
                    n_distinct(analysis_cohort$Participant_id),
                    sum(analysis_cohort$UTI_Status == "UTI", na.rm = TRUE),
                    sum(analysis_cohort$UTI_Status == "Not_UTI", na.rm = TRUE)))
}

add_check("status_map exists", !is.null(status), FILE_STATUS_MAP)
if (!is.null(status)) {
  status <- apply_manual_sample_curation(status, context = "verify_status")
  status_primary <- filter_primary_analysis(status)
  add_check("status_map has primary UTI fields",
            all(c("UTI_Status", "UTI_binary") %in% names(status)),
            paste(setdiff(c("UTI_Status", "UTI_binary"), names(status)), collapse = ", "))
  add_check("status_map has manual curation fields",
            all(c("analysis_include_primary", "analysis_exclusion_reason", "duplicate_role",
                  "genomics_expected_include") %in% names(status)),
            paste(setdiff(c("analysis_include_primary", "analysis_exclusion_reason",
                            "duplicate_role", "genomics_expected_include"), names(status)), collapse = ", "))
  add_check("status_map UTI_Status values valid",
            all(stats::na.omit(status$UTI_Status) %in% c("UTI", "Not_UTI")),
            paste(sort(unique(status$UTI_Status)), collapse = ", "))
  add_check("status_map UTI_binary values valid",
            all(stats::na.omit(status$UTI_binary) %in% c(0, 1)),
            paste(sort(unique(status$UTI_binary)), collapse = ", "))
  if ("Episode_ID" %in% names(status)) {
    dup_epi <- status %>% count(Episode_ID, name = "n") %>% filter(!is.na(Episode_ID), n > 1)
    add_check("no duplicate Episode_ID", nrow(dup_epi) == 0, sprintf("%d duplicated Episode_ID(s)", nrow(dup_epi)))
  }
  if (all(c("UTI_Status", "culture_supports_uti", "symptom_compatible_uti") %in% names(status))) {
    bad_culture <- status_primary %>% filter(UTI_Status == "UTI", !(culture_supports_uti %in% TRUE))
    bad_symptom <- status_primary %>% filter(UTI_Status == "UTI", !(symptom_compatible_uti %in% TRUE))
    add_check("all UTI rows have culture support", nrow(bad_culture) == 0, sprintf("%d bad row(s)", nrow(bad_culture)))
    add_check("all UTI rows have symptom support", nrow(bad_symptom) == 0, sprintf("%d bad row(s)", nrow(bad_symptom)))
  }
  add_check("source clinical denominator retained for attrition/QC only",
            nrow(status_primary) == 583 &&
              sum(status_primary$UTI_Status == "UTI", na.rm = TRUE) == 18 &&
              sum(status_primary$UTI_Status == "Not_UTI", na.rm = TRUE) == 565,
            sprintf("n=%d UTI=%d Not_UTI=%d",
                    nrow(status_primary),
                    sum(status_primary$UTI_Status == "UTI", na.rm = TRUE),
                    sum(status_primary$UTI_Status == "Not_UTI", na.rm = TRUE)))
  unknown_excluded <- status %>%
    filter(Participant_id == "UNKNOWN", tp_lab == "UTI-?",
           as.character(UTI_Label) == "89003",
           !(analysis_include_primary %in% TRUE))
  add_check("UNKNOWN UTI-? 89003 excluded from primary status denominator",
            nrow(unknown_excluded) == 1,
            sprintf("%d matching excluded row(s)", nrow(unknown_excluded)))
  dup_excluded <- status %>%
    filter(Participant_id == "31036", tp_lab == "UTI-2",
           as.character(UTI_Label) == "39009",
           !(analysis_include_primary %in% TRUE),
           duplicate_role == "secondary_duplicate")
  add_check("31036 UTI-2 duplicate excluded from primary status denominator",
            nrow(dup_excluded) == 1,
            sprintf("%d matching excluded row(s)", nrow(dup_excluded)))
}

poster_ok <- !is.null(poster) &&
  all(c("UTI_Status", "UTI_binary") %in% names(poster)) &&
  !status_map_is_stale(FILE_STATUS_MAP_POSTER, FILE_STATUS_MAP)
add_check("status_map_with_poster_tp current and primary", poster_ok, FILE_STATUS_MAP_POSTER)

if (!is.null(analysis_cohort) && !is.null(table02) && all(c("dataset_layer", "Infection_Status", "n_episodes") %in% names(table02))) {
  clinical_uti <- table02 %>% filter(dataset_layer == "selected_longcycler_analysis_cohort", Infection_Status == "UTI") %>% pull(n_episodes)
  expected_uti <- sum(analysis_cohort$UTI_Status == "UTI", na.rm = TRUE)
  add_check("summary table selected-cohort UTI count matches analysis cohort",
            length(clinical_uti) == 1 && clinical_uti == expected_uti,
            sprintf("table=%s expected=%d", paste(clinical_uti, collapse = ","), expected_uti))
}

if (!is.null(vf_ready)) {
  vf_ready <- apply_manual_sample_curation(vf_ready, context = "verify_vf_ready")
  vf_primary <- filter_primary_genomics(vf_ready)
  add_check("vf_analysis_ready has primary UTI fields",
            all(c("UTI_Status", "UTI_binary") %in% names(vf_ready)),
            paste(setdiff(c("UTI_Status", "UTI_binary"), names(vf_ready)), collapse = ", "))
  add_check("vf_analysis_ready has manual curation fields",
            all(c("analysis_include_primary", "analysis_exclusion_reason", "duplicate_role",
                  "genomics_expected_include") %in% names(vf_ready)),
            paste(setdiff(c("analysis_include_primary", "analysis_exclusion_reason",
                            "duplicate_role", "genomics_expected_include"), names(vf_ready)), collapse = ", "))
  if (!is.null(active_longcycler_keys) && !is.null(analysis_cohort)) {
    vf_primary <- vf_primary %>%
      mutate(
        Participant_id = as.character(.data$Participant_id),
        tp_lab = normalise_timepoint_preserve_events(.data$tp_lab)
      )
    expected_active_status <- analysis_cohort
    expected_active_total <- nrow(active_longcycler_keys)
    expected_active_uti <- sum(expected_active_status$UTI_Status == "UTI", na.rm = TRUE)
    expected_active_not_uti <- sum(expected_active_status$UTI_Status == "Not_UTI", na.rm = TRUE)
    vf_keys <- vf_primary %>% distinct(.data$Participant_id, .data$tp_lab)
    key_sets_match <- nrow(anti_join(vf_keys, active_longcycler_keys,
                                    by = c("Participant_id", "tp_lab"))) == 0 &&
      nrow(anti_join(active_longcycler_keys, vf_keys,
                    by = c("Participant_id", "tp_lab"))) == 0
    add_check("primary VF/model keys equal the active selected QC-pass Longcycler manifest",
              key_sets_match && nrow(vf_primary) == expected_active_total,
              sprintf("VF n=%d; active Longcycler n=%d", nrow(vf_primary), expected_active_total))
    add_check("active Longcycler VF/model status counts match the primary status map",
              sum(vf_primary$UTI_Status == "UTI", na.rm = TRUE) == expected_active_uti &&
                sum(vf_primary$UTI_Status == "Not_UTI", na.rm = TRUE) == expected_active_not_uti,
              sprintf("observed UTI=%d Not_UTI=%d; expected UTI=%d Not_UTI=%d",
                      sum(vf_primary$UTI_Status == "UTI", na.rm = TRUE),
                      sum(vf_primary$UTI_Status == "Not_UTI", na.rm = TRUE),
                      expected_active_uti, expected_active_not_uti))
  }
}

quarantine <- read_optional(FILE_QUARANTINED_FASTA_EXPECTATIONS)
add_check("quarantined failed/not-expected FASTA report exists",
          !is.null(quarantine), FILE_QUARANTINED_FASTA_EXPECTATIONS)
if (!is.null(quarantine)) {
  expected_quarantine <- tibble(
    Participant_id = c("110059", "110060", "20043", "30010"),
    tp_lab = c("T0", "T3", "UTI-4", "T6"),
    Isolate_ID = c("24210082601-2", "2508C147101-1", "2443C048901-1", "2543C105201-1")
  )
  quarantine_norm <- quarantine %>%
    mutate(Participant_id = as.character(Participant_id),
           tp_lab = normalise_timepoint_preserve_events(tp_lab),
           Isolate_ID = normalise_isolate_id(Isolate_ID))
  add_check("four failed/not-expected FASTA rows are quarantined",
            nrow(semi_join(expected_quarantine, quarantine_norm,
                           by = c("Participant_id", "tp_lab", "Isolate_ID"))) == 4,
            sprintf("quarantine rows=%d", nrow(quarantine_norm)))
}

if (!is.null(table10) && nrow(table10) > 0 && all(c("from_status", "to_status") %in% names(table10))) {
  bad_t10 <- table10 %>% filter(!(from_status == "Not_UTI" & to_status == "UTI"))
  add_check("table_10 contains only Not_UTI->UTI transitions",
            nrow(bad_t10) == 0,
            sprintf("%d non-target row(s)", nrow(bad_t10)))
  add_check("table_10 contains all 9 selected Longcycler Not_UTI->UTI transitions with paired evidence",
            nrow(table10) == 9L && all(table10$has_vf_pair %in% TRUE) && all(!is.na(table10$SNPs)),
            sprintf("n=%d linked=%d direct_snp=%d", nrow(table10),
                    sum(table10$has_vf_pair %in% TRUE), sum(!is.na(table10$SNPs))))
}
if (!is.null(transition_index)) {
  add_check("transition index contains 371 selected adjacent pairs and 9 Not_UTI->UTI cases",
            nrow(transition_index) == 371L &&
              sum(transition_index$is_not_uti_to_uti %in% TRUE) == 9L &&
              all(transition_index$has_vf_pair %in% TRUE),
            sprintf("n=%d target=%d linked=%d", nrow(transition_index),
                    sum(transition_index$is_not_uti_to_uti %in% TRUE),
                    sum(transition_index$has_vf_pair %in% TRUE)))
}
add_check("canonical Longcycler transition export exists", !is.null(canonical_transitions),
          file.path(DIR_RESULTS, "longitudinal", "longcycler_transitions.csv"))
if (!is.null(canonical_transitions)) {
  canonical_keys <- canonical_transitions %>%
    transmute(key = paste(as.character(.data$Participant_id),
                          normalise_timepoint_preserve_events(.data$tp_from),
                          normalise_timepoint_preserve_events(.data$tp_to), sep = "|"))
  add_check("canonical transition export satisfies 371/9/direct-SNP contract",
            nrow(canonical_transitions) == 371L &&
              sum(canonical_transitions$status_from == "Not_UTI" &
                    canonical_transitions$status_to == "UTI", na.rm = TRUE) == 9L &&
              all(!is.na(canonical_transitions$TotalSNPs)) &&
              !anyDuplicated(canonical_keys$key),
            sprintf("n=%d target=%d direct_snp=%d", nrow(canonical_transitions),
                    sum(canonical_transitions$status_from == "Not_UTI" &
                          canonical_transitions$status_to == "UTI", na.rm = TRUE),
                    sum(!is.na(canonical_transitions$TotalSNPs))))
  if (!is.null(transition_index)) {
    index_keys <- transition_index %>%
      transmute(key = paste(as.character(.data$Participant_id),
                            normalise_timepoint_preserve_events(.data$from_tp),
                            normalise_timepoint_preserve_events(.data$to_tp), sep = "|"))
    add_check("transition case index keys equal canonical transition keys",
              setequal(index_keys$key, canonical_keys$key),
              sprintf("index=%d canonical=%d", nrow(index_keys), nrow(canonical_keys)))
  }
}

current_legacy_outputs <- c(
  list.files(DIR_RESULTS, recursive = TRUE, full.names = TRUE, pattern = "(UTI_vs_ASB|ASB_vs_UTI|asb_uti)", ignore.case = TRUE),
  list.files(DIR_PLOTS, recursive = TRUE, full.names = TRUE, pattern = "(UTI_vs_ASB|ASB_vs_UTI|asb_uti)", ignore.case = TRUE)
)
current_legacy_outputs <- current_legacy_outputs[
  !str_detect(current_legacy_outputs, regex("/legacy/", ignore_case = TRUE))
]
add_check("no ASB-vs-UTI generated outputs in current result/plot folders",
          length(current_legacy_outputs) == 0,
          paste(current_legacy_outputs, collapse = "; "))

out_csv <- file.path(DIR_QC, "uti_not_uti_alignment_checks.csv")
out_txt <- file.path(DIR_QC, "uti_not_uti_alignment_checks.txt")
write_csv(checks, out_csv)
writeLines(c(
  "UTI vs Not_UTI alignment checks",
  sprintf("Generated: %s", format(Sys.time())),
  "",
  apply(checks, 1, function(x) sprintf("[%s] %s - %s", x[["status"]], x[["check"]], x[["detail"]]))
), out_txt)

failures <- checks %>% filter(status == "FAIL")
msg("Alignment checks written to %s", out_csv)
if (nrow(failures) > 0) {
  print(failures)
  stop(sprintf("%d UTI/Not_UTI alignment check(s) failed.", nrow(failures)))
}
msg("✓ All UTI/Not_UTI alignment checks passed.")

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
vf_ready <- read_optional(FILE_VF_READY)
table02 <- read_optional(file.path(DIR_RESULTS, "summary", "table_02_clinical_status_counts.csv"))
table10 <- read_optional(file.path(DIR_RESULTS, "summary", "table_10_not_uti_uti_transition_cases.csv"))

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
  add_check("primary clinical denominator expected after manual exclusions",
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

if (!is.null(status) && !is.null(table02) && all(c("dataset_layer", "Infection_Status", "n_episodes") %in% names(table02))) {
  clinical_uti <- table02 %>% filter(dataset_layer == "all_batch_clinical", Infection_Status == "UTI") %>% pull(n_episodes)
  poster_uti <- table02 %>% filter(dataset_layer == "ordered_poster_clinical", Infection_Status == "UTI") %>% pull(n_episodes)
  expected_uti <- sum(filter_primary_analysis(status)$UTI_Status == "UTI", na.rm = TRUE)
  add_check("summary table all-batch UTI count matches status_map",
            length(clinical_uti) == 1 && clinical_uti == expected_uti,
            sprintf("table=%s expected=%d", paste(clinical_uti, collapse = ","), expected_uti))
  add_check("summary table ordered-poster UTI count is not legacy",
            length(poster_uti) == 0 || poster_uti == expected_uti,
            sprintf("poster=%s expected=%d", paste(poster_uti, collapse = ","), expected_uti))
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
  add_check("primary VF/model denominator expected after manual exclusions",
            nrow(vf_primary) == 556 &&
              sum(vf_primary$UTI_Status == "UTI", na.rm = TRUE) == 17 &&
              sum(vf_primary$UTI_Status == "Not_UTI", na.rm = TRUE) == 539,
            sprintf("n=%d UTI=%d Not_UTI=%d",
                    nrow(vf_primary),
                    sum(vf_primary$UTI_Status == "UTI", na.rm = TRUE),
                    sum(vf_primary$UTI_Status == "Not_UTI", na.rm = TRUE)))
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

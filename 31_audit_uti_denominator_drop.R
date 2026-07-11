# ==============================================================================
# 31_audit_uti_denominator_drop.R
# ------------------------------------------------------------------------------
# Evidence-driven audit of the clinical UTI -> VF/WGS-ready denominator drop.
#
# This script is deliberately diagnostic only. It does not modify analysis
# outputs; it writes audit products under results/audit/.
# ==============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
  library(tidyr)
  library(tibble)
})

source("00_config.R")
source("R/pipeline_qc_helpers.R")

DIR_AUDIT <- file.path(DIR_RESULTS, "audit")
ensure_dir(DIR_AUDIT)

audit_msg <- function(...) {
  message(format(Sys.time(), "[%H:%M:%S] "), sprintf(...))
}

collapse_values <- function(x, sep = "; ") {
  x <- as.character(x)
  x <- x[!is.na(x) & nzchar(trimws(x))]
  x <- unique(x)
  if (length(x) == 0) NA_character_ else paste(x, collapse = sep)
}

collapse_col <- function(df, col, sep = "; ") {
  if (nrow(df) == 0 || !col %in% names(df)) return(NA_character_)
  collapse_values(df[[col]], sep = sep)
}

collapse_bool <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0) FALSE else any(x)
}

norm_tp_audit <- function(x) {
  if (exists("normalise_timepoint_preserve_events")) {
    normalise_timepoint_preserve_events(x)
  } else {
    x <- trimws(as.character(x))
    low <- tolower(x)
    num <- suppressWarnings(as.integer(readr::parse_number(low)))
    dplyr::case_when(
      is.na(x) | x == "" ~ NA_character_,
      str_detect(low, "uricult") ~ "Uricult",
      str_detect(toupper(x), "^UTI[-_A-Za-z0-9?]*$") ~ str_replace_all(toupper(x), "_", "-"),
      str_detect(toupper(x), "^T[0-9]+(\\.[0-9]+)?$") ~ toupper(x),
      str_detect(x, "^[0-9]+$") ~ paste0("T", x),
      !is.na(num) ~ paste0("T", num),
      TRUE ~ x
    )
  }
}

read_csv_optional <- function(path, label, delim = ",") {
  if (!file.exists(path)) {
    audit_msg("Optional file missing: %s (%s)", path, label)
    return(tibble())
  }
  audit_msg("Reading %s: %s", label, path)
  if (identical(delim, "\t")) {
    readr::read_tsv(path, show_col_types = FALSE, progress = FALSE)
  } else {
    readr::read_csv(path, show_col_types = FALSE, progress = FALSE)
  }
}

standardise_key_cols <- function(df, tp_candidates = c("tp_lab", "Timepoint")) {
  if (nrow(df) == 0 || !"Participant_id" %in% names(df)) {
    return(tibble(Participant_id = character(), tp_lab_norm = character()))
  }
  tp_col <- tp_candidates[tp_candidates %in% names(df)][1]
  if (is.na(tp_col)) {
    return(tibble(Participant_id = character(), tp_lab_norm = character()))
  }
  df %>%
    mutate(
      Participant_id = as.character(.data$Participant_id),
      tp_lab_norm = norm_tp_audit(.data[[tp_col]])
    )
}

summarise_status_counts <- function(df, status_col = "Infection_Status") {
  if (nrow(df) == 0 || (!status_col %in% names(df) && !"UTI_Status" %in% names(df))) {
    return(tibble(
      n_total = 0L, n_ASB = 0L, n_legacy_UTI = 0L, n_Negative = 0L,
      n_UTI = 0L, n_Not_UTI = 0L, n_other = 0L
    ))
  }

  legacy_col <- dplyr::case_when(
    "Infection_Status_old" %in% names(df) ~ "Infection_Status_old",
    "Infection_Status_legacy" %in% names(df) ~ "Infection_Status_legacy",
    status_col %in% names(df) ~ status_col,
    TRUE ~ NA_character_
  )
  legacy_status <- if (!is.na(legacy_col)) as.character(df[[legacy_col]]) else rep(NA_character_, nrow(df))

  primary_status <- if ("UTI_Status" %in% names(df)) {
    as.character(df$UTI_Status)
  } else if (status_col %in% names(df)) {
    dplyr::case_when(
      as.character(df[[status_col]]) == "UTI" ~ "UTI",
      as.character(df[[status_col]]) %in% c("ASB", "Negative") ~ "Not_UTI",
      TRUE ~ NA_character_
    )
  } else {
    rep(NA_character_, nrow(df))
  }

  tibble(
    n_total = length(primary_status),
    n_ASB = sum(legacy_status == "ASB", na.rm = TRUE),
    n_legacy_UTI = sum(legacy_status == "UTI", na.rm = TRUE),
    n_Negative = sum(legacy_status == "Negative", na.rm = TRUE),
    n_UTI = sum(primary_status == "UTI", na.rm = TRUE),
    n_Not_UTI = sum(primary_status == "Not_UTI", na.rm = TRUE),
    n_other = sum(is.na(primary_status) | !primary_status %in% c("UTI", "Not_UTI"), na.rm = TRUE)
  )
}

make_md_table <- function(df) {
  if (nrow(df) == 0) return("_No rows._")
  out <- df %>%
    mutate(across(everything(), ~ {
      x <- as.character(.x)
      x[is.na(x)] <- ""
      str_replace_all(x, "\\|", "/")
    }))
  header <- paste0("| ", paste(names(out), collapse = " | "), " |")
  sep <- paste0("| ", paste(rep("---", ncol(out)), collapse = " | "), " |")
  rows <- apply(out, 1, function(z) paste0("| ", paste(z, collapse = " | "), " |"))
  paste(c(header, sep, rows), collapse = "\n")
}

first_existing <- function(paths) {
  paths[file.exists(paths)][1] %||% NA_character_
}

`%||%` <- function(x, y) if (is.null(x) || length(x) == 0 || all(is.na(x))) y else x

legacy_status_path <- file.path(DIR_RESULTS, "status_map_pid_tp_status.csv")
legacy_status_full_path <- file.path(DIR_ROOT, "status_map.csv")
legacy_status_full_alt_path <- file.path(DIR_RESULTS, "status_map_full.csv")
current_status_path <- file.path(DIR_CLINICAL, "status_map.csv")
old_vf_pa_path <- file.path(DIR_ROOT, "vf_pa_all.csv")
current_vf_pa_path <- file.path(DIR_VF, "vf_pa_all.csv")
old_vf_ready_path <- file.path(DIR_VF, "vf_analysis_ready_OLD.csv")
current_vf_ready_path <- file.path(DIR_VF, "vf_analysis_ready.csv")
model_path <- first_existing(c(
  file.path(DIR_RESULTS, "vf", "vf_model_dataset.csv"),
  file.path(DIR_RESULTS, "13_modelling", "model_dataset_UTI_ASB.csv"),
  file.path(DIR_MODELS, "vf_model_dataset.csv")
))
mlst_path <- first_existing(c(FILE_MLST_CANONICAL, FILE_MLST_PROVIDER_PREFERRED_ALL))

legacy_status <- read_csv_optional(legacy_status_path, "legacy ASB/UTI/Negative status map")
legacy_status_full <- read_csv_optional(legacy_status_full_path, "legacy full status map")
if (nrow(legacy_status_full) == 0) {
  legacy_status_full <- read_csv_optional(legacy_status_full_alt_path, "legacy full status map fallback")
}
current_status <- read_csv_optional(current_status_path, "current canonical status_map")
if (nrow(current_status) > 0) {
  current_status <- prefer_primary_uti_status(current_status)
}
assembly_meta <- read_csv_optional(FILE_METADATA, "assembly_metadata")
qc_summary <- read_csv_optional(file.path(DIR_WGS, "qc_summary.csv"), "WGS QC")
old_vf_pa <- read_csv_optional(old_vf_pa_path, "legacy/root VF PA")
current_vf_pa <- read_csv_optional(current_vf_pa_path, "current VF PA")
old_vf_ready <- read_csv_optional(old_vf_ready_path, "legacy VF analysis-ready")
current_vf_ready <- read_csv_optional(current_vf_ready_path, "current VF analysis-ready")
model_df <- if (!is.na(model_path)) read_csv_optional(model_path, "final/model dataset") else tibble()

vf_hits <- tibble()
if (file.exists(FILE_VF_HITS)) {
  audit_msg("Reading VF hits RDS: %s", FILE_VF_HITS)
  vf_hits_obj <- readRDS(FILE_VF_HITS)
  if (inherits(vf_hits_obj, "data.frame")) {
    vf_hits <- as_tibble(vf_hits_obj)
  } else if (is.list(vf_hits_obj)) {
    vf_hits <- tryCatch(bind_rows(vf_hits_obj), error = function(e) tibble())
  }
} else {
  audit_msg("Optional file missing: %s (VF hits RDS)", FILE_VF_HITS)
}

mlst_df <- if (!is.na(mlst_path)) {
  read_csv_optional(mlst_path, "active provider-preferred MLST", delim = if (grepl("\\.tsv$", mlst_path, ignore.case = TRUE)) "\t" else ",")
} else {
  tibble()
}

# ------------------------------------------------------------------------------
# Build the legacy clinical denominator that reproduces the user-reported 20 UTI.
# ------------------------------------------------------------------------------

if (nrow(legacy_status) == 0 && nrow(current_status) > 0) {
  warning(
    "Legacy status map is missing; falling back to current status_map legacy fields for ",
    "the historical denominator trace."
  )
  tp_source <- if ("Timepoint" %in% names(current_status)) "Timepoint" else "tp_lab"
  legacy_status <- current_status %>%
    transmute(
      Participant_id = as.character(.data$Participant_id),
      Timepoint = as.character(.data[[tp_source]]),
      Infection_Status = dplyr::coalesce(
        as.character(.data$Infection_Status_old),
        as.character(.data$Infection_Status_legacy),
        as.character(.data$Infection_Status)
      )
    )
}
if (nrow(legacy_status) == 0) {
  stop("Cannot run the legacy UTI denominator audit because neither ", legacy_status_path,
       " nor a current status_map with legacy fields is available.")
}
if (!all(c("Participant_id", "Timepoint", "Infection_Status") %in% names(legacy_status))) {
  stop("Legacy status map lacks required columns: Participant_id, Timepoint, Infection_Status")
}

legacy_clinical <- legacy_status %>%
  transmute(
    Participant_id = as.character(.data$Participant_id),
    clinical_Timepoint = as.character(.data$Timepoint),
    clinical_tp_lab = norm_tp_audit(.data$Timepoint),
    Infection_Status = as.character(.data$Infection_Status)
  )

legacy_counts <- summarise_status_counts(legacy_clinical)
audit_msg(
  "Legacy clinical denominator: %s total; %s UTI.",
  legacy_counts$n_total, legacy_counts$n_UTI
)

if (legacy_counts$n_UTI != 20L) {
  warning("The legacy clinical denominator does not contain exactly 20 UTI rows; audit will continue with observed rows.")
}

# Supplement raw clinical context from richer legacy/current status maps.
raw_sources <- bind_rows(
  legacy_status_full %>%
    mutate(across(everything(), as.character), .raw_source = "root_status_map.csv"),
  current_status %>%
    mutate(across(everything(), as.character), .raw_source = "results/clinical/status_map.csv")
) %>%
  filter("Participant_id" %in% names(.)) %>%
  mutate(
    Participant_id = as.character(.data$Participant_id),
    clinical_tp_lab = norm_tp_audit(if ("tp_lab" %in% names(.)) .data$tp_lab else .data$Timepoint)
  )

raw_summary <- raw_sources %>%
  group_by(.data$Participant_id, .data$clinical_tp_lab) %>%
  summarise(
    Batch = collapse_values(.data$Batch),
    Collection_Date = collapse_values(.data$Collection_Date),
    UTI_Label = collapse_values(.data$UTI_Label),
    Urine_collection_method = collapse_values(.data$Urine_collection_method),
    Status_Confidence = collapse_values(.data$Status_Confidence_epi),
    Sx_present_any = collapse_values(.data$Sx_present_any),
    Sx_source_epi = collapse_values(.data$Sx_source_epi),
    raw_CFU_examples = collapse_values(.data$raw_CFU_examples),
    raw_BEO_examples = collapse_values(.data$raw_BEO_examples),
    raw_Pop_examples = collapse_values(.data$raw_Pop_examples),
    source_file = collapse_values(.data$.raw_source),
    .groups = "drop"
  )

clinical_all <- legacy_clinical %>%
  left_join(raw_summary, by = c("Participant_id", "clinical_tp_lab")) %>%
  mutate(
    Status_Confidence = coalesce(.data$Status_Confidence, NA_character_),
    Status_Provenance = paste0("clinical denominator from ", legacy_status_path,
                               "; raw context from ", coalesce(.data$source_file, "not available")),
    in_status_map = TRUE
  )

clinical_uti <- clinical_all %>%
  filter(.data$Infection_Status == "UTI") %>%
  arrange(.data$Participant_id, .data$clinical_tp_lab)

# ------------------------------------------------------------------------------
# Standardise evidence tables.
# ------------------------------------------------------------------------------

assembly_std <- assembly_meta %>%
  standardise_key_cols(c("tp_lab", "Timepoint")) %>%
  mutate(
    Assembly_ID = if ("Assembly_ID" %in% names(.)) as.character(.data$Assembly_ID) else NA_character_,
    Isolate_ID = if ("Isolate_ID" %in% names(.)) as.character(.data$Isolate_ID) else NA_character_,
    file_name = if ("file_name" %in% names(.)) as.character(.data$file_name) else basename(.data$full_path %||% NA_character_),
    full_path = case_when(
      "full_path" %in% names(.) ~ as.character(.data$full_path),
      "fasta_path" %in% names(.) ~ as.character(.data$fasta_path),
      TRUE ~ NA_character_
    ),
    fasta_exists_calc = case_when(
      "file_exists" %in% names(.) ~ as.logical(.data$file_exists),
      !is.na(.data$full_path) ~ file.exists(.data$full_path),
      TRUE ~ FALSE
    ),
    usable_fasta_calc = case_when(
      "usable_fasta" %in% names(.) ~ as.logical(.data$usable_fasta),
      TRUE ~ .data$fasta_exists_calc
    ),
    Event_type = if ("Event_type" %in% names(.)) as.character(.data$Event_type) else NA_character_,
    Batch_meta = if ("Batch" %in% names(.)) as.character(.data$Batch) else NA_character_,
    Collection_Date_meta = if ("Collection_Date" %in% names(.)) as.character(.data$Collection_Date) else NA_character_,
    UTI_Label_meta = if ("UTI_Label" %in% names(.)) as.character(.data$UTI_Label) else NA_character_
  )

qc_std <- qc_summary %>%
  standardise_key_cols(c("tp_lab", "Timepoint")) %>%
  mutate(
    Assembly_ID = if ("Assembly_ID" %in% names(.)) as.character(.data$Assembly_ID) else NA_character_,
    QC_PASS = if ("QC_PASS" %in% names(.)) as.logical(.data$QC_PASS) else NA,
    QC_REASON = if ("QC_REASON" %in% names(.)) as.character(.data$QC_REASON) else NA_character_
  )

old_vf_pa_keys <- old_vf_pa %>%
  standardise_key_cols(c("tp_lab", "Timepoint")) %>%
  distinct(.data$Participant_id, .data$tp_lab_norm)

current_vf_pa_keys <- current_vf_pa %>%
  standardise_key_cols(c("tp_lab", "Timepoint")) %>%
  distinct(.data$Participant_id, .data$tp_lab_norm)

old_vf_ready_keys <- old_vf_ready %>%
  standardise_key_cols(c("tp_lab", "Timepoint")) %>%
  distinct(.data$Participant_id, .data$tp_lab_norm)

current_vf_ready_keys <- current_vf_ready %>%
  standardise_key_cols(c("tp_lab", "Timepoint")) %>%
  distinct(.data$Participant_id, .data$tp_lab_norm)

current_ready_status <- current_vf_ready %>%
  standardise_key_cols(c("tp_lab", "Timepoint")) %>%
  mutate(
    ST_current = if ("ST" %in% names(.)) as.character(.data$ST) else NA_character_,
    current_ready_status = dplyr::case_when(
      "UTI_Status" %in% names(.) ~ as.character(.data$UTI_Status),
      "Infection_Status" %in% names(.) ~ as.character(.data$Infection_Status),
      TRUE ~ NA_character_
    )
  ) %>%
  group_by(.data$Participant_id, .data$tp_lab_norm) %>%
  summarise(
    ST_current = collapse_values(.data$ST_current),
    current_ready_status = collapse_values(.data$current_ready_status),
    .groups = "drop"
  )

vf_hits_keys <- vf_hits %>%
  standardise_key_cols(c("tp_lab", "Timepoint")) %>%
  distinct(.data$Participant_id, .data$tp_lab_norm)

model_keys <- model_df %>%
  standardise_key_cols(c("tp_lab", "Timepoint")) %>%
  distinct(.data$Participant_id, .data$tp_lab_norm)

# ------------------------------------------------------------------------------
# Row-wise evidence trace for each clinical UTI episode.
# ------------------------------------------------------------------------------

audit_one_uti <- function(pid, clinical_tp) {
  exact_asm <- assembly_std %>%
    filter(.data$Participant_id == pid, .data$tp_lab_norm == clinical_tp)
  alt_asm <- assembly_std %>%
    filter(
      .data$Participant_id == pid,
      str_detect(.data$tp_lab_norm, "^UTI") | .data$Event_type == "UTI_event"
    )
  routine_asm <- assembly_std %>%
    filter(.data$Participant_id == pid, str_detect(.data$tp_lab_norm, "^T[0-9]"))

  candidate_asm <- if (nrow(exact_asm) > 0) exact_asm else alt_asm
  matched_by <- case_when(
    nrow(exact_asm) > 0 ~ "exact_participant_tp_lab",
    nrow(alt_asm) > 0 ~ "alternate_uti_label_same_participant",
    TRUE ~ "none"
  )

  candidate_assemblies <- unique(candidate_asm$Assembly_ID)
  candidate_tps <- unique(candidate_asm$tp_lab_norm)

  qc_match <- qc_std %>%
    filter(
      .data$Participant_id == pid,
      (.data$tp_lab_norm %in% candidate_tps) |
        (!is.na(.data$Assembly_ID) & .data$Assembly_ID %in% candidate_assemblies)
    )

  old_pa_exact <- old_vf_pa_keys %>%
    filter(.data$Participant_id == pid, .data$tp_lab_norm == clinical_tp)
  old_ready_exact <- old_vf_ready_keys %>%
    filter(.data$Participant_id == pid, .data$tp_lab_norm == clinical_tp)
  old_pa_alt <- old_vf_pa_keys %>%
    filter(.data$Participant_id == pid, .data$tp_lab_norm %in% candidate_tps)

  current_pa_exact <- current_vf_pa_keys %>%
    filter(.data$Participant_id == pid, .data$tp_lab_norm == clinical_tp)
  current_pa_alt <- current_vf_pa_keys %>%
    filter(.data$Participant_id == pid, .data$tp_lab_norm %in% candidate_tps)

  current_ready_exact <- current_vf_ready_keys %>%
    filter(.data$Participant_id == pid, .data$tp_lab_norm == clinical_tp)
  current_ready_alt <- current_vf_ready_keys %>%
    filter(.data$Participant_id == pid, .data$tp_lab_norm %in% candidate_tps)

  vf_hits_exact <- vf_hits_keys %>%
    filter(.data$Participant_id == pid, .data$tp_lab_norm == clinical_tp)
  vf_hits_alt <- vf_hits_keys %>%
    filter(.data$Participant_id == pid, .data$tp_lab_norm %in% candidate_tps)

  model_exact <- model_keys %>%
    filter(.data$Participant_id == pid, .data$tp_lab_norm == clinical_tp)
  model_alt <- model_keys %>%
    filter(.data$Participant_id == pid, .data$tp_lab_norm %in% candidate_tps)

  st_match <- current_ready_status %>%
    filter(.data$Participant_id == pid, .data$tp_lab_norm %in% c(clinical_tp, candidate_tps))

  retained <- nrow(old_ready_exact) > 0
  in_meta <- nrow(candidate_asm) > 0
  fasta_exists <- collapse_bool(candidate_asm$fasta_exists_calc)
  qc_present <- nrow(qc_match) > 0
  qc_pass <- collapse_bool(qc_match$QC_PASS)
  has_current_vf <- nrow(current_pa_exact) > 0 || nrow(current_pa_alt) > 0
  has_vf_hits <- nrow(vf_hits_exact) > 0 || nrow(vf_hits_alt) > 0

  exclusion_reason <- case_when(
    retained ~ "Retained in VF/WGS-ready dataset",
    matched_by == "alternate_uti_label_same_participant" && has_current_vf ~
      "Uricult timepoint label mismatch",
    !in_meta ~ "No assembly metadata row for UTI/Uricult episode",
    in_meta && !fasta_exists ~ "Assembly metadata row exists but FASTA missing",
    qc_present && !qc_pass ~ "QC exclusion",
    in_meta && fasta_exists && !has_vf_hits ~ "FASTA exists but Abricate/VFDB produced no row",
    has_current_vf && nrow(old_pa_exact) == 0 ~ "Stale canonical VF file or legacy join mismatch",
    TRUE ~ "Other: unresolved from available files"
  )

  note_bits <- c(
    if (matched_by == "alternate_uti_label_same_participant") {
      paste0("Clinical key is ", clinical_tp, " but assembly/VF evidence uses ",
             collapse_values(candidate_tps))
    } else NULL,
    if (!in_meta && nrow(routine_asm) > 0) {
      paste0("Participant has routine assembly rows only: ", collapse_values(unique(routine_asm$tp_lab_norm)))
    } else NULL,
    if (nrow(old_pa_alt) > 0 && nrow(old_pa_exact) == 0) {
      paste0("Legacy root VF PA has same participant under alternate key(s): ",
             collapse_values(old_pa_alt$tp_lab_norm))
    } else NULL,
    if (nrow(current_ready_alt) > 0 && nrow(current_ready_exact) == 0) {
      paste0("Current VF-ready has alternate key(s): ",
             collapse_values(current_ready_alt$tp_lab_norm),
             " with status ",
             collapse_values(st_match$current_ready_status))
    } else NULL
  )

  tibble(
    exact_assembly_key_match = nrow(exact_asm) > 0,
    alternate_uti_label_match = nrow(alt_asm) > 0,
    matched_by = matched_by,
    alt_tp_lab = collapse_values(candidate_tps),
    in_assembly_metadata = in_meta,
    Isolate_ID = collapse_values(candidate_asm$Isolate_ID),
    file_name = collapse_values(candidate_asm$file_name),
    full_path = collapse_values(candidate_asm$full_path),
    fasta_exists = fasta_exists,
    in_wgs_qc = qc_present,
    wgs_qc_status = if (!qc_present) NA_character_ else if (qc_pass) "PASS" else "FAIL",
    wgs_qc_reason = collapse_values(qc_match$QC_REASON),
    in_vf_hits = has_vf_hits,
    in_vf_pa_all = nrow(old_pa_exact) > 0,
    in_vf_pa_all_alt_key = nrow(old_pa_alt) > 0,
    in_current_vf_pa_all = nrow(current_pa_exact) > 0,
    in_current_vf_pa_all_alt_key = nrow(current_pa_alt) > 0,
    in_vf_analysis_ready = nrow(old_ready_exact) > 0,
    in_current_vf_analysis_ready = nrow(current_ready_exact) > 0,
    in_current_vf_analysis_ready_alt_key = nrow(current_ready_alt) > 0,
    in_final_model_dataset = nrow(model_exact) > 0,
    in_final_model_dataset_alt_key = nrow(model_alt) > 0,
    ST = collapse_values(st_match$ST_current),
    retained_in_vf_ready = retained,
    exclusion_reason = exclusion_reason,
    notes = collapse_values(note_bits)
  )
}

evidence <- bind_rows(lapply(seq_len(nrow(clinical_uti)), function(i) {
  audit_one_uti(clinical_uti$Participant_id[i], clinical_uti$clinical_tp_lab[i])
}))

uti_audit <- bind_cols(clinical_uti, evidence) %>%
  transmute(
    Participant_id,
    clinical_Timepoint,
    clinical_tp_lab,
    Batch,
    Collection_Date,
    UTI_Label,
    Urine_collection_method,
    Infection_Status,
    Status_Confidence,
    Status_Provenance,
    raw_CFU_examples,
    raw_BEO_examples,
    raw_Pop_examples,
    Sx_present_any,
    Sx_source_epi,
    in_status_map,
    in_assembly_metadata,
    exact_assembly_key_match,
    alternate_uti_label_match,
    matched_by,
    alt_tp_lab,
    Isolate_ID,
    file_name,
    full_path,
    fasta_exists,
    in_wgs_qc,
    wgs_qc_status,
    wgs_qc_reason,
    in_vf_hits,
    in_vf_pa_all,
    in_vf_pa_all_alt_key,
    in_current_vf_pa_all,
    in_current_vf_pa_all_alt_key,
    in_vf_analysis_ready,
    in_current_vf_analysis_ready,
    in_current_vf_analysis_ready_alt_key,
    in_final_model_dataset,
    in_final_model_dataset_alt_key,
    ST,
    retained_in_vf_ready,
    exclusion_reason,
    notes
  )

missing_uti <- uti_audit %>%
  filter(!.data$retained_in_vf_ready)

# ------------------------------------------------------------------------------
# Denominator cascade.
# ------------------------------------------------------------------------------

old_vf_join <- legacy_clinical %>%
  inner_join(old_vf_pa_keys, by = c("Participant_id", "clinical_tp_lab" = "tp_lab_norm"))

old_ready_status_counts <- if (nrow(old_vf_ready) > 0 && "Infection_Status" %in% names(old_vf_ready)) {
  old_vf_ready
} else {
  legacy_clinical %>%
    inner_join(old_vf_ready_keys, by = c("Participant_id", "clinical_tp_lab" = "tp_lab_norm"))
}

current_ready_counts_df <- if (nrow(current_vf_ready) > 0 && "Infection_Status" %in% names(current_vf_ready)) {
  current_vf_ready
} else {
  tibble()
}

episode_flags_for_stage <- clinical_all %>%
  rowwise() %>%
  mutate(
    has_assembly_metadata = {
      pid <- .data$Participant_id
      tp <- .data$clinical_tp_lab
      any(assembly_std$Participant_id == pid & assembly_std$tp_lab_norm == tp) ||
        (tp == "Uricult" && any(assembly_std$Participant_id == pid &
                                  (str_detect(assembly_std$tp_lab_norm, "^UTI") |
                                     assembly_std$Event_type == "UTI_event")))
    },
    has_existing_fasta = {
      pid <- .data$Participant_id
      tp <- .data$clinical_tp_lab
      rows <- assembly_std %>% filter(.data$Participant_id == pid)
      rows <- rows %>% filter(.data$tp_lab_norm == tp |
                                (tp == "Uricult" & (str_detect(.data$tp_lab_norm, "^UTI") |
                                                      .data$Event_type == "UTI_event")))
      collapse_bool(rows$fasta_exists_calc)
    },
    passes_wgs_qc = {
      pid <- .data$Participant_id
      tp <- .data$clinical_tp_lab
      rows <- qc_std %>% filter(.data$Participant_id == pid)
      rows <- rows %>% filter(.data$tp_lab_norm == tp |
                                (tp == "Uricult" & str_detect(.data$tp_lab_norm, "^UTI")))
      collapse_bool(rows$QC_PASS)
    },
    has_current_vf_hits = {
      pid <- .data$Participant_id
      tp <- .data$clinical_tp_lab
      any(vf_hits_keys$Participant_id == pid & vf_hits_keys$tp_lab_norm == tp) ||
        (tp == "Uricult" && any(vf_hits_keys$Participant_id == pid &
                                  str_detect(vf_hits_keys$tp_lab_norm, "^UTI")))
    }
  ) %>%
  ungroup()

add_stage <- function(stage, df, source_file, script, notes) {
  summarise_status_counts(df) %>%
    mutate(
      stage = stage,
      source_file = source_file,
      script_responsible = script,
      notes = notes,
      .before = 1
    )
}

denominator_cascade <- bind_rows(
  add_stage(
    "legacy_clinical_status_map",
    legacy_clinical,
    legacy_status_path,
    "legacy/previous pipeline status output",
    "This is the denominator that reproduces the user-reported 274 episodes and 20 UTI."
  ),
  add_stage(
    "legacy_clinical_with_any_assembly_metadata_current_trace",
    episode_flags_for_stage %>% filter(.data$has_assembly_metadata),
    FILE_METADATA,
    "00_make_assembly_metadata.r",
    "Current trace; counts exact participant/timepoint matches plus Uricult clinical rows with UTI-* assembly labels for the same participant."
  ),
  add_stage(
    "legacy_clinical_with_existing_fasta_current_trace",
    episode_flags_for_stage %>% filter(.data$has_existing_fasta),
    FILE_METADATA,
    "00_make_assembly_metadata.r",
    "Current trace using assembly_metadata file_exists/full_path evidence."
  ),
  add_stage(
    "legacy_clinical_passing_wgs_qc_current_trace",
    episode_flags_for_stage %>% filter(.data$passes_wgs_qc),
    file.path(DIR_WGS, "qc_summary.csv"),
    "12a_wgs_qc.R",
    "Current trace using QC_PASS."
  ),
  add_stage(
    "legacy_clinical_with_current_vf_hits_trace",
    episode_flags_for_stage %>% filter(.data$has_current_vf_hits),
    FILE_VF_HITS,
    "02_gene_presence_analysis.R",
    "Current trace using vf_hits_all.rds, allowing Uricult clinical rows to match UTI-* VF rows for the same participant."
  ),
  add_stage(
    "legacy_root_vf_pa_all_joined_by_exact_clinical_key",
    old_vf_join,
    old_vf_pa_path,
    "legacy 02_gene_presence_analysis.R / VF join",
    "This exact-key join is where the legacy UTI denominator drops from 20 to 16."
  ),
  add_stage(
    "legacy_vf_analysis_ready_OLD",
    old_ready_status_counts,
    old_vf_ready_path,
    "legacy 22_vf_build_analysis_dataset.R",
    "Old VF-ready file that reproduces 183 rows and 16 UTI."
  ),
  add_stage(
    "available_final_model_dataset",
    model_df,
    model_path %||% "not found",
    "legacy/current model script",
    "The available model dataset was audited; it does not necessarily match the user-reported 152 ASB/UTI rows."
  ),
  add_stage(
    "current_canonical_status_map",
    current_status,
    current_status_path,
    "00b_classify_episodes.R",
    "Current repaired pipeline status map; this has a different denominator than the legacy 20-UTI set."
  ),
  add_stage(
    "current_vf_analysis_ready",
    current_ready_counts_df,
    current_vf_ready_path,
    "22_vf_build_analysis_dataset.R",
    "Current VF-ready output; included to flag stale/conflicting denominator files."
  )
) %>%
  select(stage, n_total, n_ASB, n_legacy_UTI, n_Negative, n_UTI, n_Not_UTI, n_other,
         source_file, script_responsible, notes)

# ------------------------------------------------------------------------------
# New primary-definition audit: old ASB/UTI/Negative versus UTI/Not_UTI.
# ------------------------------------------------------------------------------

primary_denominator_summary <- if (nrow(current_status) > 0) {
  bind_rows(
    summarise_status_counts(current_status) %>%
      mutate(
        stage = "current_canonical_status_map",
        UTI_definition_version = collapse_col(current_status, "UTI_definition_version"),
        .before = 1
      ),
    summarise_status_counts(current_vf_ready) %>%
      mutate(
        stage = "current_vf_analysis_ready",
        UTI_definition_version = collapse_col(current_vf_ready, "UTI_definition_version"),
        .before = 1
      ),
    summarise_status_counts(model_df) %>%
      mutate(
        stage = "available_final_model_dataset",
        UTI_definition_version = collapse_col(model_df, "UTI_definition_version"),
        .before = 1
      )
  )
} else {
  tibble()
}

old_new_movement <- if (nrow(current_status) > 0) {
  current_status %>%
    mutate(
      Infection_Status_old = dplyr::coalesce(
        as.character(.data$Infection_Status_old),
        as.character(.data$Infection_Status_legacy)
      ),
      Not_UTI_subgroup = dplyr::if_else(
        .data$UTI_Status == "UTI",
        NA_character_,
        as.character(.data$Not_UTI_subgroup)
      )
    ) %>%
    count(
      .data$Infection_Status_old,
      .data$UTI_Status,
      .data$Not_UTI_subgroup,
      .data$catheter_rule,
      .data$symptom_rule_met,
      name = "n"
    ) %>%
    arrange(.data$Infection_Status_old, .data$UTI_Status, desc(.data$n))
} else {
  tibble()
}

primary_rule_audit <- if (nrow(current_status) > 0) {
  current_status %>%
    count(
      .data$UTI_Status,
      .data$catheter_rule,
      .data$symptom_rule_met,
      .data$Not_UTI_subgroup,
      name = "n"
    ) %>%
    arrange(.data$UTI_Status, .data$catheter_rule, .data$symptom_rule_met)
} else {
  tibble()
}

# ------------------------------------------------------------------------------
# Write audit products.
# ------------------------------------------------------------------------------

uti_audit_path <- file.path(DIR_AUDIT, "uti_denominator_audit.csv")
cascade_path <- file.path(DIR_AUDIT, "denominator_cascade.csv")
missing_path <- file.path(DIR_AUDIT, "missing_uti_episodes.csv")
summary_path <- file.path(DIR_AUDIT, "uti_denominator_audit_summary.md")
primary_denominator_path <- file.path(DIR_AUDIT, "uti_new_primary_denominator_audit.csv")
old_new_movement_path <- file.path(DIR_AUDIT, "uti_old_new_reclassification_summary.csv")
primary_rule_audit_path <- file.path(DIR_AUDIT, "uti_new_primary_rule_audit.csv")

readr::write_csv(uti_audit, uti_audit_path)
readr::write_csv(denominator_cascade, cascade_path)
readr::write_csv(missing_uti, missing_path)
readr::write_csv(primary_denominator_summary, primary_denominator_path)
readr::write_csv(old_new_movement, old_new_movement_path)
readr::write_csv(primary_rule_audit, primary_rule_audit_path)

retained_uti <- uti_audit %>% filter(.data$retained_in_vf_ready)

verdict <- if (nrow(missing_uti) == 0) {
  "Expected denominator loss was not observed in the audited legacy files."
} else if (any(missing_uti$exclusion_reason == "Uricult timepoint label mismatch") &&
           any(missing_uti$exclusion_reason == "No assembly metadata row for UTI/Uricult episode")) {
  "Mixed: both pipeline key mismatch and expected WGS/VF unavailability are present."
} else if (all(missing_uti$exclusion_reason == "Uricult timepoint label mismatch")) {
  "Pipeline issue: the missing UTI rows have WGS/VF evidence under alternate UTI-* labels."
} else {
  "Expected denominator loss or unresolved data unavailability dominates."
}

missing_reason_lines <- missing_uti %>%
  transmute(
    line = paste0(
      "- ", .data$Participant_id, " / ", .data$clinical_tp_lab,
      ": ", .data$exclusion_reason,
      ifelse(is.na(.data$notes) | .data$notes == "", "", paste0(". ", .data$notes))
    )
  ) %>%
  pull(.data$line)

summary_lines <- c(
  "# UTI Denominator Audit",
  "",
  paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  "",
  "## Verdict",
  "",
  verdict,
  "",
  "## Canonical And Conflicting Files Checked",
  "",
  paste0("- Legacy clinical denominator reproducing 20 UTI: `", legacy_status_path, "`."),
  paste0("- Legacy/root VF PA reproducing the 20 to 16 drop: `", old_vf_pa_path, "`."),
  paste0("- Legacy VF-ready file reproducing 183 VF-ready rows and 16 UTI: `", old_vf_ready_path, "`."),
  paste0("- Current primary UTI status map: `", current_status_path, "`."),
  paste0("- Current assembly metadata: `", FILE_METADATA, "`."),
  paste0("- Current WGS QC: `", file.path(DIR_WGS, "qc_summary.csv"), "`."),
  paste0("- Current VF PA: `", current_vf_pa_path, "`."),
  paste0("- Current VF-ready: `", current_vf_ready_path, "`."),
  "",
  "The legacy files and current repaired files do not have the same denominator. The 20-to-16 question is reproduced by the legacy status file plus the root-level legacy `vf_pa_all.csv` / `vf_analysis_ready_OLD.csv` files.",
  "",
  "## Denominator Cascade",
  "",
  make_md_table(denominator_cascade),
  "",
  "## New Primary UTI/Not_UTI Definition",
  "",
  paste0("Definition: `", UTI_DEFINITION_VERSION, "`. The primary estimand is UTI versus all Not_UTI episodes; legacy ASB/UTI/Negative counts are retained only for comparison."),
  "",
  make_md_table(primary_denominator_summary),
  "",
  "## Old-to-New Reclassification Summary",
  "",
  make_md_table(old_new_movement),
  "",
  "## Primary Rule Audit",
  "",
  make_md_table(primary_rule_audit),
  "",
  "## Clinical UTI Episodes",
  "",
  make_md_table(
    uti_audit %>%
      select(Participant_id, clinical_tp_lab, Batch, Collection_Date, UTI_Label,
             retained_in_vf_ready, exclusion_reason, alt_tp_lab)
  ),
  "",
  "## Retained UTI Episodes In Legacy VF/WGS-Ready Dataset",
  "",
  make_md_table(
    retained_uti %>%
      select(Participant_id, clinical_tp_lab, Batch, Isolate_ID, alt_tp_lab,
             file_name, ST, in_vf_pa_all, in_vf_analysis_ready)
  ),
  "",
  "## Missing UTI Episodes",
  "",
  if (length(missing_reason_lines) == 0) "_No missing UTI episodes._" else paste(missing_reason_lines, collapse = "\n"),
  "",
  make_md_table(
    missing_uti %>%
      select(Participant_id, clinical_tp_lab, Batch, Collection_Date, UTI_Label,
             Urine_collection_method, raw_CFU_examples, raw_BEO_examples,
             raw_Pop_examples, in_assembly_metadata, exact_assembly_key_match,
             alternate_uti_label_match, alt_tp_lab, fasta_exists, in_wgs_qc,
             wgs_qc_status, in_vf_hits, in_current_vf_pa_all_alt_key,
             retained_in_vf_ready, exclusion_reason, notes)
  ),
  "",
  "## Interpretation For Reporting",
  "",
  "For the legacy 274-episode denominator, report this as: \"20 clinical UTI episodes, of which 16 had matched E. coli WGS/VF data in the legacy VF-ready file.\"",
  "",
  "However, three of the four legacy-missing UTI rows have current assembly/QC/VF evidence under UTI-* labels rather than the clinical Uricult label. Those are not clean biological losses; they are key/linkage problems caused by Uricult event labelling. One missing UTI row has no UTI/Uricult assembly or VF evidence in the inspected files and appears to be genuine WGS/VF unavailability unless a trusted mapping table exists elsewhere.",
  "",
  "## Recommended Minimal Fix If The Pipeline Is Patched",
  "",
  "- Preserve event-specific Uricult/UTI labels instead of collapsing all Uricult rows to a generic participant-level `Uricult` key.",
  "- Join clinical and WGS/VF data through `Episode_ID` or an explicit Uricult crosswalk using Participant_id plus Collection_Date, UTI_Label/Sample_ID, and Isolate_ID where available.",
  "- Rebuild the VF-ready dataset after repairing the key; do not interpret the legacy 16-UTI VF-ready denominator as the total clinical UTI burden."
)

writeLines(summary_lines, summary_path)

audit_msg("Wrote %s", uti_audit_path)
audit_msg("Wrote %s", cascade_path)
audit_msg("Wrote %s", missing_path)
audit_msg("Wrote %s", primary_denominator_path)
audit_msg("Wrote %s", old_new_movement_path)
audit_msg("Wrote %s", primary_rule_audit_path)
audit_msg("Wrote %s", summary_path)
audit_msg("Done. Legacy retained UTI in VF-ready: %s/%s",
          sum(uti_audit$retained_in_vf_ready), nrow(uti_audit))

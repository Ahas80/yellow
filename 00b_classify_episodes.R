#!/usr/bin/env Rscript
# ==============================================================================
# 00b_classify_episodes.R
# ==============================================================================
#
# GOAL:
#   Build the authoritative primary UTI status map for the pipeline.
#
# PRIMARY STATUS DEFINITION:
#   UTI_Status uses catheter-aware signs/symptoms plus culture support at the
#   lower >=10^3 CFU/mL threshold:
#
#     UTI     = culture_supports_uti AND symptom_compatible_uti
#     Not_UTI = all other episodes, with Not_UTI_subgroup retained
#
# LEGACY STATUS:
#   Infection_Status, Infection_Status_legacy, and Infection_Status_old preserve
#   the former ASB / UTI / Negative framing for comparability only.
# ==============================================================================

source("00_config.R")
source(here::here("R", "clinical_helpers.R"))

suppressPackageStartupMessages({
    library(dplyr)
    library(readr)
    library(stringr)
    library(tidyr)
})

msg("Starting 00b_classify_episodes.R")
msg("Primary UTI definition: %s", UTI_DEFINITION_VERSION)

input_rds <- file.path(DIR_CLINICAL, "intermediate", "clinical_merged.rds")
if (!file.exists(input_rds)) stop("Intermediate data not found. Run 00a_load_clean_clinical.R first.")
data_list <- readRDS(input_rds)

metadata <- data_list$metadata
sym_all_raw <- data_list$sym_all_raw
sns_12 <- data_list$sns_12

symptom_fields <- c(
    "dysuria_present",
    "urgency_present",
    "frequency_present",
    "incontinence_present",
    "pus_present",
    "flankpain_present",
    "suprapubic_pain_present",
    "fever_present",
    "rigors_present",
    "delirium_present",
    "other_sxs_present"
)

empty_symptom_episode_tbl <- function() {
    out <- tibble(Participant_id = character(), Timepoint = character())
    for (nm in symptom_fields) out[[nm]] <- logical()
    out
}

collapse_bool_episode <- function(x) {
    x <- x[!is.na(x)]
    if (!length(x)) return(NA)
    any(x)
}

collapse_chr_unique <- function(x, max_items = 8) {
    vals <- sort(unique(na.omit(as.character(x))))
    vals <- vals[nzchar(vals)]
    if (!length(vals)) return(NA_character_)
    paste(head(vals, max_items), collapse = ";")
}

collapse_catheter_rule <- function(x) {
    vals <- unique(na.omit(as.character(x)))
    vals <- vals[nzchar(vals)]
    known <- setdiff(vals, "Unknown_collection_method")
    if (!length(known)) return("Unknown_collection_method")
    if (length(unique(known)) == 1) return(unique(known))
    "Unknown_collection_method"
}

collapse_collection_norm <- function(x) {
    vals <- unique(na.omit(as.character(x)))
    vals <- vals[nzchar(vals)]
    known <- setdiff(vals, "unknown")
    if (!length(known)) return("unknown")
    if (length(unique(known)) == 1) return(unique(known))
    "mixed_ambiguous"
}

collapse_cfu_source <- function(source, culture_supports_uti) {
    source <- unique(na.omit(as.character(source)))
    if (!length(source)) return(NA_character_)
    if ("cfu_ge_1e3_lower_bound" %in% source) return("cfu_ge_1e3_lower_bound")
    beoord_sources <- c(
        "beoord_plus3_fallback_1e5",
        "beoord_plus2_fallback_1e4",
        "beoord_plus1_fallback_1e3"
    )
    matched_beoord <- intersect(beoord_sources, source)
    if (length(matched_beoord)) return(matched_beoord[1])
    if (any(culture_supports_uti %in% FALSE, na.rm = TRUE)) return(paste(source, collapse = ";"))
    paste(source, collapse = ";")
}

# ------------------------------------------------------------------------------
# 1. Parse individual S&S columns
# ------------------------------------------------------------------------------
if (is.null(sym_all_raw) || nrow(sym_all_raw) == 0) {
    sym_all_1row <- empty_symptom_episode_tbl()
    sym_all_1row$Sx_columns_detected <- logical()
} else {
    sym_std <- sym_all_raw %>%
        mutate(
            Participant_id = str_trim(as.character(Participant_id)),
            Timepoint = canon_tp(Timepoint)
        )
    sym_flags <- derive_symptom_flags(sym_std)
    sym_all_1row <- bind_cols(sym_std %>% select(Participant_id, Timepoint), sym_flags) %>%
        group_by(Participant_id, Timepoint) %>%
        summarise(across(all_of(symptom_fields), collapse_bool_episode), .groups = "drop") %>%
        distinct() %>%
        mutate(Sx_columns_detected = TRUE)
}

msg("Parsed S&S fields for %d participant-timepoint rows", nrow(sym_all_1row))

# ------------------------------------------------------------------------------
# 2. Row-level legacy and primary clinical classification inputs
# ------------------------------------------------------------------------------
meta_plus <- metadata %>%
    left_join(sns_12,
        by = c("Participant_id", "Timepoint", "Batch"),
        relationship = "many-to-one"
    ) %>%
    left_join(sym_all_1row,
        by = c("Participant_id", "Timepoint"),
        relationship = "many-to-one"
    )

for (nm in symptom_fields) {
    if (!nm %in% names(meta_plus)) meta_plus[[nm]] <- NA
}
if (!"SnS_status" %in% names(meta_plus)) meta_plus$SnS_status <- NA_integer_
if (!"Collection_Date" %in% names(meta_plus)) meta_plus$Collection_Date <- NA_character_
if (!"UTI_Label" %in% names(meta_plus)) meta_plus$UTI_Label <- NA_character_
if (!"Urine collection method" %in% names(meta_plus)) meta_plus[["Urine collection method"]] <- NA_character_
if (!"CFU_Count" %in% names(meta_plus)) meta_plus$CFU_Count <- NA_character_
if (!"Beoord" %in% names(meta_plus)) meta_plus$Beoord <- NA_character_
if (!"Population" %in% names(meta_plus)) meta_plus$Population <- NA_character_

cfu_detail <- parse_cfu_detail(meta_plus$CFU_Count)
symptom_rule_row <- derive_symptom_rule(
    catheter_rule = classify_catheter_rule(meta_plus[["Urine collection method"]]),
    dysuria_present = meta_plus$dysuria_present,
    urgency_present = meta_plus$urgency_present,
    frequency_present = meta_plus$frequency_present,
    incontinence_present = meta_plus$incontinence_present,
    pus_present = meta_plus$pus_present,
    flankpain_present = meta_plus$flankpain_present,
    fever_present = meta_plus$fever_present,
    rigors_present = meta_plus$rigors_present,
    delirium_present = meta_plus$delirium_present
)

classified <- meta_plus %>%
    bind_cols(cfu_detail %>% select(-cfu_raw)) %>%
    mutate(
        cfu_raw = as.character(CFU_Count),
        cfu_recorded = !is_unknown_clinical_text(CFU_Count),
        beoord_cat = parse_beoord_cat(Beoord),
        beoord_cfu_lower_bound = beoord_cfu_lower_bound(beoord_cat),
        cfu_cat = cfu_bucket(CFU_Count, threshold = UTI_CFU_THRESHOLD_LEGACY),
        culture_pos = case_when(
            cfu_recorded ~ (cfu_cat == ">=1e5"),
            !cfu_recorded & !is.na(beoord_cat) ~ (beoord_cat == "+++"),
            TRUE ~ NA
        ),
        sx_present_pop = population_to_sns(Population),
        Sx_present_sym = logic_or_known(
            dysuria_present, urgency_present, frequency_present, incontinence_present,
            pus_present, flankpain_present, suprapubic_pain_present, fever_present,
            rigors_present, delirium_present, other_sxs_present
        ),
        Sx_present_sym_legacy = ifelse(is.na(Sx_present_sym) & Sx_columns_detected %in% TRUE, FALSE, Sx_present_sym),
        Sx_present_final = case_when(
            !is.na(sx_present_pop) ~ sx_present_pop,
            !is.na(Sx_present_sym_legacy) ~ Sx_present_sym_legacy,
            Batch %in% c(1L, 2L) & is.na(Sx_present_sym_legacy) & !is.na(SnS_status) ~ case_when(
                SnS_status == 2L ~ TRUE,
                SnS_status == 0L ~ FALSE,
                TRUE ~ as.logical(NA)
            ),
            TRUE ~ as.logical(NA)
        ),
        Infection_Status = case_when(
            culture_pos == TRUE & Sx_present_final %in% TRUE ~ "UTI",
            culture_pos == TRUE & Sx_present_final %in% FALSE ~ "ASB",
            culture_pos == TRUE & is.na(Sx_present_final) ~ "Culture-positive, S&S unknown",
            culture_pos == FALSE ~ "Negative",
            TRUE ~ "None"
        ),
        Infection_Status_legacy = Infection_Status,
        Infection_Status_old = Infection_Status,
        Status_Confidence = case_when(
            !is.na(sx_present_pop) & cfu_recorded ~ "High",
            is.na(sx_present_pop) & !is.na(Sx_present_final) & cfu_recorded ~ "Medium",
            !cfu_recorded | is.na(Sx_present_final) ~ "Low",
            TRUE ~ "Unknown"
        ),
        Status_Provenance = case_when(
            !is.na(sx_present_pop) ~ "Population_field",
            !is.na(Sx_present_sym) ~ "SxS_columns",
            Batch %in% c(1L, 2L) & !is.na(SnS_status) ~ "SnS_status_fallback",
            TRUE ~ "Unknown"
        ),
        urine_collection_method_raw = as.character(`Urine collection method`),
        urine_collection_method_norm = normalise_urine_collection_method(urine_collection_method_raw),
        catheter_rule = classify_catheter_rule(urine_collection_method_raw),
        cfu_threshold_used_for_uti = "1e3",
        cfu_threshold_source = case_when(
            cfu_recorded & cfu_ge_1e3 %in% TRUE ~ "cfu_ge_1e3_lower_bound",
            cfu_recorded & cfu_ge_1e3 %in% FALSE ~ "cfu_below_1e3_lower_bound",
            cfu_recorded & is.na(cfu_ge_1e3) ~ "cfu_ambiguous",
            !cfu_recorded & !is.na(beoord_cfu_lower_bound) & beoord_cfu_lower_bound >= UTI_CFU_THRESHOLD_PRIMARY ~
                beoord_primary_source(beoord_cat),
            !cfu_recorded & !is.na(beoord_cfu_lower_bound) ~ "beoord_below_primary",
            TRUE ~ "missing_cfu_and_beoord"
        ),
        culture_supports_uti = case_when(
            cfu_recorded & !is.na(cfu_ge_1e3) ~ cfu_ge_1e3,
            cfu_recorded & is.na(cfu_ge_1e3) ~ as.logical(NA),
            !cfu_recorded & !is.na(beoord_cfu_lower_bound) ~ beoord_cfu_lower_bound >= UTI_CFU_THRESHOLD_PRIMARY,
            TRUE ~ as.logical(NA)
        )
    ) %>%
    bind_cols(symptom_rule_row)

classified <- bind_cols(
    classified,
    derive_uti_status(
        culture_supports_uti = classified$culture_supports_uti,
        symptom_compatible_uti = classified$symptom_compatible_uti,
        catheter_rule = classified$catheter_rule,
        cfu_threshold_source = classified$cfu_threshold_source,
        symptom_rule_met = classified$symptom_rule_met
    )
)

# ------------------------------------------------------------------------------
# 3. Episode-level collapse
# ------------------------------------------------------------------------------
episode_core <- classified %>%
    group_by(Participant_id, Timepoint) %>%
    summarise(
        Batch = paste(sort(unique(Batch)), collapse = ";"),
        Collection_Date = collapse_chr_unique(Collection_Date),
        UTI_Label = collapse_chr_unique(UTI_Label),
        Urine_collection_method = collapse_chr_unique(`Urine collection method`),
        urine_collection_method_raw = collapse_chr_unique(urine_collection_method_raw),
        urine_collection_method_norm = collapse_collection_norm(urine_collection_method_norm),
        catheter_rule = collapse_catheter_rule(catheter_rule),
        culture_pos_epi = collapse_bool_episode(culture_pos),
        cfu_recorded_any = any(ifelse(is.na(cfu_recorded), FALSE, cfu_recorded)),
        cfu_ge_1e3 = collapse_bool_episode(cfu_ge_1e3),
        cfu_ge_1e4 = collapse_bool_episode(cfu_ge_1e4),
        cfu_ge_1e5 = collapse_bool_episode(cfu_ge_1e5),
        cfu_ge_1e5_any = cfu_ge_1e5 %in% TRUE,
        beoord_plus3_any = any(beoord_cat == "+++", na.rm = TRUE),
        beoord_cfu_lower_bound = if (all(is.na(beoord_cfu_lower_bound))) NA_real_ else max(beoord_cfu_lower_bound, na.rm = TRUE),
        culture_supports_uti = collapse_bool_episode(culture_supports_uti),
        cfu_threshold_used_for_uti = "1e3",
        cfu_threshold_source = collapse_cfu_source(cfu_threshold_source, culture_supports_uti),
        cfu_raw = collapse_chr_unique(cfu_raw),
        cfu_raw_parsed = collapse_chr_unique(cfu_raw_parsed),
        cfu_lower_bound = if (all(is.na(cfu_lower_bound))) NA_real_ else min(cfu_lower_bound, na.rm = TRUE),
        cfu_upper_bound = if (all(is.na(cfu_upper_bound))) NA_real_ else max(cfu_upper_bound, na.rm = TRUE),
        beoord_cat = collapse_chr_unique(beoord_cat),
        Sx_present_any = collapse_bool_episode(Sx_present_final),
        Sx_source_epi = case_when(
            any(!is.na(sx_present_pop)) ~ "Population",
            any(!is.na(Sx_present_sym)) ~ "S&S columns",
            any(!is.na(SnS_status)) ~ "SnS_status fallback",
            TRUE ~ NA_character_
        ),
        across(all_of(symptom_fields), collapse_bool_episode),
        raw_CFU_examples = collapse_chr_unique(CFU_Count, max_items = 5),
        raw_BEO_examples = collapse_chr_unique(Beoord, max_items = 5),
        raw_Pop_examples = collapse_chr_unique(Population, max_items = 5),
        .groups = "drop"
    ) %>%
    mutate(
        Collection_Date = na_if(Collection_Date, ""),
        UTI_Label = na_if(UTI_Label, ""),
        Urine_collection_method = na_if(Urine_collection_method, ""),
        tp_lab = normalise_timepoint_preserve_events(Timepoint),
        Event_type = episode_event_type(tp_lab),
        Infection_Status = case_when(
            culture_pos_epi %in% TRUE & Sx_present_any %in% TRUE ~ "UTI",
            culture_pos_epi %in% TRUE & Sx_present_any %in% FALSE ~ "ASB",
            culture_pos_epi %in% TRUE & is.na(Sx_present_any) ~ "Culture-positive, S&S unknown",
            culture_pos_epi %in% FALSE ~ "Negative",
            TRUE ~ "None"
        ),
        Infection_Status_legacy = Infection_Status,
        Infection_Status_old = Infection_Status,
        Status_Confidence_epi = case_when(
            Sx_source_epi == "Population" & cfu_recorded_any ~ "High",
            Sx_source_epi %in% c("S&S columns", "SnS_status fallback") & cfu_recorded_any ~ "Medium",
            TRUE ~ "Low"
        ),
        UTI_definition_version = UTI_DEFINITION_VERSION
    )

episode_rule <- derive_symptom_rule(
    catheter_rule = episode_core$catheter_rule,
    dysuria_present = episode_core$dysuria_present,
    urgency_present = episode_core$urgency_present,
    frequency_present = episode_core$frequency_present,
    incontinence_present = episode_core$incontinence_present,
    pus_present = episode_core$pus_present,
    flankpain_present = episode_core$flankpain_present,
    fever_present = episode_core$fever_present,
    rigors_present = episode_core$rigors_present,
    delirium_present = episode_core$delirium_present
)

episode_tbl <- bind_cols(
    episode_core,
    episode_rule,
    derive_uti_status(
        culture_supports_uti = episode_core$culture_supports_uti,
        symptom_compatible_uti = episode_rule$symptom_compatible_uti,
        catheter_rule = episode_core$catheter_rule,
        cfu_threshold_source = episode_core$cfu_threshold_source,
        symptom_rule_met = episode_rule$symptom_rule_met
    )
) %>%
    mutate(
        Episode_ID = build_episode_id(., timepoint_col = "tp_lab", event_col = "Event_type",
                                      date_col = "Collection_Date"),
        changed_to_new_uti = UTI_Status == "UTI" & Infection_Status_old != "UTI",
        changed_due_to_lower_cfu_threshold = UTI_Status == "UTI" & !(cfu_ge_1e5 %in% TRUE) & cfu_ge_1e3 %in% TRUE,
        changed_due_to_catheter_aware_sns = UTI_Status == "UTI" & catheter_rule == "B_indwelling" &
            symptom_rule_met == "B_systemic_catheter",
        movement_category = paste0(Infection_Status_old, "_to_", UTI_Status,
                                   ifelse(is.na(Not_UTI_subgroup), "", paste0("_", Not_UTI_subgroup)))
    )

episode_tbl <- apply_manual_sample_curation(
    episode_tbl,
    context = "clinical_status_map",
    write_audit = TRUE
)

# ------------------------------------------------------------------------------
# 4. Validation and audit outputs
# ------------------------------------------------------------------------------
assert_unique_keys(
    episode_tbl,
    c("Episode_ID"),
    context = "status_map Episode_ID",
    out_path = file.path(DIR_QC, "status_map_duplicate_episode_ids.csv")
)

assert_unique_keys(
    episode_tbl,
    c("Participant_id", "tp_lab"),
    context = "status_map Participant_id + tp_lab",
    out_path = file.path(DIR_QC, "status_map_duplicate_episode_keys.csv")
)

if (any(!episode_tbl$UTI_Status %in% c("UTI", "Not_UTI"))) {
    stop("UTI_Status contains values outside UTI/Not_UTI")
}
if (any(!episode_tbl$UTI_binary %in% c(0L, 1L))) {
    stop("UTI_binary contains values outside 0/1")
}
if (any(episode_tbl$UTI_Status == "UTI" & !(episode_tbl$culture_supports_uti %in% TRUE), na.rm = TRUE)) {
    stop("Validation failed: at least one UTI row lacks culture_supports_uti == TRUE")
}
if (any(episode_tbl$UTI_Status == "UTI" & !(episode_tbl$symptom_compatible_uti %in% TRUE), na.rm = TRUE)) {
    stop("Validation failed: at least one UTI row lacks symptom_compatible_uti == TRUE")
}
if (any(episode_tbl$catheter_rule == "B_indwelling" &
        episode_tbl$symptom_rule_met %in% c("A_local_urinary_symptom", "A_flankpain_plus_systemic"), na.rm = TRUE)) {
    stop("Validation failed: catheter rows used rule A symptom logic")
}
if (any(episode_tbl$catheter_rule == "A_non_indwelling" &
        episode_tbl$symptom_rule_met == "B_systemic_catheter", na.rm = TRUE)) {
    stop("Validation failed: non-catheter rows used rule B symptom logic")
}

write_csv(episode_tbl, FILE_STATUS_MAP)
msg("Saved status_map.csv (%d rows) to %s", nrow(episode_tbl), FILE_STATUS_MAP)

legacy_comparison <- episode_tbl %>%
    select(
        Participant_id, Timepoint, tp_lab, Episode_ID, Batch,
        Infection_Status_old, Infection_Status_legacy, Infection_Status,
        UTI_Status, UTI_binary, Not_UTI_subgroup, movement_category,
        changed_to_new_uti, changed_due_to_lower_cfu_threshold,
        changed_due_to_catheter_aware_sns, UTI_definition_version,
        analysis_include_primary, analysis_exclusion_reason,
        duplicate_role, duplicate_of_participant_id, duplicate_of_tp_lab,
        allow_secondary_duplicate_qc, duplicate_use_note,
        manual_curation_applied, manual_curation_note,
        catheter_rule, symptom_rule_met, culture_supports_uti,
        cfu_ge_1e3, cfu_ge_1e5, cfu_threshold_source
    )
write_csv(legacy_comparison, FILE_STATUS_MAP_LEGACY_COMPARISON)

classification_audit <- episode_tbl %>%
    select(
        Participant_id, Timepoint, tp_lab, Episode_ID, Batch, Event_type,
        UTI_definition_version, UTI_Status, UTI_binary, Not_UTI_subgroup,
        Infection_Status_old, Infection_Status_legacy,
        analysis_include_primary, analysis_exclusion_reason,
        duplicate_role, duplicate_of_participant_id, duplicate_of_tp_lab,
        allow_secondary_duplicate_qc, duplicate_use_note,
        manual_curation_applied, manual_curation_note,
        UTI_classification_confidence, UTI_classification_reason,
        urine_collection_method_raw, urine_collection_method_norm, catheter_rule,
        symptom_compatible_uti, symptom_rule_met,
        local_urinary_symptom_any, systemic_symptom_any,
        all_of(symptom_fields), culture_supports_uti, cfu_raw, cfu_raw_parsed,
        cfu_lower_bound, cfu_upper_bound, cfu_ge_1e3, cfu_ge_1e4, cfu_ge_1e5,
        beoord_cat, cfu_threshold_used_for_uti, cfu_threshold_source,
        raw_CFU_examples, raw_BEO_examples, raw_Pop_examples
    )
write_csv(classification_audit, FILE_UTI_CLASSIFICATION_AUDIT)

movement_table <- episode_tbl %>%
    count(Infection_Status_old, UTI_Status, Not_UTI_subgroup,
          changed_to_new_uti, changed_due_to_lower_cfu_threshold,
          changed_due_to_catheter_aware_sns, name = "n") %>%
    arrange(desc(changed_to_new_uti), Infection_Status_old, UTI_Status, Not_UTI_subgroup)
write_csv(movement_table, FILE_UTI_RECLASSIFICATION_MOVEMENT)

symptom_rule_audit <- episode_tbl %>%
    count(catheter_rule, symptom_rule_met, symptom_compatible_uti, UTI_Status,
          Not_UTI_subgroup, name = "n") %>%
    arrange(catheter_rule, symptom_rule_met, UTI_Status)
write_csv(symptom_rule_audit, FILE_UTI_SYMPTOM_RULE_AUDIT)

cfu_threshold_audit <- episode_tbl %>%
    count(cfu_ge_1e3, cfu_ge_1e4, cfu_ge_1e5, culture_supports_uti,
          cfu_threshold_source, Infection_Status_old, UTI_Status,
          Not_UTI_subgroup, name = "n") %>%
    arrange(desc(culture_supports_uti), cfu_threshold_source, Infection_Status_old, UTI_Status)
write_csv(cfu_threshold_audit, FILE_UTI_CFU_THRESHOLD_AUDIT)

msg("Old status counts:")
print(table(episode_tbl$Infection_Status_old, useNA = "ifany"))
msg("New primary status counts:")
print(table(episode_tbl$UTI_Status, useNA = "ifany"))
msg("New primary status counts after manual primary exclusions:")
print(table(filter_primary_analysis(episode_tbl)$UTI_Status, useNA = "ifany"))
manual_excluded <- episode_tbl %>% filter(!(analysis_include_primary %in% TRUE))
if (nrow(manual_excluded) > 0) {
    msg("Manual primary exclusions applied: %d row(s)", nrow(manual_excluded))
    print(manual_excluded %>% count(analysis_exclusion_reason, name = "n"))
}
msg("Not_UTI subgroup counts (restricted to Not_UTI rows):")
print(table(episode_tbl$Not_UTI_subgroup[episode_tbl$UTI_Status == "Not_UTI"], useNA = "ifany"))
msg("UTI rows with Not_UTI_subgroup intentionally blank: %d",
    sum(episode_tbl$UTI_Status == "UTI" & is.na(episode_tbl$Not_UTI_subgroup), na.rm = TRUE))

append_denominator_summary(
    episode_tbl,
    "00b_classify_episodes.R",
    "status_map",
    "clinical_episode",
    input_rds,
    "Primary UTI_Status uses catheter-aware S&S and >=1e3 CFU support; legacy Infection_Status retained side-by-side"
)
append_denominator_summary(
    filter_primary_analysis(episode_tbl),
    "00b_classify_episodes.R",
    "status_map_primary_included",
    "clinical_episode",
    input_rds,
    "Primary clinical denominator after manual curation exclusions; excluded rows remain in status_map.csv for audit"
)
write_uti_attrition_outputs(FILE_STATUS_MAP)

msg("✓ Done.")

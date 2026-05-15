#!/usr/bin/env Rscript
# ==============================================================================
# 00b_classify_episodes.R
# ==============================================================================
#
# GOAL:
#   Classify each clinical episode as UTI, ASB, or Negative based on culture
#   results and symptom data.  This produces the authoritative status_map.csv
#   that every downstream analysis depends on.
#
# WHY THIS SCRIPT EXISTS:
#   The distinction between UTI and ASB is the central clinical question of
#   this project.  In elderly nursing home residents, most bacteriuria is
#   asymptomatic (ASB).  This script applies the project's operational episode-classification rules
#   systematically, handles three data sources for symptom information
#   (Population field, S&S sheets, SNS tables), and assigns confidence levels
#   for sensitivity analyses.
#
# IMPORTANT TERMINOLOGY:
#   In this script, "culture-positive" does NOT mean any bacterial growth.
#   It means a qualifying positive culture under the current operational rule:
#     - explicit CFU category interpreted as >=10^5 CFU/mL, OR
#     - Beoordeling "+++" when CFU_Count is unavailable.
#
#   This is more conservative than the broader YELLOW protocol description,
#   where uropathogen identification may occur from >10^3 CFU for urine
#   samples and >10^4 CFU for dipslides. Therefore, this script's threshold
#   should be described as the analysis threshold used for status_map.csv.
#
# CLASSIFICATION RULES:
#   UTI      = qualifying culture-positive AND UTI-related S&S present
#   ASB      = qualifying culture-positive AND UTI-related S&S absent
#   Culture-positive, S&S unknown
#            = qualifying culture-positive AND S&S cannot be determined
#   Negative = no qualifying positive culture
#
#   In this analysis script, UTI classification requires BOTH:
#     1. Qualifying culture positivity under the script's operational threshold
#        (>=10^5 CFU/mL or Beoordeling +++)
#     2. UTI-related signs/symptoms present
#
#   Bacteriuria/culture positivity alone is not classified as UTI.
#
# INPUTS:
#   - results/clinical/intermediate/clinical_merged.rds  (from 00a_)
#
# OUTPUTS:
#   - results/clinical/status_map.csv    (Participant_id × Timepoint × Status)
#   - assembly_metadata.csv              (isolate-level metadata for genomics)
#
# DOWNSTREAM:
#   → Infection_Status is used by virtually every analysis script (02–25).
#   → The GLMM in 14_ models UTI vs ASB as the binary outcome.
#   → The VF pipeline (22–25) stratifies all analyses by this status.
# ==============================================================================

source("00_config.R")
source(here::here("R", "clinical_helpers.R"))

suppressPackageStartupMessages({
    library(dplyr)
    library(readr)
    library(stringr)
    library(tidyr)
    library(Biostrings)
})

msg("Starting 00b_classify_episodes.R")

# Load Intermediate Data
input_rds <- file.path(DIR_CLINICAL, "intermediate", "clinical_merged.rds")
if (!file.exists(input_rds)) stop("Intermediate data not found. Run 00a_load_clean_clinical.R first.")
data_list <- readRDS(input_rds)

metadata <- data_list$metadata
sym_all_raw <- data_list$sym_all_raw
sns_12 <- data_list$sns_12

# 1. Process S&S Columns
# ------------------------------------------------------------------------------
if (nrow(sym_all_raw) == 0) {
    sym_all_1row <- tibble(Participant_id = character(), Timepoint = character(), Sx_present = logical())
} else {
    sym_input_cols <- setdiff(names(sym_all_raw), c("Participant_id", "Timepoint"))
    sym_bool <- sym_all_raw %>%
        mutate(
            Participant_id = str_trim(as.character(Participant_id)),
            Timepoint = canon_tp(Timepoint)
        ) %>%
        mutate(
            across(
                all_of(sym_input_cols),
                ~ {
                    x <- str_to_lower(str_trim(as.character(.x)))
                    num <- suppressWarnings(readr::parse_number(x))
                    dplyr::case_when(
                        !is.na(num) ~ num > 0,
                        x %in% c("ja", "yes", "y", "true", "1", "present", "pos", "positive") ~ TRUE,
                        x %in% c("nee", "no", "n", "false", "0", "", "absent", "neg", "negative", NA_character_) ~ FALSE,
                        TRUE ~ FALSE
                    )
                }
            )
        )
    sym_cols <- setdiff(names(sym_bool), c("Participant_id", "Timepoint"))
    sym_all_1row <- sym_bool %>%
        group_by(Participant_id, Timepoint) %>%
        summarise(Sx_present = any(rowSums(across(all_of(sym_cols), ~.x)) > 0), .groups = "drop") %>%
        distinct()
}

# 2. Classification Logic
# ------------------------------------------------------------------------------
meta_plus <- metadata %>%
    left_join(sns_12,
        by = c("Participant_id", "Timepoint", "Batch"),
        relationship = "many-to-one"
    ) %>%
    left_join(sym_all_1row %>% dplyr::rename(Sx_present_sym = Sx_present),
        by = c("Participant_id", "Timepoint"),
        relationship = "many-to-one"
    )

if (!"Sx_present_sym" %in% names(meta_plus)) meta_plus$Sx_present_sym <- NA

classified <- meta_plus %>%
    mutate(
        Beoord_clean = if ("Beoord" %in% names(.)) str_to_lower(str_trim(Beoord)) else NA_character_,

        # CFU: explicit > / ≥ / "meer dan 100000" => positive; ranges not positive
        # Use helper function cfu_bucket
        cfu_cat = if ("CFU_Count" %in% names(.)) cfu_bucket(CFU_Count) else NA_character_,
        cfu_recorded = if ("CFU_Count" %in% names(.)) (!is.na(CFU_Count) & str_trim(as.character(CFU_Count)) != "") else FALSE,
        beoord_cat = case_when(
            str_detect(ifelse(is.na(Beoord_clean), "", Beoord_clean), "\\+\\+\\+") ~ "+++",
            str_detect(ifelse(is.na(Beoord_clean), "", Beoord_clean), "\\+\\+") ~ "++",
            str_detect(ifelse(is.na(Beoord_clean), "", Beoord_clean), "\\+") ~ "+",
            TRUE ~ NA_character_
        ),

        # [UTI] Culture positivity per classification rule
        # ≥10^5 CFU/mL threshold is the standard diagnostic cutoff for UTI in clean-catch
        # midstream urine (IDSA guidelines, Nicolle et al. 2005; Hooton et al. 2010).
        # In nursing home residents and catheterized patients, this threshold may overdiagnose
        # ASB due to chronic colonization. Some guidelines suggest pyuria (≥10 WBC/HPF) as
        # additional criterion (Nicolle 2014, JAMA).
        culture_pos = case_when(
            cfu_recorded ~ (cfu_cat == ">=1e5"),
            !cfu_recorded & !is.na(beoord_cat) ~ (beoord_cat == "+++"),
            TRUE ~ NA
        ),

        # Population-derived S&S (authoritative when determinable, ALL batches)
        sx_present_pop = if ("Population" %in% names(.)) population_to_sns(Population) else as.logical(NA),

        # FINAL S&S decision (Population first, then explicit S&S columns)
        Sx_present_final = case_when(
            !is.na(sx_present_pop) ~ sx_present_pop,
            !is.na(Sx_present_sym) ~ Sx_present_sym,
            Batch %in% c(1L, 2L) & is.na(Sx_present_sym) &
                !is.na(SnS_status) ~ case_when(
                SnS_status == 2L ~ TRUE,
                SnS_status == 0L ~ FALSE,
                TRUE ~ as.logical(NA) # 1 or other -> unknown
            ),
            TRUE ~ as.logical(NA)
        ),

        # FINAL status
        Infection_Status = case_when(
            culture_pos == TRUE & Sx_present_final %in% TRUE ~ "UTI",
            culture_pos == TRUE & Sx_present_final %in% FALSE ~ "ASB",
            culture_pos == TRUE & is.na(Sx_present_final) ~ "Culture-positive, S&S unknown",
            culture_pos == FALSE ~ "Negative",
            TRUE ~ "None"
        ),

        # [EPI][STAT] Add confidence and provenance flags for sensitivity analyses
        Status_Confidence = case_when(
            # High confidence: Population field (clinician-determined) + explicit CFU
            !is.na(sx_present_pop) & cfu_recorded ~ "High",
            # Medium: S&S columns or SnS_status + explicit CFU
            is.na(sx_present_pop) & !is.na(Sx_present_final) & cfu_recorded ~ "Medium",
            # Low: Beoord-only CFU or missing S&S
            !cfu_recorded | is.na(Sx_present_final) ~ "Low",
            TRUE ~ "Unknown"
        ),
        Status_Provenance = case_when(
            !is.na(sx_present_pop) ~ "Population_field",
            !is.na(Sx_present_sym) ~ "SxS_columns",
            Batch %in% c(1L, 2L) & !is.na(SnS_status) ~ "SnS_status_fallback",
            TRUE ~ "Unknown"
        )
    )

# 3. Episode-level Collapse
# ------------------------------------------------------------------------------
# We collapse the raw clinical data down to one row per Participant_id + Timepoint.
# This is crucial because the scientific unit of interest is the clinical episode,
# not the individual isolate or assembly file.
# Without this step, participants with multiple isolates sequenced at the same
# timepoint would contribute more heavily to ASB/UTI comparisons, biasing the summaries.
if (!"Collection_Date" %in% names(classified)) classified$Collection_Date <- NA_character_
if (!"UTI_Label" %in% names(classified)) classified$UTI_Label <- NA_character_
if (!"Urine collection method" %in% names(classified)) classified[["Urine collection method"]] <- NA_character_

episode_tbl <- classified %>%
    group_by(Participant_id, Timepoint) %>%
    summarise(
        Batch = paste(sort(unique(Batch)), collapse = ";"),
        Collection_Date = paste(sort(unique(na.omit(as.character(Collection_Date)))), collapse = ";"),
        UTI_Label = paste(sort(unique(na.omit(as.character(UTI_Label)))), collapse = ";"),
        Urine_collection_method = paste(sort(unique(na.omit(as.character(`Urine collection method`)))), collapse = ";"),
        culture_pos_epi = any(culture_pos, na.rm = TRUE),
        cfu_recorded_any = any(ifelse(is.na(cfu_recorded), FALSE, cfu_recorded)),
        cfu_ge_1e5_any = any(cfu_cat == ">=1e5", na.rm = TRUE),
        beoord_plus3_any = any(beoord_cat == "+++", na.rm = TRUE),

        # tri-state collapse for S&S at episode level
        Sx_present_any = collapse_tristate(Sx_present_final),

        # Source labelling (Population first)
        from_pop = any(!is.na(sx_present_pop)),
        from_cols = any(!is.na(Sx_present_sym)),
        from_status = any(!is.na(SnS_status)),
        Sx_source_epi = case_when(
            from_pop ~ "Population",
            from_cols ~ "S&S columns",
            from_status ~ "SnS_status fallback",
            TRUE ~ NA_character_
        ),
        raw_CFU_examples = paste(head(unique(ifelse(is.na(CFU_Count), "", CFU_Count)), 5), collapse = " | "),
        raw_BEO_examples = paste(head(unique(ifelse(is.na(Beoord), "", Beoord)), 5), collapse = " | "),
        raw_Pop_examples = paste(head(unique(ifelse(is.na(Population), "", Population)), 5), collapse = " | "),
        .groups = "drop"
    ) %>%
    select(-from_pop, -from_cols, -from_status) %>%
    mutate(
        Collection_Date = na_if(Collection_Date, ""),
        UTI_Label = na_if(UTI_Label, ""),
        Urine_collection_method = na_if(Urine_collection_method, ""),
        tp_lab = normalise_timepoint_preserve_events(Timepoint),
        Event_type = episode_event_type(tp_lab),
        Infection_Status = case_when(
            culture_pos_epi & Sx_present_any %in% TRUE ~ "UTI",
            culture_pos_epi & Sx_present_any %in% FALSE ~ "ASB",
            culture_pos_epi & is.na(Sx_present_any) ~ "Culture-positive, S&S unknown",
            !culture_pos_epi ~ "Negative",
            TRUE ~ "None"
        ),
        # [EPI] Episode-level confidence (conservative: lowest confidence wins)
        Status_Confidence_epi = case_when(
            Sx_source_epi == "Population" & cfu_recorded_any ~ "High",
            Sx_source_epi %in% c("S&S columns", "SnS_status fallback") & cfu_recorded_any ~ "Medium",
            TRUE ~ "Low"
        )
    ) %>%
    mutate(
        Episode_ID = build_episode_id(., timepoint_col = "tp_lab", event_col = "Event_type",
                                      date_col = "Collection_Date")
    )

# Save Status Map
# Save Status Map
out_status <- file.path(DIR_CLINICAL, "status_map.csv")
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

write_csv(episode_tbl, out_status)
msg("Saved status_map.csv (%d rows) to %s", nrow(episode_tbl), out_status)

append_denominator_summary(
    episode_tbl,
    "00b_classify_episodes.R",
    "status_map",
    "clinical_episode",
    input_rds,
    "Episode_ID preserves routine, Uricult, and UTI-* event labels; duplicate keys are written to QC if present"
)
write_uti_attrition_outputs(out_status)

# (Section removed. FASTA Discovery & Metrics is now handled exclusively by 00_make_assembly_metadata.r)
msg("✓ Done.")

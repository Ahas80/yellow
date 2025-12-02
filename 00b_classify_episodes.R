#!/usr/bin/env Rscript
# ==============================================================================
# 00b_classify_episodes.R
# ------------------------------------------------------------------------------
# Role: [Data Prep] - Classify episodes as UTI, ASB, or Negative.
#
# Inputs:
#   - results/clinical/intermediate/clinical_merged.rds
#
# Outputs:
#   - results/clinical/status_map.csv
#   - assembly_metadata.csv (Project Root)
#
# Usage:
#   Rscript 00b_classify_episodes.R
#
# Biological/Statistical purpose:
#   - Applies clinical rules (CFU thresholds + Symptoms) to define the infection
#     status (UTI/ASB/Negative) for each episode.
#   - Assigns confidence levels for sensitivity analyses.
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
    sym_bool <- sym_all_raw %>%
        mutate(across(where(is.character), ~ str_to_lower(str_trim(.x)))) %>%
        mutate(
            across(
                -c(Participant_id, Timepoint),
                ~ {
                    x <- .x
                    num <- suppressWarnings(readr::parse_number(x))
                    dplyr::case_when(
                        !is.na(num) ~ num > 0,
                        x %in% c("ja", "yes", "y", "true", "1", "present", "pos", "positive") ~ TRUE,
                        x %in% c("nee", "no", "n", "false", "0", "", "absent", "neg", "negative", NA_character_) ~ FALSE,
                        TRUE ~ FALSE
                    )
                }
            ),
            Timepoint = canon_tp(Timepoint)
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
        relationship = "many-to-many"
    ) %>%
    left_join(sym_all_1row %>% dplyr::rename(Sx_present_sym = Sx_present),
        by = c("Participant_id", "Timepoint"),
        relationship = "many-to-many"
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

        # FINAL S&S decision (Population first)
        Sx_present_final = case_when(
            !is.na(sx_present_pop) ~ sx_present_pop,
            Batch %in% c(1L, 2L) & !is.na(Sx_present_sym) ~ Sx_present_sym,
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
            Batch %in% c(1L, 2L) & !is.na(Sx_present_sym) ~ "SxS_columns",
            Batch %in% c(1L, 2L) & !is.na(SnS_status) ~ "SnS_status_fallback",
            TRUE ~ "Unknown"
        )
    )

# 3. Episode-level Collapse
# ------------------------------------------------------------------------------
episode_tbl <- classified %>%
    group_by(Participant_id, Timepoint) %>%
    summarise(
        Batch = paste(sort(unique(Batch)), collapse = ","),
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
    )

# Save Status Map
# Save Status Map
out_status <- file.path(DIR_CLINICAL, "status_map.csv")
write_csv(episode_tbl, out_status)
msg("Saved status_map.csv (%d rows) to %s", nrow(episode_tbl), out_status)

# 4. FASTA Discovery & Metrics
# ------------------------------------------------------------------------------
# This generates assembly_metadata.csv, linking Isolates to Participants/Timepoints
# It uses the metadata we just loaded to map Isolate_ID -> Participant_id

all_fasta <- list.files(DIR_FASTAS, "\\.fasta$", full.names = TRUE)
if (length(all_fasta) == 0) {
    warning("No FASTA files found in ", DIR_FASTAS)
} else {
    assembly_tbl <- tibble(full_path = all_fasta) %>%
        mutate(
            file_name = basename(full_path),
            Isolate_ID = str_extract(str_split_fixed(file_name, "_", 4)[, 3], "[0-9A-Za-z]+-[0-9]+"),
            assembler = case_when(
                str_detect(file_name, "flye") ~ "flye",
                str_detect(file_name, "longcycler") ~ "longcycler",
                TRUE ~ "unknown"
            )
        )

    summarise_fasta <- function(fp) {
        x <- readDNAStringSet(fp)
        af <- colSums(alphabetFrequency(x, baseOnly = TRUE))
        tibble(
            num_contigs = length(x),
            total_bases = sum(width(x)),
            gc_content  = round((af["G"] + af["C"]) / sum(af) * 100, 2)
        )
    }

    msg("Calculating metrics for %d assemblies...", nrow(assembly_tbl))
    metrics_tbl <- purrr::map_dfr(assembly_tbl$full_path, summarise_fasta)

    assembly_df <- assembly_tbl %>%
        select(file_name, full_path, Isolate_ID, assembler) %>%
        bind_cols(metrics_tbl) %>%
        left_join(metadata, by = c("Isolate_ID" = "isolate_ID"), relationship = "many-to-many")

    write_csv(assembly_df, "assembly_metadata.csv")
    msg("Saved assembly_metadata.csv (%d rows)", nrow(assembly_df))
}

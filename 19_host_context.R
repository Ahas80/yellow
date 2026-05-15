#!/usr/bin/env Rscript
# ==============================================================================
# 19_host_context.R
# ==============================================================================
#
# GOAL:
#   Extract clinical and host context (symptoms, catheter use, etc.) specifically
#   for the episodes involved in phenotype switching.  This bridges the gap
#   between the genomic events (from script 18) and the clinical reality.
#
# WHY THIS SCRIPT EXISTS:
#   Even if a strain acquires a new VF gene, the transition to symptomatic
#   UTI might be driven primarily by a host factor change (e.g., insertion
#   of a urinary catheter).  We need to review the host context alongside
#   the genomic context to avoid incorrectly associating genomic changes with
#   clinical outcomes when a massive host confounder is present.
#
# ------------------------------------------------------------------------------
# Role: [Analysis] - Priority 2: Host Context Integration.
#
# Inputs:
#   - results/longitudinal/phenotype_switch_candidates.csv
#   - results/clinical/intermediate/clinical_merged.rds
#
# Outputs:
#   - results/longitudinal/host_context_table.csv
#
# Purpose:
#   - Extract "Urine collection method" (Catheter use) and Symptoms for
#     the specific episodes involved in phenotype switching.
#   - (Antibiotics data is currently missing from inputs).
# ==============================================================================

source("00_config.R")
suppressPackageStartupMessages({
    library(dplyr)
    library(readr)
    library(tidyr)
    library(stringr)
})

msg("Starting 19_host_context.R")

# 1. Load Candidates
# ------------------------------------------------------------------------------
cand_file <- file.path(DIR_RESULTS, "longitudinal", "phenotype_switch_candidates.csv")
if (!file.exists(cand_file)) stop("No candidates file found.")
candidates <- read_csv(cand_file, show_col_types = FALSE)

# 2. Load Clinical Data
# ------------------------------------------------------------------------------
rds_file <- file.path(DIR_CLINICAL, "intermediate", "clinical_merged.rds")
if (!file.exists(rds_file)) stop("Clinical RDS not found.")
data_list <- readRDS(rds_file)

# We need metadata (for collection method) and sym_all_raw (for specific symptoms)
meta <- data_list$metadata
syms <- data_list$sym_all_raw

# 3. Extract Context for Candidates
# ------------------------------------------------------------------------------
# We need to look up T1 and T2 for each candidate.
# Candidate columns: Participant_id, From_Time, To_Time, From_Status, To_Status

get_context <- function(pid, tp) {
    # 1. Collection Method
    m_row <- meta %>% filter(Participant_id == pid, Timepoint == tp)
    coll_method <- if (nrow(m_row) > 0) m_row$`Urine collection method`[1] else NA_character_

    # 2. Symptoms (Detailed)
    # sym_all_raw has columns like "S&S Fever", "S&S Dysuria", etc.
    s_row <- syms %>% filter(Participant_id == pid, Timepoint == tp)

    symptom_str <- "None"
    if (nrow(s_row) > 0) {
        # Find columns starting with "S&S" that are TRUE or "1" or "Checked"
        # The raw data might be messy. Let's look for known keywords.

        # Helper to check if a column indicates presence
        is_present <- function(x) {
            !is.na(x) & (x == "TRUE" | x == "1" | x == "Checked" | x == "Yes" | x == "x")
        }

        active_syms <- c()
        for (col in names(s_row)) {
            # Match keywords used in 00a_load_clean_clinical.R
            if (grepl("sympt|dysur|urg|urge|fever|koorts|pijn|pain|burn|S&S", col, ignore.case = TRUE)) {
                val <- s_row[[col]][1]
                if (is_present(val)) {
                    # Clean column name: Remove "S&S " prefix if present, or just use col name
                    clean_name <- gsub("S&S |S&S_", "", col)
                    active_syms <- c(active_syms, clean_name)
                }
            }
        }
        if (length(active_syms) > 0) symptom_str <- paste(active_syms, collapse = "; ")
    }

    tibble(
        Collection_Method = coll_method,
        Symptoms_Detail = symptom_str
    )
}

results <- list()

for (i in 1:nrow(candidates)) {
    row <- candidates[i, ]

    # Context for "From" timepoint
    ctx_A <- get_context(row$Participant_id, row$From_Time)
    names(ctx_A) <- paste0("From_", names(ctx_A))

    # Context for "To" timepoint
    ctx_B <- get_context(row$Participant_id, row$To_Time)
    names(ctx_B) <- paste0("To_", names(ctx_B))

    # Combine
    combined <- bind_cols(row, ctx_A, ctx_B)
    results[[i]] <- combined
}

final_table <- bind_rows(results)

# 4. Interpret Catheter Status
# ------------------------------------------------------------------------------
# CAD = Catheter à Demeure (Indwelling)
final_table <- final_table %>%
    mutate(
        From_Catheter = grepl("CAD", From_Collection_Method, ignore.case = TRUE),
        To_Catheter = grepl("CAD", To_Collection_Method, ignore.case = TRUE),
        Catheter_Change = case_when(
            From_Catheter & !To_Catheter ~ "Removed",
            !From_Catheter & To_Catheter ~ "Inserted",
            From_Catheter & To_Catheter ~ "Maintained",
            TRUE ~ "None"
        )
    )

# 5. Save
# ------------------------------------------------------------------------------
out_file <- file.path(DIR_RESULTS, "longitudinal", "host_context_table.csv")
write_csv(final_table, out_file)
msg("Saved host context to %s", out_file)

# Print Summary for Log
print(final_table %>% select(Participant_id, From_Status, To_Status, From_Collection_Method, To_Collection_Method, Catheter_Change))

msg("Done.")

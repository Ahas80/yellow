# Script: get_stratified_stats.R
# Goal: Generate side-by-side stats for >=2 TPs and >=3 TPs cohorts
# To populate the two abstract options requested.

suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(tidyr))

# Load Data - Using status_map.csv as it works for the verified stats (has cfu_recorded_any)
status_df <- read_csv("status_map.csv", show_col_types = FALSE) %>%
    mutate(Participant_id = as.character(Participant_id))

# Load Genomic Data (for STs)
mlst_df <- read_csv("results/mlst/mlst_with_meta.csv", show_col_types = FALSE) %>%
    mutate(Participant_id = as.character(Participant_id))

# Clean & Filter Clinical Data
# LOGIC:
# 1. Count Total Timepoints per Patient (Retention)
# 2. THEN Filter for episodes with CFU data for analysis
clean_df <- status_df %>%
    filter(!is.na(Infection_Status)) %>%
    group_by(Participant_id) %>%
    mutate(n_tp = n_distinct(Timepoint)) %>% # Count retention first
    ungroup() %>%
    filter(cfu_recorded_any == TRUE) # Filter for stats

# Function to calculate stats for a given threshold
calc_stats <- function(df, min_tps, label) {
    # 1. Filter Cohort
    cohort <- df %>% filter(n_tp >= min_tps)

    pids <- unique(cohort$Participant_id)
    n_pids <- length(pids)
    n_episodes <- nrow(cohort)

    # 2. Status Counts
    status_counts <- table(cohort$Infection_Status)
    status_pct <- prop.table(status_counts) * 100

    # 3. Transitions
    # Sort by TP
    trans_df <- cohort %>%
        mutate(
            tp_num_val = case_when(
                grepl("T[0-9]+", Timepoint) ~ as.numeric(gsub("T", "", Timepoint)),
                Timepoint == "Uricult" ~ 99,
                TRUE ~ 999
            )
        ) %>%
        group_by(Participant_id) %>%
        arrange(tp_num_val) %>%
        mutate(Next_Status = lead(Infection_Status)) %>%
        filter(!is.na(Next_Status))

    n_trans <- nrow(trans_df)
    asb_trans <- trans_df %>% filter(Infection_Status == "ASB")
    asb_to_asb <- sum(asb_trans$Next_Status == "ASB")
    asb_to_uti <- sum(asb_trans$Next_Status == "UTI")
    asb_to_neg <- sum(asb_trans$Next_Status %in% c("Negative", "Cult_Neg"))

    neg_trans <- trans_df %>% filter(Infection_Status %in% c("Negative", "Cult_Neg"))
    neg_to_asb <- sum(neg_trans$Next_Status == "ASB")
    neg_to_uti <- sum(neg_trans$Next_Status == "UTI")
    neg_to_neg <- sum(neg_trans$Next_Status %in% c("Negative", "Cult_Neg"))


    # 4. Genomic Stats (Subset of Cohort with MLST)
    class_inputs <- read_csv("results/class_inputs_full.csv", show_col_types = FALSE) %>%
        mutate(Participant_id = as.character(Participant_id)) %>%
        select(Participant_id, tp_lab, isolate_ID)

    cohort_isolates <- cohort %>%
        left_join(class_inputs, by = c("Participant_id", "Timepoint" = "tp_lab"), relationship = "many-to-many")

    cohort_mlst <- cohort_isolates %>%
        left_join(mlst_df, by = c("isolate_ID" = "Isolate_ID"), relationship = "many-to-many") %>%
        rename(
            Participant_id = Participant_id.x,
            Timepoint = Timepoint.x
        ) %>%
        filter(!is.na(ST)) %>%
        mutate(
            tp_num_val = case_when(
                grepl("T[0-9]+", Timepoint) ~ as.numeric(gsub("T", "", Timepoint)),
                Timepoint == "Uricult" ~ 99,
                TRUE ~ 999
            )
        )

    # Top STs
    st_counts <- table(cohort_mlst$ST)
    st_counts_df <- as.data.frame(st_counts) %>%
        arrange(desc(Freq)) %>%
        mutate(Pct = round(Freq / sum(Freq) * 100, 1))

    # Distinct STs per status
    st_in_uti <- cohort_mlst %>%
        filter(Infection_Status == "UTI") %>%
        pull(ST) %>%
        n_distinct()
    st_in_asb <- cohort_mlst %>%
        filter(Infection_Status == "ASB") %>%
        pull(ST) %>%
        n_distinct()

    # Stability (Consecutive Pairs)
    # The prompt refers to "consecutive isolate pairs".

    st_trans_df <- cohort_mlst %>%
        distinct(Participant_id, tp_num_val, ST, Infection_Status) %>%
        group_by(Participant_id) %>%
        arrange(tp_num_val) %>%
        mutate(Next_ST = lead(ST)) %>%
        filter(!is.na(Next_ST))

    n_st_pairs <- nrow(st_trans_df)
    n_stable <- sum(st_trans_df$ST == st_trans_df$Next_ST)
    pct_stable <- round(n_stable / n_st_pairs * 100, 1)

    # Switchers (Participants with at least one switch)
    switchers <- st_trans_df %>%
        group_by(Participant_id) %>%
        summarize(Any_Switch = any(ST != Next_ST))

    n_switchers <- sum(switchers$Any_Switch)
    pct_switchers <- round(n_switchers / nrow(switchers) * 100, 1)

    # Output
    cat("\n---", label, "---\n")
    cat("Participants (Clinical):", n_pids, "\n")
    cat("Episodes with CFU:", n_episodes, "\n")
    cat("Episodes with MLST:", nrow(cohort_mlst), "(", round(nrow(cohort_mlst) / n_episodes * 100, 1), "%)\n")

    cat("\nStatus Distribution:\n")
    print(status_counts)
    print(round(status_pct, 1))

    cat("\nTotal Transitions (Clinical):", n_trans, "\n")

    cat("\n--- ST STATS ---\n")
    cat("Distinct STs:", n_distinct(cohort_mlst$ST), "\n")
    cat("Distinct STs in UTI:", st_in_uti, "\n")
    cat("Distinct STs in ASB:", st_in_asb, "\n")

    cat("Top 5 STs:\n")
    print(head(st_counts_df, 5))

    cat("\nST Stability (Consecutive Pairs):\n")
    cat("Pairs:", n_st_pairs, "\n")
    cat("Stable Pairs:", n_stable, "(", pct_stable, "%)\n")

    cat("Participants with Switch:", n_switchers, "/", nrow(switchers), "(", pct_switchers, "%)\n")
}

calc_stats(clean_df, 2, "Cohort >= 2 Timepoints + CFU")
calc_stats(clean_df, 3, "Cohort >= 3 Timepoints + CFU")
calc_stats(clean_df, 4, "Cohort >= 4 Timepoints + CFU")

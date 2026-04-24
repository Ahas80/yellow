#!/usr/bin/env Rscript

# generate_abstract_stats.R
# Recompute publication-grade stats for Abstract
# STRICT adherence to file-derived counts.

suppressPackageStartupMessages({
    library(tidyverse)
    library(lubridate)
    library(glue)
})

# --- CONFIGURATION ---
FILE_STATUS <- "results/clinical/status_map.csv"
FILE_MLST <- "results/mlst/mlst_with_meta.csv"
FILE_VF <- "results/vf/vf_hits_all.rds"

ALLELE_COLS <- c("dinb", "icda", "pabb", "polb", "putp", "trpa", "trpb", "uida")

# --- HELPER FUNCTIONS ---
print_header <- function(title) {
    cat("\n# ", title, "\n\n", sep = "")
}

print_sub_header <- function(title) {
    cat("\n## ", title, "\n\n", sep = "")
}

# --- MAIN ---

print_header("SECTION 1: Data Integrity & Linkage")

# 1.1 Check files
files_to_check <- c(FILE_STATUS, FILE_MLST)
file_info <- file.info(files_to_check) %>%
    rownames_to_column("File") %>%
    select(File, mtime, size)

print_sub_header("1.1 File Timestamps")
print(knitr::kable(file_info))

if (any(is.na(file_info$size))) stop("CRITICAL: Missing required files.")

# 1.2 Load Data
df_status <- read_csv(FILE_STATUS, show_col_types = FALSE) %>%
    mutate(Participant_id = as.character(Participant_id))

df_mlst <- read_csv(FILE_MLST, show_col_types = FALSE) %>%
    mutate(Participant_id = as.character(Participant_id))

# Verify Columns
req_cols_status <- c("Participant_id", "Timepoint", "Infection_Status")
req_cols_mlst <- c("Participant_id", "Timepoint", "ST", "Isolate_ID", "has_new_allele", "ambiguous_call")

missing_status <- setdiff(req_cols_status, names(df_status))
missing_mlst <- setdiff(req_cols_mlst, names(df_mlst))

if (length(missing_status) > 0) stop(glue("Missing columns in status_map: {paste(missing_status, collapse=', ')}"))
if (length(missing_mlst) > 0) stop(glue("Missing columns in mlst: {paste(missing_mlst, collapse=', ')}"))

print_sub_header("1.2 Row Counts Loaded")
cat(glue("Status Map: {nrow(df_status)} rows\n"))
cat(glue("MLST File:  {nrow(df_mlst)} rows\n"))

# 1.3 Linkage
# Prepare MLST for joining: Deduplicate per episode
# Rule: Prefer mlst_complete (no new, no ambig) -> highest n_loci -> first row

# Check which allele cols actually exist
actual_allele_cols <- intersect(names(df_mlst), ALLELE_COLS)

df_mlst_clean <- df_mlst %>%
    mutate(
        # Clean ST column (sometimes it has (?) or other marks, usually pure integer or -)
        ST_clean = suppressWarnings(as.numeric(as.character(ST))),

        # Calculate Completeness metrics
        has_new = replace_na(has_new_allele, FALSE),
        is_ambig = replace_na(ambiguous_call, FALSE),
        mlst_complete = !has_new & !is_ambig & !is.na(ST_clean),

        # Count typed loci
        n_loci_typed = rowSums(!is.na(select(., all_of(actual_allele_cols))))
    ) %>%
    group_by(Participant_id, Timepoint) %>%
    arrange(desc(mlst_complete), desc(n_loci_typed), Isolate_ID) %>% # Deterministic sort
    mutate(rank = row_number()) %>%
    ungroup()

# Use rank 1 for linkage
df_mlst_unique <- df_mlst_clean %>% filter(rank == 1)

collapsed_count <- nrow(df_mlst_clean) - nrow(df_mlst_unique)

# JOIN
# Note: status_map is the master for "Episodes"
df_linked <- df_status %>%
    left_join(df_mlst_unique %>% select(Participant_id, Timepoint, ST, ST_clean, Isolate_ID),
        by = c("Participant_id", "Timepoint")
    )

# Stats for 1.3
n_episodes_total <- nrow(df_status)
n_asb_uti <- df_linked %>%
    filter(Infection_Status %in% c("ASB", "UTI")) %>%
    nrow()
n_linked <- df_linked %>%
    filter(!is.na(ST)) %>%
    nrow()
n_asb_uti_linked <- df_linked %>%
    filter(Infection_Status %in% c("ASB", "UTI"), !is.na(ST)) %>%
    nrow()

print_sub_header("1.3 Linkage Yield")
cat(glue(" Total Episodes (status_map): {n_episodes_total}\n"))
cat(glue(" ASB/UTI Episodes: {n_asb_uti}\n"))
cat(glue(" Episodes with ST data linked: {n_linked} ({round(n_linked/n_episodes_total*100, 1)}%)\n"))
cat(glue(" ASB/UTI Episodes with ST linked: {n_asb_uti_linked} ({round(n_asb_uti_linked/n_asb_uti*100, 1)}%)\n"))
cat(glue(" Duplicate isolates collapsed: {collapsed_count}\n"))

print_sub_header("1.4 QC & Dropout Narrative")
n_missing_st <- n_asb_uti - n_asb_uti_linked
pct_missing <- round(n_missing_st / n_asb_uti * 100, 1)

cat(glue("> \"In the primary cohort (>=2 timepoints), we identified {n_asb_uti} ASB/UTI episodes. Sequence type could not be assigned for {n_missing_st} episodes ({pct_missing}%) due to sequencing ambiguity or QC failure. These were excluded from lineage analysis, leaving {n_asb_uti_linked} typeable isolates for detailed characterization.\"\n"))


# --- SECTION 2: COHORT SIZE ---
print_header("SECTION 2: Cohort Size by Follow-up Depth")

# Determine eligibility (>= K distinct timepoints in status_map)
ppt_counts <- df_status %>%
    group_by(Participant_id) %>%
    summarise(n_timepoints = n_distinct(Timepoint))

cohort_stats <- list()

for (K in c(2, 3, 4)) {
    eligible_ids <- ppt_counts %>%
        filter(n_timepoints >= K) %>%
        pull(Participant_id)

    # Filter Data
    subset_data <- df_linked %>% filter(Participant_id %in% eligible_ids)
    subset_asb_uti <- subset_data %>% filter(Infection_Status %in% c("ASB", "UTI"))

    # Calc Stats
    n_participants <- length(eligible_ids)
    n_epi_asb_uti <- nrow(subset_asb_uti)
    n_epi_linked <- sum(!is.na(subset_asb_uti$ST))

    # Timepoints range
    tp_per_ppt <- subset_data %>%
        group_by(Participant_id) %>%
        summarise(n = n()) %>%
        pull(n)
    tp_range <- glue("{min(tp_per_ppt)}-{max(tp_per_ppt)}")

    # ASB/UTI Episodes per participant
    epi_per_ppt <- subset_asb_uti %>%
        group_by(Participant_id) %>%
        summarise(n = n())
    epi_counts_all <- tibble(Participant_id = eligible_ids) %>%
        left_join(epi_per_ppt, by = "Participant_id") %>%
        mutate(n = replace_na(n, 0)) %>%
        pull(n)

    epi_median <- median(epi_counts_all)
    epi_range <- glue("{min(epi_counts_all)}-{max(epi_counts_all)}")

    cohort_stats[[paste0("K", K)]] <- tibble(
        Recall_Depth = glue(">={K} Timepoints"),
        Participants = n_participants,
        ASB_UTI_Episodes = n_epi_asb_uti,
        Linked_Episodes = n_epi_linked,
        TP_Range = tp_range,
        Epi_Per_Ppt_Median_Range = glue("{epi_median} ({epi_range})")
    )
}

df_cohort_stats <- bind_rows(cohort_stats)
print(knitr::kable(df_cohort_stats))


# --- SECTION 3: ST DISTRIBUTION ---
print_header("SECTION 3: ST Distribution (ASB/UTI Episodes)")

for (K in c(2, 3, 4)) {
    print_sub_header(glue("K>={K} Subset"))

    eligible_ids <- ppt_counts %>%
        filter(n_timepoints >= K) %>%
        pull(Participant_id)

    # Analysis Set: Eligible Participants + ASB/UTI Status + Has ST
    analysis_set <- df_linked %>%
        filter(
            Participant_id %in% eligible_ids,
            Infection_Status %in% c("ASB", "UTI"),
            !is.na(ST)
        )

    total_linked <- nrow(analysis_set)
    total_ppts <- n_distinct(analysis_set$Participant_id)

    st_counts <- analysis_set %>%
        group_by(ST) %>%
        summarise(
            Episodes = n(),
            Participants = n_distinct(Participant_id)
        ) %>%
        arrange(desc(Episodes)) %>%
        mutate(
            Pct_Episodes = round(Episodes / total_linked * 100, 1),
            Pct_Participants = round(Participants / length(eligible_ids) * 100, 1)
        ) %>%
        head(10)

    print(knitr::kable(st_counts))

    distinct_sts <- n_distinct(analysis_set$ST)
    top5_sum <- st_counts %>%
        slice(1:5) %>%
        summarise(s = sum(Episodes)) %>%
        pull(s)
    top5_pct <- round(top5_sum / total_linked * 100, 1)

    cat(glue("\n Distinct STs observed: {distinct_sts}"))
    cat(glue("\n Top 5 STs comprise {top5_pct}% of all linked episodes.\n"))
}

# --- SECTION 4: UTI PROPORTION ---
print_header("SECTION 4: UTI Proportion by ST")

for (K in c(2, 3, 4)) {
    print_sub_header(glue("K>={K} Subset"))

    eligible_ids <- ppt_counts %>%
        filter(n_timepoints >= K) %>%
        pull(Participant_id)

    analysis_set <- df_linked %>%
        filter(
            Participant_id %in% eligible_ids,
            Infection_Status %in% c("ASB", "UTI"),
            !is.na(ST)
        )

    st_risk <- analysis_set %>%
        group_by(ST) %>%
        summarise(
            n_total = n(),
            n_UTI = sum(Infection_Status == "UTI"),
            n_ASB = sum(Infection_Status == "ASB")
        ) %>%
        filter(n_total >= 5) %>%
        mutate(
            UTI_Prop = round(n_UTI / n_total * 100, 1)
        ) %>%
        arrange(desc(UTI_Prop))

    print(knitr::kable(st_risk))

    # Abstract sentence
    if (nrow(st_risk) > 0) {
        top_st <- st_risk %>% slice(1)
        cat(glue("\n> Among STs with >=5 episodes, ST{top_st$ST} had the highest UTI proportion ({top_st$n_UTI}/{top_st$n_total} = {top_st$UTI_Prop}%), compared with... \n"))
    } else {
        cat("\n> No STs with >= 5 episodes in this subset.\n")
    }
}

# --- SECTION 5: WITHIN-PARTICIPANT STABILITY ---
print_header("SECTION 5: Within-Participant Stability")

for (K in c(2, 3, 4)) {
    print_sub_header(glue("K>={K} Subset"))

    eligible_ids <- ppt_counts %>%
        filter(n_timepoints >= K) %>%
        pull(Participant_id)

    long_data <- df_linked %>%
        filter(
            Participant_id %in% eligible_ids,
            Infection_Status %in% c("ASB", "UTI"),
            !is.na(ST)
        ) %>%
        mutate(TP_Num = as.numeric(str_extract(Timepoint, "\\d+"))) %>%
        arrange(Participant_id, TP_Num) %>%
        group_by(Participant_id) %>%
        mutate(
            Prev_ST = lag(ST),
            Is_Switch = ST != Prev_ST
        ) %>%
        summarise(
            n_episodes = n(),
            n_switches = sum(Is_Switch, na.rm = TRUE),
            is_stable = n_switches == 0 & n_episodes > 1
        )

    long_valid <- long_data %>% filter(n_episodes >= 2)
    if (nrow(long_valid) == 0) {
        cat("No participants with >=2 linked episodes in this subset.\n")
        next
    }

    n_stable <- sum(long_valid$is_stable)
    n_switcher <- sum(long_valid$n_switches > 0)
    pct_stable <- round(n_stable / nrow(long_valid) * 100, 1)
    pct_switcher <- round(n_switcher / nrow(long_valid) * 100, 1)

    cat(glue("Participants with >=2 linked episodes: {nrow(long_valid)}\n"))
    cat(glue("Stable ST: {n_stable} ({pct_stable}%)\n"))
    cat(glue("At least one switch: {n_switcher} ({pct_switcher}%)\n"))
    cat(glue("Median switches (IQR): {median(long_valid$n_switches)} ({IQR(long_valid$n_switches)})\n"))

    transitions <- df_linked %>%
        filter(
            Participant_id %in% eligible_ids,
            Infection_Status %in% c("ASB", "UTI"),
            !is.na(ST)
        ) %>%
        mutate(TP_Num = as.numeric(str_extract(Timepoint, "\\d+"))) %>%
        arrange(Participant_id, TP_Num) %>%
        group_by(Participant_id) %>%
        mutate(
            Prev_ST = lag(ST)
        ) %>%
        filter(!is.na(Prev_ST), ST != Prev_ST) %>%
        ungroup() %>%
        count(Prev_ST, ST, sort = TRUE) %>%
        head(10)

    print(knitr::kable(transitions, caption = "Top ST Switches"))
}

# --- SECTION 6: VF BURDEN ---
print_header("SECTION 6: VF Burden Analysis")

if (file.exists(FILE_VF)) {
    df_vf <- readRDS(FILE_VF)

    if (is.data.frame(df_vf)) {
        cat("Loaded VF data. Rows: ", nrow(df_vf), "\n")

        # Check for GENE and Isolate_ID columns
        if ("GENE" %in% names(df_vf) && "Isolate_ID" %in% names(df_vf)) {
            burden_df <- df_vf %>%
                group_by(Isolate_ID) %>%
                summarise(vf_burden = n_distinct(GENE)) %>%
                ungroup()

            # Join to Linked Data (All Cohorts combined for Top STs)
            df_vf_linked <- df_linked %>%
                filter(Infection_Status %in% c("ASB", "UTI"), !is.na(ST)) %>%
                inner_join(burden_df, by = "Isolate_ID")

            # Identify Top 10 STs from the linked data
            top10_sts <- df_vf_linked %>%
                count(ST, sort = TRUE) %>%
                head(10) %>%
                pull(ST)

            vf_stats <- df_vf_linked %>%
                filter(ST %in% top10_sts) %>%
                group_by(ST) %>%
                summarise(
                    N_Episodes = n(),
                    Median_Burden = median(vf_burden),
                    IQR_Burden = IQR(vf_burden),
                    UTI_Burden_Med = if (sum(Infection_Status == "UTI") > 0) median(vf_burden[Infection_Status == "UTI"]) else NA,
                    ASB_Burden_Med = if (sum(Infection_Status == "ASB") > 0) median(vf_burden[Infection_Status == "ASB"]) else NA
                )

            print(knitr::kable(vf_stats))
        } else {
            cat("Missing GENE or Isolate_ID column in VF file. Skipping.\n")
            cat("Columns found: ", paste(names(df_vf), collapse = ", "), "\n")
        }
    } else {
        cat("VF file is not a dataframe. Skipping.\n")
    }
} else {
    cat("VF File not found. Skipping Section 6.\n")
}

# --- FINAL OUTPUT ---
print_header("FINAL DELIVERABLE: Abstract Numbers")
cat("Run the script to see the generated lists above. Copy specific rows as needed.\n")
cat("Methods One-Liner:\n")
cat(glue("'We analyzed E. coli ST dynamics in {n_participants} participants with >=X follow-up episodes. STs were determined from whole-genome assemblies and linked to clinical episodes (UTI/ASB) via timepoint-matched sequence typing.'\n"))

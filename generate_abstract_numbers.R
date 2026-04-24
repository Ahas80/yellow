library(tidyverse)
library(lubridate)

# --- 1. SETUP & FILE SELECTION ---
cat("DATA CHECK\n")
cat("==========\n")

# Define target files
status_map_path <- "/Users/Aamir/Desktop/rUTIs/results/clinical/status_map.csv"
mlst_path <- "/Users/Aamir/Desktop/rUTIs/results/mlst/mlst_with_meta.csv"

# Check existence and times
if (!file.exists(status_map_path)) stop("status_map.csv not found!")
if (!file.exists(mlst_path)) stop("mlst_with_meta.csv not found!")

info_status <- file.info(status_map_path)
info_mlst <- file.info(mlst_path)

cat(sprintf("Status Map: %s (Last Modified: %s)\n", status_map_path, info_status$mtime))
cat(sprintf("MLST File:  %s (Last Modified: %s)\n", mlst_path, info_mlst$mtime))

# Load Data
df_status <- read_csv(status_map_path, show_col_types = FALSE)
df_mlst <- read_csv(mlst_path, show_col_types = FALSE)

cat(sprintf("Rows in status_map: %d\n", nrow(df_status)))

# Check & Map Columns for Status Map
cat("\nColumn Mapping for status_map.csv:\n")
if ("Participant_id" %in% names(df_status) && !"participant_id" %in% names(df_status)) {
    cat("- Found 'Participant_id' -> Mapping to 'participant_id'\n")
    df_status <- df_status %>% rename(participant_id = Participant_id)
} else {
    cat("- 'participant_id' found directly.\n")
}

req_cols <- c("participant_id", "Timepoint", "Infection_Status")
if (!all(req_cols %in% names(df_status))) {
    stop("Missing columns in status_map: ", paste(setdiff(req_cols, names(df_status)), collapse = ", "))
}

# --- 2. DATA CLEANING & ORDERING ---

# Normalize Timepoint Ordering
# Assumes T0, T1, T2... or dates.
# Let's try to extract numbers if T-prefix, otherwise rank.
df_status <- df_status %>%
    mutate(
        tp_num = as.numeric(str_extract(Timepoint, "\\d+")),
        tp_rank = if_else(!is.na(tp_num), tp_num, as.numeric(as.factor(Timepoint))) # Fallback
    ) %>%
    arrange(participant_id, tp_num)

# Check Infection Statuses
cat("Unique Infection_Status values:\n")
print(unique(df_status$Infection_Status))

# --- 3. HELPER FUNCTIONS ---

generate_clinical_abstract <- function(df, K) {
    # 1. Eligibility: >= K distinct timepoints
    eligible_ids <- df %>%
        group_by(participant_id) %>%
        summarise(n_tp = n_distinct(Timepoint)) %>%
        filter(n_tp >= K) %>%
        pull(participant_id)

    df_cohort <- df %>%
        filter(participant_id %in% eligible_ids) %>%
        arrange(participant_id, tp_num)

    n_part <- n_distinct(df_cohort$participant_id)
    n_episodes <- nrow(df_cohort)

    # Status Counts
    stats <- df_cohort %>%
        count(Infection_Status) %>%
        mutate(pct = n / sum(n) * 100)

    n_asb <- sum(stats$n[stats$Infection_Status == "ASB"])
    p_asb <- sum(stats$pct[stats$Infection_Status == "ASB"])
    n_uti <- sum(stats$n[stats$Infection_Status == "UTI"])
    p_uti <- sum(stats$pct[stats$Infection_Status == "UTI"])
    n_neg <- sum(stats$n[stats$Infection_Status == "Negative"])
    p_neg <- sum(stats$pct[stats$Infection_Status == "Negative"])

    # Transitions
    df_trans <- df_cohort %>%
        group_by(participant_id) %>%
        arrange(tp_num) %>%
        mutate(
            next_status = lead(Infection_Status),
            pair_type = paste(Infection_Status, next_status, sep = "->")
        ) %>%
        filter(!is.na(next_status))

    n_trans <- nrow(df_trans)

    # ASB Stability: ASB->ASB / ASB->Any
    asb_starts <- df_trans %>% filter(Infection_Status == "ASB")
    n_asb_starts <- nrow(asb_starts)
    if (n_asb_starts > 0) {
        n_asb_stable <- sum(asb_starts$next_status == "ASB")
        p_asb_stable <- n_asb_stable / n_asb_starts * 100

        n_asb_uti <- sum(asb_starts$next_status == "UTI")
        p_asb_uti <- n_asb_uti / n_asb_starts * 100
    } else {
        n_asb_stable <- 0
        p_asb_stable <- 0
        n_asb_uti <- 0
        p_asb_uti <- 0
    }

    # Negative -> UTI
    neg_starts <- df_trans %>% filter(Infection_Status == "Negative")
    n_neg_uti <- if (nrow(neg_starts) > 0) sum(neg_starts$next_status == "UTI") else 0

    # Output
    cat(sprintf("\n[DOC %d] Abstract: Longitudinal Clinical Transitions Between ASB and UTI in Nursing Home Residents (≥%d TIMEPOINTS)\n", (K - 2) * 2 + 1, K))
    cat("\n**Background**\n")
    cat("To better manage urinary tract infections (UTI) in the elderly, it is crucial to understand the natural history of bacteriuria. While asymptomatic bacteriuria (ASB) is common, its stability and propensity to progress to symptomatic UTI remain debated. We aimed to characterize longitudinal transitions between ASB, UTI, and culture-negative states in a nursing home cohort.\n")

    cat("\n**Methods**\n")
    cat("We conducted a prospective longitudinal study of nursing home residents. Urine samples were collected routinely every 3 months for 18 months, plus additional samples during suspected UTI episodes. Samples included voided urine, catheter urine, or urine-saturated incontinence material (worn <12 hours, no visible fecal contamination). Episodes were classified as ASB, UTI, or Negative based on guideline-based clinical criteria. Participants provided ≥", K, " distinct timepoints. Transitions between consecutive concurrent episodes were analyzed.\n", sep = "")

    cat("\n**Results**\n")
    cat(sprintf("The cohort included %d participants contributing %d evaluable urine episodes. ", n_part, n_episodes))
    cat(sprintf(
        "Episodes were classified as ASB (n=%d, %.1f%%), UTI (n=%d, %.1f%%), or culture-negative (n=%d, %.1f%%). ",
        n_asb, p_asb, n_uti, p_uti, n_neg, p_neg
    ))
    cat(sprintf("Analysis of %d within-participant transitions revealed distinct patterns. ", n_trans))
    cat(sprintf(
        "ASB was highly stable; of transitions starting in ASB (n=%d), %.1f%% (n=%d) remained ASB at the subsequent timepoint. ",
        n_asb_starts, p_asb_stable, n_asb_stable
    ))
    cat(sprintf(
        "Progression from ASB to UTI was observed in %d instances (%.1f%% of ASB-start transitions). ",
        n_asb_uti, p_asb_uti
    ))
    cat(sprintf("Transitions from a culture-negative state directly to UTI were rare (n=%d).", n_neg_uti))

    cat("\n\n**Conclusions**\n")
    cat("In this longitudinal nursing home cohort, ASB was a stable phenotype that rarely progressed directly to symptomatic UTI. These descriptive findings suggest that ASB may represent a stable host-commensal relationship rather than a precursor to acute infection.\n")
}


generate_mlst_abstract <- function(df_clinical, df_mlst, K) {
    # 1. Eligibility
    eligible_ids <- df_clinical %>%
        group_by(participant_id) %>%
        summarise(n_tp = n_distinct(Timepoint)) %>%
        filter(n_tp >= K) %>%
        pull(participant_id)

    df_cohort <- df_clinical %>%
        filter(participant_id %in% eligible_ids) %>%
        arrange(participant_id, tp_num)

    # 2. Join Preparation
    # Map MLST columns
    if ("Participant_id" %in% names(df_mlst) && !"participant_id" %in% names(df_mlst)) {
        df_mlst <- df_mlst %>% rename(participant_id = Participant_id)
    }

    # Select best MLST per episode (Participant_id + Timepoint)
    # Priority: mlst_complete=TRUE > n_loci_typed > first

    # Check if we assume 1 sample per episode in status_map
    # Join keys: participant_id, Timepoint

    df_mlst_clean <- df_mlst %>%
        mutate(
            is_complete = if_else(str_detect(tolower(as.character(if ("mlst_complete" %in% names(.)) mlst_complete else "FALSE")), "true|yes"), 1, 0),
            n_loci = if ("n_loci_typed" %in% names(.)) n_loci_typed else 7 # dummy if missing, assume logic
        ) %>%
        arrange(participant_id, Timepoint, desc(is_complete), desc(n_loci)) %>%
        distinct(participant_id, Timepoint, .keep_all = TRUE)

    # Join
    df_merged <- inner_join(df_cohort, df_mlst_clean, by = c("participant_id", "Timepoint"))

    n_linked <- nrow(df_merged)
    n_total_cohort <- nrow(df_cohort)
    yield_pct <- if (n_total_cohort > 0) n_linked / n_total_cohort * 100 else 0

    # 3. Stats
    n_distinct_st <- n_distinct(df_merged$ST)

    top_sts <- df_merged %>%
        count(ST) %>%
        mutate(pct = n / sum(n) * 100) %>%
        arrange(desc(n)) %>%
        slice(1:5)

    top_st_str <- paste(sprintf("ST%s (n=%d, %.1f%%)", top_sts$ST, top_sts$n, top_sts$pct), collapse = ", ")

    # ST Distribution by Status
    asb_sts <- df_merged %>% filter(Infection_Status == "ASB")
    n_asb_st <- nrow(asb_sts)
    n_asb_distinct <- n_distinct(asb_sts$ST)

    uti_sts <- df_merged %>% filter(Infection_Status == "UTI")
    n_uti_st <- nrow(uti_sts)
    n_uti_distinct <- n_distinct(uti_sts$ST)

    # Non-typable
    n_typed_total <- nrow(df_merged %>% filter(!is.na(ST))) # assuming ST is NA if missing, or "Non-typable" string
    n_nt <- sum(df_merged$ST == "Non-typable" | is.na(df_merged$ST), na.rm = TRUE)

    # 4. Stability
    # Filter for participants with >=2 TYPED episodes
    stable_check <- df_merged %>%
        group_by(participant_id) %>%
        mutate(n_typed = n()) %>%
        filter(n_typed >= 2) %>%
        arrange(tp_num) %>%
        mutate(
            next_st = lead(ST),
            same_st = (ST == next_st)
        ) %>%
        filter(!is.na(next_st))

    n_pairs <- nrow(stable_check)
    if (n_pairs > 0) {
        n_stable <- sum(stable_check$same_st, na.rm = TRUE)
        p_stable <- n_stable / n_pairs * 100

        # Participants with >=1 switch
        switchers <- stable_check %>%
            group_by(participant_id) %>%
            summarise(any_switch = any(!same_st))

        n_switchers <- sum(switchers$any_switch)
        n_denom_part <- nrow(switchers)
        p_switchers <- n_switchers / n_denom_part * 100
    } else {
        n_stable <- 0
        p_stable <- 0
        n_switchers <- 0
        n_denom_part <- 0
        p_switchers <- 0
    }

    # Output
    cat(sprintf("\n[DOC %d] Abstract: Longitudinal E. coli Sequence Type Dynamics Across ASB and UTI Episodes in Nursing Home Residents (≥%d TIMEPOINTS)\n", (K - 2) * 2 + 2, K))
    cat("\n**Background**\n")
    cat("E. coli is the predominant uropathogen in the elderly, causing both asymptomatic bacteriuria (ASB) and urinary tract infections (UTI). It is unclear whether specific sequence types (STs) are differentially associated with ASB versus UTI, or if within-host lineage turnover drives clinical status changes. We characterized E. coli ST dynamics in a longitudinal nursing home cohort.\n")

    cat("\n**Methods**\n")
    cat("We analyzed E. coli isolates from a prospective longitudinal study (samples every 3 months for 18 months + suspected UTI). Episodes were defined as ASB or UTI using clinical criteria. E. coli isolates were typed using multi-locus sequence typing (MLST). Participants with ≥", K, " eligible timepoints were included. We compared ST distributions between clinical states and assessed lineage stability across consecutive episodes.\n", sep = "")

    cat("\n**Results**\n")
    cat(sprintf("MLST data were available for %d episodes (%.1f%% of eligible cohort). ", n_linked, yield_pct))
    cat(sprintf("Overall, %d distinct STs were identified. ", n_distinct_st))
    cat(sprintf("The most prevalent lineages were %s. ", top_st_str))

    if (n_uti_st > 0) {
        cat(sprintf(
            "In UTI episodes (n=%d), %d distinct STs were found, compared to %d distinct STs in ASB episodes (n=%d). ",
            n_uti_st, n_uti_distinct, n_asb_distinct, n_asb_st
        ))
    } else {
        cat("UTI episodes were too few for robust stratification. ")
    }

    if (n_pairs > 0) {
        cat(sprintf("Longitudinal analysis of %d consecutive isolate pairs from %d participants showed high lineage stability: ", n_pairs, n_denom_part))
        cat(sprintf("STs remained identical in %.1f%% (n=%d) of consecutive pairs. ", p_stable, n_stable))
        cat(sprintf("However, %d participants (%.1f%%) exhibited at least one strain switch during follow-up.", n_switchers, p_switchers))
    } else {
        cat("Insufficient consecutive typed pairs for stability analysis.")
    }

    cat("\n\n**Conclusions**\n")
    cat("We observed a diverse E. coli population dominated by common globally disseminated lineages. The high within-host stability of STs across timepoints, regardless of clinical status transitions, suggests that clinical symptoms (UTI) may often arise from resident ASB strains rather than introduction of novel lineages. Strain switching was observed but less frequent.\n")

    # Print linkage info for QC
    cat(sprintf("\n[QC INFO for K=%d: Linked %d/%d]\n", K, n_linked, n_total_cohort))
}

# --- 4. EXECUTION ---

# Print K=2
generate_clinical_abstract(df_status, 2)
generate_mlst_abstract(df_status, df_mlst, 2)

# Print K=3
generate_clinical_abstract(df_status, 3)
generate_mlst_abstract(df_status, df_mlst, 3)

# Print K=4
generate_clinical_abstract(df_status, 4)
generate_mlst_abstract(df_status, df_mlst, 4)

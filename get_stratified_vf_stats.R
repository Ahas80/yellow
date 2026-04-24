suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(tidyr))

# Script: get_stratified_vf_stats.R
# Goal: Compare VF loads between ASB and UTI across >=2, >=3, and >=4 TP cohorts
# Filter: cfu_recorded_any == TRUE (Same as verified stats)

# 1. Load Data
# Clinical
status_df <- read_csv("status_map.csv", show_col_types = FALSE) %>%
    mutate(Participant_id = as.character(Participant_id))

# Linkage (Episode -> Isolate)
# DEDUPLICATION LOGIC: One isolate per episode (Participant + Timepoint)
# Priority: Has valid ST > Has Isolate ID > First available
filt_mlst <- read_csv("results/mlst/mlst_with_meta.csv", show_col_types = FALSE) %>%
    select(Isolate_ID, ST)

raw_inputs <- read_csv("results/class_inputs_full.csv", show_col_types = FALSE) %>%
    mutate(Participant_id = as.character(Participant_id)) %>%
    select(Participant_id, tp_lab, isolate_ID)

# Join to check ST status for ranking
class_inputs <- raw_inputs %>%
    left_join(filt_mlst, by = c("isolate_ID" = "Isolate_ID")) %>%
    mutate(
        has_st = !is.na(ST) & ST != "-" & ST != "NA",
        rank = case_when(
            has_st ~ 1,
            !is.na(isolate_ID) ~ 2,
            TRUE ~ 3
        )
    ) %>%
    group_by(Participant_id, tp_lab) %>%
    arrange(rank, isolate_ID) %>%
    slice(1) %>%
    ungroup() %>%
    select(Participant_id, tp_lab, isolate_ID) # Keep just the best linkage

# VF Data
vf_df <- read_csv("results/annotated_gene_table.csv", show_col_types = FALSE) %>%
    filter(Category != "Unassigned") %>% # Only count actual VFs
    distinct(Isolate_ID, Gene) # Deduplicate

# MLST Data (Reload for final join)
mlst_df <- read_csv("results/mlst/mlst_with_meta.csv", show_col_types = FALSE) %>%
    mutate(Participant_id = as.character(Participant_id)) %>%
    select(Isolate_ID, ST) %>%
    distinct(Isolate_ID, .keep_all = TRUE)

# 2. Clean & Filter Clinical Data
# LOGIC MATCHING VERIFIED STATS:
# 1. Determine eligible participants by counting n_distinct(Timepoint) (BEFORE filtering for CFU)
# 2. Restrict to episodes with cfu_recorded_any == TRUE
clean_df <- status_df %>%
    filter(!is.na(Infection_Status)) %>%
    group_by(Participant_id) %>%
    mutate(n_tp = n_distinct(Timepoint)) %>% # Count ALL timepoints
    ungroup() %>%
    filter(cfu_recorded_any == TRUE) # Then filter for analysis

# 3. Calculate VF Counts per Isolate
vf_counts <- vf_df %>%
    group_by(Isolate_ID) %>%
    summarize(n_vf = n())

# 4. Function: Calculate Stats per Stratum and RETURN them
get_stratum_metrics <- function(df, min_tps) {
    # Filter Cohort (Based on total retention)
    cohort <- df %>% filter(n_tp >= min_tps)

    # Link to Isolates (Now 1:1)
    cohort_isolates <- cohort %>%
        left_join(class_inputs, by = c("Participant_id", "Timepoint" = "tp_lab")) %>%
        filter(!is.na(isolate_ID))

    # Link to VF Counts & MLST
    cohort_vf <- cohort_isolates %>%
        left_join(vf_counts, by = c("isolate_ID" = "Isolate_ID")) %>%
        left_join(mlst_df, by = c("isolate_ID" = "Isolate_ID")) %>%
        mutate(n_vf = replace_na(n_vf, 0)) %>%
        filter(Infection_Status %in% c("ASB", "UTI"))

    # Calc Stats
    stats <- cohort_vf %>%
        group_by(Infection_Status) %>%
        summarize(
            N = n(),
            Mean = mean(n_vf),
            SD = sd(n_vf),
            Median = median(n_vf),
            Min = min(n_vf),
            Max = max(n_vf),
            .groups = "drop"
        )

    # Get values for Tests
    asb_vals <- cohort_vf$n_vf[cohort_vf$Infection_Status == "ASB"]
    uti_vals <- cohort_vf$n_vf[cohort_vf$Infection_Status == "UTI"]

    p_wilcox <- NA
    p_ttest <- NA

    if (length(asb_vals) > 0 && length(uti_vals) > 0) {
        p_wilcox <- wilcox.test(asb_vals, uti_vals)$p.value
        p_ttest <- t.test(asb_vals, uti_vals)$p.value
    }

    # Helper: Get Top ST String
    get_top_st <- function(sub_df, n_total) {
        if (nrow(sub_df) == 0) {
            return("None")
        }

        top <- sub_df %>%
            mutate(ST = case_when(
                is.na(ST) ~ "Not Sequenced",
                ST == "-" ~ "Unassigned ST",
                ST == "NA" ~ "Not Sequenced",
                TRUE ~ ST
            )) %>%
            count(ST) %>%
            mutate(prop = n / sum(n) * 100) %>%
            arrange(desc(n)) %>%
            head(3) %>%
            mutate(label = sprintf("%s (n=%d, %.1f%%)", ST, n, prop))

        return(paste(top$label, collapse = "; "))
    }

    # Format Helper
    fmt <- function(x, d = 1) sprintf(paste0("%.", d, "f"), x)

    # Extract ASB
    asb <- stats %>% filter(Infection_Status == "ASB")
    if (nrow(asb) == 0) asb <- list(N = 0, Mean = 0, SD = 0, Median = 0, Min = 0, Max = 0)

    # Extract UTI
    uti <- stats %>% filter(Infection_Status == "UTI")
    if (nrow(uti) == 0) uti <- list(N = 0, Mean = 0, SD = 0, Median = 0, Min = 0, Max = 0)

    # Extract Top STs
    asb_st_str <- get_top_st(cohort_vf %>% filter(Infection_Status == "ASB"), asb$N)
    uti_st_str <- get_top_st(cohort_vf %>% filter(Infection_Status == "UTI"), uti$N)

    return(list(
        ASB_N = as.character(asb$N),
        ASB_Mean_SD = paste0(fmt(asb$Mean), " (", fmt(asb$SD), ")"),
        ASB_Median_Range = paste0(asb$Median, " (", asb$Min, "-", asb$Max, ")"),
        ASB_TopST = asb_st_str,
        UTI_N = as.character(uti$N),
        UTI_Mean_SD = paste0(fmt(uti$Mean), " (", fmt(uti$SD), ")"),
        UTI_Median_Range = paste0(uti$Median, " (", uti$Min, "-", uti$Max, ")"),
        UTI_TopST = uti_st_str,
        P_Wilcox = if (!is.na(p_wilcox)) sprintf("%.4f", p_wilcox) else "NA",
        P_Ttest = if (!is.na(p_ttest)) sprintf("%.4f", p_ttest) else "NA"
    ))
}

# Run Analysis for all Strata
res2 <- get_stratum_metrics(clean_df, 2)
res3 <- get_stratum_metrics(clean_df, 3)
res4 <- get_stratum_metrics(clean_df, 4)

# Build Side-by-Side Table
final_table <- data.frame(
    Metric = c(
        "ASB Episodes (N)",
        "ASB VF Score: Mean (SD)",
        "ASB VF Score: Median (Range)",
        "Top ASB STs (Ranked)",
        "UTI Episodes (N)",
        "UTI VF Score: Mean (SD)",
        "UTI VF Score: Median (Range)",
        "Top UTI STs (Ranked)",
        "Comparison: Wilcoxon p-value",
        "Comparison: T-test p-value"
    ),
    `TP_ge_2` = c(
        res2$ASB_N, res2$ASB_Mean_SD, res2$ASB_Median_Range, res2$ASB_TopST,
        res2$UTI_N, res2$UTI_Mean_SD, res2$UTI_Median_Range, res2$UTI_TopST,
        res2$P_Wilcox, res2$P_Ttest
    ),
    `TP_ge_3` = c(
        res3$ASB_N, res3$ASB_Mean_SD, res3$ASB_Median_Range, res3$ASB_TopST,
        res3$UTI_N, res3$UTI_Mean_SD, res3$UTI_Median_Range, res3$UTI_TopST,
        res3$P_Wilcox, res3$P_Ttest
    ),
    `TP_ge_4` = c(
        res4$ASB_N, res4$ASB_Mean_SD, res4$ASB_Median_Range, res4$ASB_TopST,
        res4$UTI_N, res4$UTI_Mean_SD, res4$UTI_Median_Range, res4$UTI_TopST,
        res4$P_Wilcox, res4$P_Ttest
    )
)

# Rename Columns for Display
colnames(final_table) <- c("Metric", ">=2 Timepoints", ">=3 Timepoints", ">=4 Timepoints")

cat("\n--- Side-by-Side VF Statistics ---\n")
print(final_table, row.names = FALSE, right = FALSE)

# Export
write_csv(final_table, "results/stratified_vf_stats_table.csv")
cat("\nTable saved to results/stratified_vf_stats_table.csv\n")

# Script: get_extended_stats.R
# Goal: Calculate additional stats (Sample Type, VF Load) for the verified Abstract cohort

suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(tidyr))

# --- CONFIG ---
MIN_COMPLETE_EPISODES <- 1
class_inputs_path <- "results/class_inputs_full.csv"
mlst_path <- "results/mlst/mlst_with_meta.csv"
vf_path <- "results/annotated_gene_table.csv"

# --- 1. LOAD & RECONSTRUCT N=87 COHORT (Same logic as V3) ---
message("Loading Data...")
class_inputs <- read_csv(class_inputs_path, show_col_types = FALSE) %>%
    mutate(Participant_id = as.character(Participant_id), tp_num = as.character(tp_num))

mlst_df <- read_csv(mlst_path, show_col_types = FALSE) %>%
    mutate(Participant_id = as.character(Participant_id))

vf_df <- read_csv(vf_path, show_col_types = FALSE) %>%
    mutate(Participant_id = as.character(Participant_id))

# Rebuild Master DF
base_df <- class_inputs %>%
    rename(pid = Participant_id, isolate_id = isolate_ID, status_clinical = Infection_Status) %>%
    distinct(pid, tp_num, isolate_id, .keep_all = TRUE)

gen_data <- mlst_df %>%
    rename(isolate_id = Isolate_ID) %>%
    select(isolate_id, ST) %>%
    distinct(isolate_id, .keep_all = TRUE)

master_df <- base_df %>%
    left_join(gen_data, by = "isolate_id") %>%
    mutate(
        has_clinical = TRUE,
        has_isolate = !is.na(isolate_id),
        has_genome = !is.na(ST),
        is_complete_episode = has_clinical & has_isolate & has_genome
    )

# Filter
pid_stats <- master_df %>%
    group_by(pid) %>%
    summarize(n_complete = n_distinct(paste(tp_num, isolate_id)[is_complete_episode]))

valid_pids <- pid_stats %>%
    filter(n_complete >= MIN_COMPLETE_EPISODES) %>%
    pull(pid)

final_df <- master_df %>%
    filter(pid %in% valid_pids) %>%
    filter(is_complete_episode)

cat("Verified Cohort N:", n_distinct(final_df$pid), "\n")

# --- 2. SAMPLE TYPE STATS ---
message("\n--- SAMPLE TYPES ---")
# 'Urine collection method'
# Clean up labels if needed (e.g. Translate Dutch)
final_df <- final_df %>%
    mutate(SampleType = case_when(
        grepl("Spontaan", `Urine collection method`, ignore.case = T) ~ "Voided",
        grepl("Inco", `Urine collection method`, ignore.case = T) ~ "Incontinence Material",
        grepl("Katheter", `Urine collection method`, ignore.case = T) ~ "Catheter",
        TRUE ~ "Other/Unknown"
    ))

sample_stats <- final_df %>%
    distinct(pid, tp_num, status_clinical, SampleType) %>% # distinct episodes
    count(status_clinical, SampleType) %>%
    group_by(status_clinical) %>%
    mutate(pct = n / sum(n) * 100)

print(sample_stats)

# --- 3. VIRULENCE FACTOR LOAD ---
message("\n--- VIRULENCE FACTOR LOAD ---")

# Count Genes per Isolate in VF outcomes
# Using Category != "Unassigned" as proxy for VF
vf_counts <- vf_df %>%
    filter(Category != "Unassigned") %>%
    group_by(Isolate_ID) %>%
    summarize(n_vf = n_distinct(Gene))

# Join to final DF
vf_analysis <- final_df %>%
    left_join(vf_counts, by = c("isolate_id" = "Isolate_ID")) %>%
    mutate(n_vf = replace_na(n_vf, 0)) # 0 if no match

vf_summary <- vf_analysis %>%
    group_by(status_clinical) %>%
    summarize(
        Mean_VF = mean(n_vf),
        Median_VF = median(n_vf),
        SD_VF = sd(n_vf),
        N_Isolates = n()
    )

print(vf_summary)

# T-test if possible
tryCatch(
    {
        asb_vfs <- vf_analysis$n_vf[vf_analysis$status_clinical == "ASB"]
        uti_vfs <- vf_analysis$n_vf[vf_analysis$status_clinical == "UTI"]
        test <- t.test(asb_vfs, uti_vfs)
        cat("\nT-test ASB vs UTI VF Count: p-value =", test$p.value, "\n")
    },
    error = function(e) {
        cat("Test failed\n")
    }
)

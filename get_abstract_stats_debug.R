suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(readr))

# Load Data
status_df <- read_csv("status_map.csv", show_col_types = FALSE)
ci <- read_csv("results/class_inputs_full.csv", show_col_types = F) %>%
    select(Participant_id, tp_lab, isolate_ID) %>%
    mutate(Participant_id = as.character(Participant_id)) %>%
    distinct(Participant_id, tp_lab, .keep_all = TRUE)
mlst <- read_csv("results/mlst/mlst_with_meta.csv", show_col_types = F) %>%
    select(Isolate_ID, ST) %>%
    distinct(Isolate_ID, .keep_all = TRUE)

# Base Cleaning
df_base <- status_df %>%
    filter(!is.na(Infection_Status)) %>%
    group_by(Participant_id) %>%
    mutate(n_total_tp = n_distinct(Timepoint)) %>%
    ungroup() %>%
    filter(cfu_recorded_any == TRUE)

cat("\n--- Base (>=2 Total TPs) ---\n")
base2 <- df_base %>% filter(n_total_tp >= 2)
cat("Participants:", n_distinct(base2$Participant_id), "\n") # Exp: 92

# Hypothesis F: >= 2 *Genotyped* Episodes (Valid ST)
df_st <- df_base %>%
    left_join(ci, by = c("Participant_id", "Timepoint" = "tp_lab")) %>%
    left_join(mlst, by = c("isolate_ID" = "Isolate_ID")) %>%
    mutate(has_st = !is.na(ST) & ST != "-" & ST != "NA")

h_st <- df_st %>%
    filter(has_st == TRUE) %>%
    group_by(Participant_id) %>%
    filter(n() >= 2) %>% # At least 2 genotyped episodes
    ungroup()

cat("\n--- Hypothesis F: >=2 Genotyped Episodes (Valid ST) ---\n")
n_h_st <- n_distinct(h_st$Participant_id)
cat("Participants:", n_h_st, "\n")
cat("Genotyped Episodes:", nrow(h_st), "\n")

# If this is the cohort, get their TOTAL stats (including non-sequenced eps)
if (n_h_st > 0) {
    h_st_ids <- unique(h_st$Participant_id)
    h_st_all <- df_base %>% filter(Participant_id %in% h_st_ids)
    cat("Total Episodes for Genotyped Cohort:", nrow(h_st_all), "\n")
    print(table(h_st_all$Infection_Status))
} else {
    cat("Hypothesis F yielded 0 participants.\n")
}

# --- EXTENDED: CALCULATE >=4 TP Stats for Hypothesis F ---
# (Assuming Hypothesis F is close enough to 71)
if (n_h_st > 0) {
    cat("\n--- >= 4 TPs for Genotyped Cohort Logic ---\n")
    # Logic: Filter for people with >=4 Genotyped Episodes? Or >=4 Total?
    # User Abstract for >=3: "54 participants" (Matches my Base >=3).
    # Wait, my Base >=3 was 54. User was 54.
    # But User >=2 was 71. My Base >=2 was 92.

    # This implies the filter is applied ONLY to the >=2 group? Or inconsistent?
    # Let's check >= 3 Genotyped Episodes
    h_st_3 <- df_st %>%
        filter(has_st == TRUE) %>%
        group_by(Participant_id) %>%
        filter(n() >= 3) %>%
        ungroup()
    cat("Participants with >=3 Genotyped Eps:", n_distinct(h_st_3$Participant_id), "\n")
}

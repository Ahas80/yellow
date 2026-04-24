suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(readr))

# Load Data
status_df <- read_csv("status_map.csv", show_col_types = FALSE)

# Base Cleaning (My previous method)
df_base <- status_df %>%
    group_by(Participant_id) %>%
    mutate(n_total_tp = n_distinct(Timepoint)) %>%
    ungroup() %>%
    filter(cfu_recorded_any == TRUE) # Episodes with CFU data

# Hypothesis 1: My Base (>=2 Total Timepoints)
h1 <- df_base %>% filter(n_total_tp >= 2)
cat("\n--- Hypothesis 1: Base (n_total_tp >= 2) ---\n")
cat("Participants:", n_distinct(h1$Participant_id), "\n")
cat("Episodes:", nrow(h1), "\n")

# Hypothesis 2: >= 2 *Evaluated* Timepoints (with CFU)
# (i.e., not just retention, but actual data points)
h2 <- df_base %>%
    group_by(Participant_id) %>%
    filter(n() >= 2) %>%
    ungroup()
cat("\n--- Hypothesis 2: >=2 Episodes with CFU ---\n")
cat("Participants:", n_distinct(h2$Participant_id), "\n")
cat("Episodes:", nrow(h2), "\n")

# Hypothesis 3: Exclude Participants with ONLY Negative results
# (i.e. must have at least one ASB or UTI)
h3 <- h1 %>%
    group_by(Participant_id) %>%
    filter(any(Infection_Status %in% c("ASB", "UTI"))) %>%
    ungroup()
cat("\n--- Hypothesis 3: Base + At least 1 Positive (ASB/UTI) ---\n")
cat("Participants:", n_distinct(h3$Participant_id), "\n")
cat("Episodes:", nrow(h3), "\n")

# Hypothesis 4: Exclude Participants with NO Sequencing Data?
# (Requires joining with mlst or class_inputs)
# Load linkage
ci <- read_csv("results/class_inputs_full.csv", show_col_types = F) %>%
    select(Participant_id, tp_lab, isolate_ID) %>%
    mutate(Participant_id = as.character(Participant_id))

h4 <- h1 %>%
    left_join(ci, by = c("Participant_id", "Timepoint" = "tp_lab")) %>%
    group_by(Participant_id) %>%
    filter(any(!is.na(isolate_ID))) %>%
    ungroup() %>%
    distinct(Participant_id, Timepoint, .keep_all = T)

cat("\n--- Hypothesis 4: Base + At least 1 Isolate ---\n")
cat("Participants:", n_distinct(h4$Participant_id), "\n")
cat("Episodes:", nrow(h4), "\n")

# Hypothesis 5: Check the 'evaluable urine episodes' phrasing
# User text: 236 evaluable episodes.
# My Base h1 has 258. Difference = 22.
# H2 has 256.
# Let's check status dist of h2:
cat("\nStatus Distribution of H2 (>=2 CFU episodes):\n")
print(table(h2$Infection_Status))

# Hypothesis 6: Maybe user's text comes from `get_abstract_stats_v3.R` logic?
# I'll check if there's a strict filter like "culture_pos_epi"
h6 <- h1 %>% filter(culture_positive_epi == TRUE) # if column exists?
if ("culture_pos_epi" %in% names(h1)) {
    h6 <- h1 %>% filter(culture_pos_epi == TRUE)
    cat("\n--- Hypothesis 6: Culture Positive Only ---\n")
    cat("Participants:", n_distinct(h6$Participant_id), "\n")
    cat("Episodes:", nrow(h6), "\n")
}

# Check counts specifically for >=3 TPs match
# User: 54 Pts, 185 Eps.
# My Base >=3: 54 Pts, 187 Eps. (2 extra UTI).
base3 <- df_base %>% filter(n_total_tp >= 3)
cat("\n--- Base >= 3 TPs ---\n")
cat("Participants:", n_distinct(base3$Participant_id), "\n")
cat("Episodes:", nrow(base3), "\n")
print(table(base3$Infection_Status))

library(dplyr)
library(readr)

# Load Data
message("Loading Data...")
class_inputs <- read_csv("results/class_inputs_full.csv", show_col_types = FALSE)
mlst <- read_csv("results/mlst/mlst_with_meta.csv", show_col_types = FALSE)
clinical <- read_csv("results/clinical/status_map.csv", show_col_types = FALSE) %>%
    filter(cfu_recorded_any == TRUE)

# Re-establish the cohort (>=2 TPs)
cohort_2tp <- clinical %>%
    group_by(Participant_id) %>%
    filter(n_distinct(Timepoint) >= 2) %>%
    ungroup() %>%
    mutate(Participant_id = as.character(Participant_id))

# Link to isolate IDs
clinical_with_isolate <- cohort_2tp %>%
    left_join(class_inputs, by = c("Participant_id", "Timepoint" = "tp_lab"))

# Identify the 106 missing isolates
# (Episodes -> Linked IsolateID -> Not in MLST file)
missing_mlst <- clinical_with_isolate %>%
    filter(!is.na(isolate_ID)) %>%
    anti_join(mlst, by = c("isolate_ID" = "Isolate_ID"))

message(paste("Count of episodes with Isolate ID but missing MLST:", nrow(missing_mlst)))

# Analyze these missing isolates using class_inputs metadata
# class_inputs might have organism info
message("\n--- Analyzing Missing Isolates using Class Inputs ---")
print(colnames(class_inputs))



# Link back to class_inputs to get Organism
missing_with_org <- missing_mlst %>%
    left_join(class_inputs %>% select(isolate_ID, Organism), by = "isolate_ID")

message("\n--- Organism Distribution of Missing MLST Isolates ---")
print(table(missing_with_org$Organism.y, useNA = "ifany"))

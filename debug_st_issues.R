suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(readr))

# Load Data
status_df <- read_csv("status_map.csv", show_col_types = FALSE)
class_inputs <- read_csv("results/class_inputs_full.csv", show_col_types = FALSE)
mlst_df <- read_csv("results/mlst/mlst_with_meta.csv", show_col_types = FALSE)

# 1. Investigate Many-to-Many Join
# Check duplicates in class_inputs for (Participant_id, tp_lab)
dups_ci <- class_inputs %>%
    mutate(Participant_id = as.character(Participant_id)) %>%
    group_by(Participant_id, tp_lab) %>%
    summarize(n = n(), isolates = paste(unique(isolate_ID), collapse = ";")) %>%
    filter(n > 1)

cat("\n--- Duplicates in Class Inputs (P_ID + TP_Lab) ---\n")
print(head(dups_ci, 10))

# 2. Investigate Missing STs
# Recreate the cohort logic briefly
clean_df <- status_df %>%
    filter(cfu_recorded_any == TRUE, !is.na(Infection_Status)) %>%
    mutate(Participant_id = as.character(Participant_id))

# Join
joined <- clean_df %>%
    filter(n_distinct(Timepoint) >= 2) %>% # Just check >=2 cohort
    left_join(class_inputs %>% mutate(Participant_id = as.character(Participant_id)),
        by = c("Participant_id", "Timepoint" = "tp_lab")
    )

# Look at rows with NA or - ST
# Link to MLST
joined_mlst <- joined %>%
    left_join(mlst_df %>% select(Isolate_ID, ST), by = c("isolate_ID" = "Isolate_ID"))

missing_st <- joined_mlst %>%
    filter(is.na(ST) | ST == "-" | ST == "NA") %>%
    select(Participant_id, Timepoint, isolate_ID, ST, Infection_Status)

cat("\n--- Rows with Missing/Unknown ST (First 20) ---\n")
print(head(missing_st, 20))

# Check if these 'missing' isolate_IDs actually exist in MLST file
isolated_ids_missing <- unique(missing_st$isolate_ID)
present_in_mlst <- mlst_df %>% filter(Isolate_ID %in% isolated_ids_missing)

cat("\n--- Do these 'Missing' Isolate IDs exist in the MLST file? ---\n")
cat("Missing IDs count:", length(isolated_ids_missing), "\n")
cat("Found in MLST file:", nrow(present_in_mlst), "\n")
if (nrow(present_in_mlst) > 0) print(head(present_in_mlst))

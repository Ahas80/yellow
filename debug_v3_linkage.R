# Script: debug_v3_linkage.R
library(dplyr)
library(readr)

status_map_path <- "results/status_map_full.csv"
class_inputs_path <- "results/class_inputs_full.csv"

# Load and standardize
status_map <- read_csv(status_map_path, show_col_types = FALSE) %>%
    mutate(Participant_id = as.character(Participant_id), Timepoint = as.character(Timepoint))

class_inputs <- read_csv(class_inputs_path, show_col_types = FALSE) %>%
    mutate(Participant_id = as.character(Participant_id), tp_num = as.character(tp_num))

# Check keys before join
cat("--- Clinical Base (Status Map) ---\n")
print(head(status_map %>% select(Participant_id, Timepoint, Infection_Status)))

cat("\n--- Sequence Data (Class Inputs) ---\n")
print(head(class_inputs %>% select(Participant_id, tp_num, isolate_ID, Infection_Status)))

# Check intersection of Keys (PID + Timepoint)
clinical_keys <- paste(status_map$Participant_id, status_map$Timepoint)
seq_keys <- paste(class_inputs$Participant_id, class_inputs$tp_num)

cat("\nIntersection of Episode Keys (PID + TP):\n")
cat("Clinical N:", length(clinical_keys), "\n")
cat("Seq N:", length(seq_keys), "\n")
cat("Intersect N:", length(intersect(clinical_keys, seq_keys)), "\n")

# Check why Isolates are NA in join
# Try joining specifically for a known case: 120003 T0
cat("\n--- Specific Check: 120003 T0 ---\n")
c_sub <- status_map %>% filter(Participant_id == "120003", Timepoint == "0")
s_sub <- class_inputs %>% filter(Participant_id == "120003", tp_num == "0")

print(c_sub)
print(s_sub)

joined_sub <- c_sub %>%
    left_join(s_sub, by = c("Participant_id" = "Participant_id", "Timepoint" = "tp_num"))
print(joined_sub)

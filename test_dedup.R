suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(readr))

# Load
class_inputs <- read_csv("results/class_inputs_full.csv", show_col_types = FALSE) %>%
    mutate(Participant_id = as.character(Participant_id)) %>%
    select(Participant_id, tp_lab, isolate_ID)

mlst_df <- read_csv("results/mlst/mlst_with_meta.csv", show_col_types = FALSE) %>%
    select(Isolate_ID, ST)

# Check duplicates before
dups_before <- class_inputs %>%
    group_by(Participant_id, tp_lab) %>%
    filter(n() > 1) %>%
    nrow()
cat("Duplicate rows in class_inputs before:", dups_before, "\n")

# Logic: Prioritize Isolates
# 1. Join with MLST to see which have STs
# 2. Arrange by: Has ST, then Isolate ID
# 3. Slice Head
class_inputs_dedup <- class_inputs %>%
    left_join(mlst_df, by = c("isolate_ID" = "Isolate_ID")) %>%
    mutate(
        has_st = !is.na(ST) & ST != "-" & ST != "NA",
        # Create a rank for sorting: 1=Best
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
    select(Participant_id, tp_lab, isolate_ID, ST) # Keep determined ST

# Check duplicates after
dups_after <- class_inputs_dedup %>%
    group_by(Participant_id, tp_lab) %>%
    filter(n() > 1) %>%
    nrow()
cat("Duplicate rows in class_inputs after:", dups_after, "\n")

# Check specific problem case 20002 T0
cat("\nCheck 20002 T0:\n")
print(class_inputs_dedup %>% filter(Participant_id == "20002", tp_lab == "T0"))

# Check NA ST count in deduped set
cat("\nIsolates with NA/Unknown ST in deduped set:\n")
print(table(is.na(class_inputs_dedup$ST) | class_inputs_dedup$ST == "-"))

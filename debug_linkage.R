library(dplyr)
library(readr)

status_df <- read_csv("status_map.csv", show_col_types = FALSE)
class_inputs <- read_csv("results/class_inputs_full.csv", show_col_types = FALSE)

# Check UTI episodes in status_map
uti_eps <- status_df %>%
    filter(Infection_Status == "UTI", cfu_recorded_any == TRUE) %>%
    mutate(tp_clean = gsub("T", "", Timepoint))

cat("UTI Episodes in Status Map (CFU+):\n")
print(uti_eps %>% select(Participant_id, Timepoint, tp_clean) %>% head())

# Check class_inputs for these PIDs
cat("\nClass Inputs for these PIDs:\n")
pids <- unique(uti_eps$Participant_id)
ci_sub <- class_inputs %>%
    filter(Participant_id %in% pids) %>%
    select(Participant_id, tp_num, isolate_ID)
print(ci_sub)

# Try join
joined <- uti_eps %>%
    mutate(Participant_id = as.character(Participant_id)) %>%
    left_join(ci_sub %>% mutate(Participant_id = as.character(Participant_id), tp_num = as.character(tp_num)),
        by = c("Participant_id", "tp_clean" = "tp_num")
    )

cat("\nJoined UTIs:\n")
print(joined %>% select(Participant_id, Timepoint, isolate_ID))

cat("\nNon-NA Isolates in Joined:\n")
print(sum(!is.na(joined$isolate_ID)))

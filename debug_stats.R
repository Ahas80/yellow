suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(readr))

status_path <- "results/class_inputs_full.csv"
mlst_path <- "results/mlst/mlst_with_meta.csv"

status_df <- read_csv(status_path, show_col_types = FALSE) %>% rename(isolate_id = isolate_ID, pid = Participant_id)
mlst_df <- read_csv(mlst_path, show_col_types = FALSE) %>% rename(isolate_id = Isolate_ID)

complete_df <- status_df %>%
    filter(Infection_Status %in% c("ASB", "UTI")) %>%
    filter(!is.na(isolate_id)) %>%
    inner_join(mlst_df %>% select(isolate_id, ST), by = "isolate_id")

# Check overlap
asb_pids <- unique(complete_df$pid[complete_df$Infection_Status == "ASB"])
uti_pids <- unique(complete_df$pid[complete_df$Infection_Status == "UTI"])

cat("ASB PIDs:", length(asb_pids), "\n")
cat("UTI PIDs:", length(uti_pids), "\n")
intersect_pids <- intersect(asb_pids, uti_pids)
cat("Intersect PIDs:", length(intersect_pids), "\n")
print(intersect_pids)

if (length(intersect_pids) == 0) {
    cat("\nChecking if PIDs are formatted differently?\n")
    print(head(asb_pids))
    print(head(uti_pids))
}

suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(readr))

# Paths
status_path <- "results/class_inputs_full.csv"
mlst_path <- "results/mlst/mlst_with_meta.csv"

# Load
message("Loading raw files...")
status_df <- read_csv(status_path, show_col_types = FALSE)
mlst_df <- read_csv(mlst_path, show_col_types = FALSE) %>% rename(isolate_id = Isolate_ID)

cat("\n--- 1. RAW Clinical File (class_inputs_full.csv) ---\n")
cat("Total Rows:", nrow(status_df), "\n")
cat("Unique Participant_id:", n_distinct(status_df$Participant_id), "\n")
cat("Unique Participant + Timepoint (Episodes?):", nrow(status_df %>% distinct(Participant_id, tp_num)), "\n")

# Check duplicates by status
cat("\nRaw Row Counts by Status:\n")
print(table(status_df$Infection_Status))

# Check for duplicates of (Participant + Timepoint)
# Are there multiple rows for the same episode?
dupes <- status_df %>%
    group_by(Participant_id, tp_num, Infection_Status) %>%
    count() %>%
    filter(n > 1)

cat("\nEpisodes with multiple rows (duplicates?):", nrow(dupes), "\n")
if (nrow(dupes) > 0) {
    cat("Example duplicate episode:\n")
    ex_pid <- dupes$Participant_id[1]
    ex_tp <- dupes$tp_num[1]
    print(status_df %>% filter(Participant_id == ex_pid, tp_num == ex_tp))
}

cat("\n--- 2. Filter: Has Isolate ID ---\n")
status_iso <- status_df %>% filter(!is.na(isolate_ID))
cat("Rows with Isolate ID:", nrow(status_iso), "\n")
cat("Unique Isolate IDs:", n_distinct(status_iso$isolate_ID), "\n")

# Duplicate check again
dupe_iso <- status_iso %>%
    group_by(isolate_ID) %>%
    count() %>%
    filter(n > 1)
cat("Isolate IDs appearing in multiple rows:", nrow(dupe_iso), "\n")


cat("\n--- 3. Filter: Join with MLST ---\n")
# We want to know how many DISTINCT EPISODES have valid MLST.
# Unique Key for Episode: Participant_id + tp_num (or tp_lab)
# Unique Key for Isolate: isolate_ID

# Left join status to MLST to see what matches
joined <- status_iso %>%
    left_join(mlst_df %>% select(isolate_id, ST, scheme), by = c("isolate_ID" = "isolate_id"))

cat("Status Rows after Join:", nrow(joined), "\n")
joined_complete <- joined %>% filter(!is.na(ST))
cat("Rows with valid MLST match:", nrow(joined_complete), "\n")

# Now calculate the STATISTICS for the abstract
# 1. Unique Participants
n_part <- n_distinct(joined_complete$Participant_id)
# 2. Unique Episodes (ASB)
# Filter for ASB first
asb_df <- joined_complete %>% filter(Infection_Status == "ASB")
n_asb_episodes <- n_distinct(paste(asb_df$Participant_id, asb_df$tp_num))

# 3. Unique Episodes (UTI)
uti_df <- joined_complete %>% filter(Infection_Status == "UTI")
n_uti_episodes <- n_distinct(paste(uti_df$Participant_id, uti_df$tp_num))

cat("\n--- FINAL AUDITED COUNTS (Unique Episodes) ---\n")
cat("Participants:", n_part, "\n")
cat("ASB Episodes:", n_asb_episodes, "\n")
cat("UTI Episodes:", n_uti_episodes, "\n")

cat("\n--- CHECK: Why was 20034 UTI dropped? ---\n")
# 20034 UTI isolates from before: 24240065001-1, 24360035101-1
check_ids <- c("24240065001-1", "24360035101-1")
cat("Checking specific IDs in MLST file:\n")
for (id in check_ids) {
    present <- id %in% mlst_df$isolate_id
    cat(id, "Present in MLST?", present, "\n")
}

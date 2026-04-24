suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(readr))

# Load Clinical
status_df <- read_csv("results/class_inputs_full.csv", show_col_types = FALSE)

# 1. Raw intersection
pids_asb_raw <- status_df %>%
    filter(Infection_Status == "ASB") %>%
    pull(Participant_id) %>%
    unique()
pids_uti_raw <- status_df %>%
    filter(Infection_Status == "UTI") %>%
    pull(Participant_id) %>%
    unique()
intersect_raw <- intersect(pids_asb_raw, pids_uti_raw)

cat("--- RAW CLINICAL DATA ---\n")
cat("ASB PIDs:", length(pids_asb_raw), "\n")
cat("UTI PIDs:", length(pids_uti_raw), "\n")
cat("Intersection (Any ASB + Any UTI):", length(intersect_raw), "\n")
if (length(intersect_raw) > 0) cat("PIDs:", paste(head(intersect_raw), collapse = ", "), "\n")

# 2. With Isolate ID
status_iso <- status_df %>% filter(!is.na(isolate_ID))
pids_asb_iso <- status_iso %>%
    filter(Infection_Status == "ASB") %>%
    pull(Participant_id) %>%
    unique()
pids_uti_iso <- status_iso %>%
    filter(Infection_Status == "UTI") %>%
    pull(Participant_id) %>%
    unique()
intersect_iso <- intersect(pids_asb_iso, pids_uti_iso)

cat("\n--- WITH ISOLATE ID ---\n")
cat("ASB PIDs:", length(pids_asb_iso), "\n")
cat("UTI PIDs:", length(pids_uti_iso), "\n")
cat("Intersection:", length(intersect_iso), "\n")
if (length(intersect_iso) > 0) cat("PIDs:", paste(head(intersect_iso), collapse = ", "), "\n")

# 3. With MLST Match
mlst_df <- read_csv("results/mlst/mlst_with_meta.csv", show_col_types = FALSE)
# Check IDs
# status_iso$isolate_ID vs mlst_df$Isolate_ID

# Direct check of IDs in intersection
if (length(intersect_iso) > 0) {
    target_pids <- intersect_iso

    # Get isolate IDs for these people
    target_iso_ids <- status_iso %>%
        filter(Participant_id %in% target_pids) %>%
        filter(Infection_Status %in% c("ASB", "UTI")) %>%
        pull(isolate_ID)

    # Check existence in MLST
    found_in_mlst <- target_iso_ids %in% mlst_df$Isolate_ID
    cat("\n--- MLST LOOKUP ---\n")
    cat("Target Isolates from Intersection participants:", length(target_iso_ids), "\n")
    cat("Found in MLST file:", sum(found_in_mlst), "\n")

    if (sum(found_in_mlst) < length(target_iso_ids)) {
        cat("Missing examples:\n")
        print(head(target_iso_ids[!found_in_mlst]))
        cat("Sample MLST IDs:\n")
        print(head(mlst_df$Isolate_ID))
    }
}

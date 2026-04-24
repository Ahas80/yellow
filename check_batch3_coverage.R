library(tidyverse)

# Load Data
df_status <- read_csv("results/clinical/status_map.csv", show_col_types = FALSE)
df_meta <- read_csv("assembly_metadata.csv", show_col_types = FALSE)

# Identify Clinical Batch 3 Participants
batch3_ids <- df_status %>%
    filter(grepl("3", as.character(Batch))) %>% # catching "3" and "1,3" etc
    pull(Participant_id) %>%
    unique()

cat(sprintf("Number of participants with Batch 3 clinical data: %d\n", length(batch3_ids)))
print(head(batch3_ids))

# Check coverage in Assembly Metadata
meta_participants <- unique(df_meta$Participant_id)

batch3_in_meta <- intersect(batch3_ids, meta_participants)
cat(sprintf("Number of Batch 3 participants found in Assembly Metadata: %d\n", length(batch3_in_meta)))

if (length(batch3_in_meta) < length(batch3_ids)) {
    cat("MISSING: Some Batch 3 participants have no assemblies!\n")
    missing_ids <- setdiff(batch3_ids, meta_participants)
    print(head(missing_ids))
} else {
    cat("SUCCESS: All Batch 3 participants have assemblies.\n")
}

# Check if the Batch column in metadata is just missing "3"
if ("Batch" %in% names(df_meta)) {
    batch3_meta_rows <- df_meta %>% filter(Participant_id %in% batch3_ids)
    print("Batch labels in metadata for these 'Batch 3' participants:")
    print(table(batch3_meta_rows$Batch, useNA = "ifany"))
}

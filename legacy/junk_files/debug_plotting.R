source("00_config.R")
suppressPackageStartupMessages({
    library(dplyr)
    library(tidyr)
    library(readr)
    library(stringr)
})

# Helper to safely read and cast types
read_batch_debug <- function(path) {
    df <- read_csv(path, show_col_types = FALSE)

    # Fix BOM
    nms <- names(df)
    nms[1] <- gsub("^[^A-Za-z0-9]+", "", nms[1])
    names(df) <- nms

    if (!"Timepoint" %in% names(df)) {
        tp_col <- names(df)[grepl("time.*point", names(df), ignore.case = TRUE)]
        if (length(tp_col) > 0) df <- df %>% rename(Timepoint = all_of(tp_col))
    }

    df %>%
        mutate(
            Participant_id = as.character(Participant_id),
            isolate_ID = as.character(isolate_ID),
            Timepoint = as.character(Timepoint)
        ) %>%
        select(Participant_id, isolate_ID, Timepoint)
}

b1 <- read_batch_debug("data/inputs/batch1.csv")
b2 <- read_batch_debug("data/inputs/batch2.csv")
b3 <- read_batch_debug("data/inputs/batch3.csv")

meta_map <- bind_rows(b1, b2, b3)

# Helper to normalize TP
tp_norm_val <- function(x) {
    x <- as.character(x)
    is_uricult <- str_detect(x, regex("uricult", ignore_case = TRUE))
    tp_num <- suppressWarnings(as.integer(str_extract(x, "\\d+")))
    case_when(
        is_uricult ~ "Uricult",
        !is.na(tp_num) ~ paste0("T", tp_num),
        TRUE ~ "Unscheduled"
    )
}

meta_map <- meta_map %>%
    mutate(Timepoint = tp_norm_val(Timepoint))

# Load MLST
mlst <- read_tsv(FILE_MLST_ALL, show_col_types = FALSE)

# Extract Isolate_ID from filename
mlst_parsed <- mlst %>%
    mutate(
        extracted_ID = str_extract(file_name, "24[0-9A-Za-z]+-[0-9]+")
    ) %>%
    # Drop original Participant_id and Timepoint from MLST to avoid conflict
    select(-any_of(c("Participant_id", "Timepoint")))

# Join
joined_full <- mlst_parsed %>%
    inner_join(meta_map, by = c("extracted_ID" = "isolate_ID"))

message("Joined Full rows: ", nrow(joined_full))
message("Joined Full columns: ", paste(names(joined_full), collapse = ", "))

# Join with status
status <- read_csv(file.path(DIR_RESULTS, "status_map.csv"), show_col_types = FALSE) %>%
    mutate(Participant_id = as.character(Participant_id)) %>%
    mutate(tp_lab = tp_norm_val(Timepoint))

final_df <- joined_full %>%
    inner_join(status, by = c("Participant_id", "Timepoint" = "tp_lab"))

message("Final Joined rows: ", nrow(final_df))
if (nrow(final_df) > 0) {
    print(head(final_df %>% select(Participant_id, Timepoint, ST, Infection_Status)))
}

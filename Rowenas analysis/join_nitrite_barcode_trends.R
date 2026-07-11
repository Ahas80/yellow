#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(haven)
  library(readr)
  library(stringr)
  library(tidyr)
})

base_file <- file.path("outputs", "nitrate_blast", "nitrate_trend_ready_needs_nitrite.csv")
nitrite_file <- "SPSS_nitrite_results"
out_dir <- file.path("outputs", "nitrate_blast")

required_nitrite_cols <- c("participant_id", "timepoint", "urinestick_nitriet")
required_base_cols <- c(
  "Isolate_ID", "Participant_id", "tp_lab", "Timepoint", "Event_type",
  "sequencing_run", "sequencing_barcode", "nitrate_profile_01",
  "genes_present", "genes_screened", "narA_complete",
  "narZ_complete_available_refs", "nar_regulators_complete",
  "gene_CsgD", "gene_NarG", "gene_NarH", "gene_NarI", "gene_NarJ",
  "gene_NarL", "gene_NarP", "gene_NarQ", "gene_NarV", "gene_NarX",
  "gene_NarY", "gene_NarZ"
)

normalise_pid <- function(x) {
  x_chr <- trimws(as.character(x))
  is_intish <- str_detect(x_chr, "^[0-9]+([.]0+)?$")
  x_num <- suppressWarnings(as.integer(as.numeric(x_chr)))
  ifelse(is_intish & !is.na(x_num), as.character(x_num), x_chr)
}

assert_required_cols <- function(df, required, context) {
  missing <- setdiff(required, names(df))
  if (length(missing) > 0) {
    stop(context, " is missing required columns: ", paste(missing, collapse = ", "))
  }
}

summarise_nitrite <- function(df, ...) {
  df %>%
    group_by(...) %>%
    summarise(
      total_rows = n(),
      matched_binary_n = sum(.data$nitrite_binary %in% c("positive", "negative"), na.rm = TRUE),
      positive_n = sum(.data$nitrite_binary == "positive", na.rm = TRUE),
      negative_n = sum(.data$nitrite_binary == "negative", na.rm = TRUE),
      missing_or_unmatched_n = .data$total_rows - .data$matched_binary_n,
      positive_pct_among_matched = ifelse(
        .data$matched_binary_n > 0,
        round(100 * .data$positive_n / .data$matched_binary_n, 1),
        NA_real_
      ),
      distinct_participants = n_distinct(.data$Participant_id),
      distinct_isolates = n_distinct(.data$Isolate_ID),
      .groups = "drop"
    ) %>%
    arrange(desc(.data$matched_binary_n), desc(.data$positive_pct_among_matched))
}

if (!file.exists(base_file)) {
  stop("Cannot find FASTA/barcode base table: ", base_file)
}
if (!file.exists(nitrite_file)) {
  stop("Cannot find nitrite source file: ", nitrite_file)
}

message("Reading FASTA/barcode base table...")
base <- read_csv(base_file, show_col_types = FALSE)
assert_required_cols(base, required_base_cols, base_file)

message("Reading nitrite source RDS...")
nitrite_raw <- readRDS(nitrite_file)
if (!is.data.frame(nitrite_raw)) {
  stop("Nitrite source is not a data frame/tibble: ", nitrite_file)
}
if (!identical(names(nitrite_raw), required_nitrite_cols)) {
  stop(
    "Nitrite source must have exactly these columns, in order: ",
    paste(required_nitrite_cols, collapse = ", "),
    "\nFound: ", paste(names(nitrite_raw), collapse = ", ")
  )
}

nitrite_clean <- nitrite_raw %>%
  transmute(
    Participant_id = normalise_pid(.data$participant_id),
    spss_timepoint = as.character(.data$timepoint),
    spss_timepoint_integer = suppressWarnings(as.integer(.data$timepoint)),
    tp_lab = ifelse(
      !is.na(.data$spss_timepoint_integer),
      paste0("T", .data$spss_timepoint_integer - 1L),
      NA_character_
    ),
    urinestick_nitriet_code = suppressWarnings(as.numeric(.data$urinestick_nitriet)),
    urinestick_nitriet_label = as.character(haven::as_factor(.data$urinestick_nitriet)),
    nitrite_binary = case_when(
      .data$urinestick_nitriet_code == 1 |
        .data$urinestick_nitriet_label == "Positief" ~ "positive",
      .data$urinestick_nitriet_code == 2 |
        .data$urinestick_nitriet_label == "Negatief" ~ "negative",
      is.na(.data$urinestick_nitriet_code) |
        is.na(.data$urinestick_nitriet_label) |
        .data$urinestick_nitriet_label %in% c(
          "Not done", "Asked but unknown", "Not asked",
          "Not applicable", "Measurement failed"
        ) ~ "missing",
      TRUE ~ "missing"
    )
  )

duplicate_keys <- nitrite_clean %>%
  count(.data$Participant_id, .data$tp_lab, name = "n") %>%
  filter(.data$n > 1)

if (nrow(duplicate_keys) > 0) {
  write_csv(duplicate_keys, file.path(out_dir, "nitrite_duplicate_participant_timepoint_keys.csv"))
  stop(
    "Nitrite source has duplicate Participant_id + tp_lab keys. ",
    "See outputs/nitrate_blast/nitrite_duplicate_participant_timepoint_keys.csv"
  )
}

base_prepped <- base %>%
  select(-any_of(c("dipstick_nitrite", "nitrite_join_note"))) %>%
  mutate(
    Participant_id = normalise_pid(.data$Participant_id),
    sequencing_run = as.character(.data$sequencing_run),
    sequencing_barcode = str_pad(as.character(.data$sequencing_barcode), width = 2, pad = "0"),
    run_barcode = paste0(.data$sequencing_run, "_barcode", .data$sequencing_barcode),
    all_12_screened_genes_present = .data$genes_present == .data$genes_screened
  )

joined <- base_prepped %>%
  left_join(
    nitrite_clean,
    by = c("Participant_id", "tp_lab"),
    relationship = "many-to-one"
  ) %>%
  mutate(
    nitrite_join_status = case_when(
      .data$nitrite_binary %in% c("positive", "negative") ~ "matched_binary_result",
      .data$nitrite_binary == "missing" ~ "matched_missing_or_unclear_result",
      TRUE ~ "no_matching_nitrite_row"
    ),
    nitrite_positive = case_when(
      .data$nitrite_binary == "positive" ~ TRUE,
      .data$nitrite_binary == "negative" ~ FALSE,
      TRUE ~ NA
    )
  ) %>%
  relocate(
    run_barcode, nitrite_binary, nitrite_positive,
    nitrite_join_status, urinestick_nitriet_code,
    urinestick_nitriet_label, spss_timepoint,
    spss_timepoint_integer,
    .after = sequencing_barcode
  )

unmatched <- joined %>%
  filter(.data$nitrite_join_status != "matched_binary_result") %>%
  arrange(.data$Event_type, .data$Participant_id, .data$tp_lab, .data$Isolate_ID)

by_run <- summarise_nitrite(joined, .data$sequencing_run)

by_barcode_number <- summarise_nitrite(joined, .data$sequencing_barcode) %>%
  rename(barcode_number = sequencing_barcode)

by_run_barcode <- summarise_nitrite(
  joined,
  .data$sequencing_run,
  .data$sequencing_barcode,
  .data$run_barcode
)

by_profile <- summarise_nitrite(
  joined,
  .data$nitrate_profile_01,
  .data$genes_absent,
  .data$genes_present,
  .data$genes_screened,
  .data$all_12_screened_genes_present,
  .data$narA_complete,
  .data$narZ_complete_available_refs,
  .data$nar_regulators_complete
)

join_qc <- tibble(
  metric = c(
    "base_rows",
    "nitrite_source_rows",
    "nitrite_unique_participant_timepoint_keys",
    "joined_rows",
    "matched_binary_rows",
    "positive_rows",
    "negative_rows",
    "unmatched_or_missing_rows",
    "unmatched_uti_event_rows",
    "unmatched_routine_rows"
  ),
  value = c(
    nrow(base),
    nrow(nitrite_clean),
    n_distinct(paste(nitrite_clean$Participant_id, nitrite_clean$tp_lab, sep = "__")),
    nrow(joined),
    sum(joined$nitrite_binary %in% c("positive", "negative"), na.rm = TRUE),
    sum(joined$nitrite_binary == "positive", na.rm = TRUE),
    sum(joined$nitrite_binary == "negative", na.rm = TRUE),
    nrow(unmatched),
    sum(unmatched$Event_type == "UTI_event", na.rm = TRUE),
    sum(unmatched$Event_type == "Routine", na.rm = TRUE)
  )
)

write_csv(joined, file.path(out_dir, "nitrate_barcode_nitrite_joined.csv"))
write_csv(unmatched, file.path(out_dir, "nitrate_barcode_nitrite_unmatched.csv"))
write_csv(by_run, file.path(out_dir, "nitrite_by_sequencing_run.csv"))
write_csv(by_barcode_number, file.path(out_dir, "nitrite_by_barcode_number.csv"))
write_csv(by_run_barcode, file.path(out_dir, "nitrite_by_run_barcode.csv"))
write_csv(by_profile, file.path(out_dir, "nitrite_by_nitrate_gene_profile.csv"))
write_csv(join_qc, file.path(out_dir, "nitrite_barcode_join_qc.csv"))

message("Done. Wrote barcode-nitrite CSV outputs to: ", out_dir)
print(join_qc)

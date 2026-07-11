#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(purrr)
  library(readr)
  library(readxl)
  library(stringr)
  library(tidyr)
})

source("00_config.R")

input_dir <- file.path(DIR_ROOT, "data", "inputs")
qc_dir <- DIR_QC
dir.create(input_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(qc_dir, recursive = TRUE, showWarnings = FALSE)

overview_path <- FILE_OVERVIEW_XLSX
overview_sheet <- OVERVIEW_SHEET

individual_sources <- tibble::tribble(
  ~Batch, ~source_path, ~skip_rows,
  1L, file.path(input_dir, "list Ecoli YELLOW sequencing - batch 1 - uitgebreid.xlsx"), 5L,
  2L, file.path(input_dir, "list Ecoli YELLOW sequencing - batch 2 - uitgebreid - correctie.xlsx"), 3L,
  3L, file.path(input_dir, "Lijst Ecoli YELLOW sequencing - batch 3 - UPDATE (definitief).xlsx"), 3L
)

canonical_cols <- c(
  "Participant_id", "UTI_Label", "Collection_Date", "Spec", "Obj",
  "isolate_ID", "Organism", "Beoord", "CFU_Count", "Archive", "Timepoint",
  "Population", "UWI#", "Urine collection method",
  "S&S Dysuria", "S&S Urgency", "S&S Frequency", "S&S Incontinence",
  "S&S Pus", "S&S Flank_pain", "S&S Suprapubic pain", "S&S Fever",
  "S&S Chills", "S&S Delirium", "S&S Other", "No S&S"
)

format_date <- function(x) {
  if (inherits(x, "POSIXt") || inherits(x, "Date")) {
    return(format(as.Date(x), "%d/%m/%Y"))
  }
  x_chr <- as.character(x)
  suppressWarnings({
    num <- as.numeric(x_chr)
  })
  out <- x_chr
  excel_like <- !is.na(num) & str_detect(x_chr, "^[0-9]+(\\.[0-9]+)?$")
  out[excel_like] <- format(as.Date(num[excel_like], origin = "1899-12-30"), "%d/%m/%Y")
  out
}

as_chr <- function(x) {
  out <- as.character(x)
  out[is.na(out) | out == "NA"] <- NA_character_
  out
}

rename_canonical <- function(df) {
  df %>%
    rename(
      UTI_Label = any_of("UWIsticker"),
      Collection_Date = any_of("AfnDat"),
      isolate_ID = any_of("Isolaat id"),
      Organism = any_of("Organisme"),
      CFU_Count = any_of("Kiemgetal"),
      Archive = any_of("Archief"),
      Timepoint = any_of("Meetmoment"),
      Population = any_of("Populatie"),
      `Urine collection method` = any_of("Urine opvang methode"),
      `S&S Dysuria` = any_of("S&S Dysurie"),
      `S&S Urgency` = any_of("S&S Aandrang"),
      `S&S Frequency` = any_of("S&S Frequentie"),
      `S&S Incontinence` = any_of("S&S Incontinentie"),
      `S&S Flank_pain` = any_of("S&S Flankpijn"),
      `S&S Suprapubic pain` = any_of("S&S Suprapub pijn"),
      `S&S Fever` = any_of("S&S Koorts"),
      `S&S Chills` = any_of("S&S Rillingen"),
      `S&S Delirium` = any_of("S&S Delier"),
      `S&S Other` = any_of("S&S anders"),
      `No S&S` = any_of("Geen S&S")
    )
}

normalise_batch_frame <- function(df) {
  df <- rename_canonical(df)
  for (col in canonical_cols) {
    if (!col %in% names(df)) {
      df[[col]] <- NA_character_
    }
  }

  df %>%
    mutate(Collection_Date = format_date(Collection_Date)) %>%
    mutate(across(all_of(canonical_cols), as_chr)) %>%
    select(all_of(canonical_cols))
}

read_individual_source <- function(Batch, source_path, skip_rows) {
  if (!file.exists(source_path)) {
    stop("Missing expected source workbook: ", source_path)
  }
  read_excel(source_path, skip = skip_rows) %>%
    normalise_batch_frame() %>%
    mutate(Batch = Batch, .before = 1) %>%
    filter(!is.na(isolate_ID), !is.na(Participant_id)) %>%
    filter(!str_detect(Participant_id, regex("^NIET\\s+MEEGENOMEN", ignore_case = TRUE)))
}

overview_raw <- read_excel(overview_path, sheet = overview_sheet, skip = 3)
overview <- overview_raw %>%
  rename_canonical() %>%
  mutate(Batch = as.integer(Batch)) %>%
  filter(Batch %in% 1:6) %>%
  filter(!is.na(isolate_ID), !is.na(Participant_id)) %>%
  mutate(Collection_Date = format_date(Collection_Date)) %>%
  mutate(across(-Batch, as_chr))

overview_canonical <- overview
for (col in canonical_cols) {
  if (!col %in% names(overview_canonical)) {
    overview_canonical[[col]] <- NA_character_
  }
}

overview_canonical <- overview_canonical %>%
  select(Batch, all_of(canonical_cols))

individual <- bind_rows(purrr::pmap(individual_sources, read_individual_source))

source_supplements <- individual %>%
  select(Batch, isolate_ID, Population_source = Population, UWI_source = `UWI#`)

written <- overview_canonical %>%
  left_join(source_supplements, by = c("Batch", "isolate_ID")) %>%
  mutate(
    Population = coalesce(Population_source, Population),
    `UWI#` = coalesce(UWI_source, `UWI#`)
  ) %>%
  select(Batch, all_of(canonical_cols)) %>%
  arrange(Batch, Collection_Date, Participant_id, isolate_ID)

source_vs_overview <- individual %>%
  filter(Batch %in% 1:3) %>%
  select(Batch, source_Participant_id = Participant_id, source_Timepoint = Timepoint,
         source_UTI_Label = UTI_Label, isolate_ID, source_Population = Population,
         source_UWI = `UWI#`) %>%
  full_join(
    overview_canonical %>%
      filter(Batch %in% 1:3) %>%
      select(Batch, overview_Participant_id = Participant_id,
             overview_Timepoint = Timepoint, overview_UTI_Label = UTI_Label,
             isolate_ID),
    by = c("Batch", "isolate_ID")
  ) %>%
  mutate(
    comparison = case_when(
      is.na(overview_Participant_id) ~ "source_row_not_in_clean_rc_overview",
      is.na(source_Participant_id) ~ "clean_rc_overview_row_not_in_individual_workbook",
      source_Participant_id != overview_Participant_id |
        source_Timepoint != overview_Timepoint |
        coalesce(source_UTI_Label, "") != coalesce(overview_UTI_Label, "") ~
        "clean_rc_overview_corrects_source_row",
      TRUE ~ "matches"
    )
  ) %>%
  arrange(Batch, comparison, isolate_ID)

write_csv(source_vs_overview, file.path(qc_dir, "batch_input_source_vs_clean_rc_overview.csv"))

population_summary <- written %>%
  count(Batch, name = "written_rows") %>%
  left_join(
    individual %>% count(Batch, name = "individual_valid_rows"),
    by = "Batch"
  ) %>%
  mutate(
    workbook_source = case_when(
      Batch %in% 1:3 ~ individual_sources$source_path[match(Batch, individual_sources$Batch)],
      TRUE ~ overview_path
    ),
    overview_rows_for_batch = as.integer(table(factor(overview_canonical$Batch, levels = 1:6))[as.character(Batch)])
  ) %>%
  select(Batch, workbook_source, individual_valid_rows, overview_rows_for_batch, written_rows)

write_csv(population_summary, file.path(qc_dir, "batch_input_population_summary.csv"))

for (batch in sort(unique(written$Batch))) {
  out_path <- file.path(input_dir, paste0("batch", batch, ".csv"))
  written %>%
    filter(Batch == batch) %>%
    select(all_of(canonical_cols)) %>%
    write_csv(out_path, na = "")
  message("Wrote ", out_path, " from canonical workbook inputs")
}

message("Wrote audit: ", file.path(qc_dir, "batch_input_population_summary.csv"))
message("Wrote audit: ", file.path(qc_dir, "batch_input_source_vs_clean_rc_overview.csv"))
print(population_summary)
print(source_vs_overview %>% count(Batch, comparison, name = "n"))

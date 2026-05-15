#!/usr/bin/env Rscript
# ==============================================================================
# 00_make_assembly_metadata.r
# ------------------------------------------------------------------------------
# Build the authoritative assembly-level metadata table from:
#   1. OVERVIEW E.coli batch 1-6 - CLEAN.xlsx (expected isolate universe)
#   2. data/inputs/batch1.csv ... batch6.csv (metadata supplements only)
#   3. ont-yellow-routine-fastas/ (canonical FASTA discovery helper)
#
# This script deliberately does not add unlinked FASTAs as biological episodes.
# FASTAs absent from both overview and batch CSVs are reported for provenance and
# excluded from episode-level analyses until a trusted mapping is supplied.
# ==============================================================================

source("00_config.R")
suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
  library(tidyr)
  library(fs)
})

msg("Starting 00_make_assembly_metadata.r")

overview <- load_overview_spreadsheet()
if (is.null(overview)) {
  stop("Overview spreadsheet not found. This pipeline requires the master overview spreadsheet.")
}

expected_df <- overview %>%
  select(any_of(c(
    "Batch", "Participant_id", "Isolate_ID", "Timepoint", "Organism",
    "Beoordeling", "Kiemgetal"
  ))) %>%
  mutate(
    Batch = suppressWarnings(as.integer(Batch)),
    Participant_id = as.character(Participant_id),
    Isolate_ID = as.character(Isolate_ID),
    Timepoint = as.character(Timepoint),
    tp_lab = normalise_timepoint_preserve_events(Timepoint),
    Event_type = episode_event_type(tp_lab),
    overview_row_number = row_number()
  ) %>%
  filter(!is.na(Isolate_ID), Isolate_ID != "")

msg("Loaded %d expected E. coli isolates from overview spreadsheet.", nrow(expected_df))

# ------------------------------------------------------------------------------
# Metadata provenance: overview is authoritative; batch CSVs supplement rows that
# are already in the overview.
# ------------------------------------------------------------------------------
clinical_csvs <- get_available_clinical_csvs()
manifest_paths <- c(FILE_OVERVIEW_XLSX, clinical_csvs)
manifest_roles <- c("authoritative_expected_isolate_universe", rep("clinical_metadata_supplement", length(clinical_csvs)))
metadata_manifest <- tibble(
  role = manifest_roles,
  path = normalizePath(manifest_paths, winslash = "/", mustWork = FALSE),
  exists = file.exists(manifest_paths),
  file_size = ifelse(file.exists(manifest_paths), file.size(manifest_paths), NA_real_),
  modified_time = ifelse(file.exists(manifest_paths), as.character(file.info(manifest_paths)$mtime), NA_character_),
  md5 = ifelse(file.exists(manifest_paths), unname(tools::md5sum(manifest_paths)), NA_character_)
)
write_csv(metadata_manifest, file.path(DIR_QC, "metadata_input_manifest.csv"))

read_batch_supplement <- function(path, batch_id) {
  df <- read_csv(path, show_col_types = FALSE) %>% mutate(across(everything(), as.character))
  if ("isolate_ID" %in% names(df) && !"Isolate_ID" %in% names(df)) {
    df <- rename(df, Isolate_ID = isolate_ID)
  }
  if (!"Isolate_ID" %in% names(df)) {
    return(tibble())
  }

  out <- df %>%
    mutate(
      Batch_csv = batch_id,
      Participant_id_csv = if ("Participant_id" %in% names(.)) as.character(Participant_id) else NA_character_,
      Timepoint_csv = if ("Timepoint" %in% names(.)) as.character(Timepoint) else NA_character_,
      tp_lab_csv = normalise_timepoint_preserve_events(Timepoint_csv),
      Collection_Date = if ("Collection_Date" %in% names(.)) as.character(Collection_Date) else NA_character_,
      UTI_Label = if ("UTI_Label" %in% names(.)) as.character(UTI_Label) else NA_character_,
      Clinical_Beoord = if ("Beoord" %in% names(.)) as.character(Beoord) else NA_character_,
      Clinical_CFU_Count = if ("CFU_Count" %in% names(.)) as.character(CFU_Count) else NA_character_,
      Clinical_Organism = if ("Organism" %in% names(.)) as.character(Organism) else NA_character_,
      Urine_collection_method = if ("Urine collection method" %in% names(.)) as.character(`Urine collection method`) else NA_character_,
      Population = if ("Population" %in% names(.)) as.character(Population) else NA_character_,
      Spec = if ("Spec" %in% names(.)) as.character(Spec) else NA_character_,
      Obj = if ("Obj" %in% names(.)) as.character(Obj) else NA_character_,
      Archive = if ("Archive" %in% names(.)) as.character(Archive) else NA_character_,
      UWI_number = if ("UWI#" %in% names(.)) as.character(`UWI#`) else NA_character_
    ) %>%
    select(
      Isolate_ID, Batch_csv, Participant_id_csv, Timepoint_csv, tp_lab_csv,
      Collection_Date, UTI_Label, Clinical_Beoord, Clinical_CFU_Count,
      Clinical_Organism, Urine_collection_method, Population, Spec, Obj,
      Archive, UWI_number
    ) %>%
    filter(!is.na(Isolate_ID), Isolate_ID != "")

  out
}

batch_supp <- bind_rows(lapply(seq_along(clinical_csvs), function(i) {
  csv_name <- basename(clinical_csvs[i])
  batch_id <- suppressWarnings(as.integer(str_extract(csv_name, "[0-9]+")))
  read_batch_supplement(clinical_csvs[i], batch_id)
}))

if (nrow(batch_supp) > 0) {
  batch_supp <- batch_supp %>%
    group_by(Isolate_ID) %>%
    summarise(across(everything(), ~ paste(unique(na.omit(.x)), collapse = ";")), .groups = "drop") %>%
    mutate(across(everything(), ~ na_if(.x, "")))
}

crosswalk_ids <- union(expected_df$Isolate_ID, batch_supp$Isolate_ID)
crosswalk <- tibble(Isolate_ID = crosswalk_ids) %>%
  left_join(
    expected_df %>%
      transmute(
        Isolate_ID,
        in_overview = TRUE,
        overview_Batch = as.character(Batch),
        overview_Participant_id = Participant_id,
        overview_Timepoint = Timepoint,
        overview_tp_lab = tp_lab
      ),
    by = "Isolate_ID"
  ) %>%
  left_join(
    batch_supp %>%
      transmute(
        Isolate_ID,
        in_batch_csv = TRUE,
        batch_csv_Batch = Batch_csv,
        batch_csv_Participant_id = Participant_id_csv,
        batch_csv_Timepoint = Timepoint_csv,
        batch_csv_tp_lab = tp_lab_csv
      ),
    by = "Isolate_ID"
  ) %>%
  mutate(
    in_overview = coalesce(in_overview, FALSE),
    in_batch_csv = coalesce(in_batch_csv, FALSE),
    provenance_status = case_when(
      in_overview & in_batch_csv ~ "overview_and_batch_csv",
      in_overview & !in_batch_csv ~ "overview_only",
      !in_overview & in_batch_csv ~ "batch_csv_only_not_authoritative",
      TRUE ~ "unlinked"
    ),
    participant_disagree = in_overview & in_batch_csv &
      !is.na(batch_csv_Participant_id) & overview_Participant_id != batch_csv_Participant_id,
    timepoint_disagree = in_overview & in_batch_csv &
      !is.na(batch_csv_tp_lab) & overview_tp_lab != batch_csv_tp_lab
  )
write_csv(crosswalk, file.path(DIR_QC, "metadata_input_crosswalk_audit.csv"))

batch_only <- crosswalk %>% filter(provenance_status == "batch_csv_only_not_authoritative")
if (nrow(batch_only) > 0) {
  msg("WARNING: %d isolate ID(s) appear in batch CSVs but not overview; they are not added automatically.", nrow(batch_only))
}

# ------------------------------------------------------------------------------
# Canonical FASTA discovery.  The full scan is written for audit; only top-level
# candidate input FASTAs are eligible for assembly_metadata.csv.
# ------------------------------------------------------------------------------
fasta_scan_all <- discover_project_fastas(DIR_FASTAS, recursive = TRUE, include_excluded = TRUE)
write_csv(fasta_scan_all, file.path(DIR_QC, "metadata_fasta_discovery_manifest.csv"))

candidate_fastas <- fasta_scan_all %>% filter(include_in_metadata)

# Batch 5/6 names sometimes omit a unique trailing replicate suffix present in
# the overview. Repair only when stripping that suffix maps to exactly one
# expected isolate.
expected_suffix_lookup <- expected_df %>%
  mutate(file_id_without_suffix = str_remove(Isolate_ID, "-\\d+$")) %>%
  filter(file_id_without_suffix != Isolate_ID) %>%
  group_by(file_id_without_suffix) %>%
  summarise(
    matched_isolate_id = first(Isolate_ID),
    n_expected_ids = n_distinct(Isolate_ID),
    .groups = "drop"
  ) %>%
  filter(n_expected_ids == 1)

candidate_fastas <- candidate_fastas %>%
  left_join(expected_suffix_lookup, by = c("Isolate_ID" = "file_id_without_suffix")) %>%
  mutate(
    Isolate_ID_raw = Isolate_ID,
    Isolate_ID = case_when(
      !Isolate_ID %in% expected_df$Isolate_ID &
        !str_detect(Isolate_ID, "-\\d+$") &
        !is.na(matched_isolate_id) ~ matched_isolate_id,
      TRUE ~ Isolate_ID
    ),
    naming_normalisation_note = case_when(
      Isolate_ID != Isolate_ID_raw ~ paste0("mapped omitted suffix from ", Isolate_ID_raw, " to ", Isolate_ID),
      TRUE ~ NA_character_
    )
  ) %>%
  select(-matched_isolate_id, -n_expected_ids)

found_df <- candidate_fastas %>%
  mutate(
    # Assembly_ID is one row per FASTA.  Some batches contain both .fasta and
    # .fna files with the same stem, so the extension is kept in the stable ID
    # while Assembly_Base_ID preserves the old stem for Prokka/GFF matching.
    Assembly_Base_ID = str_remove(file_name, regex("\\.(fasta|fa|fna)(\\.gz)?$", ignore_case = TRUE)),
    Assembly_ID = paste0(Assembly_Base_ID, "__", str_replace_all(extension, "\\.", "_")),
    fasta_path = full_path,
    file_exists = file.exists(full_path),
    usable_fasta = file_exists
  )

found_unexpected <- found_df %>%
  filter(!Isolate_ID %in% expected_df$Isolate_ID) %>%
  mutate(
    in_batch_csv = Isolate_ID %in% batch_supp$Isolate_ID,
    fasta_source_status = case_when(
      in_batch_csv ~ "unexpected_in_batch_csv_only",
      fasta_class == "candidate_input_fasta" ~ "unexpected_unlinked",
      TRUE ~ fasta_class
    )
  )

excluded_fastas <- fasta_scan_all %>%
  filter(!include_in_metadata) %>%
  mutate(
    in_overview = Isolate_ID %in% expected_df$Isolate_ID,
    in_batch_csv = Isolate_ID %in% batch_supp$Isolate_ID,
    fasta_source_status = case_when(
      fasta_class == "legacy_or_quarantine" ~ "legacy_or_quarantine",
      fasta_class == "outside_expected_input_directory" ~ "outside_expected_input_directories",
      fasta_class == "derived_output_or_cache" ~ "derived_output_or_cache",
      TRUE ~ fasta_class
    )
  )

unlinked_report <- bind_rows(
  found_unexpected %>%
    transmute(full_path, relative_path, file_name, Isolate_ID, Assembler,
              fasta_class, fasta_source_status, in_overview = FALSE, in_batch_csv,
              exclusion_reason = "Candidate FASTA does not map to overview expected isolate universe"),
  excluded_fastas %>%
    transmute(full_path, relative_path, file_name, Isolate_ID, Assembler,
              fasta_class, fasta_source_status, in_overview, in_batch_csv,
              exclusion_reason)
)
write_csv(unlinked_report, file.path(DIR_QC, "unlinked_unexpected_fastas.csv"))

msg("Found %d candidate FASTA files in %s", nrow(candidate_fastas), DIR_FASTAS)
msg("Full recursive FASTA audit universe: %d files (%d excluded from metadata).",
    nrow(fasta_scan_all), sum(!fasta_scan_all$include_in_metadata))

# ------------------------------------------------------------------------------
# Match expected isolates to candidate FASTAs.  This intentionally expands to
# assembly-level rows when flye and longcycler are both present.
# ------------------------------------------------------------------------------
meta <- expected_df %>%
  left_join(batch_supp, by = "Isolate_ID") %>%
  left_join(
    found_df %>%
      select(
        Isolate_ID, Assembly_ID, Assembly_Base_ID, file_name, full_path, fasta_path, relative_path,
        Assembler, file_exists, usable_fasta, Isolate_ID_raw,
        naming_normalisation_note
      ),
    by = "Isolate_ID"
  ) %>%
  mutate(
    found = !is.na(file_name),
    full_path = ifelse(found, full_path, NA_character_),
    fasta_path = ifelse(found, fasta_path, NA_character_),
    file_exists = coalesce(file_exists, FALSE),
    usable_fasta = coalesce(usable_fasta, FALSE),
    Assembler = ifelse(found, Assembler, NA_character_),
    assembler = Assembler,
    Assembly_ID = ifelse(found, Assembly_ID, paste0("missing__", Isolate_ID)),
    Assembly_Base_ID = ifelse(found, Assembly_Base_ID, NA_character_),
    metadata_source_status = case_when(
      found & !is.na(Collection_Date) ~ "overview_plus_batch_csv_plus_fasta",
      found ~ "overview_plus_fasta",
      !found & !is.na(Collection_Date) ~ "overview_plus_batch_csv_no_fasta",
      TRUE ~ "overview_only_no_fasta"
    ),
    Episode_ID = build_episode_id(., timepoint_col = "tp_lab", event_col = "Event_type",
                                  date_col = "Collection_Date")
  )

missing_expected <- meta %>% filter(!found)

msg("")
msg("=== Validation Summary ===")
msg("Total expected isolates: %d", nrow(expected_df))
msg("Candidate FASTAs found: %d", nrow(candidate_fastas))
msg("Successfully matched assemblies: %d", sum(meta$found))
msg("Expected isolate rows without FASTA: %d", nrow(missing_expected))
msg("Unexpected/unlinked candidate FASTAs: %d", nrow(found_unexpected))

if (nrow(missing_expected) > 0) {
  missing_file <- file.path(DIR_QC, "00_missing_expected_assemblies.csv")
  write_csv(missing_expected, missing_file)
  msg("Wrote missing assemblies report to: %s", missing_file)

  batch_summary <- missing_expected %>% count(Batch, name = "Missing_Count")
  msg("Missing FASTAs by batch:")
  for (i in seq_len(nrow(batch_summary))) {
    msg("  Batch %s: %d missing", as.character(batch_summary$Batch[i]), batch_summary$Missing_Count[i])
  }
}

if (nrow(found_unexpected) > 0) {
  unex_file <- file.path(DIR_QC, "00_unexpected_assemblies.csv")
  write_csv(found_unexpected, unex_file)
  msg("Wrote unexpected candidate FASTA report to: %s", unex_file)
}

# ------------------------------------------------------------------------------
# Compute lightweight FASTA metrics for rows with files.
# ------------------------------------------------------------------------------
if (requireNamespace("Biostrings", quietly = TRUE)) {
  msg("Calculating metrics for %d valid assemblies...", sum(meta$found))
  summarise_fasta <- function(fp) {
    if (is.na(fp) || !file.exists(fp)) {
      return(tibble(num_contigs = NA_integer_, total_bases = NA_real_, gc_content = NA_real_))
    }
    tryCatch({
      x <- Biostrings::readDNAStringSet(fp)
      af <- colSums(Biostrings::alphabetFrequency(x, baseOnly = TRUE))
      tibble(
        num_contigs = length(x),
        total_bases = sum(Biostrings::width(x)),
        gc_content = round((af["G"] + af["C"]) / sum(af) * 100, 2)
      )
    }, error = function(e) tibble(num_contigs = NA_integer_, total_bases = NA_real_, gc_content = NA_real_))
  }

  metrics_tbl <- bind_rows(lapply(meta$full_path, summarise_fasta))
  meta <- bind_cols(meta, metrics_tbl)
} else {
  msg("Biostrings package not available. Skipping FASTA metric calculation.")
  meta <- meta %>% mutate(num_contigs = NA_integer_, total_bases = NA_real_, gc_content = NA_real_)
}

out_meta <- meta %>%
  select(
    Assembly_ID, Assembly_Base_ID, Isolate_ID, Participant_id, tp_lab, Timepoint, Event_type,
    Episode_ID, Batch, Assembler, assembler, file_name, full_path, fasta_path,
    file_exists, usable_fasta, found, metadata_source_status,
    Collection_Date, UTI_Label, Clinical_Beoord, Clinical_CFU_Count,
    Clinical_Organism, Urine_collection_method, Population, Spec, Obj,
    Archive, UWI_number, Organism, Beoordeling, Kiemgetal,
    num_contigs, total_bases, gc_content, relative_path, Isolate_ID_raw,
    naming_normalisation_note, overview_row_number
  ) %>%
  arrange(Batch, Participant_id, tp_lab, Isolate_ID, Assembler, Assembly_ID)

assert_unique_keys(
  out_meta %>% filter(found),
  c("Assembly_ID"),
  context = "assembly_metadata found assemblies",
  out_path = file.path(DIR_QC, "assembly_metadata_duplicate_assembly_ids.csv")
)

write_csv(out_meta, FILE_METADATA)
write_csv(out_meta, file.path(DIR_RESULTS, "assembly_metadata.csv"))
write_csv(out_meta, FILE_ASSEMBLY_META_ALL)

append_denominator_summary(expected_df, "00_make_assembly_metadata.r", "expected_overview",
                           "isolate", FILE_OVERVIEW_XLSX,
                           "Overview spreadsheet is authoritative expected-isolate universe")
append_denominator_summary(out_meta, "00_make_assembly_metadata.r", "assembly_metadata",
                           "assembly", FILE_METADATA,
                           "Assembly-level rows; flye/longcycler alternatives are not independent biological episodes")
append_denominator_summary(found_unexpected, "00_make_assembly_metadata.r", "unexpected_fastas",
                           "assembly", file.path(DIR_QC, "unlinked_unexpected_fastas.csv"),
                           "Unlinked candidate FASTAs excluded from episode-level analyses")

msg("✓ Master metadata generated: %s", FILE_METADATA)
msg("✓ Pipeline is ready. Downstream scripts will gracefully skip the %d missing expected isolate rows.", nrow(missing_expected))

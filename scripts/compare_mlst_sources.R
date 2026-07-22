#!/usr/bin/env Rscript

# Compare Longcycler-local MLST calls with explicitly source-labelled
# Longcycler provider SeqSphere/wgMLST calls linked to the exact selected
# Longcycler manifest key and FASTA path/SHA-256.

suppressPackageStartupMessages({
  library(dplyr)
  library(purrr)
  library(readr)
  library(stringr)
  library(tibble)
  library(tidyr)
})

if (!file.exists("00_config.R")) {
  stop("Run this script from the project root: /Users/Aamir/Desktop/rUTIs")
}
source("00_config.R")

out_dir <- file.path("results", "mlst_source_comparison")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
legacy_excluded_audit <- file.path(out_dir, "provider_mlst_non_longcycler_excluded.csv")
if (file.exists(legacy_excluded_audit)) unlink(legacy_excluded_audit)

local_mlst_path <- FILE_MLST_LOCAL_CANONICAL
canonical_path <- FILE_ANALYSIS_ASSEMBLY_MANIFEST

provider_files <- tibble(
  provider_source = c(
    "batch1_2_first_158_bsr_mlst",
    "batch1_2_repeat_34_bsr_mlst",
    "batch3_seqsphere_bsr_mlst",
    "batch4_6_bsr_mlst_longcycler"
  ),
  provider_scope = c("Batch1+2 first 158", "Batch1+2 repeat 34", "Batch3", "Batch4-6"),
  expected_batches = c("1;2", "1;2", "3", "4;5;6"),
  provider_role = rep("selected_manifest_key_linked", 4),
  provider_path = c(
    file.path(
      "archive", "cleanup_2026-05-14", "put_away_later", "zipped files",
      "Batch1+2", "Data eerste 158 isolaten", "wgMLST_Eco_BSR_+_MLST.csv"
    ),
    file.path(
      "archive", "cleanup_2026-05-14", "put_away_later", "zipped files",
      "Batch1+2", "Data laatste 34 isolaten", "herhaling_wgMLST_Eco_BSR_+_MLST.csv"
    ),
    file.path(
      "archive", "cleanup_2026-05-14", "put_away_later", "zipped files",
      "Batch3", "seqsphere", "wgMLST_Eco_BSR_+_MLST.csv"
    ),
    file.path(
      "archive", "cleanup_2026-05-14", "put_away_later", "zipped files",
      "batch4-6", "wgMLST", "wgMLST-Eco-BSR_MLST_longcycler.csv"
    )
  )
)

missing_tokens <- c(
  "", ".", "?", "-", "--", "na", "n/a", "nan", "missing", "unknown",
  "unk", "nt", "non-typable", "nontypable", "not typed"
)

normalise_st <- function(x) {
  y <- str_squish(as.character(x))
  y[is.na(x) | str_to_lower(y) %in% missing_tokens] <- NA_character_
  y
}

normalise_bool <- function(x) {
  if (is.logical(x)) return(x %in% TRUE)
  str_to_lower(str_squish(as.character(x))) %in% c("true", "t", "1", "yes", "y")
}

normalise_sample_base <- function(x) {
  y <- tools::file_path_sans_ext(basename(as.character(x)))
  y <- str_squish(y)
  y <- str_remove(y, regex("(_assembly)$", ignore_case = TRUE))
  y <- str_remove(y, regex("[-_](flye|longcycler)$", ignore_case = TRUE))
  y
}

detect_assembler <- function(x) {
  y <- tools::file_path_sans_ext(basename(as.character(x)))
  y <- str_remove(y, regex("(_assembly)$", ignore_case = TRUE))
  case_when(
    str_detect(y, regex("[-_]flye$", ignore_case = TRUE)) ~ "flye",
    str_detect(y, regex("[-_]longcycler$", ignore_case = TRUE)) ~ "longcycler",
    TRUE ~ NA_character_
  )
}

first_nonmissing <- function(x) {
  x <- x[!is.na(x) & nzchar(as.character(x))]
  if (length(x) == 0) NA_character_ else as.character(x[1])
}

collapse_set <- function(x) {
  vals <- sort(unique(as.character(x[!is.na(x) & nzchar(as.character(x))])))
  if (length(vals) == 0) NA_character_ else paste(vals, collapse = ";")
}

format_pct <- function(n, d) {
  if (is.na(d) || d == 0) return("NA")
  sprintf("%.1f%%", 100 * n / d)
}

read_provider_file <- function(provider_source, provider_scope, expected_batches, provider_role, provider_path) {
  if (!file.exists(provider_path)) {
    stop("Missing provider file: ", provider_path)
  }
  dat <- read.csv(
    provider_path,
    sep = ";",
    check.names = FALSE,
    quote = "\"",
    stringsAsFactors = FALSE
  )
  required <- c("Sample ID", "ST", "PercGoodTargets")
  missing_required <- setdiff(required, names(dat))
  if (length(missing_required) > 0) {
    stop(provider_path, " missing required columns: ", paste(missing_required, collapse = ", "))
  }

  tibble(
    provider_source = provider_source,
    provider_scope = provider_scope,
    expected_batches = expected_batches,
    provider_role = provider_role,
    provider_path = provider_path,
    provider_file = basename(provider_path),
    provider_sample_id = dat[["Sample ID"]],
    provider_assembler = detect_assembler(dat[["Sample ID"]]),
    provider_norm_id = normalise_sample_base(dat[["Sample ID"]]),
    provider_ST = normalise_st(dat[["ST"]]),
    provider_CC = if ("CC" %in% names(dat)) normalise_st(dat[["CC"]]) else NA_character_,
    provider_Profile = if ("Profile" %in% names(dat)) normalise_st(dat[["Profile"]]) else NA_character_,
    provider_PercGoodTargets = suppressWarnings(as.numeric(dat[["PercGoodTargets"]])),
    provider_qc_ge_90 = !is.na(provider_PercGoodTargets) & provider_PercGoodTargets >= 90,
    provider_qc_ge_95 = !is.na(provider_PercGoodTargets) & provider_PercGoodTargets >= 95,
    provider_has_classic_7_loci = all(c("adk", "fumC", "gyrB", "icd", "mdh", "purA", "recA") %in% names(dat)),
    provider_delimiter = ";",
    provider_n_rows_in_file = nrow(dat),
    provider_n_cols_in_file = ncol(dat)
  )
}

provider_file_inventory <- provider_files %>%
  mutate(
    exists = file.exists(provider_path),
    file_size_bytes = ifelse(exists, file.size(provider_path), NA_real_),
    modified_time = ifelse(exists, as.character(file.info(provider_path)$mtime), NA_character_)
  )
write_csv(provider_file_inventory, file.path(out_dir, "provider_file_inventory.csv"))

required_provider_missing <- provider_file_inventory %>% filter(!exists)
if (nrow(required_provider_missing) > 0) {
  stop(
    "One or more Longcycler-eligible provider source files are missing. See ",
    file.path(out_dir, "provider_file_inventory.csv")
  )
}

if (!file.exists(canonical_path)) {
  stop("Canonical assembly selection is mandatory: ", canonical_path)
}

canonical_all <- load_analysis_assemblies(canonical_path, require_files = TRUE)
required_canonical_cols <- c(
  "Batch", "Participant_id", "tp_lab", "Isolate_ID", "Assembly_ID",
  "selected_canonical", "QC_PASS", "full_path"
)
missing_canonical_cols <- setdiff(required_canonical_cols, names(canonical_all))
if (length(missing_canonical_cols) > 0) {
  stop("Canonical assembly selection lacks required column(s): ", paste(missing_canonical_cols, collapse = ", "))
}
if (!"Assembly_Base_ID" %in% names(canonical_all)) {
  canonical_all$Assembly_Base_ID <- tools::file_path_sans_ext(basename(canonical_all$full_path))
}
if (!"file_name" %in% names(canonical_all)) canonical_all$file_name <- basename(canonical_all$full_path)

canonical_all <- canonical_all %>%
  mutate(
    selected_canonical = normalise_bool(selected_canonical),
    QC_PASS = normalise_bool(QC_PASS),
    Batch = suppressWarnings(as.integer(Batch)),
    canonical_norm_id = normalise_sample_base(Assembly_Base_ID),
    canonical_assembler = str_to_lower(coalesce(
      if ("assembler" %in% names(.)) as.character(assembler) else NA_character_,
      if ("Assembler" %in% names(.)) as.character(Assembler) else NA_character_,
      detect_assembler(coalesce(as.character(full_path), as.character(Assembly_Base_ID)))
    )),
    Participant_id = as.character(Participant_id),
    tp_lab = as.character(tp_lab),
    Isolate_ID = as.character(Isolate_ID),
    full_path = normalizePath(as.character(full_path), winslash = "/", mustWork = TRUE),
    fasta_sha256 = vapply(full_path, digest::digest, character(1), algo = "sha256", file = TRUE)
  )

canonical <- canonical_all %>%
  filter(selected_canonical, QC_PASS)

if (nrow(canonical) == 0) {
  stop("No selected QC-passing canonical rows found in ", canonical_path)
}
if (any(is.na(canonical$canonical_assembler) | canonical$canonical_assembler != "longcycler")) {
  stop(
    "Selected canonical MLST denominator contains non-Longcycler or unknown assembler rows. ",
    "Rerun 12a_wgs_qc.R with Longcycler-only canonical selection."
  )
}
if (anyDuplicated(canonical[c("Participant_id", "tp_lab")])) {
  stop("Selected Longcycler canonical manifest has duplicate participant-timepoint rows.")
}
if (anyDuplicated(canonical$canonical_norm_id)) {
  stop("Selected Longcycler canonical manifest has duplicate normalized isolate IDs.")
}

canonical <- canonical %>%
  select(
    Batch, Participant_id, tp_lab, Isolate_ID, Assembly_ID, Assembly_Base_ID,
    canonical_norm_id, canonical_assembler, file_name, full_path, fasta_sha256
  )

if (!file.exists(local_mlst_path)) {
  stop("Longcycler local MLST provenance is missing: ", local_mlst_path)
}

local_mlst_raw <- read_csv(local_mlst_path, show_col_types = FALSE)
if (!"ST" %in% names(local_mlst_raw)) stop("Local MLST table lacks ST: ", local_mlst_path)
if (!"full_path" %in% names(local_mlst_raw)) stop("Local MLST table lacks full_path provenance: ", local_mlst_path)
if (!"file_name" %in% names(local_mlst_raw)) local_mlst_raw$file_name <- basename(local_mlst_raw$full_path)

local_mlst <- local_mlst_raw %>%
  mutate(
    local_norm_id = if ("Assembly_Base_ID" %in% names(.)) {
      normalise_sample_base(Assembly_Base_ID)
    } else if ("file_name" %in% names(.)) {
      normalise_sample_base(file_name)
    } else {
      normalise_sample_base(Isolate_ID)
    },
    local_assembler = str_to_lower(coalesce(
      if ("assembler" %in% names(.)) as.character(assembler) else NA_character_,
      if ("Assembler" %in% names(.)) as.character(Assembler) else NA_character_,
      detect_assembler(coalesce(as.character(full_path), as.character(file_name)))
    )),
    full_path = normalizePath(as.character(full_path), winslash = "/", mustWork = FALSE),
    local_ST = normalise_st(ST),
    local_scheme = if ("scheme" %in% names(.)) as.character(scheme) else NA_character_,
    local_mlst_complete = if ("mlst_complete" %in% names(.)) normalise_bool(mlst_complete) else NA
  )

if (any(is.na(local_mlst$local_assembler) | local_mlst$local_assembler != "longcycler")) {
  stop("Active local MLST provenance contains non-Longcycler or unknown assembler rows.")
}

local_mlst <- local_mlst %>%
  semi_join(canonical %>% select(full_path), by = "full_path")
if (nrow(local_mlst) != nrow(canonical) || !setequal(local_mlst$full_path, canonical$full_path)) {
  stop("Local MLST FASTA paths do not exactly match the selected Longcycler canonical manifest.")
}

local_mlst <- local_mlst %>%
  group_by(local_norm_id, local_assembler) %>%
  summarise(
    local_ST = first_nonmissing(local_ST),
    local_ST_set = collapse_set(local_ST),
    local_scheme = first_nonmissing(local_scheme),
    local_mlst_complete = any(local_mlst_complete %in% TRUE, na.rm = TRUE),
    .groups = "drop"
  )

canonical_local <- canonical %>%
  left_join(
    local_mlst,
    by = c("canonical_norm_id" = "local_norm_id", "canonical_assembler" = "local_assembler")
  )

if (nrow(canonical_local) != nrow(canonical)) {
  stop("Longcycler local-MLST join changed the canonical denominator.")
}

provider_raw_all <- provider_file_inventory %>%
  filter(exists) %>%
  select(provider_source, provider_scope, expected_batches, provider_role, provider_path) %>%
  pmap_dfr(read_provider_file)

provider_raw <- provider_raw_all %>%
  filter(provider_assembler == "longcycler")
if (nrow(provider_raw) == 0) stop("No explicitly source-labelled Longcycler provider MLST rows were found.")
if (any(is.na(provider_raw$provider_assembler) | provider_raw$provider_assembler != "longcycler")) {
  stop("Internal error: a non-Longcycler or missing source label survived provider filtering.")
}

canonical_lookup <- canonical %>%
  select(
    canonical_norm_id, Batch, Participant_id, tp_lab, Isolate_ID, Assembly_ID,
    Assembly_Base_ID, canonical_assembler, file_name, full_path, fasta_sha256
  )

provider_normalized <- provider_raw %>%
  inner_join(canonical_lookup, by = c("provider_norm_id" = "canonical_norm_id")) %>%
  left_join(
    canonical_local %>% select(canonical_norm_id, canonical_assembler, local_ST, local_scheme),
    by = c("provider_norm_id" = "canonical_norm_id")
  ) %>%
  mutate(
    assembler_matches_canonical = provider_assembler == "longcycler" &
      provider_assembler == canonical_assembler.x &
      provider_assembler == canonical_assembler.y,
    matched_canonical = !is.na(Isolate_ID) & assembler_matches_canonical,
    provider_ST_called = !is.na(provider_ST),
    expected_batch_match = case_when(
      !matched_canonical ~ NA,
      str_detect(expected_batches, paste0("(^|;)", Batch, "($|;)")) ~ TRUE,
      TRUE ~ FALSE
    ),
    canonical_assembler = coalesce(canonical_assembler.x, canonical_assembler.y),
    provider_match_basis = "explicit_source_longcycler_label_plus_unique_selected_manifest_normalized_key;current_manifest_path_sha256_attached_after_match",
    provider_scheme_note = "Provider rows require an explicit Longcycler source label and a unique selected-manifest key match. Current path/SHA-256 provenance is attached after that match and does not independently prove the bytes analysed externally. Provider ST is the classic seven-locus call emitted alongside SeqSphere wgMLST results."
  ) %>%
  select(
    provider_source, provider_scope, expected_batches, provider_role, provider_file, provider_path,
    provider_sample_id, provider_assembler, provider_norm_id, provider_ST,
    provider_CC, provider_Profile, provider_PercGoodTargets, provider_qc_ge_90,
    provider_qc_ge_95, provider_ST_called, provider_has_classic_7_loci,
    provider_delimiter, matched_canonical, assembler_matches_canonical, expected_batch_match, Batch,
    Participant_id, tp_lab, Isolate_ID, Assembly_ID, Assembly_Base_ID, file_name, full_path, fasta_sha256,
    canonical_assembler, local_ST, local_scheme, provider_scheme_note,
    provider_match_basis,
    provider_n_rows_in_file, provider_n_cols_in_file
  )

if (any(provider_normalized$provider_assembler != "longcycler")) {
  stop("Active provider-normalized MLST contains non-Longcycler provenance.")
}
if (anyDuplicated(provider_normalized$provider_norm_id)) {
  stop("Active provider-normalized MLST contains duplicate selected-manifest keys.")
}
if (
  nrow(provider_normalized) != nrow(canonical) ||
    !setequal(provider_normalized$provider_norm_id, canonical$canonical_norm_id)
) {
  stop("Active provider-normalized rows do not exactly cover the selected Longcycler manifest keys.")
}
matched_provider <- provider_normalized %>% filter(matched_canonical)
if (nrow(matched_provider) > 0 && any(matched_provider$canonical_assembler != "longcycler")) {
  stop("Matched provider MLST contains non-Longcycler canonical provenance.")
}

write_csv(provider_normalized, file.path(out_dir, "provider_mlst_normalized.csv"))

summarise_provider_by_id <- function(threshold = NA_real_, threshold_label = "none") {
  provider_normalized %>%
    mutate(
      provider_valid_for_threshold = provider_ST_called &
        (is.na(threshold) | provider_PercGoodTargets >= threshold)
    ) %>%
    group_by(provider_norm_id) %>%
    summarise(
      provider_qc_threshold = threshold_label,
      provider_ST_called = any(provider_valid_for_threshold, na.rm = TRUE),
      provider_ST = first_nonmissing(provider_ST[provider_valid_for_threshold]),
      provider_ST_set = collapse_set(provider_ST[provider_valid_for_threshold]),
      provider_ST_n = n_distinct(provider_ST[provider_valid_for_threshold], na.rm = TRUE),
      provider_sources = collapse_set(provider_source[provider_valid_for_threshold]),
      provider_files = collapse_set(provider_file[provider_valid_for_threshold]),
      provider_assemblers = collapse_set(provider_assembler[provider_valid_for_threshold]),
      provider_max_PercGoodTargets = suppressWarnings(max(provider_PercGoodTargets[provider_valid_for_threshold], na.rm = TRUE)),
      provider_min_PercGoodTargets = suppressWarnings(min(provider_PercGoodTargets[provider_valid_for_threshold], na.rm = TRUE)),
      n_provider_rows_for_threshold = sum(provider_valid_for_threshold, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      provider_ST_n = if_else(is.na(provider_ST), 0L, as.integer(provider_ST_n)),
      provider_max_PercGoodTargets = if_else(is.infinite(provider_max_PercGoodTargets), NA_real_, provider_max_PercGoodTargets),
      provider_min_PercGoodTargets = if_else(is.infinite(provider_min_PercGoodTargets), NA_real_, provider_min_PercGoodTargets)
    )
}

provider_by_id_all_thresholds <- bind_rows(
  summarise_provider_by_id(NA_real_, "none"),
  summarise_provider_by_id(90, ">=90"),
  summarise_provider_by_id(95, ">=95")
)

coverage_by_batch <- canonical_local %>%
  select(Batch, Participant_id, tp_lab, Isolate_ID, canonical_norm_id, local_ST) %>%
  tidyr::crossing(provider_qc_threshold = c("none", ">=90", ">=95")) %>%
  left_join(
    provider_by_id_all_thresholds,
    by = c("canonical_norm_id" = "provider_norm_id", "provider_qc_threshold")
  ) %>%
  mutate(provider_ST_called_flag = coalesce(provider_ST_called, FALSE)) %>%
  group_by(provider_qc_threshold, Batch) %>%
  summarise(
    n_canonical = n(),
    local_ST_called = sum(!is.na(local_ST)),
    provider_ST_called = sum(provider_ST_called_flag),
    both_called = sum(!is.na(local_ST) & provider_ST_called_flag),
    provider_only = sum(is.na(local_ST) & provider_ST_called_flag),
    local_only = sum(!is.na(local_ST) & !provider_ST_called_flag),
    neither = sum(is.na(local_ST) & !provider_ST_called_flag),
    local_coverage_pct = round(100 * local_ST_called / n_canonical, 1),
    provider_coverage_pct = round(100 * provider_ST_called / n_canonical, 1),
    provider_minus_local = provider_ST_called - local_ST_called,
    .groups = "drop"
  ) %>%
  arrange(match(provider_qc_threshold, c("none", ">=90", ">=95")), Batch)

write_csv(coverage_by_batch, file.path(out_dir, "mlst_coverage_by_batch.csv"))

coverage_by_qc_threshold <- coverage_by_batch %>%
  group_by(provider_qc_threshold) %>%
  summarise(
    n_canonical = sum(n_canonical),
    local_ST_called = sum(local_ST_called),
    provider_ST_called = sum(provider_ST_called),
    both_called = sum(both_called),
    provider_only = sum(provider_only),
    local_only = sum(local_only),
    neither = sum(neither),
    local_coverage_pct = round(100 * local_ST_called / n_canonical, 1),
    provider_coverage_pct = round(100 * provider_ST_called / n_canonical, 1),
    provider_minus_local = provider_ST_called - local_ST_called,
    .groups = "drop"
  ) %>%
  arrange(match(provider_qc_threshold, c("none", ">=90", ">=95")))

write_csv(coverage_by_qc_threshold, file.path(out_dir, "mlst_coverage_by_qc_threshold.csv"))

coverage_by_provider_source <- provider_normalized %>%
  filter(matched_canonical) %>%
  mutate(provider_ST_called_flag = provider_ST_called) %>%
  group_by(provider_source, provider_scope, provider_file) %>%
  summarise(
    n_provider_rows_matched = n(),
    n_canonical_ids_matched = n_distinct(provider_norm_id),
    provider_ST_called = n_distinct(provider_norm_id[provider_ST_called_flag]),
    provider_ST_qc90 = n_distinct(provider_norm_id[provider_ST_called_flag & provider_qc_ge_90]),
    provider_ST_qc95 = n_distinct(provider_norm_id[provider_ST_called_flag & provider_qc_ge_95]),
    local_ST_called_among_matched = n_distinct(provider_norm_id[!is.na(local_ST)]),
    provider_rescues_local_missing_qc95 = n_distinct(provider_norm_id[provider_ST_called_flag & provider_qc_ge_95 & is.na(local_ST)]),
    local_rescues_provider_missing_qc95 = n_distinct(provider_norm_id[!provider_ST_called_flag & !is.na(local_ST)]),
    .groups = "drop"
  )

write_csv(coverage_by_provider_source, file.path(out_dir, "mlst_coverage_by_provider_source.csv"))

default_provider_by_id <- provider_by_id_all_thresholds %>%
  filter(provider_qc_threshold == ">=95")

canonical_with_default_provider <- canonical_local %>%
  left_join(default_provider_by_id, by = c("canonical_norm_id" = "provider_norm_id")) %>%
  mutate(provider_ST_called = coalesce(provider_ST_called, FALSE))

mlst_rescue_cases <- canonical_with_default_provider %>%
  filter(is.na(local_ST), provider_ST_called) %>%
  transmute(
    Batch, Participant_id, tp_lab, Isolate_ID, Assembly_ID, Assembly_Base_ID,
    canonical_assembler, local_ST, provider_ST, provider_ST_set,
    provider_sources, provider_files, provider_assemblers,
    provider_min_PercGoodTargets, provider_max_PercGoodTargets,
    n_provider_rows_for_threshold
  ) %>%
  arrange(Batch, Participant_id, tp_lab, Isolate_ID)

write_csv(mlst_rescue_cases, file.path(out_dir, "mlst_rescue_cases.csv"))

local_only_cases <- canonical_with_default_provider %>%
  filter(!is.na(local_ST), !provider_ST_called) %>%
  transmute(
    Batch, Participant_id, tp_lab, Isolate_ID, Assembly_ID, Assembly_Base_ID,
    canonical_assembler, local_ST, local_scheme, provider_ST,
    provider_sources, provider_files
  ) %>%
  arrange(Batch, Participant_id, tp_lab, Isolate_ID)

write_csv(local_only_cases, file.path(out_dir, "mlst_local_only_cases.csv"))

neither_cases <- canonical_with_default_provider %>%
  filter(is.na(local_ST), !provider_ST_called) %>%
  transmute(
    Batch, Participant_id, tp_lab, Isolate_ID, Assembly_ID, Assembly_Base_ID,
    canonical_assembler, local_ST, provider_ST,
    provider_sources, provider_files
  ) %>%
  arrange(Batch, Participant_id, tp_lab, Isolate_ID)

write_csv(neither_cases, file.path(out_dir, "mlst_neither_source_cases.csv"))

provider_internal_conflicts <- provider_by_id_all_thresholds %>%
  filter(provider_ST_n > 1) %>%
  inner_join(
    canonical %>% select(canonical_norm_id, Batch, Participant_id, tp_lab, Isolate_ID),
    by = c("provider_norm_id" = "canonical_norm_id")
  ) %>%
  transmute(
    issue_type = "provider_internal_ST_conflict",
    severity = "review",
    provider_qc_threshold,
    Batch, Participant_id, tp_lab, Isolate_ID,
    norm_id = provider_norm_id,
    local_ST = NA_character_,
    provider_ST = provider_ST_set,
    provider_sources,
    provider_files,
    note = "Multiple provider rows/files assign different provider ST labels for the same normalized isolate."
  )

local_provider_numeric_differences <- canonical_with_default_provider %>%
  filter(!is.na(local_ST), provider_ST_called, local_ST != provider_ST) %>%
  transmute(
    issue_type = "local_provider_numeric_ST_difference",
    severity = "info",
    provider_qc_threshold = ">=95",
    Batch, Participant_id, tp_lab, Isolate_ID,
    norm_id = canonical_norm_id,
    local_ST,
    provider_ST = provider_ST_set,
    provider_sources,
    provider_files,
    note = "Numeric ST labels differ. This is not treated as an error because local and provider outputs may use different E. coli MLST schemes/nomenclature."
  )

mlst_source_discordance_audit <- bind_rows(
  provider_internal_conflicts,
  local_provider_numeric_differences
) %>%
  arrange(issue_type, Batch, Participant_id, tp_lab, Isolate_ID)

write_csv(mlst_source_discordance_audit, file.path(out_dir, "mlst_source_discordance_audit.csv"))

pair_context <- function(st_a, st_b) {
  st_a <- normalise_st(st_a)
  st_b <- normalise_st(st_b)
  case_when(
    is.na(st_a) | is.na(st_b) ~ "Missing ST evidence",
    st_a == st_b ~ "Same ST",
    TRUE ~ "Different ST"
  )
}

make_within_participant_pairs <- function(df) {
  df <- df %>%
    arrange(Participant_id, tp_lab, Isolate_ID) %>%
    group_by(Participant_id) %>%
    mutate(sample_index = row_number()) %>%
    ungroup()

  split(df, df$Participant_id) %>%
    map_dfr(function(x) {
      if (nrow(x) < 2) return(tibble())
      pair_idx <- utils::combn(seq_len(nrow(x)), 2)
      tibble(
        Participant_id = x$Participant_id[pair_idx[1, ]],
        Timepoint_A = x$tp_lab[pair_idx[1, ]],
        Isolate_ID_A = x$Isolate_ID[pair_idx[1, ]],
        norm_id_A = x$canonical_norm_id[pair_idx[1, ]],
        local_ST_A = x$local_ST[pair_idx[1, ]],
        provider_ST_A = x$provider_ST[pair_idx[1, ]],
        Timepoint_B = x$tp_lab[pair_idx[2, ]],
        Isolate_ID_B = x$Isolate_ID[pair_idx[2, ]],
        norm_id_B = x$canonical_norm_id[pair_idx[2, ]],
        local_ST_B = x$local_ST[pair_idx[2, ]],
        provider_ST_B = x$provider_ST[pair_idx[2, ]]
      )
    })
}

mlst_pair_context_comparison <- canonical_with_default_provider %>%
  transmute(
    Batch, Participant_id, tp_lab, Isolate_ID, canonical_norm_id,
    local_ST,
    provider_ST = if_else(provider_ST_called, provider_ST, NA_character_)
  ) %>%
  make_within_participant_pairs() %>%
  mutate(
    local_ST_context = pair_context(local_ST_A, local_ST_B),
    provider_ST_context = pair_context(provider_ST_A, provider_ST_B),
    provider_rescues_missing_local_context =
      local_ST_context == "Missing ST evidence" & provider_ST_context != "Missing ST evidence",
    local_rescues_missing_provider_context =
      provider_ST_context == "Missing ST evidence" & local_ST_context != "Missing ST evidence",
    changed_same_different_interpretation =
      local_ST_context != "Missing ST evidence" &
      provider_ST_context != "Missing ST evidence" &
      local_ST_context != provider_ST_context,
    scheme_caveat = "Compare same/different-ST context within source; numeric local/provider ST labels may not share nomenclature."
  ) %>%
  arrange(Participant_id, Timepoint_A, Timepoint_B, Isolate_ID_A, Isolate_ID_B)

write_csv(mlst_pair_context_comparison, file.path(out_dir, "mlst_pair_context_comparison.csv"))

pair_context_summary <- mlst_pair_context_comparison %>%
  summarise(
    n_pairs = n(),
    local_missing_context = sum(local_ST_context == "Missing ST evidence"),
    provider_missing_context = sum(provider_ST_context == "Missing ST evidence"),
    provider_rescues_missing_local_context = sum(provider_rescues_missing_local_context),
    local_rescues_missing_provider_context = sum(local_rescues_missing_provider_context),
    changed_same_different_interpretation = sum(changed_same_different_interpretation)
  )

write_csv(pair_context_summary, file.path(out_dir, "mlst_pair_context_summary.csv"))

default_overall <- coverage_by_qc_threshold %>%
  filter(provider_qc_threshold == ">=95")

default_batch <- coverage_by_batch %>%
  filter(provider_qc_threshold == ">=95") %>%
  mutate(
    line = sprintf(
      "| %s | %d | %d (%s) | %d (%s) | %+d | %d | %d | %d |",
      Batch,
      n_canonical,
      local_ST_called,
      sprintf("%.1f%%", local_coverage_pct),
      provider_ST_called,
      sprintf("%.1f%%", provider_coverage_pct),
      provider_minus_local,
      provider_only,
      local_only,
      neither
    )
  ) %>%
  pull(line)

threshold_lines <- coverage_by_qc_threshold %>%
  mutate(
    line = sprintf(
      "| %s | %d | %d (%s) | %d (%s) | %+d | %d | %d | %d |",
      provider_qc_threshold,
      n_canonical,
      local_ST_called,
      sprintf("%.1f%%", local_coverage_pct),
      provider_ST_called,
      sprintf("%.1f%%", provider_coverage_pct),
      provider_minus_local,
      provider_only,
      local_only,
      neither
    )
  ) %>%
  pull(line)

summary_lines <- c(
  "# MLST Source Coverage Comparison",
  "",
  sprintf("Generated: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  "",
  "## Main Result",
  "",
  sprintf(
    "Canonical selected-isolate denominator: **%d**.",
    default_overall$n_canonical
  ),
  sprintf(
    "Local `mlst` usable ST coverage: **%d/%d (%s)**.",
    default_overall$local_ST_called,
    default_overall$n_canonical,
    format_pct(default_overall$local_ST_called, default_overall$n_canonical)
  ),
  sprintf(
    "Provider usable ST coverage at `PercGoodTargets >= 95`: **%d/%d (%s)**.",
    default_overall$provider_ST_called,
    default_overall$n_canonical,
    format_pct(default_overall$provider_ST_called, default_overall$n_canonical)
  ),
  sprintf(
    "Provider adds **%d** more usable ST calls than local MLST at the default high-confidence threshold.",
    default_overall$provider_minus_local
  ),
  "",
  "## QC Threshold Sensitivity",
  "",
  "| Provider QC threshold | Canonical n | Local ST called | Provider ST called | Provider minus local | Provider only | Local only | Neither |",
  "|---|---:|---:|---:|---:|---:|---:|---:|",
  threshold_lines,
  "",
  "## Batch-Level Coverage at Provider `PercGoodTargets >= 95`",
  "",
  "| Batch | Canonical n | Local ST called | Provider ST called | Provider minus local | Provider only | Local only | Neither |",
  "|---|---:|---:|---:|---:|---:|---:|---:|",
  default_batch,
  "",
  "## Pair Context Summary",
  "",
  sprintf("Within-participant isolate pairs assessed: **%d**.", pair_context_summary$n_pairs),
  sprintf("Provider rescues missing local same/different-ST context for **%d** pairs.", pair_context_summary$provider_rescues_missing_local_context),
  sprintf("Local rescues missing provider same/different-ST context for **%d** pairs.", pair_context_summary$local_rescues_missing_provider_context),
  sprintf("Both sources have ST context but differ in same/different interpretation for **%d** pairs.", pair_context_summary$changed_same_different_interpretation),
  "",
  "## Caveats",
  "",
  "- Active coverage and integration use only explicitly source-labelled Longcycler provider rows with a unique selected-manifest key match. Current FASTA path/SHA-256 provenance is attached after matching; it does not independently prove the bytes analysed externally.",
  "- Provider and local ST labels are accepted into one active classic seven-locus layer only when all dual-usable calls agree; the integration gate fails on any discordance.",
  "- The primary coverage comparison is isolate-level: a usable ST is any non-missing ST after treating `\"\"`, `-`, `?`, and `NA`-like values as missing.",
  "- Provider calls with `PercGoodTargets >= 95` are treated as the default high-confidence provider calls.",
  "- This audit does not switch pipeline source files or overwrite `results/mlst/mlst_with_meta.csv`.",
  "",
  "## Output Files",
  "",
  "- `provider_mlst_normalized.csv`",
  "- `mlst_coverage_by_batch.csv`",
  "- `mlst_coverage_by_qc_threshold.csv`",
  "- `mlst_rescue_cases.csv`",
  "- `mlst_source_discordance_audit.csv`",
  "- `mlst_pair_context_comparison.csv`"
)

writeLines(summary_lines, file.path(out_dir, "summary.md"))

message("MLST source comparison complete.")
message("Outputs written to: ", out_dir)

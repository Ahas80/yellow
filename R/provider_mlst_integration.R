# ==============================================================================
# R/provider_mlst_integration.R
# ------------------------------------------------------------------------------
# Internal provider/RIVM MLST integration helper sourced by 06_MLST.R.
# It promotes Longcycler provider SeqSphere MLST to the active chromosomal ST
# source while preserving local MLST from the same Longcycler FASTA as a
# labelled fallback only.
# ==============================================================================

source("00_config.R")
source("R/pipeline_qc_helpers.R")

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
  library(tidyr)
})

provider_qc_threshold <- 95
canonical_selection_path <- file.path(DIR_QC, "canonical_assembly_selection.csv")

required_files <- c(
  FILE_MLST_LOCAL_CANONICAL,
  FILE_MLST_PROVIDER_NORMALIZED,
  canonical_selection_path
)
missing_required <- required_files[!file.exists(required_files)]
if (length(missing_required) > 0) {
  stop(
    "Missing required MLST integration input(s):\n",
    paste(missing_required, collapse = "\n"),
    "\nRun Rscript 06_MLST.R to refresh the active provider/RIVM MLST layer."
  )
}

missing_tokens <- c("", ".", "?", "-", "--", "na", "n/a", "nan", "missing", "unknown", "null")

is_missing_like <- function(x) {
  y <- str_squish(str_to_lower(as.character(x)))
  is.na(x) | is.na(y) | y %in% missing_tokens
}

boolish <- function(x, default = FALSE) {
  if (is.logical(x)) return(replace_na(x, default))
  y <- str_squish(str_to_lower(as.character(x)))
  dplyr::case_when(
    is.na(x) | y == "" ~ default,
    y %in% c("true", "t", "yes", "y", "1") ~ TRUE,
    y %in% c("false", "f", "no", "n", "0") ~ FALSE,
    TRUE ~ default
  )
}

read_mlst_table <- function(path) {
  if (grepl("\\.tsv$", path, ignore.case = TRUE)) {
    readr::read_tsv(path, show_col_types = FALSE, progress = FALSE)
  } else {
    readr::read_csv(path, show_col_types = FALSE, progress = FALSE)
  }
}

collapse_unique <- function(x) {
  vals <- sort(unique(na.omit(as.character(x))))
  vals <- vals[nzchar(vals)]
  if (length(vals) == 0) NA_character_ else paste(vals, collapse = ";")
}

prepare_local_for_preference <- function(df) {
  if (!"Isolate_ID" %in% names(df)) {
    stop("Local MLST table lacks Isolate_ID.")
  }
  if ("ST" %in% names(df)) {
    df <- df %>% rename(ST_local = ST)
  }
  if (!"ST_local" %in% names(df)) {
    df$ST_local <- NA_character_
  }
  if (!"full_path" %in% names(df)) {
    stop("Local MLST table lacks full_path provenance.")
  }
  if (!"file_name" %in% names(df)) df$file_name <- basename(df$full_path)
  df %>%
    mutate(
      Isolate_ID = as.character(Isolate_ID),
      full_path = normalizePath(as.character(full_path), winslash = "/", mustWork = FALSE),
      local_assembler = str_to_lower(coalesce(
        if ("assembler" %in% names(.)) as.character(assembler) else NA_character_,
        if ("Assembler" %in% names(.)) as.character(Assembler) else NA_character_,
        detect_assembler(coalesce(as.character(full_path), as.character(file_name)))
      )),
      ST_local = as.character(ST_local),
      local_ST_called = !is_missing_like(ST_local)
    )
}

apply_provider_preference <- function(df, provider_summary) {
  df %>%
    left_join(provider_summary, by = "Isolate_ID") %>%
    mutate(
      provider_ST_called = !is_missing_like(ST_provider),
      ST = case_when(
        provider_ST_called ~ ST_provider,
        !provider_ST_called & local_ST_called ~ ST_local,
        TRUE ~ NA_character_
      ),
      ST_source = case_when(
        provider_ST_called ~ "provider_qc95",
        provider_internal_conflict %in% TRUE & local_ST_called ~ "local_fallback_provider_conflict",
        provider_internal_conflict %in% TRUE ~ "missing_provider_conflict",
        local_ST_called ~ "local_fallback_provider_missing",
        TRUE ~ "missing"
      ),
      ST_scheme = case_when(
        ST_source == "provider_qc95" ~ "provider_seqsphere_wgMLST_Eco_BSR_MLST",
        str_starts(ST_source, "local_fallback") ~ "local_mlst",
        TRUE ~ NA_character_
      ),
      ST_numeric_comparable_to_local = FALSE
    ) %>%
    select(
      any_of(c(
        "Isolate_ID", "Participant_id", "tp_lab", "Timepoint", "Episode_ID", "Batch",
        "Assembly_ID", "Assembly_Base_ID", "Assembler", "assembler", "file_name", "full_path"
      )),
      ST, ST_source, ST_provider, ST_local, ST_scheme,
      ST_numeric_comparable_to_local, provider_PercGoodTargets,
      provider_file, provider_batch_match, provider_assembler, provider_source,
      provider_internal_conflict, provider_ST_values,
      everything()
    )
}

selection <- read_csv(canonical_selection_path, show_col_types = FALSE, progress = FALSE)
required_selection_cols <- c("Participant_id", "tp_lab", "Isolate_ID", "selected_canonical", "QC_PASS")
missing_selection_cols <- setdiff(required_selection_cols, names(selection))
if (length(missing_selection_cols) > 0) {
  stop("Canonical assembly selection lacks required column(s): ", paste(missing_selection_cols, collapse = ", "))
}
if (!"full_path" %in% names(selection) && "fasta_path" %in% names(selection)) {
  selection$full_path <- selection$fasta_path
}
if (!"full_path" %in% names(selection)) stop("Canonical assembly selection lacks full_path/fasta_path.")
if (!"file_name" %in% names(selection)) selection$file_name <- basename(selection$full_path)

selection <- selection %>%
  mutate(
    selected_canonical = boolish(selected_canonical),
    QC_PASS = boolish(QC_PASS),
    canonical_assembler = str_to_lower(coalesce(
      if ("assembler" %in% names(.)) as.character(assembler) else NA_character_,
      if ("Assembler" %in% names(.)) as.character(Assembler) else NA_character_,
      detect_assembler(coalesce(as.character(full_path), as.character(file_name)))
    )),
    full_path = normalizePath(as.character(full_path), winslash = "/", mustWork = FALSE)
  )

canonical_manifest <- selection %>%
  filter(selected_canonical %in% TRUE, QC_PASS %in% TRUE)
if (nrow(canonical_manifest) == 0) stop("Canonical selection contains no selected QC-passing assemblies.")
if (any(is.na(canonical_manifest$canonical_assembler) | canonical_manifest$canonical_assembler != "longcycler")) {
  stop("Canonical MLST denominator contains non-Longcycler or unknown assembler rows.")
}
if (anyDuplicated(canonical_manifest[c("Participant_id", "tp_lab")])) {
  stop("Canonical Longcycler MLST denominator has duplicate participant-timepoint rows.")
}
canonical_denominator <- nrow(canonical_manifest)
canonical_paths <- unique(canonical_manifest$full_path)

msg("Reading local Longcycler MLST provenance: %s", FILE_MLST_LOCAL_CANONICAL)
local_canonical <- read_csv(FILE_MLST_LOCAL_CANONICAL, show_col_types = FALSE, progress = FALSE) %>%
  prepare_local_for_preference() %>%
  filter(full_path %in% canonical_paths)

if (any(is.na(local_canonical$local_assembler) | local_canonical$local_assembler != "longcycler")) {
  stop("Local canonical MLST contains non-Longcycler or unknown assembler provenance.")
}
if (nrow(local_canonical) != canonical_denominator || !setequal(local_canonical$full_path, canonical_paths)) {
  stop(
    "Local canonical MLST does not exactly match the selected Longcycler manifest: local=",
    nrow(local_canonical), ", canonical=", canonical_denominator, "."
  )
}

if (anyDuplicated(local_canonical$Isolate_ID)) {
  dup_ids <- unique(local_canonical$Isolate_ID[duplicated(local_canonical$Isolate_ID)])
  stop("Local canonical MLST has duplicate Isolate_ID rows: ", paste(head(dup_ids, 10), collapse = ", "))
}

msg("Reading provider normalized MLST: %s", FILE_MLST_PROVIDER_NORMALIZED)
provider <- read_csv(FILE_MLST_PROVIDER_NORMALIZED, show_col_types = FALSE, progress = FALSE)
required_provider_cols <- c(
  "provider_assembler", "provider_norm_id", "provider_ST", "Isolate_ID",
  "matched_canonical", "assembler_matches_canonical"
)
missing_provider_cols <- setdiff(required_provider_cols, names(provider))
if (length(missing_provider_cols) > 0) {
  stop("Provider normalized MLST lacks required column(s): ", paste(missing_provider_cols, collapse = ", "))
}
provider <- provider %>%
  mutate(
    Isolate_ID = as.character(Isolate_ID),
    provider_assembler = str_to_lower(str_squish(as.character(provider_assembler))),
    provider_ST = as.character(provider_ST),
    provider_ST_called = if ("provider_ST_called" %in% names(.)) boolish(provider_ST_called) else !is_missing_like(provider_ST),
    provider_qc_ge_95 = if ("provider_qc_ge_95" %in% names(.)) boolish(provider_qc_ge_95) else provider_PercGoodTargets >= provider_qc_threshold,
    matched_canonical = if ("matched_canonical" %in% names(.)) boolish(matched_canonical) else !is.na(Isolate_ID),
    assembler_matches_canonical = boolish(assembler_matches_canonical),
    expected_batch_match = if ("expected_batch_match" %in% names(.)) boolish(expected_batch_match, default = NA) else NA
  )

if (any(is.na(provider$provider_assembler) | provider$provider_assembler != "longcycler")) {
  stop("Provider-normalized active input contains Flye, combined, or unknown assembler provenance.")
}
if (any(provider$matched_canonical %in% TRUE & !(provider$assembler_matches_canonical %in% TRUE))) {
  stop("Provider-normalized input marks an assembler-mismatched row as canonical.")
}

provider_qc95 <- provider %>%
  filter(
    provider_assembler == "longcycler",
    matched_canonical %in% TRUE,
    provider_qc_ge_95 %in% TRUE,
    provider_ST_called %in% TRUE,
    !is_missing_like(provider_ST)
  )

provider_conflicts <- provider_qc95 %>%
  group_by(Isolate_ID) %>%
  summarise(
    n_provider_qc95_rows = n(),
    n_provider_qc95_ST = n_distinct(provider_ST),
    provider_ST_values = collapse_unique(provider_ST),
    provider_files = collapse_unique(provider_file),
    provider_sources = collapse_unique(provider_source),
    .groups = "drop"
  ) %>%
  filter(n_provider_qc95_ST > 1)

provider_summary <- provider_qc95 %>%
  group_by(Isolate_ID) %>%
  summarise(
    n_provider_qc95_rows = n(),
    n_provider_qc95_ST = n_distinct(provider_ST),
    provider_ST_values = collapse_unique(provider_ST),
    ST_provider = if (n_distinct(provider_ST) == 1) first(provider_ST) else NA_character_,
    provider_PercGoodTargets = suppressWarnings(max(as.numeric(provider_PercGoodTargets), na.rm = TRUE)),
    provider_file = collapse_unique(provider_file),
    provider_source = collapse_unique(provider_source),
    provider_assembler = collapse_unique(provider_assembler),
    provider_batch_match = any(expected_batch_match %in% TRUE, na.rm = TRUE),
    provider_internal_conflict = n_distinct(provider_ST) > 1,
    .groups = "drop"
  ) %>%
  mutate(provider_PercGoodTargets = ifelse(is.infinite(provider_PercGoodTargets), NA_real_, provider_PercGoodTargets))

preferred_canonical <- apply_provider_preference(local_canonical, provider_summary)

if (nrow(preferred_canonical) != canonical_denominator) {
  stop(
    "Provider-preferred MLST denominator mismatch: wrote ", nrow(preferred_canonical),
    " row(s), but canonical selection has ", canonical_denominator, " selected row(s)."
  )
}

provider_primary_n <- sum(preferred_canonical$ST_source == "provider_qc95", na.rm = TRUE)
local_called_n <- sum(preferred_canonical$local_ST_called %in% TRUE, na.rm = TRUE)
preferred_called_n <- sum(!is_missing_like(preferred_canonical$ST), na.rm = TRUE)
local_fallback_n <- sum(str_starts(preferred_canonical$ST_source, "local_fallback"), na.rm = TRUE)
missing_n <- sum(preferred_canonical$ST_source %in% c("missing", "missing_provider_conflict"), na.rm = TRUE)

active_assembler <- str_to_lower(coalesce(
  if ("assembler" %in% names(preferred_canonical)) as.character(preferred_canonical$assembler) else NA_character_,
  if ("Assembler" %in% names(preferred_canonical)) as.character(preferred_canonical$Assembler) else NA_character_,
  preferred_canonical$local_assembler
))
if (any(is.na(active_assembler) | active_assembler != "longcycler")) {
  stop("Provider-preferred MLST contains non-Longcycler or missing active assembly provenance.")
}
provider_primary <- preferred_canonical %>% filter(ST_source == "provider_qc95")
if (nrow(provider_primary) > 0 && any(
  is.na(provider_primary$provider_assembler) |
    provider_primary$provider_assembler != "longcycler"
)) {
  stop("Provider-primary MLST contains Flye, combined, or missing provider provenance.")
}

if (provider_primary_n < local_called_n) {
  stop(
    "Provider QC95 MLST coverage (", provider_primary_n, ") is below local MLST coverage (",
    local_called_n, "). Review provider inputs before promoting this source."
  )
}

msg("Writing provider-preferred canonical MLST: %s", FILE_MLST_PROVIDER_PREFERRED)
write_csv(preferred_canonical, FILE_MLST_PROVIDER_PREFERRED)

if (file.exists(FILE_MLST_LOCAL_ALL)) {
  local_all <- read_mlst_table(FILE_MLST_LOCAL_ALL) %>%
    prepare_local_for_preference() %>%
    filter(full_path %in% canonical_paths)
  if (any(is.na(local_all$local_assembler) | local_all$local_assembler != "longcycler")) {
    stop("Local all/isolate MLST input contains non-Longcycler provenance.")
  }
  if (nrow(local_all) != canonical_denominator || !setequal(local_all$full_path, canonical_paths)) {
    stop("Local all/isolate MLST input does not exactly match the selected Longcycler manifest.")
  }
  preferred_all <- apply_provider_preference(local_all, provider_summary)
} else {
  preferred_all <- preferred_canonical
}

msg("Writing provider-preferred all/isolate MLST: %s", FILE_MLST_PROVIDER_PREFERRED_ALL)
write_csv(preferred_all, FILE_MLST_PROVIDER_PREFERRED_ALL)

coverage_audit <- tibble(
  audit_type = "coverage_summary",
  metric = c(
    "canonical_denominator",
    "local_usable_ST",
    "provider_qc95_usable_ST",
    "provider_preferred_usable_ST",
    "local_fallback_ST",
    "missing_ST"
  ),
  value = c(
    canonical_denominator,
    local_called_n,
    provider_primary_n,
    preferred_called_n,
    local_fallback_n,
    missing_n
  ),
  note = c(
    "Selected canonical isolates",
    "Usable local mlst ST before provider promotion",
    "Provider SeqSphere ST with PercGoodTargets >= 95",
    "Active ST after provider-primary plus labelled local fallback",
    "Rows where local ST is used because provider QC95 ST is absent or conflicted",
    "Rows where neither usable provider QC95 nor local ST is available"
  )
)

assignment_audit <- preferred_canonical %>%
  count(ST_source, name = "n") %>%
  transmute(
    audit_type = "source_assignment",
    metric = ST_source,
    value = n,
    note = "Provider-preferred canonical MLST source assignment"
  )

conflict_audit <- provider_conflicts %>%
  transmute(
    audit_type = "provider_qc95_conflict",
    metric = Isolate_ID,
    value = n_provider_qc95_ST,
    note = paste0(
      "Conflicting provider QC95 ST values: ", provider_ST_values,
      "; files: ", provider_files,
      "; sources: ", provider_sources
    )
  )

fallback_audit <- preferred_canonical %>%
  filter(str_starts(ST_source, "local_fallback") | ST_source %in% c("missing", "missing_provider_conflict")) %>%
  transmute(
    audit_type = "non_provider_primary_row",
    metric = Isolate_ID,
    value = NA_real_,
    note = paste0(
      ST_source,
      "; local_ST=", coalesce(ST_local, "NA"),
      "; provider_ST=", coalesce(ST_provider, "NA")
    )
  )

source_audit <- bind_rows(coverage_audit, assignment_audit, conflict_audit, fallback_audit)

msg("Writing MLST source audit: %s", FILE_MLST_PROVIDER_SOURCE_AUDIT)
write_csv(source_audit, FILE_MLST_PROVIDER_SOURCE_AUDIT)

msg(
  "Provider-preferred MLST complete: provider_qc95=%d, local_fallback=%d, missing=%d, active_usable=%d/%d.",
  provider_primary_n,
  local_fallback_n,
  missing_n,
  preferred_called_n,
  canonical_denominator
)

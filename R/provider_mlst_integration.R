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
canonical_selection_path <- FILE_ANALYSIS_ASSEMBLY_MANIFEST

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
  required_local_qc <- c("mlst_complete", "ambiguous_call")
  missing_local_qc <- setdiff(required_local_qc, names(df))
  if (length(missing_local_qc) > 0L) {
    stop("Local MLST table lacks QC column(s): ", paste(missing_local_qc, collapse = ", "))
  }
  if (!"file_name" %in% names(df)) df$file_name <- basename(df$full_path)
  df %>%
    mutate(
      Isolate_ID = as.character(Isolate_ID),
      full_path = normalizePath(as.character(full_path), winslash = "/", mustWork = FALSE),
      local_assembler = str_to_lower(coalesce(
        if ("assembler" %in% names(.)) as.character(assembler) else NA_character_,
        if ("Assembler" %in% names(.)) as.character(Assembler) else NA_character_
      )),
      ST_local = as.character(ST_local),
      local_mlst_complete = boolish(mlst_complete),
      local_ambiguous_call = boolish(ambiguous_call),
      local_ST_called = !is_missing_like(ST_local) &
        local_mlst_complete %in% TRUE &
        !(local_ambiguous_call %in% TRUE)
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
      ST_numeric_comparable_to_local = case_when(
        provider_ST_called & provider_has_classic_7_loci %in% TRUE ~ TRUE,
        !provider_ST_called & local_ST_called ~ TRUE,
        TRUE ~ NA
      )
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

selection <- load_analysis_assemblies(canonical_selection_path, require_files = TRUE)
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
      if ("Assembler" %in% names(.)) as.character(Assembler) else NA_character_
    )),
    full_path = normalizePath(as.character(full_path), winslash = "/", mustWork = TRUE),
    fasta_sha256 = vapply(full_path, digest::digest, character(1), algo = "sha256", file = TRUE)
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
canonical_fasta_hashes <- canonical_manifest$fasta_sha256

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
if (!"fasta_sha256" %in% names(local_canonical)) stop("Local canonical MLST lacks FASTA SHA-256 provenance.")
canonical_path_hash <- paste(canonical_manifest$full_path, canonical_manifest$fasta_sha256, sep = "\n")
local_path_hash <- paste(local_canonical$full_path, local_canonical$fasta_sha256, sep = "\n")
if (!setequal(local_path_hash, canonical_path_hash)) stop("Local canonical MLST FASTA path/hash does not exactly match the selected manifest.")

if (anyDuplicated(local_canonical$Isolate_ID)) {
  dup_ids <- unique(local_canonical$Isolate_ID[duplicated(local_canonical$Isolate_ID)])
  stop("Local canonical MLST has duplicate Isolate_ID rows: ", paste(head(dup_ids, 10), collapse = ", "))
}

msg("Reading provider normalized MLST: %s", FILE_MLST_PROVIDER_NORMALIZED)
provider <- read_csv(FILE_MLST_PROVIDER_NORMALIZED, show_col_types = FALSE, progress = FALSE)
required_provider_cols <- c(
  "provider_assembler", "provider_norm_id", "provider_ST", "Isolate_ID",
  "matched_canonical", "assembler_matches_canonical", "full_path", "fasta_sha256",
  "provider_has_classic_7_loci"
)
missing_provider_cols <- setdiff(required_provider_cols, names(provider))
if (length(missing_provider_cols) > 0) {
  stop("Provider normalized MLST lacks required column(s): ", paste(missing_provider_cols, collapse = ", "))
}
provider <- provider %>%
  mutate(
    Isolate_ID = as.character(Isolate_ID),
    provider_assembler = str_to_lower(str_squish(as.character(provider_assembler))),
    full_path = normalizePath(as.character(full_path), winslash = "/", mustWork = FALSE),
    fasta_sha256 = as.character(fasta_sha256),
    provider_ST = as.character(provider_ST),
    provider_has_classic_7_loci = boolish(provider_has_classic_7_loci),
    provider_ST_called = if ("provider_ST_called" %in% names(.)) boolish(provider_ST_called) else !is_missing_like(provider_ST),
    provider_qc_ge_95 = if ("provider_qc_ge_95" %in% names(.)) boolish(provider_qc_ge_95) else provider_PercGoodTargets >= provider_qc_threshold,
    matched_canonical = if ("matched_canonical" %in% names(.)) boolish(matched_canonical) else !is.na(Isolate_ID),
    assembler_matches_canonical = boolish(assembler_matches_canonical),
    expected_batch_match = if ("expected_batch_match" %in% names(.)) boolish(expected_batch_match, default = NA) else NA
  )

if (any(provider$matched_canonical %in% TRUE & !(provider$assembler_matches_canonical %in% TRUE))) {
  stop("Provider-normalized input marks an assembler-mismatched row as canonical.")
}
provider_path_hash <- paste(provider$full_path, provider$fasta_sha256, sep = "\n")
if (any(!provider_path_hash %in% canonical_path_hash)) {
  stop("Provider-normalized input contains provenance not tied to the current selected Longcycler FASTA path/hash.")
}

provider_qc95 <- provider %>%
  filter(
    matched_canonical %in% TRUE,
    provider_qc_ge_95 %in% TRUE,
    provider_ST_called %in% TRUE,
    !is_missing_like(provider_ST)
  )

if (nrow(provider_qc95) > 0 && any(!(provider_qc95$provider_has_classic_7_loci %in% TRUE))) {
  stop("A provider QC95 MLST call lacks evidence for the classic seven-locus scheme.")
}

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

provider_below_qc95 <- provider %>%
  filter(
    matched_canonical %in% TRUE,
    !(provider_qc_ge_95 %in% TRUE),
    provider_ST_called %in% TRUE,
    !is_missing_like(provider_ST)
  ) %>%
  group_by(Isolate_ID) %>%
  summarise(
    provider_below_qc95_n_rows = n(),
    provider_below_qc95_n_ST = n_distinct(provider_ST),
    provider_below_qc95_ST_values = collapse_unique(provider_ST),
    ST_provider_below_qc95 = if (n_distinct(provider_ST) == 1) first(provider_ST) else NA_character_,
    provider_below_qc95_PercGoodTargets = suppressWarnings(max(as.numeric(provider_PercGoodTargets), na.rm = TRUE)),
    provider_below_qc95_file = collapse_unique(provider_file),
    provider_below_qc95_has_classic_7_loci = all(provider_has_classic_7_loci %in% TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    provider_below_qc95_PercGoodTargets = ifelse(
      is.infinite(provider_below_qc95_PercGoodTargets),
      NA_real_,
      provider_below_qc95_PercGoodTargets
    )
  )

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
    provider_has_classic_7_loci = all(provider_has_classic_7_loci %in% TRUE),
    provider_batch_match = any(expected_batch_match %in% TRUE, na.rm = TRUE),
    provider_internal_conflict = n_distinct(provider_ST) > 1,
    .groups = "drop"
  ) %>%
  mutate(provider_PercGoodTargets = ifelse(is.infinite(provider_PercGoodTargets), NA_real_, provider_PercGoodTargets))

preferred_canonical <- apply_provider_preference(local_canonical, provider_summary) %>%
  left_join(provider_below_qc95, by = "Isolate_ID")

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
below_qc95_evidence_n <- sum(!is_missing_like(preferred_canonical$ST_provider_below_qc95), na.rm = TRUE)
fallback_with_below_qc95_n <- sum(
  str_starts(preferred_canonical$ST_source, "local_fallback") &
    !is_missing_like(preferred_canonical$ST_provider_below_qc95),
  na.rm = TRUE
)
missing_with_below_qc95_n <- sum(
  preferred_canonical$ST_source %in% c("missing", "missing_provider_conflict") &
    !is_missing_like(preferred_canonical$ST_provider_below_qc95),
  na.rm = TRUE
)
fallback_below_qc95 <- preferred_canonical %>%
  filter(
    str_starts(.data$ST_source, "local_fallback"),
    !is_missing_like(.data$ST_provider_below_qc95)
  )
fallback_below_qc95_discordant_n <- sum(
  as.character(fallback_below_qc95$ST_local) !=
    as.character(fallback_below_qc95$ST_provider_below_qc95),
  na.rm = TRUE
)
dual_usable <- preferred_canonical %>%
  filter(.data$ST_source == "provider_qc95", .data$local_ST_called %in% TRUE)
dual_usable_n <- nrow(dual_usable)
dual_discordant_n <- sum(
  as.character(dual_usable$ST_provider) != as.character(dual_usable$ST_local),
  na.rm = TRUE
)

active_assembler <- str_to_lower(coalesce(
  if ("assembler" %in% names(preferred_canonical)) as.character(preferred_canonical$assembler) else NA_character_,
  if ("Assembler" %in% names(preferred_canonical)) as.character(preferred_canonical$Assembler) else NA_character_,
  preferred_canonical$local_assembler
))
if (any(is.na(active_assembler) | active_assembler != "longcycler")) {
  stop("Provider-preferred MLST contains non-Longcycler or missing active assembly provenance.")
}
provider_primary <- preferred_canonical %>% filter(ST_source == "provider_qc95")
if (nrow(provider_primary) > 0 && any(provider_primary$provider_assembler != "longcycler")) {
  stop("Provider-primary MLST is not tied to the selected Longcycler manifest.")
}
if (nrow(provider_primary) > 0 && any(!(provider_primary$provider_has_classic_7_loci %in% TRUE))) {
  stop("Provider-primary MLST lacks classic seven-locus scheme evidence.")
}
if (dual_discordant_n > 0L) {
  stop(
    "Provider and local classic seven-locus ST disagree for ",
    dual_discordant_n,
    " dual-usable selected Longcycler isolate(s)."
  )
}
if (fallback_below_qc95_discordant_n > 0L) {
  stop(
    "Local fallback and retained below-QC95 provider ST disagree for ",
    fallback_below_qc95_discordant_n,
    " selected Longcycler isolate(s)."
  )
}

if (provider_primary_n == 0L) {
  stop("Provider-preferred MLST has no usable provider QC95 calls.")
}
if (preferred_called_n != provider_primary_n + local_fallback_n) {
  stop(
    "Provider-preferred usable-ST count does not equal provider-primary plus labelled local fallback."
  )
}
if (provider_primary_n + local_fallback_n + missing_n != canonical_denominator) {
  stop("Provider-preferred MLST source assignments do not sum to the canonical denominator.")
}
if (preferred_called_n < local_called_n) {
  stop(
    "Provider-preferred MLST coverage after labelled Longcycler fallback (",
    preferred_called_n, ") is below local MLST coverage (", local_called_n,
    "). Review provider inputs before promoting this source."
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
    "missing_ST",
    "below_qc95_provider_ST_evidence",
    "local_fallback_with_below_qc95_provider_ST",
    "missing_with_below_qc95_provider_ST",
    "local_fallback_below_qc95_provider_discordant_ST",
    "provider_local_dual_usable_ST",
    "provider_local_dual_discordant_ST"
  ),
  value = c(
    canonical_denominator,
    local_called_n,
    provider_primary_n,
    preferred_called_n,
    local_fallback_n,
    missing_n,
    below_qc95_evidence_n,
    fallback_with_below_qc95_n,
    missing_with_below_qc95_n,
    fallback_below_qc95_discordant_n,
    dual_usable_n,
    dual_discordant_n
  ),
  note = c(
    "Selected canonical isolates",
    "Usable local mlst ST before provider promotion",
    "Provider SeqSphere ST with PercGoodTargets >= 95",
    "Active ST after provider-primary plus labelled local fallback",
    "Rows where local ST is used because provider QC95 ST is absent or conflicted",
    "Rows where neither usable provider QC95 nor local ST is available",
    "Rows retaining a provider ST below the QC95 promotion threshold",
    "Local-fallback rows that retain concordant below-QC95 provider ST evidence",
    "Missing active-ST rows that retain below-QC95 provider ST evidence",
    "Local-fallback rows whose retained below-QC95 provider ST disagrees with the local ST",
    "Rows with both provider QC95 and usable local classic seven-locus ST",
    "Dual-usable rows whose provider and local classic seven-locus ST disagree"
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
      "; provider_qc95_ST=", coalesce(ST_provider, "NA"),
      "; provider_below_qc95_ST=", coalesce(ST_provider_below_qc95, "NA"),
      "; provider_below_qc95_good_targets=", coalesce(
        as.character(provider_below_qc95_PercGoodTargets),
        "NA"
      )
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

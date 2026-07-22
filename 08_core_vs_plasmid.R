#!/usr/bin/env Rscript
# ==============================================================================
# 08_core_vs_plasmid.R
# ==============================================================================
#
# Descriptive lineage context for the canonical gene-level PlasmidFinder calls
# produced by script 09. This script never calls ABRicate or pMLST and does not
# treat repeated episodes as independent evidence for plasmid-ST association.
# ==============================================================================

source("00_config.R")
suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
  library(scales)
  library(tibble)
})

EXPECTED_EPISODES <- 532L
pf_hits_path <- file.path(DIR_PLASMIDS, "plasmidfinder_hits_long.csv")
pf_pa_path <- file.path(DIR_PLASMIDS, "plasmidfinder_presence_absence.csv")
pf_marker <- file.path(DIR_PLASMIDS, "PLASMIDFINDER_RUN_COMPLETE.txt")

required <- c(FILE_MLST_CANONICAL, pf_hits_path, pf_pa_path, pf_marker)
if (!all(file.exists(required))) {
  stop(
    "Script 08 requires completed canonical outputs from script 09: ",
    paste(basename(required[!file.exists(required)]), collapse = ", "),
    call. = FALSE
  )
}

manifest <- load_analysis_assemblies(
  FILE_ANALYSIS_ASSEMBLY_MANIFEST, require_files = TRUE
) %>%
  mutate(
    Isolate_ID = as.character(Isolate_ID),
    Participant_id = as.character(Participant_id),
    fasta_path = normalizePath(full_path, winslash = "/", mustWork = TRUE),
    fasta_sha256 = vapply(
      fasta_path, digest::digest, character(1),
      algo = "sha256", file = TRUE, serialize = FALSE
    )
  ) %>%
  distinct(Isolate_ID, .keep_all = TRUE)
if (nrow(manifest) != EXPECTED_EPISODES || anyDuplicated(manifest$Isolate_ID)) {
  stop("Script 08 requires the exact 532 selected Longcycler episodes.", call. = FALSE)
}

mlst <- read_csv(FILE_MLST_CANONICAL, show_col_types = FALSE) %>%
  mutate(
    Isolate_ID = as.character(Isolate_ID),
    ST = as.character(ST),
    ST_source = if ("ST_source" %in% names(.)) as.character(ST_source) else NA_character_
  )
if (nrow(mlst) != EXPECTED_EPISODES ||
    anyDuplicated(mlst$Isolate_ID) ||
    !setequal(mlst$Isolate_ID, manifest$Isolate_ID)) {
  stop("Canonical MLST rows do not match the exact 532-episode cohort.", call. = FALSE)
}

mlst %>%
  filter(!is.na(ST), nzchar(ST)) %>%
  count(ST, ST_source, sort = TRUE) %>%
  mutate(pct = percent(n / sum(n))) %>%
  write_csv(file.path(DIR_MLST, "ST_core_freq.csv"))

hits <- read_csv(pf_hits_path, show_col_types = FALSE) %>%
  mutate(
    Isolate_ID = as.character(Isolate_ID),
    fasta_path = normalizePath(fasta_path, winslash = "/", mustWork = TRUE),
    fasta_sha256 = tolower(as.character(fasta_sha256))
  )
required_hit_cols <- c(
  "Isolate_ID", "fasta_path", "fasta_sha256", "GENE",
  "accession", "identity", "coverage"
)
if (!all(required_hit_cols %in% names(hits))) {
  stop("Canonical PlasmidFinder long table lacks required columns.", call. = FALSE)
}
unmatched_hits <- hits %>%
  anti_join(
    manifest %>% select(Isolate_ID, fasta_path, fasta_sha256),
    by = c("Isolate_ID", "fasta_path", "fasta_sha256")
  )
if (nrow(unmatched_hits)) {
  stop("Canonical replicon hits do not match selected FASTA path/hash keys.", call. = FALSE)
}

pa <- read_csv(pf_pa_path, show_col_types = FALSE) %>%
  mutate(Isolate_ID = as.character(Isolate_ID))
if (nrow(pa) != EXPECTED_EPISODES ||
    anyDuplicated(pa$Isolate_ID) ||
    !setequal(pa$Isolate_ID, manifest$Isolate_ID)) {
  stop("Canonical gene-level plasmid matrix does not match all 532 episodes.", call. = FALSE)
}
replicon_cols <- setdiff(names(pa), "Isolate_ID")
if (!length(replicon_cols) ||
    any(replicon_cols %in% unique(hits$accession)) ||
    any(!vapply(pa[replicon_cols], function(x) all(x %in% c(0, 1)), logical(1)))) {
  stop("Compatibility exports require a binary GENE-label matrix.", call. = FALSE)
}

compat_long <- hits %>%
  transmute(
    isolate_id = Isolate_ID,
    fasta_path,
    fasta_sha256,
    replicon = GENE,
    accession,
    identity,
    coverage,
    contig_id = SEQUENCE
  )
write_csv(compat_long, file.path(DIR_MLST, "plasmid_replicons_long.csv"))

compat_wide <- pa %>%
  arrange(match(Isolate_ID, manifest$Isolate_ID))
write_csv(compat_wide, file.path(DIR_MLST, "plasmid_replicons_wide.csv"))
write_csv(compat_wide, file.path(DIR_MLST, "plasmid_types_per_isolate.csv"))

episode_long <- pa %>%
  pivot_longer(
    cols = all_of(replicon_cols),
    names_to = "replicon",
    values_to = "present"
  ) %>%
  left_join(
    manifest %>% select(Isolate_ID, Participant_id),
    by = "Isolate_ID"
  ) %>%
  left_join(
    mlst %>% select(Isolate_ID, ST, ST_source),
    by = "Isolate_ID"
  )

episode_prevalence <- episode_long %>%
  filter(!is.na(ST), nzchar(ST)) %>%
  group_by(ST, ST_source, replicon) %>%
  summarise(
    n_ST_episodes = n(),
    n_with_replicon = sum(present == 1L),
    prevalence = n_with_replicon / n_ST_episodes,
    .groups = "drop"
  ) %>%
  arrange(desc(n_ST_episodes), ST, desc(prevalence), replicon)
write_csv(
  episode_prevalence,
  file.path(DIR_MLST, "ST_replicon_prevalence_by_episode.csv")
)

resident_profile <- episode_long %>%
  filter(!is.na(ST), nzchar(ST)) %>%
  group_by(Participant_id, ST, ST_source, replicon) %>%
  summarise(present = as.integer(any(present == 1L)), .groups = "drop")
resident_prevalence <- resident_profile %>%
  group_by(ST, ST_source, replicon) %>%
  summarise(
    n_ST_residents = n(),
    n_residents_with_replicon = sum(present == 1L),
    prevalence = n_residents_with_replicon / n_ST_residents,
    .groups = "drop"
  ) %>%
  arrange(desc(n_ST_residents), ST, desc(prevalence), replicon)
write_csv(
  resident_prevalence,
  file.path(DIR_MLST, "ST_replicon_prevalence_by_resident.csv")
)

legacy_assoc <- file.path(DIR_MLST, "ST_plasmid_associations.csv")
if (file.exists(legacy_assoc)) {
  old <- tryCatch(read_csv(legacy_assoc, show_col_types = FALSE), error = function(e) NULL)
  if (!is.null(old) && any(c("p_value", "FDR", "OR") %in% names(old))) {
    archive_dir <- file.path(DIR_MLST, "superseded")
    dir.create(archive_dir, recursive = TRUE, showWarnings = FALSE)
    content_hash <- digest::digest(
      legacy_assoc, algo = "sha256", file = TRUE, serialize = FALSE
    )
    archive_path <- file.path(
      archive_dir,
      paste0("ST_plasmid_associations_repeated_episode_fisher_",
             substr(content_hash, 1L, 16L), ".csv")
    )
    if (!file.exists(archive_path) && !file.copy(legacy_assoc, archive_path)) {
      stop("Could not archive the superseded ST-plasmid association table.", call. = FALSE)
    }
    unlink(legacy_assoc)
  }
}
writeLines(
  c(
    "ST_plasmid_associations.csv is retired.",
    "Reason: the former Fisher tests omitted plasmid-negative episodes and treated repeated episodes as independent.",
    "Use ST_replicon_prevalence_by_episode.csv and ST_replicon_prevalence_by_resident.csv for descriptive lineage context.",
    "No inferential ST-plasmid association claim is made by script 08."
  ),
  file.path(DIR_MLST, "ST_plasmid_associations_RETIRED.txt")
)

writeLines(
  c(
    "Core-versus-plasmid descriptive layer: PASS",
    paste0("episodes=", nrow(pa)),
    paste0("replicon_gene_features=", length(replicon_cols)),
    "episode_denominator_includes_successful_no_hit_profiles=TRUE",
    "inference=descriptive only; repeated episodes are not treated as independent"
  ),
  file.path(DIR_MLST, "CORE_VS_PLASMID_RUN_COMPLETE.txt")
)
msg("✓ Descriptive ST–replicon context complete.")

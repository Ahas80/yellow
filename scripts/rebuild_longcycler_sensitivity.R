#!/usr/bin/env Rscript
# ==============================================================================
# Rebuild the canonical selected-Longcycler transition export
# ==============================================================================
# The primary genomic workflow is now selected QC-passing Longcycler only.
# This script independently rebuilds the episode and adjacent-transition tables
# from the canonical analysis cohort and verifies direct pairwise provenance.
# The historical filename is retained only so existing launchers keep working.

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
  library(tidyr)
})

source("00_config.R")

outdir <- file.path(DIR_RESULTS, "longitudinal")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

as_flag <- function(x) {
  if (is.logical(x)) return(replace(x, is.na(x), FALSE))
  tolower(trimws(as.character(x))) %in% c("true", "t", "1", "yes", "y")
}

parse_date_safe <- function(x) {
  x <- as.character(x)
  ans <- suppressWarnings(as.Date(x, format = "%d/%m/%Y"))
  missing <- is.na(ans)
  if (any(missing)) ans[missing] <- suppressWarnings(as.Date(x[missing], format = "%Y-%m-%d"))
  ans
}

fallback_time_order <- function(tp) {
  tp <- as.character(tp)
  upper <- str_to_upper(tp)
  case_when(
    str_detect(upper, "^T\\d+$") ~ suppressWarnings(as.numeric(str_remove(upper, "^T"))),
    str_detect(upper, "^UTI-\\d+$") ~ 100 + suppressWarnings(as.numeric(str_remove(upper, "^UTI-"))),
    str_detect(tp, regex("uricult", ignore_case = TRUE)) ~ 99,
    TRUE ~ NA_real_
  )
}

normalise_st <- function(x) {
  x <- str_trim(as.character(x))
  x[str_to_upper(x) %in% c("", "-", "ST-", "NA", "N/A", "UNKNOWN", "UNK", "NT",
                            "NON-TYPABLE", "NONTYPABLE", "NOT TYPED")] <- NA_character_
  x
}

unordered_pair_key <- function(pid_a, tp_a, pid_b, tp_b) {
  ifelse(
    paste(pid_a, tp_a, sep = "|") <= paste(pid_b, tp_b, sep = "|"),
    paste(pid_a, tp_a, pid_b, tp_b, sep = "|"),
    paste(pid_b, tp_b, pid_a, tp_a, sep = "|")
  )
}

required_files <- c(
  cohort = FILE_ANALYSIS_CLINICAL_COHORT,
  vf_ready = FILE_VF_READY,
  vf_pa = FILE_VF_PA,
  mlst = FILE_MLST_CANONICAL,
  pairwise = file.path(DIR_STRAIN, "pairwise_metrics.csv")
)
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files)) stop("Missing required file(s): ", paste(missing_files, collapse = ", "))

longcycler_raw <- read_csv(required_files[["cohort"]], show_col_types = FALSE)
if (!all(c("Participant_id", "tp_lab", "fasta_path", "fasta_sha256") %in% names(longcycler_raw))) {
  stop("Selected Longcycler analysis cohort lacks endpoint path/hash provenance")
}
longcycler <- longcycler_raw %>%
  prefer_primary_uti_status(allow_legacy_fallback = FALSE) %>%
  mutate(
    Participant_id = as.character(Participant_id),
    tp_lab = normalise_timepoint_preserve_events(tp_lab),
    assembler = str_to_lower(coalesce(
      if ("assembler" %in% names(.)) as.character(assembler) else NA_character_,
      if ("Assembler" %in% names(.)) as.character(Assembler) else NA_character_
    )),
    fasta_path = coalesce(
      if ("fasta_path" %in% names(.)) as.character(fasta_path) else NA_character_,
      if ("full_path" %in% names(.)) as.character(full_path) else NA_character_
    ),
    fasta_sha256 = as.character(.data$fasta_sha256)
  )

if (nrow(longcycler) == 0) stop("Selected Longcycler analysis cohort is empty")
if (any(is.na(longcycler$assembler) | longcycler$assembler != ANALYSIS_ASSEMBLER)) {
  stop("Analysis cohort contains an assembler outside the selected Longcycler policy")
}
if (any(is.na(longcycler$fasta_path) | !file.exists(longcycler$fasta_path))) {
  stop("Longcycler-primary manifest contains a missing selected FASTA")
}
longcycler <- longcycler %>%
  mutate(fasta_path = normalizePath(.data$fasta_path, winslash = "/", mustWork = TRUE))
current_cohort_hashes <- file_content_sha256(longcycler$fasta_path)
if (any(is.na(longcycler$fasta_sha256) | !nzchar(longcycler$fasta_sha256)) ||
    any(is.na(current_cohort_hashes)) ||
    any(tolower(longcycler$fasta_sha256) != tolower(current_cohort_hashes))) {
  stop("Selected Longcycler cohort FASTA hashes do not match current file content")
}
if (anyDuplicated(longcycler[c("Participant_id", "tp_lab")])) {
  stop("Longcycler-primary manifest contains duplicate participant-timepoint rows")
}

vf_ready <- read_csv(required_files[["vf_ready"]], show_col_types = FALSE) %>%
  mutate(Participant_id = as.character(Participant_id),
         tp_lab = normalise_timepoint_preserve_events(tp_lab))
vf_keys <- vf_ready %>% distinct(Participant_id, tp_lab)
cohort_keys <- longcycler %>% distinct(Participant_id, tp_lab)
if (nrow(vf_keys) != nrow(vf_ready) ||
    nrow(anti_join(cohort_keys, vf_keys, by = c("Participant_id", "tp_lab"))) ||
    nrow(anti_join(vf_keys, cohort_keys, by = c("Participant_id", "tp_lab")))) {
  stop("VF-ready keys do not exactly equal the selected Longcycler analysis cohort")
}

vf_pa <- read_csv(required_files[["vf_pa"]], show_col_types = FALSE)
gene_cols <- canonical_vf_gene_cols(names(vf_pa), vf_pa_file = required_files[["vf_pa"]])

mlst <- read_csv(required_files[["mlst"]], show_col_types = FALSE) %>%
  mutate(Participant_id = as.character(Participant_id),
         tp_lab = normalise_timepoint_preserve_events(tp_lab), ST = normalise_st(ST))

pairwise <- read_csv(required_files[["pairwise"]], show_col_types = FALSE) %>%
  mutate(
    Participant_id_A = as.character(Participant_id_A),
    Participant_id_B = as.character(Participant_id_B),
    Timepoint_A = normalise_timepoint_preserve_events(Timepoint_A),
    Timepoint_B = normalise_timepoint_preserve_events(Timepoint_B),
    pair_key = unordered_pair_key(Participant_id_A, Timepoint_A, Participant_id_B, Timepoint_B)
  )

if (anyDuplicated(pairwise$pair_key)) stop("pairwise_metrics.csv contains duplicate unordered endpoint pairs")

expected_pairwise_n <- longcycler %>%
  count(Participant_id, name = "n_episodes") %>%
  summarise(n_pairs = sum(n_episodes * (n_episodes - 1L) / 2L)) %>%
  pull(n_pairs)
if (nrow(pairwise) != expected_pairwise_n) {
  stop(
    "Primary pairwise denominator does not match all within-participant combinations from the Longcycler manifest: expected ",
    expected_pairwise_n, ", observed ", nrow(pairwise)
  )
}

provenance_cols <- c(
  "Assembler_A", "Assembler_B", "Fasta_path_A", "Fasta_path_B",
  "Fasta_SHA256_A", "Fasta_SHA256_B", "dnadiff_report_path",
  "dnadiff_sidecar_path", "dnadiff_report_sha256", "dnadiff_cache_signature",
  "dnadiff_cache_status", "dnadiff_version"
)
plasmid_pair_cols <- c(
  "Replicon_Jaccard", "Replicon_Both_Empty",
  "Replicon_Profile_Available", "Replicon_Gains", "Replicon_Losses",
  "Replicon_Gain_Count", "Replicon_Loss_Count",
  "MOB_Cluster_Jaccard", "MOB_Cluster_Both_Empty",
  "MOB_Profile_Available", "MOB_Cluster_Gains", "MOB_Cluster_Losses",
  "MOB_Cluster_Gain_Count", "MOB_Cluster_Loss_Count",
  "Predicted_Plasmid_Count_A", "Predicted_Plasmid_Count_B",
  "Predicted_Plasmid_Count_Difference",
  "MOB_High_Confidence_Profile_A", "MOB_High_Confidence_Profile_B",
  "MOB_High_Confidence_Profile_Both"
)
missing_provenance <- setdiff(provenance_cols, names(pairwise))
if (length(missing_provenance)) {
  stop("Pairwise output lacks required cache-provenance columns: ", paste(missing_provenance, collapse = ", "))
}
missing_plasmid_pair <- setdiff(plasmid_pair_cols, names(pairwise))
if (length(missing_plasmid_pair)) {
  stop(
    "Pairwise output lacks required corrected plasmid columns: ",
    paste(missing_plasmid_pair, collapse = ", ")
  )
}
if (any(is.na(pairwise$Fasta_SHA256_A) | pairwise$Fasta_SHA256_A == "" |
        is.na(pairwise$Fasta_SHA256_B) | pairwise$Fasta_SHA256_B == "")) {
  stop("Pairwise output contains missing FASTA fingerprints")
}
if (any(!file.exists(pairwise$dnadiff_report_path))) stop("Pairwise output points to missing dnadiff reports")
if (any(str_to_lower(pairwise$Assembler_A) != "longcycler" |
        str_to_lower(pairwise$Assembler_B) != "longcycler")) {
  stop("Primary pairwise output contains a non-Longcycler endpoint")
}

cohort_endpoint_lookup <- longcycler %>%
  transmute(
    endpoint_key = paste(.data$Participant_id, .data$tp_lab, sep = "|"),
    cohort_path = .data$fasta_path,
    cohort_sha256 = tolower(.data$fasta_sha256)
  )
audit_pairwise_endpoint <- function(key_pid, key_tp, pair_path, pair_hash, label) {
  pair_path <- as.character(pair_path)
  pair_path_norm <- vapply(pair_path, function(p) {
    if (is.na(p) || !nzchar(p) || !file.exists(p)) return(NA_character_)
    normalizePath(p, winslash = "/", mustWork = TRUE)
  }, character(1))
  endpoint <- tibble(
    endpoint_key = paste(as.character(key_pid), normalise_timepoint_preserve_events(key_tp), sep = "|"),
    pair_path = pair_path_norm,
    pair_sha256 = tolower(as.character(pair_hash))
  ) %>%
    left_join(cohort_endpoint_lookup, by = "endpoint_key", relationship = "many-to-one")
  bad <- endpoint %>%
    filter(is.na(.data$cohort_path) | is.na(.data$pair_path) |
             .data$pair_path != .data$cohort_path |
             is.na(.data$pair_sha256) | .data$pair_sha256 != .data$cohort_sha256)
  if (nrow(bad) > 0) {
    stop("Pairwise endpoint ", label, " has ", nrow(bad),
         " key/path/hash row(s) outside the exact selected Longcycler cohort")
  }
  invisible(TRUE)
}
audit_pairwise_endpoint(pairwise$Participant_id_A, pairwise$Timepoint_A,
                        pairwise$Fasta_path_A, pairwise$Fasta_SHA256_A, "A")
audit_pairwise_endpoint(pairwise$Participant_id_B, pairwise$Timepoint_B,
                        pairwise$Fasta_path_B, pairwise$Fasta_SHA256_B, "B")

lc_rows <- vf_ready %>%
  select(-any_of(c("Assembly_ID", "assembler", "Assembler", "fasta_path", "full_path"))) %>%
  inner_join(
    longcycler %>%
      transmute(Participant_id, tp_lab, Assembly_ID,
                Isolate_ID_canonical = if ("Isolate_ID" %in% names(.)) as.character(Isolate_ID) else NA_character_,
                assembler, fasta_path),
    by = c("Participant_id", "tp_lab"), relationship = "one-to-one"
  )

if (nrow(lc_rows) != nrow(longcycler)) {
  stop(
    "Longcycler-primary manifest-to-VF join changed the denominator: manifest=",
    nrow(longcycler), ", VF-ready=", nrow(lc_rows)
  )
}
if (anyDuplicated(lc_rows[c("Participant_id", "tp_lab")])) stop("Duplicate Longcycler episode keys")
if (nrow(lc_rows) != 532L || n_distinct(lc_rows$Participant_id) != 161L ||
    sum(lc_rows$Infection_Status == "UTI", na.rm = TRUE) != 16L ||
    sum(lc_rows$Infection_Status == "Not_UTI", na.rm = TRUE) != 516L) {
  stop("Selected Longcycler cohort contract failed: expected 532 episodes, 161 residents, 16 UTI, and 516 Not_UTI")
}

lc_rows <- lc_rows %>%
  mutate(
    Collection_Date_parsed = parse_date_safe(Collection_Date),
    fallback_order = fallback_time_order(tp_lab)
  ) %>%
  group_by(Participant_id) %>%
  mutate(
    first_collection_date = if (all(is.na(Collection_Date_parsed))) as.Date(NA) else min(Collection_Date_parsed, na.rm = TRUE),
    date_order = as.numeric(Collection_Date_parsed - first_collection_date),
    time_order = coalesce(date_order, fallback_order),
    time_order_source = case_when(
      !is.na(date_order) ~ "Collection_Date",
      !is.na(fallback_order) & str_detect(tp_lab, regex("uricult|^UTI-", ignore_case = TRUE)) ~ "Fallback_event_label_order",
      !is.na(fallback_order) ~ "Fallback_routine_timepoint_order",
      TRUE ~ "Unavailable"
    )
  ) %>%
  ungroup() %>%
  arrange(Participant_id, time_order, tp_lab)

transitions <- lc_rows %>%
  group_by(Participant_id) %>%
  arrange(time_order, tp_lab, .by_group = TRUE) %>%
  mutate(
    tp_to = lead(tp_lab),
    status_to = lead(Infection_Status),
    collection_date_to = lead(Collection_Date_parsed),
    time_order_to = lead(time_order),
    time_order_source_to = lead(time_order_source),
    ST_to = lead(ST),
    assembler_to = lead(assembler),
    fasta_path_to = lead(fasta_path),
    Isolate_ID_to = lead(Isolate_ID_canonical)
  ) %>%
  ungroup() %>%
  filter(!is.na(tp_to)) %>%
  transmute(
    Participant_id,
    tp_from = tp_lab,
    tp_to,
    status_from = Infection_Status,
    status_to,
    transition_type = paste(status_from, status_to, sep = "\u2192"),
    collection_date_from = as.character(Collection_Date_parsed),
    collection_date_to = as.character(collection_date_to),
    days_between_samples = as.numeric(as.Date(collection_date_to) - as.Date(collection_date_from)),
    time_order_from = time_order,
    time_order_to,
    time_order_source = paste(time_order_source, time_order_source_to, sep = " -> "),
    comparison_scope = "adjacent in the Longcycler-primary timeline",
    ST_from = normalise_st(ST),
    ST_to = normalise_st(ST_to),
    assembler_from = assembler,
    assembler_to,
    fasta_path_from = fasta_path,
    fasta_path_to,
    Isolate_ID_from = Isolate_ID_canonical,
    Isolate_ID_to,
    pair_key = unordered_pair_key(Participant_id, tp_from, Participant_id, tp_to)
  ) %>%
  left_join(
    pairwise %>%
      select(
        pair_key, AvgIdentity, TotalSNPs, MashDistance,
        Classification, RuleUsed,
        legacy_accessory_composite_classification,
        legacy_accessory_composite_rule,
        strict_same_strain,
        snp_strain_context, st_lineage_context, pair_interpretation,
        all_of(plasmid_pair_cols),
        all_of(provenance_cols)),
    by = "pair_key"
  ) %>%
  mutate(
    same_strain_snp_threshold = strain_snp_threshold(),
    strict_same_strain = !is.na(TotalSNPs) &
      TotalSNPs <= same_strain_snp_threshold,
    Classification_is_canonical = FALSE
  )

expected_transition_n <- lc_rows %>%
  count(Participant_id, name = "n_episodes") %>%
  summarise(n_transitions = sum(pmax(n_episodes - 1L, 0L))) %>%
  pull(n_transitions)
expected_transition_participants <- lc_rows %>%
  count(Participant_id, name = "n_episodes") %>%
  filter(n_episodes >= 2L) %>%
  summarise(n = n()) %>%
  pull(n)
if (nrow(transitions) != expected_transition_n) {
  stop("Rebuilt Longcycler-primary transition count is inconsistent with the episode manifest")
}
if (n_distinct(transitions$Participant_id) != expected_transition_participants) {
  stop("Rebuilt Longcycler-primary transition participant count is inconsistent with the episode manifest")
}
if (nrow(transitions) != 371L || n_distinct(transitions$Participant_id) != 139L) {
  stop("Selected Longcycler transition contract failed: expected 371 transitions from 139 residents")
}
if (any(is.na(transitions$TotalSNPs))) stop("Fresh pairwise SNP evidence is missing for one or more Longcycler transitions")
if (any(transitions$Assembler_A != "longcycler" | transitions$Assembler_B != "longcycler")) {
  stop("Longcycler transition joined to a non-Longcycler pairwise endpoint")
}

not_uti_to_uti <- transitions %>% filter(status_from == "Not_UTI", status_to == "UTI")
if (nrow(not_uti_to_uti) != 9L) {
  stop("Selected Longcycler transition contract failed: expected 9 Not_UTI-to-UTI transitions")
}

target_checks <- tibble::tribble(
  ~Participant_id, ~tp_from, ~tp_to, ~expected_snps,
  "122006", "T2", "UTI-1", 79,
  "20011", "T5", "UTI-1", 97
) %>%
  left_join(
    transitions %>% select(Participant_id, tp_from, tp_to, observed_snps = TotalSNPs),
    by = c("Participant_id", "tp_from", "tp_to")
  ) %>%
  mutate(pass = !is.na(observed_snps) & observed_snps == expected_snps)
if (!all(target_checks$pass)) stop("Targeted Longcycler dnadiff regression check failed")

threshold_sensitivity <- tidyr::crossing(
  threshold = c(5L, 10L, 25L, 50L, 100L),
  analysis = c("all_longcycler_transitions", "not_uti_to_uti")
) %>%
  rowwise() %>%
  mutate(
    denominator = if (analysis == "all_longcycler_transitions") nrow(transitions) else nrow(not_uti_to_uti),
    n_at_or_below_threshold = if (analysis == "all_longcycler_transitions") {
      sum(transitions$TotalSNPs <= threshold, na.rm = TRUE)
    } else {
      sum(not_uti_to_uti$TotalSNPs <= threshold, na.rm = TRUE)
    },
    proportion = n_at_or_below_threshold / denominator
  ) %>%
  ungroup()

typed_lc <- lc_rows %>% mutate(ST = normalise_st(ST)) %>% filter(!is.na(ST))

count_rows <- tibble::tribble(
  ~metric, ~value, ~unit, ~source_file, ~filter, ~denominator, ~formula, ~verification_status,
  "selected_longcycler_rows", nrow(longcycler), "genome-linked episode rows", required_files[["cohort"]], "canonical analysis cohort", "selected QC-pass Longcycler rows", "nrow(cohort)", "verified",
  "selected_longcycler_participants", n_distinct(longcycler$Participant_id), "participants", required_files[["cohort"]], "canonical analysis cohort", "selected QC-pass Longcycler rows", "n_distinct(Participant_id)", "verified",
  "longcycler_uti", sum(lc_rows$Infection_Status == "UTI", na.rm = TRUE), "episode rows", required_files[["cohort"]], "Infection_Status == UTI", "selected Longcycler rows", "sum(condition)", "verified",
  "longcycler_not_uti", sum(lc_rows$Infection_Status == "Not_UTI", na.rm = TRUE), "episode rows", required_files[["cohort"]], "Infection_Status == Not_UTI", "selected Longcycler rows", "sum(condition)", "verified",
  "longcycler_mlst_typed_rows", nrow(typed_lc), "episode rows", required_files[["mlst"]], "preferred ST non-missing", "selected Longcycler rows", "sum(!is.na(ST))", "verified",
  "longcycler_distinct_st", n_distinct(typed_lc$ST), "ST labels", required_files[["mlst"]], "preferred ST non-missing", "typed selected Longcycler rows", "n_distinct(ST)", "verified",
  "vf_binary_features", length(gene_cols), "VFDB-derived binary features", required_files[["vf_pa"]], "canonical VF gene columns", "feature columns, not rows", "length(canonical_vf_gene_cols)", "verified",
  "longcycler_longitudinal_transitions", nrow(transitions), "transitions", file.path(outdir, "longcycler_transitions.csv"), "adjacent selected episodes", "selected Longcycler rows", "sum(n_i - 1) for participants with n_i >= 2", "verified",
  "longcycler_transition_participants", n_distinct(transitions$Participant_id), "participants", file.path(outdir, "longcycler_transitions.csv"), "participants with >=2 retained episodes", "371 transitions", "n_distinct(Participant_id)", "verified",
  "longcycler_strict_same_strain_transitions", sum(transitions$strict_same_strain), "transitions", file.path(outdir, "longcycler_transitions.csv"), "TotalSNPs <= 25", "371 Longcycler transitions", "sum(TotalSNPs <= 25)", "verified_after_fresh_dnadiff_rerun",
  "longcycler_not_uti_to_uti", nrow(not_uti_to_uti), "transitions", file.path(outdir, "not_uti_to_uti_transitions.csv"), "status_from == Not_UTI; status_to == UTI", "371 Longcycler transitions", "nrow(filtered transitions)", "verified",
  "longcycler_not_uti_to_uti_strict", sum(not_uti_to_uti$strict_same_strain), "transitions", file.path(outdir, "not_uti_to_uti_transitions.csv"), "Not_UTI to UTI; TotalSNPs <= 25", "9 Not_UTI-to-UTI transitions", "sum(TotalSNPs <= 25)", "verified_after_fresh_dnadiff_rerun"
)

write_csv(lc_rows, file.path(outdir, "longcycler_episode_manifest.csv"))
write_csv(transitions, file.path(outdir, "longcycler_transitions.csv"))
write_csv(not_uti_to_uti, file.path(outdir, "not_uti_to_uti_transitions.csv"))
write_csv(threshold_sensitivity, file.path(outdir, "snp_threshold_context.csv"))
write_csv(target_checks, file.path(outdir, "targeted_dnadiff_regression_checks.csv"))
write_csv(count_rows, file.path(outdir, "longcycler_denominator_counts.csv"))

summary_lines <- c(
  "Selected Longcycler transition analysis",
  sprintf("Generated: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  sprintf("Selected Longcycler rows: %d from %d participants", nrow(longcycler), n_distinct(longcycler$Participant_id)),
  sprintf("Longcycler clinical labels: %d UTI; %d Not_UTI", sum(lc_rows$Infection_Status == "UTI"), sum(lc_rows$Infection_Status == "Not_UTI")),
  sprintf("Preferred MLST calls linked to Longcycler-selected episodes: %d typed rows; %d ST labels", nrow(typed_lc), n_distinct(typed_lc$ST)),
  sprintf("VFDB-derived binary features: %d", length(gene_cols)),
  sprintf("Longcycler transitions: %d from %d participants", nrow(transitions), n_distinct(transitions$Participant_id)),
  sprintf("Operational <=%d dnadiff-SNP rule: %d/%d Longcycler transitions", strain_snp_threshold(), sum(transitions$strict_same_strain), nrow(transitions)),
  sprintf("Not_UTI-to-UTI: %d transitions; %d meet the operational <=%d dnadiff-SNP rule", nrow(not_uti_to_uti), sum(not_uti_to_uti$strict_same_strain), strain_snp_threshold()),
  "Interpretation boundary: assembly-to-assembly dnadiff SNP differences are not Parsnp core-genome SNP distances or wgMLST allele distances."
)
writeLines(summary_lines, file.path(outdir, "longcycler_transition_summary.txt"))

message("Selected Longcycler transition rebuild complete: ", outdir)

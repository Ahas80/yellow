#!/usr/bin/env Rscript
# ==============================================================================
# Rebuild the Longcycler-only sensitivity denominator and longitudinal pairs
# ==============================================================================
# This is a restricted sensitivity analysis of the executed mixed-canonical
# workflow. It does not redefine the project's primary clinical labels and it
# does not replace the mixed 556-row analysis.

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
  library(tidyr)
})

source("00_config.R")

outdir <- file.path("results", "sensitivity", "longcycler_only")
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
  canonical = file.path(DIR_QC, "canonical_assembly_selection.csv"),
  status = FILE_STATUS_MAP,
  vf_ready = FILE_VF_READY,
  vf_pa = FILE_VF_PA,
  mlst = FILE_MLST_CANONICAL,
  pairwise = file.path(DIR_STRAIN, "pairwise_metrics.csv")
)
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files)) stop("Missing required file(s): ", paste(missing_files, collapse = ", "))

canonical <- read_csv(required_files[["canonical"]], show_col_types = FALSE) %>%
  mutate(
    Participant_id = as.character(Participant_id),
    tp_lab = as.character(tp_lab),
    selected_canonical = as_flag(selected_canonical),
    QC_PASS = as_flag(QC_PASS),
    assembler = str_to_lower(coalesce(as.character(assembler), as.character(Assembler))),
    fasta_path = coalesce(as.character(fasta_path), as.character(full_path))
  )

selected <- canonical %>%
  filter(selected_canonical, QC_PASS, file.exists(fasta_path))

longcycler <- selected %>%
  filter(assembler == "longcycler") %>%
  distinct(Participant_id, tp_lab, .keep_all = TRUE)

if (nrow(selected) != 556L) stop("Expected 556 selected canonical rows; observed ", nrow(selected))
if (nrow(longcycler) != 532L) stop("Expected 532 selected Longcycler rows; observed ", nrow(longcycler))
if (n_distinct(longcycler$Participant_id) != 161L) stop("Expected 161 Longcycler participants")

status <- read_csv(required_files[["status"]], show_col_types = FALSE) %>%
  prefer_primary_uti_status() %>%
  mutate(
    Participant_id = as.character(Participant_id),
    tp_lab = as.character(tp_lab),
    analysis_include_primary = as_flag(analysis_include_primary)
  )

vf_ready <- read_csv(required_files[["vf_ready"]], show_col_types = FALSE) %>%
  mutate(Participant_id = as.character(Participant_id), tp_lab = as.character(tp_lab))

vf_pa <- read_csv(required_files[["vf_pa"]], show_col_types = FALSE)
gene_cols <- canonical_vf_gene_cols(names(vf_pa), vf_pa_file = required_files[["vf_pa"]])

mlst <- read_csv(required_files[["mlst"]], show_col_types = FALSE) %>%
  mutate(Participant_id = as.character(Participant_id), tp_lab = as.character(tp_lab), ST = normalise_st(ST))

pairwise <- read_csv(required_files[["pairwise"]], show_col_types = FALSE) %>%
  mutate(
    Participant_id_A = as.character(Participant_id_A),
    Participant_id_B = as.character(Participant_id_B),
    Timepoint_A = as.character(Timepoint_A),
    Timepoint_B = as.character(Timepoint_B),
    pair_key = unordered_pair_key(Participant_id_A, Timepoint_A, Participant_id_B, Timepoint_B)
  )

if (nrow(pairwise) != 963L) stop("Expected 963 current canonical within-participant pairs; observed ", nrow(pairwise))
if (anyDuplicated(pairwise$pair_key)) stop("pairwise_metrics.csv contains duplicate unordered endpoint pairs")

provenance_cols <- c(
  "Assembler_A", "Assembler_B", "Fasta_path_A", "Fasta_path_B",
  "Fasta_SHA256_A", "Fasta_SHA256_B", "dnadiff_report_path",
  "dnadiff_sidecar_path", "dnadiff_report_sha256", "dnadiff_cache_signature",
  "dnadiff_cache_status", "dnadiff_version"
)
missing_provenance <- setdiff(provenance_cols, names(pairwise))
if (length(missing_provenance)) {
  stop("Pairwise output lacks required cache-provenance columns: ", paste(missing_provenance, collapse = ", "))
}
if (any(is.na(pairwise$Fasta_SHA256_A) | pairwise$Fasta_SHA256_A == "" |
        is.na(pairwise$Fasta_SHA256_B) | pairwise$Fasta_SHA256_B == "")) {
  stop("Pairwise output contains missing FASTA fingerprints")
}
if (any(!file.exists(pairwise$dnadiff_report_path))) stop("Pairwise output points to missing dnadiff reports")

lc_keys <- longcycler %>% select(Participant_id, tp_lab)
lc_rows <- vf_ready %>%
  inner_join(
    longcycler %>%
      select(Participant_id, tp_lab, Assembly_ID, Isolate_ID_canonical = Isolate_ID,
             assembler, fasta_path, QC_PASS, selected_canonical),
    by = c("Participant_id", "tp_lab")
  )

if (nrow(lc_rows) != 532L) stop("Longcycler-to-VF join did not retain exactly 532 rows")
if (anyDuplicated(lc_rows[c("Participant_id", "tp_lab")])) stop("Duplicate Longcycler episode keys")

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
    comparison_scope = "adjacent among retained Longcycler-selected episodes",
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
      select(pair_key, AvgIdentity, TotalSNPs, MashDistance, Classification, RuleUsed,
             snp_strain_context, st_lineage_context, pair_interpretation,
             all_of(provenance_cols)),
    by = "pair_key"
  ) %>%
  mutate(
    same_strain_snp_threshold = strain_snp_threshold(),
    strict_same_strain = !is.na(TotalSNPs) & TotalSNPs <= same_strain_snp_threshold
  )

if (nrow(transitions) != 371L) stop("Expected 371 Longcycler transitions; observed ", nrow(transitions))
if (n_distinct(transitions$Participant_id) != 139L) stop("Expected 139 transition participants")
if (any(is.na(transitions$TotalSNPs))) stop("Fresh pairwise SNP evidence is missing for one or more Longcycler transitions")
if (any(transitions$Assembler_A != "longcycler" | transitions$Assembler_B != "longcycler")) {
  stop("Longcycler transition joined to a non-Longcycler pairwise endpoint")
}

not_uti_to_uti <- transitions %>% filter(status_from == "Not_UTI", status_to == "UTI")
if (nrow(not_uti_to_uti) != 9L) stop("Expected 9 Longcycler Not_UTI-to-UTI transitions; observed ", nrow(not_uti_to_uti))

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

clinical_primary <- status %>% filter(analysis_include_primary)
typed_lc <- lc_rows %>% mutate(ST = normalise_st(ST)) %>% filter(!is.na(ST))
mixed_transitions <- read_csv(file.path(DIR_VF, "vf_longitudinal_transitions.csv"), show_col_types = FALSE)
if ("cohort" %in% names(mixed_transitions)) mixed_transitions <- mixed_transitions %>% filter(cohort == "all")
mixed_transitions <- mixed_transitions %>%
  mutate(
    Participant_id = as.character(Participant_id),
    pair_key = unordered_pair_key(Participant_id, as.character(tp_from), Participant_id, as.character(tp_to))
  )
n_new_longcycler_adjacencies <- sum(!transitions$pair_key %in% mixed_transitions$pair_key)
if (n_new_longcycler_adjacencies != 9L) {
  stop("Expected 9 newly adjacent Longcycler pairs after removing Flye rows; observed ", n_new_longcycler_adjacencies)
}

count_rows <- tibble::tribble(
  ~metric, ~value, ~unit, ~source_file, ~filter, ~denominator, ~formula, ~verification_status,
  "clinical_records_before_primary_exclusions", nrow(status), "clinical records", required_files[["status"]], "none", "all status-map rows", "nrow(status_map)", "verified",
  "primary_clinical_episodes", nrow(clinical_primary), "clinical episodes", required_files[["status"]], "analysis_include_primary == TRUE", "eligible clinical episodes", "sum(analysis_include_primary)", "verified",
  "primary_clinical_participants", n_distinct(clinical_primary$Participant_id), "participants", required_files[["status"]], "analysis_include_primary == TRUE", "eligible clinical participants", "n_distinct(Participant_id)", "verified",
  "primary_clinical_uti", sum(clinical_primary$Infection_Status == "UTI", na.rm = TRUE), "episodes", required_files[["status"]], "analysis_include_primary == TRUE; Infection_Status == UTI", "583 primary clinical episodes", "sum(condition)", "verified",
  "primary_clinical_not_uti", sum(clinical_primary$Infection_Status == "Not_UTI", na.rm = TRUE), "episodes", required_files[["status"]], "analysis_include_primary == TRUE; Infection_Status == Not_UTI", "583 primary clinical episodes", "sum(condition)", "verified",
  "candidate_assembly_records", nrow(canonical), "assembly records", required_files[["canonical"]], "none", "all candidate assembler alternatives", "nrow(canonical selection table)", "verified",
  "mixed_canonical_rows", nrow(selected), "genome-linked episode rows", required_files[["canonical"]], "selected_canonical == TRUE; QC_PASS == TRUE; FASTA exists", "selected mixed-canonical analysis", "nrow(selected)", "verified",
  "mixed_canonical_participants", n_distinct(selected$Participant_id), "participants", required_files[["canonical"]], "selected canonical QC-pass", "556 mixed rows", "n_distinct(Participant_id)", "verified",
  "selected_longcycler_rows", nrow(longcycler), "genome-linked episode rows", required_files[["canonical"]], "selected canonical QC-pass; assembler == longcycler", "Longcycler sensitivity set", "nrow(longcycler)", "verified",
  "selected_longcycler_participants", n_distinct(longcycler$Participant_id), "participants", required_files[["canonical"]], "selected canonical QC-pass; assembler == longcycler", "532 Longcycler rows", "n_distinct(Participant_id)", "verified",
  "longcycler_uti", sum(lc_rows$Infection_Status == "UTI", na.rm = TRUE), "episode rows", required_files[["vf_ready"]], "Longcycler keys; Infection_Status == UTI", "532 Longcycler rows", "sum(condition)", "verified",
  "longcycler_not_uti", sum(lc_rows$Infection_Status == "Not_UTI", na.rm = TRUE), "episode rows", required_files[["vf_ready"]], "Longcycler keys; Infection_Status == Not_UTI", "532 Longcycler rows", "sum(condition)", "verified",
  "longcycler_mlst_typed_rows", nrow(typed_lc), "episode rows", required_files[["mlst"]], "Longcycler keys; preferred ST non-missing", "532 Longcycler rows", "sum(!is.na(ST))", "verified",
  "longcycler_distinct_st", n_distinct(typed_lc$ST), "ST labels", required_files[["mlst"]], "Longcycler keys; preferred ST non-missing", "514 typed Longcycler-selected rows", "n_distinct(ST)", "verified",
  "vf_binary_features", length(gene_cols), "VFDB-derived binary features", required_files[["vf_pa"]], "canonical VF gene columns", "feature columns, not rows", "length(canonical_vf_gene_cols)", "verified",
  "mixed_longitudinal_transitions", nrow(mixed_transitions), "transitions", file.path(DIR_VF, "vf_longitudinal_transitions.csv"), "cohort == all when present", "mixed 556-row timeline", "nrow(filtered transition table)", "verified",
  "longcycler_longitudinal_transitions", nrow(transitions), "transitions", file.path(outdir, "longcycler_transitions.csv"), "adjacent after Longcycler restriction", "532 Longcycler rows", "sum(n_i - 1) for participants with n_i >= 2", "verified",
  "longcycler_transition_participants", n_distinct(transitions$Participant_id), "participants", file.path(outdir, "longcycler_transitions.csv"), "participants with >=2 retained episodes", "371 transitions", "n_distinct(Participant_id)", "verified",
  "longcycler_transitions_already_adjacent_in_mixed", nrow(transitions) - n_new_longcycler_adjacencies, "transitions", file.path(outdir, "longcycler_transitions.csv"), "Longcycler pair key also present in mixed transition table", "371 rebuilt Longcycler transitions", "sum(pair_key %in% mixed_pair_keys)", "verified",
  "new_longcycler_adjacencies_after_restriction", n_new_longcycler_adjacencies, "transitions", file.path(outdir, "longcycler_transitions.csv"), "Longcycler pair key absent from mixed transition table", "371 rebuilt Longcycler transitions", "sum(!pair_key %in% mixed_pair_keys)", "verified",
  "longcycler_strict_same_strain_transitions", sum(transitions$strict_same_strain), "transitions", file.path(outdir, "longcycler_transitions.csv"), "TotalSNPs <= 25", "371 Longcycler transitions", "sum(TotalSNPs <= 25)", "verified_after_fresh_dnadiff_rerun",
  "longcycler_not_uti_to_uti", nrow(not_uti_to_uti), "transitions", file.path(outdir, "not_uti_to_uti_transitions.csv"), "status_from == Not_UTI; status_to == UTI", "371 Longcycler transitions", "nrow(filtered transitions)", "verified",
  "longcycler_not_uti_to_uti_strict", sum(not_uti_to_uti$strict_same_strain), "transitions", file.path(outdir, "not_uti_to_uti_transitions.csv"), "Not_UTI to UTI; TotalSNPs <= 25", "9 Not_UTI-to-UTI transitions", "sum(TotalSNPs <= 25)", "verified_after_fresh_dnadiff_rerun"
)

write_csv(lc_rows, file.path(outdir, "longcycler_episode_manifest.csv"))
write_csv(transitions, file.path(outdir, "longcycler_transitions.csv"))
write_csv(not_uti_to_uti, file.path(outdir, "not_uti_to_uti_transitions.csv"))
write_csv(threshold_sensitivity, file.path(outdir, "snp_threshold_sensitivity.csv"))
write_csv(target_checks, file.path(outdir, "targeted_dnadiff_regression_checks.csv"))
write_csv(count_rows, file.path(outdir, "denominator_counts.csv"))

summary_lines <- c(
  "Longcycler-only sensitivity analysis",
  sprintf("Generated: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  "Framing: executed mixed-canonical pipeline plus a restricted Longcycler-only sensitivity analysis.",
  sprintf("Selected mixed canonical rows: %d", nrow(selected)),
  sprintf("Selected Longcycler rows: %d from %d participants", nrow(longcycler), n_distinct(longcycler$Participant_id)),
  sprintf("Longcycler clinical labels: %d UTI; %d Not_UTI", sum(lc_rows$Infection_Status == "UTI"), sum(lc_rows$Infection_Status == "Not_UTI")),
  sprintf("Preferred MLST calls linked to Longcycler-selected episodes: %d typed rows; %d ST labels", nrow(typed_lc), n_distinct(typed_lc$ST)),
  sprintf("VFDB-derived binary features: %d", length(gene_cols)),
  sprintf("Longcycler transitions: %d from %d participants", nrow(transitions), n_distinct(transitions$Participant_id)),
  sprintf("New Longcycler-to-Longcycler adjacencies created by rebuilding after restriction: %d", n_new_longcycler_adjacencies),
  sprintf("Operational <=%d dnadiff-SNP rule: %d/%d Longcycler transitions", strain_snp_threshold(), sum(transitions$strict_same_strain), nrow(transitions)),
  sprintf("Not_UTI-to-UTI: %d transitions; %d meet the operational <=%d dnadiff-SNP rule", nrow(not_uti_to_uti), sum(not_uti_to_uti$strict_same_strain), strain_snp_threshold()),
  "Interpretation boundary: assembly-to-assembly dnadiff SNP differences are not Parsnp core-genome SNP distances or wgMLST allele distances."
)
writeLines(summary_lines, file.path(outdir, "longcycler_sensitivity_summary.txt"))

message("Longcycler-only sensitivity rebuild complete: ", outdir)

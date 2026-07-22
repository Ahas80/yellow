#!/usr/bin/env Rscript

# Build the single machine-readable claim registry used by every regenerated
# deck, handout, and codebook. Publication is fail-closed on the exact selected
# cohort and direct genomic evidence contracts.

source("00_config.R")
source("R/pipeline_qc_helpers.R")
suppressPackageStartupMessages({
  library(digest)
  library(dplyr)
  library(jsonlite)
  library(readr)
  library(stringr)
})

require_file <- function(path, label) {
  if (!file.exists(path)) stop(label, " is missing: ", path, call. = FALSE)
  path
}

normalise_path <- function(x, must_work = TRUE) {
  normalizePath(as.character(x), winslash = "/", mustWork = must_work)
}

sha256_file <- function(path) unname(digest(path, algo = "sha256", file = TRUE))

release_qc_config <- function() {
  list(
    MAX_CONTIGS = 200,
    MIN_N50 = 20000,
    MIN_GENOME_SIZE = 4.0e6,
    MAX_GENOME_SIZE = 6.0e6
  )
}

cohort_path <- require_file(FILE_ANALYSIS_CLINICAL_COHORT, "Selected clinical cohort")
manifest_path <- require_file(FILE_ANALYSIS_ASSEMBLY_MANIFEST, "Selected assembly manifest")
vf_path <- require_file(FILE_VF_PA, "VF presence/absence table")
mlst_path <- require_file(FILE_MLST_CANONICAL, "Canonical MLST table")
pair_path <- require_file(file.path(DIR_STRAIN, "pairwise_metrics.csv"), "Direct pair table")
transition_path <- require_file(file.path(DIR_RESULTS, "longitudinal", "longcycler_transitions.csv"), "Canonical transition table")
panaroo_manifest_path <- require_file(file.path(DIR_WGS_PAN, "panaroo_input_manifest.csv"), "Panaroo input manifest")
panaroo_roary_path <- require_file(file.path(DIR_WGS_PAN, "gene_presence_absence_roary.csv"), "Panaroo Roary-compatible table")
casebook_path <- require_file(file.path(DIR_RESULTS, "mechanism", "not_uti_to_uti_casebook.csv"), "Mechanism casebook")
near_miss_path <- require_file(file.path(DIR_RESULTS, "audit", "uti_not_uti_near_miss_rows.csv"), "Near-miss audit")
rq_marker <- require_file(file.path(DIR_RESULTS, "research_questions", "RUN_COMPLETE.txt"), "RQ01-RQ10 marker")
vf_method_path <- require_file(
  file.path(DIR_RESULTS, "research_questions", "_inputs", "method_manifest.csv"),
  "RQ06-RQ08 method manifest"
)
amr_profile_path <- require_file(
  file.path(DIR_RESULTS, "amr", "episode_amr_profiles.csv"),
  "Script-29 genomic AMR episode profiles"
)
amr_transition_path <- require_file(
  file.path(DIR_RESULTS, "amr", "adjacent_pair_amr_profiles_371.csv"),
  "Script-29 adjacent-pair AMR profiles"
)
amr_resident_path <- require_file(
  file.path(DIR_RESULTS, "amr", "resident_amr_profiles.csv"),
  "Script-29 genomic AMR resident profiles"
)
amr_focused_path <- require_file(
  file.path(DIR_RESULTS, "amr", "not_uti_to_uti_amr_profiles_9.csv"),
  "Script-29 focused AMR transition profiles"
)
amr_inference_path <- require_file(
  file.path(DIR_RESULTS, "amr", "longitudinal_resident_bootstrap_inference.csv"),
  "Script-29 AMR longitudinal inference"
)
amr_coverage_path <- require_file(
  file.path(DIR_RESULTS, "amr", "caller_coverage_summary.csv"),
  "Script-29 AMR caller coverage"
)
amr_table13_path <- require_file(
  file.path(DIR_RESULTS, "summary", "table_13_genomic_amr_summary.csv"),
  "Genomic AMR result summary table"
)
amr_validation_path <- require_file(
  file.path(DIR_RESULTS, "amr", "validation_checks.csv"),
  "Script-29 genomic AMR validation"
)
amr_provenance_path <- require_file(
  file.path(DIR_RESULTS, "amr", "provenance", "tool_database_versions.csv"),
  "Script-29 AMR tool/database provenance"
)
pf_run_path <- require_file(
  file.path(DIR_PLASMIDS, "plasmidfinder_run_manifest.csv"),
  "Canonical PlasmidFinder run manifest"
)
pf_database_path <- require_file(
  file.path(DIR_PLASMIDS, "plasmidfinder_database_manifest.csv"),
  "Canonical PlasmidFinder database manifest"
)
pf_pa_path <- require_file(
  file.path(DIR_PLASMIDS, "plasmidfinder_presence_absence.csv"),
  "Canonical gene-level PlasmidFinder matrix"
)
mob_root <- file.path(DIR_PLASMIDS, "mob_suite")
mob_marker_path <- require_file(
  file.path(mob_root, "RUN_COMPLETE.txt"),
  "MOB-suite completion marker"
)
mob_run_path <- require_file(
  file.path(mob_root, "run_manifest.csv"),
  "MOB-suite run manifest"
)
mob_contig_path <- require_file(
  file.path(mob_root, "contig_assignments.csv"),
  "MOB-suite contig assignments"
)
mob_plasmid_path <- require_file(
  file.path(mob_root, "plasmids_long.csv"),
  "MOB-suite predicted plasmids"
)
mob_profile_path <- require_file(
  file.path(mob_root, "episode_plasmid_profiles.csv"),
  "MOB-suite episode plasmid profiles"
)
mob_gene_location_path <- require_file(
  file.path(mob_root, "plasmid_gene_locations_long.csv"),
  "Predicted plasmid/chromosome gene locations"
)
mob_adjacent_path <- require_file(
  file.path(mob_root, "adjacent_pair_plasmid_metrics_371.csv"),
  "Adjacent-pair plasmid mechanism metrics"
)
mob_focused_path <- require_file(
  file.path(mob_root, "not_uti_to_uti_plasmid_metrics_9.csv"),
  "Focused plasmid mechanism metrics"
)
mob_mechanism_profile_path <- require_file(
  file.path(mob_root, "episode_mechanism_profiles.csv"),
  "Episode plasmid mechanism profiles"
)
mob_location_validation_path <- require_file(
  file.path(mob_root, "plasmid_gene_location_validation.csv"),
  "Predicted plasmid gene-location validation"
)
plasmid_summary_path <- require_file(
  file.path(DIR_RESULTS, "summary", "table_13b_predicted_plasmid_summary.csv"),
  "Predicted plasmid result summary table"
)

cohort <- read_csv(cohort_path, show_col_types = FALSE) %>%
  mutate(
    Participant_id = as.character(Participant_id),
    tp_lab = normalise_timepoint_preserve_events(tp_lab)
  )
manifest <- read_csv(manifest_path, show_col_types = FALSE) %>%
  mutate(
    Participant_id = as.character(Participant_id),
    tp_lab = normalise_timepoint_preserve_events(tp_lab),
    full_path = normalise_path(full_path)
  )
vf <- read_csv(vf_path, show_col_types = FALSE)
mlst <- read_csv(mlst_path, show_col_types = FALSE)
pairs <- read_csv(pair_path, show_col_types = FALSE)
transitions <- read_csv(transition_path, show_col_types = FALSE)
panaroo_manifest <- read_csv(panaroo_manifest_path, show_col_types = FALSE)
casebook <- read_csv(casebook_path, show_col_types = FALSE)
near_miss <- read_csv(near_miss_path, show_col_types = FALSE)
vf_methods <- read_csv(vf_method_path, show_col_types = FALSE) %>%
  mutate(parameter = as.character(parameter), value = as.character(value))
amr_profiles <- read_csv(amr_profile_path, show_col_types = FALSE)
amr_transitions <- read_csv(amr_transition_path, show_col_types = FALSE)
amr_residents <- read_csv(amr_resident_path, show_col_types = FALSE)
amr_focused <- read_csv(amr_focused_path, show_col_types = FALSE)
amr_inference <- read_csv(amr_inference_path, show_col_types = FALSE)
amr_coverage <- read_csv(amr_coverage_path, show_col_types = FALSE)
amr_validation <- read_csv(amr_validation_path, show_col_types = FALSE)
amr_provenance <- read_csv(amr_provenance_path, show_col_types = FALSE)
pf_run <- read_csv(pf_run_path, show_col_types = FALSE)
pf_database <- read_csv(pf_database_path, show_col_types = FALSE)
pf_pa <- read_csv(pf_pa_path, show_col_types = FALSE)
mob_run <- read_csv(mob_run_path, show_col_types = FALSE)
mob_contigs <- read_csv(mob_contig_path, show_col_types = FALSE)
mob_plasmids <- read_csv(mob_plasmid_path, show_col_types = FALSE)
mob_profiles <- read_csv(mob_profile_path, show_col_types = FALSE)
mob_gene_locations <- read_csv(mob_gene_location_path, show_col_types = FALSE)
mob_adjacent <- read_csv(mob_adjacent_path, show_col_types = FALSE)
mob_focused <- read_csv(mob_focused_path, show_col_types = FALSE)
mob_mechanism_profiles <- read_csv(
  mob_mechanism_profile_path, show_col_types = FALSE
)
mob_location_validation <- read_csv(
  mob_location_validation_path, show_col_types = FALSE
)

method_value <- function(parameter) {
  values <- vf_methods$value[vf_methods$parameter == parameter]
  if (length(values) != 1L || is.na(values) || !nzchar(values)) {
    stop("Method manifest must contain exactly one non-missing value for: ", parameter, call. = FALSE)
  }
  values[[1]]
}

method_number <- function(parameter) {
  value <- suppressWarnings(as.numeric(method_value(parameter)))
  if (length(value) != 1L || is.na(value)) {
    stop("Method manifest value is not numeric for: ", parameter, call. = FALSE)
  }
  value
}

cohort_keys <- paste(cohort$Participant_id, cohort$tp_lab, sep = "|")
manifest_keys <- paste(manifest$Participant_id, manifest$tp_lab, sep = "|")
if (nrow(cohort) != 532L || n_distinct(cohort$Participant_id) != 161L ||
    sum(cohort$UTI_Status == "UTI", na.rm = TRUE) != 16L ||
    sum(cohort$UTI_Status == "Not_UTI", na.rm = TRUE) != 516L ||
    anyDuplicated(cohort_keys) || !setequal(cohort_keys, manifest_keys)) {
  stop("Selected cohort contract failed; claim registry was not published.", call. = FALSE)
}

if (!all(c("fasta_sha256", "full_path") %in% names(manifest))) {
  stop("Selected manifest lacks content provenance.", call. = FALSE)
}
observed_manifest_hash <- vapply(manifest$full_path, sha256_file, character(1))
if (any(tolower(observed_manifest_hash) != tolower(manifest$fasta_sha256))) {
  stop("Selected manifest FASTA content changed before claim publication.", call. = FALSE)
}

if (!all(c("fasta_path", "fasta_sha256", "gff_path") %in% names(panaroo_manifest))) {
  stop("Panaroo manifest lacks exact FASTA/GFF provenance.", call. = FALSE)
}
panaroo_manifest <- panaroo_manifest %>%
  mutate(fasta_path = normalise_path(fasta_path))
selected_path_hash <- paste(manifest$full_path, tolower(manifest$fasta_sha256), sep = "\n")
panaroo_path_hash <- paste(panaroo_manifest$fasta_path, tolower(panaroo_manifest$fasta_sha256), sep = "\n")
if (nrow(panaroo_manifest) != nrow(manifest) ||
    !setequal(panaroo_path_hash, selected_path_hash)) {
  stop("Panaroo manifest does not exactly match the selected FASTA/SHA-256 set.", call. = FALSE)
}
panaroo_expected_samples <- sub(
  "\\.gff$", "", basename(as.character(panaroo_manifest$gff_path)),
  ignore.case = TRUE
)
panaroo_metadata_columns <- c(
  "Gene", "Non-unique Gene name", "Annotation", "No. isolates",
  "No. sequences", "Avg sequences per isolate", "Genome Fragment",
  "Order within Fragment", "Accessory Fragment",
  "Accessory Order with Fragment", "QC", "Min group size nuc",
  "Max group size nuc", "Avg group size nuc"
)
panaroo_roary_samples <- setdiff(
  names(read_csv(panaroo_roary_path, n_max = 0, show_col_types = FALSE)),
  panaroo_metadata_columns
)
if (length(panaroo_roary_samples) != nrow(manifest) ||
    !setequal(panaroo_roary_samples, panaroo_expected_samples)) {
  stop("Panaroo Roary-compatible sample columns do not exactly match the selected GFF set.", call. = FALSE)
}

if (nrow(vf) != 532L || nrow(mlst) != 532L || nrow(pairs) != 893L ||
    nrow(transitions) != 371L || n_distinct(transitions$Participant_id) != 139L) {
  stop("Genomic denominator contract failed; claim registry was not published.", call. = FALSE)
}
if (nrow(amr_profiles) != 532L ||
    n_distinct(amr_profiles$Participant_id) != 161L ||
    nrow(amr_residents) != 161L ||
    nrow(amr_transitions) != 371L ||
    nrow(amr_focused) != 9L ||
    nrow(amr_inference) != 3L ||
    nrow(amr_coverage) != 3L ||
    any(amr_coverage$n_complete != 532L) ||
    any(!amr_validation$pass & amr_validation$severity == "critical")) {
  stop("Genomic AMR release contract failed; claim registry was not published.",
       call. = FALSE)
}
if (nrow(pf_run) != 532L ||
    any(pf_run$call_status != "complete") ||
    sum(as_pipeline_bool(pf_run$no_hit)) != 110L ||
    nrow(pf_pa) != 532L ||
    nrow(mob_run) != 532L ||
    any(mob_run$status != "complete") ||
    nrow(mob_profiles) != 532L ||
    nrow(mob_mechanism_profiles) != 532L ||
    nrow(mob_adjacent) != 371L ||
    nrow(mob_focused) != 9L ||
    any(!mob_location_validation$pass) ||
    anyNA(mob_gene_locations$localization) ||
    anyDuplicated(mob_contigs[c("Isolate_ID", "contig_id")]) ||
    any(!mob_contigs$molecule_type %in%
          c("predicted_plasmid", "predicted_chromosome", "unassigned"))) {
  stop("Plasmid reconstruction release contract failed; claim registry was not published.",
       call. = FALSE)
}

n_to_uti <- sum(transitions$status_from == "Not_UTI" & transitions$status_to == "UTI")
n_to_uti_le25 <- sum(
  transitions$status_from == "Not_UTI" & transitions$status_to == "UTI" &
    transitions$TotalSNPs <= SAME_STRAIN_SNP_THRESHOLD
)
n_all_le25 <- sum(transitions$TotalSNPs <= SAME_STRAIN_SNP_THRESHOLD)
if (n_all_le25 != 140L || n_to_uti != 9L || n_to_uti_le25 != 5L) {
  stop("Direct transition evidence contract failed; claim registry was not published.", call. = FALSE)
}

link_col <- intersect(c("has_vf_pair", "wgs_linked"), names(casebook))
if (!length(link_col)) stop("Mechanism casebook lacks a linkage field.", call. = FALSE)
linked <- as_pipeline_bool(casebook[[link_col[[1]]]])
if (nrow(casebook) != 9L || sum(linked, na.rm = TRUE) != 9L || any(!linked | is.na(linked))) {
  stop("Mechanism casebook contract failed; claim registry was not published.", call. = FALSE)
}
if (nrow(near_miss) != 17L) stop("Near-miss audit contract failed.", call. = FALSE)

qc_config <- release_qc_config()
vfdb_min_identity <- method_number("vfdb_min_identity_pct")
vfdb_min_coverage <- method_number("vfdb_min_coverage_pct")
provider_mlst_min_good_targets <- method_number("provider_mlst_min_good_targets_pct")
provider_mlst_policy <- method_value("provider_st_policy")
if (vfdb_min_identity != 80 || vfdb_min_coverage != 80 ||
    provider_mlst_min_good_targets != 95 ||
    method_number("same_strain_snp_threshold") != SAME_STRAIN_SNP_THRESHOLD) {
  stop("Published method thresholds do not match the implemented Longcycler release contract.", call. = FALSE)
}

gene_cols <- canonical_vf_gene_cols(names(vf), vf_pa_file = vf_path)
typed_mlst <- mlst %>%
  mutate(ST = str_trim(as.character(ST))) %>%
  filter(!is.na(ST), nzchar(ST), !str_to_upper(ST) %in% c("NA", "N/A", "UNKNOWN", "NT"))

event_type <- if ("Event_type" %in% names(cohort)) as.character(cohort$Event_type) else episode_event_type(cohort$tp_lab)
event_rows <- event_type == "UTI_event"
amr_registry_longitudinal <- amr_transitions %>%
  group_by(.data$direct_snp_context) %>%
  summarise(
    pairs = n(),
    gain_or_loss = sum(.data$any_informative_acquired_gene_gain_or_loss),
    median_jaccard = median(.data$acquired_gene_jaccard),
    .groups = "drop"
  )
amr_registry_focused <- amr_focused %>%
  count(.data$focused_strain_context, name = "pairs")

plot_files <- if (dir.exists(DIR_PLOTS)) {
  sort(list.files(DIR_PLOTS, pattern = "\\.(png|pdf|svg)$", recursive = TRUE, full.names = TRUE, ignore.case = TRUE))
} else {
  character()
}
plot_files <- normalise_path(plot_files, must_work = FALSE)
historical_plot <- str_detect(
  plot_files,
  regex("/(?:legacy|archive[^/]*|backup[^/]*)/", ignore_case = TRUE)
) | str_detect(
  basename(plot_files),
  regex("(?:ASB.*UTI|UTI.*ASB)", ignore_case = TRUE)
)
plot_mtime <- file.info(plot_files)$mtime
cohort_mtime <- file.info(cohort_path)$mtime
plot_files <- plot_files[
  !historical_plot & !is.na(plot_mtime) & plot_mtime >= cohort_mtime
]
if (!length(plot_files)) {
  stop("No current-run analytical plot files were available for claim publication.", call. = FALSE)
}
if (any(str_detect(plot_files, regex("flye", ignore_case = TRUE)))) {
  stop("A retired input name remains in the active plot inventory.", call. = FALSE)
}

source_paths <- c(
  cohort = cohort_path,
  manifest = manifest_path,
  vf = vf_path,
  mlst = mlst_path,
  pairs = pair_path,
  transitions = transition_path,
  panaroo_manifest = panaroo_manifest_path,
  panaroo_roary_compatibility = panaroo_roary_path,
  casebook = casebook_path,
  near_miss = near_miss_path,
  rq_marker = rq_marker,
  vf_method_manifest = vf_method_path,
  amr_episode_profiles = amr_profile_path,
  amr_resident_profiles = amr_resident_path,
  amr_adjacent_profiles = amr_transition_path,
  amr_focused_profiles = amr_focused_path,
  amr_longitudinal_inference = amr_inference_path,
  amr_caller_coverage = amr_coverage_path,
  amr_result_summary = amr_table13_path,
  amr_validation = amr_validation_path,
  amr_provenance = amr_provenance_path,
  plasmidfinder_run = pf_run_path,
  plasmidfinder_database = pf_database_path,
  plasmidfinder_gene_matrix = pf_pa_path,
  mob_suite_marker = mob_marker_path,
  mob_suite_run = mob_run_path,
  mob_contig_assignments = mob_contig_path,
  mob_predicted_plasmids = mob_plasmid_path,
  mob_episode_profiles = mob_profile_path,
  plasmid_gene_locations = mob_gene_location_path,
  plasmid_episode_mechanism_profiles = mob_mechanism_profile_path,
  plasmid_adjacent_metrics = mob_adjacent_path,
  plasmid_focused_metrics = mob_focused_path,
  plasmid_location_validation = mob_location_validation_path,
  plasmid_result_summary = plasmid_summary_path,
  project_configuration = file.path(DIR_ROOT, "00_config.R"),
  assembly_qc_implementation = file.path(DIR_ROOT, "R", "wgs_helpers.R"),
  provider_mlst_implementation = file.path(DIR_ROOT, "R", "provider_mlst_integration.R"),
  vfdb_implementation = file.path(DIR_ROOT, "scripts", "research_questions", "run_rq06_08.R"),
  genomic_amr_implementation = file.path(DIR_ROOT, "29_vf_amr_combined_profile.R"),
  plasmidfinder_implementation = file.path(DIR_ROOT, "09_inc_plasmid_network.R"),
  mob_reconstruction_implementation = file.path(
    DIR_ROOT, "09b_mob_plasmid_reconstruction.R"
  ),
  plasmid_localization_implementation = file.path(
    DIR_ROOT, "R", "plasmid_mechanism_helpers.R"
  ),
  direct_pair_implementation = file.path(DIR_ROOT, "11_compare_strains.R"),
  core_genome_implementation = file.path(DIR_ROOT, "12b_core_snp.R"),
  pangenome_implementation = file.path(DIR_ROOT, "12c_panaroo.R")
)
source_paths <- vapply(source_paths, normalise_path, character(1))
sources <- lapply(names(source_paths), function(name) {
  list(role = name, path = source_paths[[name]], sha256 = sha256_file(source_paths[[name]]))
})

registry <- list(
  schema_version = "longcycler_release_claim_registry_v1",
  generated_at = format(Sys.time(), tz = "Europe/Amsterdam", usetz = TRUE),
  analysis_scope = list(
    assembly_policy = "selected QC-passing Longcycler only",
    clinical_phenotype = "operational UTI phenotype",
    clinical_definition_version = UTI_DEFINITION_VERSION,
    interpretation = "exploratory observational analysis; no causal claim"
  ),
  method_contract = list(
    operational_phenotype = list(
      culture_lower_bound_cfu_per_ml = UTI_CFU_THRESHOLD_PRIMARY,
      rule = "versioned operational culture-plus-compatible-symptom phenotype",
      caveat = "not a reconstruction of the full published protocol"
    ),
    assembly_qc = list(
      max_contigs = qc_config$MAX_CONTIGS,
      min_n50_bp = qc_config$MIN_N50,
      min_genome_size_bp = qc_config$MIN_GENOME_SIZE,
      max_genome_size_bp = qc_config$MAX_GENOME_SIZE,
      excluded_metrics = c("read coverage", "completeness", "contamination")
    ),
    vfdb = list(
      tool = "ABRicate",
      database = "VFDB",
      min_identity_pct = vfdb_min_identity,
      min_coverage_pct = vfdb_min_coverage,
      provenance = "SHA-bound calls from the selected Longcycler FASTA manifest"
    ),
    genomic_amr = list(
      role = "supporting mechanism and longitudinal context; not a research question",
      primary_caller = "AMRFinderPlus 4.2.7 acquired genes and known resistance mutations",
      complementary_caller = "ResFinder 4.7.2 with PointFinder",
      legacy_comparison = "ABRicate-ResFinder at 80% identity and coverage",
      amrfinder_database = amr_provenance$version_or_commit[
        amr_provenance$component == "amrfinderplus_database"
      ],
      resfinder_database_commit = amr_provenance$version_or_commit[
        amr_provenance$component == "resfinder_database"
      ],
      pointfinder_database_commit = amr_provenance$version_or_commit[
        amr_provenance$component == "pointfinder_database"
      ],
      mdfA_policy = "retained raw and in sensitivity metrics; excluded from primary acquired-gene burden, gain/loss and Jaccard",
      longitudinal_pairs = nrow(amr_transitions),
      bootstrap_replicates = 10000L,
      bootstrap_seed = 20260712L,
      interpretation = "genomic determinants and predicted phenotypes—not phenotypic AST"
    ),
    plasmids = list(
      role = "secondary assembly-based mechanism context within existing research questions",
      plasmidfinder = list(
        tool = pf_database$caller_version[[1]],
        database = "PlasmidFinder",
        scan_floor_identity_pct = 80L,
        scan_floor_coverage_pct = 60L,
        primary_identity_pct = 80L,
        primary_coverage_pct = 80L,
        feature_key = "exact GENE label; accession retained as metadata"
      ),
      mob_suite = list(
        version = "3.1.9",
        blast_version = "2.15.0",
        mash_version = "2.3",
        database_zenodo_record = "10304948",
        profiled_episodes = nrow(mob_profiles),
        predicted_plasmids = nrow(mob_plasmids),
        interpretation = paste(
          "assembly-based predicted plasmid bins and predicted gene location;",
          "not confirmed circular plasmids, transfer, transmission, or causal evidence"
        )
      )
    ),
    mlst = list(
      role = "lineage context; not pair-specific continuity proof",
      provider_min_good_targets_pct = provider_mlst_min_good_targets,
      provider_policy = provider_mlst_policy,
      fallback = "labelled local MLST from the same selected Longcycler FASTA where required"
    ),
    direct_pair_evidence = list(
      tool = "dnadiff",
      role = "primary pair-specific distance evidence",
      operational_snp_threshold = SAME_STRAIN_SNP_THRESHOLD,
      priority = "graph connectivity and MLST agreement cannot override a conflicting direct pair"
    ),
    population_context = list(
      core_genome_tool = "Parsnp",
      pangenome_tool = "Panaroo",
      role = "population context; not a substitute for direct pair evidence"
    )
  ),
  analytical_cohort = list(
    episodes = nrow(cohort),
    residents = n_distinct(cohort$Participant_id),
    operational_UTI = sum(cohort$UTI_Status == "UTI"),
    operational_Not_UTI = sum(cohort$UTI_Status == "Not_UTI")
  ),
  attrition_qc_context = list(
    label = "full clinical source retained only for attrition/QC context",
    episodes = 583L,
    residents = 166L,
    operational_UTI = 18L,
    operational_Not_UTI = 565L
  ),
  direct_pairs = list(all_within_resident = nrow(pairs)),
  adjacent_transitions = list(
    pairs = nrow(transitions),
    residents = n_distinct(transitions$Participant_id),
    operational_snp_threshold = SAME_STRAIN_SNP_THRESHOLD,
    at_or_below_threshold = n_all_le25,
    Not_UTI_to_UTI = n_to_uti,
    Not_UTI_to_UTI_at_or_below_threshold = n_to_uti_le25
  ),
  mechanism_casebook = list(cases = nrow(casebook), linked = sum(linked), missing = sum(!linked | is.na(linked))),
  near_miss_audit = list(rows = nrow(near_miss), label = "near-miss rows; not operational UTI cases"),
  selected_uti_event_genomes = list(
    genomes = sum(event_rows),
    residents = n_distinct(cohort$Participant_id[event_rows]),
    operational_UTI = sum(event_rows & cohort$UTI_Status == "UTI"),
    operational_Not_UTI = sum(event_rows & cohort$UTI_Status == "Not_UTI")
  ),
  genomic_dimensions = list(
    VFDB_binary_features = length(gene_cols),
    MLST_typed_episodes = nrow(typed_mlst),
    distinct_preferred_ST_labels = n_distinct(typed_mlst$ST),
    genomic_AMR_profiled_episodes = nrow(amr_profiles),
    genomic_AMR_profiled_residents = n_distinct(amr_profiles$Participant_id),
    informative_acquired_AMR_positive_episodes =
      sum(amr_profiles$any_informative_acquired_amr %in% TRUE),
    informative_acquired_AMR_positive_residents =
      sum(amr_residents$any_informative_acquired_amr %in% TRUE),
    known_resistance_mutation_positive_episodes =
      sum(amr_profiles$amr_mutation_count > 0, na.rm = TRUE),
    mdfA_positive_episodes = sum(amr_profiles$mdfA_detected %in% TRUE),
    AMR_adjacent_pairs_with_informative_gain_or_loss =
      sum(amr_transitions$any_informative_acquired_gene_gain_or_loss),
    AMR_focused_Not_UTI_to_UTI_pairs = nrow(amr_focused),
    PlasmidFinder_profiled_episodes = nrow(pf_pa),
    PlasmidFinder_valid_no_hit_episodes =
      sum(as_pipeline_bool(pf_run$no_hit)),
    MOB_profiled_episodes = nrow(mob_profiles),
    MOB_predicted_plasmids = nrow(mob_plasmids),
    predicted_plasmid_localized_gene_hits =
      sum(mob_gene_locations$localization == "predicted_plasmid",
          na.rm = TRUE)
  ),
  genomic_amr_results = list(
    caller_coverage = split(
      amr_coverage[c("n_expected", "n_complete", "n_generated", "n_reused")],
      amr_coverage$caller
    ),
    longitudinal_summary = split(
      amr_registry_longitudinal,
      amr_registry_longitudinal$direct_snp_context
    ),
    resident_cluster_bootstrap = amr_inference,
    focused_context_counts = amr_registry_focused
  ),
  research_questions = list(first = "RQ01", last = "RQ10", count = 10L, retired_questions = 0L),
  plot_files = unname(plot_files),
  sources = sources
)

out_path <- file.path(DIR_RESULTS, "pipeline", "longcycler_release_claim_registry.json")
dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
tmp <- tempfile(pattern = ".claim_registry_", tmpdir = dirname(out_path), fileext = ".json")
write_json(registry, tmp, pretty = TRUE, auto_unbox = TRUE, na = "null")
registry_text <- readLines(tmp, warn = FALSE)
if (any(str_detect(registry_text, regex("flye", ignore_case = TRUE)))) {
  unlink(tmp)
  stop("Claim registry content guard failed.", call. = FALSE)
}
if (!file.rename(tmp, out_path)) {
  unlink(tmp)
  stop("Could not atomically publish claim registry: ", out_path, call. = FALSE)
}
message("Longcycler release claim registry written: ", out_path)

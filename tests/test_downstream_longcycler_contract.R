suppressPackageStartupMessages({
  library(testthat)
  library(dplyr)
  library(readr)
})

project_root <- normalizePath(testthat::test_path(".."), winslash = "/", mustWork = TRUE)

test_that("downstream Longcycler scripts parse", {
  scripts <- c(
    "08_core_vs_plasmid.R", "09_inc_plasmid_network.R",
    "09b_mob_plasmid_reconstruction.R", "10_replicon_heatmap.R",
    "11_compare_strains.R", "11_compare_strains_helpers.R",
    "00c_plot_clinical_summary.R", "15_longitudinal_patterns.R",
    "21_publication_figures.R", "24_vf_longitudinal_dynamics.R",
    "28_vf_transition_case_studies.R", "30_vf_project_summary_tables.R",
    "32_uti_not_uti_diagnostic_stats.R", "33_mechanism_first_addon.R",
    "34_robustness_first_addon.R", "35_final_figure_pack.R",
    "36_statistical_sensitivity_addon.R",
    "scripts/rebuild_longcycler_sensitivity.R",
    "scripts/research_questions/plasmid_mechanism_addon_rq06_08.R",
    "scripts/verify_uti_not_uti_alignment.R",
    "scripts/create_workflow_case_count_flowchart.R",
    "scripts/audit_uti_status_count_explanation.R"
  )
  for (script in scripts) {
    expect_error(parse(file.path(project_root, script)), NA, info = script)
  }
})

test_that("complete runner can safely resume at final verification only", {
  runner <- readLines(file.path(project_root, "RUN_COMPLETE_ANALYSIS.sh"), warn = FALSE)
  expect_true(any(grepl('FINALIZE_ONLY:-0', runner, fixed = TRUE)))
  expect_true(any(grepl('Rscript scripts/verify_longcycler_only_pipeline.R --stage final', runner, fixed = TRUE)))
  expect_true(any(grepl('publish_complete_marker', runner, fixed = TRUE)))
  expect_true(any(grepl('require_completed_genomic_amr', runner, fixed = TRUE)))
  expect_true(any(grepl('REQUIRE_COMPLETED_AMR=1 Rscript scripts/preflight_complete_analysis.R', runner, fixed = TRUE)))
})

test_that("corrected plasmid profile comparison distinguishes empty from failed calls", {
  old_wd <- getwd()
  on.exit(setwd(old_wd), add = TRUE)
  setwd(project_root)
  source(file.path(project_root, "11_compare_strains_helpers.R"), local = TRUE)

  both_empty <- set_profile_metrics(character(), character())
  expect_identical(both_empty$jaccard, 1)
  expect_true(both_empty$both_empty)
  expect_true(both_empty$available)

  failed <- set_profile_metrics(character(), character(), available_a = FALSE)
  expect_true(is.na(failed$jaccard))
  expect_true(is.na(failed$both_empty))
  expect_false(failed$available)

  changed <- set_profile_metrics(c("IncFIB", "IncFIC"), c("IncFIC", "IncI1"))
  expect_equal(changed$jaccard, 1 / 3)
  expect_identical(changed$gains, "IncI1")
  expect_identical(changed$losses, "IncFIB")
})

test_that("legacy ABRicate harmonization collapses alignment duplicates deterministically", {
  old_wd <- getwd()
  on.exit(setwd(old_wd), add = TRUE)
  setwd(project_root)
  source(file.path(project_root, "R", "amr_helpers.R"), local = TRUE)

  raw_path <- tempfile(fileext = ".tsv")
  on.exit(unlink(raw_path), add = TRUE)
  write_tsv(
    tibble(
      GENE = rep("blaTEM-1B_1", 3L),
      RESISTANCE = rep("beta-lactam", 3L),
      PRODUCT = rep("TEM beta-lactamase", 3L),
      `%IDENTITY` = rep(100, 3L),
      `%COVERAGE` = c(100, 100, 99.5),
      SEQUENCE = c("contig_b", "contig_a", "contig_c"),
      START = c(20, 10, 30),
      END = c(900, 890, 910),
      ACCESSION = rep("AF427133", 3L)
    ),
    raw_path
  )
  meta <- tibble(
    Participant_id = "fixture_resident",
    tp_lab = "T1",
    Assembly_ID = "fixture_assembly"
  )
  parsed <- amr_read_abricate(raw_path, meta)

  expect_equal(nrow(read_tsv(raw_path, show_col_types = FALSE)), 3L)
  expect_equal(nrow(parsed), 2L)
  expect_identical(parsed$contig[[1L]], "contig_a")
  expect_setequal(parsed$coverage, c(99.5, 100))
})

test_that("canonical PlasmidFinder and MOB products satisfy exact cohort contracts when complete", {
  plasmid_root <- file.path(project_root, "results", "plasmids")
  pf_marker <- file.path(plasmid_root, "PLASMIDFINDER_RUN_COMPLETE.txt")
  mob_root <- file.path(plasmid_root, "mob_suite")
  mob_marker <- file.path(mob_root, "RUN_COMPLETE.txt")
  skip_if_not(file.exists(pf_marker), "Canonical PlasmidFinder cohort has not completed")

  pf_run <- read_csv(
    file.path(plasmid_root, "plasmidfinder_run_manifest.csv"),
    show_col_types = FALSE
  )
  pf_pa <- read_csv(
    file.path(plasmid_root, "plasmidfinder_presence_absence.csv"),
    show_col_types = FALSE
  )
  pf_catalog <- read_csv(
    file.path(plasmid_root, "plasmidfinder_replicon_catalog.csv"),
    show_col_types = FALSE
  )
  features <- setdiff(names(pf_pa), "Isolate_ID")
  expect_equal(nrow(pf_run), 532L)
  expect_true(all(pf_run$call_status == "complete"))
  expect_equal(sum(pf_run$no_hit %in% TRUE), 110L)
  expect_equal(nrow(pf_pa), 532L)
  expect_equal(length(features), 42L)
  expect_true(all(features %in% pf_catalog$GENE))
  expect_false(any(features %in% pf_catalog$accession))
  expect_true(all(c(
    "IncFIA_1", "IncFIB(AP001918)_1", "IncFIC(FII)_1"
  ) %in% pf_catalog$GENE[pf_catalog$accession == "AP001918"]))

  skip_if_not(file.exists(mob_marker), "MOB-suite cohort has not completed")
  mob_status <- read_csv(file.path(mob_root, "sample_status.csv"),
                         show_col_types = FALSE)
  mob_profiles <- read_csv(file.path(mob_root, "episode_plasmid_profiles.csv"),
                           show_col_types = FALSE)
  mob_contigs <- read_csv(file.path(mob_root, "contig_assignments.csv"),
                          show_col_types = FALSE)
  expect_equal(nrow(mob_status), 532L)
  expect_true(all(mob_status$status == "complete"))
  expect_equal(nrow(mob_profiles), 532L)
  expect_false(anyDuplicated(mob_contigs[c("Isolate_ID", "contig_id")]) > 0)
  expect_true(all(mob_contigs$molecule_type %in%
                    c("predicted_plasmid", "predicted_chromosome", "unassigned")))
})

test_that("plasmid localization retains exact pair denominators and explicit mappings", {
  root <- file.path(project_root, "results", "plasmids", "mob_suite")
  validation_path <- file.path(root, "plasmid_gene_location_validation.csv")
  skip_if_not(file.exists(validation_path), "Plasmid localization has not completed")

  locations <- read_csv(file.path(root, "plasmid_gene_locations_long.csv"),
                        show_col_types = FALSE)
  episodes <- read_csv(file.path(root, "episode_mechanism_profiles.csv"),
                       show_col_types = FALSE)
  adjacent <- read_csv(file.path(root, "adjacent_pair_plasmid_metrics_371.csv"),
                       show_col_types = FALSE)
  focused <- read_csv(file.path(root, "not_uti_to_uti_plasmid_metrics_9.csv"),
                      show_col_types = FALSE)
  validation <- read_csv(validation_path, show_col_types = FALSE)
  expect_equal(nrow(episodes), 532L)
  expect_equal(nrow(adjacent), 371L)
  expect_equal(nrow(focused), 9L)
  expect_false(anyNA(locations$localization))
  expect_true(all(locations$localization %in%
                    c("predicted_plasmid", "chromosome",
                      "ambiguous_or_unassigned")))
  expect_true(all(validation$pass))
})

test_that("canonical clinical cohort is exactly the selected Longcycler analysis set", {
  cohort_path <- file.path(project_root, "results", "clinical", "analysis_cohort_longcycler.csv")
  manifest_path <- file.path(project_root, "results", "qc", "analysis_assembly_manifest.csv")
  skip_if_not(file.exists(cohort_path), "Canonical Longcycler clinical cohort has not been rebuilt")
  skip_if_not(file.exists(manifest_path), "Analysis assembly manifest has not been rebuilt")

  cohort <- read_csv(cohort_path, show_col_types = FALSE) %>%
    mutate(Participant_id = as.character(.data$Participant_id), tp_lab = as.character(.data$tp_lab))
  manifest <- read_csv(manifest_path, show_col_types = FALSE) %>%
    mutate(Participant_id = as.character(.data$Participant_id), tp_lab = as.character(.data$tp_lab))

  expect_equal(nrow(cohort), 532L)
  expect_equal(n_distinct(cohort$Participant_id), 161L)
  expect_equal(sum(cohort$UTI_Status == "UTI"), 16L)
  expect_equal(sum(cohort$UTI_Status == "Not_UTI"), 516L)
  expect_false(anyDuplicated(cohort[c("Participant_id", "tp_lab")]) > 0)
  expect_setequal(
    paste(cohort$Participant_id, cohort$tp_lab, sep = "|"),
    paste(manifest$Participant_id, manifest$tp_lab, sep = "|")
  )
})

test_that("canonical VF-ready export uses exact Longcycler episode IDs and statuses", {
  cohort_path <- file.path(project_root, "results", "clinical", "analysis_cohort_longcycler.csv")
  vf_path <- file.path(project_root, "results", "vf", "vf_analysis_ready.csv")
  skip_if_not(file.exists(cohort_path), "Canonical Longcycler clinical cohort has not been rebuilt")
  skip_if_not(file.exists(vf_path), "Canonical VF-ready export has not been rebuilt")

  cohort <- read_csv(cohort_path, show_col_types = FALSE) %>%
    transmute(
      Participant_id = as.character(.data$Participant_id),
      tp_lab = as.character(.data$tp_lab),
      Episode_ID = as.character(.data$Episode_ID),
      UTI_Status = as.character(.data$UTI_Status)
    ) %>%
    arrange(.data$Participant_id, .data$tp_lab)
  vf_ready <- read_csv(vf_path, show_col_types = FALSE) %>%
    transmute(
      Participant_id = as.character(.data$Participant_id),
      tp_lab = as.character(.data$tp_lab),
      Episode_ID = as.character(.data$Episode_ID),
      UTI_Status = as.character(.data$UTI_Status)
    ) %>%
    arrange(.data$Participant_id, .data$tp_lab)

  expect_identical(vf_ready, cohort)
})

test_that("canonical Longcycler transition exports are complete and directly linked", {
  path <- file.path(project_root, "results", "longitudinal", "longcycler_transitions.csv")
  skip_if_not(file.exists(path), "Canonical Longcycler transition export has not been rebuilt")
  transitions <- read_csv(path, show_col_types = FALSE)

  expect_equal(nrow(transitions), 371L)
  expect_equal(n_distinct(transitions$Participant_id), 139L)
  expect_equal(sum(transitions$status_from == "Not_UTI" & transitions$status_to == "UTI"), 9L)
  expect_true(all(!is.na(transitions$TotalSNPs)))
  expect_true(all(c(
    "pair_key", "TotalSNPs", "AvgIdentity", "Classification", "RuleUsed",
    "snp_strain_context", "st_lineage_context", "pair_interpretation",
    "Assembler_A", "Assembler_B", "Fasta_SHA256_A", "Fasta_SHA256_B",
    "strict_same_strain"
  ) %in% names(transitions)))
  expect_identical(
    transitions$strict_same_strain,
    !is.na(transitions$TotalSNPs) & transitions$TotalSNPs <= 25
  )
  expect_true(all(!is.na(transitions$Fasta_SHA256_A) & nzchar(transitions$Fasta_SHA256_A)))
  expect_true(all(!is.na(transitions$Fasta_SHA256_B) & nzchar(transitions$Fasta_SHA256_B)))
  expect_true(all(tolower(transitions$Assembler_A) == "longcycler"))
  expect_true(all(tolower(transitions$Assembler_B) == "longcycler"))
})

test_that("plasmid hits agree with selected Longcycler path/hash episode keys", {
  hits_path <- file.path(project_root, "results", "plasmids", "plasmidfinder_hits_long.csv")
  manifest_path <- file.path(project_root, "results", "qc", "analysis_assembly_manifest.csv")
  skip_if_not(file.exists(hits_path) && file.exists(manifest_path),
              "Current plasmid hits or Longcycler manifest are unavailable")

  hits <- read_csv(hits_path, show_col_types = FALSE) %>%
    mutate(
      Participant_id = as.character(.data$Participant_id),
      tp_lab = as.character(.data$tp_lab),
      fasta_path = normalizePath(.data$fasta_path, winslash = "/", mustWork = TRUE),
      fasta_sha256 = tolower(as.character(.data$fasta_sha256))
    )
  expect_true(all(c("SEQUENCE", "GENE", "start", "end") %in% names(hits)))
  expect_true(all(is.finite(hits$start) & is.finite(hits$end)))
  expect_true(all(hits$start >= 1 & hits$end >= hits$start))
  manifest <- read_csv(manifest_path, show_col_types = FALSE) %>%
    mutate(
      Participant_id = as.character(.data$Participant_id),
      tp_lab = as.character(.data$tp_lab),
      fasta_path = normalizePath(.data$full_path, winslash = "/", mustWork = TRUE),
      fasta_sha256 = tolower(as.character(.data$fasta_sha256))
    ) %>%
    transmute(fasta_path, fasta_sha256,
              manifest_Participant_id = Participant_id,
              manifest_tp_lab = tp_lab) %>%
    distinct()

  expect_false(anyDuplicated(manifest[c("fasta_path", "fasta_sha256")]) > 0)
  checked <- hits %>%
    left_join(manifest, by = c("fasta_path", "fasta_sha256"), relationship = "many-to-one")
  expect_false(anyNA(checked$manifest_Participant_id))
  expect_identical(checked$Participant_id, checked$manifest_Participant_id)
  expect_identical(checked$tp_lab, checked$manifest_tp_lab)
})

test_that("transition case summaries retain nonmissing VF Jaccard for paired endpoints", {
  cohort_path <- file.path(project_root, "results", "clinical", "analysis_cohort_longcycler.csv")
  path <- file.path(project_root, "results", "vf", "vf_transition_case_summary.csv")
  skip_if_not(file.exists(cohort_path), "Canonical Longcycler cohort has not been rebuilt")
  skip_if_not(file.exists(path), "Transition case summary has not been rebuilt")
  skip_if(file.info(path)$mtime < file.info(cohort_path)$mtime,
          "Transition case summary predates the canonical Longcycler cohort")
  cases <- read_csv(path, show_col_types = FALSE)
  selected <- cases %>% filter(.data$has_vf_pair %in% TRUE)
  expect_true(all(!is.na(selected$vf_jaccard)))
})

test_that("provider-preferred MLST retains or improves usable Longcycler coverage", {
  path <- file.path(project_root, "results", "mlst", "mlst_provider_preferred.csv")
  skip_if_not(file.exists(path), "Provider-preferred MLST has not been rebuilt")

  mlst <- read_csv(path, show_col_types = FALSE, progress = FALSE)
  missing_like <- function(x) {
    y <- trimws(tolower(as.character(x)))
    is.na(x) | y %in% c("", ".", "?", "-", "--", "na", "n/a", "nan", "missing", "unknown", "null")
  }
  source_counts <- mlst %>% count(.data$ST_source, name = "n")
  source_n <- function(source) {
    value <- source_counts$n[match(source, source_counts$ST_source)]
    ifelse(length(value) == 0L || is.na(value), 0L, value)
  }

  provider_n <- source_n("provider_qc95")
  fallback_n <- sum(source_counts$n[grepl("^local_fallback", source_counts$ST_source)])
  missing_n <- source_n("missing") + source_n("missing_provider_conflict")
  preferred_called_n <- sum(!missing_like(mlst$ST))
  local_called_n <- sum(!missing_like(mlst$ST_local))

  expect_equal(nrow(mlst), 532L)
  expect_equal(provider_n, 472L)
  expect_equal(fallback_n, 21L)
  expect_equal(missing_n, 39L)
  expect_equal(preferred_called_n, 493L)
  expect_equal(local_called_n, 481L)
  expect_equal(provider_n + fallback_n + missing_n, nrow(mlst))
  expect_equal(preferred_called_n, provider_n + fallback_n)
  expect_gte(preferred_called_n, local_called_n)
  expect_true(all(tolower(mlst$provider_assembler[mlst$ST_source == "provider_qc95"]) == "longcycler"))
  expect_true(all(mlst$provider_has_classic_7_loci[mlst$ST_source == "provider_qc95"] %in% TRUE))
  dual <- mlst %>%
    filter(.data$ST_source == "provider_qc95", !missing_like(.data$ST_local))
  expect_equal(nrow(dual), 460L)
  expect_equal(sum(as.character(dual$ST_provider) != as.character(dual$ST_local)), 0L)
  expect_true(all(mlst$ST_numeric_comparable_to_local[!missing_like(mlst$ST)] %in% TRUE))
  fallback <- mlst %>% filter(grepl("^local_fallback", .data$ST_source))
  expect_equal(nrow(fallback), 21L)
  expect_true(all(fallback$local_mlst_complete %in% TRUE))
  expect_true(all(!(fallback$local_ambiguous_call %in% TRUE)))
  expect_equal(sum(!missing_like(mlst$ST_provider_below_qc95)), 22L)
  expect_equal(sum(
    grepl("^local_fallback", mlst$ST_source) & !missing_like(mlst$ST_provider_below_qc95)
  ), 21L)
  expect_equal(sum(
    mlst$ST_source %in% c("missing", "missing_provider_conflict") &
      !missing_like(mlst$ST_provider_below_qc95)
  ), 1L)
})

test_that("local MLST QC uses exactly the seven classic loci", {
  locus_path <- file.path(project_root, "results", "mlst", "mlst_locus_list.txt")
  local_path <- file.path(project_root, "results", "mlst", "mlst_with_meta.csv")
  skip_if_not(file.exists(locus_path) && file.exists(local_path), "Local MLST has not been rebuilt")

  loci <- tolower(readLines(locus_path, warn = FALSE))
  local <- read_csv(local_path, show_col_types = FALSE, progress = FALSE)
  expect_setequal(loci, c("adk", "fumc", "gyrb", "icd", "mdh", "pura", "reca"))
  expect_equal(length(loci), 7L)
  expect_true(all(local$n_loci_typed == 7L))
  expect_true(all(local$mlst_complete %in% TRUE))
  expect_equal(sum(local$ambiguous_call %in% TRUE), 8L)
})

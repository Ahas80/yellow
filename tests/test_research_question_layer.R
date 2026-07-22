suppressPackageStartupMessages({
  library(testthat)
  library(dplyr)
  library(readr)
})

project_root <- normalizePath(testthat::test_path(".."), winslash = "/", mustWork = TRUE)

test_that("research-question scripts parse", {
  paths <- file.path(
    project_root,
    "scripts",
    "research_questions",
    c("run_all.R", "run_rq01_05.R", "run_rq06_08.R", "run_rq09_10.R")
  )
  expect_true(all(file.exists(paths)))
  for (path in paths) expect_error(parse(path), NA)
})

test_that("RQ01 and RQ06 both implement the prespecified plasmid change models", {
  rq01 <- paste(
    readLines(
      file.path(
        project_root, "scripts", "research_questions", "run_rq01_05.R"
      ),
      warn = FALSE
    ),
    collapse = "\n"
  )
  rq06 <- paste(
    readLines(
      file.path(
        project_root, "scripts", "research_questions",
        "plasmid_mechanism_addon_rq06_08.R"
      ),
      warn = FALSE
    ),
    collapse = "\n"
  )
  for (text in c(rq01, rq06)) {
    expect_match(text, "any_replicon_profile_change", fixed = TRUE)
    expect_match(text, "any_mob_cluster_change", fixed = TRUE)
    expect_match(text, "ns(days", fixed = TRUE)
    expect_match(text, "c(10L, 25L, 50L)", fixed = TRUE)
    expect_match(text, "exclude_both_empty", fixed = TRUE)
    expect_match(text, "high_confidence_mob_profiles_only", fixed = TRUE)
    expect_match(text, "adjusted_risk_difference", fixed = TRUE)
  }
})

test_that("RQ09 uses resident-label permutation and the gated multiplier interval", {
  rq09 <- paste(
    readLines(
      file.path(
        project_root, "scripts", "research_questions", "run_rq09_10.R"
      ),
      warn = FALSE
    ),
    collapse = "\n"
  )
  expect_match(rq09, 'primary_inference = "resident-label permutation"',
               fixed = TRUE)
  expect_match(rq09, 'rexp(length(bootstrap_residents), rate = 1)',
               fixed = TRUE)
  expect_match(rq09, "9990L", fixed = TRUE)
  expect_match(
    rq09,
    "resident_exponential_multiplier_bootstrap_distribution.csv",
    fixed = TRUE
  )
  expect_false(grepl(
    "resident_cluster_bootstrap_distribution.csv", rq09, fixed = TRUE
  ))
})

test_that("research-question scripts do not read from the excluded Rowena tree", {
  paths <- list.files(
    file.path(project_root, "scripts", "research_questions"),
    pattern = "\\.R$",
    full.names = TRUE
  )
  for (path in paths) {
    lines <- readLines(path, warn = FALSE)
    prohibited <- grepl(
      "(read_|read\\.|source\\(|load\\(|file\\.path\\(|here::here\\().*Rowenas analysis",
      lines,
      ignore.case = TRUE
    )
    expect_false(any(prohibited), info = paste("Excluded input reference in", path))
  }
})

test_that("release runner fixes publication resampling contracts", {
  run_all <- readLines(
    file.path(project_root, "scripts", "research_questions", "run_all.R"),
    warn = FALSE
  )
  rq06 <- readLines(
    file.path(project_root, "scripts", "research_questions", "run_rq06_08.R"),
    warn = FALSE
  )
  rq_tail <- readLines(
    file.path(project_root, "scripts", "research_questions", "run_rq09_10.R"),
    warn = FALSE
  )
  expect_true(any(grepl('RQ_BOOTSTRAP_REPS = "10000"', run_all, fixed = TRUE)))
  expect_true(any(grepl('RQ_PERMUTATIONS = "10000"', run_all, fixed = TRUE)))
  expect_false(any(grepl("RQ_BOOT_REPS", rq06, fixed = TRUE)))
  expect_true(any(grepl('c("fresh", "generated", "reused")', rq_tail, fixed = TRUE)))
  expect_true(any(grepl("validate_dnadiff_sidecar", rq_tail, fixed = TRUE)))
  expect_false(any(grepl("provider_assembler", rq_tail, fixed = TRUE)))
})

test_that("RQ06-RQ08 VFDB pmap inputs match the reader callback", {
  rq06 <- readLines(
    file.path(project_root, "scripts", "research_questions", "run_rq06_08.R"),
    warn = FALSE
  )
  expected_mappings <- c(
    "tab_path = cache_file",
    "episode_key_value = episode_key",
    "pid = Participant_id",
    "tp = tp_lab",
    "assembly_id = Assembly_ID",
    "fasta_sha = fasta_sha256"
  )
  for (mapping in expected_mappings) {
    expect_true(any(grepl(mapping, rq06, fixed = TRUE)), info = mapping)
  }
})

test_that("RQ06 binary inference keeps the scalar exposure name out of mutate scope", {
  rq_text <- paste(
    readLines(file.path(project_root, "scripts", "research_questions", "run_rq06_08.R")),
    collapse = "\n"
  )
  expect_match(rq_text, "exposure_name <- exposure", fixed = TRUE)
  expect_match(rq_text, "n_exposed = sum(df[[exposure_name]] %in% TRUE", fixed = TRUE)
  expect_false(grepl("n_exposed = sum(df[[exposure]] %in% TRUE", rq_text, fixed = TRUE))
})

test_that("release reruns clear prior question inventories but retain validated caches", {
  run_all <- readLines(
    file.path(project_root, "scripts", "research_questions", "run_all.R"),
    warn = FALSE
  )
  expect_true(any(grepl('sprintf("RQ%02d", 1:10)', run_all, fixed = TRUE)))
  expect_true(any(grepl('c("vfdb_cache_sha256_v1", "rq09_mash")', run_all, fixed = TRUE)))
  expect_true(any(grepl("prior_shared_inputs", run_all, fixed = TRUE)))
  expect_true(any(grepl('pattern = "(rq11|rq09_11)"', run_all, fixed = TRUE)))
})

test_that("generated-content cleanup and final gate scan ignored files fail-closed", {
  paths <- c(
    file.path(project_root, "scripts", "prepare_longcycler_release.R"),
    file.path(project_root, "scripts", "verify_longcycler_only_pipeline.R")
  )
  for (path in paths) {
    lines <- readLines(path, warn = FALSE)
    expect_true(any(grepl('"--no-ignore"', lines, fixed = TRUE)), info = path)
    expect_true(any(grepl("shQuote(glob)", lines, fixed = TRUE)), info = path)
    expect_true(any(grepl("rg_status", lines, fixed = TRUE)), info = path)
    expect_true(any(grepl("forbidden_content_pattern", lines, fixed = TRUE)), info = path)
  }
})

test_that("final verifier enumerates RQ01-RQ10 statuses and accepts the published marker scope", {
  verifier <- readLines(
    file.path(project_root, "scripts", "verify_longcycler_only_pipeline.R"),
    warn = FALSE
  )
  expect_true(any(grepl(
    'file.path(rq_root, sprintf("RQ%02d", 1:10), "analysis_status.csv")',
    verifier,
    fixed = TRUE
  )))
  expect_true(any(grepl('fixed("RQ01-RQ10")', verifier, fixed = TRUE)))
  expect_true(any(grepl('fixed("RQ01--RQ10")', verifier, fixed = TRUE)))
})

test_that("inactive sensitivity results and binary VF inputs cannot survive release", {
  cleanup <- readLines(
    file.path(project_root, "scripts", "prepare_longcycler_release.R"),
    warn = FALSE
  )
  verifier <- readLines(
    file.path(project_root, "scripts", "verify_longcycler_only_pipeline.R"),
    warn = FALSE
  )
  expect_true(any(grepl('file.path("results", "sensitivity")', cleanup, fixed = TRUE)))
  expect_true(any(grepl('file.path(DIR_RESULTS, "sensitivity")', verifier, fixed = TRUE)))
  expect_true(any(grepl('require_file(FILE_VF_HITS, "VF hit RDS")', verifier, fixed = TRUE)))
  expect_true(any(grepl('"paths_outside_manifest"', verifier, fixed = TRUE)))
  expect_true(any(grepl('"not_older_than_selected_cohort"', verifier, fixed = TRUE)))
})

test_that("diagnostic reruns refresh the selected-cohort count audit", {
  diagnostics <- readLines(
    file.path(project_root, "32_uti_not_uti_diagnostic_stats.R"),
    warn = FALSE
  )
  verifier <- readLines(
    file.path(project_root, "scripts", "verify_longcycler_only_pipeline.R"),
    warn = FALSE
  )
  expect_true(any(grepl('"audit_uti_status_count_explanation.R"', diagnostics, fixed = TRUE)))
  expect_true(any(grepl('"Selected-cohort UTI/Not_UTI count audit"', verifier, fixed = TRUE)))
})

test_that("the first downstream Panaroo consumer enforces a current cleanup sweep", {
  selection_plot <- readLines(
    file.path(project_root, "13_visualise_panaroo_selection.R"),
    warn = FALSE
  )
  expect_true(any(grepl('"prepare_longcycler_release.R"', selection_plot, fixed = TRUE)))
  expect_true(any(grepl("cleanup_summary_mtime >= selected_cohort_mtime", selection_plot, fixed = TRUE)))
  expect_true(any(grepl('c(cleanup_script, "--apply")', selection_plot, fixed = TRUE)))
})

test_that("Panaroo compatibility output is manifest-bound and uses GFF sample ids", {
  panaroo <- readLines(file.path(project_root, "12c_panaroo.R"), warn = FALSE)
  mechanism <- readLines(file.path(project_root, "33_mechanism_first_addon.R"), warn = FALSE)
  verifier <- readLines(
    file.path(project_root, "scripts", "verify_longcycler_only_pipeline.R"),
    warn = FALSE
  )
  expect_true(any(grepl('"gene_presence_absence_roary.csv"', panaroo, fixed = TRUE)))
  expect_true(any(grepl('sprintf("Exit status: %d", res)', panaroo, fixed = TRUE)))
  expect_true(any(grepl('"Exit status: 0"', panaroo, fixed = TRUE)))
  expect_true(any(grepl("panaroo_sample_id = str_remove", mechanism, fixed = TRUE)))
  expect_true(any(grepl("basename(as.character(gff_path))", mechanism, fixed = TRUE)))
  expect_true(any(grepl('"roary_sample_columns_outside_manifest"', verifier, fixed = TRUE)))
  expect_true(any(grepl('"panaroo_roary_compatibility"', verifier, fixed = TRUE)))
  expect_true(any(grepl('"source_sha256_mismatches"', verifier, fixed = TRUE)))
})

test_that("source attrition context and selected analytical denominators are exact", {
  status <- read_csv(
    file.path(project_root, "results", "clinical", "status_map.csv"),
    show_col_types = FALSE
  ) %>%
    filter(.data$analysis_include_primary %in% TRUE) %>%
    mutate(
      Participant_id = as.character(.data$Participant_id),
      tp_lab = as.character(.data$tp_lab)
    )
  manifest <- read_csv(
    file.path(project_root, "results", "qc", "analysis_assembly_manifest.csv"),
    show_col_types = FALSE
  )
  analytical <- status %>%
    semi_join(
      manifest %>% transmute(
        Participant_id = as.character(.data$Participant_id),
        tp_lab = as.character(.data$tp_lab)
      ),
      by = c("Participant_id", "tp_lab")
    )
  transition_path <- file.path(project_root, "results", "longitudinal", "longcycler_transitions.csv")
  skip_if_not(file.exists(transition_path), "Canonical transition table has not been rebuilt yet")
  transitions <- read_csv(transition_path, show_col_types = FALSE)

  expect_equal(nrow(status), 583)
  expect_equal(n_distinct(status$Participant_id), 166)
  expect_equal(sum(status$UTI_Status == "UTI"), 18)
  expect_equal(sum(status$UTI_Status == "Not_UTI"), 565)

  expect_equal(nrow(manifest), 532)
  expect_equal(n_distinct(manifest$Participant_id), 161)
  expect_false(anyDuplicated(manifest[c("Participant_id", "tp_lab")]) > 0)
  assembler <- tolower(if ("assembler" %in% names(manifest)) manifest$assembler else manifest$Assembler)
  expect_true(all(assembler == "longcycler"))
  expect_true(all(manifest$QC_PASS %in% TRUE))

  expect_equal(nrow(analytical), 532)
  expect_equal(n_distinct(analytical$Participant_id), 161)
  expect_equal(sum(analytical$UTI_Status == "UTI"), 16)
  expect_equal(sum(analytical$UTI_Status == "Not_UTI"), 516)

  expect_equal(nrow(transitions), 371)
  expect_equal(n_distinct(transitions$Participant_id), 139)
  expect_true(all(!is.na(transitions$TotalSNPs)))
  expect_equal(sum(transitions$TotalSNPs <= 25), 140)
  expect_equal(sum(transitions$status_from == "Not_UTI" & transitions$status_to == "UTI"), 9)
  expect_equal(sum(
    transitions$status_from == "Not_UTI" & transitions$status_to == "UTI" &
      transitions$TotalSNPs <= 25
  ), 5)

  expected_all_pairs <- manifest %>%
    count(.data$Participant_id, name = "n") %>%
    summarise(n_pairs = sum(.data$n * (.data$n - 1) / 2)) %>%
    pull(.data$n_pairs)
  expect_equal(expected_all_pairs, 893)

  cohort_path <- file.path(project_root, "results", "clinical", "analysis_cohort_longcycler.csv")
  if (file.exists(cohort_path)) {
    cohort <- read_csv(cohort_path, show_col_types = FALSE)
    expect_equal(nrow(cohort), 532)
    expect_setequal(
      paste(cohort$Participant_id, cohort$tp_lab, sep = "||"),
      paste(manifest$Participant_id, manifest$tp_lab, sep = "||")
    )
  }
})

test_that("every current transition has direct pair evidence", {
  transition_path <- file.path(project_root, "results", "longitudinal", "longcycler_transitions.csv")
  skip_if_not(file.exists(transition_path), "Canonical transition table has not been rebuilt yet")
  transitions <- read_csv(transition_path, show_col_types = FALSE)
  expect_true(all(!is.na(transitions$pair_key) & nzchar(transitions$pair_key)))
  expect_true(all(!is.na(transitions$Fasta_SHA256_A) & nzchar(transitions$Fasta_SHA256_A)))
  expect_true(all(!is.na(transitions$Fasta_SHA256_B) & nzchar(transitions$Fasta_SHA256_B)))
  expect_true(all(tolower(transitions$Assembler_A) == "longcycler"))
  expect_true(all(tolower(transitions$Assembler_B) == "longcycler"))
  expect_true(all(transitions$dnadiff_cache_status %in% c("fresh", "generated", "reused")))
  expect_true(all(file.exists(transitions$dnadiff_sidecar_path)))
})

test_that("final research-question release satisfies all anchor contracts", {
  out_root <- file.path(project_root, "results", "research_questions")
  release_marker <- file.path(out_root, "RUN_COMPLETE.txt")
  skip_if_not(file.exists(release_marker), "Final RQ01--RQ10 release has not been run yet")

  checks <- read_csv(file.path(out_root, "final_contract_checks.csv"), show_col_types = FALSE)
  expect_true(all(checks$pass))
  expected_checks <- c(
    "source_attrition_clinical_episodes", "selected_longcycler_genomes",
    "selected_operational_uti_genomes", "selected_operational_not_uti_genomes",
    "adjacent_longcycler_pairs", "all_direct_within_resident_pairs",
    "selected_uti_event_genomes"
  )
  expect_true(all(expected_checks %in% checks$check))

  status <- read_csv(file.path(out_root, "final_question_status.csv"), show_col_types = FALSE)
  expect_equal(nrow(status), 10L)
  expect_true(all(status$status == "complete"))
})

test_that("completed RQ release contains no excluded-assembler content", {
  out_root <- file.path(project_root, "results", "research_questions")
  release_marker <- file.path(out_root, "RUN_COMPLETE.txt")
  skip_if_not(file.exists(release_marker), "Final RQ01--RQ10 release has not been run yet")

  forbidden <- paste0("(^|[^[:alnum:]])", "fl", "ye", "([^[:alnum:]]|$)")
  paths <- list.files(out_root, recursive = TRUE, full.names = TRUE, all.files = TRUE)
  expect_false(any(grepl(paste0("fl", "ye"), paths, ignore.case = TRUE)))

  text_paths <- paths[file.exists(paths) & !dir.exists(paths) &
    grepl("\\.(csv|tsv|txt|md|json)$", paths, ignore.case = TRUE)]
  violations <- vapply(text_paths, function(path) {
    any(grepl(forbidden, readLines(path, warn = FALSE), ignore.case = TRUE))
  }, logical(1))
  expect_false(any(violations), info = paste(text_paths[violations], collapse = ", "))
})

test_that("publication-facing deidentified tables contain no raw key columns", {
  out_root <- file.path(project_root, "results", "research_questions")
  tables <- list.files(
    out_root,
    pattern = "(deidentified|case_matrix|case_table).*\\.csv$",
    recursive = TRUE,
    full.names = TRUE,
    ignore.case = TRUE
  )
  prohibited <- c("Participant_id", "tp_lab", "Episode_ID", "Isolate_ID", "Assembly_ID")
  for (path in tables) {
    x <- read_csv(path, show_col_types = FALSE)
    expect_false(any(names(x) %in% prohibited), info = path)
  }
})

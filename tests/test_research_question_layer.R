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
    c("run_all.R", "run_rq01_05.R", "run_rq06_08.R", "run_rq09_11.R")
  )
  expect_true(all(file.exists(paths)))
  for (path in paths) expect_error(parse(path), NA)
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
  expect_true(any(grepl('RQ_BOOTSTRAP_REPS = "10000"', run_all, fixed = TRUE)))
  expect_true(any(grepl('RQ_PERMUTATIONS = "10000"', run_all, fixed = TRUE)))
  expect_false(any(grepl("RQ_BOOT_REPS", rq06, fixed = TRUE)))
})

test_that("authoritative clinical and genomic denominators are unchanged", {
  status <- read_csv(
    file.path(project_root, "results", "clinical", "status_map.csv"),
    show_col_types = FALSE
  ) %>% filter(.data$analysis_include_primary %in% TRUE)
  manifest <- read_csv(
    file.path(project_root, "results", "qc", "analysis_assembly_manifest.csv"),
    show_col_types = FALSE
  )
  transitions <- read_csv(
    file.path(project_root, "results", "sensitivity", "longcycler_only", "longcycler_transitions.csv"),
    show_col_types = FALSE
  )

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

  expect_equal(nrow(transitions), 371)
  expect_equal(n_distinct(transitions$Participant_id), 139)
  expect_true(all(!is.na(transitions$TotalSNPs)))
  expect_equal(sum(transitions$TotalSNPs <= 25), 140)
})

test_that("every current transition has direct pair evidence", {
  transitions <- read_csv(
    file.path(project_root, "results", "sensitivity", "longcycler_only", "longcycler_transitions.csv"),
    show_col_types = FALSE
  )
  expect_true(all(!is.na(transitions$pair_key) & nzchar(transitions$pair_key)))
  expect_true(all(!is.na(transitions$Fasta_SHA256_A) & nzchar(transitions$Fasta_SHA256_A)))
  expect_true(all(!is.na(transitions$Fasta_SHA256_B) & nzchar(transitions$Fasta_SHA256_B)))
  expect_true(all(tolower(transitions$Assembler_A) == "longcycler"))
  expect_true(all(tolower(transitions$Assembler_B) == "longcycler"))
})

test_that("final research-question release satisfies all anchor contracts", {
  out_root <- file.path(project_root, "results", "research_questions")
  release_marker <- file.path(out_root, "RUN_COMPLETE.txt")
  skip_if_not(file.exists(release_marker), "Final RQ01--RQ11 release has not been run yet")

  checks <- read_csv(file.path(out_root, "final_contract_checks.csv"), show_col_types = FALSE)
  expect_true(all(checks$pass))
  expected_checks <- c(
    "eligible_clinical_episodes", "selected_longcycler_genomes",
    "adjacent_longcycler_pairs", "all_direct_within_resident_pairs",
    "selected_uti_event_genomes", "rq11_paired_episode_keys"
  )
  expect_true(all(expected_checks %in% checks$check))

  status <- read_csv(file.path(out_root, "final_question_status.csv"), show_col_types = FALSE)
  expect_equal(nrow(status), 11L)
  expect_true(all(status$status == "complete"))
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

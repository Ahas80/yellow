library(testthat)

withr::local_dir(here::here())
source(here::here("00_config.R"))

write_selection_fasta <- function(path) {
  writeLines(c(">contig", "ACGTACGT"), path)
  normalizePath(path, winslash = "/", mustWork = TRUE)
}

test_that("Flye is never selected when Longcycler passes", {
  td <- tempfile("longcycler-selection-")
  dir.create(td)
  lc <- write_selection_fasta(file.path(td, "episode-longcycler.fasta"))
  flye <- write_selection_fasta(file.path(td, "episode-flye.fasta"))
  candidates <- tibble::tibble(
    Participant_id = c("1", "1"),
    tp_lab = c("T0", "T0"),
    Assembly_ID = c("episode-longcycler", "episode-flye"),
    Assembler = c("longcycler", "flye"),
    full_path = c(lc, flye),
    file_name = basename(c(lc, flye)),
    QC_PASS = c(TRUE, TRUE),
    N50 = c(100, 200),
    n_contigs = c(2, 1)
  )

  selected <- select_canonical_assemblies(candidates)
  expect_equal(selected$Assembler[selected$selected_canonical], "longcycler")
  expect_false(any(selected$selected_canonical & selected$Assembler == "flye"))
})

test_that("a passing Flye is not used when Longcycler fails", {
  td <- tempfile("longcycler-no-fallback-")
  dir.create(td)
  lc <- write_selection_fasta(file.path(td, "episode-longcycler.fasta"))
  flye <- write_selection_fasta(file.path(td, "episode-flye.fasta"))
  candidates <- tibble::tibble(
    Participant_id = c("1", "1"),
    tp_lab = c("T0", "T0"),
    Assembly_ID = c("episode-longcycler", "episode-flye"),
    Assembler = c("longcycler", "flye"),
    full_path = c(lc, flye),
    file_name = basename(c(lc, flye)),
    QC_PASS = c(FALSE, TRUE),
    N50 = c(100, 200),
    n_contigs = c(2, 1)
  )

  selected <- select_canonical_assemblies(candidates)
  expect_false(any(selected$selected_canonical))
  expect_match(
    selected$canonical_reason[selected$Assembler == "flye"],
    "excluded: assembler not permitted"
  )
})

test_that("GFF lookup requires an exact selected assembly key", {
  td <- tempfile("strict-gff-key-")
  dir.create(td)
  expected <- tibble::tibble(
    Assembly_ID = "episode-longcycler__fasta",
    Assembly_Base_ID = "episode-longcycler"
  )

  substring_only <- file.path(td, "prefix-episode-longcycler-extra.gff")
  writeLines("##gff-version 3", substring_only)
  unresolved <- attach_pipeline_gff_paths(expected, gff_dirs = td)
  expect_false(unresolved$gff_available)
  expect_true(is.na(unresolved$gff_path))

  exact <- file.path(td, "episode-longcycler.gff")
  writeLines("##gff-version 3", exact)
  resolved <- attach_pipeline_gff_paths(expected, gff_dirs = td)
  expect_true(resolved$gff_available)
  expect_equal(resolved$gff_path, normalizePath(exact, winslash = "/", mustWork = TRUE))
})

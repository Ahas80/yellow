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

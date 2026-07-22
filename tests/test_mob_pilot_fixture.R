suppressPackageStartupMessages({
  library(testthat)
  library(readr)
  library(stringr)
})

pilot_dir <- Sys.getenv(
  "MOB_PILOT_DIR",
  "/tmp/codex_mob_pilot_2510C119001-1-debug"
)

test_that("validated MOB-suite pilot retains its raw reconstruction contract", {
  contig_path <- file.path(pilot_dir, "contig_report.txt")
  typer_path <- file.path(pilot_dir, "mobtyper_results.txt")
  skip_if_not(
    file.exists(contig_path) && file.exists(typer_path),
    "Set MOB_PILOT_DIR to the retained 2510C119001-1 pilot output"
  )

  contigs <- read_tsv(
    contig_path, show_col_types = FALSE, col_types = cols(.default = "c")
  )
  typer <- read_tsv(
    typer_path, show_col_types = FALSE, col_types = cols(.default = "c")
  )
  expect_equal(nrow(contigs), 9L)
  expect_equal(sum(contigs$molecule_type == "chromosome"), 3L)
  expect_equal(sum(contigs$molecule_type == "plasmid"), 6L)
  expect_false(anyDuplicated(contigs$contig_id) > 0)

  incf <- typer[typer$primary_cluster_id == "AA179", , drop = FALSE]
  expect_equal(nrow(incf), 1L)
  expect_equal(as.numeric(incf$size), 147696)
  incf_types <- str_split(incf$`rep_type(s)`, ",", simplify = TRUE)
  expect_true(all(c("IncFIB", "IncFIC", "IncFII") %in% incf_types))
  expect_identical(incf$predicted_mobility, "conjugative")
})

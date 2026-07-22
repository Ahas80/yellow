suppressPackageStartupMessages({
  library(testthat)
  library(digest)
  library(dplyr)
})

project_root <- normalizePath(testthat::test_path(".."), winslash = "/",
                              mustWork = TRUE)
source(file.path(project_root, "R", "amr_helpers.R"))

test_that("AMR normalization preserves alleles and flags mdf(A)", {
  expect_identical(amr_normalise_symbol("blaTEM-1_1"), "blaTEM-1")
  expect_identical(amr_normalise_symbol("mdfA"), "mdf(A)")
  expect_identical(amr_gene_family("blaTEM-1"), "blaTEM")
  expect_true(amr_is_background("mdf(A)"))
  expect_false(amr_is_background("tet(A)"))
  expect_identical(amr_normalise_mutation("gyrA_S83L"), "gyrA:S83L")
  expect_identical(
    amr_infer_class("blaTEM-1", "BETA-LACTAM"), "Beta-lactams"
  )
})

test_that("mdf(A) inflates similarity without informative overlap", {
  including <- amr_jaccard(
    c("mdf(A)", "blaTEM"), c("mdf(A)", "tet(A)")
  )
  informative <- amr_jaccard(c("blaTEM"), c("tet(A)"))
  expect_equal(including, 1 / 3)
  expect_equal(informative, 0)
})

test_that("cache keys bind hashes versions databases and parameters", {
  base <- amr_cache_key(
    "fasta", "annotation", "caller", "1.0", "database",
    list(identity = 80, coverage = 80)
  )
  expect_identical(
    base,
    amr_cache_key(
      "fasta", "annotation", "caller", "1.0", "database",
      list(identity = 80, coverage = 80)
    )
  )
  expect_false(identical(
    base,
    amr_cache_key(
      "changed", "annotation", "caller", "1.0", "database",
      list(identity = 80, coverage = 80)
    )
  ))
  expect_false(identical(
    base,
    amr_cache_key(
      "fasta", "annotation", "caller", "1.0", "database",
      list(identity = 90, coverage = 80)
    )
  ))
})

test_that("verified AMRFinderPlus 4.2.7 schema is parsed", {
  path <- tempfile(fileext = ".tsv")
  writeLines(
    paste(
      c(
        "Name", "Protein id", "Contig id", "Start", "Stop", "Strand",
        "Element symbol", "Element name", "Scope", "Type", "Subtype",
        "Class", "Subclass", "Method", "Target length",
        "Reference sequence length", "% Coverage of reference",
        "% Identity to reference", "Alignment length",
        "Closest reference accession", "Closest reference name",
        "HMM accession", "HMM description"
      ),
      collapse = "\t"
    ),
    path
  )
  cat(
    paste(
      c(
        "sample", "protein", "contig", "1", "100", "+", "gyrA_S83L",
        "quinolone-resistant GyrA", "core", "AMR", "POINT",
        "FLUOROQUINOLONE", "FLUOROQUINOLONE", "POINTP", "100", "100",
        "100.0", "99.0", "100", "ACC", "GyrA", "NA", "NA"
      ),
      collapse = "\t"
    ),
    "\n", file = path, append = TRUE, sep = ""
  )
  meta <- tibble(
    Participant_id = "P1", tp_lab = "T1", Assembly_ID = "assembly"
  )
  parsed <- amr_read_amrfinder(path, meta)
  expect_equal(nrow(parsed), 1L)
  expect_identical(parsed$normalized_symbol, "gyrA:S83L")
  expect_identical(parsed$gene_family, "gyrA")
  expect_identical(parsed$determinant_type, "point_mutation")
  expect_identical(parsed$drug_class, "Fluoroquinolones")
})

test_that("ResFinder E. coli genomic-phenotype table is comment-aware", {
  directory <- tempfile()
  dir.create(directory)
  path <- file.path(directory, "pheno_table_escherichia_coli.txt")
  writeLines(
    c(
      "# ResFinder phenotype results for escherichia coli.",
      "#",
      paste(
        "# Antimicrobial", "Class", "WGS-predicted phenotype", "Match",
        "Genetic background", sep = "\t"
      ),
      paste("ciprofloxacin", "quinolone", "Resistant", "3", "gyrA:S83L",
            sep = "\t")
    ),
    path
  )
  meta <- tibble(
    Participant_id = "P1", tp_lab = "T1", Assembly_ID = "assembly"
  )
  parsed <- amr_read_predicted_phenotype(directory, meta)
  expect_equal(nrow(parsed), 1L)
  expect_identical(parsed$antimicrobial, "ciprofloxacin")
  expect_identical(parsed$wgs_predicted_phenotype, "Resistant")
  expect_identical(parsed$match_level, 3L)
  expect_identical(
    parsed$interpretation_scope,
    "genomic prediction—not phenotypic AST"
  )
})

test_that("cache reuse requires a matching successful completion marker", {
  directory <- tempfile()
  dir.create(directory)
  output <- file.path(directory, "output.tsv")
  stderr <- file.path(directory, "stderr.txt")
  marker <- file.path(directory, "complete.txt")
  writeLines("raw", output)
  writeLines("", stderr)
  amr_write_completion_marker(
    marker, "key",
    list(
      exit_status = 0L, started_at = "start", completed_at = "end",
      command = "caller"
    ),
    "caller", output
  )
  expect_true(amr_cache_complete(marker, "key", c(output, stderr)))
  expect_false(amr_cache_complete(marker, "other-key", c(output, stderr)))
})

test_that("script 29 and downstream integration parse", {
  paths <- file.path(
    project_root,
    c(
      "29_vf_amr_combined_profile.R", "30_vf_project_summary_tables.R",
      "33_mechanism_first_addon.R", "scripts/nursing_home_candidate_clusters.R",
      "scripts/verify_longcycler_only_pipeline.R"
    )
  )
  for (path in paths) expect_error(parse(path), NA, info = path)
  downstream_text <- paste(
    unlist(lapply(
      file.path(
        project_root,
        c(
          "33_mechanism_first_addon.R",
          "scripts/nursing_home_candidate_clusters.R"
        )
      ),
      readLines, warn = FALSE
    )),
    collapse = "\n"
  )
  expect_false(grepl("abricate_resfinder_cache", downstream_text, fixed = TRUE))
  expect_true(grepl(
    "results.*amr|DIR_RESULTS.*amr", downstream_text, perl = TRUE
  ))
})

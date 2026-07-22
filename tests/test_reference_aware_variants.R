suppressPackageStartupMessages({
  library(testthat)
  library(tibble)
})

project_root <- normalizePath(testthat::test_path(".."), winslash = "/", mustWork = TRUE)
source(file.path(project_root, "scripts", "prepare_reference_aware_variants.R"), local = FALSE)

test_that("sequence ID normalization is conservative and deterministic", {
  expect_equal(
    normalise_sequence_id(c(">contig_1 description", " contig_2\tmetadata ", "contig_3")),
    c("contig_1", "contig_2", "contig_3")
  )
  expect_error(
    contig_evidence_from_sequences(c("same one", "same two"), c("ACGT", "TGCA"), "synthetic"),
    "duplicate normalized contig IDs"
  )
})

test_that("renamed GFF contigs require a unique complete-sequence hash", {
  reference <- contig_evidence_from_sequences(
    c("reference_A", "reference_B"), c("AACCGGTT", "TTTTCCCC"), "reference"
  )
  renamed <- contig_evidence_from_sequences(
    c("gnl|X|A", "gnl|X|B"), c("AACCGGTT", "TTTTCCCC"), "gff"
  )
  mapping <- resolve_reference_to_gff_contigs(reference, renamed)
  expect_true(all(mapping$Contig_Mapping_Valid))
  expect_equal(mapping$Contig_Mapping_Method, rep("unique_sequence_sha256", 2))
  expect_equal(mapping$GFF_Contig_ID, c("gnl|X|A", "gnl|X|B"))
})

test_that("same length alone is never accepted as contig identity", {
  reference <- contig_evidence_from_sequences("reference_A", "AACCGGTT", "reference")
  gff <- contig_evidence_from_sequences("renamed_A", "TTGGCCAA", "gff")
  mapping <- resolve_reference_to_gff_contigs(reference, gff)
  expect_false(mapping$Contig_Mapping_Valid)
  expect_equal(mapping$Contig_Mapping_Method, "unmapped_no_sequence_identity")
  expect_equal(mapping$Length_Only_Candidate_Count, 1)
})

test_that("non-unique sequence matches fail closed", {
  reference <- contig_evidence_from_sequences("reference_A", "AACCGGTT", "reference")
  gff <- contig_evidence_from_sequences(
    c("renamed_A", "renamed_B"), c("AACCGGTT", "AACCGGTT"), "gff"
  )
  mapping <- resolve_reference_to_gff_contigs(reference, gff)
  expect_false(mapping$Contig_Mapping_Valid)
  expect_equal(mapping$Contig_Mapping_Method, "ambiguous_sequence_sha256")
})

test_that("show-snps keys preserve contig-local identity", {
  first <- variant_key(10, "A", "G", 12, "ref_1", "qry_1")
  second <- variant_key(10, "A", "G", 12, "ref_2", "qry_1")
  expect_false(identical(first, second))
})


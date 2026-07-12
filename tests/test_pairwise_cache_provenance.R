library(testthat)

withr::local_dir(here::here())
source(here::here("11_compare_strains_helpers.R"))

write_test_fasta <- function(path, sequence) {
  writeLines(c(">test", sequence), path)
  normalizePath(path, winslash = "/", mustWork = TRUE)
}

write_matching_test_sidecar <- function(spec) {
  payload <- list(
    schema_version = spec$schema_version,
    pair_key = spec$pair_key,
    cache_signature = spec$cache_signature,
    input_a = spec$input_a,
    input_b = spec$input_b,
    report = fasta_fingerprint(spec$report_path),
    dnadiff = list(
      executable = "test-dnadiff",
      version = "test-version",
      command = "test command",
      generated_at_utc = "2026-01-01T00:00:00Z"
    )
  )
  jsonlite::write_json(payload, spec$sidecar_path, auto_unbox = TRUE, pretty = TRUE)
  invisible(payload)
}

test_that("cache signature changes when the FASTA selected for the same sample key changes", {
  td <- tempfile("pairwise-cache-")
  dir.create(td)
  a <- write_test_fasta(file.path(td, "sample-longcycler.fasta"), paste(rep("ACGT", 50), collapse = ""))
  b <- write_test_fasta(file.path(td, "sample-flye.fasta"), paste(rep("ACGA", 50), collapse = ""))
  other <- write_test_fasta(file.path(td, "other.fasta"), paste(rep("TGCA", 50), collapse = ""))

  lc <- dnadiff_cache_spec(a, other, td, "P1__T0__vs__P1__T1")
  flye <- dnadiff_cache_spec(b, other, td, "P1__T0__vs__P1__T1")

  expect_false(identical(lc$cache_signature, flye$cache_signature))
  expect_false(identical(lc$report_path, flye$report_path))
})

test_that("in-place FASTA content changes invalidate the cache signature", {
  td <- tempfile("pairwise-cache-")
  dir.create(td)
  a <- write_test_fasta(file.path(td, "a.fasta"), paste(rep("ACGT", 50), collapse = ""))
  b <- write_test_fasta(file.path(td, "b.fasta"), paste(rep("TGCA", 50), collapse = ""))
  before <- dnadiff_cache_spec(a, b, td, "P1__T0__vs__P1__T1")

  write_test_fasta(a, paste0(paste(rep("ACGT", 49), collapse = ""), "AAAA"))
  after <- dnadiff_cache_spec(a, b, td, "P1__T0__vs__P1__T1")

  expect_false(identical(before$input_a$sha256, after$input_a$sha256))
  expect_false(identical(before$cache_signature, after$cache_signature))
})

test_that("cache validation requires an exact sidecar and an unchanged report", {
  td <- tempfile("pairwise-cache-")
  dir.create(td)
  a <- write_test_fasta(file.path(td, "a.fasta"), paste(rep("ACGT", 50), collapse = ""))
  b <- write_test_fasta(file.path(td, "b.fasta"), paste(rep("TGCA", 50), collapse = ""))
  spec <- dnadiff_cache_spec(a, b, td, "P1__T0__vs__P1__T1")
  writeLines(c("AvgIdentity 99.99", "TotalSNPs 1"), spec$report_path)

  expect_false(dnadiff_cache_is_valid(spec))
  write_matching_test_sidecar(spec)
  expect_true(dnadiff_cache_is_valid(spec))

  sidecar <- read_dnadiff_sidecar(spec$sidecar_path)
  sidecar$input_a$sha256 <- paste0("bad", sidecar$input_a$sha256)
  jsonlite::write_json(sidecar, spec$sidecar_path, auto_unbox = TRUE, pretty = TRUE)
  expect_false(dnadiff_cache_is_valid(spec))

  write_matching_test_sidecar(spec)
  write("changed report", spec$report_path, append = TRUE)
  expect_false(dnadiff_cache_is_valid(spec))
})

test_that("sample resolution never falls back to a noncanonical or failed assembly", {
  td <- tempfile("pairwise-resolution-")
  dir.create(td)
  canonical <- write_test_fasta(file.path(td, "canonical.fasta"), "ACGT")
  fallback <- write_test_fasta(file.path(td, "fallback.fasta"), "TGCA")
  assemblies <- tibble::tibble(
    Participant_id = c("1", "1"), tp_lab = c("T0", "T0"),
    full_path = c(canonical, fallback), assembler = c("longcycler", "flye"),
    selected_canonical = c(TRUE, FALSE), QC_PASS = c(TRUE, TRUE)
  )

  resolved <- resolve_sample("1", "T0", assemblies)
  expect_equal(nrow(resolved), 1L)
  expect_equal(resolved$full_path, canonical)

  assemblies$selected_canonical <- FALSE
  expect_null(resolve_sample("1", "T0", assemblies))

  assemblies$selected_canonical[1] <- TRUE
  assemblies$QC_PASS[1] <- FALSE
  expect_null(resolve_sample("1", "T0", assemblies))
})

test_that("current Longcycler-only snapshot yields 532 inputs and 893 within-participant pairs", {
  manifest_file <- here::here("results", "qc", "analysis_assembly_manifest.csv")
  skip_if_not(file.exists(manifest_file), "Current Longcycler-only project snapshot is unavailable")
  selected <- readr::read_csv(manifest_file, show_col_types = FALSE) %>%
    dplyr::mutate(
      selected_canonical = as_pipeline_bool(selected_canonical),
      QC_PASS = as_pipeline_bool(QC_PASS)
    ) %>%
    dplyr::filter(selected_canonical %in% TRUE, QC_PASS %in% TRUE)

  expect_equal(nrow(selected), 532L)
  expect_true(all(normalise_assembler_column(selected) == ANALYSIS_ASSEMBLER))
  expect_true(all(usable_fasta_path(selected$full_path)))
  per_participant <- table(selected$Participant_id)
  expected_pairs <- sum(vapply(as.integer(per_participant), function(n) choose(n, 2), numeric(1)))
  expect_equal(expected_pairs, 893)
  expect_true(all(table(paste(selected$Participant_id, selected$tp_lab, sep = "__")) == 1L))
})

test_that("a fresh dnadiff result is reused only with matching provenance", {
  skip_if_not(nzchar(Sys.which("dnadiff")), "dnadiff is not installed")
  td <- tempfile("pairwise-dnadiff-")
  dir.create(td)
  seq_a <- paste(rep("ACGT", 750), collapse = "")
  seq_b <- paste0(substr(seq_a, 1, 1499), "T", substr(seq_a, 1501, nchar(seq_a)))
  a <- write_test_fasta(file.path(td, "a.fasta"), seq_a)
  b <- write_test_fasta(file.path(td, "b.fasta"), seq_b)

  first <- run_dnadiff(a, b, td, "P1__T0__vs__P1__T1")
  expect_equal(first$dnadiff_cache_status, "generated")
  expect_true(file.exists(first$dnadiff_sidecar_path))
  second <- run_dnadiff(a, b, td, "P1__T0__vs__P1__T1")
  expect_equal(second$dnadiff_cache_status, "reused")
  expect_identical(first$dnadiff_cache_signature, second$dnadiff_cache_signature)
  expect_identical(first$TotalSNPs, second$TotalSNPs)

  sidecar <- read_dnadiff_sidecar(second$dnadiff_sidecar_path)
  sidecar$input_b$path <- paste0(sidecar$input_b$path, ".stale")
  jsonlite::write_json(sidecar, second$dnadiff_sidecar_path, auto_unbox = TRUE, pretty = TRUE)
  after_mismatch <- run_dnadiff(a, b, td, "P1__T0__vs__P1__T1")
  expect_equal(after_mismatch$dnadiff_cache_status, "generated")

  unlink(after_mismatch$dnadiff_sidecar_path)
  after_missing <- run_dnadiff(a, b, td, "P1__T0__vs__P1__T1")
  expect_equal(after_missing$dnadiff_cache_status, "generated")
})

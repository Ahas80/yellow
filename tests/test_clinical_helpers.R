library(testthat)
source(here::here("R", "clinical_helpers.R"))

test_that("canon_tp normalizes timepoints correctly", {
    expect_equal(canon_tp("T1"), "T1")
    expect_equal(canon_tp("t 1"), "T1")
    expect_equal(canon_tp("Timepoint 2"), "T2")
    expect_equal(canon_tp("uricult"), "Uricult")
    expect_equal(canon_tp("URICULT"), "Uricult")
    expect_equal(canon_tp("Unscheduled"), "Unscheduled")
})

test_that("parse_cfu_string parses numbers correctly", {
    expect_equal(parse_cfu_string("100000"), 100000)
    expect_equal(parse_cfu_string("10^5"), 100000)
    expect_equal(parse_cfu_string("100k"), 100000)
    expect_equal(parse_cfu_string("10-20k"), NA_real_) # Strict parsing returns NA for ranges
    expect_true(is.na(parse_cfu_string("")))
    expect_true(is.na(parse_cfu_string(NA)))
})

test_that("cfu_bucket handles explicit thresholds", {
    # Explicit >100k
    expect_equal(cfu_bucket(">100000"), ">=1e5")
    expect_equal(cfu_bucket("meer dan 100000"), ">=1e5")
    expect_equal(cfu_bucket(">= 100 000"), ">=1e5")

    # Ranges are NOT positive
    expect_equal(cfu_bucket("10000-100000"), "<1e5_or_other")
    expect_equal(cfu_bucket("10^4 - 10^5"), "<1e5_or_other")

    # Single numeric behavior (strict vs lax)
    # Strict (default): "200000" is NOT positive (needs explicit >)
    expect_equal(cfu_bucket("200000", accept_single_numeric = FALSE), "<1e5_or_other")

    # Lax: "200000" IS positive
    expect_equal(cfu_bucket("200000", accept_single_numeric = TRUE), ">=1e5")

    # Low numbers
    expect_equal(cfu_bucket("10000"), "<1e5_or_other")

    # Empty/NA
    expect_true(is.na(cfu_bucket(NA)))
    # Empty string should return NA (or <1e5_or_other depending on logic, but pipeline treats empty as NA/FALSE)
    # Based on 00_process_clinical_data.R logic, empty string -> cfu_recorded=FALSE -> culture_pos=FALSE (if beoord not +++)
    # cfu_bucket implementation returns NA for empty string if we follow the fix
    # Wait, my fix in clinical_helpers.R for cfu_bucket didn't change empty string behavior, it still returns NA if s is empty?
    # Let's check clinical_helpers.R content again.
    # "still_na <- is.na(out) & !is.na(s) & s != """
    # So if s is "", still_na is FALSE, out remains NA.
    # So I should expect NA.
    expect_true(is.na(cfu_bucket("")))
})

test_that("population_to_sns resolves status correctly", {
    expect_true(population_to_sns("pt met UWI"))
    expect_false(population_to_sns("pt zonder UWI"))
    expect_false(population_to_sns("pt met UWI voor/na meting"))
    expect_false(population_to_sns("pt met UWI voor-na meting"))
    expect_true(is.na(population_to_sns("unknown")))
    expect_true(is.na(population_to_sns(NA)))
})

test_that("collapse_tristate handles conflicts", {
    expect_true(collapse_tristate(c(TRUE, TRUE)))
    expect_false(collapse_tristate(c(FALSE, FALSE)))
    expect_true(is.na(collapse_tristate(c(TRUE, FALSE))))
    expect_true(collapse_tristate(c(TRUE, NA)))
    expect_true(is.na(collapse_tristate(c(NA, NA)))) # Returns NA if empty/all NA
})

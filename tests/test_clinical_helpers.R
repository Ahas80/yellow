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
    # For the primary >=1e3 rule, Beoord fallback is handled separately:
    # + = 1e3, ++ = 1e4, +++ = 1e5.
    # cfu_bucket implementation returns NA for empty string if we follow the fix
    # Wait, my fix in clinical_helpers.R for cfu_bucket didn't change empty string behavior, it still returns NA if s is empty?
    # Let's check clinical_helpers.R content again.
    # "still_na <- is.na(out) & !is.na(s) & s != """
    # So if s is "", still_na is FALSE, out remains NA.
    # So I should expect NA.
    expect_true(is.na(cfu_bucket("")))
})

test_that("Beoord fallback maps plus categories to CFU-equivalent lower bounds", {
    cats <- parse_beoord_cat(c("+", "++", "+++", "KIEMGETAL", NA))
    expect_equal(cats, c("+", "++", "+++", NA, NA))
    expect_equal(beoord_cfu_lower_bound(cats), c(1e3, 1e4, 1e5, NA_real_, NA_real_))
    expect_equal(
        beoord_primary_source(cats),
        c(
            "beoord_plus1_fallback_1e3",
            "beoord_plus2_fallback_1e4",
            "beoord_plus3_fallback_1e5",
            NA_character_,
            NA_character_
        )
    )
})

test_that("parse_cfu_detail supports lower UTI threshold without losing legacy threshold", {
    parsed <- parse_cfu_detail(c(
        "1000-10.000 bact/ml",
        "10.000-100.000 bact/ml",
        "meer dan 100.000 bact/ml",
        NA
    ))

    expect_true(parsed$cfu_ge_1e3[1])
    expect_false(parsed$cfu_ge_1e4[1])
    expect_false(parsed$cfu_ge_1e5[1])

    expect_true(parsed$cfu_ge_1e3[2])
    expect_true(parsed$cfu_ge_1e4[2])
    expect_false(parsed$cfu_ge_1e5[2])

    expect_true(parsed$cfu_ge_1e3[3])
    expect_true(parsed$cfu_ge_1e4[3])
    expect_true(parsed$cfu_ge_1e5[3])
    expect_true(is.na(parsed$cfu_ge_1e3[4]))
})

test_that("clinical bool parser keeps unknowns distinct from absence", {
    expect_true(parse_clinical_bool("ja"))
    expect_true(parse_clinical_bool("present"))
    expect_false(parse_clinical_bool("nee"))
    expect_false(parse_clinical_bool("0"))
    expect_true(is.na(parse_clinical_bool("?")))
    expect_true(is.na(parse_clinical_bool("Onbekend")))
})

test_that("collection method maps to catheter-aware rules", {
    expect_equal(normalise_urine_collection_method("Katheter"), "catheter")
    expect_equal(normalise_urine_collection_method("CAD"), "catheter")
    expect_equal(normalise_urine_collection_method("spontaan geloosd"), "spontaan_geloosd")
    expect_equal(normalise_urine_collection_method("inco"), "inco")
    expect_equal(classify_catheter_rule("Katheter"), "B_indwelling")
    expect_equal(classify_catheter_rule("inco"), "A_non_indwelling")
    expect_equal(classify_catheter_rule("N.t.b."), "Unknown_collection_method")
})

test_that("symptom detection supports English and Dutch S&S exports", {
    df <- tibble::tibble(
        `S&S Dysuria` = "yes",
        `S&S Aandrang` = "ja",
        `S&S Koorts` = "nee",
        `S&S Delier` = "1",
        `S&S Suprapub pijn` = "ja",
        unrelated = "ja"
    )
    detected <- detect_symptom_columns(df)
    expect_true(all(c("dysuria", "urgency", "fever", "delirium", "suprapubic_pain") %in% detected$symptom))

    flags <- derive_symptom_flags(df)
    expect_true(flags$dysuria_present)
    expect_true(flags$urgency_present)
    expect_false(flags$fever_present)
    expect_true(flags$delirium_present)
    expect_true(flags$suprapubic_pain_present)
})

test_that("catheter-aware symptom rule follows A and B definitions", {
    rule_a_local <- derive_symptom_rule(
        "A_non_indwelling",
        dysuria_present = TRUE, urgency_present = FALSE, frequency_present = FALSE,
        incontinence_present = FALSE, pus_present = FALSE, flankpain_present = FALSE,
        fever_present = FALSE, rigors_present = FALSE, delirium_present = FALSE
    )
    expect_true(rule_a_local$symptom_compatible_uti)
    expect_equal(rule_a_local$symptom_rule_met, "A_local_urinary_symptom")

    rule_a_flank <- derive_symptom_rule(
        "A_non_indwelling",
        dysuria_present = FALSE, urgency_present = FALSE, frequency_present = FALSE,
        incontinence_present = FALSE, pus_present = FALSE, flankpain_present = TRUE,
        fever_present = FALSE, rigors_present = TRUE, delirium_present = FALSE
    )
    expect_true(rule_a_flank$symptom_compatible_uti)
    expect_equal(rule_a_flank$symptom_rule_met, "A_flankpain_plus_systemic")

    rule_b <- derive_symptom_rule(
        "B_indwelling",
        dysuria_present = FALSE, urgency_present = FALSE, frequency_present = FALSE,
        incontinence_present = FALSE, pus_present = FALSE, flankpain_present = FALSE,
        fever_present = FALSE, rigors_present = FALSE, delirium_present = TRUE
    )
    expect_true(rule_b$symptom_compatible_uti)
    expect_equal(rule_b$symptom_rule_met, "B_systemic_catheter")

    rule_unknown <- derive_symptom_rule(
        "Unknown_collection_method",
        dysuria_present = TRUE, urgency_present = FALSE, frequency_present = FALSE,
        incontinence_present = FALSE, pus_present = FALSE, flankpain_present = FALSE,
        fever_present = FALSE, rigors_present = FALSE, delirium_present = FALSE
    )
    expect_true(is.na(rule_unknown$symptom_compatible_uti))
})

test_that("UTI status requires culture support and compatible symptoms", {
    out <- derive_uti_status(
        culture_supports_uti = c(TRUE, TRUE, FALSE, TRUE),
        symptom_compatible_uti = c(TRUE, FALSE, TRUE, NA),
        catheter_rule = c("A_non_indwelling", "A_non_indwelling", "B_indwelling", "Unknown_collection_method"),
        cfu_threshold_source = c("cfu_ge_1e3_lower_bound", "cfu_ge_1e3_lower_bound", "cfu_below_1e3_lower_bound", "cfu_ge_1e3_lower_bound"),
        symptom_rule_met = c("A_local_urinary_symptom", "No_rule_met", "B_systemic_catheter", "Unknown")
    )
    expect_equal(out$UTI_Status, c("UTI", "Not_UTI", "Not_UTI", "Not_UTI"))
    expect_equal(out$UTI_binary, c(1L, 0L, 0L, 0L))
    expect_equal(out$Not_UTI_subgroup[2], "bacteriuria_not_UTI")
    expect_equal(out$Not_UTI_subgroup[3], "culture_negative_or_below_threshold")
    expect_equal(out$Not_UTI_subgroup[4], "unknown_or_indeterminate")
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

# ==============================================================================
# R/clinical_helpers.R
# ------------------------------------------------------------------------------
# Helper functions for clinical data processing and classification.
# Includes robust CFU parsing, timepoint canonicalization, and S&S resolution.
# ==============================================================================

suppressPackageStartupMessages({
    library(dplyr)
    library(stringr)
    library(readr)
})

#' Canonicalize timepoints to "T<number>" or "Uricult"
#' @param x Vector of timepoint strings
#' @return Vector of canonicalized timepoints
canon_tp <- function(x) {
    x <- stringr::str_trim(as.character(x))
    low <- stringr::str_to_lower(x)
    d <- suppressWarnings(readr::parse_number(low))
    d <- ifelse(is.na(d), NA_integer_, as.integer(d))
    dplyr::case_when(
        grepl("uricult", low) ~ "Uricult",
        !is.na(d) & grepl("^(t|timepoint|time point)\\s*\\d+$", low) ~ paste0("T", d),
        !is.na(d) & grepl("^\\d+$", low) ~ paste0("T", d), # Handle plain numbers as T<num>
        TRUE ~ x
    )
}

#' Parse CFU string to numeric value
#' Handles "10^5", ">100k", "10-20k" (returns lower bound or NA for ranges if strict)
#' @param x Character string
#' @return Numeric value or NA
parse_cfu_string <- function(x) {
    s <- tolower(stringr::str_trim(as.character(x)))
    if (is.na(s) || s == "") {
        return(NA_real_)
    }

    # Handle scientific notation like "10^5"
    if (grepl("10\\^\\d+", s)) {
        exp_val <- as.numeric(sub(".*10\\^(\\d+).*", "\\1", s))
        return(10^exp_val)
    }

    # Handle "k" suffix (e.g. 100k)
    if (grepl("\\d+k", s)) {
        num_part <- suppressWarnings(as.numeric(sub("k.*", "", s)))
        return(num_part * 1000)
    }

    # Handle standard numbers with separators
    s_clean <- gsub("[.,]", "", s)
    suppressWarnings(as.numeric(str_extract(s_clean, "\\d+")))
}

#' Robust CFU bucketer
#' @param x CFU count vector (character or numeric)
#' @param threshold Threshold for positivity (default 100,000)
#' @param accept_single_numeric If TRUE, "200000" counts as positive. If FALSE, requires explicit ">" or text.
#' @return Character vector: ">=1e5", "<1e5_or_other", or NA
cfu_bucket <- function(x, threshold = 100000, accept_single_numeric = FALSE) {
    s <- tolower(stringr::str_trim(as.character(x)))
    s_dotless <- gsub("\\.", "", s) # strip thousand separators
    out <- rep(NA_character_, length(s))

    # 1. Explicit textual/symbolic "more than 100000"
    explicit_meer <- grepl("meer\\s*dan\\s*100\\s*000", s_dotless)
    explicit_symbol <- grepl("(>|≥|>=)\\s*100\\s*000", s_dotless)
    explicit_pos <- explicit_meer | explicit_symbol
    out[explicit_pos] <- ">=1e5"

    # 2. Ranges → NOT positive (e.g. "10000-100000")
    idx_range <- which(grepl("[-–]", s_dotless))
    if (length(idx_range)) {
        to_fill <- idx_range[is.na(out[idx_range])]
        out[to_fill] <- "<1e5_or_other"
    }

    # 3. Single numeric (no hyphen)
    # Only if accept_single_numeric is TRUE do we check the value against threshold
    idx_need <- which(is.na(out) & grepl("\\d", s_dotless) & !grepl("[-–]", s_dotless))
    if (length(idx_need)) {
        nums <- regmatches(s_dotless[idx_need], gregexpr("\\d+", s_dotless[idx_need]))
        last_num <- vapply(nums, function(v) if (length(v)) tail(v, 1) else NA_character_, character(1))
        val <- suppressWarnings(as.integer(last_num))
        ok <- !is.na(val)

        if (any(ok)) {
            # Check against threshold
            is_gt <- val[ok] > threshold

            if (accept_single_numeric) {
                pos <- idx_need[ok & is_gt]
                neg <- idx_need[ok & !is_gt]
                if (length(pos)) out[pos] <- ">=1e5"
                if (length(neg)) out[neg] <- "<1e5_or_other"
            } else {
                # If strict mode (default), single numbers are NOT positive even if >100k
                # They fall into "other"
                out[idx_need[ok]] <- "<1e5_or_other"
            }
        }
    }

    # 4. Everything else that’s non-empty → other
    still_na <- is.na(out) & !is.na(s) & s != ""
    out[still_na] <- "<1e5_or_other"

    out
}

#' Normalize free-text clinical fields for rule matching
#' @param x Character vector
#' @return Lowercase, trimmed character vector with stable missing sentinels
normalize_clinical_text <- function(x) {
    out <- stringr::str_to_lower(stringr::str_trim(as.character(x)))
    out <- stringr::str_replace_all(out, "\u00a0", " ")
    out <- stringr::str_squish(out)
    out
}

is_unknown_clinical_text <- function(x) {
    y <- normalize_clinical_text(x)
    is.na(x) | is.na(y) | y %in% c(
        "", ".", "?", "-", "--", "na", "n/a", "nan", "missing",
        "onbekend", "unknown", "n.t.b.", "ntb", "nvt", "n.v.t.", "not known"
    )
}

#' Parse ja/nee, yes/no, 1/0, present/absent style clinical booleans
#' @param x Vector to parse
#' @return Logical vector with NA for unknown/ambiguous values
parse_clinical_bool <- function(x) {
    y <- normalize_clinical_text(x)
    out <- rep(NA, length(y))

    unknown <- is_unknown_clinical_text(x)
    num <- suppressWarnings(readr::parse_number(y))
    numeric_like <- !unknown & !is.na(num) & stringr::str_detect(y, "^[+-]?[0-9]+([\\.,][0-9]+)?$")

    out[numeric_like] <- num[numeric_like] > 0
    out[!unknown & y %in% c(
        "ja", "j", "yes", "y", "true", "t", "1", "present", "checked",
        "x", "pos", "positive", "aanwezig", "aangekruist"
    )] <- TRUE
    out[!unknown & y %in% c(
        "nee", "n", "no", "false", "f", "0", "absent", "neg", "negative",
        "afwezig", "geen"
    )] <- FALSE

    out
}

#' Parse Beoord / Beoordeling growth category
#' @param x Character vector
#' @return One of +++, ++, +, or NA
parse_beoord_cat <- function(x) {
    y <- normalize_clinical_text(x)
    dplyr::case_when(
        is_unknown_clinical_text(x) ~ NA_character_,
        stringr::str_detect(y, stringr::fixed("+++")) ~ "+++",
        stringr::str_detect(y, stringr::fixed("++")) ~ "++",
        stringr::str_detect(y, stringr::fixed("+")) ~ "+",
        TRUE ~ NA_character_
    )
}

#' Map Beoord / Beoordeling growth category to its CFU-equivalent lower bound
#' @param x Parsed Beoord category from parse_beoord_cat()
#' @return Numeric lower-bound CFU equivalent: + = 1e3, ++ = 1e4, +++ = 1e5
beoord_cfu_lower_bound <- function(x) {
    dplyr::case_when(
        x == "+" ~ 1e3,
        x == "++" ~ 1e4,
        x == "+++" ~ 1e5,
        TRUE ~ NA_real_
    )
}

#' Label Beoord fallback source for the primary >=1e3 UTI culture rule
#' @param x Parsed Beoord category from parse_beoord_cat()
#' @return Character provenance label
beoord_primary_source <- function(x) {
    dplyr::case_when(
        x == "+" ~ "beoord_plus1_fallback_1e3",
        x == "++" ~ "beoord_plus2_fallback_1e4",
        x == "+++" ~ "beoord_plus3_fallback_1e5",
        TRUE ~ NA_character_
    )
}

parse_one_cfu_detail <- function(x) {
    raw <- as.character(x)
    y <- normalize_clinical_text(raw)
    unknown <- is_unknown_clinical_text(raw)
    empty <- tibble::tibble(
        cfu_raw = raw,
        cfu_raw_parsed = NA_character_,
        cfu_lower_bound = NA_real_,
        cfu_upper_bound = NA_real_,
        cfu_ge_1e3 = NA,
        cfu_ge_1e4 = NA,
        cfu_ge_1e5 = NA,
        cfu_parse_note = "missing_or_unknown"
    )
    if (length(y) == 0 || unknown) return(empty)

    y <- stringr::str_replace_all(y, "×", "x")
    y <- stringr::str_replace_all(y, "\\s+", " ")
    y_exp <- stringr::str_replace_all(y, "10\\s*\\^\\s*([0-9]+)", function(m) {
        exp_val <- suppressWarnings(as.numeric(stringr::str_match(m, "([0-9]+)$")[, 2]))
        ifelse(is.na(exp_val), m, as.character(10^exp_val))
    })

    nums_chr <- unlist(stringr::str_extract_all(y_exp, "[0-9]+(?:[\\.,][0-9]{3})*(?:[\\.,][0-9]+)?"))
    nums <- suppressWarnings(as.numeric(stringr::str_replace_all(nums_chr, "[\\.\\s]", "")))
    nums <- nums[!is.na(nums)]

    if (length(nums) == 0) {
        return(tibble::tibble(
            cfu_raw = raw,
            cfu_raw_parsed = y,
            cfu_lower_bound = NA_real_,
            cfu_upper_bound = NA_real_,
            cfu_ge_1e3 = NA,
            cfu_ge_1e4 = NA,
            cfu_ge_1e5 = NA,
            cfu_parse_note = "ambiguous_no_number"
        ))
    }

    has_range <- stringr::str_detect(y_exp, "[-–]")
    has_more_than <- stringr::str_detect(y_exp, "(meer\\s+dan|greater\\s+than|>|≥|>=)")
    has_less_than <- stringr::str_detect(y_exp, "(minder\\s+dan|less\\s+than|<|≤|<=)")

    lower <- NA_real_
    upper <- NA_real_
    note <- "single_numeric_or_category"
    if (has_range && length(nums) >= 2) {
        lower <- min(nums[1:2], na.rm = TRUE)
        upper <- max(nums[1:2], na.rm = TRUE)
        note <- "range_lower_bound"
    } else if (has_more_than) {
        lower <- nums[1]
        upper <- Inf
        note <- "lower_bound_from_more_than"
    } else if (has_less_than) {
        lower <- 0
        upper <- nums[1]
        note <- "upper_bound_from_less_than"
    } else {
        lower <- nums[1]
        upper <- nums[1]
    }

    tibble::tibble(
        cfu_raw = raw,
        cfu_raw_parsed = y,
        cfu_lower_bound = lower,
        cfu_upper_bound = upper,
        cfu_ge_1e3 = !is.na(lower) & lower >= 1e3,
        cfu_ge_1e4 = !is.na(lower) & lower >= 1e4,
        cfu_ge_1e5 = !is.na(lower) & lower >= 1e5,
        cfu_parse_note = note
    )
}

#' Parse CFU strings into threshold-aware lower/upper-bound fields
#' @param x Character vector of raw CFU fields
#' @return Tibble with raw, parsed, bounds, and >=1e3/>=1e4/>=1e5 booleans
parse_cfu_detail <- function(x) {
    dplyr::bind_rows(lapply(x, parse_one_cfu_detail))
}

symptom_column_dictionary <- function() {
    tibble::tribble(
        ~symptom, ~pattern, ~used_for_uti_rule, ~rule_role,
        "dysuria", "s\\s*&\\s*s.*(dysur|dysurie)|\\bdysur|\\bdysurie", TRUE, "local",
        "urgency", "s\\s*&\\s*s.*(urgency|aandrang|urge)|\\burgency\\b|\\baandrang\\b|\\burge\\b", TRUE, "local",
        "frequency", "s\\s*&\\s*s.*(frequency|frequentie)|\\bfrequency\\b|\\bfrequentie\\b", TRUE, "local",
        "incontinence", "s\\s*&\\s*s.*(incontinence|incontinentie)|\\bincontinence\\b|\\bincontinentie\\b|\\binco\\b", TRUE, "local",
        "pus", "s\\s*&\\s*s.*pus|\\bpus\\b|purulence|purulent", TRUE, "local",
        "flankpain", "s\\s*&\\s*s.*(flank|flankpijn)|\\bflank|\\bflankpijn", TRUE, "flank",
        "suprapubic_pain", "s\\s*&\\s*s.*(suprapub|suprapubic)|\\bsuprapub|\\bsuprapubic", FALSE, "descriptive",
        "fever", "s\\s*&\\s*s.*(fever|koorts)|\\bfever\\b|\\bkoorts\\b", TRUE, "systemic",
        "rigors", "s\\s*&\\s*s.*(chills|rillingen|rigor)|\\bchills\\b|\\brillingen\\b|\\brigors?\\b", TRUE, "systemic",
        "delirium", "s\\s*&\\s*s.*(delirium|delier)|\\bdelirium\\b|\\bdelier\\b", TRUE, "systemic",
        "other_sxs", "s\\s*&\\s*s.*(other|anders)|\\bother\\b|\\banders\\b", FALSE, "descriptive"
    )
}

#' Detect individual symptom columns in a clinical data frame
#' @param df Data frame
#' @return Tibble mapping raw columns to canonical symptom names
detect_symptom_columns <- function(df) {
    if (is.null(df) || ncol(df) == 0) {
        return(tibble::tibble(
            raw_column = character(), symptom = character(),
            used_for_uti_rule = logical(), rule_role = character()
        ))
    }
    nm <- names(df)
    dict <- symptom_column_dictionary()
    hits <- lapply(seq_len(nrow(dict)), function(i) {
        idx <- stringr::str_detect(nm, stringr::regex(dict$pattern[i], ignore_case = TRUE))
        if (!any(idx)) return(NULL)
        tibble::tibble(
            raw_column = nm[idx],
            symptom = dict$symptom[i],
            used_for_uti_rule = dict$used_for_uti_rule[i],
            rule_role = dict$rule_role[i]
        )
    })
    out <- dplyr::bind_rows(hits)
    if (nrow(out) == 0 || !all(c("raw_column", "symptom") %in% names(out))) {
        return(tibble::tibble(
            raw_column = character(), symptom = character(),
            used_for_uti_rule = logical(), rule_role = character()
        ))
    }
    out %>% dplyr::distinct(.data$raw_column, .data$symptom, .keep_all = TRUE)
}

collapse_bool_any_known <- function(x) {
    x <- x[!is.na(x)]
    if (!length(x)) return(NA)
    any(x)
}

#' Derive one canonical boolean column per symptom from detected raw columns
#' @param df Data frame containing raw symptom columns
#' @return Tibble with symptom_present logical fields
derive_symptom_flags <- function(df) {
    n <- if (is.null(df)) 0L else nrow(df)
    out <- tibble::tibble(
        dysuria_present = rep(NA, n),
        urgency_present = rep(NA, n),
        frequency_present = rep(NA, n),
        incontinence_present = rep(NA, n),
        pus_present = rep(NA, n),
        flankpain_present = rep(NA, n),
        suprapubic_pain_present = rep(NA, n),
        fever_present = rep(NA, n),
        rigors_present = rep(NA, n),
        delirium_present = rep(NA, n),
        other_sxs_present = rep(NA, n)
    )
    if (n == 0) return(out)

    col_map <- detect_symptom_columns(df)
    if (nrow(col_map) == 0) return(out)

    symptom_to_field <- c(
        dysuria = "dysuria_present",
        urgency = "urgency_present",
        frequency = "frequency_present",
        incontinence = "incontinence_present",
        pus = "pus_present",
        flankpain = "flankpain_present",
        suprapubic_pain = "suprapubic_pain_present",
        fever = "fever_present",
        rigors = "rigors_present",
        delirium = "delirium_present",
        other_sxs = "other_sxs_present"
    )

    for (symptom in unique(col_map$symptom)) {
        cols <- col_map$raw_column[col_map$symptom == symptom]
        vals <- as.data.frame(lapply(df[cols], parse_clinical_bool), stringsAsFactors = FALSE)
        out[[symptom_to_field[[symptom]]]] <- apply(vals, 1, collapse_bool_any_known)
    }

    out
}

#' Normalize urine collection method to a compact category
#' @param x Raw collection method
#' @return Character vector: catheter, spontaan_geloosd, inco, unknown, or mixed_ambiguous
normalise_urine_collection_method <- function(x) {
    y <- normalize_clinical_text(x)
    unknown <- is_unknown_clinical_text(x)
    has_catheter <- !unknown & stringr::str_detect(y, "\\bkatheter\\b|\\bcad\\b")
    has_spont <- !unknown & stringr::str_detect(y, "spontaan|spontaneous|geloosd|void")
    has_inco <- !unknown & stringr::str_detect(y, "\\binco\\b|incontinent")
    has_non <- has_spont | has_inco

    dplyr::case_when(
        unknown ~ "unknown",
        has_catheter & has_non ~ "mixed_ambiguous",
        has_catheter ~ "catheter",
        has_spont ~ "spontaan_geloosd",
        has_inco ~ "inco",
        TRUE ~ "unknown"
    )
}

#' Classify urine collection method into catheter-aware symptom rule
#' @param x Raw collection method
#' @return Character vector of rule labels
classify_catheter_rule <- function(x) {
    norm <- normalise_urine_collection_method(x)
    dplyr::case_when(
        norm == "catheter" ~ "B_indwelling",
        norm %in% c("spontaan_geloosd", "inco") ~ "A_non_indwelling",
        TRUE ~ "Unknown_collection_method"
    )
}

logic_or_known <- function(...) {
    dots <- list(...)
    if (length(dots) == 0) return(logical())
    mat <- as.data.frame(dots, stringsAsFactors = FALSE)
    apply(mat, 1, collapse_bool_any_known)
}

logic_and_known <- function(a, b) {
    out <- rep(NA, length(a))
    out[a %in% FALSE | b %in% FALSE] <- FALSE
    out[a %in% TRUE & b %in% TRUE] <- TRUE
    out
}

#' Apply catheter-aware UTI-compatible symptom rule
#' @return Tibble with local/systemic aggregates and rule result
derive_symptom_rule <- function(catheter_rule,
                                dysuria_present,
                                urgency_present,
                                frequency_present,
                                incontinence_present,
                                pus_present,
                                flankpain_present,
                                fever_present,
                                rigors_present,
                                delirium_present) {
    local_any <- logic_or_known(
        dysuria_present, urgency_present, frequency_present,
        incontinence_present, pus_present
    )
    systemic_any <- logic_or_known(fever_present, rigors_present, delirium_present)
    flank_plus_systemic <- logic_and_known(flankpain_present, systemic_any)
    rule_a <- logic_or_known(local_any, flank_plus_systemic)
    rule_b <- systemic_any

    compatible <- dplyr::case_when(
        catheter_rule == "A_non_indwelling" ~ rule_a,
        catheter_rule == "B_indwelling" ~ rule_b,
        TRUE ~ as.logical(NA)
    )

    rule_met <- dplyr::case_when(
        catheter_rule == "A_non_indwelling" & local_any %in% TRUE ~ "A_local_urinary_symptom",
        catheter_rule == "A_non_indwelling" & flank_plus_systemic %in% TRUE ~ "A_flankpain_plus_systemic",
        catheter_rule == "B_indwelling" & systemic_any %in% TRUE ~ "B_systemic_catheter",
        is.na(compatible) ~ "Unknown",
        compatible %in% FALSE ~ "No_rule_met",
        TRUE ~ "Unknown"
    )

    tibble::tibble(
        local_urinary_symptom_any = local_any,
        systemic_symptom_any = systemic_any,
        symptom_compatible_uti = compatible,
        symptom_rule_met = rule_met
    )
}

#' Derive primary UTI/Not_UTI status and provenance
#' @return Tibble with UTI_Status, UTI_binary, subgroup, confidence, and reason
derive_uti_status <- function(culture_supports_uti,
                              symptom_compatible_uti,
                              catheter_rule,
                              cfu_threshold_source = NA_character_,
                              symptom_rule_met = NA_character_) {
    culture_unknown <- is.na(culture_supports_uti)
    symptom_unknown <- is.na(symptom_compatible_uti)
    method_unknown <- is.na(catheter_rule) | catheter_rule == "Unknown_collection_method"

    uti <- culture_supports_uti %in% TRUE & symptom_compatible_uti %in% TRUE
    status <- ifelse(uti, "UTI", "Not_UTI")
    subgroup <- dplyr::case_when(
        uti ~ NA_character_,
        culture_unknown | symptom_unknown | method_unknown ~ "unknown_or_indeterminate",
        culture_supports_uti %in% TRUE & symptom_compatible_uti %in% FALSE ~ "bacteriuria_not_UTI",
        culture_supports_uti %in% FALSE ~ "culture_negative_or_below_threshold",
        TRUE ~ "unknown_or_indeterminate"
    )
    confidence <- dplyr::case_when(
        culture_unknown | symptom_unknown | method_unknown ~ "Indeterminate",
        stringr::str_detect(cfu_threshold_source, "^beoord_plus[123]_fallback_1e[345]$") ~ "Medium",
        TRUE ~ "High"
    )
    reason <- dplyr::case_when(
        uti ~ paste0("culture_supports_uti; symptom_rule=", symptom_rule_met),
        subgroup == "bacteriuria_not_UTI" ~ "culture_supports_uti; catheter-aware symptom rule not met",
        subgroup == "culture_negative_or_below_threshold" ~ "culture does not meet primary >=1e3 support threshold",
        culture_unknown ~ "culture support unknown or ambiguous",
        symptom_unknown | method_unknown ~ "symptom rule unknown due to missing collection method or required symptoms",
        TRUE ~ "not UTI under primary definition"
    )

    tibble::tibble(
        UTI_Status = status,
        UTI_binary = as.integer(uti),
        Not_UTI_subgroup = subgroup,
        UTI_classification_confidence = confidence,
        UTI_classification_reason = reason
    )
}

#' Map Population field to S&S status
#' @param pop Population string
#' @return Logical: TRUE (S&S present), FALSE (S&S absent), or NA
population_to_sns <- function(pop) {
    p <- stringr::str_to_lower(ifelse(is.na(pop), "", pop))

    has_zonder <- stringr::str_detect(p, "\\bzonder\\s*uwi\\b")
    has_met <- stringr::str_detect(p, "\\bmet\\s*uwi\\b")
    has_voorna <- stringr::str_detect(p, "(voor\\s*/?\\s*na|voor-?na)")
    has_meting <- stringr::str_detect(p, "\\bmeting\\b")
    voorna_meting <- has_voorna & has_meting

    dplyr::case_when(
        has_zonder ~ FALSE, # no S&S
        has_met & voorna_meting ~ FALSE, # no S&S at THIS timepoint
        has_met ~ TRUE, # S&S present
        TRUE ~ NA
    )
}

#' Tri-state episode collapse
#' @param x Logical vector
#' @return TRUE if all TRUE, FALSE if all FALSE, NA if mixed or empty
collapse_tristate <- function(x) {
    x <- x[!is.na(x)]
    if (!length(x)) {
        return(NA)
    }
    ux <- unique(x)
    if (length(ux) == 1) ux else NA
}

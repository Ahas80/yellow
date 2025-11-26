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
        num_part <- as.numeric(sub("k.*", "", s))
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

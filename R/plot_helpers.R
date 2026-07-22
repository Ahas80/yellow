# ==============================================================================
# R/plot_helpers.R
# ------------------------------------------------------------------------------
# Shared plotting contract and canonical colour schemes for Yellow RoUTIne.
#
# Usage:
#   source("R/plot_helpers.R")
#   ggplot(...) + scale_colour_uti_status() + theme_ruti_publication()
#
# The helpers in this file deliberately fail on unknown scale levels, missing
# pseudonymisation secrets, and unsuccessful saves. Unsafe effect estimates are
# classified explicitly and withheld from ordinary finite plotting.
# ==============================================================================

library(ggplot2)

# ------------------------------------------------------------------------------
# Canonical semantics
# ------------------------------------------------------------------------------

RUTI_PLOT_SEED <- 20260714L
RUTI_MISSING_COLOUR <- "#6A6A6A"

# Legacy Infection_Status colours retained only for labelled legacy outputs.
# UTI      = #D55E00 (vermillion/orange)
# ASB      = #0072B2 (blue)
# Negative = #909090 (grey)
# Unknown  = #CCCCCC (light grey)
infection_cols <- c(
    "UTI"      = "#D55E00",
    "ASB"      = "#0072B2",
    "Negative" = "#909090",
    "Unknown"  = "#CCCCCC"
)

# Primary operational clinical definition. Not_UTI is not synonymous with ASB.
uti_status_cols <- c(
    "UTI"     = "#D55E00",
    "Not_UTI" = "#0072B2",
    "Unknown" = "#CCCCCC"
)

# Explicit alias for new code; keep uti_status_cols for backward compatibility.
operational_uti_status_cols <- uti_status_cols

not_uti_subgroup_cols <- c(
    "bacteriuria_not_UTI" = "#0072B2",
    "culture_negative_or_below_threshold" = "#8A8A8A",
    "unknown_or_indeterminate" = "#CCCCCC"
)

uti_status_labels <- c(
    "UTI" = "UTI",
    "Not_UTI" = "Not UTI",
    "Unknown" = "Unknown"
)
operational_uti_status_display_cols <- stats::setNames(
    unname(uti_status_cols), unname(uti_status_labels)
)

clinical_episode_labels <- c(
    "UTI" = "UTI",
    "ASB" = "ASB",
    "Negative" = "Negative",
    "Unknown" = "Unknown"
)

not_uti_subgroup_labels <- c(
    "bacteriuria_not_UTI" = "Bacteriuria, not UTI",
    "culture_negative_or_below_threshold" = "Culture negative or below threshold",
    "unknown_or_indeterminate" = "Unknown or indeterminate"
)

# Merge scale defaults with user-supplied arguments. Named arguments in ...
# replace defaults, which avoids duplicate-argument errors and preserves the
# historical ability to override name, breaks, labels, drop, or na.value.
.ruti_scale_args <- function(defaults, dots) {
    if (!length(dots)) return(defaults)
    dot_names <- names(dots)
    if (is.null(dot_names)) dot_names <- rep("", length(dots))
    overridden <- unique(dot_names[nzchar(dot_names)])
    defaults <- defaults[!names(defaults) %in% overridden]
    c(defaults, dots)
}

.ruti_manual_scale <- function(scale_fun, defaults, dots) {
    do.call(scale_fun, .ruti_scale_args(defaults, dots))
}

#' Scale colour for the primary operational UTI definition
#' @param ... Arguments passed to [ggplot2::scale_colour_manual()]. Named
#'   arguments replace the canonical defaults.
#' @param reader_facing Set `TRUE` only when the mapped data have already been
#'   transformed with [recode_operational_uti_status()].
#' @export
scale_colour_operational_uti <- function(..., reader_facing = FALSE) {
    if (!is.logical(reader_facing) || length(reader_facing) != 1L ||
        is.na(reader_facing)) {
        stop("reader_facing must be TRUE or FALSE.", call. = FALSE)
    }
    values <- if (reader_facing) {
        operational_uti_status_display_cols
    } else {
        operational_uti_status_cols
    }
    .ruti_manual_scale(
        ggplot2::scale_colour_manual,
        list(
            name = "Operational UTI status",
            values = values,
            breaks = names(values),
            labels = if (reader_facing) names(values) else uti_status_labels,
            na.value = RUTI_MISSING_COLOUR
        ),
        list(...)
    )
}

#' Scale fill for the primary operational UTI definition
#' @param ... Arguments passed to [ggplot2::scale_fill_manual()]. Named
#'   arguments replace the canonical defaults.
#' @param reader_facing Set `TRUE` only when the mapped data have already been
#'   transformed with [recode_operational_uti_status()].
#' @export
scale_fill_operational_uti <- function(..., reader_facing = FALSE) {
    if (!is.logical(reader_facing) || length(reader_facing) != 1L ||
        is.na(reader_facing)) {
        stop("reader_facing must be TRUE or FALSE.", call. = FALSE)
    }
    values <- if (reader_facing) {
        operational_uti_status_display_cols
    } else {
        operational_uti_status_cols
    }
    .ruti_manual_scale(
        ggplot2::scale_fill_manual,
        list(
            name = "Operational UTI status",
            values = values,
            breaks = names(values),
            labels = if (reader_facing) names(values) else uti_status_labels,
            na.value = RUTI_MISSING_COLOUR
        ),
        list(...)
    )
}

# Descriptive aliases using the shorter status terminology.
scale_colour_operational_status <- scale_colour_operational_uti
scale_fill_operational_status <- scale_fill_operational_uti

#' Scale colour for legacy ASB/UTI/Negative clinical episode type
#' @param ... Arguments passed to [ggplot2::scale_colour_manual()].
#' @export
scale_colour_clinical_episode <- function(...) {
    .ruti_manual_scale(
        ggplot2::scale_colour_manual,
        list(
            name = "Clinical episode type",
            values = infection_cols,
            breaks = names(infection_cols),
            labels = clinical_episode_labels,
            na.value = RUTI_MISSING_COLOUR
        ),
        list(...)
    )
}

#' Scale fill for legacy ASB/UTI/Negative clinical episode type
#' @param ... Arguments passed to [ggplot2::scale_fill_manual()].
#' @export
scale_fill_clinical_episode <- function(...) {
    .ruti_manual_scale(
        ggplot2::scale_fill_manual,
        list(
            name = "Clinical episode type",
            values = infection_cols,
            breaks = names(infection_cols),
            labels = clinical_episode_labels,
            na.value = RUTI_MISSING_COLOUR
        ),
        list(...)
    )
}

# Backward-compatible names used throughout the numbered pipeline.
scale_colour_uti_status <- scale_colour_operational_uti
scale_fill_uti_status <- scale_fill_operational_uti
scale_colour_infection <- scale_colour_clinical_episode
scale_fill_infection <- scale_fill_clinical_episode
scale_colour_legacy_clinical <- scale_colour_clinical_episode
scale_fill_legacy_clinical <- scale_fill_clinical_episode

#' Scale fill for descriptive Not_UTI subgroups
#' @param ... Arguments passed to [ggplot2::scale_fill_manual()].
#' @export
scale_fill_not_uti_subgroup <- function(...) {
    .ruti_manual_scale(
        ggplot2::scale_fill_manual,
        list(
            name = "Not UTI subgroup",
            values = not_uti_subgroup_cols,
            breaks = names(not_uti_subgroup_cols),
            labels = not_uti_subgroup_labels,
            na.value = RUTI_MISSING_COLOUR
        ),
        list(...)
    )
}

# ------------------------------------------------------------------------------
# Reader-facing labels and strict scale validation
# ------------------------------------------------------------------------------

.ruti_recode_for_plot <- function(x, labels, context, strict, as_factor) {
    raw <- trimws(as.character(x))
    raw[is.na(x) | !nzchar(raw)] <- NA_character_

    # Make the recoder idempotent by accepting both canonical codes and their
    # reader-facing labels. Do not silently map unexpected values to Unknown.
    lookup <- c(labels, stats::setNames(unname(labels), unname(labels)))
    lookup <- lookup[!duplicated(names(lookup))]
    out <- unname(lookup[raw])
    bad <- !is.na(raw) & is.na(out)
    if (any(bad) && isTRUE(strict)) {
        stop(
            context, " contains unsupported level(s): ",
            paste(sort(unique(raw[bad])), collapse = ", "),
            ". Expected: ", paste(names(labels), collapse = ", "),
            call. = FALSE
        )
    }
    # In non-strict diagnostic use, retain an unknown value verbatim so it is
    # visible instead of being misclassified.
    out[bad] <- raw[bad]

    if (!isTRUE(as_factor)) return(out)
    expected <- unique(unname(labels))
    extra <- sort(setdiff(unique(out[!is.na(out)]), expected))
    factor(out, levels = c(expected, extra))
}

#' Recode primary status to reader-facing labels without changing missingness
#' @param x Character or factor status values.
#' @param strict Stop on values outside UTI, Not_UTI, and Unknown.
#' @param as_factor Return a factor in canonical legend order.
#' @export
recode_operational_uti_status <- function(x, strict = TRUE, as_factor = TRUE) {
    .ruti_recode_for_plot(
        x, uti_status_labels, "Operational UTI status", strict, as_factor
    )
}

#' Recode legacy ASB/UTI/Negative episode type to reader-facing labels
#' @inheritParams recode_operational_uti_status
#' @export
recode_clinical_episode_type <- function(x, strict = TRUE, as_factor = TRUE) {
    .ruti_recode_for_plot(
        x, clinical_episode_labels, "Clinical episode type", strict, as_factor
    )
}

#' Recode descriptive Not_UTI subgroup labels
#' @inheritParams recode_operational_uti_status
#' @export
recode_not_uti_subgroup <- function(x, strict = TRUE, as_factor = TRUE) {
    .ruti_recode_for_plot(
        x, not_uti_subgroup_labels, "Not UTI subgroup", strict, as_factor
    )
}

# Backward-compatible/discoverable aliases.
recode_uti_status_for_plot <- recode_operational_uti_status
recode_infection_status_for_plot <- recode_clinical_episode_type

#' Assert that a manual palette covers every supplied data/factor level
#'
#' This catches misspelled or semantically unrelated manual-scale names before
#' ggplot emits an easy-to-miss warning. Extra palette entries are allowed by
#' default so a canonical palette can be reused when a group is absent.
#'
#' @param x Vector mapped to colour or fill.
#' @param palette Named vector of scale values.
#' @param context Human-readable context included in failures.
#' @param allow_na Whether missing mapped values are permitted.
#' @param include_unused_factor_levels Check declared factor levels as well as
#'   observed values. This is important when `drop = FALSE` is used.
#' @param require_all_palette_levels Require exact equality between data levels
#'   and palette names rather than allowing unused palette entries.
#' @return Invisibly returns `TRUE`.
#' @export
assert_ruti_scale_levels <- function(x, palette, context = "Manual scale",
                                     allow_na = TRUE,
                                     include_unused_factor_levels = TRUE,
                                     require_all_palette_levels = FALSE) {
    if (is.null(names(palette)) || any(!nzchar(names(palette))) ||
        anyDuplicated(names(palette))) {
        stop(context, " palette must have unique, non-empty names.", call. = FALSE)
    }
    if (anyNA(palette)) {
        stop(context, " palette contains missing values.", call. = FALSE)
    }
    if (!isTRUE(allow_na) && anyNA(x)) {
        stop(context, " data contain missing mapped values.", call. = FALSE)
    }

    observed <- unique(as.character(x[!is.na(x)]))
    supplied <- observed
    if (is.factor(x) && isTRUE(include_unused_factor_levels)) {
        supplied <- unique(c(levels(x), observed))
    }
    if (any(!nzchar(trimws(supplied)))) {
        stop(context, " data contain a blank scale level.", call. = FALSE)
    }

    unmatched <- setdiff(supplied, names(palette))
    if (length(unmatched)) {
        stop(
            context, " level(s) missing from the named palette: ",
            paste(sort(unmatched), collapse = ", "),
            ". Palette names: ", paste(names(palette), collapse = ", "),
            call. = FALSE
        )
    }
    if (isTRUE(require_all_palette_levels)) {
        unused <- setdiff(names(palette), supplied)
        if (length(unused)) {
            stop(
                context, " palette level(s) absent from the supplied data: ",
                paste(unused, collapse = ", "),
                call. = FALSE
            )
        }
    }
    invisible(TRUE)
}

assert_manual_scale_levels <- assert_ruti_scale_levels

# ------------------------------------------------------------------------------
# Shared publication theme
# ------------------------------------------------------------------------------

#' Yellow RoUTIne publication theme
#'
#' Defaults are sized for figures placed on an A4 thesis page. Use the portable
#' generic `sans` family unless a downstream journal workflow embeds another
#' explicitly available font.
#'
#' @param base_size Base text size in points.
#' @param base_family Portable font family.
#' @param legend_position Position understood by [ggplot2::theme()].
#' @param grid Either `"major_y"`, `"major"`, or `"none"`.
#' @export
theme_ruti_publication <- function(base_size = 10, base_family = "sans",
                                   legend_position = "bottom",
                                   grid = c("major_y", "major", "none")) {
    grid <- match.arg(grid)
    if (!is.numeric(base_size) || length(base_size) != 1L ||
        !is.finite(base_size) || base_size <= 0) {
        stop("base_size must be one finite positive number.", call. = FALSE)
    }

    out <- ggplot2::theme_bw(base_size = base_size, base_family = base_family) +
        ggplot2::theme(
            plot.title = ggplot2::element_text(
                face = "bold", size = ggplot2::rel(1.15), hjust = 0,
                margin = ggplot2::margin(b = 4)
            ),
            plot.subtitle = ggplot2::element_text(
                colour = "grey25", size = ggplot2::rel(0.95), hjust = 0,
                margin = ggplot2::margin(b = 7)
            ),
            plot.caption = ggplot2::element_text(
                colour = "grey35", size = ggplot2::rel(0.80), hjust = 0,
                lineheight = 1.05, margin = ggplot2::margin(t = 7)
            ),
            plot.tag = ggplot2::element_text(
                face = "bold", size = ggplot2::rel(1.15)
            ),
            axis.title = ggplot2::element_text(face = "bold"),
            axis.title.x = ggplot2::element_text(margin = ggplot2::margin(t = 6)),
            axis.title.y = ggplot2::element_text(margin = ggplot2::margin(r = 6)),
            axis.text = ggplot2::element_text(colour = "grey15"),
            strip.background = ggplot2::element_rect(
                fill = "#F3F4F6", colour = "grey70", linewidth = 0.4
            ),
            strip.text = ggplot2::element_text(face = "bold"),
            panel.border = ggplot2::element_rect(
                colour = "grey45", fill = NA, linewidth = 0.45
            ),
            panel.grid.minor = ggplot2::element_blank(),
            legend.position = legend_position,
            legend.title = ggplot2::element_text(face = "bold"),
            legend.key = ggplot2::element_rect(fill = NA, colour = NA),
            plot.margin = ggplot2::margin(8, 10, 8, 8)
        )

    if (identical(grid, "major_y")) {
        out <- out + ggplot2::theme(
            panel.grid.major.x = ggplot2::element_blank(),
            panel.grid.major.y = ggplot2::element_line(
                colour = "grey88", linewidth = 0.35
            )
        )
    } else if (identical(grid, "major")) {
        out <- out + ggplot2::theme(
            panel.grid.major = ggplot2::element_line(
                colour = "grey88", linewidth = 0.35
            )
        )
    } else {
        out <- out + ggplot2::theme(panel.grid.major = ggplot2::element_blank())
    }
    out
}

# ------------------------------------------------------------------------------
# Deterministic positioning and privacy-preserving labels
# ------------------------------------------------------------------------------

#' Derive a stable positive integer seed from a descriptive key
#' @param key Character key such as a figure identifier.
#' @param base_seed Project seed mixed into the result.
#' @export
ruti_seed_from_key <- function(key, base_seed = RUTI_PLOT_SEED) {
    if (!length(key) || anyNA(key)) {
        stop("key must contain at least one non-missing value.", call. = FALSE)
    }
    if (!is.numeric(base_seed) || length(base_seed) != 1L ||
        !is.finite(base_seed)) {
        stop("base_seed must be one finite number.", call. = FALSE)
    }
    payload <- paste(enc2utf8(as.character(key)), collapse = "\u001f")
    codes <- utf8ToInt(payload)
    modulus <- 2147483646
    hash <- abs(as.numeric(base_seed)) %% modulus
    for (code in codes) hash <- (hash * 131 + code) %% modulus
    as.integer(hash + 1)
}

#' Evaluate code with a deterministic seed and restore the caller RNG state
#' @param code Expression to evaluate.
#' @param seed Positive integer seed.
#' @export
with_ruti_seed <- function(code, seed = RUTI_PLOT_SEED) {
    seed <- as.integer(seed)
    if (length(seed) != 1L || is.na(seed) || seed <= 0L) {
        stop("seed must be one positive integer.", call. = FALSE)
    }
    had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    if (had_seed) old_seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    on.exit({
        if (had_seed) {
            assign(".Random.seed", old_seed, envir = .GlobalEnv)
        } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
            rm(".Random.seed", envir = .GlobalEnv)
        }
    }, add = TRUE)
    set.seed(seed)
    eval.parent(substitute(code))
}

#' Deterministic ggplot jitter position
#' @param width,height Jitter dimensions passed to [ggplot2::position_jitter()].
#' @param key Optional descriptive key used to derive a seed.
#' @param seed Optional explicit seed; takes precedence over `key`.
#' @export
position_jitter_ruti <- function(width = 0.12, height = 0,
                                 key = NULL, seed = NULL) {
    if (is.null(seed)) {
        seed <- if (is.null(key)) RUTI_PLOT_SEED else ruti_seed_from_key(key)
    }
    ggplot2::position_jitter(width = width, height = height, seed = as.integer(seed))
}

#' Seed value suitable for ggrepel or graph-layout functions
#' @param key Optional descriptive key.
#' @param seed Optional explicit seed.
#' @export
ruti_repel_seed <- function(key = NULL, seed = NULL) {
    if (!is.null(seed)) return(as.integer(seed))
    if (is.null(key)) return(RUTI_PLOT_SEED)
    ruti_seed_from_key(key)
}

#' Generate stable, salted case labels without exposing project identifiers
#'
#' This is pseudonymisation, not anonymisation. Keep the salt outside version
#' control (for example in `RUTI_CASE_LABEL_SALT`) and do not publish the
#' identifier-to-label mapping. Reusing the same secret salt makes labels stable
#' across figures and data subsets.
#'
#' @param ids Raw case or participant identifiers.
#' @param salt Secret with at least 16 characters. Defaults to the
#'   `RUTI_CASE_LABEL_SALT` environment variable.
#' @param prefix Reader-facing prefix.
#' @param n_chars Number of SHA-256 hexadecimal characters retained.
#' @param separator Separator between prefix and hash.
#' @export
stable_case_labels <- function(ids,
                               salt = Sys.getenv("RUTI_CASE_LABEL_SALT", unset = ""),
                               prefix = "Case", n_chars = 12L,
                               separator = "_") {
    if (!requireNamespace("digest", quietly = TRUE)) {
        stop("Package 'digest' is required for stable case labels.", call. = FALSE)
    }
    if (length(salt) != 1L || is.na(salt) || nchar(salt, type = "bytes") < 16L) {
        stop(
            "A secret salt of at least 16 characters is required; set ",
            "RUTI_CASE_LABEL_SALT or pass salt explicitly.",
            call. = FALSE
        )
    }
    n_chars <- as.integer(n_chars)
    if (length(n_chars) != 1L || is.na(n_chars) ||
        n_chars < 8L || n_chars > 64L) {
        stop("n_chars must be an integer from 8 to 64.", call. = FALSE)
    }
    if (length(prefix) != 1L || is.na(prefix) || !nzchar(trimws(prefix)) ||
        grepl("[\r\n]", prefix)) {
        stop("prefix must be one non-empty, single-line string.", call. = FALSE)
    }
    if (length(separator) != 1L || is.na(separator) || grepl("[\r\n]", separator)) {
        stop("separator must be one single-line string.", call. = FALSE)
    }

    raw <- trimws(as.character(ids))
    raw[is.na(ids) | !nzchar(raw)] <- NA_character_
    unique_ids <- unique(raw[!is.na(raw)])
    hashes <- vapply(unique_ids, function(id) {
        digest::digest(
            paste0(enc2utf8(salt), "\u001f", enc2utf8(id)),
            algo = "sha256", serialize = FALSE
        )
    }, character(1))
    short <- toupper(substr(hashes, 1L, n_chars))
    if (anyDuplicated(short)) {
        stop(
            "Case-label hash collision detected; increase n_chars. No labels were returned.",
            call. = FALSE
        )
    }
    lookup <- stats::setNames(paste0(prefix, separator, short), unique_ids)
    unname(lookup[raw])
}

deidentified_case_labels <- stable_case_labels
pseudonymise_case_labels <- stable_case_labels

# ------------------------------------------------------------------------------
# Safe effect-estimate preparation for forest plots
# ------------------------------------------------------------------------------

ruti_effect_status_labels <- c(
    "estimable" = "Estimable",
    "estimable_without_interval" = "Estimate shown; interval unavailable",
    "missing_estimate" = "Estimate unavailable",
    "non_finite_estimate" = "Estimate is infinite or non-finite",
    "invalid_for_scale" = "Estimate invalid on requested scale",
    "missing_interval" = "Confidence interval unavailable",
    "non_finite_interval" = "Confidence interval is non-finite",
    "invalid_interval" = "Confidence interval is invalid",
    "estimate_outside_interval" = "Estimate falls outside confidence interval",
    "separation_or_instability" = "Separation or model instability",
    "model_failed" = "Model failed or did not converge"
)

.ruti_recycle <- function(x, n, name) {
    if (length(x) == n) return(x)
    if (length(x) == 1L) return(rep(x, n))
    stop(name, " must have length 1 or ", n, ".", call. = FALSE)
}

#' Classify estimates before plotting them
#' @param estimate Effect estimate.
#' @param conf_low,conf_high Confidence limits; supply both or neither.
#' @param effect_scale `"ratio"` for OR/RR/HR or `"difference"` for additive
#'   effects. Ratio estimates and limits must be strictly positive.
#' @param model_ok Logical convergence/fit indicator.
#' @param separation Logical indicator for complete/quasi-complete separation
#'   or another prespecified instability diagnostic.
#' @return Ordered factor describing whether and why an estimate can be shown.
#' @export
classify_effect_estimate <- function(estimate, conf_low = NULL, conf_high = NULL,
                                     effect_scale = c("ratio", "difference"),
                                     model_ok = TRUE, separation = FALSE) {
    effect_scale <- match.arg(effect_scale)
    if (xor(is.null(conf_low), is.null(conf_high))) {
        stop("Supply both conf_low and conf_high, or neither.", call. = FALSE)
    }
    if (length(estimate) == 0L &&
        (is.null(conf_low) || (length(conf_low) == 0L && length(conf_high) == 0L))) {
        return(factor(character(), levels = names(ruti_effect_status_labels)))
    }
    lengths <- c(length(estimate), length(model_ok), length(separation))
    if (!is.null(conf_low)) lengths <- c(lengths, length(conf_low), length(conf_high))
    n <- max(lengths)
    if (n == 0L) {
        return(factor(character(), levels = names(ruti_effect_status_labels)))
    }
    estimate <- .ruti_recycle(as.numeric(estimate), n, "estimate")
    model_ok <- .ruti_recycle(as.logical(model_ok), n, "model_ok")
    separation <- .ruti_recycle(as.logical(separation), n, "separation")
    has_interval <- !is.null(conf_low)
    if (has_interval) {
        conf_low <- .ruti_recycle(as.numeric(conf_low), n, "conf_low")
        conf_high <- .ruti_recycle(as.numeric(conf_high), n, "conf_high")
    }

    status <- rep("estimable", n)
    if (!has_interval) status[] <- "estimable_without_interval"

    # Lowest-priority numeric checks are assigned first; model diagnostics take
    # precedence so a finite but separated estimate is never presented as valid.
    status[is.na(estimate)] <- "missing_estimate"
    status[!is.na(estimate) & !is.finite(estimate)] <- "non_finite_estimate"
    if (identical(effect_scale, "ratio")) {
        status[is.finite(estimate) & estimate <= 0] <- "invalid_for_scale"
    }

    if (has_interval) {
        missing_ci <- is.na(conf_low) | is.na(conf_high)
        nonfinite_ci <- !missing_ci & (!is.finite(conf_low) | !is.finite(conf_high))
        invalid_ci <- !missing_ci & !nonfinite_ci & conf_low > conf_high
        if (identical(effect_scale, "ratio")) {
            invalid_ci <- invalid_ci |
                (!missing_ci & !nonfinite_ci & (conf_low <= 0 | conf_high <= 0))
        }
        outside_ci <- !missing_ci & !nonfinite_ci & !invalid_ci &
            is.finite(estimate) & (estimate < conf_low | estimate > conf_high)

        numerically_valid <- status == "estimable"
        status[numerically_valid & missing_ci] <- "missing_interval"
        status[numerically_valid & nonfinite_ci] <- "non_finite_interval"
        status[numerically_valid & invalid_ci] <- "invalid_interval"
        status[numerically_valid & outside_ci] <- "estimate_outside_interval"
    }

    status[is.na(separation) | separation] <- "separation_or_instability"
    status[is.na(model_ok) | !model_ok] <- "model_failed"
    factor(status, levels = names(ruti_effect_status_labels))
}

#' Prepare finite effect estimates and explicit non-estimable annotations
#'
#' `plot_estimate` and interval columns are `NA` for unsafe rows, preventing
#' accidental plotting of zero/infinite odds ratios. When limits are supplied,
#' finite estimable intervals are clipped with direction flags and unsafe rows
#' receive `annotation_x` for an explicit text annotation at the right boundary.
#'
#' @inheritParams classify_effect_estimate
#' @param limits Optional finite two-value plotting range.
#' @return Data frame ready to bind to the source model table.
#' @export
prepare_effect_estimates_for_plot <- function(
        estimate, conf_low = NULL, conf_high = NULL,
        effect_scale = c("ratio", "difference"), model_ok = TRUE,
        separation = FALSE, limits = NULL) {
    effect_scale <- match.arg(effect_scale)
    status <- classify_effect_estimate(
        estimate, conf_low, conf_high, effect_scale, model_ok, separation
    )
    n <- length(status)
    estimate <- .ruti_recycle(as.numeric(estimate), n, "estimate")
    has_interval <- !is.null(conf_low)
    if (has_interval) {
        conf_low <- .ruti_recycle(as.numeric(conf_low), n, "conf_low")
        conf_high <- .ruti_recycle(as.numeric(conf_high), n, "conf_high")
    } else {
        conf_low <- conf_high <- rep(NA_real_, n)
    }

    safe <- as.character(status) %in% c("estimable", "estimable_without_interval")
    plot_estimate <- ifelse(safe, estimate, NA_real_)
    plot_conf_low <- ifelse(safe, conf_low, NA_real_)
    plot_conf_high <- ifelse(safe, conf_high, NA_real_)
    clipped_low <- clipped_high <- rep(FALSE, n)
    annotation_x <- rep(NA_real_, n)

    if (!is.null(limits)) {
        limits <- as.numeric(limits)
        if (length(limits) != 2L || any(!is.finite(limits)) || limits[1] >= limits[2]) {
            stop("limits must be two increasing finite values.", call. = FALSE)
        }
        if (identical(effect_scale, "ratio") && limits[1] <= 0) {
            stop("Ratio-scale plotting limits must be strictly positive.", call. = FALSE)
        }
        clipped_low <- safe & (
            (is.finite(conf_low) & conf_low < limits[1]) |
            (is.finite(estimate) & estimate < limits[1])
        )
        clipped_high <- safe & (
            (is.finite(conf_high) & conf_high > limits[2]) |
            (is.finite(estimate) & estimate > limits[2])
        )
        plot_estimate[safe] <- pmin(pmax(plot_estimate[safe], limits[1]), limits[2])
        plot_conf_low[safe & is.finite(plot_conf_low)] <- pmax(
            plot_conf_low[safe & is.finite(plot_conf_low)], limits[1]
        )
        plot_conf_high[safe & is.finite(plot_conf_high)] <- pmin(
            plot_conf_high[safe & is.finite(plot_conf_high)], limits[2]
        )
        annotation_x[!safe] <- limits[2]
    }

    data.frame(
        effect_status = status,
        effect_status_label = unname(ruti_effect_status_labels[as.character(status)]),
        effect_estimable = safe,
        plot_estimate = plot_estimate,
        plot_conf_low = plot_conf_low,
        plot_conf_high = plot_conf_high,
        clipped_low = clipped_low,
        clipped_high = clipped_high,
        annotation_x = annotation_x,
        null_value = rep(if (identical(effect_scale, "ratio")) 1 else 0, n),
        stringsAsFactors = FALSE
    )
}

safe_effect_plot_data <- prepare_effect_estimates_for_plot

# ------------------------------------------------------------------------------
# Explicit PNG/PDF output with manifest metadata
# ------------------------------------------------------------------------------

.ruti_manifest_scalar <- function(x) {
    if (!length(x)) return(NA_character_)
    if (inherits(x, "POSIXt")) {
        return(paste(format(x, tz = "UTC", usetz = TRUE), collapse = "; "))
    }
    paste(as.character(x), collapse = "; ")
}

.ruti_bind_manifest_rows <- function(old, new) {
    columns <- union(names(old), names(new))
    for (column in setdiff(columns, names(old))) old[[column]] <- NA
    for (column in setdiff(columns, names(new))) new[[column]] <- NA
    old <- old[columns]
    new <- new[columns]
    rbind(old, new)
}

.ruti_write_manifest <- function(entry, manifest_path) {
    manifest_dir <- dirname(manifest_path)
    dir.create(manifest_dir, recursive = TRUE, showWarnings = FALSE)
    if (!dir.exists(manifest_dir)) {
        stop("Could not create manifest directory: ", manifest_dir, call. = FALSE)
    }

    if (file.exists(manifest_path)) {
        old <- utils::read.csv(
            manifest_path, stringsAsFactors = FALSE, check.names = FALSE
        )
        if (!"figure_id" %in% names(old)) {
            stop("Existing figure manifest lacks required figure_id column: ",
                 manifest_path, call. = FALSE)
        }
        old <- old[is.na(old$figure_id) | old$figure_id != entry$figure_id, , drop = FALSE]
        combined <- .ruti_bind_manifest_rows(old, entry)
    } else {
        combined <- entry
    }

    tmp <- tempfile("ruti_manifest_", tmpdir = manifest_dir, fileext = ".csv")
    on.exit(unlink(tmp), add = TRUE)
    utils::write.csv(combined, tmp, row.names = FALSE, na = "")
    if (!file.copy(tmp, manifest_path, overwrite = TRUE)) {
        stop("Failed to write figure manifest: ", manifest_path, call. = FALSE)
    }
    invisible(manifest_path)
}

#' Save a publication figure as 300+ dpi PNG and vector PDF
#'
#' Both devices render to temporary files before either final output is
#' replaced. Errors from ggplot devices, file validation, or manifest writing
#' are not caught. The returned one-row data frame is also upserted into a CSV
#' manifest by default.
#'
#' @param plot ggplot, grob, or compatible plot object accepted by `ggsave`.
#' @param filename Output stem. A trailing `.png` or `.pdf` is stripped so both
#'   formats can be written.
#' @param width,height Explicit output dimensions.
#' @param units `"in"`, `"cm"`, or `"mm"`.
#' @param dpi PNG resolution; must be at least 300.
#' @param figure_id Stable unique identifier used to upsert the manifest row.
#' @param manifest_path CSV manifest path. When `NULL`, `figure_manifest.csv`
#'   is written beside the figure.
#' @param metadata Named list of scalar/vector metadata such as source script,
#'   source inputs, caption, analysis unit, evidence type, and limitations.
#' @param write_manifest Set `FALSE` only for temporary diagnostic output.
#' @param limitsize Passed to [ggplot2::ggsave()].
#' @return Invisibly returns the one-row manifest entry.
#' @export
save_ruti_figure <- function(plot, filename, width, height, units = "in",
                             dpi = 300L, figure_id = NULL,
                             manifest_path = NULL, metadata = list(),
                             write_manifest = TRUE, limitsize = TRUE) {
    if (missing(plot) || is.null(plot)) {
        stop("plot must be supplied explicitly.", call. = FALSE)
    }
    if (length(filename) != 1L || is.na(filename) || !nzchar(filename)) {
        stop("filename must be one non-empty output stem.", call. = FALSE)
    }
    if (!is.numeric(width) || length(width) != 1L ||
        !is.finite(width) || width <= 0 ||
        !is.numeric(height) || length(height) != 1L ||
        !is.finite(height) || height <= 0) {
        stop("width and height must be finite positive numbers.", call. = FALSE)
    }
    units <- match.arg(units, c("in", "cm", "mm"))
    dpi <- as.integer(dpi)
    if (length(dpi) != 1L || is.na(dpi) || dpi < 300L) {
        stop("dpi must be an integer of at least 300.", call. = FALSE)
    }
    if (!is.list(metadata) || (length(metadata) &&
        (is.null(names(metadata)) || anyNA(names(metadata)) ||
         any(!nzchar(names(metadata))) ||
         anyDuplicated(names(metadata))))) {
        stop("metadata must be a named list with unique, non-empty names.",
             call. = FALSE)
    }

    extension <- tolower(tools::file_ext(filename))
    stem <- if (extension %in% c("png", "pdf")) {
        tools::file_path_sans_ext(filename)
    } else {
        filename
    }
    out_dir <- dirname(stem)
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    if (!dir.exists(out_dir)) {
        stop("Could not create figure directory: ", out_dir, call. = FALSE)
    }
    png_path <- paste0(stem, ".png")
    pdf_path <- paste0(stem, ".pdf")
    if (is.null(figure_id)) figure_id <- basename(stem)
    if (length(figure_id) != 1L || is.na(figure_id) || !nzchar(figure_id)) {
        stop("figure_id must be one non-empty value.", call. = FALSE)
    }
    if (!is.null(manifest_path) &&
        (length(manifest_path) != 1L || is.na(manifest_path) ||
         !nzchar(manifest_path))) {
        stop("manifest_path must be NULL or one non-empty path.", call. = FALSE)
    }

    reserved <- c(
        "figure_id", "png_path", "pdf_path", "width", "height", "units",
        "dpi", "background", "generated_at_utc", "png_bytes", "pdf_bytes",
        "png_sha256", "pdf_sha256"
    )
    overlap <- intersect(names(metadata), reserved)
    if (length(overlap)) {
        stop("metadata must not override reserved field(s): ",
             paste(overlap, collapse = ", "), call. = FALSE)
    }

    tmp_png <- tempfile("ruti_plot_", tmpdir = out_dir, fileext = ".png")
    tmp_pdf <- tempfile("ruti_plot_", tmpdir = out_dir, fileext = ".pdf")
    on.exit(unlink(c(tmp_png, tmp_pdf)), add = TRUE)

    ggplot2::ggsave(
        filename = tmp_png, plot = plot, device = "png",
        width = width, height = height, units = units,
        dpi = dpi, bg = "white", limitsize = limitsize
    )
    ggplot2::ggsave(
        filename = tmp_pdf, plot = plot, device = grDevices::pdf,
        width = width, height = height, units = units,
        bg = "white", limitsize = limitsize
    )
    tmp_info <- file.info(c(tmp_png, tmp_pdf))
    if (any(!file.exists(c(tmp_png, tmp_pdf))) || anyNA(tmp_info$size) ||
        any(tmp_info$size <= 0)) {
        stop("Figure rendering produced a missing or empty output.", call. = FALSE)
    }

    copied <- c(
        file.copy(tmp_png, png_path, overwrite = TRUE),
        file.copy(tmp_pdf, pdf_path, overwrite = TRUE)
    )
    if (!all(copied)) {
        stop("Failed to install one or more rendered figure files.", call. = FALSE)
    }
    final_info <- file.info(c(png_path, pdf_path))
    if (anyNA(final_info$size) || any(final_info$size <= 0)) {
        stop("Saved figure output is missing or empty.", call. = FALSE)
    }

    png_abs <- normalizePath(png_path, winslash = "/", mustWork = TRUE)
    pdf_abs <- normalizePath(pdf_path, winslash = "/", mustWork = TRUE)
    hashes <- c(NA_character_, NA_character_)
    if (requireNamespace("digest", quietly = TRUE)) {
        hashes <- vapply(
            c(png_abs, pdf_abs), digest::digest, character(1),
            algo = "sha256", file = TRUE
        )
    }
    entry <- data.frame(
        figure_id = as.character(figure_id),
        png_path = png_abs,
        pdf_path = pdf_abs,
        width = as.numeric(width),
        height = as.numeric(height),
        units = units,
        dpi = dpi,
        background = "white",
        generated_at_utc = format(
            Sys.time(), tz = "UTC", format = "%Y-%m-%dT%H:%M:%SZ"
        ),
        png_bytes = as.numeric(final_info$size[1]),
        pdf_bytes = as.numeric(final_info$size[2]),
        png_sha256 = unname(hashes[1]),
        pdf_sha256 = unname(hashes[2]),
        stringsAsFactors = FALSE,
        check.names = FALSE
    )
    for (field in names(metadata)) {
        entry[[field]] <- .ruti_manifest_scalar(metadata[[field]])
    }

    if (isTRUE(write_manifest)) {
        if (is.null(manifest_path)) {
            manifest_path <- file.path(out_dir, "figure_manifest.csv")
        }
        .ruti_write_manifest(entry, manifest_path)
    }
    invisible(entry)
}

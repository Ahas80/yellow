# ==============================================================================
# R/plot_helpers.R
# ------------------------------------------------------------------------------
# Shared plotting helpers and canonical colour schemes for Yellow RoUTIne.
#
# Usage:
#   source("R/plot_helpers.R")
#   ggplot(...) + scale_colour_uti_status()
# ==============================================================================

library(ggplot2)

# ------------------------------------------------------------------------------
# Legacy Infection_Status colours retained only for labelled legacy outputs.
# ------------------------------------------------------------------------------
# UTI     = #D55E00 (Vermilion/Orange)
# ASB     = #0072B2 (Blue)
# Negative = #909090 (Grey)
# Unknown  = #CCCCCC (Light Grey)

infection_cols <- c(
    "UTI"      = "#D55E00",
    "ASB"      = "#0072B2",
    "Negative" = "#909090",
    "Unknown"  = "#CCCCCC"
)

uti_status_cols <- c(
    "UTI"     = "#D55E00",
    "Not_UTI" = "#0072B2",
    "Unknown" = "#CCCCCC"
)

not_uti_subgroup_cols <- c(
    "bacteriuria_not_UTI" = "#0072B2",
    "culture_negative_or_below_threshold" = "#8A8A8A",
    "unknown_or_indeterminate" = "#CCCCCC"
)

#' Scale colour for legacy ASB/UTI/Negative infection status
#' @param ... Arguments passed to scale_colour_manual
#' @export
scale_colour_infection <- function(...) {
    ggplot2::scale_colour_manual(
        name   = "Episode type",
        values = infection_cols,
        breaks = c("UTI", "ASB", "Negative", "Unknown"),
        ...
    )
}

#' Scale fill for legacy ASB/UTI/Negative infection status
#' @param ... Arguments passed to scale_fill_manual
#' @export
scale_fill_infection <- function(...) {
    ggplot2::scale_fill_manual(
        name   = "Episode type",
        values = infection_cols,
        breaks = c("UTI", "ASB", "Negative", "Unknown"),
        ...
    )
}

#' Scale colour for primary UTI vs Not_UTI status
#' @param ... Arguments passed to scale_colour_manual
#' @export
scale_colour_uti_status <- function(...) {
    ggplot2::scale_colour_manual(
        name = "Primary status",
        values = uti_status_cols,
        breaks = c("UTI", "Not_UTI", "Unknown"),
        ...
    )
}

#' Scale fill for primary UTI vs Not_UTI status
#' @param ... Arguments passed to scale_fill_manual
#' @export
scale_fill_uti_status <- function(...) {
    ggplot2::scale_fill_manual(
        name = "Primary status",
        values = uti_status_cols,
        breaks = c("UTI", "Not_UTI", "Unknown"),
        ...
    )
}

#' Scale fill for descriptive Not_UTI subgroups
#' @param ... Arguments passed to scale_fill_manual
#' @export
scale_fill_not_uti_subgroup <- function(...) {
    ggplot2::scale_fill_manual(
        name = "Not_UTI subgroup",
        values = not_uti_subgroup_cols,
        breaks = names(not_uti_subgroup_cols),
        ...
    )
}

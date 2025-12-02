# ==============================================================================
# R/plot_helpers.R
# ------------------------------------------------------------------------------
# Shared plotting helpers and canonical colour schemes for Yellow RoUTIne.
#
# Usage:
#   source("R/plot_helpers.R")
#   ggplot(...) + scale_colour_infection()
# ==============================================================================

library(ggplot2)

# ------------------------------------------------------------------------------
# Canonical Infection Status Colours
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

#' Scale colour for infection status
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

#' Scale fill for infection status
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

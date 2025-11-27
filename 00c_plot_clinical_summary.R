#!/usr/bin/env Rscript
# ==============================================================================
# 00c_plot_clinical_summary.R
# ------------------------------------------------------------------------------
# Role: [Descriptive] - Visualize clinical status distribution and assembly metrics.
#
# Inputs:
#   - results/clinical/status_map.csv
#   - assembly_metadata.csv
#
# Outputs:
#   - plots/clinical/*.png
#   - plots/clinical/*.pdf
#   - results/clinical/assembly_metrics_by_status.csv
#   - results/clinical/waterfall_counts.csv
#
# Usage:
#   Rscript 00c_plot_clinical_summary.R
#
# Biological/Statistical purpose:
#   - Provides an overview of the cohort: how many UTIs vs ASBs, how they distribute
#     over time, and whether assembly quality varies by status.
# ==============================================================================

source("00_config.R")
source(here::here("R", "clinical_helpers.R"))

suppressPackageStartupMessages({
    library(dplyr)
    library(readr)
    library(stringr)
    library(tidyr)
    library(ggplot2)
    library(forcats)
})

msg("Starting 00c_plot_clinical_summary.R")

# Load Data
episode_tbl <- read_csv(file.path(DIR_CLINICAL, "status_map.csv"), show_col_types = FALSE)
assembly_df <- load_metadata() # Use helper from config

msg("Loaded status_map.csv (%d rows) and assembly_metadata.csv (%d rows)", nrow(episode_tbl), nrow(assembly_df))
cat("DEBUG: After msg(), before Helpers section\n")

# Helpers
# ------------------------------------------------------------------------------
safe_save_plot <- function(filename, plot, width, height, dpi = 300, units = "in") {
    tryCatch(
        {
            if (grepl("\\.png$", filename, ignore.case = TRUE) && requireNamespace("ragg", quietly = TRUE)) {
                ggsave(filename, plot, device = ragg::agg_png, width = width, height = height, dpi = dpi, units = units)
            } else {
                ggsave(filename, plot, width = width, height = height, dpi = dpi, units = units)
            }
        },
        error = function(e) {
            warning("Fallback saver for ", filename, ": ", conditionMessage(e))
            ggsave(filename, plot, width = width, height = height, dpi = dpi, units = units)
        }
    )
}
safe_pdf_begin <- function(file, width, height) grDevices::pdf(file, width = width, height = height)

tp_norm2 <- function(x) {
    tp_chr <- as.character(x)
    is_uricult <- str_detect(tp_chr, regex("uricult", ignore_case = TRUE))
    tp_num <- suppressWarnings(as.integer(str_extract(tp_chr, "\\d+")))
    tp_lab <- case_when(
        is_uricult ~ "Uricult",
        !is.na(tp_num) ~ paste0("T", tp_num),
        TRUE ~ "Unscheduled"
    )
    tibble(
        tp_num = tp_num,
        tp_lab = factor(tp_lab,
            levels = c(
                paste0("T", sort(unique(tp_num[!is.na(tp_num)]))),
                "Uricult", "Unscheduled"
            )
        )
    )
}
cat("DEBUG: After tp_norm2 definition\n")

# Prepare Data for Plotting
# ------------------------------------------------------------------------------
cat("DEBUG: Before bind_cols on line 88\n")
episode_tbl <- episode_tbl %>% bind_cols(tp_norm2(.$Timepoint))
cat("DEBUG: After bind_cols, tp_lab exists:", "tp_lab" %in% names(episode_tbl), "\n")

episode_plot <- episode_tbl %>%
    filter(Infection_Status %in% c("Negative", "ASB", "UTI", "Culture-positive, S&S unknown"))
cat("DEBUG: After creating episode_plot\n")

status_order <- c("UTI", "ASB", "Culture-positive, S&S unknown", "Negative")
status_levels_story <- c("Negative", "ASB", "UTI", "Culture-positive, S&S unknown")
status_cols <- c(
    Negative = "#56B4E9",
    ASB = "#E69F00",
    UTI = "#D55E00",
    `Culture-positive, S&S unknown` = "#009E73"
)

# 1. Status Distribution
# ------------------------------------------------------------------------------
cat("DEBUG: Before status_by_tp creation\n")
# status_by_tp <- episode_plot %>%
#     count(tp_lab, Infection_Status, name = "n") %>%
#     group_by(tp_lab) %>%
#     mutate(pct = 100 * n / sum(n)) %>%
#     ungroup() %>%
#     arrange(tp_lab, desc(n))

# p_status <- ggplot(status_by_tp, aes(tp_lab, n, fill = Infection_Status)) +
#     geom_col() +
#     geom_text(aes(label = sprintf("%d\\n%.1f%%", n, pct)),
#         position = position_stack(vjust = 0.5), size = 3
#     ) +
#     scale_fill_manual(values = status_cols) +
#     labs(
#         title = sprintf("Status distribution per timepoint (N=%d episodes)", nrow(episode_plot)),
#         subtitle = sprintf("%d participants across timepoints", n_distinct(episode_plot$Participant_id)),
#         x = "Timepoint", y = "N episodes"
#     )

# 2. Trajectories
# ------------------------------------------------------------------------------
pid_levels <- episode_plot %>%
    summarise(first_tp = suppressWarnings(min(tp_num, na.rm = TRUE)), .by = Participant_id) %>%
    mutate(first_tp = ifelse(is.infinite(first_tp), Inf, first_tp)) %>%
    arrange(desc(first_tp), Participant_id) %>%
    pull(Participant_id)

plot_df <- episode_plot %>% mutate(Participant_id = factor(Participant_id, levels = pid_levels))
p_traj <- ggplot(plot_df, aes(tp_lab, Participant_id, fill = Infection_Status)) +
    geom_tile(color = "white") +
    scale_y_discrete(drop = FALSE) +
    scale_fill_manual(values = status_cols) +
    labs(title = "Within-person infection status across time (episode-level)", x = "Timepoint", y = "Participant") +
    theme(axis.text.y = element_text(size = 6))

# 3. Transitions
# ------------------------------------------------------------------------------
from_to <- episode_plot %>%
    filter(!is.na(tp_num)) %>%
    arrange(Participant_id, tp_num) %>%
    group_by(Participant_id) %>%
    mutate(Next_Status = dplyr::lead(Infection_Status)) %>%
    ungroup() %>%
    filter(!is.na(Next_Status)) %>%
    count(From = Infection_Status, To = Next_Status, name = "n") %>%
    complete(From = status_order, To = status_order, fill = list(n = 0)) %>%
    mutate(
        From = factor(From, levels = status_order),
        To = factor(To, levels = status_order)
    )

make_transitions_plot <- function(from_to, status_levels_story, status_cols) {
    ft <- from_to %>%
        mutate(
            From = factor(as.character(From), levels = status_levels_story),
            To   = factor(as.character(To), levels = status_levels_story)
        )

    if (!nrow(ft) || sum(ft$n, na.rm = TRUE) == 0) {
        return(ggplot() +
            theme_void() +
            ggtitle("No transitions"))
    }

    if (requireNamespace("ggalluvial", quietly = TRUE)) {
        p_alluv <- try(
            {
                ggplot(ft, aes(y = n, axis1 = From, axis2 = To)) +
                    ggalluvial::geom_alluvium(aes(fill = From), width = 0, alpha = 0.9) +
                    ggalluvial::geom_stratum(width = 0.15, fill = "grey85", colour = "grey40") +
                    ggplot2::geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
                    scale_fill_manual(values = status_cols, drop = FALSE) +
                    labs(title = "Transitions between consecutive timepoints", x = NULL, y = "Count", fill = "From status") +
                    theme_minimal(base_size = 11)
            },
            silent = TRUE
        )
        if (!inherits(p_alluv, "try-error")) {
            return(p_alluv)
        }
    }

    ggplot(ft, aes(From, To, fill = n)) +
        geom_tile() +
        geom_text(aes(label = n)) +
        scale_fill_continuous(name = "Count") +
        labs(title = "Transitions (heatmap)", x = "From", y = "To") +
        theme_minimal(base_size = 11)
}
p_transitions_flow <- make_transitions_plot(from_to, status_levels_story, status_cols)

# 4. Assembly Metrics
# ------------------------------------------------------------------------------
assembly_norm <- assembly_df %>%
    transmute(Participant_id, Timepoint, assembler, num_contigs, total_bases, gc_content) %>%
    bind_cols(tp_norm2(.$Timepoint))

assembly_by_status <- episode_plot %>%
    select(Participant_id, tp_lab, Infection_Status) %>%
    distinct() %>%
    inner_join(assembly_norm, by = c("Participant_id", "tp_lab"), relationship = "many-to-many") %>%
    group_by(tp_lab, Infection_Status) %>%
    summarise(
        n_assemblies = n(),
        mean_contigs = round(mean(num_contigs, na.rm = TRUE), 2),
        median_contigs = median(num_contigs, na.rm = TRUE),
        mean_size_mb = round(mean(total_bases, na.rm = TRUE) / 1e6, 3),
        median_gc = round(median(gc_content, na.rm = TRUE), 2),
        asm_mix = paste(sort(table(assembler)), collapse = "|"),
        .groups = "drop"
    ) %>%
    arrange(tp_lab, Infection_Status)

write_csv(assembly_by_status, file.path(DIR_CLINICAL, "assembly_metrics_by_status.csv"))

p_asm_contigs <- ggplot(
    inner_join(episode_plot %>% select(Participant_id, tp_lab, Infection_Status) %>% distinct(),
        assembly_norm,
        by = c("Participant_id", "tp_lab"), relationship = "many-to-many"
    ),
    aes(Infection_Status, num_contigs)
) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = 0.2, alpha = 0.3) +
    facet_wrap(~tp_lab, scales = "free_x") +
    labs(title = "Assembly contig counts by status", x = "Status", y = "# contigs")

# 5. Waterfall
# ------------------------------------------------------------------------------
waterfall <- tibble::tibble(
    step = c("All episodes", "Culture positive", "S&S determinable", "S&S present", "Final UTIs"),
    n = c(
        nrow(episode_tbl),
        sum(episode_tbl$culture_pos_epi, na.rm = TRUE),
        sum(!is.na(episode_tbl$Sx_present_any)),
        sum(episode_tbl$Sx_present_any %in% TRUE, na.rm = TRUE),
        sum(episode_tbl$Infection_Status == "UTI")
    )
)
readr::write_csv(waterfall, file.path(DIR_CLINICAL, "waterfall_counts.csv"))
p_waterfall <- ggplot(waterfall, aes(x = factor(step, levels = step), y = n)) +
    geom_col() +
    geom_text(aes(label = n), vjust = -0.3, size = 3) +
    labs(title = "How rules narrow episodes down to UTIs", x = NULL, y = "Episodes") +
    theme_minimal(base_size = 11) +
    theme(axis.text.x = element_text(angle = 25, hjust = 1))

# 6. Save Plots
# ------------------------------------------------------------------------------
msg("Saving plots to %s", DIR_PLOTS_CLINICAL)
# safe_save_plot(file.path(DIR_PLOTS_CLINICAL, "status_distribution.png"), p_status, width = 8, height = 5)
safe_save_plot(file.path(DIR_PLOTS_CLINICAL, "trajectories_heatmap.png"), p_traj, width = 10, height = 12)
safe_save_plot(file.path(DIR_PLOTS_CLINICAL, "transitions_alluvial_or_heatmap.png"), p_transitions_flow, width = 12, height = 8)
safe_save_plot(file.path(DIR_PLOTS_CLINICAL, "assembly_contigs_boxplot.png"), p_asm_contigs, width = 10, height = 6)
safe_save_plot(file.path(DIR_PLOTS_CLINICAL, "waterfall_counts.png"), p_waterfall, width = 9, height = 5)

safe_pdf_begin(file.path(DIR_PLOTS_CLINICAL, "overview_plots.pdf"), width = 10, height = 6)
# print(p_status)
print(p_traj)
print(p_transitions_flow)
print(p_asm_contigs)
print(p_waterfall)
dev.off()

msg("00c_plot_clinical_summary.R completed.")

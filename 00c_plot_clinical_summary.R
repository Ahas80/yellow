#!/usr/bin/env Rscript
# ==============================================================================
# 00c_plot_clinical_summary.R
# ==============================================================================
#
# GOAL:
#   Visualise the clinical cohort: status distribution, assembly quality
#   metrics, waterfall/attrition counts. Produces the figures that describe
#   the study population in a thesis Results section or poster.
#
# WHY THIS SCRIPT EXISTS:
#   Before any genomic analysis, reviewers and supervisors need to see the
#   clinical landscape: how many UTI vs Not_UTI episodes, how they distribute
#   across batches, and whether assembly quality differs by primary status
#   (which could bias downstream VF detection).
#
# INPUTS:
#   - results/clinical/status_map.csv    (from 00b_classify_episodes.R)
#   - assembly_metadata.csv              (isolate-level metadata)
#
# OUTPUTS:
#   - plots/clinical/*.png, *.pdf        (status distribution figures)
#   - results/clinical/assembly_metrics_by_status.csv
#   - results/clinical/waterfall_counts.csv
# ==============================================================================
#   - Provides an overview of the cohort under the catheter-aware UTI definition.
# ==============================================================================

source("00_config.R")
source(here::here("R", "clinical_helpers.R"))
source("R/plot_helpers.R")

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
episode_tbl <- read_csv(FILE_STATUS_MAP, show_col_types = FALSE)
if (!all(c("UTI_Status", "UTI_binary") %in% names(episode_tbl))) {
    stop("status_map.csv lacks UTI_Status/UTI_binary. Rerun 00b_classify_episodes.R before plotting.")
}
assembly_df <- load_metadata() # Use helper from config

msg("Loaded status_map.csv (%d rows) and assembly_metadata.csv (%d rows)", nrow(episode_tbl), nrow(assembly_df))

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
    is_uti_event <- str_detect(tp_chr, regex("^UTI[-_ ]*[0-9]+$", ignore_case = TRUE))
    tp_num <- suppressWarnings(as.integer(str_extract(tp_chr, "\\d+")))
    tp_clean <- normalise_timepoint_preserve_events(tp_chr)
    routine_levels <- paste0("T", sort(unique(tp_num[!is.na(tp_num) & !is_uricult & !is_uti_event])))
    uti_levels <- sort(unique(tp_clean[is_uti_event]))
    tp_lab <- case_when(
        is_uricult ~ "Uricult",
        is_uti_event ~ tp_clean,
        !is.na(tp_num) ~ tp_clean,
        TRUE ~ "Unscheduled"
    )
    tibble(
        tp_num = tp_num,
        tp_lab = factor(tp_lab,
            levels = c(routine_levels, uti_levels, "Uricult", "Unscheduled")
        )
    )
}

# Prepare Data for Plotting
# ------------------------------------------------------------------------------
episode_tp_source <- if ("tp_lab" %in% names(episode_tbl)) episode_tbl$tp_lab else episode_tbl$Timepoint
episode_tbl <- episode_tbl %>%
    select(-any_of(c("tp_num", "tp_lab"))) %>%
    bind_cols(tp_norm2(episode_tp_source))

episode_plot <- episode_tbl %>%
    mutate(
        UTI_Status = as.character(UTI_Status),
        Status_plot = factor(UTI_Status, levels = c("UTI", "Not_UTI"))
    ) %>%
    filter(UTI_Status %in% c("UTI", "Not_UTI"))

status_order <- c("UTI", "Not_UTI")
status_levels_story <- c("Not_UTI", "UTI")

# 1. Status Distribution
# ------------------------------------------------------------------------------
# status_by_tp <- episode_plot %>%
#     count(tp_lab, Status_plot, name = "n") %>%
#     group_by(tp_lab) %>%
#     mutate(pct = 100 * n / sum(n)) %>%
#     ungroup() %>%
#     arrange(tp_lab, desc(n))

# p_status <- ggplot(status_by_tp, aes(tp_lab, n, fill = Status_plot)) +
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
p_traj <- ggplot(plot_df, aes(tp_lab, Participant_id, fill = Status_plot)) +
    geom_tile(color = "white") +
    scale_y_discrete(drop = FALSE) +
    scale_fill_uti_status() +
    labs(
        title = "Longitudinal UTI Status by Participant",
        subtitle = "Primary definition: catheter-aware S&S plus >=10^3 CFU culture support",
        x = "Timepoint",
        y = "Participant ID",
        fill = "Primary Status"
    ) +
    theme_minimal(base_size = 11) +
    theme(axis.text.y = element_text(size = 6))

# 3. Transitions
# ------------------------------------------------------------------------------
from_to <- episode_plot %>%
    filter(!is.na(tp_num)) %>%
    arrange(Participant_id, tp_num) %>%
    group_by(Participant_id) %>%
    mutate(Next_Status = dplyr::lead(Status_plot)) %>%
    ungroup() %>%
    filter(!is.na(Next_Status)) %>%
    count(From = Status_plot, To = Next_Status, name = "n") %>%
    complete(From = status_order, To = status_order, fill = list(n = 0)) %>%
    mutate(
        From = factor(From, levels = status_order),
        To = factor(To, levels = status_order)
    )

make_transitions_plot <- function(from_to, status_levels_story, status_cols = NULL) {
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
                    scale_fill_uti_status(drop = FALSE) +
                    labs(
                        title = "UTI Status Transitions Between Consecutive Timepoints",
                        x = "From Status",
                        y = "Count",
                        fill = "Previous Status"
                    ) +
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
        labs(
            title = "UTI Status Transitions (Heatmap)",
            x = "From Status",
            y = "To Status"
        ) +
        theme_minimal(base_size = 11)
}
p_transitions_flow <- make_transitions_plot(from_to, status_levels_story, NULL)

# 4. Assembly Metrics
# ------------------------------------------------------------------------------
assembly_norm <- assembly_df %>%
    transmute(
        Participant_id = as.character(Participant_id),
        Timepoint,
        tp_source = if ("tp_lab" %in% names(assembly_df)) as.character(tp_lab) else as.character(Timepoint),
        assembler, num_contigs, total_bases, gc_content
    ) %>%
    bind_cols(tp_norm2(.$tp_source)) %>%
    select(-tp_source)

episode_status_keys <- episode_plot %>%
    transmute(Participant_id = as.character(Participant_id), tp_lab, UTI_Status, Not_UTI_subgroup) %>%
    distinct()

assert_unique_keys(
    episode_status_keys,
    c("Participant_id", "tp_lab"),
    context = "00c clinical episode status keys",
    out_path = file.path(DIR_QC, "00c_duplicate_clinical_plot_keys.csv")
)

assembly_by_status <- episode_plot %>%
    transmute(Participant_id = as.character(Participant_id), tp_lab, UTI_Status, Not_UTI_subgroup) %>%
    distinct() %>%
    inner_join(assembly_norm, by = c("Participant_id", "tp_lab"), relationship = "one-to-many") %>%
    group_by(tp_lab, UTI_Status, Not_UTI_subgroup) %>%
    summarise(
        n_assemblies = n(),
        mean_contigs = round(mean(num_contigs, na.rm = TRUE), 2),
        median_contigs = median(num_contigs, na.rm = TRUE),
        mean_size_mb = round(mean(total_bases, na.rm = TRUE) / 1e6, 3),
        median_gc = round(median(gc_content, na.rm = TRUE), 2),
        asm_mix = paste(sort(table(assembler)), collapse = "|"),
        .groups = "drop"
    ) %>%
    arrange(tp_lab, UTI_Status, Not_UTI_subgroup)

write_csv(assembly_by_status, file.path(DIR_CLINICAL, "assembly_metrics_by_status.csv"))

assembly_plot_df <- inner_join(
    episode_status_keys,
        assembly_norm,
        by = c("Participant_id", "tp_lab"), relationship = "one-to-many"
    ) %>%
    filter(is.finite(num_contigs))

p_asm_contigs <- ggplot(
    assembly_plot_df,
    aes(UTI_Status, num_contigs, fill = UTI_Status)
) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = 0.2, alpha = 0.3) +
    facet_wrap(~tp_lab, scales = "free_x") +
    scale_fill_uti_status() +
    labs(
        title = "Assembly Quality Distribution by Primary UTI Status",
        x = "Primary UTI Status",
        y = "Number of Contigs"
    ) +
    theme_minimal(base_size = 11) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

# 5. Waterfall
# ------------------------------------------------------------------------------
waterfall <- tibble::tibble(
    step = c("All episodes", "Collection method known", "Culture supports UTI >=1e3", "Symptom-compatible", "Final UTI"),
    n = c(
        nrow(episode_tbl),
        sum(episode_tbl$catheter_rule %in% c("A_non_indwelling", "B_indwelling"), na.rm = TRUE),
        sum(episode_tbl$culture_supports_uti %in% TRUE, na.rm = TRUE),
        sum(episode_tbl$symptom_compatible_uti %in% TRUE, na.rm = TRUE),
        sum(episode_tbl$UTI_Status == "UTI", na.rm = TRUE)
    )
)
readr::write_csv(waterfall, file.path(DIR_CLINICAL, "waterfall_counts.csv"))
p_waterfall <- ggplot(waterfall, aes(x = factor(step, levels = step), y = n)) +
    geom_col() +
    geom_text(aes(label = n), vjust = -0.3, size = 3) +
    labs(
        title = "Catheter-Aware UTI Case Definition",
        subtitle = "Culture support uses the primary >=10^3 CFU/mL threshold where CFU is available",
        x = "Selection Step",
        y = "Number of Episodes"
    ) +
    theme_minimal(base_size = 11) +
    theme(axis.text.x = element_text(angle = 25, hjust = 1))

# 6. Definition diagnostics
# ------------------------------------------------------------------------------
legacy_col <- if ("Infection_Status_old" %in% names(episode_tbl)) "Infection_Status_old" else "Infection_Status"
reclassification <- episode_tbl %>%
    mutate(
        legacy_status = as.character(.data[[legacy_col]]),
        primary_status = as.character(UTI_Status),
        Not_UTI_subgroup_plot = ifelse(primary_status == "UTI", "UTI", as.character(Not_UTI_subgroup))
    ) %>%
    count(legacy_status, primary_status, Not_UTI_subgroup_plot, name = "n") %>%
    arrange(legacy_status, primary_status, Not_UTI_subgroup_plot)
write_csv(reclassification, file.path(DIR_CLINICAL, "uti_reclassification_plot_table.csv"))

p_reclass <- ggplot(reclassification, aes(x = legacy_status, y = Not_UTI_subgroup_plot, fill = n)) +
    geom_tile(colour = "white") +
    geom_text(aes(label = n), size = 3.2) +
    scale_fill_gradient(low = "#f7fbff", high = "#084081", name = "Episodes") +
    labs(
        title = "Legacy Status Reclassified Under Primary UTI Definition",
        subtitle = "Primary status is UTI vs Not_UTI; Not_UTI subgroups are descriptive",
        x = "Legacy ASB / UTI / Negative status",
        y = "Primary UTI or Not_UTI subgroup"
    ) +
    theme_minimal(base_size = 11) +
    theme(axis.text.x = element_text(angle = 25, hjust = 1))

subgroup_by_context <- episode_tbl %>%
    filter(UTI_Status == "Not_UTI") %>%
    mutate(
        Event_type_plot = if ("Event_type" %in% names(.)) as.character(Event_type) else episode_event_type(tp_lab),
        Batch_plot = if ("Batch" %in% names(.)) paste0("Batch ", Batch) else "Batch unknown"
    ) %>%
    count(Batch_plot, Event_type_plot, Not_UTI_subgroup, name = "n") %>%
    group_by(Batch_plot, Event_type_plot) %>%
    mutate(pct = 100 * n / sum(n)) %>%
    ungroup()
write_csv(subgroup_by_context, file.path(DIR_CLINICAL, "not_uti_subgroup_by_batch_event.csv"))

p_subgroup <- ggplot(subgroup_by_context,
                     aes(x = Event_type_plot, y = n, fill = Not_UTI_subgroup)) +
    geom_col(width = 0.72) +
    facet_wrap(~Batch_plot, scales = "free_y") +
    labs(
        title = "Not_UTI Subgroup Composition",
        subtitle = "Shows why the primary comparator is heterogeneous",
        x = "Sampling context",
        y = "Not_UTI episodes",
        fill = "Not_UTI subgroup"
    ) +
    theme_minimal(base_size = 10) +
    theme(axis.text.x = element_text(angle = 30, hjust = 1))

symptom_rule_counts <- episode_tbl %>%
    count(catheter_rule, symptom_rule_met, UTI_Status, name = "n") %>%
    arrange(catheter_rule, symptom_rule_met, UTI_Status)
write_csv(symptom_rule_counts, file.path(DIR_CLINICAL, "uti_symptom_rule_plot_table.csv"))

p_symptom_rule <- ggplot(symptom_rule_counts,
                         aes(x = symptom_rule_met, y = n, fill = UTI_Status)) +
    geom_col(position = "stack", width = 0.75) +
    facet_wrap(~catheter_rule, scales = "free_x") +
    scale_fill_uti_status() +
    labs(
        title = "Symptom-Rule Provenance For Primary UTI Classification",
        subtitle = "Catheter and non-catheter episodes are evaluated under different symptom rules",
        x = "Symptom rule result",
        y = "Episodes",
        fill = "Primary UTI status"
    ) +
    theme_minimal(base_size = 10) +
    theme(axis.text.x = element_text(angle = 35, hjust = 1))

cfu_plot_counts <- episode_tbl %>%
    mutate(
        cfu_1e3 = case_when(cfu_ge_1e3 %in% TRUE ~ ">=1e3", cfu_ge_1e3 %in% FALSE ~ "<1e3", TRUE ~ "unknown"),
        cfu_1e5 = case_when(cfu_ge_1e5 %in% TRUE ~ ">=1e5", cfu_ge_1e5 %in% FALSE ~ "<1e5", TRUE ~ "unknown"),
        threshold_profile = paste(cfu_1e3, cfu_1e5, sep = " / "),
        cfu_source_plot = coalesce(as.character(cfu_threshold_source), "missing_or_not_supporting")
    ) %>%
    count(threshold_profile, cfu_source_plot, UTI_Status, name = "n") %>%
    arrange(threshold_profile, cfu_source_plot, UTI_Status)
write_csv(cfu_plot_counts, file.path(DIR_CLINICAL, "uti_cfu_threshold_plot_table.csv"))

p_cfu <- ggplot(cfu_plot_counts,
                aes(x = threshold_profile, y = n, fill = UTI_Status)) +
    geom_col(width = 0.72) +
    facet_wrap(~cfu_source_plot, scales = "free_y") +
    scale_fill_uti_status() +
    labs(
        title = "CFU Threshold Provenance For Primary UTI Classification",
        subtitle = "Primary culture support uses >=10^3 CFU/mL or explicit +++ fallback when CFU is missing",
        x = "Parsed CFU threshold profile",
        y = "Episodes",
        fill = "Primary UTI status"
    ) +
    theme_minimal(base_size = 10) +
    theme(axis.text.x = element_text(angle = 30, hjust = 1))

uti_episode_inspection <- episode_tbl %>%
    filter(UTI_Status == "UTI") %>%
    select(any_of(c(
        "Participant_id", "Timepoint", "tp_lab", "Episode_ID", "Batch",
        "UTI_definition_version", "UTI_Status", "UTI_binary",
        "urine_collection_method_raw", "urine_collection_method_norm",
        "catheter_rule", "symptom_rule_met", "symptom_compatible_uti",
        "culture_supports_uti", "cfu_raw", "cfu_raw_parsed",
        "cfu_ge_1e3", "cfu_ge_1e5", "cfu_threshold_source",
        "dysuria_present", "urgency_present", "frequency_present",
        "incontinence_present", "pus_present", "flankpain_present",
        "fever_present", "rigors_present", "delirium_present",
        "UTI_classification_confidence", "UTI_classification_reason"
    ))) %>%
    arrange(Participant_id, tp_lab)
write_csv(uti_episode_inspection, file.path(DIR_CLINICAL, "uti_episode_inspection_table.csv"))

# 6. Save Plots
# ------------------------------------------------------------------------------
msg("Saving plots to %s", DIR_PLOTS_CLINICAL)
# safe_save_plot(file.path(DIR_PLOTS_CLINICAL, "status_distribution.png"), p_status, width = 8, height = 5)
safe_save_plot(file.path(DIR_PLOTS_CLINICAL, "trajectories_heatmap.png"), p_traj, width = 10, height = 12)
safe_save_plot(file.path(DIR_PLOTS_CLINICAL, "transitions_alluvial_or_heatmap.png"), p_transitions_flow, width = 12, height = 8)
safe_save_plot(file.path(DIR_PLOTS_CLINICAL, "assembly_contigs_boxplot.png"), p_asm_contigs, width = 10, height = 6)
safe_save_plot(file.path(DIR_PLOTS_CLINICAL, "waterfall_counts.png"), p_waterfall, width = 9, height = 5)
safe_save_plot(file.path(DIR_PLOTS_CLINICAL, "uti_reclassification_heatmap.png"), p_reclass, width = 8.5, height = 5.5)
safe_save_plot(file.path(DIR_PLOTS_CLINICAL, "not_uti_subgroup_by_batch_event.png"), p_subgroup, width = 11, height = 7)
safe_save_plot(file.path(DIR_PLOTS_CLINICAL, "uti_symptom_rule_provenance.png"), p_symptom_rule, width = 10, height = 6)
safe_save_plot(file.path(DIR_PLOTS_CLINICAL, "uti_cfu_threshold_provenance.png"), p_cfu, width = 11, height = 6)

safe_pdf_begin(file.path(DIR_PLOTS_CLINICAL, "overview_plots.pdf"), width = 10, height = 6)
# print(p_status)
print(p_traj)
print(p_transitions_flow)
print(p_asm_contigs)
print(p_waterfall)
print(p_reclass)
print(p_subgroup)
print(p_symptom_rule)
print(p_cfu)
dev.off()

msg("00c_plot_clinical_summary.R completed.")

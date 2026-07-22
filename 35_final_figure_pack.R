#!/usr/bin/env Rscript

# ==============================================================================
# 35_final_figure_pack.R
# ==============================================================================
# Canonical thesis-figure layer for the selected QC-passing Longcycler cohort.
# This script runs only after RQ01-RQ10 have completed. It consumes validated
# tables, derives every displayed denominator, and never reclassifies episodes.
# Legacy ASB/UTI/Negative figures remain descriptive and are not read here.
# ==============================================================================

source("00_config.R")
source("R/plot_helpers.R")
source("R/wgs_helpers.R")

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
  library(scales)
  library(patchwork)
  library(ggrepel)
  library(digest)
  library(ape)
  library(ggtree)
})

set.seed(20260714)
msg("Starting canonical thesis figure pack")

DIR_FINAL <- file.path(DIR_PLOTS, "final")
DIR_SUPP <- file.path(DIR_FINAL, "supplementary")
DIR_FIG_AUDIT <- file.path(DIR_RESULTS, "figure_audit")
ensure_dir(DIR_FINAL)
ensure_dir(DIR_SUPP)
ensure_dir(DIR_FIG_AUDIT)

paths <- c(
  source_status = file.path(DIR_CLINICAL, "status_map.csv"),
  cohort = FILE_ANALYSIS_CLINICAL_COHORT,
  qc = FILE_CANONICAL_ASSEMBLY_SELECTION,
  mlst = file.path(DIR_MLST, "mlst_provider_preferred.csv"),
  vf_ready = FILE_VF_READY,
  vf_pa = file.path(DIR_VF, "vf_pa_all.csv"),
  gene_map = file.path(DIR_VF, "gene_map.csv"),
  score_table = file.path(DIR_VF, "vf_score_table.csv"),
  score_models = file.path(DIR_RESULTS, "statistical_sensitivity", "score_glmm_sensitivity.csv"),
  score_collapsed = file.path(DIR_RESULTS, "statistical_sensitivity", "participant_collapsed_score_values.csv"),
  gene_models = file.path(DIR_RESULTS, "models", "gwas_multivariable_glmm.csv"),
  fisher = file.path(DIR_RESULTS, "models", "gwas_univariable_stats.csv"),
  timelines = file.path(DIR_RESULTS, "longitudinal", "participant_timelines.csv"),
  transitions = file.path(DIR_RESULTS, "longitudinal", "longcycler_transitions.csv"),
  casebook = file.path(DIR_RESULTS, "mechanism", "not_uti_to_uti_casebook.csv"),
  mechanism_summary = file.path(DIR_RESULTS, "mechanism", "transition_mechanism_summary.csv"),
  pairwise = file.path(DIR_RESULTS, "strain_compare", "pairwise_metrics.csv"),
  rq_marker = file.path(DIR_RESULTS, "research_questions", "RUN_COMPLETE.txt"),
  rq01_threshold = file.path(DIR_RESULTS, "research_questions", "RQ01", "threshold_sensitivity.csv"),
  rq03_cases = file.path(DIR_RESULTS, "research_questions", "RQ03", "deidentified_case_matrix.csv"),
  rq06_pairs = file.path(DIR_RESULTS, "research_questions", "RQ06", "rq06_adjacent_pair_vf_changes.csv"),
  rq07_inference = file.path(DIR_RESULTS, "research_questions", "RQ07", "rq07_event_sample_inference.csv"),
  rq07_paired = file.path(DIR_RESULTS, "research_questions", "RQ07", "rq07_nearest_paired_resident_values.csv"),
  module_changes = file.path(DIR_RESULTS, "statistical_sensitivity", "not_uti_to_uti_module_change_matrix.csv"),
  pcoa = file.path(DIR_VF, "vf_pcoa_jaccard_coordinates.csv"),
  near_miss = file.path(DIR_RESULTS, "audit", "uti_not_uti_near_miss_rows.csv"),
  leave_one = file.path(DIR_VF, "uti_not_uti_leave_one_uti_out.csv"),
  plasmid_manifest = file.path(DIR_RESULTS, "plasmids", "plasmidfinder_input_manifest.csv"),
  replicons = file.path(DIR_MLST, "plasmid_replicons_long.csv"),
  mob_profiles = file.path(
    DIR_RESULTS, "plasmids", "mob_suite",
    "episode_plasmid_profiles.csv"
  ),
  plasmid_mechanism_profiles = file.path(
    DIR_RESULTS, "plasmids", "mob_suite",
    "episode_mechanism_profiles.csv"
  ),
  plasmid_focused = file.path(
    DIR_RESULTS, "plasmids", "mob_suite",
    "not_uti_to_uti_plasmid_metrics_9.csv"
  ),
  plasmid_location_validation = file.path(
    DIR_RESULTS, "plasmids", "mob_suite",
    "plasmid_gene_location_validation.csv"
  ),
  amr = file.path(DIR_RESULTS, "amr", "episode_amr_profiles.csv"),
  tree = file.path(DIR_RESULTS, "wgs", "core", "core_genome.tree"),
  tree_map = file.path(DIR_RESULTS, "wgs", "core", "core_snp_sample_map.csv"),
  variants = file.path(DIR_FIG_AUDIT, "reference_aware_variants.csv"),
  variant_validation = file.path(DIR_FIG_AUDIT, "reference_aware_variants_validation.csv")
)

required <- setdiff(names(paths), c("variants", "variant_validation"))
missing_required <- paths[required][!file.exists(paths[required])]
if (length(missing_required)) {
  stop("Missing canonical figure input(s): ", paste(missing_required, collapse = "; "), call. = FALSE)
}
if (any(str_detect(paths[required], regex("_OLD|old_asb|legacy", ignore_case = TRUE)))) {
  stop("A canonical figure input points to a legacy result.", call. = FALSE)
}

read_current <- function(path) read_csv(path, show_col_types = FALSE)
assert_unique <- function(x, keys, label) {
  if (any(!keys %in% names(x))) stop(label, " lacks key column(s): ", paste(setdiff(keys, names(x)), collapse = ", "))
  dup <- x %>% count(across(all_of(keys)), name = ".n") %>% filter(.data$.n > 1L)
  if (nrow(dup)) stop(label, " contains ", nrow(dup), " duplicated key(s).", call. = FALSE)
  invisible(x)
}
assert_value <- function(observed, expected, label) {
  if (!identical(as.integer(observed), as.integer(expected))) {
    stop(label, ": observed ", observed, ", expected ", expected, call. = FALSE)
  }
}
status_label <- function(x) recode_operational_uti_status(x, strict = TRUE, as_factor = TRUE)
status_colours <- c("Not UTI" = "#0072B2", "UTI" = "#D55E00", "Unknown" = "#CCCCCC")
status_fill <- function(...) scale_fill_operational_status(..., reader_facing = TRUE)
status_colour <- function(...) scale_colour_operational_status(..., reader_facing = TRUE)
reader_event <- function(x) recode(as.character(x), Routine = "Routine surveillance", UTI_event = "UTI-related sampling", .default = "Other sampling")
normalise_st <- function(x) {
  out <- str_trim(as.character(x))
  out[out %in% c("", "-", "NA", "N/A", "Unknown", "unknown")] <- NA_character_
  out
}
case_labels <- function(casebook) {
  casebook %>%
    distinct(case_id, Participant_id) %>%
    arrange(.data$case_id) %>%
    mutate(Case_Label = sprintf("Case %02d", row_number()))
}
short_score <- c(
  total_vf_count_all = "All detected VF genes",
  total_vf_count_curated = "Curated VF genes",
  upec_system_count = "UPEC-associated systems",
  expec_marker_count = "ExPEC-like markers"
)
mechanism_labels <- c(
  same_strain_stable_profile = "Same strain, stable profile",
  same_strain_genomic_change = "Same strain, profile change",
  strain_replacement = "Strain replacement",
  uncertain = "Uncertain",
  missing_wgs_endpoint = "Genomic endpoint unavailable"
)
mechanism_colours <- c(
  same_strain_stable_profile = "#009E73",
  same_strain_genomic_change = "#E69F00",
  strain_replacement = "#CC79A7",
  uncertain = "#6B7280",
  missing_wgs_endpoint = "#BDBDBD"
)

cohort <- read_current(paths[["cohort"]]) %>%
  mutate(
    Participant_id = as.character(.data$Participant_id),
    tp_lab = normalise_timepoint_preserve_events(.data$tp_lab),
    status_display = status_label(.data$UTI_Status)
  )
assert_unique(cohort, c("Participant_id", "tp_lab"), "Selected Longcycler cohort")
assert_value(nrow(cohort), 532L, "Selected Longcycler episodes")
assert_value(n_distinct(cohort$Participant_id), 161L, "Selected participants")
assert_value(sum(cohort$UTI_Status == "UTI"), 16L, "Operational UTI episodes")
assert_value(sum(cohort$UTI_Status == "Not_UTI"), 516L, "Operational Not UTI episodes")

source_status <- read_current(paths[["source_status"]]) %>%
  filter(.data$analysis_include_primary %in% TRUE) %>%
  mutate(Participant_id = as.character(.data$Participant_id), status_display = status_label(.data$UTI_Status))
assert_unique(source_status, c("Participant_id", "tp_lab"), "Primary clinical source cohort")
assert_value(nrow(source_status), 583L, "Primary clinical source episodes")
assert_value(n_distinct(source_status$Participant_id), 166L, "Primary clinical source participants")

qc <- read_current(paths[["qc"]]) %>% mutate(Participant_id = as.character(.data$Participant_id))
assert_unique(qc, c("Assembly_ID"), "Longcycler candidate QC table")
assert_value(nrow(qc), 592L, "Longcycler candidates")
assert_value(sum(qc$selected_canonical %in% TRUE & qc$QC_PASS %in% TRUE), 532L, "Selected QC-passing assemblies")

mlst <- read_current(paths[["mlst"]]) %>%
  mutate(Participant_id = as.character(.data$Participant_id), tp_lab = normalise_timepoint_preserve_events(.data$tp_lab))
assert_unique(mlst, c("Participant_id", "tp_lab"), "Provider-preferred MLST")
assert_value(nrow(mlst), 532L, "MLST denominator")

vf_ready <- read_current(paths[["vf_ready"]]) %>%
  mutate(Participant_id = as.character(.data$Participant_id), tp_lab = normalise_timepoint_preserve_events(.data$tp_lab))
assert_unique(vf_ready, c("Participant_id", "tp_lab"), "VF-ready table")
assert_value(nrow(vf_ready), 532L, "VF-ready denominator")

score_table <- read_current(paths[["score_table"]]) %>%
  mutate(Participant_id = as.character(.data$Participant_id), tp_lab = normalise_timepoint_preserve_events(.data$tp_lab)) %>%
  select(-any_of(c("Infection_Status", "UTI_Status"))) %>%
  left_join(cohort %>% select(Participant_id, tp_lab, UTI_Status, status_display),
            by = c("Participant_id", "tp_lab"), relationship = "one-to-one")
assert_unique(score_table, c("Participant_id", "tp_lab"), "VF score table")
assert_value(nrow(score_table), 532L, "VF score denominator")
if (anyNA(score_table$UTI_Status)) stop("VF score table has unmatched clinical status.")

transitions <- read_current(paths[["transitions"]]) %>% mutate(Participant_id = as.character(.data$Participant_id))
assert_unique(transitions, c("Participant_id", "tp_from", "tp_to"), "Canonical adjacent transitions")
assert_value(nrow(transitions), 371L, "Adjacent transition count")
assert_value(n_distinct(transitions$Participant_id), 139L, "Participants with adjacent transitions")
assert_value(sum(transitions$status_from == "Not_UTI" & transitions$status_to == "UTI"), 9L, "Not UTI-to-UTI transitions")
assert_value(sum(transitions$status_from == "Not_UTI" & transitions$status_to == "UTI" & transitions$TotalSNPs <= SAME_STRAIN_SNP_THRESHOLD), 5L, "Not UTI-to-UTI transitions at or below 25 SNPs")

pairwise <- read_current(paths[["pairwise"]])
assert_value(nrow(pairwise), 893L, "Direct within-participant pair count")
if (any(pairwise$within_participant %in% FALSE, na.rm = TRUE)) stop("Pairwise table contains cross-participant pairs.")

casebook <- read_current(paths[["casebook"]]) %>% mutate(Participant_id = as.character(.data$Participant_id))
assert_unique(casebook, c("Participant_id", "from_tp", "to_tp"), "Not UTI-to-UTI casebook")
assert_value(nrow(casebook), 9L, "Not UTI-to-UTI casebook size")
case_key <- case_labels(casebook)
casebook <- casebook %>% left_join(case_key, by = c("case_id", "Participant_id"), relationship = "many-to-one")

rq_status_files <- file.path(DIR_RESULTS, "research_questions", sprintf("RQ%02d", 1:10), "analysis_status.csv")
if (any(!file.exists(rq_status_files))) stop("Final figures require all ten per-question completion files.")
rq_status <- bind_rows(lapply(rq_status_files, read_current))
if (nrow(rq_status) != 10L || any(rq_status$status != "complete")) stop("At least one research question is not complete.")

manifest <- tibble()
register_figure <- function(figure_id, figure_number, title, scientific_question, figure_class,
                            plot, width, height, caption, source_inputs, statistical_method,
                            caveat, unit, filters, visual_encodings, multiplicity,
                            validation_status = "validated") {
  out_dir <- if (identical(figure_class, "Main")) DIR_FINAL else DIR_SUPP
  stub <- file.path(out_dir, figure_id)
  save_result <- save_ruti_figure(
    plot = plot,
    filename = stub,
    width = width,
    height = height,
    dpi = 300,
    figure_id = figure_id,
    manifest_path = file.path(DIR_FIG_AUDIT, "save_ruti_figure_manifest.csv"),
    metadata = list(
      figure_number = figure_number,
      title = title,
      scientific_question = scientific_question,
      figure_class = figure_class,
      caption = caption,
      source_inputs = source_inputs,
      statistical_method = statistical_method,
      caveat = caveat,
      unit = unit,
      filters = filters,
      visual_encodings = visual_encodings,
      multiplicity = multiplicity,
      validation_status = validation_status
    )
  )
  manifest <<- bind_rows(manifest, tibble(
    figure_id, figure_number, title, scientific_question, figure_class,
    png_path = paste0(stub, ".png"), pdf_path = paste0(stub, ".pdf"),
    width_in = width, height_in = height, dpi = 300L,
    caption, source_inputs = paste(source_inputs, collapse = "; "),
    statistical_method, caveat, unit, filters, visual_encodings, multiplicity,
    validation_status
  ))
  invisible(save_result)
}

panel_annotation_theme <- theme(
  plot.title = element_text(face = "bold", size = 15),
  plot.subtitle = element_text(size = 10.5, colour = "grey25"),
  plot.caption = element_text(size = 8.5, colour = "grey35", hjust = 0),
  plot.margin = margin(8, 10, 8, 8)
)

# ----------------------------------------------------------------------------
# Fig01. Cohort and analytical denominators
# ----------------------------------------------------------------------------
cohort_flow <- tibble(
  stage = factor(c("Primary clinical cohort", "Selected Longcycler cohort"),
                 levels = rev(c("Primary clinical cohort", "Selected Longcycler cohort"))),
  episodes = c(nrow(source_status), nrow(cohort)),
  participants = c(n_distinct(source_status$Participant_id), n_distinct(cohort$Participant_id))
)
p01a <- cohort_flow %>%
  pivot_longer(c("episodes", "participants"), names_to = "unit", values_to = "n") %>%
  mutate(unit = recode(.data$unit, episodes = "Episodes", participants = "Participants")) %>%
  ggplot(aes(.data$n, .data$stage, colour = .data$unit)) +
  geom_segment(aes(x = 0, xend = .data$n, yend = .data$stage), linewidth = 1.1, colour = "grey78") +
  geom_point(size = 4) +
  geom_text(aes(label = .data$n), nudge_x = -25, hjust = 1, fontface = "bold", size = 3.3) +
  facet_wrap(~unit, scales = "free_x") +
  scale_colour_manual(values = c(Episodes = "#4C78A8", Participants = "#7A5195"), guide = "none") +
  scale_x_continuous(expand = expansion(mult = c(0, .14))) +
  labs(title = "Selection into the genomic cohort", x = "Count", y = NULL) +
  theme_ruti_publication(9.5)

status_counts <- bind_rows(
  source_status %>% count(.data$status_display, name = "n") %>% mutate(cohort = "Primary clinical cohort"),
  cohort %>% count(.data$status_display, name = "n") %>% mutate(cohort = "Selected Longcycler cohort")
) %>%
  group_by(.data$cohort) %>% mutate(percent = 100 * .data$n / sum(.data$n)) %>% ungroup()
p01b <- ggplot(status_counts, aes(.data$cohort, .data$n, fill = .data$status_display)) +
  geom_col(width = .65, colour = "white") +
  geom_text(aes(label = .data$n), position = position_stack(vjust = .5), colour = "white", fontface = "bold", size = 3.2) +
  status_fill(name = "Operational UTI status") +
  labs(title = "Operational-status composition", x = NULL, y = "Number of episodes") +
  theme_ruti_publication(9.5)

not_uti_subgroups <- cohort %>%
  filter(.data$UTI_Status == "Not_UTI") %>%
  mutate(
    legacy_type = factor(.data$Infection_Status_legacy, levels = c("UTI", "Negative", "ASB")),
    subgroup = factor(
      recode(.data$Infection_Status_legacy,
             ASB = "Legacy ASB",
             Negative = "Legacy Negative",
             UTI = "Legacy UTI, reclassified"),
      levels = c("Legacy UTI, reclassified", "Legacy Negative", "Legacy ASB")
    )
  ) %>%
  count(.data$legacy_type, .data$subgroup, sort = TRUE)
assert_value(sum(not_uti_subgroups$n), 516L, "Legacy composition of selected Not UTI episodes")
p01c <- ggplot(not_uti_subgroups, aes(.data$n, .data$subgroup, fill = .data$legacy_type)) +
  geom_col(width = .68) +
  geom_text(aes(label = .data$n), hjust = -0.2, fontface = "bold", size = 3.1) +
  scale_fill_clinical_episode(name = "Legacy clinical episode type", guide = "none") +
  scale_x_continuous(expand = expansion(mult = c(0, .12))) +
  labs(title = "Composition of operational Not UTI", subtitle = "Legacy episode type; descriptive context only", x = "Number of episodes", y = NULL) +
  theme_ruti_publication(9.5)

fig01 <- (p01a | p01b | p01c) +
  plot_annotation(
    title = "Cohort selection and analytical denominators",
    subtitle = "The canonical analysis contains 532 episodes from 161 participants: 16 UTI and 516 heterogeneous Not UTI episodes.",
    caption = "Counts are derived from the primary clinical status map and the selected QC-passing Longcycler manifest.",
    tag_levels = "A", theme = panel_annotation_theme
  ) & theme(legend.position = "bottom")
cap01 <- paste0(
  "Cohort selection and analytical denominators. Panel A shows episode- and participant-level attrition from 583 primary clinical episodes among 166 participants to 532 selected QC-passing Longcycler episodes among 161 participants. ",
  "Panel B shows operational UTI status before and after genomic selection; the selected cohort contains 16 UTI and 516 Not UTI episodes. Panel C provides descriptive legacy context for the 516 operational Not UTI episodes: 416 were formerly classified as ASB, 83 as Negative, and 17 as UTI. Bars encode counts; colours retain the canonical operational and legacy episode-type mappings. No inferential test is shown. Repeated episodes from the same participant are descriptive observations and are not independent."
)
register_figure("Fig01_cohort_and_denominators", "Fig01", "Cohort selection and analytical denominators",
                "Which episodes and participants enter each analytical denominator?", "Main", fig01, 15, 6.2, cap01,
                paths[c("source_status", "cohort", "qc")], "Descriptive denominator audit",
                "Not UTI is heterogeneous; episode counts are not independent participant counts.",
                "Episode and participant", "Primary clinical inclusion and selected QC-passing Longcycler assemblies",
                "Bar length and labels are counts; fill denotes operational UTI status or, in panel C, labelled legacy episode type.", "None")

# ----------------------------------------------------------------------------
# Fig02. WGS quality control
# ----------------------------------------------------------------------------
qc_plot <- qc %>%
  mutate(
    qc_state = factor(if_else(.data$QC_PASS, "Pass", "Fail"), levels = c("Pass", "Fail")),
    selection_state = factor(if_else(.data$selected_canonical, "Selected", "Not selected"), levels = c("Selected", "Not selected")),
    assembly_mb = .data$total_bp / 1e6,
    outlier_rank = rank(abs(.data$assembly_mb - 5), ties.method = "first", na.last = "keep"),
    outlier_label = NA_character_
  )
fail_labs <- qc_plot %>% filter(.data$qc_state == "Fail") %>% arrange(desc(abs(.data$assembly_mb - 5))) %>% slice_head(n = 8) %>% mutate(outlier_label = sprintf("QC outlier %02d", row_number()))
qc_plot <- qc_plot %>% select(-outlier_label) %>% left_join(fail_labs %>% select(Assembly_ID, outlier_label), by = "Assembly_ID", relationship = "one-to-one")
qc_cols <- c(Pass = "#0072B2", Fail = "#D55E00")
selection_fills <- c(Selected = "#555555", `Not selected` = "white")
p02a <- ggplot(qc_plot, aes(.data$n_contigs, .data$N50 / 1000,
                            colour = .data$qc_state, shape = .data$qc_state,
                            fill = .data$selection_state)) +
  geom_vline(xintercept = 200, linetype = "dashed", colour = "grey35") +
  geom_hline(yintercept = 20, linetype = "dashed", colour = "grey35") +
  geom_point(alpha = .78, size = 2.15, stroke = .7) +
  scale_y_log10(labels = label_number()) +
  scale_colour_manual(values = qc_cols, name = "Assembly QC") +
  scale_shape_manual(values = c(Pass = 21, Fail = 24), name = "Assembly QC") +
  scale_fill_manual(values = selection_fills, name = "Canonical selection") +
  guides(fill = guide_legend(override.aes = list(shape = 21, colour = "grey20", alpha = 1, size = 2.8))) +
  labs(title = "Contiguity metrics", subtitle = "N50 axis is logarithmic; dashed lines: <=200 contigs and N50 >=20 kb",
       x = "Number of contigs", y = "Assembly N50 (kb; logarithmic scale)") +
  theme_ruti_publication(9.5)
p02b <- ggplot(qc_plot, aes(.data$assembly_mb, .data$GC_pct,
                            colour = .data$qc_state, shape = .data$qc_state,
                            fill = .data$selection_state)) +
  geom_vline(xintercept = c(4, 6), linetype = "dashed", colour = "grey35") +
  geom_point(alpha = .78, size = 2.15, stroke = .7) +
  ggrepel::geom_text_repel(data = qc_plot %>% filter(!is.na(.data$outlier_label)), aes(label = .data$outlier_label),
                           size = 2.6, seed = 20260714, min.segment.length = 0, max.overlaps = Inf, show.legend = FALSE) +
  scale_colour_manual(values = qc_cols, name = "Assembly QC") +
  scale_shape_manual(values = c(Pass = 21, Fail = 24), name = "Assembly QC") +
  scale_fill_manual(values = selection_fills, name = "Canonical selection") +
  guides(fill = guide_legend(override.aes = list(shape = 21, colour = "grey20", alpha = 1, size = 2.8))) +
  labs(title = "Assembly size and GC content", subtitle = "Dashed lines: accepted assembly-size range 4-6 Mb; no GC threshold was imposed",
       x = "Assembly length (Mb)", y = "GC content (%)") +
  theme_ruti_publication(9.5)
fig02 <- (p02a | p02b) + plot_layout(guides = "collect") +
  plot_annotation(
    title = "Whole-genome assembly quality control",
    subtitle = sprintf("%d Longcycler candidates; %d passed all implemented thresholds; %d were selected for analysis.", nrow(qc_plot), sum(qc_plot$QC_PASS), sum(qc_plot$selected_canonical & qc_plot$QC_PASS)),
    caption = "All candidates were assembled with Longcycler. Failure labels follow an objective extreme-size rule and disclose no isolate identifier.",
    tag_levels = "A", theme = panel_annotation_theme
  ) & theme(legend.position = "bottom")
cap02 <- paste0(
  "Whole-genome assembly quality control for 592 Longcycler candidate assemblies. Panel A plots N50 on a logarithmic kilobase axis against contig count with the implemented thresholds of N50 >=20 kb and <=200 contigs. Panel B plots assembly length against GC content with the implemented 4-6 Mb assembly-length limits; GC content is shown diagnostically and was not thresholded. ",
  "Points are assemblies; colour and shape redundantly indicate pass/fail, while filled versus open symbols distinguish the 532 canonically selected assemblies from candidates not selected. Eight most extreme failing assembly lengths are labelled using deidentified QC labels. No inferential test is shown. Repeated episodes from the same participant can contribute multiple assemblies and are not independent."
)
register_figure("Fig02_wgs_quality_control", "Fig02", "Whole-genome assembly quality control",
                "Do candidate assemblies satisfy the implemented quality thresholds?", "Main", fig02, 13.5, 6.4, cap02,
                paths[c("qc")], "Descriptive quality-control thresholds",
                "GC content was diagnostic only; repeated episodes from the same participant are not independent.",
                "Assembly", "All primary-analysis Longcycler candidates",
                "Points are assemblies; colour and shape are pass/fail; symbol fill is canonical selection; dashed lines are implemented thresholds; N50 is logarithmic.", "None")

# ----------------------------------------------------------------------------
# Fig03. Sequence types and provenance
# ----------------------------------------------------------------------------
st_data <- mlst %>%
  select(Participant_id, tp_lab, ST, ST_source, ST_provider, ST_local) %>%
  left_join(cohort %>% select(Participant_id, tp_lab, UTI_Status, status_display),
            by = c("Participant_id", "tp_lab"), relationship = "one-to-one") %>%
  mutate(ST_norm = normalise_st(.data$ST))
st_keep <- st_data %>%
  count(.data$ST_norm, .data$UTI_Status, name = "n") %>%
  group_by(.data$ST_norm) %>%
  summarise(n_all = sum(.data$n), n_uti = sum(.data$n[.data$UTI_Status == "UTI"]), .groups = "drop") %>%
  filter(!is.na(.data$ST_norm), .data$n_all >= 10L | .data$n_uti >= 2L) %>% pull("ST_norm")
st_data <- st_data %>%
  mutate(ST_group = case_when(
    is.na(.data$ST_norm) ~ "Missing / non-typable",
    .data$ST_norm %in% st_keep ~ paste0("ST", .data$ST_norm),
    TRUE ~ "Other typed ST"
  ))
st_summary <- st_data %>%
  count(.data$ST_group, .data$status_display, name = "n") %>%
  group_by(.data$status_display) %>% mutate(percent = 100 * .data$n / sum(.data$n)) %>% ungroup()
st_order <- st_summary %>% group_by(.data$ST_group) %>% summarise(n = sum(.data$n), .groups = "drop") %>% arrange(.data$n) %>% pull("ST_group")
st_summary <- st_summary %>%
  mutate(
    ST_group = factor(.data$ST_group, levels = st_order),
    label_y = as.numeric(.data$ST_group) + if_else(.data$status_display == "UTI", .14, -.14)
  )
p03a <- ggplot(st_summary, aes(.data$percent, .data$ST_group, colour = .data$status_display)) +
  geom_line(aes(group = .data$ST_group), colour = "grey75", linewidth = .7) +
  geom_point(aes(shape = .data$status_display), size = 3) +
  geom_text(aes(y = .data$label_y, label = paste0("n=", .data$n)), hjust = -.25, size = 2.6, show.legend = FALSE) +
  status_colour(name = "Operational UTI status") +
  scale_shape_manual(values = c("Not UTI" = 16, "UTI" = 17), name = "Operational UTI status") +
  scale_x_continuous(labels = label_number(suffix = "%"), expand = expansion(mult = c(.02, .18))) +
  labs(title = "Sequence-type composition", subtitle = "Display rule: >=10 episodes overall or >=2 UTI episodes; remaining typed calls are grouped",
       x = "Within-status prevalence", y = NULL) + theme_ruti_publication(9)
prov_summary <- st_data %>%
  mutate(provenance = recode(.data$ST_source,
                             provider_qc95 = "Provider call (QC >=95%)",
                             local_fallback_provider_missing = "Local fallback",
                             missing = "Missing / non-typable",
                             .default = "Other / unavailable")) %>%
  count(.data$status_display, .data$provenance, name = "n") %>%
  group_by(.data$status_display) %>% mutate(percent = .data$n / sum(.data$n)) %>% ungroup()
prov_cols <- c("Provider call (QC >=95%)" = "#0072B2", "Local fallback" = "#E69F00", "Missing / non-typable" = "#BDBDBD", "Other / unavailable" = "#6B7280")
p03b <- ggplot(prov_summary, aes(.data$status_display, .data$percent, fill = .data$provenance)) +
  geom_col(width = .65, colour = "white") +
  geom_text(aes(label = as.character(.data$n)), position = position_stack(vjust = .5), size = 2.8) +
  scale_y_continuous(labels = label_percent()) + scale_fill_manual(values = prov_cols, name = "ST provenance") +
  labs(title = "Typing provenance", subtitle = "Counts are episodes, not participants", x = "Operational UTI status", y = "Episodes within status") +
  theme_ruti_publication(9)
fig03 <- (p03a | p03b) + plot_layout(widths = c(.66, .34), guides = "collect") +
  plot_annotation(
    title = "Sequence-type distribution and provenance",
    subtitle = "Provider calls are preferred; local calls are labelled as fallback rather than treated as equivalent provenance.",
    caption = "Lineage comparisons are descriptive because UTI counts are sparse and repeated episodes can share participant and lineage.",
    tag_levels = "A", theme = panel_annotation_theme
  ) & theme(legend.position = "bottom")
cap03 <- paste0(
  "Sequence-type distribution and provenance across 532 selected episodes (516 Not UTI; 16 UTI) from 161 participants. Panel A compares within-status prevalence for sequence types observed in at least 10 episodes overall or at least two UTI episodes; remaining typed calls are grouped as Other typed ST and missing/non-typable calls remain separate. Points show percentages, connecting lines aid comparison, and labels give episode counts. ",
  "Panel B shows provider-derived, local-fallback, and missing typing provenance. In Panel A, colour and shape redundantly encode operational UTI status. Provider-only calls should be used for lineage claims. The display is descriptive; repeated episodes and lineage structure preclude treating episode counts as independent evidence of UTI risk."
)
register_figure("Fig03_sequence_type_distribution", "Fig03", "Sequence-type distribution and provenance",
                "How are sequence types distributed by operational status, and where did calls originate?", "Main", fig03, 14.5, 7.2, cap03,
                paths[c("mlst", "cohort")], "Descriptive prevalence and provenance",
                "Provider/local calls have distinct provenance; lineage claims require provider-only sensitivity and participant-aware analysis.",
                "Episode", "Selected cohort; displayed ST rule >=10 episodes overall or >=2 UTI episodes",
                "Dots are within-status percentages; labels are episode counts; stacked bars encode ST provenance.", "None")

# ----------------------------------------------------------------------------
# Fig04. VF score distributions and paired participant summaries
# ----------------------------------------------------------------------------
score_long <- score_table %>%
  select(Participant_id, tp_lab, status_display, all_of(names(short_score))) %>%
  pivot_longer(all_of(names(short_score)), names_to = "score", values_to = "value") %>%
  mutate(score_label = factor(short_score[.data$score], levels = unname(short_score)),
         status_display = factor(.data$status_display, levels = c("Not UTI", "UTI")))
p04a <- ggplot(score_long, aes(.data$status_display, .data$value, fill = .data$status_display, colour = .data$status_display)) +
  geom_violin(alpha = .16, linewidth = .35, width = .82, trim = FALSE) +
  geom_boxplot(width = .24, outlier.shape = NA, alpha = .5, linewidth = .45) +
  geom_point(position = position_jitter(width = .14, height = 0, seed = 20260714), alpha = .22, size = .8) +
  stat_summary(fun = median, geom = "point", shape = 23, fill = "white", colour = "black", size = 2.5) +
  facet_wrap(~score_label, scales = "free_y", nrow = 2) +
  status_fill(name = "Operational UTI status") + status_colour(name = "Operational UTI status") +
  scale_x_discrete(labels = c("Not UTI" = "Not UTI\n516 episodes", "UTI" = "UTI\n16 episodes")) +
  labs(title = "Episode-level distributions", subtitle = "Points are episodes; diamonds are medians", x = NULL, y = "Count") +
  theme_ruti_publication(9) + theme(legend.position = "none", plot.margin = margin(8, 8, 8, 14))

collapsed <- read_current(paths[["score_collapsed"]]) %>%
  mutate(Participant_id = as.character(.data$Participant_id), status_display = status_label(.data$UTI_Status),
         status_display = factor(.data$status_display, levels = c("Not UTI", "UTI")),
         score_label = factor(short_score[.data$score], levels = unname(short_score))) %>%
  filter(.data$score %in% names(short_score)) %>%
  group_by(.data$Participant_id, .data$score) %>%
  filter(n_distinct(.data$status_display) == 2L) %>% ungroup()
n_paired <- n_distinct(collapsed$Participant_id)
p04b <- ggplot(collapsed, aes(.data$status_display, .data$participant_status_median, group = interaction(.data$Participant_id, .data$score))) +
  geom_line(colour = "grey65", linewidth = .45, alpha = .75) +
  geom_point(aes(colour = .data$status_display), size = 1.8) +
  facet_wrap(~score_label, scales = "free_y", nrow = 2) +
  status_colour(name = "Operational UTI status") +
  labs(title = "Within-participant summaries", subtitle = sprintf("Participant-status medians for %d participants observed in both states", n_paired),
       x = NULL, y = "Median count") + theme_ruti_publication(9) +
  theme(legend.position = "none", plot.margin = margin(8, 8, 8, 14))
fig04 <- (p04a | p04b) + plot_annotation(
  title = "Virulence-factor burden and prespecified score distributions",
  subtitle = "Raw episode distributions and participant-collapsed trajectories show both imbalance and within-participant structure.",
  caption = "Distributional displays are descriptive. Participant lines include only people sampled in both operational states.",
  tag_levels = "A", theme = panel_annotation_theme
)
cap04 <- paste0(
  "Virulence-factor burden and prespecified score distributions. Panel A shows all 532 episodes (516 Not UTI; 16 UTI) for all detected VF genes, curated VF genes, UPEC-associated systems, and ExPEC-like markers. Violin envelopes show density, boxes show the median and interquartile range, points are episodes, and white diamonds mark medians. ",
  "Panel B collapses repeated episodes to participant-status medians and connects the ", n_paired, " participants observed in both operational states; line segments compare observed state-specific summaries and do not imply continuous observation. No p-values are shown. Episode-level observations are repeated within participants and should not be interpreted as independent."
)
register_figure("Fig04_vf_burden", "Fig04", "Virulence-factor burden and prespecified score distributions",
                "How do VF burdens and prespecified scores vary between operational states and within participants?", "Main", fig04, 14, 8.2, cap04,
                paths[c("score_table", "score_collapsed", "cohort")], "Descriptive episode distributions and participant-status medians",
                "Episode points are repeated within participants; paired panels include only participants observed in both states.",
                "Episode and participant-status summary", "Selected cohort; four prespecified VF scores",
                "Violins show density, boxes show IQR/median, points are episodes, and lines join participant-status medians.", "None")

# ----------------------------------------------------------------------------
# Fig05. Prespecified score-level model evidence
# ----------------------------------------------------------------------------
score_models <- read_current(paths[["score_models"]]) %>%
  filter(.data$model_variant == "batch_timepoint_collapsed_st_glmm", .data$score %in% names(short_score)) %>%
  mutate(
    score_label = factor(short_score[.data$score], levels = rev(unname(short_score))),
    estimability = if_else(is.finite(.data$OR_per_1sd) & is.finite(.data$OR_lower) & is.finite(.data$OR_upper) &
                             .data$OR_per_1sd > 0 & .data$OR_lower > 0 & .data$OR_upper > 0, "Finite estimate", "Not estimable"),
    singular = str_detect(.data$model_type, fixed("Singular")),
    q_label = sprintf("BH q = %.3f", .data$q_value_BH),
    model_state = factor(
      case_when(
        .data$estimability == "Not estimable" ~ "Not estimable",
        .data$singular ~ "Finite estimate; singular fit",
        TRUE ~ "Finite estimate; non-singular fit"
      ),
      levels = c("Finite estimate; non-singular fit", "Finite estimate; singular fit", "Not estimable")
    )
  )
assert_value(nrow(score_models), 4L, "Prespecified ST-adjusted score models")
score_state_cols <- c(
  "Finite estimate; non-singular fit" = "#0072B2",
  "Finite estimate; singular fit" = "#D55E00",
  "Not estimable" = "#6B7280"
)
score_state_shapes <- c(
  "Finite estimate; non-singular fit" = 16,
  "Finite estimate; singular fit" = 17,
  "Not estimable" = 4
)
p05 <- ggplot(score_models, aes(.data$OR_per_1sd, .data$score_label)) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey35") +
  geom_errorbarh(aes(xmin = .data$OR_lower, xmax = .data$OR_upper, colour = .data$model_state), height = .14, linewidth = .8) +
  geom_point(aes(shape = .data$model_state, colour = .data$model_state), size = 3.3, stroke = .8) +
  geom_text(
    data = score_models,
    aes(x = 1.18, y = .data$score_label, label = .data$q_label),
    inherit.aes = FALSE, hjust = 0, size = 3.25, fontface = "bold"
  ) +
  scale_x_log10(limits = c(min(score_models$OR_lower) * .75, max(score_models$OR_upper) * 2.3), breaks = c(.05, .1, .25, .5, 1, 2, 4)) +
  scale_colour_manual(values = score_state_cols, name = "Estimate and fit diagnostic") +
  scale_shape_manual(values = score_state_shapes, name = "Estimate and fit diagnostic") +
  labs(
    title = "Prespecified score-level associations with operational UTI",
    subtitle = "Adjusted odds ratio per 1-SD increase; triangles flag the four singular fits and text gives BH-adjusted q-values",
    x = "Adjusted odds ratio for UTI versus Not UTI (log scale)", y = NULL,
    caption = "Outcome: operational UTI. Covariates: batch, timepoint and collapsed ST group; random intercept: participant. BH q-values span four prespecified ST-adjusted score models."
  ) + theme_ruti_publication(10) + theme(legend.position = "bottom", plot.margin = margin(8, 80, 8, 8)) + coord_cartesian(clip = "off")
cap05 <- paste0(
  "Prespecified score-level mixed-effects model estimates for 532 episodes from 161 participants (16 UTI; 516 Not UTI). Triangles are adjusted odds ratios for operational UTI per one-standard-deviation increase in each score, horizontal lines are Wald 95% confidence intervals on a logarithmic scale, and the dashed line marks OR=1. Text reports the Benjamini-Hochberg-adjusted q-value for each score. ",
  "Models include batch, timepoint, and collapsed sequence-type group as fixed effects and a participant random intercept, with Not UTI as the reference. Benjamini-Hochberg q-values were calculated across the four prespecified ST-adjusted score models. All four fits were singular and are marked accordingly, so the apparent inverse associations are exploratory and model-dependent rather than confirmatory evidence."
)
register_figure("Fig05_vf_association_evidence", "Fig05", "VF association evidence and uncertainty",
                "What do participant-aware prespecified score models estimate, and how uncertain are they?", "Main", p05, 10.5, 5.8, cap05,
                paths[c("score_models", "cohort")], "Binomial GLMM with participant random intercept; batch, timepoint and ST-group covariates",
                "All prespecified ST-adjusted fits are singular; sparse UTI data limit inference.",
                "Episode with participant random intercept", "Four prespecified standardized VF scores; selected cohort",
                "Triangles are finite adjusted ORs from singular fits; lines are 95% CIs; text gives BH q-values.", "Benjamini-Hochberg across four ST-adjusted score models")

# ----------------------------------------------------------------------------
# Fig06. Deidentified longitudinal trajectories
# ----------------------------------------------------------------------------
timelines <- read_current(paths[["timelines"]]) %>%
  mutate(Participant_id = as.character(.data$Participant_id), tp_lab = normalise_timepoint_preserve_events(.data$tp_lab)) %>%
  semi_join(case_key, by = "Participant_id") %>%
  left_join(case_key, by = "Participant_id", relationship = "many-to-one") %>%
  mutate(status_display = status_label(.data$UTI_Status), event_display = factor(reader_event(.data$Event_type),
         levels = c("Routine surveillance", "UTI-related sampling", "Other sampling")),
         Case_Label = factor(.data$Case_Label, levels = rev(case_key$Case_Label)))
assert_unique(timelines, c("Participant_id", "tp_lab"), "Switch-case timeline episodes")
transition_segments <- casebook %>%
  select(Participant_id, from_tp, to_tp, SNPs, Case_Label) %>%
  left_join(timelines %>% select(Participant_id, from_tp = tp_lab, x_from = Time_Order),
            by = c("Participant_id", "from_tp"), relationship = "one-to-one") %>%
  left_join(timelines %>% select(Participant_id, to_tp = tp_lab, x_to = Time_Order),
            by = c("Participant_id", "to_tp"), relationship = "one-to-one") %>%
  mutate(snp_support = factor(if_else(.data$SNPs <= SAME_STRAIN_SNP_THRESHOLD, "<=25 SNPs", ">25 SNPs"), levels = c("<=25 SNPs", ">25 SNPs")),
         Case_Label = factor(.data$Case_Label, levels = levels(timelines$Case_Label)))
if (anyNA(transition_segments$x_from) || anyNA(transition_segments$x_to)) stop("A casebook transition did not map to its timeline endpoints.")
p06 <- ggplot(timelines, aes(.data$Time_Order, .data$Case_Label)) +
  geom_line(aes(group = .data$Case_Label), colour = "grey78", linetype = "dotted", linewidth = .55) +
  geom_segment(data = transition_segments,
               aes(x = .data$x_from, xend = .data$x_to, y = .data$Case_Label, yend = .data$Case_Label, linetype = .data$snp_support),
               colour = "black", linewidth = 1.25, inherit.aes = FALSE) +
  geom_point(aes(fill = .data$status_display, shape = .data$event_display, size = .data$status_display), colour = "black", stroke = .45) +
  status_fill(name = "Operational UTI status", breaks = c("UTI", "Not UTI"), labels = c("UTI", "Not UTI")) +
  scale_size_manual(values = c("UTI" = 4.5, "Not UTI" = 2.8), guide = "none") +
  scale_shape_manual(values = c("Routine surveillance" = 21, "UTI-related sampling" = 24, "Other sampling" = 22), name = "Sampling event type") +
  scale_linetype_manual(values = c("<=25 SNPs" = "solid", ">25 SNPs" = "longdash"), name = "Direct genomic evidence") +
  labs(
    title = "Deidentified trajectories for participants with a Not UTI-to-UTI transition",
    subtitle = "All observed episodes are shown; UTI symbols are larger, and bold segments distinguish five intervals with <=25-SNP support",
    x = "Days since first selected episode for that participant", y = NULL,
    caption = "Dotted lines connect observed sampling events for orientation only; they do not imply continuous carriage or observation between samples."
  ) + theme_ruti_publication(10) +
  guides(fill = guide_legend(override.aes = list(shape = c(21, 21), colour = "black", size = c(4.5, 2.8)))) +
  theme(legend.position = "bottom")
cap06 <- paste0(
  "Deidentified longitudinal trajectories for the nine participants with an adjacent Not UTI-to-UTI transition. Points are selected episodes ordered by observed collection date; fill and symbol size redundantly denote operational UTI status, and shape denotes sampling event type. Dotted grey lines connect successive observed samples for orientation only. Bold segments identify the nine Not UTI-to-UTI intervals; solid segments mark the five intervals with direct same-strain support at the operational <=25-SNP reference, whereas dashed segments exceed 25 SNPs. ",
  "Case labels are stable research-facing pseudonyms. Unequal sampling density and gaps mean that connecting lines do not imply continuous observation, persistence, or causation."
)
register_figure("Fig06_longitudinal_trajectories", "Fig06", "Deidentified longitudinal trajectories",
                "How do observed clinical states evolve around Not UTI-to-UTI transitions?", "Main", p06, 12.5, 7.2, cap06,
                paths[c("timelines", "casebook", "transitions")], "Descriptive longitudinal reconstruction with direct SNP context",
                "Connections link observed samples only; five of nine transitions meet the operational <=25-SNP reference.",
                "Episode nested within participant", "Nine data-derived adjacent Not UTI-to-UTI transitions",
                "Point fill and size are status, point shape is sampling event, and bold interval linetype is direct SNP context.", "None")

# ----------------------------------------------------------------------------
# Fig07. Within-host genomic continuity and VF similarity
# ----------------------------------------------------------------------------
rq06 <- read_current(paths[["rq06_pairs"]]) %>%
  mutate(
    transition_display = case_when(
      .data$status_from == "Not_UTI" & .data$status_to == "Not_UTI" ~ "Not UTI -> Not UTI",
      .data$status_from == "Not_UTI" & .data$status_to == "UTI" ~ "Not UTI -> UTI",
      .data$status_from == "UTI" & .data$status_to == "Not_UTI" ~ "UTI -> Not UTI",
      TRUE ~ "UTI -> UTI"
    ),
    snp_context = factor(if_else(.data$TotalSNPs <= SAME_STRAIN_SNP_THRESHOLD, "<=25 SNPs", ">25 SNPs"), levels = c("<=25 SNPs", ">25 SNPs"))
  )
assert_value(nrow(rq06), 371L, "RQ06 adjacent-pair figure denominator")
transition_cols <- c("Not UTI -> Not UTI" = "#0072B2", "Not UTI -> UTI" = "#D55E00", "UTI -> Not UTI" = "#CC79A7", "UTI -> UTI" = "#E69F00")
p07a <- ggplot(rq06, aes(.data$days_between, .data$TotalSNPs + 1, colour = .data$transition_display)) +
  geom_hline(yintercept = SAME_STRAIN_SNP_THRESHOLD + 1, linetype = "dashed", colour = "grey25") +
  geom_point(aes(shape = .data$transition_display), alpha = .55, size = 1.7) +
  scale_y_log10(labels = label_number()) +
  scale_colour_manual(values = transition_cols, name = "Operational-status transition") +
  scale_shape_manual(values = c("Not UTI -> Not UTI" = 16, "Not UTI -> UTI" = 17,
                                "UTI -> Not UTI" = 15, "UTI -> UTI" = 3),
                     name = "Operational-status transition") +
  labs(title = "Adjacent-pair genomic distance", subtitle = "The y-axis displays SNP distance + 1; the dashed line is the operational 25-SNP reference",
       x = "Days between adjacent samples", y = "SNP distance + 1") +
  theme_ruti_publication(9) + theme(plot.margin = margin(8, 8, 8, 14))
p07b <- ggplot(rq06, aes(.data$snp_context, .data$vf_jaccard, fill = .data$snp_context)) +
  geom_violin(trim = FALSE, alpha = .25, linewidth = .4) +
  geom_boxplot(width = .24, outlier.shape = NA, alpha = .55) +
  geom_point(position = position_jitter(width = .13, seed = 20260714), alpha = .32, size = .9) +
  scale_fill_manual(values = c("<=25 SNPs" = "#009E73", ">25 SNPs" = "#7A5195"), guide = "none") +
  scale_y_continuous(labels = label_number(accuracy = .1)) +
  labs(title = "VF-profile similarity", subtitle = sprintf("%d pairs <=25 SNPs; %d pairs >25 SNPs", sum(rq06$close25), sum(!rq06$close25)),
       x = "Direct genomic-distance context", y = "VF Jaccard similarity") +
  coord_cartesian(ylim = c(0, 1)) + theme_ruti_publication(9) +
  theme(plot.margin = margin(8, 8, 8, 14))
fig07 <- (p07a | p07b) + plot_layout(widths = c(.58, .42), guides = "collect") +
  plot_annotation(
    title = "Within-host genomic continuity and virulence-factor similarity",
    subtitle = "All 371 adjacent direct comparisons from 139 participants are shown.",
    caption = "The <=25-SNP reference is operational rather than proof of persistence or transmission; participants contribute unequal numbers of pairs.",
    tag_levels = "A", theme = panel_annotation_theme
  ) & theme(legend.position = "bottom")
cap07 <- paste0(
  "Within-host genomic continuity and virulence-factor similarity for 371 adjacent direct comparisons from 139 participants. Panel A plots days between samples against pairwise SNP distance plus one on a logarithmic scale; colour and shape redundantly indicate the observed operational-status transition and the dashed line corresponds to 25 SNPs after the displayed +1 transformation. Panel B shows VF Jaccard similarity for 140 pairs at or below and 231 pairs above the operational 25-SNP reference; violins show density, boxes show medians and interquartile ranges, and points are comparisons. ",
  "The threshold is an operational genomic-continuity reference, not evidence of transmission. Participants contribute unequal numbers of comparisons, so points are not independent."
)
register_figure("Fig07_within_host_genomic_continuity", "Fig07", "Within-host genomic continuity and VF similarity",
                "How do adjacent direct SNP distances relate to time, clinical transition, and VF similarity?", "Main", fig07, 14, 6.5, cap07,
                paths[c("rq06_pairs", "transitions")], "Descriptive adjacent-pair comparison; operational 25-SNP reference",
                "Comparisons repeat within participants; genomic proximity is not demonstrated transmission.",
                "Adjacent pair nested within participant", "371 adjacent direct comparisons from 139 participants",
                "Points are pairs; colour is status transition; box/violin summaries stratify by SNP reference.", "None")

# ----------------------------------------------------------------------------
# Fig08. Reference-aware variant map (strictly conditional)
# ----------------------------------------------------------------------------
variant_valid <- FALSE
if (file.exists(paths[["variants"]]) && file.exists(paths[["variant_validation"]])) {
  vv <- read_current(paths[["variant_validation"]])
  variant_valid <- nrow(vv) > 0 && "Figure_Eligible" %in% names(vv) && all(vv$Figure_Eligible %in% TRUE)
}
if (variant_valid) {
  var <- read_current(paths[["variants"]])
  required_var_cols <- c("Reference_ID", "Participant_id", "From_Time", "To_Time",
                         "Reference_Cumulative_Position", "Reference_Total_Length",
                         "Reference_Contig_ID", "Reference_Contig_Offset",
                         "Reference_Contig_Length", "Gene", "Type", "Figure_Eligible")
  if (any(!required_var_cols %in% names(var))) stop("Reference-aware variant table lacks required plotting columns.")
  reference_key <- var %>% distinct(Reference_ID) %>% arrange(.data$Reference_ID) %>%
    mutate(reference_pseudonym = sprintf("Reference %02d", row_number()))
  var_case <- var %>%
    filter(.data$Figure_Eligible %in% TRUE) %>%
    mutate(Participant_id = as.character(.data$Participant_id),
           From_Time = normalise_timepoint_preserve_events(.data$From_Time),
           To_Time = normalise_timepoint_preserve_events(.data$To_Time)) %>%
    left_join(casebook %>% select(case_id, Participant_id, from_tp, to_tp, Case_Label),
              by = c("Participant_id", "From_Time" = "from_tp", "To_Time" = "to_tp"), relationship = "many-to-one") %>%
    left_join(reference_key, by = "Reference_ID", relationship = "many-to-one") %>%
    mutate(Case_Label = factor(.data$Case_Label, levels = rev(case_key$Case_Label)),
           reference_label = paste(.data$reference_pseudonym, .data$Case_Label, sep = " - "),
           reference_contig = .data$Reference_Contig_ID,
           cumulative_position = .data$Reference_Cumulative_Position,
           reference_length = .data$Reference_Total_Length,
           contig_start = .data$Reference_Contig_Offset,
           contig_end = .data$Reference_Contig_Offset + .data$Reference_Contig_Length,
           label = if_else(!is.na(.data$Gene) & nzchar(.data$Gene), .data$Gene, NA_character_),
           mutation_type = if_else(.data$Type %in% c("SNP", "Substitution"), "Substitution", as.character(.data$Type)))
  if (anyNA(var_case$Case_Label)) stop("A validated variant lacks a deidentified case label.")
  if (!identical(sort(unique(var_case$mutation_type)), "Substitution")) {
    stop("The current reference-aware variant map expects validated substitutions only.")
  }
  lab_rule <- var_case %>% filter(!is.na(.data$label)) %>% group_by(.data$case_id, .data$label) %>% slice_min(.data$cumulative_position, n = 1, with_ties = FALSE) %>% ungroup()
  contigs <- var_case %>% distinct(reference_label, reference_contig, contig_start, contig_end)
  p08 <- ggplot(var_case, aes(.data$cumulative_position / 1e6, 1)) +
    geom_vline(data = contigs, aes(xintercept = .data$contig_end / 1e6), colour = "grey88", linewidth = .25) +
    geom_hline(yintercept = 1, colour = "grey70", linewidth = .35) +
    geom_point(size = 2, alpha = .85, colour = "#4C78A8") +
    ggrepel::geom_text_repel(data = lab_rule, aes(label = .data$label), seed = 20260714, size = 2.5,
                             min.segment.length = 0, max.overlaps = Inf, direction = "both", show.legend = FALSE) +
    facet_wrap(~reference_label, scales = "free_x", ncol = 1) +
    scale_y_continuous(NULL, breaks = NULL, expand = expansion(add = c(.55, .8))) +
    labs(title = "Reference-aware within-host variant map",
         subtitle = sprintf("All %d validated variants are substitutions; coordinates are cumulative only within each exact reference", nrow(var_case)),
         x = "Reference-genome coordinate (Mb)", y = NULL,
         caption = "Vertical lines mark contig boundaries. Gene labels follow a reproducible rule: the first mapped occurrence of each annotated gene within each comparison. Reference panels are not coordinate-comparable.") +
    theme_ruti_publication(9.5) + theme(panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank())
  cap08 <- paste0(
    "Reference-aware within-host variant map for ", nrow(var_case), " validated variants from ",
    n_distinct(var_case$case_id), " Not UTI-to-UTI comparisons across ",
    n_distinct(var_case$reference_pseudonym), " exact references. All validated variants are substitutions. Points are positioned using contig lengths derived from each exact reference FASTA; cumulative coordinates are calculated only within a reference, and different references are faceted separately. Vertical grey lines mark contig boundaries, and labels identify the first mapped occurrence of each annotated gene per comparison. ",
    "Every position was checked against its named reference contig and reference hash before plotting. Reference panels are not directly coordinate-comparable, unlabelled points remain included, and candidate variants are descriptive rather than causal."
  )
  register_figure("Fig08_reference_aware_variant_map", "Fig08", "Reference-aware within-host variant map",
                  "Where do validated within-host variants fall on their exact reference assemblies?", "Main", p08, 12, 8.2, cap08,
                  paths[c("variants", "variant_validation", "casebook")], "Validated reference-contig coordinate reconstruction",
                  "Coordinates are comparable only within the same reference; candidate variants are not causal evidence.",
                  "Variant within reference-specific comparison", "Only variants passing exact reference/contig/position validation",
                  "Points are substitutions, facets are exact references, vertical lines are contig boundaries, labels follow a deterministic gene rule.", "None")
} else {
  writeLines(c(
    "Fig08 was not generated.",
    "Reason: the exact reference/contig coordinate validator did not report all checks as passing.",
    "The obsolete position maps from scripts 21 and the previous script 35 must not be used."
  ), file.path(DIR_FIG_AUDIT, "Fig08_UNVALIDATED.txt"))
}

# ----------------------------------------------------------------------------
# FigS01. Curated VF presence/absence/missing heatmap
# ----------------------------------------------------------------------------
vf_pa <- read_current(paths[["vf_pa"]]) %>% mutate(Participant_id = as.character(.data$Participant_id), tp_lab = normalise_timepoint_preserve_events(.data$tp_lab))
assert_unique(vf_pa, c("Participant_id", "tp_lab"), "VF presence/absence matrix")
gene_cols <- setdiff(names(vf_pa), c("Participant_id", "tp_lab", "Episode_ID"))
gene_prev <- tibble(gene = gene_cols, prevalence = vapply(vf_pa[gene_cols], function(z) mean(z == 1, na.rm = TRUE), numeric(1)),
                    missing_fraction = vapply(vf_pa[gene_cols], function(z) mean(is.na(z)), numeric(1)))
display_genes <- gene_prev %>% filter(.data$prevalence >= .05, .data$prevalence <= .95) %>% pull("gene")
gene_map <- read_current(paths[["gene_map"]]) %>% transmute(gene = .data$Gene, category = .data$Category)
gene_order <- gene_prev %>% filter(.data$gene %in% display_genes) %>%
  left_join(gene_map, by = "gene", relationship = "many-to-one") %>%
  mutate(category = coalesce(.data$category, "Unassigned"), category = str_replace_all(.data$category, "/", " / ")) %>%
  arrange(.data$category, desc(.data$prevalence)) %>% mutate(gene_factor = factor(.data$gene, levels = rev(.data$gene)))
episode_order <- cohort %>% select(Participant_id, tp_lab, status_display, Collection_Date) %>%
  left_join(mlst %>% select(Participant_id, tp_lab, ST), by = c("Participant_id", "tp_lab"), relationship = "one-to-one") %>%
  arrange(.data$status_display, .data$Participant_id, .data$Collection_Date, .data$ST, .data$tp_lab) %>% mutate(episode_index = row_number())
heat_long <- vf_pa %>% inner_join(episode_order %>% select(Participant_id, tp_lab, episode_index),
                                  by = c("Participant_id", "tp_lab"), relationship = "one-to-one") %>%
  pivot_longer(all_of(display_genes), names_to = "gene", values_to = "presence") %>%
  left_join(gene_order %>% select(gene, category, gene_factor), by = "gene", relationship = "many-to-one") %>%
  mutate(state = factor(case_when(is.na(.data$presence) ~ "Unavailable", .data$presence == 1 ~ "Present", TRUE ~ "Absent"),
                        levels = c("Present", "Absent", "Unavailable")))
pS01bar <- ggplot(episode_order, aes(.data$episode_index, 1, fill = .data$status_display)) + geom_raster() + status_fill(name = "Operational UTI status") +
  labs(x = NULL, y = NULL) + theme_void() + theme(legend.position = "bottom")
pS01heat <- ggplot(heat_long, aes(.data$episode_index, .data$gene_factor, fill = .data$state)) +
  geom_raster() + facet_grid(rows = vars(category), scales = "free_y", space = "free_y", switch = "y") +
  scale_fill_manual(values = c(Present = "#253494", Absent = "#F7FBFF", Unavailable = "#D55E00"), name = "VF state", drop = FALSE) +
  labs(x = "Episodes ordered by operational status, participant, collection date and sequence type", y = "Virulence-factor gene") +
  theme_ruti_publication(7) + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), strip.placement = "outside",
                                      strip.text.y.left = element_text(angle = 0, size = 6.3), legend.position = "bottom")
figS01 <- pS01bar / pS01heat + plot_layout(heights = c(.035, .965), guides = "collect") +
  plot_annotation(title = "Curated virulence-factor presence heatmap",
                  subtitle = sprintf("%d genes with 5-95%% prevalence across 532 episodes; no binary Euclidean clustering", length(display_genes)),
                  caption = "Rows are biologically ordered by functional category and prevalence. Missing/unavailable calls have a distinct state and are never converted to absence.",
                  theme = panel_annotation_theme) & theme(legend.position = "bottom")
capS01 <- paste0(
  "Curated VF heatmap for 532 selected episodes and ", length(display_genes), " genes with overall prevalence between 5% and 95%. Columns are episodes ordered by operational UTI status, participant, collection date, and sequence type; rows are genes ordered by curated functional category and decreasing prevalence. Dark blue indicates presence, pale blue absence, and orange unavailable data. The status strip uses the canonical operational colours. No clustering or inferential test is applied, avoiding inappropriate default Euclidean clustering of binary observations. Repeated participant episodes can visually dominate recurrent profiles and should not be treated as independent clusters."
)
register_figure("FigS01_vf_presence_heatmap", "FigS01", "Curated VF presence heatmap",
                "Which variably prevalent curated VF genes co-occur across episodes?", "Supplementary", figS01, 15, 15.5, capS01,
                paths[c("vf_pa", "gene_map", "cohort", "mlst")], "Descriptive binary matrix",
                "Repeated episodes can dominate patterns; ordering is biological rather than a transmission or lineage inference.",
                "Gene-by-episode cell", "Genes with 5-95% prevalence; all selected episodes",
                "Dark blue present, pale blue absent, orange unavailable; status strip uses canonical colours.", "None")

# ----------------------------------------------------------------------------
# FigS02. Exact 532-tip core-genome tree
# ----------------------------------------------------------------------------
tips_in_component <- function(tree, start, blocked) {
  n_tip <- length(tree$tip.label)
  edges <- tree$edge
  adj <- split(c(edges[, 2], edges[, 1]), c(edges[, 1], edges[, 2]))
  seen <- blocked
  queue <- start
  out <- integer()
  while (length(queue)) {
    node <- queue[[1]]; queue <- queue[-1]
    if (node %in% seen) next
    seen <- c(seen, node)
    if (node <= n_tip) out <- c(out, node)
    queue <- c(queue, setdiff(adj[[as.character(node)]], seen))
  }
  unique(out)
}
midpoint_root_display <- function(tree) {
  if (is.null(tree$edge.length)) stop("Core tree lacks branch lengths.")
  negative_n <- sum(tree$edge.length < 0, na.rm = TRUE)
  tree$edge.length[tree$edge.length < 0] <- 0
  dm <- cophenetic.phylo(tree)
  ij <- which(dm == max(dm), arr.ind = TRUE)[1, ]
  tip1 <- ij[[1]]; tip2 <- ij[[2]]
  edges <- tree$edge
  lens <- tree$edge.length
  adj <- split(c(edges[, 2], edges[, 1]), c(edges[, 1], edges[, 2]))
  edge_key <- function(a, b) which((edges[,1] == a & edges[,2] == b) | (edges[,1] == b & edges[,2] == a))[[1]]
  prev <- rep(NA_integer_, length(tree$tip.label) + tree$Nnode)
  queue <- tip1; prev[tip1] <- 0L
  while (length(queue) && is.na(prev[tip2])) {
    node <- queue[[1]]; queue <- queue[-1]
    for (nb in adj[[as.character(node)]]) if (is.na(prev[nb])) { prev[nb] <- node; queue <- c(queue, nb) }
  }
  path <- tip2
  while (tail(path, 1) != tip1) path <- c(path, prev[tail(path, 1)])
  path <- rev(path)
  path_lens <- vapply(seq_len(length(path) - 1L), function(i) lens[edge_key(path[i], path[i+1])], numeric(1))
  target <- sum(path_lens) / 2
  j <- which(cumsum(path_lens) >= target)[1]
  before <- if (j == 1L) 0 else sum(path_lens[seq_len(j - 1L)])
  left_len <- target - before
  right_len <- path_lens[j] - left_len
  left_node <- path[j]; right_node <- path[j + 1L]
  left_tips <- tips_in_component(tree, left_node, right_node)
  rooted <- ape::root(tree, outgroup = left_tips, resolve.root = TRUE)
  root_rows <- which(rooted$edge[, 1] == length(rooted$tip.label) + 1L)
  if (length(root_rows) != 2L) stop("Midpoint display root did not create two root branches.")
  descendant_has_tip <- function(child, tip) tip %in% tips_in_component(rooted, child, length(rooted$tip.label) + 1L)
  for (rr in root_rows) rooted$edge.length[rr] <- if (descendant_has_tip(rooted$edge[rr, 2], tip1)) left_len else right_len
  attr(rooted, "negative_edges_truncated_for_display") <- negative_n
  rooted
}
tree <- read.tree(paths[["tree"]])
tree_map <- read_current(paths[["tree_map"]]) %>% mutate(Participant_id = as.character(.data$Participant_id), tp_lab = normalise_timepoint_preserve_events(.data$tp_lab))
assert_value(length(tree$tip.label), 532L, "Core-genome tree tips")
assert_unique(tree_map, c("parsnp_alignment_label"), "Core-genome tree sample map")
if (!setequal(tree$tip.label, tree_map$parsnp_alignment_label)) stop("Core tree tips do not exactly match the 532-row sample map.")
transition_tip_labels <- casebook %>%
  transmute(.data$Participant_id, tp_lab = .data$to_tp, .data$Case_Label) %>%
  distinct()
assert_unique(transition_tip_labels, c("Participant_id", "tp_lab"), "Transition endpoint tip labels")
assert_value(nrow(transition_tip_labels), 9L, "Transition endpoint tip-label count")
tree_meta <- tree_map %>% transmute(label = .data$parsnp_alignment_label, .data$Participant_id, .data$tp_lab) %>%
  left_join(cohort %>% select(Participant_id, tp_lab, status_display), by = c("Participant_id", "tp_lab"), relationship = "one-to-one") %>%
  left_join(mlst %>% select(Participant_id, tp_lab, ST), by = c("Participant_id", "tp_lab"), relationship = "one-to-one") %>%
  left_join(transition_tip_labels, by = c("Participant_id", "tp_lab"), relationship = "one-to-one") %>%
  mutate(ST_plot = case_when(is.na(.data$ST) ~ "Missing ST", as.character(.data$ST) == "131" ~ "ST131", as.character(.data$ST) == "73" ~ "ST73",
                             as.character(.data$ST) == "69" ~ "ST69", TRUE ~ "Other typed ST"),
         case_tip_label = if_else(!is.na(.data$Case_Label), .data$Case_Label, NA_character_))
if (anyNA(tree_meta$status_display)) stop("Core tree metadata join left an unmatched tip.")
tree_display <- midpoint_root_display(tree)
negative_edges <- attr(tree_display, "negative_edges_truncated_for_display")
tree_base <- ggtree(tree_display, linewidth = .22) %<+% tree_meta
tree_label_data <- tree_base$data %>% filter(.data$isTip, !is.na(.data$case_tip_label))
assert_value(nrow(tree_label_data), 9L, "Labelled transition endpoints on the core tree")
pS02 <- tree_base +
  geom_tippoint(aes(colour = .data$status_display, shape = .data$ST_plot), size = .75, alpha = .8) +
  ggrepel::geom_text_repel(data = tree_label_data, aes(x = .data$x, y = .data$y, label = .data$case_tip_label),
                           inherit.aes = FALSE, direction = "y", hjust = 0,
                           nudge_x = max(tree_base$data$x) * .015, seed = 20260714,
                           size = 1.8, min.segment.length = 0, max.overlaps = Inf, colour = "grey20") +
  status_colour(name = "Operational UTI status") +
  scale_shape_manual(values = c(ST131 = 16, ST73 = 17, ST69 = 15, `Other typed ST` = 1, `Missing ST` = 4), name = "Sequence-type annotation") +
  geom_treescale(width = 100, x = 0, y = 0, fontsize = 3) +
  labs(title = "Canonical core-genome neighbour-joining tree",
       subtitle = "Exact 532-tip sample map; midpoint-rooted for display only; branch lengths are SNP-distance units",
       caption = sprintf("No branch-support values are claimed. %d negative neighbour-joining branch estimates were truncated to zero for display; the source tree is unchanged.", negative_edges)) +
  theme_tree2() + theme_ruti_publication(8) + theme(axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(), legend.position = "bottom")
capS02 <- paste0(
  "Canonical core-genome neighbour-joining tree for all 532 selected assemblies. Every tip was matched one-to-one to the exact core-SNP sample map before clinical and sequence-type metadata were added. Branch lengths are in SNP-distance units; the tree was midpoint-rooted only for display and no rooting or transmission claim is made. Tip colour denotes operational UTI status, tip shape gives a compact sequence-type annotation, and labels identify only the nine deidentified UTI endpoints of the data-derived Not UTI-to-UTI transitions. ",
  negative_edges, " negative neighbour-joining branch estimates were truncated to zero for display while the source tree remained unchanged. Branch support values were unavailable and are not shown. Repeated episodes from the same participant may cluster and are not independent."
)
register_figure("FigS02_core_genome_phylogeny", "FigS02", "Canonical core-genome phylogeny",
                "How are selected assemblies related in the core-genome SNP tree?", "Supplementary", pS02, 13, 14, capS02,
                paths[c("tree", "tree_map", "cohort", "mlst")], "Neighbour-joining tree from the canonical core-SNP distance matrix",
                "Midpoint root is display-only; genomic clustering is not demonstrated transmission; no support values available.",
                "Assembly tip", "Exact 532-tip core-SNP sample map",
                "Branch length is SNP distance; colour is status; shape is compact ST annotation; labels are deidentified cases.", "None")

# ----------------------------------------------------------------------------
# FigS03. Module gain/loss
# ----------------------------------------------------------------------------
module <- read_current(paths[["module_changes"]]) %>% mutate(Participant_id = as.character(.data$Participant_id)) %>%
  select(-any_of("case_label")) %>% left_join(case_key, by = c("case_id", "Participant_id"), relationship = "many-to-one") %>%
  mutate(Case_Label = factor(.data$Case_Label, levels = case_key$Case_Label),
         module_label = factor(.data$module_label, levels = rev(unique(.data$module_label))),
         change_display = factor(.data$change_display, levels = c("Gained", "Lost", "Stable present", "Stable absent")))
pS03 <- ggplot(module, aes(.data$Case_Label, .data$module_label, fill = .data$change_display)) +
  geom_tile(colour = "white", linewidth = .25) +
  scale_fill_manual(values = c("Gained" = "#009E73", "Lost" = "#CC79A7", "Stable present" = "#0072B2", "Stable absent" = "#F2F2F2"), name = "Module state") +
  labs(title = "VF module changes across Not UTI-to-UTI transitions", subtitle = "Nine deidentified adjacent comparisons; stable absence remains explicit",
       x = "Transition case", y = "VF module", caption = "Presence/absence changes describe observed endpoints and do not establish acquisition timing or causation.") +
  theme_ruti_publication(8.5) + theme(axis.text.x = element_text(angle = 35, hjust = 1), legend.position = "bottom")
capS03 <- paste0(
  "Virulence-factor module gain/loss matrix for ", n_distinct(module$case_id),
  " deidentified adjacent Not UTI-to-UTI comparisons. Columns are transition cases and rows are ",
  n_distinct(module$module_label),
  " curated VF modules. Green indicates gained, magenta lost, blue stable present, and pale grey stable absent between the two observed endpoints; unavailable values would remain distinct rather than be converted to absence. The analysis is descriptive, repeated episodes are nested within participants, and endpoint differences do not establish when a change occurred or whether it caused UTI."
)
register_figure("FigS03_module_gain_loss", "FigS03", "VF module gain/loss",
                "Which VF modules differ between adjacent Not UTI and UTI endpoints?", "Supplementary", pS03, 12.5, 9.5, capS03,
                paths[c("module_changes", "casebook")], "Descriptive paired presence/absence comparison",
                "Endpoint changes do not determine acquisition timing or causation.", "Module within transition pair",
                "Nine adjacent Not UTI-to-UTI transitions", "Tile fill denotes gained, lost, stable present, or stable absent.", "None")

# ----------------------------------------------------------------------------
# FigS04. Global VF PCoA
# ----------------------------------------------------------------------------
pcoa <- read_current(paths[["pcoa"]]) %>% mutate(Participant_id = as.character(.data$Participant_id), tp_lab = normalise_timepoint_preserve_events(.data$tp_lab)) %>%
  select(-any_of(c("UTI_Status", "status_display"))) %>%
  left_join(cohort %>% select(Participant_id, tp_lab, status_display), by = c("Participant_id", "tp_lab"), relationship = "one-to-one") %>%
  left_join(mlst %>% select(Participant_id, tp_lab, ST_plot = ST), by = c("Participant_id", "tp_lab"), relationship = "one-to-one") %>%
  mutate(ST_group = case_when(is.na(.data$ST_plot) ~ "Missing ST", as.character(.data$ST_plot) == "131" ~ "ST131", as.character(.data$ST_plot) == "73" ~ "ST73", TRUE ~ "Other typed ST"))
axis1 <- unique(pcoa$var_Axis1)[1]; axis2 <- unique(pcoa$var_Axis2)[1]
pS04 <- ggplot(pcoa, aes(.data$Axis1, .data$Axis2, colour = .data$status_display, shape = .data$ST_group)) +
  geom_point(aes(size = .data$status_display), alpha = .72) + status_colour(name = "Operational UTI status") +
  scale_size_manual(values = c("Not UTI" = 1.8, "UTI" = 3), guide = "none") +
  scale_shape_manual(values = c(ST131 = 16, ST73 = 17, `Other typed ST` = 1, `Missing ST` = 4), name = "Sequence-type group") +
  labs(title = "Global VF-profile principal coordinates", subtitle = "Jaccard distance on module presence/absence; descriptive and unadjusted",
       x = sprintf("PCoA axis 1 (%.1f%%)", axis1), y = sprintf("PCoA axis 2 (%.1f%%)", axis2),
       caption = str_wrap("Status colour and size redundantly encode UTI status; shape encodes the ST group. Descriptive only: no independent-episode inference.", width = 115)) +
  theme_ruti_publication(10) + theme(legend.position = "bottom")
capS04 <- "Principal-coordinate analysis of Jaccard distances among VF module presence/absence profiles for all 532 selected episodes. Each point is an episode; colour denotes operational UTI status, UTI points are also larger for accessibility, and shape provides a compact sequence-type annotation. Axes report the displayed proportion of variation. The ordination is descriptive, has no hypothesis-test overlay, and is not adjusted for repeated episodes, unequal group sizes, or lineage structure. Visual overlap or separation should therefore not be interpreted as an independent or causal UTI association."
register_figure("FigS04_vf_pcoa", "FigS04", "Global VF-profile PCoA",
                "Do global VF module profiles visually separate by operational status?", "Supplementary", pS04, 8.5, 6.5, capS04,
                paths[c("pcoa", "cohort", "mlst")], "Jaccard PCoA",
                "Unadjusted ordination; points repeat within participants and reflect lineage structure.", "Episode",
                "All selected episodes and VF module presence/absence", "Point colour and size encode status; shape is compact ST group.", "None")

# ----------------------------------------------------------------------------
# FigS05. Near-miss and leave-one-UTI diagnostics
# ----------------------------------------------------------------------------
near <- read_current(paths[["near_miss"]])
near_features <- tibble(
  criterion = c("Culture supports UTI", "Local urinary symptom", "Systemic symptom", "Symptom rule met", "Indwelling-catheter rule"),
  n = c(sum(near$culture_supports_uti %in% TRUE, na.rm = TRUE), sum(near$local_urinary_symptom_any %in% TRUE, na.rm = TRUE),
        sum(near$systemic_symptom_any %in% TRUE, na.rm = TRUE), sum(near$symptom_rule_met %in% TRUE, na.rm = TRUE),
        sum(near$catheter_rule == "B_indwelling", na.rm = TRUE))
)
pS05a <- ggplot(near_features, aes(.data$n, reorder(.data$criterion, .data$n))) + geom_col(fill = "#E69F00", width = .68) +
  geom_text(aes(label = .data$n), hjust = -.25, fontface = "bold") + scale_x_continuous(expand = expansion(mult = c(0, .15))) +
  labs(title = "Near-miss rule components", subtitle = sprintf("%d episodes remain Not UTI under the primary definition", nrow(near)), x = "Episodes meeting component", y = NULL) + theme_ruti_publication(9)
loo <- read_current(paths[["leave_one"]]) %>% filter(.data$feature_type == "score", .data$feature %in% names(short_score), .data$statistic == "mean difference") %>%
  mutate(score_label = factor(short_score[.data$feature], levels = rev(unname(short_score))))
pS05b <- ggplot(loo, aes(.data$full_effect, .data$score_label)) + geom_vline(xintercept = 0, linetype = "dashed") +
  geom_segment(aes(x = .data$min_without_one_uti, xend = .data$max_without_one_uti, yend = .data$score_label, colour = .data$direction_flip), linewidth = 1.2) +
  geom_point(size = 3, colour = "black") + scale_colour_manual(values = c(`TRUE` = "#D55E00", `FALSE` = "#0072B2"), labels = c(`TRUE` = "Direction flip", `FALSE` = "No direction flip"), name = NULL) +
  labs(title = "Leave-one-UTI-out sensitivity", subtitle = "Point: full UTI-minus-Not-UTI mean difference; line: range after omitting one UTI episode",
       x = "Mean difference", y = NULL) + theme_ruti_publication(9)
figS05 <- (pS05a | pS05b) + plot_annotation(title = "Operational-definition and sparse-case diagnostics", tag_levels = "A", theme = panel_annotation_theme) & theme(legend.position = "bottom")
capS05 <- paste0(
  "Operational-definition and sparse-case diagnostics. Panel A counts the prespecified clinical components among ", nrow(near), " near-miss episodes that remain Not UTI under the primary operational definition; criteria overlap and bars must not be summed. Panel B shows full UTI-minus-Not-UTI mean differences for four prespecified VF scores and the range obtained after omitting each of the 16 UTI episodes in turn; colour flags direction changes. These are sensitivity diagnostics, not corrected association tests. Episodes repeat within participants, and leave-one-episode analysis does not replace participant-aware modelling."
)
register_figure("FigS05_near_miss_leave_one_uti", "FigS05", "Near-miss and leave-one-UTI diagnostics",
                "How sensitive are descriptive findings to the operational rule and individual UTI episodes?", "Supplementary", figS05, 13.5, 6.1, capS05,
                paths[c("near_miss", "leave_one")], "Rule-component audit and leave-one-UTI sensitivity",
                "Criteria overlap; sensitivity ranges are not inferential confidence intervals.", "Episode and score diagnostic",
                "Prespecified near-miss rows and four score mean differences", "Bars count rule components; points are full effects; lines are leave-one ranges.", "None")

# ----------------------------------------------------------------------------
# FigS06. Transition mechanism composition
# ----------------------------------------------------------------------------
mech <- read_current(paths[["mechanism_summary"]]) %>%
  mutate(transition_display = recode(.data$transition_type, `Not_UTI->Not_UTI` = "Not UTI -> Not UTI", `Not_UTI->UTI` = "Not UTI -> UTI",
                                     `UTI->Not_UTI` = "UTI -> Not UTI", `UTI->UTI` = "UTI -> UTI"),
         mechanism_display = factor(mechanism_labels[.data$mechanism_bucket], levels = mechanism_labels)) %>%
  group_by(.data$transition_display) %>% mutate(prop = .data$n_transitions / sum(.data$n_transitions)) %>% ungroup()
pS06 <- ggplot(mech, aes(.data$transition_display, .data$prop, fill = .data$mechanism_bucket)) + geom_col(colour = "white", width = .72) +
  geom_text(aes(label = if_else(.data$n_transitions > 0, as.character(.data$n_transitions), "")),
            position = position_stack(vjust = .5), size = 2.7, fontface = "bold") +
  scale_y_continuous(labels = label_percent()) + scale_fill_manual(values = mechanism_colours, labels = mechanism_labels, name = "Evidence category") +
  labs(title = "Evidence categories across adjacent operational-status transitions", subtitle = "Counts are derived from all 371 adjacent comparisons",
       x = "Observed adjacent transition", y = "Proportion of transitions",
       caption = "Categories organise SNP, ST and VF evidence; they do not prove a biological mechanism or transmission.") +
  theme_ruti_publication(9.5) + theme(legend.position = "bottom")
capS06 <- "Evidence-category composition across 371 adjacent comparisons from 139 participants: 349 Not UTI-to-Not UTI, nine Not UTI-to-UTI, 12 UTI-to-Not UTI, and one UTI-to-UTI interval. Bar height is the within-transition-type proportion, fill combines direct SNP, sequence-type, and VF-profile context, and visible labels give transition counts. Categories are descriptive evidence summaries, not proof of replacement, within-host evolution, persistence, or transmission. Repeated comparisons from the same participant are not independent."
register_figure("FigS06_transition_mechanisms", "FigS06", "Transition mechanism composition",
                "How do evidence categories vary across observed operational-status transitions?", "Supplementary", pS06, 10.5, 6.5, capS06,
                paths[c("mechanism_summary", "transitions")], "Descriptive mechanism evidence classification",
                "Evidence categories are not causal mechanisms or transmission claims.", "Adjacent pair",
                "All 371 adjacent comparisons", "Stacked proportions show evidence-category composition; labels are counts.", "None")

# Keep the historical unnumbered compatibility alias in sync with the canonical
# numbered figure because current presentation builders still resolve release
# assets through the audited plots/final_figures inventory.
transition_alias_dir <- file.path(DIR_PLOTS, "final_figures")
dir.create(transition_alias_dir, recursive = TRUE, showWarnings = FALSE)
transition_alias_ok <- vapply(c("png", "pdf"), function(ext) {
  file.copy(
    file.path(DIR_SUPP, paste0("FigS06_transition_mechanisms.", ext)),
    file.path(transition_alias_dir, paste0("transition_mechanisms_by_transition_type.", ext)),
    overwrite = TRUE
  )
}, logical(1))
if (!all(transition_alias_ok)) stop("Failed to refresh the transition-mechanism compatibility aliases.")

# ----------------------------------------------------------------------------
# FigS07. Predicted plasmid and AMR mechanism context
# ----------------------------------------------------------------------------
plasmid_manifest <- read_current(paths[["plasmid_manifest"]]) %>% mutate(Participant_id = as.character(.data$Participant_id), tp_lab = normalise_timepoint_preserve_events(.data$tp_lab))
assert_unique(plasmid_manifest, c("Participant_id", "tp_lab"), "PlasmidFinder manifest")
assert_value(nrow(plasmid_manifest), 532L, "PlasmidFinder manifest denominator")
rep <- read_current(paths[["replicons"]]) %>% distinct(isolate_id, replicon) %>%
  left_join(mlst %>% transmute(isolate_id = .data$Isolate_ID, .data$Participant_id, .data$tp_lab), by = "isolate_id", relationship = "many-to-one")
if (anyNA(rep$Participant_id)) stop("Replicon hit did not map to the selected MLST isolate table.")
rep_prev <- rep %>% count(.data$replicon, name = "n") %>% mutate(prevalence = .data$n / 532) %>% filter(.data$prevalence >= .05) %>% pull("replicon")
rep_summary <- plasmid_manifest %>% select(Participant_id, tp_lab) %>%
  left_join(cohort %>% select(Participant_id, tp_lab, status_display), by = c("Participant_id", "tp_lab"), relationship = "one-to-one") %>%
  crossing(replicon = rep_prev) %>%
  left_join(rep %>% filter(.data$replicon %in% rep_prev) %>% transmute(.data$Participant_id, .data$tp_lab, .data$replicon, present = TRUE),
            by = c("Participant_id", "tp_lab", "replicon"), relationship = "one-to-one") %>%
  mutate(present = coalesce(.data$present, FALSE)) %>% group_by(.data$replicon, .data$status_display) %>%
  summarise(n = sum(.data$present), N = n(), percent = 100 * .data$n / .data$N, .groups = "drop")
rep_order <- rep_summary %>% group_by(.data$replicon) %>% summarise(n = sum(.data$n), .groups = "drop") %>% arrange(.data$n) %>% pull("replicon")
rep_summary <- rep_summary %>%
  mutate(replicon = factor(.data$replicon, levels = rep_order),
         label_y = as.numeric(.data$replicon) + if_else(.data$status_display == "UTI", .20, -.20),
         label_hjust = if_else(.data$status_display == "UTI" & .data$percent >= 5, 1.25, -.25))
status_order <- c("UTI", "Not UTI")
rep_denoms <- rep_summary %>% distinct(.data$status_display, .data$N)
rep_n <- setNames(rep_denoms$N[match(status_order, as.character(rep_denoms$status_display))], status_order)
if (anyNA(rep_n)) stop("Replicon status denominators are incomplete.")
rep_status_labels <- setNames(sprintf("%s (n=%d)", status_order, rep_n), status_order)
pS07a <- ggplot(rep_summary, aes(.data$percent, .data$replicon, colour = .data$status_display)) +
  geom_point(aes(shape = .data$status_display), size = 3) +
  geom_text(aes(y = .data$label_y, label = paste0("n=", .data$n), hjust = .data$label_hjust), size = 3, show.legend = FALSE) +
  status_colour(name = "Operational UTI status", breaks = status_order, labels = unname(rep_status_labels)) +
  scale_shape_manual(values = c("Not UTI" = 16, "UTI" = 17), name = "Operational UTI status",
                     breaks = status_order, labels = unname(rep_status_labels)) +
  scale_x_continuous(labels = label_number(suffix = "%"), expand = expansion(mult = c(0, .18))) +
  labs(title = "Replicon prevalence", subtitle = "Replicons present in at least 5% overall; legend gives the within-status denominators",
       x = "Within-status prevalence", y = NULL) +
  theme_ruti_publication(10) + theme(axis.text.y = element_text(size = 9.2))

amr <- read_current(paths[["amr"]]) %>% mutate(Participant_id = as.character(.data$Participant_id), tp_lab = normalise_timepoint_preserve_events(.data$tp_lab))
assert_unique(amr, c("Participant_id", "tp_lab"), "Script-29 AMR profile")
assert_value(nrow(amr), 532L, "Script-29 AMR denominator")
mob_profiles <- read_current(paths[["mob_profiles"]]) %>%
  mutate(
    Participant_id = as.character(.data$Participant_id),
    tp_lab = normalise_timepoint_preserve_events(.data$tp_lab)
  )
mechanism_profiles <- read_current(paths[["plasmid_mechanism_profiles"]]) %>%
  mutate(
    Participant_id = as.character(.data$Participant_id),
    tp_lab = normalise_timepoint_preserve_events(.data$tp_lab)
  )
location_validation <- read_current(paths[["plasmid_location_validation"]])
assert_unique(mob_profiles, c("Participant_id", "tp_lab"), "MOB episode profiles")
assert_unique(mechanism_profiles, c("Participant_id", "tp_lab"), "Plasmid mechanism profiles")
assert_value(nrow(mob_profiles), 532L, "MOB episode denominator")
assert_value(nrow(mechanism_profiles), 532L, "Plasmid mechanism denominator")
if (any(!location_validation$pass)) {
  stop("Predicted-plasmid location validation failed; FigS07 cannot be published.")
}

episode_mechanism <- cohort %>% select(Participant_id, tp_lab, status_display) %>%
  left_join(
    amr %>% select(
      Participant_id, tp_lab,
      amr_gene_count = amr_gene_count_informative,
      mdfA_detected
    ),
    by = c("Participant_id", "tp_lab"), relationship = "one-to-one"
  ) %>%
  left_join(
    mechanism_profiles %>% select(
      Participant_id, tp_lab, predicted_plasmid_count,
      plasmid_binned_informative_vf_amr_burden,
      mob_high_confidence_profile
    ),
    by = c("Participant_id", "tp_lab"), relationship = "one-to-one"
  ) %>%
  mutate(
    profile_state = if_else(
      is.na(.data$amr_gene_count) |
        is.na(.data$predicted_plasmid_count) |
        is.na(.data$plasmid_binned_informative_vf_amr_burden),
      "Unavailable", "Available"
    )
  )
if (any(episode_mechanism$profile_state != "Available")) {
  stop("All 532 episodes require validated AMR and predicted-plasmid profiles.")
}
mechanism_denoms <- episode_mechanism %>%
  count(.data$status_display, name = "n")
mechanism_n <- setNames(
  mechanism_denoms$n[
    match(status_order, as.character(mechanism_denoms$status_display))
  ],
  status_order
)
if (anyNA(mechanism_n)) stop("Predicted-plasmid status denominators are incomplete.")
mechanism_status_labels <- setNames(
  sprintf("%s\n(n=%d episodes)", status_order, mechanism_n),
  status_order
)
mechanism_long <- episode_mechanism %>%
  select(
    Participant_id, tp_lab, status_display,
    predicted_plasmid_count,
    plasmid_binned_informative_vf_amr_burden
  ) %>%
  pivot_longer(
    c(
      predicted_plasmid_count,
      plasmid_binned_informative_vf_amr_burden
    ),
    names_to = "endpoint", values_to = "value"
  ) %>%
  mutate(
    endpoint = recode(
      endpoint,
      predicted_plasmid_count = "MOB predicted plasmid count",
      plasmid_binned_informative_vf_amr_burden =
        "Plasmid-binned informative VF/AMR burden"
    )
  )
pS07b <- ggplot(
  mechanism_long,
  aes(.data$status_display, .data$value, fill = .data$status_display)
) +
  geom_boxplot(width = .45, outlier.shape = NA, alpha = .5) +
  geom_point(
    position = position_jitter(width = .12, seed = 20260714),
    alpha = .25, size = .7
  ) +
  facet_wrap(~ endpoint, scales = "free_y", ncol = 2) +
  status_fill(name = "Operational UTI status", guide = "none") +
  scale_x_discrete(labels = mechanism_status_labels) +
  labs(
    title = "Predicted plasmid burden and localized gene context",
    subtitle = "532/532 MOB profiles; VF/AMR burden counts distinct informative features placed in predicted bins",
    x = "Operational UTI status and complete-profile denominator",
    y = "Episode-level count"
  ) +
  theme_ruti_publication(9.5) +
  theme(legend.position = "none", plot.margin = margin(8, 8, 8, 14))

focused_plasmid <- read_current(paths[["plasmid_focused"]]) %>%
  arrange(.data$Participant_id, .data$tp_from, .data$tp_to) %>%
  mutate(
    case = sprintf("Case %02d", row_number()),
    `Replicon profile` = if_else(
      .data$any_replicon_profile_change, "Changed", "Stable"
    ),
    `MOB cluster profile` = if_else(
      .data$any_mob_cluster_change, "Changed", "Stable"
    ),
    `Predicted plasmid count` = if_else(
      .data$predicted_plasmid_count_difference != 0, "Changed", "Stable"
    ),
    `Plasmid-binned VF` = if_else(
      nzchar(coalesce(.data$plasmid_binned_vf_genes_gained, "")) |
        nzchar(coalesce(.data$plasmid_binned_vf_genes_lost, "")),
      "Changed", "Stable"
    ),
    `Plasmid-binned AMR` = if_else(
      nzchar(coalesce(
        .data$plasmid_binned_informative_amr_genes_gained, ""
      )) |
        nzchar(coalesce(
          .data$plasmid_binned_informative_amr_genes_lost, ""
        )),
      "Changed", "Stable"
    )
  )
assert_value(nrow(focused_plasmid), 9L, "Focused plasmid transition denominator")
focused_heat <- focused_plasmid %>%
  select(
    case, `Replicon profile`, `MOB cluster profile`,
    `Predicted plasmid count`, `Plasmid-binned VF`,
    `Plasmid-binned AMR`
  ) %>%
  pivot_longer(-case, names_to = "mechanism", values_to = "state") %>%
  mutate(
    case = factor(case, levels = rev(unique(case))),
    state = factor(state, levels = c("Stable", "Changed"))
  )
pS07c <- ggplot(focused_heat, aes(.data$mechanism, .data$case, fill = .data$state)) +
  geom_tile(colour = "white", linewidth = .7) +
  geom_text(aes(label = .data$state), size = 2.6) +
  scale_fill_manual(
    values = c(Stable = "#DCE8F2", Changed = "#D55E00"),
    name = "Predicted change"
  ) +
  labs(
    title = "Nine descriptive Not UTI-to-UTI transitions",
    subtitle = "Deidentified cases; changes are assembly-based predicted context",
    x = NULL, y = NULL
  ) +
  theme_ruti_publication(9.5) +
  theme(
    axis.text.x = element_text(angle = 25, hjust = 1),
    panel.grid = element_blank()
  )

figS07 <- (pS07a / pS07b / pS07c) +
  plot_layout(heights = c(1.2, .95, 1.05), guides = "collect") +
  plot_annotation(
    title = "Replicon detection, predicted plasmid context and localized AMR/VF changes",
    subtitle = "Marker detection, reconstructed bins and predicted gene localization are shown as distinct evidence layers",
    tag_levels = "A", theme = panel_annotation_theme
  ) &
  theme(
    legend.position = "bottom",
    legend.box = "vertical",
    legend.box.just = "center"
  )
capS07 <- paste0(
  "Predicted plasmid and genomic-AMR mechanism context. Panel A shows corrected gene-level PlasmidFinder marker prevalence for markers detected in at least 5% of 532 complete episode profiles (",
  rep_status_labels[["UTI"]], "; ", rep_status_labels[["Not UTI"]],
  "); successful no-hit episodes remain valid zero profiles. Panel B separates MOB predicted plasmid-bin count from the count of distinct pinned CGE VirulenceFinder and primary informative AMRFinderPlus features placed in predicted plasmid bins (",
  mechanism_status_labels[["UTI"]] %>% str_replace("\\n", " "), "; ",
  mechanism_status_labels[["Not UTI"]] %>% str_replace("\\n", " "),
  "). Panel C shows descriptive changes across the nine deidentified adjacent Not UTI-to-UTI transitions; no regression is fitted. Episodes and adjacent pairs are repeated within participants and are not independent. MOB results are assembly-based predictions, not confirmed circular plasmids. Same-bin placement is predicted linkage only and does not establish transfer, transmission, phenotypic susceptibility, causality, or direction."
)
register_figure(
  "FigS07_plasmid_amr_context", "FigS07",
  "Predicted plasmid and AMR mechanism context",
  "How do replicon markers, reconstructed plasmid predictions and predicted gene localization contribute distinct mechanism context?",
  "Supplementary", figS07, 8.27, 13.5, capS07,
  paths[c(
    "plasmid_manifest", "replicons", "mob_profiles",
    "plasmid_mechanism_profiles", "plasmid_focused",
    "plasmid_location_validation", "amr", "mlst", "cohort"
  )],
  "Descriptive marker prevalence, predicted burden and nine-case changes",
  "Assembly-only predictions; no circularity, physical transfer, transmission, phenotype or causality claim.",
  "Episode and adjacent transition",
  "532 complete profiles and nine deidentified Not UTI-to-UTI transitions",
  "Panel A shows marker prevalence; Panel B shows predicted bin/gene burden; Panel C shows descriptive transition changes.",
  "None"
)

# ----------------------------------------------------------------------------
# FigS08. Full prescreened gene-model evidence
# ----------------------------------------------------------------------------
gm <- read_current(paths[["gene_models"]]) %>%
  mutate(
    estimable = is.finite(.data$OR) & is.finite(.data$OR_lower) & is.finite(.data$OR_upper) & .data$OR > 0 & .data$OR_lower > 0 & .data$OR_upper > 0,
    sparse = .data$sparse_data_separation_risk %in% TRUE,
    singular = str_detect(.data$model_type, fixed("Singular")),
    display_or = if_else(.data$estimable, pmin(pmax(.data$OR, .01), 100), 1),
    display_lo = if_else(.data$estimable, pmin(pmax(.data$OR_lower, .01), 100), 1),
    display_hi = if_else(.data$estimable, pmin(pmax(.data$OR_upper, .01), 100), 1),
    feature_label = factor(.data$feature, levels = rev(.data$feature[order(.data$OR, na.last = TRUE)])),
    estimate_state = case_when(!.data$estimable ~ "Not estimable", .data$sparse ~ "Sparse/separation risk", .data$singular ~ "Singular fit", TRUE ~ "Finite model")
  )
assert_value(nrow(gm), 50L, "Prescreened gene-model result count")
gm_cols <- c("Finite model" = "#0072B2", "Singular fit" = "#E69F00", "Sparse/separation risk" = "#CC79A7", "Not estimable" = "#6B7280")
gm_shapes <- c("Finite model" = 16, "Singular fit" = 17, "Sparse/separation risk" = 18, "Not estimable" = 4)
pS08 <- ggplot(gm, aes(.data$display_or, .data$feature_label, colour = .data$estimate_state)) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey35") +
  geom_errorbarh(aes(xmin = .data$display_lo, xmax = .data$display_hi), height = .14, linewidth = .65) +
  geom_point(aes(shape = .data$estimate_state), size = 2.45, stroke = .8) +
  scale_x_log10(limits = c(.01, 100), breaks = c(.01, .1, 1, 10, 100)) +
  scale_colour_manual(values = gm_cols, name = "Model diagnostic") +
  scale_shape_manual(values = gm_shapes, name = "Model diagnostic") +
  labs(title = "Prescreened gene-level model evidence", subtitle = "All 50 fitted features; CIs are clipped at 0.01 and 100, and colour plus shape flag model diagnostics",
       x = "Adjusted odds ratio for UTI versus Not UTI (log scale)", y = "Virulence-factor gene / feature",
       caption = str_wrap("Inclusion rule from script 14: Fisher p<0.1 or the first 50 Fisher-ranked features. No model was BH-FDR significant; prescreening and sparse fits make this exploratory.", width = 118)) +
  theme_ruti_publication(9) +
  theme(
    legend.position = "bottom",
    axis.text.y = element_text(size = 8.5),
    legend.text = element_text(size = 8.3),
    plot.caption = element_text(size = 8.1, lineheight = 1.05)
  )
capS08 <- "Gene-level association evidence for all 50 features fitted after the declared exploratory prescreen (Fisher p<0.1 or the first 50 Fisher-ranked features). Points are adjusted odds ratios for operational UTI, horizontal lines are Wald 95% confidence intervals, and the logarithmic x-axis is truncated at 0.01 and 100 for display. Non-estimable or separated results are shown at the null with an x rather than plotted as finite effects; colour and shape redundantly flag finite, singular, sparse/separation-risk, and non-estimable fits. Models use participant random intercepts where GLMM fitting succeeded and adjust for timepoint and batch. Benjamini-Hochberg correction was applied across fitted features; none had FDR<0.05. This is exploratory prescreened evidence, not an independent confirmatory analysis."
register_figure("FigS08_gene_model_forest", "FigS08", "Prescreened gene-level model evidence",
                "What does the full fitted gene-level evidence show, including unstable estimates?", "Supplementary", pS08, 8.27, 11.69, capS08,
                paths[c("gene_models", "fisher")], "Prescreened gene-level GLMM/declared fallback model",
                "Prescreening, singularity, separation and repeated episodes limit inference; no FDR-significant gene.", "Episode with participant-aware model when estimable",
                "Fisher p<0.1 or first 50 Fisher-ranked features", "Points are ORs, lines are 95% CIs, and colour plus shape flag model diagnostics.",
                "Benjamini-Hochberg across fitted gene models")

# ----------------------------------------------------------------------------
# FigS09. UTI-event paired sensitivity
# ----------------------------------------------------------------------------
rq07inf <- read_current(paths[["rq07_inference"]]) %>% filter(.data$analysis == "all_UTI_event_samples")
rq07pair <- read_current(paths[["rq07_paired"]]) %>% mutate(Participant_id = as.character(.data$Participant_id)) %>%
  pivot_longer(c("curated_delta_uti_minus_not", "expec_delta_uti_minus_not"), names_to = "endpoint", values_to = "delta") %>%
  mutate(endpoint_label = recode(.data$endpoint, curated_delta_uti_minus_not = "Curated VF genes", expec_delta_uti_minus_not = "ExPEC-like markers"))
pS09 <- ggplot(rq07pair, aes(.data$delta, .data$endpoint_label)) + geom_vline(xintercept = 0, linetype = "dashed") +
  geom_point(position = position_jitter(height = .08, seed = 20260714), alpha = .65, size = 2, colour = "#4C78A8") +
  stat_summary(fun = median, geom = "point", shape = 23, fill = "#D55E00", colour = "black", size = 3.5) +
  labs(title = "Nearest within-participant event-sample sensitivity", subtitle = sprintf("%d paired participants; each point is UTI minus nearest Not UTI", n_distinct(rq07pair$Participant_id)),
       x = "Within-participant difference", y = NULL,
       caption = "Pairs are selected by temporal proximity and do not imply uninterrupted carriage. Diamonds are medians; inference remains exploratory.") + theme_ruti_publication(10)
capS09 <- paste0(
  "Nearest within-participant sensitivity analysis for ", n_distinct(rq07pair$Participant_id), " participants with an operational UTI among UTI-related samples and a temporally nearest Not UTI comparator. Points are within-participant UTI-minus-Not-UTI differences for curated VF-gene count and ExPEC-like marker count; diamonds are medians and the dashed line marks no difference. The prespecified all-event-sample analysis included 32 UTI-related samples from 29 participants (15 UTI; 17 Not UTI) with resident-cluster bootstrap uncertainty and Holm adjustment across two frozen endpoints. Temporal pairing reduces, but does not eliminate, confounding and does not establish continuous carriage."
)
register_figure("FigS09_event_sample_sensitivity", "FigS09", "UTI-event paired VF sensitivity",
                "Do nearest within-participant UTI-event comparisons support a VF burden difference?", "Supplementary", pS09, 8.5, 5.5, capS09,
                paths[c("rq07_inference", "rq07_paired")], "Nearest within-participant descriptive contrast; resident-cluster bootstrap in companion table",
                "Temporal pairing does not imply continuous carriage and includes a sparse subset.", "Participant pair",
                "Nearest available UTI and Not UTI UTI-event samples", "Points are participant differences; diamonds are medians.", "Holm across two frozen RQ07 endpoints")

# ----------------------------------------------------------------------------
# FigS10. SNP-threshold sensitivity
# ----------------------------------------------------------------------------
thr <- read_current(paths[["rq01_threshold"]])
pS10 <- ggplot(thr, aes(.data$threshold, .data$estimate)) +
  geom_ribbon(aes(ymin = .data$conf_low, ymax = .data$conf_high), fill = "#0072B2", alpha = .18) +
  geom_line(colour = "#0072B2", linewidth = .9) + geom_point(colour = "#0072B2", size = 2.6) +
  geom_vline(xintercept = SAME_STRAIN_SNP_THRESHOLD, linetype = "dashed", colour = "#D55E00") +
  scale_y_continuous(labels = label_percent()) +
  labs(title = "Sensitivity to the operational SNP reference", subtitle = "Resident-cluster bootstrap proportions across prespecified thresholds",
       x = "Pairwise SNP threshold", y = "Proportion of 371 adjacent pairs at or below threshold",
       caption = "Ribbon shows the 95% resident-cluster bootstrap interval; the orange line marks the prespecified 25-SNP reference.") + theme_ruti_publication(10)
capS10 <- "Threshold sensitivity for the proportion of 371 adjacent direct comparisons from 139 participants classified at or below each prespecified SNP threshold. Points and line show the observed proportions; the ribbon is the 95% resident-cluster bootstrap interval from 10,000 resamples, and the orange dashed line marks the operational 25-SNP reference. The analysis describes how a threshold changes classification and does not validate a universal biological boundary or imply transmission."
register_figure("FigS10_snp_threshold_sensitivity", "FigS10", "SNP-threshold sensitivity",
                "How does the classified continuity proportion change across SNP thresholds?", "Supplementary", pS10, 8.5, 5.5, capS10,
                paths[c("rq01_threshold", "transitions")], "Resident-cluster bootstrap proportion",
                "Thresholds are operational continuity references, not universal biological boundaries.", "Adjacent pair nested within participant",
                "371 adjacent comparisons; prespecified thresholds", "Line/points are observed proportions; ribbon is 95% cluster-bootstrap interval.", "None")

# ----------------------------------------------------------------------------
# Structured outputs
# ----------------------------------------------------------------------------
write_csv(manifest, file.path(DIR_FIG_AUDIT, "final_figure_manifest.csv"))

caption_lines <- c(
  "# Final thesis figure captions",
  "",
  paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  "",
  "Canonical analytical contract: 532 selected QC-passing Longcycler episodes from 161 participants (16 operational UTI; 516 heterogeneous Not UTI).",
  ""
)
for (i in seq_len(nrow(manifest))) {
  caption_lines <- c(caption_lines,
    paste0("## ", manifest$figure_number[[i]], ". ", manifest$title[[i]]), "",
    manifest$caption[[i]], "",
    paste0("Statistical/scientific method: ", manifest$statistical_method[[i]], "."),
    paste0("Required caveat: ", manifest$caveat[[i]]),
    paste0("Files: `", manifest$png_path[[i]], "`; `", manifest$pdf_path[[i]], "`."), "")
}
if (!variant_valid) {
  caption_lines <- c(caption_lines,
    "## Fig08. Reference-aware within-host variant map - withheld", "",
    "Fig08 remains unvalidated and is not part of the final set because the exact reference/contig coordinate checks did not all pass. Earlier genomic-position maps are obsolete and must not be substituted.", "")
}
writeLines(caption_lines, file.path(DIR_FIG_AUDIT, "final_figure_captions.md"))

checks <- tibble(
  check = c(
    "selected episodes", "selected participants", "operational UTI", "operational Not UTI",
    "direct pairs", "adjacent pairs", "Not UTI-to-UTI transitions", "same-strain-supported transitions",
    "MLST rows", "VF-ready rows", "tree tips", "tree-map one-to-one", "all expected figure files exist",
    "all expected figure files nonempty", "no raw participant identifiers in manifest captions"
  ),
  observed = c(
    nrow(cohort), n_distinct(cohort$Participant_id), sum(cohort$UTI_Status == "UTI"), sum(cohort$UTI_Status == "Not_UTI"),
    nrow(pairwise), nrow(transitions), sum(transitions$status_from == "Not_UTI" & transitions$status_to == "UTI"),
    sum(transitions$status_from == "Not_UTI" & transitions$status_to == "UTI" & transitions$TotalSNPs <= SAME_STRAIN_SNP_THRESHOLD),
    nrow(mlst), nrow(vf_ready), length(tree$tip.label), as.integer(setequal(tree$tip.label, tree_map$parsnp_alignment_label)),
    sum(file.exists(c(manifest$png_path, manifest$pdf_path))),
    sum(file.info(c(manifest$png_path, manifest$pdf_path))$size > 0),
    sum(str_detect(paste(manifest$caption, collapse = " "), "[0-9]{5,}"))
  ),
  expected = c(532, 161, 16, 516, 893, 371, 9, 5, 532, 532, 532, 1,
               2 * nrow(manifest), 2 * nrow(manifest), 0)
) %>% mutate(pass = .data$observed == .data$expected)
write_csv(checks, file.path(DIR_FIG_AUDIT, "final_figure_data_checks.csv"))
if (any(!checks$pass)) stop("Final figure data checks failed; see results/figure_audit/final_figure_data_checks.csv.")

msg("Canonical thesis figure pack complete: %d retained figure families (%d main, %d supplementary).",
    nrow(manifest), sum(manifest$figure_class == "Main"), sum(manifest$figure_class == "Supplementary"))

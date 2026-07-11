#!/usr/bin/env Rscript
# ==============================================================================
# 28_vf_transition_case_studies.R
# ==============================================================================
#
# GOAL:
#   Build clinical-first Not_UTI->UTI and phenotype-switch case-study tables. Every
#   ordered clinical transition is retained, including cases with missing WGS/VF,
#   and genomic/VF/module/score evidence is layered on only where available.
#
# METHOD:
#   1. Load ordered primary UTI status data, canonical VF data, module outputs,
#      score outputs, longitudinal outputs, and strain metrics.
#   2. Build a complete ordered transition index from clinical episodes.
#   3. Mark WGS/VF/module/score/strain availability for both endpoints.
#   4. Calculate VF gene, module, and score changes for WGS-linked pairs.
#   5. Classify same-strain/replacement/missing-genomics evidence.
#   6. Write case-study tables, notes, summary text, and descriptive figures.
#
# OUTPUT:
#   - results/vf/vf_transition_case_index.csv
#   - results/vf/vf_transition_case_summary.csv
#   - results/vf/vf_transition_gene_changes.csv
#   - results/vf/vf_transition_module_changes.csv
#   - results/vf/vf_transition_score_changes.csv
#   - results/vf/vf_transition_strain_context.csv
#   - results/vf/vf_transition_case_notes.csv
#   - results/vf/vf_transition_case_study_summary.txt
#   - plots/vf/vf_transition_case_timeline.png
#   - plots/vf/vf_transition_module_change_heatmap.png
#   - plots/vf/vf_transition_score_slopeplot.png
#   - plots/vf/vf_transition_gene_gain_loss_tile.png
#   - plots/vf/vf_transition_snp_vs_vf_jaccard.png
#
# ==============================================================================

source("00_config.R")
suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
  library(tibble)
  library(lubridate)
})

msg("Starting 28_vf_transition_case_studies.R")
ensure_dir(DIR_VF)
ensure_dir(DIR_PLOTS_VF)

SNP_THRESHOLD <- strain_snp_threshold()

# ==============================================================================
# 1. HELPERS
# ==============================================================================
normalise_tp_label <- function(x) {
  normalise_timepoint_preserve_events(x)
}

fallback_time_order <- function(tp) {
  tp <- as.character(tp)
  case_when(
    str_detect(tp, regex("^T\\d+$", ignore_case = TRUE)) ~
      suppressWarnings(as.numeric(str_remove(str_to_upper(tp), "^T"))),
    str_detect(tp, regex("uricult", ignore_case = TRUE)) ~ 99,
    TRUE ~ NA_real_
  )
}

parse_date_safe <- function(x) {
  x <- as.character(x)
  out <- suppressWarnings(lubridate::dmy(x))
  if (all(is.na(out))) out <- suppressWarnings(lubridate::ymd(x))
  as.Date(out)
}

normalise_st_label <- function(x) {
  x <- as.character(x)
  x <- str_trim(x)
  unknown <- c("", "-", "ST-", "NA", "N/A", "UNKNOWN", "UNK", "NT",
               "NON-TYPABLE", "NONTYPABLE", "NOT TYPED")
  x[str_to_upper(x) %in% unknown] <- NA_character_
  x
}

first_or_na <- function(x) {
  if (length(x) == 0) return(NA)
  x[[1]]
}

jaccard_from_vectors <- function(from_vec, to_vec) {
  union_size <- sum(from_vec == 1 | to_vec == 1, na.rm = TRUE)
  inter_size <- sum(from_vec == 1 & to_vec == 1, na.rm = TRUE)
  if (union_size > 0) round(inter_size / union_size, 4) else NA_real_
}

plot_theme_vf <- function(base_size = 10) {
  theme_bw(base_size = base_size) +
    theme(
      plot.caption = element_text(hjust = 0, size = base_size - 2, colour = "grey35"),
      plot.subtitle = element_text(colour = "grey25"),
      legend.position = "bottom"
    )
}

stop_if_stale <- function(target, upstream, target_label, upstream_label) {
  if (!file.exists(target) || !file.exists(upstream)) return(invisible(FALSE))
  target_mtime <- file.info(target)$mtime
  upstream_mtime <- file.info(upstream)$mtime
  if (!is.na(target_mtime) && !is.na(upstream_mtime) && target_mtime < upstream_mtime) {
    if (basename(target) == "vf_analysis_ready.csv" && basename(upstream) == "vf_pa_all.csv") {
      target_probe <- read_csv(target, show_col_types = FALSE, n_max = Inf)
      upstream_probe <- read_csv(upstream, show_col_types = FALSE, n_max = Inf)
      target_genes <- canonical_vf_gene_cols(names(target_probe), vf_pa_file = upstream)
      upstream_genes <- canonical_vf_gene_cols(names(upstream_probe), vf_pa_file = upstream)
      if (nrow(target_probe) == nrow(upstream_probe) && setequal(target_genes, upstream_genes)) {
        msg("WARNING: %s timestamp is older than %s, but row count and gene set match; continuing.",
            target_label, upstream_label)
        return(invisible(FALSE))
      }
    }
    stop(sprintf(
      "%s is older than %s. Re-run 22_vf_build_analysis_dataset.R before script 28.\n  %s: %s\n  %s: %s",
      target_label, upstream_label, target_label, format(target_mtime),
      upstream_label, format(upstream_mtime)
    ))
  }
}

# ==============================================================================
# 2. LOAD REQUIRED INPUTS
# ==============================================================================
if (!file.exists(FILE_VF_READY)) stop("Missing ", FILE_VF_READY)
stop_if_stale(FILE_VF_READY, FILE_VF_PA, "vf_analysis_ready.csv", "vf_pa_all.csv")

vf <- read_csv(FILE_VF_READY, show_col_types = FALSE) %>%
  prefer_primary_uti_status() %>%
  apply_manual_sample_curation(context = "28_vf_ready") %>%
  filter_primary_genomics() %>%
  mutate(
    Participant_id = as.character(Participant_id),
    tp_lab = normalise_tp_label(tp_lab)
  )

vf_episode_lookup <- vf %>%
  filter(!is.na(Episode_ID), !is.na(Infection_Status)) %>%
  transmute(
    Participant_id,
    Episode_ID = as.character(Episode_ID),
    vf_tp_lab = tp_lab,
    vf_bridge_applied = if ("uricult_bridge_applied" %in% names(.)) uricult_bridge_applied %in% TRUE else FALSE
  ) %>%
  distinct()

vf_episode_dupes <- vf_episode_lookup %>%
  count(Participant_id, Episode_ID, name = "n") %>%
  filter(n > 1)
if (nrow(vf_episode_dupes) > 0) {
  write_csv(vf_episode_dupes, file.path(DIR_QC, "vf_transition_duplicate_episode_lookup.csv"))
  stop("vf_analysis_ready has duplicated Participant_id + Episode_ID mappings; refusing ambiguous transition linkage.")
}

vf_key_set <- paste(vf$Participant_id, vf$tp_lab, sep = "|")

meta_cols <- c("Participant_id", "tp_lab", "Episode_ID", "Event_type", "Collection_Date",
               "Infection_Status", "Batch", "ST", "vf_count_total",
               "total_vf_count_all", "total_vf_count_curated",
               "total_vf_count_unassigned", "low_confidence_count", "n_timepoints")
cat_cols <- grep("^cat_", names(vf), value = TRUE)
gene_cols <- canonical_vf_gene_cols(names(vf))

status_file <- select_primary_status_map(
  prefer_poster = TRUE,
  require_fresh = TRUE,
  caller = "28_vf_transition_case_studies.R"
)
status_plain <- FILE_STATUS_MAP

clinical_raw <- read_csv(status_file, show_col_types = FALSE)
clinical_raw <- prefer_primary_uti_status(clinical_raw, allow_legacy_fallback = FALSE) %>%
  apply_manual_sample_curation(context = "28_clinical_status") %>%
  filter_primary_analysis()
if (!"Collection_Date" %in% names(clinical_raw) && file.exists(status_plain)) {
  status_dates <- read_csv(status_plain, show_col_types = FALSE) %>%
    apply_manual_sample_curation(context = "28_status_dates") %>%
    filter_primary_analysis() %>%
    mutate(Participant_id = as.character(Participant_id),
           Episode_ID = as.character(Episode_ID)) %>%
    select(any_of(c("Participant_id", "Episode_ID", "Collection_Date")))
  clinical_raw <- clinical_raw %>%
    mutate(Participant_id = as.character(Participant_id),
           Episode_ID = as.character(Episode_ID)) %>%
    left_join(status_dates, by = c("Participant_id", "Episode_ID"), relationship = "many-to-one")
}

clinical <- clinical_raw %>%
  mutate(
    Participant_id = as.character(Participant_id),
    Episode_ID = as.character(Episode_ID),
    tp_lab = if ("tp_lab" %in% names(.)) normalise_tp_label(tp_lab) else normalise_tp_label(Timepoint),
    Collection_Date_parsed = if ("Collection_Date" %in% names(.)) parse_date_safe(Collection_Date) else as.Date(NA),
    Plot_TP_Num_Poster = if ("Plot_TP_Num_Poster" %in% names(.)) as.numeric(Plot_TP_Num_Poster) else NA_real_,
    Plot_TP_Label_Poster = if ("Plot_TP_Label_Poster" %in% names(.)) as.character(Plot_TP_Label_Poster) else NA_character_,
    Placement_Confidence = if ("Placement_Confidence" %in% names(.)) as.character(Placement_Confidence) else NA_character_
  )

msg("Loaded clinical transition source: %s (%d rows)", basename(status_file), nrow(clinical))
msg("Loaded VF-ready data: %d rows, %d VF genes", nrow(vf), length(gene_cols))

# ==============================================================================
# 3. LOAD OPTIONAL INPUTS
# ==============================================================================
f_evol <- file.path(DIR_RESULTS, "longitudinal", "evolution_events.csv")
evol <- if (file.exists(f_evol)) {
  read_csv(f_evol, show_col_types = FALSE) %>%
    mutate(Participant_id = as.character(Participant_id),
           From_Time = normalise_tp_label(From_Time),
           To_Time = normalise_tp_label(To_Time))
} else NULL

f_pair <- file.path(DIR_STRAIN, "pairwise_metrics.csv")
pairwise <- if (file.exists(f_pair)) {
  read_csv(f_pair, show_col_types = FALSE) %>%
    mutate(
      Participant_id_A = as.character(Participant_id_A),
      Participant_id_B = as.character(Participant_id_B),
      Timepoint_A = normalise_tp_label(Timepoint_A),
      Timepoint_B = normalise_tp_label(Timepoint_B)
    ) %>%
    prepare_pairwise_for_strain_context()
} else NULL

f_switch <- file.path(DIR_RESULTS, "longitudinal", "phenotype_switch_candidates.csv")
switches <- if (file.exists(f_switch)) {
  read_csv(f_switch, show_col_types = FALSE) %>%
    mutate(Participant_id = as.character(Participant_id),
           From_Time = normalise_tp_label(From_Time),
           To_Time = normalise_tp_label(To_Time),
           switch_key = paste(Participant_id, From_Time, To_Time, sep = "|"))
} else tibble(switch_key = character())

f_modmap <- file.path(DIR_VF, "gene_module_map.csv")
f_modep <- file.path(DIR_VF, "vf_module_presence_by_episode.csv")
has_mod <- file.exists(f_modmap) && file.exists(f_modep)
mod_map <- if (has_mod) read_csv(f_modmap, show_col_types = FALSE) else tibble()
mod_ep <- if (has_mod) {
  read_csv(f_modep, show_col_types = FALSE) %>%
    prefer_primary_uti_status() %>%
    mutate(Participant_id = as.character(Participant_id), tp_lab = normalise_tp_label(tp_lab))
} else tibble()

f_scores <- file.path(DIR_VF, "vf_score_table.csv")
has_scores <- file.exists(f_scores)
score_tbl <- if (has_scores) {
  read_csv(f_scores, show_col_types = FALSE) %>%
    prefer_primary_uti_status() %>%
    mutate(Participant_id = as.character(Participant_id), tp_lab = normalise_tp_label(tp_lab))
} else tibble()

# ==============================================================================
# 4. BUILD COMPLETE CLINICAL TRANSITION INDEX
# ==============================================================================
clinical_ordered <- clinical %>%
  group_by(Participant_id) %>%
  arrange(Participant_id, Collection_Date_parsed, Plot_TP_Num_Poster, fallback_time_order(tp_lab), tp_lab, .by_group = TRUE) %>%
  mutate(
    date_order = ifelse(!is.na(Collection_Date_parsed),
                        as.numeric(Collection_Date_parsed - min(Collection_Date_parsed, na.rm = TRUE)),
                        NA_real_),
    fallback_order = fallback_time_order(tp_lab),
    time_order = coalesce(date_order, Plot_TP_Num_Poster, fallback_order),
    time_order_source = case_when(
      !is.na(date_order) ~ "Collection_Date",
      !is.na(Plot_TP_Num_Poster) ~ "Plot_TP_Num_Poster",
      !is.na(fallback_order) & tp_lab == "Uricult" ~ "Fallback_Uricult_last",
      !is.na(fallback_order) ~ "Fallback_routine_timepoint",
      TRUE ~ "Unavailable"
    )
  ) %>%
  filter(!is.na(time_order)) %>%
  arrange(Participant_id, time_order, tp_lab, .by_group = TRUE) %>%
  mutate(
    from_tp = lag(tp_lab),
    from_episode_id = lag(Episode_ID),
    from_timepoint_raw = lag(Timepoint),
    from_status = lag(Infection_Status),
    from_order = lag(time_order),
    from_order_source = lag(time_order_source),
    from_collection_date = lag(as.character(Collection_Date_parsed)),
    from_plot_label = lag(Plot_TP_Label_Poster),
    from_placement_confidence = lag(Placement_Confidence),
    to_tp = tp_lab,
    to_episode_id = Episode_ID,
    to_timepoint_raw = Timepoint,
    to_status = Infection_Status,
    to_order = time_order,
    to_order_source = time_order_source,
    to_collection_date = as.character(Collection_Date_parsed),
    to_plot_label = Plot_TP_Label_Poster,
    to_placement_confidence = Placement_Confidence,
    is_consecutive_observed = !is.na(from_tp)
  ) %>%
  ungroup() %>%
  filter(!is.na(from_tp), !is.na(from_status), !is.na(to_status)) %>%
  left_join(
    vf_episode_lookup %>%
      rename(from_vf_tp = vf_tp_lab, from_vf_bridge_applied = vf_bridge_applied),
    by = c("Participant_id", "from_episode_id" = "Episode_ID"),
    relationship = "many-to-one"
  ) %>%
  left_join(
    vf_episode_lookup %>%
      rename(to_vf_tp = vf_tp_lab, to_vf_bridge_applied = vf_bridge_applied),
    by = c("Participant_id", "to_episode_id" = "Episode_ID"),
    relationship = "many-to-one"
  ) %>%
  mutate(
    from_vf_tp = coalesce(from_vf_tp, from_tp),
    to_vf_tp = coalesce(to_vf_tp, to_tp),
    transition_type = paste(from_status, to_status, sep = "->"),
    transition_key = paste(Participant_id, from_tp, to_tp, sep = "|"),
    wgs_transition_key = paste(Participant_id, from_vf_tp, to_vf_tp, sep = "|"),
    is_not_uti_to_uti = from_status == "Not_UTI" & to_status == "UTI",
    is_uricult_transition = str_detect(from_tp, regex("uricult", ignore_case = TRUE)) |
      str_detect(to_tp, regex("uricult", ignore_case = TRUE)),
    time_order_source = paste(from_order_source, to_order_source, sep = " -> "),
    timing_caveat = case_when(
      str_detect(time_order_source, "Fallback_Uricult_last") ~ "Uricult ordered by fallback convention",
      str_detect(time_order_source, "Plot_TP_Num_Poster") ~ "Poster/display timepoint used for ordering where date unavailable",
      TRUE ~ NA_character_
    ),
    from_key = paste(Participant_id, from_vf_tp, sep = "|"),
    to_key = paste(Participant_id, to_vf_tp, sep = "|"),
    has_wgs_from = from_key %in% vf_key_set,
    has_wgs_to = to_key %in% vf_key_set,
    has_vf_pair = has_wgs_from & has_wgs_to,
    has_module_pair = if (has_mod) {
      from_key %in% paste(mod_ep$Participant_id, mod_ep$tp_lab, sep = "|") &
        to_key %in% paste(mod_ep$Participant_id, mod_ep$tp_lab, sep = "|")
    } else FALSE,
    has_score_pair = if (has_scores) {
      from_key %in% paste(score_tbl$Participant_id, score_tbl$tp_lab, sep = "|") &
        to_key %in% paste(score_tbl$Participant_id, score_tbl$tp_lab, sep = "|")
    } else FALSE,
    in_script15_switch_candidates = transition_key %in% switches$switch_key,
    case_priority = case_when(
      is_not_uti_to_uti & has_vf_pair ~ "Primary_Not_UTI_to_UTI_WGS_linked",
      is_not_uti_to_uti & !has_vf_pair ~ "Not_UTI_to_UTI_missing_genomics",
      in_script15_switch_candidates ~ "Script15_phenotype_switch_context",
      TRUE ~ "Context_transition"
    )
  ) %>%
  mutate(case_id = sprintf("case_%03d", row_number())) %>%
  mutate(
    From_Time = from_tp,
    To_Time = to_tp,
    From_Status = from_status,
    To_Status = to_status
  ) %>%
  select(case_id, Participant_id,
         From_Time, To_Time, From_Status, To_Status,
         from_tp, to_tp, from_vf_tp, to_vf_tp,
         from_episode_id, to_episode_id, from_vf_bridge_applied, to_vf_bridge_applied,
         from_order, to_order, time_order_source,
         from_collection_date, to_collection_date,
         from_plot_label, to_plot_label,
         from_placement_confidence, to_placement_confidence,
         transition_type, transition_key, wgs_transition_key,
         is_not_uti_to_uti, is_consecutive_observed,
         has_wgs_from, has_wgs_to, has_vf_pair, has_module_pair, has_score_pair,
         is_uricult_transition, timing_caveat,
         in_script15_switch_candidates, case_priority)

case_index <- clinical_ordered

msg("Indexed %d clinical transitions; Not_UTI->UTI=%d; WGS-linked Not_UTI->UTI=%d",
    nrow(case_index), sum(case_index$is_not_uti_to_uti, na.rm = TRUE),
    sum(case_index$is_not_uti_to_uti & case_index$has_vf_pair, na.rm = TRUE))

# ==============================================================================
# 5. STRAIN CONTEXT
# ==============================================================================
strain_ctx <- tibble()

for (i in seq_len(nrow(case_index))) {
  row <- case_index[i, ]
  pid <- row$Participant_id
  tp_f <- row$from_vf_tp
  tp_t <- row$to_vf_tp

  vf_from <- vf %>% filter(Participant_id == pid, tp_lab == tp_f)
  vf_to <- vf %>% filter(Participant_id == pid, tp_lab == tp_t)
  st_from <- if (nrow(vf_from) > 0) normalise_st_label(vf_from$ST[1]) else NA_character_
  st_to <- if (nrow(vf_to) > 0) normalise_st_label(vf_to$ST[1]) else NA_character_
  ctx <- lookup_strain_context(
    pairwise = pairwise,
    pid = pid,
    tp_from = tp_f,
    tp_to = tp_t,
    ST_from = st_from,
    ST_to = st_to,
    has_vf_pair = row$has_vf_pair,
    evol = evol,
    threshold = SNP_THRESHOLD
  )

  strain_ctx <- bind_rows(strain_ctx, tibble(
    case_id = row$case_id,
    Participant_id = pid,
    from_vf_tp = tp_f,
    to_vf_tp = tp_t,
    ST_from = ctx$ST_from,
    ST_to = ctx$ST_to,
    same_ST = ctx$same_ST,
    SNPs = ctx$SNPs,
    AvgIdentity = ctx$AvgIdentity,
    Pairwise_Classification = ctx$Pairwise_Classification,
    Pairwise_RuleUsed = ctx$Pairwise_RuleUsed,
    Classification = ctx$Classification,
    RuleUsed = ctx$RuleUsed,
    VF_Jaccard_pairwise = ctx$VF_Jaccard_pairwise,
    Inc_Jaccard = ctx$Inc_Jaccard,
    snp_strain_context = as.character(ctx$snp_strain_context),
    st_lineage_context = as.character(ctx$st_lineage_context),
    pair_interpretation = as.character(ctx$pair_interpretation),
    same_strain_evidence = ctx$same_strain_evidence,
    strain_context_level = as.character(ctx$strain_context_level),
    replacement_flag = ctx$replacement_flag,
    strain_context_note = ctx$strain_context_note,
    same_strain_snp_threshold = ctx$same_strain_snp_threshold
  ))
}

case_index <- case_index %>%
  left_join(strain_ctx %>% select(case_id, has_strain_context = pair_interpretation),
            by = "case_id") %>%
  mutate(has_strain_context = !is.na(has_strain_context) & has_strain_context != "Missing strain metrics")

# ==============================================================================
# 6. GENE, MODULE, AND SCORE CHANGES
# ==============================================================================
gene_changes <- tibble()
module_changes <- tibble()
score_changes <- tibble()
case_summary <- tibble()

module_pres_cols <- if (has_mod) grep("^mod_.*_present$", names(mod_ep), value = TRUE) else character()
score_names <- if (has_scores) {
  intersect(c("expec_marker_count", "upec_system_count", "upec_system_fraction",
              "total_vf_count_all", "total_vf_count_curated",
              "total_vf_count_unassigned", "low_confidence_count"), names(score_tbl))
} else character()
score_label_lookup <- c(
  expec_marker_count = "ExPEC-like marker count",
  upec_system_count = "UPEC systems present",
  upec_system_fraction = "UPEC system fraction",
  total_vf_count_all = "All VF genes",
  total_vf_count_curated = "Curated VF genes",
  total_vf_count_unassigned = "Unassigned VF genes",
  low_confidence_count = "Low-confidence VF genes"
)

for (i in seq_len(nrow(case_index))) {
  row <- case_index[i, ]
  pid <- row$Participant_id
  tp_f_clinical <- row$from_tp
  tp_t_clinical <- row$to_tp
  tp_f <- row$from_vf_tp
  tp_t <- row$to_vf_tp
  sc <- strain_ctx %>% filter(case_id == row$case_id)

  summary_row <- tibble(
    case_id = row$case_id,
    Participant_id = pid,
    transition_type = row$transition_type,
    from_tp = tp_f_clinical,
    to_tp = tp_t_clinical,
    from_vf_tp = tp_f,
    to_vf_tp = tp_t,
    from_status = row$From_Status,
    to_status = row$To_Status,
    is_not_uti_to_uti = row$is_not_uti_to_uti,
    has_wgs_from = row$has_wgs_from,
    has_wgs_to = row$has_wgs_to,
    has_vf_pair = row$has_vf_pair,
    has_module_pair = row$has_module_pair,
    has_score_pair = row$has_score_pair,
    is_uricult_transition = row$is_uricult_transition,
    in_script15_switch_candidates = row$in_script15_switch_candidates,
    timing_caveat = row$timing_caveat,
    ST_from = sc$ST_from,
    ST_to = sc$ST_to,
    same_ST = sc$same_ST,
    SNPs = sc$SNPs,
    AvgIdentity = sc$AvgIdentity,
    Pairwise_Classification = sc$Pairwise_Classification,
    Pairwise_RuleUsed = sc$Pairwise_RuleUsed,
    snp_strain_context = sc$snp_strain_context,
    st_lineage_context = sc$st_lineage_context,
    pair_interpretation = sc$pair_interpretation,
    same_strain_evidence = sc$same_strain_evidence,
    strain_context_level = sc$strain_context_level,
    replacement_flag = sc$replacement_flag,
    strain_context_note = sc$strain_context_note,
    same_strain_snp_threshold = sc$same_strain_snp_threshold,
    vf_count_from = NA_real_,
    vf_count_to = NA_real_,
    n_vf_genes_gained = NA_integer_,
    n_vf_genes_lost = NA_integer_,
    n_stable_present = NA_integer_,
    vf_jaccard = NA_real_,
    genes_gained = NA_character_,
    genes_lost = NA_character_,
    n_modules_gained = NA_integer_,
    n_modules_lost = NA_integer_,
    module_jaccard = NA_real_,
    modules_gained = NA_character_,
    modules_lost = NA_character_,
    delta_expec_marker_count = NA_real_,
    delta_upec_system_count = NA_real_,
    delta_upec_system_fraction = NA_real_,
    delta_total_vf_burden = NA_real_,
    delta_curated_vf_burden = NA_real_,
    delta_upec_module_score = NA_real_,
    case_class = NA_character_,
    interpretation_short = NA_character_,
    missing_data_note = NA_character_
  )

  if (row$has_vf_pair) {
    vf_from <- vf %>% filter(Participant_id == pid, tp_lab == tp_f) %>% slice(1)
    vf_to <- vf %>% filter(Participant_id == pid, tp_lab == tp_t) %>% slice(1)
    from_vec <- as.integer(vf_from[1, gene_cols])
    to_vec <- as.integer(vf_to[1, gene_cols])
    names(from_vec) <- gene_cols
    names(to_vec) <- gene_cols

    gained <- gene_cols[from_vec == 0 & to_vec == 1]
    lost <- gene_cols[from_vec == 1 & to_vec == 0]
    stable_present <- gene_cols[from_vec == 1 & to_vec == 1]
    vf_jaccard <- jaccard_from_vectors(from_vec, to_vec)

    changed_or_present <- unique(c(gained, lost, stable_present))
    if (length(changed_or_present) > 0) {
      gene_changes <- bind_rows(gene_changes, tibble(
        case_id = row$case_id,
        Participant_id = pid,
        Gene = changed_or_present,
        from_present = as.integer(changed_or_present %in% c(lost, stable_present)),
        to_present = as.integer(changed_or_present %in% c(gained, stable_present)),
        change_type = case_when(
          changed_or_present %in% gained ~ "gained",
          changed_or_present %in% lost ~ "lost",
          TRUE ~ "stable_present"
        )
      ) %>%
        left_join(
          mod_map %>%
            filter(primary_assignment) %>%
            select(Gene, Category, module_id, system_name, broad_module,
                   assignment_confidence, upec_score_candidate),
          by = "Gene"
        ) %>%
        mutate(
          biological_priority = case_when(
            change_type %in% c("gained", "lost") & upec_score_candidate %in% TRUE &
              assignment_confidence %in% c("High", "Moderate") ~ "high",
            change_type %in% c("gained", "lost") & assignment_confidence %in% c("High", "Moderate") ~ "medium",
            change_type %in% c("gained", "lost") ~ "low",
            TRUE ~ "stable_context"
          )
        ))
    }

    summary_row <- summary_row %>%
      mutate(
        vf_count_from = sum(from_vec, na.rm = TRUE),
        vf_count_to = sum(to_vec, na.rm = TRUE),
        n_vf_genes_gained = length(gained),
        n_vf_genes_lost = length(lost),
        n_stable_present = length(stable_present),
        vf_jaccard = vf_jaccard,
        genes_gained = paste(gained, collapse = "; "),
        genes_lost = paste(lost, collapse = "; ")
      )
  }

  if (row$has_module_pair && length(module_pres_cols) > 0) {
    me_from <- mod_ep %>% filter(Participant_id == pid, tp_lab == tp_f) %>% slice(1)
    me_to <- mod_ep %>% filter(Participant_id == pid, tp_lab == tp_t) %>% slice(1)
    from_mod <- as.integer(me_from[1, module_pres_cols])
    to_mod <- as.integer(me_to[1, module_pres_cols])
    names(from_mod) <- module_pres_cols
    names(to_mod) <- module_pres_cols
    mod_jaccard <- jaccard_from_vectors(from_mod, to_mod)

    mod_ids <- str_replace_all(str_replace_all(module_pres_cols, "^mod_", ""), "_present$", "")
    mod_gained <- mod_ids[from_mod == 0 & to_mod == 1]
    mod_lost <- mod_ids[from_mod == 1 & to_mod == 0]

    mod_change_tbl <- tibble(
      case_id = row$case_id,
      Participant_id = pid,
      module_col = module_pres_cols,
      module_id = mod_ids,
      from_present = from_mod,
      to_present = to_mod,
      change_type = case_when(
        from_mod == 0 & to_mod == 1 ~ "gained",
        from_mod == 1 & to_mod == 0 ~ "lost",
        from_mod == 1 & to_mod == 1 ~ "stable_present",
        TRUE ~ "stable_absent"
      )
    ) %>%
      left_join(
        mod_map %>%
          filter(primary_assignment) %>%
          distinct(module_id, system_name, broad_module, upec_score_candidate,
                   assignment_confidence),
        by = "module_id"
      ) %>%
      mutate(from_gene_count = NA_real_, to_gene_count = NA_real_)

    for (j in seq_len(nrow(mod_change_tbl))) {
      mod_id <- mod_change_tbl$module_id[j]
      count_col <- paste0("mod_", mod_id, "_n_genes")
      mod_change_tbl$from_gene_count[j] <- if (count_col %in% names(me_from)) me_from[[count_col]][1] else NA_real_
      mod_change_tbl$to_gene_count[j] <- if (count_col %in% names(me_to)) me_to[[count_col]][1] else NA_real_
    }

    module_changes <- bind_rows(module_changes, mod_change_tbl %>%
      mutate(
        delta_gene_count = to_gene_count - from_gene_count,
        module_change_confidence = case_when(
          change_type %in% c("gained", "lost") & assignment_confidence %in% c("High", "Moderate") ~ "high",
          change_type %in% c("gained", "lost") ~ "low",
          TRUE ~ "stable"
        )
      ))

    summary_row <- summary_row %>%
      mutate(
        n_modules_gained = length(mod_gained),
        n_modules_lost = length(mod_lost),
        module_jaccard = mod_jaccard,
        modules_gained = paste(mod_gained, collapse = "; "),
        modules_lost = paste(mod_lost, collapse = "; ")
      )
  }

  if (row$has_score_pair && length(score_names) > 0) {
    s_from <- score_tbl %>% filter(Participant_id == pid, tp_lab == tp_f) %>% slice(1)
    s_to <- score_tbl %>% filter(Participant_id == pid, tp_lab == tp_t) %>% slice(1)
    for (score_name in score_names) {
      fv <- s_from[[score_name]][1]
      tv <- s_to[[score_name]][1]
      score_changes <- bind_rows(score_changes, tibble(
        case_id = row$case_id,
        Participant_id = pid,
        score_name = score_name,
        score_label = ifelse(score_name %in% names(score_label_lookup),
                             unname(score_label_lookup[score_name]), score_name),
        from_value = fv,
        to_value = tv,
        delta = tv - fv,
        absolute_delta = abs(tv - fv),
        direction = case_when(tv > fv ~ "increased", tv < fv ~ "decreased", TRUE ~ "stable"),
        interpretation = "Descriptive transition-level supplementary endpoint change"
      ))
    }
    get_delta <- function(name_options) {
      found <- intersect(name_options, score_names)[1]
      if (is.na(found)) return(NA_real_)
      (s_to[[found]][1] - s_from[[found]][1])
    }
    summary_row <- summary_row %>%
      mutate(
        delta_expec_marker_count = get_delta("expec_marker_count"),
        delta_upec_system_count = get_delta("upec_system_count"),
        delta_upec_system_fraction = get_delta("upec_system_fraction"),
        delta_total_vf_burden = get_delta(c("total_vf_count_all", "vf_count_total")),
        delta_curated_vf_burden = get_delta(c("total_vf_count_curated", "vf_count_curated")),
        delta_upec_module_score = get_delta("upec_system_count")
      )
  }

  summary_row <- summary_row %>%
    mutate(
      missing_data_note = case_when(
        !has_wgs_from & !has_wgs_to ~ "Both transition endpoints lack WGS/VF rows",
        !has_wgs_from ~ "Starting endpoint lacks WGS/VF row",
        !has_wgs_to ~ "Endpoint lacks WGS/VF row",
        !has_module_pair ~ "WGS/VF pair present but module pair unavailable",
        !has_score_pair ~ "WGS/VF pair present but score pair unavailable",
        TRUE ~ NA_character_
      ),
      case_class = case_when(
        !has_vf_pair ~ "G: Clinical transition, missing genomic endpoint",
        pair_interpretation == "Conflict: SNP same-strain but ST differs" ~
          "H: SNP/ST conflict, manual review",
        pair_interpretation == "Replacement likely" ~ "D: Strain replacement",
        snp_strain_context == "Strong same strain" &
          coalesce(n_vf_genes_gained, 0L) == 0 & coalesce(n_vf_genes_lost, 0L) == 0 &
          coalesce(n_modules_gained, 0L) == 0 & coalesce(n_modules_lost, 0L) == 0 ~
          "A: Same strain, stable VF/module profile",
        snp_strain_context == "Strong same strain" &
          (coalesce(n_modules_gained, 0L) > 0 | coalesce(n_modules_lost, 0L) > 0) ~
          "B: Same strain, module change",
        snp_strain_context == "Strong same strain" &
          (coalesce(n_vf_genes_gained, 0L) > 0 | coalesce(n_vf_genes_lost, 0L) > 0) ~
          "C: Same strain, gene-level VF change",
        pair_interpretation == "Same lineage, not same strain by SNP" ~
          "E1: Same ST lineage, not same strain by SNP",
        pair_interpretation == "ST-consistent, SNP missing" ~
          "E2: ST-consistent, SNP missing",
        snp_strain_context == "Above same-strain SNP threshold" ~
          "E3: Above same-strain SNP threshold",
        pair_interpretation == "Missing strain metrics" ~
          "F: Missing strain metrics",
        TRUE ~ "F: Missing or uncertain SNP evidence"
      ),
      interpretation_short = case_when(
        !has_vf_pair ~ "Clinical transition retained, but VF/genomic comparison is unavailable.",
        case_class == "H: SNP/ST conflict, manual review" ~
          "SNP distance supports same strain, but ST differs; review typing and pairwise metrics manually.",
        case_class == "A: Same strain, stable VF/module profile" ~
          "Low-SNP/same-strain evidence with stable VF/module profile.",
        case_class == "D: Strain replacement" ~
          "Transition likely reflects strain replacement rather than within-strain VF evolution.",
        pair_interpretation == "Same lineage, not same strain by SNP" ~
          "Same ST lineage, but SNP distance exceeds the 25-SNP same-strain rule.",
        pair_interpretation == "ST-consistent, SNP missing" ~
          "ST is consistent across endpoints, but SNP evidence is missing and cannot prove same strain.",
        snp_strain_context == "Above same-strain SNP threshold" ~
          "SNP evidence does not meet the 25-SNP same-strain rule.",
        pair_interpretation == "Missing strain metrics" ~
          "WGS/VF endpoints are present, but pairwise same-strain metrics are unavailable.",
        str_detect(case_class, "change") ~
          "Same-strain evidence with VF/module change; review changed components.",
        TRUE ~ "Transition interpretation is uncertain."
      )
    )

  case_summary <- bind_rows(case_summary, summary_row)
}

# ==============================================================================
# 7. HUMAN-READABLE CASE NOTES
# ==============================================================================
case_notes <- case_summary %>%
  mutate(
    main_finding = interpretation_short,
    supporting_evidence = case_when(
      has_vf_pair ~ sprintf("SNP context=%s; ST context=%s; pair interpretation=%s; ST %s -> %s; SNPs=%s; VF Jaccard=%s; module Jaccard=%s; VF gained=%s lost=%s",
                            coalesce(snp_strain_context, "NA"),
                            coalesce(st_lineage_context, "NA"),
                            coalesce(pair_interpretation, "NA"),
                            coalesce(ST_from, "NA"), coalesce(ST_to, "NA"),
                            ifelse(is.na(SNPs), "NA", as.character(SNPs)),
                            ifelse(is.na(vf_jaccard), "NA", sprintf("%.3f", vf_jaccard)),
                            ifelse(is.na(module_jaccard), "NA", sprintf("%.3f", module_jaccard)),
                            ifelse(is.na(n_vf_genes_gained), "NA", as.character(n_vf_genes_gained)),
                            ifelse(is.na(n_vf_genes_lost), "NA", as.character(n_vf_genes_lost))),
      TRUE ~ "No paired VF/genomic endpoint available"
    ),
    limitations = coalesce(missing_data_note, timing_caveat, "Small case-study denominator; descriptive only"),
    recommended_use = case_when(
      is_not_uti_to_uti & has_vf_pair ~ "Primary thesis/manuscript transition case",
      is_not_uti_to_uti ~ "Missing-genomics denominator accounting",
      in_script15_switch_candidates ~ "Supplementary phenotype-switch context",
      TRUE ~ "Background longitudinal context"
    )
  ) %>%
  select(case_id, Participant_id, transition_type, case_class, main_finding,
         supporting_evidence, limitations, recommended_use)

# ==============================================================================
# 8. WRITE OUTPUT TABLES
# ==============================================================================
write_csv(case_index, file.path(DIR_VF, "vf_transition_case_index.csv"))
write_csv(case_summary, file.path(DIR_VF, "vf_transition_case_summary.csv"))
write_csv(gene_changes, file.path(DIR_VF, "vf_transition_gene_changes.csv"))
write_csv(module_changes, file.path(DIR_VF, "vf_transition_module_changes.csv"))
write_csv(score_changes, file.path(DIR_VF, "vf_transition_score_changes.csv"))
write_csv(strain_ctx, file.path(DIR_VF, "vf_transition_strain_context.csv"))
write_csv(case_notes, file.path(DIR_VF, "vf_transition_case_notes.csv"))

# ==============================================================================
# 9. SUMMARY REPORT
# ==============================================================================
txt <- character()
ta <- function(...) txt <<- c(txt, sprintf(...))

not_uti_uti_index <- case_index %>% filter(is_not_uti_to_uti)
not_uti_uti_summary <- case_summary %>% filter(is_not_uti_to_uti)

ta("=== VF TRANSITION CASE STUDY SUMMARY ===")
ta("Timestamp: %s", format(Sys.time()))
ta("Clinical transition source: %s", basename(status_file))
ta("Total ordered clinical transitions indexed: %d", nrow(case_index))
ta("Not_UTI->UTI clinical transitions indexed: %d", nrow(not_uti_uti_index))
ta("Not_UTI->UTI with both WGS/VF endpoints: %d", sum(not_uti_uti_index$has_vf_pair, na.rm = TRUE))
ta("Not_UTI->UTI with SNP-defined strong same-strain evidence: %d",
   sum(not_uti_uti_summary$has_vf_pair & not_uti_uti_summary$snp_strain_context == "Strong same strain", na.rm = TRUE))
ta("Not_UTI->UTI missing one or both WGS/VF endpoints: %d", sum(!not_uti_uti_index$has_vf_pair, na.rm = TRUE))
ta("SNP context rule: Strong same strain = 0-%d SNPs; Above same-strain SNP threshold = >%d SNPs; Missing SNP evidence = SNPs unavailable.", SNP_THRESHOLD, SNP_THRESHOLD)
ta("ST context rule: Same ST, Different ST, or Missing ST evidence. ST is secondary lineage context and does not prove same strain.")
ta("Uricult-linked Not_UTI->UTI transitions: %d", sum(not_uti_uti_index$is_uricult_transition, na.rm = TRUE))
ta("Transitions also present in script 15 phenotype-switch candidates: %d",
   sum(case_index$in_script15_switch_candidates, na.rm = TRUE))
ta("")
ta("--- Time-order sources ---")
order_counts <- case_index %>% dplyr::count(time_order_source, name = "n") %>% arrange(desc(n))
for (j in seq_len(nrow(order_counts))) ta("  %s: %d", order_counts$time_order_source[j], order_counts$n[j])
ta("")
if (nrow(not_uti_uti_summary) > 0) {
  ta("--- Not_UTI->UTI case class distribution ---")
  class_counts <- not_uti_uti_summary %>% dplyr::count(case_class, name = "n") %>% arrange(desc(n))
  for (j in seq_len(nrow(class_counts))) ta("  %s: %d", class_counts$case_class[j], class_counts$n[j])
  ta("")
  ta("--- Not_UTI->UTI WGS-linked details ---")
  detail_rows <- not_uti_uti_summary %>% filter(has_vf_pair)
  for (j in seq_len(nrow(detail_rows))) {
    r <- detail_rows[j, ]
    ta("  %s: participant %s %s->%s | SNP context=%s | ST context=%s | interpretation=%s | ST %s->%s | SNPs=%s | VF Jaccard=%s | modules gained=%s lost=%s | %s",
       r$case_id, r$Participant_id, r$from_tp, r$to_tp,
       coalesce(r$snp_strain_context, "NA"),
       coalesce(r$st_lineage_context, "NA"),
       coalesce(r$pair_interpretation, "NA"),
       coalesce(r$ST_from, "NA"), coalesce(r$ST_to, "NA"),
       ifelse(is.na(r$SNPs), "NA", as.character(r$SNPs)),
       ifelse(is.na(r$vf_jaccard), "NA", sprintf("%.3f", r$vf_jaccard)),
       ifelse(is.na(r$n_modules_gained), "NA", as.character(r$n_modules_gained)),
       ifelse(is.na(r$n_modules_lost), "NA", as.character(r$n_modules_lost)),
       r$case_class)
  }
}
ta("")
ta("CAUTION: Transition analyses are descriptive case studies, not predictive models.")
ta("Uricult ordering uses Collection_Date where available; poster/fallback ordering is reported in the index and should be cited as a caveat.")

writeLines(txt, file.path(DIR_VF, "vf_transition_case_study_summary.txt"))

# ==============================================================================
# 10. FIGURES
# ==============================================================================
not_uti_uti_ids <- not_uti_uti_index$case_id

# Clinical timeline for participants with Not_UTI->UTI transitions.
timeline_participants <- not_uti_uti_index %>% pull(Participant_id) %>% unique()
timeline_plot_data <- clinical %>%
  filter(Participant_id %in% timeline_participants) %>%
  left_join(
    vf_episode_lookup %>% select(Participant_id, Episode_ID, vf_tp_lab),
    by = c("Participant_id", "Episode_ID"),
    relationship = "many-to-one"
  ) %>%
  group_by(Participant_id) %>%
  mutate(
    date_order = ifelse(!is.na(Collection_Date_parsed),
                        as.numeric(Collection_Date_parsed - min(Collection_Date_parsed, na.rm = TRUE)),
                        NA_real_),
    fallback_order = fallback_time_order(tp_lab),
    time_order = coalesce(date_order, Plot_TP_Num_Poster, fallback_order),
    has_wgs_vf = (!is.na(vf_tp_lab) & paste(Participant_id, vf_tp_lab, sep = "|") %in% vf_key_set) |
      paste(Participant_id, tp_lab, sep = "|") %in% vf_key_set
  ) %>%
  ungroup() %>%
  filter(!is.na(time_order))

if (nrow(timeline_plot_data) > 0) {
  timeline_plot_data <- timeline_plot_data %>%
    mutate(Participant_id_plot = factor(Participant_id, levels = sort(unique(Participant_id))))
  p_timeline <- ggplot(timeline_plot_data,
                       aes(x = time_order, y = Participant_id_plot)) +
    geom_line(aes(group = Participant_id), colour = "grey70", linewidth = 0.4) +
    geom_point(aes(fill = Infection_Status, shape = has_wgs_vf), size = 3, colour = "grey20") +
    scale_shape_manual(values = c("TRUE" = 21, "FALSE" = 4)) +
    labs(
      title = "Clinical timelines for residents with Not_UTI-to-UTI transitions",
      subtitle = "Point shape indicates whether the clinical episode has linked WGS/VF data",
      x = "Ordered time within participant",
      y = "Participant",
      fill = "Primary UTI status",
      shape = "WGS/VF linked",
      caption = sprintf(
        "Data: %s, %s, and %s. Denominator: %d clinical timeline rows from %d Not_UTI-to-UTI transition participants. Level of analysis: clinical episode timeline with WGS/VF availability overlay. Uricult ordering uses Collection_Date where available; display/fallback ordering is reported in results/vf/vf_transition_case_index.csv.",
        FILE_STATUS_MAP_POSTER, FILE_STATUS_MAP, FILE_VF_READY,
        nrow(timeline_plot_data), n_distinct(timeline_plot_data$Participant_id)
      )
    ) +
    plot_theme_vf(base_size = 10)
  ggsave(file.path(DIR_PLOTS_VF, "vf_transition_case_timeline.png"),
         p_timeline, width = 11, height = max(5, length(timeline_participants) * 0.25), dpi = 300)
}

if (nrow(score_changes) > 0 && length(not_uti_uti_ids) > 0) {
  sc_plot <- score_changes %>%
    filter(case_id %in% not_uti_uti_ids) %>%
    filter(score_name %in% c("expec_marker_count", "upec_system_count",
                             "upec_system_fraction", "total_vf_count_curated"))
  if (nrow(sc_plot) > 0) {
    sc_long <- sc_plot %>%
      pivot_longer(cols = c(from_value, to_value), names_to = "endpoint", values_to = "value") %>%
      mutate(endpoint = ifelse(endpoint == "from_value", "Before", "After"))
    p_slope <- ggplot(sc_long, aes(x = endpoint, y = value, group = case_id, colour = case_id)) +
      geom_line(linewidth = 0.8) +
      geom_point(size = 2.5) +
      facet_wrap(~score_label, scales = "free_y") +
      labs(
        title = "Supplementary VF endpoint changes across Not_UTI-to-UTI transitions",
        subtitle = "Each line is one WGS/VF-linked clinical Not_UTI-to-UTI transition case",
        x = NULL,
        y = "Endpoint value",
        colour = "Case",
        caption = sprintf(
          "Data: %s, script 26 module outputs, and script 27 supplementary endpoint outputs. Denominator: %d endpoint-change rows from %d Not_UTI-to-UTI cases. Level of analysis: transition case-study pair. Descriptive only; no p-values are shown, and changes may reflect lineage/strain replacement, assembly differences, or true VF content change.",
          FILE_VF_READY, nrow(sc_plot), n_distinct(sc_plot$case_id)
        )
      ) +
      plot_theme_vf(base_size = 10) +
      theme(legend.position = "bottom")
    ggsave(file.path(DIR_PLOTS_VF, "vf_transition_score_slopeplot.png"),
           p_slope, width = 12, height = 8.5, dpi = 300)
  }
}

if (nrow(module_changes) > 0 && length(not_uti_uti_ids) > 0) {
  mc_plot <- module_changes %>%
    filter(case_id %in% not_uti_uti_ids, change_type != "stable_absent")
  if (nrow(mc_plot) > 0) {
    p_mod <- ggplot(mc_plot, aes(x = case_id, y = system_name, fill = change_type)) +
      geom_tile(colour = "white") +
      scale_fill_manual(values = c("gained" = "#c0392b", "lost" = "#2c7fb8",
                                   "stable_present" = "#41ab5d")) +
      labs(
        title = "VF module changes in Not_UTI-to-UTI transition cases",
        subtitle = "Modules are descriptive biological groupings, not validated causal virulence scores",
        x = "Case",
        y = "VF module",
        fill = "Change",
        caption = sprintf(
          "Data: %s and script 26 module definitions. Denominator: %d module-change rows from %d Not_UTI-to-UTI cases. Level of analysis: transition case-study pair. Interpret SNP same-strain context first; ST is secondary lineage context.",
          FILE_VF_READY, nrow(mc_plot), n_distinct(mc_plot$case_id)
        )
      ) +
      plot_theme_vf(base_size = 9) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    ggsave(file.path(DIR_PLOTS_VF, "vf_transition_module_change_heatmap.png"),
           p_mod, width = 11, height = max(6, n_distinct(mc_plot$system_name) * 0.25), dpi = 300)
  }
}

if (nrow(gene_changes) > 0 && length(not_uti_uti_ids) > 0) {
  gc_plot <- gene_changes %>%
    filter(case_id %in% not_uti_uti_ids, change_type %in% c("gained", "lost")) %>%
    group_by(Gene) %>%
    mutate(n_cases_changed = n_distinct(case_id)) %>%
    ungroup() %>%
    arrange(desc(biological_priority == "high"), desc(n_cases_changed), Gene) %>%
    slice_head(n = 50)
  if (nrow(gc_plot) > 0) {
    p_gene <- ggplot(gc_plot, aes(x = case_id, y = reorder(Gene, n_cases_changed), fill = change_type)) +
      geom_tile(colour = "white") +
      scale_fill_manual(values = c("gained" = "#c0392b", "lost" = "#2c7fb8")) +
      labs(
        title = "VF gene gains and losses in Not_UTI-to-UTI transition cases",
        subtitle = "Top changed genes are shown for interpretability; case-study evidence is descriptive",
        x = "Case",
        y = "VF gene",
        fill = "Change",
        caption = sprintf(
          "Data: %s and results/vf/vf_transition_gene_changes.csv. Denominator: %d changed gene rows from %d Not_UTI-to-UTI cases. Level of analysis: transition case-study pair. Gene gain/loss should be interpreted after SNP same-strain classification; ST is secondary lineage context.",
          FILE_VF_READY, nrow(gc_plot), n_distinct(gc_plot$case_id)
        )
      ) +
      plot_theme_vf(base_size = 9) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    ggsave(file.path(DIR_PLOTS_VF, "vf_transition_gene_gain_loss_tile.png"),
           p_gene, width = 10, height = max(5, n_distinct(gc_plot$Gene) * 0.22), dpi = 300)
  }
}

not_uti_uti_context <- not_uti_uti_summary %>%
  left_join(
    strain_ctx %>%
      select(case_id, VF_Jaccard_pairwise, replacement_flag, strain_context_note,
             snp_context_from_strain = snp_strain_context,
             st_context_from_strain = st_lineage_context,
             pair_interpretation_from_context = pair_interpretation),
    by = "case_id"
  ) %>%
  mutate(
    snp_context_plot = coalesce(snp_context_from_strain, snp_strain_context,
                                "Missing SNP evidence"),
    st_context_plot = coalesce(st_context_from_strain, st_lineage_context,
                               "Missing ST evidence"),
    pair_interpretation_plot = coalesce(pair_interpretation_from_context, pair_interpretation,
                                        "Missing strain metrics"),
    case_class_plot = case_when(
      !has_vf_pair ~ "Missing WGS/VF endpoint",
      is.na(case_class) | case_class == "" ~ "Unclassified WGS/VF pair",
      TRUE ~ case_class
    ),
    total_vf_gene_changes = n_vf_genes_gained + n_vf_genes_lost,
    snp_distance_plus_one = SNPs + 1
  )

strain_change_plot <- not_uti_uti_context %>%
  filter(has_vf_pair, !is.na(SNPs), !is.na(delta_total_vf_burden))

if (nrow(strain_change_plot) > 0) {
  p_strain_ctx <- ggplot(strain_change_plot,
                         aes(x = snp_distance_plus_one,
                             y = delta_total_vf_burden,
                             colour = snp_context_plot,
                             shape = st_context_plot,
                             size = total_vf_gene_changes)) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey55") +
    geom_vline(xintercept = SNP_THRESHOLD + 1, linetype = "dotted", colour = "grey40") +
    geom_point(alpha = 0.82) +
    scale_x_log10(labels = scales::comma) +
    scale_size_continuous(range = c(2.2, 7), breaks = scales::pretty_breaks(n = 4)) +
    labs(
      title = "Not_UTI-to-UTI VF changes require strain-context interpretation",
      subtitle = "Colour shows SNP same-strain context; shape shows secondary ST lineage context",
      x = "SNP distance + 1 (log scale)",
      y = "Change in total VF burden (UTI - prior Not_UTI)",
      colour = "SNP context",
      shape = "ST context",
      size = "VF genes gained + lost",
      caption = sprintf(
        "Data: %s and %s. Denominator: %d WGS/VF-linked Not_UTI-to-UTI transition cases with SNP distance and VF burden change. Level of analysis: transition case-study pair. Dashed horizontal line indicates no total burden change; dotted vertical line marks the SNP same-strain threshold (%d). ST is secondary lineage context and does not prove same strain.",
        file.path(DIR_VF, "vf_transition_case_summary.csv"),
        file.path(DIR_VF, "vf_transition_strain_context.csv"),
        nrow(strain_change_plot), SNP_THRESHOLD
      )
    ) +
    plot_theme_vf(base_size = 10)

  strain_context_plot <- file.path(DIR_PLOTS_VF, "vf_not_uti_uti_transition_strain_context.png")
  ggsave(strain_context_plot, p_strain_ctx, width = 9.5, height = 6.5, dpi = 300)
}

case_class_counts <- not_uti_uti_context %>%
  count(case_class_plot, has_vf_pair, name = "n") %>%
  mutate(
    endpoint_status = ifelse(has_vf_pair, "WGS/VF pair available", "Missing WGS/VF endpoint"),
    case_class_plot = reorder(case_class_plot, n)
  )

if (nrow(case_class_counts) > 0) {
  p_case_classes <- ggplot(case_class_counts,
                           aes(x = case_class_plot, y = n, fill = endpoint_status)) +
    geom_col(width = 0.68, colour = "white", linewidth = 0.25) +
    geom_text(aes(label = n), hjust = -0.15, size = 3.3) +
    coord_flip() +
    scale_y_continuous(expand = expansion(mult = c(0, 0.14))) +
    scale_fill_manual(values = c("WGS/VF pair available" = "#0072B2",
                                 "Missing WGS/VF endpoint" = "grey65")) +
    labs(
      title = "Not_UTI-to-UTI transition case classes and missing genomic endpoints",
      subtitle = "Class counts keep missing WGS/VF pairs visible in the transition denominator",
      x = NULL,
      y = "Number of Not_UTI-to-UTI transitions",
      fill = "Endpoint status",
      caption = sprintf(
        "Data: %s. Denominator: %d clinical Not_UTI-to-UTI transitions, including %d with WGS/VF-linked endpoints. Level of analysis: transition case classification. Case classes are descriptive and should be interpreted with SNP/ST evidence; they do not test whether VF changes cause UTI.",
        file.path(DIR_VF, "vf_transition_case_summary.csv"),
        nrow(not_uti_uti_context), sum(not_uti_uti_context$has_vf_pair, na.rm = TRUE)
      )
    ) +
    plot_theme_vf(base_size = 10)

  case_class_plot <- file.path(DIR_PLOTS_VF, "vf_not_uti_uti_transition_case_classes.png")
  ggsave(case_class_plot, p_case_classes,
         width = 9, height = max(5.2, nrow(case_class_counts) * 0.34), dpi = 300)
}

snp_plot <- case_summary %>%
  filter(is_not_uti_to_uti, has_vf_pair, !is.na(SNPs), !is.na(vf_jaccard))
if (nrow(snp_plot) > 0) {
  p_snp <- ggplot(snp_plot, aes(x = SNPs, y = vf_jaccard, colour = case_class, label = case_id)) +
    geom_vline(xintercept = SNP_THRESHOLD, linetype = "dashed", colour = "grey50") +
    geom_point(size = 3) +
    geom_text(nudge_y = 0.02, size = 3, show.legend = FALSE) +
    labs(
      title = "SNP distance versus VF similarity in Not_UTI-to-UTI transition cases",
      subtitle = sprintf("Dashed line marks SNP threshold %d used by strain interpretation", SNP_THRESHOLD),
      x = "SNP distance",
      y = "VF Jaccard similarity",
      colour = "Case class",
      caption = sprintf(
        "Data: %s plus pairwise strain metrics. Denominator: %d Not_UTI-to-UTI cases with WGS/VF pair, SNP distance, and VF Jaccard. Level of analysis: transition case-study pair. Low SNP distance plus high VF similarity supports persistence but does not establish VF causality for UTI.",
        FILE_VF_READY, nrow(snp_plot)
      )
    ) +
    plot_theme_vf(base_size = 10)
  ggsave(file.path(DIR_PLOTS_VF, "vf_transition_snp_vs_vf_jaccard.png"),
         p_snp, width = 9, height = 6.5, dpi = 300)
}

append_denominator_summary(
  case_index,
  "28_vf_transition_case_studies.R",
  "clinical_transition_index",
  "clinical_episode_transition",
  file.path(DIR_VF, "vf_transition_case_index.csv"),
  "Clinical-first transition index; Collection_Date is preferred for Uricult ordering and poster labels remain caveats"
)
write_uti_attrition_outputs()

msg("✓ 28_vf_transition_case_studies.R complete.")

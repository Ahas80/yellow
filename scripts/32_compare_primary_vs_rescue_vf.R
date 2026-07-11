#!/usr/bin/env Rscript
# ==============================================================================
# 32_compare_primary_vs_rescue_vf.R
# ------------------------------------------------------------------------------
# Compare the canonical VF-ready dataset against the optional bad-size assembly
# rescue sensitivity dataset.
# ==============================================================================

source("00_config.R")
source("R/pipeline_qc_helpers.R")

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
  library(stringr)
  library(optparse)
})

option_list <- list(
  make_option(c("--primary"),
    type = "character", default = FILE_VF_READY,
    help = "Canonical VF-ready dataset [default: %default]"
  ),
  make_option(c("--rescue"),
    type = "character", default = file.path(DIR_RESULTS, "sensitivity", "vf_rescue", "vf_analysis_ready_rescue.csv"),
    help = "Rescue sensitivity VF-ready dataset [default: %default]"
  ),
  make_option(c("--selection_file"),
    type = "character", default = file.path(DIR_RESULTS, "sensitivity", "rescue_qc", "rescue_assembly_selection.csv"),
    help = "Rescue assembly selection file [default: %default]"
  ),
  make_option(c("--out_dir"),
    type = "character", default = file.path(DIR_RESULTS, "sensitivity", "vf_rescue"),
    help = "Output directory [default: %default]"
  )
)
opt <- parse_args(OptionParser(option_list = option_list))

ensure_dir(opt$out_dir)

if (!file.exists(opt$primary)) stop("Missing primary VF-ready dataset: ", opt$primary)
if (!file.exists(opt$rescue)) stop("Missing rescue VF-ready dataset: ", opt$rescue)

as_bool <- function(x) {
  if (is.logical(x)) return(x)
  x <- trimws(tolower(as.character(x)))
  dplyr::case_when(
    x %in% c("true", "t", "1", "yes", "y") ~ TRUE,
    x %in% c("false", "f", "0", "no", "n") ~ FALSE,
    TRUE ~ NA
  )
}

parse_num <- function(x) {
  suppressWarnings(as.numeric(as.character(x)))
}

normalise_ready <- function(df) {
  df %>%
    prefer_primary_uti_status() %>%
    mutate(
      Participant_id = as.character(Participant_id),
      tp_lab = normalise_timepoint_preserve_events(tp_lab),
      Infection_Status = as.character(UTI_Status),
      Episode_ID = if ("Episode_ID" %in% names(.)) as.character(Episode_ID) else NA_character_
    )
}

primary <- read_csv(opt$primary, show_col_types = FALSE) %>% normalise_ready()
rescue <- read_csv(opt$rescue, show_col_types = FALSE) %>% normalise_ready()

selection <- if (file.exists(opt$selection_file)) {
  read_csv(opt$selection_file, show_col_types = FALSE) %>%
    mutate(
      Participant_id = as.character(Participant_id),
      tp_lab = normalise_timepoint_preserve_events(tp_lab),
      rescue_accepted = if ("rescue_accepted" %in% names(.)) as_bool(rescue_accepted) else FALSE,
      selected_canonical = if ("selected_canonical" %in% names(.)) as_bool(selected_canonical) else FALSE
    )
} else {
  tibble()
}

status_counts <- bind_rows(
  primary %>% count(dataset = "primary", Infection_Status, name = "n"),
  rescue %>% count(dataset = "rescue", Infection_Status, name = "n")
) %>%
  arrange(Infection_Status, dataset)
write_csv(status_counts, file.path(opt$out_dir, "primary_vs_rescue_status_counts.csv"))

primary_keys <- primary %>%
  distinct(Participant_id, tp_lab, .keep_all = TRUE) %>%
  select(any_of(c("Participant_id", "tp_lab", "Episode_ID", "Infection_Status", "ST", "vf_count_total", "total_vf_count_all")))

rescue_keys <- rescue %>%
  distinct(Participant_id, tp_lab, .keep_all = TRUE) %>%
  select(any_of(c("Participant_id", "tp_lab", "Episode_ID", "Infection_Status", "ST", "vf_count_total", "total_vf_count_all")))

accepted_rescue <- if (nrow(selection) > 0) {
  selection %>%
    filter(rescue_accepted %in% TRUE) %>%
    select(any_of(c(
      "Participant_id", "tp_lab", "Assembly_ID", "file_name", "full_path",
      "original_full_path", "rescue_full_path", "rescue_size_class",
      "analysis_size_class", "rescue_decision", "rescue_reason",
      "checkm2_completeness", "checkm2_contamination", "gunc_pass",
      "mlst_scheme", "mlst_st", "mash_min_distance"
    )))
} else {
  tibble(Participant_id = character(), tp_lab = character())
}

rescue_only <- rescue_keys %>%
  anti_join(primary_keys %>% select(Participant_id, tp_lab), by = c("Participant_id", "tp_lab")) %>%
  left_join(accepted_rescue, by = c("Participant_id", "tp_lab"), relationship = "many-to-many") %>%
  arrange(Infection_Status, Participant_id, tp_lab)

primary_only <- primary_keys %>%
  anti_join(rescue_keys %>% select(Participant_id, tp_lab), by = c("Participant_id", "tp_lab")) %>%
  arrange(Infection_Status, Participant_id, tp_lab)

write_csv(rescue_only, file.path(opt$out_dir, "rescued_episode_table.csv"))
write_csv(primary_only, file.path(opt$out_dir, "primary_only_episode_table.csv"))

gene_cols <- canonical_vf_gene_cols(required = FALSE)
gene_cols <- intersect(gene_cols, intersect(names(primary), names(rescue)))
if (length(gene_cols) == 0) {
  metadata_cols <- c(
    "Participant_id", "tp_lab", "Episode_ID", "Event_type", "Collection_Date",
    "Infection_Status", "Batch", "Status_Confidence_epi", "Sx_source_epi",
    "UTI_Label", "Urine_collection_method", "uricult_bridge_applied", "ST",
    "vf_count_total", "total_vf_count_all", "total_vf_count_curated",
    "total_vf_count_upec_candidate", "total_vf_count_unassigned",
    "low_confidence_count", "n_timepoints"
  )
  candidates <- setdiff(intersect(names(primary), names(rescue)), metadata_cols)
  gene_cols <- candidates[vapply(candidates, function(nm) {
    vals <- unique(na.omit(c(primary[[nm]], rescue[[nm]])))
    length(vals) > 0 && all(vals %in% c(0, 1, "0", "1"))
  }, logical(1))]
}

add_total <- function(df, total_name) {
  if ("vf_count_total" %in% names(df)) {
    df[[total_name]] <- parse_num(df$vf_count_total)
  } else if (length(gene_cols) > 0) {
    df[[total_name]] <- rowSums(as.data.frame(lapply(df[gene_cols], parse_num)), na.rm = TRUE)
  } else {
    df[[total_name]] <- NA_real_
  }
  df
}

primary_comp <- primary %>%
  add_total("primary_vf_count_total") %>%
  select(Participant_id, tp_lab, primary_vf_count_total, all_of(gene_cols))

rescue_comp <- rescue %>%
  add_total("rescue_vf_count_total") %>%
  select(Participant_id, tp_lab, rescue_vf_count_total, all_of(gene_cols))

delta <- inner_join(
  primary_comp,
  rescue_comp,
  by = c("Participant_id", "tp_lab"),
  suffix = c(".primary", ".rescue")
) %>%
  mutate(vf_count_total_delta = rescue_vf_count_total - primary_vf_count_total)

if (length(gene_cols) > 0 && nrow(delta) > 0) {
  primary_mat <- as.matrix(as.data.frame(lapply(delta[paste0(gene_cols, ".primary")], parse_num)))
  rescue_mat <- as.matrix(as.data.frame(lapply(delta[paste0(gene_cols, ".rescue")], parse_num)))
  delta$genes_added <- rowSums(rescue_mat == 1 & primary_mat != 1, na.rm = TRUE)
  delta$genes_lost <- rowSums(primary_mat == 1 & rescue_mat != 1, na.rm = TRUE)
  delta$genes_changed_total <- delta$genes_added + delta$genes_lost
} else {
  delta$genes_added <- NA_integer_
  delta$genes_lost <- NA_integer_
  delta$genes_changed_total <- NA_integer_
}

delta_out <- delta %>%
  select(Participant_id, tp_lab, primary_vf_count_total, rescue_vf_count_total,
         vf_count_total_delta, genes_added, genes_lost, genes_changed_total) %>%
  arrange(desc(abs(vf_count_total_delta)), desc(genes_changed_total), Participant_id, tp_lab)
write_csv(delta_out, file.path(opt$out_dir, "vf_delta_by_episode.csv"))

accepted_count <- if (nrow(selection) > 0 && "rescue_accepted" %in% names(selection)) {
  sum(selection$rescue_accepted %in% TRUE, na.rm = TRUE)
} else {
  NA_integer_
}

dataset_summary <- tibble(
  dataset = c("primary", "rescue"),
  rows = c(nrow(primary), nrow(rescue)),
  participants = c(n_distinct(primary$Participant_id), n_distinct(rescue$Participant_id)),
  rows_with_status = c(sum(!is.na(primary$Infection_Status)), sum(!is.na(rescue$Infection_Status))),
  median_vf_count_total = c(median(parse_num(primary$vf_count_total), na.rm = TRUE), median(parse_num(rescue$vf_count_total), na.rm = TRUE))
)
write_csv(dataset_summary, file.path(opt$out_dir, "primary_vs_rescue_dataset_summary.csv"))

largest_deltas <- delta_out %>%
  filter(!is.na(vf_count_total_delta), vf_count_total_delta != 0 | genes_changed_total > 0) %>%
  slice_head(n = 20)

report <- c(
  "# Primary vs bad-size rescue VF sensitivity comparison",
  "",
  sprintf("Generated: %s", format(Sys.time())),
  "",
  "## Dataset summary",
  "",
  paste(capture.output(print(dataset_summary, n = Inf, width = Inf)), collapse = "\n"),
  "",
  "## Status counts",
  "",
  paste(capture.output(print(status_counts, n = Inf, width = Inf)), collapse = "\n"),
  "",
  "## Rescue additions",
  "",
  sprintf("- Accepted rescue assemblies in selection file: %s", accepted_count),
  sprintf("- Rescue-only participant-timepoint rows: %d", nrow(rescue_only)),
  sprintf("- Primary-only participant-timepoint rows: %d", nrow(primary_only)),
  "",
  if (nrow(rescue_only) > 0) {
    paste(capture.output(print(rescue_only, n = Inf, width = Inf)), collapse = "\n")
  } else {
    "No rescue-only rows were added to the VF-ready dataset."
  },
  "",
  "## Largest VF burden changes among shared rows",
  "",
  if (nrow(largest_deltas) > 0) {
    paste(capture.output(print(largest_deltas, n = Inf, width = Inf)), collapse = "\n")
  } else {
    "No VF burden changes were detected among shared participant-timepoint rows."
  },
  "",
  "## Output files",
  "",
  sprintf("- Dataset summary: `%s`", file.path(opt$out_dir, "primary_vs_rescue_dataset_summary.csv")),
  sprintf("- Status counts: `%s`", file.path(opt$out_dir, "primary_vs_rescue_status_counts.csv")),
  sprintf("- Rescue-only rows: `%s`", file.path(opt$out_dir, "rescued_episode_table.csv")),
  sprintf("- Primary-only rows: `%s`", file.path(opt$out_dir, "primary_only_episode_table.csv")),
  sprintf("- VF deltas: `%s`", file.path(opt$out_dir, "vf_delta_by_episode.csv"))
)
writeLines(report, file.path(opt$out_dir, "primary_vs_rescue_summary.md"))

msg("Wrote primary vs rescue summary: %s", file.path(opt$out_dir, "primary_vs_rescue_summary.md"))

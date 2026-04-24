#!/usr/bin/env Rscript
# ==============================================================================
# 22_vf_build_analysis_dataset.R
# ------------------------------------------------------------------------------
# Role: [VF Core Data] - Build canonical episode-level and transition-level
# analysis datasets for virulence-factor work.
#
# Why this script exists:
#   - Current VF logic is spread across 02_gene_presence_analysis.R,
#     04_gene_breakdown.R, 15_longitudinal_patterns.R, 16_within_host_evolution.R,
#     and MLST/clinical outputs.
#   - Downstream VF scripts should NOT repeatedly re-implement joins.
#   - This script creates one authoritative analysis-ready dataset for all later
#     VF statistics, modeling, and reporting.
#
# Inputs:
#   - 00_config.R
#   - results/clinical/status_map.csv
#   - assembly_metadata.csv
#   - results/vf/vf_pa_all.csv
#   - results/vf/annotated_gene_table.csv (preferred, optional)
#   - results/mlst/mlst_all.tsv (optional but recommended)
#
# Outputs:
#   - results/vf/vf_episode_dataset.csv
#   - results/vf/vf_gene_long_dataset.csv
#   - results/vf/vf_transition_dataset.csv
#   - results/vf/vf_dataset_qc_summary.txt
#
# Usage:
#   Rscript 22_vf_build_analysis_dataset.R
# ==============================================================================

source("00_config.R")

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
  library(stringr)
  library(purrr)
  library(forcats)
})

msg <- function(...) message(sprintf(...))

# ------------------------------------------------------------------------------
# 1. Helpers
# ------------------------------------------------------------------------------
canon_status <- function(x) {
  x <- as.character(x)
  x <- stringr::str_trim(x)
  dplyr::case_when(
    x %in% c("UTI", "uti") ~ "UTI",
    x %in% c("ASB", "asb") ~ "ASB",
    x %in% c("Negative", "negative", "NEGATIVE") ~ "Negative",
    TRUE ~ x
  )
}

canon_tp <- function(x) {
  x <- as.character(x)
  x <- str_trim(x)
  x <- str_replace_all(x, "\\s+", "")
  x_low <- str_to_lower(x)

  case_when(
    str_detect(x_low, "^t[0-9]+$") ~ str_to_upper(x_low),
    str_detect(x_low, "^timepoint[0-9]+$") ~ str_replace(str_to_upper(x_low), "TIMEPOINT", "T"),
    x_low %in% c("uricult", "uti", "suspecteduti") ~ "Uricult",
    TRUE ~ x
  )
}

parse_tp_num <- function(tp) {
  tp <- canon_tp(tp)
  out <- suppressWarnings(as.numeric(str_extract(tp, "[0-9]+")))
  ifelse(tp == "Uricult", NA_real_, out)
}

ordered_tp_factor <- function(x) {
  ux <- unique(canon_tp(x))
  routine <- ux[str_detect(ux, "^T[0-9]+$")]
  routine <- routine[order(parse_tp_num(routine))]
  lvls <- c(routine, setdiff(c("Uricult", "Unscheduled"), routine), setdiff(ux, c(routine, "Uricult", "Unscheduled")))
  factor(canon_tp(x), levels = unique(lvls), ordered = TRUE)
}

safe_read_csv <- function(path) {
  if (!file.exists(path)) stop("Missing required file: ", path)
  readr::read_csv(path, show_col_types = FALSE)
}

first_existing <- function(paths) {
  hit <- paths[file.exists(paths)][1]
  if (length(hit) == 0 || is.na(hit)) return(NA_character_)
  hit
}

coalesce_chr <- function(...) {
  args <- list(...)
  out <- args[[1]]
  for (i in seq_along(args)[-1]) out <- dplyr::coalesce(out, args[[i]])
  out
}

# ------------------------------------------------------------------------------
# 2. Locate files
# ------------------------------------------------------------------------------
status_file <- file.path(DIR_CLINICAL, "status_map.csv")
meta_file   <- "assembly_metadata.csv"

vf_pa_file <- first_existing(c(
  file.path(DIR_RESULTS, "vf", "vf_pa_all.csv"),
  file.path(DIR_RESULTS, "vf_pa_all.csv")
))
if (is.na(vf_pa_file)) stop("Could not find vf_pa_all.csv")

annot_file <- first_existing(c(
  file.path(DIR_RESULTS, "vf", "annotated_gene_table.csv"),
  file.path(DIR_RESULTS, "annotated_gene_table.csv")
))

mlst_file <- first_existing(c(
  file.path(DIR_RESULTS, "mlst", "mlst_all.tsv"),
  file.path(DIR_RESULTS, "mlst_with_meta.csv")
))

out_dir <- file.path(DIR_RESULTS, "vf")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------------------------
# 3. Load data
# ------------------------------------------------------------------------------
status <- safe_read_csv(status_file) %>%
  mutate(
    Participant_id = as.character(Participant_id),
    Timepoint = canon_tp(coalesce_chr(.data$Timepoint, .data$tp_lab)),
    Infection_Status = canon_status(Infection_Status)
  )

meta <- safe_read_csv(meta_file) %>%
  mutate(
    Participant_id = as.character(Participant_id),
    Timepoint = canon_tp(coalesce_chr(.data$Timepoint, .data$tp_lab))
  )

vf_pa <- safe_read_csv(vf_pa_file) %>%
  mutate(
    Participant_id = as.character(Participant_id),
    Timepoint = canon_tp(coalesce_chr(.data$Timepoint, .data$tp_lab))
  )

mlst <- NULL
if (!is.na(mlst_file)) {
  if (grepl("\\.tsv$", mlst_file, ignore.case = TRUE)) {
    mlst <- readr::read_tsv(mlst_file, show_col_types = FALSE)
  } else {
    mlst <- readr::read_csv(mlst_file, show_col_types = FALSE)
  }
  mlst <- mlst %>%
    mutate(
      Participant_id = as.character(coalesce_chr(.data$Participant_id, .data$participant_id)),
      Timepoint = canon_tp(coalesce_chr(.data$Timepoint, .data$tp_lab, .data$timepoint))
    )
}

# ------------------------------------------------------------------------------
# 4. Identify gene columns and construct episode-level VF matrix
# ------------------------------------------------------------------------------
id_candidates <- c(
  "Participant_id", "Timepoint", "tp_lab", "Infection_Status", "Status",
  "Isolate", "Isolate_ID", "Assembly", "sample_id", "Episode_ID",
  "ST", "Phylogroup", "Batch", "n_genes", "vf_burden", "VF_Burden"
)

gene_cols <- setdiff(names(vf_pa), intersect(names(vf_pa), id_candidates))

if (length(gene_cols) == 0) {
  stop("No gene columns detected in vf_pa_all.csv after excluding metadata columns.")
}

vf_pa <- vf_pa %>%
  mutate(across(all_of(gene_cols), ~ as.integer(.x > 0)))

vf_episode <- vf_pa %>%
  left_join(
    status %>%
      select(any_of(c("Participant_id", "Timepoint", "Infection_Status", "Confidence", "CFU_class", "Collection_method"))) %>%
      distinct(),
    by = c("Participant_id", "Timepoint")
  ) %>%
  left_join(
    meta %>%
      select(any_of(c("Participant_id", "Timepoint", "Isolate", "Assembler", "Assembly", "GenomeSize", "N50", "Contigs", "Batch"))) %>%
      distinct(),
    by = c("Participant_id", "Timepoint")
  )

if (!is.null(mlst)) {
  mlst_keep <- intersect(names(mlst), c("Participant_id", "Timepoint", "ST", "SequenceType", "Phylogroup", "ClermonTyping"))
  vf_episode <- vf_episode %>% left_join(distinct(mlst[, mlst_keep]), by = c("Participant_id", "Timepoint"))
}

vf_episode <- vf_episode %>%
  mutate(
    tp_order = ordered_tp_factor(Timepoint),
    routine_tp_num = parse_tp_num(Timepoint),
    is_uricult = Timepoint == "Uricult",
    VF_Burden = rowSums(across(all_of(gene_cols)), na.rm = TRUE),
    n_detected_vf = VF_Burden,
    n_missing_gene_calls = rowSums(is.na(across(all_of(gene_cols))))
  )

# ------------------------------------------------------------------------------
# 5. Add category counts if annotation is available
# ------------------------------------------------------------------------------
if (!is.na(annot_file)) {
  annot <- safe_read_csv(annot_file)
  gene_name_col <- names(annot)[tolower(names(annot)) %in% c("gene", "gene_name", "vf_gene")][1]
  cat_col <- names(annot)[tolower(names(annot)) %in% c("category", "functional_category", "gene_category")][1]

  if (!is.na(gene_name_col) && !is.na(cat_col)) {
    gene_to_cat <- annot %>%
      transmute(gene = .data[[gene_name_col]], category = .data[[cat_col]]) %>%
      distinct() %>%
      filter(gene %in% gene_cols)

    if (nrow(gene_to_cat) > 0) {
      cat_counts <- vf_episode %>%
        select(Participant_id, Timepoint, all_of(gene_cols)) %>%
        pivot_longer(cols = all_of(gene_cols), names_to = "gene", values_to = "present") %>%
        left_join(gene_to_cat, by = "gene") %>%
        mutate(category = coalesce(category, "Unassigned")) %>%
        group_by(Participant_id, Timepoint, category) %>%
        summarise(category_count = sum(present, na.rm = TRUE), .groups = "drop") %>%
        tidyr::pivot_wider(names_from = category, values_from = category_count, values_fill = 0, names_prefix = "cat_")

      vf_episode <- vf_episode %>% left_join(cat_counts, by = c("Participant_id", "Timepoint"))
    }
  }
}

# ------------------------------------------------------------------------------
# 6. Long gene-level dataset
# ------------------------------------------------------------------------------
long_gene <- vf_episode %>%
  select(any_of(c("Participant_id", "Timepoint", "Infection_Status", "VF_Burden", "ST", "SequenceType", "Phylogroup")), all_of(gene_cols)) %>%
  pivot_longer(cols = all_of(gene_cols), names_to = "gene", values_to = "present")

# ------------------------------------------------------------------------------
# 7. Consecutive transition dataset
# ------------------------------------------------------------------------------
transitions <- vf_episode %>%
  arrange(Participant_id, tp_order, Timepoint) %>%
  group_by(Participant_id) %>%
  mutate(next_Timepoint = lead(Timepoint),
         next_Status = lead(Infection_Status),
         next_VF_Burden = lead(VF_Burden)) %>%
  ungroup()

pair_index <- vf_episode %>%
  arrange(Participant_id, tp_order, Timepoint) %>%
  group_by(Participant_id) %>%
  mutate(row_id = row_number(), next_row_id = lead(row_number())) %>%
  ungroup() %>%
  filter(!is.na(next_row_id)) %>%
  select(Participant_id, row_id, next_row_id, Timepoint_A = Timepoint, Timepoint_B = next_Timepoint,
         Status_A = Infection_Status, Status_B = next_Status)

vf_A <- vf_episode %>%
  arrange(Participant_id, tp_order, Timepoint) %>%
  group_by(Participant_id) %>%
  mutate(row_id = row_number()) %>%
  ungroup() %>%
  select(Participant_id, row_id, all_of(gene_cols))

vf_B <- vf_A %>% rename_with(~ paste0(.x, "_B"), all_of(gene_cols))

transition_wide <- pair_index %>%
  left_join(vf_A, by = c("Participant_id", "row_id")) %>%
  left_join(vf_B, by = c("Participant_id", "next_row_id" = "row_id"))

calc_transition_metrics <- function(df, genes) {
  A <- as.matrix(df[, genes, drop = FALSE])
  B <- as.matrix(df[, paste0(genes, "_B"), drop = FALSE])

  shared <- rowSums((A == 1) & (B == 1), na.rm = TRUE)
  gained <- rowSums((A == 0) & (B == 1), na.rm = TRUE)
  lost   <- rowSums((A == 1) & (B == 0), na.rm = TRUE)
  unionn <- rowSums((A == 1) | (B == 1), na.rm = TRUE)
  jacc   <- ifelse(unionn == 0, NA_real_, shared / unionn)

  tibble(shared_genes = shared, gained_genes = gained, lost_genes = lost, jaccard_similarity = jacc)
}

transition_metrics <- calc_transition_metrics(transition_wide, gene_cols)

transition_dataset <- bind_cols(
  transition_wide %>% select(Participant_id, Timepoint_A, Timepoint_B, Status_A, Status_B),
  transition_metrics
) %>%
  mutate(
    transition_type = paste0(Status_A, "->", Status_B),
    any_vf_change = (gained_genes + lost_genes) > 0,
    zero_vf_change = !any_vf_change
  )

# ------------------------------------------------------------------------------
# 8. Write outputs
# ------------------------------------------------------------------------------
write_csv(vf_episode, file.path(out_dir, "vf_episode_dataset.csv"))
write_csv(long_gene, file.path(out_dir, "vf_gene_long_dataset.csv"))
write_csv(transition_dataset, file.path(out_dir, "vf_transition_dataset.csv"))

qc_lines <- c(
  sprintf("VF episode rows: %d", nrow(vf_episode)),
  sprintf("VF gene columns: %d", length(gene_cols)),
  sprintf("Participants: %d", dplyr::n_distinct(vf_episode$Participant_id)),
  sprintf("Transitions: %d", nrow(transition_dataset)),
  sprintf("Statuses: %s", paste(sort(unique(na.omit(vf_episode$Infection_Status))), collapse = ", ")),
  sprintf("UTI-at-Uricult rows: %d", sum(vf_episode$Infection_Status == "UTI" & vf_episode$Timepoint == "Uricult", na.rm = TRUE))
)
writeLines(qc_lines, file.path(out_dir, "vf_dataset_qc_summary.txt"))

msg("Wrote:")
msg("  - %s", file.path(out_dir, "vf_episode_dataset.csv"))
msg("  - %s", file.path(out_dir, "vf_gene_long_dataset.csv"))
msg("  - %s", file.path(out_dir, "vf_transition_dataset.csv"))
msg("Done.")

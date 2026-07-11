#!/usr/bin/env Rscript
# ==============================================================================
# 29_vf_amr_combined_profile.R
# ==============================================================================
#
# GOAL:
#   Combine VF profiles with AMR data (if available) and plasmid replicon data.
#   If no dedicated AMR database screening exists, produce a VF + plasmid
#   exploratory analysis and a clear availability report.
#
# METHOD:
#   1. Audit whether true AMR data exist (ResFinder, CARD, AMRFinder, etc.).
#   2. Load supplementary VF endpoints from script 27 and plasmid replicon data.
#   3. Build episode-level combined VF + plasmid profiles.
#   4. Summarise replicon diversity by ST and primary UTI status.
#   5. Explore VF-plasmid co-occurrence patterns.
#
# OUTPUT:
#   - results/vf_amr/vf_amr_input_availability_report.txt
#   - results/vf_amr/vf_amr_combined_profile_table.csv
#   - results/vf_amr/vf_plasmid_combined_profile.csv
#   - results/vf_amr/vf_amr_score_summary_by_status.csv
#   - results/vf_amr/vf_amr_score_summary_by_ST.csv
#   - results/vf_amr/replicon_summary_by_status.csv
#   - results/vf_amr/replicon_summary_by_ST.csv
#   - results/vf_amr/vf_plasmid_correlation.csv
#   - plots/vf_amr/replicon_burden_by_status.png
#   - plots/vf_amr/vf_vs_replicon_scatter.png
#   - plots/vf_amr/replicon_heatmap_top_STs.png
#   - plots/vf_amr/vf_plasmid_analysis_scope.png
#
# NOTE:
#   This script does NOT invent AMR data. If true AMR results are absent,
#   only VF + plasmid replicon analysis is performed, and the report clearly
#   states that true VF+AMR integration was not done.
# ==============================================================================

source("00_config.R")
source("R/plot_helpers.R")
suppressPackageStartupMessages({
  library(dplyr); library(readr); library(tidyr); library(stringr)
  library(ggplot2); library(tibble)
})

msg("Starting 29_vf_amr_combined_profile.R")

# ==============================================================================
# 1. SETUP DIRECTORIES
# ==============================================================================
DIR_VF_AMR <- file.path(DIR_RESULTS, "vf_amr")
DIR_PLOTS_VF_AMR <- file.path(DIR_PLOTS, "vf_amr")
ensure_dir(DIR_VF_AMR)
ensure_dir(DIR_PLOTS_VF_AMR)

plot_theme_vf_amr <- function(base_size = 11) {
  theme_bw(base_size = base_size) +
    theme(
      plot.caption = element_text(hjust = 0, size = base_size - 3, colour = "grey35"),
      plot.subtitle = element_text(colour = "grey25"),
      legend.position = "bottom"
    )
}

normalise_tp_label <- function(x) {
  x <- as.character(x)
  x <- str_trim(x)
  case_when(
    str_detect(str_to_lower(x), "uricult") ~ "Uricult",
    str_detect(str_to_upper(x), "^T\\d+$") ~ str_to_upper(x),
    str_detect(x, "^\\d+$") ~ paste0("T", x),
    TRUE ~ x
  )
}

normalise_assembly_key <- function(x) {
  x <- basename(as.character(x))
  str_remove(x, "\\.(fasta|fa|fna)(\\.gz)?$")
}

extract_lab_isolate_id <- function(x) {
  str_extract(as.character(x), "\\d{4}C\\d{4,}(?:-\\d+)?|\\d{8,}(?:-\\d+)?")
}

# ==============================================================================
# 2. AUDIT AMR DATA AVAILABILITY
# ==============================================================================
report <- character()
ra <- function(...) report <<- c(report, sprintf(...))

ra("=== VF / AMR INPUT AVAILABILITY REPORT ===")
ra("Timestamp: %s", format(Sys.time()))
ra("")

# Check for dedicated AMR database results
amr_candidates <- c(
  file.path(DIR_RESULTS, "amr"),
  file.path(DIR_RESULTS, "abricate", "resfinder"),
  file.path(DIR_RESULTS, "abricate", "card"),
  file.path(DIR_RESULTS, "abricate", "ncbi"),
  file.path(DIR_RESULTS, "abricate", "amrfinder"),
  file.path(DIR_RESULTS, "amrfinder")
)

amr_found <- FALSE
amr_files <- character()
for (path in amr_candidates) {
  if (dir.exists(path)) {
    csvs <- list.files(path, pattern = "\\.csv$|\\.tsv$|\\.tab$|\\.txt$", full.names = TRUE, recursive = TRUE)
    if (length(csvs) > 0) {
      amr_found <- TRUE
      amr_files <- c(amr_files, csvs)
      ra("TRUE AMR data found at: %s (%d files)", path, length(csvs))
    }
  }
}

if (!amr_found) {
  ra("NO dedicated AMR database screening results found.")
  ra("Checked: %s", paste(amr_candidates, collapse = ", "))
  ra("")
  ra("IMPORTANT: True VF + AMR integration was NOT performed.")
  ra("The gene_map.csv from 04_gene_breakdown.R does NOT contain an AMR category.")
  ra("All VFDB-derived genes are either assigned to VF categories or 'Unassigned'.")
  ra("")
  ra("This script will produce VF + plasmid replicon exploratory analysis only.")
}

# Check for plasmid data
f_plasmid_pa   <- file.path(DIR_PLASMIDS, "plasmidfinder_presence_absence.csv")
f_plasmid_long <- file.path(DIR_PLASMIDS, "plasmidfinder_hits_long.csv")
f_replicon_long <- file.path(DIR_MLST, "plasmid_replicons_long.csv")
f_replicon_wide <- file.path(DIR_MLST, "plasmid_replicons_wide.csv")

has_plasmid_pa   <- file.exists(f_plasmid_pa)
has_plasmid_long <- file.exists(f_plasmid_long)
has_replicon_long <- file.exists(f_replicon_long)
has_replicon_wide <- file.exists(f_replicon_wide)

ra("")
ra("--- Plasmid data availability ---")
ra("plasmidfinder_presence_absence.csv: %s (%s)", has_plasmid_pa, f_plasmid_pa)
ra("plasmidfinder_hits_long.csv: %s (%s)", has_plasmid_long, f_plasmid_long)
ra("plasmid_replicons_long.csv: %s (%s)", has_replicon_long, f_replicon_long)
ra("plasmid_replicons_wide.csv: %s (%s)", has_replicon_wide, f_replicon_wide)

# Check supplementary VF endpoints
f_scores <- file.path(DIR_VF, "vf_score_table.csv")
has_scores <- file.exists(f_scores)
ra("")
ra("--- Supplementary VF endpoint data ---")
ra("vf_score_table.csv: %s", has_scores)

# Check VF analysis ready
has_vf <- file.exists(FILE_VF_READY)
ra("vf_analysis_ready.csv: %s", has_vf)

writeLines(report, file.path(DIR_VF_AMR, "vf_amr_input_availability_report.txt"))
msg("Availability report written")

if (!has_vf) {
  msg("WARNING: vf_analysis_ready.csv not found. Cannot proceed.")
  msg("✓ 29_vf_amr_combined_profile.R complete (no VF data).")
  quit(save = "no", status = 0)
}

# ==============================================================================
# 3. LOAD VF DATA
# ==============================================================================
if (file.exists(FILE_VF_PA) && file.info(FILE_VF_READY)$mtime < file.info(FILE_VF_PA)$mtime) {
  ready_probe <- read_csv(FILE_VF_READY, show_col_types = FALSE, n_max = Inf)
  pa_probe <- read_csv(FILE_VF_PA, show_col_types = FALSE, n_max = Inf)
  ready_genes <- canonical_vf_gene_cols(names(ready_probe))
  pa_genes <- canonical_vf_gene_cols(names(pa_probe))
  if (nrow(ready_probe) == nrow(pa_probe) && setequal(ready_genes, pa_genes)) {
    msg("WARNING: vf_analysis_ready.csv timestamp is older than vf_pa_all.csv, but row count and gene set match; continuing.")
  } else {
    stop(sprintf(
      "vf_analysis_ready.csv is older than vf_pa_all.csv and row/gene content differs. Re-run 22_vf_build_analysis_dataset.R before script 29.\n  vf_analysis_ready.csv: %s\n  vf_pa_all.csv: %s",
      format(file.info(FILE_VF_READY)$mtime), format(file.info(FILE_VF_PA)$mtime)
    ))
  }
}

vf <- read_csv(FILE_VF_READY, show_col_types = FALSE) %>%
  prefer_primary_uti_status() %>%
  apply_manual_sample_curation(context = "29_vf_ready") %>%
  filter_primary_genomics() %>%
  mutate(Participant_id = as.character(Participant_id))

scores <- if (has_scores) {
  read_csv(f_scores, show_col_types = FALSE) %>%
    prefer_primary_uti_status() %>%
    apply_manual_sample_curation(context = "29_vf_scores") %>%
    filter_primary_genomics() %>%
    mutate(Participant_id = as.character(Participant_id))
} else NULL

msg("VF data: %d episodes", nrow(vf))

# ==============================================================================
# 4. LOAD AND PROCESS PLASMID DATA
# ==============================================================================
if (!has_replicon_long && !has_plasmid_long) {
  msg("No plasmid replicon data available. Writing minimal outputs.")
  writeLines(c(report, "", "No plasmid data available for VF+plasmid analysis."),
             file.path(DIR_VF_AMR, "vf_amr_input_availability_report.txt"))
  msg("✓ 29_vf_amr_combined_profile.R complete (no plasmid data).")
  quit(save = "no", status = 0)
}

# Load replicon long-format data
if (has_replicon_long) {
  rep_long <- read_csv(f_replicon_long, show_col_types = FALSE)
} else {
  rep_long <- read_csv(f_plasmid_long, show_col_types = FALSE) %>%
    rename(isolate_id = Isolate_ID, replicon = GENE)
}

msg("Replicon hits: %d rows, %d unique replicons", nrow(rep_long), n_distinct(rep_long$replicon))

# Parse isolate_id to extract Participant_id and timepoint
# Format: PR00XX_barcodeYY_ZZZZZZZZZZ-N_assembler
# We need to map back to Participant_id and tp_lab using assembly_metadata
meta_file <- file.path(DIR_ROOT, "assembly_metadata.csv")
if (file.exists(meta_file)) {
  asm_meta <- read_csv(meta_file, show_col_types = FALSE) %>%
    apply_manual_sample_curation(context = "29_assembly_metadata") %>%
    filter_primary_genomics() %>%
    mutate(
      Participant_id = as.character(Participant_id),
      tp_lab = if ("tp_lab" %in% names(.)) as.character(tp_lab) else normalise_tp_label(Timepoint),
      assembly_key = normalise_assembly_key(coalesce(file_name, basename(full_path))),
      lab_isolate_id = extract_lab_isolate_id(Isolate_ID)
    )

  asm_lookup <- asm_meta %>%
    select(assembly_key, lab_isolate_id, Participant_id, tp_lab) %>%
    filter(!is.na(assembly_key)) %>%
    distinct()

  rep_mapped <- rep_long %>%
    mutate(
      assembly_key = normalise_assembly_key(isolate_id),
      lab_isolate_id = extract_lab_isolate_id(isolate_id)
    ) %>%
    left_join(asm_lookup %>% select(assembly_key, Participant_id, tp_lab) %>% distinct(),
              by = "assembly_key")

  # Fill any remaining misses by the laboratory isolate ID. This supports both
  # older numeric IDs and Batch 4-6 IDs containing a "C" token.
  if (any(is.na(rep_mapped$Participant_id))) {
    lab_lookup <- asm_lookup %>%
      filter(!is.na(lab_isolate_id)) %>%
      distinct(lab_isolate_id, Participant_id, tp_lab)
    rep_mapped <- rep_mapped %>%
      left_join(lab_lookup, by = "lab_isolate_id", suffix = c("", ".lab")) %>%
      mutate(
        Participant_id = coalesce(Participant_id, Participant_id.lab),
        tp_lab = coalesce(tp_lab, tp_lab.lab)
      ) %>%
      select(-any_of(c("Participant_id.lab", "tp_lab.lab")))
  }
} else {
  rep_mapped <- rep_long %>% mutate(Participant_id = NA_character_, tp_lab = NA_character_)
}

# If direct mapping failed, try extracting from isolate_id pattern
if (all(is.na(rep_mapped$Participant_id))) {
  msg("Attempting to extract Participant_id from isolate_id naming convention...")
  # Try to use the VF PA matrix's sample keys or pairwise_metrics for mapping
  # Load pairwise_metrics which has both Isolate_ID and Participant_id
  f_pair <- file.path(DIR_STRAIN, "pairwise_metrics.csv")
  if (file.exists(f_pair)) {
    pw <- read_csv(f_pair, show_col_types = FALSE)
    id_map <- bind_rows(
      pw %>% select(Isolate_ID = Isolate_ID_A, Participant_id = Participant_id_A,
                     Timepoint = Timepoint_A) %>% distinct(),
      pw %>% select(Isolate_ID = Isolate_ID_B, Participant_id = Participant_id_B,
                     Timepoint = Timepoint_B) %>% distinct()
    ) %>%
      distinct() %>%
      mutate(Participant_id = as.character(Participant_id),
             tp_lab = as.character(Timepoint))

    # Match: isolate_id in rep_long often starts with isolate ID
    rep_mapped <- rep_long %>%
      mutate(match_id = extract_lab_isolate_id(isolate_id)) %>%
      left_join(id_map %>% mutate(match_id = extract_lab_isolate_id(Isolate_ID)) %>%
                  distinct(match_id, .keep_all = TRUE),
                by = "match_id") %>%
      select(-match_id)
  }
}

# Summarise per participant-timepoint: count replicons, list unique types
rep_episode <- rep_mapped %>%
  filter(!is.na(Participant_id), !is.na(tp_lab)) %>%
  group_by(Participant_id, tp_lab) %>%
  summarise(
    n_replicons = n_distinct(replicon),
    replicon_list = paste(sort(unique(replicon)), collapse = "; "),
    has_IncF = any(grepl("IncF", replicon, ignore.case = TRUE)),
    has_Col = any(grepl("^Col", replicon)),
    has_IncI = any(grepl("IncI", replicon, ignore.case = TRUE)),
    has_IncB = any(grepl("IncB", replicon, ignore.case = TRUE)),
    has_IncX = any(grepl("IncX", replicon, ignore.case = TRUE)),
    .groups = "drop"
  ) %>%
  mutate(Participant_id = as.character(Participant_id))

msg("Replicon episode summary: %d episodes mapped", nrow(rep_episode))

# ==============================================================================
# 5. BUILD COMBINED VF + PLASMID PROFILE
# ==============================================================================
combined <- vf %>%
  select(Participant_id, tp_lab, Infection_Status, Batch, ST, vf_count_total)

if (has_scores && !is.null(scores)) {
  score_cols <- intersect(c("total_vf_count_curated", "expec_marker_count",
                            "upec_system_count", "upec_system_fraction",
                            "total_vf_count_unassigned", "low_confidence_count"),
                          names(scores))
  combined <- combined %>%
    left_join(scores %>% select(Participant_id, tp_lab, all_of(score_cols)),
              by = c("Participant_id", "tp_lab"))
}

combined <- combined %>%
  left_join(rep_episode, by = c("Participant_id", "tp_lab")) %>%
  mutate(
    n_replicons = replace_na(n_replicons, 0L),
    replicon_list = replace_na(replicon_list, ""),
    across(any_of(c("has_IncF", "has_Col", "has_IncI", "has_IncB", "has_IncX")),
           ~replace_na(as.logical(.x), FALSE)),
    amr_data_available = FALSE,
    true_amr_integration_performed = FALSE,
    plasmid_data_available = TRUE,
    amr_gene_count_total = NA_integer_,
    amr_class_count = NA_integer_,
    high_vf_flag = if ("vf_count_curated" %in% names(.)) {
      vf_count_curated >= median(vf_count_curated, na.rm = TRUE)
    } else {
      vf_count_total >= median(vf_count_total, na.rm = TRUE)
    },
    high_amr_flag = NA,
    vf_amr_profile_group = "AMR unavailable; VF+plasmid only"
  )

write_csv(combined, file.path(DIR_VF_AMR, "vf_plasmid_combined_profile.csv"))
write_csv(combined, file.path(DIR_VF_AMR, "vf_amr_combined_profile_table.csv"))
msg("Combined profile: %d rows", nrow(combined))

# ==============================================================================
# 6. REPLICON SUMMARIES
# ==============================================================================
# By status
rep_by_status <- combined %>%
  filter(!is.na(Infection_Status)) %>%
  group_by(Infection_Status) %>%
  summarise(
    n_episodes = n(),
    median_replicons = median(n_replicons, na.rm = TRUE),
    mean_replicons = round(mean(n_replicons, na.rm = TRUE), 2),
    q25 = quantile(n_replicons, 0.25, na.rm = TRUE),
    q75 = quantile(n_replicons, 0.75, na.rm = TRUE),
    pct_with_IncF = round(100 * mean(has_IncF, na.rm = TRUE), 1),
    pct_with_Col = round(100 * mean(has_Col, na.rm = TRUE), 1),
    .groups = "drop"
  )
write_csv(rep_by_status, file.path(DIR_VF_AMR, "replicon_summary_by_status.csv"))

vf_amr_by_status <- combined %>%
  filter(!is.na(Infection_Status)) %>%
  select(Infection_Status, any_of(c("vf_count_total", "total_vf_count_curated",
                                    "expec_marker_count", "upec_system_count",
                                    "upec_system_fraction", "total_vf_count_unassigned",
                                    "low_confidence_count",
                                    "n_replicons", "amr_gene_count_total", "amr_class_count"))) %>%
  pivot_longer(-Infection_Status, names_to = "metric", values_to = "value") %>%
  group_by(Infection_Status, metric) %>%
  summarise(
    n_episodes = sum(!is.na(value)),
    median = ifelse(all(is.na(value)), NA_real_, median(value, na.rm = TRUE)),
    q25 = ifelse(all(is.na(value)), NA_real_, quantile(value, 0.25, na.rm = TRUE)),
    q75 = ifelse(all(is.na(value)), NA_real_, quantile(value, 0.75, na.rm = TRUE)),
    mean = ifelse(all(is.na(value)), NA_real_, round(mean(value, na.rm = TRUE), 2)),
    true_amr_integration_performed = FALSE,
    note = ifelse(str_detect(first(metric), "^amr_"), "No dedicated AMR screening output detected", "VF/plasmid descriptive metric"),
    .groups = "drop"
  )
write_csv(vf_amr_by_status, file.path(DIR_VF_AMR, "vf_amr_score_summary_by_status.csv"))

# By ST
rep_by_st <- combined %>%
  filter(!is.na(ST)) %>%
  group_by(ST) %>%
  summarise(
    n_episodes = n(),
    median_replicons = median(n_replicons, na.rm = TRUE),
    mean_replicons = round(mean(n_replicons, na.rm = TRUE), 2),
    pct_with_IncF = round(100 * mean(has_IncF, na.rm = TRUE), 1),
    .groups = "drop"
  ) %>%
  arrange(desc(n_episodes))
write_csv(rep_by_st, file.path(DIR_VF_AMR, "replicon_summary_by_ST.csv"))

vf_amr_by_st <- combined %>%
  filter(!is.na(ST)) %>%
  select(ST, any_of(c("vf_count_total", "total_vf_count_curated",
                      "expec_marker_count", "upec_system_count",
                      "upec_system_fraction", "total_vf_count_unassigned",
                      "low_confidence_count",
                      "n_replicons", "amr_gene_count_total", "amr_class_count"))) %>%
  pivot_longer(-ST, names_to = "metric", values_to = "value") %>%
  group_by(ST, metric) %>%
  summarise(
    n_episodes = sum(!is.na(value)),
    median = ifelse(all(is.na(value)), NA_real_, median(value, na.rm = TRUE)),
    q25 = ifelse(all(is.na(value)), NA_real_, quantile(value, 0.25, na.rm = TRUE)),
    q75 = ifelse(all(is.na(value)), NA_real_, quantile(value, 0.75, na.rm = TRUE)),
    true_amr_integration_performed = FALSE,
    note = ifelse(str_detect(first(metric), "^amr_"), "No dedicated AMR screening output detected", "VF/plasmid descriptive metric"),
    .groups = "drop"
  ) %>%
  arrange(desc(n_episodes), ST)
write_csv(vf_amr_by_st, file.path(DIR_VF_AMR, "vf_amr_score_summary_by_ST.csv"))

profile_groups <- combined %>%
  count(vf_amr_profile_group, Infection_Status, name = "n_episodes") %>%
  mutate(
    true_amr_integration_performed = FALSE,
    note = "Profile grouping not performed because true AMR data are unavailable"
  )
write_csv(profile_groups, file.path(DIR_VF_AMR, "vf_amr_profile_groups.csv"))

# ==============================================================================
# 7. VF-PLASMID CORRELATION
# ==============================================================================
cor_results <- tibble()
vf_scores <- intersect(c("vf_count_total", "total_vf_count_curated",
                         "expec_marker_count", "upec_system_count",
                         "upec_system_fraction"), names(combined))
for (sc in vf_scores) {
  valid <- combined %>% filter(!is.na(.data[[sc]]), !is.na(n_replicons))
  if (nrow(valid) >= 10) {
    ct <- cor.test(valid[[sc]], valid$n_replicons, method = "spearman", exact = FALSE)
    cor_results <- bind_rows(cor_results, tibble(
      vf_score = sc, plasmid_metric = "n_replicons",
      method = "Spearman", rho = round(ct$estimate, 3),
      p_value = ct$p.value, n = nrow(valid)
    ))
  }
}
if (nrow(cor_results) > 0) {
  write_csv(cor_results, file.path(DIR_VF_AMR, "vf_plasmid_correlation.csv"))
}

# ==============================================================================
# 8. PLOTS
# ==============================================================================
n_uti <- sum(combined$Infection_Status == "UTI", na.rm = TRUE)
n_not_uti <- sum(combined$Infection_Status == "Not_UTI", na.rm = TRUE)

# Replicon burden by status
p1 <- ggplot(combined %>% filter(!is.na(Infection_Status)),
             aes(x = Infection_Status, y = n_replicons, fill = Infection_Status)) +
  geom_boxplot(alpha = 0.7) +
  scale_fill_uti_status() +
  labs(title = "Plasmid replicon diversity by primary UTI status",
       subtitle = sprintf("UTI n=%d, Not_UTI n=%d", n_uti, n_not_uti),
       x = NULL, y = "Number of distinct replicon types") +
  theme_minimal(base_size = 12) + theme(legend.position = "none")
ggsave(file.path(DIR_PLOTS_VF_AMR, "replicon_burden_by_status.png"), p1, width = 7, height = 5, dpi = 150)

# VF vs replicon scatter
if ("vf_count_total" %in% names(combined)) {
  p2 <- ggplot(combined %>% filter(!is.na(Infection_Status)),
               aes(x = vf_count_total, y = n_replicons, colour = Infection_Status)) +
    geom_point(alpha = 0.6, size = 2) +
    geom_smooth(method = "lm", se = FALSE, linetype = "dashed", linewidth = 0.5) +
    scale_colour_uti_status() +
    labs(title = "VF burden vs plasmid replicon diversity",
         x = "Total VF gene count", y = "Number of replicon types") +
    theme_minimal(base_size = 11) + theme(legend.position = "right")
  ggsave(file.path(DIR_PLOTS_VF_AMR, "vf_vs_replicon_scatter.png"), p2, width = 8, height = 6, dpi = 150)
}

# Replicon heatmap for top STs
top_sts <- combined %>% count(ST) %>% filter(n >= 5) %>% pull(ST)
if (length(top_sts) >= 3) {
  inc_types <- c("has_IncF","has_Col","has_IncI","has_IncB","has_IncX")
  inc_types <- intersect(inc_types, names(combined))
  if (length(inc_types) > 0) {
    heat_data <- combined %>%
      filter(ST %in% top_sts) %>%
      group_by(ST) %>%
      summarise(across(all_of(inc_types), ~round(100 * mean(., na.rm = TRUE), 1)),
                .groups = "drop") %>%
      pivot_longer(-ST, names_to = "replicon_type", values_to = "pct") %>%
      mutate(replicon_type = gsub("has_", "", replicon_type))

    p3 <- ggplot(heat_data, aes(x = replicon_type, y = ST, fill = pct)) +
      geom_tile(colour = "white") +
      geom_text(aes(label = sprintf("%.0f%%", pct)), size = 3) +
      scale_fill_gradient(low = "white", high = "steelblue") +
      labs(title = "Replicon type prevalence by ST (≥5 episodes)",
           x = "Replicon family", y = "ST", fill = "% episodes") +
      theme_minimal(base_size = 11)
    ggsave(file.path(DIR_PLOTS_VF_AMR, "replicon_heatmap_top_STs.png"), p3, width = 8, height = 6, dpi = 150)
  }
}

scope_counts <- tibble(
  metric = c(
    "VF-ready rows",
    "Rows with plasmid data flag",
    "Rows with true AMR data",
    "Rows with true VF+AMR integration"
  ),
  n = c(
    nrow(combined),
    sum(combined$plasmid_data_available %in% TRUE, na.rm = TRUE),
    sum(combined$amr_data_available %in% TRUE, na.rm = TRUE),
    sum(combined$true_amr_integration_performed %in% TRUE, na.rm = TRUE)
  ),
  interpretation = c(
    "Canonical VF/WGS-linked E. coli isolates",
    "Plasmid/mobile-context data available",
    "Dedicated AMR database output available",
    "True VF+AMR integration performed"
  )
) %>%
  mutate(metric = factor(metric, levels = metric))

p_scope <- ggplot(scope_counts, aes(x = metric, y = n, fill = interpretation)) +
  geom_col(width = 0.64, colour = "white", linewidth = 0.25) +
  geom_text(aes(label = n), vjust = -0.25, size = 3.8) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.14))) +
  scale_fill_manual(values = c(
    "Canonical VF/WGS-linked E. coli isolates" = "#0072B2",
    "Plasmid/mobile-context data available" = "#009E73",
    "Dedicated AMR database output available" = "#D55E00",
    "True VF+AMR integration performed" = "grey55"
  )) +
  labs(
    title = "Scope of VF, plasmid, and AMR data integration",
    subtitle = "This script provides VF+plasmid/mobile-context summaries; true AMR screening is absent unless AMR rows are non-zero",
    x = NULL,
    y = "Rows / isolates",
    fill = "Interpretation",
    caption = paste(
      sprintf("Data: %s and %s.",
              file.path(DIR_VF_AMR, "vf_amr_input_availability_report.txt"),
              file.path(DIR_VF_AMR, "vf_amr_combined_profile_table.csv")),
      sprintf("Denominator: n=%d VF/WGS-linked E. coli isolates from %d participants.",
              nrow(combined), n_distinct(combined$Participant_id)),
      "Level of analysis: input availability and analysis-scope diagnostic.",
      "Rows with AMR data or true VF+AMR integration equal zero when no dedicated AMR screening output is available.",
      "Do not interpret plasmid replicon or VFDB-derived summaries as true AMR analysis."
    )
  ) +
  plot_theme_vf_amr(base_size = 11) +
  theme(axis.text.x = element_text(angle = 20, hjust = 1))

ggsave(file.path(DIR_PLOTS_VF_AMR, "vf_plasmid_analysis_scope.png"),
       p_scope, width = 8.5, height = 5.4, dpi = 300)

msg("✓ 29_vf_amr_combined_profile.R complete.")

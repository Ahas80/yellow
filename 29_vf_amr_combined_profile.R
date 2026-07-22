#!/usr/bin/env Rscript
# ==============================================================================
# 29_vf_amr_combined_profile.R
# ==============================================================================
#
# GOAL:
#   Own the complete genomic-AMR layer for the exact selected Longcycler cohort,
#   then integrate validated episode profiles with VF and plasmid outputs.
#
# METHOD:
#   1. Validate 532 selected FASTAs and sequence-equivalent Prokka annotations.
#   2. Run SHA-bound ABRicate-ResFinder, AMRFinderPlus and ResFinder/PointFinder.
#   3. Harmonise calls, audit caller agreement and build episode profiles.
#   4. Analyse 371 adjacent pairs and the nine Not_UTI-to-UTI transitions.
#   5. Integrate genomic AMR with existing VF and plasmid profiles.
#
# OUTPUT:
#   - results/amr/provenance/{input,run,published_output}_manifest.csv
#   - results/amr/harmonized_determinants_long.csv
#   - results/amr/{episode,resident}_amr_profiles.csv
#   - results/amr/caller_concordance_discrepancies.csv
#   - results/amr/{gene,class,mutation}_prevalence_episode_resident.csv
#   - results/amr/resfinder_predicted_phenotypes_genomic_not_ast.csv
#   - results/amr/adjacent_pair_amr_profiles_371.csv
#   - results/amr/not_uti_to_uti_amr_profiles_9.csv
#   - results/amr/longitudinal_resident_bootstrap_inference.csv
#   - results/amr/validation_checks.csv
#   - results/plasmids/mob_suite/plasmid_gene_locations_long.csv
#   - results/plasmids/mob_suite/episode_mechanism_profiles.csv
#   - results/plasmids/mob_suite/adjacent_pair_plasmid_metrics_371.csv
#   - results/plasmids/mob_suite/not_uti_to_uti_plasmid_metrics_9.csv
#   - plots/amr/*.png (three supplementary AMR figures)
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
#   AMRFinderPlus is the primary profile caller. ResFinder/PointFinder is
#   complementary and ABRicate-ResFinder is a legacy comparison. Genomic
#   predictions are not phenotypic AST.
# ==============================================================================

source("00_config.R")
source("R/plot_helpers.R")
source("R/amr_helpers.R")
source("R/plasmid_mechanism_helpers.R")
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
DIR_AMR <- file.path(DIR_RESULTS, "amr")
DIR_PLOTS_AMR <- file.path(DIR_PLOTS, "amr")
ensure_dir(DIR_VF_AMR)
ensure_dir(DIR_PLOTS_VF_AMR)
ensure_dir(DIR_AMR)
ensure_dir(DIR_PLOTS_AMR)

msg("Running/reusing authoritative genomic-AMR calls for the exact 532-episode cohort")
amr_analysis <- run_genomic_amr_analysis(
  root = DIR_ROOT,
  output_root = DIR_AMR,
  plot_root = DIR_PLOTS_AMR
)
amr_profiles <- amr_analysis$profiles
msg(
  "Genomic AMR complete: %d episode profiles, %d adjacent pairs, %d focused transitions",
  nrow(amr_profiles), nrow(amr_analysis$transitions), nrow(amr_analysis$focused)
)

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
# 2. RECORD AMR DATA AVAILABILITY
# ==============================================================================
report <- character()
ra <- function(...) report <<- c(report, sprintf(...))

ra("=== VF / AMR INPUT AVAILABILITY REPORT ===")
ra("Timestamp: %s", format(Sys.time()))
ra("")

# Record only the authoritative script-29 results. CARD and other unrelated
# screens are deliberately outside this analysis contract.
amr_files <- c(
  file.path(DIR_AMR, "episode_amr_profiles.csv"),
  file.path(DIR_AMR, "harmonized_determinants_long.csv"),
  file.path(DIR_AMR, "adjacent_pair_amr_profiles_371.csv"),
  file.path(DIR_AMR, "not_uti_to_uti_amr_profiles_9.csv"),
  file.path(DIR_AMR, "validation_checks.csv"),
  file.path(DIR_AMR, "RUN_COMPLETE.txt")
)
if (!all(file.exists(amr_files))) {
  stop(
    "Authoritative Script-29 AMR outputs are incomplete: ",
    paste(basename(amr_files[!file.exists(amr_files)]), collapse = ", "),
    call. = FALSE
  )
}
amr_found <- TRUE
ra("Authoritative AMR output directory: %s (%d required products)",
   DIR_AMR, length(amr_files))

ra("")
ra("Authoritative genomic-AMR status: COMPLETE")
ra("Primary caller: AMRFinderPlus acquired genes and known resistance mutations")
ra("Complementary caller: ResFinder/PointFinder")
ra("Legacy comparison: ABRicate-ResFinder at >=80%% identity and coverage")
ra("Episode profiles: %d", nrow(amr_profiles))
ra("Adjacent-pair profiles: %d", nrow(amr_analysis$transitions))
ra("Focused Not_UTI-to-UTI profiles: %d", nrow(amr_analysis$focused))
ra("Interpretation: genomic prediction—not phenotypic AST")

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
  stop(
    "vf_analysis_ready.csv is required for Script 29 VF/AMR integration.",
    call. = FALSE
  )
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
  stop(
    "Plasmid replicon data are required for Script 29 VF/plasmid integration.",
    call. = FALSE
  )
}

# Load replicon long-format data
if (has_plasmid_long) {
  rep_long <- read_csv(f_plasmid_long, show_col_types = FALSE) %>%
    rename(isolate_id = Isolate_ID, replicon = GENE)
} else {
  rep_long <- read_csv(f_replicon_long, show_col_types = FALSE)
}

msg("Replicon hits: %d rows, %d unique replicons", nrow(rep_long), n_distinct(rep_long$replicon))

# Map plasmid results only by the exact selected FASTA path and content hash.
meta_file <- FILE_ANALYSIS_ASSEMBLY_MANIFEST
asm_meta <- load_analysis_assemblies(meta_file, require_files = TRUE) %>%
  mutate(
    Participant_id = as.character(Participant_id),
    tp_lab = as.character(tp_lab),
    fasta_path = normalizePath(full_path, winslash = "/", mustWork = TRUE),
    fasta_sha256 = vapply(fasta_path, digest::digest, character(1), algo = "sha256", file = TRUE)
  )
required_replicon_provenance <- c("fasta_path", "fasta_sha256")
missing_replicon_provenance <- setdiff(required_replicon_provenance, names(rep_long))
if (length(missing_replicon_provenance)) {
  stop(
    "Plasmid results lack exact FASTA provenance (", paste(missing_replicon_provenance, collapse = ", "),
    "). Rerun 08_core_vs_plasmid.R or 09_inc_plasmid_network.R from the current Longcycler manifest."
  )
}
rep_long <- rep_long %>%
  mutate(
    fasta_path = normalizePath(fasta_path, winslash = "/", mustWork = FALSE),
    fasta_sha256 = tolower(as.character(fasta_sha256))
  )
replicon_endpoint_lookup <- asm_meta %>%
  transmute(
    fasta_path,
    fasta_sha256 = tolower(as.character(fasta_sha256)),
    manifest_Participant_id = as.character(Participant_id),
    manifest_tp_lab = normalise_timepoint_preserve_events(tp_lab)
  ) %>%
  distinct()
if (anyDuplicated(replicon_endpoint_lookup[c("fasta_path", "fasta_sha256")])) {
  stop("Selected Longcycler manifest has duplicate FASTA path/hash endpoint keys.")
}
unmatched_replicon <- rep_long %>%
  anti_join(replicon_endpoint_lookup, by = c("fasta_path", "fasta_sha256"))
if (nrow(unmatched_replicon)) {
  stop("Plasmid results contain ", nrow(unmatched_replicon), " row(s) not tied to the current selected Longcycler FASTA path and SHA-256.")
}
rep_long$source_Participant_id <- if ("Participant_id" %in% names(rep_long)) {
  as.character(rep_long$Participant_id)
} else {
  NA_character_
}
rep_long$source_tp_lab <- if ("tp_lab" %in% names(rep_long)) {
  normalise_timepoint_preserve_events(rep_long$tp_lab)
} else {
  NA_character_
}
rep_mapped <- rep_long %>%
  select(-any_of(c("Participant_id", "tp_lab"))) %>%
  inner_join(
    replicon_endpoint_lookup,
    by = c("fasta_path", "fasta_sha256"),
    relationship = "many-to-one"
  ) %>%
  mutate(
    Participant_id = .data$manifest_Participant_id,
    tp_lab = .data$manifest_tp_lab
  )
source_key_disagreement <-
  (!is.na(rep_mapped$source_Participant_id) &
     rep_mapped$source_Participant_id != rep_mapped$Participant_id) |
  (!is.na(rep_mapped$source_tp_lab) &
     rep_mapped$source_tp_lab != rep_mapped$tp_lab)
if (any(source_key_disagreement)) {
  stop(
    "Plasmid result episode identifiers disagree with the exact selected Longcycler FASTA path/hash for ",
    sum(source_key_disagreement), " hit row(s)."
  )
}
rep_mapped <- rep_mapped %>%
  select(-manifest_Participant_id, -manifest_tp_lab,
         -source_Participant_id, -source_tp_lab)

if (nrow(rep_mapped) == 0) {
  msg("No plasmid-replicon rows matched the selected Longcycler assembly keys; no isolate-ID fallback is permitted.")
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
  left_join(
    amr_profiles %>%
      select(
        Participant_id, tp_lab, mdfA_detected,
        informative_acquired_genes,
        acquired_genes_sensitivity_including_mdfA,
        informative_acquired_classes, primary_known_mutations,
        amr_gene_count_informative, amr_gene_count_including_mdfA,
        amr_class_count, amr_mutation_count,
        any_informative_acquired_amr,
        any_acquired_amr_including_mdfA
      ),
    by = c("Participant_id", "tp_lab"),
    relationship = "one-to-one"
  ) %>%
  mutate(
    n_replicons = replace_na(n_replicons, 0L),
    replicon_list = replace_na(replicon_list, ""),
    across(any_of(c("has_IncF", "has_Col", "has_IncI", "has_IncB", "has_IncX")),
           ~replace_na(as.logical(.x), FALSE)),
    amr_data_available = TRUE,
    true_amr_integration_performed = TRUE,
    plasmid_data_available = TRUE,
    amr_gene_count_total = amr_gene_count_informative,
    high_vf_flag = if ("vf_count_curated" %in% names(.)) {
      vf_count_curated >= median(vf_count_curated, na.rm = TRUE)
    } else {
      vf_count_total >= median(vf_count_total, na.rm = TRUE)
    },
    high_amr_flag = amr_gene_count_informative >=
      median(amr_gene_count_informative, na.rm = TRUE),
    vf_amr_profile_group = case_when(
      high_vf_flag & high_amr_flag ~ "High VF / high informative AMR",
      high_vf_flag & !high_amr_flag ~ "High VF / lower informative AMR",
      !high_vf_flag & high_amr_flag ~ "Lower VF / high informative AMR",
      TRUE ~ "Lower VF / lower informative AMR"
    )
  )

if (nrow(combined) != 532L ||
    anyNA(combined$amr_gene_count_informative) ||
    !all(combined$amr_data_available) ||
    !all(combined$true_amr_integration_performed)) {
  stop("VF/plasmid integration did not retain all 532 validated AMR profiles.")
}

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
    true_amr_integration_performed = TRUE,
    note = ifelse(
      str_detect(first(metric), "^amr_"),
      "AMRFinderPlus informative acquired-gene profile; mdf(A) excluded from primary burden",
      "VF/plasmid descriptive metric"
    ),
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
    true_amr_integration_performed = TRUE,
    note = ifelse(
      str_detect(first(metric), "^amr_"),
      "AMRFinderPlus informative acquired-gene profile; mdf(A) excluded from primary burden",
      "VF/plasmid descriptive metric"
    ),
    .groups = "drop"
  ) %>%
  arrange(desc(n_episodes), ST)
write_csv(vf_amr_by_st, file.path(DIR_VF_AMR, "vf_amr_score_summary_by_ST.csv"))

profile_groups <- combined %>%
  count(vf_amr_profile_group, Infection_Status, name = "n_episodes") %>%
  mutate(
    true_amr_integration_performed = TRUE,
    note = "Descriptive VF/AMR groups; genomic determinants are not phenotypic AST"
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
  scale_x_discrete(labels = c(UTI = "UTI", Not_UTI = "Not UTI")) +
  theme_ruti_publication() + theme(legend.position = "none")
ggsave(file.path(DIR_PLOTS_VF_AMR, "replicon_burden_by_status.png"), p1,
       width = 7, height = 5, dpi = 300, bg = "white")

# VF vs replicon scatter
if ("vf_count_total" %in% names(combined)) {
  p2 <- ggplot(combined %>% filter(!is.na(Infection_Status)),
               aes(x = vf_count_total, y = n_replicons, colour = Infection_Status)) +
    geom_point(alpha = 0.6, size = 2) +
    geom_smooth(method = "lm", se = FALSE, linetype = "dashed", linewidth = 0.5) +
    scale_colour_uti_status() +
    labs(title = "VF burden vs plasmid replicon diversity",
         x = "Total VF gene count", y = "Number of replicon types") +
    theme_ruti_publication() + theme(legend.position = "right")
  ggsave(file.path(DIR_PLOTS_VF_AMR, "vf_vs_replicon_scatter.png"), p2,
         width = 8, height = 6, dpi = 300, bg = "white")
}

# Replicon heatmap for top STs
top_sts <- combined %>%
  filter(!is.na(ST), nzchar(as.character(ST))) %>%
  count(ST) %>%
  filter(n >= 5) %>%
  pull(ST)
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
      geom_text(aes(label = ifelse(is.na(pct), "Unavailable", sprintf("%.0f%%", pct))), size = 3) +
      scale_fill_viridis_c(option = "C", limits = c(0, 100), na.value = "#BDBDBD") +
      labs(
        title = "Replicon-family prevalence by sequence type",
        subtitle = "Sequence types represented by at least five typed episodes; missing/non-typable ST calls are excluded",
        x = "Replicon family", y = "Sequence type", fill = "Episodes (%)",
        caption = "Descriptive episode-level prevalence; repeated episodes from the same participant are not independent."
      ) +
      theme_ruti_publication()
    ggsave(file.path(DIR_PLOTS_VF_AMR, "replicon_heatmap_top_STs.png"), p3,
           width = 8, height = 6, dpi = 300, bg = "white")
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
    subtitle = "All 532 selected Longcycler episodes have validated genomic-AMR profiles",
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
      "AMRFinderPlus defines the primary genomic-AMR profile; ResFinder/PointFinder and ABRicate are complementary/legacy evidence.",
      "Genomic AMR determinants and predicted phenotypes are not measured susceptibility."
    )
  ) +
  plot_theme_vf_amr(base_size = 11) +
  theme(axis.text.x = element_text(angle = 20, hjust = 1))

ggsave(file.path(DIR_PLOTS_VF_AMR, "vf_plasmid_analysis_scope.png"),
       p_scope, width = 8.5, height = 5.4, dpi = 300)

if (tolower(Sys.getenv("AMR_ONLY", "0")) %in% c("1", "true", "yes")) {
  msg(
    "AMR_ONLY requested: genomic-AMR marker is complete; ",
    "predicted-plasmid localization is deferred until script 09b completes."
  )
} else {
  msg(
    paste0(
      "Localizing pinned CGE VF, genomic-AMR and PlasmidFinder calls ",
      "to MOB predicted contig assignments"
    )
  )
  plasmid_mechanism <- run_plasmid_gene_localization(
    root = DIR_ROOT,
    amr_analysis = amr_analysis,
    output_root = file.path(DIR_PLASMIDS, "mob_suite")
  )
  if (
    nrow(plasmid_mechanism$episodes) != 532L ||
      nrow(plasmid_mechanism$adjacent) != 371L ||
      nrow(plasmid_mechanism$focused) != 9L
  ) {
    stop("Script-29 predicted-plasmid mechanism denominator gate failed.")
  }
  msg(
    paste0(
      "Predicted-plasmid localization complete: %d calls, %d episodes, ",
      "%d adjacent pairs and %d focused transitions"
    ),
    nrow(plasmid_mechanism$locations),
    nrow(plasmid_mechanism$episodes),
    nrow(plasmid_mechanism$adjacent),
    nrow(plasmid_mechanism$focused)
  )
}

msg("✓ 29_vf_amr_combined_profile.R complete.")

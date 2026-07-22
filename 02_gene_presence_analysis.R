#!/usr/bin/env Rscript
# ==============================================================================
# 02_gene_presence_analysis.R
# ==============================================================================
#
# GOAL:
#   Profile all assembled isolates for virulence factor (VF) genes using
#   Abricate against the VFDB database, and build a binary presence/absence
#   matrix.  This is the core VF data generation step: all downstream VF
#   analysis (scripts 03, 04, 05, 14, 22–25) depends on its outputs.
#
# WHY THIS SCRIPT EXISTS:
#   Virulence factors are the bacterial genes that enable pathogenesis—
#   adhesins, toxins, iron acquisition systems, etc.  By profiling every
#   isolate’s VF repertoire, we can ask: do bacteria from UTI episodes carry
#   different VFs from those from Not_UTI episodes?
#
# KEY DESIGN DECISIONS:
#   - Uses one selected, QC-passing Longcycler assembly per
#     participant×timepoint, with no alternative assembly input.
#   - Minimum thresholds: identity ≥ 80%, coverage ≥ 80% (Abricate defaults).
#   - Results are cached per-isolate in results/vf/abricate/ so re-runs
#     skip already-processed assemblies.
#
# INPUTS:
#   - assembly_metadata.csv              (isolate-level metadata)
#   - data/assemblies/*.fasta            (long-read assemblies)
#
# OUTPUTS:
#   - results/vf/vf_hits_all.rds         (full Abricate hit table)
#   - results/vf/vf_pa_all.csv           (binary P/A matrix: participant×tp)
#   - results/vf/stats_gene_level.csv    (per-gene prevalence stats)
#   - results/vf/abricate/               (per-isolate cache)
#
# DOWNSTREAM:
#   -> 22_vf_build_analysis_dataset.R joins vf_pa_all.csv with primary UTI status
#   → 14_genotype_phenotype_model.R uses vf_hits_all.rds for GLMM testing
# ==============================================================================

# 1. Load Configuration & Libraries
source("00_config.R")
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(purrr)
  library(furrr)
  library(stringr)
  library(ggplot2)
  library(forcats)
  library(glue)
  library(optparse)
})

# 2. Parse Arguments
# ------------------------------------------------------------------------------
option_list <- list(
  make_option(c("--min_id"),
    type = "integer", default = 80,
    help = "Minimum identity percentage [default: %default]"
  ),
  make_option(c("--min_cov"),
    type = "integer", default = 80,
    help = "Minimum coverage percentage [default: %default]"
  ),
  make_option(c("--selection_file"),
    type = "character", default = FILE_ANALYSIS_ASSEMBLY_MANIFEST,
    help = "Assembly selection CSV containing full_path and a logical selection column [default: %default]"
  ),
  make_option(c("--selection_column"),
    type = "character", default = "selected_canonical",
    help = "Logical column in --selection_file used to choose assemblies [default: %default]"
  ),
  make_option(c("--out_dir"),
    type = "character", default = DIR_VF,
    help = "Output directory for VF hit/matrix files [default: %default]"
  ),
  make_option(c("--out_suffix"),
    type = "character", default = "",
    help = "Suffix appended to output basenames, e.g. rescue -> vf_pa_all_rescue.csv [default: none]"
  ),
  make_option(c("--abricate_dir"),
    type = "character", default = NA_character_,
    help = "Abricate cache directory [default: <out_dir>/abricate]"
  ),
  make_option(c("--skip_nucmer"),
    action = "store_true", default = FALSE,
    help = "Skip optional Nucmer trajectory comparisons [default: %default]"
  )
)

opt <- parse_args(OptionParser(option_list = option_list))

log_info <- function(...) message(format(Sys.time(), "[%H:%M:%S] "), ...)

ensure_dir(opt$out_dir)
suffix_clean <- str_replace_all(coalesce(opt$out_suffix, ""), "^_+", "")
suffix_part <- if (nzchar(suffix_clean)) paste0("_", suffix_clean) else ""
vf_hits_file <- if (identical(normalizePath(opt$out_dir, winslash = "/", mustWork = FALSE),
                              normalizePath(DIR_VF, winslash = "/", mustWork = FALSE)) &&
                    !nzchar(suffix_part)) FILE_VF_HITS else file.path(opt$out_dir, paste0("vf_hits_all", suffix_part, ".rds"))
vf_pa_file <- if (identical(normalizePath(opt$out_dir, winslash = "/", mustWork = FALSE),
                            normalizePath(DIR_VF, winslash = "/", mustWork = FALSE)) &&
                  !nzchar(suffix_part)) FILE_VF_PA else file.path(opt$out_dir, paste0("vf_pa_all", suffix_part, ".csv"))
stats_gene_file <- file.path(opt$out_dir, paste0("stats_gene_level", suffix_part, ".csv"))
plot_out_dir <- if (identical(normalizePath(opt$out_dir, winslash = "/", mustWork = FALSE),
                              normalizePath(DIR_VF, winslash = "/", mustWork = FALSE)) &&
                    !nzchar(suffix_part)) DIR_PLOTS_VF else file.path(opt$out_dir, "plots")
is_primary_run <- identical(normalizePath(vf_pa_file, winslash = "/", mustWork = FALSE),
                            normalizePath(FILE_VF_PA, winslash = "/", mustWork = FALSE))
stage_suffix <- if (nzchar(suffix_clean)) paste0("_", suffix_clean) else ""

log_info("Starting 02_gene_presence_analysis.R")
log_info("Thresholds: ID >= ", opt$min_id, "%, Coverage >= ", opt$min_cov, "%")
log_info("Assembly selection: ", opt$selection_file, " [", opt$selection_column, "]")
log_info("VF outputs: hits=", vf_hits_file, "; matrix=", vf_pa_file)

# 3. Load Metadata
# ------------------------------------------------------------------------------
assembly_df <- load_metadata()
assembly_df <- apply_manual_sample_curation(assembly_df, context = "vf_abricate_metadata")

curation_excluded <- assembly_df %>%
  filter(!(analysis_include_primary %in% TRUE) | !(genomics_expected_include %in% TRUE))
if (nrow(curation_excluded) > 0) {
  write_csv(curation_excluded, file.path(DIR_QC, "vf_abricate_manual_curation_excluded_rows.csv"))
  log_info("Manual curation excludes ", nrow(curation_excluded), " metadata row(s) from active VF profiling denominators.")
}
assembly_df <- filter_primary_genomics(assembly_df)

# Filter out assemblies with missing Participant_id (data quality check)
n_before <- nrow(assembly_df)
assembly_df <- assembly_df %>%
  filter(!is.na(Participant_id))
n_after <- nrow(assembly_df)
if (n_before > n_after) {
  log_info("Excluded %d assemblies with missing Participant_id", n_before - n_after)
}

# Ensure full_path exists
if (!"full_path" %in% names(assembly_df)) {
  if (!"file_name" %in% names(assembly_df)) {
    stop("Metadata lacks 'full_path' and 'file_name'.")
  }
  assembly_df <- assembly_df %>%
    mutate(
      full_path = file.path(DIR_FASTAS, file_name)
    )
}

assembly_df <- assembly_df %>%
  mutate(found = !is.na(full_path) & file.exists(full_path))

if (any(!assembly_df$found)) {
  n_missing <- sum(!assembly_df$found)
  warning("Missing ", n_missing, " FASTA files. They will be skipped.")
  assembly_df <- assembly_df %>% filter(found)
}
assembly_df <- select(assembly_df, -found)

assembly_df <- assembly_df %>%
  mutate(
    tp_lab = if ("tp_lab" %in% names(.)) normalise_timepoint_preserve_events(tp_lab) else normalise_timepoint_preserve_events(Timepoint),
    Event_type = if ("Event_type" %in% names(.)) as.character(Event_type) else episode_event_type(tp_lab),
    Collection_Date = if ("Collection_Date" %in% names(.)) as.character(Collection_Date) else NA_character_,
    Episode_ID = if ("Episode_ID" %in% names(.)) as.character(Episode_ID) else build_episode_id(., timepoint_col = "tp_lab", event_col = "Event_type", date_col = "Collection_Date"),
    tp_num = suppressWarnings(readr::parse_number(tp_lab))
  )

selection_file <- opt$selection_file
if (!file.exists(selection_file)) {
  stop(selection_file, " not found. VF profiling cannot fall back to raw assembly metadata; run 12a_wgs_qc.R first.")
}
selection <- read_csv(selection_file, show_col_types = FALSE)
if (!opt$selection_column %in% names(selection)) {
  stop("Selection file lacks requested column: ", opt$selection_column)
}
if (!"full_path" %in% names(selection)) {
  stop("Selection file lacks required full_path column: ", selection_file)
}
if (!"QC_PASS" %in% names(selection)) stop("Selection file lacks required QC_PASS column: ", selection_file)
selection <- selection %>%
  filter(.data[[opt$selection_column]] %in% TRUE, QC_PASS %in% TRUE) %>%
  mutate(full_path = normalizePath(full_path, winslash = "/", mustWork = FALSE))
selection$.analysis_assembler <- normalise_assembler_column(selection)
if (!nrow(selection)) stop("Selection file contains no selected QC-passing Longcycler rows: ", selection_file)
if (any(selection$.analysis_assembler != ANALYSIS_ASSEMBLER | is.na(selection$.analysis_assembler))) {
  stop("Selection file contains a non-Longcycler row. VF profiling does not permit assembler fallback.")
}
if (any(!file.exists(selection$full_path))) stop("Selection file contains missing FASTA paths.")
selection_dupes <- selection %>% count(Participant_id, tp_lab, name = "n") %>% filter(n != 1L)
if (nrow(selection_dupes)) stop("Selection file must contain exactly one Longcycler FASTA per participant-timepoint.")

n_before_canonical <- nrow(assembly_df)
assembly_df <- assembly_df %>%
  mutate(full_path = normalizePath(full_path, winslash = "/", mustWork = FALSE))
assembly_df_selected <- assembly_df %>%
  semi_join(selection %>% select(full_path), by = "full_path")
selection_extra <- selection %>%
  select(-.analysis_assembler) %>%
  anti_join(assembly_df %>% select(full_path), by = "full_path")
if (nrow(selection_extra) > 0) {
  missing_cols <- setdiff(names(assembly_df), names(selection_extra))
  for (col in missing_cols) selection_extra[[col]] <- NA
  selection_extra <- selection_extra %>% select(all_of(names(assembly_df)))
  assembly_df <- bind_rows(assembly_df_selected, selection_extra)
  log_info("Added ", nrow(selection_extra), " selected Longcycler row(s) directly from the analysis manifest.")
} else {
  assembly_df <- assembly_df_selected
}
if (nrow(assembly_df) != nrow(selection)) {
  stop("VF input row count does not match the Longcycler selection manifest.")
}
log_info("Using Longcycler-only assemblies for VF profiling: ", nrow(assembly_df),
         " retained from ", n_before_canonical, " candidate assembly rows; no assembler fallback.")

append_denominator_summary(
  assembly_df,
  "02_gene_presence_analysis.R",
  paste0("vf_abricate_input", stage_suffix),
  "participant_timepoint",
  selection_file,
  "Selected QC-passing Longcycler assemblies with exact manifest paths"
)

# 5. Setup Output Directories
# ------------------------------------------------------------------------------
DIR_ABRICATE <- if (!is.na(opt$abricate_dir) && nzchar(opt$abricate_dir)) opt$abricate_dir else file.path(opt$out_dir, "abricate")
ensure_dir(DIR_ABRICATE)
ensure_dir(plot_out_dir)

# 6. Run Abricate (VFDB)
# ------------------------------------------------------------------------------
if (!requireNamespace("digest", quietly = TRUE)) stop("Package 'digest' is required for content-bound VF caching.")
abricate_path <- Sys.which("abricate")
if (!nzchar(abricate_path)) stop("Abricate not found in PATH.")
abricate_version <- tryCatch(
  paste(system2(abricate_path, "--version", stdout = TRUE, stderr = TRUE), collapse = " "),
  error = function(e) "unknown"
)

read_abricate_cache <- function(path) {
  if (!file.exists(path) || file.size(path) == 0) return(tibble::tibble())
  readr::read_tsv(path, show_col_types = FALSE, progress = FALSE)
}

run_abr_cached <- function(fasta, db = "vfdb", min_cov = opt$min_cov, min_id = opt$min_id) {
  fasta <- normalizePath(fasta, winslash = "/", mustWork = TRUE)
  fasta_sha256 <- unname(digest::digest(fasta, algo = "sha256", file = TRUE))
  cache_schema <- "abricate_vfdb_sha256_v1"
  signature <- digest::digest(
    paste(cache_schema, fasta, fasta_sha256, db, min_cov, min_id, abricate_version, sep = "\n"),
    algo = "sha256", serialize = FALSE
  )
  cache_name <- paste0(basename(fasta), ".", substr(signature, 1, 20), ".tsv")
  cache <- file.path(DIR_ABRICATE, cache_name)
  sidecar <- paste0(cache, ".provenance.csv")
  expected <- tibble::tibble(
    cache_schema = cache_schema,
    cache_signature = signature,
    fasta_path = fasta,
    fasta_sha256 = fasta_sha256,
    fasta_size = as.character(file.size(fasta)),
    database = db,
    min_identity = as.character(min_id),
    min_coverage = as.character(min_cov),
    abricate_path = normalizePath(abricate_path, winslash = "/", mustWork = TRUE),
    abricate_version = abricate_version
  )

  if (file.exists(cache) && file.exists(sidecar)) {
    observed <- tryCatch(readr::read_csv(sidecar, show_col_types = FALSE, col_types = cols(.default = "c")), error = function(e) NULL)
    provenance_matches <- !is.null(observed) && nrow(observed) == 1L &&
      all(names(expected) %in% names(observed)) &&
      "result_sha256" %in% names(observed) &&
      all(vapply(names(expected), function(nm) {
        identical(as.character(observed[[nm]][1]), as.character(expected[[nm]][1]))
      }, logical(1))) &&
      identical(as.character(observed$result_sha256[1]), unname(digest::digest(cache, algo = "sha256", file = TRUE)))
    if (provenance_matches) return(read_abricate_cache(cache))
  }

  res <- processx::run(
    abricate_path,
    c("--quiet", "--db", db, "--mincov", as.character(min_cov), "--minid", as.character(min_id), fasta),
    echo = FALSE,
    error_on_status = FALSE
  )
  if (res$status != 0L) {
    stop("Abricate failed for ", fasta, " (status ", res$status, "): ", trimws(res$stderr))
  }

  writeLines(res$stdout, cache, useBytes = TRUE)
  hits <- read_abricate_cache(cache)
  write_csv(expected %>% mutate(
    result_sha256 = unname(digest::digest(cache, algo = "sha256", file = TRUE)),
    n_hits = as.character(nrow(hits)),
    result_status = ifelse(nrow(hits), "hits", "zero_hits")
  ), sidecar)
  hits
}

# Parallel Execution
future::plan(future::multisession, workers = CORES_USE)
on.exit(future::plan(sequential), add = TRUE)

message("Running Abricate on ", nrow(assembly_df), " assemblies...")
vf_results <- assembly_df %>%
  mutate(vfdb = furrr::future_map(full_path, run_abr_cached, .progress = TRUE))
vf_nonempty <- vf_results %>% filter(purrr::map_int(vfdb, nrow) > 0L)
vf_hits_all <- if (nrow(vf_nonempty) > 0) {
  tidyr::unnest(vf_nonempty, vfdb)
} else {
  assembly_df[0, , drop = FALSE] %>% mutate(GENE = character())
}

provenance_files <- list.files(DIR_ABRICATE, pattern = "\\.provenance\\.csv$", full.names = TRUE)
cache_provenance <- if (length(provenance_files)) {
  provenance_files %>%
    purrr::map_dfr(~ read_csv(.x, show_col_types = FALSE, col_types = cols(.default = "c"))) %>%
    filter(fasta_path %in% assembly_df$full_path)
} else {
  tibble::tibble()
}
write_csv(cache_provenance, file.path(opt$out_dir, paste0("vf_abricate_cache_provenance", suffix_part, ".csv")))

# Clean up column names
if (!"tp_lab" %in% names(vf_hits_all)) {
  # Re-join if lost (unlikely with mutate above, but safe)
  vf_hits_all <- vf_hits_all %>%
    left_join(select(assembly_df, full_path, tp_lab), by = "full_path")
}

# Standardize Gene Column
if (nrow(vf_hits_all) > 0) {
  gene_col <- intersect(c("GENE", "GENE_NAME", "NAME", "PRODUCT", "GENE SYMBOL"), names(vf_hits_all))[1]
  if (is.na(gene_col)) stop("No gene name column found in Abricate output.")
  if (gene_col != "GENE") vf_hits_all <- vf_hits_all %>% rename(GENE = all_of(gene_col))
}

saveRDS(vf_hits_all, vf_hits_file)
message("✓ Saved VF hits: ", vf_hits_file)

# 6. Generate Presence/Absence Matrix
# ------------------------------------------------------------------------------
# We aggregate Abricate hits to the participant-timepoint level after the
# Longcycler-only selection. Each episode is represented by one assembly.
vf_episode_base <- assembly_df %>% distinct(Participant_id, tp_lab, Episode_ID)
if (nrow(vf_hits_all) > 0) {
  vf_detected <- vf_hits_all %>%
    filter(!is.na(Participant_id), !is.na(tp_lab), !is.na(GENE)) %>%
    mutate(Episode_ID = as.character(Episode_ID)) %>%
    distinct(Participant_id, tp_lab, Episode_ID, GENE) %>%
    mutate(present = 1L) %>%
    pivot_wider(names_from = GENE, values_from = present, values_fill = 0L)
  vf_pa_all <- vf_episode_base %>%
    left_join(vf_detected, by = c("Participant_id", "tp_lab", "Episode_ID")) %>%
    mutate(across(-c(Participant_id, tp_lab, Episode_ID), ~ replace_na(as.integer(.x), 0L)))
} else {
  vf_pa_all <- vf_episode_base
}

readr::write_csv(vf_pa_all, vf_pa_file)

append_denominator_summary(
  vf_pa_all,
  "02_gene_presence_analysis.R",
  paste0("vf_pa_all", stage_suffix),
  "participant_timepoint",
  vf_pa_file,
  "Episode-level VF presence/absence matrix from canonical selected assemblies"
)
if (is_primary_run) {
  write_uti_attrition_outputs()
}

# 7. Gene Prevalence Stats
# ------------------------------------------------------------------------------
tbl_gene <- vf_hits_all %>%
  distinct(Participant_id, GENE) %>%
  dplyr::count(GENE, name = "n_participants") %>%
  arrange(desc(n_participants))

readr::write_csv(tbl_gene, stats_gene_file)

# 8. Quick Plots
# ------------------------------------------------------------------------------
top25 <- tbl_gene %>%
  slice_max(n_participants, n = 25) %>%
  mutate(GENE = fct_reorder(GENE, n_participants))

p1 <- ggplot(top25, aes(GENE, n_participants)) +
  geom_col(fill = "steelblue") +
  coord_flip() +
  labs(
    title = "Most frequently detected VFDB genes in selected E. coli assemblies",
    subtitle = "Selection criterion: top 25 genes by number of participants with at least one detection",
    y = "Participants with gene detected",
    x = NULL,
    caption = sprintf(
      "Data: %s. Denominator: participant-level presence from selected assembly VFDB calls; this plot is not stratified by UTI/Not_UTI status and shows prevalence, not virulence causality.",
      vf_pa_file
    )
  ) +
  theme_bw(base_size = 11) +
  theme(plot.caption = element_text(hjust = 0, size = 8, colour = "grey35"))
ggsave(file.path(plot_out_dir, paste0("core_bar_top25_all", suffix_part, ".png")), p1, width = 6, height = 6, dpi = 300)

p2 <- ggplot(tbl_gene, aes(n_participants)) +
  geom_histogram(binwidth = 1, fill = "grey70") +
  labs(
    title = "Distribution of VFDB gene prevalence across participants",
    subtitle = "Prevalence is counted as participants with at least one selected isolate carrying the gene",
    x = "Participants with gene detected",
    y = "Number of VFDB genes",
    caption = sprintf(
      "Data: %s. Denominator: participant-level gene prevalence from selected assembly VFDB calls. This is a descriptive input-QC plot and is not a UTI-vs-Not_UTI association test.",
      vf_pa_file
    )
  ) +
  theme_bw(base_size = 11) +
  theme(plot.caption = element_text(hjust = 0, size = 8, colour = "grey35"))
ggsave(file.path(plot_out_dir, paste0("core_histogram_all", suffix_part, ".png")), p2, width = 5, height = 4, dpi = 300)

# 9. Nucmer Trajectories (Pairwise)
# ------------------------------------------------------------------------------
ids_multi <- vf_pa_all %>%
  distinct(Participant_id, tp_lab) %>%
  mutate(tp_num = suppressWarnings(readr::parse_number(as.character(tp_lab)))) %>%
  filter(!is.na(tp_num)) %>%
  dplyr::count(Participant_id, name = "n_tp") %>%
  filter(n_tp >= 2) %>%
  pull(Participant_id)

if (opt$skip_nucmer) {
  message("Skipping Nucmer trajectories (--skip_nucmer).")
} else if (length(ids_multi) > 0) {
  DIR_NUCMER <- file.path(DIR_STRAIN, "nucmer")
  ensure_dir(DIR_NUCMER)

  if (Sys.which("nucmer") == "") {
    message("⚠ Nucmer not found. Skipping trajectories.")
  } else {
    message("Running Nucmer on ", length(ids_multi), " participants...")

    assembly_long <- assembly_df %>%
      filter(Participant_id %in% ids_multi) %>%
      filter(!is.na(tp_num))

    pair_tbl <- assembly_long %>%
      group_by(Participant_id) %>% # Assuming 1 assembler per participant usually, or group by assembler too
      arrange(tp_num, .by_group = TRUE) %>%
      mutate(
        path_A = full_path,
        path_B = lead(full_path),
        tp_A   = tp_lab,
        tp_B   = lead(tp_lab),
        # STRICT OUTPUT DIRECTORY: results/nucmer/[PID]...
        outdir = file.path(DIR_NUCMER, glue::glue("{Participant_id}_{tp_A}_vs_{tp_B}"))
      ) %>%
      filter(!is.na(path_B)) %>%
      ungroup()

    run_nucmer <- function(a, b, od) {
      ensure_dir(od)
      pref <- file.path(od, "run")

      # Capture output to avoid clutter
      system(glue::glue("nucmer --mum --prefix={shQuote(pref)} {shQuote(a)} {shQuote(b)} >/dev/null 2>&1"))
      system(glue::glue("delta-filter -1 {pref}.delta > {pref}.1delta 2>/dev/null"))
      system(glue::glue("dnadiff -p {pref}_dd -d {pref}.1delta >/dev/null 2>&1"))

      rpt <- glue::glue("{pref}_dd.report")
      if (!file.exists(rpt)) {
        return(tibble(AvgIdentity = NA_real_, TotalSnpCnt = NA_real_))
      }

      L <- readLines(rpt)
      grab <- \(k) {
        m <- grep(k, L, value = TRUE)
        if (length(m)) as.numeric(str_extract(m[1], "\\d+\\.?\\d*")) else NA_real_
      }
      tibble(AvgIdentity = grab("AvgIdentity"), TotalSnpCnt = grab("TotalSNPs"))
    }

    pair_df <- pair_tbl %>%
      mutate(res = furrr::future_pmap(list(path_A, path_B, outdir), run_nucmer, .progress = TRUE)) %>%
      unnest(res)

    # Plotting logic (simplified)
    ensure_dir(file.path(DIR_PLOTS_STRAIN, "pairwise_identity"))
    # ... (keeping plotting logic minimal for now, can be expanded)
  }
}

message("✓ Analysis complete.")

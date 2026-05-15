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
#   isolate’s VF repertoire, we can ask: do bacteria causing UTI carry
#   different VFs from those causing ASB?
#
# KEY DESIGN DECISIONS:
#   - Uses UNION logic: a gene is called “present” if detected in EITHER
#     the Flye OR Longcycler assembly for a given participant×timepoint.
#     This maximises sensitivity at the cost of potential false positives
#     from a single assembler.
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
#   → 22_vf_build_analysis_dataset.R joins vf_pa_all.csv with clinical status
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
    type = "character", default = file.path(DIR_QC, "canonical_assembly_selection.csv"),
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
if (file.exists(selection_file)) {
  selection <- read_csv(selection_file, show_col_types = FALSE)
  if (!opt$selection_column %in% names(selection)) {
    stop("Selection file lacks requested column: ", opt$selection_column)
  }
  if (!"full_path" %in% names(selection)) {
    stop("Selection file lacks required full_path column: ", selection_file)
  }
  selection <- selection %>%
    filter(.data[[opt$selection_column]] %in% TRUE) %>%
    mutate(full_path = normalizePath(full_path, winslash = "/", mustWork = FALSE))
  n_before_canonical <- nrow(assembly_df)
  assembly_df <- assembly_df %>%
    mutate(full_path = normalizePath(full_path, winslash = "/", mustWork = FALSE))
  assembly_df_selected <- assembly_df %>%
    semi_join(selection %>% select(full_path), by = "full_path")
  selection_extra <- selection %>%
    anti_join(assembly_df %>% select(full_path), by = "full_path")
  if (nrow(selection_extra) > 0) {
    missing_cols <- setdiff(names(assembly_df), names(selection_extra))
    for (col in missing_cols) selection_extra[[col]] <- NA
    selection_extra <- selection_extra %>%
      select(all_of(names(assembly_df)))
    assembly_df <- bind_rows(assembly_df_selected, selection_extra)
    log_info("Added ", nrow(selection_extra), " selected assembly row(s) directly from selection file.")
  } else {
    assembly_df <- assembly_df_selected
  }
  log_info("Using selected assemblies for VF profiling: ", nrow(assembly_df),
           " retained from ", n_before_canonical, " assembly-level rows.")
} else {
  log_info("WARNING: assembly selection file not found; VF profiling will use all usable assembly rows.")
}

append_denominator_summary(
  assembly_df,
  "02_gene_presence_analysis.R",
  paste0("vf_abricate_input", stage_suffix),
  "participant_timepoint",
  FILE_METADATA,
  paste0("Selected assemblies from ", selection_file)
)

# 5. Setup Output Directories
# ------------------------------------------------------------------------------
DIR_ABRICATE <- if (!is.na(opt$abricate_dir) && nzchar(opt$abricate_dir)) opt$abricate_dir else file.path(opt$out_dir, "abricate")
ensure_dir(DIR_ABRICATE)
ensure_dir(plot_out_dir)

# 6. Run Abricate (VFDB)
# ------------------------------------------------------------------------------
run_abr_cached <- function(fasta, db = "vfdb", min_cov = opt$min_cov, min_id = opt$min_id) {
  # Cache filename includes thresholds to avoid stale results if params change
  cache_name <- paste0(basename(fasta), ".vfdb.id", min_id, ".cov", min_cov, ".tsv")
  cache <- file.path(DIR_ABRICATE, cache_name)

  if (file.exists(cache)) {
    return(readr::read_tsv(cache, show_col_types = FALSE, progress = FALSE))
  }

  # Check if abricate is available
  if (Sys.which("abricate") == "") {
    stop("Abricate not found in PATH.")
  }

  cmd <- glue::glue("abricate --quiet --db {db} --mincov {min_cov} --minid {min_id} {shQuote(fasta)} > {shQuote(cache)}")
  exit <- system(cmd)
  if (exit != 0) warning("Abricate non-zero exit: ", basename(fasta))

  readr::read_tsv(cache, show_col_types = FALSE, progress = FALSE)
}

# Parallel Execution
future::plan(future::multisession, workers = CORES_USE)
on.exit(future::plan(sequential), add = TRUE)

safe_abr <- purrr::safely(run_abr_cached, otherwise = NULL, quiet = TRUE)

message("Running Abricate on ", nrow(assembly_df), " assemblies...")
vf_hits_all <- assembly_df %>%
  mutate(vfdb = furrr::future_map(full_path, ~ safe_abr(.x)$result, .progress = TRUE)) %>%
  filter(purrr::map_lgl(vfdb, ~ !is.null(.x) && NROW(.x) > 0)) %>%
  tidyr::unnest(vfdb)

# Clean up column names
if (!"tp_lab" %in% names(vf_hits_all)) {
  # Re-join if lost (unlikely with mutate above, but safe)
  vf_hits_all <- vf_hits_all %>%
    left_join(select(assembly_df, full_path, tp_lab), by = "full_path")
}

# Standardize Gene Column
gene_col <- intersect(c("GENE", "GENE_NAME", "NAME", "PRODUCT", "GENE SYMBOL"), names(vf_hits_all))[1]
if (is.na(gene_col)) stop("No gene name column found in Abricate output.")
vf_hits_all <- vf_hits_all %>% rename(GENE = all_of(gene_col))

saveRDS(vf_hits_all, vf_hits_file)
message("✓ Saved VF hits: ", vf_hits_file)

# 6. Generate Presence/Absence Matrix
# ------------------------------------------------------------------------------
# We aggregate Abricate hits to the participant-timepoint level after canonical
# assembly selection.  This avoids treating flye and longcycler alternatives as
# independent biological episodes.
vf_pa_all <- vf_hits_all %>%
  filter(!is.na(Participant_id), !is.na(tp_lab), !is.na(GENE)) %>%
  mutate(
    Episode_ID = if ("Episode_ID" %in% names(.)) as.character(Episode_ID) else build_episode_id(., timepoint_col = "tp_lab", event_col = "Event_type", date_col = "Collection_Date")
  ) %>%
  distinct(Participant_id, tp_lab, Episode_ID, GENE) %>%
  mutate(present = 1) %>%
  pivot_wider(names_from = GENE, values_from = present, values_fill = 0)

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
      "Data: %s. Denominator: participant-level presence from selected assembly VFDB calls; this plot is not stratified by ASB/UTI status and shows prevalence, not virulence causality.",
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
      "Data: %s. Denominator: participant-level gene prevalence from selected assembly VFDB calls. This is a descriptive input-QC plot and is not an ASB-vs-UTI association test.",
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

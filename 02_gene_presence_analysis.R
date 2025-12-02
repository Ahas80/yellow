#!/usr/bin/env Rscript
# ==============================================================================
# 02_gene_presence_analysis.R
# ------------------------------------------------------------------------------
# Role: [Descriptive/Exploratory] - Profile virulence factors (VFDB) using Abricate.
#
# Inputs:
#   - assembly_metadata.csv
#   - data/assemblies/*.fasta
#
# Outputs:
#   - results/vf/vf_hits_all.rds
#   - results/vf/vf_pa_all.csv
#   - results/vf/stats_gene_level.csv
#   - results/vf/abricate/ (cache)
#   - results/strain_compare/nucmer/ (pairwise alignment)
#
# Usage:
#   Rscript 02_gene_presence_analysis.R [--min_id 80] [--min_cov 80]
#
# Biological/Statistical purpose:
#   - Detects virulence factors to determine the pathogenic potential of isolates.
#   - Generates presence/absence matrices for downstream association testing.
#   - Performs initial pairwise genomic comparisons (Nucmer) for strain tracking.
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
  )
)

opt <- parse_args(OptionParser(option_list = option_list))

log_info <- function(...) message(format(Sys.time(), "[%H:%M:%S] "), ...)

log_info("Starting 02_gene_presence_analysis.R")
log_info("Thresholds: ID >= ", opt$min_id, "%, Coverage >= ", opt$min_cov, "%")

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
      full_path = file.path(DIR_FASTAS, file_name),
      found = file.exists(full_path)
    )

  if (any(!assembly_df$found)) {
    bad <- assembly_df %>%
      filter(!found) %>%
      pull(full_path)
    stop("Missing FASTA files:\n", paste(head(bad), collapse = "\n"))
  }
  assembly_df <- select(assembly_df, -found)
}

# 4. Helper: Normalize Timepoints
# ------------------------------------------------------------------------------
tp_norm <- function(x) {
  x <- as.character(x)
  is_uricult <- stringr::str_detect(x, stringr::regex("uricult", ignore_case = TRUE))
  tp_num <- as.integer(stringr::str_extract(x, "\\d+"))
  tp_num[is_uricult] <- NA_integer_

  tp_lab <- dplyr::case_when(
    is_uricult ~ "Uricult",
    !is.na(tp_num) ~ paste0("T", tp_num),
    TRUE ~ "Unscheduled"
  )

  # Define levels for correct ordering
  tp_levels <- c(
    paste0("T", sort(unique(tp_num[!is.na(tp_num)]))),
    "Uricult", "Unscheduled"
  )

  tibble::tibble(
    tp_lab = factor(tp_lab, levels = tp_levels),
    tp_num = tp_num
  )
}

# Apply normalization
tp <- tp_norm(assembly_df$Timepoint)
assembly_df <- assembly_df %>%
  mutate(tp_lab = tp$tp_lab, tp_num = tp$tp_num)

# 5. Setup Output Directories
# ------------------------------------------------------------------------------
DIR_ABRICATE <- file.path(DIR_VF, "abricate")
ensure_dir(DIR_ABRICATE)
ensure_dir(DIR_PLOTS_VF)

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

saveRDS(vf_hits_all, FILE_VF_HITS)
message("✓ Saved VF hits: ", FILE_VF_HITS)

# 6. Generate Presence/Absence Matrix
# ------------------------------------------------------------------------------
vf_pa_all <- vf_hits_all %>%
  filter(!is.na(Participant_id), !is.na(tp_lab), !is.na(GENE)) %>%
  distinct(Participant_id, tp_lab, GENE) %>%
  mutate(present = 1) %>%
  pivot_wider(names_from = GENE, values_from = present, values_fill = 0)

readr::write_csv(vf_pa_all, FILE_VF_PA)

# 7. Gene Prevalence Stats
# ------------------------------------------------------------------------------
tbl_gene <- vf_hits_all %>%
  distinct(Participant_id, GENE) %>%
  count(GENE, name = "n_participants") %>%
  arrange(desc(n_participants))

readr::write_csv(tbl_gene, file.path(DIR_VF, "stats_gene_level.csv"))

# 8. Quick Plots
# ------------------------------------------------------------------------------
top25 <- tbl_gene %>%
  slice_max(n_participants, n = 25) %>%
  mutate(GENE = fct_reorder(GENE, n_participants))

p1 <- ggplot(top25, aes(GENE, n_participants)) +
  geom_col(fill = "steelblue") +
  coord_flip() +
  labs(title = "Top 25 VFDB genes", y = "Number of Participants", x = NULL) +
  theme_minimal(base_size = 11)
ggsave(file.path(DIR_PLOTS_VF, "core_bar_top25_all.png"), p1, width = 6, height = 6, dpi = 300)

p2 <- ggplot(tbl_gene, aes(n_participants)) +
  geom_histogram(binwidth = 1, fill = "grey70") +
  labs(title = "VF gene prevalence", x = "Number of Participants", y = "Gene count") +
  theme_minimal(base_size = 11)
ggsave(file.path(DIR_PLOTS_VF, "core_histogram_all.png"), p2, width = 5, height = 4, dpi = 300)

# 9. Nucmer Trajectories (Pairwise)
# ------------------------------------------------------------------------------
ids_multi <- vf_pa_all %>%
  distinct(Participant_id, tp_lab) %>%
  mutate(tp_num = suppressWarnings(readr::parse_number(as.character(tp_lab)))) %>%
  filter(!is.na(tp_num)) %>%
  count(Participant_id, name = "n_tp") %>%
  filter(n_tp >= 2) %>%
  pull(Participant_id)

if (length(ids_multi) > 0) {
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

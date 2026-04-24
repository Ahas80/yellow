#!/usr/bin/env Rscript
# ==============================================================================
# 08_core_vs_plasmid.R
# ==============================================================================
#
# GOAL:
#   Compare chromosomal ST lineages with plasmid replicon types to assess
#   whether specific plasmids are associated with particular E. coli lineages.
#   This informs whether plasmid-carried VFs are lineage-linked or independently
#   mobile across different STs.
#
# ------------------------------------------------------------------------------
# Role: [Inferential-core] - Compare chromosomal STs with plasmid types.
#
# Inputs:
#   - results/mlst/mlst_all.tsv
#   - data/assemblies/*.fasta
#
# Outputs:
#   - results/mlst/ST_core_freq.csv
#   - results/mlst/pMLST_hits_long.csv OR plasmid_replicons_long.csv
#   - results/mlst/plasmid_types_per_isolate.csv
#   - results/mlst/ST_plasmid_associations.csv
#
# Usage:
#   Rscript 08_core_vs_plasmid.R
#
# Biological/Statistical purpose:
#   - Investigates associations between bacterial lineages (STs) and plasmid types.
#   - Tests for significant co-occurrence (e.g., IncF plasmids in ST131).
# ==============================================================================

# 1. Load Configuration & Libraries
source("00_config.R")
suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
  library(purrr)
  library(furrr)
  library(fs)
  library(stringr)
  library(scales)
  library(processx)
})

# 2. Configuration
# ------------------------------------------------------------------------------
# 2. Configuration
# ------------------------------------------------------------------------------
FILE_MLST_ALL <- FILE_MLST_ALL # From 00_config.R
DIR_PMLST_LOG <- file.path(DIR_MLST, "pmlst_logs")
ensure_dir(DIR_PMLST_LOG)

check_tool("mlst")

# 3. Chromosomal ST Frequencies
# ------------------------------------------------------------------------------
if (file.exists(FILE_MLST_ALL)) {
  core_tbl <- read_tsv(FILE_MLST_ALL, show_col_types = FALSE)

  if ("ST" %in% names(core_tbl)) {
    core_tbl %>%
      count(ST, sort = TRUE) %>%
      mutate(pct = percent(n / sum(n))) %>%
      write_csv(file.path(DIR_MLST, "ST_core_freq.csv"))
  }
}

# 4. pMLST Setup
# ------------------------------------------------------------------------------
# Check for plasmid schemes in mlst
schemes_out <- system2("mlst", "--list", stdout = TRUE, stderr = TRUE)
writeLines(schemes_out, file.path(DIR_PMLST_LOG, "mlst_list.txt"))

pSchemes <- schemes_out[grepl("inc|plasmid|pmlst", schemes_out, ignore.case = TRUE)]
pSchemes <- sub("\\s.*$", "", trimws(pSchemes))
pSchemes <- unique(pSchemes[nzchar(pSchemes)])

fasta_files <- dir_ls(DIR_FASTAS, glob = "*.fasta")
if (length(fasta_files) == 0) fasta_files <- dir_ls(DIR_FASTAS, glob = "*.fa")

if (length(fasta_files) == 0) stop("No FASTA files found in ", DIR_FASTAS)

# 5. Run pMLST or Fallback
# ------------------------------------------------------------------------------
if (length(pSchemes) > 0) {
  msg("Running pMLST with schemes: ", paste(pSchemes, collapse = ", "))

  run_pmlst <- function(fasta, scheme) {
    basename <- fs::path_file(fasta)
    out_csv <- fs::path(DIR_MLST, "raw", paste0(basename, ".", scheme, ".csv"))

    if (file.exists(out_csv)) {
      return(read_csv(out_csv, show_col_types = FALSE))
    }

    res <- processx::run("mlst", c("--quiet", "--threads", "1", "--scheme", scheme, "--csv", "--legacy", fasta), error_on_status = FALSE)

    if (res$status == 0 && grepl(",", res$stdout)) {
      dat <- read_csv(I(res$stdout), show_col_types = FALSE) %>%
        rename_with(tolower) %>%
        mutate(scheme = scheme, isolate_id = tools::file_path_sans_ext(basename), .before = 1)
      write_csv(dat, out_csv)
      return(dat)
    }
    return(NULL)
  }

  future::plan(future::multisession, workers = CORES_USE)
  grid <- expand_grid(file = fasta_files, scheme = pSchemes)
  pmlst_hits <- future_pmap_dfr(list(grid$file, grid$scheme), ~ run_pmlst(..1, ..2), .progress = TRUE)
  future::plan(future::sequential)

  if (nrow(pmlst_hits) > 0) {
    write_csv(pmlst_hits, file.path(DIR_MLST, "pMLST_hits_long.csv"))

    pmlst_wide <- pmlst_hits %>%
      select(Isolate_ID = isolate_id, p_scheme = scheme, st) %>%
      pivot_wider(names_from = p_scheme, values_from = st, values_fill = NA, names_prefix = "pST_")
    write_csv(pmlst_wide, file.path(DIR_MLST, "plasmid_types_per_isolate.csv"))
  }
} else {
  msg("No pMLST schemes found. Falling back to ABRicate (PlasmidFinder).")
  check_tool("abricate")

  run_pf <- function(fasta) {
    res <- processx::run("abricate", c("--quiet", "--db", "plasmidfinder", fasta), error_on_status = FALSE, stderr_to_stdout = TRUE)
    if (res$status != 0 || !nzchar(res$stdout)) {
      return(NULL)
    }

    tab <- read_tsv(I(res$stdout), show_col_types = FALSE)
    if (nrow(tab) == 0) {
      return(NULL)
    }

    tab %>%
      mutate(isolate_id = tools::file_path_sans_ext(fs::path_file(fasta))) %>%
      select(isolate_id, replicon = GENE, identity = `%IDENTITY`, coverage = `%COVERAGE`)
  }

  future::plan(future::multisession, workers = CORES_USE)
  pf_hits <- future_map_dfr(fasta_files, run_pf, .progress = TRUE)
  future::plan(future::sequential)

  if (nrow(pf_hits) > 0) {
    write_csv(pf_hits, file.path(DIR_MLST, "plasmid_replicons_long.csv"))

    pf_wide <- pf_hits %>%
      distinct(isolate_id, replicon) %>%
      mutate(present = TRUE) %>%
      pivot_wider(names_from = replicon, values_from = present, values_fill = FALSE)
    write_csv(pf_wide, file.path(DIR_MLST, "plasmid_replicons_wide.csv"))
  }
}


# 6. Statistical Association (ST vs Plasmid)
# ------------------------------------------------------------------------------
# [STAT] Test for significant associations between major STs and plasmid types
# (e.g., Is IncF significantly enriched in ST131?)

# Determine which plasmid data is available
plasmid_data <- NULL
if (exists("pmlst_wide")) {
  plasmid_data <- pmlst_wide
} else if (exists("pf_wide")) {
  plasmid_data <- pf_wide
}

if (!is.null(plasmid_data) && exists("core_tbl")) {
  msg("Running statistical tests for ST-plasmid associations...")

  # Join ST and Plasmid data
  # Ensure Isolate_ID matching (case-insensitive or exact)
  # core_tbl has Isolate_ID, plasmid_data has Isolate_ID or isolate_id

  # Normalize ID column names
  if ("isolate_id" %in% names(plasmid_data)) plasmid_data <- rename(plasmid_data, Isolate_ID = isolate_id)

  merged <- inner_join(core_tbl, plasmid_data, by = "Isolate_ID")

  if (nrow(merged) > 10) {
    # Identify top STs and top Plasmids to test (avoid testing rare things)
    top_STs <- names(sort(table(merged$ST), decreasing = TRUE))[1:min(5, n_distinct(merged$ST))]

    # Identify plasmid columns (exclude ID and ST)
    plasmid_cols <- setdiff(names(plasmid_data), c("Isolate_ID", "p_scheme", "st"))

    # Function to test one ST vs one Plasmid
    test_assoc <- function(target_st, target_plasmid) {
      # Create 2x2 table: ST vs Not-ST, Plasmid vs Not-Plasmid
      st_binary <- merged$ST == target_st
      plasmid_binary <- !is.na(merged[[target_plasmid]]) & merged[[target_plasmid]] != FALSE & merged[[target_plasmid]] != 0

      tbl <- table(st_binary, plasmid_binary)

      if (all(dim(tbl) == 2)) {
        res <- fisher.test(tbl)
        tibble(
          ST = target_st,
          Plasmid = target_plasmid,
          OR = res$estimate,
          p_value = res$p.value,
          n_cooccur = tbl[2, 2]
        )
      } else {
        NULL
      }
    }

    # Run tests
    results <- list()
    for (st in top_STs) {
      for (pl in plasmid_cols) {
        # Only test if plasmid is present in >5% of isolates
        if (mean(!is.na(merged[[pl]]) & merged[[pl]] != FALSE) > 0.05) {
          results[[length(results) + 1]] <- test_assoc(st, pl)
        }
      }
    }

    assoc_results <- bind_rows(results)

    if (nrow(assoc_results) > 0) {
      assoc_results <- assoc_results %>%
        mutate(FDR = p.adjust(p_value, method = "BH")) %>%
        arrange(p_value)

      write_csv(assoc_results, file.path(DIR_MLST, "ST_plasmid_associations.csv"))
      msg("Saved association results to ST_plasmid_associations.csv")
    }
  }
}

msg("✓ Core vs Plasmid analysis complete.")

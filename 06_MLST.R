#!/usr/bin/env Rscript
# ==============================================================================
# 06_MLST.R
# ------------------------------------------------------------------------------
# Role: [Typing] - Perform Multi-Locus Sequence Typing (MLST) on all assemblies.
#
# Inputs:
#   - assembly_metadata.csv
#   - data/assemblies/*.fasta
#
# Outputs:
#   - results/mlst/mlst_all.tsv
#   - results/mlst/mlst_matrix.csv
#   - results/mlst/mlst_qc_summary.csv
#   - results/mlst/top_STs.csv
#   - results/mlst/raw/ (cache)
#   - results/mlst/logs/
#
# Usage:
#   Rscript 06_MLST.R
#
# Biological/Statistical purpose:
#   - Assigns Sequence Types (STs) to define bacterial lineages.
#   - Evaluates ST persistence within participants over time.
# ==============================================================================

# 1. Load Configuration & Libraries
source("00_config.R")
suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
  library(purrr)
  library(furrr)
  library(processx)
  library(stringr)
  library(fs)
  library(scales)
  library(rlang)
  library(ggplot2)
})

# 2. Configuration
# ------------------------------------------------------------------------------
SCHEME <- "ecoli" # PubMLST scheme, or "auto"
DEBUG <- FALSE
THREADS <- CORES_USE # From 00_config.R

# Logging helper
msg <- function(...) message(format(Sys.time(), "%H:%M:%S "), ...)

# Check tools
check_tool("mlst")
check_tool("blastn")

# Output directories
DIR_MLST <- file.path(DIR_RESULTS, "mlst")
DIR_MLST_RAW <- file.path(DIR_MLST, "raw")
DIR_MLST_LOG <- file.path(DIR_MLST, "logs")
ensure_dir(DIR_MLST)
ensure_dir(DIR_MLST_RAW)
ensure_dir(DIR_MLST_LOG)
ensure_dir(DIR_PLOTS_MLST)

# 3. Load Assembly Metadata
# ------------------------------------------------------------------------------
if (!file.exists(FILE_METADATA)) stop("Missing ", FILE_METADATA)
assembly_df <- read_csv(FILE_METADATA, show_col_types = FALSE)

# Ensure full_path exists and is valid
if (!"full_path" %in% names(assembly_df)) {
  # Try to reconstruct if missing (should be handled by 00_make_assembly_metadata.r, but robustify here)
  if ("file_name" %in% names(assembly_df)) {
    assembly_df$full_path <- file.path(DIR_FASTAS, assembly_df$file_name)
  } else {
    stop("assembly_metadata.csv missing 'full_path' and 'file_name'")
  }
}

# Verify files exist
assembly_df$found <- file.exists(assembly_df$full_path)
if (any(!assembly_df$found)) {
  missing_files <- assembly_df$full_path[!assembly_df$found]
  warning("Missing ", length(missing_files), " FASTA files. They will be skipped.")
  assembly_df <- assembly_df %>% filter(found)
}

msg("Processing ", nrow(assembly_df), " assemblies.")

# 4. Run MLST
# ------------------------------------------------------------------------------
run_mlst <- function(fasta, scheme = SCHEME) {
  basename <- fs::path_file(fasta)
  cache <- fs::path(DIR_MLST_RAW, paste0(basename, ".mlst.csv"))
  log_out <- fs::path(DIR_MLST_LOG, paste0(basename, ".log.txt"))

  if (!fs::file_exists(cache)) {
    cmd <- c("--quiet", "--threads", "1", "--scheme", scheme, "--csv", "--legacy", fasta)

    tryCatch(
      {
        px <- processx::run("mlst", cmd, echo = DEBUG, stderr_to_stdout = TRUE, error_on_status = FALSE)
        write_lines(px$stdout, log_out)

        if (px$status != 0) {
          return(tibble(.status = "failed"))
        }
        write_lines(px$stdout, cache)
      },
      error = function(e) {
        write_lines(as.character(e), log_out, append = TRUE)
        return(tibble(.status = "error"))
      }
    )
  }

  # Read cached result
  if (!fs::file_exists(cache)) {
    return(tibble(.status = "failed"))
  }

  dat <- read_csv(cache, col_types = cols(.default = col_character()), na = c("", "?"), progress = FALSE, show_col_types = FALSE) %>%
    rename_with(tolower) %>%
    mutate(scheme = scheme, .before = 1)

  if (nrow(dat) == 0) {
    return(tibble(.status = "no-data"))
  }

  st_cols <- names(dat)[tolower(names(dat)) == "st"]
  if (!length(st_cols)) {
    return(tibble(.status = "unexpected-format"))
  }

  mutate(dat, .status = "ok")
}

# Parallel Execution
future::plan(future::multisession, workers = THREADS)
mlst_raw <- assembly_df %>%
  mutate(mlst = future_map(full_path, run_mlst, .progress = TRUE)) %>%
  unnest(mlst)
future::plan(future::sequential)

# 5. Summary & QC
# ------------------------------------------------------------------------------
summary_tbl <- mlst_raw %>%
  count(.status, name = "n") %>%
  mutate(pct = percent(n / sum(n)))
msg("Run summary:")
print(summary_tbl)

if (!any(mlst_raw$.status == "ok")) stop("No successful typings!")

# Save provenance
write_lines(capture.output(sessionInfo()), file.path(DIR_MLST, "sessionInfo.txt"))
try(write_lines(system2("mlst", "--version", stdout = TRUE), file.path(DIR_MLST, "mlst_version.txt")), silent = TRUE)

# Tidy Output
st_col <- names(mlst_raw)[tolower(names(mlst_raw)) == "st"]
mlst_tbl <- mlst_raw %>%
  filter(.status == "ok") %>%
  select(-.status) %>%
  mutate(ST = .data[[st_col[1]]]) %>%
  select(-all_of(st_col))

# Ensure Isolate_ID
if (!"Isolate_ID" %in% names(mlst_tbl)) {
  mlst_tbl$Isolate_ID <- tools::file_path_sans_ext(basename(mlst_tbl$full_path))
}
mlst_tbl <- mlst_tbl %>% relocate(Isolate_ID, ST, .before = 1)

# QC Flags
meta_cols <- unique(c("scheme", "ST", "Isolate_ID", "file_name", "full_path", names(assembly_df)))
locus_cols <- setdiff(names(mlst_tbl), meta_cols)

if (length(locus_cols) > 0) {
  write_lines(locus_cols, file.path(DIR_MLST, "mlst_locus_list.txt"))

  mlst_tbl <- mlst_tbl %>%
    mutate(across(all_of(locus_cols), as.character)) %>%
    rowwise() %>%
    mutate(
      n_loci_typed    = sum(!is.na(c_across(all_of(locus_cols))) & c_across(all_of(locus_cols)) != "" & !grepl("^(\\?|0)$", c_across(all_of(locus_cols)))),
      mlst_complete   = n_loci_typed == length(locus_cols),
      has_new_allele  = any(grepl("NEW|\\*$", c_across(all_of(locus_cols)), ignore.case = TRUE)),
      ambiguous_call  = any(grepl("[,;/]", c_across(all_of(locus_cols))))
    ) %>%
    ungroup()

  mlst_qc <- mlst_tbl %>% count(mlst_complete, has_new_allele, ambiguous_call, name = "n")
  write_csv(mlst_qc, file.path(DIR_MLST, "mlst_qc_summary.csv"))
}

# 6. Top STs
# ------------------------------------------------------------------------------
top_STs <- mlst_tbl %>%
  mutate(ST = as.character(ST)) %>%
  count(ST, sort = TRUE, name = "n_isolates") %>%
  mutate(pct = n_isolates / sum(n_isolates))
write_csv(top_STs, file.path(DIR_MLST, "top_STs.csv"))

# 7. Persistence (if applicable)
# ------------------------------------------------------------------------------
if ("Participant_id" %in% names(mlst_tbl) && "Timepoint" %in% names(mlst_tbl)) {
  tp_norm <- function(x) {
    x <- as.character(x)
    is_uricult <- str_detect(x, regex("uricult", ignore_case = TRUE))
    tp_num <- suppressWarnings(as.integer(str_extract(x, "\\d+")))
    tp_lab <- case_when(is_uricult ~ "Uricult", !is.na(tp_num) ~ paste0("T", tp_num), TRUE ~ "Unscheduled")
    tibble(tp_lab = factor(tp_lab, levels = c(paste0("T", sort(unique(tp_num[!is.na(tp_num)]))), "Uricult", "Unscheduled")))
  }

  st_persist <- mlst_tbl %>%
    select(Participant_id, Timepoint, ST, mlst_complete) %>%
    bind_cols(tp_norm(.$Timepoint)) %>%
    distinct(Participant_id, tp_lab, ST, mlst_complete) %>%
    group_by(Participant_id) %>%
    summarise(
      n_tp = n(),
      n_ST = n_distinct(ST),
      dominant_ST = names(sort(table(ST), decreasing = TRUE))[1],
      frac_domin = as.numeric(max(table(ST))) / n_tp,
      .groups = "drop"
    ) %>%
    arrange(desc(frac_domin))

  write_csv(st_persist, file.path(DIR_MLST, "ST_persistence_by_participant.csv"))

  # [REPRO] Calculate consecutive-timepoint concordance (more rigorous persistence metric)
  # Sort by participant and timepoint, then compare ST to previous timepoint
  st_concordance <- mlst_tbl %>%
    select(Participant_id, Timepoint, ST) %>%
    bind_cols(tp_norm(.$Timepoint)) %>%
    arrange(Participant_id, tp_lab) %>%
    group_by(Participant_id) %>%
    mutate(
      prev_ST = lag(ST),
      is_consecutive = !is.na(prev_ST),
      same_as_prev = ST == prev_ST
    ) %>%
    filter(is_consecutive) %>%
    ungroup() %>%
    summarise(
      n_transitions = n(),
      n_same = sum(same_as_prev, na.rm = TRUE),
      pct_concordance = n_same / n_transitions
    )

  msg(
    "Consecutive timepoint ST concordance: %.1f%% (%d/%d transitions)",
    st_concordance$pct_concordance * 100, st_concordance$n_same, st_concordance$n_transitions
  )

  write_csv(st_concordance, file.path(DIR_MLST, "ST_consecutive_concordance.csv"))

  g <- ggplot(st_persist, aes(reorder(Participant_id, frac_domin), frac_domin)) +
    geom_col() +
    coord_flip() +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
    labs(x = "Participant", y = "Fraction with dominant ST", title = "Within-Host Sequence Type Persistence") +
    theme_minimal()

  ggsave(file.path(DIR_PLOTS_MLST, "st_persistence_by_participant.png"), g, width = 7, height = 8)
}

# 8. Write Outputs
# ------------------------------------------------------------------------------
write_tsv(mlst_tbl, file.path(DIR_MLST, "mlst_all.tsv"))
write_csv(mlst_tbl, file.path(DIR_MLST, "mlst_matrix.csv"))
write_csv(summary_tbl, file.path(DIR_MLST, "log_summary.csv"))

msg("✓ MLST analysis complete.")

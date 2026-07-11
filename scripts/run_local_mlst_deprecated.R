#!/usr/bin/env Rscript
# ==============================================================================
# run_local_mlst_deprecated.R
# ==============================================================================
#
# GOAL:
#   Deprecated local PubMLST/mlst runner. This preserves the older local MLST
#   workflow for provenance and sensitivity checks only.
#
# WHY THIS SCRIPT EXISTS:
#   RIVM/provider SeqSphere MLST is now the active chromosomal ST source.
#   This script is retained so local mlst outputs can still be regenerated as
#   non-authoritative provenance and explicit fallback evidence.
#
# INPUTS:
#   - assembly_metadata.csv              (isolate-level metadata)
#   - data/assemblies/*.fasta            (assemblies to type)
#
# OUTPUTS:
#   - results/mlst/mlst_all.tsv          (deprecated local full MLST results)
#   - results/mlst/mlst_with_meta.csv    (deprecated local canonical table)
#   - results/mlst/mlst_matrix.csv       (deprecated local participant × timepoint table)
#   - results/mlst/mlst_qc_summary.csv   (typing quality summary)
#   - results/mlst/top_STs.csv           (most frequent STs)
# ==============================================================================
#
# Biological/Statistical purpose:
#   - Keeps local ST calls available for comparison/fallback only.
#   - Active analyses must use results/mlst/mlst_provider_preferred.csv.
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
invisible(check_tool("mlst"))
invisible(check_tool("blastn"))

# Output directories
DIR_MLST <- file.path(DIR_RESULTS, "mlst")
DIR_MLST_RAW <- file.path(DIR_MLST, "raw")
DIR_MLST_LOG <- file.path(DIR_MLST, "logs")
ensure_dir(DIR_MLST)
ensure_dir(DIR_MLST_RAW)
ensure_dir(DIR_MLST_LOG)
ensure_dir(DIR_PLOTS_MLST)

msg("DEPRECATED local MLST runner: outputs are provenance/fallback only, not the active analysis ST source.")

# 3. Load Assembly Metadata
# ------------------------------------------------------------------------------
if (!file.exists(FILE_METADATA)) stop("Missing ", FILE_METADATA)
assembly_df <- read_csv(FILE_METADATA, show_col_types = FALSE) %>%
  apply_manual_sample_curation(context = "mlst_metadata")

curation_excluded <- assembly_df %>%
  filter(!(analysis_include_primary %in% TRUE) | !(genomics_expected_include %in% TRUE))
if (nrow(curation_excluded) > 0) {
  write_csv(curation_excluded, file.path(DIR_QC, "mlst_manual_curation_excluded_rows.csv"))
  msg("Manual curation excludes ", nrow(curation_excluded), " metadata row(s) from active MLST denominators.")
}
assembly_df <- filter_primary_genomics(assembly_df)

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
assembly_df$found <- !is.na(assembly_df$full_path) & file.exists(assembly_df$full_path)
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

  dat <- suppressWarnings(read_csv(cache, col_types = cols(.default = col_character()), na = c("", "?"), progress = FALSE, show_col_types = FALSE)) %>%
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
  dplyr::count(.status, name = "n") %>%
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
mlst_tbl <- mlst_tbl %>%
  mutate(
    Participant_id = if ("Participant_id" %in% names(.)) as.character(Participant_id) else NA_character_,
    tp_lab = if ("tp_lab" %in% names(.)) normalise_timepoint_preserve_events(tp_lab) else if ("Timepoint" %in% names(.)) normalise_timepoint_preserve_events(Timepoint) else NA_character_,
    Event_type = if ("Event_type" %in% names(.)) as.character(Event_type) else episode_event_type(tp_lab),
    Collection_Date = if ("Collection_Date" %in% names(.)) as.character(Collection_Date) else NA_character_,
    Episode_ID = if ("Episode_ID" %in% names(.)) as.character(Episode_ID) else build_episode_id(., timepoint_col = "tp_lab", event_col = "Event_type", date_col = "Collection_Date")
  )

# QC Flags
# We evaluate typing completeness because missing or ambiguous loci can result in
# non-typable or inaccurate ST assignments. This helps downstream scripts decide 
# whether to drop an isolate from lineage-specific analyses.
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

  mlst_qc <- mlst_tbl %>% dplyr::count(mlst_complete, has_new_allele, ambiguous_call, name = "n")
  write_csv(mlst_qc, file.path(DIR_MLST, "mlst_qc_summary.csv"))
}

# 6. Top STs
# ------------------------------------------------------------------------------
top_STs <- mlst_tbl %>%
  mutate(ST = as.character(ST)) %>%
  dplyr::count(ST, sort = TRUE, name = "n_isolates") %>%
  mutate(pct = n_isolates / sum(n_isolates))
write_csv(top_STs, file.path(DIR_MLST, "top_STs.csv"))

# Canonical participant-timepoint MLST table for downstream episode-level joins.
mlst_episode_tbl <- mlst_tbl
canonical_file <- file.path(DIR_QC, "canonical_assembly_selection.csv")
if (file.exists(canonical_file) && "full_path" %in% names(mlst_episode_tbl)) {
  canonical_paths <- read_csv(canonical_file, show_col_types = FALSE) %>%
    filter(selected_canonical %in% TRUE) %>%
    mutate(full_path = normalizePath(full_path, winslash = "/", mustWork = FALSE)) %>%
    pull(full_path)
  mlst_episode_tbl <- mlst_episode_tbl %>%
    mutate(full_path = normalizePath(full_path, winslash = "/", mustWork = FALSE)) %>%
    filter(full_path %in% canonical_paths)
  msg("Canonical MLST table: ", nrow(mlst_episode_tbl), " selected assembly rows.")
} else {
  mlst_conflicts <- mlst_episode_tbl %>%
    group_by(Participant_id, tp_lab) %>%
    summarise(n_ST = n_distinct(ST[!is.na(ST)]), ST_values = paste(sort(unique(na.omit(ST))), collapse = ";"), .groups = "drop") %>%
    filter(n_ST > 1)
  if (nrow(mlst_conflicts) > 0) {
    write_csv(mlst_conflicts, file.path(DIR_QC, "mlst_duplicate_participant_timepoint_st_conflicts.csv"))
    msg("WARNING: canonical selection missing; wrote ", nrow(mlst_conflicts), " ST conflict(s).")
  }
  mlst_episode_tbl <- mlst_episode_tbl %>%
    group_by(Participant_id, tp_lab) %>%
    summarise(across(everything(), ~ first(.x)), .groups = "drop")
}

append_denominator_summary(
  mlst_tbl,
  "06_MLST.R",
  "mlst_all_assemblies",
  "assembly",
  FILE_METADATA,
  "Assembly-level MLST; assembler alternatives are retained here"
)
append_denominator_summary(
  mlst_episode_tbl,
  "06_MLST.R",
  "mlst_local_canonical_episode_table",
  "participant_timepoint",
  FILE_MLST_LOCAL_CANONICAL,
  "Local mlst pipeline canonical selected assemblies; provider-preferred integration owns FILE_MLST_CANONICAL"
)

# 7. Persistence (if applicable)
# ------------------------------------------------------------------------------
if ("Participant_id" %in% names(mlst_episode_tbl) && "tp_lab" %in% names(mlst_episode_tbl)) {
  tp_norm <- function(x) {
    tp_lab <- normalise_timepoint_preserve_events(x)
    tp_num <- suppressWarnings(readr::parse_number(tp_lab))
    tibble(tp_lab = factor(tp_lab, levels = c(paste0("T", sort(unique(tp_num[!is.na(tp_num)]))), "Uricult")))
  }

  st_persist <- mlst_episode_tbl %>%
    select(Participant_id, tp_lab, ST, mlst_complete) %>%
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

  # [REPRO] Calculate consecutive-timepoint concordance
  # This is a more rigorous persistence metric than simple "dominant ST".
  # It asks: if we look at adjacent timepoints for the same participant,
  # does the ST change? This defines "strain replacement" vs "persistence".
  st_concordance <- mlst_episode_tbl %>%
    select(Participant_id, tp_lab, ST) %>%
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

  msg(sprintf(
    "Consecutive timepoint ST concordance: %.1f%% (%d/%d transitions)",
    st_concordance$pct_concordance * 100, st_concordance$n_same, st_concordance$n_transitions
  ))

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
write_tsv(mlst_tbl, FILE_MLST_LOCAL_ALL)
write_csv(mlst_episode_tbl, FILE_MLST_MATRIX)
write_csv(mlst_episode_tbl, FILE_MLST_LOCAL_CANONICAL)
write_csv(summary_tbl, file.path(DIR_MLST, "log_summary.csv"))

write_uti_attrition_outputs()

msg("✓ MLST analysis complete.")

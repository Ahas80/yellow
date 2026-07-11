#!/usr/bin/env Rscript
# ==============================================================================
# 12d_wgs_badsize_rescue_screen.R
# ------------------------------------------------------------------------------
# Optional sensitivity screen for assemblies excluded from canonical WGS/VF use
# only because total assembly size was outside the primary E. coli QC window.
#
# This script does not change canonical QC. It creates a separate sensitivity
# selection file that keeps all primary canonical assemblies and adds only
# bad-size assemblies that pass explicit rescue evidence checks.
#
# Intended pipeline use:
#   Rscript scripts/12d_wgs_badsize_rescue_screen.R
#   Rscript 02_gene_presence_analysis.R --selection_file results/sensitivity/rescue_qc/rescue_assembly_selection.csv --selection_column selected_for_vf_sensitivity --out_dir results/sensitivity/vf_rescue --out_suffix rescue
# ==============================================================================

source("00_config.R")
source("R/wgs_helpers.R")
source("R/pipeline_qc_helpers.R")

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
  library(stringr)
  library(purrr)
  library(tibble)
  library(optparse)
})

QC_CONFIG <- get_qc_config()

option_list <- list(
  make_option(c("--canonical_file"),
    type = "character", default = file.path(DIR_QC, "canonical_assembly_selection.csv"),
    help = "Canonical assembly selection CSV from 12a_wgs_qc.R [default: %default]"
  ),
  make_option(c("--out_dir"),
    type = "character", default = file.path(DIR_RESULTS, "sensitivity", "rescue_qc"),
    help = "Output directory for rescue QC files [default: %default]"
  ),
  make_option(c("--selection_out"),
    type = "character", default = NA_character_,
    help = "Sensitivity selection CSV [default: <out_dir>/rescue_assembly_selection.csv]"
  ),
  make_option(c("--threads"),
    type = "integer", default = max(1L, min(6L, CORES_USE)),
    help = "Threads for external tools [default: %default]"
  ),
  make_option(c("--min_genome_size"),
    type = "double", default = QC_CONFIG$MIN_GENOME_SIZE,
    help = "Primary minimum genome size [default: %default]"
  ),
  make_option(c("--primary_max_genome_size"),
    type = "double", default = QC_CONFIG$MAX_GENOME_SIZE,
    help = "Primary maximum genome size [default: %default]"
  ),
  make_option(c("--max_rescue_size"),
    type = "double", default = 6.8e6,
    help = "Upper size accepted for borderline/cleaned rescue candidates [default: %default]"
  ),
  make_option(c("--extreme_large_size"),
    type = "double", default = 9.0e6,
    help = "Boundary for extreme-large assemblies [default: %default]"
  ),
  make_option(c("--min_completeness"),
    type = "double", default = QC_CONFIG$MIN_COMPLETENESS,
    help = "Minimum CheckM2 completeness for rescue [default: %default]"
  ),
  make_option(c("--max_contamination"),
    type = "double", default = QC_CONFIG$MAX_CONTAMINATION,
    help = "Maximum CheckM2 contamination for rescue [default: %default]"
  ),
  make_option(c("--max_mash_distance"),
    type = "double", default = 0.10,
    help = "Maximum Mash distance before identity is treated as discordant [default: %default]"
  ),
  make_option(c("--target_taxid"),
    type = "character", default = "562",
    help = "Target taxid passed to FCS-GX when configured [default: %default]"
  ),
  make_option(c("--checkm2_bin"),
    type = "character", default = Sys.getenv("CHECKM2_BIN", unset = NA_character_),
    help = "Path/name for checkm2 [default: PATH or CHECKM2_BIN]"
  ),
  make_option(c("--gunc_bin"),
    type = "character", default = Sys.getenv("GUNC_BIN", unset = NA_character_),
    help = "Path/name for gunc [default: PATH or GUNC_BIN]"
  ),
  make_option(c("--gunc_db"),
    type = "character", default = Sys.getenv("GUNC_DB", unset = NA_character_),
    help = "GUNC database file [default: GUNC_DB]"
  ),
  make_option(c("--quast_bin"),
    type = "character", default = Sys.getenv("QUAST_BIN", unset = NA_character_),
    help = "Path/name for quast.py [default: PATH or QUAST_BIN]"
  ),
  make_option(c("--mlst_bin"),
    type = "character", default = Sys.getenv("MLST_BIN", unset = NA_character_),
    help = "Path/name for mlst [default: PATH or MLST_BIN]"
  ),
  make_option(c("--mash_bin"),
    type = "character", default = Sys.getenv("MASH_BIN", unset = NA_character_),
    help = "Path/name for mash [default: PATH or MASH_BIN]"
  ),
  make_option(c("--fcsgx_bin"),
    type = "character", default = Sys.getenv("FCSGX_BIN", unset = NA_character_),
    help = "Path/name for FCS-GX runner, for example fcs.py [default: FCSGX_BIN]"
  ),
  make_option(c("--fcsgx_db"),
    type = "character", default = Sys.getenv("FCSGX_DB", unset = NA_character_),
    help = "FCS-GX database directory [default: FCSGX_DB]"
  ),
  make_option(c("--manifest_only"),
    action = "store_true", default = FALSE,
    help = "Write manifest/selection skeleton without running external tools"
  )
)
opt <- parse_args(OptionParser(option_list = option_list))

selection_out <- if (!is.na(opt$selection_out) && nzchar(opt$selection_out)) {
  opt$selection_out
} else {
  file.path(opt$out_dir, "rescue_assembly_selection.csv")
}

ensure_dir(opt$out_dir)
ensure_dir(dirname(selection_out))
log_dir <- file.path(opt$out_dir, "logs")
original_dir <- file.path(opt$out_dir, "candidate_fastas_original")
analysis_dir <- file.path(opt$out_dir, "candidate_fastas_screened")
fcsgx_dir <- file.path(opt$out_dir, "fcsgx")
tool_dir <- file.path(opt$out_dir, "tool_outputs")
walk(c(log_dir, original_dir, analysis_dir, fcsgx_dir, tool_dir), ensure_dir)

out_manifest <- file.path(opt$out_dir, "badsize_candidate_manifest.csv")
out_tools <- file.path(opt$out_dir, "badsize_rescue_tool_results.csv")
out_availability <- file.path(opt$out_dir, "tool_availability.csv")
out_report <- file.path(opt$out_dir, "rescue_decision_report.md")

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
  suppressWarnings(as.numeric(gsub("[^0-9.eE+-]", "", as.character(x))))
}

collapse_unique <- function(x) {
  x <- unique(as.character(x[!is.na(x) & nzchar(as.character(x))]))
  if (length(x) == 0) NA_character_ else paste(sort(x), collapse = ";")
}

strip_fasta_ext <- function(x) {
  sub("\\.(fa|fasta|fna)(\\.gz)?$", "", basename(as.character(x)), ignore.case = TRUE)
}

stage_ext <- function(x) {
  if (!is.na(x) && nzchar(x) && str_detect(x, regex("\\.gz$", ignore_case = TRUE))) ".fasta.gz" else ".fasta"
}

first_existing_col <- function(df, patterns) {
  nms <- names(df)
  for (pat in patterns) {
    hit <- nms[grepl(pat, nms, ignore.case = TRUE)]
    if (length(hit) > 0) return(hit[1])
  }
  NA_character_
}

resolve_bin <- function(preferred, candidates) {
  raw <- unique(c(preferred, candidates))
  raw <- raw[!is.na(raw) & nzchar(raw)]
  for (x in raw) {
    if (grepl("/", x) && file.exists(x)) return(normalizePath(x, winslash = "/", mustWork = FALSE))
    hit <- Sys.which(x)
    if (nzchar(hit)) return(unname(hit))
  }

  conda_dirs <- c(
    file.path(Sys.getenv("CONDA_PREFIX"), "bin"),
    file.path(Sys.getenv("HOME"), "miniconda3/envs/asm-snp-x86/bin"),
    file.path(Sys.getenv("HOME"), "miniconda3/envs/wgs/bin"),
    file.path(Sys.getenv("HOME"), "miniconda3/envs/panaroo/bin")
  )
  for (d in conda_dirs) {
    for (x in candidates) {
      p <- file.path(d, x)
      if (file.exists(p)) return(normalizePath(p, winslash = "/", mustWork = FALSE))
    }
  }
  ""
}

bin_version <- function(bin, args = c("--version")) {
  if (!nzchar(bin)) return(NA_character_)
  out <- tryCatch(system2(bin, args, stdout = TRUE, stderr = TRUE), error = function(e) character())
  if (length(out) == 0) return(NA_character_)
  paste(head(out, 2), collapse = " | ")
}

run_logged <- function(cmd, args, name, stdout = NULL) {
  log_file <- file.path(log_dir, paste0(name, ".log"))
  status <- tryCatch(
    {
      system2(cmd, args = args, stdout = if (is.null(stdout)) log_file else stdout, stderr = log_file)
    },
    error = function(e) {
      writeLines(conditionMessage(e), log_file)
      127L
    }
  )
  status <- as.integer(if (is.null(status)) 0L else status)
  tibble(tool = name, status = status, log_file = log_file)
}

classify_size <- function(total_bp) {
  dplyr::case_when(
    is.na(total_bp) ~ "unknown_size",
    total_bp < opt$min_genome_size ~ "too_small",
    total_bp <= opt$primary_max_genome_size ~ "within_primary_limit",
    total_bp <= opt$max_rescue_size ~ "borderline_large",
    total_bp < opt$extreme_large_size ~ "large",
    TRUE ~ "extreme_large"
  )
}

stage_link <- function(src, dest) {
  if (is.na(src) || !nzchar(src) || !file.exists(src)) return(NA_character_)
  if (file.exists(dest) || file.info(dest)$isdir %in% TRUE) unlink(dest, recursive = TRUE, force = TRUE)
  ok <- suppressWarnings(file.symlink(normalizePath(src, winslash = "/", mustWork = FALSE), dest))
  if (!isTRUE(ok)) ok <- suppressWarnings(file.copy(src, dest, overwrite = TRUE))
  if (!isTRUE(ok)) return(NA_character_)
  normalizePath(dest, winslash = "/", mustWork = FALSE)
}

find_cleaned_fasta <- function(sample_dir) {
  fps <- list.files(sample_dir, pattern = "\\.(fa|fasta|fna)(\\.gz)?$", recursive = TRUE, full.names = TRUE, ignore.case = TRUE)
  if (length(fps) == 0) return(NA_character_)
  preferred <- fps[str_detect(tolower(basename(fps)), "clean|filtered|decontam")]
  normalizePath(if (length(preferred) > 0) preferred[1] else fps[1], winslash = "/", mustWork = FALSE)
}

match_candidate <- function(x, lookup) {
  key <- strip_fasta_ext(x)
  idx <- match(key, lookup$analysis_stem)
  if (!is.na(idx)) return(lookup$candidate_id[idx])
  idx <- match(basename(as.character(x)), basename(lookup$analysis_staged_path))
  if (!is.na(idx)) return(lookup$candidate_id[idx])
  NA_character_
}

if (!file.exists(opt$canonical_file)) {
  stop("Missing canonical selection file: ", opt$canonical_file, ". Run 12a_wgs_qc.R first.")
}

canonical <- read_csv(opt$canonical_file, show_col_types = FALSE) %>%
  mutate(.row_id = row_number())

if (!"full_path" %in% names(canonical) && "fasta_path" %in% names(canonical)) {
  canonical$full_path <- canonical$fasta_path
}
if (!"full_path" %in% names(canonical)) {
  stop("Canonical selection lacks full_path/fasta_path column: ", opt$canonical_file)
}
if (!"QC_REASON" %in% names(canonical)) canonical$QC_REASON <- NA_character_
if (!"QC_PASS" %in% names(canonical)) canonical$QC_PASS <- FALSE
if (!"selected_canonical" %in% names(canonical)) canonical$selected_canonical <- FALSE
if (!"file_exists" %in% names(canonical)) canonical$file_exists <- file.exists(canonical$full_path)
if (!"total_bp" %in% names(canonical) && "total_bases" %in% names(canonical)) canonical$total_bp <- canonical$total_bases
if (!"total_bp" %in% names(canonical)) canonical$total_bp <- NA_real_
if (!"Participant_id" %in% names(canonical)) canonical$Participant_id <- NA_character_
if (!"tp_lab" %in% names(canonical)) canonical$tp_lab <- if ("Timepoint" %in% names(canonical)) canonical$Timepoint else NA_character_
if (!"Assembly_ID" %in% names(canonical)) canonical$Assembly_ID <- tools::file_path_sans_ext(basename(canonical$full_path))

canonical <- canonical %>%
  mutate(
    Participant_id = as.character(Participant_id),
    tp_lab = normalise_timepoint_preserve_events(tp_lab),
    QC_REASON = as.character(QC_REASON),
    QC_PASS = as_bool(QC_PASS),
    selected_canonical = as_bool(selected_canonical),
    file_exists = as_bool(file_exists),
    file_exists = if_else(is.na(file_exists), !is.na(full_path) & file.exists(full_path), file_exists),
    total_bp = parse_num(total_bp),
    original_full_path = full_path
  )

status_for_join <- if (file.exists(FILE_STATUS_MAP)) {
  read_csv(FILE_STATUS_MAP, show_col_types = FALSE) %>%
    prefer_primary_uti_status() %>%
    mutate(
      Participant_id = as.character(Participant_id),
      tp_lab = if ("tp_lab" %in% names(.)) normalise_timepoint_preserve_events(tp_lab) else normalise_timepoint_preserve_events(Timepoint)
    ) %>%
    select(any_of(c("Participant_id", "tp_lab", "Infection_Status", "UTI_Status", "UTI_binary",
                    "Not_UTI_subgroup", "Infection_Status_legacy", "Episode_ID",
                    "Collection_Date", "Urine_collection_method"))) %>%
    distinct(Participant_id, tp_lab, .keep_all = TRUE)
} else {
  tibble(Participant_id = character(), tp_lab = character())
}

primary_key_status <- canonical %>%
  group_by(Participant_id, tp_lab) %>%
  summarise(
    has_primary_selected_for_key = any(selected_canonical %in% TRUE, na.rm = TRUE),
    has_any_qc_pass_for_key = any(QC_PASS %in% TRUE, na.rm = TRUE),
    .groups = "drop"
  )

bad_candidates <- canonical %>%
  filter(str_detect(QC_REASON, fixed("BadSize"))) %>%
  left_join(primary_key_status, by = c("Participant_id", "tp_lab")) %>%
  left_join(status_for_join, by = c("Participant_id", "tp_lab")) %>%
  mutate(
    candidate_id_source = coalesce(as.character(Assembly_ID), tools::file_path_sans_ext(basename(full_path)), paste0("row_", .row_id)),
    candidate_id = paste0(sprintf("%03d", row_number()), "_", sanitize_id(candidate_id_source)),
    rescue_size_class = classify_size(total_bp),
    original_staged_path = map2_chr(
      full_path,
      candidate_id,
      ~ {
        ext <- stage_ext(.x)
        stage_link(.x, file.path(original_dir, paste0(.y, ext)))
      }
    )
  )

write_csv(bad_candidates, out_manifest)

bins <- tibble(
  tool = c("quast", "checkm2", "gunc", "mlst", "mash", "fcsgx"),
  path = c(
    resolve_bin(opt$quast_bin, c("quast.py", "quast")),
    resolve_bin(opt$checkm2_bin, c("checkm2")),
    resolve_bin(opt$gunc_bin, c("gunc")),
    resolve_bin(opt$mlst_bin, c("mlst")),
    resolve_bin(opt$mash_bin, c("mash")),
    resolve_bin(opt$fcsgx_bin, c("fcs.py", "fcs-gx", "fcs"))
  )
) %>%
  mutate(
    available = nzchar(path),
    version = map_chr(path, bin_version),
    note = case_when(
      tool == "gunc" & (!available | is.na(opt$gunc_db) | !nzchar(opt$gunc_db) | !file.exists(opt$gunc_db)) ~ "GUNC requires --gunc_db or GUNC_DB",
      tool == "fcsgx" & (!available | is.na(opt$fcsgx_db) | !nzchar(opt$fcsgx_db) | !file.exists(opt$fcsgx_db)) ~ "FCS-GX optional; requires --fcsgx_bin/FCSGX_BIN and --fcsgx_db/FCSGX_DB",
      TRUE ~ NA_character_
    ),
    available = case_when(
      tool == "gunc" ~ available & !is.na(opt$gunc_db) & nzchar(opt$gunc_db) & file.exists(opt$gunc_db),
      tool == "fcsgx" ~ available & !is.na(opt$fcsgx_db) & nzchar(opt$fcsgx_db) & file.exists(opt$fcsgx_db),
      TRUE ~ available
    )
  )
write_csv(bins, out_availability)

get_bin <- function(tool) {
  bins$path[match(tool, bins$tool)]
}
has_tool <- function(tool) {
  isTRUE(bins$available[match(tool, bins$tool)])
}

required_tools_available <- has_tool("checkm2") && has_tool("gunc")

fcsgx_results <- tibble(candidate_id = character(), fcsgx_status = character(), fcsgx_log = character(), cleaned_fasta_path = character())
if (!opt$manifest_only && nrow(bad_candidates) > 0 && has_tool("fcsgx")) {
  fcsgx_results <- bad_candidates %>%
    filter(!is.na(original_staged_path), file.exists(original_staged_path)) %>%
    mutate(res = pmap(
      list(candidate_id, original_staged_path),
      function(candidate_id, original_staged_path) {
        sample_out <- file.path(fcsgx_dir, candidate_id)
        ensure_dir(sample_out)
        fcsgx_bin <- get_bin("fcsgx")
        base_args <- c("screen", "genome", "--fasta", original_staged_path, "--out-dir", sample_out, "--gx-db", opt$fcsgx_db, "--tax-id", opt$target_taxid)
        if (str_detect(basename(fcsgx_bin), regex("\\.py$", ignore_case = TRUE))) {
          py <- resolve_bin(NA_character_, c("python3", "python"))
          if (!nzchar(py)) {
            return(tibble(candidate_id = candidate_id, fcsgx_status = "python_not_found", fcsgx_log = NA_character_, cleaned_fasta_path = NA_character_))
          }
          run <- run_logged(py, c(fcsgx_bin, base_args), paste0("fcsgx_", candidate_id))
        } else {
          run <- run_logged(fcsgx_bin, base_args, paste0("fcsgx_", candidate_id))
        }
        cleaned <- find_cleaned_fasta(sample_out)
        tibble(
          candidate_id = candidate_id,
          fcsgx_status = if_else(run$status == 0L, "completed", paste0("failed_status_", run$status)),
          fcsgx_log = run$log_file,
          cleaned_fasta_path = cleaned
        )
      }
    )) %>%
    select(res) %>%
    unnest(res)
}

candidates <- bad_candidates %>%
  left_join(fcsgx_results, by = "candidate_id") %>%
  mutate(
    fcsgx_status = coalesce(fcsgx_status, if_else(has_tool("fcsgx"), "not_run", "not_configured")),
    used_cleaned_fasta = !is.na(cleaned_fasta_path) & file.exists(cleaned_fasta_path),
    screened_source_path = if_else(used_cleaned_fasta, cleaned_fasta_path, original_staged_path),
    analysis_staged_path = map2_chr(
      screened_source_path,
      candidate_id,
      ~ {
        ext <- stage_ext(.x)
        stage_link(.x, file.path(analysis_dir, paste0(.y, ext)))
      }
    )
  )

analysis_stats <- candidates %>%
  transmute(candidate_id, stats = map(analysis_staged_path, fasta_stats)) %>%
  unnest(stats) %>%
  rename(
    analysis_n_contigs = n_contigs,
    analysis_n50 = n50,
    analysis_total_bp = total_bp
  )

candidates <- candidates %>%
  left_join(analysis_stats, by = "candidate_id") %>%
  mutate(analysis_size_class = classify_size(analysis_total_bp))

lookup <- candidates %>%
  transmute(
    candidate_id,
    analysis_staged_path,
    analysis_stem = strip_fasta_ext(analysis_staged_path)
  )
analysis_paths <- lookup$analysis_staged_path[!is.na(lookup$analysis_staged_path) & file.exists(lookup$analysis_staged_path)]

tool_runs <- tibble(tool = character(), status = integer(), log_file = character())

quast_tbl <- tibble(candidate_id = character(), quast_total_length = double(), quast_n_contigs = double(), quast_n50 = double())
if (!opt$manifest_only && length(analysis_paths) > 0 && has_tool("quast")) {
  qdir <- file.path(tool_dir, "quast")
  ensure_dir(qdir)
  tool_runs <- bind_rows(tool_runs, run_logged(get_bin("quast"), c("-o", qdir, "-t", as.character(opt$threads), analysis_paths), "quast"))
  qfile <- file.path(qdir, "report.tsv")
  if (file.exists(qfile)) {
    qraw <- read_tsv(qfile, show_col_types = FALSE)
    metric_col <- names(qraw)[1]
    qlong <- qraw %>%
      pivot_longer(cols = -all_of(metric_col), names_to = "quast_name", values_to = "value") %>%
      mutate(metric = .data[[metric_col]])
    quast_tbl <- qlong %>%
      filter(metric %in% c("# contigs", "Total length", "N50")) %>%
      mutate(candidate_id = map_chr(quast_name, match_candidate, lookup = lookup)) %>%
      select(candidate_id, metric, value) %>%
      filter(!is.na(candidate_id)) %>%
      pivot_wider(names_from = metric, values_from = value) %>%
      transmute(
        candidate_id,
        quast_total_length = parse_num(.data[["Total length"]]),
        quast_n_contigs = parse_num(.data[["# contigs"]]),
        quast_n50 = parse_num(.data[["N50"]])
      )
  }
}

checkm2_tbl <- tibble(candidate_id = character(), checkm2_completeness = double(), checkm2_contamination = double(), checkm2_model = character())
if (!opt$manifest_only && length(analysis_paths) > 0 && has_tool("checkm2")) {
  cdir <- file.path(tool_dir, "checkm2")
  ensure_dir(cdir)
  tool_runs <- bind_rows(tool_runs, run_logged(get_bin("checkm2"), c("predict", "--input", analysis_dir, "--output-directory", cdir, "--threads", as.character(opt$threads)), "checkm2"))
  cfile <- file.path(cdir, "quality_report.tsv")
  if (file.exists(cfile)) {
    craw <- read_tsv(cfile, show_col_types = FALSE)
    name_col <- first_existing_col(craw, c("^Name$", "^Bin Id$", "genome", "file"))
    comp_col <- first_existing_col(craw, c("^Completeness$", "completeness"))
    contam_col <- first_existing_col(craw, c("^Contamination$", "contamination"))
    model_col <- first_existing_col(craw, c("model"))
    if (!is.na(name_col) && !is.na(comp_col) && !is.na(contam_col)) {
      checkm2_tbl <- craw %>%
        transmute(
          candidate_id = map_chr(.data[[name_col]], match_candidate, lookup = lookup),
          checkm2_completeness = parse_num(.data[[comp_col]]),
          checkm2_contamination = parse_num(.data[[contam_col]]),
          checkm2_model = if (!is.na(model_col)) as.character(.data[[model_col]]) else NA_character_
        ) %>%
        filter(!is.na(candidate_id))
    }
  }
}

gunc_tbl <- tibble(candidate_id = character(), gunc_pass = logical(), gunc_css = double(), gunc_contamination_portion = double())
if (!opt$manifest_only && length(analysis_paths) > 0 && has_tool("gunc")) {
  gdir <- file.path(tool_dir, "gunc")
  ensure_dir(gdir)
  tool_runs <- bind_rows(tool_runs, run_logged(get_bin("gunc"), c("run", "--input_dir", analysis_dir, "--db_file", opt$gunc_db, "--out_dir", gdir, "--threads", as.character(opt$threads)), "gunc"))
  gfiles <- list.files(gdir, pattern = "maxCSS.*\\.tsv$", recursive = TRUE, full.names = TRUE)
  if (length(gfiles) > 0) {
    graw <- read_tsv(gfiles[1], show_col_types = FALSE)
    name_col <- first_existing_col(graw, c("^genome$", "^name$", "file"))
    pass_col <- first_existing_col(graw, c("pass"))
    css_col <- first_existing_col(graw, c("clade_separation_score", "css"))
    contam_col <- first_existing_col(graw, c("contamination_portion", "contamination"))
    if (!is.na(name_col)) {
      pass_values <- if (!is.na(pass_col)) {
        as_bool(graw[[pass_col]])
      } else {
        rep(NA, nrow(graw))
      }
      gunc_tbl <- graw %>%
        transmute(
          candidate_id = map_chr(.data[[name_col]], match_candidate, lookup = lookup),
          gunc_pass = pass_values,
          gunc_css = if (!is.na(css_col)) parse_num(.data[[css_col]]) else NA_real_,
          gunc_contamination_portion = if (!is.na(contam_col)) parse_num(.data[[contam_col]]) else NA_real_
        ) %>%
        filter(!is.na(candidate_id))
    }
  }
}

mlst_tbl <- tibble(candidate_id = character(), mlst_scheme = character(), mlst_st = character())
if (!opt$manifest_only && length(analysis_paths) > 0 && has_tool("mlst")) {
  mlst_out <- file.path(tool_dir, "mlst.tsv")
  tool_runs <- bind_rows(tool_runs, run_logged(get_bin("mlst"), analysis_paths, "mlst", stdout = mlst_out))
  if (file.exists(mlst_out) && file.size(mlst_out) > 0) {
    mraw <- read_tsv(mlst_out, col_names = FALSE, show_col_types = FALSE)
    if (ncol(mraw) >= 3) {
      mlst_tbl <- mraw %>%
        transmute(
          candidate_id = map_chr(X1, match_candidate, lookup = lookup),
          mlst_scheme = as.character(X2),
          mlst_st = as.character(X3)
        ) %>%
        filter(!is.na(candidate_id)) %>%
        group_by(candidate_id) %>%
        summarise(
          mlst_scheme = collapse_unique(mlst_scheme),
          mlst_st = collapse_unique(mlst_st),
          .groups = "drop"
        )
    }
  }
}

mash_tbl <- tibble(candidate_id = character(), mash_min_distance = double(), mash_best_reference = character())
if (!opt$manifest_only && length(analysis_paths) > 0 && has_tool("mash")) {
  trusted_refs <- canonical %>%
    filter(selected_canonical %in% TRUE, file_exists %in% TRUE, !is.na(full_path), file.exists(full_path)) %>%
    pull(full_path) %>%
    unique()
  if (length(trusted_refs) > 0) {
    mash_dir <- file.path(tool_dir, "mash")
    ensure_dir(mash_dir)
    ref_list <- file.path(mash_dir, "trusted_primary_fastas.txt")
    writeLines(trusted_refs, ref_list)
    sketch_prefix <- file.path(mash_dir, "trusted_primary")
    tool_runs <- bind_rows(tool_runs, run_logged(get_bin("mash"), c("sketch", "-p", as.character(opt$threads), "-o", sketch_prefix, "-l", ref_list), "mash_sketch"))
    mash_out <- file.path(mash_dir, "dist_to_trusted_primary.tsv")
    msh_file <- paste0(sketch_prefix, ".msh")
    if (file.exists(msh_file)) {
      tool_runs <- bind_rows(tool_runs, run_logged(get_bin("mash"), c("dist", msh_file, analysis_paths), "mash_dist", stdout = mash_out))
      if (file.exists(mash_out) && file.size(mash_out) > 0) {
        mraw <- read_tsv(mash_out, col_names = c("reference", "query", "distance", "p_value", "shared_hashes"), show_col_types = FALSE)
        mash_tbl <- mraw %>%
          mutate(
            candidate_id = map_chr(query, match_candidate, lookup = lookup),
            distance = parse_num(distance)
          ) %>%
          filter(!is.na(candidate_id)) %>%
          group_by(candidate_id) %>%
          slice_min(distance, n = 1, with_ties = FALSE) %>%
          ungroup() %>%
          transmute(candidate_id, mash_min_distance = distance, mash_best_reference = reference)
      }
    }
  }
}

write_csv(tool_runs, file.path(opt$out_dir, "tool_run_status.csv"))

tool_results <- candidates %>%
  select(
    .row_id, candidate_id, Participant_id, tp_lab, any_of(c("Infection_Status", "Episode_ID")),
    Assembly_ID, Assembler, file_name, original_full_path, original_staged_path,
    cleaned_fasta_path, used_cleaned_fasta, analysis_staged_path,
    total_bp, rescue_size_class, analysis_total_bp, analysis_n_contigs,
    analysis_n50, analysis_size_class, QC_PASS, QC_REASON,
    has_primary_selected_for_key, has_any_qc_pass_for_key, file_exists
  ) %>%
  left_join(quast_tbl, by = "candidate_id") %>%
  left_join(checkm2_tbl, by = "candidate_id") %>%
  left_join(gunc_tbl, by = "candidate_id") %>%
  left_join(mlst_tbl, by = "candidate_id") %>%
  left_join(mash_tbl, by = "candidate_id") %>%
  mutate(
    decision_total_bp = coalesce(quast_total_length, analysis_total_bp, total_bp),
    checkm2_pass = !is.na(checkm2_completeness) &
      !is.na(checkm2_contamination) &
      checkm2_completeness >= opt$min_completeness &
      checkm2_contamination <= opt$max_contamination,
    has_required_evidence = !is.na(checkm2_completeness) &
      !is.na(checkm2_contamination) &
      !is.na(gunc_pass),
    mlst_identity_block = map_lgl(str_split(coalesce(tolower(mlst_scheme), ""), ";"), function(vals) {
      vals <- trimws(vals)
      vals <- vals[!vals %in% c("", "-", "na", "unknown")]
      length(vals) > 0 && any(!str_detect(vals, "ecoli|escherichia"))
    }),
    mash_identity_block = !is.na(mash_min_distance) & mash_min_distance > opt$max_mash_distance,
    identity_block = mlst_identity_block | mash_identity_block,
    rescue_accepted = case_when(
      opt$manifest_only ~ FALSE,
      file_exists %in% FALSE ~ FALSE,
      has_primary_selected_for_key %in% TRUE ~ FALSE,
      !required_tools_available ~ FALSE,
      !has_required_evidence ~ FALSE,
      identity_block ~ FALSE,
      !checkm2_pass ~ FALSE,
      !(gunc_pass %in% TRUE) ~ FALSE,
      rescue_size_class == "too_small" ~ TRUE,
      rescue_size_class == "borderline_large" &
        decision_total_bp <= opt$max_rescue_size ~ TRUE,
      rescue_size_class == "large" &
        used_cleaned_fasta &
        decision_total_bp >= opt$min_genome_size &
        decision_total_bp <= opt$max_rescue_size ~ TRUE,
      rescue_size_class == "extreme_large" &
        used_cleaned_fasta &
        decision_total_bp >= opt$min_genome_size &
        decision_total_bp <= opt$max_rescue_size ~ TRUE,
      TRUE ~ FALSE
    ),
    rescue_decision = case_when(
      opt$manifest_only ~ "not_rescued_manifest_only",
      file_exists %in% FALSE ~ "not_rescued_missing_fasta",
      has_primary_selected_for_key %in% TRUE ~ "not_rescued_primary_qc_pass_exists_for_same_participant_timepoint",
      identity_block ~ "not_rescued_identity_discordance",
      !required_tools_available ~ "not_rescued_missing_required_tools_checkm2_or_gunc",
      !has_required_evidence ~ "not_rescued_missing_required_tool_results",
      !checkm2_pass ~ "not_rescued_checkm2_completeness_or_contamination_fail",
      !(gunc_pass %in% TRUE) ~ "not_rescued_gunc_chimerism_or_taxonomic_discordance_fail",
      rescue_accepted & rescue_size_class == "too_small" ~ "rescued_too_small_checkm2_complete_gunc_clean",
      rescue_accepted & rescue_size_class == "borderline_large" ~ "rescued_borderline_large_checkm2_gunc_clean",
      rescue_accepted & rescue_size_class %in% c("large", "extreme_large") & used_cleaned_fasta ~ "rescued_cleaned_large_checkm2_gunc_clean",
      rescue_size_class %in% c("large", "extreme_large") & !used_cleaned_fasta ~ "not_rescued_large_requires_cleaned_plausible_fasta",
      decision_total_bp > opt$max_rescue_size ~ "not_rescued_still_too_large_after_screening",
      TRUE ~ "not_rescued_no_rule_matched"
    ),
    rescue_reason = case_when(
      rescue_accepted ~ "Accepted only for sensitivity analysis; canonical QC remains unchanged.",
      rescue_decision == "not_rescued_primary_qc_pass_exists_for_same_participant_timepoint" ~ "A QC-passing canonical assembly already represents this participant-timepoint, so adding this FASTA would duplicate the episode rather than rescue a lost episode.",
      rescue_decision == "not_rescued_missing_required_tools_checkm2_or_gunc" ~ "CheckM2 and GUNC are required for automatic rescue decisions; configure tool paths/databases and rerun.",
      rescue_decision == "not_rescued_large_requires_cleaned_plausible_fasta" ~ "Assemblies above the borderline-large window require a cleaned/decontaminated FASTA returning to plausible E. coli size.",
      TRUE ~ rescue_decision
    )
  )

write_csv(tool_results, out_tools)

if (!"Infection_Status" %in% names(tool_results)) tool_results$Infection_Status <- NA_character_

selection_cols <- tool_results %>%
  select(
    .row_id, candidate_id, rescue_accepted, rescue_decision, rescue_reason,
    rescue_size_class, analysis_size_class, original_full_path,
    rescue_full_path = analysis_staged_path, used_cleaned_fasta,
    checkm2_completeness, checkm2_contamination, gunc_pass,
    mlst_scheme, mlst_st, mash_min_distance,
    mlst_identity_block, mash_identity_block, identity_block
  )

rescue_selection <- canonical %>%
  left_join(selection_cols, by = ".row_id") %>%
  mutate(
    rescue_accepted = coalesce(rescue_accepted, FALSE),
    selected_for_vf_sensitivity = selected_canonical %in% TRUE | rescue_accepted,
    rescue_decision = case_when(
      selected_canonical %in% TRUE ~ "primary_canonical",
      !str_detect(QC_REASON, fixed("BadSize")) ~ "not_bad_size_existing_canonical_row",
      TRUE ~ coalesce(rescue_decision, "not_rescued_not_screened")
    ),
    rescue_reason = case_when(
      selected_canonical %in% TRUE ~ "Primary canonical assembly retained from 12a_wgs_qc.R.",
      TRUE ~ rescue_reason
    ),
    full_path = if_else(rescue_accepted & !is.na(rescue_full_path), rescue_full_path, full_path),
    fasta_path = if ("fasta_path" %in% names(.)) if_else(rescue_accepted & !is.na(rescue_full_path), rescue_full_path, fasta_path) else full_path
  ) %>%
  select(-.row_id)

write_csv(rescue_selection, selection_out)

decision_summary <- tool_results %>%
  count(rescue_size_class, rescue_decision, name = "n") %>%
  arrange(rescue_size_class, desc(n), rescue_decision)

status_summary <- tool_results %>%
  mutate(Infection_Status = coalesce(Infection_Status, "Missing/unknown")) %>%
  count(Infection_Status, rescue_size_class, rescue_accepted, name = "n") %>%
  arrange(Infection_Status, rescue_size_class, desc(rescue_accepted))

accepted <- tool_results %>%
  filter(rescue_accepted %in% TRUE) %>%
  arrange(Participant_id, tp_lab, candidate_id)

report <- c(
  "# Bad-size assembly rescue sensitivity screen",
  "",
  sprintf("Generated: %s", format(Sys.time())),
  "",
  "## Rule summary",
  "",
  sprintf("- Canonical size window remains %.1f-%.1f Mb.", opt$min_genome_size / 1e6, opt$primary_max_genome_size / 1e6),
  sprintf("- Borderline-large sensitivity window is >%.1f to %.1f Mb.", opt$primary_max_genome_size / 1e6, opt$max_rescue_size / 1e6),
  sprintf("- Automatic rescue requires CheckM2 completeness >= %.1f%%, CheckM2 contamination <= %.1f%%, and GUNC pass.", opt$min_completeness, opt$max_contamination),
  "- Large/extreme-large assemblies require a cleaned/decontaminated FASTA before they can be rescued.",
  "- Assemblies for participant-timepoints already represented by a primary canonical assembly are not added as rescue rows.",
  "",
  "## Tool availability",
  "",
  paste(capture.output(print(bins, n = Inf, width = Inf)), collapse = "\n"),
  "",
  "## Bad-size decisions",
  "",
  paste(capture.output(print(decision_summary, n = Inf, width = Inf)), collapse = "\n"),
  "",
  "## Status summary",
  "",
  paste(capture.output(print(status_summary, n = Inf, width = Inf)), collapse = "\n"),
  "",
  "## Accepted rescue assemblies",
  "",
  if (nrow(accepted) > 0) {
    paste(capture.output(print(accepted %>% select(any_of(c("Participant_id", "tp_lab", "Infection_Status", "candidate_id", "rescue_size_class", "decision_total_bp", "checkm2_completeness", "checkm2_contamination", "gunc_pass", "rescue_decision"))), n = Inf, width = Inf)), collapse = "\n")
  } else {
    "No bad-size assemblies met the automatic rescue criteria in this run."
  },
  "",
  "## Output files",
  "",
  sprintf("- Candidate manifest: `%s`", out_manifest),
  sprintf("- Tool availability: `%s`", out_availability),
  sprintf("- Tool results and decisions: `%s`", out_tools),
  sprintf("- Sensitivity selection file: `%s`", selection_out)
)
writeLines(report, out_report)

msg("Bad-size candidate manifest: %s", out_manifest)
msg("Tool results and decisions: %s", out_tools)
msg("Sensitivity selection file: %s", selection_out)
msg("Decision report: %s", out_report)

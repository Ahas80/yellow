#!/usr/bin/env Rscript
# ==============================================================================
# 11_compare_strains.R
# ==============================================================================
#
# GOAL:
#   Compute pairwise strain similarity between all within-participant
#   isolates using multiple metrics (ST match, VF Jaccard, plasmid Jaccard,
#   assembly-to-assembly dnadiff SNP differences).  This is the foundation for longitudinal strain
#   tracking: are participants colonised by the same or different strains
#   at each timepoint?
#
# WHY THIS SCRIPT EXISTS:
#   A central question of this project is whether the same bacterial strain
#   persists within a host over time.  If the same ST131 is present during
#   a Not_UTI episode at T0 and UTI at T2, that is biologically very different from a new
#   ST73 arriving to cause UTI.  This script produces the pairwise_metrics
#   table that scripts 15 and 16 use to define strain persistence and
#   phenotype-switch events.
#
# INPUTS:
#   - results/qc/canonical_assembly_selection.csv
#   - results/clinical/status_map.csv
#   - results/mlst/mlst_provider_preferred.csv
#   - results/vf/vf_pa_all.csv
#   - results/plasmids/plasmidfinder_presence_absence.csv
#
# OUTPUTS:
#   - results/strain_compare/pairwise_metrics.csv  (all within-participant pairs)
#   - results/strain_compare/summary_counts.csv
#   - results/strain_compare/summary_by_participant.csv
#   - results/strain_compare/stats_within_vs_between.csv
#   - plots/strain_compare/
# ==============================================================================
#
# Usage:
#   Rscript 11_compare_strains.R --participants P001,P002 --timepoints T0,T1
#
# Biological/Statistical purpose:
#   - Determines if sequential isolates are the "Same" or "Different" strain.
#   - Uses a multi-metric approach (ANI, SNPs, VF content, Plasmid content).
#   - Quantifies strain persistence vs. replacement.
# ==============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
  library(stringr)
  library(purrr)
  library(tibble)
  library(ggplot2)
  library(scales)
})

# Load configuration
source("00_config.R")

# dependencies
if (!requireNamespace("optparse", quietly = TRUE)) {
  stop("Package 'optparse' is required. Install with install.packages('optparse').")
}

source("11_compare_strains_helpers.R")

detected_physical_cores <- suppressWarnings(parallel::detectCores(logical = FALSE))
default_pair_workers <- if (length(detected_physical_cores) != 1L || is.na(detected_physical_cores)) {
  1L
} else {
  max(1L, min(8L, as.integer(detected_physical_cores) - 1L))
}

opt_list <- list(
  optparse::make_option(c("--pairs_csv"),
    type = "character", default = NA,
    help = "CSV with Participant_id_A,Timepoint_A,Participant_id_B,Timepoint_B"
  ),
  optparse::make_option(c("--participants"),
    type = "character", default = NA,
    help = "Comma-separated participants (e.g., P001,P002)"
  ),
  optparse::make_option(c("--timepoints"),
    type = "character", default = NA,
    help = "Comma-separated timepoints (e.g., T0,T1,Uricult)"
  ),
  optparse::make_option(c("--between"),
    action = "store_true", default = FALSE,
    help = "Include between-participant pairs at matching timepoints"
  ),
  optparse::make_option(c("--min_vf_prev"),
    type = "double", default = NA,
    help = "Minimum VF prevalence (0–1) for Jaccard (default: use all)"
  ),
  optparse::make_option(c("--min_inc_prev"),
    type = "double", default = 0,
    help = "Minimum Inc prevalence (0–1) for Jaccard (default: 0)"
  ),
  optparse::make_option(c("--outdir"),
    type = "character", default = DIR_STRAIN,
    help = "Output directory (default: results/strain_compare)"
  ),
  optparse::make_option(c("--id_thresh"),
    type = "double", default = 99.9,
    help = "Identity threshold for 'Same' (default: 99.9)"
  ),
  optparse::make_option(c("--snp_thresh"),
    type = "integer", default = SAME_STRAIN_SNP_THRESHOLD,
    help = sprintf("Operational SNP threshold for 'Same' (default: %d; requires study-specific calibration)", SAME_STRAIN_SNP_THRESHOLD)
  ),
  optparse::make_option(c("--workers"),
    type = "integer", default = default_pair_workers,
    help = "Parallel pair workers on macOS/Linux (default: min(physical cores - 1, 8)); Windows uses one worker"
  )
)

opt <- optparse::parse_args(optparse::OptionParser(option_list = opt_list))

outdir <- opt$outdir
safe_dir_create(outdir)
archive_dir <- archive_legacy_strain_outputs(outdir)
if (!is.null(archive_dir)) timestamp_msg("Archived pre-provenance outputs in ", archive_dir)
# Keep the legacy key-only nucmer_cache untouched.  This versioned directory is
# the only cache read by the SHA-256 provenance workflow.
cache_dir <- file.path(outdir, "dnadiff_cache_sha256_v1")
safe_dir_create(cache_dir)

timestamp_msg("Loading core tables …")
core <- load_core_tables()

# ----- build candidate pair list ---------------------------------------------
parse_csv_pairs <- function(path) {
  x <- readr::read_csv(path, show_col_types = FALSE)
  req <- c("Participant_id_A", "Timepoint_A", "Participant_id_B", "Timepoint_B")
  miss <- setdiff(req, names(x))
  if (length(miss)) stop("pairs_csv missing columns: ", paste(miss, collapse = ", "))
  a <- tp_norm(x$Timepoint_A)
  b <- tp_norm(x$Timepoint_B)
  tibble(
    Participant_id_A = as.character(x$Participant_id_A),
    Timepoint_A = a$tp_lab,
    Participant_id_B = as.character(x$Participant_id_B),
    Timepoint_B = b$tp_lab
  )
}

parse_participants <- function(participants_str, timepoints_str, between = FALSE) {
  available_assemblies <- core$assemblies %>%
    filter(
      selected_canonical %in% TRUE,
      QC_PASS %in% TRUE,
      usable_fasta_path(full_path)
    )

  if (participants_str == "ALL") {
    pids <- available_assemblies$Participant_id %>%
      unique() %>%
      as.character()
  } else {
    pids <- strsplit(participants_str, ",")[[1]] %>%
      trimws() %>%
      unique()
  }
  if (length(pids) < 1) stop("No participants parsed from --participants")

  requested_tps <- NULL
  if (is.na(timepoints_str) || timepoints_str == "") {
    timestamp_msg("Inferring per-participant timepoints from assemblies with usable FASTA paths.")
  } else {
    requested_tps <- strsplit(timepoints_str, ",")[[1]] %>% trimws()
  }
  # within-participant all combinations (unordered)
  within <- map_dfr(pids, function(pid) {
    tps <- requested_tps
    if (is.null(tps)) {
      tps <- available_assemblies %>%
        filter(Participant_id == pid) %>%
        pull(tp_lab) %>%
        as.character() %>%
        unique()
    }
    tps <- tps[!is.na(tps) & nzchar(tps)]
    if (length(tps) < 2) {
      return(tibble())
    }
    comb <- t(combn(tps, 2))
    if (!nrow(comb)) {
      return(tibble())
    }
    tibble(
      Participant_id_A = pid, Timepoint_A = comb[, 1],
      Participant_id_B = pid, Timepoint_B = comb[, 2]
    )
  })
  if (!between) {
    return(within)
  }
  # between-participant pairs at matching timepoints
  between_df <- tibble()
  between_tps <- requested_tps
  if (is.null(between_tps)) {
    between_tps <- available_assemblies %>%
      filter(Participant_id %in% pids) %>%
      pull(tp_lab) %>%
      as.character() %>%
      unique()
  }
  between_tps <- between_tps[!is.na(between_tps) & nzchar(between_tps)]
  if (length(pids) >= 2 && length(between_tps) >= 1) {
    pcomb <- t(combn(pids, 2))
    between_df <- map_dfr(seq_len(nrow(pcomb)), function(i) {
      tibble(
        Participant_id_A = pcomb[i, 1], Timepoint_A = between_tps,
        Participant_id_B = pcomb[i, 2], Timepoint_B = between_tps
      )
    })
  }
  bind_rows(within, between_df)
}

pairs <- NULL
if (!is.na(opt$pairs_csv)) {
  timestamp_msg("Using explicit pairs from ", opt$pairs_csv)
  pairs <- parse_csv_pairs(opt$pairs_csv)
} else if (!is.na(opt$participants)) {
  timestamp_msg("Generating pairs from --participants/--timepoints …")
  pairs <- parse_participants(opt$participants, opt$timepoints, opt$between)
} else {
  stop("Provide --pairs_csv or --participants (with optional --timepoints and --between)")
}

if (!nrow(pairs)) stop("No pairs to compare after parsing input")
pairs <- pairs %>% distinct()

# ----- prevalence filtering for VF / Inc -------------------------------------
vf_wide <- core$vf_pa
vf_meta_cols <- c("Participant_id", "tp_lab", "SampleKey")
vf_gene_cols <- canonical_vf_gene_cols(names(vf_wide), vf_pa_file = FILE_VF_PA, required = FALSE)
if (!is.na(opt$min_vf_prev)) {
  prev <- colMeans((vf_wide[, vf_gene_cols, drop = FALSE] > 0))
  keep <- names(prev)[prev >= opt$min_vf_prev]
  vf_gene_cols <- intersect(vf_gene_cols, keep)
  if (!length(vf_gene_cols)) warning("min_vf_prev filtered out all VF genes; VF Jaccard will be NA")
}

inc_pa <- core$inc_pa

# ----- resolve a per-sample row for each side ---------------------------------
resolve_side <- function(pid, tp) {
  tp_chr <- as.character(tp)
  row <- resolve_sample(pid, tp_chr, core$assemblies)
  if (is.null(row)) {
    return(NULL)
  }
  row$SampleKey <- make_sample_key(pid, tp_chr)
  row
}

tools <- detect_tools()
timestamp_msg("Tools: ", paste(names(tools), unlist(tools), sep = "=", collapse = ", "))

# ----- compute per-pair metrics ----------------------------------------------
compute_pair <- function(pidA, tpA, pidB, tpB) {
  A <- resolve_side(pidA, tpA)
  B <- resolve_side(pidB, tpB)
  if (is.null(A) || is.null(B)) {
    return(NULL)
  }
  A <- as_tibble(A)
  B <- as_tibble(B)

  # MLST
  mlst_idx_A <- match(A$Isolate_ID, core$mlst$Isolate_ID)
  mlst_idx_B <- match(B$Isolate_ID, core$mlst$Isolate_ID)
  st_A <- core$mlst$ST[mlst_idx_A]
  st_B <- core$mlst$ST[mlst_idx_B]
  st_source_A <- if ("ST_source" %in% names(core$mlst)) core$mlst$ST_source[mlst_idx_A] else NA_character_
  st_source_B <- if ("ST_source" %in% names(core$mlst)) core$mlst$ST_source[mlst_idx_B] else NA_character_
  st_provider_A <- if ("ST_provider" %in% names(core$mlst)) core$mlst$ST_provider[mlst_idx_A] else NA_character_
  st_provider_B <- if ("ST_provider" %in% names(core$mlst)) core$mlst$ST_provider[mlst_idx_B] else NA_character_
  st_local_A <- if ("ST_local" %in% names(core$mlst)) core$mlst$ST_local[mlst_idx_A] else NA_character_
  st_local_B <- if ("ST_local" %in% names(core$mlst)) core$mlst$ST_local[mlst_idx_B] else NA_character_
  ST_equal <- !is.na(st_A) && !is.na(st_B) && st_A == st_B

  # dnadiff
  id_res <- empty_dnadiff_result()
  if (tools$dnadiff) {
    key <- paste0(A$SampleKey, "__vs__", B$SampleKey)
    fp_A <- list(
      path = A$full_path[[1]], sha256 = A$fasta_sha256[[1]],
      size_bytes = A$fasta_size_bytes[[1]], mtime_utc = A$fasta_mtime_utc[[1]]
    )
    fp_B <- list(
      path = B$full_path[[1]], sha256 = B$fasta_sha256[[1]],
      size_bytes = B$fasta_size_bytes[[1]], mtime_utc = B$fasta_mtime_utc[[1]]
    )
    id_res <- run_dnadiff(
      A$full_path[[1]], B$full_path[[1]], cache_dir, key,
      a_fingerprint = fp_A, b_fingerprint = fp_B
    )
  }

  # mash
  mash_d <- if (tools$mash) mash_distance(A$full_path, B$full_path) else NA_real_

  # VF Jaccard
  vf_df <- vf_wide %>% select(all_of(vf_meta_cols), all_of(vf_gene_cols))
  vf_j <- jaccard_from_wide(vf_df, row_key_cols = c("Participant_id", "tp_lab"), sampleA = A$SampleKey, sampleB = B$SampleKey)

  # Focus genes (optional)
  focus_j <- list(jaccard = NA_real_, n_int = NA_integer_, n_union = NA_integer_)
  focus_file <- c(
    file.path(DIR_VF, "diff_focus_genes_UTI_vs_Not_UTI_glmm.csv")
  )
  focus_file <- focus_file[file.exists(focus_file)][1]
  if (file.exists(focus_file)) {
    foc <- readr::read_csv(focus_file, show_col_types = FALSE) %>%
      select(FocusKey = any_of(c("FocusKey", "GENE"))) %>%
      distinct()
    if (nrow(foc)) {
      cols <- intersect(foc$FocusKey, vf_gene_cols)
      if (length(cols)) {
        vf_df_focus <- vf_wide %>% select(all_of(vf_meta_cols), all_of(cols))
        focus_j <- jaccard_from_wide(vf_df_focus, row_key_cols = c("Participant_id", "tp_lab"), sampleA = A$SampleKey, sampleB = B$SampleKey)
      }
    }
  }

  # Inc Jaccard (replicons)
  inc_j <- list(jaccard = NA_real_, n_int = NA_integer_, n_union = NA_integer_)
  if (!is.null(inc_pa)) {
    inc_join <- inc_pa
    # apply prevalence filtering if requested
    inc_cols <- setdiff(names(inc_join), "Isolate_ID")
    if (!is.na(opt$min_inc_prev) && length(inc_cols)) {
      prev <- colMeans((inc_join[, inc_cols, drop = FALSE] > 0), na.rm = TRUE)
      keep <- names(prev)[prev >= opt$min_inc_prev]
      inc_cols <- intersect(inc_cols, keep)
    }
    if (length(inc_cols)) {
      # map Isolate_ID to SampleKey
      map_df <- bind_rows(A %>% select(Isolate_ID, SampleKey), B %>% select(Isolate_ID, SampleKey)) %>% distinct()
      inc_samp <- inc_join %>%
        inner_join(map_df, by = "Isolate_ID") %>%
        select(SampleKey, all_of(inc_cols))
      # widen to include both rows for SampleKey
      inc_samp <- inc_samp %>%
        mutate(dummy = 1L) %>%
        pivot_wider(names_from = SampleKey, values_from = dummy, values_fill = 0L)
      # reconstruct a small df_wide with two rows
      inc_rows <- inc_join %>%
        filter(Isolate_ID %in% c(A$Isolate_ID, B$Isolate_ID)) %>%
        mutate(SampleKey = map_df$SampleKey[match(Isolate_ID, map_df$Isolate_ID)]) %>%
        select(SampleKey, all_of(inc_cols))
      if (nrow(inc_rows) >= 2) {
        inc_rows <- inc_rows %>% distinct(SampleKey, .keep_all = TRUE)
        # create fake pid/tp by splitting SampleKey
        split_keys <- strsplit(inc_rows$SampleKey, "__", fixed = TRUE)
        inc_rows$Participant_id <- vapply(split_keys, `[[`, character(1), 1)
        inc_rows$tp_lab <- vapply(split_keys, `[[`, character(1), 2)
        inc_j <- jaccard_from_wide(inc_rows, row_key_cols = c("Participant_id", "tp_lab"), sampleA = A$SampleKey, sampleB = B$SampleKey)
      }
    }
  }

  # pMLST equality counts
  pmlst_equal_n <- NA_integer_
  pmlst_equal_schemes <- NA_character_
  if (!is.null(core$pmlst_wide)) {
    pa <- core$pmlst_wide %>% filter(Isolate_ID == A$Isolate_ID)
    pb <- core$pmlst_wide %>% filter(Isolate_ID == B$Isolate_ID)
    if (nrow(pa) && nrow(pb)) {
      schemes <- intersect(names(pa), names(pb))
      schemes <- schemes[grepl("^pST_", schemes)]
      if (length(schemes)) {
        eq <- map_lgl(schemes, ~ {
          va <- as.character(pa[[.x]])
          vb <- as.character(pb[[.x]])
          !is.na(va) && !is.na(vb) && va == vb
        })
        pmlst_equal_n <- sum(eq)
        pmlst_equal_schemes <- paste(schemes[which(eq)], collapse = ",")
        if (pmlst_equal_schemes == "") pmlst_equal_schemes <- NA_character_
      }
    }
  }

  # build metrics row
  tibble(
    Participant_id_A = pidA,
    Timepoint_A = as.character(tpA),
    SampleKey_A = A$SampleKey,
    Isolate_ID_A = A$Isolate_ID,
    Assembler_A = as.character(A$assembler),
    Fasta_path_A = A$full_path,
    Fasta_SHA256_A = A$fasta_sha256,
    Fasta_size_bytes_A = A$fasta_size_bytes,
    Fasta_mtime_utc_A = A$fasta_mtime_utc,
    Selected_canonical_A = A$selected_canonical,
    QC_PASS_A = A$QC_PASS,
    ST_A = st_A,
    ST_source_A = st_source_A,
    ST_provider_A = st_provider_A,
    ST_local_A = st_local_A,
    Participant_id_B = pidB,
    Timepoint_B = as.character(tpB),
    SampleKey_B = B$SampleKey,
    Isolate_ID_B = B$Isolate_ID,
    Assembler_B = as.character(B$assembler),
    Fasta_path_B = B$full_path,
    Fasta_SHA256_B = B$fasta_sha256,
    Fasta_size_bytes_B = B$fasta_size_bytes,
    Fasta_mtime_utc_B = B$fasta_mtime_utc,
    Selected_canonical_B = B$selected_canonical,
    QC_PASS_B = B$QC_PASS,
    ST_B = st_B,
    ST_source_B = st_source_B,
    ST_provider_B = st_provider_B,
    ST_local_B = st_local_B,
    ST_equal = ST_equal,
    AvgIdentity = id_res$AvgIdentity,
    TotalSNPs = id_res$TotalSNPs,
    dnadiff_report_path = id_res$dnadiff_report_path,
    dnadiff_sidecar_path = id_res$dnadiff_sidecar_path,
    dnadiff_report_sha256 = id_res$dnadiff_report_sha256,
    dnadiff_cache_signature = id_res$dnadiff_cache_signature,
    dnadiff_cache_status = id_res$dnadiff_cache_status,
    dnadiff_version = id_res$dnadiff_version,
    MashDistance = mash_d,
    VF_Jaccard = vf_j$jaccard,
    VF_Overlap_A = vf_j$n_int, # store intersection count as A for shorthand
    VF_Overlap_B = vf_j$n_int, # symmetric
    VF_Union = vf_j$n_union,
    Inc_Jaccard = inc_j$jaccard,
    pST_equal_n = pmlst_equal_n,
    pST_equal_schemes = pmlst_equal_schemes
  )
}

timestamp_msg("Computing pairwise metrics for ", nrow(pairs), " pairs …")
workers <- max(1L, as.integer(opt$workers))
if (.Platform$OS.type == "windows") workers <- 1L
workers <- min(workers, nrow(pairs))
pair_indices <- seq_len(nrow(pairs))
compute_pair_index <- function(i) {
  compute_pair(
    pairs$Participant_id_A[[i]], pairs$Timepoint_A[[i]],
    pairs$Participant_id_B[[i]], pairs$Timepoint_B[[i]]
  )
}
pair_list <- if (workers > 1L) {
  timestamp_msg("Using ", workers, " forked workers; output row order remains input-pair order.")
  parallel::mclapply(pair_indices, compute_pair_index, mc.cores = workers, mc.preschedule = TRUE)
} else {
  lapply(pair_indices, compute_pair_index)
}
pair_metrics <- bind_rows(pair_list)

if (!nrow(pair_metrics)) stop("No metrics computed – check inputs and availability of assemblies")

# ----- classification ---------------------------------------------------------
apply_class <- function(...) {
  row <- list(...)
  res <- classify_pair(row, thresholds = list(id = opt$id_thresh, snps = opt$snp_thresh, vf = 0.9, inc = 0.8, vf_rel = 0.7, inc_rel = 0.7))
  c(Classification = res$Classification, RuleUsed = res$RuleUsed)
}

cls <- purrr::pmap_df(pair_metrics, apply_class)
pair_metrics <- bind_cols(pair_metrics, cls)

# within vs between label
pair_metrics <- pair_metrics %>%
  mutate(
    within_participant = Participant_id_A == Participant_id_B,
    ST_A_context = normalise_strain_st_label(ST_A),
    ST_B_context = normalise_strain_st_label(ST_B),
    snp_strain_context = case_when(
      is.na(TotalSNPs) ~ "Missing SNP evidence",
      TotalSNPs <= opt$snp_thresh ~ "Strong same strain",
      TRUE ~ "Above same-strain SNP threshold"
    ),
    st_lineage_context = case_when(
      is.na(ST_A_context) | is.na(ST_B_context) ~ "Missing ST evidence",
      ST_A_context == ST_B_context ~ "Same ST",
      TRUE ~ "Different ST"
    ),
    pair_interpretation = case_when(
      snp_strain_context == "Strong same strain" & st_lineage_context == "Different ST" ~
        "Conflict: SNP same-strain but ST differs",
      snp_strain_context == "Strong same strain" ~
        "Strong same strain",
      snp_strain_context != "Strong same strain" &
        (Classification == "Different" | st_lineage_context == "Different ST") ~
        "Replacement likely",
      snp_strain_context == "Above same-strain SNP threshold" & st_lineage_context == "Same ST" ~
        "Same lineage, not same strain by SNP",
      snp_strain_context == "Missing SNP evidence" & st_lineage_context == "Same ST" ~
        "ST-consistent, SNP missing",
      snp_strain_context == "Missing SNP evidence" & st_lineage_context == "Missing ST evidence" &
        (is.na(Classification) | Classification == "") ~
        "Missing strain metrics",
      snp_strain_context == "Above same-strain SNP threshold" ~
        "Above same-strain SNP threshold",
      snp_strain_context == "Missing SNP evidence" ~
        "Missing SNP evidence",
      TRUE ~ "Missing strain metrics"
    ),
    same_strain_evidence = pair_interpretation,
    strain_context_level = factor(
      pair_interpretation,
      levels = c("Strong same strain",
                 "Conflict: SNP same-strain but ST differs",
                 "Same lineage, not same strain by SNP",
                 "ST-consistent, SNP missing",
                 "Above same-strain SNP threshold",
                 "Missing SNP evidence",
                 "Replacement likely",
                 "Missing strain metrics")
    ),
    replacement_flag = pair_interpretation == "Replacement likely"
  ) %>%
  select(-ST_A_context, -ST_B_context)

# ----- provenance validation and write outputs -------------------------------
input_manifest <- core$assemblies %>%
  transmute(
    Participant_id, tp_lab, Isolate_ID, assembler,
    selected_canonical, QC_PASS,
    Fasta_path = full_path,
    Fasta_SHA256 = fasta_sha256,
    Fasta_size_bytes = fasta_size_bytes,
    Fasta_mtime_utc = fasta_mtime_utc
  ) %>%
  arrange(Participant_id, tp_lab)
safe_write_csv(input_manifest, file.path(outdir, "pairwise_input_manifest.csv"))

endpoint_ok <- with(
  pair_metrics,
  Selected_canonical_A %in% TRUE & QC_PASS_A %in% TRUE &
    Selected_canonical_B %in% TRUE & QC_PASS_B %in% TRUE &
    tolower(Assembler_A) == ANALYSIS_ASSEMBLER &
    tolower(Assembler_B) == ANALYSIS_ASSEMBLER
)
cache_ok <- if (tools$dnadiff) {
  vapply(seq_len(nrow(pair_metrics)), function(i) {
    spec <- dnadiff_cache_spec(
      pair_metrics$Fasta_path_A[[i]], pair_metrics$Fasta_path_B[[i]], cache_dir,
      paste0(pair_metrics$SampleKey_A[[i]], "__vs__", pair_metrics$SampleKey_B[[i]]),
      a_fingerprint = list(
        path = pair_metrics$Fasta_path_A[[i]], sha256 = pair_metrics$Fasta_SHA256_A[[i]],
        size_bytes = pair_metrics$Fasta_size_bytes_A[[i]], mtime_utc = pair_metrics$Fasta_mtime_utc_A[[i]]
      ),
      b_fingerprint = list(
        path = pair_metrics$Fasta_path_B[[i]], sha256 = pair_metrics$Fasta_SHA256_B[[i]],
        size_bytes = pair_metrics$Fasta_size_bytes_B[[i]], mtime_utc = pair_metrics$Fasta_mtime_utc_B[[i]]
      )
    )
    dnadiff_cache_is_valid(spec) &&
      identical(spec$cache_signature, pair_metrics$dnadiff_cache_signature[[i]])
  }, logical(1))
} else {
  rep(NA, nrow(pair_metrics))
}
provenance_audit <- tibble(
  metric = c(
    "selected_qc_passing_input_assemblies", "pairwise_rows",
    "noncanonical_endpoint_rows", "stale_or_missing_dnadiff_cache_rows"
  ),
  value = c(
    nrow(input_manifest), nrow(pair_metrics), sum(!endpoint_ok),
    if (tools$dnadiff) sum(!cache_ok) else NA_integer_
  ),
  requirement = c(
    "All manifest rows selected_canonical=TRUE, QC_PASS=TRUE, and assembler=longcycler",
    "One row per requested unique pair",
    "Must equal 0", "Must equal 0 when dnadiff is available"
  )
)
safe_write_csv(provenance_audit, file.path(outdir, "pairwise_provenance_audit.csv"))
if (any(!endpoint_ok)) stop("Provenance audit failed: noncanonical or QC-failing endpoint detected.")
if (tools$dnadiff && any(!cache_ok)) {
  stop("Provenance audit failed: stale, missing, or mismatched dnadiff cache detected.")
}

safe_write_csv(pair_metrics, file.path(outdir, "pairwise_metrics.csv"))

sum_counts <- pair_metrics %>%
  count(Classification, name = "n") %>%
  arrange(desc(n))
safe_write_csv(sum_counts, file.path(outdir, "summary_counts.csv"))

by_participant <- pair_metrics %>%
  mutate(Participant = ifelse(within_participant, Participant_id_A, "(between)")) %>%
  count(Participant, Classification, name = "n")
safe_write_csv(by_participant, file.path(outdir, "summary_by_participant.csv"))

# ----- plots ------------------------------------------------------------------
plots_dir <- file.path(outdir, "plots")
dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)
pair_heatmap_height <- function(n) min(24, max(4, 0.02 * n))

# Heatmaps
if (any(!is.na(pair_metrics$VF_Jaccard))) {
  g <- ggplot(pair_metrics %>% mutate(idx = row_number()), aes(x = 1, y = idx, fill = VF_Jaccard)) +
    geom_tile() +
    scale_fill_viridis_c(option = "C", na.value = "grey85") +
    labs(title = "Pairwise Virulence Factor Jaccard Similarity", x = NULL, y = "Pair index") +
    theme_minimal(base_size = 11) +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
  ggsave_safe(file.path(plots_dir, "heatmap_vf_jaccard.png"), g, width = 4, height = pair_heatmap_height(nrow(pair_metrics)))
}

if (any(!is.na(pair_metrics$Inc_Jaccard))) {
  g <- ggplot(pair_metrics %>% mutate(idx = row_number()), aes(x = 1, y = idx, fill = Inc_Jaccard)) +
    geom_tile() +
    scale_fill_viridis_c(option = "C", na.value = "grey85") +
    labs(title = "Pairwise Plasmid Replicon Jaccard Similarity", x = NULL, y = "Pair index") +
    theme_minimal(base_size = 11) +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
  ggsave_safe(file.path(plots_dir, "heatmap_inc_jaccard.png"), g, width = 4, height = pair_heatmap_height(nrow(pair_metrics)))
}

# Identity vs SNPs
if (any(!is.na(pair_metrics$AvgIdentity) & !is.na(pair_metrics$TotalSNPs))) {
  g <- ggplot(pair_metrics, aes(AvgIdentity, TotalSNPs, color = Classification, shape = within_participant)) +
    geom_point(size = 2, alpha = 0.8) +
    scale_y_continuous(transform = "log1p") +
    labs(title = "Genomic Identity vs. SNP Distance", x = "AvgIdentity (%)", y = "Total SNPs (log1p)") +
    theme_minimal(base_size = 11)
  ggsave_safe(file.path(plots_dir, "identity_vs_snps_scatter.png"), g, width = 6, height = 4)

  # Violin plot for SNP distances (Within vs Between)
  g_violin <- ggplot(pair_metrics, aes(x = ifelse(within_participant, "Within-Host", "Between-Host"), y = TotalSNPs + 1, fill = ifelse(within_participant, "Within-Host", "Between-Host"))) +
    geom_violin(alpha = 0.6) +
    geom_jitter(width = 0.2, alpha = 0.2, size = 0.5) +
    scale_y_log10() +
    labs(title = "Pairwise SNP Distances: Within-Host vs. Between-Host", x = NULL, y = "SNP Distance (log10 + 1)", fill = "Comparison") +
    theme_minimal(base_size = 11) +
    theme(legend.position = "none")
  ggsave_safe(file.path(plots_dir, "snp_distance_violin.png"), g_violin, width = 6, height = 5)
}

# Same-strain network (requires igraph + ggraph)
if (requireNamespace("igraph", quietly = TRUE) && requireNamespace("ggraph", quietly = TRUE)) {
  igraph <- asNamespace("igraph")
  ggraph <- asNamespace("ggraph")
  edges <- pair_metrics %>%
    filter(Classification == "Same") %>%
    transmute(from = SampleKey_A, to = SampleKey_B)
  if (nrow(edges)) {
    nodes <- tibble(name = unique(c(edges$from, edges$to)))
    g <- igraph$graph_from_data_frame(edges, directed = FALSE, vertices = nodes)
    p <- ggraph$ggraph(g, layout = "fr") +
      ggraph$geom_edge_link(alpha = .4) +
      ggraph$geom_node_point(size = 4, colour = "steelblue") +
      ggraph$geom_node_text(ggplot2::aes(label = name), repel = TRUE, size = 3) +
      ggplot2::theme_void() + ggplot2::ggtitle("Network of Same-Strain Isolates")
    ggsave_safe(file.path(plots_dir, "network_same_strain.png"), p, width = 7, height = 6)
  }
}

# Within-participant timeline: use successive pairs if present among pairs
within_pairs <- pair_metrics %>%
  filter(within_participant) %>%
  mutate(pid = Participant_id_A)
if (nrow(within_pairs)) {
  # attempt to order timepoints numerically when possible
  ord <- function(tp) suppressWarnings(as.integer(stringr::str_extract(tp, "\\d+")))
  within_pairs <- within_pairs %>% mutate(tpA_num = ord(Timepoint_A), tpB_num = ord(Timepoint_B))
  g <- ggplot(within_pairs, aes(x = tpA_num, xend = tpB_num, y = pid, yend = pid, color = Classification)) +
    geom_segment(linewidth = 1.2, alpha = 0.8) +
    scale_x_continuous(breaks = sort(unique(c(within_pairs$tpA_num, within_pairs$tpB_num))), limits = c(min(within_pairs$tpA_num, na.rm = TRUE) - 0.5, max(within_pairs$tpB_num, na.rm = TRUE) + 0.5)) +
    labs(title = "Longitudinal Strain Classification Timeline", x = "Numeric timepoint", y = "Participant") +
    theme_minimal(base_size = 11)
  ggsave_safe(file.path(plots_dir, "timeline_by_participant.png"), g, width = 8, height = max(4, 0.3 * dplyr::n_distinct(within_pairs$pid)))
}

# ----- statistics -------------------------------------------------------------
# helper: kappa between ST_equal and Same vs not same
compute_kappa <- function(st_equal, classification) {
  same <- classification == "Same"
  valid <- !is.na(st_equal) & !is.na(same)
  if (!any(valid)) {
    return(NA_real_)
  }
  st <- st_equal[valid]
  sm <- same[valid]
  tab <- table(st, sm)
  if (nrow(tab) != 2 || ncol(tab) != 2) {
    return(NA_real_)
  }
  n <- sum(tab)
  po <- (tab[1, 1] + tab[2, 2]) / n
  pe <- ((sum(tab[1, ]) * sum(tab[, 1])) + (sum(tab[2, ]) * sum(tab[, 2]))) / (n^2)
  if (is.nan(po) || is.nan(pe) || (1 - pe) == 0) {
    return(NA_real_)
  }
  (po - pe) / (1 - pe)
}

pair_metrics <- pair_metrics %>% mutate(Same_vs_not = ifelse(Classification == "Same", "Same", "NotSame"))

# within vs between Wilcoxon tests
wilcox_metric <- function(metric) {
  df <- pair_metrics %>%
    select(within_participant, {{ metric }}) %>%
    rename(val = {{ metric }}) %>%
    filter(!is.na(val))
  if (nrow(df) < 3) {
    return(tibble(metric = deparse(substitute(metric)), p_value = NA_real_))
  }
  p <- tryCatch(stats::wilcox.test(val ~ within_participant, data = df)$p.value, error = function(e) NA_real_)
  tibble(metric = deparse(substitute(metric)), p_value = p)
}

stats_wb <- bind_rows(
  wilcox_metric(AvgIdentity),
  wilcox_metric(TotalSNPs),
  wilcox_metric(MashDistance),
  wilcox_metric(VF_Jaccard),
  wilcox_metric(Inc_Jaccard)
) %>% mutate(p_adj = p.adjust(p_value, method = "BH"))
safe_write_csv(stats_wb, file.path(outdir, "stats_within_vs_between.csv"))

# status-based comparisons (if status_map exists)
if (!is.null(core$status_map)) {
  # annotate each sample with status and propagate back to pairs
  anno <- core$status_map
  A_anno <- anno %>% rename(Timepoint_A = tp_lab, Status_A = Infection_Status)
  B_anno <- anno %>% rename(Timepoint_B = tp_lab, Status_B = Infection_Status)
  pm <- pair_metrics %>%
    left_join(A_anno, by = c("Participant_id_A" = "Participant_id", "Timepoint_A" = "Timepoint_A")) %>%
    left_join(B_anno, by = c("Participant_id_B" = "Participant_id", "Timepoint_B" = "Timepoint_B"))

  # For simplicity: compare VF_Jaccard by concatenated status pair label
  pm <- pm %>% mutate(StatusPair = paste(Status_A, Status_B, sep = "→"))
  stat <- pm %>% filter(!is.na(VF_Jaccard), !is.na(StatusPair))
  if (nrow(stat) >= 3 && dplyr::n_distinct(stat$StatusPair) >= 2) {
    p_kw <- tryCatch(stats::kruskal.test(VF_Jaccard ~ StatusPair, data = stat)$p.value, error = function(e) NA_real_)
    safe_write_csv(tibble(metric = "VF_Jaccard", test = "Kruskal-Wallis", p_value = p_kw), file.path(outdir, "stats_by_status.csv"))
  }

  # [STAT] NOTE: The Kruskal-Wallis test above is exploratory.
  # To rigorously test "Does strain persistence predict UTI recurrence?" or
  # "Does gene content change predict status change?", use a GLMM:
  # glmer(Status_B ~ Same_Strain + VF_Jaccard + (1|Participant_id), family=binomial)
  # This requires defining the outcome (e.g., Recurrent UTI) more precisely.
}

# ----- README -----------------------------------------------------------------
cmd_line <- paste(commandArgs(trailingOnly = FALSE), collapse = " ")
readme <- c(
  "Strain comparison results",
  paste0("Generated: ", Sys.time()),
  paste0("Command: ", cmd_line),
  paste0("Identity threshold: ", opt$id_thresh, "; SNP threshold: ", opt$snp_thresh),
  "Assemblies: selected_canonical=TRUE, QC_PASS=TRUE, assembler=longcycler only; no assembler fallback",
  paste0("dnadiff cache: SHA-256 input-bound cache in ", cache_dir),
  paste0("Workers: ", workers),
  paste0("min_vf_prev: ", opt$min_vf_prev, "; min_inc_prev: ", opt$min_inc_prev),
  "See pairwise_metrics.csv, summary_counts.csv, and plots/ for outputs."
)
writeLines(readme, file.path(outdir, "README.txt"))

timestamp_msg("✓ Completed. Outputs in ", outdir)

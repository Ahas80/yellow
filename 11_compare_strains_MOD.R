#!/usr/bin/env Rscript
# =============================================================
# 11_within_person_strain_compare.R
# -------------------------------------------------------------
# Purpose:
#   Within-participant strain comparison of E. coli across time,
#   extending the Yellow RoUTIne pipeline. Compares consecutive
#   timepoints per participant and quantifies similarity via:
#     • MLST equality (ST, majority vote per PID×TP)
#     • VF gene presence/absence Jaccard distance
#     • Plasmid replicon Jaccard distance (if mappable)
#     • ANI / SNPs from nucmer/dnadiff (if available)
#
# Inputs (expected from prior scripts):
#   - results/vf_pa_all.csv                      (from 02_gene_presence_analysis.R)
#   - results/mlst/mlst_matrix.csv               (from 06_MLST.R)
#   - results/status_map.csv                     (from 01_prepare_assembly_metadata.R; optional)
#   - results/plasmidfinder_presence_absence.csv (from 09_inc_plasmid_network.R; optional)
#   - results/nucmer/*/run_dd.report             (if nucmer trajectories ran in 02)
#
# Outputs:
#   - results/within_person/pairwise_summary.csv
#   - results/within_person/persistence_by_participant.csv
#   - results/within_person/gene_changes_long.csv
#   - results/within_person/benchmark_within_vs_between.csv
#   - results/within_person/nucmer_pairwise_metrics.csv         (if parsed)
#   - results/within_person/debug/*.csv                         (conflicts, missing, parsing)
#   - results/within_person/params.json
#   - results/plots/within_person/*.png                         (distributions & trajectories)
#
# Key assumptions (edit thresholds via CLI if needed):
#   1) Timepoints are normalized to labels like "T0","T1",…,"Uricult" (as in prior scripts).
#   2) "Persistent (same strain)" default rule:
#        - Same ST AND
#        - VF Jaccard distance <= 0.05 AND
#        - If ANI/SNPs available: ANI >= 99.90% (or SNPs <= 20)
#      Otherwise "Replacement". Missing/discordant metrics → "Uncertain".
#   3) If multiple isolates exist within the same Participant×Timepoint,
#      VF presence/absence collapses by any-present; ST uses majority vote
#      (conflicts logged to /results/within_person/debug/).
#
# CLI (optional):
#   --pids=P01,P02        (restrict to these participants)
#   --tp_rule=all         ("all" = consecutive by factor order; reserved for future)
#   --vf_jaccard_max=0.05 (threshold for persistence)
#   --ani_min=99.9        (threshold for persistence if ANI present)
#   --snp_max=20          (threshold for persistence if SNP counts present)
#
# Style & behavior:
#   • Only needed libraries; startup messages suppressed.
#   • Creates results/ and subfolders; logs major steps with message().
#   • Fails gracefully with informative errors if required inputs missing.
# =============================================================

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readr); library(stringr)
  library(purrr); library(ggplot2); library(forcats); library(jsonlite)
})

# ------------------------------ 0 · helpers & dirs ----------------------------
# Normalise isolate IDs so joins always match (drops paths & .fasta/.fa/.fna[.gz])
norm_iso_id <- function(x) {
  x <- as.character(x)
  x <- basename(x)
  sub("\\.(fa|fna|fasta)(\\.gz)?$", "", x, ignore.case = TRUE)
}

base_dir <- normalizePath(getwd(), mustWork = TRUE)
out_dir  <- file.path(base_dir, "results", "within_person")
plot_dir <- file.path(base_dir, "results", "plots", "within_person")
dbg_dir  <- file.path(base_dir, "results", "within_person", "debug")
dir.create(out_dir,  recursive = TRUE, showWarnings = FALSE)
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(dbg_dir,  recursive = TRUE, showWarnings = FALSE)

msg <- function(fmt, ...) {
  cat(format(Sys.time(), "[%Y-%m-%d %H:%M:%S] "), sprintf(fmt, ...), "\n")
  flush.console()
}

safe_write_csv <- function(x, path, ...) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  tmp <- paste0(path, ".tmp")
  ok <- tryCatch({ readr::write_csv(x, tmp, ...); TRUE }, error = function(e) { message(conditionMessage(e)); FALSE })
  if (ok && file.exists(tmp)) file.rename(tmp, path)
  if (file.exists(path)) msg("wrote: %s  (nrow=%d)", normalizePath(path), nrow(x))
}
safe_save_plot <- function(filename, plot, width, height, dpi = 300, units = "in") {
  try({
    dir.create(dirname(filename), recursive = TRUE, showWarnings = FALSE)
    if (requireNamespace("ragg", quietly = TRUE)) {
      ggplot2::ggsave(filename, plot, device = ragg::agg_png,
                      width = width, height = height, dpi = dpi, units = units)
    } else {
      ggplot2::ggsave(filename, plot, width = width, height = height, dpi = dpi, units = units)
    }
    msg("saved plot: %s", normalizePath(filename))
  }, silent = TRUE)
}



# Timepoint normalization used throughout the pipeline
tp_norm <- function(x) {
  tp_chr     <- as.character(x)
  is_uricult <- str_detect(tp_chr, regex("uricult", ignore_case = TRUE))
  tp_num     <- suppressWarnings(as.integer(str_extract(tp_chr, "\\d+")))
  tp_lab     <- dplyr::case_when(
    is_uricult     ~ "Uricult",
    !is.na(tp_num) ~ paste0("T", tp_num),
    TRUE           ~ "Unscheduled"
  )
  tp_levels <- c(paste0("T", sort(unique(tp_num[!is.na(tp_num)]))),
                 "Uricult", "Unscheduled")
  tibble::tibble(tp_lab = factor(tp_lab, levels = tp_levels),
                 tp_num = tp_num)
}

# Ensure we have tp_lab (+tp_num) without duplicating columns
ensure_tp <- function(df) {
  out <- df
  if ("tp_lab" %in% names(out)) {
    # Only compute tp_num if missing, using the current tp_lab
    if (!"tp_num" %in% names(out)) {
      tp_chr <- as.character(out$tp_lab)
      # Derive numeric part; Uricult/Unscheduled -> NA
      tp_num <- suppressWarnings(as.integer(stringr::str_extract(tp_chr, "\\d+")))
      out$tp_num <- tp_num
    }
  } else if ("Timepoint" %in% names(out)) {
    out <- bind_cols(out, tp_norm(out$Timepoint))
  } else {
    stop("Data frame lacks both 'tp_lab' and 'Timepoint' fields.")
  }
  out
}

# Jaccard distance for binary named vectors (0/1); returns NA if union = 0
jaccard_distance01 <- function(a, b) {
  a <- as.integer(a > 0); b <- as.integer(b > 0)
  inter <- sum(a & b, na.rm = TRUE)
  uni   <- sum((a | b), na.rm = TRUE)
  if (uni == 0) return(NA_real_)
  1 - (inter / uni)
}

# ------------------------------ 1 · parse CLI & thresholds --------------------
args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(key, default = NULL) {
  hit <- args[grepl(paste0("^--", key, "="), args)]
  if (!length(hit)) return(default)
  sub(paste0("^--", key, "="), "", hit[1])
}

PID_FILTER <- get_arg("pids", default = NA_character_)
TP_RULE    <- get_arg("tp_rule", default = "all")
VF_JAC_MAX <- suppressWarnings(as.numeric(get_arg("vf_jaccard_max", "0.05")))
ANI_MIN    <- suppressWarnings(as.numeric(get_arg("ani_min", "99.9")))
SNP_MAX    <- suppressWarnings(as.integer(get_arg("snp_max", "20")))
if (is.na(VF_JAC_MAX)) VF_JAC_MAX <- 0.05
if (is.na(ANI_MIN))    ANI_MIN    <- 99.9
if (is.na(SNP_MAX))    SNP_MAX    <- 20L

params <- list(
  pid_filter = if (is.na(PID_FILTER)) "ALL" else PID_FILTER,
  tp_rule = TP_RULE, vf_jaccard_max = VF_JAC_MAX,
  ani_min = ANI_MIN, snp_max = SNP_MAX
)
writeLines(jsonlite::toJSON(params, pretty = TRUE, auto_unbox = TRUE),
           file.path(out_dir, "params.json"))

# ------------------------------ 2 · check required inputs ---------------------
vf_path   <- "results/vf_pa_all.csv"
mlst_path <- "results/mlst/mlst_matrix.csv"
if (!file.exists(vf_path) || !file.exists(mlst_path)) {
  if (!file.exists(vf_path))   msg("ERROR: %s not found (run 02_gene_presence_analysis.R)", vf_path)
  if (!file.exists(mlst_path)) msg("ERROR: %s not found (run 06_MLST.R)", mlst_path)
  stop("[11] Cannot proceed without required inputs.")
}
status_map_path <- "results/status_map.csv"
replicon_path   <- "results/plasmidfinder_presence_absence.csv"

# ------------------------------ 3 · load input tables -------------------------
msg("Loading VF presence/absence …")
vf <- suppressMessages(readr::read_csv(vf_path, show_col_types = FALSE))
vf <- ensure_tp(vf)
vf <- vf %>% mutate(Participant_id = as.character(Participant_id)) %>%
  arrange(Participant_id, tp_lab)

# Identify gene columns; standardize binary & collapse duplicates by any-present
gene_cols <- setdiff(names(vf), c("Participant_id", "tp_lab", "tp_num", "Timepoint"))
if (!length(gene_cols)) stop("No VF gene columns detected in vf_pa_all.csv.")
vf[gene_cols] <- lapply(vf[gene_cols], function(x) {
  x <- suppressWarnings(as.integer(x)); x[is.na(x)] <- 0L; as.integer(x > 0L)
})
vf <- vf %>%
  group_by(Participant_id, tp_lab, tp_num) %>%
  summarise(across(all_of(gene_cols), ~ as.integer(any(. > 0, na.rm = TRUE))), .groups = "drop")

msg("Loading MLST table …")
mlst <- suppressMessages(readr::read_csv(mlst_path, show_col_types = FALSE))

# If mlst_matrix.csv lacks Isolate_ID, fall back to mlst_all.tsv
if (!"Isolate_ID" %in% names(mlst) && file.exists("results/mlst/mlst_all.tsv")) {
  mlst <- suppressMessages(readr::read_tsv("results/mlst/mlst_all.tsv", show_col_types = FALSE))
}

# Normalise isolate IDs early (whether they came from csv or tsv)
if ("Isolate_ID" %in% names(mlst)) {
  mlst$Isolate_ID <- norm_iso_id(mlst$Isolate_ID)
}

# If Participant_id/Timepoint missing, join from a metadata file
if (!all(c("Participant_id","Timepoint") %in% names(mlst))) {
  meta_candidates <- c("results/assembly_metadata.csv",
                       "assembly_metadata.csv",
                       "results/metadata.csv")
  meta_path <- meta_candidates[file.exists(meta_candidates)][1]
  
  if (is.na(meta_path)) {
    # Log what we have for forensics and fail clearly
    dir.create(dbg_dir, showWarnings = FALSE, recursive = TRUE)
    readr::write_lines(paste("mlst columns:", paste(names(mlst), collapse=", ")),
                       file.path(dbg_dir, "mlst_columns_seen.txt"))
    stop(
      "MLST table is missing Participant_id/Timepoint and no mapping file found.\n",
      "Provide a CSV with columns: Isolate_ID, Participant_id, Timepoint.\n",
      "Checked: ", paste(meta_candidates, collapse = ", "), "\n",
      "See debug/mlst_columns_seen.txt"
    )
  }
  
  meta <- suppressMessages(readr::read_csv(meta_path, show_col_types = FALSE))
  
  # Normalise common variants
  if ("SampleID"       %in% names(meta)) names(meta)[names(meta)=="SampleID"]       <- "Isolate_ID"
  if ("PID"            %in% names(meta)) names(meta)[names(meta)=="PID"]            <- "Participant_id"
  if ("participant_id" %in% names(meta)) names(meta)[names(meta)=="participant_id"] <- "Participant_id"
  if ("tp"             %in% names(meta)) names(meta)[names(meta)=="tp"]             <- "Timepoint"
  if ("timepoint"      %in% names(meta)) names(meta)[names(meta)=="timepoint"]      <- "Timepoint"
  
  need <- c("Isolate_ID","Participant_id","Timepoint")
  if (!all(need %in% names(meta))) {
    dir.create(dbg_dir, showWarnings = FALSE, recursive = TRUE)
    readr::write_lines(paste("metadata columns:", paste(names(meta), collapse=", ")),
                       file.path(dbg_dir, "metadata_columns_seen.txt"))
    stop("Mapping file ", meta_path, " must contain: ", paste(need, collapse = ", "),
         "\nSee debug/metadata_columns_seen.txt")
  }
  
  # Normalise IDs on both sides before join
  meta$Isolate_ID <- norm_iso_id(meta$Isolate_ID)
  
  mlst <- mlst %>%
    dplyr::left_join(meta %>%
                       dplyr::select(Isolate_ID, Participant_id, Timepoint) %>%
                       dplyr::distinct(),
                     by = "Isolate_ID",
                     relationship = "many-to-one")
}

# Now we can produce tp_lab/tp_num safely
mlst <- ensure_tp(mlst)

# Ensure ST as character
st_col <- names(mlst)[tolower(names(mlst)) == "st"]
if (!length(st_col)) stop("MLST table lacks an 'ST' column.")
mlst$ST <- as.character(mlst[[st_col[1]]])

mlst_key <- mlst %>%
  mutate(Participant_id = as.character(Participant_id)) %>%
  select(Participant_id, tp_lab, ST) %>%
  filter(!is.na(Participant_id), !is.na(tp_lab)) %>%
  group_by(Participant_id, tp_lab) %>%
  summarise(
    ST_majority = {
      s <- na.omit(ST)
      if (!length(s)) NA_character_ else names(sort(table(s), decreasing = TRUE))[1]
    },
    ST_n_unique = n_distinct(ST, na.rm = TRUE),
    .groups = "drop"
  )

conflicts <- mlst_key %>% filter(ST_n_unique > 1)
if (nrow(conflicts)) {
  msg("⚠ Detected %d Participant×Timepoint with multiple ST calls; majority vote will be used, conflicts logged.", nrow(conflicts))
  safe_write_csv(conflicts, file.path(dbg_dir, "mlst_conflicts.csv"))
}

# Optional status map (used only for annotations/plots)
status_map <- NULL
if (file.exists(status_map_path)) {
  status_map <- suppressMessages(readr::read_csv(status_map_path, show_col_types = FALSE))
  status_map <- ensure_tp(status_map)
  status_map <- status_map %>% select(Participant_id, tp_lab, Infection_Status) %>% distinct()
}

# Optional replicon matrix (map to PID×TP via Isolate_ID if possible)
# Optional replicon matrix (map to PID×TP via Isolate_ID if possible)
rep_pid_tp <- NULL
if (file.exists(replicon_path)) {
  repl <- suppressMessages(readr::read_csv(replicon_path, show_col_types = FALSE))
  
  if ("Isolate_ID" %in% names(repl) && "Isolate_ID" %in% names(mlst)) {
    # Normalise both sides (handles .fasta vs no extension, any paths)
    repl$Isolate_ID <- norm_iso_id(repl$Isolate_ID)
    mlst$Isolate_ID <- norm_iso_id(mlst$Isolate_ID)
    
    key <- mlst %>%
      select(Isolate_ID, Participant_id, tp_lab) %>%
      distinct() %>%
      filter(!is.na(Isolate_ID), !is.na(Participant_id), !is.na(tp_lab))
    
    rep_join <- repl %>%
      mutate(Isolate_ID = as.character(Isolate_ID)) %>%
      left_join(key, by = "Isolate_ID", relationship = "many-to-one") %>%
      filter(!is.na(Participant_id), !is.na(tp_lab))
    
    rep_cols <- setdiff(names(rep_join), c("Isolate_ID","Participant_id","tp_lab"))
    if (length(rep_cols)) {
      rep_join[rep_cols] <- lapply(rep_join[rep_cols], function(x) {
        x <- suppressWarnings(as.integer(x)); x[is.na(x)] <- 0L; as.integer(x > 0L)
      })
      rep_pid_tp <- rep_join %>%
        group_by(Participant_id, tp_lab) %>%
        summarise(across(all_of(rep_cols), ~ as.integer(any(. > 0, na.rm = TRUE))), .groups = "drop")
      msg("✓ Replicon matrix mapped to Participant×Timepoint (n=%d rows).", nrow(rep_pid_tp))
    } else {
      msg("↪ Replicon file found but no replicon columns to summarise.")
    }
  } else {
    msg("↪ Replicon mapping skipped: missing 'Isolate_ID' in one or both tables.")
  }
}

# ------------------------------ 4 · optional nucmer/dnadiff parsing -----------
parse_dnadiff <- function(report_file) {
  L <- readLines(report_file, warn = FALSE)
  grab <- function(key, pattern = "\\d+\\.?\\d*") {
    m <- grep(key, L, value = TRUE)
    if (!length(m)) return(NA_real_)
    as.numeric(stringr::str_extract(m[1], pattern))
  }
  tibble::tibble(
    AvgIdentity = grab("AvgIdentity"),
    TotalSNPs   = grab("TotalSNPs"),
    RefLen      = grab("^TotalBasesRef"),
    QryLen      = grab("^TotalBasesQry")
  )
}
collect_nucmer <- function() {
  root <- file.path("results","nucmer")
  if (!dir.exists(root)) return(tibble())
  dirs <- list.dirs(root, recursive = TRUE, full.names = TRUE)
  rpt  <- file.path(dirs, "run_dd.report")
  rpt  <- rpt[file.exists(rpt)]
  if (!length(rpt)) return(tibble())
  purrr::map_dfr(rpt, function(f) {
    # Expected dir name: results/nucmer/PID_ASM_Tx_vs_Ty
    dn <- basename(dirname(f))
    parts <- strsplit(dn, "_")[[1]]
    if (length(parts) < 4) return(tibble())
    pid <- parts[1]
    assembler <- parts[2]
    ab   <- strsplit(paste(parts[3:length(parts)], collapse = "_"), "_vs_")[[1]]
    if (length(ab) != 2) return(tibble())
    tpA  <- ab[1]; tpB <- ab[2]
    parse_dnadiff(f) %>%
      mutate(Participant_id = pid, assembler = assembler, tp_A = tpA, tp_B = tpB,
             report_path = f)
  })
}
nucmer_tbl <- collect_nucmer()
if (nrow(nucmer_tbl)) {
  safe_write_csv(nucmer_tbl, file.path(out_dir, "nucmer_pairwise_metrics.csv"))
  msg("✓ Parsed %d nucmer/dnadiff reports.", nrow(nucmer_tbl))
} else {
  msg("↪ No nucmer/dnadiff reports found; ANI/SNP metrics will be unavailable.")
}

# ------------------------------ 5 · restrict PIDs if requested ----------------
if (!is.na(PID_FILTER)) {
  keep <- unlist(strsplit(PID_FILTER, "\\s*,\\s*"))
  vf        <- vf        %>% filter(Participant_id %in% keep)
  mlst_key  <- mlst_key  %>% filter(Participant_id %in% keep)
  if (!is.null(status_map)) status_map <- status_map %>% filter(Participant_id %in% keep)
  if (!is.null(rep_pid_tp)) rep_pid_tp <- rep_pid_tp %>% filter(Participant_id %in% keep)
  if (nrow(vf) == 0) stop("After --pids filter, no rows remain.")
}

# ------------------------------ 6 · build sample key table --------------------
# NOTE: avoid duplicating tp_lab — we already have it; only ensure tp_num.
tp_info <- vf %>%
  select(Participant_id, tp_lab, tp_num) %>%
  mutate(tp_num = suppressWarnings(as.integer(tp_num))) %>%
  distinct() %>%
  arrange(Participant_id, as.integer(tp_num), tp_lab)

samples <- tp_info %>%
  left_join(mlst_key %>% select(Participant_id, tp_lab, ST = ST_majority),
            by = c("Participant_id","tp_lab"))

missing_st <- samples %>% filter(is.na(ST))
if (nrow(missing_st)) safe_write_csv(missing_st, file.path(dbg_dir, "missing_ST_pid_tp.csv"))

# ------------------------------ 7 · build consecutive pairs -------------------
pairs <- samples %>%
  group_by(Participant_id) %>%
  arrange(Participant_id,
          # Order: numeric T's by tp_num; then Uricult; then Unscheduled
          is.na(tp_num), tp_num, factor(tp_lab, levels = levels(tp_lab))) %>%
  mutate(tp_next = lead(tp_lab), tp_num_next = lead(tp_num), ST_next = lead(ST)) %>%
  ungroup() %>%
  filter(!is.na(tp_next)) %>%
  transmute(Participant_id,
            tp_A = tp_lab,    tp_B = tp_next,
            tpA_num = tp_num, tpB_num = tp_num_next,
            ST_A = ST, ST_B = ST_next)

if (!nrow(pairs)) {
  msg("No consecutive pairs available. Provide participants with ≥2 timepoints.")
  safe_write_csv(samples, file.path(out_dir, "sample_key.csv"))
  quit(save = "no", status = 0)
}

# ------------------------------ 8 · compute VF Jaccard & gene changes ---------
vf_mat <- vf %>% select(Participant_id, tp_lab, all_of(gene_cols))

get_vf_vec <- function(pid, tp) {
  row <- vf_mat %>% filter(Participant_id == pid, tp_lab == tp)
  if (nrow(row) != 1L) return(setNames(rep(NA_integer_, length(gene_cols)), gene_cols))
  v <- as.integer(row[1, gene_cols, drop = TRUE])
  names(v) <- gene_cols
  v
}

pairs <- pairs %>%
  mutate(
    vfA_vec = purrr::map2(Participant_id, tp_A, \(pid, tp) get_vf_vec(pid, tp)),
    vfB_vec = purrr::map2(Participant_id, tp_B, \(pid, tp) get_vf_vec(pid, tp)),
    vf_jaccard = purrr::map2_dbl(vfA_vec, vfB_vec, jaccard_distance01),
    vf_intersect = purrr::map2(vfA_vec, vfB_vec, \(a,b) names(which(a == 1L & b == 1L))),
    vf_gain     = purrr::map2(vfA_vec, vfB_vec, \(a,b) names(which(a == 0L & b == 1L))),
    vf_loss     = purrr::map2(vfA_vec, vfB_vec, \(a,b) names(which(a == 1L & b == 0L)))
  )

gene_changes_long <- pairs %>%
  select(Participant_id, tp_A, tp_B, vf_gain, vf_loss) %>%
  tidyr::pivot_longer(c(vf_gain, vf_loss), names_to = "change_type", values_to = "GENE") %>%
  tidyr::unnest(GENE, keep_empty = TRUE) %>%
  mutate(change_type = dplyr::recode(change_type, vf_gain = "gained_in_B", vf_loss = "lost_in_B"))

safe_write_csv(gene_changes_long, file.path(out_dir, "gene_changes_long.csv"))

# ------------------------------ 9 · plasmid replicon Jaccard (optional) -------
rep_cols <- NULL
if (!is.null(rep_pid_tp)) {
  rep_cols <- setdiff(names(rep_pid_tp), c("Participant_id","tp_lab"))
  get_rep_vec <- function(pid, tp) {
    row <- rep_pid_tp %>% filter(Participant_id == pid, tp_lab == tp)
    if (nrow(row) != 1L) return(rep(NA_integer_, length(rep_cols)))
    as.integer(row[1, rep_cols, drop = TRUE])
  }
  pairs <- pairs %>%
    mutate(
      repA_vec = purrr::map2(Participant_id, tp_A, \(pid, tp) get_rep_vec(pid, tp)),
      repB_vec = purrr::map2(Participant_id, tp_B, \(pid, tp) get_rep_vec(pid, tp)),
      replicon_jaccard = if (length(rep_cols)) purrr::map2_dbl(repA_vec, repB_vec, jaccard_distance01) else NA_real_
    )
} else {
  pairs$replicon_jaccard <- NA_real_
}

# ------------------------------ 10 · merge nucmer/dnadiff metrics (optional) ---
if (exists("nucmer_tbl") && nrow(nucmer_tbl)) {
  nucmer_summ <- nucmer_tbl %>%
    group_by(Participant_id, tp_A, tp_B) %>%
    summarise(
      ANI  = suppressWarnings(max(AvgIdentity, na.rm = TRUE)),
      SNPs = suppressWarnings(min(TotalSNPs, na.rm = TRUE)),
      .groups = "drop"
    ) %>%
    mutate(
      ANI  = ifelse(is.infinite(ANI),  NA_real_, ANI),
      SNPs = ifelse(is.infinite(SNPs), NA_real_, SNPs)
    )
  pairs <- pairs %>%
    left_join(nucmer_summ, by = c("Participant_id","tp_A","tp_B"))
} else {
  pairs$ANI  <- NA_real_
  pairs$SNPs <- NA_real_
}

# ------------------------------ 11 · classify persistence vs replacement -------
pairs <- pairs %>%
  mutate(
    ST_same = !is.na(ST_A) & !is.na(ST_B) & ST_A == ST_B,
    vf_ok   = !is.na(vf_jaccard) & vf_jaccard <= VF_JAC_MAX,
    ani_ok  = is.na(ANI)  | (ANI  >= ANI_MIN),
    snp_ok  = is.na(SNPs) | (SNPs <= SNP_MAX),
    call    = dplyr::case_when(
      ST_same & vf_ok & ani_ok & snp_ok ~ "Persistent",
      !ST_same                           ~ "Replacement",
      ST_same & !vf_ok                   ~ "Replacement",
      ST_same & vf_ok & (!ani_ok | !snp_ok) ~ "Replacement",
      TRUE                               ~ "Uncertain"
    )
  )

# ------------------------------ 12 · participant-level summaries ---------------
persistence_by_participant <- pairs %>%
  group_by(Participant_id) %>%
  summarise(
    n_pairs = n(),
    n_persistent = sum(call == "Persistent", na.rm = TRUE),
    n_replacement = sum(call == "Replacement", na.rm = TRUE),
    n_uncertain = sum(call == "Uncertain", na.rm = TRUE),
    prop_persistent = ifelse(n_pairs > 0, n_persistent / n_pairs, NA_real_),
    .groups = "drop"
  )

# ------------------------------ 13 · write pairwise summary --------------------
pairwise_summary <- pairs %>%
  select(Participant_id, tp_A, tp_B, tpA_num, tpB_num,
         ST_A, ST_B, ST_same,
         vf_jaccard, replicon_jaccard, ANI, SNPs, call)

safe_write_csv(pairwise_summary, file.path(out_dir, "pairwise_summary.csv"))
safe_write_csv(persistence_by_participant, file.path(out_dir, "persistence_by_participant.csv"))

# ------------------------------ 14 · benchmark within vs between ---------------
# Build a between-participant baseline using VF Jaccard for pairs at the same
# timepoint. We compute explicit pairwise combinations to keep it stable.

vf_for_bench <- vf_mat %>% select(Participant_id, tp_lab, all_of(gene_cols))

# helper: compute all unordered pairwise Jaccard distances within a single tp_lab
.compute_between_for_tp <- function(df_tp) {
  # df_tp has columns: Participant_id, tp_lab, <genes...>
  if (nrow(df_tp) < 2) return(tibble::tibble())
  pids <- df_tp$Participant_id
  tpv  <- unique(as.character(df_tp$tp_lab))
  tpv  <- if (length(tpv)) tpv[1] else NA_character_
  
  mat  <- as.matrix(df_tp[, gene_cols, drop = FALSE])
  idx  <- t(combn(seq_len(nrow(mat)), 2))
  
  jd <- apply(idx, 1, function(ix) {
    jaccard_distance01(mat[ix[1], ], mat[ix[2], ])
  })
  
  tibble::tibble(
    kind       = "between",
    tp_lab     = tpv,
    pid_A      = pids[idx[, 1]],
    pid_B      = pids[idx[, 2]],
    vf_jaccard = as.numeric(jd)
  ) %>% filter(!is.na(vf_jaccard))
}

within_dist <- pairwise_summary %>%
  transmute(kind = "within", vf_jaccard = vf_jaccard) %>%
  filter(!is.na(vf_jaccard))

between_dist_all <- vf_for_bench %>%
  group_by(tp_lab) %>% group_split() %>%
  purrr::map_dfr(.compute_between_for_tp)

# size-match the between distribution to the number of within pairs (if needed)
set.seed(1)
n_within <- nrow(within_dist)
between_dist <- if (nrow(between_dist_all) > n_within) {
  dplyr::slice_sample(between_dist_all, n = n_within) %>% transmute(kind, vf_jaccard)
} else {
  between_dist_all %>% transmute(kind, vf_jaccard)
}

benchmark_within_vs_between <- dplyr::bind_rows(within_dist, between_dist)
safe_write_csv(benchmark_within_vs_between,
               file.path(out_dir, "benchmark_within_vs_between.csv"))

# ------------------------------ 15 · plots: VF Jaccard distributions -----------
p_vf_jacc <- ggplot(benchmark_within_vs_between,
                    aes(x = vf_jaccard, fill = kind)) +
  geom_density(alpha = 0.45, adjust = 1.2) +
  theme_bw(base_size = 12) +
  labs(title = "VF Jaccard distance: within-person vs between-person",
       x = "Jaccard distance (0 = identical, 1 = disjoint)",
       y = "Density", fill = "") +
  guides(fill = guide_legend(override.aes = list(alpha = 0.6)))

# rename existing density save (change your earlier save line)
safe_save_plot(file.path(plot_dir, "vf_jaccard_within_vs_between_density.png"),
               p_vf_jacc, width = 7, height = 4.5, dpi = 300)

## Boxplot/jitter for within vs between with Wilcoxon p-value (uses 'kind' & 'vf_jaccard')
bench_df <- readr::read_csv(file.path(out_dir, "benchmark_within_vs_between.csv"), show_col_types = FALSE)

wilcox_p <- NA_real_
if (sum(bench_df$kind == "within", na.rm = TRUE) >= 5 &&
    sum(bench_df$kind == "between", na.rm = TRUE) >= 5) {
  wilcox_p <- suppressWarnings(
    wilcox.test(bench_df$vf_jaccard[bench_df$kind == "within"],
                bench_df$vf_jaccard[bench_df$kind == "between"],
                alternative = "less")$p.value
  )
}
subtitle_txt <- if (!is.na(wilcox_p)) paste0("Wilcoxon (Within < Between) p = ", signif(wilcox_p, 3)) else "Wilcoxon p = NA"

pD <- bench_df %>%
  ggplot(aes(kind, vf_jaccard, fill = kind)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.15, alpha = 0.3) +
  labs(title = "Benchmark VF distances", subtitle = subtitle_txt, x = NULL, y = "VF Jaccard", fill = "")
safe_save_plot(file.path(plot_dir, "benchmark_within_vs_between.png"), pD, 6, 4.5)
# (optional) drop a text file like the initial script did
writeLines(paste0("wilcox_p (Within < Between) = ",
                  ifelse(is.na(wilcox_p), "NA", signif(wilcox_p, 6))),
           file.path(out_dir, "benchmark_wilcoxon.txt"))

# ------------------------------ 16 · plots: call summary -----------------------
p_call <- pairs %>%
  count(call, name = "n") %>%
  mutate(call = factor(call, levels = c("Persistent","Replacement","Uncertain"))) %>%
  ggplot(aes(call, n)) + geom_col() + theme_bw(base_size = 12) +
  labs(title = "Pairwise calls", x = NULL, y = "Count")

safe_save_plot(file.path(plot_dir, "pairwise_calls_bar.png"),
               p_call, width = 5, height = 4, dpi = 300)

# ------------------------------ 17 · optional: scatter VF vs ANI ---------------
if (any(!is.na(pairs$ANI))) {
  p_scatter <- pairs %>%
    ggplot(aes(x = vf_jaccard, y = ANI, shape = call)) +
    geom_point(size = 2) +
    theme_bw(base_size = 12) +
    labs(title = "Within-person pairs: VF Jaccard vs ANI",
         x = "VF Jaccard distance",
         y = "ANI (%)",
         shape = "Call")
  safe_save_plot(file.path(plot_dir, "vf_vs_ani_scatter.png"),
                 p_scatter, width = 6, height = 4.5, dpi = 300)
}

# A) Distribution of VF Jaccard by classification (matches initial script idea)
pA <- pairs %>%
  filter(!is.na(vf_jaccard)) %>%
  ggplot(aes(vf_jaccard, fill = call)) +
  geom_histogram(bins = 30, alpha = 0.8, position = "identity") +
  geom_vline(xintercept = VF_JAC_MAX, linetype = "dashed") +
  labs(title = "Within-person VF Jaccard distances",
       subtitle = paste0("Dashed line = threshold ", VF_JAC_MAX),
       x = "VF Jaccard distance (0 = identical, 1 = disjoint)",
       y = "Pairs", fill = "Call")
safe_save_plot(file.path(plot_dir, "vf_jaccard_by_classification.png"), pA, 8, 5)

# B) Trajectories of VF distance (consecutive TPs) — requires tp_num per TP
# Build a reliable PID×TP -> tp_num lookup (avoid relying on 'samples')
tp_lookup <- vf %>%
  dplyr::select(Participant_id, tp_lab, tp_num) %>%
  dplyr::distinct()

# B) Trajectories of VF distance (consecutive TPs)
traj_df <- pairs %>%
  dplyr::left_join(tp_lookup, by = c("Participant_id","tp_A" = "tp_lab")) %>%
  dplyr::rename(tp_num_A = tp_num) %>%
  dplyr::left_join(tp_lookup, by = c("Participant_id","tp_B" = "tp_lab")) %>%
  dplyr::rename(tp_num_B = tp_num) %>%
  dplyr::filter(!is.na(tp_num_A), !is.na(tp_num_B))

pC <- traj_df %>%
  ggplot2::ggplot(ggplot2::aes(tp_num_A, vf_jaccard, group = Participant_id, color = call)) +
  ggplot2::geom_line(alpha = 0.4) +
  ggplot2::geom_point() +
  ggplot2::scale_x_continuous(breaks = sort(unique(traj_df$tp_num_A))) +
  ggplot2::labs(title = "Trajectories of within-person VF distance (consecutive TPs)",
                x = "Numeric timepoint (A)", y = "VF Jaccard")
safe_save_plot(file.path(plot_dir, "vf_trajectory_consecutive.png"), pC, 8, 5)


# ------------------------------ 18 · write extra debug aids --------------------
# (These are handy when diagnosing mismatches or empty gene change tables.)
safe_write_csv(samples, file.path(out_dir, "sample_key.csv"))
safe_write_csv(pairs,    file.path(out_dir, "pairs_with_metrics.csv"))

msg("All done. Saved benchmark table and plots in results/within_person/ and results/plots/within_person/")

# Tiny README
readme_lines <- c(
  "Within-person strain comparison (Yellow RoUTIne):",
  paste0("- Persistence rule: ST_same & VF_Jaccard ≤ ", VF_JAC_MAX,
         if (exists("nucmer_tbl") && nrow(nucmer_tbl))
           paste0(" & ANI ≥ ", ANI_MIN, "% (or SNPs ≤ ", SNP_MAX, ")")
         else ""),
  "- Missing core metrics or discordant thresholds → 'Uncertain'.",
  "- See pairwise_summary.csv for per-pair metrics; debug/ for conflicts/missing data.",
  "- Optional: replicon Jaccard is computed only if Isolate_ID→Participant×Timepoint mapping is available."
)
writeLines(readme_lines, file.path(out_dir, "README_11_within_person.txt"))

# Verification table
verify_file <- function(p) file.exists(p) && isTRUE(file.size(p) > 0)
must_have <- c(file.path(out_dir, "pairwise_summary.csv"),
               file.path(out_dir, "persistence_by_participant.csv"),
               file.path(out_dir, "gene_changes_long.csv"))
ok_vec <- vapply(must_have, verify_file, logical(1))
verif_tbl <- tibble::tibble(file = must_have, exists_nonempty = ok_vec)
safe_write_csv(verif_tbl, file.path(out_dir, "verification_table.csv"))

# Session info snapshot
sess_txt <- capture.output(sessionInfo())
writeLines(sess_txt, file.path(out_dir, "sessionInfo.txt"))
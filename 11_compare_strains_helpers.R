#!/usr/bin/env Rscript
# =============================================================
# 11_compare_strains_helpers.R
# Shared helpers for strain comparison across participants/timepoints
# Yellow RoUTIne – production-grade utilities
# =============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
  library(stringr)
  library(purrr)
  library(tibble)
})

# ------------- timepoint normalization ---------------------------------------
tp_norm <- function(x) {
  tp_chr <- as.character(x)
  is_uricult <- stringr::str_detect(tp_chr, stringr::regex("uricult", ignore_case = TRUE))
  tp_num <- suppressWarnings(as.integer(stringr::str_extract(tp_chr, "\\d+")))
  tp_lab <- dplyr::case_when(
    is_uricult ~ "Uricult",
    !is.na(tp_num) ~ paste0("T", tp_num),
    TRUE ~ "Unscheduled"
  )
  tp_levels <- c(
    paste0("T", sort(unique(tp_num[!is.na(tp_num)]))),
    "Uricult", "Unscheduled"
  )
  tibble::tibble(tp_lab = factor(tp_lab, levels = tp_levels), tp_num = tp_num)
}

# ------------- IO helpers ----------------------------------------------------
safe_dir_create <- function(path) {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  invisible(path)
}

safe_write_csv <- function(x, path, ...) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  tmp <- paste0(path, ".tmp")
  readr::write_csv(x, tmp, ...)
  if (file.exists(path)) unlink(path)
  file.rename(tmp, path)
  invisible(path)
}

timestamp_msg <- function(...) {
  cat(format(Sys.time(), "[%Y-%m-%d %H:%M:%S] "))
  message(...)
}

# ------------- tool detection ------------------------------------------------
has_tool <- function(bin) nzchar(Sys.which(bin))

detect_tools <- function() {
  list(
    nucmer = has_tool("nucmer"),
    dnadiff = has_tool("dnadiff"),
    `delta-filter` = has_tool("delta-filter"),
    mash = has_tool("mash")
  )
}

# ------------- load core tables ---------------------------------------------
load_core_tables <- function() {
  # assemblies
  if (!file.exists("assembly_metadata.csv")) {
    stop("assembly_metadata.csv not found – run 01_prepare_assembly_metadata.R")
  }
  asm <- readr::read_csv("assembly_metadata.csv", show_col_types = FALSE)
  # ensure tp_lab
  if (!"tp_lab" %in% names(asm)) {
    if (!"Timepoint" %in% names(asm)) stop("assembly_metadata.csv lacks Timepoint column")
    asm <- dplyr::bind_cols(asm, tp_norm(asm$Timepoint))
  }
  # ensure full_path
  if (!"full_path" %in% names(asm)) {
    fasta_dir <- "ont-yellow-routine-fastas"
    if (!dir.exists(fasta_dir)) stop("FASTA dir not found: ", fasta_dir)
    if (!"file_name" %in% names(asm)) stop("assembly_metadata.csv lacks file_name for reconstructing full_path")
    asm <- asm %>% mutate(full_path = file.path(fasta_dir, file_name), found = file.exists(full_path))
    if (any(!asm$found)) {
      bad <- asm$full_path[!asm$found]
      stop("Missing FASTA files (", length(bad), "):\n", paste(bad, collapse = "\n"))
    }
    asm <- select(asm, -found)
  }

  # status map (optional)
  status_map <- NULL
  if (file.exists("results/status_map.csv")) {
    status_map <- readr::read_csv("results/status_map.csv", show_col_types = FALSE)
    if (!"tp_lab" %in% names(status_map)) {
      if (!"Timepoint" %in% names(status_map)) stop("status_map.csv must contain 'tp_lab' or 'Timepoint'")
      status_map <- dplyr::bind_cols(status_map, tp_norm(status_map$Timepoint))
    }
    status_map <- status_map %>%
      select(Participant_id, tp_lab, Infection_Status) %>%
      distinct()
  }

  # MLST
  if (!file.exists("results/mlst/mlst_all.tsv")) {
    stop("results/mlst/mlst_all.tsv not found – run 06_MLST.R")
  }
  mlst <- readr::read_tsv("results/mlst/mlst_all.tsv", show_col_types = FALSE)
  st_col <- names(mlst)[tolower(names(mlst)) == "st"]
  if (!length(st_col)) {
    alt <- names(mlst)[grepl("(?i)^sequence[_ ]?type$", names(mlst))]
    if (length(alt)) st_col <- alt[1]
  }
  if (!length(st_col)) stop("Could not find ST column in mlst_all.tsv")
  if (!"ST" %in% names(mlst)) mlst$ST <- mlst[[st_col[1]]]
  # keep only relevant columns
  mlst <- mlst %>% select(Isolate_ID = any_of(c("Isolate_ID", "isolate_id")), ST, everything())

  # pMLST (wide)
  pmlst_file <- "results/mlst/plasmid_types_per_isolate.csv"
  pmlst_wide <- NULL
  if (file.exists(pmlst_file)) {
    pmlst_wide <- readr::read_csv(pmlst_file, show_col_types = FALSE)
    # ensure Isolate_ID column case
    if (!"Isolate_ID" %in% names(pmlst_wide)) {
      iso <- names(pmlst_wide)[grepl("(?i)^isolate[_ ]?id$", names(pmlst_wide))][1]
      if (!is.na(iso)) names(pmlst_wide)[names(pmlst_wide) == iso] <- "Isolate_ID"
    }
  }

  # VF presence/absence
  if (!file.exists(FILE_VF_PA)) {
    stop(FILE_VF_PA, " not found – run 02_gene_presence_analysis.R")
  }
  vf_pa <- readr::read_csv(FILE_VF_PA, show_col_types = FALSE)
  if (!"tp_lab" %in% names(vf_pa)) {
    if (!"Timepoint" %in% names(vf_pa)) stop("vf_pa_all.csv must have 'tp_lab' or 'Timepoint'")
    vf_pa <- dplyr::bind_cols(vf_pa, tp_norm(vf_pa$Timepoint))
  }
  vf_pa <- vf_pa %>% mutate(SampleKey = paste(Participant_id, as.character(tp_lab), sep = "__"))

  # replicon presence/absence (plasmidfinder)
  inc_pa <- NULL
  if (file.exists("results/plasmidfinder_presence_absence.csv")) {
    inc_pa <- readr::read_csv("results/plasmidfinder_presence_absence.csv", show_col_types = FALSE)
    # ensure Isolate_ID column exists
    if (!"Isolate_ID" %in% names(inc_pa)) {
      iso <- names(inc_pa)[grepl("(?i)^isolate[_ ]?id$", names(inc_pa))][1]
      if (!is.na(iso)) names(inc_pa)[names(inc_pa) == iso] <- "Isolate_ID"
    }
  }

  list(
    assemblies = asm,
    status_map = status_map,
    mlst = mlst,
    vf_pa = vf_pa,
    inc_pa = inc_pa,
    pmlst_wide = pmlst_wide
  )
}

# ------------- sample resolution --------------------------------------------
resolve_sample <- function(Participant_id, tp_lab, assemblies, prefer_assembler = "flye") {
  stopifnot("Participant_id" %in% names(assemblies), "tp_lab" %in% names(assemblies))
  df <- assemblies %>% filter(
    Participant_id == !!Participant_id,
    as.character(tp_lab) == as.character(!!tp_lab)
  )
  if (!nrow(df)) {
    return(NULL)
  }
  # choose assembler: prefer -> longcycler -> any; break ties by max total_bases
  df <- df %>%
    mutate(
      assembler = as.character(assembler),
      pref_rank = dplyr::case_when(
        assembler == prefer_assembler ~ 1L,
        assembler == "longcycler" ~ 2L,
        TRUE ~ 3L
      ),
      size_rank = dplyr::desc(coalesce(as.numeric(total_bases), 0))
    ) %>%
    arrange(pref_rank, size_rank)
  df[1, , drop = FALSE]
}

# ------------- pair helpers --------------------------------------------------
make_sample_key <- function(pid, tp) paste(pid, as.character(tp), sep = "__")

# ------------- similarity metrics -------------------------------------------
parse_dnadiff_report <- function(report_path) {
  if (!file.exists(report_path)) {
    return(tibble(AvgIdentity = NA_real_, TotalSNPs = NA_real_))
  }
  L <- readLines(report_path, warn = FALSE)
  get_num <- function(key) {
    m <- grep(key, L, value = TRUE)
    if (!length(m)) {
      return(NA_real_)
    }
    as.numeric(stringr::str_extract(m[1], "\\d+\\.?\\d*"))
  }
  tibble(AvgIdentity = get_num("AvgIdentity"), TotalSNPs = get_num("TotalSNPs"))
}

run_dnadiff <- function(a_fasta, b_fasta, cache_dir, key) {
  safe_dir_create(cache_dir)
  if (!has_tool("dnadiff")) {
    return(tibble(AvgIdentity = NA_real_, TotalSNPs = NA_real_))
  }
  # sanitize key for filesystem; avoid regex escaping issues by placing '-' at end of class
  pref <- file.path(cache_dir, paste0(gsub("[^A-Za-z0-9_-]+", "_", key)))
  rpt <- paste0(pref, ".report")
  if (!file.exists(rpt)) {
    cmd <- sprintf("dnadiff -p %s %s %s", shQuote(pref), shQuote(a_fasta), shQuote(b_fasta))
    status <- system(cmd)
    if (status != 0) {
      return(tibble(AvgIdentity = NA_real_, TotalSNPs = NA_real_))
    }
  }
  parse_dnadiff_report(rpt)
}

mash_distance <- function(a_fasta, b_fasta) {
  if (!has_tool("mash")) {
    return(NA_real_)
  }
  cmd <- sprintf("mash dist %s %s", shQuote(a_fasta), shQuote(b_fasta))
  out <- suppressWarnings(system(cmd, intern = TRUE))
  if (!length(out)) {
    return(NA_real_)
  }
  # typical format: a\tb\tdist\tvar\tshared-hashes
  parts <- strsplit(out[1], "\t")[[1]]
  if (length(parts) < 3) {
    return(NA_real_)
  }
  as.numeric(parts[3])
}

booleanize_df <- function(df) {
  df[] <- lapply(df, function(col) {
    if (is.logical(col)) {
      return(as.integer(col))
    }
    if (is.numeric(col)) {
      return(as.integer(col > 0))
    }
    # character: treat "1"/"0" as ints, else 0/1
    suppressWarnings(as.integer(col)) %>% tidyr::replace_na(0L)
  })
  as.data.frame(df)
}

jaccard_from_wide <- function(df_wide, row_key_cols = c("Participant_id", "tp_lab"), sampleA, sampleB, cols_regex = NULL) {
  stopifnot(all(row_key_cols %in% names(df_wide)))
  dat <- df_wide
  if (!is.null(cols_regex)) {
    keep <- grep(cols_regex, names(dat), value = TRUE)
    dat <- dat[, c(row_key_cols, keep), drop = FALSE]
  }
  # ensure one row per sample key
  dat$SampleKey <- make_sample_key(dat[[row_key_cols[1]]], dat[[row_key_cols[2]]])
  A <- dat %>% filter(SampleKey == sampleA)
  B <- dat %>% filter(SampleKey == sampleB)
  if (!nrow(A) || !nrow(B)) {
    return(list(jaccard = NA_real_, n_int = NA_integer_, n_union = NA_integer_))
  }
  mat <- rbind(A, B) %>%
    select(-all_of(c(row_key_cols, "SampleKey"))) %>%
    booleanize_df() %>%
    as.matrix()
  if (!ncol(mat)) {
    return(list(jaccard = NA_real_, n_int = 0L, n_union = 0L))
  }
  a <- as.integer(mat[1, ] > 0)
  b <- as.integer(mat[2, ] > 0)
  inter <- sum(a & b)
  un <- sum(a | b)
  list(jaccard = if (un == 0) NA_real_ else inter / un, n_int = inter, n_union = un)
}

# ------------- classification -------------------------------------------------
classify_pair <- function(metrics, thresholds = list(id = 99.9, snps = 50, vf = 0.9, inc = 0.8, vf_rel = 0.7, inc_rel = 0.7)) {
  ST_equal <- isTRUE(metrics$ST_equal)
  has_id <- !is.na(metrics$AvgIdentity) && !is.na(metrics$TotalSNPs)
  mean_acc <- mean(c(metrics$VF_Jaccard, metrics$Inc_Jaccard), na.rm = TRUE)

  # Helper: treat NA as passing (1.0) for Jaccard (implies empty union -> identical in emptiness)
  pass_j <- function(val, thresh) is.na(val) || val >= thresh

  # Same strain criteria
  same_rule <- ST_equal && (
    (!has_id || (metrics$AvgIdentity >= thresholds$id && metrics$TotalSNPs <= thresholds$snps))
  ) && pass_j(metrics$VF_Jaccard, thresholds$vf) && pass_j(metrics$Inc_Jaccard, thresholds$inc)

  if (same_rule) {
    return(list(Classification = "Same", RuleUsed = "ST + ID/SNP + accessory"))
  }

  # Related criteria
  related_rule <- ST_equal && (
    (has_id && metrics$AvgIdentity >= (thresholds$id - 0.9)) || (!is.na(metrics$MashDistance) && metrics$MashDistance <= 0.02) || (!is.na(mean_acc) && mean_acc >= 0.7)
  )
  if (related_rule) {
    return(list(Classification = "Related", RuleUsed = "ST + partial concordance"))
  }

  # Fallbacks when identity unavailable
  if (!has_id && ST_equal && !is.na(mean_acc)) {
    if (mean_acc >= 0.9) {
      return(list(Classification = "Same", RuleUsed = "ST + accessory >=0.9 (no ID)"))
    }
    if (mean_acc >= 0.7) {
      return(list(Classification = "Related", RuleUsed = "ST + accessory >=0.7 (no ID)"))
    }
  }
  list(Classification = "Different", RuleUsed = "Default")
}

# ------------- plotting helpers ---------------------------------------------
ggsave_safe <- function(filename, plot, width = 7, height = 5, dpi = 300) {
  dir.create(dirname(filename), recursive = TRUE, showWarnings = FALSE)
  if (grepl("\\.png$", filename, ignore.case = TRUE) && requireNamespace("ragg", quietly = TRUE)) {
    ggplot2::ggsave(filename, plot, device = ragg::agg_png, width = width, height = height, dpi = dpi, units = "in")
  } else {
    ggplot2::ggsave(filename, plot, width = width, height = height, dpi = dpi, units = "in")
  }
}

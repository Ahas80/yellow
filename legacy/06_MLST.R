#!/usr/bin/env Rscript
# =============================================================
# 06_MLST.R  –  Multi-locus sequence typing for “Yellow RoUTIne”
# Works with mlst ≥ 2.19   ( --quiet --csv --legacy )
# Last updated: 9 Jul 2025
# =============================================================

# ---- 0 · configuration ------------------------------------------------------
THREADS   <- pmin(4L, max(1L, parallel::detectCores() - 1L))   # PubMLST ≤ 4
FASTA_DIR <- "/Users/Aamir/Desktop/rUTIs/ont-yellow-routine-fastas"  # ← adjust
SCHEME    <- "ecoli"                      # PubMLST scheme, or "auto"
DEBUG     <- FALSE                        # TRUE → echo every mlst command

if (DEBUG) { options(warn = 1, error = recover); msg <- message } else {
  msg <- function(...) message(format(Sys.time(), "%H:%M:%S "), ...)
}

suppressPackageStartupMessages({
  library(dplyr);  library(readr);  library(tidyr)
  library(purrr);  library(furrr);  library(processx)
  library(stringr);library(fs);     library(scales); library(rlang)
})

msg("06_MLST.R starting")

# ---------- 1 · read assembly metadata ---------------------------------------
assembly_df <- read_csv("assembly_metadata.csv", show_col_types = FALSE)

if (!"full_path" %in% names(assembly_df)) {
  assembly_df <- assembly_df %>%
    mutate(full_path = file.path(FASTA_DIR, file_name)) %>%
    mutate(found     = file.exists(full_path))
  if (any(!assembly_df$found))
    stop("Missing FASTA files:\n",
         paste(assembly_df$full_path[!assembly_df$found], collapse = "\n"))
  assembly_df <- select(assembly_df, -found)
  msg("✓ full_path rebuilt (", nrow(assembly_df), " rows)")
}

# ---------- 2 · locate mlst & sanity-check -----------------------------------
MLST_BIN <- Sys.which("mlst")
if (MLST_BIN == "") stop("`mlst` not in $PATH – install or relink Homebrew")
if (system2(MLST_BIN, "--check") != 0L)
  stop("mlst --check failed – see stderr above")
msg("✓ mlst binary = ", MLST_BIN)

# ---------- 3 · output folders ------------------------------------------------
dir_create("results/mlst/raw",  recurse = TRUE)
dir_create("results/mlst/logs", recurse = TRUE)

# ---------- 4 · helpers -------------------------------------------------------
inspect_fail_logs <- function(n = 5) {
  fails <- fs::dir_ls("results/mlst/logs", glob = "*.log.txt")
  if (!length(fails)) { message("No failing logs found."); return(invisible()) }
  cat("\nShowing first", n, "failing logs:\n")
  for (lf in head(fails, n)) {
    cat("\n======", fs::path_file(lf), "======\n")
    cat(readLines(lf, n = 40), sep = "\n")
  }
  invisible()
}

run_mlst <- function(fasta, scheme = SCHEME) {
  cache   <- fs::path("results/mlst/raw",  paste0(fs::path_file(fasta), ".mlst.csv"))
  log_out <- fs::path("results/mlst/logs", paste0(fs::path_file(fasta), ".log.txt"))
  
  if (!fs::file_exists(cache)) {
    cmd <- c("--quiet", "--threads", "1",
             "--scheme", scheme, "--csv", "--legacy", fasta)
    if (DEBUG) msg("cmd: ", MLST_BIN, " ", paste(cmd, collapse = " "))
    
    px <- processx::run(MLST_BIN, cmd,
                        echo             = DEBUG,
                        stderr_to_stdout = TRUE,
                        error_on_status  = FALSE)
    
    write_lines(px$stdout, log_out)
    if (px$status != 0) return(tibble(.status = "failed"))
    write_lines(px$stdout, cache)
  }
  
  dat <- read_csv(cache,
                  col_types = cols(.default = col_character()),
                  na        = c("", "?"),
                  progress  = FALSE) %>%
    rename_with(tolower) %>%
    mutate(scheme = scheme, .before = 1)
  
  if (nrow(dat) == 0)
    return(tibble(.status = "no-data"))
  
  st_pattern <- regex("^st(\\.\\.[0-9]+)?$", ignore_case = TRUE)
  if (!any(str_detect(names(dat), st_pattern)))
    return(tibble(.status = "unexpected-format"))
  
  mutate(dat, .status = "ok")
}

# ---------- 5 · run in parallel ----------------------------------------------
plan(multisession, workers = THREADS)

mlst_raw <- assembly_df %>%
  mutate(mlst = future_map(full_path, run_mlst, .progress = TRUE)) %>%
  unnest(mlst)

# ---------- 6 · concise overview ---------------------------------------------
summary_tbl <- mlst_raw %>%
  count(.status, name = "n") %>%
  mutate(pct = percent(n / sum(n)))

msg("------ run summary ------")
print(summary_tbl)
if (!any(mlst_raw$.status == "ok"))
  stop("No successful typings – call inspect_fail_logs()")

# ---------- 7 · final tidy-up -------------------------------------------------
st_pattern <- regex("^st(\\.\\.[0-9]+)?$", ignore_case = TRUE)
st_col     <- names(mlst_raw)[str_detect(names(mlst_raw), st_pattern)]

if (!length(st_col))
  stop("Could not find an ST column – names are:\n",
       paste(names(mlst_raw), collapse = ", "))

mlst_tbl <- mlst_raw %>%
  filter(.status == "ok") %>%
  select(-.status) %>%
  mutate(ST = .data[[st_col[1]]]) %>%        # copy → ST
  select(-all_of(st_col)) %>%                # drop original
  relocate(Isolate_ID, ST, everything())

# ---------- 8 · write outputs -------------------------------------------------
write_tsv(mlst_tbl,    "results/mlst/mlst_all.tsv")
write_csv(mlst_tbl,    "results/mlst/mlst_matrix.csv")
write_csv(summary_tbl, "results/mlst/log_summary.csv")

msg("✓ wrote ", nrow(mlst_tbl), " assemblies → results/mlst/")
msg("06_MLST.R finished ✅")

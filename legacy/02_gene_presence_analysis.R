#!/usr/bin/env Rscript
# =============================================================
# 02_gene_presence_analysis.R
# -------------------------------------------------------------
#  * Expects `assembly_df` in the workspace  – OR reads the CSV
#  * Runs Abricate (VFDB) on every FASTA (cached to /results/abricate)
#  * Produces presence/absence matrices, prevalence tables & plots
#  * Optional nucmer trajectories for participants with ≥2 time-points
# =============================================================
## ---------------- 1 ·  libraries  -------------------------------------------
suppressPackageStartupMessages({
  library(dplyr);  library(tidyr);  library(readr)
  library(purrr);  library(furrr);  library(stringr)
  library(ggplot2); library(forcats); library(glue)
})

## ---------------- 2 ·  load assembly_df  ------------------------------------
if (!exists("assembly_df")) {
  if (file.exists("assembly_metadata.csv")) {
    assembly_df <- read_csv("assembly_metadata.csv", show_col_types = FALSE)
    message("✓  loaded assembly_df from assembly_metadata.csv  (",
            nrow(assembly_df), " rows)")
  } else {
    stop("assembly_df not found – run 01_prepare_assembly_metadata.R first")
  }
}
## ---- 2b · ensure full_path exists  -----------------------------
if (!"full_path" %in% names(assembly_df)) {
  if (!"file_name" %in% names(assembly_df))
    stop("assembly_df lacks both full_path and file_name columns.")
  
  fasta_dir <- "/Users/Aamir/Desktop/rUTIs/ont-yellow-routine-fastas"   # <-- adjust
  assembly_df <- assembly_df %>%
    mutate(full_path = file.path(fasta_dir, file_name)) %>%
    mutate(found = file.exists(full_path))
  
  if (any(!assembly_df$found)) {
    bad <- assembly_df %>% filter(!found) %>% pull(full_path)
    stop("Missing FASTA files:\n", paste(bad, collapse = "\n"))
  }
  assembly_df <- select(assembly_df, -found)  # clean-up helper col
  message("✓  full_path column rebuilt from file_name (", nrow(assembly_df), " rows)")
}
## ---- 2c · normalize timepoints on assemblies (tp_lab + tp_num) --------------
# Drop-in replacement for tp_norm() in 02_gene_presence_analysis.R
# --- replace your current 2c block ---
## ---- 2c · normalize timepoints on assemblies (tp_lab + tp_num) --------------
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
  tp_lab <- factor(tp_lab,
                   levels = c(paste0("T", sort(unique(tp_num[!is.na(tp_num)]))),
                              "Uricult","Unscheduled"))
  tibble::tibble(tp_lab = tp_lab, tp_num = tp_num)
}

# Overwrite clean tp_* columns (avoid ... suffixes)
tp <- tp_norm(assembly_df$Timepoint)
assembly_df <- assembly_df %>%
  dplyr::mutate(tp_lab = tp$tp_lab, tp_num = tp$tp_num)

# Drop any leftover suffixed duplicates from previous runs
dup_tp <- grep("^tp_(lab|num)\\.\\.\\.[0-9]+$", names(assembly_df), value = TRUE)
if (length(dup_tp)) {
  assembly_df <- dplyr::select(assembly_df, -tidyselect::all_of(dup_tp))
}

# (Do NOT touch vf_pa_all here; it's not built yet)



# Drop ANY pre-existing tp_* columns (including ...##) before adding fresh ones


## ---------------- 3 ·  global output dirs  ----------------------------------
dir.create("results/abricate",   recursive = TRUE, showWarnings = FALSE)
dir.create("results/plots",      recursive = TRUE, showWarnings = FALSE)

## ---------------- 4 ·  Abricate-VFDB scan  -----------------------------------
run_abr_cached <- function(fasta, db = "vfdb", min_cov = 70, min_id = 70) {
  cache <- file.path("results/abricate", paste0(basename(fasta), ".vfdb.tsv"))
  if (file.exists(cache))
    return(readr::read_tsv(cache, show_col_types = FALSE, progress = FALSE))
  
  cmd <- glue::glue("abricate --quiet --db {db} --mincov {min_cov} --minid {min_id} {shQuote(fasta)} > {shQuote(cache)}")
  exit <- system(cmd)
  if (exit != 0) warning("Abricate non-zero exit (", basename(fasta), ")")
  readr::read_tsv(cache, show_col_types = FALSE, progress = FALSE)
}

future::plan(future::multisession, workers = max(1, parallel::detectCores() - 1))
safe_abr <- purrr::safely(run_abr_cached, otherwise = NULL, quiet = TRUE)

vf_hits_all <- assembly_df %>%
  dplyr::mutate(vfdb = furrr::future_map(full_path, ~ safe_abr(.x)$result, .progress = TRUE)) %>%
  dplyr::filter(purrr::map_lgl(vfdb, ~ !is.null(.x) && NROW(.x) > 0)) %>%
  tidyr::unnest(vfdb)

if (!"tp_lab" %in% names(vf_hits_all)) {
  tp_col <- grep("^tp_lab", names(vf_hits_all), value = TRUE)[1]
  if (length(tp_col)) vf_hits_all <- dplyr::rename(vf_hits_all, tp_lab = !!tp_col)
}

# Be resilient to different gene column names from abricate
gene_col <- intersect(c("GENE","GENE_NAME","NAME","PRODUCT","GENE SYMBOL"), names(vf_hits_all))[1]
if (is.na(gene_col)) stop("No gene name column found in Abricate output.")
vf_hits_all <- vf_hits_all %>% dplyr::rename(GENE = dplyr::all_of(gene_col))

saveRDS(vf_hits_all, "results/vf_hits_all.rds")
message("✓  vf_hits_all saved  (", nrow(vf_hits_all), " rows)")

## ---------------- 5 ·  presence/absence matrix ------------------------------
vf_pa_all <- vf_hits_all %>%
  dplyr::distinct(Participant_id, tp_lab, GENE) %>%   # keep many-to-many at sample level
  dplyr::mutate(present = 1) %>%
  tidyr::pivot_wider(names_from = GENE, values_from = present, values_fill = 0)

readr::write_csv(vf_pa_all, "results/vf_pa_all.csv")

## ---------------- 6 ·  gene-level prevalence --------------------------------
tbl_gene <- vf_hits_all %>%
  dplyr::distinct(Participant_id, GENE) %>%
  dplyr::count(GENE, name = "n_participants") %>%
  dplyr::arrange(dplyr::desc(n_participants))

readr::write_csv(tbl_gene, "results/stats_gene_level.csv")

## ---------------- 7 ·  quick cohort plots -----------------------------------
dir.create("results/plots", showWarnings = FALSE)

top25 <- tbl_gene %>% slice_max(n_participants, n = 25) %>%
  mutate(GENE = fct_reorder(GENE, n_participants))

ggplot(top25, aes(GENE, n_participants)) +
  geom_col(fill = "steelblue") + coord_flip() +
  labs(title = "Top 25 VFDB genes (entire cohort)",
       y = "Participants", x = NULL) +
  theme_minimal(base_size = 11)
ggsave("results/plots/core_bar_top25_all.png", width = 6, height = 6)

ggplot(tbl_gene, aes(n_participants)) +
  geom_histogram(binwidth = 1, fill = "grey70") +
  labs(title = "VF gene prevalence distribution",
       x = "# participants", y = "Gene count") +
  theme_minimal(base_size = 11)
ggsave("results/plots/core_histogram_all.png", width = 5, height = 4)

## ---------------- 8 ·  per-participant Heat-map & UpSet ---------------------
library(ComplexUpset)

variable_genes <- tbl_gene %>%
  dplyr::filter(dplyr::between(n_participants, 1, max(n_participants) - 1)) %>%
  dplyr::pull(GENE)

vf_pa_var <- vf_pa_all %>%
  dplyr::select(Participant_id, tp_lab, tidyselect::all_of(variable_genes))

plot_heatmap <- function(pid) {
  vf_pa_var %>%
    dplyr::filter(Participant_id == pid) %>%
    tidyr::pivot_longer(-c(Participant_id, tp_lab),
                        names_to = "GENE", values_to = "present") %>%
    ggplot2::ggplot(ggplot2::aes(GENE, tp_lab, fill = factor(present))) +
    ggplot2::geom_tile(color = "grey80") +
    ggplot2::scale_fill_manual(values = c(`0` = "white", `1` = "steelblue")) +
    ggplot2::labs(title = paste("Variable genes – participant", pid),
                  x = NULL, y = "Time-point", fill = "Present") +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1),
                   panel.grid = ggplot2::element_blank(),
                   legend.position = "none")
}

get_upset_df <- function(pid) {
  vf_pa_all %>%
    dplyr::filter(Participant_id == pid) %>%
    tidyr::pivot_longer(-c(Participant_id, tp_lab),
                        names_to = "GENE", values_to = "present") %>%
    dplyr::filter(present == 1) %>%
    tidyr::unite(GENE_TP, GENE, tp_lab, sep = "_") %>%
    dplyr::mutate(val = TRUE) %>%
    tidyr::pivot_wider(names_from = GENE_TP, values_from = val, values_fill = FALSE)
}

dir.create("results/plots/participants", showWarnings = FALSE)
purrr::walk(unique(vf_pa_all$Participant_id), \(pid) {
  ggplot2::ggsave(glue::glue("results/plots/participants/heatmap_{pid}.png"),
                  plot_heatmap(pid), width = 12, height = 6, dpi = 300)
  
  df <- get_upset_df(pid)
  ComplexUpset::upset(df,
                      intersect = setdiff(names(df), "Participant_id"),
                      min_size = 1, name = as.character(pid))
  ggplot2::ggsave(glue::glue("results/plots/participants/upset_{pid}.png"),
                  width = 12, height = 8, dpi = 300)
})

## ---------------- 9 ·  nucmer trajectories (≥2 numeric TPs) ------------------
ids_multi <- vf_pa_all %>%
  dplyr::distinct(Participant_id, tp_lab) %>%
  dplyr::mutate(tp_num = readr::parse_number(as.character(tp_lab))) %>%
  dplyr::filter(!is.na(tp_num)) %>%
  dplyr::count(Participant_id, name = "n_tp") %>%
  dplyr::filter(n_tp >= 2) %>%
  dplyr::pull(Participant_id)

if (length(ids_multi)) {
  message("→  nucmer on ", length(ids_multi), " multi-TP participants")
  
  assembly_long <- assembly_df %>%
    dplyr::filter(Participant_id %in% ids_multi) %>%
    dplyr::filter(!is.na(tp_num))   # numeric T0/T1/... only
  
  pair_tbl <- assembly_long %>%
    dplyr::group_by(Participant_id, assembler) %>%
    dplyr::arrange(tp_num, .by_group = TRUE) %>%
    dplyr::mutate(path_A = full_path,
                  path_B = dplyr::lead(full_path),
                  tp_A   = tp_lab,
                  tp_B   = dplyr::lead(tp_lab),
                  outdir = glue::glue("results/nucmer/{Participant_id}_{assembler}_{tp_A}_vs_{tp_B}")) %>%
    dplyr::filter(!is.na(path_B)) %>%
    dplyr::ungroup()
  
  run_nucmer <- function(a, b, od) {
    dir.create(od, recursive = TRUE, showWarnings = FALSE)
    pref <- file.path(od, "run")
    system(glue::glue("nucmer --mum --prefix={shQuote(pref)} {shQuote(a)} {shQuote(b)}"))
    system(glue::glue("delta-filter -1 {pref}.delta > {pref}.1delta"))
    system(glue::glue("dnadiff -p {pref}_dd -d {pref}.1delta"))
    rpt <- glue::glue("{pref}_dd.report")
    if (!file.exists(rpt)) return(tibble::tibble(AvgIdentity = NA_real_, TotalSnpCnt = NA_real_))
    L <- readLines(rpt)
    grab <- \(k) { m <- grep(k, L, value = TRUE); if (length(m)) as.numeric(stringr::str_extract(m[1], "\\d+\\.?\\d*")) else NA_real_ }
    tibble::tibble(AvgIdentity = grab("AvgIdentity"), TotalSnpCnt = grab("TotalSNPs"))
  }
  
  pair_df <- pair_tbl %>%
    dplyr::mutate(res = furrr::future_pmap(list(path_A, path_B, outdir), run_nucmer, .progress = TRUE)) %>%
    tidyr::unnest(res)
  
  dir.create("results/plots/pairwise_identity", showWarnings = FALSE)
  
  plot_pair <- function(pid, asm) {
    dat <- pair_df %>%
      dplyr::filter(Participant_id == pid, assembler == asm) %>%
      dplyr::arrange(tp_num)                 # keep temporal order
    
    if (!nrow(dat)) return(NULL)
    
    g <- ggplot2::ggplot(dat, ggplot2::aes(x = tp_A, y = AvgIdentity, colour = TotalSnpCnt))
    
    if (nrow(dat) > 1) {
      g <- g + ggplot2::geom_line(ggplot2::aes(group = 1), linewidth = 0.8)
    }
    
    g + ggplot2::geom_point(size = 3) +
      ggplot2::scale_colour_viridis_c(option = "D", name = "SNP count") +
      ggplot2::labs(title = glue::glue("Pairwise identity – P{pid} ({asm})"),
                    x = "Time-point A", y = "Avg identity (%)") +
      ggplot2::theme_minimal(base_size = 10)
  }
  
  pw <- pair_df %>% dplyr::distinct(Participant_id, assembler)
  pw %>% purrr::pwalk(\(Participant_id, assembler) {
    g <- plot_pair(Participant_id, assembler)
    if (!is.null(g))
      ggplot2::ggsave(glue::glue("results/plots/pairwise_identity/P{Participant_id}_{assembler}.png"),
                      g, width = 6, height = 4, dpi = 300)
  })
}
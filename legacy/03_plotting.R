#!/usr/bin/env Rscript
# =============================================================
# 03_plotting.R (refactored)
# -------------------------------------------------------------
# Reads:
#   results/vf_hits_all.rds
#   results/vf_pa_all.csv
#   results/stats_gene_level.csv
#   results/status_map.csv                [optional]
#   results/diff_genes_UTI_vs_ASB_fisher.csv   [optional]
#
# Writes:
#   results/stats_*csv
#   results/stats_descriptive.md          [if knitr available]
#   results/plots/*.png (and subfolders)
#   results/persistence_*.csv
#   results/permanova_UTI_vs_ASB.txt      [if created earlier]
#
# Adds richer plots:
#   - Volcano (UTI vs ASB enrichment)
#   - Status-stratified prevalence heatmap
#   - PCoA (Jaccard) by status
#   - Gene–gene co-occurrence heatmap
#   - Persistence vs transient bars per participant
#   - UpSet improvements:
#       * Persistence per participant (sets = timepoints)
#       * Status-stratified (sets = UTI/ASB/Negative/None)
#       * Family-level (sets = gene families across participants)
# =============================================================

## ---------------- 0 · setup & helpers ---------------------------------------
suppressPackageStartupMessages({
  library(dplyr);    library(tidyr);    library(readr)
  library(ggplot2);  library(forcats);  library(glue)
  library(purrr);    library(stringr)
})

has_pkg <- function(p) { requireNamespace(p, quietly = TRUE) }
use_pkg <- function(p) { if (has_pkg(p)) suppressPackageStartupMessages(library(p, character.only = TRUE)) }

plots_dir <- "results/plots"
dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)

theme_set(theme_minimal(base_size = 11))
theme_update(legend.title = element_text(hjust = 0),
             legend.text  = element_text(hjust = 0))

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

ensure_tp_lab <- function(df) {
  if ("tp_lab" %in% names(df)) return(df)
  if (!"Timepoint" %in% names(df)) stop("Need 'tp_lab' or 'Timepoint' in data frame.")
  bind_cols(df, tp_norm(df$Timepoint))
}

safe_ggsave <- function(path, plot = last_plot(), width = 7, height = 5, dpi = 300) {
  try(ggsave(path, plot = plot, width = width, height = height, dpi = dpi), silent = TRUE)
}

## ---------------- 1 · load core data ----------------------------------------
vf_hits_path <- "results/vf_hits_all.rds"
vf_pa_path   <- "results/vf_pa_all.csv"
tbl_gene_path<- "results/stats_gene_level.csv"

if (!file.exists(vf_hits_path) || !file.exists(vf_pa_path) || !file.exists(tbl_gene_path)) {
  stop("Missing required inputs. Make sure 02_gene_presence_analysis.R has run.\n",
       "- results/vf_hits_all.rds\n- results/vf_pa_all.csv\n- results/stats_gene_level.csv")
}

vf_hits_all <- readRDS(vf_hits_path)
vf_pa_all   <- read_csv(vf_pa_path, show_col_types = FALSE)
tbl_gene    <- read_csv(tbl_gene_path, show_col_types = FALSE)

vf_pa_all   <- ensure_tp_lab(vf_pa_all)

## ---------------- 2 · sample/participant/cheatsheet -------------------------
# DISTINCT genes per sample (Participant_id, Timepoint)
genes_per_sample <- vf_hits_all %>%
  distinct(Participant_id, Timepoint, GENE) %>%
  count(Participant_id, Timepoint, name = "n_genes") %>%
  bind_cols(tp_norm(.$Timepoint))

participant_tbl <- genes_per_sample %>%
  filter(!is.na(tp_num)) %>%
  arrange(Participant_id, tp_num) %>%
  group_by(Participant_id) %>%
  summarise(
    n_timepoints = n(),
    mean_genes   = mean(n_genes),
    sd_genes     = sd(n_genes),
    min_genes    = min(n_genes),
    max_genes    = max(n_genes),
    genes_T0     = dplyr::first(n_genes),
    genes_last   = dplyr::last(n_genes),
    delta        = genes_last - genes_T0,
    .groups      = "drop"
  )

cohort_tbl <- participant_tbl %>%
  summarise(
    participants       = n(),
    timepoints_total   = sum(n_timepoints),
    genes_median_T0    = median(genes_T0),
    genes_median_last  = median(genes_last),
    mean_delta         = mean(delta),
    sd_delta           = sd(delta)
  )

write_csv(genes_per_sample, "results/stats_sample_level.csv")
write_csv(participant_tbl,  "results/stats_participant_level.csv")
write_csv(cohort_tbl,       "results/stats_cohort_level.csv")
message("✓  CSV summaries written to results/")

# Markdown cheat-sheet (if knitr is available)
if (has_pkg("knitr")) {
  md <- c(
    "# Descriptive statistics\n\n",
    "*Generated on:* ", as.character(Sys.Date()), "\n\n",
    "## Cohort overview\n",
    knitr::kable(cohort_tbl, format = "markdown"), "\n\n",
    "## Participant-level summary\n",
    knitr::kable(participant_tbl, format = "markdown"), "\n\n",
    "## First 10 samples\n",
    knitr::kable(head(genes_per_sample, 10), format = "markdown")
  )
  writeLines(md, "results/stats_descriptive.md")
}

## ---------------- 3 · core/accessory quick plots ----------------------------
dir.create(plots_dir, showWarnings = FALSE)

# Top-25 genes
top25 <- tbl_gene %>% slice_max(n_participants, n = 25) %>%
  mutate(GENE = fct_reorder(GENE, n_participants))
ggplot(top25, aes(GENE, n_participants)) +
  geom_col(fill = "steelblue") + coord_flip() +
  labs(title = "Top 25 VFDB genes (entire cohort)", y = "Participants", x = NULL)
safe_ggsave(file.path(plots_dir, "core_bar_top25_all.png"), width = 6, height = 6)

# Prevalence histogram
ggplot(tbl_gene, aes(n_participants)) +
  geom_histogram(binwidth = 1, fill = "grey70") +
  labs(title = "VF gene prevalence distribution", x = "# participants", y = "Gene count")
safe_ggsave(file.path(plots_dir, "core_histogram_all.png"), width = 5, height = 4)

# Richness by timepoint
ggplot(genes_per_sample, aes(tp_lab, n_genes)) +
  geom_boxplot(outlier.shape = NA, width = 0.6) +
  geom_jitter(width = 0.15, alpha = 0.5, size = 1.5) +
  labs(title = "Per-sample VF gene richness by timepoint", x = "Timepoint", y = "# VF genes (distinct)")
safe_ggsave(file.path(plots_dir, "richness_by_timepoint.png"), width = 7, height = 4.5, dpi = 300)

# Trajectories (numeric)
traj_df <- genes_per_sample %>% filter(!is.na(tp_num))
if (nrow(traj_df) > 0) {
  ggplot(traj_df, aes(tp_num, n_genes, group = Participant_id)) +
    geom_line(alpha = 0.35) +
    geom_point(size = 1.6) +
    scale_x_continuous(breaks = sort(unique(traj_df$tp_num))) +
    labs(title = "Within-participant trajectories of VF gene richness", x = "Numeric timepoint", y = "# VF genes (distinct)")
  safe_ggsave(file.path(plots_dir, "richness_trajectories_numeric.png"), width = 7, height = 4.5, dpi = 300)
} else {
  message("↪  No numeric timepoints found; skipping trajectory plot.")
}

## ---------------- 4 · UpSet improvements ------------------------------------
if (has_pkg("ComplexUpset")) {
  use_pkg("ComplexUpset")
  # 4A) Persistence per participant (sets = timepoints; elements = genes)
  dir.create(file.path(plots_dir, "persistence"), showWarnings = FALSE)
  
  persist_upset_for_pid <- function(pid) {
    df <- vf_hits_all %>%
      filter(Participant_id == pid) %>%
      mutate(tp_lab = tp_norm(Timepoint)$tp_lab) %>%   # compute tp_lab once
      distinct(GENE, tp_lab) %>%
      mutate(val = TRUE) %>%
      pivot_wider(names_from = tp_lab, values_from = val, values_fill = FALSE) %>%
      distinct(GENE, .keep_all = TRUE)
    
    tp_sets <- setdiff(names(df), "GENE")
    if (length(tp_sets) < 2) return(NULL)
    
    ComplexUpset::upset(
      df, intersect = tp_sets, min_size = 1,
      name = paste("P", pid, "genes"),
      base_annotations = list('Intersection size' = intersection_size(text = list(size = 3)))
    )
  }
  
  persist_stats <- purrr::map_dfr(unique(vf_hits_all$Participant_id), function(pid) {
    df <- vf_hits_all %>%
      filter(Participant_id == pid) %>%
      mutate(tp_lab = tp_norm(Timepoint)$tp_lab) %>%
      distinct(GENE, tp_lab) %>%
      mutate(val = TRUE) %>%
      pivot_wider(names_from = tp_lab, values_from = val, values_fill = FALSE) %>%
      distinct(GENE, .keep_all = TRUE)
    
    tp_sets <- setdiff(names(df), "GENE")
    if (length(tp_sets) < 2) return(NULL)
    
    # ✅ across() inside mutate + select(...) is valid; rowSums() on logicals works
    df <- df %>%
      mutate(
        present_n = rowSums(select(., dplyr::all_of(tp_sets))),
        persist   = present_n == length(tp_sets)
      )
    
    readr::write_csv(df %>% arrange(desc(persist), desc(present_n)),
                     glue::glue("results/persistence_P{pid}.csv"))
    
    p <- persist_upset_for_pid(pid)
    if (!is.null(p)) safe_ggsave(glue::glue("{plots_dir}/persistence/upset_P{pid}.png"),
                                 p, width = 10, height = 6, dpi = 300)
    
    tibble::tibble(
      Participant_id   = pid,
      n_timepoints     = length(tp_sets),
      persistent_genes = sum(df$persist),
      transient_genes  = sum(!df$persist)
    )
  })
  
  if (nrow(persist_stats)) readr::write_csv(persist_stats, "results/persistence_summary.csv")
  
  # 4B) Status-stratified UpSet (sets = Infection_Status)
  if (file.exists("results/status_map.csv")) {
    status_map <- read_csv("results/status_map.csv", show_col_types = FALSE) %>%
      ensure_tp_lab() %>%
      select(Participant_id, tp_lab, Infection_Status) %>%
      distinct()
    gene_status_df <- vf_pa_all %>%
      pivot_longer(-c(Participant_id, tp_lab), names_to = "GENE", values_to = "present") %>%
      filter(present > 0) %>% select(-present) %>%
      left_join(status_map, by = c("Participant_id","tp_lab")) %>%
      filter(!is.na(Infection_Status)) %>%
      distinct(GENE, Infection_Status) %>%
      mutate(val = TRUE) %>%
      pivot_wider(names_from = Infection_Status, values_from = val, values_fill = FALSE)
    
    # optional prevalence trimming (5–95%)
    prev <- tbl_gene %>%
      mutate(p = n_participants / n_distinct(vf_pa_all$Participant_id)) %>%
      filter(p >= 0.05, p <= 0.95) %>% pull(GENE)
    gene_status_df <- gene_status_df %>% filter(GENE %in% prev)
    status_sets <- setdiff(names(gene_status_df), "GENE")
    if (length(status_sets) >= 2) {
      p <- ComplexUpset::upset(
        gene_status_df, intersect = status_sets, min_size = 1,
        name = "Genes by status",
        base_annotations = list('Intersection size' = intersection_size(text = list(size = 3)))
      )
      safe_ggsave(file.path(plots_dir, "upset_genes_by_status.png"), p, width = 10, height = 6, dpi = 300)
    }
    
    # Optional: limit to significant UTI vs ASB genes
    if (file.exists("results/diff_genes_UTI_vs_ASB_fisher.csv")) {
      sig <- read_csv("results/diff_genes_UTI_vs_ASB_fisher.csv", show_col_types = FALSE) %>%
        filter(p_adj <= 0.05) %>% pull(GENE)
      sig_df <- gene_status_df %>% filter(GENE %in% sig)
      if (nrow(sig_df)) {
        p <- ComplexUpset::upset(
          sig_df, intersect = status_sets, min_size = 1,
          name = "Sig. genes (UTI vs ASB)",
          base_annotations = list('Intersection size' = intersection_size(text = list(size = 3)))
        )
        safe_ggsave(file.path(plots_dir, "upset_sig_genes_by_status.png"), p, width = 10, height = 6, dpi = 300)
      }
    }
  }
  
  # 4C) Family-level UpSet (sets = gene families; elements = participants)
  get_family <- function(g) {
    m <- stringr::str_match(g, "^([A-Za-z]+)")[,2]
    ifelse(is.na(m), g, m)
  }
  fam_df <- vf_hits_all %>%
    mutate(FAMILY = get_family(GENE)) %>%
    distinct(Participant_id, FAMILY)
  fam_wide <- fam_df %>%
    mutate(val = TRUE) %>%
    pivot_wider(names_from = FAMILY, values_from = val, values_fill = FALSE) %>%
    distinct(Participant_id, .keep_all = TRUE)
  if (nrow(fam_wide)) {
    fam_counts <- colSums(fam_wide[, -1, drop = FALSE])
    n_p <- nrow(fam_wide)
    keep <- names(fam_counts)[fam_counts >= 0.05*n_p & fam_counts <= 0.95*n_p]
    fam_mat <- select(fam_wide, Participant_id, all_of(keep))
    if (ncol(fam_mat) > 2) {
      p <- ComplexUpset::upset(
        fam_mat, intersect = names(fam_mat)[-1], min_size = 1,
        name = "Families across participants",
        base_annotations = list('Intersection size' = intersection_size(text = list(size = 3)))
      )
      safe_ggsave(file.path(plots_dir, "upset_families_by_participant.png"), p, width = 11, height = 6, dpi = 300)
    }
  }
} else {
  message("↪  ComplexUpset not installed; skipping UpSet figures.")
}

## ---------------- 5 · Volcano: UTI vs ASB enrichment -------------------------
if (file.exists("results/diff_genes_UTI_vs_ASB_fisher.csv")) {
  if (has_pkg("ggrepel")) use_pkg("ggrepel")
  de <- read_csv("results/diff_genes_UTI_vs_ASB_fisher.csv", show_col_types = FALSE) %>%
    mutate(
      log2OR = log2(OR),
      neglog10FDR = -log10(p_adj),
      sig = p_adj <= 0.05,
      total_pos = coalesce(`UTI_TRUE`, 0) + coalesce(`ASB_TRUE`, 0)
    )
  top_labs <- bind_rows(
    de %>% arrange(p_adj) %>% slice_head(n = 20),
    de %>% arrange(desc(abs(log2OR))) %>% slice_head(n = 20)
  ) %>% distinct(GENE)
  
  g <- ggplot(de, aes(log2OR, neglog10FDR)) +
    geom_point(aes(alpha = sig, size = total_pos)) +
    geom_vline(xintercept = 0, linetype = "dashed") +
    geom_hline(yintercept = -log10(0.05), linetype = "dotted") +
    { if (has_pkg("ggrepel")) ggrepel::geom_text_repel(data = semi_join(de, top_labs, by = "GENE"),
                                                       aes(label = GENE), max.overlaps = 40, size = 3) else NULL } +
    labs(title = "UTI vs ASB: per-gene enrichment",
         x = "log2(odds ratio)  (UTI / ASB)",
         y = "-log10(FDR)")
  safe_ggsave(file.path(plots_dir, "volcano_UTI_vs_ASB.png"), g, width = 7, height = 5, dpi = 300)
} else {
  message("↪  diff_genes_UTI_vs_ASB_fisher.csv not found; skipping volcano.")
}

## ---------------- 6 · Status-stratified prevalence heatmap -------------------
if (file.exists("results/status_map.csv")) {
  status_map <- read_csv("results/status_map.csv", show_col_types = FALSE)
  status_map <- status_map %>%
    mutate(tp_lab = if (!"tp_lab" %in% names(.)) tp_norm(Timepoint)$tp_lab else tp_lab) %>%
    select(Participant_id, tp_lab, Infection_Status) %>% distinct()
  
  prev_long <- vf_pa_all %>%
    pivot_longer(-c(Participant_id, tp_lab), names_to = "GENE", values_to = "present") %>%
    left_join(status_map, by = c("Participant_id","tp_lab")) %>%
    filter(!is.na(Infection_Status)) %>%
    group_by(GENE, Infection_Status) %>%
    summarise(prev = mean(present > 0) * 100, .groups = "drop")
  
  rank40 <- prev_long %>%
    pivot_wider(names_from = Infection_Status, values_from = prev) %>%
    mutate(delta_UTI_ASB = coalesce(UTI, 0) - coalesce(ASB, 0)) %>%
    arrange(desc(abs(delta_UTI_ASB))) %>%
    slice_head(n = 40) %>% pull(GENE)
  
  heat_df <- prev_long %>%
    filter(GENE %in% rank40) %>%
    mutate(
      Infection_Status = factor(Infection_Status, levels = c("UTI","ASB","Negative","None")),
      GENE = factor(GENE, levels = rev(rank40))
    )
  
  g <- ggplot(heat_df, aes(Infection_Status, GENE, fill = prev)) +
    geom_tile() +
    scale_fill_gradient(name = "% present", low = "white", high = "steelblue") +
    labs(title = "Gene prevalence by clinical status (top 40 by |UTI–ASB|)") +
    theme(axis.text.y = element_text(size = 6))
  safe_ggsave(file.path(plots_dir, "heatmap_prevalence_by_status_top40.png"), g, width = 7.5, height = 10, dpi = 300)
} else {
  message("↪  status_map.csv not found; skipping status heatmap.")
}

## ---------------- 7 · PCoA (Jaccard) by status -------------------------------
if (has_pkg("vegan")) {
  use_pkg("vegan"); if (has_pkg("ggrepel")) use_pkg("ggrepel")
  if (file.exists("results/status_map.csv")) {
    status_map <- read_csv("results/status_map.csv", show_col_types = FALSE) %>%
      mutate(tp_lab = if (!"tp_lab" %in% names(.)) tp_norm(Timepoint)$tp_lab else tp_lab) %>%
      select(Participant_id, tp_lab, Infection_Status) %>% distinct()
    ord_data <- vf_pa_all %>% left_join(status_map, by = c("Participant_id","tp_lab"))
  } else {
    ord_data <- vf_pa_all %>% mutate(Infection_Status = NA_character_)
  }
  X <- ord_data %>% select(-Participant_id, -tp_lab, -Infection_Status) %>% as.matrix()
  X <- X > 0
  if (nrow(X) >= 3) {
    dJ <- vegan::vegdist(X, method = "jaccard", binary = TRUE)
    pcoa <- stats::cmdscale(dJ, k = 2, eig = TRUE)
    ord <- tibble::tibble(PC1 = pcoa$points[,1], PC2 = pcoa$points[,2]) %>%
      bind_cols(ord_data %>% select(Participant_id, tp_lab, Infection_Status))
    g <- ggplot(ord, aes(PC1, PC2, color = Infection_Status)) +
      geom_point(alpha = 0.85, size = 2) +
      labs(title = "PCoA (Jaccard) of VF presence/absence",
           subtitle = "Colored by Infection_Status (if available)")
    safe_ggsave(file.path(plots_dir, "pcoa_jaccard_status.png"), g, width = 7, height = 5, dpi = 300)
  } else message("↪  Too few samples for PCoA.")
} else {
  message("↪  vegan not installed; skipping PCoA.")
}

## ---------------- 8 · Gene–gene co-occurrence heatmap ------------------------
# choose variable genes (10–90% prev), top 60 by variance proxy p*(1-p)
prev_tbl <- vf_pa_all %>%
  select(-Participant_id, -tp_lab) %>%
  summarise(across(everything(), ~ mean(. > 0))) %>%
  pivot_longer(everything(), names_to = "GENE", values_to = "prev")

sel_genes <- prev_tbl %>%
  filter(prev >= 0.10, prev <= 0.90) %>%
  mutate(score = prev * (1 - prev)) %>%
  arrange(desc(score)) %>%
  slice_head(n = 60) %>% pull(GENE)

if (length(sel_genes) >= 4) {
  B <- vf_pa_all %>% select(all_of(sel_genes))
  # Prefer Jaccard via vegan if available; else fallback to binary distance
  if (has_pkg("vegan")) {
    use_pkg("vegan")
    D <- vegan::vegdist(t(B > 0), method = "jaccard", binary = TRUE)
    Sim <- as.matrix(1 - D)
  } else {
    D <- stats::dist(t(B > 0), method = "binary")
    Sim <- 1 - as.matrix(D)  # not true Jaccard, but acceptable fallback
  }
  hc  <- hclust(as.dist(1 - Sim), method = "average")
  ord <- hc$labels[hc$order]
  mat <- Sim[ord, ord]
  df  <- as.data.frame(as.table(mat), stringsAsFactors = FALSE)
  names(df) <- c("Gene1","Gene2","Jaccard")
  
  g <- ggplot(df, aes(Gene1, Gene2, fill = Jaccard)) +
    geom_tile() +
    scale_fill_gradient(low = "white", high = "steelblue") +
    labs(title = "Gene–gene co-occurrence (Jaccard similarity)") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 6),
          axis.text.y = element_text(size = 6))
  safe_ggsave(file.path(plots_dir, "heatmap_gene_cooccurrence_top60.png"), g, width = 8, height = 8, dpi = 300)
} else {
  message("↪  Not enough variable genes for co-occurrence heatmap.")
}

## ---------------- 9 · Persistence vs transient summary -----------------------
persist_tbl <- vf_hits_all %>%
  distinct(Participant_id, GENE, tp_lab = tp_norm(Timepoint)$tp_lab) %>%
  group_by(Participant_id, GENE) %>%
  summarise(n_tp = n_distinct(tp_lab), .groups = "drop")

tp_counts <- vf_hits_all %>%
  distinct(Participant_id, tp_lab = tp_norm(Timepoint)$tp_lab) %>%
  count(Participant_id, name = "n_tp") %>% filter(n_tp >= 2)

stab <- persist_tbl %>%
  semi_join(tp_counts, by = "Participant_id") %>%
  group_by(Participant_id) %>%
  summarise(
    genes_total      = n(),
    genes_persistent = sum(n_tp == max(n_tp)),
    genes_transient  = sum(n_tp <  max(n_tp)),
    frac_persistent  = round(100 * genes_persistent / genes_total, 1),
    .groups = "drop"
  ) %>% arrange(desc(frac_persistent))

write_csv(stab, "results/persistence_per_participant.csv")

g <- ggplot(stab, aes(reorder(Participant_id, frac_persistent), frac_persistent)) +
  geom_col() + coord_flip() +
  labs(title = "Share of persistent genes per participant", x = "Participant", y = "% persistent genes")
safe_ggsave(file.path(plots_dir, "bar_persistence_share_by_participant.png"), g, width = 7, height = 8, dpi = 300)

## ---------------- 10 · Done --------------------------------------------------
message("✓  All plots complete – check results/plots/ and CSVs in results/.")
#!/usr/bin/env Rscript
# ==============================================================================
# 03_plotting.R
# ==============================================================================
#
# GOAL:
#   Generate comprehensive VF-focused visualisations: gene prevalence plots,
#   heatmaps, and initial burden comparisons.  This was the original "plotting
#   engine" for VF data.  Many of its figures are now superseded by the more
#   rigorous outputs from 23_vf_cross_sectional.R and 05_gene_overview_plots.R,
#   but it remains useful for exploratory visualisation.
#
# ------------------------------------------------------------------------------
# Purpose: Generates comprehensive visualizations for gene presence/absence,
#          phylogeny, epidemiology, and transmission networks.
# Inputs:
#   - 00_config.R
#   - results/vf_hits_all.rds
#   - results/vf_pa_all.csv
#   - results/stats_gene_level.csv
#   - results/status_map.csv (optional)
#   - results/mlst/mlst_all.tsv
#   - results/wgs/parsnp/parsnp.tree (optional)
#   - results/nitrate_presence_matrix.csv (optional)
#   - results/pairwise_stats.csv (optional)
#   - data/inputs/batch*.csv (for robust metadata mapping)
# Outputs:
#   - results/stats_*.csv
#   - results/plots/ (various PNGs)
# ==============================================================================

# 1. Load Configuration & Libraries
source("00_config.R")
source("R/plot_helpers.R")
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(ggplot2)
  library(forcats)
  library(glue)
  library(purrr)
  library(stringr)
  library(tibble)
})

# Standard Color Palette
# Standard Color Palette
# Note: Infection colors are now handled by scale_colour_infection() in R/plot_helpers.R
rutis_palette <- c(
  `Culture-positive` = "#009E73",
  Other = "#999999",
  `Within-Host` = "#0072B2",
  `Between-Host` = "#CC79A7"
)

# 2. Helper Functions
# ------------------------------------------------------------------------------
has_pkg <- function(p) {
  requireNamespace(p, quietly = TRUE)
}
use_pkg <- function(p) {
  if (has_pkg(p)) suppressPackageStartupMessages(library(p, character.only = TRUE))
}

safe_ggsave <- function(filename, plot = last_plot(), width = 7, height = 5, dpi = 300) {
  path <- file.path(DIR_PLOTS, filename)
  dir_name <- dirname(path)
  if (!dir.exists(dir_name)) dir.create(dir_name, recursive = TRUE, showWarnings = FALSE)

  tryCatch(
    ggsave(path, plot = plot, width = width, height = height, dpi = dpi),
    error = function(e) warning("Failed to save plot: ", path, "\n", e$message)
  )
}

# Timepoint Normalization
tp_norm <- function(x) {
  tp_chr <- as.character(x)
  is_uricult <- str_detect(tp_chr, regex("uricult", ignore_case = TRUE))
  tp_num <- suppressWarnings(as.integer(str_extract(tp_chr, "\\d+")))

  tp_lab <- dplyr::case_when(
    is_uricult ~ "Uricult",
    !is.na(tp_num) ~ paste0("T", tp_num),
    TRUE ~ "Unscheduled"
  )

  tp_levels <- c(
    paste0("T", sort(unique(tp_num[!is.na(tp_num)]))),
    "Uricult", "Unscheduled"
  )

  tibble::tibble(
    tp_lab = factor(tp_lab, levels = tp_levels),
    tp_num = tp_num
  )
}

ensure_tp_lab <- function(df) {
  if ("tp_lab" %in% names(df)) {
    return(df)
  }
  if (!"Timepoint" %in% names(df)) stop("Need 'tp_lab' or 'Timepoint' in data frame.")
  bind_cols(df, tp_norm(df$Timepoint))
}

# Robust Metadata Parser
# Loads batch*.csv files to map isolate_ID (e.g. 2446C038801-1) to Participant_id/Timepoint
load_robust_metadata <- function() {
  # Helper to safely read and cast types
  read_batch_safe <- function(path) {
    if (!file.exists(path)) {
      return(NULL)
    }
    df <- read_csv(path, show_col_types = FALSE)

    # Fix BOM if present
    nms <- names(df)
    nms[1] <- gsub("^[^A-Za-z0-9]+", "", nms[1])
    names(df) <- nms

    # Normalize Timepoint column name
    if (!"Timepoint" %in% names(df)) {
      tp_col <- names(df)[grepl("time.*point", names(df), ignore.case = TRUE)]
      if (length(tp_col) > 0) df <- df %>% rename(Timepoint = all_of(tp_col))
    }

    if (all(c("Participant_id", "isolate_ID", "Timepoint") %in% names(df))) {
      df %>%
        mutate(
          Participant_id = as.character(Participant_id),
          isolate_ID = as.character(isolate_ID),
          Timepoint = as.character(Timepoint)
        ) %>%
        select(Participant_id, isolate_ID, Timepoint)
    } else {
      NULL
    }
  }

  b1 <- read_batch_safe("data/inputs/batch1.csv")
  b2 <- read_batch_safe("data/inputs/batch2.csv")
  b3 <- read_batch_safe("data/inputs/batch3.csv")

  bind_rows(b1, b2, b3) %>%
    mutate(Timepoint = tp_norm(Timepoint)$tp_lab)
}

# 3. Load Data
# ------------------------------------------------------------------------------
FILE_GENE_STATS <- file.path(DIR_VF, "stats_gene_level.csv")
FILE_STATUS_MAP <- file.path(DIR_RESULTS, "clinical", "status_map.csv")

if (!file.exists(FILE_VF_HITS) || !file.exists(FILE_VF_PA) || !file.exists(FILE_GENE_STATS)) {
  stop("Missing required inputs. Run 02_gene_presence_analysis.R first.")
}

vf_hits_all <- readRDS(FILE_VF_HITS)
vf_pa_all <- read_csv(FILE_VF_PA, show_col_types = FALSE)
tbl_gene <- read_csv(FILE_GENE_STATS, show_col_types = FALSE)

vf_pa_all <- ensure_tp_lab(vf_pa_all)

# Load Robust Metadata Map
meta_map <- load_robust_metadata()

# 4. Generate Statistics
# ------------------------------------------------------------------------------
# Per-sample richness
genes_per_sample <- vf_hits_all %>%
  distinct(Participant_id, Timepoint, GENE) %>%
  count(Participant_id, Timepoint, name = "n_genes") %>%
  bind_cols(tp_norm(.$Timepoint))

# Per-participant summary
participant_tbl <- genes_per_sample %>%
  filter(!is.na(tp_num) | tp_lab == "Uricult") %>%
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

# Cohort summary
cohort_tbl <- participant_tbl %>%
  summarise(
    participants       = n(),
    timepoints_total   = sum(n_timepoints),
    genes_median_T0    = median(genes_T0),
    genes_median_last  = median(genes_last),
    mean_delta         = mean(delta),
    sd_delta           = sd(delta)
  )

# Write stats
write_csv(genes_per_sample, file.path(DIR_RESULTS, "stats_sample_level.csv"))
write_csv(participant_tbl, file.path(DIR_RESULTS, "stats_participant_level.csv"))
write_csv(cohort_tbl, file.path(DIR_RESULTS, "stats_cohort_level.csv"))

# Markdown report
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
  writeLines(md, file.path(DIR_RESULTS, "stats_descriptive.md"))
}

message("✓ Statistics generated.")

# 5. Core Plots
# ------------------------------------------------------------------------------
ensure_dir(DIR_PLOTS)

# Top 25 Genes
top25 <- tbl_gene %>%
  slice_max(n_participants, n = 25) %>%
  mutate(GENE = fct_reorder(GENE, n_participants))

ggplot(top25, aes(GENE, n_participants)) +
  geom_col(fill = "steelblue") +
  coord_flip() +
  labs(title = "Top 25 Most Prevalent Virulence Factor Genes", y = "Participants", x = NULL) +
  theme_minimal()
safe_ggsave("core_bar_top25_all.png", width = 6, height = 6)

# Prevalence Histogram
ggplot(tbl_gene, aes(n_participants)) +
  geom_histogram(binwidth = 1, fill = "grey70") +
  labs(title = "Distribution of Virulence Gene Prevalence", x = "Number of Participants", y = "Count of Genes") +
  theme_minimal()
safe_ggsave("core_histogram_all.png", width = 5, height = 4)

# Richness by Timepoint
ggplot(genes_per_sample, aes(tp_lab, n_genes)) +
  geom_boxplot(outlier.shape = NA, width = 0.6) +
  geom_jitter(width = 0.15, alpha = 0.5, size = 1.5) +
  labs(title = "Virulence Gene Richness per Sample by Timepoint", x = "Timepoint", y = "Number of Detected VF Genes") +
  theme_minimal()
safe_ggsave("richness_by_timepoint.png", width = 7, height = 4.5)

# Trajectories
traj_df <- genes_per_sample %>% filter(!is.na(tp_num))
if (nrow(traj_df) > 0) {
  ggplot(traj_df, aes(tp_num, n_genes, group = Participant_id)) +
    geom_line(alpha = 0.35) +
    geom_point(size = 1.6) +
    scale_x_continuous(breaks = sort(unique(traj_df$tp_num))) +
    labs(title = "Longitudinal Changes in Virulence Gene Richness", x = "Timepoint", y = "Number of Detected VF Genes") +
    theme_minimal()
  safe_ggsave("richness_trajectories_numeric.png", width = 7, height = 4.5)
}

# 6. UpSet Plots (if ComplexUpset available)
# ------------------------------------------------------------------------------
if (has_pkg("ComplexUpset")) {
  use_pkg("ComplexUpset")
  ensure_dir(file.path(DIR_PLOTS, "persistence"))

  # 6A. Persistence per Participant
  persist_stats <- purrr::map_dfr(unique(vf_hits_all$Participant_id), function(pid) {
    df <- vf_hits_all %>%
      filter(Participant_id == pid) %>%
      mutate(tp_lab = tp_norm(Timepoint)$tp_lab) %>%
      distinct(GENE, tp_lab) %>%
      mutate(val = TRUE) %>%
      pivot_wider(names_from = tp_lab, values_from = val, values_fill = FALSE) %>%
      distinct(GENE, .keep_all = TRUE)

    tp_sets <- setdiff(names(df), "GENE")
    if (length(tp_sets) < 2) {
      return(NULL)
    }

    # Calculate persistence
    df <- df %>%
      mutate(
        present_n = rowSums(select(., all_of(tp_sets))),
        persist   = present_n == length(tp_sets)
      )

    write_csv(
      df %>% arrange(desc(persist), desc(present_n)),
      file.path(DIR_RESULTS, glue::glue("persistence_P{pid}.csv"))
    )

    # Plot
    try(
      {
        p <- ComplexUpset::upset(
          df,
          intersect = tp_sets, min_size = 1,
          name = paste("P", pid, "genes"),
          base_annotations = list("Intersection size" = intersection_size(text = list(size = 3)))
        )
        safe_ggsave(file.path("persistence", glue::glue("upset_P{pid}.png")), p, width = 10, height = 6)
      },
      silent = TRUE
    )

    tibble::tibble(
      Participant_id   = pid,
      n_timepoints     = length(tp_sets),
      persistent_genes = sum(df$persist),
      transient_genes  = sum(!df$persist)
    )
  })

  if (nrow(persist_stats)) write_csv(persist_stats, file.path(DIR_RESULTS, "persistence_summary.csv"))

  # 6B. Status-stratified UpSet
  if (file.exists(FILE_STATUS_MAP)) {
    status_map <- read_csv(FILE_STATUS_MAP, show_col_types = FALSE) %>%
      ensure_tp_lab() %>%
      select(Participant_id, tp_lab, Infection_Status) %>%
      distinct()

    gene_status_df <- vf_pa_all %>%
      pivot_longer(-c(Participant_id, tp_lab), names_to = "GENE", values_to = "present") %>%
      filter(present > 0) %>%
      select(-present) %>%
      left_join(status_map, by = c("Participant_id", "tp_lab"), relationship = "many-to-many") %>%
      filter(!is.na(Infection_Status)) %>%
      distinct(GENE, Infection_Status) %>%
      mutate(val = TRUE) %>%
      pivot_wider(names_from = Infection_Status, values_from = val, values_fill = FALSE)

    # Filter by prevalence (5-95%)
    prev <- tbl_gene %>%
      mutate(p = n_participants / n_distinct(vf_pa_all$Participant_id)) %>%
      filter(p >= 0.05, p <= 0.95) %>%
      pull(GENE)

    gene_status_df <- gene_status_df %>% filter(GENE %in% prev)
    status_sets <- setdiff(names(gene_status_df), "GENE")

    if (length(status_sets) >= 2) {
      p <- ComplexUpset::upset(
        gene_status_df,
        intersect = status_sets, min_size = 1,
        name = "Intersection of Virulence Genes by Infection Status"
      )
      safe_ggsave("upset_genes_by_status.png", p, width = 10, height = 6)
    }
  }
}

# 7. Volcano Plot (UTI vs ASB)
# ------------------------------------------------------------------------------
FILE_FISHER <- file.path(DIR_RESULTS, "diff_genes_UTI_vs_ASB_fisher.csv")
if (file.exists(FILE_FISHER)) {
  if (has_pkg("ggrepel")) use_pkg("ggrepel")
  de <- read_csv(FILE_FISHER, show_col_types = FALSE) %>%
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
    {
      if (has_pkg("ggrepel")) {
        ggrepel::geom_text_repel(
          data = semi_join(de, top_labs, by = "GENE"),
          aes(label = GENE), max.overlaps = 40, size = 3
        )
      } else {
        NULL
      }
    } +
    labs(
      title = "Genomic Features Associated with UTI vs. ASB",
      x = "log2(Odds Ratio)",
      y = "-log10(FDR)"
    ) +
    theme_minimal()

  safe_ggsave("volcano_UTI_vs_ASB.png", g, width = 7, height = 5)
}

# ==============================================================================
# SECTION 1.1: Core-genome Phylogeny
# ==============================================================================
FILE_TREE <- file.path(DIR_WGS, "parsnp", "parsnp.tree")
if (file.exists(FILE_TREE) && has_pkg("ape")) {
  ensure_dir(file.path(DIR_PLOTS, "phylogeny"))
  use_pkg("ape")

  tree <- tryCatch(read.tree(FILE_TREE), error = function(e) NULL)

  if (!is.null(tree)) {
    if (file.exists(FILE_STATUS_MAP) && !is.null(meta_map)) {
      status <- read_csv(FILE_STATUS_MAP, show_col_types = FALSE) %>% ensure_tp_lab()

      # Map status to isolate_ID for tree annotation
      # Tree tips are likely isolate_IDs (e.g. 24...) or filenames
      # We need a map: isolate_ID -> Infection_Status

      annot <- meta_map %>%
        inner_join(status, by = c("Participant_id", "Timepoint" = "tp_lab"), relationship = "many-to-many") %>%
        select(isolate_ID, Infection_Status) %>%
        distinct()

      # Map tree tips to isolate_ID
      # Tree tips might be full filenames or isolate IDs.
      # We attempt to match tip labels to isolate_ID.
      tip_labels <- tree$tip.label

      # Create a mapping table
      tip_map <- tibble(label = tip_labels) %>%
        mutate(
          # Try to extract ID from label if it looks like a filename
          extracted_id = str_extract(label, "24[0-9A-Za-z]+-[0-9]+"),
          # Fallback: if extraction fails, assume label IS the ID
          match_id = coalesce(extracted_id, label)
        ) %>%
        left_join(annot, by = c("match_id" = "isolate_ID"), relationship = "many-to-many") %>%
        select(label, Infection_Status)

      # Update annot to use the exact tip labels from the tree
      annot <- tip_map %>%
        filter(!is.na(Infection_Status)) %>%
        group_by(label) %>%
        slice(1) %>%
        ungroup()

      if (has_pkg("ggtree")) {
        use_pkg("ggtree")
        try({
          p <- ggtree(tree, layout = "rectangular") %<+% annot +
            geom_tippoint(aes(color = Infection_Status), size = 2, alpha = 0.8) +
            scale_colour_infection() +
            theme_tree2() +
            labs(title = "Core Genome Phylogeny of E. coli Isolates", color = "Infection Status")
          safe_ggsave("phylogeny/core_tree_phenotype.png", p, width = 8, height = 10)
        })
      } else {
        png(file.path(DIR_PLOTS, "phylogeny", "core_tree_base.png"), width = 800, height = 1000)
        plot(tree, show.tip.label = FALSE, main = "Core-genome Phylogeny")
        dev.off()
      }
    }
  }
}

# ==============================================================================
# SECTION 1.2: ST Distribution (ASB vs UTI)
# ==============================================================================
if (file.exists(FILE_MLST_ALL) && file.exists(FILE_STATUS_MAP) && !is.null(meta_map)) {
  ensure_dir(file.path(DIR_PLOTS, "epidemiology"))

  mlst <- read_tsv(FILE_MLST_ALL, show_col_types = FALSE)
  status <- read_csv(FILE_STATUS_MAP, show_col_types = FALSE) %>% ensure_tp_lab()

  # Extract Isolate_ID from filename
  mlst_parsed <- mlst %>%
    mutate(extracted_ID = str_extract(file_name, "24[0-9A-Za-z]+-[0-9]+")) %>%
    select(-any_of(c("Participant_id", "Timepoint"))) # Drop conflicting columns

  # Join MLST -> Meta Map -> Status
  st_df <- mlst_parsed %>%
    inner_join(meta_map, by = c("extracted_ID" = "isolate_ID")) %>%
    inner_join(status, by = c("Participant_id", "Timepoint" = "tp_lab")) %>%
    filter(!is.na(ST), !is.na(Infection_Status)) %>%
    filter(Infection_Status %in% c("ASB", "UTI"))

  if (nrow(st_df) > 0) {
    top_sts <- st_df %>%
      count(ST, sort = TRUE) %>%
      slice_head(n = 10) %>%
      pull(ST)
    st_plot_df <- st_df %>%
      mutate(ST_grouped = ifelse(ST %in% top_sts, as.character(ST), "Other"))

    g <- ggplot(st_plot_df, aes(x = Infection_Status, fill = ST_grouped)) +
      geom_bar(position = "fill") +
      scale_y_continuous(labels = scales::percent) +
      scale_fill_brewer(palette = "Set3", name = "Sequence Type") +
      labs(title = "Sequence Type Distribution by Infection Status", x = "Infection Status", y = "Proportion of Isolates") +
      theme_minimal(base_size = 20) +
      theme(
        plot.title = element_text(size = 26, face = "bold", hjust = 0.5),
        legend.position = "top",
        legend.text = element_text(size = 16),
        axis.text.x = element_text(size = 18, face = "bold"),
        axis.text.y = element_text(size = 18, face = "bold"),
        plot.margin = margin(15, 15, 15, 15)
      )
    safe_ggsave("epidemiology/st_distribution_stacked.png", g, width = 12, height = 8, dpi = 600)
  }
}

# ==============================================================================
# SECTION 1.3: Virulence Burden (ASB vs UTI)
# ==============================================================================
if (exists("vf_hits_all") && file.exists(FILE_STATUS_MAP)) {
  status <- read_csv(FILE_STATUS_MAP, show_col_types = FALSE) %>% ensure_tp_lab()

  vf_counts <- vf_hits_all %>%
    count(Participant_id, Timepoint, name = "n_vf") %>%
    # Ensure Timepoint is normalized in vf_hits_all before join
    mutate(tp_lab = tp_norm(Timepoint)$tp_lab) %>%
    left_join(status, by = c("Participant_id", "tp_lab")) %>%
    filter(Infection_Status %in% c("ASB", "UTI"))

  if (nrow(vf_counts) > 0) {
    g <- ggplot(vf_counts, aes(x = Infection_Status, y = n_vf, fill = Infection_Status)) +
      geom_boxplot(outlier.shape = NA, alpha = 0.7) +
      geom_jitter(width = 0.2, alpha = 0.4) +
      scale_fill_manual(values = rutis_palette, name = "Infection Status") +
      labs(title = "Total Virulence Gene Burden by Infection Status", y = "Total VF Genes", x = "Infection Status") +
      theme_minimal() +
      theme(legend.position = "none")

    safe_ggsave("epidemiology/vf_burden_boxplot.png", g, width = 5, height = 5)

    w_test <- wilcox.test(n_vf ~ Infection_Status, data = vf_counts)
    capture.output(w_test, file = file.path(DIR_RESULTS, "stats_vf_burden_wilcox.txt"))
  }
}

# ==============================================================================
# SECTION 1.4: Pairwise SNP Distances
# ==============================================================================
FILE_PAIR_STATS <- file.path(DIR_RESULTS, "pairwise_stats.csv")
if (file.exists(FILE_PAIR_STATS)) {
  ensure_dir(file.path(DIR_PLOTS, "genomics"))
  pairs <- read_csv(FILE_PAIR_STATS, show_col_types = FALSE)

  # Check for TotalSnpCnt (or TotalSNPs)
  snp_col <- if ("TotalSnpCnt" %in% names(pairs)) "TotalSnpCnt" else "TotalSNPs"

  if (snp_col %in% names(pairs)) {
    # Normalize column name
    pairs$TotalSNPs <- pairs[[snp_col]]

    # Check for within_participant (might be missing if not generated by 12_wgs)
    if (!"within_participant" %in% names(pairs)) {
      # Infer from Participant_id if available, or skip category
      pairs$within_participant <- FALSE # Placeholder if missing
    }

    plot_pairs <- pairs %>%
      filter(!is.na(TotalSNPs)) %>%
      mutate(Category = ifelse(within_participant, "Within-Host", "Between-Host"))

    g <- ggplot(plot_pairs, aes(x = Category, y = TotalSNPs + 1, fill = Category)) +
      geom_violin(alpha = 0.5) +
      geom_jitter(height = 0, width = 0.2, alpha = 0.1, size = 0.5) +
      scale_y_log10() +
      scale_fill_manual(values = rutis_palette, name = "Comparison Type") +
      labs(title = "Pairwise SNP Distances: Within-Host vs. Between-Host", y = "SNP Distance (log10 + 1)", x = "Comparison Type") +
      theme_minimal() +
      theme(legend.position = "none")

    safe_ggsave("genomics/snp_distance_violin.png", g, width = 6, height = 5)
  }
}

# ==============================================================================
# SECTION 2.1: Longitudinal Timelines
# ==============================================================================
if ((file.exists(FILE_STATUS_MAP) || file.exists(file.path(DIR_CLINICAL, "status_map_with_poster_tp.csv"))) && file.exists(FILE_MLST_ALL) && !is.null(meta_map)) {
  ensure_dir(file.path(DIR_PLOTS, "timelines"))

  status_file <- file.path(DIR_CLINICAL, "status_map_with_poster_tp.csv")
  if (!file.exists(status_file)) status_file <- FILE_STATUS_MAP
  status <- read_csv(status_file, show_col_types = FALSE) %>% ensure_tp_lab()
  mlst <- read_tsv(FILE_MLST_ALL, show_col_types = FALSE) %>%
    mutate(extracted_ID = str_extract(file_name, "24[0-9A-Za-z]+-[0-9]+")) %>%
    select(-any_of(c("Participant_id", "Timepoint")))

  timeline_df <- mlst %>%
    inner_join(meta_map, by = c("extracted_ID" = "isolate_ID")) %>%
    inner_join(status, by = c("Participant_id", "Timepoint" = "tp_lab")) %>%
    mutate(
      ST_Label = ifelse(is.na(ST), "Unknown", paste0("ST", ST)),
      tp_num = if("Plot_TP_Num_Poster" %in% names(.)) Plot_TP_Num_Poster else tp_norm(Timepoint)$tp_num,
      Plot_Label = if("Plot_TP_Label_Poster" %in% names(.)) Plot_TP_Label_Poster else Timepoint
    ) %>%
    filter(!is.na(tp_num)) %>%
    arrange(Participant_id, tp_num)

  if (nrow(timeline_df) > 0) {
    top_pids <- timeline_df %>%
      count(Participant_id, sort = TRUE) %>%
      slice_head(n = 20) %>%
      pull(Participant_id)

    plot_data <- timeline_df %>% filter(Participant_id %in% top_pids) %>% 
      mutate(Plot_Label = fct_reorder(Plot_Label, tp_num))

    g <- ggplot(
      plot_data,
      aes(x = Plot_Label, y = as.factor(Participant_id))
    ) +
      geom_line(color = "grey80") +
      geom_point(aes(color = Infection_Status, shape = ST_Label), size = 3) +
      scale_color_manual(values = rutis_palette) +
      labs(
        title = "Longitudinal Infection Dynamics (Top 20 Participants)",
        x = "Timepoint", y = "Participant",
        color = "Infection Status", shape = "ST"
      ) +
      theme_minimal()

    safe_ggsave("timelines/swimmer_plot_top20.png", g, width = 10, height = 8)
  }
}

# ==============================================================================
# SECTION 2.2: Network / "Constellation" Plot
# ==============================================================================
if (file.exists(FILE_PAIR_STATS) && has_pkg("igraph") && has_pkg("ggraph")) {
  use_pkg("igraph")
  use_pkg("ggraph")

  pairs <- read_csv(FILE_PAIR_STATS, show_col_types = FALSE)
  SNP_THRESHOLD <- 10

  # Check for TotalSnpCnt (or TotalSNPs)
  snp_col <- if ("TotalSnpCnt" %in% names(pairs)) "TotalSnpCnt" else "TotalSNPs"

  if (snp_col %in% names(pairs)) {
    pairs$TotalSNPs <- pairs[[snp_col]]

    # Construct edges. If SampleA/SampleB missing, extract from path_A/path_B
    if (!"SampleA" %in% names(pairs) && "path_A" %in% names(pairs)) {
      pairs <- pairs %>%
        mutate(
          SampleA = str_extract(basename(path_A), "24[0-9A-Za-z]+-[0-9]+"),
          SampleB = str_extract(basename(path_B), "24[0-9A-Za-z]+-[0-9]+")
        )
      # Fallback if regex fails (e.g. different naming convention): use full basename
      pairs$SampleA[is.na(pairs$SampleA)] <- basename(pairs$path_A[is.na(pairs$SampleA)])
      pairs$SampleB[is.na(pairs$SampleB)] <- basename(pairs$path_B[is.na(pairs$SampleB)])
    }

    edges <- pairs %>%
      filter(TotalSNPs <= SNP_THRESHOLD) %>%
      mutate(weight = (SNP_THRESHOLD + 1) - TotalSNPs) %>%
      select(from = SampleA, to = SampleB, weight)

    if (nrow(edges) > 0) {
      try({
        g_net <- graph_from_data_frame(edges, directed = FALSE)

        if (file.exists(FILE_STATUS_MAP) && !is.null(meta_map)) {
          status <- read_csv(FILE_STATUS_MAP, show_col_types = FALSE) %>% ensure_tp_lab()

          # Map SampleID (isolate_ID) -> Infection_Status
          annot <- meta_map %>%
            inner_join(status, by = c("Participant_id", "Timepoint" = "tp_lab")) %>%
            select(isolate_ID, Infection_Status) %>%
            distinct()

          # Match IDs. The edges use basenames or whatever was in path_A.
          # We might need to clean IDs.
          # For now, try direct match.
          V(g_net)$Phenotype <- annot$Infection_Status[match(V(g_net)$name, annot$isolate_ID)]
        }

        p <- ggraph(g_net, layout = "fr") +
          geom_edge_link(aes(alpha = 0.5), show.legend = FALSE) +
          geom_node_point(aes(color = Phenotype), size = 3) +
          scale_color_manual(values = rutis_palette, na.value = "grey80") +
          theme_graph() +
          labs(title = paste("Transmission Network (SNPs <=", SNP_THRESHOLD, ")"))

        safe_ggsave("epidemiology/transmission_network.png", p, width = 8, height = 8)
      })
    }
  }
}

# ==============================================================================
# SECTION 2.3: Integrated Virulence + Nitrate/Fitness Heatmap
# ==============================================================================
FILE_NITRATE <- file.path(DIR_RESULTS, "nitrate_presence_matrix.csv")
if (file.exists(FILE_NITRATE) && exists("vf_pa_all") && has_pkg("pheatmap")) {
  use_pkg("pheatmap")

  tryCatch(
    {
      nitrate <- read_csv(FILE_NITRATE, show_col_types = FALSE) %>%
        mutate(Participant_id = as.character(Participant_id))

      top_vfs <- top25$GENE[1:min(20, nrow(top25))]

      vf_subset <- vf_pa_all %>%
        mutate(Participant_id = as.character(Participant_id)) %>%
        select(Participant_id, tp_lab, all_of(top_vfs))

      nitrate_subset <- nitrate %>%
        select(Participant_id, Timepoint, any_of(c("Nar", "Nap", "Nas"))) %>%
        mutate(tp_lab = tp_norm(Timepoint)$tp_lab) %>%
        select(-Timepoint)

      combined <- vf_subset %>%
        inner_join(nitrate_subset, by = c("Participant_id", "tp_lab"))

      if (nrow(combined) > 0) {
        mat <- combined %>%
          select(-Participant_id, -tp_lab) %>%
          as.matrix()
        rownames(mat) <- paste(combined$Participant_id, combined$tp_lab, sep = "_")

        if (file.exists(FILE_STATUS_MAP)) {
          status_map <- read_csv(FILE_STATUS_MAP, show_col_types = FALSE) %>% ensure_tp_lab()
          annot_df <- status_map %>%
            mutate(ID = paste(Participant_id, tp_lab, sep = "_")) %>%
            filter(ID %in% rownames(mat)) %>%
            select(ID, Infection_Status) %>%
            column_to_rownames("ID")

          png(file.path(DIR_PLOTS, "genomics", "virulence_nitrate_heatmap.png"), width = 1000, height = 1200, res = 150)
          pheatmap(mat,
            annotation_row = annot_df,
            main = "Virulence & Nitrate Gene Profiles",
            color = c("white", "navy"), legend_breaks = c(0, 1), legend_labels = c("Absent", "Present")
          )
          dev.off()
        }
      }
    },
    error = function(e) {
      message("Skipping heatmap due to error: ", e$message)
    }
  )
}

message("✓ All plots generated.")

#!/usr/bin/env Rscript
# ==============================================================================
# 09_inc_plasmid_network.R
# ==============================================================================
#
# GOAL:
#   Detect plasmid replicons (Inc types) using PlasmidFinder, build replicon
#   co-occurrence matrices, and visualise ST–plasmid associations as networks.
#   Produces the plasmidfinder_presence_absence.csv used by scripts 10, 11, 14.
#
# ------------------------------------------------------------------------------
# Role: [Exploratory] - Analyze plasmid replicon co-occurrence and ST associations.
#
# Inputs:
#   - results/mlst/mlst_all.tsv
#   - data/assemblies/*.fasta
#
# Outputs:
#   - results/plasmids/plasmidfinder_hits_long.csv
#   - results/plasmids/plasmidfinder_presence_absence.csv
#   - plots/plasmids/replicon_cooccurrence.pdf
#   - plots/plasmids/ST_vs_replicon_network.pdf
#   - results/plasmids/abricate_cache/ (cache)
#
# Usage:
#   Rscript 09_inc_plasmid_network.R
#
# Biological/Statistical purpose:
#   - Visualizes how plasmid replicons co-occur (e.g., multi-replicon plasmids).
#   - Maps plasmid content to bacterial lineages (STs).
# ==============================================================================

# 1. Load Configuration & Libraries
source("00_config.R")
suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
  library(stringr)
  library(furrr)
  library(fs)
  library(igraph)
  library(ggraph)
  library(ggplot2)
  library(purrr)
})

# 2. Configuration
# ------------------------------------------------------------------------------
invisible(check_tool("abricate"))
DIR_CACHE <- file.path(DIR_PLASMIDS, "abricate_cache")
ensure_dir(DIR_CACHE)
ensure_dir(DIR_PLOTS_PLASMIDS)

# 3. Run ABRicate (PlasmidFinder)
# ------------------------------------------------------------------------------
run_abricate <- function(fasta) {
  basename <- fs::path_file(fasta)
  out <- fs::path(DIR_CACHE, paste0(basename, ".tab"))

  if (file.exists(out)) {
    return(out)
  }

  res <- processx::run("abricate", c("--quiet", "--threads", "1", "--db", "plasmidfinder", fasta), echo = FALSE, error_on_status = FALSE)
  if (res$status == 0) {
    write_lines(res$stdout, out)
    return(out)
  }
  return(NA_character_)
}

fasta_files <- dir_ls(DIR_FASTAS, glob = "*.fasta")
if (length(fasta_files) == 0) fasta_files <- dir_ls(DIR_FASTAS, glob = "*.fa")

future::plan(future::multisession, workers = CORES_USE)
tab_files <- future_map_chr(fasta_files, run_abricate, .progress = TRUE)
future::plan(future::sequential)

tab_files <- tab_files[!is.na(tab_files)]

# 4. Tidy Data
# ------------------------------------------------------------------------------
read_one_tab <- function(f) {
  read_tsv(f, comment = "", col_types = cols(.default = "c"), show_col_types = FALSE) %>%
    mutate(file = fs::path_file(f))
}

hits_long <- map_dfr(tab_files, read_one_tab) %>%
  mutate(
    identity   = suppressWarnings(parse_number(`%IDENTITY`)),
    coverage   = suppressWarnings(parse_number(`%COVERAGE`)),
    Isolate_ID = tools::file_path_sans_ext(file)
  ) %>%
  select(Isolate_ID, accession = ACCESSION, identity, coverage, SEQUENCE, GENE)

write_csv(hits_long, file.path(DIR_PLASMIDS, "plasmidfinder_hits_long.csv"))

# Presence/Absence Matrix
matrix <- hits_long %>%
  filter(!is.na(accession)) %>%
  distinct(Isolate_ID, accession) %>%
  mutate(value = 1L) %>%
  pivot_wider(names_from = accession, values_from = value, values_fill = 0)

write_csv(matrix, file.path(DIR_PLASMIDS, "plasmidfinder_presence_absence.csv"))

# 5. Network Analysis
# ------------------------------------------------------------------------------
# A. Replicon Co-occurrence
num_mat <- matrix %>%
  select(-Isolate_ID) %>%
  as.matrix()
keep <- which(colSums(num_mat) >= 3) # Threshold: >= 3 isolates

if (length(keep) > 1) {
  mat <- num_mat[, keep, drop = FALSE]
  coocc_edges <- t(mat) %*% mat
  diag(coocc_edges) <- 0

  # Keep only upper triangle to avoid duplicate edges in undirected graph
  coocc_edges[lower.tri(coocc_edges)] <- 0

  edges_tbl <- as.data.frame(as.table(coocc_edges)) %>%
    filter(Freq > 0) %>%
    rename(from = Var1, to = Var2, weight = Freq)

  g1 <- graph_from_data_frame(edges_tbl, directed = FALSE)

  pdf(file.path(DIR_PLOTS_PLASMIDS, "replicon_cooccurrence.pdf"), width = 7, height = 7)
  print(ggraph(g1, layout = "nicely") +
    geom_edge_link(aes(width = weight), alpha = .6) +
    geom_node_point(size = 4) +
    geom_node_text(aes(label = name), vjust = 1.5) +
    theme_void() +
    ggtitle("Plasmid Replicon Co-occurrence Network"))
  dev.off()
}

# B. ST vs Replicon Bipartite Graph
FILE_MLST <- FILE_MLST_ALL
if (file.exists(FILE_MLST)) {
  mlst <- read_tsv(FILE_MLST, show_col_types = FALSE) %>% select(Isolate_ID, ST)

  bipartite_edges <- hits_long %>%
    left_join(mlst, by = "Isolate_ID") %>%
    filter(!is.na(ST), !is.na(accession)) %>%
    distinct(ST, accession)

  if (nrow(bipartite_edges) > 0) {
    g2 <- graph_from_data_frame(bipartite_edges, directed = FALSE)
    V(g2)$type <- str_detect(V(g2)$name, "^\\d+$") # Assumes STs are numeric strings

    pdf(file.path(DIR_PLOTS_PLASMIDS, "ST_vs_replicon_network.pdf"), width = 8, height = 6)
    print(ggraph(g2, layout = "fr") +
      geom_edge_link(alpha = .4) +
      geom_node_point(aes(color = type), size = 4) +
      scale_color_manual(values = c("steelblue", "tomato"), labels = c("Replicon", "ST"), name = "") +
      geom_node_text(aes(label = name), repel = TRUE, size = 3) +
      theme_void() +
      ggtitle("Bipartite Network: Sequence Types vs. Plasmid Replicons"))
    dev.off()
  }
}

msg("✓ Plasmid network analysis complete.")

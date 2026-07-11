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
#   - results/mlst/mlst_provider_preferred.csv
#   - results/qc/canonical_assembly_selection.csv
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

selection_file <- file.path(DIR_QC, "canonical_assembly_selection.csv")
if (file.exists(selection_file)) {
  fasta_manifest <- read_csv(selection_file, show_col_types = FALSE) %>%
    mutate(
      full_path = if ("full_path" %in% names(.)) as.character(full_path) else as.character(fasta_path),
      selected_canonical = if ("selected_canonical" %in% names(.)) as_pipeline_bool(selected_canonical) else FALSE,
      file_exists = if ("file_exists" %in% names(.)) as_pipeline_bool(file_exists, default = file.exists(full_path)) else file.exists(full_path),
      Isolate_ID = if ("Isolate_ID" %in% names(.)) as.character(Isolate_ID) else tools::file_path_sans_ext(basename(full_path)),
      Assembly_ID = if ("Assembly_ID" %in% names(.)) as.character(Assembly_ID) else tools::file_path_sans_ext(basename(full_path)),
      Participant_id = if ("Participant_id" %in% names(.)) as.character(Participant_id) else NA_character_,
      tp_lab = if ("tp_lab" %in% names(.)) normalise_timepoint_preserve_events(tp_lab) else NA_character_
    ) %>%
    filter(selected_canonical %in% TRUE, file_exists %in% TRUE, !is.na(full_path), file.exists(full_path)) %>%
    distinct(full_path, .keep_all = TRUE)
  msg("Using %d canonical selected FASTA(s) for PlasmidFinder.", nrow(fasta_manifest))
} else {
  warning("Canonical assembly selection not found; falling back to top-level FASTA scan.")
  fasta_files <- dir_ls(DIR_FASTAS, glob = "*.fasta")
  if (length(fasta_files) == 0) fasta_files <- dir_ls(DIR_FASTAS, glob = "*.fa")
  fasta_manifest <- tibble::tibble(
    full_path = as.character(fasta_files),
    Isolate_ID = tools::file_path_sans_ext(basename(fasta_files)),
    Assembly_ID = tools::file_path_sans_ext(basename(fasta_files)),
    Participant_id = NA_character_,
    tp_lab = NA_character_
  )
}

if (nrow(fasta_manifest) == 0) stop("No FASTA files available for PlasmidFinder.")
write_csv(fasta_manifest, file.path(DIR_PLASMIDS, "plasmidfinder_input_manifest.csv"))

future::plan(future::multisession, workers = CORES_USE)
tab_manifest <- fasta_manifest %>%
  mutate(tab_file = future_map_chr(full_path, run_abricate, .progress = TRUE)) %>%
  filter(!is.na(tab_file), file.exists(tab_file))
future::plan(future::sequential)

# 4. Tidy Data
# ------------------------------------------------------------------------------
read_one_tab <- function(tab_file, Isolate_ID, Assembly_ID, Participant_id, tp_lab) {
  if (!file.exists(tab_file) || file.size(tab_file) == 0) return(tibble::tibble())
  read_tsv(tab_file, comment = "", col_types = cols(.default = "c"), show_col_types = FALSE) %>%
    mutate(
      Isolate_ID = Isolate_ID,
      Assembly_ID = Assembly_ID,
      Participant_id = Participant_id,
      tp_lab = tp_lab,
      file = fs::path_file(tab_file)
    )
}

hits_raw <- pmap_dfr(
  tab_manifest %>% select(tab_file, Isolate_ID, Assembly_ID, Participant_id, tp_lab),
  read_one_tab
)

hits_long <- if (nrow(hits_raw) > 0 && all(c("ACCESSION", "%IDENTITY", "%COVERAGE", "SEQUENCE", "GENE") %in% names(hits_raw))) {
  hits_raw %>%
    mutate(
      identity = suppressWarnings(parse_number(`%IDENTITY`)),
      coverage = suppressWarnings(parse_number(`%COVERAGE`))
    ) %>%
    select(Isolate_ID, Assembly_ID, Participant_id, tp_lab,
           accession = ACCESSION, identity, coverage, SEQUENCE, GENE)
} else {
  tibble::tibble(
    Isolate_ID = character(), Assembly_ID = character(),
    Participant_id = character(), tp_lab = character(),
    accession = character(), identity = numeric(), coverage = numeric(),
    SEQUENCE = character(), GENE = character()
  )
}

write_csv(hits_long, file.path(DIR_PLASMIDS, "plasmidfinder_hits_long.csv"))

# Presence/Absence Matrix
matrix <- if (nrow(hits_long) > 0) {
  hits_long %>%
    filter(!is.na(accession)) %>%
    distinct(Isolate_ID, accession) %>%
    mutate(value = 1L) %>%
    pivot_wider(names_from = accession, values_from = value, values_fill = 0)
} else {
  tibble::tibble(Isolate_ID = character())
}

write_csv(matrix, file.path(DIR_PLASMIDS, "plasmidfinder_presence_absence.csv"))

# 5. Network Analysis
# ------------------------------------------------------------------------------
# A. Replicon Co-occurrence
num_mat <- matrix %>%
  select(-Isolate_ID) %>%
  as.matrix()
keep <- if (ncol(num_mat) > 0) which(colSums(num_mat) >= 3) else integer() # Threshold: >= 3 isolates

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
FILE_MLST <- FILE_MLST_CANONICAL
if (file.exists(FILE_MLST)) {
  mlst <- read_csv(FILE_MLST, show_col_types = FALSE) %>%
    select(any_of(c("Isolate_ID", "ST", "ST_source")))

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

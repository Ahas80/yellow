#!/usr/bin/env Rscript
# =============================================================
# 09_inc_plasmid_network.R  –  Inc typing via ABRicate + networks
# Yellow RoUTIne  ·  9 Jul 2025
# =============================================================

# ---- configuration ------------------------------------------
FASTA_DIR   <- "ont-yellow-routine-fastas"      # same as 06_MLST.R
THREADS     <- max(1L, parallel::detectCores() - 1L)
ABR_BIN     <- Sys.which("abricate")
DB          <- "plasmidfinder"                  # Inc typing
CACHE_DIR   <- "results/abricate_plasmid"
dir.create(CACHE_DIR, showWarnings = FALSE, recursive = TRUE)

# output files -------------------------------------------------
csv_long    <- "results/plasmidfinder_hits_long.csv"
csv_matrix  <- "results/plasmidfinder_presence_absence.csv"
net_pdf1    <- "results/replicon_cooccurrence.pdf"
net_pdf2    <- "results/ST_vs_replicon_network.pdf"

# ---- libraries ----------------------------------------------
suppressPackageStartupMessages({
  library(dplyr);  library(readr);  library(tidyr)
  library(stringr);library(furrr);  library(fs)
  library(igraph); library(ggraph); library(ggplot2)
})

if (ABR_BIN == "") stop("ABRicate not found on PATH")

# ---- 1 · run (cached) ---------------------------------------
run_abricate <- function(fasta) {
  out <- fs::path(CACHE_DIR, paste0(fs::path_file(fasta), ".tab"))
  if (fs::file_exists(out)) return(out)
  
  cmd <- c("--quiet", "--threads", "1", "--db", DB, fasta)
  res <- processx::run(ABR_BIN, cmd, echo = FALSE, error_on_status = FALSE)
  if (res$status != 0)
    warning("ABRicate failed on ", fasta, "; skipping") else
      write_lines(res$stdout, out)
  out
}

fasta_files <- fs::dir_ls(FASTA_DIR, glob = "*.fasta")
plan(multisession, workers = THREADS)
tab_files <- future_map_chr(fasta_files, run_abricate, .progress = TRUE)

## ---- 2 · read + tidy long ---------------------------------------------------
read_one_tab <- function(f) {
  read_tsv(f,
           comment = "",                 #  <-  DON'T drop the header!
           col_types = cols(.default = "c"),
           show_col_types = FALSE) |>
    mutate(file = fs::path_file(f))
}


hits_long <- purrr::map_dfr(tab_files, read_one_tab) |>
  mutate(
    identity = suppressWarnings(parse_number(`%IDENTITY`)),
    coverage = suppressWarnings(parse_number(`%COVERAGE`)),
    Isolate_ID = stringr::str_extract(file, "[0-9A-Za-z]+-[0-9]+")
  ) |>
  select(Isolate_ID, accession = ACCESSION, identity, coverage,
         SEQUENCE, GENE, everything())          # keep a few informative cols


write_csv(hits_long, csv_long)
message("✓ hits (long)  → ", csv_long)

## ---- 3 · presence/absence matrix -------------------------------------------
matrix <- hits_long |>
  filter(!is.na(accession)) |>
  distinct(Isolate_ID, accession) |>
  mutate(value = 1L) |>
  tidyr::pivot_wider(names_from = accession,
                     values_from = value,
                     values_fill = 0)

write_csv(matrix, csv_matrix)
message("✓ presence/absence → ", csv_matrix)

## ---- 4 · replicon co-occurrence --------------------------------------------
mat <- matrix %>%                        # 1) store the numeric matrix
  select(-Isolate_ID) %>%
  as.matrix()

coocc_edges <- t(mat) %*% mat            # 2) multiply once


# … the rest of the script (ggraph code, bipartite graph, etc.) stays unchanged …
# co-occur counts
diag(coocc_edges) <- 0                    # drop self loops

edges_tbl <- as.data.frame(as.table(coocc_edges)) %>%
  filter(Freq > 0) %>%
  rename(from = Var1, to = Var2, weight = Freq)

g1 <- graph_from_data_frame(edges_tbl, directed = FALSE)
pdf(net_pdf1, width = 7, height = 7)
ggraph(g1, layout = "fr") +
  geom_edge_link(aes(width = weight), alpha = .6) +
  geom_node_point(size = 4) +
  geom_node_text(aes(label = name), vjust = 1.5) +
  theme_void() +
  ggtitle("Replicon co-occurrence (Inc types)")
dev.off()
message("✓ network plot 1 → ", net_pdf1)

# ---- 5 · ST vs replicon bipartite graph ---------------------
mlst <- read_tsv("results/mlst/mlst_all.tsv", show_col_types = FALSE) %>%
  select(Isolate_ID, ST)

bipartite_edges <- hits_long %>%
  left_join(mlst, by = "Isolate_ID") %>%
  distinct(ST, accession)

g2 <- graph_from_data_frame(bipartite_edges, directed = FALSE)
V(g2)$type <- str_detect(V(g2)$name, "^\\d+$")  # TRUE for ST nodes

pdf(net_pdf2, width = 8, height = 6)
ggraph(g2, layout = "fr") +
  geom_edge_link(alpha = .4) +
  geom_node_point(aes(color = type), size = 4) +
  scale_color_manual(values = c("steelblue", "tomato"),
                     labels = c("Replicon", "ST"),
                     name = "") +
  geom_node_text(aes(label = name), repel = TRUE, size = 3) +
  theme_void() +
  ggtitle("Chromosomal ST vs Inc replicon network")
dev.off()
message("✓ network plot 2 → ", net_pdf2)

message("09_inc_plasmid_network.R finished ✅")

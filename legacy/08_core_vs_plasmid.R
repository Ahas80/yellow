#!/usr/bin/env Rscript
# =============================================================
# 08_core_vs_plasmid_STs.R – chromosomal STs vs plasmid pMLST
# Parallel version (furrr) – 9 Jul 2025
# =============================================================

suppressPackageStartupMessages({
  library(dplyr);  library(readr);  library(tidyr)
  library(purrr);   library(furrr);  library(fs)
  library(stringr); library(scales); library(processx)
})

## ---- 0 · paths -------------------------------------------------------------
fa_dir        <- "ont-yellow-routine-fastas"
core_mlst_tsv <- "results/mlst/mlst_all.tsv"

out_core_csv  <- "results/mlst/ST_core_freq.csv"
out_pmlst_long<- "results/mlst/pMLST_hits_long.csv"
out_pmlst_wide<- "results/mlst/plasmid_types_per_isolate.csv"

## ---- 1 · chromosomal ST frequencies ---------------------------------------
core_tbl <- read_tsv(core_mlst_tsv, show_col_types = FALSE)

core_tbl %>%
  count(ST, sort = TRUE) %>%
  mutate(pct = percent(n / sum(n))) %>%
  write_csv(out_core_csv)

message("✓ chromosomal ST table → ", out_core_csv)

## ---- 2 · pMLST setup -------------------------------------------------------
MLST_BIN <- Sys.which("mlst")
if (MLST_BIN == "") stop("mlst not found on PATH")

pSchemes    <- c("ac", "f", "hi1", "hi2", "i1", "n")
fasta_files <- fs::dir_ls(fa_dir, glob = "*.fasta")

run_one <- function(fasta, scheme) {
  out <- fs::path("results/mlst/raw",
                  paste0(fs::path_file(fasta), ".", scheme, ".csv"))
  if (fs::file_exists(out))                       # cached
    return(read_csv(out, col_types = cols(.default = "c")))
  
  res <- processx::run(MLST_BIN,
                       c("--quiet", "--threads", "1",
                         "--scheme", scheme, "--csv", "--legacy", fasta),
                       error_on_status = FALSE)
  
  if (res$status != 0 || !str_detect(res$stdout, ",")) return(NULL)
  
  dat <- read_csv(I(res$stdout), col_types = cols(.default = "c")) |>
    rename_with(tolower)
  write_csv(dat, out)                             # cache
  dat
}

## ---- 3 · ETA quick test ----------------------------------------------------
test_time <- system.time(run_one(fasta_files[1], pSchemes[1]))[["elapsed"]]
est_min   <- test_time * length(fasta_files) * length(pSchemes) /
  (60 * 4)                   # 4 parallel workers below
message(sprintf("~%.1f min expected on 4 workers (%.1f s/test)",
                est_min, test_time))

## ---- 4 · parallel fan-out --------------------------------------------------
plan(multisession, workers = 4)          # ← adjust cores here

grid <- tidyr::expand_grid(file = fasta_files, scheme = pSchemes)

pmlst_hits <- grid |>
  mutate(dat = furrr::future_pmap(
    list(file, scheme), run_one,
    .progress = TRUE)) |>
  unnest(dat)                            # drops NULL automatically

if (!nrow(pmlst_hits)) {
  message("No plasmid MLST alleles detected.")
  quit(status = 0)
}

## ---- 5 · write long & wide --------------------------------------------------
write_csv(pmlst_hits, out_pmlst_long)
message("✓ pMLST long table     → ", out_pmlst_long)

pmlst_hits %>%
  select(isolate_id, p_scheme = scheme, st) %>%
  pivot_wider(names_from = p_scheme,
              values_from = st,
              values_fill = NA_character_,
              names_prefix = "pST_") %>%
  arrange(isolate_id) %>%
  write_csv(out_pmlst_wide)

message("✓ per-isolate pST wide → ", out_pmlst_wide)

## ---- 6 · console top-5 summary --------------------------------------------
cat("\n---------  pMLST summary  ---------------------------------\n")
pmlst_hits %>%
  count(scheme, st, sort = TRUE) %>%
  group_by(scheme) %>%
  slice_max(n, n = 5, with_ties = FALSE) %>%
  ungroup() %>%
  print(n = 30)

cat("\nFinished 08_core_vs_plasmid_STs.R ✅\n")

#!/usr/bin/env Rscript
# =============================================================
# 10_visualise_STs_pMLST.R  –  four quick plots
#   1) Chromosomal ST bar chart (top 20)
#   2) pMLST Inc-type bar chart (top 20 across all schemes)
#   3) Presence/absence heat-map  (Isolates × Inc schemes)
#   4) Co-occurrence network  (Inc replicons that travel together)
# Yellow RoUTIne – 9 Jul 2025
# =============================================================

suppressPackageStartupMessages({
  library(dplyr);      library(readr);  library(tidyr)
  library(ggplot2);    library(scales)
  library(igraph);     library(ggraph)
  library(ComplexHeatmap)   # install.packages("ComplexHeatmap") if missing
})

# ---------- 0 · paths --------------------------------------------------------
core_csv   <- "results/mlst/ST_core_freq.csv"
p_long_csv <- "results/mlst/pMLST_hits_long.csv"
p_wide_csv <- "results/mlst/plasmid_types_per_isolate.csv"

dir.create("results/plots", showWarnings = FALSE, recursive = TRUE)

# ---------- 1 · chromosomal ST bar-chart ------------------------------------
st_freq <- read_csv(core_csv, show_col_types = FALSE) %>% slice_max(n, n = 20)

ggplot(st_freq, aes(reorder(ST, n), n)) +
  geom_col(fill = "steelblue") +
  coord_flip() +
  labs(x = "Chromosomal ST", y = "Isolates",
       title = "Top 20 chromosomal STs") +
  theme_classic()

ggsave("results/plots/top20_core_STs.pdf", width = 7, height = 5)

# ---------- 2 · pMLST Inc-type bar-chart ------------------------------------
p_long <- read_csv(p_long_csv, show_col_types = FALSE)

p_freq <- p_long %>%
  count(scheme, st, sort = TRUE) %>%
  unite("pST", scheme, st, sep = "_") %>%
  slice_max(n, n = 20)

ggplot(p_freq, aes(reorder(pST, n), n)) +
  geom_col(fill = "tomato") +
  coord_flip() +
  labs(x = "Plasmid ST (scheme_st)", y = "Isolates",
       title = "Top 20 plasmid STs across all Inc schemes") +
  theme_classic()

ggsave("results/plots/top20_pMLST.pdf", width = 8, height = 5)

# ---------- 3 · presence/absence heat-map -----------------------------------
mat <- read_csv(p_wide_csv, show_col_types = FALSE) %>%
  column_to_rownames("Isolate_ID") %>%
  mutate(across(everything(), ~ ifelse(is.na(.x), 0, 1))) %>%
  as.matrix()

Heatmap(mat,
        name = "Inc\npresent",
        col  = c("0" = "white", "1" = "darkorchid"),
        show_row_names = FALSE,
        column_title  = "Inc schemes per isolate") |
  
  pdf("results/plots/Inc_presence_heatmap.pdf", width = 7, height = 9)
draw(Heatmap(mat, name = "Inc\npresent",
             col = c("0" = "white", "1" = "darkorchid"),
             show_row_names = FALSE,
             column_title = "Inc schemes per isolate"))
dev.off()

# ---------- 4 · replicon co-occurrence network ------------------------------
# build adjacency
adj <- t(mat) %*% mat; diag(adj) <- 0
edge_tbl <- as.data.frame(as.table(adj)) %>%
  filter(Freq > 0) %>%
  rename(from = Var1, to = Var2, weight = Freq)

g <- graph_from_data_frame(edge_tbl, directed = FALSE)

pdf("results/plots/Inc_cooccurrence_network.pdf", width = 7, height = 7)
ggraph(g, layout = "fr") +
  geom_edge_link(aes(width = weight), colour = "grey60", alpha = .5) +
  scale_edge_width(range = c(.5, 3)) +
  geom_node_point(size = 6, colour = "goldenrod") +
  geom_node_text(aes(label = name), repel = TRUE) +
  theme_void() +
  ggtitle("Inc-type co-occurrence network")
dev.off()

cat("Plots saved to results/plots/ … ✅\n")

mlst <- read_tsv("mlst_all.tsv") %>% select(Isolate_ID, ST)
merged <- pf %>% select(Isolate_ID, accession) %>% left_join(mlst)
merged %>% count(ST, accession) %>% 
  filter(n >= 5) %>%            # tweak threshold
  arrange(desc(n))


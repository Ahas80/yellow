# at the top of your script -------------------------------------------------




#!/usr/bin/env Rscript
# =============================================================
# Longitudinal FASTA comparison | participants with 3-5 timepoints
# =============================================================

# ---- LOAD LIBRARIES ----------------------------------------------------------
suppressPackageStartupMessages({
  library(dplyr);   library(stringr);   library(readr)
  library(tidyr);   library(ggplot2);   library(Biostrings)
  library(glue);    library(furrr);     library(patchwork)
})

# ---- DIRECTORY ---------------------------------------------------------------
script_dir <- getwd()          # adjust if you prefer another folder
fasta_dir  <- file.path(script_dir, "ont-yellow-routine-fastas")

# ---- LOAD METADATA -----------------------------------------------------------
batch1 <- read_csv("batch1.csv", show_col_types = FALSE)
batch2 <- read_csv("batch2.csv", show_col_types = FALSE)
common_cols <- intersect(names(batch1), names(batch2))
metadata <- bind_rows(
  mutate(batch1, across(all_of(common_cols), as.character)),
  mutate(batch2, across(all_of(common_cols), as.character))
)

# ---- FASTA SCAN & BASIC METRICS ---------------------------------------------
all_fasta_paths <- list.files(fasta_dir, "\\.fasta$", full.names = TRUE)

assembly_info <- tibble(full_path = all_fasta_paths) %>%
  mutate(file_name  = basename(full_path),
         Isolate_ID = str_extract(sapply(strsplit(file_name, "_"), `[`, 3),
                                  "[0-9A-Za-z]+-[0-9]+"),
         assembler  = case_when(str_detect(file_name, "flye")        ~ "flye",
                                str_detect(file_name, "longcycler") ~ "longcycler",
                                TRUE                                ~ "unknown"))

summarize_fasta <- function(p) {
  s <- readDNAStringSet(p)
  af <- colSums(alphabetFrequency(s, baseOnly = TRUE))
  tibble(num_contigs = length(s),
         total_bases = sum(width(s)),
         gc_content  = round((af["G"] + af["C"]) / sum(af) * 100, 2))
}

assembly_df <- assembly_info %>%
  rowwise() %>% mutate(summarize_fasta(full_path)) %>% ungroup() %>%
  left_join(metadata, by = c("Isolate_ID" = "isolate_ID"))


#added code 
###############################################################################
# >>>  FULL-COHORT EXTENSION  <<<   (add below your existing pipeline)
###############################################################################
message("── Scaling up to ALL participants & ALL time-points")

# ------------------------------------------------------------------
# 1.  ABRICATE  on every FASTA (cached)   →  vf_hits_all
# ------------------------------------------------------------------
dir.create("results/abricate", showWarnings = FALSE, recursive = TRUE)

run_abricate_cached <- function(fasta) {
  cache <- file.path("results/abricate",
                     paste0(basename(fasta), ".vfdb.tsv"))
  if (file.exists(cache)) {
    return(readr::read_tsv(cache, show_col_types = FALSE))
  }
  cmd <- glue::glue(
    "abricate --quiet --db vfdb --mincov 70 --minid 70 {shQuote(fasta)} > {shQuote(cache)}"
  )
  system(cmd)
  readr::read_tsv(cache, show_col_types = FALSE)
}

safe_abr <- purrr::safely(run_abricate_cached, otherwise = NULL, quiet = TRUE)

vf_hits_all <- assembly_df %>%                           # <-- already exists
  mutate(vfdb = furrr::future_map(full_path,
                                  ~ safe_abr(.x)$result,
                                  .progress = TRUE)) %>%
  tidyr::unnest(vfdb)                                    # long table

saveRDS(vf_hits_all, "results/vf_hits_all.rds")
message("✓   vf_hits_all saved (", nrow(vf_hits_all), " rows)")

# ------------------------------------------------------------------
# 2.  Presence/absence matrix  →  vf_pa_all
# ------------------------------------------------------------------
vf_pa_all <- vf_hits_all %>%
  distinct(Participant_id, Timepoint, GENE) %>%
  mutate(present = 1) %>%
  tidyr::pivot_wider(names_from = GENE,
                     values_from = present,
                     values_fill = 0)

readr::write_csv(vf_pa_all, "results/vf_pa_all.csv")
message("✓   vf_pa_all matrix written (",
        nrow(vf_pa_all), " samples × ",
        ncol(vf_pa_all) - 2, " genes)")

# ------------------------------------------------------------------
# 3.  Gene-level prevalence for NEW bar / histogram plots
# ------------------------------------------------------------------
tbl_gene <- vf_hits_all %>%
  distinct(Participant_id, GENE) %>%
  count(GENE, name = "n_participants") %>%
  left_join(vf_hits_all %>% count(GENE, name = "total_hits"),
            by = "GENE") %>%
  arrange(desc(n_participants))

readr::write_csv(tbl_gene, "results/stats_gene_level.csv")
message("✓   stats_gene_level.csv updated")

# ------------------------------------------------------------------
# 4.  Longitudinal subset (≥ 2 TP) for nucmer & gain/loss
# ------------------------------------------------------------------
ids_multi <- vf_hits_all %>%
  distinct(Participant_id, Timepoint) %>%
  add_count(Participant_id, name = "n_tp") %>%
  filter(n_tp >= 2) %>%
  pull(Participant_id)

if (length(ids_multi) > 0) {
  message("→  Running nucmer for ", length(ids_multi),
          " participants with ≥ 2 time-points…")
  
  assembly_long <- assembly_df %>% filter(Participant_id %in% ids_multi)
  vf_hits_long  <- vf_hits_all %>%  filter(Participant_id %in% ids_multi)
  
  ## --- re-use YOUR EXISTING nucmer code but point it at assembly_long ----
  ## (basically copy the block that builds `pair_df` / `pairwise_annotated`,
  ##  replacing `assembly_df` with `assembly_long`)
  ##
  ## Example stub:
  # comparable_tbl <- assembly_long %>% ...  # same logic as before
  # pair_df        <- comparable_tbl   %>% ...  # same logic
  # pairwise_annotated <- pair_df       %>% ...  # same logic
  
} else {
  warning("No participants with ≥ 2 time-points – nucmer skipped.")
}

# ------------------------------------------------------------------
# 5.  Plots that use the bigger tables
# ------------------------------------------------------------------
##   (A)  core/accessory bar + histogram  ---------------------------
library(ggplot2)

top25 <- tbl_gene %>%
  slice_max(n_participants, n = 25) %>%
  mutate(GENE = forcats::fct_reorder(GENE, n_participants))

ggplot(top25, aes(GENE, n_participants)) +
  geom_col(fill = "steelblue") +
  coord_flip() +
  labs(title = "Top 25 most prevalent genes (all isolates)",
       y = "Participants", x = NULL) +
  theme_minimal(base_size = 12)
ggsave("results/plots/core_bar_top25_all.png", width = 6, height = 6)

ggplot(tbl_gene, aes(n_participants)) +
  geom_histogram(binwidth = 1, fill = "grey70") +
  labs(title = "Gene prevalence across entire cohort",
       x = "# participants carrying gene", y = "Gene count") +
  theme_minimal(base_size = 11)
ggsave("results/plots/core_histogram_all.png", width = 5, height = 4)

##   (B)  Heat-map & Upset still work; just replace vf_pa with vf_pa_all
##   (reuse your functions: plot_heatmap_var, get_upset_df, etc.)

message("✓   Full-cohort extension completed")
###############################################################################
# <<<  end of scale-up block >>>  
###############################################################################




# ---- SELECT PARTICIPANTS WITH 3-5 TIMEPOINTS ---------------------------------
id_tbl <- assembly_df %>%
  filter(!is.na(Timepoint), assembler %in% c("flye", "longcycler")) %>%
  add_count(Participant_id, assembler, name = "n_tp") %>%
  filter(n_tp %in% 3:5)

selected_ids <- unique(id_tbl$Participant_id)

assembly_df <- assembly_df %>% filter(Participant_id %in% selected_ids)

# ---- RUN NUCMER / DNADIFF ----------------------------------------------------
plan(multisession, workers = max(1, parallel::detectCores() - 1))

comparable_tbl <- assembly_df %>%
  group_by(Participant_id, assembler) %>%
  arrange(Timepoint, .by_group = TRUE) %>%
  mutate(path_A  = full_path,
         path_B  = lead(full_path),
         tp_B    = lead(Timepoint),
         out_dir = glue("results/{Participant_id}_{assembler}_{Timepoint}_vs_{tp_B}")) %>%
  filter(!is.na(path_B)) %>% ungroup()

run_nucmer <- function(a, b, outdir){
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  pref <- file.path(outdir, "run")
  system(glue("nucmer --mum --prefix={shQuote(pref)} {shQuote(a)} {shQuote(b)}"))
  system(glue("delta-filter -1 {pref}.delta > {pref}.1delta"))
  system(glue("dnadiff -p {pref}_dd -d {pref}.1delta"))
  rpt <- glue("{pref}_dd.report")
  if(!file.exists(rpt)) return(tibble(AvgIdentity = NA, TotalSnpCnt = NA))
  L  <- readLines(rpt)
  grab <- function(k) { m <- grep(k, L, value = TRUE); if(length(m)) as.numeric(str_extract(m[1], "\\d+\\.?\\d*")) else NA }
  tibble(AvgIdentity = grab("AvgIdentity"), TotalSnpCnt = grab("TotalSNPs"))
}

pair_df <- comparable_tbl %>%
  mutate(result = future_pmap(list(path_A, path_B, out_dir), run_nucmer,
                              .progress = TRUE)) %>%
  unnest(result) %>%
  mutate(tp_A = str_extract(out_dir, "(?<=_)[A-Za-z0-9]+(?=_vs)"))

# ---- ASB / UTI FLAG TABLE  (one row per P × A × Timepoint) -------------------
assembly_meta_flags <- assembly_df %>% 
  mutate(UTI_flag = ifelse(`No S&S` == 1, "ASB", "UTI_suspected")) %>% 
  distinct(Participant_id, assembler, Timepoint, .keep_all = TRUE) %>% 
  select(Participant_id, assembler, Timepoint, UTI_flag)


# ---- JOIN FLAG INTO PAIRS -------------------------------------------------
pairwise_annotated <- pair_df %>% 
  left_join(assembly_meta_flags,
            by = c("Participant_id", "assembler", "tp_A" = "Timepoint")) %>% 
  mutate(UTI_flag_A = UTI_flag) %>%      # copy to a stable name
  select(-UTI_flag)                      # drop the original



# ---- PLOT FUNCTION -----------------------------------------------------------
plot_pairwise <- function(df, pid, asm){
  sub <- df %>% filter(Participant_id == pid, assembler == asm) %>% arrange(tp_A)
  if(nrow(sub) == 0) return(NULL)
  label <- paste0("P", pid, " (", asm, ")")
  
  p_id  <- ggplot(sub, aes(tp_A, AvgIdentity, group = 1)) +
    geom_line(colour = "blue") +
    geom_point(aes(color = UTI_flag_A), size = 3) +
    labs(title = paste("Avg Identity —", label),
         x = "From timepoint", y = "Identity (%)") +
    theme_minimal()
  
  p_snp <- ggplot(sub, aes(tp_A, TotalSnpCnt, group = 1)) +
    geom_line(colour = "red") +
    geom_point(aes(color = UTI_flag_A), size = 3) +
    labs(title = paste("Total SNPs —", label),
         x = "From timepoint", y = "SNP count") +
    theme_minimal()
  
  p_id / p_snp            # patchwork vertical combine
}

# ---- DISPLAY PLOTS -----------------------------------------------------------
for(pid in selected_ids){
  for(asm in c("flye", "longcycler")){
    print(plot_pairwise(pairwise_annotated, pid, asm))
  }
}

plot_dir <- file.path(script_dir, "plots_pairwise")
dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)

for (pid in selected_ids) {
  for (asm in c("flye", "longcycler")) {
    
    p <- plot_pairwise(pairwise_annotated, pid, asm)
    
    if (!is.null(p)) {
      # build a concise file name, e.g. P20034_flye.png
      f_out <- file.path(
        plot_dir,
        paste0("P", pid, "_", asm, ".png")
      )
      
      ggsave(
        filename   = f_out,
        plot       = p,
        width      = 8,       # inches
        height     = 6,       # inches
        dpi        = 300
      )
    }
  }
}


# # ----------------------------------------------------------
# # 0.  IDs we do NOT want plotted any more
# # ----------------------------------------------------------
# drop_ids <- c("20002", "60020", "100042")   # ← the 3 you said to remove
# 
# # ----------------------------------------------------------
# # 1.  Build the list of IDs to plot
# #     (after pairwise_annotated has been created)
# # ----------------------------------------------------------
# plot_ids <- setdiff(unique(pairwise_annotated$Participant_id), drop_ids)
# 
# # ----------------------------------------------------------
# # 2.  Plot-export loop  (unchanged except for the ID filter)
# # ----------------------------------------------------------
# for (id in plot_ids) {
#   for (asm in c("flye", "longcycler")) {
#     
#     p <- plot_pairwise_comparison(pairwise_annotated, id, asm)
#     
#     ## ---- save to disk ----
#     fname <- glue::glue("plots/P{id}_{asm}.png")
#     ggsave(filename = fname,
#            plot      = p,
#            width     = 8,
#            height    = 6,
#            dpi       = 300)
#   }
# }
# 
# # optional clean-up: if the old PNGs still exist, delete them
# file.remove(list.files("plots",
#                        pattern = "^P(20002|60020|100042)_(flye|longcycler)\\.png$",
#                        full.names = TRUE))
# 



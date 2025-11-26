#!/usr/bin/env Rscript

# LOAD LIBRARIES -----------------------------------------------------------------
library(dplyr)
library(stringr)
library(Biostrings)
library(glue)
library(furrr)
library(tidyr)
library(readr)
library(purrr)

# SET DIRECTORY ------------------------------------------------------------------
script_dir <- tryCatch({ dirname(normalizePath(sys.frame(1)$ofile)) }, error = function(e) { message("\u26a0\ufe0f Could not detect script path. Using current working directory."); getwd() })
setwd(script_dir)

# LOAD METADATA ------------------------------------------------------------------
batch1 <- read.csv("batch1.csv", stringsAsFactors = FALSE)
batch2 <- read.csv("batch2.csv", stringsAsFactors = FALSE)
common_cols <- intersect(names(batch1), names(batch2))
batch1 <- batch1 %>% mutate(across(all_of(common_cols), as.character))
batch2 <- batch2 %>% mutate(across(all_of(common_cols), as.character))
metadata <- bind_rows(batch1, batch2)

# SCAN FASTAs --------------------------------------------------------------------
fasta_dir <- file.path(script_dir, "ont-yellow-routine-fastas")
all_fasta_paths <- list.files(fasta_dir, pattern = "\\.fasta$", full.names = TRUE)

assembly_info <- tibble(full_path = all_fasta_paths) %>%
  mutate(file_name = basename(full_path),
         Isolate_ID = str_extract(sapply(strsplit(file_name, "_"), `[`, 3), "[0-9A-Za-z]+-[0-9]+"),
         assembler = case_when(str_detect(file_name, "flye") ~ "flye", str_detect(file_name, "longcycler") ~ "longcycler", TRUE ~ "unknown"))

summarize_fasta <- function(fasta_path) {
  seqs <- readDNAStringSet(fasta_path)
  alph_freqs <- colSums(alphabetFrequency(seqs, baseOnly = TRUE))
  gc_content <- if (sum(alph_freqs) > 0) round((alph_freqs["G"] + alph_freqs["C"]) / sum(alph_freqs) * 100, 2) else NA
  list(num_contigs = length(seqs), total_bases = sum(width(seqs)), gc_content = gc_content)
}

assembly_info_with_metrics <- assembly_info %>%
  rowwise() %>%
  mutate(metrics = list(summarize_fasta(full_path))) %>%
  unnest_wider(metrics) %>%
  ungroup()

# MERGE WITH METADATA -----------------------------------------------------------
assembly_df <- assembly_info_with_metrics %>%
  left_join(metadata, by = c("Isolate_ID" = "isolate_ID"))

# FILTER FOR PAIRWISE COMPARISONS ----------------------------------------------
plan(multisession, workers = max(1, parallel::detectCores() - 1))

comparable_tbl <- assembly_df %>%
  filter(!is.na(Timepoint), assembler %in% c("flye", "longcycler")) %>%
  group_by(Participant_id, assembler) %>%
  filter(n_distinct(Timepoint) >= 2) %>%
  arrange(Timepoint, .by_group = TRUE) %>%
  ungroup()

run_nucmer <- function(path_A, path_B, out_dir) {
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  prefix <- file.path(out_dir, "run")
  system(glue("nucmer --mum --prefix={shQuote(prefix)} {shQuote(path_A)} {shQuote(path_B)}"))
  system(glue("delta-filter -1 {prefix}.delta > {prefix}.1delta"))
  system(glue("dnadiff -p {prefix}_dd -d {prefix}.1delta"))
  report_path <- glue("{prefix}_dd.report")
  if (!file.exists(report_path)) return(tibble(AvgIdentity = NA, TotalSnpCnt = NA))
  
  lines <- readLines(report_path)
  extract_metric <- function(metric) {
    match <- grep(metric, lines, value = TRUE)
    if (length(match) == 0) return(NA)
    as.numeric(str_extract(match[1], "\\d+\\.\\d+|\\d+"))
  }
  
  tibble(
    AvgIdentity = extract_metric("AvgIdentity"),
    TotalSnpCnt = extract_metric("TotalSNPs")
  )
}

pair_df <- comparable_tbl %>%
  group_by(Participant_id, assembler) %>%
  mutate(path_A = full_path,
         path_B = lead(full_path),
         tp_B = lead(Timepoint),
         out_dir = glue("results/{Participant_id}_{assembler}_{Timepoint}_vs_{tp_B}")) %>%
  filter(!is.na(path_B)) %>%
  ungroup() %>%
  mutate(result = future_pmap(list(path_A, path_B, out_dir), run_nucmer, .progress = TRUE)) %>%
  unnest(result)

# FINAL OUTPUT ------------------------------------------------------------------
write_csv(assembly_df, "assembly_info_update1.csv")
write_csv(pair_df, "results/pairwise_stats.csv")


pairwise <- read.csv("results/pairwise_stats.csv")
dim(pairwise)
head(pairwise)

summary(pairwise$AvgIdentity)
summary(pairwise$TotalSnpCnt)
str(pairwise$AvgIdentity)
str(pairwise$TotalSnpCnt)
any(is.na(pairwise$AvgIdentity))
any(is.na(pairwise$TotalSnpCnt))
pairwise %>%
  group_by(Participant_id, assembler) %>%
  summarise(mean_identity = mean(AvgIdentity, na.rm = TRUE),
            total_snps = sum(TotalSnpCnt, na.rm = TRUE),
            .groups = "drop")
#add tp_A as we dont have the column 
pairwise <- pairwise %>% 
  mutate(
    tp_A = str_extract(out_dir, "(?<=_)[[:alnum:]]+_vs") %>% str_remove("_vs")
  )

# ============================
# FASTA COMPARISON PLOTS
# ============================

# Prep: merge ASB/UTI flags into pairwise data
# Use the fully merged assembly + metadata table
assembly_meta_flags <- assembly_info_update1 %>%
  select(Timepoint, Participant_id, assembler, No.S.S) %>%
  mutate(UTI_flag = ifelse(No.S.S == 1, "ASB", "UTI_suspected"))



# Extract tp_A from out_dir (already in your environment as pair_df)
pair_df <- pair_df %>%
  mutate(tp_A = str_extract(out_dir, "(?<=_)[A-Za-z0-9]+(?=_vs)"))

# Merge ASB/UTI flag from tp_A
assembly_meta_flags <- assembly_info_update1 %>%
  select(Timepoint, Participant_id, assembler, No.S.S) %>%
  mutate(UTI_flag = ifelse(No.S.S == 1, "ASB", "UTI_suspected"))

pairwise_annotated <- pair_df %>%
  left_join(assembly_meta_flags,
            by = c("Participant_id", "assembler", "tp_A" = "Timepoint")) %>%
  rename(UTI_flag_A = UTI_flag)

assembly_meta_flags <- assembly_info_update1 %>%
  select(Participant_id, assembler, Timepoint, No.S.S, total_bases) %>%
  mutate(
    total_bases = as.numeric(total_bases),
    UTI_flag = ifelse(No.S.S == 1, "ASB", "UTI_suspected")
  ) %>%
  group_by(Participant_id, assembler, Timepoint) %>%
  arrange(desc(total_bases)) %>%
  slice(1) %>%
  ungroup() %>%
  select(Participant_id, assembler, Timepoint, UTI_flag)



# Function to plot pairwise identity & SNPs over time
plot_pairwise_comparison <- function(df, participant, assembler) {
  df_sub <- df %>%
    filter(Participant_id == participant, assembler == assembler) %>%
    arrange(tp_A)
  
  p_id <- paste0("P", participant, " (", assembler, ")")
  
  p1 <- ggplot(df_sub, aes(x = tp_A, y = AvgIdentity, group = 1)) +
    geom_line(color = "blue") +
    geom_point(aes(color = UTI_flag_A), size = 3) +
    labs(title = paste("Avg Identity —", p_id),
         x = "From Timepoint", y = "Avg Identity (%)") +
    theme_minimal()
  
  p2 <- ggplot(df_sub, aes(x = tp_A, y = TotalSnpCnt, group = 1)) +
    geom_line(color = "red") +
    geom_point(aes(color = UTI_flag_A), size = 3) +
    labs(title = paste("Total SNPs —", p_id),
         x = "From Timepoint", y = "SNP Count") +
    theme_minimal()
  
  p1 / p2
}

# Loop over selected participants
for (id in selected_ids) {
  for (asm in c("flye", "longcycler")) {
    print(plot_pairwise_comparison(pairwise_annotated, id, asm))
  }
}


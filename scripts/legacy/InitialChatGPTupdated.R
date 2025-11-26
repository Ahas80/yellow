#!/usr/bin/env Rscript

###############################################################################
# analyze_assemblies.R
#
# This script:
#   1) Merges batch1 + batch2 metadata
#   2) Discovers all FASTA files in a folder, extracts Isolate_ID and assembler
#   3) Uses Biostrings to compute number of contigs, total bases, GC% for each
#   4) Merges FASTA metrics with metadata
#   5) Final "assembly_info_update1" includes everything
###############################################################################

# (A) LOAD NECESSARY LIBRARIES ------------------------------------------------

library(dplyr)      # For data wrangling
library(stringr)    # For string manipulation
library(Biostrings) # For readDNAStringSet(), if not installed:
# if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
# BiocManager::install("Biostrings")

# (B) PREP YOUR BATCH FILES (batch1, batch2) ----------------------------------
# Below assumes you have data frames 'batch1' and 'batch2' in your environment.
# If needed, load them:
# batch1 <- read.csv("path/to/batch1.csv")   # or readRDS, etc.
# batch2 <- read.csv("path/to/batch2.csv")

# Convert mismatched columns to character in batch1

batch1<-read.csv("batch1.csv")
batch2<-read.csv("batch2.csv")
batch1 <- batch1 %>%
  mutate(
    Participant_id = as.character(Participant_id),
    Spec = as.character(Spec),
    UWI. = as.character(UWI.),
    S.S.Dysuria = as.character(S.S.Dysuria),
    S.S.Urgency = as.character(S.S.Urgency),
    S.S.Frequency = as.character(S.S.Frequency),
    S.S.Incontinence = as.character(S.S.Incontinence),
    S.S.Pus = as.character(S.S.Pus),
    S.S.Flank_pain = as.character(S.S.Flank_pain),
    S.S.Suprapubic.pain = as.character(S.S.Suprapubic.pain),
    S.S.Fever = as.character(S.S.Fever),
    S.S.Chills = as.character(S.S.Chills),
    S.S.Delirium = as.character(S.S.Delirium),
    S.S.Other = as.character(S.S.Other),
    No.S.S = as.character(No.S.S)
  )

# Convert mismatched columns to character in batch2
batch2 <- batch2 %>%
  mutate(
    Participant_id = as.character(Participant_id),
    Spec = as.character(Spec),
    UWI. = as.character(UWI.),
    S.S.Dysuria = as.character(S.S.Dysuria),
    S.S.Urgency = as.character(S.S.Urgency),
    S.S.Frequency = as.character(S.S.Frequency),
    S.S.Incontinence = as.character(S.S.Incontinence),
    S.S.Pus = as.character(S.S.Pus),
    S.S.Flank_pain = as.character(S.S.Flank_pain),
    S.S.Suprapubic.pain = as.character(S.S.Suprapubic.pain),
    S.S.Fever = as.character(S.S.Fever),
    S.S.Chills = as.character(S.S.Chills),
    S.S.Delirium = as.character(S.S.Delirium),
    S.S.Other = as.character(S.S.Other),
    No.S.S = as.character(No.S.S)
  )

# Combine batch1 and batch2 into one metadata data frame
metadata <- bind_rows(batch1, batch2)

# (C) FIND AND SUMMARIZE ALL FASTA FILES --------------------------------------

# Define the directory path where the FASTAs live
dir_path <- "/Users/Aamir/Desktop/rUTIs/ont-yellow-routine-fastas"

# List all .fasta files in the directory
all_fasta_paths <- list.files(
  path       = dir_path, 
  pattern    = "\\.fasta$",
  full.names = TRUE
)

# 1. Create a data frame with initial assembly info:
assembly_info <- tibble(
  full_path = all_fasta_paths,
  file_name = basename(all_fasta_paths)
) %>%
  mutate(
    # Extract Isolate ID (including the -1)
    Isolate_ID = str_extract(sapply(strsplit(file_name, "_"), `[`, 3), "[0-9A-Za-z]+-[0-9]+"),
    # Identify the assembler type
    assembler = case_when(
      str_detect(file_name, "flye")        ~ "flye",
      str_detect(file_name, "longcycler") ~ "longcycler",
      TRUE                                ~ "unknown"
    )
  )

# 2. Define a function to read a FASTA and compute some metrics
summarize_fasta <- function(fasta_path) {
  seqs <- readDNAStringSet(fasta_path)
  
  # Number of contigs
  num_contigs <- length(seqs)
  
  # Total number of bases
  total_bases <- sum(width(seqs))
  
  # GC content (overall)
  alph_freqs <- colSums(alphabetFrequency(seqs, baseOnly = TRUE))
  gc_content <- (alph_freqs["G"] + alph_freqs["C"]) / sum(alph_freqs) * 100
  
  # Return a simple list
  list(
    num_contigs = num_contigs,
    total_bases = total_bases,
    gc_content  = round(gc_content, 2)
  )
}

# 3. Attach these FASTA metrics to assembly_info
assembly_info_with_metrics <- assembly_info %>%
  rowwise() %>%
  mutate(
    metrics = list(summarize_fasta(full_path))
  ) %>%
  tidyr::unnest_wider(col = metrics) %>%
  ungroup()

# (D) MERGE FASTA METRICS WITH METADATA ---------------------------------------

assembly_info_update1 <- assembly_info_with_metrics %>%
  left_join(metadata, by = c("Isolate_ID" = "isolate_ID"))

# Optionally remove the 'full_path' if you no longer need it
assembly_info_update1 <- assembly_info_update1 %>%
  select(-full_path)

# (E) CHECK FOR MISSING METADATA ---------------------------------------------

missing_metadata <- assembly_info_update1 %>% filter(is.na(Participant_id))
if (nrow(missing_metadata) > 0) {
  warning(paste("There are", nrow(missing_metadata), "assemblies without metadata."))
} else {
  message("All assemblies have corresponding metadata.")
}

# (F) FINISHED! ---------------------------------------------------------------
# assembly_info_update1 now contains:
#  - file_name
#  - Isolate_ID
#  - assembler
#  - num_contigs, total_bases, gc_content
#  - all your merged metadata columns from batch1 + batch2
#
# You can continue exploring or saving the result:

# View(assembly_info_update1)
# write.csv(assembly_info_update1, "assembly_info_update1.csv", row.names = FALSE)


# ── a) rows where the join failed ────
missing_meta <- assembly_info_update1 %>% 
  filter(is.na(Participant_id))        # or whatever your participant‑ID col is

# ── b) quick report ────
if (nrow(missing_meta) == 0) {
  message("✅ Every FASTA matched a metadata record.")
} else {
  message(glue::glue(
    "⚠️  {nrow(missing_meta)} FASTA(s) have no metadata. ",
    "Here are the first few:"
  ))
  print(
    missing_meta %>% 
      select(file_name, Isolate_ID, assembler) %>% 
      head(10)
  )
}

# ── a) build the table of "comparable" participants ────
comparable_tbl <- assembly_info_update1 %>% 
  filter(!is.na(Timepoint), assembler %in% c("flye", "longcycler")) %>% 
  group_by(Participant_id, assembler) %>% 
  filter(n_distinct(Timepoint) >= 2) %>% 
  ungroup()

# ── b) counts per participant × assembler ────
count_tbl <- comparable_tbl %>% 
  count(Participant_id, assembler, name = "n_timepoints")

# participants that appear 3-to-5 times for a given assembler
selected_ids <- count_tbl %>% 
  filter(n_timepoints %in% 3:5) %>%          # keep 3, 4 or 5 TP
  distinct(Participant_id) %>%               # drop dupes if both assemblers qualify
  pull(Participant_id)

# quick sanity-check
length(selected_ids)      # how many participants?
head(selected_ids)        # show a few



# Who still needs at least one more time‑point in Flye?
# needs_more <- assembly_info_update1 %>% 
#  filter(assembler == "flye") %>% 
#  distinct(Participant_id, Timepoint) %>% 
#  add_count(Participant_id) %>% 
#  filter(n < 2)

# Quick overview of assemblies per participant
assembly_info_update1 %>% 
  count(Participant_id, assembler) %>% 
  tidyr::pivot_wider(
    names_from = assembler,
    values_from = n,
    values_fill = 0
  )


rm(list = ls())

#Set working directory: this is where my project lives. Current project is stored on a flashdrive. 
setwd("D:/Batch 1-6 Ecoli FULL sequence YELLOW Study")

#Library packages and dependencies
library(Biostrings)
library(dplyr)
library(stringr)
library(tidyr)
library(readr)
library(seqinr)

#Load file lists
isolate_files <- list.files(
  path = "raw_fasta",            #working from D:/Ecoli Sequence NIEUW Yellow Study: enter the raw_fasta file
  pattern = "\\.fasta$|\\.fa$",  #./ dubble backslash to get a slingel backslash dot as assigned for "any" regex patern $ for ends with
  full.names = TRUE              # to get "insert_name" + "extention (.fasta)"
)

head(isolate_files) # see first few names not neccesary

gene_files <- list.files(
  path = "reference_genes",             #working from D:/Ecoli Sequence NIEUW Yellow Study: enter the reference_genes file
  pattern = "\\.fasta$|\\.fa$",         #./ dubble backslash to get a slingel backslash dot as assigned for "any" regex patern $ for ends with
  full.names = TRUE                     #with file extension
)


extract_gene <- function(genome_file, gene_file) {   # make function with two variables (being genome_file and gene_file)
  
  #readDNAStringSet is used to read DNA sequences from a file into an XStringSet object
  genome <- readDNAStringSet(genome_file)
  gene_ref <- readDNAStringSet(gene_file)
  
  # combine genome (some isolates may have multiple contigs)
  genome_seq <- paste(as.character(genome), collapse = "")
  
  # get reference gene sequence
  gene_seq <- as.character(gene_ref[[1]])
  
  # find approximate match (pattern search)
  match_pos <- matchPattern(gene_seq, genome_seq, max.mismatch = 60) 
  # the max mismatch is a variable of how presice the function should operate the mismatch of 60 is based on the median of 600*10 for a 10% mismatch
  
  # if no match found fll with NA
  if (length(match_pos) == 0) {
    return(NA)
  }
  
  # retrieve start and end number
  start <- start(match_pos)[1]
  end <- end(match_pos)[1]
  
  # paste info into one return
  extracted <- substr(genome_seq, start, end)
  
  return(extracted)
}

#make empty df because we are going to be pasting into this
results <- data.frame()

# obj counting as to see the number go up to get a rough idea of how long this is going to take (long) 
obj <- 1

# for every ref gene in reference genes folder
for (gene_file in gene_files) {
  
  #get gene name
  gene_name <- tools::file_path_sans_ext(basename(gene_file))
  
  # for every isolate in raw_fasta folder
  for (iso_file in isolate_files) {
    
    # counter goes brrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrr
    obj <- obj+1
    writeLines(paste(obj))
    
    # get isolate name per input of for loop
    iso_name <- tools::file_path_sans_ext(basename(iso_file))
    
    # use to function explained above and basically regex the ref (first for loop) to the raw_fasta file (second for loop)
    seq <- extract_gene(iso_file, gene_file)
    
    # stick the data from seq to the results 
    results <- rbind(results, data.frame(
      isolate = iso_name,
      gene = gene_name,
      sequence = seq
    ))
  }
}

# the two for loops interact by running the extract_gene function first using ref1 against all raw_fasta files
# here after ref2 is extract_gene() against all raw_fasta files
# etc

#write_rds(results, "sequence_extraction_output_ref_vs_fasta_all_batches.rds")
# read_rds()
View(results)


narG_seqs <- results %>% filter(gene == "NarG") %>% filter(!is.na(sequence))

seqs <- DNAStringSet(narG_seqs$sequence)
names(seqs) <- narG_seqs$isolate

results %>% 
  group_by(gene) %>%
  filter(!is.na(sequence)) %>%
  summarise(seq_n = n())

#------------------------------------------------------------#
# ALLELE TYPING
#------------------------------------------------------------#

#Remove  missing (NA) sequences from the table.
results_clean <- results %>%
  filter(!is.na(sequence))

#Allele calling function
assign_alleles <- function(df) {
  
  # df = one gene at a time
  
  unique_seqs <- unique(df$sequence)
  
  allele_map <- setNames(
    paste0("Allele", seq_along(unique_seqs)),
    unique_seqs
  )
  
  df$allele <- allele_map[df$sequence]
  
  return(df)
}

#Apply per gene
allele_results <- results_clean %>%
  group_by(gene) %>%
  group_modify(~ assign_alleles(.x)) %>%
  ungroup()

View(allele_results)

#write_rds(allele_results, "Allele_typing_output_allele_profile_per_gene_all_batches.rds")

#Build isolate profiles allele type # per isolate per gene
profile <- allele_results %>%
  select(isolate, gene, allele) %>%
  pivot_wider(names_from = gene, values_from = allele)

View(profile)

#write_rds(profile, "Allele_typing_output_allele_profile_per_isolate_all_batches.rds")

#Build COMPACT isolate profiles allele type # per isolate per gene
profile$genotype <- apply(profile[,-1], 1, paste, collapse = "_")

#Count genotypes
genotype_counts <- profile %>%
  dplyr::count(genotype) %>%
  arrange(desc(n))

#Arranging all genotypes with isolates with this specific genotype.
genotype_isolate_type <- profile %>%
  group_by(genotype) %>%
  summarise(
    isolates = paste(isolate, collapse = ", "),
    n = n()
  ) %>%
  arrange(desc(n))

View(genotype_isolate_type)

#write_rds(genotype_isolate_type, "Genotype_isolate_combination_overview_all_batches.rds")

#------------------------------------------------------#
#METADATA ORGANIZING
#------------------------------------------------------#
install.packages("haven")
library(haven)

install.packages("readxl")
library(readxl)

#Study timepoints SPSS integrating into metadata

spss_data <- read_sav("../YELLOW dataset stage Rowena Studie meetmomenten long format _1.sav")

metadata <- read_excel(
  "D:/Batch 1-6 Ecoli FULL sequence YELLOW Study/metadata/OVERVIEW E.coli batch 1-6 - CLEAN_RC.xlsx",
  skip = 3
)

metadata <- metadata[!is.na(metadata$isolate_id), ]

#colnames(metadata) <- c("isolate_id","batch", "participant_id", "uti_sticker", "date", "species", "object",
#                      "id_isolate", "örganism", "beoordeling", "kiemgetal", "Archive", "timepoint", "nitrite_result","infection_group",
#                       "symptom_incontinence", "symptom_pus", "symptom_flankpijn", "symptom_suprapubic",
##                     "urine_collection_method", "symptom_dysuria", "symptom_aandrang", "symptom_frequency",
#                    "symptom_fever", "symptom_chills", "symptom_delirium", "symptom_other", "symptom_none")

spss_clean <- spss_data %>%
  select(participant_id, timepoint, urinestick_nitriet)

View(spss_clean)

str(metadata)
str(spss_clean)

last_column <- names(spss_clean)[ncol(spss_clean)]


#UTI cases SPSS integrating into metadata

spss_uti <- read_sav("YELLOW dataset stage Rowena NEW dipslide long format.sav")

spss_uti_clean <- spss_uti %>%
  select(participant_id, timepoint, urinestick_nitriet)



join_data <- bind_rows(
  spss_clean |>
    mutate(
      timepoint = paste0("T", as.integer(timepoint) - 1),
      participant_id = as.character(participant_id)
    ) |>
    select(participant_id, timepoint, urinestick_nitriet),
  
  spss_uti_clean |>
    mutate(
      timepoint = paste0("UTI-", as.integer(timepoint)),
      participant_id = as.character(participant_id)
    ) |>
    select(participant_id, timepoint, urinestick_nitriet)
)

metadata <- metadata |>
  left_join(
    join_data,
    by = c("participant_id", "timepoint")
  )


#Adding the allele profiles into metadata.

add_columns_to_metadata <- function(metadata, new_df, join_col_new, cols_to_add = "all") {
  if (identical(cols_to_add, "all")) {
    cols_to_add <- setdiff(names(new_df), join_col_new)
  }
  metadata |>
    left_join(
      new_df |>
        select(all_of(c(join_col_new, cols_to_add))),
      by = c("isolate_id" = join_col_new)
    )
}

metadata <- add_columns_to_metadata(
  metadata  = metadata,
  new_df    = profile,
  join_col_new = "isolate",
  cols_to_add  = "all"                #cols_to_add  = c("col1", "col2")
)

#write_rds(metadata, "full_overview_metadata_allele_profiles")

#---------------------------------------------------------------#
#ANALYSIS: PROFILE WITH METADATA
#---------------------------------------------------------------#

#---------------------------------------------------------------#
# FIRST-PASS ALLELE VS DIPSTICK NITRITE ANALYSIS
# Run this add-on after metadata has already been created.
# Uses existing allele columns in metadata; does not rerun extraction.
#---------------------------------------------------------------#


# These are the nitrate/nitrite pathway genes included in the downstream
# allele/nitrite analysis. They should match the gene columns added to metadata
# from the allele profile table above.
selected_nitrite_genes <- c(
  "CsgD", "NarG", "NarH", "NarI", "NarJ", "NarL",
  "NarP", "NarQ", "NarV", "NarX", "NarY", "NarZ"
)

# All output files from this add-on are written here.
output_dir <- file.path("outputs", "nitrate_blast")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Work on a copy of metadata so the original metadata object stays intact.
# If a selected gene has no allele column, add it as NA so later all_of()
# selections do not fail.
metadata_for_alleles <- metadata
for (gene_name in selected_nitrite_genes) {
  if (!gene_name %in% names(metadata_for_alleles)) {
    metadata_for_alleles[[gene_name]] <- NA_character_
  }
}

# Build the main analysis table.
# urinestick_nitriet is recoded into positive/negative/missing, and the table
# keeps only the isolate identifiers, nitrite result, and selected allele columns.
nitrite_allele_table <- metadata_for_alleles %>%
  mutate(
    nitrite_code = suppressWarnings(as.numeric(urinestick_nitriet)),
    nitrite_label = as.character(haven::as_factor(urinestick_nitriet)),
    nitrite_binary = case_when(
      nitrite_code == 1 | nitrite_label == "Positief" ~ "positive",
      nitrite_code == 2 | nitrite_label == "Negatief" ~ "negative",
      TRUE ~ "missing"
    )
  ) %>%
  select(isolate_id, participant_id, timepoint, nitrite_binary, all_of(selected_nitrite_genes))

# Statistical tests use only rows with interpretable positive/negative nitrite.
analysis_rows <- nitrite_allele_table %>%
  filter(nitrite_binary %in% c("positive", "negative"))

# Descriptive genotype summary.
# Each genotype is a combined allele profile across the selected genes.
# This table shows how many isolates in each genotype group are nitrite
# positive, negative, or missing.
genotype_nitrite_summary <- metadata %>%
  mutate(
    nitrite_code = suppressWarnings(as.numeric(urinestick_nitriet)),
    nitrite_label = as.character(haven::as_factor(urinestick_nitriet)),
    nitrite_binary = case_when(
      nitrite_code == 1 | nitrite_label == "Positief" ~ "positive",
      nitrite_code == 2 | nitrite_label == "Negatief" ~ "negative",
      TRUE ~ "missing"
    )
  ) %>%
  group_by(genotype) %>%
  summarise(
    isolates_total = n(),
    nitrite_positive = sum(nitrite_binary == "positive", na.rm = TRUE),
    nitrite_negative = sum(nitrite_binary == "negative", na.rm = TRUE),
    nitrite_missing = sum(nitrite_binary == "missing", na.rm = TRUE),
    nitrite_positive_percent = if_else(
      nitrite_positive + nitrite_negative > 0,
      100 * nitrite_positive / (nitrite_positive + nitrite_negative),
      NA_real_
    ),
    .groups = "drop"
  ) %>%
  left_join(
    genotype_isolate_type %>%
      rename(genotype_group_n = n),
    by = "genotype"
  ) %>%
  arrange(desc(genotype_group_n), desc(nitrite_positive_percent))


# Main allele-distribution Fisher test.
# For each gene, this compares the full allele distribution
# (Allele1, Allele2, Allele3, absent, etc.) against nitrite positive/negative.
# simulate.p.value is used because these multi-allele contingency tables can be
# too large/sparse for exact Fisher enumeration.
allele_fisher <- function(gene_name) {
  test_data <- analysis_rows %>%
    mutate(allele = replace_na(as.character(.data[[gene_name]]), "absent"))
  
  fisher_table <- table(test_data$allele, test_data$nitrite_binary)
  
  p_value <- if (nrow(fisher_table) > 1 && ncol(fisher_table) > 1) {
    fisher.test(fisher_table, simulate.p.value = TRUE, B = 10000)$p.value
  } else {
    NA_real_
  }
  
  tibble(
    gene = gene_name,
    allele_categories = nrow(fisher_table),
    p_value = p_value
  )
}

# Apply the gene-level allele-distribution Fisher test to all selected genes,
# then use Benjamini-Hochberg correction because multiple genes are tested.
allele_vs_nitrite_fisher_bh <- bind_rows(lapply(selected_nitrite_genes, allele_fisher)) %>%
  mutate(
    bh_q_value = p.adjust(p_value, method = "BH"),
    significant_bh_0_05 = bh_q_value < 0.05
  )

# Descriptive allele counts by gene and nitrite status.
# This helps explain which allele categories drive the Fisher test.
allele_vs_nitrite_counts <- analysis_rows %>%
  select(nitrite_binary, all_of(selected_nitrite_genes)) %>%
  pivot_longer(
    cols = all_of(selected_nitrite_genes),
    names_to = "gene",
    values_to = "allele"
  ) %>%
  mutate(allele = replace_na(as.character(allele), "absent")) %>%
  dplyr::count(gene, allele, nitrite_binary, name = "n")

# QC-only presence/absence summary.
# This collapses allele categories into gene present vs absent and should not
# replace the allele-distribution analysis above.
presence_absence_qc_counts <- analysis_rows %>%
  select(nitrite_binary, all_of(selected_nitrite_genes)) %>%
  pivot_longer(
    cols = all_of(selected_nitrite_genes),
    names_to = "gene",
    values_to = "allele"
  ) %>%
  mutate(gene_present = !is.na(allele)) %>%
  dplyr::count(gene, gene_present, nitrite_binary, name = "n")

# Pull out any BH-significant gene-level allele-distribution hits for the text
# summary. If none exist, the summary says no BH q-value was below 0.05.
bh_hits <- allele_vs_nitrite_fisher_bh %>%
  filter(!is.na(bh_q_value), bh_q_value < 0.05)

# Short plain-English summary file for quick interpretation.
summary_lines <- c(
  "Allele distribution vs dipstick nitrite first-pass analysis",
  paste0("Rows with positive/negative nitrite results: ", nrow(analysis_rows)),
  "Fisher exact tests compare each gene's allele distribution with nitrite positive/negative status.",
  "BH correction was used because multiple genes were tested.",
  "Presence/absence output is QC only; allele distribution is the main analysis.",
  "Results are exploratory because participants can contribute repeated timepoints.",
  if (nrow(bh_hits) == 0) {
    "No BH q-value was below 0.05."
  } else {
    paste0("Genes with BH q-value below 0.05: ", paste(bh_hits$gene, collapse = ", "))
  }
)

# Long-format allele table: one row per isolate-gene combination.
# This makes allele counting and allele-specific tests easier.
allele_long <- analysis_rows %>%
  select(nitrite_binary, all_of(selected_nitrite_genes)) %>%
  pivot_longer(
    cols = all_of(selected_nitrite_genes),
    names_to = "gene",
    values_to = "allele"
  ) %>%
  mutate(allele = replace_na(as.character(allele), "absent"))

# Count how often each allele appears for each gene.
allele_counts_by_gene <- allele_long %>%
  dplyr::count(gene, allele, name = "n") %>%
  arrange(gene, desc(n), allele)

# Only test allele-specific associations for alleles seen more than 10 times.
# This avoids unstable tests on very rare alleles.
alleles_to_test <- allele_counts_by_gene %>%
  dplyr::filter(n > 10, allele != "absent")

# Secondary exploratory allele-specific Fisher test.
# For each common allele, test that allele vs all other alleles/absent for the
# same gene, against nitrite positive/negative.
allele_specific_fisher <- function(gene_name, allele_name, allele_n) {
  test_data <- allele_long %>%
    dplyr::filter(gene == gene_name) %>%
    dplyr::mutate(target_allele = allele == allele_name)
  allele_positive <- sum(test_data$target_allele & test_data$nitrite_binary == "positive")
  allele_negative <- sum(test_data$target_allele & test_data$nitrite_binary == "negative")
  other_positive <- sum(!test_data$target_allele & test_data$nitrite_binary == "positive")
  other_negative <- sum(!test_data$target_allele & test_data$nitrite_binary == "negative")
  fisher_table <- matrix(
    c(allele_positive, allele_negative, other_positive, other_negative),
    nrow = 2,
    byrow = TRUE
  )
  p_value <- if (all(rowSums(fisher_table) > 0)) {
    fisher.test(fisher_table)$p.value
  } else {
    NA_real_
  }
  tibble(
    gene = gene_name,
    allele = allele_name,
    allele_total = allele_n,  
    allele_positive = allele_positive,       # this allele + nitrite positive
    allele_negative = allele_negative,       # this allele + nitrite negative
    other_positive = other_positive,         # other alleles/absent + nitrite positive
    other_negative = other_negative,         # other alleles/absent + nitrite negative
    p_value = p_value
  )
}

# Run the allele-specific tests and apply BH correction across all tested alleles.
allele_specific_vs_nitrite_fisher_bh <- bind_rows(lapply(
  seq_len(nrow(alleles_to_test)),
  function(i) {
    allele_pair <- alleles_to_test[i, ]
    allele_specific_fisher(allele_pair$gene, allele_pair$allele, allele_pair$n)
  }
)) %>%
  dplyr::mutate(
    bh_q_value = p.adjust(p_value, method = "BH"),
    significant_bh_0_05 = bh_q_value < 0.05
  )

#---------------------------------------------------------------#
# SAVE OUTPUT FILES
#---------------------------------------------------------------#

write_csv(
  genotype_nitrite_summary,
  file.path(output_dir, "genotype_nitrite_summary.csv")
)

write_csv(
  allele_vs_nitrite_fisher_bh,
  file.path(output_dir, "allele_vs_nitrite_fisher_bh.csv")
)

write_csv(
  allele_vs_nitrite_counts,
  file.path(output_dir, "allele_vs_nitrite_counts.csv")
)

write_csv(
  presence_absence_qc_counts,
  file.path(output_dir, "presence_absence_qc_counts.csv")
)

write_csv(
  nitrite_allele_table,
  file.path(output_dir, "nitrite_allele_table.csv")
)

write_csv(
  allele_specific_vs_nitrite_fisher_bh,
  file.path(output_dir, "allele_specific_vs_nitrite_fisher_bh.csv")
)

write_csv(
  allele_counts_by_gene,
  file.path(output_dir, "allele_counts_by_gene.csv")
)

writeLines(
  summary_lines,
  file.path(output_dir, "nitrite_allele_association_summary.txt")
)

View(allele_specific_vs_nitrite_fisher_bh)



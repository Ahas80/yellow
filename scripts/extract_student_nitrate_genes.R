#!/usr/bin/env Rscript
# ==============================================================================
# extract_student_nitrate_genes.R
# ==============================================================================
#
# GOAL:
#   Extract the exact presence/absence data for all E. coli nitrate, nitrite, 
#   regulatory, and molybdenum cofactor genes from the massive Panaroo pangenome 
#   matrix. This creates a small, clean CSV file specifically tailored for the 
#   urinary dipstick / nitrate reduction student project.
#
# INPUT:
#   - results/wgs/pan/gene_presence_absence.csv
#
# OUTPUT:
#   - results/student_nitrate_pathway_genes.csv
# ==============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
})

# 1. Paths
panaroo_file <- "results/wgs/pan/gene_presence_absence.csv"
out_file <- "results/student_nitrate_pathway_genes.csv"

if (!file.exists(panaroo_file)) {
  stop("Cannot find Panaroo output at: ", panaroo_file)
}

message("Loading Panaroo pangenome matrix...")
df <- read_csv(panaroo_file, show_col_types = FALSE)

# 2. Define the exact genes we are hunting for
target_genes <- c(
  # Nitrate Reductase A
  "narG", "narH", "narJ", "narI", 
  # Nitrate Reductase Z
  "narZ", "narY", "narW", "narV",
  # Periplasmic Nitrate Reductase
  "napF", "napD", "napA", "napG", "napH", "napB", "napC",
  # Regulators & Receptors
  "narX", "narQ", "narL", "narP", "fnr",
  # Transporters
  "narK", "narU",
  # Nitrite Reductases
  "nirB", "nirD", "nrfA", "nrfB", "nrfC", "nrfD", "nrfE", "nrfF", "nrfG",
  # Molybdenum Cofactor Synthesis (Crucial for nitrate reduction)
  "moaA", "moaB", "moaC", "moaD", "moaE",
  "mobA", "mobB",
  "modA", "modB", "modC", "modE",
  "moeA", "moeB"
)

# 3. Filter the matrix
# We use a regex to exactly match the gene名前 or catch paralogs (e.g. narG_1, narG_2)
regex_target <- paste0("^(", paste(target_genes, collapse = "|"), ")(_.*)?$")

message("Filtering for nitrate/nitrite/molybdenum pathway genes...")
df_student <- df %>%
  filter(str_detect(Gene, regex_target))

# 4. Save and summarize
write_csv(df_student, out_file)

message("==================================================")
message("SUCCESS!")
message("Extracted ", nrow(df_student), " pathway genes from the pangenome.")
message("Saved clean dataset for the student to: ", out_file)
message("==================================================")

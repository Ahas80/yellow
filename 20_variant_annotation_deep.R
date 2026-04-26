#!/usr/bin/env Rscript
# ==============================================================================
# 20_variant_annotation_deep.R
# ==============================================================================
#
# GOAL:
#   Annotate within-host SNP variants with gene-level context by cross-
#   referencing variant positions against Prokka GFF annotations.  This adds
#   biological meaning to the raw SNP calls: is a variant in a known gene?
#   What is the gene product?  Is it intergenic?
#
# WHY THIS SCRIPT EXISTS:
#   Script 18_annotate_variants.R detects SNPs between consecutive timepoints
#   for the same participant, but only reports positions (contig + coordinate).
#   This script overlays those positions on the GFF annotation to determine
#   which genes are affected.  For the phenotype-switch candidates (ASB→UTI),
#   this is critical: if a SNP hits a known virulence gene, it could explain
#   the clinical transition.
#
# CURRENT SCOPE:
#   Currently focused on two specific participants (40001, 40004) that were
#   identified as ASB→UTI transition candidates.  Can be extended to all
#   participants by modifying the 'candidates' vector.
#
# INPUTS:
#   - results/longitudinal/annotated_snps.csv  (from 18_annotate_variants.R)
#   - assembly_metadata.csv
#   - results/prokka_prefixed_slim/*.gff       (Prokka gene annotations)
#
# OUTPUTS:
#   - results/longitudinal/variant_annotation_detailed.csv
#
# KEY DESIGN DECISIONS:
#   - Focuses on curated phenotype-switch participants to prioritise manual,
#     biologically interpretable annotation during thesis-stage deep dive.
#
# POSITION IN PIPELINE:
#   - Phase 3 deep annotation step after 18_annotate_variants.R.
#
# NOTES / LIMITATIONS:
#   - Current participant subset is intentionally narrow and not cohort-wide.
# ==============================================================================

suppressPackageStartupMessages({
    library(dplyr)
    library(readr)
    library(tidyr)
    library(stringr)
})

msg <- function(...) message(sprintf(...))

# 0. Load Config
source("00_config.R")

# Hardcode paths for debugging
# DIR_RESULTS is defined in 00_config.R

# 1. Load Data
# ------------------------------------------------------------------------------
snps_file <- file.path(DIR_RESULTS, "longitudinal", "annotated_snps.csv")
snps <- read_csv(snps_file, show_col_types = FALSE) %>%
    mutate(Participant_id = as.character(Participant_id))

meta <- read_csv("assembly_metadata.csv", show_col_types = FALSE) %>%
    mutate(Participant_id = as.character(Participant_id))

# Filter for the key candidates (ASB -> UTI)
candidates <- c("40001", "40004")
target_snps <- snps %>%
    filter(Participant_id %in% candidates) %>%
    filter(Type == "SNP")

msg("Target SNPs: %d", nrow(target_snps))

# 2. Helper: Parse GFF Manually
# ------------------------------------------------------------------------------
get_gff_data <- function(pid, tp) {
    iso_row <- meta %>% filter(Participant_id == pid, Timepoint == tp)
    if (nrow(iso_row) == 0) {
        return(NULL)
    }

    iso_id <- iso_row$file_name[1]
    # iso_id usually ends in .fasta. Remove it.
    id_clean <- sub("\\.fasta$", "", iso_id)

    # Search for GFF in DIR_PROKKA_SLIM (defined in 00_config.R)
    # Pattern: id_clean + anything + .gff
    gff_pattern <- paste0(id_clean, ".*\\.gff$")
    gff_files <- list.files(DIR_PROKKA_SLIM, pattern = gff_pattern, full.names = TRUE)

    # Fallback: Search in full Prokka dir (recursive)
    if (length(gff_files) == 0) {
        msg("Not found in slim dir. Searching full Prokka dir: %s", DIR_PROKKA)
        # Recursive search in DIR_PROKKA
        gff_files <- list.files(DIR_PROKKA, pattern = paste0(id_clean, ".*\\.gff$"), full.names = TRUE, recursive = TRUE)
    }

    if (length(gff_files) == 0) {
        msg("No GFF found for %s in %s or fallback", id_clean, DIR_PROKKA_SLIM)
        return(NULL)
    }

    gff_file <- gff_files[1]
    msg("Loading GFF: %s", basename(gff_file))

    # Manual GFF parsing
    gff <- read_tsv(gff_file, comment = "#", col_names = c("seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes"), show_col_types = FALSE)

    gff <- gff %>% filter(type %in% c("CDS", "gene"))

    extract_attr <- function(attr_col, key) {
        str_extract(attr_col, paste0("(?<=", key, "=)[^;]+"))
    }

    gff <- gff %>%
        mutate(
            gene = extract_attr(attributes, "gene"),
            product = extract_attr(attributes, "product"),
            locus_tag = extract_attr(attributes, "locus_tag")
        )

    return(gff)
}

# 3. Annotate
# ------------------------------------------------------------------------------
results <- list()

pairs <- target_snps %>%
    select(Participant_id, From_Time) %>%
    distinct()

msg("Processing %d pairs", nrow(pairs))

for (i in 1:nrow(pairs)) {
    p <- pairs[i, ]
    pid <- p$Participant_id
    t_ref <- p$From_Time

    msg("Processing Pair: %s (Ref: %s)", pid, t_ref)

    gff_df <- get_gff_data(pid, t_ref)
    if (is.null(gff_df)) next

    pair_snps <- target_snps %>%
        filter(Participant_id == pid, From_Time == t_ref)

    for (j in 1:nrow(pair_snps)) {
        pos <- pair_snps$Pos_Ref[j]

        hit <- gff_df %>%
            filter(start <= pos & end >= pos) %>%
            head(1)

        if (nrow(hit) > 0) {
            res_row <- pair_snps[j, ] %>%
                mutate(
                    Region = hit$type,
                    Gene = hit$gene,
                    Product = hit$product,
                    Locus = hit$locus_tag
                )
        } else {
            res_row <- pair_snps[j, ] %>%
                mutate(
                    Region = "Intergenic",
                    Gene = NA,
                    Product = NA,
                    Locus = NA
                )
        }
        results[[length(results) + 1]] <- res_row
    }
}

final_df <- bind_rows(results)

out_file <- file.path(DIR_RESULTS, "longitudinal", "variant_annotation_detailed.csv")
write_csv(final_df, out_file)
msg("Saved annotated variants to %s", out_file)

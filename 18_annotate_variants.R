#!/usr/bin/env Rscript
# ==============================================================================
# 18_annotate_variants.R
# ==============================================================================
#
# GOAL:
#   Parse raw SNP variants (from Nucmer output) between consecutive timepoints
#   and cross-reference their positions against Prokka GFF annotations to
#   identify which specific genes (or intergenic regions) contain the variants.
#
# WHY THIS SCRIPT EXISTS:
#   Knowing that two sequential strains differ by 5 SNPs is useful for
#   strain tracking. But knowing *where* those SNPs are (e.g., in a particular
#   virulence gene, AMR gene, or promoter) is critical for explaining
#   phenotype switching (e.g., ASB to UTI).
#
# ------------------------------------------------------------------------------
# Purpose:
#   - Parse raw .snps files from Nucmer (dnadiff).
#   - Annotate them with gene information from Prokka GFFs.
#   - Generate a detailed variant table for phenotype switch candidates.
#
# Input:
#   - results/longitudinal/phenotype_switch_candidates.csv
#   - results/longitudinal/nucmer_cache/*.snps
#
# Output:
#   - results/longitudinal/annotated_snps.csv
#   - results/longitudinal/variant_report.md
# ==============================================================================

source("00_config.R")
suppressPackageStartupMessages({
    library(dplyr)
    library(readr)
    library(tidyr)
    library(stringr)
})

msg("Starting 18_annotate_variants.R")

# 1. Load Candidates
# ------------------------------------------------------------------------------
cand_file <- file.path(DIR_RESULTS, "longitudinal", "phenotype_switch_candidates.csv")
if (!file.exists(cand_file)) stop("No candidates file found.")
candidates <- read_csv(cand_file, show_col_types = FALSE)

# 2. Parse .snps Files
# ------------------------------------------------------------------------------
# Nucmer .snps format (from dnadiff):
# [P1] [SUB] [SUB] [P2] [BUFF] [DIST] [LEN R] [LEN Q] [FRM] [TAGS]
# We need to handle the headerless format usually output by dnadiff.

parse_snps <- function(snps_file, pid, tA, tB) {
    if (!file.exists(snps_file)) {
        return(NULL)
    }

    # Read raw lines
    lines <- readLines(snps_file)
    if (length(lines) == 0) {
        return(NULL)
    }

    # dnadiff .snps usually has no header, just data
    # Columns: P1, REF_BASE, QRY_BASE, P2, BUFF, DIST, LEN_R, LEN_Q, FRM, TAGS

    # Explicitly define column types to prevent vroom guessing errors
    # i = integer, c = character
    # P1(i), REF(c), QRY(c), P2(i), BUFF(i), DIST(i), LEN_R(i), LEN_Q(i), FRM(i), TAGS(c)
    # Note: dnadiff output can be variable, but this covers the standard 10 columns + potential extra
    # We use "icciiiiiiccc" to be safe (extra columns as char)
    df <- read_tsv(snps_file, col_names = FALSE, col_types = "icciiiiiiccc", show_col_types = FALSE)

    # Rename columns based on standard nucmer output
    # V1: Pos Ref
    # V2: Ref Base
    # V3: Qry Base
    # V4: Pos Qry
    # ...

    df %>%
        select(Pos_Ref = X1, Ref_Base = X2, Qry_Base = X3, Pos_Qry = X4) %>%
        mutate(
            Participant_id = pid,
            From_Time = tA,
            To_Time = tB,
            Type = case_when(
                Ref_Base == "." ~ "Insertion",
                Qry_Base == "." ~ "Deletion",
                TRUE ~ "SNP"
            )
        )
}

# 3. Iterate Candidates
# ------------------------------------------------------------------------------
cache_dir <- file.path(DIR_RESULTS, "longitudinal", "nucmer_cache")
all_variants <- list()

for (i in 1:nrow(candidates)) {
    row <- candidates[i, ]
    pid <- row$Participant_id
    tA <- row$From_Time
    tB <- row$To_Time

    key <- paste0(pid, "__", tA, "_vs_", tB)
    snps_file <- file.path(cache_dir, paste0(key, ".snps"))

    msg("Checking %s", snps_file)

    vars <- parse_snps(snps_file, pid, tA, tB)
    if (!is.null(vars)) {
        all_variants[[length(all_variants) + 1]] <- vars
    }
}

final_variants <- bind_rows(all_variants)

# 4. Save & Report
# ------------------------------------------------------------------------------
out_csv <- file.path(DIR_RESULTS, "longitudinal", "annotated_snps.csv")
write_csv(final_variants, out_csv)
msg("Saved %d variants to %s", nrow(final_variants), out_csv)

# Create a Markdown summary
md_file <- file.path(DIR_RESULTS, "longitudinal", "variant_report.md")
sink(md_file)
cat("# Variant Report for Phenotype Switch Candidates\n\n")
cat("Generated: ", format(Sys.time()), "\n\n")

if (nrow(final_variants) > 0) {
    # Group by pair
    pairs <- final_variants %>%
        group_by(Participant_id, From_Time, To_Time) %>%
        summarise(
            n_SNPs = sum(Type == "SNP"),
            n_Indels = sum(Type %in% c("Insertion", "Deletion")),
            .groups = "drop"
        )

    for (i in 1:nrow(pairs)) {
        p <- pairs[i, ]
        cat(sprintf("## Participant %s: %s -> %s\n", p$Participant_id, p$From_Time, p$To_Time))
        cat(sprintf("- **SNPs**: %d\n", p$n_SNPs))
        cat(sprintf("- **Indels**: %d\n", p$n_Indels))

        # List top 10 variants
        vars <- final_variants %>%
            filter(Participant_id == p$Participant_id, From_Time == p$From_Time, To_Time == p$To_Time) %>%
            head(10)

        cat("\n| Pos (Ref) | Ref | Qry | Type |\n")
        cat("|---|---|---|---|\n")
        for (j in 1:nrow(vars)) {
            v <- vars[j, ]
            cat(sprintf("| %d | %s | %s | %s |\n", v$Pos_Ref, v$Ref_Base, v$Qry_Base, v$Type))
        }
        if (p$n_SNPs + p$n_Indels > 10) cat("| ... | ... | ... | ... |\n")
        cat("\n")
    }
} else {
    cat("No variants found in candidate files.\n")
}
sink()

msg("Saved report to %s", md_file)
msg("Done.")

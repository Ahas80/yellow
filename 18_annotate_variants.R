#!/usr/bin/env Rscript
# ==============================================================================
# 18_annotate_variants.R
# ==============================================================================
#
# GOAL:
#   Parse raw SNP variants (from Nucmer output) between consecutive timepoints
#   for phenotype-switch candidates.  Script 20 adds Prokka GFF gene context.
#
# WHY THIS SCRIPT EXISTS:
#   Knowing that two sequential strains differ by 5 SNPs is useful for strain
#   tracking. The raw position table produced here is the input to the deeper
#   gene-level annotation in 20_variant_annotation_deep.R.
#
# ------------------------------------------------------------------------------
# Purpose:
#   - Parse raw .snps files from Nucmer (dnadiff).
#   - Generate a variant position table for phenotype switch candidates.
#
# Input:
#   - results/longitudinal/phenotype_switch_candidates.csv
#   - results/longitudinal/evolution_events.csv (exact SHA-bound .snps paths)
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
candidates <- candidates %>%
    mutate(
        Participant_id = as.character(Participant_id),
        From_Time = normalise_timepoint_preserve_events(From_Time),
        To_Time = normalise_timepoint_preserve_events(To_Time)
    )

# 2. Parse .snps Files
# ------------------------------------------------------------------------------
# Nucmer .snps format (from dnadiff):
# [P1] [SUB] [SUB] [P2] [BUFF] [DIST] [LEN R] [LEN Q] [FRM] [TAGS]
# We need to handle the headerless format usually output by dnadiff.

parse_snps <- function(snps_file, expected_sha256, pid, tA, tB,
                       reference_fasta, reference_sha256,
                       query_fasta, query_sha256) {
    if (!file.exists(snps_file)) stop("Recorded SNP file is missing: ", snps_file)
    observed_sha256 <- unname(digest::digest(snps_file, algo = "sha256", file = TRUE))
    if (is.na(expected_sha256) || !identical(observed_sha256, expected_sha256)) {
        stop("SNP file content hash does not match evolution_events.csv: ", snps_file)
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
        select(Pos_Ref = X1, Ref_Base = X2, Qry_Base = X3, Pos_Qry = X4,
               Ref_Seqid = X11, Qry_Seqid = X12) %>%
        mutate(
            Participant_id = pid,
            From_Time = tA,
            To_Time = tB,
            SNP_Path = normalizePath(snps_file, winslash = "/", mustWork = TRUE),
            SNP_SHA256 = observed_sha256,
            Reference_FASTA_Path = reference_fasta,
            Reference_FASTA_SHA256 = reference_sha256,
            Query_FASTA_Path = query_fasta,
            Query_FASTA_SHA256 = query_sha256,
            Type = case_when(
                Ref_Base == "." ~ "Insertion",
                Qry_Base == "." ~ "Deletion",
                TRUE ~ "SNP"
            )
        )
}

# 3. Iterate Candidates
# ------------------------------------------------------------------------------
events_file <- file.path(DIR_RESULTS, "longitudinal", "evolution_events.csv")
if (!file.exists(events_file)) stop("Missing exact SNP provenance table: ", events_file, ". Run 16_within_host_evolution.R first.")
events <- read_csv(events_file, show_col_types = FALSE) %>%
    mutate(
        Participant_id = as.character(Participant_id),
        From_Time = normalise_timepoint_preserve_events(From_Time),
        To_Time = normalise_timepoint_preserve_events(To_Time)
    )
required_provenance <- c(
    "SNP_Path", "SNP_SHA256", "Fasta_Path_A", "Fasta_SHA256_A",
    "Fasta_Path_B", "Fasta_SHA256_B"
)
missing_provenance <- setdiff(required_provenance, names(events))
if (length(missing_provenance)) stop("evolution_events.csv lacks SNP provenance: ", paste(missing_provenance, collapse = ", "))
event_dupes <- events %>% count(Participant_id, From_Time, To_Time, name = "n") %>% filter(n != 1L)
if (nrow(event_dupes)) stop("evolution_events.csv must contain exactly one provenance row per candidate pair.")

pair_manifest <- candidates %>%
    select(Participant_id, From_Time, To_Time) %>%
    left_join(
        events %>% select(Participant_id, From_Time, To_Time, all_of(required_provenance)),
        by = c("Participant_id", "From_Time", "To_Time"), relationship = "one-to-one"
    )
if (any(is.na(pair_manifest$SNP_Path) | is.na(pair_manifest$SNP_SHA256))) {
    stop("One or more phenotype-switch pairs lack an exact SNP path/hash. Rerun script 16 and inspect dnadiff failures.")
}

analysis_manifest <- load_analysis_assemblies(FILE_ANALYSIS_ASSEMBLY_MANIFEST, require_files = TRUE) %>%
    mutate(
        full_path = normalizePath(full_path, winslash = "/", mustWork = TRUE),
        fasta_sha256 = vapply(full_path, digest::digest, character(1), algo = "sha256", file = TRUE)
    )
for (side in c("A", "B")) {
    path_col <- paste0("Fasta_Path_", side)
    sha_col <- paste0("Fasta_SHA256_", side)
    selected_key <- paste(pair_manifest[[path_col]], pair_manifest[[sha_col]], sep = "\n")
    current_key <- paste(analysis_manifest$full_path, analysis_manifest$fasta_sha256, sep = "\n")
    if (any(!selected_key %in% current_key)) stop("SNP provenance contains a reference/query FASTA not in the current Longcycler analysis manifest.")
}

all_variants <- list()

for (i in seq_len(nrow(pair_manifest))) {
    row <- pair_manifest[i, ]
    pid <- row$Participant_id
    tA <- row$From_Time
    tB <- row$To_Time

    snps_file <- row$SNP_Path

    msg("Checking %s", snps_file)

    vars <- parse_snps(
        snps_file, row$SNP_SHA256, pid, tA, tB,
        row$Fasta_Path_A, row$Fasta_SHA256_A,
        row$Fasta_Path_B, row$Fasta_SHA256_B
    )
    if (!is.null(vars)) {
        all_variants[[length(all_variants) + 1]] <- vars
    }
}

empty_variants <- tibble::tibble(
    Pos_Ref = integer(),
    Ref_Base = character(),
    Qry_Base = character(),
    Pos_Qry = integer(),
    Ref_Seqid = character(),
    Qry_Seqid = character(),
    Participant_id = character(),
    From_Time = character(),
    To_Time = character(),
    SNP_Path = character(),
    SNP_SHA256 = character(),
    Reference_FASTA_Path = character(),
    Reference_FASTA_SHA256 = character(),
    Query_FASTA_Path = character(),
    Query_FASTA_SHA256 = character(),
    Type = character()
)

final_variants <- if (length(all_variants) > 0) bind_rows(all_variants) else empty_variants

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

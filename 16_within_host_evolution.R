#!/usr/bin/env Rscript
# ==============================================================================
# 16_within_host_evolution.R
# ==============================================================================
#
# GOAL:
#   Deep-dive into phenotype-switch pairs (same strain, status changed).
#   For each pair: compute exact SNP distance and identify VF/plasmid gene
#   content changes.  This answers: "When the same E. coli lineage causes
#   ASB and then UTI, what genomic changes accompanied the transition?"
#
# NOTE:
#   Script 24_vf_longitudinal_dynamics.R complements this by computing VF
#   stability across ALL transitions (not just phenotype switches).
#
# ------------------------------------------------------------------------------
# Role: [Analysis] - Deep dive into "Phenotype Switch" pairs.
#
# Inputs:
#   - results/longitudinal/phenotype_switch_candidates.csv
#   - assembly_metadata.csv
#   - results/vf/vf_pa_all.csv
#   - results/plasmids/plasmidfinder_presence_absence.csv
#
# Outputs:
#   - results/longitudinal/evolution_events.csv
#   - results/longitudinal/evolution_summary.txt
#
# Purpose:
#   - For each Same-Strain pair that switched status (e.g. ASB -> UTI):
#     1. Calculate precise SNP distance (nucmer).
#     2. Identify gene content changes (VF gain/loss, Plasmid gain/loss).
#     3. Report "Chameleon" candidates: 0 SNPs but status change?
#
# KEY DESIGN DECISIONS:
#   - Restricts deep-dive to phenotype-switch candidates to prioritise biologically
#     meaningful transitions instead of all pairwise comparisons.
#
# POSITION IN PIPELINE:
#   - Phase 3 mechanistic follow-up after candidate generation in 15_.
#
# NOTES / LIMITATIONS:
#   - Candidate count is typically small; findings are hypothesis-generating.
# ==============================================================================

source("00_config.R")
source("11_compare_strains_helpers.R") # Reuse helpers

suppressPackageStartupMessages({
    library(dplyr)
    library(readr)
    library(tidyr)
    library(stringr)
})

msg("Starting 16_within_host_evolution.R")

# 1. Load Candidates
# ------------------------------------------------------------------------------
cand_file <- file.path(DIR_RESULTS, "longitudinal", "phenotype_switch_candidates.csv")
if (!file.exists(cand_file)) stop("No candidates file found. Run 15_longitudinal_patterns.R first.")
candidates <- read_csv(cand_file, show_col_types = FALSE)

msg("Loaded %d candidate pairs", nrow(candidates))

if (nrow(candidates) == 0) {
    msg("No candidates to analyze. Exiting.")
    quit(save = "no")
}

# 2. Load Genomic Data
# ------------------------------------------------------------------------------
# Assemblies
assemblies <- load_core_tables()$assemblies

# VF Data
vf_pa <- load_core_tables()$vf_pa

# Plasmid Data
inc_pa <- load_core_tables()$inc_pa

# 3. Analyze Each Pair
# ------------------------------------------------------------------------------
analyze_pair <- function(i) {
    row <- candidates[i, ]
    pid <- row$Participant_id
    tA <- row$From_Time
    tB <- row$To_Time

    msg("Analyzing %s: %s (%s) -> %s (%s)", pid, tA, row$From_Status, tB, row$To_Status)

    # Resolve samples
    sampA <- resolve_sample(pid, tA, assemblies)
    sampB <- resolve_sample(pid, tB, assemblies)

    if (is.null(sampA) || is.null(sampB)) {
        warning("Missing assembly for ", pid, " ", tA, " or ", tB)
        return(NULL)
    }

    # A. SNP Distance (re-run dnadiff to be sure or get details)
    # We use a temp cache for this specific analysis
    cache_dir <- file.path(DIR_RESULTS, "longitudinal", "nucmer_cache")
    key <- paste0(pid, "__", tA, "_vs_", tB)

    dna_res <- run_dnadiff(sampA$full_path, sampB$full_path, cache_dir, key)

    # B. VF Changes
    get_genes <- function(pa_df, p, t) {
        # pa_df has SampleKey or we filter by pid/tp
        # vf_pa from load_core_tables has SampleKey
        k <- make_sample_key(p, t)
        if (!k %in% pa_df$SampleKey) {
            return(character(0))
        }

        # Get columns that are 1
        row_vals <- pa_df %>%
            filter(SampleKey == k) %>%
            select(-any_of(c("Participant_id", "Timepoint", "tp_lab", "SampleKey")))
        names(row_vals)[which(row_vals[1, ] > 0)]
    }

    vf_A <- get_genes(vf_pa, pid, tA)
    vf_B <- get_genes(vf_pa, pid, tB)

    vf_gain <- setdiff(vf_B, vf_A)
    vf_loss <- setdiff(vf_A, vf_B)

    # C. Plasmid Changes
    # inc_pa has Isolate_ID, need to map via assemblies
    get_plasmids <- function(pa_df, iso_id) {
        if (is.null(pa_df)) {
            return(character(0))
        }
        if (!iso_id %in% pa_df$Isolate_ID) {
            return(character(0))
        }

        row_vals <- pa_df %>%
            filter(Isolate_ID == iso_id) %>%
            select(-Isolate_ID)
        names(row_vals)[which(row_vals[1, ] > 0)]
    }

    inc_A <- get_plasmids(inc_pa, sampA$Isolate_ID)
    inc_B <- get_plasmids(inc_pa, sampB$Isolate_ID)

    inc_gain <- setdiff(inc_B, inc_A)
    inc_loss <- setdiff(inc_A, inc_B)

    # D. Construct Result
    tibble(
        Participant_id = pid,
        From_Time = tA,
        To_Time = tB,
        From_Status = row$From_Status,
        To_Status = row$To_Status,
        Strain_ID = row$Strain_ID,
        SNPs = dna_res$TotalSNPs,
        AvgIdentity = dna_res$AvgIdentity,
        VF_Gain_Count = length(vf_gain),
        VF_Loss_Count = length(vf_loss),
        VF_Gained = paste(vf_gain, collapse = ";"),
        VF_Lost = paste(vf_loss, collapse = ";"),
        Plasmid_Gain_Count = length(inc_gain),
        Plasmid_Loss_Count = length(inc_loss),
        Plasmid_Gained = paste(inc_gain, collapse = ";"),
        Plasmid_Lost = paste(inc_loss, collapse = ";")
    )
}

results <- lapply(1:nrow(candidates), analyze_pair) %>% bind_rows()

# 4. Save & Summarize
# ------------------------------------------------------------------------------
out_file <- file.path(DIR_RESULTS, "longitudinal", "evolution_events.csv")
write_csv(results, out_file)
msg("Saved evolution events to %s", out_file)

# Summary Text
summary_file <- file.path(DIR_RESULTS, "longitudinal", "evolution_summary.txt")
sink(summary_file)
cat("=== Within-Host Evolution Summary ===\n")
cat("Generated: ", format(Sys.time()), "\n\n")

cat("Total Pairs Analyzed: ", nrow(results), "\n")
cat("Pairs with 0 SNPs: ", sum(results$SNPs == 0, na.rm = TRUE), "\n")
cat("Pairs with VF Changes: ", sum(results$VF_Gain_Count > 0 | results$VF_Loss_Count > 0), "\n")
cat("Pairs with Plasmid Changes: ", sum(results$Plasmid_Gain_Count > 0 | results$Plasmid_Loss_Count > 0), "\n\n")

cat("--- Top Candidates (ASB -> UTI) ---\n")
asb_uti <- results %>% filter(From_Status == "ASB", To_Status == "UTI")
if (nrow(asb_uti) > 0) {
    print(asb_uti)
} else {
    cat("None found.\n")
}

sink()
msg("Saved summary to %s", summary_file)
msg("Done.")

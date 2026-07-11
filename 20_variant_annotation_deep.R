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
#   which genes are affected.  For phenotype-switch candidates (Not_UTI -> UTI),
#   this is critical: if a SNP hits a known virulence gene, it could explain
#   the clinical transition.
#
# CURRENT SCOPE:
#   Uses current Not_UTI -> UTI transition tables when available, falling back
#   to script 15 phenotype-switch candidates.  This avoids hard-coded residents
#   and keeps the detailed SNP annotation aligned with the primary status logic.
#
# INPUTS:
#   - results/longitudinal/annotated_snps.csv  (variant positions from 18_)
#   - results/vf/vf_transition_case_index.csv or phenotype switch candidates
#   - assembly_metadata.csv
#   - results/prokka_prefixed_slim/*.gff       (Prokka gene annotations)
#
# OUTPUTS:
#   - results/longitudinal/variant_annotation_detailed.csv
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
    mutate(
        Participant_id = as.character(Participant_id),
        From_Time = normalise_timepoint_preserve_events(From_Time),
        To_Time = normalise_timepoint_preserve_events(To_Time)
    )

canonical_selection_file <- file.path(DIR_RESULTS, "qc", "canonical_assembly_selection.csv")
metadata_file <- if (file.exists(canonical_selection_file)) {
    canonical_selection_file
} else if (file.exists(FILE_METADATA)) {
    FILE_METADATA
} else {
    "assembly_metadata.csv"
}
meta <- read_csv(metadata_file, show_col_types = FALSE) %>%
    mutate(
        Participant_id = as.character(Participant_id),
        tp_lab = if ("tp_lab" %in% names(.)) normalise_timepoint_preserve_events(tp_lab) else normalise_timepoint_preserve_events(Timepoint),
        selected_canonical = if ("selected_canonical" %in% names(.)) as_pipeline_bool(selected_canonical) else FALSE,
        assembler = if ("assembler" %in% names(.)) as.character(assembler) else NA_character_,
        total_bases = if ("total_bases" %in% names(.)) suppressWarnings(as.numeric(total_bases)) else NA_real_
    )

load_transition_candidates <- function() {
    candidate_sources <- list()

    case_index_file <- file.path(DIR_VF, "vf_transition_case_index.csv")
    if (file.exists(case_index_file)) {
        candidate_sources[["vf_transition_case_index"]] <- read_csv(case_index_file, show_col_types = FALSE) %>%
            mutate(
                Participant_id = as.character(Participant_id),
                From_Time = normalise_timepoint_preserve_events(From_Time),
                To_Time = normalise_timepoint_preserve_events(To_Time),
                is_not_uti_to_uti = if ("is_not_uti_to_uti" %in% names(.)) as_pipeline_bool(is_not_uti_to_uti) else From_Status == "Not_UTI" & To_Status == "UTI"
            ) %>%
            filter(is_not_uti_to_uti %in% TRUE) %>%
            transmute(Participant_id, From_Time, To_Time, source = "vf_transition_case_index")
    }

    table10_file <- file.path(DIR_RESULTS, "summary", "table_10_not_uti_uti_transition_cases.csv")
    if (file.exists(table10_file)) {
        candidate_sources[["table_10"]] <- read_csv(table10_file, show_col_types = FALSE) %>%
            mutate(
                Participant_id = as.character(Participant_id),
                From_Time = normalise_timepoint_preserve_events(from_tp),
                To_Time = normalise_timepoint_preserve_events(to_tp)
            ) %>%
            filter(from_status == "Not_UTI", to_status == "UTI") %>%
            transmute(Participant_id, From_Time, To_Time, source = "table_10_not_uti_uti_transition_cases")
    }

    phen_file <- file.path(DIR_RESULTS, "longitudinal", "phenotype_switch_candidates.csv")
    if (file.exists(phen_file)) {
        candidate_sources[["phenotype_switch_candidates"]] <- read_csv(phen_file, show_col_types = FALSE) %>%
            mutate(
                Participant_id = as.character(Participant_id),
                From_Time = normalise_timepoint_preserve_events(From_Time),
                To_Time = normalise_timepoint_preserve_events(To_Time)
            ) %>%
            filter(From_Status == "Not_UTI", To_Status == "UTI") %>%
            transmute(Participant_id, From_Time, To_Time, source = "phenotype_switch_candidates")
    }

    if (length(candidate_sources) == 0) {
        return(tibble::tibble(
            Participant_id = character(),
            From_Time = character(),
            To_Time = character(),
            source = character()
        ))
    }

    bind_rows(candidate_sources) %>%
        distinct(Participant_id, From_Time, To_Time, .keep_all = TRUE)
}

candidate_pairs <- load_transition_candidates()
if (nrow(candidate_pairs) == 0) {
    warning("No Not_UTI -> UTI transition candidate table found; annotating all SNP pairs from annotated_snps.csv.")
    candidate_pairs <- snps %>%
        distinct(Participant_id, From_Time, To_Time) %>%
        mutate(source = "all_annotated_snp_pairs")
}

target_snps <- snps %>%
    inner_join(candidate_pairs %>% select(Participant_id, From_Time, To_Time),
               by = c("Participant_id", "From_Time", "To_Time")) %>%
    filter(Type == "SNP")

msg("Target SNPs: %d across %d Not_UTI -> UTI candidate pair(s)",
    nrow(target_snps), nrow(candidate_pairs))

empty_result <- target_snps %>%
    slice(0) %>%
    mutate(
        Region = character(),
        Gene = character(),
        Product = character(),
        Locus = character()
    )

out_file <- file.path(DIR_RESULTS, "longitudinal", "variant_annotation_detailed.csv")

if (nrow(target_snps) == 0) {
    msg("No target SNPs found for the configured candidate participants. Writing empty annotation table and exiting.")
    write_csv(empty_result, out_file)
    quit(save = "no", status = 0)
}

# 2. Helper: Parse GFF Manually
# ------------------------------------------------------------------------------
get_gff_data <- function(pid, tp) {
    tp_norm <- normalise_timepoint_preserve_events(tp)
    iso_row <- meta %>%
        filter(Participant_id == pid, tp_lab == tp_norm, !is.na(file_name)) %>%
        mutate(
            selected_rank = if_else(selected_canonical %in% TRUE, 0L, 1L),
            assembler_rank = case_when(
                tolower(assembler) == "longcycler" ~ 0L,
                tolower(assembler) == "flye" ~ 1L,
                TRUE ~ 2L
            ),
            size_rank = coalesce(total_bases, 0)
        ) %>%
        arrange(selected_rank, assembler_rank, desc(size_rank))
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

for (i in seq_len(nrow(pairs))) {
    p <- pairs[i, ]
    pid <- p$Participant_id
    t_ref <- p$From_Time

    msg("Processing Pair: %s (Ref: %s)", pid, t_ref)

    gff_df <- get_gff_data(pid, t_ref)
    if (is.null(gff_df)) next

    pair_snps <- target_snps %>%
        filter(Participant_id == pid, From_Time == t_ref)

    for (j in seq_len(nrow(pair_snps))) {
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
if (nrow(final_df) == 0) {
    final_df <- empty_result
}

write_csv(final_df, out_file)
msg("Saved annotated variants to %s", out_file)

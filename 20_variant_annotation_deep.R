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
#   Requires the canonical current Longcycler Not_UTI -> UTI case index and
#   fails closed unless its nine transition endpoints belong to the exact
#   Longcycler analysis cohort.
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

normalise_gff_sequence_id <- function(x) {
    x <- trimws(as.character(x))
    x <- sub("^>", "", x)
    sub("[[:space:]].*$", "", x)
}

extract_gff_attribute_strict <- function(attributes, key) {
    hit <- stringr::str_match(attributes, paste0("(?:^|;)", key, "=([^;]*)"))[, 2]
    ifelse(
        is.na(hit) | !nzchar(hit),
        NA_character_,
        vapply(hit, utils::URLdecode, character(1))
    )
}

read_gff_features_strict <- function(path) {
    lines <- readLines(path, warn = FALSE)
    fasta_marker <- which(lines == "##FASTA")
    if (length(fasta_marker) != 1L) {
        stop(basename(path), " must contain exactly one ##FASTA marker.")
    }
    marker <- fasta_marker[[1L]]
    feature_lines <- if (marker > 1L) lines[seq_len(marker - 1L)] else character()
    feature_lines <- feature_lines[nzchar(feature_lines) & !startsWith(feature_lines, "#")]
    fields <- strsplit(feature_lines, "\t", fixed = TRUE)
    if (any(lengths(fields) != 9L)) {
        stop(basename(path), " contains malformed non-comment GFF feature rows.")
    }
    feature_matrix <- if (length(fields)) {
        do.call(rbind, fields)
    } else {
        matrix(character(), nrow = 0L, ncol = 9L)
    }
    features <- tibble::as_tibble(feature_matrix, .name_repair = "minimal")
    names(features) <- c("seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes")
    features <- features %>%
        transmute(
            seqid = normalise_gff_sequence_id(.data$seqid),
            source = .data$source,
            type = .data$type,
            start = suppressWarnings(as.numeric(.data$start)),
            end = suppressWarnings(as.numeric(.data$end)),
            strand = .data$strand,
            attributes = .data$attributes,
            gene = extract_gff_attribute_strict(.data$attributes, "gene"),
            product = extract_gff_attribute_strict(.data$attributes, "product"),
            locus_tag = extract_gff_attribute_strict(.data$attributes, "locus_tag")
        )
    if (any(!is.finite(features$start) | !is.finite(features$end) |
            features$start < 1 | features$end < features$start)) {
        stop(basename(path), " contains invalid feature coordinates.")
    }

    fasta_lines <- if (marker < length(lines)) lines[(marker + 1L):length(lines)] else character()
    header_idx <- which(startsWith(fasta_lines, ">"))
    if (!length(header_idx)) stop(basename(path), " has no embedded FASTA records.")
    headers <- normalise_gff_sequence_id(fasta_lines[header_idx])
    if (any(!nzchar(headers)) || anyDuplicated(headers)) {
        stop(basename(path), " has missing or duplicated embedded FASTA sequence IDs.")
    }
    end_idx <- c(header_idx[-1L] - 1L, length(fasta_lines))
    contig_lengths <- vapply(seq_along(header_idx), function(i) {
        seq_lines <- if (end_idx[i] > header_idx[i]) fasta_lines[(header_idx[i] + 1L):end_idx[i]] else character()
        nchar(paste0(trimws(seq_lines), collapse = ""), type = "bytes")
    }, numeric(1))
    names(contig_lengths) <- headers
    feature_contig_lengths <- unname(contig_lengths[features$seqid])
    if (any(is.na(feature_contig_lengths)) || any(features$end > feature_contig_lengths)) {
        stop(basename(path), " contains features outside its embedded FASTA contigs.")
    }
    features
}

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
required_snp_provenance <- c(
    "Ref_Seqid", "SNP_Path", "SNP_SHA256",
    "Reference_FASTA_Path", "Reference_FASTA_SHA256"
)
missing_snp_provenance <- setdiff(required_snp_provenance, names(snps))
if (length(missing_snp_provenance)) {
    stop("annotated_snps.csv lacks exact reference provenance: ", paste(missing_snp_provenance, collapse = ", "))
}

gff_manifest_file <- file.path(DIR_WGS_PAN, "panaroo_input_manifest.csv")
if (!file.exists(gff_manifest_file)) stop("Missing strict Panaroo/GFF manifest: ", gff_manifest_file)
gff_manifest <- read_csv(gff_manifest_file, show_col_types = FALSE)
required_gff_manifest <- c("fasta_path", "fasta_sha256", "gff_path", "gff_sha256", "gff_available")
missing_gff_manifest <- setdiff(required_gff_manifest, names(gff_manifest))
if (length(missing_gff_manifest)) {
    stop("Panaroo input manifest lacks exact FASTA/GFF provenance: ", paste(missing_gff_manifest, collapse = ", "))
}
gff_manifest <- gff_manifest %>%
    mutate(
        fasta_path = normalizePath(fasta_path, winslash = "/", mustWork = FALSE),
        fasta_sha256 = as.character(fasta_sha256),
        gff_path = normalizePath(gff_path, winslash = "/", mustWork = FALSE),
        gff_sha256 = as.character(gff_sha256),
        gff_available = as_pipeline_bool(gff_available)
    )
if (any(!gff_manifest$gff_available | !file.exists(gff_manifest$gff_path))) {
    stop("Panaroo input manifest contains a missing required GFF. Rerun 12c_panaroo.R.")
}
observed_gff_sha256 <- vapply(gff_manifest$gff_path, digest::digest, character(1), algo = "sha256", file = TRUE)
if (any(is.na(gff_manifest$gff_sha256) | !nzchar(gff_manifest$gff_sha256) |
        observed_gff_sha256 != gff_manifest$gff_sha256)) {
    stop("One or more GFF files no longer match the SHA-256 recorded in the Panaroo input manifest.")
}
if (anyDuplicated(gff_manifest[c("fasta_path", "fasta_sha256")])) {
    stop("Panaroo input manifest contains duplicated selected FASTA path/SHA-256 pairs.")
}

case_index_file <- file.path(DIR_VF, "vf_transition_case_index.csv")
if (!file.exists(case_index_file)) stop("Missing canonical Longcycler transition case index: ", case_index_file)
case_index <- read_csv(case_index_file, show_col_types = FALSE)
required_case_cols <- c(
    "Participant_id", "From_Time", "To_Time", "From_Status", "To_Status",
    "is_not_uti_to_uti", "has_wgs_from", "has_wgs_to"
)
missing_case_cols <- setdiff(required_case_cols, names(case_index))
if (length(missing_case_cols)) stop("Canonical transition case index lacks: ", paste(missing_case_cols, collapse = ", "))
candidate_pairs <- case_index %>%
    mutate(
        Participant_id = as.character(Participant_id),
        From_Time = normalise_timepoint_preserve_events(From_Time),
        To_Time = normalise_timepoint_preserve_events(To_Time),
        is_not_uti_to_uti = as_pipeline_bool(is_not_uti_to_uti),
        has_wgs_from = as_pipeline_bool(has_wgs_from),
        has_wgs_to = as_pipeline_bool(has_wgs_to)
    ) %>%
    filter(
        is_not_uti_to_uti %in% TRUE,
        From_Status == "Not_UTI", To_Status == "UTI",
        has_wgs_from %in% TRUE, has_wgs_to %in% TRUE
    ) %>%
    transmute(Participant_id, From_Time, To_Time, source = "vf_transition_case_index") %>%
    distinct()
if (nrow(candidate_pairs) != 9L) {
    stop("Canonical current Longcycler Not_UTI -> UTI case index must contain exactly 9 WGS transitions; found ", nrow(candidate_pairs), ".")
}
if (anyDuplicated(candidate_pairs[c("Participant_id", "From_Time", "To_Time")])) {
    stop("Canonical Longcycler transition case index contains duplicated pair keys.")
}

if (!file.exists(FILE_ANALYSIS_CLINICAL_COHORT)) stop("Missing exact Longcycler analysis cohort: ", FILE_ANALYSIS_CLINICAL_COHORT)
analysis_cohort <- read_csv(FILE_ANALYSIS_CLINICAL_COHORT, show_col_types = FALSE) %>%
    mutate(
        Participant_id = as.character(Participant_id),
        tp_lab = normalise_timepoint_preserve_events(tp_lab)
    )
assert_analysis_assembly_manifest(
    analysis_cohort,
    context = "variant annotation Longcycler clinical cohort",
    require_selected = TRUE,
    require_qc = TRUE,
    require_files = TRUE,
    require_unique_episode = TRUE,
    require_nonempty = TRUE,
    require_assembly_id = TRUE
)
cohort_key <- paste(analysis_cohort$Participant_id, analysis_cohort$tp_lab, sep = "\n")
from_key <- paste(candidate_pairs$Participant_id, candidate_pairs$From_Time, sep = "\n")
to_key <- paste(candidate_pairs$Participant_id, candidate_pairs$To_Time, sep = "\n")
if (any(!from_key %in% cohort_key | !to_key %in% cohort_key)) {
    stop("Canonical transition case index contains endpoint(s) outside the exact Longcycler analysis cohort.")
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
        Locus = character(),
        GFF_Path = character(),
        GFF_SHA256 = character(),
        GFF_Assembly_ID = character()
    )

out_file <- file.path(DIR_RESULTS, "longitudinal", "variant_annotation_detailed.csv")

if (nrow(target_snps) == 0) {
    msg("No target SNPs found for the configured candidate participants. Writing empty annotation table and exiting.")
    write_csv(empty_result, out_file)
    quit(save = "no", status = 0)
}

# 2. Helper: Parse the exact GFF bound to a reference FASTA hash
# ------------------------------------------------------------------------------
get_gff_data <- function(reference_path, reference_sha256) {
    reference_path <- normalizePath(reference_path, winslash = "/", mustWork = TRUE)
    manifest_row <- gff_manifest %>%
        filter(fasta_sha256 == reference_sha256, fasta_path == reference_path)
    if (nrow(manifest_row) != 1L) {
        stop("Reference FASTA path/hash does not map to exactly one current Panaroo GFF: ", reference_path)
    }
    gff_file <- manifest_row$gff_path[1]
    msg("Loading GFF: %s", basename(gff_file))

    # Parse only the nine-column GFF feature section before ##FASTA.  The
    # shared strict parser rejects malformed feature rows and validates feature
    # coordinates against the embedded reference sequences; treating embedded
    # FASTA sequence lines as tabular rows previously emitted large readr
    # parsing-warning summaries.
    gff <- read_gff_features_strict(gff_file) %>%
        filter(type %in% c("CDS", "gene"))

    list(
        features = gff,
        gff_path = gff_file,
        gff_sha256 = manifest_row$gff_sha256[1],
        assembly_id = if ("Assembly_ID" %in% names(manifest_row)) as.character(manifest_row$Assembly_ID[1]) else NA_character_
    )
}

# 3. Annotate
# ------------------------------------------------------------------------------
results <- list()

pairs <- target_snps %>%
    select(Participant_id, From_Time, Reference_FASTA_Path, Reference_FASTA_SHA256) %>%
    distinct()

msg("Processing %d pairs", nrow(pairs))

for (i in seq_len(nrow(pairs))) {
    p <- pairs[i, ]
    pid <- p$Participant_id
    t_ref <- p$From_Time

    msg("Processing Pair: %s (Ref: %s)", pid, t_ref)

    gff_info <- get_gff_data(p$Reference_FASTA_Path, p$Reference_FASTA_SHA256)
    gff_df <- gff_info$features

    pair_snps <- target_snps %>%
        filter(
            Participant_id == pid, From_Time == t_ref,
            Reference_FASTA_Path == p$Reference_FASTA_Path,
            Reference_FASTA_SHA256 == p$Reference_FASTA_SHA256
        )

    for (j in seq_len(nrow(pair_snps))) {
        pos <- pair_snps$Pos_Ref[j]
        ref_seqid <- pair_snps$Ref_Seqid[j]

        hit <- gff_df %>%
            filter(seqid == ref_seqid, start <= pos, end >= pos) %>%
            mutate(.feature_priority = if_else(type == "CDS", 0L, 1L)) %>%
            arrange(.feature_priority, start, end) %>%
            slice_head(n = 1)

        if (nrow(hit) > 0) {
            res_row <- pair_snps[j, ] %>%
                mutate(
                    Region = hit$type,
                    Gene = hit$gene,
                    Product = hit$product,
                    Locus = hit$locus_tag,
                    GFF_Path = gff_info$gff_path,
                    GFF_SHA256 = gff_info$gff_sha256,
                    GFF_Assembly_ID = gff_info$assembly_id
                )
        } else {
            seqid_present <- ref_seqid %in% gff_df$seqid
            res_row <- pair_snps[j, ] %>%
                mutate(
                    Region = ifelse(seqid_present, "Intergenic", "Reference_contig_not_in_GFF"),
                    Gene = NA,
                    Product = NA,
                    Locus = NA,
                    GFF_Path = gff_info$gff_path,
                    GFF_SHA256 = gff_info$gff_sha256,
                    GFF_Assembly_ID = gff_info$assembly_id
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

#!/usr/bin/env Rscript
# ==============================================================================
# 12b_core_snp.R
# ==============================================================================
#
# GOAL:
#   Perform core genome SNP calling and alignment to quantify genetic distance
#   between isolates.  SNP distances are used by script 11 (strain comparison)
#   and script 15 (longitudinal patterns) to distinguish "same strain" from
#   "different strain" within the same participant.
#
# ------------------------------------------------------------------------------
# Role: [Phylogeny] - Perform core genome SNP calling and alignment.
#
# Inputs:
#   - results/qc/analysis_assembly_manifest.csv (validated Longcycler-only)
#   - exact FASTA paths declared by that manifest
#
# Outputs:
#   - results/wgs/core/core.aln.fasta
#   - results/wgs/core/snp_dists.tsv
#   - results/wgs/core/strain_pairs.csv
#   - results/wgs/core/core_genome.tree
#   - results/wgs/core/core_snp_sample_map.csv
#   - results/wgs/core/core_snp_command.log
#
# Usage:
#   Rscript 12b_core_snp.R
#
# Biological/Statistical purpose:
#   - Constructs a core genome alignment to infer phylogenetic relationships.
#   - Calculates pairwise SNP distances to identify "Same" vs "Different" strains.
# ==============================================================================

# 1. Load Configuration & Libraries
source("00_config.R")
source("R/wgs_helpers.R")

suppressPackageStartupMessages({
    library(dplyr)
    library(readr)
    library(fs)
    library(ape)
    library(tidyr)
})

# 2. Setup
# ------------------------------------------------------------------------------
DIR_CORE <- file.path(DIR_WGS, "core")
ensure_dir(DIR_CORE)

log_info("Starting 12b_core_snp.R")

tree_file <- file.path(DIR_CORE, "core_genome.tree")
dist_file <- file.path(DIR_CORE, "snp_dists.tsv")
pairs_file <- file.path(DIR_CORE, "strain_pairs.csv")
manifest_file <- file.path(DIR_CORE, "core_snp_input_manifest.csv")
hash_file <- file.path(DIR_CORE, "core_snp_input_manifest.hash")
stale_report <- file.path(DIR_CORE, "core_snp_staleness_report.txt")
command_log_file <- file.path(DIR_CORE, "core_snp_command.log")
sample_map_file <- file.path(DIR_CORE, "core_snp_sample_map.csv")
parsnp_out <- file.path(DIR_CORE, "parsnp_out")
aln_file <- file.path(parsnp_out, "parsnp.fasta")
xmfa_file <- file.path(parsnp_out, "parsnp.xmfa")
ggr_file <- file.path(parsnp_out, "parsnp.ggr")
force_core <- identical(Sys.getenv("FORCE_RERUN_CORE_SNP", "0"), "1")
resume_core <- identical(Sys.getenv("RESUME_CORE_SNP_ALIGNMENT", "0"), "1")

files_complete <- function(paths) {
    all(file.exists(paths) & !is.na(file.size(paths)) & file.size(paths) > 0)
}

sample_map_matches <- function(path, expected) {
    if (!files_complete(path)) return(FALSE)
    observed <- tryCatch(
        suppressMessages(read_csv(path, show_col_types = FALSE)),
        error = function(e) NULL
    )
    if (is.null(observed) || !identical(names(observed), names(expected)) ||
        nrow(observed) != nrow(expected)) {
        return(FALSE)
    }
    all(vapply(
        names(expected),
        function(column) identical(as.character(observed[[column]]), as.character(expected[[column]])),
        logical(1)
    ))
}

command_log_matches_manifest <- function(path, manifest_hash) {
    if (!files_complete(path)) return(FALSE)
    hash_anchor <- paste0("Input manifest SHA-256: ", manifest_hash)
    lines <- tryCatch(readLines(path, warn = FALSE), error = function(e) character())
    sum(lines == hash_anchor) == 1L
}

canonicalise_snp_distance_matrix <- function(dists_df, expected_labels) {
    expected_labels <- as.character(expected_labels)
    if (anyNA(expected_labels) || any(!nzchar(expected_labels)) || anyDuplicated(expected_labels)) {
        stop("Expected core-SNP labels are missing, empty, or duplicated.")
    }
    if (ncol(dists_df) < 3L || nrow(dists_df) != ncol(dists_df) - 1L) {
        stop("snp-dists output is not a square distance matrix.")
    }

    raw_labels <- as.character(dists_df[[1]])
    raw_matrix <- suppressWarnings(as.matrix(
        data.frame(lapply(dists_df[-1], as.numeric), check.names = FALSE)
    ))
    rownames(raw_matrix) <- raw_labels
    colnames(raw_matrix) <- names(dists_df)[-1]

    if (anyNA(raw_matrix) || anyDuplicated(raw_labels) ||
        !identical(rownames(raw_matrix), colnames(raw_matrix)) ||
        !isTRUE(all.equal(raw_matrix, t(raw_matrix), tolerance = 0)) ||
        any(raw_matrix < 0) || any(raw_matrix != floor(raw_matrix)) ||
        any(diag(raw_matrix) != 0)) {
        stop(
            "snp-dists output must have unique identical row/column labels and a ",
            "symmetric, non-negative integer matrix with a zero diagonal."
        )
    }

    # With `parsnp -r !`, Parsnp can retain the selected input genome and also
    # emit a synthetic copy of it suffixed `.ref`.  This is the only supported
    # label alias: it must map to one SHA-bound staged filename, and a retained
    # input/reference pair must have exactly identical distance vectors.
    reference_labels <- raw_labels[endsWith(raw_labels, ".ref")]
    if (length(reference_labels) > 1L) {
        stop("snp-dists output contains more than one Parsnp .ref label.")
    }
    canonical_raw_labels <- sub("[.]ref$", "", raw_labels)
    if (!setequal(unique(canonical_raw_labels), expected_labels)) {
        stop("Core-SNP distance labels do not map exactly to the SHA-bound staged input manifest.")
    }

    canonical_counts <- table(factor(canonical_raw_labels, levels = expected_labels))
    if (any(canonical_counts < 1L) || any(canonical_counts > 2L) ||
        sum(canonical_counts == 2L) > 1L) {
        stop("Core-SNP distance labels have unsupported multiplicity after .ref canonicalisation.")
    }

    reference_alias_collapsed <- FALSE
    if (length(reference_labels) == 1L) {
        reference_label <- reference_labels[[1]]
        staged_reference_label <- sub("[.]ref$", "", reference_label)
        if (!staged_reference_label %in% expected_labels) {
            stop("The Parsnp .ref label is not a SHA-bound staged input filename.")
        }

        if (staged_reference_label %in% raw_labels) {
            if (length(raw_labels) != length(expected_labels) + 1L ||
                canonical_counts[[staged_reference_label]] != 2L) {
                stop("Unexpected matrix dimensions for the retained Parsnp input/reference alias pair.")
            }
            ref_i <- match(reference_label, raw_labels)
            staged_i <- match(staged_reference_label, raw_labels)
            if (raw_matrix[ref_i, staged_i] != 0 ||
                !identical(unname(raw_matrix[ref_i, ]), unname(raw_matrix[staged_i, ])) ||
                !identical(unname(raw_matrix[, ref_i]), unname(raw_matrix[, staged_i]))) {
                stop(
                    "The Parsnp .ref entry and its SHA-bound staged input do not have ",
                    "identical zero-distance profiles."
                )
            }
            keep <- raw_labels != reference_label
            canonical_matrix <- raw_matrix[keep, keep, drop = FALSE]
            reference_alias_collapsed <- TRUE
        } else {
            if (length(raw_labels) != length(expected_labels)) {
                stop("Unexpected matrix dimensions when the Parsnp .ref label replaces its staged input label.")
            }
            canonical_matrix <- raw_matrix
            rownames(canonical_matrix)[rownames(canonical_matrix) == reference_label] <- staged_reference_label
            colnames(canonical_matrix)[colnames(canonical_matrix) == reference_label] <- staged_reference_label
        }
    } else {
        if (length(raw_labels) != length(expected_labels)) {
            stop("Core-SNP distance matrix dimensions do not match the SHA-bound staged input manifest.")
        }
        canonical_matrix <- raw_matrix
    }

    if (nrow(canonical_matrix) != length(expected_labels) ||
        anyDuplicated(rownames(canonical_matrix)) ||
        !setequal(rownames(canonical_matrix), expected_labels) ||
        !identical(rownames(canonical_matrix), colnames(canonical_matrix))) {
        stop("Canonical core-SNP distance matrix is not an exact one-to-one manifest mapping.")
    }

    # Publish a stable manifest-ordered matrix so every downstream row/column is
    # an exact staged filename rather than a Parsnp-specific reference alias.
    canonical_matrix <- canonical_matrix[expected_labels, expected_labels, drop = FALSE]
    canonical_df <- data.frame(
        sample = rownames(canonical_matrix),
        as.data.frame(canonical_matrix, check.names = FALSE),
        check.names = FALSE,
        row.names = NULL
    )
    list(
        data = canonical_df,
        matrix = canonical_matrix,
        reference_alias_collapsed = reference_alias_collapsed
    )
}

required_outputs <- c(aln_file, xmfa_file, dist_file, tree_file, pairs_file, command_log_file, sample_map_file)

# 3. Load the validated Longcycler-only inputs and fingerprint them before deciding
# whether an existing core-SNP directory is current.
selected_df <- load_analysis_assemblies(FILE_ANALYSIS_ASSEMBLY_MANIFEST, require_files = TRUE)
input_source <- FILE_ANALYSIS_ASSEMBLY_MANIFEST

manifest <- selected_df %>%
    transmute(
        Assembly_ID = as.character(Assembly_ID),
        Participant_id = as.character(Participant_id),
        tp_lab = normalise_timepoint_preserve_events(tp_lab),
        Assembler = normalise_assembler_column(selected_df),
        selection_policy_version = as.character(selection_policy_version),
        full_path = normalizePath(full_path, winslash = "/", mustWork = FALSE),
        file_size = file.size(full_path)
    ) %>%
    add_file_content_sha256("full_path", "fasta_sha256") %>%
    mutate(
        staged_extension = fasta_extension(full_path),
        staged_extension = ifelse(is.na(staged_extension), ".fasta", staged_extension),
        staged_fasta_name = paste0(
            sanitize_key_part(Assembly_ID), "__", substr(fasta_sha256, 1, 16),
            staged_extension
        )
    ) %>%
    select(-staged_extension) %>%
    arrange(Participant_id, tp_lab, Assembly_ID)

if (any(is.na(manifest$fasta_sha256) | nchar(manifest$fasta_sha256) != 64L)) {
    stop("Failed to create SHA-256 provenance for every core-SNP FASTA input.")
}
if (anyDuplicated(manifest$staged_fasta_name)) {
    stop("Core-SNP staged FASTA names are not unique after SHA-bound naming.")
}

temp_fasta_dir <- file.path(DIR_CORE, "temp_fastas")
sample_map <- manifest %>%
    transmute(
        parsnp_sample = stringr::str_remove(
            staged_fasta_name,
            stringr::regex("\\.(fasta|fna|fa)(\\.gz)?$", ignore_case = TRUE)
        ),
        parsnp_alignment_label = staged_fasta_name,
        parsnp_reference_label = paste0(staged_fasta_name, ".ref"),
        staged_fasta_name,
        staged_path = file.path(temp_fasta_dir, staged_fasta_name),
        Assembly_ID,
        Participant_id,
        tp_lab,
        Assembler,
        full_path,
        fasta_sha256
    )

write_csv(manifest, manifest_file)
current_hash <- hash_input_manifest(manifest)
previous_hash <- if (file.exists(hash_file)) readLines(hash_file, warn = FALSE)[1] else NA_character_
outputs_complete <- files_complete(required_outputs)

writeLines(
    c(
        "Core SNP input staleness report",
        sprintf("Generated: %s", format(Sys.time())),
        sprintf("Input source: %s", input_source),
        sprintf("Selected assemblies in current manifest: %d", nrow(manifest)),
        sprintf("Current manifest hash: %s", current_hash),
        sprintf("Previous manifest hash: %s", ifelse(is.na(previous_hash), "<none>", previous_hash)),
        sprintf("Required outputs complete: %s", outputs_complete),
        sprintf("Required outputs: %s", paste(required_outputs, collapse = "; ")),
        sprintf("FORCE_RERUN_CORE_SNP: %s", Sys.getenv("FORCE_RERUN_CORE_SNP", "0")),
        sprintf("RESUME_CORE_SNP_ALIGNMENT: %s", Sys.getenv("RESUME_CORE_SNP_ALIGNMENT", "0")),
        if (resume_core) "Status: RESUME REQUESTED - existing alignment provenance must pass before post-processing is recomputed."
        else if (outputs_complete && identical(current_hash, previous_hash) && !force_core) "Status: GREEN - complete outputs match current SHA-256 input manifest."
        else if (outputs_complete && identical(current_hash, previous_hash) && force_core) "Status: FORCED - matching complete outputs will be regenerated."
        else if (outputs_complete) "Status: RED - Core SNP outputs are stale relative to current assembly inputs."
        else if (any(file.exists(required_outputs))) "Status: INCOMPLETE - one or more required core SNP outputs are absent or empty."
        else "Status: MISSING - core SNP outputs are absent."
    ),
    stale_report
)

append_denominator_summary(
    manifest,
    "12b_core_snp.R",
    "core_snp_input_manifest",
    "participant_timepoint",
    manifest_file,
    "Fingerprint of canonical selected assemblies used to decide whether Parsnp outputs are current"
)

if (outputs_complete && identical(current_hash, previous_hash) && !force_core && !resume_core) {
    log_info("All required core-SNP outputs exist and the SHA-256 input-manifest hash matches.")
    log_info("Skipping Parsnp run.")
    quit(save = "no", status = 0)
}

if (outputs_complete && !identical(current_hash, previous_hash) && !force_core && !resume_core) {
    stop(
        "Core SNP outputs are stale relative to the Longcycler-only manifest. ",
        "Review ", stale_report, " and rerun with FORCE_RERUN_CORE_SNP=1; stale mixed-assembler outputs will not be accepted."
    )
}

# Reusing a multi-hour alignment is an explicit recovery mode, never an
# automatic fallback. It is allowed only when both surviving alignments are
# non-empty, the command log contains exactly the current manifest hash, and
# every prior sample-map field exactly matches the map derived above from that
# same SHA-bound manifest. Labels are checked again after snp-dists.
if (resume_core) {
    if (!files_complete(c(aln_file, xmfa_file))) {
        stop("RESUME_CORE_SNP_ALIGNMENT=1 requires non-empty parsnp.fasta and parsnp.xmfa files.")
    }
    if (!command_log_matches_manifest(command_log_file, current_hash)) {
        stop(
            "RESUME_CORE_SNP_ALIGNMENT=1 rejected: core_snp_command.log is not bound ",
            "exactly once to the current SHA-256 input-manifest hash."
        )
    }
    if (!sample_map_matches(sample_map_file, sample_map)) {
        stop(
            "RESUME_CORE_SNP_ALIGNMENT=1 rejected: the existing sample map does not ",
            "exactly match the current SHA/path/name-bound expected rows."
        )
    }
}
reuse_existing_alignment <- resume_core

# Check tools
has_snp_dists <- check_wgs_tool("snp-dists")
if (!has_snp_dists) stop("snp-dists is required for this module.")
snp_dists_bin <- unname(Sys.which("snp-dists"))
if (!reuse_existing_alignment) {
    has_parsnp <- check_wgs_tool("parsnp")
    if (!has_parsnp) stop("Parsnp is required for this module.")
    parsnp_bin <- unname(Sys.which("parsnp"))
}

valid_genomes <- manifest %>%
    filter(!is.na(full_path), file.exists(full_path)) %>%
    pull(full_path)

if (length(valid_genomes) < 2) {
    stop("Not enough valid genomes for core SNP analysis (Need >= 2, found ", length(valid_genomes), ")")
}

log_info("Processing ", length(valid_genomes), " selected QC-passing Longcycler genomes.")

# Once regeneration starts, the old success marker must not survive an
# interrupted or failed run. External outputs are also cleared so a failed
# command cannot be mistaken for a complete refreshed result.
unlink(c(hash_file, tree_file, dist_file, pairs_file), force = TRUE)
if (!reuse_existing_alignment) {
    unlink(c(command_log_file, sample_map_file), force = TRUE)
}

# 4. Prepare Input for Parsnp
# ------------------------------------------------------------------------------
# Parsnp takes a directory of FASTAs or a file list.
# We'll create a temporary directory with symlinks to ensure clean input.
if (reuse_existing_alignment) {
    log_info("Reusing Parsnp alignment bound to the exact current SHA-256 input manifest; rerunning post-processing only.")
} else {
    if (dir_exists(temp_fasta_dir)) dir_delete(temp_fasta_dir)
    dir_create(temp_fasta_dir)

    log_info("Staging FASTAs in ", temp_fasta_dir)
    write_csv(sample_map, sample_map_file)

    for (i in seq_len(nrow(manifest))) {
        link_name <- file.path(temp_fasta_dir, manifest$staged_fasta_name[[i]])
        fs::link_create(manifest$full_path[[i]], link_name)
    }

    # 5. Run Parsnp
    # --------------------------------------------------------------------------
    # -c: core genome alignment
    # -r !: random reference
    # -p: threads
    # -o: output dir
    # -x: Parsnp xtrafast mode in the installed Parsnp CLI
    if (dir_exists(parsnp_out)) dir_delete(parsnp_out) # Parsnp requires clean dir

    parsnp_args <- c(
        "-c",
        "-r", "!",
        "-d", temp_fasta_dir,
        "-o", parsnp_out,
        "-p", as.character(CORES_USE),
        "-x",
        "--verbose"
    )
    parsnp_capture <- tempfile("parsnp-command-")
    on.exit(unlink(parsnp_capture), add = TRUE)

    log_info("Running Parsnp...")
    res <- system2(parsnp_bin, args = parsnp_args, stdout = parsnp_capture, stderr = parsnp_capture)
    writeLines(
        c(
            sprintf("Generated: %s", format(Sys.time())),
            sprintf("Input manifest SHA-256: %s", current_hash),
            paste(c(shQuote(parsnp_bin), vapply(parsnp_args, shQuote, character(1))), collapse = " "),
            "",
            readLines(parsnp_capture, warn = FALSE)
        ),
        command_log_file
    )

    # Check if alignment exists even if Parsnp failed (e.g. at tree step)
    if (res != 0) {
        if (file.exists(xmfa_file)) {
            log_warn("Parsnp exited with error (likely RAxML failure), but alignment exists. Proceeding...")
        } else {
            log_error("Parsnp failed with exit code ", res)
            stop("Parsnp execution failed and no alignment found.")
        }
    } else {
        log_info("Parsnp completed successfully.")
    }
}

# 5b. Convert GGR to FASTA (if needed)
# ------------------------------------------------------------------------------
if (!file.exists(aln_file) && file.exists(ggr_file)) {
    log_info("Converting GGR to FASTA...")
    harvesttools_bin <- unname(Sys.which("harvesttools"))
    if (!nzchar(harvesttools_bin)) stop("harvesttools is required to convert Parsnp GGR output to FASTA.")
    conv_res <- system2(harvesttools_bin, args = c("-i", ggr_file, "-M", aln_file))
    if (conv_res != 0L) stop("harvesttools failed to convert Parsnp output (exit ", conv_res, ").")
}
if (!files_complete(c(xmfa_file, aln_file))) stop("Parsnp completed without complete XMFA and FASTA alignments.")

# 6. Run snp-dists & Build Tree
# ------------------------------------------------------------------------------
log_info("Running snp-dists...")
snp_dists_stderr <- tempfile("snp-dists-command-")
raw_dist_file <- tempfile("snp-dists-raw-", tmpdir = DIR_CORE, fileext = ".tsv")
on.exit(unlink(snp_dists_stderr), add = TRUE)
on.exit(unlink(raw_dist_file), add = TRUE)
dists_res <- system2(snp_dists_bin, args = aln_file, stdout = raw_dist_file, stderr = snp_dists_stderr)
cat(
    "\n\n", paste(c(shQuote(snp_dists_bin), shQuote(aln_file)), collapse = " "), "\n",
    paste(readLines(snp_dists_stderr, warn = FALSE), collapse = "\n"), "\n",
    file = command_log_file,
    append = TRUE,
    sep = ""
)
if (dists_res != 0L || !files_complete(raw_dist_file)) {
    stop("snp-dists failed or produced an empty distance table (exit ", dists_res, ").")
}

# 6b. Build Tree (Robust Fix for "Too Few Species")
# ------------------------------------------------------------------------------
# A deterministic NJ tree is built from the required SNP-distance matrix even
# when Parsnp's optional RAxML tree step fails on identical sequences.
log_info("Building Neighbor-Joining tree from SNP distances...")

dists_raw <- read_tsv(raw_dist_file, show_col_types = FALSE, name_repair = "minimal")
canonical_dists <- canonicalise_snp_distance_matrix(
    dists_raw,
    sample_map$parsnp_alignment_label
)
dists_df <- canonical_dists$data
dist_mat <- canonical_dists$matrix
if (canonical_dists$reference_alias_collapsed) {
    log_info("Collapsed Parsnp's validated zero-distance .ref alias to its SHA-bound staged input label.")
}
write_tsv(dists_df, dist_file)
log_info("Canonical SNP distances written to ", dist_file)

tree <- ape::nj(as.dist(dist_mat))

ape::write.tree(tree, file = tree_file)
log_info("Tree saved to ", tree_file)

# 7. Pairwise Classification
# ------------------------------------------------------------------------------
log_info("Classifying strain pairs...")

dists <- dists_df

colnames(dists)[1] <- "A"

pairs_long <- dists %>%
    pivot_longer(-A, names_to = "B", values_to = "snps") %>%
    mutate(snps = as.numeric(snps)) %>%
    filter(A < B)

same_strain_snp_threshold <- strain_snp_threshold()

pairs_classified <- pairs_long %>%
    mutate(
        call = case_when(
            snps <= same_strain_snp_threshold ~ "Same",
            snps <= 1000 ~ "Related",
            TRUE ~ "Different"
        )
    )

write_csv(pairs_classified, pairs_file)
log_info("Classified ", nrow(pairs_classified), " pairs. Saved to ", pairs_file)

log_info("Pair Summary:")
print(table(pairs_classified$call))

if (files_complete(required_outputs)) {
    writeLines(current_hash, hash_file)
    writeLines(
        c(
            "Core SNP input staleness report",
            sprintf("Generated: %s", format(Sys.time())),
            sprintf("Input source: %s", input_source),
            sprintf("Selected assemblies in current manifest: %d", nrow(manifest)),
            sprintf("Current manifest hash: %s", current_hash),
            sprintf("Previous manifest hash: %s", ifelse(is.na(previous_hash), "<none>", previous_hash)),
            sprintf("Required outputs: %s", paste(required_outputs, collapse = "; ")),
            "Status: GREEN - all required outputs were generated/refreshed for the current SHA-256 input manifest."
        ),
        stale_report
    )
} else {
    missing_outputs <- required_outputs[!file.exists(required_outputs) | is.na(file.size(required_outputs)) | file.size(required_outputs) <= 0]
    stop("Core-SNP run ended without all required outputs: ", paste(missing_outputs, collapse = ", "))
}

# Cleanup
# dir_delete(temp_fasta_dir) # Optional: keep for debugging

log_info("12b_core_snp.R complete.")

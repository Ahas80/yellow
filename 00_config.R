# ==============================================================================
# 00_config.R
# ==============================================================================
#
# GOAL:
#   Central configuration for the entire rUTIs / YELLOW RoUTIne pipeline.
#   Every numbered R script sources this file first to get consistent paths,
#   helpers, and directory setup.
#
# WHY THIS FILE EXISTS:
#   Without a single config, each script would hardcode its own paths and
#   directory-creation logic.  This file guarantees that:
#     - All scripts agree on where inputs and outputs live
#     - Output directories are created automatically
#     - Common helpers (logging, directory creation, metadata loading) are
#       available everywhere without duplication
#
# WHAT IT PROVIDES:
#   Section 1 — Global paths (DIR_ROOT, DIR_VF, DIR_MLST, etc.)
#   Section 2 — Key file paths (FILE_VF_PA, FILE_MLST_ALL, etc.)
#   Section 3 — Batch configuration (BATCH_DEFINITIONS, clinical/assembly paths)
#   Section 4 — Settings (parallel cores)
#   Section 5 — Helper functions (msg, ensure_dir, check_tool, load_metadata)
#   Section 6 — Automatic directory initialisation
#
# USAGE:
#   source("00_config.R")   # at the top of every pipeline script
# ==============================================================================

# ------------------------------------------------------------------------------
# 1. Global Paths
# ------------------------------------------------------------------------------
# Use here::here() for robust project root detection
if (!requireNamespace("here", quietly = TRUE)) {
    stop("Package 'here' is required. Please install it with install.packages('here')")
}
DIR_ROOT <- here::here()

# Input Directories
DIR_FASTAS <- file.path(DIR_ROOT, "ont-yellow-routine-fastas")

# Output Directories
DIR_RESULTS <- file.path(DIR_ROOT, "results")

# Level 1 Results
DIR_CLINICAL <- file.path(DIR_RESULTS, "clinical")
DIR_VF <- file.path(DIR_RESULTS, "vf")
DIR_MLST <- file.path(DIR_RESULTS, "mlst")
DIR_PLASMIDS <- file.path(DIR_RESULTS, "plasmids")
DIR_STRAIN <- file.path(DIR_RESULTS, "strain_compare")
DIR_WGS <- file.path(DIR_RESULTS, "wgs")
DIR_MODELS <- file.path(DIR_RESULTS, "models")
DIR_QC <- file.path(DIR_RESULTS, "qc") # General QC
DIR_REPORTS <- file.path(DIR_RESULTS, "reports")
DIR_METADATA <- file.path(DIR_RESULTS, "metadata")

# WGS Subdirectories
DIR_WGS_QC <- file.path(DIR_WGS, "qc")
DIR_WGS_CORE <- file.path(DIR_WGS, "core")
DIR_WGS_PAN <- file.path(DIR_WGS, "pan")
DIR_WGS_KMER <- file.path(DIR_WGS, "kmer")
DIR_WGS_PLASMIDS <- file.path(DIR_WGS, "plasmids")
DIR_WGS_SV <- file.path(DIR_WGS, "sv")
DIR_PROKKA <- file.path(DIR_RESULTS, "prokka")
DIR_PROKKA_SLIM <- file.path(DIR_RESULTS, "prokka_prefixed_slim")

# Plots Directories
DIR_PLOTS <- file.path(DIR_ROOT, "plots")
DIR_PLOTS_CLINICAL <- file.path(DIR_PLOTS, "clinical")
DIR_PLOTS_VF <- file.path(DIR_PLOTS, "vf")
DIR_PLOTS_MLST <- file.path(DIR_PLOTS, "mlst")
DIR_PLOTS_PLASMIDS <- file.path(DIR_PLOTS, "plasmids")
DIR_PLOTS_STRAIN <- file.path(DIR_PLOTS, "strain_compare")
DIR_PLOTS_WGS <- file.path(DIR_PLOTS, "wgs")
DIR_PLOTS_MODELS <- file.path(DIR_PLOTS, "models")

# Logs Directory
DIR_LOGS <- file.path(DIR_ROOT, "logs")


# ------------------------------------------------------------------------------
# 2. Files
# ------------------------------------------------------------------------------
FILE_METADATA <- file.path(DIR_ROOT, "assembly_metadata.csv")
FILE_ASSEMBLIES <- file.path(DIR_ROOT, "assemblies.list")
FILE_STATUS_MAP <- file.path(DIR_CLINICAL, "status_map.csv")
FILE_STATUS_MAP_POSTER <- file.path(DIR_CLINICAL, "status_map_with_poster_tp.csv")
FILE_STATUS_MAP_LEGACY_COMPARISON <- file.path(DIR_CLINICAL, "status_map_legacy_comparison.csv")
FILE_UTI_CLASSIFICATION_AUDIT <- file.path(DIR_CLINICAL, "uti_binary_classification_audit.csv")
FILE_UTI_RECLASSIFICATION_MOVEMENT <- file.path(DIR_CLINICAL, "uti_reclassification_movement_table.csv")
FILE_UTI_SYMPTOM_RULE_AUDIT <- file.path(DIR_CLINICAL, "uti_symptom_rule_audit.csv")
FILE_UTI_CFU_THRESHOLD_AUDIT <- file.path(DIR_CLINICAL, "uti_cfu_threshold_audit.csv")
FILE_UTI_DETECTED_COLUMN_MAP <- file.path(DIR_CLINICAL, "uti_detected_column_map.csv")
FILE_MLST_LOCAL_ALL <- file.path(DIR_MLST, "mlst_all.tsv")
FILE_MLST_LOCAL_CANONICAL <- file.path(DIR_MLST, "mlst_with_meta.csv")
FILE_MLST_SOURCE_COMPARISON_DIR <- file.path(DIR_RESULTS, "mlst_source_comparison")
FILE_MLST_PROVIDER_NORMALIZED <- file.path(FILE_MLST_SOURCE_COMPARISON_DIR, "provider_mlst_normalized.csv")
FILE_MLST_PROVIDER_PREFERRED <- file.path(DIR_MLST, "mlst_provider_preferred.csv")
FILE_MLST_PROVIDER_PREFERRED_ALL <- file.path(DIR_MLST, "mlst_provider_preferred_all.csv")
FILE_MLST_PROVIDER_SOURCE_AUDIT <- file.path(DIR_MLST, "mlst_provider_source_audit.csv")
FILE_MLST_ALL <- FILE_MLST_LOCAL_ALL
FILE_MLST_CANONICAL <- FILE_MLST_PROVIDER_PREFERRED
FILE_MLST_EPISODE <- FILE_MLST_CANONICAL
FILE_MLST_MATRIX <- file.path(DIR_MLST, "mlst_matrix.csv")
FILE_MLST_ISOLATE_EXPLORATORY <- file.path(DIR_MLST, "mlst_isolate_with_metadata_exploratory.csv")
FILE_VF_HITS <- file.path(DIR_VF, "vf_hits_all.rds")
FILE_VF_PA <- file.path(DIR_VF, "vf_pa_all.csv")
FILE_VF_READY <- file.path(DIR_VF, "vf_analysis_ready.csv")
FILE_VF_BINARY_UTI_READY <- file.path(DIR_VF, "vf_binary_uti_ready.csv")
FILE_UTI_BINARY_MODEL_WARNINGS <- file.path(DIR_MODELS, "uti_binary_model_warnings.txt")
FILE_MANUAL_SAMPLE_CURATION <- file.path(DIR_ROOT, "data", "manual_sample_curation.csv")
FILE_SAMPLE_CURATION_AUDIT <- file.path(DIR_QC, "manual_sample_curation_applied.csv")
FILE_QUARANTINED_FASTA_EXPECTATIONS <- file.path(DIR_QC, "quarantined_failed_or_not_expected_fastas.csv")
FILE_DUPLICATE_CULTURE_QC <- file.path(DIR_QC, "duplicate_culture_qc.csv")

# New: canonical all-batches assembly metadata (generated by 00_make_assembly_metadata.r)
FILE_ASSEMBLY_META_ALL <- file.path(DIR_METADATA, "assembly_metadata_all_batches.csv")

# Overview spreadsheet: master reference for expected isolates across all batches
# This file was provided by the lab and contains the authoritative list of all
# E. coli isolates sent for sequencing in batches 1–6 plus extras.
FILE_OVERVIEW_XLSX <- file.path(DIR_ROOT, "OVERVIEW E.coli batch 1-6 - CLEAN_RC.xlsx")

# ------------------------------------------------------------------------------
# 3. Batch Configuration
# ------------------------------------------------------------------------------
# BATCH_DEFINITIONS: Central definition of all known batches.
# To add a future batch, add a single entry here. All downstream scripts
# will automatically pick it up because they read this list.
#
# Each entry specifies:
#   - batch_id:       Integer batch number
#   - label:          Human-readable label
#   - clinical_csv:   Filename of the clinical CSV in data/inputs/ (NA if none)
#   - expected_count: Expected number of E. coli isolates (from study protocol)
#   - notes:          Free-text notes
BATCH_DEFINITIONS <- list(
    list(batch_id = 1L, label = "Batch 1", clinical_csv = "batch1.csv",
         expected_count = 96L, notes = "First sequencing batch"),
    list(batch_id = 2L, label = "Batch 2", clinical_csv = "batch2.csv",
         expected_count = 96L, notes = "Second sequencing batch"),
    list(batch_id = 3L, label = "Batch 3", clinical_csv = "batch3.csv",
         expected_count = 96L, notes = "Third sequencing batch"),
    list(batch_id = 4L, label = "Batch 4", clinical_csv = "batch4.csv",
         expected_count = 96L, notes = "Fourth sequencing batch"),
    list(batch_id = 5L, label = "Batch 5", clinical_csv = "batch5.csv",
         expected_count = 96L, notes = "Fifth sequencing batch"),
    list(batch_id = 6L, label = "Batch 6", clinical_csv = "batch6.csv",
         expected_count = 119L, notes = "Sixth batch (96 + 23 extra to reach ~600 total)")
)

# Convenience: extract the batch IDs and clinical CSV file list
BATCH_IDS <- vapply(BATCH_DEFINITIONS, `[[`, integer(1), "batch_id")
BATCH_CLINICAL_CSVS <- vapply(BATCH_DEFINITIONS, `[[`, character(1), "clinical_csv")
BATCH_EXPECTED_COUNTS <- vapply(BATCH_DEFINITIONS, `[[`, integer(1), "expected_count")

# Assembly file extensions we accept
ASSEMBLY_EXTENSIONS <- c("fasta", "fa", "fna")
ASSEMBLY_EXTENSIONS_GZ <- paste0(ASSEMBLY_EXTENSIONS, ".gz")
ASSEMBLY_PATTERN <- paste0("\\.(", paste(c(ASSEMBLY_EXTENSIONS, ASSEMBLY_EXTENSIONS_GZ),
                                          collapse = "|"), ")$")

# Primary-analysis assembly policy. Generated metadata and QC products contain
# only Longcycler assemblies; no alternative assembler is an analysis input.
KNOWN_ASSEMBLERS <- c("longcycler")
ANALYSIS_ASSEMBLER <- "longcycler"
ALLOW_ASSEMBLER_FALLBACK <- FALSE
ASSEMBLY_SELECTION_POLICY_VERSION <- "longcycler_only_qcpass_v1"
FILE_CANONICAL_ASSEMBLY_SELECTION <- file.path(DIR_QC, "canonical_assembly_selection.csv")
FILE_ANALYSIS_ASSEMBLY_MANIFEST <- file.path(DIR_QC, "analysis_assembly_manifest.csv")
FILE_ANALYSIS_CLINICAL_COHORT <- file.path(DIR_CLINICAL, "analysis_cohort_longcycler.csv")

# Overview spreadsheet sheet name
OVERVIEW_SHEET <- "Batches overzicht"

# Overview spreadsheet header row (row 4 in the file, data starts at row 5)
OVERVIEW_HEADER_ROW <- 4L
OVERVIEW_DATA_START_ROW <- 5L

# Overview spreadsheet key column names (after reading with proper header)
OVERVIEW_COL_BATCH <- "Batch"
OVERVIEW_COL_PID <- "Participant_id"
OVERVIEW_COL_ISOLATE_ID <- "Isolaat id"
OVERVIEW_COL_ORGANISM <- "Organisme"
OVERVIEW_COL_TIMEPOINT <- "Meetmoment"
OVERVIEW_COL_BEOORD <- "Beoord"
OVERVIEW_COL_KIEMGETAL <- "Kiemgetal"

# Total expected E. coli samples across all batches (protocol target)
TOTAL_EXPECTED_ECOLI <- 600L

# ------------------------------------------------------------------------------
# 4. Settings & Constants
# ------------------------------------------------------------------------------
# Parallel processing
CORES_MAX <- parallel::detectCores(logical = FALSE)
CORES_USE <- max(1, CORES_MAX - 1)

# Clinical UTI definition for the primary analysis.
UTI_DEFINITION_VERSION <- "catheter_adjusted_sns_cfu1e3_v1"
UTI_CFU_THRESHOLD_PRIMARY <- 1000
UTI_CFU_THRESHOLD_LEGACY <- 100000

# Same-strain threshold used by longitudinal WGS/VF interpretation.
# Prior YELLOW study work used <=25 SNPs as the same-strain boundary.
SAME_STRAIN_SNP_THRESHOLD <- 25L

# ------------------------------------------------------------------------------
# 5. Helper Functions
# ------------------------------------------------------------------------------

#' Log a message with timestamp
#' @param ... Message components passed to sprintf
msg <- function(...) {
    message(format(Sys.time(), "[%H:%M:%S] "), sprintf(...))
}

#' Ensure a directory exists, creating it if necessary
#' @param path Path to the directory
ensure_dir <- function(path) {
    if (!dir.exists(path)) {
        dir.create(path, recursive = TRUE, showWarnings = FALSE)
        message("Created directory: ", path)
    }
}

#' Check if a system tool is available
#' @param tool_name Name of the command line tool
#' @return TRUE if found, stops execution if not found
check_tool <- function(tool_name) {
    path <- Sys.which(tool_name)
    if (path == "") {
        stop(sprintf("Critical tool '%s' not found in PATH. Please install it.", tool_name))
    }
    return(TRUE)
}

#' Log session information and tool versions
#' @param tools Vector of tool names to check versions for (optional)
log_session_info <- function(tools = NULL) {
    msg("R Session Info:")
    print(sessionInfo())

    if (!is.null(tools)) {
        msg("External Tool Versions:")
        for (tool in tools) {
            tryCatch(
                {
                    ver <- system(paste(tool, "--version"), intern = TRUE)
                    msg("%s: %s", tool, paste(ver, collapse = " "))
                },
                error = function(e) {
                    msg("%s: Version check failed", tool)
                }
            )
        }
    }
}

#' Load Assembly Metadata
#' @return Tibble containing assembly metadata
load_metadata <- function() {
    if (!file.exists(FILE_METADATA)) {
        stop("Metadata file not found: ", FILE_METADATA, "\nRun 00_make_assembly_metadata.r first.")
    }
    readr::read_csv(FILE_METADATA, show_col_types = FALSE)
}

#' Load the overview spreadsheet (OVERVIEW E.coli batch 1-6 - CLEAN_RC.xlsx)
#' Returns a tibble with one row per expected E. coli isolate across all batches.
#' @return Tibble with columns: Batch, Participant_id, Isolate_ID, Organism, Timepoint, etc.
load_overview_spreadsheet <- function() {
    if (!file.exists(FILE_OVERVIEW_XLSX)) {
        warning("Overview spreadsheet not found: ", FILE_OVERVIEW_XLSX)
        return(NULL)
    }
    if (!requireNamespace("readxl", quietly = TRUE)) {
        warning("Package 'readxl' is required to read the overview spreadsheet. Install with install.packages('readxl')")
        return(NULL)
    }
    # Read the spreadsheet — header is in row 4, data starts row 5
    raw <- readxl::read_excel(
        FILE_OVERVIEW_XLSX,
        sheet = OVERVIEW_SHEET,
        skip = OVERVIEW_HEADER_ROW - 1,  # skip rows before the header
        col_names = TRUE
    )
    # Standardize column names for downstream use
    if (OVERVIEW_COL_ISOLATE_ID %in% names(raw)) {
        raw <- dplyr::rename(raw, Isolate_ID = !!OVERVIEW_COL_ISOLATE_ID)
    }
    if (OVERVIEW_COL_ORGANISM %in% names(raw)) {
        raw <- dplyr::rename(raw, Organism = !!OVERVIEW_COL_ORGANISM)
    }
    if (OVERVIEW_COL_TIMEPOINT %in% names(raw)) {
        raw <- dplyr::rename(raw, Timepoint = !!OVERVIEW_COL_TIMEPOINT)
    }
    if (OVERVIEW_COL_BEOORD %in% names(raw)) {
        raw <- dplyr::rename(raw, Beoordeling = !!OVERVIEW_COL_BEOORD)
    }
    if (OVERVIEW_COL_KIEMGETAL %in% names(raw)) {
        raw <- dplyr::rename(raw, Kiemgetal = !!OVERVIEW_COL_KIEMGETAL)
    }
    if (OVERVIEW_COL_PID %in% names(raw)) {
        raw$Participant_id <- as.character(raw[[OVERVIEW_COL_PID]])
    }
    if (OVERVIEW_COL_BATCH %in% names(raw)) {
        raw$Batch <- as.integer(raw[[OVERVIEW_COL_BATCH]])
    }
    # Filter to only E. coli rows (all rows in the clean file should be E. coli,
    # but guard against any non-E. coli entries)
    if ("Organism" %in% names(raw)) {
        raw <- dplyr::filter(
            raw,
            grepl("Escherichia", Organism, ignore.case = TRUE) | is.na(Organism)
        )
    }
    raw
}

#' Get list of available clinical batch CSV filenames
#' Only returns filenames for CSVs that actually exist on disk.
#' @return Character vector of available batch CSV paths
get_available_clinical_csvs <- function() {
    csvs <- BATCH_CLINICAL_CSVS[!is.na(BATCH_CLINICAL_CSVS)]
    paths <- file.path("data", "inputs", csvs)
    existing <- file.exists(paths)
    if (sum(!existing) > 0) {
        msg("⚠ Clinical CSV(s) not found: %s",
            paste(csvs[!existing], collapse = ", "))
    }
    paths[existing]
}

# Shared QA helpers for canonical FASTA discovery, episode IDs, denominator
# logging, and WGS input manifests. This is optional for legacy compatibility
# but required by the repaired source-level QA workflow.
PIPELINE_QC_HELPERS <- file.path(DIR_ROOT, "R", "pipeline_qc_helpers.R")
if (file.exists(PIPELINE_QC_HELPERS)) {
    source(PIPELINE_QC_HELPERS)
}

STRAIN_CONTEXT_HELPERS <- file.path(DIR_ROOT, "R", "strain_context_helpers.R")
if (file.exists(STRAIN_CONTEXT_HELPERS)) {
    source(STRAIN_CONTEXT_HELPERS)
}

# ------------------------------------------------------------------------------
# 6. Initialization
# ------------------------------------------------------------------------------
# Create essential output directories immediately
# Create essential output directories immediately
ensure_dir(DIR_RESULTS)
ensure_dir(DIR_PLOTS)

# Ensure Level 1 Subdirs
for (d in c(DIR_CLINICAL, DIR_VF, DIR_MLST, DIR_PLASMIDS, DIR_STRAIN, DIR_WGS,
            DIR_MODELS, DIR_QC, DIR_METADATA)) {
    ensure_dir(d)
}

# Ensure Plot Subdirs
for (d in c(DIR_PLOTS_CLINICAL, DIR_PLOTS_VF, DIR_PLOTS_MLST, DIR_PLOTS_PLASMIDS,
            DIR_PLOTS_STRAIN, DIR_PLOTS_WGS, DIR_PLOTS_MODELS)) {
    ensure_dir(d)
}

# Ensure Logs Directory
ensure_dir(DIR_LOGS)

message("Loaded configuration from 00_config.R")

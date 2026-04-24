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
#   Section 3 — Settings (parallel cores)
#   Section 4 — Helper functions (msg, ensure_dir, check_tool, load_metadata)
#   Section 5 — Automatic directory initialisation
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
FILE_MLST_ALL <- file.path(DIR_MLST, "mlst_all.tsv")
FILE_VF_HITS <- file.path(DIR_VF, "vf_hits_all.rds")
FILE_VF_PA <- file.path(DIR_VF, "vf_pa_all.csv")

# ------------------------------------------------------------------------------
# 3. Settings & Constants
# ------------------------------------------------------------------------------
# Parallel processing
CORES_MAX <- parallel::detectCores(logical = FALSE)
CORES_USE <- max(1, CORES_MAX - 1)

# ------------------------------------------------------------------------------
# 4. Helper Functions
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

# ------------------------------------------------------------------------------
# 5. Initialization
# ------------------------------------------------------------------------------
# Create essential output directories immediately
# Create essential output directories immediately
ensure_dir(DIR_RESULTS)
ensure_dir(DIR_PLOTS)

# Ensure Level 1 Subdirs
for (d in c(DIR_CLINICAL, DIR_VF, DIR_MLST, DIR_PLASMIDS, DIR_STRAIN, DIR_WGS, DIR_MODELS, DIR_QC)) {
    ensure_dir(d)
}

# Ensure Plot Subdirs
for (d in c(DIR_PLOTS_CLINICAL, DIR_PLOTS_VF, DIR_PLOTS_MLST, DIR_PLOTS_PLASMIDS, DIR_PLOTS_STRAIN, DIR_PLOTS_WGS, DIR_PLOTS_MODELS)) {
    ensure_dir(d)
}

# Ensure Logs Directory
ensure_dir(DIR_LOGS)

message("Loaded configuration from 00_config.R")

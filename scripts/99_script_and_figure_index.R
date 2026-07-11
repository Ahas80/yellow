#!/usr/bin/env Rscript
# ==============================================================================
# 99_script_and_figure_index.R
# ------------------------------------------------------------------------------
# Role: [Meta] - Generate an index of all scripts and their functions.
#
# Inputs:
#   - Root numbered R scripts and non-legacy helper scripts under scripts/
#
# Outputs:
#   - results/meta/script_and_figure_index.csv
#
# Purpose:
#   - Programmatically parses the headers of all pipeline scripts to create
#     a summary table of roles, inputs, and outputs.
# ==============================================================================

# Adjust paths for running from scripts/ subdirectory
root_dir <- ".."
config_path <- file.path(root_dir, "00_config.R")

if (!file.exists(config_path)) {
    # Fallback if running from root
    if (file.exists("00_config.R")) {
        root_dir <- "."
        config_path <- "00_config.R"
    } else {
        stop("Cannot find 00_config.R")
    }
}

source(config_path)
suppressPackageStartupMessages({
    library(dplyr)
    library(readr)
    library(stringr)
    library(fs)
})

msg("Starting 99_script_and_figure_index.R")

# 1. Identify Scripts
# ------------------------------------------------------------------------------
# Look for active root pipeline scripts and non-legacy scripts/ utilities.
# Use list.files here because fs::dir_ls regex matching differs for "." and ".."
# paths and previously produced an empty/stale index when run from the repo root.
root_scripts <- list.files(root_dir, pattern = "^[0-9]{2}.*\\.[Rr]$", full.names = TRUE)
scripts_dir <- file.path(root_dir, "scripts")
utility_scripts <- if (dir.exists(scripts_dir)) {
    list.files(
        scripts_dir,
        pattern = "\\.(R|r|sh|awk)$",
        full.names = TRUE,
        recursive = FALSE,
        ignore.case = TRUE
    )
} else {
    character()
}

scripts <- sort(unique(c(root_scripts, utility_scripts)))
msg("Found %d scripts to index.", length(scripts))

# 2. Parse Headers
# ------------------------------------------------------------------------------
parse_script <- function(f) {
    lines <- readLines(f, n = 100) # Read first 100 lines

    # Helper to extract content after a tag
    extract_tag <- function(tag) {
        # Find line with tag
        idx <- grep(paste0("^#\\s*", tag), lines)
        if (length(idx) == 0) {
            return(NA_character_)
        }

        # Extract content from that line
        content <- sub(paste0("^#\\s*", tag, "\\s*"), "", lines[idx[1]])

        # Check for continuation lines (indented comments following the tag)
        # This is a simple heuristic: look at next lines starting with "#   -" or "#     "
        j <- idx[1] + 1
        while (j <= length(lines)) {
            if (grepl("^#\\s{2,}-", lines[j]) || grepl("^#\\s{4,}", lines[j])) {
                extra <- sub("^#\\s*", "", lines[j])
                content <- paste(content, extra, sep = "; ")
                j <- j + 1
            } else {
                break
            }
        }
        return(content)
    }

    tibble(
        Script = fs::path_file(f),
        Role = extract_tag("Role:"),
        Inputs = extract_tag("Inputs:"),
        Outputs = extract_tag("Outputs:"),
        Purpose = extract_tag("Biological/Statistical purpose:")
    )
}

# 3. Build Index
# ------------------------------------------------------------------------------
index <- if (length(scripts) > 0) {
    purrr::map_dfr(scripts, parse_script)
} else {
    tibble(Script = character(), Role = character(), Inputs = character(),
           Outputs = character(), Purpose = character())
}

# Clean up fields
index <- index %>%
    mutate(
        Role = str_remove_all(Role, "\\[.*?\\]\\s*-\\s*"), # Remove [Tag] - prefix if desired, or keep
        Inputs = str_replace_all(Inputs, "; -", ";"),
        Outputs = str_replace_all(Outputs, "; -", ";")
    )

# 4. Save
# ------------------------------------------------------------------------------
out_dir <- file.path(DIR_RESULTS, "meta")
ensure_dir(out_dir)
out_file <- file.path(out_dir, "script_and_figure_index.csv")

write_csv(index, out_file)
msg("Saved index to %s", out_file)

# Print preview
print(index %>% select(Script, Role) %>% head(10))

msg("Done.")

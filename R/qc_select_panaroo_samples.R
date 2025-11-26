# R/qc_select_panaroo_samples.R

suppressPackageStartupMessages({
    library(dplyr)
    library(readr)
    library(stringr)
    library(tools)
    library(fs)
})

# Function to count CDS and Contigs using grep (fast)
get_gff_stats <- function(f) {
    # Size in MB
    sz_mb <- file.size(f) / 1024 / 1024

    # Count CDS
    n_cds <- suppressWarnings(as.integer(system2("grep", c("-c", "'\tCDS\t'", shQuote(f)), stdout = TRUE)))

    # Count Contigs (##sequence-region)
    n_contigs <- suppressWarnings(as.integer(system2("grep", c("-c", "'##sequence-region'", shQuote(f)), stdout = TRUE)))

    # Handle potential failures
    if (length(n_cds) == 0) n_cds <- NA
    if (length(n_contigs) == 0) n_contigs <- NA

    tibble(
        path = f,
        filename = basename(f),
        SampleID = tools::file_path_sans_ext(basename(f)),
        size_mb = sz_mb,
        n_cds = n_cds,
        n_contigs = n_contigs
    )
}

# Main function to select samples
select_panaroo_samples <- function(gff_dirs, metadata_file) {
    # 1. Discover GFF files
    message("Discovering GFF files in: ", paste(gff_dirs, collapse = ", "))
    gffs <- unique(unlist(lapply(
        gff_dirs,
        function(d) {
            if (dir.exists(d)) {
                list.files(d, pattern = "\\.(gff|gff3)$", full.names = TRUE, recursive = TRUE)
            } else {
                character(0)
            }
        }
    )))

    if (length(gffs) == 0) {
        warning("No GFF files found.")
        return(NULL)
    }

    message(sprintf("Found %d GFF files. Analyzing stats...", length(gffs)))

    # 2. Analyze stats
    stats_list <- lapply(gffs, get_gff_stats)
    df <- bind_rows(stats_list)

    # 3. Load Metadata
    if (!file.exists(metadata_file)) {
        stop("Metadata file not found: ", metadata_file)
    }
    message("Loading metadata from: ", metadata_file)
    meta_df <- read_csv(metadata_file, show_col_types = FALSE)

    # Prepare metadata for join
    meta_clean <- meta_df %>%
        mutate(
            SampleID = tools::file_path_sans_ext(file_name),
            Participant_id = as.character(Participant_id),
            Timepoint = as.character(Timepoint)
        ) %>%
        select(SampleID, Participant_id, Timepoint) %>%
        distinct(SampleID, .keep_all = TRUE)

    # 4. Join Metadata
    # We perform an INNER JOIN implicitly by checking if SampleID exists in metadata?
    # Or do we keep all GFFs and just mark them?
    # The user said: "Ensure that only assemblies present in the QC-selected set are used."
    # And "The script now pulls Timepoint... and Participant_id directly from assembly_info_update1.csv."
    # So if it's not in metadata, it shouldn't be in Panaroo (probably).

    df <- df %>%
        left_join(meta_clean, by = "SampleID")

    # 5. Apply Selection Logic
    # Criteria: Size <= 7MB, CDS <= 6000, Contigs <= 500
    # AND must have valid metadata (Participant_id) - implying it's a known sample

    df <- df %>%
        mutate(
            pass_size = !is.na(size_mb) & size_mb <= 15,
            pass_cds = !is.na(n_cds) & n_cds <= 6000,
            pass_contigs = !is.na(n_contigs) & n_contigs <= 500,
            has_metadata = !is.na(Participant_id),
            Selected = pass_size & pass_cds & pass_contigs & has_metadata,
            Reason = case_when(
                Selected ~ "Selected",
                !has_metadata ~ "No Metadata",
                !pass_size ~ "Size > 15MB",
                !pass_cds ~ "CDS > 6000",
                !pass_contigs ~ "Contigs > 500",
                TRUE ~ "Other"
            )
        )

    # --- SMART SUBSAMPLING (Limit to 250) ---
    # If we have too many samples, we prioritize diversity:
    # 1. Keep at least 1 sample per Participant.
    # 2. Then fill up to 250 with additional samples from participants.

    MAX_SAMPLES <- 300
    selected_df <- df %>% filter(Selected)

    if (nrow(selected_df) > MAX_SAMPLES) {
        message(sprintf("⚠️ Too many samples selected (%d). Smart subsampling to %d...", nrow(selected_df), MAX_SAMPLES))

        # Rank samples within each participant (randomly or by quality?)
        # Let's just take the first ones found (usually timepoints)
        ranked <- selected_df %>%
            group_by(Participant_id) %>%
            mutate(rank = row_number()) %>%
            ungroup() %>%
            arrange(rank, Participant_id) %>% # Pick 1st from all, then 2nd from all...
            head(n = MAX_SAMPLES)

        # Update the main df to unselect those that didn't make the cut
        kept_ids <- ranked$SampleID

        df <- df %>%
            mutate(
                Selected = Selected & (SampleID %in% kept_ids),
                Reason = ifelse(Selected, "Selected", ifelse(Reason == "Selected", "Subsampled (Limit 250)", Reason))
            )

        message(sprintf("✓ Reduced to %d samples (prioritizing unique participants).", sum(df$Selected)))
    }

    # Write summary
    message("Class of df: ", paste(class(df), collapse = ", "))
    utils::write.csv(df, "results/wgs/panaroo_selection_summary.csv", row.names = FALSE)

    # Return the full dataframe so the caller can filter/inspect it
    return(df)
}

# Load libraries
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(stringr))
# Helper for lists
`%||%` <- function(a, b) if (!is.null(a)) a else b

# Paths
status_path <- "results/class_inputs_full.csv" # CORRECTED PATH
pairwise_path <- "results/pairwise_stats.csv"
mlst_path <- "results/mlst/mlst_with_meta.csv" # CORRECTED PATH

# 1. Load Data
message("Loading data...")
status_df <- read_csv(status_path, show_col_types = FALSE)
pairwise_df <- read_csv(pairwise_path, show_col_types = FALSE)
mlst_df <- read_csv(mlst_path, show_col_types = FALSE)

# 2. Define Populations & Completeness
# Rule: Must have Status + Isolate ID + MLST
target_status <- c("ASB", "UTI")

# Check column names
# status_df has 'isolate_ID'
# mlst_df has 'Isolate_ID'
# Harmonize
status_df <- status_df %>% rename(isolate_id = isolate_ID, pid = Participant_id)
mlst_df <- mlst_df %>% rename(isolate_id = Isolate_ID)

complete_df <- status_df %>%
    filter(Infection_Status %in% target_status) %>%
    filter(!is.na(isolate_id)) %>%
    # Check if isolate has MLST. Inner join ensures MLST exists.
    inner_join(mlst_df %>% select(isolate_id, ST, scheme, everything()), by = "isolate_id") %>%
    # Keep only necessary columns
    select(pid, timepoint = tp_num, label_timepoint = tp_lab, status = Infection_Status, isolate_id, ST, everything())

# Report Population N
# Completeness criteria: Clinical (Status) + Isolate + Genomic (MLST joined)
n_participants <- n_distinct(complete_df$pid)
n_asb_episodes <- sum(complete_df$status == "ASB")
n_uti_episodes <- sum(complete_df$status == "UTI")

cat("\n--- POPULATION SUMMARY ---\n")
cat("Participants (Complete Data):", n_participants, "\n")
cat("ASB Episodes:", n_asb_episodes, "\n")
cat("UTI Episodes:", n_uti_episodes, "\n")

# 3. Comparisons (ASB vs UTI)
# We need pairs within same participant: ASB -> UTI or UTI -> ASB
# Logic: Find participants with BOTH status types.
pids_with_both <- complete_df %>%
    group_by(pid) %>%
    summarize(
        has_asb = "ASB" %in% status,
        has_uti = "UTI" %in% status
    ) %>%
    filter(has_asb & has_uti) %>%
    pull(pid)

cat("Participants with ASB and UTI episodes:", length(pids_with_both), "\n")

# For these PIDs, we look at their pairs.
# Use pairwise_stats.csv.
# We need to map pairwise records to these PIDs and Statuses.
# Pairwise has 'Isolate_ID' (Query) and 'path_B' (Ref).
# We need to map path_B to an Isolate ID or assume matching naming convention.
# However, usually pairwise includes Isolate_ID for both or filename.

# Let's try to map filenames in pairwise_df to status.
# mlst_df$file_name maps to isolate_id.
iso_lookup <- complete_df %>% select(isolate_id, pid, status, ST, timepoint)

# Filter pairwise for relevant participants
pairs_filtered <- pairwise_df %>%
    filter(Participant_id %in% pids_with_both) %>%
    rename(isolate_1_id = Isolate_ID)

# Now we need isolate_2_id.
# `pairwise_df` column `path_B` usually ends with the filename.
# Let's extract filename from path_B
pairs_filtered$file_name_2 <- basename(pairs_filtered$path_B)

# Now map file_name_2 to Isolate ID using mlst_df (which has file_name and Isolate_ID)
file_to_id <- mlst_df %>%
    select(file_name, isolate_id) %>%
    distinct()

pairs_mapped <- pairs_filtered %>%
    left_join(file_to_id, by = c("file_name_2" = "file_name")) %>%
    rename(isolate_2_id = isolate_id)

# Now add status info
pairs_status <- pairs_mapped %>%
    inner_join(iso_lookup, by = c("isolate_1_id" = "isolate_id"), suffix = c("_1", "_1_meta")) %>%
    inner_join(iso_lookup, by = c("isolate_2_id" = "isolate_id", "pid" = "pid"), suffix = c("_1", "_2"))
# Note: joining by PID ensures within-patient

# Filter for ASB vs UTI
# Valid transitions: ASB->UTI or UTI->ASB
# And we generally only care about distinct timepoints.
pairs_asb_uti <- pairs_status %>%
    filter(status_1 != status_2) %>%
    filter(status_1 %in% c("ASB", "UTI")) %>%
    filter(status_2 %in% c("ASB", "UTI"))

# Limit to unique pairs (A-B is same as B-A).
# Create a sorted key
pairs_asb_uti <- pairs_asb_uti %>%
    rowwise() %>%
    mutate(pair_key = paste(sort(c(isolate_1_id, isolate_2_id)), collapse = "-")) %>%
    ungroup() %>%
    distinct(pair_key, .keep_all = TRUE)

cat("\n--- PAIRWISE COMPARISONS (ASB <-> UTI) ---\n")
cat("Total ASB-UTI pairs analyzed:", nrow(pairs_asb_uti), "\n")

# Similarity
# Threshold: 10 SNPs.
# Columns: TotalSnpCnt
snp_threshold <- 10

pairs_scored <- pairs_asb_uti %>%
    mutate(
        is_similar = case_when(
            !is.na(TotalSnpCnt) & TotalSnpCnt <= snp_threshold ~ TRUE,
            !is.na(TotalSnpCnt) & TotalSnpCnt > snp_threshold ~ FALSE,
            # If SNP is NA, use ST
            is.na(TotalSnpCnt) & (ST_1 == ST_2) ~ TRUE,
            is.na(TotalSnpCnt) & (ST_1 != ST_2) ~ FALSE,
            TRUE ~ NA
        )
    ) %>%
    filter(!is.na(is_similar))

n_similar <- sum(pairs_scored$is_similar)
pct_similar <- mean(pairs_scored$is_similar) * 100

cat("Similar E. coli pairs:", n_similar, "(", round(pct_similar, 1), "%)\n")
cat("Dissimilar E. coli pairs:", nrow(pairs_scored) - n_similar, "\n")

# Median SNPs
cat("Median SNPs (Similar):", median(pairs_scored$TotalSnpCnt[pairs_scored$is_similar], na.rm = TRUE), "\n")
cat("Median SNPs (Dissimilar):", median(pairs_scored$TotalSnpCnt[!pairs_scored$is_similar], na.rm = TRUE), "\n")


# 4. Top E. coli Types
cat("\n--- TOP E. COLI TYPES ---\n")
# Count by Participant-Episode (avoid double counting if duplicates exist)
top_asb <- complete_df %>%
    filter(status == "ASB") %>%
    count(ST) %>%
    arrange(desc(n)) %>%
    head(5)
top_uti <- complete_df %>%
    filter(status == "UTI") %>%
    count(ST) %>%
    arrange(desc(n)) %>%
    head(5)
print("Top ASB STs:")
print(top_asb)
print("Top UTI STs:")
print(top_uti)

# 5. Allele Changes
cat("\n--- ALLELE CHANGES ---\n")
# Filter for Similar (=Persistent) pairs
persistent <- pairs_scored %>% filter(is_similar == TRUE)

# Identify allele columns. In mlst_with_meta, they are between 'scheme' and 'has_new_allele' usually.
# Or we can just grep logical names.
# Based on head: dinb,icda,pabb,polb,putp,trpa,trpb,uida (Likely Pasteur scheme?)
# I'll just exclude known columns.
known_cols <- c(
    "isolate_id", "ST", "file_name", "full_path", "assembler", "num_contigs",
    "total_bases", "gc_content", "pid", "Participant_id", "UTI_Label",
    "Collection_Date", "Spec", "Obj", "Organism", "Beoord", "CFU_Count",
    "Archive", "Timepoint", "Population", "UWI#", "Urine collection method",
    "Batch", "found", "file", "scheme", "has_new_allele", "ambiguous_call",
    "file_name_meta", "full_path_meta", "assembler_meta", "num_contigs_meta",
    "total_bases_meta", "gc_content_meta", "Participant_id_meta",
    "UTI_Label_meta", "Collection_Date_meta", "Spec_meta", "Obj_meta",
    "Organism_meta", "Beoord_meta", "CFU_Count_meta", "Archive_meta",
    "Timepoint_meta", "Population_meta", "UWI#_meta", "Urine collection method_meta",
    "Batch_meta"
)

# We select columns from mlst_df that are NOT in these.
allele_cols <- setdiff(colnames(mlst_df), known_cols)
# Remove suffix matching
allele_cols <- allele_cols[!grepl("_meta$", allele_cols)]
# Also remove id cols if any remain
allele_cols <- allele_cols[!allele_cols %in% c("isolate_id", "ST")]

cat("Identified allele columns:", paste(allele_cols, collapse = ", "), "\n")

changes_list <- list()

if (nrow(persistent) > 0) {
    for (i in 1:nrow(persistent)) {
        id1 <- persistent$isolate_1_id[i]
        id2 <- persistent$isolate_2_id[i]

        # Get rows
        r1 <- mlst_df %>% filter(isolate_id == id1)
        r2 <- mlst_df %>% filter(isolate_id == id2)

        if (nrow(r1) == 1 && nrow(r2) == 1) {
            for (col in allele_cols) {
                if (as.character(r1[[col]]) != as.character(r2[[col]])) {
                    k <- paste(col, r1[[col]], "->", r2[[col]])
                    changes_list[[k]] <- (changes_list[[k]] %||% 0) + 1
                }
            }
        }
    }
}

# Print top changes
changes_vec <- unlist(changes_list)
if (length(changes_vec) > 0) {
    changes_df <- data.frame(change = names(changes_vec), count = changes_vec) %>% arrange(desc(count))
    print(head(changes_df, 10))
} else {
    cat("No allele changes found in persistent strains.\n")
}

library(dplyr)
library(readr)
library(stringr)

# 1. Investigate ST- (Untypeable)
message("Loading MLST data...")
mlst <- read_csv("results/mlst/mlst_with_meta.csv", show_col_types = FALSE)

st_minus <- mlst %>% filter(ST == "-")

message(paste("Total ST- isolates:", nrow(st_minus)))

# Analyze reasons
st_minus_analysis <- st_minus %>%
    transmute(
        Isolate_ID,
        file_name,
        has_new_allele,
        ambiguous_call,
        # Check for missing alleles (NA or 0 or -) in standard 7 genes
        # Assuming columns are: adk, fumC, gyrB, icd, mdh, purA, recA
        # (Note: column names in file might be 'icda', 'pabb' etc from mlst software schema usually 'ecoli')
        # Let's check the schema specific columns. Based on previous view_file, they are:
        # dinb, icda, pabb, polb, putp, trpa, trpb, uida ... wait, previous view showed odd columns.
        # Let's inspect the `scheme` column and then look for gene columns.
        scheme
    )

# Let's look at the gene columns from the file structure seen previously:
# "dinb,icda,pabb,polb,putp,trpa,trpb,uida" -> These are Achtman 7 gene MLST for E. coli?
# wait, standard Achtman is: adk, fumC, gyrB, icd, mdh, purA, recA.
# The columns in line 1 of mlst_with_meta.csv were:
# "found,file,scheme,dinb,icda,pabb,polb,putp,trpa,trpb,uida"
# This looks like the Pasteur scheme (dinB, icdA, pabB, polB, putP, trpA, trpB, uidA).
# But usually verified stats mentioned lineages like ST43, ST131 etc which are usually Achtman.
# Let's check if there are other columns or if 'scheme' says 'ecoli' (Achtman) vs 'ecoli_2' (Pasteur).
# Actually, let's just dump the columns of st_minus_analysis to be safe.

message("Breakdown of ST- reasons:")
print(table(st_minus$has_new_allele, st_minus$ambiguous_call, useNA = "ifany"))

# detailed breakdown
detailed_st_minus <- st_minus %>%
    select(Isolate_ID, has_new_allele, ambiguous_call, contains("~"), contains("?")) # Look for approx or uncertain alleles

print(head(detailed_st_minus))

# 2. Investigate Unlinked Episodes
message("\nLoading Clinical Data...")
clinical <- read_csv("results/clinical/status_map.csv", show_col_types = FALSE) %>%
    filter(cfu_recorded_any == TRUE) # Only consider those that were "valid" episodes clinically

message(paste("Total Clinical Episodes (Any CFU):", nrow(clinical)))

# Load Linkage Map
# We need to see which of these clinical episodes made it into the genomic dataset
# The updated abstract said: "K >= 2: 258 episodes (ST Analysis: 181 episodes linked)"
# So we are looking for the missing (258 - 181) = 77 episodes.

# Re-create the verified cohort logic to be exact
# Group by participant, count TIMEPOINTS
cohort_2tp <- clinical %>%
    group_by(Participant_id) %>%
    filter(n_distinct(Timepoint) >= 2) %>%
    ungroup()

message(paste("Clinical Episodes (>=2 TPs):", nrow(cohort_2tp)))

# Now check linkage
# Logic from get_stratified_stats.R:
# class_inputs <- read_csv("results/class_inputs_full.csv") joined to cohort
class_inputs <- read_csv("results/class_inputs_full.csv", show_col_types = FALSE) %>%
    mutate(Participant_id = as.character(Participant_id))

# Link clinical to isolate ID
clinical_with_isolate <- cohort_2tp %>%
    left_join(class_inputs, by = c("Participant_id", "Timepoint" = "tp_lab"))

# Identify unlinked
unlinked_episodes <- clinical_with_isolate %>%
    filter(is.na(isolate_ID))

message(paste("Episodes missing isolate_ID in class_inputs:", nrow(unlinked_episodes)))

# Identify linked but missing MLST
# i.e. we have an isolate_ID but it's not in mlst_df
linked_episodes <- clinical_with_isolate %>%
    filter(!is.na(isolate_ID))

missing_mlst <- linked_episodes %>%
    anti_join(mlst, by = c("isolate_ID" = "Isolate_ID"))

message(paste("Episodes with isolate_ID but missing from MLST file:", nrow(missing_mlst)))

# Output distinctive reasons
message("\n--- Summary of Missing Links ---")
if (nrow(unlinked_episodes) > 0) {
    message("Top reasons for missing isolate_ID (checking notes/statuses):")
    # Check if there are columns indicating sample status in clinical file
    print(head(unlinked_episodes))
}

if (nrow(missing_mlst) > 0) {
    message("Isolate IDs missing MLST:")
    print(head(missing_mlst$isolate_ID))
}

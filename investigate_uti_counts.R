# Load libraries
library(dplyr)
library(readr)

df <- read_csv("results/clinical/status_map.csv", show_col_types = FALSE)

# Filter UTIs
utis <- df %>% filter(Infection_Status == "UTI")

cat("Total UTIs (including 'Still to be linked'):", nrow(utis), "\n")

# Remove 'Still to be linked' to get the "19" count likely responsible for the 'raw' count
utis_clean <- utis %>% filter(Participant_id != "Still to be linked")
cat("Total UTIs (valid IDs):", nrow(utis_clean), "\n")

# Now print ALL columns for these 19 to find the 2 that might be "missing" something
# specific columns of interest
cols_to_show <- c("Participant_id", "Timepoint", "cfu_recorded_any", "cfu_ge_1e5_any", "beoord_plus3_any", "raw_CFU_examples", "raw_BEO_examples")

cat("\n=== Detailed UTI Data ===\n")
print(utis_clean %>% select(any_of(cols_to_show)) %>% print(n = Inf))

# Check for specific missingness missingness
cat("\n=== Potential Exclusions ===\n")

# Hypothesis 1: cfu_ge_1e5_any is FALSE
excluded_h1 <- utis_clean %>% filter(cfu_ge_1e5_any == FALSE | is.na(cfu_ge_1e5_any))
cat("Count where cfu_ge_1e5_any is FALSE/NA:", nrow(excluded_h1), "\n")
if (nrow(excluded_h1) > 0) print(excluded_h1 %>% select(Participant_id, any_of(cols_to_show)))

# Hypothesis 2: beoord_plus3_any is FALSE
excluded_h2 <- utis_clean %>% filter(beoord_plus3_any == FALSE | is.na(beoord_plus3_any))
cat("Count where beoord_plus3_any is FALSE/NA:", nrow(excluded_h2), "\n")
if (nrow(excluded_h2) > 0) print(excluded_h2 %>% select(Participant_id, any_of(cols_to_show)))

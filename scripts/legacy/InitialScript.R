# Load necessary libraries
library(dplyr)    # For %>% pipe, mutate(), and case_when()
library(stringr)  # For str_extract() and str_detect()
library(Biostrings) #For readDNAStringSet(),

# Define the directory path
dir_path <- "/Users/Aamir/Desktop/rUTIs/ont-yellow-routine-fastas"

# List all .fasta files in the directory
all_fasta_paths <- list.files(
  path       = dir_path, 
  pattern    = "\\.fasta$",
  full.names = TRUE
)

# Create a data frame with assembly information
assembly_info <- tibble(
  full_path = all_fasta_paths,
  file_name = basename(all_fasta_paths)
) %>%
  mutate(
    # Extract Isolate ID (e.g., "24110099601-1") has to include the -1 
    Isolate_ID = str_extract(sapply(strsplit(file_name, "_"), `[`, 3), "[0-9A-Za-z]+-[0-9]+"),
    # Identify the assembler type
    assembler = case_when(
      str_detect(file_name, "flye")        ~ "flye",
      str_detect(file_name, "longcycler") ~ "longcycler",
      TRUE                                ~ "unknown"
    )
  )

# View the result
assembly_info

# Load both batch csv files so we can attach both together (i have converted the xlsx files to csvs be aware no additional 
# info has been added)
# eg: the unused sample in batch 2/ red study numbers etc
# had to add a new mutate as.charachter because one of the csvs were loading participant IDs as charachters and the other 
# loaded these as numbers
# Convert mismatched columns to character in batch1
# Didnt mutate_all(as.charachter) because wanted to maintain some of these integer types in case of numeric analysis later ? is this even required 
batch1<-read.csv("batch1.csv")
batch2<-read.csv("batch2.csv")

batch1 <- batch1 %>% 
  mutate(
    Participant_id = as.character(Participant_id),
    Spec = as.character(Spec),
    UWI. = as.character(UWI.),
    S.S.Dysuria = as.character(S.S.Dysuria),
    S.S.Urgency = as.character(S.S.Urgency),
    S.S.Frequency = as.character(S.S.Frequency),
    S.S.Incontinence = as.character(S.S.Incontinence),
    S.S.Pus = as.character(S.S.Pus),
    S.S.Flank_pain = as.character(S.S.Flank_pain),
    S.S.Suprapubic.pain = as.character(S.S.Suprapubic.pain),
    S.S.Fever = as.character(S.S.Fever),
    S.S.Chills = as.character(S.S.Chills),
    S.S.Delirium = as.character(S.S.Delirium),
    S.S.Other = as.character(S.S.Other),
    No.S.S = as.character(No.S.S)
  )

# Convert mismatched columns to character in batch2
batch2 <- batch2 %>% 
  mutate(
    Participant_id = as.character(Participant_id),
    Spec = as.character(Spec),
    UWI. = as.character(UWI.),
    S.S.Dysuria = as.character(S.S.Dysuria),
    S.S.Urgency = as.character(S.S.Urgency),
    S.S.Frequency = as.character(S.S.Frequency),
    S.S.Incontinence = as.character(S.S.Incontinence),
    S.S.Pus = as.character(S.S.Pus),
    S.S.Flank_pain = as.character(S.S.Flank_pain),
    S.S.Suprapubic.pain = as.character(S.S.Suprapubic.pain),
    S.S.Fever = as.character(S.S.Fever),
    S.S.Chills = as.character(S.S.Chills),
    S.S.Delirium = as.character(S.S.Delirium),
    S.S.Other = as.character(S.S.Other),
    No.S.S = as.character(No.S.S)
  )
# Check data types for batch1 in order to make above changes
str(batch1)
str(batch2)

# Combine batch1 and batch2 into one metadata data frame
metadata <- bind_rows(batch1, batch2)

# Merge with assembly_info and make a new dataframe so we know if something goes wrong the previous step is protected
assembly_info_update1 <- assembly_info %>% 
  left_join(metadata, by = c("Isolate_ID" = "isolate_ID"))
# Remove assembly_info_update1[1] (the long line for where the fasta is stored - only needed this to extract the exact isolate_ID)
assembly_info_update1 <- assembly_info_update1 %>%
  select(-full_path)
# Check for missing metadata
#If this alerts check Isolate_ID step - most probably an improper seperating of the Isolate ID from the file name 
missing_metadata <- assembly_info_update1 %>% filter(is.na(Participant_id))


if (nrow(missing_metadata) > 0) {
  warning(paste("There are", nrow(missing_metadata), "assemblies without metadata."))
} else {
  message("All assemblies have corresponding metadata.")
}




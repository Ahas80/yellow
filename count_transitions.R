suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(readr))

status_df <- read_csv("results/class_inputs_full.csv", show_col_types = FALSE)

# Transitions
pids_asb_uti <- status_df %>%
    group_by(Participant_id) %>%
    summarize(
        has_asb = any(Infection_Status == "ASB"),
        has_uti = any(Infection_Status == "UTI"),
        has_neg = any(Infection_Status == "Negative")
    )

n_asb_uti <- sum(pids_asb_uti$has_asb & pids_asb_uti$has_uti)
n_neg_uti <- sum(pids_asb_uti$has_neg & pids_asb_uti$has_uti)

cat("Participants with ASB + UTI episodes:", n_asb_uti, "\n")
cat("Participants with Negative + UTI episodes:", n_neg_uti, "\n")

# Check if they overlap (ASB+Negative+UTI)
n_all <- sum(pids_asb_uti$has_asb & pids_asb_uti$has_neg & pids_asb_uti$has_uti)
cat("Participants with All 3:", n_all, "\n")

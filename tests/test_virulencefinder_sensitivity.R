#!/usr/bin/env Rscript

# Fast contract tests for the RQ06-RQ08 CGE sensitivity inputs. The final
# analysis repeats these gates and additionally validates all model outputs.

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(jsonlite)
})

root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
output <- file.path(root, "results", "virulencefinder_cge_3_2_1")
effective <- jsonlite::fromJSON(file.path(output, "effective_config.json"), simplifyVector = TRUE)
cfg <- effective$config
cohort <- cfg$cohort

stopifnot(
  as.integer(cfg$analysis$bootstrap_reps) == 10000L,
  as.integer(cfg$analysis$seed) == 20260712L,
  identical(as.integer(cfg$analysis$snp_thresholds), c(10L, 25L, 50L)),
  as.integer(cfg$analysis$primary_snp_threshold) == 25L,
  as.numeric(cfg$analysis$confidence_level) == 0.95
)

manifest <- read_csv(cfg$paths$manifest, show_col_types = FALSE, progress = FALSE) %>%
  mutate(Participant_id = as.character(Participant_id),
         episode_key = paste(Participant_id, tp_lab, sep = "||"))
run_manifest <- read_csv(file.path(output, "run_manifest.csv"), show_col_types = FALSE, progress = FALSE)
adjacent <- read_csv(cfg$paths$adjacent_pairs, show_col_types = FALSE, progress = FALSE) %>%
  mutate(Participant_id = as.character(Participant_id))
direct <- read_csv(cfg$paths$direct_pairs, show_col_types = FALSE, progress = FALSE)
primary_episode <- read_csv(cfg$paths$primary_episode_metrics, show_col_types = FALSE, progress = FALSE) %>%
  mutate(Participant_id = as.character(Participant_id))

stopifnot(
  nrow(manifest) == 532L,
  dplyr::n_distinct(manifest$Participant_id) == 161L,
  nrow(run_manifest) == 1064L,
  all(run_manifest$status == "success"),
  nrow(adjacent) == 371L,
  dplyr::n_distinct(adjacent$Participant_id) == 139L,
  nrow(direct) == 893L,
  nrow(filter(primary_episode, Event_type == "UTI_event")) == 32L,
  dplyr::n_distinct(filter(primary_episode, Event_type == "UTI_event")$Participant_id) == 29L
)

metadata <- c("Participant_id", "tp_lab", "episode_key", "Assembly_ID", "fasta_sha256",
              "UTI_Status", "Event_type", "full_path")
for (profile in c("web_default_id90_cov60", "matched_id80_cov80")) {
  pa <- read_csv(file.path(output, paste0("presence_absence_", profile, ".csv")),
                 show_col_types = FALSE, progress = FALSE)
  genes <- setdiff(names(pa), metadata)
  stopifnot(nrow(pa) == 532L, length(genes) == 681L, !anyDuplicated(pa$episode_key))
}

# Resident bootstrap unit test: a sampled resident contributes every row from
# that resident, and duplicate sampled residents contribute duplicate clusters.
synthetic <- data.frame(Participant_id = c("A", "A", "B", "B", "B"), row_id = 1:5)
ids <- unique(synthetic$Participant_id)
groups <- split(seq_len(nrow(synthetic)), factor(synthetic$Participant_id, levels = ids))
sampled_resident_indices <- c(2L, 2L)
resampled_rows <- unlist(groups[sampled_resident_indices], use.names = FALSE)
stopifnot(identical(resampled_rows, c(3L, 4L, 5L, 3L, 4L, 5L)))

cat("VirulenceFinder sensitivity contract tests: PASS\n")

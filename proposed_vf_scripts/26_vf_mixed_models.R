#!/usr/bin/env Rscript
# ==============================================================================
# 26_vf_mixed_models.R
# ------------------------------------------------------------------------------
# Role: [VF Inferential] - Mixed-effects models for exploratory candidate VF
# features while accounting for repeated measures within participants.
#
# Inputs:
#   - results/vf/vf_episode_dataset.csv
#
# Outputs:
#   - results/vf/vf_feature_glmm_results.csv
#   - results/vf/vf_model_input_summary.csv
#
# Notes:
#   - Keep this script intentionally targeted rather than fitting 164 unstable
#     models by default.
#   - Recommended initial feature set: total VF burden, iro operon burden,
#     capsule burden, selected genes of a priori interest.
#   - Extend only after dataset size improves.
# ==============================================================================

source("00_config.R")

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
  library(purrr)
  library(lme4)
  library(broom.mixed)
})

msg <- function(...) message(sprintf(...))

in_episode <- file.path(DIR_RESULTS, "vf", "vf_episode_dataset.csv")
if (!file.exists(in_episode)) stop("Run 22_vf_build_analysis_dataset.R first.")

d <- read_csv(in_episode, show_col_types = FALSE) %>%
  filter(Infection_Status %in% c("ASB", "UTI")) %>%
  mutate(is_uti = as.integer(Infection_Status == "UTI"))

# Build a conservative candidate feature set
candidate_features <- c("VF_Burden")

iro_genes <- intersect(c("iroB", "iroC", "iroD", "iroE", "iroN"), names(d))
if (length(iro_genes) > 0) d$iro_operon_count <- rowSums(d[, iro_genes, drop = FALSE], na.rm = TRUE)
if ("iro_operon_count" %in% names(d)) candidate_features <- c(candidate_features, "iro_operon_count")

capsule_candidates <- grep("^(kps|cat_Capsule|cat_Capsule_Surface)", names(d), value = TRUE)
if (length(capsule_candidates) > 0) {
  cap_col <- capsule_candidates[1]
  if (cap_col != "capsule_feature") d$capsule_feature <- d[[cap_col]]
  candidate_features <- c(candidate_features, "capsule_feature")
}

model_summary <- tibble(
  n_rows = nrow(d),
  n_participants = dplyr::n_distinct(d$Participant_id),
  uti_rows = sum(d$is_uti == 1, na.rm = TRUE),
  asb_rows = sum(d$is_uti == 0, na.rm = TRUE)
)
write_csv(model_summary, file.path(DIR_RESULTS, "vf", "vf_model_input_summary.csv"))

fit_one <- function(feature) {
  fml <- as.formula(paste0("is_uti ~ ", feature, " + (1|Participant_id)"))
  mod <- tryCatch(
    glmer(fml, data = d, family = binomial(), control = glmerControl(optimizer = "bobyqa")),
    error = function(e) NULL
  )
  if (is.null(mod)) return(tibble(feature = feature, term = feature, estimate = NA_real_, std.error = NA_real_, statistic = NA_real_, p.value = NA_real_, note = "model_failed"))

  broom.mixed::tidy(mod, effects = "fixed") %>%
    filter(term == feature) %>%
    mutate(feature = feature, odds_ratio = exp(estimate), note = NA_character_) %>%
    select(feature, term, estimate, std.error, statistic, p.value, odds_ratio, note)
}

res <- map_dfr(unique(candidate_features), fit_one)
write_csv(res, file.path(DIR_RESULTS, "vf", "vf_feature_glmm_results.csv"))
msg("Done.")

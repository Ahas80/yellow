# Plasmid-mechanism extension for run_rq06_08.R.
#
# This file is sourced by the numbered RQ runner after its authoritative
# 893-pair and 371-adjacent-pair universes have been validated. It deliberately
# adds only the prespecified RQ06, RQ07 and RQ08 endpoints.

messagef("RQ06-RQ08 plasmid extension: validating predicted-plasmid inputs")

path_mob_profiles <- file.path(
  root, "results", "plasmids", "mob_suite", "episode_plasmid_profiles.csv"
)
path_mechanism_profiles <- file.path(
  root, "results", "plasmids", "mob_suite", "episode_mechanism_profiles.csv"
)
path_mob_marker <- file.path(
  root, "results", "plasmids", "mob_suite", "RUN_COMPLETE.txt"
)
for (p in c(path_mob_profiles, path_mechanism_profiles, path_mob_marker)) {
  if (!file.exists(p)) {
    block(
      "plasmid_rq_input_exists", "complete script-09b/29 plasmid layer", p,
      "Plasmid RQ integration cannot use partial or marker-only outputs"
    )
  }
}

required_pair_plasmid <- c(
  "Replicon_Jaccard", "Replicon_Both_Empty",
  "Replicon_Profile_Available", "MOB_Cluster_Jaccard",
  "MOB_Cluster_Both_Empty", "MOB_Profile_Available",
  "Predicted_Plasmid_Count_A", "Predicted_Plasmid_Count_B",
  "MOB_High_Confidence_Profile_Both"
)
require_columns(pair_dat, required_pair_plasmid, "plasmid_pair_metrics")

mechanism_profiles <- read_csv(
  path_mechanism_profiles, show_col_types = FALSE, progress = FALSE
)
require_columns(
  mechanism_profiles,
  c(
    "Participant_id", "tp_lab", "episode_key", "predicted_plasmid_count",
    "plasmid_binned_informative_vf_amr_burden"
  ),
  "episode_mechanism_profiles"
)
mechanism_profiles <- mechanism_profiles %>%
  mutate(
    Participant_id = as.character(Participant_id),
    tp_lab = as.character(tp_lab),
    episode_key = as.character(episode_key),
    predicted_plasmid_count = as.numeric(predicted_plasmid_count),
    plasmid_binned_informative_vf_amr_burden =
      as.numeric(plasmid_binned_informative_vf_amr_burden)
  )
if (
  nrow(mechanism_profiles) != EXPECTED_EPISODES ||
    anyDuplicated(mechanism_profiles$episode_key) ||
    !setequal(mechanism_profiles$episode_key, episode$episode_key)
) {
  block(
    "episode_mechanism_profile_universe", EXPECTED_EPISODES,
    nrow(mechanism_profiles),
    "Mechanism profiles do not equal the authoritative episode universe"
  )
}
record_check(
  "episode_mechanism_profile_universe", "PASS", EXPECTED_EPISODES,
  nrow(mechanism_profiles),
  "Predicted-plasmid endpoints are complete for the selected cohort"
)

adjacent <- adjacent %>%
  mutate(
    any_replicon_profile_change =
      Replicon_Profile_Available & Replicon_Jaccard < 1,
    any_mob_cluster_change =
      MOB_Profile_Available & MOB_Cluster_Jaccard < 1
  )
if (
  sum(adjacent$Replicon_Profile_Available) != EXPECTED_ADJACENT_PAIRS ||
    sum(adjacent$MOB_Profile_Available) != EXPECTED_ADJACENT_PAIRS
) {
  block(
    "adjacent_plasmid_call_availability",
    paste(EXPECTED_ADJACENT_PAIRS, "complete calls for both methods"),
    paste(
      sum(adjacent$Replicon_Profile_Available), "replicon;",
      sum(adjacent$MOB_Profile_Available), "MOB"
    ),
    "Failed plasmid calls must remain missing and cannot enter inference"
  )
}

plasmid_binary_standardized <- function(df, outcome, threshold) {
  z <- df %>%
    mutate(
      .outcome = as.integer(.data[[outcome]]),
      .exposure = as.integer(TotalSNPs <= threshold)
    ) %>%
    filter(
      !is.na(.outcome), !is.na(.exposure), is.finite(days_between)
    )
  if (
    nrow(z) < 20L ||
      n_distinct(z$.outcome) != 2L ||
      n_distinct(z$.exposure) != 2L
  ) {
    return(c(log_or = NA_real_, adjusted_rd = NA_real_))
  }
  fit <- tryCatch(
    suppressWarnings(glm(
      .outcome ~ .exposure + splines::ns(days_between, df = 3),
      family = binomial(), data = z, control = glm.control(maxit = 100)
    )),
    error = function(e) NULL
  )
  if (is.null(fit) || !is.finite(coef(fit)[[".exposure"]])) {
    return(c(log_or = NA_real_, adjusted_rd = NA_real_))
  }
  nd1 <- z
  nd0 <- z
  nd1$.exposure <- 1L
  nd0$.exposure <- 0L
  p1 <- suppressWarnings(predict(fit, newdata = nd1, type = "response"))
  p0 <- suppressWarnings(predict(fit, newdata = nd0, type = "response"))
  c(
    log_or = unname(coef(fit)[[".exposure"]]),
    adjusted_rd = mean(p1 - p0, na.rm = TRUE)
  )
}

plasmid_binary_inference <- function(
    df, outcome, threshold, analysis, seed_offset) {
  outcome_name <- outcome
  z <- df %>%
    filter(
      !is.na(.data[[outcome_name]]), is.finite(days_between),
      !is.na(TotalSNPs)
    )
  outcome_event_count <- sum(z[[outcome_name]] %in% TRUE)
  point <- plasmid_binary_standardized(z, outcome_name, threshold)
  draws <- bootstrap_vector(
    z,
    function(b) plasmid_binary_standardized(b, outcome_name, threshold),
    c("log_or", "adjusted_rd"),
    seed = RQ_SEED + seed_offset
  )
  bind_rows(
    bootstrap_interval(
      exp(draws[, "log_or"]), exp(point[["log_or"]])
    ) %>% mutate(estimand = "odds_ratio", .before = 1),
    bootstrap_interval(
      draws[, "adjusted_rd"], point[["adjusted_rd"]]
    ) %>% mutate(estimand = "adjusted_risk_difference", .before = 1)
  ) %>%
    mutate(
      analysis = analysis,
      outcome = outcome_name,
      exposure = paste0("direct_SNP_le_", threshold),
      n_pairs = nrow(z),
      n_residents = n_distinct(z$Participant_id),
      n_close = sum(z$TotalSNPs <= threshold),
      outcome_events = outcome_event_count,
      model = paste0(
        outcome_name, " ~ direct SNP<=", threshold,
        " + natural spline(days, df=3)"
      ),
      interpretation =
        "assembly-based plasmid predictions; not transfer or transmission evidence",
      .before = 1
    )
}

rq06_plasmid_primary <- bind_rows(
  plasmid_binary_inference(
    adjacent, "any_replicon_profile_change", 25L,
    "primary_replicon_change", 1600L
  ),
  plasmid_binary_inference(
    adjacent, "any_mob_cluster_change", 25L,
    "primary_mob_cluster_change", 1610L
  )
)

rq06_plasmid_thresholds <- map_dfr(
  c("any_replicon_profile_change", "any_mob_cluster_change"),
  function(outcome) {
    map_dfr(c(10L, 25L, 50L), function(threshold) {
      plasmid_binary_inference(
        adjacent, outcome, threshold,
        paste0("threshold_", outcome), 1700L + threshold +
          ifelse(outcome == "any_mob_cluster_change", 100L, 0L)
      )
    })
  }
)

rq06_plasmid_sensitivity <- bind_rows(
  plasmid_binary_inference(
    adjacent %>% filter(!Replicon_Both_Empty),
    "any_replicon_profile_change", 25L,
    "exclude_both_empty_replicon_profiles", 1900L
  ),
  plasmid_binary_inference(
    adjacent %>% filter(!MOB_Cluster_Both_Empty),
    "any_mob_cluster_change", 25L,
    "exclude_both_empty_mob_profiles", 1910L
  ),
  plasmid_binary_inference(
    adjacent %>% filter(MOB_High_Confidence_Profile_Both),
    "any_mob_cluster_change", 25L,
    "high_confidence_mob_profiles_only", 1920L
  )
)

atomic_write_csv(
  rq06_plasmid_primary,
  file.path(dir_rq06, "rq06_plasmid_primary_inference.csv")
)
atomic_write_csv(
  rq06_plasmid_thresholds,
  file.path(dir_rq06, "rq06_plasmid_snp_threshold_sensitivity.csv")
)
atomic_write_csv(
  rq06_plasmid_sensitivity,
  file.path(dir_rq06, "rq06_plasmid_profile_sensitivity.csv")
)
atomic_write_csv(
  adjacent %>%
    select(
      Participant_id, pair_id, key_from, key_to, days_between, TotalSNPs,
      Replicon_Jaccard, Replicon_Both_Empty,
      any_replicon_profile_change, MOB_Cluster_Jaccard,
      MOB_Cluster_Both_Empty, any_mob_cluster_change,
      Predicted_Plasmid_Count_A, Predicted_Plasmid_Count_B,
      MOB_High_Confidence_Profile_Both
    ),
  file.path(dir_rq06, "rq06_adjacent_plasmid_changes_371.csv")
)

# RQ07: only two exploratory endpoints, with resident-aware descriptive
# intervals and no broad feature tests.
event_plasmid <- event_samples %>%
  select(Participant_id, tp_lab, episode_key, collection_date, UTI_Status) %>%
  left_join(
    mechanism_profiles %>%
      select(
        episode_key, predicted_plasmid_count,
        plasmid_binned_informative_vf_amr_burden
      ),
    by = "episode_key", relationship = "one-to-one"
  )
if (anyNA(event_plasmid$predicted_plasmid_count) ||
    anyNA(event_plasmid$plasmid_binned_informative_vf_amr_burden)) {
  block(
    "rq07_plasmid_endpoints_complete", "all UTI-event samples", "missing",
    "At least one UTI-event episode lacks a predicted-plasmid endpoint"
  )
}

descriptive_median_difference <- function(
    df, endpoint, analysis, seed_offset) {
  statistic <- function(z) {
    median(z[[endpoint]][z$UTI_Status == "UTI"], na.rm = TRUE) -
      median(z[[endpoint]][z$UTI_Status == "Not_UTI"], na.rm = TRUE)
  }
  point <- statistic(df)
  draws <- bootstrap_stat(df, statistic, seed = RQ_SEED + seed_offset)
  bootstrap_interval(draws, point) %>%
    mutate(
      analysis = analysis,
      endpoint = endpoint,
      estimand = "median_UTI_minus_median_Not_UTI",
      n_rows = nrow(df),
      n_residents = n_distinct(df$Participant_id),
      n_uti = sum(df$UTI_Status == "UTI"),
      n_not_uti = sum(df$UTI_Status == "Not_UTI"),
      inferential_scope =
        "resident-bootstrap descriptive interval; no hypothesis test",
      .before = 1
    )
}

rq07_plasmid_event <- bind_rows(
  descriptive_median_difference(
    event_plasmid, "predicted_plasmid_count",
    "all_UTI_event_samples", 2000L
  ),
  descriptive_median_difference(
    event_plasmid, "plasmid_binned_informative_vf_amr_burden",
    "all_UTI_event_samples", 2010L
  )
)

paired_plasmid <- nearest_pairs %>%
  select(Participant_id, uti_key, not_uti_key, absolute_days) %>%
  left_join(
    mechanism_profiles %>%
      select(
        uti_key = episode_key,
        uti_plasmid_count = predicted_plasmid_count,
        uti_plasmid_gene_burden =
          plasmid_binned_informative_vf_amr_burden
      ),
    by = "uti_key", relationship = "many-to-one"
  ) %>%
  left_join(
    mechanism_profiles %>%
      select(
        not_uti_key = episode_key,
        not_uti_plasmid_count = predicted_plasmid_count,
        not_uti_plasmid_gene_burden =
          plasmid_binned_informative_vf_amr_burden
      ),
    by = "not_uti_key", relationship = "many-to-one"
  ) %>%
  mutate(
    predicted_plasmid_count_delta =
      uti_plasmid_count - not_uti_plasmid_count,
    plasmid_binned_informative_vf_amr_burden_delta =
      uti_plasmid_gene_burden - not_uti_plasmid_gene_burden
  )

paired_descriptive <- function(df, delta_col, endpoint, seed_offset) {
  statistic <- function(z) median(z[[delta_col]], na.rm = TRUE)
  point <- statistic(df)
  draws <- bootstrap_stat(df, statistic, seed = RQ_SEED + seed_offset)
  bootstrap_interval(draws, point) %>%
    mutate(
      analysis = "nearest_within_resident_pair",
      endpoint = endpoint,
      estimand = "median_within_resident_UTI_minus_Not_UTI",
      n_pairs = nrow(df),
      n_positive = sum(df[[delta_col]] > 0),
      n_zero = sum(df[[delta_col]] == 0),
      n_negative = sum(df[[delta_col]] < 0),
      inferential_scope =
        "resident-paired descriptive interval; no hypothesis test",
      .before = 1
    )
}

rq07_plasmid_paired <- bind_rows(
  paired_descriptive(
    paired_plasmid, "predicted_plasmid_count_delta",
    "predicted_plasmid_count", 2020L
  ),
  paired_descriptive(
    paired_plasmid,
    "plasmid_binned_informative_vf_amr_burden_delta",
    "plasmid_binned_informative_vf_amr_burden", 2030L
  )
)
atomic_write_csv(
  event_plasmid, file.path(dir_rq07, "rq07_plasmid_event_samples.csv")
)
atomic_write_csv(
  rq07_plasmid_event,
  file.path(dir_rq07, "rq07_plasmid_event_descriptive_differences.csv")
)
atomic_write_csv(
  paired_plasmid,
  file.path(dir_rq07, "rq07_nearest_paired_plasmid_values.csv")
)
atomic_write_csv(
  rq07_plasmid_paired,
  file.path(dir_rq07, "rq07_nearest_paired_plasmid_descriptive.csv")
)

# RQ08: plasmid similarities as proxies for direct SNP <=25, and a prespecified
# paired delta-AUC gate before augmenting the current LOOR composite.
rq08_plasmid_auc <- bind_rows(
  auc_inference(
    pair_dat %>% filter(Replicon_Profile_Available),
    "Replicon_Jaccard", "higher_is_close",
    "Corrected replicon Jaccard", 2100L
  ),
  auc_inference(
    pair_dat %>% filter(MOB_Profile_Available),
    "MOB_Cluster_Jaccard", "higher_is_close",
    "MOB primary-cluster Jaccard", 2110L
  )
)
atomic_write_csv(
  rq08_plasmid_auc,
  file.path(dir_rq08, "rq08_plasmid_similarity_auc.csv")
)

fit_loor_plasmid_composite <- function(
    df, plasmid_columns, model_name, reference_threshold = 25L) {
  z <- df %>%
    mutate(
      close_reference = as.integer(TotalSNPs <= reference_threshold),
      same_st_num = as.integer(same_provider_ST)
    ) %>%
    filter(
      provider_ST_both, is.finite(MashDistance),
      is.finite(fresh_vf_jaccard),
      if_all(all_of(plasmid_columns), is.finite)
    )
  ids <- unique(z$Participant_id)
  bind_rows(lapply(ids, function(held_out) {
    train <- z %>% filter(Participant_id != held_out)
    test <- z %>% filter(Participant_id == held_out)
    scale_columns <- c(
      "MashDistance", "fresh_vf_jaccard", plasmid_columns
    )
    scaled_names <- paste0("z_", seq_along(scale_columns))
    for (j in seq_along(scale_columns)) {
      mu <- mean(train[[scale_columns[[j]]]])
      sigma <- sd(train[[scale_columns[[j]]]])
      if (!is.finite(sigma) || sigma == 0) sigma <- 1
      train[[scaled_names[[j]]]] <-
        (train[[scale_columns[[j]]]] - mu) / sigma
      test[[scaled_names[[j]]]] <-
        (test[[scale_columns[[j]]]] - mu) / sigma
    }
    formula <- reformulate(
      c("same_st_num", scaled_names), response = "close_reference"
    )
    fit <- tryCatch(
      suppressWarnings(glm(
        formula, data = train, family = binomial(),
        control = glm.control(maxit = 100)
      )),
      error = function(e) NULL
    )
    pred <- if (is.null(fit)) {
      rep(NA_real_, nrow(test))
    } else {
      suppressWarnings(as.numeric(
        predict(fit, newdata = test, type = "response")
      ))
    }
    test %>%
      transmute(
        Participant_id, pair_id, reference_threshold,
        outcome = close_reference, predicted_probability = pred,
        model = model_name
      )
  }))
}

augmentation_specs <- list(
  list(
    model = "current_plus_replicon_jaccard",
    columns = "Replicon_Jaccard"
  ),
  list(
    model = "current_plus_mob_cluster_jaccard",
    columns = "MOB_Cluster_Jaccard"
  ),
  list(
    model = "current_plus_replicon_and_mob_jaccard",
    columns = c("Replicon_Jaccard", "MOB_Cluster_Jaccard")
  )
)

augmented_predictions <- map_dfr(augmentation_specs, function(spec) {
  fit_loor_plasmid_composite(
    provider_pairs, spec$columns, spec$model, 25L
  )
})
atomic_write_csv(
  augmented_predictions,
  file.path(dir_rq08, "rq08_plasmid_augmented_loor_predictions.csv")
)

paired_auc_delta <- function(aug, model_name, seed_offset) {
  paired <- aug %>%
    filter(model == model_name, is.finite(predicted_probability)) %>%
    rename(augmented_probability = predicted_probability) %>%
    inner_join(
      composite_primary %>%
        select(
          Participant_id, pair_id,
          base_outcome = outcome,
          base_probability = predicted_probability
        ),
      by = c("Participant_id", "pair_id"), relationship = "one-to-one"
    ) %>%
    filter(
      outcome == base_outcome, is.finite(base_probability),
      is.finite(augmented_probability)
    )
  statistic <- function(z) {
    auc_rank(z$outcome, z$augmented_probability) -
      auc_rank(z$outcome, z$base_probability)
  }
  point <- statistic(paired)
  draws <- bootstrap_stat(
    paired, statistic, seed = RQ_SEED + seed_offset
  )
  bootstrap_interval(draws, point) %>%
    mutate(
      model = model_name,
      base_model =
        "current LOOR provider-ST + Mash + fresh-VF composite",
      augmented_auc = auc_rank(
        paired$outcome, paired$augmented_probability
      ),
      base_auc_same_pairs = auc_rank(
        paired$outcome, paired$base_probability
      ),
      n_pairs = nrow(paired),
      n_residents = n_distinct(paired$Participant_id),
      estimand = "paired_delta_auc_augmented_minus_current",
      improvement_supported =
        is.finite(ci_lower) & ci_lower > 0,
      .before = 1
    )
}

rq08_augmentation <- map2_dfr(
  augmentation_specs, seq_along(augmentation_specs),
  function(spec, i) {
    paired_auc_delta(
      augmented_predictions, spec$model, 2200L + 10L * i
    )
  }
)
atomic_write_csv(
  rq08_augmentation,
  file.path(dir_rq08, "rq08_plasmid_augmented_delta_auc.csv")
)

eligible_augmentation <- rq08_augmentation %>%
  filter(improvement_supported) %>%
  arrange(desc(estimate), desc(augmented_auc), model)
selected_model <- if (nrow(eligible_augmentation)) {
  eligible_augmentation$model[[1]]
} else {
  "current_loor_composite_no_plasmid_metric"
}
selected_predictions <- if (
  selected_model == "current_loor_composite_no_plasmid_metric"
) {
  composite_primary %>%
    mutate(model = selected_model)
} else {
  augmented_predictions %>%
    filter(model == selected_model)
}
atomic_write_csv(
  selected_predictions,
  file.path(dir_rq08, "rq08_selected_loor_composite_predictions.csv")
)
atomic_write_csv(
  tibble(
    selected_model = selected_model,
    plasmid_metric_added =
      selected_model != "current_loor_composite_no_plasmid_metric",
    decision_rule =
      "add only when resident-bootstrap paired delta-AUC 95% interval is entirely above zero",
    interpretation =
      "prediction support only; no transfer, transmission or causal inference"
  ),
  file.path(dir_rq08, "rq08_plasmid_composite_selection.csv")
)

plasmid_roc <- bind_rows(
  roc_coordinates(
    pair_dat$TotalSNPs <= 25, pair_dat$Replicon_Jaccard,
    "Corrected replicon Jaccard"
  ),
  roc_coordinates(
    pair_dat$TotalSNPs <= 25, pair_dat$MOB_Cluster_Jaccard,
    "MOB primary-cluster Jaccard"
  )
)
atomic_write_csv(
  plasmid_roc, file.path(dir_rq08, "rq08_plasmid_roc_coordinates.csv")
)
p_rq08_plasmid <- ggplot(
  plasmid_roc, aes(false_positive_rate, sensitivity, colour = predictor)
) +
  geom_abline(
    slope = 1, intercept = 0, linetype = "dashed", colour = "grey65"
  ) +
  geom_path(linewidth = 0.9) +
  coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
  labs(
    title = "Plasmid-profile similarity versus direct SNP closeness",
    subtitle = "Operational reference: direct assembly SNP distance <=25",
    x = "False-positive rate", y = "Sensitivity", colour = "Predictor",
    caption = paste(
      "Assembly-based predicted plasmid context.",
      "Resident-bootstrap AUCs and paired delta-AUCs are tabulated."
    )
  ) +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom")
atomic_ggsave(
  p_rq08_plasmid,
  file.path(dir_rq08, "rq08_plasmid_similarity_roc.png"), 7, 6
)

record_check(
  "rq06_plasmid_outputs", "PASS",
  "371 complete adjacent pairs", nrow(adjacent),
  "Replicon and MOB change outcomes only"
)
record_check(
  "rq07_plasmid_outputs", "PASS",
  "two prespecified exploratory endpoints", 2,
  "Resident-aware descriptive intervals; no broad feature tests"
)
record_check(
  "rq08_plasmid_outputs", "PASS",
  "resident-bootstrap AUC and paired delta-AUC gate", selected_model,
  "Plasmid metrics enter the selected composite only if improvement is supported"
)

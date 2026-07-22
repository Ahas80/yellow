#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(digest)
  library(dplyr)
  library(ggplot2)
  library(readr)
  library(splines)
  library(stringr)
  library(tidyr)
})

project_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
if (!file.exists(file.path(project_root, "00_config.R"))) {
  stop("Run this script from the rUTIs project root.", call. = FALSE)
}
source(file.path(project_root, "R", "pipeline_qc_helpers.R"))

seed <- 20260712L
n_boot <- suppressWarnings(as.integer(Sys.getenv("RQ_BOOTSTRAP_REPS", "10000")))
if (is.na(n_boot) || n_boot < 100L) stop("RQ_BOOTSTRAP_REPS must be at least 100.", call. = FALSE)

out_root <- file.path(project_root, "results", "research_questions")
dir.create(out_root, recursive = TRUE, showWarnings = FALSE)

status_path <- file.path(project_root, "results", "clinical", "status_map.csv")
manifest_path <- file.path(project_root, "results", "qc", "analysis_assembly_manifest.csv")
transition_path <- file.path(project_root, "results", "longitudinal", "longcycler_transitions.csv")
input_paths <- c(status_path, manifest_path, transition_path)

if (any(grepl("Rowenas analysis", input_paths, fixed = TRUE))) {
  stop("Excluded input tree detected.", call. = FALSE)
}
if (any(!file.exists(input_paths))) {
  stop("Missing required input(s): ", paste(input_paths[!file.exists(input_paths)], collapse = ", "), call. = FALSE)
}

atomic_write_csv <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  tmp <- tempfile(pattern = paste0(".", basename(path), "."), tmpdir = dirname(path), fileext = ".csv")
  readr::write_csv(x, tmp, na = "")
  if (!file.rename(tmp, path)) stop("Could not atomically publish ", path, call. = FALSE)
  invisible(path)
}

atomic_write_lines <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  tmp <- tempfile(pattern = paste0(".", basename(path), "."), tmpdir = dirname(path), fileext = ".txt")
  writeLines(x, tmp, useBytes = TRUE)
  if (!file.rename(tmp, path)) stop("Could not atomically publish ", path, call. = FALSE)
  invisible(path)
}

atomic_ggsave <- function(path, plot, width = 7, height = 5, dpi = 300) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  tmp <- tempfile(pattern = paste0(".", tools::file_path_sans_ext(basename(path)), "."),
                  tmpdir = dirname(path), fileext = paste0(".", tools::file_ext(path)))
  ggplot2::ggsave(tmp, plot, width = width, height = height, dpi = dpi)
  if (!file.rename(tmp, path)) stop("Could not atomically publish ", path, call. = FALSE)
  invisible(path)
}

as_flag <- function(x) {
  if (is.logical(x)) return(replace(x, is.na(x), FALSE))
  tolower(trimws(as.character(x))) %in% c("true", "t", "1", "yes", "y")
}

git_head <- function() {
  x <- tryCatch(system2("git", c("rev-parse", "HEAD"), stdout = TRUE, stderr = FALSE),
                error = function(e) NA_character_)
  if (length(x)) x[[1]] else NA_character_
}

input_manifest <- tibble(
  input_role = c("clinical_status", "selected_longcycler_manifest", "direct_adjacent_transitions"),
  path = normalizePath(input_paths, winslash = "/", mustWork = TRUE),
  sha256 = vapply(input_paths, digest::digest, character(1), file = TRUE, algo = "sha256"),
  bytes = as.numeric(file.info(input_paths)$size)
)

status <- read_csv(status_path, show_col_types = FALSE) %>%
  mutate(
    Participant_id = as.character(.data$Participant_id),
    tp_lab = as.character(.data$tp_lab),
    analysis_include_primary = as_flag(.data$analysis_include_primary)
  ) %>%
  filter(.data$analysis_include_primary)

manifest <- read_csv(manifest_path, show_col_types = FALSE) %>%
  mutate(
    Participant_id = as.character(.data$Participant_id),
    tp_lab = as.character(.data$tp_lab)
  )

transitions <- read_csv(transition_path, show_col_types = FALSE) %>%
  mutate(
    Participant_id = as.character(.data$Participant_id),
    tp_from = as.character(.data$tp_from),
    tp_to = as.character(.data$tp_to),
    days_between_samples = as.numeric(.data$days_between_samples),
    TotalSNPs = as.numeric(.data$TotalSNPs),
    Replicon_Jaccard = as.numeric(.data$Replicon_Jaccard),
    Replicon_Both_Empty = as_flag(.data$Replicon_Both_Empty),
    Replicon_Profile_Available =
      as_flag(.data$Replicon_Profile_Available),
    MOB_Cluster_Jaccard = as.numeric(.data$MOB_Cluster_Jaccard),
    MOB_Cluster_Both_Empty = as_flag(.data$MOB_Cluster_Both_Empty),
    MOB_Profile_Available = as_flag(.data$MOB_Profile_Available),
    MOB_High_Confidence_Profile_Both =
      as_flag(.data$MOB_High_Confidence_Profile_Both)
  )

stopifnot(
  nrow(status) == 583L,
  n_distinct(status$Participant_id) == 166L,
  sum(status$UTI_Status == "UTI") == 18L,
  sum(status$UTI_Status == "Not_UTI") == 565L,
  nrow(manifest) == 532L,
  n_distinct(manifest$Participant_id) == 161L,
  !anyDuplicated(manifest[c("Participant_id", "tp_lab")]),
  nrow(transitions) == 371L,
  n_distinct(transitions$Participant_id) == 139L,
  all(!is.na(transitions$TotalSNPs)),
  all(!is.na(transitions$days_between_samples)),
  sum(transitions$TotalSNPs <= 25) == 140L
)
required_plasmid_transition_columns <- c(
  "Replicon_Jaccard", "Replicon_Both_Empty",
  "Replicon_Profile_Available", "MOB_Cluster_Jaccard",
  "MOB_Cluster_Both_Empty", "MOB_Profile_Available",
  "MOB_High_Confidence_Profile_Both"
)
if (
  !all(required_plasmid_transition_columns %in% names(transitions)) ||
    sum(transitions$Replicon_Profile_Available) != 371L ||
    sum(transitions$MOB_Profile_Available) != 371L ||
    anyNA(transitions$Replicon_Jaccard) ||
    anyNA(transitions$MOB_Cluster_Jaccard)
) {
  stop(
    "RQ01 requires complete corrected replicon and MOB profiles for all ",
    "371 adjacent pairs.", call. = FALSE
  )
}

analysis_status <- status %>%
  semi_join(
    manifest %>% distinct(.data$Participant_id, .data$tp_lab),
    by = c("Participant_id", "tp_lab")
  )
stopifnot(
  nrow(analysis_status) == 532L,
  n_distinct(analysis_status$Participant_id) == 161L,
  sum(analysis_status$UTI_Status == "UTI") == 16L,
  sum(analysis_status$UTI_Status == "Not_UTI") == 516L,
  !anyDuplicated(analysis_status[c("Participant_id", "tp_lab")])
)

manifest_assembler <- tolower(if ("assembler" %in% names(manifest)) manifest$assembler else manifest$Assembler)
if (any(manifest_assembler != "longcycler" | is.na(manifest_assembler))) {
  stop("The authoritative analysis manifest is not Longcycler-only.", call. = FALSE)
}

episode_context <- analysis_status %>%
  select(Participant_id, tp_lab, Event_type, UTI_Status)

transitions <- transitions %>%
  left_join(
    episode_context %>% rename(tp_from = tp_lab, event_from = Event_type,
                               verified_status_from = UTI_Status),
    by = c("Participant_id", "tp_from"), relationship = "many-to-one"
  ) %>%
  left_join(
    episode_context %>% rename(tp_to = tp_lab, event_to = Event_type,
                               verified_status_to = UTI_Status),
    by = c("Participant_id", "tp_to"), relationship = "many-to-one"
  ) %>%
  mutate(
    close10 = .data$TotalSNPs <= 10,
    close25 = .data$TotalSNPs <= 25,
    close50 = .data$TotalSNPs <= 50,
    event_interval = .data$event_from == "UTI_event" | .data$event_to == "UTI_event",
    interval_type = case_when(
      .data$event_from == "Routine" & .data$event_to == "Routine" ~ "Routine->Routine",
      .data$event_from == "Routine" & .data$event_to == "UTI_event" ~ "Routine->UTI_event",
      .data$event_from == "UTI_event" & .data$event_to == "Routine" ~ "UTI_event->Routine",
      .data$event_from == "UTI_event" & .data$event_to == "UTI_event" ~ "UTI_event->UTI_event",
      TRUE ~ "Unknown"
    ),
    any_replicon_profile_change =
      .data$Replicon_Profile_Available &
      .data$Replicon_Jaccard < 1,
    any_mob_cluster_change =
      .data$MOB_Profile_Available &
      .data$MOB_Cluster_Jaccard < 1
  )

if (any(is.na(transitions$event_from) | is.na(transitions$event_to))) {
  stop("At least one transition endpoint lacks verified event context.", call. = FALSE)
}
if (any(transitions$status_from != transitions$verified_status_from |
        transitions$status_to != transitions$verified_status_to)) {
  stop("Transition clinical status disagrees with current status_map.csv.", call. = FALSE)
}
stopifnot(
  sum(transitions$interval_type == "Routine->Routine") == 322L,
  sum(transitions$event_interval) == 49L,
  sum(transitions$status_from == "Not_UTI" & transitions$status_to == "UTI") == 9L,
  sum(transitions$status_from == "Not_UTI" & transitions$status_to == "UTI" & transitions$close25) == 5L
)

cluster_weights <- function(data, B, seed_offset = 0L) {
  ids <- unique(data$Participant_id)
  set.seed(seed + seed_offset)
  replicate(B, tabulate(sample.int(length(ids), length(ids), replace = TRUE), nbins = length(ids))) %>%
    matrix(nrow = length(ids), ncol = B, dimnames = list(ids, NULL))
}

weighted_prop_boot <- function(data, outcome, B = n_boot, seed_offset = 0L) {
  ids <- unique(data$Participant_id)
  by_id <- data %>%
    group_by(.data$Participant_id) %>%
    summarise(num = sum(.data[[outcome]] %in% TRUE), den = n(), .groups = "drop") %>%
    slice(match(ids, .data$Participant_id))
  w <- cluster_weights(data, B, seed_offset)
  reps <- colSums(w * by_id$num) / colSums(w * by_id$den)
  tibble(
    estimate = mean(data[[outcome]] %in% TRUE),
    conf_low = unname(quantile(reps, 0.025, na.rm = TRUE)),
    conf_high = unname(quantile(reps, 0.975, na.rm = TRUE)),
    bootstrap_reps = B,
    bootstrap_success = sum(is.finite(reps)),
    bootstrap_failures = B - sum(is.finite(reps))
  )
}

bootstrap_glm_stat <- function(data, fit_stat, B = n_boot, seed_offset = 0L) {
  ids <- unique(data$Participant_id)
  w <- cluster_weights(data, B, seed_offset)
  out <- vector("list", B)
  for (b in seq_len(B)) {
    row_w <- w[match(data$Participant_id, ids), b]
    out[[b]] <- tryCatch(fit_stat(data, row_w), error = function(e) NULL)
  }
  template <- fit_stat(data, rep(1, nrow(data)))
  mat <- matrix(NA_real_, nrow = B, ncol = length(template),
                dimnames = list(NULL, names(template)))
  for (b in seq_len(B)) {
    if (!is.null(out[[b]]) && length(out[[b]]) == ncol(mat)) mat[b, ] <- out[[b]]
  }
  list(point = template, replicates = mat)
}

summarise_boot_vector <- function(point, reps) {
  tibble(
    metric = names(point),
    estimate = as.numeric(point),
    conf_low = apply(reps, 2, quantile, probs = 0.025, na.rm = TRUE),
    conf_high = apply(reps, 2, quantile, probs = 0.975, na.rm = TRUE),
    bootstrap_success = colSums(is.finite(reps)),
    bootstrap_failures = nrow(reps) - colSums(is.finite(reps)),
    bootstrap_reps = nrow(reps)
  )
}

write_common <- function(rq_dir, analysis_name, denominator_rows, denominator_residents,
                         extra = tibble()) {
  atomic_write_csv(input_manifest, file.path(rq_dir, "input_manifest.csv"))
  provenance <- bind_rows(
    tibble(
      field = c("analysis", "git_head", "generated_at", "seed", "bootstrap_reps",
                "denominator_rows", "denominator_residents", "snp_rule"),
      value = c(
        analysis_name,
        git_head(),
        format(Sys.time(), tz = "Europe/Amsterdam", usetz = TRUE),
        as.character(seed),
        as.character(n_boot),
        as.character(denominator_rows),
        as.character(denominator_residents),
        "Direct pairwise DNAdiff <=25 SNPs; operational reference, not biological gold standard"
      )
    ),
    extra
  )
  atomic_write_csv(provenance, file.path(rq_dir, "provenance.csv"))
}

write_analysis_status <- function(rq_dir, research_question, denominator_detail) {
  atomic_write_csv(
    tibble(
      research_question = research_question,
      status = "complete",
      reason = "analysis_completed",
      detail = denominator_detail,
      run_timestamp_utc = format(Sys.time(), tz = "UTC", usetz = TRUE)
    ),
    file.path(rq_dir, "analysis_status.csv")
  )
}

theme_rq <- function() {
  theme_bw(base_size = 11) +
    theme(plot.caption = element_text(hjust = 0, size = 8, colour = "grey35"),
          legend.position = "bottom")
}

# -----------------------------------------------------------------------------
# RQ01 -- Genomic continuity over time
# -----------------------------------------------------------------------------
rq <- file.path(out_root, "RQ01")
dir.create(rq, recursive = TRUE, showWarnings = FALSE)

rq1_primary <- weighted_prop_boot(transitions, "close25", seed_offset = 101L) %>%
  mutate(n_pairs = nrow(transitions), n_residents = n_distinct(transitions$Participant_id),
         n_close = sum(transitions$close25), threshold = 25L, .before = 1)
atomic_write_csv(rq1_primary, file.path(rq, "primary_estimate.csv"))

rq1_thresholds <- bind_rows(lapply(seq_along(c(10L, 25L, 50L)), function(i) {
  threshold <- c(10L, 25L, 50L)[i]
  col <- paste0("close", threshold)
  weighted_prop_boot(transitions, col, seed_offset = 110L + i) %>%
    mutate(threshold = threshold, n_pairs = nrow(transitions), n_close = sum(transitions[[col]]), .before = 1)
}))
atomic_write_csv(rq1_thresholds, file.path(rq, "threshold_sensitivity.csv"))

rq1_subsets <- bind_rows(
  weighted_prop_boot(filter(transitions, interval_type == "Routine->Routine"), "close25", seed_offset = 121L) %>%
    mutate(subset = "routine_to_routine", n_pairs = sum(transitions$interval_type == "Routine->Routine"), .before = 1),
  weighted_prop_boot(filter(transitions, days_between_samples <= 365), "close25", seed_offset = 122L) %>%
    mutate(subset = "interval_le365_days", n_pairs = sum(transitions$days_between_samples <= 365), .before = 1)
)
atomic_write_csv(rq1_subsets, file.path(rq, "subset_sensitivity.csv"))

prediction_days <- c(30, 90, 180, 365)
rq1_fit_stat <- function(d, w) {
  fit <- suppressWarnings(glm(close25 ~ ns(days_between_samples, df = 3),
                              data = d, weights = w, family = binomial()))
  pred <- predict(fit, newdata = data.frame(days_between_samples = prediction_days), type = "response")
  names(pred) <- paste0("day_", prediction_days)
  c(overall = weighted.mean(d$close25, w), pred)
}
rq1_boot <- bootstrap_glm_stat(transitions, rq1_fit_stat, seed_offset = 130L)
rq1_model <- summarise_boot_vector(rq1_boot$point, rq1_boot$replicates)
atomic_write_csv(rq1_model, file.path(rq, "time_model_estimates.csv"))

grid_days <- seq(min(transitions$days_between_samples), max(transitions$days_between_samples), length.out = 200)
rq1_base_fit <- glm(close25 ~ ns(days_between_samples, df = 3), data = transitions, family = binomial())
atomic_write_csv(
  tibble(
    model = "close25 ~ ns(days_between_samples, df = 3)",
    input_rows = nrow(transitions),
    complete_case_rows = sum(complete.cases(transitions[c("close25", "days_between_samples")])),
    point_model_converged = isTRUE(rq1_base_fit$converged),
    bootstrap_reps = n_boot,
    bootstrap_failures_max_across_estimands = max(n_boot - colSums(is.finite(rq1_boot$replicates)))
  ),
  file.path(rq, "model_diagnostics.csv")
)
rq1_curve <- tibble(days_between_samples = grid_days,
                    predicted_probability = predict(rq1_base_fit,
                                                    newdata = data.frame(days_between_samples = grid_days),
                                                    type = "response"))
atomic_write_csv(rq1_curve, file.path(rq, "time_prediction_curve.csv"))

rq1_context <- transitions %>%
  mutate(
    ST_comparable = !is.na(.data$ST_from) & !is.na(.data$ST_to) &
      nzchar(as.character(.data$ST_from)) & nzchar(as.character(.data$ST_to)),
    same_ST = .data$ST_comparable & as.character(.data$ST_from) == as.character(.data$ST_to)
  ) %>%
  group_by(.data$close25) %>%
  summarise(
    n_pairs = n(),
    n_ST_comparable = sum(.data$ST_comparable),
    n_same_ST = sum(.data$same_ST),
    median_Mash_distance = median(.data$MashDistance, na.rm = TRUE),
    q1_Mash_distance = quantile(.data$MashDistance, .25, na.rm = TRUE),
    q3_Mash_distance = quantile(.data$MashDistance, .75, na.rm = TRUE),
    .groups = "drop"
  )
atomic_write_csv(rq1_context, file.path(rq, "st_mash_context.csv"))

rq1_plasmid_fit_stat <- function(d, w, outcome, threshold) {
  failure <- c(odds_ratio = NA_real_, adjusted_risk_difference = NA_real_)
  z <- d %>%
    mutate(
      .outcome = as.integer(.data[[outcome]]),
      .close = as.integer(.data$TotalSNPs <= threshold),
      .weight = as.numeric(w)
    ) %>%
    filter(
      !is.na(.data$.outcome), !is.na(.data$.close),
      is.finite(.data$days_between_samples),
      is.finite(.data$.weight), .data$.weight > 0
    )
  if (
    nrow(z) < 20L ||
      n_distinct(z$.outcome) != 2L ||
      n_distinct(z$.close) != 2L
  ) {
    return(failure)
  }
  fit <- tryCatch(
    suppressWarnings(glm(
      .outcome ~ .close + ns(days_between_samples, df = 3),
      data = z, weights = .weight, family = binomial(),
      control = glm.control(maxit = 100)
    )),
    error = function(e) NULL
  )
  if (
    is.null(fit) || !".close" %in% names(coef(fit)) ||
      !is.finite(coef(fit)[[".close"]])
  ) {
    return(failure)
  }
  close <- z
  distant <- z
  close$.close <- 1L
  distant$.close <- 0L
  p_close <- suppressWarnings(
    predict(fit, newdata = close, type = "response")
  )
  p_distant <- suppressWarnings(
    predict(fit, newdata = distant, type = "response")
  )
  c(
    odds_ratio = exp(unname(coef(fit)[[".close"]])),
    adjusted_risk_difference = weighted.mean(
      p_close - p_distant, z$.weight, na.rm = TRUE
    )
  )
}

rq1_plasmid_inference <- function(
    data, outcome, threshold, analysis, seed_offset) {
  outcome_name <- outcome
  outcome_event_count <- sum(data[[outcome_name]] %in% TRUE)
  boot <- bootstrap_glm_stat(
    data,
    function(d, w) {
      rq1_plasmid_fit_stat(d, w, outcome = outcome, threshold = threshold)
    },
    seed_offset = seed_offset
  )
  summarise_boot_vector(boot$point, boot$replicates) %>%
    mutate(
      analysis = analysis,
      outcome = outcome_name,
      exposure = paste0("direct_SNP_le_", threshold),
      n_pairs = nrow(data),
      n_residents = n_distinct(data$Participant_id),
      n_close = sum(data$TotalSNPs <= threshold),
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

rq1_plasmid_primary <- bind_rows(
  rq1_plasmid_inference(
    transitions, "any_replicon_profile_change", 25L,
    "primary_replicon_change", 1600L
  ),
  rq1_plasmid_inference(
    transitions, "any_mob_cluster_change", 25L,
    "primary_mob_cluster_change", 1610L
  )
)
rq1_plasmid_thresholds <- bind_rows(lapply(
  c("any_replicon_profile_change", "any_mob_cluster_change"),
  function(outcome) {
    bind_rows(lapply(c(10L, 25L, 50L), function(threshold) {
      rq1_plasmid_inference(
        transitions, outcome, threshold,
        paste0("threshold_", outcome),
        1700L + threshold +
          ifelse(outcome == "any_mob_cluster_change", 100L, 0L)
      )
    }))
  }
))
rq1_plasmid_sensitivity <- bind_rows(
  rq1_plasmid_inference(
    filter(transitions, !.data$Replicon_Both_Empty),
    "any_replicon_profile_change", 25L,
    "exclude_both_empty_replicon_profiles", 1900L
  ),
  rq1_plasmid_inference(
    filter(transitions, !.data$MOB_Cluster_Both_Empty),
    "any_mob_cluster_change", 25L,
    "exclude_both_empty_mob_profiles", 1910L
  ),
  rq1_plasmid_inference(
    filter(transitions, .data$MOB_High_Confidence_Profile_Both),
    "any_mob_cluster_change", 25L,
    "high_confidence_mob_profiles_only", 1920L
  )
)
atomic_write_csv(
  rq1_plasmid_primary,
  file.path(rq, "plasmid_change_primary_inference.csv")
)
atomic_write_csv(
  rq1_plasmid_thresholds,
  file.path(rq, "plasmid_change_snp_threshold_sensitivity.csv")
)
atomic_write_csv(
  rq1_plasmid_sensitivity,
  file.path(rq, "plasmid_change_profile_sensitivity.csv")
)

p_rq1 <- ggplot(rq1_curve, aes(.data$days_between_samples, .data$predicted_probability)) +
  geom_rug(data = transitions, aes(x = .data$days_between_samples), inherit.aes = FALSE, alpha = 0.18) +
  geom_line(linewidth = 1, colour = "#1B6CA8") +
  geom_point(data = filter(rq1_model, grepl("^day_", .data$metric)) %>%
               mutate(days_between_samples = as.numeric(str_remove(.data$metric, "day_")),
                      predicted_probability = .data$estimate), size = 2.2) +
  geom_errorbar(data = filter(rq1_model, grepl("^day_", .data$metric)) %>%
                  mutate(days_between_samples = as.numeric(str_remove(.data$metric, "day_")),
                         predicted_probability = .data$estimate),
                aes(ymin = .data$conf_low, ymax = .data$conf_high), width = 8) +
  scale_y_continuous(labels = scales::label_percent(accuracy = 1), limits = c(0, 1)) +
  labs(title = "Operational genomic continuity over elapsed time",
       x = "Days between consecutive observed isolates", y = "Predicted probability of <=25 DNAdiff SNPs",
       caption = "Direct Longcycler endpoint comparisons; spline is descriptive and resident-cluster uncertainty is reported in the accompanying table.") +
  theme_rq()
atomic_ggsave(file.path(rq, "continuity_over_time.png"), p_rq1, width = 7.5, height = 5)

p_rq1_snp <- ggplot(transitions, aes(x = .data$TotalSNPs)) +
  geom_histogram(bins = 60, fill = "#6B7A8F", colour = "white", linewidth = .15) +
  geom_vline(xintercept = 25, linetype = 2, colour = "#A23B3B") +
  scale_x_continuous(trans = scales::pseudo_log_trans(base = 10)) +
  labs(
    title = "Distribution of direct DNAdiff SNP distances",
    x = "DNAdiff SNPs (log scale)", y = "Consecutive isolate pairs",
    caption = "The dashed line marks the operational 25-SNP reference; the distribution is not evidence of a universal biological cutoff."
  ) +
  theme_rq()
atomic_ggsave(file.path(rq, "snp_distribution.png"), p_rq1_snp, width = 7.5, height = 5)

atomic_write_csv(transitions %>%
                   summarise(n_pairs = n(), n_residents = n_distinct(.data$Participant_id),
                             min_days = min(.data$days_between_samples),
                             q1_days = quantile(.data$days_between_samples, .25),
                             median_days = median(.data$days_between_samples),
                             q3_days = quantile(.data$days_between_samples, .75),
                             max_days = max(.data$days_between_samples)),
                 file.path(rq, "denominator.csv"))
atomic_write_lines(c(
  "RQ01 interpretation boundary",
  "The analysis estimates genomic continuity among consecutively recovered, sequenced E. coli isolates.",
  "It does not demonstrate uninterrupted bladder colonisation, clearance, relapse, or a validated biological strain boundary."
), file.path(rq, "interpretation.txt"))
write_common(rq, "RQ01_genomic_continuity_over_time", nrow(transitions), n_distinct(transitions$Participant_id))
write_analysis_status(rq, "RQ01", "371 adjacent direct-pair comparisons from 139 residents")

# -----------------------------------------------------------------------------
# RQ02 -- Routine surveillance versus UTI-event sampling
# -----------------------------------------------------------------------------
rq <- file.path(out_root, "RQ02")
dir.create(rq, recursive = TRUE, showWarnings = FALSE)

rq2_fit_stat <- function(d, w, outcome = "close25") {
  f <- as.formula(paste0(outcome, " ~ event_interval + ns(days_between_samples, df = 3)"))
  fit <- suppressWarnings(glm(f, data = d, weights = w, family = binomial()))
  d0 <- d; d0$event_interval <- FALSE
  d1 <- d; d1$event_interval <- TRUE
  p0 <- weighted.mean(predict(fit, newdata = d0, type = "response"), w)
  p1 <- weighted.mean(predict(fit, newdata = d1, type = "response"), w)
  c(routine_probability = p0, event_probability = p1,
    risk_difference = p1 - p0, risk_ratio = ifelse(p0 > 0, p1 / p0, NA_real_))
}

rq2_results <- list()
for (i in seq_along(c(10L, 25L, 50L))) {
  threshold <- c(10L, 25L, 50L)[i]
  outcome <- paste0("close", threshold)
  boot <- bootstrap_glm_stat(
    transitions,
    function(d, w) rq2_fit_stat(d, w, outcome = outcome),
    seed_offset = 200L + i
  )
  rq2_results[[i]] <- summarise_boot_vector(boot$point, boot$replicates) %>%
    mutate(threshold = threshold, .before = 1)
}
rq2_estimates <- bind_rows(rq2_results)
atomic_write_csv(rq2_estimates, file.path(rq, "adjusted_event_comparison.csv"))

rq2_diagnostics <- bind_rows(lapply(c(10L, 25L, 50L), function(threshold) {
  outcome <- paste0("close", threshold)
  fit <- suppressWarnings(glm(
    as.formula(paste0(outcome, " ~ event_interval + ns(days_between_samples, df = 3)")),
    data = transitions,
    family = binomial()
  ))
  result_rows <- filter(rq2_estimates, .data$threshold == threshold)
  tibble(
    threshold = threshold,
    input_rows = nrow(transitions),
    complete_case_rows = sum(complete.cases(transitions[c(outcome, "event_interval", "days_between_samples")])),
    point_model_converged = isTRUE(fit$converged),
    bootstrap_reps = n_boot,
    bootstrap_failures_max_across_estimands = max(result_rows$bootstrap_failures)
  )
}))
atomic_write_csv(rq2_diagnostics, file.path(rq, "model_diagnostics.csv"))

rq2_excluding_event_event <- filter(transitions, interval_type != "UTI_event->UTI_event")
boot_rq2_sens <- bootstrap_glm_stat(rq2_excluding_event_event, rq2_fit_stat, seed_offset = 210L)
atomic_write_csv(summarise_boot_vector(boot_rq2_sens$point, boot_rq2_sens$replicates),
                 file.path(rq, "exclude_event_to_event_sensitivity.csv"))

rq2_direction_stat <- function(d, w, outcome = "close25") {
  levels_kept <- c("Routine->Routine", "Routine->UTI_event", "UTI_event->Routine")
  d <- d[d$interval_type %in% levels_kept, , drop = FALSE]
  d$interval_type <- factor(d$interval_type, levels = levels_kept)
  fit <- suppressWarnings(glm(
    as.formula(paste0(outcome, " ~ interval_type + ns(days_between_samples, df = 3)")),
    data = d, weights = w, family = binomial()
  ))
  probabilities <- vapply(levels_kept, function(level) {
    nd <- d
    nd$interval_type <- factor(level, levels = levels_kept)
    weighted.mean(predict(fit, newdata = nd, type = "response"), w)
  }, numeric(1))
  names(probabilities) <- paste0("probability_", make.names(levels_kept))
  c(
    probabilities,
    rd_routine_to_event_vs_routine = probabilities[[2]] - probabilities[[1]],
    rd_event_to_routine_vs_routine = probabilities[[3]] - probabilities[[1]]
  )
}

rq2_directional <- bind_rows(lapply(seq_along(c(10L, 25L, 50L)), function(i) {
  threshold <- c(10L, 25L, 50L)[[i]]
  outcome <- paste0("close", threshold)
  boot <- bootstrap_glm_stat(
    rq2_excluding_event_event,
    function(d, w) rq2_direction_stat(d, w, outcome),
    seed_offset = 220L + i
  )
  summarise_boot_vector(boot$point, boot$replicates) %>%
    mutate(threshold = threshold, .before = 1)
}))
atomic_write_csv(rq2_directional, file.path(rq, "directional_interval_sensitivity.csv"))

rq2_counts <- transitions %>%
  count(.data$interval_type, .data$close25, name = "n") %>%
  group_by(.data$interval_type) %>%
  mutate(total = sum(.data$n), percent = 100 * .data$n / .data$total) %>%
  ungroup()
atomic_write_csv(rq2_counts, file.path(rq, "interval_type_counts.csv"))

p_rq2 <- rq2_counts %>% filter(.data$close25) %>%
  ggplot(aes(x = reorder(.data$interval_type, .data$percent), y = .data$percent)) +
  geom_col(fill = "#3A7D44") + coord_flip() +
  labs(title = "Observed genomic continuity by sampling-interval type",
       x = NULL, y = "Pairs at <=25 DNAdiff SNPs (%)",
       caption = "Unadjusted percentages; elapsed-time-adjusted estimates and resident-bootstrap intervals are in adjusted_event_comparison.csv.") +
  theme_rq()
atomic_ggsave(file.path(rq, "continuity_by_interval_type.png"), p_rq2, width = 7, height = 4.8)

atomic_write_lines(c(
  "RQ02 interpretation boundary",
  "This analysis compares sampling contexts, not an effect of clinical UTI.",
  "UTI-event sampling and operational phenotype overlap strongly, so causal language is prohibited."
), file.path(rq, "interpretation.txt"))
write_common(rq, "RQ02_event_interval_vs_routine", nrow(transitions), n_distinct(transitions$Participant_id))
write_analysis_status(rq, "RQ02", "371 adjacent intervals: 322 routine-to-routine and 49 event-involved")

# -----------------------------------------------------------------------------
# RQ03 -- The isolate preceding operational UTI
# -----------------------------------------------------------------------------
rq <- file.path(out_root, "RQ03")
dir.create(rq, recursive = TRUE, showWarnings = FALSE)

rq3_plasmid_path <- file.path(
  project_root, "results", "plasmids", "mob_suite",
  "not_uti_to_uti_plasmid_metrics_9.csv"
)
if (!file.exists(rq3_plasmid_path)) {
  stop(
    "RQ03 requires the validated nine-case predicted-plasmid table from ",
    "numbered Script 29."
  )
}
rq3_plasmid <- read_csv(
  rq3_plasmid_path, show_col_types = FALSE, progress = FALSE
) %>%
  mutate(
    Participant_id = as.character(Participant_id),
    tp_from = normalise_timepoint_preserve_events(tp_from),
    tp_to = normalise_timepoint_preserve_events(tp_to)
  )
if (
  nrow(rq3_plasmid) != 9L ||
    anyDuplicated(rq3_plasmid[c("Participant_id", "tp_from", "tp_to")])
) {
  stop("RQ03 predicted-plasmid table is not exactly nine unique transitions.")
}

rq3 <- transitions %>%
  filter(.data$status_from == "Not_UTI", .data$status_to == "UTI") %>%
  left_join(
    rq3_plasmid,
    by = c("Participant_id", "tp_from", "tp_to"),
    relationship = "one-to-one",
    suffix = c("", "_plasmid")
  ) %>%
  arrange(.data$collection_date_to, .data$days_between_samples, .data$TotalSNPs) %>%
  mutate(Case_ID = sprintf("Case_%02d", row_number()))
stopifnot(
  nrow(rq3) == 9L, sum(rq3$close25) == 5L,
  all(rq3$TotalSNPs == rq3$TotalSNPs_plasmid)
)

bt <- binom.test(sum(rq3$close25), nrow(rq3))
rq3_primary <- tibble(
  threshold = 25L, n_transitions = nrow(rq3), n_close = sum(rq3$close25),
  proportion = unname(bt$estimate), conf_low = bt$conf.int[[1]], conf_high = bt$conf.int[[2]],
  interval_method = "Clopper-Pearson exact"
)
atomic_write_csv(rq3_primary, file.path(rq, "primary_exact_estimate.csv"))

exact_row <- function(x, subset, threshold) {
  b <- binom.test(sum(x), length(x))
  tibble(
    subset = subset,
    threshold = threshold,
    denominator = length(x),
    n_close = sum(x),
    proportion = as.numeric(b$estimate),
    conf_low = as.numeric(b$conf.int[1]),
    conf_high = as.numeric(b$conf.int[2])
  )
}
rq3_sensitivity <- bind_rows(
  lapply(c(10L, 25L, 50L), function(threshold) {
    exact_row(rq3$TotalSNPs <= threshold, "all", threshold)
  }),
  exact_row(rq3$close25[rq3$event_from == "Routine"], "preceding_routine", 25L),
  exact_row(rq3$close25[rq3$days_between_samples <= 180], "interval_le180_days", 25L)
)
atomic_write_csv(rq3_sensitivity, file.path(rq, "sensitivity_estimates.csv"))

rq3_case_table <- rq3 %>%
  transmute(
    .data$Case_ID,
    days_between_samples = .data$days_between_samples,
    ST_concordant = !is.na(.data$ST_from) & !is.na(.data$ST_to) & as.character(.data$ST_from) == as.character(.data$ST_to),
    DNAdiff_SNPs = .data$TotalSNPs, Mash_distance = .data$MashDistance,
    ANI_percent = .data$AvgIdentity,
    pair_interpretation = .data$pair_interpretation,
    legacy_accessory_composite_classification = .data$Classification,
    at_or_below_25_SNPs = .data$close25,
    corrected_replicon_jaccard = .data$replicon_jaccard,
    corrected_replicon_both_empty = .data$replicon_both_empty,
    replicons_gained = coalesce(.data$replicons_gained, ""),
    replicons_lost = coalesce(.data$replicons_lost, ""),
    mob_cluster_jaccard = .data$mob_cluster_jaccard,
    mob_cluster_both_empty = .data$mob_cluster_both_empty,
    mob_clusters_gained = coalesce(.data$mob_clusters_gained, ""),
    mob_clusters_lost = coalesce(.data$mob_clusters_lost, ""),
    predicted_plasmid_count_from = .data$predicted_plasmid_count_from,
    predicted_plasmid_count_to = .data$predicted_plasmid_count_to,
    predicted_plasmid_count_difference =
      .data$predicted_plasmid_count_difference,
    plasmid_binned_vf_genes_gained =
      coalesce(.data$plasmid_binned_vf_genes_gained, ""),
    plasmid_binned_vf_genes_lost =
      coalesce(.data$plasmid_binned_vf_genes_lost, ""),
    plasmid_binned_informative_amr_genes_gained = coalesce(
      .data$plasmid_binned_informative_amr_genes_gained, ""
    ),
    plasmid_binned_informative_amr_genes_lost = coalesce(
      .data$plasmid_binned_informative_amr_genes_lost, ""
    ),
    mob_high_confidence_profiles_both =
      .data$mob_high_confidence_profiles_both,
    interpretation_scope =
      "descriptive assembly-based plasmid predictions; no transfer or causal claim"
  )
atomic_write_csv(rq3_case_table, file.path(rq, "deidentified_case_matrix.csv"))
atomic_write_csv(
  rq3_case_table %>%
    select(
      Case_ID, DNAdiff_SNPs, pair_interpretation,
      starts_with("corrected_replicon"),
      starts_with("replicons_"), starts_with("mob_"),
      starts_with("predicted_plasmid"),
      starts_with("plasmid_binned"), interpretation_scope
    ),
  file.path(rq, "deidentified_plasmid_mechanism_table.csv")
)

p_rq3 <- ggplot(rq3_case_table, aes(x = reorder(.data$Case_ID, .data$DNAdiff_SNPs),
                                    y = .data$DNAdiff_SNPs, fill = .data$at_or_below_25_SNPs)) +
  geom_col() + geom_hline(yintercept = 25, linetype = 2) + coord_flip() +
  scale_y_log10() + scale_fill_manual(values = c(`TRUE` = "#287D3C", `FALSE` = "#A23B3B")) +
  labs(title = "Direct genomic distance in Not_UTI-to-UTI transitions",
       x = NULL, y = "DNAdiff SNPs (log scale)", fill = "<=25 SNPs",
       caption = "Nine de-identified adjacent transitions; the 25-SNP line is an operational reference.") +
  theme_rq()
atomic_ggsave(file.path(rq, "deidentified_case_snp_distances.png"), p_rq3, width = 7, height = 5)

rq3_plasmid_plot <- rq3_case_table %>%
  transmute(
    Case_ID,
    `Replicon profile` = if_else(
      nzchar(replicons_gained) | nzchar(replicons_lost),
      "Changed", "Stable"
    ),
    `MOB cluster profile` = if_else(
      nzchar(mob_clusters_gained) | nzchar(mob_clusters_lost),
      "Changed", "Stable"
    ),
    `Predicted plasmid count` = if_else(
      predicted_plasmid_count_difference != 0,
      "Changed", "Stable"
    ),
    `Plasmid-binned VF` = if_else(
      nzchar(plasmid_binned_vf_genes_gained) |
        nzchar(plasmid_binned_vf_genes_lost),
      "Changed", "Stable"
    ),
    `Plasmid-binned AMR` = if_else(
      nzchar(plasmid_binned_informative_amr_genes_gained) |
        nzchar(plasmid_binned_informative_amr_genes_lost),
      "Changed", "Stable"
    )
  ) %>%
  pivot_longer(-Case_ID, names_to = "mechanism", values_to = "state") %>%
  mutate(
    Case_ID = factor(Case_ID, levels = rev(unique(Case_ID))),
    state = factor(state, levels = c("Stable", "Changed"))
  )
p_rq3_plasmid <- ggplot(
  rq3_plasmid_plot, aes(mechanism, Case_ID, fill = state)
) +
  geom_tile(colour = "white") +
  geom_text(aes(label = state), size = 2.7) +
  scale_fill_manual(values = c(Stable = "#DCE8F2", Changed = "#D55E00")) +
  labs(
    title = "Predicted plasmid mechanism context in nine transitions",
    subtitle = "Deidentified descriptive cases; no regression is fitted",
    x = NULL, y = NULL, fill = "Predicted state",
    caption = paste(
      "MOB bins and same-bin VF/AMR placement are assembly-based predictions;",
      "they do not establish transfer, transmission or causality."
    )
  ) +
  theme_rq() +
  theme(axis.text.x = element_text(angle = 25, hjust = 1))
atomic_ggsave(
  file.path(rq, "deidentified_plasmid_mechanism.png"),
  p_rq3_plasmid, width = 8.5, height = 5.5
)

atomic_write_lines(c(
  "RQ03 interpretation boundary",
  "The result describes whether the later UTI-labelled isolate was operationally close to the preceding isolate.",
  "It does not show that bacteriuria progressed to, caused, or protected against UTI. No regression is fitted to nine cases."
), file.path(rq, "interpretation.txt"))
write_common(rq, "RQ03_pre_UTI_isolate_relatedness", nrow(rq3), n_distinct(rq3$Participant_id),
             tibble(field = "privacy", value = "Publication case table uses generated Case_ID only"))
write_analysis_status(rq, "RQ03", "9 adjacent Not_UTI-to-UTI comparisons from 9 residents")

# -----------------------------------------------------------------------------
# RQ04 -- Robustness of operational UTI ascertainment
# -----------------------------------------------------------------------------
rq <- file.path(out_root, "RQ04")
dir.create(rq, recursive = TRUE, showWarnings = FALSE)

phenotype <- analysis_status %>%
  mutate(
    current_symptom = .data$symptom_compatible_uti %in% TRUE,
    expanded_symptom = case_when(
      .data$catheter_rule == "A_non_indwelling" ~
        (.data$symptom_compatible_uti %in% TRUE | .data$suprapubic_pain_present %in% TRUE),
      .data$catheter_rule == "B_indwelling" ~ .data$symptom_compatible_uti %in% TRUE,
      TRUE ~ FALSE
    ),
    culture_1e3 = .data$cfu_ge_1e3 %in% TRUE,
    culture_1e4 = .data$cfu_ge_1e4 %in% TRUE,
    culture_1e5 = .data$cfu_ge_1e5 %in% TRUE
  )

rule_grid <- tidyr::crossing(
  symptom_definition = c("current", "expanded_suprapubic"),
  cfu_threshold = c("1e3", "1e4", "1e5")
)
rule_values <- list()
for (i in seq_len(nrow(rule_grid))) {
  sx <- if (rule_grid$symptom_definition[[i]] == "current") phenotype$current_symptom else phenotype$expanded_symptom
  cx <- phenotype[[paste0("culture_", rule_grid$cfu_threshold[[i]])]]
  rule_values[[i]] <- sx & cx
}
names(rule_values) <- paste(rule_grid$symptom_definition, rule_grid$cfu_threshold, sep = "_")

rq4_summary <- bind_rows(lapply(seq_len(nrow(rule_grid)), function(i) {
  x <- rule_values[[i]]
  tibble(
    symptom_definition = rule_grid$symptom_definition[[i]],
    cfu_threshold = rule_grid$cfu_threshold[[i]],
    n_uti_episodes = sum(x),
    n_uti_residents = n_distinct(phenotype$Participant_id[x]),
    overlap_with_primary_uti = sum(x & phenotype$UTI_Status == "UTI"),
    newly_uti_vs_primary = sum(x & phenotype$UTI_Status != "UTI"),
    no_longer_uti_vs_primary = sum(!x & phenotype$UTI_Status == "UTI")
  )
}))
stopifnot(identical(rq4_summary$n_uti_episodes, c(16L, 16L, 16L, 18L, 18L, 18L)))
atomic_write_csv(rq4_summary, file.path(rq, "six_rule_summary.csv"))

rq4_reclass <- bind_rows(lapply(names(rule_values), function(rule) {
  tibble(primary_status = phenotype$UTI_Status,
         alternative_status = ifelse(rule_values[[rule]], "UTI", "Not_UTI")) %>%
    count(.data$primary_status, .data$alternative_status, name = "n") %>%
    mutate(rule = rule, .before = 1)
}))
atomic_write_csv(rq4_reclass, file.path(rq, "reclassification_matrix.csv"))

movement_mask <- Reduce(`|`, lapply(rule_values, function(x) x != (phenotype$UTI_Status == "UTI")))
movement_rows <- phenotype[movement_mask, , drop = FALSE] %>%
  arrange(.data$Collection_Date, .data$tp_lab) %>%
  mutate(Case_ID = sprintf("PhenotypeCase_%02d", row_number()))
rq4_movements <- bind_rows(lapply(names(rule_values), function(rule) {
  idx <- match(movement_rows$Episode_ID, phenotype$Episode_ID)
  alternative <- ifelse(rule_values[[rule]][idx], "UTI", "Not_UTI")
  tibble(
    Case_ID = movement_rows$Case_ID,
    rule = rule,
    primary_status = movement_rows$UTI_Status,
    alternative_status = alternative,
    collection_method = movement_rows$urine_collection_method_norm,
    suprapubic_pain_recorded = movement_rows$suprapubic_pain_present,
    cfu_lower_bound = movement_rows$cfu_lower_bound
  ) %>% filter(.data$primary_status != .data$alternative_status)
}))
atomic_write_csv(rq4_movements, file.path(rq, "deidentified_episode_movements.csv"))

rq4_strata <- bind_rows(lapply(names(rule_values), function(rule) {
  phenotype %>%
    mutate(alternative_uti = rule_values[[rule]]) %>%
    count(.data$urine_collection_method_norm, .data$alternative_uti, name = "n") %>%
    mutate(rule = rule, .before = 1)
}))
atomic_write_csv(rq4_strata, file.path(rq, "collection_method_strata.csv"))

p_rq4 <- ggplot(rq4_summary,
                 aes(x = .data$cfu_threshold, y = .data$n_uti_episodes,
                     fill = .data$symptom_definition, group = .data$symptom_definition)) +
  geom_col(position = position_dodge(width = .75), width = .7) +
  geom_text(aes(label = .data$n_uti_episodes), position = position_dodge(width = .75), vjust = -0.35) +
  coord_cartesian(ylim = c(0, max(rq4_summary$n_uti_episodes) + 4)) +
  labs(title = "Operational UTI count under six observable rule combinations",
       x = "Culture threshold", y = "Episodes classified UTI", fill = "Symptom rule",
       caption = "Deterministic sensitivity analysis only; none of the alternatives is labelled the full protocol phenotype.") +
  theme_rq()
atomic_ggsave(file.path(rq, "six_rule_case_counts.png"), p_rq4, width = 7, height = 5)

atomic_write_lines(c(
  "RQ04 interpretation boundary",
  "This deterministic sensitivity analysis measures dependence on observable CFU and symptom choices.",
  "It cannot reconstruct the full protocol phenotype because recent symptom onset and alternative infectious focus are unavailable."
), file.path(rq, "interpretation.txt"))
write_common(rq, "RQ04_operational_phenotype_robustness", nrow(phenotype), n_distinct(phenotype$Participant_id),
             tibble(field = "case_table_privacy", value = "Episode movements use generated PhenotypeCase_ID"))
write_analysis_status(rq, "RQ04", "532 Longcycler-linked episodes from 161 residents")

# -----------------------------------------------------------------------------
# RQ05 -- Selection into the genomic analysis
# -----------------------------------------------------------------------------
rq <- file.path(out_root, "RQ05")
dir.create(rq, recursive = TRUE, showWarnings = FALSE)

selected_keys <- manifest %>% distinct(.data$Participant_id, .data$tp_lab) %>% mutate(selected_wgs = TRUE)
selection <- status %>%
  left_join(selected_keys, by = c("Participant_id", "tp_lab"), relationship = "one-to-one") %>%
  mutate(selected_wgs = coalesce(.data$selected_wgs, FALSE))
stopifnot(sum(selection$selected_wgs) == 532L,
          sum(selection$selected_wgs & selection$UTI_Status == "UTI") == 16L,
          sum(selection$selected_wgs & selection$UTI_Status == "Not_UTI") == 516L)

rq5_point <- function(d, w) {
  p_uti <- weighted.mean(d$selected_wgs[d$UTI_Status == "UTI"], w[d$UTI_Status == "UTI"])
  p_not <- weighted.mean(d$selected_wgs[d$UTI_Status == "Not_UTI"], w[d$UTI_Status == "Not_UTI"])
  c(uti_selection_probability = p_uti, not_uti_selection_probability = p_not,
    risk_difference = p_uti - p_not, risk_ratio = ifelse(p_not > 0, p_uti / p_not, NA_real_))
}
rq5_boot <- bootstrap_glm_stat(selection, rq5_point, seed_offset = 501L)
rq5_estimates <- summarise_boot_vector(rq5_boot$point, rq5_boot$replicates)
atomic_write_csv(rq5_estimates, file.path(rq, "selection_effect_estimates.csv"))

rq5_counts <- selection %>%
  count(.data$UTI_Status, .data$selected_wgs, name = "n") %>%
  group_by(.data$UTI_Status) %>%
  mutate(status_total = sum(.data$n), percent = 100 * .data$n / .data$status_total) %>%
  ungroup()
atomic_write_csv(rq5_counts, file.path(rq, "selection_counts_by_status.csv"))

fisher_table <- table(selection$UTI_Status, selection$selected_wgs)
ft <- fisher.test(fisher_table)
atomic_write_csv(tibble(
  odds_ratio = unname(ft$estimate), conf_low = ft$conf.int[[1]], conf_high = ft$conf.int[[2]],
  p_value = ft$p.value, method = ft$method
), file.path(rq, "fisher_sensitivity.csv"))

strata_vars <- c("Batch", "Event_type", "urine_collection_method_norm")
rq5_strata <- bind_rows(lapply(strata_vars, function(v) {
  selection %>%
    mutate(stratum = as.character(.data[[v]])) %>%
    count(.data$stratum, .data$selected_wgs, name = "n") %>%
    group_by(.data$stratum) %>%
    mutate(total = sum(.data$n), percent = 100 * .data$n / .data$total) %>%
    ungroup() %>%
    mutate(variable = v, .before = 1)
}))
atomic_write_csv(rq5_strata, file.path(rq, "selection_descriptive_strata.csv"))

p_rq5 <- rq5_counts %>% filter(.data$selected_wgs) %>%
  ggplot(aes(x = .data$UTI_Status, y = .data$percent, fill = .data$UTI_Status)) +
  geom_col(width = .65, show.legend = FALSE) +
  geom_text(aes(label = sprintf("%d/%d", .data$n, .data$status_total)), vjust = -0.4) +
  coord_cartesian(ylim = c(0, 105)) +
  labs(title = "Selection of a QC-passing Longcycler genome by operational status",
       x = NULL, y = "Episodes with selected WGS (%)",
       caption = "Descriptive selection audit; resident-bootstrap risk contrasts are reported separately.") +
  theme_rq()
atomic_ggsave(file.path(rq, "wgs_selection_by_status.png"), p_rq5, width = 6.5, height = 5)

atomic_write_lines(c(
  "RQ05 interpretation boundary",
  "This is a non-analytical attrition and selection-quality audit.",
  "The full clinical source cohort appears here only to account for selection into the 532-episode Longcycler-linked analysis cohort.",
  "Batch, event type, and collection method summaries are missingness diagnostics, not causal adjustment models."
), file.path(rq, "interpretation.txt"))
write_common(
  rq,
  "RQ05_non_analytical_selection_attrition_audit",
  nrow(selection),
  n_distinct(selection$Participant_id),
  tibble(field = "analysis_role", value = "attrition/QC only; not an analytical denominator")
)
write_analysis_status(rq, "RQ05", "attrition/QC audit: 583 source episodes to 532 Longcycler-linked episodes")

message("RQ01--RQ05 analyses complete.")

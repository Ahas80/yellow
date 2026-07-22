#!/usr/bin/env Rscript

# Replay the three model-heavy canonical modules with warning conditions emitted
# immediately, aggregate the exact condition messages, and retain a reviewable
# classification. This is an audit utility: because the modules refresh model
# outputs, the canonical final figure pack must be regenerated afterwards.

root <- normalizePath(".", winslash = "/", mustWork = TRUE)
if (!file.exists(file.path(root, "RUN_COMPLETE_ANALYSIS.sh"))) {
  stop("Run scripts/replay_model_warnings.R from the repository root.", call. = FALSE)
}

suppressPackageStartupMessages({
  library(digest)
  library(dplyr)
  library(readr)
  library(tibble)
})

audit_dir <- file.path(root, "results", "figure_audit")
dir.create(audit_dir, recursive = TRUE, showWarnings = FALSE)

modules <- tibble(
  module = c(
    "Focused gene models",
    "Formal genotype-phenotype models",
    "Statistical sensitivity models"
  ),
  script = file.path(root, c(
    "04_gene_breakdown.R",
    "14_genotype_phenotype_model.R",
    "36_statistical_sensitivity_addon.R"
  ))
)
if (any(!file.exists(modules$script))) {
  stop("One or more warning-replay modules are missing.", call. = FALSE)
}

normalise_condition <- function(value) {
  value <- gsub("[\r\n]+", " ", as.character(value))
  value <- gsub("[[:space:]]+", " ", value)
  trimws(value)
}

classify_condition <- function(message) {
  if (grepl("boundary.*singular|singular fit|isSingular", message, ignore.case = TRUE)) {
    return("model_singularity")
  }
  if (grepl(
    "converg|Hessian|eigenvalue|scaled gradient|max\\|grad\\||nearly unidentifiable",
    message, ignore.case = TRUE
  )) {
    return("model_convergence_diagnostic")
  }
  if (grepl(
    "separation|fitted probabilities numerically 0 or 1|algorithm did not converge|infinite",
    message, ignore.case = TRUE
  )) {
    return("separation_or_nonfinite_effect")
  }
  if (grepl("rank[- ]deficien|not full rank|aliased", message, ignore.case = TRUE)) {
    return("rank_deficiency")
  }
  if (grepl("NaNs produced|numerically 0|essentially perfect fit", message, ignore.case = TRUE)) {
    return("numerical_instability")
  }
  if (grepl("deprecated|deprecation|lifecycle", message, ignore.case = TRUE)) {
    return("prohibited_deprecation")
  }
  if (grepl(
    "removed [0-9]+ rows?|missing values?|non[- ]?finite|outside the scale range|manual scale|many-to-many",
    message, ignore.case = TRUE
  )) {
    return("prohibited_plot_or_join_warning")
  }
  "unclassified_requires_review"
}

condition_rows <- list()
module_rows <- list()
condition_index <- 0L

record_condition <- function(module, script, type, condition) {
  condition_index <<- condition_index + 1L
  call_text <- conditionCall(condition)
  condition_rows[[condition_index]] <<- tibble(
    module = module,
    script = script,
    condition_type = type,
    classification = classify_condition(conditionMessage(condition)),
    message = normalise_condition(conditionMessage(condition)),
    call = if (is.null(call_text)) "" else normalise_condition(paste(deparse(call_text), collapse = " "))
  )
}

for (i in seq_len(nrow(modules))) {
  module <- modules$module[[i]]
  script <- modules$script[[i]]
  message("Replaying warnings for ", basename(script), " ...")
  start <- Sys.time()
  error_message <- ""
  status <- tryCatch(
    withCallingHandlers(
      {
        sys.source(script, envir = new.env(parent = globalenv()), chdir = FALSE)
        TRUE
      },
      warning = function(w) {
        record_condition(module, script, "warning", w)
        invokeRestart("muffleWarning")
      },
      message = function(m) {
        text <- conditionMessage(m)
        if (grepl(
          "boundary.*singular|singular fit|isSingular|converg|Hessian|eigenvalue|scaled gradient|nearly unidentifiable",
          text, ignore.case = TRUE
        )) {
          record_condition(module, script, "diagnostic_message", m)
          invokeRestart("muffleMessage")
        }
      }
    ),
    error = function(e) {
      error_message <<- normalise_condition(conditionMessage(e))
      FALSE
    }
  )
  module_rows[[i]] <- tibble(
    module = module,
    script = script,
    script_sha256 = unname(digest::digest(script, algo = "sha256", file = TRUE)),
    completed = isTRUE(status),
    elapsed_seconds = as.numeric(difftime(Sys.time(), start, units = "secs")),
    error = error_message
  )
  if (!isTRUE(status)) break
}

conditions <- if (length(condition_rows)) bind_rows(condition_rows) else tibble(
  module = character(), script = character(), condition_type = character(),
  classification = character(), message = character(), call = character()
)
modules_run <- bind_rows(module_rows)

condition_summary <- conditions %>%
  count(.data$module, .data$script, .data$condition_type, .data$classification,
        .data$message, .data$call, name = "n") %>%
  arrange(.data$module, .data$condition_type, .data$classification, desc(.data$n), .data$message)

write_csv(condition_summary, file.path(audit_dir, "model_warning_replay.csv"))
write_csv(modules_run, file.path(audit_dir, "model_warning_replay_modules.csv"))

unresolved <- condition_summary %>%
  filter(.data$classification %in% c(
    "prohibited_deprecation", "prohibited_plot_or_join_warning", "unclassified_requires_review"
  ))
all_completed <- nrow(modules_run) == nrow(modules) && all(modules_run$completed)

summary_lines <- c(
  "MODEL WARNING REPLAY",
  "====================",
  paste0("Generated: ", format(Sys.time(), tz = "Europe/Amsterdam", usetz = TRUE)),
  paste0("Modules completed: ", sum(modules_run$completed), "/", nrow(modules)),
  paste0("Distinct condition records: ", nrow(condition_summary)),
  paste0("Total condition occurrences: ", sum(condition_summary$n)),
  paste0("Unresolved/prohibited distinct records: ", nrow(unresolved)),
  "",
  "Classification totals:",
  if (nrow(condition_summary)) {
    condition_summary %>%
      group_by(.data$condition_type, .data$classification) %>%
      summarise(n = sum(.data$n), .groups = "drop") %>%
      transmute(line = paste0("- ", .data$condition_type, " / ", .data$classification, ": ", .data$n)) %>%
      pull("line")
  } else {
    "- No warning or model-diagnostic conditions captured."
  },
  "",
  "Interpretation:",
  "- This replay distinguishes model singularity, convergence, separation and numerical-instability diagnostics from plotting, join and deprecation defects.",
  "- Replaying these modules refreshes model outputs; regenerate the final figure pack and its validation/visual-QA derivatives afterwards."
)
writeLines(summary_lines, file.path(audit_dir, "model_warning_replay_summary.txt"), useBytes = TRUE)

if (!all_completed) {
  stop("At least one model-warning replay module failed; see model_warning_replay_modules.csv.", call. = FALSE)
}
if (nrow(unresolved)) {
  stop("Unresolved or prohibited warning classes remain; review model_warning_replay.csv.", call. = FALSE)
}

message("Model-warning replay completed with ", sum(condition_summary$n),
        " classified condition occurrence(s) and no unresolved class.")

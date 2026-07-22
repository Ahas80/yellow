#!/usr/bin/env Rscript

# Consolidate the machine-verifiable figure, visual-QA, and pipeline release
# evidence into one human-readable report.  This script deliberately derives
# status from current evidence; it never creates or repairs pipeline markers.

root <- normalizePath(".", winslash = "/", mustWork = TRUE)
if (!file.exists(file.path(root, "RUN_COMPLETE_ANALYSIS.sh"))) {
  stop("Run this script from the repository root.", call. = FALSE)
}

audit_dir <- file.path(root, "results", "figure_audit")
dir.create(audit_dir, recursive = TRUE, showWarnings = FALSE)
output_path <- file.path(audit_dir, "validation_results.txt")

rel <- function(path) {
  path <- normalizePath(path, winslash = "/", mustWork = FALSE)
  prefix <- paste0(root, "/")
  ifelse(startsWith(path, prefix), substring(path, nchar(prefix) + 1L), path)
}

fmt_time <- function(value) {
  if (!length(value) || is.na(value)) return("unavailable")
  format(value, "%Y-%m-%d %H:%M:%S %Z")
}

file_summary <- function(path) {
  if (!file.exists(path)) return(paste0(rel(path), ": MISSING"))
  info <- file.info(path)
  paste0(
    rel(path), ": ", format(info$size, scientific = FALSE, trim = TRUE),
    " bytes; modified ", fmt_time(info$mtime)
  )
}

read_csv_evidence <- function(path) {
  if (!file.exists(path)) {
    return(list(data = NULL, error = paste0("missing file: ", rel(path))))
  }
  value <- tryCatch(
    utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE),
    error = function(e) e
  )
  if (inherits(value, "error")) {
    return(list(data = NULL, error = conditionMessage(value)))
  }
  list(data = value, error = NULL)
}

as_flag <- function(x) {
  if (is.logical(x)) return(!is.na(x) & x)
  tolower(trimws(as.character(x))) %in% c("true", "t", "1", "yes", "pass", "passed")
}

marker_value <- function(lines, label) {
  hit <- grep(paste0("^", label, "[[:space:]]*:"), lines, value = TRUE)
  if (!length(hit)) return(NA_character_)
  trimws(sub(paste0("^", label, "[[:space:]]*:"), "", hit[[1L]]))
}

nonempty_files <- function(paths) {
  paths <- as.character(paths)
  good_path <- !is.na(paths) & nzchar(trimws(paths))
  exists <- rep(FALSE, length(paths))
  size <- rep(NA_real_, length(paths))
  exists[good_path] <- file.exists(paths[good_path])
  if (any(exists)) size[exists] <- file.info(paths[exists])$size
  good_path & exists & is.finite(size) & size > 0
}

evidence_paths <- list(
  manifest = file.path(audit_dir, "final_figure_manifest.csv"),
  automated = file.path(audit_dir, "automated_validation_checks.csv"),
  automated_summary = file.path(audit_dir, "automated_validation_results.txt"),
  data_checks = file.path(audit_dir, "final_figure_data_checks.csv"),
  visual_qa = file.path(audit_dir, "visual_qa_results.csv"),
  visual_qa_summary = file.path(audit_dir, "visual_qa_generation_summary.txt"),
  visual_review = file.path(audit_dir, "visual_qa_review.csv"),
  model_warning_replay = file.path(audit_dir, "model_warning_replay.csv"),
  model_warning_modules = file.path(audit_dir, "model_warning_replay_modules.csv"),
  model_warning_summary = file.path(audit_dir, "model_warning_replay_summary.txt"),
  palette = file.path(audit_dir, "palette_accessibility_checks.csv"),
  colour_separation = file.path(audit_dir, "operational_status_colour_separation.csv"),
  final_verification = file.path(root, "results", "qc", "longcycler_only_pipeline_verification.csv"),
  final_verification_summary = file.path(root, "results", "qc", "longcycler_only_pipeline_verification.txt"),
  rq_marker = file.path(root, "results", "research_questions", "RUN_COMPLETE.txt"),
  complete_marker = file.path(root, "results", "pipeline", "RUN_COMPLETE.txt"),
  failed_marker = file.path(root, "results", "pipeline", "RUN_FAILED.txt")
)

manifest_e <- read_csv_evidence(evidence_paths$manifest)
automated_e <- read_csv_evidence(evidence_paths$automated)
data_e <- read_csv_evidence(evidence_paths$data_checks)
visual_e <- read_csv_evidence(evidence_paths$visual_qa)
visual_review_e <- read_csv_evidence(evidence_paths$visual_review)
model_warning_e <- read_csv_evidence(evidence_paths$model_warning_replay)
model_warning_modules_e <- read_csv_evidence(evidence_paths$model_warning_modules)
palette_e <- read_csv_evidence(evidence_paths$palette)
separation_e <- read_csv_evidence(evidence_paths$colour_separation)
verification_e <- read_csv_evidence(evidence_paths$final_verification)

gate <- list()
gate_detail <- list()
set_gate <- function(name, pass, detail) {
  gate[[name]] <<- isTRUE(pass)
  gate_detail[[name]] <<- as.character(detail)
}

manifest <- manifest_e$data
manifest_ok <- !is.null(manifest) &&
  all(c("figure_id", "figure_class", "png_path", "pdf_path", "validation_status") %in% names(manifest)) &&
  nrow(manifest) == 18L && !anyDuplicated(manifest$figure_id) &&
  sum(manifest$figure_class == "Main", na.rm = TRUE) == 8L &&
  sum(manifest$figure_class == "Supplementary", na.rm = TRUE) == 10L &&
  all(tolower(trimws(manifest$validation_status)) == "validated") &&
  all(nonempty_files(manifest$png_path)) && all(nonempty_files(manifest$pdf_path))
set_gate(
  "Canonical figure manifest and files", manifest_ok,
  if (is.null(manifest)) manifest_e$error else paste0(
    nrow(manifest), " families; ", sum(manifest$figure_class == "Main", na.rm = TRUE),
    " main; ", sum(manifest$figure_class == "Supplementary", na.rm = TRUE),
    " supplementary; ", sum(nonempty_files(c(manifest$png_path, manifest$pdf_path))),
    "/", 2L * nrow(manifest), " nonempty PNG/PDF files"
  )
)

automated <- automated_e$data
automated_ok <- !is.null(automated) && all(c("pass", "severity") %in% names(automated)) &&
  nrow(automated) > 0L && all(as_flag(automated$pass))
automated_failures <- if (!is.null(automated) && "pass" %in% names(automated)) {
  sum(!as_flag(automated$pass))
} else {
  NA_integer_
}
set_gate(
  "Automated final-figure checks", automated_ok,
  if (is.null(automated)) automated_e$error else paste0(
    sum(as_flag(automated$pass)), "/", nrow(automated), " passed; ",
    automated_failures, " failed"
  )
)

data_checks <- data_e$data
data_ok <- !is.null(data_checks) && all(c("check", "pass") %in% names(data_checks)) &&
  nrow(data_checks) > 0L && all(as_flag(data_checks$pass))
set_gate(
  "Independent plotted-data anchors", data_ok,
  if (is.null(data_checks)) data_e$error else paste0(
    sum(as_flag(data_checks$pass)), "/", nrow(data_checks), " passed"
  )
)

visual <- visual_e$data
visual_required <- c(
  "figure_id", "pdf_render_path", "a4_preview_path", "greyscale_path",
  "deutan_path", "protan_path", "all_derivatives_nonempty"
)
visual_paths <- character()
visual_fresh <- logical()
if (!is.null(visual) && all(visual_required %in% names(visual))) {
  visual_paths <- unlist(visual[c(
    "pdf_render_path", "a4_preview_path", "greyscale_path", "deutan_path", "protan_path"
  )], use.names = FALSE)
  if (!is.null(manifest) && all(c("figure_id", "png_path", "pdf_path") %in% names(manifest))) {
    match_idx <- match(visual$figure_id, manifest$figure_id)
    visual_fresh <- vapply(seq_len(nrow(visual)), function(i) {
      j <- match_idx[[i]]
      if (is.na(j)) return(FALSE)
      source_paths <- c(manifest$png_path[[j]], manifest$pdf_path[[j]])
      derivative_paths <- unlist(visual[i, c(
        "pdf_render_path", "a4_preview_path", "greyscale_path", "deutan_path", "protan_path"
      )], use.names = FALSE)
      if (!all(nonempty_files(c(source_paths, derivative_paths)))) return(FALSE)
      newest_source <- max(file.info(source_paths)$mtime)
      oldest_derivative <- min(file.info(derivative_paths)$mtime)
      as.numeric(difftime(oldest_derivative, newest_source, units = "secs")) >= -2
    }, logical(1))
  }
}
visual_ok <- !is.null(visual) && all(visual_required %in% names(visual)) &&
  nrow(visual) == 18L && !anyDuplicated(visual$figure_id) &&
  !is.null(manifest) && setequal(visual$figure_id, manifest$figure_id) &&
  all(as_flag(visual$all_derivatives_nonempty)) &&
  length(visual_paths) == 90L && all(nonempty_files(visual_paths)) &&
  length(visual_fresh) == 18L && all(visual_fresh)
set_gate(
  "Reproducible visual-QA derivatives", visual_ok,
  if (is.null(visual)) visual_e$error else paste0(
    nrow(visual), " figure families; ", sum(nonempty_files(visual_paths)), "/", length(visual_paths),
    " derivatives nonempty; ", sum(visual_fresh), "/", length(visual_fresh),
    " families current relative to source figures"
  )
)

palette <- palette_e$data
separation <- separation_e$data
palette_recorded <- !is.null(palette) && nrow(palette) > 0L &&
  all(c("semantic_level", "hex", "contrast_on_white") %in% names(palette))
palette_contrast <- if (palette_recorded) {
  suppressWarnings(as.numeric(palette$contrast_on_white))
} else {
  numeric()
}
palette_recorded <- palette_recorded && length(palette_contrast) == nrow(palette) &&
  all(is.finite(palette_contrast))
separation_distance <- if (!is.null(separation) && "ciede2000_distance" %in% names(separation)) {
  suppressWarnings(as.numeric(separation$ciede2000_distance))
} else {
  numeric()
}
separation_ok <- !is.null(separation) && nrow(separation) >= 3L &&
  "ciede2000_distance" %in% names(separation) &&
  length(separation_distance) == nrow(separation) &&
  all(is.finite(separation_distance)) && all(separation_distance >= 10)
set_gate(
  "Semantic-colour accessibility evidence", palette_recorded && separation_ok,
  if (is.null(palette) || is.null(separation)) {
    paste(na.omit(c(palette_e$error, separation_e$error)), collapse = "; ")
  } else if (!palette_recorded || !separation_ok) {
    "palette or colour-separation evidence is malformed or below the declared CIEDE2000 threshold"
  } else {
    paste0(
      nrow(palette), " semantic colours recorded; minimum operational UTI/Not UTI CIEDE2000 distance = ",
      format(min(separation_distance), digits = 4)
    )
  }
)

visual_review <- visual_review_e$data
visual_review_status <- "NOT RECORDED (this is distinct from derivative generation)"
if (!is.null(visual_review)) {
  review_id_col <- intersect(c("figure_id", "id"), names(visual_review))
  review_pass_col <- intersect(c("pass", "review_pass", "visual_qa_pass"), names(visual_review))
  review_status_col <- intersect(c("status", "review_status", "visual_qa_status"), names(visual_review))
  review_pass <- NULL
  if (length(review_pass_col)) {
    review_pass <- as_flag(visual_review[[review_pass_col[[1L]]]])
  } else if (length(review_status_col)) {
    review_pass <- tolower(trimws(visual_review[[review_status_col[[1L]]]])) %in%
      c("pass", "passed", "acceptable", "approved", "complete", "completed")
  } else {
    component_status_cols <- grep("_status$", names(visual_review), value = TRUE)
    if (length(component_status_cols)) {
      review_pass <- apply(visual_review[component_status_cols], 1L, function(x) {
        all(tolower(trimws(as.character(x))) %in%
              c("pass", "passed", "acceptable", "approved", "complete", "completed"))
      })
    }
  }
  review_ids_ok <- length(review_id_col) && !is.null(manifest) &&
    nrow(visual_review) == nrow(manifest) &&
    setequal(visual_review[[review_id_col[[1L]]]], manifest$figure_id)
  review_fresh <- FALSE
  if (file.exists(evidence_paths$visual_review) && length(visual_paths) &&
      all(nonempty_files(visual_paths))) {
    review_fresh <- as.numeric(difftime(
      file.info(evidence_paths$visual_review)$mtime,
      max(file.info(visual_paths)$mtime),
      units = "secs"
    )) >= -2
  }
  if (!is.null(review_pass) && review_ids_ok) {
    review_ok <- length(review_pass) == nrow(visual_review) && all(review_pass) && review_fresh
    visual_review_status <- paste0(
      sum(review_pass), "/", length(review_pass), " figure families recorded as passing",
      if (all(review_pass)) "" else "; REVIEW FAILURES PRESENT",
      "; review table is ", if (review_fresh) "current" else "STALE",
      " relative to the visual-QA derivatives"
    )
    set_gate("Recorded manual visual review", review_ok, visual_review_status)
  } else {
    visual_review_status <- "PRESENT BUT NOT MACHINE-INTERPRETABLE as a complete 18-family review"
    set_gate("Recorded manual visual review", FALSE, visual_review_status)
  }
}

# The replay is an audit-only supplement because it intentionally reruns and
# refreshes model outputs. When present, however, it must be complete,
# classified, and hash-bound to the current three model scripts.
model_warning <- model_warning_e$data
model_warning_modules <- model_warning_modules_e$data
model_warning_status <- "not recorded"
if (file.exists(evidence_paths$model_warning_replay) ||
    file.exists(evidence_paths$model_warning_modules)) {
  required_module_cols <- c("module", "script", "script_sha256", "completed")
  required_warning_cols <- c("condition_type", "classification", "message", "n")
  module_shape_ok <- !is.null(model_warning_modules) &&
    all(required_module_cols %in% names(model_warning_modules)) &&
    nrow(model_warning_modules) == 3L && all(as_flag(model_warning_modules$completed))
  warning_shape_ok <- !is.null(model_warning) &&
    all(required_warning_cols %in% names(model_warning))
  scripts_current <- FALSE
  if (module_shape_ok && requireNamespace("digest", quietly = TRUE)) {
    scripts_current <- all(file.exists(model_warning_modules$script)) && all(vapply(
      seq_len(nrow(model_warning_modules)),
      function(i) identical(
        unname(digest::digest(model_warning_modules$script[[i]], algo = "sha256", file = TRUE)),
        as.character(model_warning_modules$script_sha256[[i]])
      ),
      logical(1)
    ))
  }
  prohibited_classes <- c(
    "prohibited_deprecation", "prohibited_plot_or_join_warning",
    "unclassified_requires_review"
  )
  unresolved_n <- if (warning_shape_ok) {
    sum(model_warning$classification %in% prohibited_classes)
  } else {
    NA_integer_
  }
  replay_ok <- module_shape_ok && warning_shape_ok && scripts_current &&
    isTRUE(unresolved_n == 0L)
  total_conditions <- if (warning_shape_ok) sum(suppressWarnings(as.numeric(model_warning$n)), na.rm = TRUE) else NA_real_
  model_warning_status <- paste0(
    if (module_shape_ok) "3/3 modules completed" else "module evidence incomplete",
    "; current script hashes ", if (scripts_current) "match" else "do not match",
    "; exact condition occurrences ", if (is.finite(total_conditions)) format(total_conditions, trim = TRUE) else "unavailable",
    "; unresolved/prohibited distinct records ", if (!is.na(unresolved_n)) unresolved_n else "unavailable"
  )
  set_gate("Exact model-warning replay", replay_ok, model_warning_status)
}

verification <- verification_e$data
verification_ok <- !is.null(verification) && all(c("stage", "pass") %in% names(verification)) &&
  nrow(verification) > 0L && all(verification$stage == "final") && all(as_flag(verification$pass))
set_gate(
  "Final Longcycler-only release verification", verification_ok,
  if (is.null(verification)) verification_e$error else paste0(
    sum(as_flag(verification$pass)), "/", nrow(verification), " final-stage checks passed"
  )
)

rq_lines <- if (file.exists(evidence_paths$rq_marker)) {
  readLines(evidence_paths$rq_marker, warn = FALSE)
} else {
  character()
}
rq_ok <- length(rq_lines) > 0L && any(grepl("PASS", rq_lines, fixed = TRUE)) &&
  any(grepl("RQ01-+RQ10", rq_lines))
set_gate(
  "RQ01-RQ10 completion marker", rq_ok,
  if (length(rq_lines)) paste(rq_lines, collapse = " | ") else "completion marker missing"
)

complete_lines <- if (file.exists(evidence_paths$complete_marker)) {
  readLines(evidence_paths$complete_marker, warn = FALSE)
} else {
  character()
}
failed_lines <- if (file.exists(evidence_paths$failed_marker)) {
  readLines(evidence_paths$failed_marker, warn = FALSE)
} else {
  character()
}
marker_run_id <- marker_value(complete_lines, "Run ID")
environment_run_id <- trimws(Sys.getenv("RUN_ID", unset = ""))
run_id_ok <- !nzchar(environment_run_id) ||
  (!is.na(marker_run_id) && identical(marker_run_id, environment_run_id))

# Human visual review is deliberately excluded from completion-marker
# freshness: it is performed after computational completion and must instead
# be newer than the QA derivatives (checked above).
freshness_inputs <- unlist(evidence_paths[c(
  "manifest", "automated", "data_checks", "visual_qa", "palette", "colour_separation",
  "final_verification", "rq_marker"
)], use.names = FALSE)
freshness_inputs <- freshness_inputs[file.exists(freshness_inputs)]
marker_fresh <- FALSE
freshness_detail <- "completion marker missing"
if (file.exists(evidence_paths$complete_marker) && length(freshness_inputs)) {
  marker_time <- file.info(evidence_paths$complete_marker)$mtime
  newest_evidence <- max(file.info(freshness_inputs)$mtime)
  marker_fresh <- as.numeric(difftime(marker_time, newest_evidence, units = "secs")) >= -2
  freshness_detail <- paste0(
    "marker modified ", fmt_time(marker_time), "; newest required evidence modified ",
    fmt_time(newest_evidence)
  )
}
complete_marker_ok <- length(complete_lines) > 0L &&
  any(grepl("^Complete analysis: PASS$", complete_lines)) &&
  !length(failed_lines) && run_id_ok && marker_fresh
set_gate(
  "Current pipeline completion marker", complete_marker_ok,
  paste0(
    if (length(complete_lines)) paste(complete_lines, collapse = " | ") else "marker missing",
    "; failed marker ", if (length(failed_lines)) "PRESENT" else "absent",
    "; environment/marker run ID ", if (run_id_ok) "consistent" else "MISMATCH",
    "; ", freshness_detail
  )
)

explicit_log_path <- trimws(Sys.getenv("PIPELINE_LOG_PATH", unset = ""))
log_path <- ""
log_provenance <- "not supplied"
if (nzchar(explicit_log_path)) {
  log_path <- if (startsWith(explicit_log_path, "/")) explicit_log_path else file.path(root, explicit_log_path)
  log_provenance <- "PIPELINE_LOG_PATH"
} else if (!is.na(marker_run_id)) {
  candidate <- file.path(root, "logs", paste0("complete_analysis_", marker_run_id, ".log"))
  if (file.exists(candidate)) {
    log_path <- candidate
    log_provenance <- "inferred from completion-marker Run ID"
  }
}

log_lines <- character()
warning_lines <- character()
warning_classes <- character()
prohibited_warning_lines <- character()
if (nzchar(log_path) && file.exists(log_path)) {
  log_lines <- readLines(log_path, warn = FALSE)
  warning_pattern <- paste0(
    "^[[:space:]]*(\\[[^]]+\\][[:space:]]*){0,2}(⚠️?[[:space:]]*)?",
    "(Warning( message)?s?([[:space:]]+in[^:]*)?:|In addition: Warning|",
    "There (was|were) .*warnings?|\\[(WARN|WARNING)\\]|WARN(ING)?[[:space:]]*:)"
  )
  warning_idx <- grep(warning_pattern, log_lines, ignore.case = TRUE)
  if (length(warning_idx)) {
    warning_lines <- unique(vapply(warning_idx, function(i) {
      header_has_detail_block <- grepl(
        "(Warning( message)?s?:|In addition: Warning)[[:space:]]*$",
        log_lines[[i]], ignore.case = TRUE
      )
      if (!header_has_detail_block) return(trimws(log_lines[[i]]))
      upper <- min(length(log_lines), i + 200L)
      candidate <- log_lines[i:upper]
      if (length(candidate) > 1L) {
        blank <- which(!nzchar(trimws(candidate[-1L])))
        if (length(blank)) candidate <- candidate[seq_len(blank[[1L]])]
      }
      paste(trimws(candidate), collapse = " ")
    }, character(1)))
  }
  classify_warning <- function(value) {
    if (grepl("participants appear in multiple batches|expected if the same individuals", value, ignore.case = TRUE)) {
      "expected longitudinal sampling notice"
    } else if (grepl("manual scale|shared levels|palette|colour scale|color scale", value, ignore.case = TRUE)) {
      "scale/palette"
    } else if (grepl("many-to-many|cardinality|join relationship", value, ignore.case = TRUE)) {
      "join/cardinality"
    } else if (grepl("removed [0-9]+ rows?|missing values?|non[- ]?finite|outside the scale range", value, ignore.case = TRUE)) {
      "removed/missing/nonfinite"
    } else if (grepl("singular|converg|separation|Hessian|boundary fit", value, ignore.case = TRUE)) {
      "model diagnostics"
    } else if (grepl("deprecated|deprecation", value, ignore.case = TRUE)) {
      "deprecated syntax"
    } else if (grepl("built under R version|package", value, ignore.case = TRUE)) {
      "package/runtime"
    } else {
      "other"
    }
  }
  warning_classes <- vapply(warning_lines, classify_warning, character(1))
  prohibited_pattern <- paste(
    "insufficient values in manual scale|no shared levels|unknown .*scale|",
    "removed [0-9]+ rows?|many-to-many|non[- ]?finite.*(effect|estimate)|",
    "clipp(ed|ing).*annotation|unlabeled data points|too many overlaps|failed figure|figure.*failed|stale hard-coded count|",
    "deprecated as of ggplot2|deprecated ggplot|deprecated in tidyselect|",
    "one or more parsing issues",
    sep = ""
  )
  prohibited_warning_lines <- warning_lines[
    grepl(prohibited_pattern, warning_lines, ignore.case = TRUE, perl = TRUE)
  ]
}

if (nzchar(explicit_log_path)) {
  set_gate(
    "Requested pipeline run log is readable",
    nzchar(log_path) && file.exists(log_path),
    if (nzchar(log_path) && file.exists(log_path)) file_summary(log_path) else
      paste0("PIPELINE_LOG_PATH does not resolve to a readable file: ", explicit_log_path)
  )
  if (nzchar(log_path) && file.exists(log_path)) {
    set_gate(
      "No prohibited warning pattern in pipeline log",
      !length(prohibited_warning_lines),
      if (length(prohibited_warning_lines)) {
        paste0(length(prohibited_warning_lines), " prohibited warning record(s) detected")
      } else {
        paste0(length(warning_lines), " warning-like record(s) scanned; none matched prohibited patterns")
      }
    )
  }
}

all_machine_gates_pass <- length(gate) > 0L && all(unlist(gate, use.names = FALSE))
overall_status <- if (length(failed_lines)) {
  "FAIL — a pipeline failure marker is present"
} else if (!all_machine_gates_pass) {
  "INCOMPLETE OR FAIL — one or more required machine gates did not pass"
} else if (length(warning_lines)) {
  "PASS WITH RECORDED WARNING-LIKE LOG LINES — machine gates passed; warnings require review"
} else if (nzchar(log_path) && file.exists(log_path)) {
  "PASS — all required machine gates passed and no warning-like lines were detected in the run log"
} else {
  "PASS — all required machine gates passed; no specific run log was available for warning extraction"
}

failed_checks <- character()
if (!is.null(automated) && all(c("pass", "scope", "check", "observed", "expected") %in% names(automated))) {
  bad <- automated[!as_flag(automated$pass), , drop = FALSE]
  if (nrow(bad)) {
    failed_checks <- apply(bad, 1L, function(x) paste0(
      x[["scope"]], " — ", x[["check"]], ": observed ", x[["observed"]],
      "; expected ", x[["expected"]]
    ))
  }
}

data_lines <- character()
if (!is.null(data_checks) && all(c("check", "observed", "expected", "pass") %in% names(data_checks))) {
  data_lines <- apply(data_checks, 1L, function(x) paste0(
    "- ", x[["check"]], ": observed ", x[["observed"]],
    "; expected ", x[["expected"]], "; ",
    if (as_flag(x[["pass"]])) "PASS" else "FAIL"
  ))
}

gate_lines <- vapply(names(gate), function(name) paste0(
  "- ", name, ": ", if (gate[[name]]) "PASS" else "FAIL", " — ", gate_detail[[name]]
), character(1))

source_lines <- vapply(evidence_paths, file_summary, character(1))
session_lines <- capture.output(utils::sessionInfo())

report <- c(
  "FIGURE AUDIT VALIDATION RESULTS",
  "===============================",
  paste0("Generated: ", fmt_time(Sys.time())),
  paste0("Repository: ", root),
  paste0("Run ID from environment: ", if (nzchar(environment_run_id)) environment_run_id else "not supplied"),
  paste0("Run ID from completion marker: ", if (!is.na(marker_run_id)) marker_run_id else "unavailable"),
  "",
  paste0("OVERALL STATUS: ", overall_status),
  "",
  "Machine gates",
  "-------------",
  gate_lines,
  "",
  "Final figure release",
  "--------------------",
  if (is.null(manifest)) {
    paste0("- Manifest unavailable: ", manifest_e$error)
  } else {
    c(
      paste0("- Figure families: ", nrow(manifest), " (", sum(manifest$figure_class == "Main"),
             " main; ", sum(manifest$figure_class == "Supplementary"), " supplementary)"),
      paste0("- PNG/PDF files: ", sum(nonempty_files(c(manifest$png_path, manifest$pdf_path))),
             "/", 2L * nrow(manifest), " exist and are nonempty"),
      paste0("- Manifest validation labels: ", paste(sort(unique(manifest$validation_status)), collapse = ", "))
    )
  },
  paste0("- Automated checks: ", if (is.null(automated)) "unavailable" else paste0(
    sum(as_flag(automated$pass)), "/", nrow(automated), " passed"
  )),
  if (length(failed_checks)) c("- Failed automated checks:", paste0("  * ", failed_checks)) else
    "- Failed automated checks: none",
  "",
  "Independently recomputed displayed-data anchors",
  "-----------------------------------------------",
  if (length(data_lines)) data_lines else paste0("- Unavailable: ", data_e$error),
  "",
  "Exact model-warning replay",
  "--------------------------",
  paste0("- Status: ", model_warning_status),
  if (!is.null(model_warning) && all(c("condition_type", "classification", "n") %in% names(model_warning))) {
    aggregate(
      suppressWarnings(as.numeric(model_warning$n)),
      by = list(type = model_warning$condition_type, class = model_warning$classification),
      FUN = sum,
      na.rm = TRUE
    ) |>
      transform(line = paste0("- ", type, " / ", class, ": ", x, " occurrence(s)")) |>
      (`[[`)("line")
  } else {
    "- No machine-readable replay table was present."
  },
  "",
  "Visual and accessibility QA",
  "---------------------------",
  if (is.null(visual)) {
    paste0("- Derivative-generation table unavailable: ", visual_e$error)
  } else {
    c(
      paste0("- Figure families represented: ", nrow(visual)),
      paste0("- Nonempty PDF/A4/greyscale/deutan/protan derivatives: ",
             sum(nonempty_files(visual_paths)), "/", length(visual_paths)),
      paste0("- Families whose derivatives are current relative to their source figure: ",
             sum(visual_fresh), "/", length(visual_fresh))
    )
  },
  paste0("- Manual full-size/thesis-size review: ", visual_review_status),
  if (is.null(palette)) paste0("- Palette table unavailable: ", palette_e$error) else paste0(
    "- Palette entries recorded: ", nrow(palette),
    "; contrast ratios on white range ",
    if (length(palette_contrast) && any(is.finite(palette_contrast))) {
      paste0(format(min(palette_contrast, na.rm = TRUE), digits = 4), "–",
             format(max(palette_contrast, na.rm = TRUE), digits = 4))
    } else {
      "unavailable"
    }
  ),
  if (is.null(separation)) paste0("- Colour-separation table unavailable: ", separation_e$error) else paste0(
    "- Operational UTI/Not UTI CIEDE2000 distances: ",
    if (length(separation_distance) == nrow(separation) &&
        "simulation" %in% names(separation) && all(is.finite(separation_distance))) {
      paste(paste0(separation$simulation, "=", format(separation_distance, digits = 4)), collapse = "; ")
    } else {
      "unavailable because the table is malformed"
    }
  ),
  "- Contrast ratios are recorded as evidence, not represented as WCAG text-conformance claims; redundant shape/label encodings remain necessary.",
  "",
  "Pipeline markers and warning review",
  "-----------------------------------",
  paste0("- Completion marker: ", if (length(complete_lines)) paste(complete_lines, collapse = " | ") else "MISSING"),
  paste0("- Failure marker: ", if (length(failed_lines)) paste(failed_lines, collapse = " | ") else "absent"),
  paste0("- Run log: ", if (nzchar(log_path)) paste0(rel(log_path), " (", log_provenance, ")") else "not available"),
  if (nzchar(log_path) && !file.exists(log_path)) "- Run-log warning extraction: requested log is missing" else
    if (length(log_lines)) paste0("- Run-log lines scanned: ", length(log_lines)) else
      "- Run-log warning extraction: not performed because no readable run log was available",
  paste0("- Warning-like lines detected by anchored log scan: ", length(warning_lines)),
  if (length(warning_lines)) c(
    "- Warning-like records (deduplicated and heuristically classified; these require human review):",
    paste0("  * [", head(warning_classes, 200L), "] ", head(warning_lines, 200L)),
    if (length(warning_lines) > 200L) paste0("  * ... ", length(warning_lines) - 200L, " additional unique lines omitted")
  ) else NULL,
  paste0("- Prohibited warning-pattern records: ", length(prohibited_warning_lines)),
  "",
  "Evidence files",
  "--------------",
  paste0("- ", source_lines),
  "",
  "R session information",
  "---------------------",
  session_lines,
  "",
  "Interpretation",
  "--------------",
  "A PASS above applies only to the machine gates explicitly listed in this report. Derivative generation is not a substitute for human visual inspection; when a manual review table is absent, that limitation is stated rather than inferred away. Warning-like log extraction is intentionally conservative and requires review of the cited run log."
)

tmp <- tempfile(pattern = "validation_results_", tmpdir = audit_dir, fileext = ".txt")
writeLines(report, tmp, useBytes = TRUE)
if (!file.rename(tmp, output_path)) {
  unlink(tmp)
  stop("Could not atomically publish ", rel(output_path), call. = FALSE)
}

message("Wrote ", rel(output_path), ": ", overall_status)
if (!all_machine_gates_pass || length(failed_lines)) quit(save = "no", status = 1L)

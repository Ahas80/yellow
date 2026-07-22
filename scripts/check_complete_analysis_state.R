#!/usr/bin/env Rscript

# Evidence-first monitor for the scheduled Longcycler-only continuation.
# It never declares completion itself; it records process, marker, log-growth,
# and explicitly scoped file-change evidence for the supervising heartbeat.

root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
if (!file.exists(file.path(root, "RUN_COMPLETE_ANALYSIS.sh"))) {
  stop("Run this monitor from the rUTIs project root.")
}

pipeline_dir <- file.path(root, "results", "pipeline")
dir.create(pipeline_dir, recursive = TRUE, showWarnings = FALSE)
state_path <- file.path(pipeline_dir, "heartbeat_monitor_state.rds")
snapshot_path <- file.path(pipeline_dir, "heartbeat_monitor_files.csv")
complete_path <- file.path(pipeline_dir, "RUN_COMPLETE.txt")
failed_path <- file.path(pipeline_dir, "RUN_FAILED.txt")

previous <- if (file.exists(state_path)) {
  tryCatch(readRDS(state_path), error = function(e) NULL)
} else {
  NULL
}

logs <- list.files(
  file.path(root, "logs"),
  pattern = "^complete_analysis_[0-9T+]+\\.log$",
  full.names = TRUE
)
if (length(logs) > 0L) {
  log_info <- file.info(logs)
  latest_log <- logs[[which.max(log_info$mtime)]]
  latest_log_size <- unname(file.size(latest_log))
  latest_log_mtime <- unname(file.info(latest_log)$mtime)
} else {
  latest_log <- NA_character_
  latest_log_size <- NA_real_
  latest_log_mtime <- as.POSIXct(NA)
}

previous_log_size <- if (!is.null(previous) && identical(previous$latest_log, latest_log)) {
  previous$latest_log_size
} else {
  NA_real_
}
log_delta_bytes <- if (is.na(latest_log_size) || is.na(previous_log_size)) {
  NA_real_
} else {
  latest_log_size - previous_log_size
}

ps_output <- tryCatch(
  system2(
    "ps",
    c("-Ao", "pid=,ppid=,etime=,%cpu=,%mem=,command="),
    stdout = TRUE,
    stderr = TRUE
  ),
  error = function(e) structure(conditionMessage(e), status = 1L)
)
ps_status <- attr(ps_output, "status")
if (is.null(ps_status)) ps_status <- 0L
process_check_ok <- identical(as.integer(ps_status), 0L)
process_pattern <- paste(
  c(
    "RUN_IN_TERMINAL\\.sh", "RUN_COMPLETE_ANALYSIS\\.sh", "Rscript",
    "panaroo", "parsnp", "abricate", "mlst"
  ),
  collapse = "|"
)
active_processes <- if (process_check_ok) {
  ps_output[
    grepl(process_pattern, ps_output, ignore.case = TRUE) &
      !grepl("check_complete_analysis_state\\.R", ps_output, fixed = FALSE) &
      !grepl("ps -Ao", ps_output, fixed = TRUE)
  ]
} else {
  character()
}
live_process <- length(active_processes) > 0L

scan_roots <- file.path(
  root,
  c(
    "results/clinical", "results/qc", "results/vf", "results/mlst",
    "results/mlst_source_comparison", "results/strain_compare",
    "results/models", "results/longitudinal", "results/mechanism",
    "results/robustness", "results/research_questions", "results/pipeline",
    "results/summary", "results/audit", "plots/clinical", "plots/vf",
    "plots/mlst", "plots/plasmids", "plots/publication", "plots/wgs"
  )
)
scan_roots <- scan_roots[dir.exists(scan_roots)]
active_files <- unique(unlist(lapply(scan_roots, function(path) {
  list.files(path, recursive = TRUE, full.names = TRUE, all.files = TRUE, no.. = TRUE)
})))
active_files <- active_files[file.exists(active_files) & !dir.exists(active_files)]
active_files <- setdiff(normalizePath(active_files, winslash = "/", mustWork = FALSE), c(
  normalizePath(state_path, winslash = "/", mustWork = FALSE),
  normalizePath(snapshot_path, winslash = "/", mustWork = FALSE)
))

if (length(active_files) > 0L) {
  active_info <- file.info(active_files)
  current_files <- data.frame(
    path = active_files,
    size_bytes = as.numeric(active_info$size),
    mtime_epoch = as.numeric(active_info$mtime),
    stringsAsFactors = FALSE
  )
} else {
  current_files <- data.frame(
    path = character(), size_bytes = numeric(), mtime_epoch = numeric(),
    stringsAsFactors = FALSE
  )
}

prior_files <- if (file.exists(snapshot_path)) {
  tryCatch(read.csv(snapshot_path, stringsAsFactors = FALSE), error = function(e) NULL)
} else {
  NULL
}
if (is.null(prior_files) || !all(c("path", "size_bytes", "mtime_epoch") %in% names(prior_files))) {
  changed <- rep(TRUE, nrow(current_files))
  deleted_count <- 0L
  snapshot_baseline <- TRUE
} else {
  prior_index <- match(current_files$path, prior_files$path)
  changed <- is.na(prior_index) |
    current_files$size_bytes != prior_files$size_bytes[prior_index] |
    abs(current_files$mtime_epoch - prior_files$mtime_epoch[prior_index]) > 0.001
  deleted_count <- sum(!prior_files$path %in% current_files$path)
  snapshot_baseline <- FALSE
}
changed[is.na(changed)] <- TRUE
changed_file_count <- sum(changed)
changed_file_bytes <- sum(current_files$size_bytes[changed], na.rm = TRUE)

read_log_tail <- function(path, max_bytes = 65536L, max_lines = 20L) {
  if (is.na(path) || !file.exists(path)) return(character())
  size <- file.size(path)
  con <- file(path, open = "rb")
  on.exit(close(con), add = TRUE)
  seek(con, where = max(0, size - max_bytes), origin = "start")
  raw <- readBin(con, what = "raw", n = min(max_bytes, size))
  text <- iconv(rawToChar(raw), from = "UTF-8", to = "UTF-8", sub = "")
  if (is.na(text)) return(character())
  lines <- unlist(strsplit(text, "[\\r\\n]+", perl = TRUE), use.names = FALSE)
  lines <- trimws(lines)
  lines <- lines[nzchar(lines) & !grepl("^Progress:", lines)]
  tail(lines, max_lines)
}
meaningful_tail <- read_log_tail(latest_log)

complete_present <- file.exists(complete_path) && file.size(complete_path) > 0L
failed_present <- file.exists(failed_path) && file.size(failed_path) > 0L
log_advanced <- !is.na(log_delta_bytes) && log_delta_bytes > 0

status <- if (failed_present) {
  if (live_process) "FAILED_RECOVERY_PROCESS_PRESENT" else "FAILED_STOPPED"
} else if (complete_present) {
  "COMPLETE_MARKER_PRESENT"
} else if (live_process && (log_advanced || length(active_processes) > 0L)) {
  "RUNNING_CONFIRMED"
} else if (!process_check_ok) {
  "UNKNOWN_PROCESS_CHECK_FAILED"
} else {
  "STOPPED_WITHOUT_COMPLETION"
}

cat(sprintf("monitor_time=%s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S %z")))
cat(sprintf("status=%s\n", status))
cat(sprintf("complete_marker=%s\n", complete_present))
cat(sprintf("failed_marker=%s\n", failed_present))
cat(sprintf("process_check_ok=%s\n", process_check_ok))
cat(sprintf("live_process_count=%d\n", length(active_processes)))
cat(sprintf("latest_log=%s\n", ifelse(is.na(latest_log), "<none>", latest_log)))
cat(sprintf("latest_log_size_bytes=%s\n", ifelse(is.na(latest_log_size), "NA", latest_log_size)))
cat(sprintf("latest_log_mtime=%s\n", ifelse(is.na(latest_log_mtime), "NA", format(latest_log_mtime, "%Y-%m-%d %H:%M:%S %z"))))
cat(sprintf("log_delta_bytes_since_previous_check=%s\n", ifelse(is.na(log_delta_bytes), "NA", log_delta_bytes)))
cat(sprintf("file_snapshot_is_baseline=%s\n", snapshot_baseline))
cat(sprintf("changed_files_in_active_analysis_scope=%d\n", changed_file_count))
cat(sprintf("changed_file_bytes_in_active_analysis_scope=%.0f\n", changed_file_bytes))
cat(sprintf("deleted_files_in_active_analysis_scope=%d\n", deleted_count))
if (length(active_processes) > 0L) {
  cat("active_processes:\n")
  cat(paste0("  ", active_processes, collapse = "\n"), "\n", sep = "")
}
if (length(meaningful_tail) > 0L) {
  cat("latest_meaningful_log_lines:\n")
  cat(paste0("  ", meaningful_tail, collapse = "\n"), "\n", sep = "")
}
if (failed_present) {
  cat("failure_marker_contents:\n")
  cat(paste0("  ", readLines(failed_path, warn = FALSE), collapse = "\n"), "\n", sep = "")
}

write.csv(current_files, snapshot_path, row.names = FALSE, na = "")
saveRDS(
  list(
    checked_at = Sys.time(), latest_log = latest_log,
    latest_log_size = latest_log_size, status = status
  ),
  state_path
)

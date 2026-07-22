#!/usr/bin/env Rscript

# Restore the RQ01-RQ10 completion marker only from a fully validated persisted
# release. This avoids rerunning the expensive research-question analyses when
# a downstream cleanup accidentally removes the marker but leaves the complete
# audited outputs intact.

root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
rq_root <- file.path(root, "results", "research_questions")
expected_rq <- sprintf("RQ%02d", 1:10)

stop_release <- function(message) stop(message, call. = FALSE)
read_required_csv <- function(path) {
  if (!file.exists(path)) stop_release(paste("Missing required RQ release file:", path))
  read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
}

status_paths <- file.path(rq_root, expected_rq, "analysis_status.csv")
if (any(!file.exists(status_paths))) {
  stop_release("One or more RQ01-RQ10 analysis_status.csv files are missing.")
}
status_rows <- do.call(rbind, lapply(status_paths, read_required_csv))
if (nrow(status_rows) != 10L ||
    !identical(as.character(status_rows$research_question), expected_rq) ||
    any(as.character(status_rows$status) != "complete")) {
  stop_release("RQ01-RQ10 status files do not form one complete ordered release.")
}

final_status <- read_required_csv(file.path(rq_root, "final_question_status.csv"))
if (nrow(final_status) != 10L ||
    !identical(as.character(final_status$research_question), expected_rq) ||
    any(as.character(final_status$status) != "complete") ||
    any(!file.exists(as.character(final_status$status_file)))) {
  stop_release("final_question_status.csv does not validate all RQ01-RQ10 outputs.")
}

contracts <- read_required_csv(file.path(rq_root, "final_contract_checks.csv"))
contract_pass <- tolower(trimws(as.character(contracts$pass))) %in% c("true", "t", "1")
if (!nrow(contracts) || any(!contract_pass)) {
  stop_release("One or more final RQ denominator contracts failed.")
}

execution <- read_required_csv(file.path(rq_root, "_inputs", "execution_status.csv"))
if (nrow(execution) != 1L ||
    as.character(execution$state[[1]]) != "COMPLETE" ||
    as.integer(execution$boot_reps[[1]]) != 10000L) {
  stop_release("RQ06-RQ08 execution status is not a complete 10,000-bootstrap release.")
}

retired <- list.files(
  rq_root,
  pattern = "rq11|rq09_11",
  recursive = TRUE,
  full.names = TRUE,
  ignore.case = TRUE,
  all.files = TRUE,
  include.dirs = TRUE,
  no.. = TRUE
)
if (length(retired)) stop_release("Retired RQ11 artifacts remain in the active release.")

marker <- file.path(rq_root, "RUN_COMPLETE.txt")
tmp <- tempfile(pattern = ".RUN_COMPLETE.", tmpdir = rq_root, fileext = ".txt")
writeLines(c(
  "Research-question analysis runner: PASS",
  sprintf("Restored after persisted-release validation: %s", format(Sys.time(), tz = "Europe/Amsterdam", usetz = TRUE)),
  "All RQ01-RQ10 modules are complete.",
  "All final denominator contracts pass.",
  "RQ06-RQ08 execution records 10,000 bootstrap replicates.",
  "No retired RQ11 artifact is present."
), tmp, useBytes = TRUE)
if (!file.rename(tmp, marker)) stop_release("Could not atomically restore the RQ completion marker.")

message("Validated persisted RQ01-RQ10 release and restored: ", marker)

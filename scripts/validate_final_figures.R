#!/usr/bin/env Rscript

# Programmatic acceptance checks for the canonical thesis figure pack.

source("00_config.R")

suppressPackageStartupMessages({
  library(digest)
  library(dplyr)
  library(png)
  library(readr)
  library(stringr)
  library(tibble)
  library(tidyr)
})

options(warn = 2)

audit_dir <- file.path(DIR_RESULTS, "figure_audit")
manifest_path <- file.path(audit_dir, "final_figure_manifest.csv")
save_manifest_path <- file.path(audit_dir, "save_ruti_figure_manifest.csv")
if (!file.exists(manifest_path) || !file.exists(save_manifest_path)) {
  stop("Final figure manifests are missing; run 35_final_figure_pack.R first.", call. = FALSE)
}

manifest <- read_csv(manifest_path, show_col_types = FALSE)
save_manifest <- read_csv(save_manifest_path, show_col_types = FALSE)

expected_ids <- c(
  sprintf("Fig%02d_%s", 1:8, c(
    "cohort_and_denominators", "wgs_quality_control", "sequence_type_distribution",
    "vf_burden", "vf_association_evidence", "longitudinal_trajectories",
    "within_host_genomic_continuity", "reference_aware_variant_map"
  )),
  sprintf("FigS%02d_%s", 1:10, c(
    "vf_presence_heatmap", "core_genome_phylogeny", "module_gain_loss", "vf_pcoa",
    "near_miss_leave_one_uti", "transition_mechanisms", "plasmid_amr_context",
    "gene_model_forest", "event_sample_sensitivity", "snp_threshold_sensitivity"
  ))
)
expected_plot_titles <- setNames(c(
  "Cohort selection and analytical denominators",
  "Whole-genome assembly quality control",
  "Sequence-type distribution and provenance",
  "Virulence-factor burden and prespecified score distributions",
  "Prespecified score-level associations with operational UTI",
  "Deidentified trajectories for participants with a Not UTI-to-UTI transition",
  "Within-host genomic continuity and virulence-factor similarity",
  "Reference-aware within-host variant map",
  "Curated virulence-factor presence heatmap",
  "Canonical core-genome neighbour-joining tree",
  "VF module changes across Not UTI-to-UTI transitions",
  "Global VF-profile principal coordinates",
  "Operational-definition and sparse-case diagnostics",
  "Evidence categories across adjacent operational-status transitions",
  "Replicon detection, predicted plasmid context and localized AMR/VF changes",
  "Prescreened gene-level model evidence",
  "Nearest within-participant event-sample sensitivity",
  "Sensitivity to the operational SNP reference"
), expected_ids)

checks <- list()
add_check <- function(scope, check, observed, expected, pass, severity = "critical", detail = "") {
  checks[[length(checks) + 1L]] <<- tibble(
    scope = as.character(scope), check = as.character(check), observed = as.character(observed),
    expected = as.character(expected), pass = isTRUE(pass), severity = severity, detail = as.character(detail)
  )
}

add_check("manifest", "figure family count", nrow(manifest), 18L, nrow(manifest) == 18L)
add_check("manifest", "expected figure IDs", length(intersect(manifest$figure_id, expected_ids)), length(expected_ids),
          setequal(manifest$figure_id, expected_ids))
add_check("manifest", "unique figure IDs", anyDuplicated(manifest$figure_id), 0L, !anyDuplicated(manifest$figure_id))
add_check("manifest", "eight main figures", sum(manifest$figure_class == "Main"), 8L,
          sum(manifest$figure_class == "Main") == 8L)
add_check("manifest", "ten supplementary figures", sum(manifest$figure_class == "Supplementary"), 10L,
          sum(manifest$figure_class == "Supplementary") == 10L)

required_meta <- c(
  "figure_id", "figure_number", "title", "scientific_question", "figure_class", "png_path", "pdf_path",
  "width_in", "height_in", "dpi", "caption", "source_inputs", "statistical_method", "caveat", "unit",
  "filters", "visual_encodings", "multiplicity", "validation_status"
)
add_check("manifest", "required metadata columns", length(intersect(required_meta, names(manifest))), length(required_meta),
          all(required_meta %in% names(manifest)))
for (column in intersect(required_meta, names(manifest))) {
  value <- manifest[[column]]
  if (is.character(value)) {
    ok <- all(!is.na(value) & nzchar(trimws(value)))
    add_check("manifest", paste0("nonblank ", column), sum(!is.na(value) & nzchar(trimws(value))), nrow(manifest), ok)
  }
}

pdfinfo_bin <- unname(Sys.which("pdfinfo"))
if (!nzchar(pdfinfo_bin)) stop("pdfinfo is required for PDF page validation.", call. = FALSE)
pdftotext_bin <- unname(Sys.which("pdftotext"))
if (!nzchar(pdftotext_bin)) {
  candidate <- file.path(
    path.expand("~/.cache/codex-runtimes/codex-primary-runtime/dependencies/native/poppler/poppler/bin"),
    "pdftotext"
  )
  if (file.exists(candidate)) pdftotext_bin <- candidate
}
if (!nzchar(pdftotext_bin) || !file.exists(pdftotext_bin)) {
  stop("pdftotext is required for reader-facing label and deidentification checks.", call. = FALSE)
}

extract_pdf_text <- function(path) {
  out <- tempfile(fileext = ".txt")
  on.exit(unlink(out), add = TRUE)
  status <- system2(pdftotext_bin, c("-layout", shQuote(path), shQuote(out)), stdout = TRUE, stderr = TRUE)
  exit_status <- attr(status, "status")
  if (!is.null(exit_status) && exit_status != 0L) stop("pdftotext failed for ", path, call. = FALSE)
  paste(readLines(out, warn = FALSE), collapse = "\n")
}

parse_pdfinfo <- function(path) {
  out <- system2(pdfinfo_bin, shQuote(path), stdout = TRUE, stderr = TRUE)
  exit_status <- attr(out, "status")
  if (!is.null(exit_status) && exit_status != 0L) stop("pdfinfo failed for ", path, call. = FALSE)
  pages <- as.integer(str_match(out[str_detect(out, "^Pages:")][1], "Pages:\\s+([0-9]+)")[, 2])
  size_line <- out[str_detect(out, "^Page size:")][1]
  size <- as.numeric(str_match(size_line, "Page size:\\s+([0-9.]+) x ([0-9.]+) pts")[, 2:3])
  list(pages = pages, width_pt = size[[1]], height_pt = size[[2]], raw = out)
}

cohort <- read_csv(FILE_ANALYSIS_CLINICAL_COHORT, show_col_types = FALSE)
participant_tokens <- unique(as.character(cohort$Participant_id))
isolate_tokens <- unique(as.character(cohort$Isolate_ID))
raw_tokens <- c(participant_tokens, isolate_tokens)
raw_tokens <- raw_tokens[!is.na(raw_tokens) & nzchar(raw_tokens)]

internal_patterns <- c(
  "\\bNot_UTI\\b", "\\bUTI_Status\\b", "\\bInfection_Status\\b", "\\btp_lab\\b",
  "\\bvf_count_total\\b", "\\bcat_Iron_acquisition\\b", "\\bneglog10FDR\\b", "\\bPos_Ref\\b"
)

for (i in seq_len(nrow(manifest))) {
  id <- manifest$figure_id[[i]]
  png_path <- manifest$png_path[[i]]
  pdf_path <- manifest$pdf_path[[i]]
  expected_width_px <- round(manifest$width_in[[i]] * manifest$dpi[[i]])
  expected_height_px <- round(manifest$height_in[[i]] * manifest$dpi[[i]])

  for (path in c(png_path, pdf_path)) {
    exists <- file.exists(path)
    add_check(id, paste0(tools::file_ext(path), " exists"), exists, TRUE, exists)
    bytes <- if (exists) file.info(path)$size else NA_real_
    add_check(id, paste0(tools::file_ext(path), " nonempty"), bytes, ">0", exists && is.finite(bytes) && bytes > 0)
  }
  if (!file.exists(png_path) || !file.exists(pdf_path)) next

  png_obj <- readPNG(png_path, info = TRUE)
  png_info <- attr(png_obj, "info")
  png_dim <- png_info$dim
  add_check(id, "PNG width", png_dim[[1]], expected_width_px, png_dim[[1]] == expected_width_px)
  add_check(id, "PNG height", png_dim[[2]], expected_height_px, png_dim[[2]] == expected_height_px)
  add_check(id, "PNG horizontal DPI", round(png_info$dpi[[1]], 3), manifest$dpi[[i]],
            abs(png_info$dpi[[1]] - manifest$dpi[[i]]) < .01)
  add_check(id, "PNG vertical DPI", round(png_info$dpi[[2]], 3), manifest$dpi[[i]],
            abs(png_info$dpi[[2]] - manifest$dpi[[i]]) < .01)

  pdf_sig <- readBin(pdf_path, what = "raw", n = 5L)
  add_check(id, "PDF signature", rawToChar(pdf_sig), "%PDF-", identical(rawToChar(pdf_sig), "%PDF-"))
  pdf_info <- parse_pdfinfo(pdf_path)
  add_check(id, "PDF page count", pdf_info$pages, 1L, identical(pdf_info$pages, 1L))
  add_check(id, "PDF width", round(pdf_info$width_pt, 2), round(manifest$width_in[[i]] * 72, 2),
            abs(pdf_info$width_pt - manifest$width_in[[i]] * 72) <= 1)
  add_check(id, "PDF height", round(pdf_info$height_pt, 2), round(manifest$height_in[[i]] * 72, 2),
            abs(pdf_info$height_pt - manifest$height_in[[i]] * 72) <= 1)

  save_row <- save_manifest %>% filter(.data$figure_id == id)
  add_check(id, "one save-manifest row", nrow(save_row), 1L, nrow(save_row) == 1L)
  if (nrow(save_row) == 1L) {
    png_hash <- digest(png_path, algo = "sha256", file = TRUE)
    pdf_hash <- digest(pdf_path, algo = "sha256", file = TRUE)
    add_check(id, "PNG SHA-256", png_hash, save_row$png_sha256[[1]], identical(png_hash, save_row$png_sha256[[1]]))
    add_check(id, "PDF SHA-256", pdf_hash, save_row$pdf_sha256[[1]], identical(pdf_hash, save_row$pdf_sha256[[1]]))
  }

  pdf_text <- extract_pdf_text(pdf_path)
  pdf_text_normalized <- pdf_text %>%
    str_replace_all("[\u2212\u2010\u2011\u2012\u2013\u2014]", "-") %>%
    str_squish()
  internal_hits <- internal_patterns[vapply(
    internal_patterns,
    function(pattern) str_detect(pdf_text, regex(pattern)),
    logical(1)
  )]
  add_check(id, "no internal reader labels", length(internal_hits), 0L, length(internal_hits) == 0L,
            detail = paste(internal_hits, collapse = "; "))
  na_hit <- str_detect(pdf_text, regex("(^|[^[:alnum:]_])NA([^[:alnum:]_]|$)"))
  add_check(id, "no accidental NA legend text", na_hit, FALSE, !na_hit)
  raw_hits <- raw_tokens[vapply(raw_tokens, function(token) str_detect(pdf_text, fixed(token)), logical(1))]
  add_check(id, "no raw participant/isolate identifiers", length(raw_hits), 0L, length(raw_hits) == 0L,
            detail = paste(head(raw_hits, 10), collapse = "; "))
  title_present <- str_detect(pdf_text_normalized, fixed(expected_plot_titles[[id]]))
  add_check(id, "reader-facing title present in PDF", title_present, TRUE, title_present)
}

pairwise <- read_csv(file.path(DIR_RESULTS, "strain_compare", "pairwise_metrics.csv"), show_col_types = FALSE)
transitions <- read_csv(file.path(DIR_RESULTS, "longitudinal", "longcycler_transitions.csv"), show_col_types = FALSE)
anchors <- tibble(
  name = c("episodes", "participants", "UTI", "Not UTI", "direct pairs", "adjacent pairs",
           "Not UTI-to-UTI", "Not UTI-to-UTI <=25 SNPs"),
  observed = c(nrow(cohort), n_distinct(cohort$Participant_id), sum(cohort$UTI_Status == "UTI"),
               sum(cohort$UTI_Status == "Not_UTI"), nrow(pairwise), nrow(transitions),
               sum(transitions$status_from == "Not_UTI" & transitions$status_to == "UTI"),
               sum(transitions$status_from == "Not_UTI" & transitions$status_to == "UTI" & transitions$TotalSNPs <= 25)),
  expected = c(532, 161, 16, 516, 893, 371, 9, 5)
)
for (i in seq_len(nrow(anchors))) {
  add_check("source anchors", anchors$name[[i]], anchors$observed[[i]], anchors$expected[[i]],
            anchors$observed[[i]] == anchors$expected[[i]])
}

caption_required <- c("unit", "filter", "sample", "model", "repeat", "adjust")
caption_blob <- str_to_lower(paste(manifest$caption, manifest$statistical_method, manifest$caveat,
                                   manifest$unit, manifest$filters, manifest$multiplicity))
caption_terms <- list(
  unit = "episode|participant|pair|assembly|variant|gene|module|tip",
  filter = "selected|filter|prevalence|prescreen|candidate|transition|criteri|threshold|nearest|available",
  sample = "[0-9]+",
  model = "model|descriptive|bootstrap|pcoa|neighbour|quality|presence|denominator",
  repeated = "repeat|nested|participant|within-host",
  adjust = "benjamini|holm|none|no inferential|no p-value|not applied"
)
names(caption_terms)[names(caption_terms) == "repeated"] <- "repeat"
for (term in caption_required) {
  present <- str_detect(caption_blob, regex(caption_terms[[term]], ignore_case = TRUE))
  add_check("captions", paste0("all captions contain ", term, " context"), sum(present), nrow(manifest), all(present),
            severity = if (term %in% c("unit", "repeat")) "critical" else "major")
}

result <- bind_rows(checks)
dir.create(audit_dir, recursive = TRUE, showWarnings = FALSE)
write_csv(result, file.path(audit_dir, "automated_validation_checks.csv"))

failed <- result %>% filter(!.data$pass)
summary_lines <- c(
  "FINAL FIGURE AUTOMATED VALIDATION",
  "=================================",
  paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  paste0("Figure families: ", nrow(manifest), " (", sum(manifest$figure_class == "Main"), " main; ",
         sum(manifest$figure_class == "Supplementary"), " supplementary)"),
  paste0("Files checked: ", 2L * nrow(manifest)),
  paste0("Checks passed: ", sum(result$pass), "/", nrow(result)),
  paste0("Critical failures: ", sum(!result$pass & result$severity == "critical")),
  paste0("Other failures: ", sum(!result$pass & result$severity != "critical")),
  "",
  "Independently recomputed anchors:",
  paste0("- ", anchors$name, ": ", anchors$observed, " (expected ", anchors$expected, ")"),
  "",
  if (nrow(failed)) c("FAILED CHECKS:", capture.output(print(failed, n = Inf, width = Inf))) else "All automated checks passed.",
  "",
  "Visual inspection and full-pipeline status are appended to validation_results.txt by the repository audit workflow."
)
writeLines(summary_lines, file.path(audit_dir, "automated_validation_results.txt"))

if (nrow(failed)) {
  stop("Final figure automated validation failed; see results/figure_audit/automated_validation_checks.csv.", call. = FALSE)
}
message("Final figure automated validation passed: ", nrow(manifest), " families; ", nrow(result), " checks.")

#!/usr/bin/env Rscript

# Final, post-regeneration gate for the active Longcycler-only communication
# release. The analytical pipeline has its own gate; this script verifies that
# the canonical decks, handouts, codebook and current release documentation
# were rebuilt from the audited claim registry and contain no retired claims.

suppressPackageStartupMessages({
  library(digest)
  library(jsonlite)
  library(readr)
  library(stringr)
  library(tibble)
})

root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
if (!file.exists(file.path(root, "00_config.R"))) {
  stop("Run from the rUTIs project root.", call. = FALSE)
}

registry_path <- file.path(root, "results", "pipeline", "longcycler_release_claim_registry.json")
report_csv <- file.path(root, "results", "qc", "longcycler_release_deliverables_verification.csv")
report_txt <- file.path(root, "results", "qc", "longcycler_release_deliverables_verification.txt")
dir.create(dirname(report_csv), recursive = TRUE, showWarnings = FALSE)

bundled_pdftotext <- "/Users/Aamir/.cache/codex-runtimes/codex-primary-runtime/dependencies/native/poppler/poppler/bin/pdftotext"
pdftotext_bin <- Sys.getenv("PDFTOTEXT_BIN", unset = bundled_pdftotext)
if (!nzchar(pdftotext_bin) || !file.exists(pdftotext_bin) || file.access(pdftotext_bin, mode = 1L) != 0L) {
  pdftotext_bin <- unname(Sys.which("pdftotext"))
}

v3 <- "outputs/manual-20260527-current-review/presentations/ruti-longitudinal-vf-pipeline-review-v3"
v4 <- "outputs/manual-20260527-current-review/presentations/ruti-longitudinal-vf-pipeline-review-v4"
v5 <- "outputs/manual-20260527-current-review/presentations/ruti-longitudinal-vf-pipeline-review-v5"

canonical_decks <- c(
  "outputs/lecturer_methodology_pack/rUTI_complete_methodology_for_lecturer.pptx",
  "outputs/longcycler_only_methods_summary/Longcycler_only_methods_summary.pptx",
  "outputs/longcycler_only_methods_summary/Longcycler_only_analysis_flowchart.pptx",
  "outputs/longcycler_only_methods_summary/Longcycler_only_methods_summary_with_flowchart.pptx",
  "outputs/manual-20260526-ruti-onboarding/presentations/ruti-clinical-genomic-onboarding/output/ruti-clinical-genomic-onboarding-review.pptx",
  "outputs/manual-20260527-current-review/presentations/ruti-current-results-scientific-review/output/rUTI_Current_Results_Scientific_Review_2026-05-27.pptx",
  "outputs/manual-20260527-current-review/presentations/ruti-current-results-scientific-review-v2/output/rUTI_Current_Results_VF_Focused_Review_2026-05-27.pptx",
  file.path(v3, "output", "Longitudinal_Urinary_Ecoli_VF_Pipeline_Review_2026-05-27.pptx"),
  file.path(v4, "output", "Longitudinal_Urinary_Ecoli_VF_Pipeline_Review_With_Onboarding_2026-05-28.pptx"),
  file.path(v5, "output", "Longitudinal_Urinary_Ecoli_VF_Pipeline_Review_2026-05-28.pptx"),
  file.path(v5, "output", "Longitudinal_Urinary_Ecoli_VF_Pipeline_Review_Compact_Onboarding_2026-05-28.pptx")
)
canonical_deck_slide_counts <- c(8L, 3L, 1L, 4L, 15L, 22L, 22L, 22L, 25L, 17L, 17L)

canonical_companions <- c(
  "outputs/lecturer_methodology_pack/layperson_to_technical_talking_points.md",
  "outputs/lecturer_methodology_pack/methodology_audit_findings.md",
  "outputs/lecturer_methodology_pack/numbered_R_script_methods_register.csv",
  "outputs/lecturer_methodology_pack/numbered_R_script_methods_register.md",
  "outputs/lecturer_methodology_pack/presentation_number_provenance.csv",
  "outputs/lecturer_methodology_pack/rUTI_methods_section_audited.md",
  "outputs/longcycler_only_methods_summary/Longcycler_only_methods_handout.md",
  "outputs/longcycler_only_methods_summary/Longcycler_only_methods_layperson_explanation.md",
  "outputs/longcycler_only_methods_summary/Longcycler_only_methods_talking_points.md",
  "outputs/longcycler_only_methods_summary/flowchart_counts.csv",
  "outputs/longcycler_only_methods_summary/longcycler_only_methods_counts.csv",
  file.path(v3, "output", "Longitudinal_Urinary_Ecoli_VF_Pipeline_Detailed_Presenter_Script_2026-05-28.docx"),
  file.path(v3, "output", "Longitudinal_Urinary_Ecoli_VF_Pipeline_Detailed_Presenter_Script_2026-05-28.md"),
  file.path(v4, "output", "Longitudinal_Urinary_Ecoli_VF_Pipeline_Presenter_Guide_With_Onboarding_2026-05-28.docx"),
  file.path(v4, "output", "Longitudinal_Urinary_Ecoli_VF_Pipeline_Presenter_Guide_With_Onboarding_2026-05-28.md"),
  file.path(v5, "output", "Longitudinal_Urinary_Ecoli_VF_Pipeline_Presenter_Guide_Compact_Onboarding_2026-05-28.docx"),
  file.path(v5, "output", "Longitudinal_Urinary_Ecoli_VF_Pipeline_Presenter_Guide_Compact_Onboarding_2026-05-28.md"),
  "outputs/codebooks/vf_analysis_ready_lay_codebook.docx",
  "outputs/codebooks/vf_analysis_ready_lay_codebook.pdf"
)

release_docs <- c(
  "FOLDER_MAP.md",
  "CODE_REVIEW_RECONCILIATION_README.md",
  "docs/LECTURER_README.md",
  "docs/PIPELINE_FAILURE_LOG.md",
  "docs/VF_abstract_draft.md",
  "docs/VF_merge_diagnostics.md",
  "docs/VF_verification_report.md",
  "docs/figures/timepoint_vs_isolate_clarification.md",
  "docs/legacy_asb_uti_docs/VF_provenance_map.md",
  "docs/workflow_case_count_flowchart.md",
  "docs/workflow_flowchart.md"
)

checks <- list()
add_check <- function(component, file, metric, observed, requirement, pass) {
  checks[[length(checks) + 1L]] <<- tibble(
    checked_at = as.character(Sys.time()),
    component = component,
    file = as.character(file),
    metric = metric,
    observed = as.character(observed),
    requirement = as.character(requirement),
    pass = isTRUE(pass)
  )
}

xml_unescape <- function(x) {
  x <- str_replace_all(x, fixed("&amp;"), "&")
  x <- str_replace_all(x, fixed("&lt;"), "<")
  x <- str_replace_all(x, fixed("&gt;"), ">")
  x <- str_replace_all(x, fixed("&quot;"), '"')
  x <- str_replace_all(x, fixed("&apos;"), "'")
  str_squish(x)
}

extract_tag_text <- function(xml, tag) {
  pattern <- paste0("(?s)<", tag, "(?:\\s[^>]*)?>(.*?)</", tag, ">")
  matches <- str_match_all(xml, regex(pattern))[[1]]
  if (!nrow(matches)) return(character())
  xml_unescape(str_replace_all(matches[, 2], "<[^>]+>", " "))
}

read_zip_member <- function(archive, member, member_size) {
  connection <- unz(archive, member, open = "rb")
  on.exit(close(connection), add = TRUE)
  rawToChar(readBin(connection, what = "raw", n = as.integer(member_size)))
}

office_visible_text <- function(path, extension) {
  listing <- utils::unzip(path, list = TRUE)
  if (extension == "pptx") {
    members <- listing[str_detect(listing$Name, "^ppt/(slides/slide[0-9]+|notesSlides/notesSlide[0-9]+)\\.xml$"), , drop = FALSE]
    tag <- "a:t"
  } else {
    members <- listing[str_detect(listing$Name, "^word/(document|header[0-9]+|footer[0-9]+|footnotes|endnotes|comments)\\.xml$"), , drop = FALSE]
    tag <- "w:t"
  }
  if (!nrow(members)) stop("No visible-text XML members found in ", path, call. = FALSE)
  paste(unlist(Map(
    function(member, member_size) extract_tag_text(read_zip_member(path, member, member_size), tag),
    members$Name,
    members$Length
  )), collapse = "\n")
}

visible_text <- function(path) {
  extension <- str_to_lower(tools::file_ext(path))
  if (extension %in% c("pptx", "docx")) return(office_visible_text(path, extension))
  if (extension == "pdf") {
    if (!nzchar(pdftotext_bin) || !file.exists(pdftotext_bin)) {
      stop("pdftotext is required to verify PDF text.", call. = FALSE)
    }
    output <- system2(pdftotext_bin, c("-layout", path, "-"), stdout = TRUE, stderr = TRUE)
    if (!is.null(attr(output, "status")) && attr(output, "status") != 0L) {
      stop("pdftotext failed for ", path, call. = FALSE)
    }
    return(paste(output, collapse = "\n"))
  }
  paste(readLines(path, warn = FALSE, encoding = "UTF-8"), collapse = "\n")
}

all_required <- c(canonical_decks, canonical_companions, release_docs)
absolute_required <- file.path(root, all_required)
registry_sha256 <- NA_character_
for (index in seq_along(all_required)) {
  path <- absolute_required[[index]]
  exists_nonempty <- file.exists(path) && isTRUE(file.info(path)$size > 0)
  add_check("inventory", all_required[[index]], "exists_nonempty", exists_nonempty, TRUE, exists_nonempty)
}

for (index in seq_along(canonical_decks)) {
  relative <- canonical_decks[[index]]
  path <- file.path(root, relative)
  if (!file.exists(path)) next
  listing <- tryCatch(utils::unzip(path, list = TRUE)$Name, error = function(error) character())
  observed_slides <- sum(str_detect(listing, "^ppt/slides/slide[0-9]+\\.xml$"))
  add_check("deck_structure", relative, "slide_count", observed_slides, canonical_deck_slide_counts[[index]], observed_slides == canonical_deck_slide_counts[[index]])
}

if (!file.exists(registry_path)) {
  add_check("claim_registry", registry_path, "exists", FALSE, TRUE, FALSE)
  registry <- NULL
} else {
  registry <- tryCatch(fromJSON(registry_path, simplifyVector = TRUE), error = function(error) error)
  add_check("claim_registry", registry_path, "readable_json", !inherits(registry, "error"), TRUE, !inherits(registry, "error"))
  if (!inherits(registry, "error")) {
    anchors <- c(
      episodes = registry$analytical_cohort$episodes,
      residents = registry$analytical_cohort$residents,
      operational_UTI = registry$analytical_cohort$operational_UTI,
      operational_Not_UTI = registry$analytical_cohort$operational_Not_UTI,
      source_episodes = registry$attrition_qc_context$episodes,
      source_residents = registry$attrition_qc_context$residents,
      source_operational_UTI = registry$attrition_qc_context$operational_UTI,
      source_operational_Not_UTI = registry$attrition_qc_context$operational_Not_UTI,
      direct_pairs = registry$direct_pairs$all_within_resident,
      adjacent_pairs = registry$adjacent_transitions$pairs,
      adjacent_residents = registry$adjacent_transitions$residents,
      adjacent_threshold = registry$adjacent_transitions$operational_snp_threshold,
      adjacent_le25 = registry$adjacent_transitions$at_or_below_threshold,
      focused_pairs = registry$adjacent_transitions$Not_UTI_to_UTI,
      focused_le25 = registry$adjacent_transitions$Not_UTI_to_UTI_at_or_below_threshold,
      casebook_cases = registry$mechanism_casebook$cases,
      casebook_linked = registry$mechanism_casebook$linked,
      casebook_missing = registry$mechanism_casebook$missing,
      near_miss = registry$near_miss_audit$rows,
      rq_count = registry$research_questions$count
    )
    required_anchors <- c(
      532L, 161L, 16L, 516L,
      583L, 166L, 18L, 565L,
      893L, 371L, 139L, 25L, 140L, 9L, 5L,
      9L, 9L, 0L, 17L, 10L
    )
    add_check("claim_registry", registry_path, "release_anchors", paste(anchors, collapse = "/"), paste(required_anchors, collapse = "/"), identical(as.integer(anchors), required_anchors))
    registry_identity_ok <- identical(as.character(registry$schema_version), "longcycler_release_claim_registry_v1") &&
      identical(as.character(registry$analysis_scope$assembly_policy), "selected QC-passing Longcycler only") &&
      identical(as.character(registry$analysis_scope$clinical_phenotype), "operational UTI phenotype") &&
      identical(as.character(registry$analysis_scope$interpretation), "exploratory observational analysis; no causal claim") &&
      identical(as.character(registry$attrition_qc_context$label), "full clinical source retained only for attrition/QC context") &&
      identical(as.character(registry$near_miss_audit$label), "near-miss rows; not operational UTI cases") &&
      identical(as.character(registry$research_questions$first), "RQ01") &&
      identical(as.character(registry$research_questions$last), "RQ10") &&
      identical(as.integer(registry$research_questions$retired_questions), 0L)
    add_check("claim_registry", registry_path, "identity_and_scope", registry_identity_ok, TRUE, registry_identity_ok)
    excluded_metrics <- as.character(registry$method_contract$assembly_qc$excluded_metrics)
    methods_ok <- identical(as.integer(registry$method_contract$operational_phenotype$culture_lower_bound_cfu_per_ml), 1000L) &&
      identical(as.character(registry$method_contract$operational_phenotype$rule), "versioned operational culture-plus-compatible-symptom phenotype") &&
      identical(as.character(registry$method_contract$operational_phenotype$caveat), "not a reconstruction of the full published protocol") &&
      identical(as.integer(registry$method_contract$assembly_qc$max_contigs), 200L) &&
      identical(as.integer(registry$method_contract$assembly_qc$min_n50_bp), 20000L) &&
      identical(as.integer(registry$method_contract$assembly_qc$min_genome_size_bp), 4000000L) &&
      identical(as.integer(registry$method_contract$assembly_qc$max_genome_size_bp), 6000000L) &&
      length(excluded_metrics) == 3L &&
      setequal(excluded_metrics, c("read coverage", "completeness", "contamination")) &&
      identical(as.character(registry$method_contract$vfdb$tool), "ABRicate") &&
      identical(as.character(registry$method_contract$vfdb$database), "VFDB") &&
      identical(as.integer(registry$method_contract$vfdb$min_identity_pct), 80L) &&
      identical(as.integer(registry$method_contract$vfdb$min_coverage_pct), 80L) &&
      identical(as.character(registry$method_contract$vfdb$provenance), "SHA-bound calls from the selected Longcycler FASTA manifest") &&
      identical(as.character(registry$method_contract$mlst$role), "lineage context; not pair-specific continuity proof") &&
      identical(as.integer(registry$method_contract$mlst$provider_min_good_targets_pct), 95L) &&
      identical(as.character(registry$method_contract$mlst$provider_policy), "provider_qc95 call key/path-linked to the selected Longcycler episode; local fallback excluded") &&
      identical(as.character(registry$method_contract$mlst$fallback), "labelled local MLST from the same selected Longcycler FASTA where required") &&
      identical(as.character(registry$method_contract$direct_pair_evidence$tool), "dnadiff") &&
      identical(as.character(registry$method_contract$direct_pair_evidence$role), "primary pair-specific distance evidence") &&
      identical(as.integer(registry$method_contract$direct_pair_evidence$operational_snp_threshold), 25L) &&
      identical(as.character(registry$method_contract$direct_pair_evidence$priority), "graph connectivity and MLST agreement cannot override a conflicting direct pair") &&
      identical(as.character(registry$method_contract$population_context$core_genome_tool), "Parsnp") &&
      identical(as.character(registry$method_contract$population_context$pangenome_tool), "Panaroo") &&
      identical(as.character(registry$method_contract$population_context$role), "population context; not a substitute for direct pair evidence")
    add_check("claim_registry", registry_path, "method_contract", methods_ok, TRUE, methods_ok)
    registry_sha256 <- unname(digest(registry_path, algo = "sha256", file = TRUE))
    add_check("claim_registry", registry_path, "sha256", registry_sha256, "64 lowercase hexadecimal characters", str_detect(registry_sha256, "^[0-9a-f]{64}$"))
  }
}

forbidden_token <- str_c("fl", "ye")
stale_patterns <- c(
  retired_input = paste0("(?i)", forbidden_token),
  retired_question = "(?i)\\bRQ11\\b",
  stale_556 = "\\b556\\b",
  stale_394 = "\\b394\\b",
  stale_539_not_uti = "(?i)\\b539\\s+(?:operational\\s+)?Not_UTI\\b",
  stale_162_people = "(?i)\\b162\\s+(?:participants?|residents?)\\b",
  stale_17_uti = "(?i)\\b17\\s+(?:operational\\s+)?UTI\\b",
  stale_116 = "\\b116\\b",
  stale_7_of_9 = "(?i)\\b7\\s*(?:/|of)\\s*9\\b"
)

texts <- list()
for (relative in all_required[file.exists(file.path(root, all_required))]) {
  path <- file.path(root, relative)
  text <- tryCatch(visible_text(path), error = function(error) error)
  readable <- !inherits(text, "error")
  add_check("content", relative, "visible_text_readable", readable, TRUE, readable)
  if (!readable) next
  texts[[relative]] <- text
  for (name in names(stale_patterns)) {
    hits <- str_count(text, regex(stale_patterns[[name]]))
    add_check("content", relative, name, hits, 0L, hits == 0L)
  }
  has_source_counts <- str_detect(text, regex("\\b(?:583|166|565)\\b"))
  source_labelled <- str_detect(
    text,
    regex("attrition(?:\\s*/\\s*|\\s+and\\s+)QC", ignore_case = TRUE)
  )
  add_check("content", relative, "full_source_counts_labelled", if (has_source_counts) source_labelled else "not present", "labelled attrition/QC when present", !has_source_counts || source_labelled)
}

registry_bound_text_artifacts <- c(
  canonical_companions[str_detect(canonical_companions, regex("Presenter_(?:Script|Guide)|codebook", ignore_case = TRUE))],
  release_docs
)
if (!is.na(registry_sha256)) {
  for (relative in registry_bound_text_artifacts[file.exists(file.path(root, registry_bound_text_artifacts))]) {
    text <- texts[[relative]]
    if (is.null(text)) next
    has_exact_registry_sha <- str_detect(str_remove_all(text, "\\s+"), fixed(registry_sha256))
    add_check("provenance", relative, "exact_claim_registry_sha256", has_exact_registry_sha, registry_sha256, has_exact_registry_sha)
  }
}

for (relative in canonical_decks[file.exists(file.path(root, canonical_decks))]) {
  text <- texts[[relative]]
  if (is.null(text)) next
  required_visible <- c(
    Longcycler = "(?i)Longcycler",
    episodes_532 = "\\b532\\b",
    residents_161 = "\\b161\\b",
    operational_UTI_16 = "\\b16\\b",
    operational_Not_UTI_516 = "\\b516\\b",
    UTI_label = "(?i)\\bUTI\\b",
    Not_UTI_label = "(?i)\\bNot_UTI\\b"
  )
  missing <- names(required_visible)[!vapply(required_visible, function(pattern) str_detect(text, regex(pattern)), logical(1))]
  add_check("deck_claims", relative, "visible_release_anchors", if (length(missing)) paste(missing, collapse = ",") else "all present", paste(names(required_visible), collapse = ","), !length(missing))
}

codebook_relative <- "outputs/codebooks/vf_analysis_ready_lay_codebook.docx"
if (!is.null(texts[[codebook_relative]])) {
  codebook_ok <- str_detect(texts[[codebook_relative]], regex("\\b532\\b")) && str_detect(texts[[codebook_relative]], regex("Longcycler", ignore_case = TRUE))
  add_check("codebook", codebook_relative, "selected_cohort_anchor", codebook_ok, TRUE, codebook_ok)
}

if (file.exists(registry_path)) {
  registry_time <- file.info(registry_path)$mtime
  for (relative in c(canonical_decks, canonical_companions, release_docs)) {
    path <- file.path(root, relative)
    if (!file.exists(path)) next
    current <- file.info(path)$mtime >= registry_time
    add_check("freshness", relative, "not_older_than_claim_registry", current, TRUE, current)
  }
}

result <- if (length(checks)) do.call(rbind, checks) else tibble()
write_csv(result, report_csv)
failed <- result[!result$pass, , drop = FALSE]
summary <- c(
  "Longcycler-only release deliverable verification",
  paste0("Generated: ", Sys.time()),
  paste0("Canonical decks: ", length(canonical_decks)),
  paste0("Canonical companions: ", length(canonical_companions)),
  paste0("Current release documents: ", length(release_docs)),
  paste0("Checks passed: ", sum(result$pass), "/", nrow(result)),
  if (nrow(failed)) {
    paste0("FAILED: ", failed$component, " / ", failed$file, " / ", failed$metric, " observed=", failed$observed, " required=", failed$requirement)
  } else {
    "PASS: canonical analytical communication artifacts are current and Longcycler-only."
  }
)
writeLines(summary, report_txt)

if (nrow(failed)) {
  stop("Longcycler release deliverable verification failed. See ", report_csv, call. = FALSE)
}
message("Longcycler release deliverable verification passed: ", nrow(result), " checks.")

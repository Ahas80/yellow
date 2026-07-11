#!/usr/bin/env Rscript

source("00_config.R")

suppressPackageStartupMessages({
  library(dplyr)
  library(purrr)
  library(readr)
  library(stringr)
  library(tidyr)
})

out_dir <- file.path("results", "audit")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

missing_tokens <- c(
  "", ".", "?", "-", "--", "na", "n/a", "nan", "missing",
  "unknown", "onbekend", "n.t.b.", "ntb", "nvt", "n.v.t.",
  "niet te koppelen", "nog te koppelen", "still to be linked",
  "not linked", "not known"
)

is_missing_like <- function(x) {
  if (is.logical(x)) {
    return(is.na(x))
  }
  if (is.numeric(x)) {
    return(is.na(x))
  }
  y <- str_squish(str_to_lower(as.character(x)))
  is.na(x) | is.na(y) | y %in% missing_tokens
}

read_if_exists <- function(path) {
  if (!file.exists(path)) {
    return(NULL)
  }
  read_csv(path, show_col_types = FALSE)
}

scope_vec <- function(df, expr = NULL) {
  if (is.null(expr)) {
    rep(TRUE, nrow(df))
  } else {
    out <- eval(substitute(expr), df, parent.frame())
    out %in% TRUE
  }
}

examples_for <- function(df, col, miss) {
  if (!col %in% names(df) || !any(miss)) {
    return(NA_character_)
  }
  key_cols <- intersect(
    c("Participant_id", "tp_lab", "Timepoint", "Isolate_ID", "Assembly_ID", "Episode_ID"),
    names(df)
  )
  if (length(key_cols) == 0) {
    key_cols <- col
  }
  ex <- df[miss, key_cols, drop = FALSE] %>%
    head(6) %>%
    mutate(across(everything(), as.character))
  if (nrow(ex) == 0) {
    NA_character_
  } else {
    apply(ex, 1, function(row) paste(paste(names(row), row, sep = "="), collapse = "; ")) %>%
      paste(collapse = " | ")
  }
}

field_summary <- function(df, layer, file, column, scope_name = "all rows",
                          scope = rep(TRUE, nrow(df)), severity = "info",
                          interpretation = "") {
  n_scope <- sum(scope, na.rm = TRUE)
  if (!column %in% names(df)) {
    return(tibble(
      layer, file, column, scope = scope_name, n_scope,
      missing_like_n = NA_integer_, missing_like_pct = NA_real_,
      severity = "absent_column",
      interpretation = paste("Column is absent.", interpretation),
      examples = NA_character_
    ))
  }
  miss <- rep(FALSE, nrow(df))
  miss[scope] <- is_missing_like(df[[column]][scope])
  tibble(
    layer, file, column, scope = scope_name, n_scope,
    missing_like_n = sum(miss, na.rm = TRUE),
    missing_like_pct = if (n_scope > 0) round(100 * sum(miss, na.rm = TRUE) / n_scope, 1) else NA_real_,
    severity,
    interpretation,
    examples = examples_for(df, column, miss)
  )
}

status_for_count <- function(n, severity) {
  case_when(
    is.na(n) ~ "CHECK",
    n == 0 ~ "OK",
    severity %in% c("blocker", "warning") ~ "CHECK",
    TRUE ~ "INFO"
  )
}

summaries <- list()
notes <- list()

batch_paths <- file.path("data", "inputs", paste0("batch", 1:6, ".csv"))
batch_inputs <- map2(batch_paths, seq_along(batch_paths), function(path, batch_id) {
  x <- read_if_exists(path)
  if (is.null(x)) {
    return(NULL)
  }
  x %>%
    mutate(across(everything(), as.character)) %>%
    mutate(Batch = batch_id, .before = 1)
}) %>% compact() %>% bind_rows()

if (nrow(batch_inputs) > 0) {
  batch_inputs <- apply_manual_sample_curation(batch_inputs, context = "metadata_audit_batch_inputs")
  batch_inputs_all <- batch_inputs
  batch_inputs <- filter_primary_analysis(batch_inputs)
  file <- "data/inputs/batch1.csv..batch6.csv"
  uti_scope <- str_detect(as.character(batch_inputs$Timepoint), "^UTI-")
  symptom_cols <- intersect(
    c("S&S Dysuria", "S&S Urgency", "S&S Frequency", "S&S Incontinence",
      "S&S Pus", "S&S Flank_pain", "S&S Suprapubic pain", "S&S Fever",
      "S&S Chills", "S&S Delirium", "S&S Other", "No S&S"),
    names(batch_inputs)
  )
  symptom_missing_all <- if (length(symptom_cols) > 0) {
    apply(batch_inputs[symptom_cols], 1, function(row) all(is_missing_like(row)))
  } else {
    rep(TRUE, nrow(batch_inputs))
  }

  summaries$batch_inputs <- bind_rows(
    field_summary(batch_inputs, "batch_inputs", file, "Participant_id", severity = "blocker",
                  interpretation = "Participant key used by every downstream join."),
    field_summary(batch_inputs, "batch_inputs", file, "isolate_ID", severity = "blocker",
                  interpretation = "Isolate key used to connect clinical rows to assemblies."),
    field_summary(batch_inputs, "batch_inputs", file, "Collection_Date", severity = "warning",
                  interpretation = "Needed for timeline and episode interpretation."),
    field_summary(batch_inputs, "batch_inputs", file, "Timepoint", severity = "blocker",
                  interpretation = "Used to derive tp_lab and episode keys."),
    field_summary(batch_inputs, "batch_inputs", file, "Organism", severity = "warning",
                  interpretation = "Expected organism metadata from the workbook."),
    field_summary(batch_inputs, "batch_inputs", file, "Beoord", severity = "warning",
                  interpretation = "Growth/category field used as CFU fallback when CFU text is missing."),
    field_summary(batch_inputs, "batch_inputs", file, "CFU_Count", severity = "warning",
                  interpretation = "Culture burden text used for primary culture support."),
    field_summary(batch_inputs, "batch_inputs", file, "Urine collection method", severity = "warning",
                  interpretation = "Controls catheter-aware symptom logic."),
    field_summary(batch_inputs, "batch_inputs", file, "UTI_Label", "UTI-* rows only", uti_scope,
                  severity = "warning", interpretation = "Expected for UTI-labelled episodes, not routine T rows."),
    field_summary(batch_inputs, "batch_inputs", file, "Population", severity = "info",
                  interpretation = "Supplemental cohort descriptor; not used for primary UTI/Not_UTI status."),
    field_summary(batch_inputs, "batch_inputs", file, "UWI#", severity = "info",
                  interpretation = "Supplemental UWI count; not used for primary UTI/Not_UTI status."),
    tibble(
      layer = "batch_inputs", file, column = "all_S&S_columns",
      scope = "all rows", n_scope = nrow(batch_inputs),
      missing_like_n = sum(symptom_missing_all),
      missing_like_pct = round(100 * mean(symptom_missing_all), 1),
      severity = "warning",
      interpretation = "Rows where every detected symptom/No-S&S field is blank or missing-like.",
      examples = examples_for(batch_inputs, "Participant_id", symptom_missing_all)
    )
  )
}

assembly <- read_if_exists("assembly_metadata.csv")
if (!is.null(assembly)) {
  assembly <- apply_manual_sample_curation(assembly, context = "metadata_audit_assembly")
  assembly_all <- assembly
  assembly <- filter_primary_genomics(assembly)
  file <- "assembly_metadata.csv"
  found_scope <- if ("found" %in% names(assembly)) assembly$found %in% TRUE else rep(TRUE, nrow(assembly))
  uti_scope <- str_detect(as.character(assembly$tp_lab), "^UTI-")
  summaries$assembly <- bind_rows(
    field_summary(assembly, "assembly_metadata", file, "Participant_id", severity = "blocker"),
    field_summary(assembly, "assembly_metadata", file, "tp_lab", severity = "blocker"),
    field_summary(assembly, "assembly_metadata", file, "Isolate_ID", severity = "blocker"),
    field_summary(assembly, "assembly_metadata", file, "Batch", severity = "warning"),
    field_summary(assembly, "assembly_metadata", file, "Collection_Date", severity = "warning"),
    field_summary(assembly, "assembly_metadata", file, "Clinical_Organism", severity = "warning"),
    field_summary(assembly, "assembly_metadata", file, "Clinical_Beoord", severity = "warning"),
    field_summary(assembly, "assembly_metadata", file, "Clinical_CFU_Count", severity = "warning"),
    field_summary(assembly, "assembly_metadata", file, "Urine_collection_method", severity = "warning"),
    field_summary(assembly, "assembly_metadata", file, "UTI_Label", "UTI-* rows only", uti_scope,
                  severity = "warning", interpretation = "Expected for UTI-labelled metadata rows."),
    field_summary(assembly, "assembly_metadata", file, "file_name", severity = "info",
                  interpretation = "Blank on expected isolate rows whose FASTA is absent."),
    field_summary(assembly, "assembly_metadata", file, "full_path", severity = "info",
                  interpretation = "Blank on expected isolate rows whose FASTA is absent."),
    field_summary(assembly, "assembly_metadata", file, "file_name", "found FASTA rows only", found_scope,
                  severity = "blocker", interpretation = "Should be populated for all found FASTA rows."),
    field_summary(assembly, "assembly_metadata", file, "full_path", "found FASTA rows only", found_scope,
                  severity = "blocker", interpretation = "Should be populated for all found FASTA rows."),
    field_summary(assembly, "assembly_metadata", file, "num_contigs", "found FASTA rows only", found_scope,
                  severity = "warning"),
    field_summary(assembly, "assembly_metadata", file, "total_bases", "found FASTA rows only", found_scope,
                  severity = "warning"),
    field_summary(assembly, "assembly_metadata", file, "gc_content", "found FASTA rows only", found_scope,
                  severity = "warning"),
    field_summary(assembly, "assembly_metadata", file, "Population", severity = "info",
                  interpretation = "Supplemental participant population field; not always present in the RC overview.")
  )

  missing_expected <- assembly %>%
    filter(!(found %in% TRUE) | !(file_exists %in% TRUE) | !(usable_fasta %in% TRUE)) %>%
    select(any_of(c(
      "Participant_id", "tp_lab", "Timepoint", "Isolate_ID", "Batch",
      "Collection_Date", "Clinical_Beoord", "Clinical_CFU_Count",
      "Urine_collection_method", "metadata_source_status", "found",
      "file_exists", "usable_fasta", "genomics_expected_include",
      "genomics_exclusion_reason"
    )))
  write_csv(missing_expected, file.path(out_dir, "metadata_expected_rows_missing_fasta.csv"))

  quarantined_sequences <- assembly_all %>%
    filter(!(genomics_expected_include %in% TRUE)) %>%
    select(any_of(c(
      "Participant_id", "tp_lab", "Timepoint", "Isolate_ID", "Batch",
      "Collection_Date", "Clinical_Beoord", "Clinical_CFU_Count",
      "Urine_collection_method", "found", "file_exists", "usable_fasta",
      "genomics_expected_include", "genomics_exclusion_reason",
      "manual_curation_note"
    )))
  write_csv(quarantined_sequences, file.path(out_dir, "metadata_quarantined_failed_sequences.csv"))
}

canonical <- read_if_exists("results/qc/canonical_assembly_selection.csv")
if (!is.null(canonical)) {
  canonical <- apply_manual_sample_curation(canonical, context = "metadata_audit_canonical") %>%
    filter_primary_genomics()
  file <- "results/qc/canonical_assembly_selection.csv"
  selected_scope <- if ("selected_canonical" %in% names(canonical)) canonical$selected_canonical %in% TRUE else rep(TRUE, nrow(canonical))
  summaries$canonical <- bind_rows(
    field_summary(canonical, "canonical_assembly_selection", file, "Participant_id", "selected canonical rows", selected_scope, "blocker"),
    field_summary(canonical, "canonical_assembly_selection", file, "tp_lab", "selected canonical rows", selected_scope, "blocker"),
    field_summary(canonical, "canonical_assembly_selection", file, "Isolate_ID", "selected canonical rows", selected_scope, "blocker"),
    field_summary(canonical, "canonical_assembly_selection", file, "full_path", "selected canonical rows", selected_scope, "blocker"),
    field_summary(canonical, "canonical_assembly_selection", file, "QC_PASS", "selected canonical rows", selected_scope, "blocker"),
    field_summary(canonical, "canonical_assembly_selection", file, "QC_REASON", "all rows", rep(TRUE, nrow(canonical)), "info",
                  "Blank is expected for QC-passing rows; populated values explain excluded rows.")
  )
}

status <- read_if_exists("results/clinical/status_map.csv")
if (!is.null(status)) {
  status <- apply_manual_sample_curation(status, context = "metadata_audit_status") %>%
    filter_primary_analysis()
  file <- "results/clinical/status_map.csv"
  not_uti_scope <- status$UTI_Status == "Not_UTI"
  uti_scope <- status$UTI_Status == "UTI"
  uti_event_scope <- str_detect(as.character(status$tp_lab), "^UTI-")
  summaries$status <- bind_rows(
    field_summary(status, "clinical_status", file, "Participant_id", severity = "blocker"),
    field_summary(status, "clinical_status", file, "tp_lab", severity = "blocker"),
    field_summary(status, "clinical_status", file, "Collection_Date", severity = "warning"),
    field_summary(status, "clinical_status", file, "UTI_Status", severity = "blocker"),
    field_summary(status, "clinical_status", file, "UTI_binary", severity = "blocker"),
    field_summary(status, "clinical_status", file, "Not_UTI_subgroup", "Not_UTI rows only", not_uti_scope, "blocker",
                  "Should explain every current Not_UTI row."),
    field_summary(status, "clinical_status", file, "Not_UTI_subgroup", "UTI rows only", uti_scope, "info",
                  "Blank is expected for primary UTI rows."),
    field_summary(status, "clinical_status", file, "Urine_collection_method", severity = "warning"),
    field_summary(status, "clinical_status", file, "urine_collection_method_norm", severity = "warning"),
    field_summary(status, "clinical_status", file, "culture_supports_uti", severity = "warning"),
    field_summary(status, "clinical_status", file, "cfu_threshold_source", severity = "warning"),
    field_summary(status, "clinical_status", file, "symptom_compatible_uti", severity = "warning",
                  interpretation = "NA means symptom rule could not be evaluated, usually because collection method is unknown."),
    field_summary(status, "clinical_status", file, "UTI_Label", "UTI-* rows only", uti_event_scope, "warning")
  )

  clinical_incomplete <- status %>%
    filter(
      is_missing_like(UTI_Status) |
        is_missing_like(UTI_binary) |
        is_missing_like(Participant_id) |
        is_missing_like(tp_lab) |
        Not_UTI_subgroup == "unknown_or_indeterminate" |
        is.na(symptom_compatible_uti) |
        urine_collection_method_norm == "unknown"
    ) %>%
    select(any_of(c(
      "Participant_id", "tp_lab", "Collection_Date", "UTI_Status",
      "UTI_binary", "Not_UTI_subgroup", "Urine_collection_method",
      "urine_collection_method_norm", "catheter_rule", "culture_supports_uti",
      "cfu_threshold_source", "symptom_compatible_uti",
      "UTI_classification_reason", "Episode_ID"
    )))
  write_csv(clinical_incomplete, file.path(out_dir, "clinical_status_incomplete_or_indeterminate_rows.csv"))
}

vf <- read_if_exists("results/vf/vf_analysis_ready.csv")
if (!is.null(vf)) {
  vf <- apply_manual_sample_curation(vf, context = "metadata_audit_vf") %>%
    filter_primary_genomics()
  file <- "results/vf/vf_analysis_ready.csv"
  not_uti_scope <- vf$UTI_Status == "Not_UTI"
  uti_scope <- vf$UTI_Status == "UTI"
  uti_event_scope <- str_detect(as.character(vf$tp_lab), "^UTI-")
  summaries$vf <- bind_rows(
    field_summary(vf, "vf_analysis_ready", file, "Participant_id", severity = "blocker"),
    field_summary(vf, "vf_analysis_ready", file, "tp_lab", severity = "blocker"),
    field_summary(vf, "vf_analysis_ready", file, "UTI_Status", severity = "blocker"),
    field_summary(vf, "vf_analysis_ready", file, "UTI_binary", severity = "blocker"),
    field_summary(vf, "vf_analysis_ready", file, "Not_UTI_subgroup", "Not_UTI rows only", not_uti_scope, "blocker"),
    field_summary(vf, "vf_analysis_ready", file, "Not_UTI_subgroup", "UTI rows only", uti_scope, "info",
                  "Blank is expected for primary UTI rows."),
    field_summary(vf, "vf_analysis_ready", file, "Collection_Date", severity = "warning"),
    field_summary(vf, "vf_analysis_ready", file, "Urine_collection_method", severity = "warning"),
    field_summary(vf, "vf_analysis_ready", file, "culture_supports_uti", severity = "warning"),
    field_summary(vf, "vf_analysis_ready", file, "symptom_compatible_uti", severity = "warning"),
    field_summary(vf, "vf_analysis_ready", file, "ST", severity = "warning",
                  interpretation = "Missing-like ST means neither provider QC95 nor labelled local fallback returned a usable sequence type."),
    field_summary(vf, "vf_analysis_ready", file, "ST_source", severity = "warning",
                  interpretation = "ST_source records whether the active ST came from provider QC95, labelled local fallback, or remains missing."),
    field_summary(vf, "vf_analysis_ready", file, "UTI_Label", "UTI-* rows only", uti_event_scope, "warning")
  )

  vf_incomplete <- vf %>%
    filter(
      is_missing_like(UTI_Status) |
        is_missing_like(UTI_binary) |
        is_missing_like(ST) |
        is.na(symptom_compatible_uti) |
        Not_UTI_subgroup == "unknown_or_indeterminate"
    ) %>%
    select(any_of(c(
      "Participant_id", "tp_lab", "Collection_Date", "UTI_Status",
      "UTI_binary", "Not_UTI_subgroup", "ST", "ST_source", "ST_provider", "ST_local", "Urine_collection_method",
      "culture_supports_uti", "symptom_compatible_uti", "Episode_ID"
    )))
  write_csv(vf_incomplete, file.path(out_dir, "vf_ready_incomplete_metadata_rows.csv"))
}

mlst <- read_if_exists(FILE_MLST_CANONICAL)
if (!is.null(mlst)) {
  mlst <- apply_manual_sample_curation(mlst, context = "metadata_audit_mlst") %>%
    filter_primary_genomics()
  file <- FILE_MLST_CANONICAL
  summaries$mlst <- bind_rows(
    field_summary(mlst, "mlst", file, "Participant_id", severity = "blocker"),
    field_summary(mlst, "mlst", file, "tp_lab", severity = "blocker"),
    field_summary(mlst, "mlst", file, "Isolate_ID", severity = "blocker"),
    field_summary(mlst, "mlst", file, "ST", severity = "warning",
                  interpretation = "Missing-like ST indicates provider-preferred MLST still lacks usable sequence type."),
    field_summary(mlst, "mlst", file, "ST_source", severity = "warning",
                  interpretation = "Active MLST must preserve provider/local provenance."),
    field_summary(mlst, "mlst", file, "n_loci_typed", severity = "warning"),
    field_summary(mlst, "mlst", file, "mlst_complete", severity = "warning")
  )

  mlst_incomplete <- mlst %>%
    filter(
      is_missing_like(ST) |
        is.na(n_loci_typed) |
        n_loci_typed == 0
    ) %>%
    select(any_of(c(
      "Participant_id", "tp_lab", "Isolate_ID", "Assembly_ID", "ST", "ST_source", "ST_provider", "ST_local",
      "n_loci_typed", "mlst_complete", "has_new_allele",
      "ambiguous_call", "Collection_Date", "UTI_Label"
    )))
  write_csv(mlst_incomplete, file.path(out_dir, "mlst_incomplete_rows.csv"))
}

model <- read_if_exists("results/models/model_dataset_denominator.csv")
if (!is.null(model)) {
  model <- apply_manual_sample_curation(model, context = "metadata_audit_model") %>%
    filter_primary_genomics()
  file <- "results/models/model_dataset_denominator.csv"
  summaries$model <- bind_rows(
    field_summary(model, "model_denominator", file, "Participant_id", severity = "blocker"),
    field_summary(model, "model_denominator", file, "Timepoint", severity = "blocker"),
    field_summary(model, "model_denominator", file, "UTI_Status", severity = "blocker"),
    field_summary(model, "model_denominator", file, "UTI_binary", severity = "blocker"),
    field_summary(model, "model_denominator", file, "Not_UTI_subgroup", "Not_UTI rows only", model$UTI_Status == "Not_UTI", "blocker"),
    field_summary(model, "model_denominator", file, "model_interpretation", severity = "warning")
  )
}

gene_gap <- read_if_exists("results/vf/vf_gene_annotation_gap_report.csv")
module_map <- read_if_exists("results/vf/gene_module_map.csv")
if (!is.null(gene_gap)) {
  notes$gene_gap <- gene_gap
}
if (!is.null(module_map)) {
  summaries$gene_modules <- bind_rows(
    field_summary(module_map, "gene_module_metadata", "results/vf/gene_module_map.csv", "Gene", severity = "blocker"),
    field_summary(module_map, "gene_module_metadata", "results/vf/gene_module_map.csv", "module_id", severity = "warning"),
    field_summary(module_map, "gene_module_metadata", "results/vf/gene_module_map.csv", "assignment_confidence", severity = "warning"),
    field_summary(module_map, "gene_module_metadata", "results/vf/gene_module_map.csv", "Category", "genes in VF matrix", module_map$in_vf_matrix %in% TRUE, "warning"),
    field_summary(module_map, "gene_module_metadata", "results/vf/gene_module_map.csv", "Subcategory", "genes in VF matrix", module_map$in_vf_matrix %in% TRUE, "warning")
  )
}

summary_tbl <- bind_rows(summaries) %>%
  mutate(status = status_for_count(missing_like_n, severity)) %>%
  arrange(factor(status, levels = c("CHECK", "INFO", "OK")),
          factor(severity, levels = c("blocker", "warning", "info", "absent_column")),
          layer, column)

write_csv(summary_tbl, file.path(out_dir, "metadata_missing_like_by_layer.csv"))

priority_tbl <- summary_tbl %>%
  filter(status == "CHECK", !is.na(missing_like_n), missing_like_n > 0) %>%
  arrange(factor(severity, levels = c("blocker", "warning", "info")),
          desc(missing_like_n))

md_table <- function(df) {
  if (nrow(df) == 0) {
    return("_None._")
  }
  df <- as.data.frame(df, stringsAsFactors = FALSE)
  df[] <- lapply(df, function(x) ifelse(is.na(x), "", as.character(x)))
  header <- paste0("| ", paste(names(df), collapse = " | "), " |")
  sep <- paste0("| ", paste(rep("---", ncol(df)), collapse = " | "), " |")
  rows <- apply(df, 1, function(x) paste0("| ", paste(x, collapse = " | "), " |"))
  paste(c(header, sep, rows), collapse = "\n")
}

compact_priority <- priority_tbl %>%
  transmute(
    layer, column, scope, missing_like_n, missing_like_pct,
    severity, interpretation
  )

clinical_unknown_n <- if (!is.null(status) && "Not_UTI_subgroup" %in% names(status)) {
  sum(status$Not_UTI_subgroup == "unknown_or_indeterminate", na.rm = TRUE)
} else NA_integer_
clinical_status_missing_n <- if (!is.null(status) && "UTI_Status" %in% names(status)) {
  sum(is_missing_like(status$UTI_Status))
} else NA_integer_
vf_status_missing_n <- if (!is.null(vf) && "UTI_Status" %in% names(vf)) {
  sum(is_missing_like(vf$UTI_Status))
} else NA_integer_
vf_st_missing_n <- if (!is.null(vf) && "ST" %in% names(vf)) {
  sum(is_missing_like(vf$ST))
} else NA_integer_
missing_fasta_n <- if (!is.null(assembly) && "found" %in% names(assembly)) {
  sum(!(assembly$found %in% TRUE), na.rm = TRUE)
} else NA_integer_
quarantined_fasta_n <- if (exists("quarantined_sequences")) {
  nrow(quarantined_sequences)
} else 0L
mlst_incomplete_n <- if (file.exists(file.path(out_dir, "mlst_incomplete_rows.csv"))) {
  nrow(read_csv(file.path(out_dir, "mlst_incomplete_rows.csv"), show_col_types = FALSE))
} else NA_integer_

metric_value <- function(df, metric, default = NA_real_) {
  if (is.null(df)) {
    return(default)
  }
  if (metric %in% names(df) && nrow(df) > 0) {
    return(df[[metric]][[1]])
  }
  if (!all(c("metric", "value") %in% names(df))) {
    return(default)
  }
  val <- df$value[df$metric == metric]
  if (length(val) == 0) default else val[[1]]
}

gene_gap_text <- if (!is.null(gene_gap) && nrow(gene_gap) > 0) {
  vf_gene_columns <- metric_value(gene_gap, "vf_gene_columns")
  genes_in_matrix_not_gene_map <- metric_value(gene_gap, "genes_in_matrix_not_gene_map")
  unassigned_genes_in_matrix <- metric_value(gene_gap, "unassigned_genes_in_matrix")
  if (is.na(unassigned_genes_in_matrix) && !is.null(module_map) && all(c("in_vf_matrix", "assignment_confidence") %in% names(module_map))) {
    sum(module_map$in_vf_matrix %in% TRUE & module_map$assignment_confidence == "Unassigned", na.rm = TRUE)
  }
  low_confidence_genes_in_matrix <- metric_value(gene_gap, "low_confidence_genes_in_matrix")
  if (is.na(low_confidence_genes_in_matrix) && !is.null(module_map) && all(c("in_vf_matrix", "assignment_confidence") %in% names(module_map))) {
    sum(module_map$in_vf_matrix %in% TRUE & module_map$assignment_confidence == "Low", na.rm = TRUE)
  }
  audit_flags_need_review <- metric_value(gene_gap, "audit_flags_need_review")
  if (is.na(audit_flags_need_review) && !is.null(module_map) && all(c("in_vf_matrix", "needs_manual_review") %in% names(module_map))) {
    sum(module_map$in_vf_matrix %in% TRUE & module_map$needs_manual_review %in% TRUE, na.rm = TRUE)
  }
  sprintf(
    "- VF gene annotation layer: %d VF gene columns; %d genes in the matrix are not in `gene_map.csv`; %d matrix genes are unassigned; %d are low confidence; %d audit flags need review.",
    vf_gene_columns,
    genes_in_matrix_not_gene_map,
    unassigned_genes_in_matrix,
    low_confidence_genes_in_matrix,
    audit_flags_need_review
  )
} else {
  "- VF gene annotation layer: gap report not found."
}

report <- c(
  "# Metadata Completeness Audit",
  "",
  "This audit treats `NA`, blank strings, `.`, `?`, `-`, `N/A`, `unknown`, `onbekend`, `n.t.b.`, `niet te koppelen`, and similar placeholders as missing-like values.",
  "",
  "## Bottom line",
  "",
  sprintf("- Primary clinical status is complete: %d missing `UTI_Status` rows in `status_map.csv`.", clinical_status_missing_n),
  sprintf("- VF-ready clinical status is complete: %d missing `UTI_Status` rows in `vf_analysis_ready.csv`.", vf_status_missing_n),
  sprintf("- Clinical indeterminate layer remains: %d status-map rows have `Not_UTI_subgroup == unknown_or_indeterminate`.", clinical_unknown_n),
  sprintf("- Assembly metadata has %d expected isolate row(s) without a found/usable FASTA.", missing_fasta_n),
  sprintf("- Separate quarantine list has %d failed/not-expected sequence row(s) removed from active genomics completeness checks.", quarantined_fasta_n),
  sprintf("- VF-ready MLST layer has %d missing-like `ST` value(s).", vf_st_missing_n),
  sprintf("- MLST detailed output has %d incomplete/non-typable row(s).", mlst_incomplete_n),
  gene_gap_text,
  "",
  "## Priority Checks",
  "",
  md_table(compact_priority),
  "",
  "## Interpretation",
  "",
  "- The UTI/Not_UTI labels themselves are no longer missing in either the clinical or VF-ready denominator.",
  "- The remaining incomplete pieces are not one single status problem. They live in separate layers: expected isolates without FASTA, unknown/indeterminate clinical evidence, MLST typing gaps, and VF gene annotation/module metadata.",
  "- Some blanks are expected by design: `UTI_Label` is normally blank for routine `T*` visits, and `Not_UTI_subgroup` is blank for primary `UTI` rows.",
  "- Rows in `clinical_status_incomplete_or_indeterminate_rows.csv`, `metadata_expected_rows_missing_fasta.csv`, `vf_ready_incomplete_metadata_rows.csv`, and `mlst_incomplete_rows.csv` are the concrete active places to inspect before trying to eliminate remaining unknown-like outputs.",
  "- Rows in `metadata_quarantined_failed_sequences.csv` are intentionally outside the active genomics denominator and should not be interpreted as unresolved active missing FASTA problems.",
  "",
  "## Output files",
  "",
  "- `metadata_missing_like_by_layer.csv`: column-by-column missing-like summary across key metadata layers.",
  "- `clinical_status_incomplete_or_indeterminate_rows.csv`: clinical rows with missing/indeterminate status-rule evidence.",
  "- `metadata_expected_rows_missing_fasta.csv`: expected overview rows that do not have a found/usable FASTA.",
  "- `metadata_quarantined_failed_sequences.csv`: failed/not-expected sequence rows removed from active genomics denominator checks.",
  "- `vf_ready_incomplete_metadata_rows.csv`: VF-ready rows with missing ST, status, or indeterminate clinical metadata.",
  "- `mlst_incomplete_rows.csv`: MLST rows with missing/non-typable ST or incomplete locus calls."
)

write_lines(report, file.path(out_dir, "metadata_completeness_audit.md"))

message("Wrote metadata completeness audit to ", normalizePath(out_dir))
message("Priority checks:")
print(compact_priority)

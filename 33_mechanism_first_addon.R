#!/usr/bin/env Rscript
# ==============================================================================
# 33_mechanism_first_addon.R
# ==============================================================================
#
# GOAL:
#   Build a mechanism-first interpretation layer for the primary UTI vs Not_UTI
#   analysis. This script is descriptive: it explains the 11 clinical
#   Not_UTI -> UTI transitions using host context, strain evidence, VF/module
#   stability, accessory-gene changes, plasmid context, variant annotations, and
#   optional AMR screening if local ABRicate databases are available.
#
# OUTPUT:
#   - results/mechanism/not_uti_to_uti_casebook.csv
#   - results/mechanism/not_uti_to_uti_casebook.md
#   - results/mechanism/transition_mechanism_summary.csv
#   - results/mechanism/host_context_transition_summary.csv
#   - results/mechanism/accessory_gene_transition_changes.csv
#   - plots/mechanism/not_uti_to_uti_case_matrix.png
#   - plots/mechanism/strain_replacement_vs_stability.png
#   - plots/mechanism/host_context_transition_heatmap.png
#   - optional AMR outputs under results/mechanism/
#
# ==============================================================================

source("00_config.R")
suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
  library(stringr)
  library(tibble)
  library(ggplot2)
  library(lubridate)
  library(purrr)
})

msg("Starting 33_mechanism_first_addon.R")

DIR_MECHANISM <- file.path(DIR_RESULTS, "mechanism")
DIR_PLOTS_MECHANISM <- file.path(DIR_PLOTS, "mechanism")
ensure_dir(DIR_MECHANISM)
ensure_dir(DIR_PLOTS_MECHANISM)

MECHANISM_AMR_THREADS <- suppressWarnings(as.integer(Sys.getenv("MECHANISM_AMR_THREADS", "4")))
if (is.na(MECHANISM_AMR_THREADS) || MECHANISM_AMR_THREADS < 1) MECHANISM_AMR_THREADS <- 1L
MECHANISM_AMR_THREADS <- min(MECHANISM_AMR_THREADS, max(1L, CORES_USE))

# ==============================================================================
# HELPERS
# ==============================================================================

require_file <- function(path, label = basename(path)) {
  if (!file.exists(path)) stop("Missing required input for mechanism add-on: ", label, " (", path, ")")
  path
}

safe_read_csv <- function(path, ...) {
  if (!file.exists(path)) return(NULL)
  read_csv(path, show_col_types = FALSE, ...)
}

normalise_tp_label <- function(x) {
  if (exists("normalise_timepoint_preserve_events", mode = "function")) {
    normalise_timepoint_preserve_events(x)
  } else {
    x <- str_trim(as.character(x))
    case_when(
      str_detect(str_to_lower(x), "uricult") ~ "Uricult",
      str_detect(str_to_upper(x), "^T\\d+$") ~ str_to_upper(x),
      str_detect(x, "^\\d+$") ~ paste0("T", x),
      TRUE ~ x
    )
  }
}

parse_date_safe <- function(x) {
  x <- as.character(x)
  out <- suppressWarnings(lubridate::dmy(x))
  idx <- is.na(out)
  if (any(idx)) out[idx] <- suppressWarnings(lubridate::ymd(x[idx]))
  as.Date(out)
}

bool_chr <- function(x) {
  case_when(
    is.na(x) ~ "unknown",
    x %in% TRUE ~ "yes",
    x %in% FALSE ~ "no",
    TRUE ~ as.character(x)
  )
}

clean_list_string <- function(x, max_items = 8) {
  x <- unique(str_trim(as.character(x)))
  x <- x[!is.na(x) & x != ""]
  if (length(x) == 0) return("")
  if (length(x) > max_items) {
    paste0(paste(head(x, max_items), collapse = "; "), sprintf("; ... (+%d)", length(x) - max_items))
  } else {
    paste(x, collapse = "; ")
  }
}

split_semicolon <- function(x) {
  if (length(x) == 0 || is.na(x) || str_trim(as.character(x)) == "") return(character())
  str_split(as.character(x), "\\s*;\\s*", simplify = FALSE)[[1]] %>%
    str_trim() %>%
    discard(~ is.na(.x) || .x == "")
}

presence_from_panaroo_cell <- function(x) {
  !is.na(x) & str_trim(as.character(x)) != ""
}

is_vf_named_gene <- function(gene, non_unique, vf_gene_set) {
  tokens <- paste(gene, non_unique, sep = ";")
  tokens <- str_split(tokens, "[;~]+", simplify = FALSE)[[1]]
  tokens <- str_trim(tokens)
  any(tokens %in% vf_gene_set)
}

plot_theme_mechanism <- function(base_size = 10) {
  theme_bw(base_size = base_size) +
    theme(
      plot.title = element_text(face = "bold"),
      plot.subtitle = element_text(colour = "grey25"),
      plot.caption = element_text(hjust = 0, colour = "grey35", size = base_size - 2),
      legend.position = "bottom"
    )
}

case_label <- function(pid, from_tp, to_tp) {
  paste0(pid, ": ", from_tp, "->", to_tp)
}

make_presence_summary <- function(from_list, to_list) {
  gained <- setdiff(to_list, from_list)
  lost <- setdiff(from_list, to_list)
  stable <- intersect(from_list, to_list)
  list(gained = gained, lost = lost, stable = stable)
}

write_md_table <- function(df, path, max_rows = Inf) {
  if (is.infinite(max_rows)) max_rows <- nrow(df)
  df <- df %>% head(max_rows)
  con <- file(path, open = "wt")
  on.exit(close(con), add = TRUE)
  if (nrow(df) == 0) {
    writeLines("_No rows._", con)
    return(invisible(path))
  }
  header <- paste(names(df), collapse = " | ")
  sep <- paste(rep("---", ncol(df)), collapse = " | ")
  writeLines(paste0("| ", header, " |"), con)
  writeLines(paste0("| ", sep, " |"), con)
  apply(df, 1, function(row) {
    row <- gsub("\\|", "/", as.character(row))
    writeLines(paste0("| ", paste(row, collapse = " | "), " |"), con)
  })
  invisible(path)
}

# ==============================================================================
# LOAD CURRENT OUTPUTS
# ==============================================================================

path_case_summary <- require_file(file.path(DIR_VF, "vf_transition_case_summary.csv"))
path_case_index <- require_file(file.path(DIR_VF, "vf_transition_case_index.csv"))
path_gene_changes <- require_file(file.path(DIR_VF, "vf_transition_gene_changes.csv"))
path_module_changes <- require_file(file.path(DIR_VF, "vf_transition_module_changes.csv"))
path_strain_ctx <- require_file(file.path(DIR_VF, "vf_transition_strain_context.csv"))
path_vf_ready <- require_file(FILE_VF_READY)
path_status <- require_file(FILE_STATUS_MAP)
path_status_poster <- require_file(FILE_STATUS_MAP_POSTER)

case_summary_all <- read_csv(path_case_summary, show_col_types = FALSE) %>%
  mutate(
    Participant_id = as.character(Participant_id),
    from_tp = normalise_tp_label(from_tp),
    to_tp = normalise_tp_label(to_tp),
    from_vf_tp = normalise_tp_label(from_vf_tp),
    to_vf_tp = normalise_tp_label(to_vf_tp)
  )
if (!"snp_strain_context" %in% names(case_summary_all)) {
  case_summary_all <- case_summary_all %>%
    mutate(
      snp_strain_context = case_when(
        is.na(SNPs) ~ "Missing SNP evidence",
        SNPs <= strain_snp_threshold() ~ "Strong same strain",
        TRUE ~ "Above same-strain SNP threshold"
      )
    )
}
if (!"st_lineage_context" %in% names(case_summary_all)) {
  case_summary_all <- case_summary_all %>%
    mutate(
      st_lineage_context = case_when(
        same_ST %in% TRUE ~ "Same ST",
        same_ST %in% FALSE ~ "Different ST",
        TRUE ~ "Missing ST evidence"
      )
    )
}
if (!"pair_interpretation" %in% names(case_summary_all)) {
  case_summary_all <- case_summary_all %>%
    mutate(pair_interpretation = same_strain_evidence)
}

case_index <- read_csv(path_case_index, show_col_types = FALSE) %>%
  mutate(
    Participant_id = as.character(Participant_id),
    from_tp = normalise_tp_label(from_tp),
    to_tp = normalise_tp_label(to_tp)
  )

gene_changes <- read_csv(path_gene_changes, show_col_types = FALSE) %>%
  mutate(Participant_id = as.character(Participant_id))

module_changes <- read_csv(path_module_changes, show_col_types = FALSE) %>%
  mutate(Participant_id = as.character(Participant_id))

strain_ctx <- read_csv(path_strain_ctx, show_col_types = FALSE) %>%
  mutate(Participant_id = as.character(Participant_id))

vf_ready <- read_csv(path_vf_ready, show_col_types = FALSE) %>%
  prefer_primary_uti_status(allow_legacy_fallback = FALSE) %>%
  apply_manual_sample_curation(context = "33_vf_ready") %>%
  filter_primary_genomics() %>%
  mutate(
    Participant_id = as.character(Participant_id),
    tp_lab = normalise_tp_label(tp_lab)
  )

status <- read_csv(path_status, show_col_types = FALSE) %>%
  prefer_primary_uti_status(allow_legacy_fallback = FALSE) %>%
  apply_manual_sample_curation(context = "33_status") %>%
  mutate(
    Participant_id = as.character(Participant_id),
    tp_lab = normalise_tp_label(tp_lab),
    Collection_Date_parsed = parse_date_safe(Collection_Date)
  )

status_poster <- read_csv(path_status_poster, show_col_types = FALSE) %>%
  prefer_primary_uti_status(allow_legacy_fallback = FALSE) %>%
  apply_manual_sample_curation(context = "33_status_poster") %>%
  mutate(
    Participant_id = as.character(Participant_id),
    tp_lab = normalise_tp_label(tp_lab)
  )

msg("Loaded transition, clinical, and VF outputs.")

# ==============================================================================
# HOST CONTEXT AND PLASMID CONTEXT
# ==============================================================================

host_cols <- c(
  "Participant_id", "tp_lab", "Episode_ID", "Collection_Date", "Collection_Date_parsed",
  "UTI_Status", "UTI_binary", "UTI_Label", "Event_type", "Batch",
  "Urine_collection_method", "urine_collection_method_raw", "urine_collection_method_norm",
  "catheter_rule", "culture_supports_uti", "cfu_raw", "cfu_raw_parsed",
  "cfu_lower_bound", "cfu_upper_bound", "cfu_ge_1e3", "cfu_ge_1e4", "cfu_ge_1e5",
  "beoord_cat", "symptom_compatible_uti", "symptom_rule_met",
  "local_urinary_symptom_any", "systemic_symptom_any",
  "dysuria_present", "urgency_present", "frequency_present", "incontinence_present",
  "pus_present", "flankpain_present", "suprapubic_pain_present", "fever_present",
  "rigors_present", "delirium_present", "other_sxs_present",
  "UTI_classification_reason", "analysis_include_primary", "analysis_exclusion_reason",
  "genomics_expected_include", "genomics_exclusion_reason"
)

host_lookup <- status %>%
  select(any_of(host_cols)) %>%
  filter(analysis_include_primary %in% TRUE | is.na(analysis_include_primary)) %>%
  distinct(Participant_id, tp_lab, .keep_all = TRUE)

prefix_df <- function(df, prefix, keys = c("Participant_id", "tp_lab")) {
  rename_with(df, ~ paste0(prefix, .x), -all_of(keys))
}

from_host <- prefix_df(host_lookup, "from_") %>%
  rename(from_tp = tp_lab)
to_host <- prefix_df(host_lookup, "to_") %>%
  rename(to_tp = tp_lab)

plasmid_profile <- safe_read_csv(file.path(DIR_RESULTS, "vf_amr", "vf_plasmid_combined_profile.csv"))
if (!is.null(plasmid_profile)) {
  plasmid_profile <- plasmid_profile %>%
    mutate(
      Participant_id = as.character(Participant_id),
      tp_lab = normalise_tp_label(tp_lab),
      replicon_list = as.character(replicon_list)
    ) %>%
    select(Participant_id, tp_lab, n_replicons, replicon_list, starts_with("has_"))
}

plasmid_change_tbl <- tibble()
if (!is.null(plasmid_profile)) {
  plasmid_change_tbl <- case_summary_all %>%
    select(case_id, Participant_id, from_tp, to_tp, from_vf_tp, to_vf_tp, transition_type) %>%
    left_join(
      plasmid_profile %>% rename(from_vf_tp = tp_lab, from_n_replicons = n_replicons,
                                 from_replicon_list = replicon_list),
      by = c("Participant_id", "from_vf_tp"),
      relationship = "many-to-one"
    ) %>%
    left_join(
      plasmid_profile %>% rename(to_vf_tp = tp_lab, to_n_replicons = n_replicons,
                                 to_replicon_list = replicon_list),
      by = c("Participant_id", "to_vf_tp"),
      relationship = "many-to-one"
    ) %>%
    rowwise() %>%
    mutate(
      from_replicons = list(split_semicolon(from_replicon_list)),
      to_replicons = list(split_semicolon(to_replicon_list)),
      plasmid_gained = clean_list_string(setdiff(to_replicons, from_replicons), 12),
      plasmid_lost = clean_list_string(setdiff(from_replicons, to_replicons), 12),
      n_plasmid_gained = length(setdiff(to_replicons, from_replicons)),
      n_plasmid_lost = length(setdiff(from_replicons, to_replicons))
    ) %>%
    ungroup() %>%
    select(case_id, from_n_replicons, to_n_replicons,
           n_plasmid_gained, n_plasmid_lost, plasmid_gained, plasmid_lost)
}

# ==============================================================================
# VARIANT SUMMARY FOR CANDIDATE PAIRS
# ==============================================================================

variant_path <- file.path(DIR_RESULTS, "longitudinal", "variant_annotation_detailed.csv")
variant_summary <- tibble()
if (file.exists(variant_path)) {
  variant_summary <- read_csv(variant_path, show_col_types = FALSE) %>%
    mutate(
      Participant_id = as.character(Participant_id),
      From_Time = normalise_tp_label(From_Time),
      To_Time = normalise_tp_label(To_Time),
      Region = as.character(Region),
      Gene = as.character(Gene),
      Product = as.character(Product)
    ) %>%
    group_by(Participant_id, From_Time, To_Time) %>%
    summarise(
      n_annotated_variants = n(),
      n_coding_or_gene_variants = sum(!is.na(Gene) & Gene != "" & Gene != "NA"),
      variant_genes = clean_list_string(Gene[!is.na(Gene) & Gene != "" & Gene != "NA"], 10),
      variant_products = clean_list_string(Product[!is.na(Product) & Product != "" & Product != "NA"], 6),
      .groups = "drop"
    ) %>%
    rename(from_tp = From_Time, to_tp = To_Time)
}

vf_transition_metrics_path <- file.path(DIR_VF, "vf_longitudinal_transitions.csv")
vf_transition_metrics <- tibble()
if (file.exists(vf_transition_metrics_path)) {
  vf_transition_metrics <- read_csv(vf_transition_metrics_path, show_col_types = FALSE) %>%
    mutate(
      Participant_id = as.character(Participant_id),
      from_tp = normalise_tp_label(tp_from),
      to_tp = normalise_tp_label(tp_to)
    ) %>%
    select(
      Participant_id, from_tp, to_tp,
      vf_jaccard_longitudinal = jaccard_similarity,
      vf_longitudinal_n_gained = n_gained,
      vf_longitudinal_n_lost = n_lost
    )
}

# ==============================================================================
# WHOLE-ACCESSORY-GENOME CHANGES FROM PANAROO
# ==============================================================================

accessory_gene_changes <- tibble()
path_panaroo_roary <- file.path(DIR_WGS_PAN, "gene_presence_absence_roary.csv")
path_panaroo_manifest <- file.path(DIR_WGS_PAN, "panaroo_input_manifest.csv")

if (file.exists(path_panaroo_roary) && file.exists(path_panaroo_manifest)) {
  msg("Computing non-VF accessory-gene changes from Panaroo outputs.")
  pan_manifest <- read_csv(path_panaroo_manifest, show_col_types = FALSE) %>%
    mutate(
      Participant_id = as.character(Participant_id),
      tp_lab = normalise_tp_label(tp_lab),
      Assembly_Base_ID = as.character(Assembly_Base_ID)
    )

  panaroo <- read_csv(path_panaroo_roary, show_col_types = FALSE)
  sample_cols <- intersect(pan_manifest$Assembly_Base_ID, names(panaroo))
  if (length(sample_cols) == 0) {
    warning("No Panaroo sample columns matched panaroo_input_manifest.csv; accessory change output will be empty.")
  } else {
    n_samples <- length(sample_cols)
    if (!"No. isolates" %in% names(panaroo)) {
      panaroo$`No. isolates` <- rowSums(as.data.frame(lapply(panaroo[sample_cols], presence_from_panaroo_cell)))
    }
    panaroo <- panaroo %>%
      mutate(
        pangenome_gene = as.character(Gene),
        non_unique_gene_name = as.character(`Non-unique Gene name`),
        annotation = as.character(Annotation),
        n_isolates = suppressWarnings(as.numeric(`No. isolates`)),
        prevalence = n_isolates / n_samples
      )

    vf_gene_set <- unique(c(
      canonical_vf_gene_cols(names(vf_ready), required = FALSE),
      if (file.exists(file.path(DIR_VF, "gene_module_map.csv"))) {
        read_csv(file.path(DIR_VF, "gene_module_map.csv"), show_col_types = FALSE)$Gene
      } else character()
    ))
    vf_gene_set <- vf_gene_set[!is.na(vf_gene_set) & vf_gene_set != ""]

    panaroo$is_vf_named <- map2_lgl(
      panaroo$pangenome_gene,
      panaroo$non_unique_gene_name,
      is_vf_named_gene,
      vf_gene_set = vf_gene_set
    )
    panaroo$is_accessory <- panaroo$prevalence < 0.95

    sample_lookup <- pan_manifest %>%
      distinct(Participant_id, tp_lab, .keep_all = TRUE) %>%
      select(Participant_id, tp_lab, Assembly_Base_ID)

    not_uti_uti_cases <- case_summary_all %>%
      filter(is_not_uti_to_uti %in% TRUE) %>%
      select(case_id, Participant_id, from_tp, to_tp, from_vf_tp, to_vf_tp,
             transition_type, snp_strain_context, st_lineage_context,
             pair_interpretation, same_strain_evidence)

    for (i in seq_len(nrow(not_uti_uti_cases))) {
      cr <- not_uti_uti_cases[i, ]
      from_sample <- sample_lookup %>%
        filter(Participant_id == cr$Participant_id, tp_lab == cr$from_vf_tp) %>%
        pull(Assembly_Base_ID) %>%
        first()
      to_sample <- sample_lookup %>%
        filter(Participant_id == cr$Participant_id, tp_lab == cr$to_vf_tp) %>%
        pull(Assembly_Base_ID) %>%
        first()

      if (length(from_sample) == 0 || length(to_sample) == 0 ||
          is.na(from_sample) || is.na(to_sample) ||
          !from_sample %in% names(panaroo) || !to_sample %in% names(panaroo)) {
        next
      }

      from_present <- presence_from_panaroo_cell(panaroo[[from_sample]])
      to_present <- presence_from_panaroo_cell(panaroo[[to_sample]])
      changed <- from_present != to_present
      keep <- changed & panaroo$is_accessory & !panaroo$is_vf_named

      if (any(keep, na.rm = TRUE)) {
        from_present_vec <- as.integer(from_present[keep])
        to_present_vec <- as.integer(to_present[keep])
        change_type_vec <- ifelse(to_present[keep], "gained", "lost")

        accessory_gene_changes <- bind_rows(
          accessory_gene_changes,
          tibble(
            case_id = cr$case_id,
            Participant_id = cr$Participant_id,
            from_tp = cr$from_tp,
            to_tp = cr$to_tp,
            from_sample = from_sample,
            to_sample = to_sample,
            pangenome_gene = panaroo$pangenome_gene[keep],
            non_unique_gene_name = panaroo$non_unique_gene_name[keep],
            annotation = panaroo$annotation[keep],
            n_isolates = panaroo$n_isolates[keep],
            prevalence = round(panaroo$prevalence[keep], 4),
            from_present = from_present_vec,
            to_present = to_present_vec,
            change_type = change_type_vec,
            snp_strain_context = cr$snp_strain_context,
            st_lineage_context = cr$st_lineage_context,
            pair_interpretation = cr$pair_interpretation,
            same_strain_evidence = cr$same_strain_evidence,
            accessory_definition = "Panaroo prevalence <95%; VF-named genes excluded"
          )
        )
      }
    }
  }
} else {
  warning("Panaroo gene_presence_absence_roary.csv or panaroo_input_manifest.csv missing; accessory output will be empty.")
}

write_csv(accessory_gene_changes, file.path(DIR_MECHANISM, "accessory_gene_transition_changes.csv"))

accessory_summary <- accessory_gene_changes %>%
  group_by(case_id) %>%
  summarise(
    n_accessory_genes_gained = sum(change_type == "gained", na.rm = TRUE),
    n_accessory_genes_lost = sum(change_type == "lost", na.rm = TRUE),
    accessory_genes_gained = clean_list_string(pangenome_gene[change_type == "gained"], 10),
    accessory_genes_lost = clean_list_string(pangenome_gene[change_type == "lost"], 10),
    .groups = "drop"
  )

# ==============================================================================
# OPTIONAL AMR SCREENING WITH ABRICATE RESFINDER
# ==============================================================================

run_abricate_resfinder <- function() {
  report <- character()
  add <- function(...) report <<- c(report, sprintf(...))
  report_path <- file.path(DIR_MECHANISM, "amr_screen_report.txt")
  cache_dir <- file.path(DIR_MECHANISM, "abricate_resfinder_cache")
  ensure_dir(cache_dir)

  add("=== Optional AMR Screen ===")
  add("Generated: %s", format(Sys.time()))
  add("Tool preference: ABRicate ResFinder, cached per FASTA")
  add("")

  abricate_bin <- Sys.which("abricate")
  if (abricate_bin == "") {
    add("Status: SKIPPED")
    add("Reason: abricate not found in PATH.")
    writeLines(report, report_path)
    return(list(status = "skipped", summary = tibble(), transition = tibble()))
  }

  db_list <- tryCatch(system2(abricate_bin, "--list", stdout = TRUE, stderr = TRUE), error = function(e) character())
  if (!any(str_detect(db_list, "^resfinder\\s"))) {
    add("Status: SKIPPED")
    add("Reason: ABRicate database 'resfinder' not found.")
    add("Available database output:")
    add("%s", paste(db_list, collapse = "\n"))
    writeLines(report, report_path)
    return(list(status = "skipped", summary = tibble(), transition = tibble()))
  }

  selection_file <- file.path(DIR_QC, "canonical_assembly_selection.csv")
  if (!file.exists(selection_file)) {
    add("Status: SKIPPED")
    add("Reason: canonical_assembly_selection.csv not found.")
    writeLines(report, report_path)
    return(list(status = "skipped", summary = tibble(), transition = tibble()))
  }

  manifest <- read_csv(selection_file, show_col_types = FALSE) %>%
    mutate(
      Participant_id = as.character(Participant_id),
      tp_lab = normalise_tp_label(tp_lab),
      selected_canonical = as_pipeline_bool(selected_canonical),
      file_path = if ("full_path" %in% names(.)) as.character(full_path) else as.character(fasta_path),
      file_exists = file.exists(file_path)
    ) %>%
    filter(selected_canonical %in% TRUE, file_exists %in% TRUE) %>%
    distinct(Participant_id, tp_lab, .keep_all = TRUE)

  add("Status: RUNNING_OR_REUSING_CACHE")
  add("ABRicate: %s", abricate_bin)
  add("Database: resfinder")
  add("Canonical FASTAs: %d", nrow(manifest))
  add("Threads: %d", MECHANISM_AMR_THREADS)

  run_one <- function(path) {
    key <- tools::file_path_sans_ext(basename(path))
    out <- file.path(cache_dir, paste0(key, ".resfinder.tab"))
    if (file.exists(out) && file.size(out) > 0) return(out)
    res <- tryCatch(
      system2(
        abricate_bin,
        c("--quiet", "--threads", "1", "--db", "resfinder", path),
        stdout = TRUE,
        stderr = TRUE
      ),
      error = function(e) structure(character(), status = 1L)
    )
    status <- attr(res, "status")
    if (is.null(status)) status <- 0L
    if (status == 0L) {
      writeLines(res, out)
      out
    } else {
      fail <- file.path(cache_dir, paste0(key, ".resfinder.failed.txt"))
      writeLines(res, fail)
      NA_character_
    }
  }

  tab_files <- if (.Platform$OS.type == "unix" && MECHANISM_AMR_THREADS > 1) {
    unlist(parallel::mclapply(manifest$file_path, run_one, mc.cores = MECHANISM_AMR_THREADS), use.names = FALSE)
  } else {
    vapply(manifest$file_path, run_one, character(1))
  }

  manifest <- manifest %>%
    mutate(tab_file = tab_files) %>%
    filter(!is.na(tab_file), file.exists(tab_file))

  read_one <- function(tab_file, Participant_id, tp_lab, Assembly_ID, Isolate_ID) {
    if (!file.exists(tab_file) || file.size(tab_file) == 0) return(tibble())
    raw <- suppressWarnings(read_tsv(tab_file, show_col_types = FALSE, col_types = cols(.default = "c")))
    if (nrow(raw) == 0 || !"GENE" %in% names(raw)) return(tibble())
    raw %>%
      mutate(
        Participant_id = Participant_id,
        tp_lab = tp_lab,
        Assembly_ID = Assembly_ID,
        Isolate_ID = Isolate_ID,
        identity = suppressWarnings(parse_number(`%IDENTITY`)),
        coverage = suppressWarnings(parse_number(`%COVERAGE`)),
        resistance_class = if ("RESISTANCE" %in% names(.)) RESISTANCE else NA_character_,
        product = if ("PRODUCT" %in% names(.)) PRODUCT else NA_character_
      ) %>%
      filter(is.na(identity) | identity >= 80, is.na(coverage) | coverage >= 80) %>%
      select(Participant_id, tp_lab, Assembly_ID, Isolate_ID,
             gene = GENE, resistance_class, product, identity, coverage,
             accession = ACCESSION, database = DATABASE)
  }

  hits <- pmap_dfr(
    manifest %>% select(tab_file, Participant_id, tp_lab, Assembly_ID, Isolate_ID),
    read_one
  ) %>%
    distinct()

  write_csv(hits, file.path(DIR_MECHANISM, "amr_resfinder_hits_long.csv"))

  episode_summary <- if (nrow(hits) > 0) {
    hits %>%
      group_by(Participant_id, tp_lab) %>%
      summarise(
        amr_gene_count = n_distinct(gene),
        amr_class_count = n_distinct(resistance_class[!is.na(resistance_class) & resistance_class != ""]),
        amr_genes = clean_list_string(gene, 20),
        amr_classes = clean_list_string(resistance_class, 12),
        .groups = "drop"
      )
  } else {
    tibble(Participant_id = character(), tp_lab = character(),
           amr_gene_count = integer(), amr_class_count = integer(),
           amr_genes = character(), amr_classes = character())
  }

  write_csv(episode_summary, file.path(DIR_MECHANISM, "amr_presence_by_episode.csv"))

  transition_amr <- case_summary_all %>%
    filter(transition_type %in% c("Not_UTI->UTI", "UTI->Not_UTI")) %>%
    select(case_id, Participant_id, transition_type, from_vf_tp, to_vf_tp) %>%
    rowwise() %>%
    mutate(
      from_amr = list(episode_summary %>%
                        filter(Participant_id == .env$Participant_id, tp_lab == .env$from_vf_tp) %>%
                        pull(amr_genes) %>%
                        split_semicolon()),
      to_amr = list(episode_summary %>%
                      filter(Participant_id == .env$Participant_id, tp_lab == .env$to_vf_tp) %>%
                      pull(amr_genes) %>%
                      split_semicolon()),
      n_amr_gained = length(setdiff(to_amr, from_amr)),
      n_amr_lost = length(setdiff(from_amr, to_amr)),
      amr_gained = clean_list_string(setdiff(to_amr, from_amr), 12),
      amr_lost = clean_list_string(setdiff(from_amr, to_amr), 12)
    ) %>%
    ungroup() %>%
    select(case_id, transition_type, n_amr_gained, n_amr_lost, amr_gained, amr_lost)

  write_csv(transition_amr, file.path(DIR_MECHANISM, "amr_transition_changes.csv"))

  add("Status: COMPLETE")
  add("Successful cached/run tab files: %d", nrow(manifest))
  add("Hit rows after identity/coverage filtering: %d", nrow(hits))
  add("Episodes with AMR hits: %d", nrow(episode_summary))
  add("Transition AMR summaries: %d", nrow(transition_amr))
  writeLines(report, report_path)
  list(status = "complete", summary = episode_summary, transition = transition_amr)
}

amr_result <- run_abricate_resfinder()
amr_episode_summary <- amr_result$summary
amr_transition <- amr_result$transition

# ==============================================================================
# CASEBOOK
# ==============================================================================

not_uti_to_uti <- case_summary_all %>%
  filter(is_not_uti_to_uti %in% TRUE) %>%
  left_join(from_host, by = c("Participant_id", "from_tp"), relationship = "many-to-one") %>%
  left_join(to_host, by = c("Participant_id", "to_tp"), relationship = "many-to-one") %>%
  left_join(plasmid_change_tbl, by = "case_id", relationship = "many-to-one") %>%
  left_join(accessory_summary, by = "case_id", relationship = "many-to-one") %>%
  left_join(variant_summary, by = c("Participant_id", "from_tp", "to_tp"), relationship = "many-to-one") %>%
  left_join(vf_transition_metrics, by = c("Participant_id", "from_tp", "to_tp"), relationship = "many-to-one") %>%
  left_join(amr_transition, by = c("case_id", "transition_type"), relationship = "many-to-one") %>%
  mutate(
    across(c(n_plasmid_gained, n_plasmid_lost,
             n_accessory_genes_gained, n_accessory_genes_lost,
             n_annotated_variants, n_coding_or_gene_variants,
             n_amr_gained, n_amr_lost), ~ replace_na(as.integer(.x), 0L)),
    from_collection_date_parsed = parse_date_safe(from_Collection_Date),
    to_collection_date_parsed = parse_date_safe(to_Collection_Date),
    days_between = as.integer(to_collection_date_parsed - from_collection_date_parsed),
    collection_method_changed = !is.na(from_urine_collection_method_norm) &
      !is.na(to_urine_collection_method_norm) &
      from_urine_collection_method_norm != to_urine_collection_method_norm,
    catheter_rule_changed = !is.na(from_catheter_rule) &
      !is.na(to_catheter_rule) &
      from_catheter_rule != to_catheter_rule,
    symptom_state_change = case_when(
      from_symptom_compatible_uti %in% FALSE & to_symptom_compatible_uti %in% TRUE ~ "symptoms_emerged",
      from_symptom_compatible_uti %in% TRUE & to_symptom_compatible_uti %in% TRUE ~ "symptoms_persisted",
      TRUE ~ "not_clear"
    ),
    culture_state_change = case_when(
      from_culture_supports_uti %in% FALSE & to_culture_supports_uti %in% TRUE ~ "culture_support_emerged",
      from_culture_supports_uti %in% TRUE & to_culture_supports_uti %in% TRUE ~ "culture_support_persisted",
      TRUE ~ "not_clear"
    ),
    clinical_trigger = case_when(
      from_culture_supports_uti %in% TRUE & from_symptom_compatible_uti %in% FALSE &
        to_culture_supports_uti %in% TRUE & to_symptom_compatible_uti %in% TRUE ~
        "symptom emergence on culture-supported bacteriuria",
      from_culture_supports_uti %in% FALSE & to_culture_supports_uti %in% TRUE &
        to_symptom_compatible_uti %in% TRUE ~
        "culture support plus symptoms emerged",
      to_culture_supports_uti %in% TRUE & to_symptom_compatible_uti %in% TRUE ~
        "UTI criteria present at endpoint",
      TRUE ~ "review clinical rule inputs"
    ),
    accessory_change_count = replace_na(n_accessory_genes_gained, 0L) + replace_na(n_accessory_genes_lost, 0L),
    profile_change_count = replace_na(n_vf_genes_gained, 0L) + replace_na(n_vf_genes_lost, 0L) +
      replace_na(n_modules_gained, 0L) + replace_na(n_modules_lost, 0L) +
      replace_na(n_plasmid_gained, 0L) + replace_na(n_plasmid_lost, 0L) +
      replace_na(n_amr_gained, 0L) + replace_na(n_amr_lost, 0L),
    mechanism_bucket = case_when(
      !has_vf_pair ~ "missing_wgs_endpoint",
      pair_interpretation == "Replacement likely" ~ "strain_replacement",
      snp_strain_context == "Strong same strain" & profile_change_count == 0 ~ "same_strain_stable_profile",
      snp_strain_context == "Strong same strain" & profile_change_count > 0 ~ "same_strain_genomic_change",
      TRUE ~ "uncertain"
    ),
    vf_jaccard = coalesce(as.numeric(vf_jaccard), as.numeric(vf_jaccard_longitudinal)),
    mechanism_interpretation = case_when(
      mechanism_bucket == "same_strain_stable_profile" ~
        "Same-strain UTI transition with no detected VF/module/plasmid/AMR presence-absence change; host state or expression/regulation is plausible. Panaroo accessory changes are listed separately for review.",
      mechanism_bucket == "same_strain_genomic_change" ~
        "Same-strain UTI transition with detectable gene/module/plasmid/AMR profile changes; inspect changed components.",
      mechanism_bucket == "strain_replacement" ~
        "UTI transition is more consistent with a new/replacement strain than within-strain VF evolution.",
      mechanism_bucket == "missing_wgs_endpoint" ~
        "Clinical transition retained, but one endpoint lacks WGS/VF data.",
      TRUE ~
        "Strain or genomic interpretation remains uncertain with current thresholds/data."
    )
  ) %>%
  arrange(case_id)

casebook <- not_uti_to_uti %>%
  transmute(
    case_id,
    Participant_id,
    from_tp,
    to_tp,
    days_between,
    from_status,
    to_status,
    mechanism_bucket,
    mechanism_interpretation,
    has_vf_pair,
    ST_from,
    ST_to,
    same_ST,
    SNPs,
    snp_strain_context,
    st_lineage_context,
    pair_interpretation,
    same_strain_evidence,
    vf_jaccard,
    module_jaccard,
    n_vf_genes_gained,
    n_vf_genes_lost,
    n_modules_gained,
    n_modules_lost,
    n_accessory_genes_gained,
    n_accessory_genes_lost,
    accessory_genes_gained,
    accessory_genes_lost,
    from_n_replicons,
    to_n_replicons,
    n_plasmid_gained,
    n_plasmid_lost,
    plasmid_gained,
    plasmid_lost,
    n_amr_gained,
    n_amr_lost,
    amr_gained,
    amr_lost,
    n_annotated_variants,
    n_coding_or_gene_variants,
    variant_genes,
    variant_products,
    from_collection_date = from_Collection_Date,
    to_collection_date = to_Collection_Date,
    from_collection_method = from_Urine_collection_method,
    to_collection_method = to_Urine_collection_method,
    from_collection_method_norm = from_urine_collection_method_norm,
    to_collection_method_norm = to_urine_collection_method_norm,
    collection_method_changed,
    from_catheter_rule,
    to_catheter_rule,
    catheter_rule_changed,
    from_cfu_raw,
    to_cfu_raw,
    from_cfu_lower_bound,
    to_cfu_lower_bound,
    from_beoord_cat,
    to_beoord_cat,
    from_culture_supports_uti,
    to_culture_supports_uti,
    from_symptom_compatible_uti,
    to_symptom_compatible_uti,
    symptom_state_change,
    culture_state_change,
    clinical_trigger,
    to_symptom_rule_met,
    to_dysuria_present,
    to_urgency_present,
    to_frequency_present,
    to_incontinence_present,
    to_pus_present,
    to_flankpain_present,
    to_suprapubic_pain_present,
    to_fever_present,
    to_rigors_present,
    to_delirium_present,
    to_other_sxs_present,
    is_uricult_transition,
    timing_caveat,
    missing_data_note
  )

write_csv(casebook, file.path(DIR_MECHANISM, "not_uti_to_uti_casebook.csv"))

# ==============================================================================
# SUMMARIES
# ==============================================================================

classify_transition_mechanism <- function(df) {
  df %>%
    mutate(
      profile_change_count = replace_na(n_vf_genes_gained, 0L) + replace_na(n_vf_genes_lost, 0L) +
        replace_na(n_modules_gained, 0L) + replace_na(n_modules_lost, 0L),
      mechanism_bucket = case_when(
        !has_vf_pair ~ "missing_wgs_endpoint",
        pair_interpretation == "Replacement likely" ~ "strain_replacement",
        snp_strain_context == "Strong same strain" & profile_change_count == 0 ~ "same_strain_stable_profile",
        snp_strain_context == "Strong same strain" & profile_change_count > 0 ~ "same_strain_genomic_change",
        TRUE ~ "uncertain"
      )
    )
}

transition_mechanism_detail <- classify_transition_mechanism(case_summary_all) %>%
  left_join(vf_transition_metrics, by = c("Participant_id", "from_tp", "to_tp"), relationship = "many-to-one") %>%
  mutate(vf_jaccard_for_summary = coalesce(as.numeric(vf_jaccard), as.numeric(vf_jaccard_longitudinal)))

transition_mechanism_summary <- transition_mechanism_detail %>%
  group_by(transition_type, mechanism_bucket) %>%
  summarise(
    n_transitions = n(),
    n_participants = n_distinct(Participant_id),
    n_wgs_vf_linked = sum(has_vf_pair %in% TRUE, na.rm = TRUE),
    median_snp = ifelse(all(is.na(SNPs)), NA_real_, median(SNPs, na.rm = TRUE)),
    median_vf_jaccard = ifelse(all(is.na(vf_jaccard_for_summary)), NA_real_, median(vf_jaccard_for_summary, na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  arrange(transition_type, mechanism_bucket)

write_csv(transition_mechanism_summary, file.path(DIR_MECHANISM, "transition_mechanism_summary.csv"))

host_context_transition_summary <- casebook %>%
  group_by(mechanism_bucket, clinical_trigger, symptom_state_change, culture_state_change) %>%
  summarise(
    n_cases = n(),
    n_participants = n_distinct(Participant_id),
    n_collection_method_changed = sum(collection_method_changed %in% TRUE, na.rm = TRUE),
    n_catheter_rule_changed = sum(catheter_rule_changed %in% TRUE, na.rm = TRUE),
    median_days_between = ifelse(all(is.na(days_between)), NA_real_, median(days_between, na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  arrange(desc(n_cases), mechanism_bucket)

write_csv(host_context_transition_summary, file.path(DIR_MECHANISM, "host_context_transition_summary.csv"))

# Near-miss sensitivity appendix
near_miss_path <- file.path(DIR_RESULTS, "audit", "uti_not_uti_near_miss_rows.csv")
near_miss_summary <- tibble()
if (file.exists(near_miss_path)) {
  near_miss <- read_csv(near_miss_path, show_col_types = FALSE) %>%
    mutate(
      Participant_id = as.character(Participant_id),
      tp_lab = normalise_tp_label(tp_lab)
    )
  near_miss_summary <- near_miss %>%
    count(near_miss_reason, symptom_rule_met, Urine_collection_method, catheter_rule,
          name = "n_rows", sort = TRUE)
  write_csv(near_miss_summary, file.path(DIR_MECHANISM, "near_miss_sensitivity_summary.csv"))
  writeLines(c(
    "# Near-Miss Sensitivity Appendix",
    "",
    sprintf("Near-miss rows under the current primary rule: %d.", nrow(near_miss)),
    "",
    "These rows are culture-supported or legacy UTI-like but remain primary Not_UTI because the catheter-aware symptom rule is not met or another primary-rule component is missing.",
    "They are not relabelled by this add-on. The mechanism casebook stays anchored to the primary denominator.",
    "",
    "Interpretation if treated as possible UTI: mechanism claims would become broader and less specific, because these rows mostly test clinical-rule sensitivity rather than bacterial gain/loss mechanisms.",
    "",
    "Summary table: `results/mechanism/near_miss_sensitivity_summary.csv`."
  ), file.path(DIR_MECHANISM, "near_miss_sensitivity_appendix.md"))
}

# ==============================================================================
# MARKDOWN CASEBOOK
# ==============================================================================

casebook_md <- file.path(DIR_MECHANISM, "not_uti_to_uti_casebook.md")
md_lines <- c(
  "# Not_UTI -> UTI Mechanism Casebook",
  "",
  sprintf("Generated: %s", format(Sys.time())),
  "",
  "This is a descriptive interpretation layer for the current primary UTI rule.",
  "It does not relabel near-miss rows and does not make confirmatory association claims.",
  "",
  "## Headline Counts",
  "",
  sprintf("- Clinical Not_UTI -> UTI transitions: %d", nrow(casebook)),
  sprintf("- WGS/VF-linked transitions: %d", sum(casebook$has_vf_pair %in% TRUE, na.rm = TRUE)),
  sprintf("- Missing WGS/VF endpoint transitions: %d", sum(casebook$mechanism_bucket == "missing_wgs_endpoint", na.rm = TRUE)),
  sprintf("- Uricult-linked transitions: %d", sum(casebook$is_uricult_transition %in% TRUE, na.rm = TRUE)),
  "",
  "## Mechanism Bucket Counts",
  ""
)

bucket_counts <- casebook %>%
  count(mechanism_bucket, name = "n_cases") %>%
  arrange(desc(n_cases), mechanism_bucket)

tmp_bucket_md <- tempfile(fileext = ".md")
write_md_table(bucket_counts, tmp_bucket_md)
md_lines <- c(md_lines, readLines(tmp_bucket_md), "", "## Case Summaries", "")

case_md_table <- casebook %>%
  transmute(
    case = case_label(Participant_id, from_tp, to_tp),
    mechanism_bucket,
    days_between,
    ST = paste0(coalesce(as.character(ST_from), "NA"), "->", coalesce(as.character(ST_to), "NA")),
    SNPs,
    same_strain_evidence,
    VF_delta = paste0("+", n_vf_genes_gained, "/-", n_vf_genes_lost),
    accessory_delta = paste0("+", n_accessory_genes_gained, "/-", n_accessory_genes_lost),
    plasmid_delta = paste0("+", n_plasmid_gained, "/-", n_plasmid_lost),
    AMR_delta = paste0("+", n_amr_gained, "/-", n_amr_lost),
    clinical_trigger
  )
tmp_case_md <- tempfile(fileext = ".md")
write_md_table(case_md_table, tmp_case_md)
md_lines <- c(md_lines, readLines(tmp_case_md), "", "## Interpretation Guardrails", "",
              "- Treat same-strain stable-profile cases as evidence against simple VF gain/loss explanations, not as evidence against regulation, expression, inoculum, or host-state effects.",
              "- Panaroo accessory-gene changes are shown as exploratory leads and are not used alone to reclassify stable VF/module cases, because low-SNP pairs can still show annotation/presence-absence noise.",
              "- Treat replacement cases separately from within-strain evolution.",
              "- AMR output is optional and should be read only if `amr_screen_report.txt` says COMPLETE.",
              "- UTI n remains sparse; these outputs are hypothesis-generating.")

writeLines(md_lines, casebook_md)

# ==============================================================================
# PLOTS
# ==============================================================================

plot_cases <- casebook %>%
  mutate(case = factor(case_label(Participant_id, from_tp, to_tp),
                       levels = rev(case_label(Participant_id, from_tp, to_tp)))) %>%
  transmute(
    case,
    `WGS/VF linked` = ifelse(has_vf_pair, "yes", "no"),
    `SNP strain context` = case_when(
      snp_strain_context == "Strong same strain" ~ "strong",
      snp_strain_context == "Above same-strain SNP threshold" ~ "above_threshold",
      snp_strain_context == "Missing SNP evidence" ~ "missing",
      TRUE ~ "uncertain"
    ),
    `ST lineage context` = case_when(
      st_lineage_context == "Same ST" ~ "same_ST",
      st_lineage_context == "Different ST" ~ "different_ST",
      TRUE ~ "missing_ST"
    ),
    `VF or module change` = ifelse((n_vf_genes_gained + n_vf_genes_lost + n_modules_gained + n_modules_lost) > 0, "yes", "no"),
    `Accessory change` = ifelse((n_accessory_genes_gained + n_accessory_genes_lost) > 0, "yes", "no"),
    `Plasmid change` = ifelse((n_plasmid_gained + n_plasmid_lost) > 0, "yes", "no"),
    `AMR change` = ifelse((n_amr_gained + n_amr_lost) > 0, "yes", "no"),
    `Host trigger` = str_trunc(clinical_trigger, 32)
  ) %>%
  pivot_longer(-case, names_to = "feature", values_to = "value") %>%
  mutate(feature = factor(feature, levels = unique(feature)))

p_case_matrix <- ggplot(plot_cases, aes(feature, case, fill = value)) +
  geom_tile(colour = "white", linewidth = 0.4) +
  scale_fill_manual(
    values = c(
      yes = "#2C7A7B", no = "#D8DEE9", strong = "#2F855A",
      above_threshold = "#D69E2E", replacement = "#C05621",
      same_ST = "#009E73", different_ST = "#D55E00", missing_ST = "#718096",
      missing = "#718096", uncertain = "#B7791F",
      `symptom emergence on culture-supported bacteriuria` = "#805AD5",
      `culture support plus symptoms emerged` = "#3182CE",
      `UTI criteria present at endpoint` = "#4A5568",
      `review clinical rule inputs` = "#E53E3E"
    ),
    na.value = "grey90"
  ) +
  labs(
    title = "Not_UTI -> UTI mechanism case matrix",
    x = NULL, y = NULL, fill = NULL,
    caption = "Primary UTI rule only; descriptive mechanism categories."
  ) +
  plot_theme_mechanism(9) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1))

ggsave(file.path(DIR_PLOTS_MECHANISM, "not_uti_to_uti_case_matrix.png"),
       p_case_matrix, width = 11, height = 6.5, dpi = 300)

plot_mech <- transition_mechanism_summary %>%
  mutate(
    transition_type = factor(transition_type, levels = c("Not_UTI->Not_UTI", "Not_UTI->UTI", "UTI->Not_UTI", "UTI->UTI"))
  )

p_strain <- ggplot(plot_mech, aes(transition_type, n_transitions, fill = mechanism_bucket)) +
  geom_col(position = "stack", colour = "white", linewidth = 0.2) +
  scale_fill_manual(values = c(
    same_strain_stable_profile = "#2F855A",
    same_strain_genomic_change = "#2B6CB0",
    strain_replacement = "#C05621",
    uncertain = "#B7791F",
    missing_wgs_endpoint = "#718096"
  )) +
  labs(
    title = "Strain stability and replacement across clinical transitions",
    x = NULL, y = "Transitions", fill = "Mechanism bucket",
    caption = "All ordered primary-status transitions from script 28."
  ) +
  plot_theme_mechanism(10)

ggsave(file.path(DIR_PLOTS_MECHANISM, "strain_replacement_vs_stability.png"),
       p_strain, width = 9, height = 5.5, dpi = 300)

symptom_cols <- c(
  to_dysuria_present = "Dysuria",
  to_urgency_present = "Urgency",
  to_frequency_present = "Frequency",
  to_incontinence_present = "Incontinence",
  to_pus_present = "Pus",
  to_flankpain_present = "Flank pain",
  to_suprapubic_pain_present = "Suprapubic pain",
  to_fever_present = "Fever",
  to_rigors_present = "Rigors",
  to_delirium_present = "Delirium",
  to_other_sxs_present = "Other"
)

plot_host <- casebook %>%
  mutate(case = factor(case_label(Participant_id, from_tp, to_tp),
                       levels = rev(case_label(Participant_id, from_tp, to_tp)))) %>%
  select(case, all_of(names(symptom_cols)), collection_method_changed, catheter_rule_changed) %>%
  rename(`Collection changed` = collection_method_changed,
         `Catheter rule changed` = catheter_rule_changed) %>%
  pivot_longer(-case, names_to = "feature", values_to = "present") %>%
  mutate(
    feature = recode(feature, !!!symptom_cols),
    present = bool_chr(present),
    feature = factor(feature, levels = c(unname(symptom_cols), "Collection changed", "Catheter rule changed"))
  )

p_host <- ggplot(plot_host, aes(feature, case, fill = present)) +
  geom_tile(colour = "white", linewidth = 0.4) +
  scale_fill_manual(values = c(yes = "#2C7A7B", no = "#E2E8F0", unknown = "#A0AEC0"), na.value = "#A0AEC0") +
  labs(
    title = "Host and symptom context at UTI endpoint",
    x = NULL, y = NULL, fill = NULL,
    caption = "Endpoint symptom flags and whether collection/catheter context changed from the preceding Not_UTI row."
  ) +
  plot_theme_mechanism(9) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1))

ggsave(file.path(DIR_PLOTS_MECHANISM, "host_context_transition_heatmap.png"),
       p_host, width = 11, height = 6.5, dpi = 300)

# ==============================================================================
# VALIDATION
# ==============================================================================

validation <- tibble(
  check = c(
    "casebook has exactly 11 clinical Not_UTI->UTI transitions",
    "casebook flags exactly 1 missing WGS/VF endpoint",
    "casebook has exactly 10 WGS/VF-linked Not_UTI->UTI transitions",
    "casebook has 0 Uricult-linked Not_UTI->UTI transitions",
    "vf_analysis_ready primary status has no missing UTI_Status",
    "legacy ASB-vs-UTI OLD files were not used as inputs"
  ),
  status = c(
    ifelse(nrow(casebook) == 11, "PASS", "FAIL"),
    ifelse(sum(!casebook$has_vf_pair, na.rm = TRUE) == 1, "PASS", "FAIL"),
    ifelse(sum(casebook$has_vf_pair %in% TRUE, na.rm = TRUE) == 10, "PASS", "FAIL"),
    ifelse(sum(casebook$is_uricult_transition %in% TRUE, na.rm = TRUE) == 0, "PASS", "FAIL"),
    ifelse(sum(is.na(vf_ready$UTI_Status)) == 0, "PASS", "FAIL"),
    "PASS"
  ),
  detail = c(
    sprintf("n=%d", nrow(casebook)),
    sprintf("missing=%d", sum(!casebook$has_vf_pair, na.rm = TRUE)),
    sprintf("linked=%d", sum(casebook$has_vf_pair %in% TRUE, na.rm = TRUE)),
    sprintf("uricult=%d", sum(casebook$is_uricult_transition %in% TRUE, na.rm = TRUE)),
    sprintf("missing_status=%d", sum(is.na(vf_ready$UTI_Status))),
    "Inputs restricted to current status_map, vf_transition, vf_ready, Panaroo, plasmid, variant, and optional ABRicate outputs."
  )
)

write_csv(validation, file.path(DIR_MECHANISM, "mechanism_validation_checks.csv"))
if (any(validation$status != "PASS")) {
  print(validation)
  stop("Mechanism add-on validation failed.")
}

msg("Mechanism add-on complete.")
msg("Outputs written to %s and %s", DIR_MECHANISM, DIR_PLOTS_MECHANISM)

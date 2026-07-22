#!/usr/bin/env Rscript

# Genetics-first candidate-cluster screen for nursing-home linkage.
#
# This script:
#   1. Reconciles the emailed 556-row VF file with the current 532-isolate
#      QC-passing Longcycler cohort.
#   2. Screens cross-participant isolate pairs using ST, 227-gene VF Jaccard
#      similarity, and the all-isolate core-SNP matrix.
#   3. Directly validates exact-VF/core-SNP<=25 pairs through the existing
#      provenance-aware 11_compare_strains.R workflow.
#   4. Adds validated script-29 primary genomic-AMR profiles as context only.
#   5. Writes flat CSV/JSON outputs consumed by the final workbook and DOCX
#      builders.

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
  library(tidyr)
  library(purrr)
  library(tibble)
})

source("00_config.R")
source("R/pipeline_qc_helpers.R")
source("11_compare_strains_helpers.R")

parse_args <- function(args) {
  out <- list(
    output_dir = file.path("outputs", "nursing_home_candidate_clusters_20260717"),
    external_vf = NA_character_,
    stage = "all",
    workers = 8L
  )
  for (arg in args) {
    if (!startsWith(arg, "--") || !grepl("=", arg, fixed = TRUE)) next
    bits <- strsplit(sub("^--", "", arg), "=", fixed = TRUE)[[1]]
    key <- bits[[1]]
    value <- paste(bits[-1], collapse = "=")
    if (key %in% names(out)) out[[key]] <- value
  }
  out$workers <- max(1L, as.integer(out$workers))
  out$stage <- tolower(out$stage)
  out
}

opt <- parse_args(commandArgs(trailingOnly = TRUE))
if (!opt$stage %in% c("prepare", "finalize", "all")) {
  stop("--stage must be prepare, finalize, or all")
}
if (is.na(opt$external_vf) || !file.exists(opt$external_vf)) {
  stop("Provide an existing --external_vf=/path/to/vf_analysis_ready.csv")
}

dir.create(opt$output_dir, recursive = TRUE, showWarnings = FALSE)
direct_dir <- file.path(opt$output_dir, "direct_pairwise")
dir.create(direct_dir, recursive = TRUE, showWarnings = FALSE)

path_preliminary <- file.path(opt$output_dir, "preliminary_candidate_pairs.csv")
path_direct_input <- file.path(opt$output_dir, "direct_pairs_input.csv")
path_direct_metrics <- file.path(opt$output_dir, "direct_pairwise_exact.csv")
path_candidate_pairs <- file.path(opt$output_dir, "candidate_pairs.csv")
path_cluster_summary <- file.path(opt$output_dir, "cluster_summary.csv")
path_isolate_profiles <- file.path(opt$output_dir, "isolate_profiles.csv")
path_excluded <- file.path(opt$output_dir, "excluded_isolates.csv")
path_guide <- file.path(opt$output_dir, "interpretation_guide.csv")
path_summary_json <- file.path(opt$output_dir, "analysis_summary.json")

key2 <- function(pid, tp) paste(as.character(pid), as.character(tp), sep = "__")

normalise_bool <- function(x) {
  if (is.logical(x)) return(replace(x, is.na(x), FALSE))
  tolower(trimws(as.character(x))) %in% c("true", "t", "1", "yes", "y")
}

parse_collection_date <- function(x) {
  x <- as.character(x)
  d <- suppressWarnings(as.Date(x, format = "%d/%m/%Y"))
  missing <- is.na(d)
  d[missing] <- suppressWarnings(as.Date(x[missing], format = "%Y-%m-%d"))
  d
}

clean_text <- function(x) {
  x <- as.character(x)
  x[is.na(x)] <- ""
  x
}

collapse_sorted <- function(x, sep = "; ") {
  x <- sort(unique(clean_text(x)))
  x <- x[nzchar(x)]
  if (!length(x)) "" else paste(x, collapse = sep)
}

set_jaccard <- function(a, b) {
  a <- sort(unique(a[nzchar(a)]))
  b <- sort(unique(b[nzchar(b)]))
  un <- union(a, b)
  if (!length(un)) return(NA_real_)
  length(intersect(a, b)) / length(un)
}

safe_min <- function(x) if (all(is.na(x))) NA_real_ else min(x, na.rm = TRUE)
safe_max <- function(x) if (all(is.na(x))) NA_real_ else max(x, na.rm = TRUE)

disjoint_components <- function(nodes, edges_a, edges_b, prefix) {
  nodes <- sort(unique(as.character(nodes)))
  if (!length(nodes)) return(tibble(node = character(), component_id = character()))
  parent <- setNames(nodes, nodes)
  find_root <- function(x) {
    while (parent[[x]] != x) {
      parent[[x]] <<- parent[[parent[[x]]]]
      x <- parent[[x]]
    }
    x
  }
  union_nodes <- function(a, b) {
    ra <- find_root(a)
    rb <- find_root(b)
    if (ra != rb) parent[[rb]] <<- ra
  }
  if (length(edges_a)) {
    for (i in seq_along(edges_a)) union_nodes(as.character(edges_a[[i]]), as.character(edges_b[[i]]))
  }
  roots <- vapply(nodes, find_root, character(1))
  root_levels <- sort(unique(roots))
  tibble(
    node = nodes,
    component_id = sprintf("%s-%03d", prefix, match(roots, root_levels))
  )
}

read_core_inputs <- function() {
  vf <- readr::read_csv(FILE_VF_READY, show_col_types = FALSE) %>%
    mutate(
      Participant_id = as.character(Participant_id),
      tp_lab = as.character(tp_lab),
      sample_key = key2(Participant_id, tp_lab),
      Collection_Date_parsed = parse_collection_date(Collection_Date),
      ST = as.character(ST),
      ST_reliable = ST_source == "provider_qc95" &
        !is.na(provider_PercGoodTargets) & provider_PercGoodTargets >= 95
    )
  gene_cols <- canonical_vf_gene_cols(names(vf), vf_pa_file = FILE_VF_PA)
  if (length(gene_cols) != 227L) {
    stop("Expected 227 canonical VF gene columns; found ", length(gene_cols))
  }
  if (nrow(vf) != 532L || n_distinct(vf$sample_key) != 532L) {
    stop("Current VF-ready cohort must contain 532 unique participant/timepoint rows.")
  }

  asm <- readr::read_csv(FILE_ANALYSIS_ASSEMBLY_MANIFEST, show_col_types = FALSE) %>%
    mutate(
      Participant_id = as.character(Participant_id),
      tp_lab = as.character(tp_lab),
      sample_key = key2(Participant_id, tp_lab)
    ) %>%
    filter(normalise_bool(selected_canonical), normalise_bool(QC_PASS)) %>%
    select(
      sample_key, Isolate_ID, Assembly_ID, assembler, full_path, fasta_sha256,
      QC_PASS, QC_REASON, selected_canonical, canonical_reason
    )
  if (nrow(asm) != 532L || n_distinct(asm$sample_key) != 532L) {
    stop("Analysis assembly manifest must contain 532 unique selected QC-passing rows.")
  }

  sample_map <- readr::read_csv(
    file.path(DIR_WGS_CORE, "core_snp_sample_map.csv"),
    show_col_types = FALSE
  ) %>%
    mutate(
      Participant_id = as.character(Participant_id),
      tp_lab = as.character(tp_lab),
      sample_key = key2(Participant_id, tp_lab)
    ) %>%
    select(sample_key, parsnp_alignment_label)
  if (nrow(sample_map) != 532L || n_distinct(sample_map$sample_key) != 532L) {
    stop("Core-SNP sample map must contain 532 unique participant/timepoint rows.")
  }

  meta <- vf %>%
    left_join(asm, by = "sample_key") %>%
    left_join(sample_map, by = "sample_key")
  if (anyNA(meta$Isolate_ID) || anyNA(meta$parsnp_alignment_label)) {
    stop("Failed to map every VF-ready row to the selected assembly and core-SNP label.")
  }

  list(vf = vf, meta = meta, gene_cols = gene_cols)
}

build_preliminary_pairs <- function(inputs) {
  meta <- inputs$meta
  gene_cols <- inputs$gene_cols
  core_pairs <- readr::read_csv(
    file.path(DIR_WGS_CORE, "strain_pairs.csv"),
    show_col_types = FALSE
  )
  idx_a <- match(core_pairs$A, meta$parsnp_alignment_label)
  idx_b <- match(core_pairs$B, meta$parsnp_alignment_label)
  if (anyNA(idx_a) || anyNA(idx_b)) stop("Unmapped core-SNP pair label(s).")

  vf_matrix <- as.matrix(meta[, gene_cols, drop = FALSE])
  storage.mode(vf_matrix) <- "integer"
  intersection_n <- rowSums(
    vf_matrix[idx_a, , drop = FALSE] * vf_matrix[idx_b, , drop = FALSE],
    na.rm = TRUE
  )
  union_n <- rowSums(
    (vf_matrix[idx_a, , drop = FALSE] + vf_matrix[idx_b, , drop = FALSE]) > 0,
    na.rm = TRUE
  )
  vf_jaccard <- ifelse(union_n > 0, intersection_n / union_n, NA_real_)
  same_st <- !is.na(meta$ST[idx_a]) & !is.na(meta$ST[idx_b]) &
    nzchar(meta$ST[idx_a]) & nzchar(meta$ST[idx_b]) &
    meta$ST[idx_a] == meta$ST[idx_b]
  cross_participant <- meta$Participant_id[idx_a] != meta$Participant_id[idx_b]

  pair_tbl <- tibble(
    pair_id = paste(meta$sample_key[idx_a], meta$sample_key[idx_b], sep = "__vs__"),
    sample_key_A = meta$sample_key[idx_a],
    sample_key_B = meta$sample_key[idx_b],
    Participant_id_A = meta$Participant_id[idx_a],
    Timepoint_A = meta$tp_lab[idx_a],
    Participant_id_B = meta$Participant_id[idx_b],
    Timepoint_B = meta$tp_lab[idx_b],
    preliminary_core_snps = as.numeric(core_pairs$snps),
    preliminary_core_call = as.character(core_pairs$call),
    ST_A = meta$ST[idx_a],
    ST_B = meta$ST[idx_b],
    ST_source_A = meta$ST_source[idx_a],
    ST_source_B = meta$ST_source[idx_b],
    ST_reliable_A = meta$ST_reliable[idx_a],
    ST_reliable_B = meta$ST_reliable[idx_b],
    provider_PercGoodTargets_A = meta$provider_PercGoodTargets[idx_a],
    provider_PercGoodTargets_B = meta$provider_PercGoodTargets[idx_b],
    vf_jaccard = vf_jaccard,
    vf_intersection_n = intersection_n,
    vf_union_n = union_n,
    vf_identical = !is.na(vf_jaccard) & abs(vf_jaccard - 1) < 1e-12,
    meta_index_A = idx_a,
    meta_index_B = idx_b
  ) %>%
    filter(
      cross_participant,
      same_st,
      !is.na(vf_jaccard),
      vf_jaccard >= 0.90,
      preliminary_core_snps <= 25
    ) %>%
    arrange(preliminary_core_snps, desc(vf_jaccard), ST_A, pair_id)

  if (!nrow(pair_tbl)) stop("No preliminary candidate pairs met the screen.")
  pair_tbl
}

prepare_stage <- function(inputs) {
  pair_tbl <- build_preliminary_pairs(inputs)
  write_csv(
    pair_tbl %>% select(-meta_index_A, -meta_index_B),
    path_preliminary,
    na = ""
  )
  direct_input <- pair_tbl %>%
    filter(vf_identical) %>%
    transmute(
      Participant_id_A, Timepoint_A,
      Participant_id_B, Timepoint_B
    ) %>%
    distinct()
  write_csv(direct_input, path_direct_input, na = "")
  message(
    "Prepared ", nrow(pair_tbl), " preliminary pairs and ",
    nrow(direct_input), " exact-VF pairs for direct validation."
  )
  invisible(pair_tbl)
}

run_direct_validation <- function() {
  if (!file.exists(path_direct_input)) stop("Missing ", path_direct_input)
  direct_input <- readr::read_csv(path_direct_input, show_col_types = FALSE) %>%
    mutate(
      Participant_id_A = as.character(Participant_id_A),
      Timepoint_A = as.character(Timepoint_A),
      Participant_id_B = as.character(Participant_id_B),
      Timepoint_B = as.character(Timepoint_B)
    )
  asm <- readr::read_csv(FILE_ANALYSIS_ASSEMBLY_MANIFEST, show_col_types = FALSE) %>%
    mutate(
      Participant_id = as.character(Participant_id),
      tp_lab = as.character(tp_lab),
      full_path = normalizePath(full_path, winslash = "/", mustWork = TRUE),
      selected_canonical = normalise_bool(selected_canonical),
      QC_PASS = normalise_bool(QC_PASS)
    ) %>%
    filter(selected_canonical, QC_PASS)
  file_info <- file.info(asm$full_path)
  asm <- asm %>%
    mutate(
      fasta_size_bytes = unname(as.numeric(file_info$size)),
      fasta_mtime_utc = format(
        file_info$mtime,
        "%Y-%m-%dT%H:%M:%OS6Z",
        tz = "UTC"
      )
    )
  if (nrow(asm) != 532L || anyNA(asm$fasta_sha256)) {
    stop("Direct validation requires 532 selected assemblies with SHA-256 fingerprints.")
  }

  cache_dir <- file.path(direct_dir, "dnadiff_cache_sha256_v1")
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  resolve_exact <- function(pid, tp) {
    row <- asm %>%
      filter(Participant_id == as.character(pid), tp_lab == as.character(tp))
    if (nrow(row) != 1L) {
      stop("Expected one exact assembly for ", pid, "__", tp, "; found ", nrow(row))
    }
    row
  }
  compute_one <- function(i) {
    row <- direct_input[i, , drop = FALSE]
    a <- resolve_exact(row$Participant_id_A[[1]], row$Timepoint_A[[1]])
    b <- resolve_exact(row$Participant_id_B[[1]], row$Timepoint_B[[1]])
    key <- paste0(
      key2(row$Participant_id_A[[1]], row$Timepoint_A[[1]]),
      "__vs__",
      key2(row$Participant_id_B[[1]], row$Timepoint_B[[1]])
    )
    fp_a <- list(
      path = a$full_path[[1]],
      sha256 = a$fasta_sha256[[1]],
      size_bytes = a$fasta_size_bytes[[1]],
      mtime_utc = a$fasta_mtime_utc[[1]]
    )
    fp_b <- list(
      path = b$full_path[[1]],
      sha256 = b$fasta_sha256[[1]],
      size_bytes = b$fasta_size_bytes[[1]],
      mtime_utc = b$fasta_mtime_utc[[1]]
    )
    result <- run_dnadiff(
      a$full_path[[1]], b$full_path[[1]], cache_dir, key,
      a_fingerprint = fp_a, b_fingerprint = fp_b
    )
    bind_cols(row, result)
  }
  workers <- min(opt$workers, nrow(direct_input))
  message(
    "Running exact-label direct validation for ", nrow(direct_input),
    " pairs with ", workers, " workers."
  )
  pair_list <- if (.Platform$OS.type != "windows" && workers > 1L) {
    parallel::mclapply(
      seq_len(nrow(direct_input)),
      compute_one,
      mc.cores = workers,
      mc.preschedule = TRUE
    )
  } else {
    lapply(seq_len(nrow(direct_input)), compute_one)
  }
  exact_direct <- bind_rows(pair_list)
  if (nrow(exact_direct) != nrow(direct_input)) {
    stop("Exact-label direct comparison did not return one row per requested pair.")
  }
  write_csv(exact_direct, path_direct_metrics, na = "")
  if (any(is.na(exact_direct$TotalSNPs))) {
    warning(sum(is.na(exact_direct$TotalSNPs)), " direct comparison(s) lack SNP results.")
  }
}

normalise_amr_gene <- function(x) {
  str_replace(as.character(x), "_[0-9]+$", "")
}

resfinder_broad_class <- function(drug_name) {
  drug <- str_to_lower(clean_text(drug_name))
  case_when(
    str_detect(
      drug,
      "amoxicillin|ampicillin|aztreonam|cefepime|cefotaxime|cefoxitin|ceftazidime|ceftriaxone|cephalothin|piperacillin|ticarcillin"
    ) ~ "Beta-lactams",
    str_detect(drug, "gentamicin|streptomycin|tobramycin") ~ "Aminoglycosides",
    str_detect(drug, "doxycycline|minocycline|tetracycline") ~ "Tetracyclines",
    str_detect(drug, "chloramphenicol|florfenicol") ~ "Phenicols",
    str_detect(drug, "sulfamethoxazole") ~ "Sulfonamides",
    str_detect(drug, "trimethoprim") ~ "Trimethoprim",
    str_detect(drug, "ciprofloxacin") ~ "Fluoroquinolones",
    str_detect(drug, "azithromycin|erythromycin|spiramycin|telithromycin") ~
      "Macrolides/ketolides",
    str_detect(drug, "lincomycin") ~ "Lincosamides",
    str_detect(drug, "fosfomycin") ~ "Fosfomycin",
    nzchar(drug) ~ "Other/unspecified",
    TRUE ~ ""
  )
}

build_amr_profiles <- function(meta) {
  profile_path <- file.path(DIR_RESULTS, "amr", "episode_amr_profiles.csv")
  if (!file.exists(profile_path)) {
    stop("Validated script-29 AMR profiles are missing: ", profile_path)
  }
  episode <- readr::read_csv(profile_path, show_col_types = FALSE) %>%
    mutate(
      Participant_id = as.character(Participant_id),
      tp_lab = as.character(tp_lab),
      sample_key = key2(Participant_id, tp_lab),
      amr_genes = strsplit(
        replace_na(as.character(informative_acquired_genes), ""),
        ";", fixed = TRUE
      ),
      amr_classes = strsplit(
        replace_na(as.character(informative_acquired_classes), ""),
        ";", fixed = TRUE
      ),
      amr_genes = map(amr_genes, ~ sort(unique(.x[nzchar(.x)]))),
      amr_classes = map(amr_classes, ~ sort(unique(.x[nzchar(.x)])))
    )
  if (nrow(episode) != 532L || anyDuplicated(episode$sample_key)) {
    stop("Script-29 AMR profile must contain 532 unique episode keys.")
  }
  gene_sets <- setNames(episode$amr_genes, episode$sample_key)
  class_sets <- setNames(episode$amr_classes, episode$sample_key)
  profiles <- meta %>%
    transmute(sample_key) %>%
    mutate(
      amr_genes_all = map(sample_key, ~ sort(unique(gene_sets[[.x]] %||% character()))),
      amr_genes_excluding_mdfA = amr_genes_all,
      mdfA_detected = map_lgl(
        sample_key,
        ~ episode$mdfA_detected[match(.x, episode$sample_key)] %in% TRUE
      ),
      amr_classes = map(sample_key, ~ sort(unique(class_sets[[.x]] %||% character())))
    )
  profiles
}

`%||%` <- function(x, y) if (is.null(x)) y else x

vf_diff_for_pair <- function(meta, gene_cols, idx_a, idx_b) {
  va <- as.integer(unlist(meta[idx_a, gene_cols, drop = FALSE], use.names = FALSE))
  vb <- as.integer(unlist(meta[idx_b, gene_cols, drop = FALSE], use.names = FALSE))
  names(va) <- gene_cols
  names(vb) <- gene_cols
  list(
    only_a = names(va)[va > 0 & vb <= 0],
    only_b = names(vb)[vb > 0 & va <= 0]
  )
}

add_pair_profiles <- function(pair_tbl, inputs, amr_profiles) {
  meta <- inputs$meta
  gene_cols <- inputs$gene_cols
  idx_a <- match(pair_tbl$sample_key_A, meta$sample_key)
  idx_b <- match(pair_tbl$sample_key_B, meta$sample_key)
  amr_idx_a <- match(pair_tbl$sample_key_A, amr_profiles$sample_key)
  amr_idx_b <- match(pair_tbl$sample_key_B, amr_profiles$sample_key)

  vf_diff <- map2(idx_a, idx_b, ~ vf_diff_for_pair(meta, gene_cols, .x, .y))
  amr_a <- amr_profiles$amr_genes_excluding_mdfA[amr_idx_a]
  amr_b <- amr_profiles$amr_genes_excluding_mdfA[amr_idx_b]
  amr_all_a <- amr_profiles$amr_genes_all[amr_idx_a]
  amr_all_b <- amr_profiles$amr_genes_all[amr_idx_b]
  class_a <- amr_profiles$amr_classes[amr_idx_a]
  class_b <- amr_profiles$amr_classes[amr_idx_b]

  pair_tbl %>%
    mutate(
      Episode_ID_A = meta$Episode_ID[idx_a],
      Episode_ID_B = meta$Episode_ID[idx_b],
      Isolate_ID_A = meta$Isolate_ID[idx_a],
      Isolate_ID_B = meta$Isolate_ID[idx_b],
      Collection_Date_A = meta$Collection_Date_parsed[idx_a],
      Collection_Date_B = meta$Collection_Date_parsed[idx_b],
      collection_date_gap_days = abs(as.integer(Collection_Date_A - Collection_Date_B)),
      Event_type_A = meta$Event_type[idx_a],
      Event_type_B = meta$Event_type[idx_b],
      UTI_Status_A = meta$UTI_Status[idx_a],
      UTI_Status_B = meta$UTI_Status[idx_b],
      Primary_Status_A = meta$Primary_Status[idx_a],
      Primary_Status_B = meta$Primary_Status[idx_b],
      vf_discordant_gene_n = map_int(vf_diff, ~ length(.x$only_a) + length(.x$only_b)),
      vf_genes_only_A = map_chr(vf_diff, ~ collapse_sorted(.x$only_a)),
      vf_genes_only_B = map_chr(vf_diff, ~ collapse_sorted(.x$only_b)),
      amr_jaccard_excluding_mdfA = map2_dbl(amr_a, amr_b, set_jaccard),
      shared_amr_genes_excluding_mdfA = map2_chr(
        amr_a, amr_b, ~ collapse_sorted(intersect(.x, .y))
      ),
      amr_genes_only_A_excluding_mdfA = map2_chr(
        amr_a, amr_b, ~ collapse_sorted(setdiff(.x, .y))
      ),
      amr_genes_only_B_excluding_mdfA = map2_chr(
        amr_a, amr_b, ~ collapse_sorted(setdiff(.y, .x))
      ),
      all_amr_genes_A = map_chr(amr_all_a, collapse_sorted),
      all_amr_genes_B = map_chr(amr_all_b, collapse_sorted),
      mdfA_detected_A = amr_profiles$mdfA_detected[amr_idx_a],
      mdfA_detected_B = amr_profiles$mdfA_detected[amr_idx_b],
      shared_amr_classes = map2_chr(
        class_a, class_b, ~ collapse_sorted(intersect(.x, .y))
      ),
      all_amr_classes_A = map_chr(class_a, collapse_sorted),
      all_amr_classes_B = map_chr(class_b, collapse_sorted)
    )
}

assign_groups <- function(pair_tbl) {
  eligible_screen <- pair_tbl$priority_category != "Same lineage, not the same strain"
  screen_edges <- pair_tbl[eligible_screen, , drop = FALSE]
  screen_nodes <- unique(c(screen_edges$sample_key_A, screen_edges$sample_key_B))
  screen_components <- disjoint_components(
    screen_nodes,
    screen_edges$sample_key_A,
    screen_edges$sample_key_B,
    "SCREEN"
  )

  high_edges <- pair_tbl %>% filter(priority_category == "High-priority candidate")
  high_nodes <- unique(c(high_edges$sample_key_A, high_edges$sample_key_B))
  high_components <- disjoint_components(
    high_nodes,
    high_edges$sample_key_A,
    high_edges$sample_key_B,
    "CGC"
  )

  pair_tbl %>%
    mutate(
      screen_group_A = screen_components$component_id[
        match(sample_key_A, screen_components$node)
      ],
      screen_group_B = screen_components$component_id[
        match(sample_key_B, screen_components$node)
      ],
      screen_group_id = if_else(
        !is.na(screen_group_A) & screen_group_A == screen_group_B,
        screen_group_A,
        NA_character_
      ),
      candidate_cluster_A = high_components$component_id[
        match(sample_key_A, high_components$node)
      ],
      candidate_cluster_B = high_components$component_id[
        match(sample_key_B, high_components$node)
      ],
      candidate_cluster_id = if_else(
        !is.na(candidate_cluster_A) & candidate_cluster_A == candidate_cluster_B,
        candidate_cluster_A,
        NA_character_
      )
    ) %>%
    select(-screen_group_A, -screen_group_B,
           -candidate_cluster_A, -candidate_cluster_B)
}

build_cluster_summary <- function(pair_tbl, inputs, amr_profiles) {
  meta <- inputs$meta
  high <- pair_tbl %>% filter(priority_category == "High-priority candidate")
  cluster_ids <- sort(unique(na.omit(high$candidate_cluster_id)))
  map_dfr(cluster_ids, function(cid) {
    edges <- high %>% filter(candidate_cluster_id == cid)
    members <- sort(unique(c(edges$sample_key_A, edges$sample_key_B)))
    midx <- match(members, meta$sample_key)
    aidx <- match(members, amr_profiles$sample_key)
    amr_sets <- amr_profiles$amr_genes_excluding_mdfA[aidx]
    class_sets <- amr_profiles$amr_classes[aidx]
    shared_amr <- if (length(amr_sets)) Reduce(intersect, amr_sets) else character()
    union_amr <- if (length(amr_sets)) Reduce(union, amr_sets) else character()
    shared_classes <- if (length(class_sets)) Reduce(intersect, class_sets) else character()
    union_classes <- if (length(class_sets)) Reduce(union, class_sets) else character()
    dates <- meta$Collection_Date_parsed[midx]
    participant_counts <- table(meta$Participant_id[midx])
    n_all_pairs <- length(members) * (length(members) - 1) / 2
    n_within_participant_pairs <- sum(
      participant_counts * (participant_counts - 1) / 2
    )
    n_possible_cross_participant_pairs <- n_all_pairs - n_within_participant_pairs
    tibble(
      candidate_cluster_id = cid,
      interpretation = "Candidate genomic cluster; requires nursing-home, ward and temporal confirmation.",
      n_isolates = length(members),
      n_participants = n_distinct(meta$Participant_id[midx]),
      Participant_ids = collapse_sorted(meta$Participant_id[midx]),
      Episode_IDs = collapse_sorted(meta$Episode_ID[midx]),
      Isolate_IDs = collapse_sorted(meta$Isolate_ID[midx]),
      ST = collapse_sorted(meta$ST[midx]),
      earliest_collection_date = if (all(is.na(dates))) as.Date(NA) else min(dates, na.rm = TRUE),
      latest_collection_date = if (all(is.na(dates))) as.Date(NA) else max(dates, na.rm = TRUE),
      collection_span_days = if (all(is.na(dates))) NA_integer_ else
        as.integer(max(dates, na.rm = TRUE) - min(dates, na.rm = TRUE)),
      high_priority_pair_n = nrow(edges),
      possible_cross_participant_pair_n = n_possible_cross_participant_pairs,
      chain_linked_flag = nrow(edges) < n_possible_cross_participant_pairs,
      direct_snp_min = safe_min(edges$direct_total_snps),
      direct_snp_max = safe_max(edges$direct_total_snps),
      preliminary_core_snp_min = safe_min(edges$preliminary_core_snps),
      preliminary_core_snp_max = safe_max(edges$preliminary_core_snps),
      vf_jaccard_min = safe_min(edges$vf_jaccard),
      vf_jaccard_max = safe_max(edges$vf_jaccard),
      shared_amr_genes_all_members_excluding_mdfA = collapse_sorted(shared_amr),
      variable_amr_genes_excluding_mdfA = collapse_sorted(setdiff(union_amr, shared_amr)),
      mdfA_detected_in_all_members = all(amr_profiles$mdfA_detected[aidx]),
      shared_amr_classes_all_members = collapse_sorted(shared_classes),
      resistance_classes_present = collapse_sorted(union_classes),
      reporting_note = ifelse(
        nrow(edges) < n_possible_cross_participant_pairs,
        "Chain-linked component: review all pair rows before interpreting the group.",
        "All cross-participant isolate pairs in this component meet the high-priority rule."
      )
    )
  })
}

build_isolate_profiles <- function(pair_tbl, inputs, amr_profiles) {
  meta <- inputs$meta
  high_nodes <- pair_tbl %>%
    filter(priority_category == "High-priority candidate", !is.na(candidate_cluster_id)) %>%
    select(candidate_cluster_id, sample_key_A, sample_key_B) %>%
    pivot_longer(c(sample_key_A, sample_key_B), values_to = "sample_key") %>%
    select(candidate_cluster_id, sample_key) %>%
    distinct()
  screen_nodes <- pair_tbl %>%
    filter(priority_category != "Same lineage, not the same strain",
           !is.na(screen_group_id)) %>%
    select(screen_group_id, sample_key_A, sample_key_B) %>%
    pivot_longer(c(sample_key_A, sample_key_B), values_to = "sample_key") %>%
    select(screen_group_id, sample_key) %>%
    distinct()

  meta %>%
    left_join(amr_profiles, by = "sample_key") %>%
    left_join(high_nodes, by = "sample_key") %>%
    left_join(screen_nodes, by = "sample_key") %>%
    transmute(
      candidate_cluster_id,
      screen_group_id,
      Participant_id,
      tp_lab,
      Episode_ID,
      Isolate_ID,
      Collection_Date = Collection_Date_parsed,
      Event_type,
      UTI_Status,
      Primary_Status,
      QC_PASS = as.logical(QC_PASS),
      QC_REASON,
      ST,
      ST_source,
      provider_PercGoodTargets,
      ST_reliable,
      vf_gene_count = rowSums(across(all_of(inputs$gene_cols)) > 0, na.rm = TRUE),
      amr_gene_count_excluding_mdfA = map_int(amr_genes_excluding_mdfA, length),
      amr_genes_excluding_mdfA = map_chr(amr_genes_excluding_mdfA, collapse_sorted),
      mdfA_detected,
      all_amr_genes = map_chr(amr_genes_all, collapse_sorted),
      resistance_classes = map_chr(amr_classes, collapse_sorted),
      phenotype_available = FALSE,
      phenotype_note = "No phenotypic susceptibility or antibiogram data were available."
    ) %>%
    arrange(candidate_cluster_id, screen_group_id, ST, Participant_id, Collection_Date)
}

build_excluded_isolates <- function(inputs) {
  external <- readr::read_csv(opt$external_vf, show_col_types = FALSE) %>%
    mutate(
      Participant_id = as.character(Participant_id),
      tp_lab = as.character(tp_lab),
      sample_key = key2(Participant_id, tp_lab)
    )
  current_keys <- inputs$meta$sample_key
  excluded <- external %>% filter(!sample_key %in% current_keys)
  if (nrow(excluded) != 24L) {
    stop("Expected 24 emailed-file rows outside the current cohort; found ", nrow(excluded))
  }

  canonical <- readr::read_csv(
    file.path(DIR_QC, "canonical_assembly_selection.csv"),
    show_col_types = FALSE
  ) %>%
    mutate(
      Participant_id = as.character(Participant_id),
      tp_lab = as.character(tp_lab),
      sample_key = key2(Participant_id, tp_lab)
    ) %>%
    filter(sample_key %in% excluded$sample_key) %>%
    group_by(sample_key) %>%
    summarise(
      current_manifest_rows = n(),
      any_selected_canonical = any(normalise_bool(selected_canonical)),
      any_QC_PASS = any(normalise_bool(QC_PASS)),
      QC_REASON_current = collapse_sorted(QC_REASON),
      canonical_reason_current = collapse_sorted(canonical_reason),
      .groups = "drop"
    )

  excluded %>%
    left_join(canonical, by = "sample_key") %>%
    transmute(
      Participant_id,
      tp_lab,
      Episode_ID,
      Collection_Date = parse_collection_date(Collection_Date),
      Event_type,
      ST,
      ST_source,
      provider_PercGoodTargets,
      emailed_analysis_include_primary = analysis_include_primary,
      emailed_genomics_expected_include = genomics_expected_include,
      current_manifest_rows = replace_na(current_manifest_rows, 0L),
      any_selected_canonical = replace_na(any_selected_canonical, FALSE),
      any_QC_PASS = replace_na(any_QC_PASS, FALSE),
      current_QC_reason = if_else(
        nzchar(clean_text(QC_REASON_current)),
        QC_REASON_current,
        "No matching row in the current canonical assembly manifest"
      ),
      current_selection_reason = if_else(
        nzchar(clean_text(canonical_reason_current)),
        canonical_reason_current,
        "No selected QC-passing assembly in the current cohort"
      ),
      analysis_decision = "Excluded from the 532-isolate cluster screen"
    ) %>%
    arrange(Participant_id, tp_lab)
}

build_interpretation_guide <- function() {
  tribble(
    ~section, ~term, ~plain_language_meaning, ~recommended_wording, ~avoid_wording,
    "Purpose", "Candidate genomic cluster",
    "A group of isolates with strong genomic similarity that should be checked against nursing-home, ward and timing information.",
    "Candidate genomic cluster; requires epidemiological confirmation.",
    "Confirmed outbreak or proven transmission.",
    "Evidence", "ST",
    "Sequence type describes a broad bacterial lineage. The same ST is useful for screening but is not strain-level proof.",
    "The isolates share ST [number], supporting membership of the same lineage.",
    "They are the same strain because they have the same ST.",
    "Evidence", "VF Jaccard similarity",
    "The proportion of detected virulence-factor genes shared by two isolates. 1.00 means identical detected profiles across the 227 genes.",
    "The isolates had identical/highly similar detected virulence-factor profiles.",
    "The virulence profile proves transmission.",
    "Evidence", "Preliminary core-SNP distance",
    "A study-wide screening distance from the common core-genome alignment. It is used to shortlist pairs and is not the final pair-specific result.",
    "The pair passed the preliminary core-SNP screen.",
    "This screening distance alone confirms the same strain.",
    "Evidence", "Direct SNP distance",
    "Pair-specific whole-assembly comparison. The project uses <=25 SNPs as an operational same-strain threshold.",
    "The pair was genetically closely related under the project's operational <=25-SNP rule.",
    "The <=25-SNP result proves who transmitted to whom.",
    "Resistance", "Primary genomic AMR gene",
    "An acquired resistance-associated gene in the script-29 AMRFinderPlus primary profile; mdf(A) is excluded from informative burden and similarity.",
    "The genome contains [gene], which is associated with [class]; phenotypic susceptibility was not available for confirmation.",
    "The isolate is clinically resistant to [drug].",
    "Resistance", "AMR-profile similarity",
    "Jaccard similarity of informative AMRFinderPlus acquired genes; the near-ubiquitous mdf(A) marker is excluded.",
    "The isolates shared a similar detected resistance-gene profile.",
    "The shared AMR profile defines the cluster.",
    "Limitations", "Nursing-home linkage",
    "Facility and ward variables are absent from these genomic files. The colleague must link Participant_id to her local metadata.",
    "Please check whether cluster members overlap by nursing home, ward and time.",
    "This is a nursing-home cluster.",
    "Limitations", "Phenotype",
    "No MIC, disk-diffusion or susceptible/intermediate/resistant results were found.",
    "These are genotypic AMR findings only.",
    "Genotypic detection is confirmed clinical resistance.",
    "Limitations", "Plasmids",
    "Plasmid replicon screening is not included in this resistance-gene-only report and a shared replicon would not prove an identical plasmid.",
    "Plasmid structure was not assessed in this report.",
    "The isolates share the same plasmid.",
    "Methods", "High-priority rule",
    "Different participants, same reliable provider ST, identical 227-gene VF profile, preliminary core SNPs <=25, and direct SNPs <=25.",
    "High-priority candidate genomic relationship.",
    "Confirmed transmission.",
    "Methods", "Moderate-priority rule",
    "Different participants, same reliable provider ST, VF Jaccard >=0.95, and preliminary core SNPs <=25; direct validation is pending or supporting information is incomplete.",
    "Moderate-priority candidate for epidemiological linkage.",
    "Same strain.",
    "Methods", "Possible related lineage",
    "Different participants with the same ST, VF Jaccard 0.90-0.949, and preliminary core SNPs <=25, or incomplete ST-quality support.",
    "Possible related lineage; check facility, ward and timing.",
    "Transmission cluster.",
    "Source", "Active VF dataset",
    FILE_VF_READY,
    "Current 532-isolate QC-selected VF dataset.",
    "",
    "Source", "Direct pair method",
    "11_compare_strains.R using dnadiff and the project's operational 25-SNP threshold.",
    "Pair-specific validation method.",
    "",
    "Source", "AMR data",
    file.path(DIR_RESULTS, "amr", "episode_amr_profiles.csv"),
    "Validated AMRFinderPlus primary profiles; ResFinder/PointFinder is complementary evidence.",
    ""
  )
}

finalize_stage <- function(inputs) {
  if (!file.exists(path_preliminary)) stop("Missing ", path_preliminary)
  if (!file.exists(path_direct_metrics)) stop("Missing ", path_direct_metrics)

  preliminary <- readr::read_csv(path_preliminary, show_col_types = FALSE) %>%
    mutate(
      Participant_id_A = as.character(Participant_id_A),
      Timepoint_A = as.character(Timepoint_A),
      Participant_id_B = as.character(Participant_id_B),
      Timepoint_B = as.character(Timepoint_B)
    )
  direct <- readr::read_csv(path_direct_metrics, show_col_types = FALSE) %>%
    transmute(
      Participant_id_A = as.character(Participant_id_A),
      Timepoint_A = as.character(Timepoint_A),
      Participant_id_B = as.character(Participant_id_B),
      Timepoint_B = as.character(Timepoint_B),
      direct_avg_identity = AvgIdentity,
      direct_total_snps = TotalSNPs,
      direct_cache_status = dnadiff_cache_status,
      direct_report_path = dnadiff_report_path,
      direct_report_sha256 = dnadiff_report_sha256
    )

  pair_tbl <- preliminary %>%
    left_join(
      direct,
      by = c(
        "Participant_id_A", "Timepoint_A",
        "Participant_id_B", "Timepoint_B"
      )
    ) %>%
    mutate(
      priority_category = case_when(
        !is.na(direct_total_snps) & direct_total_snps > 25 ~
          "Same lineage, not the same strain",
        vf_identical &
          ST_reliable_A %in% TRUE & ST_reliable_B %in% TRUE &
          !is.na(direct_total_snps) & direct_total_snps <= 25 ~
          "High-priority candidate",
        vf_jaccard >= 0.95 &
          ST_reliable_A %in% TRUE & ST_reliable_B %in% TRUE ~
          "Moderate-priority candidate",
        vf_identical & !is.na(direct_total_snps) & direct_total_snps <= 25 ~
          "Moderate-priority candidate",
        TRUE ~ "Possible related lineage"
      ),
      priority_reason = case_when(
        priority_category == "High-priority candidate" ~
          "Same reliable ST; exact 227-gene VF profile; core SNP <=25; direct SNP <=25.",
        priority_category == "Moderate-priority candidate" & vf_identical ~
          "Exact VF and direct SNP support, but provider ST-quality support is incomplete.",
        priority_category == "Moderate-priority candidate" ~
          "Same reliable ST; VF Jaccard >=0.95; core SNP <=25; direct validation pending.",
        priority_category == "Same lineage, not the same strain" ~
          "Same ST and exact VF profile, but direct SNP distance is above 25.",
        TRUE ~
          "Same ST; VF Jaccard 0.90-0.949 or incomplete ST-quality support; core SNP <=25."
      ),
      direct_validation_status = case_when(
        !vf_identical ~ "Not run: reserved for exact-VF preliminary pairs",
        is.na(direct_total_snps) ~ "Attempted but unavailable",
        direct_total_snps <= 25 ~ "Passed project operational <=25-SNP rule",
        TRUE ~ "Did not pass project operational <=25-SNP rule"
      )
    )

  amr_profiles <- build_amr_profiles(inputs$meta)
  pair_tbl <- pair_tbl %>%
    add_pair_profiles(inputs, amr_profiles) %>%
    assign_groups() %>%
    select(
      candidate_cluster_id, screen_group_id,
      priority_category, priority_reason,
      pair_id, sample_key_A, sample_key_B,
      Participant_id_A, Timepoint_A, Episode_ID_A, Isolate_ID_A,
      Collection_Date_A, Event_type_A, UTI_Status_A, Primary_Status_A,
      Participant_id_B, Timepoint_B, Episode_ID_B, Isolate_ID_B,
      Collection_Date_B, Event_type_B, UTI_Status_B, Primary_Status_B,
      collection_date_gap_days,
      ST_A, ST_source_A, ST_reliable_A, provider_PercGoodTargets_A,
      ST_B, ST_source_B, ST_reliable_B, provider_PercGoodTargets_B,
      preliminary_core_snps, preliminary_core_call,
      direct_total_snps, direct_avg_identity, direct_validation_status,
      direct_cache_status, direct_report_path, direct_report_sha256,
      vf_jaccard, vf_identical, vf_intersection_n, vf_union_n,
      vf_discordant_gene_n, vf_genes_only_A, vf_genes_only_B,
      amr_jaccard_excluding_mdfA,
      shared_amr_genes_excluding_mdfA,
      amr_genes_only_A_excluding_mdfA,
      amr_genes_only_B_excluding_mdfA,
      all_amr_genes_A, all_amr_genes_B,
      mdfA_detected_A, mdfA_detected_B,
      shared_amr_classes, all_amr_classes_A, all_amr_classes_B
    ) %>%
    arrange(
      factor(
        priority_category,
        levels = c(
          "High-priority candidate",
          "Moderate-priority candidate",
          "Possible related lineage",
          "Same lineage, not the same strain"
        )
      ),
      candidate_cluster_id,
      screen_group_id,
      direct_total_snps,
      preliminary_core_snps,
      desc(vf_jaccard)
    )

  clusters <- build_cluster_summary(pair_tbl, inputs, amr_profiles)
  isolates <- build_isolate_profiles(pair_tbl, inputs, amr_profiles)
  excluded <- build_excluded_isolates(inputs)
  guide <- build_interpretation_guide()

  write_csv(pair_tbl, path_candidate_pairs, na = "")
  write_csv(clusters, path_cluster_summary, na = "")
  write_csv(isolates, path_isolate_profiles, na = "")
  write_csv(excluded, path_excluded, na = "")
  write_csv(guide, path_guide, na = "")

  summary <- list(
    generated_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
    active_isolates = nrow(inputs$meta),
    active_participants = n_distinct(inputs$meta$Participant_id),
    emailed_rows = nrow(readr::read_csv(opt$external_vf, show_col_types = FALSE)),
    excluded_emailed_rows = nrow(excluded),
    vf_gene_count = length(inputs$gene_cols),
    preliminary_candidate_pairs = nrow(pair_tbl),
    direct_pairs_attempted = sum(pair_tbl$vf_identical),
    direct_pairs_with_snp_result = sum(!is.na(pair_tbl$direct_total_snps)),
    priority_counts = as.list(table(pair_tbl$priority_category)),
    high_priority_cluster_count = nrow(clusters),
    high_priority_cluster_isolate_count = n_distinct(
      c(
        pair_tbl$sample_key_A[pair_tbl$priority_category == "High-priority candidate"],
        pair_tbl$sample_key_B[pair_tbl$priority_category == "High-priority candidate"]
      )
    ),
    amr_screen_method = "ABRicate ResFinder; retained identity >=80% and coverage >=80%",
    phenotype_available = FALSE,
    phenotype_note = "No phenotypic susceptibility or antibiogram data were available.",
    caution = "Candidate genomic clusters require nursing-home, ward and temporal confirmation."
  )
  jsonlite::write_json(summary, path_summary_json, auto_unbox = TRUE, pretty = TRUE)
  message(
    "Finalized ", nrow(pair_tbl), " candidate pairs and ",
    nrow(clusters), " high-priority candidate genomic clusters."
  )
}

inputs <- read_core_inputs()

if (opt$stage %in% c("prepare", "all")) {
  prepare_stage(inputs)
}

if (opt$stage == "all") {
  run_direct_validation()
}

if (opt$stage %in% c("finalize", "all")) {
  finalize_stage(inputs)
}

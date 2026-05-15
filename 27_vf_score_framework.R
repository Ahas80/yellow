#!/usr/bin/env Rscript
# ==============================================================================
# 27_vf_score_framework.R
# ==============================================================================
#
# GOAL:
#   Compare multiple VF scoring approaches: total burden, curated burden,
#   module count, UPEC-associated scores, and exploratory PCA/ordination.
#
# METHOD:
#   1. Load canonical VF dataset (script 22) and module outputs (script 26).
#   2. Calculate score variants per episode.
#   3. Summarise by Infection_Status and ST.
#   4. Run exploratory comparisons and PCA.
#
# OUTPUT:
#   - results/vf/vf_score_table.csv
#   - results/vf/vf_module_score_table.csv
#   - results/vf/vf_upec_score_components.csv
#   - results/vf/vf_score_summary_by_status.csv
#   - results/vf/vf_score_summary_by_ST.csv
#   - results/vf/vf_score_summary_by_status_within_ST.csv
#   - results/vf/vf_score_correlations.csv
#   - results/vf/vf_score_tests_exploratory.csv
#   - results/vf/vf_pca_coordinates.csv
#   - results/vf/vf_pca_loadings.csv
#   - results/vf/vf_pcoa_jaccard_coordinates.csv [if supported]
#   - results/vf/vf_score_framework_summary.txt
#   - plots/vf/vf_scores_by_status.png
#   - plots/vf/vf_scores_by_ST.png
#   - plots/vf/vf_score_correlation_heatmap.png
#   - plots/vf/vf_pca_status.png
#   - plots/vf/vf_pca_ST.png
#
# CAUTION:
#   UTI count is small (~16). All comparisons are exploratory/descriptive.
#   Do not interpret as validated predictive scores.
# ==============================================================================

source("00_config.R")
suppressPackageStartupMessages({
  library(dplyr); library(readr); library(tidyr); library(stringr)
  library(ggplot2); library(tibble)
})

msg("Starting 27_vf_score_framework.R")

plot_theme_vf <- function(base_size = 11) {
  theme_bw(base_size = base_size) +
    theme(
      plot.caption = element_text(hjust = 0, size = base_size - 3, colour = "grey35"),
      plot.subtitle = element_text(colour = "grey25"),
      legend.position = "bottom"
    )
}

normalise_st_label <- function(x) {
  x <- as.character(x)
  x <- str_trim(x)
  unknown <- c("", "-", "ST-", "NA", "N/A", "UNKNOWN", "UNK", "NT",
               "NON-TYPABLE", "NONTYPABLE", "NOT TYPED")
  x[str_to_upper(x) %in% unknown] <- NA_character_
  x
}

# ==============================================================================
# 1. PATHS
# ==============================================================================
FILE_MODULE_MAP  <- file.path(DIR_VF, "gene_module_map.csv")
FILE_MODULE_EP   <- file.path(DIR_VF, "vf_module_presence_by_episode.csv")
FILE_VF_PA_RAW    <- FILE_VF_PA

ensure_dir(DIR_VF); ensure_dir(DIR_PLOTS_VF)

# ==============================================================================
# 2. LOAD DATA
# ==============================================================================
if (!file.exists(FILE_VF_READY)) stop("Missing ", FILE_VF_READY)
if (!file.exists(FILE_MODULE_MAP)) stop("Missing ", FILE_MODULE_MAP, ". Run 26_vf_define_gene_modules.R first.")
if (!file.exists(FILE_MODULE_EP)) stop("Missing ", FILE_MODULE_EP, ". Run 26_vf_define_gene_modules.R first.")

stop_if_stale <- function(target, upstream, target_label, upstream_label) {
  if (!file.exists(target) || !file.exists(upstream)) return(invisible(FALSE))
  target_mtime <- file.info(target)$mtime
  upstream_mtime <- file.info(upstream)$mtime
  if (!is.na(target_mtime) && !is.na(upstream_mtime) && target_mtime < upstream_mtime) {
    if (basename(target) == "vf_analysis_ready.csv" && basename(upstream) == "vf_pa_all.csv") {
      target_probe <- read_csv(target, show_col_types = FALSE, n_max = Inf)
      upstream_probe <- read_csv(upstream, show_col_types = FALSE, n_max = Inf)
      target_genes <- canonical_vf_gene_cols(names(target_probe), vf_pa_file = upstream)
      upstream_genes <- canonical_vf_gene_cols(names(upstream_probe), vf_pa_file = upstream)
      if (nrow(target_probe) == nrow(upstream_probe) && setequal(target_genes, upstream_genes)) {
        msg("WARNING: %s timestamp is older than %s, but row count and gene set match; continuing.",
            target_label, upstream_label)
        return(invisible(FALSE))
      }
    }
    stop(sprintf(
      "%s is older than %s. Re-run 22_vf_build_analysis_dataset.R and 26_vf_define_gene_modules.R before script 27.\n  %s: %s\n  %s: %s",
      target_label, upstream_label, target_label, format(target_mtime),
      upstream_label, format(upstream_mtime)
    ))
  }
}

stop_if_stale(FILE_VF_READY, FILE_VF_PA_RAW, "vf_analysis_ready.csv", "vf_pa_all.csv")
stop_if_stale(FILE_MODULE_EP, FILE_VF_READY, "vf_module_presence_by_episode.csv", "vf_analysis_ready.csv")

vf <- read_csv(FILE_VF_READY, show_col_types = FALSE) %>%
  mutate(Participant_id = as.character(Participant_id),
         tp_lab = normalise_timepoint_preserve_events(tp_lab),
         ST = if ("ST" %in% names(.)) normalise_st_label(ST) else NA_character_)
mod_map <- read_csv(FILE_MODULE_MAP, show_col_types = FALSE)
mod_ep  <- read_csv(FILE_MODULE_EP, show_col_types = FALSE) %>%
  mutate(Participant_id = as.character(Participant_id))

gene_cols <- canonical_vf_gene_cols(names(vf), vf_pa_file = FILE_VF_PA_RAW)
meta_cols <- intersect(
  c(
    "Participant_id", "tp_lab", "Episode_ID", "Collection_Date",
    "Infection_Status", "Batch", "Status_Confidence_epi", "Sx_source_epi",
    "UTI_Label", "Urine_collection_method", "ST", "vf_count_total",
    "total_vf_count_all", "total_vf_count_curated",
    "total_vf_count_upec_candidate", "total_vf_count_unassigned",
    "low_confidence_count", "is_ecoli", "n_timepoints"
  ),
  names(vf)
)

n_asb <- sum(vf$Infection_Status == "ASB", na.rm = TRUE)
n_uti <- sum(vf$Infection_Status == "UTI", na.rm = TRUE)
n_neg <- sum(vf$Infection_Status == "Negative", na.rm = TRUE)
msg("Loaded: %d episodes (ASB=%d, UTI=%d, Neg=%d), %d gene cols",
    nrow(vf), n_asb, n_uti, n_neg, length(gene_cols))

# ==============================================================================
# 3. CALCULATE SCORES
# ==============================================================================
# Curated genes: High/Moderate confidence, primary, not unassigned
curated_genes <- mod_map %>%
  filter(in_vf_matrix, primary_assignment,
         assignment_confidence %in% c("High", "Moderate"),
         module_id != "unassigned") %>%
  pull(Gene) %>% intersect(gene_cols)

upec_genes <- mod_map %>%
  filter(in_vf_matrix, primary_assignment, upec_score_candidate,
         assignment_confidence %in% c("High", "Moderate")) %>%
  pull(Gene) %>% intersect(gene_cols)

# Unassigned genes
unassigned_genes <- mod_map %>%
  filter(in_vf_matrix, primary_assignment, module_id == "unassigned") %>%
  pull(Gene) %>% intersect(gene_cols)

# Low confidence genes
low_conf_genes <- mod_map %>%
  filter(in_vf_matrix, primary_assignment,
         assignment_confidence == "Low") %>%
  pull(Gene) %>% intersect(gene_cols)

msg("Score components: curated=%d, UPEC=%d, unassigned=%d, low_conf=%d",
    length(curated_genes), length(upec_genes), length(unassigned_genes), length(low_conf_genes))

# Module-level scores from mod_ep
pres_cols <- grep("_present$", names(mod_ep), value = TRUE)

# UPEC modules
upec_modules <- mod_map %>%
  filter(upec_score_candidate, module_id != "unassigned",
         assignment_confidence %in% c("High","Moderate")) %>%
  pull(module_id) %>% unique()
upec_pres_cols <- paste0("mod_", upec_modules, "_present")
upec_pres_cols <- intersect(upec_pres_cols, names(mod_ep))

scores <- vf %>%
  select(any_of(meta_cols)) %>%
  mutate(
    total_vf_count_all = if ("total_vf_count_all" %in% names(.)) total_vf_count_all else vf_count_total,
    vf_count_curated = rowSums(vf[, curated_genes, drop = FALSE], na.rm = TRUE),
    vf_count_unassigned = rowSums(vf[, unassigned_genes, drop = FALSE], na.rm = TRUE),
    vf_count_low_conf = rowSums(vf[, low_conf_genes, drop = FALSE], na.rm = TRUE),
    upec_gene_score = rowSums(vf[, upec_genes, drop = FALSE], na.rm = TRUE),
    total_vf_count_curated = vf_count_curated,
    total_vf_count_upec_candidate = upec_gene_score,
    total_vf_count_unassigned = vf_count_unassigned,
    low_confidence_count = vf_count_low_conf
  ) %>%
  left_join(
    mod_ep %>% select(Participant_id, tp_lab, all_of(pres_cols),
                      n_modules_present, n_upec_modules_present),
    by = c("Participant_id", "tp_lab")
  ) %>%
  mutate(
    upec_n_possible_genes = length(upec_genes),
    upec_gene_fraction = ifelse(upec_n_possible_genes > 0,
                                round(upec_gene_score / upec_n_possible_genes, 3), NA_real_),
    upec_n_possible_modules = length(upec_modules),
    upec_module_fraction = ifelse(upec_n_possible_modules > 0,
                                  round(n_upec_modules_present / upec_n_possible_modules, 3), NA_real_),
    # Compatibility names used by later scripts and thesis tables.
    upec_gene_score_unweighted = upec_gene_score,
    upec_module_score_unweighted = n_upec_modules_present,
    module_count_present = n_modules_present,
    module_count_curated = n_modules_present,
    upec_score_fraction = upec_module_fraction
  )

msg("Score table: %d rows × %d columns", nrow(scores), ncol(scores))

# ==============================================================================
# 3b. MODULE SCORE TABLE
# ==============================================================================
count_cols <- grep("_n_genes$", names(mod_ep), value = TRUE)
module_scores <- mod_ep %>%
  select(any_of(meta_cols), all_of(pres_cols), all_of(count_cols),
         n_modules_present, n_upec_modules_present) %>%
  mutate(
    module_count_present = n_modules_present,
    upec_module_score_unweighted = n_upec_modules_present
  )

module_lookup <- mod_map %>%
  filter(primary_assignment, in_vf_matrix, module_id != "unassigned") %>%
  distinct(module_id, broad_module)

for (bm in sort(unique(module_lookup$broad_module))) {
  bm_modules <- module_lookup %>% filter(broad_module == bm) %>% pull(module_id)
  bm_cols <- paste0("mod_", bm_modules, "_n_genes")
  bm_cols <- intersect(bm_cols, names(module_scores))
  safe_bm <- str_replace_all(str_to_lower(bm), "[^a-z0-9]+", "_")
  safe_bm <- str_replace_all(safe_bm, "^_|_$", "")
  out_col <- paste0("module_score_", safe_bm)
  module_scores[[out_col]] <- if (length(bm_cols) > 0) {
    rowSums(module_scores[, bm_cols, drop = FALSE], na.rm = TRUE)
  } else {
    0
  }
}

# ==============================================================================
# 4. UPEC SCORE COMPONENTS TABLE
# ==============================================================================
upec_components <- mod_map %>%
  filter(upec_score_candidate, in_vf_matrix, primary_assignment) %>%
  mutate(
    detected_in_n = sapply(Gene, function(g) if (g %in% gene_cols) sum(vf[[g]], na.rm = TRUE) else 0L),
    detected_pct = round(100 * detected_in_n / nrow(vf), 1)
  ) %>%
  select(Gene, module_id, system_name, broad_module, assignment_confidence,
         upec_score_candidate, detected_in_n, detected_pct, notes) %>%
  arrange(broad_module, module_id, Gene)

# ==============================================================================
# 5. SUMMARISE SCORES BY STATUS
# ==============================================================================
score_names <- c("total_vf_count_all", "total_vf_count_curated",
                 "total_vf_count_upec_candidate", "total_vf_count_unassigned",
                 "low_confidence_count", "vf_count_total", "vf_count_curated",
                 "vf_count_unassigned", "upec_gene_score_unweighted",
                 "module_count_present", "upec_module_score_unweighted",
                 "upec_score_fraction")
score_names <- intersect(score_names, names(scores))

summarise_scores <- function(df, group_col) {
  results <- tibble()
  for (sc in score_names) {
    grp <- df %>%
      filter(!is.na(.data[[group_col]]), !is.na(.data[[sc]])) %>%
      group_by(.data[[group_col]]) %>%
      summarise(
        score_name = sc,
        n_episodes = n(),
        n_participants = n_distinct(Participant_id),
        median = median(.data[[sc]], na.rm = TRUE),
        q25 = quantile(.data[[sc]], 0.25, na.rm = TRUE),
        q75 = quantile(.data[[sc]], 0.75, na.rm = TRUE),
        mean = round(mean(.data[[sc]], na.rm = TRUE), 2),
        sd = round(sd(.data[[sc]], na.rm = TRUE), 2),
        min = min(.data[[sc]], na.rm = TRUE),
        max = max(.data[[sc]], na.rm = TRUE),
        .groups = "drop"
      )
    names(grp)[1] <- group_col
    results <- bind_rows(results, grp)
  }
  results
}

by_status <- summarise_scores(scores, "Infection_Status")
by_st     <- summarise_scores(scores, "ST")

common_sts_for_status <- scores %>%
  filter(!is.na(ST), !is.na(Infection_Status)) %>%
  dplyr::count(ST, Infection_Status, name = "n_status") %>%
  group_by(ST) %>%
  summarise(
    n_episodes = sum(n_status),
    n_statuses = n_distinct(Infection_Status),
    has_ASB = any(Infection_Status == "ASB"),
    has_UTI = any(Infection_Status == "UTI"),
    .groups = "drop"
  ) %>%
  filter(n_episodes >= 5, has_ASB, has_UTI) %>%
  pull(ST)

by_status_within_st <- tibble()
if (length(common_sts_for_status) > 0) {
  by_status_within_st <- scores %>%
    filter(ST %in% common_sts_for_status, !is.na(Infection_Status)) %>%
    select(Participant_id, ST, Infection_Status, all_of(score_names)) %>%
    pivot_longer(cols = all_of(score_names), names_to = "score_name", values_to = "value") %>%
    group_by(ST, Infection_Status, score_name) %>%
    summarise(
      n_episodes = n(),
      n_participants = n_distinct(Participant_id),
      median = median(value, na.rm = TRUE),
      q25 = quantile(value, 0.25, na.rm = TRUE),
      q75 = quantile(value, 0.75, na.rm = TRUE),
      mean = round(mean(value, na.rm = TRUE), 2),
      sd = round(sd(value, na.rm = TRUE), 2),
      min = min(value, na.rm = TRUE),
      max = max(value, na.rm = TRUE),
      .groups = "drop"
    )
}

# ==============================================================================
# 6. EXPLORATORY TESTS (ASB vs UTI)
# ==============================================================================
test_results <- tibble()
asb_uti <- scores %>% filter(Infection_Status %in% c("ASB", "UTI"))

for (sc in score_names) {
  x_asb <- asb_uti %>% filter(Infection_Status == "ASB") %>% pull(!!sym(sc))
  x_uti <- asb_uti %>% filter(Infection_Status == "UTI") %>% pull(!!sym(sc))
  if (length(x_asb) >= 3 && length(x_uti) >= 3) {
    wt <- tryCatch(wilcox.test(x_uti, x_asb, exact = FALSE), error = function(e) NULL)
    test_results <- bind_rows(test_results, tibble(
      comparison = "UTI vs ASB",
      score_name = sc,
      test = "Wilcoxon rank-sum",
      n_ASB = length(x_asb), n_UTI = length(x_uti),
      median_ASB = median(x_asb, na.rm = TRUE),
      median_UTI = median(x_uti, na.rm = TRUE),
      p_value = if (!is.null(wt)) wt$p.value else NA_real_,
      note = "Exploratory; repeated measures not accounted for"
    ))
  }
}
if (nrow(test_results) > 0) {
  test_results <- test_results %>%
    mutate(p_adj_BH = p.adjust(p_value, method = "BH"))
}

# ==============================================================================
# 7. SCORE CORRELATIONS
# ==============================================================================
score_mat <- scores[, score_names, drop = FALSE]
cor_results <- tibble()
pairs <- combn(score_names, 2, simplify = FALSE)
for (pr in pairs) {
  valid <- is.finite(score_mat[[pr[1]]]) & is.finite(score_mat[[pr[2]]])
  x <- score_mat[[pr[1]]][valid]
  y <- score_mat[[pr[2]]][valid]
  if (length(x) < 3 || length(unique(x)) < 2 || length(unique(y)) < 2) {
    cor_results <- bind_rows(cor_results, tibble(
      score_x = pr[1], score_y = pr[2],
      method = "Spearman", rho = NA_real_, p_value = NA_real_,
      n = length(x), note = "Skipped: insufficient variation"
    ))
    next
  }
  ct <- cor.test(x, y, method = "spearman", exact = FALSE)
  cor_results <- bind_rows(cor_results, tibble(
    score_x = pr[1], score_y = pr[2],
    method = "Spearman", rho = round(ct$estimate, 3),
    p_value = ct$p.value, n = length(x), note = "Exploratory"
  ))
}
if (nrow(cor_results) > 0) {
  cor_results <- cor_results %>%
    mutate(p_adj_BH = ifelse(is.na(p_value), NA_real_, p.adjust(p_value, method = "BH")))
}

# ==============================================================================
# 8. PCA ON MODULE PROFILES
# ==============================================================================
pca_coords <- tibble()
pca_loadings <- tibble()
pcoa_coords <- tibble()
pca_ok <- FALSE
pcoa_ok <- FALSE

mod_pres_cols <- grep("^mod_.*_present$", names(mod_ep), value = TRUE)
if (length(mod_pres_cols) >= 3) {
  pca_mat <- mod_ep[, mod_pres_cols, drop = FALSE]
  # Remove zero-variance columns
  var_ok <- apply(pca_mat, 2, var, na.rm = TRUE) > 0
  pca_mat <- pca_mat[, var_ok, drop = FALSE]

  if (ncol(pca_mat) >= 3) {
    pc <- prcomp(pca_mat, center = TRUE, scale. = FALSE)
    ve <- round(100 * summary(pc)$importance[2, 1:min(3, ncol(pc$x))], 1)
    pca_coords <- bind_cols(
      mod_ep %>% select(Participant_id, tp_lab, Infection_Status, ST),
      as_tibble(pc$x[, 1:min(3, ncol(pc$x))])
    ) %>%
      mutate(var_PC1 = ve[1], var_PC2 = ve[2])
    load_cols <- paste0("PC", seq_len(min(3, ncol(pc$rotation))))
    pca_loadings <- as_tibble(pc$rotation[, seq_along(load_cols), drop = FALSE], rownames = "feature") %>%
      setNames(c("feature", paste0(load_cols, "_loading"))) %>%
      mutate(
        module_id = str_replace_all(str_replace_all(feature, "^mod_", ""), "_present$", "")
      ) %>%
      left_join(
        mod_map %>%
          filter(primary_assignment, in_vf_matrix) %>%
          distinct(module_id, system_name, broad_module),
        by = "module_id"
      )
    pca_ok <- TRUE
    msg("PCA: %d features, PC1=%.1f%%, PC2=%.1f%%", ncol(pca_mat), ve[1], ve[2])
  }
}
if (!pca_ok) msg("PCA skipped: insufficient variable modules")

# Optional PCoA on binary module presence using Jaccard distance. This is a
# sensitivity ordination for binary profiles, not a formal hypothesis test.
if (length(mod_pres_cols) >= 3) {
  pcoa_mat <- as.matrix(mod_ep[, mod_pres_cols, drop = FALSE])
  var_ok <- apply(pcoa_mat, 2, var, na.rm = TRUE) > 0
  pcoa_mat <- pcoa_mat[, var_ok, drop = FALSE]
  if (ncol(pcoa_mat) >= 3 && nrow(pcoa_mat) >= 3) {
    pcoa_mat[is.na(pcoa_mat)] <- 0
    inter <- pcoa_mat %*% t(pcoa_mat)
    row_totals <- rowSums(pcoa_mat)
    union <- outer(row_totals, row_totals, "+") - inter
    dist_mat <- 1 - inter / union
    dist_mat[union == 0] <- 0
    diag(dist_mat) <- 0
    pcj <- cmdscale(as.dist(dist_mat), k = 3, eig = TRUE)
    if (!is.null(pcj$points) && ncol(pcj$points) >= 2) {
      eig_pos <- pcj$eig[pcj$eig > 0]
      var_axes <- if (length(eig_pos) >= 2) round(100 * pcj$eig[1:2] / sum(eig_pos), 1) else c(NA_real_, NA_real_)
      point_mat <- pcj$points[, 1:min(3, ncol(pcj$points)), drop = FALSE]
      colnames(point_mat) <- paste0("Axis", seq_len(ncol(point_mat)))
      pcoa_coords <- bind_cols(
        mod_ep %>% select(Participant_id, tp_lab, Infection_Status, ST),
        as_tibble(point_mat)
      )
      pcoa_coords <- pcoa_coords %>% mutate(var_Axis1 = var_axes[1], var_Axis2 = var_axes[2])
      pcoa_ok <- TRUE
      msg("Jaccard PCoA: %d features, Axis1=%.1f%%, Axis2=%.1f%%", ncol(pcoa_mat), var_axes[1], var_axes[2])
    }
  }
}
if (!pcoa_ok) msg("Jaccard PCoA skipped: insufficient variable modules")

# ==============================================================================
# 9. WRITE OUTPUTS
# ==============================================================================
write_csv(scores, file.path(DIR_VF, "vf_score_table.csv"))
write_csv(module_scores, file.path(DIR_VF, "vf_module_score_table.csv"))
write_csv(upec_components, file.path(DIR_VF, "vf_upec_score_components.csv"))
write_csv(by_status, file.path(DIR_VF, "vf_score_summary_by_status.csv"))
write_csv(by_st, file.path(DIR_VF, "vf_score_summary_by_ST.csv"))
write_csv(by_status_within_st, file.path(DIR_VF, "vf_score_summary_by_status_within_ST.csv"))
write_csv(test_results, file.path(DIR_VF, "vf_score_tests_exploratory.csv"))
write_csv(cor_results, file.path(DIR_VF, "vf_score_correlations.csv"))
if (pca_ok) write_csv(pca_coords, file.path(DIR_VF, "vf_pca_coordinates.csv"))
if (pca_ok) write_csv(pca_loadings, file.path(DIR_VF, "vf_pca_loadings.csv"))
if (pcoa_ok) write_csv(pcoa_coords, file.path(DIR_VF, "vf_pcoa_jaccard_coordinates.csv"))

# Summary text
summ <- character()
sa <- function(...) summ <<- c(summ, sprintf(...))
sa("=== VF SCORE FRAMEWORK SUMMARY ===")
sa("Timestamp: %s", format(Sys.time()))
sa("Episodes: %d (ASB=%d, UTI=%d, Neg=%d)", nrow(scores), n_asb, n_uti, n_neg)
sa("Curated genes: %d, UPEC genes: %d", length(curated_genes), length(upec_genes))
sa("Unassigned genes: %d, low-confidence genes: %d", length(unassigned_genes), length(low_conf_genes))
sa("UPEC modules: %d", length(upec_modules))
sa("Within-ST status summaries: %d STs with >=5 episodes and both ASB/UTI",
   length(common_sts_for_status))
sa("PCA performed: %s", pca_ok)
sa("Jaccard PCoA performed: %s", pcoa_ok)
sa("")
sa("--- Score medians by status ---")
for (sc in score_names) {
  row_asb <- by_status %>% filter(Infection_Status == "ASB", score_name == sc)
  row_uti <- by_status %>% filter(Infection_Status == "UTI", score_name == sc)
  if (nrow(row_asb) > 0 && nrow(row_uti) > 0) {
    sa("  %s: ASB median=%.0f, UTI median=%.0f", sc, row_asb$median, row_uti$median)
  }
}
sa("")
sa("CAUTION: UTI n=%d. All comparisons are exploratory.", n_uti)
if (length(unassigned_genes) / max(1, length(gene_cols)) > 0.25) {
  sa("WARNING: Unassigned genes are %.1f%% of the VF matrix; interpret total burden separately from curated/UPEC-candidate scores.",
     100 * length(unassigned_genes) / max(1, length(gene_cols)))
}
sa("Repeated measures from same participants are not adjusted.")
writeLines(summ, file.path(DIR_VF, "vf_score_framework_summary.txt"))

cat(
  sprintf("\nScore warning: unassigned genes represent %.1f%% of the VF matrix; curated/UPEC candidate scores are separated from total counts.\n",
          100 * length(unassigned_genes) / max(1, length(gene_cols))),
  file = file.path(DIR_VF, "vf_gene_annotation_gap_report.txt"),
  append = TRUE
)

append_denominator_summary(
  scores,
  "27_vf_score_framework.R",
  "vf_score_table",
  "participant_timepoint",
  file.path(DIR_VF, "vf_score_table.csv"),
  "Score table separates all, curated, UPEC-candidate, unassigned, and low-confidence VF burden"
)

# ==============================================================================
# 10. PLOTS
# ==============================================================================
# Scores by status
score_long <- scores %>%
  filter(Infection_Status %in% c("ASB","UTI","Negative")) %>%
  pivot_longer(cols = all_of(score_names), names_to = "score", values_to = "value")

p1 <- ggplot(score_long, aes(x = Infection_Status, y = value, fill = Infection_Status)) +
  geom_boxplot(outlier.size = 0.8, alpha = 0.7) +
  facet_wrap(~score, scales = "free_y") +
  labs(
    title = "Exploratory VF score distributions by clinical status",
    subtitle = sprintf("ASB n=%d, UTI n=%d, Negative n=%d; scores are descriptive, not validated predictors", n_asb, n_uti, n_neg),
    x = "Clinical status",
    y = "Score value",
    caption = sprintf(
      "Data: %s and module outputs from script 26. Denominator: %d VF/WGS-linked E. coli isolates from %d participants. Level of analysis: isolate-level score summary. Residents may contribute repeated isolates; score comparisons are exploratory unless modelled with Participant_id clustering. UTI n=%d is small and ST/lineage may confound score-status differences.",
      FILE_VF_READY, nrow(vf), n_distinct(vf$Participant_id), n_uti
    )
  ) +
  plot_theme_vf(base_size = 11) + theme(legend.position = "none")
ggsave(file.path(DIR_PLOTS_VF, "vf_scores_by_status.png"), p1, width = 12, height = 8.5, dpi = 300)

score_effects <- test_results %>%
  filter(comparison == "UTI vs ASB") %>%
  mutate(
    median_difference_UTI_minus_ASB = median_UTI - median_ASB,
    score_label = str_replace_all(score_name, "_", " "),
    q_label = sprintf("q=%.2g", p_adj_BH),
    result_status = ifelse(!is.na(p_adj_BH) & p_adj_BH < 0.05,
                           "BH q < 0.05", "Not BH-significant")
  ) %>%
  arrange(median_difference_UTI_minus_ASB)

if (nrow(score_effects) > 0) {
  p_score_effect <- ggplot(score_effects,
                           aes(x = median_difference_UTI_minus_ASB,
                               y = reorder(score_label, median_difference_UTI_minus_ASB),
                               fill = result_status)) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey45") +
    geom_col(width = 0.68, colour = "white", linewidth = 0.2) +
    geom_text(aes(label = q_label),
              hjust = ifelse(score_effects$median_difference_UTI_minus_ASB >= 0, -0.08, 1.08),
              size = 3, colour = "grey20") +
    scale_fill_manual(values = c("BH q < 0.05" = "#D55E00", "Not BH-significant" = "grey55")) +
    labs(
      title = "Exploratory ASB-UTI VF score differences",
      subtitle = sprintf("Median UTI minus ASB score differences; ASB n=%d, UTI n=%d", n_asb, n_uti),
      x = "Median difference (UTI - ASB)",
      y = "VF score",
      fill = "Adjusted result",
      caption = sprintf(
        "Data: %s. Denominator: ASB n=%d and UTI n=%d VF/WGS-linked isolates. Level of analysis: isolate-level score comparison. Wilcoxon tests and BH-adjusted q-values are exploratory because residents may contribute repeated isolates and ST/lineage, timepoint, batch, and event-type structure are not modelled here. Scores are descriptive summaries, not validated UTI predictors.",
        file.path(DIR_VF, "vf_score_tests_exploratory.csv"), n_asb, n_uti
      )
    ) +
    plot_theme_vf(base_size = 11) +
    theme(legend.position = "bottom")

  ggsave(file.path(DIR_PLOTS_VF, "vf_score_effect_summary_asb_uti.png"),
         p_score_effect, width = 8.8, height = max(5.2, nrow(score_effects) * 0.34), dpi = 300)
}

# Scores by top STs
top_sts <- scores %>% dplyr::count(ST) %>% filter(n >= 5) %>% pull(ST)
if (length(top_sts) >= 3) {
  st_long <- scores %>%
    filter(ST %in% top_sts) %>%
    pivot_longer(cols = all_of(score_names), names_to = "score", values_to = "value")
  p2 <- ggplot(st_long, aes(x = reorder(ST, value, FUN = median), y = value)) +
    geom_boxplot(fill = "steelblue", alpha = 0.6, outlier.size = 0.5) +
    facet_wrap(~score, scales = "free_y") +
    coord_flip() +
    labs(
      title = "Exploratory VF score distributions by E. coli sequence type",
      subtitle = "STs with at least five VF-ready isolates are shown",
      x = "Sequence type",
      y = "Score value",
      caption = sprintf(
        "Data: %s. Denominator: common STs among %d VF-ready isolates. Level of analysis: isolate-level lineage diagnostic. Residents may contribute repeated isolates; ST-associated score differences are expected because VF content is lineage structured.",
        FILE_VF_READY, nrow(vf)
      )
    ) +
    plot_theme_vf(base_size = 10)
  ggsave(file.path(DIR_PLOTS_VF, "vf_scores_by_ST.png"), p2, width = 12, height = 8.5, dpi = 300)
}

# Correlation heatmap
if (nrow(cor_results) > 0) {
  all_scores <- unique(c(cor_results$score_x, cor_results$score_y))
  cor_full <- bind_rows(
    cor_results %>% select(score_x, score_y, rho),
    cor_results %>% transmute(score_x = score_y, score_y = score_x, rho),
    tibble(score_x = all_scores, score_y = all_scores, rho = 1.0)
  )
  p3 <- ggplot(cor_full, aes(x = score_x, y = score_y, fill = rho)) +
    geom_tile() + geom_text(aes(label = rho), size = 3) +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
    labs(
      title = "Spearman correlations among exploratory VF scores",
      subtitle = "Correlation reflects shared gene/module content and should not be interpreted as clinical prediction",
      x = NULL,
      y = NULL,
      caption = sprintf("Data: %s. Level of analysis: score-score correlation across VF-ready isolates; repeated residents and lineage structure are not modelled.", FILE_VF_READY)
    ) +
    plot_theme_vf(base_size = 9) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  ggsave(file.path(DIR_PLOTS_VF, "vf_score_correlation_heatmap.png"), p3, width = 8, height = 7.5, dpi = 300)
}

# PCA plots
if (pca_ok && nrow(pca_coords) > 0) {
  p4 <- ggplot(pca_coords, aes(x = PC1, y = PC2, colour = Infection_Status)) +
    geom_point(alpha = 0.7, size = 2) +
    labs(
      title = "Exploratory PCA of VF module profiles by clinical status",
      subtitle = "Ordination is descriptive and does not adjust for repeated residents or ST/lineage",
      x = sprintf("PC1 (%.1f%%)", pca_coords$var_PC1[1]),
      y = sprintf("PC2 (%.1f%%)", pca_coords$var_PC2[1]),
      caption = sprintf("Data: %s and script 26 module matrix. Denominator: VF-ready isolates with module profiles; UTI n=%d is small.", FILE_VF_READY, n_uti)
    ) +
    plot_theme_vf(base_size = 11)
  ggsave(file.path(DIR_PLOTS_VF, "vf_pca_status.png"), p4, width = 8, height = 6.5, dpi = 300)

  p5 <- ggplot(pca_coords %>% mutate(ST_label = ifelse(ST %in% top_sts, ST, "Other")),
               aes(x = PC1, y = PC2, colour = ST_label)) +
    geom_point(alpha = 0.7, size = 2) +
    labs(
      title = "Exploratory PCA of VF module profiles by sequence type",
      subtitle = "Lineage clustering is a confounding diagnostic for VF-status analyses",
      x = sprintf("PC1 (%.1f%%)", pca_coords$var_PC1[1]),
      y = sprintf("PC2 (%.1f%%)", pca_coords$var_PC2[1]),
      caption = sprintf("Data: %s and script 26 module matrix. ST labels are grouped as Other for sparse lineages.", FILE_VF_READY)
    ) +
    plot_theme_vf(base_size = 11) +
    theme(legend.position = "right")
  ggsave(file.path(DIR_PLOTS_VF, "vf_pca_ST.png"), p5, width = 10, height = 6.5, dpi = 300)
}

if (pcoa_ok && nrow(pcoa_coords) > 0) {
  p6 <- ggplot(pcoa_coords, aes(x = Axis1, y = Axis2, colour = Infection_Status)) +
    geom_point(alpha = 0.7, size = 2) +
    labs(
      title = "Exploratory Jaccard PCoA of VF module profiles by clinical status",
      subtitle = "Binary module-profile ordination; descriptive only",
      x = sprintf("Axis 1 (%.1f%%)", pcoa_coords$var_Axis1[1]),
      y = sprintf("Axis 2 (%.1f%%)", pcoa_coords$var_Axis2[1]),
      caption = sprintf("Data: %s and script 26 module matrix. Residents may contribute repeated isolates; UTI n=%d is small.", FILE_VF_READY, n_uti)
    ) +
    plot_theme_vf(base_size = 11)
  ggsave(file.path(DIR_PLOTS_VF, "vf_pcoa_jaccard_status.png"), p6, width = 8, height = 6.5, dpi = 300)

  p7 <- ggplot(pcoa_coords %>% mutate(ST_label = ifelse(ST %in% top_sts, ST, "Other")),
               aes(x = Axis1, y = Axis2, colour = ST_label)) +
    geom_point(alpha = 0.7, size = 2) +
    labs(
      title = "Exploratory Jaccard PCoA of VF module profiles by sequence type",
      subtitle = "Lineage structure in VF modules should be considered before interpreting status associations",
      x = sprintf("Axis 1 (%.1f%%)", pcoa_coords$var_Axis1[1]),
      y = sprintf("Axis 2 (%.1f%%)", pcoa_coords$var_Axis2[1]),
      caption = sprintf("Data: %s and script 26 module matrix. ST labels are grouped as Other for sparse lineages.", FILE_VF_READY)
    ) +
    plot_theme_vf(base_size = 11) +
    theme(legend.position = "right")
  ggsave(file.path(DIR_PLOTS_VF, "vf_pcoa_jaccard_ST.png"), p7, width = 10, height = 6.5, dpi = 300)
}

msg("✓ 27_vf_score_framework.R complete.")

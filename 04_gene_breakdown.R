#!/usr/bin/env Rscript
# ==============================================================================
# 04_gene_breakdown.R
# ==============================================================================
#
# GOAL:
#   Annotate detected VF genes with functional categories (e.g., Adhesion,
#   Iron Acquisition, Toxins), perform gene-level GLMM tests on focus genes,
#   and produce the gene_map.csv that the VF pipeline (22–25) depends on.
#
# WHY THIS SCRIPT EXISTS:
#   Raw Abricate output (from 02_) gives gene names but no biological context.
#   This script maps each gene to a VF functional category using the VFDB
#   classification, enabling category-level summaries ("How many adhesion
#   genes does this isolate carry?") that are more biologically interpretable
#   than individual gene lists.
#
#   It also runs mixed-effects logistic regression (GLMM) on selected
#   "focus genes" (e.g., fimH, papGII, hlyA, cnf1) to test their association
#   with primary UTI vs Not_UTI status, accounting for within-participant
#   correlation.
#
# INPUTS:
#   - results/vf/vf_hits_all.rds          (from 02_gene_presence_analysis.R)
#   - results/clinical/status_map.csv     (from 00b_classify_episodes.R)
#
# OUTPUTS:
#   - results/vf/annotated_gene_table.csv       (gene-level annotation)
#   - results/vf/gene_map.csv                   (Gene → Category mapping)
#   - results/vf/per_sample_category_counts.csv  (category counts per isolate)
#   - results/vf/diff_focus_genes_UTI_vs_Not_UTI_glmm.csv  (focus gene GLMM)
#   - plots/vf/                                 (category and gene plots)
# ==============================================================================
#
# Usage:
#   Rscript 04_gene_breakdown.R
#
# Biological/Statistical purpose:
#   - Annotates detected genes with functional categories.
#   - Tests whether specific virulence factors (e.g., adhesins, toxins) are
#     associated with UTI vs Not_UTI using GLMMs to account for repeated measures.
# ==============================================================================

# 1. Load Configuration & Libraries
source("00_config.R")
source("R/pipeline_qc_helpers.R")
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(tibble)
  library(ggplot2)
  library(stringr)
  library(purrr)
  library(forcats)
})

# 2. Load Data
# ------------------------------------------------------------------------------
if (!file.exists(FILE_VF_HITS)) stop("Missing ", FILE_VF_HITS)
vf_hits_all <- readRDS(FILE_VF_HITS)

ensure_dir(DIR_PLOTS_VF)

# 3. Gene Mapping & Annotation
# ------------------------------------------------------------------------------
FILE_GENE_MAP <- file.path(DIR_VF, "gene_map.csv")

if (file.exists(FILE_GENE_MAP)) {
  gene_map <- read_csv(FILE_GENE_MAP, show_col_types = FALSE) %>%
    mutate(
      Gene = as.character(Gene),
      Category = coalesce(as.character(Category), "Unassigned"),
      Subcategory = coalesce(as.character(Subcategory), Category)
    )
} else {
  message("↪ Creating starter gene map...")
  genes_seen <- vf_hits_all %>%
    rename(Gene = any_of(c("GENE", "gene", "Gene"))) %>%
    mutate(Gene = trimws(as.character(Gene))) %>%
    distinct(Gene) %>%
    filter(!is.na(Gene), Gene != "")

  # Heuristic categorization
  gene_map <- genes_seen %>%
    mutate(
      Category = case_when(
        str_detect(Gene, regex("^(fim|fml|pil|foc|sfa|pap|afa|dra|cfa)", TRUE)) ~ "Adhesion/Fimbriae",
        str_detect(Gene, regex("^(kps|kfi|neu|ugd|rmpA|caps|wzx|wzy)", TRUE)) ~ "Capsule/Surface",
        str_detect(Gene, regex("^(iut|iuc|iro|irp|fyuA|chu|fep|ent|fec|ybt)", TRUE)) ~ "Iron acquisition",
        str_detect(Gene, regex("^(hly|cnf|sat|vat|cdt|astA|subAB|stx|lt|st)", TRUE)) ~ "Toxins",
        str_detect(Gene, regex("^(omp|iss|ibe|tra|usp|malX)", TRUE)) ~ "Invasion/Evasion",
        str_detect(Gene, regex("^(bla|qnr|aac|aph|aad|erm|cat|tet|sul|dfr|mcr|gyrA|parC)", TRUE)) ~ "AMR",
        TRUE ~ "Unassigned"
      ),
      Subcategory = Category
    )

  write_csv(gene_map, FILE_GENE_MAP)
  message("✓ Starter map written to ", FILE_GENE_MAP)
}

# 4. Standardize Hits Table
# ------------------------------------------------------------------------------
hits <- vf_hits_all
gene_col <- intersect(c("GENE", "gene", "Gene"), names(hits))[1]
if (is.na(gene_col)) stop("No gene column found.")
names(hits)[names(hits) == gene_col] <- "Gene"

# Timepoint normalization
tp_norm <- function(x) {
  tp_chr <- as.character(x)
  is_uricult <- str_detect(tp_chr, regex("uricult", ignore_case = TRUE))
  tp_num <- suppressWarnings(as.integer(str_extract(tp_chr, "\\d+")))
  tp_lab <- case_when(
    is_uricult ~ "Uricult",
    !is.na(tp_num) ~ paste0("T", tp_num),
    TRUE ~ "Unscheduled"
  )
  tibble(tp_lab = tp_lab, tp_num = tp_num)
}

if (!"tp_lab" %in% names(hits)) {
  hits <- bind_cols(hits, tp_norm(hits$Timepoint))
}

hits$Gene <- trimws(hits$Gene)

# 5. Apply Aliases
# ------------------------------------------------------------------------------
alias_rules <- tribble(
  ~pattern, ~Gene,
  "^bla[-_ ]?ctx[-_ ]?m[-_ ]?15(\\.[0-9]+)?$", "blaCTX-M-15",
  "^bla[-_ ]?ctx[-_ ]?m[-_ ]?[0-9]+[a-zA-Z]?$", "blaCTX-M",
  "^bla[-_ ]?oxa[-_ ]?48$", "blaOXA-48",
  "^aac\\(6'\\)[-_ ]?ib[-_ ]?cr[0-9]*$", "aac(6')-Ib-cr",
  "^qnrA[0-9]*$", "qnrA",
  "^qnrB[0-9]*$", "qnrB",
  "^qnrS[0-9]*$", "qnrS",
  "^dfrA[0-9]+.*$", "dfrA",
  "^mcr[-_ ]?([1-9]|10)(\\.[0-9]+)?$", "mcr-\\1"
)

apply_aliases <- function(df, gene_col = "Gene") {
  x <- as.character(df[[gene_col]])
  seen <- rep(FALSE, length(x))
  for (i in seq_len(nrow(alias_rules))) {
    pat <- alias_rules$pattern[i]
    repl <- alias_rules$Gene[i]
    idx <- !seen & grepl(pat, x, ignore.case = TRUE, perl = TRUE)
    if (any(idx, na.rm = TRUE)) {
      x[idx] <- gsub(pat, repl, x[idx], ignore.case = TRUE, perl = TRUE)
      seen[idx] <- TRUE
    }
  }
  df[[gene_col]] <- x
  df
}

hits <- apply_aliases(hits, "Gene")

# Annotate
annotated <- hits %>%
  left_join(gene_map, by = "Gene") %>%
  mutate(
    Category    = coalesce(Category, "Unassigned"),
    Subcategory = coalesce(Subcategory, "Unassigned")
  )

write_csv(annotated, file.path(DIR_VF, "annotated_gene_table.csv"))

# Category Summary
cat_summary <- annotated %>%
  group_by(Participant_id, Timepoint, Category) %>%
  summarise(n_genes = n_distinct(Gene), .groups = "drop") %>%
  pivot_wider(names_from = Category, values_from = n_genes, values_fill = 0)

write_csv(cat_summary, file.path(DIR_VF, "per_sample_category_counts.csv"))

# 6. Nitrate Analysis
# ------------------------------------------------------------------------------
canon <- function(x) str_to_lower(gsub("[^a-z0-9]+", "", x))
nitrate_sys <- list(
  Nar = c("narG", "narH", "narJ", "narI"),
  Nap = c("napF", "napD", "napA", "napG", "napH", "napB", "napC"),
  Nas = c("nasA", "nasB", "nasG", "nasH", "nasC")
)

hits_nc <- hits %>% mutate(Gene_clean = canon(Gene))
nitrate_long <- hits_nc %>%
  filter(Gene_clean %in% canon(unlist(nitrate_sys))) %>%
  mutate(System = case_when(
    Gene_clean %in% canon(nitrate_sys$Nar) ~ "Nar",
    Gene_clean %in% canon(nitrate_sys$Nap) ~ "Nap",
    Gene_clean %in% canon(nitrate_sys$Nas) ~ "Nas"
  )) %>%
  distinct(Participant_id, Timepoint, System)

if (nrow(nitrate_long)) {
  nitrate_mat <- nitrate_long %>%
    mutate(present = 1) %>%
    pivot_wider(names_from = System, values_from = present, values_fill = 0)

  for (col in c("Nar", "Nap", "Nas")) if (!col %in% names(nitrate_mat)) nitrate_mat[[col]] <- 0

  write_csv(nitrate_mat, file.path(DIR_VF, "nitrate_presence_matrix.csv"))

  # Plot
  if (any(colSums(nitrate_mat[, c("Nar", "Nap", "Nas"), drop = FALSE]) > 0)) {
    if (requireNamespace("ComplexUpset", quietly = TRUE)) {
      ups_df <- nitrate_mat %>%
        unite(Sample, Participant_id, Timepoint, sep = "_") %>%
        mutate(across(c(Nar, Nap, Nas), as.logical))
      p <- ComplexUpset::upset(ups_df, intersect = c("Nar", "Nap", "Nas"), name = "Nitrate System")
      ggsave(file.path(DIR_PLOTS_VF, "nitrate_upset.png"), p, width = 6, height = 4)
    } else {
      message("⚠ Skipping Nitrate UpSet plot: Package 'ComplexUpset' not installed.")
    }
  }
}

# 7. Focus Genes Analysis (GLMM)
# ------------------------------------------------------------------------------
FILE_STATUS_MAP <- file.path(DIR_CLINICAL, "status_map.csv")

if (file.exists(FILE_STATUS_MAP)) {
  status_map <- read_csv(FILE_STATUS_MAP, show_col_types = FALSE) %>%
    prefer_primary_uti_status()
  if (!"tp_lab" %in% names(status_map)) {
    status_map <- bind_cols(status_map, tp_norm(status_map$Timepoint))
  }
  status_map <- status_map %>%
    select(any_of(c("Participant_id", "tp_lab", "Infection_Status", "UTI_Status",
                    "UTI_binary", "Not_UTI_subgroup", "Batch"))) %>%
    distinct()
  if (!"Batch" %in% names(status_map)) status_map$Batch <- NA_character_

  # Focus gene list (subset for brevity in example, expand as needed)
  focus_genes <- c("afa/draBC", "hlyA", "cnf1", "iroN", "fyuA", "papA", "papC", "iutA", "kpsM II")

  # Filter hits for focus genes
  hits_focus <- hits %>% filter(Gene %in% focus_genes)

  # [STAT] NOTE: Fisher tests below use isolate-level data → pseudoreplication
  # (Multiple timepoints per participant treated as independent)
  # Results are EXPLORATORY ONLY. Use GLMM results (below) for inference.

  if (nrow(hits_focus) > 0) {
    focus_pa <- hits_focus %>%
      distinct(Participant_id, tp_lab, Gene) %>%
      mutate(present = 1L) %>%
      pivot_wider(names_from = Gene, values_from = present, values_fill = 0L)

    samples_focus <- status_map %>% left_join(focus_pa, by = c("Participant_id", "tp_lab"))

    # GLMM Analysis (UTI vs Not_UTI)
    focus_cols <- intersect(names(samples_focus), focus_genes)
    if (length(focus_cols) > 0) {
      if (!requireNamespace("lme4", quietly = TRUE)) {
        warning("lme4 package not found. Skipping GLMM analysis.")
      } else {
        samples_focus <- samples_focus %>%
          mutate(across(all_of(focus_cols), ~ replace_na(., 0L))) %>%
          filter(UTI_Status %in% c("UTI", "Not_UTI")) %>%
          mutate(
            Outcome = as.integer(UTI_binary),
            Participant_id = as.factor(Participant_id),
            # [STAT] Timepoint as numeric for trend test (if defined)
            tp_num = suppressWarnings(as.integer(str_extract(tp_lab, "\\d+"))),
            # [STAT] Batch as factor (if exists in status_map)
            Batch = if ("Batch" %in% names(.)) as.factor(Batch) else NA_character_
          )

        # [STAT] Enhanced GLMM with proper adjustment for confounders
        # [LIMITATION] Clinical covariates (Age, Catheter, Diabetes, Antibiotics) are currently
        # missing from the dataset. Their absence is a potential source of unmeasured confounding.
        # Future versions should include them in the fixed effects if available.
        run_glmm <- function(gene_col) {
          # Helper to format result
          format_res <- function(mod, type, converged) {
            coefs <- summary(mod)$coefficients
            if (gene_col %in% rownames(coefs)) {
              gene_prev <- mean(samples_focus[[gene_col]], na.rm = TRUE)
              tibble(
                Gene = gene_col,
                OR = exp(coefs[gene_col, "Estimate"]),
                CI_lower = exp(coefs[gene_col, "Estimate"] - 1.96 * coefs[gene_col, "Std. Error"]),
                CI_upper = exp(coefs[gene_col, "Estimate"] + 1.96 * coefs[gene_col, "Std. Error"]),
                p = coefs[gene_col, "Pr(>|z|)"],
                Converged = converged,
                Adjusted_for = paste(covariate_terms, collapse = "+"),
                Prevalence = gene_prev,
                Power_Flag = ifelse(gene_prev >= 0.10, "Adequate", "Underpowered (<10%)"),
                Role = ifelse(gene_prev >= 0.10, "Inferential-core", "Exploratory"),
                Model_Type = type
              )
            } else {
              tibble(Gene = gene_col, OR = NA, p = NA, Converged = FALSE, Model_Type = type)
            }
          }

          # Build formula based on available covariates
          # Priority: adjust for Timepoint (temporal trend) and Batch (batch effects)
          has_tp <- sum(!is.na(samples_focus$tp_num)) > 10
          has_batch <- "Batch" %in% names(samples_focus) && n_distinct(samples_focus$Batch, na.rm = TRUE) > 1

          covariate_terms <- c()
          if (has_tp) covariate_terms <- c(covariate_terms, "tp_num")
          if (has_batch) covariate_terms <- c(covariate_terms, "Batch")

          fixed_part <- paste(c(gene_col, covariate_terms), collapse = " + ")

          tryCatch(
            {
              # 1. Try GLMM
              fmla_glmm <- as.formula(paste("Outcome ~", fixed_part, "+ (1|Participant_id)"))
              m <- lme4::glmer(fmla_glmm,
                data = samples_focus, family = binomial,
                control = lme4::glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1e5))
              )

              is_singular <- isSingular(m)
              if (is_singular) {
                format_res(m, "GLMM (Singular)", TRUE)
              } else {
                format_res(m, "GLMM", TRUE)
              }
            },
            error = function(e) {
              # 2. Fallback to GLM
              tryCatch(
                {
                  fmla_glm <- as.formula(paste("Outcome ~", fixed_part))
                  m_glm <- glm(fmla_glm, data = samples_focus, family = binomial)
                  format_res(m_glm, "GLM_Fallback", TRUE)
                },
                error = function(e2) {
                  tibble(Gene = gene_col, OR = NA, p = NA, Converged = FALSE, Model_Type = "Failed")
                }
              )
            }
          )
        }

        gene_enrichment <- purrr::map_dfr(focus_cols, run_glmm) %>%
          mutate(p_adj = p.adjust(p, method = "BH")) %>%
          arrange(p_adj)

        write_csv(gene_enrichment, file.path(DIR_VF, "diff_focus_genes_UTI_vs_Not_UTI_glmm.csv"))
        message("✓ GLMM results written to diff_focus_genes_UTI_vs_Not_UTI_glmm.csv")
      }
    }
  }
}

message("✓ Gene breakdown complete.")

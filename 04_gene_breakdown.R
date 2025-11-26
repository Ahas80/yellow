#!/usr/bin/env Rscript
# ==============================================================================
# 04_gene_breakdown.R
# ------------------------------------------------------------------------------
# Role: [Inferential-core] - Detailed gene analysis, annotation, and focus-gene statistics.
#
# Inputs:
#   - results/vf/vf_hits_all.rds
#   - results/clinical/status_map.csv
#
# Outputs:
#   - results/vf/annotated_gene_table.csv
#   - results/vf/per_sample_category_counts.csv
#   - results/vf/nitrate_presence_matrix.csv
#   - results/vf/diff_focus_genes_UTI_vs_ASB_glmm.csv
#   - plots/vf/
#
# Usage:
#   Rscript 04_gene_breakdown.R
#
# Biological/Statistical purpose:
#   - Annotates detected genes with functional categories.
#   - Tests whether specific virulence factors (e.g., adhesins, toxins) are
#     associated with UTI vs ASB using GLMMs to account for repeated measures.
# ==============================================================================

# 1. Load Configuration & Libraries
source("00_config.R")
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(tibble)
  library(ggplot2)
  library(ComplexUpset)
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
    ups_df <- nitrate_mat %>%
      unite(Sample, Participant_id, Timepoint, sep = "_") %>%
      mutate(across(c(Nar, Nap, Nas), as.logical))
    p <- ComplexUpset::upset(ups_df, intersect = c("Nar", "Nap", "Nas"), name = "Nitrate System")
    ggsave(file.path(DIR_PLOTS_VF, "nitrate_upset.png"), p, width = 6, height = 4)
  }
}

# 7. Focus Genes Analysis (GLMM)
# ------------------------------------------------------------------------------
FILE_STATUS_MAP <- file.path(DIR_CLINICAL, "status_map.csv")

if (file.exists(FILE_STATUS_MAP)) {
  status_map <- read_csv(FILE_STATUS_MAP, show_col_types = FALSE)
  if (!"tp_lab" %in% names(status_map)) {
    status_map <- bind_cols(status_map, tp_norm(status_map$Timepoint))
  }
  status_map <- status_map %>%
    select(Participant_id, tp_lab, Infection_Status) %>%
    distinct()

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

    # GLMM Analysis (UTI vs ASB)
    focus_cols <- intersect(names(samples_focus), focus_genes)
    if (length(focus_cols) > 0) {
      if (!requireNamespace("lme4", quietly = TRUE)) {
        warning("lme4 package not found. Skipping GLMM analysis.")
      } else {
        samples_focus <- samples_focus %>%
          mutate(across(all_of(focus_cols), ~ replace_na(., 0L))) %>%
          filter(Infection_Status %in% c("UTI", "ASB")) %>%
          mutate(
            Outcome = ifelse(Infection_Status == "UTI", 1, 0),
            Participant_id = as.factor(Participant_id),
            # [STAT] Timepoint as numeric for trend test (if defined)
            tp_num = suppressWarnings(as.integer(str_extract(tp_lab, "\\d+"))),
            # [STAT] Batch as factor (if exists in status_map)
            Batch = if ("Batch" %in% names(.)) as.factor(Batch) else NA_character_
          )

        # [STAT] Enhanced GLMM with proper adjustment for confounders
        run_glmm <- function(gene_col) {
          # Build formula based on available covariates
          # Priority: adjust for Timepoint (temporal trend) and Batch (batch effects)
          has_tp <- sum(!is.na(samples_focus$tp_num)) > 10
          has_batch <- "Batch" %in% names(samples_focus) && n_distinct(samples_focus$Batch, na.rm = TRUE) > 1

          covariate_terms <- c()
          if (has_tp) covariate_terms <- c(covariate_terms, "tp_num")
          if (has_batch) covariate_terms <- c(covariate_terms, "Batch")

          fixed_part <- paste(c(gene_col, covariate_terms), collapse = " + ")
          fmla <- as.formula(paste("Outcome ~", fixed_part, "+ (1|Participant_id)"))

          tryCatch(
            {
              # Use glmer for mixed effects
              m <- lme4::glmer(fmla,
                data = samples_focus, family = binomial,
                control = lme4::glmerControl(optimizer = "bobyqa")
              )
              coefs <- summary(m)$coefficients
              if (gene_col %in% rownames(coefs)) {
                # Calculate gene prevalence for power flag
                gene_prev <- mean(samples_focus[[gene_col]], na.rm = TRUE)
                tibble(
                  Gene = gene_col,
                  OR = exp(coefs[gene_col, "Estimate"]),
                  CI_lower = exp(coefs[gene_col, "Estimate"] - 1.96 * coefs[gene_col, "Std. Error"]),
                  CI_upper = exp(coefs[gene_col, "Estimate"] + 1.96 * coefs[gene_col, "Std. Error"]),
                  p = coefs[gene_col, "Pr(>|z|)"],
                  Converged = TRUE,
                  Adjusted_for = paste(covariate_terms, collapse = "+"),
                  Prevalence = gene_prev,
                  Power_Flag = ifelse(gene_prev >= 0.10, "Adequate", "Underpowered (<10%)"),
                  Role = ifelse(gene_prev >= 0.10, "Inferential-core", "Exploratory")
                )
              } else {
                tibble(Gene = gene_col, OR = NA, p = NA, Converged = FALSE)
              }
            },
            error = function(e) {
              # Fallback to simple glm if singular fit or other error (often happens with rare genes)
              tryCatch(
                {
                  m_glm <- glm(as.formula(paste("Outcome ~", gene_col)), data = samples_focus, family = binomial)
                  coefs <- summary(m_glm)$coefficients
                  if (gene_col %in% rownames(coefs)) {
                    tibble(
                      Gene = gene_col,
                      OR = exp(coefs[gene_col, "Estimate"]),
                      p = coefs[gene_col, "Pr(>|z|)"],
                      Converged = "GLM_Fallback"
                    )
                  } else {
                    tibble(Gene = gene_col, OR = NA, p = NA, Converged = FALSE)
                  }
                },
                error = function(e2) {
                  tibble(Gene = gene_col, OR = NA, p = NA, Converged = FALSE)
                }
              )
            }
          )
        }

        gene_enrichment <- purrr::map_dfr(focus_cols, run_glmm) %>%
          mutate(p_adj = p.adjust(p, method = "BH")) %>%
          arrange(p_adj)

        write_csv(gene_enrichment, file.path(DIR_VF, "diff_focus_genes_UTI_vs_ASB_glmm.csv"))
        message("✓ GLMM results written to diff_focus_genes_UTI_vs_ASB_glmm.csv")
      }
    }
  }
}

message("✓ Gene breakdown complete.")

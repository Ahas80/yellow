#!/usr/bin/env Rscript
# ==============================================================================
# 14_genotype_phenotype_model.R
# ------------------------------------------------------------------------------
# Role: [Modelling] - Test genotype-phenotype associations (GWAS).
#
# Inputs:
#   - results/clinical/status_map.csv
#   - results/vf/vf_hits_all.csv
#   - results/plasmids/plasmidfinder_presence_absence.csv
#   - results/mlst/mlst_with_meta.csv
#   - assembly_metadata.csv
#
# Outputs:
#   - results/models/gwas_univariable_stats.csv
#   - results/models/gwas_multivariable_glmm.csv
#   - results/models/plots/volcano_plot_UTI_vs_ASB.png
#   - results/models/plots/forest_plot_top_hits.png
#   - results/models/plots/heatmap_top_discriminators.png
#
# Usage:
#   Rscript 14_genotype_phenotype_model.R
#
# Biological/Statistical purpose:
#   - Identifies bacterial genomic features (VFs, plasmids, lineages) associated
#     with symptomatic UTI vs ASB, adjusting for repeated measures (GLMM).
# ==============================================================================

Sys.setlocale("LC_ALL", "en_US.UTF-8")
options(bitmapType = "cairo")

# 1. Load Configuration & Libraries
# ------------------------------------------------------------------------------
source("00_config.R")

suppressPackageStartupMessages({
    library(tidyverse)
    library(optparse)
    library(lme4)
    library(broom.mixed)
    library(ggrepel)
    library(pheatmap)
    library(future)
    library(future.apply)
})

# 2. CLI Options
# ------------------------------------------------------------------------------
option_list <- list(
    make_option(c("-t", "--threads"),
        type = "integer", default = 10,
        help = "Number of parallel threads [default: %default]"
    ),
    make_option("--min_prev",
        type = "double", default = 0.05,
        help = "Minimum feature prevalence to retain [default: %default]"
    ),
    make_option("--max_prev",
        type = "double", default = 0.95,
        help = "Maximum feature prevalence to retain [default: %default]"
    ),
    make_option("--run_rf",
        action = "store_true", default = FALSE,
        help = "Run Random Forest variable importance analysis"
    ),
    make_option("--simple-glm",
        action = "store_true", default = FALSE,
        help = "Use standard GLM instead of GLMM (no random effects). Use when you have severe class imbalance or singular fit warnings."
    ),
    make_option("--fdr_thresh",
        type = "double", default = 0.05,
        help = "FDR threshold for reporting top hits [default: %default]"
    ),
    make_option("--outdir",
        type = "character", default = "results/models",
        help = "Output directory [default: %default]"
    )
)

opt <- parse_args(OptionParser(option_list = option_list))
CORES_USE <- min(opt$threads, parallel::detectCores() - 1)

msg <- function(fmt, ...) {
    cat(sprintf("[%s] %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), sprintf(fmt, ...)))
}

msg("Starting genotype-phenotype association analysis")
msg(
    "Threads: %d | Min/Max prevalence: %.2f-%.2f | FDR threshold: %.3f | Model: %s",
    CORES_USE, opt$min_prev, opt$max_prev, opt$fdr_thresh,
    ifelse(opt$`simple-glm`, "Standard GLM", "Mixed-effects GLMM")
)

# 3. Setup Directories
# ------------------------------------------------------------------------------
DIR_OUT <- opt$outdir
DIR_PLOTS <- file.path(DIR_OUT, "plots")
dir.create(DIR_OUT, recursive = TRUE, showWarnings = FALSE)
dir.create(DIR_PLOTS, recursive = TRUE, showWarnings = FALSE)

# 4. Load Data
# ------------------------------------------------------------------------------
msg("Loading input datasets...")

# Clinical status
FILE_STATUS <- file.path(DIR_CLINICAL, "status_map.csv")
if (!file.exists(FILE_STATUS)) stop("Missing ", FILE_STATUS)
status <- read_csv(FILE_STATUS, show_col_types = FALSE) %>%
    filter(Infection_Status %in% c("UTI", "ASB")) %>%
    mutate(Outcome = as.integer(Infection_Status == "UTI"))

msg(
    "  Loaded %d samples with UTI/ASB status (%d UTI, %d ASB)",
    nrow(status), sum(status$Outcome == 1), sum(status$Outcome == 0)
)

# VF/AMR matrix - Use hits file and create presence/absence matrix
FILE_VF <- file.path(DIR_VF, "vf_hits_all.csv")
if (!file.exists(FILE_VF)) stop("Missing ", FILE_VF)
vf_hits <- read_csv(FILE_VF, show_col_types = FALSE)

# Create presence/absence: 1 if GENE present for that Participant/Timepoint
# Need to extract participant and timepoint from the hits
vf_raw <- vf_hits %>%
    filter(!is.na(Participant_id), !is.na(GENE)) %>%
    distinct(Participant_id, Timepoint, GENE) %>%
    mutate(
        present = 1,
        # Normalize timepoint to match status_map format
        Timepoint = case_when(
            str_detect(Timepoint, "(?i)uricult") ~ "Uricult",
            str_detect(Timepoint, "^T\\d+$") ~ Timepoint,
            str_detect(Timepoint, "\\d+") ~ paste0("T", str_extract(Timepoint, "\\d+")),
            TRUE ~ Timepoint
        )
    ) %>%
    pivot_wider(
        id_cols = c(Participant_id, Timepoint),
        names_from = GENE,
        values_from = present,
        values_fill = 0
    )

msg(
    "  Loaded VF data: %d participant-timepoints, %d genes",
    nrow(vf_raw), ncol(vf_raw) - 2
)

# Plasmid replicons
FILE_PLASMID <- file.path(DIR_PLASMIDS, "plasmidfinder_presence_absence.csv")
plasmid_raw <- NULL
if (file.exists(FILE_PLASMID)) {
    plasmid_raw <- read_csv(FILE_PLASMID, show_col_types = FALSE)
    msg(
        "  Loaded plasmid data: %d isolates, %d replicons",
        nrow(plasmid_raw), ncol(plasmid_raw) - 1
    )
}

# MLST
FILE_MLST <- file.path(DIR_MLST, "mlst_with_meta.csv")
mlst <- NULL
if (file.exists(FILE_MLST)) {
    mlst <- read_csv(FILE_MLST, show_col_types = FALSE) %>%
        select(Isolate_ID, ST) %>%
        distinct()
    msg("  Loaded MLST data: %d isolates", nrow(mlst))
}

# Assembly metadata (to link Isolate_ID with Participant_id/Timepoint)
# FILE_METADATA is defined in 00_config.R
if (!file.exists(FILE_METADATA)) stop("Missing ", FILE_METADATA)
metadata <- read_csv(FILE_METADATA, show_col_types = FALSE) %>%
    select(Isolate_ID = file_name, Participant_id, Timepoint) %>%
    distinct()

# ============================= DATA HARMONIZATION =============================
msg("Harmonizing datasets...")

# VF data is already in Participant_id/Timepoint format from pivot above
vf <- vf_raw %>%
    # Ensure character types for join keys
    mutate(
        Participant_id = as.character(Participant_id),
        Timepoint = as.character(Timepoint)
    )

msg("  VF data ready: %d participant-timepoints", nrow(vf))

# Ensure status has character types for join
status <- status %>%
    mutate(
        Participant_id = as.character(Participant_id),
        Timepoint = as.character(Timepoint)
    )

# Join with status
data_vf <- status %>%
    inner_join(vf, by = c("Participant_id", "Timepoint"))

msg("  After VF join: %d samples", nrow(data_vf))

# Plasmid data: Need to link Isolate_ID -> Participant_id/Timepoint
data_plasmid <- NULL
if (!is.null(plasmid_raw)) {
    # Link plasmid Isolate_ID to Participant_id/Timepoint via metadata
    plasmid_linked <- plasmid_raw %>%
        inner_join(metadata, by = "Isolate_ID") %>%
        select(-Isolate_ID)

    # Normalize Timepoint format to match status map
    # Extract timepoint from format like "T0", "T1", "Uricult"
    plasmid_linked <- plasmid_linked %>%
        mutate(Timepoint = case_when(
            str_detect(Timepoint, "(?i)uricult") ~ "Uricult",
            str_detect(Timepoint, "T\\d+") ~ Timepoint,
            TRUE ~ Timepoint
        ))

    # For each Participant_id x Timepoint, take max across isolates (multiple assemblies)
    rep_cols <- setdiff(names(plasmid_linked), c("Participant_id", "Timepoint"))
    plasmid_agg <- plasmid_linked %>%
        group_by(Participant_id, Timepoint) %>%
        summarise(across(all_of(rep_cols), ~ max(as.numeric(.), na.rm = TRUE)), .groups = "drop") %>%
        mutate(across(all_of(rep_cols), ~ ifelse(is.finite(.), ., 0)))

    # Join with status
    data_plasmid <- status %>%
        inner_join(plasmid_agg, by = c("Participant_id", "Timepoint"))

    msg("  After plasmid join: %d samples", nrow(data_plasmid))
}

# Combine VF + Plasmid features
if (!is.null(data_plasmid)) {
    # Merge on status columns
    data_merged <- data_vf %>%
        left_join(
            data_plasmid %>% select(-Infection_Status, -Outcome),
            by = c("Participant_id", "Timepoint")
        )
} else {
    data_merged <- data_vf
}

msg(
    "  Final merged dataset: %d samples, %d features",
    nrow(data_merged), ncol(data_merged) - 4
)

# Feature filtering: Remove low/high prevalence
feature_cols <- setdiff(names(data_merged), c("Participant_id", "Timepoint", "Infection_Status", "Outcome"))
prevalences <- data_merged %>%
    select(all_of(feature_cols)) %>%
    summarise(across(everything(), ~ mean(. > 0, na.rm = TRUE))) %>%
    pivot_longer(everything(), names_to = "feature", values_to = "prevalence")

features_keep <- prevalences %>%
    filter(prevalence >= opt$min_prev, prevalence <= opt$max_prev) %>%
    pull(feature)

msg(
    "  Features retained after prevalence filter: %d (from %d total)",
    length(features_keep), length(feature_cols)
)

data_final <- data_merged %>%
    select(Participant_id, Timepoint, Infection_Status, Outcome, all_of(features_keep))

# ========================= UNIVARIABLE ASSOCIATION ============================
msg("Running univariable association tests...")

# [STAT] NOTE: Fisher tests use episode-level data but ignore participant clustering
# This violates independence → inflated Type I error (false positives)
# Results are EXPLORATORY screening only. Rely on GLMM for inference.

run_univar <- function(feature, data) {
    tryCatch(
        {
            # Create 2x2 contingency table
            tab <- table(data[[feature]] > 0, data$Outcome)

            # Fisher's exact test
            test <- fisher.test(tab)

            # Calculate prevalence in each group
            prev_asb <- mean(data[[feature]][data$Outcome == 0] > 0, na.rm = TRUE)
            prev_uti <- mean(data[[feature]][data$Outcome == 1] > 0, na.rm = TRUE)

            tibble(
                feature = feature,
                prev_ASB = prev_asb,
                prev_UTI = prev_uti,
                OR = test$estimate,
                CI_lower = test$conf.int[1],
                CI_upper = test$conf.int[2],
                p_value = test$p.value,
                test_type = "Fisher"
            )
        },
        error = function(e) {
            tibble(
                feature = feature,
                prev_ASB = NA_real_,
                prev_UTI = NA_real_,
                OR = NA_real_,
                CI_lower = NA_real_,
                CI_upper = NA_real_,
                p_value = NA_real_,
                test_type = "Failed"
            )
        }
    )
}

# Parallel execution
plan(multisession, workers = CORES_USE)
univar_results <- future.apply::future_lapply(features_keep, run_univar,
    data = data_final,
    future.seed = TRUE
) %>%
    bind_rows()
plan(sequential)

# FDR correction
univar_results <- univar_results %>%
    filter(!is.na(p_value)) %>%
    mutate(
        FDR = p.adjust(p_value, method = "BH"),
        log2_OR = log2(OR),
        neg_log10_p = -log10(p_value)
    ) %>%
    arrange(p_value)

msg(
    "  Completed %d univariable tests (%d significant at FDR < %.3f)",
    nrow(univar_results), sum(univar_results$FDR < opt$fdr_thresh, na.rm = TRUE), opt$fdr_thresh
)

# Save univariable results
write_csv(univar_results, file.path(DIR_OUT, "gwas_univariable_stats.csv"))

# ======================== MULTIVARIABLE MODELLING =================================
if (opt$`simple-glm`) {
    msg("Running multivariable standard GLM analyses (no random effects)...")
} else {
    msg("Running multivariable GLMM analyses (with participant random intercepts)...")
}

# Select features for GLMM: those with p < 0.1 in univariable (or top 50)
glmm_features <- univar_results %>%
    filter(p_value < 0.1 | row_number() <= 50) %>%
    pull(feature)

msg("  Selected %d features for GLMM", length(glmm_features))

# Simple GLM (no random effects)
run_simple_glm <- function(feature, data) {
    tryCatch(
        {
            # Build formula
            # Include Timepoint if more than one level exists
            if (n_distinct(data$Timepoint) > 1) {
                fml <- as.formula(sprintf("Outcome ~ `%s` + Timepoint", feature))
            } else {
                fml <- as.formula(sprintf("Outcome ~ `%s`", feature))
            }

            # Fit model
            mod <- glm(fml, data = data, family = binomial(link = "logit"))

            # Extract coefficients
            tidy_mod <- broom::tidy(mod, conf.int = TRUE) %>%
                filter(term == paste0("`", feature, "`") | term == feature) %>%
                mutate(
                    feature = feature,
                    OR = exp(estimate),
                    OR_lower = exp(conf.low),
                    OR_upper = exp(conf.high),
                    converged = TRUE,
                    model_type = "GLM"
                ) %>%
                select(feature, estimate, std.error, OR, OR_lower, OR_upper, p.value, converged, model_type)

            tidy_mod
        },
        error = function(e) {
            tibble(
                feature = feature,
                estimate = NA_real_,
                std.error = NA_real_,
                OR = NA_real_,
                OR_lower = NA_real_,
                OR_upper = NA_real_,
                p.value = NA_real_,
                converged = FALSE,
                model_type = "GLM"
            )
        }
    )
}

# Mixed-effects GLMM (with participant random intercepts)
run_glmm <- function(feature, data) {
    tryCatch(
        {
            # Build formula
            # [STAT] Adjust for Timepoint (temporal trend) and Batch (batch effects) if available
            has_batch <- "Batch" %in% names(data) && n_distinct(data$Batch, na.rm = TRUE) > 1

            covariate_terms <- c()
            if (n_distinct(data$Timepoint) > 1) covariate_terms <- c(covariate_terms, "Timepoint")
            if (has_batch) covariate_terms <- c(covariate_terms, "Batch")

            fixed_part <- paste(c(paste0("`", feature, "`"), covariate_terms), collapse = " + ")
            fml <- as.formula(sprintf("%s + (1 | Participant_id)", fixed_part))
            fml <- reformulate(
                termlabels = c(paste0("`", feature, "`"), covariate_terms, "(1 | Participant_id)"),
                response = "Outcome"
            )

            # Fit model
            mod <- glmer(fml,
                data = data, family = binomial(link = "logit"),
                control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1e5))
            )

            # Check for singular fit
            is_singular <- isSingular(mod)

            # Extract fixed effects
            # We only want the feature coefficient, not Timepoint or Intercept
            tidy_mod <- tidy(mod, effects = "fixed", conf.int = TRUE) %>%
                filter(term == paste0("`", feature, "`") | term == feature) %>%
                mutate(
                    feature = feature,
                    OR = exp(estimate),
                    OR_lower = exp(conf.low),
                    OR_upper = exp(conf.high),
                    converged = !is_singular,
                    model_type = "GLMM"
                ) %>%
                select(feature, estimate, std.error, OR, OR_lower, OR_upper, p.value, converged, model_type)

            tidy_mod
        },
        error = function(e) {
            tibble(
                feature = feature,
                estimate = NA_real_,
                std.error = NA_real_,
                OR = NA_real_,
                OR_lower = NA_real_,
                OR_upper = NA_real_,
                p.value = NA_real_,
                converged = FALSE,
                model_type = "GLMM"
            )
        }
    )
}

# Parallel execution - choose model type
plan(multisession, workers = CORES_USE)
if (opt$`simple-glm`) {
    glmm_results <- future.apply::future_lapply(glmm_features, run_simple_glm,
        data = data_final,
        future.seed = TRUE
    ) %>%
        bind_rows()
} else {
    glmm_results <- future.apply::future_lapply(glmm_features, run_glmm,
        data = data_final,
        future.seed = TRUE
    ) %>%
        bind_rows()
}
plan(sequential)

# FDR correction
glmm_results <- glmm_results %>%
    filter(!is.na(p.value)) %>%
    mutate(FDR = p.adjust(p.value, method = "BH")) %>%
    arrange(p.value)

if (opt$`simple-glm`) {
    msg(
        "  Completed %d GLM models (%d converged, %d significant at FDR < %.3f)",
        length(glmm_features),
        sum(glmm_results$converged, na.rm = TRUE),
        sum(glmm_results$FDR < opt$fdr_thresh, na.rm = TRUE),
        opt$fdr_thresh
    )
} else {
    msg(
        "  Completed %d GLMM models (%d converged, %d significant at FDR < %.3f)",
        length(glmm_features),
        sum(glmm_results$converged, na.rm = TRUE),
        sum(glmm_results$FDR < opt$fdr_thresh, na.rm = TRUE),
        opt$fdr_thresh
    )
}

# Save GLMM results
msg("Generating plots...")

# Volcano plot
top_hits <- univar_results %>%
    filter(FDR < 0.1) %>%
    slice_min(p_value, n = 20)

p_volcano <- ggplot(univar_results, aes(log2_OR, neg_log10_p)) +
    geom_point(aes(color = FDR < opt$fdr_thresh), alpha = 0.6, size = 2) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
    geom_vline(xintercept = 0, linetype = "solid", color = "grey50") +
    geom_text_repel(
        data = top_hits,
        aes(label = feature),
        size = 3,
        max.overlaps = 15,
        box.padding = 0.5
    ) +
    scale_color_manual(values = c("TRUE" = "red", "FALSE" = "grey40")) +
    labs(
        title = "Genotype-Phenotype Association: UTI vs ASB",
        subtitle = sprintf("%d features tested (Fisher's exact)", nrow(univar_results)),
        x = "Log2 Odds Ratio",
        y = "-Log10 P-value",
        color = sprintf("FDR < %.2f", opt$fdr_thresh)
    ) +
    theme_minimal(base_size = 11) +
    theme(legend.position = "top")

ggsave(file.path(DIR_PLOTS, "volcano_plot_UTI_vs_ASB.png"),
    p_volcano,
    width = 10, height = 8, dpi = 300
)

# Forest plot (top GLMM hits)
top_glmm <- glmm_results %>%
    filter(converged, FDR < 0.2) %>%
    slice_min(FDR, n = 20)

if (nrow(top_glmm) > 0) {
    p_forest <- ggplot(top_glmm, aes(x = OR, y = reorder(feature, OR))) +
        geom_vline(xintercept = 1, linetype = "dashed", color = "grey50") +
        geom_point(size = 3, color = "steelblue") +
        geom_errorbarh(aes(xmin = OR_lower, xmax = OR_upper), height = 0.2) +
        geom_text(aes(label = sprintf("%.2f", OR)), hjust = -0.5, size = 3) +
        scale_x_log10() +
        labs(
            title = "Top GLMM Associations (UTI vs ASB)",
            subtitle = sprintf("%d features with FDR < 0.2", nrow(top_glmm)),
            x = "Odds Ratio (95% CI)",
            y = NULL
        ) +
        theme_minimal(base_size = 11)

    ggsave(file.path(DIR_PLOTS, "forest_plot_top_hits.png"),
        p_forest,
        width = 8, height = max(6, nrow(top_glmm) * 0.3), dpi = 300
    )
}

# Heatmap of top discriminators
top_features <- glmm_results %>%
    filter(converged, !is.na(FDR)) %>%
    slice_min(FDR, n = 30) %>%
    pull(feature)

if (length(top_features) > 1) {
    heatmap_data <- data_final %>%
        select(Participant_id, Infection_Status, all_of(top_features)) %>%
        arrange(Infection_Status, Participant_id) %>%
        column_to_rownames("Participant_id") %>%
        select(-Infection_Status)

    annotation <- data_final %>%
        select(Participant_id, Infection_Status) %>%
        distinct() %>%
        column_to_rownames("Participant_id")

    pheatmap(
        t(as.matrix(heatmap_data)),
        annotation_col = annotation,
        cluster_rows = TRUE,
        cluster_cols = FALSE,
        show_colnames = FALSE,
        color = c("white", "darkblue"),
        border_color = NA,
        filename = file.path(DIR_PLOTS, "heatmap_top_discriminators.png"),
        width = 10,
        height = max(6, length(top_features) * 0.2)
    )
}

# ======================== OPTIONAL: RANDOM FOREST =============================
if (opt$run_rf && requireNamespace("randomForest", quietly = TRUE)) {
    msg("Running Random Forest variable importance analysis...")

    rf_data <- data_final %>%
        select(Outcome, all_of(features_keep)) %>%
        mutate(Outcome = factor(Outcome, levels = c(0, 1), labels = c("ASB", "UTI"))) %>%
        na.omit()

    set.seed(42)
    rf_model <- randomForest::randomForest(
        Outcome ~ .,
        data = rf_data,
        importance = TRUE,
        ntree = 500
    )

    rf_importance <- randomForest::importance(rf_model) %>%
        as.data.frame() %>%
        rownames_to_column("feature") %>%
        arrange(desc(MeanDecreaseGini))

    write_csv(rf_importance, file.path(DIR_OUT, "rf_variable_importance.csv"))
    msg("  Random Forest complete. OOB error: %.2f%%", rf_model$err.rate[500, "OOB"] * 100)
}

# ============================= SUMMARY ========================================
msg("Analysis complete!")
msg("Output directory: %s", normalizePath(DIR_OUT))
msg("")
msg("Summary:")
msg(
    "  Total samples: %d (%d UTI, %d ASB)",
    nrow(data_final),
    sum(data_final$Outcome == 1),
    sum(data_final$Outcome == 0)
)
msg("  Features tested: %d", length(features_keep))
msg(
    "  Univariable significant (FDR < %.2f): %d",
    opt$fdr_thresh, sum(univar_results$FDR < opt$fdr_thresh, na.rm = TRUE)
)
msg(
    "  GLMM significant (FDR < %.2f): %d",
    opt$fdr_thresh, sum(glmm_results$FDR < opt$fdr_thresh, na.rm = TRUE)
)
msg("")

# Print top 10 hits
msg("Top 10 associations (GLMM):")
top10 <- glmm_results %>%
    filter(converged) %>%
    slice_min(FDR, n = 10)

if (nrow(top10) > 0) {
    for (i in 1:nrow(top10)) {
        msg(
            "  %2d. %-20s OR=%.2f (%.2f-%.2f), p=%.2e, FDR=%.2e",
            i,
            top10$feature[i],
            top10$OR[i],
            top10$OR_lower[i],
            top10$OR_upper[i],
            top10$p.value[i],
            top10$FDR[i]
        )
    }
} else {
    msg("  (No significant associations found)")
}

msg("Done.")

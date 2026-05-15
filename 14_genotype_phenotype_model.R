#!/usr/bin/env Rscript
# ==============================================================================
# 14_genotype_phenotype_model.R
# ==============================================================================
#
# GOAL:
#   Perform genome-wide association testing (GWAS-style) to identify bacterial
#   genomic features (VF genes, plasmid replicons, ST lineage) associated with
#   symptomatic UTI vs asymptomatic bacteriuria (ASB).
#
# WHY THIS SCRIPT EXISTS:
#   Scripts 03 and 23 provide descriptive and exploratory Fisher-based
#   comparisons, but these VIOLATE the independence assumption because
#   participants contribute multiple episodes.  This script is the formal
#   inferential engine: it uses Generalised Linear Mixed Models (GLMM) with
#   a random intercept for Participant_id — (1|Participant_id) — to properly
#   account for within-participant correlation.
#
#   This is the ONLY script in the pipeline that provides statistically
#   valid p-values for VF–status associations.  All other enrichment tests
#   (in 23_, 04_) are labelled EXPLORATORY.
#
# STATISTICAL APPROACH:
#   For each genomic feature with prevalence between 5–95%:
#     1. Univariable Fisher exact test (screening)
#     2. GLMM: Outcome ~ Feature + (1|Participant_id),
#        family = binomial, where Outcome = 1 for UTI, 0 for ASB
#     3. Benjamini-Hochberg FDR correction across all features
#   Falls back to standard GLM if GLMM fails to converge (singular fit).
#
# INPUTS:
#   - results/clinical/status_map.csv
#   - results/vf/vf_hits_all.rds
#   - results/plasmids/plasmidfinder_presence_absence.csv
#   - results/mlst/mlst_with_meta.csv
#   - assembly_metadata.csv
#
# OUTPUTS:
#   - results/models/gwas_univariable_stats.csv   (Fisher + prevalence)
#   - results/models/gwas_multivariable_glmm.csv  (GLMM coefficients, FDR)
#   - results/models/plots/volcano_plot_UTI_vs_ASB.png
#   - results/models/plots/forest_plot_top_hits.png
#   - results/models/plots/heatmap_top_discriminators.png
#
# RELATIONSHIP TO OTHER VF SCRIPTS:
#   → 23_vf_cross_sectional.R runs exploratory Fisher tests (screening only)
#   → 25_vf_lineage_vf_interaction.R checks whether ST should be added as
#     a covariate in this GLMM
#   → This script is the definitive statistical test for VF–UTI associations
# ==============================================================================

invisible(Sys.setlocale("LC_ALL", "en_US.UTF-8"))
if (identical(Sys.info()[["sysname"]], "Darwin")) {
    options(bitmapType = "quartz")
}

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
    library(broom)
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

normalize_analysis_timepoint <- function(x) {
    normalise_timepoint_preserve_events(x)
}

save_png_device <- function(filename, width, height, dpi = 300, draw_fun) {
    dir.create(dirname(filename), recursive = TRUE, showWarnings = FALSE)
    opened <- FALSE
    tryCatch(
        {
            if (requireNamespace("ragg", quietly = TRUE)) {
                ragg::agg_png(filename, width = width, height = height, units = "in", res = dpi)
            } else {
                grDevices::png(filename, width = width, height = height, units = "in", res = dpi)
            }
            opened <- TRUE
            draw_fun()
        },
        finally = {
            if (opened && grDevices::dev.cur() > 1) {
                grDevices::dev.off()
            }
        }
    )
}

plot_theme_vf_model <- function(base_size = 11) {
    theme_bw(base_size = base_size) +
        theme(
            plot.caption = element_text(hjust = 0, size = base_size - 3, colour = "grey35"),
            plot.subtitle = element_text(colour = "grey25"),
            legend.position = "bottom"
        )
}

# 3. Setup Directories
# ------------------------------------------------------------------------------
DIR_OUT <- opt$outdir
DIR_PLOTS <- file.path(DIR_OUT, "plots")
dir.create(DIR_OUT, recursive = TRUE, showWarnings = FALSE)
dir.create(DIR_PLOTS, recursive = TRUE, showWarnings = FALSE)
dir.create(DIR_PLOTS_VF, recursive = TRUE, showWarnings = FALSE)

# 4. Load Data
# ------------------------------------------------------------------------------
msg("Loading input datasets...")

# Clinical status
FILE_STATUS <- FILE_STATUS_MAP
if (!file.exists(FILE_STATUS)) stop("Missing ", FILE_STATUS)
status <- read_csv(FILE_STATUS, show_col_types = FALSE) %>%
    filter(Infection_Status %in% c("UTI", "ASB")) %>%
    mutate(
        Participant_id = as.character(Participant_id),
        Timepoint = if ("tp_lab" %in% names(.)) normalize_analysis_timepoint(tp_lab) else normalize_analysis_timepoint(Timepoint),
        Event_type = if ("Event_type" %in% names(.)) as.character(Event_type) else episode_event_type(Timepoint),
        Collection_Date = if ("Collection_Date" %in% names(.)) as.character(Collection_Date) else NA_character_,
        Episode_ID = if ("Episode_ID" %in% names(.)) as.character(Episode_ID) else build_episode_id(., timepoint_col = "Timepoint", event_col = "Event_type", date_col = "Collection_Date"),
        Batch = if ("Batch" %in% names(.)) as.factor(Batch) else NA,
        Outcome = as.integer(Infection_Status == "UTI")
    ) %>%
    select(Participant_id, Timepoint, Episode_ID, Batch, Infection_Status, Outcome)

msg(
    "  Loaded %d samples with UTI/ASB status (%d UTI, %d ASB)",
    nrow(status), sum(status$Outcome == 1), sum(status$Outcome == 0)
)

# VF/AMR matrix - prefer canonical episode-level P/A matrix from script 02.
if (!file.exists(FILE_VF_PA)) stop("Missing ", FILE_VF_PA, ". Run 02_gene_presence_analysis.R first.")
vf_raw <- read_csv(FILE_VF_PA, show_col_types = FALSE) %>%
    mutate(
        Participant_id = as.character(Participant_id),
        Timepoint = normalize_analysis_timepoint(tp_lab)
    )
vf_gene_cols <- canonical_vf_gene_cols(names(vf_raw), vf_pa_file = FILE_VF_PA)
vf_raw <- vf_raw %>%
    select(Participant_id, Timepoint, all_of(vf_gene_cols))

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
FILE_MLST <- FILE_MLST_CANONICAL
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
    select(any_of(c("Assembly_ID", "Isolate_ID", "file_name", "Participant_id", "Timepoint", "tp_lab", "full_path", "usable_fasta"))) %>%
    mutate(
        Participant_id = as.character(Participant_id),
        Timepoint = if ("tp_lab" %in% names(.)) normalize_analysis_timepoint(tp_lab) else normalize_analysis_timepoint(Timepoint)
    )

canonical_file <- file.path(DIR_QC, "canonical_assembly_selection.csv")
if (file.exists(canonical_file)) {
    metadata <- read_csv(canonical_file, show_col_types = FALSE) %>%
        filter(selected_canonical %in% TRUE, QC_PASS %in% TRUE) %>%
        mutate(
            Participant_id = as.character(Participant_id),
            Timepoint = normalize_analysis_timepoint(tp_lab)
        ) %>%
        select(any_of(c("Assembly_ID", "Isolate_ID", "file_name", "Participant_id", "Timepoint", "full_path")))
}
metadata <- metadata %>% distinct(Participant_id, Timepoint, .keep_all = TRUE)

# ============================= DATA HARMONIZATION =============================
msg("Harmonizing datasets...")

analysis_keys <- metadata %>%
    filter(!is.na(full_path), file.exists(full_path)) %>%
    distinct(Participant_id, Timepoint)

status_without_assembly <- status %>%
    anti_join(analysis_keys, by = c("Participant_id", "Timepoint"))

if (nrow(status_without_assembly) > 0) {
    msg(
        "  Pre-bridge direct-key diagnostic: %d UTI/ASB sample(s) do not match assembly keys directly (%d UTI, %d ASB)",
        nrow(status_without_assembly),
        sum(status_without_assembly$Outcome == 1),
        sum(status_without_assembly$Outcome == 0)
    )
}

status <- status %>%
    inner_join(analysis_keys, by = c("Participant_id", "Timepoint"), relationship = "one-to-one")

msg(
    "  Direct-key diagnostic retained %d UTI/ASB samples with usable assemblies (%d UTI, %d ASB)",
    nrow(status), sum(status$Outcome == 1), sum(status$Outcome == 0)
)

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
vf_feature_cols <- vf_gene_cols
data_vf <- status %>%
    left_join(vf, by = c("Participant_id", "Timepoint"), relationship = "one-to-one") %>%
    mutate(across(all_of(vf_feature_cols), ~ tidyr::replace_na(as.integer(.), 0L)))

msg("  Direct-key VF join diagnostic: %d samples", nrow(data_vf))

# Plasmid data: Need to link Isolate_ID -> Participant_id/Timepoint
data_plasmid <- NULL
if (!is.null(plasmid_raw)) {
    # Link plasmid Isolate_ID to Participant_id/Timepoint via metadata
    metadata_plasmid <- metadata %>%
        mutate(plasmid_sample_id = if ("file_name" %in% names(.)) file_name else Isolate_ID) %>%
        select(plasmid_sample_id, Participant_id, Timepoint, full_path) %>%
        distinct()
    plasmid_linked <- plasmid_raw %>%
        rename(plasmid_sample_id = Isolate_ID) %>%
        inner_join(metadata_plasmid, by = "plasmid_sample_id", relationship = "many-to-one") %>%
        select(-plasmid_sample_id, -full_path)

    # Normalize Timepoint format to match status map
    # Extract timepoint from format like "T0", "T1", "Uricult"
    plasmid_linked <- plasmid_linked %>%
        mutate(Timepoint = normalize_analysis_timepoint(Timepoint))

    # For each Participant_id x Timepoint, take max across isolates (multiple assemblies)
    rep_cols <- setdiff(names(plasmid_linked), c("Participant_id", "Timepoint"))
    plasmid_agg <- plasmid_linked %>%
        group_by(Participant_id, Timepoint) %>%
        summarise(across(all_of(rep_cols), ~ max(as.numeric(.), na.rm = TRUE)), .groups = "drop") %>%
        mutate(across(all_of(rep_cols), ~ ifelse(is.finite(.), ., 0)))

    # Join with status
    data_plasmid <- status %>%
        left_join(plasmid_agg, by = c("Participant_id", "Timepoint"), relationship = "one-to-one") %>%
        mutate(across(all_of(rep_cols), ~ tidyr::replace_na(as.integer(.), 0L)))

    msg("  Direct-key plasmid join diagnostic: %d samples", nrow(data_plasmid))
}

# Combine VF + Plasmid features
if (!is.null(data_plasmid)) {
    # Merge on status columns
    data_merged <- data_vf %>%
        left_join(
            data_plasmid %>% select(Participant_id, Timepoint, all_of(rep_cols)),
            by = c("Participant_id", "Timepoint"),
            relationship = "one-to-one"
        )
} else {
    data_merged <- data_vf
}

# Prefer the canonical VF-ready episode table from script 22 when available.
# This table carries the audited Uricult bridge, so Uricult clinical episodes
# mapped to UTI-N WGS rows are retained by clinical Episode_ID instead of being
# lost by a raw Participant_id + Timepoint join.
if (file.exists(FILE_VF_READY)) {
    if (file.info(FILE_VF_READY)$mtime < file.info(FILE_VF_PA)$mtime) {
        stop(
            "vf_analysis_ready.csv is older than vf_pa_all.csv. ",
            "Run Rscript 22_vf_build_analysis_dataset.R before 14_genotype_phenotype_model.R."
        )
    }

    vf_ready_model <- read_csv(FILE_VF_READY, show_col_types = FALSE) %>%
        mutate(
            Participant_id = as.character(Participant_id),
            Timepoint = normalize_analysis_timepoint(tp_lab),
            Infection_Status = as.character(Infection_Status),
            Outcome = as.integer(Infection_Status == "UTI"),
            Batch = if ("Batch" %in% names(.)) as.factor(Batch) else NA
        ) %>%
        filter(Infection_Status %in% c("UTI", "ASB"))

    vf_feature_cols <- canonical_vf_gene_cols(names(vf_ready_model), vf_pa_file = FILE_VF_PA)
    model_meta_cols <- intersect(
        c(
            "Participant_id", "Timepoint", "tp_lab", "Episode_ID", "Collection_Date",
            "Batch", "Infection_Status", "Outcome", "ST", "uricult_bridge_applied"
        ),
        names(vf_ready_model)
    )

    data_merged <- vf_ready_model %>%
        select(all_of(model_meta_cols), all_of(vf_feature_cols))

    if (exists("plasmid_agg") && exists("rep_cols") && length(rep_cols) > 0) {
        data_merged <- data_merged %>%
            left_join(plasmid_agg, by = c("Participant_id", "Timepoint"), relationship = "one-to-one") %>%
            mutate(across(all_of(rep_cols), ~ tidyr::replace_na(as.integer(.), 0L)))
    }

    msg(
        "  Using canonical vf_analysis_ready.csv for modelling: %d samples (%d UTI, %d ASB; %d Uricult-bridged)",
        nrow(data_merged),
        sum(data_merged$Outcome == 1, na.rm = TRUE),
        sum(data_merged$Outcome == 0, na.rm = TRUE),
        if ("uricult_bridge_applied" %in% names(data_merged)) sum(data_merged$uricult_bridge_applied %in% TRUE, na.rm = TRUE) else 0L
    )
}

base_cols <- intersect(
    c(
        "Participant_id", "Timepoint", "tp_lab", "Episode_ID", "Collection_Date",
        "Batch", "Infection_Status", "Outcome", "ST", "uricult_bridge_applied"
    ),
    names(data_merged)
)
feature_cols <- setdiff(names(data_merged), base_cols)

msg(
    "  Final merged dataset: %d samples, %d features",
    nrow(data_merged), length(feature_cols)
)

# Feature filtering: Remove low/high prevalence
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
    select(all_of(base_cols), all_of(features_keep))

model_denom <- data_final %>%
    select(any_of(c("Participant_id", "Timepoint", "Episode_ID", "Batch", "Infection_Status", "Outcome"))) %>%
    mutate(model_interpretation = case_when(
        sum(Outcome == 1, na.rm = TRUE) < 10 ~ "not_interpretable_sparse_UTI",
        sum(Outcome == 1, na.rm = TRUE) < 20 ~ "exploratory_underpowered_UTI",
        TRUE ~ "standard"
    ))
write_csv(model_denom, file.path(DIR_OUT, "model_dataset_denominator.csv"))
append_denominator_summary(
    data_final,
    "14_genotype_phenotype_model.R",
    "model_dataset",
    "participant_timepoint",
    file.path(DIR_OUT, "model_dataset_denominator.csv"),
    "ASB/UTI model dataset built from canonical vf_analysis_ready.csv so audited Uricult bridge rows are retained"
)

n_uti_model <- sum(data_final$Outcome == 1, na.rm = TRUE)
model_warnings <- c(
    "Genotype-phenotype model interpretation warnings",
    sprintf("Generated: %s", format(Sys.time())),
    sprintf("Model dataset: %d samples (%d UTI, %d ASB)", nrow(data_final), n_uti_model, sum(data_final$Outcome == 0, na.rm = TRUE))
)
if (n_uti_model < 20) {
    model_warnings <- c(model_warnings, "WARNING: UTI < 20. Association models are exploratory only and underpowered for definitive UTI inference.")
}
if (n_uti_model < 10) {
    model_warnings <- c(model_warnings, "RED: UTI < 10. Association testing is not interpretable; use descriptive summaries only.")
}

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
    future.seed = TRUE,
    future.packages = c("dplyr", "tibble")
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
            covariate_terms <- c()
            if (n_distinct(data$Timepoint) > 1) covariate_terms <- c(covariate_terms, "Timepoint")
            if ("Batch" %in% names(data) && n_distinct(data$Batch, na.rm = TRUE) > 1) covariate_terms <- c(covariate_terms, "Batch")
            fixed_part <- paste(c(paste0("`", feature, "`"), covariate_terms), collapse = " + ")
            fml <- as.formula(sprintf("Outcome ~ %s", fixed_part))

            # Fit model
            mod <- stats::glm(fml, data = data, family = binomial(link = "logit"))

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
# Mixed-effects GLMM (with participant random intercepts)
run_glmm <- function(feature, data) {
    # Helper to format result
    format_res <- function(mod, type, converged) {
        tidy_mod <- if (inherits(mod, "merMod")) {
            broom.mixed::tidy(mod, conf.int = TRUE, effects = "fixed")
        } else {
            broom::tidy(mod, conf.int = TRUE)
        }
        if (!"p.value" %in% names(tidy_mod)) {
            tidy_mod <- tidy_mod %>%
                mutate(p.value = if ("statistic" %in% names(.)) {
                    2 * pnorm(abs(.data$statistic), lower.tail = FALSE)
                } else {
                    NA_real_
                })
        }
        tidy_mod %>%
            filter(term == paste0("`", feature, "`") | term == feature) %>%
            mutate(
                feature = feature,
                OR = exp(estimate),
                OR_lower = exp(conf.low),
                OR_upper = exp(conf.high),
                converged = converged,
                model_type = type
            ) %>%
            select(feature, estimate, std.error, OR, OR_lower, OR_upper, p.value, converged, model_type)
    }

    # Empty result helper
    empty_res <- function(type = "GLMM") {
        tibble(
            feature = feature,
            estimate = NA_real_,
            std.error = NA_real_,
            OR = NA_real_,
            OR_lower = NA_real_,
            OR_upper = NA_real_,
            p.value = NA_real_,
            converged = FALSE,
            model_type = type
        )
    }

    tryCatch(
        {
            # Build formula
            # [STAT] NOTE ON CONFOUNDING:
            # We adjust for Timepoint and Batch if available. 
            # Bacterial lineage (ST) is NOT automatically included as a covariate here
            # because adding a factor with many levels can cause severe convergence issues.
            # If script 25_vf_lineage_vf_interaction.R indicates significant lineage 
            # confounding, you should manually add `ST` to the covariate_terms here 
            # (or run a sensitivity analysis on a single dominant ST).
            has_batch <- "Batch" %in% names(data) && n_distinct(data$Batch, na.rm = TRUE) > 1
            covariate_terms <- c()
            if (n_distinct(data$Timepoint) > 1) covariate_terms <- c(covariate_terms, "Timepoint")
            if (has_batch) covariate_terms <- c(covariate_terms, "Batch")

            fixed_part <- paste(c(paste0("`", feature, "`"), covariate_terms), collapse = " + ")

            # 1. Try GLMM
            fml_glmm <- as.formula(sprintf("Outcome ~ %s + (1 | Participant_id)", fixed_part))

            mod <- lme4::glmer(fml_glmm,
                data = data, family = binomial(link = "logit"),
                control = lme4::glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))
            )

            is_singular <- lme4::isSingular(mod)
            if (is_singular) {
                # Singular fit often means random effect variance is ~0.
                # We can treat this as converged for practical purposes or fallback.
                # Here we'll accept it but maybe flag it? For now, treat as GLMM.
                format_res(mod, "GLMM (Singular)", TRUE)
            } else {
                format_res(mod, "GLMM", TRUE)
            }
        },
        error = function(e) {
            # 2. Fallback to GLM
            tryCatch(
                {
                    fml_glm <- as.formula(sprintf("Outcome ~ %s", fixed_part))
                    mod_glm <- stats::glm(fml_glm, data = data, family = binomial(link = "logit"))
                    format_res(mod_glm, "GLM_Fallback", TRUE)
                },
                error = function(e2) {
                    empty_res("Failed")
                }
            )
        }
    )
}

# Parallel execution - choose model type
plan(multisession, workers = CORES_USE)
if (opt$`simple-glm`) {
    glmm_results <- future.apply::future_lapply(glmm_features, run_simple_glm,
        data = data_final,
        future.seed = TRUE,
        future.packages = c("broom", "dplyr", "tibble")
    ) %>%
        bind_rows()
} else {
    glmm_results <- future.apply::future_lapply(glmm_features, run_glmm,
        data = data_final,
        future.seed = TRUE,
        future.packages = c("broom", "broom.mixed", "dplyr", "lme4", "tibble")
    ) %>%
        bind_rows()
}
plan(sequential)

# FDR correction
glmm_results <- glmm_results %>%
    filter(!is.na(p.value)) %>%
    mutate(
        FDR = p.adjust(p.value, method = "BH"),
        sparse_data_separation_risk = !is.finite(OR) | !is.finite(OR_lower) | !is.finite(OR_upper) |
            OR_upper > 100 | OR_lower < 0.01 | std.error > 5,
        interpretation = case_when(
            n_uti_model < 10 ~ "not_interpretable_sparse_UTI",
            n_uti_model < 20 ~ "exploratory_underpowered_UTI",
            str_detect(model_type, "Singular") ~ "exploratory_singular_fit",
            sparse_data_separation_risk ~ "exploratory_sparse_data_separation_risk",
            TRUE ~ "exploratory_screening"
        )
    ) %>%
    arrange(p.value)

if (any(str_detect(glmm_results$model_type, "Singular"), na.rm = TRUE)) {
    model_warnings <- c(model_warnings, sprintf("WARNING: %d GLMM result(s) had singular random-effect fits.",
                                                sum(str_detect(glmm_results$model_type, "Singular"), na.rm = TRUE)))
}
if (any(glmm_results$sparse_data_separation_risk, na.rm = TRUE)) {
    model_warnings <- c(model_warnings, sprintf("WARNING: %d model result(s) have extreme OR/CI or large SE consistent with sparse-data/separation risk.",
                                                sum(glmm_results$sparse_data_separation_risk, na.rm = TRUE)))
}
writeLines(model_warnings, file.path(DIR_OUT, "model_interpretation_warnings.txt"))

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
write_csv(glmm_results, file.path(DIR_OUT, "gwas_multivariable_glmm.csv"))
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
    scale_color_manual(values = c("TRUE" = "red", "FALSE" = "grey40"), labels = c("Not Significant", sprintf("FDR < %.2f", opt$fdr_thresh))) +
    labs(
        title = "Exploratory Genotype-Phenotype Association: UTI vs. ASB",
        subtitle = sprintf("%d features tested (Fisher's exact); UTI n=%d", nrow(univar_results), n_uti_model),
        x = "Log2 Odds Ratio",
        y = "-Log10 P-value",
        color = "Significance"
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
    slice_min(FDR, n = 20, with_ties = FALSE)

if (nrow(top_glmm) > 0) {
    p_forest <- ggplot(top_glmm, aes(x = OR, y = reorder(feature, OR))) +
        geom_vline(xintercept = 1, linetype = "dashed", color = "grey50") +
        geom_point(size = 3, color = "steelblue") +
        geom_errorbarh(aes(xmin = OR_lower, xmax = OR_upper), height = 0.2) +
        geom_text(aes(label = sprintf("%.2f", OR)), hjust = -0.5, size = 3) +
        scale_x_log10() +
        labs(
            title = "Exploratory Top GLMM Associations (UTI vs. ASB)",
            subtitle = sprintf("%d features with FDR < 0.2; sparse data flags retained", nrow(top_glmm)),
            x = "Odds Ratio (95% CI)",
            y = NULL
        ) +
        theme_minimal(base_size = 11)

    ggsave(file.path(DIR_PLOTS, "forest_plot_top_hits.png"),
        p_forest,
        width = 8, height = max(6, nrow(top_glmm) * 0.3), dpi = 300
    )
}

# VF interpretation bridge: univariable Fisher screening versus participant-aware
# model output. This deliberately lives in plots/vf because its purpose is to
# explain VF interpretation, even though it reuses model outputs from this script.
bridge_top_n <- 12
bridge_fisher <- univar_results %>%
    filter(!is.na(p_value)) %>%
    arrange(p_value) %>%
    slice_head(n = bridge_top_n) %>%
    transmute(
        feature,
        fisher_OR = OR,
        fisher_lower = CI_lower,
        fisher_upper = CI_upper,
        fisher_p = p_value,
        fisher_FDR = FDR
    )

bridge_model <- glmm_results %>%
    transmute(
        feature,
        model_OR = OR,
        model_lower = OR_lower,
        model_upper = OR_upper,
        model_p = p.value,
        model_FDR = FDR,
        model_type,
        sparse_data_separation_risk
    )

bridge_join <- bridge_fisher %>%
    left_join(bridge_model, by = "feature") %>%
    mutate(
        feature = factor(feature, levels = rev(feature)),
        model_status = case_when(
            is.na(model_OR) ~ "Not modelled",
            sparse_data_separation_risk %in% TRUE ~ "Model sparse/separation risk",
            !is.na(model_FDR) & model_FDR < opt$fdr_thresh ~ "Model FDR-significant",
            TRUE ~ "Model not FDR-significant"
        ),
        fisher_status = ifelse(fisher_p < 0.05, "Nominal Fisher p < 0.05", "Fisher p >= 0.05")
    )

bridge_plot_data <- bind_rows(
    bridge_join %>%
        transmute(
            feature,
            evidence_layer = "Univariable Fisher screen",
            estimate = fisher_OR,
            lower = fisher_lower,
            upper = fisher_upper,
            status = fisher_status,
            label = sprintf("p=%.2g; q=%.2g", fisher_p, fisher_FDR)
        ),
    bridge_join %>%
        filter(!is.na(model_OR)) %>%
        transmute(
            feature,
            evidence_layer = "Participant-aware model",
            estimate = model_OR,
            lower = model_lower,
            upper = model_upper,
            status = model_status,
            label = sprintf("q=%.2g%s", model_FDR,
                            ifelse(sparse_data_separation_risk %in% TRUE, "; sparse", ""))
        )
) %>%
    mutate(
        evidence_layer = factor(evidence_layer,
                                levels = c("Univariable Fisher screen", "Participant-aware model")),
        estimate = ifelse(is.finite(estimate) & estimate > 0, estimate, NA_real_),
        lower = ifelse(is.finite(lower) & lower > 0, lower, NA_real_),
        upper = ifelse(is.finite(upper) & upper > 0, upper, NA_real_)
    ) %>%
    filter(!is.na(estimate))

if (nrow(bridge_plot_data) > 0) {
    n_model_sig <- sum(glmm_results$FDR < opt$fdr_thresh, na.rm = TRUE)
    n_sparse <- sum(glmm_results$sparse_data_separation_risk %in% TRUE, na.rm = TRUE)
    p_bridge <- ggplot(bridge_plot_data,
                       aes(x = estimate, y = feature, xmin = lower, xmax = upper, colour = status)) +
        geom_vline(xintercept = 1, linetype = "dashed", colour = "grey50") +
        geom_errorbarh(height = 0.18, linewidth = 0.45, na.rm = TRUE) +
        geom_point(size = 2.6, na.rm = TRUE) +
        geom_text(aes(label = label), hjust = -0.05, size = 2.8, show.legend = FALSE, na.rm = TRUE) +
        facet_wrap(~evidence_layer, nrow = 1) +
        scale_x_log10() +
        coord_cartesian(xlim = c(0.03, 80), clip = "off") +
        scale_colour_manual(values = c(
            "Nominal Fisher p < 0.05" = "#D55E00",
            "Fisher p >= 0.05" = "grey45",
            "Model FDR-significant" = "#009E73",
            "Model not FDR-significant" = "#0072B2",
            "Model sparse/separation risk" = "#CC79A7",
            "Not modelled" = "grey70"
        )) +
        labs(
            title = "Exploratory VF gene screening versus participant-aware model evidence",
            subtitle = sprintf(
                "Top %d Fisher-ranked features; model results have %d FDR-significant hits and %d sparse/separation flags",
                bridge_top_n, n_model_sig, n_sparse
            ),
            x = "Odds ratio on log scale",
            y = "VF gene / feature",
            colour = "Evidence status",
            caption = paste(
                sprintf("Data: %s and %s.",
                        file.path(DIR_OUT, "gwas_univariable_stats.csv"),
                        file.path(DIR_OUT, "gwas_multivariable_glmm.csv")),
                sprintf("Denominator: ASB/UTI model dataset n=%d, including UTI n=%d.", nrow(data_final), n_uti_model),
                "Level of analysis: gene/feature-level association screen.",
                "Fisher tests are exploratory and do not account for repeated resident isolates.",
                "Participant-aware models include a Participant_id random effect where GLMM fitting is possible, but sparse or singular fits remain underpowered.",
                "No cross-sectional VF comparison should be interpreted causally without considering ST/lineage, timepoint, batch, and event-type structure."
            )
        ) +
        plot_theme_vf_model(base_size = 10) +
        theme(
            plot.margin = margin(5.5, 36, 5.5, 5.5),
            strip.background = element_rect(fill = "grey92", colour = "grey70")
        )

    ggsave(file.path(DIR_PLOTS_VF, "vf_gene_screening_vs_model_evidence.png"),
           p_bridge, width = 11, height = 6.5, dpi = 300)
}

# Heatmap of top discriminators
top_features <- glmm_results %>%
    filter(converged, !is.na(FDR)) %>%
    slice_min(FDR, n = 30, with_ties = FALSE) %>%
    pull(feature)

if (length(top_features) > 1) {
    # Create unique Sample_ID to avoid duplicate row names
    heatmap_data <- data_final %>%
        mutate(Sample_ID = ifelse(!is.na(Timepoint),
            paste(Participant_id, Timepoint, sep = "_"),
            make.unique(as.character(Participant_id))
        )) %>%
        select(Sample_ID, Infection_Status, all_of(top_features)) %>%
        arrange(Infection_Status, Sample_ID) %>%
        column_to_rownames("Sample_ID") %>%
        select(-Infection_Status)

    annotation <- data_final %>%
        mutate(Sample_ID = ifelse(!is.na(Timepoint),
            paste(Participant_id, Timepoint, sep = "_"),
            make.unique(as.character(Participant_id))
        )) %>%
        select(Sample_ID, Infection_Status) %>%
        distinct() %>%
        column_to_rownames("Sample_ID")

    save_png_device(
        file.path(DIR_PLOTS, "heatmap_top_discriminators.png"),
        width = 10,
        height = max(6, length(top_features) * 0.2),
        draw_fun = function() {
            pheatmap(
                t(as.matrix(heatmap_data)),
                annotation_col = annotation,
                cluster_rows = TRUE,
                cluster_cols = FALSE,
                show_colnames = FALSE,
                color = c("white", "darkblue"),
                border_color = NA,
                main = "Exploratory Top Discriminatory Features (UTI vs. ASB)"
            )
        }
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
    slice_min(FDR, n = 10, with_ties = FALSE)

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

write_uti_attrition_outputs()
msg("Done.")

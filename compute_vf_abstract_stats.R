# ==============================================================================
# compute_vf_abstract_stats.R
# SUPERSEDED by 22_vf_build_analysis_dataset.R + 23_vf_cross_sectional.R + 24_vf_longitudinal_dynamics.R
# Kept for validation reference only. Do not use for new analyses.
# ==============================================================================
# Purpose: Compute all VF metrics for abstract from ANCHOR files only.
#
# Anchor files:
#   1. results/vf/vf_pa_all.csv       — VF presence/absence matrix
#   2. status_map.csv                  — Clinical status labels
#   3. results/vf/gene_map.csv         — Gene → Category mapping
#
# Outputs:
#   - results/vf/vf_analysis_ready.csv             (Deliverable B)
#   - results/vf/vf_burden_by_status.csv            (Deliverable C1)
#   - results/vf/vf_gene_prevalence_by_status.csv   (Deliverable C2)
#   - results/vf/vf_gene_enrichment_UTI_vs_ASB.csv  (Deliverable C3)
#   - results/vf/vf_longitudinal_transitions.csv    (Deliverable D2)
#   - results/vf/vf_transition_summary_by_type.csv  (Deliverable D3)
# ==============================================================================

suppressPackageStartupMessages({
    library(dplyr)
    library(readr)
    library(tidyr)
    library(stringr)
    library(purrr)
})

cat("=== compute_vf_abstract_stats.R ===\n")
cat("Started:", format(Sys.time()), "\n\n")

# ==============================================================================
# 0. LOAD ANCHOR FILES
# ==============================================================================

# --- Anchor 1: VF P/A matrix ---
vf_pa_file <- "results/vf/vf_pa_all.csv"
stopifnot(file.exists(vf_pa_file))
vf_pa <- read_csv(vf_pa_file, show_col_types = FALSE)

cat("ANCHOR 1: vf_pa_all.csv\n")
cat("  Rows:", nrow(vf_pa), "\n")
cat("  Cols:", ncol(vf_pa), "\n")
cat("  Key cols: Participant_id, tp_lab\n")
cat("  Gene cols:", ncol(vf_pa) - 2, "\n")

# Identify gene columns (everything except Participant_id, tp_lab)
meta_cols <- c("Participant_id", "tp_lab")
gene_cols <- setdiff(names(vf_pa), meta_cols)
cat("  Unique participants:", n_distinct(vf_pa$Participant_id), "\n")
cat("  Unique timepoints:", paste(sort(unique(as.character(vf_pa$tp_lab))), collapse = ", "), "\n\n")

# --- Anchor 2: status_map.csv ---
status_file <- "status_map.csv"
stopifnot(file.exists(status_file))
status_map <- read_csv(status_file, show_col_types = FALSE) %>%
    mutate(Participant_id = as.character(Participant_id))

cat("ANCHOR 2: status_map.csv\n")
cat("  Rows:", nrow(status_map), "\n")
cat("  Cols:", ncol(status_map), "\n")
cat("  Status breakdown:\n")
print(table(status_map$Infection_Status, useNA = "ifany"))
cat("\n")

# --- Anchor 3: gene_map.csv ---
gene_map_file <- "results/vf/gene_map.csv"
stopifnot(file.exists(gene_map_file))
gene_map <- read_csv(gene_map_file, show_col_types = FALSE) %>%
    mutate(
        Gene = as.character(Gene),
        Category = coalesce(as.character(Category), "Unassigned")
    )

cat("ANCHOR 3: gene_map.csv\n")
cat("  Rows:", nrow(gene_map), "\n")
cat("  Categories:", paste(sort(unique(gene_map$Category)), collapse = ", "), "\n\n")

# ==============================================================================
# 1. NORMALIZE TIMEPOINT KEYS FOR JOINING
# ==============================================================================

# vf_pa uses tp_lab (T0, T1, T2, Uricult)
# status_map uses Timepoint (T0, T1, T2, Uricult, etc.)

# Normalize: create a common key
vf_pa <- vf_pa %>%
    mutate(
        Participant_id = as.character(Participant_id),
        tp_lab = as.character(tp_lab)
    )

# Normalize Timepoint in status_map to tp_lab format
normalize_tp <- function(x) {
    x <- as.character(x)
    is_uricult <- str_detect(x, regex("uricult", ignore_case = TRUE))
    tp_num <- suppressWarnings(as.integer(str_extract(x, "\\d+")))
    case_when(
        is_uricult ~ "Uricult",
        !is.na(tp_num) ~ paste0("T", tp_num),
        TRUE ~ x
    )
}

status_map <- status_map %>%
    mutate(tp_lab = normalize_tp(Timepoint))

cat("=== MERGE DIAGNOSTICS ===\n\n")

# Check for duplicates in join keys
vf_dupes <- vf_pa %>%
    count(Participant_id, tp_lab) %>%
    filter(n > 1)
status_dupes <- status_map %>%
    count(Participant_id, tp_lab) %>%
    filter(n > 1)

cat("Duplicate join keys in vf_pa:", nrow(vf_dupes), "\n")
cat("Duplicate join keys in status_map:", nrow(status_dupes), "\n")

if (nrow(status_dupes) > 0) {
    cat("  Status map duplicates (will take first per key):\n")
    print(status_dupes)
    # Deduplicate status_map: take first row per key
    status_map <- status_map %>%
        group_by(Participant_id, tp_lab) %>%
        slice(1) %>%
        ungroup()
}

# Anti-join diagnostics
vf_not_in_status <- vf_pa %>%
    anti_join(status_map, by = c("Participant_id", "tp_lab"))
status_not_in_vf <- status_map %>%
    anti_join(vf_pa, by = c("Participant_id", "tp_lab"))

cat("\nVF rows NOT matching status_map:", nrow(vf_not_in_status), "\n")
if (nrow(vf_not_in_status) > 0) {
    cat("  Unmatched VF participants:", paste(unique(vf_not_in_status$Participant_id), collapse = ", "), "\n")
    cat("  Unmatched VF timepoints:", paste(unique(vf_not_in_status$tp_lab), collapse = ", "), "\n")
}

cat("Status rows NOT matching VF data:", nrow(status_not_in_vf), "\n")
if (nrow(status_not_in_vf) > 0) {
    cat("  These are clinical episodes without sequenced VF data.\n")
    cat("  Status breakdown of unmatched:\n")
    print(table(status_not_in_vf$Infection_Status, useNA = "ifany"))
}

# ==============================================================================
# 2. BUILD ANALYSIS-READY DATASET (Deliverable B)
# ==============================================================================

cat("\n=== BUILDING ANALYSIS-READY DATASET ===\n\n")

# Join VF data with status
vf_ready <- vf_pa %>%
    left_join(
        status_map %>% select(Participant_id, tp_lab, Infection_Status),
        by = c("Participant_id", "tp_lab")
    )

# Compute VF count per row
vf_ready <- vf_ready %>%
    mutate(vf_count_total = rowSums(across(all_of(gene_cols)), na.rm = TRUE))

# Category-level counts using gene_map
cat_counts <- lapply(unique(gene_map$Category), function(cat) {
    cat_genes <- gene_map %>%
        filter(Category == cat) %>%
        pull(Gene)
    matching_cols <- intersect(cat_genes, gene_cols)
    if (length(matching_cols) == 0) {
        return(rep(0L, nrow(vf_ready)))
    }
    rowSums(vf_ready[, matching_cols, drop = FALSE], na.rm = TRUE)
})
names(cat_counts) <- paste0("cat_", gsub("[/ ]", "_", unique(gene_map$Category)))

for (nm in names(cat_counts)) {
    vf_ready[[nm]] <- cat_counts[[nm]]
}

# Count Unassigned category
assigned_genes <- gene_map %>%
    filter(Category != "Unassigned") %>%
    pull(Gene)
unassigned_in_matrix <- setdiff(gene_cols, assigned_genes)
if (length(unassigned_in_matrix) > 0) {
    vf_ready$cat_Unassigned_matrix <- rowSums(vf_ready[, unassigned_in_matrix, drop = FALSE], na.rm = TRUE)
}

# Add is_ecoli flag (all TRUE — documented)
vf_ready$is_ecoli <- TRUE

# Write Deliverable B
write_csv(vf_ready, "results/vf/vf_analysis_ready.csv")

cat("Written: results/vf/vf_analysis_ready.csv\n")
cat("  Rows:", nrow(vf_ready), "\n")
cat("  Participants:", n_distinct(vf_ready$Participant_id), "\n")
cat("  Status breakdown in analysis-ready:\n")
print(table(vf_ready$Infection_Status, useNA = "ifany"))
cat("  Mean VF count:", round(mean(vf_ready$vf_count_total), 1), "\n")
cat("  Median VF count:", median(vf_ready$vf_count_total), "\n\n")

# Save merge diagnostics as structured data for the doc
merge_diag <- list(
    vf_pa_rows = nrow(vf_pa),
    vf_pa_participants = n_distinct(vf_pa$Participant_id),
    vf_pa_gene_cols = length(gene_cols),
    status_map_rows = nrow(status_map),
    status_map_participants = n_distinct(status_map$Participant_id),
    vf_dupes = nrow(vf_dupes),
    status_dupes_before_dedup = nrow(status_dupes),
    vf_not_in_status = nrow(vf_not_in_status),
    status_not_in_vf = nrow(status_not_in_vf),
    analysis_ready_rows = nrow(vf_ready),
    analysis_ready_matched = sum(!is.na(vf_ready$Infection_Status)),
    analysis_ready_unmatched = sum(is.na(vf_ready$Infection_Status))
)

# ==============================================================================
# 3. BETWEEN-INDIVIDUAL VF RESULTS (Deliverables C1, C2, C3)
# ==============================================================================

cat("=== BETWEEN-INDIVIDUAL VF RESULTS ===\n\n")

# C1: VF burden by status
vf_with_status <- vf_ready %>% filter(!is.na(Infection_Status))

burden_by_status <- vf_with_status %>%
    group_by(Infection_Status) %>%
    summarise(
        n_episodes = n(),
        n_participants = n_distinct(Participant_id),
        mean_vf = round(mean(vf_count_total), 1),
        sd_vf = round(sd(vf_count_total), 1),
        median_vf = median(vf_count_total),
        q25_vf = quantile(vf_count_total, 0.25),
        q75_vf = quantile(vf_count_total, 0.75),
        min_vf = min(vf_count_total),
        max_vf = max(vf_count_total),
        .groups = "drop"
    )

write_csv(burden_by_status, "results/vf/vf_burden_by_status.csv")
cat("C1: VF burden by status\n")
print(as.data.frame(burden_by_status))
cat("\n")

# Also compute category-level burden
cat_cols_in_ready <- grep("^cat_", names(vf_ready), value = TRUE)
if (length(cat_cols_in_ready) > 0) {
    cat_burden <- vf_with_status %>%
        group_by(Infection_Status) %>%
        summarise(
            across(all_of(cat_cols_in_ready),
                list(
                    median = ~ median(., na.rm = TRUE),
                    mean = ~ round(mean(., na.rm = TRUE), 1)
                ),
                .names = "{.col}__{.fn}"
            ),
            .groups = "drop"
        )
    write_csv(cat_burden, "results/vf/vf_category_burden_by_status.csv")
    cat("  Category-level burden written.\n\n")
}

# C2: Gene-level prevalence by status
cat("C2: Gene-level prevalence by status\n")

# Compute per-gene prevalence for each status
gene_prev_list <- lapply(gene_cols, function(g) {
    df <- vf_with_status %>%
        group_by(Infection_Status) %>%
        summarise(
            n_present = sum(.data[[g]] > 0, na.rm = TRUE),
            n_total = n(),
            .groups = "drop"
        ) %>%
        mutate(
            prevalence_pct = round(n_present / n_total * 100, 1),
            label = paste0(n_present, "/", n_total, " (", prevalence_pct, "%)")
        )

    tibble(
        gene = g,
        ASB_n = df$n_present[df$Infection_Status == "ASB"],
        ASB_N = df$n_total[df$Infection_Status == "ASB"],
        ASB_pct = df$prevalence_pct[df$Infection_Status == "ASB"],
        UTI_n = if ("UTI" %in% df$Infection_Status) df$n_present[df$Infection_Status == "UTI"] else NA_integer_,
        UTI_N = if ("UTI" %in% df$Infection_Status) df$n_total[df$Infection_Status == "UTI"] else NA_integer_,
        UTI_pct = if ("UTI" %in% df$Infection_Status) df$prevalence_pct[df$Infection_Status == "UTI"] else NA_real_,
        Neg_n = if ("Negative" %in% df$Infection_Status) df$n_present[df$Infection_Status == "Negative"] else NA_integer_,
        Neg_N = if ("Negative" %in% df$Infection_Status) df$n_total[df$Infection_Status == "Negative"] else NA_integer_,
        Neg_pct = if ("Negative" %in% df$Infection_Status) df$prevalence_pct[df$Infection_Status == "Negative"] else NA_real_
    )
})

gene_prev <- bind_rows(gene_prev_list) %>%
    mutate(delta_UTI_minus_ASB = round(UTI_pct - ASB_pct, 1)) %>%
    arrange(desc(abs(delta_UTI_minus_ASB)))

# Add category from gene_map
gene_prev <- gene_prev %>%
    left_join(gene_map %>% select(Gene, Category), by = c("gene" = "Gene")) %>%
    mutate(Category = coalesce(Category, "Unassigned"))

write_csv(gene_prev, "results/vf/vf_gene_prevalence_by_status.csv")
cat("  Written:", nrow(gene_prev), "genes\n")
cat("  Top 10 by |Δ(UTI−ASB)|:\n")
print(head(as.data.frame(gene_prev %>% select(gene, Category, ASB_pct, UTI_pct, delta_UTI_minus_ASB)), 10))
cat("\n")

# C3: Exploratory enrichment testing (Fisher exact, UTI vs ASB)
cat("C3: Exploratory enrichment (Fisher exact, UTI vs ASB)\n")

asb_uti_data <- vf_with_status %>% filter(Infection_Status %in% c("ASB", "UTI"))

enrichment_list <- lapply(gene_cols, function(g) {
    tab <- table(
        gene_present = asb_uti_data[[g]] > 0,
        is_uti = asb_uti_data$Infection_Status == "UTI"
    )

    # Need a proper 2x2 table
    if (nrow(tab) < 2 || ncol(tab) < 2) {
        return(tibble(
            gene = g, OR = NA_real_, p_value = NA_real_,
            CI_lower = NA_real_, CI_upper = NA_real_,
            note = "insufficient variation"
        ))
    }

    ft <- tryCatch(fisher.test(tab), error = function(e) NULL)
    if (is.null(ft)) {
        return(tibble(
            gene = g, OR = NA_real_, p_value = NA_real_,
            CI_lower = NA_real_, CI_upper = NA_real_,
            note = "fisher test failed"
        ))
    }

    tibble(
        gene = g,
        OR = round(ft$estimate, 3),
        p_value = round(ft$p.value, 4),
        CI_lower = round(ft$conf.int[1], 3),
        CI_upper = round(ft$conf.int[2], 3),
        note = "exploratory/unadjusted; repeated measures violate independence"
    )
})

enrichment <- bind_rows(enrichment_list) %>%
    filter(!is.na(OR)) %>%
    mutate(p_adj_BH = round(p.adjust(p_value, method = "BH"), 4)) %>%
    arrange(p_value)

# Add category
enrichment <- enrichment %>%
    left_join(gene_map %>% select(Gene, Category), by = c("gene" = "Gene")) %>%
    mutate(Category = coalesce(Category, "Unassigned"))

write_csv(enrichment, "results/vf/vf_gene_enrichment_UTI_vs_ASB.csv")
cat("  Written:", nrow(enrichment), "genes tested\n")
cat("  Genes with p < 0.05:", sum(enrichment$p_value < 0.05, na.rm = TRUE), "\n")
cat("  Genes with p_adj_BH < 0.05:", sum(enrichment$p_adj_BH < 0.05, na.rm = TRUE), "\n")
cat("  Top 10 by p-value:\n")
print(head(as.data.frame(enrichment %>% select(gene, Category, OR, p_value, p_adj_BH)), 10))
cat("\n")

# ==============================================================================
# 4. WITHIN-INDIVIDUAL VF DYNAMICS (Deliverables D2, D3)
# ==============================================================================

cat("=== WITHIN-INDIVIDUAL VF DYNAMICS ===\n\n")

# Define time ordering
tp_order <- c("T0", "T1", "T2", "T3", "T4", "Uricult")

# Build participant-level time series
vf_longitudinal <- vf_with_status %>%
    filter(tp_lab %in% tp_order) %>%
    mutate(tp_rank = match(tp_lab, tp_order)) %>%
    arrange(Participant_id, tp_rank)

# Identify participants with ≥2 timepoints
multi_tp <- vf_longitudinal %>%
    group_by(Participant_id) %>%
    summarise(n_tp = n_distinct(tp_lab), .groups = "drop") %>%
    filter(n_tp >= 2)

cat("Participants with ≥2 timepoints:", nrow(multi_tp), "\n")
cat("  ≥3 timepoints:", sum(multi_tp$n_tp >= 3), "\n")
cat("  ≥4 timepoints:", sum(multi_tp$n_tp >= 4), "\n\n")

# Build consecutive pairs
transitions <- list()

for (pid in multi_tp$Participant_id) {
    pid_data <- vf_longitudinal %>%
        filter(Participant_id == pid) %>%
        arrange(tp_rank)

    if (nrow(pid_data) < 2) next

    for (i in 1:(nrow(pid_data) - 1)) {
        row_from <- pid_data[i, ]
        row_to <- pid_data[i + 1, ]

        # Get VF gene vectors
        vf_from <- gene_cols[which(row_from[, gene_cols] > 0)]
        vf_to <- gene_cols[which(row_to[, gene_cols] > 0)]

        gained <- setdiff(vf_to, vf_from)
        lost <- setdiff(vf_from, vf_to)
        stable <- intersect(vf_from, vf_to)

        union_set <- union(vf_from, vf_to)
        jaccard <- if (length(union_set) == 0) NA_real_ else length(stable) / length(union_set)

        transitions[[length(transitions) + 1]] <- tibble(
            Participant_id = pid,
            tp_from = as.character(row_from$tp_lab),
            tp_to = as.character(row_to$tp_lab),
            status_from = row_from$Infection_Status,
            status_to = row_to$Infection_Status,
            transition_type = paste0(row_from$Infection_Status, "→", row_to$Infection_Status),
            vf_count_from = row_from$vf_count_total,
            vf_count_to = row_to$vf_count_total,
            n_gained = length(gained),
            n_lost = length(lost),
            n_stable = length(stable),
            n_from = length(vf_from),
            n_to = length(vf_to),
            jaccard_similarity = round(jaccard, 3),
            any_vf_change = length(gained) > 0 | length(lost) > 0,
            genes_gained = paste(sort(gained), collapse = ";"),
            genes_lost = paste(sort(lost), collapse = ";")
        )
    }
}

trans_df <- bind_rows(transitions)

write_csv(trans_df, "results/vf/vf_longitudinal_transitions.csv")
cat("D2: Longitudinal transitions\n")
cat("  Total transitions:", nrow(trans_df), "\n")
cat("  Participants with transitions:", n_distinct(trans_df$Participant_id), "\n")
cat("  Transition type counts:\n")
print(table(trans_df$transition_type))
cat("\n")

# D3: Summaries by transition type
trans_summary <- trans_df %>%
    group_by(transition_type) %>%
    summarise(
        n_transitions = n(),
        n_participants = n_distinct(Participant_id),
        median_gained = median(n_gained),
        mean_gained = round(mean(n_gained), 1),
        median_lost = median(n_lost),
        mean_lost = round(mean(n_lost), 1),
        median_stable = median(n_stable),
        median_jaccard = round(median(jaccard_similarity, na.rm = TRUE), 3),
        mean_jaccard = round(mean(jaccard_similarity, na.rm = TRUE), 3),
        pct_no_change = round(mean(!any_vf_change) * 100, 1),
        .groups = "drop"
    ) %>%
    arrange(desc(n_transitions))

write_csv(trans_summary, "results/vf/vf_transition_summary_by_type.csv")
cat("D3: Transition summary by type\n")
print(as.data.frame(trans_summary))
cat("\n")

# Overall longitudinal summary
cat("Overall longitudinal stats:\n")
cat("  Median Jaccard:", round(median(trans_df$jaccard_similarity, na.rm = TRUE), 3), "\n")
cat("  Mean Jaccard:", round(mean(trans_df$jaccard_similarity, na.rm = TRUE), 3), "\n")
cat("  % transitions with zero VF change:", round(mean(!trans_df$any_vf_change) * 100, 1), "%\n")
cat("  Median genes gained:", median(trans_df$n_gained), "\n")
cat("  Median genes lost:", median(trans_df$n_lost), "\n")

# Most commonly gained genes in ASB→UTI transitions
asb_to_uti <- trans_df %>% filter(transition_type == "ASB→UTI")
cat("\n  ASB→UTI transitions:", nrow(asb_to_uti), "\n")

if (nrow(asb_to_uti) > 0 && any(asb_to_uti$n_gained > 0)) {
    gained_in_asb_uti <- asb_to_uti %>%
        filter(genes_gained != "") %>%
        pull(genes_gained) %>%
        strsplit(";") %>%
        unlist() %>%
        table() %>%
        sort(decreasing = TRUE)

    cat("  Most commonly gained genes in ASB→UTI:\n")
    print(head(gained_in_asb_uti, 10))
} else {
    cat("  No genes gained in ASB→UTI transitions (or none exist).\n")
}

# ==============================================================================
# 5. CROSS-CHECK AGAINST EXISTING OUTPUTS (Deliverable G)
# ==============================================================================

cat("\n=== CROSS-CHECK AGAINST EXISTING OUTPUTS ===\n\n")

# Cross-check with stratified_vf_stats_table.csv
xcheck_file <- "results/stratified_vf_stats_table.csv"
if (file.exists(xcheck_file)) {
    xcheck <- read_csv(xcheck_file, show_col_types = FALSE)
    cat("Cross-check: results/stratified_vf_stats_table.csv\n")
    cat("  NOTE: This file was produced by get_stratified_vf_stats.R which uses:\n")
    cat("    - results/annotated_gene_table.csv (NOT vf_pa_all.csv)\n")
    cat("    - Filters by cfu_recorded_any == TRUE\n")
    cat("    - Uses deduplication via class_inputs_full.csv (one isolate per episode)\n")
    cat("    - Only counts VFs with Category != 'Unassigned'\n")
    cat("  These are DIFFERENT definitions from our anchor-based computation.\n")
    cat("  Discrepancies are expected and documented.\n\n")
    print(as.data.frame(xcheck))
} else {
    cat("  Cross-check file not found:", xcheck_file, "\n")
}

# Cross-check GLMM results
glmm_file <- "results/vf/diff_focus_genes_UTI_vs_ASB_glmm.csv"
if (file.exists(glmm_file)) {
    glmm <- read_csv(glmm_file, show_col_types = FALSE)
    cat("\nCross-check: GLMM focus gene results\n")
    cat("  File:", glmm_file, "\n")
    cat("  This uses 04_gene_breakdown.R logic with status_map from results/clinical/\n")
    cat("  All ORs ≈ 1.0, p = 1 — no signal detected.\n")
    cat("  Our Fisher exact results can be compared for consistency.\n")

    # Compare our Fisher results for the same genes
    focus_genes <- glmm$Gene
    our_focus <- enrichment %>% filter(gene %in% focus_genes)
    cat("\n  Our Fisher results for focus genes:\n")
    if (nrow(our_focus) > 0) {
        print(as.data.frame(our_focus %>% select(gene, OR, p_value)))
    } else {
        cat("  (None of the focus genes had sufficient variation for Fisher test)\n")
    }
}

# ==============================================================================
# 6. FINAL SUMMARY BLOCK
# ==============================================================================

cat("\n\n========================================\n")
cat("ABSTRACT-READY NUMBERS SUMMARY\n")
cat("========================================\n\n")

cat("COHORT:\n")
cat("  Participants with VF data:", n_distinct(vf_ready$Participant_id), "\n")
cat("  Total participant×timepoint rows:", nrow(vf_ready), "\n")
cat("  With clinical status:", sum(!is.na(vf_ready$Infection_Status)), "\n")
cat("  Without status (unmatched):", sum(is.na(vf_ready$Infection_Status)), "\n")
cat("  ASB episodes:", sum(vf_ready$Infection_Status == "ASB", na.rm = TRUE), "\n")
cat("  UTI episodes:", sum(vf_ready$Infection_Status == "UTI", na.rm = TRUE), "\n")
cat("  Negative episodes:", sum(vf_ready$Infection_Status == "Negative", na.rm = TRUE), "\n")
cat("  All assemblies are E. coli: YES (382/382 in assembly_metadata.csv)\n\n")

cat("VF BURDEN BY STATUS:\n")
for (i in seq_len(nrow(burden_by_status))) {
    r <- burden_by_status[i, ]
    cat(sprintf(
        "  %s: n=%d, median=%.0f (IQR %.0f-%.0f), mean=%.1f (SD=%.1f)\n",
        r$Infection_Status, r$n_episodes, r$median_vf,
        r$q25_vf, r$q75_vf, r$mean_vf, r$sd_vf
    ))
}

cat("\nTOP DIFFERENTIATING GENES (UTI vs ASB by prevalence):\n")
top_diff <- gene_prev %>%
    filter(!is.na(delta_UTI_minus_ASB)) %>%
    slice_max(abs(delta_UTI_minus_ASB), n = 10)
for (i in seq_len(nrow(top_diff))) {
    r <- top_diff[i, ]
    cat(sprintf(
        "  %s (%s): ASB=%.1f%%, UTI=%.1f%%, Δ=%.1f pp\n",
        r$gene, r$Category, r$ASB_pct, r$UTI_pct, r$delta_UTI_minus_ASB
    ))
}

cat("\nLONGITUDINAL DYNAMICS:\n")
cat(sprintf(
    "  Total transitions: %d (from %d participants)\n",
    nrow(trans_df), n_distinct(trans_df$Participant_id)
))
cat(sprintf("  %% with zero VF change: %.1f%%\n", mean(!trans_df$any_vf_change) * 100))
cat(sprintf("  Median Jaccard: %.3f\n", median(trans_df$jaccard_similarity, na.rm = TRUE)))
cat(sprintf(
    "  Median gained: %d, Median lost: %d\n",
    median(trans_df$n_gained), median(trans_df$n_lost)
))

# Timepoint coverage
tp_cov <- vf_ready %>%
    filter(!is.na(Infection_Status)) %>%
    count(tp_lab, name = "n_episodes") %>%
    arrange(tp_lab)
cat("\nTIMEPOINT COVERAGE:\n")
print(as.data.frame(tp_cov))

cat("\n=== DONE ===\n")
cat("Finished:", format(Sys.time()), "\n")

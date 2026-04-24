# ==============================================================================
# compute_vf_stratified_by_depth.R
# SUPERSEDED by 23_vf_cross_sectional.R + 24_vf_longitudinal_dynamics.R
# Kept for validation reference only. Do not use for new analyses.
# ==============================================================================
# Stratify ALL VF metrics by participant timepoint depth:
#   ≥2, ≥3, ≥4 timepoints
#
# For each cohort: burden, prevalence, enrichment, longitudinal transitions
# ==============================================================================

suppressPackageStartupMessages({
    library(dplyr)
    library(readr)
    library(tidyr)
    library(stringr)
    library(purrr)
})

cat("=== compute_vf_stratified_by_depth.R ===\n")
cat("Started:", format(Sys.time()), "\n\n")

# ---- Load analysis-ready dataset ----
vf_ready <- read_csv("results/vf/vf_analysis_ready.csv", show_col_types = FALSE) %>%
    mutate(Participant_id = as.character(Participant_id))

gene_map <- read_csv("results/vf/gene_map.csv", show_col_types = FALSE) %>%
    mutate(Gene = as.character(Gene), Category = coalesce(as.character(Category), "Unassigned"))

# Identify gene columns
meta_cols <- c("Participant_id", "tp_lab", "Infection_Status", "vf_count_total", "is_ecoli")
cat_cols <- grep("^cat_", names(vf_ready), value = TRUE)
gene_cols <- setdiff(names(vf_ready), c(meta_cols, cat_cols))

cat("Loaded vf_analysis_ready.csv:", nrow(vf_ready), "rows,", n_distinct(vf_ready$Participant_id), "participants\n")
cat("Gene cols:", length(gene_cols), "\n\n")

# ---- Determine participant timepoint counts ----
# Count clinical timepoints per participant (in the VF-available data)
tp_order <- c("T0", "T1", "T2", "T3", "T4", "Uricult")

participant_tp_counts <- vf_ready %>%
    filter(!is.na(Infection_Status)) %>%
    group_by(Participant_id) %>%
    summarise(n_tp = n_distinct(tp_lab), .groups = "drop")

cat("Participant timepoint distribution:\n")
print(table(participant_tp_counts$n_tp))
cat("\n")

# ---- Master function: compute all metrics for a cohort ----
compute_cohort_metrics <- function(min_tp, vf_data, gene_cols, gene_map, tp_order) {
    cohort_label <- paste0(">=", min_tp, " timepoints")
    cat("\n============================================================\n")
    cat("COHORT:", cohort_label, "\n")
    cat("============================================================\n\n")

    # Filter to participants with >= min_tp timepoints
    eligible <- participant_tp_counts %>% filter(n_tp >= min_tp)
    cohort <- vf_data %>%
        filter(Participant_id %in% eligible$Participant_id, !is.na(Infection_Status))

    cat("Participants:", nrow(eligible), "\n")
    cat("Episodes:", nrow(cohort), "\n")
    cat("Status breakdown:\n")
    print(table(cohort$Infection_Status))
    cat("\n")

    # ---- C1: VF Burden ----
    burden <- cohort %>%
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
        ) %>%
        mutate(cohort = cohort_label)

    cat("VF BURDEN:\n")
    for (i in seq_len(nrow(burden))) {
        r <- burden[i, ]
        cat(sprintf(
            "  %s: n=%d, median=%.0f (IQR %.0f-%.0f), mean=%.1f (SD=%.1f)\n",
            r$Infection_Status, r$n_episodes, r$median_vf,
            r$q25_vf, r$q75_vf, r$mean_vf, r$sd_vf
        ))
    }
    cat("\n")

    # ---- C2: Gene Prevalence ----
    gene_prev_list <- lapply(gene_cols, function(g) {
        df <- cohort %>%
            group_by(Infection_Status) %>%
            summarise(
                n_present = sum(.data[[g]] > 0, na.rm = TRUE),
                n_total = n(), .groups = "drop"
            ) %>%
            mutate(pct = round(n_present / n_total * 100, 1))

        asb <- df %>% filter(Infection_Status == "ASB")
        uti <- df %>% filter(Infection_Status == "UTI")
        neg <- df %>% filter(Infection_Status == "Negative")

        tibble(
            gene = g,
            ASB_n = if (nrow(asb)) asb$n_present else NA_integer_,
            ASB_N = if (nrow(asb)) asb$n_total else NA_integer_,
            ASB_pct = if (nrow(asb)) asb$pct else NA_real_,
            UTI_n = if (nrow(uti)) uti$n_present else NA_integer_,
            UTI_N = if (nrow(uti)) uti$n_total else NA_integer_,
            UTI_pct = if (nrow(uti)) uti$pct else NA_real_,
            Neg_n = if (nrow(neg)) neg$n_present else NA_integer_,
            Neg_N = if (nrow(neg)) neg$n_total else NA_integer_,
            Neg_pct = if (nrow(neg)) neg$pct else NA_real_
        )
    })

    gene_prev <- bind_rows(gene_prev_list) %>%
        mutate(delta_UTI_minus_ASB = round(UTI_pct - ASB_pct, 1)) %>%
        left_join(gene_map %>% select(Gene, Category), by = c("gene" = "Gene")) %>%
        mutate(
            Category = coalesce(Category, "Unassigned"),
            cohort = cohort_label
        ) %>%
        arrange(desc(abs(delta_UTI_minus_ASB)))

    cat("TOP 10 GENES BY |Δ(UTI−ASB)|:\n")
    top10 <- head(gene_prev, 10)
    for (i in seq_len(nrow(top10))) {
        r <- top10[i, ]
        cat(sprintf(
            "  %s (%s): ASB=%.1f%%, UTI=%.1f%%, Δ=%+.1f pp\n",
            r$gene, r$Category, r$ASB_pct,
            ifelse(is.na(r$UTI_pct), 0, r$UTI_pct),
            ifelse(is.na(r$delta_UTI_minus_ASB), 0, r$delta_UTI_minus_ASB)
        ))
    }
    cat("\n")

    # ---- C3: Enrichment (Fisher, UTI vs ASB) ----
    asb_uti <- cohort %>% filter(Infection_Status %in% c("ASB", "UTI"))
    n_asb <- sum(asb_uti$Infection_Status == "ASB")
    n_uti <- sum(asb_uti$Infection_Status == "UTI")

    enrichment <- NULL
    if (n_uti >= 2) {
        enr_list <- lapply(gene_cols, function(g) {
            tab <- table(
                gene_present = asb_uti[[g]] > 0,
                is_uti = asb_uti$Infection_Status == "UTI"
            )
            if (nrow(tab) < 2 || ncol(tab) < 2) {
                return(NULL)
            }
            ft <- tryCatch(fisher.test(tab), error = function(e) NULL)
            if (is.null(ft)) {
                return(NULL)
            }
            tibble(
                gene = g, OR = round(ft$estimate, 3),
                p_value = round(ft$p.value, 4),
                CI_lower = round(ft$conf.int[1], 3),
                CI_upper = round(ft$conf.int[2], 3)
            )
        })
        enrichment <- bind_rows(enr_list) %>%
            mutate(
                p_adj_BH = round(p.adjust(p_value, method = "BH"), 4),
                cohort = cohort_label
            ) %>%
            left_join(gene_map %>% select(Gene, Category), by = c("gene" = "Gene")) %>%
            mutate(Category = coalesce(Category, "Unassigned")) %>%
            arrange(p_value)

        cat("ENRICHMENT (Fisher, UTI vs ASB, n_ASB=", n_asb, ", n_UTI=", n_uti, "):\n")
        cat("  Genes with p < 0.05:", sum(enrichment$p_value < 0.05, na.rm = TRUE), "\n")
        cat("  Genes with p_adj < 0.05:", sum(enrichment$p_adj_BH < 0.05, na.rm = TRUE), "\n")
        cat("  Top 5:\n")
        print(head(as.data.frame(enrichment %>% select(gene, Category, OR, p_value, p_adj_BH)), 5))
        cat("\n")
    } else {
        cat("ENRICHMENT: Skipped (only", n_uti, "UTI episodes)\n\n")
    }

    # ---- D: Longitudinal transitions ----
    eligible_long <- participant_tp_counts %>% filter(n_tp >= max(min_tp, 2))
    long_data <- cohort %>%
        filter(
            Participant_id %in% eligible_long$Participant_id,
            tp_lab %in% tp_order
        ) %>%
        mutate(tp_rank = match(tp_lab, tp_order)) %>%
        arrange(Participant_id, tp_rank)

    transitions <- list()
    for (pid in unique(long_data$Participant_id)) {
        pid_data <- long_data %>%
            filter(Participant_id == pid) %>%
            arrange(tp_rank)
        if (nrow(pid_data) < 2) next
        for (i in 1:(nrow(pid_data) - 1)) {
            rf <- pid_data[i, ]
            rt <- pid_data[i + 1, ]
            vf_from <- gene_cols[which(rf[, gene_cols] > 0)]
            vf_to <- gene_cols[which(rt[, gene_cols] > 0)]
            gained <- setdiff(vf_to, vf_from)
            lost <- setdiff(vf_from, vf_to)
            stable <- intersect(vf_from, vf_to)
            union_s <- union(vf_from, vf_to)
            jac <- if (length(union_s) == 0) NA_real_ else length(stable) / length(union_s)
            transitions[[length(transitions) + 1]] <- tibble(
                Participant_id = pid,
                tp_from = as.character(rf$tp_lab), tp_to = as.character(rt$tp_lab),
                status_from = rf$Infection_Status, status_to = rt$Infection_Status,
                transition_type = paste0(rf$Infection_Status, "→", rt$Infection_Status),
                vf_count_from = rf$vf_count_total, vf_count_to = rt$vf_count_total,
                n_gained = length(gained), n_lost = length(lost),
                n_stable = length(stable), jaccard = round(jac, 3),
                any_change = length(gained) > 0 | length(lost) > 0,
                genes_gained = paste(sort(gained), collapse = ";"),
                genes_lost = paste(sort(lost), collapse = ";")
            )
        }
    }
    trans_df <- bind_rows(transitions) %>% mutate(cohort = cohort_label)

    if (nrow(trans_df) > 0) {
        trans_summ <- trans_df %>%
            group_by(transition_type) %>%
            summarise(
                n_transitions = n(), n_participants = n_distinct(Participant_id),
                median_gained = median(n_gained), mean_gained = round(mean(n_gained), 1),
                median_lost = median(n_lost), mean_lost = round(mean(n_lost), 1),
                median_stable = median(n_stable),
                median_jaccard = round(median(jaccard, na.rm = TRUE), 3),
                mean_jaccard = round(mean(jaccard, na.rm = TRUE), 3),
                pct_no_change = round(mean(!any_change) * 100, 1),
                .groups = "drop"
            ) %>%
            mutate(cohort = cohort_label) %>%
            arrange(desc(n_transitions))

        cat("LONGITUDINAL TRANSITIONS:\n")
        cat("  Total transitions:", nrow(trans_df), "\n")
        cat("  Participants:", n_distinct(trans_df$Participant_id), "\n")
        cat("  % zero change:", round(mean(!trans_df$any_change) * 100, 1), "%\n")
        cat("  Median Jaccard:", round(median(trans_df$jaccard, na.rm = TRUE), 3), "\n")
        cat("  Transition types:\n")
        print(as.data.frame(trans_summ %>% select(transition_type, n_transitions, median_jaccard, pct_no_change)))

        # ASB→UTI gains
        a2u <- trans_df %>% filter(transition_type == "ASB→UTI")
        if (nrow(a2u) > 0) {
            cat("\n  ASB→UTI transitions:", nrow(a2u), "\n")
            if (any(a2u$n_gained > 0)) {
                top_gained <- a2u %>%
                    filter(genes_gained != "") %>%
                    pull(genes_gained) %>%
                    strsplit(";") %>%
                    unlist() %>%
                    table() %>%
                    sort(decreasing = TRUE)
                cat("  Top gained genes:\n")
                print(head(top_gained, 8))
            }
        }
    } else {
        trans_summ <- tibble()
        cat("LONGITUDINAL: No transitions (participants have only 1 timepoint in cohort)\n")
    }

    cat("\n")

    return(list(
        burden = burden,
        gene_prev = gene_prev,
        enrichment = enrichment,
        transitions = trans_df,
        trans_summary = trans_summ
    ))
}

# ---- Run for each cohort ----
results <- list()
for (min_tp in c(2, 3, 4)) {
    results[[paste0("tp", min_tp)]] <- compute_cohort_metrics(
        min_tp = min_tp,
        vf_data = vf_ready,
        gene_cols = gene_cols,
        gene_map = gene_map,
        tp_order = tp_order
    )
}

# ---- Combine and write outputs ----
cat("\n=== WRITING COMBINED OUTPUTS ===\n\n")

# Burden
all_burden <- bind_rows(lapply(results, `[[`, "burden"))
write_csv(all_burden, "results/vf/vf_burden_by_status_stratified.csv")
cat("Written: vf_burden_by_status_stratified.csv\n")

# Prevalence
all_prev <- bind_rows(lapply(results, `[[`, "gene_prev"))
write_csv(all_prev, "results/vf/vf_gene_prevalence_stratified.csv")
cat("Written: vf_gene_prevalence_stratified.csv\n")

# Enrichment
all_enr <- bind_rows(compact(lapply(results, `[[`, "enrichment")))
if (nrow(all_enr) > 0) {
    write_csv(all_enr, "results/vf/vf_enrichment_stratified.csv")
    cat("Written: vf_enrichment_stratified.csv\n")
}

# Transitions
all_trans <- bind_rows(lapply(results, `[[`, "transitions"))
write_csv(all_trans, "results/vf/vf_transitions_stratified.csv")
cat("Written: vf_transitions_stratified.csv\n")

# Transition summaries
all_tsumm <- bind_rows(lapply(results, `[[`, "trans_summary"))
write_csv(all_tsumm, "results/vf/vf_transition_summary_stratified.csv")
cat("Written: vf_transition_summary_stratified.csv\n")

# ---- Print side-by-side comparison table ----
cat("\n\n================================================================\n")
cat("SIDE-BY-SIDE COMPARISON: VF BURDEN BY STATUS × COHORT\n")
cat("================================================================\n\n")

comparison <- all_burden %>%
    mutate(label = sprintf(
        "%d episodes, median=%.0f (IQR %.0f-%.0f), mean=%.1f±%.1f",
        n_episodes, median_vf, q25_vf, q75_vf, mean_vf, sd_vf
    )) %>%
    select(cohort, Infection_Status, label) %>%
    pivot_wider(names_from = cohort, values_from = label)

print(as.data.frame(comparison))

cat("\n\n================================================================\n")
cat("SIDE-BY-SIDE COMPARISON: LONGITUDINAL SUMMARY × COHORT\n")
cat("================================================================\n\n")

if (nrow(all_tsumm) > 0) {
    long_compare <- all_tsumm %>%
        mutate(label = sprintf(
            "n=%d, Jac=%.3f, %%noΔ=%.0f%%",
            n_transitions, median_jaccard, pct_no_change
        )) %>%
        select(cohort, transition_type, label) %>%
        pivot_wider(names_from = cohort, values_from = label, values_fill = "—")
    print(as.data.frame(long_compare))
}

cat("\n=== DONE ===\n")
cat("Finished:", format(Sys.time()), "\n")

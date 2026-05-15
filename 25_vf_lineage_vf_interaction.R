#!/usr/bin/env Rscript
# ==============================================================================
# 25_vf_lineage_vf_interaction.R
# ==============================================================================
#
# GOAL:
#   Determine whether observed VF burden differences between ASB and UTI
#   are driven by clinical status or by the underlying bacterial lineage (ST).
#
# WHY THIS SCRIPT IS CRITICAL:
#   This is the study's most important confounding check.  Consider this
#   scenario:
#     - ST131 naturally carries 90 VF genes.
#     - ST73 naturally carries 70 VF genes.
#     - If ST131 happens to cause more UTIs in our cohort, then a naive
#       comparison would show "UTI episodes have more VFs" — but this is
#       entirely a lineage artefact, not a genuine VF→UTI association.
#
#   Without this analysis, any VF-status association claim is suspect.
#
#   The script addresses three sequential questions:
#     Q1: Does VF burden vary across STs?
#         (Kruskal-Wallis test)
#     Q2: Within each major ST, is there a VF burden difference between
#         ASB and UTI?
#         (Within-ST Wilcoxon tests — small sample sizes expected)
#     Q3: Is ST composition different between ASB and UTI?
#         (Fisher exact on ST×Status contingency table)
#
#   Interpretation:
#     Q1=YES + Q3=YES → Lineage IS a confounder.  Add ST as a covariate
#                        in 14_genotype_phenotype_model.R.
#     Q1=YES + Q3=NO  → STs carry different VFs but aren't differentially
#                        distributed.  Less confounding concern.
#     Q1=NO           → Lineage is unlikely to confound VF-status results.
#
# RELATIONSHIP TO EXISTING SCRIPTS:
#   - 17_lineage_analysis.R computes UTI *risk* per ST (what proportion of
#     each ST's episodes are UTI?).  It does NOT examine VF burden per ST.
#     This script fills that gap.
#   - 14_genotype_phenotype_model.R runs GLMM on individual VF genes.
#     If this script detects confounding, the recommendation is to add ST
#     as a covariate there, NOT to build a separate modelling script.
#
# INPUTS:
#   - results/vf/vf_analysis_ready.csv   (from 22_vf_build_analysis_dataset.R)
#   - results/vf/gene_map.csv            (for category-level medians per ST)
#
# OUTPUTS (all in results/vf/):
#   - vf_burden_by_st.csv                  VF burden per ST (top STs only)
#   - vf_burden_by_st_and_status.csv       VF burden per ST × Status
#   - vf_lineage_confounding_summary.txt   Human-readable confounding report
#
# PLOTS (in plots/vf/):
#   - vf_burden_by_st.png                  Boxplot of VF burden per top STs
#   - vf_burden_st_x_status.png            Faceted: ASB vs UTI within each ST
# ==============================================================================

source("00_config.R")
source("R/plot_helpers.R")  # Canonical infection status colours
suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
})

msg("Starting 25_vf_lineage_vf_interaction.R")

# ==============================================================================
# VF visualisation shared helpers
# ==============================================================================

STATUS_LEVELS <- c("ASB", "UTI", "Negative", "Culture-positive/S&S unknown", "Unknown")

status_for_plot <- function(x) {
  x <- as.character(x)
  x[is.na(x) | x == ""] <- "Unknown"
  x[!x %in% STATUS_LEVELS] <- "Culture-positive/S&S unknown"
  factor(x, levels = STATUS_LEVELS)
}

normalise_st_label <- function(x) {
  x <- as.character(x)
  x <- str_trim(x)
  unknown <- c("", "-", "ST-", "NA", "N/A", "UNKNOWN", "UNK", "NT",
               "NON-TYPABLE", "NONTYPABLE", "NOT TYPED")
  x[str_to_upper(x) %in% unknown] <- NA_character_
  x
}

status_count_text <- function(data) {
  counts <- data %>%
    mutate(.status = status_for_plot(Infection_Status)) %>%
    count(.status, name = "n") %>%
    filter(n > 0) %>%
    mutate(label = paste0(as.character(.status), " n=", n))
  paste(counts$label, collapse = "; ")
}

vf_caption <- function(input_file, data, analysis_unit,
                       p_value_note = "Tests and plots are diagnostic/exploratory and do not establish causality.",
                       extra_note = NULL) {
  repeated <- if ("Participant_id" %in% names(data)) any(table(data$Participant_id) > 1) else FALSE
  n_uti <- if ("Infection_Status" %in% names(data)) sum(data$Infection_Status == "UTI", na.rm = TRUE) else NA_integer_
  paste(
    sprintf("Data: %s.", input_file),
    sprintf("Denominator: n=%d VF/WGS-linked E. coli isolates from %d participants%s.",
            nrow(data),
            if ("Participant_id" %in% names(data)) n_distinct(data$Participant_id) else NA_integer_,
            if ("Infection_Status" %in% names(data)) paste0(" (", status_count_text(data), ")") else ""),
    sprintf("Level of analysis: %s.", analysis_unit),
    if (repeated) "Residents may contribute repeated isolates; isolate-level summaries can be pseudoreplicated." else "",
    p_value_note,
    sprintf("UTI denominator is small (n=%d).", n_uti),
    "Lineage/ST, batch, timepoint, and event-driven sampling structure should be considered when interpreting VF associations.",
    extra_note %||% "",
    sep = " "
  ) %>% str_squish()
}

plot_theme_vf <- function(base_size = 11) {
  theme_bw(base_size = base_size) +
    theme(
      plot.caption = element_text(hjust = 0, size = base_size - 3, colour = "grey35"),
      plot.subtitle = element_text(colour = "grey25"),
      legend.position = "bottom"
    )
}

timepoint_display_order <- function(x) {
  ux <- unique(as.character(x))
  ux <- ux[!is.na(ux) & ux != ""]
  routine <- ux[str_detect(str_to_upper(ux), "^T\\d+$")]
  routine <- routine[order(suppressWarnings(as.numeric(str_remove(str_to_upper(routine), "^T"))))]
  uti <- ux[str_detect(str_to_upper(ux), "^UTI-\\d+$")]
  uti <- uti[order(suppressWarnings(as.numeric(str_remove(str_to_upper(uti), "^UTI-"))))]
  uricult <- ux[str_detect(ux, regex("uricult", ignore_case = TRUE))]
  other <- setdiff(ux, c(routine, uricult, uti))
  c(routine, uricult, uti, sort(other))
}

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
      "%s is older than %s. Re-run 22_vf_build_analysis_dataset.R before script 25 so confounding plots use the current VF matrix.\n  %s: %s\n  %s: %s",
      target_label, upstream_label, target_label, format(target_mtime),
      upstream_label, format(upstream_mtime)
    ))
  }
}

# ==============================================================================
# 1. LOAD DATA
# ==============================================================================

ready_file <- FILE_VF_READY
if (!file.exists(ready_file)) stop("Missing ", ready_file, ". Run 22_vf_build_analysis_dataset.R first.")
stop_if_stale(ready_file, FILE_VF_PA, "vf_analysis_ready.csv", "vf_pa_all.csv")
vf_ready <- read_csv(ready_file, show_col_types = FALSE) %>%
  mutate(Participant_id = as.character(Participant_id),
         ST = normalise_st_label(ST),
         Infection_Status = as.character(Infection_Status))

gene_map <- read_csv(file.path(DIR_VF, "gene_map.csv"), show_col_types = FALSE) %>%
  mutate(Gene = as.character(Gene),
         Category = coalesce(as.character(Category), "Unassigned"))

cat_cols <- grep("^cat_", names(vf_ready), value = TRUE)

msg("Loaded: %d rows, %d with ST, %d distinct STs",
    nrow(vf_ready), sum(!is.na(vf_ready$ST)), n_distinct(vf_ready$ST, na.rm = TRUE))
msg("VF-ready denominator by status: %s", status_count_text(vf_ready))

# ==============================================================================
# 2. CHECK ST DATA AVAILABILITY
# ==============================================================================
# If MLST data was not joined in 22_ (e.g., mlst_with_meta.csv was missing),
# this script cannot run.  Fail gracefully with a clear message.

has_st <- vf_ready %>% filter(!is.na(ST), !is.na(Infection_Status))

if (nrow(has_st) < 10) {
  msg("WARNING: Only %d episodes have both ST and Infection_Status.", nrow(has_st))
  msg("Lineage confounding analysis will be limited.")
}

if (nrow(has_st) == 0) {
  msg("ERROR: No episodes with ST data.  Cannot run lineage analysis.")
  msg("Ensure 22_vf_build_analysis_dataset.R successfully joined MLST data.")
  msg("Writing empty outputs and exiting.")

  writeLines("No ST data available for lineage–VF analysis.",
             file.path(DIR_VF, "vf_lineage_confounding_summary.txt"))
  quit(save = "no")
}

# ==============================================================================
# 3. VF BURDEN BY SEQUENCE TYPE
# ==============================================================================
# We focus on STs with sufficient data (≥5 episodes) to make meaningful
# comparisons.  With <5 episodes, a single outlier dominates the summary.

st_counts <- has_st %>%
  count(ST, name = "n_episodes") %>%
  arrange(desc(n_episodes))

min_episodes <- 5
top_sts <- st_counts %>% filter(n_episodes >= min_episodes) %>% pull(ST)

msg("STs with ≥%d episodes: %d (out of %d total)",
    min_episodes, length(top_sts), nrow(st_counts))

# If no ST has ≥5 episodes, relax the threshold to ≥3 to still produce output
if (length(top_sts) == 0) {
  msg("WARNING: No ST has ≥%d episodes.  Relaxing threshold to ≥3.", min_episodes)
  min_episodes <- 3
  top_sts <- st_counts %>% filter(n_episodes >= min_episodes) %>% pull(ST)
}

# Compute VF burden descriptive statistics for each top ST.
# This table directly answers Q1: do different STs carry different VF loads?
burden_by_st <- has_st %>%
  filter(ST %in% top_sts) %>%
  group_by(ST) %>%
  summarise(
    n_episodes     = n(),
    n_participants = n_distinct(Participant_id),
    mean_vf   = round(mean(vf_count_total), 1),
    sd_vf     = round(sd(vf_count_total), 1),
    median_vf = median(vf_count_total),
    q25_vf    = quantile(vf_count_total, 0.25),
    q75_vf    = quantile(vf_count_total, 0.75),
    .groups   = "drop"
  ) %>%
  arrange(desc(median_vf))

# Also add category-level medians per ST (e.g., how many Iron Acquisition
# genes does ST131 typically carry vs ST73?)
for (cc in cat_cols) {
  cat_med <- has_st %>%
    filter(ST %in% top_sts) %>%
    group_by(ST) %>%
    summarise(!!paste0(cc, "_median") := median(.data[[cc]], na.rm = TRUE),
              .groups = "drop")
  burden_by_st <- burden_by_st %>% left_join(cat_med, by = "ST")
}

write_csv(burden_by_st, file.path(DIR_VF, "vf_burden_by_st.csv"))

# Compute VF burden per ST × Status (ASB/UTI only).
# This answers Q2: within the same ST, does VF burden differ by status?
asb_uti_st <- has_st %>%
  filter(ST %in% top_sts, Infection_Status %in% c("ASB", "UTI"))

burden_st_status <- asb_uti_st %>%
  group_by(ST, Infection_Status) %>%
  summarise(
    n_episodes = n(),
    mean_vf    = round(mean(vf_count_total), 1),
    median_vf  = median(vf_count_total),
    .groups    = "drop"
  ) %>%
  arrange(ST, Infection_Status)

write_csv(burden_st_status, file.path(DIR_VF, "vf_burden_by_st_and_status.csv"))

# ==============================================================================
# 4. STATISTICAL TESTS FOR CONFOUNDING
# ==============================================================================
# We apply the three-question framework described in the header.

summary_lines <- character()
sl <- function(...) summary_lines <<- c(summary_lines, sprintf(...))

sl("=== VF–LINEAGE CONFOUNDING ASSESSMENT ===")
sl("Generated: %s", format(Sys.time()))
sl("")

# ------------------------------------------------------------------
# Q1: Does VF burden differ significantly across STs?
# ------------------------------------------------------------------
# Kruskal-Wallis is a non-parametric test for differences in medians
# across >2 groups (distribution-free alternative to one-way ANOVA).
# If significant, it means STs carry genuinely different VF loads,
# which is a prerequisite for lineage confounding.
if (length(top_sts) >= 2) {
  kw_data <- has_st %>% filter(ST %in% top_sts)
  kw <- kruskal.test(vf_count_total ~ factor(ST), data = kw_data)
  sl("Q1: Does VF burden differ across STs?")
  sl("  Kruskal–Wallis chi-sq = %.2f, df = %d, p = %.4f",
     kw$statistic, kw$parameter, kw$p.value)
  if (kw$p.value < 0.05) {
    sl("  → YES: Significant VF burden variation across STs.")
    sl("  → This means lineage is a potential confounder.")
  } else {
    sl("  → NO: VF burden does not significantly vary across STs.")
  }
  sl("")
}

# ------------------------------------------------------------------
# Q2: Within each major ST, does VF burden differ between ASB and UTI?
# ------------------------------------------------------------------
# This is the key within-lineage test.  If VF burden is the SAME for
# ASB and UTI episodes of the same ST, then VF differences observed
# in the overall cohort are likely driven by ST composition, not by
# VF content per se.
#
# NOTE: Small sample sizes within STs are expected.  Many STs will have
# too few UTI episodes for a valid test.  This is a known limitation.
sl("Q2: Within each top ST, does VF burden differ between ASB and UTI?")
for (st in top_sts) {
  st_data <- asb_uti_st %>% filter(ST == st)
  n_asb <- sum(st_data$Infection_Status == "ASB")
  n_uti <- sum(st_data$Infection_Status == "UTI")

  if (n_asb >= 2 && n_uti >= 2) {
    # Wilcoxon rank-sum (Mann-Whitney U) test: non-parametric comparison
    # of two independent groups.
    wt <- tryCatch(
      wilcox.test(vf_count_total ~ Infection_Status, data = st_data),
      error = function(e) NULL
    )
    if (!is.null(wt)) {
      sl("  ST%s (ASB=%d, UTI=%d): Wilcoxon p = %.4f, median ASB=%.0f, UTI=%.0f",
         st, n_asb, n_uti, wt$p.value,
         median(st_data$vf_count_total[st_data$Infection_Status == "ASB"]),
         median(st_data$vf_count_total[st_data$Infection_Status == "UTI"]))
    }
  } else {
    sl("  ST%s (ASB=%d, UTI=%d): Too few in one group for within-ST test.", st, n_asb, n_uti)
  }
}

# ------------------------------------------------------------------
# Q3: Does ST composition differ between ASB and UTI?
# ------------------------------------------------------------------
# This tests whether certain STs are over-represented in UTI vs ASB.
# If yes, and Q1 is also yes, then lineage is a LIKELY CONFOUNDER:
# UTI episodes might have higher VF counts simply because they tend
# to harbour VF-heavy STs, not because VFs cause UTI.
sl("")
sl("Q3: Does ST composition differ between ASB and UTI?")
st_status_tab <- asb_uti_st %>%
  filter(ST %in% top_sts) %>%
  count(ST, Infection_Status) %>%
  pivot_wider(names_from = Infection_Status, values_from = n, values_fill = 0)

if (nrow(st_status_tab) >= 2 && ncol(st_status_tab) >= 3) {
  mat <- as.matrix(st_status_tab %>% select(-ST))
  # Use simulated p-values because the contingency table may be sparse
  ft <- tryCatch(fisher.test(mat, simulate.p.value = TRUE, B = 10000),
                 error = function(e) NULL)
  if (!is.null(ft)) {
    sl("  Fisher exact (simulated): p = %.4f", ft$p.value)
    if (ft$p.value < 0.05) {
      sl("  → YES: ST composition differs between ASB and UTI.")
      sl("  → Combined with Q1, this means lineage is a LIKELY CONFOUNDER.")
    } else {
      sl("  → NO: ST composition is not significantly different.")
    }
  }
}

# ------------------------------------------------------------------
# Interpretation guide
# ------------------------------------------------------------------
sl("")
sl("INTERPRETATION GUIDANCE:")
sl("  If Q1=YES and Q3=YES → VF–status differences may be ST-driven artefacts.")
sl("  → Add ST as covariate in 14_genotype_phenotype_model.R.")
sl("  If Q1=YES and Q3=NO  → STs carry different VFs but aren't over-represented")
sl("    in either status.  Less confounding concern.")
sl("  If Q1=NO  → Lineage is unlikely to confound VF–status results.")

# ==============================================================================
# 5. PLOTS
# ==============================================================================

ensure_dir(DIR_PLOTS_VF)

# PLOT 1: VF burden boxplot by ST
#   Shows the range of VF gene counts for each major ST.
#   Large differences here (Q1=YES) mean that lineage choice alone
#   substantially determines VF load.
if (length(top_sts) >= 2) {
  plot_st <- has_st %>%
    filter(ST %in% top_sts) %>%
    mutate(ST_label = paste0("ST", ST))

  p_st <- ggplot(plot_st, aes(x = reorder(ST_label, -vf_count_total, FUN = median),
                                y = vf_count_total)) +
    geom_boxplot(outlier.shape = NA, width = 0.5, fill = "grey90") +
    geom_jitter(width = 0.15, alpha = 0.5, size = 1.5) +
    labs(
      title = "Virulence factor burden varies by E. coli sequence type",
      subtitle = sprintf("STs with at least %d VF-ready isolates; VF burden = detected VF genes per isolate", min_episodes),
      x = "Sequence type",
      y = "Detected VF genes per isolate",
      caption = vf_caption(
        ready_file, plot_st, "isolate-level diagnostic by ST",
        p_value_note = "Kruskal-Wallis testing in the summary file is exploratory and does not account for repeated resident isolates.",
        extra_note = "Sparse STs are omitted from this burden plot to avoid unstable visual summaries."
      )
    ) +
    plot_theme_vf(base_size = 11) +
    theme(axis.text.x = element_text(angle = 30, hjust = 1))

  ggsave(file.path(DIR_PLOTS_VF, "vf_burden_by_st.png"), p_st,
         width = max(7, length(top_sts) * 0.8), height = 5.6, dpi = 300)
  ggsave(file.path(DIR_PLOTS_VF, "vf_burden_by_top_st.png"), p_st,
         width = max(7, length(top_sts) * 0.8), height = 5.6, dpi = 300)
}

# PLOT 2: VF burden by ST × Status (ASB vs UTI, faceted)
#   This is the visual complement to the Q2 tests above.
#   For each ST that has both ASB and UTI episodes, show side-by-side
#   boxplots.  If boxes overlap heavily within each ST, the VF-status
#   association is likely an artefact of ST composition.
if (nrow(asb_uti_st) > 0 && length(top_sts) >= 2) {
  plot_st_status <- asb_uti_st %>%
    mutate(ST_label = paste0("ST", ST))

  # Only facet STs that have at least 1 episode of each status
  sts_both <- plot_st_status %>%
    count(ST_label, Infection_Status) %>%
    group_by(ST_label) %>%
    filter(n() >= 2) %>%
    pull(ST_label) %>% unique()

  if (length(sts_both) >= 1) {
    p_stx <- ggplot(plot_st_status %>% filter(ST_label %in% sts_both),
                    aes(x = Infection_Status, y = vf_count_total,
                        fill = Infection_Status)) +
      geom_boxplot(outlier.shape = NA, width = 0.5, alpha = 0.7) +
      geom_jitter(width = 0.15, alpha = 0.4, size = 1.2) +
      scale_fill_infection() +
      facet_wrap(~ST_label, scales = "free_x") +
      labs(
        title = "ASB-UTI VF burden contrasts within E. coli sequence types",
        subtitle = "Only STs with at least one ASB and one UTI VF-ready isolate are shown",
        x = "Clinical status",
        y = "Detected VF genes per isolate",
        caption = vf_caption(
          ready_file, plot_st_status %>% filter(ST_label %in% sts_both),
          "isolate-level within-ST ASB-UTI diagnostic",
          p_value_note = "Within-ST Wilcoxon tests in the summary file are exploratory and often underpowered.",
          extra_note = "This plot is a confounding diagnostic; overlapping distributions argue against interpreting naive ASB-UTI contrasts as causal."
        )
      ) +
      plot_theme_vf(base_size = 10) +
      theme(legend.position = "none",
            strip.text = element_text(face = "bold"))

    ggsave(file.path(DIR_PLOTS_VF, "vf_burden_st_x_status.png"), p_stx,
           width = max(6, length(sts_both) * 2.5),
           height = max(4, ceiling(length(sts_both) / 3) * 3),
           dpi = 300)
  }
}

# =============================================================================
# VF visualisation module 05: Confounding diagnostics
# =============================================================================

status_st_data <- vf_ready %>%
  filter(!is.na(Infection_Status), Infection_Status %in% c("ASB", "UTI", "Negative")) %>%
  mutate(
    Infection_Status = status_for_plot(Infection_Status),
    ST_group = ifelse(is.na(ST), "Missing/non-typable ST", paste0("ST", ST))
  )

if (nrow(status_st_data) > 0) {
  common_st_groups <- status_st_data %>%
    count(ST_group, sort = TRUE) %>%
    filter(ST_group != "Missing/non-typable ST") %>%
    slice_head(n = 12) %>%
    pull(ST_group)

  st_comp <- status_st_data %>%
    mutate(ST_group = case_when(
      ST_group %in% common_st_groups ~ ST_group,
      ST_group == "Missing/non-typable ST" ~ ST_group,
      TRUE ~ "Other STs"
    )) %>%
    count(Infection_Status, ST_group, name = "n") %>%
    group_by(Infection_Status) %>%
    mutate(prop = n / sum(n)) %>%
    ungroup()

  p_st_comp <- ggplot(st_comp, aes(x = Infection_Status, y = prop, fill = ST_group)) +
    geom_col(width = 0.65, colour = "white", linewidth = 0.2) +
    scale_y_continuous(labels = scales::percent) +
    labs(
      title = "Clinical status distribution across E. coli sequence types",
      subtitle = "Top STs are shown individually; sparse STs are grouped as Other STs",
      x = "Clinical status",
      y = "Within-status proportion of VF-ready isolates",
      fill = "Sequence type",
      caption = vf_caption(
        ready_file, status_st_data, "isolate-level ST composition diagnostic",
        p_value_note = "ST composition tests are exploratory diagnostics and do not establish that lineage causes clinical status.",
        extra_note = "Missing/non-typable STs are shown separately from Other STs."
      )
    ) +
    plot_theme_vf(base_size = 11)

  ggsave(file.path(DIR_PLOTS_VF, "vf_st_composition_by_status.png"),
         p_st_comp, width = 8.5, height = 5.6, dpi = 300)
}

if ("Batch" %in% names(vf_ready)) {
  batch_data <- vf_ready %>%
    filter(!is.na(Infection_Status), Infection_Status %in% c("ASB", "UTI", "Negative")) %>%
    mutate(
      Infection_Status = status_for_plot(Infection_Status),
      Batch = ifelse(is.na(Batch) | Batch == "", "Missing batch", paste0("Batch ", Batch))
    ) %>%
    count(Infection_Status, Batch, name = "n") %>%
    group_by(Infection_Status) %>%
    mutate(prop = n / sum(n)) %>%
    ungroup()

  p_batch <- ggplot(batch_data, aes(x = Infection_Status, y = prop, fill = Batch)) +
    geom_col(width = 0.65, colour = "white", linewidth = 0.2) +
    scale_y_continuous(labels = scales::percent) +
    labs(
      title = "Batch structure across VF-ready clinical states",
      subtitle = "Batch imbalance should be considered before interpreting ASB-UTI VF contrasts",
      x = "Clinical status",
      y = "Within-status proportion of isolates",
      fill = "Batch",
      caption = vf_caption(
        ready_file, vf_ready %>% filter(!is.na(Infection_Status)),
        "isolate-level batch composition diagnostic",
        extra_note = "This plot checks study design and processing structure, not biological causation."
      )
    ) +
    plot_theme_vf(base_size = 11)

  ggsave(file.path(DIR_PLOTS_VF, "vf_batch_by_status.png"),
         p_batch, width = 7.5, height = 5.4, dpi = 300)
} else {
  msg("Skipping batch-by-status diagnostic: Batch column is absent.")
}

if ("Event_type" %in% names(vf_ready)) {
  event_data <- vf_ready %>%
    filter(!is.na(Infection_Status), Infection_Status %in% c("ASB", "UTI", "Negative")) %>%
    mutate(
      Infection_Status = status_for_plot(Infection_Status),
      Event_type = ifelse(is.na(Event_type) | Event_type == "", "Missing event type", Event_type)
    ) %>%
    count(Infection_Status, Event_type, name = "n") %>%
    group_by(Infection_Status) %>%
    mutate(prop = n / sum(n)) %>%
    ungroup()

  p_event <- ggplot(event_data, aes(x = Infection_Status, y = prop, fill = Event_type)) +
    geom_col(width = 0.65, colour = "white", linewidth = 0.2) +
    scale_y_continuous(labels = scales::percent) +
    labs(
      title = "Sampling context across VF-ready clinical states",
      subtitle = "Routine and suspected-UTI event sampling are not interchangeable denominators",
      x = "Clinical status",
      y = "Within-status proportion of isolates",
      fill = "Event type",
      caption = vf_caption(
        ready_file, vf_ready %>% filter(!is.na(Infection_Status)),
        "isolate-level event-type composition diagnostic",
        extra_note = "Uricult/UTI-labelled WGS rows are event-driven sampling contexts and should not be interpreted as routine timepoints."
      )
    ) +
    plot_theme_vf(base_size = 11)

  ggsave(file.path(DIR_PLOTS_VF, "vf_event_type_by_status.png"),
         p_event, width = 7.5, height = 5.4, dpi = 300)
} else {
  msg("Skipping event-type-by-status diagnostic: Event_type column is absent.")
}

if (all(c("tp_lab", "Event_type", "Infection_Status") %in% names(vf_ready))) {
  tp_levels <- timepoint_display_order(vf_ready$tp_lab)
  status_tp_event <- vf_ready %>%
    filter(!is.na(Infection_Status), Infection_Status %in% c("ASB", "UTI", "Negative")) %>%
    mutate(
      Infection_Status = status_for_plot(Infection_Status),
      tp_lab = factor(ifelse(is.na(tp_lab) | tp_lab == "", "Missing timepoint", tp_lab),
                      levels = c(tp_levels, "Missing timepoint")),
      Event_type = ifelse(is.na(Event_type) | Event_type == "", "Missing event type", Event_type)
    ) %>%
    count(Event_type, Infection_Status, tp_lab, name = "n") %>%
    filter(n > 0)

  if (nrow(status_tp_event) > 0) {
    p_tp_event <- ggplot(status_tp_event, aes(x = tp_lab, y = Infection_Status, fill = n)) +
      geom_tile(colour = "white", linewidth = 0.35) +
      geom_text(aes(label = n), size = 3, colour = "grey10") +
      facet_wrap(~Event_type, ncol = 1, scales = "free_x") +
      scale_fill_gradient(low = "white", high = "#0072B2") +
      labs(
        title = "Clinical status, timepoint, and event context in VF-ready isolates",
        subtitle = "UTI-labelled VF/WGS rows are concentrated in event-driven timepoint labels rather than routine sampling labels",
        x = "Categorical timepoint / WGS event label",
        y = "Clinical status",
        fill = "n isolates",
        caption = vf_caption(
          ready_file,
          vf_ready %>% filter(!is.na(Infection_Status)),
          "isolate-level timepoint and event-context diagnostic",
          extra_note = "Timepoint labels are ordered for display only; Uricult/UTI-N labels are not treated as numeric timepoints."
        )
      ) +
      plot_theme_vf(base_size = 11) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))

    ggsave(file.path(DIR_PLOTS_VF, "vf_status_timepoint_event_tile.png"),
           p_tp_event, width = 9.5, height = 6.5, dpi = 300)
  }
} else {
  msg("Skipping status-timepoint-event diagnostic: tp_lab, Event_type, or Infection_Status column is absent.")
}

qc_bias_file <- file.path(DIR_QC, "qc_selection_bias_by_status.csv")
if (file.exists(qc_bias_file)) {
  qc_bias <- read_csv(qc_bias_file, show_col_types = FALSE)
  selection_flags <- intersect(
    c("qc_pass", "selected_canonical", "gff_available", "included_in_current_panaroo"),
    names(qc_bias)
  )

  if (length(selection_flags) > 0 && all(c("Infection_Status", "n") %in% names(qc_bias))) {
    qc_rates <- qc_bias %>%
      filter(!is.na(Infection_Status), Infection_Status %in% c("ASB", "UTI", "Negative")) %>%
      mutate(Infection_Status = status_for_plot(Infection_Status)) %>%
      pivot_longer(all_of(selection_flags), names_to = "selection_step", values_to = "included") %>%
      group_by(Infection_Status, selection_step) %>%
      summarise(
        n_included = sum(n[included %in% TRUE], na.rm = TRUE),
        n_total = sum(n, na.rm = TRUE),
        prop = ifelse(n_total > 0, n_included / n_total, NA_real_),
        .groups = "drop"
      ) %>%
      mutate(
        selection_step = recode(selection_step,
                                qc_pass = "QC pass",
                                selected_canonical = "Canonical assembly selected",
                                gff_available = "GFF available",
                                included_in_current_panaroo = "Included in current Panaroo",
                                .default = selection_step),
        selection_step = factor(selection_step, levels = c(
          "QC pass", "Canonical assembly selected", "GFF available", "Included in current Panaroo"
        )),
        label = sprintf("%d/%d", n_included, n_total)
      )

    if (nrow(qc_rates) > 0) {
      p_qc <- ggplot(qc_rates, aes(x = selection_step, y = Infection_Status, fill = prop)) +
        geom_tile(colour = "white", linewidth = 0.35) +
        geom_text(aes(label = label), size = 3.1, colour = "grey10") +
        scale_fill_gradient(labels = scales::percent, low = "white", high = "#009E73",
                            na.value = "grey90", limits = c(0, 1)) +
        labs(
          title = "WGS/QC selection structure by clinical status",
          subtitle = "Status-specific inclusion rates help assess selection bias before interpreting VF associations",
          x = "Selection / availability step",
          y = "Clinical status",
          fill = "Included",
          caption = paste(
            sprintf("Data: %s.", qc_bias_file),
            "Level of analysis: QC and selection-bias diagnostic by clinical status.",
            "The figure shows weighted inclusion counts from the QC selection table, not a biological VF association.",
            "Unequal inclusion across status groups can affect VF-ready denominators and downstream ASB-UTI interpretation."
          )
        ) +
        plot_theme_vf(base_size = 11) +
        theme(axis.text.x = element_text(angle = 25, hjust = 1))

      ggsave(file.path(DIR_PLOTS_VF, "vf_qc_selection_by_status.png"),
             p_qc, width = 8.5, height = 5.4, dpi = 300)
    }
  } else {
    msg("Skipping QC selection diagnostic: required QC selection columns are absent.")
  }
} else {
  msg("Skipping QC selection diagnostic: %s is absent.", qc_bias_file)
}

status_map <- if (file.exists(FILE_STATUS_MAP)) {
  read_csv(FILE_STATUS_MAP, show_col_types = FALSE) %>%
    mutate(Participant_id = as.character(Participant_id))
} else {
  NULL
}
vf_pa <- if (file.exists(FILE_VF_PA)) {
  read_csv(FILE_VF_PA, show_col_types = FALSE) %>%
    mutate(Participant_id = as.character(Participant_id))
} else {
  NULL
}

flow_rows <- bind_rows(
  if (!is.null(status_map)) {
    status_map %>%
      filter(!is.na(Infection_Status)) %>%
      count(stage = "Clinical status map", Infection_Status, name = "n")
  },
  if (!is.null(vf_pa)) {
    tibble(stage = "Raw VF P/A matrix", Infection_Status = "All VF rows before clinical join", n = nrow(vf_pa))
  },
  vf_ready %>%
    mutate(Infection_Status = ifelse(is.na(Infection_Status), "Missing clinical status", Infection_Status)) %>%
    count(stage = "Canonical VF-ready table", Infection_Status, name = "n"),
  vf_ready %>%
    filter(Infection_Status %in% c("ASB", "UTI")) %>%
    count(stage = "ASB/UTI VF subset", Infection_Status, name = "n")
) %>%
  mutate(stage = factor(stage, levels = c("Clinical status map", "Raw VF P/A matrix",
                                          "Canonical VF-ready table", "ASB/UTI VF subset")))

if (nrow(flow_rows) > 0) {
  p_flow <- ggplot(flow_rows, aes(x = stage, y = n, fill = Infection_Status)) +
    geom_col(width = 0.68) +
    geom_text(aes(label = n), position = position_stack(vjust = 0.5),
              size = 3, colour = "white", check_overlap = TRUE) +
    labs(
      title = "Clinical-to-genomic denominator flow for VF analysis",
      subtitle = "Attrition from clinical episodes to VF-ready and ASB/UTI analysis denominators is shown explicitly",
      x = NULL,
      y = "Rows / episodes",
      fill = "Status or layer",
      caption = paste(
        sprintf("Data: %s, %s, and %s.", FILE_STATUS_MAP, FILE_VF_PA, ready_file),
        "Level of analysis: denominator-flow diagnostic.",
        "Do not hide denominator attrition when interpreting ASB-UTI VF analyses.",
        "Clinical UTI rows may require Uricult-to-UTI-N harmonisation before becoming VF-ready rows."
      )
    ) +
    plot_theme_vf(base_size = 11) +
    theme(axis.text.x = element_text(angle = 20, hjust = 1))

  ggsave(file.path(DIR_PLOTS_VF, "vf_denominator_flow.png"),
         p_flow, width = 8.5, height = 5.8, dpi = 300)
}

if (!is.null(status_map) && "uricult_bridge_applied" %in% names(vf_ready)) {
  bridge_diag <- tibble(
    metric = c("Clinical UTI episodes", "VF-ready UTI rows", "VF-ready UTI rows via Uricult bridge"),
    n = c(
      sum(status_map$Infection_Status == "UTI", na.rm = TRUE),
      sum(vf_ready$Infection_Status == "UTI", na.rm = TRUE),
      sum(vf_ready$Infection_Status == "UTI" & vf_ready$uricult_bridge_applied %in% TRUE, na.rm = TRUE)
    )
  )

  p_bridge <- ggplot(bridge_diag, aes(x = metric, y = n)) +
    geom_col(fill = "#0072B2", width = 0.6) +
    geom_text(aes(label = n), vjust = -0.25, size = 4) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
    labs(
      title = "Uricult clinical events require harmonisation with UTI-labelled WGS isolates",
      subtitle = "Diagnostic count of clinical UTI episodes, VF-ready UTI rows, and Uricult-bridge-linked UTI rows",
      x = NULL,
      y = "Episodes / rows",
      caption = paste(
        sprintf("Data: %s and %s.", FILE_STATUS_MAP, ready_file),
        "Level of analysis: join/denominator diagnostic.",
        "Uricult clinical rows and UTI-N WGS rows are not identical labels; bridge assumptions must be cited.",
        "This plot is not a biological association test."
      )
    ) +
    plot_theme_vf(base_size = 11) +
    theme(axis.text.x = element_text(angle = 20, hjust = 1))

  ggsave(file.path(DIR_PLOTS_VF, "vf_uricult_join_diagnostic.png"),
         p_bridge, width = 7.5, height = 5.2, dpi = 300)
} else {
  msg("Skipping Uricult join diagnostic: status map or uricult_bridge_applied column is absent.")
}

# ==============================================================================
# 6. WRITE SUMMARY
# ==============================================================================

writeLines(summary_lines, file.path(DIR_VF, "vf_lineage_confounding_summary.txt"))
msg("Summary written to %s", file.path(DIR_VF, "vf_lineage_confounding_summary.txt"))

# Also print to console
cat(paste(summary_lines, collapse = "\n"), "\n")

msg("✓ 25_vf_lineage_vf_interaction.R complete.")

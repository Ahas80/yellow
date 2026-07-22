#!/usr/bin/env Rscript

# Compare the two CGE VirulenceFinder threshold profiles with the completed
# primary ABRicate/VFDB analysis. This older diagnostic uses 52 exact/case-
# normalized symbol overlaps. Nine of those are now known to have ambiguous
# reference scope, so the tiered comparison is the scientifically preferred
# figure. Raw database-wide burdens remain deliberately incomparable.

options(stringsAsFactors = FALSE, warn = 1)

required <- c("dplyr", "tidyr", "readr", "ggplot2", "ggrepel", "patchwork", "scales")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Missing R package(s): ", paste(missing, collapse = ", "))

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(ggplot2)
  library(ggrepel)
  library(patchwork)
})

root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
batch <- file.path(root, "results", "virulencefinder_cge_3_2_1")
primary <- file.path(root, "results", "research_questions", "_inputs", "vf_presence_absence_532.csv")
out <- file.path(batch, "concordance", "three_method_comparison")
dir.create(out, recursive = TRUE, showWarnings = FALSE)

atomic_write_csv <- function(x, path) {
  tmp <- tempfile(pattern = paste0(basename(path), "."), tmpdir = dirname(path), fileext = ".tmp")
  readr::write_csv(x, tmp, na = "")
  if (file.exists(path) && !file.remove(path)) stop("Could not replace ", path)
  if (!file.rename(tmp, path)) stop("Could not publish ", path)
}

atomic_write_lines <- function(x, path) {
  tmp <- tempfile(pattern = paste0(basename(path), "."), tmpdir = dirname(path), fileext = ".tmp")
  writeLines(x, tmp, useBytes = TRUE)
  if (file.exists(path) && !file.remove(path)) stop("Could not replace ", path)
  if (!file.rename(tmp, path)) stop("Could not publish ", path)
}

web <- read_csv(file.path(batch, "presence_absence_web_default_id90_cov60.csv"),
                show_col_types = FALSE, progress = FALSE)
matched <- read_csv(file.path(batch, "presence_absence_matched_id80_cov80.csv"),
                    show_col_types = FALSE, progress = FALSE)
vfdb <- read_csv(primary, show_col_types = FALSE, progress = FALSE)
crosswalk <- read_csv(file.path(batch, "concordance", "cge_to_primary_vfdb_gene_crosswalk.csv"),
                      show_col_types = FALSE, progress = FALSE)

metadata_cge <- c("Participant_id", "tp_lab", "episode_key", "Assembly_ID", "fasta_sha256",
                  "UTI_Status", "Event_type", "full_path")
metadata_vfdb <- c("Participant_id", "tp_lab", "episode_key", "Assembly_ID", "fasta_sha256")
cge_genes <- setdiff(names(web), metadata_cge)
vfdb_genes <- setdiff(names(vfdb), metadata_vfdb)
shared <- crosswalk %>%
  filter(tolower(as.character(eligible_for_concordance)) %in% c("true", "t", "1")) %>%
  arrange(cge_gene_family)

stopifnot(
  nrow(web) == 532L, nrow(matched) == 532L, nrow(vfdb) == 532L,
  length(cge_genes) == 681L, length(vfdb_genes) == 227L, nrow(shared) == 52L,
  setequal(web$episode_key, matched$episode_key), setequal(web$episode_key, vfdb$episode_key),
  !anyDuplicated(shared$cge_gene_family), !anyDuplicated(shared$primary_vfdb_gene)
)

# Align every matrix to the primary episode-key order before comparing cells.
web <- web[match(vfdb$episode_key, web$episode_key), ]
matched <- matched[match(vfdb$episode_key, matched$episode_key), ]

as_binary_matrix <- function(data, columns) {
  value <- as.matrix(data[, columns, drop = FALSE])
  storage.mode(value) <- "integer"
  value
}

web_all <- as_binary_matrix(web, cge_genes)
matched_all <- as_binary_matrix(matched, cge_genes)
web_shared <- as_binary_matrix(web, shared$cge_gene_family)
matched_shared <- as_binary_matrix(matched, shared$cge_gene_family)
vfdb_shared <- as_binary_matrix(vfdb, shared$primary_vfdb_gene)

pair_metrics <- function(a, b) {
  intersection <- rowSums(a == 1L & b == 1L)
  union <- rowSums(a == 1L | b == 1L)
  jaccard <- ifelse(union > 0L, intersection / union, 1)
  list(
    median_jaccard = median(jaccard), mean_jaccard = mean(jaccard),
    binary_cell_agreement = mean(a == b), exact_profile_n = sum(rowSums(a != b) == 0L),
    exact_profile_pct = mean(rowSums(a != b) == 0L),
    burden_spearman = suppressWarnings(cor(rowSums(a), rowSums(b), method = "spearman")),
    median_burden_a = median(rowSums(a)), median_burden_b = median(rowSums(b))
  )
}

make_summary <- function(label, universe, a, b) {
  metric <- pair_metrics(a, b)
  tibble(
    comparison = label, comparison_universe = universe,
    n_assemblies = nrow(a), n_gene_families = ncol(a),
    median_jaccard = metric$median_jaccard, mean_jaccard = metric$mean_jaccard,
    binary_cell_agreement = metric$binary_cell_agreement,
    exact_profile_n = metric$exact_profile_n, exact_profile_pct = metric$exact_profile_pct,
    burden_spearman = metric$burden_spearman,
    median_burden_first = metric$median_burden_a,
    median_burden_second = metric$median_burden_b
  )
}

comparison_summary <- bind_rows(
  make_summary("Web 90/60 vs primary VFDB", "52 name-overlap candidates; 9 biologically ambiguous",
               web_shared, vfdb_shared),
  make_summary("Matched 80/80 vs primary VFDB", "52 name-overlap candidates; 9 biologically ambiguous",
               matched_shared, vfdb_shared),
  make_summary("Web 90/60 vs matched 80/80", "52 name-overlap candidates; 9 biologically ambiguous",
               web_shared, matched_shared)
)
atomic_write_csv(comparison_summary, file.path(out, "comparison_summary.csv"))

threshold_full_summary <- make_summary(
  "Web 90/60 vs matched 80/80", "681-family pinned CGE universe", web_all, matched_all
)
atomic_write_csv(threshold_full_summary, file.path(out, "threshold_full_universe_summary.csv"))

shared_prevalence <- shared %>%
  transmute(
    gene_family = cge_gene_family, primary_vfdb_gene,
    mapping_type,
    primary_vfdb_prevalence = colMeans(vfdb_shared),
    web_90_60_prevalence = colMeans(web_shared),
    matched_80_80_prevalence = colMeans(matched_shared)
  ) %>%
  mutate(
    web_minus_primary = web_90_60_prevalence - primary_vfdb_prevalence,
    matched_minus_primary = matched_80_80_prevalence - primary_vfdb_prevalence,
    matched_minus_web = matched_80_80_prevalence - web_90_60_prevalence,
    maximum_absolute_primary_difference = pmax(abs(web_minus_primary), abs(matched_minus_primary))
  ) %>%
  arrange(desc(maximum_absolute_primary_difference), gene_family)
atomic_write_csv(shared_prevalence, file.path(out, "shared_gene_prevalence_three_methods.csv"))

threshold_prevalence <- tibble(
  gene_family = cge_genes,
  web_90_60_prevalence = colMeans(web_all),
  matched_80_80_prevalence = colMeans(matched_all)
) %>%
  mutate(
    matched_minus_web = matched_80_80_prevalence - web_90_60_prevalence,
    absolute_difference = abs(matched_minus_web),
    direction = case_when(
      matched_minus_web > 0 ~ "Higher at matched 80/80",
      matched_minus_web < 0 ~ "Higher at web 90/60",
      TRUE ~ "No difference"
    )
  ) %>%
  arrange(desc(absolute_difference), gene_family)
atomic_write_csv(threshold_prevalence, file.path(out, "threshold_gene_prevalence_differences.csv"))

# Explicitly labelled post-hoc scope check for the largest formal cross-database
# discrepancy. This does not replace the prespecified exact-symbol comparison.
kpsm_family_columns <- grep("^kpsM", cge_genes, value = TRUE)
primary_kpsm <- as.integer(vfdb$kpsM > 0)
capsule_scope_sensitivity <- bind_rows(lapply(
  list(`CGE web 90/60` = web, `CGE matched 80/80` = matched),
  function(data) {
    exact <- as.integer(data$kpsM > 0)
    prefix <- as.integer(rowSums(as_binary_matrix(data, kpsm_family_columns)) > 0)
    bind_rows(
      tibble(scope = "prespecified exact kpsM", cge_present = exact),
      tibble(scope = "exploratory any kpsM-prefix family", cge_present = prefix)
    ) %>%
      group_by(scope) %>%
      summarise(
        cge_present_n = sum(cge_present), primary_kpsM_present_n = sum(primary_kpsm),
        both_present_n = sum(cge_present == 1L & primary_kpsm == 1L),
        cge_only_n = sum(cge_present == 1L & primary_kpsm == 0L),
        primary_only_n = sum(cge_present == 0L & primary_kpsm == 1L),
        both_absent_n = sum(cge_present == 0L & primary_kpsm == 0L),
        agreement = mean(cge_present == primary_kpsm), .groups = "drop"
      )
  }
), .id = "cge_profile")
atomic_write_csv(capsule_scope_sensitivity, file.path(out, "capsule_kpsM_scope_sensitivity.csv"))

universe <- tibble(
  definition = factor(
    c("CGE pinned reference families", "Primary VFDB retained features", "Name-overlap candidates"),
    levels = rev(c("CGE pinned reference families", "Primary VFDB retained features", "Name-overlap candidates"))
  ),
  n = c(length(cge_genes), length(vfdb_genes), nrow(shared)),
  detail = c(
    paste0(sum(colSums(web_all | matched_all) > 0), " observed"),
    paste0(sum(colSums(as_binary_matrix(vfdb, vfdb_genes)) > 0), " observed"),
    "cross-database overlap"
  )
)

palette <- c(
  "Primary ABRicate/VFDB" = "#4D4D4D",
  "CGE web 90/60" = "#0072B2",
  "CGE matched 80/80" = "#D55E00",
  "Higher at matched 80/80" = "#D55E00",
  "Higher at web 90/60" = "#0072B2"
)

p_universe <- ggplot(universe, aes(n, definition)) +
  geom_col(width = 0.62, fill = "#6B7280") +
  geom_text(aes(label = paste0(n, "  (", detail, ")")), hjust = -0.02, size = 3.2) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.34))) +
  labs(title = "Feature definitions are not interchangeable", x = "Number of feature columns", y = NULL,
       caption = "The 52 symbol overlaps are candidates, not 52 proven one-to-one biological mappings.") +
  theme_minimal(base_size = 10) +
  theme(panel.grid.major.y = element_blank(), panel.grid.minor = element_blank())

scatter <- shared_prevalence %>%
  select(gene_family, primary_vfdb_gene, primary_vfdb_prevalence,
         `CGE web 90/60` = web_90_60_prevalence,
         `CGE matched 80/80` = matched_80_80_prevalence,
         maximum_absolute_primary_difference) %>%
  pivot_longer(starts_with("CGE"), names_to = "method", values_to = "cge_prevalence")
label_genes <- shared_prevalence %>% slice_max(maximum_absolute_primary_difference, n = 6, with_ties = FALSE) %>%
  pull(gene_family)

p_scatter <- ggplot(scatter, aes(primary_vfdb_prevalence, cge_prevalence, colour = method, shape = method)) +
  geom_abline(slope = 1, intercept = 0, colour = "grey65", linewidth = 0.5) +
  geom_point(alpha = 0.72, size = 2) +
  ggrepel::geom_text_repel(
    data = filter(scatter, gene_family %in% label_genes, method == "CGE matched 80/80"),
    aes(label = ifelse(gene_family %in% c("kpsM", "sfaE", "focC"),
                       paste0(gene_family, "*"), gene_family)),
    colour = "grey20", size = 2.8, max.overlaps = Inf,
    box.padding = 0.25, point.padding = 0.15, min.segment.length = 0
  ) +
  scale_colour_manual(values = palette) +
  scale_shape_manual(values = c("CGE web 90/60" = 16, "CGE matched 80/80" = 17)) +
  scale_x_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0, 1)) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0, 1)) +
  coord_equal() +
  labs(title = "Prevalence across 52 name-overlap candidates",
       x = "Primary ABRicate/VFDB prevalence", y = "CGE VirulenceFinder prevalence",
       colour = NULL, shape = NULL,
       caption = "Asterisks flag method-specific reference scope; differences are not biological gains/losses.") +
  theme_minimal(base_size = 10) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank())

top_threshold <- threshold_prevalence %>%
  filter(absolute_difference > 0) %>%
  slice_max(absolute_difference, n = 12, with_ties = FALSE) %>%
  arrange(matched_minus_web) %>%
  mutate(gene_family = factor(gene_family, levels = gene_family))

p_threshold <- ggplot(top_threshold, aes(matched_minus_web, gene_family, fill = direction)) +
  geom_col(width = 0.66) +
  geom_vline(xintercept = 0, colour = "grey55", linewidth = 0.45) +
  geom_text(aes(label = scales::percent(matched_minus_web, accuracy = 0.1)),
            hjust = ifelse(top_threshold$matched_minus_web >= 0, -0.08, 1.08), size = 3) +
  scale_fill_manual(values = palette) +
  scale_x_continuous(labels = scales::percent_format(accuracy = 1),
                     expand = expansion(mult = c(0.22, 0.2))) +
  labs(title = "Largest threshold-sensitive CGE families",
       subtitle = "Matched 80/80 prevalence minus web 90/60 prevalence",
       x = "Prevalence-point difference", y = NULL, fill = NULL,
       caption = "Positive and negative values confirm that the profiles are complementary, not nested.") +
  theme_minimal(base_size = 10) +
  theme(legend.position = "bottom", panel.grid.major.y = element_blank(), panel.grid.minor = element_blank())

agreement_plot <- comparison_summary %>%
  select(comparison, comparison_universe, median_jaccard, exact_profile_pct,
         exact_profile_n) %>%
  pivot_longer(c(median_jaccard, exact_profile_pct), names_to = "metric", values_to = "value") %>%
  mutate(
    metric = recode(metric, median_jaccard = "Median isolate Jaccard",
                    exact_profile_pct = "Exact shared profile"),
    comparison = factor(comparison, levels = rev(comparison_summary$comparison))
  )

p_agreement <- ggplot(agreement_plot, aes(value, comparison, fill = metric)) +
  geom_col(position = position_dodge(width = 0.68), width = 0.6) +
  geom_text(aes(label = scales::percent(value, accuracy = 0.1)),
            position = position_dodge(width = 0.68), hjust = -0.08, size = 3) +
  scale_fill_manual(values = c("Median isolate Jaccard" = "#009E73",
                               "Exact shared profile" = "#CC79A7")) +
  scale_x_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0, 1.08)) +
  labs(title = "Agreement across the same 532 assemblies", x = "Agreement", y = NULL, fill = NULL,
       caption = paste0(
         "Exact whole-profile matches: web-VFDB ", comparison_summary$exact_profile_n[1], "/532; ",
         "matched-VFDB ", comparison_summary$exact_profile_n[2], "/532; ",
         "web-matched ", comparison_summary$exact_profile_n[3], "/532."
       )) +
  theme_minimal(base_size = 10) +
  theme(legend.position = "bottom", panel.grid.major.y = element_blank(), panel.grid.minor = element_blank())

figure <- (p_universe | p_agreement) / (p_scatter | p_threshold) +
  plot_annotation(
    title = "Historical name-only diagnostic: CGE versus primary VFDB",
    subtitle = "532 assemblies; 52 symbol overlaps, including nine biologically ambiguous mappings",
    caption = "Use the tiered comparison for inference. This name-only view is retained for transparency; raw database-wide burdens are not compared.",
    theme = theme(plot.title = element_text(face = "bold", size = 15),
                  plot.subtitle = element_text(size = 10), plot.caption = element_text(size = 9))
  )

ggsave(file.path(out, "virulence_factor_comparison.png"), figure,
       width = 13, height = 10, units = "in", dpi = 300, bg = "white")
ggsave(file.path(out, "virulence_factor_comparison.pdf"), figure,
       width = 13, height = 10, units = "in", bg = "white")

atomic_write_lines(c(
  "# Historical name-only three-method comparison",
  "",
  "This comparison uses all 532 selected assemblies.",
  "",
  "- Threshold comparison: both CGE profiles are compared over the full pinned 681-family CGE universe.",
  "- This historical diagnostic uses 52 exact/case-normalized symbol overlaps; nine are now classified as biologically ambiguous.",
  "- Use `../tiered_cross_database_comparison/tiered_virulence_comparison.*` as the preferred cross-database figure.",
  "- Raw total burden across CGE and VFDB is intentionally not compared.",
  "- `kpsM` is retained as an exact symbol match but requires biological caution because CGE also separates multiple capsule-family variants such as `kpsMII`.",
  "- Differences indicate method/database sensitivity and must not automatically be interpreted as biological acquisition or loss.",
  "",
  "Generated files:",
  "",
  "- `comparison_summary.csv`",
  "- `threshold_full_universe_summary.csv`",
  "- `shared_gene_prevalence_three_methods.csv`",
  "- `threshold_gene_prevalence_differences.csv`",
  "- `capsule_kpsM_scope_sensitivity.csv` (explicitly exploratory)",
  "- `virulence_factor_comparison.png`",
  "- `virulence_factor_comparison.pdf`"
), file.path(out, "README.md"))

message("Three-method virulence-factor comparison: PASS")

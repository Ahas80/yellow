#!/usr/bin/env Rscript

# Publication-style visualization of the conservative CGE/VFDB crosswalk.
#
# USER-TUNABLE DISPLAY SETTINGS
# -----------------------------
TOP_DIRECT_TARGETS <- 12L # Number of direct pairs shown in the difference panel.
FIGURE_WIDTH_IN <- 13     # Output width; does not change analytical results.
FIGURE_HEIGHT_IN <- 10    # Output height; does not change analytical results.
FIGURE_DPI <- 300L        # Raster resolution for the PNG.

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
out <- file.path(
  root, "results", "virulencefinder_cge_3_2_1", "concordance",
  "tiered_cross_database_comparison"
)

atomic_write_csv <- function(x, path) {
  tmp <- tempfile(pattern = paste0(basename(path), "."), tmpdir = dirname(path), fileext = ".tmp")
  readr::write_csv(x, tmp, na = "")
  if (file.exists(path) && !file.remove(path)) stop("Could not replace ", path)
  if (!file.rename(tmp, path)) stop("Could not publish ", path)
}

scope <- read_csv(file.path(out, "mapping_scope_summary.csv"), show_col_types = FALSE)
direct <- read_csv(file.path(out, "strict_direct_gene_concordance.csv"), show_col_types = FALSE)
profile <- read_csv(file.path(out, "strict_direct_profile_summary.csv"), show_col_types = FALSE)
systems <- read_csv(file.path(out, "system_level_concordance.csv"), show_col_types = FALSE)
ambiguous <- read_csv(file.path(out, "ambiguous_exact_name_diagnostics.csv"), show_col_types = FALSE)

stopifnot(
  nrow(scope) == 8L,
  nrow(direct) == 96L,
  length(unique(paste(direct$primary_feature, direct$cge_family))) == 48L,
  nrow(profile) == 2L,
  nrow(systems) == 14L,
  length(unique(systems$system_id)) == 7L,
  nrow(ambiguous) == 18L
)

profile_labels <- c(
  web_default_id90_cov60 = "CGE web 90/60",
  matched_id80_cov80 = "CGE matched 80/80"
)
method_colours <- c(
  "Primary ABRicate/VFDB" = "#4D4D4D",
  "CGE web 90/60" = "#0072B2",
  "CGE matched 80/80" = "#D55E00"
)
scope_colours <- c(
  strict_direct = "#009E73",
  ambiguous_exact_name = "#CC79A7",
  aggregate_only = "#E69F00",
  not_comparable = "#9AA0A6"
)

scope_plot <- scope %>%
  mutate(
    method = recode(
      method,
      primary_ABRicate_VFDB = "ABRicate/VFDB\n227 features",
      CGE_VirulenceFinder = "CGE VirulenceFinder\n681 families"
    ),
    mapping_category = factor(
      mapping_category,
      levels = c("strict_direct", "ambiguous_exact_name", "aggregate_only", "not_comparable"),
      labels = c("Strict direct", "Ambiguous name", "Aggregate only", "Not comparable")
    )
  )

p_scope <- ggplot(scope_plot, aes(proportion, method, fill = mapping_category)) +
  geom_col(width = 0.62, colour = "white", linewidth = 0.25) +
  geom_text(
    aes(label = ifelse(proportion >= 0.06, paste0(n, "\n", scales::percent(proportion, accuracy = 0.1)), "")),
    position = position_stack(vjust = 0.5), size = 3, colour = "white", lineheight = 0.9
  ) +
  scale_fill_manual(values = setNames(scope_colours, levels(scope_plot$mapping_category))) +
  scale_x_continuous(labels = scales::percent_format(accuracy = 1), expand = c(0, 0)) +
  labs(
    title = "Only a minority supports direct gene-wise comparison",
    x = "Share of each method-specific target universe", y = NULL, fill = NULL,
    caption = "Unmapped means not comparable; it does not mean biologically absent."
  ) +
  theme_minimal(base_size = 10) +
  theme(
    legend.position = "bottom", panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank()
  )

direct_rank <- direct %>%
  group_by(primary_feature, cge_family) %>%
  summarise(max_abs_difference = max(abs(cge_minus_primary_prevalence)), .groups = "drop") %>%
  slice_max(max_abs_difference, n = TOP_DIRECT_TARGETS, with_ties = FALSE)

direct_top <- direct %>%
  semi_join(direct_rank, by = c("primary_feature", "cge_family")) %>%
  mutate(
    profile = recode(profile, !!!profile_labels),
    target = ifelse(primary_feature == cge_family, primary_feature,
                    paste0(primary_feature, " <-> ", cge_family))
  ) %>%
  arrange(cge_minus_primary_prevalence) %>%
  mutate(target = factor(target, levels = unique(target)))

p_direct <- ggplot(
  direct_top,
  aes(cge_minus_primary_prevalence, target, colour = profile, shape = profile)
) +
  geom_vline(xintercept = 0, colour = "grey65", linewidth = 0.45) +
  geom_segment(aes(x = 0, xend = cge_minus_primary_prevalence, yend = target),
               colour = "grey78", linewidth = 0.45) +
  geom_point(size = 2.4, position = position_dodge(width = 0.45)) +
  scale_colour_manual(values = method_colours[names(method_colours) != "Primary ABRicate/VFDB"]) +
  scale_shape_manual(values = c("CGE web 90/60" = 16, "CGE matched 80/80" = 17)) +
  scale_x_continuous(labels = scales::percent_format(accuracy = 1),
                     expand = expansion(mult = c(0.18, 0.18))) +
  labs(
    title = "Largest differences among 48 strict direct pairs",
    subtitle = "CGE prevalence minus primary VFDB prevalence",
    x = "Prevalence-point difference", y = NULL, colour = NULL, shape = NULL,
    caption = "Direct pairs are sequence/annotation-supported and provisionally curated."
  ) +
  theme_minimal(base_size = 10) +
  theme(
    legend.position = "bottom", panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank()
  )

system_prevalence <- systems %>%
  select(profile, system_id, system_label, primary_prevalence, cge_prevalence) %>%
  mutate(profile = recode(profile, !!!profile_labels)) %>%
  pivot_longer(c(primary_prevalence, cge_prevalence), names_to = "source", values_to = "prevalence") %>%
  mutate(
    method = ifelse(source == "primary_prevalence", "Primary ABRicate/VFDB", profile),
    system_label = recode(
      system_label,
      `Intimin (any CGE eae subtype)` = "Intimin / eae",
      `PapA major pilin (any CGE papA allele)` = "PapA major pilin",
      `Capsule KpsM (any CGE kpsM family)` = "Capsule / kpsM",
      `NleB effector (either primary variant)` = "NleB effector",
      `ChuA/ShuA outer-membrane heme receptor` = "ChuA/ShuA receptor",
      `Afa/Dr/DAF adherence system` = "Afa/Dr/DAF system",
      `Sfa/Foc fimbrial system` = "Sfa/Foc system"
    )
  ) %>%
  distinct(system_id, system_label, method, prevalence) %>%
  mutate(
    system_label = factor(
      system_label,
      levels = rev(c("Intimin / eae", "NleB effector", "ChuA/ShuA receptor",
                     "Capsule / kpsM", "Afa/Dr/DAF system", "Sfa/Foc system",
                     "PapA major pilin"))
    )
  )

p_system <- ggplot(system_prevalence, aes(prevalence, system_label, colour = method, shape = method)) +
  geom_line(aes(group = system_id), colour = "grey78", linewidth = 0.5) +
  geom_point(size = 2.5) +
  scale_colour_manual(values = method_colours) +
  scale_shape_manual(values = c(
    "Primary ABRicate/VFDB" = 15, "CGE web 90/60" = 16, "CGE matched 80/80" = 17
  )) +
  scale_x_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0, 1)) +
  labs(
    title = "Seven aggregate endpoints compare system presence only",
    x = "Assembly prevalence", y = NULL, colour = NULL, shape = NULL,
    caption = "Member counts and raw burden are deliberately not compared."
  ) +
  theme_minimal(base_size = 10) +
  theme(
    legend.position = "bottom", panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank()
  )

agreement <- profile %>%
  transmute(
    profile = recode(profile, !!!profile_labels),
    `Mean episode Jaccard` = mean_episode_jaccard,
    `Exact 48-target profile` = exact_profile_proportion
  ) %>%
  pivot_longer(-profile, names_to = "metric", values_to = "value")

p_agreement <- ggplot(agreement, aes(value, profile, fill = metric)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.62) +
  geom_text(
    aes(label = scales::percent(value, accuracy = 0.1)),
    position = position_dodge(width = 0.7), hjust = -0.08, size = 3
  ) +
  scale_fill_manual(values = c("Mean episode Jaccard" = "#009E73", "Exact 48-target profile" = "#56B4E9")) +
  scale_x_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0, 1.08)) +
  labs(
    title = "Strict direct agreement is high",
    x = "Agreement across 532 assemblies", y = NULL, fill = NULL,
    caption = "Nine ambiguous same-name pairs are excluded from these metrics."
  ) +
  theme_minimal(base_size = 10) +
  theme(
    legend.position = "bottom", panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank()
  )

figure <- (p_scope | p_agreement) / (p_direct | p_system) +
  plot_annotation(
    title = "Tiered comparison of primary VFDB features and CGE gene families",
    subtitle = "532 QC-passed canonical Longcycler assemblies; direct, aggregate and non-comparable targets remain separate",
    caption = paste0(
      "Exploratory harmonization only. Primary ABRicate/VFDB remains authoritative; ",
      "differences are method/database detections, not biological gain or loss."
    ),
    theme = theme(
      plot.title = element_text(face = "bold", size = 15),
      plot.subtitle = element_text(size = 10), plot.caption = element_text(size = 9)
    )
  )

ggsave(file.path(out, "tiered_virulence_comparison.png"), figure,
       width = FIGURE_WIDTH_IN, height = FIGURE_HEIGHT_IN, units = "in",
       dpi = FIGURE_DPI, bg = "white")
ggsave(file.path(out, "tiered_virulence_comparison.pdf"), figure,
       width = FIGURE_WIDTH_IN, height = FIGURE_HEIGHT_IN, units = "in", bg = "white")

key_findings <- bind_rows(
  profile %>% transmute(
    comparison_unit = "strict_direct_profile", item = profile,
    primary_prevalence = NA_real_, cge_prevalence = NA_real_,
    mean_episode_jaccard, exact_profile_proportion,
    interpretation = "48 strict direct pairs only"
  ),
  systems %>% filter(profile == "matched_id80_cov80") %>% transmute(
    comparison_unit = "system_presence", item = system_label,
    primary_prevalence, cge_prevalence,
    mean_episode_jaccard = NA_real_, exact_profile_proportion = overall_agreement,
    interpretation = "Binary any-member endpoint; no component-count comparison"
  )
)
atomic_write_csv(key_findings, file.path(out, "comparison_key_findings.csv"))

message("Tiered virulence comparison visualization: PASS")

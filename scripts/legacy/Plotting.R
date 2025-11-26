###############################################################################
# 0 ·  HOUSE-KEEPING  ----------------------------------------------------------
###############################################################################
Sys.setenv(PATH = paste(Sys.getenv("PATH"), "/opt/homebrew/bin", sep = ":"))

pkgs <- c("dplyr", "tidyr", "purrr", "readr", "ggplot2",
          "ComplexUpset", "tibble", "stringr")
invisible(lapply(pkgs, require, character.only = TRUE))

theme_set(theme_minimal(base_size = 10))  # consistent theme across plots

plots_dir <- file.path(getwd(), "plots")
if (!dir.exists(plots_dir)) dir.create(plots_dir)



###############################################################################
# 3a ·  DESCRIPTIVE-STATISTICS  ───────────────────────────────────────────────
#       ▸ One-liner tables saved as CSV
#       ▸ Human-readable summary saved as   stats_descriptive.md
###############################################################################

# --- sample-level stats -------------------------------------------------------
genes_per_sample <- vf_hits %>%
  count(Participant_id, Timepoint, name = "n_genes")        # genes / sample

samples_tbl <- genes_per_sample %>%
  arrange(Participant_id, Timepoint)

# --- participant-level stats --------------------------------------------------
#uricult is ignored
participant_tbl <- samples_tbl %>%
  mutate(tp_num = parse_number(Timepoint)) %>%   # "T0" → 0, "Uricult" → NA
  filter(!is.na(tp_num)) %>%                     # drop the Uricult rows
  arrange(Participant_id, tp_num) %>%            # make sure rows are ordered
  group_by(Participant_id) %>%
  summarise(
    n_timepoints = n(),                          # 3, 4, 5 …
    mean_genes   = mean(n_genes),
    sd_genes     = sd(n_genes),
    min_genes    = min(n_genes),
    max_genes    = max(n_genes),
    first_tp     = dplyr::first(tp_num),         # numeric 0,1,2…
    last_tp      = dplyr::last(tp_num),
    genes_T0     = dplyr::first(n_genes),        # n_genes at earliest TP
    genes_last   = dplyr::last(n_genes),         # n_genes at latest  TP
    delta        = genes_last - genes_T0,
    .groups      = "drop"
  )




# --- cohort-level stats -------------------------------------------------------
cohort_tbl <- participant_tbl %>%
  summarise(
    participants       = n(),
    timepoints_total   = sum(n_timepoints),
    genes_median_T0    = median(genes_T0),
    genes_median_last  = median(genes_last),
    mean_delta         = mean(delta),
    sd_delta           = sd(delta)
  )

# --- write CSVs ---------------------------------------------------------------
write_csv(samples_tbl,      "stats_sample_level.csv")
write_csv(participant_tbl,  "stats_participant_level.csv")
write_csv(cohort_tbl,       "stats_cohort_level.csv")
message("✓ wrote sample / participant / cohort stats")

# --- Markdown summary ---------------------------------------------------------
md <- c(
  "# Descriptive statistics  \n",
  "Generated on: ", Sys.Date(), "  \n\n",
  
  "## Cohort-wide overview\n",
  knitr::kable(cohort_tbl, format = "markdown"), "\n\n",
  
  "## Participant-level summary\n",
  knitr::kable(participant_tbl, format = "markdown"), "\n\n",
  
  "## Per-sample detail (first 10 rows)\n",
  knitr::kable(head(samples_tbl, 10), format = "markdown"), "\n"
)

writeLines(md, "stats_descriptive.md")
message("✓ wrote stats_descriptive.md  (easy to copy-paste into a report)")


###############################################################################
# 3a ·  QUICK DESCRIPTIVE STATISTICS
###############################################################################
genes_per_sample <- vf_hits %>%
  count(Participant_id, Timepoint, name = "n_genes")

tbl_core <- genes_per_sample %>%
  group_by(Participant_id) %>%
  summarise(
    N_samples  = n(),
    Genes_T0   = n_genes[Timepoint == min(Timepoint)],
    Genes_last = n_genes[Timepoint == max(Timepoint)],
    .groups    = "drop"
  )

write_csv(tbl_core, "summary_by_participant.csv")
message("✓ wrote summary_by_participant.csv")


###############################################################################
# 4 ·  GAIN / LOSS TABLE
###############################################################################
gain_loss_tbl <- vf_hits %>%
  select(Participant_id, assembler, Timepoint, GENE) %>%
  group_by(Participant_id, assembler, Timepoint) %>%
  summarise(genes = list(unique(GENE)), .groups = "drop") %>%
  arrange(Participant_id, assembler, Timepoint) %>%
  group_by(Participant_id, assembler) %>%
  mutate(prev_genes = lag(genes),
         from_tp    = lag(Timepoint),
         gained     = map2(genes, prev_genes, ~ setdiff(.x, .y)),
         lost       = map2(prev_genes, genes, ~ setdiff(.x, .y))) %>%
  ungroup() %>%
  filter(!map_lgl(prev_genes, is.null)) %>%
  select(Participant_id, assembler, from_tp, to_tp = Timepoint, gained, lost)

write_csv(gain_loss_tbl, "vfdb_gain_loss_by_participant.csv")
message("✓ wrote vfdb_gain_loss_by_participant.csv")

vf_pa <- vf_hits %>%
  distinct(Participant_id, Timepoint, GENE) %>%   # de-duplicate hits
  mutate(present = 1) %>%
  pivot_wider(
    names_from  = GENE,
    values_from = present,
    values_fill = 0
  )

# quick sanity check
print(dim(vf_pa))          # rows = samples, cols = 2 + #genes
head(vf_pa[, 1:10])
###############################################################################
# 5 ·  VARIABLE-GENE HEATMAPS
###############################################################################
variable_genes <- vf_hits %>%
  distinct(Participant_id, GENE) %>%
  count(GENE, name = "n_participants") %>%
  filter(between(n_participants, 1, n_distinct(vf_hits$Participant_id) - 1)) %>%
  pull(GENE)

vf_pa_var <- vf_pa %>%
  select(Participant_id, Timepoint, all_of(variable_genes))

plot_heatmap_var <- function(id) {
  df <- vf_pa_var %>%
    filter(Participant_id == id) %>%
    pivot_longer(-c(Participant_id, Timepoint),
                 names_to = "GENE", values_to = "present")
  
  total_genes <- df %>% filter(present == 1) %>% nrow()
  
  ggplot(df, aes(GENE, Timepoint, fill = factor(present))) +
    geom_tile(color = "grey80") +
    scale_fill_manual(values = c(`0` = "white", `1` = "steelblue")) +
    labs(title = paste("Variable VF genes – Participant", id),
         subtitle = paste(total_genes, "total gene-timepoint hits"),
         x = NULL, y = "Timepoint", fill = "Present") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1),
          panel.grid = element_blank())
}


###############################################################################
# 6 ·  COMPLEX UPSET DATA + PLOTTERF
###############################################################################
get_upset_df <- function(id) {
  vf_pa %>%
    filter(Participant_id == id) %>%
    pivot_longer(-c(Participant_id, Timepoint),
                 names_to = "GENE", values_to = "present") %>%
    filter(present == 1) %>%
    unite(GENE_tp, GENE, Timepoint, sep = "_") %>%
    mutate(val = TRUE) %>%
    pivot_wider(names_from = GENE_tp, values_from = val, values_fill = FALSE) %>%
    mutate(across(-Participant_id, as.logical))
}

show_upset <- function(id) {
  df <- get_upset_df(id)
  n_sets <- ncol(df) - 1
  
  ComplexUpset::upset(df,
                      intersect = colnames(df)[-1],
                      name = paste("Participant", id),
                      base_annotations = list(
                        'Intersection size' = intersection_size(
                          text = list(size = 3),
                          bar_number_threshold = 0
                        )
                      ),
                      themes = upset_modify_themes(
                        list('intersections_matrix'=theme(axis.text.y=element_text(size=7)))
                      )) +
    labs(title = paste("ComplexUpSet – Participant", id),
         subtitle = paste(n_sets, "gene-timepoint combinations"))
}


###############################################################################
# 7 ·  SAVE TO PLOTS/
###############################################################################
for (id in sort(unique(vf_pa$Participant_id))) {
  message("Rendering heat-map / ComplexUpSet for ", id)
  
  ggsave(file.path(plots_dir, sprintf("heatmap_%s.png", id)),
         plot = plot_heatmap_var(id),
         width = 12, height = 6, dpi = 300)
  
  ggsave(file.path(plots_dir, sprintf("upset_%s.png", id)),
         plot = show_upset(id),
         width = 12, height = 8, dpi = 300)
}

message("✓  Annotated plots saved to ", plots_dir)

###############################################################################
# 8 ·  PER-PARTICIPANT  ↔  PAIRWISE NUCLEOTIDE-IDENTITY PLOTS
#     – one graphic per Participant × Assembler
#       • X-axis  :           From-time-point (tp_A)
#       • Y-axis  :           % identity to the **next** time-point (tp_B)
#       • Point   :           colour-coded by ASB / UTI suspicion
#       • Facet   :           none (1 plot ⇒ 1 PNG)   ─ keeps slides clean
###############################################################################

# ---- OUTPUT DIRECTORY --------------------------------------------------------
pair_plot_dir <- file.path(plots_dir, "pairwise_identity")
if (!dir.exists(pair_plot_dir)) dir.create(pair_plot_dir, recursive = TRUE)

# ---- PLOTTER -----------------------------------------------------------------
plot_pair_id <- function(df, pid, asm) {
  
  sub <- df |>
    filter(Participant_id == pid,
           assembler      == asm) |>
    arrange(tp_A)
  
  if (nrow(sub) == 0) return(NULL)
  
  ggplot(sub,
         aes(x = tp_A,
             y = AvgIdentity,
             group = 1,
             colour = UTI_flag_A)) +
    geom_line(linewidth = 0.8) +        # ← updated (was size = 0.8)
    geom_point(size = 3) +              # points can still use `size`
    scale_y_continuous(limits = c(min(sub$AvgIdentity, na.rm = TRUE) - 0.2,
                                  max(sub$AvgIdentity, na.rm = TRUE) + 0.2),
                       expand = c(0, 0)) +
    scale_colour_manual(values = c(ASB           = "#3182bd",
                                   UTI_suspected = "#de2d26"),
                        name = "Clinical flag") +
    labs(title  = glue::glue("Pairwise identity across time – P{pid} ({asm})"),
         x = "Time-point A (vs A+1)",
         y = "Avg nucleotide identity (%)") +
    theme_minimal(base_size = 10) +
    theme(panel.grid.minor = element_blank())
}
# ---- LOOP & SAVE -------------------------------------------------------------
for (pid in selected_ids) {
  for (asm in c("flye", "longcycler")) {
    
    g <- plot_pair_id(pairwise_annotated, pid, asm)
    
    if (is.null(g)) next
    
    ggsave(
      filename = file.path(pair_plot_dir,
                           glue::glue("pairwise_id_P{pid}_{asm}.png")),
      plot     = g,
      width    = 6,
      height   = 4,
      dpi      = 300
    )
    message("✓  pairwise identity → ",
            glue::glue("pairwise_id_P{pid}_{asm}.png"))
  }
}

###############################################################################
# END  Section 8
###############################################################################
write_csv(tbl_core, "stats_participant_level.csv")
write_csv(genes_per_sample, "stats_sample_level.csv")
###############################################################################
# Gene-level prevalence table  ->  tbl_gene
###############################################################################
library(dplyr)

tbl_gene <- vf_hits %>%                             # long table of hits
  distinct(Participant_id, Timepoint, GENE) %>%     # one row per sample hit
  group_by(GENE) %>%
  summarise(
    n_participants = n_distinct(Participant_id),    # in how many participants
    total_hits     = n(),                           # total sample-level hits
    .groups = "drop"
  ) %>%
  arrange(desc(n_participants))

write_csv(tbl_gene, "stats_gene_level.csv")
message("✓ wrote stats_gene_level.csv")

# Cohort-wide summaries
stats_cohort <- tibble(
  n_participants   = n_distinct(vf_pa$Participant_id),
  median_samples   = median(tbl_core$N_samples),
  median_genes_t0  = median(tbl_core$Genes_T0,   na.rm = TRUE),
  median_genes_end = median(tbl_core$Genes_last, na.rm = TRUE)
)
write_csv(stats_cohort, "stats_cohort_level.csv")


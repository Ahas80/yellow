#!/usr/bin/env Rscript
# ==============================================================================
# 26_vf_define_gene_modules.R
# ==============================================================================
#
# GOAL:
#   Curate detected VF genes into biologically interpretable modules/systems
#   and create episode-level module presence/count tables for downstream
#   VF scoring, longitudinal transition analysis, and lineage-aware modelling.
#
# METHOD:
#   1. Load canonical VF analysis dataset from script 22.
#   2. Load current gene_map from script 04.
#   3. Apply conservative, reproducible gene-to-module assignment rules.
#   4. Preserve all unassigned and ambiguous genes in a review table.
#   5. Collapse gene-level presence/absence into module-level presence/counts.
#
# OUTPUT:
#   - results/vf/gene_module_map.csv
#   - results/vf/vf_module_presence_by_episode.csv
#   - results/vf/vf_module_summary.csv
#   - results/vf/vf_module_assignment_audit.csv
#   - results/vf/gene_module_unassigned_review.csv
#   - results/vf/vf_module_qc_report.txt
#   - results/vf/vf_module_definition_notes.md
#   - plots/vf/module_gene_counts.png
#   - plots/vf/module_prevalence_by_status.png
#
# NOTE:
#   This is a curation and data-structuring step. It does not test whether
#   modules are associated with UTI. Module definitions should be interpreted
#   as an analysis framework to be validated in downstream scripts.
#   Module presence is based on the existing union-called VF gene matrix.
#   This script does not change gene detection thresholds or assembler logic.
# ==============================================================================

source("00_config.R")
source("R/plot_helpers.R")
suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
})

msg("Starting 26_vf_define_gene_modules.R")

# ==============================================================================
# VF visualisation shared helpers
# ==============================================================================

STATUS_LEVELS <- c("UTI", "Not_UTI", "Unknown")

status_for_plot <- function(x) {
  x <- as.character(x)
  x[is.na(x) | x == ""] <- "Unknown"
  x[!x %in% STATUS_LEVELS] <- "Unknown"
  factor(x, levels = STATUS_LEVELS)
}

plot_theme_vf <- function(base_size = 11) {
  theme_bw(base_size = base_size) +
    theme(
      plot.caption = element_text(hjust = 0, size = base_size - 3, colour = "grey35"),
      plot.subtitle = element_text(colour = "grey25"),
      legend.position = "bottom"
    )
}

format_category_label <- function(x) {
  recode(
    x,
    "Adhesion/Fimbriae" = "Adhesion / fimbriae",
    "Iron acquisition" = "Iron acquisition",
    "Toxins" = "Toxins",
    "Capsule/Surface" = "Capsule / surface structures",
    "Invasion/Evasion" = "Invasion / immune evasion",
    "Not in gene_map" = "Unassigned VFDB hits not in gene map",
    "Unassigned" = "Unassigned VFDB hits",
    .default = as.character(x)
  )
}

# ==============================================================================
# 1. DEFINE INPUT AND OUTPUT PATHS
# ==============================================================================
FILE_GENE_MAP   <- file.path(DIR_VF, "gene_map.csv")
FILE_VF_PA_RAW   <- FILE_VF_PA

OUT_MODULE_MAP     <- file.path(DIR_VF, "gene_module_map.csv")
OUT_MODULE_EP      <- file.path(DIR_VF, "vf_module_presence_by_episode.csv")
OUT_MODULE_SUMMARY <- file.path(DIR_VF, "vf_module_summary.csv")
OUT_ASSIGN_AUDIT   <- file.path(DIR_VF, "vf_module_assignment_audit.csv")
OUT_UNASSIGNED     <- file.path(DIR_VF, "gene_module_unassigned_review.csv")
OUT_CATEGORY_DEFS   <- file.path(DIR_VF, "vf_category_definitions.csv")
OUT_QC_REPORT      <- file.path(DIR_VF, "vf_module_qc_report.txt")
OUT_NOTES          <- file.path(DIR_VF, "vf_module_definition_notes.md")
OUT_PLOT_COUNTS    <- file.path(DIR_PLOTS_VF, "module_gene_counts.png")
OUT_PLOT_PREV      <- file.path(DIR_PLOTS_VF, "module_prevalence_by_status.png")
OUT_PLOT_CATEGORY_COMPOSITION <- file.path(DIR_PLOTS_VF, "vf_category_composition_by_status.png")
OUT_PLOT_ASSIGNMENT_CONFIDENCE <- file.path(DIR_PLOTS_VF, "vf_module_assignment_confidence.png")

ensure_dir(DIR_VF)
ensure_dir(DIR_PLOTS_VF)

# ==============================================================================
# 2. LOAD INPUTS
# ==============================================================================
if (!file.exists(FILE_VF_READY)) stop("Missing ", FILE_VF_READY, ". Run 22_vf_build_analysis_dataset.R first.")
if (!file.exists(FILE_GENE_MAP)) stop("Missing ", FILE_GENE_MAP, ". Run 04_gene_breakdown.R first.")

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
      "%s is older than %s. Re-run 22_vf_build_analysis_dataset.R before script 26 so module outputs use the current VF matrix.\n  %s: %s\n  %s: %s",
      target_label, upstream_label, target_label, format(target_mtime),
      upstream_label, format(upstream_mtime)
    ))
  }
}

stop_if_stale(FILE_VF_READY, FILE_VF_PA_RAW, "vf_analysis_ready.csv", "vf_pa_all.csv")

vf_ready <- read_csv(FILE_VF_READY, show_col_types = FALSE) %>%
  prefer_primary_uti_status() %>%
  apply_manual_sample_curation(context = "26_vf_ready") %>%
  filter_primary_genomics() %>%
  mutate(Participant_id = as.character(Participant_id), tp_lab = normalise_timepoint_preserve_events(tp_lab))

gene_map <- read_csv(FILE_GENE_MAP, show_col_types = FALSE) %>%
  mutate(Gene = as.character(Gene), Category = coalesce(as.character(Category), "Unassigned"))

gene_cols <- canonical_vf_gene_cols(names(vf_ready), vf_pa_file = FILE_VF_PA_RAW)
raw_gene_cols <- canonical_vf_gene_cols(vf_pa_file = FILE_VF_PA_RAW)
meta_cols <- intersect(
  c(
    "Participant_id", "tp_lab", "Episode_ID", "Collection_Date",
    "Infection_Status", "Batch", "Status_Confidence_epi", "Sx_source_epi",
    "UTI_Label", "Urine_collection_method", "ST", "vf_count_total",
    "total_vf_count_all", "total_vf_count_curated",
    "total_vf_count_upec_candidate", "total_vf_count_unassigned",
    "low_confidence_count", "is_ecoli", "n_timepoints",
    "analysis_include_primary", "analysis_exclusion_reason",
    "duplicate_role", "duplicate_of_participant_id", "duplicate_of_tp_lab",
    "allow_secondary_duplicate_qc", "duplicate_use_note",
    "genomics_expected_include", "genomics_exclusion_reason"
  ),
  names(vf_ready)
)

msg("Loaded VF dataset: %d episodes, %d gene columns, %d participants",
    nrow(vf_ready), length(gene_cols), n_distinct(vf_ready$Participant_id))
msg("Loaded gene_map: %d genes", nrow(gene_map))
if (length(raw_gene_cols) > 0) {
  msg("Raw VF P/A matrix contains %d gene columns", length(raw_gene_cols))
}

# ==============================================================================
# 3. DEFINE MODULE ASSIGNMENT RULES
# ==============================================================================
# Each rule: list(module_id, system_name, broad_module, regex_pattern,
#                 confidence, upec_candidate, notes)
# Rules are applied in order; first match wins for primary assignment.

module_rules <- list(
  # --- Adhesion modules ---
  list("adhesion_type1_fimbriae", "Type 1 fimbriae (fim)", "Adhesion",
       "^fim[A-Z]$|^fim[A-Z][0-9]?$|^fimB$|^fimE$|^fimI$",
       "High", TRUE, "Mannose-sensitive adhesion; key UPEC colonisation factor"),
  list("adhesion_p_fimbriae", "P fimbriae (pap)", "Adhesion",
       "^pap[A-Z]$|^papX$",
       "High", TRUE, "Pyelonephritis-associated; classic UPEC marker"),
  list("adhesion_s_fimbriae", "S fimbriae (sfa)", "Adhesion",
       "^sfa[A-Z]$|^sfaX$|^sfaY$",
       "High", TRUE, "S fimbriae; sialic acid binding"),
  list("adhesion_f1c_fimbriae", "F1C fimbriae (foc)", "Adhesion",
       "^foc[A-Z]$",
       "High", TRUE, "F1C fimbriae; related to S fimbriae"),
  list("adhesion_afa_dr", "Afa/Dr adhesins", "Adhesion",
       "^afa[A-Z]|^dra[A-Z]|^draE|^draP$|^daaF$",
       "High", TRUE, "Afa/Dr family; associated with cystitis and pyelonephritis"),
  list("adhesion_ecp", "E. coli common pilus (ecp)", "Adhesion",
       "ecpR$|ecpA$|ecpB$|ecpC$|ecpD$|ecpE$|^ykgK|^yag[ZYWXV]",
       "High", FALSE, "Common pilus; widespread in E. coli"),
  list("adhesion_fae", "F4 fimbriae (fae)", "Adhesion",
       "^fae[A-Z]$",
       "Moderate", FALSE, "F4 fimbriae; typically animal-associated"),
  list("adhesion_lpf", "Long polar fimbriae (lpf)", "Adhesion",
       "^lpf[ABC]$",
       "Moderate", FALSE, "Long polar fimbriae"),
  list("adhesion_other", "Other adhesins (fdeC, cswA, sinH)", "Adhesion",
       "^fdeC$|^cswA$|^sinH$",
       "Moderate", FALSE, "Miscellaneous adhesins"),

  # --- Biofilm ---
  list("biofilm_curli", "Curli (csg)", "Biofilm",
       "^csg[ABCDEFG]$",
       "High", FALSE, "Curli fibres; biofilm formation and surface adhesion"),

  # --- Iron acquisition modules ---
  list("iron_aerobactin", "Aerobactin (iuc/iut)", "Iron acquisition",
       "^iuc[ABCD]$|^iutA$",
       "High", TRUE, "Aerobactin siderophore; key UPEC iron acquisition"),
  list("iron_yersiniabactin", "Yersiniabactin (ybt/irp/fyuA)", "Iron acquisition",
       "^ybt[AEPQSTUX]$|^irp[12]$|^fyuA$",
       "High", TRUE, "Yersiniabactin; HPI-associated siderophore"),
  list("iron_salmochelin", "Salmochelin (iro)", "Iron acquisition",
       "^iro[BCDNE]$",
       "High", TRUE, "Salmochelin; glucosylated enterobactin"),
  list("iron_enterobactin", "Enterobactin (ent/fep)", "Iron acquisition",
       "^ent[ABCDEF]$|^entS$|^fep[ABCDG]$|^fes$",
       "High", FALSE, "Enterobactin; conserved core E. coli siderophore"),
  list("iron_heme_chu", "Heme uptake (chu)", "Iron acquisition",
       "^chu[ASTUWXYV]$",
       "High", TRUE, "Heme/hemin uptake system"),
  list("iron_heme_shu", "Heme uptake (shu)", "Iron acquisition",
       "^shu[ASTXY]$",
       "Moderate", TRUE, "Shigella heme uptake homologue"),

  # --- Toxins ---
  list("toxin_hemolysin", "Alpha-hemolysin (hly)", "Toxin",
       "^hly[ABCD]$|^hlyC$",
       "High", TRUE, "Alpha-hemolysin; pore-forming toxin"),
  list("toxin_cnf", "CNF1 (cnf1)", "Toxin",
       "^cnf1$",
       "High", TRUE, "Cytotoxic necrotising factor 1"),
  list("toxin_sat", "Sat autotransporter", "Toxin",
       "^sat$",
       "High", TRUE, "Secreted autotransporter toxin; UPEC-associated"),
  list("toxin_vat", "Vat autotransporter", "Toxin",
       "^vat$",
       "High", TRUE, "Vacuolating autotransporter toxin"),
  list("toxin_east1", "EAST1 (astA)", "Toxin",
       "^astA$",
       "Moderate", FALSE, "Enteroaggregative heat-stable toxin"),
  list("toxin_pic", "Pic mucinase", "Toxin",
       "^pic$",
       "Moderate", TRUE, "Mucinase autotransporter"),
  list("toxin_senB", "SenB enterotoxin", "Toxin",
       "^senB$",
       "Moderate", FALSE, "Plasmid-encoded enterotoxin"),
  list("toxin_espP", "EspP serine protease", "Toxin",
       "^espP$",
       "Moderate", FALSE, "EHEC serine protease homologue"),
  list("toxin_sigA", "SigA autotransporter", "Toxin",
       "^sigA$",
       "Moderate", FALSE, "Shigella IgA protease homologue"),
  list("toxin_pet", "Pet autotransporter", "Toxin",
       "^pet$",
       "Low", FALSE, "Plasmid-encoded toxin"),

  # --- Capsule / Surface ---
  list("capsule_kps", "Group 2 capsule (kps)", "Capsule/Surface",
       "^kps[DMT]$",
       "High", TRUE, "Group 2 capsule; serum resistance"),
  list("surface_lps", "LPS/O-antigen (gtr/fcl/gmd)", "Capsule/Surface",
       "^gtr[AB]$|^gtrII$|^fcl$|^gmd$",
       "Moderate", FALSE, "LPS modification / O-antigen biosynthesis"),

  # --- Invasion/Evasion ---
  list("invasion_ibeA", "IbeA invasin", "Invasion/Evasion",
       "^ibeA$",
       "High", TRUE, "Invasion of brain endothelial cells"),
  list("invasion_ompA", "OmpA outer membrane protein", "Invasion/Evasion",
       "^ompA$",
       "High", FALSE, "Major OMP; multiple roles including evasion"),
  list("invasion_tcpC", "TcpC TIR-domain protein", "Invasion/Evasion",
       "^tcpC$",
       "High", TRUE, "TLR signalling inhibitor; immune evasion"),
  list("invasion_aslA", "AslA arylsulfatase", "Invasion/Evasion",
       "^aslA$",
       "Moderate", FALSE, "Arylsulfatase; potential immune modulation"),
  list("invasion_pla", "Pla protease", "Invasion/Evasion",
       "^pla$",
       "Low", FALSE, "Plasminogen activator"),

  # --- Secretion systems ---
  list("secretion_t2ss", "Type II secretion (gsp)", "Secretion system",
       "^gsp[CDEFGHIJKLM]$",
       "High", FALSE, "General secretion pathway"),

  # --- Flagella / Motility ---
  list("motility_flagella", "Flagella/chemotaxis (fli/flg/flh/che/mot)", "Motility",
       "^fli[AGIMPQZ]$|^flg[BCDGHI]$|^flh[ACD]$|^che[BRWYZ]$|^mot[A]$",
       "High", FALSE, "Flagellar motility and chemotaxis"),

  # --- T6SS ---
  list("secretion_t6ss", "Type VI secretion (hcp/hsi/vip)", "Secretion system",
       "^hcp|^hsi[BC]|^vipB$|^vipA$|^mglB$|^vipB/mglB$|hsiB1/vipA|hsiC1/vipB",
       "Moderate", FALSE, "Type VI secretion system components"),

  # --- Prophage/LEE-associated esp genes ---
  list("prophage_esp", "Prophage/LEE effectors (esp)", "Prophage/mobile",
       "^esp[XLRY][0-9]?$",
       "Low", FALSE, "LEE/prophage-associated effector-like genes"),

  # --- Miscellaneous ---
  list("misc_metabolism", "Metabolic/housekeeping (kdsA, galU, etc.)", "Metabolism",
       "^kdsA$|^galU$|^rfaD$|^lpx[CD]$|^gmhA|^lpcA$|^luxS$|^sodB$",
       "Low", FALSE, "Metabolic genes sometimes found in VFDB; not classic VF"),
  list("misc_IlpA", "IlpA lipoprotein", "Other",
       "^IlpA$",
       "Low", FALSE, "Lipoprotein; function in UPEC unclear"),
  list("misc_htpB", "HtpB chaperonin", "Other",
       "^htpB$",
       "Low", FALSE, "Heat-shock protein; Legionella VF homologue")
)

msg("Defined %d module assignment rules", length(module_rules))

# ==============================================================================
# 4. APPLY GENE-TO-MODULE CURATION
# ==============================================================================
# Build the module map by applying rules to all genes in gene_map
all_genes <- unique(c(gene_map$Gene, gene_cols))

assignments <- tibble()

for (g in all_genes) {
  matched <- FALSE
  g_clean <- tolower(trimws(g))

  for (rule in module_rules) {
    pattern <- rule[[4]]
    # Test against both original and cleaned gene name
    if (grepl(pattern, g, ignore.case = TRUE) || grepl(pattern, g_clean)) {
      assignments <- bind_rows(assignments, tibble(
        Gene = g,
        Gene_clean = g_clean,
        module_id = rule[[1]],
        system_name = rule[[2]],
        broad_module = rule[[3]],
        assignment_rule = pattern,
        assignment_confidence = rule[[5]],
        upec_score_candidate = rule[[6]],
        notes = rule[[7]],
        primary_assignment = TRUE
      ))
      matched <- TRUE
      break  # first match wins
    }
  }

  if (!matched) {
    assignments <- bind_rows(assignments, tibble(
      Gene = g,
      Gene_clean = g_clean,
      module_id = "unassigned",
      system_name = "Unassigned",
      broad_module = "Unassigned",
      assignment_rule = NA_character_,
      assignment_confidence = "Unassigned",
      upec_score_candidate = FALSE,
      notes = "No rule matched",
      primary_assignment = TRUE
    ))
  }
}

# Add existing category from gene_map
assignments <- assignments %>%
  left_join(gene_map %>% select(Gene, Category, Subcategory), by = "Gene") %>%
  mutate(
    Category = coalesce(Category, "Not in gene_map"),
    Subcategory = coalesce(Subcategory, "Not in gene_map"),
    in_vf_matrix = Gene %in% gene_cols,
    in_raw_vf_pa = if (length(raw_gene_cols) > 0) Gene %in% raw_gene_cols else NA,
    needs_manual_review = assignment_confidence %in% c("Low", "Unassigned"),
    gene_family = str_extract(Gene_clean, "^[a-z]+"),
    module_score_weight = 1L
  )

msg("Module assignments: %d genes total, %d in VF matrix",
    nrow(assignments), sum(assignments$in_vf_matrix))

# ==============================================================================
# 5. COMPUTE EPISODE-LEVEL MODULE PRESENCE
# ==============================================================================
# Only use genes that are (a) in the VF matrix AND (b) have primary assignment
active_assignments <- assignments %>%
  filter(in_vf_matrix, primary_assignment)

modules_detected <- active_assignments %>%
  filter(module_id != "unassigned") %>%
  pull(module_id) %>%
  unique() %>%
  sort()

msg("Modules with detected genes: %d", length(modules_detected))

# Build module presence and gene count per episode
module_ep <- vf_ready %>% select(all_of(meta_cols))

for (mod in modules_detected) {
  mod_genes <- active_assignments %>% filter(module_id == mod) %>% pull(Gene)
  matching  <- intersect(mod_genes, gene_cols)
  col_pres  <- paste0("mod_", mod, "_present")
  col_count <- paste0("mod_", mod, "_n_genes")
  if (length(matching) > 0) {
    gene_sums <- rowSums(vf_ready[, matching, drop = FALSE], na.rm = TRUE)
    module_ep[[col_pres]]  <- as.integer(gene_sums > 0)
    module_ep[[col_count]] <- as.integer(gene_sums)
  } else {
    module_ep[[col_pres]]  <- 0L
    module_ep[[col_count]] <- 0L
  }
}

# Compute summary columns
pres_cols <- grep("_present$", names(module_ep), value = TRUE)

# Calculate n_modules_present
module_ep$n_modules_present <- rowSums(module_ep[, pres_cols, drop = FALSE], na.rm = TRUE)

# Calculate n_upec_modules_present
upec_mods <- active_assignments %>%
  filter(upec_score_candidate, module_id != "unassigned") %>%
  pull(module_id) %>% unique()
upec_pres_cols <- pres_cols[gsub("^mod_|_present$", "", pres_cols) %in% upec_mods]

if (length(upec_pres_cols) > 0) {
  module_ep$n_upec_modules_present <- rowSums(module_ep[, upec_pres_cols, drop = FALSE], na.rm = TRUE)
} else {
  module_ep$n_upec_modules_present <- 0L
}

msg("Module episode matrix: %d rows × %d columns", nrow(module_ep), ncol(module_ep))

# ==============================================================================
# 6. BUILD MODULE SUMMARY TABLE
# ==============================================================================
mod_summary <- tibble()
for (mod in modules_detected) {
  info <- active_assignments %>% filter(module_id == mod) %>% slice(1)
  col_pres <- paste0("mod_", mod, "_present")
  col_count <- paste0("mod_", mod, "_n_genes")
  n_genes <- active_assignments %>% filter(module_id == mod) %>% nrow()

  ep_present <- module_ep[[col_pres]]
  status_vec <- module_ep$Infection_Status

  mod_summary <- bind_rows(mod_summary, tibble(
    module_id = mod,
    system_name = info$system_name,
    broad_module = info$broad_module,
    n_genes_in_module = n_genes,
    n_episodes_present = sum(ep_present, na.rm = TRUE),
    pct_episodes_present = round(100 * sum(ep_present, na.rm = TRUE) / nrow(module_ep), 1),
    n_Not_UTI_present = sum(ep_present[status_vec == "Not_UTI"], na.rm = TRUE),
    n_UTI_present  = sum(ep_present[status_vec == "UTI"], na.rm = TRUE),
    upec_score_candidate = info$upec_score_candidate,
    assignment_confidence = info$assignment_confidence
  ))
}

# Category definitions make the biological grouping denominator explicit for
# figure captions, manuscripts, and manual review.
category_definitions <- assignments %>%
  mutate(
    category_label = format_category_label(Category),
    category_definition = case_when(
      category_label == "Adhesion / fimbriae" ~ "Adhesins, fimbriae, and colonisation-associated surface structures.",
      category_label == "Iron acquisition" ~ "Siderophore and heme/iron uptake systems.",
      category_label == "Toxins" ~ "Secreted or cell-associated toxin genes.",
      category_label == "Capsule / surface structures" ~ "Capsule, surface polysaccharide, or outer-surface-associated genes.",
      category_label == "Invasion / immune evasion" ~ "Genes plausibly involved in invasion, serum resistance, or immune evasion.",
      category_label == "Unassigned VFDB hits not in gene map" ~ "Detected VFDB matrix genes absent from the curated gene_map input.",
      category_label == "Unassigned VFDB hits" ~ "Detected VFDB genes retained in the matrix but not assigned to a curated biological category.",
      TRUE ~ "Descriptive VFDB-derived category from gene_map."
    )
  ) %>%
  group_by(category_label, category_definition) %>%
  summarise(
    n_genes = n(),
    n_genes_in_vf_matrix = sum(in_vf_matrix, na.rm = TRUE),
    n_upec_score_candidates = sum(upec_score_candidate, na.rm = TRUE),
    example_genes = paste(head(sort(unique(Gene)), 12), collapse = ";"),
    interpretation_note = "Descriptive biological grouping only; not a validated causal virulence score.",
    .groups = "drop"
  ) %>%
  arrange(category_label)

# ==============================================================================
# 7. BUILD UNASSIGNED REVIEW TABLE
# ==============================================================================
unassigned_genes <- assignments %>%
  filter(module_id == "unassigned", in_vf_matrix)

if (nrow(unassigned_genes) > 0) {
  unassigned_review <- unassigned_genes %>%
    rowwise() %>%
    mutate(
      detected_in_n_episodes = if (Gene %in% gene_cols) sum(vf_ready[[Gene]], na.rm = TRUE) else 0L,
      detected_in_pct_episodes = round(100 * detected_in_n_episodes / nrow(vf_ready), 1)
    ) %>%
    ungroup() %>%
    select(Gene, Gene_clean, Category, gene_family,
           detected_in_n_episodes, detected_in_pct_episodes,
           notes) %>%
    arrange(desc(detected_in_n_episodes))
} else {
  unassigned_review <- tibble(Gene = character(), notes = character())
}

# ==============================================================================
# 8. BUILD ASSIGNMENT AUDIT TABLE
# ==============================================================================
# This table is intentionally broader than gene_module_map.csv. It records
# whether each gene is present in the raw VF P/A matrix, the canonical
# analysis-ready matrix, and the curated map, so stale or partial curation is
# visible instead of silently changing downstream denominators.
assignment_audit <- assignments %>%
  mutate(
    in_gene_map = Gene %in% gene_map$Gene,
    in_current_raw_vf_pa = if (length(raw_gene_cols) > 0) Gene %in% raw_gene_cols else NA,
    in_current_vf_ready = Gene %in% gene_cols,
    detected_in_n_episodes = sapply(Gene, function(g) {
      if (g %in% gene_cols) sum(vf_ready[[g]], na.rm = TRUE) else 0L
    }),
    detected_in_pct_episodes = round(100 * detected_in_n_episodes / nrow(vf_ready), 1),
    audit_flag = case_when(
      !in_gene_map & in_current_vf_ready ~ "In VF-ready matrix but absent from gene_map",
      !is.na(in_current_raw_vf_pa) & in_current_raw_vf_pa & !in_current_vf_ready ~ "In raw VF matrix but absent from vf_analysis_ready",
      module_id == "unassigned" ~ "Unassigned module",
      needs_manual_review ~ "Manual review recommended",
      TRUE ~ "OK"
    )
  ) %>%
  arrange(desc(audit_flag != "OK"), audit_flag, module_id, Gene)

# ==============================================================================
# 9. QC REPORT
# ==============================================================================
qc <- character()
qc_add <- function(...) qc <<- c(qc, sprintf(...))

qc_add("=== VF MODULE QC REPORT ===")
qc_add("Timestamp: %s", format(Sys.time()))
qc_add("")
qc_add("--- INPUT ---")
qc_add("Episodes: %d", nrow(vf_ready))
qc_add("Participants: %d", n_distinct(vf_ready$Participant_id))
qc_add("UTI: %d, Not_UTI: %d",
       sum(vf_ready$Infection_Status == "UTI", na.rm = TRUE),
       sum(vf_ready$Infection_Status == "Not_UTI", na.rm = TRUE))
qc_add("VF gene columns in matrix: %d", length(gene_cols))
if (length(raw_gene_cols) > 0) {
  qc_add("Raw VF P/A gene columns: %d", length(raw_gene_cols))
  qc_add("Genes in raw VF P/A but not in vf_analysis_ready: %d",
         length(setdiff(raw_gene_cols, gene_cols)))
  qc_add("Genes in raw VF P/A but not in gene_map: %d",
         length(setdiff(raw_gene_cols, gene_map$Gene)))
}
qc_add("Genes in gene_map: %d", nrow(gene_map))
qc_add("")
qc_add("--- ASSIGNMENT ---")
qc_add("Total genes processed: %d", nrow(assignments))
qc_add("Genes in VF matrix: %d", sum(assignments$in_vf_matrix))
conf_tbl <- table(assignments$assignment_confidence[assignments$in_vf_matrix])
for (nm in names(conf_tbl)) qc_add("  %s: %d", nm, conf_tbl[[nm]])
qc_add("Modules created: %d", length(modules_detected))
qc_add("Unassigned genes (in matrix): %d", nrow(unassigned_review))
qc_add("Assignment audit flags needing review: %d",
       sum(assignment_audit$audit_flag != "OK", na.rm = TRUE))
qc_add("UPEC-candidate modules: %d",
       sum(mod_summary$upec_score_candidate, na.rm = TRUE))
qc_add("WARNING: Unassigned/low-confidence genes are retained separately and should not be interpreted as curated UPEC-relevant genes.")
qc_add("")
qc_add("--- OUTPUT ---")
qc_add("Module episode matrix: %d rows", nrow(module_ep))
qc_add("Module presence columns: %d", length(pres_cols))

# Duplicate checks
dup_ep <- module_ep %>% dplyr::count(Participant_id, tp_lab) %>% filter(n > 1)
qc_add("Duplicate Participant_id × tp_lab in module output: %d", nrow(dup_ep))

dup_map <- assignments %>% filter(primary_assignment, in_vf_matrix) %>%
  dplyr::count(Gene) %>% filter(n > 1)
qc_add("Genes with multiple primary assignments: %d", nrow(dup_map))

writeLines(qc, OUT_QC_REPORT)
msg("QC report written to %s", OUT_QC_REPORT)

gap_report <- assignment_audit %>%
  summarise(
    vf_gene_columns = sum(in_current_vf_ready, na.rm = TRUE),
    genes_in_gene_map = sum(in_gene_map, na.rm = TRUE),
    genes_in_matrix_not_gene_map = sum(!in_gene_map & in_current_vf_ready, na.rm = TRUE),
    unassigned_genes_in_matrix = sum(module_id == "unassigned" & in_current_vf_ready, na.rm = TRUE),
    low_confidence_genes_in_matrix = sum(assignment_confidence == "Low" & in_current_vf_ready, na.rm = TRUE),
    audit_flags_need_review = sum(audit_flag != "OK", na.rm = TRUE)
  )
write_csv(gap_report, file.path(DIR_VF, "vf_gene_annotation_gap_report.csv"))
writeLines(
  c(
    "VF gene annotation gap report",
    sprintf("Generated: %s", format(Sys.time())),
    sprintf("VF matrix gene columns: %d", gap_report$vf_gene_columns),
    sprintf("Genes in matrix not in gene_map: %d", gap_report$genes_in_matrix_not_gene_map),
    sprintf("Unassigned genes in matrix: %d", gap_report$unassigned_genes_in_matrix),
    sprintf("Low-confidence genes in matrix: %d", gap_report$low_confidence_genes_in_matrix),
    sprintf("Assignment audit flags needing review: %d", gap_report$audit_flags_need_review),
    "Interpretation: total_vf_count_all includes all detected VFDB calls. Curated and UPEC-candidate counts should be used for biological interpretation."
  ),
  file.path(DIR_VF, "vf_gene_annotation_gap_report.txt")
)

# ==============================================================================
# 10. WRITE OUTPUTS
# ==============================================================================
notes_md <- c(
  "# VF Module Definition Notes",
  "",
  sprintf("Generated: %s", format(Sys.time())),
  "",
  "This file documents the rule-based biological module framework used by `26_vf_define_gene_modules.R`.",
  "Module presence is derived from the canonical selected-assembly VF gene matrix; this script does not change Abricate thresholds.",
  "",
  "## Interpretation",
  "",
  "- Modules are descriptive curation units, not validated predictors of UTI.",
  "- `upec_score_candidate` marks classic/plausible UPEC-associated systems for downstream exploratory scoring.",
  "- AMR/mobile-like or metabolic modules are retained for transparency but are not treated as true AMR analysis.",
  "- Genes with no defensible rule are retained as `unassigned` and excluded from default curated scores.",
  "",
  "## Assignment Rules",
  "",
  paste0(
    "- `", vapply(module_rules, `[[`, character(1), 1), "`: ",
    vapply(module_rules, `[[`, character(1), 2),
    " | broad module: ", vapply(module_rules, `[[`, character(1), 3),
    " | confidence: ", vapply(module_rules, `[[`, character(1), 5),
    " | UPEC candidate: ", vapply(module_rules, `[[`, logical(1), 6),
    " | regex: `", vapply(module_rules, `[[`, character(1), 4), "`"
  ),
  "",
  "## Review Outputs",
  "",
  "- `vf_module_assignment_audit.csv` records map/raw-matrix/canonical-matrix representation for every processed gene.",
  "- `gene_module_unassigned_review.csv` lists detected genes that remain unassigned and need manual review if they are biologically important."
)

write_csv(assignments, OUT_MODULE_MAP)
write_csv(module_ep, OUT_MODULE_EP)
write_csv(mod_summary, OUT_MODULE_SUMMARY)
write_csv(assignment_audit, OUT_ASSIGN_AUDIT)
write_csv(unassigned_review, OUT_UNASSIGNED)
write_csv(category_definitions, OUT_CATEGORY_DEFS)
writeLines(notes_md, OUT_NOTES)

msg("Wrote gene_module_map.csv (%d rows)", nrow(assignments))
msg("Wrote vf_module_presence_by_episode.csv (%d rows)", nrow(module_ep))
msg("Wrote vf_module_summary.csv (%d modules)", nrow(mod_summary))
msg("Wrote vf_module_assignment_audit.csv (%d rows)", nrow(assignment_audit))
msg("Wrote gene_module_unassigned_review.csv (%d genes)", nrow(unassigned_review))
msg("Wrote vf_category_definitions.csv (%d categories)", nrow(category_definitions))
msg("Wrote vf_module_definition_notes.md")

append_denominator_summary(
  module_ep,
  "26_vf_define_gene_modules.R",
  "vf_module_presence_by_episode",
  "participant_timepoint",
  OUT_MODULE_EP,
  "Module matrix keeps unassigned/low-confidence genes separate from curated UPEC-candidate systems"
)

# ==============================================================================
# 11. PLOTS
# ==============================================================================

# --- Module gene counts bar plot ---
if (nrow(mod_summary) > 0) {
  p1 <- ggplot(mod_summary, aes(x = reorder(system_name, n_genes_in_module),
                                 y = n_genes_in_module, fill = broad_module)) +
    geom_col() +
    coord_flip() +
    labs(
      title = "Virulence factor genes per curated biological module",
      subtitle = "Module definitions are descriptive curation units, not validated UTI-causality scores",
      x = NULL,
      y = "Number of VF genes assigned to module",
      fill = "Broad module",
      caption = sprintf(
        "Data: %s and %s. Denominator: %d detected/curated VF genes assigned across %d modules. Level of analysis: gene-to-module curation. Categories and modules are descriptive biological groupings and should not be interpreted as validated causal virulence scores.",
        FILE_GENE_MAP, FILE_VF_READY, nrow(assignments), nrow(mod_summary)
      )
    ) +
    plot_theme_vf(base_size = 11)
  ggsave(OUT_PLOT_COUNTS, p1, width = 10, height = max(6, nrow(mod_summary) * 0.3), dpi = 300)
  msg("Saved %s", OUT_PLOT_COUNTS)
}

# --- Module prevalence by status ---
if (nrow(mod_summary) > 0) {
  prev_long <- mod_summary %>%
    select(module_id, system_name, broad_module, n_Not_UTI_present, n_UTI_present) %>%
    pivot_longer(cols = starts_with("n_"), names_to = "status", values_to = "n_present") %>%
    mutate(
      status = case_when(
        status == "n_Not_UTI_present" ~ "Not_UTI",
        status == "n_UTI_present" ~ "UTI",
      ),
      denom = case_when(
        status == "Not_UTI" ~ sum(vf_ready$Infection_Status == "Not_UTI", na.rm = TRUE),
        status == "UTI" ~ sum(vf_ready$Infection_Status == "UTI", na.rm = TRUE),
      ),
      pct = round(100 * n_present / denom, 1)
    ) %>%
    mutate(status = factor(status, levels = c("UTI", "Not_UTI")))

  p2 <- ggplot(prev_long, aes(x = reorder(system_name, pct), y = pct, fill = status)) +
    geom_col(position = "dodge") +
    coord_flip() +
    scale_fill_uti_status() +
    labs(
      title = "Virulence factor module prevalence across primary UTI status",
      subtitle = sprintf("Presence = at least one detected gene in module; UTI n=%d, Not_UTI n=%d",
                         sum(vf_ready$Infection_Status == "UTI", na.rm = TRUE),
                         sum(vf_ready$Infection_Status == "Not_UTI", na.rm = TRUE)),
      x = NULL,
      y = "Isolates with module present",
      fill = "Primary UTI status",
      caption = sprintf(
        "Data: %s. Denominator: %d VF/WGS-linked E. coli isolates from %d participants. Level of analysis: isolate-level module presence. Residents may contribute repeated isolates; module prevalence is descriptive and should not be read as causal evidence. UTI n=%d is small, Not_UTI is heterogeneous, and ST/lineage may confound module-status contrasts.",
        FILE_VF_READY, nrow(vf_ready), n_distinct(vf_ready$Participant_id),
        sum(vf_ready$Infection_Status == "UTI", na.rm = TRUE)
      )
    ) +
    plot_theme_vf(base_size = 11)
  ggsave(OUT_PLOT_PREV, p2, width = 12, height = max(7, nrow(mod_summary) * 0.3), dpi = 300)
  msg("Saved %s", OUT_PLOT_PREV)
}

# --- Category composition by status ---
cat_cols_ready <- grep("^cat_", names(vf_ready), value = TRUE)
if (length(cat_cols_ready) > 0) {
  cat_comp <- vf_ready %>%
    filter(!is.na(Infection_Status), Infection_Status %in% c("UTI", "Not_UTI")) %>%
    mutate(Infection_Status = status_for_plot(Infection_Status)) %>%
    select(Participant_id, Infection_Status, all_of(cat_cols_ready)) %>%
    pivot_longer(cols = all_of(cat_cols_ready),
                 names_to = "category_col", values_to = "gene_count") %>%
    mutate(category_label = recode(
      category_col,
      "cat_Adhesion_Fimbriae" = "Adhesion / fimbriae",
      "cat_Iron_acquisition" = "Iron acquisition",
      "cat_Toxins" = "Toxins",
      "cat_Capsule_Surface" = "Capsule / surface structures",
      "cat_Invasion_Evasion" = "Invasion / immune evasion",
      "cat_Unassigned" = "Unassigned VFDB hits",
      "cat_Unassigned_matrix" = "Unassigned VFDB hits not in gene map",
      .default = gsub("_", " ", gsub("^cat_", "", category_col))
    )) %>%
    group_by(Infection_Status, category_label) %>%
    summarise(mean_count = mean(gene_count, na.rm = TRUE), .groups = "drop") %>%
    group_by(Infection_Status) %>%
    mutate(prop_of_mean_category_burden = mean_count / sum(mean_count)) %>%
    ungroup()

  p_cat_comp <- ggplot(cat_comp, aes(x = Infection_Status,
                                     y = prop_of_mean_category_burden,
                                     fill = category_label)) +
    geom_col(width = 0.65, colour = "white", linewidth = 0.2) +
    scale_y_continuous(labels = scales::percent) +
    labs(
      title = "Virulence factor category profiles across primary UTI status",
      subtitle = "Category burden represents the number of detected genes per isolate within each curated VF category",
      x = "Primary UTI status",
      y = "Share of mean category burden",
      fill = "VF category",
      caption = sprintf(
        "Data: %s. Denominator: %d VF/WGS-linked E. coli isolates from %d participants. Level of analysis: isolate-level category composition summarized by primary UTI status. Residents may contribute repeated isolates. Categories are descriptive biological groupings and should not be interpreted as validated causal virulence scores; UTI n=%d is small and ST/lineage may confound interpretation.",
        FILE_VF_READY, nrow(vf_ready), n_distinct(vf_ready$Participant_id),
        sum(vf_ready$Infection_Status == "UTI", na.rm = TRUE)
      )
    ) +
    plot_theme_vf(base_size = 11)

  ggsave(OUT_PLOT_CATEGORY_COMPOSITION, p_cat_comp, width = 8, height = 5.8, dpi = 300)
  msg("Saved %s", OUT_PLOT_CATEGORY_COMPOSITION)
}

# --- Module assignment confidence and gene-map coverage diagnostic ---
if (nrow(assignment_audit) > 0) {
  assignment_confidence_plot <- assignment_audit %>%
    mutate(
      assignment_confidence = factor(assignment_confidence,
                                     levels = c("High", "Moderate", "Low", "Unassigned")),
      gene_map_status = case_when(
        in_current_vf_ready %in% TRUE & !(in_gene_map %in% TRUE) ~ "Detected in VF matrix; absent from gene_map",
        in_gene_map %in% TRUE ~ "Represented in gene_map",
        TRUE ~ "Not represented in gene_map"
      )
    ) %>%
    count(assignment_confidence, gene_map_status, name = "n") %>%
    filter(!is.na(assignment_confidence), n > 0)

  if (nrow(assignment_confidence_plot) > 0) {
    gap_text <- if (exists("gap_report") && nrow(gap_report) > 0) {
      sprintf(
        "Current VF matrix: %d genes; %d absent from gene_map; %d unassigned; %d low-confidence",
        gap_report$vf_gene_columns[1],
        gap_report$genes_in_matrix_not_gene_map[1],
        gap_report$unassigned_genes_in_matrix[1],
        gap_report$low_confidence_genes_in_matrix[1]
      )
    } else {
      "Assignment confidence is summarized from vf_module_assignment_audit.csv"
    }

    p_assign <- ggplot(assignment_confidence_plot,
                       aes(x = assignment_confidence, y = n, fill = gene_map_status)) +
      geom_col(width = 0.68, colour = "white", linewidth = 0.25) +
      geom_text(aes(label = n), position = position_stack(vjust = 0.5),
                size = 3.2, colour = "white", check_overlap = TRUE) +
      scale_fill_manual(values = c(
        "Represented in gene_map" = "#0072B2",
        "Detected in VF matrix; absent from gene_map" = "#D55E00",
        "Not represented in gene_map" = "grey60"
      )) +
      labs(
        title = "VF module assignment confidence and gene-map coverage",
        subtitle = gap_text,
        x = "Module assignment confidence",
        y = "Number of VF genes / assignment rows",
        fill = "Gene-map coverage",
        caption = sprintf(
          "Data: %s and %s. Level of analysis: gene-to-module curation diagnostic. This plot audits the annotation basis for VF modules; low-confidence, unassigned, and gene-map-absent VFDB hits should not be interpreted as validated causal virulence modules.",
          OUT_ASSIGN_AUDIT, file.path(DIR_VF, "vf_gene_annotation_gap_report.csv")
        )
      ) +
      plot_theme_vf(base_size = 11)

    ggsave(OUT_PLOT_ASSIGNMENT_CONFIDENCE, p_assign, width = 8.5, height = 5.4, dpi = 300)
    msg("Saved %s", OUT_PLOT_ASSIGNMENT_CONFIDENCE)
  }
}

# ==============================================================================
# 12. FINAL CONSOLE SUMMARY
# ==============================================================================
cat("\n")
for (line in qc) cat(line, "\n")
cat("\n")
msg("✓ 26_vf_define_gene_modules.R complete.")

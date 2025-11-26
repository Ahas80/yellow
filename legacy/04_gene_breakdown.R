#!/usr/bin/env Rscript
# =============================================================
# 04_gene_breakdown.R  (updated: adds focus-gene Fisher test)
# -------------------------------------------------------------
# Reads : results/vf_hits_all.rds   (long VFDB-hit table from 02)
# Writes:
#   results/annotated_gene_table.csv
#   results/per_sample_category_counts.csv
#   results/nitrate_presence_matrix.csv
#   results/plots/nitrate_upset.png
#   results/plots/vf_focus_prevalence_by_status.png    (if status_map exists)
#   results/vf_focus_prevalence_by_status.csv          (if status_map exists)
#   results/diff_focus_genes_UTI_vs_ASB_fisher.csv     (NEW)
# Notes:
#  • Safe GENE→Gene rename; adds tp_lab only if missing
#  • DISTINCT per-sample counts (no assembler double-counts)
#  • Optional status-stratified prevalence & Fisher test for *focus* genes
#  • ggplot2 3.5+ legend alignment (no deprecated args)
# =============================================================

## ---------- 1 · libraries, dirs, theme --------------------------------------
suppressPackageStartupMessages({
  library(dplyr);   library(tidyr);      library(readr)
  library(tibble);  library(ggplot2);    library(ComplexUpset)
  library(stringr); library(purrr);      library(forcats)
})

dir.create("results/plots", recursive = TRUE, showWarnings = FALSE)
# --- PATCH: load or create a gene_map -----------------------------------------
# Where to load/save the mapping table
gene_map_path <- "results/gene_map.csv"

# If a mapping already exists on disk, use it; otherwise create a sensible default
if (file.exists(gene_map_path)) {
  gene_map <- readr::read_csv(gene_map_path, show_col_types = FALSE) %>%
    dplyr::mutate(
      Gene = as.character(Gene),
      Category = dplyr::coalesce(as.character(Category), "Unassigned"),
      Subcategory = dplyr::coalesce(as.character(Subcategory), Category)
    )
} else {
  message("↪  No results/gene_map.csv found; creating a starter map…")
  
  # Build a starter mapping from genes seen in vf_hits_all
  genes_seen <- vf_hits_all %>%
    dplyr::rename(Gene = dplyr::any_of(c("GENE","gene","Gene"))) %>%
    dplyr::mutate(Gene = trimws(as.character(Gene))) %>%
    dplyr::distinct(Gene) %>%
    dplyr::filter(!is.na(Gene), Gene != "")
  
  # Heuristic categories (edit later in the CSV if you prefer)
  gene_map <- genes_seen %>%
    dplyr::mutate(
      Category = dplyr::case_when(
        stringr::str_detect(Gene, regex("^(fim|fml|pil|foc|sfa|pap|afa|dra|cfa)", TRUE)) ~ "Adhesion/Fimbriae",
        stringr::str_detect(Gene, regex("^(kps|kfi|neu|ugd|rmpA|caps|wzx|wzy)", TRUE))    ~ "Capsule/Surface",
        stringr::str_detect(Gene, regex("^(iut|iuc|iro|irp|fyuA|chu|fep|ent|fec|ybt)", TRUE)) ~ "Iron acquisition/Siderophores",
        stringr::str_detect(Gene, regex("^(hly|cnf|sat|vat|cdt|astA|subAB|stx|lt|st)", TRUE)) ~ "Toxins",
        stringr::str_detect(Gene, regex("^(omp|iss|ibe|tra|usp|malX)", TRUE))            ~ "Invasion/Evasion",
        stringr::str_detect(Gene, regex("^(bla|qnr|aac|aph|aad|erm|cat|tet|sul|dfr|mcr|gyrA|parC)", TRUE)) ~ "AMR",
        TRUE ~ "Unassigned"
      ),
      Subcategory = Category
    )
  
  readr::write_csv(gene_map, gene_map_path)
  message("✓  Starter mapping written to ", gene_map_path,
          " — review/edit categories there any time and re-run.")
}
# -------------------------------------------------------------------------------
# ggplot2 3.5+ friendly legend alignment
theme_set(theme_minimal(base_size = 11))
theme_update(legend.title = element_text(hjust = 0),
             legend.text  = element_text(hjust = 0))

## ---------- 2 · helpers ------------------------------------------------------
# Uniform timepoint labels (‘T0’, ‘T1’, …, ‘Uricult’, ‘Unscheduled’)
tp_norm <- function(x) {
  tp_chr     <- as.character(x)
  is_uricult <- str_detect(tp_chr, regex("uricult", ignore_case = TRUE))
  tp_num     <- suppressWarnings(as.integer(str_extract(tp_chr, "\\d+")))
  tp_lab     <- dplyr::case_when(
    is_uricult     ~ "Uricult",
    !is.na(tp_num) ~ paste0("T", tp_num),
    TRUE           ~ "Unscheduled"
  )
  tp_levels <- c(paste0("T", sort(unique(tp_num[!is.na(tp_num)]))),
                 "Uricult", "Unscheduled")
  tibble(tp_lab = factor(tp_lab, levels = tp_levels), tp_num = tp_num)
}

# Canonicalise gene strings (case/space/punct insensitive)
canon <- function(x) str_to_lower(gsub("[^a-z0-9]+", "", x))

## ---------- 3 · load VFDB hits ----------------------------------------------
vf_hits_all <- readRDS("results/vf_hits_all.rds")   # from script 02

# --- PATCH: standardize 'Gene', add tp_lab once, map messy names via regex ----

# 1) Ensure we have a 'Gene' column and add tp_lab only if missing
hits <- vf_hits_all
gene_col <- intersect(c("GENE","gene","Gene"), names(hits))[1]
if (is.na(gene_col)) stop("Could not find a gene column among: GENE/gene/Gene")
names(hits)[names(hits) == gene_col] <- "Gene"

# add tp_lab/tp_num if not already present
if (!"tp_lab" %in% names(hits)) {
  hits <- dplyr::bind_cols(hits, tp_norm(hits$Timepoint))
}

# trim whitespace in Gene
hits$Gene <- trimws(hits$Gene)

# 2) Regex alias rules (first match wins)
alias_rules <- tibble::tribble(
  ~pattern,                                   ~Gene,
  "^bla[-_ ]?ctx[-_ ]?m[-_ ]?15(\\.[0-9]+)?$",  "blaCTX-M-15",
  "^bla[-_ ]?ctx[-_ ]?m[-_ ]?[0-9]+[a-zA-Z]?$", "blaCTX-M",
  "^bla[-_ ]?oxa[-_ ]?48$",                     "blaOXA-48",
  "^bla[-_ ]?oxa[-_ ]?181$",                    "blaOXA-181",
  "^bla[-_ ]?oxa[-_ ]?232$",                    "blaOXA-232",
  "^bla[-_ ]?oxa[-_ ]?244$",                    "blaOXA-244",
  "^aac\\(6'\\)[-_ ]?ib[-_ ]?cr[0-9]*$",        "aac(6')-Ib-cr",
  "^aac\\(6'\\)[-_ ]?ib($|[-_ ].*)",            "aac(6')-Ib",
  "^aac\\(3\\).*$",                              "aac(3)",
  "^qnrA[0-9]*$",                               "qnrA",
  "^qnrB[0-9]*$",                               "qnrB",
  "^qnrS[0-9]*$",                               "qnrS",
  "^dfrA[-_ ]?1(\\b|\\D|$)",                    "dfrA1",
  "^dfrA[-_ ]?17(\\b|\\D|$)",                   "dfrA17",
  "^dfrA[0-9]+.*$",                              "dfrA",
  "^mcr[-_ ]?([1-9]|10)(\\.[0-9]+)?$",          "mcr-\\1",
  "^gyrA(\\b|[_-].*)$",                         "gyrA",
  "^parC(\\b|[_-].*)$",                         "parC",
  "^sfa[/_-]?focde$",                           "sfa/focDE",
  "^afa[/_-]?drabc$",                           "afa/draBC",
  "^kpsm[-_ ]?ii$",                              "kpsM II",
  "^kpsmt[-_ ]?ii$",                             "kpsMT II",
  "^repfia.*$",                                  "repFIA",
  "^repfii.*$",                                  "repFII",
  "^is1(\\b|[_-].*)$",                           "IS1",
  "^is3(\\b|[_-].*)$",                           "IS3",
  "^is91(\\b|[_-].*)$",                          "IS91",
  "^iscr1(\\b|[_-].*)$",                         "ISCR1",
  "^tn21(\\b|[_-].*)$",                          "Tn21"
)

apply_aliases <- function(df, gene_col = "Gene") {
  x <- as.character(df[[gene_col]])
  seen <- rep(FALSE, length(x))
  for (i in seq_len(nrow(alias_rules))) {
    pat  <- alias_rules$pattern[i]
    repl <- alias_rules$Gene[i]
    idx  <- !seen & grepl(pat, x, ignore.case = TRUE, perl = TRUE)
    if (any(idx, na.rm = TRUE)) {
      x[idx] <- gsub(pat, repl, x[idx], ignore.case = TRUE, perl = TRUE)
      seen[idx] <- TRUE
    }
  }
  df[[gene_col]] <- x
  df
}

# 3) Apply aliases, then join to your basic gene_map (exact names)
hits <- apply_aliases(hits, "Gene")

annotated <- hits %>%
  dplyr::left_join(gene_map, by = "Gene") %>%
  dplyr::mutate(
    Category    = dplyr::coalesce(Category,    "Unassigned"),
    Subcategory = dplyr::coalesce(Subcategory, "Unassigned")
  )

# from here your script can continue as before, e.g.:
# readr::write_csv(annotated, "results/annotated_gene_table.csv")
# ... etc.
## ---------- 4 · (your existing) gene annotation ------------------------------
# =========================
# Canonical gene annotation
# =========================
gene_map <- tibble::tribble(
  ~Gene,          ~Category,                ~Subcategory,
  
  # ---- Virulence (focus + core) ----
  "afa/draBC",    "Virulence factor",       "Adhesin (Dr-binding fimbriae)",
  "bmaE",         "Virulence factor",       "Adhesin (M hemagglutinin)",
  "cvaC",         "Virulence factor",       "Bacteriocin (Colicin V production)",
  "f17A",         "Virulence factor",       "Adhesin (F17 fimbriae)",
  "fimG",         "Virulence factor",       "Type 1 fimbriae (assembly)",
  "fimH",         "Virulence factor",       "Adhesin (Type 1 fimbriae)",
  "focG",         "Virulence factor",       "Adhesin (F1C fimbriae)",
  "hlyA",         "Virulence factor",       "Toxin (alpha-hemolysin)",
  "ibeA",         "Virulence factor",       "Invasin (brain endothelium)",
  "iha",          "Virulence factor",       "Adhesin (iron-regulated)",
  "iroN",         "Virulence factor",       "Siderophore receptor (iron acquisition)",
  "kpsM II",      "Virulence factor",       "Capsule synthesis (Group II)",
  "kpsMT II",     "Virulence factor",       "Capsule transport (Group II)",
  "malX",         "Virulence factor",       "PAI-associated PTS",
  "neuA",         "Virulence factor",       "Capsule (sialic acid metabolism)",
  "ompT",         "Virulence factor",       "Outer membrane protease",
  "papA",         "Virulence factor",       "Adhesin (P fimbriae)",
  "papC",         "Virulence factor",       "P fimbriae assembly (usher)",
  "papEF",        "Virulence factor",       "P fimbriae accessory",
  "pic",          "Virulence factor",       "Protease (colonization)",
  "sfa/focDE",    "Virulence factor",       "Adhesin (S/F1C fimbriae)",
  "sfaS",         "Virulence factor",       "Adhesin (S fimbriae)",
  "sat",          "Virulence factor",       "Toxin (autotransporter)",
  "tcpC",         "Virulence factor",       "Immune modulation (TLR antagonist)",
  "traT",         "Virulence factor",       "Serum resistance (surface exclusion)",
  "tsh",          "Virulence factor",       "Adhesin (temperature-sensitive hemagglutinin)",
  "usp",          "Virulence factor",       "Uropathogenic-specific protein",
  "vat",          "Virulence factor",       "Toxin (vacuolating autotransporter)",
  "yfcV",         "Virulence factor",       "Adhesin (fimbrial-like protein)",
  "chuA",         "Virulence factor",       "Heme uptake receptor",
  "fyuA",         "Virulence factor",       "Yersiniabactin receptor",
  "iucD",         "Virulence factor",       "Aerobactin synthesis",
  
  # ---- Antibiotic resistance ----
  "blaOXA-48",    "Antibiotic resistance",  "Beta-lactamase (OXA-48-like)",
  "blaOXA-181",   "Antibiotic resistance",  "Beta-lactamase (OXA-48 variant)",
  "blaOXA-232",   "Antibiotic resistance",  "Beta-lactamase (OXA-48 variant)",
  "blaOXA-244",   "Antibiotic resistance",  "Beta-lactamase (OXA-48 variant)",
  "blaCTX-M-15",  "Antibiotic resistance",  "ESBL (cephalosporins)",
  "blaCTX-M",     "Antibiotic resistance",  "ESBL (CTX-M family)",
  "blaTEM",       "Antibiotic resistance",  "Beta-lactamase",
  "blaSHV",       "Antibiotic resistance",  "Beta-lactamase",
  "aac(6')-Ib-cr","Antibiotic resistance",  "PMQR / aminoglycoside",
  "aac(6')-Ib",   "Antibiotic resistance",  "Aminoglycoside acetyltransferase",
  "aac(3)",       "Antibiotic resistance",  "Aminoglycoside acetyltransferase (AAC3 family)",
  "aadA",         "Antibiotic resistance",  "Aminoglycoside adenyltransferase",
  "aph(3')-Ia",   "Antibiotic resistance",  "Aminoglycoside phosphotransferase",
  "aph(3')",      "Antibiotic resistance",  "Aminoglycoside phosphotransferase (family)",
  "mph(A)",       "Antibiotic resistance",  "Macrolide phosphotransferase",
  "dfrA1",        "Antibiotic resistance",  "Trimethoprim (DHFR A1)",
  "dfrA17",       "Antibiotic resistance",  "Trimethoprim (DHFR A17)",
  "dfrA",         "Antibiotic resistance",  "Trimethoprim (DHFR family)",
  "sul1",         "Antibiotic resistance",  "Sulfonamide",
  "sul2",         "Antibiotic resistance",  "Sulfonamide",
  "sul3",         "Antibiotic resistance",  "Sulfonamide",
  "tetA",         "Antibiotic resistance",  "Tetracycline efflux",
  "tetB",         "Antibiotic resistance",  "Tetracycline efflux",
  "catA1",        "Antibiotic resistance",  "Chloramphenicol acetyltransferase",
  "floR",         "Antibiotic resistance",  "Phenicols",
  "qnrA",         "Antibiotic resistance",  "PMQR (qnrA)",
  "qnrB",         "Antibiotic resistance",  "PMQR (qnrB)",
  "qnrS",         "Antibiotic resistance",  "PMQR (qnrS)",
  "gyrA",         "Antibiotic resistance",  "Quinolone (QRDR site)",
  "parC",         "Antibiotic resistance",  "Quinolone (QRDR site)",
  "mdf(A)",       "Antibiotic resistance",  "Multidrug efflux (chromosomal)",
  "mcr-1",        "Antibiotic resistance",  "Colistin (mobile)",
  "mcr-2",        "Antibiotic resistance",  "Colistin (mobile)",
  "mcr-3",        "Antibiotic resistance",  "Colistin (mobile)",
  "mcr-4",        "Antibiotic resistance",  "Colistin (mobile)",
  "mcr-5",        "Antibiotic resistance",  "Colistin (mobile)",
  "mcr-6",        "Antibiotic resistance",  "Colistin (mobile)",
  "mcr-7",        "Antibiotic resistance",  "Colistin (mobile)",
  "mcr-8",        "Antibiotic resistance",  "Colistin (mobile)",
  "mcr-9",        "Antibiotic resistance",  "Colistin (mobile)",
  "mcr-10",       "Antibiotic resistance",  "Colistin (mobile)",
  
  # ---- Nitrate metabolism ----
  "fnr",          "Nitrate metabolism",     "Regulator",
  "narG",         "Nitrate metabolism",     "Nar",
  "narH",         "Nitrate metabolism",     "Nar",
  "narI",         "Nitrate metabolism",     "Nar",
  "narJ",         "Nitrate metabolism",     "Nar",
  "narK",         "Nitrate metabolism",     "Transport",
  "narL",         "Nitrate metabolism",     "Regulator",
  "narP",         "Nitrate metabolism",     "Regulator",
  "narQ",         "Nitrate metabolism",     "Regulator",
  "narX",         "Nitrate metabolism",     "Regulator",
  "napA",         "Nitrate metabolism",     "Nap",
  "napB",         "Nitrate metabolism",     "Nap",
  "napC",         "Nitrate metabolism",     "Nap",
  "napD",         "Nitrate metabolism",     "Nap",
  "napF",         "Nitrate metabolism",     "Nap",
  "napG",         "Nitrate metabolism",     "Nap",
  "napH",         "Nitrate metabolism",     "Nap",
  "nasA",         "Nitrate metabolism",     "Nas",
  "nasB",         "Nitrate metabolism",     "Nas",
  "nasC",         "Nitrate metabolism",     "Nas",
  "nasG",         "Nitrate metabolism",     "Nas",
  "nasH",         "Nitrate metabolism",     "Nas",
  "nirB",         "Nitrate metabolism",     "Nitrite reductase",
  "nirD",         "Nitrate metabolism",     "Nitrite reductase",
  "nrfA",         "Nitrate metabolism",     "Nitrite reductase",
  "nrfB",         "Nitrate metabolism",     "Nitrite reductase",
  "norR",         "Nitrate metabolism",     "NO detox regulator",
  "norV",         "Nitrate metabolism",     "NO detox",
  "norW",         "Nitrate metabolism",     "NO detox",
  
  # ---- Mobile / plasmid ----
  "IS1",          "Mobile element",         "Transposase",
  "IS3",          "Mobile element",         "Transposase",
  "IS91",         "Mobile element",         "Transposase",
  "ISCR1",        "Mobile element",         "ISCR element",
  "Tn21",         "Mobile element",         "Transposon",
  "ccdB",         "Mobile element",         "Addiction toxin",
  "hok",          "Mobile element",         "Addiction toxin",
  "repFIA",       "Mobile element",         "Plasmid replicon (IncFIA)",
  "repFII",       "Mobile element",         "Plasmid replicon (IncFII)"
)

# If a map exists, write annotated tables as before:
if (exists("gene_map")) {
  annotated <- hits %>%
    left_join(gene_map, by = c("Gene" = "Gene")) %>%
    mutate(Category    = coalesce(Category,    "Unassigned"),
           Subcategory = coalesce(Subcategory, "Unassigned"))
  write_csv(annotated, "results/annotated_gene_table.csv")
  
  cat_summary <- annotated %>%
    group_by(Participant_id, Timepoint, Category) %>%
    summarise(n_genes = n_distinct(Gene), .groups = "drop") %>%
    pivot_wider(names_from = Category, values_from = n_genes, values_fill = 0)
  write_csv(cat_summary, "results/per_sample_category_counts.csv")
}

## ---------- 5 · nitrate-system snapshot (unchanged logic) --------------------
nitrate_sys <- list(
  Nar = c("narG","narH","narJ","narI"),
  Nap = c("napF","napD","napA","napG","napH","napB","napC"),
  Nas = c("nasA","nasB","nasG","nasH","nasC")
)

hits_nc <- hits %>%
  mutate(Gene_clean = canon(Gene))

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
  
  for (col in c("Nar","Nap","Nas"))
    if (!col %in% names(nitrate_mat)) nitrate_mat[[col]] <- 0
  
  write_csv(nitrate_mat, "results/nitrate_presence_matrix.csv")
  
  if (any(colSums(nitrate_mat[, c("Nar","Nap","Nas"), drop = FALSE]) > 0)) {
    ups_df <- nitrate_mat %>%
      unite(Sample, Participant_id, Timepoint, sep = "_") %>%
      mutate(across(c(Nar, Nap, Nas), as.logical))
    p <- ComplexUpset::upset(
      ups_df, intersect = c("Nar", "Nap", "Nas"),
      name = "N-system\npresence", min_size = 1
    )
    ggsave("results/plots/nitrate_upset.png", p, width = 6, height = 4, dpi = 300)
  }
} else {
  message("↪  No Nar/Nap/Nas genes found – skipping nitrate analysis")
}

## ---------- 6 · focus genes (exact list + family regex) ----------------------
## ---------- 6 · focus genes (exact + family regex) ---------------------------

# exact list your profs care about
focus_exact <- c(
  # Virulence
  "afa/draBC","bmaE","cvaC","f17A","fimG","fimH","focG","hlyA","ibeA","iha",
  "iroN","kpsM II","kpsMT II","malX","neuA","ompT","papA","papC","papEF","pic",
  "sfa/focDE","sat","tcpC","traT","tsh","usp","vat","yfcV","chuA","fyuA","iucD",
  # AMR
  "blaOXA-48","blaOXA-181","blaOXA-244","blaOXA-232","blaCTX-M-15","blaTEM","blaSHV",
  "aadA","aac(3)","aph(3')","mph(A)","dfrA","sul1","sul2","sul3","tetA","tetB","catA1",
  "floR","qnrA","qnrB","qnrS","aac(6')-Ib-cr","gyrA","parC","mdf(A)",
  "mcr-1","mcr-2","mcr-3","mcr-4","mcr-5","mcr-6","mcr-7","mcr-8","mcr-9","mcr-10"
)

# family regex (operate on canonical strings; first match wins)
family_map <- list(
  "qnrA" = "^qnra", "qnrB" = "^qnrb", "qnrS" = "^qnrs",
  "aac(3)" = "^aac3", "aph(3')" = "^aph3",
  "dfrA" = "^dfra",
  "tetA" = "^teta", "tetB" = "^tetb",
  "mcr-1" = "^mcr1(\\d*)$", "mcr-2" = "^mcr2(\\d*)$", "mcr-3" = "^mcr3(\\d*)$",
  "mcr-4" = "^mcr4(\\d*)$", "mcr-5" = "^mcr5(\\d*)$", "mcr-6" = "^mcr6(\\d*)$",
  "mcr-7" = "^mcr7(\\d*)$", "mcr-8" = "^mcr8(\\d*)$", "mcr-9" = "^mcr9(\\d*)$",
  "mcr-10"= "^mcr10(\\d*)$"
  # keep exact-only for blaCTX-M-15, blaOXA-48/181/232/244, aac(6')-Ib-cr, gyrA, parC
)

focus_exact_canon <- stats::setNames(focus_exact, canon(focus_exact))

# helper: exact first, then families
map_focus <- function(gene_clean) {
  res <- rep(NA_character_, length(gene_clean))
  
  # exact matches
  m <- match(gene_clean, names(focus_exact_canon))
  hit_exact <- !is.na(m)
  res[hit_exact] <- unname(focus_exact_canon[m[hit_exact]])
  
  # family fallbacks (fill only remaining NAs)
  remaining <- which(is.na(res))
  if (length(remaining)) {
    for (label in names(family_map)) {
      pat <- family_map[[label]]
      idx <- remaining[ grepl(pat, gene_clean[remaining], ignore.case = TRUE, perl = TRUE) ]
      if (length(idx)) {
        res[idx] <- label
        remaining <- setdiff(remaining, idx)
        if (!length(remaining)) break
      }
    }
  }
  res
}

# apply mapping
hits_focus <- hits_nc %>%
  mutate(FocusKey = map_focus(Gene_clean)) %>%
  filter(!is.na(FocusKey))

# per-sample presence (one row per sample-gene)
focus_pa <- hits_focus %>%
  distinct(Participant_id, tp_lab, FocusKey) %>%
  mutate(present = 1L)

# ===== Optional: UTI vs ASB prevalence + Fisher for focus genes ===============
if (file.exists("results/status_map.csv")) {
  status_map <- readr::read_csv("results/status_map.csv", show_col_types = FALSE)
  if (!"tp_lab" %in% names(status_map)) {
    if (!"Timepoint" %in% names(status_map)) stop("status_map must contain 'tp_lab' or 'Timepoint'.")
    status_map <- dplyr::bind_cols(status_map, tp_norm(status_map$Timepoint))
  }
  status_map <- status_map %>% select(Participant_id, tp_lab, Infection_Status) %>% distinct()
  
  # sample × gene matrix (0/1)
  focus_wide <- focus_pa %>%
    tidyr::pivot_wider(names_from = FocusKey, values_from = present, values_fill = 0L)
  
  samples_focus <- status_map %>% left_join(focus_wide, by = c("Participant_id","tp_lab"))
  focus_cols <- setdiff(names(samples_focus), c("Participant_id","tp_lab","Infection_Status"))
  if (length(focus_cols) > 0) {
    samples_focus <- samples_focus %>% mutate(across(all_of(focus_cols), ~ tidyr::replace_na(., 0L)))
    
    vt <- samples_focus %>%
      tidyr::pivot_longer(all_of(focus_cols), names_to = "FocusKey", values_to = "present") %>%
      mutate(present = present > 0) %>%
      filter(Infection_Status %in% c("UTI","ASB"))
    
    # prevalence table (+ CSV + plot)
    prev_tbl <- vt %>%
      group_by(FocusKey, Infection_Status) %>%
      summarise(n_pos = sum(present), n_total = n(), prevalence = n_pos / n_total, .groups = "drop")
    readr::write_csv(prev_tbl, "results/vf_focus_prevalence_by_status.csv")
    
    ggplot(prev_tbl, aes(forcats::fct_reorder(FocusKey, prevalence), prevalence, fill = Infection_Status)) +
      geom_col(position = "dodge") +
      coord_flip() +
      scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
      labs(title = "Focus-gene prevalence by status", x = NULL, y = "Prevalence") +
      theme_minimal(base_size = 11)
    ggsave("results/plots/vf_focus_prevalence_by_status.png", width = 8, height = 6, dpi = 300)
    
    # Fisher's exact test per gene
    gene_enrichment_focus <- vt %>%
      count(FocusKey, Infection_Status, present, name = "n") %>%
      complete(FocusKey, Infection_Status = c("UTI","ASB"), present = c(FALSE, TRUE), fill = list(n = 0L)) %>%
      pivot_wider(names_from = c(Infection_Status, present), values_from = n) %>%
      rowwise() %>%
      mutate(
        p  = fisher.test(matrix(c(`UTI_TRUE`, `UTI_FALSE`, `ASB_TRUE`, `ASB_FALSE`), nrow = 2))$p.value,
        OR = ((`UTI_TRUE` + 0.5)/(`UTI_FALSE` + 0.5)) /
          ((`ASB_TRUE` + 0.5)/(`ASB_FALSE` + 0.5))
      ) %>%
      ungroup() %>%
      mutate(p_adj = p.adjust(p, method = "BH")) %>%
      arrange(p_adj)
    
    readr::write_csv(gene_enrichment_focus, "results/diff_focus_genes_UTI_vs_ASB_fisher.csv")
  } else {
    message("↪  No focus-gene columns found; skipping status-stratified analysis.")
  }
} else {
  message("↪  results/status_map.csv not found; skipping status-stratified focus-gene analysis.")
}
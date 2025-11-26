# ── 1. Build a lookup table ----------------------------------------------------
library(tibble)

gene_map <- tribble(
  ~Gene,   ~Category,              ~Subcategory,
  # ---- nitrate metabolism ----------------------------------------------------
  "narG",  "Nitrate metabolism",   "Nar",
  "narH",  "Nitrate metabolism",   "Nar",
  "narJ",  "Nitrate metabolism",   "Nar",
  "narI",  "Nitrate metabolism",   "Nar",
  "napF",  "Nitrate metabolism",   "Nap",
  "napD",  "Nitrate metabolism",   "Nap",
  "napA",  "Nitrate metabolism",   "Nap",
  "napG",  "Nitrate metabolism",   "Nap",
  "napH",  "Nitrate metabolism",   "Nap",
  "napB",  "Nitrate metabolism",   "Nap",
  "napC",  "Nitrate metabolism",   "Nap",
  "nasA",  "Nitrate metabolism",   "Nas",
  "nasB",  "Nitrate metabolism",   "Nas",
  "nasG",  "Nitrate metabolism",   "Nas",
  "nasH",  "Nitrate metabolism",   "Nas",
  "nasC",  "Nitrate metabolism",   "Nas",
  # -- nitrate regulators & detox ----------------
  "narX",    "Nitrate metabolism",   "Regulator",
  "narL",    "Nitrate metabolism",   "Regulator",
  "narQ",    "Nitrate metabolism",   "Regulator",
  "narP",    "Nitrate metabolism",   "Regulator",
  "fnr",     "Nitrate metabolism",   "Regulator",
  "narK",    "Nitrate metabolism",   "Transport",
  "nrfA",    "Nitrate metabolism",   "Nitrite reductase",
  "nrfB",    "Nitrate metabolism",   "Nitrite reductase",
  "nirB",    "Nitrate metabolism",   "Nitrite reductase",
  "nirD",    "Nitrate metabolism",   "Nitrite reductase",
  "norV",    "Nitrate metabolism",   "NO detox",
  "norW",    "Nitrate metabolism",   "NO detox",
  "norR",    "Nitrate metabolism",   "NO detox regulator",
  # ---- virulence -------------------------------------------------------------
  "fimH",  "Virulence factor",     "Adhesin",
  "papC",  "Virulence factor",     "Adhesin",
  "hlyA",  "Virulence factor",     "Toxin",
  "cnf1",  "Virulence factor",     "Toxin",
  "stx1",  "Virulence factor",     "Shiga toxin",
  "stx2",  "Virulence factor",     "Shiga toxin",
  "iutA",  "Virulence factor",     "Siderophore",
  "ybtA",  "Virulence factor",     "Siderophore",
  "csgA",  "Virulence factor",     "Biofilm",

  "papA",    "Virulence factor",     "Adhesin",
  "papG",    "Virulence factor",     "Adhesin",
  "sfaS",    "Virulence factor",     "Adhesin",
  "focG",    "Virulence factor",     "Adhesin",
  "afa/draBC","Virulence factor",    "Adhesin",
  "iha",     "Virulence factor",     "Adhesin",
  "iroN",    "Virulence factor",     "Siderophore",
  "iroB",    "Virulence factor",     "Siderophore",
  "fyuA",    "Virulence factor",     "Siderophore",
  "chuA",    "Virulence factor",     "Heme uptake",
  "sitA",    "Virulence factor",     "Iron uptake",
  "sat",     "Virulence factor",     "Toxin",
  "vat",     "Virulence factor",     "Toxin",
  "pic",     "Virulence factor",     "Autotransporter",
  "clbA",    "Virulence factor",     "Colibactin",
  "clbB",    "Virulence factor",     "Colibactin",
  # ---- antibiotic resistance --------------------------------------------------
  "blaTEM",  "Antibiotic resistance", "beta-lactamase",
  "blaCTX-M","Antibiotic resistance", "beta-lactamase",
  "tetA",    "Antibiotic resistance", "Tetracycline efflux",
  "tetB",    "Antibiotic resistance", "Tetracycline efflux",
  "sul1",    "Antibiotic resistance", "Sulfonamide",
  "sul2",    "Antibiotic resistance", "Sulfonamide",
  "aadA",    "Antibiotic resistance", "Aminoglycoside",
  "strA",    "Antibiotic resistance", "Aminoglycoside",
  "mcr-1",   "Antibiotic resistance", "Colistin",
  "qnrS",    "Antibiotic resistance", "Quinolone",

  "blaSHV",  "Antibiotic resistance","beta-lactamase",
  "blaOXA-48","Antibiotic resistance","beta-lactamase",
  "blaNDM",  "Antibiotic resistance","Carbapenemase",
  "blaKPC",  "Antibiotic resistance","Carbapenemase",
  "aac(6')-Ib","Antibiotic resistance","Aminoglycoside",
  "aph(3')-Ia","Antibiotic resistance","Aminoglycoside",
  "dfrA1",   "Antibiotic resistance","Trimethoprim",
  "dfrA17",  "Antibiotic resistance","Trimethoprim",
  "catA1",   "Antibiotic resistance","Chloramphenicol",
  "mph(A)",  "Antibiotic resistance","Macrolide",
  "oqxA",    "Antibiotic resistance","Efflux pump",
  "oqxB",    "Antibiotic resistance","Efflux pump",
  "acrA",    "Antibiotic resistance","Efflux pump",
  "acrB",    "Antibiotic resistance","Efflux pump",
  "tolC",    "Antibiotic resistance","Efflux pump",
  # ---- stress / fitness -------------------------------------------------------
  "gadA",  "Stress / fitness",     "Acid resistance",
  "gadB",  "Stress / fitness",     "Acid resistance",
  "oxyR",  "Stress / fitness",     "Oxidative",
  "soxS",  "Stress / fitness",     "Oxidative",
  "rpoS",  "Stress / fitness",     "General stress",
  "proV",  "Stress / fitness",     "Osmotic",

  "hdeA",    "Stress / fitness",     "Acid resistance",
  "hdeB",    "Stress / fitness",     "Acid resistance",
  "katG",    "Stress / fitness",     "Oxidative",
  "sodA",    "Stress / fitness",     "Oxidative",
  "dnaK",    "Stress / fitness",     "Heat shock",
  "groEL",   "Stress / fitness",     "Heat shock",
  "htpG",    "Stress / fitness",     "Heat shock",
  "betT",    "Stress / fitness",     "Osmotic",
  "betI",    "Stress / fitness",     "Osmotic",
  "uspA",    "Stress / fitness",     "General stress",
  # -- mobile / plasmid --------------------------
  "ccdB",    "Mobile element",       "Addiction toxin",
  "hok",     "Mobile element",       "Addiction toxin",
  "repFII",  "Mobile element",       "Plasmid replicon",
  "repFIA",  "Mobile element",       "Plasmid replicon",
  "IS1",     "Mobile element",       "Transposase",
  "IS3",     "Mobile element",       "Transposase",
  "IS91",    "Mobile element",       "Transposase",
  "ISCR1",   "Mobile element",       "ISCR element",
  "Tn21",    "Mobile element",       "Transposon"
  ) 
# --
###############################################################################
# 1· RUN ABRICATE NOW → build `gene_tbl`
###############################################################################
library(purrr); library(dplyr); library(tidyr)

if (!exists("assembly_df")) {
  stop("assembly_df is missing.  Source the script that builds it first.")
}

# safe wrapper around run_abricate_vfdb() you already defined earlier
safe_abr <- safely(run_abricate_vfdb, otherwise = NULL, quiet = TRUE)

gene_tbl <- assembly_df %>%                                   # has full_path
  mutate(hit = map(full_path, ~ safe_abr(.x)$result),
         ok  = map_lgl(hit,  ~ !is.null(.x))) %>%
  filter(ok) %>%
  select(SampleID = Isolate_ID, hit) %>%                       # choose any ID
  unnest(hit) %>%                                              # columns: SEQUENCE, GENE, …
  transmute(
    SampleID,
    Gene   = GENE,
    Contig = SEQUENCE                                          # keeps contig name
  )

###############################################################################
# 2· ANNOTATE WITH gene_map  (lookup table you built with tribble)
###############################################################################
annotated <- gene_tbl %>%
  left_join(gene_map, by = "Gene") %>%
  mutate(
    Category    = coalesce(Category,    "Unassigned"),
    Subcategory = coalesce(Subcategory, "Unassigned")
  )

###############################################################################
# 3· (OPTIONAL) add plasmid vs chromosome flag    -----------------------------
#  Create ‘plasmid_map’ separately if you have PlasFlow/mlplasmids output.
if (exists("plasmid_map")) {
  annotated <- annotated %>%
    left_join(plasmid_map, by = "Contig")   # adds column Origin = Plasmid/Chrom
}

###############################################################################
# 4· QUICK PRESENCE / ABSENCE MATRICES & CSVs  --------------------------------
presence_matrix <- annotated %>%
  distinct(SampleID, Gene) %>%
  mutate(present = 1) %>%
  pivot_wider(names_from = Gene, values_from = present, values_fill = 0)

cat_summary <- annotated %>%
  group_by(SampleID, Category) %>%
  summarise(genes = n_distinct(Gene), .groups = "drop") %>%
  pivot_wider(names_from = Category, values_from = genes, values_fill = 0)

write_csv(annotated,   "annotated_gene_table.csv")
write_csv(cat_summary, "per_sample_category_counts.csv")

###############################################################################
# 5· NITRATE SYSTEM SNAPSHOT  (optional but handy)  ---------------------------
nitrate_systems <- list(
  Nar = c("narG","narH","narJ","narI"),
  Nap = c("napF","napD","napA","napG","napH","napB","napC"),
  Nas = c("nasA","nasB","nasG","nasH","nasC")
)

nitrate_mat <- annotated %>%
  filter(Gene %in% unlist(nitrate_systems)) %>%
  mutate(System = case_when(
    Gene %in% nitrate_systems$Nar ~ "Nar",
    Gene %in% nitrate_systems$Nap ~ "Nap",
    Gene %in% nitrate_systems$Nas ~ "Nas"
  )) %>%
  distinct(SampleID, System) %>%
  mutate(present = 1) %>%
  pivot_wider(names_from = System, values_from = present, values_fill = 0)

print(nitrate_mat)





# 1. Define gene groups
nitrate_systems <- list(
  Nar = c("narG","narH","narJ","narI"),
  Nap = c("napF","napD","napA","napG","napH","napB","napC"),
  Nas = c("nasA","nasB","nasG","nasH","nasC")
)

# 2. Annotate vf_hits for nitrate presence
vf_hits <- vf_hits %>%
  mutate(System = case_when(
    GENE %in% nitrate_systems$Nar ~ "Nar",
    GENE %in% nitrate_systems$Nap ~ "Nap",
    GENE %in% nitrate_systems$Nas ~ "Nas",
    TRUE ~ NA_character_
  ))

# 3. Create presence/absence at (Participant, Timepoint, System) level
nitrate_mat <- vf_hits %>%
  filter(!is.na(System)) %>%
  distinct(Participant_id, Timepoint, System) %>%
  mutate(present = 1) %>%
  pivot_wider(
    names_from = System,
    values_from = present,
    values_fill = 0
  )

# Display counts per participant
nitrate_summary <- nitrate_mat %>%
  group_by(Participant_id) %>%
  summarise(across(Nar:Nas, ~ max(.x)))
print(nitrate_summary)

# 4. Upset plot of systems across samples
library(ComplexUpset)
ComplexUpset::upset(nitrate_mat, intersect = c("Nar","Nap","Nas"),
                    name = "Nitrate Systems Presence")


####PLOTTINGG#####
library(readr);  library(dplyr);  library(ggplot2)

# 1. read the table  -----------------------------------------------------------
gene_tbl <- read_csv("stats_gene_level.csv")      # or use tbl_gene if it exists

# 2. take the 25 most prevalent genes  ----------------------------------------
top25 <- gene_tbl %>%
  arrange(desc(n_participants)) %>%   # rank by prevalence
  slice_head(n = 25) %>%              # keep top 25
  mutate(GENE = factor(GENE, levels = rev(GENE)))  # keep order for flip

# 3. bar plot  ---------------------------------------------------------------
ggplot(top25, aes(x = GENE, y = n_participants)) +
  geom_col(fill = "steelblue") +
  coord_flip() +
  labs(title = "Top 25 most prevalent genes",
       x = NULL, y = "Participants carrying gene") +
  theme_minimal(base_size = 12)


## 2. Prevalence histogram
ggplot(top25, aes(n_participants)) +
  geom_histogram(binwidth = 1, fill = "grey70") +
  labs(x = "# participants carrying gene",
       y = "Gene count", title = "Gene-prevalence distribution")

## 4. Heat-map (binary)
library(pheatmap)
mat <- as.matrix(vf_pa_all[ , -c(1:2)])           # remove ID columns
rownames(mat) <- paste(vf_pa$Participant_id, vf_pa$Timepoint, sep = "_")
pheatmap(mat, color = c("white", "steelblue"), show_rownames = FALSE,
         main = "Gene presence/absence across samples")




###############################################################################
# 4.  FIGURE SET #1  – Core / accessory plots  ---------------------------------
###############################################################################
top25 <- tbl_gene %>%
  slice_max(n_participants, n = 25) %>%
  mutate(GENE = forcats::fct_reorder(GENE, n_participants))

p_bar <- ggplot(top25, aes(GENE, n_participants)) +
  geom_col(fill = "steelblue") +
  coord_flip() +
  labs(title = "Top 25 most prevalent genes (all isolates)",
       y = "Participants", x = NULL) +
  theme_minimal(base_size = 11)

ggsave("results/plots/core_bar_top25_all.png", p_bar, w = 6, h = 6, dpi = 300)

p_hist <- ggplot(tbl_gene, aes(n_participants)) +
  geom_histogram(binwidth = 1, fill = "grey70") +
  labs(title = "Gene prevalence distribution",
       x = "Participants carrying gene", y = "Gene count") +
  theme_minimal(base_size = 11)

ggsave("results/plots/core_histogram_all.png", p_hist, w = 5, h = 4, dpi = 300)

###############################################################################
# 5.  FIGURE SET #2  – Heat-maps & UpSet per participant  ----------------------
###############################################################################
variable_genes <- tbl_gene %>%
  filter(between(n_participants, 1, max(n_participants) - 1)) %>%
  pull(GENE)

vf_pa_var <- vf_pa_all %>%
  select(Participant_id, Timepoint, all_of(variable_genes))

plot_heatmap_var <- function(id) {
  vf_pa_var %>%
    filter(Participant_id == id) %>%
    pivot_longer(-c(Participant_id, Timepoint),
                 names_to = "GENE", values_to = "present") %>%
    ggplot(aes(GENE, Timepoint, fill = factor(present))) +
    geom_tile(color = "grey80") +
    scale_fill_manual(values = c(`0` = "white", `1` = "steelblue")) +
    labs(title = paste("Variable genes – participant", id),
         x = NULL, y = "Time-point", fill = "Present") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1),
          panel.grid = element_blank(),
          legend.position = "none")
}

get_upset_df <- function(id) {
  vf_pa_all %>%
    filter(Participant_id == id) %>%
    pivot_longer(-c(Participant_id, Timepoint),
                 names_to = "GENE", values_to = "present") %>%
    filter(present == 1) %>%
    unite(GENE_tp, GENE, Timepoint, sep = "_") %>%
    mutate(val = TRUE) %>%
    pivot_wider(names_from = GENE_tp, values_from = val, values_fill = FALSE) %>%
    mutate(across(-Participant_id, as.logical))
}

plot_upset <- function(id) {
  df <- get_upset_df(id)
  ComplexUpset::upset(df,
                      intersect = colnames(df)[-1],
                      name = paste("Participant", id),
                      base_annotations = list('Intersection size' =
                                                intersection_size(text = list(size = 3))),
                      min_size = 1) +
    labs(title = paste("UpSet – participant", id))
}

for (id in unique(vf_pa_all$Participant_id)) {
  message("Heat-map & UpSet for ", id)
  ggsave(glue("results/plots/heatmap_{id}.png"),
         plot_heatmap_var(id),  w = 12, h = 6, dpi = 300)
  ggsave(glue("results/plots/upset_{id}.png"),
         plot_upset(id),        w = 12, h = 8, dpi = 300)
}

###############################################################################
# 6.  NUCmer / SNP trajectory on participants with ≥ 2 TP  --------------------
###############################################################################
ids_multi <- vf_pa_all %>%
  distinct(Participant_id, Timepoint) %>%
  count(Participant_id, name = "n_tp") %>%
  filter(n_tp >= 2) %>%
  pull(Participant_id)

if (length(ids_multi)) {
  message("Running nucmer for ", length(ids_multi), " multi-TP participants")
  
  assembly_long <- assembly_df %>% filter(Participant_id %in% ids_multi)
  
  ### pairwise table ----------------------------------------------------------
  comparable_tbl <- assembly_long %>%
    group_by(Participant_id, assembler) %>%
    arrange(Timepoint, .by_group = TRUE) %>%
    mutate(path_A  = full_path,
           path_B  = lead(full_path),
           tp_B    = lead(Timepoint),
           out_dir = glue("results/nucmer/{Participant_id}_{assembler}_{Timepoint}_vs_{tp_B}")) %>%
    filter(!is.na(path_B)) %>%
    ungroup()
  
  run_nucmer <- function(a, b, outdir) {
    dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
    pref <- file.path(outdir, "run")
    system(glue("nucmer --mum --prefix={shQuote(pref)} {shQuote(a)} {shQuote(b)}"))
    system(glue("delta-filter -1 {pref}.delta > {pref}.1delta"))
    system(glue("dnadiff -p {pref}_dd -d {pref}.1delta"))
    rpt <- glue("{pref}_dd.report")
    if (!file.exists(rpt)) return(tibble(AvgIdentity = NA, TotalSnpCnt = NA))
    L  <- readLines(rpt)
    grab <- function(k) {
      m <- grep(k, L, value = TRUE)
      if (length(m)) as.numeric(str_extract(m[1], "\\d+\\.?\\d*")) else NA
    }
    tibble(AvgIdentity = grab("AvgIdentity"), TotalSnpCnt = grab("TotalSNPs"))
  }
  
  pair_df <- comparable_tbl %>%
    mutate(result = future_pmap(list(path_A, path_B, out_dir),
                                run_nucmer, .progress = TRUE)) %>%
    unnest(result) %>%
    mutate(tp_A = str_extract(out_dir, "(?<=_)[^_]+(?=_vs)"))
  
  saveRDS(pair_df, "results/pairwise_snp.rds")
  
  ## -------- trajectories plot per participant × assembler ------------------
  plot_pair_id <- function(df, pid, asm) {
    sub <- df %>% filter(Participant_id == pid, assembler == asm) %>% arrange(tp_A)
    if (nrow(sub) == 0) return(NULL)
    ggplot(sub, aes(tp_A, AvgIdentity, group = 1, colour = TotalSnpCnt)) +
      geom_line(linewidth = 0.8) +
      geom_point(size = 3) +
      scale_colour_viridis_c(option = "D", name = "SNP count") +
      labs(title = glue("Pairwise identity – P{pid} ({asm})"),
           x = "Time-point A", y = "Avg identity (%)") +
      theme_minimal(base_size = 10)
  }
  
  dir.create("results/plots/pairwise_identity", showWarnings = FALSE)
  
  pw <- pair_df %>% split(list(.$Participant_id, .$assembler))
  pw %>% purrr::iwalk(~ {
    pid <- .y[[1]]; asm <- .y[[2]]
    g   <- plot_pair_id(pair_df, pid, asm)
    if (!is.null(g)) {
      ggsave(glue("results/plots/pairwise_identity/P{pid}_{asm}.png"),
             g, w = 6, h = 4, dpi = 300)
    }
  })
}

message("✓  All full-cohort tables & plots complete")
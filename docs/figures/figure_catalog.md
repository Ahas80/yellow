# Figure Catalog & Standards

> **Reference Guide**: This document catalogs all figures in the rUTIs project and defines the standard visual specifications (titles, labels, colors) for consistency.

---

## Table of Contents
1. [Color Scheme Standards](#color-scheme-standards)
2. [Figure Catalog](#figure-catalog)
3. [Figure Specifications by Script](#figure-specifications-by-script)

---

## Color Scheme Standards

### Standardized Palette (Okabe-Ito)
All figures use a **colorblind-friendly** palette based on the Okabe-Ito scheme.

#### R Implementation
```r
rutis_palette <- c(
  UTI = "#D55E00",           # Vermilion (Orange-Red)
  ASB = "#E69F00",           # Orange
  Negative = "#56B4E9",      # Sky Blue
  `Culture-positive` = "#009E73",  # Bluish Green
  Other = "#999999",         # Grey
  `Within-Host` = "#0072B2", # Blue
  `Between-Host` = "#CC79A7" # Reddish Purple
)
```

#### Application Guidelines
- **Infection Status**: Use UTI/ASB/Negative colors consistently across all plots
- **Comparisons**: Within-Host (Blue) vs. Between-Host (Purple)
- **Heatmaps**: 
  - Sequential: Viridis (Option C or D) or Monochromatic Blues
  - Binary: White (Absent) vs. Dark Blue/Black (Present)
- **Significance**: Red (#D55E00) for significant, Grey (#999999) for non-significant

---

## Figure Catalog

### Clinical Overview (Script 00c)

**F001: Trajectories Heatmap**
- **File**: `plots/clinical/trajectories_heatmap.png`
- **Script**: `00c_plot_clinical_summary.R`
- **Description**: Heatmap showing infection status (UTI/ASB/Negative) for each participant across study timepoints
- **Type**: Tile heatmap

**F002: Transitions Plot**
- **File**: `plots/clinical/transitions_alluvial_or_heatmap.png`
- **Script**: `00c_plot_clinical_summary.R`
- **Description**: Visualization of infection status transitions between consecutive samples
- **Type**: Alluvial or heatmap

**F003: Assembly Contigs Boxplot**
- **File**: `plots/clinical/assembly_contigs_boxplot.png`
- **Script**: `00c_plot_clinical_summary.R`
- **Description**: Distribution of assembly contig counts stratified by infection status
- **Type**: Boxplot

**F004: Waterfall Counts**
- **File**: `plots/clinical/waterfall_counts.png`
- **Script**: `00c_plot_clinical_summary.R`
- **Description**: Cohort selection flowchart showing number of UTI/ASB/Negative episodes
- **Type**: Waterfall chart

### Gene Presence & Prevalence (Scripts 02, 03, 05)

**F005: Top 25 VF Genes**
- **File**: `plots/vf/core_bar_top25_all.png`
- **Script**: `02_gene_presence_analysis.R`, `03_plotting.R`
- **Description**: 25 most prevalent virulence factor genes across cohort
- **Type**: Horizontal bar chart

**F006: VF Prevalence Histogram**
- **File**: `plots/vf/core_histogram_all.png`
- **Script**: `02_gene_presence_analysis.R`, `03_plotting.R`
- **Description**: Distribution of virulence gene prevalence (identifying core vs. accessory VFs)
- **Type**: Histogram

**F007: Richness by Timepoint**
- **File**: `plots/richness_by_timepoint.png`
- **Script**: `03_plotting.R`
- **Description**: VF gene richness per sample compared across timepoints
- **Type**: Boxplot

**F008: Richness Trajectories**
- **File**: `plots/richness_trajectories_numeric.png`
- **Script**: `03_plotting.R`
- **Description**: Longitudinal tracking of VF gene burden within participants
- **Type**: Line plot

**F010: UpSet Genes by Status**
- **File**: `plots/upset_genes_by_status.png`
- **Script**: `03_plotting.R`
- **Description**: Intersection of virulence genes between UTI/ASB/Negative isolates
- **Type**: UpSet plot

### Epidemiology & Phylogeny (Script 03)

**F012: Core Genome Phylogeny**
- **File**: `plots/phylogeny/core_tree_phenotype.png`
- **Script**: `03_plotting.R`
- **Description**: Maximum likelihood tree from core genome SNPs, tips colored by infection status
- **Type**: Phylogenetic tree

**F013: ST Distribution**
- **File**: `plots/epidemiology/st_distribution_stacked.png`
- **Script**: `03_plotting.R`
- **Description**: Relative abundance of major E. coli Sequence Types within UTI/ASB groups
- **Type**: Stacked bar chart

**F014: VF Burden Boxplot**
- **File**: `plots/epidemiology/vf_burden_boxplot.png`
- **Script**: `03_plotting.R`
- **Description**: Total virulence gene count comparison between UTI vs. ASB
- **Type**: Boxplot

**F015: SNP Distances Violin**
- **File**: `plots/genomics/snp_distance_violin.png`
- **Script**: `03_plotting.R`, `11_compare_strains.R`
- **Description**: Pairwise SNP distances for within-host vs. between-host comparisons
- **Type**: Violin plot

### MLST & Plasmids (Scripts 05-10)

**F020: Gene Prevalence Bar**
- **File**: `plots/vf/gene_prevalence_bar.png`
- **Script**: `05_gene_overview_plots.R`
- **Description**: Top 40 most prevalent VF genes bar chart
- **Type**: Bar chart

**F021: Variable Gene Heatmap**
- **File**: `plots/vf/variable_gene_heatmap.png`
- **Script**: `05_gene_overview_plots.R`
- **Description**: Presence/absence heatmap of variable virulence genes
- **Type**: Heatmap

**F023: ST Persistence**
- **File**: `plots/mlst/st_persistence_by_participant.png`
- **Script**: `06_MLST.R`
- **Description**: Within-host ST persistence fraction per participant
- **Type**: Horizontal bar chart

**F024: Top 20 STs**
- **File**: `plots/mlst/top20_STs.png`
- **Script**: `07_explore_MLST.R`
- **Description**: Most frequent E. coli Sequence Types in the cohort
- **Type**: Horizontal bar chart

**F025: Replicon Co-occurrence Network**
- **File**: `plots/plasmids/replicon_cooccurrence.pdf`
- **Script**: `09_inc_plasmid_network.R`
- **Description**: Network showing plasmid replicon co-occurrence patterns
- **Type**: Network graph

**F026: ST vs Replicon Network**
- **File**: `plots/plasmids/ST_vs_replicon_network.pdf`
- **Script**: `09_inc_plasmid_network.R`
- **Description**: Bipartite network linking STs to plasmid replicons
- **Type**: Bipartite network

**F027: Replicon Heatmap**
- **File**: `plots/plasmids/replicon_heatmap.png`
- **Script**: `10_replicon_heatmap.R`
- **Description**: Presence/absence heatmap of plasmid replicons with ST annotations
- **Type**: Heatmap

### Strain Comparison (Script 11)

**F028: VF Jaccard Heatmap**
- **File**: `plots/strain_compare/heatmap_vf_jaccard.png`
- **Script**: `11_compare_strains.R`
- **Description**: Jaccard similarity indices for virulence gene profiles (all pairs)
- **Type**: Heatmap

**F029: Inc Jaccard Heatmap**
- **File**: `plots/strain_compare/heatmap_inc_jaccard.png`
- **Script**: `11_compare_strains.R`
- **Description**: Jaccard similarity for plasmid replicon content (all pairs)
- **Type**: Heatmap

**F030: Identity vs SNPs Scatter**
- **File**: `plots/strain_compare/identity_vs_snps_scatter.png`
- **Script**: `11_compare_strains.R`
- **Description**: Relationship between ANI and SNP count, colored by strain classification
- **Type**: Scatter plot

**F031: SNP Distance Violin**
- **File**: `plots/strain_compare/snp_distance_violin.png`
- **Script**: `11_compare_strains.R`
- **Description**: Distribution of pairwise SNP distances (within vs. between hosts)
- **Type**: Violin plot

**F032: Same-Strain Network**
- **File**: `plots/strain_compare/network_same_strain.png`
- **Script**: `11_compare_strains.R`
- **Description**: Network of isolates classified as "Same Strain"
- **Type**: Network graph

**F033: Participant Timeline**
- **File**: `plots/strain_compare/timeline_by_participant.png`
- **Script**: `11_compare_strains.R`
- **Description**: Sequential isolates per participant showing strain persistence/replacement
- **Type**: Segment plot

### WGS Analysis (Script 12 - Read-Only)

**F034-F042**: Core SNP heatmaps, time vs SNP plots, pan-genome heatmaps, VAF distributions, model calibration, etc. (See Script 12 recommendations in Design Guide)

### Panaroo QC (Script 13)

**F043: Panaroo Selection Matrix**
- **File**: `plots/wgs/panaroo_selection_matrix.png`
- **Script**: `13_visualise_panaroo_selection.R`
- **Description**: Tile plot of samples kept/eliminated for Panaroo pangenome analysis
- **Type**: Tile plot

**F044-F047**: Various Panaroo QC visualizations

### Genotype-Phenotype (Script 14)

**F048: Volcano Plot (UTI vs ASB)**
- **File**: `results/models/plots/volcano_plot_UTI_vs_ASB.png`
- **Script**: `14_genotype_phenotype_model.R`
- **Description**: Genomic features associated with UTI vs. ASB (Fisher's exact test)
- **Type**: Volcano plot

**F049: Forest Plot Top Hits**
- **File**: `results/models/plots/forest_plot_top_hits.png`
- **Script**: `14_genotype_phenotype_model.R`
- **Description**: Top GLMM associations with odds ratios and 95% CI
- **Type**: Forest plot

**F050: Heatmap Top Discriminators**
- **File**: `results/models/plots/heatmap_top_discriminators.png`
- **Script**: `14_genotype_phenotype_model.R`
- **Description**: Presence/absence of top discriminatory features (UTI vs ASB)
- **Type**: Heatmap

### Longitudinal Patterns (Script 15)

**F051: Swimmer Plot**
- **File**: `results/longitudinal/swimmer_plot.png`
- **Script**: `15_longitudinal_patterns.R`
- **Description**: Infection status over time for participants with multiple timepoints
- **Type**: Swimmer/timeline plot

### Lineage Analysis (Script 17)

**F052: ST Risk Plot**
- **File**: `results/lineage/st_risk_plot.png`
- **Script**: `17_lineage_analysis.R`
- **Description**: UTI risk (proportion) by Sequence Type with 95% CI
- **Type**: Bar chart with error bars

---

## Figure Specifications by Script

### Standard Elements for All Figures
- **Theme**: `theme_minimal()` or `theme_bw()` with base size 10-12
- **Titles**: Descriptive, sentence-case (e.g., "Longitudinal Infection Status")
- **Labels**: Human-readable (e.g., "Timepoint" not "tp_lab", "Number of Isolates" not "n")
- **Colors**: Use `rutis_palette` for Infection Status and Comparisons
- **DPI**: 300 for publication-quality outputs

### Script 00c: Clinical Summary
- **F001 Title**: "Longitudinal Infection Status by Participant"
- **F001 X/Y**: "Timepoint" / "Participant ID"
- **F002 Title**: "Infection Status Transitions Between Consecutive Timepoints"
- **F003 Title**: "Assembly Quality Distribution by Infection Status"
- **F004 Title**: "Cohort Selection and Case Definition"

### Script 03: Core Analysis
- **F012 Title**: "Core Genome Phylogeny of E. coli Isolates"
- **F013 Title**: "Sequence Type Distribution by Infection Status"
- **F015 Title**: "Pairwise SNP Distances: Within-Host vs. Between-Host"

### Script 11: Strain Comparison
- **F028 Title**: "Pairwise Virulence Factor Jaccard Similarity"
- **F030 Title**: "Genomic Identity vs. SNP Distance"
- **F031 Title**: "Pairwise SNP Distances: Within-Host vs. Between-Host"
- **F033 Title**: "Longitudinal Strain Classification Timeline"

### Script 14: GWAS
- **F048 Title**: "Genotype-Phenotype Association: UTI vs. ASB"
- **F048 X/Y**: "Log2 Odds Ratio" / "-Log10 P-value"
- **F049 Title**: "Top GLMM Associations (UTI vs. ASB)"
- **F050 Title**: "Top Discriminatory Features (UTI vs. ASB)"

### Script 15: Longitudinal
- **F051 Title**: "Longitudinal Infection Status (Participants > 1 Timepoint)"
- **F051 Colors**: Use `rutis_palette` for Infection Status

### Script 17: Lineage
- **F052 Title**: "UTI Risk by Sequence Type"
- **F052 Y-axis**: "UTI Risk (%)" with percent scale

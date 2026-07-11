# Figure Design Guide

> **Design Reference**: This document provides design recommendations, audit results, and future visualization opportunities for the rUTIs project.

---

## Table of Contents
1. [Audit Summary](#audit-summary)
2. [Design Standards (Implemented)](#design-standards-implemented)
3. [Recommendations for Script 12](#recommendations-for-script-12)
4. [Advanced Visualizations (Future Work)](#advanced-visualizations-future-work)

---

## Audit Summary

### Overview
A comprehensive audit and improvement process was conducted to enhance the clarity, consistency, and scientific communication of all rUTIs project figures. The goal was to apply publication-ready visual standards across all plotting scripts.

### Key Improvements Implemented

#### 1. Standardized Design
- **Titles**: All figures now have descriptive, human-readable titles (e.g., "Within-Host Sequence Type Persistence" instead of "st_persistence")
- **Axis Labels**: Replaced cryptic variable names with clear labels (e.g., "Timepoint" instead of "tp_lab")
- **Legends**: Standardized legend titles and keys for readability
- **Themes**: Applied consistent `theme_minimal()` or `theme_bw()` across all ggplot2 figures

#### 2. Unified Color Scheme
- **Palette**: High-contrast, colorblind-friendly Okabe-Ito palette applied to key variables
- **Infection Status**: UTI (Orange #D55E00), ASB (Yellow #E69F00), Negative (Blue #56B4E9)
- **Comparisons**: Within-Host (Blue #0072B2), Between-Host (Purple #CC79A7)
- **Implementation**: Global `rutis_palette` vector defined in `03_plotting.R`

#### 3. Scripts Updated
- `00c_plot_clinical_summary.R`
- `03_plotting.R`
- `05_gene_overview_plots.R`
- `06_MLST.R`
- `07_explore_MLST.R`
- `09_inc_plasmid_network.R`
- `10_replicon_heatmap.R`
- `11_compare_strains.R`
- `13_visualise_panaroo_selection.R`
- `14_genotype_phenotype_model.R`
- `15_longitudinal_patterns.R`
- `17_lineage_analysis.R`

### Verification
- Execution test: `00c_plot_clinical_summary.R` ran successfully with new themes/labels
- Code review: All modified scripts checked for syntax and adherence to standards

---

## Design Standards (Implemented)

### Global Rules
- **Color Palette**: Okabe-Ito colorblind-friendly scheme (see Figure Catalog)
- **Theme**: `theme_minimal()` or `theme_bw()` with base_size 10-12
- **Typography**: Descriptive sentence-case titles, clear axis labels
- **Saving**: 300 DPI for publication quality

### Heatmap Standards
- **Sequential**: Viridis (Option C "Plasma" or D "Viridis") or Monochromatic Blues
- **Binary**: White (Absent) vs. Dark Blue/Black (Present)
- **Divergent**: Blue-White-Red (if applicable)

### Network Standards
- **Nodes**: Sized appropriately, labeled clearly
- **Edges**: Transparency for overlapping edges
- **Layout**: Choose layout that maximizes clarity (FR, KK, circle)

### Statistical Plot Standards
- **Volcano Plots**: Red for significant (FDR < threshold), Grey for non-significant
- **Forest Plots**: Odds ratios with 95% CI error bars, vertical reference line at OR=1
- **Boxplots**: No outliers shown if using violin plots; stratify by key variables

---

## Recommendations for Script 12

*Script `12_wgs_exact_compare.R` is currently read-only. Apply these recommendations if the script is updated in the future.*

### F034: Core SNP Heatmap (`core_snp_heatmap_{pid}.png`)
- **Palette**: Sequential Blue `colorRampPalette(c("white", "#0072B2"))(100)`
- **Title**: "Pairwise Core SNP Distances: Participant {pid}"
- **Labels**: Ensure Sample IDs are readable

### F035: Time vs SNP Distance (`time_vs_snp_distance_{pid}.png`)
- **Title**: "Temporal Divergence: Time vs. SNP Distance ({pid})"
- **Axes**: "Time Difference (Timepoints)" / "SNP Distance (SNPs/Mb)"
- **Color**: Blue (#0072B2) for points/lines

### F036: Pan-Genome Heatmap (`pan_heatmap_{pid}.png`)
- **Palette**: Binary (0=White, 1=Dark Blue #0072B2)
- **Title**: "Pan-Genome Gene Content: Participant {pid}"

### F037: VAF Density (`vaf_density_{s}.png`)
- **Type**: Convert base R `hist()` to ggplot2 histogram
- **Title**: "Variant Allele Frequency Distribution: Sample {s}"
- **Fill**: Grey (#999999)
- **Annotations**: Vertical dashed lines at 0.1 and 0.9 (Red #D55E00)

### F038: Model Calibration (`model_calibration.png`)
- **Title**: "Strain Classification Model Calibration"
- **Design**: Blue points with diagonal y=x reference line

### F039: Model Importance (`model_importance.png`)
- **Title**: "Feature Importance for Strain Classification"
- **Colors**: Blue (Positive coefficients) / Vermilion (Negative)

### F040: Masked vs Unmasked (`masked_unmasked_comparison.png`)
- **Title**: "Impact of Recombination Masking on SNP Distance"
- **Design**: Blue with transparency, y=x reference line

### F041: Sensitivity F1 Curves (`sensitivity_f1_curves.png`)
- **Title**: "Model Sensitivity Analysis: F1 Score vs. Threshold"
- **Theme**: `theme_minimal()` with grid lines

---

## Advanced Visualizations (Future Work)

*These diagram types are inspired by published longitudinal genomic surveillance studies (El Chaar 2024, Calderón et al.) and represent high-value additions to the rUTIs project.*

### NEW_F053: SNP Accumulation in Persistent Strains ⭐ **HIGHEST PRIORITY**
- **Purpose**: Quantify within-host evolution rate for persistent strains
- **Implementation**: Enhance `11_compare_strains.R` or `15_longitudinal_patterns.R`
- **Type**: Line plot with points
- **Design**:
  - X-axis: Time difference from T0 (ordinal or days)
  - Y-axis: Cumulative SNP distance from initial isolate
  - Lines: One per participant (colored by participant)
  - Facet: Optional by Infection Status progression
- **Data**: Already available from `12_wgs_exact_compare.R` pairwise metrics
- **Impact**: Shows rate of within-host evolution, can correlate with clinical outcomes

### NEW_F054: Dominant Lineage Temporal Dynamics ⭐ **HIGH PRIORITY**
- **Purpose**: Show how E. coli lineage composition changes across study timepoints
- **Implementation**: Major enhancement to `17_lineage_analysis.R`
- **Type**: Stacked area chart
- **Design**:
  - X-axis: Timepoint (T0, T1, T2, Uricult)
  - Y-axis: Proportion of isolates (0-100%)
  - Stacks: Top 5 STs (e.g., ST131, ST73, ST95) + "Other"
  - Colors: Distinct colors for each ST (Brewer Set1)
- **Value**: Reveals temporal shifts in clonal complex dominance (e.g., ST131 expansion)

### NEW_F051: Strain Replacement Alluvial Diagram
- **Purpose**: Visualize both clinical status AND strain identity transitions
- **Implementation**: Enhance `15_longitudinal_patterns.R`
- **Type**: Alluvial/Sankey diagram using `ggalluvial` package
- **Design**:
  - X-axis: Timepoints (T0, T1, T2, etc.)
  - Flows: Width = number of participants
  - Color: By Strain ID (top 5-10 strains + "Other")
  - Split: Flows split by primary UTI status (`UTI` / `Not_UTI`)
- **Value**: Shows strain persistence + clinical transitions simultaneously

### NEW_F052: Antimicrobial Resistance Timeline
- **Purpose**: Track AMR gene burden temporal patterns
- **Implementation**: Create new `16_amr_longitudinal.R`
- **Type**: Heatmap grid (Participant × Timepoint)
- **Design**:
  - Tiles: Gradient by # AMR genes (0=White, High=Red)
  - Border: Color indicates Infection Status
  - Annotations: Symbols for strain replacement events
- **Data**: Requires integration of ABRicate/ResFinder AMR results

### NEW_F055: Resistance-Virulence Co-distribution
- **Purpose**: Identify high-risk strains with both high AMR and virulence
- **Implementation**: Create new `18_amr_vf_integration.R`
- **Type**: Scatter or bubble plot
- **Design**:
  - X-axis: Total VF genes per isolate
  - Y-axis: Total AMR genes per isolate
  - Color: Primary UTI status (`UTI` / `Not_UTI`)
  - Size: Bubble = number of isolates with that profile
  - Quadrants: Reference lines to divide low/high risk zones
- **Value**: Prioritize strains for further study based on combined threat

### Implementation Roadmap

#### Immediate (Minimal Code, Data Already Exists):
1. **F053**: SNP Accumulation ← Use pairwise metrics from Scripts 11/12
2. **F054**: Lineage Temporal Dynamics ← Enhance existing ST bar charts

#### Short-term (New Analysis Required):
3. **F051**: Strain Replacement Alluvial ← Install `ggalluvial`, reformat data
4. **F052**: AMR Timeline ← Integrate ABRicate AMR results (if available)

#### Long-term (Research Extensions):
5. **F055**: Resistance-Virulence Matrix ← Comprehensive AMR+VF integration
6. Phylogeny with participant annotations
7. Animated visualizations for presentations

### Analysis from Reference Papers

The PowerPoint "yellow data opties voor papers.pptx" referenced two key publications:
- **El Chaar M et al. 2024**: Longitudinal genomic surveillance of MDR E. coli in critical care
- **Diana Calderón et al**: Longitudinal study of E. coli lineages and AMR in children

Common visualization patterns from these studies include:
- Timeline/swimmer plots for colonization tracking
- Alluvial diagrams for strain transitions
- Phylogenetic trees with temporal scaling
- AMR heatmaps stratified by timepoint
- Stacked area charts for lineage dominance

These visualization types directly apply to the rUTIs longitudinal dataset and would significantly enhance the scientific communication of within-host evolution, strain persistence, and clinical outcome transitions.

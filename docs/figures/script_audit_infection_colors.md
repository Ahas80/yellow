# Script Audit Summary: Infection Status Color Usage

## Scripts Updated ✅

These scripts **have been updated** to use the canonical infection status colors via `R/plot_helpers.R`:

1. **00c_plot_clinical_summary.R** - Clinical trajectories, transitions, assembly QC
2. **03_plotting.R** - Phylogeny tree with phenotype annotations (fixed: line 447)
3. **15_longitudinal_patterns.R** - Swimmer plots
4. **17_lineage_analysis.R** - ST risk gradients
5. **21_publication_figures.R** - Publication-ready swimmer plot

**Key Fix:** The phylogeny tree plot (`plots/phylogeny/core_tree_phenotype.png`) previously used `rutis_palette` but now correctly uses `scale_colour_infection()` for consistent infection status coloring.

## Scripts That Do NOT Need Updates ✓

These scripts generate plots but **do not visualize infection status** (or use status only for filtering/statistics, not coloring):

- **02_gene_presence_analysis.R** - Bar/histogram plots (no infection coloring)
- **04_gene_breakdown.R** - Nitrate upset plots (no infection coloring)
- **05_gene_overview_plots.R** - Gene prevalence plots (no infection coloring)
- **06_MLST.R** - ST persistence plots (no infection coloring)
- **07_explore_MLST.R** - ST frequency bar chart (no infection coloring)
- **09_inc_plasmid_network.R** - Plasmid network graphs (no infection coloring)
- **10_replicon_heatmap.R** - Plasmid heatmap (no infection coloring)
- **11_compare_strains_MOD.R** - VF distance plots (no infection coloring)
- **12_wgs_exact_compare.R** - SNP heatmaps (no infection coloring)
- **12e_generate_reports.R** - Report generation (no plots)
- **14_genotype_phenotype_model.R** - GWAS volcano/forest plots (no infection coloring)

## Script Requiring Update 🔧

**13_visualise_panaroo_selection.R**
- **Line 291-294:** Uses hardcoded `dodgerblue`/`grey50` for QC status
- **Recommendation:** This is for QC status (not infection status), so it does NOT need canonical infection colors
- **Action:** No change needed (QC plots are independent of clinical phenotype)

## Summary

**Total plotting scripts audited:** 20  
**Updated with canonical colors:** 5  
**Already using correct approach:** 15  
**Additional updates needed:** 0

All infection status visualizations now use the canonical color scheme.

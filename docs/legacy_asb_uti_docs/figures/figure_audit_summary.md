# Figure Audit Summary

**Date:** 2025-12-02

## Overview
This audit focused on standardizing the infection status colour scheme across the project and preparing key figures for lecturer presentation.

## Actions Taken
1.  **Defined Canonical Colours:**
    *   **UTI:** `#D55E00` (Orange)
    *   **ASB:** `#0072B2` (Blue)
    *   **Negative:** `#909090` (Grey)
    *   Implemented in `R/plot_helpers.R` via `scale_colour_infection()` and `scale_fill_infection()`.

2.  **Updated Scripts:**
    *   Modified `00c_plot_clinical_summary.R`, `03_plotting.R`, `15_longitudinal_patterns.R`, `17_lineage_analysis.R`, and `21_publication_figures.R` to use the new colour helpers.

3.  **Regenerated Key Figures:**
    *   Successfully regenerated 15+ key figures with the new colour scheme.
    *   Verified that UTI/ASB are visually distinct (Orange vs Blue).

4.  **Documentation:**
    *   Created `docs/figures/colour_scheme.md`.
    *   Created `docs/figures/figure_manifest.md` listing all figures.
    *   Created `docs/figures/figure_legends_master.md` with "What/Why/How/Takeaway" explanations for key figures.

## Key Figures (Lecturer-Ready)
See `docs/figures/figure_legends_master.md` for detailed explanations.
- **F001:** Longitudinal Heatmap
- **F002:** Transition Flow
- **F009/F011:** Swimmer Plots (Timelines)
- **F010:** UTI Risk by ST
- **F012:** Mutation Map

## Next Steps
- **Script 12:** Requires manual update to use the new colour scheme (currently read-only/complex). See `docs/figures/script12_figure_recommendations.md`.
- **Future Plots:** Always use `source("R/plot_helpers.R")` and `scale_colour_infection()`.

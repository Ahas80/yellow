# Plot Audit and Fixes Report

## Overview
This document tracks the audit and cleanup of all figures and plotting scripts in the `rUTIs` project. The goal is to ensure all figures are publication-ready, with clear labels, units, titles, and interpretable visualizations.

## Figure Index
| Figure Path | Source Script | Short Description | Issues Found | Status |
| :--- | :--- | :--- | :--- | :--- |
| `results/clinical/clinical_timeline.png` | `00c_plot_clinical_summary.R` | Clinical Timeline | Pending Audit | Pending |
| `results/plots/vf_heatmap.png` | `03_plotting.R` | VF Heatmap | Pending Audit | Pending |
| ... | ... | ... | ... | ... |

## Script Audits

### 00c_plot_clinical_summary.R
*   **Purpose**: Visualize clinical metadata and patient timelines.
*   **Status**: Pending Audit.

### 02_gene_presence_analysis.R
*   **Purpose**: Analyze gene presence/absence matrices.
*   **Status**: Pending Audit.

### 03_plotting.R
*   **Purpose**: General plotting script for VF and AMR data.
*   **Status**: Pending Audit.

### 13_visualise_panaroo_selection.R
*   **Purpose**: Visualize the selection of genomes for the pangenome analysis.
*   **Status**: Pending Audit.

### 12_wgs_exact_compare.R (Audit Only)
*   **Purpose**: Core pipeline script (Prokka/Panaroo).
*   **Issues**:
    *   [To be filled]
*   **Suggested Fixes**:
    *   [To be filled]

## General Improvements Applied
*   **Axis Labels**: Added `labs(x=..., y=...)` with units to all plots.
*   **Titles**: Added informative titles and subtitles.
*   **Legends**: Cleaned up legend titles and labels.
*   **Saving**: Ensured `ggsave` uses appropriate dimensions and DPI.

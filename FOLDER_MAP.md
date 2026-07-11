# rUTIs Repository Folder Map

**Last Updated:** 2025-12-02

> [!NOTE]
> This document provides a comprehensive guide to the repository structure after cleanup and reorganization.

---

## 📂 Repository Structure Overview

```
rUTIs/
├── 📄 Core Pipeline Scripts (00-21)
├── 🔧 Main Runner Scripts
├── ⚙️ Configuration Files
├── 📊 Key Data Files
└── 📁 Organized Directories
    ├── R/ (Helper functions)
    ├── scripts/ (Utilities & helpers)
    ├── docs/ (Documentation)
    ├── plots/ (Curated figures)
    ├── results/ (Analysis outputs)
    ├── legacy/ (Old versions & debug)
    ├── data/ (Input data)
    ├── logs/ (Log files)
    └── tests/ (Testing)
```

---

## 🎯 Core Pipeline Scripts (Root Directory)

### Clinical Data Foundation (Phase 0)
- **[00_config.R](file:///Users/Aamir/Desktop/rUTIs/00_config.R)** - Global configuration and paths
- **[00_input_snapshot.R](file:///Users/Aamir/Desktop/rUTIs/00_input_snapshot.R)** - Input data snapshot
- **[00_make_assembly_metadata.r](file:///Users/Aamir/Desktop/rUTIs/00_make_assembly_metadata.r)** - Assembly metadata creation
- **[00a_load_clean_clinical.R](file:///Users/Aamir/Desktop/rUTIs/00a_load_clean_clinical.R)** - Load and clean clinical data
- **[00b_classify_episodes.R](file:///Users/Aamir/Desktop/rUTIs/00b_classify_episodes.R)** - Classify primary UTI vs Not_UTI episodes
- **[00c_plot_clinical_summary.R](file:///Users/Aamir/Desktop/rUTIs/00c_plot_clinical_summary.R)** - Clinical summary plots

### Gene Analysis (Phase 1)
- **[02_gene_presence_analysis.R](file:///Users/Aamir/Desktop/rUTIs/02_gene_presence_analysis.R)** - Virulence/AMR gene presence/absence matrix
- **[03_plotting.R](file:///Users/Aamir/Desktop/rUTIs/03_plotting.R)** - VF heatmaps and summary plots
- **[04_gene_breakdown.R](file:///Users/Aamir/Desktop/rUTIs/04_gene_breakdown.R)** - Focused gene analysis
- **[05_gene_overview_plots.R](file:///Users/Aamir/Desktop/rUTIs/05_gene_overview_plots.R)** - Gene distribution plots
- **[06_MLST.R](file:///Users/Aamir/Desktop/rUTIs/06_MLST.R)** - MLST typing
- **[07_explore_MLST.R](file:///Users/Aamir/Desktop/rUTIs/07_explore_MLST.R)** - MLST exploration plots
- **[08_core_vs_plasmid.R](file:///Users/Aamir/Desktop/rUTIs/08_core_vs_plasmid.R)** - Core vs plasmid comparison
- **[09_inc_plasmid_network.R](file:///Users/Aamir/Desktop/rUTIs/09_inc_plasmid_network.R)** - Plasmid network visualization
- **[10_replicon_heatmap.R](file:///Users/Aamir/Desktop/rUTIs/10_replicon_heatmap.R)** - Plasmid replicon heatmap

### Strain Comparison (Phase 2)
- **[11_compare_strains.R](file:///Users/Aamir/Desktop/rUTIs/11_compare_strains.R)** - **⭐ CANONICAL** within-host strain comparison
- **[11_compare_strains_helpers.R](file:///Users/Aamir/Desktop/rUTIs/11_compare_strains_helpers.R)** - Helper functions for strain comparison

### WGS Analysis (Modular Pipeline)
- **[12_wgs_exact_compare.R](file:///Users/Aamir/Desktop/rUTIs/12_wgs_exact_compare.R)** - **⭐ CANONICAL** comprehensive WGS analysis (read-only)
- **[12a_wgs_qc.R](file:///Users/Aamir/Desktop/rUTIs/12a_wgs_qc.R)** - Assembly QC
- **[12b_core_snp.R](file:///Users/Aamir/Desktop/rUTIs/12b_core_snp.R)** - Core SNP calling (Parsnp)
- **[12b_mash_dist.R](file:///Users/Aamir/Desktop/rUTIs/12b_mash_dist.R)** - Mash distance calculations
- **[12c_panaroo.R](file:///Users/Aamir/Desktop/rUTIs/12c_panaroo.R)** - Pangenome analysis
- **[12d_plasmid_analysis.R](file:///Users/Aamir/Desktop/rUTIs/12d_plasmid_analysis.R)** - Plasmid analysis
- **[12d_sv_analysis.R](file:///Users/Aamir/Desktop/rUTIs/12d_sv_analysis.R)** - Structural variant analysis
- **[12e_generate_reports.R](file:///Users/Aamir/Desktop/rUTIs/12e_generate_reports.R)** - Report generation
- **[13_visualise_panaroo_selection.R](file:///Users/Aamir/Desktop/rUTIs/13_visualise_panaroo_selection.R)** - Panaroo selection visualization

### Genotype-Phenotype & Lineage (Phase 3)
- **[14_genotype_phenotype_model.R](file:///Users/Aamir/Desktop/rUTIs/14_genotype_phenotype_model.R)** - GWAS for UTI-associated genes
- **[17_lineage_analysis.R](file:///Users/Aamir/Desktop/rUTIs/17_lineage_analysis.R)** - Lineage risk analysis

### Longitudinal & Evolution (Phase 4)
- **[15_longitudinal_patterns.R](file:///Users/Aamir/Desktop/rUTIs/15_longitudinal_patterns.R)** - Patient timelines
- **[16_within_host_evolution.R](file:///Users/Aamir/Desktop/rUTIs/16_within_host_evolution.R)** - Within-host evolution
- **[18_annotate_variants.R](file:///Users/Aamir/Desktop/rUTIs/18_annotate_variants.R)** - Variant annotation
- **[19_host_context.R](file:///Users/Aamir/Desktop/rUTIs/19_host_context.R)** - Host context analysis
- **[20_variant_annotation_deep.R](file:///Users/Aamir/Desktop/rUTIs/20_variant_annotation_deep.R)** - Deep variant annotation
- **[21_publication_figures.R](file:///Users/Aamir/Desktop/rUTIs/21_publication_figures.R)** - Publication-ready figures

---

## 🔧 Main Runner Scripts

- **[RUN_COMPLETE_ANALYSIS.sh](file:///Users/Aamir/Desktop/rUTIs/RUN_COMPLETE_ANALYSIS.sh)** - **⭐ PRIMARY** Complete pipeline execution (Phases 0-4)
- **[RUN_IN_TERMINAL.sh](file:///Users/Aamir/Desktop/rUTIs/RUN_IN_TERMINAL.sh)** - Alternative runner
- **[run_wgs.sh](file:///Users/Aamir/Desktop/rUTIs/run_wgs.sh)** - Legacy WGS runner (calls 12_wgs_exact_compare.R)
- **[run_12.sh](file:///Users/Aamir/Desktop/rUTIs/run_12.sh)** - Legacy WGS runner alternative
- **[11_cleanup.sh](file:///Users/Aamir/Desktop/rUTIs/11_cleanup.sh)** - Cleanup utility
- **[13_deep_clean.sh](file:///Users/Aamir/Desktop/rUTIs/13_deep_clean.sh)** - Deep cleanup

---

## ⚙️ Environment & Configuration

- **[env-annot.yml](file:///Users/Aamir/Desktop/rUTIs/env-annot.yml)** - Annotation environment
- **[env-wgs.yml](file:///Users/Aamir/Desktop/rUTIs/env-wgs.yml)** - WGS tools environment
- **[environment-rutis-core.yml](file:///Users/Aamir/Desktop/rUTIs/environment-rutis-core.yml)** - Core R environment
- **.gitignore** - Git ignore rules (updated to allow plots/)
- **rUTIs.Rproj** - RStudio project file

---

## 📊 Key Data Files (Root)

> [!TIP]
> These are generated outputs or required inputs for the pipeline

- **assemblies.list** - List of assembly files
- **assembly_metadata.csv** - Assembly metadata (generated by 00_make_assembly_metadata.r)
- **status_map.csv** - Clinical status map (PRIMARY output from Phase 0)
- **status_map_full.csv** - Full status map
- **vf_pa_all.csv** - Virulence factor presence/absence matrix
- **mlst_all.tsv** - MLST typing results

---

## 📁 R/ - Helper Functions

Contains reusable R functions used across multiple scripts:

- **[R/wgs_helpers.R](file:///Users/Aamir/Desktop/rUTIs/R/wgs_helpers.R)** - WGS analysis helpers
- **[R/clinical_helpers.R](file:///Users/Aamir/Desktop/rUTIs/R/clinical_helpers.R)** - Clinical data helpers
- **[R/plot_helpers.R](file:///Users/Aamir/Desktop/rUTIs/R/plot_helpers.R)** - Plotting utilities
- **[R/qc_select_panaroo_samples.R](file:///Users/Aamir/Desktop/rUTIs/R/qc_select_panaroo_samples.R)** - Panaroo sample QC
- **[R/wgs_report_template.Rmd](file:///Users/Aamir/Desktop/rUTIs/R/wgs_report_template.Rmd)** - Report template

---

## 📁 scripts/ - Utilities & Shell Scripts

Utility scripts and helpers organized by type:

### Main Utilities
- Various analysis and utility scripts

### AWK Scripts
- **[scripts/awk/inject_safe.awk](file:///Users/Aamir/Desktop/rUTIs/scripts/awk/inject_safe.awk)** - Safe code injection
- **[scripts/awk/patch_progress.awk](file:///Users/Aamir/Desktop/rUTIs/scripts/awk/patch_progress.awk)** - Progress patching
- **[scripts/awk/clean_gff.awk](file:///Users/Aamir/Desktop/rUTIs/scripts/awk/clean_gff.awk)** - GFF cleaning

### Cleanup Scripts
- **[scripts/cleanup_plan.sh](file:///Users/Aamir/Desktop/rUTIs/scripts/cleanup_plan.sh)** - Cleanup planning
- **[scripts/cleanup_root_clutter.sh](file:///Users/Aamir/Desktop/rUTIs/scripts/cleanup_root_clutter.sh)** - Root cleanup

---

## 📁 docs/ - Documentation

Comprehensive project documentation:

### Main Documentation
- **[docs/LECTURER_README.md](file:///Users/Aamir/Desktop/rUTIs/docs/LECTURER_README.md)** - Overview for lecturers
- **[docs/DOCUMENTATION_INDEX.md](file:///Users/Aamir/Desktop/rUTIs/docs/DOCUMENTATION_INDEX.md)** - Index of all docs
- **[docs/FINAL_SUMMARY.md](file:///Users/Aamir/Desktop/rUTIs/docs/FINAL_SUMMARY.md)** - Project summary
- **[docs/YELLOW_RoUTIne_R_Script_Breakdown.txt](file:///Users/Aamir/Desktop/rUTIs/docs/YELLOW_RoUTIne_R_Script_Breakdown.txt)** - Script descriptions

### Technical Documentation
- **[docs/pipeline_architecture.md](file:///Users/Aamir/Desktop/rUTIs/docs/pipeline_architecture.md)** - Pipeline architecture
- **[docs/script_and_figure_guide.md](file:///Users/Aamir/Desktop/rUTIs/docs/script_and_figure_guide.md)** - Scripts and figures
- **[docs/analysis_outputs_guide.md](file:///Users/Aamir/Desktop/rUTIs/docs/analysis_outputs_guide.md)** - Output guide
- **[docs/PIPELINE_FAILURE_LOG.md](file:///Users/Aamir/Desktop/rUTIs/docs/PIPELINE_FAILURE_LOG.md)** - Failure tracking

### Scientific Documentation
- **[docs/research_outcomes.md](file:///Users/Aamir/Desktop/rUTIs/docs/research_outcomes.md)** - Research findings
- **[docs/research_roadmap.md](file:///Users/Aamir/Desktop/rUTIs/docs/research_roadmap.md)** - Future plans
- **[docs/methods_draft.md](file:///Users/Aamir/Desktop/rUTIs/docs/methods_draft.md)** - Methods draft

### Subdirectories
- **docs/figures/** - Figure recommendations and designs
- **docs/papers/** - Related papers
- **docs/slides/** - Presentation materials

---

## 📁 plots/ - Curated Figures for Lecturer

> [!IMPORTANT]
> This directory contains only the KEY, PUBLICATION-READY figures intended for presentation to your lecturer. All diagnostic/intermediate plots are in subdirectories or legacy/.

### Structure
```
plots/
├── publication/         ⭐ Main publication figures
├── clinical/            Clinical summary plots
├── genomics/            Genomic analysis plots
├── phylogeny/           Phylogenetic trees
├── wgs/                 WGS QC and selection plots
├── mlst/                MLST analysis
├── plasmids/            Plasmid networks and heatmaps
├── vf/                  Virulence factor analysis
├── models/              Statistical model outputs
├── timelines/           Patient timelines
├── epidemiology/        Epidemiological plots
├── persistence/         Persistence patterns
├── strain_compare/      Within-host comparisons
├── pairwise_identity/   Pairwise comparisons
└── participant_specific/ Per-participant diagnostic plots
```

### Key Figures in plots/publication/
1. **Fig1_Swimmer_Plot.png** - Patient episode timelines
2. **Fig2_Mutation_Map.png** - Mutation patterns

### Important Clinical Plots in plots/clinical/
- **trajectories_heatmap.png** - Episode trajectories
- **waterfall_counts.png** - Episode counts by participant
- **transitions_alluvial_or_heatmap.png** - UTI/ASB transitions

### Core Phylogeny in plots/phylogeny/
- **core_tree_phenotype.png** - Core SNP tree colored by primary UTI status where available

---

## 📁 legacy/ - Old Versions & Debug Material

> [!CAUTION]
> These files are kept for reference but are NOT used in the current pipeline.

### Subdirectories
- **legacy/11_compare_variants/** - Old versions of 11_compare_strains.R
  - `11_compare_strains_MOD.R`
- **legacy/12_wgs_variants/** - Old versions of 12_wgs scripts
  - `12_wgs_exact_compare.RY`
  - `12_wgs_exact_compare_clean.R`
  - `12_wgs_exact_compare_clean.R.prepatch*` (multiple)
  - `12_wgs_runner.R`
  - `12e_report.R`
- **legacy/debug_scripts/** - Debug and experimental scripts
  - `debug_panaroo*.sh` (3 variants)
- **legacy/temp_dev_files/** - Temporary development files
  - History files (.Rapp.history, .Rhistory)
  - Test scripts (test_metadata_join.R)
  - Temp scripts (temp_generate_pairs*.R)
  - Utility scripts (verify_*.R, repair_prokka.R, install_r_packages.R)
- **legacy/installers/** - Large installer files
  - `miniforge_x86.sh`
  - `igraph_*.tgz`
  - `rUTIs_backup_*.tar.gz`
  - `Rplots.pdf`
  - Old log files

---

## 📁 Other Directories

### data/
Input data files (see `.gitignore` - most data is ignored for size)

### results/
Pipeline output files organized by analysis type

### logs/
Log files from pipeline execution

### tests/
Testing scripts and fixtures

### tools/
Analysis tools and utilities
- **tools/used_files.json** - Dependency map of all used files

### backups/
Backup files

### gff/
GFF annotation files

### not-assemblies/, ont-yellow-routine-fastas/, extras/
Data directories (see `.gitignore`)

---

## 🔍 Canonical Scripts Reference

> [!NOTE]
> If you're unsure which script to use, these are the canonical versions:

| Family | Canonical Script | Location | Used By |
|--------|-----------------|----------|---------|
| **11** | `11_compare_strains.R` | Root | RUN_COMPLETE_ANALYSIS.sh |
| **12** | `12_wgs_exact_compare.R` | Root | run_wgs.sh, run_12.sh |
| **12 modules** | `12a-12e` scripts | Root | RUN_COMPLETE_ANALYSIS.sh |

**Legacy variants** are in `legacy/11_compare_variants/` and `legacy/12_wgs_variants/`

---

## 🚀 Quick Start

1. **To run the complete pipeline:**
   ```bash
   bash RUN_COMPLETE_ANALYSIS.sh
   ```

2. **To view main figures:**
   - Navigate to `plots/publication/`
   - See also `plots/clinical/`, `plots/phylogeny/`, and `plots/wgs/`

3. **To understand outputs:**
   - Read `docs/analysis_outputs_guide.md`
   - Check `docs/script_and_figure_guide.md`

4. **For your lecturer:**
   - Start with `docs/LECTURER_README.md`
   - Main figures in `plots/publication/`
   - Research outcomes in `docs/research_outcomes.md`

---

**Repository Cleanup completed:** 2025-12-02  
**Maintained by:** Repository cleanup automation

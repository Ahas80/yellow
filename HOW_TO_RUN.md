# How to Run the rUTIs Pipeline

This guide provides instructions for running the `rUTIs` Whole Genome Sequencing (WGS) analysis pipeline.

## Prerequisites

1.  **R (>= 4.0)** with required packages (`tidyverse`, `optparse`, `lme4`, `pheatmap`, etc.).
2.  **Conda** environment with bioinformatics tools (`abricate`, `mlst`, `plasmidfinder`, `parsnp`, `panaroo`, `snp-dists`).
    *   Activate environment: `conda activate wgs_env` (or similar).

## Execution Order

The pipeline is designed to be run sequentially. Scripts are numbered to indicate the order.

### Phase 1: Clinical Data & Setup
1.  `Rscript 00a_load_clean_clinical.R` - Load and clean clinical metadata.
2.  `Rscript 00b_classify_episodes.R` - Classify episodes (UTI vs ASB).
3.  `Rscript 00c_plot_clinical_summary.R` - Generate clinical summary plots.

### Phase 2: Bacterial Genomics (VF, MLST, Plasmids)
4.  `Rscript 02_gene_presence_analysis.R` - Run ABRicate for VFs and AMR genes.
5.  `Rscript 04_gene_breakdown.R` - Analyze gene prevalence and associations.
6.  `Rscript 05_gene_overview_plots.R` - Visualize gene presence.
7.  `Rscript 06_MLST.R` - Run MLST typing.
8.  `Rscript 07_explore_MLST.R` - Explore MLST distribution.
9.  `Rscript 08_core_vs_plasmid.R` - Compare core STs with plasmid types.
10. `Rscript 09_inc_plasmid_network.R` - Analyze plasmid networks.
11. `Rscript 10_replicon_heatmap.R` - Generate plasmid replicon heatmaps.

### Phase 3: WGS Core Pipeline
12. `Rscript 12a_wgs_qc.R` - Perform assembly QC.
13. `Rscript 12b_core_snp.R` - Call core SNPs and calculate distances.
14. `Rscript 12c_panaroo.R` - Run Pangenome analysis (Panaroo).
15. `Rscript 13_visualise_panaroo_selection.R` - Visualize QC and selection bias.

### Phase 4: Integration & Reporting
16. `Rscript 11_compare_strains.R` - Compare specific strain pairs.
17. `Rscript 14_genotype_phenotype_model.R` - Run GWAS (Genotype-Phenotype Association).
18. `Rscript 12e_generate_reports.R` - Generate per-participant PDF reports.

## Helper Shell Scripts

In addition to the R scripts, the following shell scripts are available for maintenance and alternative analyses:

| Script | Purpose | When to Use |
| :--- | :--- | :--- |
| `11_cleanup.sh` | Standard cleanup. Moves plots to `docs/` and legacy files to `results/legacy/`. | Run this **after** completing the full pipeline to tidy up the project folder. |
| `13_deep_clean.sh` | Deep cleanup. Archives logs, temp files, and old scripts to `results/cleanup_<date>/`. | Use when you want to "reset" the workspace or archive a previous run's artifacts. |
| `run_wgs.sh` | **Legacy** runner for the old monolithic WGS script (`12_wgs_exact_compare.R`). | Use only if you need to reproduce results from the old single-script pipeline. For new runs, use the modular R scripts (12a-12e). |
| `01_kmer_mash.sh` | Standalone Mash analysis. | Use for quick K-mer based clustering independent of the main R pipeline. |
| `02_core_snp.sh` | Standalone Parsnp analysis. | Use if you prefer Parsnp over the pipeline's `minimap2`/`bcftools` approach for core SNPs. |

## Configuration

All file paths and constants are defined in `00_config.R`. Modify this file to change input directories, QC thresholds, or output locations.

## Troubleshooting

- **Missing Files**: Check `FOLDER_MAP.md` to see where files should be.
- **Logs**: Check `logs/debug/` for detailed error messages.
- **Environment**: Ensure all external tools (`abricate`, `parsnp`, etc.) are in your PATH.

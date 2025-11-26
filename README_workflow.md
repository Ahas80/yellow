# rUTIs Pipeline Workflow

This guide details the execution order, environment requirements, and expected outputs for the refactored rUTIs pipeline.

## Environment Setup

1.  **R Environment**: Ensure you have R installed with the following packages:
    *   `tidyverse` (dplyr, readr, tidyr, stringr, ggplot2, purrr, forcats)
    *   `fs`, `glue`, `Biostrings`
    *   `furrr` (for parallel processing)
    *   `ComplexHeatmap` or `pheatmap` (for heatmaps)
    *   `optparse`, `processx`

2.  **Conda Environment** (for shell scripts):
    *   Ensure you have a conda environment (e.g., `asm-snp-x86` or similar) with:
        *   `mash`
        *   `parsnp`, `harvesttools`
        *   `snp-dists`
        *   `abricate`, `nucmer` (mummer)
        *   `mlst`, `blastn`

## Execution Order

Run the scripts in the following order. Always run from the project root (`Desktop/rUTIs/`).

### Phase 1: Foundation & Metadata

1.  **`00_process_clinical_data.R`** (formerly `01_prepare_assembly_metadata.R`)
    *   *Purpose*: Processes `batch*.csv` clinical data and creates the ID mapping.
    *   *Input*: `batch1.csv`, `batch2.csv`, `batch3.csv`.
    *   *Output*: `results/assembly_pid_timepoint_map.csv`, `results/debug/`.

2.  **`00_make_assembly_metadata.r`**
    *   *Purpose*: Combines FASTA files with clinical metadata to create the master index.
    *   *Input*: `ont-yellow-routine-fastas/`, `results/assembly_pid_timepoint_map.csv`.
    *   *Output*: `assembly_metadata.csv` (root), `results/assembly_metadata.csv`.

3.  **`01_kmer_mash.sh`**
    *   *Purpose*: Runs Mash for K-mer distance estimation.
    *   *Input*: `assembly_metadata.csv`, FASTA files.
    *   *Output*: `results/kmer/mash_distance_matrix.csv`.

4.  **`02_core_snp.sh`**
    *   *Purpose*: Runs Parsnp for core SNP alignment.
    *   *Input*: FASTA files.
    *   *Output*: `results/parsnp/parsnp.tree`, `results/parsnp/parsnp.ggr`.

5.  **`02_gene_presence_analysis.R`**
    *   *Purpose*: Runs Abricate (VFDB) and Nucmer for gene/plasmid detection.
    *   *Input*: FASTA files.
    *   *Output*: `results/abricate/`, `results/nucmer/`.

### Phase 2: Core Analysis

6.  **`03_plotting.R`**
    *   *Purpose*: Generates high-level descriptive plots (richness, trajectories).
    *   *Output*: `results/plots/`.

7.  **`04_gene_breakdown.R`**
    *   *Purpose*: Detailed gene analysis (nitrate, focus genes).
    *   *Output*: `results/annotated_gene_table.csv`, `results/nitrate_presence_matrix.csv`.

8.  **`05_gene_overview_plots.R`**
    *   *Purpose*: Heatmaps and prevalence plots for genes.
    *   *Output*: `results/plots/variable_gene_heatmap.pdf`.

9.  **`06_MLST.R`**
    *   *Purpose*: Runs MLST typing.
    *   *Output*: `results/mlst/mlst_all.tsv`.

10. **`07_explore_MLST.R`**
    *   *Purpose*: Visualizes MLST results (ST frequency).
    *   *Output*: `results/mlst/top20_STs.pdf`.

11. **`08_core_vs_plasmid.R`**
    *   *Purpose*: Compares chromosomal vs. plasmid STs/replicons.
    *   *Output*: `results/plasmid_replicons_wide.csv`.

12. **`09_inc_plasmid_network.R`**
    *   *Purpose*: Network analysis of plasmid replicons.
    *   *Output*: `results/replicon_cooccurrence.pdf`.

### Phase 3: Deep Dive & Cleanup

13. **`10.R`**
    *   *Purpose*: Plasmid replicon heatmap.
    *   *Output*: `results/replicon_heatmap.pdf`.

14. **`11_cleanup.sh`**
    *   *Purpose*: Moves legacy files and plots to `docs/` and `results/legacy/`.
    *   *Run*: `bash 11_cleanup.sh`

15. **`12_wgs_exact_compare.R`**
    *   *Purpose*: Detailed pairwise WGS comparison (SNPs, SVs).
    *   *Run*: `Rscript 12_wgs_exact_compare.R` (can take a long time).

## Output Structure

*   `results/`
    *   `plots/`: All visualizations.
    *   `mlst/`: MLST tables and QC.
    *   `kmer/`: Mash results.
    *   `parsnp/`: Core SNP trees.
    *   `wgs/`: Deep dive outputs (from script 12).
    *   `legacy/`: Old files moved by cleanup script.
*   `docs/`: Papers and slides.

## Troubleshooting

*   **Missing Metadata**: If `00_make_assembly_metadata.r` warns about rows needing manual fill, check `results/assembly_metadata_TODO_fill.csv`.
*   **Tool Not Found**: Ensure your conda environment is active (`conda activate asm-snp-x86`).
*   **Permission Denied**: Run `chmod +x *.sh` to make shell scripts executable.

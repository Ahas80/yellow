# rUTIs Pipeline Architecture & Data Flow

This document explains how every script in the project relates to each other, why they are necessary, and exactly what they produce.

---

## **Phase 1: Clinical Data Foundation**
*Goal: Turn raw hospital CSVs into a clean, usable dataset for analysis.*

### `00_config.R`
*   **Role**: The "Brain".
*   **Why**: Defines global paths, constants, and helper functions used by **every** other script. Ensures consistency.
*   **Input**: None.
*   **Output**: Global R environment variables (`DIR_RESULTS`, `CORES_USE`, etc.).

### `00a_load_clean_clinical.R`
*   **Role**: The "Loader".
*   **Why**: Raw data comes in messy batches (`batch1.csv`, `batch2.csv`). This script merges them and fixes column names.
*   **Input**: `data/inputs/batch*.csv`.
*   **Output**: `results/clinical/intermediate/clinical_merged.rds` (Raw merged data).

### `00b_classify_episodes.R`
*   **Role**: The "Judge".
*   **Why**: We need to decide if a patient has a **UTI**, **ASB**, or is **Negative**. This script applies the clinical rules (Symptoms + CFU count).
*   **Input**: `clinical_merged.rds`.
*   **Output**: `status_map.csv` (The "Master Key": maps every `Participant_id` + `Timepoint` to a clinical status).

### `00c_plot_clinical_summary.R`
*   **Role**: The "Reporter".
*   **Why**: Visualizes the cohort. How many UTIs? How many ASB?
*   **Input**: `status_map.csv`.
*   **Output**: `plots/clinical/cohort_summary.png`.

---

## **Phase 2: WGS Pipeline (The Heavy Lifting)**
*Goal: Process raw DNA sequencing data into annotated genomes.*

### `12_wgs_exact_compare.R` (and `12_wgs_runner.R`)
*   **Role**: The "Factory".
*   **Why**: This is the main bioinformatics pipeline. It wraps multiple external tools (minimap2, samtools, mash, panaroo).
*   **Input**: Raw FASTQ reads (`data/reads`) or Assemblies (`ont-yellow-routine-fastas`).
*   **Output**:
    *   `results/wgs/core/`: Core SNP alignments.
    *   `results/wgs/kmer/`: MASH distances (genomic similarity).
    *   `results/wgs/pangenome/`: Gene presence/absence matrix (`gene_data.csv`).

### `02_gene_presence_analysis.R`
*   **Role**: The "Librarian".
*   **Why**: Takes the raw gene data and organizes it into a clean "Presence/Absence" matrix for statistical testing.
*   **Input**: `results/wgs/pangenome/gene_data.csv`.
*   **Output**: `results/vf/vf_pa_all.csv` (Rows = Isolates, Cols = Genes).

### `03_plotting.R`
*   **Role**: The "Visualizer".
*   **Why**: Creates heatmaps and summary plots of gene presence across samples.
*   **Input**: `vf_pa_all.csv`.
*   **Output**: Various plots in `plots/vf/`.

### `06_MLST.R`
*   **Role**: The "Taxonomist".
*   **Why**: Determines the bacterial lineage (Sequence Type) of each isolate using 7 housekeeping genes.
*   **Input**: Assemblies.
*   **Output**: `results/mlst/mlst_all.tsv`.

---

## **Phase 3: Comparative Genomics**
*Goal: Compare strains to find "Same" vs "Different" and identify risk factors.*

### `11_compare_strains.R` (and `11_compare_strains_helpers.R`)
*   **Role**: The "Detective".
*   **Why**: Compares every isolate against every other isolate from the same patient.
*   **How**: Uses ANI (Average Nucleotide Identity), SNP distance, VF Jaccard, Plasmid Jaccard.
*   **Input**:
    *   `status_map.csv` (Clinical info).
    *   `mash_distance.csv` (Genomic distance).
    *   `vf_pa_all.csv` (Gene content).
*   **Output**: `results/strain_compare/pairwise_metrics.csv` (Classifies pairs as "Same Strain" or "Different").

### `14_genotype_phenotype_model.R`
*   **Role**: The "Statistician".
*   **Why**: Performs GWAS (Genome-Wide Association Study). Tests if specific genes (e.g., *lpfA*) are more common in UTI vs ASB.
*   **Input**: `status_map.csv` + `vf_pa_all.csv`.
*   **Output**: 
    *   `results/models/volcano_plot.png`.
    *   `gwas_multivariable_glmm.csv` (Statistical results).

### `17_lineage_analysis.R`
*   **Role**: The "Epidemiologist".
*   **Why**: Checks if specific bacterial lineages (Sequence Types like ST131) are more dangerous.
*   **Input**: `results/mlst/mlst_all.tsv` + `status_map.csv`.
*   **Output**: `results/lineage/st_risk_profile.csv`.

---

## **Phase 4: Longitudinal & Mechanism (The New Roadmap)**
*Goal: Track evolution over time and explain phenotype switches.*

### `15_longitudinal_patterns.R`
*   **Role**: The "Time Traveler".
*   **Why**: Stitches individual episodes into patient timelines. Identifies "Phenotype Switches" (ASB → UTI).
*   **How**: Uses graph-based clustering to assign global "Strain IDs" based on "Same" classification.
*   **Input**: `pairwise_metrics.csv` (from script 11).
*   **Output**: 
    *   `results/longitudinal/participant_timelines.csv`.
    *   `phenotype_switch_candidates.csv`.
    *   `swimmer_plot.png`.

### `16_within_host_evolution.R`
*   **Role**: The "Microscope".
*   **Why**: Zooms in on the "Switch" candidates. Calculates exact SNP distance and gene gain/loss between T1 and T2.
*   **How**: Uses `nucmer` (from MUMmer4) to align assemblies and count SNPs.
*   **Input**: `phenotype_switch_candidates.csv`.
*   **Output**: 
    *   `evolution_events.csv` (List of specific changes).
    *   `evolution_summary.txt`.

### `18_annotate_variants.R`
*   **Role**: The "Translator".
*   **Why**: `16` gives us SNP positions. This script parses the raw `.snps` files from nucmer.
*   **Input**: `.snps` files (from `nucmer`).
*   **Output**: `annotated_snps.csv` (Position, Ref Base, Alt Base, Type = SNP/Insertion/Deletion).

### `19_host_context.R`
*   **Role**: The "Clinician".
*   **Why**: Checks if the patient had a catheter or antibiotics during the switch.
*   **Input**: 
    *   `phenotype_switch_candidates.csv`.
    *   `clinical_merged.rds`.
*   **Output**: `host_context_table.csv` (Catheter status, symptoms for each switch).

### `20_variant_annotation_deep.R`
*   **Role**: The "Biologist".
*   **Why**: Deep dive into specific mutations. Maps SNPs to GFF files to find gene products (e.g., "Sigma 70").
*   **How**: Parses Prokka GFF files to match SNP positions to genes.
*   **Input**: 
    *   `annotated_snps.csv`.
    *   GFF files from `results/prokka_prefixed_slim/`.
*   **Output**: `variant_annotation_detailed.csv` (SNP position → Gene name → Product).

### `21_publication_figures.R`
*   **Role**: The "Artist".
*   **Why**: Generates the final, polished figures for the paper.
*   **Input**: All results above.
*   **Output**: 
    *   `plots/publication/Fig1_Swimmer_Plot.png`.
    *   `plots/publication/Fig2_Mutation_Map.png`.

---

## **Summary of Data Flow**

```
1. Raw Clinical Data → 00a/00b → status_map.csv (The Master Key)
                                      ↓
2. Raw DNA → 12_wgs → Gene Matrices (vf_pa_all.csv, MASH, etc.)
                                      ↓
3. status_map + Gene Matrices → 11_compare → pairwise_metrics.csv
                                                    ↓
4. pairwise_metrics → 15_longitudinal → Timelines & phenotype_switch_candidates.csv
                                                    ↓
5. phenotype_switch_candidates → 16/18/20 → Mechanism (SNPs → Genes → rpoD, lpxL)
                                                    ↓
6. Everything → 21_figures → Publication Figures
```

---

## **Key Files to Understand the Project**

1.  **`status_map.csv`**: Every participant-timepoint combination and their clinical status (UTI/ASB/Negative).
2.  **`vf_pa_all.csv`**: Every gene and which samples have it (binary matrix).
3.  **`pairwise_metrics.csv`**: Every pair of samples compared (genomic distance + classification as "Same" or "Different").
4.  **`participant_timelines.csv`**: Longitudinal view showing how strains persist or switch phenotypes.
5.  **`variant_annotation_detailed.csv`**: The exact mutations (e.g., *rpoD* SNP) that occurred during phenotype switches.

---

## **Why This Architecture?**

*   **Modularity**: Each script does ONE thing well. You can re-run step 15 without re-running the WGS pipeline.
*   **Reproducibility**: All paths are in `00_config.R`. Anyone can run this on their machine.
*   **Layered Analysis**: Static → Comparative → Longitudinal → Mechanistic. Each layer builds on the previous.

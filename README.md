# rUTIs Pipeline - Yellow Cohort Analysis

**Escherichia coli** genomic and clinical analysis for recurrent UTI (rUTIs) study.

---

## 🔧 Setup (First Time Only)

### 1. Install Conda Environment
```bash
# Create environment with bioinformatics tools
conda create -n yellow-wgs-x86 -c bioconda -c conda-forge \
    parsnp panaroo snp-dists abricate mlst plasmidfinder mummer prokka openjdk

# Activate it
conda activate yellow-wgs-x86
```

### 2. Install R Packages
**CRITICAL**: Run this before the first pipeline execution:
```bash
# Install all required R packages (takes 5-10 minutes)
Rscript install_r_packages.R
```

This installs 38 CRAN packages (including `ggraph`, `tidyverse`, `lme4`, etc.) and 4 Bioconductor packages.

**Verify installation**:
```r
# In R, check key packages load:
library(ggraph)     # Network plots
library(tidyverse)  # Data manipulation
library(lme4)       # Mixed models
library(ComplexUpset)  # UpSet plots
```

If any package is missing, install individually:
```bash
Rscript -e "install.packages('ggraph', repos='https://cloud.r-project.org')"
```

---

## Quick Start

```bash
# 1. Activate conda environment
conda activate yellow-wgs-x86

# 2. Run COMPLETE pipeline (all phases)
bash RUN_COMPLETE_ANALYSIS.sh

# OR run phases individually (see below)
```

---

## 🔄 Complete Analysis Workflow

The pipeline is divided into **4 phases**. Run them in order:

### **Phase 0: Clinical Data Foundation** (Required First)
```bash
Rscript 00a_load_clean_clinical.R
Rscript 00b_classify_episodes.R
Rscript 00c_plot_clinical_summary.R
```
**Output**: `results/clinical/status_map.csv` (Master clinical status table)

---

### **Phase 1: WGS Processing** (Long-running)
```bash
# QC and core SNPs
Rscript 12a_wgs_qc.R
Rscript 12b_core_snp.R        # ~30-60 mins

# Pangenome
Rscript 12c_panaroo.R         # ~15 mins
Rscript 13_visualise_panaroo_selection.R

# Gene presence/absence
Rscript 02_gene_presence_analysis.R
Rscript 06_MLST.R
```
**Output**: 
- `results/wgs/qc_summary.csv`
- `results/wgs/core/` (SNP alignments)
- `results/wgs/pan/gene_data.csv` (Pangenome)
- `results/vf/vf_pa_all.csv` (Gene matrix)
- `results/mlst/mlst_all.tsv` (Lineages)

---

### **Phase 2: Comparative Genomics**
```bash
# Within-host comparisons
Rscript 11_compare_strains.R --participants ALL

# GWAS (find UTI-associated genes)
Rscript 14_genotype_phenotype_model.R

# Lineage risk
Rscript 17_lineage_analysis.R
```
**Output**:
- `results/strain_compare/pairwise_metrics.csv` (Same vs Different strains)
- `results/models/volcano_plot.png` (GWAS results)
- `results/lineage/st_risk_profile.csv`

---

### **Phase 3: Longitudinal & Mechanism** (NEW - Research Roadmap)
```bash
# Reconstruct patient timelines
Rscript 15_longitudinal_patterns.R

# Identify mutations in phenotype switches
Rscript 16_within_host_evolution.R
Rscript 18_annotate_variants.R
Rscript 20_variant_annotation_deep.R

# Check host factors (catheter, antibiotics)
Rscript 19_host_context.R

# Generate publication figures
Rscript 21_publication_figures.R
```
**Output**:
- `results/longitudinal/participant_timelines.csv` (Timeline data)
- `results/longitudinal/phenotype_switch_candidates.csv` (ASB→UTI switches)
- `results/longitudinal/variant_annotation_detailed.csv` (Mutated genes)
- `plots/publication/*.png` (Figures 1-2)

---

## 📊 Key Results

**Original Finding (Phase 2)**:
- 🔬 **Top Hit**: *Long polar fimbriae* (`lpfA`, `lpfB`) strongly associated with symptomatic UTI (OR=7.6, FDR<0.001)

**NEW Findings (Phase 3 - Longitudinal)**:
- 🧬 **"The Chameleon Effect"**: Found 2 cases where the *exact same strain* caused ASB at T1 and UTI at T2
- 🧫 **Mechanism**: Phenotype switches involved minimal genomic changes (10-12 SNPs), NOT gene acquisition
- 🎯 **Key Mutations**:
  - Participant 40004: `rpoD` (Sigma 70) mutation → global transcriptional shift
  - Participant 40001: `lpxL` (Lipid A biosynthesis) mutation → endotoxin modification
- ❌ **Lineage Risk**: No ST was significantly more dangerous (challenges "high-risk clone" approach)

---

## 📂 Essential Documentation

| Document | Purpose |
|:---------|:--------|
| **`docs/pipeline_architecture.md`** | How scripts relate to each other |
| **`docs/research_outcomes.md`** | Comprehensive findings summary |
| **`docs/methods_draft.md`** | Methods section for manuscript |
| **`results/KEY_FINDINGS.md`** | Executive summary (Lpf association) |
| **`results/ANALYSIS_README.md`** | Navigate results files |
| **`FINAL_SUMMARY.md`** | General audience summary |

---

## Project Structure

| Directory | Purpose |
|:----------|:--------|
| `data/inputs/` | Raw clinical CSVs (batch1-3) |
| `results/clinical/` | Processed clinical data |
| `results/vf/` | Virulence factor matrices |
| `results/wgs/` | WGS outputs (SNPs, pangenome) |
| `results/strain_compare/` | Within-host comparisons |
| `results/models/` | GWAS results |
| `results/longitudinal/` | **NEW** - Timelines and evolution |
| `results/lineage/` | **NEW** - ST risk analysis |
| `plots/publication/` | **NEW** - Manuscript figures |
| `logs/` | Execution logs |
| `ont-yellow-routine-fastas/` | Assembly FASTA files |

---

## solved : Participant ID Logic - deprecated R file so this is solved 

**Problem**: Legacy pipeline (`12_wgs_exact_compare.R`) treated batch IDs (e.g., `PR0010`) as single participants.

**Solution**: Use modular pipeline (`12a-12e`) which correctly reads `Participant_id` from `assembly_metadata.csv`.

**Verification**:
```bash
# Should show ~87 unique participants (not 5 batch IDs)
cut -d',' -f8 results/wgs/qc_summary.csv | tail -n +2 | sort -u | wc -l
```

---

## Configuration

All paths in `00_config.R`:

```r
DIR_RESULTS <- "results"
DIR_PLOTS <- "plots"
CORES_USE <- 10  # Adjust for your system
```

---

## Dependencies

### R Packages (38 CRAN + 4 Bioconductor)
**Install via**: `Rscript install_r_packages.R`

**Key packages**: 
- Data: `tidyverse`, `data.table`, `vroom`
- Modeling: `lme4`, `broom.mixed`
- Visualization: `ggplot2`, `ggraph`, `pheatmap`, `ComplexUpset`, `ggtree`
- Parallel: `future`, `future.apply`, `furrr`

### Conda Tools
**Install via**:
```bash
conda create -n yellow-wgs-x86 -c bioconda -c conda-forge \
    parsnp panaroo snp-dists abricate mlst plasmidfinder mummer prokka openjdk
```

**Tools**: `parsnp` (SNPs), `panaroo` (pangenome), `abricate` (AMR/VF), `mlst` (typing), `prokka` (annotation)

---

## Troubleshooting

### "No GFF files found" (Script 20)
- **Cause**: GFF files in `results/prokka_*` directories
- **Fix**: Script auto-searches multiple directories
- **Note**: Missing GFFs for some assemblies is normal (uses fallback)

### Parsnp slow/hanging
- **Cause**: Large dataset (361 genomes)
- **Solution**: Expected, allow 30-60 mins

### "phenotype_switch_candidates.csv is empty"
- **Cause**: Need to re-run `11_compare_strains.R` after fixing Jaccard NA handling
- **Fix**: Run `Rscript 11_compare_strains.R --participants ALL`

---

## Citation & Contact

**Project**: Yellow rUTIs Cohort, E. coli Genomic Analysis  
**Date**: 2025-11-28  
**Status**: Initial analysis complete - now to figure out exactly what different E colis are and how they are different.

## 📝 To-Do & Future Directions

### 1. Investigate Within-Host E. coli Variation
**Question:** Are the *E. coli* isolates in a single participant different across timepoints (T0, T1, T2, Uricult)? - especially between UTIs vs ASBs (Uricult vs Timepoint X)

**Postulate & Approach:**
To determine if we are seeing the same strain persisting or new strains entering (strain replacement), we can use the following results:

1.  **Genomic Distance (SNPs)**:
    *   **Source:** `results/strain_compare/pairwise_metrics.csv`
    *   **Method:** Filter for `within_participant == TRUE`.
    *   **Thresholds:** TBD
        *   **< 10 SNPs**: Likely the **same strain** (persistent infection)? 
        *   **> 1000 SNPs**: Likely a **different strain** (re-infection/replacement)?
        *   **10-100 SNPs**: Grey area (potential within-host evolution)?

2.  **Gene Content (Pangenome)**:
    *   **Source:** `results/wgs/pan/gene_data.csv` (or `presence_absence.Rtab`)
    *   **Method:** Calculate Jaccard distance between isolates from the same patient.
    *   **Expectation:** Same strains should share >95% of their accessory genome.

3.  **Visual Check**:
    *   **Source:** `plots/publication/Fig1_Swimmer_Plot.png`
    *   **Method:** Look for patients with multiple timepoints. If the "Strain ID" (if annotated) changes, it's a replacement.

**Next Steps:**
- [ ] Create a boxplot of "Pairwise SNP Distance" for all within-patient pairs.
- [ ] Flag any patient with >0 but <50 SNPs for manual review (evolution candidates).

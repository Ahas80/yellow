# rUTIs Pipeline - Yellow Cohort Analysis

**Escherichia coli** genomic and clinical analysis for recurrent UTI (rUTIs) study.

---

## Quick Start

```bash
# 1. Activate conda environment
conda activate asm-snp-x86  # or your WGS env

# 2. Run modular WGS pipeline
Rscript 12a_wgs_qc.R          # QC assemblies
Rscript 12b_core_snp.R        # Core genome SNPs  
Rscript 12c_panaroo.R         # Pangenome
Rscript 13_visualise_panaroo_selection.R  # Selection plots

# 3. Integration & modeling
Rscript 11_compare_strains.R --participants ALL
Rscript 14_genotype_phenotype_model.R

# See docs/run_plan_wgs_restart.md for detailed execution plan
```

---

## Project Structure

| Directory | Purpose |
|:----------|:--------|
| `data/` | Clinical metadata, raw inputs |
| `results/` | Computational outputs (VF, MLST, GWAS, WGS) |
| `plots/` | Visualizations |
| `logs/` | Execution logs |
| `R/` | Helper scripts (`wgs_helpers.R`, etc.) |
| `docs/` | Documentation, run plans, status reports |
| `ont-yellow-routine-fastas/` | Assembly FASTA files |

---

## Key Results

📊 **Sample the data**: 87 participants, 361 QC-passed genomes  
🔬 **Top Finding**: *Long polar fimbriae* (`lpfA`, `lpfB`) strongly associated with symptomatic UTI

**Essential Documentation**:
- **`results/KEY_FINDINGS.md`**: Executive summary of Lpf association findings
- **`results/ANALYSIS_README.md`**: Navigation guide for all results files
- **`docs/wgs_restart_status_pass2.md`**: Latest WGS pipeline execution status

---

## Pipeline Phases

### Phase 0: Clinical (00a-00c)
1. `00a_load_clean_clinical.R` - Load/clean clinical data
2. `00b_classify_episodes.R` - Classify UTI vs ASB episodes
3. `00c_plot_clinical_summary.R` - Clinical plots

Outputs: `results/clinical/status_map.csv`

### Phase 1: Genomics (02-10)
4. `02_gene_presence_analysis.R` - Virulence factors (Abricate)
5. `04_gene_breakdown.R` - Focus gene analysis
6. `06_MLST.R` - Multi-locus sequence typing
7. `08_core_vs_plasmid.R` - Core vs plasmid comparison
8. `09_inc_plasmid_network.R` - Replicon network analysis

Outputs: `results/vf/vf_pa_all.csv`, `results/mlst/`, `results/plasmids/`

### Phase 2: WGS Core (12a-12e) **[MODULAR PATH]**
9. `12a_wgs_qc.R` - Assembly QC
10. `12b_core_snp.R` - Core genome SNP calling (Parsnp)
11. `12c_panaroo.R` - Pangenome analysis
12. `13_visualise_panaroo_selection.R` - QC/selection visualization

Outputs: `results/wgs/qc_summary.csv`, `results/wgs/core/`, `results/wgs/pan/`

⚠️ **Legacy**: `12_wgs_exact_compare.R` is deprecated (batch ID bug). Use modular 12a-12e instead.

### Phase 3: Integration (11, 14)
13. `11_compare_strains.R` - Within-host strain comparison
14. `14_genotype_phenotype_model.R` - Genotype-phenotype GWAS

Outputs: `results/strain_compare/`, `results/models/gwas_*.csv`

---

##  Critical Fix: Participant ID Logic

**Problem**: Legacy pipeline treated batch IDs (e.g., `PR0010`) as single participants, causing invalid cross-patient comparisons.

**Solution**: Modular pipeline (12a-12e) correctly uses `Participant_id` from `assembly_metadata.csv`.

**Verification**:
```bash
# Check QC summary has true participant IDs
cut -d',' -f8 results/wgs/qc_summary.csv | sort -u | wc -l
# Should return ~87, not 5

# Ensure no batch IDs in participant columns
grep "PR0010" results/strain_compare/pairwise_metrics.csv results/models/*.csv
# Should return nothing
```

---

## Configuration

All paths configured in `00_config.R`. Key settings:

```r
DIR_RESULTS <- "results"
DIR_PLOTS <- "results/plots"
DIR_WGS <- "results/wgs"
CORES_USE <- 10  # Adjust based on your system
```

---

## Dependencies

**R Packages**: `tidyverse`, `lme4`, `broom.mixed`, `optparse`, `pheatmap`, `furrr`, `seqinr`, `fs`

**Conda Tools**: `parsnp`, `panaroo`, `snp-dists`, `abricate`, `mlst`, `plasmidfinder`, `nucmer`

Install conda env:
```bash
conda create -n asm-snp-x86 -c bioconda parsnp panaroo snp-dists abricate mlst plasmidfinder mummer
```

---

## Troubleshooting

### Parsnp slow/hanging
- **Cause**: Large dataset (361 genomes)
- **Solution**: Expected, allow 30-60 mins for completion

### Tools not found in PATH
- **Fix**: Added in `R/wgs_helpers.R` (lines 17-25) - conda paths prepended to PATH

### file_symlink error
- **Fix**: Updated `12b_core_snp.R` to use `fs::link_create()`

See `docs/wgs_restart_status_pass2.md` for recent execution logs and fixes.

---

## Contact & Citation

**Project**: Yellow rUTIs Cohort, E. coli Genomic Analysis  
**Date**: 2025-11-27  
**Status**: Active - WGS pipeline completion in progress

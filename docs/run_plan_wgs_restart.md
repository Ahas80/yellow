# WGS Pipeline Restart Run Plan

## 1. Dependency Analysis

| Script | Changed? | Inputs | Outputs | Action | Reason |
| :--- | :--- | :--- | :--- | :--- | :--- |
| `00a-00c` (Clinical) | No | Raw Clinical Data | `status_map.csv` | **Skip** | Clinical data is stable. |
| `02_gene_presence` | No | Assemblies | `vf_pa_all.csv` | **Skip** | VF data is stable. |
| `12a_wgs_qc.R` | **Yes** (Logic check) | `assembly_metadata.csv` | `qc_summary.csv` | **RUN** | Must ensure `qc_summary.csv` has correct `Participant_id` from metadata. |
| `12b_core_snp.R` | No | `qc_summary.csv` | `strain_pairs.csv` | **RUN** | Depends on `qc_summary.csv`. Needs to use correct PIDs for labeling. |
| `12c_panaroo.R` | No | `qc_summary.csv` | `gene_presence.csv` | **RUN** | Depends on `qc_summary.csv` (valid genomes list). |
| `11_compare_strains.R` | **Yes** (Logic check) | `strain_pairs.csv`, Metadata | `pairwise_metrics.csv` | **RUN** | **CRITICAL**. Must use correct PIDs for within-host classification. |
| `14_genotype_phenotype` | **Yes** (Logic check) | `vf_pa_all.csv`, Metadata | `gwas_stats.csv` | **RUN** | **CRITICAL**. Must use correct PIDs for GLMM random effects. |
| `12e_reports.R` | No | All of above | PDFs | **RUN** | Needs to generate reports for correct participants. |

## 2. Participant Logic Verification

- **`R/wgs_helpers.R`**: `discover_samples` MUST join `assembly_metadata.csv` correctly.
    - *Status*: Verifying with `test_metadata_join.R`.
- **`12a_wgs_qc.R`**: Uses `meta_df` from `assembly_metadata.csv`.
    - *Status*: Safe, provided metadata is loaded correctly.
- **`11_compare_strains.R`**: Uses `core$assemblies` (from `discover_samples`).
    - *Status*: Safe if `discover_samples` is fixed.

## 3. Execution Order

1.  **Fix `R/wgs_helpers.R`** (if needed) to ensure metadata join works.
2.  **Run `12a_wgs_qc.R`** -> `results/wgs/qc_summary.csv`.
3.  **Run `12b_core_snp.R`** -> `results/wgs/core/strain_pairs.csv`.
4.  **Run `12c_panaroo.R`** -> `results/wgs/pan/gene_presence_absence.csv`.
5.  **Run `11_compare_strains.R`** -> `results/strain_compare/pairwise_metrics.csv`.
6.  **Run `14_genotype_phenotype_model.R`** -> `results/models/`.
7.  **Run `12e_generate_reports.R`** -> `results/reports/`.

## 4. Sanity Checks

- **QC**: `qc_summary.csv` should have ~87 unique participants, not 5.
- **Pairs**: `pairwise_metrics.csv` should show within-host pairs matching `Participant_id`.
- **GWAS**: `gwas_multivariable_glmm.csv` should run without "singular fit" warnings due to bad grouping.

---

# Pass 2: Complete Restart (After Kill)

## Context

**Date**: 2025-11-27 05:41 AM

**Previous session**:
- ✅ **Verified**: `discover_samples()` in `R/wgs_helpers.R` correctly uses Participant_id from `assembly_metadata.csv`
- ✅ **Fixed**: PATH handling for conda tools (parsnp, etc.)
- ✅ **Fixed**: `file_symlink` → `fs::link_create` in `12b_core_snp.R`
- ✅ **Completed**: `12a_wgs_qc.R` (361 genomes passed QC, correct Participant IDs)
- ⚠️ **Interrupted**: `12b_core_snp.R` killed mid-Parsnp alignment

## Dependency Analysis (Pass 2)

| Script | Last Modified | Key Inputs | Key Outputs | Status | Action |
|:-------|:-------------|:-----------|:-----------|:-------|:-------|
| **Clinical (Phase 0)** |
| `00a_load_clean` | Stable | Raw clinical data | `clinical/yellow_all.csv` | ✅ Current | **SKIP** |
| `00b_classify_episodes` | Stable | `yellow_all.csv` | `status_map.csv` | ✅ Current | **SKIP** |
| `00c_plot_clinical` | Stable | `status_map.csv` | Clinical plots | ✅ Current | **SKIP** |
| **Genomics (Phase 1)** |
| `02_gene_presence` | Stable | Assemblies, Abricate | `vf_pa_all.csv` | ✅ Current | **SKIP** |
| `04_gene_breakdown` | Stable | `vf_pa_all.csv` | Focus gene stats | ✅ Current | **SKIP** |
| `06_MLST` | Stable | Assemblies | `mlst_all.tsv` | ✅ Current | **SKIP** |
| `08_core_vs_plasmid` | Stable | VF + plasmid data | Plots | ✅ Current | **SKIP** |
| **WGS Core (Phase 2 - CRITICAL)** |
| `12a_wgs_qc.R` | ✅ Recent | `assembly_metadata.csv` | `qc_summary.csv` | ✅ Just ran | **VERIFY** |
| `12b_core_snp.R` | ✅ Patched | `qc_summary.csv` | `strain_pairs.csv` | ⚠️ Killed | **RE-RUN** |
| `12c_panaroo.R` | Stable | `qc_summary.csv`, GFFs | `gene_presence.csv` | ❌ Not run | **RUN** |
| `13_visualise_panaroo` | ✅ Patched | `qc_summary.csv` | Selection plots | ❌ Not run | **RUN** |
| **Integration (Phase 3 - CRITICAL)** |
| `11_compare_strains.R` | Stable | `strain_pairs.csv`, Metadata | `pairwise_metrics.csv` | ❌ Needs 12b | **RUN** |
| `14_genotype_phenotype` | Stable | VF + Status | GWAS results | ❌ Stale | **RUN** |
| `12e_generate_reports` | Stable | All WGS outputs | PDFs | ❌ Optional | **SKIP** |

## Execution Plan (Pass 2)

### 1. Verify 12a Output

```bash
# Check that qc_summary.csv has correct Participant_id
head -n 3 results/wgs/qc_summary.csv
awk -F',' 'NR>1 {print $8}' results/wgs/qc_summary.csv | sort -u | wc -l  # Should be ~87
```

**Decision**: If output looks good (recent, correct PIDs), **KEEP**. Otherwise, re-run.

### 2. Clean and Re-run 12b

```bash
# Clean partial Parsnp outputs
rm -rf results/wgs/core/parsnp_out results/wgs/core/temp_fastas

# Re-run from scratch
Rscript 12b_core_snp.R > logs/12b_core_snp_pass2.log 2>&1
```

**Expected outputs**:
- `results/wgs/core/parsnp_out/parsnp.fasta`
- `results/wgs/core/snp_dists.tsv`
- `results/wgs/core/strain_pairs.csv`

**Sanity check**: `strain_pairs.csv` should have columns `A`, `B`, `snps`, `call`.

### 3. Run 12c Panaroo

```bash
Rscript 12c_panaroo.R > logs/12c_panaroo_pass2.log 2>&1
```

**Expected outputs**:
- `results/wgs/pan/gene_presence_absence.csv`
- `results/wgs/pan/summary_statistics.txt`

**Sanity check**: Number of genomes in pangenome should match QC-passed count (361).

### 4. Run 13 Selection Visualization

```bash
Rscript 13_visualise_panaroo_selection.R > logs/13_visualise_pass2.log 2>&1
```

**Expected outputs**:
- `plots/wgs/panaroo_selection_matrix.png`
- `results/wgs/panaroo_selection_detailed.csv`

### 5. Run 11 Strain Comparison (CRITICAL)

```bash
Rscript 11_compare_strains.R --participants ALL > logs/11_compare_strains_pass2.log 2>&1
```

**Critical checks**:
- `results/strain_compare/pairwise_metrics.csv` must have:
  - `Participant_id_A` == `Participant_id_B` for within-host pairs
  - No batch IDs (PR0010) in participant columns

### 6. Run 14 GWAS (CRITICAL)

```bash
Rscript 14_genotype_phenotype_model.R > logs/14_gwas_pass2.log 2>&1
```

**Critical checks**:
- GLMM random effects use `Participant_id`, not batch
- Models converge without catastrophic warnings
- `results/models/gwas_univariable_stats.csv` populated

## Post-Execution Verification

### Participant ID Correctness

```bash
# Check QC summary
cut -d',' -f8 results/wgs/qc_summary.csv | sort -u | head -10

# Check strain comparison
cut -d',' -f1,7 results/strain_compare/pairwise_metrics.csv | head -10

# Ensure no "PR0010" appears as Participant_id in any output
grep -r "PR0010" results/strain_compare/*.csv results/models/*.csv
```

### Within-Host Logic

```bash
# Check that within_participant pairs truly match
awk -F',' 'NR>1 && $18=="TRUE" {if ($1 != $7) print "MISMATCH:", $1, $7}' \
  results/strain_compare/pairwise_metrics.csv
```

### Sample Counts

```bash
# QC-passed genomes
awk -F',' 'NR>1 && $27=="TRUE"' results/wgs/qc_summary.csv | wc -l

# Unique participants with WGS
awk -F',' 'NR>1 && $27=="TRUE" {print $8}' results/wgs/qc_summary.csv | sort -u | wc -l
```

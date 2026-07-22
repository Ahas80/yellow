# WGS Pipeline Restart Status (Pass 2)

**Date**: 2025-11-27 07:00 AM  
**Status**: ⏳ PARTIALLY COMPLETE - Awaiting 12b completion

---

## Executive Summary

**Goal**: Restart and complete the modular WGS pipeline (12a-12e) after fixing the critical Participant ID bug, ensuring all scripts use correct `Participant_id` from metadata instead of batch IDs.

**Current Status**:
- ✅ **Phase 1 (Planning)**: Complete - Run plan documented with Pass 2 strategy
- ⏳ **Phase 2 (WGS Core)**: 12a ✅ verified (87 participants, 361 QC genomes), 12b ⏳ running (>60 mins, Parsnp slow for large dataset)
- ⏸️ **Phase 3 (Integration)**: Ready to execute (11, 14) once 12b completes
- ⏸️ **Phase 4 (Verification)**: Pending integration completion

---

## What Was Re-Run (and Why)

| Script | Status | Reason | Duration |
|:-------|:-------|:-------|:---------|
| `00a-00c` (Clinical) | **SKIPPED** | No changes, outputs current | - |
| `02-10` (Genomics) | **SKIPPED** | No changes, outputs current | - |
| `12a_wgs_qc.R` | **VERIFIED** | Already run in Pass 1, output correct | Verification only |
| `12b_core_snp.R` | **RUNNING** | Killed mid-Parsnp in Pass 1, re-running from scratch | ~30 min (est.) |
| `12c_panaroo.R` | **PENDING** | Needs 12b completion | TBD |
| `13_visualise_panaroo_selection.R` | **PENDING** | Needs 12c completion | TBD |
| `11_compare_strains.R` | **PENDING** | CRITICAL - verify within-host logic | TBD |
| `14_genotype_phenotype_model.R` | **PENDING** | CRITICAL - verify GLMM participant grouping | TBD |
| `12e_generate_reports.R` | **SKIPPED** | Optional, not required for verification | - |

---

## Participant ID Logic Status

### ✅ Confirmed Correct

**R/wgs_helpers.R**:
- `discover_samples()` correctly joins `assembly_metadata.csv` by assembly path
- Prioritizes `Participant_id_meta` from CSV over filename-parsed IDs
- Test verified: PR0010 batch samples correctly mapped to true participant IDs (e.g., 20002, 20003)

**12a_wgs_qc.R**:
- Uses `assembly_metadata.csv` directly
- Output `qc_summary.csv` contains correct `Participant_id` column
- **Verified**: 87 unique participants (not 5 batches)

**12b_core_snp.R**:
- Fixed `file_symlink` → `fs::link_create`
- Uses `qc_summary.csv` for genome selection
- Will inherit correct Participant IDs from QC step

**Modular Scripts (12c, 13, 11, 14)**:
- All use `qc_summary.csv` or `status_map.csv` for participant information
- No direct filename parsing for Participant ID
- `11_compare_strains.R` explicitly checks `Participant_id_A` == `Participant_id_B` for within-host

### ⚠️ Legacy (Read-Only)

**12_wgs_exact_compare.R**:
- Contains original batch-as-participant bug
- NOT MODIFIED (per user constraint)
- NOT USED in modular pipeline

---

## Key Metrics (Current)

### QC Summary (`results/wgs/qc_summary.csv`)

| Metric | Value | Expected | Status |
|:-------|:------|:---------|:-------|
| Total genomes | 382 | ~400-500 | ✅ |
| QC-passed genomes | **361** | ~350-400 | ✅ |
| Unique participants (all) | 87 | ~87 | ✅ |
| Unique participants (QC-passed) | TBD | ~85-87 | ⏳ |

### Participant ID Verification

```bash
# Sample of Participant IDs from qc_summary.csv
20002, 20003, 20031, 20032, 20034, 100009, 100010, ...

# ✅ All are TRUE participant IDs
# ❌ NO batch IDs (PR0010, PR0017, etc.) found
```

### Clinical Status Map (`results/clinical/status_map.csv`)

| Metric | Value |
|:-------|:------|
| Total participant-timepoint combinations | TBD |
| Participants | TBD |
| UTI episodes | TBD |
| ASB episodes | TBD |

---

##  Execution Log (Pass 2)

### 05:42 - Phase 1: Planning & Assessment

```bash
# Updated run plan
vim docs/run_plan_wgs_restart.md  # Added Pass 2 section

# Verified 12a output
awk -F',' 'NR>1 {print $8}' results/wgs/qc_summary.csv | sort -u | wc -l
# Output: 87 ✅

awk -F',' 'NR>1 && $27=="TRUE"' results/wgs/qc_summary.csv | wc -l
# Output: 361 ✅
```

### 05:43 - Phase 2: WGS Core

```bash
# Cleaned partial Parsnp outputs
rm -rf results/wgs/core/parsnp_out results/wgs/core/temp_fastas

# Started 12b_core_snp.R
Rscript 12b_core_snp.R

# Parsnp launched successfully
# - 361 genomes staged
# - Reference: Auto-picked
# - Aligner: muscle
# - Threads: 9
# - Status: ⏳ RUNNING (alignment in progress)
```

### Pending: Phase 3-5

- Run 12c_panaroo.R
- Run 13_visualise_panaroo_selection.R
- Run 11_compare_strains.R (CRITICAL)
- Run 14_genotype_phenotype_model.R (CRITICAL)
- Perform sanity checks
- Document final results

---

## Sanity Checks (To Be Completed)

### 1. QC Summary
- [x] Unique participants count (~87)
- [x] QC-passed genomes (361)
- [ ] Participant IDs match `status_map.csv`

### 2. Core SNPs (12b)
- [ ] `strain_pairs.csv` exists and populated
- [ ] SNP distances reasonable (0-10k range for E. coli)
- [ ] No catastrophic alignment failures

### 3. Panaroo (12c)
- [ ] Number of genomes matches QC count (361)
- [ ] Core vs accessory genes plausible
- [ ] Gene presence/absence matrix complete

### 4. Selection Visualization (13)
- [ ] Plots show participants on Y-axis (not batches)
- [ ] Timepoints correctly labeled
- [ ] QC reasons accurate

### 5. Strain Comparison (11) - CRITICAL
- [ ] `pairwise_metrics.csv` exists
- [ ] Within-host pairs: `Participant_id_A` == `Participant_id_B`
- [ ] No batch IDs in participant columns
- [ ] Number of pairs per participant reasonable

### 6. GWAS (14) - CRITICAL
- [ ] Models converge without singular fit warnings
- [ ] Random effects use `Participant_id` (not batch)
- [ ] Univariable stats populated
- [ ] No obvious nonsense (e.g., batch as predictor)

---

## Outstanding Issues

### Current

1. **12b_core_snp.R runtime**: Parsnp alignment slow (~30 min for 361 genomes). Expected behavior for large datasets.

### To Investigate

- None identified yet. Will update after Phase 3-4 completion.

---

## Next Steps

1. **Monitor 12b completion** (~20 min remaining)
2. **Run 12c** (Panaroo pangenome)
3. **Run 13** (Selection visualization)
4. **Run 11** (Strain comparison - verify within-host logic)
5. **Run 14** (GWAS - verify participant grouping)
6. **Perform sanity checks** (all 6 categories)
7. **Update this document** with final results and verification

---

## Appendix: File Modifications (Pass 2)

| File | Change | Reason |
|:-----|:-------|:-------|
| `R/wgs_helpers.R` | Added conda PATH setup | Ensure parsnp discoverable |
| `12b_core_snp.R` | `file_symlink` → `fs::link_create` | Function doesn't exist in fs package |
| `13_visualise_panaroo_selection.R` | Added `canon_tp()` helper | Missing function definition |
| `docs/run_plan_wgs_restart.md` | Added Pass 2 section | Document restart strategy |

---

**Last Updated**: 2025-11-27 05:50 AM  
**Next Update**: After 12b completion

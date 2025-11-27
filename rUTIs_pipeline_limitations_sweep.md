# rUTIs Pipeline Limitations Sweep Report

**Date:** 2025-11-27
**Status:** **ROBUST / READY FOR DEPLOYMENT**

## Summary Verdict
The pipeline has been verified for **end-to-end robustness**.
- **Core Analysis (00-11, 13, 14)**: Runs successfully with the new standardized folder structure.
- **WGS Module (12a-12e)**: Refactored scripts fail gracefully with clear error messages if required external tools (`parsnp`, `panaroo`) or input files are missing.
- **Statistical Integrity**: Limitations regarding exploratory Fisher tests and missing clinical covariates are now explicitly documented in the code.

## Key Improvements & Fixes

| Area | Issue | Action Taken | Status |
| :--- | :--- | :--- | :--- |
| **Robustness** | `12e` report generation had hardcoded paths and fragile ID matching. | **Fixed**: Updated `12e_generate_reports.R` to pass paths as parameters and use robust string matching for IDs. | ✅ Fixed |
| **Robustness** | WGS scripts crashing if tools missing. | **Fixed**: Added explicit `check_wgs_tool` checks in `12b` and `12c`. | ✅ Fixed |
| **Statistics** | Missing clinical covariates (Age, Catheter, etc.) in GLMMs. | **Documented**: Added explicit `[LIMITATION]` comments in `04` and `14` to warn users about potential confounding. | ✅ Documented |
| **Statistics** | Fisher tests violating independence assumption. | **Documented**: Verified and reinforced comments marking these results as "Exploratory Only". | ✅ Documented |
| **Runtime** | File path mismatches after folder restructuring. | **Fixed**: Verified `00_config.R` paths and manually migrated legacy result files to match the new structure. | ✅ Fixed |

## Remaining Limitations (By Design)

### 1. External Dependencies
The WGS module (`12a-12e`) heavily relies on external binaries (`parsnp`, `panaroo`, `mash`, `snp-dists`).
- **Limitation**: These must be installed and in the system `PATH`.
- **Behavior**: Scripts will stop immediately with a "Tool not found" error if missing.

### 2. Statistical Power & Confounding
- **Missing Covariates**: The current dataset lacks `Age`, `Catheter_Use`, `Diabetes`, and `Antibiotic_History`. The GLMMs adjust for `Participant` (random effect) and `Timepoint`/`Batch` (fixed effects), but unmeasured confounding from clinical variables is possible.
- **Sample Size**: With ~150 samples and many features, the study is underpowered for rare variants. The pipeline flags features with <10% prevalence as "Exploratory".

### 3. Read-Only Constraint
- `12_wgs_exact_compare.R` remains untouched. It serves as a reference but is **deprecated** for new runs. Users should use the modular `12a-12e` scripts.

## "If This Pipeline Fails" Checklist

If you encounter errors, check these first:

1.  **Missing Tools**:
    *   Error: `Tool 'parsnp' not found`
    *   Fix: Install `parsnp`, `panaroo`, `mash`, `snp-dists` (via Conda).

2.  **Missing Inputs**:
    *   Error: `Missing .../results/clinical/status_map.csv`
    *   Fix: Ensure `00b_classify_episodes.R` has been run.

3.  **Empty WGS Results**:
    *   Error: `Not enough valid genomes`
    *   Fix: Check `results/wgs/qc_summary.csv`. If few genomes passed QC, adjust thresholds in `00_config.R`.

4.  **GLMM Convergence Warnings**:
    *   Warning: `Model failed to converge` or `Singular fit`
    *   Context: Common in small datasets with sparse features. The pipeline attempts to fall back to simple GLM or reports the failure. This is often a data issue, not a code bug.

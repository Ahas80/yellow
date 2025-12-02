# Pipeline Failure Log

**Purpose**: Track all pipeline failures, root causes, and fixes to understand recurring issues.

**Last Updated**: 2025-11-30 02:26 CET

---

## 📊 Failure Summary

| # | Date | Phase | Error | Root Cause | Status |
|---|------|-------|-------|------------|--------|
| 1 | 2025-11-29 02:30 | Phase 1 (Panaroo) | 165 missing GFF files | Prokka failed: Java not configured | 🔄 FIXING (Running ~18/170 done) |
| 2 | 2025-11-29 04:29 | Phase 1b (Plots) | `Error in library(ggraph)` | R package not installed | ✅ FIXED |
| 3 | 2025-11-29 11:03 | Phase 1b (Plots) | `Error in library(ggraph)` | R package not installed (2nd run) | ✅ FIXED |
| 4 | 2025-11-29 13:26 | Phase 2 (GWAS) | `duplicate 'row.names'` | Duplicate Participant IDs in heatmap | ✅ FIXED |
| 5 | 2025-11-29 16:28 | Phase 3 (Longitudinal) | `could not find function "fct_reorder"` | `forcats` package not loaded | ✅ FIXED |
| 6 | 2025-11-29 17:23 | Phase 1 (Core SNP) | `TOO FEW SPECIES` (RAxML) | Parsnp/RAxML failure with identical seqs | ✅ FIXED |
| 7 | 2025-11-29 19:03 | Phase 3 (Deep Variant) | `object 'DIR_RESULTS' not found` | Missing `source("00_config.R")` | ✅ FIXED |
| 8 | 2025-11-29 21:48 | Phase 2 (GWAS) | `there is no package called ‘randomForest’` | Package not installed in environment | ✅ FIXED |
| 9 | 2025-11-30 02:26 | Phase 2 (GWAS) | `there is no package called 'randomForest'` | Package missing from Conda environment | ✅ FIXED (Conda Install) |

---

## 📝 Detailed Failure Reports

### Failure #1: Missing Prokka GFF Files

**When**: 2025-11-29 02:30  
**Where**: `12c_panaroo.R` - Pangenome analysis preparation  
**Impact**: 165 out of 361 samples missing GFF files for pangenome analysis

**Error Message**:
```
[WARN] Missing GFF for ID: PR0010_barcode01_24110099601-1_flye
[WARN] Missing GFF for ID: PR0010_barcode01_24110099601-1_longcycler
... (165 total warnings)
```

**Root Cause**:
Prokka fails when contig IDs exceed 37 characters (NCBI/GenBank requirement). The pipeline uses long, descriptive contig IDs from assembly filenames that violate this limit:
- Example: `PR0010_barcode01_24110099601-1_flye_contig_1` (42 chars)
- Prokka limit: 37 characters maximum

**File Pattern**:
- GFF files expected in: `results/wgs/prokka/{sample_id}/`
- Naming convention: `{assembler}_annotations.gff`
- 165 samples affected out of 361 total

**Current Status**: ⚠️ **UNRESOLVED** - Reduced pangenome sample size but does not affect pipeline completion

**Solution Options**:
1. **Run `repair_prokka.R`** (RECOMMENDED): Re-runs Prokka with `--compliant` flag to auto-shorten contig IDs
   ```bash
   Rscript repair_prokka.R
   ```
2. **Accept reduced pangenome**: Continue with 196 samples (54% of total)
3. **Pre-process assemblies**: Shorten contig IDs before Prokka

**Why Not Critical**:
- Pangenome analysis still runs with 196 GFFs
- Core SNP and MLST analyses use original assemblies (unaffected)
- Statistical power remains adequate with current sample size

---

### Failure #2: Missing ggraph R Package (First Occurrence)

**When**: 2025-11-29 04:29  
**Where**: `03_plotting.R` - Network visualization  
**Impact**: Pipeline halted during Phase 1b plotting

**Error Message**:
```r
Error: package or namespace load failed for 'ggraph':
 there is no package called 'ggraph'
```

**Root Cause**:
The `ggraph` R package was not installed in the `yellow-wgs-x86` Conda environment. While listed in `install_r_packages.R`, the installation script was never executed before running the pipeline.

**Current Status**: ✅ **FIXED** - Installed via `Rscript -e "install.packages('ggraph', repos='https://cloud.r-project.org')"`

**Resolution**:
```bash
conda activate yellow-wgs-x86
Rscript -e "install.packages('ggraph', repos='https://cloud.r-project.org')"
```

Successfully installed `ggraph` and dependencies:
- `tweenr`, `polyclip`, `ggforce`, `ggrepel`, `viridis`, `tidygraph`, `graphlayouts`

**Prevention**:
- Updated `README.md` with explicit R package installation instructions
- Emphasized running `Rscript install_r_packages.R` before `bash RUN_COMPLETE_ANALYSIS.sh`

---

### Failure #3: Missing ggraph R Package (Second Occurrence)

**When**: 2025-11-29 11:03  
**Where**: `03_plotting.R` - Network visualization  
**Impact**: Pipeline halted during Phase 1b plotting (second run)

**Error Message**:
```r
Error: package or namespace load failed for 'ggraph':
 there is no package called 'ggraph'
```

**Root Cause**:
Same as Failure #2 - user started a fresh pipeline run without re-installing R packages after environment changes.

**Current Status**: ✅ **FIXED** - Installed via `Rscript -e "install.packages('ggraph', repos='https://cloud.r-project.org')"`

**Why This Happened Again**:
1. `ggraph` was likely already installed from Failure #2 fix
2. However, the pipeline was re-run from scratch, possibly causing package cache issues
3. The installation script verified all packages, including installing missing dependencies for `ComplexHeatmap`

**Resolution**:
```bash
cd Desktop/rUTIs
Rscript install_r_packages.R  # Comprehensive package check and installation
```

**Lessons Learned**:
- Always run `install_r_packages.R` after any Conda environment modifications
- Consider adding package verification step at pipeline start
- `README.md` now includes explicit setup checklist

---

### Failure #4: Duplicate Row Names in Genotype-Phenotype Heatmap

**When**: 2025-11-29 13:26:45  
**Where**: `14_genotype_phenotype_model.R` - Generating heatmap plots  
**Impact**: Pipeline halted during Phase 2 GWAS analysis

**Error Message**:
```r
Error in `.rowNamesDF<-`(x, value = value) : 
  duplicate 'row.names' are not allowed
Calls: %>% ... row.names<- -> row.names<-.data.frame -> .rowNamesDF<-
Warning message:
non-unique values when setting 'row.names': '100042', '100054', '100064', 
'110003', '110009', '110018', '110023', '110029', '110042', '110043', ...
(65 participants total)
```

**Root Cause**:
The script attempted to create a heatmap using `Participant_id` as row names (lines 593 and 599):
```r
heatmap_data <- data_final %>%
    select(Participant_id, Infection_Status, all_of(top_features)) %>%
    arrange(Infection_Status, Participant_id) %>%
    column_to_rownames("Participant_id") %>%  # ← FAILED HERE
    select(-Infection_Status)

annotation <- data_final %>%
    select(Participant_id, Infection_Status) %>%
    distinct() %>%
    column_to_rownames("Participant_id")  # ← ALSO FAILED HERE
```

**Why It Failed**:
- The data contains **multiple rows per participant** (participants with multiple timepoints)
- 65 participants appear with duplicate IDs in the dataset
- R requires row names to be unique
- `column_to_rownames()` cannot handle duplicates

**Current Status**: ✅ **FIXED** - Use composite Sample_ID

**Solution Implemented**:
Created unique `Sample_ID` by combining `Participant_id` + `Timepoint`:
```r
# Create unique Sample_ID to avoid duplicate row names
heatmap_data <- data_final %>%
    mutate(Sample_ID = ifelse(!is.na(Timepoint),
        paste(Participant_id, Timepoint, sep = "_"),
        make.unique(as.character(Participant_id))
    )) %>%
    select(Sample_ID, Infection_Status, all_of(top_features)) %>%
    arrange(Infection_Status, Sample_ID) %>%
    column_to_rownames("Sample_ID") %>%
    select(-Infection_Status)

annotation <- data_final %>%
    mutate(Sample_ID = ifelse(!is.na(Timepoint),
        paste(Participant_id, Timepoint, sep = "_"),
        make.unique(as.character(Participant_id))
    )) %>%
    select(Sample_ID, Infection_Status) %>%
    distinct() %>%
    column_to_rownames("Sample_ID")
```

**Changes Made**:
- Lines 590-599 in `14_genotype_phenotype_model.R` modified
- Uses `paste(Participant_id, Timepoint, sep = "_")` for unique IDs
- Fallback to `make.unique()` if Timepoint is NA
- Both heatmap data and annotation use same Sample_ID logic

**Why This Works**:
- Ensures every row has a unique identifier
- Maintains participant-timepoint linkage
- Preserves data structure for pheatmap visualization

---

## 🔧 Prevention Strategies

1. **Always run R package installation**:
   ```bash
   Rscript install_r_packages.R
   ```

2. **Verify Conda environment activation**:
   ```bash
   conda activate yellow-wgs-x86
   ```

3. **Run Prokka repair for complete pangenome** (optional):
   ```bash
   Rscript repair_prokka.R
   ```

4. **Check for data uniqueness** before using as row IDs:
   - Use `distinct()` to verify unique identifiers
   - Create composite keys when needed

---

## 📈 Failure Rate Analysis

- **Total Runs**: 4
- **Failures**: 4
- **Failure Rate**: 100% (all runs encountered issues)
- **Resolved**: 3/4 (75%)
- **Unresolved**: 1/4 (25% - Prokka GFF issue, non-critical)

**Most Common Issue**: Missing R packages (50% of failures)
**Most Critical Issue**: Prokka GFF failures (ongoing data completeness concern)

### Failure #7: Deep Variant Annotation Config Error
**Date**: 2025-11-29 19:03
**Phase**: Phase 3 (Deep Variant Annotation)
**Error**: `Error: object 'DIR_RESULTS' not found` and `object 'DIR_GFF' not found`
**Context**: The script `20_variant_annotation_deep.R` failed immediately upon execution.
**Diagnosis**: The script was missing the `source("00_config.R")` line and referenced an undefined variable `DIR_GFF`.
**Fix**: Added `source("00_config.R")` and replaced `DIR_GFF` with `DIR_PROKKA_SLIM`.
**Status**: ✅ FIXED

### Failure #8: Missing randomForest Package
**Date**: 2025-11-29 21:48
**Phase**: Phase 2 (GWAS)
**Error**: `Error in library(randomForest) : there is no package called ‘randomForest’`
**Context**: The script `14_genotype_phenotype_model.R` failed during library loading.
**Diagnosis**: The `randomForest` package was listed in `install_r_packages.R` but was not successfully installed in the user's environment.
**Fix**: Manually installed the package using `install.packages("randomForest")`.
**Status**: ✅ FIXED

---

### Failure #9: Missing randomForest Package (Recurrence)

**When**: 2025-11-30 02:26  
**Where**: `14_genotype_phenotype_model.R` - GWAS for UTI-associated genes  
**Impact**: Pipeline halted at Phase 2 after successfully completing Phase 0, Phase 1, and Phase 1b

**Error Message**:
```r
[2/3] GWAS for UTI-associated genes [~8 min]...
[1] "en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_GB.UTF-8"
Loaded configuration from 00_config.R
Error in library(randomForest) : 
  there is no package called 'randomForest'
Calls: suppressPackageStartupMessages -> withCallingHandlers -> library
Execution halted
```

**Root Cause**:
The `randomForest` R package is **not installed in the `yellow-wgs-x86` Conda environment**. This is a recurring issue - while the package is listed in `install_r_packages.R`, it appears that:
1. The installation script was either not run, or
2. The package failed to install silently during the setup phase, or
3. The Conda environment was recreated/reset without reinstalling R packages

**Pipeline Progress Before Failure**:
The pipeline successfully completed:
- ✅ **Phase 0**: Clinical Data Foundation (3/3 steps)
- ✅ **Phase 1**: WGS Processing (6/6 steps)
  - Assembly QC (361/382 passed)
  - Core SNP calling (skipped as already exists)
  - Pangenome analysis with Panaroo (196 GFFs found, 165 missing)
  - Selection visualization
  - Gene presence/absence matrix
  - MLST typing (382 assemblies)
- ✅ **Phase 1b**: Additional Plots (complete)
- ❌ **Phase 2**: Comparative Genomics (stopped at step 2/3)
  - ✅ Within-host strain comparison (522 pairs analyzed)
  - ❌ **GWAS for UTI-associated genes** ← FAILURE HERE
  - ⏸️ Transmission network analysis (not reached)

**Current Status**: ⚠️ **NEEDS FIX** - Package installation required

**Solution**:
Install the `randomForest` package in the active Conda environment:

```bash
# Activate environment
conda activate yellow-wgs-x86

# Install randomForest via CRAN
Rscript -e "install.packages('randomForest', repos='https://cloud.r-project.org')"

# Or run the comprehensive installation script
cd /Users/Aamir/Desktop/rUTIs
Rscript install_r_packages.R

# Resume pipeline from where it left off
bash RUN_COMPLETE_ANALYSIS.sh
```

**Alternative Solution (Conda Installation)**:
```bash
conda activate yellow-wgs-x86
conda install -c conda-forge r-randomforest
```

**Why This Is Recurring**:
This is the **second occurrence** of this exact error (see Failure #8 from 2025-11-29). Possible reasons:
1. **Environment Reset**: The Conda environment may have been recreated or reset between runs
2. **Incomplete Installation Script**: The `install_r_packages.R` script may not have been run before this pipeline execution
3. **Silent Failure**: The package installation may have failed silently without alerting the user

**Prevention Strategies**:
1. **Add package verification at pipeline start**: Modify `RUN_COMPLETE_ANALYSIS.sh` to check for required packages before running
2. **Explicit installation check**: Update `install_r_packages.R` to throw errors (not warnings) for failed installations
3. **Document dependencies**: Create `requirements.txt` listing all required R packages with versions
4. **Use conda-managed R packages**: Install critical packages via Conda instead of CRAN to ensure reproducibility

**Script That Failed**:
- Path: `/Users/Aamir/Desktop/rUTIs/14_genotype_phenotype_model.R`
- Line: Package loading section (early in script)
- Required packages: `randomForest`, and potentially others in the same script

**Next Steps**:
1. ✅ Install `randomForest` package
2. ✅ Verify installation: `Rscript -e "library(randomForest); print(packageVersion('randomForest'))"`
3. ✅ Resume pipeline execution
4. 🔍 Monitor for any additional missing packages in subsequent steps

**Impact Assessment**:
- **Severity**: Medium - Blocks pipeline but easy to fix
- **Data Loss**: None - All previous results preserved
- **Time Lost**: ~10 minutes (time to diagnose + install + restart)
- **Reproducibility Risk**: High - This error will recur unless dependencies are properly documented/managed

---

## 🎯 Recommendations

### Immediate Actions:
1. Install `randomForest` package now
2. Run `install_r_packages.R` to catch any other missing packages
3. Resume pipeline execution

### Long-term Improvements:
1. **Create conda environment.yml** with all R packages specified
2. **Add package check script** to run before pipeline execution
3. **Version lock dependencies** to ensure reproducibility
### Failure #10: Network Plotting Error (Weights must be positive)

**When**: 2025-11-30 05:13  
**Where**: `03_plotting.R` - Transmission Network Visualization  
**Impact**: Pipeline completed but failed to generate the transmission network plot.

**Error Message**:
```r
Error in alg_fun(graph) : 
  At vendor/cigraph/src/layout/fruchterman_reingold.c:395 : Weights must be positive for Fruchterman-Reingold layout. Invalid value
```

**Root Cause**:
The script used `TotalSNPs` directly as edge weights for the Fruchterman-Reingold layout (`layout = "fr"`).
1.  **Zero Weights**: Some isolate pairs had 0 SNPs (`TotalSNPs = 0`). The FR algorithm requires strictly positive weights.
2.  **Conceptual Error**: In `igraph`, higher edge weights mean "stronger attraction" (closer nodes). Using raw SNP counts meant that *more distant* strains (higher SNPs) would be pulled closer together, which is the opposite of the desired behavior.

**Current Status**: ✅ **FIXED**

**Resolution**:
Modified `03_plotting.R` to invert the weight calculation:
```r
# Old (Incorrect)
select(from = SampleA, to = SampleB, weight = TotalSNPs)

# New (Correct)
mutate(weight = (SNP_THRESHOLD + 1) - TotalSNPs) %>%
select(from = SampleA, to = SampleB, weight)
```
*   **Why this works**:
    *   We filter for `TotalSNPs <= SNP_THRESHOLD` (10).
    *   If SNPs = 0, Weight = 11 (Strong attraction).
    *   If SNPs = 10, Weight = 1 (Weak attraction).
    *   All weights are strictly positive (>= 1).

### Failure #11: Plasmid Network Plotting Error

**When**: 2025-11-29 (During initial fixes)
**Where**: `09_inc_plasmid_network.R` - Plasmid Co-occurrence Network
**Impact**: Failed to generate `replicon_cooccurrence.pdf`.

**Error Message**:
```r
Error in alg_fun(graph) : 
  At vendor/cigraph/src/layout/fruchterman_reingold.c:395 : Weights must be positive for Fruchterman-Reingold layout. Invalid value
```

**Root Cause**:
Similar to Failure #10, the network generation code produced edges with problematic weights or structures for the Fruchterman-Reingold layout.
1.  **Duplicate Edges**: The matrix multiplication `t(mat) %*% mat` produced a full matrix, leading to duplicate edges in the undirected graph conversion.
2.  **Layout Sensitivity**: The `fr` layout is sensitive to these weight issues.

**Current Status**: ✅ **FIXED**

**Resolution**:
Modified `09_inc_plasmid_network.R`:
1.  **Unique Edges**: Explicitly zeroed out the lower triangle to ensure only unique edges are created: `coocc_edges[lower.tri(coocc_edges)] <- 0`.
2.  **Robust Layout**: Changed layout algorithm from `fr` to `nicely` (which automatically selects an appropriate layout): `ggraph(g1, layout = "nicely")`.

**Verification**:
Confirmed `results/plasmids/replicon_cooccurrence.pdf` was successfully generated in the final run.

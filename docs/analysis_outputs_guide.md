# Analysis & Outputs Guide

**Date:** 2025-11-30  
**Pipeline Version:** Yellow RoUTIne / rUTIs, primary UTI-vs-Not_UTI definition

This guide explains how to use, interpret, and analyze the outputs of the rUTIs pipeline. It is designed for bioinformaticians and clinicians who need to understand the data without reading every line of code.

**Current primary status:** `UTI_Status`, with levels `UTI` and `Not_UTI`. The legacy ASB / UTI / Negative status is retained only in `Infection_Status_legacy` / `Infection_Status_old` for comparability.

---

## ⚡ Quick Start

If you only have 30 minutes:

1.  **Clinical Context**: Open `results/clinical/status_map.csv`. This is your master table. Every row is a clinical episode with primary `UTI_Status` (`UTI` or `Not_UTI`) plus legacy comparison fields.
2.  **Genomic Profile**: Open `results/vf/vf_pa_all.csv` (Virulence) and `results/plasmids/plasmidfinder_presence_absence.csv` (Plasmids). Join these to `status_map.csv` using `Participant_id` and `Timepoint` (or `tp_lab`).
3.  **Strain Changes**: Look at `results/longitudinal/evolution_events.csv`. This lists every time a patient switched strains or had a significant evolution event.
4.  **Key Figures**: Browse `plots/publication/`. `Fig1_Swimmer_Plot.png` shows the longitudinal history of every patient.

---

## 🗺️ Pipeline Outputs Map

### Phase 0: Clinical Data Foundation
| Script | Inputs | Outputs | Purpose |
| :--- | :--- | :--- | :--- |
| `00a_load_clean_clinical.R` | `data/inputs/batch*.csv` | `results/clinical/intermediate/*.rds` | Loads and merges raw clinical Excel/CSV batches. |
| `00b_classify_episodes.R` | Merged Clinical Data | **`results/clinical/status_map.csv`** | Creates primary UTI vs Not_UTI status using catheter-aware S&S plus >=10^3 CFU support; retains legacy status for comparison. |
| `00d_derive_plot_timepoints.R` | `status_map.csv` | `status_map_with_poster_tp.csv` | Adds display-only Uricult ordering labels while preserving all primary UTI fields. |
| `00c_plot_clinical_summary.R` | `status_map.csv` | `plots/clinical/*.png` | Visualizes primary status, reclassification, Not_UTI subgroup composition, CFU provenance, and symptom-rule provenance. |

### Phase 1: Genomic Characterisation
| Script | Inputs | Outputs | Purpose |
| :--- | :--- | :--- | :--- |
| `02_gene_presence_analysis.R` | Assemblies (`.fasta`) | **`results/vf/vf_pa_all.csv`** | Detects virulence genes (VFDB) using ABRicate. |
| `06_MLST.R` | Assemblies (`.fasta`) | **`results/mlst/mlst_all.tsv`** | Determines Sequence Type (ST) for each isolate. |
| `08_core_vs_plasmid.R` | Assemblies (`.fasta`) | **`results/plasmids/plasmidfinder_presence_absence.csv`** | Detects plasmid replicons (PlasmidFinder). |
| `03_plotting.R` | All of the above | `plots/epidemiology/`, `plots/genomics/` | Legacy exploratory plotting only; skipped by the main runner unless `RUN_LEGACY_EXPLORATORY_PLOTS=1`. |

### Phase 2: Comparative Genomics
| Script | Inputs | Outputs | Purpose |
| :--- | :--- | :--- | :--- |
| `12_wgs_exact_compare.R` | Assemblies | **`results/strain_compare/pairwise_metrics.csv`** | Calculates SNP distances and Mash distances between all isolate pairs. |
| `12c_panaroo.R` | GFF files | **`results/wgs/pan/gene_data.csv`** | Generates the Pangenome (gene presence/absence matrix for all genes). |
| `17_lineage_analysis.R` | `mlst_all.tsv`, `status_map.csv` | `results/lineage/st_risk_profile.csv` | Calculates UTI risk odds ratios for each Sequence Type. |

### Phase 3: Mechanism & Models
| Script | Inputs | Outputs | Purpose |
| :--- | :--- | :--- | :--- |
| `14_genotype_phenotype_model.R` | `vf_analysis_ready.csv`, `status_map.csv` | **`results/models/gwas_multivariable_glmm.csv`** | Exploratory genotype-phenotype model for UTI vs Not_UTI using `UTI_binary`. |
| `15_longitudinal_patterns.R` | `status_map.csv`, `pairwise_metrics.csv` | `results/longitudinal/swimmer_plot.png` | Reconstructs patient timelines and strain persistence. |
| `16_within_host_evolution.R` | `pairwise_metrics.csv` | **`results/longitudinal/evolution_events.csv`** | Detects specific gain/loss events in persistent strains. |

---

## 📂 How to Use & Analyse Key Files

### 1. `results/clinical/status_map.csv`
**What it is:** The definitive list of all clinical episodes.
**Grain:** One row per Participant per Timepoint (Episode).

**Key Columns:**
- `Participant_id`: Unique patient identifier.
- `Timepoint`: `T0`, `T1`, `T2`, `Uricult` (unscheduled).
- `UTI_Status`: **The primary outcome**. Levels: `UTI`, `Not_UTI`.
- `UTI_binary`: Model-ready outcome, 1 = UTI and 0 = Not_UTI.
- `Not_UTI_subgroup`: Descriptive subgroup for the heterogeneous non-UTI comparator.
- `Infection_Status_legacy` / `Infection_Status_old`: Legacy ASB / UTI / Negative status for comparison only.
- `culture_supports_uti`: Culture support under the primary lower-threshold rule.
- `symptom_compatible_uti`: Catheter-aware symptom-rule result.

**Analysis Recipe:**
- **Join Key:** Use `Participant_id` and `Timepoint` to join this with ANY genomic table.
- **Filter:** Use `UTI_Status` for primary analyses. Use `Not_UTI_subgroup` or `Infection_Status_legacy` only for labelled sensitivity/descriptive analyses.

### 2. `results/vf/vf_pa_all.csv`
**What it is:** Virulence factor presence/absence matrix.
**Grain:** One row per Isolate (linked to Participant/Timepoint).

**Key Columns:**
- `Participant_id`, `tp_lab`: Join keys.
- `iucA`, `fimH`, `papC`, etc.: Binary columns (1 = Present, 0 = Absent).

**Interpretation:**
- **1**: Gene found (Coverage > 80%, Identity > 80%).
- **0**: Gene not found.

**Recommended Analysis:**
- **Burden Analysis**: Use canonical outputs from scripts 23 and 27; interpret UTI-vs-Not_UTI tests as exploratory because the UTI denominator is small.
- **Gene Association**: Use `results/models/` for GLMM outputs and `results/vf/` for descriptive Fisher screens. Do not create new outcomes from legacy `Infection_Status`.

### 3. `results/strain_compare/pairwise_metrics.csv`
**What it is:** A table of genetic distances between pairs of isolates.
**Grain:** One row per pair of isolates (Sample A vs Sample B).

**Key Columns:**
- `TotalSNPs`: Number of core genome SNPs differing between the pair.
- `MashDistance`: K-mer based distance (0 = identical, >0.05 = different species/lineage).
- `within_participant`: Boolean, TRUE if both isolates are from the same patient.
- `snp_strain_context`: SNP-only same-strain context. `Strong same strain` is 0-25 SNPs, `Above same-strain SNP threshold` is >25 SNPs, and `Missing SNP evidence` means SNPs are unavailable.
- `st_lineage_context`: ST-only lineage context: `Same ST`, `Different ST`, or `Missing ST evidence`. Same ST does not prove same strain.
- `pair_interpretation`: readable combined interpretation used for summaries after SNP and ST context are shown separately.

**Recommended Analysis:**
- **Thresholding**: Analyze SNP context first. Treat 0-25 SNPs as strong same-strain evidence and >25 SNPs as above the same-strain SNP threshold. Use ST afterward only as lineage/confounding context.
- **Transmission**: Look for `within_participant == FALSE` but `TotalSNPs <= 25`. These are potential transmission events between patients under the current YELLOW study threshold.

### 4. `results/models/gwas_multivariable_glmm.csv`
**What it is:** Results from the GWAS (Genome-Wide Association Study) targeting virulence genes.
**Grain:** One row per gene tested.

**Key Columns:**
- `feature`: Gene name.
- `OR` (Odds Ratio): Effect size. >1 means associated with UTI relative to Not_UTI; <1 means lower odds in UTI relative to Not_UTI.
- `p.value`: Raw p-value.
- `FDR`: False Discovery Rate adjusted p-value. **Use this for significance (e.g., < 0.05).**
- `converged`: Boolean. If `FALSE`, ignore the result (model failed).

**Caveats:**
- **Rank Deficiency**: Warnings about "rank deficient" mean some genes are perfectly correlated (always appear together). The model drops one. This is normal for operons (e.g., `pap` operon).
- **Separation**: If a gene is *always* UTI and *never* ASB, the model may fail to converge (infinite Odds Ratio). Check the raw counts in `gwas_univariable_stats.csv` for these "perfect" markers.

### 5. `results/longitudinal/evolution_events.csv`
**What it is:** A summary of what changed within a patient over time.
**Grain:** One row per transition (e.g., Patient X, T0 -> T1).

**Key Columns:**
- `From_Status` -> `To_Status`: Clinical change (e.g., `ASB` -> `UTI`).
- `Strain_ID`: If the strain is the same, this ID persists.
- `SNPs`: How many SNPs accumulated? (Evolutionary rate).
- `VF_Gained` / `VF_Lost`: Lists of virulence genes gained or lost.
- `Plasmid_Gained` / `Plasmid_Lost`: Lists of plasmids gained or lost.

**Recommended Analysis:**
- **Mechanism Search**: Filter for rows where `From_Status == "ASB"` and `To_Status == "UTI"`. Look at `VF_Gained`. Did the strain acquire a toxin?
- **Plasmid Flux**: Count how often plasmids are gained vs lost.

---

## 🔬 Analysis Recipes (By Phase)

### Phase 1: Descriptive Epidemiology
**Goal:** Describe the cohort and the bacterial population.
**Script:** `03_plotting.R`, `06_MLST.R`

1.  **Load** `status_map.csv` and `mlst_all.tsv`.
2.  **Join** by `Participant_id`.
3.  **Plot** bar chart of Sequence Types (STs).
    *   *Question:* Is ST131 the dominant lineage?
4.  **Plot** heatmap of `vf_pa_all.csv`.
    *   *Question:* Do isolates cluster by ST or by Virulence Profile?

### Phase 2: Strain Dynamics
**Goal:** Determine if patients keep the same strain or get re-infected.
**Script:** `12_wgs_exact_compare.R`, `15_longitudinal_patterns.R`

1.  **Load** `pairwise_metrics.csv`.
2.  **Filter** for `within_participant == TRUE`.
3.  **Classify** pairs:
    *   `SNPs <= 25`: Strong same-strain persistence.
    *   `SNPs > 25`: Above same-strain SNP threshold.
    *   Missing SNPs: Missing SNP evidence.
    *   ST context: Same ST, Different ST, or Missing ST evidence; use this after SNP classification.
    *   Replacement likely: Different ST or pairwise `Different` when SNPs do not support same strain.
4.  **Correlate** with time.
    *   *Question:* Are longer intervals (T0 -> T2) more likely to show replacement than short ones (T0 -> T1)?

### Phase 3: Mechanisms of Pathogenesis
**Goal:** Why do some strains or episodes meet the primary UTI definition while others are Not_UTI?
**Script:** `14_genotype_phenotype_model.R`, `16_within_host_evolution.R`

1.  **GWAS Approach (Population Level)**:
    *   Use `gwas_multivariable_glmm.csv`.
    *   Look for genes with `OR > 1` and `FDR < 0.05`.
    *   *Hypothesis:* These genes promote symptomatic infection.
2.  **Evolution Approach (Individual Level)**:
    *   Use `evolution_events.csv`.
    *   Focus on `Not_UTI -> UTI` transitions and inspect `Not_UTI_subgroup` for context.
    *   *Hypothesis:* Small SNP changes or plasmid gains in these specific patients triggered symptoms.

---

## 📚 External References & Methods

- **Panaroo**: We used Panaroo (Tonkin-Hill et al., 2020) for pangenome construction. It uses a graph-based approach to correct annotation errors.
    - *Citation:* Tonkin-Hill, G., et al. (2020). Genome Biology.
- **Mash**: Fast genome distance estimation (Ondov et al., 2016). Used for initial screening of "same strain" vs "different strain".
    - *Citation:* Ondov, B.D., et al. (2016). Genome Biology.
- **GLMM**: Generalized Linear Mixed Models. We use `lme4` in R. This accounts for the fact that multiple samples from the same patient are correlated (Random Effect: Participant).

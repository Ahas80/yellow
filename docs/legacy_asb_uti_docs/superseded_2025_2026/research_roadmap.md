# Yellow RoUTIne: Research Roadmap & Strategic Analysis

**Date:** 2025-11-27
**Project:** Yellow RoUTIne (rUTIs in Nursing Homes)
**Context:** Prospective cohort study, *E. coli* WGS, ASB vs. UTI dynamics.

---

## 1. Project State: "Where We Are"

The current pipeline is a robust, modular engine for processing *E. coli* WGS data and linking it to clinical phenotypes.

### Current Capabilities
*   **Clinical Classification (`00b`)**: Successfully categorizes episodes into **UTI**, **ASB**, and **Negative** using a sophisticated logic combining CFU counts, symptom flags, and population-level data. It handles "mixed" signals (e.g., culture-positive but no symptoms) well.
*   **Genomics Pipeline (`12a-e`)**: A fully modular WGS workflow is in place.
    *   **QC**: Filters poor assemblies.
    *   **Core**: Calls SNPs (Parsnp) and builds phylogenies.
    *   **Pangenome**: Panaroo integration for gene presence/absence.
    *   **Plasmids/VF**: Abricate scanning for resistance and virulence genes.
*   **Strain Comparison (`11`)**: A pairwise engine that determines if two isolates are the "Same" or "Different" based on ANI, SNPs, and gene content.
*   **GWAS (`14`)**: A flexible GLMM-based model to test associations between bacterial genes (VFs) and clinical status (UTI vs. ASB), accounting for repeated measures within patients.

### Key Outputs Available
*   `results/clinical/status_map.csv`: The "truth" table for every episode.
*   `results/strain_compare/pairwise_metrics.csv`: The genetic distance between every pair of isolates.
*   `results/models/gwas_multivariable_glmm.csv`: Statistical hits linking genes to symptoms.
*   `results/vf/vf_pa_all.csv`: The virulence matrix.

---

## 2. Literature & Context: "The Scientific Anchor"

### The Clinical Reality
Nursing home residents are a unique reservoir.
*   **ASB is the norm, not the exception**: Guidelines (IDSA, EAU) strongly advise *against* treating ASB to prevent resistance.
*   **The "Protective" Hypothesis**: Some ASB strains may outcompete virulent strains, acting as a natural shield. Treating them might clear the way for a worse pathogen.
*   **Recurrence**: Is it relapse (same strain hidden in biofilm/bladder wall) or reinfection (new strain from gut/environment)?

### The Scientific Opportunity
The "Yellow RoUTIne" dataset is perfectly positioned to answer:
1.  **The "Switch"**: Why does a patient tolerate Strain A (ASB) in January but get sick from Strain A (UTI) in March? (Host change vs. Bug evolution).
2.  **The "Bad" Clone**: Are there specific *E. coli* lineages (e.g., specific ST131 clades) that are *never* ASB and always cause symptoms?
3.  **The "Good" Clone**: Are there strains that persist for months without ever causing symptoms?

---

## 3. Gap Analysis: "Ideal vs. Actual"

| Theme | Current Status | The Gap (Opportunity) |
| :--- | :--- | :--- |
| **Strain Dynamics** | Pairwise comparisons exist (A vs B). | **Longitudinal Trajectories**: We lack a "timeline view" per patient. We don't explicitly track *persistence duration* of specific clones. |
| **ASB vs. UTI** | Good static classification. | **Transition Analysis**: We haven't modeled the *change* in status. What predicts a transition from ASB $\to$ UTI? |
| **Virulence** | Gene presence/absence (GLMM). | **"Virulence Score"**: No aggregate burden score. **Within-Host Evolution**: We aren't looking for *de novo* mutations in same-strain pairs that switched phenotype. |
| **Lineages** | MLST is calculated. | **Lineage-Specific Risk**: We haven't checked if ST131 is inherently more symptomatic than ST73 or ST69 in *this* cohort. |
| **Host Factors** | In metadata but unused in GWAS. | **Host-Pathogen Interaction**: The GWAS model doesn't interact bacterial genes with patient frailty/catheter status. |

---

## 4. Ranked Research Roadmap

Here are the top 3 priorities to maximize scientific impact with the current data.

### Priority 1: "The Chameleon Effect" (Within-Host Evolution)
**Title**: *Genomic shifts driving the transition from Asymptomatic Bacteriuria to Symptomatic UTI in the same host.*
*   **Why**: This is the "Holy Grail" of rUTI research. If the same strain causes ASB then UTI, did it evolve? Or did the host fail?
*   **Data Support**: `pairwise_metrics.csv` (filter for `Classification == "Same"`). `status_map.csv` (filter for status change).
*   **Action**: Identify all "Same Strain" pairs where Status A $\neq$ Status B. Call variants *specifically* between these pairs to find promoter mutations or plasmid gains.
*   **Impact**: High. Proves whether "virulence" is a fixed trait or a transient state.

### Priority 2: "The Protective Shield" (Longitudinal Persistence)
**Title**: *Does long-term ASB colonization protect against symptomatic recurrence?*
*   **Why**: Directly informs antibiotic stewardship. If persistent ASB prevents new virulent strains from entering, treating it is harmful.
*   **Data Support**: `pairwise_metrics.csv` (persistence). `status_map.csv` (outcomes).
*   **Action**: Calculate "Duration of Colonization" for every strain. Correlate duration with frequency of UTI episodes.
*   **Impact**: High. Clinical policy relevance.

### Priority 3: "The Usual Suspects" (Lineage-Specific Virulence)
**Title**: *Virulence potential is lineage-dependent: Defining "High-Risk" clones in the elderly.*
*   **Why**: Not all *E. coli* are equal. If ST131-H30 is 90% UTI and ST73 is 90% ASB, we can risk-stratify patients based on admission swabs.
*   **Data Support**: `mlst_with_meta.csv`, `gwas_multivariable_glmm.csv`.
*   **Action**: Aggregate UTI rate by ST. Run GWAS *within* major STs to remove population structure bias.
*   **Impact**: Moderate/High. Refines the "virulence factor" story.

---

## 5. Suggested Implementation Sketches

### Sketch A: The "Episode Master" & Transition Table (For Priority 1 & 2)
We need a new data structure that moves beyond "pairs" to "timelines".

**New Script**: `15_longitudinal_patterns.R`
1.  **Input**: `status_map.csv`, `pairwise_metrics.csv`.
2.  **Logic**:
    *   Group by `Participant_id`.
    *   Order episodes by time.
    *   Use `pairwise_metrics` to assign a global `Strain_ID` (e.g., "Strain_A", "Strain_B") to every isolate in the timeline. (If Isolate 1 == Isolate 2, they share Strain_ID).
3.  **Output Table**: `participant_timelines.csv`
    *   `Participant`, `Timepoint`, `Date`, `Status`, `Strain_ID`, `Duration_Carriage_Days`.
4.  **Analysis**:
    *   Count "ASB $\to$ UTI" transitions (Same Strain).
    *   Count "ASB $\to$ UTI" replacements (New Strain).

### Sketch B: The "Chameleon" Variant Caller (For Priority 1)
**New Script**: `16_within_host_evolution.R`
1.  **Input**: The list of "Same Strain, Different Status" pairs from Sketch A.
2.  **Logic**:
    *   Loop through these specific pairs.
    *   Run `snippy` or `nucmer` (high sensitivity) on Assembly A vs Assembly B.
    *   Filter for non-synonymous SNPs or intergenic SNPs (promoters).
    *   Check for plasmid gain/loss (using `plasmidfinder` diffs).
3.  **Output**: `evolutionary_events.csv` (Participant, Strain, EventType, GeneAffected).

### Sketch C: Lineage-Adjusted GWAS (For Priority 3)
**Enhance**: `14_genotype_phenotype_model.R`
1.  **Logic**:
    *   Add `ST` (Sequence Type) as a random effect or fixed covariate in the GLMM: `Outcome ~ Gene + (1|Participant) + (1|ST)`.
    *   Alternatively, run the model *separately* for the top 3 STs (e.g., "GWAS in ST131 only").
2.  **Output**: `gwas_adjusted_by_lineage.csv`.

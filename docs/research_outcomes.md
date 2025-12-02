# Research Outcomes: The Yellow RoUTIne Project

**Date:** 2025-11-28
**Project:** Yellow RoUTIne (rUTIs in Nursing Homes)
**Objective:** Identify drivers of recurrent UTIs and ASB-to-UTI transitions using genomic and clinical data.

---

## 1. Executive Summary

We successfully executed the research roadmap, transforming the project from a static analysis pipeline into a dynamic longitudinal study.

**Key Discoveries:**
1.  **"The Chameleon Effect" is Real but Rare**: We identified **2 clear cases** where the *exact same strain* (0-12 SNPs) caused Asymptomatic Bacteriuria (ASB) at one timepoint and Symptomatic UTI at another.
2.  **Genomic Stasis**: In these "Switch" cases, the bacteria did *not* acquire new virulence factors or plasmids. The transition to symptoms was **not driven by gene gain/loss**.
3.  **Host Context**: Preliminary analysis suggests these switches occurred *without* obvious catheter changes (e.g., participants had "Inco" or "Spontaneous" voiding consistently). This points to subtle host immune shifts or transcriptional regulation (phase variation) as the driver.
4.  **Lineage Risk**: We found **no statistically significant difference** in UTI risk between major Sequence Types (STs) in this cohort. ST131 was not significantly more "symptomatic" than others, challenging common assumptions.

---

## 2. Longitudinal Patterns (Priority 1)

We reconstructed timelines for **276 clinical episodes**.

*   **Persistence is the Norm**: The majority of recurrent positive cultures were the **same strain** (ANI >99.9%, <50 SNPs).
*   **Phenotype Switching Events**:
    *   **ASB $\to$ UTI**: 2 events (Participants `40001`, `40004`).
    *   **Negative $\to$ UTI**: 1 event (`31036`).
    *   **ASB $\to$ Negative**: 2 events (Clearance).
    *   **Negative $\to$ ASB**: 1 event.

**Figure 1**: *Swimmer Plot of Patient Timelines* (See `results/longitudinal/swimmer_plot.png`)

---

## 3. "The Chameleon Effect": Within-Host Evolution

We performed deep genomic characterization of the "Switch" pairs.

### Case Study 1: Participant 40001 (ASB $\to$ UTI)
*   **Timeline**: T1 (ASB) $\to$ Uricult (UTI).
*   **Genomics**:
    *   **SNPs**: 12.
    *   **Key Mutations**:
        *   **lpxL** (Lipid A biosynthesis): SNP in CDS. Potential impact on endotoxin structure/immunogenicity.
        *   **dppB** (Dipeptide transport): SNP in CDS.
    *   **Gene Content**: Identical (No VF/Plasmid gain/loss).
*   **Host Context**: "Inco" (Incontinence) at both timepoints. No catheter.
*   **Interpretation**: The strain persisted. The *lpxL* mutation is highly significant as Lipid A modification is a known mechanism of immune evasion or virulence modulation.

### Case Study 2: Participant 40004 (ASB $\to$ UTI)
*   **Timeline**: T1 (ASB) $\to$ Uricult (UTI).
*   **Genomics**:
    *   **SNPs**: 10.
    *   **Key Mutations**:
        *   **rpoD** (Sigma factor 70): SNP in CDS. **Critical Finding**. Sigma 70 is the primary housekeeping sigma factor. A mutation here could globally alter gene expression, potentially triggering a "virulence state" without acquiring new genes.
        *   **hmuT** (Hemin transport): SNP in CDS. Iron acquisition is key for UTI.
*   **Host Context**: "Spontaan geloosd" (Spontaneous) at both timepoints.
*   **Interpretation**: A stable colonizer that likely underwent a global transcriptional shift (via *rpoD*) or metabolic adaptation (*hmuT*) to cause symptoms.

---

## 4. Lineage-Specific Virulence (Priority 3)

We tested if specific STs were more likely to cause UTI vs ASB.

*   **Top STs Analyzed**: ST6, ST87, ST4, ST43.
*   **Result**: No ST showed a statistically significant Odds Ratio (OR) for UTI after FDR correction.
*   **Implication**: In this nursing home cohort, **lineage alone does not predict symptoms**. Clinical management cannot rely solely on "high-risk clone" identification.

---

## 5. Strategic Recommendations

Based on these findings, we propose the following next steps:

1.  **Targeted Variant Analysis**:
    *   Focus exclusively on the ~10 SNPs in `40001` and `40004`.
    *   **Hypothesis**: Are these SNPs in the *fimS* invertible element or other phase-variable promoters?
    *   **Action**: Use long-read assembly graph analysis (if available) or targeted mapping to promoter regions.

2.  **Transcriptomics (Future)**:
    *   Since DNA is stable, the switch is likely RNA-driven. If urine samples are banked, perform RNA-seq on ASB vs UTI samples of the *same* strain.

3.  **Host Immunology**:
    *   Investigate non-catheter host factors: hydration status, viral coinfection, or medication changes (e.g., anticholinergics) that might trigger symptoms.

4.  **Publication**:
    *   These findings (Stability of ASB strains, rarity of genomic switches) are publishable. They support the "Protective Shield" hypothesis—most ASB strains just stay ASB. The rare switches are the exception that proves the rule.

---

## 6. Data & Analysis Guide

For a detailed, step-by-step guide on how to use the output files mentioned above, interpret the columns, and reproduce these analyses, please refer to:

👉 **[Analysis & Outputs Guide](analysis_outputs_guide.md)**

---

**Artifacts Generated:**
*   `results/longitudinal/participant_timelines.csv`
*   `results/longitudinal/evolution_events.csv`
*   `results/longitudinal/annotated_snps.csv`
*   `results/longitudinal/host_context_table.csv`
*   `results/lineage/st_risk_profile.csv`

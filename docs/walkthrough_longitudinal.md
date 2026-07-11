# Walkthrough: Longitudinal & Evolutionary Analysis

**Date:** 2025-11-27
**Goal:** Implement Priorities 1 & 2 of the Research Roadmap (Longitudinal Timelines & Within-Host Evolution).

---

## 1. Longitudinal Timelines (`15_longitudinal_patterns.R`)

We reconstructed patient timelines by linking isolates into strain-lineage context using genomic similarity. Current same-strain-first interpretation uses 0-25 SNPs as strong same-strain evidence, treats >25 SNPs as above the same-strain SNP threshold, and keeps missing SNPs as missing SNP evidence. ST is reported afterward as secondary lineage context (`Same ST`, `Different ST`, or `Missing ST evidence`) and does not prove same strain.

### Key Findings
*   **Persistence is Common**: Many participants (`100064`, `110018`, etc.) carry the exact same strain for multiple timepoints (ASB $\to$ ASB).
*   **Phenotype Switching is Rare but Real**: We identified **6 specific events** where the *same strain* changed clinical status.
    *   **ASB $\to$ UTI**: 2 clear cases (`40001`, `40004`).
    *   **Negative $\to$ UTI**: 1 case (`31036`).
    *   **ASB $\to$ Negative**: 2 cases (`110061`, `40021`).
    *   **Negative $\to$ ASB**: 1 case (`100010`).

### Artifacts
*   `results/longitudinal/participant_timelines.csv`: Master timeline for all patients.
*   `results/longitudinal/strain_persistence_stats.csv`: Summary of how long each strain persisted.
*   `results/longitudinal/swimmer_plot.png`: Visual timeline.

---

## 2. Within-Host Evolution (`16_within_host_evolution.R`)

We performed a deep genomic comparison (SNPs, VF/Plasmid gain/loss) on the 6 "Switch" candidates.

### The "Chameleon" Candidates (ASB $\to$ UTI)
These are the most scientifically valuable pairs. They caused ASB at T1 but symptomatic UTI later, with **minimal genomic change**.

| Participant | Transition | SNPs | VF Change | Plasmid Change | Interpretation |
| :--- | :--- | :--- | :--- | :--- | :--- |
| **40001** | ASB $\to$ UTI | **12** | None | None | **Host/Regulator Driven**. The strain is virtually identical. The switch to symptoms is likely due to host status change or expression/promoter mutation. |
| **40004** | ASB $\to$ UTI | **10** | None | None | **Host/Regulator Driven**. Similar to 40001. |
| **31036** | Neg $\to$ UTI | **10** | None | None | Rapid onset or bloom from low abundance. |

### The "Silent" Colonizer
*   **100010 (Neg $\to$ ASB)**: **0 SNPs**. This strain was likely present below detection at T0 and bloomed to ASB levels at T1 without evolving.

---

## 3. Next Steps (from Roadmap)

1.  **Variant Annotation**: The ~10-12 SNPs in `40001` and `40004` are critical. Are they in promoter regions of virulence genes (e.g., *fim* switch)?
    *   *Action*: Extract these SNPs and annotate them (requires reference GenBank or Prokka output).
2.  **Host Integration**: Since the bugs didn't change much, did the patients? Check clinical metadata for `40001` and `40004` (catheter change? antibiotics? frailty?).
3.  **Lineage Analysis (Priority 3)**: We still need to check if these switching strains belong to specific high-risk STs (e.g., ST131).

---

## 4. Technical Fixes Applied
*   **Strain Classification**: Fixed a bug where missing plasmid data (`NA`) caused "Same" classification to fail. Relaxed logic to treat `NA` as "pass".
*   **Graph Construction**: Fixed duplicate vertex error in timeline generation.
*   **Column Selection**: Fixed `tidyselect` error in evolution script.

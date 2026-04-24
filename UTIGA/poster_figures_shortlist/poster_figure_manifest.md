# UTIGA Poster Figures — Shortlist Manifest

This document outlines the top 5 recommended figures for the UTIGA poster, ready for placement. All figures have been generated specifically using verified, deduplicated, episode-level data.

## Top 5 Poster Figure Recommendations

---

### 1. Cohort Flow Diagram
*   **Filename**: `fig1_cohort_flow.png`
*   **What it shows**: CONSORT-style exclusion flowchart showing participant and episode counts at each analysis step.
*   **Source script**: `scripts/poster_02_flow_diagram.R`
*   **Status**: ✅ **Final-ready**.
*   **Caption Text**: *Participant and episode flow through the YELLOW RoUTIne cohort. Episodes without CFU data or untypeable MLST results ('-') were excluded from the primary longitudinal stability analysis.*

---

### 2. Overall Top ST Distribution
*   **Filename**: `fig2_st_distribution.png`
*   **What it shows**: Horizontal bar chart of the top Sequence Types, restricted to typeable isolates in the ≥2-timepoint cohort. Contains exact percentages.
*   **Source script**: `scripts/poster_03_st_bar.R`
*   **Status**: ✅ **Final-ready**.
*   **Caption Text**: *Distribution of the most prevalent E. coli sequence types among nursing home residents with ≥2 follow-up timepoints. Only typeable isolate episodes are shown. Globally disseminated lineages ST43, ST4, and ST6 dominate the community landscape.*

---

### 3. Within-Host ST Stability (Headline Finding)
*   **Filename**: `fig3_stability.png`
*   **What it shows**: Side-by-side grouped bars comparing pair-level ST stability in the ≥2 timepoint vs ≥3 timepoint cohorts.
*   **Source script**: `scripts/poster_04_stability.R`
*   **Status**: ✅ **Final-ready**. Ensures that "pair-level" consecutive continuity matches abstract logic.
*   **Caption Text**: *Within-host E. coli strain stability over time. Stability was assessed across consecutive sampling pairs for residents within the primary (≥2 TP) and extended follow-up (≥3 TP) sub-cohorts, demonstrating high baseline colonization persistence.*

---

### 4. ASB vs UTI Lineage Comparison
*   **Filename**: `fig4_st_by_status.png`
*   *(Alternate candidate: `fig4_st_distribution_stacked_alt.png`)*
*   **What it shows**: Grouped bar chart comparing the proportional composition of STs in ASB vs UTI episodes.
*   **Source script**: `scripts/poster_06_asb_vs_uti.R`
*   **Status**: ✅ **Final-ready**. (Use the `_alt` version if you prefer stacked bars; use the primary if you prefer grouped side-by-side).
*   **Caption Text**: *Distribution of top E. coli sequence types stratified by clinical episode presentation (ASB vs. UTI). While sample numbers for UTI are limited, lineages like ST10 and ST6 appear proportionately overrepresented in symptomatic episodes.*

---

### 5. Participant Trajectories
*   **Filename**: `fig5_trajectories.png`
*   **What it shows**: Curated swimmer/trajectory plots for individuals with ≥3 timepoints, visually demonstrating both stable carriage and strain-switching events during UTI.
*   **Source script**: `scripts/poster_05_trajectories.R`
*   **Status**: ✅ **Final-ready**. 
*   **Caption Text**: *Longitudinal ST trajectories for selected E. coli carriers (≥3 timepoints). Shape denotes clinical status (▲ = ASB, ● = UTI); color denotes Sequence Type. Illustrates individual-level dynamics ranging from single-strain persistence to lineage turnover.*

---

## Poster Storyline Integration
For the most cohesive flow, arrange the poster Left-to-Right:
**Methods/Cohort**: Figure 1 (Flow Diagram)
**Result Section 1 / Population Level**: Figure 2 (Distribution) & Figure 4 (ASB vs UTI)
**Result Section 2 / Within-Host Level**: Figure 3 (Stability) & Figure 5 (Trajectories)

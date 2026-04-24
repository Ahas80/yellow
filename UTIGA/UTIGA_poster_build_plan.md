# UTIGA Poster Build Plan

## The Core Concept
The single strong scientific story is **high within-host persistence despite population-level diversity**. Your poster must visually scream this. Do not clutter it with tables.

---

## Poster Structure (3-Column Layout)

### 1. Header & Take-Home Box (Top)
*   **Title:** *Persistent Pathogens: High Within-Host Lineage Stability of Escherichia coli in Nursing Home Residents*
*   **Authors & Affiliations:** Standard format under the title.
*   **Take-Home Box (Center-Top, bright background):** 
    *   *Text:* "Within-host E. coli lineages are remarkably stable over 18 months (~80% stability), despite massive strain diversity across the nursing home population."

### 2. Left Column: Background & Methods (25% Width)
*   **Background (3 bullets max):**
    *   *E. coli* causes both ASB and UTI in older adults.
    *   Unclear if lineage turnover inside the host drives the clinical transition to UTI.
*   **Methods:** Brief bullet points detailing the YELLOW RoUTIne cohort, 3-month sampling, Uricults during UTIs, and MLST typing.
*   **Key Numbers Panel (Text/Icons, NO TABLE):**
    *   83 Participants (≥2 timepoints)
    *   181 Isolates Sequenced
    *   38 Distinct Sequence Types (STs) identified
*   **Figure:** `fig1_cohort_flow.png`
    *   *Size:* Small (bottom of column).
    *   *Caption:* "**Figure 1: Cohort Inclusion Flow.** Longitudinal sampling strategy capturing routine 3-monthly visits and interim suspected UTI episodes."

### 3. Center Column: The Main Result (50% Width)
This is the heart of the poster. It proves the 80% stability claim.

*   **Top Center Figure:** `fig5_trajectories.png`
    *   *Size:* Huge Centerpiece (occupies top 2/3 of the column).
    *   *Header:* "Longitudinal Dynamics Show Persistent Lineages Through UTIs"
    *   *Caption:* "**Figure 2: Illustrative Longitudinal Lineage Trajectories.** A curated subset of patient timelines demonstrating both persistent colonization (same colour across timepoints) and minority lineage switching events. Across the entire ≥3-timepoint cohort, 12 residents experienced a UTI; sequence types (STs) remained strictly consistent through the UTI in 66% (n=8) of cases, while only 4 residents displayed an ST switch."
*   **Bottom Center Figure:** `fig3_stability.png`
    *   *Size:* Medium (bottom 1/3 of the column).
    *   *Header:* "80.2% Overall Isolate Pair Stability"
    *   *Caption:* "**Figure 3: Pair-Level Lineage Stability.** Among 96 consecutive isolate pairs from the ≥2-timepoint sub-cohort, STs remained identical in 77 cases (80.2%), demonstrating tight within-host persistence. Stability increased slightly to 81.8% in residents sampled at ≥3 timepoints."

### 4. Right Column: Population Context & Conclusion (25% Width)
This section proves the diversity claim.

*   **Top Right Figure:** `fig4_st_by_status.png`
    *   *Size:* Medium.
    *   *Header:* "Lineage Diversity is High, Led by Globally Disseminated STs"
    *   *Caption:* "**Figure 4: Sequence Type Distribution by Episode Status.** Across the population, 38 distinct STs were identified. ST43 (11.6%), ST4 (7.7%), and ST6 (6.6%) dominated overall. In suspected UTI episodes specifically, ST6 and ST10 were the most common isolates recovered."
*   **Scientific Caveat (Small text box):** 17.1% of isolates (n=31) were untypeable due to sequencing ambiguity and excluded from lineage turnover analysis.
*   **Conclusion (Bottom Right):**
    *   *E. coli* populations in nursing homes are diverse but dominated by a few lineages.
    *   Colonization is persistent: individuals retain strains for many months.
    *   Most suspected UTIs are caused by the patient’s "resident" strain, not a novel acquired lineage.

---

## Exactly What I Should Place on the Poster

### 1. The Final Figures to Use
*   **`UTIGA/poster_figures/fig5_trajectories.png`** (Centerpiece)
*   **`UTIGA/poster_figures/fig3_stability.png`** (Center bottom, supporting metric)
*   **`UTIGA/poster_figures/fig4_st_by_status.png`** (Right column, showing diversity)
*   **`UTIGA/poster_figures/fig1_cohort_flow.png`** (Left column, small methodology)
*   *(Drop `fig0_all_swimmer.png` and `fig2_st_distribution.png` entirely).*

### 2. The Final Tables to Use or Not Use
*   **DO NOT USE ANY TABLES.** 
*   Tables stall visual momentum. Use the "Key Numbers Panel" text format (outlined in Left Column) instead of a grid.

### 3. The Exact Poster Section Where Each Item Goes
*   **Top/Header:** Take-Home Message Box.
*   **Left Column:** Background, Methods, Key Numbers Box, `fig1_cohort_flow.png`.
*   **Center Column:** `fig5_trajectories.png` (Top), `fig3_stability.png` (Bottom).
*   **Right Column:** `fig4_st_by_status.png` (Top), Caveat Box, Conclusions Box (Bottom).

### 4. Any Scripts I Need to Run Before Building
*   None. Your existing scripts (`poster_04_stability.R`, `poster_05_trajectories.R`, `poster_06_asb_vs_uti.R`) have already generated the required poster figures in the `poster_figures` folder. The counts have been extracted and verified against the abstract. You are ready to drag and drop these images directly into PowerPoint or Canva.

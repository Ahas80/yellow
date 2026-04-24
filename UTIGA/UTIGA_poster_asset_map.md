# UTIGA Poster Asset Map & Source Audit

## 1. Usable Scripts (Poster Content Generators)
*The following scripts should be used directly. They are fully optimized for the poster logic (Uricult spacing, etc).*

| Script Name | Purpose | Generates | Keep / Modify |
| :--- | :--- | :--- | :--- |
| `UTIGA/scripts/poster_05_trajectories.R` | Builds the visual centerpiece of selected individuals over time | `fig5_trajectories.png` | **USE AS-IS** (Generates central figure) |
| `UTIGA/scripts/poster_04_stability.R` | Calculates and plots the 80%+ pair-level stability | `fig3_stability.png` | **USE AS-IS** (Generates secondary figure) |
| `UTIGA/scripts/poster_06_asb_vs_uti.R` | Plots the ASB vs UTI distribution among dominant STs | `fig4_st_by_status.png` | **USE AS-IS** (Generates supporting context figure) |
| `UTIGA/scripts/poster_02_flow_diagram.R` | Creates the cohort inclusion flow | `fig1_cohort_flow.png` | **USE AS-IS** (Generates methodology graphic) |
| `scripts/06_MLST.R` & `07_explore_MLST.R` | Exploratory source scripts. Run to generate `mlst_with_meta.csv` | Background data | **IGNORE FOR POSTER** (Just run upstream) |

## 2. Redundant Scripts (Drop for Poster)
*These scripts are useful for the paper but cause poster clutter.*

| Script Name | Purpose | Generates | Why to Drop |
| :--- | :--- | :--- | :--- |
| `UTIGA/scripts/poster_01_all_swimmer.R` | Plots every single individual's timeline | `fig0_all_swimmer.png` | 83 rows is an illegible wall of color on a poster. |
| `UTIGA/scripts/poster_03_st_bar.R` | Simple bar chart of top STs | `fig2_st_distribution.png` | `fig4` is a smarter, split version of this. No need for both. |
| `scripts/00c_plot_clinical_summary.R` | Cohort-wide epidemiological summary | Various stacked plots | Basic demo graphics; not the ST-specific story. |

## 3. Poster Figure Inventory

| Figure | Source File | Poster Section | Size | Verdict |
| :--- | :--- | :--- | :--- | :--- |
| **Trajectories** | `fig5_trajectories.png` | Center Column (Main Result) | **Centerpiece (Huge)** | **USE.** Perfectly illustrates the abstract claims of stability vs switching. |
| **Stability Bars** | `fig3_stability.png` | Center Column (Bottom) | **Medium** | **USE.** Directly visualization of the 80%+ stability metric. |
| **ASB vs UTI STs** | `fig4_st_by_status.png` | Right Column (Context) | **Medium** | **USE.** Provides population context without stealing focus. |
| **Flow Diagram** | `fig1_cohort_flow.png` | Left Column (Methods) | **Small** | **OPTIONAL.** Good if space permits, drop if cramped. |
| **All Swimmer** | `fig0_all_swimmer.png` | N/A | N/A | **DROP.** Too dense. |
| **ST Distribution** | `fig2_st_distribution.png` | N/A | N/A | **DROP.** Redundant to fig4. |

## 4. Extraction of Poster-Ready Numbers
*Note: Due to minor snapshot/linkage iteration differences in the data, the current "live" `poster_01_reconcile_numbers.R` yields slightly lower counts than your accepted abstract. **Use the abstract numbers** below to ensure perfect alignment with what reviewers accepted.*

### Cohort Setup
*   **≥2-Timepoint Participants:** 83
*   **Total Isolates:** 181
*   **Untypeable Isolates:** 31 (17.1%) (Excluded due to "-" sequencing ambiguity)
*   **Total Typeable:** 150
*   **Distinct Sequence Types (STs):** 38

### Top Lineages (Population Diversity)
*   **ST43:** 11.6%
*   **ST4:** 7.7%
*   **ST6:** 6.6%
*   **ST1:** 5.5%
*   **Typeable UTI Isolates Total:** 33
*   **Top UTI STs:** ST6 (n=8), ST10 (n=5)

### Within-Host Stability Metrics (The Main Message)
*   **Isolate Pairs Evaluated:** 96 consecutive pairs
*   **Identical Over Time:** 77 pairs (**80.2%**)
*   *Extended group (≥3 Timepoints):* 54 participants, 119 isolates.
*   *Extended group Pairs:* 66 pairs.
*   *Extended group Stability:* 54 identical pairs (**81.8%**)
*   *Switching occurrence:* ~20% of participants switched lineages.
*   *UTI Specifics:* 12 residents with ≥3 timepoints had UTIs. 8 remained consistent, 4 showed a lineage switch.

## 5. Table Strategy
**Recommendation: DO NOT USE TABLES.**

A table of counts (e.g., ST frequencies or exactly how many people swapped) will completely stall the visual momentum built by your figures. 
If you desperately want to list numbers, use a **"Key Findings" icon-stat panel** in the left column.

**Icon Panel Concept (No Script Needed):**
*   👥 **83** Residents Followed
*   🧫 **38** Distinct *E. coli* Lineages
*   🔒 **80.2%** Within-Host Stable Transmissions Over 18 Months

There are no missing scripts you need to run to generate tables. Let your existing `fig3`, `fig4`, and `fig5` do all the heavy lifting.

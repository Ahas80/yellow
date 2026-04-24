# UTIGA Conference Poster Strategy & Blueprint

## A. Executive Recommendation
**The core scientific story for UTIGA:** *Escherichia coli* colonization in nursing home residents is highly diverse across the population but displays massive stability (~80%) within individual hosts, even across 18 months, regardless of whether they develop a UTI or remain asymptomatically colonized (ASB). 

**The Visual Strategy:** Avoid cluttering the poster with generic demographic charts. Make the story about **longitudinal persistence**. Your trajectory plot (`fig5`) should be the centerpiece, supported by the pair-level stability bar chart (`fig3`). The message must visually scream: "People keep the same *E. coli* strain for a very long time, even when they get sick."

## B. Best Poster Title Options
1. **Persistent Pathogens:** High Within-Host Lineage Stability of *Escherichia coli* in Nursing Home Residents During ASB and UTI
2. **Longitudinal Dynamics of *E. coli* Sequence Types in Older Adults:** High Within-Host Stability Despite Population-Level Diversity
3. **Within-Host Lineage Stability of *E. coli* Over 18 Months:** A Genomic Analysis of the YELLOW RoUTIne Cohort
4. *E. coli* Sequence Type Dynamics in the YELLOW RoUTIne Cohort: Dominance and Persistence of Globally Disseminated Lineages
5. **From ASB to UTI:** ~80% Within-Host Strain Stability Over 18 Months in a Dutch Nursing Home Cohort

*(Recommendation: Go with Title 1 or 2. They emphasize the most interesting numerical finding right in the headline).*

## C. Recommended Poster Layout
**Format:** 3-Column Billboard Style (Wide orientation, e.g., 48" x 36")

*   **Top Banner:** Title, Authors, Affiliations.
*   **Column 1 (Left - 25% width): Context & Setup**
    *   *Background:* 3 bullet points max.
    *   *Methods & Cohort Flow:* Brief text + `fig1_cohort_flow.png` (small) to establish the robust methodology.
*   **Column 2 (Center - 50% width): The Main Event (Results)**
    *   *Key Finding Box:* Massive font explicitly stating the 80% stability.
    *   *Visual Centerpiece:* `fig5_trajectories.png` showing the 10-15 selected illustrative cases.
    *   *Supporting Finding:* `fig3_stability.png` (the pair-level stability) positioned below or beside it.
*   **Column 3 (Right - 25% width): Population Context & Walkaway**
    *   *Population Diversity:* `fig4_st_by_status.png` downsized to show that while individuals are stable, the cohort is diverse.
    *   *Conclusions:* 3 clear take-home bullet points.
    *   *QR Code:* Link to the pre-print or GitHub repo if available.

## D. Recommended Figures
| Rank | File Path | Keep / Drop | Why | Required Edits | Suggested Poster Size |
| :--- | :--- | :--- | :--- | :--- | :--- |
| **1** | `UTIGA/poster_figures/fig5_trajectories.png` | **USE** | It visually proves the "within-host stability vs switching" claim perfectly by showing real patient timelines. | Ensure y-axis patient IDs aren't too small. Reduce the legend box size. | **Large (Centerpiece)** |
| **2** | `UTIGA/poster_figures/fig3_stability.png` | **USE** | Directly visualizes the 80.2% and 81.8% abstract numbers. | Increase the font size of the percentages. | **Medium (Below Centerpiece)** |
| **3** | `UTIGA/poster_figures/fig4_st_by_status.png` | **USE** | Better than fig2 because it groups the dominant lineages (ST6, ST10) by ASB vs UTI, hitting your abstract note. | Simplify labels. Remove the x-axis "%" grid lines to increase whitespace. | **Small/Medium (Right Col.)** |
| **4** | `UTIGA/poster_figures/fig1_cohort_flow.png` | **USE (Optional)** | Validates the YELLOW pipeline visually. | Shrink significantly. | **Small (Left Col.)** |
| **5** | `UTIGA/poster_figures/fig0_all_swimmer.png` | **DROP** | 83 individuals is far too dense for a conference poster. | None. Do not use. | N/A |
| **6** | `UTIGA/poster_figures/fig2_st_distribution.png` | **DROP** | Redundant if you use `fig4_st_by_status.png`. | None. Do not use. | N/A |

## E. Poster Text Draft

### Take-Home Message Box (Top Center)
> **Within-host *E. coli* lineages are remarkably stable over 18 months (~80% stability), despite massive strain diversity across the nursing home population.**

### Background
*   *Escherichia coli* is the dominant uropathogen in older adults, causing both asymptomatic bacteriuria (ASB) and symptomatic urinary tract infections (UTI).
*   It is currently unclear if lineage turnover inside the host drives the transition from ASB to UTI.
*   **Objective:** Characterize longitudinal *E. coli* Sequence Type (ST) dynamics to assess lineage stability during routine colonization and suspected UTI episodes.

### Methods
*   **Cohort:** YELLOW roUTIne nursing home cohort (Older adults, Netherlands).
*   **Sampling:** Routine urine cultures every 3 months for up to 18 months + interim Uricults during suspected UTIs. 
*   **Sequencing & Typing:** Whole Genome Sequencing (WGS) + Multilocus Sequence Typing (MLST) on isolates from episodes clinically classified as ASB or UTI.

### Results
**1. Exceptional Within-Host Stability** (Pair with Fig 3 & 5)
*   Among 96 consecutive isolate pairs, **STs remained identical in 80.2%** (77/96), demonstrating tight within-host persistence.
*   In the extended tracking group (≥3 timepoints; 54 participants), stability increased slightly to 81.8%.
*   **UTI Dynamics:** 12 residents with ≥3 timepoints experienced a UTI. In 8 cases (66%), the ST remained strictly consistent throughout follow-up, while 4 residents showed an ST switch.

**2. Population-Level Diversity** (Pair with Fig 4)
*   **High Diversity:** In the ≥2-timepoint sub-cohort (83 participants, 181 isolates), 38 distinct sequence types (STs) were identified.
*   **Dominant Lineages:** Most prevalent STs were globally disseminated lineages: ST43 (11.6%), ST4 (7.7%), and ST6 (6.6%).
*   **UTI-Specific:** Among 33 typeable UTI isolates, ST6 (n=8) and ST10 (n=5) were the most frequently observed. *(17.1% of isolates were untypeable due to sequencing ambiguity and excluded).*

### Conclusions
*   *E. coli* populations in nursing homes are highly diverse but dominated by a few globally successful lineages.
*   **Colonization is persistent:** Individuals tend to retain the same *E. coli* lineage over many months (~80% within-host stability).
*   Most suspected UTIs are caused by the patient’s persistently colonizing "resident" strain rather than a newly acquired lineage.

## F. Figure-by-Figure Edit Instructions
1.  **`fig3_stability.png` (Stability Bar Plot)**
    *   **Edit:** The subtitle "Consecutive isolate pairs from ASB/UTI episodes" should be larger. 
    *   **Edit:** Make the text labels inside the bars (e.g., "77 (80.2%)") white, bold, and massive (like +4pt size up). 
    *   **Edit:** Remove the dark vertical axis line entirely.
2.  **`fig5_trajectories.png` (Trajectories Plot)**
    *   **Edit:** The legend for shapes (triangle = ASB, circle = UTI) needs to be positioned closer to the top left of the actual plotting area to avoid blank space on the right. 
    *   **Edit:** Increase the thickness of the horizontal grey connecting lines slightly to visually group the patient's timeline better from a distance.
3.  **`fig4_st_by_status.png` (ASB vs UTI Bar Plot)**
    *   **Edit:** The x-axis is cluttered. Remove the minor grid lines.
    *   **Edit:** Use the exact same green/red or blue/red hex codes used in Fig 3 and 5 to ensure visual consistency across the poster (e.g., ASB = Blue/Yellow, UTI = Red). Right now, the scripts might be using different palettes. Red `#EF5350` for UTI is great, make sure `fig5` uses the exact same red for UTI circles.

## G. Missing Figure(s) to Create
No *new* complex figures are needed. Your scripts (`poster_04_stability.R` and `poster_05_trajectories.R`) perfectly map to the abstract claims. 

*Optional Addition:* A **Graphical Abstract Icon** showing a generic "Patient -> Time 1 (Blue Bug) -> Time 2 (Blue Bug)" vs "Patient -> Time 1 (Blue Bug) -> Time 2 (Red Bug)" would make the 80% stability metric instantly understandable to someone walking by in 3 seconds.

## H. Risks / Scientific Cautions
1.  **Untypeable Isolates (17.1%):** Be prepared for questions about the 31 untypeable isolates. Someone *will* ask if the switching rate is artificially deflated because the messy/untypeable ones were excluded. Have a 1-sentence defense ready (e.g. "Exclusion was strictly due to WGS short/long read ambiguity, not biological absence, and was distributed randomly across patients").
2.  **"Causation" vs "Association":** You wrote "4 residents showed a switch in STs" during a UTI. Be careful not to claim the switch *caused* the UTI on the poster, only that a UTI *coincided* with an ST turnover event in a minority of cases.
3.  **ASB vs UTI classification:** Because the poster doesn't rigorously define clinical UTI vs ASB, someone will poke at how you distinguished them in elderly patients (where symptoms are notoriously vague). Mention "suspected UTI based on nursing home physician assessment / Uricult request" if asked.

## I. Final "Ready-to-build" Assembly Plan

### Canva/PowerPoint Assembly Order Step-by-Step
1.  **Set up the canvas:** Set slide size to your conference specs (e.g., 48" x 36"). Add a clean, light background (pure white or very light grey/blue).
2.  **Header:** Place the Title, Authors, and Logos at the top. Use dark blue or black for the title, sans-serif font (Arial, Helvetica, or Roboto).
3.  **Central Highlight Box:** Immediately below the header in the center, create a colored box (e.g., soft blue or yellow) and paste the "Take-Home Message." Making this the visual anchor ensures everyone reads your main point.
4.  **Columns:** Draw 3 invisible or faintly bordered columns. 
5.  **Left Column:** Paste Background text, Methods text, and drop in `fig1_cohort_flow.png`. Keep text bulleted.
6.  **Middle Column:** Drop in `fig5_trajectories.png` (make it huge) and `fig3_stability.png` right below it. This is the heart of the poster. 
7.  **Right Column:** Paste the Results text (Population Diversity section), drop in `fig4_st_by_status.png`, and add the Conclusions text box. 

### Final Submission Checklist
*   [ ] **Color consistency:** Are the colors for "UTI" and "ASB" identical across Fig 3, Fig 4, and Fig 5?
*   [ ] **Font check:** Print a test page on standard A4 paper at 100% scale. Place it on the floor. If you can't read the axis labels while standing over it, the font is too small.
*   [ ] **Numbers match:** Verify the abstract numbers (83 participants, 80.2% stability, 38 STs) exactly match what is printed on the poster panels.
*   [ ] **Vector export:** Export final as a flattened high-res PDF to preserve vector text crispness at printing scale. Do not export as PNG/JPEG.

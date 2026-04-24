# UTIGA Conference Poster — Complete Analysis & Figure Plan

## Overview

This document audits the YELLOW RoUTIne / rUTI project files to build a conference poster based on the submitted MLST abstract. It maps every abstract claim to its source file, recommends which figures to include, flags data inconsistencies, and provides an execution plan.

---

## STEP 1 — Project Audit Summary

### Core MLST Scripts

| Script | Purpose | Key Outputs |
|--------|---------|-------------|
| [06_MLST.R](file:///Users/Aamir/Desktop/rUTIs/06_MLST.R) | Run MLST typing on all assemblies | `results/mlst/mlst_all.tsv`, `top_STs.csv`, `ST_consecutive_concordance.csv`, `ST_persistence_by_participant.csv` |
| [07_explore_MLST.R](file:///Users/Aamir/Desktop/rUTIs/07_explore_MLST.R) | Plot top 20 STs, join metadata | `results/mlst/ST_frequencies.csv`, `mlst_with_meta.csv`, `plots/mlst/top20_STs.png` |
| [15_longitudinal_patterns.R](file:///Users/Aamir/Desktop/rUTIs/15_longitudinal_patterns.R) | Swimmer plot, strain assignments, transitions | `results/longitudinal/swimmer_plot.png`, `participant_timelines.csv`, `transitions.csv` |
| [get_stratified_stats.R](file:///Users/Aamir/Desktop/rUTIs/get_stratified_stats.R) | Stratified stats for ≥2/≥3/≥4 TP cohorts | Console output (not auto-saved) |
| [generate_abstract_stats.R](file:///Users/Aamir/Desktop/rUTIs/generate_abstract_stats.R) | Full abstract number generation | [abstract_stats_output.md](file:///Users/Aamir/Desktop/rUTIs/abstract_stats_output.md) |

### Key Data Files

| File | Role | Rows |
|------|------|------|
| [status_map.csv](file:///Users/Aamir/Desktop/rUTIs/status_map.csv) | Clinical master: participant × timepoint × infection status | 277 |
| [mlst_with_meta.csv](file:///Users/Aamir/Desktop/rUTIs/results/mlst/mlst_with_meta.csv) | MLST calls joined with assembly metadata | 382 |
| [ST_frequencies.csv](file:///Users/Aamir/Desktop/rUTIs/results/mlst/ST_frequencies.csv) | ST distribution (deduplicated, all isolates) | 39 STs |
| [ST_persistence_by_participant.csv](file:///Users/Aamir/Desktop/rUTIs/results/mlst/ST_persistence_by_participant.csv) | Per-participant dominant ST and stability | 88 rows |
| [st_risk_profile.csv](file:///Users/Aamir/Desktop/rUTIs/results/lineage/st_risk_profile.csv) | UTI proportion by ST with Fisher's exact test | 8 STs |

### Existing Figures (Reviewed)

| File | Usable? | Notes |
|------|---------|-------|
| [top20_STs.png](file:///Users/Aamir/Desktop/rUTIs/plots/mlst/top20_STs.png) | ⚠️ Needs adaptation | Uses full cohort (all 382 isolates); "-" (untypeable) is the largest bar. Needs filtering to ≥2-TP subset and excluding untypeable |
| [swimmer_plot.png](file:///Users/Aamir/Desktop/rUTIs/results/longitudinal/swimmer_plot.png) | ⚠️ Needs adaptation | Shows infection status over time but NOT ST identity; good starting point but needs ST overlay |
| [st_distribution_stacked.png](file:///Users/Aamir/Desktop/rUTIs/plots/epidemiology/st_distribution_stacked.png) | ✅ Reusable as-is | ASB vs UTI ST composition — excellent poster candidate |
| [st_risk_plot.png](file:///Users/Aamir/Desktop/rUTIs/results/lineage/st_risk_plot.png) | ✅ Reusable | UTI risk by ST with confidence intervals |
| [swimmer_plot_top20.png](file:///Users/Aamir/Desktop/rUTIs/plots/timelines/swimmer_plot_top20.png) | ❌ Poor quality | Greyscale, missing participants, hard to read — needs full redo |
| [st_persistence_by_participant.png](file:///Users/Aamir/Desktop/rUTIs/results/plots/st_persistence_by_participant.png) | ❌ Broken | Rendering failure (solid grey block), axis labels show "NA" |
| [waterfall_counts.png](file:///Users/Aamir/Desktop/rUTIs/results/plots/waterfall_counts.png) | ✅ Reusable with edits | Episode classification funnel — useful for flow diagram |

---

## STEP 2 — Source-of-Truth Mapping

> [!WARNING]
> **Critical discrepancies exist between the submitted abstract and the most recent script outputs.** The abstract uses numbers from an earlier `get_stratified_stats.R` run. The latest `abstract_stats_output.md` (from `generate_abstract_stats.R`) shows different denominators due to changed linkage/filtering logic.

| Poster Element | Abstract Claim | Best Source File | Supporting Script | Why Best Source | Caveats |
|---|---|---|---|---|---|
| **≥2-TP cohort size** | 83 participants, 181 isolates | `abstract_stats_output.md` says **92 participants, 150 linked** | `generate_abstract_stats.R` | Most recent, most rigorous linkage logic | **MISMATCH**: 83 vs 92 participants; 181 vs 150 isolates. Likely different filtering (83 may exclude Neg-only participants, 181 may count all ASB+UTI episodes vs only linked) |
| **Number of distinct STs** | 38 STs | [ST_frequencies.csv](file:///Users/Aamir/Desktop/rUTIs/results/mlst/ST_frequencies.csv) shows **39 STs** (incl. "-") = **38 typeable** | `07_explore_MLST.R` | Direct file count | 38 matches if "-" excluded; but this is full cohort, not ≥2-TP subset. `abstract_stats_output.md` says 34 STs for ≥2-TP |
| **Untypeable isolates** | 31 (17.1%) | [ST_frequencies.csv](file:///Users/Aamir/Desktop/rUTIs/results/mlst/ST_frequencies.csv): 34 "-" of 191 = 17.8% | `07_explore_MLST.R` | Full cohort deduplicated count | **MISMATCH**: 31 vs 34, 17.1% vs 17.8%. The 31 may come from the ≥2-TP subset |
| **UTI isolates** | 33 typeable UTI episodes | Not directly in any CSV | `get_stratified_stats.R` console output | Needs re-running | Abstract says "thirty-three typeable isolates from UTI episodes" — needs verification |
| **Top STs overall** | ST43 (11.6%), ST4 (7.7%), ST6 (6.6%), ST1 (5.5%) | [abstract_stats_output.md](file:///Users/Aamir/Desktop/rUTIs/abstract_stats_output.md): ST43=13.3%, ST4=8.0%, ST6=7.3%, ST1=4.0% | `generate_abstract_stats.R` | ≥2 TP subset, typeable only | **MISMATCH** vs abstract percentages. These shift depending on whether "-" is excluded from denominator |
| **Top UTI STs** | ST6 (n=8), ST10 (n=5) | [st_risk_profile.csv](file:///Users/Aamir/Desktop/rUTIs/results/lineage/st_risk_profile.csv): ST6 has 3 UTI, not 8 | `17_lineage_analysis.R` | **SIGNIFICANT MISMATCH**: ST6 n=3 UTI, not 8. Abstract number may be outdated or from different subset |
| **Stability ≥2 TP** | 80.2% (77/96 pairs) | [abstract_stats_output.md](file:///Users/Aamir/Desktop/rUTIs/abstract_stats_output.md): **86.2% (56/65 stable participants)** | `generate_abstract_stats.R` | Different metric: abstract counts *pairs*, output counts *participants*. `get_stratified_stats.R` may be closer | Scripts define stability differently (participant-level vs pair-level) |
| **Stability ≥3 TP** | 81.8% (54/66 pairs) | `abstract_stats_output.md`: **86.7% (39/45 stable participants)** | `generate_abstract_stats.R` | Same metric mismatch as above | The 80.2% / 81.8% pair-level figures likely came from `get_stratified_stats.R` |
| **Switching ~20%** | ~20% of participants | `abstract_stats_output.md`: **13.8%** (9/65 for ≥2 TP) | `generate_abstract_stats.R` | **MISMATCH**: 13.8% vs ~20%. Depends on definition and whether "-" STs count as switches |
| **≥3-TP UTI residents** | 12 residents with UTI, 8 stable, 4 switched | No pre-computed file | Needs new analysis | — | Not directly available in any existing output |

> [!IMPORTANT]
> **The abstract numbers and the latest script outputs use different definitions.** The most likely explanation:
> - The abstract's "181 isolates" counts all ASB/UTI **episodes** (not just those with linked STs)
> - The abstract's "80.2% (77/96)" counts **consecutive isolate pairs**, while `generate_abstract_stats.R` counts **participants**
> - The abstract's "~20%" switching may include "-" → typed-ST changes as switches
>
> **Recommendation**: Re-run `get_stratified_stats.R` (which uses pair-level counting) and reconcile. The poster should use whichever definition you choose, applied consistently.

---

## STEP 3 — Figure Selection

### Figure 1: Cohort / Analysis Flow Diagram
- **Include**: ✅ **Must-have** (Priority 1)
- **Why useful**: Anchors the poster, shows sample attrition, ensures transparency about untypeable isolates
- **Scientific question**: How were participants and isolates selected for analysis?
- **Type**: Descriptive (CONSORT-style flow)
- **Poster-worthy**: Yes — every clinical poster needs this
- **Input files**: [status_map.csv](file:///Users/Aamir/Desktop/rUTIs/status_map.csv), [abstract_stats_output.md](file:///Users/Aamir/Desktop/rUTIs/abstract_stats_output.md)
- **Existing**: [waterfall_counts.png](file:///Users/Aamir/Desktop/rUTIs/results/plots/waterfall_counts.png) is a starting point but needs a full CONSORT-style redesign
- **Status**: 🔨 **Needs new figure** — create as a clean flowchart
- **Proposed output**: `poster_figures/fig1_cohort_flow.png`
- **Plot type**: Box-and-arrow CONSORT diagram
- **Message**: From the YELLOW RoUTIne cohort, 92 participants with ≥2 timepoints contributed 150 ST-linked E. coli episodes for longitudinal analysis
- **Caption**: *Participant and isolate flow through the YELLOW RoUTIne cohort. Episodes without CFU data or untypeable MLST results were excluded from ST stability analysis.*

---

### Figure 2: Overall ST Distribution (Top 10–15 STs)
- **Include**: ✅ **Must-have** (Priority 1)
- **Why useful**: Shows the diversity and dominance of key lineages
- **Scientific question**: Which E. coli lineages dominate in NH residents?
- **Type**: Descriptive
- **Poster-worthy**: Yes — core result
- **Input files**: [mlst_with_meta.csv](file:///Users/Aamir/Desktop/rUTIs/results/mlst/mlst_with_meta.csv), [status_map.csv](file:///Users/Aamir/Desktop/rUTIs/status_map.csv)
- **Existing**: [top20_STs.png](file:///Users/Aamir/Desktop/rUTIs/plots/mlst/top20_STs.png) — includes "-" and uses full cohort. **Not directly usable.**
- **Status**: 🔨 **Needs regeneration** — filter to ≥2-TP cohort, exclude untypeable, show top 10–15 STs with percentages
- **Proposed output**: `poster_figures/fig2_st_distribution.png`
- **Plot type**: Horizontal bar chart with percentage labels, colored by known global lineage status
- **Message**: E. coli populations were diverse (34 STs) but dominated by globally disseminated lineages ST43 (13.3%), ST4 (8.0%), and ST6 (7.3%)
- **Caption**: *Distribution of E. coli sequence types among typeable isolates in nursing home residents with ≥2 follow-up timepoints (n=150 episodes, 34 distinct STs). Top 10 STs shown.*

---

### Figure 3: Within-Host ST Stability (Stacked Bar or Donut)
- **Include**: ✅ **Must-have** (Priority 1)
- **Why useful**: The central finding — ~80% stability is the take-home message
- **Scientific question**: How stable are E. coli lineages within individual NH residents over time?
- **Type**: Descriptive
- **Poster-worthy**: Yes — the headline result
- **Input files**: [mlst_with_meta.csv](file:///Users/Aamir/Desktop/rUTIs/results/mlst/mlst_with_meta.csv), [status_map.csv](file:///Users/Aamir/Desktop/rUTIs/status_map.csv)
- **Existing**: `st_persistence_by_participant.png` is broken. No clean version exists.
- **Status**: 🔨 **Needs new figure**
- **Proposed output**: `poster_figures/fig3_stability.png`
- **Plot type**: Side-by-side stacked bar or paired donut chart comparing ≥2-TP and ≥3-TP cohorts (stable vs switched)
- **Message**: ~80% of consecutive E. coli isolate pairs maintained the same ST, regardless of follow-up duration
- **Caption**: *Within-host E. coli ST stability across consecutive sampling timepoints. Stability was consistent between primary (≥2 TP) and extended follow-up (≥3 TP) sub-cohorts.*

---

### Figure 4: ST Distribution by Infection Status (ASB vs UTI)
- **Include**: ✅ **Good-to-have** (Priority 2)
- **Why useful**: Directly addresses whether certain STs are associated with UTI
- **Scientific question**: Do the dominant STs differ between ASB and UTI episodes?
- **Type**: Descriptive (inferential interpretation discussed cautiously)
- **Poster-worthy**: Yes — addresses the core research question
- **Input files**: [mlst_with_meta.csv](file:///Users/Aamir/Desktop/rUTIs/results/mlst/mlst_with_meta.csv), [status_map.csv](file:///Users/Aamir/Desktop/rUTIs/status_map.csv)
- **Existing**: [st_distribution_stacked.png](file:///Users/Aamir/Desktop/rUTIs/plots/epidemiology/st_distribution_stacked.png) — **reusable** (shows ASB vs UTI side-by-side)
- **Status**: ✅ **Can reuse** — but may benefit from poster-quality styling
- **Proposed output**: `poster_figures/fig4_st_by_status.png`
- **Plot type**: Stacked 100% bar chart (ASB vs UTI) or grouped bar chart
- **Message**: ST composition differed qualitatively between ASB and UTI, with ST6 and ST10 over-represented in UTI episodes, but sample sizes preclude formal comparison
- **Caption**: *Proportional distribution of E. coli STs in ASB vs UTI episodes. Due to the low number of UTI episodes (n≈7–33 depending on subset), differences should be interpreted cautiously.*

---

### Figure 5: Participant-Level ST Trajectories (Selected Illustrative Cases)
- **Include**: ✅ **Good-to-have** (Priority 2)
- **Why useful**: Makes longitudinal dynamics tangible — shows real patient journeys
- **Scientific question**: What do within-host ST trajectories look like in individual residents?
- **Type**: Descriptive
- **Poster-worthy**: Yes — provides the "human story" for a poster
- **Input files**: [mlst_with_meta.csv](file:///Users/Aamir/Desktop/rUTIs/results/mlst/mlst_with_meta.csv), [status_map.csv](file:///Users/Aamir/Desktop/rUTIs/status_map.csv), [participant_timelines.csv](file:///Users/Aamir/Desktop/rUTIs/results/longitudinal/participant_timelines.csv)
- **Existing**: [swimmer_plot.png](file:///Users/Aamir/Desktop/rUTIs/results/longitudinal/swimmer_plot.png) — shows infection status but NOT STs. [swimmer_plot_top20.png](file:///Users/Aamir/Desktop/rUTIs/plots/timelines/swimmer_plot_top20.png) — shows STs but poor quality
- **Status**: 🔨 **Needs new figure** — select 8–12 illustrative participants (mix of stable, switchers, UTI)
- **Proposed output**: `poster_figures/fig5_trajectories.png`
- **Plot type**: Swimmer / timeline plot with ST identity as color and infection status as shape, for selected participants with ≥3 TPs
- **Message**: Most residents maintained the same E. coli lineage across all timepoints, regardless of clinical status transitions
- **Caption**: *Selected participant ST trajectories over longitudinal follow-up. Shape denotes infection status (▲ ASB, ● UTI); color denotes ST. Most participants show stable colonization with a single ST.*

---

### Figure 6: Summary Table (Optional)
- **Include**: ⬜ **Optional** (Priority 3)
- **Why useful**: Allows precise numbers to be read directly on the poster
- **Scientific question**: Summary of key cohort and analysis metrics
- **Type**: Descriptive
- **Poster-worthy**: Only if space permits — can be placed in Methods or as a small panel
- **Input files**: [abstract_stats_output.md](file:///Users/Aamir/Desktop/rUTIs/abstract_stats_output.md)
- **Status**: 🔨 **Needs creation** — formatted as a clean 2-column table
- **Proposed output**: `poster_figures/table1_summary.png` (rendered as image) or typeset directly in poster software
- **Message**: Key cohort characteristics at a glance
- **Caption**: *Cohort characteristics and ST analysis summary for the ≥2-timepoint and ≥3-timepoint sub-cohorts.*

> [!TIP]
> **Recommended final poster layout (5 panels):** Flow diagram → ST distribution bar → Stability comparison → ASB vs UTI stacked bars → Trajectory plot. Drop the summary table if space is tight.

---

## STEP 4 — New Analyses Needed

### 4A. Reconcile abstract numbers (CRITICAL)
- Re-run `get_stratified_stats.R` and capture output to a file
- Compare pair-level stability (abstract uses this) vs participant-level stability
- Decide which definition to use on the poster and apply consistently
- Verify the "33 UTI typeable isolates", "ST6 (n=8) in UTI" claims

### 4B. UTI-only ST distribution for ≥2-TP cohort
- Filter `mlst_with_meta.csv` × `status_map.csv` to UTI episodes only
- Count STs in UTI episodes
- **Already partially available** in [abstract_stats_output.md](file:///Users/Aamir/Desktop/rUTIs/abstract_stats_output.md) Section 4

### 4C. ≥3-TP residents with UTI: ST stability breakdown
- The abstract claims "12 residents with ≥3 TPs experienced UTI; 8 stable, 4 switched"
- This is **not** in any existing output file
- Needs a small new query on `mlst_with_meta.csv` + `status_map.csv`
- Quick to do — filter ≥3-TP participants who have at least one UTI episode, then check ST consistency

### 4D. Trajectory plot data preparation
- Select illustrative participants for Figure 5
- Prioritize: (a) participants with ≥3 TPs and stable ST, (b) participants who switched, (c) participants with UTI
- Build a tidy data frame: Participant × Timepoint × ST × Infection_Status

### 4E. Flow diagram data
- Count totals at each step: YELLOW cohort → ≥2 TPs → with ST data → typeable → analysis set
- Can be computed from `status_map.csv` and `abstract_stats_output.md`

---

## STEP 5 — Execution Plan

### A. Figures reusable immediately
1. [st_distribution_stacked.png](file:///Users/Aamir/Desktop/rUTIs/plots/epidemiology/st_distribution_stacked.png) — ASB vs UTI ST comparison (Fig 4)
2. [st_risk_plot.png](file:///Users/Aamir/Desktop/rUTIs/results/lineage/st_risk_plot.png) — UTI risk by ST (optional supplementary)

### B. Figures to regenerate from existing scripts
1. **Top STs bar chart** — re-run `07_explore_MLST.R` after filtering `mlst_with_meta.csv` to ≥2-TP cohort, excluding "-"
2. **Swimmer plot** — adapt `15_longitudinal_patterns.R` to overlay ST identity as color

### C. Figures to create from scratch
1. **CONSORT-style flow diagram** — new R script using `DiagrammeR` or `ggplot2` boxes
2. **Within-host stability comparison plot** — new script: stacked bars or donuts for ≥2-TP vs ≥3-TP
3. **Selected trajectory plot** — new script: curated participant swimmer/timeline

### D. Scripts to create or edit

| # | Script | Action | Purpose |
|---|--------|--------|---------|
| 1 | `poster_01_reconcile_numbers.R` | **New** | Re-run stratified stats, capture to file, verify abstract claims |
| 2 | `poster_02_flow_diagram.R` | **New** | Generate CONSORT flow diagram |
| 3 | `poster_03_st_bar.R` | **New** | Generate filtered top-ST bar chart for ≥2-TP cohort |
| 4 | `poster_04_stability.R` | **New** | Generate stability comparison figure (pair-level, ≥2-TP vs ≥3-TP) |
| 5 | `poster_05_trajectories.R` | **New** | Generate selected participant trajectory plot |
| 6 | `poster_06_asb_vs_uti.R` | **New** (optional) | Poster-quality version of ASB vs UTI ST stacked bar |

### E. Order of Operations
1. Run `poster_01_reconcile_numbers.R` → verify all abstract numbers
2. Run `poster_02_flow_diagram.R` → Figure 1
3. Run `poster_03_st_bar.R` → Figure 2
4. Run `poster_04_stability.R` → Figure 3
5. Run `poster_06_asb_vs_uti.R` → Figure 4 (or reuse existing)
6. Run `poster_05_trajectories.R` → Figure 5

---

## STEP 6 — Deliverables

### Ranked Figure List

| Rank | Figure | Priority | Status | Action |
|------|--------|----------|--------|--------|
| 1 | Flow diagram | Must-have | ❌ Not yet | Create `poster_02_flow_diagram.R` |
| 2 | Overall ST distribution | Must-have | ⚠️ Needs redo | Create `poster_03_st_bar.R` |
| 3 | Within-host stability | Must-have | ❌ Not yet | Create `poster_04_stability.R` |
| 4 | ASB vs UTI ST comparison | Good-to-have | ✅ Reusable | Polish existing or reuse as-is |
| 5 | Participant trajectories | Good-to-have | ❌ Not yet | Create `poster_05_trajectories.R` |
| 6 | Summary table | Optional | ❌ Not yet | Create in poster software |

### Poster Storyline (left-to-right reading order)

1. **Background** → E. coli lineage diversity in NH residents is poorly characterized; unclear if specific STs cause UTI
2. **Methods** → YELLOW RoUTIne cohort, 18-month follow-up, WGS + MLST → **Flow diagram** (Fig 1)
3. **Results Panel 1** → E. coli populations are diverse but dominated by globally disseminated lineages → **ST distribution** (Fig 2)
4. **Results Panel 2** → Within-host colonization is remarkably stable (~80% pair-level concordance) → **Stability figure** (Fig 3)
5. **Results Panel 3** → ST composition differs qualitatively between ASB and UTI → **ASB vs UTI** (Fig 4)
6. **Results Panel 4** → Individual trajectories illustrate stable colonization with occasional switching → **Trajectories** (Fig 5)
7. **Conclusions** → Stability suggests symptomatic UTI typically arises from resident commensal strains, not novel introductions. Future work: characterize untypeable isolates, explore ST–status associations with larger UTI sample

### Final Recommended 3–5 Visuals

> **If limited to 3**: Flow diagram + ST distribution + Stability comparison
>
> **If 4 panels**: Add ASB vs UTI stacked bars
>
> **If 5 panels**: Add selected trajectories

---

## Verification Plan

### Automated
- Run `poster_01_reconcile_numbers.R` and compare output against abstract claims line-by-line
- Each figure script should print the exact N to console so figures can be verified

### Manual
- Visually inspect each generated figure for readability at poster print size (A0)
- Cross-check figure labels against the reconciled numbers
- Verify color schemes are colorblind-friendly

---

## Proposed Changes

### Scripts (all NEW files in project root)

#### [NEW] [poster_01_reconcile_numbers.R](file:///Users/Aamir/Desktop/rUTIs/poster_01_reconcile_numbers.R)
Reconcile abstract claims against current data. Outputs a line-by-line verification report.

#### [NEW] [poster_02_flow_diagram.R](file:///Users/Aamir/Desktop/rUTIs/poster_02_flow_diagram.R)
Generate CONSORT-style flow diagram for the poster.

#### [NEW] [poster_03_st_bar.R](file:///Users/Aamir/Desktop/rUTIs/poster_03_st_bar.R)
Generate filtered top-ST bar chart (≥2-TP cohort, typeable only).

#### [NEW] [poster_04_stability.R](file:///Users/Aamir/Desktop/rUTIs/poster_04_stability.R)
Generate within-host stability comparison (pair-level, ≥2-TP vs ≥3-TP).

#### [NEW] [poster_05_trajectories.R](file:///Users/Aamir/Desktop/rUTIs/poster_05_trajectories.R)
Generate curated participant trajectory plot.

#### [NEW] poster_figures/ directory
All poster figures will be output to `poster_figures/` to keep them separate from pipeline outputs.

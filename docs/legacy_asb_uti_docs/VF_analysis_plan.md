# VF Analysis Plan: Directions, Figures, and Models

## A) Ranked Next-Step Directions

### 1. iro Operon Deep Dive — The Strongest Candidate (★★★★★)

**Hypothesis**: Salmochelin (iroB–N) carriage distinguishes UTI-causing from ASB-colonizing *E. coli*.

| Field | Detail |
|-------|--------|
| Biological plausibility | ★★★★★ — Salmochelin is a stealth siderophore that evades lipocalin-2; well-established uropathogenic virulence factor |
| Statistical feasibility | ★★★★ — 5-gene operon behaves as one unit (all-or-nothing), so can be scored as a single binary variable |
| Confounding risk | Medium — Must control for ST (iro is phylogenetically structured) and for Uricult timepoint |
| Inputs | `vf_analysis_ready.csv` (iroB–N columns), `status_map.csv`, MLST results if available |
| Supportive signal | iro+ isolates significantly overrepresented in UTI vs ASB, even after adjusting for ST or within-person |
| Key plot | Mosaic plot or grouped bar: iro+ proportion by status; within-person paired slopegraph for ASB→UTI transitions |

---

### 2. Within-Person Paired VF Comparison (ASB vs UTI from Same Individual) (★★★★★)

**Hypothesis**: When the same participant has both ASB and UTI episodes, VF profiles differ at the UTI event.

| Field | Detail |
|-------|--------|
| Biological plausibility | ★★★★★ — The gold standard: within-person controls eliminate between-person confounding |
| Statistical feasibility | ★★★ — Only ~12 participants contribute ASB→UTI pairs; paired tests (McNemar, signed rank) |
| Confounding risk | Low — self-controlled; Uricult timing remains a concern |
| Inputs | `vf_longitudinal_transitions.csv` (ASB→UTI subset + ASB→ASB for comparison) |
| Supportive signal | VF genes gained more often in ASB→UTI than ASB→ASB transitions |
| Key plot | Paired dot plots (one dot per participant, ASB vs UTI VF count connected by line); gained-gene frequency barplot |

---

### 3. VF Burden as Total Score — Participant-Level Comparison (★★★★)

**Hypothesis**: UTI-associated episodes carry more VF genes overall than ASB episodes.

| Field | Detail |
|-------|--------|
| Biological plausibility | ★★★ — Mixed evidence in literature; individual genes may matter more than total count |
| Statistical feasibility | ★★★★ — Simple summary; aggregate to participant-level mean to avoid repeats |
| Confounding risk | Medium — Must show participant-aggregated version alongside episode-level |
| Inputs | `vf_analysis_ready.csv` |
| Supportive signal | Mean or median VF count higher for UTI (currently: not significant) |
| Key plot | Raincloud/violin plot by status; participant-level beeswarm with repeated-measure lines |

---

### 4. Category-Level VF Composition by Status (★★★★)

**Hypothesis**: Even if total VF burden is similar, the *composition* (relative proportion of iron/adhesion/toxin genes) differs by status.

| Field | Detail |
|-------|--------|
| Biological plausibility | ★★★★ — Functional group comparison is biologically interpretable |
| Statistical feasibility | ★★★★ — Category counts from `gene_map.csv`; stacked bar or radar chart |
| Confounding risk | Low for descriptive; medium for inferential |
| Inputs | `vf_analysis_ready.csv` (cat_ columns), `gene_map.csv` |
| Supportive signal | Iron acquisition or toxin proportion elevated in UTI relative to ASB |
| Key plot | Stacked bar chart by status; heatmap of category × status |

---

### 5. VF Profile Clustering / Dimensionality Reduction (★★★)

**Hypothesis**: Distinct VF profile types exist that correlate with clinical outcome.

| Field | Detail |
|-------|--------|
| Biological plausibility | ★★★ — ExPEC pathotype concept (well-known VF combinations) |
| Statistical feasibility | ★★★ — PCA or UMAP on binary VF matrix; color by status; test cluster-status association |
| Confounding risk | Medium — ST will dominate the first axes; must disentangle from status |
| Inputs | `vf_pa_all.csv`, status labels |
| Supportive signal | UTI episodes cluster together or in specific PCA regions |
| Key plot | PCA/UMAP scatter with status colors; loadings plot for top contributing genes |

---

### 6. Timepoint-Controlled Enrichment (Uricult-Only ASB vs UTI) (★★★★)

**Hypothesis**: If we restrict to Uricult timepoint only, do ASB and UTI differ in VF?

| Field | Detail |
|-------|--------|
| Biological plausibility | ★★★ — Removes timepoint confound but reduces sample size to 17 at Uricult |
| Statistical feasibility | ★★ — Only 16 UTI + 0-1 ASB at Uricult; may not be feasible |
| Confounding risk | None (same timepoint) — but only if ASB cases exist at Uricult |
| Inputs | `vf_analysis_ready.csv` filtered to `tp_lab == "Uricult"` |
| Supportive signal | Same VF differences persist when timepoint is held constant |
| Key plot | Side-by-side with full cohort comparison to show whether effect changes after restriction |

---

### 7. Longitudinal VF Stability by Microbiological Context (★★★)

**Hypothesis**: VF profiles change more when the ST changes (strain replacement) vs when the same ST persists.

| Field | Detail |
|-------|--------|
| Biological plausibility | ★★★★ — Strain replacement vs within-strain evolution is the key distinction |
| Statistical feasibility | ★★ — Requires linked MLST data for transitions (may or may not exist for all pairs) |
| Confounding risk | Low — purely descriptive |
| Inputs | `vf_longitudinal_transitions.csv`, MLST data, strain comparison results |
| Supportive signal | Jaccard drops significantly when ST changes |
| Key plot | Jaccard boxplot stratified by same-ST vs different-ST |

---

### 8. Gene Co-Occurrence Network (★★)

**Hypothesis**: VF genes form co-occurring modules (e.g., pathogenicity island blocks) that differ by status.

| Field | Detail |
|-------|--------|
| Biological plausibility | ★★★★ — Known PAIs co-carry genes |
| Statistical feasibility | ★★ — Requires careful correlation analysis; φ-coefficient matrix |
| Confounding risk | Medium — Co-occurrence may reflect clonal structure |
| Inputs | `vf_pa_all.csv` (binary matrix only) |
| Supportive signal | Gene modules differentially present by status |
| Key plot | Network graph or clustered heatmap of gene-gene co-occurrence |

---

## B) Figure Set Blueprint

### Figure 1: VF Burden Distribution by Clinical Status

**What this proves**: Overall VF arsenal size does not significantly differ between ASB and UTI episodes.

| Panel | Content | Cohort | x | y | Annotations |
|-------|---------|--------|---|---|------------|
| A | Raincloud plot (violin + jitter + boxplot) | All 183 episodes | Infection_Status | vf_count_total | n per group, median lines, Wilcoxon p (labelled exploratory) |
| B | Participant-aggregated beeswarm | One mean per participant per status (avoid pseudo-replication) | Status | Mean VF count | n participants per status, lines connecting participants with both ASB+UTI |
| C | Faceted by timepoint (T0, T1, T2, Uricult) | All episodes | Status within facet | vf_count | n per cell |

---

### Figure 2: Top Differentially Prevalent VF Genes (ASB vs UTI)

**What this proves**: A small number of genes (astA, iro operon, pic) show notable prevalence differences, but none survive multiple testing correction.

| Panel | Content | Cohort | Annotations |
|-------|---------|--------|------------|
| A | Lollipop/dot plot: top 20 genes by |Δ(UTI−ASB)| | ASB+UTI only (n=152) | Prevalence %, n/N, Δ with 95% CI |
| B | Forest plot: Fisher exact ORs + 95% CI for top 15 genes | ASB+UTI | OR + CI, dotted line at OR=1, BH-adjusted p annotation |
| C | Highlighted panel for iro operon genes showing co-inheritance | All statuses | Correlation matrix or co-occurrence heatmap of iro genes |

---

### Figure 3: VF Category Composition by Status

**What this proves**: The functional composition of VF arsenals (iron/adhesion/toxin distribution) is broadly similar across statuses, with subtle iron-acquisition enrichment in UTI.

| Panel | Content | Cohort |
|-------|---------|--------|
| A | Stacked bar chart: mean category counts by status | All with status (n=183) |
| B | Heatmap: category × status, color = median gene count | Same |
| C | Per-category violin or boxplot (one per category, faceted), by status | Same |

---

### Figure 4: Within-Person VF Stability (Longitudinal)

**What this proves**: VF profiles are remarkably stable within individuals over time; limited VF change accompanies ASB→UTI transitions.

| Panel | Content | Cohort |
|-------|---------|--------|
| A | Histogram of Jaccard similarity, by transition type (ASB→ASB vs ASB→UTI vs other) | 96 transitions |
| B | Paired slopegraph: VF count at t_from vs t_to, colored by transition type | Same |
| C | Barplot: top 10 most frequently gained genes in ASB→UTI (with gene name + count/12) | ASB→UTI only (n=12) |
| D | Heatmap: genes gained/lost per ASB→UTI transition (12 rows × top variable genes cols) | ASB→UTI only |

---

### Figure 5: Confounding Check — Status × Timepoint and Adjusted Analysis

**What this proves**: UTI status is confounded with Uricult timepoint; the VF story must be interpreted within this constraint.

| Panel | Content |
|-------|---------|
| A | Tile/balloon plot: Status × Timepoint cross-tabulation (all 183 episodes) |
| B | VF burden by status at scheduled visits only (T0+T1+T2, excluding Uricult): ASB vs Negative |
| C | VF burden at Uricult only: UTI vs Negative (if any Neg exist at Uricult) |
| D | Within-person comparison: participants with both scheduled + Uricult data (paired VF count) |

---

## C) Model Plan

### Model 1: Mixed-Effects Logistic Regression per Gene (Primary)

```
Outcome ~ gene_present + (1 | Participant_id)
```

- **Outcome**: Binary, 1 = UTI, 0 = ASB (exclude Negative)
- **Predictor**: Gene presence (one model per gene, or gene-set)
- **Random effect**: Participant intercept (accounts for repeated measures)
- **Confounders to explore**: Timepoint (but collinear with outcome); Batch if available
- **Limitation**: UTI n=16 is very small for GLMM convergence; expect singularity warnings
- **Report**: OR + 95% CI + p; flag any singular fits; compare with Fisher as sensitivity check

### Model 2: VF Score as Continuous Outcome (GEE)

```
vf_count ~ Infection_Status + (1 | Participant_id), family = gaussian, corstr = "exchangeable"
```

- **Outcome**: VF count (continuous)
- **GEE** with exchangeable correlation handles repeated measures robustly without convergence issues
- **Variants**: Replace `vf_count` with iron-acquisition count, or iro binary
- **Report**: β ± SE, 95% CI, p; interpret as adjusted mean difference in VF count

### Model 3: Dimensionality-Reduced VF Profile (PCA + Mixed Model)

```
PC1 ~ Infection_Status + (1 | Participant_id)
```

- Compute PCA on binary VF matrix
- Test whether first 2-3 PCs differ by Infection_Status using mixed model
- **Advantage**: Reduces 164 tests to 2-3; avoids multiple testing burden
- **Report**: Variance explained, loadings of top genes, association with status

---

## D) Script File Paths (to be created)

| Script | Purpose | Output dir |
|--------|---------|-----------|
| `scripts/build_vf_analysis_table.R` | Build unified analysis table | `results/vf/` |
| `scripts/visualize_vf_01_burden.R` | Figure 1 | `plots/vf/` |
| `scripts/visualize_vf_02_gene_prevalence.R` | Figure 2 | `plots/vf/` |
| `scripts/visualize_vf_03_category_profiles.R` | Figure 3 | `plots/vf/` |
| `scripts/visualize_vf_04_longitudinal_stability.R` | Figure 4 | `plots/vf/` |
| `scripts/visualize_vf_05_confounding_checks.R` | Figure 5 | `plots/vf/` |

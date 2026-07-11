# Stratified VF Analysis: ≥2, ≥3, ≥4 Timepoints

> All numbers from [compute_vf_stratified_by_depth.R](file:///Users/Aamir/Desktop/rUTIs/compute_vf_stratified_by_depth.R)

---

## Cohort Sizes Across Strata

| | ≥2 Timepoints | ≥3 Timepoints | ≥4 Timepoints |
|---|---|---|---|
| **Participants** | 83 | 8 | 5 |
| **Total episodes** | 179 | 29 | 20 |
| ASB episodes | 135 | 18 | 13 |
| UTI episodes | 13 | 7 | 4 |
| Negative episodes | 31 | 4 | 3 |

> [!NOTE]
> The ≥3 and ≥4 cohorts are dramatically smaller — only 8 and 5 participants respectively — so all numbers from those strata have wide confidence intervals and should be interpreted cautiously.

---

## VF Burden by Status × Cohort

| Status | Metric | ≥2 Timepoints | ≥3 Timepoints | ≥4 Timepoints |
|--------|--------|--------------|--------------|--------------|
| **ASB** | n | 135 | 18 | 13 |
| | Median (IQR) | 80 (70–92) | 76 (68–84) | 77 (71–84) |
| | Mean ± SD | 79.6 ± 16.2 | 75.4 ± 12.8 | 74.8 ± 13.0 |
| **UTI** | n | 13 | 7 | 4 |
| | Median (IQR) | 84 (71–99) | 77 (74–98) | 86 (76–99) |
| | Mean ± SD | 82.8 ± 21.3 | 85.4 ± 16.9 | 88.8 ± 18.6 |
| **Negative** | n | 31 | 4 | 3 |
| | Median (IQR) | 88 (77–100) | 97 (94–98) | 97 (90–97) |
| | Mean ± SD | 85.7 ± 20.1 | 94.2 ± 7.6 | 92.3 ± 8.1 |

**What this means**: The overall VF burden (total number of virulence genes per isolate) stays comparable across ASB and UTI at all timepoint depths. As we narrow to the most deeply followed participants (≥4), UTI episodes trend slightly higher (median 86 vs 77 for ASB), but the sample size is too small (n=4 UTI) to draw conclusions.

---

## Top Differentiating Genes: UTI vs ASB

### ≥2 Timepoints (n_ASB=135, n_UTI=13)

| Gene | Category | ASB% | UTI% | Δ | Fisher p |
|------|----------|------|------|---|---------|
| sfaF | Adhesion/Fimbriae | 25.2% | 46.2% | +21.0 | 0.114 |
| iroB–N (×5) | Iron acquisition | 41.5% | 61.5% | +20.0 | — |
| sfaE/G | Adhesion/Fimbriae | 26.7% | 46.2% | +19.5 | — |
| pic | Unassigned | 18.5% | 0.0% | −18.5 | 0.126 |
| astA | Toxins | 6.7% | 23.1% | +16.4 | 0.074 |

### ≥3 Timepoints (n_ASB=18, n_UTI=7)

| Gene | Category | ASB% | UTI% | Δ | Fisher p |
|------|----------|------|------|---|---------|
| iroB–N (×5) | Iron acquisition | 16.7% | 42.9% | **+26.2** | 0.299 |
| sfaB/C/D/E/F | Adhesion/Fimbriae | 11.1% | 28.6% | +17.5 | — |
| espR4 | Unassigned | 0.0% | 14.3% | +14.3 | 0.280 |

### ≥4 Timepoints (n_ASB=13, n_UTI=4)

| Gene | Category | ASB% | UTI% | Δ | Fisher p |
|------|----------|------|------|---|---------|
| **iroB–N (×5)** | **Iron acquisition** | **15.4%** | **50.0%** | **+34.6** | **0.219** |
| irp2, ybtA/P/Q/S/X | Iron acquisition | 76.9% | 100.0% | +23.1 | — |
| espR4 | Unassigned | 0.0% | 25.0% | +25.0 | — |

> [!IMPORTANT]
> **The iro operon (salmochelin) prevalence gap grows with deeper follow-up**: +15pp at all participants → +20pp at ≥2tp → +26pp at ≥3tp → **+35pp at ≥4tp**. This is the most consistent signal across all strata, though none reach significance given the small UTI n's.

---

## Longitudinal VF Dynamics × Cohort

### Overview

| Metric | ≥2 Timepoints | ≥3 Timepoints | ≥4 Timepoints |
|--------|--------------|--------------|--------------|
| Participants | 83 | 8 | 5 |
| Total transitions | 96 | 21 | 15 |
| Median Jaccard | 1.000 | 1.000 | 1.000 |
| % zero VF change | **72.9%** | **66.7%** | **53.3%** |

**What this means**: More VF changes are detected with deeper follow-up. At ≥2 timepoints, 73% of transitions show no change at all; by ≥4 timepoints, this drops to 53%. This makes biological sense — more sampling occasions means more opportunity to capture genuine VF variation.

### By Transition Type

| Transition | ≥2tp: n (Jac, %noΔ) | ≥3tp: n (Jac, %noΔ) | ≥4tp: n (Jac, %noΔ) |
|-----------|---------------------|---------------------|---------------------|
| ASB→ASB | 61 (1.000, 82%) | 10 (1.000, 70%) | 8 (1.000, 62%) |
| ASB→UTI | 12 (1.000, 58%) | 6 (1.000, 67%) | 4 (0.942, 50%) |
| Neg→Neg | 9 (0.989, 44%) | — | — |
| Neg→ASB | 7 (1.000, 57%) | 2 (0.717, 50%) | 2 (0.717, 50%) |
| ASB→Neg | 6 (1.000, 67%) | 2 (0.717, 50%) | 1 (0.433, 0%) |
| Neg→UTI | 1 (1.000, 100%) | 1 (1.000, 100%) | — |

**Key observation**: In the ≥4 timepoint cohort, the **ASB→UTI Jaccard drops to 0.942** — meaning that among deeply followed participants who transitioned to UTI, there is measurable VF gene content change. This contrasts with the ≥2tp cohort where the median Jaccard for ASB→UTI is 1.000 (mostly no change).

### Most Gained Genes in ASB→UTI Transitions

| | ≥2tp (12 transitions) | ≥3tp (6 transitions) | ≥4tp (4 transitions) |
|---|---|---|---|
| **Top gained** | iroB/C/D/E/N (4 each) | afaA/B/C, draD/P (1 each) | afaA/B/C, draD/P (1 each) |
| **2nd gained** | kpsD/M (3 each) | espL1, sfaC/D (1 each) | espL1, sfaC/D (1 each) |
| **3rd gained** | chuA/T/X (2 each) | — | — |

**What this means**: In the broader ≥2tp cohort, iron acquisition (iro) and capsule (kps) genes dominate. In the narrower ≥3/≥4tp cohorts, adhesion genes (afa, dra) and S fimbriae (sfa) appear — suggesting that different VF mechanisms may be gained in different individuals transitioning to UTI. The small numbers make it impossible to identify a single universal mechanism.

---

## Summary: What Changes With Deeper Follow-Up?

| Finding | ≥2tp | ≥3tp | ≥4tp | Trend |
|---------|------|------|------|-------|
| ASB vs UTI burden similar? | Yes | Yes | Yes (slight UTI ↑) | **Stable** |
| iro operon UTI–ASB gap | +20pp | +26pp | **+35pp** | **↑ Growing** |
| Any gene reaches significance? | astA p=0.07 | None | None | Underpowered |
| % transitions with zero change | 73% | 67% | 53% | **↓ More change seen** |
| ASB→UTI Jaccard | 1.000 | 1.000 | 0.942 | **↓ Slightly less stable** |

> [!IMPORTANT]
> **The core story strengthens with deeper follow-up**: The iro operon prevalence gap between UTI and ASB *widens*, and longitudinal VF stability *decreases* — suggesting that among the most-followed participants, genuine VF dynamics are present. However, the ≥3tp and ≥4tp cohorts are very small (8 and 5 participants), so this trend needs confirmation in larger datasets.

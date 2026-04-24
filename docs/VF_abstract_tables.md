# VF Abstract-Ready Tables

> All numbers computed from anchor files by [compute_vf_abstract_stats.R](file:///Users/Aamir/Desktop/rUTIs/compute_vf_abstract_stats.R). Source CSVs cited in [VF_traceability_log.md](file:///Users/Aamir/Desktop/rUTIs/docs/VF_traceability_log.md).

---

## Table 1: VF Analysis Cohort

| Metric | Value |
|--------|-------|
| Participants with VF data | 87 |
| Participant×Timepoint rows with VF | 183 |
| All *E. coli* | Yes (382/382 assemblies) |
| Clinical episodes in `status_map.csv` | 276 |
| Clinical episodes WITH matched VF data | 183 (66.3%) |
| Clinical episodes WITHOUT VF data | 93 (33.7%) |

### Status Breakdown (VF-available episodes only)

| Infection Status | n episodes | n participants |
|-----------------|-----------|---------------|
| ASB | 136 | 75 |
| UTI | 16 | 16 |
| Negative | 31 | 21 |
| **Total** | **183** | **87** |

### Timepoint Coverage

| Timepoints per participant | n participants |
|---------------------------|---------------|
| ≥2 | 83 |
| ≥3 | 8 |
| ≥4 | 5 |

### Episodes per Timepoint

| Timepoint | ASB | UTI | Negative | Total |
|-----------|-----|-----|----------|-------|
| T0 | 64 | 0 | 16 | 80 |
| T1 | 64 | 0 | 13 | 77 |
| T2 | 8 | 0 | 1 | 9 |
| Uricult | 0 | 16 | 1 | 17 |

> [!IMPORTANT]
> All 16 UTI episodes occur at the Uricult timepoint. T0/T1 are exclusively ASB + Negative. This is a design feature of the study: Uricult samples are collected when UTI is suspected.

> Source: [vf_analysis_ready.csv](file:///Users/Aamir/Desktop/rUTIs/results/vf/vf_analysis_ready.csv), [vf_burden_by_status.csv](file:///Users/Aamir/Desktop/rUTIs/results/vf/vf_burden_by_status.csv)

---

## Table 2: VF Burden by Clinical Status

| Status | n | Mean VF (SD) | Median VF (IQR) | Range |
|--------|---|-------------|-----------------|-------|
| ASB | 136 | 79.3 (16.5) | 80 (70–92) | 32–109 |
| UTI | 16 | 81.3 (19.6) | 80.5 (70–98) | 32–112 |
| Negative | 31 | 85.7 (20.1) | 88 (77–100) | 43–121 |

> VF count = number of distinct VFDB genes detected per participant-timepoint (out of 164 possible).
> Source: [vf_burden_by_status.csv](file:///Users/Aamir/Desktop/rUTIs/results/vf/vf_burden_by_status.csv)

---

## Table 3: Top VF Genes by Status (ranked by |Δ(UTI−ASB)|)

| Gene | Category | ASB (n/136, %) | UTI (n/16, %) | Neg (n/31, %) | Δ(UTI−ASB) |
|------|----------|---------------|--------------|--------------|------------|
| astA | Toxins | 9/136 (6.6%) | 4/16 (25.0%) | 3/31 (9.7%) | +18.4 |
| pic | Unassigned | 25/136 (18.4%) | 0/16 (0.0%) | 2/31 (6.5%) | −18.4 |
| iroN | Iron acquisition | 56/136 (41.2%) | 9/16 (56.2%) | 6/31 (19.4%) | +15.0 |
| iroE | Iron acquisition | 56/136 (41.2%) | 9/16 (56.2%) | 6/31 (19.4%) | +15.0 |
| iroD | Iron acquisition | 56/136 (41.2%) | 9/16 (56.2%) | 6/31 (19.4%) | +15.0 |
| iroC | Iron acquisition | 56/136 (41.2%) | 9/16 (56.2%) | 6/31 (19.4%) | +15.0 |
| iroB | Iron acquisition | 56/136 (41.2%) | 9/16 (56.2%) | 6/31 (19.4%) | +15.0 |
| papE | Adhesion/Fimbriae | 28/136 (20.6%) | 1/16 (6.2%) | 4/31 (12.9%) | −14.4 |
| vat | Toxins | 76/136 (55.9%) | 11/16 (68.8%) | 10/31 (32.3%) | +12.9 |
| sfaF | Adhesion/Fimbriae | 34/136 (25.0%) | 6/16 (37.5%) | 2/31 (6.5%) | +12.5 |

> Fisher exact test (exploratory, unadjusted): astA OR=4.63, p=0.033; pic OR=0, p=0.075. No genes significant after BH correction.
> Source: [vf_gene_prevalence_by_status.csv](file:///Users/Aamir/Desktop/rUTIs/results/vf/vf_gene_prevalence_by_status.csv), [vf_gene_enrichment_UTI_vs_ASB.csv](file:///Users/Aamir/Desktop/rUTIs/results/vf/vf_gene_enrichment_UTI_vs_ASB.csv)

---

## Table 4: Within-Person VF Dynamics (Longitudinal Transitions)

### 4a: Transition Counts

| Transition | n | n participants | Median Jaccard | % No VF Change |
|-----------|---|---------------|---------------|----------------|
| ASB→ASB | 61 | 58 | 1.000 | 82.0% |
| ASB→UTI | 12 | 12 | 1.000 | 58.3% |
| Negative→Negative | 9 | 9 | 0.989 | 44.4% |
| Negative→ASB | 7 | 7 | 1.000 | 57.1% |
| ASB→Negative | 6 | 6 | 1.000 | 66.7% |
| Negative→UTI | 1 | 1 | 1.000 | 100.0% |
| **Total** | **96** | **83** | **1.000** | **72.9%** |

### 4b: VF Changes per Transition

| Transition | Median Gained | Mean Gained | Median Lost | Mean Lost | Median Stable |
|-----------|-------------|-----------|------------|---------|-------------|
| ASB→ASB | 0 | 2.3 | 0 | 0.9 | 77 |
| ASB→UTI | 0 | 9.6 | 0 | 3.2 | 74.5 |
| Negative→Negative | 0 | 12.8 | 0 | 3.8 | 86 |
| Negative→ASB | 0 | 9.0 | 0 | 11.9 | 79 |
| ASB→Negative | 0 | 13.0 | 0 | 1.3 | 89.5 |

### 4c: Most Commonly Gained Genes in ASB→UTI Transitions (5/12 had gains)

| Gene | Times Gained (out of 12 ASB→UTI) | Category |
|------|--------------------------------|----------|
| iroB | 4 | Iron acquisition |
| iroC | 4 | Iron acquisition |
| iroD | 4 | Iron acquisition |
| iroE | 4 | Iron acquisition |
| iroN | 4 | Iron acquisition |
| kpsD | 3 | Capsule/Surface |
| kpsM | 3 | Capsule/Surface |
| chuA | 2 | Iron acquisition |
| chuT | 2 | Iron acquisition |
| chuX | 2 | Iron acquisition |

> Source: [vf_longitudinal_transitions.csv](file:///Users/Aamir/Desktop/rUTIs/results/vf/vf_longitudinal_transitions.csv), [vf_transition_summary_by_type.csv](file:///Users/Aamir/Desktop/rUTIs/results/vf/vf_transition_summary_by_type.csv)

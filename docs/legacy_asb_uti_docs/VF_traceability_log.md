# VF Traceability Log

> Every statistic reported in the abstract tables and draft is traced here to its source file, join logic, filters, and code origin.

## Methodology

- **Primary anchor files**: `results/vf/vf_pa_all.csv`, `status_map.csv`, `results/vf/gene_map.csv`
- **Computation script**: [compute_vf_abstract_stats.R](file:///Users/Aamir/Desktop/rUTIs/compute_vf_abstract_stats.R)
- **All joins documented in**: [VF_merge_diagnostics.md](file:///Users/Aamir/Desktop/rUTIs/docs/VF_merge_diagnostics.md)

---

## Statistic: Cohort Size (87 participants, 183 rows)

| Field | Value |
|-------|-------|
| File(s) | `results/vf/vf_pa_all.csv` |
| Code | `n_distinct(vf_pa$Participant_id)` → 87; `nrow(vf_pa)` → 183 |
| Join | None (direct from anchor) |
| Filters | None |
| Type | **Primary** |
| Output CSV | `results/vf/vf_analysis_ready.csv` (183 rows) |

---

## Statistic: Status Breakdown (136 ASB, 16 UTI, 31 Negative)

| Field | Value |
|-------|-------|
| File(s) | `results/vf/vf_pa_all.csv` + `status_map.csv` |
| Join | `left_join(vf_pa, status_map, by = c("Participant_id", "tp_lab"))` |
| Match rate | 183/183 (100%) VF rows matched; 0 unmatched |
| Filters | `!is.na(Infection_Status)` → 183 rows (all matched) |
| Code | `table(vf_ready$Infection_Status)` |
| Type | **Primary** |
| Output CSV | `results/vf/vf_analysis_ready.csv` column `Infection_Status` |

---

## Statistic: VF Burden — ASB median=80 (IQR 70–92)

| Field | Value |
|-------|-------|
| File(s) | `results/vf/vf_analysis_ready.csv` |
| Code | Section 3, `summarise(median_vf = median(vf_count_total), ...)` grouped by `Infection_Status` |
| Filters | `Infection_Status == "ASB"` → n=136 |
| Definition | `vf_count_total = rowSums(across(all_of(gene_cols)))` — count of all 164 VF genes present |
| Type | **Primary** |
| Output CSV | `results/vf/vf_burden_by_status.csv` row 1 |

---

## Statistic: VF Burden — UTI median=80.5 (IQR 70–98)

| Field | Value |
|-------|-------|
| File(s) | `results/vf/vf_analysis_ready.csv` |
| Code | Same as above, `Infection_Status == "UTI"` → n=16 |
| Type | **Primary** |
| Output CSV | `results/vf/vf_burden_by_status.csv` row 3 |

---

## Statistic: astA prevalence 25.0% UTI vs 6.6% ASB, OR=4.63, p=0.033

| Field | Value |
|-------|-------|
| File(s) | `results/vf/vf_analysis_ready.csv` |
| Code — prevalence | Section C2, per-gene `sum(.data[[g]] > 0) / n()` grouped by `Infection_Status` |
| Code — Fisher | Section C3, `fisher.test(table(gene_present, is_uti))` |
| Filters | `Infection_Status %in% c("ASB", "UTI")` → n=152 |
| BH correction | `p.adjust(p_value, method = "BH")` → p_adj = 1.000 (not significant) |
| Type | **Primary** |
| Output CSV | `results/vf/vf_gene_prevalence_by_status.csv` (gene=astA), `results/vf/vf_gene_enrichment_UTI_vs_ASB.csv` (gene=astA) |

---

## Statistic: 96 transitions, 72.9% zero VF change, median Jaccard=1.000

| Field | Value |
|-------|-------|
| File(s) | `results/vf/vf_analysis_ready.csv` |
| Code | Section 4, consecutive-pair loop for participants with ≥2 timepoints |
| Join | Uses `vf_with_status` (already joined with status) |
| Filters | `tp_lab %in% c("T0","T1","T2","T3","T4","Uricult")`, participants with ≥2 timepoints |
| Ordering | T0 < T1 < T2 < T3 < T4 < Uricult |
| Pair construction | Consecutive only (t_i → t_{i+1}) |
| Jaccard | `length(intersect) / length(union)` of VF gene sets |
| Type | **Primary** |
| Output CSV | `results/vf/vf_longitudinal_transitions.csv` (96 rows) |

---

## Statistic: 12 ASB→UTI transitions, iroB/C/D/E/N gained 4/12

| Field | Value |
|-------|-------|
| File(s) | `results/vf/vf_longitudinal_transitions.csv` |
| Code | `filter(transition_type == "ASB→UTI")` → 12 rows; `genes_gained` column parsed |
| Filters | transition_type exactly "ASB→UTI" |
| Type | **Primary** |
| Output CSV | `results/vf/vf_longitudinal_transitions.csv` (filtered), `results/vf/vf_transition_summary_by_type.csv` |

---

## Statistic: Timepoint coverage (T0=80, T1=77, T2=9, Uricult=17)

| Field | Value |
|-------|-------|
| File(s) | `results/vf/vf_analysis_ready.csv` |
| Code | `table(d$tp_lab)` |
| Filters | None (all 183 rows) |
| Type | **Primary** |

---

## Statistic: All assemblies E. coli (382/382)

| Field | Value |
|-------|-------|
| File(s) | `assembly_metadata.csv` |
| Code | `cut -d',' -f13 assembly_metadata.csv | sort | uniq -c` |
| Type | **Cross-check** (not from anchor files, but confirms species assumption) |

---

## Cross-Check Sources

| Existing output | Compared against | Concordance | Discrepancy explanation |
|----------------|-----------------|-------------|----------------------|
| `results/stratified_vf_stats_table.csv` | Our burden stats | Different magnitudes | They exclude "Unassigned" category genes; we count all 164 |
| `results/vf/diff_focus_genes_UTI_vs_ASB_glmm.csv` | Our Fisher results | Both show no signal | GLMM degenerate (OR=1, p=1); Fisher more interpretable |
| `results/vf/stats_gene_level.csv` | Our gene prevalence | Consistent | Same underlying data |

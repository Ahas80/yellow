# VF Merge Diagnostics

## Source Files

| File | Rows | Unique Participants | Key Columns |
|------|------|-------------------|-------------|
| `results/vf/vf_pa_all.csv` | 183 | 87 | Participant_id, tp_lab |
| `status_map.csv` | 276 | 97 | Participant_id, Timepoint |

## Duplicate Checks

| Dataset | Duplicate join keys (Participant_id × tp_lab) |
|---------|-----------------------------------------------|
| `vf_pa_all.csv` | **0** — no duplicates |
| `status_map.csv` | **0** — no duplicates after tp normalization |

## Join Results

| Metric | Count |
|--------|-------|
| VF rows matched to status | **183** (100%) |
| VF rows NOT in status_map | **0** |
| Status rows NOT in VF data | **93** |
| → These are clinical episodes without sequenced VF data | |

### Unmatched Status Rows (clinical episodes without VF data)

These 93 clinical episodes exist in `status_map.csv` but have no corresponding VF data — either no assembly was available or the participant was not sequenced at that timepoint.

| Status | Count |
|--------|-------|
| ASB | 77 |
| Negative | 9 |
| UTI | 6 |
| Culture-positive | 1 |

## Final Analysis-Ready Dataset

| Metric | Value |
|--------|-------|
| Output file | `results/vf/vf_analysis_ready.csv` |
| Total rows | **183** |
| Rows with clinical status | **183** (100% match) |
| Rows without status | **0** |
| Unique participants | **87** |
| Species | All *E. coli* (382/382 in assembly_metadata.csv) |

### Status Breakdown in Analysis-Ready Dataset

| Infection_Status | n episodes |
|-----------------|-----------|
| ASB | 136 |
| Negative | 31 |
| UTI | 16 |

### Timepoints Present

| Timepoint | In VF data | In status_map |
|-----------|-----------|--------------|
| T0 | ✅ | ✅ |
| T1 | ✅ | ✅ |
| T2 | ✅ | ✅ |
| Uricult | ✅ | ✅ |

## Methodology Note

- **Join key normalization**: `status_map.csv` uses "Timepoint" column; `vf_pa_all.csv` uses "tp_lab". Both were normalized using regex: Uricult-like → "Uricult", digit-containing → "T{n}".
- **Deduplication in VF matrix**: The pipeline's `02_gene_presence_analysis.R` constructs `vf_pa_all.csv` by `distinct(Participant_id, tp_lab, GENE)` before pivoting — so a VF gene is marked present if detected in *any* assembler (flye or longcycler) for that participant-timepoint pair.
- **No VF data loss**: All 183 VF rows matched a clinical status; the 93 unmatched status rows represent clinical episodes for which no sequenced isolate/assembly exists.

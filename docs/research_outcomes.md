# Research Outcomes: Current UTI vs Not_UTI Analysis

**Current primary definition:** `catheter_adjusted_sns_cfu1e3_v1`  
**Primary contrast:** `UTI` vs `Not_UTI`

This document supersedes older ASB-vs-UTI summaries. Legacy ASB / UTI / Negative labels are retained only for comparability and reclassification auditing.

## Executive Summary

The current pipeline asks whether clinical UTI episodes differ from all non-UTI episodes under a catheter-aware signs/symptoms rule and lower culture-support threshold.

Current clinical denominator after classification:

- `status_map.csv`: 579 clinical episodes.
- Primary status: 12 `UTI`, 567 `Not_UTI`.
- Legacy comparison: 30 old `UTI`, 450 old `ASB`, 99 old `Negative`.
- `Not_UTI` is heterogeneous and includes bacteriuria-not-UTI, culture-negative/below-threshold episodes, and indeterminate episodes.

The main interpretation change is important: this is no longer a pure ASB-vs-UTI bacterial virulence contrast. It is a broader clinical contrast between UTI and all non-UTI episodes.

## Current Outputs To Use

- Clinical truth table: `results/clinical/status_map.csv`
- Display-ordered clinical table: `results/clinical/status_map_with_poster_tp.csv`
- Reclassification audit: `results/clinical/status_map_legacy_comparison.csv`
- CFU/provenance audits: `results/clinical/uti_cfu_threshold_audit.csv`
- Symptom-rule audit: `results/clinical/uti_symptom_rule_audit.csv`
- VF-ready table: `results/vf/vf_analysis_ready.csv`
- Binary VF-ready subset: `results/vf/vf_binary_uti_ready.csv`
- Model denominator: `results/models/model_dataset_denominator.csv`
- Key summary: `results/summary/final_key_results_summary.md`

## Figure Interpretation

Use current figures whose titles or filenames say `UTI`, `Not_UTI`, primary status, catheter-aware S&S, or `>=10^3 CFU`.

Do not use ASB-named plots or tables as primary outputs. They should either be absent from current output folders or archived under `results/legacy/old_asb_uti_outputs/` and `plots/legacy/old_asb_uti_outputs/`.

Most useful current clinical diagnostics:

- `plots/clinical/uti_reclassification_heatmap.png`
- `plots/clinical/not_uti_subgroup_by_batch_event.png`
- `plots/clinical/uti_symptom_rule_provenance.png`
- `plots/clinical/uti_cfu_threshold_provenance.png`
- `plots/clinical/waterfall_counts.png`

## Scientific Caveats

- The primary comparator `Not_UTI` is heterogeneous.
- The lower CFU threshold increases sensitivity but may reduce specificity.
- Catheter and non-catheter residents use different symptom rules.
- Missing or ambiguous collection method can affect classification.
- Uricult/event-driven samples may have different sampling bias than routine samples.
- VF and genotype-phenotype models are underpowered because the VF-ready primary UTI count is small.
- Repeated measures and lineage/ST structure remain important in all modelling.

## Current Recommendation

Interpret primary analyses as exploratory UTI-vs-Not_UTI evidence. Use `Not_UTI_subgroup` and legacy ASB/UTI/Negative outputs only as labelled descriptive or sensitivity analyses.

# Current Figure Legends

These legends describe the current primary analysis: `UTI` vs `Not_UTI` under the catheter-aware S&S plus >=10^3 CFU/mL definition. Legacy ASB-vs-UTI legends are obsolete unless explicitly used for a labelled legacy sensitivity analysis.

## F001 – Longitudinal Primary UTI Status Heatmap

Shows participants across timepoints, with cells coloured by primary `UTI_Status` (`UTI` or `Not_UTI`). This is a cohort overview, not an ASB-vs-UTI comparison.

## F002 – Primary Status Transitions

Shows transitions between `Not_UTI` and `UTI` across consecutive ordered episodes. Transition counts should be interpreted with Uricult/date-ordering caveats from the transition index.

## F003 – Catheter-Aware UTI Definition Waterfall

Shows how the denominator narrows from all clinical episodes to episodes with known collection method, culture support at the lower threshold, symptom-compatible episodes, and final primary UTI.

## F004 – Legacy-To-Primary Reclassification Heatmap

Shows how old ASB / UTI / Negative labels map into the current primary `UTI` and `Not_UTI` subgroup framework. This figure is for interpretability and denominator auditing.

## F005 – Not_UTI Subgroup Composition

Shows the heterogeneous makeup of the `Not_UTI` comparator by batch and sampling context. This helps distinguish bacteriuria-not-UTI from culture-negative/below-threshold or indeterminate rows.

## F006 – Symptom-Rule Provenance

Shows which catheter-aware symptom rule was applied. Catheter episodes use systemic symptoms only; non-catheter episodes use local urinary symptoms or flank pain plus systemic signs.

## F007 – CFU Threshold Provenance

Shows CFU threshold parsing and fallback sources. The primary definition uses >=10^3 CFU/mL where parsed, while >=10^5 is retained for legacy comparison.

## F008/F009 – VF Burden And Gene Prevalence By Primary Status

These figures compare UTI to Not_UTI and are exploratory because the UTI denominator is small. `Not_UTI_subgroup` should be checked before making biological claims.

## F010 – Genotype-Phenotype Volcano Plot

Shows exploratory association results for `UTI_binary` from `vf_analysis_ready.csv`. Odds ratios compare UTI to Not_UTI, not UTI to ASB.

## F011/F012 – Publication Longitudinal Figures

Show primary-status trajectories and candidate switch events. These should be interpreted as descriptive longitudinal evidence, not proof that VF changes cause UTI.

# Longcycler-only methods summary

## Release contract

All analytical tables, figures, models and timelines use the same selected QC-passing Longcycler cohort: **532 genome-linked episodes from 161 residents**, comprising **16 operational UTI** and **516 operational Not_UTI** episodes. The full clinical source (583 episodes; 166 residents; 18 UTI; 565 Not_UTI) is retained only as explicitly labelled attrition/QC context.

The UTI label is the current versioned **operational phenotype**. It is not presented as a reconstruction of the full published protocol. Not_UTI is not a healthy-control category.

## Analytical units

| Unit | Count | Interpretation |
|---|---:|---|
| Selected genome-linked episodes | 532 | One selected QC-passing Longcycler assembly per exact episode key. |
| Residents | 161 | Distinct residents represented in the analytical cohort. |
| Direct within-resident pairs | 893 | Every unordered pair has two selected Longcycler endpoints. |
| Adjacent transitions | 371 | Rebuilt in clinical time order from 139 residents. |
| Adjacent transitions at ≤25 SNPs | 140 | Operational direct-pair support boundary. |
| Not_UTI→UTI transitions | 9 | Focused longitudinal denominator. |
| Not_UTI→UTI at ≤25 SNPs | 5 | Direct evidence subset. |
| Mechanism casebook | 9 cases; 9 linked; 0 missing | Complete focused traceability set. |
| Near-miss audit | 17 | Separate audit rows; not operational UTI cases. |
| VFDB feature space | 227 | Binary features measured on selected episodes. |
| MLST lineage context | 493 typed; 76 ST labels | Lineage context, not direct same-strain proof. |

## Methods wording

Genome analyses were restricted to selected canonical Longcycler assemblies that passed the implemented assembly-QC screen. The exact selected assembly keys were joined to the versioned operational clinical phenotype, VFDB presence/absence calls, MLST lineage calls and content-bound direct genomic comparisons. Direct within-resident pair evidence was treated as primary: graph connectivity or ST agreement could not override a direct distance above the operational ≤25-SNP support boundary. Longitudinal transitions were rebuilt from the retained 532-episode cohort rather than filtered from an earlier timeline.

## Interpretation boundary

The release supports descriptive VF, lineage, direct-pair, population-context and longitudinal summaries. It does not support causal UTI mechanisms, antibiotic effects, demographic or host-factor claims, named-mutation conclusions, or conversion between distinct distance methods. RQ01–RQ10 form the numbered release layer.

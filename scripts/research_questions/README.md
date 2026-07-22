# Research-question analysis layer

This directory implements the active RQ01--RQ10 release menu. All analytical
outputs use the exact selected Longcycler cohort and are written below
`results/research_questions/`.

## Non-negotiable contracts

- The entire `Rowenas analysis/` tree is out of scope and must never be used as
  an input.
- Analytical episodes are the 532 selected, QC-passing Longcycler-linked rows
  from 161 residents in `results/clinical/analysis_cohort_longcycler.csv`:
  16 operational UTI and 516 operational Not_UTI episodes.
- The 583-episode clinical source is used only for clearly labelled attrition
  and QC context; it is not an analytical denominator.
- The assembly manifest must contain the same 532 participant-timepoint keys,
  exact FASTA paths, and no assembler fallback.
- The longitudinal universes are exactly 893 direct within-resident pairs and
  371 adjacent pairs rebuilt after cohort restriction. The single canonical
  adjacent-pair source is `results/longitudinal/longcycler_transitions.csv`.
- Longitudinal relatedness uses direct endpoint comparisons. Graph-component
  strain identifiers are never used as the primary outcome.
- The `<=25` DNAdiff SNP rule is an operational reference, not a validated
  biological gold standard.
- DNAdiff evidence must retain the selected endpoint paths and hashes and pass
  report-hash, signature, and sidecar validation. Validated reused caches are
  accepted; reconstructed legacy paths and unvalidated fallbacks are not.
- Key-linked provider ST calls are accepted by exact selected FASTA path; a
  provider assembler field is neither required nor published.
- Bootstrap resampling is performed at resident level, with seed `20260712` and
  10,000 replicates for final results.
- Research-facing case tables use generated case labels and never expose raw
  participant identifiers.

## Entry points

- `run_rq01_05.R`: verified clinical and direct-transition questions.
- `run_rq06_08.R`: VF stability, event-matched VF summaries, and relatedness
  surrogate performance.
- `run_rq09_10.R`: resident-specific clustering and adjacent provider-ST
  turnover with hash-bound direct-pair validation.
- `run_all.R`: sequential runner and final contract verification.

Use `RQ_BOOTSTRAP_REPS` to reduce bootstrap replicates for development only.
Final publication outputs must be generated with the default value of 10,000.

## Interpretation boundary

All estimates are conditional on an *E. coli*-positive episode and, for genomic
questions, successful sequencing and QC. The implemented `UTI_Status` is an
operational phenotype. These analyses do not estimate cohort UTI incidence,
the effect of asymptomatic bacteriuria on future UTI, transmission, treatment
effects, or causal effects of individual genomic features.

# Research-question analysis layer

This directory implements the prespecified RQ01--RQ11 menu without changing the
numbered clinical/genomic pipeline or its existing outputs. All generated files
are written below `results/research_questions/`.

## Non-negotiable contracts

- The entire `Rowenas analysis/` tree is out of scope and must never be used as
  an input.
- Clinical analyses start from the 583 primary-eligible episodes in
  `results/clinical/status_map.csv`.
- Primary genomic analyses start from the 532 selected, QC-passing Longcycler
  rows in `results/qc/analysis_assembly_manifest.csv`.
- Longitudinal relatedness uses direct endpoint comparisons. Graph-component
  strain identifiers are never used as the primary outcome.
- The `<=25` DNAdiff SNP rule is an operational reference, not a validated
  biological gold standard.
- Bootstrap resampling is performed at resident level (or isolate level for the
  explicitly paired assembler analysis), with seed `20260712` and 10,000
  replicates for final results.
- Research-facing case tables use generated case labels and never expose raw
  participant identifiers.

## Entry points

- `run_rq01_05.R`: verified clinical and direct-transition questions.
- `run_rq06_08.R`: VF stability, event-matched VF summaries, and relatedness
  surrogate performance.
- `run_rq09_11.R`: resident-specific clustering, ST turnover, and paired
  assembler concordance.
- `run_all.R`: sequential runner and final contract verification.

Use `RQ_BOOTSTRAP_REPS` to reduce bootstrap replicates for development only.
Final publication outputs must be generated with the default value of 10,000.

## Interpretation boundary

All estimates are conditional on an *E. coli*-positive episode and, for genomic
questions, successful sequencing and QC. The implemented `UTI_Status` is an
operational phenotype. These analyses do not estimate cohort UTI incidence,
the effect of asymptomatic bacteriuria on future UTI, transmission, treatment
effects, or causal effects of individual genomic features.

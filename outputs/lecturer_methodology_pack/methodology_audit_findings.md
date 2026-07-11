# Methodology audit findings

Date: 10 July 2026  
Purpose: internal scientific-validity and reproducibility review for lecturer discussion

## Executive finding

The project contains a coherent clinical-to-genomic analysis structure, but the defensible description is narrower than several older decks and handouts imply.

The 556-row mixed canonical dataset is the analysis that was executed. It contains 532 selected Longcycler assemblies and 24 selected Flye fallback assemblies. The 532-row Longcycler-only dataset is a useful single-assembler sensitivity analysis, not yet a replacement description of the entire executed pipeline.

Clinical, assembly-selection, MLST-linkage, VF-feature and transition denominators can be reproduced from current generated files. An earlier path-only `dnadiff` cache was unsafe because it did not prove that a report belonged to the current endpoint FASTAs. That historical output remains withdrawn and archived as audit evidence. The corrected cache-safe workflow has now recomputed the current canonical comparisons and verifies every reused report against exact FASTA paths and SHA-256 hashes.

The current provenance-validated longitudinal results are:

- mixed adjacent threshold-supported transitions: 138 of 394;
- Longcycler adjacent threshold-supported transitions: 140 of 371; and
- Longcycler Not_UTI-to-UTI threshold-supported transitions: five of nine.

This audit concerns reproducibility and scientific interpretation. It does not allege intent or misconduct.

## 1. Factual reproducibility

### Reproducible denominator anchors

| Finding | Current evidence | Audit status |
|---|---|---|
| 585 clinical records before primary exclusions | Current clinical status output | Reproduced |
| 583 primary-included episodes from 166 participants | `analysis_include_primary == TRUE` | Reproduced |
| 18 UTI and 565 Not_UTI included clinical episodes | Current primary `UTI_Status` | Reproduced |
| Two excluded records | One secondary duplicate and one unknown-participant record | Reproduced |
| 1,303 discovered candidate FASTAs | `results/qc/metadata_fasta_discovery_manifest.csv` | Reproduced; these are files, not independent episodes |
| 1,299 assembly-metadata rows | `results/assembly_metadata.csv` | Reproduced as a reconciliation boundary: 1,295 rows have matched, existing FASTAs and four are audit-only rows without a FASTA; eight discovered FASTAs are absent from the metadata table |
| 1,291 primary assembly candidates representing 579 episode keys | `results/qc/canonical_assembly_selection.csv` | Reproduced after primary/genomics curation; these are candidate alternatives, not independent episodes |
| 1,211 QC-passing assembly records | Implemented QC flags | Reproduced |
| 556 selected canonical rows from 162 participants | `selected_canonical == TRUE` and QC pass | Reproduced |
| Mixed selected composition: 532 Longcycler and 24 Flye fallback rows | Selected assembler field | Reproduced |
| Mixed clinical composition: 17 UTI and 539 Not_UTI rows | Clinical status joined to selected rows | Reproduced |
| 532 Longcycler-selected rows from 161 participants | Selected canonical, QC-pass and assembler Longcycler | Reproduced |
| Longcycler clinical composition: 16 UTI and 516 Not_UTI rows | Clinical status joined to sensitivity rows | Reproduced |
| 514 typed Longcycler-selected episode rows across 80 ST labels | Preferred MLST table joined by episode key | Reproduced |
| 227 VFDB-derived binary features | Analysis-ready VF matrix | Reproduced |
| 394 mixed adjacent transitions from 144 participants | Full-cohort transition output | Reproduced |
| 138 of 394 mixed transitions at no more than 25 `dnadiff` SNPs | Current provenance-validated pairwise output linked to mixed transitions | Reproduced |
| 371 rebuilt Longcycler transitions from 139 participants | Timeline rebuilt after Longcycler restriction | Reproduced |
| 140 of 371 Longcycler transitions at no more than 25 `dnadiff` SNPs | Current provenance-validated pairwise output linked to rebuilt transitions | Reproduced |
| Nine Longcycler Not_UTI-to-UTI transitions, five at no more than 25 `dnadiff` SNPs | Rebuilt clinical transition labels and current pairwise output | Reproduced |

### Filters and formula that must accompany the numbers

- Clinical counts use the current primary status and exclude records where `analysis_include_primary` is not true.
- The 1,303-to-1,299 step is metadata reconciliation rather than a four-file exclusion: eight discovered FASTAs are not represented in the metadata table, while four metadata rows are retained for audit visibility despite lacking a FASTA.
- Canonical genomic counts require `selected_canonical == TRUE`, QC pass and an existing selected FASTA.
- Longcycler sensitivity counts additionally require the selected assembler to be Longcycler after case normalisation.
- Transition outputs must be restricted to `cohort == "all"` where cohort-stratified copies coexist.
- Within each participant, `n` retained episodes produce `n - 1` adjacent retained-episode transitions.
- Rows, participants, assembly records, features, ST labels, all within-participant pairs and adjacent transitions are different units.

### Clinical definition

The current primary UTI label requires culture support plus compatible catheter-aware symptoms. Culture support uses at least 10³ CFU/mL where quantitative data are available.

Not_UTI is not a healthy or bacteria-free control category. The current included Not_UTI rows are bacteriuria episodes that do not meet the combined UTI definition. This heterogeneity must be stated when interpreting comparisons.

### Implemented assembly QC

The current QC script applies:

- FASTA readability;
- genome size from 4 to 6 megabases;
- no more than 200 contigs; and
- N50 of at least 20 kilobases.

Completeness and contamination thresholds appear as configuration constants but are not applied in the current QC implementation. They must not be reported as executed criteria without code changes and rerunning.

## 2. Pairwise-cache reproducibility incident and completed repair

### What was found

The earlier pairwise analysis generated 1,018 rows and reused existing `dnadiff` reports using a cache key based on participant/timepoint identifiers rather than the exact FASTA contents.

- Only 14 reports matched both current selected canonical FASTA paths.
- 949 rows had current selected files available but the cached report referred to different assembly inputs.
- A further 55 rows included at least one endpoint without a selected canonical assembly because the resolver fell back to another usable FASTA.
- In total, 1,004 of the 1,018 earlier rows were invalid for the current canonical-input claim.
- Within the 371 Longcycler-only adjacent transitions, only four cached reports matched both current Longcycler inputs; 367 did not.

Therefore, every prior SNP-derived headline count based on that historical output was withdrawn. The old cache and outputs are preserved only as audit evidence; they are not sources for the current presentation values.

### Why it happened

The former helper logic:

- named cached reports from participant/timepoint keys only;
- reused an existing report without verifying current FASTA paths or content hashes; and
- allowed fallback from a missing canonical endpoint to another usable assembly.

When canonical selection changed, the cache filename did not change. Old reports could therefore be interpreted as if they represented current inputs.

### Implemented repair and verification

The corrected comparison workflow now does the following:

1. Resolves only QC-passing selected canonical endpoints for the main mixed analysis.
2. Computes a content hash for each selected FASTA.
3. Includes both endpoint hashes in a deterministic pair signature.
4. Stores a sidecar recording paths, hashes, sizes, modification times, software version, command and generation time.
5. Reuses a report only when the sidecar exactly matches both current inputs.
6. Records endpoint assembler, path, hash, report and signature in the output.
7. Fails the provenance audit if any report is stale or any endpoint is noncanonical.

With the current 556 selected canonical rows, the workflow produces exactly 963 unordered within-participant canonical pairs. All endpoints are current QC-passing selected canonical assemblies, all reports have matching provenance sidecars and no report with a mismatched signature is reused. This is an all-pair denominator, not the adjacent-transition denominator.

### Current provenance-validated results

The lecturer materials should use the following values, derived only after joining the current pairwise output to the relevant rebuilt adjacent-transition tables:

| Current result | Reproduced value |
|---|---|
| Mixed adjacent transitions at or below the operational boundary | 138 of 394 |
| Longcycler adjacent transitions at or below the operational boundary | 140 of 371 |
| Longcycler Not_UTI-to-UTI transitions at or below the operational boundary | Five of nine |

The historical cached summary must not be used to replace or supplement these values.

## 3. Scientific interpretation

### Executed mixed analysis versus Longcycler sensitivity

Factual statement:

> The executed canonical analysis uses 556 selected rows: 532 Longcycler assemblies and 24 Flye fallback assemblies.

Supported interpretation:

> Restricting to the 532 selected Longcycler rows is a high-retention sensitivity analysis that reduces one source of technical heterogeneity.

Unsupported interpretation:

> Longcycler-only is inherently true, unbiased, universally superior or guaranteed to turn every remaining difference into a biological difference.

The sensitivity exclusion can itself produce selection if Flye fallback episodes differ from retained episodes. Both analyses should be reported.

### MLST provenance and interpretation

In the Longcycler sensitivity set, 514 episode rows link to usable preferred MLST calls and contain 80 ST labels. The typed subset comprises 509 provider calls meeting the provider QC rule and five labelled local fallback calls; 18 sensitivity rows are untyped.

These calls are linked to Longcycler-selected episode keys, but they should not be described as MLST generated exclusively from Longcycler assemblies. Provider provenance frequently references Flye and/or Longcycler inputs.

MLST supports lineage context. It cannot independently establish same-strain persistence, particularly for common STs.

### Pairwise `dnadiff`, Parsnp and wgMLST

`pairwise_metrics.csv` is based on direct assembly-to-assembly MUMmer `dnadiff` comparisons, supplemented by Mash and other pair metrics. It is not the Parsnp core-genome distance table.

Parsnp and `snp-dists` form a separate multi-genome core-alignment branch. wgMLST is another distinct method based on allele differences under a named locus scheme. These approaches can provide directionally compatible evidence, but their numerical distances and thresholds are not interchangeable.

The at-most-25 `dnadiff` SNP boundary is a predefined project-level operational rule for strong same-strain support. No citation or empirical calibration is recorded in the current code. It must not be described as a universal, validated or inherently conservative biological cutoff.

### Longitudinal interpretation

The 371 Longcycler transitions are adjacent among retained Longcycler episodes in pipeline time order. Of these, 362 were already Longcycler-to-Longcycler adjacencies in the mixed timeline and nine become newly adjacent after Flye removal. This is why the timeline must be rebuilt rather than filtered.

Consequently, the 371-pair Longcycler table is not a simple subset of the 394-pair mixed table. Its 140 threshold-supported pairs can exceed the mixed table's 138 because the sensitivity restriction changes which sampled episodes are adjacent; the two numerators must be interpreted with their own denominators.

Adjacency among retained observations does not prove continuous carriage between visits. Unobserved intermediate strains, unequal elapsed time and sampling gaps remain possible.

### VF interpretation

The 227 columns are VFDB-derived detected/not-detected features. They are not all validated UPEC virulence determinants. Current curation diagnostics identify:

- 64 features absent from the curated gene map;
- 57 unassigned features; and
- 15 low-confidence mappings.

Feature presence does not establish expression, biological activity or disease causation. Apparent gain or loss means a change in detected presence/absence between assemblies; it is not proof of biological acquisition or deletion.

### Formal modelling and causality

The UTI outcome is sparse and participants contribute repeated observations. Episode-level Fisher tests, Wilcoxon tests and unclustered lineage comparisons should be treated as exploratory. Formal models need transparent handling of participant dependence, singular fits, separation, multiple testing, missingness and potential overfitting.

The current analyses do not justify causal statements that antibiotics caused replacement, carriage caused UTI, colonisation was protective, or host factors drove an observed genomic transition. Such claims require appropriate exposure data, temporal ordering, confounder control and a causal design.

## 4. Claim-linked presentation boundary

The lecturer presentation names a script only when it supplies an input, transformation or result explicitly discussed in a slide. Its headline methodology chain is therefore limited to:

- `00a_load_clean_clinical.R` and `00b_classify_episodes.R` for clinical cleaning and classification;
- `00_make_assembly_metadata.r` and `12a_wgs_qc.R` for assembly reconciliation, implemented QC and canonical selection;
- `02_gene_presence_analysis.R` and `06_MLST.R` for VFDB-derived features and lineage context;
- `11_compare_strains.R` plus its helper for provenance-verified `dnadiff` and Mash comparisons;
- `22_vf_build_analysis_dataset.R` for episode-level clinical/genomic integration;
- `24_vf_longitudinal_dynamics.R` and `scripts/rebuild_longcycler_sensitivity.R` for mixed and Longcycler-only timelines;
- `12b_core_snp.R` and `12c_panaroo.R` only where the separate core-genome and pangenome context is discussed; and
- `14_genotype_phenotype_model.R` only in the statistical-caution appendix.

Other scripts can be described briefly when they clarify an upstream assumption, supporting diagnostic or limitation, but they are not presented as result-generating stages unless their outputs support a displayed claim. The separate comprehensive register remains the place for script-by-script audit detail.

Claims about antibiotic resistance, demographics, host factors, wgMLST allele distances, exact mutation locations or causal mechanisms remain outside the headline methodology because the required data or validated analysis is not established here. The exploratory association model is not confirmatory: sparse UTI outcomes, repeated participants, data-dependent screening, singular fits, separation, multiple testing and overfitting all require explicit diagnostics and cautious interpretation.

## 5. Missing provenance

The current R scripts cannot supply all information needed for a complete sequencing and bioinformatics methods section. The following must be recovered from laboratory records, command logs, provider manifests or environment records.

### Sequencing and read processing

- Oxford Nanopore platform/device and flow-cell type;
- library-preparation kit and barcode kit;
- basecaller, version, model and basecalling mode;
- demultiplexing and adapter-trimming software/settings;
- read-quality, read-length and coverage summaries; and
- contamination screening before assembly.

### Assembly

- Longcycler and Flye versions;
- complete commands and parameters;
- input-read selection;
- polishing software, versions, number of rounds and read type;
- assembly-generation date and batch; and
- empirical comparison or rationale for assembler choice for each downstream task.

### Typing and annotation

- ABRicate version and VFDB release/date;
- SeqSphere/provider software version and MLST scheme version;
- definition and provenance of the provider `PercGoodTargets >= 95` rule;
- local MLST software/database version;
- Prokka and Panaroo versions and complete commands.

### Relatedness analysis

- MUMmer/`dnadiff`, Mash, Parsnp and `snp-dists` versions;
- reference-selection provenance for Parsnp;
- exact comparison commands and exclusion rules;
- citation or study-specific calibration for the operational 25-SNP boundary; and
- any recombination handling used before biological interpretation.

### Computational reproducibility

- R version and package versions/session information;
- operating system/container or environment manifest;
- immutable input manifest with hashes;
- code revision used for each final output; and
- a complete run log connecting source files to tables and figures.

### Runner boundary

`RUN_COMPLETE_ANALYSIS.sh` assumes existing FASTAs and `assembly_metadata.csv`. It does not run sequencing, basecalling, read QC, Longcycler/Flye assembly or the metadata-builder script. The methods must not call it an end-to-end raw-read pipeline unless those upstream stages are integrated.

## 6. Wording corrections for all lecturer materials

| Replace | With |
|---|---|
| Framing the Longcycler-only subset as the sole primary denominator | “The 556-row mixed canonical analysis was executed; 532 Longcycler-selected rows form a single-assembler sensitivity analysis.” |
| “Longcycler-only supports reproducible methods” | “Longcycler-only reduces one source of technical heterogeneity and tests robustness to assembler restriction.” |
| “core-genome SNPs” for `pairwise_metrics.csv` | “assembly-to-assembly `dnadiff` SNP differences” |
| “the strict/conservative same-strain rule” | “a predefined operational at-most-25 `dnadiff` SNP boundary for strong support that still requires justification/calibration” |
| Implying that MLST was generated only from Longcycler | “preferred MLST calls linked to Longcycler-selected episode rows” |
| “consecutive transitions” without qualification | “adjacent among retained episodes in pipeline time order” |
| “227 virulence genes” | “227 VFDB-derived detected/not-detected gene features” |
| “gene acquisition/deletion” from presence/absence change | “detected feature gain/loss between retained assemblies” |
| “Not_UTI controls” or “healthy episodes” | “bacteriuria episodes not meeting the combined UTI definition” |
| “mechanism” or “causal driver” | “descriptive pattern,” “case-prioritisation hypothesis” or “association,” as appropriate |

The older four-slide Longcycler presentation, its handout, its detailed talking points and its layperson explanation are superseded for lecturer use because they contain the withdrawn cache-derived results and the incorrect primary-analysis framing. The separate `results/thesis_audit/hamdi_thesis_discussion_brief_today.md` and `.docx` are also superseded until regenerated, because they still contain the pre-repair strict-transition count; the current scientific audit and calculation-filter reconciliation are the authoritative Hamdi-audit outputs.

## 7. Lecturer decisions required

1. Confirm the reporting hierarchy: executed mixed canonical analysis, followed by Longcycler-only sensitivity analysis.
2. Decide whether the operational SNP boundary can be justified from literature and current data or requires a threshold-sensitivity presentation rather than a binary label.
3. Approve exclusion of unvalidated persistence, mutation-location, antibiotic-resistance and causal claims from headline results.
4. Decide which sparse association models remain exploratory supplementary analyses.
5. Recover and approve the missing sequencing, assembly, database and software provenance before thesis submission.

## Audit conclusion

The strongest currently reproducible parts of the methodology are the clinical classification, implemented assembly-QC/canonical selection, selected-row denominators, VFDB-derived feature matrix, preferred MLST linkage, provenance-verified pairwise comparison and reconstructed longitudinal adjacency. Scientific conclusions should remain descriptive and lineage-aware, with clear separation between factual reproduction, interpretation, sensitivity analysis and causal speculation.
